// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file testDeviceTracker.cxx
/// \brief A timeframe on a device against the host tracker's, pass by pass
///        and iteration by iteration.
///
/// Built once per backend. The host tracker and the device run one plan in
/// lockstep: the host's slice observer runs the same slice on the device,
/// its iteration observer the device's iteration end, each with its own used
/// mask; every pass's output, every iteration's selection, proposals, tracks
/// and mask, and the timeframe's tracks are compared: byte for byte in the
/// deterministic build, up to the last bits of the float columns with fast
/// math. The device also runs a timeframe
/// alone, and its lifecycle and failures are checked. Needs the device.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK device tracker
#include <boost/test/unit_test.hpp>
#include <numbers>

#include "DeviceTestFixture.h"
#include "DeviceTestKernels.h"
#include "FixedPointReferences.h"
#include "ToyBarrel8Fixture.h"

#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceTracker.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Study/StageComparison.h"
#include "GenericTracking/Tracking/Tracker.h"
#include "GenericTracking/Tracking/VertexingPass.h"
#include "UniformMaterial.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk
{
namespace
{

using namespace device_test;

using namespace toybarrel8;
using Host = Tracker<Det, UniformMaterial>;

/// Twelve tracks on the sensors over 1.1 us, the even ones stiff (2 to 3
/// GeV), the odd ones soft (0.4 to 0.7 GeV): the first iteration takes the
/// stiff ones, and the second, down to 0.3 GeV, works around their rows.
constexpr std::size_t kTracks = 12U;
constexpr double kDurationNS = 1100.;

/// The layers from which the soft tracks get twins (`hits`).
constexpr LayerID kFirstTwinLayer = 5;

/// With \a twins, a second hit 5 um beside each soft track's on the layers
/// from `kFirstTwinLayer` out, on the same sensor and at the same time: two
/// rows for the follow to choose between (its lookahead), where an iteration
/// whose links stop below those layers leaves them to it.
std::vector<Hit> hits(bool twins = false)
{
  std::vector<Hit> result;
  for (std::size_t track = 0; track < kTracks; ++track) {
    const auto k = static_cast<double>(track);
    const bool stiff = track % 2U == 0U;
    const double pt = stiff ? 2. + (0.2 * static_cast<double>(track % 5U)) : 0.4 + (0.06 * static_cast<double>(track % 5U));
    const auto one = sensorHelixHits(std::fmod(0.3 + (0.9 * k), 2 * std::numbers::pi_v<double>), 0.03 * std::cos(1.3 * k), pt, 37. + (83. * k));
    result.insert(result.end(), one.begin(), one.end());
    if (!twins || stiff) {
      continue;
    }
    for (const auto& hit : one) {
      if (hit.point.layer < kFirstTwinLayer) {
        continue;
      }
      Hit twin = hit;
      twin.point = moved(hit.point, 5.e-4, 5.e-4);
      const Scalar global[3]{twin.point.x, twin.point.y, twin.point.z};
      Scalar local[3]{};
      surfaceToLocal(trackingGeometry().surface(hit.surface), global, local);
      twin.loc0 = local[0];
      twin.loc1 = local[1];
      result.push_back(twin);
    }
  }
  return result;
}

std::vector<IterationSettings> iterations(VertexSource source)
{
  auto first = iterationSettings();
  auto second = iterationSettings();
  second.windows.minPtGeV = 0.3;
  first.vertexSource = source;
  second.vertexSource = source;
  return {first, second};
}

/// One vertex per track, at its start time on the beam axis: a point box.
std::vector<VertexSeed> trackVertices()
{
  std::vector<VertexSeed> vertices;
  vertices.reserve(kTracks);
  for (std::size_t track = 0; track < kTracks; ++track) {
    const auto timePS = static_cast<std::uint32_t>(37'000U + (83'000U * track));
    vertices.push_back({.timePS = timePS, .widthPS = 0U, .z = 0.F, .varZ = 4.e-4F, .varTransverse = 1.e-6F});
  }
  return vertices;
}

/// The device's output is the host's whatever its schedule, in the
/// deterministic build byte for byte. With fast math, the same
/// records, whose float columns may differ in their last bits.
template <template <template <class> class> class S>
bool matchesHost(const StageOutput<S>& host, const StageOutput<S>& device)
{
  const auto comparison = study::compareStageOutputs(host, device);
#ifdef GPUCA_DETERMINISTIC_MODE
  return comparison.identical();
#else
  constexpr double kTolerance = 1.e-4;
  return comparison.sameElements && comparison.elementsWithOtherCount == 0U && comparison.elementsReordered == 0U && std::ranges::all_of(comparison.maxAbsDifference, [&](double difference) { return difference <= kTolerance; });
#endif
}

struct Totals {
  std::size_t iterations = 0U;
  std::size_t slices = 0U;
  std::size_t tracklets = 0U;
  std::size_t triplets = 0U;
  std::size_t neighbours = 0U;
  std::size_t chained = 0U;
  std::size_t roads = 0U;
  std::size_t records = 0U;
  std::size_t gathered = 0U;
  std::size_t tracks = 0U;
  std::size_t rounds = 0U;
  std::size_t extended = 0U;     ///< tracks the follow added clusters to, on the device
  std::size_t maskedSlices = 0U; ///< slices that read a used mask with exclusive rows
  bool missed = false;
};

/// A configuration of both trackers.
struct Options {
  VertexSource source = VertexSource::kBeamLine;
  std::int32_t width = 13;
  std::uint32_t sharedInnermost = 0U;
  bool extension = true;
  /// The second iteration's kink allowance, whose follow runs in lane groups
  /// on the device, its links below
  /// `kFirstTwinLayer` and the soft tracks' twins there for the follow
  /// (`hits`), and the first's end decided again:
  /// never both in one iteration.
  double kinkSigma = 1.;
  bool followEnd = false;
  FitReference reference = FitReference::kFilteredState;
};

std::string describe(const Options& options) { return std::string{options.source == VertexSource::kBeamLine ? "beam line" : "vertex list"} + ", width " + std::to_string(options.width) + ", shared " + std::to_string(options.sharedInnermost) + (options.extension ? ", follow" : ", no follow") + (options.kinkSigma > 1. ? ", kink allowance " + std::to_string(options.kinkSigma) + " sigma" : "") + (options.followEnd ? ", end again" : ""); }

std::unique_ptr<Host> makeHost(const Options& options)
{
  const std::array<double, kLayers> material{};
  static const UniformMaterial none{};
  auto settings = iterations(options.source);
  for (auto& iteration : settings) {
    iteration.fitReference = options.reference;
    iteration.sharedInnermost = options.sharedInnermost;
    iteration.follow = options.extension;
  }
  settings.front().followEnd = options.followEnd;
  settings.back().kinkSigma = options.kinkSigma;
  if (options.kinkSigma > 1.) {
    settings.back().layerMask = (std::uint64_t{1} << kFirstTwinLayer) - 1U;
  }
  auto tracker = std::make_unique<Host>(settings, Fixture::kBeamConstraint, scatteringTable<Det>(material), std::make_shared<BoundedMemoryResource>(), 2U, 2);
  tracker->setTrackModel(trackingGeometry(), none, Propagator{});
  return tracker;
}

std::vector<IterationParameters<Det>> parametersOf(const Host& host)
{
  std::vector<IterationParameters<Det>> parameters;
  parameters.reserve(host.iterationCount());
  for (std::size_t iteration = 0; iteration < host.iterationCount(); ++iteration) {
    parameters.push_back(host.iteration(iteration));
  }
  return parameters;
}

void configure(DeviceTracker<Det, UniformMaterial>& device, const Host& host) { BOOST_REQUIRE((device.configureRun(parametersOf(host), trackingGeometry(), UniformMaterial{}, Propagator{}) == StageStatus::kComplete)); }

[[maybe_unused]] bool sameProposal(const FollowProposal<maxChainLayers<Det>>& lhs, const FollowProposal<maxChainLayers<Det>>& rhs) { return lhs.outwardCount == rhs.outwardCount && lhs.inwardCount == rhs.inwardCount && std::ranges::equal(lhs.outward, rhs.outward) && std::ranges::equal(lhs.inward, rhs.inward) && lhs.replacedEnd == rhs.replacedEnd; }

/// The same for tracks: byte for byte in the deterministic build, the states field
/// by field. With fast math a chi2's last bits may reorder ties, so the same
/// tracks -- rows, count, iteration -- in any order, each with the same time
/// box: its rows', not the fit's.
template <typename View>
bool sameTracks(const View& host, std::size_t hostRows, const View& device, std::size_t deviceRows)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return study::compareColumns<TrackOf<maxChainLayers<Det>>::RecordT>(host, hostRows, device, deviceRows).identical();
#else
  const auto rowsOf = [](const View& view, std::size_t rows) {
    std::vector<std::vector<std::uint64_t>> result;
    result.reserve(rows);
    for (std::size_t i = 0; i < rows; ++i) {
      std::vector<std::uint64_t> track{view.clusters[i].index, view.clusters[i].index + view.count[i]};
      track.push_back(view.iteration[i]);
      track.push_back(view.timePS[i]);
      track.push_back(view.widthPS[i]);
      result.push_back(std::move(track));
    }
    std::ranges::sort(result);
    return result;
  };
  if (hostRows != deviceRows || rowsOf(host, hostRows) != rowsOf(device, deviceRows)) {
    return false;
  }
  // Check states explicitly: a non-arithmetic memlayout column has no scalar
  // maxAbsDifference in StageComparison, so matching identities is not a
  // numerical check of its parameters and covariance.
  const auto order = [](const View& view, std::size_t rows) {
    std::vector<std::pair<std::vector<std::uint32_t>, std::size_t>> sorted;
    for (std::size_t i = 0; i < rows; ++i) {
      sorted.emplace_back(std::vector<std::uint32_t>{view.clusters[i].index, view.clusters[i].index + view.count[i]}, i);
    }
    std::ranges::sort(sorted);
    return sorted;
  };
  const auto lhs = order(host, hostRows);
  const auto rhs = order(device, deviceRows);
  const auto closeState = [](const BoundTrackState& a, const BoundTrackState& b) {
    constexpr Scalar tolerance = 1.e-4F;
    if (a.surface != b.surface || std::abs(a.surfacePhi - b.surfacePhi) > tolerance) {
      return false;
    }
    for (int i = 0; i < kNBoundParams; ++i) {
      Scalar delta = a.par[i] - b.par[i];
      if (i == kPhi) {
        o2::math_utils::detail::bringToPMPi<Scalar>(delta);
      }
      if (!(std::abs(delta) <= tolerance)) {
        return false;
      }
    }
    return std::ranges::equal(a.cov, b.cov, [=](Scalar x, Scalar y) { return std::abs(x - y) <= tolerance; });
  };
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    if (!closeState(host.state[lhs[i].second], device.state[rhs[i].second]) || !closeState(host.outwardState[lhs[i].second], device.outwardState[rhs[i].second])) {
      return false;
    }
  }
  return true;
#endif
}

/// The same indices, in the same order in the deterministic build, in any
/// order with fast math.
bool sameIndices(std::span<const std::uint32_t> host, std::span<const std::uint32_t> device)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return std::ranges::equal(host, device);
#else
  std::vector<std::uint32_t> lhs{host.begin(), host.end()};
  std::vector<std::uint32_t> rhs{device.begin(), device.end()};
  std::ranges::sort(lhs);
  std::ranges::sort(rhs);
  return lhs == rhs;
#endif
}

/// The proposals, one per track: equal in the deterministic build; with fast
/// math, in the tracks' order, which may differ, so the clusters they add.
template <typename Proposals>
bool sameProposals(const Proposals& host, const std::vector<FollowProposal<maxChainLayers<Det>>>& device)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return std::ranges::equal(host, device, sameProposal);
#else
  const auto added = [](const auto& proposals) {
    std::size_t sum = 0U;
    for (const auto& proposal : proposals) {
      sum += std::size_t{proposal.outwardCount} + proposal.inwardCount;
    }
    return sum;
  };
  return host.size() == device.size() && added(host) == added(device);
#endif
}

/// The host tracker's whole timeframe and the device's, in lockstep over one
/// plan: every slice compared pass by pass, every iteration's selection,
/// proposals, tracks and used mask, then the timeframe's tracks.
Totals compare(DeviceTracker<Det, UniformMaterial>& device, const Options& options)
{
  const Fixture fixture{hits(options.kinkSigma > 1.), kDurationNS};
  const auto vertexList = trackVertices();
  const std::span<const VertexSeed> vertices = options.source == VertexSource::kVertexList ? std::span<const VertexSeed>{vertexList} : std::span<const VertexSeed>{};
  const auto tracker = makeHost(options);
  configure(device, *tracker);
  const auto plan = tracker->planTimeframe(*fixture.frame, vertices, options.width, 0.);
  BOOST_REQUIRE((device.beginTimeframe(*fixture.frame, plan) == StageStatus::kComplete));

  Totals totals;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  bool open = false;
  std::size_t openIteration = 0U;
  std::size_t nextSlice = 0U;
  tracker->setIterationObserver([&](const Host&, const Host::IterationObservation& iteration) {
    BOOST_REQUIRE(open);
    BOOST_TEST(openIteration == iteration.iteration);
    DeviceIterationDownload<Det> onDevice{pool};
    BOOST_REQUIRE((device.endIteration(&onDevice) == StageStatus::kComplete));
    open = false;
    const auto& selection = *iteration.selection;
    const auto& report = device.lastIteration();
    BOOST_TEST_CONTEXT("iteration " << iteration.iteration)
    {
      // The records of every slice, gathered in slice order: the host's appended ones.
      const auto records = study::compareColumns<SelectionRecordOf<maxChainLayers<Det>>::RecordT>(selection.records(), selection.recordCount(), std::as_const(onDevice.records).view(), onDevice.records.size());
#ifdef GPUCA_DETERMINISTIC_MODE
      BOOST_TEST(records.identical());
      BOOST_TEST(std::ranges::equal(selection.ranked(), onDevice.ranked));
      BOOST_TEST(report.rounds == selection.rounds());
#else
      BOOST_TEST(records.sameElements);
#endif
      BOOST_TEST(report.records == selection.recordCount());
      BOOST_TEST(report.duplicates == selection.duplicates());
      BOOST_TEST(report.mismatches == 0U);
      BOOST_TEST(sameIndices(selection.selected(), onDevice.selected));
      BOOST_TEST(report.tracks == selection.selected().size());
      BOOST_TEST(std::ranges::equal(iteration.used, onDevice.used));
      BOOST_TEST(sameProposals(iteration.proposals, onDevice.proposals));
      BOOST_TEST(onDevice.proposals.size() == (options.extension ? onDevice.selected.size() : 0U));
      const auto hostTracks = std::as_const(*iteration.tracks).view() + static_cast<std::ptrdiff_t>(iteration.trackOffset);
      BOOST_TEST(sameTracks(hostTracks, iteration.tracks->size() - iteration.trackOffset, std::as_const(onDevice.tracks).view(), onDevice.tracks.size()));
    }
    totals.gathered += onDevice.records.size();
    totals.tracks += onDevice.tracks.size();
    totals.rounds += report.rounds;
    totals.extended += static_cast<std::size_t>(std::ranges::count_if(onDevice.proposals, [](const auto& proposal) { return proposal.outwardCount + proposal.inwardCount != 0U; }));
    ++totals.iterations;
  });
  tracker->setSliceObserver([&](const Host& host, const Host::SliceObservation& slice) {
    if (!open || openIteration != slice.iteration) {
      device.beginIteration(slice.iteration);
      open = true;
      openIteration = slice.iteration;
      nextSlice = 0U;
    }
    BOOST_REQUIRE(nextSlice < plan.slices.size());
    BOOST_TEST(slice.slice == &plan.slices[nextSlice].slice);
    DeviceSliceDownload<Det> onDevice{pool};
    BOOST_REQUIRE((device.runSlice(nextSlice, &onDevice) == StageStatus::kComplete));
    ++nextSlice;
    const auto& result = *slice.result;
    BOOST_TEST_CONTEXT("iteration " << slice.iteration << ", slice from block " << slice.slice->coreBlocks.begin)
    {
      BOOST_TEST(matchesHost(*result.tracklets, onDevice.tracklets));
      BOOST_TEST(matchesHost(*result.triplets, onDevice.triplets));
      BOOST_TEST(matchesHost(*result.neighbours, onDevice.neighbours));
      BOOST_TEST(std::ranges::equal(chainedTriplets(*result.neighbours), onDevice.chained));
      BOOST_TEST(std::ranges::equal(slice.roads->successorAbove(), onDevice.successorAbove));
      BOOST_TEST(std::ranges::equal(slice.roads->roads(), onDevice.roads));
      BOOST_TEST(std::ranges::equal(slice.roads->starts(), onDevice.starts));
      BOOST_TEST(std::ranges::equal(slice.roads->roadOffsets(), onDevice.roadOffsets));
      BOOST_TEST(matchesHost(*slice.records, onDevice.records));
    }
    const auto& report = device.lastSlice();
    BOOST_TEST(report.tracklets == result.tracklets->size());
    BOOST_TEST(report.triplets == result.triplets->size());
    BOOST_TEST(report.neighbours == result.neighbours->size());
    BOOST_TEST(report.readbacks >= 6U);
    ++totals.slices;
    totals.tracklets += report.tracklets;
    totals.triplets += report.triplets;
    totals.neighbours += report.neighbours;
    totals.chained += report.chained;
    totals.roads += report.roads;
    totals.records += report.records;
    totals.missed = totals.missed || report.stagingPredictionMissed;
    totals.maskedSlices += std::ranges::any_of(host.used(), [](std::uint8_t state) { return state == kExclusive; }) ? 1U : 0U;
  });
  const auto result = tracker->run(*fixture.frame, plan);
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  BOOST_TEST(result.tracks->size() > 0U);
  BOOST_TEST(totals.iterations == tracker->iterationCount());
  // The timeframe's tracks: every iteration's, in order.
  typename TrackRefit<Det>::Tracks tracks;
  tracks.setMemoryResource(pool.get());
  device.downloadTracks(tracks);
  BOOST_TEST(sameTracks(std::as_const(*result.tracks).view(), result.tracks->size(), std::as_const(tracks).view(), tracks.size()));
  return totals;
}

std::vector<Options> configurations()
{
  std::vector<Options> result;
  for (const auto source : {VertexSource::kBeamLine, VertexSource::kVertexList}) {
    for (const std::int32_t width : {5, 13, 1000}) {
      result.push_back({.source = source, .width = width});
    }
    result.push_back({.source = source, .width = 13, .sharedInnermost = 2U});
    result.push_back({.source = source, .width = 13, .extension = false});
    result.push_back({.source = source, .width = 13, .kinkSigma = 5., .followEnd = true});
  }
  return result;
}

} // namespace

BOOST_AUTO_TEST_CASE(EverySliceAndIterationMatchesTheHostTracker)
{
  auto* device = runtime();
  const auto tracker = kernels().makeToyBarrel8Tracker(*device, block(kArenaBytes));
  BOOST_TEST(tracker->arenaBytes() == kArenaBytes);
  for (const auto reference : {FitReference::kFilteredState, FitReference::kCircleOutward, FitReference::kIndependentBoth}) {
    for (auto options : configurations()) {
      options.reference = reference;
      BOOST_TEST_CONTEXT(describe(options) << ", reference " << static_cast<int>(reference))
      {
        const auto totals = compare(*tracker, options);
        BOOST_TEST(totals.triplets > 0U);
        BOOST_TEST(totals.neighbours > 0U);
        BOOST_TEST(totals.chained > 0U);
        BOOST_TEST(totals.roads > 0U);
        BOOST_TEST(totals.records > 0U);
        BOOST_TEST(totals.gathered == totals.records);
        BOOST_TEST(totals.tracks == kTracks);
        BOOST_TEST(totals.maskedSlices > 0U);
        if (options.kinkSigma > 1.) {
          BOOST_TEST(totals.extended > 0U);
        }
        BOOST_TEST_MESSAGE(device->backend() << " " << describe(options) << ": " << totals.slices << " slices, " << totals.tracklets << " tracklets, " << totals.triplets << " triplets, " << totals.neighbours << " neighbours, " << totals.roads << " roads, " << totals.records << " records, " << totals.rounds << " rounds, " << totals.tracks << " tracks, " << totals.extended << " extended");
      }
    }
  }
}

/// The flat runs in chunks of 16 pairs: sources and first tracklets split
/// across chunks, every chunk's elements re-walked past the others', the
/// chunks' outputs assembled; the slices are the host's as before.
BOOST_AUTO_TEST_CASE(ChunkedFlatRunsMatchTheHostTracker)
{
  auto* device = runtime();
  const auto tracker = kernels().makeToyBarrel8Tracker(*device, block(kArenaBytes));
  tracker->setFlatChunkPairs(16U);
  for (const auto& options : configurations()) {
    BOOST_TEST_CONTEXT(describe(options))
    {
      const auto totals = compare(*tracker, options);
      BOOST_TEST(totals.triplets > 0U);
      BOOST_TEST(totals.tracks == kTracks);
    }
  }
}

/// The staged stages' staging at one record -- the neighbours, the road
/// fit's records and, where a wave is 64 lanes wide, the tracklet stage's
/// lane groups; the flat runs (the triplet stage, CUDA's tracklet stage) list
/// their candidates exactly and stage nothing: each overflows, retries at the
/// exact size and keeps it; the outputs do not change.
BOOST_AUTO_TEST_CASE(MispredictedStagingGivesTheSameSlices)
{
  auto* device = runtime();
  const auto tracker = kernels().makeToyBarrel8Tracker(*device, block(kArenaBytes));
  tracker->setStagingCapacity(1U);
  const auto totals = compare(*tracker, {.width = 1000});
  BOOST_TEST(totals.missed);
  BOOST_TEST(totals.neighbours > 0U);
  BOOST_TEST(totals.records > 0U);
}

BOOST_AUTO_TEST_CASE(AnEdgelessTripletGraphNeedsNoOverflowReadback)
{
  auto rows = sensorHelixHits(0.3, 0.03, 2., 37.);
  rows.resize(3U);
  const Fixture fixture{rows, kDurationNS};
  auto settings = iterationSettings();
  settings.minClusters = 3;
  settings.follow = false;
  const std::array<double, kLayers> material{};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  Host host{settings, Fixture::kBeamConstraint, scatteringTable<Det>(material), pool, 2U, 2};
  const UniformMaterial none{};
  host.setTrackModel(trackingGeometry(), none, Propagator{});
  const auto plan = host.planTimeframe(*fixture.frame, {}, 1000, 0.);
  auto* device = runtime();
  const auto tracker = kernels().makeToyBarrel8Tracker(*device, block(kArenaBytes));
  configure(*tracker, host);
  tracker->setStagingCapacity(16U);
  BOOST_REQUIRE((tracker->beginTimeframe(*fixture.frame, plan) == StageStatus::kComplete));
  tracker->beginIteration(0U);
  DeviceSliceDownload<Det> output{pool};
  BOOST_REQUIRE((tracker->runSlice(0U, &output) == StageStatus::kComplete));
  const auto& report = tracker->lastSlice();
  BOOST_REQUIRE(report.triplets == 1U);
  BOOST_TEST(report.chained == 0U);
  // Four stages, the tracklet launch's list, permutation, the chained count and (with a download)
  // the chained list, starts, offsets: 10; a flat run
  // reads back three -- the list's length, the selected count, the
  // limit flag -- for the staging's one: the tracklets' and the triplets'; and the
  // count of the triplets kept on a path.
  BOOST_TEST(report.readbacks == 15U);
  BOOST_REQUIRE(output.roads.size() == 1U);
  BOOST_TEST(output.roads[0] == 1U);
  BOOST_TEST(output.successorAbove[0] == 0U);
  BOOST_TEST((tracker->endIteration(nullptr) == StageStatus::kComplete));
}

/// trackTimeframe alone, with no observer driving it: the host's tracks.
BOOST_AUTO_TEST_CASE(ARunAloneGivesTheHostsTracks)
{
  auto* device = runtime();
  const auto tracker = kernels().makeToyBarrel8Tracker(*device, block(kArenaBytes));
  const auto vertexList = trackVertices();
  const auto pool = std::make_shared<BoundedMemoryResource>();
  for (const auto& options : configurations()) {
    BOOST_TEST_CONTEXT(describe(options))
    {
      const Fixture fixture{hits(options.kinkSigma > 1.), kDurationNS};
      const std::span<const VertexSeed> vertices = options.source == VertexSource::kVertexList ? std::span<const VertexSeed>{vertexList} : std::span<const VertexSeed>{};
      const auto host = makeHost(options);
      configure(*tracker, *host);
      const auto plan = host->planTimeframe(*fixture.frame, vertices, options.width, 0.);
      const auto expected = host->run(*fixture.frame, plan);
      BOOST_REQUIRE(expected.status == StageStatus::kComplete);
      typename TrackRefit<Det>::Tracks tracks;
      tracks.setMemoryResource(pool.get());
      typename DeviceTracker<Det, UniformMaterial>::TimeframeReport report;
      BOOST_REQUIRE((tracker->trackTimeframe(*fixture.frame, plan, plan.parameters.size(), false, tracks, report) == StageStatus::kComplete));
      tracker->endTimeframe();
      BOOST_TEST(tracks.size() == kTracks);
      BOOST_TEST(sameTracks(std::as_const(*expected.tracks).view(), expected.tracks->size(), std::as_const(tracks).view(), tracks.size()));
    }
  }
}

namespace
{

/// A configured tracker, the fixture's timeframe and its plan, for the
/// lifecycle tests.
struct Configured {
  explicit Configured(std::size_t arenaBytes = kArenaBytes)
    : runtime(device_test::runtime()),
      tracker(kernels().makeToyBarrel8Tracker(*runtime, block(arenaBytes))),
      fixture(hits(), kDurationNS),
      host(makeHost({}))
  {
    parameters = parametersOf(*host);
    plan = host->planTimeframe(*fixture.frame, {}, 13, 0.);
  }

  [[nodiscard]] StageStatus configure() { return tracker->configureRun(parameters, trackingGeometry(), UniformMaterial{}, Propagator{}); }

  SliceScope slices;
  DeviceRuntime* runtime;
  std::unique_ptr<DeviceTracker<Det, UniformMaterial>> tracker;
  Fixture fixture;
  std::unique_ptr<Host> host;
  std::vector<IterationParameters<Det>> parameters;
  TimeframePlan<Det> plan;
};

} // namespace

/// Every call out of order throws.
BOOST_AUTO_TEST_CASE(CallsOutOfOrderAreRefused)
{
  Configured configured;
  auto& tracker = *configured.tracker;
  const auto& frame = *configured.fixture.frame;
  BOOST_CHECK_THROW(static_cast<void>(tracker.beginTimeframe(frame, configured.plan)), std::logic_error);
  BOOST_REQUIRE((configured.configure() == StageStatus::kComplete));
  BOOST_CHECK_THROW(tracker.beginIteration(0U), std::logic_error);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  typename TrackRefit<Det>::Tracks tracks;
  tracks.setMemoryResource(pool.get());
  BOOST_CHECK_THROW(tracker.downloadTracks(tracks), std::logic_error);
  BOOST_REQUIRE((tracker.beginTimeframe(frame, configured.plan) == StageStatus::kComplete));
  BOOST_CHECK_THROW(static_cast<void>(tracker.runSlice(0U, nullptr)), std::logic_error);
  BOOST_CHECK_THROW(static_cast<void>(tracker.endIteration(nullptr)), std::logic_error);
  BOOST_CHECK_THROW(tracker.beginIteration(configured.parameters.size()), std::invalid_argument);
  tracker.beginIteration(0U);
  BOOST_CHECK_THROW(tracker.beginIteration(1U), std::logic_error);
  BOOST_CHECK_THROW(static_cast<void>(tracker.runSlice(configured.plan.slices.size(), nullptr)), std::logic_error);
  BOOST_CHECK_THROW(tracker.downloadTracks(tracks), std::logic_error);
  BOOST_TEST((tracker.runSlice(0U, nullptr) == StageStatus::kComplete));
  BOOST_TEST((tracker.endIteration(nullptr) == StageStatus::kComplete));
  BOOST_CHECK_THROW(static_cast<void>(tracker.endIteration(nullptr)), std::logic_error);
  tracker.downloadTracks(tracks);
  // A new timeframe closes an open iteration, and the tracks.
  tracker.beginIteration(1U);
  BOOST_REQUIRE((tracker.beginTimeframe(frame, configured.plan) == StageStatus::kComplete));
  BOOST_CHECK_THROW(static_cast<void>(tracker.endIteration(nullptr)), std::logic_error);
  BOOST_TEST(tracker.arenaUsedBytes(ArenaEnd::kScratch) == 0U);
  tracker.downloadTracks(tracks);
  BOOST_TEST(tracks.size() == 0U);
}

/// A plan of another frame or another run is refused before anything
/// changes.
BOOST_AUTO_TEST_CASE(APlanOfAnotherFrameOrRunIsRefused)
{
  Configured configured;
  auto& tracker = *configured.tracker;
  BOOST_REQUIRE((configured.configure() == StageStatus::kComplete));
  BOOST_REQUIRE((tracker.beginTimeframe(*configured.fixture.frame, configured.plan) == StageStatus::kComplete));
  const auto base = tracker.arenaUsedBytes(ArenaEnd::kBase);
  const Fixture longer{hits(), 2. * kDurationNS};
  BOOST_CHECK_THROW(static_cast<void>(tracker.beginTimeframe(*longer.frame, configured.plan)), std::invalid_argument);
  const auto other = makeHost({.source = VertexSource::kVertexList});
  const auto otherPlan = other->planTimeframe(*configured.fixture.frame, {}, 13, 0.);
  BOOST_CHECK_THROW(static_cast<void>(tracker.beginTimeframe(*configured.fixture.frame, otherPlan)), std::invalid_argument);
  BOOST_TEST(tracker.arenaUsedBytes(ArenaEnd::kBase) == base);
  // The timeframe is still there.
  tracker.beginIteration(0U);
  BOOST_TEST((tracker.runSlice(0U, nullptr) == StageStatus::kComplete));
  BOOST_TEST((tracker.endIteration(nullptr) == StageStatus::kComplete));
}

/// A slice refused for memory fails its iteration: the next slice is
/// refused, endIteration returns the failure and leaves the scratch empty,
/// and the next iteration begins.
BOOST_AUTO_TEST_CASE(AFailedSliceFailsItsIteration)
{
  std::size_t timeframeBytes = 0U;
  {
    Configured probe;
    BOOST_REQUIRE((probe.configure() == StageStatus::kComplete));
    BOOST_REQUIRE((probe.tracker->beginTimeframe(*probe.fixture.frame, probe.plan) == StageStatus::kComplete));
    timeframeBytes = probe.tracker->arenaUsedBytes(ArenaEnd::kBase);
  }
  // Room for the run and the timeframe, and not for a slice's staging.
  Configured configured{timeframeBytes + std::size_t{64} * 1024U};
  auto& tracker = *configured.tracker;
  BOOST_REQUIRE((configured.configure() == StageStatus::kComplete));
  BOOST_REQUIRE((tracker.beginTimeframe(*configured.fixture.frame, configured.plan) == StageStatus::kComplete));
  const auto base = tracker.arenaUsedBytes(ArenaEnd::kBase);
  tracker.beginIteration(0U);
  BOOST_TEST((tracker.runSlice(0U, nullptr) == StageStatus::kOutOfMemory));
  BOOST_TEST(tracker.arenaUsedBytes(ArenaEnd::kBase) == base);
  BOOST_CHECK_THROW(static_cast<void>(tracker.runSlice(1U, nullptr)), std::logic_error);
  BOOST_TEST((tracker.endIteration(nullptr) == StageStatus::kOutOfMemory));
  BOOST_TEST(tracker.arenaUsedBytes(ArenaEnd::kScratch) == 0U);
  BOOST_TEST(tracker.arenaUsedBytes(ArenaEnd::kBase) == base);
  tracker.beginIteration(1U);
  BOOST_TEST((tracker.endIteration(nullptr) == StageStatus::kComplete));

  // A run the arena cannot hold is refused, and leaves no run configured.
  Configured tiny{4096U};
  BOOST_TEST((tiny.configure() == StageStatus::kOutOfMemory));
  BOOST_CHECK_THROW(static_cast<void>(tiny.tracker->beginTimeframe(*tiny.fixture.frame, tiny.plan)), std::logic_error);
}

namespace
{

/// The tracks, and 30,000 hits on the outermost layer after them in time:
/// no layer's tracklets start there, and nothing reaches them in time, so
/// they only make the timeframe's rows -- the used mask, the rounds' slots
/// -- larger than anything a slice needs.
std::vector<Hit> hitsWithLonelyRows()
{
  constexpr std::size_t kLonely = 30000U;
  auto result = hits();
  result.reserve(result.size() + kLonely);
  for (std::size_t i = 0; i < kLonely; ++i) {
    const auto k = static_cast<double>(i);
    result.push_back({
      .point = helixAt(std::fmod(0.37 * k, 2 * std::numbers::pi_v<double>), 0.02 * std::sin(k), 5., static_cast<LayerID>(kLayers - 1U)),
      .timeNS = 1000. + (99. * k / static_cast<double>(kLonely)),
    });
  }
  return result;
}

} // namespace

/// An iteration whose end the arena refuses -- in the gather, the
/// selection, the refit's or the extension's storage, the track chunk --
/// fails and leaves the used mask and the tracks as they were:
/// nothing kept, and the next iteration
/// runs on a free mask, as a host run of that iteration alone does.
namespace
{

/// Two collisions of eight tracks each, 500 ns apart: one vertex each for
/// the vertexing pass.
constexpr std::size_t kTracksPerCollision = 8U;
constexpr std::array<double, 2> kCollisionsNS{37., 537.};
constexpr double kCollisionDurationNS = 600.;

std::vector<Hit> collisionHits()
{
  std::vector<Hit> result;
  for (std::size_t track = 0; track < kTracksPerCollision * kCollisionsNS.size(); ++track) {
    const auto k = static_cast<double>(track);
    const double pt = 0.5 + (0.3 * static_cast<double>(track % kTracksPerCollision));
    const auto one = sensorHelixHits(std::fmod(0.3 + (2.1 * k), 2 * std::numbers::pi_v<double>), 0.03 * std::cos(1.3 * k), pt, kCollisionsNS[track / kTracksPerCollision]);
    result.insert(result.end(), one.begin(), one.end());
  }
  return result;
}

/// The pass: the beam line, the three innermost layers, no holes, 0.3 GeV.
IterationSettings passSettings()
{
  auto settings = iterationSettings();
  settings.windows.minPtGeV = 0.3;
  settings.minClusters = 3;
  settings.maxKinks = 0;
  settings.skipLinks = false;
  settings.layerMask = 0b111U;
  settings.follow = false;
  return settings;
}

} // namespace

/// The device's rerun of the frames without a vertex: the two collisions
/// and, in a frame of its own, two tracks, a vertex only at the lower floor;
/// the device's list the host's.
BOOST_AUTO_TEST_CASE(TheVertexerRerunsAFrameWithoutAVertexAsTheHost)
{
  auto* device = runtime();
  const auto tracker = kernels().makeToyBarrel8Tracker(*device, block(kArenaBytes));
  auto hits = collisionHits();
  for (const double pt : {0.8, 1.1}) {
    const auto one = sensorHelixHits(pt, 0.02 * pt, pt, 287.);
    hits.insert(hits.end(), one.begin(), one.end());
  }
  const Fixture fixture{hits, kCollisionDurationNS};
  const std::array<double, kLayers> material{};
  const auto scattering = scatteringTable<Det>(material);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const VertexerSettings settings{.aloneMinContributors = 2U};

  VertexingPass<Det, UniformMaterial> pass{passSettings(), settings, Fixture::kBeamConstraint, scattering, pool, 2U, 2};
  const auto expected = pass.find(*fixture.frame, 13, 0);
  BOOST_REQUIRE(expected.status == StageStatus::kComplete);
  BOOST_REQUIRE(expected.found.vertices.size() == kCollisionsNS.size() + 1U);

  const auto host = makeHost({.source = VertexSource::kVertexList});
  auto parameters = parametersOf(*host);
  parameters.push_back(pass.passParameters());
  std::vector<SearchWindowParams> windows;
  windows.reserve(parameters.size());
  for (std::size_t iteration = 0; iteration < host->iterationCount(); ++iteration) {
    windows.push_back(iterations(VertexSource::kVertexList)[iteration].windows);
  }
  windows.push_back(passSettings().windows);
  BOOST_REQUIRE((tracker->configureRun(parameters, trackingGeometry(), UniformMaterial{}, Propagator{}) == StageStatus::kComplete));
  const auto empty = planTimeframe<Det>(*fixture.frame, parameters, windows, host->beamLine(), {}, 13, 0);
  BOOST_REQUIRE((tracker->beginTimeframe(*fixture.frame, empty) == StageStatus::kComplete));
  VertexerResult found;
  BOOST_REQUIRE((tracker->findVertices(host->iterationCount(), settings, Fixture::kBeamConstraint, found) == StageStatus::kComplete));
  BOOST_REQUIRE(found.vertices.size() == expected.found.vertices.size());
  BOOST_TEST(found.rounds == expected.found.rounds);
  for (std::size_t v = 0; v < found.vertices.size(); ++v) {
    BOOST_TEST_CONTEXT("vertex " << v)
    {
      BOOST_TEST(found.vertices[v].timePS == expected.found.vertices[v].timePS);
      BOOST_TEST(found.vertices[v].widthPS == expected.found.vertices[v].widthPS);
      BOOST_TEST(found.contributors[v] == expected.found.contributors[v]);
      BOOST_TEST(std::abs(found.vertices[v].z - expected.found.vertices[v].z) < 1.e-5F);
    }
  }
}

/// The vertexing pass on the device, configured after the tracking's
/// iterations: the host's vertex list from it, and with the list -- the plan
/// replaced, the measurements kept -- the host's tracks; out of order it is
/// refused.
BOOST_AUTO_TEST_CASE(TheVertexingPassGivesTheHostsListAndTracks)
{
  auto* device = runtime();
  const auto tracker = kernels().makeToyBarrel8Tracker(*device, block(kArenaBytes));
  const Fixture fixture{collisionHits(), kCollisionDurationNS};
  const std::array<double, kLayers> material{};
  const auto scattering = scatteringTable<Det>(material);
  const auto pool = std::make_shared<BoundedMemoryResource>();

  // The host: the pass's list, and the tracking from it.
  VertexingPass<Det, UniformMaterial> pass{passSettings(), VertexerSettings{}, Fixture::kBeamConstraint, scattering, pool, 2U, 2};
  const auto expected = pass.find(*fixture.frame, 13, 0, true);
  BOOST_REQUIRE(expected.status == StageStatus::kComplete);
  BOOST_REQUIRE(expected.found.vertices.size() == kCollisionsNS.size());
  const auto host = makeHost({.source = VertexSource::kVertexList});
  const auto hostPlan = host->planTimeframe(*fixture.frame, expected.found.vertices, 13, 0.);
  const auto hostRun = host->run(*fixture.frame, hostPlan);
  BOOST_REQUIRE(hostRun.status == StageStatus::kComplete);
  BOOST_REQUIRE(hostRun.tracks->size() == kTracksPerCollision * kCollisionsNS.size());

  // The device: the tracking's iterations, then the pass.
  auto parameters = parametersOf(*host);
  parameters.push_back(pass.passParameters());
  std::vector<SearchWindowParams> windows;
  windows.reserve(parameters.size());
  for (std::size_t iteration = 0; iteration < host->iterationCount(); ++iteration) {
    windows.push_back(iterations(VertexSource::kVertexList)[iteration].windows);
  }
  windows.push_back(passSettings().windows);
  BOOST_REQUIRE((tracker->configureRun(parameters, trackingGeometry(), UniformMaterial{}, Propagator{}) == StageStatus::kComplete));
  const auto empty = planTimeframe<Det>(*fixture.frame, parameters, windows, host->beamLine(), {}, 13, 0);
  BOOST_REQUIRE((tracker->beginTimeframe(*fixture.frame, empty) == StageStatus::kComplete));
  VertexerResult found;
  std::vector<TripletClaim> claims;
  BOOST_REQUIRE((tracker->findVertices(host->iterationCount(), VertexerSettings{}, Fixture::kBeamConstraint, found, &claims) == StageStatus::kComplete));
  BOOST_REQUIRE(found.vertices.size() == expected.found.vertices.size());
  BOOST_TEST(found.candidates == expected.found.candidates);
  for (std::size_t v = 0; v < found.vertices.size(); ++v) {
    BOOST_TEST_CONTEXT("vertex " << v)
    {
      BOOST_TEST(found.vertices[v].timePS == expected.found.vertices[v].timePS);
      BOOST_TEST(found.vertices[v].widthPS == expected.found.vertices[v].widthPS);
      BOOST_TEST(found.contributors[v] == expected.found.contributors[v]);
      BOOST_TEST(std::abs(found.vertices[v].z - expected.found.vertices[v].z) < 1.e-5F);
    }
  }

  // Contributor identities must survive both vertex sorting and slice offsets.
  auto expectedClaims = expected.claims;
  std::erase_if(expectedClaims, [](const TripletClaim& claim) { return claim.vertex < 0; });
  const auto order = [](const TripletClaim& claim) { return std::pair{claim.vertex, claim.rows}; };
  std::ranges::sort(expectedClaims, {}, order);
  std::ranges::sort(claims, {}, order);
  BOOST_REQUIRE(claims.size() == expectedClaims.size());
  for (std::size_t c = 0; c < claims.size(); ++c) {
    BOOST_TEST(claims[c].vertex == expectedClaims[c].vertex);
    BOOST_TEST(claims[c].rows == expectedClaims[c].rows, boost::test_tools::per_element());
  }

  // The tracking's plan with the device's list replaces the empty one.
  const auto trackingPlan = planTimeframe<Det>(*fixture.frame, parameters, windows, host->beamLine(), found.vertices, 13, 0);
  BOOST_REQUIRE((tracker->replacePlan(trackingPlan) == StageStatus::kComplete));
  for (std::size_t iteration = 0; iteration < host->iterationCount(); ++iteration) {
    tracker->beginIteration(iteration);
    for (std::size_t slice = 0; slice < trackingPlan.slices.size(); ++slice) {
      BOOST_REQUIRE((tracker->runSlice(slice, nullptr) == StageStatus::kComplete));
    }
    BOOST_REQUIRE((tracker->endIteration(nullptr) == StageStatus::kComplete));
  }
  typename TrackRefit<Det>::Tracks tracks;
  tracks.setMemoryResource(pool.get());
  tracker->downloadTracks(tracks);
  BOOST_TEST(tracks.size() == hostRun.tracks->size());
  BOOST_TEST(sameTracks(std::as_const(*hostRun.tracks).view(), hostRun.tracks->size(), std::as_const(tracks).view(), tracks.size()));

  // Out of order, or with another run's plan: refused.
  tracker->beginIteration(0U);
  BOOST_CHECK_THROW(static_cast<void>(tracker->findVertices(host->iterationCount(), VertexerSettings{}, Fixture::kBeamConstraint, found)), std::logic_error);
  BOOST_CHECK_THROW(static_cast<void>(tracker->replacePlan(trackingPlan)), std::logic_error);
  static_cast<void>(tracker->endIteration(nullptr));
  BOOST_CHECK_THROW(static_cast<void>(tracker->replacePlan(hostPlan)), std::invalid_argument);
  BOOST_CHECK_THROW(static_cast<void>(tracker->findVertices(host->iterationCount() + 1U, VertexerSettings{}, Fixture::kBeamConstraint, found)), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(AFailedEndLeavesTheMaskAndTheTracks)
{
  auto* device = runtime();
  const Fixture fixture{hitsWithLonelyRows(), kDurationNS};
  const auto host = makeHost({});
  const auto parameters = parametersOf(*host);
  const auto plan = host->planTimeframe(*fixture.frame, {}, 13, 0.);
  const auto pool = std::make_shared<BoundedMemoryResource>();

  // The second iteration alone on the host: its mask from a free one.
  auto alone = std::make_unique<Host>(std::vector<IterationSettings>{iterations(VertexSource::kBeamLine)[1]}, Fixture::kBeamConstraint, scatteringTable<Det>(std::array<double, kLayers>{}), std::make_shared<BoundedMemoryResource>(), 2U, 2);
  static const UniformMaterial none{};
  alone->setTrackModel(trackingGeometry(), none, Propagator{});
  std::vector<std::size_t> aloneTracklets;
  alone->setSliceObserver([&](const Host&, const Host::SliceObservation& slice) { aloneTracklets.push_back(slice.result->tracklets->size()); });
  const auto aloneResult = alone->run(*fixture.frame, alone->planTimeframe(*fixture.frame, {}, 13, 0.));
  BOOST_REQUIRE(aloneResult.status == StageStatus::kComplete);
  const std::vector<std::uint8_t> aloneMask{alone->used().begin(), alone->used().end()};

  // What the first iteration's slices need, and what its end needs besides.
  std::vector<std::size_t> tracklets;
  const auto runSlices = [&](DeviceTracker<Det, UniformMaterial>& tracker, std::size_t iteration) {
    tracklets.clear();
    tracker.beginIteration(iteration);
    for (std::size_t slice = 0; slice < plan.slices.size(); ++slice) {
      if (tracker.runSlice(slice, nullptr) != StageStatus::kComplete) {
        return false;
      }
      tracklets.push_back(tracker.lastSlice().tracklets);
    }
    return true;
  };
  std::size_t slicePeak = 0U;
  std::size_t endPeak = 0U;
  {
    const auto probe = kernels().makeToyBarrel8Tracker(*device, block(kArenaBytes));
    probe->setStagingCapacity(1U);
    BOOST_REQUIRE((probe->configureRun(parameters, trackingGeometry(), UniformMaterial{}, Propagator{}) == StageStatus::kComplete));
    BOOST_REQUIRE((probe->beginTimeframe(*fixture.frame, plan) == StageStatus::kComplete));
    BOOST_REQUIRE(runSlices(*probe, 0U));
    slicePeak = probe->arenaPeakBytes();
    BOOST_REQUIRE((probe->endIteration(nullptr) == StageStatus::kComplete));
    endPeak = probe->arenaPeakBytes();
    // After a committed first iteration the second's tracklets are fewer:
    // the check below tells a free mask from a committed one.
    BOOST_REQUIRE(runSlices(*probe, 1U));
    BOOST_TEST(tracklets != aloneTracklets);
    BOOST_REQUIRE((probe->endIteration(nullptr) == StageStatus::kComplete));
  }
  BOOST_REQUIRE(endPeak > slicePeak);
  constexpr std::size_t kAlignment = 256U;
  const std::size_t step = std::max(kAlignment, ((endPeak - slicePeak) / 48U) / kAlignment * kAlignment);
  std::size_t failedEnds = 0U;
  std::size_t checkedAfter = 0U;
  for (std::size_t bytes = slicePeak; bytes < endPeak; bytes += step) {
    const SliceScope slices;
    const auto tracker = kernels().makeToyBarrel8Tracker(*device, block(bytes));
    tracker->setStagingCapacity(1U);
    BOOST_REQUIRE((tracker->configureRun(parameters, trackingGeometry(), UniformMaterial{}, Propagator{}) == StageStatus::kComplete));
    BOOST_REQUIRE((tracker->beginTimeframe(*fixture.frame, plan) == StageStatus::kComplete));
    const auto scratchBefore = tracker->arenaUsedBytes(ArenaEnd::kScratch);
    const auto baseBefore = tracker->arenaUsedBytes(ArenaEnd::kBase);
    if (!runSlices(*tracker, 0U)) {
      static_cast<void>(tracker->endIteration(nullptr));
      continue;
    }
    const auto status = tracker->endIteration(nullptr);
    if (status == StageStatus::kComplete) {
      continue;
    }
    BOOST_TEST_CONTEXT("arena of " << bytes << " bytes")
    {
      BOOST_TEST((status == StageStatus::kOutOfMemory));
      ++failedEnds;
      BOOST_TEST(tracker->arenaUsedBytes(ArenaEnd::kScratch) == scratchBefore);
      BOOST_TEST(tracker->arenaUsedBytes(ArenaEnd::kBase) == baseBefore);
      typename TrackRefit<Det>::Tracks tracks;
      tracks.setMemoryResource(pool.get());
      tracker->downloadTracks(tracks);
      BOOST_TEST(tracks.size() == 0U);
      // The next iteration reads a free mask: its tracklets are those of the
      // iteration alone, which would lose the first iteration's rows to a
      // committed mask. Its end, if the arena has room, leaves the mask the
      // iteration alone does.
      if (runSlices(*tracker, 1U)) {
        BOOST_TEST(tracklets == aloneTracklets, boost::test_tools::per_element());
        ++checkedAfter;
        DeviceIterationDownload<Det> onDevice{pool};
        if (tracker->endIteration(&onDevice) == StageStatus::kComplete) {
          BOOST_TEST(onDevice.used == aloneMask, boost::test_tools::per_element());
          BOOST_TEST(onDevice.tracks.size() == aloneResult.tracks->size());
        }
      } else {
        static_cast<void>(tracker->endIteration(nullptr));
      }
    }
  }
  BOOST_TEST(failedEnds > 0U);
  BOOST_TEST(checkedAfter > 0U);
  BOOST_TEST_MESSAGE(device->backend() << ": ends need " << slicePeak << " to " << endPeak << " bytes; " << failedEnds << " arenas failed the end, " << checkedAfter << " then ran the next iteration");
}

} // namespace o2::gtrk
