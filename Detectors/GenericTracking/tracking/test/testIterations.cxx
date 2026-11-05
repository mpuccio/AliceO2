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

/// \file testIterations.cxx
/// \brief The timeframe's tracking, iteration after iteration over every
///        slice, the selection between them and the refit: each track found
///        once, by the iteration it belongs to, whatever the slicing.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Iterations
#include <boost/test/unit_test.hpp>
#include "CommonConstants/LHCConstants.h"
#include <numbers>

#include "GenericTracking/Study/Digest.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/SlicePlan.h"
#include "GenericTracking/Tracking/Tracker.h"

#include "ToyBarrel8Fixture.h"
#include "UniformMaterial.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

namespace o2::gtrk
{

namespace
{

using namespace toybarrel8;

/// Twelve tracks over 1.1 us: the even ones stiff, 2 to 3 GeV, the odd ones
/// soft, 0.4 to 0.7 GeV -- below the first iteration's 1 GeV, above the
/// second's 0.3.
constexpr std::size_t kTracks = 12U;
constexpr double kDurationNS = 1100.;

bool stiff(std::size_t track) { return track % 2U == 0U; }

std::vector<Hit> hits()
{
  std::vector<Hit> result;
  for (std::size_t track = 0; track < kTracks; ++track) {
    const auto k = static_cast<double>(track);
    const double pt = stiff(track) ? 2. + (0.2 * static_cast<double>(track % 5U)) : 0.4 + (0.06 * static_cast<double>(track % 5U));
    const auto one = sensorHelixHits(std::fmod(0.3 + (0.9 * k), 2 * std::numbers::pi_v<double>), 0.03 * std::cos(1.3 * k), pt, 37. + (83. * k));
    BOOST_REQUIRE(one.size() == kLayers);
    result.insert(result.end(), one.begin(), one.end());
  }
  return result;
}

std::vector<IterationSettings> iterations()
{
  auto first = iterationSettings();
  auto second = iterationSettings();
  second.windows.minPtGeV = 0.3;
  return {first, second};
}

/// A track as the test compares it: its iteration, its rows and its chi2's
/// bits.
struct Found {
  std::uint8_t iteration = 0U;
  std::vector<std::uint32_t> rows;
  std::uint32_t chi2 = 0U;
  bool operator==(const Found&) const = default;
};

struct Outcome {
  std::vector<Found> tracks;
  std::vector<std::uint8_t> used;
  std::vector<std::size_t> records;
};

auto rows(const std::vector<Found>& tracks)
{
  std::vector<std::vector<std::uint32_t>> result;
  result.reserve(tracks.size());
  for (const auto& track : tracks) {
    result.push_back(track.rows);
  }
  std::ranges::sort(result);
  return result;
}

std::unique_ptr<Tracker<Det, UniformMaterial>> makeTracker(std::vector<IterationSettings> settings, bool follow = true, std::uint32_t shared = 0U, int threads = 1)
{
  for (auto& iteration : settings) {
    iteration.follow = follow;
    iteration.sharedInnermost = shared;
  }
  const std::array<double, kLayers> material{};
  static const UniformMaterial none{};
  auto tracker = std::make_unique<Tracker<Det, UniformMaterial>>(settings, Fixture::kBeamConstraint, scatteringTable<Det>(material),
                                                                 std::make_shared<BoundedMemoryResource>(), 2U, threads);
  tracker->setTrackModel(trackingGeometry(), none, Propagator{});
  return tracker;
}

std::vector<Found> found(const Tracker<Det, UniformMaterial>::TimeframeResult& result)
{
  std::vector<Found> tracks;
  tracks.reserve(result.tracks->size());
  const auto view = result.tracks->view();
  for (std::size_t i = 0; i < result.tracks->size(); ++i) {
    tracks.push_back({.iteration = view.iteration[i],
                      .rows = {view.clusters[i].index, view.clusters[i].index + view.count[i]},
                      .chi2 = std::bit_cast<std::uint32_t>(view.chi2[i])});
  }
  return tracks;
}

Outcome track(const Fixture& fixture, std::int32_t widthBlocks, int threads)
{
  auto owned = makeTracker(iterations(), false, 0U, threads);
  auto& tracker = *owned;
  const auto result = tracker.run(*fixture.frame, tracker.planTimeframe(*fixture.frame, {}, widthBlocks, 0.));
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  BOOST_REQUIRE(result.tracks != nullptr);
  Outcome outcome;
  for (const auto& summary : result.iterations) {
    BOOST_TEST(summary.mismatches == 0U);
    BOOST_TEST(summary.tracksFollowed == 0U);
    BOOST_TEST(summary.clustersFollowedOutward == 0U);
    BOOST_TEST(summary.clustersFollowedInward == 0U);
    outcome.records.push_back(summary.records);
  }
  outcome.tracks = found(result);
  outcome.used.assign(tracker.used().begin(), tracker.used().end());
  return outcome;
}

} // namespace

BOOST_AUTO_TEST_CASE(EachIterationFindsItsTracksOnce)
{
  const Fixture fixture{hits(), kDurationNS};
  const auto whole = track(fixture, fixture.geometry.blockCount(), 1);
  BOOST_REQUIRE(whole.tracks.size() == kTracks);
  std::vector<bool> found(kTracks, false);
  for (const auto& t : whole.tracks) {
    BOOST_REQUIRE(t.rows.size() == kLayers);
    const auto index = fixture.frame->measurement(t.rows.front()).clusterIndex / kLayers;
    BOOST_TEST_CONTEXT("track " << index)
    {
      // The stiff tracks by the first iteration, the soft ones by the second,
      // each whole and once.
      BOOST_TEST(t.iteration == (stiff(index) ? 0U : 1U));
      BOOST_TEST(!found[index]);
      for (std::size_t layer = 0; layer < kLayers; ++layer) {
        BOOST_TEST(t.rows[layer] == fixture.row((kLayers * index) + layer));
      }
    }
    found[index] = true;
  }
  // The first iteration's rows are exclusive when the second runs: it fits
  // only the soft tracks' roads.
  BOOST_TEST(whole.records[0] == kTracks / 2U);
  BOOST_TEST(whole.records[1] == kTracks / 2U);
  BOOST_TEST(std::ranges::all_of(whole.used, [](std::uint8_t state) { return state == kExclusive; }));
}

/// A bounded pool ends a run with kOutOfMemory wherever it runs out -- never
/// with an exception out of `run` -- and leaves a tracker that runs again and
/// gives back every byte when it goes. Every
/// limit from a few bytes to a whole run's peak, in steps finer than the
/// stages' staging, so that each allocation of a run is the one refused at
/// some limit: the road fit's and the selection's setup included.
BOOST_AUTO_TEST_CASE(AnyMemoryLimitEndsARunWithAStatus)
{
  const Fixture fixture{hits(), kDurationNS};
  const std::array<double, kLayers> material{};
  const UniformMaterial none{};
  constexpr std::int32_t kWidth = 3;
  const auto runOn = [&](const std::shared_ptr<BoundedMemoryResource>& pool, int threads, std::vector<StageStatus>& statuses) {
    Tracker<Det, UniformMaterial> tracker{iterations(), Fixture::kBeamConstraint, scatteringTable<Det>(material), pool, 2U, threads};
    tracker.setTrackModel(trackingGeometry(), none, Propagator{});
    for (int run = 0; run < 2; ++run) {
      Tracker<Det, UniformMaterial>::TimeframeResult result;
      BOOST_CHECK_NO_THROW(result = tracker.run(*fixture.frame, tracker.planTimeframe(*fixture.frame, {}, kWidth, 0.)));
      statuses.push_back(result.status);
    }
  };
  std::size_t peak = 0U;
  {
    const auto unbounded = std::make_shared<BoundedMemoryResource>();
    std::vector<StageStatus> statuses;
    runOn(unbounded, 1, statuses);
    BOOST_REQUIRE((statuses == std::vector<StageStatus>{StageStatus::kComplete, StageStatus::kComplete}));
    BOOST_TEST(unbounded->getUsedMemory() == 0U);
    peak = unbounded->getPeakMemory();
  }
  BOOST_REQUIRE(peak > 0U);
  constexpr std::size_t kSteps = 1000U;
  std::size_t failed = 0U;
  for (std::size_t step = 0; step <= kSteps; ++step) {
    const std::size_t limit = 1U + (peak * step / kSteps);
    for (const int threads : {1, 3}) {
      const auto bounded = std::make_shared<BoundedMemoryResource>(limit);
      std::vector<StageStatus> statuses;
      runOn(bounded, threads, statuses);
      BOOST_TEST_CONTEXT("limit " << limit << " of " << peak << " bytes, " << threads << " threads")
      {
        for (const auto status : statuses) {
          BOOST_TEST((status == StageStatus::kComplete || status == StageStatus::kOutOfMemory));
          failed += status == StageStatus::kOutOfMemory ? 1U : 0U;
        }
        BOOST_TEST(bounded->getUsedMemory() == 0U);
      }
    }
  }
  BOOST_TEST(failed != 0U);
}

BOOST_AUTO_TEST_CASE(TheTracksDoNotDependOnTheSlicingOrTheThreads)
{
  const Fixture fixture{hits(), kDurationNS};
  const auto whole = track(fixture, fixture.geometry.blockCount(), 1);
  for (const std::int32_t width : {1, 2, 3, 7}) {
    for (const int threads : {1, 4}) {
      const auto sliced = track(fixture, width, threads);
      BOOST_TEST_CONTEXT("slices of " << width << " blocks, " << threads << " threads")
      {
        BOOST_TEST((sliced.tracks == whole.tracks));
        BOOST_TEST(sliced.used == whole.used, boost::test_tools::per_element());
      }
    }
  }
}

/// The observer sees every slice of every iteration once, completed, with the
/// used mask it read: all free in the first iteration, the first iteration's
/// tracks exclusive in the second. Observing changes nothing.
BOOST_AUTO_TEST_CASE(TheSliceObserverSeesEverySliceWithItsUsedMask)
{
  const Fixture fixture{hits(), kDurationNS};
  constexpr std::int32_t kWidth = 3;
  const auto unobserved = track(fixture, kWidth, 1);

  const std::array<double, kLayers> material{};
  Tracker<Det, UniformMaterial> tracker{
    iterations(), Fixture::kBeamConstraint, scatteringTable<Det>(material), std::make_shared<BoundedMemoryResource>(), 2U, 1};
  const UniformMaterial none{};
  tracker.setTrackModel(trackingGeometry(), none, Propagator{});
  std::vector<std::size_t> slicesPerIteration(2U, 0U);
  std::vector<std::size_t> exclusivePerIteration(2U, 0U);
  std::vector<std::size_t> recordsPerIteration(2U, 0U);
  std::int32_t nextBlock = 0;
  tracker.setSliceObserver(
    [&](const Tracker<Det, UniformMaterial>& observed, const Tracker<Det, UniformMaterial>::SliceObservation& slice) {
      BOOST_TEST(&observed == &tracker);
      BOOST_REQUIRE(slice.iteration < 2U);
      BOOST_TEST(slice.frame == fixture.frame.get());
      BOOST_TEST((slice.result->status == StageStatus::kComplete));
      BOOST_TEST(slice.result->tracklets != nullptr);
      // After the roads: the slice's road counts and records are there.
      BOOST_REQUIRE(slice.roads != nullptr);
      BOOST_REQUIRE(slice.records != nullptr);
      BOOST_TEST(slice.roads->roads().size() == slice.result->triplets->size());
      BOOST_TEST(slice.records->elementCount() == slice.roads->roadCount());
      recordsPerIteration[slice.iteration] += slice.records->size();
      BOOST_TEST(slice.slice->coreBlocks.begin == nextBlock);
      nextBlock = slice.slice->coreBlocks.end == fixture.geometry.blockCount() ? 0 : slice.slice->coreBlocks.end;
      ++slicesPerIteration[slice.iteration];
      exclusivePerIteration[slice.iteration] = static_cast<std::size_t>(std::ranges::count(observed.used(), std::uint8_t{kExclusive}));
    });
  const auto result = tracker.run(*fixture.frame, tracker.planTimeframe(*fixture.frame, {}, kWidth, 0.));
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  const auto slices = static_cast<std::size_t>((fixture.geometry.blockCount() + kWidth - 1) / kWidth);
  BOOST_TEST(slicesPerIteration[0] == slices);
  BOOST_TEST(slicesPerIteration[1] == slices);
  BOOST_TEST(exclusivePerIteration[0] == 0U);
  BOOST_TEST(exclusivePerIteration[1] == (kTracks / 2U) * kLayers);
  BOOST_TEST((found(result) == unobserved.tracks));
  BOOST_TEST(recordsPerIteration[0] > 0U);
  BOOST_TEST(recordsPerIteration[0] == result.iterations[0].records);
  BOOST_TEST(recordsPerIteration[1] == result.iterations[1].records);
}

/// The iteration observer sees every iteration once, completed: its
/// selection holds every slice's records in slice order, and what it
/// accepted; the used mask is the one the selection and follow left.
BOOST_AUTO_TEST_CASE(TheIterationObserverSeesEveryIterationsSelection)
{
  const Fixture fixture{hits(), kDurationNS};
  const std::array<double, kLayers> material{};
  Tracker<Det, UniformMaterial> tracker{
    iterations(), Fixture::kBeamConstraint, scatteringTable<Det>(material), std::make_shared<BoundedMemoryResource>(), 2U, 1};
  const UniformMaterial none{};
  tracker.setTrackModel(trackingGeometry(), none, Propagator{});
  std::vector<std::uint32_t> firstCluster;
  tracker.setSliceObserver([&](const Tracker<Det, UniformMaterial>&, const Tracker<Det, UniformMaterial>::SliceObservation& slice) {
    for (std::uint32_t record = 0; record < slice.records->size(); ++record) {
      firstCluster.push_back(slice.records->view().clusters[record].index[0]);
    }
  });
  std::vector<std::size_t> seen;
  tracker.setIterationObserver(
    [&](const Tracker<Det, UniformMaterial>& observed, const Tracker<Det, UniformMaterial>::IterationObservation& iteration) {
      seen.push_back(iteration.iteration);
      BOOST_REQUIRE(iteration.selection != nullptr && iteration.summary != nullptr);
      if (iteration.selection == nullptr || iteration.summary == nullptr) {
        return; // for the static analysis: BOOST_REQUIRE has already failed
      }
      BOOST_TEST(iteration.summary->records == iteration.selection->recordCount());
      BOOST_TEST(iteration.summary->tracks == iteration.selection->selected().size());
      // The appended records are the slices', in slice order.
      BOOST_REQUIRE(iteration.selection->recordCount() == firstCluster.size());
      for (std::uint32_t record = 0; record < iteration.selection->recordCount(); ++record) {
        BOOST_TEST(iteration.selection->records().clusters[record].index[0] == firstCluster[record]);
      }
      firstCluster.clear();
      const auto exclusive = std::ranges::count(observed.used(), std::uint8_t{kExclusive});
      BOOST_TEST(static_cast<std::size_t>(exclusive) == (iteration.iteration + 1U) * (kTracks / 2U) * kLayers);
    });
  const auto result = tracker.run(*fixture.frame, tracker.planTimeframe(*fixture.frame, {}, 3, 0.));
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  BOOST_TEST((seen == std::vector<std::size_t>{0U, 1U}));
}

namespace
{

std::vector<IterationSettings> iterations(VertexSource source)
{
  auto settings = iterations();
  for (auto& one : settings) {
    one.vertexSource = source;
  }
  return settings;
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

} // namespace

/// A plan run is the run that plans for itself, beam line and vertex list.
BOOST_AUTO_TEST_CASE(APlannedRunGivesThePlanningRunsTracks)
{
  const Fixture fixture{hits(), kDurationNS};
  const auto vertices = trackVertices();
  for (const auto source : {VertexSource::kBeamLine, VertexSource::kVertexList}) {
    BOOST_TEST_CONTEXT((source == VertexSource::kBeamLine ? "beam line" : "vertex list"))
    {
      const auto settings = iterations(source);
      auto planning = makeTracker(settings);
      const auto expected = planning->run(*fixture.frame, planning->planTimeframe(*fixture.frame, vertices, 3, 0.));
      BOOST_REQUIRE(expected.status == StageStatus::kComplete);
      auto planned = makeTracker(settings);
      const auto plan = planned->planTimeframe(*fixture.frame, vertices, 3, 0.);
      BOOST_TEST(plan.beamLine.varZ == planned->beamLine().varZ);
      BOOST_TEST(plan.vertices.size() == vertices.size());
      const auto result = planned->run(*fixture.frame, plan);
      BOOST_REQUIRE(result.status == StageStatus::kComplete);
      BOOST_TEST((found(result) == found(expected)));
      BOOST_TEST(found(result).size() == kTracks);
    }
  }
}

/// A run resumed from the checkpoint of its first iteration, by another tracker
/// of the same configuration, is the whole run: the same digest -- tracks,
/// used mask and summaries -- and tracks.
BOOST_AUTO_TEST_CASE(AResumedRunIsTheWholeRun)
{
  const Fixture fixture{hits(), kDurationNS};
  const auto settings = iterations();
  auto whole = makeTracker(settings);
  const auto plan = whole->planTimeframe(*fixture.frame, {}, 3, 0.);
  const auto expected = whole->run(*fixture.frame, plan);
  BOOST_REQUIRE(expected.status == StageStatus::kComplete);
  const auto expectedDigest = study::trackingDigest(expected, whole->used());
  // A plan is its tracker's: each plans its own.
  auto first = makeTracker(settings);
  const auto part = first->run(*fixture.frame, first->planTimeframe(*fixture.frame, {}, 3, 0.), 1U);
  BOOST_REQUIRE(part.status == StageStatus::kComplete);
  BOOST_TEST(part.iterations.size() == 1U);
  const auto checkpoint = first->checkpoint(part);
  BOOST_TEST(checkpoint.tracks.size() == kTracks / 2U);
  auto resumed = makeTracker(settings);
  const auto result = resumed->resume(*fixture.frame, resumed->planTimeframe(*fixture.frame, {}, 3, 0.), checkpoint);
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  BOOST_TEST(study::trackingDigest(result, resumed->used()) == expectedDigest);
  BOOST_TEST((found(result) == found(expected)));
}

/// The plan holds the frame's bounds and one slice list for every iteration:
/// the halo, the gate around each core, the source maps; and each
/// iteration's vertex windows.
BOOST_AUTO_TEST_CASE(ASlicePlanHoldsTheSlicesBoundsAndGate)
{
  const Fixture fixture{hits(), kDurationNS};
  const auto tracker = makeTracker(iterations());
  const auto plan = tracker->planTimeframe(*fixture.frame, {}, 3, 0.);
  const auto& frame = *fixture.frame;
  BOOST_TEST(plan.slices.size() == static_cast<std::size_t>((fixture.geometry.blockCount() + 2) / 3));
  // The toys' blocks start at 0 and are whole picoseconds.
  const auto reach = sliceReach(frameBounds(frame), fixture.geometry).gatePS();
  const auto width = fixture.geometry.blockWidth();
  BOOST_REQUIRE(fixture.geometry.partsPerPS() == 1);
  for (const auto& slicePlan : plan.slices) {
    const auto& slice = slicePlan.slice;
    BOOST_TEST(slicePlan.gateEarliestPS == (static_cast<std::int64_t>(slice.coreBlocks.begin) * width) - reach);
    BOOST_TEST(slicePlan.gateLatestPS == (static_cast<std::int64_t>(slice.coreBlocks.end) * width) + reach);
    BOOST_TEST(slicePlan.sources.size() == TimeFrameScratch<Det>::makeSources(slice.halo).size());
    const auto inputs = slicePlan.inputs(tracker->iteration(0U), plan.layerMaxWidthPS);
    BOOST_TEST(inputs.gate.core == slice.core.data());
    BOOST_TEST(inputs.gate.earliestPS == slicePlan.gateEarliestPS);
  }
}

/// A plan of another frame, another number of iterations or other
/// parameters is refused before anything runs.
BOOST_AUTO_TEST_CASE(APlanOfAnotherFrameOrConfigurationIsRefused)
{
  const Fixture fixture{hits(), kDurationNS};
  const Fixture longer{hits(), 2. * kDurationNS};
  const auto tracker = makeTracker(iterations());
  const auto plan = tracker->planTimeframe(*fixture.frame, {}, 3, 0.);
  BOOST_CHECK_THROW(static_cast<void>(tracker->run(*longer.frame, plan)), std::invalid_argument);

  const auto one = makeTracker({iterations().front()});
  BOOST_CHECK_THROW(static_cast<void>(one->run(*fixture.frame, plan)), std::invalid_argument);

  auto other = iterations();
  other[1].windows.minPtGeV = 0.4;
  const auto otherTracker = makeTracker(other);
  BOOST_CHECK_THROW(static_cast<void>(otherTracker->run(*fixture.frame, plan)), std::invalid_argument);
  BOOST_CHECK_THROW(plan.check(*fixture.frame, std::span<const IterationParameters<Det>>{}), std::invalid_argument);
  BOOST_CHECK_NO_THROW(plan.check(*fixture.frame, plan.parameters));
}

/// The iteration observation carries the mask after the iteration, the
/// tracks with the iteration's offset, and the follow's proposals -- none
/// without the follow, never an earlier iteration's.
BOOST_AUTO_TEST_CASE(TheIterationObservationCarriesMaskTracksAndProposals)
{
  const Fixture fixture{hits(), kDurationNS};
  for (const bool follow : {true, false}) {
    BOOST_TEST_CONTEXT("follow " << follow)
    {
      auto tracker = makeTracker(iterations(), follow);
      std::vector<std::uint32_t> offsets;
      std::vector<std::size_t> proposals;
      tracker->setIterationObserver(
        [&](const Tracker<Det, UniformMaterial>& observed, const Tracker<Det, UniformMaterial>::IterationObservation& iteration) {
          BOOST_TEST(iteration.used.data() == observed.used().data());
          BOOST_REQUIRE(iteration.tracks != nullptr);
          BOOST_TEST(iteration.tracks->size() == iteration.trackOffset + iteration.summary->tracks);
          offsets.push_back(iteration.trackOffset);
          proposals.push_back(iteration.proposals.size());
        });
      const auto result = tracker->run(*fixture.frame, tracker->planTimeframe(*fixture.frame, {}, 3, 0.));
      BOOST_REQUIRE(result.status == StageStatus::kComplete);
      BOOST_TEST((offsets == std::vector<std::uint32_t>{0U, static_cast<std::uint32_t>(kTracks / 2U)}));
      const std::vector<std::size_t> expected = follow ? std::vector<std::size_t>{kTracks / 2U, kTracks / 2U} : std::vector<std::size_t>{0U, 0U};
      BOOST_TEST(proposals == expected, boost::test_tools::per_element());
    }
  }
}

/// Seeds on the inner layers, the follow over every layer:
/// the iterations' roads join layers
/// 0-4, their follow every layer, and each track comes out whole. A track
/// shorter than an iteration's length after the follow is dropped, and its
/// rows go back to the next iteration.
BOOST_AUTO_TEST_CASE(InnerSeedsFollowOverEveryLayerAndShortTracksGiveTheirRowsBack)
{
  const Fixture fixture{hits(), kDurationNS};
  const auto everyLayerTracker = makeTracker(iterations(), false);
  const auto everyLayer = everyLayerTracker->run(*fixture.frame, everyLayerTracker->planTimeframe(*fixture.frame, {}, 3, 0.));
  BOOST_REQUIRE(everyLayer.status == StageStatus::kComplete);
  const auto expected = found(everyLayer);
  BOOST_REQUIRE(expected.size() == kTracks);

  constexpr int kWhole = static_cast<int>(kLayers);
  const auto seeded = [](std::array<int, 2> minClustersAfterFollow) {
    auto settings = iterations();
    for (std::size_t i = 0; i < settings.size(); ++i) {
      settings[i].layerMask = 0b11111U;
      settings[i].skipLinks = false;
      settings[i].minClustersAfterFollow = minClustersAfterFollow[i];
    }
    return settings;
  };
  {
    auto tracker = makeTracker(seeded({0, 0}));
    const auto result = tracker->run(*fixture.frame, tracker->planTimeframe(*fixture.frame, {}, 3, 0.));
    BOOST_REQUIRE(result.status == StageStatus::kComplete);
    const auto tracks = found(result);
    BOOST_TEST(tracks.size() == kTracks);
    BOOST_TEST((rows(tracks) == rows(expected)));
    BOOST_TEST(result.iterations[0].tracksFollowed == kTracks / 2U);
    BOOST_TEST(result.iterations[0].clustersFollowedOutward == 3U * kTracks / 2U);
    BOOST_TEST(result.iterations[0].tracksDropped == 0U);
    BOOST_TEST(std::ranges::all_of(tracker->used(), [](std::uint8_t state) { return state == kExclusive; }));
  }
  {
    // The first iteration asks more than a track can hold: its six stiff
    // tracks are dropped, and the second, at 0.3 GeV, finds all twelve.
    auto tracker = makeTracker(seeded({kWhole + 1, kWhole}));
    std::vector<std::size_t> observed;
    tracker->setIterationObserver(
      [&](const Tracker<Det, UniformMaterial>&, const Tracker<Det, UniformMaterial>::IterationObservation& iteration) {
        BOOST_TEST(iteration.tracks->size() == iteration.trackOffset + iteration.summary->tracks - iteration.summary->tracksDropped);
        BOOST_TEST(iteration.proposals.size() == iteration.tracks->size() - iteration.trackOffset);
        observed.push_back(iteration.tracks->size());
      });
    const auto result = tracker->run(*fixture.frame, tracker->planTimeframe(*fixture.frame, {}, 3, 0.));
    BOOST_REQUIRE(result.status == StageStatus::kComplete);
    BOOST_TEST(result.iterations[0].tracksDropped == kTracks / 2U);
    BOOST_TEST(result.iterations[1].tracksDropped == 0U);
    BOOST_TEST((observed == std::vector<std::size_t>{0U, kTracks}));
    const auto tracks = found(result);
    BOOST_TEST((rows(tracks) == rows(expected)));
    BOOST_TEST(std::ranges::all_of(tracks, [](const Found& track) { return track.iteration == 1U; }));
    BOOST_TEST(std::ranges::all_of(tracker->used(), [](std::uint8_t state) { return state == kExclusive; }));
  }
  {
    // Both ask too much: nothing is kept, and every row is free again.
    auto tracker = makeTracker(seeded({kWhole + 1, kWhole + 1}));
    const auto result = tracker->run(*fixture.frame, tracker->planTimeframe(*fixture.frame, {}, 3, 0.));
    BOOST_REQUIRE(result.status == StageStatus::kComplete);
    BOOST_TEST(result.tracks->size() == 0U);
    BOOST_TEST(std::ranges::all_of(tracker->used(), [](std::uint8_t state) { return state == kFree; }));
  }
}

/// Seeds on the outer layers, updated with the beam line, the follow over
/// every layer inward from the constrained state: each track comes out whole, as the iterations over every layer
/// find it.
BOOST_AUTO_TEST_CASE(OuterSeedsOnTheBeamLineFollowInwardWhole)
{
  const Fixture fixture{hits(), kDurationNS};
  const auto everyLayerTracker = makeTracker(iterations(), false);
  const auto everyLayer = everyLayerTracker->run(*fixture.frame, everyLayerTracker->planTimeframe(*fixture.frame, {}, 3, 0.));
  BOOST_REQUIRE(everyLayer.status == StageStatus::kComplete);

  const auto expected = rows(found(everyLayer));
  BOOST_REQUIRE(expected.size() == kTracks);
  auto settings = iterations();
  for (auto& iteration : settings) {
    iteration.layerMask = ((1U << kLayers) - 1U) & ~0b111U;
    iteration.skipLinks = false;
    iteration.beamConstraint = true;
  }
  auto tracker = makeTracker(settings);
  const auto result = tracker->run(*fixture.frame, tracker->planTimeframe(*fixture.frame, {}, 3, 0.));
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  BOOST_TEST((rows(found(result)) == expected));
  BOOST_TEST(result.iterations[0].clustersFollowedInward == 3U * kTracks / 2U);
  BOOST_TEST(std::ranges::all_of(tracker->used(), [](std::uint8_t state) { return state == kExclusive; }));
}

/// The recovery pass: a second iteration with the
/// kink allowance, following the first's tracks again (followEarlier), and a
/// length after the follow for its own tracks. On tracks without kinks it
/// finds what the first iteration did, every row held once: the earlier
/// tracks' follow takes nothing more and gives nothing back.
BOOST_AUTO_TEST_CASE(TheRecoveryPassKeepsWholeTracksWhole)
{
  const Fixture fixture{hits(), kDurationNS};

  const auto everyLayerTracker = makeTracker(iterations(), false);
  const auto expected = rows(found(everyLayerTracker->run(*fixture.frame, everyLayerTracker->planTimeframe(*fixture.frame, {}, 3, 0.))));
  BOOST_REQUIRE(expected.size() == kTracks);
  auto settings = iterations();
  for (auto& iteration : settings) {
    iteration.layerMask = ((1U << kLayers) - 1U) & ~0b111U;
    iteration.skipLinks = false;
    iteration.beamConstraint = true;
  }
  settings[1].kinkSigma = 5.;
  settings[1].maxKinks = 2;
  settings[1].followEarlier = true;
  settings[1].minClustersAfterFollow = static_cast<int>(kLayers);
  auto tracker = makeTracker(settings);
  const auto result = tracker->run(*fixture.frame, tracker->planTimeframe(*fixture.frame, {}, 3, 0.));
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  BOOST_TEST((rows(found(result)) == expected));
  BOOST_TEST(result.iterations[1].tracksDropped == 0U);
  BOOST_TEST(std::ranges::all_of(tracker->used(), [](std::uint8_t state) { return state == kExclusive; }));
}

/// The recovery pass on one graph (one graph for both rounds): the second
/// iteration takes the first one's graph, built with its
/// kink allowance, and finds the tracks the two iterations with a graph each
/// find -- their iterations, rows and chi2s -- building no slice of its own.
/// A second iteration whose graph settings differ, or without an allowance,
/// is refused.
BOOST_AUTO_TEST_CASE(TheRecoveryPassOnOneGraphFindsTheTwoGraphsTracks)
{
  const Fixture fixture{hits(), kDurationNS};
  auto settings = iterations();
  for (auto& iteration : settings) {
    iteration.windows = settings[0].windows;
    iteration.layerMask = ((1U << kLayers) - 1U) & ~0b111U;
    iteration.skipLinks = false;
    iteration.beamConstraint = true;
  }
  settings[1].kinkSigma = 5.;
  settings[1].maxKinks = 2;
  settings[1].followEarlier = true;
  settings[1].minClustersAfterFollow = static_cast<int>(kLayers);
  auto separate = makeTracker(settings);
  const auto twoGraphs = separate->run(*fixture.frame, separate->planTimeframe(*fixture.frame, {}, 3, 0.));
  BOOST_REQUIRE(twoGraphs.status == StageStatus::kComplete);
  settings[1].sharedGraph = true;
  auto tracker = makeTracker(settings);
  const auto oneGraph = tracker->run(*fixture.frame, tracker->planTimeframe(*fixture.frame, {}, 3, 0.));
  BOOST_REQUIRE(oneGraph.status == StageStatus::kComplete);
  // Both at the first iteration's floor: the stiffer half of the tracks.
  BOOST_TEST(found(oneGraph).size() == kTracks / 2U);
  BOOST_TEST((found(oneGraph) == found(twoGraphs)));
  BOOST_TEST(
    (std::vector<std::uint8_t>(tracker->used().begin(), tracker->used().end()) == std::vector<std::uint8_t>(separate->used().begin(), separate->used().end())));
  BOOST_TEST(oneGraph.iterations[1].slices == 0U);
  BOOST_TEST(oneGraph.iterations[1].roads == 0U);

  auto otherFloor = settings;
  otherFloor[1].windows.minPtGeV = 0.3;
  BOOST_CHECK_THROW(makeTracker(otherFloor), std::invalid_argument);
  auto noAllowance = settings;
  noAllowance[1].kinkSigma = 1.;
  BOOST_CHECK_THROW(makeTracker(noAllowance), std::invalid_argument);
  auto first = settings;
  first[0].sharedGraph = true;
  BOOST_CHECK_THROW(makeTracker(first), std::invalid_argument);
}

namespace
{

/// The toy with ALICE 3's two frame lengths:
/// layers 0-2 read out every 22 BC, layers 3-7 every 81 BC, both from BC 0,
/// widened by \a marginBC. The blocks are two of the longer frames, as
/// GridConfig.h's `gridParams` makes them.
constexpr std::int32_t kInnerFrameBC = 22;
constexpr std::int32_t kOuterFrameBC = 81;
constexpr std::int64_t kPartsPerPS = kPSPerBCDenominator;
constexpr std::int64_t kPartsPerBC = kPSPerBCNumerator;
/// A BC [ns], O2's LHCBunchSpacingNS.
constexpr double kBunchNS = o2::constants::lhc::LHCBunchSpacingNS;
/// The tracks' crossings. The inner frame of those at 160, 170, 326 and 650
/// BC straddles an outer frame's edge; for the last three the edge is a
/// block's (162, 324, 648), so their inner and outer rows start in
/// neighbouring blocks.
constexpr std::array<std::int64_t, 12> kCrossingBC{5, 40, 100, 160, 170, 250, 326, 400, 480, 520, 650, 700};
constexpr double kFramedDurationNS = 820. * kBunchNS;

std::int32_t frameLengthBC(LayerID layer) { return layer < 3 ? kInnerFrameBC : kOuterFrameBC; }

GridParams twoLengthGrid()
{
  GridParams params;
  params.partsPerPS = kPartsPerPS;
  params.blockWidth = kPartsPerBC * 2 * kOuterFrameBC;
  params.blockOrigin = -kPartsPerPS;
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    params.frames[layer] = {.lengthBC = frameLengthBC(static_cast<LayerID>(layer)), .firstStartBC = 0};
  }
  return params;
}

/// The box of the frame \a hit falls in on its layer, widened by \a marginBC.
TimeEstPS frameOf(const Hit& hit, std::int64_t marginBC)
{
  const auto bc = static_cast<std::int64_t>(std::floor(hit.timeNS / kBunchNS));
  const auto length = frameLengthBC(hit.point.layer);
  const auto start = (bc / length) * length;
  return frameBox(start, start + length, marginBC);
}

std::vector<Hit> framedHits()
{
  std::vector<Hit> result;
  for (std::size_t track = 0; track < kCrossingBC.size(); ++track) {
    const auto k = static_cast<double>(track);
    const double pt = stiff(track) ? 2. + (0.2 * static_cast<double>(track % 5U)) : 0.4 + (0.06 * static_cast<double>(track % 5U));
    const double startNS = (static_cast<double>(kCrossingBC[track]) * kBunchNS) + 5.;
    const auto one = sensorHelixHits(std::fmod(0.3 + (0.9 * k), 2 * std::numbers::pi_v<double>), 0.03 * std::cos(1.3 * k), pt, startNS);
    BOOST_REQUIRE(one.size() == kLayers);
    result.insert(result.end(), one.begin(), one.end());
  }
  return result;
}

struct FramedRun {
  std::vector<Found> tracks;
  std::vector<TimeEstPS> times; ///< by track index, from the clusters
  std::int64_t haloPS = 0;
};

FramedRun trackFramed(const Fixture& fixture, std::int32_t widthBlocks)
{
  auto tracker = makeTracker(iterations(), false);
  const auto result = tracker->run(*fixture.frame, tracker->planTimeframe(*fixture.frame, {}, widthBlocks, 0.));
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  FramedRun run{.tracks = found(result), .times = std::vector<TimeEstPS>(kCrossingBC.size()), .haloPS = result.iterations.front().haloPS};
  const auto view = result.tracks->view();
  for (std::size_t i = 0; i < result.tracks->size(); ++i) {
    const auto index = fixture.frame->measurement(view.clusters[i].index[0]).clusterIndex / kLayers;
    run.times[index] = {.start = view.timePS[i], .width = view.widthPS[i]};
  }
  return run;
}

} // namespace

/// Layers of two frame lengths, no margin: only the inner -> outer link
/// reaches, an inner frame forward and an outer one backward, so the gate is
/// one of each and the halo one outer frame more; the slices give the
/// one-slice tracks,
/// those whose inner and outer rows start on two sides of a slice's edge
/// among them; a track's time is its frames' intersection, 8 BC where the
/// inner frame straddles an outer edge.
BOOST_AUTO_TEST_CASE(TwoFrameLengthsSliceWithOneOuterFrameOfHalo)
{
  const auto hits = framedHits();
  const Fixture fixture{hits, kFramedDurationNS, twoLengthGrid(), [](const Hit& hit) { return frameOf(hit, 0); }};
  const auto bounds = frameBounds(*fixture.frame);
  BOOST_TEST(std::ranges::all_of(bounds.rowsAreFrames, [](bool framed) { return framed; }));

  // The crossing at 170 BC: inner frame [154, 176) in block 0, outer frame
  // [162, 243) in block 1.
  const auto& geometry = fixture.geometry;
  BOOST_TEST(geometry.blockOfPS(frameOf(hits[Fixture::logicalIndex(4U, 2U)], 0).start) == 0);
  BOOST_TEST(geometry.blockOfPS(frameOf(hits[Fixture::logicalIndex(4U, 3U)], 0).start) == 1);

  const auto whole = trackFramed(fixture, geometry.blockCount());
  BOOST_REQUIRE(whole.tracks.size() == kCrossingBC.size());
  for (std::size_t track = 0; track < kCrossingBC.size(); ++track) {
    BOOST_TEST_CONTEXT("track at " << kCrossingBC[track] << " BC")
    {
      const auto bc = kCrossingBC[track];
      const auto inner = (bc / kInnerFrameBC) * kInnerFrameBC;
      const auto outer = (bc / kOuterFrameBC) * kOuterFrameBC;
      const auto lowBC = std::max(inner, outer);
      const auto highBC = std::min(inner + kInnerFrameBC, outer + kOuterFrameBC);
      BOOST_TEST(whole.times[track].start == psOfBC(lowBC));
      // The ends are each rounded down from BC: within a picosecond of the BC interval's.
      BOOST_TEST(std::abs(static_cast<std::int64_t>(whole.times[track].upper()) - psOfBC(highBC)) <= 1);
    }
  }
  // At 160 BC: inner [154, 176), outer [81, 162), the time [154, 162).
  BOOST_TEST(std::abs(static_cast<std::int64_t>(whole.times[3].width) - psOfBC(8)) <= 1);
  for (const std::int32_t width : {1, 2, 3}) {
    BOOST_TEST_CONTEXT("slices of " << width << " blocks")
    {
      const auto sliced = trackFramed(fixture, width);
      BOOST_TEST((sliced.tracks == whole.tracks));
    }
  }
}

/// The same toy with its boxes widened by a margin while the grid still
/// declares the frames -- as a dump keeps the boxes its producer widened,
/// whatever margin the tracking is configured with: the rows are not their
/// frames, every link
/// reaches the wider box, and the slices still give the one-slice tracks.
BOOST_AUTO_TEST_CASE(TwoFrameLengthsWithAMarginReachEveryLink)
{
  constexpr std::int64_t kMarginBC = 4;
  const Fixture fixture{framedHits(), kFramedDurationNS, twoLengthGrid(), [](const Hit& hit) { return frameOf(hit, kMarginBC); }};
  const auto bounds = frameBounds(*fixture.frame);
  BOOST_TEST(std::ranges::none_of(bounds.rowsAreFrames, [](bool framed) { return framed; }));

  const auto whole = trackFramed(fixture, fixture.geometry.blockCount());
  BOOST_REQUIRE(whole.tracks.size() == kCrossingBC.size());
  for (const std::int32_t width : {1, 2, 3}) {
    BOOST_TEST_CONTEXT("slices of " << width << " blocks") { BOOST_TEST((trackFramed(fixture, width).tracks == whole.tracks)); }
  }
}

} // namespace o2::gtrk
