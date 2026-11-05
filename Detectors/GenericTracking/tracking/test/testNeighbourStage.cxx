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

/// \file testNeighbourStage.cxx
/// \brief The pulling neighbour stage against a brute-force scan of every
///        triplet pair sharing a tracklet, and the level sweep against levels
///        worked out by hand, on the deep ToyBarrel8.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK NeighbourStage
#include <boost/test/unit_test.hpp>

#include "ToyBarrel8.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/NeighbourStage.h"
#include "GenericTracking/Tracking/Tracker.h"
#include "GenericTracking/Tracking/TrackletPaths.h"
#include "GenericTracking/Tracking/TripletEndIndex.h"
#include "GenericTracking/Tracking/TripletPaths.h"
#include "GenericTracking/Tracking/TripletPermutation.h"

#include "Oracles.h"
#include "StageOutputComparison.h"
#include "ToyBarrel8Fixture.h"
#include "SerialReferences.h"
#include "FixedPointReferences.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <vector>

namespace o2::gtrk
{

namespace
{

using namespace toybarrel8;
using NeighbourOutput = StageOutput<NeighbourT>;
constexpr std::size_t kTracks = kHelixTracks;
/// The last track drifts in time, its layers in the toys' ticks 0, 1, 1, 2, 2,
/// 3, 3, 4 (GridTestSupport.h): every pair and every triplet of it shares a
/// time, but every other `a-d` span is two ticks -- equal boxes cannot fail
/// every `a-d` with every triplet passing.
constexpr std::size_t kDriftingTrack = 4U;
constexpr std::array<int, 8> kDriftTicks{0, 1, 1, 2, 2, 3, 3, 4};
constexpr std::size_t kNoisePerLayer = 3U;
/// Hits on layer 4 next to a track, moved along z and r-phi [cm]: branches
/// whose triplets all pass, and whose neighbours pass, fail the curvature
/// agreement, or fail only the fit -- which sees the theta kinks and one
/// shared curvature together.
struct Branch {
  std::size_t track = 0U;
  double dz = 0.;
  double dRPhi = 0.;
};
constexpr std::array<Branch, 10> kBranches{{{.track = 0, .dz = 0.008, .dRPhi = 0.003},
                                            {.track = 0, .dz = 0.012, .dRPhi = 0.005},
                                            {.track = 1, .dz = 0.010, .dRPhi = 0.004},
                                            {.track = 1, .dz = 0.016, .dRPhi = 0.002},
                                            {.track = 2, .dz = 0.006, .dRPhi = 0.006},
                                            {.track = 2, .dz = 0.004, .dRPhi = 0.008},
                                            {.track = 3, .dz = 0.012, .dRPhi = 0.004},
                                            {.track = 3, .dz = 0.008, .dRPhi = 0.007},
                                            {.track = 2, .dRPhi = 0.012},
                                            {.track = 3, .dRPhi = 0.016}}};
constexpr LayerID kBranchLayer = 4;

/// The neighbour fixture: kTracks through all eight layers, the last one
/// drifting in time, noise, and the branch hits.
std::vector<Hit> neighbourHits()
{
  std::vector<Hit> hits;
  for (std::size_t track = 0; track < kTracks; ++track) {
    for (std::size_t layer = 0; layer < kLayers; ++layer) {
      const double drift = track == kDriftingTrack ? test::kToyTickNS * static_cast<double>(kDriftTicks[layer]) : 0.;
      hits.push_back({.point = helixPoint(track, static_cast<LayerID>(layer)), .timeNS = 10. + drift});
    }
  }
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    const double radius = layerDistance<Det>(static_cast<LayerID>(layer));
    for (std::size_t noise = 0; noise < kNoisePerLayer; ++noise) {
      const double phi = 0.55 + (1.85 * static_cast<double>(noise)) + (0.23 * static_cast<double>(layer));
      hits.push_back({.point = {.layer = static_cast<LayerID>(layer),
                                .x = static_cast<Scalar>(radius * std::cos(phi)),
                                .y = static_cast<Scalar>(radius * std::sin(phi)),
                                .z = static_cast<Scalar>((noise % 2U == 0U ? 0.65 : -0.55) + (0.03 * static_cast<double>(layer))),
                                .r = static_cast<Scalar>(radius),
                                .phi = static_cast<Scalar>(phi)}});
    }
  }
  for (const auto& branch : kBranches) {
    hits.push_back({.point = moved(helixPoint(branch.track, kBranchLayer), branch.dz, branch.dRPhi)});
  }
  return hits;
}

/// A writable copy of \a triplets, in the same arena order and back at level
/// 1: the sweep raises levels in place, and the tracker hands its triplets
/// out const, already swept.
void copyTriplets(const Triplets& triplets, Triplets& copy)
{
  const auto view = triplets.view();
  copy.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  copy.prepare(triplets.size());
  const auto status = runSerial(copy, [&](Triplets::Index element, auto& emitter) {
    emitter.emit(Triplet{.first = view.first[element],
                         .second = view.second[element],
                         .curvature = view.curvature[element],
                         .varCurvature = view.varCurvature[element],
                         .level = 1U,
                         .kinked = view.kinked[element]});
  });
  BOOST_REQUIRE(status == StageStatus::kComplete);
}

struct Levels {
  std::vector<std::uint32_t> level;
  std::uint32_t sweeps = 0U;
};

/// The level sweep over a copy of \a triplets: serial when \a threads is 0,
/// else on the CPU backend in chunks of \a chunk.
Levels sweepLevels(const Triplets& triplets, const NeighbourOutput& neighbours, int threads = 0, std::uint32_t chunk = 1U)
{
  Triplets copy;
  copyTriplets(triplets, copy);
  Levels result;
  if (threads == 0) {
    result.sweeps = fixedPointLevels(
      [](std::uint32_t elements, const LevelUpdate& update) {
        for (std::uint32_t element = 0; element < elements; ++element) {
          update(element);
        }
      },
      neighbours, copy);
  } else {
    detail::TaskArena cpuArena{threads};
    CpuStage<NeighbourT> cpu{cpuArena};
    cpu.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    result.sweeps = fixedPointLevels([&](std::uint32_t elements, const LevelUpdate& update) { cpu.forEachElement(elements, chunk, update); }, neighbours, copy);
  }
  const auto view = copy.view();
  result.level.assign(view.level, view.level + copy.size());
  return result;
}

/// Every level is one above its highest predecessor's, or 1 without one.
void checkFixedPoint(const Levels& levels, const NeighbourOutput& neighbours)
{
  const auto* const from = neighbours.view().from;
  for (std::uint32_t triplet = 0; triplet < levels.level.size(); ++triplet) {
    std::uint32_t highest = 0U;
    for (auto edge = neighbours.offsets()[triplet]; edge < neighbours.offsets()[triplet + 1U]; ++edge) {
      highest = std::max(highest, levels.level[from[edge]]);
    }
    BOOST_TEST_CONTEXT("triplet " << triplet) { BOOST_TEST(levels.level[triplet] == highest + 1U); }
  }
}

/// The same levels from every schedule: serial, and 1, 2 and 4 threads in
/// chunks of 1 to 16. The number of sweeps may differ; the levels may not.
void checkScheduleInvariance(const Triplets& triplets, const NeighbourOutput& neighbours, const Levels& serial)
{
  for (const int threads : {1, 2, 4}) {
    for (const std::uint32_t chunk : {1U, 2U, 3U, 7U, 16U}) {
      const auto parallel = sweepLevels(triplets, neighbours, threads, chunk);
      BOOST_TEST_CONTEXT("threads " << threads << ", chunk " << chunk) { BOOST_TEST(parallel.level == serial.level, boost::test_tools::per_element()); }
    }
  }
}

/// The tracker's own neighbours and levels are the ones the stages give run
/// by hand.
void checkTrackerOutput(const Triplets& trackerTriplets, const NeighbourOutput& trackerNeighbours, const NeighbourOutput& neighbours, const Levels& levels)
{
  // The independent full-fit reference establishes levels first. The fused
  // tracker keeps precisely the edges that its road walk can use.
  NeighbourOutput longest;
  longest.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  longest.prepare(trackerTriplets.size());
  const auto* const from = neighbours.view().from;
  const auto offsets = neighbours.offsets();
  const auto status = runSerial(longest, [&](std::uint32_t self, auto& emitter) {
    const auto end = offsets[self + 1U];
    for (auto edge = offsets[self]; edge < end; ++edge) {
      if (levels.level[from[edge]] + 1U == levels.level[self]) {
        emitter.emit(Neighbour{.from = from[edge]});
      }
    }
  });
  BOOST_REQUIRE(status == StageStatus::kComplete);
  BOOST_TEST(sameStageOutput(trackerNeighbours, longest));
  const auto view = trackerTriplets.view();
  BOOST_TEST(std::vector<std::uint32_t>(view.level, view.level + trackerTriplets.size()) == levels.level, boost::test_tools::per_element());
}

/// The triplet path bounds as the tracker builds them, at \a minClusters and
/// \a parameters' road ends.
struct PathBounds {
  TrackletPaths<Det> trackletPaths;
  TripletPermutation<Det> permutation;
  TripletPaths<Det> tripletPaths;

  TrackletPathBounds build(IterationParameters<Det> parameters, std::uint32_t minClusters, const LayerRowMap<Det>& sources, const TrackletOutput& tracklets,
                           const Triplets& triplets)
  {
    parameters.minClusters = minClusters;
    const auto pool = std::make_shared<BoundedMemoryResource>();
    trackletPaths.setMemoryPool(pool);
    permutation.setMemoryPool(pool);
    tripletPaths.setMemoryPool(pool);
    trackletPaths.build(parameters, sources, tracklets, SerialForEach{});
    BOOST_REQUIRE(permutation.build(parameters, sources, tracklets, trackletPaths.bounds(), SerialForEach{}, 256U) == StageStatus::kComplete);
    tripletPaths.build(parameters, sources, tracklets, permutation, placedTriplets(triplets), SerialForEach{});
    return tripletPaths.bounds();
  }
};

/// Tracklets, triplets, the end index and the neighbours of \a fixture's
/// whole timeframe, with the tracker's path bounds; \a tracker owns the
/// tracklets and triplets.
void buildNeighbours(const Fixture& fixture, const IterationParameters<Det>& parameters, Controller& tracker, TripletEndIndex& ending,
                     NeighbourOutput& neighbours, const TrackletOutput*& tracklets, const Triplets*& triplets)
{
  const Frame::TimeSlice slice{*fixture.frame, fixture.geometry.allBlocks(), 0};
  BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
  const auto result = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  tracklets = result.tracklets;
  triplets = result.triplets;
  ending.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  BOOST_REQUIRE(ending.build(tracklets->size(), *triplets, SerialForEach{}, std::numeric_limits<std::uint32_t>::max()) == StageStatus::kComplete);
  PathBounds paths;
  const auto bounds = paths.build(parameters, parameters.minClusters, tracker.sources(), *tracklets, *triplets);
  const NeighbourStage<Det> stage{parameters, *fixture.frame, *tracklets, *triplets, ending, bounds};
  neighbours.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  neighbours.prepare(triplets->size());
  BOOST_REQUIRE(runSerial(neighbours, stage) == StageStatus::kComplete);
}

} // namespace

BOOST_AUTO_TEST_CASE(NeighboursMatchBruteForceAndEveryTrueChainIsLinked)
{
  const Fixture fixture{neighbourHits()};
  const Frame::TimeSlice slice{*fixture.frame, fixture.geometry.allBlocks(), 0};
  const std::array<double, kLayers> material{};
  const auto scattering = scatteringTable<Det>(material);
  const auto parameters = iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scattering);
  Controller tracker{iterationSettings(), Fixture::kBeamConstraint, scattering, std::make_shared<BoundedMemoryResource>(), 2U, 2};
  BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
  const auto result = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  const auto& tracklets = *result.tracklets;
  const auto& triplets = *result.triplets;
  BOOST_REQUIRE(triplets.size() > 0U);

  TripletEndIndex ending;
  ending.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  BOOST_REQUIRE(ending.build(tracklets.size(), triplets, SerialForEach{}, std::numeric_limits<std::uint32_t>::max()) == StageStatus::kComplete);
  // At 3 clusters no triplet path is too short: the stage's own cut alone.
  PathBounds paths;
  const auto bounds = paths.build(parameters, 3U, tracker.sources(), tracklets, triplets);
  const NeighbourStage<Det> stage{parameters, *fixture.frame, tracklets, triplets, ending, bounds};
  NeighbourOutput serial;
  serial.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  serial.prepare(triplets.size());
  BOOST_REQUIRE(runSerial(serial, stage) == StageStatus::kComplete);

  // Brute force: every pair T1, T2 with T1.second == T2.first, the whole cut,
  // predecessors in arena order.
  const auto trackletView = tracklets.view();
  const auto tripletView = triplets.view();
  const auto neighbourView = serial.view();
  std::size_t candidates = 0U;
  std::size_t curvatureRejected = 0U;
  std::size_t timeRejected = 0U;
  std::size_t fitRejected = 0U;
  for (std::uint32_t to = 0; to < triplets.size(); ++to) {
    std::vector<std::uint32_t> expected;
    const auto shared = tripletView.first[to];
    const auto last = tripletView.second[to];
    for (std::uint32_t from = 0; from < triplets.size(); ++from) {
      if (tripletView.second[from] != shared) {
        continue;
      }
      ++candidates;
      const auto first = tripletView.first[from];
      const auto a = timedPointAt(fixture.frame->measurements(), trackletView.source[first]);
      const auto b = timedPointAt(fixture.frame->measurements(), trackletView.source[shared]);
      const auto c = timedPointAt(fixture.frame->measurements(), trackletView.target[shared]);
      const auto d = timedPointAt(fixture.frame->measurements(), trackletView.target[last]);
      const LayerID secondLayer = parameters.linkSources[trackletView.link[shared]];
      const LayerID thirdLayer = parameters.linkTargets[trackletView.link[shared]];
      if (!neighbourCurvatureCompatible(parameters, tripletView.curvature[from], tripletView.varCurvature[from], tripletView.curvature[to],
                                        tripletView.varCurvature[to])) {
        ++curvatureRejected;
        continue;
      }
      if (!neighbourTimeCompatible(a, d)) {
        ++timeRejected;
        continue;
      }
      TripletFitResult fit;
      if (!fitAdjacentTriplets(a.point, b.point, c.point, d.point, secondLayer, thirdLayer, parameters.scattering, parameters.pairCut, fit) ||
          !(fit.chi2 <= parameters.neighbourChi2Max)) {
        ++fitRejected;
        continue;
      }
      BOOST_TEST(neighbourCompatible(parameters, secondLayer, thirdLayer, a, b, c, d, fit));
      expected.push_back(from);
    }
    std::vector<std::uint32_t> emitted(neighbourView.from + serial.offsets()[to], neighbourView.from + serial.offsets()[to + 1U]);
    BOOST_TEST_CONTEXT("triplet " << to) { BOOST_TEST(emitted == expected, boost::test_tools::per_element()); }
  }
  BOOST_TEST_MESSAGE("triplets " << triplets.size() << ", candidates " << candidates << ", neighbours " << serial.size() << ", rejected by curvature "
                                 << curvatureRejected << ", time " << timeRejected << ", fit " << fitRejected);
  BOOST_TEST(curvatureRejected > 0U);
  BOOST_TEST(timeRejected > 0U);
  BOOST_TEST(fitRejected > 0U);

  // Each track crosses all eight layers: six true triplets, linked by five
  // true neighbours -- except the drifting one, whose links 1, 3 and 5 span
  // two ticks `a-d`.
  for (std::size_t track = 0; track < kTracks; ++track) {
    std::array<std::uint32_t, kLayers - 2U> chain{};
    for (std::size_t layer = 0; layer + 2U < kLayers; ++layer) {
      chain[layer] = findTriplet(tracklets, triplets, fixture.row(Fixture::logicalIndex(track, layer)), fixture.row(Fixture::logicalIndex(track, layer + 1U)),
                                 fixture.row(Fixture::logicalIndex(track, layer + 2U)));
      BOOST_REQUIRE_MESSAGE(chain[layer] < triplets.size(), "track " << track << " misses the triplet from layer " << layer);
    }
    for (std::size_t link = 1; link < chain.size(); ++link) {
      const auto* begin = neighbourView.from + serial.offsets()[chain[link]];
      const auto* end = neighbourView.from + serial.offsets()[chain[link] + 1U];
      const bool linked = std::find(begin, end, chain[link - 1U]) != end;
      BOOST_TEST(linked == (track != kDriftingTrack || link % 2U == 0U), "track " << track << " link " << link);
    }
  }

  for (const int threads : {1, 2, 4}) {
    for (const NeighbourOutput::Index chunk : {1U, 2U, 3U, 7U, 16U}) {
      NeighbourOutput parallel;
      const auto parallelPool = std::make_shared<BoundedMemoryResource>();
      parallel.setMemoryPool(parallelPool);
      parallel.prepare(triplets.size());
      detail::TaskArena cpuArena{threads};
      CpuStage<NeighbourT> cpu{cpuArena};
      cpu.setMemoryPool(parallelPool);
      BOOST_REQUIRE(cpu.run(parallel, chunk, stage) == StageStatus::kComplete);
      BOOST_TEST(sameStageOutput(parallel, serial), "threads " << threads << ", chunk " << chunk);
    }
  }

  // The levels over this graph -- branches, noise, a broken track -- are a
  // fixed point, and every schedule reaches the same one.
  const auto levels = sweepLevels(triplets, serial);
  BOOST_TEST_MESSAGE("chained " << chainedTriplets(serial).size() << ", sweeps " << levels.sweeps);
  checkFixedPoint(levels, serial);
  checkScheduleInvariance(triplets, serial, levels);

  // The tracker runs the same stages itself.
  BOOST_REQUIRE(result.neighbours != nullptr);
  checkTrackerOutput(triplets, *result.neighbours, serial, levels);
}

namespace
{

struct Expected {
  std::array<std::uint32_t, 3> hits{};
  std::uint32_t level = 0U;
};

} // namespace

BOOST_AUTO_TEST_CASE(LevelsOnTheDeepBarrelMatchTheirHandCount)
{
  using namespace toybarrel8::levels;
  const Fixture fixture{levelHits()};
  const std::array<double, kLayers> material{};
  const auto scattering = scatteringTable<Det>(material);
  const auto parameters = iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scattering);
  Controller tracker{iterationSettings(), Fixture::kBeamConstraint, scattering, std::make_shared<BoundedMemoryResource>(), 2U, 2};
  TripletEndIndex ending;
  NeighbourOutput neighbours;
  const StageOutput<TrackletT>* tracklets = nullptr;
  const Triplets* triplets = nullptr;
  buildNeighbours(fixture, parameters, tracker, ending, neighbours, tracklets, triplets);

  // Track 0: six triplets in a row, levels 1 to 6, and a branch through the
  // second hit on layer 3 that runs alongside from level 2 to 4 and rejoins
  // at (4, 5, 6), which has two predecessors at level 4. Track 1 lost its
  // layer-3 hit, and links join only neighbouring layers: its chain breaks
  // and restarts at level 1. Its first piece, (9, 10, 11), and track 2's
  // single triplet are on no path of four clusters and not made.
  const std::array<Expected, 11> expected{{
    {.hits = {0, 1, 2}, .level = 1},
    {.hits = {1, 2, 3}, .level = 2},
    {.hits = {1, 2, kSecondOnLayer3}, .level = 2},
    {.hits = {2, 3, 4}, .level = 3},
    {.hits = {2, kSecondOnLayer3, 4}, .level = 3},
    {.hits = {3, 4, 5}, .level = 4},
    {.hits = {kSecondOnLayer3, 4, 5}, .level = 4},
    {.hits = {4, 5, 6}, .level = 5},
    {.hits = {5, 6, 7}, .level = 6},
    {.hits = {kTrack1[4], kTrack1[5], kTrack1[6]}, .level = 1},
    {.hits = {kTrack1[5], kTrack1[6], kTrack1[7]}, .level = 2},
  }};
  BOOST_REQUIRE(triplets->size() == expected.size());
  const auto levels = sweepLevels(*triplets, neighbours);
  // The tracker's own sweep left the same levels on its triplets.
  const auto swept = triplets->view();
  BOOST_TEST(std::vector<std::uint32_t>(swept.level, swept.level + triplets->size()) == levels.level, boost::test_tools::per_element());
  for (const auto& triplet : expected) {
    const auto index = findTriplet(*tracklets, *triplets, fixture.row(triplet.hits[0]), fixture.row(triplet.hits[1]), fixture.row(triplet.hits[2]));
    BOOST_TEST_CONTEXT("triplet " << triplet.hits[0] << " " << triplet.hits[1] << " " << triplet.hits[2])
    {
      BOOST_REQUIRE(index < triplets->size());
      BOOST_TEST(levels.level[index] == triplet.level);
    }
  }
  // Only triplets with predecessors are swept: eight on track 0, one on track
  // 1. Levels up to 6 need at most five sweeps that change one.
  BOOST_TEST(chainedTriplets(neighbours).size() == 9U);
  BOOST_TEST(levels.sweeps >= 1U);
  BOOST_TEST(levels.sweeps <= static_cast<std::uint32_t>(maxChainLayers<Det>) - 3U);
  checkFixedPoint(levels, neighbours);
  checkScheduleInvariance(*triplets, neighbours, levels);
}

} // namespace o2::gtrk
