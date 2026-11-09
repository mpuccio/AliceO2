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

/// \file testRoads.cxx
/// \brief Where roads start and how many leave each triplet, against a hand
///        count and a brute-force walk, on the deep ToyBarrel8.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Roads
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/RoadCount.h"

#include "ToyBarrel8Fixture.h"
#include "SerialReferences.h"
#include "FixedPointReferences.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <numbers>
#include <set>
#include <span>
#include <vector>

namespace o2::gtrk
{

namespace
{

using namespace toybarrel8;

/// A road as its triplets, innermost first.
using Road = std::vector<std::uint32_t>;

struct Counted {
  std::vector<std::uint32_t> successorAbove;
  std::vector<std::uint32_t> roads;
  std::vector<std::uint32_t> starts;
  std::vector<std::uint32_t> roadOffsets;
  std::vector<Road> enumerated; ///< every road by its index
};

constexpr std::size_t kMaxRoadTriplets = maxChainLayers<Det> - 2;

/// Every road by its index through the walk, one element per road.
template <typename ForEach>
std::vector<Road> enumerate(ForEach forEach, const RoadWalk& walk, RoadCount::Index roads)
{
  std::vector<std::array<std::uint32_t, kMaxRoadTriplets>> triplets(roads);
  std::vector<std::uint32_t> lengths(roads);
  forEach(roads, [&](RoadCount::Index road) { lengths[road] = walk(road, triplets[road].data()); });
  std::vector<Road> result(roads);
  for (RoadCount::Index road = 0; road < roads; ++road) {
    result[road].assign(triplets[road].begin(), triplets[road].begin() + lengths[road]);
  }
  return result;
}

/// The road count of one tracker run, serial when \a threads is 0, else on
/// the CPU backend in chunks of \a chunk.
template <typename Result>
Counted countRoads(const Result& run, const IterationParameters<Det>& parameters, std::span<const MeasurementRange> core, std::uint32_t minStartLevel, int threads = 0, RoadCount::Index chunk = 1U, std::uint64_t endLayers = 0U)
{
  const RoadStartRule rule{.minLevel = minStartLevel, .endLayers = endLayers, .linkTargets = parameters.linkTargets};
  RoadCount count;
  count.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  const auto serial = [](RoadCount::Index elements, const auto& body) {
    for (RoadCount::Index element = 0; element < elements; ++element) {
      body(element);
    }
  };
  detail::TaskArena cpuArena{std::max(threads, 1)};
  CpuStage<NeighbourT> cpu{cpuArena};
  cpu.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  const auto parallel = [&](RoadCount::Index elements, const auto& body) { cpu.forEachElement(elements, chunk, body); };
  const auto status = threads == 0 ? fixedPointRoads(count, serial, parameters.linkSources, core.data(), rule, *run.tracklets, *run.triplets, *run.neighbours) : fixedPointRoads(count, parallel, parameters.linkSources, core.data(), rule, *run.tracklets, *run.triplets, *run.neighbours);
  BOOST_REQUIRE(status == StageStatus::kComplete);
  // An independently scheduled ordered pass uses the final graph, rather
  // than the fitter's generation. It must give the reference's numbering.
  RoadCount ordered;
  ordered.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  const auto compareOrdered = [&](const auto& forEach) {
    ordered.prepareOrdered(forEach, run.triplets->size());
    const auto state = ordered.ordered(parameters.linkSources, core.data());
    const auto triplets = run.triplets->view();
    const auto offsets = run.neighbours->offsets();
    const auto* const from = run.neighbours->view().from;
    const auto highest = *std::max_element(triplets.level, triplets.level + run.triplets->size());
    // Replaying the complete ordered pass models a staging retry. Counts are
    // overwritten; successor marks and overflow must remain idempotent.
    for (int attempt = 0; attempt < 2; ++attempt) {
      for (std::uint32_t level = 1U; level <= highest; ++level) {
        forEach(run.triplets->size(), [&](std::uint32_t element) {
          if (triplets.level[element] != level) {
            return;
          }
          state.base(element, run.tracklets->view(), triplets);
          if (level > 1U) {
            std::uint64_t sum = 0U;
            for (auto edge = offsets[element]; edge < offsets[element + 1U]; ++edge) {
              if (triplets.level[from[edge]] + 1U == level) {
                state.add(from[edge], sum);
              }
            }
            state.finish(element, sum);
          }
        });
      }
      const auto orderedStatus = ordered.runPrepared(forEach, rule, *run.tracklets, *run.triplets, chunk);
      BOOST_REQUIRE(orderedStatus == StageStatus::kComplete);
      BOOST_TEST(std::ranges::equal(ordered.successorAbove(), count.successorAbove()));
      BOOST_TEST(std::ranges::equal(ordered.roads(), count.roads()));
      BOOST_TEST(std::ranges::equal(ordered.starts(), count.starts()));
      BOOST_TEST(std::ranges::equal(ordered.roadOffsets(), count.roadOffsets()));
      BOOST_TEST((enumerate(forEach, ordered.walk(*run.triplets, *run.neighbours), ordered.roadCount()) == enumerate(forEach, count.walk(*run.triplets, *run.neighbours), count.roadCount())));
    }
  };
  if (run.triplets->size() != 0U) {
    if (threads == 0) {
      compareOrdered(serial);
    } else {
      compareOrdered(parallel);
    }
  }
  const auto walk = count.walk(*run.triplets, *run.neighbours);
  return {
    .successorAbove = {count.successorAbove().begin(), count.successorAbove().end()},
    .roads = {count.roads().begin(), count.roads().end()},
    .starts = {count.starts().begin(), count.starts().end()},
    .roadOffsets = {count.roadOffsets().begin(), count.roadOffsets().end()},
    .enumerated = threads == 0 ? enumerate(serial, walk, count.roadCount()) : enumerate(parallel, walk, count.roadCount()),
  };
}

/// By brute force: the successor flags from every edge, the roads by a
/// depth-first walk down through predecessors one level below to owned
/// level-1 triplets, the starts by their definition; and, with every level
/// allowed to start, that the walks from the starts cover every triplet with
/// a road.
template <typename Result>
void checkAgainstBruteForce(const Result& run, const IterationParameters<Det>& parameters, std::span<const MeasurementRange> core, std::uint32_t minStartLevel, const Counted& counted)
{
  const auto tracklets = run.tracklets->view();
  const auto triplets = run.triplets->view();
  const auto& neighbours = *run.neighbours;
  const auto* const from = neighbours.view().from;
  const auto count = run.triplets->size();
  std::vector<std::uint32_t> successorAbove(count, 0U);
  for (std::uint32_t to = 0; to < count; ++to) {
    for (auto edge = neighbours.offsets()[to]; edge < neighbours.offsets()[to + 1U]; ++edge) {
      successorAbove[from[edge]] |= triplets.level[from[edge]] + 1U == triplets.level[to] ? 1U : 0U;
    }
  }
  // Recursive on purpose: the reference walks every road, the stage does not.
  const auto walk = [&](const auto& self, std::uint32_t triplet, std::vector<bool>* visited) -> std::uint64_t { // NOLINT(misc-no-recursion)
    if (visited != nullptr) {
      (*visited)[triplet] = true;
    }
    if (triplets.level[triplet] == 1U) {
      const auto first = triplets.first[triplet];
      const auto owned = core[parameters.linkSources[tracklets.link[first]]];
      return tracklets.source[first] >= owned.begin && tracklets.source[first] < owned.end ? 1U : 0U;
    }
    std::uint64_t roads = 0U;
    for (auto edge = neighbours.offsets()[triplet]; edge < neighbours.offsets()[triplet + 1U]; ++edge) {
      if (triplets.level[from[edge]] + 1U == triplets.level[triplet]) {
        roads += self(self, from[edge], visited);
      }
    }
    return roads;
  };
  std::vector<std::uint32_t> starts;
  std::vector<bool> onRoad(count, false);
  for (std::uint32_t triplet = 0; triplet < count; ++triplet) {
    const auto roads = walk(walk, triplet, nullptr);
    BOOST_TEST_CONTEXT("triplet " << triplet)
    {
      BOOST_TEST(counted.successorAbove[triplet] == successorAbove[triplet]);
      BOOST_TEST(counted.roads[triplet] == roads);
    }
    if (triplets.level[triplet] >= minStartLevel && successorAbove[triplet] == 0U && roads != 0U) {
      starts.push_back(triplet);
      walk(walk, triplet, &onRoad);
    }
  }
  BOOST_TEST(counted.starts == starts, boost::test_tools::per_element());

  // The roads by a depth-first walk from each start in turn, predecessors in
  // the neighbour stage's order: the walk by index must give the same roads
  // in the same order, each once.
  std::vector<Road> roads;
  Road path;
  const auto descend = [&](const auto& self, std::uint32_t triplet) -> void { // NOLINT(misc-no-recursion)
    path.push_back(triplet);
    if (triplets.level[triplet] == 1U) {
      if (walk(walk, triplet, nullptr) != 0U) {
        roads.emplace_back(path.rbegin(), path.rend());
      }
    } else {
      for (auto edge = neighbours.offsets()[triplet]; edge < neighbours.offsets()[triplet + 1U]; ++edge) {
        if (triplets.level[from[edge]] + 1U == triplets.level[triplet]) {
          self(self, from[edge]);
        }
      }
    }
    path.pop_back();
  };
  for (const auto start : starts) {
    descend(descend, start);
  }
  BOOST_REQUIRE(counted.enumerated.size() == roads.size());
  BOOST_TEST(counted.roadOffsets.back() == roads.size());
  for (std::size_t road = 0; road < roads.size(); ++road) {
    BOOST_TEST_CONTEXT("road " << road) { BOOST_TEST(counted.enumerated[road] == roads[road], boost::test_tools::per_element()); }
  }
  if (minStartLevel == 1U) {
    for (std::uint32_t triplet = 0; triplet < count; ++triplet) {
      BOOST_TEST((onRoad[triplet] || walk(walk, triplet, nullptr) == 0U), "triplet " << triplet << " lies on no road");
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(RoadsOnTheDeepBarrelMatchTheirHandCount)
{
  using namespace toybarrel8::levels;
  const Fixture fixture{levelHits()};
  const std::array<double, kLayers> material{};
  const auto scattering = scatteringTable<Det>(material);
  const auto parameters = iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scattering);
  Controller tracker{iterationSettings(), Fixture::kBeamConstraint, scattering, std::make_shared<BoundedMemoryResource>(), 2U, 2};
  const Frame::TimeSlice slice{*fixture.frame, fixture.geometry.allBlocks(), 0};
  BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
  const auto run = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
  BOOST_REQUIRE(run.status == StageStatus::kComplete);
  BOOST_REQUIRE(parameters.minRoadStartLevel == 2U);
  const auto counted = countRoads(run, parameters, slice.core, parameters.minRoadStartLevel);

  // Track 0 runs through all eight layers with a second hit on layer 3: the
  // branch doubles the roads from (4, 5, 6) up, and (5, 6, 7), level 6, is
  // its one start. Track 1 lost its layer-3 hit: (12, 13, 14) -> (13, 14, 15)
  // at level 2; its (9, 10, 11) is on no path of four clusters, and neither is
  // track 2's lone triplet, so neither is made.
  struct Expected {
    std::array<std::uint32_t, 3> hits{};
    std::uint32_t successorAbove = 0U;
    std::uint32_t roads = 0U;
    bool start = false;
  };
  const std::array<Expected, 11> expected{{
    {.hits = {0, 1, 2}, .successorAbove = 1, .roads = 1},
    {.hits = {1, 2, 3}, .successorAbove = 1, .roads = 1},
    {.hits = {1, 2, kSecondOnLayer3}, .successorAbove = 1, .roads = 1},
    {.hits = {2, 3, 4}, .successorAbove = 1, .roads = 1},
    {.hits = {2, kSecondOnLayer3, 4}, .successorAbove = 1, .roads = 1},
    {.hits = {3, 4, 5}, .successorAbove = 1, .roads = 1},
    {.hits = {kSecondOnLayer3, 4, 5}, .successorAbove = 1, .roads = 1},
    {.hits = {4, 5, 6}, .successorAbove = 1, .roads = 2},
    {.hits = {5, 6, 7}, .successorAbove = 0, .roads = 2, .start = true},
    {.hits = {kTrack1[4], kTrack1[5], kTrack1[6]}, .successorAbove = 1, .roads = 1},
    {.hits = {kTrack1[5], kTrack1[6], kTrack1[7]}, .successorAbove = 0, .roads = 1, .start = true},
  }};
  BOOST_REQUIRE(run.triplets->size() == expected.size());
  std::vector<std::uint32_t> starts;
  std::uint32_t fullTrackStart = run.triplets->size();
  for (const auto& triplet : expected) {
    const auto index = findTriplet(*run.tracklets, *run.triplets, fixture.row(triplet.hits[0]), fixture.row(triplet.hits[1]), fixture.row(triplet.hits[2]));
    BOOST_TEST_CONTEXT("triplet " << triplet.hits[0] << " " << triplet.hits[1] << " " << triplet.hits[2])
    {
      BOOST_REQUIRE(index < run.triplets->size());
      BOOST_TEST(counted.successorAbove[index] == triplet.successorAbove);
      BOOST_TEST(counted.roads[index] == triplet.roads);
    }
    if (triplet.start) {
      starts.push_back(index);
    }
    if (triplet.hits[0] == 5U) {
      fullTrackStart = index;
    }
  }
  std::ranges::sort(starts);
  BOOST_TEST(counted.starts == starts, boost::test_tools::per_element());

  // Seven clusters or more: only track 0's start, level 6, remains; nine is
  // beyond its eight.
  auto seven = iterationSettings();
  seven.minClusters = 7;
  const auto sevenParameters = iterationParameters<Det>(seven, Fixture::kBeamConstraint, scattering);
  BOOST_TEST(sevenParameters.minRoadStartLevel == 5U);
  const auto longOnly = countRoads(run, sevenParameters, slice.core, sevenParameters.minRoadStartLevel);
  BOOST_TEST(longOnly.starts == std::vector<std::uint32_t>{fullTrackStart}, boost::test_tools::per_element());
  BOOST_TEST(longOnly.roads == counted.roads, boost::test_tools::per_element());
  BOOST_TEST(countRoads(run, parameters, slice.core, 7U).starts.empty());

  // Three roads: the two of track 0, which part at levels 2 to 4 -- one
  // through each hit on layer 3 -- and meet again below and above; then
  // track 1's.
  const auto& roads = counted.enumerated;
  BOOST_REQUIRE(roads.size() == 3U);
  BOOST_REQUIRE(roads[0].size() == 6U);
  BOOST_REQUIRE(roads[1].size() == 6U);
  BOOST_TEST(roads[2].size() == 2U);
  BOOST_TEST(roads[0][0] == roads[1][0]);
  for (std::size_t level = 1; level < 4U; ++level) {
    BOOST_TEST(roads[0][level] != roads[1][level]);
  }
  BOOST_TEST(roads[0][4] == roads[1][4]);
  BOOST_TEST(roads[0][5] == fullTrackStart);
  BOOST_TEST(roads[1][5] == fullTrackStart);
  BOOST_TEST(roads[2][1] == starts.back());
  checkAgainstBruteForce(run, parameters, slice.core, 1U, countRoads(run, parameters, slice.core, 1U));

  // With an end-layer mask a road starts at every triplet ending on one of its
  // layers, extended or not: on layer 3, track 0's two prefixes (1, 2, 3) and
  // (1, 2, 3'), a road of four clusters each; with layer 7 as well, the two
  // starts above and their three roads besides.
  std::vector<std::uint32_t> onThree;
  for (const auto third : {3U, kSecondOnLayer3}) {
    onThree.push_back(findTriplet(*run.tracklets, *run.triplets, fixture.row(1U), fixture.row(2U), fixture.row(third)));
  }
  std::ranges::sort(onThree);
  const auto three = countRoads(run, parameters, slice.core, 2U, 0, 1U, std::uint64_t{1U} << 3U);
  BOOST_TEST(three.starts == onThree, boost::test_tools::per_element());
  BOOST_TEST(three.roads == counted.roads, boost::test_tools::per_element());
  BOOST_REQUIRE(three.enumerated.size() == 2U);
  BOOST_TEST(three.enumerated[0].size() == 2U);
  BOOST_TEST(three.enumerated[1].size() == 2U);
  auto threeAndSeven = onThree;
  threeAndSeven.insert(threeAndSeven.end(), starts.begin(), starts.end());
  std::ranges::sort(threeAndSeven);
  const auto both = countRoads(run, parameters, slice.core, 2U, 0, 1U, (std::uint64_t{1U} << 3U) | (std::uint64_t{1U} << 7U));
  BOOST_TEST(both.starts == threeAndSeven, boost::test_tools::per_element());
  BOOST_TEST(both.enumerated.size() == 5U);
  // A mask of layers no triplet ends on starts nothing.
  BOOST_TEST(countRoads(run, parameters, slice.core, 2U, 0, 1U, std::uint64_t{1U} << 1U).starts.empty());
}

BOOST_AUTO_TEST_CASE(RoadsMatchABruteForceWalkAndOnlyTheCoresAreCounted)
{
  const Fixture fixture{sliceHits(), kSliceDurationNS};
  const std::array<double, kLayers> material{};
  const auto scattering = scatteringTable<Det>(material);
  const auto parameters = iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scattering);
  Controller tracker{iterationSettings(), Fixture::kBeamConstraint, scattering, std::make_shared<BoundedMemoryResource>(), 2U, 2};
  const std::int64_t haloPS = sliceReach(frameBounds(*fixture.frame), fixture.geometry).haloPS();
  const auto blockCount = fixture.geometry.blockCount();

  // The whole timeframe, then three slices whose cores own parts of it.
  std::size_t owned = 0U;
  std::size_t total = 0U;
  for (const auto& [begin, end] : std::array<std::array<std::int32_t, 2>, 4>{{{0, blockCount}, {0, 20}, {20, 35}, {35, blockCount}}}) {
    const Frame::TimeSlice slice{*fixture.frame, {.begin = begin, .end = end}, haloPS};
    BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
    const auto run = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
    BOOST_REQUIRE(run.status == StageStatus::kComplete);
    BOOST_TEST_CONTEXT("core blocks " << begin << "-" << end)
    {
      for (const std::uint32_t minStartLevel : {1U, 2U, 5U}) {
        const auto counted = countRoads(run, parameters, slice.core, minStartLevel);
        checkAgainstBruteForce(run, parameters, slice.core, minStartLevel, counted);
        for (const int threads : {1, 2, 4}) {
          for (const RoadCount::Index chunk : {1U, 3U, 16U}) {
            const auto parallel = countRoads(run, parameters, slice.core, minStartLevel, threads, chunk);
            BOOST_TEST(parallel.successorAbove == counted.successorAbove, boost::test_tools::per_element());
            BOOST_TEST(parallel.roads == counted.roads, boost::test_tools::per_element());
            BOOST_TEST(parallel.starts == counted.starts, boost::test_tools::per_element());
            BOOST_TEST((parallel.enumerated == counted.enumerated));
          }
        }
      }
      // The roads of the minimum length, summed over the starts.
      const auto counted = countRoads(run, parameters, slice.core, parameters.minRoadStartLevel);
      std::size_t roads = 0U;
      for (const auto start : counted.starts) {
        roads += counted.roads[start];
      }
      (begin == 0 && end == blockCount ? total : owned) += roads;
    }
  }
  // The three cores tile the timeframe, and each road is owned by one.
  BOOST_TEST(total > 0U);
  BOOST_TEST(owned == total);
}

namespace
{

/// The slicing fixture's tracks, a third of them from layer 1 and a third
/// from layer 2 on, so their roads start on those layers.
std::vector<Hit> staggeredHits()
{
  std::vector<Hit> hits;
  for (std::size_t track = 0; track < kSliceTracks; ++track) {
    const double phi0 = std::fmod(0.3 + (0.9 * static_cast<double>(track)), 2 * std::numbers::pi_v<double>);
    const double slope = (track % 2U == 0U ? 0.03 : -0.02) + (0.002 * static_cast<double>(track));
    const double pt = 1.2 + (0.3 * static_cast<double>(track % 4U));
    const double start = 37. + (83. * static_cast<double>(track));
    for (std::size_t layer = track % 3U; layer < kLayers; ++layer) {
      hits.push_back({.point = helixAt(phi0, slope, pt, static_cast<LayerID>(layer)), .timeNS = start});
    }
  }
  return hits;
}

} // namespace

BOOST_AUTO_TEST_CASE(AStartLayerMaskKeepsExactlyTheRoadsStartingOnIt)
{
  const Fixture fixture{staggeredHits(), kSliceDurationNS};
  const std::array<double, kLayers> material{};
  const auto scattering = scatteringTable<Det>(material);
  const auto parameters = iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scattering);
  Controller tracker{iterationSettings(), Fixture::kBeamConstraint, scattering, std::make_shared<BoundedMemoryResource>(), 2U, 2};
  const std::int64_t haloPS = sliceReach(frameBounds(*fixture.frame), fixture.geometry).haloPS();
  const Frame::TimeSlice slice{*fixture.frame, {.begin = 0, .end = fixture.geometry.blockCount()}, haloPS};
  BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
  const auto run = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
  BOOST_REQUIRE(run.status == StageStatus::kComplete);
  const auto all = countRoads(run, parameters, slice.core, parameters.minRoadStartLevel);
  const auto tracklets = run.tracklets->view();
  const auto triplets = run.triplets->view();
  const auto startLayer = [&](const Road& road) { return parameters.linkSources[tracklets.link[triplets.first[road.front()]]]; };
  std::set<LayerID> starts;
  for (const auto& road : all.enumerated) {
    starts.insert(startLayer(road));
  }
  BOOST_REQUIRE((starts == std::set<LayerID>{0, 1, 2}));
  for (const std::uint64_t mask : {std::uint64_t{1U}, std::uint64_t{1U} << 2U, std::uint64_t{6U}}) {
    const auto cores = roadCores(std::span<const MeasurementRange, kLayers>{slice.core}, mask);
    const auto masked = countRoads(run, parameters, cores, parameters.minRoadStartLevel);
    checkAgainstBruteForce(run, parameters, cores, parameters.minRoadStartLevel, masked);
    std::vector<Road> expected;
    std::ranges::copy_if(all.enumerated, std::back_inserter(expected), [&](const Road& road) { return ((mask >> startLayer(road)) & 1U) != 0U; });
    BOOST_TEST_CONTEXT("mask " << mask)
    {
      BOOST_TEST(!expected.empty());
      BOOST_TEST((masked.enumerated == expected));
    }
  }
}

namespace
{

/// A road as its clusters, innermost first.
using RoadClusters = std::vector<std::uint32_t>;

/// Every road the slices of \a widthBlocks blocks own, the first
/// \a firstWidth wide, as clusters, sorted.
std::vector<RoadClusters> roadsOfSlicing(const Fixture& fixture, std::int32_t widthBlocks, std::int32_t firstWidth)
{
  const std::array<double, kLayers> material{};
  const auto scattering = scatteringTable<Det>(material);
  const auto parameters = iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scattering);
  Controller tracker{iterationSettings(), Fixture::kBeamConstraint, scattering, std::make_shared<BoundedMemoryResource>(), 2U, 2};
  const std::int64_t haloPS = sliceReach(frameBounds(*fixture.frame), fixture.geometry).haloPS();
  const auto blockCount = fixture.geometry.blockCount();
  std::vector<RoadClusters> result;
  for (std::int32_t begin = 0; begin < blockCount;) {
    const auto end = std::min(blockCount, begin + (begin == 0 ? firstWidth : widthBlocks));
    const Frame::TimeSlice slice{*fixture.frame, {.begin = begin, .end = end}, haloPS};
    BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
    const auto run = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
    BOOST_REQUIRE(run.status == StageStatus::kComplete);
    const auto counted = countRoads(run, parameters, slice.core, parameters.minRoadStartLevel);
    const auto tracklets = run.tracklets->view();
    const auto triplets = run.triplets->view();
    const auto cluster = [&](std::uint32_t row) { return fixture.frame->measurement(row).clusterIndex; };
    for (const auto& road : counted.enumerated) {
      const auto first = triplets.first[road.front()];
      RoadClusters clusters{cluster(tracklets.source[first]), cluster(tracklets.target[first])};
      for (const auto triplet : road) {
        clusters.push_back(cluster(tracklets.target[triplets.second[triplet]]));
      }
      result.push_back(clusters);
    }
    begin = end;
  }
  std::ranges::sort(result);
  return result;
}

} // namespace

BOOST_AUTO_TEST_CASE(TheRoadsDoNotDependOnTheSlicing)
{
  const Fixture fixture{sliceHits(), kSliceDurationNS};
  const auto blockCount = fixture.geometry.blockCount();
  const auto whole = roadsOfSlicing(fixture, blockCount, blockCount);
  BOOST_TEST((std::ranges::adjacent_find(whole) == whole.end()), "a road owned twice");
  // Every track: one road through all eight layers, two where it is branched.
  BOOST_TEST(whole.size() == kSliceTracks + ((kSliceTracks + 2U) / 3U));
  for (const auto& [width, first] : std::array<std::array<std::int32_t, 2>, 5>{{{1, 1}, {2, 2}, {2, 1}, {3, 2}, {7, 4}}}) {
    BOOST_TEST_CONTEXT("slices of " << width << " blocks, the first " << first) { BOOST_TEST((roadsOfSlicing(fixture, width, first) == whole)); }
  }
}

BOOST_AUTO_TEST_CASE(OrderedCountsSaturateAndKeepOverflowAcrossReplay)
{
  constexpr std::uint32_t kNodes = 66U;
  std::array<std::uint32_t, kNodes> roads{};
  std::array<std::uint32_t, kNodes> successors{};
  std::uint32_t overflow = 0U;
  const OrderedRoadCounts state{.roads = roads.data(), .successorAbove = successors.data(), .overflow = &overflow};
  for (int attempt = 0; attempt < 2; ++attempt) {
    roads[0] = roads[1] = 1U;
    for (std::uint32_t node = 2U; node < kNodes; ++node) {
      std::uint64_t sum = 0U;
      const auto below = ((node / 2U) - 1U) * 2U;
      state.add(below, sum);
      state.add(below + 1U, sum);
      state.finish(node, sum);
      BOOST_TEST(roads[node] == (node < 64U ? std::uint32_t{1U} << (node / 2U) : std::numeric_limits<std::uint32_t>::max()));
    }
    BOOST_TEST(overflow == 1U);
    BOOST_TEST(std::all_of(successors.begin(), successors.end() - 2, [](auto flag) { return flag == 1U; }));
  }
  // A subsequent small sum cannot clear an earlier overflowing node's flag.
  state.finish(64U, 1U);
  BOOST_TEST(overflow == 1U);
}

BOOST_AUTO_TEST_CASE(OrderedNumberingChecksEmptyIndividualAndTotalOverflow)
{
  const std::array<double, kLayers> material{};
  const auto parameters = iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scatteringTable<Det>(material));
  const auto serial = [](std::uint32_t elements, const auto& body) {
    for (std::uint32_t element = 0U; element < elements; ++element) {
      body(element);
    }
  };
  const auto pool = std::make_shared<BoundedMemoryResource>();
  StageOutput<TrackletT> tracklets;
  tracklets.setMemoryPool(pool);
  tracklets.prepare(1U);
  const auto madeTracklets = runSerial(tracklets, [](auto, auto& emitter) { emitter.emit(Tracklet{.source = 0U, .target = 0U, .link = 0U, .lambda = 0.F, .phi = 0.F}); });
  BOOST_REQUIRE(madeTracklets == StageStatus::kComplete);
  std::array<MeasurementRange, kLayers> core{};
  core[parameters.linkSources[0]] = {.begin = 0U, .end = 1U};
  for (const std::uint32_t depth : {0U, 31U, 32U, 33U}) {
    StageOutput<TripletT> triplets;
    triplets.setMemoryPool(pool);
    triplets.prepare(2U * depth);
    const auto madeTriplets = runSerial(triplets, [](auto element, auto& emitter) { emitter.emit(Triplet{.first = 0U, .second = 0U, .curvature = 0.F, .varCurvature = 0.F, .level = element / 2U + 1U, .kinked = 0U}); });
    BOOST_REQUIRE(madeTriplets == StageStatus::kComplete);
    RoadCount counter;
    counter.setMemoryPool(pool);
    counter.prepareOrdered(serial, triplets.size());
    const auto state = counter.ordered(parameters.linkSources, core.data());
    for (std::uint32_t node = 0U; node < triplets.size(); ++node) {
      state.base(node, tracklets.view(), triplets.view());
      if (node >= 2U) {
        std::uint64_t sum = 0U;
        const auto below = ((node / 2U) - 1U) * 2U;
        state.add(below, sum);
        state.add(below + 1U, sum);
        state.finish(node, sum);
      }
    }
    const auto status = counter.runPrepared(serial, roadStartRule(parameters), tracklets, triplets, 256U);
    BOOST_TEST((status == (depth >= 32U ? StageStatus::kIndexOverflow : StageStatus::kComplete)));
    if (depth == 0U) {
      BOOST_TEST(counter.roadCount() == 0U);
      BOOST_TEST(counter.roadOffsets().size() == 1U);
    }
    // No starts under this mask: individual overflow must still be reported.
    const RoadStartRule noStarts{.minLevel = 2U, .endLayers = std::uint64_t{1U} << 7U, .linkTargets = parameters.linkTargets};
    const auto maskedStatus = counter.runPrepared(serial, noStarts, tracklets, triplets, 256U);
    BOOST_TEST((maskedStatus == (depth == 33U ? StageStatus::kIndexOverflow : StageStatus::kComplete)));
  }
}

} // namespace o2::gtrk
