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

/// \file testRoadGlue.cxx
/// \brief A slice's road counts, starts and offsets on a device against
///        RoadCount on generated triplet graphs.
///
/// Built once per backend. Needs the device.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK road glue
#include <boost/test/unit_test.hpp>

#include "DeviceTestFixture.h"
#include "DeviceTestKernels.h"
#include "ToyBarrel8Fixture.h"

#include "GenericTracking/Tracking/RoadCount.h"
#include "SerialReferences.h"
#include "FixedPointReferences.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <span>
#include <vector>

namespace o2::gtrk
{
namespace
{

using namespace device_test;
using namespace toybarrel8;
constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);
constexpr std::uint32_t kRowsPerLayer = 1000U;

const auto kSerial = [](std::uint32_t elements, const auto& body) {
  for (std::uint32_t element = 0; element < elements; ++element) {
    body(element);
  }
};

IterationParameters<Det> parameters()
{
  const std::array<double, kLayers> material{};
  return iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scatteringTable<Det>(material));
}

/// The first half of every layer's rows is the slice's core.
std::array<MeasurementRange, kLayerCount> core()
{
  std::array<MeasurementRange, kLayerCount> result{};
  for (std::size_t layer = 0; layer < kLayerCount; ++layer) {
    const auto begin = static_cast<std::uint32_t>(layer) * kRowsPerLayer;
    result[layer] = {.begin = begin, .end = begin + (kRowsPerLayer / 2U)};
  }
  return result;
}

/// A triplet graph: the tracklets, the triplets on them, each triplet's
/// predecessors -- all earlier in the arena, so the graph is acyclic -- and
/// the levels the host's sweep gives them.
struct Graph {
  std::unique_ptr<StageOutput<TrackletT>> tracklets;
  std::unique_ptr<StageOutput<TripletT>> triplets;
  std::unique_ptr<StageOutput<NeighbourT>> neighbours;
  std::vector<std::uint32_t> chained;
};

Graph graph(const IterationParameters<Det>& iteration, const std::vector<Tracklet>& tracklets, const std::vector<std::uint32_t>& firstTracklet, const std::vector<std::vector<std::uint32_t>>& predecessors)
{
  Graph result;
  result.tracklets = recordsOf<TrackletT>(tracklets);
  std::vector<Triplet> triplets(firstTracklet.size());
  for (std::size_t i = 0; i < triplets.size(); ++i) {
    triplets[i] = Triplet{.first = firstTracklet[i], .second = firstTracklet[i], .curvature = 0.F, .varCurvature = 0.F, .level = 1U, .kinked = 0U};
  }
  result.triplets = recordsOf<TripletT>(triplets);
  result.neighbours = std::make_unique<StageOutput<NeighbourT>>();
  result.neighbours->setMemoryPool(std::make_shared<BoundedMemoryResource>());
  result.neighbours->prepare(static_cast<std::uint32_t>(predecessors.size()));
  BOOST_REQUIRE(runSerial(*result.neighbours, [&](std::uint32_t element, auto& emitter) {
                  for (const auto from : predecessors[element]) {
                    emitter.emit(Neighbour{.from = from});
                  }
                }) == StageStatus::kComplete);
  result.chained = chainedTriplets(*result.neighbours);
  static_cast<void>(fixedPointLevels(kSerial, *result.neighbours, *result.triplets));
  static_cast<void>(iteration);
  return result;
}

/// The host's RoadCount and the device's glue on \a graph: the same status,
/// and on success the same flags, counts, starts and offsets.
StageStatus compare(const IterationParameters<Det>& iteration, const Graph& input)
{
  const auto cores = core();
  RoadCount host;
  host.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  const auto hostStatus = fixedPointRoads(host, kSerial, iteration.linkSources, cores.data(), roadStartRule(iteration), *input.tracklets, *input.triplets, *input.neighbours);
  BOOST_TEST_MESSAGE("host: " << host.starts().size() << " starts, " << host.roadCount() << " roads, status " << static_cast<int>(hostStatus));
  Device device;
  RoadGlueResult ordered;
  const auto orderedStatus = kernels().orderedRoadGlue(*device.runtime, device.stream, *device.arena, iteration, *input.tracklets, *input.triplets, *input.neighbours, cores, ordered);
  BOOST_TEST((orderedStatus == hostStatus));
  BOOST_TEST(ordered.chained == input.chained.size());
  if (hostStatus == StageStatus::kComplete && orderedStatus == StageStatus::kComplete) {
    BOOST_TEST(std::ranges::equal(ordered.successorAbove, host.successorAbove()));
    BOOST_TEST(std::ranges::equal(ordered.roads, host.roads()));
    BOOST_TEST(std::ranges::equal(ordered.starts, host.starts()));
    BOOST_TEST(std::ranges::equal(ordered.offsets, host.roadOffsets()));
  }
  BOOST_TEST(device.arena->used() == 0U);
  return hostStatus;
}

/// A ladder of \a levels levels, two triplets each, every one a predecessor of
/// both above it: the level-k triplets have 2^(k-1) roads each.
Graph ladder(const IterationParameters<Det>& iteration, std::uint32_t levels)
{
  const auto coreRow = static_cast<std::uint32_t>(iteration.linkSources[0]) * kRowsPerLayer;
  const std::size_t tripletCount = std::size_t{2U} * levels;
  std::vector<Tracklet> tracklets(tripletCount);
  std::vector<std::uint32_t> first(tripletCount);
  std::vector<std::vector<std::uint32_t>> predecessors(tripletCount);
  for (std::uint32_t triplet = 0; triplet < 2U * levels; ++triplet) {
    tracklets[triplet] = Tracklet{.source = coreRow, .target = coreRow, .link = 0U, .lambda = 0.F, .phi = 0.F};
    first[triplet] = triplet;
    if (triplet >= 2U) {
      const std::uint32_t below = ((triplet / 2U) - 1U) * 2U;
      predecessors[triplet] = {below, below + 1U};
    }
  }
  return graph(iteration, tracklets, first, predecessors);
}

} // namespace

/// Random acyclic graphs: tracklets inside and outside the core, triplets
/// with up to three predecessors.
BOOST_AUTO_TEST_CASE(RandomGraphsGiveTheHostsRoads)
{
  const auto iteration = parameters();
  std::mt19937 random{21U};
  for (const std::uint32_t tripletCount : {0U, 1U, 50U, 3000U, 40000U}) {
    constexpr std::uint32_t kTracklets = 5000U;
    std::vector<Tracklet> tracklets(kTracklets);
    for (auto& tracklet : tracklets) {
      const auto link = static_cast<std::uint16_t>(static_cast<std::uint32_t>(random()) % static_cast<std::uint32_t>(nLinks<Det>));
      const auto layer = static_cast<std::uint32_t>(iteration.linkSources[link]);
      const auto source = (layer * kRowsPerLayer) + (static_cast<std::uint32_t>(random()) % kRowsPerLayer);
      tracklet = Tracklet{.source = source, .target = source, .link = link, .lambda = 0.F, .phi = 0.F};
    }
    std::vector<std::uint32_t> first(tripletCount);
    std::vector<std::vector<std::uint32_t>> predecessors(tripletCount);
    for (std::uint32_t triplet = 0; triplet < tripletCount; ++triplet) {
      first[triplet] = static_cast<std::uint32_t>(random()) % kTracklets;
      const auto count = triplet == 0U ? 0U : static_cast<std::uint32_t>(random()) % 4U;
      for (std::uint32_t i = 0; i < count; ++i) {
        predecessors[triplet].push_back(static_cast<std::uint32_t>(random()) % triplet);
      }
    }
    const auto input = graph(iteration, tracklets, first, predecessors);
    BOOST_TEST_CONTEXT("triplets " << tripletCount) { BOOST_TEST((compare(iteration, input) == StageStatus::kComplete)); }
  }
}

/// Counts past 32 bits: a triplet's is kIndexOverflow, and so is a total of
/// the starts' that is, each count still fitting.
BOOST_AUTO_TEST_CASE(OverflowsAreStatusesAsOnTheHost)
{
  const auto iteration = parameters();
  BOOST_TEST((compare(iteration, ladder(iteration, 31U)) == StageStatus::kComplete));      // two starts of 2^30
  BOOST_TEST((compare(iteration, ladder(iteration, 32U)) == StageStatus::kIndexOverflow)); // two of 2^31: the total
  BOOST_TEST((compare(iteration, ladder(iteration, 33U)) == StageStatus::kIndexOverflow)); // 2^32 each
  auto masked = iteration;
  masked.roadEndLayers = std::uint64_t{1U} << 7U;
  BOOST_TEST((compare(masked, ladder(masked, 32U)) == StageStatus::kComplete));      // excluded starts: no total overflow
  BOOST_TEST((compare(masked, ladder(masked, 33U)) == StageStatus::kIndexOverflow)); // individual overflow remains sticky
}

/// Only level-1 triplets, below the start level: no roads, one offset.
BOOST_AUTO_TEST_CASE(NoStartsGiveOneOffset)
{
  const auto iteration = parameters();
  BOOST_REQUIRE(iteration.minRoadStartLevel > 1U);
  const auto coreRow = static_cast<std::uint32_t>(iteration.linkSources[0]) * kRowsPerLayer;
  const std::vector<Tracklet> tracklets(4U, Tracklet{.source = coreRow, .target = coreRow, .link = 0U, .lambda = 0.F, .phi = 0.F});
  const auto input = graph(iteration, tracklets, {0U, 1U, 2U, 3U}, std::vector<std::vector<std::uint32_t>>(4U));
  BOOST_TEST((compare(iteration, input) == StageStatus::kComplete));
  RoadCount host;
  host.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  const auto cores = core();
  BOOST_REQUIRE(fixedPointRoads(host, kSerial, iteration.linkSources, cores.data(), roadStartRule(iteration), *input.tracklets, *input.triplets, *input.neighbours) == StageStatus::kComplete);
  BOOST_TEST(host.starts().empty());
  BOOST_TEST(host.roadOffsets().size() == 1U);
}

} // namespace o2::gtrk
