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

/// \file testSliceGlue.cxx
/// \brief The passes between a slice's stages on a device against the host
///        classes, on generated inputs.
///
/// Built once per backend. Needs the device.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK slice glue
#include <boost/test/unit_test.hpp>

#include "DeviceTestFixture.h"
#include "DeviceTestKernels.h"
#include "ToyBarrel8Fixture.h"

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/TrackletPaths.h"
#include "GenericTracking/Tracking/TripletEndIndex.h"
#include "GenericTracking/Tracking/TripletPaths.h"
#include "GenericTracking/Tracking/TripletPermutation.h"
#include "FixedPointReferences.h"
#include "SerialReferences.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace o2::gtrk
{
namespace
{

using namespace device_test;

using namespace toybarrel8;
using Sources = LayerRowMap<Det>;
constexpr std::size_t kLayerCount = Sources::kLayerCount;

constexpr std::uint32_t kLayerStride = 100000U; ///< rows of layer l start at l * kLayerStride

/// A source map of \a total rows spread over every layer, the first layers
/// taking the remainder.
Sources sourceMap(std::uint32_t total)
{
  std::array<MeasurementRange, kLayerCount> ranges{};
  std::array<bool, kLayerCount> presented{};
  for (std::size_t layer = 0; layer < kLayerCount; ++layer) {
    const auto rows = static_cast<std::uint32_t>((total / kLayerCount) + (layer < total % kLayerCount ? 1U : 0U));
    const auto begin = static_cast<std::uint32_t>(layer) * kLayerStride;
    ranges[layer] = {.begin = begin, .end = begin + rows};
    presented[layer] = true;
  }
  return Sources{ranges, presented};
}

IterationParameters<Det> parameters()
{
  const std::array<double, kLayers> material{};
  return iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scatteringTable<Det>(material));
}

} // namespace

/// Tracklets from every source to rows inside the target layer's source range
/// or just past it, interleaved, over several source counts and lengths: the
/// device's compaction with the host's path bounds is the host's permutation
/// (the path bounds included).
BOOST_AUTO_TEST_CASE(ThePermutationIsTheHostOrder)
{
  Device device;
  auto iteration = parameters();
  std::mt19937 random{2024U};
  for (const std::uint32_t sourceCount : {37U, 256U, 1000U, 1024U}) {
    const auto sources = sourceMap(sourceCount);
    constexpr std::uint32_t kTracklets = 100000U;
    const std::uint32_t perSource = kTracklets / sources.size();
    auto pool = std::make_shared<BoundedMemoryResource>();
    StageOutput<TrackletT> tracklets;
    tracklets.setMemoryPool(pool);
    tracklets.prepare(sources.size());
    const auto generate = [&](std::uint32_t element, auto& emitter) {
      const auto row = sources.row(element);
      const auto outBegin = iteration.outLinkOffsets[row.layer];
      const auto outCount = iteration.outLinkOffsets[static_cast<std::size_t>(row.layer) + 1U] - outBegin;
      if (outCount == 0U) {
        return;
      }
      for (std::uint32_t i = 0; i < perSource; ++i) {
        const auto link = static_cast<std::uint16_t>(iteration.outLinkIndices[outBegin + (static_cast<std::uint32_t>(random()) % outCount)]);
        const auto range = sources.ranges[iteration.linkTargets[link]];
        // Every other tracklet one row past its range: ineligible.
        const std::uint32_t span = std::max(range.size(), 1U);
        const std::uint32_t target = (i % 2U == 0U) ? range.begin + (static_cast<std::uint32_t>(random()) % span) : range.end + (static_cast<std::uint32_t>(random()) % 7U);
        emitter.emit(Tracklet{.source = row.row, .target = target, .link = link, .lambda = 0.F, .phi = 0.F});
      }
    };
    BOOST_REQUIRE(runSerial(tracklets, generate) == StageStatus::kComplete);
    for (const std::uint32_t length : {3U, 5U, 8U}) {
      iteration.minClusters = length;
      TrackletPaths<Det> paths;
      paths.setMemoryPool(pool);
      paths.build(iteration, sources, tracklets, SerialForEach{});
      TripletPermutation<Det> host;
      host.setMemoryPool(pool);
      BOOST_REQUIRE(host.build(iteration, sources, tracklets, paths.bounds(), SerialForEach{}, 256U) == StageStatus::kComplete);
      std::vector<std::uint32_t> onDevice;
      BOOST_REQUIRE((kernels().permutation(*device.runtime, device.stream, *device.arena, iteration, sources, tracklets, std::span{paths.bounds().prefix, sources.size()}, std::span{paths.bounds().suffix, sources.size()}, onDevice) == StageStatus::kComplete));
      BOOST_TEST_CONTEXT("sources " << sourceCount << " length " << length)
      {
        BOOST_TEST(host.size() < tracklets.size());
        BOOST_TEST(std::ranges::equal(onDevice, host.indices()));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TheEndIndexIsTheHostIndex)
{
  Device device;
  std::mt19937 random{7U};
  for (const std::uint32_t trackletCount : {1U, 64U, 1000U, 65536U}) {
    constexpr std::uint32_t kTriplets = 50000U;
    std::vector<Triplet> records(kTriplets);
    for (std::uint32_t i = 0; i < kTriplets; ++i) {
      records[i] = Triplet{
        .first = static_cast<std::uint32_t>(random()) % trackletCount,
        .second = static_cast<std::uint32_t>(random()) % trackletCount,
        .curvature = 0.F,
        .varCurvature = 0.F,
        .level = 1U,
        .kinked = 0U,
      };
    }
    const auto triplets = recordsOf<TripletT>(records);
    TripletEndIndex host;
    host.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    BOOST_REQUIRE(host.build(trackletCount, *triplets, SerialForEach{}, std::numeric_limits<std::uint32_t>::max()) == StageStatus::kComplete);
    std::vector<std::uint32_t> offsets;
    std::vector<std::uint32_t> indices;
    BOOST_REQUIRE((kernels().endIndex(*device.runtime, device.stream, *device.arena, *triplets, trackletCount, offsets, indices) == StageStatus::kComplete));
    BOOST_TEST_CONTEXT("tracklets " << trackletCount)
    {
      BOOST_TEST(std::ranges::equal(offsets, std::span<const std::uint32_t>{host.offsets(), std::size_t{trackletCount} + 1U}));
      BOOST_TEST(std::ranges::equal(indices, std::span<const std::uint32_t>{host.triplets(), host.offsets()[trackletCount]}));
    }
  }
}

BOOST_AUTO_TEST_CASE(TheChainedTripletsAreTheHostOnes)
{
  Device device;
  std::mt19937 random{11U};
  constexpr std::uint32_t kTriplets = 70000U;
  std::vector<std::vector<Neighbour>> perTriplet(kTriplets);
  for (auto& predecessors : perTriplet) {
    const auto count = static_cast<std::uint32_t>(random()) % 4U == 0U ? static_cast<std::uint32_t>(random()) % 3U + 1U : 0U;
    for (std::uint32_t i = 0; i < count; ++i) {
      predecessors.push_back(Neighbour{.from = static_cast<std::uint32_t>(random()) % kTriplets});
    }
  }
  StageOutput<NeighbourT> neighbours;
  neighbours.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  neighbours.prepare(kTriplets);
  BOOST_REQUIRE(runSerial(neighbours, [&](std::uint32_t element, auto& emitter) {
                  for (const auto& neighbour : perTriplet[element]) {
                    emitter.emit(neighbour);
                  }
                }) == StageStatus::kComplete);
  const auto host = chainedTriplets(neighbours);
  std::vector<std::uint32_t> onDevice;
  BOOST_REQUIRE((kernels().chained(*device.runtime, device.stream, *device.arena, neighbours.offsets(), onDevice) == StageStatus::kComplete));
  BOOST_TEST(host.size() > 0U);
  BOOST_TEST(std::ranges::equal(onDevice, host));
}

// Actual layered DAGs, including empty graphs, terminal rows absent from the
// source map, and irregular source ranges. Compare every bound byte with the
// independently checked CPU implementation, plus stable early compaction.
BOOST_AUTO_TEST_CASE(ThePathBoundsAndPrunedPermutationMatchTheHost)
{
  Device device;
  auto pool = std::make_shared<BoundedMemoryResource>();
  std::mt19937 random{9102026U};
  for (std::uint32_t graph = 0U; graph < 12U; ++graph) {
    auto iteration = parameters();
    std::array<MeasurementRange, kLayerCount> ranges{};
    std::array<bool, kLayerCount> presented{};
    for (std::size_t layer = 0U; layer < kLayerCount; ++layer) {
      const auto begin = static_cast<std::uint32_t>(layer) * 3U;
      ranges[layer] = {.begin = begin, .end = begin + static_cast<std::uint32_t>(random() % 4U)};
      presented[layer] = layer + 1U < kLayerCount;
    }
    const Sources sources{ranges, presented};
    StageOutput<TrackletT> tracklets;
    tracklets.setMemoryPool(pool);
    tracklets.prepare(sources.size());
    std::vector<std::vector<Tracklet>> outgoing(sources.size());
    for (std::uint32_t element = 0U; element < sources.size(); ++element) {
      const auto row = sources.row(element);
      for (std::uint32_t target = (row.layer + 1U) * 3U; target < (row.layer + 2U) * 3U; ++target) {
        if (graph != 0U && random() % 3U != 0U) {
          outgoing[element].push_back(Tracklet{.source = row.row, .target = target, .link = static_cast<std::uint16_t>(row.layer), .lambda = 0.F, .phi = 0.F});
        }
      }
    }
    BOOST_REQUIRE(runSerial(tracklets, [&](std::uint32_t element, auto& emitter) {
                    for (const auto& edge : outgoing[element]) {
                      emitter.emit(edge);
                    }
                  }) == StageStatus::kComplete);
    // The first tracklets whose middle has tracklets: at 3 clusters and any
    // end, no path is too short.
    iteration.minClusters = 3U;
    iteration.roadEndLayers = 0U;
    TrackletPaths<Det> firsts;
    firsts.setMemoryPool(pool);
    firsts.build(iteration, sources, tracklets, SerialForEach{});
    TripletPermutation<Det> permutation;
    permutation.setMemoryPool(pool);
    BOOST_REQUIRE(permutation.build(iteration, sources, tracklets, firsts.bounds(), SerialForEach{}, 256U) == StageStatus::kComplete);
    StageOutput<TripletT> triplets;
    triplets.setMemoryPool(pool);
    triplets.prepare(permutation.size());
    std::vector<std::vector<Triplet>> accepted(permutation.size());
    const auto view = tracklets.view();
    for (std::uint32_t element = 0U; element < permutation.size(); ++element) {
      const auto first = permutation.indices()[element];
      const auto middle = sources.element(iteration.linkTargets[view.link[first]], view.target[first]);
      for (auto second = tracklets.offsets()[middle]; second < tracklets.offsets()[middle + 1U]; ++second) {
        if (random() % 3U != 0U) {
          accepted[element].push_back(Triplet{.first = first, .second = second, .curvature = 0.F, .varCurvature = 0.F, .level = 1U, .kinked = 0U});
        }
      }
    }
    BOOST_REQUIRE(runSerial(triplets, [&](std::uint32_t element, auto& emitter) {
                    for (const auto& edge : accepted[element]) {
                      emitter.emit(edge);
                    }
                  }) == StageStatus::kComplete);
    TrackletPaths<Det> measurements;
    TripletPaths<Det> acceptedPaths;
    TripletPermutation<Det> pruned;
    measurements.setMemoryPool(pool);
    acceptedPaths.setMemoryPool(pool);
    pruned.setMemoryPool(pool);
    for (const std::uint64_t mask : {0ULL, 1ULL << 3U, 1ULL << 7U, (1ULL << 4U) | (1ULL << 6U)}) {
      for (const auto length : {4U, 7U, 9U}) {
        iteration.roadEndLayers = mask;
        iteration.minClusters = length;
        measurements.build(iteration, sources, tracklets, SerialForEach{});
        acceptedPaths.build(iteration, sources, tracklets, permutation, placedTriplets(triplets), SerialForEach{});
        BOOST_REQUIRE(pruned.build(iteration, sources, tracklets, measurements.bounds(), SerialForEach{}, 256U) == StageStatus::kComplete);
        PathBoundsResult result;
        BOOST_REQUIRE(kernels().pathBounds(*device.runtime, device.stream, *device.arena, iteration, sources, tracklets, triplets, permutation.indices(), result) == StageStatus::kComplete);
        BOOST_TEST_CONTEXT("graph " << graph << " mask " << mask << " length " << length)
        {
          BOOST_TEST(std::ranges::equal(result.measurementPrefix, std::span{measurements.bounds().prefix, sources.size()}));
          BOOST_TEST(std::ranges::equal(result.measurementSuffix, std::span{measurements.bounds().suffix, sources.size()}));
          BOOST_TEST(std::ranges::equal(result.tripletPrefix, std::span{acceptedPaths.bounds().prefix, tracklets.size()}));
          BOOST_TEST(std::ranges::equal(result.tripletSuffix, std::span{acceptedPaths.bounds().suffix, tracklets.size()}));
          BOOST_TEST(std::ranges::equal(result.permutation, pruned.indices()));
        }
      }
    }
  }
  BOOST_TEST(device.arena->used(ArenaEnd::kScratch) == 0U);
}

} // namespace o2::gtrk
