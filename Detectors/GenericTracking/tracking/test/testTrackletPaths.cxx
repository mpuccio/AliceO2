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

/// \file testTrackletPaths.cxx
/// \brief The path bounds against every path enumerated, serial against
///        parallel, and the sweeps in the layout's order, not the layers'.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TrackletPaths
#include <boost/test/unit_test.hpp>

#include "ToyBarrel8.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/Neighbour.h"
#include "GenericTracking/Tracking/TrackletPaths.h"
#include "GenericTracking/Tracking/TripletPaths.h"
#include "SerialReferences.h"

#include <array>
#include <cstdint>
#include <random>
#include <vector>

namespace o2::gtrk
{
namespace
{
using Det = detectors::ToyBarrel8;
using Index = StageOutput<TrackletT>::Index;
constexpr Index kRows = 24U;
constexpr Index kPerLayer = 3U;

// Enumerate every path, without dynamic programming, and mark every triplet
// on paths of the requested length terminating on an allowed end layer.
void checkGraph(std::mt19937& random)
{
  std::array<MeasurementRange, 8> ranges{};
  std::array<bool, 8> presented{};
  for (Index layer = 0U; layer < 8U; ++layer) {
    ranges[layer] = {.begin = layer * kPerLayer, .end = (layer + 1U) * kPerLayer};
    presented[layer] = layer < 7U;
  }
  LayerRowMap<Det> sources{ranges, presented};
  std::array<std::vector<Index>, kRows> outgoing;
  for (Index source = 0U; source < 21U; ++source) {
    for (Index target = (source / kPerLayer + 1U) * kPerLayer; target < (source / kPerLayer + 2U) * kPerLayer; ++target) {
      if (random() % 4U == 0U) {
        outgoing[source].push_back(target);
      }
    }
  }
  StageOutput<TrackletT> tracklets;
  tracklets.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  tracklets.prepare(sources.size());
  BOOST_REQUIRE(runSerial(tracklets, [&](Index element, auto& emit) {
                  const auto row = sources.row(element);
                  for (const auto target : outgoing[row.row]) {
                    emit.emit(Tracklet{.source = row.row, .target = target, .link = static_cast<std::uint16_t>(row.layer), .lambda = 0.F, .phi = 0.F});
                  }
                }) == StageStatus::kComplete);
  TrackletPaths<Det> paths;
  paths.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  for (std::uint32_t length = 4U; length <= 9U; ++length) {
    for (const std::uint64_t mask : {std::uint64_t{0U}, std::uint64_t{1U << 3U}, std::uint64_t{1U << 7U}, std::uint64_t{(1U << 4U) | (1U << 6U)}}) {
      IterationParameters<Det> parameters;
      parameters.minClusters = length;
      parameters.roadEndLayers = mask;
      for (std::size_t link = 0U; link < Det::layerLinks.size(); ++link) {
        parameters.linkSources[link] = Det::layerLinks[link].from;
        parameters.linkTargets[link] = Det::layerLinks[link].to;
      }
      paths.build(parameters, sources, tracklets, SerialForEach{});
      const auto bounds = paths.bounds();
      if (length == 7U) {
        const std::vector<std::uint8_t> prefix(bounds.prefix, bounds.prefix + sources.size());
        const std::vector<std::uint8_t> suffix(bounds.suffix, bounds.suffix + sources.size());
        detail::TaskArena cpuArena{8};
        CpuStage<NeighbourT> cpu{cpuArena};
        cpu.setMemoryPool(std::make_shared<BoundedMemoryResource>());
        const auto parallel = [&](Index count, const auto& body) { cpu.forEachElement(count, 1U, body); };
        paths.build(parameters, sources, tracklets, parallel);
        for (Index element = 0U; element < sources.size(); ++element) {
          BOOST_TEST(paths.bounds().prefix[element] == prefix[element]);
          BOOST_TEST(paths.bounds().suffix[element] == suffix[element]);
        }
      }
      std::array<bool, std::size_t{kRows} * kRows * kRows> kept{};
      std::vector<Index> chain;
      const auto visit = [&](const auto& self, Index row) -> void { // NOLINT(misc-no-recursion)
        chain.push_back(row);
        if (chain.size() >= length && (mask == 0U || ((mask >> (row / kPerLayer)) & 1U) != 0U)) {
          for (std::size_t i = 0U; i + 2U < chain.size(); ++i) {
            kept[(chain[i] * kRows + chain[i + 1U]) * kRows + chain[i + 2U]] = true;
          }
        }
        for (const auto target : outgoing[row]) {
          self(self, target);
        }
        chain.pop_back();
      };
      for (Index row = 0U; row < kRows; ++row) {
        visit(visit, row);
      }
      for (Index a = 0U; a < kRows; ++a) {
        for (const Index b : outgoing[a]) {
          for (const Index c : outgoing[b]) {
            // a has tracklets, so it is a source; c may be on the last layer.
            const std::uint32_t prefix = bounds.prefix[sources.element(static_cast<LayerID>(a / kPerLayer), a)];
            const auto suffix = measurementSuffix(bounds, sources, mask, static_cast<LayerID>(c / kPerLayer), c);
            const bool possible = suffix != 0U && prefix + 1U + suffix >= length;
            BOOST_TEST(possible == kept[(a * kRows + b) * kRows + c]);
          }
        }
      }
    }
  }
}
} // namespace
BOOST_AUTO_TEST_CASE(ActualPathsMatchExhaustiveEnumeration)
{
  std::mt19937 random{47U};
  for (unsigned graph = 0U; graph < 100U; ++graph) {
    checkGraph(random);
  }
}
} // namespace o2::gtrk

namespace o2::gtrk
{
BOOST_AUTO_TEST_CASE(TripletPathsMatchExhaustiveOverlapPaths)
{
  std::mt19937 random{4711U};
  for (unsigned graph = 0U; graph < 50U; ++graph) {
    std::array<MeasurementRange, 8> ranges{};
    std::array<bool, 8> presented{};
    for (Index layer = 0U; layer < 8U; ++layer) {
      ranges[layer] = {.begin = layer * kPerLayer, .end = (layer + 1U) * kPerLayer};
      presented[layer] = layer < 7U;
    }
    const LayerRowMap<Det> sources{ranges, presented};
    IterationParameters<Det> parameters;
    for (std::size_t link = 0U; link < Det::layerLinks.size(); ++link) {
      parameters.linkSources[link] = Det::layerLinks[link].from;
      parameters.linkTargets[link] = Det::layerLinks[link].to;
    }
    StageOutput<TrackletT> tracklets;
    tracklets.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    tracklets.prepare(sources.size());
    const auto generateTracklets = [&](Index element, auto& emit) {
      const auto row = sources.row(element);
      for (Index target = (std::uint32_t{row.layer} + 1U) * kPerLayer; target < (std::uint32_t{row.layer} + 2U) * kPerLayer; ++target) {
        if (random() % 3U != 0U) {
          emit.emit(Tracklet{.source = row.row, .target = target, .link = static_cast<std::uint16_t>(row.layer), .lambda = 0.F, .phi = 0.F});
        }
      }
    };
    BOOST_REQUIRE(runSerial(tracklets, generateTracklets) == StageStatus::kComplete);
    // The first tracklets whose middle has tracklets: at 3 clusters and any
    // end, no path is too short.
    parameters.minClusters = 3U;
    TrackletPaths<Det> firsts;
    firsts.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    firsts.build(parameters, sources, tracklets, SerialForEach{});
    TripletPermutation<Det> permutation;
    permutation.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    BOOST_REQUIRE(permutation.build(parameters, sources, tracklets, firsts.bounds(), SerialForEach{}, 256U) == StageStatus::kComplete);
    const auto tv = tracklets.view();
    std::vector<std::vector<Index>> outgoing(tracklets.size());
    StageOutput<TripletT> triplets;
    triplets.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    triplets.prepare(permutation.size());
    const auto generateTriplets = [&](Index element, auto& emit) {
      const auto first = permutation.indices()[element];
      for (Index second = 0U; second < tracklets.size(); ++second) {
        if (tv.target[first] == tv.source[second] && random() % 3U != 0U) {
          outgoing[first].push_back(second);
          emit.emit(Triplet{.first = first, .second = second, .curvature = 0.F, .varCurvature = 1.F, .level = 1U, .kinked = 0U});
        }
      }
    };
    BOOST_REQUIRE(runSerial(triplets, generateTriplets) == StageStatus::kComplete);
    TripletPaths<Det> paths;
    paths.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    const std::size_t n = tracklets.size();
    for (std::uint32_t length = 4U; length <= 9U; ++length) {
      for (const std::uint64_t mask : {std::uint64_t{0U}, std::uint64_t{1U << 3U}, std::uint64_t{1U << 7U}, std::uint64_t{(1U << 4U) | (1U << 6U)}}) {
        parameters.minClusters = length;
        parameters.roadEndLayers = mask;
        paths.build(parameters, sources, tracklets, permutation, placedTriplets(triplets), SerialForEach{});
        const auto bounds = paths.bounds();
        if (length == 7U) {
          const std::vector<std::uint8_t> prefix(bounds.prefix, bounds.prefix + n);
          const std::vector<std::uint8_t> suffix(bounds.suffix, bounds.suffix + n);
          detail::TaskArena cpuArena{8};
          CpuStage<NeighbourT> cpu{cpuArena};
          cpu.setMemoryPool(std::make_shared<BoundedMemoryResource>());
          const auto parallel = [&](Index count, const auto& body) { cpu.forEachElement(count, 1U, body); };
          paths.build(parameters, sources, tracklets, permutation, placedTriplets(triplets), parallel);
          for (Index row = 0U; row < tracklets.size(); ++row) {
            BOOST_TEST(paths.bounds().prefix[row] == prefix[row]);
            BOOST_TEST(paths.bounds().suffix[row] == suffix[row]);
          }
        }
        std::vector<bool> kept(n * n * n, false);
        std::vector<Index> chain;
        const auto visit = [&](const auto& self, Index row) -> void { // NOLINT(misc-no-recursion)
          chain.push_back(row);
          const auto endLayer = parameters.linkTargets[tv.link[row]];
          if (chain.size() + 1U >= length && (mask == 0U || ((mask >> endLayer) & 1U) != 0U)) {
            for (std::size_t i = 0U; i + 2U < chain.size(); ++i) {
              kept[(chain[i] * n + chain[i + 1U]) * n + chain[i + 2U]] = true;
            }
          }
          for (const auto target : outgoing[row]) {
            self(self, target);
          }
          chain.pop_back();
        };
        for (Index row = 0U; row < tracklets.size(); ++row) {
          visit(visit, row);
        }
        for (Index a = 0U; a < tracklets.size(); ++a) {
          for (const Index b : outgoing[a]) {
            for (const Index c : outgoing[b]) {
              const bool possible = bounds.suffix[c] != 0U && std::uint32_t{bounds.prefix[a]} + 1U + bounds.suffix[c] >= length;
              BOOST_TEST(possible == kept[(a * n + b) * n + c]);
            }
          }
        }
      }
    }
  }
}
} // namespace o2::gtrk

namespace o2::gtrk
{
namespace
{
// A numeric layer order deliberately opposite to the topological order.
struct ReverseBarrel8 : detectors::ToyBarrel8 {
  static constexpr std::array<LayerLink, 7> layerLinks{{{.from = 7, .to = 6, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false}, {.from = 6, .to = 5, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false}, {.from = 5, .to = 4, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false}, {.from = 4, .to = 3, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false}, {.from = 3, .to = 2, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false}, {.from = 2, .to = 1, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false}, {.from = 1, .to = 0, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false}}};
};
} // namespace
BOOST_AUTO_TEST_CASE(PathSweepsFollowTopologyInsteadOfNumericLayerOrder)
{
  std::array<MeasurementRange, 8> ranges{};
  std::array<bool, 8> presented{};
  for (Index layer = 0U; layer < 8U; ++layer) {
    ranges[layer] = {.begin = layer, .end = layer + 1U};
    presented[layer] = layer != 0U;
  }
  const LayerRowMap<ReverseBarrel8> sources{ranges, presented};
  IterationParameters<ReverseBarrel8> parameters;
  parameters.roadEndLayers = 1U; // The end of the chain is layer 0.
  parameters.minClusters = 8U;
  for (std::size_t link = 0U; link < ReverseBarrel8::layerLinks.size(); ++link) {
    parameters.linkSources[link] = ReverseBarrel8::layerLinks[link].from;
    parameters.linkTargets[link] = ReverseBarrel8::layerLinks[link].to;
  }
  StageOutput<TrackletT> tracklets;
  tracklets.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  tracklets.prepare(sources.size());
  BOOST_REQUIRE(runSerial(tracklets, [&](Index element, auto& emit) {
                  const auto source = sources.row(element).row;
                  emit.emit(Tracklet{.source = source, .target = source - 1U, .link = static_cast<std::uint16_t>(7U - source), .lambda = 0.F, .phi = 0.F});
                }) == StageStatus::kComplete);
  TrackletPaths<ReverseBarrel8> paths;
  paths.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  paths.build(parameters, sources, tracklets, SerialForEach{});
  // Row 0, the chain's end, has no tracklets and is no source.
  BOOST_TEST(paths.bounds().prefix[sources.element(1, 1U)] == 7U);
  BOOST_TEST(paths.bounds().suffix[sources.element(7, 7U)] == 8U);
  BOOST_TEST(measurementSuffix(paths.bounds(), sources, parameters.roadEndLayers, 0, 0U) == 1U);
  TripletPermutation<ReverseBarrel8> permutation;
  permutation.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  BOOST_REQUIRE(permutation.build(parameters, sources, tracklets, paths.bounds(), SerialForEach{}, 256U) == StageStatus::kComplete);
  BOOST_TEST(permutation.size() == 6U);
  StageOutput<TripletT> triplets;
  triplets.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  triplets.prepare(permutation.size());
  BOOST_REQUIRE(runSerial(triplets, [&](Index element, auto& emit) {
                  const auto first = permutation.indices()[element];
                  emit.emit(Triplet{.first = first, .second = first - 1U, .curvature = 0.F, .varCurvature = 1.F, .level = 1U, .kinked = 0U});
                }) == StageStatus::kComplete);
  TripletPaths<ReverseBarrel8> tripletPaths;
  tripletPaths.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  tripletPaths.build(parameters, sources, tracklets, permutation, placedTriplets(triplets), SerialForEach{});
  BOOST_TEST(tripletPaths.bounds().prefix[0] == 8U);
  BOOST_TEST(tripletPaths.bounds().suffix[6] == 7U);
}
} // namespace o2::gtrk
