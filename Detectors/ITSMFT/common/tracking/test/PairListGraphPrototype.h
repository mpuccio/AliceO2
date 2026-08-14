// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_TEST_PAIRLISTGRAPHPROTOTYPE_H_
#define ALICEO2_ITSMFT_TRACKING_TEST_PAIRLISTGRAPHPROTOTYPE_H_

#include <algorithm>
#include <cstdint>
#include <optional>
#include <type_traits>
#include <vector>

#include "ITSMFTTracking/SurfaceGraph.h"

namespace o2::itsmft::tracking::test
{

// This record is intentionally the only representation of an immediate edge
// in the prototype. It is not a production graph authority.
struct PairListBasePair {
  uint16_t fromIndex{0};
  uint16_t toIndex{0};
};

static_assert(std::is_standard_layout_v<PairListBasePair> && std::is_trivially_copyable_v<PairListBasePair>);
static_assert(sizeof(PairListBasePair) == 4);

struct PairListComponentInput {
  std::vector<SurfaceId> activeSurfaces;
  std::vector<PairListBasePair> basePairs;
};

struct PairListHolePolicy {
  int maxSkipped{0};
  SurfaceMask skippableSurfaceMask{};
};

struct PairListGraphInput {
  std::vector<PairListComponentInput> components;
  PairListHolePolicy holePolicy{};
  SurfaceMask seedingMask{};
};

enum class PairListGraphError : uint8_t {
  None,
  EmptyComponent,
  NegativeHoleBudget,
  InvalidSurfaceId,
  DuplicateActiveSurface,
  HoleMaskOutsideActive,
  SeedingMaskOutsideActive,
  InvalidPair,
  DuplicatePair,
};

struct PairListComponentSequences {
  std::vector<LinkId> links;
  std::vector<CellTopologyId> cells;
  std::vector<CellTopologyId> scheduledCells;
  std::vector<CellTopologyId> roadStartCells;
};

struct PairListGraph {
  std::vector<SurfaceLink> links;
  // Kept separately to make skipped-surface witnesses an explicit parity
  // output, even though each witness is also in SurfaceLink.
  std::vector<SurfaceMask> skippedWitnesses;
  std::vector<SurfaceCellTopology> cells;
  std::vector<uint32_t> cellsByFirstLinkOffsets;
  std::vector<CellTopologyId> cellsByFirstLink;
  std::vector<CellTopologyId> scheduledCells;
  std::vector<CellTopologyId> roadStartCells;
  std::vector<PairListComponentSequences> components;
};

struct PairListGraphBuildResult {
  std::optional<PairListGraph> graph;
  PairListGraphError error{PairListGraphError::None};

  bool ok() const noexcept { return graph.has_value(); }
};

inline PairListGraphBuildResult derivePairListGraph(const PairListGraphInput& input)
{
  PairListGraphBuildResult result;
  PairListGraph graph;
  graph.components.resize(input.components.size());

  SurfaceMask activeAcrossComponents;
  std::vector<std::vector<bool>> immediate(input.components.size());
  std::vector<SurfaceMask> componentHoleMasks(input.components.size());
  std::vector<SurfaceMask> componentSeedingMasks(input.components.size());
  for (size_t componentIndex = 0; componentIndex < input.components.size(); ++componentIndex) {
    const auto& component = input.components[componentIndex];
    if (component.activeSurfaces.empty()) {
      result.error = PairListGraphError::EmptyComponent;
      return result;
    }
    if (input.holePolicy.maxSkipped < 0) {
      result.error = PairListGraphError::NegativeHoleBudget;
      return result;
    }
    SurfaceMask active;
    for (const auto surface : component.activeSurfaces) {
      if (!surface.isValid() || surface.value() >= MaxLayoutSurfaces) {
        result.error = PairListGraphError::InvalidSurfaceId;
        return result;
      }
      if (active.has(surface) || activeAcrossComponents.has(surface)) {
        result.error = PairListGraphError::DuplicateActiveSurface;
        return result;
      }
      active.set(surface);
      activeAcrossComponents.set(surface);
    }
    componentHoleMasks[componentIndex] = input.holePolicy.skippableSurfaceMask & active;
    componentSeedingMasks[componentIndex] = input.seedingMask & active;
    immediate[componentIndex].assign(component.activeSurfaces.size() - 1, false);
    for (const auto pair : component.basePairs) {
      if (pair.fromIndex >= component.activeSurfaces.size() || pair.toIndex >= component.activeSurfaces.size() ||
          pair.fromIndex >= pair.toIndex || pair.toIndex != pair.fromIndex + 1) {
        result.error = PairListGraphError::InvalidPair;
        return result;
      }
      if (immediate[componentIndex][pair.fromIndex]) {
        result.error = PairListGraphError::DuplicatePair;
        return result;
      }
      immediate[componentIndex][pair.fromIndex] = true;
    }
  }
  if (!input.holePolicy.skippableSurfaceMask.isSubsetOf(activeAcrossComponents)) {
    result.error = PairListGraphError::HoleMaskOutsideActive;
    return result;
  }
  if (!input.seedingMask.isSubsetOf(activeAcrossComponents)) {
    result.error = PairListGraphError::SeedingMaskOutsideActive;
    return result;
  }

  struct LinkInfo {
    size_t component{0};
    uint16_t fromRank{0};
    uint16_t toRank{0};
  };
  std::vector<LinkInfo> linkInfo;
  for (size_t componentIndex = 0; componentIndex < input.components.size(); ++componentIndex) {
    const auto& component = input.components[componentIndex];
    for (uint16_t fromRank = 0; fromRank + 1 < component.activeSurfaces.size(); ++fromRank) {
      bool allEdgesPresent = true;
      for (uint16_t toRank = fromRank + 1; toRank < component.activeSurfaces.size(); ++toRank) {
        allEdgesPresent = allEdgesPresent && immediate[componentIndex][toRank - 1];
        if (!allEdgesPresent) {
          break;
        }
        SurfaceMask skipped;
        for (uint16_t skippedRank = fromRank + 1; skippedRank < toRank; ++skippedRank) {
          skipped.set(component.activeSurfaces[skippedRank]);
        }
        if (skipped.count() > input.holePolicy.maxSkipped ||
            !skipped.isSubsetOf(componentHoleMasks[componentIndex])) {
          continue;
        }
        graph.links.push_back(SurfaceLink{component.activeSurfaces[fromRank],
                                                      component.activeSurfaces[toRank], skipped, 0});
        graph.skippedWitnesses.push_back(skipped);
        graph.components[componentIndex].links.push_back(
          LinkId{static_cast<uint16_t>(graph.links.size() - 1)});
        linkInfo.push_back(LinkInfo{componentIndex, fromRank, toRank});
      }
    }
  }

  struct CellInfo {
    size_t component{0};
    uint16_t targetRank{0};
  };
  std::vector<CellInfo> cellInfo;
  for (size_t first = 0; first < graph.links.size(); ++first) {
    for (size_t second = 0; second < graph.links.size(); ++second) {
      const auto& firstLink = graph.links[first];
      const auto& secondLink = graph.links[second];
      if (firstLink.to != secondLink.from ||
          (firstLink.skippedSurfaces | secondLink.skippedSurfaces).count() >
            input.holePolicy.maxSkipped) {
        continue;
      }
      const auto componentIndex = linkInfo[first].component;
      if (componentIndex != linkInfo[second].component || firstLink.from == secondLink.to) {
        continue;
      }
      const auto id = CellTopologyId{static_cast<uint16_t>(graph.cells.size())};
      SurfaceMask hits;
      hits.set(firstLink.from);
      hits.set(firstLink.to);
      hits.set(secondLink.to);
      graph.cells.push_back(SurfaceCellTopology{LinkId{static_cast<uint16_t>(first)},
                                                LinkId{static_cast<uint16_t>(second)}, hits});
      cellInfo.push_back(CellInfo{componentIndex, linkInfo[second].toRank});
      graph.components[componentIndex].cells.push_back(id);
      if (componentSeedingMasks[componentIndex].has(secondLink.to)) {
        graph.roadStartCells.push_back(id);
        graph.components[componentIndex].roadStartCells.push_back(id);
      }
    }
  }

  graph.cellsByFirstLinkOffsets.assign(graph.links.size() + 1, 0);
  for (const auto& cell : graph.cells) {
    ++graph.cellsByFirstLinkOffsets[cell.firstLink.value() + 1];
  }
  for (size_t i = 1; i < graph.cellsByFirstLinkOffsets.size(); ++i) {
    graph.cellsByFirstLinkOffsets[i] += graph.cellsByFirstLinkOffsets[i - 1];
  }
  graph.cellsByFirstLink.resize(graph.cells.size());
  auto cursor = graph.cellsByFirstLinkOffsets;
  for (size_t id = 0; id < graph.cells.size(); ++id) {
    graph.cellsByFirstLink[cursor[graph.cells[id].firstLink.value()]++] = CellTopologyId{static_cast<uint16_t>(id)};
  }

  std::vector<CellTopologyId> sorted;
  sorted.reserve(graph.cells.size());
  for (uint16_t id = 0; id < graph.cells.size(); ++id) {
    sorted.push_back(CellTopologyId{id});
  }
  std::sort(sorted.begin(), sorted.end(), [&](CellTopologyId lhs, CellTopologyId rhs) {
    const auto& left = cellInfo[lhs.value()];
    const auto& right = cellInfo[rhs.value()];
    if (left.component != right.component) {
      return left.component < right.component;
    }
    return left.targetRank != right.targetRank ? left.targetRank < right.targetRank : lhs < rhs;
  });
  graph.scheduledCells = sorted;
  for (const auto id : sorted) {
    graph.components[cellInfo[id.value()].component].scheduledCells.push_back(id);
  }

  result.graph.emplace(std::move(graph));
  return result;
}

} // namespace o2::itsmft::tracking::test

#endif
