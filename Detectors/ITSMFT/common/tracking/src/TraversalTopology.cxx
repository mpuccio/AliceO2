// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.

#include "ITSMFTTracking/TraversalTopology.h"
#include "ITSMFTTracking/Configuration.h"

#include <algorithm>
#include <limits>

#include <gsl/span>

namespace o2::itsmft::tracking
{

namespace
{
uint16_t componentForPosition(gsl::span<const uint16_t> componentOffsets, uint16_t position) noexcept
{
  const auto upper = std::upper_bound(componentOffsets.begin(), componentOffsets.end(), position);
  return static_cast<uint16_t>(std::distance(componentOffsets.begin(), upper) - 1);
}

LayerMask skippedBetween(uint16_t fromPosition, uint16_t toPosition) noexcept
{
  return LayerMask::skipped(fromPosition, toPosition);
}
} // namespace

TraversalTopologyBuildResult deriveTraversalTopology(const SurfaceLayout& layout,
                                                     const o2::itsmft::IterationParameters& parameters)
{
  TraversalTopologyBuildResult result;
  if (!layout.valid()) {
    result.error = TraversalTopologyError::InvalidLayout;
    return result;
  }

  const auto orderedSurfaces = layout.getOrderedSurfaces();
  if (parameters.NLayers != 0 && parameters.NLayers != orderedSurfaces.size()) {
    result.error = TraversalTopologyError::LayerCountMismatch;
    return result;
  }
  if (parameters.MaxHoles < 0) {
    result.error = TraversalTopologyError::NegativeMaxHoles;
    return result;
  }
  const auto seedingLayers = parameters.SeedingLayers;
  const auto roadStartLayers = parameters.StartLayerMask;
  const auto disabledLayers = parameters.InactiveLayerMask;

  TraversalTopology topology;
  topology.orderedSurfaces.assign(orderedSurfaces.begin(), orderedSurfaces.end());
  topology.surfacePositionById.assign(MaxLayoutSurfaces, -1);
  for (uint16_t position = 0; position < orderedSurfaces.size(); ++position) {
    topology.surfacePositionById[orderedSurfaces[position].value()] = static_cast<int16_t>(position);
    if (!disabledLayers.has(position)) {
      topology.activeLayers.set(position);
      topology.activeSurfaceList.push_back(orderedSurfaces[position]);
    }
  }
  if (topology.activeSurfaceList.empty()) {
    result.error = TraversalTopologyError::NoActiveSurfaces;
    return result;
  }
  topology.seedingLayers = seedingLayers.empty() ? topology.activeLayers : (seedingLayers & topology.activeLayers);

  const auto componentOffsets = layout.getComponentOffsets();
  const auto holeLayers = layout.getHoleLayers();
  const auto componentOf = [componentOffsets](uint16_t position) {
    return componentForPosition(componentOffsets, position);
  };

  for (uint16_t fromPosition = 0; fromPosition + 1 < orderedSurfaces.size(); ++fromPosition) {
    if (!topology.seedingLayers.has(fromPosition)) {
      continue;
    }
    for (uint16_t toPosition = fromPosition + 1; toPosition < orderedSurfaces.size(); ++toPosition) {
      if (!topology.seedingLayers.has(toPosition) ||
          componentOf(fromPosition) != componentOf(toPosition)) {
        continue;
      }
      const auto skipped = skippedBetween(fromPosition, toPosition) & topology.seedingLayers;
      if (skipped.count() > parameters.MaxHoles || !skipped.isSubsetOf(holeLayers)) {
        continue;
      }
      if (topology.edges.size() >= MaxLayoutEdges) {
        result.error = TraversalTopologyError::TooManyEdges;
        return result;
      }
      topology.edges.push_back(Edge{orderedSurfaces[fromPosition], orderedSurfaces[toPosition]});
    }
  }

  for (uint32_t first = 0; first < topology.edges.size(); ++first) {
    for (uint32_t second = 0; second < topology.edges.size(); ++second) {
      const auto& firstEdge = topology.edges[first];
      const auto& secondEdge = topology.edges[second];
      if (firstEdge.to != secondEdge.from || firstEdge.from == secondEdge.to) {
        continue;
      }
      const auto firstFrom = topology.surfacePositionById[firstEdge.from.value()];
      const auto firstTo = topology.surfacePositionById[firstEdge.to.value()];
      const auto secondFrom = topology.surfacePositionById[secondEdge.from.value()];
      const auto secondTo = topology.surfacePositionById[secondEdge.to.value()];
      const auto skipped = (skippedBetween(static_cast<uint16_t>(firstFrom), static_cast<uint16_t>(firstTo)) |
                            skippedBetween(static_cast<uint16_t>(secondFrom), static_cast<uint16_t>(secondTo))) &
                           topology.seedingLayers;
      if (skipped.count() > parameters.MaxHoles || !skipped.isSubsetOf(holeLayers)) {
        continue;
      }
      if (topology.paths.size() >= MaxLayoutPaths) {
        result.error = TraversalTopologyError::TooManyPaths;
        return result;
      }
      topology.paths.push_back(CellPath{EdgeId{static_cast<uint16_t>(first)}, EdgeId{static_cast<uint16_t>(second)}});
    }
  }

  topology.pathsByFirstEdgeOffsets.assign(topology.edges.size() + 1, 0);
  for (const auto& path : topology.paths) {
    ++topology.pathsByFirstEdgeOffsets[path.first.value() + 1];
  }
  for (size_t offset = 1; offset < topology.pathsByFirstEdgeOffsets.size(); ++offset) {
    topology.pathsByFirstEdgeOffsets[offset] += topology.pathsByFirstEdgeOffsets[offset - 1];
  }
  topology.pathsByFirstEdge.resize(topology.paths.size());
  auto cursor = topology.pathsByFirstEdgeOffsets;
  for (uint32_t path = 0; path < topology.paths.size(); ++path) {
    topology.pathsByFirstEdge[cursor[topology.paths[path].first.value()]++] = CellPathId{static_cast<uint16_t>(path)};
  }

  topology.scheduledPaths.reserve(topology.paths.size());
  for (uint32_t path = 0; path < topology.paths.size(); ++path) {
    topology.scheduledPaths.push_back(CellPathId{static_cast<uint16_t>(path)});
  }
  const auto pathOrder = [&](CellPathId lhs, CellPathId rhs) {
    const auto lhsTarget = topology.edges[topology.paths[lhs.value()].second.value()].to;
    const auto rhsTarget = topology.edges[topology.paths[rhs.value()].second.value()].to;
    const auto lhsPosition = topology.surfacePositionById[lhsTarget.value()];
    const auto rhsPosition = topology.surfacePositionById[rhsTarget.value()];
    return lhsPosition != rhsPosition ? lhsPosition < rhsPosition : lhs < rhs;
  };
  std::sort(topology.scheduledPaths.begin(), topology.scheduledPaths.end(), pathOrder);

  topology.roadStartPaths.reserve(topology.paths.size());
  for (const auto path : topology.scheduledPaths) {
    const auto target = topology.edges[topology.paths[path.value()].second.value()].to;
    if (roadStartLayers.has(topology.surfacePositionById[target.value()])) {
      topology.roadStartPaths.push_back(path);
    }
  }
  topology.roadStartComponentOffsets.push_back(0);
  uint16_t previousComponent = std::numeric_limits<uint16_t>::max();
  for (uint32_t index = 0; index < topology.roadStartPaths.size(); ++index) {
    const auto path = topology.roadStartPaths[index];
    const auto target = topology.edges[topology.paths[path.value()].second.value()].to;
    const auto component = componentOf(static_cast<uint16_t>(topology.surfacePositionById[target.value()]));
    if (component != previousComponent && index != 0) {
      topology.roadStartComponentOffsets.push_back(index);
    }
    previousComponent = component;
  }
  topology.roadStartComponentOffsets.push_back(static_cast<uint32_t>(topology.roadStartPaths.size()));

  result.topology.emplace(std::move(topology));
  return result;
}

} // namespace o2::itsmft::tracking
