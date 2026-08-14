// Copyright 2019-2026 CERN and copyright holders of ALICE O2.

#include "ITSMFTTracking/SurfaceGraphBuilder.h"

namespace o2::itsmft::tracking
{

namespace
{
struct BuiltLink {
  LinkId id{};
  SurfaceId from{};
  SurfaceId to{};
  SurfaceMask skipped{};
};
} // namespace

SurfaceGraphDefinition makeSurfaceChain(gsl::span<const SurfaceId> orderedSurfaces,
                                        int maxHoles,
                                        SurfaceMask holeSurfaces,
                                        SurfaceMask seedingSurfaces)
{
  SurfaceGraphDefinition definition;
  definition.orderedSurfaces.assign(orderedSurfaces.begin(), orderedSurfaces.end());
  definition.basePairs.reserve(orderedSurfaces.empty() ? 0 : orderedSurfaces.size() - 1);
  for (uint16_t index = 0; index + 1 < orderedSurfaces.size(); ++index) {
    definition.basePairs.push_back(SurfaceAdjacencyPair{index, static_cast<uint16_t>(index + 1)});
  }
  definition.maxHoles = maxHoles;
  definition.holeSurfaces = holeSurfaces;
  definition.seedingSurfaces = seedingSurfaces;
  return definition;
}

SurfaceGraphBuildResult SurfaceGraphBuilder::build() const
{
  SurfaceGraphBuildResult result;
  if (mCatalog.nSurfaces > MaxLayoutSurfaces) {
    result.error = SurfaceGraphBuildError::TopologyRejected;
    result.topologyError = SurfaceGraphTopologyError::InvalidSurfaceCount;
    return result;
  }
  if (mDefinition.maxHoles < 0) {
    result.error = SurfaceGraphBuildError::NegativeMaxHoles;
    return result;
  }
  SurfaceMask activeSurfaces;
  for (const auto id : mDefinition.orderedSurfaces) {
    if (!mCatalog.hasSurface(id)) {
      result.error = SurfaceGraphBuildError::InvalidSurfaceId;
      return result;
    }
    if (activeSurfaces.has(id)) {
      result.error = SurfaceGraphBuildError::DuplicateSurface;
      return result;
    }
    activeSurfaces.set(id);
  }
  if (!mDefinition.holeSurfaces.isSubsetOf(activeSurfaces)) {
    result.error = SurfaceGraphBuildError::HoleSurfacesOutsideGraph;
    return result;
  }
  if (!mDefinition.seedingSurfaces.isSubsetOf(activeSurfaces)) {
    result.error = SurfaceGraphBuildError::SeedingSurfacesOutsideGraph;
    return result;
  }

  std::vector<bool> immediate(mDefinition.orderedSurfaces.empty() ? 0 : mDefinition.orderedSurfaces.size() - 1, false);
  for (const auto pair : mDefinition.basePairs) {
    if (pair.fromIndex >= mDefinition.orderedSurfaces.size() || pair.toIndex >= mDefinition.orderedSurfaces.size() ||
        pair.toIndex != pair.fromIndex + 1) {
      result.error = SurfaceGraphBuildError::InvalidBasePair;
      return result;
    }
    if (immediate[pair.fromIndex]) {
      result.error = SurfaceGraphBuildError::DuplicateBasePair;
      return result;
    }
    immediate[pair.fromIndex] = true;
  }

  SurfaceGraph graph{gsl::span<const SurfaceDescriptor>{mCatalog.surfaces, mCatalog.nSurfaces}, mDefinition.seedingSurfaces};
  graph.setOrderedSurfaces(mDefinition.orderedSurfaces);
  std::vector<BuiltLink> links;
  for (size_t posFrom = 0; posFrom < mDefinition.orderedSurfaces.size(); ++posFrom) {
    bool connected = true;
    for (size_t posTo = posFrom + 1; posTo < mDefinition.orderedSurfaces.size(); ++posTo) {
      connected = connected && immediate[posTo - 1];
      if (!connected) {
        break;
      }
      SurfaceMask skipped;
      for (size_t k = posFrom + 1; k < posTo; ++k) {
        skipped.set(mDefinition.orderedSurfaces[k]);
      }
      if (skipped.count() > mDefinition.maxHoles || !skipped.isSubsetOf(mDefinition.holeSurfaces)) {
        continue;
      }
      const auto from = mDefinition.orderedSurfaces[posFrom];
      const auto to = mDefinition.orderedSurfaces[posTo];
      const auto id = graph.addLink(SurfaceLink{from, to, skipped, 0});
      if (!id.isValid()) {
        result.error = SurfaceGraphBuildError::TopologyRejected;
        result.topologyError = graph.getTopologyError();
        return result;
      }
      links.push_back(BuiltLink{id, from, to, skipped});
    }
  }
  for (const auto& first : links) {
    for (const auto& second : links) {
      if (first.to != second.from || (first.skipped | second.skipped).count() > mDefinition.maxHoles) {
        continue;
      }
      if (!graph.addCell(first.id, second.id).isValid()) {
        result.error = SurfaceGraphBuildError::TopologyRejected;
        result.topologyError = graph.getTopologyError();
        return result;
      }
    }
  }
  if (!graph.finalize()) {
    result.error = SurfaceGraphBuildError::TopologyRejected;
    result.topologyError = graph.getTopologyError();
    result.graphError = graph.getError();
    return result;
  }
  if (!graph.valid()) {
    result.error = SurfaceGraphBuildError::GraphRejected;
    result.topologyError = graph.getTopologyError();
    result.graphError = graph.getError();
    return result;
  }
  result.graph.emplace(std::move(graph));
  return result;
}

SurfaceGraphBatchResult buildSurfaceGraphs(SurfaceCatalogView catalog,
                                           gsl::span<const SurfaceId> orderedSurfaces,
                                           gsl::span<const o2::itsmft::TrackingParameters> trackingParameters)
{
  SurfaceGraphBatchResult result;
  result.graphs.reserve(trackingParameters.size());
  for (size_t iteration = 0; iteration < trackingParameters.size(); ++iteration) {
    const auto& parameters = trackingParameters[iteration];
    if (parameters.NLayers < 0 || static_cast<size_t>(parameters.NLayers) > orderedSurfaces.size()) {
      result.error = SurfaceGraphBuildError::GraphRejected;
      result.failedIteration = iteration;
      return result;
    }
    const auto definition = makeSurfaceChain(
      gsl::span<const SurfaceId>{orderedSurfaces.data(), static_cast<size_t>(parameters.NLayers)}, parameters.MaxHoles,
      positionalSurfaceMask(parameters.HoleLayerMask, orderedSurfaces, static_cast<uint32_t>(parameters.NLayers)),
      positionalSurfaceMask(parameters.StartLayerMask, orderedSurfaces, static_cast<uint32_t>(parameters.NLayers)));
    const auto graphResult = SurfaceGraphBuilder{catalog, definition}.build();
    if (!graphResult.ok()) {
      result.error = SurfaceGraphBuildError::GraphRejected;
      result.detail = graphResult.error;
      result.failedIteration = iteration;
      result.topologyError = graphResult.topologyError;
      result.graphError = graphResult.graphError;
      result.graphs.clear();
      return result;
    }
    result.graphs.push_back(std::move(*graphResult.graph));
  }
  return result;
}

} // namespace o2::itsmft::tracking
