// Copyright 2019-2026 CERN and copyright holders of ALICE O2.

#include "ITSMFTTracking/SurfaceGraphBuilder.h"

namespace o2::itsmft::tracking
{

namespace
{
struct BuiltTransition {
  TransitionId id{};
  SurfaceId from{};
  SurfaceId to{};
  SurfaceMask skipped{};
};
} // namespace

SurfaceGraphBuildResult SurfaceGraphBuilder::build() const
{
  SurfaceGraphBuildResult result;
  if (mCatalog.nSurfaces > MaxLayoutSurfaces) {
    result.error = SurfaceGraphBuildError::TopologyRejected;
    result.topologyError = SurfaceGraphTopologyError::InvalidSurfaceCount;
    return result;
  }
  SurfaceMask seenAcrossSubgraphs;
  SurfaceMask combinedSeeding;
  std::vector<SurfaceId> graphOrder;
  for (const auto& subgraph : mSubgraphs) {
    if (subgraph.orderedSurfaces.empty()) {
      result.error = SurfaceGraphBuildError::EmptySubgraph;
      return result;
    }
    if (subgraph.maxHoles < 0) {
      result.error = SurfaceGraphBuildError::NegativeMaxHoles;
      return result;
    }
    SurfaceMask subgraphSurfaces;
    std::optional<SurfaceKind> expectedKind;
    for (const auto id : subgraph.orderedSurfaces) {
      if (!mCatalog.hasSurface(id)) {
        result.error = SurfaceGraphBuildError::InvalidSubgraphSurfaceId;
        return result;
      }
      if (subgraphSurfaces.has(id)) {
        result.error = SurfaceGraphBuildError::DuplicateSurfaceInSubgraph;
        return result;
      }
      if (seenAcrossSubgraphs.has(id)) {
        result.error = SurfaceGraphBuildError::SurfaceDuplicatedAcrossSubgraphs;
        return result;
      }
      const auto kind = mCatalog.getSurface(id).kind;
      if (!expectedKind) {
        expectedKind = kind;
      } else if (*expectedKind != kind) {
        result.error = SurfaceGraphBuildError::GraphRejected;
        result.graphError = SurfaceGraphError::PolicySurfaceKindMismatch;
        return result;
      }
      subgraphSurfaces.set(id);
      graphOrder.push_back(id);
    }
    seenAcrossSubgraphs |= subgraphSurfaces;
    if (!subgraph.holeSurfaces.isSubsetOf(subgraphSurfaces)) {
      result.error = SurfaceGraphBuildError::HoleSurfacesOutsideSubgraph;
      return result;
    }
    if (!subgraph.seedingSurfaces.isSubsetOf(subgraphSurfaces)) {
      result.error = SurfaceGraphBuildError::SeedingSurfacesOutsideSubgraph;
      return result;
    }
    combinedSeeding |= subgraph.seedingSurfaces;
  }

  SurfaceGraph graph{gsl::span<const SurfaceDescriptor>{mCatalog.surfaces, mCatalog.nSurfaces}, combinedSeeding};
  graph.setOrderedSurfaces(std::move(graphOrder));
  for (const auto& subgraph : mSubgraphs) {
    std::vector<BuiltTransition> transitions;
    const auto& ordered = subgraph.orderedSurfaces;
    for (size_t posFrom = 0; posFrom < ordered.size(); ++posFrom) {
      for (size_t posTo = posFrom + 1; posTo < ordered.size(); ++posTo) {
        SurfaceMask skipped;
        for (size_t k = posFrom + 1; k < posTo; ++k) {
          skipped.set(ordered[k]);
        }
        if (skipped.count() > subgraph.maxHoles || !skipped.isSubsetOf(subgraph.holeSurfaces)) {
          continue;
        }
        const auto id = graph.addTransition(SurfaceTransition{ordered[posFrom], ordered[posTo], skipped, 0});
        if (!id.isValid()) {
          result.error = SurfaceGraphBuildError::TopologyRejected;
          result.topologyError = graph.getTopologyError();
          return result;
        }
        transitions.push_back(BuiltTransition{id, ordered[posFrom], ordered[posTo], skipped});
      }
    }
    for (const auto& first : transitions) {
      for (const auto& second : transitions) {
        if (first.to != second.from) {
          continue;
        }
        if ((first.skipped | second.skipped).count() > subgraph.maxHoles) {
          continue;
        }
        if (!graph.addCell(first.id, second.id).isValid()) {
          result.error = SurfaceGraphBuildError::TopologyRejected;
          result.topologyError = graph.getTopologyError();
          return result;
        }
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
    const auto activeEnd = orderedSurfaces.begin() + parameters.NLayers;
    SurfaceGraphSubgraph subgraph;
    subgraph.orderedSurfaces.assign(orderedSurfaces.begin(), activeEnd);
    subgraph.maxHoles = parameters.MaxHoles;
    subgraph.holeSurfaces = positionalSurfaceMask(parameters.HoleLayerMask, orderedSurfaces, static_cast<uint32_t>(parameters.NLayers));
    subgraph.seedingSurfaces = positionalSurfaceMask(parameters.StartLayerMask, orderedSurfaces, static_cast<uint32_t>(parameters.NLayers));
    SurfaceGraphBuilder builder{catalog};
    const auto graphResult = builder.addSubgraph(std::move(subgraph)).build();
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
