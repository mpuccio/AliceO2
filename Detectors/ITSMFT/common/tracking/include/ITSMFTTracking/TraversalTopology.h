// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
#ifndef ALICEO2_ITSMFT_TRACKING_TRAVERSALTOPOLOGY_H_
#define ALICEO2_ITSMFT_TRACKING_TRAVERSALTOPOLOGY_H_

#include <cstdint>

#include <gsl/span>

#include "ITSMFTTracking/IdTypes.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMask.h"

namespace o2::itsmft::tracking
{

struct Edge {
  SurfaceId from{};
  SurfaceId to{};
  SurfaceMask skippedSurfaces{};
  uint16_t flags{0};
};

struct CellPath {
  EdgeId first{};
  EdgeId second{};
};

struct TopologyRange {
  uint32_t firstEntry{0};
  uint32_t entries{0};

  uint32_t getFirstEntry() const noexcept { return firstEntry; }
  uint32_t getEntries() const noexcept { return entries; }
  uint32_t getEntriesBound() const noexcept { return firstEntry + entries; }
};

struct TraversalTopologyView {
  SurfaceCatalogView catalog{};
  gsl::span<const SurfaceId> orderedSurfaces{};
  gsl::span<const Edge> edges{};
  gsl::span<const CellPath> paths{};
  gsl::span<const uint32_t> pathsByFirstEdgeOffsets{};
  gsl::span<const CellTopologyId> pathsByFirstEdge{};
  gsl::span<const CellTopologyId> scheduledPaths{};
  gsl::span<const CellTopologyId> roadStartPaths{};
  gsl::span<const uint32_t> roadStartComponentOffsets{};
  SurfaceMask seedingSurfaces{};

  const SurfaceDescriptor& getSurface(SurfaceId id) const { return catalog.surfaces[catalog.getSurfaceIndex(id)]; }
  const Edge& getEdge(EdgeId id) const { return edges[id.value()]; }
  const CellPath& getPath(CellTopologyId id) const { return paths[id.value()]; }
  TopologyRange getPathsStartingWithEdge(EdgeId edge) const
  {
    const auto index = edge.value();
    return {pathsByFirstEdgeOffsets[index], pathsByFirstEdgeOffsets[index + 1] - pathsByFirstEdgeOffsets[index]};
  }
};

static_assert(sizeof(CellPath) == 4);

} // namespace o2::itsmft::tracking

#endif
