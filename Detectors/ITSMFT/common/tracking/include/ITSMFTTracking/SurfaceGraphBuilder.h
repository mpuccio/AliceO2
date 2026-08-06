// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEGRAPHBUILDER_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEGRAPHBUILDER_H_

#include <cstddef>
#include <cstdint>

#ifndef GPUCA_GPUCODE
#include <optional>
#include <utility>
#include <vector>

#include <gsl/span>

#include "ITSMFTTracking/SurfaceCatalogView.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/SurfaceId.h"
#include "ITSMFTTracking/SurfaceMask.h"

namespace o2::itsmft::tracking
{

struct SurfaceGraphSubgraph {
  std::vector<SurfaceId> orderedSurfaces;
  int maxHoles{0};
  SurfaceMask holeSurfaces{};
  SurfaceMask seedingSurfaces{};
};

enum class SurfaceGraphBuildError : uint8_t {
  None,
  EmptySubgraph,
  InvalidSubgraphSurfaceId,
  DuplicateSurfaceInSubgraph,
  SurfaceDuplicatedAcrossSubgraphs,
  NegativeMaxHoles,
  HoleSurfacesOutsideSubgraph,
  SeedingSurfacesOutsideSubgraph,
  TopologyRejected,
  GraphRejected
};

struct SurfaceGraphBuildResult {
  std::optional<SurfaceGraph> graph{};
  SurfaceGraphBuildError error{SurfaceGraphBuildError::None};
  SurfaceGraphTopologyError topologyError{SurfaceGraphTopologyError::None};
  SurfaceGraphError graphError{SurfaceGraphError::None};

  bool ok() const noexcept { return graph.has_value(); }
};

class SurfaceGraphBuilder
{
 public:
  explicit SurfaceGraphBuilder(SurfaceCatalogView catalog) : mCatalog{catalog} {}

  SurfaceGraphBuilder& addSubgraph(SurfaceGraphSubgraph subgraph)
  {
    mSubgraphs.push_back(std::move(subgraph));
    return *this;
  }

  SurfaceGraphBuildResult build() const;

 private:
  SurfaceCatalogView mCatalog;
  std::vector<SurfaceGraphSubgraph> mSubgraphs;
};

struct SurfaceGraphBatchResult {
  SurfaceGraphBuildError error{SurfaceGraphBuildError::None};
  SurfaceGraphBuildError detail{SurfaceGraphBuildError::None};
  size_t failedIteration{static_cast<size_t>(-1)};
  SurfaceGraphTopologyError topologyError{SurfaceGraphTopologyError::None};
  SurfaceGraphError graphError{SurfaceGraphError::None};
  std::vector<SurfaceGraph> graphs{};

  bool ok() const noexcept { return error == SurfaceGraphBuildError::None; }
};

SurfaceGraphBatchResult buildSurfaceGraphs(SurfaceCatalogView catalog,
                                           gsl::span<const SurfaceId> orderedSurfaces,
                                           gsl::span<const o2::itsmft::TrackingParameters> trackingParameters);

} // namespace o2::itsmft::tracking
#endif // GPUCA_GPUCODE

#endif
