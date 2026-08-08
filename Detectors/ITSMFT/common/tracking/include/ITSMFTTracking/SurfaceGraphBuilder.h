// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEGRAPHBUILDER_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEGRAPHBUILDER_H_

#include <cstddef>
#include <cstdint>
#include <type_traits>

#ifndef GPUCA_GPUCODE
#include <optional>
#include <utility>
#include <vector>

#include <gsl/span>

#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/SurfaceId.h"
#include "ITSMFTTracking/SurfaceMask.h"

namespace o2::itsmft::tracking
{

struct SurfaceAdjacencyPair {
  uint16_t fromIndex{0};
  uint16_t toIndex{0};
};

static_assert(std::is_standard_layout_v<SurfaceAdjacencyPair> && std::is_trivially_copyable_v<SurfaceAdjacencyPair>);
static_assert(sizeof(SurfaceAdjacencyPair) == 4);

struct SurfaceGraphDefinition {
  std::vector<SurfaceId> orderedSurfaces;
  std::vector<SurfaceAdjacencyPair> basePairs;
  int maxHoles{0};
  SurfaceMask holeSurfaces{};
  SurfaceMask seedingSurfaces{};
};

SurfaceGraphDefinition makeSurfaceChain(gsl::span<const SurfaceId> orderedSurfaces,
                                        int maxHoles = 0,
                                        SurfaceMask holeSurfaces = {},
                                        SurfaceMask seedingSurfaces = {});

enum class SurfaceGraphBuildError : uint8_t {
  None,
  InvalidSurfaceId,
  DuplicateSurface,
  InvalidBasePair,
  DuplicateBasePair,
  NegativeMaxHoles,
  HoleSurfacesOutsideGraph,
  SeedingSurfacesOutsideGraph,
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
  SurfaceGraphBuilder(SurfaceCatalogView catalog, SurfaceGraphDefinition definition)
    : mCatalog{catalog}, mDefinition{std::move(definition)}
  {
  }

  SurfaceGraphBuildResult build() const;

 private:
  SurfaceCatalogView mCatalog;
  SurfaceGraphDefinition mDefinition;
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
