// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORLAYOUTBUILDER_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORLAYOUTBUILDER_H_

#include <cstdint>

#ifndef GPUCA_GPUCODE
#include <optional>
#include <vector>

#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/SparseTrackingTopology.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceId.h"
#include "ITSMFTTracking/SurfaceMask.h"
#include "ITSMFTTracking/TransitionPolicy.h"

namespace o2::itsmft::tracking
{

/// One explicit traversal order over a subset of the catalog surfaces.
///
/// `orderedSurfaces` fixes the adjacency used to derive transitions and cell
/// holes: skipped surfaces and hole counts are computed from *positions* in
/// this span, never from a numeric comparison of the SurfaceId values it
/// contains. A catalog is free to assign global ids in an order that does not
/// follow the physical traversal (see testDetectorLayoutBuilder.cxx's
/// non-monotonic chain case).
struct DetectorLayoutSubgraph {
  std::vector<SurfaceId> orderedSurfaces;
  int maxHoles{0};
  SurfaceMask holeSurfaces{};    // subset of orderedSurfaces allowed to be skipped as a hole
  SurfaceMask seedingSurfaces{}; // subset of orderedSurfaces usable as a seed
  TransitionPolicyTag policyTag{TransitionPolicyTag::Invalid};
};

/// Failures the builder itself detects before (or instead of) delegating to
/// SparseTrackingTopology / DetectorLayout. TopologyRejected and
/// LayoutRejected mean the failure is one those types already know how to
/// diagnose; consult DetectorLayoutBuildResult::topologyError /
/// ::layoutError respectively.
enum class DetectorLayoutBuildError : uint8_t {
  None,
  EmptySubgraph,
  InvalidSubgraphSurfaceId,
  DuplicateSurfaceInSubgraph,
  SurfaceDuplicatedAcrossSubgraphs,
  NegativeMaxHoles,
  HoleSurfacesOutsideSubgraph,
  SeedingSurfacesOutsideSubgraph,
  TopologyRejected,
  LayoutRejected
};

/// Outcome of a DetectorLayoutBuilder::build() call. A failed build never
/// carries a constructed-but-invalid DetectorLayout: `layout` is populated
/// only when the build fully succeeds, so `ok()` and `layout.has_value()`
/// always agree. The three error fields are diagnostics, not a discriminated
/// union: only the field matching `error` (or, for `TopologyRejected` /
/// `LayoutRejected`, the correspondingly named field) is meaningful.
struct DetectorLayoutBuildResult {
  std::optional<DetectorLayout> layout{};
  DetectorLayoutBuildError error{DetectorLayoutBuildError::None};
  TopologyBuildError topologyError{TopologyBuildError::None};
  DetectorLayoutError layoutError{DetectorLayoutError::None};

  bool ok() const noexcept { return layout.has_value(); }
};

/// Builds a DetectorLayout from a dense surface catalog and one or more
/// explicit traversal orders over it.
///
/// The catalog is always the complete surface set: its size is the topology
/// surface count regardless of how many surfaces any individual subgraph
/// activates, and global SurfaceIds are never renumbered or shrunk to fit a
/// subgraph. Multiple subgraphs may be added to build a single disconnected
/// layout in one call (e.g. an ITS-like cylinder stack alongside an MFT-like
/// disk stack); their surfaces must be disjoint.
class DetectorLayoutBuilder
{
 public:
  explicit DetectorLayoutBuilder(std::vector<SurfaceDescriptor> catalog) : mCatalog{std::move(catalog)} {}

  DetectorLayoutBuilder& addSubgraph(DetectorLayoutSubgraph subgraph)
  {
    mSubgraphs.push_back(std::move(subgraph));
    return *this;
  }

  DetectorLayoutBuildResult build() const;

 private:
  std::vector<SurfaceDescriptor> mCatalog;
  std::vector<DetectorLayoutSubgraph> mSubgraphs;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTORLAYOUTBUILDER_H_ */
