// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORLAYOUTSET_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORLAYOUTSET_H_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#ifndef GPUCA_GPUCODE
#include <gsl/span>
#include <optional>
#include <utility>
#include <vector>

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"
#include "ITSMFTTracking/LayerMask.h"
#include "ITSMFTTracking/SurfaceCatalogView.h"

namespace o2::itsmft::tracking
{

struct DetectorLayoutIterationConfiguration {
  uint32_t activeCount{0};
  int maxHoles{0};
  LayerMask holeLayerMask{};
  LayerMask startLayerMask{};

  bool operator==(const DetectorLayoutIterationConfiguration& other) const noexcept
  {
    return activeCount == other.activeCount && maxHoles == other.maxHoles &&
           holeLayerMask.value() == other.holeLayerMask.value() &&
           startLayerMask.value() == other.startLayerMask.value();
  }
};

// No geometryEpoch/catalogRequest: this key no longer identifies a runtime
// provider query to re-issue or compare against (Gate 4 B2 Slice 2 -- the
// catalog is now a static, process-lifetime SurfaceCatalogView with nothing
// to invalidate; see StaticDetectorCatalogs.h). A DetectorLayoutSet is built
// exactly once by its owner (buildDetectorLayoutSet() below) and never
// rebuilt/compared for currency.
struct DetectorLayoutConfigurationKey {
  std::vector<SurfaceId> orderedSurfaces{};
  std::vector<DetectorLayoutIterationConfiguration> iterations{};
  TransitionPolicyTag policyTag{TransitionPolicyTag::Invalid};

  bool operator==(const DetectorLayoutConfigurationKey& other) const noexcept
  {
    return orderedSurfaces == other.orderedSurfaces && iterations == other.iterations && policyTag == other.policyTag;
  }
};

// Owns only the per-iteration layouts/topology and this configuration key;
// borrows the shared surface catalog (SurfaceCatalogView) rather than owning
// a copy of it. The full-catalogue cylinder/disk masks are still derived and
// cached once here. Iteration-specific DetectorLayout objects hold no copy
// of the catalog either; they are validated against a borrowed span at
// construction and combined with it only when a DetectorLayoutView is
// assembled (getLayoutView() below, the sole production assembly point).
class DetectorLayoutSet
{
 public:
  DetectorLayoutSet(DetectorLayoutConfigurationKey key, SurfaceCatalogView catalog,
                    std::vector<DetectorLayout> layouts) noexcept
    : mConfigurationKey{std::move(key)}, mCatalog{catalog}, mLayouts{std::move(layouts)}
  {
    const auto masks = computeSurfaceKindMasks(gsl::span<const SurfaceDescriptor>{mCatalog.surfaces, mCatalog.nSurfaces});
    mCylinderSurfaces = masks.first;
    mDiskSurfaces = masks.second;
  }

  DetectorLayoutSet(DetectorLayoutSet&&) noexcept = default;
  DetectorLayoutSet& operator=(DetectorLayoutSet&&) noexcept = default;

  const DetectorLayoutConfigurationKey& getConfigurationKey() const noexcept { return mConfigurationKey; }
  SurfaceCatalogView getSurfaceCatalog() const noexcept { return mCatalog; }
  SurfaceMask getCylinderSurfaces() const noexcept { return mCylinderSurfaces; }
  SurfaceMask getDiskSurfaces() const noexcept { return mDiskSurfaces; }
  size_t size() const noexcept { return mLayouts.size(); }
  const std::vector<DetectorLayout>& getLayouts() const noexcept { return mLayouts; }
  const DetectorLayout* getLayout(size_t iteration) const noexcept
  {
    return iteration < mLayouts.size() ? &mLayouts[iteration] : nullptr;
  }
  DetectorLayoutView getLayoutView(size_t iteration) const noexcept
  {
    const auto* layout = getLayout(iteration);
    return layout ? layout->getView(gsl::span<const SurfaceDescriptor>{mCatalog.surfaces, mCatalog.nSurfaces}, mCylinderSurfaces, mDiskSurfaces) : DetectorLayoutView{};
  }

 private:
  DetectorLayoutConfigurationKey mConfigurationKey;
  SurfaceCatalogView mCatalog;
  SurfaceMask mCylinderSurfaces{};
  SurfaceMask mDiskSurfaces{};
  std::vector<DetectorLayout> mLayouts;
};

static_assert(std::is_nothrow_move_constructible_v<DetectorLayoutSet>);

enum class DetectorLayoutSetBuildError : uint8_t {
  None,
  MissingProvider,
  CatalogProviderFailure,
  InvalidCatalog,
  InvalidActiveCount,
  LayoutBuilderFailure
};

enum class DetectorSurfaceCatalogValidationError : uint8_t {
  None,
  InvalidDetector,
  InvalidFirstSurface,
  EmptyDetector,
  TooManySurfaces,
  SizeMismatch,
  NonDenseGlobalSurfaceIds,
  DetectorMismatch,
  DetectorSurfaceIndexOutOfRange,
  DuplicateDetectorSurfaceIndex,
  MissingDetectorSurfaceIndex,
  InvalidMaterial
};

// catalogError/catalogValidationError are diagnostics inherited from the
// pre-Slice-2 runtime-provider build path; buildDetectorLayoutSet() below
// never produces MissingProvider/CatalogProviderFailure/InvalidCatalog or a
// non-None catalogValidationError (there is no provider to fail and no
// runtime catalog to validate -- the static catalog it borrows is already
// proven valid at compile time, see SurfaceSpec.h), so both stay at their
// default (::None) on every real call. `rebuilt` is always true on success:
// there is no currency concept left to make a rebuild conditional on (see
// DetectorLayoutConfigurationKey's own doc).
struct DetectorLayoutSetBuildResult {
  DetectorLayoutSetBuildError error{DetectorLayoutSetBuildError::None};
  DetectorSurfaceCatalogError catalogError{DetectorSurfaceCatalogError::None};
  DetectorSurfaceCatalogValidationError catalogValidationError{DetectorSurfaceCatalogValidationError::None};
  size_t failedIteration{std::numeric_limits<size_t>::max()};
  DetectorLayoutBuildError layoutBuildError{DetectorLayoutBuildError::None};
  TopologyBuildError topologyError{TopologyBuildError::None};
  DetectorLayoutError layoutError{DetectorLayoutError::None};
  bool rebuilt{false};
  // Populated iff ok(). The built product itself -- absent from earlier
  // (TimeFrame-owned) build results, which stored it as TimeFrame state
  // instead; this result now carries it directly since the owner
  // (ITSMFTTrackingInterface) has nowhere else to receive it from.
  std::optional<DetectorLayoutSet> layout{};

  bool ok() const noexcept { return error == DetectorLayoutSetBuildError::None; }
};

// Builds one DetectorLayoutSet, once, from a borrowed SurfaceCatalogView
// (process-lifetime static storage, see StaticDetectorCatalogs.h) and this
// detector's per-iteration TrackingParameters. One DetectorLayoutSubgraph
// per iteration, built through the same (unmodified) DetectorLayoutBuilder
// every caller already used; `orderedSurfaces` fixes the legacy-layer
// traversal order exactly as it did for the pre-Slice-2 runtime path. Pure
// function of its arguments: no TimeFrame, no provider, nothing to
// invalidate or re-issue -- the intended single caller (ITSMFTTrackingInterface
// ::initialiseTracker()) calls this exactly once per interface instance.
DetectorLayoutSetBuildResult buildDetectorLayoutSet(SurfaceCatalogView catalog,
                                                    gsl::span<const SurfaceId> orderedSurfaces,
                                                    TransitionPolicyTag policyTag,
                                                    gsl::span<const o2::itsmft::TrackingParameters> trackingParameters);

} // namespace o2::itsmft::tracking
#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTORLAYOUTSET_H_ */
