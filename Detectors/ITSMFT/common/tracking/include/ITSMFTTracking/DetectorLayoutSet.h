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
#include <utility>
#include <vector>

#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"
#include "ITSMFTTracking/LayerMask.h"
#include "ITSMFTTracking/NominalSurfaceMaterial.h"

namespace o2::itsmft::tracking
{

using DetectorGeometryEpoch = uint64_t;
inline constexpr DetectorGeometryEpoch InitialDetectorGeometryEpoch = 1;

constexpr DetectorGeometryEpoch nextDetectorGeometryEpoch(DetectorGeometryEpoch epoch) noexcept
{
  return epoch == std::numeric_limits<DetectorGeometryEpoch>::max() ? InitialDetectorGeometryEpoch : epoch + 1;
}

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

struct DetectorLayoutConfigurationKey {
  DetectorGeometryEpoch geometryEpoch{InitialDetectorGeometryEpoch};
  DetectorSurfaceCatalogRequest catalogRequest{};
  std::vector<SurfaceId> orderedSurfaces{};
  std::vector<DetectorLayoutIterationConfiguration> iterations{};
  TransitionPolicyTag policyTag{TransitionPolicyTag::Invalid};

  bool operator==(const DetectorLayoutConfigurationKey& other) const noexcept
  {
    return geometryEpoch == other.geometryEpoch && catalogRequest == other.catalogRequest && orderedSurfaces == other.orderedSurfaces &&
           iterations == other.iterations && policyTag == other.policyTag;
  }
};

// Owns the complete shared detector description exactly once: the dense
// surface catalog, its parallel nominal-material budgets (same index space),
// and the full-catalogue cylinder/disk masks derived from it. Iteration-
// specific DetectorLayout objects hold no copy of any of these; they are
// validated against a borrowed view at construction and combined with the
// shared arrays only when a DetectorLayoutView is assembled (getLayoutView()
// below, the sole production assembly point).
//
// Precondition: nominalMaterial.size() == catalog.size() (parallel, same-
// index arrays). This is not re-validated here -- the caller (currently
// TimeFrame::ensureDetectorLayouts()) is responsible for supplying already
// size-matched, already-validated data; see NominalSurfaceMaterialEntry-based
// validation added alongside the caller.
class DetectorLayoutSet
{
 public:
  DetectorLayoutSet(DetectorLayoutConfigurationKey key, std::vector<SurfaceDescriptor> catalog,
                    std::vector<NominalSurfaceMaterialBudget> nominalMaterial, std::vector<DetectorLayout> layouts) noexcept
    : mConfigurationKey{std::move(key)}, mCatalog{std::move(catalog)}, mNominalMaterial{std::move(nominalMaterial)}, mLayouts{std::move(layouts)}
  {
    const auto masks = computeSurfaceKindMasks(mCatalog);
    mCylinderSurfaces = masks.first;
    mDiskSurfaces = masks.second;
  }

  DetectorLayoutSet(DetectorLayoutSet&&) noexcept = default;
  DetectorLayoutSet& operator=(DetectorLayoutSet&&) noexcept = default;

  const DetectorLayoutConfigurationKey& getConfigurationKey() const noexcept { return mConfigurationKey; }
  const std::vector<SurfaceDescriptor>& getSurfaceCatalog() const noexcept { return mCatalog; }
  const std::vector<NominalSurfaceMaterialBudget>& getNominalMaterial() const noexcept { return mNominalMaterial; }
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
    return layout ? layout->getView(mCatalog, mCylinderSurfaces, mDiskSurfaces) : DetectorLayoutView{};
  }

 private:
  DetectorLayoutConfigurationKey mConfigurationKey;
  std::vector<SurfaceDescriptor> mCatalog;
  std::vector<NominalSurfaceMaterialBudget> mNominalMaterial;
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
  MissingDetectorSurfaceIndex
};

struct DetectorLayoutSetBuildResult {
  DetectorLayoutSetBuildError error{DetectorLayoutSetBuildError::None};
  DetectorSurfaceCatalogError catalogError{DetectorSurfaceCatalogError::None};
  DetectorSurfaceCatalogValidationError catalogValidationError{DetectorSurfaceCatalogValidationError::None};
  size_t failedIteration{std::numeric_limits<size_t>::max()};
  DetectorLayoutBuildError layoutBuildError{DetectorLayoutBuildError::None};
  TopologyBuildError topologyError{TopologyBuildError::None};
  DetectorLayoutError layoutError{DetectorLayoutError::None};
  bool rebuilt{false};

  bool ok() const noexcept { return error == DetectorLayoutSetBuildError::None; }
};

} // namespace o2::itsmft::tracking
#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTORLAYOUTSET_H_ */
