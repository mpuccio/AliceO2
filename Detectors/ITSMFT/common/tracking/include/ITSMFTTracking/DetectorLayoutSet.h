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

namespace o2::itsmft::tracking
{

using DetectorGeometryEpoch = uint64_t;
inline constexpr DetectorGeometryEpoch InitialDetectorGeometryEpoch = 1;

constexpr DetectorGeometryEpoch nextDetectorGeometryEpoch(DetectorGeometryEpoch epoch) noexcept
{
  return epoch == std::numeric_limits<DetectorGeometryEpoch>::max() ? InitialDetectorGeometryEpoch : epoch + 1;
}

// Independent versioning for the shared nominal-material content, distinct
// from DetectorGeometryEpoch: a material-only or geometry-only rebuild
// changes exactly one of the two. Owned/incremented by whichever caller
// detects the *material content* changed (TimeFrame::invalidateNominalMaterial());
// DetectorLayoutSet/ensureDetectorLayouts() only compare and store it.
// Changing material content without bumping this epoch is a caller-contract
// violation the framework cannot independently detect -- symmetric to the
// existing DetectorGeometryEpoch contract.
using MaterialCatalogEpoch = uint64_t;
inline constexpr MaterialCatalogEpoch InitialMaterialCatalogEpoch = 1;

constexpr MaterialCatalogEpoch nextMaterialCatalogEpoch(MaterialCatalogEpoch epoch) noexcept
{
  return epoch == std::numeric_limits<MaterialCatalogEpoch>::max() ? InitialMaterialCatalogEpoch : epoch + 1;
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
  MaterialCatalogEpoch materialEpoch{InitialMaterialCatalogEpoch};
  DetectorSurfaceCatalogRequest catalogRequest{};
  std::vector<SurfaceId> orderedSurfaces{};
  std::vector<DetectorLayoutIterationConfiguration> iterations{};
  TransitionPolicyTag policyTag{TransitionPolicyTag::Invalid};

  bool operator==(const DetectorLayoutConfigurationKey& other) const noexcept
  {
    return geometryEpoch == other.geometryEpoch && materialEpoch == other.materialEpoch && catalogRequest == other.catalogRequest &&
           orderedSurfaces == other.orderedSurfaces && iterations == other.iterations && policyTag == other.policyTag;
  }
};

// Owns the complete shared detector description exactly once: the dense
// surface catalog (each entry carrying its own nominal material, see
// SurfaceDescriptor::material) and the full-catalogue cylinder/disk masks
// derived from it. Iteration-specific DetectorLayout objects hold no copy of
// either; they are validated against a borrowed view at construction and
// combined with the shared catalog only when a DetectorLayoutView is
// assembled (getLayoutView() below, the sole production assembly point).
class DetectorLayoutSet
{
 public:
  DetectorLayoutSet(DetectorLayoutConfigurationKey key, std::vector<SurfaceDescriptor> catalog,
                    std::vector<DetectorLayout> layouts) noexcept
    : mConfigurationKey{std::move(key)}, mCatalog{std::move(catalog)}, mLayouts{std::move(layouts)}
  {
    const auto masks = computeSurfaceKindMasks(mCatalog);
    mCylinderSurfaces = masks.first;
    mDiskSurfaces = masks.second;
  }

  DetectorLayoutSet(DetectorLayoutSet&&) noexcept = default;
  DetectorLayoutSet& operator=(DetectorLayoutSet&&) noexcept = default;

  const DetectorLayoutConfigurationKey& getConfigurationKey() const noexcept { return mConfigurationKey; }
  const std::vector<SurfaceDescriptor>& getSurfaceCatalog() const noexcept { return mCatalog; }
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
  InvalidMaterial,
  InvalidActiveCount,
  LayoutBuilderFailure
};

// Diagnostic detail for DetectorLayoutSetBuildError::InvalidMaterial.
// Checked in this order (see TimeFrame::ensureDetectorLayouts()): a nonzero
// material epoch, then positional identity-bearing entry validation against
// the already-dense, already-validated surface catalog, then per-entry
// budget finiteness/non-negativity.
enum class DetectorLayoutMaterialValidationError : uint8_t {
  None,
  InvalidMaterialEpoch,
  SizeMismatch,
  SurfaceIdMismatch,
  InvalidBudget
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
  DetectorLayoutMaterialValidationError materialValidationError{DetectorLayoutMaterialValidationError::None};
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
