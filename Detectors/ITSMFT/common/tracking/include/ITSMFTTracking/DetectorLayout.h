// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORLAYOUT_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORLAYOUT_H_

#include <cstdint>
#include <type_traits>

#ifndef GPUCA_GPUCODE
#include <utility>
#include <vector>

#include <gsl/gsl>
#endif

#include "ITSMFTTracking/SparseTrackingTopology.h"
#include "ITSMFTTracking/SurfaceCatalogView.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/TransitionPolicy.h"

namespace o2::itsmft::tracking
{

struct DetectorLayoutView {
  const SurfaceDescriptor* surfaces{nullptr};
  uint32_t nSurfaces{0};
  SurfaceMask cylinderSurfaces{};
  SurfaceMask diskSurfaces{};
  SparseTrackingTopologyView topology{};

  GPUhdi() const SurfaceDescriptor& getSurface(SurfaceId id) const { return surfaces[id.value()]; }
  // Narrowed, catalog-only view: geometry identity alone, with no topology,
  // masks or transition-policy dependency.
  GPUhdi() SurfaceCatalogView getSurfaceCatalogView() const noexcept { return SurfaceCatalogView{surfaces, nSurfaces}; }
};

static_assert(std::is_standard_layout_v<DetectorLayoutView>);
static_assert(std::is_trivially_copyable_v<DetectorLayoutView>);

#ifndef GPUCA_GPUCODE

enum class DetectorLayoutError : uint8_t {
  None,
  TooManySurfaces,
  NonDenseSurfaceIds,
  TopologySurfaceCountMismatch,
  TopologyNotFinalized,
  MixedSurfaceTransition,
  PolicySurfaceKindMismatch
};

// Full-catalogue cylinder/disk masks: a pure function of a dense surface
// catalog's `.kind` values, independent of any iteration/subgraph. Callers
// that own a shared catalog (DetectorLayoutSet) should compute this once and
// cache the result rather than calling it per iteration/per view.
inline std::pair<SurfaceMask, SurfaceMask> computeSurfaceKindMasks(gsl::span<const SurfaceDescriptor> surfaces) noexcept
{
  SurfaceMask cylinderSurfaces{};
  SurfaceMask diskSurfaces{};
  for (const auto& surface : surfaces) {
    if (surface.kind == SurfaceKind::Cylinder) {
      cylinderSurfaces.set(surface.id);
    } else {
      diskSurfaces.set(surface.id);
    }
  }
  return {cylinderSurfaces, diskSurfaces};
}

// Owns only iteration-specific topology and validation state. Surface
// geometry, nominal material and the full-catalogue cylinder/disk masks are
// not owned here -- they live once in the shared detector description
// (DetectorLayoutSet) that every iteration's DetectorLayout is validated
// against. `surfaces` is borrowed for construction/validation only and is
// never retained past the constructor call.
class DetectorLayout
{
 public:
  DetectorLayout(gsl::span<const SurfaceDescriptor> surfaces, SparseTrackingTopology topology)
    : mTopology{std::move(topology)}
  {
    validate(surfaces);
  }

  // Assembles a DetectorLayoutView from this iteration's topology plus a
  // caller-supplied shared surface catalog and its precomputed full-catalogue
  // masks. Precondition: `surfaces`/`cylinderSurfaces`/`diskSurfaces` are the
  // exact shared description this DetectorLayout was validated against
  // (DetectorLayoutSet::getLayoutView() is the sole production caller;
  // behavior is unspecified for a mismatched catalog).
  DetectorLayoutView getView(gsl::span<const SurfaceDescriptor> surfaces, SurfaceMask cylinderSurfaces, SurfaceMask diskSurfaces) const noexcept
  {
    return valid() ? DetectorLayoutView{surfaces.data(), static_cast<uint32_t>(surfaces.size()), cylinderSurfaces, diskSurfaces, mTopology.getView()}
                   : DetectorLayoutView{};
  }

  bool valid() const noexcept { return mError == DetectorLayoutError::None; }
  DetectorLayoutError getError() const noexcept { return mError; }
  const auto& getTopology() const noexcept { return mTopology; }

 private:
  void validate(gsl::span<const SurfaceDescriptor> surfaces)
  {
    if (surfaces.size() > MaxLayoutSurfaces) {
      mError = DetectorLayoutError::TooManySurfaces;
      return;
    }
    if (mTopology.getSurfaceCount() != surfaces.size()) {
      mError = DetectorLayoutError::TopologySurfaceCountMismatch;
      return;
    }
    if (!mTopology.isFinalized()) {
      mError = DetectorLayoutError::TopologyNotFinalized;
      return;
    }
    for (uint16_t i = 0; i < surfaces.size(); ++i) {
      if (surfaces[i].id != SurfaceId{i}) {
        mError = DetectorLayoutError::NonDenseSurfaceIds;
        return;
      }
    }
    for (const auto& transition : mTopology.getTransitions()) {
      const auto fromKind = surfaces[transition.from.value()].kind;
      const auto toKind = surfaces[transition.to.value()].kind;
      if (fromKind != toKind) {
        mError = DetectorLayoutError::MixedSurfaceTransition;
        return;
      }
      if (!isSurfaceKindCompatible(transition.policyTag, fromKind)) {
        mError = DetectorLayoutError::PolicySurfaceKindMismatch;
        return;
      }
    }
  }

  SparseTrackingTopology mTopology;
  DetectorLayoutError mError{DetectorLayoutError::None};
};

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTORLAYOUT_H_ */
