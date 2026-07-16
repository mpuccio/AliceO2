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
#endif

#include "ITSMFTTracking/SparseTrackingTopology.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"

namespace o2::itsmft::tracking
{

struct DetectorLayoutView {
  const SurfaceDescriptor* surfaces{nullptr};
  uint32_t nSurfaces{0};
  SurfaceMask cylinderSurfaces{};
  SurfaceMask diskSurfaces{};
  SparseTrackingTopologyView topology{};

  GPUhdi() const SurfaceDescriptor& getSurface(SurfaceId id) const { return surfaces[id.value()]; }
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

class DetectorLayout
{
 public:
  DetectorLayout(std::vector<SurfaceDescriptor> surfaces, SparseTrackingTopology topology)
    : mSurfaces{std::move(surfaces)}, mTopology{std::move(topology)}
  {
    validate();
  }

  DetectorLayoutView getView() const noexcept
  {
    return valid() ? DetectorLayoutView{mSurfaces.data(), static_cast<uint32_t>(mSurfaces.size()), mCylinderSurfaces, mDiskSurfaces, mTopology.getView()}
                   : DetectorLayoutView{};
  }

  bool valid() const noexcept { return mError == DetectorLayoutError::None; }
  DetectorLayoutError getError() const noexcept { return mError; }
  const auto& getSurfaces() const noexcept { return mSurfaces; }
  const auto& getTopology() const noexcept { return mTopology; }
  SurfaceMask getCylinderSurfaces() const noexcept { return mCylinderSurfaces; }
  SurfaceMask getDiskSurfaces() const noexcept { return mDiskSurfaces; }

 private:
  void validate()
  {
    if (mSurfaces.size() > MaxLayoutSurfaces) {
      mError = DetectorLayoutError::TooManySurfaces;
      return;
    }
    if (mTopology.getSurfaceCount() != mSurfaces.size()) {
      mError = DetectorLayoutError::TopologySurfaceCountMismatch;
      return;
    }
    if (!mTopology.isFinalized()) {
      mError = DetectorLayoutError::TopologyNotFinalized;
      return;
    }
    for (uint16_t i = 0; i < mSurfaces.size(); ++i) {
      const auto& surface = mSurfaces[i];
      if (surface.id != SurfaceId{i}) {
        mError = DetectorLayoutError::NonDenseSurfaceIds;
        return;
      }
      if (surface.kind == SurfaceKind::Cylinder) {
        mCylinderSurfaces.set(surface.id);
      } else {
        mDiskSurfaces.set(surface.id);
      }
    }
    for (const auto& transition : mTopology.getTransitions()) {
      const auto fromKind = mSurfaces[transition.from.value()].kind;
      const auto toKind = mSurfaces[transition.to.value()].kind;
      if (fromKind != toKind) {
        mError = DetectorLayoutError::MixedSurfaceTransition;
        return;
      }
      const bool matchesPolicy =
        (transition.policyTag == TransitionPolicyTag::CylinderCylinder && fromKind == SurfaceKind::Cylinder) ||
        (transition.policyTag == TransitionPolicyTag::DiskDisk && fromKind == SurfaceKind::Disk);
      if (!matchesPolicy) {
        mError = DetectorLayoutError::PolicySurfaceKindMismatch;
        return;
      }
    }
  }

  std::vector<SurfaceDescriptor> mSurfaces;
  SparseTrackingTopology mTopology;
  SurfaceMask mCylinderSurfaces{};
  SurfaceMask mDiskSurfaces{};
  DetectorLayoutError mError{DetectorLayoutError::None};
};

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTORLAYOUT_H_ */
