// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYSTATE_H_
#define ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYSTATE_H_

#include <cstddef>
#include <type_traits>

#include "GPUCommonDef.h"
#include "GPUCommonMath.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/TransitionPolicy.h"

namespace o2::itsmft::tracking
{

/// True for finite `x` (excludes +-Inf and NaN, both of which set every
/// exponent bit of an IEEE-754 single-precision value). GPUCommonMath::Finite()/
/// IsNaN() are a float-determinism switch (TPC compression) that unconditionally
/// return "finite"/"not NaN" outside GPUCA_DETERMINISTIC_MODE, so they cannot
/// back a real bounds check here; this reads the bit pattern directly through
/// the same device-portable reinterpretation GPUCommonMath itself uses.
GPUhdi() bool isFiniteParam(float x) noexcept
{
  return (o2::gpu::GPUCommonMath::Float2UIntReint(x) & 0x7f800000u) != 0x7f800000u;
}

/// Device-compatible cylinder-cylinder policy parameter boundary. Every cut
/// is a named, bounds-checked field; an untyped flat-float block is not an
/// accepted substitute for this ABI. Defaults mirror the current production
/// barrel defaults in TrackingParameters (Configuration.h); this struct is a
/// separate contract and does not yet feed production tracking.
struct CylinderCylinderPolicyParams {
  float trackletMinPt{0.3f};
  float cellDeltaTanLambdaSigma{0.007f};
  float nSigmaCut{5.f};
  float maxChi2ClusterAttachment{60.f};
  float maxChi2NDF{30.f};
  float pvResolution{1.e-2f};

  GPUhdi() bool isValid() const noexcept
  {
    return isFiniteParam(trackletMinPt) && trackletMinPt > 0.f &&
           isFiniteParam(cellDeltaTanLambdaSigma) && cellDeltaTanLambdaSigma > 0.f &&
           isFiniteParam(nSigmaCut) && nSigmaCut > 0.f &&
           isFiniteParam(maxChi2ClusterAttachment) && maxChi2ClusterAttachment > 0.f &&
           isFiniteParam(maxChi2NDF) && maxChi2NDF > 0.f &&
           isFiniteParam(pvResolution) && pvResolution >= 0.f;
  }
};

/// Device-compatible disk-disk policy parameter boundary. `trackletMinPt` and
/// `cellDeltaTanLambdaSigma` are required here, not only in the barrel
/// struct, because the disk path reads TrackletMinPt during MFT projection
/// and applies CellDeltaTanLambdaSigma ahead of the detector-specific
/// cell-building branch (TrackerTraits::computeLayerTracklets/computeLayerCells).
/// Defaults mirror resetDetectorDefaults(..., MFT) (Configuration.cxx):
/// TrackletMinPt and CellDeltaTanLambdaSigma are left at the TrackingParameters
/// struct defaults (MFT does not override them), while TrackletMinAbsX is
/// explicitly set to 0.05f.
struct DiskDiskPolicyParams {
  float trackletMinPt{0.3f};
  float cellDeltaTanLambdaSigma{0.007f};
  float cellRoadRCut{0.05f};
  float trackletMinAbsX{0.05f};
  float nSigmaCut{5.f};
  float maxChi2ClusterAttachment{60.f};
  float maxChi2NDF{30.f};

  GPUhdi() bool isValid() const noexcept
  {
    return isFiniteParam(trackletMinPt) && trackletMinPt > 0.f &&
           isFiniteParam(cellDeltaTanLambdaSigma) && cellDeltaTanLambdaSigma > 0.f &&
           isFiniteParam(cellRoadRCut) && cellRoadRCut > 0.f &&
           isFiniteParam(trackletMinAbsX) && trackletMinAbsX >= 0.f &&
           isFiniteParam(nSigmaCut) && nSigmaCut > 0.f &&
           isFiniteParam(maxChi2ClusterAttachment) && maxChi2ClusterAttachment > 0.f &&
           isFiniteParam(maxChi2NDF) && maxChi2NDF > 0.f;
  }
};

// Device-facing ABI lock: standard-layout/trivially-copyable, exact size and
// alignment, and a fixed byte offset for every field. Any reordering,
// insertion, or width change is a breaking ABI change and must fail here
// first.
static_assert(std::is_standard_layout_v<CylinderCylinderPolicyParams> && std::is_trivially_copyable_v<CylinderCylinderPolicyParams>);
static_assert(sizeof(CylinderCylinderPolicyParams) == 24);
static_assert(alignof(CylinderCylinderPolicyParams) == alignof(float));
static_assert(offsetof(CylinderCylinderPolicyParams, trackletMinPt) == 0);
static_assert(offsetof(CylinderCylinderPolicyParams, cellDeltaTanLambdaSigma) == 4);
static_assert(offsetof(CylinderCylinderPolicyParams, nSigmaCut) == 8);
static_assert(offsetof(CylinderCylinderPolicyParams, maxChi2ClusterAttachment) == 12);
static_assert(offsetof(CylinderCylinderPolicyParams, maxChi2NDF) == 16);
static_assert(offsetof(CylinderCylinderPolicyParams, pvResolution) == 20);

static_assert(std::is_standard_layout_v<DiskDiskPolicyParams> && std::is_trivially_copyable_v<DiskDiskPolicyParams>);
static_assert(sizeof(DiskDiskPolicyParams) == 28);
static_assert(alignof(DiskDiskPolicyParams) == alignof(float));
static_assert(offsetof(DiskDiskPolicyParams, trackletMinPt) == 0);
static_assert(offsetof(DiskDiskPolicyParams, cellDeltaTanLambdaSigma) == 4);
static_assert(offsetof(DiskDiskPolicyParams, cellRoadRCut) == 8);
static_assert(offsetof(DiskDiskPolicyParams, trackletMinAbsX) == 12);
static_assert(offsetof(DiskDiskPolicyParams, nSigmaCut) == 16);
static_assert(offsetof(DiskDiskPolicyParams, maxChi2ClusterAttachment) == 20);
static_assert(offsetof(DiskDiskPolicyParams, maxChi2NDF) == 24);

/// Stage-B activation compatibility hypothesis, used only at initial
/// Cell-state construction (native buildCellSeed<Tag>, TransitionPolicyOperations.h);
/// never re-set afterward. Centrally named so every call site shares one
/// definition rather than repeating a magic charge/PID literal.
inline constexpr uint8_t kCompatibilityAbsCharge = 1;
inline const o2::track::PID kCompatibilityPID = o2::track::PID::Pion;

/// Per-tag policy/state boundary: the derived state family, the surface kind
/// every transition carrying this tag must have, the Stage-A track state used
/// by cell/road construction, and the tag-specific bounds-checked parameter
/// block. One instantiation per accepted TransitionPolicyTag; instantiating
/// the primary template for an unsupported tag is a compile error rather than
/// a silent fallback.
template <TransitionPolicyTag Tag>
struct TransitionPolicyTraits;

template <>
struct TransitionPolicyTraits<TransitionPolicyTag::CylinderCylinder> {
  static constexpr TransitionPolicyTag Tag = TransitionPolicyTag::CylinderCylinder;
  static constexpr StateFamily Family = StateFamily::Barrel;
  static constexpr SurfaceKind ExpectedSurfaceKind = SurfaceKind::Cylinder;
  using SeedState = SurfaceKinematicState;
  using Params = CylinderCylinderPolicyParams;
};

template <>
struct TransitionPolicyTraits<TransitionPolicyTag::DiskDisk> {
  static constexpr TransitionPolicyTag Tag = TransitionPolicyTag::DiskDisk;
  static constexpr StateFamily Family = StateFamily::Forward;
  static constexpr SurfaceKind ExpectedSurfaceKind = SurfaceKind::Disk;
  using SeedState = SurfaceKinematicState;
  using Params = DiskDiskPolicyParams;
};

static_assert(TransitionPolicyTraits<TransitionPolicyTag::CylinderCylinder>::Family ==
              stateFamilyOf(TransitionPolicyTag::CylinderCylinder));
static_assert(TransitionPolicyTraits<TransitionPolicyTag::DiskDisk>::Family ==
              stateFamilyOf(TransitionPolicyTag::DiskDisk));
static_assert(isSurfaceKindCompatible(TransitionPolicyTag::CylinderCylinder,
                                      TransitionPolicyTraits<TransitionPolicyTag::CylinderCylinder>::ExpectedSurfaceKind));
static_assert(isSurfaceKindCompatible(TransitionPolicyTag::DiskDisk,
                                      TransitionPolicyTraits<TransitionPolicyTag::DiskDisk>::ExpectedSurfaceKind));
static_assert(!isSurfaceKindCompatible(TransitionPolicyTag::CylinderCylinder, SurfaceKind::Disk));
static_assert(!isSurfaceKindCompatible(TransitionPolicyTag::DiskDisk, SurfaceKind::Cylinder));

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYSTATE_H_ */
