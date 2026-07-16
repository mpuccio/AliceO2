// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYSTATE_H_
#define ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYSTATE_H_

#include <type_traits>

#include "GPUCommonDef.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/TransitionPolicy.h"

namespace o2::itsmft::tracking
{

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

  GPUhdi() constexpr bool isValid() const noexcept
  {
    return trackletMinPt > 0.f && cellDeltaTanLambdaSigma > 0.f && nSigmaCut > 0.f &&
           maxChi2ClusterAttachment > 0.f && maxChi2NDF > 0.f;
  }
};

/// Device-compatible disk-disk policy parameter boundary. Defaults mirror the
/// current production MFT-style defaults in TrackingParameters.
struct DiskDiskPolicyParams {
  float cellRoadRCut{0.05f};
  float trackletMinAbsX{0.f};
  float nSigmaCut{5.f};
  float maxChi2ClusterAttachment{60.f};
  float maxChi2NDF{30.f};

  GPUhdi() constexpr bool isValid() const noexcept
  {
    return cellRoadRCut > 0.f && trackletMinAbsX >= 0.f && nSigmaCut > 0.f &&
           maxChi2ClusterAttachment > 0.f && maxChi2NDF > 0.f;
  }
};

static_assert(std::is_standard_layout_v<CylinderCylinderPolicyParams> && std::is_trivially_copyable_v<CylinderCylinderPolicyParams>);
static_assert(std::is_standard_layout_v<DiskDiskPolicyParams> && std::is_trivially_copyable_v<DiskDiskPolicyParams>);

/// True if a transition carrying `tag` may legally connect surfaces of `kind`.
/// This mirrors, but does not replace, DetectorLayout's own construction-time
/// policy/surface-kind validation; it exists so the dispatch boundary and its
/// tests can state the same compatibility rule without depending on layout
/// construction succeeding first.
GPUhdi() constexpr bool isSurfaceKindCompatible(TransitionPolicyTag tag, SurfaceKind kind) noexcept
{
  switch (tag) {
    case TransitionPolicyTag::CylinderCylinder:
      return kind == SurfaceKind::Cylinder;
    case TransitionPolicyTag::DiskDisk:
      return kind == SurfaceKind::Disk;
    case TransitionPolicyTag::Invalid:
      return false;
  }
  return false;
}

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
  using SeedState = o2::track::TrackParCovF;
  using Params = CylinderCylinderPolicyParams;
};

template <>
struct TransitionPolicyTraits<TransitionPolicyTag::DiskDisk> {
  static constexpr TransitionPolicyTag Tag = TransitionPolicyTag::DiskDisk;
  static constexpr StateFamily Family = StateFamily::Forward;
  static constexpr SurfaceKind ExpectedSurfaceKind = SurfaceKind::Disk;
  using SeedState = o2::track::TrackParCovFwd;
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
