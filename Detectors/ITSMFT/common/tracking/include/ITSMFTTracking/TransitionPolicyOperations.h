// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYOPERATIONS_H_
#define ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYOPERATIONS_H_

#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/TransitionPolicyState.h"

#ifndef GPUCA_GPUCODE
#include "DetectorsBase/Propagator.h"

namespace o2::its
{
struct TrackingFrameInfo;
}
#endif

namespace o2::itsmft::tracking
{

/// D007 policy-boundary operation (Architecture.md 10): true if `nextCell`
/// may extend a road built on `currentCell`. cellsAreCompatible and attachHit
/// are now migrated to this boundary; projectSearchWindow, buildCellSeed, and
/// finalRefit remain TrackerTraits-internal until later slices. Dispatch is a
/// compile-time tag switch -- one specialization per Tag, selected by the
/// caller once outside any candidate/neighbour loop -- never a runtime detector
/// branch inside this operation or its callers' hot loops. Instantiating the
/// primary template for an unsupported tag is a compile error rather than a
/// silent fallback (mirrors TransitionPolicyTraits).
///
/// Host-only for this CPU migration slice: the disk (DiskDisk) specialization
/// is defined in TransitionPolicyOperations.cxx against the existing MFT
/// forward-state helpers (TrackParCovFwd::propagateToZhelix/Linear and
/// friends), which are themselves host-only. This operation must not be
/// declared or assumed device-compatible until a device-capable forward
/// state/propagator exists; do not add GPU qualifiers here without first
/// making that dependency chain device-compatible. The public header
/// deliberately does not include MFTFwdTrackHelpers.h, so the MFT-specific
/// helper stays behind this implementation boundary rather than leaking MFT
/// constants/TimeFrame/MFTCATrack dependencies into a common policy header.
template <TransitionPolicyTag Tag>
bool cellsAreCompatible(const CellSeedTpl<typename TransitionPolicyTraits<Tag>::SeedState>& currentCell,
                        const CellSeedTpl<typename TransitionPolicyTraits<Tag>::SeedState>& nextCell,
                        float bz,
                        const typename TransitionPolicyTraits<Tag>::Params& params);

/// Barrel formula: rotate/propagate `nextCell` into `currentCell`'s frame and
/// accept if the predicted chi2 stays within the bound. `nextCell` is never
/// mutated -- the propagated state is a local copy -- so callers may pass a
/// stored cell directly without pre-copying it for this call. Defined out of
/// line in TransitionPolicyOperations.cxx.
template <>
bool cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(
  const CellSeedTpl<o2::track::TrackParCovF>& currentCell,
  const CellSeedTpl<o2::track::TrackParCovF>& nextCell,
  float bz,
  const CylinderCylinderPolicyParams& params);

/// Disk formula: delegates to the existing MFT forward-state helper, which
/// already performs its own internal copy/propagation; neither input cell is
/// mutated. Defined out of line in TransitionPolicyOperations.cxx, the only
/// translation unit permitted to include MFTFwdTrackHelpers.h on behalf of
/// this policy operation.
template <>
bool cellsAreCompatible<TransitionPolicyTag::DiskDisk>(
  const CellSeedTpl<o2::track::TrackParCovFwd>& currentCell,
  const CellSeedTpl<o2::track::TrackParCovFwd>& nextCell,
  float bz,
  const DiskDiskPolicyParams& params);

#ifndef GPUCA_GPUCODE

/// D007 attach-hit policy operation. The typed family state and parameter
/// block are selected once by the caller's outer policy dispatch. `xOverX0`
/// is the already-selected material budget for the hit surface; no legacy
/// TrackingParameters object or detector identity crosses this boundary.
/// `o2::its::TrackingFrameInfo` is a temporary Gate 3 compatibility boundary,
/// not a detector-neutral measurement contract; production normalized loading
/// must eventually supply SurfaceMeasurement directly.
///
/// Host-only: CylinderCylinder calls the host Propagator singleton and
/// DiskDisk calls the host forward-state propagation/update chain. The two
/// specializations deliberately preserve their existing failure mutation
/// contracts: the cylinder state is modified as each successful step is
/// applied, while the disk state and chi2 are committed only after the full
/// attachment succeeds.
template <TransitionPolicyTag Tag>
bool attachHit(typename TransitionPolicyTraits<Tag>::SeedState& state,
               const o2::its::TrackingFrameInfo& hit,
               float xOverX0,
               o2::base::PropagatorF::MatCorrType corrType,
               float bz,
               float& chi2,
               const typename TransitionPolicyTraits<Tag>::Params& params);

template <>
bool attachHit<TransitionPolicyTag::CylinderCylinder>(
  o2::track::TrackParCovF& state,
  const o2::its::TrackingFrameInfo& hit,
  float xOverX0,
  o2::base::PropagatorF::MatCorrType corrType,
  float bz,
  float& chi2,
  const CylinderCylinderPolicyParams& params);

template <>
bool attachHit<TransitionPolicyTag::DiskDisk>(
  o2::track::TrackParCovFwd& state,
  const o2::its::TrackingFrameInfo& hit,
  float xOverX0,
  o2::base::PropagatorF::MatCorrType corrType,
  float bz,
  float& chi2,
  const DiskDiskPolicyParams& params);

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYOPERATIONS_H_ */
