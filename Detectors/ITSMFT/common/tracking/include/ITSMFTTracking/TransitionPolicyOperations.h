// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYOPERATIONS_H_
#define ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYOPERATIONS_H_

#include "GPUCommonDef.h"
#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/TransitionPolicyState.h"

namespace o2::itsmft::tracking
{

/// D007 policy-boundary operation (Architecture.md 10): true if `nextCell`
/// may extend a road built on `currentCell`. This is the first of the five
/// policy operations listed there (projectSearchWindow, buildCellSeed,
/// cellsAreCompatible, attachHit, finalRefit); the other four remain
/// TrackerTraits-internal until later slices. Dispatch is a compile-time tag
/// switch -- one specialization per Tag, selected by the caller once outside
/// any candidate/neighbour loop -- never a runtime detector branch inside
/// this operation or its callers' hot loops. Instantiating the primary
/// template for an unsupported tag is a compile error rather than a silent
/// fallback (mirrors TransitionPolicyTraits).
template <TransitionPolicyTag Tag>
GPUhdi() bool cellsAreCompatible(const CellSeedTpl<typename TransitionPolicyTraits<Tag>::SeedState>& currentCell,
                                 const CellSeedTpl<typename TransitionPolicyTraits<Tag>::SeedState>& nextCell,
                                 float bz,
                                 const typename TransitionPolicyTraits<Tag>::Params& params) noexcept;

/// Barrel formula: rotate/propagate `nextCell` into `currentCell`'s frame and
/// accept if the predicted chi2 stays within the bound. `nextCell` is never
/// mutated -- the propagated state is a local copy -- so callers may pass a
/// stored cell directly without pre-copying it for this call.
template <>
GPUhdi() bool cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(
  const CellSeedTpl<o2::track::TrackParCovF>& currentCell,
  const CellSeedTpl<o2::track::TrackParCovF>& nextCell,
  float bz,
  const CylinderCylinderPolicyParams& params) noexcept
{
  auto propagated = nextCell;
  if (!propagated.rotate(currentCell.getAlpha()) || !propagated.propagateTo(currentCell.getX(), bz)) {
    return false;
  }
  return currentCell.getPredictedChi2(propagated) <= params.maxChi2ClusterAttachment;
}

/// Disk formula: delegates to the existing MFT forward-state helper, which
/// already performs its own internal copy/propagation; neither input cell is
/// mutated.
template <>
GPUhdi() bool cellsAreCompatible<TransitionPolicyTag::DiskDisk>(
  const CellSeedTpl<o2::track::TrackParCovFwd>& currentCell,
  const CellSeedTpl<o2::track::TrackParCovFwd>& nextCell,
  float bz,
  const DiskDiskPolicyParams& params) noexcept
{
  return detail::mftFwdCellsAreCompatible(currentCell, nextCell, bz, params.maxChi2ClusterAttachment);
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYOPERATIONS_H_ */
