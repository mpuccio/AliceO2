// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file TransitionPolicyOperations.cxx
/// \brief Out-of-line transition-policy operation specializations
///
/// Only this translation unit may include MFTFwdTrackHelpers.h on behalf of
/// the D007 policy operation boundary, so the common public policy header
/// stays free of MFT-specific constants/TimeFrame/MFTCATrack dependencies.

#include "ITSMFTTracking/TransitionPolicyOperations.h"

#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking
{

template <>
bool cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(
  const CellSeedTpl<o2::track::TrackParCovF>& currentCell,
  const CellSeedTpl<o2::track::TrackParCovF>& nextCell,
  float bz,
  const CylinderCylinderPolicyParams& params)
{
  auto propagated = nextCell;
  if (!propagated.rotate(currentCell.getAlpha()) || !propagated.propagateTo(currentCell.getX(), bz)) {
    return false;
  }
  return currentCell.getPredictedChi2(propagated) <= params.maxChi2ClusterAttachment;
}

template <>
bool cellsAreCompatible<TransitionPolicyTag::DiskDisk>(
  const CellSeedTpl<o2::track::TrackParCovFwd>& currentCell,
  const CellSeedTpl<o2::track::TrackParCovFwd>& nextCell,
  float bz,
  const DiskDiskPolicyParams& params)
{
  return detail::mftFwdCellsAreCompatible(currentCell, nextCell, bz, params.maxChi2ClusterAttachment);
}

template <>
bool attachHit<TransitionPolicyTag::CylinderCylinder>(
  o2::track::TrackParCovF& state,
  const o2::its::TrackingFrameInfo& hit,
  float xOverX0,
  o2::base::PropagatorF::MatCorrType corrType,
  float bz,
  float& chi2,
  const CylinderCylinderPolicyParams& params)
{
  if (!state.rotate(hit.alphaTrackingFrame)) {
    return false;
  }
  const auto* propagator = o2::base::Propagator::Instance();
  if (!propagator->propagateToX(state, hit.xTrackingFrame, bz,
                                o2::base::PropagatorImpl<float>::MAX_SIN_PHI,
                                o2::base::PropagatorImpl<float>::MAX_STEP, corrType)) {
    return false;
  }
  if (corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE &&
      !state.correctForMaterial(xOverX0, xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho, true)) {
    return false;
  }
  const float predictedChi2 = state.getPredictedChi2Quiet(hit.positionTrackingFrame, hit.covarianceTrackingFrame);
  if (predictedChi2 > params.maxChi2ClusterAttachment || predictedChi2 < 0.f) {
    return false;
  }
  chi2 += predictedChi2;
  return state.o2::track::TrackParCov::update(hit.positionTrackingFrame, hit.covarianceTrackingFrame);
}

template <>
bool attachHit<TransitionPolicyTag::DiskDisk>(
  o2::track::TrackParCovFwd& state,
  const o2::its::TrackingFrameInfo& hit,
  float xOverX0,
  o2::base::PropagatorF::MatCorrType,
  float bz,
  float& chi2,
  const DiskDiskPolicyParams& params)
{
  auto updatedState = state;
  float updatedChi2 = chi2;
  if (!detail::mftFwdAttachCluster(updatedState, hit.zCoordinate, hit.xCoordinate, hit.yCoordinate,
                                   hit.covarianceTrackingFrame[0], hit.covarianceTrackingFrame[2],
                                   xOverX0, bz, params.maxChi2ClusterAttachment, updatedChi2, true)) {
    return false;
  }
  state = updatedState;
  chi2 = updatedChi2;
  return true;
}

} // namespace o2::itsmft::tracking
