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

#include <algorithm>
#include <cmath>

#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Constants.h"
#include "ITStracking/TrackHelpers.h"

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

template <>
bool buildCellSeed<TransitionPolicyTag::CylinderCylinder>(
  const o2::its::Cluster& clusterInner,
  const o2::its::Cluster& clusterMiddle,
  const o2::its::Cluster& /*clusterOuter*/,
  const o2::its::TrackingFrameInfo& hitInner,
  const o2::its::TrackingFrameInfo& hitMiddle,
  const o2::its::TrackingFrameInfo& hitOuter,
  const std::array<float, 3>& xOverX0,
  float bz,
  o2::track::TrackParCovF& outState,
  float& chi2,
  const CylinderCylinderPolicyParams& params)
{
  // Matches TrackerTraits::computeLayerCells' barrel branch exactly: the
  // outer surface only enters through hitOuter inside buildTrackSeed, so
  // xOverX0[2] (outer) is never read here.
  auto track = o2::its::track::buildTrackSeed(clusterInner, clusterMiddle, hitOuter, bz);
  float localChi2{0.f};

  const std::array<const o2::its::TrackingFrameInfo*, 2> steps{&hitMiddle, &hitInner};
  const std::array<float, 2> stepsXOverX0{xOverX0[1], xOverX0[0]};
  bool good = false;
  for (int step = 0; step < 2; ++step) {
    const bool isLast = (step == 1);
    const auto& trackingHit = *steps[step];

    if (!track.rotate(trackingHit.alphaTrackingFrame)) {
      good = false;
      break;
    }
    if (!track.propagateTo(trackingHit.xTrackingFrame, bz)) {
      good = false;
      break;
    }
    const float x0 = stepsXOverX0[step];
    if (!track.correctForMaterial(x0, x0 * o2::its::constants::Radl * o2::its::constants::Rho, true)) {
      good = false;
      break;
    }
    const float predChi2 = track.getPredictedChi2Quiet(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame);
    if (isLast && predChi2 > params.maxChi2ClusterAttachment) {
      good = false;
      break;
    }
    if (!track.o2::track::TrackParCov::update(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)) {
      good = false;
      break;
    }
    localChi2 += predChi2;
    good = isLast;
  }

  if (!good) {
    return false;
  }
  outState = track;
  chi2 = localChi2;
  return true;
}

template <>
bool buildCellSeed<TransitionPolicyTag::DiskDisk>(
  const o2::its::Cluster& clusterInner,
  const o2::its::Cluster& clusterMiddle,
  const o2::its::Cluster& clusterOuter,
  const o2::its::TrackingFrameInfo& hitInner,
  const o2::its::TrackingFrameInfo& hitMiddle,
  const o2::its::TrackingFrameInfo& hitOuter,
  const std::array<float, 3>& xOverX0,
  float bz,
  o2::track::TrackParCovFwd& outState,
  float& chi2,
  const DiskDiskPolicyParams& params)
{
  // Matches detail::mftFwdFitCellClusters exactly, reading its three
  // clusters/hits directly instead of through a TimeFrame. The geometric
  // road pre-cut (detail::validateMFTCellClusters) is intentionally not
  // repeated here -- it depends on nominal layer position, not on these
  // measurements, and remains a TrackerTraits-owned guard.
  if (clusterInner.zCoordinate <= clusterOuter.zCoordinate + 1.e-6f) {
    return false;
  }

  const float dxTan = clusterMiddle.xCoordinate - clusterInner.xCoordinate;
  const float dyTan = clusterMiddle.yCoordinate - clusterInner.yCoordinate;
  const float dzTan = clusterMiddle.zCoordinate - clusterInner.zCoordinate;
  const float drTan = std::sqrt(dxTan * dxTan + dyTan * dyTan);
  const float dxPhi = clusterOuter.xCoordinate - clusterInner.xCoordinate;
  const float dyPhi = clusterOuter.yCoordinate - clusterInner.yCoordinate;
  const float dzPhi = clusterOuter.zCoordinate - clusterInner.zCoordinate;
  const float drPhi = std::sqrt(dxPhi * dxPhi + dyPhi * dyPhi);
  if (drTan < 1.e-6f || std::abs(dzTan) < 1.e-6f || drPhi < 1.e-6f || std::abs(dzPhi) < 1.e-6f) {
    return false;
  }

  const float minPt = params.trackletMinPt;
  const float invQPt = (minPt > 0.f) ? 1.f / minPt : 0.f;
  float tanl{0.f};
  float phi{0.f};
  if (std::abs(bz) > 0.01f) {
    tanl = -std::abs(dzTan) / drTan;
    phi = std::atan2(dyPhi, dxPhi);
    if (std::abs(tanl) > 1.e-6f) {
      const float k = std::abs(o2::constants::math::B2C * bz);
      const float hz = (bz > 0.f) ? 1.f : -1.f;
      phi -= 0.5f * hz * invQPt * dzPhi * k / tanl;
    }
  } else {
    tanl = -std::abs(dzPhi) / drPhi;
    phi = std::atan2(dyPhi, dxPhi);
  }

  ROOT::Math::SVector<double, 5> seedParams{clusterOuter.xCoordinate, clusterOuter.yCoordinate, phi, tanl, invQPt};
  ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>> seedCov{};
  seedCov(0, 0) = hitOuter.covarianceTrackingFrame[0] > 0.f ? hitOuter.covarianceTrackingFrame[0] : 1.f;
  seedCov(1, 1) = hitOuter.covarianceTrackingFrame[2] > 0.f ? hitOuter.covarianceTrackingFrame[2] : 1.f;
  seedCov(2, 2) = seedCov(3, 3) = 1.;
  const double qptSigma = std::clamp(static_cast<double>(std::abs(invQPt)), 1., 10.);
  seedCov(4, 4) = qptSigma * qptSigma;

  o2::track::TrackParCovFwd track{clusterOuter.zCoordinate, seedParams, seedCov, 0.};
  float localChi2{0.f};

  const std::array<const o2::its::TrackingFrameInfo*, 3> steps{&hitOuter, &hitMiddle, &hitInner};
  const std::array<float, 3> stepsXOverX0{xOverX0[2], xOverX0[1], xOverX0[0]};
  for (int step = 0; step < 3; ++step) {
    const bool isLast = (step == 2);
    const auto& tfInfo = *steps[step];
    if (!detail::mftFwdAttachCluster(track, tfInfo.zCoordinate, tfInfo.xCoordinate, tfInfo.yCoordinate,
                                     tfInfo.covarianceTrackingFrame[0], tfInfo.covarianceTrackingFrame[2],
                                     stepsXOverX0[step], bz, params.maxChi2ClusterAttachment, localChi2, isLast)) {
      return false;
    }
  }

  outState = track;
  chi2 = localChi2;
  return true;
}

} // namespace o2::itsmft::tracking
