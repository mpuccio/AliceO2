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
#include <array>
#include <cmath>

#include "CommonConstants/MathConstants.h"
#include "DataFormatsITS/Vertex.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Constants.h"
#include "ITStracking/MathUtils.h"
#include "ITStracking/TrackHelpers.h"
#include "MFTTracking/Constants.h"

namespace o2::itsmft::tracking
{

bool TrackletSearchWindow<TransitionPolicyTag::CylinderCylinder>::acceptCandidate(
  const o2::its::Cluster& source,
  const o2::its::Cluster& target,
  float& tanLambdaOut) const
{
  const float deltaZ = o2::gpu::CAMath::Abs((tanLambda * (target.radius - source.radius)) + source.zCoordinate - target.zCoordinate);
  if (deltaZ / sigmaZ < nSigmaCut &&
      o2::its::math_utils::isPhiDifferenceBelow(source.phi, target.phi, phiCut)) {
    const float acceptedTanLambda = (source.zCoordinate - target.zCoordinate) / (source.radius - target.radius);
    tanLambdaOut = acceptedTanLambda;
    return true;
  }
  return false;
}

bool TrackletSearchWindow<TransitionPolicyTag::DiskDisk>::acceptCandidate(
  const o2::its::Cluster& source,
  const o2::its::Cluster& target,
  float& tanLambdaOut) const
{
  const float dx = target.xCoordinate - xProj;
  const float dy = target.yCoordinate - yProj;
  const float invSigmaX2 = (sigmaX > 0.f) ? 1.f / (sigmaX * sigmaX) : 0.f;
  const float invSigmaY2 = (sigmaY > 0.f) ? 1.f / (sigmaY * sigmaY) : 0.f;
  const float transChi2 = dx * dx * invSigmaX2 + dy * dy * invSigmaY2;
  if (transChi2 < o2::its::math_utils::Sq(nSigmaCut) && std::abs(meanDeltaZ) > 1.e-6f) {
    const float acceptedTanLambda = (source.zCoordinate - target.zCoordinate) / meanDeltaZ;
    tanLambdaOut = acceptedTanLambda;
    return true;
  }
  return false;
}

template <TransitionPolicyTag Tag, int NLayers>
bool projectSearchWindow(const o2::its::Cluster& source,
                         const o2::its::TrackingFrameInfo& sourceHit,
                         const o2::its::Vertex& vertex,
                         const TrackletProjectionState<Tag>& transitionState,
                         float bz,
                         const o2::itsmft::IndexTableUtils<NLayers>& indexUtils,
                         const typename TransitionPolicyTraits<Tag>::Params& params,
                         TrackletSearchWindow<Tag>& out)
{
  if constexpr (Tag == TransitionPolicyTag::CylinderCylinder) {
    const float inverseR0 = 1.f / source.radius;
    const float resolution = o2::gpu::CAMath::Sqrt(o2::its::math_utils::Sq(transitionState.sourcePositionResolution) +
                                                   o2::its::math_utils::Sq(params.pvResolution) / float(vertex.getNContributors()));
    const float tanLambda = (source.zCoordinate - vertex.getZ()) * inverseR0;
    const float zAtTargetMinR = tanLambda * (transitionState.targetMinR - source.radius) + source.zCoordinate;
    const float zAtTargetMaxR = tanLambda * (transitionState.targetMaxR - source.radius) + source.zCoordinate;
    const float sqInvDeltaZ0 = 1.f / (o2::its::math_utils::Sq(source.zCoordinate - vertex.getZ()) + o2::its::constants::Tolerance);
    const float sigmaZ = o2::gpu::CAMath::Sqrt((o2::its::math_utils::Sq(resolution) * o2::its::math_utils::Sq(tanLambda) *
                                                ((o2::its::math_utils::Sq(inverseR0) + sqInvDeltaZ0) * o2::its::math_utils::Sq(transitionState.meanDeltaR) + 1.f)) +
                                               o2::its::math_utils::Sq(transitionState.meanDeltaR * transitionState.transitionMSAngle));
    const auto bins = o2::itsmft::getBinsPhiZ(source.phi, transitionState.toLayer,
                                              zAtTargetMinR, zAtTargetMaxR,
                                              sigmaZ * params.nSigmaCut, transitionState.transitionPhiCut,
                                              indexUtils);
    if (bins.x < 0) {
      return false;
    }
    out = {bins, tanLambda, sigmaZ, transitionState.transitionPhiCut, params.nSigmaCut};
    return true;
  } else if constexpr (Tag == TransitionPolicyTag::DiskDisk) {
    float xProj = 0.f;
    float yProj = 0.f;
    detail::mftTrackletProject(source.xCoordinate, source.yCoordinate, source.zCoordinate,
                               vertex.getX(), vertex.getY(), vertex.getZ(),
                               transitionState.fromLayer, transitionState.toLayer, bz, params.trackletMinPt,
                               xProj, yProj);
    float sigmaX = 0.f;
    float sigmaY = 0.f;
    detail::mftTrackletSigmaXY(source.xCoordinate, source.yCoordinate,
                               vertex.getX(), vertex.getY(), vertex.getZ(),
                               sourceHit.covarianceTrackingFrame[0], sourceHit.covarianceTrackingFrame[2],
                               vertex.getSigmaX2(), vertex.getSigmaY2(), vertex.getSigmaZ2(),
                               transitionState.fromLayer, transitionState.toLayer,
                               transitionState.sourceReferenceRadius, transitionState.meanDeltaZ,
                               transitionState.transitionMSAngle, transitionState.transitionBendingAngle,
                               xProj, yProj, sigmaX, sigmaY);

    const float zSpread = params.nSigmaCut * vertex.getSigmaZ();
    const float zVtxMin = vertex.getZ() - zSpread;
    const float zVtxMax = vertex.getZ() + zSpread;
    const float absZFrom = std::abs(transitionState.fromZ);
    const float absZTo = std::abs(transitionState.toZ);
    const float denomMin = zVtxMax + absZFrom;
    const float denomMax = absZFrom + zVtxMin;
    float radialRangeMin = (std::abs(denomMin) > 1.e-6f) ? source.radius * (zVtxMax + absZTo) / denomMin : source.radius;
    float radialRangeMax = (std::abs(denomMax) > 1.e-6f) ? source.radius * (absZTo + zVtxMin) / denomMax : source.radius;
    if (radialRangeMin > radialRangeMax) {
      const float tmp = radialRangeMin;
      radialRangeMin = radialRangeMax;
      radialRangeMax = tmp;
    }
    const auto bins = o2::itsmft::getBinsRectClusterAtProj<NLayers>(xProj, yProj, transitionState.toLayer,
                                                                    radialRangeMin, radialRangeMax,
                                                                    sigmaX * params.nSigmaCut, sigmaY * params.nSigmaCut,
                                                                    indexUtils);
    if (bins.x < 0) {
      return false;
    }
    out = {bins, xProj, yProj, sigmaX, sigmaY, transitionState.meanDeltaZ, params.nSigmaCut};
    return true;
  } else {
    static_assert(Tag != Tag, "Unsupported transition policy tag");
  }
}

template bool projectSearchWindow<TransitionPolicyTag::CylinderCylinder, 7>(
  const o2::its::Cluster&, const o2::its::TrackingFrameInfo&, const o2::its::Vertex&,
  const TrackletProjectionState<TransitionPolicyTag::CylinderCylinder>&, float,
  const o2::itsmft::IndexTableUtils<7>&, const CylinderCylinderPolicyParams&,
  TrackletSearchWindow<TransitionPolicyTag::CylinderCylinder>&);

template bool projectSearchWindow<TransitionPolicyTag::DiskDisk, 10>(
  const o2::its::Cluster&, const o2::its::TrackingFrameInfo&, const o2::its::Vertex&,
  const TrackletProjectionState<TransitionPolicyTag::DiskDisk>&, float,
  const o2::itsmft::IndexTableUtils<10>&, const DiskDiskPolicyParams&,
  TrackletSearchWindow<TransitionPolicyTag::DiskDisk>&);

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
  // road pre-cut (passesCellRoadPrecut<DiskDisk>) is intentionally not
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

template <>
float layerMultipleScatteringAngle<TransitionPolicyTag::CylinderCylinder>(
  const LayerScatteringInputs<TransitionPolicyTag::CylinderCylinder>& inputs, float trackletMinPt)
{
  // Unchanged from the frozen ITS expression:
  // math_utils::MSangle(0.14f, trkParam.TrackletMinPt, trkParam.LayerxX0[iLayer]).
  return o2::its::math_utils::MSangle(0.14f, trackletMinPt, inputs.layerxX0);
}

template <>
float layerMultipleScatteringAngle<TransitionPolicyTag::DiskDisk>(
  const LayerScatteringInputs<TransitionPolicyTag::DiskDisk>& inputs, float trackletMinPt)
{
  // Same formula as the legacy detail::mftLayerMSAngle(), except zLayer/rRef
  // are supplied explicitly by the caller instead of being derived here from
  // mftLayerZ()/LayerZCoordinate() -- see the header doc on this
  // specialization and bindLegacyMFTReferenceCoordinates() below.
  const float invP = 1.f / trackletMinPt;
  const float zLayer = inputs.referenceCoordinate;
  const float rRef = inputs.layerRadius;
  const float tanlRef = (std::abs(rRef) > 1e-6f) ? zLayer / rRef : 0.f;
  const float absTanl = std::abs(tanlRef);
  const float cscLambda = (absTanl > 1e-6f) ? std::sqrt(1.f + tanlRef * tanlRef) / absTanl : 1e6f;
  return 0.0136f * invP * std::sqrt(inputs.layerxX0 * cscLambda);
}

namespace
{
// Time-boxed Gate 3 compatibility values: the legacy nominal half-disk z
// coordinates already used by mftLayerMSAngle() today, preserved bit-for-bit
// so layerMultipleScatteringAngle<DiskDisk> reproduces the accepted 91-track
// / hash 826dc653cd936a472929c600c97c140b baseline. Deliberately NOT
// SurfaceDescriptor::referenceCoordinate -- see the header doc on
// bindLegacyMFTReferenceCoordinates(). static constexpr storage duration:
// initialized at compile time, lives for the process lifetime, so a span
// over it is valid indefinitely and needs no per-iteration staging.
static constexpr std::array<float, o2::mft::constants::mft::LayersNumber> kLegacyMFTReferenceCoordinate =
  o2::mft::constants::mft::LayerZCoordinate();
} // namespace

DiskDiskReferenceCoordinateView bindLegacyMFTReferenceCoordinates() noexcept
{
  return DiskDiskReferenceCoordinateView{gsl::span<const float>(kLegacyMFTReferenceCoordinate)};
}

template <>
float clampTransitionCurvature<TransitionPolicyTag::CylinderCylinder>(float oneOverR, float r2) noexcept
{
  // Preserves the legacy double-promoted comparison verbatim (frozen ITS
  // TimeFrame.cxx / common-CA non-MFT branch): `0.5` is a double literal.
  return (0.5 * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
}

template <>
float clampTransitionCurvature<TransitionPolicyTag::DiskDisk>(float oneOverR, float r2) noexcept
{
  // Preserves the legacy float-only comparison verbatim (common-CA MFT
  // branch): `0.5f` stays in float.
  return (0.5f * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
}

TransitionScatteringBendingPrep prepareTransitionScatteringAndBending(
  gsl::span<const float> perLayerMSAngle, int fromLayer, int toLayer,
  float r1, float r2, float clampedOneOverR, float res1, float res2) noexcept
{
  float ms2 = 0.f;
  for (int layer = fromLayer; layer < toLayer; ++layer) {
    ms2 += o2::its::math_utils::Sq(perLayerMSAngle[layer]);
  }
  const float msAngle = o2::gpu::CAMath::Sqrt(ms2);
  const float cosTheta1half = o2::gpu::CAMath::Sqrt(1.f - o2::its::math_utils::Sq(0.5f * r1 * clampedOneOverR));
  const float cosTheta2half = o2::gpu::CAMath::Sqrt(1.f - o2::its::math_utils::Sq(0.5f * r2 * clampedOneOverR));
  const float x = (r2 * cosTheta1half) - (r1 * cosTheta2half);
  const float delta = o2::gpu::CAMath::Sqrt(1.f / (1.f - 0.25f * o2::its::math_utils::Sq(x * clampedOneOverR)) *
                                            (o2::its::math_utils::Sq((0.25f * r1 * r2 * o2::its::math_utils::Sq(clampedOneOverR) / cosTheta2half) + cosTheta1half) * o2::its::math_utils::Sq(res1) +
                                             o2::its::math_utils::Sq((0.25f * r1 * r2 * o2::its::math_utils::Sq(clampedOneOverR) / cosTheta1half) + cosTheta2half) * o2::its::math_utils::Sq(res2)));
  const float phiCut = o2::gpu::CAMath::Min(o2::gpu::CAMath::ASin(0.5f * x * clampedOneOverR) + 2.f * msAngle + delta, o2::constants::math::PI * 0.5f);
  return TransitionScatteringBendingPrep{msAngle, phiCut};
}

} // namespace o2::itsmft::tracking
