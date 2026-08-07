// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file TrackletFinding.cxx
/// \brief Out-of-line tracklet operation helpers
///
/// Only this translation unit may include MFTFwdTrackHelpers.h on behalf of
/// the D007 operation boundary, so the common public headers
/// stays free of MFT-specific constants, TimeFrame, and typed-output
/// dependencies.

#include "ITSMFTTracking/detail/TrackletFinding.h"
#include "ITSMFTTracking/detail/CellFinding.h"

#include <algorithm>
#include <array>
#include <cmath>

#include "CommonConstants/MathConstants.h"
#include "DataFormatsITS/Vertex.h"
#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Constants.h"
#include "ITStracking/MathUtils.h"
#include "ITStracking/TrackHelpers.h"
#include "MFTTracking/Constants.h"

namespace o2::itsmft::tracking
{

bool CylinderTrackletSearchWindow::acceptCandidate(
  const SurfaceMeasurement& sourceMeasurement,
  const o2::its::Cluster& sourceLocator,
  const SurfaceMeasurement& targetMeasurement,
  const o2::its::Cluster& targetLocator,
  float& tanLambdaOut) const
{
  const float deltaZ = o2::gpu::CAMath::Abs((tanLambda * (targetLocator.radius - sourceLocator.radius)) + sourceMeasurement.global.z - targetMeasurement.global.z);
  if (deltaZ / sigmaZ < nSigmaCut &&
      o2::its::math_utils::isPhiDifferenceBelow(sourceLocator.phi, targetLocator.phi, phiCut)) {
    const float acceptedTanLambda = (sourceMeasurement.global.z - targetMeasurement.global.z) / (sourceLocator.radius - targetLocator.radius);
    tanLambdaOut = acceptedTanLambda;
    return true;
  }
  return false;
}

bool DiskTrackletSearchWindow::acceptCandidate(
  const SurfaceMeasurement& sourceMeasurement,
  const SurfaceMeasurement& targetMeasurement,
  float& tanLambdaOut) const
{
  const float dx = targetMeasurement.global.x - xProj;
  const float dy = targetMeasurement.global.y - yProj;
  const float invSigmaX2 = (sigmaX > 0.f) ? 1.f / (sigmaX * sigmaX) : 0.f;
  const float invSigmaY2 = (sigmaY > 0.f) ? 1.f / (sigmaY * sigmaY) : 0.f;
  const float transChi2 = dx * dx * invSigmaX2 + dy * dy * invSigmaY2;
  if (transChi2 < o2::its::math_utils::Sq(nSigmaCut) && std::abs(meanDeltaZ) > 1.e-6f) {
    const float acceptedTanLambda = (sourceMeasurement.global.z - targetMeasurement.global.z) / meanDeltaZ;
    tanLambdaOut = acceptedTanLambda;
    return true;
  }
  return false;
}

bool projectCylinderSearchWindow(const SurfaceMeasurement& sourceMeasurement,
                                 const o2::its::Cluster& sourceLocator,
                                 const o2::its::Vertex& vertex,
                                 const CylinderTrackletProjectionState& transitionState,
                                 float /*bz*/, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                                 const TrackingKernelParameters& params,
                                 CylinderTrackletSearchWindow& out)
{
  const float inverseR0 = 1.f / sourceLocator.radius;
  const float resolution = o2::gpu::CAMath::Sqrt(o2::its::math_utils::Sq(transitionState.sourcePositionResolution) +
                                                 o2::its::math_utils::Sq(params.pvResolution) / float(vertex.getNContributors()));
  const float tanLambda = (sourceMeasurement.global.z - vertex.getZ()) * inverseR0;
  const float zAtTargetMinR = tanLambda * (transitionState.targetMinR - sourceLocator.radius) + sourceMeasurement.global.z;
  const float zAtTargetMaxR = tanLambda * (transitionState.targetMaxR - sourceLocator.radius) + sourceMeasurement.global.z;
  const float sqInvDeltaZ0 = 1.f / (o2::its::math_utils::Sq(sourceMeasurement.global.z - vertex.getZ()) + o2::its::constants::Tolerance);
  const float sigmaZ = o2::gpu::CAMath::Sqrt((o2::its::math_utils::Sq(resolution) * o2::its::math_utils::Sq(tanLambda) *
                                              ((o2::its::math_utils::Sq(inverseR0) + sqInvDeltaZ0) * o2::its::math_utils::Sq(transitionState.meanDeltaR) + 1.f)) +
                                             o2::its::math_utils::Sq(transitionState.meanDeltaR * transitionState.transitionMSAngle));
  const auto bins = o2::itsmft::getBinsPhiZ(sourceLocator.phi, transitionState.toLayer,
                                            zAtTargetMinR, zAtTargetMaxR,
                                            sigmaZ * params.nSigmaCut, transitionState.transitionPhiCut,
                                            indexUtils);
  if (bins.x < 0) {
    return false;
  }
  out = {bins, tanLambda, sigmaZ, transitionState.transitionPhiCut, params.nSigmaCut};
  return true;
}

bool projectDiskSearchWindow(const SurfaceMeasurement& sourceMeasurement,
                             const o2::its::Cluster& sourceLocator,
                             const o2::its::Vertex& vertex,
                             const DiskTrackletProjectionState& transitionState,
                             float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                             const TrackingKernelParameters& params,
                             DiskTrackletSearchWindow& out)
{
  float xProj = 0.f;
  float yProj = 0.f;
  detail::mftTrackletProject(sourceMeasurement.global.x, sourceMeasurement.global.y, sourceMeasurement.global.z,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             transitionState.fromLayer, transitionState.toLayer, bz, params.trackletMinPt,
                             xProj, yProj);
  float sigmaX = 0.f;
  float sigmaY = 0.f;
  detail::mftTrackletSigmaXY(sourceMeasurement.global.x, sourceMeasurement.global.y,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             sourceMeasurement.covariance.uu, sourceMeasurement.covariance.vv,
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
  float radialRangeMin = (std::abs(denomMin) > 1.e-6f) ? sourceLocator.radius * (zVtxMax + absZTo) / denomMin : sourceLocator.radius;
  float radialRangeMax = (std::abs(denomMax) > 1.e-6f) ? sourceLocator.radius * (absZTo + zVtxMin) / denomMax : sourceLocator.radius;
  if (radialRangeMin > radialRangeMax) {
    const float tmp = radialRangeMin;
    radialRangeMin = radialRangeMax;
    radialRangeMax = tmp;
  }
  const auto bins = o2::itsmft::getBinsRectClusterAtProj(xProj, yProj, transitionState.toLayer,
                                                         radialRangeMin, radialRangeMax,
                                                         sigmaX * params.nSigmaCut, sigmaY * params.nSigmaCut,
                                                         indexUtils);
  if (bins.x < 0) {
    return false;
  }
  out = {bins, xProj, yProj, sigmaX, sigmaY, transitionState.meanDeltaZ, params.nSigmaCut};
  return true;
}

// --- Native cylinder/disk cell leaves ---
//
// The cell leaves are composed entirely from the existing
// barrel::/forward:: primitives (BarrelSurfaceStateOperations.h/
// ForwardSurfaceStateOperations.h) and the shared PID/absCharge-aware
// material kernel (MaterialPhysics.h). See TrackletFinding.h for
// the per-operation contract documentation.

namespace
{

// "The accepted forward model" (Stage-B kickoff): reproduces
// detail::mftFwdPropagateToZ's own field-magnitude dispatch exactly --
// forward::propagate<Helix> when |bz| > 0.01f, otherwise
// forward::propagate<Linear> -- the same threshold and the same two models
// the legacy CA-construction path (detail::mftFwdAttachCluster, used by
// today's disk cell construction already uses.
// PropagationModel::Optimized (params: helix, errors: quadratic) is the
// separate model used only by the MFT final-track refit
// (TrackParCovFwd::propagateToZ, MFTTracking/TrackFitter.cxx) and is
// deliberately not used here.
bool forwardPropagateAcceptedModel(SurfaceKinematicState& state, float targetZ, float bz,
                                   OperationFailureReason& reason) noexcept
{
  if (std::abs(bz) > 0.01f) {
    return forward::propagate<forward::PropagationModel::Helix>(state, targetZ, bz, reason);
  }
  return forward::propagate<forward::PropagationModel::Linear>(state, targetZ, bz, reason);
}

} // namespace

bool buildCylinderCellSeed(
  const SurfaceMeasurement& measurementInner,
  const SurfaceMeasurement& measurementMiddle,
  const SurfaceMeasurement& measurementOuter,
  const std::array<NominalSurfaceMaterial, 3>& material,
  float bz,
  uint8_t absCharge,
  o2::track::PID pid,
  SurfaceKinematicState& outState,
  float& chi2,
  const TrackingKernelParameters& params,
  OperationFailureReason& reason) noexcept
{
  SurfaceKinematicState scratch{};
  if (!barrel::buildSeed(measurementInner, measurementMiddle, measurementOuter, bz, absCharge, pid, scratch, reason)) {
    return false;
  }

  float localChi2{0.f};
  const std::array<const SurfaceMeasurement*, 2> steps{&measurementMiddle, &measurementInner};
  const std::array<NominalSurfaceMaterial, 2> stepsMaterial{material[1], material[0]};
  for (int step = 0; step < 2; ++step) {
    const bool isLast = (step == 1);
    const auto& measurement = *steps[step];

    if (!barrel::rotate(scratch, measurement.frame.frameAngle, reason)) {
      return false;
    }
    if (!barrel::propagate(scratch, measurement.frame.q, bz, reason)) {
      return false;
    }
    const auto& stepMaterial = stepsMaterial[step];
    const auto materialResult = barrel::correctForMaterial(
      scratch, material::IntegratedMaterialBudget{stepMaterial.xOverX0, stepMaterial.arealDensityGPerCm2},
      material::MaterialTraversalDirection::OppositeMomentum);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    float predChi2{0.f};
    if (!barrel::predictedChi2(scratch, measurement, predChi2, reason)) {
      return false;
    }
    if (isLast && predChi2 > params.maxChi2ClusterAttachment) {
      reason = OperationFailureReason::PredictedChi2Failure;
      return false;
    }
    float updateChi2{0.f};
    if (!barrel::update(scratch, measurement, updateChi2, reason)) {
      return false;
    }
    localChi2 += updateChi2;
  }

  outState = scratch;
  chi2 = localChi2;
  return true;
}

bool buildDiskCellSeed(
  const SurfaceMeasurement& measurementInner,
  const SurfaceMeasurement& measurementMiddle,
  const SurfaceMeasurement& measurementOuter,
  const std::array<NominalSurfaceMaterial, 3>& material,
  float bz,
  uint8_t absCharge,
  o2::track::PID pid,
  SurfaceKinematicState& outState,
  float& chi2,
  const TrackingKernelParameters& params,
  OperationFailureReason& reason) noexcept
{
  SurfaceKinematicState scratch{};
  if (!forward::buildSeed(measurementInner, measurementMiddle, measurementOuter, bz, params.trackletMinPt,
                          absCharge, pid, scratch, reason)) {
    return false;
  }

  float localChi2{0.f};
  const std::array<const SurfaceMeasurement*, 3> steps{&measurementOuter, &measurementMiddle, &measurementInner};
  const std::array<NominalSurfaceMaterial, 3> stepsMaterial{material[2], material[1], material[0]};
  for (int step = 0; step < 3; ++step) {
    const bool isLast = (step == 2);
    const auto& measurement = *steps[step];

    if (!forwardPropagateAcceptedModel(scratch, measurement.frame.q, bz, reason)) {
      return false;
    }
    const auto& stepMaterial = stepsMaterial[step];
    const auto materialResult = forward::correctForMaterial(
      scratch, material::IntegratedMaterialBudget{stepMaterial.xOverX0, stepMaterial.arealDensityGPerCm2},
      material::MaterialTraversalDirection::OppositeMomentum);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    float predChi2{0.f};
    if (!forward::predictedChi2(scratch, measurement, predChi2, reason)) {
      return false;
    }
    if (isLast && predChi2 > params.maxChi2ClusterAttachment) {
      reason = OperationFailureReason::PredictedChi2Failure;
      return false;
    }
    float updateChi2{0.f};
    if (!forward::update(scratch, measurement, updateChi2, reason)) {
      return false;
    }
    localChi2 += updateChi2;
  }

  outState = scratch;
  chi2 = localChi2;
  return true;
}

bool attachCylinderHit(
  SurfaceKinematicState& state,
  const SurfaceMeasurement& measurement,
  const NominalSurfaceMaterial& material,
  float bz,
  float& chi2,
  const TrackingKernelParameters& params,
  OperationFailureReason& reason) noexcept
{
  SurfaceKinematicState scratch = state;
  float scratchChi2 = chi2;

  if (!barrel::rotate(scratch, measurement.frame.frameAngle, reason)) {
    return false;
  }
  if (!barrel::propagate(scratch, measurement.frame.q, bz, reason)) {
    return false;
  }
  const auto materialResult = barrel::correctForMaterial(
    scratch, material::IntegratedMaterialBudget{material.xOverX0, material.arealDensityGPerCm2},
    material::MaterialTraversalDirection::OppositeMomentum);
  if (!materialResult.ok()) {
    reason = OperationFailureReason::MaterialFailure;
    return false;
  }
  float predChi2{0.f};
  if (!barrel::predictedChi2(scratch, measurement, predChi2, reason)) {
    return false;
  }
  if (predChi2 > params.maxChi2ClusterAttachment || predChi2 < 0.f) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  float updateChi2{0.f};
  if (!barrel::update(scratch, measurement, updateChi2, reason)) {
    return false;
  }
  scratchChi2 += updateChi2;

  state = scratch;
  chi2 = scratchChi2;
  return true;
}

bool attachDiskHit(
  SurfaceKinematicState& state,
  const SurfaceMeasurement& measurement,
  const NominalSurfaceMaterial& material,
  float bz,
  float& chi2,
  const TrackingKernelParameters& params,
  OperationFailureReason& reason) noexcept
{
  SurfaceKinematicState scratch = state;
  float scratchChi2 = chi2;

  if (!forwardPropagateAcceptedModel(scratch, measurement.frame.q, bz, reason)) {
    return false;
  }
  const auto materialResult = forward::correctForMaterial(
    scratch, material::IntegratedMaterialBudget{material.xOverX0, material.arealDensityGPerCm2},
    material::MaterialTraversalDirection::OppositeMomentum);
  if (!materialResult.ok()) {
    reason = OperationFailureReason::MaterialFailure;
    return false;
  }
  float predChi2{0.f};
  if (!forward::predictedChi2(scratch, measurement, predChi2, reason)) {
    return false;
  }
  if (predChi2 > params.maxChi2ClusterAttachment) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  float updateChi2{0.f};
  if (!forward::update(scratch, measurement, updateChi2, reason)) {
    return false;
  }
  scratchChi2 += updateChi2;

  state = scratch;
  chi2 = scratchChi2;
  return true;
}

bool cellsCylinderAreCompatible(
  const SurfaceKinematicState& current,
  const SurfaceKinematicState& next,
  int /*currentSecondClusterIndex*/,
  int /*nextFirstClusterIndex*/,
  float bz,
  const TrackingKernelParameters& params) noexcept
{
  SurfaceKinematicState scratch = next;
  OperationFailureReason reason{};
  if (!barrel::rotate(scratch, current.alpha, reason) ||
      !barrel::propagate(scratch, current.referenceCoordinate, bz, reason)) {
    return false;
  }
  float chi2{0.f};
  if (!barrel::stateChi2(current, scratch, chi2, reason)) {
    return false;
  }
  return chi2 <= params.maxChi2ClusterAttachment;
}

bool cellsDiskAreCompatible(
  const SurfaceKinematicState& current,
  const SurfaceKinematicState& next,
  int currentSecondClusterIndex,
  int nextFirstClusterIndex,
  float bz,
  const TrackingKernelParameters& params) noexcept
{
  // Temporary Gate-3 compatibility input (see the header doc on the primary
  // template): checked first, exactly mirroring
  // detail::mftFwdCellsAreCompatible's own precedence.
  if (currentSecondClusterIndex != nextFirstClusterIndex) {
    return false;
  }
  SurfaceKinematicState scratch = next;
  OperationFailureReason reason{};
  if (!forwardPropagateAcceptedModel(scratch, current.referenceCoordinate, bz, reason)) {
    return false;
  }
  float chi2{0.f};
  if (!forward::stateChi2(current, scratch, chi2, reason)) {
    return false;
  }
  return chi2 <= params.maxChi2ClusterAttachment;
}

float cylinderLayerMultipleScatteringAngle(
  const CylinderLayerScatteringInputs& inputs, float trackletMinPt)
{
  // Unchanged from the frozen ITS expression:
  // math_utils::MSangle(0.14f, trkParam.TrackletMinPt, trkParam.LayerxX0[iLayer]).
  return o2::its::math_utils::MSangle(0.14f, trackletMinPt, inputs.layerxX0);
}

float diskLayerMultipleScatteringAngle(
  const DiskLayerScatteringInputs& inputs, float trackletMinPt)
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
// so layerMultipleScatteringAngle<Disk> reproduces the accepted 91-track
// / hash 826dc653cd936a472929c600c97c140b baseline. Deliberately NOT
// SurfaceDescriptor::referenceCoordinate -- see the header doc on
// bindLegacyMFTReferenceCoordinates(). static constexpr storage duration:
// initialized at compile time, lives for the process lifetime, so a span
// over it is valid indefinitely and needs no per-iteration staging.
static constexpr std::array<float, o2::mft::constants::mft::LayersNumber> kLegacyMFTReferenceCoordinate =
  o2::mft::constants::mft::LayerZCoordinate();
} // namespace

DiskReferenceCoordinateView bindLegacyMFTReferenceCoordinates() noexcept
{
  return DiskReferenceCoordinateView{gsl::span<const float>(kLegacyMFTReferenceCoordinate)};
}

float clampCylinderTransitionCurvature(float oneOverR, float r2) noexcept
{
  // Preserves the legacy double-promoted comparison verbatim (frozen ITS
  // TimeFrame.cxx / common-CA non-MFT branch): `0.5` is a double literal.
  return (0.5 * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
}

float clampDiskTransitionCurvature(float oneOverR, float r2) noexcept
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
