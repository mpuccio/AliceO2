// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "ITSMFTTracking/Propagator.h"

#include <cmath>

#include "ITSMFTTracking/detail/SurfaceStateOperations.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

namespace o2::itsmft::tracking
{

namespace
{

// Remove tiny negative diagonal values caused by floating-point cancellation
// during covariance transport. Larger negative values remain errors.
void clampNegligibleCovarianceNoise(SurfaceTrackState& state) noexcept
{
  constexpr float kNoiseFloor = 1.e-3f;
  for (uint8_t i = 0; i < 5; ++i) {
    const uint8_t index = packedCovarianceIndex(i, i);
    if (state.covariance[index] < 0.f && state.covariance[index] > -kNoiseFloor) {
      state.covariance[index] = 0.f;
    }
  }
}

// Apply outCov = J * inCov * J^T to a packed-symmetric 5x5 covariance.
void congruenceTransform(const float (&inCov)[15], const float (&jacobian)[5][5], float (&outCov)[15]) noexcept
{
  float full[5][5];
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t col = 0; col < 5; ++col) {
      full[row][col] = inCov[packedCovarianceIndex(row, col)];
    }
  }
  float tmp[5][5];
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t col = 0; col < 5; ++col) {
      float sum = 0.f;
      for (uint8_t k = 0; k < 5; ++k) {
        sum += jacobian[row][k] * full[k][col];
      }
      tmp[row][col] = sum;
    }
  }
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t col = 0; col <= row; ++col) {
      float sum = 0.f;
      for (uint8_t k = 0; k < 5; ++k) {
        sum += tmp[row][k] * jacobian[col][k];
      }
      outCov[packedCovarianceIndex(row, col)] = sum;
    }
  }
}

// Convert Barrel (bY, bZ, Snp, Tgl, Q2Pt) to Forward
// (X, Y, Phi, Tanl, InvQPt) at the state's current point.
bool barrelToForward(SurfaceTrackState& state, OperationFailureReason& reason) noexcept
{
  const float snp = state.parameters[2];
  if (!(std::abs(snp) < 1.f)) {
    reason = OperationFailureReason::SurfaceKindConversionFailure;
    return false;
  }
  const float csA = std::cos(state.alpha);
  const float snA = std::sin(state.alpha);
  const float csp = std::sqrt((1.f - snp) * (1.f + snp));
  const float bX = state.referenceCoordinate;
  const float bY = state.parameters[0];

  const float xGlo = bX * csA - bY * snA;
  const float yGlo = bX * snA + bY * csA;
  const float zGlo = state.parameters[1];
  float phi = std::remainder(state.alpha + std::asin(snp), o2::constants::math::TwoPI);
  // Match the library's (-pi, pi] angle convention.
  if (phi <= -o2::constants::math::PI) {
    phi += o2::constants::math::TwoPI;
  }

  const float jacobian[5][5] = {
    {-snA, 0.f, 0.f, 0.f, 0.f},
    {csA, 0.f, 0.f, 0.f, 0.f},
    {0.f, 0.f, 1.f / csp, 0.f, 0.f},
    {0.f, 0.f, 0.f, 1.f, 0.f},
    {0.f, 0.f, 0.f, 0.f, 1.f}};
  float newCov[15];
  congruenceTransform(state.covariance, jacobian, newCov);

  const float newParameters[5] = {xGlo, yGlo, phi, state.parameters[3], state.parameters[4]};
  for (uint8_t i = 0; i < 5; ++i) {
    state.parameters[i] = newParameters[i];
  }
  for (uint8_t i = 0; i < 15; ++i) {
    state.covariance[i] = newCov[i];
  }
  state.referenceCoordinate = zGlo;
  state.alpha = 0.f;
  state.kind = SurfaceKind::Disk;
  return true;
}

// Convert Forward (X, Y, Phi, Tanl, InvQPt) to Barrel
// (bY, bZ, Snp, Tgl, Q2Pt) at the current nominal point. The target alpha
// is fixed while constructing the linearized Jacobian. Since Forward has no
// variance for its reference z, the newly freed bZ uses kCZ2max.
bool forwardToBarrel(SurfaceTrackState& state, OperationFailureReason& reason) noexcept
{
  const float x = state.parameters[0];
  const float y = state.parameters[1];
  const float r = std::sqrt(x * x + y * y);
  if (!(r > 1.e-6f)) {
    reason = OperationFailureReason::SurfaceKindConversionFailure;
    return false;
  }
  const float alpha = std::atan2(y, x);
  const float csA = std::cos(alpha);
  const float snA = std::sin(alpha);
  const float phi = state.parameters[2];
  const float csp = std::cos(phi - alpha);
  const float snp = std::sin(phi - alpha);
  if (!(std::abs(snp) < 1.f)) {
    reason = OperationFailureReason::SurfaceKindConversionFailure;
    return false;
  }

  const float bX = x * csA + y * snA;
  const float bY = -x * snA + y * csA;
  const float bZ = state.referenceCoordinate;

  const float jacobian[5][5] = {
    {-snA, csA, 0.f, 0.f, 0.f},
    {0.f, 0.f, 0.f, 0.f, 0.f},
    {0.f, 0.f, csp, 0.f, 0.f},
    {0.f, 0.f, 0.f, 1.f, 0.f},
    {0.f, 0.f, 0.f, 0.f, 1.f}};
  float newCov[15];
  congruenceTransform(state.covariance, jacobian, newCov);
  newCov[packedCovarianceIndex(1, 1)] += o2::track::kCZ2max;

  const float newParameters[5] = {bY, bZ, snp, state.parameters[3], state.parameters[4]};
  for (uint8_t i = 0; i < 5; ++i) {
    state.parameters[i] = newParameters[i];
  }
  for (uint8_t i = 0; i < 15; ++i) {
    state.covariance[i] = newCov[i];
  }
  state.referenceCoordinate = bX;
  state.alpha = alpha;
  state.kind = SurfaceKind::Cylinder;
  return true;
}

} // namespace

bool Propagator::attachMeasurement(SurfaceTrackState& state, const SurfaceMeasurement& measurement,
                                   NominalSurfaceMaterial materialBudget, float bz,
                                   material::MaterialTraversalDirection direction,
                                   bool chi2GateEnabled, float maxChi2, float& chi2,
                                   OperationFailureReason& reason) noexcept
{
  if (chi2GateEnabled && maxChi2 < 0.f) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }

  SurfaceTrackState scratch = state;
  float predictedChi2 = 0.f;
  float updateChi2 = 0.f;
  const material::IntegratedMaterialBudget integratedMaterial{materialBudget.xOverX0, materialBudget.arealDensityGPerCm2};
  if (scratch.kind == SurfaceKind::Cylinder) {
    if (!detail::barrel::rotate(scratch, measurement.frame.frameAngle, reason) ||
        !detail::barrel::propagate(scratch, measurement.frame.q, bz, reason)) {
      return false;
    }
    const auto materialResult = detail::barrel::correctForMaterial(scratch, integratedMaterial, direction);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    if (!detail::barrel::predictedChi2(scratch, measurement, predictedChi2, reason)) {
      return false;
    }
    if (predictedChi2 < 0.f || (chi2GateEnabled && predictedChi2 > maxChi2)) {
      reason = OperationFailureReason::PredictedChi2Failure;
      return false;
    }
    if (!detail::barrel::update(scratch, measurement, updateChi2, reason)) {
      return false;
    }
  } else if (scratch.kind == SurfaceKind::Disk) {
    if (!propagateToReference(scratch, measurement.frame.q, bz, reason)) {
      return false;
    }
    const auto materialResult = detail::forward::correctForMaterial(scratch, integratedMaterial, direction);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    if (!detail::forward::predictedChi2(scratch, measurement, predictedChi2, reason)) {
      return false;
    }
    if (predictedChi2 < 0.f || (chi2GateEnabled && predictedChi2 > maxChi2)) {
      reason = OperationFailureReason::PredictedChi2Failure;
      return false;
    }
    if (!detail::forward::update(scratch, measurement, updateChi2, reason)) {
      return false;
    }
  } else {
    reason = OperationFailureReason::SourceSurfaceKindMismatch;
    return false;
  }
  const float updatedChi2 = chi2 + updateChi2;
  state = scratch;
  chi2 = updatedChi2;
  return true;
}

bool Propagator::stateChi2(const SurfaceTrackState& reference, const SurfaceTrackState& candidate,
                           float& chi2, OperationFailureReason& reason) noexcept
{
  if (reference.kind != candidate.kind) {
    reason = OperationFailureReason::SourceSurfaceKindMismatch;
    return false;
  }
  if (reference.kind == SurfaceKind::Cylinder) {
    return detail::barrel::stateChi2(reference, candidate, chi2, reason);
  }
  if (reference.kind == SurfaceKind::Disk) {
    return detail::forward::stateChi2(reference, candidate, chi2, reason);
  }
  reason = OperationFailureReason::SourceSurfaceKindMismatch;
  return false;
}

bool Propagator::propagateToReference(SurfaceTrackState& state, float targetReferenceCoordinate, float bz,
                                      OperationFailureReason& reason) noexcept
{
  if (state.kind == SurfaceKind::Cylinder) {
    return detail::barrel::propagate(state, targetReferenceCoordinate, bz, reason);
  }
  if (state.kind == SurfaceKind::Disk) {
    return detail::forward::propagate(state, targetReferenceCoordinate, bz, reason);
  }
  reason = OperationFailureReason::SourceSurfaceKindMismatch;
  return false;
}

bool Propagator::propagateToReference(SurfaceTrackState& state, SurfaceTrackParameters& linRef,
                                      float targetReferenceCoordinate, float bz,
                                      OperationFailureReason& reason) noexcept
{
  if (state.kind != linRef.kind) {
    reason = OperationFailureReason::SourceSurfaceKindMismatch;
    return false;
  }
  if (state.kind == SurfaceKind::Cylinder) {
    return detail::barrel::propagate(state, linRef, targetReferenceCoordinate, bz, reason);
  }
  if (state.kind == SurfaceKind::Disk) {
    return detail::forward::propagate(state, linRef, targetReferenceCoordinate, bz, reason);
  }
  reason = OperationFailureReason::SourceSurfaceKindMismatch;
  return false;
}

bool Propagator::convertKind(SurfaceTrackState& state, SurfaceKind targetKind,
                             OperationFailureReason& reason) noexcept
{
  if (targetKind != SurfaceKind::Cylinder && targetKind != SurfaceKind::Disk) {
    reason = OperationFailureReason::SurfaceKindConversionFailure;
    return false;
  }
  if (state.kind != SurfaceKind::Cylinder && state.kind != SurfaceKind::Disk) {
    reason = OperationFailureReason::SourceSurfaceKindMismatch;
    return false;
  }
  if (state.kind == targetKind) {
    return true;
  }
  if (targetKind == SurfaceKind::Disk) {
    return barrelToForward(state, reason);
  }
  return forwardToBarrel(state, reason);
}

bool Propagator::propagateToMeasurement(SurfaceTrackState& state, SurfaceTrackParameters& linRef,
                                        const SurfaceDescriptor& targetSurface, const SurfaceMeasurement& targetMeasurement,
                                        float bz, material::MaterialTraversalDirection direction,
                                        bool chi2GateEnabled, float maxChi2, float& chi2,
                                        bool shiftReferenceToMeasurement, OperationFailureReason& reason) noexcept
{
  if (chi2 < 0.f) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  if (chi2GateEnabled && maxChi2 < 0.f) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }

  const SurfaceKind targetKind = targetSurface.kind;
  if (targetKind == SurfaceKind::Undefined) {
    reason = OperationFailureReason::SurfaceKindConversionFailure;
    return false;
  }

  SurfaceTrackState scratchState = state;
  SurfaceTrackParameters scratchRef = linRef;

  if (scratchState.kind != targetKind) {
    if (!convertKind(scratchState, targetKind, reason)) {
      return false;
    }
    // Changing parameter conventions is also a relinearization boundary.
    // The conversion Jacobian is evaluated at scratchState, so begin the
    // target-kind propagation from that same point.
    scratchRef = SurfaceTrackParameters{scratchState};
  }

  const material::IntegratedMaterialBudget materialBudget{targetSurface.material.xOverX0, targetSurface.material.arealDensityGPerCm2};
  float scratchChi2 = chi2;
  float predChi2 = 0.f;
  float updateChi2 = 0.f;

  if (targetKind == SurfaceKind::Cylinder) {
    if (!detail::barrel::rotate(scratchState, scratchRef, targetMeasurement.frame.frameAngle, bz, reason)) {
      return false;
    }
    if (!detail::barrel::propagate(scratchState, scratchRef, targetMeasurement.frame.q, bz, reason)) {
      return false;
    }
    clampNegligibleCovarianceNoise(scratchState);
    const auto materialResult = detail::barrel::correctForMaterial(scratchState, materialBudget, direction);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    if (!detail::barrel::predictedChi2(scratchState, targetMeasurement, predChi2, reason)) {
      return false;
    }
  } else {
    if (!Propagator::propagateToReference(scratchState, scratchRef, targetMeasurement.frame.q, bz, reason)) {
      return false;
    }
    clampNegligibleCovarianceNoise(scratchState);
    const auto materialResult = detail::forward::correctForMaterial(scratchState, materialBudget, direction);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    if (!detail::forward::predictedChi2(scratchState, targetMeasurement, predChi2, reason)) {
      return false;
    }
  }

  if (predChi2 < 0.f) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  if (chi2GateEnabled && predChi2 > maxChi2) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }

  if (targetKind == SurfaceKind::Cylinder) {
    if (!detail::barrel::update(scratchState, targetMeasurement, updateChi2, reason)) {
      return false;
    }
  } else {
    if (!detail::forward::update(scratchState, targetMeasurement, updateChi2, reason)) {
      return false;
    }
  }
  scratchChi2 += updateChi2;
  if (scratchChi2 < 0.f) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }

  if (shiftReferenceToMeasurement) {
    if (targetKind == SurfaceKind::Cylinder) {
      if (!detail::barrel::shiftReferenceToMeasurement(scratchRef, targetMeasurement, reason)) {
        return false;
      }
    } else {
      if (!detail::forward::shiftReferenceToMeasurement(scratchRef, targetMeasurement, reason)) {
        return false;
      }
    }
  }

  state = scratchState;
  linRef = scratchRef;
  chi2 = scratchChi2;
  return true;
}

} // namespace o2::itsmft::tracking
