// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"

#include <cmath>
#include <cstdint>

#include "CommonConstants/MathConstants.h"
#include "GPUROOTSMatrixFwd.h"
#include <Math/SMatrix.h>

namespace o2::itsmft::tracking::barrel
{
namespace
{

using Matrix5 = float[5][5];

// Packed symmetric 5x5 covariance for the same-family stateChi2 combined
// matrix. o2::math_utils::MatRepSym<float, 5>::offset(row, column) matches
// packedCovarianceIndex bit-for-bit (both are row*(row+1)/2+column for
// row >= column), so the combined covariance is formed directly in packed
// storage without an intermediate dense unpack.
using CombinedCovariance = o2::math_utils::SMatrix<float, 5, 5, o2::math_utils::MatRepSym<float, 5>>;
static_assert(o2::math_utils::MatRepSym<float, 5>::kSize == 15, "packed symmetric 5x5 representation must hold exactly 15 floats");
static_assert(sizeof(CombinedCovariance) == 15 * sizeof(float), "combined covariance must occupy exactly 15 floats");

bool finiteState(const SurfaceKinematicState& state) noexcept
{
  if (!std::isfinite(state.referenceCoordinate) || !std::isfinite(state.alpha)) {
    return false;
  }
  for (const float value : state.parameters) {
    if (!std::isfinite(value)) {
      return false;
    }
  }
  for (const float value : state.covariance) {
    if (!std::isfinite(value)) {
      return false;
    }
  }
  return true;
}

bool finiteMeasurement(const SurfaceMeasurement& measurement) noexcept
{
  return std::isfinite(measurement.frame.u) && std::isfinite(measurement.frame.v) &&
         std::isfinite(measurement.covariance.uu) && std::isfinite(measurement.covariance.uv) &&
         std::isfinite(measurement.covariance.vv);
}

bool validateSource(const SurfaceKinematicState& state, OperationFailureReason& reason) noexcept
{
  if (state.family != StateFamily::Barrel) {
    reason = OperationFailureReason::SourceFamilyMismatch;
    return false;
  }
  if (!finiteState(state)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }
  return true;
}

void unpackCovariance(const SurfaceKinematicState& state, Matrix5& covariance) noexcept
{
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column < 5; ++column) {
      covariance[row][column] = state.covariance[packedCovarianceIndex(row, column)];
    }
  }
}

void packCovariance(const Matrix5& covariance, SurfaceKinematicState& state) noexcept
{
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = covariance[row][column];
    }
  }
}

void identity(Matrix5& matrix) noexcept
{
  for (uint8_t i = 0; i < 5; ++i) {
    matrix[i][i] = 1.f;
  }
}

void transportCovariance(SurfaceKinematicState& state, const Matrix5& jacobian) noexcept
{
  Matrix5 covariance{};
  Matrix5 product{};
  Matrix5 transported{};
  unpackCovariance(state, covariance);
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column < 5; ++column) {
      for (uint8_t inner = 0; inner < 5; ++inner) {
        product[row][column] += jacobian[row][inner] * covariance[inner][column];
      }
    }
  }
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column < 5; ++column) {
      for (uint8_t inner = 0; inner < 5; ++inner) {
        transported[row][column] += product[row][inner] * jacobian[column][inner];
      }
    }
  }
  packCovariance(transported, state);
}

bool commit(SurfaceKinematicState& destination, const SurfaceKinematicState& scratch, OperationFailureReason& reason) noexcept
{
  if (!finiteState(scratch)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }
  destination = scratch;
  return true;
}

bool residualInverse(const SurfaceKinematicState& state, const SurfaceMeasurement& measurement,
                     float& inverse00, float& inverse01, float& inverse11,
                     OperationFailureReason& reason, OperationFailureReason failure) noexcept
{
  const float s00 = state.covariance[packedCovarianceIndex(0, 0)] + measurement.covariance.uu;
  const float s01 = state.covariance[packedCovarianceIndex(1, 0)] + measurement.covariance.uv;
  const float s11 = state.covariance[packedCovarianceIndex(1, 1)] + measurement.covariance.vv;
  const float determinant = s00 * s11 - s01 * s01;
  if (!std::isfinite(determinant) || determinant == 0.f) {
    reason = std::isfinite(determinant) ? OperationFailureReason::InvalidCovariance : failure;
    return false;
  }
  const float inverseDeterminant = 1.f / determinant;
  inverse00 = s11 * inverseDeterminant;
  inverse01 = -s01 * inverseDeterminant;
  inverse11 = s00 * inverseDeterminant;
  if (!std::isfinite(inverse00) || !std::isfinite(inverse01) || !std::isfinite(inverse11)) {
    reason = failure;
    return false;
  }
  return true;
}

} // namespace

bool rotate(SurfaceKinematicState& state, float targetAlpha, OperationFailureReason& reason) noexcept
{
  if (!validateSource(state, reason) || !std::isfinite(targetAlpha)) {
    if (state.family == StateFamily::Barrel && finiteState(state)) {
      reason = OperationFailureReason::NonFiniteInput;
    }
    return false;
  }
  SurfaceKinematicState scratch = state;
  const float canonicalTargetAlpha = std::remainder(targetAlpha, 2.f * o2::constants::math::PI);
  const float delta = std::remainder(canonicalTargetAlpha - scratch.alpha, 2.f * o2::constants::math::PI);
  const float sine = std::sin(delta);
  const float cosine = std::cos(delta);
  const float snp = scratch.parameters[2];
  if (std::abs(snp) >= 1.f) {
    reason = OperationFailureReason::RotationFailure;
    return false;
  }
  const float csp = std::sqrt((1.f - snp) * (1.f + snp));
  const float rotatedCosine = csp * cosine + snp * sine;
  const float rotatedSnp = snp * cosine - csp * sine;
  if (rotatedCosine < 0.f || std::abs(rotatedSnp) >= 1.f || csp == 0.f) {
    reason = OperationFailureReason::RotationFailure;
    return false;
  }
  const float x = scratch.referenceCoordinate;
  const float y = scratch.parameters[0];
  scratch.referenceCoordinate = x * cosine + y * sine;
  scratch.parameters[0] = -x * sine + y * cosine;
  scratch.parameters[2] = rotatedSnp;
  scratch.alpha = canonicalTargetAlpha;
  const float ratio = cosine + snp / csp * sine;
  scratch.covariance[packedCovarianceIndex(0, 0)] *= cosine * cosine;
  scratch.covariance[packedCovarianceIndex(1, 0)] *= cosine;
  scratch.covariance[packedCovarianceIndex(2, 0)] *= cosine * ratio;
  scratch.covariance[packedCovarianceIndex(2, 1)] *= ratio;
  scratch.covariance[packedCovarianceIndex(2, 2)] *= ratio * ratio;
  scratch.covariance[packedCovarianceIndex(3, 0)] *= cosine;
  scratch.covariance[packedCovarianceIndex(3, 2)] *= ratio;
  scratch.covariance[packedCovarianceIndex(4, 0)] *= cosine;
  scratch.covariance[packedCovarianceIndex(4, 2)] *= ratio;
  return commit(state, scratch, reason);
}

bool propagate(SurfaceKinematicState& state, float targetX, float bz, OperationFailureReason& reason) noexcept
{
  if (!validateSource(state, reason) || !std::isfinite(targetX) || !std::isfinite(bz)) {
    if (state.family == StateFamily::Barrel && finiteState(state)) {
      reason = OperationFailureReason::NonFiniteInput;
    }
    return false;
  }
  SurfaceKinematicState scratch = state;
  const float dx = targetX - scratch.referenceCoordinate;
  if (dx == 0.f) {
    scratch.referenceCoordinate = targetX;
    return commit(state, scratch, reason);
  }
  const float snp = scratch.parameters[2];
  const float curvature = scratch.absCharge == 0 ? 0.f : scratch.parameters[4] * bz * o2::constants::math::B2C;
  const float propagatedSnp = snp + curvature * dx;
  if (std::abs(snp) >= 1.f || std::abs(propagatedSnp) >= 1.f) {
    reason = OperationFailureReason::UnreachableTarget;
    return false;
  }
  const float csp = std::sqrt((1.f - snp) * (1.f + snp));
  const float propagatedCsp = std::sqrt((1.f - propagatedSnp) * (1.f + propagatedSnp));
  if (csp == 0.f || propagatedCsp == 0.f) {
    reason = OperationFailureReason::UnreachableTarget;
    return false;
  }
  const float reciprocalCosines = 1.f / (csp + propagatedCsp);
  const float dyOverDx = (snp + propagatedSnp) * reciprocalCosines;
  const float x2r = curvature * dx;
  const bool arcZ = std::abs(x2r) > 0.05f;
  float dz = 0.f;
  if (arcZ) {
    const float argument = csp * propagatedSnp - propagatedCsp * snp;
    if (std::abs(argument) > 1.f || curvature == 0.f) {
      reason = OperationFailureReason::PropagationFailure;
      return false;
    }
    float angle = std::asin(argument);
    if (snp * snp + propagatedSnp * propagatedSnp > 1.f && snp * propagatedSnp < 0.f) {
      angle = propagatedSnp > 0.f ? o2::constants::math::PI - angle : -o2::constants::math::PI - angle;
    }
    dz = scratch.parameters[3] / curvature * angle;
  } else {
    dz = dx * (propagatedCsp + propagatedSnp * dyOverDx) * scratch.parameters[3];
  }
  scratch.referenceCoordinate = targetX;
  scratch.parameters[0] += dx * dyOverDx;
  scratch.parameters[1] += dz;
  scratch.parameters[2] = propagatedSnp;

  const float propagatedCspInverse = 1.f / propagatedCsp;
  const float dxOverCosines = dx * reciprocalCosines;
  const float hh = dxOverCosines * propagatedCspInverse * (1.f + csp * propagatedCsp + snp * propagatedSnp);
  const float jj = dx * (dyOverDx - propagatedSnp * propagatedCspInverse);
  Matrix5 jacobian{};
  identity(jacobian);
  jacobian[0][2] = hh / csp;
  jacobian[0][4] = hh * dxOverCosines * bz * o2::constants::math::B2C;
  jacobian[1][2] = scratch.parameters[3] * (jacobian[0][2] * propagatedSnp + jj);
  jacobian[1][3] = dx * (propagatedCsp + propagatedSnp * dyOverDx);
  jacobian[1][4] = scratch.parameters[3] * (jacobian[0][4] * propagatedSnp + jj * dx * bz * o2::constants::math::B2C);
  jacobian[2][4] = dx * bz * o2::constants::math::B2C;
  transportCovariance(scratch, jacobian);
  return commit(state, scratch, reason);
}

bool predictedChi2(const SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
                   OperationFailureReason& reason) noexcept
{
  if (!validateSource(state, reason) || !finiteMeasurement(measurement)) {
    if (state.family == StateFamily::Barrel && finiteState(state)) {
      reason = OperationFailureReason::NonFiniteInput;
    }
    return false;
  }
  float inverse00 = 0.f;
  float inverse01 = 0.f;
  float inverse11 = 0.f;
  if (!residualInverse(state, measurement, inverse00, inverse01, inverse11, reason,
                       OperationFailureReason::PredictedChi2Failure)) {
    return false;
  }
  const float residualY = measurement.frame.u - state.parameters[0];
  const float residualZ = measurement.frame.v - state.parameters[1];
  const float scratchChi2 = residualY * (inverse00 * residualY + inverse01 * residualZ) +
                            residualZ * (inverse01 * residualY + inverse11 * residualZ);
  if (!std::isfinite(scratchChi2)) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  chi2 = scratchChi2;
  return true;
}

bool update(SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
            OperationFailureReason& reason) noexcept
{
  if (!validateSource(state, reason) || !finiteMeasurement(measurement)) {
    if (state.family == StateFamily::Barrel && finiteState(state)) {
      reason = OperationFailureReason::NonFiniteInput;
    }
    return false;
  }
  float inverse00 = 0.f;
  float inverse01 = 0.f;
  float inverse11 = 0.f;
  if (!residualInverse(state, measurement, inverse00, inverse01, inverse11, reason, OperationFailureReason::UpdateFailure)) {
    return false;
  }
  Matrix5 covariance{};
  Matrix5 updatedCovariance{};
  float gain[5][2]{};
  unpackCovariance(state, covariance);
  const float residual[2] = {measurement.frame.u - state.parameters[0], measurement.frame.v - state.parameters[1]};
  SurfaceKinematicState scratch = state;
  for (uint8_t row = 0; row < 5; ++row) {
    gain[row][0] = covariance[row][0] * inverse00 + covariance[row][1] * inverse01;
    gain[row][1] = covariance[row][0] * inverse01 + covariance[row][1] * inverse11;
    scratch.parameters[row] += gain[row][0] * residual[0] + gain[row][1] * residual[1];
    for (uint8_t column = 0; column < 5; ++column) {
      updatedCovariance[row][column] = covariance[row][column] - gain[row][0] * covariance[0][column] - gain[row][1] * covariance[1][column];
    }
  }
  packCovariance(updatedCovariance, scratch);
  const float scratchChi2 = residual[0] * (inverse00 * residual[0] + inverse01 * residual[1]) +
                            residual[1] * (inverse01 * residual[0] + inverse11 * residual[1]);
  if (!std::isfinite(scratchChi2) || !finiteState(scratch)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }
  state = scratch;
  chi2 = scratchChi2;
  return true;
}

bool stateChi2(const SurfaceKinematicState& reference, const SurfaceKinematicState& candidate, float& chi2,
               OperationFailureReason& reason) noexcept
{
  if (reference.family != StateFamily::Barrel || candidate.family != StateFamily::Barrel) {
    reason = OperationFailureReason::SourceFamilyMismatch;
    return false;
  }
  if (!finiteState(reference) || !finiteState(candidate)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }
  if (std::abs(reference.alpha - candidate.alpha) > o2::constants::math::Epsilon) {
    reason = OperationFailureReason::AlphaMismatch;
    return false;
  }
  if (std::abs(reference.referenceCoordinate - candidate.referenceCoordinate) > o2::constants::math::Epsilon) {
    reason = OperationFailureReason::ReferenceCoordinateMismatch;
    return false;
  }

  CombinedCovariance combined;
  float* packed = combined.Array();
  for (uint8_t i = 0; i < 15; ++i) {
    packed[i] = reference.covariance[i] + candidate.covariance[i];
  }
  if (!combined.Invert()) {
    reason = OperationFailureReason::InvalidCovariance;
    return false;
  }

  float diff[5];
  for (uint8_t i = 0; i < 5; ++i) {
    diff[i] = reference.parameters[i] - candidate.parameters[i];
  }
  float chi2diag = 0.f;
  float chi2ndiag = 0.f;
  for (uint8_t i = 0; i < 5; ++i) {
    chi2diag += diff[i] * diff[i] * packed[packedCovarianceIndex(i, i)];
    for (uint8_t j = 0; j < i; ++j) {
      chi2ndiag += diff[i] * diff[j] * packed[packedCovarianceIndex(i, j)];
    }
  }
  const float scratchChi2 = chi2diag + 2.f * chi2ndiag;
  if (!std::isfinite(scratchChi2)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }
  chi2 = scratchChi2;
  return true;
}

} // namespace o2::itsmft::tracking::barrel
