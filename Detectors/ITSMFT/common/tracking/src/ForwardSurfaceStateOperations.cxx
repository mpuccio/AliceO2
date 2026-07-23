// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "CommonConstants/MathConstants.h"
#include "GPUROOTSMatrixFwd.h"
#include <Math/SMatrix.h>

#ifndef GPUCA_GPUCODE
#include "ITStracking/Cluster.h"
#endif

namespace o2::itsmft::tracking::forward
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
  return std::isfinite(measurement.global.x) && std::isfinite(measurement.global.y) &&
         std::isfinite(measurement.covariance.uu) && std::isfinite(measurement.covariance.uv) &&
         std::isfinite(measurement.covariance.vv);
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

void identity(Matrix5& matrix) noexcept
{
  for (uint8_t i = 0; i < 5; ++i) {
    matrix[i][i] = 1.f;
  }
}

bool validateSource(const SurfaceKinematicState& state, OperationFailureReason& reason) noexcept
{
  if (state.family != StateFamily::Forward) {
    reason = OperationFailureReason::SourceFamilyMismatch;
    return false;
  }
  if (!finiteState(state)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }
  return true;
}

bool commitPropagation(SurfaceKinematicState& destination, SurfaceKinematicState& scratch,
                       OperationFailureReason& reason) noexcept
{
  if (!finiteState(scratch)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }
  destination = scratch;
  return true;
}

bool propagateLinear(SurfaceKinematicState& state, float targetZ, OperationFailureReason& reason) noexcept
{
  const float dz = targetZ - state.referenceCoordinate;
  const float tanl = state.parameters[3];
  if (tanl == 0.f && dz != 0.f) {
    reason = OperationFailureReason::UnreachableTarget;
    return false;
  }
  if (dz == 0.f) {
    return true;
  }
  const float inverseTanl = 1.f / tanl;
  const float n = dz * inverseTanl;
  const float m = n * inverseTanl;
  const float sinPhi = std::sin(state.parameters[2]);
  const float cosPhi = std::cos(state.parameters[2]);
  state.parameters[0] += n * cosPhi;
  state.parameters[1] += n * sinPhi;
  state.referenceCoordinate = targetZ;

  Matrix5 jacobian{};
  identity(jacobian);
  jacobian[0][2] = -n * sinPhi;
  jacobian[0][3] = -m * cosPhi;
  jacobian[1][2] = n * cosPhi;
  jacobian[1][3] = -m * sinPhi;
  transportCovariance(state, jacobian);
  return true;
}

bool propagateQuadratic(SurfaceKinematicState& state, float targetZ, float bz,
                        OperationFailureReason& reason) noexcept
{
  const float dz = targetZ - state.referenceCoordinate;
  const float tanl = state.parameters[3];
  if (tanl == 0.f && dz != 0.f) {
    reason = OperationFailureReason::UnreachableTarget;
    return false;
  }
  if (dz == 0.f) {
    return true;
  }
  const float inverseTanl = 1.f / tanl;
  const float inverseQPt = state.parameters[4];
  const float sinPhi = std::sin(state.parameters[2]);
  const float cosPhi = std::cos(state.parameters[2]);
  const float fieldSign = std::copysign(1.f, bz);
  const float k = std::abs(o2::constants::math::B2C * bz);
  const float n = dz * inverseTanl;
  const float m = n * inverseTanl;
  const float theta = -inverseQPt * dz * k * inverseTanl;

  state.parameters[0] += n * cosPhi - 0.5f * n * theta * fieldSign * sinPhi;
  state.parameters[1] += n * sinPhi + 0.5f * n * theta * fieldSign * cosPhi;
  state.parameters[2] += fieldSign * theta;
  state.referenceCoordinate = targetZ;

  Matrix5 jacobian{};
  identity(jacobian);
  jacobian[0][2] = -0.5f * n * theta * fieldSign * cosPhi - n * sinPhi;
  jacobian[0][3] = fieldSign * m * theta * sinPhi - m * cosPhi;
  jacobian[0][4] = 0.5f * k * m * fieldSign * dz * sinPhi;
  jacobian[1][2] = -0.5f * n * theta * fieldSign * sinPhi + n * cosPhi;
  jacobian[1][3] = -fieldSign * m * theta * cosPhi - m * sinPhi;
  jacobian[1][4] = -0.5f * k * m * fieldSign * dz * cosPhi;
  jacobian[2][3] = -fieldSign * theta * inverseTanl;
  jacobian[2][4] = -fieldSign * k * n;
  transportCovariance(state, jacobian);
  return true;
}

bool propagateHelixParameters(SurfaceKinematicState& state, float targetZ, float bz,
                              OperationFailureReason& reason) noexcept
{
  const float dz = targetZ - state.referenceCoordinate;
  if (dz == 0.f) {
    return true;
  }
  const float tanl = state.parameters[3];
  const float inverseQPt = state.parameters[4];
  if (tanl == 0.f) {
    reason = OperationFailureReason::UnreachableTarget;
    return false;
  }
  if (bz == 0.f || inverseQPt == 0.f) {
    reason = OperationFailureReason::PropagationFailure;
    return false;
  }
  const float inverseTanl = 1.f / tanl;
  const float qPt = 1.f / inverseQPt;
  const float phi = state.parameters[2];
  const float sinPhi = std::sin(phi);
  const float cosPhi = std::cos(phi);
  const float k = std::abs(o2::constants::math::B2C * bz);
  const float inverseK = 1.f / k;
  const float theta = -inverseQPt * dz * k * inverseTanl;
  const float sinTheta = std::sin(theta);
  const float cosTheta = std::cos(theta);
  const float fieldSign = std::copysign(1.f, bz);
  const float y = sinPhi * qPt * inverseK;
  const float x = cosPhi * qPt * inverseK;
  state.parameters[0] += fieldSign * (y - y * cosTheta) - x * sinTheta;
  state.parameters[1] += fieldSign * (-x + x * cosTheta) - y * sinTheta;
  state.parameters[2] += fieldSign * theta;
  state.referenceCoordinate = targetZ;
  return true;
}

bool propagateHelix(SurfaceKinematicState& state, float targetZ, float bz,
                    OperationFailureReason& reason) noexcept
{
  const float originalZ = state.referenceCoordinate;
  const float dz = targetZ - originalZ;
  if (dz == 0.f) {
    return true;
  }
  const float phi = state.parameters[2];
  const float tanl = state.parameters[3];
  const float inverseQPt = state.parameters[4];
  if (!propagateHelixParameters(state, targetZ, bz, reason)) {
    return false;
  }
  const float inverseTanl = 1.f / tanl;
  const float qPt = 1.f / inverseQPt;
  const float sinPhi = std::sin(phi);
  const float cosPhi = std::cos(phi);
  const float k = std::abs(o2::constants::math::B2C * bz);
  const float inverseK = 1.f / k;
  const float theta = -inverseQPt * dz * k * inverseTanl;
  const float sinTheta = std::sin(theta);
  const float cosTheta = std::cos(theta);
  const float fieldSign = std::copysign(1.f, bz);
  const float n = dz * inverseTanl;
  const float m = n * inverseTanl;
  const float o = sinTheta * cosPhi;
  const float p = sinPhi * cosTheta;
  const float r = sinPhi * sinTheta;
  const float s = cosPhi * cosTheta;
  const float y = sinPhi * qPt * inverseK;
  const float x = cosPhi * qPt * inverseK;
  const float t = qPt * cosTheta;
  const float u = qPt * sinTheta;
  const float v = qPt;
  const float nn = dz * inverseTanl * qPt;

  Matrix5 jacobian{};
  identity(jacobian);
  jacobian[0][2] = fieldSign * x - fieldSign * x * cosTheta + y * sinTheta;
  jacobian[0][3] = fieldSign * r * m - s * m;
  jacobian[0][4] = -fieldSign * nn * r + fieldSign * t * y - fieldSign * v * y + nn * s + u * x;
  jacobian[1][2] = fieldSign * y - fieldSign * y * cosTheta - x * sinTheta;
  jacobian[1][3] = -fieldSign * o * m - p * m;
  jacobian[1][4] = fieldSign * nn * o - fieldSign * t * x + fieldSign * v * x + nn * p + u * y;
  jacobian[2][3] = -fieldSign * theta * inverseTanl;
  jacobian[2][4] = -fieldSign * k * n;
  transportCovariance(state, jacobian);
  return true;
}

template <PropagationModel Model>
bool propagateBound(SurfaceKinematicState& destination, float targetZ, float bz,
                    OperationFailureReason& reason) noexcept
{
  if (!validateSource(destination, reason) || !std::isfinite(targetZ) || !std::isfinite(bz)) {
    if (destination.family == StateFamily::Forward && finiteState(destination)) {
      reason = OperationFailureReason::NonFiniteInput;
    }
    return false;
  }
  SurfaceKinematicState scratch = destination;
  bool success = false;
  if constexpr (Model == PropagationModel::Linear) {
    success = propagateLinear(scratch, targetZ, reason);
  } else if constexpr (Model == PropagationModel::Quadratic) {
    success = propagateQuadratic(scratch, targetZ, bz, reason);
  } else if constexpr (Model == PropagationModel::Helix) {
    success = propagateHelix(scratch, targetZ, bz, reason);
  } else {
    if (bz == 0.f) {
      success = propagateLinear(scratch, targetZ, reason);
    } else {
      const SurfaceKinematicState original = scratch;
      success = propagateHelixParameters(scratch, targetZ, bz, reason);
      if (success) {
        SurfaceKinematicState covarianceScratch = original;
        success = propagateQuadratic(covarianceScratch, targetZ, bz, reason);
        if (success) {
          for (uint8_t i = 0; i < 15; ++i) {
            scratch.covariance[i] = covarianceScratch.covariance[i];
          }
        }
      }
    }
  }
  return success && commitPropagation(destination, scratch, reason);
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

template <>
bool propagate<PropagationModel::Linear>(SurfaceKinematicState& state, float targetZ, float bz,
                                         OperationFailureReason& reason) noexcept
{
  return propagateBound<PropagationModel::Linear>(state, targetZ, bz, reason);
}

template <>
bool propagate<PropagationModel::Quadratic>(SurfaceKinematicState& state, float targetZ, float bz,
                                            OperationFailureReason& reason) noexcept
{
  return propagateBound<PropagationModel::Quadratic>(state, targetZ, bz, reason);
}

template <>
bool propagate<PropagationModel::Helix>(SurfaceKinematicState& state, float targetZ, float bz,
                                        OperationFailureReason& reason) noexcept
{
  return propagateBound<PropagationModel::Helix>(state, targetZ, bz, reason);
}

template <>
bool propagate<PropagationModel::Optimized>(SurfaceKinematicState& state, float targetZ, float bz,
                                            OperationFailureReason& reason) noexcept
{
  return propagateBound<PropagationModel::Optimized>(state, targetZ, bz, reason);
}

bool predictedChi2(const SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
                   OperationFailureReason& reason) noexcept
{
  if (!validateSource(state, reason) || !finiteMeasurement(measurement)) {
    if (state.family == StateFamily::Forward && finiteState(state)) {
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
  const float residualX = measurement.global.x - state.parameters[0];
  const float residualY = measurement.global.y - state.parameters[1];
  const float scratchChi2 = residualX * (inverse00 * residualX + inverse01 * residualY) +
                            residualY * (inverse01 * residualX + inverse11 * residualY);
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
    if (state.family == StateFamily::Forward && finiteState(state)) {
      reason = OperationFailureReason::NonFiniteInput;
    }
    return false;
  }
  float inverse00 = 0.f;
  float inverse01 = 0.f;
  float inverse11 = 0.f;
  if (!residualInverse(state, measurement, inverse00, inverse01, inverse11, reason,
                       OperationFailureReason::UpdateFailure)) {
    return false;
  }

  Matrix5 covariance{};
  Matrix5 updatedCovariance{};
  float gain[5][2]{};
  unpackCovariance(state, covariance);
  const float residual[2] = {measurement.global.x - state.parameters[0], measurement.global.y - state.parameters[1]};
  SurfaceKinematicState scratch = state;
  for (uint8_t row = 0; row < 5; ++row) {
    gain[row][0] = covariance[row][0] * inverse00 + covariance[row][1] * inverse01;
    gain[row][1] = covariance[row][0] * inverse01 + covariance[row][1] * inverse11;
    scratch.parameters[row] += gain[row][0] * residual[0] + gain[row][1] * residual[1];
  }
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column < 5; ++column) {
      updatedCovariance[row][column] = covariance[row][column] -
                                       gain[row][0] * covariance[0][column] -
                                       gain[row][1] * covariance[1][column];
    }
  }
  packCovariance(updatedCovariance, scratch);
  const float scratchChi2 = residual[0] * (inverse00 * residual[0] + inverse01 * residual[1]) +
                            residual[1] * (inverse01 * residual[0] + inverse11 * residual[1]);
  if (!finiteState(scratch) || !std::isfinite(scratchChi2)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }
  state = scratch;
  chi2 = scratchChi2;
  return true;
}

bool correctForMaterial(SurfaceKinematicState& state, float xOverX0, OperationFailureReason& reason) noexcept
{
  if (!validateSource(state, reason) || !std::isfinite(xOverX0)) {
    if (state.family == StateFamily::Forward && finiteState(state)) {
      reason = OperationFailureReason::NonFiniteInput;
    }
    return false;
  }
  if (xOverX0 == 0.f) {
    return true;
  }
  const float tanl = state.parameters[3];
  if (tanl == 0.f) {
    reason = OperationFailureReason::MaterialFailure;
    return false;
  }
  const float inverseQPt = state.parameters[4];
  const float onePlusTanl2 = 1.f + tanl * tanl;
  const float inverseMomentum = std::abs(inverseQPt) / std::sqrt(onePlusTanl2);
  const float pathLengthOverX0 = xOverX0 * std::abs(std::sqrt(onePlusTanl2) / tanl);
  const float theta2 = highlandTheta2(inverseMomentum, pathLengthOverX0);
  SurfaceKinematicState scratch = state;
  scratch.covariance[packedCovarianceIndex(2, 2)] += theta2 * onePlusTanl2;
  scratch.covariance[packedCovarianceIndex(3, 3)] += theta2 * onePlusTanl2 * onePlusTanl2;
  scratch.covariance[packedCovarianceIndex(4, 4)] += theta2 * tanl * tanl * inverseQPt * inverseQPt;
  if (!finiteState(scratch)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }
  state = scratch;
  return true;
}

bool stateChi2(const SurfaceKinematicState& reference, const SurfaceKinematicState& candidate, float& chi2,
               OperationFailureReason& reason) noexcept
{
  if (reference.family != StateFamily::Forward || candidate.family != StateFamily::Forward) {
    reason = OperationFailureReason::SourceFamilyMismatch;
    return false;
  }
  if (!finiteState(reference) || !finiteState(candidate)) {
    reason = OperationFailureReason::NonFiniteInput;
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

  // Direct differences of (X, Y, Phi, Tanl, InvQPt): no phi wrapping. See
  // ITSMFTTracking/ForwardSurfaceStateOperations.h for the raw-phi rationale.
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

#ifndef GPUCA_GPUCODE

bool buildSeed(const o2::its::Cluster& clusterInner, const o2::its::Cluster& clusterMiddle, const o2::its::Cluster& clusterOuter,
               const o2::its::TrackingFrameInfo& hitOuter, float bz, float trackletMinPt,
               uint8_t absCharge, o2::track::PID pid,
               SurfaceKinematicState& outState, OperationFailureReason& reason) noexcept
{
  if (!std::isfinite(clusterInner.xCoordinate) || !std::isfinite(clusterInner.yCoordinate) || !std::isfinite(clusterInner.zCoordinate) ||
      !std::isfinite(clusterMiddle.xCoordinate) || !std::isfinite(clusterMiddle.yCoordinate) || !std::isfinite(clusterMiddle.zCoordinate) ||
      !std::isfinite(clusterOuter.xCoordinate) || !std::isfinite(clusterOuter.yCoordinate) || !std::isfinite(clusterOuter.zCoordinate) ||
      !std::isfinite(hitOuter.covarianceTrackingFrame[0]) || !std::isfinite(hitOuter.covarianceTrackingFrame[2]) ||
      !std::isfinite(bz) || !std::isfinite(trackletMinPt)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }

  // Strict boundary, transcribed verbatim from buildCellSeed<DiskDisk>
  // (TransitionPolicyOperations.cxx) / detail::mftFwdFitCellClusters
  // (MFTFwdTrackHelpers.h): established hard rejections, not non-finite-
  // output artifacts, so they are reported through the dedicated
  // SeedGeometryDegenerate reason.
  if (clusterInner.zCoordinate <= clusterOuter.zCoordinate + 1.e-6f) {
    reason = OperationFailureReason::SeedGeometryDegenerate;
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
    reason = OperationFailureReason::SeedGeometryDegenerate;
    return false;
  }

  // trackletMinPt<=0 fallback preserved verbatim (established behavior, not
  // re-validated here): invQPt becomes 0 rather than being rejected.
  const float invQPt = (trackletMinPt > 0.f) ? 1.f / trackletMinPt : 0.f;
  float tanl = 0.f;
  float phi = 0.f;
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

  SurfaceKinematicState scratch{};
  scratch.referenceCoordinate = clusterOuter.zCoordinate;
  scratch.alpha = 0.f;
  scratch.parameters[0] = clusterOuter.xCoordinate;
  scratch.parameters[1] = clusterOuter.yCoordinate;
  scratch.parameters[2] = phi;
  scratch.parameters[3] = tanl;
  scratch.parameters[4] = invQPt;
  // Only the diagonal is populated, matching the legacy seedCov (a
  // default-constructed, i.e. zero, SMatrix55Sym with just five diagonal
  // assignments) -- every off-diagonal packed entry stays at its
  // value-initialized 0.f.
  scratch.covariance[packedCovarianceIndex(0, 0)] = hitOuter.covarianceTrackingFrame[0] > 0.f ? hitOuter.covarianceTrackingFrame[0] : 1.f;
  scratch.covariance[packedCovarianceIndex(1, 1)] = hitOuter.covarianceTrackingFrame[2] > 0.f ? hitOuter.covarianceTrackingFrame[2] : 1.f;
  scratch.covariance[packedCovarianceIndex(2, 2)] = 1.f;
  scratch.covariance[packedCovarianceIndex(3, 3)] = 1.f;
  const float qptSigma = std::clamp(std::abs(invQPt), 1.f, 10.f);
  scratch.covariance[packedCovarianceIndex(4, 4)] = qptSigma * qptSigma;
  scratch.family = StateFamily::Forward;
  scratch.flags = 0;
  scratch.absCharge = absCharge;
  scratch.pid = pid;

  if (!finiteState(scratch)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }
  outState = scratch;
  return true;
}

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking::forward
