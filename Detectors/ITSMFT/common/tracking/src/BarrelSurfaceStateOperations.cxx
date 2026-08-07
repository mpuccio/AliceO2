// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "CommonConstants/MathConstants.h"
#include "GPUROOTSMatrixFwd.h"
#include <Math/SMatrix.h>

// TrackParametrization.h is included solely to reuse its public
// kCY2max/kCZ2max/kCSnp2max/kCTgl2max/kC1Pt2max/kMostProbablePt constants,
// matching the exact covariance/curvature conventions
// o2::its::track::buildTrackSeed (ITStracking/TrackHelpers.h) uses and the
// covariance-validity upper range sanitizeCovariance() (SurfaceKinematicState.h)
// enforces below -- no narrower public header declares them (the same reuse
// pattern MaterialPhysics.cxx already uses for its own constants). Not part
// of this translation unit's public interface, and nothing here constructs
// or references a TrackParametrization/TrackParCov object. Unconditional
// (not GPUCA_GPUCODE-guarded): unlike ITStracking/MathUtils.h below (needed
// only by the host-only buildSeed), these constants are also needed by the
// device-visible rotate/propagate/update operations' sanitization call.
#include "ReconstructionDataFormats/TrackParametrization.h"

#ifndef GPUCA_GPUCODE
#include "ITStracking/MathUtils.h"
#endif

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

// Upper bound for sanitizeCovariance() (SurfaceKinematicState.h), in (Y, Z,
// Snp, Tgl, Q2Pt) slot order -- the exact five constants
// FamilyMaterialOperations.cxx's barrel covariance-range limiting already
// uses, so post-propagate/rotate/update sanitization and the material-step
// range clamp enforce the identical upper bound.
constexpr float kBarrelMaxDiagonal[5] = {o2::track::kCY2max, o2::track::kCZ2max, o2::track::kCSnp2max,
                                         o2::track::kCTgl2max, o2::track::kC1Pt2max};

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

// Shared commit point for rotate() and propagate() (non-linRef) above: every
// exit from either function funnels through here, so sanitizing the
// covariance-validity invariant (ADR 0008) unconditionally here -- once --
// covers both operations and every one of their exit paths (including the
// dx==0 trivial-step early return, where covariance is unchanged and
// sanitization is a no-op) without relying on each call site to remember it.
bool commit(SurfaceKinematicState& destination, SurfaceKinematicState& scratch, OperationFailureReason& reason) noexcept
{
  if (!finiteState(scratch)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }
  sanitizeCovariance(scratch, kBarrelMaxDiagonal);
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
  // ADR 0008: the naive/non-Joseph-form Kalman covariance subtraction above
  // can reveal an already out-of-bounds correlation (introduced upstream by
  // a large-step propagate) as a small negative diagonal; sanitize
  // unconditionally before committing so no caller ever observes it.
  sanitizeCovariance(scratch, kBarrelMaxDiagonal);
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

#ifndef GPUCA_GPUCODE

namespace
{

// Shared closed-form seed construction, transcribed verbatim from
// o2::its::track::buildTrackSeed (ITStracking/TrackHelpers.h) with its
// `cluster1`/`cluster2`/`tf3`/`reverse` arguments renamed to
// `clusterA`/`clusterB`/`frameMeasurement`/`sign` respectively.
bool buildSeedImpl(const SurfaceMeasurement& clusterA, const SurfaceMeasurement& clusterB,
                   const SurfaceMeasurement& frameMeasurement, float bz, float sign,
                   uint8_t absCharge, o2::track::PID pid,
                   SurfaceKinematicState& outState, OperationFailureReason& reason) noexcept
{
  if (!std::isfinite(clusterA.global.x) || !std::isfinite(clusterA.global.y) || !std::isfinite(clusterA.global.z) ||
      !std::isfinite(clusterB.global.x) || !std::isfinite(clusterB.global.y) || !std::isfinite(clusterB.global.z) ||
      !std::isfinite(frameMeasurement.frame.q) || !std::isfinite(frameMeasurement.frame.frameAngle) ||
      !std::isfinite(frameMeasurement.frame.u) || !std::isfinite(frameMeasurement.frame.v) ||
      !std::isfinite(frameMeasurement.covariance.uu) || !std::isfinite(frameMeasurement.covariance.uv) ||
      !std::isfinite(frameMeasurement.covariance.vv) || !std::isfinite(bz)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }

  const float cosAlpha = std::cos(frameMeasurement.frame.frameAngle);
  const float sinAlpha = std::sin(frameMeasurement.frame.frameAngle);
  const float x1 = (clusterA.global.x * cosAlpha) + (clusterA.global.y * sinAlpha);
  const float y1 = (-clusterA.global.x * sinAlpha) + (clusterA.global.y * cosAlpha);
  const float x2 = (clusterB.global.x * cosAlpha) + (clusterB.global.y * sinAlpha);
  const float y2 = (-clusterB.global.x * sinAlpha) + (clusterB.global.y * cosAlpha);
  const float x3 = frameMeasurement.frame.q;
  const float y3 = frameMeasurement.frame.u;

  float snp = 0.f;
  float q2pt = 0.f;
  float q2pt2 = 0.f;
  if (std::abs(bz) < 0.01f) { // zero field
    const float dx = x3 - x1;
    const float dy = y3 - y1;
    snp = sign * dy / std::hypot(dx, dy);
    q2pt = 1.f / o2::track::kMostProbablePt;
    q2pt2 = 1.f;
  } else {
    const float crv = o2::its::math_utils::computeCurvature(x3, y3, x2, y2, x1, y1);
    snp = sign * crv * (x3 - o2::its::math_utils::computeCurvatureCentreX(x3, y3, x2, y2, x1, y1));
    q2pt = sign * crv / (bz * o2::constants::math::B2C);
    q2pt2 = crv * crv;
  }
  const float tgl = -0.5f * sign * (o2::its::math_utils::computeTanDipAngle(x1, y1, x2, y2, clusterA.global.z, clusterB.global.z) + o2::its::math_utils::computeTanDipAngle(x2, y2, x3, y3, clusterB.global.z, frameMeasurement.frame.v));
  const float sg2q2pt = o2::track::kC1Pt2max * std::clamp(q2pt2, 0.0005f, 1.0f);

  SurfaceKinematicState scratch{};
  scratch.referenceCoordinate = x3;
  scratch.alpha = frameMeasurement.frame.frameAngle;
  scratch.parameters[0] = y3;
  scratch.parameters[1] = frameMeasurement.frame.v;
  scratch.parameters[2] = snp;
  scratch.parameters[3] = tgl;
  scratch.parameters[4] = q2pt;
  // Packed layout matches o2::track::TrackParametrizationWithError's own
  // packed-symmetric storage index-for-index (both are row*(row+1)/2+column
  // for row >= column), so the legacy flat covariance initializer list
  // transcribes directly, entry for entry, with no reshuffling.
  scratch.covariance[packedCovarianceIndex(0, 0)] = frameMeasurement.covariance.uu;
  scratch.covariance[packedCovarianceIndex(1, 0)] = frameMeasurement.covariance.uv;
  scratch.covariance[packedCovarianceIndex(1, 1)] = frameMeasurement.covariance.vv;
  scratch.covariance[packedCovarianceIndex(2, 0)] = 0.f;
  scratch.covariance[packedCovarianceIndex(2, 1)] = 0.f;
  scratch.covariance[packedCovarianceIndex(2, 2)] = o2::track::kCSnp2max;
  scratch.covariance[packedCovarianceIndex(3, 0)] = 0.f;
  scratch.covariance[packedCovarianceIndex(3, 1)] = 0.f;
  scratch.covariance[packedCovarianceIndex(3, 2)] = 0.f;
  scratch.covariance[packedCovarianceIndex(3, 3)] = o2::track::kCTgl2max;
  scratch.covariance[packedCovarianceIndex(4, 0)] = 0.f;
  scratch.covariance[packedCovarianceIndex(4, 1)] = 0.f;
  scratch.covariance[packedCovarianceIndex(4, 2)] = 0.f;
  scratch.covariance[packedCovarianceIndex(4, 3)] = 0.f;
  scratch.covariance[packedCovarianceIndex(4, 4)] = sg2q2pt;
  scratch.family = StateFamily::Barrel;
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

// Non-covariance propagation of a SurfaceLinearizationReference, transcribed
// from o2::track::TrackParametrization<float>::propagateParamTo(xk, b)
// (DataFormats/Reconstruction/src/TrackParametrization.cxx) -- the exact
// parameter-space formula the non-linRef `propagate` above already
// transcribes for SurfaceKinematicState, applied here to the covariance-free
// reference. `stateAbsCharge` substitutes for the legacy reference's own
// absCharge field (a SurfaceLinearizationReference carries none -- see the
// paired type in SurfaceKinematicState.h): the reference and its paired state
// always describe the same particle hypothesis.
bool propagateReferenceParams(SurfaceLinearizationReference& ref, uint8_t stateAbsCharge, float targetX, float bz,
                              OperationFailureReason& reason) noexcept
{
  const float dx = targetX - ref.referenceCoordinate;
  if (dx == 0.f) {
    ref.referenceCoordinate = targetX;
    return true;
  }
  const float snp = ref.parameters[2];
  const float curvature = stateAbsCharge == 0 ? 0.f : ref.parameters[4] * bz * o2::constants::math::B2C;
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
    dz = ref.parameters[3] / curvature * angle;
  } else {
    dz = dx * (propagatedCsp + propagatedSnp * dyOverDx) * ref.parameters[3];
  }
  ref.referenceCoordinate = targetX;
  ref.parameters[0] += dx * dyOverDx;
  ref.parameters[1] += dz;
  ref.parameters[2] = propagatedSnp;
  return true;
}

bool finiteLinRef(const SurfaceLinearizationReference& ref) noexcept
{
  if (!std::isfinite(ref.referenceCoordinate) || !std::isfinite(ref.alpha)) {
    return false;
  }
  for (const float value : ref.parameters) {
    if (!std::isfinite(value)) {
      return false;
    }
  }
  return true;
}

} // namespace

bool buildSeed(const SurfaceMeasurement& measurementInner, const SurfaceMeasurement& measurementMiddle,
               const SurfaceMeasurement& measurementOuter, float bz,
               uint8_t absCharge, o2::track::PID pid,
               SurfaceKinematicState& outState, OperationFailureReason& reason) noexcept
{
  return buildSeedImpl(measurementInner, measurementMiddle, measurementOuter, bz, 1.f, absCharge, pid, outState, reason);
}

bool rotate(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef, float targetAlpha, float bz,
            OperationFailureReason& reason) noexcept
{
  if (!validateSource(state, reason)) {
    return false;
  }
  if (linRef.family != StateFamily::Barrel) {
    reason = OperationFailureReason::SourceFamilyMismatch;
    return false;
  }
  if (!finiteLinRef(linRef) || !std::isfinite(targetAlpha) || !std::isfinite(bz)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }
  // Fitted-state/linRef pairing precondition: parameters may legitimately
  // differ (that is the entire purpose of a linearization reference), but
  // the anchor and frame may not. makeLinearizationReference and every
  // successful paired rotate/propagate establish referenceCoordinate/alpha
  // identically (bit-for-bit, not merely within tolerance), so this is an
  // exact comparison, not an epsilon check.
  if (state.referenceCoordinate != linRef.referenceCoordinate) {
    reason = OperationFailureReason::ReferenceCoordinateMismatch;
    return false;
  }
  if (state.alpha != linRef.alpha) {
    reason = OperationFailureReason::AlphaMismatch;
    return false;
  }
  const float stateSnp = state.parameters[2];
  if (std::abs(stateSnp) >= 1.f) {
    reason = OperationFailureReason::RotationFailure;
    return false;
  }

  SurfaceKinematicState scratchState = state;
  SurfaceLinearizationReference scratchRef = linRef;

  const float canonicalAlpha = std::remainder(targetAlpha, 2.f * o2::constants::math::PI);

  // Rotate the reference (linRef1.rotateParam(alpha, ca, sa)): gated by the
  // reference's own pre-rotation snp, not state's.
  const float refSnpBefore = scratchRef.parameters[2];
  if (std::abs(refSnpBefore) >= 1.f) {
    reason = OperationFailureReason::RotationFailure;
    return false;
  }
  const float delta = std::remainder(canonicalAlpha - scratchRef.alpha, 2.f * o2::constants::math::PI);
  const float sa = std::sin(delta);
  const float ca = std::cos(delta);
  const float refCsp0 = std::sqrt((1.f - refSnpBefore) * (1.f + refSnpBefore));
  if (refCsp0 * ca + refSnpBefore * sa < 0.f) {
    reason = OperationFailureReason::RotationFailure;
    return false;
  }
  const float refSnpRotated = refSnpBefore * ca - refCsp0 * sa;
  if (std::abs(refSnpRotated) >= 1.f) {
    reason = OperationFailureReason::RotationFailure;
    return false;
  }
  const float refXOld = scratchRef.referenceCoordinate;
  const float refYOld = scratchRef.parameters[0];
  scratchRef.alpha = canonicalAlpha;
  scratchRef.referenceCoordinate = refXOld * ca + refYOld * sa;
  scratchRef.parameters[0] = -refXOld * sa + refYOld * ca;
  scratchRef.parameters[2] = refSnpRotated;

  // trackX = state's own (pre-rotation) X,Y rotated by the reference's delta.
  const float trackX = scratchState.referenceCoordinate * ca + scratchState.parameters[0] * sa;

  if (!propagateReferenceParams(scratchRef, state.absCharge, trackX, bz, reason)) {
    reason = OperationFailureReason::RotationFailure;
    return false;
  }

  // Rotate state itself, gated by state's own snp (already checked above)
  // and its own post-rotation validity.
  const float csp = std::sqrt((1.f - stateSnp) * (1.f + stateSnp));
  if (csp * ca + stateSnp * sa < 0.f) {
    reason = OperationFailureReason::RotationFailure;
    return false;
  }
  const float updatedSnp = stateSnp * ca - csp * sa;
  if (std::abs(updatedSnp) >= 1.f) {
    reason = OperationFailureReason::RotationFailure;
    return false;
  }
  const float stateXOld = scratchState.referenceCoordinate;
  const float stateYOld = scratchState.parameters[0];
  scratchState.parameters[0] = -stateXOld * sa + stateYOld * ca;
  scratchState.referenceCoordinate = trackX;
  scratchState.parameters[2] = updatedSnp;
  scratchState.alpha = canonicalAlpha;

  // Covariance: Jacobian evaluated at the reference (cspRef0/cspRef1), not
  // at state's own snp -- the defining difference from the non-linRef
  // rotate() above. cspRef1 is computed algebraically (ca*cspRef0 +
  // sa*snpRef0), matching the legacy formula exactly rather than
  // re-deriving it via sqrt(1-refSnpRotated^2).
  const float cspRef1 = ca * refCsp0 + sa * refSnpBefore;
  if (cspRef1 == 0.f) {
    reason = OperationFailureReason::RotationFailure;
    return false;
  }
  const float rr = cspRef1 / refCsp0;

  // "Extra row" of the lower triangle, computed from the covariance values
  // as they stand *before* the plane-rotation multiplies below (matching
  // the legacy evaluation order exactly).
  const float cXSigY = scratchState.covariance[packedCovarianceIndex(0, 0)] * ca * sa;
  const float cXSigZ = scratchState.covariance[packedCovarianceIndex(1, 0)] * sa;
  const float cXSigSnp = scratchState.covariance[packedCovarianceIndex(2, 0)] * rr * sa;
  const float cXSigTgl = scratchState.covariance[packedCovarianceIndex(3, 0)] * sa;
  const float cXSigQ2Pt = scratchState.covariance[packedCovarianceIndex(4, 0)] * sa;
  const float cSigX2 = scratchState.covariance[packedCovarianceIndex(0, 0)] * sa * sa;

  scratchState.covariance[packedCovarianceIndex(0, 0)] *= ca * ca;
  scratchState.covariance[packedCovarianceIndex(1, 0)] *= ca;
  scratchState.covariance[packedCovarianceIndex(2, 0)] *= ca * rr;
  scratchState.covariance[packedCovarianceIndex(2, 1)] *= rr;
  scratchState.covariance[packedCovarianceIndex(2, 2)] *= rr * rr;
  scratchState.covariance[packedCovarianceIndex(3, 0)] *= ca;
  scratchState.covariance[packedCovarianceIndex(3, 2)] *= rr;
  scratchState.covariance[packedCovarianceIndex(4, 0)] *= ca;
  scratchState.covariance[packedCovarianceIndex(4, 2)] *= rr;

  const float cspRef1Inv = 1.f / cspRef1;
  const float j3 = -refSnpRotated * cspRef1Inv;
  const float j4 = -scratchRef.parameters[3] * cspRef1Inv;
  const float j5 = state.absCharge != 0 ? scratchRef.parameters[4] * bz * o2::constants::math::B2C : 0.f;

  const float hXSigY = cXSigY + cSigX2 * j3;
  const float hXSigZ = cXSigZ + cSigX2 * j4;
  const float hXSigSnp = cXSigSnp + cSigX2 * j5;

  scratchState.covariance[packedCovarianceIndex(0, 0)] += j3 * (cXSigY + hXSigY);
  scratchState.covariance[packedCovarianceIndex(1, 1)] += j4 * (cXSigZ + hXSigZ);
  scratchState.covariance[packedCovarianceIndex(2, 0)] += cXSigSnp * j3 + hXSigY * j5;
  scratchState.covariance[packedCovarianceIndex(2, 2)] += j5 * (cXSigSnp + hXSigSnp);
  scratchState.covariance[packedCovarianceIndex(3, 1)] += cXSigTgl * j4;
  scratchState.covariance[packedCovarianceIndex(4, 0)] += cXSigQ2Pt * j3;
  scratchState.covariance[packedCovarianceIndex(4, 2)] += cXSigQ2Pt * j5;

  scratchState.covariance[packedCovarianceIndex(1, 0)] += cXSigZ * j3 + hXSigY * j4;
  scratchState.covariance[packedCovarianceIndex(2, 1)] += cXSigSnp * j4 + hXSigZ * j5;
  scratchState.covariance[packedCovarianceIndex(3, 0)] += cXSigTgl * j3;
  scratchState.covariance[packedCovarianceIndex(3, 2)] += cXSigTgl * j5;
  scratchState.covariance[packedCovarianceIndex(4, 1)] += cXSigQ2Pt * j4;

  if (!finiteState(scratchState) || !finiteLinRef(scratchRef)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }
  sanitizeCovariance(scratchState, kBarrelMaxDiagonal);
  state = scratchState;
  linRef = scratchRef;
  return true;
}

bool propagate(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef, float targetX, float bz,
               OperationFailureReason& reason) noexcept
{
  if (!validateSource(state, reason)) {
    return false;
  }
  if (linRef.family != StateFamily::Barrel) {
    reason = OperationFailureReason::SourceFamilyMismatch;
    return false;
  }
  if (!finiteLinRef(linRef) || !std::isfinite(targetX) || !std::isfinite(bz)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }
  // Fitted-state/linRef pairing precondition -- see the identical check in
  // rotate() above for the rationale (exact comparison; parameters may
  // differ, anchor/frame may not).
  if (state.referenceCoordinate != linRef.referenceCoordinate) {
    reason = OperationFailureReason::ReferenceCoordinateMismatch;
    return false;
  }
  if (state.alpha != linRef.alpha) {
    reason = OperationFailureReason::AlphaMismatch;
    return false;
  }

  const float effectiveBz = state.absCharge == 0 ? 0.f : bz;
  const float dx = targetX - state.referenceCoordinate;
  if (std::abs(dx) < o2::constants::math::Almost0) {
    SurfaceKinematicState scratchState = state;
    SurfaceLinearizationReference scratchRef = linRef;
    scratchState.referenceCoordinate = targetX;
    scratchRef.referenceCoordinate = targetX;
    state = scratchState;
    linRef = scratchRef;
    return true;
  }

  SurfaceLinearizationReference scratchRef = linRef;
  const float snpRef0 = scratchRef.parameters[2];
  const float cspRef0 = std::sqrt((1.f - snpRef0) * (1.f + snpRef0));
  const float tglRef0 = scratchRef.parameters[3];

  if (!propagateReferenceParams(scratchRef, state.absCharge, targetX, effectiveBz, reason)) {
    return false;
  }
  const float snpRef1 = scratchRef.parameters[2];
  const float cspRef1 = std::sqrt((1.f - snpRef1) * (1.f + snpRef1));
  if (cspRef0 == 0.f || cspRef1 == 0.f) {
    reason = OperationFailureReason::PropagationFailure;
    return false;
  }

  const float kb = effectiveBz * o2::constants::math::B2C;
  const float cspRef0Inv = 1.f / cspRef0;
  const float cspRef1Inv = 1.f / cspRef1;
  const float cc = cspRef0 + cspRef1;
  const float ccInv = 1.f / cc;
  const float dy2dx = (snpRef0 + snpRef1) * ccInv;
  const float dxccInv = dx * ccInv;
  const float hh = dxccInv * cspRef1Inv * (1.f + cspRef0 * cspRef1 + snpRef0 * snpRef1);
  const float jj = dx * (dy2dx - snpRef1 * cspRef1Inv);

  const float f02 = hh * cspRef0Inv;
  const float f04 = hh * dxccInv * kb;
  const float f24 = dx * kb;
  const float f12 = tglRef0 * (f02 * snpRef1 + jj);
  const float f13 = dx * (cspRef1 + snpRef1 * dy2dx);
  const float f14 = tglRef0 * (f04 * snpRef1 + jj * f24);

  float diff[5];
  for (uint8_t i = 0; i < 5; ++i) {
    diff[i] = state.parameters[i] - linRef.parameters[i];
  }
  const float snpUpd = snpRef1 + diff[2] + f24 * diff[4];
  if (std::abs(snpUpd) >= 1.f) {
    reason = OperationFailureReason::PropagationFailure;
    return false;
  }

  SurfaceKinematicState scratchState = state;
  scratchState.referenceCoordinate = targetX;
  scratchState.parameters[0] = scratchRef.parameters[0] + diff[0] + f02 * diff[2] + f04 * diff[4];
  scratchState.parameters[1] = scratchRef.parameters[1] + diff[1] + f13 * diff[3] + f14 * diff[4];
  scratchState.parameters[2] = snpUpd;
  scratchState.parameters[3] = scratchRef.parameters[3] + diff[3];
  scratchState.parameters[4] = scratchRef.parameters[4] + diff[4];

  const float c00 = state.covariance[packedCovarianceIndex(0, 0)];
  const float c10 = state.covariance[packedCovarianceIndex(1, 0)];
  const float c11 = state.covariance[packedCovarianceIndex(1, 1)];
  const float c20 = state.covariance[packedCovarianceIndex(2, 0)];
  const float c21 = state.covariance[packedCovarianceIndex(2, 1)];
  const float c22 = state.covariance[packedCovarianceIndex(2, 2)];
  const float c30 = state.covariance[packedCovarianceIndex(3, 0)];
  const float c31 = state.covariance[packedCovarianceIndex(3, 1)];
  const float c32 = state.covariance[packedCovarianceIndex(3, 2)];
  const float c33 = state.covariance[packedCovarianceIndex(3, 3)];
  const float c40 = state.covariance[packedCovarianceIndex(4, 0)];
  const float c41 = state.covariance[packedCovarianceIndex(4, 1)];
  const float c42 = state.covariance[packedCovarianceIndex(4, 2)];
  const float c43 = state.covariance[packedCovarianceIndex(4, 3)];
  const float c44 = state.covariance[packedCovarianceIndex(4, 4)];

  const float b00 = f02 * c20 + f04 * c40;
  const float b01 = f12 * c20 + f14 * c40 + f13 * c30;
  const float b02 = f24 * c40;
  const float b10 = f02 * c21 + f04 * c41;
  const float b11 = f12 * c21 + f14 * c41 + f13 * c31;
  const float b12 = f24 * c41;
  const float b20 = f02 * c22 + f04 * c42;
  const float b21 = f12 * c22 + f14 * c42 + f13 * c32;
  const float b22 = f24 * c42;
  const float b40 = f02 * c42 + f04 * c44;
  const float b41 = f12 * c42 + f14 * c44 + f13 * c43;
  const float b42 = f24 * c44;
  const float b30 = f02 * c32 + f04 * c43;
  const float b31 = f12 * c32 + f14 * c43 + f13 * c33;
  const float b32 = f24 * c43;

  const float a00 = f02 * b20 + f04 * b40;
  const float a01 = f02 * b21 + f04 * b41;
  const float a02 = f02 * b22 + f04 * b42;
  const float a11 = f12 * b21 + f14 * b41 + f13 * b31;
  const float a12 = f12 * b22 + f14 * b42 + f13 * b32;
  const float a22 = f24 * b42;

  scratchState.covariance[packedCovarianceIndex(0, 0)] = c00 + b00 + b00 + a00;
  scratchState.covariance[packedCovarianceIndex(1, 0)] = c10 + b10 + b01 + a01;
  scratchState.covariance[packedCovarianceIndex(2, 0)] = c20 + b20 + b02 + a02;
  scratchState.covariance[packedCovarianceIndex(3, 0)] = c30 + b30;
  scratchState.covariance[packedCovarianceIndex(4, 0)] = c40 + b40;
  scratchState.covariance[packedCovarianceIndex(1, 1)] = c11 + b11 + b11 + a11;
  scratchState.covariance[packedCovarianceIndex(2, 1)] = c21 + b21 + b12 + a12;
  scratchState.covariance[packedCovarianceIndex(3, 1)] = c31 + b31;
  scratchState.covariance[packedCovarianceIndex(4, 1)] = c41 + b41;
  scratchState.covariance[packedCovarianceIndex(2, 2)] = c22 + b22 + b22 + a22;
  scratchState.covariance[packedCovarianceIndex(3, 2)] = c32 + b32;
  scratchState.covariance[packedCovarianceIndex(4, 2)] = c42 + b42;
  scratchState.covariance[packedCovarianceIndex(3, 3)] = c33;
  scratchState.covariance[packedCovarianceIndex(4, 3)] = c43;
  scratchState.covariance[packedCovarianceIndex(4, 4)] = c44;

  if (!finiteState(scratchState) || !finiteLinRef(scratchRef)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }
  // ADR 0008: a large single-step Jacobian transport can produce an
  // off-diagonal term large enough that the matrix is no longer
  // positive-semi-definite even though every individual diagonal still
  // looks valid; sanitize unconditionally so the next operation (typically
  // a measurement update) never receives an already-invalid covariance.
  sanitizeCovariance(scratchState, kBarrelMaxDiagonal);
  state = scratchState;
  linRef = scratchRef;
  return true;
}

bool shiftReferenceToMeasurement(SurfaceLinearizationReference& linRef, const SurfaceMeasurement& measurement,
                                 OperationFailureReason& reason) noexcept
{
  if (linRef.family != StateFamily::Barrel) {
    reason = OperationFailureReason::SourceFamilyMismatch;
    return false;
  }
  if (!std::isfinite(measurement.frame.u) || !std::isfinite(measurement.frame.v)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }
  SurfaceLinearizationReference scratch = linRef;
  scratch.parameters[0] = measurement.frame.u;
  scratch.parameters[1] = measurement.frame.v;
  if (!finiteLinRef(scratch)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }
  linRef = scratch;
  return true;
}

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking::barrel
