// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file Propagator.cxx
/// \brief M5d generic, descriptor-driven SurfaceKinematicState propagator
///

#include "ITSMFTTracking/Propagator.h"

#include <cmath>

#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

namespace o2::itsmft::tracking
{

namespace
{

// Minimum accepted-hit count (within one Propagator::driveRefitLeg leg) at
// which the maxChi2 gate activates -- the same threshold and rationale as
// detail::kRefitLegChi2GateMinAcceptedHits (TransitionPolicyOperations.h),
// duplicated here (rather than included) so this generic Propagator has no
// dependency on that Tag-templated legacy-adjacent detail header at all.
constexpr uint32_t kChi2GateMinAcceptedHits = 3;

// Reuses the exact |bz| > 0.01f Helix/Linear threshold already established
// as "the accepted forward model" by detail::forwardPropagateAcceptedModel
// (TransitionPolicyOperations.cxx, anonymous namespace -- not exported).
// Duplicated here, narrowly, rather than exported from that legacy-adjacent
// translation unit: this is exactly the kind of small, family-local leaf
// choice ADR 0007 decision 10 expects to live at the operation boundary, not
// shared orchestration.
bool forwardPropagateAcceptedModel(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef,
                                   float targetZ, float bz, OperationFailureReason& reason) noexcept
{
  if (std::abs(bz) > 0.01f) {
    return forward::propagate<forward::PropagationModel::Helix>(state, linRef, targetZ, bz, reason);
  }
  return forward::propagate<forward::PropagationModel::Linear>(state, linRef, targetZ, bz, reason);
}

// A congruence transform (rotate/propagate's own Jacobian-based covariance
// transport) that starts from a "loose ceiling" reset diagonal (see
// NativeRefitDriver.h's resetCovarianceForRefit) and crosses a large step at
// an extreme dip angle can leave a diagonal entry that is mathematically
// exactly zero as an infinitesimally negative float (catastrophic
// cancellation inside the Jacobian congruence, not a sign of a genuinely
// corrupted covariance) -- observed for a near-longitudinal forward
// trajectory. correctForMaterialImpl's own covarianceDiagonalsNonNegative()
// check is an exact `< 0.f` comparison with no tolerance, so this noise must
// be floored before it reaches that check. Deliberately narrow: only
// sub-float-precision noise (magnitude below this threshold, itself already
// tiny next to the smallest loose-ceiling constant these states are ever
// reset to, o2::track::kCSnp2max == 1) is floored to zero; anything larger
// still fails downstream exactly as before, so a genuine covariance-transport
// defect is not masked.
void clampNegligibleCovarianceNoise(SurfaceKinematicState& state) noexcept
{
  constexpr float kNoiseFloor = 1.e-3f;
  for (uint8_t i = 0; i < 5; ++i) {
    const uint8_t index = packedCovarianceIndex(i, i);
    if (state.covariance[index] < 0.f && state.covariance[index] > -kNoiseFloor) {
      state.covariance[index] = 0.f;
    }
  }
}

bool isFinite5(const float (&values)[5]) noexcept
{
  for (float value : values) {
    if (!std::isfinite(value)) {
      return false;
    }
  }
  return true;
}

bool isFiniteCov(const float (&cov)[15]) noexcept
{
  for (float value : cov) {
    if (!std::isfinite(value)) {
      return false;
    }
  }
  return true;
}

// Congruence transform outCov = J * inCov * J^T for packed-symmetric 5x5
// covariances, J a dense 5x5 Jacobian. Shared by both conversion directions
// below.
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

// Barrel (bY, bZ, Snp, Tgl, Q2Pt) @ (referenceCoordinate=bX, alpha) ->
// Forward (X, Y, Phi, Tanl, InvQPt) @ (referenceCoordinate=Z). A real
// coordinate/covariance transform at the state's current point: bX/alpha
// select the global point (Xg, Yg, Zg) the state already sits at, and that
// same point becomes Forward's (X, Y) parameters / Z reference. bZ's own
// fitted variance is discarded as it becomes the new fixed reference
// (Propagator.h's class doc explains why); every other quantity is a real,
// linearized (first-order Jacobian) reparametrization, not a relabeling.
bool barrelToForward(SurfaceKinematicState& state, SurfaceLinearizationReference* linRef,
                     OperationFailureReason& reason) noexcept
{
  if (!isFinite5(state.parameters) || !isFiniteCov(state.covariance) ||
      !std::isfinite(state.referenceCoordinate) || !std::isfinite(state.alpha)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }
  const float snp = state.parameters[2];
  if (!(std::abs(snp) < 1.f)) {
    reason = OperationFailureReason::FamilyConversionFailure;
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
  float phi = state.alpha + std::asin(snp);
  // Canonicalize to (-pi, pi], matching this library's own alpha convention
  // elsewhere (barrel::rotate).
  while (phi > o2::constants::math::PI) {
    phi -= o2::constants::math::TwoPI;
  }
  while (phi <= -o2::constants::math::PI) {
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
  if (!isFinite5(newParameters) || !isFiniteCov(newCov) || !std::isfinite(zGlo)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }

  if (linRef != nullptr) {
    if (linRef->family != state.family || linRef->referenceCoordinate != state.referenceCoordinate || linRef->alpha != state.alpha) {
      reason = OperationFailureReason::FamilyConversionFailure;
      return false;
    }
    const float linSnp = linRef->parameters[2];
    if (!(std::abs(linSnp) < 1.f)) {
      reason = OperationFailureReason::FamilyConversionFailure;
      return false;
    }
    const float linBX = linRef->referenceCoordinate;
    const float linBY = linRef->parameters[0];
    const float linXGlo = linBX * csA - linBY * snA;
    const float linYGlo = linBX * snA + linBY * csA;
    const float linZGlo = linRef->parameters[1];
    float linPhi = linRef->alpha + std::asin(linSnp);
    while (linPhi > o2::constants::math::PI) {
      linPhi -= o2::constants::math::TwoPI;
    }
    while (linPhi <= -o2::constants::math::PI) {
      linPhi += o2::constants::math::TwoPI;
    }
    const float newLinParameters[5] = {linXGlo, linYGlo, linPhi, linRef->parameters[3], linRef->parameters[4]};
    if (!isFinite5(newLinParameters) || !std::isfinite(linZGlo)) {
      reason = OperationFailureReason::NonFiniteOutput;
      return false;
    }
    for (uint8_t i = 0; i < 5; ++i) {
      linRef->parameters[i] = newLinParameters[i];
    }
    linRef->referenceCoordinate = linZGlo;
    linRef->alpha = 0.f;
    linRef->family = StateFamily::Forward;
  }

  for (uint8_t i = 0; i < 5; ++i) {
    state.parameters[i] = newParameters[i];
  }
  for (uint8_t i = 0; i < 15; ++i) {
    state.covariance[i] = newCov[i];
  }
  state.referenceCoordinate = zGlo;
  state.alpha = 0.f;
  state.family = StateFamily::Forward;
  return true;
}

// Forward (X, Y, Phi, Tanl, InvQPt) @ (referenceCoordinate=Z) -> Barrel (bY,
// bZ, Snp, Tgl, Q2Pt) @ (referenceCoordinate=bX, alpha). alpha is chosen
// once as atan2(Y, X) at the current nominal point and then frozen for the
// linearized Jacobian (the same "linearize at an externally supplied,
// fixed target" technique barrel::rotate/propagate already use for their
// own targetAlpha argument) -- not a self-referential function of the
// varying parameters. Z's own fitted variance has no Forward-side
// representation to transport (Forward's referenceCoordinate carries no
// variance slot), so the newly freed bZ parameter is assigned the same
// loose, uninformative ceiling this library already uses for an
// analogous "we do not actually know this coordinate's uncertainty here"
// situation (o2::track::kCZ2max).
bool forwardToBarrel(SurfaceKinematicState& state, SurfaceLinearizationReference* linRef,
                     OperationFailureReason& reason) noexcept
{
  if (!isFinite5(state.parameters) || !isFiniteCov(state.covariance) || !std::isfinite(state.referenceCoordinate)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }
  const float x = state.parameters[0];
  const float y = state.parameters[1];
  const float r = std::sqrt(x * x + y * y);
  if (!(r > 1.e-6f)) {
    reason = OperationFailureReason::FamilyConversionFailure;
    return false;
  }
  const float alpha = std::atan2(y, x);
  const float csA = std::cos(alpha);
  const float snA = std::sin(alpha);
  const float phi = state.parameters[2];
  const float csp = std::cos(phi - alpha);
  const float snp = std::sin(phi - alpha);
  if (!(std::abs(snp) < 1.f)) {
    reason = OperationFailureReason::FamilyConversionFailure;
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
  if (!isFinite5(newParameters) || !isFiniteCov(newCov) || !std::isfinite(bX) || !std::isfinite(alpha)) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }

  if (linRef != nullptr) {
    if (linRef->family != state.family || linRef->referenceCoordinate != state.referenceCoordinate) {
      reason = OperationFailureReason::FamilyConversionFailure;
      return false;
    }
    const float linX = linRef->parameters[0];
    const float linY = linRef->parameters[1];
    const float linPhi = linRef->parameters[2];
    const float linBY = -linX * snA + linY * csA;
    const float linBX = linX * csA + linY * snA;
    const float linSnp = std::sin(linPhi - alpha);
    const float newLinParameters[5] = {linBY, linRef->referenceCoordinate, linSnp, linRef->parameters[3], linRef->parameters[4]};
    if (!isFinite5(newLinParameters) || !std::isfinite(linBX)) {
      reason = OperationFailureReason::NonFiniteOutput;
      return false;
    }
    for (uint8_t i = 0; i < 5; ++i) {
      linRef->parameters[i] = newLinParameters[i];
    }
    linRef->referenceCoordinate = linBX;
    linRef->alpha = alpha;
    linRef->family = StateFamily::Barrel;
  }

  for (uint8_t i = 0; i < 5; ++i) {
    state.parameters[i] = newParameters[i];
  }
  for (uint8_t i = 0; i < 15; ++i) {
    state.covariance[i] = newCov[i];
  }
  state.referenceCoordinate = bX;
  state.alpha = alpha;
  state.family = StateFamily::Barrel;
  return true;
}

} // namespace

bool Propagator::convertFamily(SurfaceKinematicState& state, SurfaceLinearizationReference* linRef,
                               StateFamily targetFamily, OperationFailureReason& reason) noexcept
{
  if (targetFamily != StateFamily::Barrel && targetFamily != StateFamily::Forward) {
    reason = OperationFailureReason::FamilyConversionFailure;
    return false;
  }
  if (state.family != StateFamily::Barrel && state.family != StateFamily::Forward) {
    reason = OperationFailureReason::SourceFamilyMismatch;
    return false;
  }
  if (state.family == targetFamily) {
    return true;
  }
  if (targetFamily == StateFamily::Forward) {
    return barrelToForward(state, linRef, reason);
  }
  return forwardToBarrel(state, linRef, reason);
}

bool Propagator::propagateToMeasurement(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef,
                                        const SurfaceDescriptor& targetSurface, const SurfaceMeasurement& targetMeasurement,
                                        float bz, material::MaterialTraversalDirection direction,
                                        bool chi2GateEnabled, float maxChi2, float& chi2,
                                        bool shiftReferenceToMeasurement, OperationFailureReason& reason) noexcept
{
  if (!std::isfinite(chi2)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }
  if (chi2 < 0.f) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  if (chi2GateEnabled) {
    if (!std::isfinite(maxChi2)) {
      reason = OperationFailureReason::NonFiniteInput;
      return false;
    }
    if (maxChi2 < 0.f) {
      reason = OperationFailureReason::PredictedChi2Failure;
      return false;
    }
  }

  const StateFamily targetFamily = stateFamilyOf(targetSurface.kind);
  if (targetFamily == StateFamily::Invalid) {
    reason = OperationFailureReason::FamilyConversionFailure;
    return false;
  }

  SurfaceKinematicState scratchState = state;
  SurfaceLinearizationReference scratchRef = linRef;

  if (scratchState.family != targetFamily) {
    if (!convertFamily(scratchState, &scratchRef, targetFamily, reason)) {
      return false;
    }
  }

  const material::IntegratedMaterialBudget materialBudget{targetSurface.material.xOverX0, targetSurface.material.arealDensityGPerCm2};
  float scratchChi2 = chi2;
  float predChi2 = 0.f;
  float updateChi2 = 0.f;

  if (targetFamily == StateFamily::Barrel) {
    if (!barrel::rotate(scratchState, scratchRef, targetMeasurement.frame.frameAngle, bz, reason)) {
      return false;
    }
    if (!barrel::propagate(scratchState, scratchRef, targetMeasurement.frame.q, bz, reason)) {
      return false;
    }
    clampNegligibleCovarianceNoise(scratchState);
    const auto materialResult = barrel::correctForMaterial(scratchState, materialBudget, direction);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    if (!barrel::predictedChi2(scratchState, targetMeasurement, predChi2, reason)) {
      return false;
    }
  } else {
    if (!forwardPropagateAcceptedModel(scratchState, scratchRef, targetMeasurement.frame.q, bz, reason)) {
      return false;
    }
    clampNegligibleCovarianceNoise(scratchState);
    const auto materialResult = forward::correctForMaterial(scratchState, materialBudget, direction);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    if (!forward::predictedChi2(scratchState, targetMeasurement, predChi2, reason)) {
      return false;
    }
  }

  if (predChi2 < 0.f || !std::isfinite(predChi2)) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  if (chi2GateEnabled && predChi2 > maxChi2) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }

  if (targetFamily == StateFamily::Barrel) {
    if (!barrel::update(scratchState, targetMeasurement, updateChi2, reason)) {
      return false;
    }
  } else {
    if (!forward::update(scratchState, targetMeasurement, updateChi2, reason)) {
      return false;
    }
  }
  scratchChi2 += updateChi2;
  if (!std::isfinite(scratchChi2) || scratchChi2 < 0.f) {
    reason = OperationFailureReason::NonFiniteOutput;
    return false;
  }

  if (shiftReferenceToMeasurement) {
    if (targetFamily == StateFamily::Barrel) {
      if (!barrel::shiftReferenceToMeasurement(scratchRef, targetMeasurement, reason)) {
        return false;
      }
    } else {
      if (!forward::shiftReferenceToMeasurement(scratchRef, targetMeasurement, reason)) {
        return false;
      }
    }
  }

  state = scratchState;
  linRef = scratchRef;
  chi2 = scratchChi2;
  return true;
}

bool Propagator::driveRefitLeg(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef,
                               float& chi2, uint32_t& acceptedHitCount,
                               gsl::span<const SurfaceMeasurement> orderedSlots, SurfaceCatalogView surfaceCatalog,
                               float bz, material::MaterialTraversalDirection direction,
                               bool shiftReferenceToMeasurement, float maxChi2, OperationFailureReason& reason) noexcept
{
  if (!std::isfinite(chi2)) {
    reason = OperationFailureReason::NonFiniteInput;
    return false;
  }
  if (chi2 < 0.f) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }

  SurfaceKinematicState scratchState = state;
  SurfaceLinearizationReference scratchLinRef = linRef;
  float scratchChi2 = chi2;
  uint32_t scratchAcceptedHitCount = 0;

  for (const auto& measurement : orderedSlots) {
    if (!measurement.cluster.isValid()) {
      continue;
    }

    if (!measurement.surface.isValid()) {
      reason = OperationFailureReason::InvalidSurfaceCatalogAssociation;
      return false;
    }
    if (!(surfaceCatalog.nSurfaces == 0 || surfaceCatalog.surfaces != nullptr)) {
      reason = OperationFailureReason::InvalidSurfaceCatalogAssociation;
      return false;
    }
    if (!(measurement.surface.value() < surfaceCatalog.nSurfaces)) {
      reason = OperationFailureReason::InvalidSurfaceCatalogAssociation;
      return false;
    }
    const SurfaceDescriptor& descriptor = surfaceCatalog.getSurface(measurement.surface);
    if (!(descriptor.id == measurement.surface)) {
      reason = OperationFailureReason::InvalidSurfaceCatalogAssociation;
      return false;
    }

    const bool chi2GateEnabled = scratchAcceptedHitCount >= kChi2GateMinAcceptedHits;
    if (!propagateToMeasurement(scratchState, scratchLinRef, descriptor, measurement, bz, direction,
                                chi2GateEnabled, maxChi2, scratchChi2, shiftReferenceToMeasurement, reason)) {
      return false;
    }
    ++scratchAcceptedHitCount;
  }

  state = scratchState;
  linRef = scratchLinRef;
  chi2 = scratchChi2;
  acceptedHitCount = scratchAcceptedHitCount;
  return true;
}

} // namespace o2::itsmft::tracking
