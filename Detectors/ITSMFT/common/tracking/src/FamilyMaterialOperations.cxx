// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Defines both barrel::correctForMaterial(state, material, direction) and
// forward::correctForMaterial(state, material, direction): the PID/absCharge-
// aware composite family operations built on top of the detector-neutral
// scalar kernel in MaterialPhysics.h. Both overloads share the complete
// preflight-validation/momentum-derivation/scratch-and-commit orchestration
// below; only the family-specific extra kinematics check and covariance
// projection formula differ between them.
//
// This translation unit is host-only, does not construct or delegate through
// TrackParCovF or TrackParCovFwd, and includes TrackParametrization.h solely
// to reuse its public kCY2max/kCZ2max/kCSnp2max/kCTgl2max/kC1Pt2max constants
// for the retained barrel covariance-range handling (no narrower public
// header declares them; the same reuse pattern is already used by
// MaterialPhysics.cxx for its own constants).

#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/MaterialPhysics.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

namespace o2::itsmft::tracking
{
namespace
{

bool allStateFloatsFinite(const SurfaceKinematicState& state) noexcept
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

bool covarianceDiagonalsNonNegative(const SurfaceKinematicState& state) noexcept
{
  for (uint8_t i = 0; i < 5; ++i) {
    if (state.covariance[packedCovarianceIndex(i, i)] < 0.f) {
      return false;
    }
  }
  return true;
}

// Physical-momentum derivation shared by both families. u = slot 4, t = slot
// 3. Rejects overflow, underflow to zero, and non-finite derived values.
bool derivePhysicalMomentum(const SurfaceKinematicState& state, float& momentumGeV) noexcept
{
  const float t = state.parameters[3];
  const float u = state.parameters[4];
  const float absU = std::abs(u);
  const float pT = (state.absCharge == 0) ? (1.f / absU) : (static_cast<float>(state.absCharge) / absU);
  if (!std::isfinite(pT) || pT <= 0.f) {
    return false;
  }
  const float p = pT * std::sqrt(1.f + t * t);
  if (!std::isfinite(p) || p <= 0.f) {
    return false;
  }
  momentumGeV = p;
  return true;
}

material::MaterialOperationResult makePreflightFailure(material::MaterialFailureReason reason) noexcept
{
  material::MaterialOperationResult result{};
  result.momentumBeforeGeV = 0.f;
  result.momentumAfterGeV = 0.f;
  result.signedEnergyChangeGeV = 0.f;
  result.highlandTheta2Rad2 = 0.f;
  result.relativeInverseMomentumVariance = 0.f;
  result.energyLossSubsteps = 0;
  result.flags = material::MaterialOperationFlags::None;
  result.failure = reason;
  result.reserved = 0;
  return result;
}

material::MaterialOperationResult makeProjectionFailure(const material::MaterialOperationResult& scalarResult,
                                                        material::MaterialFailureReason reason) noexcept
{
  material::MaterialOperationResult result{};
  result.momentumBeforeGeV = scalarResult.momentumBeforeGeV;
  result.momentumAfterGeV = 0.f;
  result.signedEnergyChangeGeV = 0.f;
  result.highlandTheta2Rad2 = 0.f;
  result.relativeInverseMomentumVariance = 0.f;
  result.energyLossSubsteps = 0;
  result.flags = material::MaterialOperationFlags::None;
  result.failure = reason;
  result.reserved = 0;
  return result;
}

// Translation-unit-private packed-array reproduction of the retained
// TrackParametrizationWithError<float>::checkCovariance() range handling: for
// each diagonal, take its absolute value, and if it exceeds the retained
// maximum, clamp it and scale every off-diagonal entry involving that
// parameter by sqrt(max/diagonal). Diagonals are processed in the retained
// (Y, Z, Snp, Tgl, Q2Pt) order so cumulative scaling of shared off-diagonal
// entries (e.g. the Tgl/Q2Pt cross term) matches the retained sequential
// behavior bit-for-bit. No legacy track object is constructed; this operates
// directly on the packed float covariance array. This is a faithful
// reproduction with no intentional difference from the retained
// implementation.
void limitBarrelCovariance(SurfaceKinematicState& scratch) noexcept
{
  constexpr float kMaxValues[5] = {o2::track::kCY2max, o2::track::kCZ2max, o2::track::kCSnp2max,
                                   o2::track::kCTgl2max, o2::track::kC1Pt2max};
  auto& c = scratch.covariance;
  for (uint8_t i = 0; i < 5; ++i) {
    const uint8_t diagIndex = packedCovarianceIndex(i, i);
    c[diagIndex] = std::abs(c[diagIndex]);
    if (c[diagIndex] > kMaxValues[i]) {
      const float scale = std::sqrt(kMaxValues[i] / c[diagIndex]);
      c[diagIndex] = kMaxValues[i];
      for (uint8_t j = 0; j < 5; ++j) {
        if (j == i) {
          continue;
        }
        c[packedCovarianceIndex(std::max(i, j), std::min(i, j))] *= scale;
      }
    }
  }
}

// Shared preflight validation, steps 1-6 of the required order. Step 3's
// family-specific extra check (barrel |Snp|<1 / forward alpha==0) is
// supplied by the caller; the shared slope-finite/slot4-finite-nonzero part
// of step 3 is applied here for both families.
template <typename FamilyKinematicsCheck>
bool preflightValidate(const SurfaceKinematicState& state, StateFamily expectedFamily, FamilyKinematicsCheck&& familyCheck,
                       material::MaterialFailureReason& failure) noexcept
{
  if (state.family != expectedFamily) {
    failure = material::MaterialFailureReason::SourceFamilyMismatch;
    return false;
  }
  if (!allStateFloatsFinite(state)) {
    failure = material::MaterialFailureReason::NonFiniteState;
    return false;
  }
  const float t = state.parameters[3];
  const float u = state.parameters[4];
  if (!familyCheck(state) || !std::isfinite(t) || !std::isfinite(u) || u == 0.f) {
    failure = material::MaterialFailureReason::InvalidStateKinematics;
    return false;
  }
  if (state.pid.getID() >= o2::track::PID::NIDsTot) {
    failure = material::MaterialFailureReason::InvalidPID;
    return false;
  }
  if (state.absCharge != 0 && state.pid.getMass() == 0.f) {
    failure = material::MaterialFailureReason::ChargedMasslessPID;
    return false;
  }
  if (!covarianceDiagonalsNonNegative(state)) {
    failure = material::MaterialFailureReason::InvalidCovariance;
    return false;
  }
  return true;
}

// Complete transactional operation shared by both families (Slice 2
// "Transactional result contract"): validate, derive physical momentum,
// scratch-copy, invoke the scalar kernel, project covariance on scratch
// only, validate the projected scratch, and commit exactly once.
// projectCovariance may additionally apply family-specific covariance range
// handling (barrel only); it must not touch state.parameters[4], which this
// function updates uniformly for both families after projection.
template <typename FamilyKinematicsCheck, typename ProjectCovariance>
material::MaterialOperationResult correctForMaterialImpl(SurfaceKinematicState& state, StateFamily expectedFamily,
                                                         material::IntegratedMaterialBudget materialBudget,
                                                         material::MaterialTraversalDirection direction,
                                                         FamilyKinematicsCheck&& familyCheck,
                                                         ProjectCovariance&& projectCovariance) noexcept
{
  material::MaterialFailureReason failure{};
  if (!preflightValidate(state, expectedFamily, familyCheck, failure)) {
    return makePreflightFailure(failure);
  }

  float momentumBeforeGeV = 0.f;
  if (!derivePhysicalMomentum(state, momentumBeforeGeV)) {
    return makePreflightFailure(material::MaterialFailureReason::InvalidStateKinematics);
  }

  SurfaceKinematicState scratch = state;
  const auto scalarResult = material::calculateMaterialPhysics(momentumBeforeGeV, scratch.pid, scratch.absCharge, direction, materialBudget);
  if (!scalarResult.ok()) {
    return scalarResult;
  }

  const float tBefore = scratch.parameters[3];
  const float kBefore = scratch.parameters[4];
  projectCovariance(scratch, scalarResult, tBefore, kBefore);

  // Compute the ratio before multiplying: when momentumBeforeGeV and
  // momentumAfterGeV are bit-identical (zero material, or neutral states),
  // x/x is exactly 1.0f for any finite nonzero x, so kAfter == kBefore
  // bit-for-bit. Multiplying first then dividing would round-trip through an
  // intermediate product and generally fail to recover kBefore exactly,
  // breaking the required byte-identical invariant.
  const float kAfter = kBefore * (scalarResult.momentumBeforeGeV / scalarResult.momentumAfterGeV);
  scratch.parameters[4] = kAfter;

  if (!allStateFloatsFinite(scratch)) {
    return makeProjectionFailure(scalarResult, material::MaterialFailureReason::NonFiniteResult);
  }
  if (!covarianceDiagonalsNonNegative(scratch)) {
    return makeProjectionFailure(scalarResult, material::MaterialFailureReason::InvalidCovariance);
  }

  state = scratch;
  return scalarResult;
}

} // namespace
} // namespace o2::itsmft::tracking

namespace o2::itsmft::tracking::barrel
{

material::MaterialOperationResult correctForMaterial(SurfaceKinematicState& state, material::IntegratedMaterialBudget materialBudget,
                                                     material::MaterialTraversalDirection direction) noexcept
{
  auto familyCheck = [](const SurfaceKinematicState& s) noexcept {
    return std::abs(s.parameters[2]) < 1.f;
  };
  // Barrel parameters are (Y, Z, Snp, Tgl, Q2Pt). The accepted Jacobian
  // requires q/pT unconditionally in slots 13/14, fixing the retained
  // TrackParametrizationWithError<float>::correctForMaterial() unit-charge
  // conditional that omits q/pT there (see the module doc comment).
  auto projectCovariance = [](SurfaceKinematicState& scratch, const material::MaterialOperationResult& scalarResult,
                              float t, float k) noexcept {
    const float A = 1.f + t * t;
    const float snp = scratch.parameters[2];
    const float c2 = 1.f - snp * snp;
    const float h = scalarResult.highlandTheta2Rad2;
    const float R = scalarResult.relativeInverseMomentumVariance;
    scratch.covariance[packedCovarianceIndex(2, 2)] += h * A * c2;
    scratch.covariance[packedCovarianceIndex(3, 3)] += h * A * A;
    scratch.covariance[packedCovarianceIndex(4, 3)] += h * A * t * k;
    scratch.covariance[packedCovarianceIndex(4, 4)] += h * (t * k) * (t * k) + k * k * R;
    limitBarrelCovariance(scratch);
  };
  return correctForMaterialImpl(state, StateFamily::Barrel, materialBudget, direction, familyCheck, projectCovariance);
}

} // namespace o2::itsmft::tracking::barrel

namespace o2::itsmft::tracking::forward
{

material::MaterialOperationResult correctForMaterial(SurfaceKinematicState& state, material::IntegratedMaterialBudget materialBudget,
                                                     material::MaterialTraversalDirection direction) noexcept
{
  auto familyCheck = [](const SurfaceKinematicState& s) noexcept {
    return s.alpha == 0.f;
  };
  // Forward parameters are (X, Y, Phi, Tanl, Q2Pt); unlike barrel there is no
  // cos(phi)-like factor on the angular diagonal term, and forward does not
  // inherit barrel-only covariance range limiting. Slot 13 (the Q2Pt/Tanl
  // cross term) and the k^2*R straggling contribution to slot 14 are new
  // physics: the legacy TrackParCovFwd::addMCSEffect() never populates
  // slot 13 and has no charge/PID/energy-loss awareness at all.
  auto projectCovariance = [](SurfaceKinematicState& scratch, const material::MaterialOperationResult& scalarResult,
                              float t, float k) noexcept {
    const float A = 1.f + t * t;
    const float h = scalarResult.highlandTheta2Rad2;
    const float R = scalarResult.relativeInverseMomentumVariance;
    scratch.covariance[packedCovarianceIndex(2, 2)] += h * A;
    scratch.covariance[packedCovarianceIndex(3, 3)] += h * A * A;
    scratch.covariance[packedCovarianceIndex(4, 3)] += h * A * t * k;
    scratch.covariance[packedCovarianceIndex(4, 4)] += h * (t * k) * (t * k) + k * k * R;
  };
  return correctForMaterialImpl(state, StateFamily::Forward, materialBudget, direction, familyCheck, projectCovariance);
}

} // namespace o2::itsmft::tracking::forward
