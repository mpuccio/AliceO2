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

// Defines both detail::barrel::correctForMaterial(state, material, direction) and
// detail::forward::correctForMaterial(state, material, direction): the PID/absCharge-
// aware composite cylinder/disk operations built on the detector-neutral
// scalar kernel in MaterialPhysics.h. Both overloads share the complete
// preflight-validation/momentum-derivation/scratch-and-commit orchestration
// below; only the coordinate-specific kinematics check and covariance
// projection formula differ between them.
//
// This translation unit is host-only, does not construct or delegate through
// TrackParCovF or TrackParCovFwd, and includes TrackParametrization.h solely
// to reuse its public kCY2max/kCZ2max/kCSnp2max/kCTgl2max/kC1Pt2max constants
// for the retained barrel covariance-range handling (no narrower public
// header declares them; the same reuse pattern is already used by
// MaterialPhysics.cxx for its own constants).

#include "ITSMFTTracking/detail/SurfaceStateOperations.h"
#include "ITSMFTTracking/MaterialPhysics.h"

#include <cmath>
#include <cstdint>

#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

namespace o2::itsmft::tracking
{
namespace
{

bool covarianceDiagonalsNonNegative(const SurfaceTrackState& state) noexcept
{
  for (uint8_t i = 0; i < 5; ++i) {
    if (state.covariance[packedCovarianceIndex(i, i)] < 0.f) {
      return false;
    }
  }
  return true;
}

// Physical-momentum derivation shared by both coordinate conventions. u = slot 4, t = slot 3.
bool derivePhysicalMomentum(const SurfaceTrackState& state, float& momentumGeV) noexcept
{
  const float t = state.parameters[3];
  const float u = state.parameters[4];
  const float absU = std::abs(u);
  const float pT = (state.absCharge == 0) ? (1.f / absU) : (static_cast<float>(state.absCharge) / absU);
  if (pT <= 0.f) {
    return false;
  }
  const float p = pT * std::sqrt(1.f + t * t);
  if (p <= 0.f) {
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

// Barrel covariance-range upper bound, in (Y, Z, Snp, Tgl, Q2Pt) slot order:
// the retained TrackParametrizationWithError<float>::checkCovariance()
// range-clamp values, and the same five constants
// PropagatorBarrelOperations.cxx's post-propagate/rotate/update
// sanitization (ADR 0008) enforces.
constexpr float kBarrelMaxDiagonal[5] = {o2::track::kCY2max, o2::track::kCZ2max, o2::track::kCSnp2max,
                                         o2::track::kCTgl2max, o2::track::kC1Pt2max};

// Thin wrapper over the shared, detector-neutral sanitizeCovariance()
// (SurfaceTrackState.h): abs()'s each diagonal and, if it still exceeds
// the retained maximum, clamps it and rescales every off-diagonal entry
// involving that parameter by sqrt(max/diagonal). No legacy track object is
// constructed; this operates directly on the packed float covariance array.
// Formerly a private reimplementation of this exact behavior; now delegates to
// the one shared implementation also used by the barrel state operations'
// own post-propagate/rotate/update sanitization, with no behavioral change.
void limitBarrelCovariance(SurfaceTrackState& scratch) noexcept
{
  sanitizeCovariance(scratch, kBarrelMaxDiagonal);
}

// Shared preflight validation, steps 1-6 of the required order. Step 3's
// kind-specific extra check (barrel |Snp|<1 / forward alpha==0) is
// supplied by the caller; the shared slot-4-nonzero part of step 3 is applied
// here for both kinds.
template <typename FamilyKinematicsCheck>
bool preflightValidate(const SurfaceTrackState& state, SurfaceKind expectedFamily, FamilyKinematicsCheck&& familyCheck,
                       material::MaterialFailureReason& failure) noexcept
{
  if (state.kind != expectedFamily) {
    failure = material::MaterialFailureReason::SourceSurfaceKindMismatch;
    return false;
  }
  const float u = state.parameters[4];
  if (!familyCheck(state) || u == 0.f) {
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

// Complete transactional operation shared by cylinder and disk states (Slice 2
// "Transactional result contract"): validate, derive physical momentum,
// scratch-copy, invoke the scalar kernel, project covariance on scratch
// only, validate the projected scratch, and commit exactly once.
// projectCovariance may additionally apply cylinder-specific covariance range
// handling (barrel only); it must not touch state.parameters[4], which this
// function updates uniformly for both kinds after projection.
//
// Unconditional no-op contract: once the scalar kernel succeeds, absCharge
// == 0 or an exactly-{0,0} materialBudget returns the scalar result
// immediately, before projectCovariance (and any barrel covariance-range
// limiting it applies) or the slot-4 update ever run. This holds even when
// the source state's barrel covariance diagonals already exceed the
// retained checkCovariance limits: those diagonals must not be silently
// clamped by an operation that has no material to apply.
template <typename FamilyKinematicsCheck, typename ProjectCovariance>
material::MaterialOperationResult correctForMaterialImpl(SurfaceTrackState& state, SurfaceKind expectedFamily,
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

  SurfaceTrackState scratch = state;
  const auto scalarResult = material::calculateMaterialPhysics(momentumBeforeGeV, scratch.pid, scratch.absCharge, direction, materialBudget);
  if (!scalarResult.ok()) {
    return scalarResult;
  }

  const bool isNoopMaterial = (materialBudget.xOverX0 == 0.f && materialBudget.arealDensityGPerCm2 == 0.f);
  if (scratch.absCharge == 0 || isNoopMaterial) {
    return scalarResult;
  }

  const float tBefore = scratch.parameters[3];
  const float kBefore = scratch.parameters[4];
  projectCovariance(scratch, scalarResult, tBefore, kBefore);

  // The equality branch preserves the exact no-op invariant for the
  // MCS-only-with-unchanged-momentum case (xOverX0 > 0, arealDensity == 0):
  // x == y implies kAfter == kBefore bit-for-bit with no division rounding.
  // The nonzero-change branch keeps the accepted/legacy left-to-right
  // arithmetic (multiply, then divide) rather than dividing the momenta
  // first, which would prematurely underflow for extreme momentum ratios
  // and would not reproduce the retained nonzero-material rounding.
  const float kAfter = (scalarResult.momentumBeforeGeV == scalarResult.momentumAfterGeV)
                         ? kBefore
                         : (kBefore * scalarResult.momentumBeforeGeV) / scalarResult.momentumAfterGeV;
  scratch.parameters[4] = kAfter;

  // Complete post-projection domain validation: the projected state must
  // still satisfy every kind/kinematics precondition the source state was
  // required to satisfy, and physical momentum must still be re-derivable.
  if (scratch.parameters[4] == 0.f || !familyCheck(scratch)) {
    return makeProjectionFailure(scalarResult, material::MaterialFailureReason::InvalidStateKinematics);
  }
  float momentumAfterDerived = 0.f;
  if (!derivePhysicalMomentum(scratch, momentumAfterDerived)) {
    return makeProjectionFailure(scalarResult, material::MaterialFailureReason::InvalidStateKinematics);
  }
  if (!covarianceDiagonalsNonNegative(scratch)) {
    return makeProjectionFailure(scalarResult, material::MaterialFailureReason::InvalidCovariance);
  }

  state = scratch;
  return scalarResult;
}

} // namespace
} // namespace o2::itsmft::tracking

namespace o2::itsmft::tracking::detail::barrel
{

material::MaterialOperationResult correctForMaterial(SurfaceTrackState& state, material::IntegratedMaterialBudget materialBudget,
                                                     material::MaterialTraversalDirection direction) noexcept
{
  auto familyCheck = [](const SurfaceTrackState& s) noexcept {
    return std::abs(s.parameters[2]) < 1.f;
  };
  // Barrel parameters are (Y, Z, Snp, Tgl, Q2Pt). The accepted Jacobian
  // requires q/pT unconditionally in slots 13/14, fixing the retained
  // TrackParametrizationWithError<float>::correctForMaterial() unit-charge
  // conditional that omits q/pT there (see the module doc comment).
  auto projectCovariance = [](SurfaceTrackState& scratch, const material::MaterialOperationResult& scalarResult,
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
  return correctForMaterialImpl(state, SurfaceKind::Cylinder, materialBudget, direction, familyCheck, projectCovariance);
}

} // namespace o2::itsmft::tracking::detail::barrel

namespace o2::itsmft::tracking::detail::forward
{

material::MaterialOperationResult correctForMaterial(SurfaceTrackState& state, material::IntegratedMaterialBudget materialBudget,
                                                     material::MaterialTraversalDirection direction) noexcept
{
  auto familyCheck = [](const SurfaceTrackState& s) noexcept {
    return s.alpha == 0.f;
  };
  // Forward parameters are (X, Y, Phi, Tanl, Q2Pt); unlike barrel there is no
  // cos(phi)-like factor on the angular diagonal term, and forward does not
  // inherit barrel-only covariance range limiting. Slot 13 (the Q2Pt/Tanl
  // cross term) and the k^2*R straggling contribution to slot 14 are new
  // physics: the legacy TrackParCovFwd::addMCSEffect() never populates
  // slot 13 and has no charge/PID/energy-loss awareness at all.
  auto projectCovariance = [](SurfaceTrackState& scratch, const material::MaterialOperationResult& scalarResult,
                              float t, float k) noexcept {
    const float A = 1.f + t * t;
    const float h = scalarResult.highlandTheta2Rad2;
    const float R = scalarResult.relativeInverseMomentumVariance;
    scratch.covariance[packedCovarianceIndex(2, 2)] += h * A;
    scratch.covariance[packedCovarianceIndex(3, 3)] += h * A * A;
    scratch.covariance[packedCovarianceIndex(4, 3)] += h * A * t * k;
    scratch.covariance[packedCovarianceIndex(4, 4)] += h * (t * k) * (t * k) + k * k * R;
  };
  return correctForMaterialImpl(state, SurfaceKind::Disk, materialBudget, direction, familyCheck, projectCovariance);
}

} // namespace o2::itsmft::tracking::detail::forward
