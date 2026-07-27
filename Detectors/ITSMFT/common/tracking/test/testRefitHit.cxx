// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFTRefitHit
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <cstring>
#include <limits>

#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/TransitionPolicyOperations.h"

/// Stage-B refit-primitive slice: focused coverage for refitHit<Tag>
/// (TransitionPolicyOperations.h/.cxx). This is deliberately unwired: there
/// is no production call site anywhere in this slice, and attachHit<Tag> is
/// not touched by any of this.
using namespace o2::itsmft::tracking;

namespace
{

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

// --- Barrel fixtures --------------------------------------------------

SurfaceKinematicState barrelState(uint8_t absCharge = 1, o2::track::PID pid = o2::track::PID::Pion)
{
  SurfaceKinematicState state{};
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.2f;
  state.parameters[3] = -0.35f;
  state.parameters[4] = 0.8f;
  state.referenceCoordinate = 4.f;
  state.alpha = 0.3f;
  state.family = StateFamily::Barrel;
  state.absCharge = absCharge;
  state.pid = pid;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
    }
  }
  return state;
}

SurfaceLinearizationReference barrelLinRef(const SurfaceKinematicState& state)
{
  SurfaceLinearizationReference ref{};
  BOOST_REQUIRE(makeLinearizationReference(state, ref));
  return ref;
}

SurfaceMeasurement barrelMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.frame.q = 2.5f;
  measurement.frame.frameAngle = 0.3f; // same alpha as barrelState(): no rotation needed
  measurement.frame.u = 0.8f;
  measurement.frame.v = -0.45f;
  measurement.covariance = {0.04f, 0.012f, 0.09f};
  return measurement;
}

constexpr float BarrelBz = 5.f;

NominalSurfaceMaterial barrelMaterial() { return NominalSurfaceMaterial{0.01f, 0.001f}; }

// Independent re-transcription of refitHit<CylinderCylinder>'s own
// documented call sequence (rotate, linRef-propagate, direction-aware
// material, chi2 gate, update, optional reference shift).
bool replayRefitHitBarrel(SurfaceKinematicState state, SurfaceLinearizationReference linRef,
                          const SurfaceMeasurement& measurement, const NominalSurfaceMaterial& material, float bz,
                          material::MaterialTraversalDirection direction, bool gateEnabled, float maxChi2, float chi2In,
                          bool shiftRef, SurfaceKinematicState& outState, SurfaceLinearizationReference& outRef,
                          float& outChi2, OperationFailureReason& reason)
{
  if (!barrel::rotate(state, linRef, measurement.frame.frameAngle, bz, reason)) {
    return false;
  }
  if (!barrel::propagate(state, linRef, measurement.frame.q, bz, reason)) {
    return false;
  }
  const auto materialResult = barrel::correctForMaterial(
    state, material::IntegratedMaterialBudget{material.xOverX0, material.arealDensityGPerCm2}, direction);
  if (!materialResult.ok()) {
    reason = OperationFailureReason::MaterialFailure;
    return false;
  }
  float predChi2 = 0.f;
  if (!barrel::predictedChi2(state, measurement, predChi2, reason)) {
    return false;
  }
  if (predChi2 < 0.f || !std::isfinite(predChi2)) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  if (gateEnabled && predChi2 > maxChi2) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  float updateChi2 = 0.f;
  if (!barrel::update(state, measurement, updateChi2, reason)) {
    return false;
  }
  float chi2 = chi2In + updateChi2;
  if (shiftRef) {
    if (!barrel::shiftReferenceToMeasurement(linRef, measurement, reason)) {
      return false;
    }
  }
  outState = state;
  outRef = linRef;
  outChi2 = chi2;
  return true;
}

// --- Disk fixtures ------------------------------------------------------

SurfaceKinematicState diskState(uint8_t absCharge = 1, o2::track::PID pid = o2::track::PID::Pion)
{
  SurfaceKinematicState state{};
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.35f;
  state.parameters[3] = -2.5f;
  state.parameters[4] = 0.8f;
  state.referenceCoordinate = -45.f;
  state.family = StateFamily::Forward;
  state.absCharge = absCharge;
  state.pid = pid;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
    }
  }
  return state;
}

SurfaceLinearizationReference diskLinRef(const SurfaceKinematicState& state)
{
  SurfaceLinearizationReference ref{};
  BOOST_REQUIRE(makeLinearizationReference(state, ref));
  return ref;
}

SurfaceMeasurement diskMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.global = {0.8f, -0.45f, -50.f};
  measurement.frame.q = -50.f;
  measurement.frame.u = 0.8f;
  measurement.frame.v = -0.45f;
  measurement.covariance = {0.04f, 0.f, 0.09f};
  return measurement;
}

constexpr float DiskBz = 5.f;

NominalSurfaceMaterial diskMaterial() { return NominalSurfaceMaterial{0.02f, 0.002f}; }

bool replayForwardPropagateAcceptedModelWithLinRef(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef,
                                                   float targetZ, float bz, OperationFailureReason& reason)
{
  if (std::abs(bz) > 0.01f) {
    return forward::propagate<forward::PropagationModel::Helix>(state, linRef, targetZ, bz, reason);
  }
  return forward::propagate<forward::PropagationModel::Linear>(state, linRef, targetZ, bz, reason);
}

bool replayRefitHitDisk(SurfaceKinematicState state, SurfaceLinearizationReference linRef,
                        const SurfaceMeasurement& measurement, const NominalSurfaceMaterial& material, float bz,
                        material::MaterialTraversalDirection direction, bool gateEnabled, float maxChi2, float chi2In,
                        bool shiftRef, SurfaceKinematicState& outState, SurfaceLinearizationReference& outRef,
                        float& outChi2, OperationFailureReason& reason)
{
  if (!replayForwardPropagateAcceptedModelWithLinRef(state, linRef, measurement.frame.q, bz, reason)) {
    return false;
  }
  const auto materialResult = forward::correctForMaterial(
    state, material::IntegratedMaterialBudget{material.xOverX0, material.arealDensityGPerCm2}, direction);
  if (!materialResult.ok()) {
    reason = OperationFailureReason::MaterialFailure;
    return false;
  }
  float predChi2 = 0.f;
  if (!forward::predictedChi2(state, measurement, predChi2, reason)) {
    return false;
  }
  if (predChi2 < 0.f || !std::isfinite(predChi2)) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  if (gateEnabled && predChi2 > maxChi2) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  float updateChi2 = 0.f;
  if (!forward::update(state, measurement, updateChi2, reason)) {
    return false;
  }
  float chi2 = chi2In + updateChi2;
  if (shiftRef) {
    if (!forward::shiftReferenceToMeasurement(linRef, measurement, reason)) {
      return false;
    }
  }
  outState = state;
  outRef = linRef;
  outChi2 = chi2;
  return true;
}

} // namespace

// ===========================================================================
// refitHit<CylinderCylinder>
// ===========================================================================

BOOST_AUTO_TEST_CASE(RefitHitBarrelMatchesIndependentReplay)
{
  const auto state = barrelState();
  const auto linRef = barrelLinRef(state);
  const auto measurement = barrelMeasurement();
  const auto material = barrelMaterial();

  SurfaceKinematicState expectedState{};
  SurfaceLinearizationReference expectedRef{};
  float expectedChi2 = 0.f;
  OperationFailureReason expectedReason{};
  BOOST_REQUIRE(replayRefitHitBarrel(state, linRef, measurement, material, BarrelBz,
                                     material::MaterialTraversalDirection::AlongMomentum, false, 0.f, 0.5f, true,
                                     expectedState, expectedRef, expectedChi2, expectedReason));

  auto actualState = state;
  auto actualRef = linRef;
  float actualChi2 = 0.5f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::CylinderCylinder>(
    actualState, actualRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::AlongMomentum,
    false, 0.f, actualChi2, true, reason));

  BOOST_CHECK(bitEqual(actualState, expectedState));
  BOOST_CHECK(bitEqual(actualRef, expectedRef));
  BOOST_CHECK_EQUAL(actualChi2, expectedChi2);
}

// The Kalman update step (shared identically by both directions) dominates
// the overall parameters[4] shift from its pre-fit value, so the material
// direction's own contribution (energy loss vs gain, order 1e-6 for this
// fixture's thin material) must be isolated by comparing the two full
// results directly against each other, not each one's own delta from the
// pre-fit value.
BOOST_AUTO_TEST_CASE(RefitHitBarrelDirectionAffectsQ2PtOppositely)
{
  const auto state = barrelState();
  const auto linRef = barrelLinRef(state);
  const auto measurement = barrelMeasurement();
  const auto material = barrelMaterial();

  auto alongState = state;
  auto alongRef = linRef;
  float alongChi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::CylinderCylinder>(
    alongState, alongRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::AlongMomentum,
    false, 0.f, alongChi2, false, reason));

  auto oppositeState = state;
  auto oppositeRef = linRef;
  float oppositeChi2 = 0.f;
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::CylinderCylinder>(
    oppositeState, oppositeRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::OppositeMomentum,
    false, 0.f, oppositeChi2, false, reason));

  // AlongMomentum loses energy (momentum decreases, |Q2Pt| increases);
  // OppositeMomentum gains energy (momentum increases, |Q2Pt| decreases) --
  // for this fixture's positive-sign Q2Pt, Along must end up strictly
  // larger than Opposite.
  BOOST_CHECK_GT(alongState.parameters[4], oppositeState.parameters[4]);
}

BOOST_AUTO_TEST_CASE(RefitHitBarrelAbsChargeZeroIsDirectionIndependent)
{
  const auto state = barrelState(0, o2::track::PID::Pion);
  const auto linRef = barrelLinRef(state);
  const auto measurement = barrelMeasurement();
  const auto material = barrelMaterial();

  auto alongState = state;
  auto alongRef = linRef;
  float alongChi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::CylinderCylinder>(
    alongState, alongRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::AlongMomentum,
    false, 0.f, alongChi2, false, reason));

  auto oppositeState = state;
  auto oppositeRef = linRef;
  float oppositeChi2 = 0.f;
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::CylinderCylinder>(
    oppositeState, oppositeRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::OppositeMomentum,
    false, 0.f, oppositeChi2, false, reason));

  BOOST_CHECK(bitEqual(alongState, oppositeState));
}

BOOST_AUTO_TEST_CASE(RefitHitBarrelPidAffectsResult)
{
  const auto measurement = barrelMeasurement();
  const auto material = barrelMaterial();

  const auto pionState = barrelState(1, o2::track::PID::Pion);
  auto pionResult = pionState;
  auto pionRef = barrelLinRef(pionState);
  float pionChi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::CylinderCylinder>(
    pionResult, pionRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::OppositeMomentum,
    false, 0.f, pionChi2, false, reason));

  const auto protonState = barrelState(1, o2::track::PID::Proton);
  auto protonResult = protonState;
  auto protonRef = barrelLinRef(protonState);
  float protonChi2 = 0.f;
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::CylinderCylinder>(
    protonResult, protonRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::OppositeMomentum,
    false, 0.f, protonChi2, false, reason));

  BOOST_CHECK_NE(pionResult.parameters[4], protonResult.parameters[4]);
}

BOOST_AUTO_TEST_CASE(RefitHitBarrelGateDisabledAcceptsLargeChi2)
{
  const auto state = barrelState();
  auto linRef = barrelLinRef(state);
  auto measurement = barrelMeasurement();
  measurement.frame.u = 500.f; // huge residual
  const auto material = barrelMaterial();

  auto resultState = state;
  float chi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_CHECK(refitHit<TransitionPolicyTag::CylinderCylinder>(
    resultState, linRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::AlongMomentum,
    false, 1.f, chi2, false, reason));
}

BOOST_AUTO_TEST_CASE(RefitHitBarrelGateAcceptedAndRejectedBoundary)
{
  const auto state = barrelState();
  const auto linRef = barrelLinRef(state);
  const auto measurement = barrelMeasurement();
  const auto material = barrelMaterial();

  // Discover the natural predChi2 for this fixture via the independent
  // replay, then set maxChi2 at exactly that value (inclusive accept) and
  // just below it (reject).
  SurfaceKinematicState replState{};
  SurfaceLinearizationReference replRef{};
  float replChi2 = 0.f;
  OperationFailureReason replReason{};
  BOOST_REQUIRE(replayRefitHitBarrel(state, linRef, measurement, material, BarrelBz,
                                     material::MaterialTraversalDirection::AlongMomentum, false, 0.f, 0.f, false,
                                     replState, replRef, replChi2, replReason));
  const float naturalChi2 = replChi2; // gate disabled above -> chi2In=0, so this equals the update chi2 increment

  auto acceptedState = state;
  auto acceptedRef = linRef;
  float acceptedChi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_CHECK(refitHit<TransitionPolicyTag::CylinderCylinder>(
    acceptedState, acceptedRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::AlongMomentum,
    true, naturalChi2 + 1.f, acceptedChi2, false, reason));

  auto rejectedState = state;
  auto rejectedRef = linRef;
  float rejectedChi2 = 0.f;
  BOOST_CHECK(!refitHit<TransitionPolicyTag::CylinderCylinder>(
    rejectedState, rejectedRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::AlongMomentum,
    true, -1.f, rejectedChi2, false, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(rejectedState, state));
  BOOST_CHECK(bitEqual(rejectedRef, linRef));
}

BOOST_AUTO_TEST_CASE(RefitHitBarrelRejectsNegativeChi2RegardlessOfGate)
{
  const auto state = barrelState();
  const auto linRef = barrelLinRef(state);
  const auto measurement = barrelMeasurement();
  const auto material = barrelMaterial();

  // Deterministically construct a measurement with residualZ == 0 (by
  // querying the exact post-rotate/propagate/material state) and a large
  // residualY combined with a deliberately invalid negative measurement
  // variance, forcing predictedChi2 negative.
  SurfaceKinematicState probe = state;
  SurfaceLinearizationReference probeRef = linRef;
  OperationFailureReason probeReason{};
  BOOST_REQUIRE(barrel::rotate(probe, probeRef, measurement.frame.frameAngle, BarrelBz, probeReason));
  BOOST_REQUIRE(barrel::propagate(probe, probeRef, measurement.frame.q, BarrelBz, probeReason));
  BOOST_REQUIRE(barrel::correctForMaterial(probe, material::IntegratedMaterialBudget{material.xOverX0, material.arealDensityGPerCm2},
                                           material::MaterialTraversalDirection::AlongMomentum)
                  .ok());

  auto badMeasurement = measurement;
  badMeasurement.frame.v = probe.parameters[1];
  badMeasurement.frame.u = probe.parameters[0] + 50.f;
  badMeasurement.covariance.uu = -0.5f;
  badMeasurement.covariance.uv = 0.f;
  badMeasurement.covariance.vv = 0.05f;

  auto resultState = state;
  auto resultRef = linRef;
  float chi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_CHECK(!refitHit<TransitionPolicyTag::CylinderCylinder>(
    resultState, resultRef, badMeasurement, material, BarrelBz, material::MaterialTraversalDirection::AlongMomentum,
    /*gateEnabled=*/true, /*maxChi2=*/1.e9f, chi2, false, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(resultState, state));
  BOOST_CHECK(bitEqual(resultRef, linRef));
  BOOST_CHECK_EQUAL(chi2, 0.f);
}

BOOST_AUTO_TEST_CASE(RefitHitBarrelNonMutationOnEveryFailureStage)
{
  const auto material = barrelMaterial();

  // Rotation failure.
  {
    auto state = barrelState();
    state.parameters[2] = 0.9999999f;
    auto linRef = barrelLinRef(state);
    auto measurement = barrelMeasurement();
    measurement.frame.frameAngle = state.alpha + o2::constants::math::PI; // half-turn: cos(delta)+snp*sin(delta)/csp < 0 with snp ~ 1
    const auto stateBefore = state;
    const auto refBefore = linRef;
    float chi2 = 0.7f;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    BOOST_CHECK(!refitHit<TransitionPolicyTag::CylinderCylinder>(
      state, linRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::AlongMomentum, false, 0.f,
      chi2, false, reason));
    BOOST_CHECK(bitEqual(state, stateBefore));
    BOOST_CHECK(bitEqual(linRef, refBefore));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  }

  // Non-finite target (propagation stage).
  {
    auto state = barrelState();
    auto linRef = barrelLinRef(state);
    auto measurement = barrelMeasurement();
    measurement.frame.q = std::numeric_limits<float>::quiet_NaN();
    const auto stateBefore = state;
    const auto refBefore = linRef;
    float chi2 = 0.3f;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    BOOST_CHECK(!refitHit<TransitionPolicyTag::CylinderCylinder>(
      state, linRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::AlongMomentum, false, 0.f,
      chi2, false, reason));
    BOOST_CHECK(bitEqual(state, stateBefore));
    BOOST_CHECK(bitEqual(linRef, refBefore));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  }

  // Material failure (non-finite material budget).
  {
    auto state = barrelState();
    auto linRef = barrelLinRef(state);
    auto measurement = barrelMeasurement();
    const NominalSurfaceMaterial badMaterial{std::numeric_limits<float>::quiet_NaN(), 0.f};
    const auto stateBefore = state;
    const auto refBefore = linRef;
    float chi2 = 0.1f;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    BOOST_CHECK(!refitHit<TransitionPolicyTag::CylinderCylinder>(
      state, linRef, measurement, badMaterial, BarrelBz, material::MaterialTraversalDirection::AlongMomentum, false, 0.f,
      chi2, false, reason));
    BOOST_CHECK(reason == OperationFailureReason::MaterialFailure);
    BOOST_CHECK(bitEqual(state, stateBefore));
    BOOST_CHECK(bitEqual(linRef, refBefore));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  }
}

BOOST_AUTO_TEST_CASE(RefitHitBarrelShiftReferenceOnOff)
{
  const auto state = barrelState();
  const auto linRef = barrelLinRef(state);
  const auto measurement = barrelMeasurement();
  const auto material = barrelMaterial();

  auto onState = state;
  auto onRef = linRef;
  float onChi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::CylinderCylinder>(
    onState, onRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::AlongMomentum, false, 0.f,
    onChi2, true, reason));
  BOOST_CHECK_EQUAL(onRef.parameters[0], measurement.frame.u);
  BOOST_CHECK_EQUAL(onRef.parameters[1], measurement.frame.v);

  auto offState = state;
  auto offRef = linRef;
  float offChi2 = 0.f;
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::CylinderCylinder>(
    offState, offRef, measurement, material, BarrelBz, material::MaterialTraversalDirection::AlongMomentum, false, 0.f,
    offChi2, false, reason));
  BOOST_CHECK_NE(offRef.parameters[0], measurement.frame.u);
  BOOST_CHECK(bitEqual(onState, offState)); // shift never touches the fitted state
}

// ===========================================================================
// refitHit<DiskDisk>: same orchestration contract as CylinderCylinder above.
// ===========================================================================

BOOST_AUTO_TEST_CASE(RefitHitDiskMatchesIndependentReplay)
{
  const auto state = diskState();
  const auto linRef = diskLinRef(state);
  const auto measurement = diskMeasurement();
  const auto material = diskMaterial();

  SurfaceKinematicState expectedState{};
  SurfaceLinearizationReference expectedRef{};
  float expectedChi2 = 0.f;
  OperationFailureReason expectedReason{};
  BOOST_REQUIRE(replayRefitHitDisk(state, linRef, measurement, material, DiskBz,
                                   material::MaterialTraversalDirection::OppositeMomentum, false, 0.f, 0.25f, true,
                                   expectedState, expectedRef, expectedChi2, expectedReason));

  auto actualState = state;
  auto actualRef = linRef;
  float actualChi2 = 0.25f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::DiskDisk>(
    actualState, actualRef, measurement, material, DiskBz, material::MaterialTraversalDirection::OppositeMomentum,
    false, 0.f, actualChi2, true, reason));

  BOOST_CHECK(bitEqual(actualState, expectedState));
  BOOST_CHECK(bitEqual(actualRef, expectedRef));
  BOOST_CHECK_EQUAL(actualChi2, expectedChi2);
}

BOOST_AUTO_TEST_CASE(RefitHitDiskDirectionAffectsQ2PtOppositely)
{
  const auto state = diskState();
  const auto linRef = diskLinRef(state);
  const auto measurement = diskMeasurement();
  const auto material = diskMaterial();

  auto alongState = state;
  auto alongRef = linRef;
  float alongChi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::DiskDisk>(
    alongState, alongRef, measurement, material, DiskBz, material::MaterialTraversalDirection::AlongMomentum, false,
    0.f, alongChi2, false, reason));

  auto oppositeState = state;
  auto oppositeRef = linRef;
  float oppositeChi2 = 0.f;
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::DiskDisk>(
    oppositeState, oppositeRef, measurement, material, DiskBz, material::MaterialTraversalDirection::OppositeMomentum,
    false, 0.f, oppositeChi2, false, reason));

  // Same isolation rationale as RefitHitBarrelDirectionAffectsQ2PtOppositely
  // above: the shared Kalman update dominates the delta from the pre-fit
  // value, so compare the two full results against each other directly.
  BOOST_CHECK_GT(alongState.parameters[4], oppositeState.parameters[4]);
}

BOOST_AUTO_TEST_CASE(RefitHitDiskGateAcceptedAndRejectedBoundary)
{
  const auto state = diskState();
  const auto linRef = diskLinRef(state);
  const auto measurement = diskMeasurement();
  const auto material = diskMaterial();

  SurfaceKinematicState replState{};
  SurfaceLinearizationReference replRef{};
  float replChi2 = 0.f;
  OperationFailureReason replReason{};
  BOOST_REQUIRE(replayRefitHitDisk(state, linRef, measurement, material, DiskBz,
                                   material::MaterialTraversalDirection::OppositeMomentum, false, 0.f, 0.f, false,
                                   replState, replRef, replChi2, replReason));
  const float naturalChi2 = replChi2;

  auto acceptedState = state;
  auto acceptedRef = linRef;
  float acceptedChi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_CHECK(refitHit<TransitionPolicyTag::DiskDisk>(
    acceptedState, acceptedRef, measurement, material, DiskBz, material::MaterialTraversalDirection::OppositeMomentum,
    true, naturalChi2 + 1.f, acceptedChi2, false, reason));

  auto rejectedState = state;
  auto rejectedRef = linRef;
  float rejectedChi2 = 0.f;
  BOOST_CHECK(!refitHit<TransitionPolicyTag::DiskDisk>(
    rejectedState, rejectedRef, measurement, material, DiskBz, material::MaterialTraversalDirection::OppositeMomentum,
    true, -1.f, rejectedChi2, false, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(rejectedState, state));
  BOOST_CHECK(bitEqual(rejectedRef, linRef));
}

// Unlike attachHit<DiskDisk> (no `< 0.f` guard at all), refitHit<DiskDisk>
// unconditionally rejects a negative predicted chi2 -- the same uniform
// rule as CylinderCylinder above, with no family-specific branch.
BOOST_AUTO_TEST_CASE(RefitHitDiskRejectsNegativeChi2RegardlessOfGate)
{
  const auto state = diskState();
  const auto linRef = diskLinRef(state);
  const auto measurement = diskMeasurement();
  const auto material = diskMaterial();

  SurfaceKinematicState probe = state;
  SurfaceLinearizationReference probeRef = linRef;
  OperationFailureReason probeReason{};
  BOOST_REQUIRE(replayForwardPropagateAcceptedModelWithLinRef(probe, probeRef, measurement.frame.q, DiskBz, probeReason));
  BOOST_REQUIRE(forward::correctForMaterial(probe, material::IntegratedMaterialBudget{material.xOverX0, material.arealDensityGPerCm2},
                                            material::MaterialTraversalDirection::OppositeMomentum)
                  .ok());

  auto badMeasurement = measurement;
  badMeasurement.global.y = probe.parameters[1];
  badMeasurement.frame.v = probe.parameters[1];
  badMeasurement.global.x = probe.parameters[0] + 50.f;
  badMeasurement.frame.u = probe.parameters[0] + 50.f;
  badMeasurement.covariance.uu = -0.5f;
  badMeasurement.covariance.uv = 0.f;
  badMeasurement.covariance.vv = 0.05f;

  auto resultState = state;
  auto resultRef = linRef;
  float chi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_CHECK(!refitHit<TransitionPolicyTag::DiskDisk>(
    resultState, resultRef, badMeasurement, material, DiskBz, material::MaterialTraversalDirection::OppositeMomentum,
    /*gateEnabled=*/true, /*maxChi2=*/1.e9f, chi2, false, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(resultState, state));
  BOOST_CHECK(bitEqual(resultRef, linRef));
  BOOST_CHECK_EQUAL(chi2, 0.f);
}

BOOST_AUTO_TEST_CASE(RefitHitDiskNonMutationOnMaterialFailure)
{
  auto state = diskState();
  auto linRef = diskLinRef(state);
  auto measurement = diskMeasurement();
  const NominalSurfaceMaterial badMaterial{std::numeric_limits<float>::quiet_NaN(), 0.f};
  const auto stateBefore = state;
  const auto refBefore = linRef;
  float chi2 = 0.2f;
  const float chi2Before = chi2;
  OperationFailureReason reason{};
  BOOST_CHECK(!refitHit<TransitionPolicyTag::DiskDisk>(
    state, linRef, measurement, badMaterial, DiskBz, material::MaterialTraversalDirection::AlongMomentum, false, 0.f,
    chi2, false, reason));
  BOOST_CHECK(reason == OperationFailureReason::MaterialFailure);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, refBefore));
  BOOST_CHECK_EQUAL(chi2, chi2Before);
}

BOOST_AUTO_TEST_CASE(RefitHitDiskShiftReferenceOnOff)
{
  const auto state = diskState();
  const auto linRef = diskLinRef(state);
  const auto measurement = diskMeasurement();
  const auto material = diskMaterial();

  auto onState = state;
  auto onRef = linRef;
  float onChi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::DiskDisk>(
    onState, onRef, measurement, material, DiskBz, material::MaterialTraversalDirection::OppositeMomentum, false, 0.f,
    onChi2, true, reason));
  BOOST_CHECK_EQUAL(onRef.parameters[0], measurement.frame.u);
  BOOST_CHECK_EQUAL(onRef.parameters[1], measurement.frame.v);

  auto offState = state;
  auto offRef = linRef;
  float offChi2 = 0.f;
  BOOST_REQUIRE(refitHit<TransitionPolicyTag::DiskDisk>(
    offState, offRef, measurement, material, DiskBz, material::MaterialTraversalDirection::OppositeMomentum, false, 0.f,
    offChi2, false, reason));
  BOOST_CHECK_NE(offRef.parameters[0], measurement.frame.u);
  BOOST_CHECK(bitEqual(onState, offState));
}
