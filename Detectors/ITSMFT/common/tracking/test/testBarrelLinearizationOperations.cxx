// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFTBarrelLinearizationOperations
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/detail/SurfaceKinematicStateLegacyAdapters.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"

namespace
{
using namespace o2::itsmft::tracking;
using namespace o2::itsmft::tracking::barrel;

constexpr float AbsoluteTolerance = 3.e-5f;
constexpr float RelativeTolerance = 3.e-4f;

struct Drift {
  float absolute{0.f};
  float relative{0.f};
};

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

void checkClose(float value, float oracle, Drift& drift)
{
  const float difference = std::abs(value - oracle);
  drift.absolute = std::max(drift.absolute, difference);
  if (oracle != 0.f) {
    drift.relative = std::max(drift.relative, difference / std::abs(oracle));
  }
  BOOST_CHECK_LE(difference, AbsoluteTolerance + RelativeTolerance * std::abs(oracle));
}

SurfaceKinematicState makeState()
{
  SurfaceKinematicState state{};
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.2f;
  state.parameters[3] = -0.35f;
  state.parameters[4] = 0.8f;
  state.referenceCoordinate = 4.f;
  state.alpha = 0.3f;
  state.kind = SurfaceKind::Cylinder;
  state.flags = 0x5a;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
    }
  }
  return state;
}

SurfaceLinearizationReference linRefFromState(const SurfaceKinematicState& state)
{
  SurfaceLinearizationReference ref{};
  BOOST_REQUIRE(makeLinearizationReference(state, ref));
  return ref;
}

// A linearization reference deliberately offset from `state`'s own
// parameters, used to prove reference-driven and self-driven propagation
// differ.
SurfaceLinearizationReference perturbedLinRef(const SurfaceKinematicState& state)
{
  auto ref = linRefFromState(state);
  ref.parameters[0] += 0.08f;
  ref.parameters[1] -= 0.05f;
  ref.parameters[2] += 0.01f;
  ref.parameters[3] += 0.015f;
  ref.parameters[4] -= 0.01f;
  return ref;
}

bool exportLinRef(const SurfaceLinearizationReference& ref, uint8_t absCharge, o2::track::TrackParF& out)
{
  if (ref.kind != SurfaceKind::Cylinder) {
    return false;
  }
  const o2::track::TrackParF::params_t params{ref.parameters[0], ref.parameters[1], ref.parameters[2], ref.parameters[3], ref.parameters[4]};
  out = o2::track::TrackParF(ref.referenceCoordinate, ref.alpha, params, absCharge);
  return true;
}

void compareLinRefWithOracle(const SurfaceLinearizationReference& ref, const o2::track::TrackParF& oracle, Drift& drift)
{
  BOOST_CHECK(ref.kind == SurfaceKind::Cylinder);
  checkClose(ref.referenceCoordinate, oracle.getX(), drift);
  checkClose(ref.alpha, oracle.getAlpha(), drift);
  for (uint8_t i = 0; i < 5; ++i) {
    checkClose(ref.parameters[i], oracle.getParam(i), drift);
  }
}

void compareStateWithOracle(const SurfaceKinematicState& state, const o2::track::TrackParCovF& oracle, Drift& drift)
{
  checkClose(state.referenceCoordinate, oracle.getX(), drift);
  checkClose(state.alpha, oracle.getAlpha(), drift);
  for (uint8_t i = 0; i < 5; ++i) {
    checkClose(state.parameters[i], oracle.getParam(i), drift);
  }
  for (uint8_t i = 0; i < 15; ++i) {
    checkClose(state.covariance[i], oracle.getCov()[i], drift);
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(RotateWithLinRefMatchesRetainedLegacyOracle)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);
  const auto stateBefore = state;
  const auto linRefBefore = linRef;

  o2::track::TrackParCovF oracleState{};
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(state, oracleState));
  o2::track::TrackParF oracleRef{};
  BOOST_REQUIRE(exportLinRef(linRef, state.absCharge, oracleRef));

  const float bz = -5.f;
  const float targetAlpha = 4.f * o2::constants::math::PI + 0.55f;
  BOOST_REQUIRE(oracleState.rotate(targetAlpha, oracleRef, bz));

  OperationFailureReason reason{};
  BOOST_REQUIRE(rotate(state, linRef, targetAlpha, bz, reason));

  Drift drift{};
  compareStateWithOracle(state, oracleState, drift);
  compareLinRefWithOracle(linRef, oracleRef, drift);
  // Metadata (family/flags/absCharge/pid) never touched by rotate.
  BOOST_CHECK(state.kind == stateBefore.kind);
  BOOST_CHECK_EQUAL(state.flags, stateBefore.flags);
  BOOST_CHECK_EQUAL(state.absCharge, stateBefore.absCharge);
}

BOOST_AUTO_TEST_CASE(PropagateWithLinRefMatchesRetainedLegacyOracle)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);

  o2::track::TrackParCovF oracleState{};
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(state, oracleState));
  o2::track::TrackParF oracleRef{};
  BOOST_REQUIRE(exportLinRef(linRef, state.absCharge, oracleRef));

  const float bz = -5.f;
  const float targetX = 6.5f;
  BOOST_REQUIRE(oracleState.propagateTo(targetX, oracleRef, bz));

  OperationFailureReason reason{};
  BOOST_REQUIRE(propagate(state, linRef, targetX, bz, reason));

  Drift drift{};
  compareStateWithOracle(state, oracleState, drift);
  compareLinRefWithOracle(linRef, oracleRef, drift);
}

BOOST_AUTO_TEST_CASE(PropagateWithLinRefMatchesOracleZeroField)
{
  auto state = makeState();
  state.absCharge = 0;
  auto linRef = linRefFromState(state);

  o2::track::TrackParCovF oracleState{};
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(state, oracleState));
  o2::track::TrackParF oracleRef{};
  BOOST_REQUIRE(exportLinRef(linRef, state.absCharge, oracleRef));

  const float bz = 0.f;
  const float targetX = 5.2f;
  BOOST_REQUIRE(oracleState.propagateTo(targetX, oracleRef, bz));

  OperationFailureReason reason{};
  BOOST_REQUIRE(propagate(state, linRef, targetX, bz, reason));

  Drift drift{};
  compareStateWithOracle(state, oracleState, drift);
  compareLinRefWithOracle(linRef, oracleRef, drift);
}

BOOST_AUTO_TEST_CASE(PropagateWithLinRefTrivialStepUpdatesOnlyReferenceCoordinate)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);
  const auto stateBefore = state;
  const auto linRefBefore = linRef;

  OperationFailureReason reason{};
  BOOST_REQUIRE(propagate(state, linRef, state.referenceCoordinate, -5.f, reason));
  BOOST_CHECK_EQUAL(state.referenceCoordinate, stateBefore.referenceCoordinate);
  BOOST_CHECK_EQUAL(linRef.referenceCoordinate, linRefBefore.referenceCoordinate);
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_EQUAL(state.parameters[i], stateBefore.parameters[i]);
    BOOST_CHECK_EQUAL(linRef.parameters[i], linRefBefore.parameters[i]);
  }
}

// Reference-driven propagation (Jacobian evaluated at a perturbed linRef)
// must differ from self-driven propagation (Jacobian evaluated at state's
// own snp via the non-linRef propagate) of the same starting state.
BOOST_AUTO_TEST_CASE(ReferenceDrivenPropagationDiffersFromSelfDrivenPropagation)
{
  const auto baseState = makeState();
  const float bz = -5.f;
  const float targetX = 6.5f;

  auto selfDriven = baseState;
  OperationFailureReason reason{};
  BOOST_REQUIRE(barrel::propagate(selfDriven, targetX, bz, reason));

  auto referenceDrivenState = baseState;
  auto referenceDrivenRef = perturbedLinRef(baseState);
  BOOST_REQUIRE(propagate(referenceDrivenState, referenceDrivenRef, targetX, bz, reason));

  bool differs = false;
  for (uint8_t i = 0; i < 15; ++i) {
    if (std::abs(selfDriven.covariance[i] - referenceDrivenState.covariance[i]) > 1.e-6f) {
      differs = true;
      break;
    }
  }
  BOOST_CHECK(differs);

  // Cross-check the reference-driven result against the real oracle with
  // the SAME perturbed linRef, proving the divergence above is expected
  // (not a bug), not merely "some difference exists".
  auto oracleStateBefore = baseState;
  o2::track::TrackParCovF oracleState{};
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(oracleStateBefore, oracleState));
  o2::track::TrackParF oracleRef{};
  BOOST_REQUIRE(exportLinRef(perturbedLinRef(baseState), baseState.absCharge, oracleRef));
  BOOST_REQUIRE(oracleState.propagateTo(targetX, oracleRef, bz));
  Drift drift{};
  compareStateWithOracle(referenceDrivenState, oracleState, drift);
}

BOOST_AUTO_TEST_CASE(RotateWithLinRefRejectsFamilyMismatchAndPreservesBytes)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);
  linRef.kind = SurfaceKind::Disk;
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!rotate(state, linRef, 0.5f, -5.f, reason));
  BOOST_CHECK(reason == OperationFailureReason::SourceSurfaceKindMismatch);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(RotateWithLinRefRejectsInvalidSnpAndPreservesBytes)
{
  auto state = makeState();
  state.parameters[2] = 0.9999999f; // snp ~ 1: local direction nearly perpendicular to X
  auto linRef = linRefFromState(state);
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  // A half-turn (delta = pi) makes cos(delta) = -1, so
  // csp*cos(delta) + snp*sin(delta) ~= -csp < 0 for any snp close to 1,
  // reliably violating the post-rotation local-cosine precondition.
  const float targetAlpha = state.alpha + o2::constants::math::PI;
  BOOST_CHECK(!rotate(state, linRef, targetAlpha, -5.f, reason));
  BOOST_CHECK(reason == OperationFailureReason::RotationFailure);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(PropagateWithLinRefRejectsNonFiniteInputAndPreservesBytes)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!propagate(state, linRef, std::numeric_limits<float>::quiet_NaN(), -5.f, reason));
  BOOST_CHECK(reason == OperationFailureReason::NonFiniteInput);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(PropagateWithLinRefRejectsUnreachableTargetAndPreservesBytes)
{
  auto state = makeState();
  state.parameters[2] = 0.9999999f;
  auto linRef = linRefFromState(state);
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!propagate(state, linRef, 50.f, -5.f, reason));
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(ShiftReferenceToMeasurementSetsOnlyYAndZ)
{
  auto linRef = linRefFromState(makeState());
  const auto before = linRef;
  SurfaceMeasurement measurement{};
  measurement.frame.u = 1.23f;
  measurement.frame.v = -4.56f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(shiftReferenceToMeasurement(linRef, measurement, reason));
  BOOST_CHECK_EQUAL(linRef.parameters[0], measurement.frame.u);
  BOOST_CHECK_EQUAL(linRef.parameters[1], measurement.frame.v);
  BOOST_CHECK_EQUAL(linRef.parameters[2], before.parameters[2]);
  BOOST_CHECK_EQUAL(linRef.parameters[3], before.parameters[3]);
  BOOST_CHECK_EQUAL(linRef.parameters[4], before.parameters[4]);
  BOOST_CHECK_EQUAL(linRef.referenceCoordinate, before.referenceCoordinate);
  BOOST_CHECK_EQUAL(linRef.alpha, before.alpha);
}

BOOST_AUTO_TEST_CASE(ShiftReferenceToMeasurementRejectsNonFiniteAndPreservesBytes)
{
  auto linRef = linRefFromState(makeState());
  const auto before = linRef;
  SurfaceMeasurement measurement{};
  measurement.frame.u = std::numeric_limits<float>::quiet_NaN();
  measurement.frame.v = 1.f;
  OperationFailureReason reason{};
  BOOST_CHECK(!shiftReferenceToMeasurement(linRef, measurement, reason));
  BOOST_CHECK(reason == OperationFailureReason::NonFiniteInput);
  BOOST_CHECK(bitEqual(linRef, before));
}

BOOST_AUTO_TEST_CASE(ShiftReferenceToMeasurementRejectsFamilyMismatch)
{
  auto linRef = linRefFromState(makeState());
  linRef.kind = SurfaceKind::Disk;
  const auto before = linRef;
  SurfaceMeasurement measurement{};
  measurement.frame.u = 1.f;
  measurement.frame.v = 2.f;
  OperationFailureReason reason{};
  BOOST_CHECK(!shiftReferenceToMeasurement(linRef, measurement, reason));
  BOOST_CHECK(bitEqual(linRef, before));
}

// --- Fitted-state/linRef pairing preconditions (hardening) -----------------
// Parameters may legitimately differ between `state` and `linRef` -- that is
// the entire purpose of a linearization reference -- but referenceCoordinate
// and (for Barrel) alpha must match exactly.

BOOST_AUTO_TEST_CASE(RotateWithLinRefRejectsReferenceCoordinateMismatchAndPreservesBytes)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);
  linRef.referenceCoordinate += 0.5f;
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!rotate(state, linRef, 0.5f, -5.f, reason));
  BOOST_CHECK(reason == OperationFailureReason::ReferenceCoordinateMismatch);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(RotateWithLinRefRejectsAlphaMismatchAndPreservesBytes)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);
  linRef.alpha += 0.1f;
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!rotate(state, linRef, 0.5f, -5.f, reason));
  BOOST_CHECK(reason == OperationFailureReason::AlphaMismatch);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(PropagateWithLinRefRejectsReferenceCoordinateMismatchAndPreservesBytes)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);
  linRef.referenceCoordinate += 0.5f;
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!propagate(state, linRef, 6.5f, -5.f, reason));
  BOOST_CHECK(reason == OperationFailureReason::ReferenceCoordinateMismatch);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(PropagateWithLinRefRejectsAlphaMismatchAndPreservesBytes)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);
  linRef.alpha += 0.1f;
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!propagate(state, linRef, 6.5f, -5.f, reason));
  BOOST_CHECK(reason == OperationFailureReason::AlphaMismatch);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}
