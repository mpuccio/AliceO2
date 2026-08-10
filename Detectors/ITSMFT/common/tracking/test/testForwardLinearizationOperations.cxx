// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFTForwardLinearizationOperations
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"

namespace
{
using namespace o2::itsmft::tracking;
using namespace o2::itsmft::tracking::forward;

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

constexpr float AbsoluteTolerance = 3.e-4f;
constexpr float RelativeTolerance = 3.e-4f;

struct Drift {
  float maximumAbsolute{0.f};
  float maximumRelative{0.f};
};

// Module-scope drift accumulators (item 4 hardening): every numerical
// cross-check below folds its observed drift into one of these, and the
// final BOOST_AUTO_TEST_CASE at the bottom of the file reports the maxima
// -- Boost.Test executes cases in registration order within one
// translation unit, so this reliably runs last.
Drift selfDrivenReductionParameterDrift;
Drift selfDrivenReductionCovarianceDrift;
Drift jacobianFiniteDifferenceDrift;

void checkClose(float value, float oracle, Drift& drift)
{
  const float difference = std::abs(value - oracle);
  drift.maximumAbsolute = std::max(drift.maximumAbsolute, difference);
  if (oracle != 0.f) {
    drift.maximumRelative = std::max(drift.maximumRelative, difference / std::abs(oracle));
  }
  BOOST_CHECK_LE(difference, AbsoluteTolerance + RelativeTolerance * std::abs(oracle));
}

SurfaceKinematicState makeState()
{
  SurfaceKinematicState state{};
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.35f;
  state.parameters[3] = -2.5f;
  state.parameters[4] = 0.8f;
  state.referenceCoordinate = -45.f;
  state.family = StateFamily::Forward;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] =
        row == column ? 0.01f * (row + 1) : 0.0002f * (1 + packedCovarianceIndex(row, column));
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

SurfaceLinearizationReference perturbedLinRef(const SurfaceKinematicState& state)
{
  auto ref = linRefFromState(state);
  ref.parameters[0] -= 0.09f;
  ref.parameters[1] += 0.06f;
  ref.parameters[2] -= 0.02f;
  ref.parameters[3] += 0.015f;
  ref.parameters[4] -= 0.01f;
  return ref;
}

// Independent analytic float oracle for the Linear model, written
// separately from ForwardSurfaceStateOperations.cxx's referencePropagateLinear
// (different structure/variable names, same underlying secant-line physics)
// so this test does not merely echo the production formula back at itself.
void computeExpectedLinear(const SurfaceKinematicState& state, const SurfaceLinearizationReference& ref, float targetZ,
                           SurfaceKinematicState& expectedState, SurfaceLinearizationReference& expectedRef)
{
  const float step = targetZ - ref.referenceCoordinate;
  const float cotTanl = 1.f / ref.parameters[3];
  const float alongZ = step * cotTanl;
  const float alongZ2 = alongZ * cotTanl;
  const float sinPhi = std::sin(ref.parameters[2]);
  const float cosPhi = std::cos(ref.parameters[2]);

  expectedRef = ref;
  expectedRef.referenceCoordinate = targetZ;
  expectedRef.parameters[0] = ref.parameters[0] + alongZ * cosPhi;
  expectedRef.parameters[1] = ref.parameters[1] + alongZ * sinPhi;

  float jacobian[5][5] = {};
  for (uint8_t i = 0; i < 5; ++i) {
    jacobian[i][i] = 1.f;
  }
  jacobian[0][2] = -alongZ * sinPhi;
  jacobian[0][3] = -alongZ2 * cosPhi;
  jacobian[1][2] = alongZ * cosPhi;
  jacobian[1][3] = -alongZ2 * sinPhi;

  float diff[5];
  for (uint8_t i = 0; i < 5; ++i) {
    diff[i] = state.parameters[i] - ref.parameters[i];
  }

  expectedState = state;
  expectedState.referenceCoordinate = targetZ;
  for (uint8_t row = 0; row < 5; ++row) {
    float value = expectedRef.parameters[row];
    for (uint8_t column = 0; column < 5; ++column) {
      value += jacobian[row][column] * diff[column];
    }
    expectedState.parameters[row] = value;
  }
}

// State == linRef reduction: when the reference is an exact, unperturbed
// copy of the fitted state, the linRef-aware propagation must agree with the
// independently oracle-tested non-linRef propagation for both parameters and
// covariance -- not merely to first order, but to the tolerance floor of
// this file's shared float arithmetic (the two paths compute the position
// update and Jacobian from identical starting values and feed the same
// covariance transform).
void checkSelfDrivenReduction(float bz, float targetZ)
{
  const auto baseState = makeState();

  auto linRefDriven = baseState;
  auto linRef = linRefFromState(baseState);
  OperationFailureReason reason{};
  BOOST_REQUIRE(Propagator::propagateForward(linRefDriven, linRef, targetZ, bz, reason));

  auto selfDriven = baseState;
  BOOST_REQUIRE(Propagator::propagateForward(selfDriven, targetZ, bz, reason));

  checkClose(linRefDriven.referenceCoordinate, selfDriven.referenceCoordinate, selfDrivenReductionParameterDrift);
  for (uint8_t i = 0; i < 5; ++i) {
    checkClose(linRefDriven.parameters[i], selfDriven.parameters[i], selfDrivenReductionParameterDrift);
  }
  for (uint8_t i = 0; i < 15; ++i) {
    checkClose(linRefDriven.covariance[i], selfDriven.covariance[i], selfDrivenReductionCovarianceDrift);
  }
}

// Finite-difference Jacobian check (item 4 hardening): perturbs exactly one
// of the five parameters of `state` away from an otherwise-identical
// `linRef`, and compares the linRef-aware propagation parameter
// output against a first-order Taylor prediction built from an
// INDEPENDENT central-difference numerical derivative -- computed purely
// from two calls to the already-accepted, non-linRef propagation
// (never this slice's own analytic Jacobian formulas). This validates the
// analytic Jacobian against a numerically-derived one.
void checkJacobianAgainstFiniteDifference(float bz, float targetZ, uint8_t paramIndex)
{
  const auto baseState = makeState();
  const float epsilon = 1.e-3f;

  auto plusSelf = baseState;
  plusSelf.parameters[paramIndex] += epsilon;
  auto minusSelf = baseState;
  minusSelf.parameters[paramIndex] -= epsilon;
  OperationFailureReason reason{};
  BOOST_REQUIRE(Propagator::propagateForward(plusSelf, targetZ, bz, reason));
  BOOST_REQUIRE(Propagator::propagateForward(minusSelf, targetZ, bz, reason));

  auto baseSelf = baseState;
  BOOST_REQUIRE(Propagator::propagateForward(baseSelf, targetZ, bz, reason));

  float numericalDerivative[5];
  for (uint8_t row = 0; row < 5; ++row) {
    numericalDerivative[row] = (plusSelf.parameters[row] - minusSelf.parameters[row]) / (2.f * epsilon);
  }

  auto perturbedState = baseState;
  perturbedState.parameters[paramIndex] += epsilon;
  auto linRef = linRefFromState(baseState);
  BOOST_REQUIRE(Propagator::propagateForward(perturbedState, linRef, targetZ, bz, reason));

  for (uint8_t row = 0; row < 5; ++row) {
    const float predicted = baseSelf.parameters[row] + numericalDerivative[row] * epsilon;
    checkClose(perturbedState.parameters[row], predicted, jacobianFiniteDifferenceDrift);
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(LowFieldPropagationMatchesIndependentAnalyticOracle)
{
  const auto baseState = makeState();
  auto state = baseState;
  state.parameters[0] += 0.05f;
  state.parameters[2] -= 0.01f;
  state.parameters[3] += 0.02f;
  auto linRef = linRefFromState(baseState);

  SurfaceKinematicState expectedState{};
  SurfaceLinearizationReference expectedRef{};
  computeExpectedLinear(state, linRef, -80.f, expectedState, expectedRef);

  OperationFailureReason reason{};
  BOOST_REQUIRE(Propagator::propagateForward(state, linRef, -80.f, 0.f, reason));

  BOOST_CHECK_CLOSE(state.referenceCoordinate, expectedState.referenceCoordinate, 1.e-3f);
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_SMALL(state.parameters[i] - expectedState.parameters[i], 1.e-4f);
    BOOST_CHECK_SMALL(linRef.parameters[i] - expectedRef.parameters[i], 1.e-4f);
  }
}

BOOST_AUTO_TEST_CASE(LowFieldPropagationMatchesOracleAtZeroField)
{
  const auto baseState = makeState();
  auto state = baseState;
  state.parameters[1] -= 0.03f;
  auto linRef = linRefFromState(baseState);

  SurfaceKinematicState expectedState{};
  SurfaceLinearizationReference expectedRef{};
  computeExpectedLinear(state, linRef, -60.f, expectedState, expectedRef);

  OperationFailureReason reason{};
  BOOST_REQUIRE(Propagator::propagateForward(state, linRef, -60.f, 0.f, reason));
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_SMALL(state.parameters[i] - expectedState.parameters[i], 1.e-4f);
  }
}

// First-order consistency: for a small state/linRef offset, the
// linRef-aware propagation parameter output must agree with the
// already-accepted, independently oracle-tested non-linRef propagation
// applied directly to the same starting state, since the linRef-aware
// transform is exactly a first-order (Jacobian) correction around the
// reference. This is an independent numerical cross-check: it never
// constructs or reuses this slice's own Jacobian formulas.
void checkFirstOrderConsistency(float bz, float targetZ)
{
  const auto baseState = makeState();
  auto state = baseState;
  const float epsilon = 0.01f;
  state.parameters[0] += epsilon;
  state.parameters[2] += 0.4f * epsilon;
  state.parameters[3] += 0.3f * epsilon;
  auto linRef = linRefFromState(baseState);

  auto referenceDrivenState = state;
  auto referenceDrivenLinRef = linRef;
  OperationFailureReason reason{};
  BOOST_REQUIRE(Propagator::propagateForward(referenceDrivenState, referenceDrivenLinRef, targetZ, bz, reason));

  auto selfDrivenState = state;
  BOOST_REQUIRE(Propagator::propagateForward(selfDrivenState, targetZ, bz, reason));

  const float tolerance = 5.e-3f;
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_SMALL(referenceDrivenState.parameters[i] - selfDrivenState.parameters[i], tolerance);
  }
}

BOOST_AUTO_TEST_CASE(AcceptedPropagationFirstOrderConsistency)
{
  checkFirstOrderConsistency(0.5f, -80.f);
  checkFirstOrderConsistency(-0.5f, -80.f);
  checkFirstOrderConsistency(0.f, -60.f);
}

BOOST_AUTO_TEST_CASE(AcceptedPropagationSelfDrivenReduction)
{
  checkSelfDrivenReduction(0.5f, -80.f);
  checkSelfDrivenReduction(-0.5f, -80.f);
  checkSelfDrivenReduction(0.f, -60.f);
}

BOOST_AUTO_TEST_CASE(AcceptedPropagationJacobianMatchesFiniteDifference)
{
  for (uint8_t param = 0; param < 5; ++param) {
    checkJacobianAgainstFiniteDifference(0.5f, -80.f, param);
    checkJacobianAgainstFiniteDifference(-0.5f, -80.f, param);
  }
}

// Reference-driven propagation with a perturbed linRef must differ from
// self-driven propagation of the same starting state (Jacobian evaluated
// at different points).
BOOST_AUTO_TEST_CASE(ReferenceDrivenPropagationDiffersFromSelfDrivenPropagation)
{
  const auto baseState = makeState();
  const float bz = 0.5f;
  const float targetZ = -80.f;

  auto selfDriven = baseState;
  OperationFailureReason reason{};
  BOOST_REQUIRE(Propagator::propagateForward(selfDriven, targetZ, bz, reason));

  auto referenceDrivenState = baseState;
  auto referenceDrivenRef = perturbedLinRef(baseState);
  BOOST_REQUIRE(Propagator::propagateForward(referenceDrivenState, referenceDrivenRef, targetZ, bz, reason));

  bool differs = false;
  for (uint8_t i = 0; i < 15; ++i) {
    if (std::abs(selfDriven.covariance[i] - referenceDrivenState.covariance[i]) > 1.e-6f) {
      differs = true;
      break;
    }
  }
  BOOST_CHECK(differs);
}

BOOST_AUTO_TEST_CASE(PropagateWithLinRefRejectsFamilyMismatchAndPreservesBytes)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);
  linRef.family = StateFamily::Barrel;
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!Propagator::propagateForward(state, linRef, -60.f, 0.f, reason));
  BOOST_CHECK(reason == OperationFailureReason::SourceFamilyMismatch);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

// Fitted-state/linRef pairing precondition (hardening): parameters may
// legitimately differ (the entire purpose of a linearization reference),
// but referenceCoordinate must match exactly. No alpha check for Forward
// (always 0/unused).
BOOST_AUTO_TEST_CASE(PropagateWithLinRefRejectsReferenceCoordinateMismatchAndPreservesBytes)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);
  linRef.referenceCoordinate += 1.f;
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!Propagator::propagateForward(state, linRef, -60.f, 0.f, reason));
  BOOST_CHECK(reason == OperationFailureReason::ReferenceCoordinateMismatch);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(FieldPropagationRejectsZeroCurvatureAndPreservesBytes)
{
  auto state = makeState();
  state.parameters[4] = 0.f;
  auto linRef = linRefFromState(state);
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!Propagator::propagateForward(state, linRef, -60.f, 0.5f, reason));
  BOOST_CHECK(reason == OperationFailureReason::PropagationFailure);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(PropagationRejectsZeroTanlAndPreservesBytes)
{
  auto state = makeState();
  state.parameters[3] = 0.f;
  auto linRef = linRefFromState(state);
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!Propagator::propagateForward(state, linRef, -60.f, 0.f, reason));
  BOOST_CHECK(reason == OperationFailureReason::UnreachableTarget);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(ShiftReferenceToMeasurementSetsOnlyXAndY)
{
  auto linRef = linRefFromState(makeState());
  const auto before = linRef;
  SurfaceMeasurement measurement{};
  measurement.frame.u = 12.3f;
  measurement.frame.v = -45.6f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(shiftReferenceToMeasurement(linRef, measurement, reason));
  BOOST_CHECK_EQUAL(linRef.parameters[0], measurement.frame.u);
  BOOST_CHECK_EQUAL(linRef.parameters[1], measurement.frame.v);
  BOOST_CHECK_EQUAL(linRef.parameters[2], before.parameters[2]);
  BOOST_CHECK_EQUAL(linRef.parameters[3], before.parameters[3]);
  BOOST_CHECK_EQUAL(linRef.parameters[4], before.parameters[4]);
  BOOST_CHECK_EQUAL(linRef.referenceCoordinate, before.referenceCoordinate);
}

BOOST_AUTO_TEST_CASE(ShiftReferenceToMeasurementRejectsNonFiniteAndPreservesBytes)
{
  auto linRef = linRefFromState(makeState());
  const auto before = linRef;
  SurfaceMeasurement measurement{};
  measurement.frame.u = 1.f;
  measurement.frame.v = std::numeric_limits<float>::infinity();
  OperationFailureReason reason{};
  BOOST_CHECK(!shiftReferenceToMeasurement(linRef, measurement, reason));
  BOOST_CHECK(bitEqual(linRef, before));
}

// buildSeed fixtures, matching testForwardSurfaceStateOperations.cxx's own
// conventions (z-ordered, well-separated hits).
SurfaceMeasurement makeInnerMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.frame = {-46.f, 0.5f, 0.3f, 0.f};
  measurement.frame.q = -46.f;
  measurement.covariance = {0.02f, 0.001f, 0.03f};
  return measurement;
}

SurfaceMeasurement makeMiddleMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.frame = {-65.f, 1.1f, 0.9f, 0.f};
  measurement.frame.q = -65.f;
  return measurement;
}

SurfaceMeasurement makeOuterMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.frame = {-90.f, 2.0f, 1.8f, 0.f};
  measurement.frame.q = -90.f;
  measurement.covariance = {0.04f, 0.002f, 0.05f};
  return measurement;
}

// A second transverse handedness: negates the y-component of all three
// global positions (z-ordering, hence the SeedGeometryDegenerate boundary,
// is untouched), giving a geometrically distinct candidate whose estimated
// signed q/pT has the opposite sense from the un-mirrored fixture above.
SurfaceMeasurement mirrorTransverse(SurfaceMeasurement measurement)
{
  measurement.frame.v = -measurement.frame.v;
  return measurement;
}

// The mirrored fixture is a geometrically distinct candidate (opposite
// transverse handedness): its outer-seed phi (the atan2(dyPhi, dxPhi)
// direction estimate, sign-sensitive to handedness) must differ from the
// un-mirrored fixture's, confirming the two handedness fixtures above are
// not accidentally equivalent. parameters[4] (invQPt) is deliberately not
// used as the discriminator here: forward::buildSeed's established
// `(trackletMinPt > 0.f) ? 1.f/trackletMinPt : 0.f` compatibility fallback
// (see the header doc) makes it a fixed configured magnitude, not a
// geometry-derived estimate, so it is identical for both fixtures by
// construction.
BOOST_AUTO_TEST_CASE(MirroredFixtureIsGeometricallyDistinctFromUnmirrored)
{
  SurfaceKinematicState plain{};
  OperationFailureReason reasonPlain{};
  BOOST_REQUIRE(forward::buildSeed(makeInnerMeasurement(), makeMiddleMeasurement(), makeOuterMeasurement(), 0.5f, 0.3f,
                                   1, o2::track::PID::Pion, plain, reasonPlain));

  SurfaceKinematicState mirrored{};
  OperationFailureReason reasonMirrored{};
  BOOST_REQUIRE(forward::buildSeed(mirrorTransverse(makeInnerMeasurement()), mirrorTransverse(makeMiddleMeasurement()),
                                   mirrorTransverse(makeOuterMeasurement()), 0.5f, 0.3f, 1, o2::track::PID::Pion,
                                   mirrored, reasonMirrored));

  BOOST_CHECK_NE(plain.parameters[2], mirrored.parameters[2]);
}

// Reports the observed maximum absolute/relative drift accumulated by the
// item 4 numerical cross-checks above. Relies on Boost.Test's default
// registration-order execution within one translation unit to run last.
BOOST_AUTO_TEST_CASE(ReportNonlinearForwardNumericalDrift)
{
  BOOST_TEST_MESSAGE("self-driven reduction max drift: parameters abs="
                     << selfDrivenReductionParameterDrift.maximumAbsolute
                     << " rel=" << selfDrivenReductionParameterDrift.maximumRelative
                     << ", covariance abs=" << selfDrivenReductionCovarianceDrift.maximumAbsolute
                     << " rel=" << selfDrivenReductionCovarianceDrift.maximumRelative);
  BOOST_TEST_MESSAGE("finite-difference Jacobian max drift: abs=" << jacobianFiniteDifferenceDrift.maximumAbsolute
                                                                  << " rel=" << jacobianFiniteDifferenceDrift.maximumRelative);
}
