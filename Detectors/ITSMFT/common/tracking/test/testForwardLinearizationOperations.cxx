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
#include "ITSMFTTracking/SurfaceLinearizationReference.h"
#include "ITSMFTTracking/TransitionPolicy.h"

namespace
{
using namespace o2::itsmft::tracking;
using namespace o2::itsmft::tracking::forward;

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
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

} // namespace

BOOST_AUTO_TEST_CASE(PropagateLinearMatchesIndependentAnalyticOracle)
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
  BOOST_REQUIRE((propagate<PropagationModel::Linear>(state, linRef, -80.f, 0.f, reason)));

  BOOST_CHECK_CLOSE(state.referenceCoordinate, expectedState.referenceCoordinate, 1.e-3f);
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_SMALL(state.parameters[i] - expectedState.parameters[i], 1.e-4f);
    BOOST_CHECK_SMALL(linRef.parameters[i] - expectedRef.parameters[i], 1.e-4f);
  }
}

BOOST_AUTO_TEST_CASE(PropagateLinearMatchesOracleZeroField)
{
  const auto baseState = makeState();
  auto state = baseState;
  state.parameters[1] -= 0.03f;
  auto linRef = linRefFromState(baseState);

  SurfaceKinematicState expectedState{};
  SurfaceLinearizationReference expectedRef{};
  computeExpectedLinear(state, linRef, -60.f, expectedState, expectedRef);

  OperationFailureReason reason{};
  BOOST_REQUIRE((propagate<PropagationModel::Linear>(state, linRef, -60.f, 0.f, reason)));
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_SMALL(state.parameters[i] - expectedState.parameters[i], 1.e-4f);
  }
}

// First-order consistency: for a small state/linRef offset, the
// linRef-aware propagate<Model> parameter output must agree with the
// already-accepted, independently oracle-tested non-linRef propagate<Model>
// applied directly to the same starting state, since the linRef-aware
// transform is exactly a first-order (Jacobian) correction around the
// reference. This is an independent numerical cross-check: it never
// constructs or reuses this slice's own Jacobian formulas.
template <PropagationModel Model>
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
  BOOST_REQUIRE((propagate<Model>(referenceDrivenState, referenceDrivenLinRef, targetZ, bz, reason)));

  auto selfDrivenState = state;
  BOOST_REQUIRE((propagate<Model>(selfDrivenState, targetZ, bz, reason)));

  const float tolerance = 5.e-3f;
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_SMALL(referenceDrivenState.parameters[i] - selfDrivenState.parameters[i], tolerance);
  }
}

BOOST_AUTO_TEST_CASE(PropagateQuadraticFirstOrderConsistencyNonzeroField)
{
  checkFirstOrderConsistency<PropagationModel::Quadratic>(0.5f, -80.f);
}

BOOST_AUTO_TEST_CASE(PropagateHelixFirstOrderConsistencyNonzeroField)
{
  checkFirstOrderConsistency<PropagationModel::Helix>(0.5f, -80.f);
}

BOOST_AUTO_TEST_CASE(PropagateOptimizedFirstOrderConsistencyNonzeroField)
{
  checkFirstOrderConsistency<PropagationModel::Optimized>(0.5f, -80.f);
}

BOOST_AUTO_TEST_CASE(PropagateOptimizedFirstOrderConsistencyZeroField)
{
  checkFirstOrderConsistency<PropagationModel::Optimized>(0.f, -60.f);
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
  BOOST_REQUIRE((propagate<PropagationModel::Helix>(selfDriven, targetZ, bz, reason)));

  auto referenceDrivenState = baseState;
  auto referenceDrivenRef = perturbedLinRef(baseState);
  BOOST_REQUIRE((propagate<PropagationModel::Helix>(referenceDrivenState, referenceDrivenRef, targetZ, bz, reason)));

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
  BOOST_CHECK(!(propagate<PropagationModel::Linear>(state, linRef, -60.f, 0.f, reason)));
  BOOST_CHECK(reason == OperationFailureReason::SourceFamilyMismatch);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(PropagateHelixRejectsZeroFieldAndPreservesBytes)
{
  auto state = makeState();
  auto linRef = linRefFromState(state);
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!(propagate<PropagationModel::Helix>(state, linRef, -60.f, 0.f, reason)));
  BOOST_CHECK(reason == OperationFailureReason::PropagationFailure);
  BOOST_CHECK(bitEqual(state, stateBefore));
  BOOST_CHECK(bitEqual(linRef, linRefBefore));
}

BOOST_AUTO_TEST_CASE(PropagateLinearRejectsZeroTanlAndPreservesBytes)
{
  auto state = makeState();
  state.parameters[3] = 0.f;
  auto linRef = linRefFromState(state);
  const auto stateBefore = state;
  const auto linRefBefore = linRef;
  OperationFailureReason reason{};
  BOOST_CHECK(!(propagate<PropagationModel::Linear>(state, linRef, -60.f, 0.f, reason)));
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
  measurement.global = {0.5f, 0.3f, -46.f};
  measurement.frame.q = -46.f;
  measurement.covariance = {0.02f, 0.001f, 0.03f};
  return measurement;
}

SurfaceMeasurement makeMiddleMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.global = {1.1f, 0.9f, -65.f};
  measurement.frame.q = -65.f;
  return measurement;
}

SurfaceMeasurement makeOuterMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.global = {2.0f, 1.8f, -90.f};
  measurement.frame.q = -90.f;
  measurement.covariance = {0.04f, 0.002f, 0.05f};
  return measurement;
}

BOOST_AUTO_TEST_CASE(AnchoredBuildSeedOuterIsByteCompatibleWithExistingBuildSeed)
{
  const auto inner = makeInnerMeasurement();
  const auto middle = makeMiddleMeasurement();
  const auto outer = makeOuterMeasurement();
  const float bz = 0.5f;
  const float trackletMinPt = 0.3f;

  SurfaceKinematicState plain{};
  OperationFailureReason reasonPlain{};
  BOOST_REQUIRE(forward::buildSeed(inner, middle, outer, bz, trackletMinPt, 1, o2::track::PID::Pion, plain, reasonPlain));

  SurfaceKinematicState anchored{};
  OperationFailureReason reasonAnchored{};
  BOOST_REQUIRE(forward::buildSeed(SeedAnchor::Outer, inner, middle, outer, bz, trackletMinPt, 1, o2::track::PID::Pion, anchored, reasonAnchored));

  BOOST_CHECK(bitEqual(plain, anchored));
}

// SeedAnchor::Inner: same physical hits, anchored at the inner measurement's
// own frame/reference/covariance instead of the outer's, with the identical
// (anchor-symmetric) direction estimate as Outer -- correct signed q/pT and
// endpoint frame per the header contract.
BOOST_AUTO_TEST_CASE(AnchoredBuildSeedInnerUsesInnerFrameWithSameDirectionAsOuter)
{
  const auto inner = makeInnerMeasurement();
  const auto middle = makeMiddleMeasurement();
  const auto outer = makeOuterMeasurement();
  const float bz = 0.5f;
  const float trackletMinPt = 0.3f;

  SurfaceKinematicState outerAnchored{};
  OperationFailureReason reasonOuter{};
  BOOST_REQUIRE(forward::buildSeed(SeedAnchor::Outer, inner, middle, outer, bz, trackletMinPt, 1, o2::track::PID::Pion, outerAnchored, reasonOuter));

  SurfaceKinematicState innerAnchored{};
  OperationFailureReason reasonInner{};
  BOOST_REQUIRE(forward::buildSeed(SeedAnchor::Inner, inner, middle, outer, bz, trackletMinPt, 1, o2::track::PID::Pion, innerAnchored, reasonInner));

  BOOST_CHECK(innerAnchored.family == StateFamily::Forward);
  BOOST_CHECK_EQUAL(innerAnchored.referenceCoordinate, inner.frame.q);
  BOOST_CHECK_EQUAL(innerAnchored.parameters[0], inner.global.x);
  BOOST_CHECK_EQUAL(innerAnchored.parameters[1], inner.global.y);
  BOOST_CHECK_EQUAL(innerAnchored.covariance[packedCovarianceIndex(0, 0)], inner.covariance.uu);
  BOOST_CHECK_EQUAL(innerAnchored.covariance[packedCovarianceIndex(1, 1)], inner.covariance.vv);

  // Same direction estimate (phi/tanl/invQPt -- the signed q/pT) as Outer,
  // by the documented anchor-symmetry of the closed-form geometry estimate.
  BOOST_CHECK_CLOSE(innerAnchored.parameters[2], outerAnchored.parameters[2], 1.e-3f);
  BOOST_CHECK_CLOSE(innerAnchored.parameters[3], outerAnchored.parameters[3], 1.e-3f);
  BOOST_CHECK_CLOSE(innerAnchored.parameters[4], outerAnchored.parameters[4], 1.e-3f);
}

BOOST_AUTO_TEST_CASE(AnchoredBuildSeedRejectsInvalidAnchorTransactionally)
{
  const auto inner = makeInnerMeasurement();
  const auto middle = makeMiddleMeasurement();
  const auto outer = makeOuterMeasurement();
  auto outState = makeState();
  const auto before = outState;
  OperationFailureReason reason{};
  BOOST_CHECK(!forward::buildSeed(static_cast<SeedAnchor>(2), inner, middle, outer, 0.5f, 0.3f, 1, o2::track::PID::Pion, outState, reason));
  BOOST_CHECK(bitEqual(outState, before));
}
