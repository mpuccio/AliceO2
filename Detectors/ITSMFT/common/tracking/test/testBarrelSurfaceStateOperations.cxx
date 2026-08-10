// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFTBarrelSurfaceStateOperations
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
#include "ITStracking/Cluster.h"
// buildTrackSeed is the retained legacy oracle for buildSeed's test coverage
// only (see BuildSeed* test cases below) -- production never calls it or
// constructs the TrackParCov it returns (BarrelSurfaceStateOperations.cxx).
#include "ITStracking/TrackHelpers.h"

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
  state.family = StateFamily::Barrel;
  state.flags = 0x5a;
  state.absCharge = 1;
  state.pid = o2::track::PID::Kaon;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
    }
  }
  return state;
}

SurfaceKinematicState makeCandidateState()
{
  auto state = makeState();
  state.parameters[0] += 0.05f;
  state.parameters[1] -= 0.03f;
  state.parameters[2] += 0.01f;
  state.parameters[3] += 0.02f;
  state.parameters[4] -= 0.015f;
  state.flags = 0x3c; // Deliberately different non-family metadata: stateChi2 must not care.
  state.absCharge = 2;
  state.pid = o2::track::PID::Pion;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.02f * (row + 1) : 0.0001f * (row + column + 1);
    }
  }
  return state;
}

SurfaceMeasurement makeMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.frame = {456.f, 0.8f, -0.45f, 0.f};
  measurement.covariance = {0.04f, 0.012f, 0.09f};
  return measurement;
}

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

// buildSeed fixtures: three well-separated, non-collinear points so
// computeCurvature/computeCurvatureCentreX/computeTanDipAngle stay away from
// their own degenerate fallbacks in the "success" oracle tests below.
// "Mirror" negates the transverse offsets to flip the sign of the estimated
// curvature/q2pt, matching an opposite-sign candidate at the same field.
o2::its::Cluster makeInnerCluster(bool mirror = false)
{
  const float sign = mirror ? -1.f : 1.f;
  return o2::its::Cluster{1.8f, sign * 0.10f, -0.60f, 0};
}

o2::its::Cluster makeMiddleCluster(bool mirror = false)
{
  const float sign = mirror ? -1.f : 1.f;
  return o2::its::Cluster{3.0f, sign * 0.35f, 0.10f, 1};
}

o2::its::TrackingFrameInfo makeOuterHit(bool mirror = false)
{
  // xCoordinate/yCoordinate/zCoordinate are irrelevant to buildTrackSeed's
  // formula (it reads only xTrackingFrame/alphaTrackingFrame/
  // positionTrackingFrame/covarianceTrackingFrame) and are set to 0 here.
  const float sign = mirror ? -1.f : 1.f;
  return o2::its::TrackingFrameInfo{0.f, 0.f, 0.f, 4.4f, 0.15f, {sign * 0.5f, -0.8f}, {0.02f, 0.001f, 0.03f}};
}

void checkBuildSeedMetadata(const SurfaceKinematicState& state, uint8_t expectedAbsCharge, o2::track::PID expectedPid)
{
  BOOST_CHECK(state.family == StateFamily::Barrel);
  BOOST_CHECK_EQUAL(state.flags, uint8_t{0});
  BOOST_CHECK_EQUAL(state.absCharge, expectedAbsCharge);
  BOOST_CHECK_EQUAL(static_cast<uint8_t>(state.pid), static_cast<uint8_t>(expectedPid));
}

GlobalPoint3F globalPointFromCluster(const o2::its::Cluster& cluster)
{
  return {cluster.xCoordinate, cluster.yCoordinate, cluster.zCoordinate};
}

// Test-local field-mapping helper: builds the SurfaceMeasurement fields
// barrel::buildSeed reads from the outer tracking-frame hit (Cylinder field
// mapping: reference coordinate -> measurement.frame.q, frame angle ->
// measurement.frame.frameAngle, measured local coordinates ->
// measurement.frame.u/v, measured covariance -> measurement.covariance).
SurfaceMeasurement measurementFromOuterHit(const o2::its::TrackingFrameInfo& hit)
{
  SurfaceMeasurement measurement{};
  measurement.frame.q = hit.xTrackingFrame;
  measurement.frame.frameAngle = hit.alphaTrackingFrame;
  measurement.frame.u = hit.positionTrackingFrame[0];
  measurement.frame.v = hit.positionTrackingFrame[1];
  measurement.covariance.uu = hit.covarianceTrackingFrame[0];
  measurement.covariance.uv = hit.covarianceTrackingFrame[1];
  measurement.covariance.vv = hit.covarianceTrackingFrame[2];
  return measurement;
}

void checkBuildSeedFailurePreservesBytes(const o2::its::Cluster& clusterInner, const o2::its::Cluster& clusterMiddle,
                                         const o2::its::TrackingFrameInfo& hitOuter, float bz,
                                         OperationFailureReason expected)
{
  auto outState = makeState(); // deliberately non-default sentinel pattern
  const auto before = outState;
  OperationFailureReason reason{};
  BOOST_CHECK(!barrel::buildSeed(globalPointFromCluster(clusterInner), globalPointFromCluster(clusterMiddle),
                                 measurementFromOuterHit(hitOuter), bz, 1, o2::track::PID::Pion, outState, reason));
  BOOST_CHECK(reason == expected);
  BOOST_CHECK(bitEqual(outState, before));
}

void checkStateChi2FailurePreservesBytes(SurfaceKinematicState reference, SurfaceKinematicState candidate, float chi2,
                                         OperationFailureReason expected)
{
  const auto referenceBefore = reference;
  const auto candidateBefore = candidate;
  const float chi2Before = chi2;
  OperationFailureReason reason{};
  BOOST_CHECK(!stateChi2(reference, candidate, chi2, reason));
  BOOST_CHECK_EQUAL(static_cast<uint8_t>(reason), static_cast<uint8_t>(expected));
  BOOST_CHECK(bitEqual(reference, referenceBefore));
  BOOST_CHECK(bitEqual(candidate, candidateBefore));
  BOOST_CHECK_EQUAL(chi2, chi2Before);
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

void checkMetadata(const SurfaceKinematicState& state, const SurfaceKinematicState& before)
{
  BOOST_CHECK_EQUAL(static_cast<uint8_t>(state.family), static_cast<uint8_t>(before.family));
  BOOST_CHECK_EQUAL(state.flags, before.flags);
  BOOST_CHECK_EQUAL(state.absCharge, before.absCharge);
  BOOST_CHECK_EQUAL(static_cast<uint8_t>(state.pid), static_cast<uint8_t>(before.pid));
}

void compareWithOracle(const SurfaceKinematicState& state, const SurfaceKinematicState& before,
                       const o2::track::TrackParCovF& oracle, Drift& drift)
{
  checkClose(state.referenceCoordinate, oracle.getX(), drift);
  checkClose(state.alpha, oracle.getAlpha(), drift);
  for (uint8_t i = 0; i < 5; ++i) {
    checkClose(state.parameters[i], oracle.getParam(i), drift);
  }
  for (uint8_t i = 0; i < 15; ++i) {
    checkClose(state.covariance[i], oracle.getCov()[i], drift);
  }
  checkMetadata(state, before);
}

float directPredictedChi2(const SurfaceKinematicState& state, const SurfaceMeasurement& measurement)
{
  const float s00 = state.covariance[packedCovarianceIndex(0, 0)] + measurement.covariance.uu;
  const float s01 = state.covariance[packedCovarianceIndex(1, 0)] + measurement.covariance.uv;
  const float s11 = state.covariance[packedCovarianceIndex(1, 1)] + measurement.covariance.vv;
  const float determinant = s00 * s11 - s01 * s01;
  const float residualY = measurement.frame.u - state.parameters[0];
  const float residualZ = measurement.frame.v - state.parameters[1];
  return (s11 * residualY * residualY - 2.f * s01 * residualY * residualZ + s00 * residualZ * residualZ) / determinant;
}

void directUpdate(const SurfaceKinematicState& before, const SurfaceMeasurement& measurement,
                  SurfaceKinematicState& expected, float& expectedChi2)
{
  expected = before;
  const float s00 = before.covariance[0] + measurement.covariance.uu;
  const float s01 = before.covariance[1] + measurement.covariance.uv;
  const float s11 = before.covariance[2] + measurement.covariance.vv;
  const float determinant = s00 * s11 - s01 * s01;
  const float inverse00 = s11 / determinant;
  const float inverse01 = -s01 / determinant;
  const float inverse11 = s00 / determinant;
  float covariance[5][5]{};
  float updated[5][5]{};
  const float residual[2] = {measurement.frame.u - before.parameters[0], measurement.frame.v - before.parameters[1]};
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column < 5; ++column) {
      covariance[row][column] = before.covariance[packedCovarianceIndex(row, column)];
    }
  }
  for (uint8_t row = 0; row < 5; ++row) {
    const float gain0 = covariance[row][0] * inverse00 + covariance[row][1] * inverse01;
    const float gain1 = covariance[row][0] * inverse01 + covariance[row][1] * inverse11;
    expected.parameters[row] += gain0 * residual[0] + gain1 * residual[1];
    for (uint8_t column = 0; column < 5; ++column) {
      updated[row][column] = covariance[row][column] - gain0 * covariance[0][column] - gain1 * covariance[1][column];
    }
  }
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      expected.covariance[packedCovarianceIndex(row, column)] = updated[row][column];
    }
  }
  expectedChi2 = directPredictedChi2(before, measurement);
}

template <typename Callable>
void checkStateFailurePreservesBytes(SurfaceKinematicState state, float chi2, OperationFailureReason expected, Callable&& operation)
{
  const auto before = state;
  const float beforeChi2 = chi2;
  OperationFailureReason reason{};
  BOOST_CHECK(!operation(state, chi2, reason));
  BOOST_CHECK_EQUAL(static_cast<uint8_t>(reason), static_cast<uint8_t>(expected));
  BOOST_CHECK(bitEqual(state, before));
  BOOST_CHECK_EQUAL(chi2, beforeChi2);
}

} // namespace

BOOST_AUTO_TEST_CASE(RotationMatchesCompleteRetainedLegacyOracleAndCanonicalizesAlpha)
{
  auto state = makeState();
  const auto before = state;
  o2::track::TrackParCovF oracle{};
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(state, oracle));
  const float targetAlpha = 4.f * o2::constants::math::PI + 0.55f;
  BOOST_REQUIRE(oracle.rotate(targetAlpha));
  OperationFailureReason reason{};
  BOOST_REQUIRE(rotate(state, targetAlpha, reason));
  BOOST_CHECK_GE(state.alpha, -o2::constants::math::PI);
  BOOST_CHECK_LE(state.alpha, o2::constants::math::PI);
  Drift drift{};
  compareWithOracle(state, before, oracle, drift);
  BOOST_TEST_MESSAGE("rotate retained-oracle max drift: abs=" << drift.absolute << " rel=" << drift.relative);
}

BOOST_AUTO_TEST_CASE(PropagationMatchesCompleteRetainedLegacyOracleAtBoundaries)
{
  for (const float targetX : {2.5f, 4.f, 6.5f}) { // inward, no-op, outward
    auto state = makeState();
    const auto before = state;
    o2::track::TrackParCovF oracle{};
    BOOST_REQUIRE(legacy::exportBarrelTrackParCov(state, oracle));
    BOOST_REQUIRE(oracle.propagateTo(targetX, 5.f));
    OperationFailureReason reason{};
    BOOST_REQUIRE(propagate(state, targetX, 5.f, reason));
    Drift drift{};
    compareWithOracle(state, before, oracle, drift);
    BOOST_TEST_MESSAGE("propagate retained-oracle max drift: abs=" << drift.absolute << " rel=" << drift.relative);
  }
  auto zeroField = makeState();
  auto oracle = o2::track::TrackParCovF{};
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(zeroField, oracle));
  BOOST_REQUIRE(oracle.propagateTo(6.5f, 0.f));
  OperationFailureReason reason{};
  BOOST_REQUIRE(propagate(zeroField, 6.5f, 0.f, reason));
  Drift drift{};
  compareWithOracle(zeroField, makeState(), oracle, drift);
}

BOOST_AUTO_TEST_CASE(PredictedChi2AndUpdateMatchDirectAndRetainedOracles)
{
  const auto before = makeState();
  const auto measurement = makeMeasurement();
  float chi2 = -1.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(predictedChi2(before, measurement, chi2, reason));
  BOOST_CHECK_CLOSE_FRACTION(chi2, directPredictedChi2(before, measurement), 2.e-6f);

  o2::track::TrackParCovF oracle{};
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(before, oracle));
  const o2::track::TrackParCovF::dim2_t position{measurement.frame.u, measurement.frame.v};
  const o2::track::TrackParCovF::dim3_t covariance{measurement.covariance.uu, measurement.covariance.uv, measurement.covariance.vv};
  const float oracleChi2 = oracle.getPredictedChi2Quiet(position, covariance);
  BOOST_CHECK_CLOSE_FRACTION(chi2, oracleChi2, 2.e-5f);
  BOOST_REQUIRE(oracle.update(position, covariance));

  SurfaceKinematicState expected{};
  float expectedChi2 = 0.f;
  directUpdate(before, measurement, expected, expectedChi2);
  auto state = before;
  BOOST_REQUIRE(update(state, measurement, chi2, reason));
  BOOST_CHECK_CLOSE_FRACTION(chi2, expectedChi2, 2.e-6f);
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_CLOSE_FRACTION(state.parameters[i], expected.parameters[i], 2.e-6f);
  }
  for (uint8_t i = 0; i < 15; ++i) {
    BOOST_CHECK_CLOSE_FRACTION(state.covariance[i], expected.covariance[i], 3.e-6f);
  }
  Drift drift{};
  compareWithOracle(state, before, oracle, drift);
  BOOST_TEST_MESSAGE("update retained-oracle max drift: abs=" << drift.absolute << " rel=" << drift.relative);
}

BOOST_AUTO_TEST_CASE(FailuresAreTransactionalAndAnalyticallyGuaranteed)
{
  constexpr float InitialChi2 = 123.f;
  auto wrongFamily = makeState();
  wrongFamily.family = StateFamily::Forward;
  checkStateFailurePreservesBytes(wrongFamily, InitialChi2, OperationFailureReason::SourceFamilyMismatch,
                                  [](auto& state, auto&, auto& reason) { return rotate(state, 0.5f, reason); });
  auto nonFiniteState = makeState();
  nonFiniteState.parameters[0] = std::numeric_limits<float>::infinity();
  checkStateFailurePreservesBytes(nonFiniteState, InitialChi2, OperationFailureReason::NonFiniteInput,
                                  [](auto& state, auto&, auto& reason) { return propagate(state, 5.f, 0.f, reason); });
  auto nonFiniteTarget = makeState();
  checkStateFailurePreservesBytes(nonFiniteTarget, InitialChi2, OperationFailureReason::NonFiniteInput,
                                  [](auto& state, auto&, auto& reason) { return rotate(state, std::numeric_limits<float>::infinity(), reason); });
  auto nonFiniteBz = makeState();
  checkStateFailurePreservesBytes(nonFiniteBz, InitialChi2, OperationFailureReason::NonFiniteInput,
                                  [](auto& state, auto&, auto& reason) { return propagate(state, 5.f, std::numeric_limits<float>::infinity(), reason); });

  auto rotationFailure = makeState();
  const float delta = o2::constants::math::PI;
  const float csp = std::sqrt(1.f - rotationFailure.parameters[2] * rotationFailure.parameters[2]);
  BOOST_REQUIRE_LT(csp * std::cos(delta) + rotationFailure.parameters[2] * std::sin(delta), 0.f);
  checkStateFailurePreservesBytes(rotationFailure, InitialChi2, OperationFailureReason::RotationFailure,
                                  [delta](auto& state, auto&, auto& reason) { return rotate(state, state.alpha + delta, reason); });

  auto propagationFailure = makeState();
  constexpr float bz = 5.f;
  const float targetX = propagationFailure.referenceCoordinate +
                        (1.1f - propagationFailure.parameters[2]) / (propagationFailure.parameters[4] * bz * o2::constants::math::B2C);
  const float propagatedSnp = propagationFailure.parameters[2] + propagationFailure.parameters[4] * bz * o2::constants::math::B2C *
                                                                   (targetX - propagationFailure.referenceCoordinate);
  BOOST_REQUIRE_GE(std::abs(propagatedSnp), 1.f);
  checkStateFailurePreservesBytes(propagationFailure, InitialChi2, OperationFailureReason::UnreachableTarget,
                                  [targetX](auto& state, auto&, auto& reason) { return propagate(state, targetX, bz, reason); });

  auto singular = makeState();
  auto measurement = makeMeasurement();
  measurement.covariance = {-singular.covariance[0], -singular.covariance[1], -singular.covariance[2]};
  float chi2 = InitialChi2;
  const auto before = singular;
  OperationFailureReason reason{};
  BOOST_CHECK(!predictedChi2(singular, measurement, chi2, reason));
  BOOST_CHECK_EQUAL(static_cast<uint8_t>(reason), static_cast<uint8_t>(OperationFailureReason::InvalidCovariance));
  BOOST_CHECK_EQUAL(chi2, InitialChi2);
  BOOST_CHECK(bitEqual(singular, before));
  BOOST_CHECK(!update(singular, measurement, chi2, reason));
  BOOST_CHECK_EQUAL(static_cast<uint8_t>(reason), static_cast<uint8_t>(OperationFailureReason::InvalidCovariance));
  BOOST_CHECK_EQUAL(chi2, InitialChi2);
  BOOST_CHECK(bitEqual(singular, before));

  auto nonFiniteOutput = makeState();
  nonFiniteOutput.referenceCoordinate = std::numeric_limits<float>::max();
  nonFiniteOutput.parameters[0] = std::numeric_limits<float>::max();
  checkStateFailurePreservesBytes(nonFiniteOutput, InitialChi2, OperationFailureReason::NonFiniteOutput,
                                  [](auto& state, auto&, auto& reason) { return rotate(state, state.alpha + o2::constants::math::PIQuarter, reason); });
}

BOOST_AUTO_TEST_CASE(RepeatedOperationChainIsByteDeterministic)
{
  const auto measurement = makeMeasurement();
  auto first = makeState();
  auto second = first;
  for (auto* state : {&first, &second}) {
    OperationFailureReason reason{};
    float chi2 = 0.f;
    BOOST_REQUIRE(rotate(*state, 0.55f, reason));
    BOOST_REQUIRE(propagate(*state, 5.5f, 5.f, reason));
    BOOST_REQUIRE(update(*state, measurement, chi2, reason));
  }
  BOOST_CHECK(bitEqual(first, second));
}

BOOST_AUTO_TEST_CASE(StateChi2MatchesRetainedLegacyOracleAndHandReferences)
{
  const auto reference = makeState();
  const auto candidate = makeCandidateState();
  float chi2 = -1.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(stateChi2(reference, candidate, chi2, reason));

  o2::track::TrackParCovF referenceOracle{};
  o2::track::TrackParCovF candidateOracle{};
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(reference, referenceOracle));
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(candidate, candidateOracle));
  // Retained-oracle drift characterization only: TrackParCovF::getPredictedChi2
  // inverts the combined covariance in double precision and does not prove
  // physical correctness of this float-native primitive.
  const float oracleChi2 = referenceOracle.getPredictedChi2(candidateOracle);
  Drift drift{};
  checkClose(chi2, oracleChi2, drift);
  BOOST_TEST_MESSAGE("stateChi2 retained-oracle max drift: abs=" << drift.absolute << " rel=" << drift.relative);

  // Independent hand reference #1: diagonal combined covariance, split evenly
  // between reference and candidate to also exercise the packed-storage sum.
  auto diagonalReference = makeState();
  auto diagonalCandidate = makeCandidateState();
  for (float& value : diagonalReference.covariance) {
    value = 0.f;
  }
  for (float& value : diagonalCandidate.covariance) {
    value = 0.f;
  }
  const float diagShare[5] = {0.02f, 0.015f, 0.05f, 0.03f, 0.04f};
  for (uint8_t i = 0; i < 5; ++i) {
    diagonalReference.covariance[packedCovarianceIndex(i, i)] = diagShare[i];
    diagonalCandidate.covariance[packedCovarianceIndex(i, i)] = diagShare[i];
  }
  float diagonalChi2 = -1.f;
  BOOST_REQUIRE(stateChi2(diagonalReference, diagonalCandidate, diagonalChi2, reason));
  float diff[5];
  for (uint8_t i = 0; i < 5; ++i) {
    diff[i] = diagonalReference.parameters[i] - diagonalCandidate.parameters[i];
  }
  float handDiagonalChi2 = 0.f;
  for (uint8_t i = 0; i < 5; ++i) {
    handDiagonalChi2 += diff[i] * diff[i] / (2.f * diagShare[i]);
  }
  BOOST_CHECK_CLOSE_FRACTION(diagonalChi2, handDiagonalChi2, 2.e-5f);

  // Independent hand reference #2: rows (0,1) correlated, rows 2/3/4 diagonal;
  // the combined covariance is block-diagonal so a closed-form 2x2 inverse
  // plus per-row division verifies the packed symmetric access independently
  // of the production Bunch-Kaufman inversion.
  auto blockReference = makeState();
  auto blockCandidate = makeCandidateState();
  for (float& value : blockReference.covariance) {
    value = 0.f;
  }
  for (float& value : blockCandidate.covariance) {
    value = 0.f;
  }
  constexpr float d0 = 0.05f, d1 = 0.04f, c01 = 0.012f, d2 = 0.02f, d3 = 0.03f, d4 = 0.025f;
  blockReference.covariance[packedCovarianceIndex(0, 0)] = 0.6f * d0;
  blockCandidate.covariance[packedCovarianceIndex(0, 0)] = 0.4f * d0;
  blockReference.covariance[packedCovarianceIndex(1, 1)] = 0.6f * d1;
  blockCandidate.covariance[packedCovarianceIndex(1, 1)] = 0.4f * d1;
  blockReference.covariance[packedCovarianceIndex(1, 0)] = 0.6f * c01;
  blockCandidate.covariance[packedCovarianceIndex(1, 0)] = 0.4f * c01;
  blockReference.covariance[packedCovarianceIndex(2, 2)] = 0.6f * d2;
  blockCandidate.covariance[packedCovarianceIndex(2, 2)] = 0.4f * d2;
  blockReference.covariance[packedCovarianceIndex(3, 3)] = 0.6f * d3;
  blockCandidate.covariance[packedCovarianceIndex(3, 3)] = 0.4f * d3;
  blockReference.covariance[packedCovarianceIndex(4, 4)] = 0.6f * d4;
  blockCandidate.covariance[packedCovarianceIndex(4, 4)] = 0.4f * d4;
  float blockChi2 = -1.f;
  BOOST_REQUIRE(stateChi2(blockReference, blockCandidate, blockChi2, reason));
  for (uint8_t i = 0; i < 5; ++i) {
    diff[i] = blockReference.parameters[i] - blockCandidate.parameters[i];
  }
  const float det = d0 * d1 - c01 * c01;
  const float inverse00 = d1 / det;
  const float inverse01 = -c01 / det;
  const float inverse11 = d0 / det;
  const float handBlockChi2 = diff[0] * (inverse00 * diff[0] + inverse01 * diff[1]) +
                              diff[1] * (inverse01 * diff[0] + inverse11 * diff[1]) +
                              diff[2] * diff[2] / d2 + diff[3] * diff[3] / d3 + diff[4] * diff[4] / d4;
  BOOST_CHECK_CLOSE_FRACTION(blockChi2, handBlockChi2, 3.e-5f);

  // Ill-conditioned but invertible characterization: a near-perfectly
  // correlated (0,1) block (determinant close to, but not exactly, zero).
  auto illConditionedReference = makeState();
  auto illConditionedCandidate = makeCandidateState();
  for (float& value : illConditionedCandidate.covariance) {
    value = 0.f;
  }
  constexpr float illD0 = 1.f, illD1 = 1.f, illC01 = 0.999999f;
  for (float& value : illConditionedReference.covariance) {
    value = 0.f;
  }
  illConditionedReference.covariance[packedCovarianceIndex(0, 0)] = illD0;
  illConditionedReference.covariance[packedCovarianceIndex(1, 1)] = illD1;
  illConditionedReference.covariance[packedCovarianceIndex(1, 0)] = illC01;
  illConditionedReference.covariance[packedCovarianceIndex(2, 2)] = 0.02f;
  illConditionedReference.covariance[packedCovarianceIndex(3, 3)] = 0.03f;
  illConditionedReference.covariance[packedCovarianceIndex(4, 4)] = 0.04f;
  float illConditionedChi2 = -1.f;
  BOOST_REQUIRE(stateChi2(illConditionedReference, illConditionedCandidate, illConditionedChi2, reason));
  BOOST_CHECK(std::isfinite(illConditionedChi2));
  for (uint8_t i = 0; i < 5; ++i) {
    diff[i] = illConditionedReference.parameters[i] - illConditionedCandidate.parameters[i];
  }
  const float illDet = illD0 * illD1 - illC01 * illC01;
  const float illInverse00 = illD1 / illDet;
  const float illInverse01 = -illC01 / illDet;
  const float illInverse11 = illD0 / illDet;
  const float handIllConditionedChi2 = diff[0] * (illInverse00 * diff[0] + illInverse01 * diff[1]) +
                                       diff[1] * (illInverse01 * diff[0] + illInverse11 * diff[1]) +
                                       diff[2] * diff[2] / 0.02f + diff[3] * diff[3] / 0.03f + diff[4] * diff[4] / 0.04f;
  BOOST_TEST_MESSAGE("stateChi2 ill-conditioned characterization: production=" << illConditionedChi2
                                                                               << " hand=" << handIllConditionedChi2
                                                                               << " determinant=" << illDet);
  BOOST_CHECK_CLOSE_FRACTION(illConditionedChi2, handIllConditionedChi2, 5.e-3f); // Loose: characterizes conditioning, not a proof.
}

BOOST_AUTO_TEST_CASE(StateChi2RejectsInvalidInputsAndPreservesState)
{
  constexpr float InitialChi2 = 77.f;

  auto wrongFamilyCandidate = makeCandidateState();
  wrongFamilyCandidate.family = StateFamily::Forward;
  checkStateChi2FailurePreservesBytes(makeState(), wrongFamilyCandidate, InitialChi2, OperationFailureReason::SourceFamilyMismatch);

  auto wrongFamilyReference = makeState();
  wrongFamilyReference.family = StateFamily::Forward;
  checkStateChi2FailurePreservesBytes(wrongFamilyReference, makeCandidateState(), InitialChi2, OperationFailureReason::SourceFamilyMismatch);

  auto nonFiniteReference = makeState();
  nonFiniteReference.parameters[2] = std::numeric_limits<float>::quiet_NaN();
  checkStateChi2FailurePreservesBytes(nonFiniteReference, makeCandidateState(), InitialChi2, OperationFailureReason::NonFiniteInput);

  auto nonFiniteCandidate = makeCandidateState();
  nonFiniteCandidate.covariance[packedCovarianceIndex(3, 3)] = std::numeric_limits<float>::infinity();
  checkStateChi2FailurePreservesBytes(makeState(), nonFiniteCandidate, InitialChi2, OperationFailureReason::NonFiniteInput);

  auto alphaMismatch = makeCandidateState();
  alphaMismatch.alpha += 0.01f;
  checkStateChi2FailurePreservesBytes(makeState(), alphaMismatch, InitialChi2, OperationFailureReason::AlphaMismatch);

  auto referenceXMismatch = makeCandidateState();
  referenceXMismatch.referenceCoordinate += 0.5f;
  checkStateChi2FailurePreservesBytes(makeState(), referenceXMismatch, InitialChi2, OperationFailureReason::ReferenceCoordinateMismatch);

  auto singularReference = makeState();
  auto singularCandidate = makeCandidateState();
  for (float& value : singularReference.covariance) {
    value = 0.f;
  }
  for (float& value : singularCandidate.covariance) {
    value = 0.f;
  }
  checkStateChi2FailurePreservesBytes(singularReference, singularCandidate, InitialChi2, OperationFailureReason::InvalidCovariance);

  auto nonFiniteOutputReference = makeState();
  nonFiniteOutputReference.parameters[0] = std::numeric_limits<float>::max();
  auto nonFiniteOutputCandidate = makeCandidateState();
  nonFiniteOutputCandidate.parameters[0] = -std::numeric_limits<float>::max();
  checkStateChi2FailurePreservesBytes(nonFiniteOutputReference, nonFiniteOutputCandidate, InitialChi2, OperationFailureReason::NonFiniteOutput);
}

BOOST_AUTO_TEST_CASE(StateChi2PreservesInputsOnSuccessAndIsByteDeterministic)
{
  const auto reference = makeState();
  const auto candidate = makeCandidateState();
  const auto referenceBefore = reference;
  const auto candidateBefore = candidate;
  float firstChi2 = -1.f;
  float secondChi2 = -2.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(stateChi2(reference, candidate, firstChi2, reason));
  BOOST_REQUIRE(stateChi2(reference, candidate, secondChi2, reason));
  BOOST_CHECK(bitEqual(reference, referenceBefore));
  BOOST_CHECK(bitEqual(candidate, candidateBefore));
  BOOST_CHECK_EQUAL(firstChi2, secondChi2);
}

// --- barrel::buildSeed (Stage-B Slice A) ---

BOOST_AUTO_TEST_CASE(BuildSeedMatchesRetainedLegacyOracleNonzeroField)
{
  for (const bool mirror : {false, true}) {
    const auto clusterInner = makeInnerCluster(mirror);
    const auto clusterMiddle = makeMiddleCluster(mirror);
    const auto hitOuter = makeOuterHit(mirror);
    const float bz = 0.5f;
    const auto oracle = o2::its::track::buildTrackSeed(clusterInner, clusterMiddle, hitOuter, bz);

    SurfaceKinematicState outState{};
    OperationFailureReason reason{};
    BOOST_REQUIRE(barrel::buildSeed(globalPointFromCluster(clusterInner), globalPointFromCluster(clusterMiddle),
                                    measurementFromOuterHit(hitOuter), bz, 1, o2::track::PID::Pion, outState, reason));

    Drift drift{};
    checkClose(outState.referenceCoordinate, oracle.getX(), drift);
    checkClose(outState.alpha, oracle.getAlpha(), drift);
    for (uint8_t i = 0; i < 5; ++i) {
      checkClose(outState.parameters[i], oracle.getParam(i), drift);
    }
    for (uint8_t i = 0; i < 15; ++i) {
      checkClose(outState.covariance[i], oracle.getCov()[i], drift);
    }
    checkBuildSeedMetadata(outState, 1, o2::track::PID::Pion);
    BOOST_TEST_MESSAGE("buildSeed nonzero-field (mirror=" << mirror << ") retained-oracle max drift: abs=" << drift.absolute << " rel=" << drift.relative);
  }
}

BOOST_AUTO_TEST_CASE(BuildSeedMatchesRetainedLegacyOracleZeroField)
{
  const auto clusterInner = makeInnerCluster();
  const auto clusterMiddle = makeMiddleCluster();
  const auto hitOuter = makeOuterHit();
  const float bz = 0.005f; // |bz| < 0.01f -> zero-field branch
  const auto oracle = o2::its::track::buildTrackSeed(clusterInner, clusterMiddle, hitOuter, bz);

  SurfaceKinematicState outState{};
  OperationFailureReason reason{};
  BOOST_REQUIRE(barrel::buildSeed(globalPointFromCluster(clusterInner), globalPointFromCluster(clusterMiddle),
                                  measurementFromOuterHit(hitOuter), bz, 2, o2::track::PID::Kaon, outState, reason));

  Drift drift{};
  checkClose(outState.referenceCoordinate, oracle.getX(), drift);
  checkClose(outState.alpha, oracle.getAlpha(), drift);
  for (uint8_t i = 0; i < 5; ++i) {
    checkClose(outState.parameters[i], oracle.getParam(i), drift);
  }
  for (uint8_t i = 0; i < 15; ++i) {
    checkClose(outState.covariance[i], oracle.getCov()[i], drift);
  }
  checkBuildSeedMetadata(outState, 2, o2::track::PID::Kaon);
  BOOST_TEST_MESSAGE("buildSeed zero-field retained-oracle max drift: abs=" << drift.absolute << " rel=" << drift.relative);
}

BOOST_AUTO_TEST_CASE(BuildSeedCurvatureSignFlipsWithMirroredGeometry)
{
  // Dedicated alpha=0 fixture (not makeOuterHit's alpha=0.15f): with a zero
  // rotation, clusterInner.y/clusterMiddle.y/hitOuter's Y3 map identically
  // onto computeCurvature's y1/y2/y3, so negating all three is an exact
  // reflection about the local X axis. Under that reflection the retained
  // formula's signed area (hence crv and q2pt) flips sign while
  // computeCurvatureCentreX -- and therefore snp -- transforms consistently
  // with it; this is not true for a mirror applied before a nonzero-alpha
  // rotation (rotation mixes x and y), which is why the oracle-match
  // fixtures above use their own alpha and are not reused here.
  const float bz = 0.5f;
  const o2::its::Cluster clusterInner{1.8f, 0.10f, -0.60f, 0};
  const o2::its::Cluster clusterMiddle{3.0f, 0.35f, 0.10f, 1};
  const o2::its::TrackingFrameInfo hitOuter{0.f, 0.f, 0.f, 4.4f, 0.f, {0.5f, -0.8f}, {0.02f, 0.001f, 0.03f}};
  const o2::its::Cluster clusterInnerMirrored{1.8f, -0.10f, -0.60f, 0};
  const o2::its::Cluster clusterMiddleMirrored{3.0f, -0.35f, 0.10f, 1};
  const o2::its::TrackingFrameInfo hitOuterMirrored{0.f, 0.f, 0.f, 4.4f, 0.f, {-0.5f, -0.8f}, {0.02f, 0.001f, 0.03f}};

  SurfaceKinematicState plain{};
  SurfaceKinematicState mirrored{};
  OperationFailureReason reason{};
  BOOST_REQUIRE(barrel::buildSeed(globalPointFromCluster(clusterInner), globalPointFromCluster(clusterMiddle),
                                  measurementFromOuterHit(hitOuter), bz, 1, o2::track::PID::Pion, plain, reason));
  BOOST_REQUIRE(barrel::buildSeed(globalPointFromCluster(clusterInnerMirrored), globalPointFromCluster(clusterMiddleMirrored),
                                  measurementFromOuterHit(hitOuterMirrored), bz,
                                  1, o2::track::PID::Pion, mirrored, reason));
  BOOST_CHECK_LT(plain.parameters[4] * mirrored.parameters[4], 0.f);
  BOOST_CHECK_LT(plain.parameters[2] * mirrored.parameters[2], 0.f);
}

BOOST_AUTO_TEST_CASE(BuildSeedRejectsNonFiniteInput)
{
  const auto clusterInner = makeInnerCluster();
  const auto clusterMiddle = makeMiddleCluster();
  const auto hitOuter = makeOuterHit();
  const float bz = 0.5f;

  auto badInner = clusterInner;
  badInner.xCoordinate = std::numeric_limits<float>::quiet_NaN();
  checkBuildSeedFailurePreservesBytes(badInner, clusterMiddle, hitOuter, bz, OperationFailureReason::NonFiniteInput);

  auto badMiddle = clusterMiddle;
  badMiddle.zCoordinate = std::numeric_limits<float>::infinity();
  checkBuildSeedFailurePreservesBytes(clusterInner, badMiddle, hitOuter, bz, OperationFailureReason::NonFiniteInput);

  auto badHit = hitOuter;
  badHit.alphaTrackingFrame = std::numeric_limits<float>::quiet_NaN();
  checkBuildSeedFailurePreservesBytes(clusterInner, clusterMiddle, badHit, bz, OperationFailureReason::NonFiniteInput);

  checkBuildSeedFailurePreservesBytes(clusterInner, clusterMiddle, hitOuter, std::numeric_limits<float>::quiet_NaN(),
                                      OperationFailureReason::NonFiniteInput);
}

BOOST_AUTO_TEST_CASE(BuildSeedDegenerateZeroFieldGeometryRejectsViaNonFiniteOutput)
{
  // alpha=0 makes the rotation identity; clusterInner is chosen to coincide
  // exactly with (x3, y3) in that frame, forcing dx=dy=0 in the zero-field
  // branch (0.f/hypot(0,0) == NaN) -- an inherent numeric singularity of the
  // retained formula itself, not a new rejection added by this operation.
  const o2::its::TrackingFrameInfo hitOuter{0.f, 0.f, 0.f, 4.0f, 0.f, {0.5f, -0.8f}, {0.02f, 0.001f, 0.03f}};
  const o2::its::Cluster clusterInner{4.0f, 0.5f, -0.6f, 0};
  const auto clusterMiddle = makeMiddleCluster();
  checkBuildSeedFailurePreservesBytes(clusterInner, clusterMiddle, hitOuter, 0.f, OperationFailureReason::NonFiniteOutput);
}

BOOST_AUTO_TEST_CASE(BuildSeedIsByteDeterministic)
{
  const auto clusterInner = makeInnerCluster();
  const auto clusterMiddle = makeMiddleCluster();
  const auto hitOuter = makeOuterHit();
  const float bz = -0.35f;

  SurfaceKinematicState firstState{};
  SurfaceKinematicState secondState{};
  OperationFailureReason reason{};
  BOOST_REQUIRE(barrel::buildSeed(globalPointFromCluster(clusterInner), globalPointFromCluster(clusterMiddle),
                                  measurementFromOuterHit(hitOuter), bz, 1, o2::track::PID::Pion, firstState, reason));
  BOOST_REQUIRE(barrel::buildSeed(globalPointFromCluster(clusterInner), globalPointFromCluster(clusterMiddle),
                                  measurementFromOuterHit(hitOuter), bz, 1, o2::track::PID::Pion, secondState, reason));
  BOOST_CHECK(bitEqual(firstState, secondState));
}
