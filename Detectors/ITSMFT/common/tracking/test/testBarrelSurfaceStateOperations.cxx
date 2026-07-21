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
#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"

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

SurfaceMeasurement makeMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.global = {99.f, -77.f, 123.f}; // Deliberately inconsistent with barrel frame u/v.
  measurement.frame = {456.f, 0.8f, -0.45f, 0.f};
  measurement.covariance = {0.04f, 0.012f, 0.09f};
  return measurement;
}

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
