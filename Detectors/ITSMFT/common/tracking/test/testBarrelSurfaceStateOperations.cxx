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

#include <cmath>
#include <cstring>
#include <limits>

#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"

namespace
{
using namespace o2::itsmft::tracking;
using namespace o2::itsmft::tracking::barrel;

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
  state.absCharge = 1;
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
  measurement.global = {99.f, -77.f, 0.f}; // Deliberately not barrel u/v.
  measurement.frame = {123.f, 0.8f, -0.45f, 0.f};
  measurement.covariance = {0.04f, 0.012f, 0.09f};
  return measurement;
}

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

float directPredictedChi2(const SurfaceKinematicState& state, const SurfaceMeasurement& measurement)
{
  const float s00 = state.covariance[packedCovarianceIndex(0, 0)] + measurement.covariance.uu;
  const float s01 = state.covariance[packedCovarianceIndex(1, 0)] + measurement.covariance.uv;
  const float s11 = state.covariance[packedCovarianceIndex(1, 1)] + measurement.covariance.vv;
  const float determinant = s00 * s11 - s01 * s01;
  const float dy = measurement.frame.u - state.parameters[0];
  const float dz = measurement.frame.v - state.parameters[1];
  return (s11 * dy * dy - 2.f * s01 * dy * dz + s00 * dz * dz) / determinant;
}

template <typename Callable>
void checkStateFailurePreservesBytes(SurfaceKinematicState state, OperationFailureReason expected, Callable&& operation)
{
  const auto before = state;
  OperationFailureReason reason{};
  BOOST_CHECK(!operation(state, reason));
  BOOST_CHECK_EQUAL(static_cast<uint8_t>(reason), static_cast<uint8_t>(expected));
  BOOST_CHECK(bitEqual(state, before));
}

} // namespace

BOOST_AUTO_TEST_CASE(RotationTransformsFrameAndMatchesRetainedLegacyOracle)
{
  auto state = makeState();
  o2::track::TrackParCovF oracle{};
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(state, oracle));
  constexpr float targetAlpha = 0.55f;
  BOOST_REQUIRE(oracle.rotate(targetAlpha));
  OperationFailureReason reason{};
  BOOST_REQUIRE(rotate(state, targetAlpha, reason));
  BOOST_CHECK_EQUAL(state.alpha, targetAlpha);
  BOOST_CHECK_CLOSE_FRACTION(state.referenceCoordinate, oracle.getX(), 2.e-5f);
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_CLOSE_FRACTION(state.parameters[i], oracle.getParam(i), 2.e-5f);
  }
}

BOOST_AUTO_TEST_CASE(PropagationUsesTargetXAndMatchesRetainedLegacyOracle)
{
  auto state = makeState();
  o2::track::TrackParCovF oracle{};
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(state, oracle));
  constexpr float targetX = 6.5f;
  constexpr float bz = 5.f;
  BOOST_REQUIRE(oracle.propagateTo(targetX, bz));
  OperationFailureReason reason{};
  BOOST_REQUIRE(propagate(state, targetX, bz, reason));
  BOOST_CHECK_EQUAL(state.referenceCoordinate, targetX);
  BOOST_CHECK_EQUAL(state.alpha, 0.3f);
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_CLOSE_FRACTION(state.parameters[i], oracle.getParam(i), 2.e-4f);
  }
}

BOOST_AUTO_TEST_CASE(PredictedChi2AndUpdateUseFrameUV)
{
  const auto before = makeState();
  const auto measurement = makeMeasurement();
  float chi2 = -1.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(predictedChi2(before, measurement, chi2, reason));
  BOOST_CHECK_CLOSE_FRACTION(chi2, directPredictedChi2(before, measurement), 2.e-6f);

  auto state = before;
  BOOST_REQUIRE(update(state, measurement, chi2, reason));
  BOOST_CHECK_CLOSE_FRACTION(chi2, directPredictedChi2(before, measurement), 2.e-6f);
  BOOST_CHECK_NE(state.parameters[0], before.parameters[0]);
  BOOST_CHECK_NE(state.parameters[1], before.parameters[1]);
}

BOOST_AUTO_TEST_CASE(FailuresPreserveStateAndChi2)
{
  auto wrongFamily = makeState();
  wrongFamily.family = StateFamily::Forward;
  checkStateFailurePreservesBytes(wrongFamily, OperationFailureReason::SourceFamilyMismatch,
                                  [](auto& state, auto& reason) { return rotate(state, 0.5f, reason); });

  auto rotateFailure = makeState();
  rotateFailure.parameters[2] = 0.99f;
  checkStateFailurePreservesBytes(rotateFailure, OperationFailureReason::RotationFailure,
                                  [](auto& state, auto& reason) { return rotate(state, 2.f, reason); });

  auto propagationFailure = makeState();
  checkStateFailurePreservesBytes(propagationFailure, OperationFailureReason::UnreachableTarget,
                                  [](auto& state, auto& reason) { return propagate(state, 1000.f, 5.f, reason); });

  auto state = makeState();
  const auto before = state;
  auto measurement = makeMeasurement();
  measurement.covariance = {-before.covariance[packedCovarianceIndex(0, 0)], 0.f,
                            -before.covariance[packedCovarianceIndex(1, 1)]};
  float chi2 = 123.f;
  OperationFailureReason reason{};
  BOOST_CHECK(!update(state, measurement, chi2, reason));
  BOOST_CHECK_EQUAL(static_cast<uint8_t>(reason), static_cast<uint8_t>(OperationFailureReason::InvalidCovariance));
  BOOST_CHECK_EQUAL(chi2, 123.f);
  BOOST_CHECK(bitEqual(state, before));

  measurement = makeMeasurement();
  measurement.frame.u = std::numeric_limits<float>::infinity();
  BOOST_CHECK(!predictedChi2(state, measurement, chi2, reason));
  BOOST_CHECK_EQUAL(chi2, 123.f);
}
