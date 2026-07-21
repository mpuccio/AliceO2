// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFTForwardSurfaceStateOperations
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"

namespace
{
using namespace o2::itsmft::tracking;

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
  state.absCharge = 3;
  state.pid = o2::track::PID::Kaon;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] =
        row == column ? 0.01f * (row + 1) : 0.0002f * (1 + packedCovarianceIndex(row, column));
    }
  }
  return state;
}

SurfaceMeasurement makeMeasurement()
{
  SurfaceMeasurement measurement{};
  measurement.global = {0.8f, -0.45f, -50.f};
  measurement.frame = {999.f, 777.f, 555.f, 1.2f}; // Deliberately not disk u/v.
  measurement.covariance = {0.04f, 0.012f, 0.09f};
  return measurement;
}

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

o2::track::TrackParCovFwd makeOracle(const SurfaceKinematicState& state)
{
  o2::track::SMatrix5 parameters{};
  o2::track::SMatrix55Sym covariance{};
  for (uint8_t i = 0; i < 5; ++i) {
    parameters(i) = state.parameters[i];
  }
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      covariance(row, column) = state.covariance[packedCovarianceIndex(row, column)];
    }
  }
  return {state.referenceCoordinate, parameters, covariance, 0.};
}

struct OracleDrift {
  double maximumAbsolute{0.};
  double maximumRelative{0.};
};

OracleDrift parameterDrift;
OracleDrift covarianceDrift;
OracleDrift referenceDrift;
OracleDrift chi2Drift;

void checkClose(float value, double oracle, double absoluteTolerance, double relativeTolerance, OracleDrift& drift)
{
  const double difference = std::abs(static_cast<double>(value) - oracle);
  drift.maximumAbsolute = std::max(drift.maximumAbsolute, difference);
  if (oracle != 0.) {
    drift.maximumRelative = std::max(drift.maximumRelative, difference / std::abs(oracle));
  }
  BOOST_CHECK_LE(difference, absoluteTolerance + relativeTolerance * std::abs(oracle));
}

void compareWithRetainedLegacyOracle(const SurfaceKinematicState& state, const o2::track::TrackParCovFwd& oracle)
{
  // This tolerance comparison characterizes float drift against the retained
  // legacy double oracle. It does not establish physical correctness.
  const double oracleParameters[5] = {oracle.getX(), oracle.getY(), oracle.getPhi(), oracle.getTanl(), oracle.getInvQPt()};
  for (uint8_t i = 0; i < 5; ++i) {
    checkClose(state.parameters[i], oracleParameters[i], 3.e-5, 3.e-5, parameterDrift);
  }
  checkClose(state.referenceCoordinate, oracle.getZ(), 1.e-6, 1.e-6, referenceDrift);
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      checkClose(state.covariance[packedCovarianceIndex(row, column)], oracle.getCovariances()(row, column), 3.e-6, 2.e-4, covarianceDrift);
    }
  }
}

template <ForwardPropagationModel Model>
void comparePropagationWithOracle(float targetZ, float bz)
{
  auto state = makeState();
  auto oracle = makeOracle(state);
  if constexpr (Model == ForwardPropagationModel::Linear) {
    oracle.propagateToZlinear(targetZ);
  } else if constexpr (Model == ForwardPropagationModel::Quadratic) {
    oracle.propagateToZquadratic(targetZ, bz);
  } else if constexpr (Model == ForwardPropagationModel::Helix) {
    oracle.propagateToZhelix(targetZ, bz);
  } else {
    oracle.propagateToZ(targetZ, bz);
  }
  OperationFailureReason reason{};
  BOOST_REQUIRE(propagateToDisk<Model>(state, targetZ, bz, reason));
  compareWithRetainedLegacyOracle(state, oracle);
}

float directPredictedChi2(const SurfaceKinematicState& state, const SurfaceMeasurement& measurement)
{
  const float s00 = state.covariance[packedCovarianceIndex(0, 0)] + measurement.covariance.uu;
  const float s01 = state.covariance[packedCovarianceIndex(1, 0)] + measurement.covariance.uv;
  const float s11 = state.covariance[packedCovarianceIndex(1, 1)] + measurement.covariance.vv;
  const float determinant = s00 * s11 - s01 * s01;
  const float dx = measurement.global.x - state.parameters[0];
  const float dy = measurement.global.y - state.parameters[1];
  return (s11 * dx * dx - 2.f * s01 * dx * dy + s00 * dy * dy) / determinant;
}

template <typename Callable>
void checkStateFailurePreservesBytes(SurfaceKinematicState state, OperationFailureReason expected, Callable&& operation)
{
  const auto before = state;
  OperationFailureReason reason{};
  BOOST_CHECK(!operation(state, reason));
  BOOST_CHECK(reason == expected);
  BOOST_CHECK(bitEqual(state, before));
}

} // namespace

BOOST_AUTO_TEST_CASE(LinearPropagationMatchesHandCalculationAndPackedCovariance)
{
  auto state = makeState();
  for (float& value : state.covariance) {
    value = 0.f;
  }
  state.covariance[packedCovarianceIndex(2, 2)] = 0.25f;
  const auto before = state;
  constexpr float targetZ = -50.f;
  const float n = (targetZ - before.referenceCoordinate) / before.parameters[3];
  OperationFailureReason reason{};
  BOOST_REQUIRE(propagateToDisk<ForwardPropagationModel::Linear>(state, targetZ, 5.f, reason));
  BOOST_CHECK_CLOSE_FRACTION(state.parameters[0], before.parameters[0] + n * std::cos(before.parameters[2]), 2.e-6f);
  BOOST_CHECK_CLOSE_FRACTION(state.parameters[1], before.parameters[1] + n * std::sin(before.parameters[2]), 2.e-6f);
  BOOST_CHECK_EQUAL(state.referenceCoordinate, targetZ);
  BOOST_CHECK_CLOSE_FRACTION(state.covariance[packedCovarianceIndex(0, 0)],
                             n * n * std::sin(before.parameters[2]) * std::sin(before.parameters[2]) * 0.25f, 3.e-6f);
  BOOST_CHECK_CLOSE_FRACTION(state.covariance[packedCovarianceIndex(1, 0)],
                             -n * n * std::sin(before.parameters[2]) * std::cos(before.parameters[2]) * 0.25f, 3.e-6f);
  BOOST_CHECK_CLOSE_FRACTION(state.covariance[packedCovarianceIndex(1, 1)],
                             n * n * std::cos(before.parameters[2]) * std::cos(before.parameters[2]) * 0.25f, 3.e-6f);
}

BOOST_AUTO_TEST_CASE(ZeroFieldOptimizedIsExactlyLinearInBothDirections)
{
  for (const float targetZ : {-60.f, -35.f}) {
    auto linear = makeState();
    auto optimized = linear;
    OperationFailureReason reason{};
    BOOST_REQUIRE(propagateToDisk<ForwardPropagationModel::Linear>(linear, targetZ, 0.f, reason));
    BOOST_REQUIRE(propagateToDisk<ForwardPropagationModel::Optimized>(optimized, targetZ, 0.f, reason));
    BOOST_CHECK(bitEqual(linear, optimized));
  }
}

BOOST_AUTO_TEST_CASE(NonzeroFieldModelsMatchRetainedLegacyOracle)
{
  comparePropagationWithOracle<ForwardPropagationModel::Linear>(-52.f, 5.f);
  comparePropagationWithOracle<ForwardPropagationModel::Quadratic>(-52.f, 5.f);
  comparePropagationWithOracle<ForwardPropagationModel::Helix>(-52.f, 5.f);
  comparePropagationWithOracle<ForwardPropagationModel::Optimized>(-52.f, 5.f);
  comparePropagationWithOracle<ForwardPropagationModel::Helix>(-38.f, -5.f);
}

BOOST_AUTO_TEST_CASE(NearZeroNonzeroFieldRetainsOptimizedHelixSelection)
{
  comparePropagationWithOracle<ForwardPropagationModel::Optimized>(-45.01f, 1.e-4f);
}

BOOST_AUTO_TEST_CASE(PredictedChi2UsesDiskGlobalXYAndNonzeroUV)
{
  const auto state = makeState();
  const auto measurement = makeMeasurement();
  float chi2 = -1.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(predictedChi2(state, measurement, chi2, reason));
  BOOST_CHECK_CLOSE_FRACTION(chi2, directPredictedChi2(state, measurement), 2.e-6f);
  BOOST_CHECK_LT(chi2, 100.f); // Would be enormous if frame.u/frame.v were used.
}

BOOST_AUTO_TEST_CASE(UpdateMatchesDirectTwoDimensionalKalmanReference)
{
  auto state = makeState();
  const auto before = state;
  const auto measurement = makeMeasurement();
  const float s00 = before.covariance[packedCovarianceIndex(0, 0)] + measurement.covariance.uu;
  const float s01 = before.covariance[packedCovarianceIndex(1, 0)] + measurement.covariance.uv;
  const float s11 = before.covariance[packedCovarianceIndex(1, 1)] + measurement.covariance.vv;
  const float determinant = s00 * s11 - s01 * s01;
  const float inv00 = s11 / determinant;
  const float inv01 = -s01 / determinant;
  const float inv11 = s00 / determinant;
  const float residual[2] = {measurement.global.x - before.parameters[0], measurement.global.y - before.parameters[1]};
  float expectedParameters[5]{};
  float expectedCovariance[5][5]{};
  for (uint8_t row = 0; row < 5; ++row) {
    const float c0 = before.covariance[packedCovarianceIndex(row, 0)];
    const float c1 = before.covariance[packedCovarianceIndex(row, 1)];
    const float gain0 = c0 * inv00 + c1 * inv01;
    const float gain1 = c0 * inv01 + c1 * inv11;
    expectedParameters[row] = before.parameters[row] + gain0 * residual[0] + gain1 * residual[1];
    for (uint8_t column = 0; column < 5; ++column) {
      expectedCovariance[row][column] = before.covariance[packedCovarianceIndex(row, column)] -
                                        gain0 * before.covariance[packedCovarianceIndex(0, column)] -
                                        gain1 * before.covariance[packedCovarianceIndex(1, column)];
    }
  }
  float chi2 = -1.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(update(state, measurement, chi2, reason));
  BOOST_CHECK_CLOSE_FRACTION(chi2, directPredictedChi2(before, measurement), 2.e-6f);
  for (uint8_t row = 0; row < 5; ++row) {
    BOOST_CHECK_CLOSE_FRACTION(state.parameters[row], expectedParameters[row], 2.e-6f);
    for (uint8_t column = 0; column <= row; ++column) {
      BOOST_CHECK_CLOSE_FRACTION(state.covariance[packedCovarianceIndex(row, column)], expectedCovariance[row][column], 3.e-6f);
    }
  }
}

BOOST_AUTO_TEST_CASE(DiagonalMeasurementUpdateMatchesRetainedLegacyOracle)
{
  auto state = makeState();
  auto oracle = makeOracle(state);
  auto measurement = makeMeasurement();
  measurement.covariance.uv = 0.f;
  const std::array<float, 2> position = {measurement.global.x, measurement.global.y};
  const std::array<float, 2> covariance = {measurement.covariance.uu, measurement.covariance.vv};
  BOOST_REQUIRE(oracle.update(position, covariance));
  float chi2 = -1.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(update(state, measurement, chi2, reason));
  // Tolerance comparison against a retained legacy oracle; this does not
  // claim that the oracle proves physical correctness of the float port.
  compareWithRetainedLegacyOracle(state, oracle);
  checkClose(chi2, oracle.getTrackChi2(), 3.e-6, 2.e-4, chi2Drift);
}

BOOST_AUTO_TEST_CASE(MaterialMatchesHandCalculationAndRetainedLegacyOracle)
{
  auto state = makeState();
  const auto before = state;
  auto oracle = makeOracle(state);
  constexpr float xOverX0 = 0.018f;
  oracle.addMCSEffect(xOverX0);
  OperationFailureReason reason{};
  BOOST_REQUIRE(correctForMaterial(state, xOverX0, reason));
  const float a = 1.f + before.parameters[3] * before.parameters[3];
  const float inverseMomentum = std::abs(before.parameters[4]) / std::sqrt(a);
  const float path = xOverX0 * std::abs(std::sqrt(a) / before.parameters[3]);
  const float theta2 = highlandTheta2(inverseMomentum, path);
  BOOST_CHECK_CLOSE_FRACTION(state.covariance[packedCovarianceIndex(2, 2)], before.covariance[packedCovarianceIndex(2, 2)] + theta2 * a, 2.e-6f);
  BOOST_CHECK_CLOSE_FRACTION(state.covariance[packedCovarianceIndex(3, 3)], before.covariance[packedCovarianceIndex(3, 3)] + theta2 * a * a, 2.e-6f);
  BOOST_CHECK_CLOSE_FRACTION(state.covariance[packedCovarianceIndex(4, 4)], before.covariance[packedCovarianceIndex(4, 4)] + theta2 * before.parameters[3] * before.parameters[3] * before.parameters[4] * before.parameters[4], 2.e-6f);
  BOOST_CHECK_EQUAL(state.absCharge, before.absCharge);
  BOOST_CHECK_EQUAL(state.pid.getID(), before.pid.getID());
  compareWithRetainedLegacyOracle(state, oracle);
}

BOOST_AUTO_TEST_CASE(FailuresPreserveEveryDestinationByte)
{
  auto wrongFamily = makeState();
  wrongFamily.family = StateFamily::Barrel;
  checkStateFailurePreservesBytes(wrongFamily, OperationFailureReason::SourceFamilyMismatch,
                                  [](auto& state, auto& reason) { return propagateToDisk<ForwardPropagationModel::Linear>(state, -50.f, 0.f, reason); });

  auto nonFinite = makeState();
  nonFinite.parameters[0] = std::numeric_limits<float>::quiet_NaN();
  checkStateFailurePreservesBytes(nonFinite, OperationFailureReason::NonFiniteInput,
                                  [](auto& state, auto& reason) { return correctForMaterial(state, 0.01f, reason); });

  auto unreachable = makeState();
  unreachable.parameters[3] = 0.f;
  checkStateFailurePreservesBytes(unreachable, OperationFailureReason::UnreachableTarget,
                                  [](auto& state, auto& reason) { return propagateToDisk<ForwardPropagationModel::Quadratic>(state, -50.f, 5.f, reason); });

  checkStateFailurePreservesBytes(makeState(), OperationFailureReason::PropagationFailure,
                                  [](auto& state, auto& reason) { return propagateToDisk<ForwardPropagationModel::Helix>(state, -50.f, 0.f, reason); });

  auto material = makeState();
  material.parameters[3] = 0.f;
  checkStateFailurePreservesBytes(material, OperationFailureReason::MaterialFailure,
                                  [](auto& state, auto& reason) { return correctForMaterial(state, 0.01f, reason); });

  SurfaceMeasurement singular{};
  auto singularState = makeState();
  for (float& value : singularState.covariance) {
    value = 0.f;
  }
  float chi2 = 123.f;
  const float chi2Before = chi2;
  OperationFailureReason reason{};
  BOOST_CHECK(!predictedChi2(singularState, singular, chi2, reason));
  BOOST_CHECK(reason == OperationFailureReason::InvalidCovariance);
  BOOST_CHECK_EQUAL(chi2, chi2Before);
  const auto updateBefore = singularState;
  BOOST_CHECK(!update(singularState, singular, chi2, reason));
  BOOST_CHECK(reason == OperationFailureReason::InvalidCovariance);
  BOOST_CHECK(bitEqual(singularState, updateBefore));
  BOOST_CHECK_EQUAL(chi2, chi2Before);

  auto invalidMeasurement = makeMeasurement();
  invalidMeasurement.covariance.uv = std::numeric_limits<float>::infinity();
  const auto stateBefore = singularState;
  BOOST_CHECK(!update(singularState, invalidMeasurement, chi2, reason));
  BOOST_CHECK(reason == OperationFailureReason::NonFiniteInput);
  BOOST_CHECK(bitEqual(singularState, stateBefore));
  BOOST_CHECK_EQUAL(chi2, chi2Before);

  auto nonFiniteOutput = makeState();
  nonFiniteOutput.parameters[4] = std::numeric_limits<float>::max();
  checkStateFailurePreservesBytes(nonFiniteOutput, OperationFailureReason::NonFiniteOutput,
                                  [](auto& state, auto& reason) { return correctForMaterial(state, std::numeric_limits<float>::max(), reason); });

  auto failedChi2 = makeState();
  failedChi2.covariance[packedCovarianceIndex(0, 0)] = std::numeric_limits<float>::max();
  failedChi2.covariance[packedCovarianceIndex(1, 1)] = std::numeric_limits<float>::max();
  const auto failedChi2Before = failedChi2;
  auto finiteMeasurement = makeMeasurement();
  BOOST_CHECK(!predictedChi2(failedChi2, finiteMeasurement, chi2, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(failedChi2, failedChi2Before));
  BOOST_CHECK_EQUAL(chi2, chi2Before);
  BOOST_CHECK(!update(failedChi2, finiteMeasurement, chi2, reason));
  BOOST_CHECK(reason == OperationFailureReason::UpdateFailure);
  BOOST_CHECK(bitEqual(failedChi2, failedChi2Before));
  BOOST_CHECK_EQUAL(chi2, chi2Before);

  checkStateFailurePreservesBytes(makeState(), OperationFailureReason::NonFiniteInput,
                                  [](auto& state, auto& reason) { return propagateToDisk<ForwardPropagationModel::Linear>(state, std::numeric_limits<float>::infinity(), 0.f, reason); });
}

BOOST_AUTO_TEST_CASE(RepeatedMultiStepChainsAreByteDeterministicAndCharacterizeOracleDrift)
{
  auto runChain = [] {
    auto state = makeState();
    OperationFailureReason reason{};
    float chi2 = 0.f;
    for (int step = 0; step < 5; ++step) {
      const float z = -48.f - 2.f * step;
      BOOST_REQUIRE(propagateToDisk<ForwardPropagationModel::Optimized>(state, z, 5.f, reason));
      BOOST_REQUIRE(correctForMaterial(state, 0.004f + 0.001f * step, reason));
      auto measurement = makeMeasurement();
      measurement.global.x = state.parameters[0] + 0.01f * (step + 1);
      measurement.global.y = state.parameters[1] - 0.015f * (step + 1);
      float increment = 0.f;
      BOOST_REQUIRE(update(state, measurement, increment, reason));
      chi2 += increment;
    }
    return std::pair{state, chi2};
  };
  const auto first = runChain();
  const auto second = runChain();
  BOOST_CHECK(bitEqual(first.first, second.first));
  BOOST_CHECK_EQUAL(first.second, second.second);

  // Retained-oracle drift characterization for the propagation/material part
  // of the same repeated chain; uv-aware updates have no legacy counterpart.
  auto state = makeState();
  auto oracle = makeOracle(state);
  OperationFailureReason reason{};
  for (int step = 0; step < 5; ++step) {
    const float z = -48.f - 2.f * step;
    BOOST_REQUIRE(propagateToDisk<ForwardPropagationModel::Optimized>(state, z, 5.f, reason));
    BOOST_REQUIRE(correctForMaterial(state, 0.004f + 0.001f * step, reason));
    oracle.propagateToZ(z, 5.f);
    oracle.addMCSEffect(0.004f + 0.001f * step);
  }
  compareWithRetainedLegacyOracle(state, oracle);
}

BOOST_AUTO_TEST_CASE(NewProductionFilesHaveNoLegacyForwardDependency)
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const std::array<std::string, 2> productionFiles = {
    testDirectory + "/../include/ITSMFTTracking/ForwardSurfaceStateOperations.h",
    testDirectory + "/../src/ForwardSurfaceStateOperations.cxx"};
  const std::array<std::string, 4> forbidden = {
    "SurfaceKinematicStateLegacyAdapters.h", "TrackFwd.h", "TrackParCovFwd", "TrackParFwd"};
  for (const auto& file : productionFiles) {
    std::ifstream input{file};
    BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << file);
    std::ostringstream buffer;
    buffer << input.rdbuf();
    for (const auto& token : forbidden) {
      BOOST_CHECK_MESSAGE(buffer.str().find(token) == std::string::npos, file << " contains forbidden token " << token);
    }
  }
  BOOST_TEST_MESSAGE("retained-oracle maximum absolute/relative drift: parameters="
                     << parameterDrift.maximumAbsolute << "/" << parameterDrift.maximumRelative
                     << ", covariance=" << covarianceDrift.maximumAbsolute << "/" << covarianceDrift.maximumRelative
                     << ", reference=" << referenceDrift.maximumAbsolute << "/" << referenceDrift.maximumRelative
                     << ", update chi2=" << chi2Drift.maximumAbsolute << "/" << chi2Drift.maximumRelative);
}
