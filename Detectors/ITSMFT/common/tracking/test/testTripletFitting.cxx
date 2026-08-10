// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFTTrackingTripletFitting
#include <boost/test/unit_test.hpp>

#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>

#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/TripletFitting.h"

using namespace o2::itsmft::tracking;

namespace
{

constexpr double Radius = 50.;
constexpr double TanLambda = 0.4;

SymmetricCovariance3D makeCovariance()
{
  // Positive definite, non-axis-aligned covariance in cm^2.
  return {4.e-6, 0.8e-6, -0.4e-6, 3.e-6, 0.3e-6, 5.e-6};
}

std::array<TripletFitObservation, 3> makeHelixObservations()
{
  const std::array<double, 3> angles{0.1, 0.16, 0.25};
  std::array<TripletFitObservation, 3> observations{};
  for (std::size_t i = 0; i < observations.size(); ++i) {
    observations[i].position = {3. + Radius * std::cos(angles[i]),
                                -2. + Radius * std::sin(angles[i]),
                                1.5 + Radius * angles[i] * TanLambda};
    observations[i].covariance = makeCovariance();
  }
  return observations;
}

SymmetricCovariance3D rotateCovarianceAroundZ(const SymmetricCovariance3D& covariance,
                                              double angle)
{
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  return {
    cosine * cosine * covariance.xx - 2. * sine * cosine * covariance.xy + sine * sine * covariance.yy,
    sine * cosine * covariance.xx + (cosine * cosine - sine * sine) * covariance.xy -
      sine * cosine * covariance.yy,
    cosine * covariance.xz - sine * covariance.yz,
    sine * sine * covariance.xx + 2. * sine * cosine * covariance.xy + cosine * cosine * covariance.yy,
    sine * covariance.xz + cosine * covariance.yz,
    covariance.zz};
}

void checkClose(double actual, double expected, double relativeTolerance, double absoluteTolerance = 0.)
{
  BOOST_CHECK_SMALL(actual - expected,
                    std::max(absoluteTolerance, relativeTolerance * std::max(std::abs(actual), std::abs(expected))));
}

std::string readFile(const std::filesystem::path& path)
{
  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path.string());
  return {std::istreambuf_iterator<char>{input}, {}};
}

} // namespace

BOOST_AUTO_TEST_CASE(CylinderObservationLiftsMeasuredUVIntoGlobalCovariance)
{
  SurfaceDescriptor surface{};
  surface.kind = SurfaceKind::Cylinder;
  surface.referenceCoordinate = 5.f;
  SurfaceMeasurement measurement{};
  measurement.frame = {4.f, 3.f, -2.f, static_cast<float>(std::atan2(4., 3.))};
  measurement.covariance = {4.f, 1.5f, 9.f};

  TripletFitObservation observation{};
  BOOST_REQUIRE(makeCylinderTripletFitObservation(surface, measurement, observation));
  const double sine = std::sin(measurement.frame.frameAngle);
  const double cosine = std::cos(measurement.frame.frameAngle);
  checkClose(observation.position[0], measurement.frame.q * cosine - measurement.frame.u * sine, 1.e-14);
  checkClose(observation.position[1], measurement.frame.q * sine + measurement.frame.u * cosine, 1.e-14);
  BOOST_CHECK_EQUAL(observation.position[2], -2.);
  checkClose(observation.covariance.xx, sine * sine * 4., 1.e-14);
  checkClose(observation.covariance.xy, -sine * cosine * 4., 1.e-14);
  checkClose(observation.covariance.xz, -sine * 1.5, 1.e-14);
  checkClose(observation.covariance.yy, cosine * cosine * 4., 1.e-14);
  checkClose(observation.covariance.yz, cosine * 1.5, 1.e-14);
  BOOST_CHECK_EQUAL(observation.covariance.zz, 9.);
}

BOOST_AUTO_TEST_CASE(DiskObservationUsesGlobalXYAndExactSurfaceZ)
{
  SurfaceDescriptor surface{};
  surface.kind = SurfaceKind::Disk;
  surface.referenceCoordinate = -46.f;
  SurfaceMeasurement measurement{};
  measurement.global = {2.f, -3.f, 999.f};
  measurement.covariance = {4.f, -1.f, 9.f};

  TripletFitObservation observation{};
  BOOST_REQUIRE(makeDiskTripletFitObservation(surface, measurement, observation));
  BOOST_CHECK_EQUAL(observation.position[0], 2.);
  BOOST_CHECK_EQUAL(observation.position[1], -3.);
  BOOST_CHECK_EQUAL(observation.position[2], -46.);
  BOOST_CHECK_EQUAL(observation.covariance.xx, 4.);
  BOOST_CHECK_EQUAL(observation.covariance.xy, -1.);
  BOOST_CHECK_EQUAL(observation.covariance.yy, 9.);
  BOOST_CHECK_EQUAL(observation.covariance.xz, 0.);
  BOOST_CHECK_EQUAL(observation.covariance.yz, 0.);
  BOOST_CHECK_EQUAL(observation.covariance.zz, 0.);
}

BOOST_AUTO_TEST_CASE(ExactHelixHasExpectedCurvatureAndZeroLocalQuality)
{
  const auto observations = makeHelixObservations();
  LocalTripletFitResult result{};
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(observations, {1.e-8}, result));

  const double expectedTransverseCurvature = 1. / Radius;
  const double expectedSinTheta = 1. / std::sqrt(1. + TanLambda * TanLambda);
  const double expectedCurvature = expectedTransverseCurvature * expectedSinTheta;
  checkClose(result.referenceTransverseCurvature, expectedTransverseCurvature, 2.e-13);
  checkClose(result.referenceSinTheta, expectedSinTheta, 3.e-4);
  checkClose(result.curvature, expectedCurvature, 3.e-4);
  BOOST_CHECK_SMALL(result.chi2, 1.e-18);
  BOOST_CHECK_GT(result.curvatureVariance, 0.);

  const double expectedPt = std::abs(static_cast<double>(o2::constants::math::B2C) * 5. * Radius);
  checkClose(fittedTripletTransverseMomentum(result, 5.), expectedPt, 4.e-4);
}

BOOST_AUTO_TEST_CASE(GlobalRotationAndTranslationLeaveFitInvariant)
{
  const auto original = makeHelixObservations();
  auto transformed = original;
  const double angle = 1.234;
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  for (auto& observation : transformed) {
    const double x = observation.position[0];
    const double y = observation.position[1];
    observation.position = {cosine * x - sine * y + 17.,
                            sine * x + cosine * y - 8.,
                            observation.position[2] + 23.};
    observation.covariance = rotateCovarianceAroundZ(observation.covariance, angle);
  }

  LocalTripletFitResult first{};
  LocalTripletFitResult second{};
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(original, {2.e-7}, first));
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(transformed, {2.e-7}, second));
  checkClose(second.curvature, first.curvature, 2.e-13);
  checkClose(second.curvatureVariance, first.curvatureVariance, 2.e-12);
  checkClose(second.chi2, first.chi2, 1.e-9, 1.e-20);
  checkClose(second.gammaThetaPhi, first.gammaThetaPhi, 3.e-12);
}

BOOST_AUTO_TEST_CASE(FullCrossCovarianceContributesToTheFit)
{
  const auto withCrossTerms = makeHelixObservations();
  auto diagonal = withCrossTerms;
  for (auto& observation : diagonal) {
    observation.covariance.xy = 0.;
    observation.covariance.xz = 0.;
    observation.covariance.yz = 0.;
  }

  LocalTripletFitResult full{};
  LocalTripletFitResult reduced{};
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(withCrossTerms, {3.e-8}, full));
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(diagonal, {3.e-8}, reduced));
  BOOST_CHECK_NE(full.gammaThetaPhi, reduced.gammaThetaPhi);
  BOOST_CHECK_NE(full.curvatureVariance, reduced.curvatureVariance);
}

BOOST_AUTO_TEST_CASE(StraightTripletUsesTheRemovableZeroBendingLimit)
{
  const SymmetricCovariance3D covariance{1.e-6, 0., 0., 1.e-6, 0., 1.e-6};
  const std::array<TripletFitObservation, 3> observations{{
    {{1., 2., 3.}, covariance},
    {{2., 2., 3.5}, covariance},
    {{4., 2., 4.5}, covariance},
  }};
  LocalTripletFitResult result{};
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(observations, {1.e-8}, result));
  BOOST_CHECK_SMALL(result.referenceTransverseCurvature, 1.e-15);
  BOOST_CHECK_SMALL(result.curvature, 1.e-15);
  BOOST_CHECK_SMALL(result.chi2, 1.e-18);
  BOOST_CHECK(std::isinf(fittedTripletTransverseMomentum(result, 5.)));
}

BOOST_AUTO_TEST_CASE(NonHelicalKinkHasPositiveQualityAndMSReducesIt)
{
  auto observations = makeHelixObservations();
  observations[2].position[2] += 0.04;
  LocalTripletFitResult measurementOnly{};
  LocalTripletFitResult withScattering{};
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(observations, {0.}, measurementOnly));
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(observations, {2.e-4}, withScattering));
  BOOST_CHECK_GT(measurementOnly.chi2, 0.);
  BOOST_CHECK_LT(withScattering.chi2, measurementOnly.chi2);
  BOOST_CHECK_GT(withScattering.curvatureVariance, measurementOnly.curvatureVariance);
}

BOOST_AUTO_TEST_CASE(InvalidInputsFailTransactionally)
{
  const LocalTripletFitResult sentinel{1., 2., 3., 4., 0.5, 6., 7., 8., 9., 10., 11., 12.};
  auto observations = makeHelixObservations();
  LocalTripletFitResult result = sentinel;

  observations[1].covariance.xy = 1.;
  BOOST_CHECK(!fitLocalTripletUniformSolenoid(observations, {0.}, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations = makeHelixObservations();
  observations[2].position[0] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK(!fitLocalTripletUniformSolenoid(observations, {0.}, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations = makeHelixObservations();
  observations[1].position = observations[0].position;
  BOOST_CHECK(!fitLocalTripletUniformSolenoid(observations, {0.}, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations = makeHelixObservations();
  BOOST_CHECK(!fitLocalTripletUniformSolenoid(observations, {-1.}, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);
}

BOOST_AUTO_TEST_CASE(StandaloneFitHasNoFamilyDispatchOrTrackerTraitsCallSite)
{
  const auto trackingRoot = std::filesystem::path{__FILE__}.parent_path().parent_path();
  const auto source = readFile(trackingRoot / "src/TripletFitting.cxx");
  const auto begin = source.find("bool fitLocalTripletUniformSolenoid(");
  const auto end = source.find("double fittedTripletTransverseMomentum(", begin);
  BOOST_REQUIRE_NE(begin, std::string::npos);
  BOOST_REQUIRE_NE(end, std::string::npos);
  const auto fitBody = source.substr(begin, end - begin);
  BOOST_CHECK_EQUAL(fitBody.find("SurfaceKind"), std::string::npos);
  BOOST_CHECK_EQUAL(fitBody.find("detectorId"), std::string::npos);
  BOOST_CHECK_EQUAL(fitBody.find("sourceId"), std::string::npos);

  BOOST_CHECK_EQUAL(source.find("fitLocalTripletWithMaterial"), std::string::npos);
  BOOST_CHECK_EQUAL(source.find("calculateMaterialPhysics"), std::string::npos);
  BOOST_CHECK_EQUAL(source.find("MaximumIterations"), std::string::npos);

  const auto trackerTraits = readFile(trackingRoot / "src/TrackerTraits.cxx");
  BOOST_CHECK_EQUAL(trackerTraits.find("fitLocalTripletUniformSolenoid"), std::string::npos);
  BOOST_CHECK_EQUAL(trackerTraits.find("TripletFitting.h"), std::string::npos);
}

BOOST_AUTO_TEST_CASE(CharacterizeStandaloneHostCost)
{
  const auto observations = makeHelixObservations();
  constexpr int Repetitions = 20000;
  double checksum = 0.;
  const auto start = std::chrono::steady_clock::now();
  for (int iteration = 0; iteration < Repetitions; ++iteration) {
    LocalTripletFitResult result{};
    BOOST_REQUIRE(fitLocalTripletUniformSolenoid(observations, {2.e-7}, result));
    checksum += result.curvature + result.chi2;
  }
  const auto elapsed = std::chrono::steady_clock::now() - start;
  const double nanosecondsPerFit =
    std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count() /
    static_cast<double>(Repetitions);
  BOOST_TEST_MESSAGE("local triplet fit host cost: " << nanosecondsPerFit << " ns/fit; checksum=" << checksum);
  BOOST_CHECK_GT(checksum, 0.);
}
