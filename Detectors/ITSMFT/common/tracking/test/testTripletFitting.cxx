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

std::array<TripletFitObservation, 4> makeAdjacentHelixObservations()
{
  const std::array<double, 4> angles{0.1, 0.16, 0.25, 0.33};
  std::array<TripletFitObservation, 4> observations{};
  for (std::size_t i = 0; i < observations.size(); ++i) {
    observations[i].position = {3. + Radius * std::cos(angles[i]),
                                -2. + Radius * std::sin(angles[i]),
                                1.5 + Radius * angles[i] * TanLambda};
    observations[i].covariance = makeCovariance();
  }
  return observations;
}

std::array<TripletFitFactor, 2> fitAdjacentFactors(
  const std::array<TripletFitObservation, 4>& observations)
{
  const std::array<TripletFitObservation, 3> first{
    observations[0], observations[1], observations[2]};
  const std::array<TripletFitObservation, 3> second{
    observations[1], observations[2], observations[3]};
  std::array<TripletFitFactor, 2> factors{};
  BOOST_REQUIRE(makeTripletFitFactor(first, factors[0]));
  BOOST_REQUIRE(makeTripletFitFactor(second, factors[1]));
  return factors;
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

double factorCovariance(const TripletFitFactor& factor,
                        const std::array<TripletFitObservation, 3>& observations,
                        bool leftTheta, bool rightTheta)
{
  double covariance = 0.;
  for (std::size_t hit = 0; hit < observations.size(); ++hit) {
    const auto& left = leftTheta ? factor.h[hit].theta : factor.h[hit].phi;
    const auto& right = rightTheta ? factor.h[hit].theta : factor.h[hit].phi;
    const auto& v = observations[hit].covariance;
    covariance += left[0] * (v.xx * right[0] + v.xy * right[1] + v.xz * right[2]) +
                  left[1] * (v.xy * right[0] + v.yy * right[1] + v.yz * right[2]) +
                  left[2] * (v.xz * right[0] + v.yz * right[1] + v.zz * right[2]);
  }
  return covariance;
}

std::string readFile(const std::filesystem::path& path)
{
  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path.string());
  return {std::istreambuf_iterator<char>{input}, {}};
}

} // namespace

BOOST_AUTO_TEST_CASE(GlobalObservationCopiesDecodedGeometryAndCovariance)
{
  GlobalMeasurement measurement{};
  measurement.position = {4.f, 3.f, -2.f};
  measurement.radius = 5.f;
  measurement.covariance = {4.f, -1.f, 1.5f, 9.f, -2.f, 16.f};

  TripletFitObservation observation{};
  BOOST_REQUIRE(makeTripletFitObservation(measurement, observation));
  BOOST_CHECK_EQUAL(observation.position[0], 4.);
  BOOST_CHECK_EQUAL(observation.position[1], 3.);
  BOOST_CHECK_EQUAL(observation.position[2], -2.);
  BOOST_CHECK_EQUAL(observation.covariance.xx, 4.);
  BOOST_CHECK_EQUAL(observation.covariance.xy, -1.);
  BOOST_CHECK_EQUAL(observation.covariance.xz, 1.5);
  BOOST_CHECK_EQUAL(observation.covariance.yy, 9.);
  BOOST_CHECK_EQUAL(observation.covariance.yz, -2.);
  BOOST_CHECK_EQUAL(observation.covariance.zz, 16.);
}

BOOST_AUTO_TEST_CASE(GlobalObservationAcceptsOnlyFloatRoundoffAtPSDBoundary)
{
  GlobalMeasurement measurement{};
  measurement.position = {1.f, 2.f, 3.f};
  measurement.covariance = {3.9999948e-10f, -1.9999946e-7f, 0.f,
                            9.999959e-5f, 0.f, 1.e-4f};
  const double minor = static_cast<double>(measurement.covariance.xx) * measurement.covariance.yy -
                       static_cast<double>(measurement.covariance.xy) * measurement.covariance.xy;
  BOOST_REQUIRE_LT(minor, 0.);

  TripletFitObservation observation{};
  BOOST_CHECK(makeTripletFitObservation(measurement, observation));

  measurement.covariance.xy *= 1.01f;
  BOOST_CHECK(!makeTripletFitObservation(measurement, observation));
}

BOOST_AUTO_TEST_CASE(ExactHelixProducesAConsistentFactor)
{
  const auto observations = makeHelixObservations();
  TripletFitFactor factor{};
  BOOST_REQUIRE(makeTripletFitFactor(observations, factor));
  BOOST_REQUIRE(factor.isValid());
  const double referenceCurvature = -static_cast<double>(factor.psi.phi) / factor.rho.phi;
  const double expectedCurvature = (1. / Radius) / std::sqrt(1. + TanLambda * TanLambda);
  checkClose(referenceCurvature, expectedCurvature, 4.e-4);
  BOOST_CHECK_SMALL(static_cast<double>(factor.psi.theta) +
                      static_cast<double>(factor.rho.theta) * referenceCurvature,
                    1.e-8);
  BOOST_CHECK_GT(factorCovariance(factor, observations, true, true), 0.);
  BOOST_CHECK_NE(factorCovariance(factor, observations, true, false), 0.);
}

BOOST_AUTO_TEST_CASE(AdjacentFactorsImplementEquation19ClosedForm)
{
  const SymmetricCovariance3D exact{};
  const std::array<TripletFitObservation, 4> observations{{
    {{{0., 0., 0.}}, exact},
    {{{1., 0., 0.}}, exact},
    {{{2., 0., 0.}}, exact},
    {{{3., 0., 0.}}, exact},
  }};
  std::array<TripletFitFactor, 2> factors{};
  factors[0].psi = {1.f, 2.f};
  factors[0].rho = {1.f, 1.f};
  factors[1].psi = {3.f, 4.f};
  factors[1].rho = {1.f, 1.f};

  AdjacentTripletFitResult result{};
  BOOST_REQUIRE(fitAdjacentTripletFactors(factors[0], factors[1], observations, {4., 9.}, result));
  const double rhoKpsi = 1. / 4. + 2. / 4. + 3. / 9. + 4. / 9.;
  const double rhoKrho = 1. / 4. + 1. / 4. + 1. / 9. + 1. / 9.;
  const double psiKpsi = 1. / 4. + 4. / 4. + 9. / 9. + 16. / 9.;
  checkClose(result.curvature, -rhoKpsi / rhoKrho, 1.e-14);
  checkClose(result.curvatureVariance, 1. / rhoKrho, 1.e-14);
  checkClose(result.chi2, psiKpsi - rhoKpsi * rhoKpsi / rhoKrho, 1.e-14);
}

BOOST_AUTO_TEST_CASE(AdjacentFactorsRetainSharedHitCrossCovariance)
{
  const SymmetricCovariance3D exact{};
  std::array<TripletFitObservation, 4> observations{{
    {{{0., 0., 0.}}, exact},
    {{{1., 0., 0.}}, {2., 0.5, 0., 3., 0., 0.}},
    {{{2., 0., 0.}}, exact},
    {{{3., 0., 0.}}, exact},
  }};
  std::array<TripletFitFactor, 2> factors{};
  factors[0].rho.phi = 1.f;
  factors[1].rho.phi = 1.f;
  factors[0].h[1].theta = {1.f, 2.f, 0.f};
  factors[0].h[1].phi = {-1.f, 1.f, 0.f};
  factors[1].h[0].theta = {3.f, -2.f, 0.f};
  factors[1].h[0].phi = {2.f, 4.f, 0.f};

  AdjacentTripletFitResult correlated{};
  BOOST_REQUIRE(fitAdjacentTripletFactors(factors[0], factors[1], observations,
                                          {100., 100.}, correlated));

  // Move the second factor's identical covariance contribution from shared
  // hit 1 to private hit 3. Diagonal blocks stay equal; only H V H^T's
  // cross-triplet block disappears.
  auto independentFactors = factors;
  auto independentObservations = observations;
  independentFactors[1].h[2] = independentFactors[1].h[0];
  independentFactors[1].h[0] = {};
  independentObservations[3].covariance = observations[1].covariance;
  AdjacentTripletFitResult independent{};
  BOOST_REQUIRE(fitAdjacentTripletFactors(independentFactors[0], independentFactors[1],
                                          independentObservations, {100., 100.}, independent));
  BOOST_CHECK_NE(correlated.curvatureVariance, independent.curvatureVariance);
}

BOOST_AUTO_TEST_CASE(AdjacentFactorsApplySpaceAngleMSGeometry)
{
  const SymmetricCovariance3D exact{};
  const std::array<TripletFitObservation, 4> observations{{
    {{{0., 0., 0.}}, exact},
    {{{1., 0., 1.}}, exact},
    {{{2., 0., 2.}}, exact},
    {{{3., 0., 3.}}, exact},
  }};
  std::array<TripletFitFactor, 2> factors{};
  factors[0].rho = {1.f, 1.f};
  factors[1].rho = {1.f, 1.f};
  AdjacentTripletFitResult result{};
  BOOST_REQUIRE(fitAdjacentTripletFactors(factors[0], factors[1], observations, {4., 9.}, result));
  const double expectedPrecision = 1. / 4. + 1. / 8. + 1. / 9. + 1. / 18.;
  checkClose(result.curvatureVariance, 1. / expectedPrecision, 1.e-14);
}

BOOST_AUTO_TEST_CASE(AdjacentExactHelixHasCommonCurvatureAndZeroQuality)
{
  const auto observations = makeAdjacentHelixObservations();
  const std::array<double, 2> angularVariance{1.e-8, 2.e-8};
  const auto factors = fitAdjacentFactors(observations);
  AdjacentTripletFitResult result{};
  BOOST_REQUIRE(fitAdjacentTripletFactors(factors[0], factors[1], observations, angularVariance, result));
  const double expectedCurvature = (1. / Radius) / std::sqrt(1. + TanLambda * TanLambda);
  checkClose(result.curvature, expectedCurvature, 4.e-4);
  // The persisted factors use floats; their round-trip leaves only this
  // numerical residue in an otherwise exactly common-curvature helix.
  BOOST_CHECK_SMALL(result.chi2, 1.e-9);
  BOOST_CHECK_GT(result.curvatureVariance, 0.);
}

BOOST_AUTO_TEST_CASE(AdjacentFactorFitIsRotationInvariant)
{
  const auto original = makeAdjacentHelixObservations();
  auto rotated = original;
  const double angle = 0.73;
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  for (auto& observation : rotated) {
    const double x = observation.position[0];
    const double y = observation.position[1];
    observation.position[0] = cosine * x - sine * y;
    observation.position[1] = sine * x + cosine * y;
    observation.covariance = rotateCovarianceAroundZ(observation.covariance, angle);
  }
  const std::array<double, 2> angularVariance{2.e-8, 3.e-8};
  const auto originalFactors = fitAdjacentFactors(original);
  const auto rotatedFactors = fitAdjacentFactors(rotated);
  AdjacentTripletFitResult first{};
  AdjacentTripletFitResult second{};
  BOOST_REQUIRE(fitAdjacentTripletFactors(originalFactors[0], originalFactors[1], original, angularVariance, first));
  BOOST_REQUIRE(fitAdjacentTripletFactors(rotatedFactors[0], rotatedFactors[1], rotated, angularVariance, second));
  checkClose(second.curvature, first.curvature, 2.e-7);
  checkClose(second.curvatureVariance, first.curvatureVariance, 2.e-6);
  checkClose(second.chi2, first.chi2, 1.e-5, 1.e-14);
}

BOOST_AUTO_TEST_CASE(AdjacentFactorInvalidInputsFailTransactionally)
{
  AdjacentTripletFitResult sentinel{};
  sentinel.curvature = 1.;
  sentinel.curvatureVariance = 2.;
  sentinel.chi2 = 3.;
  auto observations = makeAdjacentHelixObservations();
  const auto factors = fitAdjacentFactors(observations);
  AdjacentTripletFitResult result = sentinel;

  auto invalidFactors = factors;
  invalidFactors[0] = {};
  BOOST_CHECK(!fitAdjacentTripletFactors(invalidFactors[0], invalidFactors[1], observations, {1.e-8, 1.e-8}, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations[1].covariance.xy = 1.;
  BOOST_CHECK(!fitAdjacentTripletFactors(factors[0], factors[1], observations, {1.e-8, 1.e-8}, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations = makeAdjacentHelixObservations();
  BOOST_CHECK(!fitAdjacentTripletFactors(factors[0], factors[1], observations, {-1., 1.e-8}, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);
}

BOOST_AUTO_TEST_CASE(StraightTripletUsesTheRemovableZeroBendingLimit)
{
  const SymmetricCovariance3D covariance{1.e-6, 0., 0., 1.e-6, 0., 1.e-6};
  const std::array<TripletFitObservation, 3> observations{{
    {{1., 2., 3.}, covariance},
    {{2., 2., 3.5}, covariance},
    {{4., 2., 4.5}, covariance},
  }};
  TripletFitFactor factor{};
  BOOST_REQUIRE(makeTripletFitFactor(observations, factor));
  BOOST_REQUIRE(factor.isValid());
  BOOST_CHECK_SMALL(-static_cast<double>(factor.psi.phi) / factor.rho.phi, 1.e-15);
}

BOOST_AUTO_TEST_CASE(FactorConstructionFailsTransactionally)
{
  TripletFitFactor sentinel{};
  sentinel.psi = {1.f, 2.f};
  sentinel.rho = {3.f, 4.f};
  auto observations = makeHelixObservations();
  TripletFitFactor result = sentinel;

  observations[1].covariance.xy = 1.;
  BOOST_CHECK(!makeTripletFitFactor(observations, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations = makeHelixObservations();
  observations[2].position[0] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK(!makeTripletFitFactor(observations, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations = makeHelixObservations();
  observations[1].position = observations[0].position;
  BOOST_CHECK(!makeTripletFitFactor(observations, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);
}

BOOST_AUTO_TEST_CASE(FitAndTrackerTraitsCallSiteHaveNoFamilyDispatch)
{
  const auto trackingRoot = std::filesystem::path{__FILE__}.parent_path().parent_path();
  const auto source = readFile(trackingRoot / "src/TripletFitting.cxx");
  const auto begin = source.find("bool makeTripletFitFactor(");
  const auto end = source.find("bool fitAdjacentTripletFactors(", begin);
  BOOST_REQUIRE_NE(begin, std::string::npos);
  BOOST_REQUIRE_NE(end, std::string::npos);
  const auto fitBody = source.substr(begin, end - begin);
  BOOST_CHECK_EQUAL(fitBody.find("SurfaceKind"), std::string::npos);
  BOOST_CHECK_EQUAL(fitBody.find("detectorId"), std::string::npos);
  BOOST_CHECK_EQUAL(fitBody.find("sourceId"), std::string::npos);

  BOOST_CHECK_EQUAL(source.find("fitLocalTripletWithMaterial"), std::string::npos);
  BOOST_CHECK_EQUAL(source.find("calculateMaterialPhysics"), std::string::npos);
  BOOST_CHECK_EQUAL(source.find("MaximumIterations"), std::string::npos);
  BOOST_CHECK_EQUAL(source.find("TripletCompatibilityState"), std::string::npos);
  BOOST_CHECK_EQUAL(source.find("makeTripletCompatibilityState"), std::string::npos);
  BOOST_CHECK_EQUAL(source.find("fitLocalTripletUniformSolenoid"), std::string::npos);
  BOOST_CHECK_EQUAL(source.find("TripletSegmentEstimate"), std::string::npos);

  const auto trackerTraits = readFile(trackingRoot / "src/TrackerTraits.cxx");
  BOOST_CHECK_NE(trackerTraits.find("makeTripletFitFactor"), std::string::npos);
  BOOST_CHECK_NE(trackerTraits.find("fitAdjacentTripletFactors"), std::string::npos);
  BOOST_CHECK_NE(trackerTraits.find("TripletFitting.h"), std::string::npos);

  const auto candidateSource = readFile(trackingRoot / "src/CandidateFinding.cxx");
  const auto candidateHeader = readFile(trackingRoot / "include/ITSMFTTracking/detail/CandidateFinding.h");
  for (const auto retired : {"cellsCylinderAreCompatible", "cellsDiskAreCompatible"}) {
    BOOST_CHECK_EQUAL(candidateSource.find(retired), std::string::npos);
    BOOST_CHECK_EQUAL(candidateHeader.find(retired), std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(CharacterizeFactorConstructionHostCost)
{
  const auto observations = makeHelixObservations();
  constexpr int Repetitions = 20000;
  double checksum = 0.;
  const auto start = std::chrono::steady_clock::now();
  for (int iteration = 0; iteration < Repetitions; ++iteration) {
    TripletFitFactor factor{};
    BOOST_REQUIRE(makeTripletFitFactor(observations, factor));
    checksum += factor.psi.theta + factor.psi.phi + factor.rho.theta + factor.rho.phi;
  }
  const auto elapsed = std::chrono::steady_clock::now() - start;
  const double nanosecondsPerFit =
    std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count() /
    static_cast<double>(Repetitions);
  BOOST_TEST_MESSAGE("triplet-factor construction host cost: " << nanosecondsPerFit << " ns/factor; checksum=" << checksum);
  BOOST_CHECK_NE(checksum, 0.);
}

BOOST_AUTO_TEST_CASE(CharacterizeAdjacentFactorHostCost)
{
  const auto observations = makeAdjacentHelixObservations();
  const std::array<double, 2> angularVariance{2.e-7, 3.e-7};
  const auto factors = fitAdjacentFactors(observations);
  constexpr int Repetitions = 20000;
  double checksum = 0.;
  const auto start = std::chrono::steady_clock::now();
  for (int iteration = 0; iteration < Repetitions; ++iteration) {
    AdjacentTripletFitResult result{};
    BOOST_REQUIRE(fitAdjacentTripletFactors(factors[0], factors[1], observations, angularVariance, result));
    checksum += result.curvature + result.chi2;
  }
  const auto elapsed = std::chrono::steady_clock::now() - start;
  const double nanosecondsPerFit =
    std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count() /
    static_cast<double>(Repetitions);
  BOOST_TEST_MESSAGE("adjacent triplet-factor fit host cost: " << nanosecondsPerFit << " ns/fit; checksum=" << checksum);
  BOOST_CHECK_GT(checksum, 0.);
}
