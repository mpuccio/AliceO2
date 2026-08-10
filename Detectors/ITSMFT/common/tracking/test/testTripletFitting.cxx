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
  const std::array<TripletFitObservation, 4>& observations,
  const std::array<double, 2>& angularVariance)
{
  const std::array<TripletFitObservation, 3> first{
    observations[0], observations[1], observations[2]};
  const std::array<TripletFitObservation, 3> second{
    observations[1], observations[2], observations[3]};
  LocalTripletFitResult firstFit{};
  LocalTripletFitResult secondFit{};
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(first, {angularVariance[0]}, 5., firstFit));
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(second, {angularVariance[1]}, 5., secondFit));
  return {firstFit.factor, secondFit.factor};
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

BOOST_AUTO_TEST_CASE(ExactHelixHasExpectedCurvatureAndZeroLocalQuality)
{
  const auto observations = makeHelixObservations();
  LocalTripletFitResult result{};
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(observations, {1.e-8}, 5., result));

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
  checkClose(fittedTripletTransverseMomentum(result.inner), expectedPt, 4.e-4);
  BOOST_CHECK_GT(fittedTripletTransverseMomentumVariance(result.inner), 0.);
  checkClose(result.inner.tanLambda, TanLambda, 4.e-4);
  checkClose(result.outer.tanLambda, TanLambda, 4.e-4);
  BOOST_CHECK_GT(result.inner.covariance.qOverPtQOverPt, 0.);
  BOOST_CHECK_GT(result.inner.covariance.phiPhi, 0.);
  BOOST_CHECK_GT(result.inner.covariance.tanLambdaTanLambda, 0.);

  const auto& factor = result.factor;
  BOOST_REQUIRE(factor.isValid());
  checkClose(factor.psi.theta, result.thetaTilde, 1.e-7);
  checkClose(factor.psi.phi, result.phiTilde, 1.e-7);
  checkClose(factor.rho.theta, result.rhoTheta, 1.e-7);
  checkClose(factor.rho.phi, result.rhoPhi, 1.e-7);

  // Eq. (19): the persisted H reproduces the measurement part of the local
  // kink covariance, including every per-hit xyz cross term.
  checkClose(factorCovariance(factor, observations, true, true),
             result.gammaThetaTheta - 1.e-8, 2.e-6, 1.e-15);
  checkClose(factorCovariance(factor, observations, true, false),
             result.gammaThetaPhi, 2.e-6, 1.e-15);
  const double phiMSVariance = 1.e-8 / (result.referenceSinTheta * result.referenceSinTheta);
  checkClose(factorCovariance(factor, observations, false, false),
             result.gammaPhiPhi - phiMSVariance, 2.e-6, 1.e-15);
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
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(original, {2.e-7}, 5., first));
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(transformed, {2.e-7}, 5., second));
  checkClose(second.curvature, first.curvature, 2.e-13);
  checkClose(second.curvatureVariance, first.curvatureVariance, 2.e-12);
  checkClose(second.chi2, first.chi2, 1.e-9, 1.e-20);
  checkClose(second.gammaThetaPhi, first.gammaThetaPhi, 3.e-12);
  checkClose(second.inner.qOverPt, first.inner.qOverPt, 2.e-13);
  checkClose(second.inner.phi - first.inner.phi, angle, 2.e-13);
  checkClose(second.inner.tanLambda, first.inner.tanLambda, 2.e-13);
  checkClose(second.inner.covariance.qOverPtQOverPt,
             first.inner.covariance.qOverPtQOverPt, 3.e-11);
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
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(withCrossTerms, {3.e-8}, 5., full));
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(diagonal, {3.e-8}, 5., reduced));
  BOOST_CHECK_NE(full.gammaThetaPhi, reduced.gammaThetaPhi);
  BOOST_CHECK_NE(full.curvatureVariance, reduced.curvatureVariance);
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
  BOOST_REQUIRE(fitAdjacentTripletFactors(factors, observations, {4., 9.}, result));
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

  AdjacentTripletFitResult result{};
  BOOST_REQUIRE(fitAdjacentTripletFactors(factors, observations, {100., 100.}, result));
  const auto contraction = [&](const std::array<float, 3>& left,
                               const std::array<float, 3>& right) {
    const auto& v = observations[1].covariance;
    return left[0] * (v.xx * right[0] + v.xy * right[1]) +
           left[1] * (v.xy * right[0] + v.yy * right[1]);
  };
  checkClose(result.covariance[0][2],
             contraction(factors[0].h[1].theta, factors[1].h[0].theta), 1.e-14);
  checkClose(result.covariance[0][3],
             contraction(factors[0].h[1].theta, factors[1].h[0].phi), 1.e-14);
  checkClose(result.covariance[1][2],
             contraction(factors[0].h[1].phi, factors[1].h[0].theta), 1.e-14);
  checkClose(result.covariance[1][3],
             contraction(factors[0].h[1].phi, factors[1].h[0].phi), 1.e-14);
  BOOST_CHECK_EQUAL(result.covariance[0][2], result.covariance[2][0]);
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
  BOOST_REQUIRE(fitAdjacentTripletFactors(factors, observations, {4., 9.}, result));
  checkClose(result.covariance[0][0], 4., 1.e-14);
  checkClose(result.covariance[1][1], 8., 1.e-14);
  checkClose(result.covariance[2][2], 9., 1.e-14);
  checkClose(result.covariance[3][3], 18., 1.e-14);
}

BOOST_AUTO_TEST_CASE(AdjacentExactHelixHasCommonCurvatureAndZeroQuality)
{
  const auto observations = makeAdjacentHelixObservations();
  const std::array<double, 2> angularVariance{1.e-8, 2.e-8};
  const auto factors = fitAdjacentFactors(observations, angularVariance);
  AdjacentTripletFitResult result{};
  BOOST_REQUIRE(fitAdjacentTripletFactors(factors, observations, angularVariance, result));
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
  const auto originalFactors = fitAdjacentFactors(original, angularVariance);
  const auto rotatedFactors = fitAdjacentFactors(rotated, angularVariance);
  AdjacentTripletFitResult first{};
  AdjacentTripletFitResult second{};
  BOOST_REQUIRE(fitAdjacentTripletFactors(originalFactors, original, angularVariance, first));
  BOOST_REQUIRE(fitAdjacentTripletFactors(rotatedFactors, rotated, angularVariance, second));
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
  sentinel.covariance[0][0] = 4.;
  auto observations = makeAdjacentHelixObservations();
  const auto factors = fitAdjacentFactors(observations, {1.e-8, 1.e-8});
  AdjacentTripletFitResult result = sentinel;

  auto invalidFactors = factors;
  invalidFactors[0] = {};
  BOOST_CHECK(!fitAdjacentTripletFactors(invalidFactors, observations, {1.e-8, 1.e-8}, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations[1].covariance.xy = 1.;
  BOOST_CHECK(!fitAdjacentTripletFactors(factors, observations, {1.e-8, 1.e-8}, result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations = makeAdjacentHelixObservations();
  BOOST_CHECK(!fitAdjacentTripletFactors(factors, observations, {-1., 1.e-8}, result));
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
  LocalTripletFitResult result{};
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(observations, {1.e-8}, 5., result));
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
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(observations, {0.}, 5., measurementOnly));
  BOOST_REQUIRE(fitLocalTripletUniformSolenoid(observations, {2.e-4}, 5., withScattering));
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
  BOOST_CHECK(!fitLocalTripletUniformSolenoid(observations, {0.}, 5., result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations = makeHelixObservations();
  observations[2].position[0] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK(!fitLocalTripletUniformSolenoid(observations, {0.}, 5., result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations = makeHelixObservations();
  observations[1].position = observations[0].position;
  BOOST_CHECK(!fitLocalTripletUniformSolenoid(observations, {0.}, 5., result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations = makeHelixObservations();
  BOOST_CHECK(!fitLocalTripletUniformSolenoid(observations, {-1.}, 5., result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);

  observations = makeHelixObservations();
  BOOST_CHECK(!fitLocalTripletUniformSolenoid(observations, {0.}, 0., result));
  BOOST_CHECK_EQUAL(std::memcmp(&result, &sentinel, sizeof(result)), 0);
}

BOOST_AUTO_TEST_CASE(FitAndTrackerTraitsCallSiteHaveNoFamilyDispatch)
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
  BOOST_CHECK_EQUAL(source.find("TripletCompatibilityState"), std::string::npos);
  BOOST_CHECK_EQUAL(source.find("makeTripletCompatibilityState"), std::string::npos);

  const auto trackerTraits = readFile(trackingRoot / "src/TrackerTraits.cxx");
  BOOST_CHECK_NE(trackerTraits.find("fitLocalTripletUniformSolenoid"), std::string::npos);
  BOOST_CHECK_NE(trackerTraits.find("fitAdjacentTripletFactors"), std::string::npos);
  BOOST_CHECK_NE(trackerTraits.find("TripletFitting.h"), std::string::npos);

  const auto candidateSource = readFile(trackingRoot / "src/CandidateFinding.cxx");
  const auto candidateHeader = readFile(trackingRoot / "include/ITSMFTTracking/detail/CandidateFinding.h");
  for (const auto retired : {"cellsCylinderAreCompatible", "cellsDiskAreCompatible"}) {
    BOOST_CHECK_EQUAL(candidateSource.find(retired), std::string::npos);
    BOOST_CHECK_EQUAL(candidateHeader.find(retired), std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(CharacterizeStandaloneHostCost)
{
  const auto observations = makeHelixObservations();
  constexpr int Repetitions = 20000;
  double checksum = 0.;
  const auto start = std::chrono::steady_clock::now();
  for (int iteration = 0; iteration < Repetitions; ++iteration) {
    LocalTripletFitResult result{};
    BOOST_REQUIRE(fitLocalTripletUniformSolenoid(observations, {2.e-7}, 5., result));
    checksum += result.curvature + result.chi2;
  }
  const auto elapsed = std::chrono::steady_clock::now() - start;
  const double nanosecondsPerFit =
    std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count() /
    static_cast<double>(Repetitions);
  BOOST_TEST_MESSAGE("local triplet fit host cost: " << nanosecondsPerFit << " ns/fit; checksum=" << checksum);
  BOOST_CHECK_GT(checksum, 0.);
}

BOOST_AUTO_TEST_CASE(CharacterizeAdjacentFactorHostCost)
{
  const auto observations = makeAdjacentHelixObservations();
  const std::array<double, 2> angularVariance{2.e-7, 3.e-7};
  const auto factors = fitAdjacentFactors(observations, angularVariance);
  constexpr int Repetitions = 20000;
  double checksum = 0.;
  const auto start = std::chrono::steady_clock::now();
  for (int iteration = 0; iteration < Repetitions; ++iteration) {
    AdjacentTripletFitResult result{};
    BOOST_REQUIRE(fitAdjacentTripletFactors(factors, observations, angularVariance, result));
    checksum += result.curvature + result.chi2;
  }
  const auto elapsed = std::chrono::steady_clock::now() - start;
  const double nanosecondsPerFit =
    std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count() /
    static_cast<double>(Repetitions);
  BOOST_TEST_MESSAGE("adjacent triplet-factor fit host cost: " << nanosecondsPerFit << " ns/fit; checksum=" << checksum);
  BOOST_CHECK_GT(checksum, 0.);
}
