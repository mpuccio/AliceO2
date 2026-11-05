// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file testSeedCuts.cxx
/// \brief Triplet and neighbour acceptance and rejection on generated helices.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK SeedCuts
#include <boost/test/unit_test.hpp>

#include "ToyBarrel4Disc8.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LayerMaterial.h"
#include "GenericTracking/Tracking/TripletCut.h"
#include "GenericTracking/Tracking/NeighbourCut.h"

#include "GeneratedHelixFixture.h"
#include "Oracles.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <random>

namespace o2::gtrk
{
namespace
{
using detectors::ToyBarrel4Disc8;
using Helix = GeneratedHelixFixture<ToyBarrel4Disc8>;
using State = Helix::State;
using Errors = Helix::Errors;
constexpr int kDraws = 2000;
constexpr double kKept = 0.99;
} // namespace

namespace triplets
{
namespace
{
struct Case {
  std::array<LayerID, 3> layers{};
  double slope = 0.;
  Errors errors;
};

struct Rates {
  int lambda = 0;
  int time = 0;
  int curvature = 0;
  int all = 0;
};

/// Generate kDraws triplets on the case's layers at \a ptGeV, with the
/// Highland kink at the middle layer at that momentum and Gaussian errors, and
/// count how many pass each cut of an iteration at \a minPtGeV. \a zShift and
/// \a timeShift move the third measurement.
Rates rates(const Case& c, double ptGeV, double minPtGeV, double zShift = 0., double timeShift = 0.)
{
  const auto parameters = Helix::parametersAt(minPtGeV);
  const double curvature = 1. / curvatureRadius(ptGeV, 20.);
  const auto secondLink = linkOf<ToyBarrel4Disc8>(c.layers[1], c.layers[2]);
  const TripletTargetBounds target{
    .varZ = static_cast<Scalar>(c.errors.varZ), .varR = static_cast<Scalar>(c.errors.varR), .varRPhi = static_cast<Scalar>(c.errors.varRPhi)};

  // A fixed seed: the rates are compared against fixed thresholds.
  std::mt19937_64 generator(20260926U); // NOLINT(cert-msc51-cpp)
  Rates result;
  for (int draw = 0; draw < kDraws; ++draw) {
    const State origin{.phi = 0.3, .slope = c.slope};
    const State a = Helix::toLayer(origin, curvature, c.layers[0]);
    const State b = Helix::toLayer(a, curvature, c.layers[1]);
    State third = Helix::toLayer(Helix::scatter(b, c.layers[1], ptGeV, parameters, generator), curvature, c.layers[2]);
    third.z += zShift;

    const auto first = Helix::measure(a, c.errors, 0., generator);
    const auto middle = Helix::measure(b, c.errors, 0., generator);
    const auto last = Helix::measure(third, c.errors, timeShift, generator);
    const auto lambdaFirst = std::atan2(middle.point.z - first.point.z, pairChord(first.point, middle.point).chord);
    const auto lambdaSecond = std::atan2(last.point.z - middle.point.z, pairChord(middle.point, last.point).chord);
    const Scalar halfWidth = tripletLambdaHalfWidth(parameters, secondLink, c.layers[1], first.point, middle.point, target);
    const bool lambdaPasses = std::abs(lambdaSecond - lambdaFirst) <= halfWidth;
    const bool timePasses = tripletTimeCompatible(first, last);
    const bool curvaturePasses = tripletCurvatureCompatible(parameters, c.layers[1], first.point, middle.point, last.point);
    TripletFitResult fit;
    const bool allPass = lambdaPasses && tripletCompatible(parameters, c.layers[1], first, middle, last, fit);
    result.lambda += lambdaPasses ? 1 : 0;
    result.time += timePasses ? 1 : 0;
    result.curvature += curvaturePasses ? 1 : 0;
    result.all += allPass ? 1 : 0;
  }
  return result;
}

const std::array<Case, 5> kCases{{
  {.layers = {0, 1, 2}, .slope = 0.3, .errors = Helix::kPixelErrors}, // barrel
  {.layers = {1, 2, 3}, .slope = 1., .errors = Helix::kPixelErrors},  // steep barrel
  {.layers = {1, 2, 3}, .slope = 0.6, .errors = Helix::kLongZErrors}, // barrel, poor z
  {.layers = {8, 9, 10}, .slope = 3., .errors = Helix::kDiscErrors},  // discs
  {.layers = {0, 8, 9}, .slope = 3., .errors = Helix::kPixelErrors},  // barrel into the discs
}};

} // namespace

BOOST_AUTO_TEST_CASE(TrueTripletsPassEveryCut)
{
  for (const auto& c : kCases) {
    for (const auto& [pt, minPt] : std::array<std::array<double, 2>, 4>{{{1., 1.}, {1.5, 1.}, {5., 1.}, {0.3, 0.3}}}) {
      const auto kept = rates(c, pt, minPt);
      BOOST_TEST_CONTEXT("layers " << int(c.layers[0]) << " " << int(c.layers[1]) << " " << int(c.layers[2]) << " pT " << pt << " min " << minPt)
      {
        BOOST_TEST(kept.lambda >= kKept * kDraws);
        BOOST_TEST(kept.time == kDraws);
        BOOST_TEST(kept.curvature >= kKept * kDraws);
        BOOST_TEST(kept.all >= kKept * kDraws);
      }
    }
  }
}

namespace
{

/// Triplets whose middle point kinked by \a kinkTheta0 times the Highland
/// theta0 at \a ptGeV (in the dip and in azimuth alike), under an iteration
/// at \a minPtGeV with the kink allowance \a kinkSigma: how many pass the
/// nominal cuts, how many only with the allowance.
struct KinkedRates {
  int nominal = 0;
  int kinked = 0;
};

KinkedRates kinkedRates(const Case& c, double ptGeV, double minPtGeV, double kinkTheta0, double kinkSigma)
{
  const auto parameters = Helix::parametersAt(minPtGeV, kinkSigma);
  const double curvature = 1. / curvatureRadius(ptGeV, 20.);
  std::mt19937_64 generator(20261004U); // NOLINT(cert-msc51-cpp)
  KinkedRates result;
  for (int draw = 0; draw < kDraws; ++draw) {
    const State origin{.phi = 0.3, .slope = c.slope};
    const State a = Helix::toLayer(origin, curvature, c.layers[0]);
    const State b = Helix::toLayer(a, curvature, c.layers[1]);
    // The hard kink: theta0 of the middle layer at the track's momentum, as
    // Helix::scatter draws it, times kinkTheta0, on both angles.
    const double cosine = 1. / std::sqrt(1. + (b.slope * b.slope));
    const auto k = static_cast<double>(parameters.scattering.layers[c.layers[1]].k);
    const double theta0 = std::sqrt(k / cosine) * cosine / ptGeV;
    State kinkedState = b;
    const double sign = (draw % 2 == 0) ? 1. : -1.;
    kinkedState.slope = std::tan(std::atan(b.slope) + (sign * kinkTheta0 * theta0));
    kinkedState.phi += sign * kinkTheta0 * theta0 / cosine;
    const State third = Helix::toLayer(kinkedState, curvature, c.layers[2]);
    const auto first = Helix::measure(a, c.errors, 0., generator);
    const auto middle = Helix::measure(b, c.errors, 0., generator);
    const auto last = Helix::measure(third, c.errors, 0., generator);
    TripletFitResult fit;
    KinkTerms terms;
    bool kinked = false;
    if (tripletCompatible(parameters, c.layers[1], first, middle, last, fit, terms, kinked)) {
      result.nominal += kinked ? 0 : 1;
      result.kinked += kinked ? 1 : 0;
    }
  }
  return result;
}

} // namespace

BOOST_AUTO_TEST_CASE(TheKinkAllowanceKeepsAHardKinkAndNothingElse)
{
  // A 4 theta0 kink at the middle point at pT_min, on both angles: the
  // nominal cuts keep few of them (the curvature test is one-sided, a kink
  // can straighten the circle), an allowance of 2 keeps the rest; without
  // one nothing more survives, and the nominal count is the same either way.
  const auto without = kinkedRates(kCases[0], 1., 1., 4., 1.);
  BOOST_TEST(without.nominal <= 0.15 * kDraws);
  BOOST_TEST(without.kinked == 0);
  const auto with = kinkedRates(kCases[0], 1., 1., 4., 2.);
  BOOST_TEST(with.nominal == without.nominal);
  BOOST_TEST(with.nominal + with.kinked >= 0.95 * kDraws);
  // An 8 theta0 kink is beyond an allowance of 2 as well.
  const auto beyond = kinkedRates(kCases[0], 1., 1., 8., 2.);
  BOOST_TEST(beyond.kinked <= 0.02 * kDraws);
  // True triplets without a kink are never flagged: the nominal cuts come first.
  const auto straight = kinkedRates(kCases[0], 1., 1., 0., 2.);
  BOOST_TEST(straight.nominal >= kKept * kDraws);
  BOOST_TEST(straight.kinked <= 0.01 * kDraws);
}

BOOST_AUTO_TEST_CASE(ASofterTrackFailsTheCurvatureCut)
{
  const auto kept = rates(kCases[0], 0.5, 1.);
  BOOST_TEST(kept.curvature <= 0.05 * kDraws);
}

BOOST_AUTO_TEST_CASE(AThirdMeasurementOffTheHelixFailsTheDipAngleCut)
{
  const auto kept = rates(kCases[0], 1., 1., 0.1);
  BOOST_TEST(kept.lambda <= 0.05 * kDraws);
}

BOOST_AUTO_TEST_CASE(AThirdMeasurementTicksLateFailsTheTimeCut)
{
  const auto kept = rates(kCases[0], 1., 1., 0., 75.);
  BOOST_TEST(kept.time == 0);
}

} // namespace triplets

namespace neighbours
{
namespace
{
struct Case {
  std::array<LayerID, 4> layers{};
  double slope = 0.;
  Errors errors;
};

/// How a draw departs from one helix after `c`.
struct Departure {
  double softPtGeV = 0.; ///< `d` from a track of this pT through `c`, if set
  double zShift = 0.;    ///< `d` moved along z [cm]
  double timeShift = 0.;
};

struct Rates {
  int curvature = 0;
  int time = 0;
  int fit = 0;
  int all = 0;
  double chi2Sum = 0.; ///< over the draws the adjacent fit converged on
  int fitted = 0;
};

/// Generate kDraws four-measurement chains on the case's layers at \a ptGeV,
/// with the Highland kinks at `b` and `c` and Gaussian errors, and count how
/// many pass each cut of an iteration at \a minPtGeV.
Rates rates(const Case& c, double ptGeV, double minPtGeV, const Departure& departure = {})
{
  const auto parameters = Helix::parametersAt(minPtGeV);
  const double curvature = 1. / curvatureRadius(ptGeV, 20.);

  // A fixed seed: the rates are compared against fixed thresholds.
  std::mt19937_64 generator(20260926U); // NOLINT(cert-msc51-cpp)
  Rates result;
  for (int draw = 0; draw < kDraws; ++draw) {
    const State origin{.phi = 0.3, .slope = c.slope};
    const State a = Helix::toLayer(origin, curvature, c.layers[0]);
    const State b = Helix::toLayer(a, curvature, c.layers[1]);
    const State third = Helix::toLayer(Helix::scatter(b, c.layers[1], ptGeV, parameters, generator), curvature, c.layers[2]);
    const double lastCurvature = departure.softPtGeV > 0. ? 1. / curvatureRadius(departure.softPtGeV, 20.) : curvature;
    State fourth = Helix::toLayer(Helix::scatter(third, c.layers[2], ptGeV, parameters, generator), lastCurvature, c.layers[3]);
    fourth.z += departure.zShift;

    const std::array<TimedPoint, 4> m{Helix::measure(a, c.errors, 0., generator), Helix::measure(b, c.errors, 0., generator),
                                      Helix::measure(third, c.errors, 0., generator), Helix::measure(fourth, c.errors, departure.timeShift, generator)};
    TripletFitResult first;
    TripletFitResult second;
    const bool fitted = fitTriplet(m[0].point, m[1].point, m[2].point, c.layers[1], parameters.scattering, parameters.pairCut, first) &&
                        fitTriplet(m[1].point, m[2].point, m[3].point, c.layers[2], parameters.scattering, parameters.pairCut, second);
    const bool curvaturePasses = fitted && neighbourCurvatureCompatible(parameters, first.curvature, first.varCurvature, second.curvature, second.varCurvature);
    const bool timePasses = neighbourTimeCompatible(m[0], m[3]);
    TripletFitResult adjacent;
    const bool adjacentFitted =
      fitAdjacentTriplets(m[0].point, m[1].point, m[2].point, m[3].point, c.layers[1], c.layers[2], parameters.scattering, parameters.pairCut, adjacent);
    TripletFitResult fit;
    const bool allPass = curvaturePasses && neighbourCompatible(parameters, c.layers[1], c.layers[2], m[0], m[1], m[2], m[3], fit);
    result.curvature += curvaturePasses ? 1 : 0;
    result.time += timePasses ? 1 : 0;
    result.fit += adjacentFitted && adjacent.chi2 <= parameters.neighbourChi2Max ? 1 : 0;
    result.all += allPass ? 1 : 0;
    if (adjacentFitted) {
      result.chi2Sum += static_cast<double>(adjacent.chi2);
      ++result.fitted;
    }
  }
  return result;
}

const std::array<Case, 5> kCases{{
  {.layers = {0, 1, 2, 3}, .slope = 0.3, .errors = Helix::kPixelErrors}, // barrel
  {.layers = {0, 1, 2, 3}, .slope = 1., .errors = Helix::kPixelErrors},  // steep barrel
  {.layers = {0, 1, 2, 3}, .slope = 0.6, .errors = Helix::kLongZErrors}, // barrel, poor z
  {.layers = {8, 9, 10, 11}, .slope = 3., .errors = Helix::kDiscErrors}, // discs
  {.layers = {0, 8, 9, 10}, .slope = 3., .errors = Helix::kPixelErrors}, // barrel into the discs
}};

} // namespace

BOOST_AUTO_TEST_CASE(TrueNeighboursPassEveryCut)
{
  for (const auto& c : kCases) {
    for (const auto& [pt, minPt] : std::array<std::array<double, 2>, 4>{{{1., 1.}, {1.5, 1.}, {5., 1.}, {0.3, 0.3}}}) {
      const auto kept = rates(c, pt, minPt);
      BOOST_TEST_CONTEXT("layers " << int(c.layers[0]) << " " << int(c.layers[1]) << " " << int(c.layers[2]) << " " << int(c.layers[3]) << " slope " << c.slope
                                   << " pT " << pt << " min " << minPt)
      {
        BOOST_TEST(kept.curvature >= kKept * kDraws);
        BOOST_TEST(kept.time == kDraws);
        BOOST_TEST(kept.fit >= kKept * kDraws);
        BOOST_TEST(kept.all >= kKept * kDraws);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TheAdjacentChi2HasThreeDegreesOfFreedomAtPtMin)
{
  // At pT_min the fit's scattering is the generator's, so its chi2 follows
  // the 3-dof distribution: mean 3, and that mean's error over kDraws is
  // sqrt(6 / kDraws) = 0.055.
  for (const auto& c : kCases) {
    const auto kept = rates(c, 1., 1.);
    const double mean = kept.chi2Sum / kept.fitted;
    BOOST_TEST_CONTEXT("layers " << int(c.layers[0]) << " " << int(c.layers[1]) << " " << int(c.layers[2]) << " " << int(c.layers[3]) << " slope " << c.slope)
    {
      BOOST_TEST_MESSAGE("mean chi2 " << mean);
      BOOST_TEST(kept.fitted == kDraws);
      BOOST_TEST(mean >= 2.7);
      BOOST_TEST(mean <= 3.3);
    }
  }
}

BOOST_AUTO_TEST_CASE(ASofterTrackThroughTheSharedPairFailsTheCurvatureCut)
{
  const auto kept = rates(kCases[0], 5., 1., {.softPtGeV = 0.5});
  BOOST_TEST(kept.curvature <= 0.05 * kDraws);
}

BOOST_AUTO_TEST_CASE(AFourthMeasurementOffTheHelixFailsTheFit)
{
  const auto kept = rates(kCases[0], 1., 1., {.zShift = 0.1});
  BOOST_TEST(kept.fit <= 0.05 * kDraws);
}

BOOST_AUTO_TEST_CASE(AFourthMeasurementTicksLateFailsTheTimeCut)
{
  const auto kept = rates(kCases[0], 1., 1., {.timeShift = 75.});
  BOOST_TEST(kept.time == 0);
}

} // namespace neighbours
} // namespace o2::gtrk
