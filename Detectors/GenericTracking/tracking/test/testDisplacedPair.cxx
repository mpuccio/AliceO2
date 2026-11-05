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

/// \file testDisplacedPair.cxx
/// \brief The displaced pair test's geometry (DisplacedPair.h): displaced
///        helices' pairs on the ITS's radii pass at their own impact and z,
///        the test agrees with a scan of explicitly built circles, at an
///        impact bound of 0 it is the pair cut's transverse test, and the
///        degenerate cases are their own.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK displaced pair
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/DisplacedPair.h"
#include "GenericTracking/Tracking/PairCut.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <random>
#include <utility>

namespace o2::gtrk
{

namespace
{

using displaced::displacedPair;

constexpr std::array<double, 7> kRadii{2.3, 3.1, 3.9, 19.6, 24.5, 34.4, 39.3};
constexpr double kPi = std::numbers::pi;
constexpr double kInfinity = std::numeric_limits<double>::infinity();

/// The curvature of a track of \a ptGeV in 0.5 T [1/cm].
double curvatureAt(double ptGeV) { return 0.3 * 0.5 / (100. * ptGeV); }

/// A helix from (x0, y0, z0) in direction phi with signed curvature k
/// (positive counterclockwise) and dip slope tanLambda.
struct Helix {
  double x0, y0, z0, phi, k, tanLambda;
};

struct Crossing {
  double x, y, z;
};

/// The helix's first crossing of the radius \a r, if any.
std::optional<Crossing> crossing(const Helix& h, double r)
{
  const double radius = 1. / std::abs(h.k);
  const double sense = h.k > 0. ? 1. : -1.;
  // The centre to the left of the direction for a counterclockwise one.
  const double cx = h.x0 - (sense * radius * std::sin(h.phi));
  const double cy = h.y0 + (sense * radius * std::cos(h.phi));
  const double centre = std::hypot(cx, cy);
  const double cosine = ((r * r) - (centre * centre) - (radius * radius)) / (2. * radius * centre);
  if (std::abs(cosine) > 1.) {
    return std::nullopt;
  }
  const double start = std::atan2(h.y0 - cy, h.x0 - cx);
  const double toward = std::atan2(cy, cx);
  double best = kInfinity;
  for (const double side : {1., -1.}) {
    const double theta = toward + (side * std::acos(cosine));
    const double turn = std::fmod((sense * (theta - start)) + (4. * kPi), 2. * kPi);
    best = std::min(best, turn);
  }
  const double theta = start + (sense * best);
  return Crossing{.x = cx + (radius * std::cos(theta)), .y = cy + (radius * std::sin(theta)), .z = h.z0 + (h.tanLambda * radius * best)};
}

/// One circle through the pair, built explicitly: whether its arc from the
/// source to the target moves outward all along (sampled), its signed impact
/// (`|C| - R` counterclockwise, its negative clockwise, the line's distance
/// on its left), and the z of its closest approach, by the turn from the
/// arc's middle in (-pi, pi]; none for a circle centred on the beam line.
struct Circle {
  bool outward = false;
  double impact = 0.;
  std::optional<double> z;
};

Circle circleThrough(double x1, double y1, double z1, double x2, double y2, double z2, double k)
{
  const double mx = 0.5 * (x1 + x2);
  const double my = 0.5 * (y1 + y2);
  const double chord = std::hypot(x2 - x1, y2 - y1);
  const double ux = (x2 - x1) / chord;
  const double uy = (y2 - y1) / chord;
  const double h = 0.5 * chord;
  Circle circle;
  constexpr int kSamples = 32;
  auto outward = true;
  double previous = std::hypot(x1, y1);
  if (exactlyEqual(k, 0.)) {
    for (int i = 1; i <= kSamples; ++i) {
      const double t = chord * i / kSamples;
      const double r = std::hypot(x1 + (t * ux), y1 + (t * uy));
      outward = outward && r >= previous - 1.e-12;
      previous = r;
    }
    circle.outward = outward;
    circle.impact = (-uy * mx) + (ux * my);
    const double foot = -((x1 * ux) + (y1 * uy));
    circle.z = z1 + (foot * (z2 - z1) / chord);
    return circle;
  }
  const double radius = 1. / std::abs(k);
  const double sense = k > 0. ? 1. : -1.;
  const double offset = std::sqrt(std::max(0., (radius * radius) - (h * h)));
  const double cx = mx + (sense * offset * -uy);
  const double cy = my + (sense * offset * ux);
  const double theta1 = std::atan2(y1 - cy, x1 - cx);
  const double turn = std::remainder(sense * (std::atan2(y2 - cy, x2 - cx) - theta1), 2. * kPi);
  for (int i = 1; i <= kSamples; ++i) {
    const double theta = theta1 + (sense * turn * i / kSamples);
    const double r = std::hypot(cx + (radius * std::cos(theta)), cy + (radius * std::sin(theta)));
    outward = outward && r >= previous - 1.e-12;
    previous = r;
  }
  circle.outward = outward && turn > 0.;
  const double centre = std::hypot(cx, cy);
  circle.impact = sense * (centre - radius);
  if (centre > 0.) {
    const double middle = theta1 + (0.5 * sense * turn);
    const double phi = std::remainder(sense * (std::atan2(-cy, -cx) - middle), 2. * kPi);
    circle.z = (0.5 * (z1 + z2)) + (phi * (z2 - z1) / turn);
  }
  return circle;
}

/// The scan's verdict: whether some curvature on a fine grid of [-k_max,
/// k_max] gives an outward circle within the impact and z bounds; the z
/// range of the circles within the impact bound into \a range.
bool scanned(double x1, double y1, double z1, double x2, double y2, double z2, double maxCurvature, double maxImpact, double zLow, double zHigh,
             std::pair<double, double>* range = nullptr)
{
  const double h = 0.5 * std::hypot(x2 - x1, y2 - y1);
  const double bound = std::min(maxCurvature, 0.999999 / h);
  constexpr int kSteps = 4000;
  bool found = false;
  std::pair<double, double> zs{kInfinity, -kInfinity};
  for (int step = 0; step <= kSteps; ++step) {
    const double k = (step == kSteps / 2) ? 0. : -bound + (2. * bound * step / kSteps);
    const auto circle = circleThrough(x1, y1, z1, x2, y2, z2, k);
    if (!circle.outward || std::abs(circle.impact) > maxImpact || !circle.z) {
      continue;
    }
    zs = {std::min(zs.first, *circle.z), std::max(zs.second, *circle.z)};
    found = found || (*circle.z >= zLow && *circle.z <= zHigh);
  }
  if (range != nullptr) {
    *range = zs;
  }
  return found;
}

} // namespace

BOOST_AUTO_TEST_CASE(ADisplacedHelixsPairsPassAtTheirOwnImpactAndZ)
{
  // Displaced helices from 0.05 to 3 GeV, both charges, born anywhere inside
  // 18 cm; their pairs on adjacent layers and across one, tested at exactly
  // their own curvature as the bound, their own impact and z. Taken the other
  // way round -- target first -- no pair moves outward.
  std::mt19937_64 rng{29};
  std::uniform_real_distribution<double> unit{0., 1.};
  int pairs = 0;
  for (int trial = 0; trial < 20000; ++trial) {
    const double pt = 0.05 * std::exp(unit(rng) * std::log(60.));
    const double production = 18. * unit(rng);
    const double azimuth = 2. * kPi * unit(rng);
    const Helix h{.x0 = production * std::cos(azimuth),
                  .y0 = production * std::sin(azimuth),
                  .z0 = -10. + (20. * unit(rng)),
                  .phi = 2. * kPi * unit(rng),
                  .k = (unit(rng) < 0.5 ? -1. : 1.) * curvatureAt(pt),
                  .tanLambda = -1. + (2. * unit(rng))};
    std::size_t first = 0;
    while (first < kRadii.size() && kRadii[first] <= production) {
      ++first;
    }
    for (std::size_t layer = first; layer + 1U < kRadii.size(); ++layer) {
      for (std::size_t skip = 1; skip <= 2 && layer + skip < kRadii.size(); ++skip) {
        const auto a = crossing(h, kRadii[layer]);
        const auto b = crossing(h, kRadii[layer + skip]);
        if (!a || !b) {
          continue;
        }
        const auto own = circleThrough(a->x, a->y, a->z, b->x, b->y, b->z, h.k);
        if (!own.outward || !own.z) {
          continue; // a looper's arc past its turn: not an outward pair
        }
        ++pairs;
        const auto exact = displacedPair(a->x, a->y, a->z, b->x, b->y, b->z, std::abs(h.k), std::abs(own.impact) + 1.e-9, *own.z - 1.e-6, *own.z + 1.e-6);
        BOOST_TEST((exact.transverse && exact.longitudinal));
        const auto single =
          displacedPair<float>(static_cast<float>(a->x), static_cast<float>(a->y), static_cast<float>(a->z), static_cast<float>(b->x), static_cast<float>(b->y),
                               static_cast<float>(b->z), static_cast<float>(std::abs(h.k) * 1.001), static_cast<float>(std::abs(own.impact) + 1.e-3),
                               static_cast<float>(*own.z - 1.e-2), static_cast<float>(*own.z + 1.e-2));
        BOOST_TEST((single.transverse && single.longitudinal));
        BOOST_TEST(!displacedPair(b->x, b->y, b->z, a->x, a->y, a->z, std::abs(h.k), 100., -kInfinity, kInfinity).transverse);
      }
    }
  }
  BOOST_TEST(pairs > 30000);
}

BOOST_AUTO_TEST_CASE(TheTestAgreesWithAScanOfExplicitCircles)
{
  // Random pairs on the ITS's radii, random bounds, every third z bound put
  // just past or just short of the scan's z range: wherever the scan finds a
  // circle with a margin to spare the test passes, and wherever the test
  // passes the scan finds one with a margin's tolerance -- the scan's grid
  // (an impact step of at most ~0.01 cm) and the bisection's apart.
  std::mt19937_64 rng{70};
  std::uniform_real_distribution<double> unit{0., 1.};
  constexpr double kImpactMargin = 0.02; // [cm]
  constexpr double kZMargin = 0.1;       // [cm]
  int passed = 0;
  int failed = 0;
  int excludingZero = 0;
  for (int trial = 0; trial < 1500; ++trial) {
    const auto inner = static_cast<std::size_t>(unit(rng) * 6.);
    const std::size_t outer = std::min<std::size_t>(inner + 1U + static_cast<std::size_t>(unit(rng) * 2.), kRadii.size() - 1U);
    const double phi1 = 2. * kPi * unit(rng);
    const double phi2 = phi1 + std::normal_distribution<double>{0., 0.3}(rng);
    const double x1 = kRadii[inner] * std::cos(phi1);
    const double y1 = kRadii[inner] * std::sin(phi1);
    const double x2 = kRadii[outer] * std::cos(phi2);
    const double y2 = kRadii[outer] * std::sin(phi2);
    const double z1 = -10. + (20. * unit(rng));
    const double z2 = -15. + (30. * unit(rng));
    const double maxCurvature = curvatureAt(0.05 + unit(rng));
    const double maxImpact = 5. * unit(rng) * unit(rng);
    double zLow = -10. + (20. * unit(rng));
    double zHigh = zLow + (15. * unit(rng));
    std::pair<double, double> range;
    scanned(x1, y1, z1, x2, y2, z2, maxCurvature, maxImpact, -kInfinity, kInfinity, &range);
    if (trial % 3 == 0 && range.first <= range.second) {
      // A z bound ending just short of the range's least z or starting just past its most.
      const double gap = (unit(rng) - 0.5) * 4. * kZMargin;
      if (trial % 2 == 0) {
        zHigh = range.first + gap;
        zLow = zHigh - 5.;
      } else {
        zLow = range.second + gap;
        zHigh = zLow + 5.;
      }
    }
    if (!circleThrough(x1, y1, z1, x2, y2, z2, 0.).outward) {
      ++excludingZero;
    }
    const auto result = displacedPair(x1, y1, z1, x2, y2, z2, maxCurvature, maxImpact, zLow, zHigh);
    const bool both = result.transverse && result.longitudinal;
    if (scanned(x1, y1, z1, x2, y2, z2, maxCurvature, std::max(0., maxImpact - kImpactMargin), zLow + kZMargin, zHigh - kZMargin)) {
      BOOST_TEST(both);
    }
    if (both) {
      BOOST_TEST(scanned(x1, y1, z1, x2, y2, z2, maxCurvature, maxImpact + kImpactMargin, zLow - kZMargin, zHigh + kZMargin));
    }
    (both ? passed : failed) += 1;
  }
  // Both verdicts are exercised, and pairs whose outward curvatures leave out the line.
  BOOST_TEST(passed > 100);
  BOOST_TEST(failed > 100);
  BOOST_TEST(excludingZero > 10);
}

BOOST_AUTO_TEST_CASE(AtZeroImpactItIsThePairCutsTransverseTest)
{
  // The pair cut's transverse test without its errors, in a field: the
  // circle through the beam point at most k_max, no turning point between
  // the two. Pairs within a hair of either edge are left out: the two
  // compute them apart.
  std::mt19937_64 rng{6};
  std::uniform_real_distribution<double> unit{0., 1.};
  int compared = 0;
  int admitted = 0;
  for (int trial = 0; trial < 20000; ++trial) {
    const auto inner = static_cast<std::size_t>(unit(rng) * 6.);
    const std::size_t outer = inner + 1U;
    const double phi1 = 2. * kPi * unit(rng);
    const double phi2 = phi1 + std::normal_distribution<double>{0., 0.2}(rng);
    const PairPoint source{.x = static_cast<Scalar>(kRadii[inner] * std::cos(phi1)),
                           .y = static_cast<Scalar>(kRadii[inner] * std::sin(phi1)),
                           .z = 0.F,
                           .r = static_cast<Scalar>(kRadii[inner])};
    const PairPoint target{.x = static_cast<Scalar>(kRadii[outer] * std::cos(phi2)),
                           .y = static_cast<Scalar>(kRadii[outer] * std::sin(phi2)),
                           .z = 0.F,
                           .r = static_cast<Scalar>(kRadii[outer])};
    const double maxCurvature = curvatureAt(0.05 + unit(rng));
    PairCut cut{};
    cut.twoRMin = static_cast<Scalar>(2. / maxCurvature);
    cut.nSigma = 0.F;
    const auto terms = transversePairTerms(cut, VertexSeed{}, source, 0.F, target);
    const double sx = source.x;
    const double sy = source.y;
    const double sr = source.r;
    const double tx = target.x;
    const double ty = target.y;
    const double tr = target.r;
    const double cross = (sx * ty) - (tx * sy);
    const double chord = std::hypot(tx - sx, ty - sy);
    const double k0 = 2. * std::abs(cross) / (sr * tr * chord);
    const double dot = (sx * tx) + (sy * ty);
    const double low = sr * sr;
    const double high = tr * tr;
    if (std::abs(k0 - maxCurvature) < 1.e-4 * maxCurvature || std::abs(dot - low) < 1.e-4 * low || std::abs(dot - high) < 1.e-4 * high) {
      continue;
    }
    ++compared;
    const bool today = transverseCompatible(cut, terms);
    const auto result = displacedPair<double>(sx, sy, 0., tx, ty, 0., maxCurvature, 0., -kInfinity, kInfinity);
    BOOST_TEST(result.transverse == today);
    admitted += today ? 1 : 0;
  }
  BOOST_TEST(compared > 19000);
  BOOST_TEST(admitted > 1000);
  BOOST_TEST(admitted < compared - 1000);
}

BOOST_AUTO_TEST_CASE(TheDegenerateCasesAreTheirOwn)
{
  // Coincident points define no circle.
  BOOST_TEST(!displacedPair(3., 0., 0., 3., 0., 1., 0.1, 10., -kInfinity, kInfinity).transverse);
  // No field: the line alone, here 1 cm from the beam line, going outward.
  const auto line = displacedPair(2., 1., 0., 30., 1., 0., 0., 1.001, -kInfinity, kInfinity);
  BOOST_TEST((line.transverse && line.longitudinal));
  BOOST_TEST(!displacedPair(2., 1., 0., 30., 1., 0., 0., 0.999, -kInfinity, kInfinity).transverse);
  // The same line travelled inward does not go out from the source.
  BOOST_TEST(!displacedPair(30., 1., 0., 2., 1., 0., 0., 2., -kInfinity, kInfinity).transverse);
  // The line's closest approach at x = 0: z there from the two points' slope, 1 + 2 * (0 - 2) / 28.
  const double z = 1. + (2. * (0. - 2.) / 28.);
  BOOST_TEST(displacedPair(2., 1., 1., 30., 1., 3., 0., 2., z - 1.e-6, z + 1.e-6).longitudinal);
  BOOST_TEST(!displacedPair(2., 1., 1., 30., 1., 3., 0., 2., z + 1.e-3, z + 1.).longitudinal);
  // A line passing the beam line between the two: inward, then outward.
  BOOST_TEST(!displacedPair(-3., 0.5, 0., 3.9, 0.5, 0., 0.001, 10., -kInfinity, kInfinity).transverse);
}

BOOST_AUTO_TEST_CASE(SmearedPairsPassAtThreeSigmaOfTheirOwnImpactAndZ)
{
  // The displaced helices' pairs again, their measurements smeared by 5 um
  // across and along the beam: at 3 sigma of `displacedPairSigmas` they pass
  // at their own impact and z (a bound of zero width there) all but rarely,
  // and the sigmas bound the spread of the impact and z at the true
  // curvature -- pulls of an rms at most 1, not far below it.
  std::mt19937_64 rng{91};
  std::uniform_real_distribution<double> unit{0., 1.};
  std::normal_distribution<double> gauss{0., 1.};
  constexpr double kSigma = 5.e-4; // [cm]
  int pairs = 0;
  int passing = 0;
  double impactPull2 = 0.;
  double zPull2 = 0.;
  for (int trial = 0; trial < 20000; ++trial) {
    const double pt = 0.05 * std::exp(unit(rng) * std::log(60.));
    const double production = 18. * unit(rng);
    const double azimuth = 2. * kPi * unit(rng);
    const Helix h{.x0 = production * std::cos(azimuth),
                  .y0 = production * std::sin(azimuth),
                  .z0 = -10. + (20. * unit(rng)),
                  .phi = 2. * kPi * unit(rng),
                  .k = (unit(rng) < 0.5 ? -1. : 1.) * curvatureAt(pt),
                  .tanLambda = -1. + (2. * unit(rng))};
    std::size_t layer = 0;
    while (layer < kRadii.size() && kRadii[layer] <= production) {
      ++layer;
    }
    if (layer + 1U >= kRadii.size()) {
      continue;
    }
    const auto a = crossing(h, kRadii[layer]);
    const auto b = crossing(h, kRadii[layer + 1U]);
    if (!a || !b) {
      continue;
    }
    const auto own = circleThrough(a->x, a->y, a->z, b->x, b->y, b->z, h.k);
    if (!own.outward || !own.z) {
      continue;
    }
    ++pairs;
    // Smeared across the radius (r phi) and in z.
    const auto smear = [&](const Crossing& c) {
      const double r = std::hypot(c.x, c.y);
      const double d = kSigma * gauss(rng);
      return Crossing{.x = c.x - (d * c.y / r), .y = c.y + (d * c.x / r), .z = c.z + (kSigma * gauss(rng))};
    };
    const auto sa = smear(*a);
    const auto sb = smear(*b);
    displaced::PairErrors<double> errors;
    errors.sourceVarRPhi = errors.targetVarRPhi = errors.sourceVarZ = errors.targetVarZ = kSigma * kSigma;
    const auto result = displaced::displacedPairWithin(sa.x, sa.y, sa.z, sb.x, sb.y, sb.z, errors, 3., std::abs(h.k), std::abs(own.impact), *own.z, *own.z);
    passing += (result.transverse && result.longitudinal) ? 1 : 0;
    const auto sigmas = displaced::displacedPairSigmas(sa.x, sa.y, sa.z, sb.x, sb.y, sb.z, errors, std::abs(h.k), std::abs(own.impact));
    const auto measured = circleThrough(sa.x, sa.y, sa.z, sb.x, sb.y, sb.z, h.k);
    if (measured.z) {
      impactPull2 += std::pow((measured.impact - own.impact) / sigmas.impact, 2);
      zPull2 += std::pow((*measured.z - *own.z) / sigmas.z, 2);
    }
  }
  BOOST_REQUIRE(pairs > 10000);
  BOOST_TEST(passing >= static_cast<int>(0.997 * pairs));
  const double impactRms = std::sqrt(impactPull2 / pairs);
  const double zRms = std::sqrt(zPull2 / pairs);
  BOOST_TEST(impactRms <= 1.);
  BOOST_TEST(impactRms >= 0.2);
  BOOST_TEST(zRms <= 1.);
  BOOST_TEST(zRms >= 0.2);
  BOOST_TEST_MESSAGE("pull rms: impact " << impactRms << ", z " << zRms << "; passing " << passing << " of " << pairs);
}

BOOST_AUTO_TEST_CASE(TheScatteringAndTheBeamWidenTheImpactAsTheirVariances)
{
  displaced::PairErrors<double> errors;
  const auto none = displaced::displacedPairSigmas(2.3, 0.1, 0., 3.1, 0.2, 1., errors, 0.02, 0.5);
  BOOST_TEST(none.impact == 0.);
  BOOST_TEST(none.z == 0.);
  errors.sourceScatter = 4.e-6;
  errors.beamVarTransverse = 9.e-6;
  const auto both = displaced::displacedPairSigmas(2.3, 0.1, 0., 3.1, 0.2, 1., errors, 0.02, 0.5);
  BOOST_TEST(both.impact * both.impact == 13.e-6, boost::test_tools::tolerance(1.e-12));
  // The kink along the 3D path, the dip (1 / 0.806) squared over the chord's.
  const double chord2 = (0.8 * 0.8) + (0.1 * 0.1);
  BOOST_TEST(both.z * both.z == (1. + (1. / chord2)) * 4.e-6, boost::test_tools::tolerance(1.e-12));
}

BOOST_AUTO_TEST_CASE(TheReviewsCounterexamplesAreRight)
{
  // Counterexamples. (1) Near the circle's farthest point the arc's branch
  // must not hang on a vanishing product's sign:
  // every outward arc here has its closest approach before the source, at
  // z <= 0, so a bound of [1, 4] admits none -- in double and in float.
  BOOST_TEST(!displacedPair(34.4, 0., 0., 0., 39.3, 10., 0.03, 35., 1., 4.).longitudinal);
  BOOST_TEST(!displacedPair<float>(8.6197029F, 33.3025633F, 0.F, 39.296894F, 0.4940891F, 10.F, 0.03F, 30.F, 1.F, 4.F).longitudinal);
  // (2) On a line through the beam line, opposite sides: only the two
  // diameter circles move outward, each at impact 34.4 and closest-approach z 0.
  const auto diameter = displacedPair(-34.4, 0., 0., 39.3, 0., 10., 0.03, 35., -1., 1.);
  BOOST_TEST((diameter.transverse && diameter.longitudinal));
  BOOST_TEST(!displacedPair(-34.4, 0., 0., 39.3, 0., 10., 0.03, 34., -1., 1.).transverse);
  // (3) The circle centred on the beam line has no closest approach: no z
  // within a bound to admit.
  BOOST_TEST(!displacedPair(20.64, 27.52, 0., -20.64, 27.52, 0., 0.03, 35., 1., 2.).longitudinal);
}

} // namespace o2::gtrk
