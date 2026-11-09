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

/// \file testPropagationPrecision.cxx
/// \brief The 1 um propagation budget, checked against the helix itself.
///
/// testPropagator.cxx checks that propagation does the right thing. This checks
/// how well it does it, over a sweep wide enough to reach the conditioning the
/// detector can produce: grazing crossings, half turns, straight tracks and
/// |z| up to 300 cm, where a float's last bit is already 0.3 um.
///
/// The check is closure, not agreement with a second implementation of the
/// same inversion. Propagation reports where it landed and how long it took to
/// get there; the helix is then evaluated independently at that time of flight,
/// and the point it gives must lie on the target plane, at the local
/// coordinates that were stored, with the azimuth that was stored. Nothing is
/// inverted to check it, so the check itself is well conditioned everywhere the
/// solve is not. What it deliberately does not police is *which* crossing was
/// chosen: that is the root policy, and testPropagator.cxx owns it.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Propagation Precision
#include <boost/test/unit_test.hpp>

#include "PropagationPrecisionFixture.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <random>

namespace o2::gtrk
{
namespace
{
using namespace propagation_precision;
} // namespace

/// A million steps over the accepted geometry, both charges and both fields.
BOOST_AUTO_TEST_CASE(PropagationClosesOnItsOwnHelixWithinOneMicron)
{
  std::mt19937_64 rng(20260928);
  std::uniform_real_distribution<double> uniform(0., 1.);
  const auto sample = [&](double low, double high) { return low + (high - low) * uniform(rng); };
  constexpr int trials = 1000000;
  Worst worst;

  for (int trial = 0; trial < trials; ++trial) {
    const double phi = sample(-std::numbers::pi, std::numbers::pi);
    const double theta = sample(0.01, std::numbers::pi - 0.01);
    const double pt = std::exp(sample(std::log(0.03), std::log(10000.)));
    const double charge = trial % 2 == 0 ? 1. : -1.;
    const double radius = sample(0., 90.);
    const double azimuth = sample(-std::numbers::pi, std::numbers::pi);
    const auto from = plane(radius * std::cos(azimuth), radius * std::sin(azimuth), sample(-300., 300.), sample(-std::numbers::pi, std::numbers::pi), trial % 3 == 0, 1);

    BoundTrackState initial;
    initial.surface = from.id;
    initial.par[kLoc0] = static_cast<float>(sample(-1.5, 1.5));
    initial.par[kLoc1] = static_cast<float>(sample(-3.5, 3.5));
    initial.par[kPhi] = static_cast<float>(phi);
    initial.par[kTheta] = static_cast<float>(theta);
    initial.par[kQOverP] = static_cast<float>(charge * std::sin(theta) / pt);

    Propagator prop;
    prop.bzKiloGauss = trial % 5 == 0 ? -20.F : 20.F;
    if (trial % 17 == 0) {
      prop.bzKiloGauss = 0.F; // a straight track, and the second root at infinity
    }

    // Put the target on a point of the track, so that the sweep asks for
    // crossings that exist rather than for mostly misses.
    const double path = sample(0.001, 300.) * (trial % 7 == 0 ? -1. : 1.);
    const auto direction = path > 0. ? PropagationDirection::kForward : PropagationDirection::kBackward;
    const double turn = static_cast<double>(initial.qOverP()) * static_cast<double>(prop.bzKiloGauss) * static_cast<double>(constants::math::B2C) * path;
    if (std::abs(turn) > 3.1) {
      continue; // beyond a half turn the root policy, not the arithmetic, decides
    }
    ExactScalar origin[3]{};
    detail::boundToGlobal(from, static_cast<ExactScalar>(initial.loc0()), static_cast<ExactScalar>(initial.loc1()), origin);
    const double transverse = path * std::sin(static_cast<double>(initial.theta()));
    const double along = transverse * sinc(turn);
    const double across = transverse * versine(turn);
    const double endX = origin[0] + along * std::cos(static_cast<double>(initial.phi())) - across * std::sin(static_cast<double>(initial.phi()));
    const double endY = origin[1] + along * std::sin(static_cast<double>(initial.phi())) + across * std::cos(static_cast<double>(initial.phi()));
    const double endZ = origin[2] + path * std::cos(static_cast<double>(initial.theta()));
    auto to = plane(endX, endY, endZ, sample(-std::numbers::pi, std::numbers::pi), trial % 3 == 1, 2);
    // Slide the target's origin inside its own plane. The plane is unchanged,
    // but the crossing is no longer at its centre -- and the offset along the
    // line is what a grazing crossing amplifies, so a sweep that always lands
    // on the origin never asks the question.
    const double offset[2] = {sample(-5., 5.), sample(-5., 5.)};
    for (int i = 0; i < 3; ++i) {
      to.origin[i] += static_cast<float>(offset[0] * static_cast<double>(to.axis0[i]) + offset[1] * static_cast<double>(to.axis1[i]));
    }

    auto state = initial;
    const auto result = propagate(prop, state, from, to, direction);
    if (result.status != PropagationStatus::kOk && result.status != PropagationStatus::kOutOfBounds) {
      continue;
    }
    worst.add(closure(prop, initial, state, from, to, result.pathLength));
  }

  worst.report("sweep");
  BOOST_TEST(worst.checked > trials / 2);
  BOOST_TEST(worst.value.position < kMicron);
  BOOST_TEST(worst.value.plane < kMicron);
  BOOST_TEST(worst.value.azimuth < 1.e-6);
}

/// The conditioning the sweep only reaches by accident: either side of the
/// incidence cut, either side of a half turn, and a vanishing curvature, all
/// from a far corner of the detector.
BOOST_AUTO_TEST_CASE(GrazingHalfTurnAndStraightLimitsStayWithinOneMicron)
{
  Worst worst;
  const Propagator prop;
  const auto from = plane(80., -70., 290., 0.3, false, 1);

  for (const double phi : {-3.1, -1.7, -0.01, 0., 0.6, 2.7, 3.1}) {
    for (const double theta : {0.05, 0.3, 0.5 * std::numbers::pi, 2.9}) {
      for (const double turn : {-std::numbers::pi + 1.e-6, -2.8, -0.5, -1.e-6, 0., 1.e-9, 1.e-4, 0.499999, 0.500001, 2.8, std::numbers::pi - 1.e-6}) {
        for (const double incidence : {0.00999, 0.01001, 0.011, 0.02, 0.1, 0.7}) {
          if (incidence > std::sin(theta)) {
            continue;
          }
          // The offset of the crossing from the target's centre, which a
          // grazing crossing multiplies by up to a hundred.
          for (const double offset : {0., -4., 4.}) {
            BoundTrackState initial;
            initial.surface = from.id;
            initial.par[kPhi] = static_cast<float>(phi);
            initial.par[kTheta] = static_cast<float>(theta);
            constexpr double path = 100.;
            initial.par[kQOverP] = static_cast<float>(turn / (path * static_cast<double>(prop.bzKiloGauss) * static_cast<double>(constants::math::B2C)));

            // Place the target on the track again, now with its normal turned
            // away from the arrival direction by exactly the wanted incidence.
            const double actualTurn = static_cast<double>(initial.qOverP()) * static_cast<double>(prop.bzKiloGauss) * static_cast<double>(constants::math::B2C) * path;
            const double middle = static_cast<double>(initial.phi()) + 0.5 * actualTurn;
            const double chord = path * std::sin(static_cast<double>(initial.theta())) * sinc(0.5 * actualTurn);
            const double normal = static_cast<double>(initial.phi()) + actualTurn + std::acos(incidence / std::sin(theta));
            auto to = plane(80. + chord * std::cos(middle), -70. + chord * std::sin(middle), 290. + path * std::cos(static_cast<double>(initial.theta())), normal, false, 2);
            for (int i = 0; i < 3; ++i) {
              to.origin[i] += static_cast<float>(offset * static_cast<double>(to.axis0[i]));
            }

            auto state = initial;
            const auto result = propagate(prop, state, from, to, PropagationDirection::kForward);
            BOOST_TEST_CONTEXT("phi " << phi << ", theta " << theta << ", turn " << turn << ", incidence " << incidence << ", offset " << offset)
            {
              // Whether a crossing this shallow is reported at all, and which of
              // the two it is, is the root policy's business, not this test's.
              if (result.status == PropagationStatus::kOk || result.status == PropagationStatus::kOutOfBounds) {
                const auto next = closure(prop, initial, state, from, to, result.pathLength);
                worst.add(next);
                BOOST_TEST(next.position < kMicron);
                BOOST_TEST(next.plane < kMicron);
              }
            }
          }
        }
      }
    }
  }

  worst.report("boundary");
  BOOST_TEST(worst.checked > 1500);
  BOOST_TEST(worst.value.azimuth < 1.e-6);
}

} // namespace o2::gtrk
