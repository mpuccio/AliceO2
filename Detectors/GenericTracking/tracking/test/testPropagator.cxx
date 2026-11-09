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

/// \file testPropagator.cxx
/// \brief Closed-form helix propagation between sensor planes.
///
/// The budget under test: 1 um in position and 1 urad in
/// direction, against a 14 um single-pixel cluster RMS, held from 300 MeV to
/// 100 GeV and out to |z| = 152 cm.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Propagator
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/CovarianceTransport.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "Oracles.h"

#include "CommonConstants/PhysicsConstants.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include <cstring>
#include <iterator>
#include <limits>

namespace o2::gtrk
{

constexpr double kPosTol = 1e-4; // cm -- 1 um
constexpr double kAngTol = 1e-6; // rad

/// A plane, built by hand so a test does not depend on a layout's numbers.
static TrackingSurface planeAt(double x, double y, double z, double nxIn, double nyIn, double nzIn, SensorID id = 0)
{
  TrackingSurface s{};
  s.origin[0] = static_cast<Scalar>(x);
  s.origin[1] = static_cast<Scalar>(y);
  s.origin[2] = static_cast<Scalar>(z);
  s.normal[0] = static_cast<Scalar>(nxIn);
  s.normal[1] = static_cast<Scalar>(nyIn);
  s.normal[2] = static_cast<Scalar>(nzIn);
  // Any orthonormal pair spanning the plane will do for these tests.
  if (std::abs(nzIn) > 0.5) { // a disc
    s.axis0[0] = 1.F;
    s.axis1[1] = 1.F;
  } else { // a barrel-like plane: loc1 along the beam
    s.axis0[0] = static_cast<Scalar>(-nyIn);
    s.axis0[1] = static_cast<Scalar>(nxIn);
    s.axis1[2] = 1.F;
  }
  s.bounds = SurfaceBounds{.halfX = 1e6F, .halfY = 1e6F};
  s.id = id;
  return s;
}

/// A state sitting at the origin of \a s, heading at (phi, theta) with a given
/// signed q/p.
static BoundTrackState stateOn(const TrackingSurface& s, double phi, double theta, double qOverP)
{
  BoundTrackState st;
  st.par[kPhi] = static_cast<Scalar>(phi);
  st.par[kTheta] = static_cast<Scalar>(theta);
  st.par[kQOverP] = static_cast<Scalar>(qOverP);
  st.surface = s.id;
  return st;
}

/// A backward step at low pT, on a target plane rotated off radial.
///
/// The transverse solve is a quadratic and both roots are crossings of the
/// target's *infinite* plane. The root of smaller magnitude is the stable one
/// and is the sensor's crossing while the turn is small, but a soft track
/// turning tens of degrees between layers swaps them. Taking the small one then
/// would point the step *forward* out of a state asked to go back, and the
/// track would be rejected and lost: a selected track whose refit fails.
///
/// What holds: a backward step is never refused for running forwards. What
/// is not: the small root can also be a backward crossing of the extension
/// rather than the sensor's, which nothing reports -- the orientations near
/// -40 deg below are exactly that, and the assertions do not claim otherwise.
BOOST_AUTO_TEST_CASE(ABackwardStepAtLowMomentumIsNeverRefusedForRunningForwards)
{
  // R = 50 cm, so 0.3 GeV at 2 T: the circle through the origin with centre
  // (0, R) reaches radius r at a turn of 2 asin(r / 2R).
  constexpr double kR = 50.;
  const Propagator prop;
  const double kappa = 1. / kR;
  const auto qOverP = kappa / (static_cast<double>(o2::constants::math::B2C) * static_cast<double>(prop.bzKiloGauss));

  const auto at = [](double psi) { return std::pair<double, double>{kR * std::sin(psi), kR * (1. - std::cos(psi))}; };
  const auto turnAt = [](double r) { return 2. * std::asin(r / (2. * kR)); };

  const double psiFrom = turnAt(82.);
  const double psiTo = turnAt(66.);
  const auto from = at(psiFrom);
  const auto to = at(psiTo);
  const double truePath = (psiTo - psiFrom) / kappa; // negative: it is behind

  int reachedTheSensor = 0;
  for (int degrees = -80; degrees <= 80; degrees += 5) {
    // The plane always passes through the true crossing; only its orientation
    // about that point changes, and with it which root carries the crossing.
    const double radial = std::atan2(to.second, to.first) + (static_cast<double>(degrees) * M_PI / 180.);
    const auto source = planeAt(from.first, from.second, 0., std::cos(radial), std::sin(radial), 0., 1);
    const auto target = planeAt(to.first, to.second, 0., std::cos(radial), std::sin(radial), 0., 2);
    auto st = stateOn(source, psiFrom, M_PI / 2., qOverP);

    const auto res = propagate(prop, st, source, target, PropagationDirection::kBackward);
    BOOST_TEST_CONTEXT("plane rotated " << degrees << " deg off radial")
    {
      BOOST_REQUIRE(res.status == PropagationStatus::kOk || res.status == PropagationStatus::kOutOfBounds);
      BOOST_TEST(static_cast<double>(res.pathLength) < 0.); // it went back, as asked

      // Whichever root it is, it is a real crossing: on the target's plane by
      // construction of loc0/loc1, and on the circle of the state it started
      // from. loc0 runs around the beam, loc1 along it.
      const double x = static_cast<double>(target.origin[0]) + (static_cast<double>(st.loc0()) * static_cast<double>(target.axis0[0]));
      const double y = static_cast<double>(target.origin[1]) + (static_cast<double>(st.loc0()) * static_cast<double>(target.axis0[1]));
      BOOST_TEST(std::abs(std::hypot(x, y - kR) - kR) < 1e-2);

      reachedTheSensor += std::abs(static_cast<double>(res.pathLength) - truePath) < 1e-3 ? 1 : 0;
    }
  }
  // Most orientations put the sensor's own crossing on the chosen root; the
  // handful that do not are the open half described above.
  BOOST_TEST(reachedTheSensor >= 28);
}

BOOST_AUTO_TEST_CASE(SincAndVersineMatchTheirSeriesAndTheirClosedForm)
{
  // sinc never cancels, so the only thing to check is that the guarded branch
  // and the quotient agree across the crossover.
  for (const double x : {1e-6, 1e-5, 1e-4, 1.0001e-4, 1e-3, 0.1, 1.0, 3.0}) {
    const double direct = std::sin(x) / x;
    BOOST_TEST(std::abs(sinc<double>(x) - direct) < 1e-12);
    // versine = (1 - cos x)/x, evaluated without ever forming 1 - cos x.
    const double naive = (1. - std::cos(x)) / x;
    BOOST_TEST(std::abs(versine<double>(x) - naive) < 1e-9 * std::max(1., naive));
  }
  BOOST_TEST(std::abs(static_cast<double>(versine<float>(1e-3F)) - (1. - std::cos(1e-3)) / 1e-3) < 1e-9);
  BOOST_TEST(std::abs(sinc<double>(0.) - 1.) < 1e-15);
  BOOST_TEST(std::abs(versine<double>(0.)) < 1e-15);
}

BOOST_AUTO_TEST_CASE(AStraightTrackCrossesAPlaneWhereGeometrySaysItShould)
{
  // Zero field: the straight line, which the closed form reaches without a
  // separate code path.
  const Propagator prop{.bzKiloGauss = 0.F};

  const TrackingSurface start = planeAt(2., 0., 0., 1., 0., 0., 7);
  const TrackingSurface target = planeAt(10., 0., 0., 1., 0., 0., 9);

  // Heading out at 30 degrees in phi and 60 in theta.
  BoundTrackState st = stateOn(start, 30. * M_PI / 180., 60. * M_PI / 180., 1.);
  const auto res = propagate(prop, st, start, target, PropagationDirection::kForward);

  BOOST_TEST((res.status == PropagationStatus::kOk));
  BOOST_TEST(st.surface == 9U);

  // Straight line: the transverse step is 8 cm along x, so the path follows.
  const double sT = 8. / std::cos(30. * M_PI / 180.);
  const double path = sT / std::sin(60. * M_PI / 180.);
  BOOST_TEST(std::abs(static_cast<double>(res.pathLength) - path) < kPosTol);

  // Angles are untouched by a zero field.
  BOOST_TEST(std::abs(static_cast<double>(st.phi()) - 30. * M_PI / 180.) < kAngTol);
  BOOST_TEST(std::abs(static_cast<double>(st.theta()) - 60. * M_PI / 180.) < kAngTol);
}

BOOST_AUTO_TEST_CASE(ADiscIsReachedByADivision)
{
  const Propagator prop;
  const TrackingSurface start = planeAt(0., 0., 0., 0., 0., 1., 1);
  const TrackingSurface target = planeAt(0., 0., 152., 0., 0., 1., 2); // representative forward plane

  BoundTrackState st = stateOn(start, 0.4, 0.3, 1. / 5.);
  const auto res = propagate(prop, st, start, target, PropagationDirection::kForward);

  BOOST_TEST((res.status == PropagationStatus::kOk));
  BOOST_TEST(std::abs(static_cast<double>(res.pathLength) - 152. / std::cos(0.3)) < kPosTol);

  // The intersection lies on the plane: loc[2] is zero by construction, and
  // this is where that gets checked rather than at run time.
  Scalar g[3] = {};
  surfaceToGlobal(target, st.par[kLoc0], st.par[kLoc1], g);
  Scalar back[3] = {};
  surfaceToLocal(target, g, back);
  BOOST_TEST(std::abs(static_cast<double>(back[2])) < kPosTol);
}

BOOST_AUTO_TEST_CASE(ForwardThenBackwardReturnsTheOriginalState)
{
  // The round trip, over the momentum range the budget has to hold across
  // and for both charges.
  const Propagator prop;

  const TrackingSurface inner = planeAt(3.4, 0., 0., 1., 0., 0., 11);
  for (const double radius : {7., 17.2}) {
    const TrackingSurface outer = planeAt(radius, 0., 0., 1., 0., 0., 12);
    const std::vector<double> momenta = radius == 7. ? std::vector<double>{0.05, 0.1, 0.5} : std::vector<double>{0.3, 0.5, 1., 10., 100.};

    for (const double pT : momenta) {
      for (const double charge : {1., -1.}) {
        for (const double theta : {0.6, 1.2, 1.5708}) {
          const double qOverP = charge / (pT / std::sin(theta));

          BoundTrackState st = stateOn(inner, radius == 7. ? 0. : 0.05, theta, qOverP);
          const BoundTrackState original = st;

          const auto out = propagate(prop, st, inner, outer, PropagationDirection::kForward);
          BOOST_TEST((out.status == PropagationStatus::kOk));
          BOOST_TEST(out.pathLength > 0.F);

          const auto back = propagate(prop, st, outer, inner, PropagationDirection::kBackward);
          BOOST_TEST((back.status == PropagationStatus::kOk));
          BOOST_TEST(back.pathLength < 0.F);

          BOOST_TEST(std::abs(static_cast<double>(st.loc0() - original.loc0())) < kPosTol);
          BOOST_TEST(std::abs(static_cast<double>(st.loc1() - original.loc1())) < kPosTol);
          BOOST_TEST(std::abs(static_cast<double>(st.phi() - original.phi())) < kAngTol);
          BOOST_TEST(std::abs(static_cast<double>(st.theta() - original.theta())) < kAngTol);
          // The two path lengths cancel.
          BOOST_TEST(std::abs(static_cast<double>(out.pathLength + back.pathLength)) < kPosTol);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(AFiftyMeVTrackReachesLayerOneAndNoFurther)
{
  // Not an arbitrary threshold: at 50 MeV the radius is 8.34 cm, so from
  // r = 3.4 the trajectory turns back at x = 11.74 cm. Layer 1 at 7.0 is
  // inside that, layer 3 at 17.2 is not -- and the status has to say so
  // rather than return a nonsense intersection.
  const Propagator prop;

  const TrackingSurface l0 = planeAt(3.4, 0., 0., 1., 0., 0., 31);
  const TrackingSurface l1 = planeAt(7.0, 0., 0., 1., 0., 0., 32);
  const TrackingSurface l3 = planeAt(17.2, 0., 0., 1., 0., 0., 33);

  BoundTrackState reaches = stateOn(l0, 0., 1.5708, 1. / 0.05);
  BOOST_TEST((propagate(prop, reaches, l0, l1, PropagationDirection::kForward).status == PropagationStatus::kOk));

  BoundTrackState curls = stateOn(l0, 0., 1.5708, 1. / 0.05);
  const auto res = propagate(prop, curls, l0, l3, PropagationDirection::kForward);
  BOOST_TEST((res.status == PropagationStatus::kNoIntersection));
  // The state is untouched when the target is unreachable.
  BOOST_TEST(std::abs(static_cast<double>(curls.loc0())) < 1e-9);
}

BOOST_AUTO_TEST_CASE(ChargeSignsBendOppositeWays)
{
  const Propagator prop;
  const TrackingSurface inner = planeAt(3.4, 0., 0., 1., 0., 0., 1);
  const TrackingSurface outer = planeAt(17.2, 0., 0., 1., 0., 0., 2);

  BoundTrackState pos = stateOn(inner, 0., 1.5708, 1.0);
  BoundTrackState neg = stateOn(inner, 0., 1.5708, -1.0);

  BOOST_TEST((propagate(prop, pos, inner, outer, PropagationDirection::kForward).status == PropagationStatus::kOk));
  BOOST_TEST((propagate(prop, neg, inner, outer, PropagationDirection::kForward).status == PropagationStatus::kOk));

  // Equal and opposite deflection, and the sign comes from B2C rather than
  // from anything written here.
  BOOST_TEST(std::abs(static_cast<double>(pos.phi() + neg.phi())) < kAngTol);
  BOOST_TEST(static_cast<double>(pos.phi()) * static_cast<double>(neg.phi()) < 0.);
}

/// The same backward steps with the other crossing rules (`CrossingRule`),
/// the target a sensor of finite width: where
/// the small root is the extension's backward crossing, both rules take the
/// sensor's, and where the default already reaches the sensor they leave its
/// state bit for bit. The Jacobian of the step they change agrees with the
/// finite differences of the same rule.
BOOST_AUTO_TEST_CASE(TheOtherCrossingRulesTakeTheSensorsCrossing)
{
  constexpr double kR = 50.;
  const Propagator prop;
  const double kappa = 1. / kR;
  const auto qOverP = kappa / (static_cast<double>(o2::constants::math::B2C) * static_cast<double>(prop.bzKiloGauss));

  const auto at = [](double psi) { return std::pair<double, double>{kR * std::sin(psi), kR * (1. - std::cos(psi))}; };
  const auto turnAt = [](double r) { return 2. * std::asin(r / (2. * kR)); };

  const double psiFrom = turnAt(82.);
  const double psiTo = turnAt(66.);
  const auto from = at(psiFrom);
  const auto to = at(psiTo);
  const double truePath = (psiTo - psiFrom) / kappa;

  int changed = 0;
  for (int degrees = -80; degrees <= 80; degrees += 5) {
    const double radial = std::atan2(to.second, to.first) + (static_cast<double>(degrees) * M_PI / 180.);
    const auto source = planeAt(from.first, from.second, 0., std::cos(radial), std::sin(radial), 0., 1);
    auto target = planeAt(to.first, to.second, 0., std::cos(radial), std::sin(radial), 0., 2);
    target.bounds = SurfaceBounds{.halfX = 3.F, .halfY = 10.F};
    const auto st = stateOn(source, psiFrom, M_PI / 2., qOverP);

    auto first = st;
    auto nearer = st;
    auto inside = st;
    const auto rFirst = propagate(prop, first, source, target, PropagationDirection::kBackward);
    const auto rNearer = propagate<CrossingRule::kNearerCentre>(prop, nearer, source, target, PropagationDirection::kBackward);
    const auto rInside = propagate<CrossingRule::kInsideFirst>(prop, inside, source, target, PropagationDirection::kBackward);
    BOOST_TEST_CONTEXT("plane rotated " << degrees << " deg off radial")
    {
      BOOST_TEST((rNearer.status == PropagationStatus::kOk));
      BOOST_TEST((rInside.status == PropagationStatus::kOk));
      BOOST_TEST(std::abs(static_cast<double>(rNearer.pathLength) - truePath) < 1e-3);
      BOOST_TEST(std::abs(static_cast<double>(rInside.pathLength) - truePath) < 1e-3);
      if (std::abs(static_cast<double>(rFirst.pathLength) - truePath) < 1e-3) {
        // Bit identity is intentional; these records have no padding.
        // NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
        BOOST_TEST(std::memcmp(&first, &nearer, sizeof(BoundTrackState)) == 0);
        // NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
        BOOST_TEST(std::memcmp(&first, &inside, sizeof(BoundTrackState)) == 0);
        continue;
      }
      ++changed;
      // The step a rule changes is still a derivative the fit can take.
      Scalar step[kNBoundParams] = {};
      defaultDifferenceSteps(step);
      Scalar reference[kJacobianSize] = {};
      BOOST_TEST_REQUIRE((jacobianFiniteDifference<CrossingRule::kNearerCentre>(prop, st, source, target, PropagationDirection::kBackward, step, reference) == PropagationStatus::kOk));
      SolenoidJacobian solenoid;
      BOOST_TEST_REQUIRE((jacobianAnalytic<Scalar, CrossingRule::kNearerCentre>(prop, st, source, target, PropagationDirection::kBackward, solenoid) == PropagationStatus::kOk));
      Scalar analytic[kJacobianSize] = {};
      expandJacobian(solenoid, analytic);
      for (int i = 0; i < kJacobianSize; ++i) {
        BOOST_TEST(std::abs(static_cast<double>(analytic[i]) - static_cast<double>(reference[i])) < 1e-2 + 5e-3 * std::abs(static_cast<double>(reference[i])));
      }
    }
  }
  // The sweep reaches the case the rules are for.
  BOOST_TEST(changed >= 2);
}

BOOST_AUTO_TEST_CASE(AskingForTheWrongSideIsRefused)
{
  const Propagator prop;
  const TrackingSurface inner = planeAt(3.4, 0., 0., 1., 0., 0., 1);
  const TrackingSurface outer = planeAt(17.2, 0., 0., 1., 0., 0., 2);

  BoundTrackState st = stateOn(inner, 0., 1.5708, 1.);
  const auto res = propagate(prop, st, inner, outer, PropagationDirection::kBackward);
  BOOST_TEST((res.status == PropagationStatus::kWrongDirection));
  // ...and the state was not touched.
  BOOST_TEST(std::abs(static_cast<double>(st.loc0())) < 1e-9);
}

BOOST_AUTO_TEST_CASE(AGrazingIncidenceIsRejected)
{
  // A track running almost parallel to a disc: the bound representation is
  // singular there, so the covariance would be meaningless.
  const Propagator prop;
  const TrackingSurface start = planeAt(0., 0., 0., 0., 0., 1., 1);
  const TrackingSurface target = planeAt(0., 0., 10., 0., 0., 1., 2);

  BoundTrackState st = stateOn(start, 0., 1.5708 - 0.002, 1. / 10.); // theta ~ 90 deg
  const auto res = propagate(prop, st, start, target, PropagationDirection::kForward);
  BOOST_TEST((res.status == PropagationStatus::kGrazing));
}

BOOST_AUTO_TEST_CASE(LeavingTheActiveAreaIsReportedButStillPropagated)
{
  const Propagator prop;
  const TrackingSurface inner = planeAt(3.4, 0., 0., 1., 0., 0., 1);
  TrackingSurface outer = planeAt(17.2, 0., 0., 1., 0., 0., 2);
  outer.bounds = SurfaceBounds{.halfX = 0.84F, .halfY = 3.6F};

  BoundTrackState st = stateOn(inner, 0.5, 1.5708, 1. / 100.); // steep in phi, so it misses
  const auto res = propagate(prop, st, inner, outer, PropagationDirection::kForward);

  BOOST_TEST((res.status == PropagationStatus::kOutOfBounds));
  // Reported, not enforced: the intersection is there for a caller that wants it.
  BOOST_TEST(res.pathLength > 0.F);
  BOOST_TEST(std::abs(static_cast<double>(st.loc0())) > 0.84);
}

BOOST_AUTO_TEST_CASE(ThetaAndMomentumAreInvariantInAUniformField)
{
  // Not a tautology: it is why neither is recomputed in the propagation, and
  // a future field-map propagator will have to stop relying on it.
  const Propagator prop;
  const TrackingSurface inner = planeAt(3.4, 0., 0., 1., 0., 0., 1);
  const TrackingSurface outer = planeAt(17.2, 0., 0., 1., 0., 0., 2);

  BoundTrackState st = stateOn(inner, 0.2, 0.7, 1. / 3.);
  const BoundTrackState before = st;
  BOOST_TEST((propagate(prop, st, inner, outer, PropagationDirection::kForward).status == PropagationStatus::kOk));

  BOOST_TEST(std::abs(static_cast<double>(st.theta() - before.theta())) < 1e-9);
  BOOST_TEST(std::abs(static_cast<double>(st.qOverP() - before.qOverP())) < 1e-9);
}

BOOST_AUTO_TEST_CASE(AFlippedNormalChangesOnlyTheSignOfLoc0)
{
  // An analytic layout's discs all face +z and its staves away from the beam.
  // The MFT's faces alternate: our normal is -y of ALICE's local frame, which is +z on odd
  // faces and -z on even ones. A
  // plane with its normal and loc0 axis both reversed is the same plane, still
  // right-handed, so a track must cross it at the same point, after the same
  // path and time, with loc0 negated -- and only loc0's covariance terms may
  // change, by their sign. That holds for the target and for the source.
  const Propagator prop{.bzKiloGauss = 5.F};
  const auto flip = [](TrackingSurface s) {
    for (int i = 0; i < 3; ++i) {
      s.normal[i] = -s.normal[i];
      s.axis0[i] = -s.axis0[i];
    }
    return s;
  };
  struct Case {
    TrackingSurface from;
    TrackingSurface to;
    double theta;
    PropagationDirection dir;
  };
  const TrackingSurface disc0 = planeAt(0., 0., -45.29, 0., 0., 1., 1);
  const TrackingSurface disc1 = planeAt(0., 0., -53.81, 0., 0., 1., 2);
  const TrackingSurface stave0 = planeAt(2.33, 0., 0., 1., 0., 0., 3);
  const TrackingSurface stave1 = planeAt(19.59, 0., 0., 1., 0., 0., 4);
  const std::array<Case, 4> cases{{
    {.from = disc0, .to = disc1, .theta = 3.0, .dir = PropagationDirection::kForward},
    {.from = disc1, .to = disc0, .theta = 3.0, .dir = PropagationDirection::kBackward},
    {.from = stave0, .to = stave1, .theta = 1.2, .dir = PropagationDirection::kForward},
    {.from = stave1, .to = stave0, .theta = 1.2, .dir = PropagationDirection::kBackward},
  }};

  Scalar start[kBoundCovSize] = {};
  for (int i = 0; i < kNBoundParams; ++i) {
    start[covIndex(i, i)] = 1e-2F * static_cast<Scalar>(i + 1);
  }
  start[covIndex(kLoc0, kPhi)] = 2e-3F;
  start[covIndex(kLoc1, kTheta)] = 1e-3F;

  for (const auto& c : cases) {
    for (const double qOverP : {1. / 0.5, -1. / 0.5, 1. / 5.}) {
      for (const bool flipSource : {false, true}) {
        BoundTrackState st = stateOn(c.from, 0.3, c.theta, qOverP);
        st.par[kLoc0] = 0.1F;
        st.par[kLoc1] = -0.2F;
        const TrackingSurface from = flipSource ? flip(c.from) : c.from;
        BoundTrackState seen = st;
        Scalar seenCov[kBoundCovSize] = {};
        std::ranges::copy(start, std::begin(seenCov));
        if (flipSource) {
          // The same point and covariance, written in the flipped source.
          seen.par[kLoc0] = -seen.par[kLoc0];
          for (int j = 0; j < kNBoundParams; ++j) {
            if (j != kLoc0) {
              seenCov[covIndex(kLoc0, j)] = -seenCov[covIndex(kLoc0, j)];
            }
          }
        }

        BoundTrackState plain = seen;
        BoundTrackState flipped = seen;
        const TrackingSurface to = flip(c.to);
        const auto a = propagate(prop, plain, from, c.to, c.dir);
        const auto b = propagate(prop, flipped, from, to, c.dir);
        BOOST_TEST_REQUIRE((a.status == PropagationStatus::kOk));
        BOOST_TEST_REQUIRE((b.status == PropagationStatus::kOk));
        BOOST_TEST(std::abs(static_cast<double>(a.pathLength - b.pathLength)) < kPosTol);
        BOOST_TEST(std::abs(static_cast<double>(plain.loc0() + flipped.loc0())) < kPosTol);
        BOOST_TEST(std::abs(static_cast<double>(plain.loc1() - flipped.loc1())) < kPosTol);
        BOOST_TEST(std::abs(static_cast<double>(plain.phi() - flipped.phi())) < kAngTol);
        BOOST_TEST(std::abs(static_cast<double>(plain.theta() - flipped.theta())) < kAngTol);
        BOOST_TEST(std::abs(static_cast<double>(plain.qOverP() - flipped.qOverP())) < kAngTol);

        SolenoidJacobian ja{};
        SolenoidJacobian jb{};
        BOOST_TEST_REQUIRE((jacobianAnalytic(prop, seen, from, c.to, c.dir, ja) == PropagationStatus::kOk));
        BOOST_TEST_REQUIRE((jacobianAnalytic(prop, seen, from, to, c.dir, jb) == PropagationStatus::kOk));
        Scalar da[kJacobianSize] = {};
        Scalar db[kJacobianSize] = {};
        expandJacobian(ja, da);
        expandJacobian(jb, db);
        Scalar ca[kBoundCovSize] = {};
        Scalar cb[kBoundCovSize] = {};
        std::ranges::copy(seenCov, std::begin(ca));
        std::ranges::copy(seenCov, std::begin(cb));
        transportCovarianceDense(da, ca);
        transportCovarianceDense(db, cb);
        for (int i = 0; i < kNBoundParams; ++i) {
          for (int j = 0; j <= i; ++j) {
            const double sign = (i == kLoc0) != (j == kLoc0) ? -1. : 1.;
            const double expected = sign * static_cast<double>(ca[covIndex(i, j)]);
            const auto got = static_cast<double>(cb[covIndex(i, j)]);
            BOOST_TEST(std::abs(got - expected) <= 1e-4 * (1. + std::abs(expected)));
          }
        }
      }
    }
  }
}

/// A state with a NaN parameter is reported as kNotFinite and left untouched,
/// in every build: the finiteness test is GTRK's own, where O2's
/// `GPUCommonMath::Finite` is a constant outside the deterministic mode.
BOOST_AUTO_TEST_CASE(ANaNParameterIsReportedNotFinite)
{
  const Propagator prop;
  const auto from = planeAt(3., 0., 0., 1., 0., 0., 1);
  const auto barrel = planeAt(7., 0., 0., 1., 0., 0., 2);
  const auto disc = planeAt(0., 0., 12., 0., 0., 1., 3);
  for (const auto* to : {&barrel, &disc}) {
    for (const int index : {kPhi, kTheta, kQOverP}) {
      auto state = stateOn(from, 0.1, 1.0, 0.5);
      state.par[index] = std::numeric_limits<Scalar>::quiet_NaN();
      const auto before = state;
      const auto result = propagate(prop, state, from, *to, PropagationDirection::kForward);
      BOOST_TEST_CONTEXT("target " << to->id << " parameter " << index)
      {
        BOOST_TEST((result.status == PropagationStatus::kNotFinite));
        BOOST_TEST(state.surface == before.surface);
        BOOST_TEST(state.loc0() == before.loc0());
      }
    }
  }
}

} // namespace o2::gtrk
