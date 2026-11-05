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

/// \file testO2Parameters.cxx
/// \brief A bound state as O2's parameters: the same point, momentum and
///        charge, and the covariance's Jacobian against finite differences of
///        the propagator onto O2's plane -- on tilted planes, in every azimuth,
///        moving outward and inward, at small sin(theta).

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK O2Parameters
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Workflow/O2Parameters.h"

#include <array>
#include <cmath>
#include <numbers>
#include <random>

namespace o2::gtrk::workflow
{
namespace
{

constexpr double kBz = 5.; // the ITS's and the MFT's 0.5 T [kG]

void cross(const double a[3], const double b[3], double out[3])
{
  out[0] = (a[1] * b[2]) - (a[2] * b[1]);
  out[1] = (a[2] * b[0]) - (a[0] * b[2]);
  out[2] = (a[0] * b[1]) - (a[1] * b[0]);
}

/// A plane through \a origin with normal \a n (unit) and loc1 along the
/// component of \a up orthogonal to it.
TrackingSurface plane(const std::array<double, 3>& origin, const double n[3], const double up[3])
{
  const double dot = (up[0] * n[0]) + (up[1] * n[1]) + (up[2] * n[2]);
  double a1[3] = {up[0] - (dot * n[0]), up[1] - (dot * n[1]), up[2] - (dot * n[2])};
  const double norm = std::hypot(a1[0], a1[1], a1[2]);
  for (auto& v : a1) {
    v /= norm;
  }
  double a0[3];
  cross(a1, n, a0); // axis0 x axis1 = normal
  TrackingSurface s{};
  for (int i = 0; i < 3; ++i) {
    s.origin[i] = static_cast<Scalar>(origin[i]);
    s.axis0[i] = static_cast<Scalar>(a0[i]);
    s.axis1[i] = static_cast<Scalar>(a1[i]);
    s.normal[i] = static_cast<Scalar>(n[i]);
  }
  s.bounds = SurfaceBounds{.halfX = 1e6F, .halfY = 1e6F};
  s.id = 1;
  return s;
}

/// A barrel sensor at azimuth \a alpha, radius \a r, its normal tilted by
/// \a tiltPhi in azimuth and \a tiltZ towards z: alignment.
TrackingSurface barrelSensor(double alpha, double r, double z, double tiltPhi, double tiltZ)
{
  const double n[3] = {std::cos(alpha + tiltPhi) * std::cos(tiltZ), std::sin(alpha + tiltPhi) * std::cos(tiltZ), std::sin(tiltZ)};
  const double up[3] = {0., 0., 1.};
  return plane({r * std::cos(alpha), r * std::sin(alpha), z}, n, up);
}

/// An MFT-like sensor at \a z, its normal +z or -z, tilted by \a tilt.
TrackingSurface discSensor(double x, double y, double z, double sign, double tilt)
{
  const double n[3] = {std::sin(tilt), 0., sign * std::cos(tilt)};
  const double up[3] = {0., 1., 0.};
  return plane({x, y, z}, n, up);
}

BoundTrackState state(double loc0, double loc1, double phi, double theta, double qOverP)
{
  BoundTrackState s;
  s.par[kLoc0] = static_cast<Scalar>(loc0);
  s.par[kLoc1] = static_cast<Scalar>(loc1);
  s.par[kPhi] = static_cast<Scalar>(phi);
  s.par[kTheta] = static_cast<Scalar>(theta);
  s.par[kQOverP] = static_cast<Scalar>(qOverP);
  for (int i = 0; i < kNBoundParams; ++i) {
    s.variance(i) = 1e-4F;
  }
  s.surface = 1;
  return s;
}

void globalOf(const BoundTrackState& s, const TrackingSurface& p, double r[3], double mom[3])
{
  const double pAbs = 1. / std::abs(static_cast<double>(s.qOverP()));
  for (int i = 0; i < 3; ++i) {
    r[i] = static_cast<double>(p.origin[i]) + (static_cast<double>(s.loc0()) * p.axis0[i]) + (static_cast<double>(s.loc1()) * p.axis1[i]);
  }
  mom[0] = pAbs * std::cos(s.phi()) * std::sin(s.theta());
  mom[1] = pAbs * std::sin(s.phi()) * std::sin(s.theta());
  mom[2] = pAbs * std::cos(s.theta());
}

double cosIncidenceOf(const BoundTrackState& s, const TrackingSurface& p)
{
  return (p.normal[0] * std::cos(s.phi()) * std::sin(s.theta())) + (p.normal[1] * std::sin(s.phi()) * std::sin(s.theta())) + (p.normal[2] * std::cos(s.theta()));
}

/// The state propagated onto \a to: forward, else backward (the points are on
/// both sides at a perturbation's distance).
bool onto(BoundTrackState& s, const TrackingSurface& from, const TrackingSurface& to)
{
  Propagator prop;
  prop.bzKiloGauss = static_cast<Scalar>(kBz);
  auto copy = s;
  if (propagate(prop, copy, from, to, PropagationDirection::kForward).status == PropagationStatus::kOk) {
    s = copy;
    return true;
  }
  copy = s;
  if (propagate(prop, copy, from, to, PropagationDirection::kBackward).status == PropagationStatus::kOk) {
    s = copy;
    return true;
  }
  return false;
}

/// The O2 parameters of a state already on O2's plane \a target.
std::array<double, 5> barrelParams(const BoundTrackState& s, const TrackingSurface& target, double alpha)
{
  const auto t = toTrackParCov(s, target, alpha, kBz);
  return {t.getY(), t.getZ(), t.getSnp(), t.getTgl(), t.getQ2Pt()};
}

std::array<double, 5> forwardParams(const BoundTrackState& s, const TrackingSurface& target)
{
  const auto t = toTrackParCovFwd(s, target, kBz, 0.);
  return {t.getX(), t.getY(), t.getPhi(), t.getTanl(), t.getInvQPt()};
}

/// The Jacobian's columns by central differences: each parameter moved by
/// +-h, propagated from \a from onto \a target, written as O2's.
template <typename Params>
detail::Jacobian finiteDifferences(const BoundTrackState& s, const TrackingSurface& from, const TrackingSurface& target, Params params, bool& ok)
{
  const std::array<double, kNBoundParams> steps{2e-3, 2e-3, 2e-4, 2e-4, 1e-3 * std::abs(static_cast<double>(s.qOverP()))};
  detail::Jacobian jac{};
  ok = true;
  for (int k = 0; k < kNBoundParams; ++k) {
    auto plus = s;
    auto minus = s;
    plus.par[k] += static_cast<Scalar>(steps[k]);
    minus.par[k] -= static_cast<Scalar>(steps[k]);
    const double h = static_cast<double>(plus.par[k]) - static_cast<double>(minus.par[k]); // as stored
    ok = ok && onto(plus, from, target) && onto(minus, from, target);
    if (!ok) {
      return jac;
    }
    const auto p = params(plus);
    const auto m = params(minus);
    for (int i = 0; i < kNBoundParams; ++i) {
      double d = p[i] - m[i];
      if (i == 2 && d > std::numbers::pi) { // a forward phi across its seam
        d -= 2. * std::numbers::pi;
      } else if (i == 2 && d < -std::numbers::pi) {
        d += 2. * std::numbers::pi;
      }
      jac[i][k] = d / h;
    }
  }
  return jac;
}

void requireClose(const detail::Jacobian& analytic, const detail::Jacobian& numeric)
{
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int k = 0; k < kNBoundParams; ++k) {
      BOOST_TEST_CONTEXT("d" << i << "/d" << k << ": analytic " << analytic[i][k] << ", numeric " << numeric[i][k])
      {
        BOOST_TEST(std::abs(analytic[i][k] - numeric[i][k]) <= 2e-4 + (1e-2 * std::abs(analytic[i][k])));
      }
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(BarrelParametersAreTheStatesPointMomentumAndCharge)
{
  std::mt19937 random{4711U};
  std::uniform_real_distribution<double> uniform{0., 1.};
  unsigned inward = 0U;
  unsigned tested = 0U;
  while (tested < 2000U) {
    const double sensorAlpha = (2. * uniform(random) - 1.) * std::numbers::pi;
    const auto p = barrelSensor(sensorAlpha, 2. + (40. * uniform(random)), 30. * (uniform(random) - 0.5), 0.2 * (uniform(random) - 0.5), 0.05 * (uniform(random) - 0.5));
    const auto s = state(1.5 * (uniform(random) - 0.5), 10. * (uniform(random) - 0.5), (2. * uniform(random) - 1.) * std::numbers::pi, 0.1 + (2.94 * uniform(random)),
                         (uniform(random) < 0.5 ? -1. : 1.) * (0.1 + (10. * uniform(random))));
    if (std::abs(cosIncidenceOf(s, p)) < 0.05) {
      continue;
    }
    ++tested;
    const double alpha = frameAlpha(s, sensorAlpha);
    BOOST_TEST(alpha >= -std::numbers::pi);
    BOOST_TEST(alpha < std::numbers::pi);
    inward += std::abs(std::remainder(alpha - sensorAlpha, 2. * std::numbers::pi)) > 1. ? 1U : 0U;
    const auto t = toTrackParCov(s, p, alpha, kBz);
    double r[3];
    double mom[3];
    globalOf(s, p, r, mom);
    const auto xyz = t.getXYZGlo();
    std::array<float, 3> pxyz{};
    BOOST_TEST(t.getPxPyPzGlo(pxyz));
    BOOST_TEST(std::abs(xyz.X() - r[0]) < 1e-4);
    BOOST_TEST(std::abs(xyz.Y() - r[1]) < 1e-4);
    BOOST_TEST(std::abs(xyz.Z() - r[2]) < 1e-4);
    for (int i = 0; i < 3; ++i) {
      BOOST_TEST(std::abs(pxyz[i] - mom[i]) <= 1e-5 * (1. + std::hypot(mom[0], mom[1], mom[2])));
    }
    BOOST_TEST(t.getSign() == (s.qOverP() > 0.F ? 1 : -1));
  }
  // Both branches were taken: the frame turns for states moving inward.
  BOOST_TEST(inward > 100U);
}

BOOST_AUTO_TEST_CASE(ForwardParametersAreTheStatesPointMomentumAndCharge)
{
  std::mt19937 random{4712U};
  std::uniform_real_distribution<double> uniform{0., 1.};
  for (int n = 0; n < 2000; ++n) {
    const auto p = discSensor(20. * (uniform(random) - 0.5), 20. * (uniform(random) - 0.5), -45. - (30. * uniform(random)), uniform(random) < 0.5 ? -1. : 1.,
                              0.02 * (uniform(random) - 0.5));
    // Towards the muon arm, sin(theta) down to 0.05; some moving back.
    const double theta = uniform(random) < 0.9 ? std::numbers::pi - 0.05 - (0.6 * uniform(random)) : 0.05 + (0.6 * uniform(random));
    const auto s = state(1.5 * (uniform(random) - 0.5), 0.7 * (uniform(random) - 0.5), (2. * uniform(random) - 1.) * std::numbers::pi, theta,
                         (uniform(random) < 0.5 ? -1. : 1.) * (0.1 + (10. * uniform(random))));
    const auto t = toTrackParCovFwd(s, p, kBz, 3.5);
    double r[3];
    double mom[3];
    globalOf(s, p, r, mom);
    BOOST_TEST(std::abs(t.getX() - r[0]) < 1e-9);
    BOOST_TEST(std::abs(t.getY() - r[1]) < 1e-9);
    BOOST_TEST(std::abs(t.getZ() - r[2]) < 1e-9);
    BOOST_TEST(std::abs(t.getPx() - mom[0]) <= 1e-7 * std::hypot(mom[0], mom[1], mom[2]));
    BOOST_TEST(std::abs(t.getPy() - mom[1]) <= 1e-7 * std::hypot(mom[0], mom[1], mom[2]));
    BOOST_TEST(std::abs(t.getPz() - mom[2]) <= 1e-7 * std::hypot(mom[0], mom[1], mom[2]));
    BOOST_TEST(t.getCharge() == (s.qOverP() > 0.F ? 1 : -1));
    BOOST_TEST(t.getTrackChi2() == 3.5);
  }
}

BOOST_AUTO_TEST_CASE(BarrelJacobianIsThePropagatorsOntoTheFrame)
{
  std::mt19937 random{4713U};
  std::uniform_real_distribution<double> uniform{0., 1.};
  unsigned tested = 0U;
  for (int n = 0; n < 4000 && tested < 400U; ++n) {
    const double sensorAlpha = (2. * uniform(random) - 1.) * std::numbers::pi;
    const auto p = barrelSensor(sensorAlpha, 3. + (40. * uniform(random)), 20. * (uniform(random) - 0.5), 1.2 * (uniform(random) - 0.5), 0.6 * (uniform(random) - 0.5));
    const auto s = state(1.5 * (uniform(random) - 0.5), 10. * (uniform(random) - 0.5), (2. * uniform(random) - 1.) * std::numbers::pi, 0.3 + (2.54 * uniform(random)),
                         (uniform(random) < 0.5 ? -1. : 1.) * (0.1 + (8. * uniform(random))));
    if (std::abs(cosIncidenceOf(s, p)) < 0.2) {
      continue;
    }
    const double alpha = frameAlpha(s, sensorAlpha);
    if (std::cos(static_cast<double>(s.phi()) - alpha) < 0.2) {
      continue;
    }
    double r[3];
    double mom[3];
    globalOf(s, p, r, mom);
    const double e[3] = {std::cos(alpha), std::sin(alpha), 0.};
    const double up[3] = {0., 0., 1.};
    const auto target = plane({r[0], r[1], r[2]}, e, up);
    bool ok = false;
    const auto numeric = finiteDifferences(s, p, target, [&](const BoundTrackState& at) { return barrelParams(at, target, alpha); }, ok);
    if (!ok) {
      continue;
    }
    ++tested;
    BOOST_TEST_CONTEXT("sensor alpha " << sensorAlpha << ", frame alpha " << alpha << ", phi " << s.phi() << ", theta " << s.theta() << ", q/p " << s.qOverP())
    {
      requireClose(detail::barrelJacobian(s, detail::shiftOnto(s, p, e, kBz), alpha), numeric);
    }
  }
  BOOST_TEST(tested >= 400U);
}

BOOST_AUTO_TEST_CASE(ForwardJacobianIsThePropagatorsOntoFixedZ)
{
  std::mt19937 random{4714U};
  std::uniform_real_distribution<double> uniform{0., 1.};
  unsigned tested = 0U;
  for (int n = 0; n < 4000 && tested < 400U; ++n) {
    const auto p = discSensor(20. * (uniform(random) - 0.5), 20. * (uniform(random) - 0.5), -45. - (30. * uniform(random)), uniform(random) < 0.5 ? -1. : 1.,
                              0.6 * (uniform(random) - 0.5));
    const double theta = std::numbers::pi - 0.1 - (0.6 * uniform(random));
    const auto s = state(1.5 * (uniform(random) - 0.5), 0.7 * (uniform(random) - 0.5), (2. * uniform(random) - 1.) * std::numbers::pi, theta,
                         (uniform(random) < 0.5 ? -1. : 1.) * (0.1 + (8. * uniform(random))));
    double r[3];
    double mom[3];
    globalOf(s, p, r, mom);
    const double e[3] = {0., 0., 1.};
    const double up[3] = {0., 1., 0.};
    const auto target = plane({r[0], r[1], r[2]}, e, up);
    bool ok = false;
    const auto numeric = finiteDifferences(s, p, target, [&](const BoundTrackState& at) { return forwardParams(at, target); }, ok);
    if (!ok) {
      continue;
    }
    ++tested;
    BOOST_TEST_CONTEXT("phi " << s.phi() << ", theta " << s.theta() << ", q/p " << s.qOverP())
    {
      requireClose(detail::forwardJacobian(s, detail::shiftOnto(s, p, e, kBz)), numeric);
    }
  }
  BOOST_TEST(tested >= 400U);
}

} // namespace o2::gtrk::workflow
