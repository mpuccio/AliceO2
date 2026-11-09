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
#ifndef GTRK_TEST_PROPAGATIONPRECISIONFIXTURE_H_
#define GTRK_TEST_PROPAGATIONPRECISIONFIXTURE_H_
#include "GenericTracking/Tracking/Propagator.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
namespace o2::gtrk::propagation_precision
{

constexpr double kMicron = 1.e-4; ///< [cm]

inline TrackingSurface plane(double x, double y, double z, double angle, bool disc, SensorID id)
{
  TrackingSurface result;
  result.id = id;
  result.origin[0] = static_cast<float>(x);
  result.origin[1] = static_cast<float>(y);
  result.origin[2] = static_cast<float>(z);
  if (disc) {
    result.normal[2] = 1.F;
    result.axis0[0] = static_cast<float>(std::cos(angle));
    result.axis0[1] = static_cast<float>(std::sin(angle));
    result.axis1[0] = -result.axis0[1];
    result.axis1[1] = result.axis0[0];
  } else {
    result.normal[0] = static_cast<float>(std::cos(angle));
    result.normal[1] = static_cast<float>(std::sin(angle));
    result.axis0[0] = -result.normal[1];
    result.axis0[1] = result.normal[0];
    result.axis1[2] = 1.F;
  }
  result.bounds = {.halfX = 1.e5F, .halfY = 1.e5F};
  return result;
}

/// How far the stored answer is from the helix it claims to be on.
struct Closure {
  double position = 0.; ///< stored local position against the helix point [cm]
  double plane = 0.;    ///< that point's distance from the target plane [cm]
  double azimuth = 0.;  ///< stored azimuth against the accumulated turn [rad]
};

/// Re-evaluate the helix in double at the path length the propagation
/// reported; its float rounding is 1e-7 relative, 30 nm over 3 m.
inline Closure closure(const Propagator& prop, const BoundTrackState& before, const BoundTrackState& after, const TrackingSurface& from, const TrackingSurface& to, Scalar reportedPath)
{
  const auto path = static_cast<double>(reportedPath);
  const auto phi = static_cast<double>(before.phi());
  const auto theta = static_cast<double>(before.theta());
  // The turn per unit path carries no theta: q/pT's sin(theta) cancels against
  // the transverse part of the step.
  const double turn = static_cast<double>(before.qOverP()) * static_cast<double>(prop.bzKiloGauss) * static_cast<double>(constants::math::B2C) * path;
  const double transverse = path * std::sin(theta);
  const double along = transverse * sinc(turn);
  const double across = transverse * versine(turn);
  const double step[3] = {along * std::cos(phi) - across * std::sin(phi), along * std::sin(phi) + across * std::cos(phi), path * std::cos(theta)};

  double local[2] = {0., 0.};
  double residual = 0.;
  for (int i = 0; i < 3; ++i) {
    const double delta = (static_cast<double>(from.origin[i]) - static_cast<double>(to.origin[i])) + static_cast<double>(before.loc0()) * static_cast<double>(from.axis0[i]) + static_cast<double>(before.loc1()) * static_cast<double>(from.axis1[i]) + step[i];
    local[0] += delta * static_cast<double>(to.axis0[i]);
    local[1] += delta * static_cast<double>(to.axis1[i]);
    residual += delta * static_cast<double>(to.normal[i]);
  }

  Closure out;
  out.position = std::hypot(local[0] - static_cast<double>(after.loc0()), local[1] - static_cast<double>(after.loc1()));
  out.plane = std::abs(residual);
  out.azimuth = std::abs(std::remainder(phi + turn - static_cast<double>(after.phi()), 2. * std::numbers::pi));
  return out;
}

struct Worst {
  Closure value;
  int checked = 0;

  void add(const Closure& next)
  {
    value.position = std::max(value.position, next.position);
    value.plane = std::max(value.plane, next.plane);
    value.azimuth = std::max(value.azimuth, next.azimuth);
    ++checked;
  }

  void report(const char* label) const { std::cout << label << ": checked=" << checked << " position um=" << value.position / kMicron << " plane um=" << value.plane / kMicron << " azimuth urad=" << value.azimuth * 1.e6 << '\n'; }
};

} // namespace o2::gtrk::propagation_precision
#endif
