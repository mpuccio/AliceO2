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

#ifndef GTRK_TOF_TRANSPORT_H_
#define GTRK_TOF_TRANSPORT_H_

#include "GenericTracking/TOF/Matcher.h"
#include "GenericTracking/Tracking/TrackFitter.h"
#include "CommonConstants/PhysicsConstants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace o2::gtrk::tof
{
inline constexpr std::array<double, kHypotheses> kMassGeV{
  o2::constants::physics::MassElectron, o2::constants::physics::MassMuon,
  o2::constants::physics::MassPionCharged, o2::constants::physics::MassKaonCharged,
  o2::constants::physics::MassProton};
inline constexpr double kLightSpeedCMNS = 29.9792458;

struct TransportResult {
  BoundTrackState state;
  double flightNS = 0.; ///< signed along physical momentum
  double pathCM = 0.;
  bool ok = false;
  bool inside = false;
};

/// Bounded helix steps, with the existing LUT energy loss/scattering model.
/// Intermediate planes are normal to the helix; the target is the real sensor.
/// This host implementation is a reference for future indexed/device matching.
template <typename Material>
TransportResult transport(const BoundTrackState& initial, TrackingSurface from, const TrackingSurface& target,
                          Propagator prop, const Material& material, PropagationDirection direction,
                          const Parameters& parameters, bool withCovariance = false)
{
  TransportResult output{.state = initial};
  if (!(parameters.maxStepCM > 0.) || !(parameters.maxPathCM > 0.)) {
    throw std::invalid_argument("TOF propagation requires positive step and path limits");
  }
  CovarianceFactors factors;
  if (withCovariance) {
    detail::factorUD(initial.cov, factors.upper, factors.diagonal);
    if (!positiveDiagonal(factors)) {
      return output;
    }
  }
  const auto inverseBeta = [&](const BoundTrackState& state) {
    const double mq = static_cast<double>(prop.massGeV) * state.qOverP();
    return std::sqrt(1. + mq * mq);
  };
  const std::size_t limit = static_cast<std::size_t>(std::ceil(parameters.maxPathCM / parameters.maxStepCM)) + 2;
  for (std::size_t stepIndex = 0; stepIndex < limit; ++stepIndex) {
    auto probe = output.state;
    const auto remaining = propagate(prop, probe, from, target, direction);
    if (remaining.status != PropagationStatus::kOk && remaining.status != PropagationStatus::kOutOfBounds) {
      return output;
    }
    const double path = remaining.pathLength;
    const bool last = std::abs(path) <= parameters.maxStepCM;
    if (std::abs(output.pathCM) + std::abs(path) > parameters.maxPathCM) {
      return output;
    }
    auto to = target;
    if (!last) {
      // A plane through the exact helix point at a bounded path length.
      const double length = std::copysign(parameters.maxStepCM, path);
      const double phi = output.state.phi();
      const double theta = output.state.theta();
      const double turn = output.state.qOverP() * prop.bzKiloGauss * o2::constants::math::B2C * length;
      const double sinc = std::abs(turn) < 1e-5 ? 1. - turn * turn / 6. : std::sin(turn) / turn;
      const double versine = std::abs(turn) < 1e-5 ? turn / 2. - turn * turn * turn / 24. : (1. - std::cos(turn)) / turn;
      Scalar position[3];
      surfaceToGlobal(from, output.state.loc0(), output.state.loc1(), position);
      to = {};
      to.origin[0] = static_cast<float>(position[0] + length * std::sin(theta) * (std::cos(phi) * sinc - std::sin(phi) * versine));
      to.origin[1] = static_cast<float>(position[1] + length * std::sin(theta) * (std::sin(phi) * sinc + std::cos(phi) * versine));
      to.origin[2] = static_cast<float>(position[2] + length * std::cos(theta));
      const double arrival = phi + turn;
      to.axis0[0] = static_cast<float>(-std::sin(arrival));
      to.axis0[1] = static_cast<float>(std::cos(arrival));
      to.axis1[0] = static_cast<float>(-std::cos(theta) * std::cos(arrival));
      to.axis1[1] = static_cast<float>(-std::cos(theta) * std::sin(arrival));
      to.axis1[2] = static_cast<float>(std::sin(theta));
      to.normal[0] = static_cast<float>(std::sin(theta) * std::cos(arrival));
      to.normal[1] = static_cast<float>(std::sin(theta) * std::sin(arrival));
      to.normal[2] = static_cast<float>(std::cos(theta));
    }
    auto next = output.state;
    SolenoidJacobian jacobian;
    PropagationResult propagated;
    if (jacobianAnalytic(prop, output.state, from, to, direction, jacobian, &next, &propagated) != PropagationStatus::kOk) {
      return output;
    }
    Scalar start[3], end[3];
    surfaceToGlobal(from, output.state.loc0(), output.state.loc1(), start);
    surfaceToGlobal(to, next.loc0(), next.loc1(), end);
    const auto budget = detail::stepMaterial(material, start, end);
    const double beforeBeta = inverseBeta(output.state);
    if (withCovariance) {
      const auto betaP = betaMomentum(1.F / std::abs(next.qOverP()), prop.massGeV);
      const auto before = scatteringNoise(highlandTheta0(budget.nearerX2X0, betaP), output.state.theta());
      if (detail::transportStep<false>(prop, budget, betaP, before, direction, jacobian, next, factors, 1.F) != PropagationStatus::kOk) {
        return output;
      }
      detail::covarianceFromUD(factors.upper, factors.diagonal, next.cov);
    } else {
      Scalar variance;
      if (!detail::applyEnergyLoss(prop, budget.xRho, direction, next, variance)) {
        return output;
      }
    }
    // Trapezoid integration of 1/beta over a short material step.
    output.flightNS += propagated.pathLength * (beforeBeta + inverseBeta(next)) / (2. * kLightSpeedCMNS);
    output.pathCM += propagated.pathLength;
    output.state = next;
    from = to;
    if (last) {
      output.ok = std::isfinite(output.flightNS);
      output.inside = target.bounds.isInside(next.loc0(), next.loc1());
      return output;
    }
  }
  return output;
}

template <typename Material>
Prediction predict(const Track& track, const Hit& hit, Propagator prop, const Material& material, const Parameters& parameters = {})
{
  Prediction prediction;
  // Inner TOF can precede the innermost fitted point of a displaced/short
  // track. Choose the nearer crossing, retaining the signed flight time.
  auto direction = PropagationDirection::kForward;
  if (hit.layer == Layer::Inner) {
    auto forward = track.state, backward = track.state;
    const auto a = propagate(prop, forward, track.surface, hit.surface, PropagationDirection::kForward);
    const auto b = propagate(prop, backward, track.surface, hit.surface, PropagationDirection::kBackward);
    const bool aOK = a.status == PropagationStatus::kOk || a.status == PropagationStatus::kOutOfBounds;
    const bool bOK = b.status == PropagationStatus::kOk || b.status == PropagationStatus::kOutOfBounds;
    if (bOK && (!aOK || std::abs(b.pathLength) < std::abs(a.pathLength))) {
      direction = PropagationDirection::kBackward;
    }
  }
  const bool outer = hit.layer == Layer::Outer;
  const auto spatial = transport(outer ? track.outwardState : track.state,
                                 outer ? track.outwardSurface : track.surface, hit.surface,
                                 prop, material, outer ? PropagationDirection::kForward : direction, parameters, true);
  if (!spatial.ok || !spatial.inside) {
    return prediction;
  }
  const double a = spatial.state.variance(kLoc0) + hit.varLoc0;
  const double b = spatial.state.covAt(kLoc0, kLoc1);
  const double d = spatial.state.variance(kLoc1) + hit.varLoc1;
  const double det = a * d - b * b;
  const double x = hit.loc0 - spatial.state.loc0();
  const double y = hit.loc1 - spatial.state.loc1();
  if (!(a > 0. && d > 0. && det > 0.)) {
    return prediction;
  }
  prediction.spatialChi2 = static_cast<float>((d * x * x - 2. * b * x * y + a * y * y) / det);
  if (!(prediction.spatialChi2 <= parameters.maxSpatialChi2)) {
    return prediction;
  }
  for (std::size_t mass = 0; mass < kHypotheses; ++mass) {
    prop.massGeV = static_cast<float>(kMassGeV[mass]);
    const auto central = transport(track.state, track.surface, hit.surface, prop, material, direction, parameters);
    if (!central.ok || !central.inside) {
      continue;
    }
    prediction.flightNS[mass] = central.flightNS;
    bool valid = true;
    for (int parameter = 0; parameter < kNBoundParams; ++parameter) {
      auto high = track.state, low = track.state;
      const double epsilon = parameter < 2 ? 0.001 : parameter < 4 ? 0.0001
                                                                   : std::max(1e-5, std::abs(static_cast<double>(track.state.qOverP())) * 0.001);
      high.par[parameter] += static_cast<float>(epsilon);
      low.par[parameter] -= static_cast<float>(epsilon);
      const auto hi = transport(high, track.surface, hit.surface, prop, material, direction, parameters);
      const auto lo = transport(low, track.surface, hit.surface, prop, material, direction, parameters);
      if (!hi.ok || !lo.ok) {
        valid = false;
        break;
      }
      prediction.derivative[mass][parameter] = (hi.flightNS - lo.flightNS) / (high.par[parameter] - low.par[parameter]);
    }
    prediction.valid[mass] = valid;
  }
  return prediction;
}
} // namespace o2::gtrk::tof
#endif
