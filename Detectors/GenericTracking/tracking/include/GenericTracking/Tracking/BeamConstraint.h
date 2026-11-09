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

/// \file BeamConstraint.h
/// \brief The beam line as a measurement of a state inside the tracker.
///
/// A road on the outer layers measures the curvature over a short lever arm;
/// the beam line, for a primary, adds the longest one. The beam point is
/// taken as a measurement of the state on the road's innermost surface, of
/// one coordinate: the transverse distance of closest approach, d0. Its
/// prediction is the state carried to the plane through the beam line normal
/// to the track's direction there (loc0 along the transverse normal: d0;
/// loc1 along the beam: z0), and its variance is the beam's width squared
/// plus the scattering between: a kink theta_k at distance s_k from the beam
/// moves the d0 the state predicts by theta_k s_k, so the sum
/// `sum_k theta_k^2 s_k^2` over the layers and the passive scatterers inside,
/// and the surface's own layer, as the pair windows take it
/// (`sourceScattering`), at the track's own momentum. The update is the
/// Kalman update of that one coordinate, made on the surface: the scattering
/// between is in the measurement's variance, not in the state's.
///
/// z0 is left out: the beam spot's length, centimetres, constrains nothing
/// the road does not.

#ifndef GTRK_TRACKING_BEAMCONSTRAINT_H_
#define GTRK_TRACKING_BEAMCONSTRAINT_H_

#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/CovarianceTransport.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/Trig.h"

#include "GPUCommonMath.h"
#include "MathUtils/detail/trigonometric.h"

namespace o2::gtrk
{

/// What the beam-line update gave: whether it could be made, and its chi2,
/// one degree of freedom.
struct BeamLineUpdate {
  bool ok = false;
  Scalar chi2 = 0.F;
  Scalar residual = 0.F;   ///< the beam's d0 less the prediction's [cm]
  Scalar variance = 0.F;   ///< the residual's [cm^2]
  Scalar scattering = 0.F; ///< of which the scattering between [cm^2]
};

/// The direction [rad] of a bound state's helix at its transverse closest
/// approach to the z axis, without a virtual beam-line surface.
///
/// With A = r.dot(direction) and B = r.dot(leftNormal), closest approach obeys
///
///   kappa A cos(psi) + (1 + kappa B) sin(psi) = 0,
///
/// solved for the turn psi; at zero curvature it is 0.
template <typename Prop>
[[nodiscard]] GPUdi() Scalar perigeePhi(const Prop& prop, const BoundTrackState& state, const TrackingSurface& surface)
{
  ExactScalar position[3]{};
  detail::boundToGlobal(surface, static_cast<ExactScalar>(state.loc0()), static_cast<ExactScalar>(state.loc1()), position);

  const auto phi = static_cast<ExactScalar>(state.phi());
  const auto theta = static_cast<ExactScalar>(state.theta());
  ExactScalar sinPhi = 0.;
  ExactScalar cosPhi = 0.;
  ExactScalar sinTheta = 0.;
  ExactScalar cosTheta = 0.;
  sincosOf(phi, sinPhi, cosPhi);
  sincosOf(theta, sinTheta, cosTheta);

  const ExactScalar curvature = pointPropagator(prop, surface, state).curvature(static_cast<ExactScalar>(state.qOverP()), sinTheta);
  const ExactScalar along = position[0] * cosPhi + position[1] * sinPhi;
  const ExactScalar across = -position[0] * sinPhi + position[1] * cosPhi;
  const ExactScalar numerator = -curvature * along;
  const ExactScalar denominator = 1. + curvature * across;
  const ExactScalar turn = atan2Of(numerator, denominator);
  ExactScalar phi0 = phi + turn;
  math_utils::detail::bringToPMPi<ExactScalar>(phi0);
  return static_cast<Scalar>(phi0);
}

/// The plane through the beam line at (\a x, \a y) normal to the transverse
/// direction \a phi: loc0 the signed transverse distance along
/// (-sin phi, cos phi), at closest approach the track's d0, loc1 z. It is no
/// sensor: its
/// bounds are empty, and the update takes the propagation's out-of-bounds
/// status as a crossing.
GPUhdi() GPUAlwaysInline() TrackingSurface beamLinePlane(Scalar x, Scalar y, Scalar phi)
{
  Scalar sinPhi = 0.F;
  Scalar cosPhi = 0.F;
  sincosOf(phi, sinPhi, cosPhi);
  TrackingSurface plane;
  plane.origin[0] = x;
  plane.origin[1] = y;
  plane.axis0[0] = -sinPhi;
  plane.axis0[1] = cosPhi;
  plane.axis1[2] = 1.F;
  plane.normal[0] = cosPhi;
  plane.normal[1] = sinPhi;
  return plane;
}

/// \a state, on \a surface of layer \a layer, updated with the iteration's
/// beam line; \a innermost is the state's measurement as the pair cut reads
/// it (the corrected beam frame), whose line from the track's z0 sets the
/// scattering's path. The state is left as it was if the update cannot be
/// made: the propagation to the plane failed, or the residual's variance is
/// not positive. \a scatteringScale multiplies the scattering term: 1 but
/// for the kink allowance.
template <DetectorLayout Det>
GPUdi() BeamLineUpdate updateWithBeamLine(const Propagator& prop, const IterationParameters<Det>& parameters, const TrackingSurface& surface, LayerID layer, const PairPoint& innermost, BoundTrackState& state, Scalar scatteringScale = 1.F)
{
  using Math = o2::gpu::GPUCommonMath;
  BeamLineUpdate result;
  const auto plane = beamLinePlane(parameters.beamX, parameters.beamY, perigeePhi(prop, state, surface));
  SolenoidJacobian jacobian;
  BoundTrackState predicted;
  const auto status = jacobianAnalytic(prop, state, surface, plane, PropagationDirection::kBackward, jacobian, &predicted);
  if (status != PropagationStatus::kOk && status != PropagationStatus::kOutOfBounds) {
    return result;
  }

  Scalar sinTheta = 0.F;
  Scalar cosTheta = 0.F;
  sincosOf(state.theta(), sinTheta, cosTheta);
  const Scalar betaPT = betaMomentum(1.F / Math::Abs(state.qOverP()), prop.massGeV) * Math::Abs(sinTheta);
  PairCut ownMomentum;
  ownMomentum.invBetaMomentum2 = 1.F / (betaPT * betaPT);
  const VertexSeed origin{.z = predicted.loc1()};
  const Scalar scattering = scatteringScale * sourceScattering<Det>(parameters.scattering, ownMomentum, origin, layer, innermost);

  // H is the Jacobian's loc0 row; u = C H^T, s = H C H^T + V. In double: the
  // update subtracts u u^T / s from a covariance whose entries span orders of
  // magnitude.
  ExactScalar u[kNBoundParams]{};
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j < kNBoundParams; ++j) {
      u[i] += static_cast<ExactScalar>(state.covAt(i, j)) * static_cast<ExactScalar>(jacobian.loc0[j]);
    }
  }
  ExactScalar variance = static_cast<ExactScalar>(parameters.beamVarTransverse) + static_cast<ExactScalar>(scattering);
  for (int j = 0; j < kNBoundParams; ++j) {
    variance += static_cast<ExactScalar>(jacobian.loc0[j]) * u[j];
  }
  if (!(variance > 0.)) {
    return result;
  }
  const auto residual = -static_cast<ExactScalar>(predicted.loc0());
  for (int i = 0; i < kNBoundParams; ++i) {
    state.par[i] = static_cast<Scalar>(static_cast<ExactScalar>(state.par[i]) + (u[i] * residual / variance));
    for (int j = 0; j <= i; ++j) {
      state.covAt(i, j) = static_cast<Scalar>(static_cast<ExactScalar>(state.covAt(i, j)) - (u[i] * u[j] / variance));
    }
  }
  result.ok = true;
  result.chi2 = static_cast<Scalar>(residual * residual / variance);
  result.residual = static_cast<Scalar>(residual);
  result.variance = static_cast<Scalar>(variance);
  result.scattering = scattering;
  return result;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_BEAMCONSTRAINT_H_
