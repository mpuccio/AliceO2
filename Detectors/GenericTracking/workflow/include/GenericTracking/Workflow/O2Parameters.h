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

/// \file O2Parameters.h
/// \brief A bound track state as O2's barrel (TrackParCov) or forward
///        (TrackParCovFwd) parameters, at the state's own point.
///
/// No propagation: the O2 parameters are taken on the plane through the
/// state's point with O2's orientation -- normal (cos a, sin a, 0) for the
/// barrel, z for the forward ones. Only the covariance sees the change of
/// plane: a shift in the source plane moves the point along the track to the
/// target plane by s = -(e . A) / (e . t), over which the direction turns by
/// (q/p) Bz B2C s. The terms vanish where the two planes are parallel.

#ifndef GTRK_WORKFLOW_O2PARAMETERS_H_
#define GTRK_WORKFLOW_O2PARAMETERS_H_

#include "CommonConstants/MathConstants.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include <array>
#include <cmath>
#include <numbers>

namespace o2::gtrk::workflow
{

namespace detail
{

using Jacobian = std::array<std::array<double, kNBoundParams>, kNBoundParams>;

/// J C J^T of the state's covariance, as O2's packed lower triangle (row by
/// row, (i, j) at i (i + 1) / 2 + j): TrackParCov's order, and the order
/// SMatrix's symmetric representation takes.
inline std::array<double, kBoundCovSize> transported(const Jacobian& jac, const BoundTrackState& state)
{
  std::array<double, kBoundCovSize> result{};
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = 0.;
      for (int k = 0; k < kNBoundParams; ++k) {
        for (int l = 0; l < kNBoundParams; ++l) {
          sum += jac[i][k] * static_cast<double>(state.covAt(k, l)) * jac[j][l];
        }
      }
      result[(i * (i + 1) / 2) + j] = sum;
    }
  }
  return result;
}

/// What both maps share: the point, the direction, and the derivatives of
/// the point and the azimuth along a shift of loc0 and loc1 that keeps the
/// point on the plane of normal \a e.
struct Shift {
  double r[3];
  double t[3];
  double dr[2][3];
  double dPhi[2];
};

inline Shift shiftOnto(const BoundTrackState& state, const TrackingSurface& plane, const double e[3], double bzKiloGauss)
{
  Shift s{};
  const double phi = state.phi();
  const double theta = state.theta();
  s.t[0] = std::cos(phi) * std::sin(theta);
  s.t[1] = std::sin(phi) * std::sin(theta);
  s.t[2] = std::cos(theta);
  const double et = (e[0] * s.t[0]) + (e[1] * s.t[1]) + (e[2] * s.t[2]);
  // d phi / d s of the helix, GenericTracking's (CovarianceTransport.h).
  const double turn = static_cast<double>(state.qOverP()) * bzKiloGauss * constants::math::B2C;
  for (int i = 0; i < 3; ++i) {
    s.r[i] = static_cast<double>(plane.origin[i]) + (static_cast<double>(state.loc0()) * plane.axis0[i]) + (static_cast<double>(state.loc1()) * plane.axis1[i]);
  }
  for (int k = 0; k < 2; ++k) {
    const Scalar* axis = k == 0 ? plane.axis0 : plane.axis1;
    const double path = -((e[0] * axis[0]) + (e[1] * axis[1]) + (e[2] * axis[2])) / et;
    for (int i = 0; i < 3; ++i) {
      s.dr[k][i] = axis[i] + (s.t[i] * path);
    }
    s.dPhi[k] = turn * path;
  }
  return s;
}

/// d (y, z, snp, tgl, q/pt) / d (loc0, loc1, phi, theta, q/p) in the frame of
/// angle \a alpha, \a s the shift onto its plane.
inline Jacobian barrelJacobian(const BoundTrackState& state, const Shift& s, double alpha)
{
  const double ca = std::cos(alpha);
  const double sa = std::sin(alpha);
  const double theta = state.theta();
  const double sinTheta = std::sin(theta);
  const double cosLocal = std::cos(static_cast<double>(state.phi()) - alpha);
  Jacobian jac{};
  for (int k = 0; k < 2; ++k) {
    jac[0][k] = (-sa * s.dr[k][0]) + (ca * s.dr[k][1]);
    jac[1][k] = s.dr[k][2];
    jac[2][k] = cosLocal * s.dPhi[k];
  }
  jac[2][kPhi] = cosLocal;
  jac[3][kTheta] = -1. / (sinTheta * sinTheta);
  jac[4][kTheta] = -static_cast<double>(state.qOverP()) * std::cos(theta) / (sinTheta * sinTheta);
  jac[4][kQOverP] = 1. / sinTheta;
  return jac;
}

/// d (x, y, phi, tanl, q/pt) / d (loc0, loc1, phi, theta, q/p) at fixed z,
/// \a s the shift onto that plane.
inline Jacobian forwardJacobian(const BoundTrackState& state, const Shift& s)
{
  const double theta = state.theta();
  const double sinTheta = std::sin(theta);
  Jacobian jac{};
  for (int k = 0; k < 2; ++k) {
    jac[0][k] = s.dr[k][0];
    jac[1][k] = s.dr[k][1];
    jac[2][k] = s.dPhi[k];
  }
  jac[2][kPhi] = 1.;
  jac[3][kTheta] = -1. / (sinTheta * sinTheta);
  jac[4][kTheta] = -static_cast<double>(state.qOverP()) * std::cos(theta) / (sinTheta * sinTheta);
  jac[4][kQOverP] = 1. / sinTheta;
  return jac;
}

} // namespace detail

/// The barrel frame a state is written in: its sensor's O2 tracking-frame
/// angle \a sensorAlpha, turned by pi when the state's momentum points against
/// that frame's x (a state moving inward through its sensor), since O2's
/// parameters are the branch with momentum along +x. In [-pi, pi).
[[nodiscard]] inline double frameAlpha(const BoundTrackState& state, double sensorAlpha)
{
  const double alpha = std::cos(static_cast<double>(state.phi()) - sensorAlpha) < 0. ? sensorAlpha + std::numbers::pi : sensorAlpha;
  return alpha - (2. * std::numbers::pi * std::floor((alpha + std::numbers::pi) / (2. * std::numbers::pi)));
}

/// \a state, bound to \a plane, as O2's barrel parameters in the frame of
/// angle \a alpha (`frameAlpha`): x and y the point's coordinates in that
/// frame, z, snp = sin(phi - alpha), tgl = cot(theta), q/pt = (q/p) /
/// sin(theta). The field \a bzKiloGauss enters the covariance only.
[[nodiscard]] inline o2::track::TrackParCov toTrackParCov(const BoundTrackState& state, const TrackingSurface& plane, double alpha, double bzKiloGauss)
{
  const double ca = std::cos(alpha);
  const double sa = std::sin(alpha);
  const double e[3] = {ca, sa, 0.};
  const auto s = detail::shiftOnto(state, plane, e, bzKiloGauss);
  const double phi = state.phi();
  const double theta = state.theta();
  const double sinTheta = std::sin(theta);
  const double qOverP = state.qOverP();
  const auto cov = detail::transported(detail::barrelJacobian(state, s, alpha), state);
  o2::track::TrackParCov::covMat_t o2Cov{};
  for (int i = 0; i < kBoundCovSize; ++i) {
    o2Cov[i] = static_cast<float>(cov[i]);
  }
  const o2::track::TrackParCov::params_t par{static_cast<float>((-sa * s.r[0]) + (ca * s.r[1])), static_cast<float>(s.r[2]), static_cast<float>(std::sin(phi - alpha)), static_cast<float>(std::cos(theta) / sinTheta), static_cast<float>(qOverP / sinTheta)};
  return {static_cast<float>((ca * s.r[0]) + (sa * s.r[1])), static_cast<float>(alpha), par, o2Cov};
}

/// \a state, bound to \a plane, as O2's forward parameters at the point's z:
/// x, y, phi, tanl = cot(theta), q/pt = (q/p) / sin(theta), \a chi2 the
/// track's. The field \a bzKiloGauss enters the covariance only.
[[nodiscard]] inline o2::track::TrackParCovFwd toTrackParCovFwd(const BoundTrackState& state, const TrackingSurface& plane, double bzKiloGauss, double chi2)
{
  const double e[3] = {0., 0., 1.};
  const auto s = detail::shiftOnto(state, plane, e, bzKiloGauss);
  const double theta = state.theta();
  const double sinTheta = std::sin(theta);
  const double qOverP = state.qOverP();
  const auto cov = detail::transported(detail::forwardJacobian(state, s), state);
  const o2::track::SMatrix5 par{s.r[0], s.r[1], static_cast<double>(state.phi()), std::cos(theta) / sinTheta, qOverP / sinTheta};
  // The lower triangle, row by row.
  return {s.r[2], par, o2::track::SMatrix55Sym(cov.begin(), cov.end(), true, true), chi2};
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_O2PARAMETERS_H_
