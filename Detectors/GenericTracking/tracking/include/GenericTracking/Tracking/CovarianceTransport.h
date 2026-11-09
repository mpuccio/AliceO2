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

/// \file CovarianceTransport.h
/// \brief The analytic Jacobian of a uniform-solenoid step, and the predicted
///        covariance F (P + Q) F^T + Q on the covariance's U D U^T factors.

#ifndef GTRK_TRACKING_COVARIANCETRANSPORT_H_
#define GTRK_TRACKING_COVARIANCETRANSPORT_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/Trig.h"

#include "CommonConstants/MathConstants.h"

#include <type_traits>

namespace o2::gtrk
{

/// The non-trivial part of a uniform-solenoid Jacobian: fifteen numbers, not
/// twenty-five.
///
/// In a purely longitudinal field theta and q/p are conserved, so their rows
/// are unit vectors. Those ten entries are constants.
///
/// They are omitted as a *type*, not as a convention. A field with transverse
/// components turns theta, so a field-map propagator's Jacobian is dense --
/// and handing one to a transport that assumes these zeros would be silently
/// wrong in the covariance while looking perfectly healthy. Making it
/// unrepresentable is the cheapest enforcement available.
struct SolenoidJacobian {
  Scalar loc0[kNBoundParams]{}; ///< d(loc0') / d(loc0, loc1, phi, theta, q/p)
  Scalar loc1[kNBoundParams]{};
  Scalar phi[kNBoundParams]{};
};

/// A step's process noise, as variances, in the order a filter step adds it:
/// the multiple scattering of the material nearer the surface left, before
/// the transport, and that of the rest after it (TrackFitter.h,
/// StepMaterial); then the energy-loss straggling on q/p.
struct StepNoise {
  Scalar phiBefore = 0.F;   ///< theta0^2 / sin^2(theta) of the nearer half
  Scalar thetaBefore = 0.F; ///< theta0^2 of the nearer half
  Scalar phiAfter = 0.F;    ///< the same, the farther half
  Scalar thetaAfter = 0.F;
  Scalar qOverPAfter = 0.F; ///< the straggling's q/p variance
};

/// The predicted covariance F (P + Q_before) F^T + Q_after of a step, on
/// \a factors in place: P's factors in, the prediction's out, P never formed
/// (Thornton's half of the Bierman-Thornton filter). False if \a factors are
/// not those of a positive-definite covariance, or a pivot of the result is
/// not positive (NaN included).
///
/// The prediction is a sum of weighted outer products, P' = sum_k w_k y_k
/// y_k^T. Each y_k, a *column*, is a five-vector with its variance w_k:
///
///   - 0-4, the transported factors: F u_k, u_k U's k-th column, with D_k;
///   - 5, 6, the nearer half's noise on phi and theta, transported: F e_phi
///     and F e_theta, F's phi and theta columns;
///   - 7, 8, the farther half's: e_phi and e_theta;
///   - the straggling, e_q/p.
///
/// Modified weighted Gram-Schmidt turns them into the five of P' = U' D' U'^T,
/// from q/p, the last row, to loc0: each row's D' is its weighted norm, and
/// the rows before it are made orthogonal to it (U' the projections). Every
/// D' is a sum of positive terms, so a small conditional variance -- loc1
/// given theta, once a precise measurement has fixed loc1 and a diffuse
/// theta is carried over a long lever arm -- comes out to float's relative
/// precision on itself, where forming P' would round it against the
/// variance of loc1 alone.
///
/// The structure is used: theta and q/p are conserved, so their rows of F U
/// are U's own, almost empty, and the two are reduced in closed form. q/p's
/// straggling column and the transported factor column 4 are both
/// proportional to that column after the q/p row, so they merge into one
/// column of weight D_4 s / (D_4 + s); without straggling it drops out.
GPUhdi() bool transportFactors(const SolenoidJacobian& jac, const StepNoise& noise, CovarianceFactors& factors)
{
  using detail::strictUpperIndex;
  if (!positiveDiagonal(factors)) {
    return false;
  }
  const Scalar* rows[3] = {jac.loc0, jac.loc1, jac.phi};
  // The three dense rows (loc0, loc1, phi) of the columns; column 7 is e_phi,
  // nonzero in phi's row only until the phi row fills the others.
  constexpr int kColumns = 9;
  Scalar y[3][kColumns]{};
  Scalar w[kColumns]{};
  for (int r = 0; r < 3; ++r) {
    for (int k = 0; k < kNBoundParams; ++k) {
      Scalar value = rows[r][k];
      for (int l = 0; l < k; ++l) {
        value += rows[r][l] * factors.upper[strictUpperIndex(l, k)];
      }
      y[r][k] = value;
    }
    y[r][5] = rows[r][kPhi];
    y[r][6] = rows[r][kTheta];
  }
  y[kPhi][7] = 1.F;
  for (int k = 0; k < kNBoundParams; ++k) {
    w[k] = factors.diagonal[k];
  }
  w[5] = noise.phiBefore;
  w[6] = noise.thetaBefore;
  w[7] = noise.phiAfter;
  w[8] = noise.thetaAfter;

  // q/p: its row is e_q/p in column 4 and the straggling's. Its rows above
  // project by D_4 / D_4', and column 4 keeps the merged weight.
  const Scalar thetaQOverP = factors.upper[strictUpperIndex(kTheta, kQOverP)];
  const Scalar dQOverP = factors.diagonal[kQOverP] + noise.qOverPAfter;
  const Scalar kept = factors.diagonal[kQOverP] / dQOverP;
  w[kQOverP] = noise.qOverPAfter * kept;
  factors.diagonal[kQOverP] = dQOverP;
  for (int r = 0; r < 3; ++r) {
    factors.upper[strictUpperIndex(r, kQOverP)] = y[r][kQOverP] * kept;
  }
  factors.upper[strictUpperIndex(kTheta, kQOverP)] = thetaQOverP * kept;

  // theta: its row is 1 in columns 3, 6 and 8 and U's theta-q/p entry in 4.
  const Scalar weighted4 = thetaQOverP * w[kQOverP];
  const Scalar dTheta = w[kTheta] + (thetaQOverP * weighted4) + w[6] + w[8];
  if (!(dTheta > 0.F)) {
    return false;
  }
  factors.diagonal[kTheta] = dTheta;
  const Scalar inverseTheta = 1.F / dTheta;
  for (int r = 0; r < 3; ++r) {
    const Scalar u = ((y[r][kTheta] * w[kTheta]) + (y[r][kQOverP] * weighted4) + (y[r][6] * w[6])) * inverseTheta;
    factors.upper[strictUpperIndex(r, kTheta)] = u;
    y[r][kTheta] -= u;
    y[r][kQOverP] -= u * thetaQOverP;
    y[r][6] -= u;
    y[r][8] = -u;
  }

  // phi, loc1, loc0: dense rows.
  for (int j = 2; j >= 0; --j) {
    Scalar weighted[kColumns];
    Scalar dj = 0.F;
    for (int k = 0; k < kColumns; ++k) {
      weighted[k] = w[k] * y[j][k];
      dj += y[j][k] * weighted[k];
    }
    if (!(dj > 0.F)) {
      return false;
    }
    factors.diagonal[j] = dj;
    const Scalar inverse = 1.F / dj;
    for (int i = 0; i < j; ++i) {
      Scalar dot = 0.F;
      for (int k = 0; k < kColumns; ++k) {
        dot += y[i][k] * weighted[k];
      }
      const Scalar u = dot * inverse;
      factors.upper[strictUpperIndex(i, j)] = u;
      for (int k = 0; k < kColumns; ++k) {
        y[i][k] -= u * y[j][k];
      }
    }
  }
  return true;
}

/// Analytic Jacobian of a propagation, evaluated at a **reference** state.
///
/// The reference, not the filtered state. A Kalman filter that linearises
/// about its own estimate is an *extended* filter, and its linearisation point
/// carries whatever error the filter currently has; linearising about a
/// reference trajectory - O2 does this with the `linRef` overload of
/// `PropagatorImpl::PropagateToXBxByBz` - keeps F and Q independent of the
/// measurements. For a two-pass fit that is a correctness property rather than
/// a refinement: the outward and inward passes then see *identical* F and Q
/// for the same physical step, which is what lets them be combined.
///
/// Unlike O2 we do not also transport the state through the reference. O2
/// writes `state = ref + F (state - ref)` because its exact propagation is a
/// stepped field integration and doing it once for the reference is a real
/// saving. Ours is a closed form, so linearising the state would trade
/// accuracy for a speed-up we do not need.
///
/// The derivation, for a uniform Bz. With `b = Bz * B2C` and `kp = (q/p) b`,
/// the turn angle is `psi = kp * s` - **independent of theta**, because the
/// 1/sin(theta) in the curvature cancels the sin(theta) in the transverse path.
/// That is what keeps this tractable. Writing `u = (cos phi, sin phi, 0)`,
/// `w = (-sin phi, cos phi, 0)`, `Sh = s sinc(psi)` and `Ch = s versine(psi)`,
/// the displacement at fixed path length is
///
///     dr = sin(theta) (Sh u + Ch w) + s cos(theta) zhat
///
/// and the five geometric derivatives at fixed s follow by inspection. The
/// path length itself moves with the source parameters, because the target is
/// a surface rather than a fixed s; differentiating `n . (r(s) - O) = 0` gives
/// `ds/dx = -(n . dr/dx) / (n . D)`, which is the term that couples every row.
///
/// \a T controls the transient derivative arithmetic. The propagated path,
/// input state and stored Jacobian are Scalar; ExactScalar is retained as a
/// host-side validation oracle while the device path uses the Scalar default.
template <typename T = Scalar, CrossingRule Rule = CrossingRule::kFirstAhead>
GPUdi() PropagationStatus jacobianAnalytic(const Propagator& prop, const BoundTrackState& reference, const TrackingSurface& from, const TrackingSurface& to, PropagationDirection dir, SolenoidJacobian& jac, BoundTrackState* propagatedState = nullptr, PropagationResult* propagatedResult = nullptr)
{
  using Math = o2::gpu::GPUCommonMath;

  // The reference is propagated exactly; its path length and arrival
  // direction are what everything below is evaluated at.
  BoundTrackState at = reference;
  detail::PropagationTrig direction;
  const auto res = propagate<Rule>(prop, at, from, to, dir, &direction);
  if (res.status != PropagationStatus::kOk && res.status != PropagationStatus::kOutOfBounds) {
    return res.status;
  }

  const auto s = static_cast<T>(res.pathLength);
  const auto phi = static_cast<T>(reference.phi());
  const auto qOverP = static_cast<T>(reference.qOverP());

  // The solve's own FP64 direction, narrowed: the Jacobian is evaluated where
  // the propagation went.
  const T sinPhi = static_cast<T>(direction.sinPhi);
  const T cosPhi = static_cast<T>(direction.cosPhi);
  const T sinTheta = static_cast<T>(direction.sinTheta);
  const T cosTheta = static_cast<T>(direction.cosTheta);

  const T bField = static_cast<T>(prop.bzKiloGauss) * static_cast<T>(constants::math::B2C);
  const T kp = qOverP * bField;
  const T psi = kp * s;

  const T u[3] = {cosPhi, sinPhi, T(0)};
  const T w[3] = {-sinPhi, cosPhi, T(0)};
  const T sh = s * sinc<T>(psi);
  const T ch = s * versine<T>(psi);

  T sinPhi1 = T(0);
  T cosPhi1 = T(0);
  if constexpr (std::is_same_v<T, Scalar>) {
    // Propagation already rotated this direction for its incidence test.
    // Reuse it rather than reevaluating the rounded arrival angle, on the
    // host and on the devices alike, so that the deterministic build's
    // backends agree. The higher-precision derivative oracle evaluates the
    // angle itself.
    sinPhi1 = direction.sinPhiArrival;
    cosPhi1 = direction.cosPhiArrival;
  } else {
    sincosOf(phi + psi, sinPhi1, cosPhi1);
  }
  const T dvec[3] = {sinTheta * cosPhi1, sinTheta * sinPhi1, cosTheta};

  T n1[3] = {};
  T u1[3] = {};
  T v1[3] = {};
  for (int i = 0; i < 3; ++i) {
    n1[i] = static_cast<T>(to.normal[i]);
    u1[i] = static_cast<T>(to.axis0[i]);
    v1[i] = static_cast<T>(to.axis1[i]);
  }

  T nD = T(0);
  T uD = T(0);
  T vD = T(0);
  for (int i = 0; i < 3; ++i) {
    nD += n1[i] * dvec[i];
    uD += u1[i] * dvec[i];
    vD += v1[i] * dvec[i];
  }
  if (Math::Abs(nD) < static_cast<T>(prop.minAbsCosIncidence)) {
    return PropagationStatus::kGrazing;
  }

  // Derivatives of the arrival point with respect to each source parameter,
  // taken at fixed path length.
  T dr[kNBoundParams][3] = {};
  for (int i = 0; i < 3; ++i) {
    dr[kLoc0][i] = static_cast<T>(from.axis0[i]);
    dr[kLoc1][i] = static_cast<T>(from.axis1[i]);
    dr[kPhi][i] = sinTheta * (sh * w[i] - ch * u[i]);
    dr[kTheta][i] = cosTheta * (sh * u[i] + ch * w[i]);
  }
  dr[kTheta][2] -= s * sinTheta;

  // d/d(q/p) enters only through kp, so it is dSh/dkp and dCh/dkp times b.
  const T scale = sinTheta * bField * s * s;
  for (int i = 0; i < 3; ++i) {
    dr[kQOverP][i] = scale * (dSinc<T>(psi) * u[i] + dVersine<T>(psi) * w[i]);
  }

  for (int x = 0; x < kNBoundParams; ++x) {
    T ndr = T(0);
    T udr = T(0);
    T vdr = T(0);
    for (int i = 0; i < 3; ++i) {
      ndr += n1[i] * dr[x][i];
      udr += u1[i] * dr[x][i];
      vdr += v1[i] * dr[x][i];
    }
    const T dsdx = -ndr / nD;

    jac.loc0[x] = static_cast<Scalar>(udr + uD * dsdx);
    jac.loc1[x] = static_cast<Scalar>(vdr + vD * dsdx);

    T dphi = kp * dsdx;
    if (x == kPhi) {
      dphi += T(1);
    }
    if (x == kQOverP) {
      dphi += s * bField;
    }
    jac.phi[x] = static_cast<Scalar>(dphi);
    // The theta and q/p rows are not stored: see SolenoidJacobian for why
    // they are absent from the type rather than filled with constants.
  }
  // The fit needs this exact crossing as well. Publish it only after the
  // Jacobian succeeds; propagatedState may alias reference, now fully read.
  if (propagatedState != nullptr) {
    *propagatedState = at;
  }
  if (propagatedResult != nullptr) {
    *propagatedResult = res;
  }
  return PropagationStatus::kOk;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_COVARIANCETRANSPORT_H_
