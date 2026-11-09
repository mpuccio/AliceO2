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

/// \file Oracles.h
/// \brief What the tests hold the production path to and nothing else calls.
///
///  - the covariance transport's oracles: the dense transport, the step's
///    prediction in double, the finite-difference Jacobian, which
///    `testCovarianceTransport` and `testPropagator` hold the structured,
///    factored and analytic versions to;
///  - the joint double-precision measurement update, `testMeasurementUpdate`'s
///    reference for the factorised float one;
///  - the three-point seed, `testTrackFitter`'s;
///  - the point-only triplet and neighbour overloads only the tests call
///    (Validation.h holds those the studies call).

#ifndef GTRK_TEST_ORACLES_H_
#define GTRK_TEST_ORACLES_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/CovarianceTransport.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/MeasurementUpdate.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/Seed.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/NeighbourCut.h"
#include "GenericTracking/Tracking/TripletCut.h"
#include "GenericTracking/Tracking/TripletFit.h"
#include "GenericTracking/Tracking/PairCut.h"

#include "MathUtils/detail/trigonometric.h"

#include <cmath>
#include <cstddef>
#include <numbers>

namespace o2::gtrk
{

// ---- the covariance transport's oracles ----

/// Row-major 5x5: jac[i * kNBoundParams + j] is d(out_i) / d(in_j).
inline constexpr int kJacobianSize = kNBoundParams * kNBoundParams;

/// Write a SolenoidJacobian out as a dense 5x5, filling in the constants.
/// For comparison against the finite-difference reference and for QA; the
/// transport never needs it.
GPUhdi() GPUAlwaysInline() void expandJacobian(const SolenoidJacobian& jac, Scalar dense[kJacobianSize])
{
  for (int i = 0; i < kJacobianSize; ++i) {
    dense[i] = 0.F;
  }
  for (int j = 0; j < kNBoundParams; ++j) {
    dense[kLoc0 * kNBoundParams + j] = jac.loc0[j];
    dense[kLoc1 * kNBoundParams + j] = jac.loc1[j];
    dense[kPhi * kNBoundParams + j] = jac.phi[j];
  }
  dense[kTheta * kNBoundParams + kTheta] = 1.F;
  dense[kQOverP * kNBoundParams + kQOverP] = 1.F;
}

/// Transport a packed covariance through a general dense Jacobian, in place.
///
/// Assuming nothing about the field.
/// This is the reference the factored transport is checked
/// against, and it is what a field-map propagator will need.
GPUhdi() GPUAlwaysInline() void transportCovarianceDense(const Scalar jac[kJacobianSize], Scalar cov[kBoundCovSize])
{
  Scalar m[kJacobianSize] = {}; // m = F P
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j < kNBoundParams; ++j) {
      Scalar sum = 0.F;
      for (int k = 0; k < kNBoundParams; ++k) {
        sum += jac[i * kNBoundParams + k] * cov[covIndex(k, j)];
      }
      m[i * kNBoundParams + j] = sum;
    }
  }

  Scalar out[kBoundCovSize] = {}; // out = m F^T, lower triangle only
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j <= i; ++j) {
      Scalar sum = 0.F;
      for (int k = 0; k < kNBoundParams; ++k) {
        sum += m[i * kNBoundParams + k] * jac[j * kNBoundParams + k];
      }
      out[covIndex(i, j)] = sum;
    }
  }

  for (int i = 0; i < kBoundCovSize; ++i) {
    cov[i] = out[i];
  }
}

/// A step's predicted covariance, F (P + Q_before) F^T + Q_after, formed in
/// double from the float inputs: the oracle `transportFactors` is held to.
/// \a predicted is packed as `cov`.
GPUhdi() GPUAlwaysInline() void predictCovarianceExact(const SolenoidJacobian& jac, const StepNoise& noise, const Scalar (&cov)[kBoundCovSize], ExactScalar (&predicted)[kBoundCovSize])
{
  Scalar dense[kJacobianSize];
  expandJacobian(jac, dense);
  ExactScalar p[kNBoundParams][kNBoundParams];
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j < kNBoundParams; ++j) {
      p[i][j] = static_cast<ExactScalar>(cov[covIndex(i, j)]);
    }
  }
  p[kPhi][kPhi] += static_cast<ExactScalar>(noise.phiBefore);
  p[kTheta][kTheta] += static_cast<ExactScalar>(noise.thetaBefore);
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j <= i; ++j) {
      ExactScalar sum = 0.;
      for (int k = 0; k < kNBoundParams; ++k) {
        for (int l = 0; l < kNBoundParams; ++l) {
          sum += static_cast<ExactScalar>(dense[(i * kNBoundParams) + k]) * p[k][l] * static_cast<ExactScalar>(dense[(j * kNBoundParams) + l]);
        }
      }
      predicted[covIndex(i, j)] = sum;
    }
  }
  predicted[covIndex(kPhi, kPhi)] += static_cast<ExactScalar>(noise.phiAfter);
  predicted[covIndex(kTheta, kTheta)] += static_cast<ExactScalar>(noise.thetaAfter);
  predicted[covIndex(kQOverP, kQOverP)] += static_cast<ExactScalar>(noise.qOverPAfter);
}

/// Finite-difference step sizes.
///
/// They depend on nothing: the position and angle steps are set by the float
/// noise floor of the output, and the q/p step by the geometry rather than by
/// the track, for the reason below.
///
/// Not round numbers: with a float state the propagated output carries about
/// 1e-6 cm of rounding, so a central difference has noise 1e-6 / 2h and
/// truncation of order h^2. The two cross near h ~ 1e-2 in a length, which is
/// where these sit. The resulting Jacobian is good to roughly 1e-4 relative --
/// ample to catch a wrong sign, a missing term or a transposed index in an
/// analytic replacement, and nowhere near good enough to *be* one.
GPUhdi() GPUAlwaysInline() void defaultDifferenceSteps(Scalar step[kNBoundParams])
{
  step[kLoc0] = 1e-2F;  // cm
  step[kLoc1] = 1e-2F;  // cm
  step[kPhi] = 1e-3F;   // rad
  step[kTheta] = 1e-3F; // rad

  // Absolute, not a fraction of q/p. The turn angle psi = (q/p) * Bz * B2C * s
  // is *linear* in q/p, so a relative step shrinks with momentum until the
  // perturbation moves the arrival point by less than a float can represent,
  // and the q/p column is noise.
  //
  // A fixed 1e-2 gives d(psi) between 2e-4 rad on the shortest barrel step and
  // 9e-3 rad on the longest disc one -- linear enough for a central difference
  // everywhere, and a few hundred times above the float noise floor at both
  // ends.
  step[kQOverP] = 1e-2F;
}

/// Central-difference Jacobian of a propagation, d(target state)/d(source state).
///
/// The *reference*, not the production path: it costs ten extra
/// propagations per step. Returns kOk only if all ten perturbed
/// propagations reached the plane. Out-of-bounds is accepted throughout: a
/// derivative does not care whether the intersection landed on active silicon,
/// and insisting it did would make the reference unusable near a sensor edge.
/// \a Rule as `propagate`'s.
template <CrossingRule Rule = CrossingRule::kFirstAhead>
GPUdi() PropagationStatus jacobianFiniteDifference(const Propagator& prop, const BoundTrackState& state, const TrackingSurface& from, const TrackingSurface& to, PropagationDirection dir, const Scalar step[kNBoundParams], Scalar jac[kJacobianSize])
{
  for (int j = 0; j < kNBoundParams; ++j) {
    BoundTrackState plus = state;
    BoundTrackState minus = state;
    plus.par[j] += step[j];
    minus.par[j] -= step[j];

    const auto rp = propagate<Rule>(prop, plus, from, to, dir);
    const auto rm = propagate<Rule>(prop, minus, from, to, dir);
    if (rp.status != PropagationStatus::kOk && rp.status != PropagationStatus::kOutOfBounds) {
      return rp.status;
    }
    if (rm.status != PropagationStatus::kOk && rm.status != PropagationStatus::kOutOfBounds) {
      return rm.status;
    }

    const Scalar inv = 1.F / (2.F * step[j]);
    for (int i = 0; i < kNBoundParams; ++i) {
      ExactScalar diff = static_cast<ExactScalar>(plus.par[i]) - static_cast<ExactScalar>(minus.par[i]);
      if (i == kPhi) {
        // phi lives on a circle: the difference of two normalised angles is
        // not the angle between them until it is unwrapped, and getting this
        // wrong shows up only for tracks that happen to straddle +-pi.
        math_utils::detail::bringToPMPi<ExactScalar>(diff);
      }
      jac[i * kNBoundParams + j] = static_cast<Scalar>(diff * static_cast<ExactScalar>(inv));
    }
  }
  return PropagationStatus::kOk;
}

// ---- the measurement update's oracle ----

namespace detail
{

/// Put the direction back in phi in [-pi, pi] and theta in [0, pi], on the
/// covariance itself: the factored update's `canonicalizeDirection` for the
/// oracle. Reflecting theta changes the sign of its local coordinate, so its
/// covariance row and column change sign with it.
GPUhdi() GPUAlwaysInline() void canonicalizeDirection(BoundTrackState& state)
{
  math_utils::detail::bringToPMPi<Scalar>(state.par[kTheta]);
  if (state.theta() < 0.F) {
    state.par[kTheta] = -state.theta();
    state.par[kPhi] += std::numbers::pi_v<Scalar>;
    for (int i = 0; i < kNBoundParams; ++i) {
      if (i != kTheta) {
        state.covAt(kTheta, i) = -state.covAt(kTheta, i);
      }
    }
  }
  math_utils::detail::bringToPMPi<Scalar>(state.par[kPhi]);
}

/// The chi2 of one update, from its innovation covariance `[[a b] [b d]]`
/// over (loc0, loc1) and its residuals. In \a T, the update's arithmetic.
template <typename T>
GPUhdi() GPUAlwaysInline() void spatialChi2(T a, T b, T d, T residualLoc0, T residualLoc1, MeasurementUpdateResult& result)
{
  const T det = a * d - b * b;
  result.chi2 = static_cast<Scalar>((d * residualLoc0 * residualLoc0 - T(2) * b * residualLoc0 * residualLoc1 + a * residualLoc1 * residualLoc1) / det);
  // The comparison fails for NaN. The chi2's own sign is not the test: a
  // near-zero residual's quadratic form can round to a tiny negative number
  // under a sound covariance, and that is a perfect update, not a broken one.
  result.positiveDefinite = det > T(0);
}

} // namespace detail

/// The factored update on the covariance itself: factored, updated, written
/// back. `updateWithMeasurementExact` below is the oracle the tests hold it to.
[[nodiscard]] GPUhdi() GPUAlwaysInline() MeasurementUpdateResult updateWithMeasurement(BoundTrackState& state, const Measurement& measurement)
{
  auto factors = factorCovariance(state.cov);
  const auto result = updateWithMeasurement(state, factors, measurement);
  covarianceFromFactors(factors, state.cov);
  return result;
}

/// Joint update of a bound state with (loc0, loc1), in double: the oracle the
/// factorised float `updateWithMeasurement` is held to.
///
/// The measurement covariance has one spatial cross term, so the two
/// components are updated together. Positive measurement variances and a
/// positive-definite covariance are construction invariants.
///
/// **The arithmetic is double; only the stored covariance is float.** The
/// update has to produce a posterior far smaller than its prior, and a direct
/// covariance subtraction cannot do that in float: the first step of an
/// inward pass over a ten-layer barrel
/// runs 20 cm from a diffuse 0.1 rad in phi, arrives with a 9 cm^2 position
/// variance against a 2.5e-6 cm^2 cluster, and a float update loses much of
/// the result to rounding. The stored float values are each accurate to
/// their own magnitude, which is all the next transport and update need.
[[nodiscard]] GPUhdi() GPUAlwaysInline() MeasurementUpdateResult updateWithMeasurementExact(BoundTrackState& state, const Measurement& measurement)
{
  MeasurementUpdateResult result;
  const ExactScalar residual[2] = {
    static_cast<ExactScalar>(measurement.loc0 - state.loc0()),
    static_cast<ExactScalar>(measurement.loc1 - state.loc1()),
  };
  result.residualLoc0 = static_cast<Scalar>(residual[0]);
  result.residualLoc1 = static_cast<Scalar>(residual[1]);

  const auto cov = [&state](int i, int j) { return static_cast<ExactScalar>(state.covAt(i, j)); };
  const ExactScalar innovation[2][2] = {
    {cov(kLoc0, kLoc0) + static_cast<ExactScalar>(measurement.varLoc0), cov(kLoc0, kLoc1) + static_cast<ExactScalar>(measurement.covLoc01)},
    {cov(kLoc1, kLoc0) + static_cast<ExactScalar>(measurement.covLoc01), cov(kLoc1, kLoc1) + static_cast<ExactScalar>(measurement.varLoc1)},
  };
  const ExactScalar a = innovation[0][0];
  const ExactScalar b = innovation[0][1];
  const ExactScalar d = innovation[1][1];
  const ExactScalar invDet = 1. / (a * d - b * b);
  const ExactScalar inverse[2][2] = {
    {d * invDet, -b * invDet},
    {-b * invDet, a * invDet},
  };

  constexpr int measured[2] = {kLoc0, kLoc1};
  ExactScalar gain[kNBoundParams][2]{};
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int column = 0; column < 2; ++column) {
      for (int row = 0; row < 2; ++row) {
        gain[i][column] += cov(i, measured[row]) * inverse[row][column];
      }
    }
  }

  for (int i = 0; i < kNBoundParams; ++i) {
    ExactScalar correction = 0.;
    for (int column = 0; column < 2; ++column) {
      correction += gain[i][column] * residual[column];
    }
    state.par[i] += static_cast<Scalar>(correction);
  }

  // Joseph form, written with the innovation covariance S:
  //
  // P' = (I - KH) P (I - KH)^T + K R K^T
  //    = P - KHP - PH^T K^T + K S K^T.
  //
  // This is algebraically the same update as P - KHP, but it does not obtain
  // a small posterior by subtracting one nearly equal float matrix: the
  // shorter form can leave the covariance indefinite, which shows as negative
  // chi2 in low-pT fits.
  ExactScalar updated[kBoundCovSize]{};
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j <= i; ++j) {
      ExactScalar value = cov(i, j);
      for (int measuredIndex = 0; measuredIndex < 2; ++measuredIndex) {
        value -= gain[i][measuredIndex] * cov(measured[measuredIndex], j);
        value -= cov(i, measured[measuredIndex]) * gain[j][measuredIndex];
      }
      for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 2; ++column) {
          value += gain[i][row] * innovation[row][column] * gain[j][column];
        }
      }
      updated[covIndex(i, j)] = value;
    }
  }
  for (int i = 0; i < kBoundCovSize; ++i) {
    state.cov[i] = static_cast<Scalar>(updated[i]);
  }

  detail::spatialChi2(a, b, d, residual[0], residual[1], result);

  detail::canonicalizeDirection(state);
  return result;
}

// ---- the three-point seed ----

/// The prior on a three-point seed, whose q/p is known to a factor of a few.
GPUhdi() GPUAlwaysInline() void setDiffuseSeedCovariance(BoundTrackState& seed) { setDiffuseCovariance(seed, 3.F * o2::gpu::GPUCommonMath::Abs(seed.qOverP())); }

/// Construct a state on the first of three ordered measurement surfaces.
///
/// The signed transverse curvature is twice the oriented triangle area over
/// the product of its side lengths. This remains zero for a straight track;
/// no straight-line seed branch is needed. The z slope is fitted against the
/// two transverse arc lengths from the first point.
///
/// \a T controls only the arithmetic used to construct the seed. The stored
/// state remains Scalar; ExactScalar is useful as a host-side validation
/// oracle while the default keeps the device path in Scalar.
template <typename T = Scalar, typename Prop>
GPUdi() BoundTrackState makeThreePointSeed(const Prop& prop, const Measurement& first, const Measurement& second, const Measurement& third, const TrackingSurface& firstSurface, const TrackingSurface& secondSurface, const TrackingSurface& thirdSurface)
{
  Scalar firstGlobalScalar[3]{};
  Scalar secondGlobalScalar[3]{};
  Scalar thirdGlobalScalar[3]{};
  surfaceToGlobal(firstSurface, first.loc0, first.loc1, firstGlobalScalar);
  surfaceToGlobal(secondSurface, second.loc0, second.loc1, secondGlobalScalar);
  surfaceToGlobal(thirdSurface, third.loc0, third.loc1, thirdGlobalScalar);

  T firstGlobal[3]{};
  T secondGlobal[3]{};
  T thirdGlobal[3]{};
  for (int i = 0; i < 3; ++i) {
    firstGlobal[i] = static_cast<T>(firstGlobalScalar[i]);
    secondGlobal[i] = static_cast<T>(secondGlobalScalar[i]);
    thirdGlobal[i] = static_cast<T>(thirdGlobalScalar[i]);
  }

  const T dx12 = secondGlobal[0] - firstGlobal[0];
  const T dy12 = secondGlobal[1] - firstGlobal[1];
  const T dx13 = thirdGlobal[0] - firstGlobal[0];
  const T dy13 = thirdGlobal[1] - firstGlobal[1];
  const T dx23 = thirdGlobal[0] - secondGlobal[0];
  const T dy23 = thirdGlobal[1] - secondGlobal[1];
  const T chord12 = std::sqrt(dx12 * dx12 + dy12 * dy12);
  const T chord13 = std::sqrt(dx13 * dx13 + dy13 * dy13);
  const T chord23 = std::sqrt(dx23 * dx23 + dy23 * dy23);
  const T cross = dx12 * dy13 - dy12 * dx13;
  const T curvature = T(2) * cross / (chord12 * chord13 * chord23);

  const T halfTurnSine12 = o2::gpu::GPUCommonMath::Clamp(T(0.5) * curvature * chord12, T(-1), T(1));
  T phi = atan2Of(dy12, dx12) - asinOf(halfTurnSine12);
  math_utils::detail::bringToPMPi<T>(phi);

  const T arc12 = chord12 * detail::asinOverArgument(halfTurnSine12);
  const T arc13 = detail::transverseArc(curvature, dx13, dy13);
  const T dz12 = secondGlobal[2] - firstGlobal[2];
  const T dz13 = thirdGlobal[2] - firstGlobal[2];
  const T cotTheta = (arc12 * dz12 + arc13 * dz13) / (arc12 * arc12 + arc13 * arc13);
  const T theta = atan2Of(T(1), cotTheta);
  const T sinTheta = T(1) / std::sqrt(T(1) + cotTheta * cotTheta);
  const auto& seedProp = pointPropagator(prop, secondSurface, second);
  const T field = static_cast<T>(seedProp.bzKiloGauss) * static_cast<T>(constants::math::B2C);

  BoundTrackState seed;
  seed.par[kLoc0] = first.loc0;
  seed.par[kLoc1] = first.loc1;
  seed.par[kPhi] = static_cast<Scalar>(phi);
  seed.par[kTheta] = static_cast<Scalar>(theta);
  seed.par[kQOverP] = static_cast<Scalar>(curvature * sinTheta / field);
  seed.surface = first.surface;
  seed.surfacePhi = firstSurface.surfacePhi;
  setDiffuseSeedCovariance(seed);
  return seed;
}

// ---- the triplet and neighbour fits from the points alone ----

/// The triplet fit from the points, the kink terms discarded.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool fitTriplet(const PairPoint& first, const PairPoint& middle, const PairPoint& third, LayerID middleLayer, const ScatteringTable<Det>& table, const PairCut& cut, TripletFitResult& result)
{
  KinkTerms terms;
  return fitTriplet(first, middle, third, middleLayer, table, cut, result, terms);
}

/// The adjacent fit from the four points and given variances at `b` and `c`,
/// the kink terms computed here.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool fitAdjacentTriplets(const PairPoint& a, const PairPoint& b, const PairPoint& c, const PairPoint& d, Scalar varianceSecond, Scalar varianceThird, TripletFitResult& result)
{
  KinkTerms first;
  KinkTerms second;
  if (!makeKinkTerms(a, b, c, first) || !makeKinkTerms(b, c, d, second)) {
    return false;
  }
  return fitAdjacentTriplets(first, second, a, b, c, d, varianceSecond, varianceThird, result);
}

/// The adjacent fit from the four points, with the scattering at `b` and `c`
/// at pT_min from the pair cut's table, each at the polar angle of its own
/// triplet; the kink terms computed here.
template <DetectorLayout Det>
[[nodiscard]] inline bool fitAdjacentTriplets(const PairPoint& a, const PairPoint& b, const PairPoint& c, const PairPoint& d, LayerID secondLayer, LayerID thirdLayer, const ScatteringTable<Det>& table, const PairCut& cut, TripletFitResult& result)
{
  KinkTerms first;
  KinkTerms second;
  if (!makeKinkTerms(a, b, c, first) || !makeKinkTerms(b, c, d, second)) {
    return false;
  }
  return fitAdjacentTriplets(first, second, a, b, c, d, secondLayer, thirdLayer, table, cut, result);
}

/// The triplet cuts after the dip angle, the kink terms and the kink flag
/// discarded.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool tripletCompatible(const IterationParameters<Det>& parameters, LayerID middleLayer, const TimedPoint& first, const TimedPoint& middle, const TimedPoint& third, TripletFitResult& fit)
{
  KinkTerms terms;
  bool kinked = false;
  return tripletCompatible(parameters, middleLayer, first, middle, third, fit, terms, kinked);
}

/// The neighbour cuts after the curvature agreement, the stage's own
/// (NeighbourCut.h), the two triplets' kink terms made there from the points,
/// no kink allowance.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool neighbourCompatible(const IterationParameters<Det>& parameters, LayerID secondLayer, LayerID thirdLayer, const TimedPoint& a, const TimedPoint& b, const TimedPoint& c, const TimedPoint& d, TripletFitResult& fit)
{
  KinkTerms bcd;
  bool made = false;
  bool valid = false;
  return neighbourCompatible(parameters, secondLayer, thirdLayer, a, b, c, d, bcd, made, valid, fit, 1.F, 1.F);
}

} // namespace o2::gtrk

#endif // GTRK_TEST_ORACLES_H_
