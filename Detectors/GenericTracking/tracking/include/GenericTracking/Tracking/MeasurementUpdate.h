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

/// \file MeasurementUpdate.h
/// \brief Kalman update of a bound track state with one cluster measurement.

#ifndef GTRK_TRACKING_MEASUREMENTUPDATE_H_
#define GTRK_TRACKING_MEASUREMENTUPDATE_H_

#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/Measurement.h"

#include "MathUtils/detail/trigonometric.h"

#include <numbers>

namespace o2::gtrk
{

/// Innovations and the fit-quality increment produced by one measurement.
///
/// A measurement's time is its readout frame's box, which carries no
/// information about where in the frame the particle passed: the update takes
/// the two local coordinates alone.
struct MeasurementUpdateResult {
  Scalar residualLoc0 = 0.F;
  Scalar residualLoc1 = 0.F;
  Scalar chi2 = 0.F; ///< of (loc0, loc1), two degrees of freedom
  /// Whether the innovation covariance was positive definite. False, and the
  /// chi2 above is meaningless, when the predicted covariance was not
  /// positive definite in float (`updateUsable`).
  bool positiveDefinite = false;
};

namespace detail
{

/// Put the direction back in phi in [-pi, pi] and theta in [0, pi], the
/// covariance as factors: reflecting theta is S P S with S = diag(1, 1, 1,
/// -1, 1), whose factors are S U S and D -- U's entries in theta's row and
/// column change sign. Phi wrapping and the accompanying half-turn are
/// constant shifts and leave the covariance alone.
GPUhdi() GPUAlwaysInline() void canonicalizeDirection(BoundTrackState& state, CovarianceFactors& factors)
{
  math_utils::detail::bringToPMPi<Scalar>(state.par[kTheta]);
  if (state.theta() < 0.F) {
    state.par[kTheta] = -state.theta();
    state.par[kPhi] += std::numbers::pi_v<Scalar>;
    for (int i = 0; i < kTheta; ++i) {
      factors.upper[strictUpperIndex(i, kTheta)] = -factors.upper[strictUpperIndex(i, kTheta)];
    }
    for (int j = kTheta + 1; j < kNBoundParams; ++j) {
      factors.upper[strictUpperIndex(kTheta, j)] = -factors.upper[strictUpperIndex(kTheta, j)];
    }
  }
  math_utils::detail::bringToPMPi<Scalar>(state.par[kPhi]);
}

GPUhdi() GPUAlwaysInline() void applySpatialGain(BoundTrackState& state, const Scalar (&gain)[kNBoundParams], Scalar residual)
{
  for (int i = 0; i < kNBoundParams; ++i) {
    state.par[i] += gain[i] * residual;
  }
}

} // namespace detail

/// A covariance's factors, and the covariance from them (CovarianceFactors).
[[nodiscard]] GPUhdi() GPUAlwaysInline() CovarianceFactors factorCovariance(const Scalar (&cov)[kBoundCovSize])
{
  CovarianceFactors factors;
  detail::factorUD(cov, factors.upper, factors.diagonal);
  return factors;
}

GPUhdi() GPUAlwaysInline() void covarianceFromFactors(const CovarianceFactors& factors, Scalar (&cov)[kBoundCovSize])
{
  detail::covarianceFromUD(factors.upper, factors.diagonal, cov);
}

/// Update of a bound state with (loc0, loc1), its covariance as factors: the
/// factors themselves are updated (Bierman), and `state.cov` is not touched.
///
/// The spatial measurement covariance is first whitened into two independent
/// scalar observations: R = L D L^T, and applying L^-1 to measurement and
/// model removes covLoc01 without rotating the track state. Sequential scalar
/// updates are then algebraically the joint update, and the joint chi2 is the
/// sum of the two scalar ones, each residual squared over its innovation
/// variance. A precise measurement produces its small posterior variance by
/// multiplication and division instead of subtracting `K H P` from `P`.
///
/// With \a Determinant, \a determinant gets the innovation covariance's
/// determinant over (loc0, loc1), det(H P H^T + R) [cm^4]: the product of the
/// two innovations.
template <bool Determinant = false>
[[nodiscard]] GPUhdi() GPUAlwaysInline() MeasurementUpdateResult updateWithMeasurement(BoundTrackState& state, CovarianceFactors& factors,
                                                                                       const Measurement& measurement,
                                                                                       [[maybe_unused]] Scalar* determinant = nullptr)
{
  MeasurementUpdateResult result;
  result.residualLoc0 = measurement.loc0 - state.loc0();
  result.residualLoc1 = measurement.loc1 - state.loc1();

  Scalar observation[kNBoundParams]{};
  Scalar gain[kNBoundParams]{};
  observation[kLoc0] = 1.F;
  const Scalar residual0 = result.residualLoc0;
  const Scalar innovation0 = detail::updateUD(factors.upper, factors.diagonal, observation, measurement.varLoc0, gain);
  detail::applySpatialGain(state, gain, residual0);

  const Scalar correlation = measurement.covLoc01 / measurement.varLoc0;
  const Scalar conditionalLoc1Variance = measurement.varLoc1 - correlation * measurement.covLoc01;
  observation[kLoc0] = -correlation;
  observation[kLoc1] = 1.F;
  const Scalar residual1 = (measurement.loc1 - correlation * measurement.loc0) - (state.loc1() - correlation * state.loc0());
  const Scalar innovation1 = detail::updateUD(factors.upper, factors.diagonal, observation, conditionalLoc1Variance, gain);
  detail::applySpatialGain(state, gain, residual1);

  result.chi2 = (residual0 * residual0 / innovation0) + (residual1 * residual1 / innovation1);
  // The whitening has determinant one and the second observation is
  // conditioned on the first: their innovations multiply to det V.
  if constexpr (Determinant) {
    *determinant = innovation0 * innovation1;
  }
  // The comparisons fail for NaN. The chi2's own sign is not the test: a
  // near-zero residual's quadratic form can round to a tiny negative number
  // under a sound covariance, and that is a perfect update, not a broken one.
  result.positiveDefinite = innovation0 > 0.F && innovation1 > 0.F;
  detail::canonicalizeDirection(state, factors);
  return result;
}

/// Whether \a update was a measurement of the track the posterior \a factors
/// describe.
///
/// A predicted covariance that is not positive definite -- a prior's float P
/// can be, factors that start from one carry it -- gives an innovation that is
/// not positive, a chi2 that is negative or NaN, or a posterior with a
/// variance that is not positive (NaN included). Either way the update is not
/// a measurement of this track: a filter pass fails on it, and a road fit
/// that let a negative chi2 through its cut would rank it first. Positive D is a positive-definite posterior, every
/// variance a sum of positive terms. The comparisons are written so that NaN
/// fails them.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool updateUsable(const MeasurementUpdateResult& update, const CovarianceFactors& factors)
{
  // NaN, spelled without a self-comparison or a `Finite` that is device-only
  // (and a constant under fast math). An infinite chi2 passes here and fails
  // every cut it meets; a tiny negative one from rounding is a perfect update.
  const bool chi2IsNaN = !(update.chi2 >= 0.F) && !(update.chi2 < 0.F);
  if (!update.positiveDefinite || chi2IsNaN) {
    return false;
  }
  return positiveDiagonal(factors);
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_MEASUREMENTUPDATE_H_
