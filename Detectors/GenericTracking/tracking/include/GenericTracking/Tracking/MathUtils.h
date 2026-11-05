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

/// \file MathUtils.h
/// \brief The tracking's numeric helpers: stable sinc, versine and asin(x)/x
///        forms, the U-D covariance factorisation, the chi2 quantile, the
///        curvature radius.

#ifndef GTRK_TRACKING_MATHUTILS_H_
#define GTRK_TRACKING_MATHUTILS_H_

#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Trig.h"

#include "GPUCommonMath.h"
#include "CommonConstants/MathConstants.h"
#include "MathUtils/detail/basicMath.h"

#include <cmath>
#include <numbers>

namespace o2::gtrk
{

/// sin(x)/x, with its value at zero.
///
/// The series and quotient agree within the scalar precision at the branch.
/// The quotient itself does not cancel, but needs the series at its removable
/// singularity.
template <typename T>
GPUhdi() GPUAlwaysInline() T sinc(T x)
{
  if (o2::gpu::GPUCommonMath::Abs(x) < T(1e-4)) {
    return T(1) - x * x / T(6);
  }
  return sinOf(x) / x;
}

/// (1 - cos(x)) / x without forming the cancelling subtraction.
template <typename T>
GPUhdi() GPUAlwaysInline() T versine(T x)
{
  const T h = sinc<T>(x / T(2));
  return T(0.5) * x * h * h;
}

/// d/dx of sinc, with a series where cos(x) - sinc(x) cancels.
template <typename T>
GPUhdi() GPUAlwaysInline() T dSinc(T x)
{
  if (o2::gpu::GPUCommonMath::Abs(x) < T(0.05)) {
    const T x2 = x * x;
    return -x / T(3) + x * x2 / T(30) - x * x2 * x2 / T(840);
  }
  return (cosOf(x) - sinc<T>(x)) / x;
}

/// d/dx of versine, with a series around zero.
template <typename T>
GPUhdi() GPUAlwaysInline() T dVersine(T x)
{
  if (o2::gpu::GPUCommonMath::Abs(x) < T(0.05)) {
    const T x2 = x * x;
    return T(0.5) - x2 / T(8) + x2 * x2 / T(144) - x2 * x2 * x2 / T(5760);
  }
  return (sinOf(x) - versine<T>(x)) / x;
}

namespace detail
{

/// Column-major strict-upper index. The caller supplies first < second.
GPUhdi() GPUAlwaysInline() constexpr int strictUpperIndex(int first, int second) { return (second * (second - 1)) / 2 + first; }

/// Factor a packed symmetric positive-definite matrix as U D U^T. U is unit
/// upper triangular, so only its strict upper triangle is stored. Working
/// from the last column towards the first keeps differently scaled state
/// coordinates in D instead of mixing their scales in a square root.
template <typename T>
GPUhdi() GPUAlwaysInline() void factorUD(const T (&matrix)[kBoundCovSize], T (&upper)[kBoundUpperSize], T (&diagonal)[kNBoundParams])
{
  for (int column = kNBoundParams - 1; column >= 0; --column) {
    T pivot = matrix[covIndex(column, column)];
    for (int k = column + 1; k < kNBoundParams; ++k) {
      const T value = upper[strictUpperIndex(column, k)];
      pivot -= value * value * diagonal[k];
    }
    diagonal[column] = pivot;

    for (int row = 0; row < column; ++row) {
      T value = matrix[covIndex(row, column)];
      for (int k = column + 1; k < kNBoundParams; ++k) {
        value -= upper[strictUpperIndex(row, k)] * upper[strictUpperIndex(column, k)] * diagonal[k];
      }
      upper[strictUpperIndex(row, column)] = value / pivot;
    }
  }
}

/// Bierman's scalar measurement update of a U D U^T covariance.
///
/// The posterior is formed by updating its factors rather than subtracting
/// K H P from P. This is the important numerical property when a precise
/// cluster follows a diffuse track prior: the small posterior variance is
/// produced by multiplication and division, never as the difference of two
/// nearly equal floats. The gain updates the state mean; the return value
/// is the innovation variance, `H P H^T + variance`, for the chi2.
template <typename T>
GPUhdi() GPUAlwaysInline() T updateUD(T (&upper)[kBoundUpperSize], T (&diagonal)[kNBoundParams], const T (&measurement)[kNBoundParams], T variance,
                                      T (&gain)[kNBoundParams])
{
  T projected[kNBoundParams]{};
  T weighted[kNBoundParams]{};
  for (int column = 0; column < kNBoundParams; ++column) {
    projected[column] = measurement[column];
    for (int row = 0; row < column; ++row) {
      projected[column] += upper[strictUpperIndex(row, column)] * measurement[row];
    }
    weighted[column] = diagonal[column] * projected[column];
  }

  T numerator[kNBoundParams]{};
  T innovation = variance;
  for (int column = 0; column < kNBoundParams; ++column) {
    const T previousInnovation = innovation;
    innovation += projected[column] * weighted[column];
    const T factor = projected[column] / previousInnovation;
    diagonal[column] *= previousInnovation / innovation;

    for (int row = 0; row < column; ++row) {
      const int index = strictUpperIndex(row, column);
      const T value = upper[index];
      upper[index] = value - factor * numerator[row];
      numerator[row] += weighted[column] * value;
    }
    numerator[column] += weighted[column];
  }

  for (int row = 0; row < kNBoundParams; ++row) {
    gain[row] = numerator[row] / innovation;
  }
  return innovation;
}

/// Reconstruct a packed symmetric covariance from its U D U^T factors.
template <typename T>
GPUhdi() GPUAlwaysInline() void covarianceFromUD(const T (&upper)[kBoundUpperSize], const T (&diagonal)[kNBoundParams], T (&matrix)[kBoundCovSize])
{
  for (int row = 0; row < kNBoundParams; ++row) {
    for (int col = 0; col <= row; ++col) {
      T value = diagonal[row];
      if (col < row) {
        value *= upper[strictUpperIndex(col, row)];
      }
      for (int k = row + 1; k < kNBoundParams; ++k) {
        value += upper[strictUpperIndex(row, k)] * diagonal[k] * upper[strictUpperIndex(col, k)];
      }
      matrix[covIndex(row, col)] = value;
    }
  }
}

/// asin(x) / x, evaluated without a singularity or an unnecessary asin near zero.
template <typename T>
GPUhdi() GPUAlwaysInline() T asinOverArgument(T x)
{
  if (o2::gpu::GPUCommonMath::Abs(x) < T(0.05)) {
    const T x2 = x * x;
    return T(1) + x2 * (T(1) / T(6) + x2 * T(3) / T(40));
  }
  return asinOf(x) / x;
}

template <typename T>
struct AsinOverArgument {
  T value;
  T derivative;
};

/// `asinOverArgument(x)` and its derivative with one asin; the value is
/// `asinOverArgument(x)` to the bit.
template <typename T>
GPUhdi() GPUAlwaysInline() AsinOverArgument<T> asinOverArgumentWithDerivative(T x)
{
  if (o2::gpu::GPUCommonMath::Abs(x) < T(0.05)) {
    const T x2 = x * x;
    return {.value = T(1) + x2 * (T(1) / T(6) + x2 * T(3) / T(40)), .derivative = x * (T(1) / T(3) + x2 * (T(3) / T(10) + x2 * T(15) / T(56)))};
  }
  const T arcsine = asinOf(x);
  return {.value = arcsine / x, .derivative = ((x / o2::math_utils::detail::sqrt(T(1) - (x * x))) - arcsine) / (x * x)};
}

/// Signed-curvature chord length converted to positive transverse arc length.
template <typename T>
GPUhdi() GPUAlwaysInline() T transverseArc(T curvature, T dx, T dy)
{
  const T chord = o2::math_utils::detail::sqrt(dx * dx + dy * dy);
  // A chord of one diameter is half a turn; rounding must not push past it.
  const T halfTurnSine = o2::gpu::GPUCommonMath::Clamp(T(0.5) * curvature * chord, T(-1), T(1));
  return chord * asinOverArgument(halfTurnSine);
}

/// The regularised lower incomplete gamma function P(a, x): its series below
/// x = a + 1, and its continued fraction above (Numerical Recipes, 6.2).
[[nodiscard]] inline double regularizedGammaP(double a, double x)
{
  if (!(x > 0.)) {
    return 0.;
  }
  const double logPrefactor = (a * std::log(x)) - x - std::lgamma(a);
  if (x < a + 1.) {
    double term = 1. / a;
    double sum = term;
    for (int n = 1; n < 1000 && std::abs(term) > std::abs(sum) * 1.e-16; ++n) {
      term *= x / (a + n);
      sum += term;
    }
    return sum * std::exp(logPrefactor);
  }
  constexpr double kTiny = 1.e-300;
  double b = x + 1. - a;
  double c = 1. / kTiny;
  double d = 1. / b;
  double h = d;
  for (int n = 1; n < 1000; ++n) {
    const double an = -n * (n - a);
    b += 2.;
    d = (an * d) + b;
    d = std::abs(d) < kTiny ? kTiny : d;
    c = b + (an / c);
    c = std::abs(c) < kTiny ? kTiny : c;
    d = 1. / d;
    const double step = d * c;
    h *= step;
    if (std::abs(step - 1.) < 1.e-16) {
      break;
    }
  }
  return 1. - (std::exp(logPrefactor) * h);
}

} // namespace detail

/// The chi2 with ndf degrees of freedom that a fit passes at the confidence
/// of nSigma of a Gaussian, two-sided. The distribution P(ndf/2, x/2) is
/// inverted by bisection. ROOT's `chisquared_quantile` has this contract; not
/// in a header the device units read.
[[nodiscard]] inline double chiSquareQuantile(int ndf, double nSigma)
{
  const double confidence = std::erf(nSigma / std::numbers::sqrt2);
  const double half = 0.5 * static_cast<double>(ndf);
  const auto distribution = [half](double x) { return detail::regularizedGammaP(half, 0.5 * x); };
  double low = 0.;
  double high = 1.;
  while (distribution(high) < confidence) {
    low = high;
    high *= 2.;
  }
  for (int i = 0; i < 100; ++i) {
    const double middle = 0.5 * (low + high);
    (distribution(middle) < confidence ? low : high) = middle;
  }
  return 0.5 * (low + high);
}

/// Transverse radius of curvature [cm] for a particle with the supplied
/// transverse momentum in a solenoid field.
///
/// B2C is negative by convention, so only the magnitude enters the radius.
template <typename T>
GPUhdi() GPUAlwaysInline() T curvatureRadius(T minPtGeV, T bzKiloGauss)
{
  const T perPt = o2::gpu::GPUCommonMath::Abs(static_cast<T>(o2::constants::math::B2C) * bzKiloGauss);
  return perPt > T(0) ? minPtGeV / perPt : T(0);
}

/// Transverse arc from the origin to a point at the supplied radius.
template <typename T>
GPUhdi() GPUAlwaysInline() T transverseArcAtRadius(T curvature, T radius)
{
  return curvature > T(0) ? detail::transverseArc(curvature, radius, T(0)) : radius;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_MATHUTILS_H_
