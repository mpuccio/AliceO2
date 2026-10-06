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

#ifndef ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_
#define ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/GlobalMeasurement.h"

namespace o2::itsmft::tracking
{

struct TripletKinkVector {
  float theta{0.f};
  float phi{0.f};
};

// Theta and phi rows of the hit-coordinate Jacobian H.
struct TripletHitJacobian {
  std::array<float, 3> theta{};
  std::array<float, 3> phi{};
};

// Linearized local-triplet factor from Eq. (19) of the General Triplet Track
// Fit. H is evaluated at kappaRef = -Psi_phi / rho_phi and hit slot i maps to
// Triplet::getClusterReference(i). Measurement and MS covariances are added
// when adjacent triplets are compared.
struct TripletFitFactor {
  TripletKinkVector psi{};
  TripletKinkVector rho{};
  std::array<TripletHitJacobian, 3> h{};

  GPUhdi() bool isValid() const noexcept
  {
    return rho.phi != 0.f;
  }
};

static_assert(std::is_standard_layout_v<TripletFitFactor>);
static_assert(std::is_trivially_copyable_v<TripletFitFactor>);
static_assert(sizeof(TripletFitFactor) == 88);

struct AdjacentTripletFitResult {
  float curvature{0.f};
  float curvatureVariance{0.f};
  float chi2{0.f};
};

// Shared host/device triplet geometry and automatic derivatives.
namespace triplet_detail
{
constexpr int NCoordinates = 9;
// Coordinates are numbered x0,y0,x1,y1,x2,y2,z0,z1,z2: transverse quantities
// then depend on [0,6) only and z differences on [6,9) only.
GPUhdi() constexpr int coordinateOf(int hit, int coordinate) noexcept
{
  return coordinate < 2 ? 2 * hit + coordinate : 6 + hit;
}

// Single-precision value with the derivatives it can have: those with respect
// to the coordinates [Begin, End); every other derivative is zero and is
// neither stored nor computed. Each stored derivative is computed with the
// same operations as by a dense dual number, so the results are identical;
// one with respect to a coordinate that neither operand depends on is zero.
template <int Begin, int End>
struct DualNumber {
  static constexpr int First = Begin, Last = End;
  float value{0.f};
  std::array<float, (End > Begin ? End - Begin : 0)> derivative{};

  GPUhdi() static constexpr bool has(int i) noexcept { return i >= Begin && i < End; }
  GPUhdi() float at(int i) const noexcept { return derivative[i - Begin]; }
};
using Constant = DualNumber<0, 0>;

GPUhdi() constexpr Constant constant(float value) noexcept { return Constant{value}; }

template <int I>
GPUhdi() DualNumber<I, I + 1> variable(float value) noexcept
{
  DualNumber<I, I + 1> result{value};
  result.derivative[0] = 1.f;
  return result;
}

// The coordinates on which a result of two operands can depend.
template <typename L, typename R>
using Joint = DualNumber<(L::Last <= L::First) ? R::First : ((R::Last <= R::First) ? L::First : (L::First < R::First ? L::First : R::First)),
                         (L::Last <= L::First) ? R::Last : ((R::Last <= R::First) ? L::Last : (L::Last > R::Last ? L::Last : R::Last))>;

template <typename L, typename R>
GPUhdi() Joint<L, R> operator+(const L& lhs, const R& rhs) noexcept
{
  using T = Joint<L, R>;
  T result{lhs.value + rhs.value};
  for (int i = T::First; i < T::Last; ++i) {
    result.derivative[i - T::First] = L::has(i) && R::has(i) ? lhs.at(i) + rhs.at(i) : (L::has(i) ? lhs.at(i) : (R::has(i) ? rhs.at(i) : 0.f));
  }
  return result;
}

template <typename L, typename R>
GPUhdi() Joint<L, R> operator-(const L& lhs, const R& rhs) noexcept
{
  using T = Joint<L, R>;
  T result{lhs.value - rhs.value};
  for (int i = T::First; i < T::Last; ++i) {
    result.derivative[i - T::First] = L::has(i) && R::has(i) ? lhs.at(i) - rhs.at(i) : (L::has(i) ? lhs.at(i) : (R::has(i) ? -rhs.at(i) : 0.f));
  }
  return result;
}

template <int B, int E>
GPUhdi() DualNumber<B, E> operator-(const DualNumber<B, E>& value) noexcept
{
  DualNumber<B, E> result{-value.value};
  for (int i = 0; i < E - B; ++i) {
    result.derivative[i] = -value.derivative[i];
  }
  return result;
}

template <typename L, typename R>
GPUhdi() Joint<L, R> operator*(const L& lhs, const R& rhs) noexcept
{
  using T = Joint<L, R>;
  T result{lhs.value * rhs.value};
  for (int i = T::First; i < T::Last; ++i) {
    result.derivative[i - T::First] = L::has(i) && R::has(i) ? lhs.at(i) * rhs.value + lhs.value * rhs.at(i)
                                                             : (L::has(i) ? lhs.at(i) * rhs.value : (R::has(i) ? lhs.value * rhs.at(i) : 0.f));
  }
  return result;
}

template <typename L, typename R>
GPUhdi() Joint<L, R> operator/(const L& lhs, const R& rhs) noexcept
{
  using T = Joint<L, R>;
  const float inverse = 1.f / rhs.value;
  T result{lhs.value * inverse};
  for (int i = T::First; i < T::Last; ++i) {
    result.derivative[i - T::First] = L::has(i) && R::has(i) ? (lhs.at(i) - result.value * rhs.at(i)) * inverse
                                                             : (L::has(i) ? lhs.at(i) * inverse : (R::has(i) ? (-(result.value * rhs.at(i))) * inverse : 0.f));
  }
  return result;
}

template <int B, int E>
GPUhdi() DualNumber<B, E> squareRoot(const DualNumber<B, E>& argument) noexcept
{
  const float root = std::sqrt(argument.value);
  DualNumber<B, E> result{root};
  const float scale = 0.5f / root;
  for (int i = 0; i < E - B; ++i) {
    result.derivative[i] = scale * argument.derivative[i];
  }
  return result;
}

template <int B, int E>
GPUhdi() DualNumber<B, E> arcSine(const DualNumber<B, E>& argument) noexcept
{
  DualNumber<B, E> result{std::asin(argument.value)};
  const float scale = 1.f / std::sqrt(1.f - argument.value * argument.value);
  for (int i = 0; i < E - B; ++i) {
    result.derivative[i] = scale * argument.derivative[i];
  }
  return result;
}

template <typename Y, typename X>
GPUhdi() Joint<Y, X> arcTangent2(const Y& y, const X& x) noexcept
{
  using T = Joint<Y, X>;
  T result{std::atan2(y.value, x.value)};
  const float denominator = x.value * x.value + y.value * y.value;
  for (int i = T::First; i < T::Last; ++i) {
    result.derivative[i - T::First] = Y::has(i) && X::has(i) ? (x.value * y.at(i) - y.value * x.at(i)) / denominator
                                                             : (Y::has(i) ? (x.value * y.at(i)) / denominator : (X::has(i) ? (-(y.value * x.at(i))) / denominator : 0.f));
  }
  return result;
}

// The same value with derivatives over all coordinates.
template <typename V>
GPUhdi() DualNumber<0, NCoordinates> full(const V& value) noexcept
{
  DualNumber<0, NCoordinates> result{value.value};
  for (int i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = V::has(i) ? value.at(i) : 0.f;
  }
  return result;
}

template <typename Angle, typename Length, typename Cotangent, typename Sine, typename Cosine, typename Index>
struct SegmentGeometry {
  bool valid{false};
  Angle bendingAngle;
  Length transverseArcLength;
  Cotangent cotangentTheta;
  Sine sineTheta;
  Cosine cosineTheta;
  Index index;
  Index oneMinusIndex;
};

template <typename Curvature, typename Chord, typename DeltaZ>
GPUhdi() auto makeSegmentGeometry(const Curvature& transverseCurvature, const Chord& chordLength, const DeltaZ& deltaZ) noexcept
{
  const auto halfSine = constant(0.5f) * transverseCurvature * chordLength;
  const auto halfSine2 = halfSine * halfSine;
  const auto halfSine4 = halfSine2 * halfSine2;
  std::remove_const_t<decltype(halfSine)> asinOverArgument;
  std::remove_const_t<decltype(halfSine)> oneMinusAngleCotangent;
  const auto halfAngle = arcSine(halfSine);
  if (std::abs(halfSine.value) < 0.05f) {
    asinOverArgument = constant(1.f) + halfSine2 * constant(1.f / 6.f) + halfSine4 * constant(3.f / 40.f);
    oneMinusAngleCotangent = halfSine2 * (constant(1.f / 3.f) + halfSine2 * (constant(2.f / 15.f) + halfSine2 * constant(8.f / 105.f)));
  } else {
    asinOverArgument = halfAngle / halfSine;
    oneMinusAngleCotangent = constant(1.f) - halfAngle * squareRoot(constant(1.f) - halfSine2) / halfSine;
  }

  const auto bendingAngle = constant(2.f) * halfAngle;
  const auto transverseArcLength = chordLength * asinOverArgument;
  const auto cotangentTheta = deltaZ / transverseArcLength;
  const auto sineTheta = constant(1.f) / squareRoot(constant(1.f) + cotangentTheta * cotangentTheta);
  const auto cosineTheta = cotangentTheta * sineTheta;
  // sin²(theta) + cos²(theta) = 1. Keep 1-index explicitly to avoid
  // subtracting nearly equal floats and amplifying their error by 1/curvature.
  const auto delta = oneMinusAngleCotangent * sineTheta * sineTheta;
  const auto index = constant(1.f) / (constant(1.f) - delta);
  const auto oneMinusIndex = -delta * index;
  const bool valid = std::abs(halfSine.value) < 1.f && !(transverseArcLength.value <= 0.f || sineTheta.value <= 0.f || index.value <= 0.f);
  return SegmentGeometry<std::remove_const_t<decltype(bendingAngle)>, std::remove_const_t<decltype(transverseArcLength)>,
                         std::remove_const_t<decltype(cotangentTheta)>, std::remove_const_t<decltype(sineTheta)>,
                         std::remove_const_t<decltype(cosineTheta)>, std::remove_const_t<decltype(index)>>{
    valid, bendingAngle, transverseArcLength, cotangentTheta, sineTheta, cosineTheta, index, oneMinusIndex};
}

struct TripletGeometry {
  DualNumber<0, NCoordinates> phiTilde;
  DualNumber<0, NCoordinates> thetaTilde;
  DualNumber<0, NCoordinates> rhoPhi;
  DualNumber<0, NCoordinates> rhoTheta;
};

GPUhdi() bool makeTripletGeometry(const std::array<GlobalMeasurement, 3>& measurements,
                                  TripletGeometry& result) noexcept
{
  const auto x0 = variable<coordinateOf(0, 0)>(measurements[0].x);
  const auto y0 = variable<coordinateOf(0, 1)>(measurements[0].y);
  const auto z0 = variable<coordinateOf(0, 2)>(measurements[0].z);
  const auto x1 = variable<coordinateOf(1, 0)>(measurements[1].x);
  const auto y1 = variable<coordinateOf(1, 1)>(measurements[1].y);
  const auto z1 = variable<coordinateOf(1, 2)>(measurements[1].z);
  const auto x2 = variable<coordinateOf(2, 0)>(measurements[2].x);
  const auto y2 = variable<coordinateOf(2, 1)>(measurements[2].y);
  const auto z2 = variable<coordinateOf(2, 2)>(measurements[2].z);

  const auto dx01 = x1 - x0;
  const auto dy01 = y1 - y0;
  const auto dz01 = z1 - z0;
  const auto dx12 = x2 - x1;
  const auto dy12 = y2 - y1;
  const auto dz12 = z2 - z1;
  const auto dx02 = x2 - x0;
  const auto dy02 = y2 - y0;
  const auto length01 = squareRoot(dx01 * dx01 + dy01 * dy01);
  const auto length12 = squareRoot(dx12 * dx12 + dy12 * dy12);
  const auto length02 = squareRoot(dx02 * dx02 + dy02 * dy02);
  if (length01.value <= 0.f ||
      length12.value <= 0.f || length02.value <= 0.f) {
    return false;
  }

  const auto cross = dx01 * dy12 - dy01 * dx12;
  const auto transverseCurvature = constant(2.f) * cross / (length01 * length12 * length02);
  const auto firstSegment = makeSegmentGeometry(transverseCurvature, length01, dz01);
  if (!firstSegment.valid) {
    return false;
  }
  const auto secondSegment = makeSegmentGeometry(transverseCurvature, length12, dz12);
  if (!secondSegment.valid) {
    return false;
  }

  const auto theta01 = arcTangent2(firstSegment.transverseArcLength, dz01);
  const auto theta12 = arcTangent2(secondSegment.transverseArcLength, dz12);
  const auto phiTilde = constant(0.5f) *
                        (firstSegment.bendingAngle * firstSegment.index +
                         secondSegment.bendingAngle * secondSegment.index);
  const auto thetaTilde = theta12 - theta01 +
                          secondSegment.oneMinusIndex * secondSegment.cotangentTheta -
                          firstSegment.oneMinusIndex * firstSegment.cotangentTheta;
  const auto rhoPhi = constant(-0.5f) *
                      (firstSegment.transverseArcLength * firstSegment.index / firstSegment.sineTheta +
                       secondSegment.transverseArcLength * secondSegment.index / secondSegment.sineTheta);

  DualNumber<0, NCoordinates> rhoTheta;
  const float maximumHalfSine = 0.5f * std::abs(transverseCurvature.value) *
                                std::max(length01.value, length12.value);
  if (maximumHalfSine < 1.e-4f) {
    rhoTheta = full(transverseCurvature *
                    (length12 * length12 * secondSegment.cosineTheta -
                     length01 * length01 * firstSegment.cosineTheta) /
                    constant(12.f));
  } else {
    rhoTheta = full((firstSegment.oneMinusIndex * firstSegment.cotangentTheta / firstSegment.sineTheta -
                     secondSegment.oneMinusIndex * secondSegment.cotangentTheta / secondSegment.sineTheta) /
                    transverseCurvature);
  }

  if (rhoPhi.value == 0.f) {
    return false;
  }
  result = {full(phiTilde), full(thetaTilde), full(rhoPhi), rhoTheta};
  return true;
}

} // namespace triplet_detail
GPUhdi() bool makeTripletFitFactor(
  const std::array<GlobalMeasurement, 3>& measurements,
  TripletFitFactor& result) noexcept
{
  triplet_detail::TripletGeometry geometry;
  if (!triplet_detail::makeTripletGeometry(measurements, geometry)) {
    return false;
  }
  const float kappaReference = -geometry.phiTilde.value / geometry.rhoPhi.value;
  TripletFitFactor scratch{
    {static_cast<float>(geometry.thetaTilde.value), static_cast<float>(geometry.phiTilde.value)},
    {static_cast<float>(geometry.rhoTheta.value), static_cast<float>(geometry.rhoPhi.value)},
    {}};

  for (int hit = 0; hit < 3; ++hit) {
    for (int coordinate = 0; coordinate < 3; ++coordinate) {
      const int index = triplet_detail::coordinateOf(hit, coordinate);
      const float gradientTheta = geometry.thetaTilde.derivative[index] +
                                  kappaReference * geometry.rhoTheta.derivative[index];
      const float gradientPhi = geometry.phiTilde.derivative[index] +
                                kappaReference * geometry.rhoPhi.derivative[index];
      scratch.h[hit].theta[coordinate] = gradientTheta;
      scratch.h[hit].phi[coordinate] = gradientPhi;
    }
  }
  if (!scratch.isValid()) {
    return false;
  }
  result = scratch;
  return true;
}

// Shared CPU/GPU adjacent-triplet fit, including correlated hit uncertainties.
namespace adjacent_detail
{

constexpr std::size_t NAdjacentKinks = 4;

using KinkVector = std::array<float, NAdjacentKinks>;
using KinkCovariance = std::array<std::array<float, NAdjacentKinks>, NAdjacentKinks>;

GPUhdi() float covarianceContraction(const std::array<float, 3>& left,
                            const GlobalCovariance3F& covariance,
                            const std::array<float, 3>& right) noexcept
{
  return left[0] * (covariance.xx * right[0] + covariance.xy * right[1] + covariance.xz * right[2]) +
         left[1] * (covariance.xy * right[0] + covariance.yy * right[1] + covariance.yz * right[2]) +
         left[2] * (covariance.xz * right[0] + covariance.yz * right[1] + covariance.zz * right[2]);
}

GPUhdi() bool choleskyDecompose(const KinkCovariance& covariance,
                       KinkCovariance& lower) noexcept
{
  for (std::size_t row = 0; row < NAdjacentKinks; ++row) {
    for (std::size_t column = 0; column <= row; ++column) {
      float value = covariance[row][column];
      for (std::size_t k = 0; k < column; ++k) {
        value -= lower[row][k] * lower[column][k];
      }
      if (row == column) {
        if (value <= 0.f) {
          return false;
        }
        lower[row][column] = std::sqrt(value);
      } else {
        lower[row][column] = value / lower[column][column];
      }
    }
  }
  return true;
}

GPUhdi() bool choleskySolve(const KinkCovariance& lower, const KinkVector& right,
                   KinkVector& solution) noexcept
{
  KinkVector intermediate{};
  for (std::size_t row = 0; row < NAdjacentKinks; ++row) {
    float value = right[row];
    for (std::size_t column = 0; column < row; ++column) {
      value -= lower[row][column] * intermediate[column];
    }
    intermediate[row] = value / lower[row][row];
  }
  for (int row = static_cast<int>(NAdjacentKinks) - 1; row >= 0; --row) {
    float value = intermediate[row];
    for (std::size_t column = static_cast<std::size_t>(row) + 1;
         column < NAdjacentKinks; ++column) {
      value -= lower[column][row] * solution[column];
    }
    solution[row] = value / lower[row][row];
  }
  return true;
}

GPUhdi() float dotProduct(const KinkVector& left, const KinkVector& right) noexcept
{
  float result = 0.f;
  for (std::size_t i = 0; i < NAdjacentKinks; ++i) {
    result += left[i] * right[i];
  }
  return result;
}

GPUhdi() bool referenceSinTheta(const GlobalMeasurement& first,
                       const GlobalMeasurement& third,
                       float& sineTheta) noexcept
{
  const float dx = third.x - first.x;
  const float dy = third.y - first.y;
  const float dz = third.z - first.z;
  const float transverse = std::hypot(dx, dy);
  const float length = std::hypot(transverse, dz);
  sineTheta = transverse / length;
  return sineTheta > 0.f && sineTheta <= 1.f;
}

} // namespace adjacent_detail

GPUhdi() bool fitAdjacentTripletFactors(
  const TripletFitFactor& firstFactor,
  const TripletFitFactor& secondFactor,
  const std::array<GlobalMeasurement, 4>& measurements,
  const std::array<float, 2>& angularVariance,
  AdjacentTripletFitResult& result) noexcept
{
  using namespace adjacent_detail;
  std::array<float, 2> sineTheta{};
  if (!referenceSinTheta(measurements[0], measurements[2], sineTheta[0]) ||
      !referenceSinTheta(measurements[1], measurements[3], sineTheta[1])) {
    return false;
  }

  const KinkVector psi{
    firstFactor.psi.theta, firstFactor.psi.phi,
    secondFactor.psi.theta, secondFactor.psi.phi};
  const KinkVector rho{
    firstFactor.rho.theta, firstFactor.rho.phi,
    secondFactor.rho.theta, secondFactor.rho.phi};
  KinkCovariance covariance{};
  covariance[0][0] = angularVariance[0];
  covariance[1][1] = angularVariance[0] / (sineTheta[0] * sineTheta[0]);
  covariance[2][2] = angularVariance[1];
  covariance[3][3] = angularVariance[1] / (sineTheta[1] * sineTheta[1]);

  // Build H for four unique hits. The factors use slots (0,1,2) and (1,2,3),
  // so shared hits contribute to the cross-triplet covariance.
  std::array<std::array<std::array<float, 3>, NAdjacentKinks>, 4> gradients{};
  for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
    for (std::size_t hit = 0; hit < 3; ++hit) {
      gradients[hit][0][coordinate] = firstFactor.h[hit].theta[coordinate];
      gradients[hit][1][coordinate] = firstFactor.h[hit].phi[coordinate];
      gradients[hit + 1][2][coordinate] = secondFactor.h[hit].theta[coordinate];
      gradients[hit + 1][3][coordinate] = secondFactor.h[hit].phi[coordinate];
    }
  }
  for (std::size_t hit = 0; hit < measurements.size(); ++hit) {
    for (std::size_t row = 0; row < NAdjacentKinks; ++row) {
      for (std::size_t column = 0; column <= row; ++column) {
        const float contribution = covarianceContraction(
          gradients[hit][row], measurements[hit].covariance, gradients[hit][column]);
        covariance[row][column] += contribution;
        if (row != column) {
          covariance[column][row] += contribution;
        }
      }
    }
  }

  KinkCovariance lower{};
  KinkVector precisionPsi{};
  KinkVector precisionRho{};
  if (!choleskyDecompose(covariance, lower) ||
      !choleskySolve(lower, psi, precisionPsi) ||
      !choleskySolve(lower, rho, precisionRho)) {
    return false;
  }
  const float rhoPrecisionPsi = dotProduct(rho, precisionPsi);
  const float rhoPrecisionRho = dotProduct(rho, precisionRho);
  if (rhoPrecisionRho <= 0.f) {
    return false;
  }

  const float curvature = -rhoPrecisionPsi / rhoPrecisionRho;
  const float curvatureVariance = 1.f / rhoPrecisionRho;
  // Evaluate the minimized quadratic form from the residual. Subtracting
  // psi^T V^-1 psi and the fitted curvature term loses the small chi2 of
  // nearly exact helices, even with accurately constructed triplet factors.
  KinkVector residual{};
  KinkVector precisionResidual{};
  for (std::size_t i = 0; i < NAdjacentKinks; ++i) {
    residual[i] = std::fma(curvature, rho[i], psi[i]);
  }
  if (!choleskySolve(lower, residual, precisionResidual)) {
    return false;
  }
  const float chi2 = dotProduct(residual, precisionResidual);
  if (curvatureVariance <= 0.f || chi2 < 0.f) {
    return false;
  }

  result = {curvature, curvatureVariance, chi2};
  return true;
}


} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_TRIPLETFITTING_H_
