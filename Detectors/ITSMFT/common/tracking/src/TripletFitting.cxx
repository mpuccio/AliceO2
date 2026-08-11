// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/TripletFitting.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace o2::itsmft::tracking
{
namespace
{

constexpr std::size_t NMeasurementCoordinates = 9;
constexpr std::size_t NCoordinates = NMeasurementCoordinates;
constexpr std::size_t NAdjacentKinks = 4;

using KinkVector = std::array<double, NAdjacentKinks>;
using KinkCovariance = std::array<std::array<double, NAdjacentKinks>, NAdjacentKinks>;

struct DualNumber {
  double value{0.};
  std::array<double, NCoordinates> derivative{};

  static DualNumber variable(double val, std::size_t index) noexcept
  {
    DualNumber result{val};
    result.derivative[index] = 1.;
    return result;
  }
};

DualNumber operator+(const DualNumber& lhs, const DualNumber& rhs) noexcept
{
  DualNumber result{lhs.value + rhs.value};
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = lhs.derivative[i] + rhs.derivative[i];
  }
  return result;
}

DualNumber operator-(const DualNumber& lhs, const DualNumber& rhs) noexcept
{
  DualNumber result{lhs.value - rhs.value};
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = lhs.derivative[i] - rhs.derivative[i];
  }
  return result;
}

DualNumber operator-(const DualNumber& value) noexcept
{
  DualNumber result{-value.value};
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = -value.derivative[i];
  }
  return result;
}

DualNumber operator*(const DualNumber& lhs, const DualNumber& rhs) noexcept
{
  DualNumber result{lhs.value * rhs.value};
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = lhs.derivative[i] * rhs.value + lhs.value * rhs.derivative[i];
  }
  return result;
}

DualNumber operator/(const DualNumber& lhs, const DualNumber& rhs) noexcept
{
  const double inverse = 1. / rhs.value;
  DualNumber result{lhs.value * inverse};
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = (lhs.derivative[i] - result.value * rhs.derivative[i]) * inverse;
  }
  return result;
}

DualNumber squareRoot(const DualNumber& argument) noexcept
{
  const double root = std::sqrt(argument.value);
  DualNumber result{root};
  const double scale = 0.5 / root;
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = scale * argument.derivative[i];
  }
  return result;
}

DualNumber arcSine(const DualNumber& argument) noexcept
{
  DualNumber result{std::asin(argument.value)};
  const double scale = 1. / std::sqrt(1. - argument.value * argument.value);
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = scale * argument.derivative[i];
  }
  return result;
}

DualNumber arcTangent2(const DualNumber& y, const DualNumber& x) noexcept
{
  DualNumber result{std::atan2(y.value, x.value)};
  const double denominator = x.value * x.value + y.value * y.value;
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = (x.value * y.derivative[i] - y.value * x.derivative[i]) / denominator;
  }
  return result;
}

bool covarianceIsPositiveSemidefinite(const SymmetricCovariance3D& covariance) noexcept
{
  const std::array<double, 6> elements{covariance.xx, covariance.xy, covariance.xz,
                                       covariance.yy, covariance.yz, covariance.zz};
  if (!std::all_of(elements.begin(), elements.end(), [](double value) { return std::isfinite(value); }) ||
      covariance.xx < 0. || covariance.yy < 0. || covariance.zz < 0.) {
    return false;
  }

  const std::array<double, 3> minors{
    covariance.xx * covariance.yy - covariance.xy * covariance.xy,
    covariance.xx * covariance.zz - covariance.xz * covariance.xz,
    covariance.yy * covariance.zz - covariance.yz * covariance.yz};
  for (const double minor : minors) {
    const double scale = std::max({std::abs(covariance.xx * covariance.yy),
                                   std::abs(covariance.xx * covariance.zz),
                                   std::abs(covariance.yy * covariance.zz),
                                   covariance.xy * covariance.xy,
                                   covariance.xz * covariance.xz,
                                   covariance.yz * covariance.yz});
    const double tolerance = 16. * std::numeric_limits<float>::epsilon() * scale;
    if (minor < -tolerance) {
      return false;
    }
  }

  const double positiveTerm = covariance.xx * covariance.yy * covariance.zz +
                              2. * covariance.xy * covariance.xz * covariance.yz;
  const double negativeTerm = covariance.xx * covariance.yz * covariance.yz +
                              covariance.yy * covariance.xz * covariance.xz +
                              covariance.zz * covariance.xy * covariance.xy;
  const double determinant = positiveTerm - negativeTerm;
  const double tolerance = 32. * std::numeric_limits<float>::epsilon() *
                           std::max(std::abs(positiveTerm), std::abs(negativeTerm));
  return determinant >= -tolerance;
}

bool observationIsValid(const TripletFitObservation& observation) noexcept
{
  return std::all_of(observation.position.begin(), observation.position.end(),
                     [](double value) { return std::isfinite(value); }) &&
         covarianceIsPositiveSemidefinite(observation.covariance);
}

struct SegmentGeometry {
  DualNumber bendingAngle;
  DualNumber transverseArcLength;
  DualNumber cotangentTheta;
  DualNumber sineTheta;
  DualNumber cosineTheta;
  DualNumber index;
};

bool makeSegmentGeometry(const DualNumber& transverseCurvature, const DualNumber& chordLength,
                         const DualNumber& deltaZ, SegmentGeometry& result) noexcept
{
  const DualNumber halfSine = DualNumber{0.5} * transverseCurvature * chordLength;
  if (!std::isfinite(halfSine.value) || std::abs(halfSine.value) >= 1.) {
    return false;
  }

  const DualNumber halfSine2 = halfSine * halfSine;
  const DualNumber halfSine4 = halfSine2 * halfSine2;
  DualNumber asinOverArgument;
  DualNumber angleCotangent;
  if (std::abs(halfSine.value) < 1.e-4) {
    asinOverArgument = DualNumber{1.} + halfSine2 * DualNumber{1. / 6.} + halfSine4 * DualNumber{3. / 40.};
    angleCotangent = DualNumber{1.} - halfSine2 * DualNumber{1. / 3.} - halfSine4 * DualNumber{2. / 15.};
  } else {
    const DualNumber halfAngle = arcSine(halfSine);
    asinOverArgument = halfAngle / halfSine;
    angleCotangent = halfAngle * squareRoot(DualNumber{1.} - halfSine2) / halfSine;
  }

  const DualNumber bendingAngle = DualNumber{2.} * arcSine(halfSine);
  const DualNumber transverseArcLength = chordLength * asinOverArgument;
  const DualNumber cotangentTheta = deltaZ / transverseArcLength;
  const DualNumber sineTheta = DualNumber{1.} / squareRoot(DualNumber{1.} + cotangentTheta * cotangentTheta);
  const DualNumber cosineTheta = cotangentTheta * sineTheta;
  const DualNumber index = DualNumber{1.} /
                    (angleCotangent * sineTheta * sineTheta + cosineTheta * cosineTheta);
  if (!std::isfinite(transverseArcLength.value) || transverseArcLength.value <= 0. ||
      !std::isfinite(sineTheta.value) || sineTheta.value <= 0. ||
      !std::isfinite(index.value) || index.value <= 0.) {
    return false;
  }
  result = {bendingAngle, transverseArcLength, cotangentTheta, sineTheta, cosineTheta, index};
  return true;
}

struct TripletGeometry {
  DualNumber phiTilde;
  DualNumber thetaTilde;
  DualNumber rhoPhi;
  DualNumber rhoTheta;
};

bool makeTripletGeometry(const std::array<TripletFitObservation, 3>& observations,
                         TripletGeometry& result) noexcept
{
  std::array<std::array<DualNumber, 3>, 3> point{};
  for (std::size_t hit = 0; hit < observations.size(); ++hit) {
    for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
      const std::size_t index = 3 * hit + coordinate;
      point[hit][coordinate] = DualNumber::variable(observations[hit].position[coordinate], index);
    }
  }

  const DualNumber dx01 = point[1][0] - point[0][0];
  const DualNumber dy01 = point[1][1] - point[0][1];
  const DualNumber dz01 = point[1][2] - point[0][2];
  const DualNumber dx12 = point[2][0] - point[1][0];
  const DualNumber dy12 = point[2][1] - point[1][1];
  const DualNumber dz12 = point[2][2] - point[1][2];
  const DualNumber dx02 = point[2][0] - point[0][0];
  const DualNumber dy02 = point[2][1] - point[0][1];
  const DualNumber length01 = squareRoot(dx01 * dx01 + dy01 * dy01);
  const DualNumber length12 = squareRoot(dx12 * dx12 + dy12 * dy12);
  const DualNumber length02 = squareRoot(dx02 * dx02 + dy02 * dy02);
  if (!std::isfinite(length01.value) || !std::isfinite(length12.value) ||
      !std::isfinite(length02.value) || length01.value <= 0. ||
      length12.value <= 0. || length02.value <= 0.) {
    return false;
  }

  const DualNumber cross = dx01 * dy12 - dy01 * dx12;
  const DualNumber transverseCurvature = DualNumber{2.} * cross / (length01 * length12 * length02);
  SegmentGeometry firstSegment;
  SegmentGeometry secondSegment;
  if (!makeSegmentGeometry(transverseCurvature, length01, dz01, firstSegment) ||
      !makeSegmentGeometry(transverseCurvature, length12, dz12, secondSegment)) {
    return false;
  }

  const DualNumber theta01 = arcTangent2(firstSegment.transverseArcLength, dz01);
  const DualNumber theta12 = arcTangent2(secondSegment.transverseArcLength, dz12);
  const DualNumber phiTilde = DualNumber{0.5} *
                       (firstSegment.bendingAngle * firstSegment.index +
                        secondSegment.bendingAngle * secondSegment.index);
  const DualNumber thetaTilde = theta12 - theta01 +
                         (DualNumber{1.} - secondSegment.index) * secondSegment.cotangentTheta -
                         (DualNumber{1.} - firstSegment.index) * firstSegment.cotangentTheta;
  const DualNumber rhoPhi = DualNumber{-0.5} *
                     (firstSegment.transverseArcLength * firstSegment.index / firstSegment.sineTheta +
                      secondSegment.transverseArcLength * secondSegment.index / secondSegment.sineTheta);

  DualNumber rhoTheta;
  const double maximumHalfSine = 0.5 * std::abs(transverseCurvature.value) *
                                 std::max(length01.value, length12.value);
  if (maximumHalfSine < 1.e-4) {
    rhoTheta = transverseCurvature *
               (length12 * length12 * secondSegment.cosineTheta -
                length01 * length01 * firstSegment.cosineTheta) /
               DualNumber{12.};
  } else {
    rhoTheta = ((DualNumber{1.} - firstSegment.index) * firstSegment.cotangentTheta / firstSegment.sineTheta -
                (DualNumber{1.} - secondSegment.index) * secondSegment.cotangentTheta / secondSegment.sineTheta) /
               transverseCurvature;
  }

  if (!std::isfinite(phiTilde.value) || !std::isfinite(thetaTilde.value) ||
      !std::isfinite(rhoPhi.value) || rhoPhi.value == 0. ||
      !std::isfinite(rhoTheta.value)) {
    return false;
  }
  result = {phiTilde, thetaTilde, rhoPhi, rhoTheta};
  return true;
}

double covarianceContraction(const std::array<double, 3>& left,
                             const SymmetricCovariance3D& covariance,
                             const std::array<double, 3>& right) noexcept
{
  return left[0] * (covariance.xx * right[0] + covariance.xy * right[1] + covariance.xz * right[2]) +
         left[1] * (covariance.xy * right[0] + covariance.yy * right[1] + covariance.yz * right[2]) +
         left[2] * (covariance.xz * right[0] + covariance.yz * right[1] + covariance.zz * right[2]);
}

bool choleskyDecompose(const KinkCovariance& covariance,
                       KinkCovariance& lower) noexcept
{
  for (std::size_t row = 0; row < NAdjacentKinks; ++row) {
    for (std::size_t column = 0; column <= row; ++column) {
      double value = covariance[row][column];
      for (std::size_t k = 0; k < column; ++k) {
        value -= lower[row][k] * lower[column][k];
      }
      if (!std::isfinite(value)) {
        return false;
      }
      if (row == column) {
        if (value <= 0.) {
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

bool choleskySolve(const KinkCovariance& lower, const KinkVector& right,
                   KinkVector& solution) noexcept
{
  KinkVector intermediate{};
  for (std::size_t row = 0; row < NAdjacentKinks; ++row) {
    double value = right[row];
    for (std::size_t column = 0; column < row; ++column) {
      value -= lower[row][column] * intermediate[column];
    }
    intermediate[row] = value / lower[row][row];
    if (!std::isfinite(intermediate[row])) {
      return false;
    }
  }
  for (int row = static_cast<int>(NAdjacentKinks) - 1; row >= 0; --row) {
    double value = intermediate[row];
    for (std::size_t column = static_cast<std::size_t>(row) + 1;
         column < NAdjacentKinks; ++column) {
      value -= lower[column][row] * solution[column];
    }
    solution[row] = value / lower[row][row];
    if (!std::isfinite(solution[row])) {
      return false;
    }
  }
  return true;
}

double dotProduct(const KinkVector& left, const KinkVector& right) noexcept
{
  double result = 0.;
  for (std::size_t i = 0; i < NAdjacentKinks; ++i) {
    result += left[i] * right[i];
  }
  return result;
}

bool referenceSinTheta(const TripletFitObservation& first,
                       const TripletFitObservation& third,
                       double& sineTheta) noexcept
{
  const double dx = third.position[0] - first.position[0];
  const double dy = third.position[1] - first.position[1];
  const double dz = third.position[2] - first.position[2];
  const double transverse = std::hypot(dx, dy);
  const double length = std::hypot(transverse, dz);
  if (!std::isfinite(transverse) || !std::isfinite(length) ||
      transverse <= 0. || length <= 0.) {
    return false;
  }
  sineTheta = transverse / length;
  return std::isfinite(sineTheta) && sineTheta > 0. && sineTheta <= 1.;
}

} // namespace

bool makeTripletFitObservation(const GlobalMeasurement& measurement,
                               TripletFitObservation& observation) noexcept
{
  const TripletFitObservation scratch{
    {measurement.position.x, measurement.position.y, measurement.position.z},
    {measurement.covariance.xx, measurement.covariance.xy, measurement.covariance.xz,
     measurement.covariance.yy, measurement.covariance.yz, measurement.covariance.zz}};
  if (!observationIsValid(scratch)) {
    return false;
  }
  observation = scratch;
  return true;
}

bool makeTripletFitFactor(
  const std::array<TripletFitObservation, 3>& observations,
  TripletFitFactor& result) noexcept
{
  if (!std::all_of(observations.begin(), observations.end(), observationIsValid)) {
    return false;
  }

  TripletGeometry geometry;
  if (!makeTripletGeometry(observations, geometry)) {
    return false;
  }
  const double kappaReference = -geometry.phiTilde.value / geometry.rhoPhi.value;
  if (!std::isfinite(kappaReference)) {
    return false;
  }

  const auto finiteFloat = [](double value) {
    return std::isfinite(value) &&
           std::abs(value) <= static_cast<double>(std::numeric_limits<float>::max());
  };
  const std::array<double, 4> factorParameters{
    geometry.thetaTilde.value, geometry.phiTilde.value,
    geometry.rhoTheta.value, geometry.rhoPhi.value};
  if (!std::all_of(factorParameters.begin(), factorParameters.end(), finiteFloat)) {
    return false;
  }
  TripletFitFactor scratch{
    {static_cast<float>(geometry.thetaTilde.value), static_cast<float>(geometry.phiTilde.value)},
    {static_cast<float>(geometry.rhoTheta.value), static_cast<float>(geometry.rhoPhi.value)},
    {}};

  for (std::size_t hit = 0; hit < observations.size(); ++hit) {
    for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
      const std::size_t index = 3 * hit + coordinate;
      const double gradientTheta = geometry.thetaTilde.derivative[index] +
                                   kappaReference * geometry.rhoTheta.derivative[index];
      const double gradientPhi = geometry.phiTilde.derivative[index] +
                                 kappaReference * geometry.rhoPhi.derivative[index];
      if (!finiteFloat(gradientTheta) || !finiteFloat(gradientPhi)) {
        return false;
      }
      scratch.h[hit].theta[coordinate] = static_cast<float>(gradientTheta);
      scratch.h[hit].phi[coordinate] = static_cast<float>(gradientPhi);
    }
  }
  if (!scratch.isValid()) {
    return false;
  }
  result = scratch;
  return true;
}

bool fitAdjacentTripletFactors(
  const TripletFitFactor& firstFactor,
  const TripletFitFactor& secondFactor,
  const std::array<TripletFitObservation, 4>& observations,
  const std::array<double, 2>& angularVariance,
  AdjacentTripletFitResult& result) noexcept
{
  if (!firstFactor.isValid() || !secondFactor.isValid() ||
      !std::all_of(observations.begin(), observations.end(), observationIsValid) ||
      !std::all_of(angularVariance.begin(), angularVariance.end(),
                   [](double value) { return std::isfinite(value) && value >= 0.; })) {
    return false;
  }

  std::array<double, 2> sineTheta{};
  if (!referenceSinTheta(observations[0], observations[2], sineTheta[0]) ||
      !referenceSinTheta(observations[1], observations[3], sineTheta[1])) {
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

  // H for four unique hits. Adjacent factors use hit slots (0,1,2) and
  // (1,2,3), so the shared-hit rows generate the cross-triplet covariance.
  std::array<std::array<std::array<double, 3>, NAdjacentKinks>, 4> gradients{};
  for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
    for (std::size_t hit = 0; hit < 3; ++hit) {
      gradients[hit][0][coordinate] = firstFactor.h[hit].theta[coordinate];
      gradients[hit][1][coordinate] = firstFactor.h[hit].phi[coordinate];
      gradients[hit + 1][2][coordinate] = secondFactor.h[hit].theta[coordinate];
      gradients[hit + 1][3][coordinate] = secondFactor.h[hit].phi[coordinate];
    }
  }
  for (std::size_t hit = 0; hit < observations.size(); ++hit) {
    for (std::size_t row = 0; row < NAdjacentKinks; ++row) {
      for (std::size_t column = 0; column <= row; ++column) {
        const double contribution = covarianceContraction(
          gradients[hit][row], observations[hit].covariance, gradients[hit][column]);
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
  const double rhoPrecisionPsi = dotProduct(rho, precisionPsi);
  const double rhoPrecisionRho = dotProduct(rho, precisionRho);
  const double psiPrecisionPsi = dotProduct(psi, precisionPsi);
  if (!std::isfinite(rhoPrecisionPsi) || !std::isfinite(rhoPrecisionRho) ||
      rhoPrecisionRho <= 0. || !std::isfinite(psiPrecisionPsi)) {
    return false;
  }

  const double curvature = -rhoPrecisionPsi / rhoPrecisionRho;
  const double curvatureVariance = 1. / rhoPrecisionRho;
  const double removedCurvatureTerm = rhoPrecisionPsi * rhoPrecisionPsi / rhoPrecisionRho;
  double chi2 = psiPrecisionPsi - removedCurvatureTerm;
  const double chi2Tolerance = 128. * std::numeric_limits<double>::epsilon() *
                               std::max(std::abs(psiPrecisionPsi), std::abs(removedCurvatureTerm));
  if (chi2 < 0. && chi2 >= -chi2Tolerance) {
    chi2 = 0.;
  }
  if (!std::isfinite(curvature) || !std::isfinite(curvatureVariance) ||
      curvatureVariance <= 0. || !std::isfinite(chi2) || chi2 < 0.) {
    return false;
  }

  result = {curvature, curvatureVariance, chi2};
  return true;
}

} // namespace o2::itsmft::tracking
