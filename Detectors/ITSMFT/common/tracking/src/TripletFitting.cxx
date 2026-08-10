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

#include "CommonConstants/MathConstants.h"

namespace o2::itsmft::tracking
{
namespace
{

constexpr std::size_t NMeasurementCoordinates = 9;
constexpr std::size_t CurvatureNoiseCoordinate = NMeasurementCoordinates;
constexpr std::size_t NCoordinates = NMeasurementCoordinates + 1;
constexpr std::size_t NAdjacentKinks = 4;

using KinkVector = std::array<double, NAdjacentKinks>;
using KinkCovariance = std::array<std::array<double, NAdjacentKinks>, NAdjacentKinks>;

struct Jet {
  double value{0.};
  std::array<double, NCoordinates> derivative{};

  static Jet variable(double value, std::size_t index) noexcept
  {
    Jet result{value};
    result.derivative[index] = 1.;
    return result;
  }
};

Jet operator+(const Jet& lhs, const Jet& rhs) noexcept
{
  Jet result{lhs.value + rhs.value};
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = lhs.derivative[i] + rhs.derivative[i];
  }
  return result;
}

Jet operator-(const Jet& lhs, const Jet& rhs) noexcept
{
  Jet result{lhs.value - rhs.value};
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = lhs.derivative[i] - rhs.derivative[i];
  }
  return result;
}

Jet operator-(const Jet& value) noexcept
{
  Jet result{-value.value};
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = -value.derivative[i];
  }
  return result;
}

Jet operator*(const Jet& lhs, const Jet& rhs) noexcept
{
  Jet result{lhs.value * rhs.value};
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = lhs.derivative[i] * rhs.value + lhs.value * rhs.derivative[i];
  }
  return result;
}

Jet operator/(const Jet& lhs, const Jet& rhs) noexcept
{
  const double inverse = 1. / rhs.value;
  Jet result{lhs.value * inverse};
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = (lhs.derivative[i] - result.value * rhs.derivative[i]) * inverse;
  }
  return result;
}

Jet squareRoot(const Jet& argument) noexcept
{
  const double root = std::sqrt(argument.value);
  Jet result{root};
  const double scale = 0.5 / root;
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = scale * argument.derivative[i];
  }
  return result;
}

Jet arcSine(const Jet& argument) noexcept
{
  Jet result{std::asin(argument.value)};
  const double scale = 1. / std::sqrt(1. - argument.value * argument.value);
  for (std::size_t i = 0; i < NCoordinates; ++i) {
    result.derivative[i] = scale * argument.derivative[i];
  }
  return result;
}

Jet arcTangent2(const Jet& y, const Jet& x) noexcept
{
  Jet result{std::atan2(y.value, x.value)};
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
  Jet bendingAngle;
  Jet transverseArcLength;
  Jet cotangentTheta;
  Jet sineTheta;
  Jet cosineTheta;
  Jet index;
};

bool makeSegmentGeometry(const Jet& transverseCurvature, const Jet& chordLength,
                         const Jet& deltaZ, SegmentGeometry& result) noexcept
{
  const Jet halfSine = Jet{0.5} * transverseCurvature * chordLength;
  if (!std::isfinite(halfSine.value) || std::abs(halfSine.value) >= 1.) {
    return false;
  }

  const Jet halfSine2 = halfSine * halfSine;
  const Jet halfSine4 = halfSine2 * halfSine2;
  Jet asinOverArgument;
  Jet angleCotangent;
  if (std::abs(halfSine.value) < 1.e-4) {
    asinOverArgument = Jet{1.} + halfSine2 * Jet{1. / 6.} + halfSine4 * Jet{3. / 40.};
    angleCotangent = Jet{1.} - halfSine2 * Jet{1. / 3.} - halfSine4 * Jet{2. / 15.};
  } else {
    const Jet halfAngle = arcSine(halfSine);
    asinOverArgument = halfAngle / halfSine;
    angleCotangent = halfAngle * squareRoot(Jet{1.} - halfSine2) / halfSine;
  }

  const Jet bendingAngle = Jet{2.} * arcSine(halfSine);
  const Jet transverseArcLength = chordLength * asinOverArgument;
  const Jet cotangentTheta = deltaZ / transverseArcLength;
  const Jet sineTheta = Jet{1.} / squareRoot(Jet{1.} + cotangentTheta * cotangentTheta);
  const Jet cosineTheta = cotangentTheta * sineTheta;
  const Jet index = Jet{1.} /
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
  Jet phiTilde;
  Jet thetaTilde;
  Jet rhoPhi;
  Jet rhoTheta;
  Jet chordAzimuth01;
  Jet chordAzimuth12;
  Jet chordLength01;
  Jet chordLength12;
  Jet deltaZ01;
  Jet deltaZ12;
  Jet referenceSinThetaJet;
  double referenceTransverseCurvature{0.};
  double referenceSinTheta{0.};
};

bool makeTripletGeometry(const std::array<TripletFitObservation, 3>& observations,
                         TripletGeometry& result) noexcept
{
  std::array<std::array<Jet, 3>, 3> point{};
  for (std::size_t hit = 0; hit < observations.size(); ++hit) {
    for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
      const std::size_t index = 3 * hit + coordinate;
      point[hit][coordinate] = Jet::variable(observations[hit].position[coordinate], index);
    }
  }

  const Jet dx01 = point[1][0] - point[0][0];
  const Jet dy01 = point[1][1] - point[0][1];
  const Jet dz01 = point[1][2] - point[0][2];
  const Jet dx12 = point[2][0] - point[1][0];
  const Jet dy12 = point[2][1] - point[1][1];
  const Jet dz12 = point[2][2] - point[1][2];
  const Jet dx02 = point[2][0] - point[0][0];
  const Jet dy02 = point[2][1] - point[0][1];
  const Jet dz02 = point[2][2] - point[0][2];

  const Jet length01 = squareRoot(dx01 * dx01 + dy01 * dy01);
  const Jet length12 = squareRoot(dx12 * dx12 + dy12 * dy12);
  const Jet length02 = squareRoot(dx02 * dx02 + dy02 * dy02);
  if (!std::isfinite(length01.value) || !std::isfinite(length12.value) ||
      !std::isfinite(length02.value) || length01.value <= 0. ||
      length12.value <= 0. || length02.value <= 0.) {
    return false;
  }

  const Jet cross = dx01 * dy12 - dy01 * dx12;
  const Jet chordAzimuth01 = arcTangent2(dy01, dx01);
  const Jet chordAzimuth12 = arcTangent2(dy12, dx12);
  const Jet transverseCurvature = Jet{2.} * cross / (length01 * length12 * length02);
  SegmentGeometry firstSegment;
  SegmentGeometry secondSegment;
  if (!makeSegmentGeometry(transverseCurvature, length01, dz01, firstSegment) ||
      !makeSegmentGeometry(transverseCurvature, length12, dz12, secondSegment)) {
    return false;
  }

  const Jet theta01 = arcTangent2(firstSegment.transverseArcLength, dz01);
  const Jet theta12 = arcTangent2(secondSegment.transverseArcLength, dz12);
  const Jet phiTilde = Jet{0.5} *
                       (firstSegment.bendingAngle * firstSegment.index +
                        secondSegment.bendingAngle * secondSegment.index);
  const Jet thetaTilde = theta12 - theta01 +
                         (Jet{1.} - secondSegment.index) * secondSegment.cotangentTheta -
                         (Jet{1.} - firstSegment.index) * firstSegment.cotangentTheta;
  const Jet rhoPhi = Jet{-0.5} *
                     (firstSegment.transverseArcLength * firstSegment.index / firstSegment.sineTheta +
                      secondSegment.transverseArcLength * secondSegment.index / secondSegment.sineTheta);

  Jet rhoTheta;
  const double maximumHalfSine = 0.5 * std::abs(transverseCurvature.value) *
                                 std::max(length01.value, length12.value);
  if (maximumHalfSine < 1.e-4) {
    rhoTheta = transverseCurvature *
               (length12 * length12 * secondSegment.cosineTheta -
                length01 * length01 * firstSegment.cosineTheta) /
               Jet{12.};
  } else {
    rhoTheta = ((Jet{1.} - firstSegment.index) * firstSegment.cotangentTheta / firstSegment.sineTheta -
                (Jet{1.} - secondSegment.index) * secondSegment.cotangentTheta / secondSegment.sineTheta) /
               transverseCurvature;
  }

  const Jet endpointLength = squareRoot(length02 * length02 + dz02 * dz02);
  const Jet referenceSinTheta = length02 / endpointLength;
  const double sineTheta = referenceSinTheta.value;
  if (!std::isfinite(phiTilde.value) || !std::isfinite(thetaTilde.value) ||
      !std::isfinite(rhoPhi.value) || rhoPhi.value == 0. ||
      !std::isfinite(rhoTheta.value) || !std::isfinite(sineTheta) || sineTheta <= 0.) {
    return false;
  }
  result = {phiTilde, thetaTilde, rhoPhi, rhoTheta,
            chordAzimuth01, chordAzimuth12, length01, length12, dz01, dz12,
            referenceSinTheta,
            transverseCurvature.value, sineTheta};
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

bool choleskySolve(const KinkCovariance& covariance, const KinkVector& right,
                   KinkVector& solution) noexcept
{
  KinkCovariance lower{};
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

double jetCovariance(const Jet& left, const Jet& right,
                     const std::array<TripletFitObservation, 3>& observations,
                     double independentCurvatureVariance) noexcept
{
  double covariance = left.derivative[CurvatureNoiseCoordinate] *
                      right.derivative[CurvatureNoiseCoordinate] *
                      independentCurvatureVariance;
  for (std::size_t hit = 0; hit < observations.size(); ++hit) {
    std::array<double, 3> leftGradient{};
    std::array<double, 3> rightGradient{};
    for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
      const auto index = 3 * hit + coordinate;
      leftGradient[coordinate] = left.derivative[index];
      rightGradient[coordinate] = right.derivative[index];
    }
    covariance += covarianceContraction(leftGradient, observations[hit].covariance, rightGradient);
  }
  return covariance;
}

bool makeSegmentEstimate(const Jet& qOverPt, const Jet& phi, const Jet& tanLambda,
                         const std::array<TripletFitObservation, 3>& observations,
                         double independentCurvatureVariance,
                         TripletSegmentEstimate& estimate) noexcept
{
  const TripletParameterCovariance covariance{
    jetCovariance(qOverPt, qOverPt, observations, independentCurvatureVariance),
    jetCovariance(qOverPt, phi, observations, independentCurvatureVariance),
    jetCovariance(qOverPt, tanLambda, observations, independentCurvatureVariance),
    jetCovariance(phi, phi, observations, independentCurvatureVariance),
    jetCovariance(phi, tanLambda, observations, independentCurvatureVariance),
    jetCovariance(tanLambda, tanLambda, observations, independentCurvatureVariance)};
  const std::array<double, 9> values{
    qOverPt.value, phi.value, tanLambda.value,
    covariance.qOverPtQOverPt, covariance.qOverPtPhi,
    covariance.qOverPtTanLambda, covariance.phiPhi,
    covariance.phiTanLambda, covariance.tanLambdaTanLambda};
  if (!std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); }) ||
      covariance.qOverPtQOverPt < 0. || covariance.phiPhi < 0. ||
      covariance.tanLambdaTanLambda < 0.) {
    return false;
  }
  estimate = {qOverPt.value, phi.value, tanLambda.value, covariance};
  return true;
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

bool fitLocalTripletUniformSolenoid(
  const std::array<TripletFitObservation, 3>& observations,
  const TripletFitProcessNoise& processNoise,
  double bz,
  LocalTripletFitResult& result) noexcept
{
  if (!std::all_of(observations.begin(), observations.end(), observationIsValid) ||
      !std::isfinite(processNoise.angularVariance) || processNoise.angularVariance < 0. ||
      !std::isfinite(bz) || bz == 0.) {
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
  TripletFitFactor factor{
    {static_cast<float>(geometry.thetaTilde.value), static_cast<float>(geometry.phiTilde.value)},
    {static_cast<float>(geometry.rhoTheta.value), static_cast<float>(geometry.rhoPhi.value)},
    {}};

  double gammaThetaTheta = processNoise.angularVariance;
  double gammaThetaPhi = 0.;
  double gammaPhiPhi = processNoise.angularVariance /
                       (geometry.referenceSinTheta * geometry.referenceSinTheta);
  for (std::size_t hit = 0; hit < observations.size(); ++hit) {
    std::array<double, 3> gradientTheta{};
    std::array<double, 3> gradientPhi{};
    for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
      const std::size_t index = 3 * hit + coordinate;
      gradientTheta[coordinate] = geometry.thetaTilde.derivative[index] +
                                  kappaReference * geometry.rhoTheta.derivative[index];
      gradientPhi[coordinate] = geometry.phiTilde.derivative[index] +
                                kappaReference * geometry.rhoPhi.derivative[index];
      if (!finiteFloat(gradientTheta[coordinate]) || !finiteFloat(gradientPhi[coordinate])) {
        return false;
      }
      factor.h[hit].theta[coordinate] = static_cast<float>(gradientTheta[coordinate]);
      factor.h[hit].phi[coordinate] = static_cast<float>(gradientPhi[coordinate]);
    }
    gammaThetaTheta += covarianceContraction(gradientTheta, observations[hit].covariance, gradientTheta);
    gammaThetaPhi += covarianceContraction(gradientTheta, observations[hit].covariance, gradientPhi);
    gammaPhiPhi += covarianceContraction(gradientPhi, observations[hit].covariance, gradientPhi);
  }
  if (!factor.isValid()) {
    return false;
  }

  const double covarianceDeterminant = gammaThetaTheta * gammaPhiPhi - gammaThetaPhi * gammaThetaPhi;
  const double rhoTheta = geometry.rhoTheta.value;
  const double rhoPhi = geometry.rhoPhi.value;
  const double denominator = rhoTheta * rhoTheta * gammaPhiPhi +
                             rhoPhi * rhoPhi * gammaThetaTheta -
                             2. * rhoTheta * rhoPhi * gammaThetaPhi;
  if (!std::isfinite(gammaThetaTheta) || !std::isfinite(gammaThetaPhi) ||
      !std::isfinite(gammaPhiPhi) || gammaThetaTheta < 0. || gammaPhiPhi < 0. ||
      !std::isfinite(covarianceDeterminant) || covarianceDeterminant <= 0. ||
      !std::isfinite(denominator) || denominator <= 0.) {
    return false;
  }

  const double thetaTilde = geometry.thetaTilde.value;
  const double phiTilde = geometry.phiTilde.value;
  const double curvature =
    -(thetaTilde * rhoTheta * gammaPhiPhi +
      phiTilde * rhoPhi * gammaThetaTheta -
      gammaThetaPhi * (phiTilde * rhoTheta + thetaTilde * rhoPhi)) /
    denominator;
  const double curvatureVariance = covarianceDeterminant / denominator;
  const double qualityResidual = thetaTilde * rhoPhi - phiTilde * rhoTheta;
  const double chi2 = qualityResidual * qualityResidual / denominator;
  if (!std::isfinite(curvature) || !std::isfinite(curvatureVariance) || curvatureVariance <= 0. ||
      !std::isfinite(chi2) || chi2 < 0.) {
    return false;
  }

  const Jet curvatureJet =
    -(geometry.thetaTilde * geometry.rhoTheta * Jet{gammaPhiPhi} +
      geometry.phiTilde * geometry.rhoPhi * Jet{gammaThetaTheta} -
      Jet{gammaThetaPhi} * (geometry.phiTilde * geometry.rhoTheta +
                            geometry.thetaTilde * geometry.rhoPhi)) /
    (geometry.rhoTheta * geometry.rhoTheta * Jet{gammaPhiPhi} +
     geometry.rhoPhi * geometry.rhoPhi * Jet{gammaThetaTheta} -
     Jet{2. * gammaThetaPhi} * geometry.rhoTheta * geometry.rhoPhi);
  if (!std::isfinite(curvatureJet.value)) {
    return false;
  }
  Jet curvatureWithProcess = curvatureJet;
  curvatureWithProcess.derivative[CurvatureNoiseCoordinate] = 1.;
  const double measurementCurvatureVariance = jetCovariance(curvatureJet, curvatureJet, observations, 0.);
  if (!std::isfinite(measurementCurvatureVariance) || measurementCurvatureVariance < 0.) {
    return false;
  }
  // The local fit variance includes the supplied MS model. Its part not
  // generated by the hit-coordinate Jacobian is represented as one latent
  // curvature uncertainty and propagated into q/pT and both segment slopes.
  const double independentCurvatureVariance =
    std::max(0., curvatureVariance - measurementCurvatureVariance);

  const Jet magneticScale{static_cast<double>(o2::constants::math::B2C) * bz};
  const Jet qOverPt = curvatureWithProcess /
                      (magneticScale * geometry.referenceSinThetaJet);
  SegmentGeometry firstFittedSegment;
  SegmentGeometry secondFittedSegment;
  if (!makeSegmentGeometry(curvatureWithProcess, geometry.chordLength01,
                           geometry.deltaZ01, firstFittedSegment) ||
      !makeSegmentGeometry(curvatureWithProcess, geometry.chordLength12,
                           geometry.deltaZ12, secondFittedSegment)) {
    return false;
  }

  TripletSegmentEstimate inner;
  TripletSegmentEstimate outer;
  if (!makeSegmentEstimate(qOverPt, geometry.chordAzimuth01,
                           firstFittedSegment.cotangentTheta, observations,
                           independentCurvatureVariance, inner) ||
      !makeSegmentEstimate(qOverPt, geometry.chordAzimuth12,
                           secondFittedSegment.cotangentTheta, observations,
                           independentCurvatureVariance, outer)) {
    return false;
  }

  result = {curvature,
            curvatureVariance,
            chi2,
            geometry.referenceTransverseCurvature,
            geometry.referenceSinTheta,
            phiTilde,
            thetaTilde,
            rhoPhi,
            rhoTheta,
            gammaThetaTheta,
            gammaThetaPhi,
            gammaPhiPhi,
            inner,
            outer,
            factor};
  return true;
}

bool fitAdjacentTripletFactors(
  const std::array<TripletFitFactor, 2>& factors,
  const std::array<TripletFitObservation, 4>& observations,
  const std::array<double, 2>& angularVariance,
  AdjacentTripletFitResult& result) noexcept
{
  if (!factors[0].isValid() || !factors[1].isValid() ||
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
    factors[0].psi.theta, factors[0].psi.phi,
    factors[1].psi.theta, factors[1].psi.phi};
  const KinkVector rho{
    factors[0].rho.theta, factors[0].rho.phi,
    factors[1].rho.theta, factors[1].rho.phi};
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
      gradients[hit][0][coordinate] = factors[0].h[hit].theta[coordinate];
      gradients[hit][1][coordinate] = factors[0].h[hit].phi[coordinate];
      gradients[hit + 1][2][coordinate] = factors[1].h[hit].theta[coordinate];
      gradients[hit + 1][3][coordinate] = factors[1].h[hit].phi[coordinate];
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

  KinkVector precisionPsi{};
  KinkVector precisionRho{};
  if (!choleskySolve(covariance, psi, precisionPsi) ||
      !choleskySolve(covariance, rho, precisionRho)) {
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

  result = {curvature, curvatureVariance, chi2, covariance};
  return true;
}

double fittedTripletTransverseMomentum(const LocalTripletFitResult& result,
                                       double bz) noexcept
{
  return fittedTripletTransverseMomentum(result, bz, 1);
}

double fittedTripletTransverseMomentum(const TripletSegmentEstimate& estimate,
                                       uint8_t absCharge) noexcept
{
  if (!std::isfinite(estimate.qOverPt) || absCharge == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (estimate.qOverPt == 0.) {
    return std::numeric_limits<double>::infinity();
  }
  return static_cast<double>(absCharge) / std::abs(estimate.qOverPt);
}

double fittedTripletTransverseMomentumVariance(const TripletSegmentEstimate& estimate,
                                               uint8_t absCharge) noexcept
{
  if (!std::isfinite(estimate.qOverPt) || estimate.qOverPt == 0. ||
      !std::isfinite(estimate.covariance.qOverPtQOverPt) ||
      estimate.covariance.qOverPtQOverPt < 0. || absCharge == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double qOverPt2 = estimate.qOverPt * estimate.qOverPt;
  return static_cast<double>(absCharge) * static_cast<double>(absCharge) *
         estimate.covariance.qOverPtQOverPt / (qOverPt2 * qOverPt2);
}

double fittedTripletTransverseMomentum(const LocalTripletFitResult& result,
                                       double bz, uint8_t absCharge) noexcept
{
  if (!std::isfinite(result.curvature) || !std::isfinite(result.referenceSinTheta) ||
      result.referenceSinTheta <= 0. || result.referenceSinTheta > 1. ||
      !std::isfinite(bz) || bz == 0. || absCharge == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (result.curvature == 0.) {
    return std::numeric_limits<double>::infinity();
  }
  return static_cast<double>(absCharge) *
         std::abs(static_cast<double>(o2::constants::math::B2C) * bz *
                  result.referenceSinTheta / result.curvature);
}

} // namespace o2::itsmft::tracking
