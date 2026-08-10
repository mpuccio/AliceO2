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

constexpr std::size_t NCoordinates = 9;

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
    const double tolerance = 32. * std::numeric_limits<double>::epsilon() * scale;
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
  const double tolerance = 64. * std::numeric_limits<double>::epsilon() *
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

  const double endpointLength = std::hypot(length02.value, dz02.value);
  const double sineTheta = length02.value / endpointLength;
  if (!std::isfinite(phiTilde.value) || !std::isfinite(thetaTilde.value) ||
      !std::isfinite(rhoPhi.value) || rhoPhi.value == 0. ||
      !std::isfinite(rhoTheta.value) || !std::isfinite(sineTheta) || sineTheta <= 0.) {
    return false;
  }
  result = {phiTilde, thetaTilde, rhoPhi, rhoTheta,
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

} // namespace

bool makeCylinderTripletFitObservation(const SurfaceDescriptor& surface,
                                       const SurfaceMeasurement& measurement,
                                       TripletFitObservation& observation) noexcept
{
  if (surface.kind != SurfaceKind::Cylinder ||
      !std::isfinite(surface.referenceCoordinate) || surface.referenceCoordinate <= 0.f ||
      !std::isfinite(measurement.frame.q) || !std::isfinite(measurement.frame.u) ||
      !std::isfinite(measurement.frame.v) || !std::isfinite(measurement.frame.frameAngle)) {
    return false;
  }
  const double sine = std::sin(measurement.frame.frameAngle);
  const double cosine = std::cos(measurement.frame.frameAngle);
  const double varianceU = measurement.covariance.uu;
  const double covarianceUV = measurement.covariance.uv;
  const double varianceV = measurement.covariance.vv;
  const TripletFitObservation scratch{
    {measurement.frame.q * cosine - measurement.frame.u * sine,
     measurement.frame.q * sine + measurement.frame.u * cosine,
     measurement.frame.v},
    {sine * sine * varianceU,
     -sine * cosine * varianceU,
     -sine * covarianceUV,
     cosine * cosine * varianceU,
     cosine * covarianceUV,
     varianceV}};
  if (!observationIsValid(scratch)) {
    return false;
  }
  observation = scratch;
  return true;
}

bool makeDiskTripletFitObservation(const SurfaceDescriptor& surface,
                                   const SurfaceMeasurement& measurement,
                                   TripletFitObservation& observation) noexcept
{
  if (surface.kind != SurfaceKind::Disk || !std::isfinite(surface.referenceCoordinate)) {
    return false;
  }
  const TripletFitObservation scratch{
    {measurement.global.x, measurement.global.y, surface.referenceCoordinate},
    {measurement.covariance.uu, measurement.covariance.uv, 0.,
     measurement.covariance.vv, 0., 0.}};
  if (!observationIsValid(scratch)) {
    return false;
  }
  observation = scratch;
  return true;
}

bool makeTripletFitObservation(const SurfaceDescriptor& surface,
                               const SurfaceMeasurement& measurement,
                               TripletFitObservation& observation) noexcept
{
  using Builder = bool (*)(const SurfaceDescriptor&, const SurfaceMeasurement&,
                           TripletFitObservation&) noexcept;
  static constexpr std::array<Builder, 2> builders{
    makeCylinderTripletFitObservation,
    makeDiskTripletFitObservation};
  const auto kindIndex = static_cast<std::size_t>(surface.kind);
  return kindIndex < builders.size() && builders[kindIndex](surface, measurement, observation);
}

bool fitLocalTripletUniformSolenoid(
  const std::array<TripletFitObservation, 3>& observations,
  const TripletFitProcessNoise& processNoise,
  LocalTripletFitResult& result) noexcept
{
  if (!std::all_of(observations.begin(), observations.end(), observationIsValid) ||
      !std::isfinite(processNoise.angularVariance) || processNoise.angularVariance < 0.) {
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
    }
    gammaThetaTheta += covarianceContraction(gradientTheta, observations[hit].covariance, gradientTheta);
    gammaThetaPhi += covarianceContraction(gradientTheta, observations[hit].covariance, gradientPhi);
    gammaPhiPhi += covarianceContraction(gradientPhi, observations[hit].covariance, gradientPhi);
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
            gammaPhiPhi};
  return true;
}

double fittedTripletTransverseMomentum(const LocalTripletFitResult& result,
                                       double bz) noexcept
{
  return fittedTripletTransverseMomentum(result, bz, 1);
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
