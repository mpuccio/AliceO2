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
#include "ITSMFTTracking/MaterialPhysics.h"

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

bool materialIsValid(const NominalSurfaceMaterial& material) noexcept
{
  return std::isfinite(material.xOverX0) && material.xOverX0 >= 0.f &&
         std::isfinite(material.arealDensityGPerCm2) && material.arealDensityGPerCm2 >= 0.f;
}

bool makeMiddleDirection(const std::array<TripletFitObservation, 3>& observations,
                         double transverseCurvature,
                         std::array<double, 3>& direction) noexcept
{
  std::array<std::array<double, 3>, 2> segmentDirections{};
  for (std::size_t segment = 0; segment < 2; ++segment) {
    const double dx = observations[segment + 1].position[0] - observations[segment].position[0];
    const double dy = observations[segment + 1].position[1] - observations[segment].position[1];
    const double dz = observations[segment + 1].position[2] - observations[segment].position[2];
    const double chord = std::hypot(dx, dy);
    const double halfSine = 0.5 * transverseCurvature * chord;
    if (!std::isfinite(chord) || chord <= 0. || !std::isfinite(halfSine) || std::abs(halfSine) >= 1.) {
      return false;
    }

    double halfAngle{0.};
    double arcOverChord{1.};
    if (std::abs(halfSine) < 1.e-4) {
      const double halfSine2 = halfSine * halfSine;
      halfAngle = halfSine * (1. + halfSine2 / 6. + 3. * halfSine2 * halfSine2 / 40.);
      arcOverChord = 1. + halfSine2 / 6. + 3. * halfSine2 * halfSine2 / 40.;
    } else {
      halfAngle = std::asin(halfSine);
      arcOverChord = halfAngle / halfSine;
    }
    const double transverseArc = chord * arcOverChord;
    const double tanLambda = dz / transverseArc;
    const double inverse3DScale = 1. / std::sqrt(1. + tanLambda * tanLambda);
    const double rotation = segment == 0 ? halfAngle : -halfAngle;
    const double cosine = std::cos(rotation);
    const double sine = std::sin(rotation);
    const double unitX = dx / chord;
    const double unitY = dy / chord;
    segmentDirections[segment] = {
      (cosine * unitX - sine * unitY) * inverse3DScale,
      (sine * unitX + cosine * unitY) * inverse3DScale,
      tanLambda * inverse3DScale};
  }

  const std::array<double, 3> sum{
    segmentDirections[0][0] + segmentDirections[1][0],
    segmentDirections[0][1] + segmentDirections[1][1],
    segmentDirections[0][2] + segmentDirections[1][2]};
  const double norm = std::sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);
  if (!std::isfinite(norm) || norm <= 0.) {
    return false;
  }
  direction = {sum[0] / norm, sum[1] / norm, sum[2] / norm};
  return true;
}

bool normalize(const std::array<double, 3>& input, std::array<double, 3>& output) noexcept
{
  const double norm = std::sqrt(input[0] * input[0] + input[1] * input[1] + input[2] * input[2]);
  if (!std::isfinite(norm) || norm <= 0.) {
    return false;
  }
  output = {input[0] / norm, input[1] / norm, input[2] / norm};
  return true;
}

bool fittedMomentum(const LocalTripletFitResult& fit, double bz, uint8_t absCharge,
                    double& momentum, double& transverseMomentum) noexcept
{
  if (absCharge == 0 || !std::isfinite(bz) || bz == 0. ||
      !std::isfinite(fit.curvature) || !std::isfinite(fit.referenceSinTheta) ||
      fit.referenceSinTheta <= 0. || fit.referenceSinTheta > 1.) {
    return false;
  }
  if (fit.curvature == 0.) {
    momentum = std::numeric_limits<double>::infinity();
    transverseMomentum = std::numeric_limits<double>::infinity();
    return true;
  }
  momentum = static_cast<double>(absCharge) *
             std::abs(static_cast<double>(o2::constants::math::B2C) * bz / fit.curvature);
  transverseMomentum = momentum * fit.referenceSinTheta;
  return std::isfinite(momentum) && momentum > 0. &&
         std::isfinite(transverseMomentum) && transverseMomentum > 0.;
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

bool makeCylinderTripletFitMaterial(const SurfaceDescriptor& surface,
                                    const SurfaceMeasurement& measurement,
                                    TripletFitMaterial& material) noexcept
{
  if (surface.kind != SurfaceKind::Cylinder || !materialIsValid(surface.material) ||
      !std::isfinite(measurement.frame.frameAngle)) {
    return false;
  }
  const TripletFitMaterial scratch{
    {std::cos(measurement.frame.frameAngle), std::sin(measurement.frame.frameAngle), 0.},
    surface.material};
  material = scratch;
  return true;
}

bool makeDiskTripletFitMaterial(const SurfaceDescriptor& surface,
                                const SurfaceMeasurement&,
                                TripletFitMaterial& material) noexcept
{
  if (surface.kind != SurfaceKind::Disk || !materialIsValid(surface.material)) {
    return false;
  }
  material = {{0., 0., 1.}, surface.material};
  return true;
}

bool makeTripletFitMaterial(const SurfaceDescriptor& surface,
                            const SurfaceMeasurement& measurement,
                            TripletFitMaterial& material) noexcept
{
  using Builder = bool (*)(const SurfaceDescriptor&, const SurfaceMeasurement&,
                           TripletFitMaterial&) noexcept;
  static constexpr std::array<Builder, 2> builders{
    makeCylinderTripletFitMaterial,
    makeDiskTripletFitMaterial};
  const auto kindIndex = static_cast<std::size_t>(surface.kind);
  return kindIndex < builders.size() && builders[kindIndex](surface, measurement, material);
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

MaterialTripletFitStatus fitLocalTripletWithMaterial(
  const std::array<TripletFitObservation, 3>& observations,
  const TripletFitMaterial& middleMaterial,
  double bz, uint8_t absCharge, o2::track::PID pid,
  MaterialAwareTripletFitResult& result) noexcept
{
  if (!materialIsValid(middleMaterial.nominal) || absCharge == 0 ||
      !std::isfinite(bz) || bz == 0. || pid.getID() >= o2::track::PID::NIDsTot ||
      (pid.getMass() == 0.f && absCharge != 0)) {
    return MaterialTripletFitStatus::InvalidInput;
  }

  std::array<double, 3> normal{};
  if (!normalize(middleMaterial.unitNormal, normal)) {
    return MaterialTripletFitStatus::InvalidInput;
  }

  LocalTripletFitResult fit{};
  if (!fitLocalTripletUniformSolenoid(observations, {0.}, fit)) {
    return MaterialTripletFitStatus::LocalFitFailure;
  }
  std::array<double, 3> middleDirection{};
  if (!makeMiddleDirection(observations, fit.referenceTransverseCurvature, middleDirection)) {
    return MaterialTripletFitStatus::LocalFitFailure;
  }
  const double incidenceCosine = std::abs(middleDirection[0] * normal[0] +
                                          middleDirection[1] * normal[1] +
                                          middleDirection[2] * normal[2]);
  if (!std::isfinite(incidenceCosine) || incidenceCosine <= 0.) {
    return MaterialTripletFitStatus::MaterialEvaluationFailure;
  }
  const double pathXOverX0 = static_cast<double>(middleMaterial.nominal.xOverX0) / incidenceCosine;
  const double pathArealDensity = static_cast<double>(middleMaterial.nominal.arealDensityGPerCm2) / incidenceCosine;
  if (!std::isfinite(pathXOverX0) || pathXOverX0 > std::numeric_limits<float>::max() ||
      !std::isfinite(pathArealDensity) || pathArealDensity > std::numeric_limits<float>::max()) {
    return MaterialTripletFitStatus::MaterialEvaluationFailure;
  }
  const material::IntegratedMaterialBudget pathMaterial{
    static_cast<float>(pathXOverX0), static_cast<float>(pathArealDensity)};

  double momentum{0.};
  double transverseMomentum{0.};
  if (!fittedMomentum(fit, bz, absCharge, momentum, transverseMomentum)) {
    return MaterialTripletFitStatus::MomentumUnresolved;
  }
  if (std::isinf(momentum)) {
    result = {fit, momentum, transverseMomentum, 0., incidenceCosine, 1};
    return MaterialTripletFitStatus::Success;
  }

  constexpr uint8_t DirectIterations = 8;
  constexpr uint8_t MaximumIterations = 128;
  constexpr double RelativeFixedPointTolerance = 4.76837158203125e-7;
  double angularVariance{0.};
  for (uint8_t iteration = 1; iteration <= MaximumIterations; ++iteration) {
    if (momentum > std::numeric_limits<float>::max()) {
      return MaterialTripletFitStatus::MaterialEvaluationFailure;
    }
    const auto materialResult = material::calculateMaterialPhysics(
      static_cast<float>(momentum), pid, absCharge,
      material::MaterialTraversalDirection::AlongMomentum, pathMaterial);
    if (materialResult.failure == material::MaterialFailureReason::MomentumBelowMinimum) {
      return MaterialTripletFitStatus::MomentumBelowMaterialRange;
    }
    if (materialResult.failure == material::MaterialFailureReason::ExcessiveScattering) {
      return MaterialTripletFitStatus::ExcessiveScattering;
    }
    if (materialResult.failure == material::MaterialFailureReason::StoppedInMaterial) {
      return MaterialTripletFitStatus::StoppedInMaterial;
    }
    if (!materialResult.ok() || !std::isfinite(materialResult.highlandTheta2Rad2) ||
        materialResult.highlandTheta2Rad2 < 0.f) {
      return MaterialTripletFitStatus::MaterialEvaluationFailure;
    }

    const double targetAngularVariance = materialResult.highlandTheta2Rad2;
    if (targetAngularVariance == angularVariance ||
        (angularVariance > 0. &&
         std::abs(targetAngularVariance - angularVariance) <=
           RelativeFixedPointTolerance * std::max(targetAngularVariance, angularVariance))) {
      result = {fit, momentum, transverseMomentum, angularVariance,
                incidenceCosine, static_cast<uint8_t>(std::max<int>(1, iteration - 1))};
      return MaterialTripletFitStatus::Success;
    }
    const double nextAngularVariance = iteration <= DirectIterations
                                         ? targetAngularVariance
                                         : 0.5 * (angularVariance + targetAngularVariance);
    LocalTripletFitResult nextFit{};
    if (!fitLocalTripletUniformSolenoid(observations, {nextAngularVariance}, nextFit)) {
      return MaterialTripletFitStatus::LocalFitFailure;
    }
    double nextMomentum{0.};
    double nextTransverseMomentum{0.};
    if (!fittedMomentum(nextFit, bz, absCharge, nextMomentum, nextTransverseMomentum)) {
      return MaterialTripletFitStatus::MomentumUnresolved;
    }
    if (std::isinf(nextMomentum) ||
        (iteration <= DirectIterations &&
         std::abs(nextMomentum - momentum) <= RelativeFixedPointTolerance * std::max(nextMomentum, momentum))) {
      result = {nextFit, nextMomentum, nextTransverseMomentum,
                nextAngularVariance, incidenceCosine, iteration};
      return MaterialTripletFitStatus::Success;
    }
    fit = nextFit;
    momentum = nextMomentum;
    transverseMomentum = nextTransverseMomentum;
    angularVariance = nextAngularVariance;
  }
  return MaterialTripletFitStatus::NoConvergence;
}

} // namespace o2::itsmft::tracking
