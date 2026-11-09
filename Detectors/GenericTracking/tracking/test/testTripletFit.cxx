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

/// \file testTripletFit.cxx
/// \brief The copied ITS triplet geometry, the kinks' gradient by hand, the
///        covariance and the fit.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Triplet Fit
#include <boost/test/unit_test.hpp>

#include "ToyBarrel4Disc8.h"
#include "Oracles.h"
#include "GenericTracking/Tracking/Triplet.h"
#include "GenericTracking/Tracking/TripletFit.h"
#include "GenericTracking/Study/Validation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <cstddef>
#include <cstdint>

namespace o2::gtrk
{
namespace
{

using detectors::ToyBarrel4Disc8;

constexpr std::size_t kCoordinates = 9U;

// The reference below is ITS's makeTripletFitFactor geometry, kept in its
// original nine-derivative form so it is independent of GTRK's gradient,
// and evaluated in double: the float port is compared with the true value,
// not with a second float evaluation whose rounding differs.
struct ReferenceDual {
  double value = 0.;
  std::array<double, kCoordinates> derivative{};

  ReferenceDual() = default;
  explicit ReferenceDual(double valueIn) : value(valueIn) {}

  static ReferenceDual variable(double value, std::size_t coordinate)
  {
    ReferenceDual result{value};
    result.derivative[coordinate] = 1.;
    return result;
  }
};

ReferenceDual operator+(const ReferenceDual& left, const ReferenceDual& right)
{
  ReferenceDual result{left.value + right.value};
  for (std::size_t i = 0; i < kCoordinates; ++i) {
    result.derivative[i] = left.derivative[i] + right.derivative[i];
  }
  return result;
}

ReferenceDual operator-(const ReferenceDual& left, const ReferenceDual& right)
{
  ReferenceDual result{left.value - right.value};
  for (std::size_t i = 0; i < kCoordinates; ++i) {
    result.derivative[i] = left.derivative[i] - right.derivative[i];
  }
  return result;
}

ReferenceDual operator*(const ReferenceDual& left, const ReferenceDual& right)
{
  ReferenceDual result{left.value * right.value};
  for (std::size_t i = 0; i < kCoordinates; ++i) {
    result.derivative[i] = left.derivative[i] * right.value + left.value * right.derivative[i];
  }
  return result;
}

ReferenceDual operator/(const ReferenceDual& left, const ReferenceDual& right)
{
  const double inverse = 1. / right.value;
  ReferenceDual result{left.value * inverse};
  for (std::size_t i = 0; i < kCoordinates; ++i) {
    result.derivative[i] = (left.derivative[i] - result.value * right.derivative[i]) * inverse;
  }
  return result;
}

ReferenceDual referenceSqrt(const ReferenceDual& argument)
{
  const double root = std::sqrt(argument.value);
  ReferenceDual result{root};
  for (std::size_t i = 0; i < kCoordinates; ++i) {
    result.derivative[i] = 0.5 * argument.derivative[i] / root;
  }
  return result;
}

ReferenceDual referenceAsin(const ReferenceDual& argument)
{
  ReferenceDual result{std::asin(argument.value)};
  const double scale = 1. / std::sqrt(1. - argument.value * argument.value);
  for (std::size_t i = 0; i < kCoordinates; ++i) {
    result.derivative[i] = scale * argument.derivative[i];
  }
  return result;
}

ReferenceDual referenceAtan2(const ReferenceDual& y, const ReferenceDual& x)
{
  ReferenceDual result{std::atan2(y.value, x.value)};
  const double denominator = x.value * x.value + y.value * y.value;
  for (std::size_t i = 0; i < kCoordinates; ++i) {
    result.derivative[i] = (x.value * y.derivative[i] - y.value * x.derivative[i]) / denominator;
  }
  return result;
}

struct ReferenceSegment {
  ReferenceDual bendingAngle;
  ReferenceDual transverseArcLength;
  ReferenceDual cotangentTheta;
  ReferenceDual sineTheta;
  ReferenceDual cosineTheta;
  ReferenceDual index;
};

ReferenceSegment referenceSegment(const ReferenceDual& curvature, const ReferenceDual& chord, const ReferenceDual& deltaZ)
{
  const ReferenceDual halfSine = ReferenceDual{0.5} * curvature * chord;
  const ReferenceDual halfSine2 = halfSine * halfSine;
  const ReferenceDual halfSine4 = halfSine2 * halfSine2;
  const ReferenceDual halfAngle = referenceAsin(halfSine);
  ReferenceDual asinOverArgument;
  ReferenceDual angleCotangent;
  if (std::abs(halfSine.value) < 1.e-4) {
    asinOverArgument = ReferenceDual{1.} + halfSine2 * ReferenceDual{1. / 6.} + halfSine4 * ReferenceDual{3. / 40.};
    angleCotangent = ReferenceDual{1.} - halfSine2 * ReferenceDual{1. / 3.} - halfSine4 * ReferenceDual{2. / 15.};
  } else {
    asinOverArgument = halfAngle / halfSine;
    angleCotangent = halfAngle * referenceSqrt(ReferenceDual{1.} - halfSine2) / halfSine;
  }
  const ReferenceDual arc = chord * asinOverArgument;
  const ReferenceDual cotangent = deltaZ / arc;
  const ReferenceDual sine = ReferenceDual{1.} / referenceSqrt(ReferenceDual{1.} + cotangent * cotangent);
  const ReferenceDual cosine = cotangent * sine;
  const ReferenceDual index = ReferenceDual{1.} / (angleCotangent * sine * sine + cosine * cosine);
  return {
    .bendingAngle = ReferenceDual{2.} * halfAngle,
    .transverseArcLength = arc,
    .cotangentTheta = cotangent,
    .sineTheta = sine,
    .cosineTheta = cosine,
    .index = index,
  };
}

struct ReferenceGeometry {
  ReferenceDual theta;
  ReferenceDual phi;
  ReferenceDual rhoTheta;
  ReferenceDual rhoPhi;
};

ReferenceGeometry referenceGeometry(const std::array<PairPoint, 3>& measurements)
{
  std::array<std::array<ReferenceDual, 3>, 3> point{};
  for (std::size_t hit = 0; hit < measurements.size(); ++hit) {
    const double position[3]{static_cast<double>(measurements[hit].x), static_cast<double>(measurements[hit].y), static_cast<double>(measurements[hit].z)};
    for (std::size_t coordinate = 0; coordinate < 3U; ++coordinate) {
      point[hit][coordinate] = ReferenceDual::variable(position[coordinate], 3U * hit + coordinate);
    }
  }
  const ReferenceDual dx01 = point[1][0] - point[0][0];
  const ReferenceDual dy01 = point[1][1] - point[0][1];
  const ReferenceDual dz01 = point[1][2] - point[0][2];
  const ReferenceDual dx12 = point[2][0] - point[1][0];
  const ReferenceDual dy12 = point[2][1] - point[1][1];
  const ReferenceDual dz12 = point[2][2] - point[1][2];
  const ReferenceDual dx02 = point[2][0] - point[0][0];
  const ReferenceDual dy02 = point[2][1] - point[0][1];
  const ReferenceDual length01 = referenceSqrt(dx01 * dx01 + dy01 * dy01);
  const ReferenceDual length12 = referenceSqrt(dx12 * dx12 + dy12 * dy12);
  const ReferenceDual length02 = referenceSqrt(dx02 * dx02 + dy02 * dy02);
  const ReferenceDual transverseCurvature = ReferenceDual{2.} * (dx01 * dy12 - dy01 * dx12) / (length01 * length12 * length02);
  const auto first = referenceSegment(transverseCurvature, length01, dz01);
  const auto second = referenceSegment(transverseCurvature, length12, dz12);
  const ReferenceDual theta01 = referenceAtan2(first.transverseArcLength, dz01);
  const ReferenceDual theta12 = referenceAtan2(second.transverseArcLength, dz12);
  const ReferenceDual phi = ReferenceDual{0.5} * (first.bendingAngle * first.index + second.bendingAngle * second.index);
  const ReferenceDual theta = theta12 - theta01 + (ReferenceDual{1.} - second.index) * second.cotangentTheta - (ReferenceDual{1.} - first.index) * first.cotangentTheta;
  const ReferenceDual rhoPhi = ReferenceDual{-0.5} * (first.transverseArcLength * first.index / first.sineTheta + second.transverseArcLength * second.index / second.sineTheta);
  ReferenceDual rhoTheta;
  const double maximumHalfSine = 0.5 * std::abs(transverseCurvature.value) * std::max(length01.value, length12.value);
  if (maximumHalfSine < 1.e-4) {
    rhoTheta = transverseCurvature * (length12 * length12 * second.cosineTheta - length01 * length01 * first.cosineTheta) / ReferenceDual{12.};
  } else {
    rhoTheta = ((ReferenceDual{1.} - first.index) * first.cotangentTheta / first.sineTheta - (ReferenceDual{1.} - second.index) * second.cotangentTheta / second.sineTheta) / transverseCurvature;
  }
  return {.theta = theta, .phi = phi, .rhoTheta = rhoTheta, .rhoPhi = rhoPhi};
}

struct ReferenceFactor {
  TripletKinkVector psi;
  TripletKinkVector rho;
  TripletKinkCovariance covariance;
};

ReferenceFactor referenceFactor(const std::array<PairPoint, 3>& points, Scalar angularVariance)
{
  const auto geometry = referenceGeometry(points);
  const Scalar kappaReference = -static_cast<Scalar>(geometry.phi.value) / static_cast<Scalar>(geometry.rhoPhi.value);
  const Scalar dx = points[2].x - points[0].x;
  const Scalar dy = points[2].y - points[0].y;
  const Scalar dz = points[2].z - points[0].z;
  const Scalar transverse2 = dx * dx + dy * dy;
  const Scalar sineTheta2 = transverse2 / (transverse2 + dz * dz);
  ReferenceFactor result{
    .psi = {.theta = static_cast<Scalar>(geometry.theta.value), .phi = static_cast<Scalar>(geometry.phi.value)},
    .rho = {.theta = static_cast<Scalar>(geometry.rhoTheta.value), .phi = static_cast<Scalar>(geometry.rhoPhi.value)},
    .covariance = {.thetaTheta = angularVariance, .thetaPhi = 0.F, .phiPhi = angularVariance / sineTheta2},
  };
  for (std::size_t hit = 0; hit < points.size(); ++hit) {
    const Scalar inverseR = 1.F / points[hit].r;
    const Scalar radialX = points[hit].x * inverseR;
    const Scalar radialY = points[hit].y * inverseR;
    const Scalar tangentX = -points[hit].y * inverseR;
    const Scalar tangentY = points[hit].x * inverseR;
    const std::size_t coordinate = 3U * hit;
    const Scalar thetaX = static_cast<Scalar>(geometry.theta.derivative[coordinate]) + kappaReference * static_cast<Scalar>(geometry.rhoTheta.derivative[coordinate]);
    const Scalar thetaY = static_cast<Scalar>(geometry.theta.derivative[coordinate + 1U]) + kappaReference * static_cast<Scalar>(geometry.rhoTheta.derivative[coordinate + 1U]);
    const Scalar thetaZ = static_cast<Scalar>(geometry.theta.derivative[coordinate + 2U]) + kappaReference * static_cast<Scalar>(geometry.rhoTheta.derivative[coordinate + 2U]);
    const Scalar phiX = static_cast<Scalar>(geometry.phi.derivative[coordinate]) + kappaReference * static_cast<Scalar>(geometry.rhoPhi.derivative[coordinate]);
    const Scalar phiY = static_cast<Scalar>(geometry.phi.derivative[coordinate + 1U]) + kappaReference * static_cast<Scalar>(geometry.rhoPhi.derivative[coordinate + 1U]);
    const Scalar phiZ = static_cast<Scalar>(geometry.phi.derivative[coordinate + 2U]) + kappaReference * static_cast<Scalar>(geometry.rhoPhi.derivative[coordinate + 2U]);
    const Scalar thetaRadial = thetaX * radialX + thetaY * radialY;
    const Scalar phiRadial = phiX * radialX + phiY * radialY;
    const Scalar thetaTangent = thetaX * tangentX + thetaY * tangentY;
    const Scalar phiTangent = phiX * tangentX + phiY * tangentY;
    const Scalar variance[3]{points[hit].varR, points[hit].varRPhi, points[hit].varZ};
    const Scalar thetaDerivative[3]{thetaRadial, thetaTangent, thetaZ};
    const Scalar phiDerivative[3]{phiRadial, phiTangent, phiZ};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      result.covariance.thetaTheta += variance[axis] * thetaDerivative[axis] * thetaDerivative[axis];
      result.covariance.thetaPhi += variance[axis] * thetaDerivative[axis] * phiDerivative[axis];
      result.covariance.phiPhi += variance[axis] * phiDerivative[axis] * phiDerivative[axis];
    }
  }
  return result;
}

struct ReferenceFit {
  double curvature = 0.;
  double varCurvature = 0.;
  double chi2 = 0.;
};

/// ITS's fitAdjacentTripletFactors, in double: each triplet's nine-coordinate
/// gradients at its own reference curvature, projected on the measurement
/// axes and summed into one 4x4 kink covariance over the four hits, then
/// solved by Gauss-Jordan elimination.
ReferenceFit referenceAdjacentFit(const std::array<PairPoint, 4>& points, double varianceSecond, double varianceThird)
{
  const std::array<ReferenceGeometry, 2> geometry{referenceGeometry({points[0], points[1], points[2]}), referenceGeometry({points[1], points[2], points[3]})};
  const std::array<double, 4> psi{geometry[0].theta.value, geometry[0].phi.value, geometry[1].theta.value, geometry[1].phi.value};
  const std::array<double, 4> rho{geometry[0].rhoTheta.value, geometry[0].rhoPhi.value, geometry[1].rhoTheta.value, geometry[1].rhoPhi.value};
  const auto sine2 = [&](std::size_t from, std::size_t to) {
    const double dx = static_cast<double>(points[to].x) - static_cast<double>(points[from].x);
    const double dy = static_cast<double>(points[to].y) - static_cast<double>(points[from].y);
    const double dz = static_cast<double>(points[to].z) - static_cast<double>(points[from].z);
    return (dx * dx + dy * dy) / (dx * dx + dy * dy + dz * dz);
  };
  std::array<std::array<double, 4>, 4> covariance{};
  covariance[0][0] = varianceSecond;
  covariance[1][1] = varianceSecond / sine2(0, 2);
  covariance[2][2] = varianceThird;
  covariance[3][3] = varianceThird / sine2(1, 3);
  for (std::size_t hit = 0; hit < points.size(); ++hit) {
    const double r = std::hypot(static_cast<double>(points[hit].x), static_cast<double>(points[hit].y));
    const double radialX = static_cast<double>(points[hit].x) / r;
    const double radialY = static_cast<double>(points[hit].y) / r;
    // gradient[kink][axis], axes radial, r-phi, z
    std::array<std::array<double, 3>, 4> gradient{};
    for (std::size_t triplet = 0; triplet < 2U; ++triplet) {
      if (hit < triplet || hit > triplet + 2U) {
        continue;
      }
      const auto& g = geometry[triplet];
      const double kappa = -g.phi.value / g.rhoPhi.value;
      const std::size_t coordinate = 3U * (hit - triplet);
      for (std::size_t angle = 0; angle < 2U; ++angle) {
        const auto& kink = angle == 0U ? g.theta : g.phi;
        const auto& slope = angle == 0U ? g.rhoTheta : g.rhoPhi;
        const double dX = kink.derivative[coordinate] + kappa * slope.derivative[coordinate];
        const double dY = kink.derivative[coordinate + 1U] + kappa * slope.derivative[coordinate + 1U];
        const double dZ = kink.derivative[coordinate + 2U] + kappa * slope.derivative[coordinate + 2U];
        gradient[2U * triplet + angle] = {dX * radialX + dY * radialY, -dX * radialY + dY * radialX, dZ};
      }
    }
    const std::array<double, 3> variance{static_cast<double>(points[hit].varR), static_cast<double>(points[hit].varRPhi), static_cast<double>(points[hit].varZ)};
    for (std::size_t row = 0; row < 4U; ++row) {
      for (std::size_t column = 0; column < 4U; ++column) {
        for (std::size_t axis = 0; axis < 3U; ++axis) {
          covariance[row][column] += variance[axis] * gradient[row][axis] * gradient[column][axis];
        }
      }
    }
  }
  // [V | psi | rho] reduced to [1 | V^-1 psi | V^-1 rho]
  std::array<std::array<double, 6>, 4> augmented{};
  for (std::size_t row = 0; row < 4U; ++row) {
    for (std::size_t column = 0; column < 4U; ++column) {
      augmented[row][column] = covariance[row][column];
    }
    augmented[row][4] = psi[row];
    augmented[row][5] = rho[row];
  }
  for (std::size_t pivot = 0; pivot < 4U; ++pivot) {
    const double scale = augmented[pivot][pivot];
    for (auto& entry : augmented[pivot]) {
      entry /= scale;
    }
    for (std::size_t row = 0; row < 4U; ++row) {
      if (row != pivot) {
        const double factor = augmented[row][pivot];
        for (std::size_t column = 0; column < 6U; ++column) {
          augmented[row][column] -= factor * augmented[pivot][column];
        }
      }
    }
  }
  double psiPsi = 0.;
  double rhoPsi = 0.;
  double rhoRho = 0.;
  for (std::size_t row = 0; row < 4U; ++row) {
    psiPsi += psi[row] * augmented[row][4];
    rhoPsi += rho[row] * augmented[row][4];
    rhoRho += rho[row] * augmented[row][5];
  }
  return {.curvature = -rhoPsi / rhoRho, .varCurvature = 1. / rhoRho, .chi2 = psiPsi - (rhoPsi * rhoPsi / rhoRho)};
}

PairPoint helixPoint(double curvature, double arc, double cotangentTheta)
{
  const double turn = curvature * arc;
  const double x = std::sin(turn) / curvature;
  const double y = (1. - std::cos(turn)) / curvature;
  return {
    .x = static_cast<Scalar>(x),
    .y = static_cast<Scalar>(y),
    .z = static_cast<Scalar>(cotangentTheta * arc),
    .r = static_cast<Scalar>(std::hypot(x, y)),
    .varZ = 4.e-6F,
    .varR = 2.e-6F,
    .varRPhi = 3.e-6F,
  };
}

std::array<PairPoint, 3> helixPoints(double curvature, double cotangentTheta, const std::array<double, 3>& arcs) { return {helixPoint(curvature, arcs[0], cotangentTheta), helixPoint(curvature, arcs[1], cotangentTheta), helixPoint(curvature, arcs[2], cotangentTheta)}; }

void checkClose(Scalar actual, Scalar expected, Scalar relativeTolerance, Scalar absoluteTolerance = 0.F)
{
  const Scalar scale = std::max(std::abs(expected), absoluteTolerance);
  BOOST_TEST(std::abs(actual - expected) <= relativeTolerance * scale + absoluteTolerance);
}

// An exact helix gives chi2 = 0 up to the float cancellation of psi' V^-1 psi
// against the removed curvature term, which are of order 1e4 here; any real
// cut is of order n^2 = 9.
constexpr Scalar kExactHelixChi2 = 1.e-2F;

/// The fit's curvature is the triplet fit's own parameter: the 3D curvature,
/// proportional to q/p, which is the transverse curvature times sin(theta).
void checkHelix(double ptGeV, const std::array<double, 3>& arcs, double cotangentTheta)
{
  constexpr double kBzKiloGauss = 5.;
  constexpr double kB2C = 0.000299792458;
  const double curvature = kB2C * kBzKiloGauss / ptGeV;
  const double sineTheta = 1. / std::sqrt(1. + (cotangentTheta * cotangentTheta));
  const auto points = helixPoints(curvature, cotangentTheta, arcs);
  TripletFitFactor factor;
  BOOST_REQUIRE(study::makeTripletFitFactor(points[0], points[1], points[2], 2.e-6F, factor));
  TripletFitResult fit;
  BOOST_REQUIRE(fitTriplet(factor, fit));
  BOOST_TEST_CONTEXT("pT " << ptGeV << " arcs " << arcs[0] << ',' << arcs[1] << ',' << arcs[2])
  {
    checkClose(fit.curvature, static_cast<Scalar>(curvature * sineTheta), 3.e-3F, 2.e-7F);
    BOOST_TEST(fit.varCurvature > 0.F);
    BOOST_TEST(fit.chi2 < kExactHelixChi2);
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(KinkGradientMatchesTheCopiedItsFactor)
{
  auto points = helixPoints(0.012, 0.7, {2.1, 5.4, 9.8});
  points[0].varR = 1.e-6F;
  points[0].varRPhi = 2.e-6F;
  points[0].varZ = 3.e-6F;
  points[1].varR = 4.e-6F;
  points[1].varRPhi = 5.e-6F;
  points[1].varZ = 6.e-6F;
  points[2].varR = 7.e-6F;
  points[2].varRPhi = 8.e-6F;
  points[2].varZ = 9.e-6F;
  constexpr Scalar kAngularVariance = 1.7e-5F;
  const auto reference = referenceFactor(points, kAngularVariance);
  const auto geometry = kinkGradient(points[0], points[1], points[2]).kinks;
  BOOST_REQUIRE(geometry.valid);
  checkClose(geometry.psi.phi, reference.psi.phi, 2.e-5F, 1.e-7F);
  checkClose(geometry.rho.phi, reference.rho.phi, 2.e-5F, 1.e-7F);
  // The theta kink is a difference of two dip angles near one radian, and its
  // slope a difference of `1 - index`, about (kappa d/2)^2/3 here: both cancel.
  // Float resolves the kink, psi + kappa rho, to a few ulp of those angles.
  constexpr Scalar kThetaKinkResolution = 5.e-7F;
  constexpr Scalar kCurvature = 0.012F;
  checkClose(geometry.psi.theta, reference.psi.theta, 0.F, kThetaKinkResolution);
  checkClose(geometry.rho.theta, reference.rho.theta, 0.F, kThetaKinkResolution / kCurvature);

  TripletFitFactor factor;
  BOOST_REQUIRE(study::makeTripletFitFactor(points[0], points[1], points[2], kAngularVariance, factor));
  checkClose(factor.covariance.thetaTheta, reference.covariance.thetaTheta, 5.e-5F, 1.e-9F);
  checkClose(factor.covariance.thetaPhi, reference.covariance.thetaPhi, 5.e-5F, 1.e-9F);
  checkClose(factor.covariance.phiPhi, reference.covariance.phiPhi, 5.e-5F, 1.e-9F);

  const TripletFitFactor referenceFitFactor{.psi = reference.psi, .rho = reference.rho, .covariance = reference.covariance};
  TripletFitResult actualFit;
  TripletFitResult referenceFit;
  BOOST_REQUIRE(fitTriplet(factor, actualFit));
  BOOST_REQUIRE(fitTriplet(referenceFitFactor, referenceFit));
  checkClose(actualFit.curvature, referenceFit.curvature, 5.e-5F, 1.e-8F);
  checkClose(actualFit.varCurvature, referenceFit.varCurvature, 5.e-5F, 1.e-10F);
  // Both are an exact helix's chi2: zero up to the same float cancellation.
  checkClose(actualFit.chi2, referenceFit.chi2, 0.F, kExactHelixChi2);
}

BOOST_AUTO_TEST_CASE(HelixCurvatureAcrossLayerKindsAndMomenta)
{
  const std::array<std::array<double, 3>, 3> geometries{{{{2., 4., 7.}}, {{4., 8., 12.}}, {{2., 6., 11.}}}};
  const double cotangents[3]{0.25, 1.5, 0.8};
  const double momenta[3]{0.1, 1., 10.};
  for (std::size_t geometry = 0; geometry < geometries.size(); ++geometry) {
    for (const double momentum : momenta) {
      checkHelix(momentum, geometries[geometry], cotangents[geometry]);
    }
  }
}

BOOST_AUTO_TEST_CASE(Chi2RespondsToAMiddleHitDisplacement)
{
  auto points = helixPoints(0.004, 0.6, {2., 5., 9.});
  TripletFitFactor exactFactor;
  BOOST_REQUIRE(study::makeTripletFitFactor(points[0], points[1], points[2], 1.e-6F, exactFactor));
  TripletFitResult exact;
  BOOST_REQUIRE(fitTriplet(exactFactor, exact));
  points[1].z += 0.05F;
  TripletFitFactor displacedFactor;
  BOOST_REQUIRE(study::makeTripletFitFactor(points[0], points[1], points[2], 1.e-6F, displacedFactor));
  TripletFitResult displaced;
  BOOST_REQUIRE(fitTriplet(displacedFactor, displaced));
  BOOST_TEST(exact.chi2 < kExactHelixChi2);
  BOOST_TEST(displaced.chi2 > exact.chi2 + 1.F);
}

// The independent copied ITS geometry and double Gauss-Jordan solve above
// checks a range of noisy helices, anisotropic errors and scattering scales.
BOOST_AUTO_TEST_CASE(AdjacentFitMatchesTheDoubleOracleAcrossScales)
{
  std::mt19937 random{501010U};
  std::uniform_real_distribution<double> unit{0., 1.};
  for (int sample = 0; sample < 300; ++sample) {
    const double curvature = (sample % 2 == 0 ? 1. : -1.) * (0.002 + 0.013 * unit(random));
    const double cotangent = 0.1 + 7. * unit(random);
    double arc = 2. + 5. * unit(random);
    std::array<PairPoint, 4> points{};
    for (auto& point : points) {
      arc += 2. + 6. * unit(random);
      point = helixPoint(curvature, arc, cotangent);
      point.x += static_cast<Scalar>(0.01 * (unit(random) - 0.5));
      point.y += static_cast<Scalar>(0.01 * (unit(random) - 0.5));
      point.z += static_cast<Scalar>(0.03 * (unit(random) - 0.5));
      point.r = std::hypot(point.x, point.y);
      point.varR = static_cast<Scalar>(1.e-7 + 5.e-6 * unit(random));
      point.varRPhi = static_cast<Scalar>(1.e-7 + 8.e-6 * unit(random));
      point.varZ = static_cast<Scalar>(1.e-7 + 2.e-5 * unit(random));
    }
    const auto second = static_cast<Scalar>(std::pow(10., -8. + 4. * unit(random)));
    const auto third = static_cast<Scalar>(std::pow(10., -8. + 4. * unit(random)));
    const auto expected = referenceAdjacentFit(points, second, third);
    TripletFitResult result;
    BOOST_TEST_CONTEXT("sample " << sample)
    {
      BOOST_REQUIRE(fitAdjacentTriplets(points[0], points[1], points[2], points[3], second, third, result));
      checkClose(result.curvature, static_cast<Scalar>(expected.curvature), 3.e-3F, 1.e-7F);
      checkClose(result.varCurvature, static_cast<Scalar>(expected.varCurvature), 3.e-3F);
      checkClose(result.chi2, static_cast<Scalar>(expected.chi2), 3.e-3F, kExactHelixChi2);
    }
  }
}

BOOST_AUTO_TEST_CASE(AdjacentFitOfAnExactHelix)
{
  constexpr double kCurvature = 0.008;
  constexpr double kCotangent = 1.2;
  const std::array<PairPoint, 4> points{helixPoint(kCurvature, 3., kCotangent), helixPoint(kCurvature, 7., kCotangent), helixPoint(kCurvature, 12., kCotangent), helixPoint(kCurvature, 18., kCotangent)};
  TripletFitResult adjacent;
  BOOST_REQUIRE(fitAdjacentTriplets(points[0], points[1], points[2], points[3], 2.e-6F, 2.e-6F, adjacent));
  TripletFitResult first;
  TripletFitResult second;
  TripletFitFactor factor;
  BOOST_REQUIRE(study::makeTripletFitFactor(points[0], points[1], points[2], 2.e-6F, factor));
  BOOST_REQUIRE(fitTriplet(factor, first));
  BOOST_REQUIRE(study::makeTripletFitFactor(points[1], points[2], points[3], 2.e-6F, factor));
  BOOST_REQUIRE(fitTriplet(factor, second));
  const double sineTheta = 1. / std::sqrt(1. + (kCotangent * kCotangent));
  checkClose(adjacent.curvature, static_cast<Scalar>(kCurvature * sineTheta), 3.e-3F);
  BOOST_TEST(adjacent.chi2 < kExactHelixChi2);
  // Two kinks constrain the curvature better than either alone.
  BOOST_TEST(adjacent.varCurvature < std::min(first.varCurvature, second.varCurvature));
}

} // namespace o2::gtrk
