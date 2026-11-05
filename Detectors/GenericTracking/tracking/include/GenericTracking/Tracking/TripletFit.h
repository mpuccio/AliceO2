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

/// \file TripletFit.h
/// \brief The General Triplet Track Fit factor for one GTRK triplet.
///
/// Adapted from O2 ITSMFTTracking/TripletFitting.h and TripletFitting.cxx.
/// The measurement covariance of the kinks needs their gradient along the
/// nine measurement axes of the three points; `kinkGradient` differentiates
/// the kinks by hand, in reverse through the few intermediates they are made
/// of, where ITS stores a nine-coordinate Jacobian.

#ifndef GTRK_TRACKING_TRIPLETFIT_H_
#define GTRK_TRACKING_TRIPLETFIT_H_

#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/Trig.h"
#include "GenericTracking/Tracking/Triplet.h"

#include "GPUCommonMath.h"

#include <cstddef>
#include <limits>

namespace o2::gtrk
{

struct TripletKinkVector {
  Scalar theta = 0.F;
  Scalar phi = 0.F;
};

/// The measurement axes of a triplet: radial, r-phi and z of each of its
/// three points, direction `3 slot + axis`.
inline constexpr std::size_t kTripletAxes = 9;

/// What the fits read of a triplet `a b c` beyond its points' errors: the
/// kinks `psi` at `b` and their slopes `rho` in the curvature, sin(theta) from
/// `a` to `c`, and the gradient of `psi + kappa rho` along the measurement
/// axes, kappa the triplet's reference curvature. The triplet fit and the
/// neighbour fit of two adjacent triplets each make them from the points
/// where a fit needs them, not stored (92 bytes per triplet).
struct KinkTerms {
  TripletKinkVector psi{};
  TripletKinkVector rho{};
  Scalar sineTheta = 0.F;
  Scalar theta[kTripletAxes]{};
  Scalar phi[kTripletAxes]{};
};

/// The kinks `psi` at the middle point, theta and phi, and their slopes `rho`
/// in the curvature.
struct TripletKinks {
  TripletKinkVector psi{};
  TripletKinkVector rho{};
  bool valid = false;
};

struct TripletKinkCovariance {
  Scalar thetaTheta = 0.F;
  Scalar thetaPhi = 0.F;
  Scalar phiPhi = 0.F;
};

struct TripletFitFactor {
  TripletKinkVector psi{};
  TripletKinkVector rho{};
  TripletKinkCovariance covariance{};

  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool isValid() const { return o2::gpu::GPUCommonMath::Abs(rho.phi) > 0.F; }
};

/// The triplet fit's curvature is its own parameter, the 3D curvature: it is
/// proportional to q/p, and the transverse curvature is it over sin(theta).
/// Adjacent triplets share it, which is why it is the one stored.
struct TripletFitResult {
  Scalar curvature = 0.F;
  Scalar varCurvature = 0.F;
  Scalar chi2 = 0.F;
};

/// A triplet's transverse chords, a -> b, b -> c and a -> c: made once per
/// candidate for the curvature test, the kinks and the reference dip,
/// the same expressions each made.
struct TripletChords {
  Scalar length01 = 0.F;
  Scalar length12 = 0.F;
  Scalar length02 = 0.F;
};

[[nodiscard]] GPUhdi() GPUAlwaysInline() TripletChords tripletChords(const PairPoint& first, const PairPoint& middle, const PairPoint& third)
{
  const Scalar dx01 = middle.x - first.x;
  const Scalar dy01 = middle.y - first.y;
  const Scalar dx12 = third.x - middle.x;
  const Scalar dy12 = third.y - middle.y;
  const Scalar dx02 = third.x - first.x;
  const Scalar dy02 = third.y - first.y;
  return {.length01 = o2::gpu::CAMath::Hypot(dx01, dy01), .length12 = o2::gpu::CAMath::Hypot(dx12, dy12), .length02 = o2::gpu::CAMath::Hypot(dx02, dy02)};
}

/// The kinks, the reference curvature `kappa = -psi.phi / rho.phi`, and the
/// gradient of `psi.theta + kappa rho.theta` and `psi.phi + kappa rho.phi`
/// along the measurement axes: what the covariance of the kinks is made of.
struct KinkGradient {
  TripletKinks kinks;
  Scalar kappa = 0.F;
  Scalar theta[kTripletAxes]{};
  Scalar phi[kTripletAxes]{};
};

namespace detail
{

/// An adjoint of the two combinations at once.
struct KinkAdjoint {
  Scalar theta = 0.F;
  Scalar phi = 0.F;
};

[[nodiscard]] GPUhdi() GPUAlwaysInline() KinkAdjoint operator*(const KinkAdjoint& adjoint, Scalar factor)
{
  return {.theta = adjoint.theta * factor, .phi = adjoint.phi * factor};
}

GPUhdi() GPUAlwaysInline() KinkAdjoint& operator+=(KinkAdjoint& adjoint, const KinkAdjoint& other)
{
  adjoint.theta += other.theta;
  adjoint.phi += other.phi;
  return adjoint;
}

GPUhdi() GPUAlwaysInline() KinkAdjoint& operator-=(KinkAdjoint& adjoint, const KinkAdjoint& other)
{
  adjoint.theta -= other.theta;
  adjoint.phi -= other.phi;
  return adjoint;
}

/// One segment of the helix between two points: its chord, delta z and the
/// values the kinks read, and the slopes of its half angle, arc factor and
/// angle cotangent in the half sine, for the reverse pass.
struct SegmentForward {
  Scalar chord = 0.F;
  Scalar deltaZ = 0.F;
  Scalar halfAngle = 0.F;
  Scalar asinOverArgument = 0.F;
  Scalar angleCotangent = 0.F;
  Scalar bendingAngle = 0.F;
  Scalar arc = 0.F;
  Scalar cotangent = 0.F;
  Scalar sine = 0.F;
  Scalar cosine = 0.F;
  Scalar index = 0.F;
  Scalar theta = 0.F;
  Scalar halfAngleSlope = 0.F;        ///< d halfAngle / d halfSine
  Scalar asinOverArgumentSlope = 0.F; ///< d asinOverArgument / d halfSine
  Scalar angleCotangentSlope = 0.F;   ///< d angleCotangent / d halfSine
  bool valid = false;
};

[[nodiscard]] GPUhdi() SegmentForward segmentForward(Scalar curvature, Scalar chord, Scalar deltaZ)
{
  SegmentForward segment{.chord = chord, .deltaZ = deltaZ};
  const Scalar halfSine = 0.5F * curvature * chord;
  if (o2::gpu::GPUCommonMath::Abs(halfSine) >= 1.F) {
    return segment;
  }
  const Scalar halfSine2 = halfSine * halfSine;
  const Scalar halfSine4 = halfSine2 * halfSine2;
  const Scalar halfAngle = asinOf(halfSine);
  const Scalar root = o2::gpu::GPUCommonMath::Sqrt(1.F - halfSine2);
  if (o2::gpu::GPUCommonMath::Abs(halfSine) < 1.e-4F) {
    segment.asinOverArgument = 1.F + halfSine2 * (1.F / 6.F) + halfSine4 * (3.F / 40.F);
    segment.angleCotangent = 1.F - halfSine2 * (1.F / 3.F) - halfSine4 * (2.F / 15.F);
    segment.asinOverArgumentSlope = halfSine * (1.F / 3.F) + halfSine2 * halfSine * (3.F / 10.F);
    segment.angleCotangentSlope = -(halfSine * (2.F / 3.F)) - halfSine2 * halfSine * (8.F / 15.F);
  } else {
    segment.asinOverArgument = halfAngle / halfSine;
    segment.angleCotangent = halfAngle * root / halfSine;
    segment.asinOverArgumentSlope = (1.F / root - segment.asinOverArgument) / halfSine;
    segment.angleCotangentSlope = (1.F - segment.angleCotangent) / halfSine - halfAngle / root;
  }
  segment.halfAngle = halfAngle;
  segment.halfAngleSlope = 1.F / root;
  segment.bendingAngle = 2.F * halfAngle;
  segment.arc = chord * segment.asinOverArgument;
  segment.cotangent = deltaZ / segment.arc;
  segment.sine = 1.F / o2::gpu::GPUCommonMath::Sqrt(1.F + segment.cotangent * segment.cotangent);
  segment.cosine = segment.cotangent * segment.sine;
  segment.index = 1.F / (segment.angleCotangent * segment.sine * segment.sine + segment.cosine * segment.cosine);
  segment.theta = atan2Of(segment.arc, deltaZ);
  segment.valid = segment.arc > 0.F && segment.sine > 0.F && segment.index > 0.F;
  return segment;
}

/// The adjoints a segment's outputs received.
struct SegmentAdjoint {
  KinkAdjoint bendingAngle;
  KinkAdjoint arc;
  KinkAdjoint cotangent;
  KinkAdjoint sine;
  KinkAdjoint cosine;
  KinkAdjoint index;
  KinkAdjoint theta;
};

/// A segment's adjoints carried back to its curvature, chord and delta z.
GPUhdi() void segmentBackward(const SegmentForward& s, SegmentAdjoint adjoint, Scalar curvature, KinkAdjoint& curvatureAdjoint, KinkAdjoint& chordAdjoint,
                              KinkAdjoint& deltaZAdjoint)
{
  // theta = atan2(arc, deltaZ)
  const Scalar radius2 = s.arc * s.arc + s.deltaZ * s.deltaZ;
  adjoint.arc += adjoint.theta * (s.deltaZ / radius2);
  deltaZAdjoint -= adjoint.theta * (s.arc / radius2);
  // index = 1 / (angleCotangent sine^2 + cosine^2)
  const KinkAdjoint denominator = adjoint.index * (-(s.index * s.index));
  const KinkAdjoint angleCotangent = denominator * (s.sine * s.sine);
  adjoint.sine += denominator * (2.F * s.angleCotangent * s.sine);
  adjoint.cosine += denominator * (2.F * s.cosine);
  // cosine = cotangent sine
  adjoint.cotangent += adjoint.cosine * s.sine;
  adjoint.sine += adjoint.cosine * s.cotangent;
  // sine = (1 + cotangent^2)^(-1/2)
  adjoint.cotangent -= adjoint.sine * (s.cotangent * s.sine * s.sine * s.sine);
  // cotangent = deltaZ / arc
  deltaZAdjoint += adjoint.cotangent * (1.F / s.arc);
  adjoint.arc -= adjoint.cotangent * (s.cotangent / s.arc);
  // arc = chord asinOverArgument
  chordAdjoint += adjoint.arc * s.asinOverArgument;
  const KinkAdjoint asinOverArgument = adjoint.arc * s.chord;
  // bendingAngle = 2 halfAngle, and the three functions of the half sine
  KinkAdjoint halfSine = adjoint.bendingAngle * (2.F * s.halfAngleSlope);
  halfSine += asinOverArgument * s.asinOverArgumentSlope;
  halfSine += angleCotangent * s.angleCotangentSlope;
  // halfSine = curvature chord / 2
  curvatureAdjoint += halfSine * (0.5F * s.chord);
  chordAdjoint += halfSine * (0.5F * curvature);
}

/// A chord's adjoint carried to its two coordinate differences.
GPUhdi() GPUAlwaysInline() void chordBackward(const KinkAdjoint& chordAdjoint, Scalar dx, Scalar dy, Scalar chord, KinkAdjoint& dxAdjoint,
                                              KinkAdjoint& dyAdjoint)
{
  dxAdjoint += chordAdjoint * (dx / chord);
  dyAdjoint += chordAdjoint * (dy / chord);
}

} // namespace detail

/// The kinks and their gradient: the forward pass computes the kinks, the
/// reverse pass carries the adjoints of `psi + kappa rho` back through the
/// chords, the curvature and each segment's half sine, arc, cotangent, sine,
/// cosine and index to the three points, and projects them on each point's
/// radial, r-phi and z axes. Invalid where the kinks are not defined.
/// \a chords are the points' (`tripletChords`).
[[nodiscard]] GPUhdi() GPUAlwaysInline() KinkGradient kinkGradient(const TripletChords& chords, const PairPoint& first, const PairPoint& middle,
                                                                   const PairPoint& third)
{
  KinkGradient result;
  const Scalar dx01 = middle.x - first.x;
  const Scalar dy01 = middle.y - first.y;
  const Scalar dz01 = middle.z - first.z;
  const Scalar dx12 = third.x - middle.x;
  const Scalar dy12 = third.y - middle.y;
  const Scalar dz12 = third.z - middle.z;
  const Scalar dx02 = third.x - first.x;
  const Scalar dy02 = third.y - first.y;
  const Scalar length01 = chords.length01;
  const Scalar length12 = chords.length12;
  const Scalar length02 = chords.length02;
  if (!(length01 > 0.F) || !(length12 > 0.F) || !(length02 > 0.F)) {
    return result;
  }
  const Scalar cross = dx01 * dy12 - dy01 * dx12;
  const Scalar lengths = length01 * length12 * length02;
  const Scalar curvature = 2.F * cross / lengths;
  const auto s1 = detail::segmentForward(curvature, length01, dz01);
  const auto s2 = detail::segmentForward(curvature, length12, dz12);
  if (!s1.valid || !s2.valid) {
    return result;
  }

  const Scalar psiPhi = 0.5F * (s1.bendingAngle * s1.index + s2.bendingAngle * s2.index);
  const Scalar psiTheta = s2.theta - s1.theta + (1.F - s2.index) * s2.cotangent - (1.F - s1.index) * s1.cotangent;
  const Scalar rhoPhi = -0.5F * (s1.arc * s1.index / s1.sine + s2.arc * s2.index / s2.sine);
  const Scalar maximumLength = length01 > length12 ? length01 : length12;
  // Below this the slope's closed form cancels; its expansion to first order.
  const bool straight = 0.5F * o2::gpu::GPUCommonMath::Abs(curvature) * maximumLength < 1.e-4F;
  Scalar rhoTheta = 0.F;
  if (straight) {
    rhoTheta = curvature * (length12 * length12 * s2.cosine - length01 * length01 * s1.cosine) / 12.F;
  } else {
    rhoTheta = ((1.F - s1.index) * s1.cotangent / s1.sine - (1.F - s2.index) * s2.cotangent / s2.sine) / curvature;
  }
  if (!(o2::gpu::GPUCommonMath::Abs(rhoPhi) > 0.F)) {
    return result;
  }
  result.kinks = {.psi = {.theta = psiTheta, .phi = psiPhi}, .rho = {.theta = rhoTheta, .phi = rhoPhi}, .valid = true};
  const Scalar kappa = -psiPhi / rhoPhi;
  result.kappa = kappa;

  // The adjoints of psi.theta + kappa rho.theta and psi.phi + kappa rho.phi.
  using Adjoint = detail::KinkAdjoint;
  const Adjoint psiThetaAdjoint{.theta = 1.F, .phi = 0.F};
  const Adjoint psiPhiAdjoint{.theta = 0.F, .phi = 1.F};
  const Adjoint rhoThetaAdjoint{.theta = kappa, .phi = 0.F};
  const Adjoint rhoPhiAdjoint{.theta = 0.F, .phi = kappa};
  detail::SegmentAdjoint a1;
  detail::SegmentAdjoint a2;
  Adjoint curvatureAdjoint;
  Adjoint length01Adjoint;
  Adjoint length12Adjoint;

  // rhoPhi = -(arc index / sine)_1 / 2 - (arc index / sine)_2 / 2
  const Adjoint half = rhoPhiAdjoint * -0.5F;
  a1.arc += half * (s1.index / s1.sine);
  a1.index += half * (s1.arc / s1.sine);
  a1.sine -= half * (s1.arc * s1.index / (s1.sine * s1.sine));
  a2.arc += half * (s2.index / s2.sine);
  a2.index += half * (s2.arc / s2.sine);
  a2.sine -= half * (s2.arc * s2.index / (s2.sine * s2.sine));
  // psiPhi = (bendingAngle index)_1 / 2 + (bendingAngle index)_2 / 2
  a1.bendingAngle += psiPhiAdjoint * (0.5F * s1.index);
  a1.index += psiPhiAdjoint * (0.5F * s1.bendingAngle);
  a2.bendingAngle += psiPhiAdjoint * (0.5F * s2.index);
  a2.index += psiPhiAdjoint * (0.5F * s2.bendingAngle);
  // psiTheta = theta_2 - theta_1 + (1 - index_2) cotangent_2 - (1 - index_1) cotangent_1
  a2.theta += psiThetaAdjoint;
  a1.theta -= psiThetaAdjoint;
  a2.index -= psiThetaAdjoint * s2.cotangent;
  a2.cotangent += psiThetaAdjoint * (1.F - s2.index);
  a1.index += psiThetaAdjoint * s1.cotangent;
  a1.cotangent -= psiThetaAdjoint * (1.F - s1.index);
  if (straight) {
    // rhoTheta = curvature (length12^2 cosine_2 - length01^2 cosine_1) / 12
    const Scalar difference = length12 * length12 * s2.cosine - length01 * length01 * s1.cosine;
    curvatureAdjoint += rhoThetaAdjoint * (difference / 12.F);
    const Adjoint differenceAdjoint = rhoThetaAdjoint * (curvature / 12.F);
    length12Adjoint += differenceAdjoint * (2.F * length12 * s2.cosine);
    a2.cosine += differenceAdjoint * (length12 * length12);
    length01Adjoint -= differenceAdjoint * (2.F * length01 * s1.cosine);
    a1.cosine -= differenceAdjoint * (length01 * length01);
  } else {
    // rhoTheta = ((1 - index_1) cotangent_1 / sine_1 - (1 - index_2) cotangent_2 / sine_2) / curvature
    const Adjoint numerator = rhoThetaAdjoint * (1.F / curvature);
    curvatureAdjoint -= rhoThetaAdjoint * (rhoTheta / curvature);
    a1.index -= numerator * (s1.cotangent / s1.sine);
    a1.cotangent += numerator * ((1.F - s1.index) / s1.sine);
    a1.sine -= numerator * ((1.F - s1.index) * s1.cotangent / (s1.sine * s1.sine));
    a2.index += numerator * (s2.cotangent / s2.sine);
    a2.cotangent -= numerator * ((1.F - s2.index) / s2.sine);
    a2.sine += numerator * ((1.F - s2.index) * s2.cotangent / (s2.sine * s2.sine));
  }

  Adjoint dz01Adjoint;
  Adjoint dz12Adjoint;
  detail::segmentBackward(s1, a1, curvature, curvatureAdjoint, length01Adjoint, dz01Adjoint);
  detail::segmentBackward(s2, a2, curvature, curvatureAdjoint, length12Adjoint, dz12Adjoint);

  // curvature = 2 cross / (length01 length12 length02)
  const Adjoint crossAdjoint = curvatureAdjoint * (2.F / lengths);
  const Adjoint lengthsAdjoint = curvatureAdjoint * (-curvature / lengths);
  length01Adjoint += lengthsAdjoint * (length12 * length02);
  length12Adjoint += lengthsAdjoint * (length01 * length02);
  const Adjoint length02Adjoint = lengthsAdjoint * (length01 * length12);
  // cross = dx01 dy12 - dy01 dx12
  Adjoint dx01Adjoint = crossAdjoint * dy12;
  Adjoint dy12Adjoint = crossAdjoint * dx01;
  Adjoint dy01Adjoint = crossAdjoint * (-dx12);
  Adjoint dx12Adjoint = crossAdjoint * (-dy01);
  Adjoint dx02Adjoint;
  Adjoint dy02Adjoint;
  detail::chordBackward(length01Adjoint, dx01, dy01, length01, dx01Adjoint, dy01Adjoint);
  detail::chordBackward(length12Adjoint, dx12, dy12, length12, dx12Adjoint, dy12Adjoint);
  detail::chordBackward(length02Adjoint, dx02, dy02, length02, dx02Adjoint, dy02Adjoint);

  // The differences to the points, the points to their measurement axes.
  Adjoint xAdjoint[3];
  Adjoint yAdjoint[3];
  Adjoint zAdjoint[3];
  xAdjoint[0] -= dx01Adjoint;
  xAdjoint[0] -= dx02Adjoint;
  xAdjoint[1] += dx01Adjoint;
  xAdjoint[1] -= dx12Adjoint;
  xAdjoint[2] += dx12Adjoint;
  xAdjoint[2] += dx02Adjoint;
  yAdjoint[0] -= dy01Adjoint;
  yAdjoint[0] -= dy02Adjoint;
  yAdjoint[1] += dy01Adjoint;
  yAdjoint[1] -= dy12Adjoint;
  yAdjoint[2] += dy12Adjoint;
  yAdjoint[2] += dy02Adjoint;
  zAdjoint[0] -= dz01Adjoint;
  zAdjoint[1] += dz01Adjoint;
  zAdjoint[1] -= dz12Adjoint;
  zAdjoint[2] += dz12Adjoint;
  const PairPoint* points[3]{&first, &middle, &third};
  for (int slot = 0; slot < 3; ++slot) {
    const Scalar inverseR = 1.F / points[slot]->r;
    const Scalar radialX = points[slot]->x * inverseR;
    const Scalar radialY = points[slot]->y * inverseR;
    const int radial = 3 * slot;
    result.theta[radial] = xAdjoint[slot].theta * radialX + yAdjoint[slot].theta * radialY;
    result.phi[radial] = xAdjoint[slot].phi * radialX + yAdjoint[slot].phi * radialY;
    result.theta[radial + 1] = yAdjoint[slot].theta * radialX - xAdjoint[slot].theta * radialY;
    result.phi[radial + 1] = yAdjoint[slot].phi * radialX - xAdjoint[slot].phi * radialY;
    result.theta[radial + 2] = zAdjoint[slot].theta;
    result.phi[radial + 2] = zAdjoint[slot].phi;
  }
  return result;
}

/// The same, the chords made here.
[[nodiscard]] GPUhdi() GPUAlwaysInline() KinkGradient kinkGradient(const PairPoint& first, const PairPoint& middle, const PairPoint& third)
{
  return kinkGradient(tripletChords(first, middle, third), first, middle, third);
}

namespace detail
{

[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar axisVariance(const PairPoint& point, int axis)
{
  if (axis == 0) {
    return point.varR;
  }
  if (axis == 1) {
    return point.varRPhi;
  }
  return point.varZ;
}

/// The chord a -> c's sin(theta), its transverse length \a transverse.
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar referenceSinTheta(Scalar transverse, const PairPoint& first, const PairPoint& third)
{
  const Scalar dz = third.z - first.z;
  return transverse / o2::gpu::CAMath::Hypot(transverse, dz);
}

[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar referenceSinTheta(const PairPoint& first, const PairPoint& third)
{
  const Scalar dx = third.x - first.x;
  const Scalar dy = third.y - first.y;
  return referenceSinTheta(o2::gpu::CAMath::Hypot(dx, dy), first, third);
}

GPUhdi() GPUAlwaysInline() void addCovariance(TripletKinkCovariance& covariance, Scalar variance, Scalar thetaDerivative, Scalar phiDerivative)
{
  covariance.thetaTheta += variance * thetaDerivative * thetaDerivative;
  covariance.thetaPhi += variance * thetaDerivative * phiDerivative;
  covariance.phiPhi += variance * phiDerivative * phiDerivative;
}

/// The curvature that minimises `(psi + kappa rho)' V^-1 (psi + kappa rho)`,
/// its variance and the chi2 left, from the three contractions with `V^-1`.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool curvatureFit(Scalar psiPrecisionPsi, Scalar rhoPrecisionPsi, Scalar rhoPrecisionRho, TripletFitResult& result)
{
  if (!(rhoPrecisionRho > 0.F)) {
    return false;
  }
  const Scalar curvature = -rhoPrecisionPsi / rhoPrecisionRho;
  const Scalar varCurvature = 1.F / rhoPrecisionRho;
  const Scalar removedCurvatureTerm = rhoPrecisionPsi * rhoPrecisionPsi / rhoPrecisionRho;
  Scalar chi2 = psiPrecisionPsi - removedCurvatureTerm;
  const Scalar scale = o2::gpu::GPUCommonMath::Max(o2::gpu::GPUCommonMath::Abs(psiPrecisionPsi), o2::gpu::GPUCommonMath::Abs(removedCurvatureTerm));
  const Scalar tolerance = 128.F * std::numeric_limits<Scalar>::epsilon() * scale;
  if (chi2 < 0.F && chi2 >= -tolerance) {
    chi2 = 0.F;
  }
  if (!(varCurvature > 0.F) || chi2 < 0.F) {
    return false;
  }
  result = {.curvature = curvature, .varCurvature = varCurvature, .chi2 = chi2};
  return true;
}

} // namespace detail

/// The space-angle scattering variance at the middle layer, at pT_min and for
/// a track at polar angle theta: `k` times the path factor over `p^2`.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar middleScattering(const ScatteringTable<Det>& table, const PairCut& cut, LayerID middleLayer, Scalar sineTheta)
{
  if (!(sineTheta > 0.F)) {
    return 0.F;
  }
  const Scalar material = table.layers[static_cast<std::size_t>(middleLayer)].k;
  if (Detector<Det>::isBarrelLayer(middleLayer)) {
    return material * sineTheta * cut.invBetaMomentum2;
  }
  const Scalar cosineTheta = o2::gpu::GPUCommonMath::Sqrt(o2::gpu::GPUCommonMath::Max(0.F, 1.F - sineTheta * sineTheta));
  return cosineTheta > 0.F ? material * sineTheta * sineTheta * cut.invBetaMomentum2 / cosineTheta : 0.F;
}

/// A passive scatterer's space-angle variance at pT_min, for a track at polar
/// angle theta: as `middleScattering` takes a layer's.
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar passiveScattering(const PassiveScatterer& passive, const PairCut& cut, Scalar sineTheta)
{
  if (!(sineTheta > 0.F)) {
    return 0.F;
  }
  if (passive.barrel) {
    return passive.at.k * sineTheta * cut.invBetaMomentum2;
  }
  const Scalar cosineTheta = o2::gpu::GPUCommonMath::Sqrt(o2::gpu::GPUCommonMath::Max(0.F, 1.F - sineTheta * sineTheta));
  return cosineTheta > 0.F ? passive.at.k * sineTheta * sineTheta * cut.invBetaMomentum2 / cosineTheta : 0.F;
}

/// Calls \a visit(variance, t) for each passive scatterer the straight segment
/// from \a from to \a to crosses, with its `passiveScattering` and the
/// fraction `t` of the way it is crossed at: in radius for a cylinder, in z
/// for a disc.
///
/// \a passives: the ones to try, a bit each, in index order; every passive by
/// default. The table's `into` and `outOf` of the segment's layers hold every
/// passive it can cross: the same visits.
template <DetectorLayout Det, typename Visit>
GPUhdi() GPUAlwaysInline() void visitPassivesBetween(const ScatteringTable<Det>& table, const PairCut& cut, const PairPoint& from, const PairPoint& to,
                                                     Scalar sineTheta, Visit visit, std::uint64_t passives = ScatteringTable<Det>::kAllPassives)
{
  if constexpr (ScatteringTable<Det>::kPassiveCount > 0) {
    for (std::uint64_t left = passives & ScatteringTable<Det>::kAllPassives; left != 0U; left &= left - 1U) {
      const auto& passive = table.passives[o2::gpu::CAMath::Ctz64(left)];
      const Scalar begin = passive.barrel ? from.r : from.z;
      const Scalar end = passive.barrel ? to.r : to.z;
      const Scalar span = end - begin;
      if (!(o2::gpu::GPUCommonMath::Abs(span) > 0.F)) {
        continue;
      }
      const Scalar t = (passive.at.position - begin) / span;
      if (t <= 0.F || t >= 1.F) {
        continue;
      }
      const Scalar along = passive.barrel ? from.z + (t * (to.z - from.z)) : from.r + (t * (to.r - from.r));
      if (along < passive.at.low || along > passive.at.high) {
        continue;
      }
      visit(passiveScattering(passive, cut, sineTheta), t);
    }
  }
}

/// The variance of the kink at a triplet's middle point: the middle layer's
/// scattering and, as ALICE 3's pipe between its layers needs, every passive
/// scatterer between the points. A kink `theta_p` a fraction `t` of the way
/// from `a` to `b` turns the angle between the chords at `b` by `theta_p t`;
/// one from `b` to `c` by `theta_p (1 - t)`. \a middleScale multiplies the
/// middle layer's variance alone: the kink allowance of a kinked triplet.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar tripletScattering(const ScatteringTable<Det>& table, const PairCut& cut, LayerID middleLayer,
                                                                  const PairPoint& first, const PairPoint& middle, const PairPoint& third, Scalar sineTheta,
                                                                  Scalar middleScale = 1.F)
{
  Scalar variance = middleScattering(table, cut, middleLayer, sineTheta) * middleScale;
  visitPassivesBetween(table, cut, first, middle, sineTheta, [&](Scalar passive, Scalar t) { variance += passive * t * t; }, table.into(middleLayer));
  visitPassivesBetween(
    table, cut, middle, third, sineTheta, [&](Scalar passive, Scalar t) { variance += passive * (1.F - t) * (1.F - t); }, table.outOf(middleLayer));
  return variance;
}

/// A triplet's kink terms from its \a chords (`tripletChords`); false where
/// its kinks are not defined.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool makeKinkTerms(const TripletChords& chords, const PairPoint& first, const PairPoint& middle,
                                                            const PairPoint& third, KinkTerms& terms)
{
  const auto gradient = kinkGradient(chords, first, middle, third);
  const Scalar sineTheta = detail::referenceSinTheta(chords.length02, first, third);
  if (!gradient.kinks.valid || !(sineTheta > 0.F)) {
    return false;
  }
  terms.psi = gradient.kinks.psi;
  terms.rho = gradient.kinks.rho;
  terms.sineTheta = sineTheta;
  for (std::size_t direction = 0; direction < kTripletAxes; ++direction) {
    terms.theta[direction] = gradient.theta[direction];
    terms.phi[direction] = gradient.phi[direction];
  }
  return true;
}

/// The same, the chords made here.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool makeKinkTerms(const PairPoint& first, const PairPoint& middle, const PairPoint& third, KinkTerms& terms)
{
  return makeKinkTerms(tripletChords(first, middle, third), first, middle, third, terms);
}

/// The fit factor of a triplet from its kink terms and its points' errors.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool makeTripletFitFactor(const KinkTerms& terms, const PairPoint& first, const PairPoint& middle,
                                                                   const PairPoint& third, Scalar angularVariance, TripletFitFactor& result)
{
  const Scalar sineTheta = terms.sineTheta;
  TripletFitFactor scratch{
    .psi = terms.psi, .rho = terms.rho, .covariance = {.thetaTheta = angularVariance, .thetaPhi = 0.F, .phiPhi = angularVariance / (sineTheta * sineTheta)}};
  for (int hit = 0; hit < 3; ++hit) {
    const PairPoint* point = &third;
    if (hit == 0) {
      point = &first;
    } else if (hit == 1) {
      point = &middle;
    }
    for (int axis = 0; axis < 3; ++axis) {
      const int direction = (3 * hit) + axis;
      detail::addCovariance(scratch.covariance, detail::axisVariance(*point, axis), terms.theta[direction], terms.phi[direction]);
    }
  }
  if (!scratch.isValid()) {
    return false;
  }
  result = scratch;
  return true;
}

[[nodiscard]] GPUhdi() GPUAlwaysInline() bool fitTriplet(const TripletFitFactor& factor, TripletFitResult& result)
{
  const Scalar determinant = factor.covariance.thetaTheta * factor.covariance.phiPhi - factor.covariance.thetaPhi * factor.covariance.thetaPhi;
  if (!(determinant > 0.F)) {
    return false;
  }
  const Scalar inverseThetaTheta = factor.covariance.phiPhi / determinant;
  const Scalar inverseThetaPhi = -factor.covariance.thetaPhi / determinant;
  const Scalar inversePhiPhi = factor.covariance.thetaTheta / determinant;
  const Scalar precisionPsiTheta = inverseThetaTheta * factor.psi.theta + inverseThetaPhi * factor.psi.phi;
  const Scalar precisionPsiPhi = inverseThetaPhi * factor.psi.theta + inversePhiPhi * factor.psi.phi;
  const Scalar precisionRhoTheta = inverseThetaTheta * factor.rho.theta + inverseThetaPhi * factor.rho.phi;
  const Scalar precisionRhoPhi = inverseThetaPhi * factor.rho.theta + inversePhiPhi * factor.rho.phi;
  const Scalar rhoPrecisionPsi = factor.rho.theta * precisionPsiTheta + factor.rho.phi * precisionPsiPhi;
  const Scalar rhoPrecisionRho = factor.rho.theta * precisionRhoTheta + factor.rho.phi * precisionRhoPhi;
  const Scalar psiPrecisionPsi = factor.psi.theta * precisionPsiTheta + factor.psi.phi * precisionPsiPhi;
  return detail::curvatureFit(psiPrecisionPsi, rhoPrecisionPsi, rhoPrecisionRho, result);
}

/// The triplet fit from its kink \a terms (`makeKinkTerms`), with the
/// scattering at the middle layer at pT_min from the pair cut's table, times
/// \a middleScale.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool fitTriplet(const KinkTerms& terms, const PairPoint& first, const PairPoint& middle, const PairPoint& third,
                                                         LayerID middleLayer, const ScatteringTable<Det>& table, const PairCut& cut, TripletFitResult& result,
                                                         Scalar middleScale = 1.F)
{
  TripletFitFactor factor;
  return makeTripletFitFactor(terms, first, middle, third, tripletScattering(table, cut, middleLayer, first, middle, third, terms.sineTheta, middleScale),
                              factor) &&
         fitTriplet(factor, result);
}

/// The same with the kink terms made here; \a terms receives them.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool fitTriplet(const PairPoint& first, const PairPoint& middle, const PairPoint& third, LayerID middleLayer,
                                                         const ScatteringTable<Det>& table, const PairCut& cut, TripletFitResult& result, KinkTerms& terms,
                                                         Scalar middleScale = 1.F)
{
  return makeKinkTerms(first, middle, third, terms) && fitTriplet(terms, first, middle, third, middleLayer, table, cut, result, middleScale);
}

namespace detail
{

constexpr int kAdjacentKinks = 4;

/// Factor the symmetric positive definite \a covariance, given by its lower
/// triangle, in place: the lower triangle becomes the Cholesky factor.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool choleskyDecompose(Scalar (&covariance)[kAdjacentKinks][kAdjacentKinks])
{
  for (int row = 0; row < kAdjacentKinks; ++row) {
    for (int column = 0; column <= row; ++column) {
      Scalar value = covariance[row][column];
      for (int k = 0; k < column; ++k) {
        value -= covariance[row][k] * covariance[column][k];
      }
      if (row == column) {
        if (!(value > 0.F)) {
          return false;
        }
        covariance[row][row] = o2::gpu::GPUCommonMath::Sqrt(value);
      } else {
        covariance[row][column] = value / covariance[column][column];
      }
    }
  }
  return true;
}

/// Whiten a vector with L^-1. The curvature fit needs only contractions with
/// C^-1 = L^-T L^-1, so backward substitution is unnecessary.
GPUhdi() GPUAlwaysInline() void choleskyWhiten(const Scalar (&lower)[kAdjacentKinks][kAdjacentKinks], const Scalar (&right)[kAdjacentKinks],
                                               Scalar (&whitened)[kAdjacentKinks])
{
  for (int row = 0; row < kAdjacentKinks; ++row) {
    Scalar value = right[row];
    for (int column = 0; column < row; ++column) {
      value -= lower[row][column] * whitened[column];
    }
    whitened[row] = value / lower[row][row];
  }
}

[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar dot(const Scalar (&left)[kAdjacentKinks], const Scalar (&right)[kAdjacentKinks])
{
  Scalar result = 0.F;
  for (int i = 0; i < kAdjacentKinks; ++i) {
    result += left[i] * right[i];
  }
  return result;
}

} // namespace detail

/// Adjacent triplets `a b c` and `b c d` fitted with one curvature:
/// ITSMFTTracking's `fitAdjacentTripletFactors`. The four kinks -- theta and phi at `b` and at
/// `c` -- are correlated through the shared measurements `b` and `c`; their
/// covariance is the scattering at `b` and `c`, \a varianceSecond and \a
/// varianceThird as space-angle variances, plus the four measurements' errors
/// carried by each triplet's kink gradient, its slope taken at its own
/// reference curvature as in makeTripletFitFactor. The chi2 has three
/// degrees of freedom. \a first and \a second are the kink terms of `a b c`
/// and `b c d`. \a covarianceSecondThird is the space-angle covariance of
/// the kinks at `b` and `c` that one scattering between them gives both
/// (`adjacentScattering`).
[[nodiscard]] GPUhdi() bool fitAdjacentTriplets(const KinkTerms& first, const KinkTerms& second, const PairPoint& a, const PairPoint& b, const PairPoint& c,
                                                const PairPoint& d, Scalar varianceSecond, Scalar varianceThird, TripletFitResult& result,
                                                Scalar covarianceSecondThird = 0.F)
{
  const Scalar sineFirst = first.sineTheta;
  const Scalar sineSecond = second.sineTheta;
  const Scalar psi[detail::kAdjacentKinks]{first.psi.theta, first.psi.phi, second.psi.theta, second.psi.phi};
  const Scalar rho[detail::kAdjacentKinks]{first.rho.theta, first.rho.phi, second.rho.theta, second.rho.phi};

  Scalar covariance[detail::kAdjacentKinks][detail::kAdjacentKinks]{};
  covariance[0][0] = varianceSecond;
  covariance[1][1] = varianceSecond / (sineFirst * sineFirst);
  covariance[2][2] = varianceThird;
  covariance[3][3] = varianceThird / (sineSecond * sineSecond);
  covariance[2][0] = covarianceSecondThird;
  covariance[3][1] = covarianceSecondThird / (sineFirst * sineSecond);
  for (int hit = 0; hit < 4; ++hit) {
    const PairPoint* point = &d;
    if (hit == 0) {
      point = &a;
    } else if (hit == 1) {
      point = &b;
    } else if (hit == 2) {
      point = &c;
    }
    for (int axis = 0; axis < 3; ++axis) {
      Scalar gradient[detail::kAdjacentKinks]{};
      if (hit < 3) {
        const int direction = (3 * hit) + axis;
        gradient[0] = first.theta[direction];
        gradient[1] = first.phi[direction];
      }
      if (hit > 0) {
        const int direction = (3 * (hit - 1)) + axis;
        gradient[2] = second.theta[direction];
        gradient[3] = second.phi[direction];
      }
      const Scalar variance = detail::axisVariance(*point, axis);
      for (int row = 0; row < detail::kAdjacentKinks; ++row) {
        for (int column = 0; column <= row; ++column) {
          covariance[row][column] += variance * gradient[row] * gradient[column];
        }
      }
    }
  }

  if (!detail::choleskyDecompose(covariance)) {
    return false;
  }
  Scalar whitenedPsi[detail::kAdjacentKinks]{};
  Scalar whitenedRho[detail::kAdjacentKinks]{};
  detail::choleskyWhiten(covariance, psi, whitenedPsi);
  detail::choleskyWhiten(covariance, rho, whitenedRho);
  return detail::curvatureFit(detail::dot(whitenedPsi, whitenedPsi), detail::dot(whitenedRho, whitenedPsi), detail::dot(whitenedRho, whitenedRho), result);
}

/// The kink variances at `b` and `c` of adjacent triplets `a b c` and `b c d`
/// and their covariance: each middle layer's scattering, at the polar angle
/// of its own triplet, and the passive scatterers between the points. One
/// between `b` and `c`, a fraction `t` of the way, turns the kink at `b` by
/// `theta_p (1 - t)` and the one at `c` by `theta_p t`: both, together.
struct AdjacentScattering {
  Scalar second = 0.F;
  Scalar third = 0.F;
  Scalar covariance = 0.F;
};

/// \a secondScale and \a thirdScale multiply the layers' own scattering at `b`
/// and `c`: the kink allowance of a kinked triplet (1 otherwise).
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() AdjacentScattering adjacentScattering(const ScatteringTable<Det>& table, const PairCut& cut, LayerID secondLayer,
                                                                               LayerID thirdLayer, const PairPoint& a, const PairPoint& b, const PairPoint& c,
                                                                               const PairPoint& d, Scalar sineFirst, Scalar sineSecond,
                                                                               Scalar secondScale = 1.F, Scalar thirdScale = 1.F)
{
  AdjacentScattering result{.second = middleScattering(table, cut, secondLayer, sineFirst) * secondScale,
                            .third = middleScattering(table, cut, thirdLayer, sineSecond) * thirdScale};
  visitPassivesBetween(table, cut, a, b, sineFirst, [&](Scalar passive, Scalar t) { result.second += passive * t * t; }, table.into(secondLayer));
  visitPassivesBetween(
    table, cut, b, c, sineFirst,
    [&](Scalar passive, Scalar t) {
      const Scalar atSecond = passive * (1.F - t) * (1.F - t);
      result.second += atSecond;
      // The same kink seen at the second triplet's polar angle.
      const Scalar scale = sineFirst > 0.F ? sineSecond / sineFirst : 1.F;
      result.third += passive * scale * t * t;
      result.covariance += passive * o2::gpu::GPUCommonMath::Sqrt(scale) * t * (1.F - t);
    },
    table.outOf(secondLayer) & table.into(thirdLayer));
  visitPassivesBetween(
    table, cut, c, d, sineSecond, [&](Scalar passive, Scalar t) { result.third += passive * (1.F - t) * (1.F - t); }, table.outOf(thirdLayer));
  return result;
}

/// The same, with the scattering at `b` and `c` at pT_min from the pair
/// cut's table, each at the polar angle of its own triplet.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool fitAdjacentTriplets(const KinkTerms& first, const KinkTerms& second, const PairPoint& a, const PairPoint& b,
                                                                  const PairPoint& c, const PairPoint& d, LayerID secondLayer, LayerID thirdLayer,
                                                                  const ScatteringTable<Det>& table, const PairCut& cut, TripletFitResult& result,
                                                                  Scalar secondScale = 1.F, Scalar thirdScale = 1.F)
{
  const auto scattering = adjacentScattering(table, cut, secondLayer, thirdLayer, a, b, c, d, first.sineTheta, second.sineTheta, secondScale, thirdScale);
  return fitAdjacentTriplets(first, second, a, b, c, d, scattering.second, scattering.third, result, scattering.covariance);
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRIPLETFIT_H_
