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

/// \file DisplacedPair.h
/// \brief Whether a pair of measurements can belong to a track that passes
///        within a distance of the beam line: the pair test of a displaced track, its geometry
///        alone, without the measurements' errors.
///
/// The circles through two points in the transverse plane are one family,
/// by their signed curvature `k` -- positive turning counterclockwise on the
/// way from the source to the target, the straight line at 0, the circle on
/// the chord as diameter at `+-1/h`, `h` the half chord. A track above pT_min
/// is one with `|k| <= k_max`; travelling outward, one whose distance to the
/// beam line grows from the source to the target. On that interval, the
/// signed distance of the circle's closest approach to the beam line (its
/// impact, `circleImpact`) is monotone in `k` (its stationary points are the
/// outward boundaries), the outward set is one interval, and the closest
/// approach's z, from the line through both measurements in arc length and z
/// (`closestApproachZ`), has at most one extremum, a numerical property,
/// not a proven one. So the transverse test is exact at the interval's two
/// ends, and the longitudinal one at the ends of the part within the impact
/// bound and at the one extremum, found by a fixed number of steps. One case
/// apart: a point on the line through the beam line and the chord's midpoint
/// moves outward only on the two diameter circles, tested one by one.
///
/// The beam line is the frame's origin, as for the pair cut (PairCut.h). At
/// an impact bound of 0 and in a field the transverse test is the pair cut's: the
/// circle through the beam point, its curvature at most k_max, no turning
/// point between the two (`transverseCompatible` without its errors). With no
/// field (`maxCurvature` 0) it is the straight line's, which the pair cut does
/// not test.

#ifndef GTRK_TRACKING_DISPLACEDPAIR_H_
#define GTRK_TRACKING_DISPLACEDPAIR_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/Trig.h"

#include "GPUCommonMath.h"

namespace o2::gtrk::displaced
{

/// atan(x) / x, without the singularity at 0.
template <typename T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() T atanOverArgument(T x)
{
  if (o2::gpu::CAMath::Abs(x) < T(0.05)) {
    const T x2 = x * x;
    return T(1) - (x2 * ((T(1) / T(3)) - (x2 * T(0.2))));
  }
  return atanOf(x) / x;
}

/// The circles through a source and a target: the chord's midpoint `M`, its
/// unit direction `u` from the source, its left normal `n`, the half chord.
template <typename T>
struct CircleFamily {
  T sourceX, sourceY, targetX, targetY;
  T midX, midY;
  T ux, uy;
  T nx, ny;
  T half;
  T midSquare; ///< |M|^2
  T midAlong;  ///< M . u
  T midNormal; ///< M . n
};

/// The family through (\a x1, \a y1) and (\a x2, \a y2); a half chord of 0
/// for coincident points, which no test admits.
template <typename T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() CircleFamily<T> circleFamily(T x1, T y1, T x2, T y2)
{
  const T dx = x2 - x1;
  const T dy = y2 - y1;
  const T chord = o2::math_utils::detail::sqrt<T>((dx * dx) + (dy * dy));
  const T inverse = chord > T(0) ? T(1) / chord : T(0);
  const T ux = dx * inverse;
  const T uy = dy * inverse;
  const T midX = T(0.5) * (x1 + x2);
  const T midY = T(0.5) * (y1 + y2);
  return {.sourceX = x1,
          .sourceY = y1,
          .targetX = x2,
          .targetY = y2,
          .midX = midX,
          .midY = midY,
          .ux = ux,
          .uy = uy,
          .nx = -uy,
          .ny = ux,
          .half = T(0.5) * chord,
          .midSquare = (midX * midX) + (midY * midY),
          .midAlong = (ux * midX) + (uy * midY),
          .midNormal = (-uy * midX) + (ux * midY)};
}

/// `sqrt(1 - (k h)^2)`: the cosine of half the turn between the two.
template <typename T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() T halfTurnCosine(const CircleFamily<T>& family, T k)
{
  const T kh = k * family.half;
  const T square = T(1) - (kh * kh);
  return square > T(0) ? o2::math_utils::detail::sqrt<T>(square) : T(0);
}

/// The signed distance from the beam line to the circle of curvature \a k:
/// `|C| - R` for a counterclockwise one, its negative for a clockwise one,
/// the line's signed distance `M . n` at 0, in one smooth expression.
template <typename T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() T circleImpact(const CircleFamily<T>& family, T k)
{
  const T s = halfTurnCosine(family, k);
  const T numerator = (k * (family.midSquare - (family.half * family.half))) + (T(2) * s * family.midNormal);
  const T centre2 = (k * k * family.midSquare) + (T(2) * k * s * family.midNormal) + (s * s);
  return numerator / (T(1) + o2::math_utils::detail::sqrt<T>(centre2 > T(0) ? centre2 : T(0)));
}

/// A closed interval of curvatures; empty if `low > high`.
template <typename T>
struct CurvatureInterval {
  T low;
  T high;
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool empty() const { return !(low <= high); }
};

/// The curvatures for which a pair moves outward at its point \a x, \a y,
/// going from the source to the target: `P . T(k) >= 0`, the tangent `T(k)`,
/// `P . T = k (P x M) + s (P . u)`, one root at most. With `P x M` zero and
/// `P . u` negative only the diameter circles, `s` = 0, move outward there:
/// \a diameterOnly, and the interval is all of [-1/h, 1/h].
template <typename T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() CurvatureInterval<T> outwardAt(const CircleFamily<T>& family, T x, T y, bool& diameterOnly)
{
  const T a = (x * family.midY) - (y * family.midX);
  const T b = (x * family.ux) + (y * family.uy);
  const T edge = family.half > T(0) ? T(1) / family.half : T(0);
  diameterOnly = false;
  if (exactlyEqual(a, T(0))) {
    diameterOnly = b < T(0);
    return {-edge, edge};
  }
  const T root = (a > T(0) ? -b : b) / o2::math_utils::detail::sqrt<T>((a * a) + (b * b * family.half * family.half));
  return a > T(0) ? CurvatureInterval<T>{root, edge} : CurvatureInterval<T>{-edge, root};
}

/// The curvatures at most \a maxCurvature in size on which the pair moves
/// outward at both its points; \a diameterOnly as `outwardAt`'s, for either.
template <typename T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() CurvatureInterval<T> outwardCurvatures(const CircleFamily<T>& family, T maxCurvature, bool& diameterOnly)
{
  diameterOnly = false;
  if (!(family.half > T(0))) {
    return {T(1), T(0)};
  }
  const T edge = T(1) / family.half;
  const T bound = maxCurvature < edge ? maxCurvature : edge;
  bool sourceDiameter = false;
  bool targetDiameter = false;
  const auto source = outwardAt(family, family.sourceX, family.sourceY, sourceDiameter);
  const auto target = outwardAt(family, family.targetX, family.targetY, targetDiameter);
  diameterOnly = sourceDiameter || targetDiameter;
  const T low = o2::gpu::CAMath::Max(-bound, o2::gpu::CAMath::Max(source.low, target.low));
  const T high = o2::gpu::CAMath::Min(bound, o2::gpu::CAMath::Min(source.high, target.high));
  return {low, high};
}

/// The z of the closest approach to the beam line of the helix of
/// curvature \a k through the two measurements, at \a z1 and \a z2: the
/// line through both in arc length and z, there. In the chord's frame, `M`
/// at `(m, v)` and `alpha = asin(k h)` half the turn between them, the
/// closest approach lies `atan2(m sin(alpha), h cos(alpha) + v sin(alpha))`
/// before the arc's middle, in turn: the branch from the geometry, not from
/// the sign of a vanishing product. \a valid false for the circle centred on
/// the beam line, which has no closest approach.
template <typename T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() T closestApproachZ(const CircleFamily<T>& family, T k, T z1, T z2, bool& valid)
{
  T sine = k * family.half;
  if (sine > T(1)) {
    sine = T(1);
  } else if (sine < T(-1)) {
    sine = T(-1);
  }
  const T cosine = halfTurnCosine(family, k);
  const T across = family.midAlong * sine;
  const T along = (family.half * cosine) + (family.midNormal * sine);
  valid = !(exactlyEqual(across, T(0)) && exactlyEqual(along, T(0)));
  if (!valid) {
    return T(0);
  }
  // atan2(across, along) / alpha, smooth through k = 0, where along is h.
  T ratio = T(0);
  if (along > T(0)) {
    // sin(alpha) / alpha = 1 / asinOverArgument(k h).
    ratio = atanOverArgument(across / along) * family.midAlong / (along * o2::gtrk::detail::asinOverArgument(sine));
  } else {
    ratio = atan2Of(across, along) / asinOf(sine);
  }
  return (T(0.5) * (z1 + z2)) - (T(0.5) * (z2 - z1) * ratio);
}

/// What the displaced pair test found.
struct DisplacedPairResult {
  bool transverse = false;   ///< a circle above pT_min, outward, within the impact bound
  bool longitudinal = false; ///< of those, a helix whose closest approach is within the z bound
};

/// Bisection steps of the impact's end points and golden-section steps of
/// the closest approach's z extremum: fixed, so every backend takes the same.
inline constexpr int kBisectionSteps = 40;
inline constexpr int kGoldenSteps = 40;

/// The two diameter circles one by one (`outwardAt`'s \a diameterOnly).
template <typename T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() DisplacedPairResult diameterCircles(const CircleFamily<T>& family, const CurvatureInterval<T>& outward, T z1, T z2,
                                                                             T maxImpact, T zLow, T zHigh)
{
  DisplacedPairResult result;
  const T edge = T(1) / family.half;
  for (int side = 0; side < 2; ++side) {
    const T k = side == 0 ? -edge : edge;
    if (k < outward.low || k > outward.high || o2::gpu::CAMath::Abs(circleImpact(family, k)) > maxImpact) {
      continue;
    }
    result.transverse = true;
    bool valid = false;
    const T z = closestApproachZ(family, k, z1, z2, valid);
    result.longitudinal = result.longitudinal || (valid && z >= zLow && z <= zHigh);
  }
  return result;
}

/// Whether some helix through the two measurements -- (\a x1, \a y1, \a z1)
/// the source, (\a x2, \a y2, \a z2) the target, in the beam line's frame --
/// has a curvature of at most \a maxCurvature, moves outward from the one to
/// the other, passes within \a maxImpact of the beam line, and does so at a z
/// within [\a zLow, \a zHigh].
template <typename T>
[[nodiscard]] GPUhdi() DisplacedPairResult displacedPair(T x1, T y1, T z1, T x2, T y2, T z2, T maxCurvature, T maxImpact, T zLow, T zHigh)
{
  DisplacedPairResult result;
  const auto family = circleFamily(x1, y1, x2, y2);
  bool diameterOnly = false;
  const auto outward = outwardCurvatures(family, maxCurvature, diameterOnly);
  if (outward.empty()) {
    return result;
  }
  if (diameterOnly) {
    return diameterCircles(family, outward, z1, z2, maxImpact, zLow, zHigh);
  }
  const T lowImpact = circleImpact(family, outward.low);
  const T highImpact = circleImpact(family, outward.high);
  const bool rising = lowImpact <= highImpact;
  const T least = rising ? lowImpact : highImpact;
  const T most = rising ? highImpact : lowImpact;
  result.transverse = least <= maxImpact && most >= -maxImpact;
  if (!result.transverse) {
    return result;
  }
  // The part of the interval within the impact bound: the impact is monotone,
  // so the curvature where it crosses each bound, by bisection.
  const auto crossing = [&](T target) {
    T low = outward.low;
    T high = outward.high;
    for (int step = 0; step < kBisectionSteps; ++step) {
      const T middle = T(0.5) * (low + high);
      const bool below = (circleImpact(family, middle) < target) == rising;
      (below ? low : high) = middle;
    }
    return T(0.5) * (low + high);
  };
  T first = outward.low;
  T last = outward.high;
  if (least < -maxImpact) {
    (rising ? first : last) = crossing(-maxImpact);
  }
  if (most > maxImpact) {
    (rising ? last : first) = crossing(maxImpact);
  }
  if (first > last) {
    const T swap = first;
    first = last;
    last = swap;
  }
  bool firstValid = false;
  bool lastValid = false;
  const T firstZ = closestApproachZ(family, first, z1, z2, firstValid);
  const T lastZ = closestApproachZ(family, last, z1, z2, lastValid);
  if (!firstValid || !lastValid) {
    // A circle centred on the beam line at an end has no closest approach,
    // so no z to test: refused, as there is none to admit.
    return result;
  }
  const auto inside = [&](T z) { return z >= zLow && z <= zHigh; };
  if (inside(firstZ) || inside(lastZ) || ((firstZ < zLow) != (lastZ < zLow))) {
    result.longitudinal = true;
    return result;
  }
  // Both ends on one side of the bound: the one extremum between them, a
  // minimum if they are above it, a maximum if below, by golden section.
  const T sign = firstZ > zHigh ? T(1) : T(-1);
  const T ratio = T(0.6180339887498949);
  T low = first;
  T high = last;
  T left = high - (ratio * (high - low));
  T right = low + (ratio * (high - low));
  bool valid = false;
  T leftZ = sign * closestApproachZ(family, left, z1, z2, valid);
  T rightZ = sign * closestApproachZ(family, right, z1, z2, valid);
  for (int step = 0; step < kGoldenSteps; ++step) {
    if (leftZ < rightZ) {
      high = right;
      right = left;
      rightZ = leftZ;
      left = high - (ratio * (high - low));
      leftZ = sign * closestApproachZ(family, left, z1, z2, valid);
    } else {
      low = left;
      left = right;
      leftZ = rightZ;
      right = low + (ratio * (high - low));
      rightZ = sign * closestApproachZ(family, right, z1, z2, valid);
    }
  }
  // Both ends above the bound: the least z must reach it; both below, the most.
  const T extremum = sign * (leftZ < rightZ ? leftZ : rightZ);
  result.longitudinal = sign > T(0) ? extremum <= zHigh : extremum >= zLow;
  return result;
}

/// The errors a displaced pair test takes: each measurement's, the
/// scattering, the beam's width [cm^2].
template <typename T>
struct PairErrors {
  T sourceVarR = T(0);
  T sourceVarRPhi = T(0);
  T sourceVarZ = T(0);
  T targetVarR = T(0);
  T targetVarRPhi = T(0);
  T targetVarZ = T(0);
  /// `sum theta_k^2 s_k^2` over the layers up to and including the source's,
  /// at pT_min (`sourceScattering`): the closest approach's displacement by a
  /// kink inside the source, its lever arm taken from the vertex.
  T sourceScatter = T(0);
  /// The target's displacement by the passives between the two
  /// (`betweenScattering`).
  T betweenScatter = T(0);
  T beamVarTransverse = T(0);
};

/// One sigma each of the impact and of the closest approach's z.
template <typename T>
struct PairSigmas {
  T impact;
  T z;
};

/// One sigma each of the impact and of the closest approach's z, bounds that
/// hold for every circle the test may take, so that widening the impact and
/// z bounds by them keeps the test's monotone and one-extremum structure
/// (a sigma per curvature would not). At a fixed
/// curvature the closest approach is the two points' extrapolation, with
/// `s_i` their arcs from it: `delta d0 = (s_2 delta_1 - s_1 delta_2) / (s_2 -
/// s_1)`, `delta` across the track, and z the same in z. The arc between them
/// is at least the chord `c`; an arc from the closest approach at most its
/// chord, `r_i + D`, times the largest arc-to-chord ratio at k_max,
/// `asin(x) / x` at `x = k_max (r_i + D) / 2` (pi / 2 at half a turn). Across
/// the track a point's error is at most `varR + varRPhi`. The transverse
/// errors move the arcs, which the dip `|z_2 - z_1| / c` turns into z; a
/// kink inside the source moves the closest approach by `sourceScatter`'s
/// lever, along the track's 3D path in z.
template <typename T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() PairSigmas<T> displacedPairSigmas(T x1, T y1, T z1, T x2, T y2, T z2, const PairErrors<T>& errors, T maxCurvature,
                                                                           T maxImpact)
{
  const T dx = x2 - x1;
  const T dy = y2 - y1;
  const T chord2 = (dx * dx) + (dy * dy);
  if (!(chord2 > T(0))) {
    return {T(0), T(0)};
  }
  const T r1 = o2::math_utils::detail::sqrt<T>((x1 * x1) + (y1 * y1)) + maxImpact;
  const T r2 = o2::math_utils::detail::sqrt<T>((x2 * x2) + (y2 * y2)) + maxImpact;
  const auto ratio = [&](T reach) {
    const T x = T(0.5) * maxCurvature * reach;
    return o2::gtrk::detail::asinOverArgument(x < T(1) ? x : T(1));
  };
  const T lever1 = ratio(r2) * r2; // the source's error's lever: the target's arc
  const T lever2 = ratio(r1) * r1;
  const T across1 = errors.sourceVarR + errors.sourceVarRPhi;
  const T across2 = errors.targetVarR + errors.targetVarRPhi + errors.betweenScatter;
  const T position = ((lever1 * lever1 * across1) + (lever2 * lever2 * across2)) / chord2;
  const T dz = z2 - z1;
  const T dip2 = dz * dz / chord2;
  const T longitudinal = ((lever1 * lever1 * errors.sourceVarZ) + (lever2 * lever2 * errors.targetVarZ)) / chord2;
  return {.impact = o2::math_utils::detail::sqrt<T>(position + errors.sourceScatter + errors.beamVarTransverse),
          .z = o2::math_utils::detail::sqrt<T>(longitudinal + (dip2 * position) + ((T(1) + dip2) * errors.sourceScatter))};
}

/// `displacedPair` with its impact and z bounds widened by \a nSigma of
/// `displacedPairSigmas`.
template <typename T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() DisplacedPairResult displacedPairWithin(T x1, T y1, T z1, T x2, T y2, T z2, const PairErrors<T>& errors, T nSigma,
                                                                                 T maxCurvature, T maxImpact, T zLow, T zHigh)
{
  const auto sigmas = displacedPairSigmas(x1, y1, z1, x2, y2, z2, errors, maxCurvature, maxImpact);
  const T impact = maxImpact + (nSigma * sigmas.impact);
  const T z = nSigma * sigmas.z;
  return displacedPair(x1, y1, z1, x2, y2, z2, maxCurvature, impact, zLow - z, zHigh + z);
}

} // namespace o2::gtrk::displaced

#endif // GTRK_TRACKING_DISPLACEDPAIR_H_
