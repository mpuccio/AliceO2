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

/// \file SearchWindow.h
/// \brief How wide a candidate window is, derived rather than tuned.
///
/// Every number here follows from the minimum pT of the current iteration, the
/// layers' envelopes measured from the loaded geometry (Base/LayerEnvelopes.h),
/// the longitudinal vertex constraint and
/// the Highland term the track model already carries. Nothing is measured off
/// truth chains, which is what lets a grid granularity be derived at all.
///
/// The split that makes this usable on a device: a window depends only on the
/// layer pair, iteration and collision vertex, never on the seed. Only its
/// centre depends on the seed, and that is one multiply. The transcendental
/// work -- asin, the Highland logarithm -- happens once per pair, iteration
/// and vertex on the host, into a `LayerPairWindow` small enough to upload.
/// The per-seed query is arithmetic.

#ifndef GTRK_TRACKING_SEARCHWINDOW_H_
#define GTRK_TRACKING_SEARCHWINDOW_H_

#include "GenericTracking/Base/LayerEnvelopes.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/DisplacedPair.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/MathUtils.h"

#include "CommonConstants/PhysicsConstants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace o2::gtrk
{

/// The iteration-wide physics inputs. A collision's longitudinal constraint
/// is separate because one timeframe can contain many vertices at different z.
struct SearchWindowParams {
  // This is host-side configuration. Window derivation evaluates in double
  // and narrows once when it creates the device-facing LayerPairWindow.
  /// Minimum transverse momentum of the current CA iteration [GeV/c].  The
  /// curvature term and the Highland term both scale as 1/p, so an early
  /// high-pT pass has a window an order of magnitude narrower on both counts.
  double minPtGeV = 0.1;

  /// Solenoid field along z [kG], the unit `o2::constants::math::B2C` and the
  /// rest of O2 are written in.  2 T is 20 kG.
  double bzKiloGauss = 20.;

  /// The mass hypothesis [GeV] (`Propagator::massGeV`, positive): the softest
  /// momentum the windows admit scatters as 1 / (beta p).
  double massGeV = o2::constants::physics::MassPionCharged;

  /// Radiation lengths crossed between two layers: the uniform material that
  /// sizes the grid's cells (`uniformSourceKink`); the pair queries take the
  /// scattering table's per-layer values (`layerPairWindow` in PairCut.h).
  /// 1 % X0 is the order of a pixel layer with its support.
  double materialX2X0 = 0.01;

  /// Cluster position resolution entering the window [cm]. One number for
  /// every layer.
  double positionSigma = 0.005;

  /// How many sigma of the Gaussian terms -- vertex, scattering, resolution --
  /// the window admits.  The curvature term is a bound, not a sigma, and is
  /// not scaled by this.
  double nSigma = 3.;

  /// A displaced iteration's bounds: the
  /// largest distance its tracks' circles may pass from the beam line, and
  /// how far the z of that closest approach may lie from the beam spot's
  /// centre [cm]. 0, the default, is a prompt iteration, whose pair test and
  /// windows are the beam point's; a displaced one tests a pair with
  /// `displacedPairWithin` and its windows hold what that admits.
  double impactMax = 0.;
  double impactZMax = 0.;

  [[nodiscard]] bool operator==(const SearchWindowParams&) const = default;
};

/// The origin shared by tracks from one collision.
///
/// Before vertexing this is the broad beam-spot constraint. Afterwards there
/// is one of these per reconstructed collision, with that collision's z and
/// uncertainty. It is deliberately not part of SearchWindowParams.
///
/// Transverse positions are in the frame of the corrected beam line, whose
/// position and width come from the MeanVertex calibration of the run. The producer subtracts the transverse position from
/// every measurement (Workflow/BeamLine.h); the constraint carries the widths
/// and the longitudinal position.
struct VertexConstraint {
  double z = 0.;                  ///< [cm]
  double sigmaZ = 6.;             ///< [cm]
  double sigmaTransverse = 0.002; ///< beam width across the beam line [cm]
  /// The beam line's transverse position in the fit's frame, the detector's
  /// [cm]: the finding works about it, in the corrected frame, and reads
  /// neither (BeamLine.h); the beam-line update does (BeamConstraint.h).
  double x = 0.;
  double y = 0.;
};

/// Where a layer is [cm], read off the envelopes measured from the geometry
/// loaded in this process (Base/LayerEnvelopes.h). Before a Geometry<Det>
/// exists they throw std::logic_error.
///
/// `layerDistance`: the layer's position (`LayerEnvelope::position`, the mean
/// of its sensors' radii or z), unsigned: the radius of a barrel layer, the |z|
/// of a disc. `discZ`: a disc's signed z,
/// which the lever arm needs; 0 on a barrel layer. `secondMin`/`secondMax`:
/// the range of the second binning coordinate, z on a barrel layer, r on a
/// disc. `radiusMin`/`radiusMax`, `zMin`/`zMax`: the ranges of r and z the
/// layer's measurements span.
template <DetectorLayout Det>
[[nodiscard]] inline double layerDistance(LayerID layer)
{
  return std::abs(PlacedEnvelopes<Det>::get(layer).position);
}

template <DetectorLayout Det>
[[nodiscard]] inline double discZ(LayerID layer)
{
  return Detector<Det>::isBarrelLayer(layer) ? 0. : PlacedEnvelopes<Det>::get(layer).position;
}

template <DetectorLayout Det>
[[nodiscard]] inline double secondMin(LayerID layer)
{
  return Detector<Det>::isBarrelLayer(layer) ? PlacedEnvelopes<Det>::get(layer).zMin : PlacedEnvelopes<Det>::get(layer).rMin;
}

template <DetectorLayout Det>
[[nodiscard]] inline double secondMax(LayerID layer)
{
  return Detector<Det>::isBarrelLayer(layer) ? PlacedEnvelopes<Det>::get(layer).zMax : PlacedEnvelopes<Det>::get(layer).rMax;
}

template <DetectorLayout Det>
[[nodiscard]] inline double radiusMin(LayerID layer)
{
  return PlacedEnvelopes<Det>::get(layer).rMin;
}

template <DetectorLayout Det>
[[nodiscard]] inline double radiusMax(LayerID layer)
{
  return PlacedEnvelopes<Det>::get(layer).rMax;
}

template <DetectorLayout Det>
[[nodiscard]] inline double zMin(LayerID layer)
{
  return PlacedEnvelopes<Det>::get(layer).zMin;
}

template <DetectorLayout Det>
[[nodiscard]] inline double zMax(LayerID layer)
{
  return PlacedEnvelopes<Det>::get(layer).zMax;
}

/// A range of values, and the ratio of two positions spanning ranges.
struct Interval {
  double low = 0.;
  double high = 0.;
};

/// The range of `(to - origin) / (from - origin)` for `to` and `from` within
/// their intervals. Neither interval may contain `origin`.
[[nodiscard]] inline Interval ratioRange(Interval to, Interval from, double origin)
{
  const double corners[4] = {(to.low - origin) / (from.low - origin), (to.low - origin) / (from.high - origin), (to.high - origin) / (from.low - origin), (to.high - origin) / (from.high - origin)};
  return {.low = *std::min_element(corners, corners + 4), .high = *std::max_element(corners, corners + 4)};
}

/// `ratioRange` over every origin in \a origin: the union of its ends'.
[[nodiscard]] inline Interval ratioRange(Interval to, Interval from, Interval origin)
{
  const auto low = ratioRange(to, from, origin.low);
  const auto high = ratioRange(to, from, origin.high);
  return {.low = std::min(low.low, high.low), .high = std::max(low.high, high.high)};
}

/// The largest magnitude in \a range.
[[nodiscard]] inline double largestAbs(Interval range) { return std::max(std::abs(range.low), std::abs(range.high)); }

/// The distances between \a a and \a b: the nearest, zero where they
/// overlap, and the farthest, between their ends.
[[nodiscard]] inline double nearestApart(Interval a, Interval b) { return std::max({a.low - b.high, b.low - a.high, 0.}); }
[[nodiscard]] inline double farthestApart(Interval a, Interval b) { return std::max({std::abs(b.low - a.low), std::abs(b.low - a.high), std::abs(b.high - a.low), std::abs(b.high - a.high)}); }

/// The z the vertex constraint reaches: `nSigma` of its width either side.
[[nodiscard]] inline Interval vertexReach(const SearchWindowParams& params, const VertexConstraint& vertex)
{
  const double reach = params.nSigma * vertex.sigmaZ;
  return {.low = vertex.z - reach, .high = vertex.z + reach};
}

/// A layer's extents as intervals: its transverse radii, its z, its second
/// coordinate (`radiusMin`, `zMin`, `secondMin` and their maxima).
template <DetectorLayout Det>
[[nodiscard]] inline Interval radiusRange(LayerID layer)
{
  return {.low = radiusMin<Det>(layer), .high = radiusMax<Det>(layer)};
}
template <DetectorLayout Det>
[[nodiscard]] inline Interval zRange(LayerID layer)
{
  return {.low = zMin<Det>(layer), .high = zMax<Det>(layer)};
}
template <DetectorLayout Det>
[[nodiscard]] inline Interval secondRange(LayerID layer)
{
  return {.low = secondMin<Det>(layer), .high = secondMax<Det>(layer)};
}

/// Whether a position lies strictly inside what layers \a a and \a b span
/// together: a radius (\a barrel) between their least and largest radius, a
/// z between their least and largest z.
template <DetectorLayout Det>
[[nodiscard]] inline bool positionBetweenLayers(bool barrel, double position, LayerID a, LayerID b)
{
  const Interval spanA = barrel ? radiusRange<Det>(a) : zRange<Det>(a);
  const Interval spanB = barrel ? radiusRange<Det>(b) : zRange<Det>(b);
  return position > std::min(spanA.low, spanB.low) && position < std::max(spanA.high, spanB.high);
}

/// The largest transverse curvature the windows admit [1/cm]: `1/R` at
/// pT_min; zero without a field, unbounded without a pT cut.
[[nodiscard]] inline double maxCurvature(const SearchWindowParams& params)
{
  if (exactlyEqual(params.bzKiloGauss, 0.)) {
    return 0.;
  }
  const double radius = curvatureRadius(params.minPtGeV, params.bzKiloGauss);
  return radius > 0. ? 1. / radius : std::numeric_limits<double>::infinity();
}

/// 1 / (beta p)^2 [GeV^-2] at the softest momentum the windows admit, pT_min:
/// what their scattering scales as, and an upper bound over every momentum
/// above it.
[[nodiscard]] inline double invBetaMomentumMin2(const SearchWindowParams& params)
{
  const double pt2 = params.minPtGeV * params.minPtGeV;
  return (1. + (params.massGeV * params.massGeV / pt2)) / pt2;
}

/// `1 / (beta pT_min)^2` over `1 / pT_min^2`: what the mass hypothesis
/// multiplies the scattering at pT_min by.
[[nodiscard]] inline double betaMomentumFactor(const SearchWindowParams& params) { return invBetaMomentumMin2(params) * params.minPtGeV * params.minPtGeV; }

/// The uniform stand-in for a source's transverse kink [rad]: `materialX2X0`
/// at normal incidence in the seed's layer alone. It sizes the grid's cells;
/// the pair queries take the pair cut's own table (`sourceKinkBound` in
/// `PairCut.h`).
[[nodiscard]] inline double uniformSourceKink(const SearchWindowParams& params) { return static_cast<double>(highlandTheta0(static_cast<Scalar>(params.materialX2X0), betaMomentum(static_cast<Scalar>(params.minPtGeV), static_cast<Scalar>(params.massGeV)))); }

/// How far the scattering before the seed moves the target's second
/// coordinate, in one standard deviation [cm], for vertices in \a vertexZ and
/// \a sourceKink bounding `sqrt(S)/r_i` over the seed's layer.
///
/// The longitudinal pair test carries the scattering as `(ds^2 + dz^2) S` in
/// the variance of `D = s_i (z_j - z_v) - s_j (z_i - z_v)`, with `S` the
/// source's scattering sum. Moving the target by `delta` changes `D` by
/// `s_i delta` on a barrel, in z, and by `|z_i - z_v| delta` on a disc, in r.
/// One standard deviation is therefore `hypot(ds, dz) sqrt(S)` over that
/// factor. The step `ds` along the transverse arc is at most the target's arc
/// less the seed's radius, with the pair cut's curvature bound; `dz` is at
/// most the pair's z extent and `ds` times the steepest slope from the vertex.
/// On a barrel `sqrt(S)/s_i <= sqrt(S)/r_i`. On a disc the factor is `r_i/|z_i -
/// z_v|` from a disc seed; from a barrel seed, `s_i/|z_i - z_v| = s_j/|z_j -
/// z_v|` for a pair the pitch test can pass, bounded at the target.
///
/// \a betweenKink bounds `sqrt(T)`, the target's displacement by the passive
/// scatterers between the two layers (`betweenKinkBound` in `PairCut.h`). It enters the variance of `D` as
/// `(s_i^2 + (z_i - z_v)^2) T`, which moves the target's second coordinate by
/// `sqrt(T) hypot(s_i, z_i - z_v)` over the same factor `s_i` or `|z_i - z_v|`.
template <DetectorLayout Det>
[[nodiscard]] inline double secondScattering(LayerID from, LayerID to, const SearchWindowParams& params, Interval vertexZ, double sourceKink, double betweenKink = 0.)
{
  const double rFromLow = radiusMin<Det>(from);
  const double rToHigh = radiusMax<Det>(to);
  if (!(rFromLow > 0.) || !(rToHigh > 0.)) {
    return std::numeric_limits<double>::infinity();
  }
  const double curvature = std::min(maxCurvature(params), 2. / rToHigh);
  const double targetArc = transverseArcAtRadius(curvature, rToHigh);
  const double arcStep = std::max(targetArc - rFromLow, 0.);
  const Interval zFrom = zRange<Det>(from);
  const Interval zTo = zRange<Det>(to);
  const double zStep = std::max(std::abs(zTo.high - zFrom.low), std::abs(zTo.low - zFrom.high));
  const double sourceFar = farthestApart(vertexZ, zFrom);
  const double step = std::hypot(arcStep, std::min(zStep, arcStep * sourceFar / rFromLow));
  if (Detector<Det>::isBarrelLayer(to)) {
    return std::hypot(sourceKink * step, betweenKink * std::hypot(1., sourceFar / rFromLow));
  }
  const double near = nearestApart(vertexZ, Detector<Det>::isBarrelLayer(from) ? zTo : zFrom);
  if (!(near > 0.)) {
    return std::numeric_limits<double>::infinity();
  }
  const double lever = Detector<Det>::isBarrelLayer(from) ? targetArc : radiusMax<Det>(from);
  return std::hypot(sourceKink * step * lever / near, betweenKink * std::hypot(1., lever / near));
}

/// Everything a query needs about one ordered layer pair, computed once per
/// pair, CA iteration and vertex constraint.
///
/// For a same-kind pair, `secondScale` and `secondOffset` are the
/// straight-line extrapolation from the collision vertex, and
/// `secondHalfWidth` is what its z uncertainty, scattering and the cluster
/// resolution do to that estimate.
///
/// A barrel <-> disc pair has no such line: the target's coordinate goes as
/// the reciprocal of the seed's. Its window is evaluated per seed from the
/// fields after `kind`, and `secondHalfWidth` is the padding added to it.
struct LayerPairWindow {
  Scalar phiHalfWidth = 0.F;    ///< [rad], symmetric about the seed's phi
  Scalar secondScale = 0.F;     ///< d(second on target) / d(second on seed)
  Scalar secondOffset = 0.F;    ///< [cm], the constant term of that extrapolation
  Scalar secondHalfWidth = 0.F; ///< [cm], about secondScale * seed + secondOffset

  PairKind kind = PairKind::kBarrel;
  Scalar targetLow = 0.F; ///< [cm] the target disc's z, or the target barrel layer's radii
  Scalar targetHigh = 0.F;
  Scalar vertexLow = 0.F; ///< [cm] the vertex's z, n sigma either side
  Scalar vertexHigh = 0.F;
  Scalar seedSpread = 0.F;      ///< [cm] n sigma of the seed's own second coordinate
  Scalar targetSecondMin = 0.F; ///< [cm] the target's extent in its second coordinate
  Scalar targetSecondMax = 0.F;
};

/// Host-derived terms shared by every per-vertex second-coordinate interval
/// of one layer pair. The interval centre still follows the individual vertex
/// on the device.
struct VertexPairWindow {
  Scalar secondScale = 0.F;    ///< barrel radius ratio
  Scalar sourcePosition = 0.F; ///< source disc z for a disc pair
  Scalar targetPosition = 0.F; ///< target disc z for a disc pair
  Scalar vertexLever = 0.F;    ///< second-coordinate shift per vertex-z shift
  Scalar secondVariance = 0.F; ///< scattering and measurement variance [cm^2]
  Scalar secondPadding = 0.F;  ///< centre, helix and rounding reach [cm]
};

/// Narrow a nonnegative host-derived bound without making it smaller.
[[nodiscard]] inline Scalar conservativeScalarBound(double value)
{
  auto narrowed = static_cast<Scalar>(value);
  if (static_cast<double>(narrowed) < value) {
    narrowed = std::nextafter(narrowed, std::numeric_limits<Scalar>::infinity());
  }
  return narrowed;
}

/// The half-width in phi: the widest azimuth the transverse pair test accepts.
///
/// The test admits a pair when
///
/// ```
/// sin(dphi) <= d / (2 R_min) + n sqrt(res^2 + (d/r_j)^2 kink^2)
/// d^2 = r_i^2 + r_j^2 - 2 r_i r_j cos(dphi),    res^2 = sigma^2 (1/r_i^2 + 1/r_j^2)
/// cos(dphi) >= min(r_i, r_j) / max(r_i, r_j)
/// ```
///
/// with \a sourceKink one standard deviation of the track's transverse
/// direction at the seed [rad] and the cluster resolution `positionSigma`.
/// The last line is the test's no-turning-point condition. The chord grows
/// with dphi, so no closed form bounds the accepted set; it is scanned from
/// the quarter turn down over the pair's radius ranges. Between two discs the target radius follows the
/// seed's through the z ratio `rho`, pulled in by the helix to
/// `r_i sin(rho a)/sin(a)` for the seed's half-turn `a`.
///
/// \a betweenKink bounds `sqrt(T)` as for `secondScattering`: the transverse
/// test's variance gains `r_i^2 T`, which in the accepted `sin(dphi)` is
/// `T / r_j^2`.
template <DetectorLayout Det>
[[nodiscard]] inline double phiWindow(LayerID from, LayerID to, const SearchWindowParams& params, const VertexConstraint& vertex, double sourceKink, double betweenKink = 0.)
{
  constexpr double kQuarter = 0.5 * std::numbers::pi_v<double>;
  constexpr double kStep = 1.e-3;
  constexpr int kSteps = static_cast<int>(kQuarter / kStep) + 1;
  constexpr int kSamples = 5;
  const double radius = curvatureRadius(params.minPtGeV, params.bzKiloGauss);
  const double fromLow = radiusMin<Det>(from);
  const double fromHigh = radiusMax<Det>(from);
  if (!(radius > 0.) || !(fromLow > 0.)) {
    return kQuarter;
  }
  const double twoR = 2. * radius;

  const bool discs = !Detector<Det>::isBarrelLayer(from) && !Detector<Det>::isBarrelLayer(to);
  const auto rho = ratioRange(zRange<Det>(to), zRange<Det>(from), vertexReach(params, vertex));
  const double rhoLow = rho.low;
  const double rhoHigh = std::max(rho.high, 0.);

  const double variance = params.positionSigma * params.positionSigma;
  const auto accepted = [&](double rFrom, double rTo, double phi) {
    if (std::cos(phi) * std::max(rFrom, rTo) < std::min(rFrom, rTo)) {
      return false;
    }
    const double chord = std::sqrt(std::max((rFrom * rFrom) + (rTo * rTo) - (2. * rFrom * rTo * std::cos(phi)), 0.));
    const double lever = chord / rTo;
    const double gaussian = std::sqrt((variance * ((1. / (rFrom * rFrom)) + (1. / (rTo * rTo)))) + (lever * lever * sourceKink * sourceKink) + (betweenKink * betweenKink / (rTo * rTo)));
    return std::sin(phi) <= (chord / twoR) + (params.nSigma * gaussian);
  };

  double widest = 0.;
  for (int i = 0; i < kSamples; ++i) {
    const double rFrom = fromLow + ((fromHigh - fromLow) * i / (kSamples - 1));
    double toLow = radiusMin<Det>(to);
    double toHigh = radiusMax<Det>(to);
    if (discs) {
      const double turn = std::min(std::asin(std::min(rFrom / twoR, 1.)), 0.5 * std::numbers::pi_v<double> / std::max(rhoLow, 1.));
      toLow = std::max(toLow, turn > 0. ? rFrom * std::sin(rhoLow * turn) / std::sin(turn) : rhoLow * rFrom);
      toHigh = std::min(toHigh, rhoHigh * rFrom);
      if (toHigh < toLow) {
        continue;
      }
    }
    for (int j = 0; j < kSamples; ++j) {
      const double rTo = toLow + ((toHigh - toLow) * j / (kSamples - 1));
      if (!(rTo > 0.)) {
        continue;
      }
      for (int step = kSteps; step >= 0 && kStep * (step + 1) > widest; --step) {
        if (accepted(rFrom, rTo, kStep * step)) {
          widest = kStep * (step + 1);
          break;
        }
      }
    }
  }
  return std::min(widest, kQuarter);
}

/// The lever arm of the beam-line constraint on the target layer's second
/// coordinate, and the window it implies around a straight extrapolation.
///
/// ```
/// barrel -> barrel   z_j = z_v + (z_i - z_v) * r_j/r_i      dz_j/dz_v = r_j/r_i - 1
/// disc   -> disc     r_j = r_i * (z_j - z_v)/(z_i - z_v)    dr_j/dz_v = r_i (z_j - z_i)/(z_i - z_v)^2
/// ```
///
/// The two behave oppositely: `r_j/r_i - 1` is about 1 on an inner barrel
/// pair, so the spot's spread passes on nearly whole and the constraint is
/// nearly worthless there; between two discs the lever arm in z is long and
/// the same spot costs under a millimetre.
///
/// For a same-kind pair (a crossing pair takes `crossingPairWindow`).
/// \a seedSecond is the seed's own second coordinate.  Barrel-to-barrel does
/// not use it -- the derivative is independent of z_i -- but disc-to-disc does,
/// and passing the target-side extreme radius there is the conservative choice
/// a precomputed per-pair window has to make. \a sourceKink is as for
/// `secondScattering`.
template <DetectorLayout Det>
[[nodiscard]] inline double secondWindow(LayerID from, LayerID to, double seedSecond, const SearchWindowParams& params, const VertexConstraint& vertexConstraint, double sourceKink, double betweenKink = 0.)
{
  const double fullExtent = secondMax<Det>(to) - secondMin<Det>(to);
  const Interval vertexZ = vertexReach(params, vertexConstraint);

  // Taken over the layers' measured ranges, where each term is largest.
  if (Detector<Det>::isBarrelLayer(from) && Detector<Det>::isBarrelLayer(to)) {
    const auto ratio = ratioRange(radiusRange<Det>(to), radiusRange<Det>(from), 0.);
    const double vertexUncertainty = std::max(std::abs(ratio.high - 1.), std::abs(ratio.low - 1.)) * vertexConstraint.sigmaZ;
    const double scattering = secondScattering<Det>(from, to, params, vertexZ, sourceKink, betweenKink);
    const double resolution = params.positionSigma * ratio.high;
    return params.nSigma * std::hypot(vertexUncertainty, scattering, resolution);
  }

  // A disc pair.
  const Interval zFrom = zRange<Det>(from);
  const Interval zTo = zRange<Det>(to);
  const double fromVertex = std::min(std::abs(zFrom.low - vertexConstraint.z), std::abs(zFrom.high - vertexConstraint.z));
  if (!(fromVertex > 0.)) {
    return fullExtent;
  }
  const double gap = std::max(std::abs(zTo.high - zFrom.low), std::abs(zTo.low - zFrom.high));
  const auto ratio = ratioRange(zTo, zFrom, vertexConstraint.z);
  const double vertexUncertainty = std::abs(seedSecond) * gap / (fromVertex * fromVertex) * vertexConstraint.sigmaZ;
  const double scattering = secondScattering<Det>(from, to, params, vertexZ, sourceKink, betweenKink);
  const double resolution = params.positionSigma * std::max(std::abs(ratio.low), std::abs(ratio.high));
  return params.nSigma * std::hypot(vertexUncertainty, scattering, resolution);
}

/// Extra second-coordinate reach when z is linear in transverse arc rather
/// than radius. The query remains centred on its cheap straight extrapolation.
///
/// With the arc `s = r f(k r/2)`, `f(x) = asin(x)/x`, a helix from the vertex
/// through the seed at `r_i` lands off the straight line through the two by
///
/// ```
/// target on a barrel, in z    |z_i - z_v| (r_j/r_i) (f_j/f_i - 1)
/// target on a disc, in r      r_j (f_j/f_i - 1)
/// ```
///
/// Every link leads outwards in radius, and both grow with `r_j` and with the
/// curvature and shrink with `r_i`: the smallest seed radius and the largest
/// target radius bound them. The pair cut holds the curvature to pT_min's, and
/// a circle through the beam point does not reach past its diameter, so the
/// curvature is at most `min(1/R_min, 2/r_j)`.
///
/// Between two discs the seed radius is not free either: the arc to it is the
/// target's over `rho = (z_j - z_v)/(z_i - z_v)`. The offset `rho r_i - r_j`
/// grows with `rho`, so it is bounded at the largest target radius and `rho`.
template <DetectorLayout Det>
[[nodiscard]] inline double helixSecondSpread(LayerID from, LayerID to, const SearchWindowParams& params, const VertexConstraint& vertex)
{
  if (exactlyEqual(params.bzKiloGauss, 0.)) {
    return 0.;
  }
  const double curvatureMax = maxCurvature(params);
  const Interval vertexZ = vertexReach(params, vertex);

  if (!Detector<Det>::isBarrelLayer(from) && !Detector<Det>::isBarrelLayer(to)) {
    const double rho = std::max(0., ratioRange(zRange<Det>(to), zRange<Det>(from), vertexZ).high);
    const double targetRadius = radiusMax<Det>(to);
    if (!(rho > 1.) || !(targetRadius > 0.)) {
      return 0.;
    }
    const double curvature = std::min(curvatureMax, 2. / targetRadius);
    const double sourceArc = transverseArcAtRadius(curvature, targetRadius) / rho;
    const double sourceRadius = 2. * std::sin(0.5 * curvature * sourceArc) / curvature;
    return (rho * sourceRadius) - targetRadius;
  }

  const double sourceRadius = radiusMin<Det>(from);
  const double targetRadius = radiusMax<Det>(to);
  if (!(sourceRadius > 0.) || !(targetRadius > sourceRadius)) {
    return 0.;
  }
  const double curvature = std::min(curvatureMax, 2. / targetRadius);
  const double arcExcess = (transverseArcAtRadius(curvature, targetRadius) / targetRadius) / (transverseArcAtRadius(curvature, sourceRadius) / sourceRadius) - 1.;
  if (!Detector<Det>::isBarrelLayer(to)) {
    return targetRadius * arcExcess;
  }
  const double lever = farthestApart(vertexZ, zRange<Det>(from));
  return lever * (targetRadius / sourceRadius) * arcExcess;
}

/// The straight-line extrapolation factor for the second coordinate about the
/// supplied vertex, for a same-kind pair (a crossing pair's centre is the
/// target layer's midpoint instead -- see `secondOffset`).
template <DetectorLayout Det>
[[nodiscard]] inline double secondScale(LayerID from, LayerID to, double vertexZ)
{
  if (Detector<Det>::isBarrelLayer(from) && Detector<Det>::isBarrelLayer(to)) {
    const double rFrom = layerDistance<Det>(from);
    return rFrom > 0. ? layerDistance<Det>(to) / rFrom : 0.;
  }
  const double fromVertex = discZ<Det>(from) - vertexZ;
  return std::abs(fromVertex) > 0. ? (discZ<Det>(to) - vertexZ) / fromVertex : 0.;
}

/// The constant term of the same extrapolation.
///
/// A barrel pair uses `(1 - scale) * vertexZ`; a disc pair needs no offset
/// because its vertex dependence is already in the scale. It carries the
/// target layer's midpoint for the undetermined barrel-to-disc case, where
/// the window is the layer's whole extent and a disc's r extent is not centred
/// on zero the way a barrel's z extent is.
template <DetectorLayout Det>
[[nodiscard]] inline double secondOffset(LayerID from, LayerID to, double vertexZ)
{
  if (Detector<Det>::isBarrelLayer(from) && Detector<Det>::isBarrelLayer(to)) {
    return (1. - secondScale<Det>(from, to, vertexZ)) * vertexZ;
  }
  if (!Detector<Det>::isBarrelLayer(from) && !Detector<Det>::isBarrelLayer(to)) {
    return 0.;
  }
  return 0.5 * (secondMin<Det>(to) + secondMax<Det>(to));
}

/// How far the straight line from the vertex through an actual seed lands
/// from the window's centre, which uses the layers' nominal positions [cm].
///
/// Barrel: the line reaches the target at `z_v + (z_i - z_v) r_j/r_i` with the
/// measured radii, against `scale` from the nominal ones. Discs: at
/// `r_i (z_j - z_v)/(z_i - z_v)` with the measured z. The spread is added to
/// the half-width, so the cells picked still contain every row the exact pair
/// cuts can accept. For a same-kind pair.
template <DetectorLayout Det>
[[nodiscard]] inline double secondCentreSpread(LayerID from, LayerID to, double vertexZ)
{
  const double scale = secondScale<Det>(from, to, vertexZ);
  if (Detector<Det>::isBarrelLayer(from) && Detector<Det>::isBarrelLayer(to)) {
    const auto ratio = ratioRange(radiusRange<Det>(to), radiusRange<Det>(from), 0.);
    const double lever = std::max(std::abs(secondMin<Det>(from) - vertexZ), std::abs(secondMax<Det>(from) - vertexZ));
    return lever * std::max(ratio.high - scale, scale - ratio.low);
  }
  const auto ratio = ratioRange(zRange<Det>(to), zRange<Det>(from), vertexZ);
  return secondMax<Det>(from) * std::max(ratio.high - scale, scale - ratio.low);
}

struct VertexSecondWindow {
  double vertexLever = 0.;
  double variance = 0.;
  double padding = 0.;
};

/// The same-kind second-coordinate uncertainty terms, bounded over the z
/// range the vertices an iteration uses can take.
template <DetectorLayout Det>
[[nodiscard]] inline VertexSecondWindow vertexSecondWindow(LayerID from, LayerID to, double seedSecond, const SearchWindowParams& params, Interval vertexZ, double sourceKink, double betweenKink = 0.)
{
  const double scattering = secondScattering<Det>(from, to, params, vertexZ, sourceKink, betweenKink);
  if (Detector<Det>::isBarrelLayer(from) && Detector<Det>::isBarrelLayer(to)) {
    const auto ratio = ratioRange(radiusRange<Det>(to), radiusRange<Det>(from), 0.);
    const double vertexLever = std::max(std::abs(ratio.high - 1.), std::abs(ratio.low - 1.));
    const double resolution = params.positionSigma * ratio.high;
    return {.vertexLever = vertexLever, .variance = (scattering * scattering) + (resolution * resolution)};
  }
  const Interval zFrom = zRange<Det>(from);
  if (vertexZ.high >= zFrom.low && vertexZ.low <= zFrom.high) {
    return {.padding = secondMax<Det>(to) - secondMin<Det>(to)};
  }
  const double fromVertex = nearestApart(vertexZ, zFrom);
  const Interval zTo = zRange<Det>(to);
  const double gap = std::max(std::abs(zTo.high - zFrom.low), std::abs(zTo.low - zFrom.high));
  const double ratio = largestAbs(ratioRange(zTo, zFrom, vertexZ));
  const double vertexLever = std::abs(seedSecond) * gap / (fromVertex * fromVertex);
  const double resolution = params.positionSigma * ratio;
  return {.vertexLever = vertexLever, .variance = (scattering * scattering) + (resolution * resolution)};
}

/// Terms shared by the per-vertex interval of one source and layer pair.
/// Vertex-dependent terms are bounded once on the host over \a vertexZ, the
/// range the vertices an iteration uses can take; the device computes only the
/// centre and, for disc pairs, one division. \a sourceKink is as for
/// `secondScattering`.
template <DetectorLayout Det>
[[nodiscard]] inline VertexPairWindow vertexPairWindow(LayerID from, LayerID to, const SearchWindowParams& params, Interval vertexZ, double maxSigmaZ, double sourceKink, double betweenKink = 0.)
{
  const double centre = 0.5 * (vertexZ.low + vertexZ.high);
  const double reach = (0.5 * (vertexZ.high - vertexZ.low)) + (params.nSigma * maxSigmaZ);
  const VertexConstraint envelope{.z = centre, .sigmaZ = params.nSigma > 0. ? reach / params.nSigma : 0.};
  const double helix = helixSecondSpread<Det>(from, to, params, envelope);
  const double boundaryRounding = 4. * static_cast<double>(std::numeric_limits<Scalar>::epsilon()) * std::max(std::abs(secondMin<Det>(to)), std::abs(secondMax<Det>(to)));

  if (Detector<Det>::isBarrelLayer(from) != Detector<Det>::isBarrelLayer(to)) {
    const double scattering = secondScattering<Det>(from, to, params, vertexZ, sourceKink, betweenKink);
    const double padding = (params.nSigma * std::hypot(scattering, params.positionSigma)) + helix + boundaryRounding;
    return {.secondPadding = conservativeScalarBound(padding)};
  }

  const double seedSecond = secondMax<Det>(from);
  const auto uncertainty = vertexSecondWindow<Det>(from, to, seedSecond, params, vertexZ, sourceKink, betweenKink);
  const double centreSpread = std::max(secondCentreSpread<Det>(from, to, vertexZ.low), secondCentreSpread<Det>(from, to, vertexZ.high));
  double rounding = 0.;
  VertexPairWindow result;
  if (Detector<Det>::isBarrelLayer(from)) {
    const double scale = secondScale<Det>(from, to, 0.);
    result.secondScale = static_cast<Scalar>(scale);
    const double lever = farthestApart(vertexZ, secondRange<Det>(from));
    rounding = lever * std::abs(scale - static_cast<double>(result.secondScale));
    const double vertexLimit = largestAbs(vertexZ);
    rounding += 4. * static_cast<double>(std::numeric_limits<Scalar>::epsilon()) * (vertexLimit + (lever * std::abs(scale)));
  } else {
    const double source = discZ<Det>(from);
    const double target = discZ<Det>(to);
    result.sourcePosition = static_cast<Scalar>(source);
    result.targetPosition = static_cast<Scalar>(target);
    for (const double z : {vertexZ.low, vertexZ.high}) {
      for (const double seed : {secondMin<Det>(from), secondMax<Det>(from)}) {
        const auto deviceZ = static_cast<Scalar>(z);
        const auto deviceSeed = static_cast<Scalar>(seed);
        const Scalar device = deviceSeed * (result.targetPosition - deviceZ) / (result.sourcePosition - deviceZ);
        const double exact = static_cast<double>(deviceSeed) * (target - static_cast<double>(deviceZ)) / (source - static_cast<double>(deviceZ));
        rounding = std::max(rounding, std::abs(static_cast<double>(device) - exact));
      }
    }
    const double ratioLimit = largestAbs(ratioRange({.low = target, .high = target}, {.low = source, .high = source}, vertexZ));
    rounding += 8. * static_cast<double>(std::numeric_limits<Scalar>::epsilon()) * secondMax<Det>(from) * ratioLimit;
  }
  result.vertexLever = conservativeScalarBound(uncertainty.vertexLever);
  result.secondVariance = conservativeScalarBound(uncertainty.variance);
  result.secondPadding = conservativeScalarBound(uncertainty.padding + centreSpread + helix + rounding + boundaryRounding);
  return result;
}

/// The per-pair part of a barrel <-> disc window.
///
/// The exact longitudinal test puts the vertex `(0, z_v)` and both
/// measurements on one line in (r, z). Across kinds the target's coordinate
/// is the reciprocal of the seed's:
///
/// ```
/// barrel -> disc (at Z)    r_j = r_i (Z - z_v) / (z_i - z_v)
/// disc -> barrel (at R)    z_j = z_v + (z_i - z_v) R / r_i
/// ```
///
/// so no line in the seed's coordinate bounds it. Each is monotonic in `z_v`,
/// in the target's position and in the seed's coordinate separately, and the
/// query takes the extremes over the corners of those ranges for each seed.
/// What is fixed per pair is stored here: the target's z or radius range, the
/// vertex range, the seed's resolution, and a padding for scattering over the
/// longest lever arm between the two layers and for the target's resolution.
template <DetectorLayout Det>
[[nodiscard]] inline LayerPairWindow crossingPairWindow(LayerID from, LayerID to, const SearchWindowParams& params, const VertexConstraint& vertex, double sourceKink, double betweenKink = 0.)
{
  const bool toDisc = !Detector<Det>::isBarrelLayer(to);
  const Interval vertexZ = vertexReach(params, vertex);
  const double scattering = secondScattering<Det>(from, to, params, vertexZ, sourceKink, betweenKink);
  const double padding = (params.nSigma * std::hypot(scattering, params.positionSigma)) + helixSecondSpread<Det>(from, to, params, vertex);

  return {
    .phiHalfWidth = conservativeScalarBound(phiWindow<Det>(from, to, params, vertex, sourceKink, betweenKink)),
    .secondScale = 0.F,
    .secondOffset = static_cast<Scalar>(secondOffset<Det>(from, to, vertex.z)),
    .secondHalfWidth = conservativeScalarBound(padding),
    .kind = toDisc ? PairKind::kBarrelToDisc : PairKind::kDiscToBarrel,
    .targetLow = static_cast<Scalar>(toDisc ? zMin<Det>(to) : radiusMin<Det>(to)),
    .targetHigh = static_cast<Scalar>(toDisc ? zMax<Det>(to) : radiusMax<Det>(to)),
    .vertexLow = static_cast<Scalar>(vertexZ.low),
    .vertexHigh = static_cast<Scalar>(vertexZ.high),
    .seedSpread = conservativeScalarBound(params.nSigma * params.positionSigma),
    .targetSecondMin = static_cast<Scalar>(secondMin<Det>(to)),
    .targetSecondMax = static_cast<Scalar>(secondMax<Det>(to)),
  };
}

/// The half-width in phi of a displaced iteration (`SearchWindowParams::
/// impactMax`): the widest azimuth between
/// the two points that the displaced pair test's transverse part admits
/// (`displacedPairWithin`), scanned as `phiWindow` scans -- from half a turn
/// down, over the layers' radius ranges -- with the window's uniform errors:
/// `positionSigma` at both points, \a sourceKink and \a betweenKink as for
/// `phiWindow` (a source kink of `theta` is `S = (theta r_i)^2`), the beam
/// spot's width. The test is symmetric in the sign of the turn.
template <DetectorLayout Det>
[[nodiscard]] inline double displacedPhiWindow(LayerID from, LayerID to, const SearchWindowParams& params, const VertexConstraint& vertex, double sourceKink, double betweenKink = 0.)
{
  constexpr double kHalf = std::numbers::pi_v<double>;
  constexpr double kStep = 1.e-3;
  constexpr int kSteps = static_cast<int>(kHalf / kStep) + 1;
  constexpr int kSamples = 5;
  const double curvatureMax = maxCurvature(params);
  const double fromLow = radiusMin<Det>(from);
  const double fromHigh = radiusMax<Det>(from);
  const double toLow = radiusMin<Det>(to);
  const double toHigh = radiusMax<Det>(to);
  if (!(fromLow > 0.) || !(toLow > 0.)) {
    return kHalf;
  }
  displaced::PairErrors<double> errors;
  errors.sourceVarRPhi = errors.targetVarRPhi = params.positionSigma * params.positionSigma;
  errors.betweenScatter = betweenKink * betweenKink;
  errors.beamVarTransverse = vertex.sigmaTransverse * vertex.sigmaTransverse;
  constexpr double kInfinity = std::numeric_limits<double>::infinity();
  double widest = 0.;
  for (int i = 0; i < kSamples; ++i) {
    const double rFrom = fromLow + ((fromHigh - fromLow) * i / (kSamples - 1));
    errors.sourceScatter = (sourceKink * rFrom) * (sourceKink * rFrom);
    for (int j = 0; j < kSamples; ++j) {
      const double rTo = toLow + ((toHigh - toLow) * j / (kSamples - 1));
      for (int step = kSteps; step >= 0 && kStep * (step + 1) > widest; --step) {
        const double phi = kStep * step;
        const auto result = displaced::displacedPairWithin(rFrom, 0., 0., rTo * std::cos(phi), rTo * std::sin(phi), 0., errors, params.nSigma, curvatureMax, params.impactMax, -kInfinity, kInfinity);
        if (result.transverse) {
          widest = kStep * (step + 1);
          break;
        }
      }
    }
  }
  return std::min(widest, kHalf);
}

/// A displaced iteration's window (`displacedPhiWindow`): in the second
/// coordinate the target's whole extent, about its middle, as the closest
/// approach's z bound is the pair test's to apply, not a line's from the
/// beam spot.
template <DetectorLayout Det>
[[nodiscard]] inline LayerPairWindow displacedPairWindow(LayerID from, LayerID to, const SearchWindowParams& params, const VertexConstraint& vertex, double sourceKink, double betweenKink = 0.)
{
  const double low = secondMin<Det>(to);
  const double high = secondMax<Det>(to);
  const auto middle = static_cast<Scalar>(0.5 * (low + high));
  const double reach = std::max(std::abs(high - static_cast<double>(middle)), std::abs(static_cast<double>(middle) - low));
  return {
    .phiHalfWidth = conservativeScalarBound(displacedPhiWindow<Det>(from, to, params, vertex, sourceKink, betweenKink)),
    .secondScale = 0.F,
    .secondOffset = middle,
    .secondHalfWidth = conservativeScalarBound(reach + params.positionSigma),
    .kind = Detector<Det>::isBarrelLayer(to) ? PairKind::kBarrel : PairKind::kDisc,
  };
}

/// Assemble the per-pair window an iteration queries with.
///
/// \a seedSecond is only read for a disc-to-disc pair; the caller passes the
/// widest radius the seeding disc carries so that the stored half-width bounds
/// every seed on it. \a sourceKink and \a betweenKink are as for
/// `secondScattering`; the pair queries pass the pair cut's
/// (`layerPairWindow(..., table)` in `PairCut.h`).
template <DetectorLayout Det>
[[nodiscard]] inline LayerPairWindow layerPairWindow(LayerID from, LayerID to, double seedSecond, const SearchWindowParams& params, const VertexConstraint& vertex, double sourceKink, double betweenKink = 0.)
{
  if (params.impactMax > 0.) {
    return displacedPairWindow<Det>(from, to, params, vertex, sourceKink, betweenKink);
  }
  if (Detector<Det>::isBarrelLayer(from) != Detector<Det>::isBarrelLayer(to)) {
    return crossingPairWindow<Det>(from, to, params, vertex, sourceKink, betweenKink);
  }
  const double scale = secondScale<Det>(from, to, vertex.z);
  const double offset = secondOffset<Det>(from, to, vertex.z);
  const auto deviceScale = static_cast<Scalar>(scale);
  const auto deviceOffset = static_cast<Scalar>(offset);

  const double seedLimit = std::max(std::abs(secondMin<Det>(from)), std::abs(secondMax<Det>(from)));
  const double centreRounding = (seedLimit * std::abs(scale - static_cast<double>(deviceScale))) + std::abs(offset - static_cast<double>(deviceOffset));

  return {
    .phiHalfWidth = conservativeScalarBound(phiWindow<Det>(from, to, params, vertex, sourceKink, betweenKink)),
    .secondScale = deviceScale,
    .secondOffset = deviceOffset,
    .secondHalfWidth = conservativeScalarBound(secondWindow<Det>(from, to, seedSecond, params, vertex, sourceKink, betweenKink) + secondCentreSpread<Det>(from, to, vertex.z) + helixSecondSpread<Det>(from, to, params, vertex) + centreRounding),
    .kind = Detector<Det>::isBarrelLayer(from) ? PairKind::kBarrel : PairKind::kDisc,
  };
}

/// The same with the uniform kink, taking the seeding layer's own outer edge
/// as the conservative seed radius: what sizes the grid's cells, and a
/// whole-layer query.
template <DetectorLayout Det>
[[nodiscard]] inline LayerPairWindow layerPairWindow(LayerID from, LayerID to, const SearchWindowParams& params, const VertexConstraint& vertex)
{
  return layerPairWindow<Det>(from, to, secondMax<Det>(from), params, vertex, uniformSourceKink(params));
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SEARCHWINDOW_H_
