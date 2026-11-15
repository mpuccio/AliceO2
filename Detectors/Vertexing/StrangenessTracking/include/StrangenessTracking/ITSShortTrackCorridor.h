// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. Distributed under the GNU GPL v3. CERN privileges and immunities are
// preserved.

#ifndef O2_STRANGENESS_TRACKING_ITSSHORTTRACKCORRIDOR_H
#define O2_STRANGENESS_TRACKING_ITSSHORTTRACKCORRIDOR_H

#include <algorithm>
#include <cmath>

namespace o2::strangeness_tracking::shorttrack
{
// Positions/radii/uncertainties in cm; angles in rad; curvature in 1/cm.
// sigmaXY bounds the Euclidean XY uncertainty (sqrt(trace(Cxy)) for a PV),
// sigmaZ is its longitudinal SD. Sums of SDs cover unknown correlations.
struct Point {
  double x = 0., y = 0., z = 0., sigmaXY = 0., sigmaZ = 0.;
};
struct Layer {
  double rMin = 0., rMax = 0., sigmaXY = 0., sigmaZ = 0., halfLength = 0.;
};
struct Settings {
  double maxAbsCurvature = 0., nSigma = 0.;
  // One-sigma projected scattering angles, including only material before the
  // outer and middle measurements respectively. The caller selects mass/pT.
  double scatteringAngle = 0., innerScatteringAngle = 0.;
};
struct Corridor {
  double phi = 0., z = 0., dPhi = 0., dZ = 0.;
  bool empty = false, fullLayer = false;
  double curvatureMin = 0., curvatureMax = 0.;
  bool contains(double hitPhi, double hitZ) const
  {
    return !empty &&
           std::abs(std::remainder(hitPhi - phi, 2. * M_PI)) <= dPhi &&
           std::abs(hitZ - z) <= dZ;
  }
};

namespace detail
{
struct Range {
  double lo, hi;
};
inline Range intersect(Range a, Range b)
{
  return {std::max(a.lo, b.lo), std::min(a.hi, b.hi)};
}
inline Range product(Range a, Range b)
{
  return {std::min({a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi}),
          std::max({a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi})};
}
inline Range square(Range a)
{
  return {a.lo <= 0. && a.hi >= 0. ? 0. : std::min(a.lo * a.lo, a.hi * a.hi),
          std::max(a.lo * a.lo, a.hi * a.hi)};
}
inline double maxAbs(Range a)
{
  return std::max(std::abs(a.lo), std::abs(a.hi));
}
inline double minAbs(Range a)
{
  return a.lo <= 0. && a.hi >= 0. ? 0.
                                  : std::min(std::abs(a.lo), std::abs(a.hi));
}
inline double arcCorrection(double x)
{
  return x == 0. ? 1. : std::asin(std::min(1., x)) / x;
}
inline Corridor full(const Layer& layer, Range curvature)
{
  return {0., 0., M_PI, layer.halfLength,
          false, true, curvature.lo, curvature.hi};
}
inline Corridor none(Range curvature)
{
  Corridor result;
  result.empty = true;
  result.curvatureMin = curvature.lo;
  result.curvatureMax = curvature.hi;
  return result;
}
// y(x) of the outgoing minor circle through (0,0) and (d,0). Rationalizing
// the two circle roots makes k=0 continuous without huge circle centres.
inline Range circleY(Range x, Range curvature, double d)
{
  const double maxK = maxAbs(curvature), minK = minAbs(curvature);
  const double furthest =
    std::max(std::abs(x.lo - .5 * d), std::abs(x.hi - .5 * d));
  const double nearest =
    x.lo <= .5 * d && x.hi >= .5 * d
      ? 0.
      : std::min(std::abs(x.lo - .5 * d), std::abs(x.hi - .5 * d));
  const double denominatorMin =
    std::sqrt(1. - .25 * maxK * maxK * d * d) +
    std::sqrt(1. - maxK * maxK * furthest * furthest);
  const double denominatorMax = std::sqrt(1. - .25 * minK * minK * d * d) +
                                std::sqrt(1. - minK * minK * nearest * nearest);
  const double f0 = x.lo * (d - x.lo), f1 = x.hi * (d - x.hi);
  const Range shape{std::min(f0, f1), x.lo <= .5 * d && x.hi >= .5 * d
                                        ? .25 * d * d
                                        : std::max(f0, f1)};
  return product(product({-curvature.hi, -curvature.lo}, shape),
                 {1. / denominatorMax, 1. / denominatorMin});
}

inline Corridor project(const Point& pv, const Point& outer, const Layer& layer,
                        const Settings& settings, Range curvature,
                        double scatteringAngle)
{
  if (curvature.lo > curvature.hi) {
    return none(curvature);
  }
  const double dx = outer.x - pv.x, dy = outer.y - pv.y, d = std::hypot(dx, dy);
  const double k = maxAbs(curvature), bend = .5 * k * d;
  // A PV outside the requested cylinder, a vanishing anchor or a near-turning
  // minor arc has no safe outgoing graph x(y). Preserve search support.
  if (!(d > 1.e-6) || !(layer.rMin > std::hypot(pv.x, pv.y)) || bend >= .95 ||
      !(layer.rMax >= layer.rMin)) {
    return full(layer, curvature);
  }
  const double ex = dx / d, ey = dy / d;
  const double pe = pv.x * ex + pv.y * ey, pn = -pv.x * ey + pv.y * ex;
  const double cosBend = std::sqrt(1. - bend * bend);
  const double arcLength = d * arcCorrection(bend);
  // Thin-scatterer bridge displacement is at most theta*arcLength/4.
  // Endpoint perturbations also rotate the chord; 2/cosBend bounds their
  // displacement on this non-turning circle graph. This is a search envelope,
  // not an additional hit/PV measurement or a fitted covariance.
  const double endpointError = settings.nSigma * (pv.sigmaXY + outer.sigmaXY);
  const double scatter = settings.nSigma * scatteringAngle * arcLength * .25;
  const double xyError =
    2. * endpointError / cosBend + settings.nSigma * layer.sigmaXY;
  const double rLo = std::max(0., layer.rMin - xyError),
               rHi = layer.rMax + xyError;
  if (!(rLo > std::hypot(pv.x, pv.y))) {
    return full(layer, curvature);
  }
  Range x{0., d}, y{0., 0.};
  // Interval contraction retains every curvature in the signed family and
  // every target radius. No sampling in curvature, phi or sensor radius.
  for (int iteration = 0; iteration < 8; ++iteration) {
    y = circleY(x, curvature, d);
    const auto radialY2 = square({y.lo + pn, y.hi + pn});
    if (radialY2.lo > rHi * rHi) {
      return none(curvature);
    }
    const Range radialX{std::sqrt(std::max(0., rLo * rLo - radialY2.hi)) - pe,
                        std::sqrt(std::max(0., rHi * rHi - radialY2.lo)) - pe};
    x = intersect(x, radialX);
    if (x.lo > x.hi) {
      return none(curvature);
    }
  }
  y = circleY(x, curvature, d);
  const double xm = .5 * (x.lo + x.hi), ym = .5 * (y.lo + y.hi);
  const double gx = pv.x + ex * xm - ey * ym, gy = pv.y + ey * xm + ex * ym;
  const double r = std::hypot(gx, gy);
  const double halfDiagonal =
    .5 * std::hypot(x.hi - x.lo, y.hi - y.lo) + xyError;
  // Scattering follows the fitter's Gaussian, linearized material model;
  // this is not a bound on arbitrary finite angular kicks. It is projected
  // at the fixed cylinder. Adding its full normal
  // displacement to the annulus as well would count the same uncertainty
  // twice and manufacture a full-phi search even for a radial crossing.
  const double maxSin =
    k * std::max(std::abs(x.lo - .5 * d), std::abs(x.hi - .5 * d));
  const auto tangentY = product(curvature, {x.lo - .5 * d, x.hi - .5 * d});
  const auto radialAlongX =
    product({pe + x.lo, pe + x.hi}, {std::sqrt(1. - maxSin * maxSin), 1.});
  const auto radialAlongY = product({pn + y.lo, pn + y.hi}, tangentY);
  // Endpoint uncertainty also perturbs the radial direction and tangent.
  const double incidence =
    std::min(1., (radialAlongX.lo + radialAlongY.lo) / rHi -
                   2. * endpointError / (rLo * cosBend));
  if (!(incidence > 1.e-3)) {
    return full(layer, curvature);
  }
  const double scatteringSpan = scatter / (layer.rMin * incidence);
  const double dPhi = halfDiagonal >= r || scatteringSpan >= 1.
                        ? M_PI
                        : std::min(M_PI, std::asin(halfDiagonal / r) +
                                           std::asin(scatteringSpan));
  const auto x2 = square(x), y2 = square(y);
  const double chordLo = std::sqrt(x2.lo + y2.lo),
               chordHi = std::min(d, std::sqrt(x2.hi + y2.hi));
  // z interpolates in transverse arc length, not cylindrical radius/chord x.
  const double minK = minAbs(curvature);
  const double fLo = std::clamp(
    chordLo / d * arcCorrection(.5 * minK * chordLo) / arcCorrection(bend),
    0., 1.);
  const double fHi = std::clamp(chordHi / d * arcCorrection(.5 * k * chordHi) /
                                  arcCorrection(.5 * minK * d),
                                0., 1.);
  const auto zOffset = product({fLo, fHi}, {outer.z - pv.z, outer.z - pv.z});
  const double slope = std::abs(outer.z - pv.z) / d;
  const double zError =
    settings.nSigma *
      ((1. - fLo) * pv.sigmaZ + fHi * outer.sigmaZ + layer.sigmaZ) +
    2. * endpointError * slope / cosBend +
    scatter * (1. + slope * slope +
               slope * std::sqrt(std::max(0., 1. - incidence * incidence)) /
                 incidence);
  const double zLo = pv.z + zOffset.lo - zError,
               zHi = pv.z + zOffset.hi + zError;
  // Rounding allowance protects exact post-bin gates at interval boundaries.
  constexpr double rounding = 1.e-6;
  return {std::atan2(gy, gx),
          .5 * (zLo + zHi),
          std::min(M_PI, dPhi + rounding),
          .5 * (zHi - zLo) + rounding,
          false,
          false,
          curvature.lo,
          curvature.hi};
}

inline Range inferCurvature(const Point& pv, const Point& outer,
                            const Point& middle, const Settings& settings)
{
  const Range floor{-settings.maxAbsCurvature, settings.maxAbsCurvature};
  const double dx = outer.x - pv.x, dy = outer.y - pv.y, d = std::hypot(dx, dy);
  const double r1 = std::hypot(middle.x - pv.x, middle.y - pv.y),
               r12 = std::hypot(outer.x - middle.x, outer.y - middle.y);
  const double endpointError = settings.nSigma * (pv.sigmaXY + outer.sigmaXY);
  const double bend = .5 * settings.maxAbsCurvature * d;
  if (!(d > endpointError) || bend >= .95) {
    return floor;
  }
  const double scatter = settings.nSigma * settings.scatteringAngle * d *
                         arcCorrection(bend) * .25;
  const double e1 = settings.nSigma * (pv.sigmaXY + middle.sigmaXY) + scatter;
  const double e12 =
    settings.nSigma * (middle.sigmaXY + outer.sigmaXY) + scatter;
  if (!(r1 > e1) || !(r12 > e12)) {
    return floor;
  }
  const double transverse =
    ((middle.y - pv.y) * dx - (middle.x - pv.x) * dy) / d;
  const double rotation = 2. * r1 * std::sin(.5 * std::asin(endpointError / d));
  const double error = e1 + rotation;
  // Signed three-point curvature. Interval products retain all combinations
  // of endpoint uncertainties; this is not an independent curvature fit.
  const Range numerator{-2. * (transverse + error), -2. * (transverse - error)};
  const Range inverseDenominator{1. / ((r1 + e1) * (r12 + e12)),
                                 1. / ((r1 - e1) * (r12 - e12))};
  return intersect(floor, product(numerator, inverseDenominator));
}
} // namespace detail

inline Corridor makeL1Corridor(const Point& pv, const Point& outer,
                               const Layer& layer, const Settings& settings)
{
  return detail::project(pv, outer, layer, settings,
                         {-settings.maxAbsCurvature, settings.maxAbsCurvature},
                         settings.scatteringAngle);
}

inline Corridor makeL0Corridor(const Point& pv, const Point& outer,
                               const Point& middle, const Layer& layer,
                               const Settings& settings)
{
  const auto curvature = detail::inferCurvature(pv, outer, middle, settings);
  auto result = detail::project(pv, outer, layer, settings, curvature,
                                settings.scatteringAngle);
  if (result.empty) {
    return result;
  }
  // Each bridge independently includes its endpoint uncertainties and the
  // intervening material. Intersect L1 first and the full-floor parent last:
  // interval-to-angular envelopes need not themselves nest when k contracts.
  const Corridor bridges[]{detail::project(pv, middle, layer, settings,
                                           curvature,
                                           settings.innerScatteringAngle),
                           makeL1Corridor(pv, outer, layer, settings)};
  for (const auto& bridge : bridges) {
    if (bridge.empty) {
      return detail::none(curvature);
    }
    if (result.dPhi >= M_PI) {
      result.phi = bridge.phi;
      result.dPhi = bridge.dPhi;
    } else if (bridge.dPhi < M_PI) {
      const double offset = std::remainder(result.phi - bridge.phi, 2. * M_PI);
      double phiLo = M_PI, phiHi = -M_PI;
      for (int turn = -1; turn <= 1; ++turn) {
        const double centre = offset + turn * 2. * M_PI;
        const double lo = std::max(-bridge.dPhi, centre - result.dPhi);
        const double hi = std::min(bridge.dPhi, centre + result.dPhi);
        if (lo <= hi) {
          phiLo = std::min(phiLo, lo);
          phiHi = std::max(phiHi, hi);
        }
      }
      if (phiLo > phiHi) {
        return detail::none(curvature);
      }
      // A disconnected circular overlap is conservatively enclosed WITHIN
      // this bridge's arc; no component is dropped at the +/-pi boundary.
      result.phi = std::remainder(bridge.phi + .5 * (phiLo + phiHi), 2. * M_PI);
      result.dPhi = .5 * (phiHi - phiLo);
    }
    const double zLo = std::max(result.z - result.dZ, bridge.z - bridge.dZ);
    const double zHi = std::min(result.z + result.dZ, bridge.z + bridge.dZ);
    if (zLo > zHi) {
      return detail::none(curvature);
    }
    result.fullLayer = result.fullLayer && bridge.fullLayer;
    result.z = .5 * (zLo + zHi);
    result.dZ = .5 * (zHi - zLo);
  }
  return result;
}
} // namespace o2::strangeness_tracking::shorttrack
#endif
