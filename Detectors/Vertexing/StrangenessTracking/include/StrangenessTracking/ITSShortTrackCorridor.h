// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. Distributed under the GNU GPL v3. CERN privileges and immunities are
// preserved.

#ifndef O2_STRANGENESS_TRACKING_ITSSHORTTRACKCORRIDOR_H
#define O2_STRANGENESS_TRACKING_ITSSHORTTRACKCORRIDOR_H

#include <algorithm>
#include <array>
#include <limits>
#include <vector>
#include <cmath>

namespace o2::strangeness_tracking::shorttrack
{
// Positions/radii/uncertainties in cm; angles in rad; curvature in 1/cm.
// sigmaXY bounds the Euclidean XY uncertainty (sqrt(trace(Cxy)) for a PV),
// sigmaZ is its longitudinal SD. Sums of SDs cover unknown correlations.
struct Point {
  double x = 0., y = 0., z = 0., sigmaXY = 0., sigmaZ = 0.;
  std::array<double, 6> covariance{}; // XX, XY, YY, XZ, YZ, ZZ; may be PSD/rank two
  bool hasCovariance = false;
};
struct Layer {
  double rMin = 0., rMax = 0., sigmaXY = 0., sigmaZ = 0., halfLength = 0.;
};
struct Settings {
  double maxAbsCurvature = 0., nSigma = 0.;
  // One-sigma projected scattering angles, including only material before the
  // outer and middle measurements respectively. The caller selects mass/pT.
  double scatteringAngle = 0., innerScatteringAngle = 0.;
  double curvatureToQPt = 0., mass = .13957039;
  std::array<double, 2> layerxX0{{.005, .005}};
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

} // namespace detail

inline Corridor makeL1Corridor(const Point& pv, const Point& outer,
                               const Layer& layer, const Settings& settings)
{
  return detail::project(pv, outer, layer, settings,
                         {-settings.maxAbsCurvature, settings.maxAbsCurvature},
                         settings.scatteringAngle);
}

struct ConditionalDiagnostics {
  bool valid = false, tailGuardUsed = false;
  int hypotheses = 0, profileComponents = 0;
  double nominalDPhi = 0., nominalDZ = 0.;
  double curvatureMean = 0., curvatureVariance = 0., qPt2 = 0., qPt4 = 0.;
  double thetaXY2L0 = 0., thetaXY2L1 = 0., positiveWeight = 0., negativeWeight = 0.;
};

namespace conditional_detail
{
struct Prediction {
  double phi = 0., z = 0., s = 0., length = 0., incidence = 0., tangent = 0.;
  bool valid = false;
};
// Exact minor-circle/cylinder intersection in the displaced PV chord frame.
// Multiplying the circle equation by k avoids a centre at infinity at k=0.
inline Prediction predict(const Point& pv, const Point& outer, double k, double radius)
{
  Prediction out;
  const double dx = outer.x - pv.x, dy = outer.y - pv.y, d = std::hypot(dx, dy);
  if (!(d > 1.e-6) || !(radius > std::hypot(pv.x, pv.y)) || !(std::abs(.5 * k * d) < .95)) {
    return out;
  }
  const double ex = dx / d, ey = dy / d, bend = .5 * k * d;
  const double cs = std::sqrt(1. - bend * bend);
  const double vx = .5 * k * (pv.x + outer.x) - cs * ey;
  const double vy = .5 * k * (pv.y + outer.y) + cs * ex;
  const double rhs = .5 * k * (radius * radius + pv.x * outer.x + pv.y * outer.y) + cs * (-pv.x * ey + pv.y * ex);
  const double argument = rhs / (radius * std::hypot(vx, vy));
  if (!(std::abs(argument) <= 1.)) {
    return out;
  }
  const double direction = std::atan2(vy, vx), offset = std::acos(argument);
  const double length = d * detail::arcCorrection(std::abs(bend));
  for (int sign : {-1, 1}) {
    const double phi = direction + sign * offset;
    const double x = radius * std::cos(phi) - pv.x, y = radius * std::sin(phi) - pv.y;
    const double along = x * ex + y * ey, transverse = -x * ey + y * ex;
    const double chord = std::hypot(x, y);
    if (along < -1.e-7 || along > d + 1.e-7 || k * transverse > 1.e-7 || chord > d + 1.e-7) {
      continue;
    }
    const double fraction = chord / d * detail::arcCorrection(.5 * std::abs(k) * chord) / detail::arcCorrection(std::abs(bend));
    const double tangent = std::atan2(dy, dx) + 2. * std::asin(bend) * (fraction - .5);
    const double incidence = std::cos(tangent - phi);
    if (!(incidence > 1.e-3)) {
      continue;
    }
    return {std::remainder(phi, 2. * M_PI), pv.z + fraction * (outer.z - pv.z), fraction * length, length, incidence, tangent, true};
  }
  return out;
}
using Matrix = std::array<std::array<double, 2>, 2>;
struct Joint {
  Matrix zero{}, one{}, cross{};
  bool valid = false;
};
inline std::array<double, 6> covariance(const Point& p)
{
  if (p.hasCovariance) {
    return p.covariance;
  }
  // 3*diag dominates any PSD 3D covariance with these diagonal bounds.
  return {3. * p.sigmaXY * p.sigmaXY, 0., 3. * p.sigmaXY * p.sigmaXY, 0., 0., 3. * p.sigmaZ * p.sigmaZ};
}
inline double covElement(const std::array<double, 6>& c, int i, int j)
{
  constexpr int index[3][3]{{0, 1, 3}, {1, 2, 4}, {3, 4, 5}};
  return c[index[i][j]];
}
inline Joint endpointCovariance(const Point& pv, const Point& outer, double k, double r0, double r1)
{
  double jac[2][2][6]{};
  for (int parameter = 0; parameter < 6; ++parameter) {
    Point lowPV = pv, highPV = pv, lowOuter = outer, highOuter = outer;
    Point& low = parameter < 3 ? lowPV : lowOuter;
    Point& high = parameter < 3 ? highPV : highOuter;
    double* lowCoordinate[] = {&low.x, &low.y, &low.z};
    double* highCoordinate[] = {&high.x, &high.y, &high.z};
    const int coordinate = parameter % 3;
    const double step = 1.e-4 * (1. + std::abs(*lowCoordinate[coordinate]));
    *lowCoordinate[coordinate] -= step;
    *highCoordinate[coordinate] += step;
    for (int layer = 0; layer < 2; ++layer) {
      const double radius = layer == 0 ? r0 : r1;
      const auto lowPrediction = predict(lowPV, lowOuter, k, radius);
      const auto highPrediction = predict(highPV, highOuter, k, radius);
      if (!lowPrediction.valid || !highPrediction.valid) {
        return {};
      }
      jac[layer][0][parameter] = radius * std::remainder(highPrediction.phi - lowPrediction.phi, 2. * M_PI) / (2. * step);
      jac[layer][1][parameter] = (highPrediction.z - lowPrediction.z) / (2. * step);
    }
  }
  Joint result;
  const auto cp = covariance(pv), co = covariance(outer);
  for (int row = 0; row < 2; ++row) {
    for (int col = 0; col < 2; ++col) {
      for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
          if ((i < 3) != (j < 3)) {
            continue;
          }
          const double c = covElement(i < 3 ? cp : co, i % 3, j % 3);
          result.zero[row][col] += jac[0][row][i] * c * jac[0][col][j];
          result.one[row][col] += jac[1][row][i] * c * jac[1][col][j];
          result.cross[row][col] += jac[0][row][i] * c * jac[1][col][j];
        }
      }
    }
  }
  result.valid = true;
  return result;
}
// Residuals are compared at the MEASURED cylinder radius. Radial hit
// errors therefore move both the measurement and the predicted crossing.
inline Matrix measurementCovariance(const Point& p, const Prediction& prediction, double dip)
{
  const auto c = covariance(p);
  const double phi = std::atan2(p.y, p.x), sn = std::sin(phi), cs = std::cos(phi);
  const double tanIncidence = std::sin(prediction.tangent - prediction.phi) / prediction.incidence;
  const double jacobian[2][3]{{-sn - tanIncidence * cs, cs - tanIncidence * sn, 0.},
                              {-dip * cs / prediction.incidence, -dip * sn / prediction.incidence, 1.}};
  Matrix result{};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
          result[i][j] += jacobian[i][row] * covElement(c, row, col) * jacobian[j][col];
        }
      }
    }
  }
  return result;
}
inline double bridge(double s, double u, double length) { return std::min(s, u) - s * u / length; }
// Projected azimuth-kick variance at EACH signed-curvature hypothesis. The
// pT-only beta and dip bound are conservative approximations, not a plug-in
// momentum from the measured sagitta. The cylinder incidence approximates the
// nominal material path; no energy loss or configured material LUT is modeled.
inline double materialVariance(double qPt2, double dipBound, double incidence, double xx0, double mass)
{
  return .0136 * .0136 * qPt2 * (1. + mass * mass * qPt2) * xx0 * std::sqrt(1. + dipBound * dipBound) / incidence;
}
inline void addProcess(Joint& joint, const Prediction& zero, const Prediction& one, double variance0, double variance1, double dipBound, double dip, double kick0Location)
{
  const Prediction positions[]{zero, one};
  double a[2], b[2];
  for (int i = 0; i < 2; ++i) {
    a[i] = 1. / positions[i].incidence;
    b[i] = -dip * std::sin(positions[i].tangent - positions[i].phi) / positions[i].incidence;
  }
  for (int kick = 0; kick < 2; ++kick) {
    const double location = kick == 0 ? kick0Location : one.s;
    const double variance = kick == 0 ? variance0 : variance1;
    const double g0 = bridge(zero.s, location, zero.length), g1 = bridge(one.s, location, one.length);
    const double green[]{g0, g1};
    Matrix* blocks[]{&joint.zero, &joint.one, &joint.cross};
    constexpr int left[]{0, 1, 0}, right[]{0, 1, 1};
    for (int block = 0; block < 3; ++block) {
      const int i = left[block], j = right[block];
      const double scale = green[i] * green[j] * variance;
      auto& c = *blocks[block];
      c[0][0] += scale * a[i] * a[j];
      c[0][1] += scale * a[i] * b[j];
      c[1][0] += scale * b[i] * a[j];
      c[1][1] += scale * (b[i] * b[j] + 1. + dipBound * dipBound);
    }
  }
}
struct Hypothesis {
  double k = 0., logWeight = 0., chi2 = 0., qPt2 = 0., qPt4 = 0., theta0 = 0., theta1 = 0.;
  std::array<double, 3> y{}, z{}, varY{}, varZ{};
  bool valid = false;
};
inline Hypothesis evaluateHypothesis(double k, double integrationWeight, const Point& pv, const Point& outer, const Point& middle,
                                     const Layer& layer, const Settings& settings, const std::array<Joint, 3>& endpoints,
                                     double dipBound, double phiReference)
{
  Hypothesis out;
  const double r1 = std::hypot(middle.x, middle.y);
  const auto one = predict(pv, outer, k, r1);
  if (!one.valid) {
    return out;
  }
  const double residual[]{r1 * std::remainder(std::atan2(middle.y, middle.x) - one.phi, 2. * M_PI), middle.z - one.z};
  const double qPt2 = k * k * settings.curvatureToQPt * settings.curvatureToQPt;
  const double dip = (outer.z - pv.z) / one.length;
  const auto rm = measurementCovariance(middle, one, dip);
  const double radius[]{layer.rMin, .5 * (layer.rMin + layer.rMax), layer.rMax};
  const auto nominalZero = predict(pv, outer, k, radius[1]);
  if (!nominalZero.valid) {
    return out;
  }
  double determinant = 0., chi2 = 0.;
  for (int radial = 0; radial < 3; ++radial) {
    const auto zero = predict(pv, outer, k, radius[radial]);
    if (!zero.valid) {
      return out;
    }
    const double theta0 = materialVariance(qPt2, dipBound, nominalZero.incidence, settings.layerxX0[0], settings.mass);
    const double theta1 = materialVariance(qPt2, dipBound, one.incidence, settings.layerxX0[1], settings.mass);
    auto joint = endpoints[radial];
    addProcess(joint, zero, one, theta0, theta1, dipBound, dip, nominalZero.s);
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        joint.one[i][j] += rm[i][j];
      }
    }
    const auto& v = joint.one;
    determinant = v[0][0] * v[1][1] - v[0][1] * v[1][0];
    if (!(v[0][0] > 0.) || !(v[1][1] > 0.) || !(determinant > 0.)) {
      return out;
    }
    const Matrix inverse{{{v[1][1] / determinant, -v[0][1] / determinant}, {-v[1][0] / determinant, v[0][0] / determinant}}};
    Matrix gain{};
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        gain[i][j] = joint.cross[i][0] * inverse[0][j] + joint.cross[i][1] * inverse[1][j];
      }
    }
    out.y[radial] = radius[radial] * std::remainder(zero.phi - phiReference, 2. * M_PI) + gain[0][0] * residual[0] + gain[0][1] * residual[1];
    out.z[radial] = zero.z + gain[1][0] * residual[0] + gain[1][1] * residual[1];
    // Target measurement noise enters ONCE, after conditioning the latent
    // position. 3*diag bounds unknown 3D target correlations and is projected
    // at its measured radius, just like the actual middle measurement.
    out.varY[radial] = std::max(0., joint.zero[0][0] - gain[0][0] * joint.cross[0][0] - gain[0][1] * joint.cross[0][1]) + 3. * layer.sigmaXY * layer.sigmaXY / (zero.incidence * zero.incidence);
    out.varZ[radial] = std::max(0., joint.zero[1][1] - gain[1][0] * joint.cross[1][0] - gain[1][1] * joint.cross[1][1]) + 3. * (layer.sigmaZ * layer.sigmaZ + dip * dip * layer.sigmaXY * layer.sigmaXY / (zero.incidence * zero.incidence));
    if (radial == 1) {
      chi2 = residual[0] * (inverse[0][0] * residual[0] + inverse[0][1] * residual[1]) + residual[1] * (inverse[1][0] * residual[0] + inverse[1][1] * residual[1]);
      out.chi2 = chi2;
      out.theta0 = theta0;
      out.theta1 = theta1;
      out.logWeight = std::log(integrationWeight) - .5 * (chi2 + std::log(determinant));
    }
  }
  out.k = k;
  out.qPt2 = qPt2;
  out.qPt4 = qPt2 * qPt2;
  out.valid = true;
  return out;
}
} // namespace conditional_detail
inline Corridor makeL0Corridor(const Point& pv, const Point& outer, const Point& middle, const Layer& layer,
                               const Settings& settings, ConditionalDiagnostics* diagnostic = nullptr, int quadratureRefinement = 1)
{
  using namespace conditional_detail;
  if (diagnostic) {
    *diagnostic = {};
  }
  const auto fallback = makeL1Corridor(pv, outer, layer, settings);
  const double d = std::hypot(outer.x - pv.x, outer.y - pv.y), maximum = settings.maxAbsCurvature;
  const double endpointError = settings.nSigma * (pv.sigmaXY + outer.sigmaXY);
  if (!(maximum > 0.) || !(std::abs(settings.curvatureToQPt) > 0.) || !(d > endpointError) ||
      !(layer.rMin > std::hypot(pv.x, pv.y) + endpointError) || .5 * maximum * d >= .95) {
    return fallback;
  }
  const double distance1 = std::hypot(middle.x - pv.x, middle.y - pv.y), distance12 = std::hypot(outer.x - middle.x, outer.y - middle.y);
  if (!(distance1 > 1.e-6) || !(distance12 > 1.e-6)) {
    return fallback;
  }
  const double transverse = ((middle.y - pv.y) * (outer.x - pv.x) - (middle.x - pv.x) * (outer.y - pv.y)) / d;
  const double measuredCurvature = std::clamp(-2. * transverse / (distance1 * distance12), -maximum, maximum);
  const double sigmaK = std::max(1.e-10 * maximum, 2. * (middle.sigmaXY + pv.sigmaXY + outer.sigmaXY) / (distance1 * distance12));
  const double radius[]{layer.rMin, .5 * (layer.rMin + layer.rMax), layer.rMax};
  const double r1 = std::hypot(middle.x, middle.y);
  std::array<Joint, 3> endpoints;
  for (int radial = 0; radial < 3; ++radial) {
    // Local-EKF endpoint Jacobians are linearized once at the measured
    // three-anchor seed. Every curvature node still has its OWN process noise.
    endpoints[radial] = endpointCovariance(pv, outer, measuredCurvature, radius[radial], r1);
    if (!endpoints[radial].valid) {
      return fallback;
    }
  }
  const double dipBound = (std::abs(outer.z - pv.z) + settings.nSigma * (pv.sigmaZ + outer.sigmaZ)) / (d - endpointError);
  const double phiReference = std::atan2(outer.y - pv.y, outer.x - pv.x);
  // Flat signed-curvature search prior, not a physical pT spectrum. A log
  // mesh resolves soft modes; measurement-scale cells resolve arbitrarily
  // narrow near-straight likelihoods. GL4 cell weights integrate d(curvature).
  std::vector<double> boundaries{-maximum, 0., maximum};
  for (int exponent = 1; exponent <= 20; ++exponent) {
    boundaries.push_back(std::ldexp(maximum, -exponent));
    boundaries.push_back(-std::ldexp(maximum, -exponent));
  }
  for (double scale : {.25, .5, 1., 2., 4., 8., 16.}) {
    for (int sign : {-1, 1}) {
      boundaries.push_back(std::clamp(measuredCurvature + sign * scale * sigmaK, -maximum, maximum));
      boundaries.push_back(std::clamp(sign * scale * sigmaK, -maximum, maximum));
    }
  }
  std::sort(boundaries.begin(), boundaries.end());
  boundaries.erase(std::unique(boundaries.begin(), boundaries.end()), boundaries.end());
  constexpr double nodes[]{-.8611363115940526, -.3399810435848563, .3399810435848563, .8611363115940526};
  constexpr double weights[]{.3478548451374539, .6521451548625461, .6521451548625461, .3478548451374539};
  std::vector<Hypothesis> hypotheses;
  double maximumLogWeight = std::numeric_limits<double>::lowest();
  const int subdivisions = std::max(1, quadratureRefinement);
  for (size_t cell = 1; cell < boundaries.size(); ++cell) {
    const double width = (boundaries[cell] - boundaries[cell - 1]) / subdivisions;
    if (!(width > 0.)) {
      continue;
    }
    for (int division = 0; division < subdivisions; ++division) {
      const double centre = boundaries[cell - 1] + (division + .5) * width;
      for (int node = 0; node < 4; ++node) {
        const double curvature = centre + .5 * width * nodes[node];
        auto hypothesis = evaluateHypothesis(curvature, .5 * width * weights[node], pv, outer, middle, layer, settings, endpoints, dipBound, phiReference);
        if (!hypothesis.valid) {
          return fallback;
        }
        maximumLogWeight = std::max(maximumLogWeight, hypothesis.logWeight);
        hypotheses.push_back(hypothesis);
      }
    }
  }
  if (hypotheses.empty()) {
    return fallback;
  }
  double normalizer = 0., meanK = 0., secondK = 0., q2 = 0., q4 = 0., theta0 = 0., theta1 = 0., positive = 0., negative = 0.;
  std::array<double, 3> meanPhi{}, meanZ{}, secondPhi{}, secondZ{};
  for (const auto& hypothesis : hypotheses) {
    // Log-sum-exp protects both a narrow straight mode and broad soft tails.
    const double weight = std::exp(hypothesis.logWeight - maximumLogWeight);
    normalizer += weight;
    meanK += weight * hypothesis.k;
    secondK += weight * hypothesis.k * hypothesis.k;
    q2 += weight * hypothesis.qPt2;
    q4 += weight * hypothesis.qPt4;
    theta0 += weight * hypothesis.theta0;
    theta1 += weight * hypothesis.theta1;
    positive += hypothesis.k > 0. ? weight : 0.;
    negative += hypothesis.k < 0. ? weight : 0.;
    for (int radial = 0; radial < 3; ++radial) {
      const double phi = hypothesis.y[radial] / radius[radial];
      meanPhi[radial] += weight * phi;
      meanZ[radial] += weight * hypothesis.z[radial];
      secondPhi[radial] += weight * (hypothesis.varY[radial] / (radius[radial] * radius[radial]) + phi * phi);
      secondZ[radial] += weight * (hypothesis.varZ[radial] + hypothesis.z[radial] * hypothesis.z[radial]);
    }
  }
  if (!(normalizer > 0.)) {
    return fallback;
  }
  double phiLo = std::numeric_limits<double>::max(), phiHi = std::numeric_limits<double>::lowest();
  double zLo = std::numeric_limits<double>::max(), zHi = std::numeric_limits<double>::lowest();
  for (int radial = 0; radial < 3; ++radial) {
    const double phi = meanPhi[radial] / normalizer, z = meanZ[radial] / normalizer;
    const double sigmaPhi = std::sqrt(std::max(0., secondPhi[radial] / normalizer - phi * phi));
    const double sigmaZ = std::sqrt(std::max(0., secondZ[radial] / normalizer - z * z));
    phiLo = std::min(phiLo, phi - settings.nSigma * sigmaPhi);
    phiHi = std::max(phiHi, phi + settings.nSigma * sigmaPhi);
    zLo = std::min(zLo, z - settings.nSigma * sigmaZ);
    zHi = std::max(zHi, z + settings.nSigma * sigmaZ);
  }
  const double nominalPhiLo = phiLo, nominalPhiHi = phiHi, nominalZLo = zLo, nominalZHi = zHi;
  int profileComponents = 0;
  auto retainProfile = [&](const Hypothesis& hypothesis) {
    // A low posterior density is not a reason to erase a physically compatible
    // soft branch. Use the existing spatial nSigma compatibility for the L1
    // residual; the determinant/flat-prior weights still define the moments.
    if (!hypothesis.valid || hypothesis.chi2 > settings.nSigma * settings.nSigma) {
      return;
    }
    ++profileComponents;
    for (int radial = 0; radial < 3; ++radial) {
      const double phi = hypothesis.y[radial] / radius[radial];
      const double halfPhi = settings.nSigma * std::sqrt(hypothesis.varY[radial]) / radius[radial];
      const double halfZ = settings.nSigma * std::sqrt(hypothesis.varZ[radial]);
      phiLo = std::min(phiLo, phi - halfPhi);
      phiHi = std::max(phiHi, phi + halfPhi);
      zLo = std::min(zLo, hypothesis.z[radial] - halfZ);
      zHi = std::max(zHi, hypothesis.z[radial] + halfZ);
    }
  };
  for (const auto& hypothesis : hypotheses) {
    retainProfile(hypothesis);
  }
  // GL nodes do not include the prior endpoints or exactly k=0. Include those
  // as support-only profiles, without giving point masses to the flat prior.
  for (double curvature : {-maximum, 0., maximum}) {
    retainProfile(evaluateHypothesis(curvature, 1., pv, outer, middle, layer, settings, endpoints, dipBound, phiReference));
  }
  const bool tailGuardUsed = phiLo < nominalPhiLo || phiHi > nominalPhiHi || zLo < nominalZLo || zHi > nominalZHi;
  // Radial shell support: hull the three conditional radial predictions and
  // add a nominal-circle derivative guard between the shell samples. Process
  // propagation and covariance are local Gaussian approximations; this is not
  // an enclosure of arbitrary finite kicks or a Gaussian-tail guarantee.
  const double bend = .5 * maximum * d, cosine = std::sqrt(1. - bend * bend);
  const double sagitta = .5 * d * bend / (1. + cosine);
  const double transverseBound = (std::hypot(pv.x, pv.y) + endpointError + sagitta) / layer.rMin;
  if (!(transverseBound < 1.)) {
    return fallback;
  }
  const double angle = std::asin(bend) + std::asin(transverseBound);
  if (!(angle < .5 * M_PI)) {
    return fallback;
  }
  const double gap = .25 * (layer.rMax - layer.rMin);
  const double radialPhi = gap * std::tan(angle) / layer.rMin;
  const double radialZ = gap * dipBound / std::cos(angle);
  constexpr double rounding = 1.e-6;
  phiLo -= radialPhi + rounding;
  phiHi += radialPhi + rounding;
  zLo -= radialZ + rounding;
  zHi += radialZ + rounding;
  if (diagnostic) {
    diagnostic->valid = true;
    diagnostic->tailGuardUsed = tailGuardUsed;
    diagnostic->profileComponents = profileComponents;
    diagnostic->nominalDPhi = std::min(M_PI, .5 * (nominalPhiHi - nominalPhiLo) + radialPhi + rounding);
    diagnostic->nominalDZ = .5 * (nominalZHi - nominalZLo) + radialZ + rounding;
    diagnostic->hypotheses = hypotheses.size();
    diagnostic->curvatureMean = meanK / normalizer;
    diagnostic->curvatureVariance = std::max(0., secondK / normalizer - diagnostic->curvatureMean * diagnostic->curvatureMean);
    diagnostic->qPt2 = q2 / normalizer;
    diagnostic->qPt4 = q4 / normalizer;
    diagnostic->thetaXY2L0 = theta0 / normalizer;
    diagnostic->thetaXY2L1 = theta1 / normalizer;
    diagnostic->positiveWeight = positive / normalizer;
    diagnostic->negativeWeight = negative / normalizer;
  }
  // Return the hull of the mixture window AND separately compatible signed
  // conditional profiles. Neither a mixture RMS nor this Gaussian process
  // approximation establishes an all-tail or arbitrary finite-kick guarantee.
  return {std::remainder(phiReference + .5 * (phiLo + phiHi), 2. * M_PI), .5 * (zLo + zHi),
          std::min(M_PI, .5 * (phiHi - phiLo)), .5 * (zHi - zLo), false, false, -maximum, maximum};
}

} // namespace o2::strangeness_tracking::shorttrack
#endif
