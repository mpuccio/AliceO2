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

/// \file Seed.h
/// \brief Whole-chain seed helices and the diffuse priors for the uniform-solenoid track model.

#ifndef GTRK_TRACKING_SEED_H_
#define GTRK_TRACKING_SEED_H_

#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/Trig.h"

#include "GPUCommonMath.h"

#include <cstddef>
#include <numbers>
#include <utility>

namespace o2::gtrk
{

namespace detail
{

/// Global position of \a m relative to \a origin, in the requested arithmetic type.
template <typename T>
GPUdi() void relativePosition(const TrackingSurface& surface, const Measurement& m, const T origin[3], T out[3])
{
  Scalar global[3]{};
  surfaceToGlobal(surface, m.loc0, m.loc1, global);
  for (int i = 0; i < 3; ++i) {
    out[i] = static_cast<T>(global[i]) - origin[i];
  }
}

} // namespace detail

/// The q/p prior on a whole-chain seed: the wider of kChainQOverPSpread
/// times |q/p| and kChainPriorOverResolution times the resolution the chain's
/// measurements alone give.
///
/// Diffuse has to mean wide against what the chain determines, and a chain
/// determines q/p very differently with its length. Three points fix it to
/// tens of percent, and there a factor of a few is what diffuse means. Ten
/// layers over a metre fix it to a few percent where scattering limits it,
/// under one at high pT -- and a factor of three is then not diffuse but out
/// of range: an inward pass opens with a 20 cm step, the q/p uncertainty
/// becomes a position variance of hundreds of cm^2 against a 2.5e-6 cm^2
/// cluster, and the float covariance cannot come back from it. 0.3 keeps the
/// prior ten or more times wider than a scattering-limited resolution; the
/// second term keeps short chains as diffuse as a three-point seed. Without
/// any material ten points fix q/p to 1e-4, which no prior both
/// non-informative and representable in float reaches, but no track crosses
/// no material.
inline constexpr Scalar kChainQOverPSpread = 0.3F;
inline constexpr Scalar kChainPriorOverResolution = 30.F;

/// A direction's diffuse prior: a three-point seed fixes the direction to
/// milliradians; 0.1 rad is two orders looser.
inline constexpr Scalar kDiffuseAngleVariance = 1e-2F;

/// Set the deliberately diffuse diagonal covariance used by both fit passes.
///
/// The seed parameters come from the chain's own measurements, but their
/// covariance is not estimated from those same points. The measurements are
/// consumed by the filter with this broad prior instead.
///
/// **Diffuse means a hundred times the measurement, not a million**: a q/p
/// variance of 1e4 cannot be reduced to 2e-2 in float -- P - K H P returns
/// negative variances. So each prior is sized against what the measurements
/// determine, and q/p relative to the seed's own estimate: a fixed large width
/// stops being representable, a fixed small one stops being diffuse for a
/// 50 MeV track.
///
/// \a qOverPSigma is that width, and it belongs to where the estimate came
/// from: see kChainQOverPSpread.
GPUhdi() GPUAlwaysInline() void setDiffuseCovariance(BoundTrackState& seed, Scalar qOverPSigma)
{
  for (auto& value : seed.cov) {
    value = 0.F;
  }
  // A sensor is about a centimetre, so 1 mm is generous against a cluster
  // resolution of tens of microns.
  seed.variance(kLoc0) = 1e-2F;
  seed.variance(kLoc1) = 1e-2F;
  seed.variance(kPhi) = kDiffuseAngleVariance;
  seed.variance(kTheta) = kDiffuseAngleVariance;
  seed.variance(kQOverP) = o2::gpu::GPUCommonMath::Max(qOverPSigma * qOverPSigma, 1e-2F);
}

/// Free a state's direction, as a kink in the fit does:
/// `phi` and `theta` get the diffuse prior and lose their covariances with
/// every other parameter; position and q/p keep theirs.
GPUhdi() GPUAlwaysInline() void freeDirection(BoundTrackState& state)
{
  state.decorrelate(kPhi, kTheta);
  state.variance(kPhi) = kDiffuseAngleVariance;
  state.variance(kTheta) = kDiffuseAngleVariance;
}

/// A whole-chain seed and the q/p width of its prior, which both fit passes
/// need.
struct ChainSeed {
  BoundTrackState state;
  Scalar qOverPSigma = 0.F;
};

/// Construct a state on the first surface of an ordered chain from a helix
/// fitted through all \a n of its measurements.
///
/// This is the fit's outward seed, and its q/p width that of both passes'
/// priors. The filter linearises each step about its own filtered state
/// (TrackFitter.h, `filterStep`), so the seed only has to start it on the
/// right helix: a three-point seed's curvature, fixed over the first few
/// centimetres, is not that once the chain is long.
///
/// Transverse: Karimaki's circle fit (NIM A 305 (1991) 187) about the first
/// measurement. It is closed form -- sums over the points and one atan2 -- and
/// stays regular as the curvature vanishes, so a straight track needs no
/// branch. Its phi is defined modulo pi and is resolved towards the second
/// measurement; its curvature is positive clockwise, ours counter-clockwise.
/// Longitudinal: z against transverse arc length from the closest approach to
/// the first measurement, least squares in intercept and slope.
///
/// Both fits are unweighted: cluster-level accuracy is all a starting point
/// needs.
///
/// \a T controls the transient fit arithmetic. The returned state remains
/// Scalar, and the device path uses the Scalar default.
///
/// \a measurements is anything indexed by `[i]` giving a `Measurement`: an
/// array, or rows read from the frame's columns.
template <typename T = Scalar, typename Prop, DetectorLayout Det, typename Chain>
GPUdi() ChainSeed makeChainSeed(const Prop& prop, const TrackingGeometry<Det>& geometry, const Chain& measurements, int n)
{
  // Everything is relative to the first measurement: the circle is fitted
  // about it, and sums of large absolute coordinates would cancel.
  const auto& firstSurface = geometry.planeOf(measurements[0]);
  T origin[3]{};
  {
    Scalar global[3]{};
    surfaceToGlobal(firstSurface, measurements[0].loc0, measurements[0].loc1, global);
    for (int i = 0; i < 3; ++i) {
      origin[i] = static_cast<T>(global[i]);
    }
  }

  // Online centred moments avoid the cancellation in E[x^2] - E[x]^2.
  // The r2 covariances are rescaled after the pass so all quantities entering
  // q1 and q2 are of comparable size even for a metre-long chain.
  T meanX = T(0);
  T meanY = T(0);
  T meanR2 = T(0);
  T cXX = T(0);
  T cXY = T(0);
  T cYY = T(0);
  T cXR = T(0);
  T cYR = T(0);
  T cRR = T(0);
  T transverseScale = T(0);
  T second[3]{};
#ifndef GPUCA_GPUCODE
  // The longitudinal pass reads these same immutable points again. Cache
  // their first transform on the CPU; avoid a larger device register frame.
  constexpr auto kCachePoints = static_cast<std::size_t>(Detector<Det>::nLayers);
  T cachedPositions[kCachePoints][3];
  T cachedVarLoc0[kCachePoints];
  const bool cachePoints = n >= 0 && std::cmp_less_equal(n, kCachePoints);
#endif
  for (int i = 0; i < n; ++i) {
    T p[3]{};
#ifdef GPUCA_GPUCODE
    detail::relativePosition<T>(geometry.planeOf(measurements[i]), measurements[i], origin, p);
#else
    const auto measurement = measurements[i];
    detail::relativePosition<T>(geometry.planeOf(measurement), measurement, origin, p);
    if (cachePoints) {
      for (int axis = 0; axis < 3; ++axis) {
        cachedPositions[i][axis] = p[axis];
      }
      cachedVarLoc0[i] = static_cast<T>(measurement.varLoc0);
    }
#endif
    if (i == 1) {
      second[0] = p[0];
      second[1] = p[1];
    }
    const T r2 = p[0] * p[0] + p[1] * p[1];
    transverseScale = o2::gpu::GPUCommonMath::Max(transverseScale, o2::math_utils::detail::sqrt(r2));

    const T count = static_cast<T>(i + 1);
    const T dx = p[0] - meanX;
    const T dy = p[1] - meanY;
    const T dr = r2 - meanR2;
    meanX += dx / count;
    meanY += dy / count;
    meanR2 += dr / count;
    cXX += dx * (p[0] - meanX);
    cXY += dx * (p[1] - meanY);
    cYY += dy * (p[1] - meanY);
    cXR += dx * (r2 - meanR2);
    cYR += dy * (r2 - meanR2);
    cRR += dr * (r2 - meanR2);
  }

  const T invScale = T(1) / transverseScale;
  const T invScale2 = invScale * invScale;
  const T scaledMeanX = meanX * invScale;
  const T scaledMeanY = meanY * invScale;
  const T scaledMeanR2 = meanR2 * invScale2;
  const T scaledCXX = cXX * invScale2;
  const T scaledCXY = cXY * invScale2;
  const T scaledCYY = cYY * invScale2;
  const T scaledCXR = cXR * invScale2 * invScale;
  const T scaledCYR = cYR * invScale2 * invScale;
  const T scaledCRR = cRR * invScale2 * invScale2;

  const T q1 = scaledCRR * scaledCXY - scaledCXR * scaledCYR;
  const T q2 = scaledCRR * (scaledCXX - scaledCYY) - scaledCXR * scaledCXR + scaledCYR * scaledCYR;
  T phi = T(0.5) * atan2Of(T(2) * q1, q2);
  T sinPhi = T(0);
  T cosPhi = T(0);
  sincosOf(phi, sinPhi, cosPhi);
  const T scaledKappa = (sinPhi * scaledCXR - cosPhi * scaledCYR) / scaledCRR;
  const T delta = -scaledKappa * scaledMeanR2 + sinPhi * scaledMeanX - cosPhi * scaledMeanY;
  const T root = o2::math_utils::detail::sqrt(o2::gpu::GPUCommonMath::Max(T(0), T(1) - T(4) * delta * scaledKappa));
  T rho = T(2) * scaledKappa * invScale / root;
  T dca = T(2) * delta * transverseScale / (T(1) + root);
  if (cosPhi * second[0] + sinPhi * second[1] < T(0)) {
    sinPhi = -sinPhi;
    cosPhi = -cosPhi;
    phi += phi > T(0) ? -std::numbers::pi_v<T> : std::numbers::pi_v<T>;
    rho = -rho;
    dca = -dca;
  }
  const T curvature = -rho;

  // The closest approach to the first measurement, where the arc starts.
  const T start[2] = {dca * sinPhi, -dca * cosPhi};

  // Arc lengths accumulate point to point, so a chain that turns through more
  // than half a circle still unwinds monotonically.
  T arc = T(0);
  T previous[2] = {start[0], start[1]};
  T meanS = T(0);
  T meanZ = T(0);
  T cSS = T(0);
  T cSZ = T(0);
  T sumVarLoc0 = T(0);
  for (int i = 0; i < n; ++i) {
    T p[3]{};
#ifndef GPUCA_GPUCODE
    if (cachePoints) {
      sumVarLoc0 += cachedVarLoc0[i];
      for (int axis = 0; axis < 3; ++axis) {
        p[axis] = cachedPositions[i][axis];
      }
    } else
#endif
    {
      sumVarLoc0 += static_cast<T>(measurements[i].varLoc0);
      detail::relativePosition<T>(geometry.planeOf(measurements[i]), measurements[i], origin, p);
    }
    arc += detail::transverseArc<T>(curvature, p[0] - previous[0], p[1] - previous[1]);
    previous[0] = p[0];
    previous[1] = p[1];

    const T count = static_cast<T>(i + 1);
    const T ds = arc - meanS;
    const T dz = p[2] - meanZ;
    meanS += ds / count;
    meanZ += dz / count;
    cSS += ds * (arc - meanS);
    cSZ += ds * (p[2] - meanZ);
  }
  const T cotTheta = cSZ / cSS;
  const T zStart = meanZ - cotTheta * meanS;

  const T theta = atan2Of(T(1), cotTheta);
  const T thetaNorm = o2::math_utils::detail::sqrt(T(1) + cotTheta * cotTheta);
  const T sinTheta = T(1) / thetaNorm;
  const T cosTheta = cotTheta / thetaNorm;
  // The chain's propagator is asked at each use, not held in a local.
  const T field = static_cast<T>(stepPropagator(prop, geometry, measurements, 0, n - 1).bzKiloGauss) * static_cast<T>(constants::math::B2C);

  // The closest approach lies within a cluster residual of the first surface
  // but not on it. Step along the track onto the plane; over microns the
  // straight line and the helix agree to far below float precision.
  const T direction[3] = {sinTheta * cosPhi, sinTheta * sinPhi, cosTheta};
  T pointFromSurface[3]{};
  T height = T(0);
  T towards = T(0);
  for (int i = 0; i < 3; ++i) {
    pointFromSurface[i] = static_cast<T>(measurements[0].loc0) * static_cast<T>(firstSurface.axis0[i]) + static_cast<T>(measurements[0].loc1) * static_cast<T>(firstSurface.axis1[i]);
    if (i < 2) {
      pointFromSurface[i] += start[i];
    } else {
      pointFromSurface[i] += zStart;
    }
    height += pointFromSurface[i] * static_cast<T>(firstSurface.normal[i]);
    towards += direction[i] * static_cast<T>(firstSurface.normal[i]);
  }
  if (o2::gpu::GPUCommonMath::Abs(towards) > static_cast<T>(stepPropagator(prop, geometry, measurements, 0, n - 1).minAbsCosIncidence)) {
    for (int i = 0; i < 3; ++i) {
      pointFromSurface[i] -= height / towards * direction[i];
    }
  }
  T local[2]{};
  for (int i = 0; i < 3; ++i) {
    local[0] += pointFromSurface[i] * static_cast<T>(firstSurface.axis0[i]);
    local[1] += pointFromSurface[i] * static_cast<T>(firstSurface.axis1[i]);
  }

  ChainSeed seed;
  auto& state = seed.state;
  state.par[kLoc0] = static_cast<Scalar>(local[0]);
  state.par[kLoc1] = static_cast<Scalar>(local[1]);
  state.par[kPhi] = static_cast<Scalar>(phi);
  state.par[kTheta] = static_cast<Scalar>(theta);
  state.par[kQOverP] = static_cast<Scalar>(curvature * sinTheta / field);
  state.surface = measurements[0].surface;
  state.surfacePhi = firstSurface.surfacePhi;

  // The measurement-limited curvature resolution over the chain's own arc,
  // Gluckstern's equal-spacing form; scattering only widens it.
  const T inv = T(1) / static_cast<T>(n);
  const T curvatureSigma = o2::math_utils::detail::sqrt(sumVarLoc0 * inv) * o2::math_utils::detail::sqrt(T(720) / (static_cast<T>(n) + T(4))) / (arc * arc);
  seed.qOverPSigma = static_cast<Scalar>(o2::gpu::GPUCommonMath::Max(static_cast<T>(kChainQOverPSpread) * o2::gpu::GPUCommonMath::Abs(curvature), static_cast<T>(kChainPriorOverResolution) * curvatureSigma) * sinTheta / field);
  setDiffuseCovariance(state, seed.qOverPSigma);
  return seed;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SEED_H_
