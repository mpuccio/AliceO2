// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file CandidateFinding.cxx
/// \brief Tracklet, cell, and candidate-extension operations
///
/// This is the only translation unit that includes MFTFwdTrackHelpers.h for
/// D007, keeping public headers free of MFT-specific and typed-output dependencies.

#include "ITSMFTTracking/detail/CandidateFinding.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "CommonConstants/MathConstants.h"
#include "DataFormatsITS/Vertex.h"
#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/detail/DirectionCompatibility.h"
#include "ITSMFTTracking/detail/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Constants.h"
#include "ITStracking/MathUtils.h"
#include "ITStracking/TrackHelpers.h"

namespace o2::itsmft::tracking
{

namespace
{
bool acceptNormalizedResidual(const std::array<float, 2>& residual,
                              const float (&variance)[3], float nSigmaCut) noexcept
{
  const float determinant = variance[0] * variance[2] - o2::its::math_utils::Sq(variance[1]);
  if (!(determinant > 0.f)) {
    return false;
  }
  const float chi2 = (residual[0] * residual[0] * variance[2] -
                      2.f * residual[0] * residual[1] * variance[1] +
                      residual[1] * residual[1] * variance[0]) /
                     determinant;
  return std::isfinite(chi2) && chi2 < o2::its::math_utils::Sq(nSigmaCut);
}

bool cylinderTrackletCoordinates(const GlobalMeasurement& targetMeasurement,
                                 const o2::its::Cluster& targetLocator,
                                 std::array<float, 2>& coordinates) noexcept
{
  coordinates = {targetMeasurement.position.z, targetLocator.phi};
  return true;
}

bool diskTrackletCoordinates(const GlobalMeasurement& targetMeasurement,
                             std::array<float, 2>& coordinates) noexcept
{
  coordinates = {targetMeasurement.position.x, targetMeasurement.position.y};
  return true;
}

bool acceptTrackletCoordinates(const TrackletSearchWindow& window,
                               const std::array<float, 2>& coordinates,
                               bool periodicSecondCoordinate,
                               float nSigmaCut) noexcept
{
  std::array<float, 2> residual{window.prediction[0] - coordinates[0],
                                window.prediction[1] - coordinates[1]};
  if (periodicSecondCoordinate) {
    residual[1] = static_cast<float>(std::remainder(static_cast<double>(residual[1]),
                                                    static_cast<double>(o2::constants::math::TwoPI)));
  }
  return acceptNormalizedResidual(residual, window.variance, nSigmaCut);
}
} // namespace

bool projectCylinderSearchWindow(const GlobalMeasurement& sourceMeasurement,
                                 const o2::its::Cluster& sourceLocator,
                                 const o2::its::Vertex& vertex,
                                 const TrackletProjectionCache& transitionCache,
                                 float /*bz*/, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                                 const TrackingKernelParameters& params,
                                 TrackletSearchWindow& out)
{
  const float inverseR0 = 1.f / sourceLocator.radius;
  const float resolution = o2::gpu::CAMath::Sqrt(o2::its::math_utils::Sq(transitionCache.sourcePositionResolution) +
                                                 o2::its::math_utils::Sq(params.pvResolution) / float(vertex.getNContributors()));
  const float tanLambda = (sourceMeasurement.position.z - vertex.getZ()) * inverseR0;
  const float targetMeanRadius = 0.5f * (transitionCache.targetMinR + transitionCache.targetMaxR);
  const float zAtTargetMeanR = tanLambda * (targetMeanRadius - sourceLocator.radius) + sourceMeasurement.position.z;
  const float sqInvDeltaZ0 = 1.f / (o2::its::math_utils::Sq(sourceMeasurement.position.z - vertex.getZ()) + o2::its::constants::Tolerance);
  const float targetRadialVariance = o2::its::math_utils::Sq(transitionCache.targetMaxR - transitionCache.targetMinR) / 12.f;
  const float sigmaZ = o2::gpu::CAMath::Sqrt((o2::its::math_utils::Sq(resolution) * o2::its::math_utils::Sq(tanLambda) *
                                              ((o2::its::math_utils::Sq(inverseR0) + sqInvDeltaZ0) * o2::its::math_utils::Sq(transitionCache.toRadius - transitionCache.fromRadius) + 1.f)) +
                                             o2::its::math_utils::Sq((transitionCache.toRadius - transitionCache.fromRadius) * transitionCache.transitionMSAngle) +
                                             o2::its::math_utils::Sq(tanLambda) * targetRadialVariance);
  const auto bins = o2::itsmft::getBinsPhiZ(sourceLocator.phi, transitionCache.toLayer,
                                            zAtTargetMeanR, zAtTargetMeanR,
                                            sigmaZ * params.nSigmaCut, transitionCache.transitionPhiCut,
                                            indexUtils);
  if (bins.x < 0) {
    return false;
  }
  const float phiSigma = transitionCache.transitionPhiCut / params.nSigmaCut;
  out = {bins, {zAtTargetMeanR, sourceLocator.phi}, {o2::its::math_utils::Sq(sigmaZ), 0.f, o2::its::math_utils::Sq(phiSigma)}};
  return true;
}

bool projectDiskSearchWindow(const GlobalMeasurement& sourceMeasurement,
                             const o2::its::Cluster& sourceLocator,
                             const o2::its::Vertex& vertex,
                             const TrackletProjectionCache& transitionCache,
                             float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                             const TrackingKernelParameters& params,
                             TrackletSearchWindow& out)
{
  if (!transitionCache.hasReferenceCoordinates) {
    return false;
  }
  if (!(std::isfinite(transitionCache.targetMinZ) && std::isfinite(transitionCache.targetMaxZ) &&
        transitionCache.targetMinZ <= transitionCache.targetMaxZ)) {
    return false;
  }
  const float targetMeanZ = 0.5f * (transitionCache.targetMinZ + transitionCache.targetMaxZ);
  float xProj = 0.f;
  float yProj = 0.f;
  detail::mftTrackletProject(sourceMeasurement.position.x, sourceMeasurement.position.y, sourceMeasurement.position.z,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             transitionCache.fromReferenceCoordinate, targetMeanZ, bz, params.trackletMinPt,
                             xProj, yProj);
  float sigmaX = 0.f;
  float sigmaY = 0.f;
  detail::mftTrackletSigmaXY(sourceMeasurement.position.x, sourceMeasurement.position.y,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             sourceMeasurement.covariance.xx, sourceMeasurement.covariance.yy,
                             vertex.getSigmaX2(), vertex.getSigmaY2(), vertex.getSigmaZ2(),
                             transitionCache.fromReferenceCoordinate, targetMeanZ,
                             transitionCache.fromRadius, targetMeanZ - transitionCache.fromReferenceCoordinate,
                             transitionCache.transitionMSAngle, transitionCache.transitionPhiCut,
                             xProj, yProj, sigmaX, sigmaY);

  float xAtMinZ = 0.f;
  float yAtMinZ = 0.f;
  float xAtMaxZ = 0.f;
  float yAtMaxZ = 0.f;
  detail::mftTrackletProject(sourceMeasurement.position.x, sourceMeasurement.position.y, sourceMeasurement.position.z,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             transitionCache.fromReferenceCoordinate, transitionCache.targetMinZ, bz, params.trackletMinPt,
                             xAtMinZ, yAtMinZ);
  detail::mftTrackletProject(sourceMeasurement.position.x, sourceMeasurement.position.y, sourceMeasurement.position.z,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             transitionCache.fromReferenceCoordinate, transitionCache.targetMaxZ, bz, params.trackletMinPt,
                             xAtMaxZ, yAtMaxZ);
  const float targetZVarianceScale = 1.f / 12.f;
  const float deltaX = xAtMaxZ - xAtMinZ;
  const float deltaY = yAtMaxZ - yAtMinZ;
  const float varianceX = o2::its::math_utils::Sq(sigmaX) + targetZVarianceScale * o2::its::math_utils::Sq(deltaX);
  const float covarianceXY = targetZVarianceScale * deltaX * deltaY;
  const float varianceY = o2::its::math_utils::Sq(sigmaY) + targetZVarianceScale * o2::its::math_utils::Sq(deltaY);

  const float zSpread = params.nSigmaCut * vertex.getSigmaZ();
  const float zVtxMin = vertex.getZ() - zSpread;
  const float zVtxMax = vertex.getZ() + zSpread;
  const float absZFrom = std::abs(transitionCache.fromReferenceCoordinate);
  const float absZTo = std::abs(targetMeanZ);
  const float denomMin = zVtxMax + absZFrom;
  const float denomMax = absZFrom + zVtxMin;
  float radialRangeMin = (std::abs(denomMin) > 1.e-6f) ? sourceLocator.radius * (zVtxMax + absZTo) / denomMin : sourceLocator.radius;
  float radialRangeMax = (std::abs(denomMax) > 1.e-6f) ? sourceLocator.radius * (absZTo + zVtxMin) / denomMax : sourceLocator.radius;
  if (radialRangeMin > radialRangeMax) {
    const float tmp = radialRangeMin;
    radialRangeMin = radialRangeMax;
    radialRangeMax = tmp;
  }
  const auto bins = o2::itsmft::getBinsRectClusterAtProj(xProj, yProj, transitionCache.toLayer,
                                                         radialRangeMin, radialRangeMax,
                                                         std::sqrt(varianceX) * params.nSigmaCut, std::sqrt(varianceY) * params.nSigmaCut,
                                                         indexUtils);
  if (bins.x < 0) {
    return false;
  }
  out = {bins, {xProj, yProj}, {varianceX, covarianceXY, varianceY}};
  return true;
}

bool bindTrackletProjectionCache(
  int fromLayer, int toLayer,
  gsl::span<const float> layerRadii,
  gsl::span<const float> diskReferenceZ,
  float targetMinR, float targetMaxR, float targetMinZ, float targetMaxZ,
  float sourcePositionResolution,
  float transitionMSAngle, float transitionPhiCut,
  TrackletProjectionCache& out) noexcept
{
  if (fromLayer < 0 || toLayer < 0 ||
      static_cast<size_t>(fromLayer) >= layerRadii.size() ||
      static_cast<size_t>(toLayer) >= layerRadii.size()) {
    return false;
  }
  const bool hasReferenceCoordinates = static_cast<size_t>(fromLayer) < diskReferenceZ.size();
  const auto referenceCoordinate = [&](int layer) {
    return static_cast<size_t>(layer) < diskReferenceZ.size() ? diskReferenceZ[layer] : 0.f;
  };
  out = {fromLayer, toLayer,
         layerRadii[fromLayer], layerRadii[toLayer],
         targetMinR, targetMaxR, targetMinZ, targetMaxZ, sourcePositionResolution,
         referenceCoordinate(fromLayer),
         transitionMSAngle, transitionPhiCut, hasReferenceCoordinates};
  return true;
}

bool projectTrackletSearchWindow(
  const GlobalMeasurement& sourceMeasurement,
  const o2::its::Cluster& sourceLocator,
  const o2::its::Vertex& vertex,
  SurfaceKind kind,
  const TrackletProjectionCache& transitionCache,
  float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
  const TrackingKernelParameters& params,
  TrackletSearchWindow& out)
{
  switch (kind) {
    case SurfaceKind::Cylinder:
      return projectCylinderSearchWindow(sourceMeasurement, sourceLocator, vertex, transitionCache, bz, indexUtils, params, out);
    case SurfaceKind::Disk:
      return projectDiskSearchWindow(sourceMeasurement, sourceLocator, vertex, transitionCache, bz, indexUtils, params, out);
    case SurfaceKind::Undefined:
      return false;
  }
  return false;
}

bool acceptTrackletCandidate(
  const TrackletSearchWindow& window,
  const GlobalMeasurement& sourceMeasurement,
  const o2::its::Cluster& sourceLocator,
  const GlobalMeasurement& targetMeasurement,
  const o2::its::Cluster& targetLocator,
  SurfaceKind kind, float nSigmaCut, float& tanLambdaOut) noexcept
{
  std::array<float, 2> residual{};
  bool validCoordinates = false;
  switch (kind) {
    case SurfaceKind::Cylinder:
      validCoordinates = cylinderTrackletCoordinates(targetMeasurement, targetLocator, residual);
      if (!validCoordinates || !acceptTrackletCoordinates(window, residual, true, nSigmaCut)) {
        return false;
      }
      break;
    case SurfaceKind::Disk:
      validCoordinates = diskTrackletCoordinates(targetMeasurement, residual);
      if (!validCoordinates || !acceptTrackletCoordinates(window, residual, false, nSigmaCut)) {
        return false;
      }
      break;
    case SurfaceKind::Undefined:
      return false;
  }
  const float deltaR = sourceMeasurement.radius - targetMeasurement.radius;
  if (!(std::abs(deltaR) > 1.e-6f)) {
    return false;
  }
  tanLambdaOut = (sourceMeasurement.position.z - targetMeasurement.position.z) / deltaR;
  return true;
}

// Native cylinder/disk cell operations use the barrel/forward primitives and
// the shared PID/absCharge-aware material kernel; see CandidateFinding.h for contracts.

namespace
{

bool covarianceIsPositiveSemidefinite(double varianceFirst,
                                      double covariance,
                                      double varianceSecond) noexcept
{
  return std::isfinite(varianceFirst) && std::isfinite(covariance) &&
         std::isfinite(varianceSecond) && varianceFirst >= 0. &&
         varianceSecond >= 0. &&
         varianceFirst * varianceSecond - covariance * covariance >= 0.;
}

bool covarianceIsPositiveSemidefiniteWithinRoundoff(double varianceFirst,
                                                    double covariance,
                                                    double varianceSecond) noexcept
{
  if (!std::isfinite(varianceFirst) || !std::isfinite(covariance) ||
      !std::isfinite(varianceSecond) || varianceFirst < 0. || varianceSecond < 0.) {
    return false;
  }
  const double diagonalProduct = varianceFirst * varianceSecond;
  const double covarianceSquared = covariance * covariance;
  const double roundoffTolerance = 16. * std::numeric_limits<float>::epsilon() *
                                   std::max(diagonalProduct, covarianceSquared);
  return diagonalProduct - covarianceSquared >= -roundoffTolerance;
}

bool covarianceIsPositiveSemidefinite(const GlobalCovariance3F& covariance) noexcept
{
  const double xx = covariance.xx;
  const double xy = covariance.xy;
  const double xz = covariance.xz;
  const double yy = covariance.yy;
  const double yz = covariance.yz;
  const double zz = covariance.zz;
  if (!covarianceIsPositiveSemidefiniteWithinRoundoff(xx, xy, yy) ||
      !covarianceIsPositiveSemidefiniteWithinRoundoff(xx, xz, zz) ||
      !covarianceIsPositiveSemidefiniteWithinRoundoff(yy, yz, zz)) {
    return false;
  }
  const double determinant = xx * yy * zz + 2. * xy * xz * yz -
                             xx * yz * yz - yy * xz * xz - zz * xy * xy;
  const double scale = std::max({std::abs(xx * yy * zz),
                                 std::abs(2. * xy * xz * yz),
                                 std::abs(xx * yz * yz),
                                 std::abs(yy * xz * xz),
                                 std::abs(zz * xy * xy)});
  return std::isfinite(determinant) &&
         determinant >= -32. * std::numeric_limits<float>::epsilon() * scale;
}

bool sanitizeProjectedCovariance(double projectionScale,
                                 double& varianceFirst,
                                 double& covariance,
                                 double varianceSecond) noexcept
{
  if (!std::isfinite(projectionScale) || projectionScale < 0. ||
      !std::isfinite(varianceFirst) || !std::isfinite(covariance) ||
      !std::isfinite(varianceSecond) || varianceSecond < 0.) {
    return false;
  }
  const double varianceTolerance = 16. * std::numeric_limits<float>::epsilon() *
                                   std::max(projectionScale, std::abs(varianceFirst));
  if (varianceFirst < -varianceTolerance) {
    return false;
  }
  varianceFirst = std::max(0., varianceFirst);

  const double diagonalProduct = varianceFirst * varianceSecond;
  const double covarianceSquared = covariance * covariance;
  const double determinantTolerance = 16. * std::numeric_limits<float>::epsilon() *
                                      std::max(diagonalProduct, covarianceSquared);
  if (diagonalProduct - covarianceSquared < -determinantTolerance) {
    return false;
  }
  if (covarianceSquared > diagonalProduct) {
    covariance = diagonalProduct == 0. ? 0. : std::copysign(std::sqrt(diagonalProduct), covariance);
  }
  return true;
}

bool observationIsValid(const DirectionObservation& observation) noexcept
{
  return std::isfinite(observation.r) && observation.r > 0. &&
         std::isfinite(observation.z) &&
         covarianceIsPositiveSemidefinite(observation.varianceR,
                                          observation.covarianceRZ,
                                          observation.varianceZ);
}

bool observationIsValid(const TransverseDirectionObservation& observation) noexcept
{
  return std::isfinite(observation.x) && std::isfinite(observation.y) &&
         covarianceIsPositiveSemidefiniteWithinRoundoff(observation.varianceX,
                                                        observation.covarianceXY,
                                                        observation.varianceY);
}

} // namespace

bool makeTransverseDirectionObservation(
  const GlobalMeasurement& measurement,
  TransverseDirectionObservation& observation) noexcept
{
  if (!covarianceIsPositiveSemidefinite(measurement.covariance)) {
    return false;
  }
  const TransverseDirectionObservation scratch{
    measurement.position.x, measurement.position.y,
    measurement.covariance.xx, measurement.covariance.xy,
    measurement.covariance.yy};
  if (!observationIsValid(scratch)) {
    return false;
  }
  observation = scratch;
  return true;
}

bool trackletDirectionsAreTransverselyCompatible(
  const std::array<TransverseDirectionObservation, 3>& observations,
  float firstPhi, float secondPhi,
  const DirectionProcessNoise& processNoise,
  float bz, float trackletMinPt, float nSigmaCut,
  TransverseDirectionCompatibility& compatibility) noexcept
{
  if (!std::isfinite(firstPhi) || !std::isfinite(secondPhi) ||
      !std::isfinite(bz) || !std::isfinite(trackletMinPt) || trackletMinPt <= 0.f ||
      !std::isfinite(nSigmaCut) || nSigmaCut <= 0.f ||
      !std::isfinite(processNoise.angularVariance) || processNoise.angularVariance < 0. ||
      !observationIsValid(observations[0]) ||
      !observationIsValid(observations[1]) ||
      !observationIsValid(observations[2])) {
    return false;
  }

  const double deltaX01 = observations[0].x - observations[1].x;
  const double deltaY01 = observations[0].y - observations[1].y;
  const double deltaX12 = observations[1].x - observations[2].x;
  const double deltaY12 = observations[1].y - observations[2].y;
  const double lengthSquared01 = deltaX01 * deltaX01 + deltaY01 * deltaY01;
  const double lengthSquared12 = deltaX12 * deltaX12 + deltaY12 * deltaY12;
  if (!std::isfinite(lengthSquared01) || !std::isfinite(lengthSquared12) ||
      lengthSquared01 <= 0. || lengthSquared12 <= 0.) {
    return false;
  }
  const double length01 = std::sqrt(lengthSquared01);
  const double length12 = std::sqrt(lengthSquared12);
  // Profile signed curvature over the physical TrackletMinPt interval;
  // each asin is a chord's half central angle.
  const double maximumCurvature = std::min({std::abs(static_cast<double>(o2::constants::math::B2C) * bz) / trackletMinPt,
                                            2. / length01,
                                            2. / length12});
  const double maximumBending =
    std::asin(std::clamp(0.5 * maximumCurvature * length01, 0., 1.)) +
    std::asin(std::clamp(0.5 * maximumCurvature * length12, 0., 1.));
  const double deltaPhi = std::remainder(static_cast<double>(firstPhi) - secondPhi,
                                         2. * static_cast<double>(o2::constants::math::PI));

  const std::array<std::array<double, 2>, 3> derivatives{{
    {-deltaY01 / lengthSquared01, deltaX01 / lengthSquared01},
    {deltaY01 / lengthSquared01 + deltaY12 / lengthSquared12,
     -deltaX01 / lengthSquared01 - deltaX12 / lengthSquared12},
    {-deltaY12 / lengthSquared12, deltaX12 / lengthSquared12},
  }};
  // The middle-hit derivative correlates both chords; the outgoing kick
  // directly changes their relative azimuth.
  double variance = processNoise.angularVariance;
  for (std::size_t i = 0; i < observations.size(); ++i) {
    const double derivativeX = derivatives[i][0];
    const double derivativeY = derivatives[i][1];
    variance += derivativeX * derivativeX * observations[i].varianceX +
                2. * derivativeX * derivativeY * observations[i].covarianceXY +
                derivativeY * derivativeY * observations[i].varianceY;
  }
  const double excess = std::max(0., std::abs(deltaPhi) - maximumBending);
  if (!std::isfinite(deltaPhi) || !std::isfinite(maximumBending) ||
      !std::isfinite(variance) || variance <= 0.) {
    return false;
  }
  const double chi2 = excess * excess / variance;
  if (!std::isfinite(chi2)) {
    return false;
  }

  compatibility = TransverseDirectionCompatibility{deltaPhi, maximumBending, variance, chi2};
  return chi2 < static_cast<double>(nSigmaCut) * nSigmaCut;
}

bool makeDirectionObservation(const GlobalMeasurement& measurement,
                              DirectionObservation& observation) noexcept
{
  const double x = measurement.position.x;
  const double y = measurement.position.y;
  const double radius = measurement.radius;
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(radius) || radius <= 0. ||
      !covarianceIsPositiveSemidefinite(measurement.covariance)) {
    return false;
  }
  const double inverseRadius = 1. / radius;
  const double xxTerm = x * x * measurement.covariance.xx;
  const double xyTerm = 2. * x * y * measurement.covariance.xy;
  const double yyTerm = y * y * measurement.covariance.yy;
  const double inverseRadiusSquared = inverseRadius * inverseRadius;
  double varianceR = (xxTerm + xyTerm + yyTerm) * inverseRadiusSquared;
  double covarianceRZ = (x * measurement.covariance.xz +
                         y * measurement.covariance.yz) *
                        inverseRadius;
  const double projectionScale = (std::abs(xxTerm) + std::abs(xyTerm) + std::abs(yyTerm)) *
                                 inverseRadiusSquared;
  if (!sanitizeProjectedCovariance(projectionScale, varianceR, covarianceRZ,
                                   measurement.covariance.zz)) {
    return false;
  }
  const DirectionObservation scratch{radius, measurement.position.z,
                                     varianceR, covarianceRZ,
                                     measurement.covariance.zz};
  if (!observationIsValid(scratch)) {
    return false;
  }
  observation = scratch;
  return true;
}

bool cellDirectionsAreCompatible(const std::array<DirectionObservation, 3>& observations,
                                 const DirectionProcessNoise& processNoise,
                                 float nSigmaCut,
                                 CellDirectionCompatibility& compatibility) noexcept
{
  if (!std::isfinite(nSigmaCut) || nSigmaCut <= 0.f ||
      !std::isfinite(processNoise.angularVariance) || processNoise.angularVariance < 0. ||
      !observationIsValid(observations[0]) ||
      !observationIsValid(observations[1]) ||
      !observationIsValid(observations[2])) {
    return false;
  }

  const auto& first = observations[0];
  const auto& middle = observations[1];
  const auto& last = observations[2];
  const double deltaZ01 = first.z - middle.z;
  const double deltaZ12 = middle.z - last.z;
  const double deltaR01 = first.r - middle.r;
  const double deltaR12 = middle.r - last.r;
  const double residual = deltaZ01 * deltaR12 - deltaZ12 * deltaR01;
  const std::array<std::array<double, 2>, 3> derivatives{{
    {-deltaZ12, deltaR12},
    {deltaZ01 + deltaZ12, -(deltaR01 + deltaR12)},
    {-deltaZ01, deltaR01},
  }};
  double variance = 0.;
  for (std::size_t i = 0; i < observations.size(); ++i) {
    const double derivativeR = derivatives[i][0];
    const double derivativeZ = derivatives[i][1];
    variance += derivativeR * derivativeR * observations[i].varianceR +
                2. * derivativeR * derivativeZ * observations[i].covarianceRZ +
                derivativeZ * derivativeZ * observations[i].varianceZ;
  }
  const double incomingR = middle.r - first.r;
  const double incomingZ = middle.z - first.z;
  const double outgoingR = last.r - middle.r;
  const double outgoingZ = last.z - middle.z;
  const double segmentDotProduct = incomingR * outgoingR + incomingZ * outgoingZ;
  // A middle-point kick displaces the outer point by L*n*dTheta, so
  // dK/dTheta is the segment-vector dot product.
  variance += segmentDotProduct * segmentDotProduct * processNoise.angularVariance;
  if (!std::isfinite(residual) || !std::isfinite(variance) || variance <= 0.) {
    return false;
  }
  const double chi2 = residual * residual / variance;
  if (!std::isfinite(chi2)) {
    return false;
  }

  compatibility = CellDirectionCompatibility{residual, variance, chi2};
  const double threshold = static_cast<double>(nSigmaCut) * nSigmaCut;
  return chi2 < threshold;
}

bool buildCylinderCellSeed(
  const GlobalMeasurement& globalInner,
  const GlobalMeasurement& globalMiddle,
  const SurfaceMeasurement& measurementInner,
  const SurfaceMeasurement& measurementMiddle,
  const SurfaceMeasurement& measurementOuter,
  const std::array<NominalSurfaceMaterial, 3>& material,
  float bz,
  uint8_t absCharge,
  o2::track::PID pid,
  SurfaceKinematicState& outState,
  float& chi2,
  const TrackingKernelParameters& params,
  OperationFailureReason& reason) noexcept
{
  SurfaceKinematicState scratch{};
  if (!barrel::buildSeed(globalInner.position, globalMiddle.position, measurementOuter, bz, absCharge, pid, scratch, reason)) {
    return false;
  }

  float localChi2{0.f};
  const std::array<const SurfaceMeasurement*, 2> steps{&measurementMiddle, &measurementInner};
  const std::array<NominalSurfaceMaterial, 2> stepsMaterial{material[1], material[0]};
  for (int step = 0; step < 2; ++step) {
    const bool isLast = (step == 1);
    const auto& measurement = *steps[step];

    if (!barrel::rotate(scratch, measurement.frame.frameAngle, reason)) {
      return false;
    }
    if (!barrel::propagate(scratch, measurement.frame.q, bz, reason)) {
      return false;
    }
    const auto& stepMaterial = stepsMaterial[step];
    const auto materialResult = barrel::correctForMaterial(
      scratch, material::IntegratedMaterialBudget{stepMaterial.xOverX0, stepMaterial.arealDensityGPerCm2},
      material::MaterialTraversalDirection::OppositeMomentum);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    float predChi2{0.f};
    if (!barrel::predictedChi2(scratch, measurement, predChi2, reason)) {
      return false;
    }
    if (isLast && predChi2 > params.maxChi2ClusterAttachment) {
      reason = OperationFailureReason::PredictedChi2Failure;
      return false;
    }
    float updateChi2{0.f};
    if (!barrel::update(scratch, measurement, updateChi2, reason)) {
      return false;
    }
    localChi2 += updateChi2;
  }

  outState = scratch;
  chi2 = localChi2;
  return true;
}

bool buildDiskCellSeed(
  const SurfaceMeasurement& measurementInner,
  const SurfaceMeasurement& measurementMiddle,
  const SurfaceMeasurement& measurementOuter,
  const std::array<NominalSurfaceMaterial, 3>& material,
  float bz,
  uint8_t absCharge,
  o2::track::PID pid,
  SurfaceKinematicState& outState,
  float& chi2,
  const TrackingKernelParameters& params,
  OperationFailureReason& reason) noexcept
{
  SurfaceKinematicState scratch{};
  if (!forward::buildSeed(measurementInner, measurementMiddle, measurementOuter, bz, params.trackletMinPt,
                          absCharge, pid, scratch, reason)) {
    return false;
  }

  float localChi2{0.f};
  const std::array<const SurfaceMeasurement*, 3> steps{&measurementOuter, &measurementMiddle, &measurementInner};
  const std::array<NominalSurfaceMaterial, 3> stepsMaterial{material[2], material[1], material[0]};
  for (int step = 0; step < 3; ++step) {
    const bool isLast = (step == 2);
    const auto& measurement = *steps[step];

    if (!Propagator::propagateForward(scratch, measurement.frame.q, bz, reason)) {
      return false;
    }
    const auto& stepMaterial = stepsMaterial[step];
    const auto materialResult = forward::correctForMaterial(
      scratch, material::IntegratedMaterialBudget{stepMaterial.xOverX0, stepMaterial.arealDensityGPerCm2},
      material::MaterialTraversalDirection::OppositeMomentum);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    float predChi2{0.f};
    if (!forward::predictedChi2(scratch, measurement, predChi2, reason)) {
      return false;
    }
    if (isLast && predChi2 > params.maxChi2ClusterAttachment) {
      reason = OperationFailureReason::PredictedChi2Failure;
      return false;
    }
    float updateChi2{0.f};
    if (!forward::update(scratch, measurement, updateChi2, reason)) {
      return false;
    }
    localChi2 += updateChi2;
  }

  outState = scratch;
  chi2 = localChi2;
  return true;
}

bool buildCellSeed(
  SurfaceKind kind,
  const GlobalMeasurement& globalInner,
  const GlobalMeasurement& globalMiddle,
  const GlobalMeasurement& globalOuter,
  const SurfaceMeasurement& measurementInner,
  const SurfaceMeasurement& measurementMiddle,
  const SurfaceMeasurement& measurementOuter,
  const std::array<NominalSurfaceMaterial, 3>& material,
  float bz,
  uint8_t absCharge,
  o2::track::PID pid,
  SurfaceKinematicState& outState,
  float& chi2,
  const TrackingKernelParameters& params,
  OperationFailureReason& reason) noexcept
{
  // Generic orchestration resolves all global observations; disk seeding uses
  // only local measurements.
  (void)globalOuter;
  switch (kind) {
    case SurfaceKind::Undefined:
      reason = OperationFailureReason::SourceSurfaceKindMismatch;
      return false;
    case SurfaceKind::Cylinder:
      return buildCylinderCellSeed(globalInner, globalMiddle, measurementInner, measurementMiddle, measurementOuter,
                                   material, bz, absCharge, pid, outState, chi2, params, reason);
    case SurfaceKind::Disk:
      return buildDiskCellSeed(measurementInner, measurementMiddle, measurementOuter,
                               material, bz, absCharge, pid, outState, chi2, params, reason);
  }
  return false;
}

bool attachCylinderHit(
  SurfaceKinematicState& state,
  const SurfaceMeasurement& measurement,
  const NominalSurfaceMaterial& material,
  float bz,
  float& chi2,
  const TrackingKernelParameters& params,
  OperationFailureReason& reason) noexcept
{
  SurfaceKinematicState scratch = state;
  float scratchChi2 = chi2;

  if (!barrel::rotate(scratch, measurement.frame.frameAngle, reason)) {
    return false;
  }
  if (!barrel::propagate(scratch, measurement.frame.q, bz, reason)) {
    return false;
  }
  const auto materialResult = barrel::correctForMaterial(
    scratch, material::IntegratedMaterialBudget{material.xOverX0, material.arealDensityGPerCm2},
    material::MaterialTraversalDirection::OppositeMomentum);
  if (!materialResult.ok()) {
    reason = OperationFailureReason::MaterialFailure;
    return false;
  }
  float predChi2{0.f};
  if (!barrel::predictedChi2(scratch, measurement, predChi2, reason)) {
    return false;
  }
  if (predChi2 > params.maxChi2ClusterAttachment || predChi2 < 0.f) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  float updateChi2{0.f};
  if (!barrel::update(scratch, measurement, updateChi2, reason)) {
    return false;
  }
  scratchChi2 += updateChi2;

  state = scratch;
  chi2 = scratchChi2;
  return true;
}

bool attachDiskHit(
  SurfaceKinematicState& state,
  const SurfaceMeasurement& measurement,
  const NominalSurfaceMaterial& material,
  float bz,
  float& chi2,
  const TrackingKernelParameters& params,
  OperationFailureReason& reason) noexcept
{
  SurfaceKinematicState scratch = state;
  float scratchChi2 = chi2;

  if (!Propagator::propagateForward(scratch, measurement.frame.q, bz, reason)) {
    return false;
  }
  const auto materialResult = forward::correctForMaterial(
    scratch, material::IntegratedMaterialBudget{material.xOverX0, material.arealDensityGPerCm2},
    material::MaterialTraversalDirection::OppositeMomentum);
  if (!materialResult.ok()) {
    reason = OperationFailureReason::MaterialFailure;
    return false;
  }
  float predChi2{0.f};
  if (!forward::predictedChi2(scratch, measurement, predChi2, reason)) {
    return false;
  }
  if (predChi2 > params.maxChi2ClusterAttachment) {
    reason = OperationFailureReason::PredictedChi2Failure;
    return false;
  }
  float updateChi2{0.f};
  if (!forward::update(scratch, measurement, updateChi2, reason)) {
    return false;
  }
  scratchChi2 += updateChi2;

  state = scratch;
  chi2 = scratchChi2;
  return true;
}

} // namespace o2::itsmft::tracking
