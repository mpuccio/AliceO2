// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file CandidateFinding.cxx
/// \brief Out-of-line tracklet, cell, and candidate-extension operations
///
/// Only this translation unit may include MFTFwdTrackHelpers.h on behalf of
/// the D007 operation boundary, so the common public headers
/// stays free of MFT-specific constants, TimeFrame, and typed-output
/// dependencies.

#include "ITSMFTTracking/detail/CandidateFinding.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>

#include "CommonConstants/MathConstants.h"
#include "DataFormatsITS/Vertex.h"
#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/detail/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Constants.h"
#include "ITStracking/MathUtils.h"
#include "ITStracking/TrackHelpers.h"
#include "MFTTracking/Constants.h"

namespace o2::itsmft::tracking
{

bool CylinderTrackletSearchWindow::acceptCandidate(
  const SurfaceMeasurement& sourceMeasurement,
  const o2::its::Cluster& sourceLocator,
  const SurfaceMeasurement& targetMeasurement,
  const o2::its::Cluster& targetLocator,
  float& tanLambdaOut) const
{
  const float deltaZ = o2::gpu::CAMath::Abs((tanLambda * (targetLocator.radius - sourceLocator.radius)) + sourceMeasurement.global.z - targetMeasurement.global.z);
  if (deltaZ / sigmaZ < nSigmaCut &&
      o2::its::math_utils::isPhiDifferenceBelow(sourceLocator.phi, targetLocator.phi, phiCut)) {
    const float acceptedTanLambda = (sourceMeasurement.global.z - targetMeasurement.global.z) / (sourceLocator.radius - targetLocator.radius);
    tanLambdaOut = acceptedTanLambda;
    return true;
  }
  return false;
}

bool DiskTrackletSearchWindow::acceptCandidate(
  const SurfaceMeasurement& sourceMeasurement,
  const SurfaceMeasurement& targetMeasurement,
  float& tanLambdaOut) const
{
  const float dx = targetMeasurement.global.x - xProj;
  const float dy = targetMeasurement.global.y - yProj;
  const float invSigmaX2 = (sigmaX > 0.f) ? 1.f / (sigmaX * sigmaX) : 0.f;
  const float invSigmaY2 = (sigmaY > 0.f) ? 1.f / (sigmaY * sigmaY) : 0.f;
  const float transChi2 = dx * dx * invSigmaX2 + dy * dy * invSigmaY2;
  const float sourceRadius = std::hypot(sourceMeasurement.global.x, sourceMeasurement.global.y);
  const float targetRadius = std::hypot(targetMeasurement.global.x, targetMeasurement.global.y);
  const float deltaR = sourceRadius - targetRadius;
  if (transChi2 < o2::its::math_utils::Sq(nSigmaCut) && std::abs(deltaR) > 1.e-6f) {
    const float acceptedTanLambda = (sourceMeasurement.global.z - targetMeasurement.global.z) / deltaR;
    tanLambdaOut = acceptedTanLambda;
    return true;
  }
  return false;
}

bool projectCylinderSearchWindow(const SurfaceMeasurement& sourceMeasurement,
                                 const o2::its::Cluster& sourceLocator,
                                 const o2::its::Vertex& vertex,
                                 const CylinderTrackletProjectionState& transitionState,
                                 float /*bz*/, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                                 const TrackingKernelParameters& params,
                                 CylinderTrackletSearchWindow& out)
{
  const float inverseR0 = 1.f / sourceLocator.radius;
  const float resolution = o2::gpu::CAMath::Sqrt(o2::its::math_utils::Sq(transitionState.sourcePositionResolution) +
                                                 o2::its::math_utils::Sq(params.pvResolution) / float(vertex.getNContributors()));
  const float tanLambda = (sourceMeasurement.global.z - vertex.getZ()) * inverseR0;
  const float zAtTargetMinR = tanLambda * (transitionState.targetMinR - sourceLocator.radius) + sourceMeasurement.global.z;
  const float zAtTargetMaxR = tanLambda * (transitionState.targetMaxR - sourceLocator.radius) + sourceMeasurement.global.z;
  const float sqInvDeltaZ0 = 1.f / (o2::its::math_utils::Sq(sourceMeasurement.global.z - vertex.getZ()) + o2::its::constants::Tolerance);
  const float sigmaZ = o2::gpu::CAMath::Sqrt((o2::its::math_utils::Sq(resolution) * o2::its::math_utils::Sq(tanLambda) *
                                              ((o2::its::math_utils::Sq(inverseR0) + sqInvDeltaZ0) * o2::its::math_utils::Sq(transitionState.meanDeltaR) + 1.f)) +
                                             o2::its::math_utils::Sq(transitionState.meanDeltaR * transitionState.transitionMSAngle));
  const auto bins = o2::itsmft::getBinsPhiZ(sourceLocator.phi, transitionState.toLayer,
                                            zAtTargetMinR, zAtTargetMaxR,
                                            sigmaZ * params.nSigmaCut, transitionState.transitionPhiCut,
                                            indexUtils);
  if (bins.x < 0) {
    return false;
  }
  out = {bins, tanLambda, sigmaZ, transitionState.transitionPhiCut, params.nSigmaCut};
  return true;
}

bool projectDiskSearchWindow(const SurfaceMeasurement& sourceMeasurement,
                             const o2::its::Cluster& sourceLocator,
                             const o2::its::Vertex& vertex,
                             const DiskTrackletProjectionState& transitionState,
                             float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
                             const TrackingKernelParameters& params,
                             DiskTrackletSearchWindow& out)
{
  float xProj = 0.f;
  float yProj = 0.f;
  detail::mftTrackletProject(sourceMeasurement.global.x, sourceMeasurement.global.y, sourceMeasurement.global.z,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             transitionState.fromZ, transitionState.toZ, bz, params.trackletMinPt,
                             xProj, yProj);
  float sigmaX = 0.f;
  float sigmaY = 0.f;
  detail::mftTrackletSigmaXY(sourceMeasurement.global.x, sourceMeasurement.global.y,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             sourceMeasurement.covariance.uu, sourceMeasurement.covariance.vv,
                             vertex.getSigmaX2(), vertex.getSigmaY2(), vertex.getSigmaZ2(),
                             transitionState.fromZ, transitionState.toZ,
                             transitionState.sourceReferenceRadius, transitionState.meanDeltaZ,
                             transitionState.transitionMSAngle, transitionState.transitionBendingAngle,
                             xProj, yProj, sigmaX, sigmaY);

  const float zSpread = params.nSigmaCut * vertex.getSigmaZ();
  const float zVtxMin = vertex.getZ() - zSpread;
  const float zVtxMax = vertex.getZ() + zSpread;
  const float absZFrom = std::abs(transitionState.fromZ);
  const float absZTo = std::abs(transitionState.toZ);
  const float denomMin = zVtxMax + absZFrom;
  const float denomMax = absZFrom + zVtxMin;
  float radialRangeMin = (std::abs(denomMin) > 1.e-6f) ? sourceLocator.radius * (zVtxMax + absZTo) / denomMin : sourceLocator.radius;
  float radialRangeMax = (std::abs(denomMax) > 1.e-6f) ? sourceLocator.radius * (absZTo + zVtxMin) / denomMax : sourceLocator.radius;
  if (radialRangeMin > radialRangeMax) {
    const float tmp = radialRangeMin;
    radialRangeMin = radialRangeMax;
    radialRangeMax = tmp;
  }
  const auto bins = o2::itsmft::getBinsRectClusterAtProj(xProj, yProj, transitionState.toLayer,
                                                         radialRangeMin, radialRangeMax,
                                                         sigmaX * params.nSigmaCut, sigmaY * params.nSigmaCut,
                                                         indexUtils);
  if (bins.x < 0) {
    return false;
  }
  out = {bins, xProj, yProj, sigmaX, sigmaY, params.nSigmaCut};
  return true;
}

bool bindTrackletProjectionState(
  SurfaceKind kind, int fromLayer, int toLayer,
  gsl::span<const float> layerRadii,
  gsl::span<const float> diskReferenceZ,
  float targetMinR, float targetMaxR,
  float sourcePositionResolution,
  float transitionMSAngle, float transitionBendingAngle,
  TrackletProjectionState& out) noexcept
{
  if (fromLayer < 0 || toLayer < 0 ||
      static_cast<size_t>(fromLayer) >= layerRadii.size() ||
      static_cast<size_t>(toLayer) >= layerRadii.size()) {
    return false;
  }
  switch (kind) {
    case SurfaceKind::Cylinder:
      out = CylinderTrackletProjectionState{fromLayer, toLayer,
                                            layerRadii[toLayer] - layerRadii[fromLayer],
                                            targetMinR, targetMaxR, sourcePositionResolution,
                                            transitionMSAngle, transitionBendingAngle};
      return true;
    case SurfaceKind::Disk:
      if (static_cast<size_t>(fromLayer) >= diskReferenceZ.size() ||
          static_cast<size_t>(toLayer) >= diskReferenceZ.size()) {
        return false;
      }
      out = DiskTrackletProjectionState{fromLayer, toLayer,
                                        diskReferenceZ[fromLayer], diskReferenceZ[toLayer],
                                        diskReferenceZ[toLayer] - diskReferenceZ[fromLayer],
                                        layerRadii[fromLayer], transitionMSAngle,
                                        transitionBendingAngle};
      return true;
  }
  return false;
}

bool projectTrackletSearchWindow(
  const SurfaceMeasurement& sourceMeasurement,
  const o2::its::Cluster& sourceLocator,
  const o2::its::Vertex& vertex,
  const TrackletProjectionState& transitionState,
  float bz, const o2::itsmft::IndexTableUtilsCore& indexUtils,
  const TrackingKernelParameters& params,
  TrackletSearchWindow& out)
{
  return std::visit(
    [&](const auto& state) {
      using State = std::decay_t<decltype(state)>;
      if constexpr (std::is_same_v<State, CylinderTrackletProjectionState>) {
        CylinderTrackletSearchWindow window{};
        if (!projectCylinderSearchWindow(sourceMeasurement, sourceLocator, vertex, state, bz, indexUtils, params, window)) {
          return false;
        }
        out = window;
      } else {
        DiskTrackletSearchWindow window{};
        if (!projectDiskSearchWindow(sourceMeasurement, sourceLocator, vertex, state, bz, indexUtils, params, window)) {
          return false;
        }
        out = window;
      }
      return true;
    },
    transitionState);
}

int4 trackletSearchBins(const TrackletSearchWindow& window) noexcept
{
  return std::visit([](const auto& value) { return value.bins; }, window);
}

int trackletSearchRowCount(const TrackletSearchWindow& window,
                           const o2::itsmft::IndexTableUtilsCore& indexUtils) noexcept
{
  return std::visit(
    [&](const auto& value) {
      int count = value.bins.w - value.bins.y + 1;
      if constexpr (std::is_same_v<std::decay_t<decltype(value)>, CylinderTrackletSearchWindow>) {
        if (count < 0) {
          count += indexUtils.getNrowBins();
        }
      }
      return std::max(0, count);
    },
    window);
}

int trackletSearchRowBin(const TrackletSearchWindow& window, int offset,
                         const o2::itsmft::IndexTableUtilsCore& indexUtils) noexcept
{
  return std::visit(
    [&](const auto& value) {
      int row = value.bins.y + offset;
      if constexpr (std::is_same_v<std::decay_t<decltype(value)>, CylinderTrackletSearchWindow>) {
        return row % indexUtils.getNrowBins();
      } else {
        return row < indexUtils.getNrowBins() ? row : -1;
      }
    },
    window);
}

bool acceptTrackletCandidate(
  const TrackletSearchWindow& window,
  const SurfaceMeasurement& sourceMeasurement,
  const o2::its::Cluster& sourceLocator,
  const SurfaceMeasurement& targetMeasurement,
  const o2::its::Cluster& targetLocator,
  float& tanLambdaOut) noexcept
{
  return std::visit(
    [&](const auto& value) {
      using Window = std::decay_t<decltype(value)>;
      if constexpr (std::is_same_v<Window, CylinderTrackletSearchWindow>) {
        return value.acceptCandidate(sourceMeasurement, sourceLocator, targetMeasurement, targetLocator, tanLambdaOut);
      } else {
        return value.acceptCandidate(sourceMeasurement, targetMeasurement, tanLambdaOut);
      }
    },
    window);
}

// --- Native cylinder/disk cell leaves ---
//
// The cell leaves are composed entirely from the existing
// barrel::/forward:: primitives (BarrelSurfaceStateOperations.h/
// ForwardSurfaceStateOperations.h) and the shared PID/absCharge-aware
// material kernel (MaterialPhysics.h). See CandidateFinding.h for
// the per-operation contract documentation.

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
  const double roundoffTolerance = 16. * std::numeric_limits<double>::epsilon() *
                                   std::max(diagonalProduct, covarianceSquared);
  return diagonalProduct - covarianceSquared >= -roundoffTolerance;
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

bool makeCylinderTransverseDirectionObservation(
  const SurfaceDescriptor& surface,
  const SurfaceMeasurement& measurement,
  TransverseDirectionObservation& observation) noexcept
{
  if (surface.kind != SurfaceKind::Cylinder ||
      !std::isfinite(surface.referenceCoordinate) || surface.referenceCoordinate <= 0.f ||
      !std::isfinite(measurement.frame.q) || !std::isfinite(measurement.frame.u) ||
      !std::isfinite(measurement.frame.frameAngle) ||
      !covarianceIsPositiveSemidefinite(measurement.covariance.uu,
                                        measurement.covariance.uv,
                                        measurement.covariance.vv)) {
    return false;
  }
  const double sine = std::sin(measurement.frame.frameAngle);
  const double cosine = std::cos(measurement.frame.frameAngle);
  const double varianceU = measurement.covariance.uu;
  const TransverseDirectionObservation scratch{
    measurement.frame.q * cosine - measurement.frame.u * sine,
    measurement.frame.q * sine + measurement.frame.u * cosine,
    sine * sine * varianceU,
    -sine * cosine * varianceU,
    cosine * cosine * varianceU};
  if (!observationIsValid(scratch)) {
    return false;
  }
  observation = scratch;
  return true;
}

bool makeDiskTransverseDirectionObservation(
  const SurfaceDescriptor& surface,
  const SurfaceMeasurement& measurement,
  TransverseDirectionObservation& observation) noexcept
{
  const TransverseDirectionObservation scratch{
    measurement.global.x, measurement.global.y,
    measurement.covariance.uu, measurement.covariance.uv,
    measurement.covariance.vv};
  if (surface.kind != SurfaceKind::Disk ||
      !std::isfinite(surface.referenceCoordinate) ||
      !observationIsValid(scratch)) {
    return false;
  }
  observation = scratch;
  return true;
}

bool makeTransverseDirectionObservation(
  const SurfaceDescriptor& surface,
  const SurfaceMeasurement& measurement,
  TransverseDirectionObservation& observation) noexcept
{
  using Builder = bool (*)(const SurfaceDescriptor&, const SurfaceMeasurement&,
                           TransverseDirectionObservation&) noexcept;
  static constexpr std::array<Builder, 2> builders{
    makeCylinderTransverseDirectionObservation,
    makeDiskTransverseDirectionObservation};
  const auto kindIndex = static_cast<std::size_t>(surface.kind);
  return kindIndex < builders.size() && builders[kindIndex](surface, measurement, observation);
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
  // Profile the unknown signed curvature over the physical TrackletMinPt
  // interval. Each asin is the chord's half central angle on one circle.
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
  // The middle-hit derivative retains the correlation between both chords;
  // the outgoing transition kick changes their relative azimuth directly.
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

bool makeCylinderDirectionObservation(const SurfaceDescriptor& surface,
                                      const SurfaceMeasurement& measurement,
                                      DirectionObservation& observation) noexcept
{
  const double q = measurement.frame.q;
  const double u = measurement.frame.u;
  const double radius2 = q * q + u * u;
  if (surface.kind != SurfaceKind::Cylinder ||
      !std::isfinite(surface.referenceCoordinate) || surface.referenceCoordinate <= 0.f ||
      !std::isfinite(q) || !std::isfinite(u) || radius2 <= 0. ||
      !std::isfinite(measurement.frame.v) ||
      !covarianceIsPositiveSemidefinite(measurement.covariance.uu,
                                        measurement.covariance.uv,
                                        measurement.covariance.vv)) {
    return false;
  }
  const double radius = std::sqrt(radius2);
  const double radialDerivativeU = u / radius;
  const DirectionObservation scratch{radius,
                                     measurement.frame.v,
                                     radialDerivativeU * radialDerivativeU * measurement.covariance.uu,
                                     radialDerivativeU * measurement.covariance.uv,
                                     measurement.covariance.vv};
  if (!observationIsValid(scratch)) {
    return false;
  }
  observation = scratch;
  return true;
}

bool makeDiskDirectionObservation(const SurfaceDescriptor& surface,
                                  const SurfaceMeasurement& measurement,
                                  DirectionObservation& observation) noexcept
{
  const double x = measurement.global.x;
  const double y = measurement.global.y;
  const double varianceX = measurement.covariance.uu;
  const double covarianceXY = measurement.covariance.uv;
  const double varianceY = measurement.covariance.vv;
  const double radius2 = x * x + y * y;
  if (surface.kind != SurfaceKind::Disk ||
      !std::isfinite(surface.referenceCoordinate) ||
      !std::isfinite(x) || !std::isfinite(y) || radius2 <= 0. ||
      !covarianceIsPositiveSemidefinite(varianceX, covarianceXY, varianceY)) {
    return false;
  }
  const double radius = std::sqrt(radius2);
  const double varianceR = (x * x * varianceX + 2. * x * y * covarianceXY +
                            y * y * varianceY) /
                           radius2;
  const DirectionObservation scratch{radius, surface.referenceCoordinate,
                                     varianceR, 0., 0.};
  if (!observationIsValid(scratch)) {
    return false;
  }
  observation = scratch;
  return true;
}

bool makeDirectionObservation(const SurfaceDescriptor& surface,
                              const SurfaceMeasurement& measurement,
                              DirectionObservation& observation) noexcept
{
  using Builder = bool (*)(const SurfaceDescriptor&, const SurfaceMeasurement&,
                           DirectionObservation&) noexcept;
  static constexpr std::array<Builder, 2> builders{
    makeCylinderDirectionObservation,
    makeDiskDirectionObservation};
  const auto kindIndex = static_cast<std::size_t>(surface.kind);
  return kindIndex < builders.size() && builders[kindIndex](surface, measurement, observation);
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
  // A middle-point kick displaces the outer point by L*n*dTheta, hence
  // dK/dTheta is the scalar product of the two segment vectors.
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
  if (!barrel::buildSeed(measurementInner, measurementMiddle, measurementOuter, bz, absCharge, pid, scratch, reason)) {
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
  switch (kind) {
    case SurfaceKind::Cylinder:
      return buildCylinderCellSeed(measurementInner, measurementMiddle, measurementOuter,
                                   material, bz, absCharge, pid, outState, chi2, params, reason);
    case SurfaceKind::Disk:
      return buildDiskCellSeed(measurementInner, measurementMiddle, measurementOuter,
                               material, bz, absCharge, pid, outState, chi2, params, reason);
  }
  reason = OperationFailureReason::SourceFamilyMismatch;
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

bool cellsCylinderAreCompatible(
  const SurfaceKinematicState& current,
  const SurfaceKinematicState& next,
  int /*currentSecondClusterIndex*/,
  int /*nextFirstClusterIndex*/,
  float bz,
  const TrackingKernelParameters& params) noexcept
{
  SurfaceKinematicState scratch = next;
  OperationFailureReason reason{};
  if (!barrel::rotate(scratch, current.alpha, reason) ||
      !barrel::propagate(scratch, current.referenceCoordinate, bz, reason)) {
    return false;
  }
  float chi2{0.f};
  if (!barrel::stateChi2(current, scratch, chi2, reason)) {
    return false;
  }
  return chi2 <= params.maxChi2ClusterAttachment;
}

bool cellsDiskAreCompatible(
  const SurfaceKinematicState& current,
  const SurfaceKinematicState& next,
  int currentSecondClusterIndex,
  int nextFirstClusterIndex,
  float bz,
  const TrackingKernelParameters& params) noexcept
{
  // Temporary Gate-3 compatibility input (see the header doc on the primary
  // template): checked first, exactly mirroring
  // detail::mftFwdCellsAreCompatible's own precedence.
  if (currentSecondClusterIndex != nextFirstClusterIndex) {
    return false;
  }
  SurfaceKinematicState scratch = next;
  OperationFailureReason reason{};
  if (!Propagator::propagateForward(scratch, current.referenceCoordinate, bz, reason)) {
    return false;
  }
  float chi2{0.f};
  if (!forward::stateChi2(current, scratch, chi2, reason)) {
    return false;
  }
  return chi2 <= params.maxChi2ClusterAttachment;
}

float cylinderLayerMultipleScatteringAngle(
  const CylinderLayerScatteringInputs& inputs, float trackletMinPt)
{
  // Keep the established ITS expression:
  // math_utils::MSangle(0.14f, trkParam.TrackletMinPt, trkParam.LayerxX0[iLayer]).
  return o2::its::math_utils::MSangle(0.14f, trackletMinPt, inputs.layerxX0);
}

float diskLayerMultipleScatteringAngle(
  const DiskLayerScatteringInputs& inputs, float trackletMinPt)
{
  // zLayer and rRef are supplied explicitly by the caller.
  const float invP = 1.f / trackletMinPt;
  const float zLayer = inputs.referenceCoordinate;
  const float rRef = inputs.layerRadius;
  const float tanlRef = (std::abs(rRef) > 1e-6f) ? zLayer / rRef : 0.f;
  const float absTanl = std::abs(tanlRef);
  const float cscLambda = (absTanl > 1e-6f) ? std::sqrt(1.f + tanlRef * tanlRef) / absTanl : 1e6f;
  return 0.0136f * invP * std::sqrt(inputs.layerxX0 * cscLambda);
}

namespace
{
// Static nominal half-disk z coordinates used by the disk scattering helper;
// compile-time storage gives returned spans process lifetime.
static constexpr std::array<float, o2::mft::constants::mft::LayersNumber> kLegacyMFTReferenceCoordinate =
  o2::mft::constants::mft::LayerZCoordinate();
} // namespace

DiskReferenceCoordinateView bindLegacyMFTReferenceCoordinates() noexcept
{
  return DiskReferenceCoordinateView{gsl::span<const float>(kLegacyMFTReferenceCoordinate)};
}

float clampCylinderTransitionCurvature(float oneOverR, float r2) noexcept
{
  // Preserves the legacy double-promoted comparison verbatim (frozen ITS
  // TimeFrame.cxx / common-CA non-MFT branch): `0.5` is a double literal.
  return (0.5 * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
}

float clampDiskTransitionCurvature(float oneOverR, float r2) noexcept
{
  // Preserves the legacy float-only comparison verbatim (common-CA MFT
  // branch): `0.5f` stays in float.
  return (0.5f * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
}

TransitionScatteringBendingPrep prepareTransitionScatteringAndBending(
  gsl::span<const float> perLayerMSAngle, int fromLayer, int toLayer,
  float r1, float r2, float clampedOneOverR, float res1, float res2) noexcept
{
  float ms2 = 0.f;
  for (int layer = fromLayer; layer < toLayer; ++layer) {
    ms2 += o2::its::math_utils::Sq(perLayerMSAngle[layer]);
  }
  const float msAngle = o2::gpu::CAMath::Sqrt(ms2);
  const float cosTheta1half = o2::gpu::CAMath::Sqrt(1.f - o2::its::math_utils::Sq(0.5f * r1 * clampedOneOverR));
  const float cosTheta2half = o2::gpu::CAMath::Sqrt(1.f - o2::its::math_utils::Sq(0.5f * r2 * clampedOneOverR));
  const float x = (r2 * cosTheta1half) - (r1 * cosTheta2half);
  const float delta = o2::gpu::CAMath::Sqrt(1.f / (1.f - 0.25f * o2::its::math_utils::Sq(x * clampedOneOverR)) *
                                            (o2::its::math_utils::Sq((0.25f * r1 * r2 * o2::its::math_utils::Sq(clampedOneOverR) / cosTheta2half) + cosTheta1half) * o2::its::math_utils::Sq(res1) +
                                             o2::its::math_utils::Sq((0.25f * r1 * r2 * o2::its::math_utils::Sq(clampedOneOverR) / cosTheta1half) + cosTheta2half) * o2::its::math_utils::Sq(res2)));
  const float phiCut = o2::gpu::CAMath::Min(o2::gpu::CAMath::ASin(0.5f * x * clampedOneOverR) + 2.f * msAngle + delta, o2::constants::math::PI * 0.5f);
  return TransitionScatteringBendingPrep{msAngle, phiCut};
}

} // namespace o2::itsmft::tracking
