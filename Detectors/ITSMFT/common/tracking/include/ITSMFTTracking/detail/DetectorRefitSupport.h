// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORREFITSUPPORT_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORREFITSUPPORT_H_

#ifndef GPUCA_GPUCODE

#include <cmath>

#include "DataFormatsCalibration/MeanVertexObject.h"
#include "ITSMFTTracking/GenericTrack.h"
#include "ITSMFTTracking/RefitDriver.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking::detail
{

inline void configureITSBeamPosition(TimeFrame& frame,
                                     const TrackingParameters& params,
                                     const o2::dataformats::MeanVertexObject* meanVertex,
                                     bool overrideBeamEstimation)
{
  const float systErrY2 = params.SystError2Row.empty() ? 0.f : params.SystError2Row[0];
  const float layerRes = params.LayerResolution.empty() ? 0.f : params.LayerResolution[0];
  if (overrideBeamEstimation && meanVertex != nullptr) {
    frame.setBeamPosition(meanVertex->getX(), meanVertex->getY(), meanVertex->getSigmaY2(), layerRes, systErrY2);
  } else if (params.UseDiamond) {
    frame.setBeamPosition(params.Diamond[0], params.Diamond[1], params.DiamondCov[3], layerRes, systErrY2);
  }
}

inline void configureMFTBeamPosition(TimeFrame& frame, const TrackingParameters& params)
{
  const float systErrY2 = params.SystError2Row.empty() ? 0.f : params.SystError2Row[0];
  const float layerRes = params.LayerResolution.empty() ? 0.f : params.LayerResolution[0];
  frame.setBeamPosition(params.Diamond[0], params.Diamond[1], params.DiamondCov[3], layerRes, systErrY2);
}

inline bool fillCandidateKinematics(TrackingCandidate& candidate) noexcept
{
  auto& state = candidate.track.innerState;
  if (!state.hasRecognizedKind() || !std::isfinite(state.parameters[4])) {
    return false;
  }
  candidate.charge = state.parameters[4] < 0.f ? -1 : 1;
  if (state.kind == SurfaceKind::Cylinder) {
    if (std::abs(state.parameters[2]) >= 1.f) {
      return false;
    }
    candidate.phi = std::asin(state.parameters[2]) + state.alpha;
    candidate.eta = std::asinh(state.parameters[3]);
  } else {
    candidate.phi = state.parameters[2];
    candidate.eta = std::asinh(state.parameters[3]);
  }
  return std::isfinite(candidate.phi) && std::isfinite(candidate.eta);
}

inline bool refitSurfaceSeed(const TrackSeed& seed,
                             const TrackingParameters& params,
                             float bz,
                             TimeFrameScratch& scratch,
                             gsl::span<const gsl::span<const GlobalMeasurement>> layerGlobals,
                             gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                             SurfaceCatalogView surfaceCatalog,
                             gsl::span<const SurfaceId> orderedSurfaces,
                             TrackingCandidate& candidate)
{
  if (layerGlobals.size() != layerMeasurements.size() || orderedSurfaces.size() != layerMeasurements.size()) {
    return false;
  }
  bool hasMeasurement = false;
  for (int position = 0; position < static_cast<int>(layerMeasurements.size()); ++position) {
    const int clusterIndex = seed.getCluster(position);
    if (clusterIndex == o2::its::constants::UnusedIndex) {
      continue;
    }
    if (clusterIndex < 0 || clusterIndex >= static_cast<int>(layerMeasurements[position].size()) ||
        clusterIndex >= static_cast<int>(layerGlobals[position].size())) {
      return false;
    }
    const auto& global = layerGlobals[position][clusterIndex];
    const auto& measurement = layerMeasurements[position][clusterIndex];
    const auto expectedSource = scratch.getSurfaceSource(position);
    if (!scratch.hasClusterExternalIndex(position, clusterIndex)) {
      return false;
    }
    const int externalIndex = scratch.getClusterExternalIndex(position, clusterIndex);
    if (!expectedSource || global.surface != orderedSurfaces[position] || !surfaceCatalog.hasSurface(global.surface) ||
        !global.cluster.isValid() || global.cluster.source != *expectedSource || externalIndex < 0 ||
        global.cluster.index != static_cast<uint32_t>(externalIndex) ||
        !std::isfinite(measurement.frame.u) || !std::isfinite(measurement.frame.v) || !std::isfinite(measurement.frame.q) ||
        !std::isfinite(measurement.covariance.uu) || !std::isfinite(measurement.covariance.uv) || !std::isfinite(measurement.covariance.vv) ||
        measurement.covariance.uu < 0.f || measurement.covariance.vv < 0.f) {
      return false;
    }
    hasMeasurement = true;
  }
  if (!hasMeasurement) {
    return false;
  }
  SurfaceKinematicState paramIn{};
  SurfaceKinematicState paramOut{};
  float chi2 = 0.f;
  OperationFailureReason reason{};
  if (!fitTrackSeedLegs(seed, layerGlobals, layerMeasurements, surfaceCatalog, bz,
                        params.ShiftRefToCluster, params.MaxChi2ClusterAttachment, params.MaxChi2NDF,
                        params.RepeatRefitOut, gsl::span<const float>(params.MinPt),
                        paramIn, paramOut, chi2, reason)) {
    return false;
  }
  candidate.seed = seed;
  candidate.track.innerState = paramIn;
  candidate.track.outerState = paramOut;
  candidate.track.chi2 = chi2;
  return fillCandidateKinematics(candidate);
}

} // namespace o2::itsmft::tracking::detail

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_DETECTORREFITSUPPORT_H_
