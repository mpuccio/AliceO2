// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORREFITSUPPORT_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORREFITSUPPORT_H_

#ifndef GPUCA_GPUCODE

#include <cmath>

#include "DataFormatsCalibration/MeanVertexObject.h"
#include "ITSMFTTracking/CommonTrack.h"
#include "ITSMFTTracking/NativeRefitDriver.h"
#include "ITSMFTTracking/detail/SurfaceKinematicStateLegacyAdapters.h"
#include "ITSMFTTracking/detail/MFTFwdTrackHelpers.h"
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

// Detector edges produce the detector-neutral TrackingCandidate consumed by
// TrackerTraits; typed accepted tracks stay outside the generic core.
inline bool fillCandidateKinematics(TrackingCandidate& candidate) noexcept
{
  if (candidate.track.innerState.family == StateFamily::Barrel) {
    o2::track::TrackParCovF param{};
    if (!legacy::exportBarrelTrackParCov(candidate.track.innerState, param)) {
      return false;
    }
    candidate.phi = param.getPhi();
    candidate.eta = param.getEta();
    candidate.charge = param.getSign();
    return true;
  }
  if (candidate.track.innerState.family == StateFamily::Forward) {
    o2::track::TrackParCovFwd param{};
    if (!legacy::exportLegacyForwardTrackParCov(candidate.track.innerState, param)) {
      return false;
    }
    candidate.phi = static_cast<float>(param.getPhi());
    candidate.eta = static_cast<float>(param.getEta());
    candidate.charge = param.getCharge();
    return true;
  }
  return false;
}

inline bool refitITSSeed(const TrackSeed& seed,
                         const TrackingParameters& params,
                         float bz,
                         SurfaceTrackingScratch&,
                         gsl::span<const gsl::span<const GlobalMeasurement>> layerGlobals,
                         gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                         SurfaceCatalogView surfaceCatalog,
                         ClusterSourceId,
                         TrackingCandidate& candidate)
{
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

inline bool refitMFTSeed(const TrackSeed& seed,
                         const TrackingParameters& params,
                         float bz,
                         SurfaceTrackingScratch& scratch,
                         gsl::span<const gsl::span<const GlobalMeasurement>> layerGlobals,
                         gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                         SurfaceCatalogView surfaceCatalog,
                         ClusterSourceId expectedSource,
                         TrackingCandidate& candidate)
{
  SurfaceKinematicState paramIn{};
  SurfaceKinematicState paramOut{};
  float chi2 = 0.f;
  if (!refitTrackFwd(seed, scratch, params, bz, layerGlobals, layerMeasurements, surfaceCatalog, expectedSource, paramIn, paramOut, chi2)) {
    return false;
  }
  candidate.seed = seed;
  candidate.track.innerState = paramIn;
  candidate.track.outerState = paramOut;
  candidate.track.chi2 = chi2;
  return fillCandidateKinematics(candidate);
}

inline bool refitSurfaceSeed(const TrackSeed& seed,
                             const TrackingParameters& params,
                             float bz,
                             SurfaceTrackingScratch& scratch,
                             gsl::span<const gsl::span<const GlobalMeasurement>> layerGlobals,
                             gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                             SurfaceCatalogView surfaceCatalog,
                             ClusterSourceId expectedSource,
                             TrackingCandidate& candidate)
{
  for (int position = 0; position < static_cast<int>(layerMeasurements.size()); ++position) {
    const int clusterIndex = seed.getCluster(position);
    if (clusterIndex == o2::its::constants::UnusedIndex) {
      continue;
    }
    if (clusterIndex < 0 || clusterIndex >= static_cast<int>(layerMeasurements[position].size())) {
      return false;
    }
    const auto surface = layerGlobals[position][clusterIndex].surface;
    if (!surfaceCatalog.hasSurface(surface)) {
      return false;
    }
    return surfaceCatalog.getSurface(surface).kind == SurfaceKind::Cylinder
             ? refitITSSeed(seed, params, bz, scratch, layerGlobals, layerMeasurements, surfaceCatalog, expectedSource, candidate)
             : refitMFTSeed(seed, params, bz, scratch, layerGlobals, layerMeasurements, surfaceCatalog, expectedSource, candidate);
  }
  return false;
}

} // namespace o2::itsmft::tracking::detail

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_DETECTORREFITSUPPORT_H_
