// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORTRACKINGOPERATIONADAPTERSUPPORT_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORTRACKINGOPERATIONADAPTERSUPPORT_H_

#ifndef GPUCA_GPUCODE

#include <cmath>

#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "ITSMFTTracking/NativeRefitDriver.h"
#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"
#include "ITSMFTTracking/TrackingOperationAdapter.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking::detail
{

template <o2::detectors::DetID::ID DetId>
void configureAdapterBeamPosition(TimeFrame& frame,
                                  const TrackingParameters& params,
                                  const o2::dataformats::MeanVertexObject* meanVertex,
                                  bool overrideBeamEstimation)
{
  const float systErrY2 = params.SystError2Row.empty() ? 0.f : params.SystError2Row[0];
  const float layerRes = params.LayerResolution.empty() ? 0.f : params.LayerResolution[0];
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    frame.setBeamPosition(params.Diamond[0], params.Diamond[1], params.DiamondCov[3], layerRes, systErrY2);
  } else if (overrideBeamEstimation && meanVertex != nullptr) {
    frame.setBeamPosition(meanVertex->getX(), meanVertex->getY(), meanVertex->getSigmaY2(), layerRes, systErrY2);
  } else if (params.UseDiamond) {
    frame.setBeamPosition(params.Diamond[0], params.Diamond[1], params.DiamondCov[3], layerRes, systErrY2);
  }
}

// These helpers are the only detector-specific operation code used by the
// two application adapters. They produce the same detector-neutral
// TrackingCandidate consumed by TrackerTraits; no typed accepted track or
// detector traits object crosses into the generic core.
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
                         gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                         SurfaceCatalogView surfaceCatalog,
                         TrackingCandidate& candidate)
{
  SurfaceKinematicState paramIn{};
  SurfaceKinematicState paramOut{};
  float chi2 = 0.f;
  OperationFailureReason reason{};
  if (!fitTrackSeedLegs(seed, layerMeasurements, surfaceCatalog, bz,
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
                         const SurfaceTrackingScratch& scratch,
                         gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                         SurfaceCatalogView surfaceCatalog,
                         ClusterSourceId expectedSource,
                         TrackingCandidate& candidate)
{
  SurfaceKinematicState paramIn{};
  SurfaceKinematicState paramOut{};
  float chi2 = 0.f;
  if (!refitTrackFwd(seed, scratch, params, bz, layerMeasurements, surfaceCatalog, expectedSource, paramIn, paramOut, chi2)) {
    return false;
  }
  candidate.seed = seed;
  candidate.track.innerState = paramIn;
  candidate.track.outerState = paramOut;
  candidate.track.chi2 = chi2;
  return fillCandidateKinematics(candidate);
}

template <o2::detectors::DetID::ID DetId>
bool refitDetectorSeed(const TrackSeed& seed,
                       const TrackingParameters& params,
                       float bz,
                       SurfaceTrackingScratch& scratch,
                       gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                       SurfaceCatalogView surfaceCatalog,
                       ClusterSourceId expectedSource,
                       TrackingCandidate& candidate)
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    return refitMFTSeed(seed, params, bz, scratch, layerMeasurements, surfaceCatalog, expectedSource, candidate);
  } else {
    return refitITSSeed(seed, params, bz, layerMeasurements, surfaceCatalog, candidate);
  }
}

} // namespace o2::itsmft::tracking::detail

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_DETECTORTRACKINGOPERATIONADAPTERSUPPORT_H_
