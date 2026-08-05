// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORTRACKINGOPERATIONADAPTERSUPPORT_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORTRACKINGOPERATIONADAPTERSUPPORT_H_

#ifndef GPUCA_GPUCODE

#include "ITSMFTTracking/CommonTrackShadow.h"
#include "ITSMFTTracking/DetectorTraits.h"
#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"
#include "ITSMFTTracking/TrackingOperationAdapter.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking::detail
{

// These helpers are deliberately kept at the ITS/MFT adapter edge. They
// translate between typed legacy refit/output values and the detector-neutral
// TrackingCandidate consumed by TrackerTraits; the generic core never sees
// DetectorTraits, typed output tracks, or publication sidecars.
template <int NLayers>
bool importCandidateStates(const typename DetectorTraits<NLayers>::TrackType& track,
                           TrackingCandidate& candidate) noexcept
{
  if constexpr (DetectorTraits<NLayers>::DetId == o2::detectors::DetID::ITS) {
    if (!legacy::importBarrelTrackParCov(track.getParamIn(), candidate.innerState) ||
        !legacy::importBarrelTrackParCov(track.getParamOut(), candidate.outerState)) {
      return false;
    }
  } else {
    if (!legacy::importLegacyForwardTrackParCov(track.getTrack(), candidate.innerState) ||
        !legacy::importLegacyForwardTrackParCov(track.getTrack().getOutParam(), candidate.outerState)) {
      return false;
    }
  }
  candidate.chi2 = track.getChi2();
  candidate.phi = track.getPhi();
  candidate.eta = track.getEta();
  if constexpr (DetectorTraits<NLayers>::DetId == o2::detectors::DetID::ITS) {
    candidate.charge = track.getSign();
  } else {
    candidate.charge = track.getCharge();
  }
  return true;
}

template <int NLayers>
bool exportCandidateTrack(const TrackingCandidate& candidate,
                          SurfaceTrackingScratch& scratch,
                          typename DetectorTraits<NLayers>::TrackType& track) noexcept
{
  if constexpr (DetectorTraits<NLayers>::DetId == o2::detectors::DetID::ITS) {
    o2::track::TrackParCovF paramIn{};
    o2::track::TrackParCovF paramOut{};
    if (!legacy::exportBarrelTrackParCov(candidate.innerState, paramIn) ||
        !legacy::exportBarrelTrackParCov(candidate.outerState, paramOut)) {
      return false;
    }
    track.getParamIn() = paramIn;
    track.getParamOut() = paramOut;
    track.setChi2(candidate.chi2);
    track.getTimeStamp() = candidate.timestamp;
    for (int layer = 0; layer < static_cast<int>(candidate.seed.getSurfaceMask().count()); ++layer) {
      const int cluster = candidate.seed.getCluster(layer);
      if (cluster != o2::its::constants::UnusedIndex) {
        track.setExternalClusterIndex(layer, cluster, true);
      }
    }
  } else {
    o2::track::TrackParCovFwd paramIn{};
    o2::track::TrackParCovFwd paramOut{};
    if (!legacy::exportLegacyForwardTrackParCov(candidate.innerState, paramIn) ||
        !legacy::exportLegacyForwardTrackParCov(candidate.outerState, paramOut)) {
      return false;
    }
    static_cast<o2::track::TrackParCovFwd&>(track.getTrack()) = paramIn;
    track.getTrack().setOutParam(paramOut);
    track.getTrack().setTrackChi2(candidate.chi2);
    track.getTrack().setCA(true);
    track.getTrack().setNumberOfPoints(candidate.seed.getSurfaceMask().count());
    track.setPattern(0);
    for (int layer = 0; layer < static_cast<int>(candidate.seed.getSurfaceMask().count()); ++layer) {
      const int cluster = candidate.seed.getCluster(layer);
      track.setClusterIndex(layer, cluster);
      if (cluster != o2::its::constants::UnusedIndex) {
        track.setClusterSize(layer, scratch.getClusterSize(layer, cluster));
      }
    }
    track.setSeedPattern(candidate.seed.getHitLayerMask().value());
    track.getTimeStamp() = candidate.timestamp;
  }
  return true;
}

template <int NLayers>
bool makeCandidateShadow(const TrackingCandidate& candidate,
                         gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                         CommonTrackShadowRecord& record) noexcept
{
  record = {};
  record.track.innerState = candidate.innerState;
  record.track.outerState = candidate.outerState;
  record.track.chi2 = candidate.chi2;
  const auto timestamp = candidate.timestamp.getTimeStamp();
  const auto timestampError = candidate.timestamp.getTimeStampError();
  record.track.timestamp = {static_cast<TFBC>(timestamp - timestampError), static_cast<TFBC>(timestamp + timestampError)};
  record.references.reserve(layerMeasurements.size());
  for (std::size_t layer = 0; layer < layerMeasurements.size(); ++layer) {
    const int localIndex = candidate.getClusterIndex(static_cast<int>(layer));
    if (localIndex == o2::its::constants::UnusedIndex) {
      continue;
    }
    if (localIndex < 0 || static_cast<size_t>(localIndex) >= layerMeasurements[layer].size()) {
      return false;
    }
    const auto& measurement = layerMeasurements[layer][localIndex];
    const TrackClusterReference reference{measurement.surface, SurfaceMeasurementIndex{static_cast<uint32_t>(localIndex)}};
    if (!reference.surface.isValid() || measurement.surface != reference.surface || !measurement.cluster.isValid()) {
      return false;
    }
    record.references.push_back(reference);
    record.track.hitSurfaces.set(reference.surface);
  }
  return !record.references.empty();
}

} // namespace o2::itsmft::tracking::detail

#endif // !GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_DETECTORTRACKINGOPERATIONADAPTERSUPPORT_H_
