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
///
/// \file DetectorTraits.cxx
/// \brief Detector-specific refit and load hooks for shared CA tracking
///

#include "ITSMFTTracking/DetectorTraits.h"

#include "Framework/Logger.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/NativeRefitDriver.h"
#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"

namespace o2::itsmft::tracking
{

namespace
{
// M5d: the barrel/ITS branch of DetectorTraits::refitSeed. Replaces the
// frozen o2::its::track::fitTrack/refitTrack/refitTrackSeed chain
// (ITStracking/TrackHelpers.h) with the shared, descriptor-driven native
// driver (fitTrackSeedLegs, NativeRefitDriver.h) that also serves the
// forward/MFT branch below -- this is the intentional, approved physics
// departure recorded in doc/decisions/0008-native-refit-activation.md.
// legacy::exportBarrelTrackParCov is retained: it is a plain data-format
// adapter into TrackITSExt's own TrackParCovF-based storage, not a fitting
// algorithm, and every detector output type in this library is still
// exported this way.
template <int NLayers>
bool refitSeedITS(const typename DetectorTraits<NLayers>::TrackSeedN& seed,
                  o2::its::TrackITSExt& track,
                  const TrackingParameters& params,
                  float bz,
                  const LayerMeasurementSpans<NLayers>& layerMeasurements,
                  SurfaceCatalogView surfaceCatalog)
{
  SurfaceKinematicState paramIn{};
  SurfaceKinematicState paramOut{};
  float chi2 = 0.f;
  OperationFailureReason reason{};
  if (!fitTrackSeedLegs<NLayers>(seed, layerMeasurements, surfaceCatalog, bz,
                                 params.ShiftRefToCluster, params.MaxChi2ClusterAttachment, params.MaxChi2NDF,
                                 params.RepeatRefitOut, gsl::span<const float>(params.MinPt),
                                 paramIn, paramOut, chi2, reason)) {
    return false;
  }

  o2::track::TrackParCovF legacyParamIn{};
  o2::track::TrackParCovF legacyParamOut{};
  if (!legacy::exportBarrelTrackParCov(paramIn, legacyParamIn) ||
      !legacy::exportBarrelTrackParCov(paramOut, legacyParamOut)) {
    return false;
  }
  o2::its::TrackITSExt scratch{};
  scratch.getParamIn() = legacyParamIn;
  scratch.getParamOut() = legacyParamOut;
  scratch.setChi2(chi2);
  scratch.getTimeStamp() = seed.getTimeStamp().makeSymmetrical();
  for (int layer = 0; layer < NLayers; ++layer) {
    const int clsIdx = seed.getCluster(layer);
    if (clsIdx != o2::its::constants::UnusedIndex) {
      scratch.setExternalClusterIndex(layer, clsIdx, true);
    }
  }
  track = scratch;
  return true;
}
} // namespace

template <int NLayers, typename ScratchT>
bool DetectorTraits<NLayers, ScratchT>::refitSeed(const TrackSeedN& seed,
                                                  TrackType& track,
                                                  const TrackingParameters& params,
                                                  float bz,
                                                  ScratchN& scratch,
                                                  const LayerMeasurementSpans<NLayers>& layerMeasurements,
                                                  SurfaceCatalogView surfaceCatalog,
                                                  ClusterSourceId expectedSource)
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    return refitTrackFwd(seed, track, scratch, params, bz, layerMeasurements, surfaceCatalog, expectedSource);
  } else {
    return refitSeedITS<NLayers>(seed, track, params, bz, layerMeasurements, surfaceCatalog);
  }
}

template <int NLayers, typename ScratchT>
void DetectorTraits<NLayers, ScratchT>::copySeedPatternToTrack(TrackType& track, const TrackSeedN& seed) noexcept
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    track.setSeedPattern(seed.getHitLayerMask().value());
  }
}

template <int NLayers, typename ScratchT>
void DetectorTraits<NLayers, ScratchT>::clearTransientLayerPattern(TrackType& track) noexcept
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    track.clearExtendedLayerPattern();
  }
}

template <int NLayers, typename ScratchT>
bool DetectorTraits<NLayers, ScratchT>::haveSamePolarity(const TrackType& a, const TrackType& b) noexcept
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    return a.getCharge() == b.getCharge();
  } else {
    return a.getSign() == b.getSign();
  }
}

template <o2::detectors::DetID::ID DetId, int NLayers>
void TrackingLoadPolicy<DetId, NLayers>::configureBeamPosition(TimeFrame& frame,
                                                               const TrackingParameters& p,
                                                               const o2::dataformats::MeanVertexObject* meanVertex,
                                                               bool overrideBeamEstimation)
{
  const float systErrY2 = p.SystError2Row.empty() ? 0.f : p.SystError2Row[0];
  const float layerRes = p.LayerResolution.empty() ? 0.f : p.LayerResolution[0];

  if constexpr (DetId == o2::detectors::DetID::MFT) {
    frame.setBeamPosition(p.Diamond[0], p.Diamond[1], p.DiamondCov[3], layerRes, systErrY2);
    LOGP(info, "MFT CA vertex seed from diamond: x={:.4f} y={:.4f} z={:.4f}",
         p.Diamond[0], p.Diamond[1], p.Diamond[2]);
  } else if (overrideBeamEstimation && meanVertex != nullptr) {
    // ITS common-CA must not consult TrackerParamRef<ITS>::get() (the frozen
    // legacy ITSCATrackerParam's overrideBeamEstimation); the common workflow
    // constructor/configuration (overrideBeamEstimation, p.UseDiamond -- both
    // sourced from ITSCommonCATrackerParam via getTrackingParameters()) is
    // the sole owner of the selected static-diamond constraint.
    frame.setBeamPosition(meanVertex->getX(), meanVertex->getY(), meanVertex->getSigmaY2(), layerRes, systErrY2);
    LOGP(info, "ITS CA beam position from MeanVertex: x={:.4f} y={:.4f}", meanVertex->getX(), meanVertex->getY());
  } else if (p.UseDiamond) {
    frame.setBeamPosition(p.Diamond[0], p.Diamond[1], p.DiamondCov[3], layerRes, systErrY2);
    LOGP(info, "ITS CA beam position from diamond: x={:.4f} y={:.4f} z={:.4f}",
         p.Diamond[0], p.Diamond[1], p.Diamond[2]);
  }
}

// M6d/M6e2: DetectorTraits<7>/<10> (default ScratchT) are still used by any
// remaining bare-default LegacyTrackerScratch<NLayers>-backed callers (e.g.
// testTrackingInterfaceLoadFailureContract); DetectorTraits<10,
// SurfaceTrackingScratch> (M6d) and DetectorTraits<7, SurfaceTrackingScratch>
// (M6e2) are what TrackerTraits<NLayers, SurfaceTrackingScratch,
// SurfacePlanBinding>::acceptTracks() (via findRoadsForPolicy()) needs for
// MFT and ITS respectively, now that both live common-CA paths use
// SurfaceTrackingScratch.
template struct DetectorTraits<7>;
template struct DetectorTraits<10>;
template struct DetectorTraits<7, SurfaceTrackingScratch>;
template struct DetectorTraits<10, SurfaceTrackingScratch>;
template struct TrackingLoadPolicy<o2::detectors::DetID::ITS, 7>;
template struct TrackingLoadPolicy<o2::detectors::DetID::MFT, 10>;

} // namespace o2::itsmft::tracking
