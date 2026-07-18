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
#include "ITStracking/TrackHelpers.h"

namespace o2::itsmft::tracking
{

namespace
{
template <int NLayers>
bool refitSeedITS(const typename DetectorTraits<NLayers>::TrackSeedN& seed,
                  o2::its::TrackITSExt& track,
                  const TrackingParameters& params,
                  float bz,
                  const o2::its::TrackingFrameInfo* const tfInfos[NLayers],
                  const o2::its::Cluster* const unsortedClusters[NLayers],
                  const o2::base::PropagatorImpl<float>* propagator)
{
  o2::its::TrackSeed<NLayers> itsSeed;
  static_cast<o2::track::TrackParCovF&>(itsSeed) = static_cast<const o2::track::TrackParCovF&>(seed);
  itsSeed.setHitLayerMask(o2::its::LayerMask{seed.getHitLayerMask().value()});
  itsSeed.setFirstTrackletIndex(seed.getFirstTrackletIndex());
  itsSeed.setSecondTrackletIndex(seed.getSecondTrackletIndex());
  itsSeed.setChi2(seed.getChi2());
  itsSeed.setLevel(seed.getLevel());
  itsSeed.getTimeStamp() = seed.getTimeStamp();
  for (int iLayer{0}; iLayer < NLayers; ++iLayer) {
    itsSeed.getClusters()[iLayer] = seed.getCluster(iLayer);
  }
  const o2::its::track::TrackFitContext<NLayers> fitCtx{
    tfInfos, params.LayerxX0.data(), params.NLayers, bz,
    params.MaxChi2ClusterAttachment, params.MaxChi2NDF,
    propagator, params.CorrType, params.ShiftRefToCluster, params.RepeatRefitOut};
  o2::its::TrackITSInternal<NLayers> internalTrack;
  if (!o2::its::track::refitTrackSeed<NLayers>(itsSeed,
                                               internalTrack,
                                               fitCtx,
                                               unsortedClusters,
                                               params.LayerRadii.data(),
                                               params.MinPt.data(),
                                               params.ReseedIfShorter)) {
    return false;
  }
  track = o2::its::makeTrackITSExt(internalTrack);
  return true;
}
} // namespace

template <int NLayers>
bool DetectorTraits<NLayers>::refitSeed(const TrackSeedN& seed,
                                         TrackType& track,
                                         const TrackingParameters& params,
                                         float bz,
                                         TimeFrameN& tf,
                                         const o2::its::TrackingFrameInfo* const tfInfos[NLayers],
                                         const o2::its::Cluster* const unsortedClusters[NLayers],
                                         const o2::base::PropagatorImpl<float>* propagator)
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    return refitTrackFwd(seed, track, tf, params, bz);
  } else {
    return refitSeedITS<NLayers>(seed, track, params, bz, tfInfos, unsortedClusters, propagator);
  }
}

template <int NLayers>
void DetectorTraits<NLayers>::configureIndexTableUtils(IndexTableUtils<NLayers>& utils, const TrackingParameters& params)
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    constexpr float defaultYMin{-20.f};
    constexpr float defaultYMax{20.f};
    const bool hasRowRange = params.IndexRowMax != 0.f;
    const float rowMin = hasRowRange ? params.IndexRowMin : defaultYMin;
    const float rowMax = hasRowRange ? params.IndexRowMax : defaultYMax;
    utils.setTrackingParametersXY(params, rowMin, rowMax);
  } else {
    utils.setTrackingParameters(params);
  }
}

template <o2::detectors::DetID::ID DetId, int NLayers>
void TrackingLoadPolicy<DetId, NLayers>::configureBeamPosition(TimeFrame<NLayers>& tf,
                                                                 const TrackingParameters& p,
                                                                 const o2::dataformats::MeanVertexObject* meanVertex,
                                                                 bool overrideBeamEstimation)
{
  const auto& tc = TrackerParamRef<DetId>::get();
  const float systErrY2 = p.SystError2Row.empty() ? 0.f : p.SystError2Row[0];
  const float layerRes = p.LayerResolution.empty() ? 0.f : p.LayerResolution[0];

  if constexpr (DetId == o2::detectors::DetID::MFT) {
    tf.setBeamPosition(p.Diamond[0], p.Diamond[1], p.DiamondCov[3], layerRes, systErrY2);
    LOGP(info, "MFT CA vertex seed from diamond: x={:.4f} y={:.4f} z={:.4f}",
         p.Diamond[0], p.Diamond[1], p.Diamond[2]);
  } else if ((overrideBeamEstimation || tc.overrideBeamEstimation) && meanVertex != nullptr) {
    tf.setBeamPosition(meanVertex->getX(), meanVertex->getY(), meanVertex->getSigmaY2(), layerRes, systErrY2);
    LOGP(info, "ITS CA beam position from MeanVertex: x={:.4f} y={:.4f}", meanVertex->getX(), meanVertex->getY());
  } else if (p.UseDiamond) {
    tf.setBeamPosition(p.Diamond[0], p.Diamond[1], p.DiamondCov[3], layerRes, systErrY2);
    LOGP(info, "ITS CA beam position from diamond: x={:.4f} y={:.4f} z={:.4f}",
         p.Diamond[0], p.Diamond[1], p.Diamond[2]);
  }
}

template struct DetectorTraits<7>;
template struct DetectorTraits<10>;
template struct TrackingLoadPolicy<o2::detectors::DetID::ITS, 7>;
template struct TrackingLoadPolicy<o2::detectors::DetID::MFT, 10>;

} // namespace o2::itsmft::tracking
