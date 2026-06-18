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

#include <algorithm>

#include "Framework/Logger.h"
#include "ITSMFTTracking/MFTForwardRefit.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/TrackingParamRef.h"
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
  // Common TrackSeed matches o2::its::TrackSeed layout (hole-layer mask aside); refit uses ITS helper.
  const auto& itsSeed = reinterpret_cast<const o2::its::TrackSeed<NLayers>&>(seed);
  return o2::its::track::refitTrack<NLayers>(itsSeed,
                                             track,
                                             params.MaxChi2ClusterAttachment,
                                             params.MaxChi2NDF,
                                             bz,
                                             tfInfos,
                                             unsortedClusters,
                                             params.LayerxX0.data(),
                                             params.LayerRadii.data(),
                                             params.MinPt.data(),
                                             propagator,
                                             params.CorrType,
                                             params.ReseedIfShorter,
                                             params.ShiftRefToCluster,
                                             params.RepeatRefitOut);
}
} // namespace

template <int NLayers>
bool DetectorTraits<NLayers>::fitCellClusters(const int clusId[3],
                                             const int hitLayers[3],
                                             const TimeFrameN& tf,
                                             const TrackingParameters& params,
                                             float bz,
                                             const o2::base::PropagatorImpl<float>* propagator,
                                             o2::track::TrackParCovF& track,
                                             float& chi2)
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    o2::track::TrackParCovFwd fwdTrack;
    if (!detail::mftFwdFitCellClusters(hitLayers, clusId, tf, params, bz, fwdTrack, chi2)) {
      return false;
    }
    track = detail::mftFwdToBarrelTrack(fwdTrack);
    return true;
  } else {
    const auto& cluster1_glo = tf.getUnsortedClusters()[hitLayers[0]][clusId[0]];
    const auto& cluster2_glo = tf.getUnsortedClusters()[hitLayers[1]][clusId[1]];
    const auto& cluster3_tf = tf.getTrackingFrameInfoOnLayer(hitLayers[2])[clusId[2]];
    track = o2::its::track::buildTrackSeed(cluster1_glo, cluster2_glo, cluster3_tf, bz);
    chi2 = 0.f;
    for (int iC{2}; iC--;) {
      const int hitLayer = hitLayers[iC];
      const TrackingFrameInfo& trackingHit = tf.getTrackingFrameInfoOnLayer(hitLayer)[clusId[iC]];

      if (!track.rotate(trackingHit.alphaTrackingFrame)) {
        return false;
      }

      if (!track.propagateTo(trackingHit.xTrackingFrame, bz)) {
        return false;
      }

      if (!track.correctForMaterial(params.LayerxX0[hitLayer], params.LayerxX0[hitLayer] * o2::its::constants::Radl * o2::its::constants::Rho, true)) {
        return false;
      }

      const auto predChi2{track.getPredictedChi2Quiet(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)};
      if (!iC && predChi2 > params.MaxChi2ClusterAttachment) {
        return false;
      }

      if (!track.o2::track::TrackParCov::update(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)) {
        return false;
      }

      if (!iC) {
        chi2 += predChi2;
      }
    }
    return true;
  }
}

template <int NLayers>
bool DetectorTraits<NLayers>::attachNeighbourToSeed(TrackSeedN& seed,
                                                     int neighbourLayer,
                                                     int neighbourCluster,
                                                     const TimeFrameN& tf,
                                                     const TrackingParameters& params,
                                                     float bz,
                                                     const o2::base::PropagatorImpl<float>* propagator)
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    o2::track::TrackParCovFwd fwdTrack;
    float chi2 = seed.getChi2();
    if (!detail::mftFwdAttachNeighbourToSeed(seed, neighbourLayer, neighbourCluster, tf, params, bz, fwdTrack, chi2)) {
      return false;
    }
    static_cast<o2::track::TrackParCovFwd&>(seed) = fwdTrack;
    seed.setChi2(chi2);
    return true;
  } else {
    const auto& trHit = tf.getTrackingFrameInfoOnLayer(neighbourLayer)[neighbourCluster];

    if (!seed.rotate(trHit.alphaTrackingFrame)) {
      return false;
    }

    if (!propagator->propagateToX(seed, trHit.xTrackingFrame, bz, o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, params.CorrType)) {
      return false;
    }

    if (params.CorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE) {
      if (!seed.correctForMaterial(params.LayerxX0[neighbourLayer], params.LayerxX0[neighbourLayer] * o2::its::constants::Radl * o2::its::constants::Rho, true)) {
        return false;
      }
    }

    auto predChi2{seed.getPredictedChi2Quiet(trHit.positionTrackingFrame, trHit.covarianceTrackingFrame)};
    if ((predChi2 > params.MaxChi2ClusterAttachment) || predChi2 < 0.f) {
      return false;
    }
    seed.setChi2(seed.getChi2() + predChi2);
    return seed.o2::track::TrackParCov::update(trHit.positionTrackingFrame, trHit.covarianceTrackingFrame);
  }
}

template <int NLayers>
bool DetectorTraits<NLayers>::cellsAreCompatible(const CellSeedN& currentCell,
                                                  const CellSeedN& nextCell,
                                                  const TimeFrameN& tf,
                                                  const TrackingParameters& params,
                                                  float bz)
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    return detail::mftFwdCellsAreCompatible(currentCell, nextCell, bz, params.MaxChi2ClusterAttachment);
  } else {
    auto nextCellSeed = nextCell;
    if (!nextCellSeed.rotate(currentCell.getAlpha()) ||
        !nextCellSeed.propagateTo(currentCell.getX(), bz)) {
      return false;
    }
    return currentCell.getPredictedChi2(nextCellSeed) <= params.MaxChi2ClusterAttachment;
  }
}

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
void DetectorTraits<NLayers>::sortRefittedTracks(bounded_vector<TrackType>& tracks)
{
  // Same ordering as o2::its::track::isBetter (longer track, then lower chi2).
  std::sort(tracks.begin(), tracks.end(), [](const TrackType& a, const TrackType& b) {
    const auto ncla = a.getNumberOfClusters();
    const auto nclb = b.getNumberOfClusters();
    return (ncla == nclb) ? (a.getChi2() < b.getChi2()) : ncla > nclb;
  });
}

template <int NLayers>
void DetectorTraits<NLayers>::finalizeAcceptedTrack(TrackType& track)
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    track.setUserField(0);
    track.getParamOut().setUserField(0);
  }
}

template <int NLayers>
bool DetectorTraits<NLayers>::sameTrackSign(const TrackType& t1, const TrackType& t2)
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    return t1.getCharge() == t2.getCharge();
  } else {
    return t1.getSign() == t2.getSign();
  }
}

template <int NLayers>
bool DetectorTraits<NLayers>::validateMFTCellClusters(const o2::its::Cluster& c0, int layer0,
                                                      const o2::its::Cluster& c1, int layer1,
                                                      const o2::its::Cluster& c2, int layer2,
                                                      float r2Cut)
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    return detail::mftDistanceToSeedSquared(c0, c2, c1) < r2Cut * detail::mftConicalRoadR2Scale(layer0, layer1) &&
           detail::mftDistanceToSeedSquared(c0, c1, c2) < r2Cut * detail::mftConicalRoadR2Scale(layer0, layer2) &&
           detail::mftDistanceToSeedSquared(c1, c2, c0) < r2Cut * detail::mftConicalRoadR2Scale(layer1, layer0);
  } else {
    return false;
  }
}

template <int NLayers>
bool DetectorTraits<NLayers>::mftCellsConnect(const o2::its::Cluster& cEnd, const o2::its::Cluster& cStart, float r2Cut)
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    const float dx = cEnd.xCoordinate - cStart.xCoordinate;
    const float dy = cEnd.yCoordinate - cStart.yCoordinate;
    return dx * dx + dy * dy <= r2Cut;
  } else {
    return false;
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
