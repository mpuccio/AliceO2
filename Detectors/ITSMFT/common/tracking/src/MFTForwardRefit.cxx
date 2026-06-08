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
/// \file MFTForwardRefit.cxx
/// \brief
///

#include "ITSMFTTracking/MFTForwardRefit.h"

#include <algorithm>
#include <cmath>

#include "Framework/Logger.h"
#include "ITStracking/Constants.h"
#include "MFTTracking/Cluster.h"
#include "MFTTracking/Constants.h"
#include "MFTTracking/MFTTrackingParam.h"
#include "MFTTracking/TrackCA.h"
#include "MFTTracking/TrackFitter.h"

namespace o2::itsmft::tracking
{

namespace
{
template <typename TrackLTFType>
bool buildTrackLTF(const TrackSeed<constants::MFTNLayers>& seed,
                   const TimeFrame<constants::MFTNLayers>& tf,
                   const TrackingParameters& params,
                   TrackLTFType& ltf)
{
  ltf = TrackLTFType(true);
  const auto& unsorted = tf.getUnsortedClusters();
  const auto hitMask = seed.getHitLayerMask();

  for (int layer = 0; layer < constants::MFTNLayers; ++layer) {
    if (!hitMask.has(layer)) {
      continue;
    }
    const int clIdx = seed.getCluster(layer);
    if (clIdx == o2::its::constants::UnusedIndex) {
      continue;
    }
    if (clIdx >= static_cast<int>(unsorted[layer].size())) {
      LOGP(warn, "MFT CA forward refit: invalid cluster index {} on layer {}", clIdx, layer);
      return false;
    }
    const auto& cluster = unsorted[layer][clIdx];
    const auto& tfInfo = tf.getClusterTrackingFrameInfo(layer, cluster);
    o2::mft::Cluster mftCluster{
      tfInfo.xCoordinate,
      tfInfo.yCoordinate,
      tfInfo.zCoordinate,
      cluster.phi,
      cluster.radius,
      clIdx,
      0,
      tfInfo.covarianceTrackingFrame[0],
      tfInfo.covarianceTrackingFrame[2],
      0};
    const int extIdx = tf.getClusterExternalIndex(layer, clIdx);
    const int clsSize = tf.getClusterSize(0, extIdx);
    ltf.setPoint(mftCluster, layer, clIdx, {}, extIdx, clsSize);
  }

  if (ltf.getNumberOfPoints() < params.MinTrackLength) {
    return false;
  }
  ltf.sort();
  return true;
}

void copyClusterRefs(const TrackSeed<constants::MFTNLayers>& seed,
                     const TimeFrame<constants::MFTNLayers>& tf,
                     MFTCATrack& track)
{
  track.setPattern(0);
  const auto hitMask = seed.getHitLayerMask();
  for (int layer = 0; layer < constants::MFTNLayers; ++layer) {
    if (!hitMask.has(layer)) {
      track.setClusterIndex(layer, o2::its::constants::UnusedIndex);
      continue;
    }
    const int clIdx = seed.getCluster(layer);
    track.setClusterIndex(layer, clIdx);
    const int extIdx = tf.getClusterExternalIndex(layer, clIdx);
    track.setClusterSize(layer, tf.getClusterSize(0, extIdx));
  }
}

template <typename TrackLTFType>
bool fitTrackLTF(TrackLTFType& ltf, float bz)
{
  const auto& mftParam = o2::mft::MFTTrackingParam::Instance();
  o2::mft::TrackFitter<TrackLTFType> fitter;
  fitter.setBz(bz);
  fitter.setMFTRadLength(mftParam.MFTRadLength);
  fitter.setTrackModel(mftParam.trackmodel);
  fitter.setAlignResiduals(mftParam.alignResidual);

  TrackLTFType outward = ltf;
  if (!fitter.initTrack(ltf) || !fitter.fit(ltf)) {
    return false;
  }
  if (!fitter.initTrack(outward, true) || !fitter.fit(outward, true)) {
    return false;
  }
  ltf.setOutParam(outward);
  return true;
}

template <typename TrackLTFType>
bool refitTrackFwdImpl(const TrackSeed<constants::MFTNLayers>& seed,
                       MFTCATrack& track,
                       const TimeFrame<constants::MFTNLayers>& tf,
                       const TrackingParameters& params,
                       float bz)
{
  TrackLTFType ltf;
  if (!buildTrackLTF(seed, tf, params, ltf)) {
    return false;
  }

  if (!fitTrackLTF(ltf, bz)) {
    return false;
  }

  const int nCl = ltf.getNumberOfPoints();
  if (nCl < params.MinTrackLength) {
    return false;
  }
  const float minPt = params.MinPt[constants::MFTNLayers - nCl];
  if (ltf.getPt() < minPt) {
    return false;
  }
  const float chi2ndf = static_cast<float>(ltf.getTrackChi2() / std::max(1, 2 * nCl - 5));
  if (chi2ndf > params.MaxChi2NDF) {
    return false;
  }

  if (params.TrackletMinAbsX > 0.f) {
    if (std::abs(ltf.getX()) < params.TrackletMinAbsX) {
      return false;
    }
    const auto hitMask = seed.getHitLayerMask();
    for (int layer = 0; layer < constants::MFTNLayers; ++layer) {
      if (!hitMask.has(layer)) {
        continue;
      }
      const int clIdx = seed.getCluster(layer);
      if (clIdx == o2::its::constants::UnusedIndex) {
        continue;
      }
      if (std::abs(tf.getUnsortedClusters()[layer][clIdx].xCoordinate) < params.TrackletMinAbsX) {
        return false;
      }
    }
  }

  auto& mftTr = track.getTrack();
  mftTr = static_cast<const o2::mft::TrackMFT&>(ltf);
  mftTr.setCA(true);
  copyClusterRefs(seed, tf, track);
  return true;
}
} // namespace

bool refitTrackFwd(const TrackSeed<constants::MFTNLayers>& seed,
                                          MFTCATrack& track,
                                          const TimeFrame<constants::MFTNLayers>& tf,
                                          const TrackingParameters& params,
                                          float bz)
{
  const auto& mftParam = o2::mft::MFTTrackingParam::Instance();
  const bool fieldOff = mftParam.forceZeroField || std::abs(bz) < 1e-6f;
  if (fieldOff) {
    return refitTrackFwdImpl<o2::mft::TrackLTFL>(seed, track, tf, params, 0.f);
  }
  return refitTrackFwdImpl<o2::mft::TrackLTF>(seed, track, tf, params, bz);
}

} // namespace o2::itsmft::tracking
