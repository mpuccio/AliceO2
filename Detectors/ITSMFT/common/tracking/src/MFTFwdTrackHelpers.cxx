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
/// \file MFTFwdTrackHelpers.cxx
/// \brief MFT CA final forward Kalman refit (legacy TrackFitter)
///

#include "ITSMFTTracking/MFTFwdTrackHelpers.h"

#include <cmath>

#include "Framework/Logger.h"
#include "ITStracking/Constants.h"
#include "MFTTracking/Cluster.h"
#include "MFTTracking/Constants.h"
#include "MFTTracking/MFTTrackingParam.h"
#include "MFTTracking/TrackFitter.h"

namespace o2::itsmft::tracking
{

namespace
{
template <typename TrackLTFType>
bool refitTrackFwdImpl(const TrackSeedN<o2::mft::constants::mft::LayersNumber>& seed,
                       MFTCATrack& track,
                       const TimeFrame<o2::mft::constants::mft::LayersNumber>& tf,
                       const TrackingParameters& params,
                       float bz)
{
  TrackLTFType ltf(true);
  const auto& unsorted = tf.getUnsortedClusters();
  const auto hitMask = seed.getHitLayerMask();

  for (int layer = 0; layer < o2::mft::constants::mft::LayersNumber; ++layer) {
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
      tfInfo.xCoordinate, tfInfo.yCoordinate, tfInfo.zCoordinate,
      cluster.phi, cluster.radius, clIdx, 0,
      tfInfo.covarianceTrackingFrame[0], tfInfo.covarianceTrackingFrame[2], 0};
    const int extIdx = tf.getClusterExternalIndex(layer, clIdx);
    ltf.setPoint(mftCluster, layer, clIdx, {}, extIdx, tf.getClusterSize(0, extIdx));
  }

  if (ltf.getNumberOfPoints() < params.MinTrackLength) {
    return false;
  }
  ltf.sort();

  const auto& mftParam = o2::mft::MFTTrackingParam::Instance();
  o2::mft::TrackFitter<TrackLTFType> fitter;
  fitter.setBz(bz);
  fitter.setMFTRadLength(mftParam.MFTRadLength);
  fitter.setTrackModel(mftParam.trackmodel);
  fitter.setAlignResiduals(mftParam.alignResidual);

  TrackLTFType outward = ltf;
  if (!fitter.initTrack(ltf) || !fitter.fit(ltf) ||
      !fitter.initTrack(outward, true) || !fitter.fit(outward, true)) {
    return false;
  }
  ltf.setOutParam(outward);

  const int nCl = ltf.getNumberOfPoints();
  if (nCl < params.MinTrackLength) {
    return false;
  }
  if (ltf.getPt() < params.MinPt[o2::mft::constants::mft::LayersNumber - nCl]) {
    return false;
  }
  if (static_cast<float>(ltf.getTrackChi2() / std::max(1, 2 * nCl - 5)) > params.MaxChi2NDF) {
    return false;
  }

  if (params.TrackletMinAbsX > 0.f) {
    if (std::abs(ltf.getX()) < params.TrackletMinAbsX) {
      return false;
    }
    for (int layer = 0; layer < o2::mft::constants::mft::LayersNumber; ++layer) {
      if (!hitMask.has(layer)) {
        continue;
      }
      const int clIdx = seed.getCluster(layer);
      if (clIdx != o2::its::constants::UnusedIndex &&
          std::abs(unsorted[layer][clIdx].xCoordinate) < params.TrackletMinAbsX) {
        return false;
      }
    }
  }

  auto& mftTr = track.getTrack();
  mftTr = static_cast<const o2::mft::TrackMFT&>(ltf);
  mftTr.setCA(true);

  track.setPattern(0);
  for (int layer = 0; layer < o2::mft::constants::mft::LayersNumber; ++layer) {
    if (!hitMask.has(layer)) {
      track.setClusterIndex(layer, o2::its::constants::UnusedIndex);
      continue;
    }
    const int clIdx = seed.getCluster(layer);
    track.setClusterIndex(layer, clIdx);
    const int extIdx = tf.getClusterExternalIndex(layer, clIdx);
    track.setClusterSize(layer, tf.getClusterSize(0, extIdx));
  }
  return true;
}
} // namespace

bool refitTrackFwd(const TrackSeedN<o2::mft::constants::mft::LayersNumber>& seed,
                   MFTCATrack& track,
                   const TimeFrame<o2::mft::constants::mft::LayersNumber>& tf,
                   const TrackingParameters& params,
                   float bz)
{
  const auto& mftParam = o2::mft::MFTTrackingParam::Instance();
  if (mftParam.forceZeroField || std::abs(bz) < 1e-6f) {
    return refitTrackFwdImpl<o2::mft::TrackLTFL>(seed, track, tf, params, 0.f);
  }
  return refitTrackFwdImpl<o2::mft::TrackLTF>(seed, track, tf, params, bz);
}

} // namespace o2::itsmft::tracking
