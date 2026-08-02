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
                       const LegacyTrackerScratch<o2::mft::constants::mft::LayersNumber>& tf,
                       const TrackingParameters& params,
                       float bz,
                       const LayerMeasurementSpans<o2::mft::constants::mft::LayersNumber>& layerMeasurements,
                       ClusterSourceId expectedSource)
{
  TrackLTFType ltf(true);
  const auto hitMask = seed.getHitLayerMask();

  for (int layer = 0; layer < o2::mft::constants::mft::LayersNumber; ++layer) {
    if (!hitMask.has(layer)) {
      continue;
    }
    const int clIdx = seed.getCluster(layer);
    if (clIdx == o2::its::constants::UnusedIndex) {
      continue;
    }
    if (clIdx < 0 || clIdx >= static_cast<int>(layerMeasurements[layer].size())) {
      LOGP(warn, "MFT CA forward refit: invalid cluster index {} on layer {}", clIdx, layer);
      return false;
    }
    const auto& measurement = layerMeasurements[layer][clIdx];
    const int extIdx = tf.getClusterExternalIndex(layer, clIdx);
    // Defensive re-check of the ClusterRef identity contract that
    // TrackerTraits::initialiseTimeFrame() already established for every
    // entry of mLayerMeasurements (NormalizedMeasurementMismatch). Checked
    // again here because a failure at this final-refit boundary must fail
    // only this one seed (return false) rather than the
    // TraversalException/dropped-TF path that guards the bulk validation.
    if (!measurement.cluster.isValid() || measurement.cluster.source != expectedSource ||
        extIdx < 0 || static_cast<uint32_t>(extIdx) != measurement.cluster.index) {
      LOGP(warn, "MFT CA forward refit: normalized measurement identity mismatch on layer {} clIdx {}", layer, clIdx);
      return false;
    }
    if (!std::isfinite(measurement.global.x) || !std::isfinite(measurement.global.y) || !std::isfinite(measurement.global.z) ||
        !std::isfinite(measurement.covariance.uu) || !std::isfinite(measurement.covariance.vv) ||
        measurement.covariance.uu < 0.f || measurement.covariance.vv < 0.f) {
      LOGP(warn, "MFT CA forward refit: invalid normalized measurement on layer {} clIdx {}", layer, clIdx);
      return false;
    }
    // o2::mft::Cluster's phi/radius constructor arguments are not read by
    // TrackLTF::setPoint (MFTTracking/TrackCA.h), which only consumes
    // x/y/z and sigmaX2/sigmaY2 (via BaseCluster); pass deterministic
    // neutral values rather than reading legacy Cluster::phi/radius or
    // recomputing them from the normalized coordinates. This fitter is
    // diagonal-only: measurement.covariance.uv is intentionally unused,
    // and no off-diagonal approximation is substituted for it.
    o2::mft::Cluster mftCluster{
      measurement.global.x, measurement.global.y, measurement.global.z,
      0.f, 0.f, clIdx, 0,
      measurement.covariance.uu, measurement.covariance.vv, 0};
    ltf.setPoint(mftCluster, layer, clIdx, {}, extIdx, tf.getClusterSize(layer, clIdx));
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
          std::abs(layerMeasurements[layer][clIdx].global.x) < params.TrackletMinAbsX) {
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
    track.setClusterSize(layer, tf.getClusterSize(layer, clIdx));
  }
  return true;
}
} // namespace

bool refitTrackFwd(const TrackSeedN<o2::mft::constants::mft::LayersNumber>& seed,
                   MFTCATrack& track,
                   const LegacyTrackerScratch<o2::mft::constants::mft::LayersNumber>& tf,
                   const TrackingParameters& params,
                   float bz,
                   const LayerMeasurementSpans<o2::mft::constants::mft::LayersNumber>& layerMeasurements,
                   ClusterSourceId expectedSource)
{
  const auto& mftParam = o2::mft::MFTTrackingParam::Instance();
  if (mftParam.forceZeroField || std::abs(bz) < 1e-6f) {
    return refitTrackFwdImpl<o2::mft::TrackLTFL>(seed, track, tf, params, 0.f, layerMeasurements, expectedSource);
  }
  return refitTrackFwdImpl<o2::mft::TrackLTF>(seed, track, tf, params, bz, layerMeasurements, expectedSource);
}

} // namespace o2::itsmft::tracking
