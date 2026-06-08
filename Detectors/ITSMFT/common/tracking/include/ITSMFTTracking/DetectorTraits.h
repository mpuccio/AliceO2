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
/// \file DetectorTraits.h
/// \brief Minimal detector-specific hooks for shared CA tracking (ITS vs MFT)
///

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTOR_TRAITS_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTOR_TRAITS_H_

#include "DetectorsBase/Propagator.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "ITSMFTTracking/CATrackTypes.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/Constants.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Tracklet.h"
#include "MFTTracking/Constants.h"

#include <cmath>
#include <limits>

namespace o2::itsmft::tracking
{

namespace detail
{
inline float mftDistanceToSeedSquared(const o2::its::Cluster& c1, const o2::its::Cluster& c2, const o2::its::Cluster& c)
{
  const float dxSeed = c2.xCoordinate - c1.xCoordinate;
  const float dySeed = c2.yCoordinate - c1.yCoordinate;
  const float dzSeed = c2.zCoordinate - c1.zCoordinate;
  if (std::abs(dzSeed) < 1e-9f) {
    return std::numeric_limits<float>::max();
  }
  const float invdzSeed = (c.zCoordinate - c1.zCoordinate) / dzSeed;
  const float xSeed = c1.xCoordinate + dxSeed * invdzSeed;
  const float ySeed = c1.yCoordinate + dySeed * invdzSeed;
  const float dx = c.xCoordinate - xSeed;
  const float dy = c.yCoordinate - ySeed;
  return dx * dx + dy * dy;
}

inline float mftLayerZ(int layer)
{
  return o2::mft::constants::mft::LayerZCoordinate()[layer];
}

inline bool mftPassesMinAbsX(const o2::its::Cluster& cluster, float minAbsX)
{
  return minAbsX <= 0.f || std::abs(cluster.xCoordinate) >= minAbsX;
}

/// Expected cluster z on toLayer for a straight line from the vertex through cluster on fromLayer.
inline float mftExpectedZAtLayer(float clusterZ, int fromLayer, int toLayer, float pvZ)
{
  const float zFrom = mftLayerZ(fromLayer);
  const float dzLayers = mftLayerZ(toLayer) - zFrom;
  const float denom = zFrom - pvZ;
  if (std::abs(denom) < 1e-6f) {
    return clusterZ + dzLayers;
  }
  return clusterZ + dzLayers * ((clusterZ - pvZ) / denom);
}

/// ITS-style z-resolution estimate with layer-to-layer dz instead of delta-R (MFT forward disks).
inline float mftTrackletSigmaZ(float clusterZ,
                               float clusterRadius,
                               float pvZ,
                               float resolution,
                               float meanDeltaZ,
                               float msAngle)
{
  const float inverseR0 = 1.f / clusterRadius;
  const float tanLambda = (clusterZ - pvZ) * inverseR0;
  const float dz0 = clusterZ - pvZ;
  const float sqInvDeltaZ0 = 1.f / (dz0 * dz0 + o2::its::constants::Tolerance);
  return std::sqrt((resolution * resolution * tanLambda * tanLambda *
                    ((inverseR0 * inverseR0 + sqInvDeltaZ0) * meanDeltaZ * meanDeltaZ + 1.f)) +
                   (meanDeltaZ * msAngle) * (meanDeltaZ * msAngle));
}

/// Cone-projected transverse distance squared (replaces phi-difference cut for MFT).
template <int NLayers>
inline float mftTrackletTransverseDist2(const o2::its::Cluster& current,
                                        const o2::its::Cluster& next,
                                        int fromLayer,
                                        int toLayer)
{
  float xProj = 0.f;
  float yProj = 0.f;
  mftConeProject<NLayers>(current, fromLayer, toLayer, xProj, yProj);
  const float dx = next.xCoordinate - xProj;
  const float dy = next.yCoordinate - yProj;
  return dx * dx + dy * dy;
}
} // namespace detail

/// Per-detector differences in refit, track acceptance, and index-table setup.
/// Everything else stays in TrackerTraits and matches ITS line-for-line.
template <int NLayers>
struct DetectorTraits {
  using TrackType = CATrackType<NLayers>;
  using TrackSeedN = TrackSeed<NLayers>;
  using TimeFrameN = TimeFrame<NLayers>;
  static constexpr o2::detectors::DetID::ID DetId = constants::detIdFromNLayers<NLayers>();

  static bool refitSeed(const TrackSeedN& seed,
                        TrackType& track,
                        const TrackingParameters& params,
                        float bz,
                        TimeFrameN& tf,
                        const o2::its::TrackingFrameInfo* const tfInfos[NLayers],
                        const o2::its::Cluster* const unsortedClusters[NLayers],
                        const o2::base::PropagatorImpl<float>* propagator);

  static void sortRefittedTracks(bounded_vector<TrackType>& tracks);
  static void finalizeAcceptedTrack(TrackType& track);
  static bool sameTrackSign(const TrackType& t1, const TrackType& t2);

  static void configureIndexTableUtils(IndexTableUtils<NLayers>& utils, const TrackingParameters& params);

  static bool validateMFTCellClusters(const o2::its::Cluster& c0, const o2::its::Cluster& c1, const o2::its::Cluster& c2, float r2Cut);
  static bool mftCellsConnect(const o2::its::Cluster& cEnd, const o2::its::Cluster& cStart, float r2Cut);
};

template <o2::detectors::DetID::ID DetId, int NLayers>
struct TrackingLoadPolicy {
  static void configureBeamPosition(TimeFrame<NLayers>& tf,
                                    const TrackingParameters& p,
                                    const o2::dataformats::MeanVertexObject* meanVertex,
                                    bool overrideBeamEstimation);
};

template <int NLayers>
using TrackingLoadPolicyN = TrackingLoadPolicy<constants::detIdFromNLayers<NLayers>(), NLayers>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTOR_TRAITS_H_ */
