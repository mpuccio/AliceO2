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
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/Cell.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Tracklet.h"

namespace o2::itsmft::tracking
{

/// Per-detector differences in refit, track acceptance, and index-table setup.
/// Everything else stays in TrackerTraits and matches ITS line-for-line.
template <int NLayers>
struct DetectorTraits {
  using TrackType = CATrackType<NLayers>;
  using TrackSeedN = o2::itsmft::tracking::TrackSeedN<NLayers>;
  using CellSeedN = o2::itsmft::tracking::CellSeedN<NLayers>;
  using TimeFrameN = TimeFrame<NLayers>;
  static constexpr o2::detectors::DetID::ID DetId = detIdFromNLayers<NLayers>();

  static bool refitSeed(const TrackSeedN& seed,
                        TrackType& track,
                        const TrackingParameters& params,
                        float bz,
                        TimeFrameN& tf,
                        const o2::its::TrackingFrameInfo* const tfInfos[NLayers],
                        const o2::its::Cluster* const unsortedClusters[NLayers],
                        const o2::base::PropagatorImpl<float>* propagator);

  static void configureIndexTableUtils(IndexTableUtils<NLayers>& utils, const TrackingParameters& params);

};

template <o2::detectors::DetID::ID DetId, int NLayers>
struct TrackingLoadPolicy {
  static void configureBeamPosition(TimeFrame<NLayers>& tf,
                                    const TrackingParameters& p,
                                    const o2::dataformats::MeanVertexObject* meanVertex,
                                    bool overrideBeamEstimation);
};

template <int NLayers>
using TrackingLoadPolicyN = TrackingLoadPolicy<detIdFromNLayers<NLayers>(), NLayers>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTOR_TRAITS_H_ */
