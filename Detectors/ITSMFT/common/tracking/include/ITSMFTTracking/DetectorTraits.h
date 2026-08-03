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

#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/LegacyTrackerScratch.h"
#include "ITSMFTTracking/SurfaceCatalogView.h"
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
  using ScratchN = LegacyTrackerScratch<NLayers>;
  static constexpr o2::detectors::DetID::ID DetId = detIdFromNLayers<NLayers>();

  // M5d: both branches now go through the shared, descriptor-driven native
  // refit (Propagator.h / NativeRefitDriver.h) instead of a frozen legacy
  // Kalman fitter -- see doc/decisions/0008-native-refit-activation.md. Both
  // branches read layerMeasurements/surfaceCatalog; neither reads a raw
  // o2::its::TrackingFrameInfo/Cluster array or an o2::base::Propagator
  // instance any longer. `scratch` is retained only for cluster
  // external-index/size bookkeeping (output metadata, never a physical
  // read) -- see MFTFwdTrackHelpers.h. `surfaceCatalog` is this iteration's
  // TrackerTraits::mTraversalLayout.getSurfaceCatalogView() (ADR 0001
  // nominal material, resolved the same way every other native operation in
  // this library already resolves it).
  //
  // expectedSource (Gate 4 C2 source-identity correction): the
  // ClusterSourceId this call's caller (TrackerTraits::findRoadsForPolicy())
  // resolved once for this invocation -- mBinding->getSource() when a
  // DetectorTraversalBinding is adopted, ClusterSourceId{0} otherwise. Same
  // MFT-only-consumer shape as layerMeasurements: forwarded to
  // refitTrackFwd()'s final-refit identity re-check; the barrel/ITS branch
  // (refitSeedITS) neither receives nor needs it, since it never reads a
  // normalized SurfaceMeasurement at all.
  static bool refitSeed(const TrackSeedN& seed,
                        TrackType& track,
                        const TrackingParameters& params,
                        float bz,
                        ScratchN& scratch,
                        const LayerMeasurementSpans<NLayers>& layerMeasurements,
                        SurfaceCatalogView surfaceCatalog,
                        ClusterSourceId expectedSource);

  static void copySeedPatternToTrack(TrackType& track, const TrackSeedN& seed) noexcept;
  static void clearTransientLayerPattern(TrackType& track) noexcept;
  static bool haveSamePolarity(const TrackType& a, const TrackType& b) noexcept;
};

template <o2::detectors::DetID::ID DetId, int NLayers>
struct TrackingLoadPolicy {
  // Beam position is TimeFrame-owned (shared, detector-neutral) state --
  // this takes the non-templated TimeFrame directly, not a per-detector
  // LegacyTrackerScratch<NLayers>.
  static void configureBeamPosition(TimeFrame& frame,
                                    const TrackingParameters& p,
                                    const o2::dataformats::MeanVertexObject* meanVertex,
                                    bool overrideBeamEstimation);
};

template <int NLayers>
using TrackingLoadPolicyN = TrackingLoadPolicy<detIdFromNLayers<NLayers>(), NLayers>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTOR_TRAITS_H_ */
