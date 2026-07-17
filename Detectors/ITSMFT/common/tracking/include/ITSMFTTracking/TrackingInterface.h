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
/// \file TrackingInterface.h
/// \brief Shared DPL-facing CA tracking interface for ITS and MFT
///

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKINGINTERFACE_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKINGINTERFACE_H_

#include <limits>
#include <memory>
#include <vector>

#include <gsl/span>

#include "CommonDataFormat/IRFrame.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/DetectorTraits.h"
#ifndef GPUCA_GPUCODE
#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"
#endif
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/ROFLookupTables.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

namespace o2::itsmft::tracking
{

template <int NLayers>
class ITSMFTTrackingInterface
{
 public:
  static_assert(NLayers == ITSNLayers || NLayers == o2::mft::constants::mft::LayersNumber,
                "ITSMFTTrackingInterface supports ITS (7) and MFT (10) layer counts only");
  static constexpr o2::detectors::DetID::ID DetId = detIdFromNLayers<NLayers>();

  using TimeFrameN = TimeFrame<NLayers>;
  using TrackerN = Tracker<NLayers>;
  using TrackerTraitsN = TrackerTraits<NLayers>;
  using ROFOverlapTableN = o2::its::ROFOverlapTable<NLayers>;
  using ROFVertexLookupTableN = o2::its::ROFVertexLookupTable<NLayers>;
  using ROFMaskTableN = o2::its::ROFMaskTable<NLayers>;
  using BoundedMemoryResourceN = BoundedMemoryResource;

  ITSMFTTrackingInterface(bool useMC, o2::itsmft::TrackingMode::Type mode, bool overrideBeamEst);
#ifndef GPUCA_GPUCODE
  ITSMFTTrackingInterface(bool useMC, o2::itsmft::TrackingMode::Type mode, bool overrideBeamEst,
                          std::unique_ptr<DetectorSurfaceCatalogProvider> catalogProvider);

  void setDetectorSurfaceCatalogProvider(std::unique_ptr<DetectorSurfaceCatalogProvider> catalogProvider)
  {
    mTimeFrame.invalidateDetectorLayouts();
    mDetectorSurfaceCatalogProvider = std::move(catalogProvider);
  }
  DetectorLayoutSetBuildResult configureDetectorLayouts(gsl::span<const SurfaceId> orderedSurfaces,
                                                        TransitionPolicyTag policyTag)
  {
    return mTimeFrame.ensureDetectorLayouts(mDetectorSurfaceCatalogProvider.get(), orderedSurfaces, policyTag, mTrackParams);
  }
#endif

  void setTrackingMode(o2::itsmft::TrackingMode::Type mode) { mTrackingMode = mode; }
  void setClusterDictionary(const o2::itsmft::TopologyDictionary* dict) { mDict = dict; }
  void setMeanVertex(const o2::dataformats::MeanVertexObject* v) { mMeanVertex = v; }

  void initialise();

  /// Phase 2 (TimeFrame load) + Phase 3 (CA tracking). Returns tracker elapsed ms or -1 on failure.
  /// Does not wipe the TimeFrame; call clearTimeFrame() after outputs are extracted.
  float processTimeFrame(gsl::span<const o2::itsmft::ROFRecord> rofs,
                         gsl::span<const o2::itsmft::CompClusterExt> clusters,
                         gsl::span<const unsigned char> patterns,
                         const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels,
                         gsl::span<const o2::dataformats::IRFrame> irFrames = {});

  void clearTimeFrame() { mTimeFrame.wipe(); }

  TimeFrameN& getTimeFrame() { return mTimeFrame; }
  const TimeFrameN& getTimeFrame() const { return mTimeFrame; }
  const std::vector<o2::itsmft::TrackingParameters>& getTrackingParameters() const { return mTrackParams; }
  bool isActive() const { return !mTrackParams.empty(); }

 protected:
  virtual void onTimeFrameLoaded() {}
  virtual void onTrackingFinished(float elapsedMs) {}

 private:
  void resolveTrackingParameters();
  void initialiseMemoryPool();
  void initialiseTracker();
  void loadTimeFrame(gsl::span<const o2::itsmft::ROFRecord> rofs,
                     gsl::span<const o2::itsmft::CompClusterExt> clusters,
                     gsl::span<const unsigned char> patterns,
                     const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels,
                     gsl::span<const o2::dataformats::IRFrame> irFrames);
  float runTracking();
  void configureROFLookupTables();
  void configureBeamPosition();
  void configureTrackingTopology();
  void configureROFMask(gsl::span<const o2::itsmft::ROFRecord> rofs,
                        gsl::span<const o2::dataformats::IRFrame> irFrames);
  void validateROFInput(gsl::span<const o2::itsmft::ROFRecord> rofs) const;

  bool mUseMC = false;
  bool mOverrideBeamEstimation = false;
  o2::itsmft::TrackingMode::Type mTrackingMode = o2::itsmft::TrackingMode::Unset;
  std::vector<o2::itsmft::TrackingParameters> mTrackParams;
  std::shared_ptr<BoundedMemoryResourceN> mMemoryPool;
  std::unique_ptr<TrackerTraitsN> mTrackerTraits;
  std::unique_ptr<TrackerN> mTracker;
#ifndef GPUCA_GPUCODE
  std::unique_ptr<DetectorSurfaceCatalogProvider> mDetectorSurfaceCatalogProvider;
#endif
  const o2::itsmft::TopologyDictionary* mDict = nullptr;
  const o2::dataformats::MeanVertexObject* mMeanVertex = nullptr;
  TimeFrameN mTimeFrame;
  int mMFTROFrameLengthInBC = 0;
  bool mMFTTriggered = false;
};

using ITSMFTTrackingInterfaceITS = ITSMFTTrackingInterface<ITSNLayers>;
using ITSMFTTrackingInterfaceMFT = ITSMFTTrackingInterface<o2::mft::constants::mft::LayersNumber>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKINGINTERFACE_H_ */
