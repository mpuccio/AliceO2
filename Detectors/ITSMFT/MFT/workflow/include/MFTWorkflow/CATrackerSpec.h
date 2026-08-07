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

/// @file   CATrackerSpec.h

#ifndef O2_MFT_CATRACKERSPEC_H_
#define O2_MFT_CATRACKERSPEC_H_

#include <memory>
#include <optional>

#include "DetectorsBase/GRPGeomHelper.h"
#include "CommonDataFormat/IRFrame.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "ITSMFTTracking/CommonTrackOutputAdapter.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/detail/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TrackingOperationAdapter.h"
#include "ITSMFTTracking/ROFViews.h"
#include "ITStracking/ROFLookupTables.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::mft
{

/// Gate 3 failure-contract publication decision for CATrackerDPL::run(),
/// factored out as a pure function of (tracker active, tracking result) so
/// the publish-vs-skip contract can be exercised by a unit test without a
/// DPL ProcessingContext.
enum class CATrackerPublicationAction {
  PublishInactiveEmpty, ///< tracker not configured/active: publish the existing echoed-empty outputs
  PublishActiveResult,  ///< active tracking returned a non-dropped result (including valid-empty input): publish tracker outputs
  SkipDroppedTimeFrame, ///< active tracking recoverably dropped this TF: publish nothing, keep the device running
};

CATrackerPublicationAction decideCATrackerPublicationAction(bool trackerActive, o2::itsmft::tracking::TrackingOutcome outcome) noexcept;

/// MFT CA tracker DPL task. Owns the standalone TimeFrame and composes the
/// workflow-owned input/timing/publication edge with Tracker.
class CATrackerDPL : public o2::framework::Task
{
 public:
  CATrackerDPL(std::shared_ptr<o2::base::GRPGeomRequest> gr,
               bool useMC,
               o2::itsmft::TrackingMode::Type trMode);
  ~CATrackerDPL() override = default;

  void init(framework::InitContext& ic) final;
  void run(framework::ProcessingContext& pc) final;
  void finaliseCCDB(framework::ConcreteDataMatcher& matcher, void* obj) final;

 private:
  void updateTimeDependentParams(framework::ProcessingContext& pc);
  void configureROFViews(gsl::span<const o2::itsmft::ROFRecord> rofs,
                         gsl::span<const o2::dataformats::IRFrame> irFrames);
  void invalidatePublication() noexcept;
  void initialiseTracking();
  o2::itsmft::tracking::TrackingOutcome processTimeFrame(
    gsl::span<const o2::itsmft::ROFRecord> rofs,
    gsl::span<const o2::itsmft::CompClusterExt> clusters,
    gsl::span<const unsigned char> patterns,
    const o2::dataformats::MCTruthContainer<MCCompLabel>* labels,
    gsl::span<const o2::dataformats::IRFrame> irFrames);
  void resetEvent() noexcept;
  bool isActive() const noexcept { return mFrame.isConfigured() && !mFrame.getTrackingParameters().empty(); }
  const o2::itsmft::tracking::SurfaceTrackingScratch& getScratch() const noexcept
  {
    return mFrame.getWorkspace(o2::itsmft::tracking::ClusterSourceId{0});
  }
  o2::itsmft::tracking::SurfaceTrackingScratch& getScratch() noexcept
  {
    return mFrame.getWorkspace(o2::itsmft::tracking::ClusterSourceId{0});
  }

  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  bool mUseMC = false;
  bool mTrackingInitialised = false;
  o2::itsmft::TrackingMode::Type mTrackingMode = o2::itsmft::TrackingMode::Unset;
  o2::itsmft::tracking::TimeFrame mFrame;
  std::unique_ptr<o2::itsmft::tracking::TrackerTraits> mTrackerTraits;
  std::unique_ptr<o2::itsmft::tracking::Tracker> mTracker;
  std::unique_ptr<o2::itsmft::tracking::TrackingOperationAdapter> mOperationAdapter;
  std::unique_ptr<o2::itsmft::tracking::ClusterDecoder> mClusterDecoder;
  const o2::itsmft::TopologyDictionary* mDictionary = nullptr;
  o2::itsmft::tracking::DetectorPublicationAdapter<o2::itsmft::tracking::MFTNLayers> mPublicationAdapter;
  o2::itsmft::tracking::MFTPublicationCompatibility mCompatibility;
  o2::its::ROFOverlapTable<o2::itsmft::tracking::MFTNLayers> mROFOverlapTable;
  o2::its::ROFVertexLookupTable<o2::itsmft::tracking::MFTNLayers> mROFVertexLookupTable;
  o2::its::ROFMaskTable<o2::itsmft::tracking::MFTNLayers> mMultiplicityMask;
  o2::its::ROFMaskTable<o2::itsmft::tracking::MFTNLayers> mUPCMask;
  std::optional<o2::itsmft::tracking::ClockTimingPublicationView> mPublicationClock;
  int mMFTROFrameLengthInBC = 0;
};

o2::framework::DataProcessorSpec getCATrackerSpec(bool useMC, bool useGeom, bool useIRFrames,
                                                  o2::itsmft::TrackingMode::Type trMode);

} // namespace o2::mft

#endif // O2_MFT_CATRACKERSPEC_H_
