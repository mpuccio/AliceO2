// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file CombinedCATrackerSpec.h
/// \brief Combined ITS+MFT common-CA tracker DPL device. It consumes both
///        detector input sets and publishes their GenericTrack output streams
///        through GenericTrackOutputAdapter.h.
///
/// This DPL task owns one TimeFrame, one Tracker, their application
/// configuration, and the detector-specific publication context. trackFrame()
/// composes the workflow-owned source inputs,
/// loadTimeFrameSources(), Tracker::run(), and publication.

#ifndef ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_COMBINEDCATRACKERSPEC_H_
#define ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_COMBINEDCATRACKERSPEC_H_

#include <array>
#include <memory>
#include <optional>
#include <vector>

#include <gsl/span>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/GenericTrackOutputAdapter.h"
#include "ITSMFTTracking/detail/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITStracking/ROFLookupTables.h"

namespace o2::itsmft::combined
{

/// Combined ITS+MFT tracker DPL task. Owns one TimeFrame, one Tracker, one
/// application plan, and detector-specific publication context.
class CombinedCATrackerDPL : public o2::framework::Task
{
 public:
  CombinedCATrackerDPL(std::shared_ptr<o2::base::GRPGeomRequest> gr, bool useMC);
  ~CombinedCATrackerDPL() override = default;

  void init(o2::framework::InitContext& ic) final;
  void run(o2::framework::ProcessingContext& pc) final;
  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final;

 private:
  void updateTimeDependentParams(o2::framework::ProcessingContext& pc);
  // Built exactly once, lazily, the first time updateTimeDependentParams()
  // observes a valid magnetic field (TrackingMode::getTrackingParameters()
  // for MFT reads o2::base::Propagator::Instance()->getNominalBz(), same
  // ordering constraint used by standalone workflow initialisation).
  void buildParticipantsOnce();

  // Composes the atomic load and, once that has committed, the tracking phase
  // into a whole-event all-or-nothing contract. Any non-success -- a load
  // failure, a non-Success tracking outcome, or an exception from tracking
  // -- performs exactly one whole reset and leaves the
  // workflow-owned publication/timing state invalidated. A successful return
  // leaves the detector sidecars and publication exports populated for
  // run() to stage.
  o2::itsmft::tracking::TrackingOutcome trackFrame(const o2::itsmft::tracking::ClusterSourceInput& itsSource,
                                                   const o2::itsmft::tracking::ClusterSourceInput& mftSource,
                                                   const o2::InteractionRecord& origin);
  std::optional<o2::itsmft::tracking::LoadSourcesResult> validateSources(
    const o2::itsmft::tracking::ClusterSourceInput& itsSource,
    const o2::itsmft::tracking::ClusterSourceInput& mftSource) const noexcept;
  o2::itsmft::tracking::SurfaceCatalogView catalogView() const noexcept;
  std::optional<bool> dropTFUponFailureFor(o2::itsmft::tracking::ClusterSourceId source) const noexcept;
  void configureRofTables(const o2::itsmft::tracking::ClusterSourceInput& itsSource,
                          const o2::itsmft::tracking::ClusterSourceInput& mftSource);
  void clearRofViews() noexcept;
  void clearPublicationSidecars() noexcept;
  void invalidatePublication() noexcept;
  void markPublicationValid() noexcept;
  std::optional<o2::itsmft::tracking::GenericTrackPublicationExport> getITSPublicationExport() const;
  std::optional<o2::itsmft::tracking::GenericTrackPublicationExport> getMFTPublicationExport() const;

  const o2::itsmft::tracking::SurfaceTrackingScratch& getITSScratch() const noexcept { return mFrame.getWorkspace(); }
  const o2::itsmft::tracking::SurfaceTrackingScratch& getMFTScratch() const noexcept { return mFrame.getWorkspace(); }
  const o2::itsmft::tracking::ITSSharedClusterCompatibility& getITSSharedClusterCompatibility() const noexcept
  {
    return mITSCompatibility;
  }
  const o2::itsmft::tracking::MFTPublicationCompatibility& getMFTPublicationCompatibility() const noexcept
  {
    return mMFTCompatibility;
  }
  gsl::span<const o2::itsmft::tracking::SurfaceId> getITSOrderedSurfaces() const noexcept
  {
    const auto* binding = mFrame.getBinding(0);
    return binding == nullptr ? gsl::span<const o2::itsmft::tracking::SurfaceId>{}
                              : binding->getOrderedSurfaces().first(o2::itsmft::tracking::ITSNLayers);
  }
  gsl::span<const o2::itsmft::tracking::SurfaceId> getMFTOrderedSurfaces() const noexcept
  {
    const auto* binding = mFrame.getBinding(0);
    return binding == nullptr ? gsl::span<const o2::itsmft::tracking::SurfaceId>{}
                              : binding->getOrderedSurfaces().subspan(o2::itsmft::tracking::ITSNLayers,
                                                                      o2::itsmft::tracking::MFTNLayers);
  }

  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  bool mUseMC = false;
  bool mParticipantsBuilt = false;
  bool mContinuousReadoutChecked = false;
  const o2::itsmft::TopologyDictionary* mITSDict = nullptr;
  const o2::itsmft::TopologyDictionary* mMFTDict = nullptr;
  std::unique_ptr<o2::itsmft::tracking::ClusterDecoder> mITSDecoder;
  std::unique_ptr<o2::itsmft::tracking::ClusterDecoder> mMFTDecoder;
  o2::itsmft::tracking::TimeFrame mFrame;
  std::unique_ptr<o2::itsmft::tracking::Tracker> mTracker;
  std::unique_ptr<o2::itsmft::tracking::TrackerTraits> mTraits;
  o2::itsmft::tracking::DetectorPublicationAdapter<o2::itsmft::tracking::ITSNLayers> mITSPublicationAdapter;
  o2::itsmft::tracking::DetectorPublicationAdapter<o2::itsmft::tracking::MFTNLayers> mMFTPublicationAdapter;
  o2::itsmft::tracking::ITSSharedClusterCompatibility mITSCompatibility;
  o2::itsmft::tracking::MFTPublicationCompatibility mMFTCompatibility;
  o2::its::ROFOverlapTable<o2::itsmft::tracking::ITSNLayers> mITSROFOverlapTable;
  o2::its::ROFVertexLookupTable<o2::itsmft::tracking::ITSNLayers> mITSROFVertexLookupTable;
  o2::its::ROFMaskTable<o2::itsmft::tracking::ITSNLayers> mITSMultiplicityMask;
  o2::its::ROFMaskTable<o2::itsmft::tracking::ITSNLayers> mITSUPCMask;
  o2::its::ROFOverlapTable<o2::itsmft::tracking::MFTNLayers> mMFTROFOverlapTable;
  o2::its::ROFVertexLookupTable<o2::itsmft::tracking::MFTNLayers> mMFTROFVertexLookupTable;
  o2::its::ROFMaskTable<o2::itsmft::tracking::MFTNLayers> mMFTMultiplicityMask;
  o2::its::ROFMaskTable<o2::itsmft::tracking::MFTNLayers> mMFTUPCMask;
  std::optional<o2::itsmft::tracking::ClockTimingPublicationView> mITSClock;
  std::optional<o2::itsmft::tracking::ClockTimingPublicationView> mMFTClock;
  bool mPublicationValid = false;
};

/// Clusters are already aligned; the task requests nominal per-detector
/// geometry and does not expose an aligned-geometry switch.
o2::framework::DataProcessorSpec getCombinedCATrackerSpec(bool useMC);

} // namespace o2::itsmft::combined

#endif // ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_COMBINEDCATRACKERSPEC_H_
