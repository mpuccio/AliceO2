// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file CombinedCATrackerSpec.h
/// \brief Gate 4 C4/M3: opt-in combined ITS+MFT common-CA tracker DPL
///        device. Consuming ITS and MFT cluster/pattern/ROF/dictionary
///        (/label) inputs in one task invocation and publishing both
///        detectors' existing CommonTrack OutputSpec streams (ITS exact, MFT
///        float-projected) through the already-approved
///        CommonTrackOutputAdapter.h staging functions. Lives in its own
///        isolated library/executable -- no link-graph overlap with
///        ITSCAWorkflow, MFTWorkflow, or either frozen legacy reco workflow
///        (see this directory's CMakeLists.txt).
///
/// This DPL task owns the shared TimeFrame, the two non-owning Tracker
/// components, their application configuration, and the event publication
/// context. trackFrame() composes the workflow-owned source inputs,
/// MultiSourceTimeFrameLoader::load(), Tracker::run(), and publication.

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
#include "ITSMFTTracking/CommonTrackOutputAdapter.h"
#include "ITSMFTTracking/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingOperationAdapter.h"
#include "ITStracking/ROFLookupTables.h"

namespace o2::itsmft::combined
{

/// Combined ITS+MFT opt-in common-CA tracker DPL task. Owns the one shared
/// TimeFrame, the two Tracker components, the combined application plan, and
/// the event publication context; delegates
/// all CommonTrack->detector-track output staging to
/// CommonTrackOutputAdapter.h. ITS/MFT source-position, invocation order, failure-
/// classification, and publication state live here in this workflow task --
/// never inside Tracker/TimeFrame/MultiSourceTimeFrameLoader, which stay
/// detector-neutral.
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
  // failure, a non-Success tracking outcome, or an exception from either
  // leg's track() -- performs exactly one whole reset and leaves the
  // workflow-owned publication/timing state invalidated. A successful return
  // leaves the two participant sidecars and publication exports populated for
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
  std::optional<o2::itsmft::tracking::CommonTrackPublicationExport> getITSPublicationExport() const;
  std::optional<o2::itsmft::tracking::CommonTrackPublicationExport> getMFTPublicationExport() const;

  const o2::itsmft::tracking::SurfaceTrackingScratch& getITSScratch() const noexcept { return mFrame.getWorkspace(o2::itsmft::tracking::ClusterSourceId{0}); }
  const o2::itsmft::tracking::SurfaceTrackingScratch& getMFTScratch() const noexcept { return mFrame.getWorkspace(o2::itsmft::tracking::ClusterSourceId{1}); }
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
    const auto* binding = mFrame.getBinding(0, o2::itsmft::tracking::ClusterSourceId{0});
    return binding == nullptr ? gsl::span<const o2::itsmft::tracking::SurfaceId>{} : binding->getOrderedSurfaces();
  }
  gsl::span<const o2::itsmft::tracking::SurfaceId> getMFTOrderedSurfaces() const noexcept
  {
    const auto* binding = mFrame.getBinding(0, o2::itsmft::tracking::ClusterSourceId{1});
    return binding == nullptr ? gsl::span<const o2::itsmft::tracking::SurfaceId>{} : binding->getOrderedSurfaces();
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
  std::unique_ptr<o2::itsmft::tracking::Tracker> mITSTracker;
  std::unique_ptr<o2::itsmft::tracking::Tracker> mMFTTracker;
  std::unique_ptr<o2::itsmft::tracking::TrackerTraits> mITSTraits;
  std::unique_ptr<o2::itsmft::tracking::TrackerTraits> mMFTTraits;
  std::unique_ptr<o2::itsmft::tracking::TrackingOperationAdapter> mITSOperationAdapter;
  std::unique_ptr<o2::itsmft::tracking::TrackingOperationAdapter> mMFTOperationAdapter;
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

/// useGeom is deliberately not a parameter: clusters entering this tracker
/// are already aligned, and this workflow never selects
/// o2::base::GRPGeomRequest::GeomRequest::Aligned -- it always requests
/// GeomRequest::None plus explicit ITS/Config/Geometry and MFT/Config/Geometry
/// condition objects (the same nominal/ideal per-detector GeometryTGeo
/// snapshot both single-detector opt-in workflows already default to).
o2::framework::DataProcessorSpec getCombinedCATrackerSpec(bool useMC);

} // namespace o2::itsmft::combined

#endif // ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_COMBINEDCATRACKERSPEC_H_
