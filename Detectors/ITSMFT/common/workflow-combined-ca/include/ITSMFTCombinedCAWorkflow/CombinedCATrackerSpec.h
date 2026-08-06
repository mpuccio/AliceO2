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
/// M6g: this DPL task is the application owner of the shared TimeFrame, the
/// two concrete plan-driven participants, their combined application plan,
/// the event publication context, and the TrackingEngine. The common tracking
/// library contains no ITS+MFT coordinator or event-loop state. run()'s own
/// trackFrame() composes, in order: the workflow-owned atomic load bindings,
/// MultiSourceTimeFrameLoader::loadEvent(), TrackingEngine::executeEvent(),
/// and publication staging.

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
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/ClockTimingPublicationView.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/SurfacePlanTrackingParticipant.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingEngine.h"
#include "ITSMFTTracking/TrackingInterface.h"
#include "ITSMFTTracking/TrackingParticipant.h"

namespace o2::itsmft::combined
{

/// Combined ITS+MFT opt-in common-CA tracker DPL task. Owns the one shared
/// TimeFrame, both concrete plan-driven participants, the combined application
/// plan, the event publication context, and the one TrackingEngine; delegates
/// all CommonTrack->detector-track output staging to
/// CommonTrackOutputAdapter.h. ITS/MFT source-position, schedule, failure-
/// classification, and publication state live here in this workflow task --
/// never inside TrackingEngine/TrackingParticipant/TimeFrame/
/// MultiSourceTimeFrameLoader, which stay detector-neutral.
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
  // ordering constraint ITSMFTTrackingInterface::initialise() already has).
  void buildParticipantsOnce();

  // Composes the atomic load and, once that has committed, the tracking phase
  // into a whole-event all-or-nothing contract. Any non-success -- a load
  // failure, a non-Success tracking outcome, or an exception from either
  // leg's track() -- performs exactly one whole reset and leaves the
  // workflow-owned publication/timing state invalidated. A successful return
  // leaves the two participant sidecars and publication exports populated for
  // run() to stage.
  o2::itsmft::tracking::ParticipantOutcome trackFrame(const o2::itsmft::tracking::ClusterSourceInput& itsSource,
                                                      const o2::itsmft::tracking::ClusterSourceInput& mftSource,
                                                      const o2::InteractionRecord& origin);

  gsl::span<o2::itsmft::tracking::TrackingParticipant* const> schedule() noexcept { return mSchedule; }
  std::optional<o2::itsmft::tracking::LoadSourcesResult> validateSources(
    const o2::itsmft::tracking::ClusterSourceInput& itsSource,
    const o2::itsmft::tracking::ClusterSourceInput& mftSource) const noexcept;
  std::array<o2::itsmft::tracking::MultiSourceTimeFrameLoader::AtomicLoadBinding, 2> loadBindings(
    const o2::itsmft::tracking::ClusterSourceInput& itsSource,
    const o2::itsmft::tracking::ClusterSourceInput& mftSource) noexcept;
  o2::itsmft::tracking::SurfaceCatalogView catalogView() const noexcept;
  std::optional<bool> dropTFUponFailureFor(o2::itsmft::tracking::ClusterSourceId source) const noexcept;
  void configureRofTables(const o2::itsmft::tracking::ClusterSourceInput& itsSource,
                          const o2::itsmft::tracking::ClusterSourceInput& mftSource);
  void clearPublicationSidecars() noexcept;
  void invalidatePublication() noexcept;
  void markPublicationValid() noexcept;
  std::optional<o2::itsmft::tracking::CommonTrackPublicationExport> getITSPublicationExport() const;
  std::optional<o2::itsmft::tracking::CommonTrackPublicationExport> getMFTPublicationExport() const;

  const o2::itsmft::tracking::SurfaceTrackingScratch& getITSScratch() const noexcept { return mITSParticipant->getScratch(); }
  const o2::itsmft::tracking::SurfaceTrackingScratch& getMFTScratch() const noexcept { return mMFTParticipant->getScratch(); }
  const o2::itsmft::tracking::ITSSharedClusterCompatibility& getITSSharedClusterCompatibility() const noexcept
  {
    return *mITSParticipant->getITSSharedClusterCompatibility();
  }
  const o2::itsmft::tracking::MFTPublicationCompatibility& getMFTPublicationCompatibility() const noexcept
  {
    return *mMFTParticipant->getMFTPublicationCompatibility();
  }
  gsl::span<const o2::itsmft::tracking::SurfaceId> getITSOrderedSurfaces() const noexcept
  {
    return mITSParticipant->ownedSurfaces();
  }
  gsl::span<const o2::itsmft::tracking::SurfaceId> getMFTOrderedSurfaces() const noexcept
  {
    return mMFTParticipant->ownedSurfaces();
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
  std::unique_ptr<o2::itsmft::tracking::SurfacePlanTrackingParticipantITS> mITSParticipant;
  std::unique_ptr<o2::itsmft::tracking::SurfacePlanTrackingParticipantMFT> mMFTParticipant;
  std::array<o2::itsmft::tracking::TrackingParticipant*, 2> mSchedule{};
  std::optional<o2::itsmft::tracking::ClockTimingPublicationView> mITSClock;
  std::optional<o2::itsmft::tracking::ClockTimingPublicationView> mMFTClock;
  bool mPublicationValid = false;
  o2::itsmft::tracking::TrackingEngine mEngine;
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
