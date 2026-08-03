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
/// M3 (GenericTrackingEngineMigration.md; ADR 0007) deletes the temporary
/// combined coordinator that used to stand between this device and its
/// owners: this device is now the sole owner of the one shared TimeFrame,
/// the one ITSMFTLegacyParticipantSet (the ITS/MFT application-layer
/// factory M2c introduced), and the one TrackingEngine.
/// run()'s own trackFrame() composes, in order: the set's atomic load
/// bindings, MultiSourceTimeFrameLoader::loadEvent(), TrackingEngine::
/// executeEvent(), and the set's publication exports/sidecars -- there is no
/// coordinator class standing between this device and those three owners.

#ifndef ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_COMBINEDCATRACKERSPEC_H_
#define ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_COMBINEDCATRACKERSPEC_H_

#include <memory>

#include <vector>

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ITSMFTLegacyParticipantSet.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingEngine.h"
#include "ITSMFTTracking/TrackingParticipant.h"

namespace o2::itsmft::combined
{

/// Combined ITS+MFT opt-in common-CA tracker DPL task. Owns the one shared
/// TimeFrame, the one ITSMFTLegacyParticipantSet both detectors' tracking
/// runs against, and the one TrackingEngine; delegates all CommonTrack->
/// detector-track output staging to CommonTrackOutputAdapter.h. ITS/MFT
/// source-position and failure-classification knowledge lives here (this
/// application adapter) and inside ITSMFTLegacyParticipantSet -- never
/// inside TrackingEngine/TrackingParticipant/TimeFrame/
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

  // Composes the atomic load (via mParticipants' own source-qualified
  // bindings and MultiSourceTimeFrameLoader::loadEvent()) and, once that has
  // committed, the tracking phase (TrackingEngine::executeEvent() over
  // mParticipants.schedule()) into the same whole-event all-or-nothing
  // contract this workflow's pre-M3 combined coordinator used to apply: any
  // non-success -- a load failure, a non-Success tracking outcome, or an
  // exception from either leg's track() -- performs exactly one whole
  // reset (TrackingEngine::resetEvent() on a load failure; executeEvent()
  // already applies it internally on a tracking failure) and leaves the
  // publication/timing bridge invalidated. A successful return leaves
  // mParticipants' publication exports/sidecars populated for run() to
  // stage. This is application-adapter composition, not a new core owner:
  // the fixed ITS=0/MFT=1 source-position mapping and per-detector
  // DropTFUponFailure classification both still live only inside
  // mParticipants, never here or in TrackingEngine.
  o2::itsmft::tracking::ParticipantOutcome trackFrame(const o2::itsmft::tracking::ClusterSourceInput& itsSource,
                                                      const o2::itsmft::tracking::ClusterSourceInput& mftSource,
                                                      const o2::InteractionRecord& origin);

  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  bool mUseMC = false;
  bool mParticipantsBuilt = false;
  bool mContinuousReadoutChecked = false;
  const o2::itsmft::TopologyDictionary* mITSDict = nullptr;
  const o2::itsmft::TopologyDictionary* mMFTDict = nullptr;
  std::unique_ptr<o2::itsmft::tracking::ClusterDecoder> mITSDecoder;
  std::unique_ptr<o2::itsmft::tracking::ClusterDecoder> mMFTDecoder;
  o2::itsmft::tracking::TimeFrame mFrame;
  std::unique_ptr<o2::itsmft::tracking::ITSMFTLegacyParticipantSet> mParticipants;
  o2::itsmft::tracking::TrackingEngine mEngine;
  // The single-iteration TrackingParameters mParticipants was built with,
  // retained so run() can derive each detector's per-TF ROFTimingConfig
  // from TrackingParameters::AddTimeError without asking the set for its
  // own construction-time arguments back.
  o2::itsmft::TrackingParameters mITSTrackingParams;
  o2::itsmft::TrackingParameters mMFTTrackingParams;
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
