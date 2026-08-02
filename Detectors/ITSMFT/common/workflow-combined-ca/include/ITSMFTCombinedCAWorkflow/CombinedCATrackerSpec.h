// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file CombinedCATrackerSpec.h
/// \brief Gate 4 C4: opt-in combined ITS+MFT common-CA tracker DPL device.
///        Wraps o2::itsmft::tracking::CombinedTimeFrameCoordinator, consuming
///        ITS and MFT cluster/pattern/ROF/dictionary(/label) inputs in one
///        task invocation and publishing both detectors' existing
///        CommonTrack OutputSpec streams (ITS exact, MFT float-projected)
///        through the already-approved CommonTrackOutputAdapter.h staging
///        functions. Lives in its own isolated library/executable -- no
///        link-graph overlap with ITSCAWorkflow, MFTWorkflow, or either
///        frozen legacy reco workflow (see this directory's CMakeLists.txt).

#ifndef ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_COMBINEDCATRACKERSPEC_H_
#define ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_COMBINEDCATRACKERSPEC_H_

#include <memory>

#include <vector>

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/CombinedTimeFrameCoordinator.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/TimeFrame.h"

namespace o2::itsmft::combined
{

/// Combined ITS+MFT opt-in common-CA tracker DPL task. Owns the one shared
/// TimeFrame and the one CombinedTimeFrameCoordinator both detectors'
/// tracking runs against; delegates all tracking to the coordinator and all
/// CommonTrack->detector-track output staging to CommonTrackOutputAdapter.h.
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
  void buildCoordinatorOnce();

  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  bool mUseMC = false;
  bool mCoordinatorBuilt = false;
  bool mContinuousReadoutChecked = false;
  const o2::itsmft::TopologyDictionary* mITSDict = nullptr;
  const o2::itsmft::TopologyDictionary* mMFTDict = nullptr;
  std::unique_ptr<o2::itsmft::tracking::ClusterDecoder> mITSDecoder;
  std::unique_ptr<o2::itsmft::tracking::ClusterDecoder> mMFTDecoder;
  std::unique_ptr<o2::itsmft::tracking::CombinedTimeFrameCoordinator> mCoordinator;
  o2::itsmft::tracking::TimeFrame mFrame;
  // The single-iteration TrackingParameters the coordinator was built with,
  // retained so run() can derive each detector's per-TF ROFTimingConfig
  // from TrackingParameters::AddTimeError without asking the coordinator
  // for its own construction-time arguments back.
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
