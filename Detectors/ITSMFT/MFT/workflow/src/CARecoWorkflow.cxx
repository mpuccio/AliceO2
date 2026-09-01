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

/// @file CARecoWorkflow.cxx

#include "MFTWorkflow/CARecoWorkflow.h"

#include "GlobalTrackingWorkflowReaders/IRFrameReaderSpec.h"
#include "ITSMFTCAWriter/MFTCATrackWriterSpec.h"
#include "ITSMFTWorkflow/ClustererSpec.h"
#include "ITSMFTWorkflow/ClusterWriterSpec.h"
#include "ITSMFTWorkflow/DigitReaderSpec.h"
#include "MFTTracking/MFTTrackingParam.h"
#include "MFTWorkflow/CATrackerSpec.h"
#include "MFTWorkflow/MFTAssessmentSpec.h"
#include "MFTWorkflow/TracksToRecordsSpec.h"

namespace o2::mft::ca_reco_workflow
{

framework::WorkflowSpec getWorkflow(
  bool useMC,
  bool doStag,
  bool useGeom,
  bool useIRFrames,
  bool upstreamDigits,
  bool upstreamClusters,
  bool clrofOnly,
  bool disableRootOutput,
  bool runAssessment,
  bool processGen,
  bool runTracking,
  o2::itsmft::TrackingMode::Type trackingMode,
  bool runTracks2Records)
{
  framework::WorkflowSpec specs;

  if (!(upstreamDigits || upstreamClusters)) {
    specs.emplace_back(o2::itsmft::getMFTDigitReaderSpec(useMC, doStag, false, true, "mftdigits.root"));
    const auto& trackingParam = MFTTrackingParam::Instance();
    if (useIRFrames || trackingParam.irFramesOnly) {
      specs.emplace_back(o2::globaltracking::getIRFrameReaderSpec("ITS", 0, "its-irframe-reader", "o2_its_irframe.root"));
    }
  }
  if (!upstreamClusters) {
    specs.emplace_back(o2::itsmft::getMFTClustererSpec(useMC, doStag));
  }
  if (!disableRootOutput || clrofOnly) {
    specs.emplace_back(o2::itsmft::getMFTClusterWriterSpec(useMC, doStag, clrofOnly));
  }

  if (runTracking) {
    specs.emplace_back(o2::mft::getCATrackerSpec(useMC, useGeom, useIRFrames, trackingMode));
    if (!disableRootOutput) {
      specs.emplace_back(o2::mft::getTrackWriterSpec(useMC, true));
    }
    if (runAssessment) {
      specs.emplace_back(o2::mft::getMFTAssessmentSpec(useMC, useGeom, processGen));
    }
    if (runTracks2Records) {
      specs.emplace_back(o2::mft::getTracksToRecordsSpec());
    }
  }
  return specs;
}

} // namespace o2::mft::ca_reco_workflow
