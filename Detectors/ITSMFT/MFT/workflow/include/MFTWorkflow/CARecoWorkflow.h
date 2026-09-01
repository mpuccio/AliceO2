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

#ifndef O2_MFT_CARECOWORKFLOW_H_
#define O2_MFT_CARECOWORKFLOW_H_

/// @file CARecoWorkflow.h

#include "Framework/WorkflowSpec.h"
#include "ITSMFTTracking/Configuration.h"

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
  bool runTracks2Records);

} // namespace o2::mft::ca_reco_workflow

#endif // O2_MFT_CARECOWORKFLOW_H_
