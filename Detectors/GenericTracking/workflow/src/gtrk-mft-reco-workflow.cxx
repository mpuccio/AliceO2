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

/// \file gtrk-mft-reco-workflow.cxx
/// \brief The MFT's tracking: O2's clusters in, O2's MFT tracks out, as
///        o2-mft-reco-workflow --clusters-from-upstream.
///
/// O2's cluster reader (unless `--clusters-from-upstream`), the converter into
/// our clusters, the measurement producer, the tracker (GTRK/MFT...), the
/// adaptor into O2's formats (unless `--disable-adaptor`) and O2's track
/// writer (unless `--disable-root-output`, or without the adaptor). With
/// `--vertices-from-upstream` the tracker takes the native vertex seeds
/// (GTRK/MFTVERTICES) from the workflow before it in the pipe.
/// With the tracking inside o2-gpu-reco-workflow (`--output-type
/// gtrk-mft-tracks`), this workflow twice: `--tracking-in-gpu-reco` before it
/// in the pipe (the converter and the producer), `--adaptor-only` after it.

#include "RecoWorkflow.h"

#include "CommonUtils/ConfigurableParam.h"
#include "Framework/CallbacksPolicy.h"
#include "Framework/CompletionPolicyHelpers.h"

using namespace o2::framework;

// DPL looks these up by symbol before runDataProcessing.h is included, so
// they have to keep external linkage.
// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::CallbacksPolicy>& policies) { o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies); }

// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::CompletionPolicy>& policies)
{
  // O2's writer in timeframe order, as o2-mft-reco-workflow has it.
  policies.push_back(CompletionPolicyHelpers::consumeWhenAllOrdered(".*(?:MFT|mft).*[W,w]riter.*"));
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions) { o2::gtrk::workflow::addRecoOptions(workflowOptions, true); }

#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Workflow/AdaptorSpec.h"

#include "ITSMFTCAWriter/MFTCATrackWriterSpec.h"
#include "ITSMFTWorkflow/ClusterReaderSpec.h"

#include <string>

#include "Framework/runDataProcessing.h" // IWYU pragma: keep

WorkflowSpec defineDataProcessing(const ConfigContext& context)
{
  o2::conf::ConfigurableParam::updateFromString(context.options().get<std::string>("configKeyValues"));

  namespace gw = o2::gtrk::workflow;
  using Det = o2::gtrk::detectors::MFT;
  const bool useMC = !context.options().get<bool>("disable-mc");
  const bool staggered = gw::staggeredInput<Det>(context);
  auto specs = gw::trackingChain<Det>(context, {o2::itsmft::getMFTClusterReaderSpec(useMC, staggered, true, false)});
  // The adaptor reads the tracks: after the GPU reconstruction when the tracking runs inside it.
  if (!context.options().get<bool>("disable-adaptor") && !gw::trackingInGPUReco(context)) {
    specs.push_back(gw::getAdaptorSpec<Det>(useMC, staggered));
    if (!context.options().get<bool>("disable-root-output")) {
      specs.push_back(o2::mft::getTrackWriterSpec(useMC));
    }
  }

  const o2::raw::HBFUtilsInitializer hbfIni(context, specs);

  return specs;
}
