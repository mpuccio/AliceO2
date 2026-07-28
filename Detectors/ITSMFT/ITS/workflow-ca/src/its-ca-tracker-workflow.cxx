// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   its-ca-tracker-workflow.cxx
/// \brief  Opt-in ITS common-CA tracker workflow (Gate 3 workflow-onboarding
///         Slice 2). Runs Tracker<7>/ITSMFTTrackingInterface<7> on real ITS
///         cluster inputs, Sync mode only, tracker-only outputs. The frozen
///         legacy o2-its-reco-workflow is untouched; this is a separate
///         executable built from an isolated library with no link-graph
///         overlap with ITSWorkflow.

#include <tuple>
#include <vector>

#include "CommonUtils/ConfigurableParam.h"
#include "DetectorsRaw/HBFUtilsInitializer.h"
#include "Framework/CallbacksPolicy.h"
#include "Framework/CompletionPolicyHelpers.h"
#include "Framework/ConfigParamSpec.h"
#include "ITSCAWorkflow/CATrackerSpec.h"
#include "ITSCAWorkflow/ConfigPreflight.h"
#include "ITSCAWorkflow/TrackWriterSpec.h"
#include "ITSMFTTracking/Configuration.h"

using namespace o2::framework;

void customize(std::vector<o2::framework::CallbacksPolicy>& policies)
{
  o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies);
}

void customize(std::vector<o2::framework::CompletionPolicy>& policies)
{
  policies.push_back(CompletionPolicyHelpers::consumeWhenAllOrdered(".*(?:ITS|its).*[W,w]riter.*"));
}

void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{"disable-mc", VariantType::Bool, false, {"disable MC labels"}});
  workflowOptions.push_back(ConfigParamSpec{"disable-root-output", VariantType::Bool, false, {"do not write output root files"}});
  workflowOptions.push_back(ConfigParamSpec{"use-geom", VariantType::Bool, false, {"use geometry from the global geometry manager"}});
  workflowOptions.push_back(ConfigParamSpec{"tracking-mode", VariantType::String, "sync", {"only 'sync' is supported by this opt-in workflow"}});
  workflowOptions.push_back(ConfigParamSpec{"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings (e.g. ITSCommonCATrackerParam.useDiamond=true)"}});
  o2::raw::HBFUtilsInitializer::addConfigOption(workflowOptions);
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{
  // Preflight order, all before any DataProcessorSpec/device is constructed:
  //   1. reject a legacy ITSCATrackerParam.* --configKeyValues override, then
  //      apply the (accepted) string;
  //   2. reject every tracking-mode other than Sync;
  //   3. reject a configuration with no reproducible vertex/beam constraint.
  o2::its::ca::applyConfigKeyValuesOrFatal(config.options().get<std::string>("configKeyValues"));

  const auto trMode = o2::itsmft::TrackingMode::fromString(config.options().get<std::string>("tracking-mode"));
  o2::its::ca::requireSyncTrackingModeOrFatal(trMode);
  o2::its::ca::requireDiamondVertexConstraintOrFatal();

  const bool useMC = !config.options().get<bool>("disable-mc");
  const bool disableRootOutput = config.options().get<bool>("disable-root-output");
  const bool useGeom = config.options().get<bool>("use-geom");

  WorkflowSpec specs;
  specs.emplace_back(o2::its::ca::getCATrackerSpec(useMC, useGeom, trMode));
  if (!disableRootOutput) {
    specs.emplace_back(o2::its::ca::getTrackWriterSpec(useMC));
  }

  o2::raw::HBFUtilsInitializer hbfIni(config, specs);

  return specs;
}
