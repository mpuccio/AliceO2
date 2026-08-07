// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file itsmft-combined-ca-tracker-workflow.cxx
/// \brief Combined ITS+MFT common-CA tracker workflow. Sync mode only, with
///        tracker outputs and no IR-frame masking.

#include <string>
#include <vector>

#include "CommonUtils/ConfigurableParam.h"
#include "DetectorsRaw/HBFUtilsInitializer.h"
#include "Framework/CallbacksPolicy.h"
#include "Framework/CompletionPolicyHelpers.h"
#include "Framework/ConfigParamSpec.h"
#include "ITSMFTCAWriter/ITSCATrackWriterSpec.h"
#include "ITSMFTCAWriter/MFTCATrackWriterSpec.h"
#include "ITSMFTCombinedCAWorkflow/CombinedCATrackerSpec.h"
#include "ITSMFTCombinedCAWorkflow/ConfigPreflight.h"
#include "ITSMFTTracking/Configuration.h"

using namespace o2::framework;

void customize(std::vector<o2::framework::CallbacksPolicy>& policies)
{
  o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies);
}

void customize(std::vector<o2::framework::CompletionPolicy>& policies)
{
  policies.push_back(CompletionPolicyHelpers::consumeWhenAllOrdered(".*(?:ITSMFT|itsmft).*[W,w]riter.*"));
}

void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{"disable-mc", VariantType::Bool, false, {"disable MC labels"}});
  workflowOptions.push_back(ConfigParamSpec{"disable-root-output", VariantType::Bool, false, {"do not write output root files"}});
  workflowOptions.push_back(ConfigParamSpec{"tracking-mode", VariantType::String, "sync", {"only 'sync' is supported by this opt-in workflow"}});
  workflowOptions.push_back(ConfigParamSpec{"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings (must include "
                                                                                         "ITSMFTCombinedCATrackerParam.enabled=true;ITSCommonCATrackerParam.useDiamond=true)"}});
  o2::raw::HBFUtilsInitializer::addConfigOption(workflowOptions);
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{
  // Apply configuration and enforce all preflight checks before constructing
  // any DataProcessorSpec/device.
  o2::conf::ConfigurableParam::updateFromString(config.options().get<std::string>("configKeyValues"));
  o2::itsmft::combined::requireCombinedTrackingEnabledOrFatal();

  const auto trMode = o2::itsmft::TrackingMode::fromString(config.options().get<std::string>("tracking-mode"));
  o2::itsmft::combined::requireSyncTrackingModeOrFatal(trMode);
  o2::itsmft::combined::requireDiamondVertexConstraintOrFatal();
  o2::itsmft::combined::requireNoMFTIRFrameConfigOrFatal();

  const bool useMC = !config.options().get<bool>("disable-mc");
  const bool disableRootOutput = config.options().get<bool>("disable-root-output");

  WorkflowSpec specs;
  specs.emplace_back(o2::itsmft::combined::getCombinedCATrackerSpec(useMC));
  if (!disableRootOutput) {
    // Reuse the shared writer specs; the combined MFT route always writes the
    // seed-pattern branch produced by the CA adapter.
    specs.emplace_back(o2::its::ca::getTrackWriterSpec(useMC));
    specs.emplace_back(o2::mft::getTrackWriterSpec(useMC, true));
  }

  o2::raw::HBFUtilsInitializer hbfIni(config, specs);

  return specs;
}
