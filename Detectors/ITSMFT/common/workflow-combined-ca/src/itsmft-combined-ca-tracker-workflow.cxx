// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file itsmft-combined-ca-tracker-workflow.cxx
/// \brief Gate 4 C4: opt-in combined ITS+MFT common-CA tracker workflow.
///        Sync mode only, tracker-only outputs (no vertex products). No
///        --use-geom (never selects an aligned-geometry route) and no
///        --use-irframes (the combined coordinator has no IR-frame/trigger
///        masking; see ConfigPreflight.h). The frozen legacy
///        o2-its-reco-workflow/o2-mft-reco-workflow and the single-detector
///        opt-in o2-its-ca-tracker-workflow/o2-mft-ca-tracker-workflow are
///        all untouched; this is a separate executable built from an
///        isolated library with no link-graph overlap with any of them.

#include <string>
#include <vector>

#include "CommonUtils/ConfigurableParam.h"
#include "DetectorsRaw/HBFUtilsInitializer.h"
#include "Framework/CallbacksPolicy.h"
#include "Framework/CompletionPolicyHelpers.h"
#include "Framework/ConfigParamSpec.h"
#include "ITSMFTCombinedCAWorkflow/CombinedCATrackerSpec.h"
#include "ITSMFTCombinedCAWorkflow/ConfigPreflight.h"
#include "ITSMFTCombinedCAWorkflow/ITSCommonTrackWriterSpec.h"
#include "ITSMFTCombinedCAWorkflow/MFTCommonTrackWriterSpec.h"
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
  // Preflight order, all before any DataProcessorSpec/device is constructed:
  //   1. apply --configKeyValues;
  //   2. reject unless explicitly opted in;
  //   3. reject every tracking-mode other than Sync;
  //   4. reject a configuration with no reproducible ITS vertex/beam constraint;
  //   5. reject a configuration requesting MFT IR-frame/triggered masking
  //      this coordinator cannot apply.
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
    specs.emplace_back(o2::itsmft::combined::getITSCommonTrackWriterSpec(useMC));
    specs.emplace_back(o2::itsmft::combined::getMFTCommonTrackWriterSpec(useMC));
  }

  o2::raw::HBFUtilsInitializer hbfIni(config, specs);

  return specs;
}
