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

/// @file   mft-ca-tracker-workflow.cxx

#include <vector>
#include <utility>

#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsITSMFT/DPLAlpideParamInitializer.h"
#include "DetectorsRaw/HBFUtilsInitializer.h"
#include "Framework/CallbacksPolicy.h"
#include "Framework/CompletionPolicyHelpers.h"
#include "Framework/ConfigParamSpec.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "MFTWorkflow/CATrackerSpec.h"
#include "ITSMFTCAWriter/MFTCATrackWriterSpec.h"

using namespace o2::framework;

void customize(std::vector<o2::framework::CallbacksPolicy>& policies)
{
  o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies);
}

void customize(std::vector<o2::framework::CompletionPolicy>& policies)
{
  policies.push_back(CompletionPolicyHelpers::consumeWhenAllOrdered(".*(?:MFT|mft).*[W,w]riter.*"));
}

void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{"disable-mc", VariantType::Bool, false, {"disable MC labels"}});
  workflowOptions.push_back(ConfigParamSpec{"disable-root-output", VariantType::Bool, false, {"do not write output root files"}});
  workflowOptions.push_back(ConfigParamSpec{"use-geom", VariantType::Bool, false, {"use geometry from the global geometry manager"}});
  workflowOptions.push_back(ConfigParamSpec{"use-irframes", VariantType::Bool, false, {"consume ITS IR frames"}});
  workflowOptions.push_back(ConfigParamSpec{"tracking-mode", VariantType::String, "sync", {"sync,async,cosmics,unset,off"}});
  workflowOptions.push_back(ConfigParamSpec{"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings (e.g. MFTCATrackerParam.nIterations=4;MFTAlpideParam.roFrameLengthInBC=594)"}});
  o2::itsmft::DPLAlpideParamInitializer::addMFTConfigOption(workflowOptions);
  o2::raw::HBFUtilsInitializer::addConfigOption(workflowOptions);
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{
  std::ignore = o2::itsmft::TrackerParamConfig<o2::detectors::DetID::MFT>::Instance();
  o2::conf::ConfigurableParam::updateFromString(config.options().get<std::string>("configKeyValues"));

  const bool useMC = !config.options().get<bool>("disable-mc");
  const bool disableRootOutput = config.options().get<bool>("disable-root-output");
  const bool useGeom = config.options().get<bool>("use-geom");
  const bool useIRFrames = config.options().get<bool>("use-irframes");
  const auto trMode = o2::itsmft::TrackingMode::fromString(config.options().get<std::string>("tracking-mode"));

  WorkflowSpec specs;
  specs.emplace_back(o2::mft::getCATrackerSpec(useMC, useGeom, useIRFrames, trMode));
  if (!disableRootOutput) {
    specs.emplace_back(o2::mft::getTrackWriterSpec(useMC, true));
  }

  o2::raw::HBFUtilsInitializer hbfIni(config, specs);

  return specs;
}
