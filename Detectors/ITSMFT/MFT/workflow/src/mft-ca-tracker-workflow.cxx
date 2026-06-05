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

#include "CommonUtils/ConfigurableParam.h"
#include "Framework/ConfigParamSpec.h"
#include "MFTWorkflow/CATrackerSpec.h"

using namespace o2::framework;

void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{"disable-mc", VariantType::Bool, false, {"disable MC labels"}});
  workflowOptions.push_back(ConfigParamSpec{"use-geom", VariantType::Bool, false, {"use geometry from the global geometry manager"}});
  workflowOptions.push_back(ConfigParamSpec{"use-irframes", VariantType::Bool, false, {"consume ITS IR frames"}});
  workflowOptions.push_back(ConfigParamSpec{"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings (e.g. MFTAlpideParam.roFrameLengthInBC=594)"}});
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{
  o2::conf::ConfigurableParam::updateFromString(config.options().get<std::string>("configKeyValues"));

  const bool useMC = !config.options().get<bool>("disable-mc");
  const bool useGeom = config.options().get<bool>("use-geom");
  const bool useIRFrames = config.options().get<bool>("use-irframes");

  WorkflowSpec specs;
  specs.emplace_back(o2::mft::getCATrackerSpec(useMC, useGeom, useIRFrames));
  return specs;
}
