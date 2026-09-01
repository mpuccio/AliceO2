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

/// @file mft-ca-reco-workflow.cxx

#include "MFTWorkflow/CARecoWorkflow.h"

#include <string>

#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsITSMFT/DPLAlpideParamInitializer.h"
#include "DetectorsRaw/HBFUtilsInitializer.h"
#include "Framework/CallbacksPolicy.h"
#include "Framework/CompletionPolicyHelpers.h"
#include "ITSMFTTracking/TrackingConfigParam.h"

using namespace o2::framework;

void customize(std::vector<o2::framework::CallbacksPolicy>& policies)
{
  o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies);
}

void customize(std::vector<o2::framework::CompletionPolicy>& policies)
{
  policies.push_back(CompletionPolicyHelpers::consumeWhenAllOrdered(".*(?:MFT|mft).*[W,w]riter.*"));
}

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<o2::framework::ConfigParamSpec> options{
    {"digits-from-upstream", o2::framework::VariantType::Bool, false, {"digits will be provided from upstream, skip digits reader"}},
    {"clusters-from-upstream", o2::framework::VariantType::Bool, false, {"clusters will be provided from upstream, skip clusterizer"}},
    {"disable-root-output", o2::framework::VariantType::Bool, false, {"do not write output root files"}},
    {"disable-mc", o2::framework::VariantType::Bool, false, {"disable MC propagation even if available"}},
    {"disable-tracking", o2::framework::VariantType::Bool, false, {"disable tracking step"}},
    {"run-assessment", o2::framework::VariantType::Bool, false, {"run MFT assessment workflow"}},
    {"disable-process-gen", o2::framework::VariantType::Bool, false, {"disable processing of all generated tracks (depends on --run-assessment)"}},
    {"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings"}},
    {"nThreads", VariantType::Int, 1, {"Number of CA tracker threads"}},
    {"use-full-geometry", o2::framework::VariantType::Bool, false, {"use full geometry instead of the light-weight MFT part"}},
    {"use-irframes", o2::framework::VariantType::Bool, false, {"consume ITS IR frames"}},
    {"tracking-mode", VariantType::String, "sync", {"sync,async,cosmics,unset,off"}},
    {"run-tracks2records", o2::framework::VariantType::Bool, false, {"run MFT alignment tracks to records workflow"}},
    {"cluster-rof-branch-only", o2::framework::VariantType::Bool, false, {"writer will store only ClustersROF branch"}}};
  o2::raw::HBFUtilsInitializer::addConfigOption(options);
  o2::itsmft::DPLAlpideParamInitializer::addMFTConfigOption(options);
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& configContext)
{
  auto trackerConfig = std::string{"MFTCATrackerParam.nThreads="} +
                       std::to_string(configContext.options().get<int>("nThreads"));
  const auto configKeyValues = configContext.options().get<std::string>("configKeyValues");
  if (!configKeyValues.empty()) {
    trackerConfig += ';';
    trackerConfig += configKeyValues;
  }
  o2::conf::ConfigurableParam::updateFromString(trackerConfig);
  o2::conf::ConfigurableParam::writeINI("o2mftcarecoflow_configuration.ini");

  const bool useMC = !configContext.options().get<bool>("disable-mc");
  const bool upstreamDigits = configContext.options().get<bool>("digits-from-upstream");
  const bool upstreamClusters = configContext.options().get<bool>("clusters-from-upstream");
  const bool disableRootOutput = configContext.options().get<bool>("disable-root-output");
  const bool runAssessment = configContext.options().get<bool>("run-assessment");
  const bool processGen = !configContext.options().get<bool>("disable-process-gen");
  const bool runTracking = !configContext.options().get<bool>("disable-tracking");
  const bool runTracks2Records = configContext.options().get<bool>("run-tracks2records");
  const bool useGeom = configContext.options().get<bool>("use-full-geometry");
  const bool useIRFrames = configContext.options().get<bool>("use-irframes");
  const bool doStag = o2::itsmft::DPLAlpideParamInitializer::isMFTStaggeringEnabled(configContext);
  const bool clrofOnly = configContext.options().get<bool>("cluster-rof-branch-only");
  const auto trackingMode = o2::itsmft::TrackingMode::fromString(configContext.options().get<std::string>("tracking-mode"));

  auto workflow = o2::mft::ca_reco_workflow::getWorkflow(
    useMC,
    doStag,
    useGeom,
    useIRFrames,
    upstreamDigits,
    upstreamClusters,
    clrofOnly,
    disableRootOutput,
    runAssessment,
    processGen,
    runTracking,
    trackingMode,
    runTracks2Records);

  o2::raw::HBFUtilsInitializer hbfInitializer(configContext, workflow);
  return workflow;
}
