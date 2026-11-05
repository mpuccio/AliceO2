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

/// \file RecoWorkflow.h
/// \brief What the detectors' reco workflows share (o2-gtrk-its-, -mft- and
///        -trk-reco-workflow): their options and the tracking chain.

#ifndef GTRK_WORKFLOW_RECOWORKFLOW_H_
#define GTRK_WORKFLOW_RECOWORKFLOW_H_

#include "GenericTracking/Workflow/ClusterConverterSpec.h"
#include "GenericTracking/Workflow/MeasurementProducerSpec.h"
#include "GenericTracking/Workflow/TrackerSpec.h"

#include "DetectorsRaw/HBFUtilsInitializer.h"
#include "Framework/ConfigContext.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/WorkflowSpec.h"

#include <vector>

namespace o2::gtrk::workflow
{

/// Every reco workflow's options; with \a o2Format the adaptor's too
/// (`--disable-adaptor`, `--disable-root-output`), for a detector with O2 track formats.
inline void addRecoOptions(std::vector<o2::framework::ConfigParamSpec>& options, bool o2Format)
{
  using o2::framework::ConfigParamSpec;
  using o2::framework::VariantType;
  options.push_back(ConfigParamSpec{"disable-mc", VariantType::Bool, false, {"do not propagate MC labels"}});
  options.push_back(ConfigParamSpec{"configKeyValues", VariantType::String, "", {"semicolon separated key=value strings"}});
  options.push_back(
    ConfigParamSpec{"clusters-from-upstream", VariantType::Bool, false, {"the clusters from the workflow before this one in the pipe, not O2's readers"}});
  options.push_back(
    ConfigParamSpec{"vertices-from-upstream", VariantType::Bool, false, {"the native vertex seeds (GTRK/<DET>VERTICES) from the workflow before this one in the pipe"}});
  if (o2Format) {
    options.push_back(ConfigParamSpec{"disable-adaptor", VariantType::Bool, false, {"only the native tracks and vertices (GTRK/<DET>...), not O2's formats nor their files"}});
    options.push_back(ConfigParamSpec{"disable-root-output", VariantType::Bool, false, {"O2's formats without O2's writers"}});
  }
  o2::raw::HBFUtilsInitializer::addConfigOption(options, o2::raw::HBFUtilsInitializer::HBFUSrc);
}

/// Detector \a Det's tracking: its \a readers (unless `--clusters-from-upstream`),
/// the converter into our clusters, the measurement producer and the tracker.
template <DetectorLayout Det>
[[nodiscard]] o2::framework::WorkflowSpec trackingChain(const o2::framework::ConfigContext& context, const o2::framework::WorkflowSpec& readers)
{
  const bool useMC = !context.options().get<bool>("disable-mc");
  o2::framework::WorkflowSpec specs;
  if (!context.options().get<bool>("clusters-from-upstream")) {
    specs.insert(specs.end(), readers.begin(), readers.end());
  }
  specs.push_back(getClusterConverterSpec<Det>(useMC));
  specs.push_back(getMeasurementProducerSpec<Det>());
  specs.push_back(getTrackerSpec<Det>(context.options().get<bool>("vertices-from-upstream")));
  return specs;
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_RECOWORKFLOW_H_
