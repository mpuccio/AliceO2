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

#include <stdexcept>
#include <vector>

namespace o2::gtrk::workflow
{

/// Every reco workflow's options; with \a o2Format the adaptor's too
/// (`--disable-adaptor`, `--disable-root-output`), for a detector with O2 track formats,
/// and the two halves of a tracking inside o2-gpu-reco-workflow: `--tracking-in-gpu-reco`
/// before it in the pipe, `--adaptor-only` after it.
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
    options.push_back(ConfigParamSpec{"tracking-in-gpu-reco", VariantType::Bool, false,
                                      {"the tracking runs inside o2-gpu-reco-workflow (--output-type gtrk-<det>-tracks), on the TPC's GPU: only the converter and "
                                       "the producer, before it in the pipe; the adaptor after it (--adaptor-only)"}});
    options.push_back(ConfigParamSpec{"adaptor-only", VariantType::Bool, false,
                                      {"only the adaptor (and O2's writers), after o2-gpu-reco-workflow in the pipe: its tracks, the measurements and clusters from upstream"}});
  }
  o2::raw::HBFUtilsInitializer::addConfigOption(options, o2::raw::HBFUtilsInitializer::HBFUSrc);
}

/// Whether the tracking runs inside o2-gpu-reco-workflow, this workflow
/// before it (`--tracking-in-gpu-reco`), or this is the adaptor after it
/// (`--adaptor-only`). False for a detector without the options.
[[nodiscard]] inline bool trackingInGPUReco(const o2::framework::ConfigContext& context)
{
  return context.options().hasOption("tracking-in-gpu-reco") && context.options().get<bool>("tracking-in-gpu-reco");
}
[[nodiscard]] inline bool adaptorOnly(const o2::framework::ConfigContext& context)
{
  return context.options().hasOption("adaptor-only") && context.options().get<bool>("adaptor-only");
}

/// Detector \a Det's tracking: its \a readers (unless `--clusters-from-upstream`),
/// the converter into our clusters, the measurement producer and the tracker
/// (unless `--tracking-in-gpu-reco`); nothing with `--adaptor-only`.
template <DetectorLayout Det>
[[nodiscard]] o2::framework::WorkflowSpec trackingChain(const o2::framework::ConfigContext& context, const o2::framework::WorkflowSpec& readers)
{
  const bool useMC = !context.options().get<bool>("disable-mc");
  o2::framework::WorkflowSpec specs;
  if (adaptorOnly(context)) {
    if (trackingInGPUReco(context)) {
      throw std::invalid_argument("--tracking-in-gpu-reco and --adaptor-only are the two halves of one tracking: one per workflow");
    }
    return specs;
  }
  if (!context.options().get<bool>("clusters-from-upstream")) {
    specs.insert(specs.end(), readers.begin(), readers.end());
  }
  specs.push_back(getClusterConverterSpec<Det>(useMC));
  specs.push_back(getMeasurementProducerSpec<Det>());
  const bool inGPUReco = trackingInGPUReco(context);
  if (inGPUReco && context.options().get<bool>("vertices-from-upstream")) {
    throw std::invalid_argument("--tracking-in-gpu-reco takes no native vertex seeds (--vertices-from-upstream)");
  }
  if (!inGPUReco) {
    specs.push_back(getTrackerSpec<Det>(context.options().get<bool>("vertices-from-upstream")));
  }
  return specs;
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_RECOWORKFLOW_H_
