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

/// \file StudyWorkflow.h
/// \brief What the detectors' study workflows share (o2-gtrk-its-, -mft- and
///        -trk-study-workflow): their options and the chain with the study's sink.
///
/// The detector's devices and messages are its own (gtrk-<det>-...,
/// GTRK/<DET>...), each device pipelined on its own (`--pipeline
/// gtrk-its-measurement-producer:2`) and configured on its own
/// (`--gtrk-mft-measurement-sink "--qa-output ..."`). With
/// `--clusters-from-upstream` the clusters come from an O2 workflow before it
/// in the pipe (`o2-its-reco-workflow --disable-tracking`), without it from
/// O2's cluster readers; run it in the timeframe's directory of an O2DPG
/// production. With `--mc-vertices` the MC vertices
/// go to the sink too, for the iterations that GTRKReco<DET>.iterationVertexSource
/// seeds from them; with GTRKReco<DET>.vertexer the sink finds its own, and
/// takes one or the other.
///
/// With `--compare-tracks` the chain also runs the tracking device
/// (gtrk-<det>-tracker, workflow/TrackerSpec.h) on the same measurements and
/// seeds, and the sink compares the device's tracks and vertices with its own
/// tracker's, every timeframe; with `--tracker-digests` the sink logs the
/// device's digests instead.

#ifndef GTRK_STUDY_STUDYWORKFLOW_H_
#define GTRK_STUDY_STUDYWORKFLOW_H_

#include "GenericTracking/Study/MCVertexReaderSpec.h"
#include "GenericTracking/Study/MeasurementSinkSpec.h"
#include "GenericTracking/Workflow/ClusterConverterSpec.h"
#include "GenericTracking/Workflow/MeasurementProducerSpec.h"
#include "GenericTracking/Workflow/TrackerSpec.h"

#include "DetectorsRaw/HBFUtilsInitializer.h"
#include "Framework/ConfigContext.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/WorkflowSpec.h"

#include <stdexcept>
#include <vector>

namespace o2::gtrk::study
{

/// Every study workflow's options.
inline void addStudyOptions(std::vector<o2::framework::ConfigParamSpec>& options)
{
  using o2::framework::ConfigParamSpec;
  using o2::framework::VariantType;
  options.push_back(ConfigParamSpec{"disable-mc", VariantType::Bool, false, {"do not propagate MC labels"}});
  options.push_back(ConfigParamSpec{"configKeyValues", VariantType::String, "", {"semicolon separated key=value strings"}});
  options.push_back(
    ConfigParamSpec{"clusters-from-upstream", VariantType::Bool, false, {"the clusters from the workflow before this one in the pipe, not O2's readers"}});
  options.push_back(ConfigParamSpec{"mc-vertices", VariantType::Bool, false, {"the MC vertices as the vertex list (any detector with an MC collision context)"}});
  options.push_back(ConfigParamSpec{
    "compare-tracks", VariantType::Bool, false, {"also run the tracking device; the sink compares its tracks and vertices with its own tracker's"}});
  options.push_back(ConfigParamSpec{
    "tracker-digests", VariantType::Bool, false, {"also run the tracking device; the sink logs its tracks' and vertices' digests, as the bench prints a device's"}});
  // Use the HBFUtils instance updated by --configKeyValues. This also makes
  // the timer advance TimingInfo::firstTForbit by nHBFPerTF for every entry.
  o2::raw::HBFUtilsInitializer::addConfigOption(options, o2::raw::HBFUtilsInitializer::HBFUSrc);
}

/// Detector \a Det's study: its \a readers (unless `--clusters-from-upstream`),
/// the converter, the producer and the sink, with the tracking device and the
/// MC vertex reader when asked for.
template <DetectorLayout Det>
[[nodiscard]] o2::framework::WorkflowSpec studyChain(const o2::framework::ConfigContext& context, const o2::framework::WorkflowSpec& readers)
{
  namespace gw = o2::gtrk::workflow;
  const bool useMC = !context.options().get<bool>("disable-mc");
  const bool mcVertices = context.options().get<bool>("mc-vertices");
  const bool compareTracks = context.options().get<bool>("compare-tracks");
  const bool trackerDigests = context.options().get<bool>("tracker-digests");
  if (compareTracks && trackerDigests) {
    throw std::invalid_argument("--compare-tracks and --tracker-digests: one or the other");
  }
  const auto products = compareTracks ? TrackerProducts::kCompare : (trackerDigests ? TrackerProducts::kDigests : TrackerProducts::kNone);
  o2::framework::WorkflowSpec specs;
  if (!context.options().get<bool>("clusters-from-upstream")) {
    specs.insert(specs.end(), readers.begin(), readers.end());
  }
  specs.push_back(gw::getClusterConverterSpec<Det>(useMC));
  specs.push_back(gw::getMeasurementProducerSpec<Det>());
  specs.push_back(getMeasurementSinkSpec<Det>(useMC, mcVertices, products));
  if (products != TrackerProducts::kNone) {
    specs.push_back(gw::getTrackerSpec<Det>(mcVertices));
  }
  if (mcVertices) {
    specs.push_back(getMCVertexReaderSpec<Det>());
  }
  return specs;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_STUDYWORKFLOW_H_
