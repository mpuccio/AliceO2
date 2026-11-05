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

/// \file TrackerDevice.h
/// \brief The tracking device's task, instantiated by each detector's unit
///        (TrackerITS.cxx, TrackerMFT.cxx, TrackerTRKFT3.cxx).
///
/// Per timeframe: the conditions; the trackers built on the first timeframe
/// and when the shipped beam line changes; the input's checks; the vertices
/// (the input seeds or the own vertexer's); the tracking; the products. The
/// timeframe ends on every way out (`endTimeframeOnExit`), and every product
/// is published, empty when there is nothing: an empty input, or a timeframe
/// dropped out of memory with `dropTFUponFailure`. It computes no truth and
/// observes no stage.

#ifndef GTRK_WORKFLOW_TRACKERDEVICE_H_
#define GTRK_WORKFLOW_TRACKERDEVICE_H_

#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Workflow/RecoParam.h"
#include "GenericTracking/Workflow/TrackerSpec.h"
#include "GenericTracking/Workflow/TrackingInterface.h"

#include "Framework/ConcreteDataMatcher.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DeviceSpec.h"
#include "Framework/EndOfStreamContext.h"
#include "Framework/InitContext.h"
#include "Framework/Logger.h"
#include "Framework/ProcessingContext.h"
#include "Framework/Task.h"
#include "Framework/TimingInfo.h"
#include "Framework/TypeTraits.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk::workflow
{

template <DetectorLayout Det>
class TrackerTask : public o2::framework::Task
{
 public:
  using Record = TrackRecord<Det>;
  static_assert(o2::framework::is_messageable<Record>::value, "a track record goes as its bytes: trivially copyable, no pointer");
  static_assert(o2::framework::is_messageable<VertexSeed>::value, "a vertex seed goes as its bytes");

  using RunConditions = Conditions<Det, true>;

  TrackerTask(bool withVertices, std::shared_ptr<RunConditions> conditions)
    : mWithVertices(withVertices), mInterface(withVertices), mConditions(std::move(conditions))
  {
  }

  void init(o2::framework::InitContext& ic) final
  {
    mConditions->init(ic);
    // The configuration, once per pipeline, as O2's ITS tracking prints its own.
    if (ic.services().get<const o2::framework::DeviceSpec>().inputTimesliceId == 0) {
      RecoParamOf<Det>::type::Instance().printKeyValues(true, true);
    }
    mInterface.initDevice({.backend = ic.options().get<std::string>("gtrk-gpu-backend"),
                           .device = ic.options().get<int>("gtrk-gpu-device"),
                           .arenaMB = static_cast<std::size_t>(ic.options().get<int>("gtrk-gpu-arena-mb"))});
  }

  /// The GPU released here: the task itself is never destroyed at the exit.
  void stop() final { mInterface.releaseDevice(); }

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final { mConditions->finaliseCCDB(matcher, obj); }

  void run(o2::framework::ProcessingContext& pc) final
  {
    using o2::framework::Output;
    mConditions->checkUpdates(pc);
    const auto orbit = pc.services().get<o2::framework::TimingInfo>().firstTForbit;
    // Every product, empty until filled: published on every way out.
    auto& tracks = pc.outputs().make<std::vector<Record>>(Output{gDataOrigin, description<Det>(kTracks), 0});
    auto& vertices = pc.outputs().make<std::vector<VertexSeed>>(Output{gDataOrigin, description<Det>(kVertexSeeds), 0});
    auto& contributors = pc.outputs().make<std::vector<std::uint32_t>>(Output{gDataOrigin, description<Det>(kVertexContributors), 0});
    const auto ended = mInterface.endTimeframeOnExit();

    const auto input = readTimeframe<Det>(pc, mWithVertices);
    mInterface.configureFor(*mConditions, input.shipped);
    if (const char* invalid = TrackingInterface<Det>::checkInput(input.columns, input.offsets, static_cast<std::size_t>(mInterface.grid().bucketCount()), input.seeds)) {
      LOGP(fatal, "{}", invalid);
      return;
    }
    if (input.columns.loc0.empty()) {
      LOGP(info, "{}: the timeframe at orbit {}: no measurements, empty products", Det::name, orbit);
      return;
    }

    const TimeFrame<Det> frame{measurementView(input.columns), mInterface.grid(), input.offsets};
    const auto found = mInterface.vertices(frame, input.seeds, false);
    if (found.dropped) {
      return;
    }
    if (found.found != nullptr) {
      vertices.assign(found.found->vertices.begin(), found.found->vertices.end());
      contributors.assign(found.found->contributors.begin(), found.found->contributors.end());
    }
    const auto tracked = mInterface.track(frame, found.seeds);
    if (tracked.dropped) {
      return;
    }
    const auto& result = *tracked.result.tracks;
    tracks.reserve(result.size());
    std::size_t failed = 0U;
    for (std::size_t track = 0; track < result.size(); ++track) {
      tracks.push_back(result.read(track));
      failed += tracks.back().status == PropagationStatus::kOk ? 0U : 1U;
    }
    LOGP(info, "{}: the timeframe at orbit {}: {} tracks ({} with a failed refit), {} vertices, {:.3f} s tracking, {:.3f} s vertexing", Det::name, orbit,
         tracks.size(), failed, vertices.size(), tracked.seconds, found.seconds);
  }

  void endOfStream(o2::framework::EndOfStreamContext& /*ec*/) final
  {
    // A validation run fails on this line: no timeframe may be dropped.
    if (const auto dropped = mInterface.droppedTimeframes(); dropped != 0U) {
      LOGP(error, "{}: {} timeframes dropped out of memory (GTRKReco<DET>.dropTFUponFailure)", Det::name, dropped);
    }
  }

 private:
  bool mWithVertices = false;
  TrackingInterface<Det> mInterface;
  std::shared_ptr<RunConditions> mConditions; ///< the run's field, geometry, MeanVertex and material LUT (Conditions.h)
};

template <DetectorLayout Det>
o2::framework::DataProcessorSpec getTrackerSpec(bool withVertices)
{
  using namespace o2::framework;
  auto inputs = measurementInputs<Det>(withVertices);
  auto conditions = std::make_shared<Conditions<Det, true>>(inputs);
  return DataProcessorSpec{
    .name = deviceName<Det>("tracker"),
    .inputs = inputs,
    .outputs = Outputs{{gDataOrigin, description<Det>(kTracks), 0, Lifetime::Timeframe},
                       {gDataOrigin, description<Det>(kVertexSeeds), 0, Lifetime::Timeframe},
                       {gDataOrigin, description<Det>(kVertexContributors), 0, Lifetime::Timeframe}},
    .algorithm = AlgorithmSpec{adaptFromTask<TrackerTask<Det>>(withVertices, conditions)},
    .options = Options{{"gtrk-gpu-backend", VariantType::String, "", {"the GPU backend to track on, cuda or hip; the host if empty"}},
                       {"gtrk-gpu-device", VariantType::Int, 0, {"the GPU backend's device index"}},
                       {"gtrk-gpu-arena-mb", VariantType::Int, 0, {"the device tracker's arena in O2's GPU pool [MB]; required with a GPU backend"}}}};
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_TRACKERDEVICE_H_
