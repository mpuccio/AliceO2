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
/// dropped out of memory with `dropTFUponFailure` (its vertices too). It computes no truth and
/// observes no stage.

#ifndef GTRK_WORKFLOW_TRACKERDEVICE_H_
#define GTRK_WORKFLOW_TRACKERDEVICE_H_

#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Workflow/Conditions.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/GPURecoTracking.h"
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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace o2::gtrk::workflow
{

/// The tracking device's products made in \a pc's outputs, empty, published on
/// every way out: (tracks, vertices, contributors), DPL's vectors.
template <DetectorLayout Det>
[[nodiscard]] auto makeTrackerProducts(o2::framework::ProcessingContext& pc)
{
  using o2::framework::Output;
  static_assert(o2::framework::is_messageable<TrackRecord<Det>>::value, "a track record goes as its bytes: trivially copyable, no pointer");
  static_assert(o2::framework::is_messageable<VertexSeed>::value, "a vertex seed goes as its bytes");
  auto& tracks = pc.outputs().make<std::vector<TrackRecord<Det>>>(Output{gDataOrigin, description<Det>(kTracks), 0});
  auto& vertices = pc.outputs().make<std::vector<VertexSeed>>(Output{gDataOrigin, description<Det>(kVertexSeeds), 0});
  auto& contributors = pc.outputs().make<std::vector<std::uint32_t>>(Output{gDataOrigin, description<Det>(kVertexContributors), 0});
  return std::tie(tracks, vertices, contributors);
}

/// The products' output specs.
template <DetectorLayout Det>
[[nodiscard]] std::vector<o2::framework::OutputSpec> trackerOutputs()
{
  using o2::framework::Lifetime;
  return {{gDataOrigin, description<Det>(kTracks), 0, Lifetime::Timeframe},
          {gDataOrigin, description<Det>(kVertexSeeds), 0, Lifetime::Timeframe},
          {gDataOrigin, description<Det>(kVertexContributors), 0, Lifetime::Timeframe}};
}

/// One timeframe of \a Det's tracking, on \a interface's host or device
/// (its block set): the conditions, the input's checks, the vertices, the
/// tracks, the products. The tracking device's run, and the GPU
/// reconstruction's (GPURecoTracking.h).
template <DetectorLayout Det, class RunConditions>
void trackTimeframe(o2::framework::ProcessingContext& pc, TrackingInterface<Det>& interface, RunConditions& conditions, bool withVertices)
{
  conditions.checkUpdates(pc);
  const auto orbit = pc.services().get<o2::framework::TimingInfo>().firstTForbit;
  auto [tracks, vertices, contributors] = makeTrackerProducts<Det>(pc);
  const auto ended = interface.endTimeframeOnExit();

  const auto input = readTimeframe<Det>(pc, withVertices);
  interface.configureFor(conditions, input.shipped);
  if (const char* invalid = TrackingInterface<Det>::checkInput(input.columns, input.offsets, static_cast<std::size_t>(interface.grid().bucketCount()), input.seeds)) {
    LOGP(fatal, "{}", invalid);
    return;
  }
  if (input.columns.loc0.empty()) {
    LOGP(info, "{}: the timeframe at orbit {}: no measurements, empty products", Det::name, orbit);
    return;
  }

  const TimeFrame<Det> frame{measurementView(input.columns), interface.grid(), input.offsets};
  const auto found = interface.vertices(frame, input.seeds, false);
  if (found.dropped) {
    return;
  }
  const auto tracked = interface.track(frame, found.seeds);
  if (tracked.dropped) {
    return;
  }
  // The vertices with the tracks only: a dropped timeframe publishes nothing.
  if (found.found != nullptr) {
    vertices.assign(found.found->vertices.begin(), found.found->vertices.end());
    contributors.assign(found.found->contributors.begin(), found.found->contributors.end());
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

/// At the end of the stream: the timeframes dropped, an error line.
template <DetectorLayout Det>
void reportDropped(const TrackingInterface<Det>& interface)
{
  // A validation run fails on this line: no timeframe may be dropped.
  if (const auto dropped = interface.droppedTimeframes(); dropped != 0U) {
    LOGP(error, "{}: {} timeframes dropped out of memory (GTRKReco<DET>.dropTFUponFailure)", Det::name, dropped);
  }
}

/// The configuration, once per pipeline, as O2's ITS tracking prints its own.
template <DetectorLayout Det>
void printConfiguration(o2::framework::InitContext& ic)
{
  if (ic.services().get<const o2::framework::DeviceSpec>().inputTimesliceId == 0) {
    RecoParamOf<Det>::type::Instance().printKeyValues(true, true);
  }
}

template <DetectorLayout Det>
class TrackerTask : public o2::framework::Task
{
 public:
  using RunConditions = Conditions<Det, true>;

  TrackerTask(bool withVertices, std::shared_ptr<RunConditions> conditions)
    : mWithVertices(withVertices), mInterface(withVertices), mConditions(std::move(conditions))
  {
  }

  void init(o2::framework::InitContext& ic) final
  {
    mConditions->init(ic);
    printConfiguration<Det>(ic);
    mInterface.initDevice({.backend = ic.options().get<std::string>("gtrk-gpu-backend"),
                           .device = ic.options().get<int>("gtrk-gpu-device"),
                           .arenaMB = static_cast<std::size_t>(ic.options().get<int>("gtrk-gpu-arena-mb"))});
  }

  /// The GPU released here: the task itself is never destroyed at the exit.
  void stop() final { mInterface.releaseDevice(); }

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final { mConditions->finaliseCCDB(matcher, obj); }

  void run(o2::framework::ProcessingContext& pc) final { trackTimeframe<Det>(pc, mInterface, *mConditions, mWithVertices); }

  void endOfStream(o2::framework::EndOfStreamContext& /*ec*/) final { reportDropped<Det>(mInterface); }

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
    .outputs = trackerOutputs<Det>(),
    .algorithm = AlgorithmSpec{adaptFromTask<TrackerTask<Det>>(withVertices, conditions)},
    .options = Options{{"gtrk-gpu-backend", VariantType::String, "", {"the GPU backend to track on, cuda or hip; the host if empty"}},
                       {"gtrk-gpu-device", VariantType::Int, 0, {"the GPU backend's device index"}},
                       {"gtrk-gpu-arena-mb", VariantType::Int, 0, {"the device tracker's arena in O2's GPU pool [MB]; required with a GPU backend"}}}};
}

/// \a Det's tracking inside the GPU reconstruction (GPURecoTracking.h): the
/// tracking device's timeframe on the reconstruction's GPU, in the pool's
/// free memory after the TPC, borrowed per timeframe.
template <DetectorLayout Det>
class GPURecoTrackingOf final : public GPURecoTracking
{
 public:
  using RunConditions = Conditions<Det, true, ConditionsHolder::kGPUReco>;

  /// No native seeds: the vertex list, where an iteration needs one, is the
  /// own vertexer's.
  GPURecoTrackingOf() : mInterface(false) {}

  void appendInputs(std::vector<o2::framework::InputSpec>& inputs) final
  {
    auto own = measurementInputs<Det>(false);
    mConditions = std::make_unique<RunConditions>(own);
    inputs.insert(inputs.end(), own.begin(), own.end());
  }

  void appendOutputs(std::vector<o2::framework::OutputSpec>& outputs) final
  {
    const auto own = trackerOutputs<Det>();
    outputs.insert(outputs.end(), own.begin(), own.end());
  }

  void init(o2::framework::InitContext& ic, o2::gpu::GPUReconstruction& reconstruction) final
  {
    mConditions->init(ic);
    printConfiguration<Det>(ic);
    mInterface.initDevice(reconstruction);
  }

  bool finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final { return mConditions->finaliseCCDB(matcher, obj); }

  void run(o2::framework::ProcessingContext& pc, bool reconstructed) final
  {
    if (!reconstructed) {
      static_cast<void>(makeTrackerProducts<Det>(pc));
      return;
    }
    // Given back after the timeframe's end, which detaches the tracker.
    const auto borrowed = mInterface.context().borrow(tag());
    mInterface.useBlock(borrowed->block());
    trackTimeframe<Det>(pc, mInterface, *mConditions, false);
    LOGP(info, "{}: the GPU reconstruction's pool lent {:.1f} MB, the arena's peak {:.1f} MB", Det::name, static_cast<double>(borrowed->block().bytes) / 1.e6,
         static_cast<double>(mInterface.lastArenaPeakBytes()) / 1.e6);
  }

  void stop() final { mInterface.releaseDevice(); }

  void endOfStream() final { reportDropped<Det>(mInterface); }

 private:
  /// The pool's stack tag of the borrow, "GTRK" and the detector's origin
  /// as O2's qStr2Tag packs a tag: the first character lowest.
  [[nodiscard]] static std::uint64_t tag()
  {
    const std::string name = "GTRK" + detectorOrigin<Det>().template as<std::string>();
    std::uint64_t packed = 0U;
    for (std::size_t i = std::min<std::size_t>(name.size(), 8U); i-- > 0U;) {
      packed = (packed << 8U) | static_cast<std::uint8_t>(name[i]);
    }
    return packed;
  }

  TrackingInterface<Det> mInterface;
  std::unique_ptr<RunConditions> mConditions; ///< the detector's own conditions; the run-wide ones are the GPU reconstruction's
};

template <DetectorLayout Det>
std::unique_ptr<GPURecoTracking> makeGPURecoTrackingOf()
{
  return std::make_unique<GPURecoTrackingOf<Det>>();
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_TRACKERDEVICE_H_
