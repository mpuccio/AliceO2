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

// Test-only CUDA kernels, hipified verbatim for HIP.

#include <cuda_runtime.h>

#include "DeviceTestKernels.h"
#include "DeviceAlgorithms.cuh"
#include "DeviceTracker.cuh"
#include "FlatStage.cuh"
#include "GpuStage.cuh"
#include "RoadGlue.cuh"
#include "SelectionGlue.cuh"
#include "SliceGlue.cuh"

#include "GenericTracking/Tracking/GPU/MaterialPort.h"
#include "GenericTracking/Tracking/Trig.h"

#include "DetectorsBase/MatLayerCylSet.h"
#include "UniformMaterial.h"

#include <cstdint>
#include <stdexcept>
#include <string>

#ifndef GTRK_BACKEND_NAMESPACE
#error "GTRK_BACKEND_NAMESPACE must name this test kernel's private namespace"
#endif

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
{
namespace
{

GPUg() void trigKernel(const double* ys, const double* xs, TrigRow* rows, std::uint32_t count)
{
  for (std::uint32_t i = (blockIdx.x * blockDim.x) + threadIdx.x; i < count; i += blockDim.x * gridDim.x) {
    rows[i] = trigRow(ys[i], xs[i]);
  }
}

cudaStream_t streamOf(DeviceRuntime& runtime, DeviceStream& stream) { return static_cast<cudaStream_t>(runtime.nativeStream(stream)); }

void runTrig(DeviceRuntime& runtime, DeviceStream& stream, const double* ys, const double* xs, TrigRow* rows, std::uint32_t count)
{
  if (count == 0U) {
    return;
  }
  const auto shape = gpu::launchFor(count);
  trigKernel<<<shape.blocks, shape.threads, 0U, streamOf(runtime, stream)>>>(ys, xs, rows, count);
  checkLaunch();
}

StageStatus runSynthetic(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, synthetic::Kind kind, std::uint32_t elements,
                         std::size_t stagingCapacity, StageOutput<TrackletT>& result, SyntheticRunInfo& info)
{
  GpuStage<TrackletT> stage{runtime, arena, stream};
  DeviceStageOutput<TrackletT> output;
  if (!output.prepare(arena, elements)) {
    return StageStatus::kOutOfMemory;
  }
  stage.prepareStaging(stagingCapacity);
  const StageStatus status = stage.run(output, synthetic::Generator{.count = {.kind = kind}});
  info.predictionMissed = stage.stagingPredictionMissed();
  info.capacityAfter = stage.stagingCapacity();
  info.scratchAfter = arena.used(ArenaEnd::kScratch);
  if (status == StageStatus::kComplete) {
    download(output, result, stream);
    stream.synchronize();
  }
  return status;
}

/// The elements a listed run generates: every one that emits, and every
/// third of the rest.
struct ListedSynthetic {
  synthetic::Count count;
  GPUdi() bool operator()(std::uint32_t element) const { return count(element) != 0U || element % 3U == 0U; }
};

StageStatus runSyntheticListed(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, synthetic::Kind kind, std::uint32_t elements,
                               std::size_t stagingCapacity, StageOutput<TrackletT>& result, SyntheticRunInfo& info)
{
  GpuStage<TrackletT> stage{runtime, arena, stream};
  DeviceStageOutput<TrackletT> output;
  if (!output.prepare(arena, elements)) {
    return StageStatus::kOutOfMemory;
  }
  stage.prepareStaging(stagingCapacity);
  auto scratch = arena.scope(ArenaEnd::kScratch);
  DeviceSpan<std::uint32_t> listed;
  std::uint32_t listedCount = 0U;
  StageStatus status =
    algorithms::selectIndices(runtime, arena, stream, elements, ListedSynthetic{.count = {.kind = kind}}, listed, listedCount, ArenaEnd::kScratch);
  if (status != StageStatus::kComplete) {
    return status;
  }
  info.listed = listedCount;
  status = stage.run(output, synthetic::Generator{.count = {.kind = kind}}, listed.data(), listedCount);
  info.predictionMissed = stage.stagingPredictionMissed();
  if (status == StageStatus::kComplete) {
    download(output, result, stream);
    stream.synchronize();
  }
  return status;
}

/// The flat form's synthetic stage:
/// element e's candidates are the rows e * 1000 + k, its context the element
/// itself, which the test checks against the row, so a pair tested with
/// another element's context is caught.
struct FlatSynthetic {
  synthetic::Count counter;
  synthetic::Accept accept;
  std::uint32_t* contexts = nullptr;
  std::uint32_t sourceBegin = 0U;
  [[nodiscard]] GPUdi() std::uint32_t count(std::uint32_t element) const { return counter(element); }
  template <typename Sink>
  GPUdi() void fill(std::uint32_t element, std::uint32_t position, Sink& sink) const
  {
    contexts[position - sourceBegin] = element;
    const std::uint32_t candidates = counter(element);
    for (std::uint32_t k = 0; k < candidates; ++k) {
      sink.put((element * 1000U) + k, position);
    }
  }
  [[nodiscard]] GPUdi() std::uint32_t sourceOf(const flat::CandidatePair& pair) const { return pair.key; }
  [[nodiscard]] GPUdi() bool test(std::uint32_t, const flat::CandidatePair& pair) const
  {
    // The pair's candidate index from its element's context: a pair tested
    // with another element's context lands outside that element's count.
    const std::uint32_t element = contexts[pair.key - sourceBegin];
    const std::uint32_t k = pair.row - (element * 1000U);
    return k < counter(element) && accept(pair.row);
  }
  [[nodiscard]] GPUdi() Tracklet emit(std::uint32_t, const flat::CandidatePair& pair, std::uint32_t element) const
  {
    return synthetic::syntheticRecord(element, pair.row - (element * 1000U));
  }
  [[nodiscard]] bool beginChunk(DeviceArena& arena, std::uint32_t begin, std::uint32_t sources, std::uint32_t)
  {
    sourceBegin = begin;
    DeviceSpan<std::uint32_t> storage;
    if (!storage.allocate(arena, std::max<std::size_t>(sources, 1U), ArenaEnd::kScratch)) {
      return false;
    }
    contexts = storage.data();
    return true;
  }
  [[nodiscard]] std::size_t contextBytes(std::uint32_t sources) const { return std::size_t{sources} * sizeof(std::uint32_t); }
  [[nodiscard]] std::size_t pairBytes() const { return 0U; }
};

/// Every element: the flat run over all of them, as a list.
struct AllElements {
  GPUdi() bool operator()(std::uint32_t) const { return true; }
};

StageStatus runFlatSynthetic(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, synthetic::Kind kind, std::uint32_t divisor,
                             std::uint32_t elements, bool listedRun, std::size_t maxChunkPairs, StageOutput<TrackletT>& result, FlatRunInfo& info)
{
  DeviceStageOutput<TrackletT> output;
  if (!output.prepare(arena, elements)) {
    return StageStatus::kOutOfMemory;
  }
  auto scratch = arena.scope(ArenaEnd::kScratch);
  DeviceSpan<std::uint32_t> listed;
  std::uint32_t listedCount = 0U;
  StageStatus status = listedRun ? algorithms::selectIndices(runtime, arena, stream, elements, ListedSynthetic{.count = {.kind = kind}}, listed, listedCount,
                                                             ArenaEnd::kScratch)
                                 : algorithms::selectIndices(runtime, arena, stream, elements, AllElements{}, listed, listedCount, ArenaEnd::kScratch);
  if (status != StageStatus::kComplete) {
    return status;
  }
  info.listed = listedCount;
  status = flat::runFlat(runtime, arena, stream, listedCount > 0U ? listed.data() : nullptr, listedCount,
                         FlatSynthetic{.counter = {.kind = kind}, .accept = {.divisor = divisor}}, output, info.report,
                         flat::FlatOptions{.maxChunkPairs = maxChunkPairs});
  if (status == StageStatus::kComplete) {
    download(output, result, stream);
    stream.synchronize();
  }
  return status;
}

struct PropagateProbe {
  const PropagationProbe* probes;
  PropagationAnswer* answers;
  GPUdi() void operator()(std::uint32_t row) const
  {
    const auto& probe = probes[row];
    auto state = probe.before;
    const auto result = propagate(probe.prop, state, probe.from, probe.to, probe.direction);
    answers[row] = {.after = state, .result = result};
  }
};
StageStatus runPropagation(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, std::span<const PropagationProbe> probes,
                           std::vector<PropagationAnswer>& answers)
{
  auto scope = arena.scope();
  DeviceSpan<PropagationProbe> inputs;
  DeviceSpan<PropagationAnswer> outputs;
  if (!inputs.allocate(arena, probes.size()) || !outputs.allocate(arena, probes.size())) {
    return StageStatus::kOutOfMemory;
  }
  upload(probes.data(), inputs, probes.size(), stream);
  algorithms::forEach(runtime, stream, static_cast<std::uint32_t>(probes.size()), PropagateProbe{inputs.data(), outputs.data()});
  answers.resize(probes.size());
  download(outputs, answers.data(), probes.size(), stream);
  stream.synchronize();
  return StageStatus::kComplete;
}

std::unique_ptr<DeviceTracker<detectors::ToyBarrel8, UniformMaterial>> makeToyBarrel8Tracker(DeviceRuntime& runtime, DeviceBlock arena)
{
  return std::make_unique<DeviceTrackerImpl<detectors::ToyBarrel8, UniformMaterial>>(runtime, arena);
}

using Toy = detectors::ToyBarrel8;
constexpr std::size_t kToyLayers = static_cast<std::size_t>(Detector<Toy>::nLayers);

StageStatus runPermutation(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, const IterationParameters<Toy>& parameters,
                           const LayerRowMap<Toy>& sources, const StageOutput<TrackletT>& tracklets, std::span<const std::uint8_t> prefix,
                           std::span<const std::uint8_t> suffix, std::vector<std::uint32_t>& result)
{
  auto scope = arena.scope();
  DeviceObject<IterationParameters<Toy>> deviceParameters;
  DeviceObject<LayerRowMap<Toy>> deviceSources;
  DeviceColumns<TrackletT> deviceTracklets;
  DeviceSpan<std::uint32_t> offsets;
  glue::DevicePathBounds paths;
  const auto rows = static_cast<std::uint32_t>(prefix.size());
  if (!deviceParameters.allocate(arena) || !deviceSources.allocate(arena) || !deviceTracklets.allocate(arena, tracklets.size()) ||
      !offsets.allocate(arena, tracklets.elementCount() + 1U) || !paths.allocate(arena, rows)) {
    return StageStatus::kOutOfMemory;
  }
  uploadObject(parameters, deviceParameters, stream);
  uploadObject(sources, deviceSources, stream);
  upload<TrackletT>(tracklets.view(), tracklets.size(), deviceTracklets, stream);
  upload(tracklets.offsets().data(), offsets, tracklets.elementCount() + 1U, stream);
  upload(prefix.data(), paths.prefix, rows, stream);
  upload(suffix.data(), paths.suffix, rows, stream);
  DeviceSpan<std::uint32_t> permutation;
  std::uint32_t size = 0U;
  const auto status = glue::permutation<Toy>(runtime, arena, stream, std::as_const(deviceTracklets).view(), tracklets.size(), deviceParameters.data(),
                                             deviceSources.data(), permutation, size, offsets.data(), paths.view(parameters.minClusters));
  if (status == StageStatus::kComplete) {
    result.resize(size);
    download(permutation, result.data(), size, stream);
  }
  stream.synchronize();
  return status;
}

StageStatus runPathBounds(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, const IterationParameters<Toy>& parameters,
                          const LayerRowMap<Toy>& sources, const StageOutput<TrackletT>& tracklets, const StageOutput<TripletT>& triplets,
                          std::span<const std::uint32_t> permutation, PathBoundsResult& result)
{
  auto scope = arena.scope();
  DeviceObject<IterationParameters<Toy>> deviceParameters;
  DeviceObject<LayerRowMap<Toy>> deviceSources;
  DeviceColumns<TrackletT> deviceTracklets;
  DeviceColumns<TripletT> deviceTriplets;
  DeviceSpan<std::uint32_t> trackletOffsets;
  DeviceSpan<std::uint32_t> tripletOffsets;
  DeviceSpan<std::uint32_t> firsts;
  if (!deviceParameters.allocate(arena) || !deviceSources.allocate(arena) || !deviceTracklets.allocate(arena, tracklets.size()) ||
      !deviceTriplets.allocate(arena, triplets.size()) || !trackletOffsets.allocate(arena, tracklets.elementCount() + 1U) ||
      !tripletOffsets.allocate(arena, triplets.elementCount() + 1U) || !firsts.allocate(arena, permutation.size())) {
    return StageStatus::kOutOfMemory;
  }
  uploadObject(parameters, deviceParameters, stream);
  uploadObject(sources, deviceSources, stream);
  upload<TrackletT>(tracklets.view(), tracklets.size(), deviceTracklets, stream);
  upload<TripletT>(triplets.view(), triplets.size(), deviceTriplets, stream);
  upload(tracklets.offsets().data(), trackletOffsets, tracklets.elementCount() + 1U, stream);
  upload(triplets.offsets().data(), tripletOffsets, triplets.elementCount() + 1U, stream);
  upload(permutation.data(), firsts, permutation.size(), stream);
  glue::DevicePathBounds measurements;
  glue::DevicePathBounds accepted;
  auto status = glue::trackletPaths<Toy>(runtime, arena, stream, deviceParameters.data(), sources, std::as_const(deviceTracklets).view(),
                                         trackletOffsets.data(), measurements);
  if (status != StageStatus::kComplete) {
    return status;
  }
  DeviceSpan<std::uint32_t> pruned;
  std::uint32_t size = 0U;
  status = glue::permutation<Toy>(runtime, arena, stream, std::as_const(deviceTracklets).view(), tracklets.size(), deviceParameters.data(),
                                  deviceSources.data(), pruned, size, trackletOffsets.data(), measurements.view(parameters.minClusters));
  if (status != StageStatus::kComplete) {
    return status;
  }
  status = glue::tripletPaths<Toy>(runtime, arena, stream, deviceParameters.data(), sources, std::as_const(deviceTracklets).view(), tracklets.size(),
                                   trackletOffsets.data(), std::as_const(deviceTriplets).view(), tripletOffsets.data(), firsts.data(), permutation.size(),
                                   accepted);
  if (status == StageStatus::kComplete) {
    result.measurementPrefix.resize(sources.size());
    result.measurementSuffix.resize(sources.size());
    result.tripletPrefix.resize(tracklets.size());
    result.tripletSuffix.resize(tracklets.size());
    result.permutation.resize(size);
    download(measurements.prefix, result.measurementPrefix.data(), sources.size(), stream);
    download(measurements.suffix, result.measurementSuffix.data(), sources.size(), stream);
    download(accepted.prefix, result.tripletPrefix.data(), tracklets.size(), stream);
    download(accepted.suffix, result.tripletSuffix.data(), tracklets.size(), stream);
    download(pruned, result.permutation.data(), size, stream);
  }
  stream.synchronize();
  return status;
}

StageStatus runEndIndex(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, const StageOutput<TripletT>& triplets, std::uint32_t trackletCount,
                        std::vector<std::uint32_t>& offsets, std::vector<std::uint32_t>& indices)
{
  auto scope = arena.scope();
  DeviceColumns<TripletT> deviceTriplets;
  if (!deviceTriplets.allocate(arena, triplets.size())) {
    return StageStatus::kOutOfMemory;
  }
  upload<TripletT>(triplets.view(), triplets.size(), deviceTriplets, stream);
  DeviceSpan<std::uint32_t> deviceOffsets;
  DeviceSpan<std::uint32_t> deviceIndices;
  const auto status =
    glue::endIndex(runtime, arena, stream, std::as_const(deviceTriplets).view(), triplets.size(), trackletCount, deviceOffsets, deviceIndices);
  if (status == StageStatus::kComplete) {
    offsets.resize(std::size_t{trackletCount} + 1U);
    indices.resize(triplets.size());
    download(deviceOffsets, offsets.data(), offsets.size(), stream);
    download(deviceIndices, indices.data(), indices.size(), stream);
  }
  stream.synchronize();
  return status;
}

StageStatus runChained(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, std::span<const std::uint32_t> neighbourOffsets,
                       std::vector<std::uint32_t>& result)
{
  auto scope = arena.scope();
  DeviceSpan<std::uint32_t> deviceOffsets;
  if (neighbourOffsets.empty()) {
    throw std::invalid_argument("neighbour offsets have one entry more than the triplets");
  }
  if (!deviceOffsets.allocate(arena, neighbourOffsets.size())) {
    return StageStatus::kOutOfMemory;
  }
  upload(neighbourOffsets.data(), deviceOffsets, neighbourOffsets.size(), stream);
  DeviceSpan<std::uint32_t> chained;
  std::uint32_t size = 0U;
  const auto triplets = static_cast<std::uint32_t>(neighbourOffsets.size() - 1U);
  const auto status = glue::chained(runtime, arena, stream, deviceOffsets.data(), triplets, chained, size);
  if (status == StageStatus::kComplete) {
    result.resize(size);
    download(chained, result.data(), size, stream);
  }
  stream.synchronize();
  return status;
}

struct SelectByMask {
  const std::uint8_t* mask = nullptr;
  GPUdi() bool operator()(std::uint32_t index) const { return mask[index] != 0U; }
};

StageStatus runSelectIndices(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, std::span<const std::uint8_t> mask,
                             std::vector<std::uint32_t>& result)
{
  auto scope = arena.scope();
  DeviceSpan<std::uint8_t> deviceMask;
  if (!deviceMask.allocate(arena, mask.size())) {
    return StageStatus::kOutOfMemory;
  }
  upload(mask.data(), deviceMask, mask.size(), stream);
  DeviceSpan<std::uint32_t> selected;
  std::uint32_t count = 0U;
  const auto status =
    algorithms::selectIndices(runtime, arena, stream, static_cast<std::uint32_t>(mask.size()), SelectByMask{.mask = deviceMask.data()}, selected, count);
  if (status == StageStatus::kComplete) {
    result.resize(count);
    download(selected, result.data(), count, stream);
  }
  stream.synchronize();
  return status;
}

template <typename T>
StageStatus runScan(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, std::span<const T> input, std::vector<T>& result)
{
  auto scope = arena.scope();
  DeviceSpan<T> deviceInput;
  DeviceSpan<T> deviceOutput;
  if (!deviceInput.allocate(arena, input.size()) || !deviceOutput.allocate(arena, input.size())) {
    return StageStatus::kOutOfMemory;
  }
  upload(input.data(), deviceInput, input.size(), stream);
  const auto status = algorithms::exclusiveScan(runtime, arena, stream, std::as_const(deviceInput).data(), deviceOutput.data(), input.size());
  if (status == StageStatus::kComplete) {
    result.resize(input.size());
    download(deviceOutput, result.data(), input.size(), stream);
  }
  stream.synchronize();
  return status;
}

StageStatus runScan32(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, std::span<const std::uint32_t> input,
                      std::vector<std::uint32_t>& result)
{
  return runScan(runtime, stream, arena, input, result);
}

StageStatus runScan64(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, std::span<const std::uint64_t> input,
                      std::vector<std::uint64_t>& result)
{
  return runScan(runtime, stream, arena, input, result);
}

/// One pass of the chain: value i rises to value i-1 plus one.
struct ChainStep {
  std::uint32_t* values = nullptr;
  std::uint32_t threshold = 0U;
  std::uint32_t* flags = nullptr;
  GPUdi() void operator()(std::uint32_t element) const
  {
    const std::uint32_t index = element + 1U;
    const std::uint32_t previous = *static_cast<const volatile std::uint32_t*>(values + index - 1U);
    if (previous == 0U) {
      return; // not reached yet
    }
    const std::uint32_t raised = previous + 1U;
    if (*static_cast<const volatile std::uint32_t*>(values + index) < raised) {
      o2::gpu::CAMath::AtomicMax(values + index, raised);
      o2::gpu::CAMath::AtomicMax(flags, 1U);
    }
    if (raised > threshold) {
      o2::gpu::CAMath::AtomicMax(flags + 1, 1U);
    }
  }
};

StageStatus runSweepChain(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, std::uint32_t length, std::uint32_t threshold,
                          std::uint32_t maxChangingSweeps, std::vector<std::uint32_t>& values, std::array<std::uint32_t, 2>& flags,
                          std::uint32_t& changingSweeps)
{
  if (length == 0U) {
    throw std::invalid_argument("a chain has at least one value");
  }
  auto scope = arena.scope();
  DeviceSpan<std::uint32_t> deviceValues;
  if (!deviceValues.allocate(arena, length)) {
    return StageStatus::kOutOfMemory;
  }
  runtime.fill(deviceValues.data(), 0U, std::size_t{length} * sizeof(std::uint32_t), stream);
  const std::uint32_t one = 1U;
  upload(&one, deviceValues, 1U, stream);
  const auto swept = algorithms::sweepToFixedPoint<2U>(runtime, arena, stream, maxChangingSweeps, [&](std::uint32_t* deviceFlags) {
    algorithms::forEach(runtime, stream, length - 1U, ChainStep{.values = deviceValues.data(), .threshold = threshold, .flags = deviceFlags});
  });
  flags = swept.flags;
  changingSweeps = swept.changingSweeps;
  values.resize(length);
  download(deviceValues, values.data(), length, stream);
  stream.synchronize();
  return swept.status;
}

StageStatus runCopyColumns(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, const StageOutput<TrackletT>& source, std::size_t sourceOffset,
                           std::size_t destinationRows, std::size_t destinationOffset, std::size_t rows, std::vector<Tracklet>& destination)
{
  auto scope = arena.scope();
  DeviceColumns<TrackletT> from;
  DeviceColumns<TrackletT> to;
  if (!from.allocate(arena, source.size()) || !to.allocate(arena, destinationRows)) {
    return StageStatus::kOutOfMemory;
  }
  upload<TrackletT>(source.view(), source.size(), from, stream);
  to.view().apply([&](auto*... column) { (runtime.fill(column, 0xFFU, destinationRows * sizeof(*column), stream), ...); });
  copyColumns(from, sourceOffset, to, destinationOffset, rows, stream);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  BoundedColumns<TrackletT> host;
  host.setMemoryResource(pool.get());
  download(to, host, destinationRows, stream);
  stream.synchronize();
  destination.resize(destinationRows);
  const auto view = host.view();
  for (std::size_t row = 0; row < destinationRows; ++row) {
    destination[row] =
      Tracklet{.source = view.source[row], .target = view.target[row], .link = view.link[row], .lambda = view.lambda[row], .phi = view.phi[row]};
  }
  return StageStatus::kComplete;
}

bool runGatherChunks(DeviceRuntime& /*runtime*/, DeviceStream& stream, DeviceArena& arena, const StageOutput<TrackletT>& source,
                     std::span<const std::size_t> chunkRows, std::vector<Tracklet>& gathered)
{
  auto base = arena.scope();
  DeviceColumns<TrackletT> from;
  if (!from.allocate(arena, source.size())) {
    return false;
  }
  upload<TrackletT>(source.view(), source.size(), from, stream);
  auto chunksScope = arena.scope(ArenaEnd::kScratch);
  DeviceColumnChunks<TrackletT> chunks;
  std::size_t at = 0U;
  for (const std::size_t rows : chunkRows) {
    if (!chunks.append(arena, ArenaEnd::kScratch, from, at, rows, stream)) {
      return false;
    }
    at += rows;
  }
  DeviceColumns<TrackletT> to;
  if (!to.allocate(arena, chunks.rows())) {
    return false;
  }
  chunks.gather(to, 0U, stream);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  BoundedColumns<TrackletT> host;
  host.setMemoryResource(pool.get());
  download(to, host, chunks.rows(), stream);
  stream.synchronize();
  gathered.resize(chunks.rows());
  const auto view = host.view();
  for (std::size_t row = 0; row < gathered.size(); ++row) {
    gathered[row] = Tracklet{.source = view.source[row], .target = view.target[row], .link = view.link[row], .lambda = view.lambda[row], .phi = view.phi[row]};
  }
  return true;
}

GPUg() void readToyFlatObject(const ToyFlatObject* object, std::int32_t* values)
{
  for (std::int32_t index = (blockIdx.x * blockDim.x) + threadIdx.x; index < object->count(); index += blockDim.x * gridDim.x) {
    values[index] = object->value(index);
  }
}

bool runFlatObject(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, const ToyFlatObject& host, std::vector<std::int32_t>& readBack)
{
  auto scope = arena.scope();
  DeviceFlatObject<ToyFlatObject> object;
  if (!object.upload(arena, host, stream)) {
    return false;
  }
  const auto count = static_cast<std::size_t>(host.count());
  DeviceSpan<std::int32_t> read;
  if (!read.allocate(arena, count)) {
    return false;
  }
  readToyFlatObject<<<4U, 64U, 0U, nativeStream(runtime, stream)>>>(object.data(), read.data());
  checkLaunch();
  readBack.resize(count);
  download(read, readBack.data(), count, stream);
  stream.synchronize();
  return true;
}

struct CountAtLevel {
  TrackletConstView tracklets;
  TripletConstView triplets;
  const std::uint32_t* offsets = nullptr;
  const std::uint32_t* from = nullptr;
  OrderedRoadCounts counts;
  std::uint32_t level = 0U;

  GPUdi() void operator()(std::uint32_t element) const
  {
    if (triplets.level[element] != level) {
      return;
    }
    counts.base(element, tracklets, triplets);
    if (level > 1U) {
      std::uint64_t sum = 0U;
      for (auto edge = offsets[element]; edge < offsets[element + 1U]; ++edge) {
        if (triplets.level[from[edge]] + 1U == level) {
          counts.add(from[edge], sum);
        }
      }
      counts.finish(element, sum);
    }
  }
};

StageStatus runOrderedRoadGlue(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, const IterationParameters<Toy>& parameters,
                               const StageOutput<TrackletT>& tracklets, const StageOutput<TripletT>& triplets, const StageOutput<NeighbourT>& neighbours,
                               const std::array<MeasurementRange, kToyLayers>& core, RoadGlueResult& result)
{
  auto scope = arena.scope();
  const auto tripletCount = triplets.size();
  DeviceObject<IterationParameters<Toy>> deviceParameters;
  DeviceColumns<TrackletT> deviceTracklets;
  DeviceColumns<TripletT> deviceTriplets;
  DeviceSpan<std::uint32_t> offsets;
  DeviceColumns<NeighbourT> from;
  DeviceSpan<MeasurementRange> deviceCore;
  if (!deviceParameters.allocate(arena) || !deviceTracklets.allocate(arena, tracklets.size()) || !deviceTriplets.allocate(arena, tripletCount) ||
      !offsets.allocate(arena, neighbours.offsets().size()) || !from.allocate(arena, neighbours.size()) || !deviceCore.allocate(arena, kToyLayers)) {
    return StageStatus::kOutOfMemory;
  }
  uploadObject(parameters, deviceParameters, stream);
  upload<TrackletT>(tracklets.view(), tracklets.size(), deviceTracklets, stream);
  upload<TripletT>(triplets.view(), tripletCount, deviceTriplets, stream);
  upload(neighbours.offsets().data(), offsets, neighbours.offsets().size(), stream);
  upload<NeighbourT>(neighbours.view(), neighbours.size(), from, stream);
  upload(core.data(), deviceCore, kToyLayers, stream);

  road_glue::OrderedCounts counts;
  auto status = counts.prepare(runtime, arena, stream, tripletCount);
  if (status != StageStatus::kComplete) {
    return status;
  }
  // Deliberately use supplied final levels, not neighbour fitting or the
  // production layer boundaries. Replay models a complete staging retry.
  std::uint32_t highest = 0U;
  for (std::uint32_t element = 0U; element < tripletCount; ++element) {
    highest = std::max(highest, triplets.view().level[element]);
  }
  for (int attempt = 0; attempt < 2; ++attempt) {
    for (std::uint32_t level = 1U; level <= highest; ++level) {
      algorithms::forEach(runtime, stream, tripletCount,
                          CountAtLevel{std::as_const(deviceTracklets).view(), std::as_const(deviceTriplets).view(), offsets.data(),
                                       std::as_const(from).view().from, counts.view(road_glue::linkSourcesOf(deviceParameters.data()), deviceCore.data()), level});
    }
  }
  status = counts.summarize(runtime, stream, std::as_const(deviceTriplets).view(), tripletCount, result.chained);
  const DeviceSpan<std::uint32_t>& successorAbove = counts.successorAbove;
  const DeviceSpan<std::uint32_t>& roads = counts.roads;
  if (status != StageStatus::kComplete) {
    stream.synchronize();
    return status;
  }
  DeviceSpan<std::uint32_t> starts;
  std::uint32_t startCount = 0U;
  status = road_glue::roadStarts(runtime, arena, stream, std::as_const(deviceTracklets).view(), std::as_const(deviceTriplets).view(), tripletCount,
                                 successorAbove.data(), roads.data(), road_glue::roadStartRuleOf(deviceParameters.data(), parameters), starts, startCount);
  if (status != StageStatus::kComplete) {
    stream.synchronize();
    return status;
  }
  DeviceSpan<std::uint32_t> roadOffsets;
  std::uint32_t total = 0U;
  status = road_glue::roadOffsets(runtime, arena, stream, starts.data(), startCount, roads.data(), roadOffsets, total);
  if (status != StageStatus::kComplete) {
    stream.synchronize();
    return status;
  }
  result.successorAbove.resize(tripletCount);
  result.roads.resize(tripletCount);
  result.starts.resize(startCount);
  result.offsets.resize(std::size_t{startCount} + 1U);
  download(successorAbove, result.successorAbove.data(), tripletCount, stream);
  download(roads, result.roads.data(), tripletCount, stream);
  download(starts, result.starts.data(), startCount, stream);
  download(roadOffsets, result.offsets.data(), result.offsets.size(), stream);
  stream.synchronize();
  return StageStatus::kComplete;
}

StageStatus runSelectionGlue(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, const ToyBarrel8Records& records,
                             std::span<const std::uint8_t> mask, std::uint32_t sharedInnermost, SelectionGlueResult& result)
{
  result = {};
  auto scope = arena.scope();
  DeviceColumns<SelectionRecordOf<kToyBarrel8Clusters>::RecordT> deviceRecords;
  DeviceSpan<std::uint8_t> deviceMask;
  const auto count = static_cast<std::uint32_t>(records.size());
  if (!deviceRecords.allocate(arena, count) || !deviceMask.allocate(arena, mask.size())) {
    return StageStatus::kOutOfMemory;
  }
  upload<SelectionRecordOf<kToyBarrel8Clusters>::RecordT>(records.view(), count, deviceRecords, stream);
  upload(mask.data(), deviceMask, mask.size(), stream);
  selection_glue::Selected selected;
  const StageStatus status = selection_glue::select<kToyBarrel8Clusters>(runtime, arena, stream, std::as_const(deviceRecords).view(), count, deviceMask.data(),
                                                                         static_cast<std::uint32_t>(mask.size()), sharedInnermost, selected);
  stream.synchronize();
  result.scratchAfter = arena.used(ArenaEnd::kScratch);
  if (status != StageStatus::kComplete) {
    return status;
  }
  result.ranked.resize(selected.unique);
  result.selected.resize(selected.tracks);
  result.mask.resize(mask.size());
  if (selected.unique != 0U) {
    download(selected.ranked, result.ranked.data(), selected.unique, stream);
  }
  if (selected.tracks != 0U) {
    download(selected.selected, result.selected.data(), selected.tracks, stream);
  }
  download(deviceMask, result.mask.data(), mask.size(), stream);
  stream.synchronize();
  result.duplicates = selected.duplicates;
  result.mismatches = selected.mismatches;
  result.rounds = selected.rounds;
  result.readbacks = selected.readbacks;
  return StageStatus::kComplete;
}

GPUg() void matBudgetKernel(const o2::base::MatLayerCylSet* lut, const float* segments, std::uint32_t count, float* result)
{
  for (std::uint32_t i = (blockIdx.x * blockDim.x) + threadIdx.x; i < count; i += blockDim.x * gridDim.x) {
    const float* s = segments + (6U * i);
    const auto budget = lut->getMatBudget(s[0], s[1], s[2], s[3], s[4], s[5]);
    result[(3U * i) + 0U] = budget.meanRho;
    result[(3U * i) + 1U] = budget.meanX2X0;
    result[(3U * i) + 2U] = budget.length;
  }
}

StageStatus runMatBudgets(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, const o2::base::MatLayerCylSet& lut,
                          std::span<const std::array<float, 6>> segments, std::vector<std::array<float, 3>>& result)
{
  auto scope = arena.scope();
  DeviceFlatObject<o2::base::MatLayerCylSet> device;
  DeviceSpan<float> deviceSegments;
  DeviceSpan<float> deviceResult;
  const auto count = static_cast<std::uint32_t>(segments.size());
  if (!device.upload(arena, lut, stream) || !deviceSegments.allocate(arena, 6U * std::size_t{count}) ||
      !deviceResult.allocate(arena, 3U * std::size_t{count})) {
    return StageStatus::kOutOfMemory;
  }
  upload(segments.data()->data(), deviceSegments, 6U * std::size_t{count}, stream);
  const auto shape = launchFor(count);
  matBudgetKernel<<<shape.blocks, shape.threads, 0U, nativeStream(runtime, stream)>>>(device.data(), deviceSegments.data(), count, deviceResult.data());
  checkLaunch();
  result.resize(count);
  download(deviceResult, result.data()->data(), 3U * std::size_t{count}, stream);
  stream.synchronize();
  return StageStatus::kComplete;
}

/// `lane_group_script::run` for one item, its group's lanes in step.
struct LaneGroupBody {
  std::uint64_t* folded = nullptr;
  template <typename Group>
  GPUdi() void together(std::uint32_t item, const Group& group) const
  {
    const auto words = lane_group_script::run(group, item);
    const std::size_t at = 2U * ((std::size_t{item} * Group::kLanes) + group.lane());
    folded[at] = words[0];
    folded[at + 1U] = words[1];
  }
};

StageStatus runLaneGroups(DeviceRuntime& runtime, DeviceStream& stream, DeviceArena& arena, std::uint32_t lanes, std::uint32_t items,
                          std::vector<std::uint64_t>& result)
{
  auto scope = arena.scope();
  DeviceSpan<std::uint64_t> folded;
  const std::size_t words = 2U * std::size_t{items} * lanes;
  if (!folded.allocate(arena, words)) {
    return StageStatus::kOutOfMemory;
  }
  const LaneGroupBody body{.folded = folded.data()};
  switch (lanes) {
    case 2U:
      algorithms::forEachGroup<2U>(runtime, stream, items, body);
      break;
    case 4U:
      algorithms::forEachGroup<4U>(runtime, stream, items, body);
      break;
    case 8U:
      algorithms::forEachGroup<8U>(runtime, stream, items, body);
      break;
    case 16U:
      algorithms::forEachGroup<16U>(runtime, stream, items, body);
      break;
    case 32U:
      algorithms::forEachGroup<32U>(runtime, stream, items, body);
      break;
    case 64U:
      if constexpr (kWarpSize == 64) {
        algorithms::forEachGroup<64U>(runtime, stream, items, body);
        break;
      }
      [[fallthrough]];
    default:
      throw std::invalid_argument("runLaneGroups: no group of " + std::to_string(lanes) + " lanes in a warp of " + std::to_string(kWarpSize));
  }
  result.resize(words);
  download(folded, result.data(), words, stream);
  stream.synchronize();
  return StageStatus::kComplete;
}

} // namespace
} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE

namespace o2::gtrk
{
#ifdef GTRK_BACKEND_HIP
const DeviceTestEntry* testKernelsHIP()
{
  static const DeviceTestEntry entry{.trig = &gpu::hip::runTrig,
                                     .runSynthetic = &gpu::hip::runSynthetic,
                                     .runSyntheticListed = &gpu::hip::runSyntheticListed,
                                     .runFlatSynthetic = &gpu::hip::runFlatSynthetic,
                                     .makeToyBarrel8Tracker = &gpu::hip::makeToyBarrel8Tracker,
                                     .permutation = &gpu::hip::runPermutation,
                                     .pathBounds = &gpu::hip::runPathBounds,
                                     .endIndex = &gpu::hip::runEndIndex,
                                     .chained = &gpu::hip::runChained,
                                     .selectIndices = &gpu::hip::runSelectIndices,
                                     .scan32 = &gpu::hip::runScan32,
                                     .scan64 = &gpu::hip::runScan64,
                                     .sweepChain = &gpu::hip::runSweepChain,
                                     .copyColumns = &gpu::hip::runCopyColumns,
                                     .gatherChunks = &gpu::hip::runGatherChunks,
                                     .flatObject = &gpu::hip::runFlatObject,
                                     .matBudgets = &gpu::hip::runMatBudgets,
                                     .selectionGlue = &gpu::hip::runSelectionGlue,
                                     .propagation = &gpu::hip::runPropagation,
                                     .laneGroups = &gpu::hip::runLaneGroups,
                                     .orderedRoadGlue = &gpu::hip::runOrderedRoadGlue};
  return &entry;
}
#else
const DeviceTestEntry* testKernelsCUDA()
{
  static const DeviceTestEntry entry{.trig = &gpu::cuda::runTrig,
                                     .runSynthetic = &gpu::cuda::runSynthetic,
                                     .runSyntheticListed = &gpu::cuda::runSyntheticListed,
                                     .runFlatSynthetic = &gpu::cuda::runFlatSynthetic,
                                     .makeToyBarrel8Tracker = &gpu::cuda::makeToyBarrel8Tracker,
                                     .permutation = &gpu::cuda::runPermutation,
                                     .pathBounds = &gpu::cuda::runPathBounds,
                                     .endIndex = &gpu::cuda::runEndIndex,
                                     .chained = &gpu::cuda::runChained,
                                     .selectIndices = &gpu::cuda::runSelectIndices,
                                     .scan32 = &gpu::cuda::runScan32,
                                     .scan64 = &gpu::cuda::runScan64,
                                     .sweepChain = &gpu::cuda::runSweepChain,
                                     .copyColumns = &gpu::cuda::runCopyColumns,
                                     .gatherChunks = &gpu::cuda::runGatherChunks,
                                     .flatObject = &gpu::cuda::runFlatObject,
                                     .matBudgets = &gpu::cuda::runMatBudgets,
                                     .selectionGlue = &gpu::cuda::runSelectionGlue,
                                     .propagation = &gpu::cuda::runPropagation,
                                     .laneGroups = &gpu::cuda::runLaneGroups,
                                     .orderedRoadGlue = &gpu::cuda::runOrderedRoadGlue};
  return &entry;
}
#endif
} // namespace o2::gtrk
