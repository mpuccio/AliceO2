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

/// \file GpuStage.cuh
/// \brief One stage's generate, scan and place on a device.
///
/// generate  one thread per element runs the stage's per-element function
///           with a DeviceEmitter: every emission takes a staging slot from a
///           64-bit atomic counter and records its element and its ordinal
///           within the element; the element's count goes to count[e].
/// scan      an exclusive scan of count over n + 1 entries, count[n] = 0, so
///           that offset[n] is the total.
/// place     one thread per staged record: record k of element e goes to
///           offset[e] + k.
///
/// The atomic decides where a record waits, never where it is placed, so the
/// output is the serial reference's (the tests' test/SerialReferences.h) whatever the
/// device's schedule.
///
/// Included only by the .cu sources of a backend library or its test library,
/// in the backend's namespace, so its instantiations stay private to that
/// library.

#ifndef GTRK_TRACKING_GPU_GPUSTAGE_CUH_
#define GTRK_TRACKING_GPU_GPUSTAGE_CUH_

#include "DeviceAlgorithms.cuh"
#include "Platform.cuh"

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStageOutput.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Base/Lanes.h"
#include "GenericTracking/Tracking/StageOutput.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
{

namespace stage_detail
{

/// Where a generator's records wait for placement, as a kernel sees it.
template <template <template <class> class> class S>
struct StagingState {
  typename DeviceColumns<S>::View records;
  typename DeviceStageOutput<S>::Index* element = nullptr;
  typename DeviceStageOutput<S>::Ordinal* ordinal = nullptr;
  Counter* counter = nullptr;
  std::uint32_t* ordinalOverflow = nullptr;
  Counter capacity = 0U;
};

/// Record \a ordinal of \a element into the staging: counted before any
/// check, as the serial reference does, so that past capacity the counter is
/// the exact size a retry needs.
template <template <template <class> class> class S>
GPUdi() void stageRecord(const StagingState<S>& staging, typename DeviceStageOutput<S>::Index element, typename DeviceStageOutput<S>::Index ordinal,
                         const typename DeviceStageOutput<S>::Record& record)
{
  const Counter slot = atomicAdd(staging.counter, Counter{1U});
  if (ordinal >= DeviceStageOutput<S>::kMaxPerElement) {
    *staging.ordinalOverflow = 1U; // every writer writes the same value
    return;
  }
  if (slot >= staging.capacity) {
    return;
  }
  staging.records.apply([&](auto*... column) { record.apply([&](const auto&... value) { ((column[slot] = value), ...); }); });
  staging.element[slot] = element;
  staging.ordinal[slot] = static_cast<typename DeviceStageOutput<S>::Ordinal>(ordinal);
}

/// What an element's generator emits into on a device: the counterpart of
/// the serial reference's and CpuStage's Emitter.
template <template <template <class> class> class S>
class DeviceEmitter
{
 public:
  using Index = typename DeviceStageOutput<S>::Index;
  using Ordinal = typename DeviceStageOutput<S>::Ordinal;
  using Record = typename DeviceStageOutput<S>::Record;

  GPUdi() DeviceEmitter(const StagingState<S>& staging, Index element) : mStaging(staging), mElement(element) {}

  GPUdi() void emit(const Record& record) { stageRecord(mStaging, mElement, mCount++, record); }

  [[nodiscard]] GPUdi() Index count() const { return mCount; }

 private:
  StagingState<S> mStaging;
  Index mElement = 0U;
  Index mCount = 0U;
};

/// What the lanes of one element emit into:
/// the element's ordinals come from one counter the lanes share, so they are
/// 0 to count - 1 in the order the lanes reach it. Only for a stage whose
/// records are put in an order of their own afterwards.
template <template <template <class> class> class S>
class GroupEmitter
{
 public:
  using Index = typename DeviceStageOutput<S>::Index;
  using Record = typename DeviceStageOutput<S>::Record;

  GPUdi() GroupEmitter(const StagingState<S>& staging, Index element, Index* count) : mStaging(staging), mElement(element), mCount(count) {}

  GPUdi() void emit(const Record& record) { stageRecord(mStaging, mElement, atomicAdd(mCount, Index{1U}), record); }

 private:
  StagingState<S> mStaging;
  Index mElement = 0U;
  Index* mCount = nullptr;
};

/// One element generated by one thread: its records staged, its count written.
template <template <template <class> class> class S, typename Generate>
GPUdi() void generateOne(const Generate& generate, const StagingState<S>& staging, typename DeviceStageOutput<S>::Index element,
                         typename DeviceStageOutput<S>::Index* counts)
{
  DeviceEmitter<S> emitter{staging, element};
  generate(element, emitter);
  counts[element] = emitter.count();
}

/// One element generated by a group of `kLanes` consecutive threads
/// (`LaneGroups`), the group's count in \a emitted: every thread of the
/// block reaches both barriers, a lane without an element (\a valid false)
/// too.
template <template <template <class> class> class S, typename Generate>
GPUdi() void generateGrouped(const Generate& generate, const StagingState<S>& staging, typename DeviceStageOutput<S>::Index element, bool valid,
                             typename DeviceStageOutput<S>::Index* counts, typename DeviceStageOutput<S>::Index* emitted)
{
  constexpr std::uint32_t kLanes = Generate::kLanes;
  const std::uint32_t lane = threadIdx.x % kLanes;
  if (lane == 0U) {
    *emitted = 0U;
  }
  __syncthreads();
  if (valid) {
    GroupEmitter<S> emitter{staging, element, emitted};
    generate.body(element, emitter, Lanes{.lane = lane, .count = kLanes});
  }
  __syncthreads();
  if (lane == 0U && valid) {
    counts[element] = *emitted;
  }
}

template <template <template <class> class> class S, typename Generate>
GPUg() void GTRK_KERNEL_BOUNDS
  generateKernel(Generate generate, StagingState<S> staging, typename DeviceStageOutput<S>::Index elements, typename DeviceStageOutput<S>::Index* counts)
{
  using Index = typename DeviceStageOutput<S>::Index;
  const Index element = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (element < elements) {
    generateOne(generate, staging, element, counts);
  }
}

/// `generateKernel` over the \a listedCount elements \a listed names only,
/// one thread each; every other element's
/// count is set to 0 before.
template <template <template <class> class> class S, typename Generate>
GPUg() void GTRK_KERNEL_BOUNDS generateListedKernel(Generate generate, StagingState<S> staging, const typename DeviceStageOutput<S>::Index* listed,
                                                    typename DeviceStageOutput<S>::Index listedCount, typename DeviceStageOutput<S>::Index* counts)
{
  using Index = typename DeviceStageOutput<S>::Index;
  const Index i = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (i >= listedCount) {
    return;
  }
  generateOne(generate, staging, listed[i], counts);
}

/// \a generate's `kLanes` lanes per listed element, `body(element, emitter,
/// lanes)` (`LaneGroups`): one element per group of consecutive threads, its
/// count in shared memory. The launch covers every listed element, and every
/// thread of the block reaches both barriers.
template <template <template <class> class> class S, typename Generate>
GPUg() void GTRK_KERNEL_BOUNDS generateGroupListedKernel(Generate generate, StagingState<S> staging, const typename DeviceStageOutput<S>::Index* listed,
                                                         typename DeviceStageOutput<S>::Index listedCount, typename DeviceStageOutput<S>::Index* counts)
{
  using Index = typename DeviceStageOutput<S>::Index;
  constexpr std::uint32_t kLanes = Generate::kLanes;
  __shared__ Index emitted[static_cast<std::uint32_t>(kBlockSize) / kLanes];
  const std::uint32_t group = threadIdx.x / kLanes;
  const std::uint64_t i = ((std::uint64_t{blockIdx.x} * blockDim.x) + threadIdx.x) / kLanes;
  generateGrouped(generate, staging, i < listedCount ? listed[i] : Index{0U}, i < listedCount, counts, emitted + group);
}

template <template <template <class> class> class S>
GPUg() void GTRK_KERNEL_BOUNDS placeKernel(typename DeviceColumns<S>::ConstView staged, const typename DeviceStageOutput<S>::Index* element,
                                           const typename DeviceStageOutput<S>::Ordinal* ordinal, const typename DeviceStageOutput<S>::Index* offsets,
                                           typename DeviceColumns<S>::View placed, typename DeviceStageOutput<S>::Index records)
{
  using Index = typename DeviceStageOutput<S>::Index;
  const Index slot = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (slot < records) {
    const Index to = offsets[element[slot]] + static_cast<Index>(ordinal[slot]);
    staged.apply([&](const auto*... from) { placed.apply([&](auto*... into) { ((into[to] = from[slot]), ...); }); });
  }
}

// Layer-dependent generation appends to one staging arena across launches;
// one count readback, scan and canonical placement follow all layers.
template <template <template <class> class> class S, typename Generate>
GPUg() void GTRK_KERNEL_BOUNDS generateLayerKernel(Generate generate, StagingState<S> staging, const std::uint32_t* boundaries, std::uint32_t layer,
                                                   typename DeviceStageOutput<S>::Index* counts)
{
  using Index = typename DeviceStageOutput<S>::Index;
  const Index begin = boundaries[layer];
  const Index end = boundaries[layer + 1U];
  const Index element = begin + (blockIdx.x * blockDim.x) + threadIdx.x;
  if (element < end) {
    generateOne(generate, staging, element, counts);
  }
}

/// A generator run layer by layer: its `body` over the elements between each
/// of `layerOrder`'s `layerBoundaries` (`Tracker::LayeredNeighbours`).
template <typename Generate>
concept LayeredGenerate = requires(const Generate& generate) {
  generate.layerBoundaries;
  generate.layerOrder;
  generate.body;
};

} // namespace stage_detail

/// A generator whose every element `L` lanes walk together:
/// `body(element, emitter, lanes)`, with the
/// element's lanes emitting through one `GroupEmitter`. For a stage whose
/// records of one element are put in an order of their own afterwards.
template <typename Body, std::uint32_t L>
struct LaneGroups {
  static_assert(static_cast<std::uint32_t>(kWarpSize) % L == 0U && static_cast<std::uint32_t>(kBlockSize) % L == 0U, "an element's lanes are part of one warp");
  static constexpr std::uint32_t kLanes = L;
  Body body;
};

/// The device backend's counterpart of CpuStage<S>: the same run and
/// statuses, on one stream and one arena.
template <template <template <class> class> class S>
class GpuStage
{
 public:
  using Output = DeviceStageOutput<S>;
  using Index = typename Output::Index;
  using Ordinal = typename Output::Ordinal;

  GpuStage(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream) : mRuntime(&runtime), mArena(&arena), mStream(&stream) {}

  /// The staging the next run starts with: the estimator's prediction.
  void prepareStaging(std::size_t records) { mCapacity = records; }
  [[nodiscard]] std::size_t stagingCapacity() const noexcept { return mCapacity; }
  /// Whether the last run's staging overflowed and was retried at the exact size.
  [[nodiscard]] bool stagingPredictionMissed() const noexcept { return mMissed; }

  /// Run \a generate(element, emitter) for every element of \a output and place
  /// the result canonically.
  template <typename Generate>
  StageStatus run(Output& output, const Generate& generate)
  {
    return runOver<false>(output, generate, nullptr, 0U);
  }

  /// `run` with only the \a listedCount elements \a listed names generated,
  /// every other element's count 0: the
  /// same records, placed the same, where the others emit nothing. \a listed
  /// is device memory that outlives the call.
  template <typename Generate>
  StageStatus run(Output& output, const Generate& generate, const Index* listed, Index listedCount)
  {
    return runOver<true>(output, generate, listed, listedCount);
  }

 private:
  /// `run` over every element, or over the \a Listed ones: chosen at compile
  /// time, so that a stage never run over a list has no listed kernel and one
  /// run over a list no unlisted one.
  template <bool Listed, typename Generate>
  StageStatus runOver(Output& output, const Generate& generate, const Index* listed, Index listedCount)
  {
    mMissed = false;
    mStaged = 0U;
    const Index elements = output.elementCount();
    if (elements >= static_cast<Index>(std::numeric_limits<int>::max())) {
      return StageStatus::kIndexOverflow; // the scan's item count is an int on CUB and hipCUB
    }
    if (elements == 0U) {
      if (!output.allocateOutput(*mArena, 0U)) {
        return StageStatus::kOutOfMemory;
      }
      mRuntime->fill(output.offsets().data(), 0U, sizeof(Index), *mStream);
      output.setSize(0U);
      return StageStatus::kComplete;
    }
    for (int attempt = 0; attempt < 2; ++attempt) {
      const StageStatus status = attemptRun<Listed>(output, generate, listed, listedCount);
      if (status != StageStatus::kStagingOverflow) {
        return status;
      }
      // Rung 2 of the overflow ladder: grow once to the exact size and keep it.
      mMissed = true;
      mCapacity = static_cast<std::size_t>(mStaged);
    }
    return StageStatus::kStagingOverflow;
  }

  /// One generate, and on success the scan and the placement. Staging and scan
  /// storage are on the scratch end and released on return; the work already
  /// queued on the stream keeps using them, and the next allocation there is
  /// only used by work queued after it on the same stream.
  template <bool Listed, typename Generate>
  StageStatus attemptRun(Output& output, const Generate& generate, const Index* listed, Index listedCount)
  {
    const Index elements = output.elementCount();
    auto scratch = mArena->scope(ArenaEnd::kScratch);
    DeviceSpan<Counter> counter;
    DeviceSpan<std::uint32_t> ordinalOverflow;
    DeviceColumns<S> staged;
    DeviceSpan<Index> element;
    DeviceSpan<Ordinal> ordinal;
    if (!counter.allocate(*mArena, 1U, ArenaEnd::kScratch) || !ordinalOverflow.allocate(*mArena, 1U, ArenaEnd::kScratch) ||
        !staged.allocate(*mArena, mCapacity, ArenaEnd::kScratch) || !element.allocate(*mArena, mCapacity, ArenaEnd::kScratch) ||
        !ordinal.allocate(*mArena, mCapacity, ArenaEnd::kScratch)) {
      return StageStatus::kOutOfMemory;
    }
    auto* counts = output.counts().data();
    mRuntime->fill(counter.data(), 0U, sizeof(Counter), *mStream);
    mRuntime->fill(ordinalOverflow.data(), 0U, sizeof(std::uint32_t), *mStream);
    mRuntime->fill(counts + elements, 0U, sizeof(Index), *mStream);

    const stage_detail::StagingState<S> state{.records = staged.view(),
                                              .element = element.data(),
                                              .ordinal = ordinal.data(),
                                              .counter = counter.data(),
                                              .ordinalOverflow = ordinalOverflow.data(),
                                              .capacity = static_cast<Counter>(mCapacity)};
    const NativeStream stream = nativeStream(*mRuntime, *mStream);
    constexpr bool kLayered = stage_detail::LayeredGenerate<Generate>;
    constexpr bool kGrouped = requires { Generate::kLanes; };
    if constexpr (Listed) {
      static_assert(!kLayered, "a layered generate runs over every element");
      mRuntime->fill(counts, 0U, std::size_t{elements} * sizeof(Index), *mStream);
      if (listedCount != 0U) {
        if constexpr (kGrouped) {
          const std::uint64_t threads = std::uint64_t{listedCount} * Generate::kLanes;
          if (threads > std::numeric_limits<std::uint32_t>::max()) {
            return StageStatus::kIndexOverflow;
          }
          const auto generateShape = launchFor(static_cast<std::uint32_t>(threads));
          stage_detail::generateGroupListedKernel<S>
            <<<generateShape.blocks, generateShape.threads, 0U, stream>>>(generate, state, listed, listedCount, counts);
        } else {
          const auto generateShape = launchFor(listedCount);
          stage_detail::generateListedKernel<S><<<generateShape.blocks, generateShape.threads, 0U, stream>>>(generate, state, listed, listedCount, counts);
        }
        checkLaunch();
      }
    } else if constexpr (kLayered) {
      // One thread per possible element (launchFor): capped grids make lanes
      // process later elements alone after divergence.
      const auto generateShape = launchFor(elements);
      for (const auto layer : generate.layerOrder) {
        stage_detail::generateLayerKernel<S>
          <<<generateShape.blocks, generateShape.threads, 0U, stream>>>(generate.body, state, generate.layerBoundaries, static_cast<Index>(layer), counts);
        checkLaunch();
      }
    } else {
      static_assert(!kGrouped, "lane groups run over a list");
      const auto generateShape = launchFor(elements);
      stage_detail::generateKernel<S><<<generateShape.blocks, generateShape.threads, 0U, stream>>>(generate, state, elements, counts);
      checkLaunch();
    }

    // The one read-back of a stage: the counter is the total, since every
    // emission counted in both its element and here.
    Counter emitted = 0U;
    std::uint32_t overflowed = 0U;
    download(counter, &emitted, 1U, *mStream);
    download(ordinalOverflow, &overflowed, 1U, *mStream);
    mStream->synchronize();
    mStaged = emitted;

    // Failures in the serial reference's order.
    if (overflowed != 0U) {
      return StageStatus::kOrdinalOverflow;
    }
    if (emitted > std::numeric_limits<Index>::max()) {
      return StageStatus::kIndexOverflow;
    }
    if (emitted > mCapacity) {
      return StageStatus::kStagingOverflow;
    }
    // The scan before the output: a refusal of either leaves the base end as
    // the run found it, counts and offsets only. run() keeps elements below
    // INT_MAX, the scan's item count.
    const auto total = static_cast<Index>(emitted);
    const StageStatus scanned =
      algorithms::exclusiveScan(*mRuntime, *mArena, *mStream, std::as_const(output.counts()).data(), output.offsets().data(), std::size_t{elements} + 1U);
    if (scanned != StageStatus::kComplete || !output.allocateOutput(*mArena, total)) {
      return StageStatus::kOutOfMemory;
    }

    if (total != 0U) {
      const auto placeShape = launchFor(total);
      stage_detail::placeKernel<S><<<placeShape.blocks, placeShape.threads, 0U, stream>>>(std::as_const(staged).view(), element.data(), ordinal.data(),
                                                                                          output.offsets().data(), output.view(), total);
      checkLaunch();
    }
    output.setSize(total);
    return StageStatus::kComplete;
  }

  DeviceRuntime* mRuntime = nullptr;
  DeviceArena* mArena = nullptr;
  DeviceStream* mStream = nullptr;
  std::size_t mCapacity = 1024U;
  std::uint64_t mStaged = 0U;
  bool mMissed = false;
};

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE

#endif // GTRK_TRACKING_GPU_GPUSTAGE_CUH_
