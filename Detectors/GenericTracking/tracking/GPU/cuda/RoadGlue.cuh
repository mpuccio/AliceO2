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

/// \file RoadGlue.cuh
/// \brief A slice's road counts, starts and numbering on a device.
///
/// RoadCount's passes, with its functors as they are:
///
///   OrderedCounts  the counts and successor flags the ordered neighbour
///                  generation finalises (OrderedRoadCounts), their flags'
///                  fill and the overflow and chained-count read-back: the
///                  tracker's path
///   roadStarts     the triplets RoadStartRule admits, in order
///   roadOffsets    each start's first road, a 64-bit scan whose total past
///                  32 bits is kIndexOverflow, narrowed for RoadWalk
///
/// The shared machinery is DeviceAlgorithms.cuh's; what is here is the
/// road-specific predicate, the widening and the narrowing. Outputs go to the
/// base end in the caller's scope. Included by the .cu that instantiates a
/// tracker, in the backend's namespace.

#ifndef GTRK_TRACKING_GPU_ROADGLUE_CUH_
#define GTRK_TRACKING_GPU_ROADGLUE_CUH_

#include "DeviceAlgorithms.cuh"
#include "Platform.cuh"

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/MeasurementRange.h"
#include "GenericTracking/Tracking/RoadCount.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/Triplet.h"

#include "GPUCommonAlgorithm.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::road_glue
{

using Index = std::uint32_t;

namespace detail
{
namespace
{

/// RoadStartRule over device arrays.
struct IsStart {
  TrackletConstView tracklets;
  TripletConstView triplets;
  const std::uint32_t* successorAbove = nullptr;
  const std::uint32_t* roads = nullptr;
  RoadStartRule rule;
  GPUdi() bool operator()(Index triplet) const { return rule(tracklets, triplets, successorAbove, roads, triplet); }
};

/// Each start's road count in 64 bits, and a zero after the last, so the
/// scan's last entry is the total.
struct WidenStartCounts {
  const Index* starts = nullptr;
  Index startCount = 0U;
  const std::uint32_t* roads = nullptr;
  std::uint64_t* counts = nullptr;
  GPUdi() void operator()(Index start) const { counts[start] = start < startCount ? std::uint64_t{roads[starts[start]]} : std::uint64_t{0U}; }
};

struct NarrowOffsets {
  const std::uint64_t* wide = nullptr;
  Index* offsets = nullptr;
  GPUdi() void operator()(Index at) const { offsets[at] = static_cast<Index>(wide[at]); }
};

/// The chained count -- the triplets above level 1, those with predecessors --
/// without their indices: one sum per block, one atomic add each.
GPUg() void GTRK_KERNEL_BOUNDS summarizeLevelsKernel(const std::uint32_t* levels, Index count, std::uint32_t* chainedCount)
{
  const auto element = blockIdx.x * blockDim.x + threadIdx.x;
  const std::uint32_t sum = work_group_count(element < count && levels[element] > 1U);
  if (threadIdx.x == 0U) {
    atomicAdd(chainedCount, sum);
  }
}

} // namespace
} // namespace detail

/// The address of \a parameters' link-to-source-layer table, from a device
/// pointer to them, without reading device memory on the host.
template <DetectorLayout Det>
[[nodiscard]] const LayerID* linkSourcesOf(const IterationParameters<Det>* parameters)
{
  static_assert(std::is_standard_layout_v<IterationParameters<Det>>, "offsetof needs a standard-layout type");
  return reinterpret_cast<const LayerID*>(reinterpret_cast<const std::byte*>(parameters) + offsetof(IterationParameters<Det>, linkSources));
}

/// The address of \a parameters' link-to-target-layer table, from a device
/// pointer to them, without reading device memory on the host.
template <DetectorLayout Det>
[[nodiscard]] const LayerID* linkTargetsOf(const IterationParameters<Det>* parameters)
{
  static_assert(std::is_standard_layout_v<IterationParameters<Det>>, "offsetof needs a standard-layout type");
  return reinterpret_cast<const LayerID*>(reinterpret_cast<const std::byte*>(parameters) + offsetof(IterationParameters<Det>, linkTargets));
}

/// The start rule of the iteration \a parameters point to on the device,
/// \a host the same parameters on the host: its link targets read in place.
template <DetectorLayout Det>
[[nodiscard]] RoadStartRule roadStartRuleOf(const IterationParameters<Det>* parameters, const IterationParameters<Det>& host)
{
  RoadStartRule rule = roadStartRule(host);
  rule.linkTargets = linkTargetsOf(parameters);
  return rule;
}

/// The road counts and successor flags the ordered neighbour generation
/// finalises (OrderedRoadCounts.h): only the flags and the two status words
/// are zeroed first; every element writes its count.
struct OrderedCounts {
  DeviceSpan<std::uint32_t> successorAbove;
  DeviceSpan<std::uint32_t> roads;
  DeviceSpan<std::uint32_t> overflow;

  [[nodiscard]] StageStatus prepare(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, Index count)
  {
    if (!successorAbove.allocate(arena, count) || !roads.allocate(arena, count) || !overflow.allocate(arena, 2U)) {
      return StageStatus::kOutOfMemory;
    }
    runtime.fill(successorAbove.data(), 0U, std::size_t{count} * sizeof(std::uint32_t), stream);
    runtime.fill(overflow.data(), 0U, 2U * sizeof(std::uint32_t), stream);
    return StageStatus::kComplete;
  }

  [[nodiscard]] OrderedRoadCounts view(const LayerID* linkSources, const MeasurementRange* core) { return {.roads = roads.data(), .successorAbove = successorAbove.data(), .overflow = overflow.data(), .linkSources = linkSources, .core = core}; }

  /// After the neighbour generation, its staging retries included, on its
  /// stream: the chained count and the sticky overflow flag in one read-back.
  [[nodiscard]] StageStatus summarize(DeviceRuntime& runtime, DeviceStream& stream, TripletConstView triplets, Index count, Index& chainedCount)
  {
    chainedCount = 0U;
    if (count == 0U) {
      return StageStatus::kComplete;
    }
    const auto shape = launchFor(count);
    detail::summarizeLevelsKernel<<<shape.blocks, shape.threads, 0U, nativeStream(runtime, stream)>>>(triplets.level, count, overflow.data() + 1U);
    checkLaunch();
    std::array<std::uint32_t, 2U> flags{};
    runtime.copyToHost(flags.data(), overflow.data(), sizeof(flags), stream);
    stream.synchronize();
    chainedCount = flags[1];
    return flags[0] != 0U ? StageStatus::kIndexOverflow : StageStatus::kComplete;
  }
};

/// The triplets roads start at, in arena order; \a startCount read back.
[[nodiscard]] inline StageStatus roadStarts(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, TrackletConstView tracklets, TripletConstView triplets, Index tripletCount, const std::uint32_t* successorAbove, const std::uint32_t* roads, const RoadStartRule& rule, DeviceSpan<Index>& starts, Index& startCount) { return algorithms::selectIndices(runtime, arena, stream, tripletCount, detail::IsStart{.tracklets = tracklets, .triplets = triplets, .successorAbove = successorAbove, .roads = roads, .rule = rule}, starts, startCount); }

/// Per start, the index of its first road, and the total after the last
/// (`startCount + 1` offsets); \a total read back. A total past 32 bits is
/// kIndexOverflow, as on the host.
[[nodiscard]] inline StageStatus roadOffsets(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, const Index* starts, Index startCount, const std::uint32_t* roads, DeviceSpan<Index>& offsets, Index& total)
{
  total = 0U;
  const std::size_t items = std::size_t{startCount} + 1U;
  if (!offsets.allocate(arena, items)) {
    return StageStatus::kOutOfMemory;
  }
  auto scratch = arena.scope(ArenaEnd::kScratch);
  DeviceSpan<std::uint64_t> counts;
  DeviceSpan<std::uint64_t> wide;
  if (!counts.allocate(arena, items, ArenaEnd::kScratch) || !wide.allocate(arena, items, ArenaEnd::kScratch)) {
    return StageStatus::kOutOfMemory;
  }
  algorithms::forEach(runtime, stream, static_cast<Index>(items), detail::WidenStartCounts{.starts = starts, .startCount = startCount, .roads = roads, .counts = counts.data()});
  const StageStatus scanned = algorithms::exclusiveScan(runtime, arena, stream, std::as_const(counts).data(), wide.data(), items);
  if (scanned != StageStatus::kComplete) {
    return scanned;
  }
  std::uint64_t wideTotal = 0U;
  runtime.copyToHost(&wideTotal, wide.data() + startCount, sizeof(std::uint64_t), stream);
  stream.synchronize();
  if (wideTotal > std::numeric_limits<Index>::max()) {
    return StageStatus::kIndexOverflow;
  }
  algorithms::forEach(runtime, stream, static_cast<Index>(items), detail::NarrowOffsets{.wide = std::as_const(wide).data(), .offsets = offsets.data()});
  total = static_cast<Index>(wideTotal);
  return StageStatus::kComplete;
}

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::road_glue

#endif // GTRK_TRACKING_GPU_ROADGLUE_CUH_
