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

/// \file SliceGlue.cuh
/// \brief The passes between a slice's stages on a device.
///
/// Each gives exactly the host class's result:
///
///   permutation   TripletPermutation: the eligible first tracklets, by their
///                 middle measurement's source element or in arena order
///   endIndex      TripletEndIndex: the triplets by their second tracklet
///   chained       the triplets with predecessors, in arena order: the slice
///                 download's, which the tests compare with `chainedTriplets`
///
/// Outputs go to the base end in the caller's scope, temporaries to the
/// scratch end, released on return. A pass that needs a count reads it back
/// and synchronises; the others are only queued. The shared parts are
/// DeviceAlgorithms.cuh's. Included by the .cu that
/// instantiates a tracker, in the backend's namespace.

#ifndef GTRK_TRACKING_GPU_SLICEGLUE_CUH_
#define GTRK_TRACKING_GPU_SLICEGLUE_CUH_

#include "DeviceAlgorithms.cuh"
#include "GenericTracking/Tracking/TrackletPaths.h"
#include "Platform.cuh"

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStageOutput.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LayerRowMap.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/TrackletSort.h"
#include "GenericTracking/Tracking/Triplet.h"
#include "GenericTracking/Tracking/TripletPermutation.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::glue
{

using Index = std::uint32_t;

// Unnamed, so that two sources of one library may include this header.
namespace detail
{
namespace
{

/// Whether a tracklet's middle measurement is a source of the slice with
/// tracklets of its own, on a tracklet path of the iteration's minClusters:
/// the triplet stage's first tracklets, as TripletPermutation::build.
template <DetectorLayout Det>
struct EligibleFirst {
  TrackletConstView tracklets;
  const IterationParameters<Det>* parameters = nullptr;
  const LayerRowMap<Det>* sources = nullptr;
  const Index* offsets = nullptr;
  TrackletPathBounds paths;
  GPUdi() bool operator()(Index first) const { return firstTrackletEligible(*parameters, *sources, tracklets, offsets, paths, first); }
};

/// The end index's sort keys: each triplet's second tracklet, the triplet as
/// its value, and the tracklet's count of them.
struct EndIndexKeys {
  TripletConstView triplets;
  Index* keys = nullptr;
  Index* values = nullptr;
  Index* perTracklet = nullptr;
  GPUdi() void operator()(Index triplet) const
  {
    const Index second = triplets.second[triplet];
    keys[triplet] = second;
    values[triplet] = triplet;
    atomicAdd(perTracklet + second, 1U);
  }
};

/// Whether a triplet has predecessors: its neighbour range is not empty.
struct HasPredecessors {
  const Index* neighbourOffsets = nullptr;
  GPUdi() bool operator()(Index triplet) const { return neighbourOffsets[triplet + 1U] != neighbourOffsets[triplet]; }
};

/// The tracklet sort's first pass: each tracklet's target as its key, its
/// row as its value.
struct TrackletTargetKeys {
  const Index* target = nullptr;
  Index* keys = nullptr;
  Index* slots = nullptr;
  GPUdi() void operator()(Index i) const
  {
    keys[i] = target[i];
    slots[i] = i;
  }
};

/// Its second: the link and lambda of the row each position holds after the
/// first pass.
struct TrackletLinkLambdaKeys {
  const std::uint16_t* link = nullptr;
  const Scalar* lambda = nullptr;
  const Index* slots = nullptr;
  std::uint64_t* keys = nullptr;
  GPUdi() void operator()(Index i) const
  {
    const Index row = slots[i];
    keys[i] = trackletLinkLambdaKey(link[row], lambda[row]);
  }
};

} // namespace
} // namespace detail

/// Each source's tracklets in (link, lambda, target) order, in place:
/// every source a segment of the library's
/// segmented sort, stable, first by target, then by link and lambda, so that
/// equal links and lambdas keep the targets' order -- `trackletBefore`'s order
/// for any segment length and any ties. Then each column gathered in that
/// order through one buffer and copied back. Every buffer lives on the
/// arena's scratch end for this call; a refusal is kOutOfMemory, with nothing
/// written to the tracklets (the copies back are the last step). Queued on
/// \a stream.
[[nodiscard]] inline StageStatus sortTracklets(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, DeviceStageOutput<TrackletT>& tracklets)
{
  const Index count = tracklets.size();
  if (count == 0U) {
    return StageStatus::kComplete;
  }
  const Index segments = tracklets.elementCount();
  const Index* offsets = std::as_const(tracklets).offsets().data();
  auto scratch = arena.scope(ArenaEnd::kScratch);
  std::array<DeviceSpan<Index>, 2> slots;
  std::array<DeviceSpan<Index>, 2> targetKeys;
  std::array<DeviceSpan<std::uint64_t>, 2> linkLambdaKeys;
  DeviceSpan<Index> column;
  for (auto* span : {&slots[0], &slots[1], &targetKeys[0], &targetKeys[1], &column}) {
    if (!span->allocate(arena, count, ArenaEnd::kScratch)) {
      return StageStatus::kOutOfMemory;
    }
  }
  for (auto& span : linkLambdaKeys) {
    if (!span.allocate(arena, count, ArenaEnd::kScratch)) {
      return StageStatus::kOutOfMemory;
    }
  }
  const auto view = tracklets.view();

  algorithms::forEach(runtime, stream, count, detail::TrackletTargetKeys{.target = view.target, .keys = targetKeys[0].data(), .slots = slots[0].data()});
  cubapi::DoubleBuffer<Index> firstKeys{targetKeys[0].data(), targetKeys[1].data()};
  cubapi::DoubleBuffer<Index> values{slots[0].data(), slots[1].data()};
  StageStatus status = algorithms::segmentedStableSortPairs(runtime, arena, stream, firstKeys, values, count, segments, offsets, offsets + 1);
  if (status != StageStatus::kComplete) {
    return status;
  }

  algorithms::forEach(runtime, stream, count,
                      detail::TrackletLinkLambdaKeys{.link = view.link, .lambda = view.lambda, .slots = values.Current(), .keys = linkLambdaKeys[0].data()});
  cubapi::DoubleBuffer<std::uint64_t> secondKeys{linkLambdaKeys[0].data(), linkLambdaKeys[1].data()};
  status = algorithms::segmentedStableSortPairs(runtime, arena, stream, secondKeys, values, count, segments, offsets, offsets + 1);
  if (status != StageStatus::kComplete) {
    return status;
  }

  const Index* sorted = values.Current();
  view.apply([&](auto*... columns) {
    const auto gather = [&](auto* data) {
      using T = std::remove_pointer_t<decltype(data)>;
      static_assert(sizeof(T) <= sizeof(Index));
      auto* buffer = reinterpret_cast<T*>(column.data());
      algorithms::forEach(runtime, stream, count, algorithms::Gather<T>{.from = data, .indices = sorted, .to = buffer});
      runtime.copyOnDevice(data, buffer, std::size_t{count} * sizeof(T), stream);
    };
    (gather(columns), ...);
  });
  return StageStatus::kComplete;
}

/// The triplet stage's first tracklets, as TripletPermutation::build: of the
/// \a count tracklets, the indices of those whose middle measurement is a
/// source with tracklets (\a offsets) and that are on a path of \a paths'
/// length, in arena order. \a permutation gets \a count entries, of which the
/// first \a size are the order.
template <DetectorLayout Det>
[[nodiscard]] StageStatus permutation(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, TrackletConstView tracklets, Index count,
                                      const IterationParameters<Det>* parameters, const LayerRowMap<Det>* sources, DeviceSpan<Index>& permutation, Index& size,
                                      const Index* offsets, TrackletPathBounds paths)
{
  return algorithms::selectIndices(
    runtime, arena, stream, count,
    detail::EligibleFirst<Det>{.tracklets = tracklets, .parameters = parameters, .sources = sources, .offsets = offsets, .paths = paths}, permutation, size);
}

/// The triplets by their second tracklet: \a offsets, `trackletCount + 1`
/// of them, delimit the triplet indices in \a indices, in arena order within
/// each tracklet.
[[nodiscard]] inline StageStatus endIndex(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, TripletConstView triplets, Index tripletCount,
                                          Index trackletCount, DeviceSpan<Index>& offsets, DeviceSpan<Index>& indices)
{
  const std::size_t offsetCount = std::size_t{trackletCount} + 1U;
  if (!offsets.allocate(arena, offsetCount) || !indices.allocate(arena, tripletCount)) {
    return StageStatus::kOutOfMemory;
  }
  if (tripletCount == 0U) {
    runtime.fill(offsets.data(), 0U, offsetCount * sizeof(Index), stream);
    return StageStatus::kComplete;
  }
  auto scratch = arena.scope(ArenaEnd::kScratch);
  DeviceSpan<Index> keys;
  DeviceSpan<Index> sortedKeys;
  DeviceSpan<Index> values;
  DeviceSpan<Index> perTracklet;
  if (!keys.allocate(arena, tripletCount, ArenaEnd::kScratch) || !sortedKeys.allocate(arena, tripletCount, ArenaEnd::kScratch) ||
      !values.allocate(arena, tripletCount, ArenaEnd::kScratch) || !perTracklet.allocate(arena, offsetCount, ArenaEnd::kScratch)) {
    return StageStatus::kOutOfMemory;
  }
  runtime.fill(perTracklet.data(), 0U, offsetCount * sizeof(Index), stream);
  algorithms::forEach(runtime, stream, tripletCount, detail::EndIndexKeys{triplets, keys.data(), values.data(), perTracklet.data()});
  const StageStatus sorted = algorithms::sortPairs(runtime, arena, stream, std::as_const(keys).data(), sortedKeys.data(), std::as_const(values).data(),
                                                   indices.data(), tripletCount, algorithms::keyBits(trackletCount));
  if (sorted != StageStatus::kComplete) {
    return sorted;
  }
  // The last count is zero, so the scan's last offset is the total.
  return algorithms::exclusiveScan(runtime, arena, stream, std::as_const(perTracklet).data(), offsets.data(), offsetCount);
}

/// The triplets whose predecessor range is not empty, in arena order, from
/// the neighbour stage's `tripletCount + 1` offsets.
[[nodiscard]] inline StageStatus chained(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, const Index* neighbourOffsets, Index tripletCount,
                                         DeviceSpan<Index>& chained, Index& size)
{
  return algorithms::selectIndices(runtime, arena, stream, tripletCount, detail::HasPredecessors{.neighbourOffsets = neighbourOffsets}, chained, size);
}

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::glue

#endif // GTRK_TRACKING_GPU_SLICEGLUE_CUH_
