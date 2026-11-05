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

/// \file PathGlue.cuh
/// \brief TrackletPaths and TripletPaths on a device: layer sweeps in stream
///        order.

#ifndef GTRK_TRACKING_GPU_PATHGLUE_CUH_
#define GTRK_TRACKING_GPU_PATHGLUE_CUH_

#include "DeviceAlgorithms.cuh"
#include "GenericTracking/Tracking/GPU/DeviceStageOutput.h"
#include "GenericTracking/Tracking/TrackletPaths.h"
#include "GenericTracking/Tracking/Triplet.h"

#include "GPUCommonAlgorithm.h"

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::glue
{

/// The path bounds of a slice on the device, in the caller's scope: per
/// element its prefix and suffix, and each layer's first permuted element and
/// first kept triplet.
struct DevicePathBounds {
  DeviceSpan<std::uint8_t> prefix;
  DeviceSpan<std::uint8_t> suffix;
  DeviceSpan<std::uint32_t> layerElements; ///< per layer, the first permuted element whose tracklet starts on it or after
  DeviceSpan<std::uint32_t> layerTriplets; ///< canonical triplet ranges for ordered neighbour generation
  [[nodiscard]] TrackletPathBounds view(std::uint32_t minClusters) const { return {prefix.data(), suffix.data(), minClusters}; }
  [[nodiscard]] bool allocate(DeviceArena& arena, std::uint32_t count) { return prefix.allocate(arena, count) && suffix.allocate(arena, count); }
};

namespace path_detail
{
using Index = std::uint32_t;

struct InitializePrefix {
  Index* prefix;
  GPUdi() void operator()(Index element) const { prefix[element] = 1U; }
};

struct PackPrefix {
  const Index* prefix;
  std::uint8_t* packed;
  GPUdi() void operator()(Index element) const { packed[element] = static_cast<std::uint8_t>(prefix[element]); }
};

/// One layer of TrackletPaths::build's forward (prefix) or reverse (suffix)
/// sweep: one element per source element of \a layer, `local` within it.
template <DetectorLayout Det, bool Forward>
struct MeasurementPaths {
  TrackletConstView tracklets;
  const Index* offsets = nullptr;
  LayerRowMap<Det> sources;
  const IterationParameters<Det>* parameters = nullptr;
  LayerID layer = 0;
  Index* prefix = nullptr;
  std::uint8_t* suffix = nullptr;
  GPUdi() void operator()(Index local) const
  {
    const Index element = sources.offsets[layer] + local;
    if constexpr (Forward) {
      raiseTargetPrefixes(*parameters, sources, tracklets, offsets, element, prefix);
    } else {
      const TrackletPathBounds partial{.prefix = nullptr, .suffix = suffix, .minClusters = 0U};
      suffix[element] = static_cast<std::uint8_t>(measurementPathSuffix(partial, *parameters, sources, tracklets, offsets, layer, element));
    }
  }
};

template <DetectorLayout Det>
struct InitializeTriplets {
  TrackletConstView tracklets;
  const IterationParameters<Det>* parameters;
  Index* prefix; ///< the scratch's wide prefix, for the forward sweep's atomic maxima
  std::uint8_t* suffix;
  Index trackletCount;
  const Index* trackletOffsets;
  const Index* permutation;
  Index firstCount;
  std::array<Index, LayerRowMap<Det>::kLayerCount + 1U> sourceOffsets;
  Index* firstRanges;
  GPUdi() void operator()(Index edge) const
  {
    if (edge < sourceOffsets.size()) {
      const Index first = trackletOffsets[sourceOffsets[edge]];
      firstRanges[edge] = o2::gpu::CAAlgo::partitionPoint(Index{0U}, firstCount, [&](Index index) { return permutation[index] < first; });
    }
    if (edge >= trackletCount) {
      return;
    }
    prefix[edge] = 2U;
    suffix[edge] = parameters->roadEndLayers == 0U || ((parameters->roadEndLayers >> parameters->linkTargets[tracklets.link[edge]]) & 1U) != 0U ? 1U : 0U;
  }
};

// Dynamic ranges are read on the device, without downloading layer counts.
// Both sweeps over the layer's permuted elements, one thread each: forward,
// each triplet raises its second tracklet's prefix (wide atomics; every
// triplet ending on the element's tracklet starts on an earlier layer of the
// order); reverse, one writer per element (the end index is not built yet).
template <bool Forward>
GPUg() void GTRK_KERNEL_BOUNDS tripletPathsKernel(TripletConstView triplets, const Index* permutation, const Index* tripletOffsets, const Index* firstRanges,
                                                  LayerID layer, Index* prefix, std::uint8_t* suffix)
{
  const Index begin = firstRanges[layer];
  const Index end = firstRanges[static_cast<std::size_t>(layer) + 1U];
  const Index stride = blockDim.x * gridDim.x;
  for (Index local = blockIdx.x * blockDim.x + threadIdx.x; local < end - begin; local += stride) {
    const Index element = begin + local;
    const Index first = permutation[element];
    if constexpr (Forward) {
      const Index length = prefix[first] + 1U;
      for (Index edge = tripletOffsets[element]; edge < tripletOffsets[element + 1U]; ++edge) {
        atomicMax(prefix + triplets.second[edge], length);
      }
    } else {
      unsigned int length = suffix[first];
      for (Index edge = tripletOffsets[element]; edge < tripletOffsets[element + 1U]; ++edge) {
        const unsigned int after = suffix[triplets.second[edge]];
        if (after != 0U && length < after + 1U) {
          length = after + 1U;
        }
      }
      suffix[first] = static_cast<std::uint8_t>(length);
    }
  }
}

/// The triplets of each permuted element on a triplet path
/// (`TripletPaths::onPath`).
struct CountOnPath {
  TripletConstView triplets;
  const Index* offsets;
  TrackletPathBounds paths;
  Index elements;
  Index* counts;
  GPUdi() void operator()(Index element) const
  {
    Index kept = 0U;
    if (element < elements) {
      for (Index edge = offsets[element]; edge < offsets[element + 1U]; ++edge) {
        kept += paths.onPath(triplets.first[edge], triplets.second[edge]) ? 1U : 0U;
      }
    }
    counts[element] = kept;
  }
};

struct OnPath {
  TripletConstView triplets;
  TrackletPathBounds paths;
  GPUdi() bool operator()(Index edge) const { return paths.onPath(triplets.first[edge], triplets.second[edge]); }
};

/// to[i] = from[indices[i]], every triplet column.
struct GatherTriplets {
  TripletConstView from;
  const Index* indices;
  DeviceColumns<TripletT>::View to;
  GPUdi() void operator()(Index i) const
  {
    const Index edge = indices[i];
    to.first[i] = from.first[edge];
    to.second[i] = from.second[edge];
    to.curvature[i] = from.curvature[edge];
    to.varCurvature[i] = from.varCurvature[edge];
    to.level[i] = from.level[edge];
    to.kinked[i] = from.kinked[edge];
  }
};

struct LayerTriplets {
  const Index* layerElements;
  const Index* offsets;
  Index* layerTriplets;
  GPUdi() void operator()(Index layer) const { layerTriplets[layer] = offsets[layerElements[layer]]; }
};
} // namespace path_detail

/// TrackletPaths::build on the device: per source element of the slice,
/// stream-ordered layer sweeps, no host synchronisation or transfer between
/// layers. The prefix's maxima are 32-bit atomics on the scratch end, packed
/// to bytes at the end.
template <DetectorLayout Det>
[[nodiscard]] StageStatus trackletPaths(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, const IterationParameters<Det>* parameters,
                                        const LayerRowMap<Det>& sources, TrackletConstView tracklets, const std::uint32_t* offsets, DevicePathBounds& bounds)
{
  static_assert(maxChainLayers<Det> <= 255U);
  const auto count = sources.size();
  if (!bounds.allocate(arena, count)) {
    return StageStatus::kOutOfMemory;
  }
  if (count == 0U) {
    return StageStatus::kComplete;
  }
  auto scratch = arena.scope(ArenaEnd::kScratch);
  DeviceSpan<std::uint32_t> prefix;
  if (!prefix.allocate(arena, count, ArenaEnd::kScratch)) {
    return StageStatus::kOutOfMemory;
  }
  algorithms::forEach(runtime, stream, count, path_detail::InitializePrefix{prefix.data()});
  const auto sweep = [&]<bool Forward>(LayerID layer) {
    const auto size = sources.ranges[layer].size();
    if (size == 0U) {
      return;
    }
    algorithms::forEach(runtime, stream, size,
                        path_detail::MeasurementPaths<Det, Forward>{tracklets, offsets, sources, parameters, layer, prefix.data(), bounds.suffix.data()});
  };
  for (const auto layer : layerDepthGroups<Det>.order) {
    sweep.template operator()<true>(layer);
  }
  for (auto order = layerDepthGroups<Det>.order.size(); order-- > 0U;) {
    sweep.template operator()<false>(layerDepthGroups<Det>.order[order]);
  }
  algorithms::forEach(runtime, stream, count, path_detail::PackPrefix{prefix.data(), bounds.prefix.data()});
  return StageStatus::kComplete;
}

/// TripletPaths::build on the device, from the triplets alone, before their
/// end index: per tracklet its prefix and suffix, and each layer's range of
/// \a permutation (the triplet stage's first tracklets, in tracklet order) in
/// `layerElements`.
template <DetectorLayout Det>
[[nodiscard]] StageStatus tripletPaths(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, const IterationParameters<Det>* parameters,
                                       const LayerRowMap<Det>& sources, TrackletConstView tracklets, std::uint32_t trackletCount,
                                       const std::uint32_t* trackletOffsets, TripletConstView triplets, const std::uint32_t* tripletOffsets,
                                       const std::uint32_t* permutation, std::uint32_t firstCount, DevicePathBounds& bounds)
{
  static_assert(maxChainLayers<Det> <= 255U);
  constexpr auto kBoundaries = static_cast<std::uint32_t>(LayerRowMap<Det>::kLayerCount + 1U);
  if (!bounds.allocate(arena, trackletCount) || !bounds.layerElements.allocate(arena, kBoundaries) || !bounds.layerTriplets.allocate(arena, kBoundaries)) {
    return StageStatus::kOutOfMemory;
  }
  if (trackletCount == 0U) {
    return StageStatus::kComplete;
  }
  auto scratch = arena.scope(ArenaEnd::kScratch);
  DeviceSpan<std::uint32_t> prefix;
  if (!prefix.allocate(arena, trackletCount, ArenaEnd::kScratch)) {
    return StageStatus::kOutOfMemory;
  }
  // Each layer's range of the permutation, worked out in the initialisation's
  // launch, once per boundary.
  algorithms::forEach(runtime, stream, std::max(trackletCount, kBoundaries),
                      path_detail::InitializeTriplets<Det>{tracklets, parameters, prefix.data(), bounds.suffix.data(), trackletCount, trackletOffsets,
                                                           permutation, firstCount, sources.offsets, bounds.layerElements.data()});
  // Grid-stride loops cover skewed layer populations without a count readback.
  const auto shape = launchFor(std::max(1U, firstCount / static_cast<std::uint32_t>(LayerRowMap<Det>::kLayerCount)));
  const auto sweep = [&]<bool Forward>(LayerID layer) {
    if (sources.offsets[layer] == sources.offsets[static_cast<std::size_t>(layer) + 1U]) {
      return;
    }
    path_detail::tripletPathsKernel<Forward><<<shape.blocks, shape.threads, 0U, nativeStream(runtime, stream)>>>(
      triplets, permutation, tripletOffsets, bounds.layerElements.data(), layer, prefix.data(), bounds.suffix.data());
    checkLaunch();
  };
  for (const auto layer : layerDepthGroups<Det>.order) {
    sweep.template operator()<true>(layer);
  }
  algorithms::forEach(runtime, stream, trackletCount, path_detail::PackPrefix{prefix.data(), bounds.prefix.data()});
  for (auto order = layerDepthGroups<Det>.order.size(); order-- > 0U;) {
    sweep.template operator()<false>(layerDepthGroups<Det>.order[order]);
  }
  return StageStatus::kComplete;
}

/// \a triplets kept where they are on a triplet path of \a paths, in arena
/// order, into \a kept with their elements' counts and offsets, and the
/// bounds' layer ranges of triplets moved to them: the host's filter before
/// the place (`Tracker::OnPathTriplets`). One readback, the kept count.
[[nodiscard]] inline StageStatus onPathTriplets(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, TrackletPathBounds paths,
                                                const DeviceStageOutput<TripletT>& triplets, DeviceStageOutput<TripletT>& kept, DevicePathBounds& bounds)
{
  const auto elements = triplets.elementCount();
  if (!kept.prepare(arena, elements)) {
    return StageStatus::kOutOfMemory;
  }
  const TripletConstView view = triplets.view();
  algorithms::forEach(runtime, stream, elements + 1U,
                      path_detail::CountOnPath{view, triplets.offsets().data(), paths, elements, kept.counts().data()});
  auto status = algorithms::exclusiveScan(runtime, arena, stream, std::as_const(kept.counts()).data(), kept.offsets().data(), std::size_t{elements} + 1U);
  if (status != StageStatus::kComplete) {
    return status;
  }
  auto scratch = arena.scope(ArenaEnd::kScratch);
  DeviceSpan<std::uint32_t> indices;
  std::uint32_t total = 0U;
  status = algorithms::selectIndices(runtime, arena, stream, triplets.size(), path_detail::OnPath{view, paths}, indices, total, ArenaEnd::kScratch);
  if (status != StageStatus::kComplete) {
    return status;
  }
  if (!kept.allocateOutput(arena, total)) {
    return StageStatus::kOutOfMemory;
  }
  algorithms::forEach(runtime, stream, total, path_detail::GatherTriplets{view, indices.data(), kept.view()});
  kept.setSize(total);
  algorithms::forEach(runtime, stream, static_cast<std::uint32_t>(bounds.layerElements.size()),
                      path_detail::LayerTriplets{bounds.layerElements.data(), kept.offsets().data(), bounds.layerTriplets.data()});
  return StageStatus::kComplete;
}
} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::glue

#endif // GTRK_TRACKING_GPU_PATHGLUE_CUH_
