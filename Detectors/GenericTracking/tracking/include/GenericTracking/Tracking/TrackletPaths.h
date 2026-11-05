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

/// \file TrackletPaths.h
/// \brief The longest tracklet paths through and from each measurement, which
///        bound the chains a triplet can be on.

#ifndef GTRK_TRACKING_TRACKLETPATHS_H_
#define GTRK_TRACKING_TRACKLETPATHS_H_

#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LayerRowMap.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/Tracklet.h"

#include "GPUCommonMath.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

namespace o2::gtrk
{

/// Upper bounds from the actual tracklet DAG, ignoring triplet/neighbour cuts.
/// A triplet a -> b -> c can be on a chain the
/// iteration uses only when a longest prefix ending at a, these two links, and
/// a suffix from c to an allowed road end together contain at least the
/// iteration's minClusters measurements. TrackletPaths indexes its bounds by
/// the slice's source element (LayerRowMap::element), TripletPaths by tracklet.
struct TrackletPathBounds {
  const std::uint8_t* prefix = nullptr;
  const std::uint8_t* suffix = nullptr;
  std::uint32_t minClusters = 0U;

  /// Whether a chain through \a before's prefix and \a after's suffix can
  /// reach minClusters.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool onPath(std::uint32_t before, std::uint32_t after) const
  {
    return suffix[after] != 0U && std::uint32_t{prefix[before]} + suffix[after] >= minClusters;
  }
};

/// Whether a road may end on \a layer (`IterationParameters::roadEndLayers`, 0 for any).
GPUhdi() GPUAlwaysInline() bool roadMayEnd(std::uint64_t endLayers, LayerID layer) { return endLayers == 0U || ((endLayers >> layer) & 1U) != 0U; }

/// The measurement suffix of \a row on \a layer: its source element's, or,
/// for a row that is no source of the slice and so has no tracklets, one
/// cluster where a road may end and none elsewhere.
template <DetectorLayout Det>
GPUhdi() GPUAlwaysInline() std::uint32_t measurementSuffix(const TrackletPathBounds& paths, const LayerRowMap<Det>& sources, std::uint64_t endLayers,
                                                           LayerID layer, std::uint32_t row)
{
  if (sources.contains(layer, row)) {
    return paths.suffix[sources.element(layer, row)];
  }
  return roadMayEnd(endLayers, layer) ? 1U : 0U;
}

/// The forward sweep's step from source element \a element: each of its
/// tracklets raises its target's prefix to the element's plus one, by an
/// atomic maximum (the sweep runs in parallel). \a Prefix is a byte on the
/// host and a 32-bit word on a device, which has no byte atomics.
template <DetectorLayout Det, typename Prefix>
GPUdi() GPUAlwaysInline() void raiseTargetPrefixes(const IterationParameters<Det>& parameters, const LayerRowMap<Det>& sources, TrackletConstView tracklets,
                                                   const std::uint32_t* offsets, std::uint32_t element, Prefix* prefix)
{
  const auto length = static_cast<Prefix>(prefix[element] + 1U);
  for (std::uint32_t edge = offsets[element]; edge < offsets[element + 1U]; ++edge) {
    const LayerID targetLayer = parameters.linkTargets[tracklets.link[edge]];
    if (sources.contains(targetLayer, tracklets.target[edge])) {
      o2::gpu::CAMath::AtomicMax(prefix + sources.element(targetLayer, tracklets.target[edge]), length);
    }
  }
}

/// The reverse sweep's suffix of source element \a element on \a layer: one
/// cluster where a road may end, or one more than the longest suffix of its
/// tracklets' targets, which \a partial holds final.
template <DetectorLayout Det>
GPUhdi() GPUAlwaysInline() std::uint32_t measurementPathSuffix(const TrackletPathBounds& partial, const IterationParameters<Det>& parameters,
                                                               const LayerRowMap<Det>& sources, TrackletConstView tracklets, const std::uint32_t* offsets,
                                                               LayerID layer, std::uint32_t element)
{
  std::uint32_t length = roadMayEnd(parameters.roadEndLayers, layer) ? 1U : 0U;
  for (std::uint32_t edge = offsets[element]; edge < offsets[element + 1U]; ++edge) {
    const std::uint32_t after =
      measurementSuffix(partial, sources, parameters.roadEndLayers, parameters.linkTargets[tracklets.link[edge]], tracklets.target[edge]);
    if (after != 0U && length < after + 1U) {
      length = after + 1U;
    }
  }
  return length;
}

/// `body(layer, index)` for every index of the layers of depth group \a group
/// of the link DAG (LayerLinks.h, `layerDepthGroups`), each layer's indices in
/// [rangeBegin(layer), rangeEnd(layer)), in one \a forEach over all of them.
/// No link joins two layers of one depth and every link into them comes from
/// a smaller depth, so a sweep's invariant -- a layer's predecessors final
/// before it -- holds depth by depth, forward or back.
template <DetectorLayout Det, typename ForEach, typename Begin, typename End, typename Body>
void forEachOfDepth(const ForEach& forEach, std::size_t group, const Begin& rangeBegin, const End& rangeEnd, const Body& body)
{
  using Index = StageOutput<TrackletT>::Index;
  constexpr auto groups = layerDepthGroups<Det>;
  std::array<LayerID, LayerDepthGroups<Det>::kLayers> layers{};
  std::array<Index, LayerDepthGroups<Det>::kLayers> bases{};
  std::array<Index, LayerDepthGroups<Det>::kLayers + 1U> starts{};
  std::size_t count = 0U;
  for (auto at = groups.groups[group]; at < groups.groups[group + 1U]; ++at) {
    layers[count] = groups.order[at];
    bases[count] = rangeBegin(layers[count]);
    starts[count + 1U] = starts[count] + (rangeEnd(layers[count]) - bases[count]);
    ++count;
  }
  forEach(starts[count], [&](Index index) {
    std::size_t i = 0U;
    while (index >= starts[i + 1U]) {
      ++i;
    }
    body(layers[i], static_cast<Index>(bases[i] + (index - starts[i])));
  });
}

/// The bounds of a slice, built by two topological sweeps over its sources;
/// no fits and no new physics threshold. Removing a path that fails this
/// necessary condition cannot remove an accepted road. Bounds deliberately
/// ignore time-box intersections and core ownership.
template <DetectorLayout Det>
class TrackletPaths
{
 public:
  using Index = StageOutput<TrackletT>::Index;
  static_assert(maxChainLayers<Det> <= 255U);

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mPool = std::move(pool);
    clearStorage();
  }

  void clearStorage() noexcept { deepVectorClear(mPool.get(), mPrefix, mSuffix); }

  /// \a forEach may run in parallel: the prefix's maxima are atomic.
  template <typename ForEach>
  void build(const IterationParameters<Det>& parameters, const LayerRowMap<Det>& sources, const StageOutput<TrackletT>& tracklets, const ForEach& forEach)
  {
    // Kept from slice to slice and filled in parallel, not reallocated:
    // the prefix's floor of one measurement.
    resizeBoundedVector(mPrefix, sources.size());
    resizeBoundedVector(mSuffix, sources.size());
    forEach(static_cast<Index>(sources.size()), [&](Index element) { mPrefix[element] = std::uint8_t{1U}; });
    mMinClusters = parameters.minClusters;
    const auto view = tracklets.view();
    const auto offsets = tracklets.offsets();
    // The layers of one depth of the link DAG in one pass (`forEachOfDepth`).
    constexpr auto groups = layerDepthGroups<Det>;
    const auto elementBegin = [&](LayerID layer) { return static_cast<Index>(sources.offsets[layer]); };
    const auto elementEnd = [&](LayerID layer) { return static_cast<Index>(sources.offsets[layer] + sources.ranges[layer].size()); };
    for (std::size_t group = 0U; group < groups.count; ++group) {
      forEachOfDepth<Det>(forEach, group, elementBegin, elementEnd,
                          [&](LayerID /*layer*/, Index element) { raiseTargetPrefixes(parameters, sources, view, offsets.data(), element, mPrefix.data()); });
    }
    const TrackletPathBounds partial{.prefix = mPrefix.data(), .suffix = mSuffix.data(), .minClusters = mMinClusters};
    for (auto group = groups.count; group-- > 0U;) {
      forEachOfDepth<Det>(forEach, group, elementBegin, elementEnd, [&](LayerID layer, Index element) {
        mSuffix[element] = static_cast<std::uint8_t>(measurementPathSuffix(partial, parameters, sources, view, offsets.data(), layer, element));
      });
    }
  }

  [[nodiscard]] TrackletPathBounds bounds() const { return {.prefix = mPrefix.data(), .suffix = mSuffix.data(), .minClusters = mMinClusters}; }

 private:
  std::shared_ptr<BoundedMemoryResource> mPool;
  BoundedVector<std::uint8_t> mPrefix;
  BoundedVector<std::uint8_t> mSuffix;
  std::uint32_t mMinClusters = 0U;
};

} // namespace o2::gtrk
#endif // GTRK_TRACKING_TRACKLETPATHS_H_
