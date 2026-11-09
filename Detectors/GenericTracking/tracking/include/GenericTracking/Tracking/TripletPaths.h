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

/// \file TripletPaths.h
/// \brief The longest triplet paths through and from each tracklet, which bound
///        the chains a neighbour can be on.

#ifndef GTRK_TRACKING_TRIPLETPATHS_H_
#define GTRK_TRACKING_TRIPLETPATHS_H_

#include "GenericTracking/Tracking/TrackletPaths.h"
#include "GenericTracking/Tracking/Triplet.h"
#include "GenericTracking/Tracking/TripletPermutation.h"

#include "GPUCommonAlgorithm.h"

#include <array>
#include <atomic>
#include <span>

namespace o2::gtrk
{
/// A necessary condition on the triplets, as `TrackletPaths` is on the
/// tracklets: the DAG whose vertices are tracklets
/// and whose edges are accepted triplets, before the expensive neighbour fit.
/// Prefix counts measurements up to a tracklet's target; suffix counts from
/// that target to an allowed road end, including that target. Built from the
/// triplets the triplet stage staged, before they are placed: the stage keeps
/// only those on a path (`onPath`).
template <DetectorLayout Det>
class TripletPaths
{
 public:
  using Index = StageOutput<TrackletT>::Index;
  static constexpr std::size_t kLayers = LayerRowMap<Det>::kLayerCount;
  static_assert(maxChainLayers<Det> <= 255U);

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mPool = std::move(pool);
    clearStorage();
  }

  void clearStorage() noexcept { deepVectorClear(mPool.get(), mPrefix, mSuffix); }

  /// \a forEachTriplet(begin, end, body) calls body(element, first, second)
  /// for every triplet of the permuted elements in [begin, end), each
  /// element's on one thread. The prefix by the link DAG's depth groups of
  /// the first tracklets' source layers, each triplet raising its second
  /// tracklet's (an atomic maximum): every triplet ending on a tracklet starts
  /// at a smaller depth, as every link raises the depth, so the prefix of an
  /// element's tracklet is final when its group runs; the suffix in reverse,
  /// one writer per element.
  template <typename ForEachTriplet, typename ForEach>
  void build(const IterationParameters<Det>& parameters, const LayerRowMap<Det>& sources, const StageOutput<TrackletT>& tracklets, const TripletPermutation<Det>& permutation, const ForEachTriplet& forEachTriplet, const ForEach& forEach)
  {
    mMinClusters = parameters.minClusters;
    // Kept from slice to slice and filled in parallel, not reallocated:
    // the prefix's floor of two measurements
    // and the suffix's of one, which with road-end layers is whether the
    // tracklet's own target may end a road.
    resizeBoundedVector(mPrefix, tracklets.size());
    resizeBoundedVector(mSuffix, tracklets.size());
    if (parameters.roadEndLayers != 0U) {
      const auto view = tracklets.view();
      std::array<std::uint8_t, nLinks<Det>> allowed{};
      for (std::size_t link = 0U; link < nLinks<Det>; ++link) {
        allowed[link] = ((parameters.roadEndLayers >> parameters.linkTargets[link]) & 1U) != 0U ? 1U : 0U;
      }
      forEach(tracklets.size(), [&](Index edge) {
        mPrefix[edge] = std::uint8_t{2U};
        mSuffix[edge] = allowed[view.link[edge]];
      });
    } else {
      forEach(tracklets.size(), [&](Index edge) {
        mPrefix[edge] = std::uint8_t{2U};
        mSuffix[edge] = std::uint8_t{1U};
      });
    }
    const auto indices = permutation.indices();
    for (std::size_t layer = 0U; layer <= kLayers; ++layer) {
      const auto first = tracklets.offsets()[sources.offsets[layer]];
      mElements[layer] = o2::gpu::CAAlgo::partitionPoint(Index{0U}, permutation.size(), [&](Index element) { return indices[element] < first; });
    }
    constexpr auto groups = layerDepthGroups<Det>;
    const auto forEachOfGroup = [&](std::size_t group, const auto& body) {
      for (auto at = groups.groups[group]; at < groups.groups[group + 1U]; ++at) {
        const auto layer = static_cast<std::size_t>(groups.order[at]);
        forEachTriplet(mElements[layer], mElements[layer + 1U], body);
      }
    };
    for (std::size_t group = 0U; group < groups.count; ++group) {
      forEachOfGroup(group, [&](Index /*element*/, Index first, Index second) {
        // A host loop: O2's AtomicMax is device code in a GPU unit.
        const auto length = static_cast<std::uint8_t>(mPrefix[first] + 1U);
        const std::atomic_ref<std::uint8_t> prefix{mPrefix[second]};
        auto current = prefix.load(std::memory_order_relaxed);
        while (current < length && !prefix.compare_exchange_weak(current, length, std::memory_order_relaxed)) {
        }
      });
    }
    for (auto group = groups.count; group-- > 0U;) {
      forEachOfGroup(group, [&](Index /*element*/, Index first, Index second) {
        const auto after = mSuffix[second];
        if (after != 0U && mSuffix[first] < after + 1U) {
          mSuffix[first] = static_cast<std::uint8_t>(after + 1U);
        }
      });
    }
  }

  /// Whether triplet (\a first, \a second) is on a triplet path of the
  /// iteration's minClusters: the neighbour stage's test of its element
  /// (NeighbourStage.h), which refuses every other triplet as element and as
  /// predecessor.
  [[nodiscard]] bool onPath(Index first, Index second) const { return bounds().onPath(first, second); }

  [[nodiscard]] TrackletPathBounds bounds() const { return {.prefix = mPrefix.data(), .suffix = mSuffix.data(), .minClusters = mMinClusters}; }

  /// Per layer, the first permuted element whose tracklet starts on it or
  /// after, and the total at the end: the layers' ranges of the permutation,
  /// as the last `build` found them.
  [[nodiscard]] std::span<const Index, kLayers + 1U> layerElements() const { return mElements; }

 private:
  std::shared_ptr<BoundedMemoryResource> mPool;
  BoundedVector<std::uint8_t> mPrefix;
  BoundedVector<std::uint8_t> mSuffix;
  std::array<Index, kLayers + 1U> mElements{};
  std::uint32_t mMinClusters = 0U;
};
} // namespace o2::gtrk
#endif // GTRK_TRACKING_TRIPLETPATHS_H_
