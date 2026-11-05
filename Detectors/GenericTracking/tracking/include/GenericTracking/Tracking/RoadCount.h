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

/// \file RoadCount.h
/// \brief Where roads start and how many leave each triplet.
///
/// A road is a chain of triplets `T_1, ..., T_L`, `T_i` a predecessor of
/// `T_{i+1}` and at level `i`. It starts at a triplet without a successor
/// exactly one level above, and descends through predecessors exactly one
/// level below, so every triplet lies on one; an iteration with an end-layer
/// mask starts them instead at every triplet ending on one of its layers
/// (`RoadStartRule`). The number of roads through a
/// triplet down to level 1 is the sum over those predecessors; the base case
/// counts only a level-1 triplet whose first measurement -- the road's
/// innermost -- is in the slice's core, so only the roads the core owns are
/// counted.

#ifndef GTRK_TRACKING_ROADCOUNT_H_
#define GTRK_TRACKING_ROADCOUNT_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Base/HostAlgorithms.h"
#include "GenericTracking/Tracking/MeasurementRange.h"
#include "GenericTracking/Tracking/Neighbour.h"
#include "GenericTracking/Tracking/OrderedRoadCounts.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/Profile.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/Triplet.h"

#include "GPUCommonAlgorithm.h"

#include <array>
#include <cstdint>
#include <memory>
#include <span>

namespace o2::gtrk
{

/// The core ranges the road counting reads (`OrderedRoadCounts::base`):
/// \a core with the layers a road may not start on emptied
/// (`IterationParameters::roadStartLayers`, 0 for any), so a road whose
/// innermost measurement lies there is none, as one starting outside the core.
template <std::size_t N>
[[nodiscard]] inline std::array<MeasurementRange, N> roadCores(std::span<const MeasurementRange, N> core, std::uint64_t startLayers)
{
  std::array<MeasurementRange, N> cores{};
  for (std::size_t layer = 0U; layer < N; ++layer) {
    cores[layer] = startLayers == 0U || ((startLayers >> layer) & 1U) != 0U ? core[layer] : MeasurementRange{core[layer].begin, core[layer].begin};
  }
  return cores;
}

/// Which triplets roads start at: at level `minLevel` or higher, with a road,
/// and -- without an end-layer mask -- without a successor exactly one level
/// above, or -- with one -- with their last measurement on a layer of the
/// mask, whatever extends them (`IterationSettings::roadEndLayers`).
struct RoadStartRule {
  std::uint32_t minLevel = 2U;
  std::uint64_t endLayers = 0U;         ///< a bit per layer; 0 for any
  const LayerID* linkTargets = nullptr; ///< per link; read only with a mask

  GPUhdi() GPUAlwaysInline() bool operator()(TrackletConstView tracklets, TripletConstView triplets, const std::uint32_t* successorAbove,
                                             const std::uint32_t* roads, StageOutput<TripletT>::Index triplet) const
  {
    if (triplets.level[triplet] < minLevel || roads[triplet] == 0U) {
      return false;
    }
    if (endLayers == 0U) {
      return successorAbove[triplet] == 0U;
    }
    const LayerID end = linkTargets[tracklets.link[triplets.second[triplet]]];
    return ((endLayers >> end) & 1U) != 0U;
  }
};

/// An iteration's rule, reading its link targets on the host.
template <DetectorLayout Det>
[[nodiscard]] inline RoadStartRule roadStartRule(const IterationParameters<Det>& parameters)
{
  return {.minLevel = parameters.minRoadStartLevel, .endLayers = parameters.roadEndLayers, .linkTargets = parameters.linkTargets};
}

/// One road by its index: the start by a binary search in the starts' road
/// offsets, then down level by level through the predecessors exactly one
/// level below, in the neighbour stage's order, subtracting the roads of each
/// one passed over until the index falls inside one (as GBTS's
/// `gbts_fill_path_store`). The walk's steps are the road's levels, each a
/// scan of one triplet's predecessors; road `r` is the same road on every
/// backend.
class RoadWalk
{
 public:
  using Index = StageOutput<TripletT>::Index;

  RoadWalk(const Index* starts, const Index* roadOffsets, Index startCount, const Index* neighbourOffsets, const Index* from, const std::uint32_t* level,
           const std::uint32_t* roads)
    : mStarts(starts), mRoadOffsets(roadOffsets), mStartCount(startCount), mNeighbourOffsets(neighbourOffsets), mFrom(from), mLevel(level), mRoads(roads)
  {
  }

  /// Write road \a road's triplets to \a triplets, innermost first:
  /// `triplets[i]` is the one at level `i + 1`. \a triplets holds the start's
  /// level, at most `maxChainLayers - 2`. Returns that level.
  GPUhdi() std::uint32_t operator()(Index road, Index* triplets) const
  {
    const Index first = o2::gpu::CAAlgo::partitionPoint(Index{0U}, mStartCount, [&](Index start) { return mRoadOffsets[start + 1U] <= road; });
    Index current = mStarts[first];
    std::uint32_t remainder = road - mRoadOffsets[first];
    const std::uint32_t length = mLevel[current];
    triplets[length - 1U] = current;
    for (std::uint32_t level = length; level > 1U; --level) {
      for (Index edge = mNeighbourOffsets[current]; edge < mNeighbourOffsets[current + 1U]; ++edge) {
        const Index predecessor = mFrom[edge];
        if (mLevel[predecessor] + 1U != level) {
          continue;
        }
        if (remainder < mRoads[predecessor]) {
          current = predecessor;
          break;
        }
        remainder -= mRoads[predecessor];
      }
      triplets[level - 2U] = current;
    }
    return length;
  }

 private:
  const Index* mStarts = nullptr;
  const Index* mRoadOffsets = nullptr;
  Index mStartCount = 0U;
  const Index* mNeighbourOffsets = nullptr;
  const Index* mFrom = nullptr;
  const std::uint32_t* mLevel = nullptr;
  const std::uint32_t* mRoads = nullptr;
};

/// The successor flags, the road counts, the start triplets of one slice's
/// roads and their numbering.
class RoadCount
{
 public:
  using Index = StageOutput<TripletT>::Index;

  RoadCount() = default;
  ~RoadCount() = default;
  RoadCount(const RoadCount&) = delete;
  RoadCount& operator=(const RoadCount&) = delete;
  RoadCount(RoadCount&&) = delete;
  RoadCount& operator=(RoadCount&&) = delete;

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool.swap(pool);
    clearStorage();
  }

  /// Storage for counts finalized by ordered neighbour generation. Every
  /// triplet writes its count before a successor can read it; only the shared
  /// flags need an initial fill. Failures are handled by the caller's guard.
  template <typename ForEach>
  void prepareOrdered(ForEach forEach, Index count)
  {
    const profile::Zone zone{profile::Phase::kRoadFills};
    resizeBoundedVector(mSuccessorAbove, count);
    resizeBoundedVector(mRoads, count);
    mOrderedOverflow = 0U;
    forEach(count, [&](Index triplet) { mSuccessorAbove[triplet] = 0U; });
  }

  [[nodiscard]] OrderedRoadCounts ordered(const LayerID* linkSources, const MeasurementRange* core)
  {
    return {.roads = mRoads.data(), .successorAbove = mSuccessorAbove.data(), .overflow = &mOrderedOverflow, .linkSources = linkSources, .core = core};
  }

  /// Counts and successor flags were finalized during neighbour generation.
  /// Select and number the starts from the counts the neighbour generation
  /// finalized.
  template <typename ForEach>
  [[nodiscard]] StageStatus runPrepared(ForEach forEach, const RoadStartRule& startRule, const StageOutput<TrackletT>& tracklets,
                                        const StageOutput<TripletT>& triplets, Index chunkSize)
  {
    return guardOutOfMemory([&] {
      if (mOrderedOverflow != 0U) {
        return StageStatus::kIndexOverflow;
      }
      return number(forEach, startRule, tracklets, triplets, chunkSize);
    });
  }

  void clearStorage() noexcept { deepVectorClear(mMemoryPool.get(), mSuccessorAbove, mRoads, mStarts, mRoadOffsets); }

  /// Per triplet: 1 if a successor sits exactly one level above it.
  [[nodiscard]] std::span<const std::uint32_t> successorAbove() const { return mSuccessorAbove; }
  /// Per triplet: the roads the core owns from it down to level 1.
  [[nodiscard]] std::span<const std::uint32_t> roads() const { return mRoads; }
  /// The triplets roads start at, in arena order.
  [[nodiscard]] std::span<const Index> starts() const { return mStarts; }
  /// Per start, the index of its first road, and the total at the end.
  [[nodiscard]] std::span<const Index> roadOffsets() const { return mRoadOffsets; }
  /// The number of roads the slice owns.
  [[nodiscard]] Index roadCount() const { return mRoadOffsets.empty() ? 0U : mRoadOffsets.back(); }

  /// The walk to any of them, over the triplets and neighbours counted.
  [[nodiscard]] RoadWalk walk(const StageOutput<TripletT>& triplets, const StageOutput<NeighbourT>& neighbours) const
  {
    return {mStarts.data(), mRoadOffsets.data(), static_cast<Index>(mStarts.size()), neighbours.offsets().data(), neighbours.view().from,
            triplets.view().level, mRoads.data()};
  }

 private:
  template <typename ForEach>
  [[nodiscard]] StageStatus number(ForEach forEach, const RoadStartRule& startRule, const StageOutput<TrackletT>& tracklets,
                                   const StageOutput<TripletT>& triplets, Index chunkSize)
  {
    const profile::Zone zone{profile::Phase::kRoadStarts};
    const auto view = triplets.view();
    const auto trackletView = tracklets.view();
    host_algorithms::compactIndices(
      forEach, triplets.size(), chunkSize, [&](Index triplet) { return startRule(trackletView, view, mSuccessorAbove.data(), mRoads.data(), triplet); },
      mStarts);
    const auto starts = static_cast<Index>(mStarts.size());
    resizeBoundedVector(mRoadOffsets, std::size_t{starts} + 1U);
    if (!host_algorithms::exclusiveScan(forEach, starts, chunkSize, [&](Index start) { return mRoads[mStarts[start]]; }, mRoadOffsets, mMemoryPool.get())) {
      return StageStatus::kIndexOverflow;
    }
    return StageStatus::kComplete;
  }

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  BoundedVector<std::uint32_t> mSuccessorAbove;
  BoundedVector<std::uint32_t> mRoads;
  BoundedVector<Index> mStarts;
  BoundedVector<Index> mRoadOffsets;
  Index mOrderedOverflow = 0U;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_ROADCOUNT_H_
