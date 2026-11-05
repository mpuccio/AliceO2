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

/// \file FixedPointReferences.h
/// \brief The levels and the road counts by sweeps to their fixed point, the
///        independent path the tests hold the ordered one to.
///
/// The tracker writes both during the neighbour stage, a stratum at a time
/// (`NeighbourStage` with its ready levels and `OrderedRoadCounts`). Here
/// they come from the whole neighbour graph instead: every triplet's level
/// from its predecessors' until nothing rises, and
/// the roads through it from those exactly one level below, the same way.
/// The sweep functors are the device's as well
/// (DeviceTestKernels.cu).

#ifndef GTRK_TEST_FIXEDPOINTREFERENCES_H_
#define GTRK_TEST_FIXEDPOINTREFERENCES_H_

#include "GenericTracking/Tracking/Neighbour.h"
#include "GenericTracking/Tracking/RoadCount.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/Triplet.h"

#include "GPUCommonMath.h"

#include <cstdint>
#include <vector>

namespace o2::gtrk
{

/// One element per triplet: the base case, one road at level 1 if the
/// triplet's first measurement is in the core, none otherwise, and none yet
/// above level 1.
class RoadBase
{
 public:
  using Index = StageOutput<TripletT>::Index;

  RoadBase(TrackletConstView tracklets, TripletConstView triplets, const LayerID* linkSources, const MeasurementRange* core, std::uint32_t* roads)
    : mTracklets(tracklets), mTriplets(triplets), mLinkSources(linkSources), mCore(core), mRoads(roads)
  {
  }

  GPUhdi() void operator()(Index triplet) const
  {
    mRoads[triplet] = mTriplets.level[triplet] == 1U && roadOwned(mTracklets, mTriplets, mLinkSources, mCore, triplet) ? 1U : 0U;
  }

 private:
  TrackletConstView mTracklets;
  TripletConstView mTriplets;
  const LayerID* mLinkSources = nullptr;
  const MeasurementRange* mCore = nullptr;
  std::uint32_t* mRoads = nullptr;
};

/// The triplets with predecessors, in arena order: the elements the sweeps
/// visit.
inline std::vector<StageOutput<TripletT>::Index> chainedTriplets(const StageOutput<NeighbourT>& neighbours)
{
  const auto offsets = neighbours.offsets();
  std::vector<StageOutput<TripletT>::Index> chained;
  for (StageOutput<TripletT>::Index triplet = 0U; triplet + 1U < offsets.size(); ++triplet) {
    if (offsets[triplet + 1U] != offsets[triplet]) {
      chained.push_back(triplet);
    }
  }
  return chained;
}

/// One element of a level sweep: a triplet with predecessors pulls their
/// highest level and raises its own to one above it.
///
/// Levels are updated in place, one buffer, with an atomic max: a sweep may
/// read a predecessor already raised in the same sweep, so how many sweeps
/// the levels take depends on the schedule, but the fixed point -- a longest
/// path, a maximum -- does not. Only the triplet itself raises its level;
/// the atomic makes the store safe against the concurrent reads.
class LevelUpdate
{
 public:
  using Index = StageOutput<TripletT>::Index;
  using Level = std::uint32_t;

  LevelUpdate(const Index* chained, const Index* neighbourOffsets, const Index* from, Level* level, Index* changed)
    : mChained(chained), mNeighbourOffsets(neighbourOffsets), mFrom(from), mLevel(level), mChanged(changed)
  {
  }

  GPUdi() void operator()(Index element) const
  {
    const Index self = mChained[element];
    const Index begin = mNeighbourOffsets[self];
    const Index end = mNeighbourOffsets[self + 1U];
    Level highest = 0U;
    for (Index neighbour = begin; neighbour < end; ++neighbour) {
      const Level level = load(mFrom[neighbour]);
      highest = level > highest ? level : highest;
    }
    const Level raised = highest + 1U;
    if (raised > load(self)) {
      o2::gpu::CAMath::AtomicMax(mLevel + self, raised);
      o2::gpu::CAMath::AtomicMax(mChanged, Index{1U});
    }
  }

 private:
  /// A level another element may be raising (`CAMath::AtomicLoad`).
  [[nodiscard]] GPUdi() Level load(Index triplet) const { return o2::gpu::CAMath::AtomicLoad(mLevel + triplet); }

  const Index* mChained = nullptr;
  const Index* mNeighbourOffsets = nullptr;
  const Index* mFrom = nullptr;
  Level* mLevel = nullptr;
  Index* mChanged = nullptr;
};

/// One element per triplet `S`: every predecessor exactly one level below
/// gets its flag. Many successors may set one flag; the atomic max makes the
/// concurrent stores of one value safe.
class SuccessorMark
{
 public:
  using Index = StageOutput<TripletT>::Index;

  SuccessorMark(const Index* neighbourOffsets, const Index* from, const std::uint32_t* level, std::uint32_t* successorAbove)
    : mNeighbourOffsets(neighbourOffsets), mFrom(from), mLevel(level), mSuccessorAbove(successorAbove)
  {
  }

  GPUdi() void operator()(Index successor) const
  {
    const std::uint32_t level = mLevel[successor];
    for (Index edge = mNeighbourOffsets[successor]; edge < mNeighbourOffsets[successor + 1U]; ++edge) {
      const Index predecessor = mFrom[edge];
      if (mLevel[predecessor] + 1U == level) {
        o2::gpu::CAMath::AtomicMax(mSuccessorAbove + predecessor, std::uint32_t{1U});
      }
    }
  }

 private:
  const Index* mNeighbourOffsets = nullptr;
  const Index* mFrom = nullptr;
  const std::uint32_t* mLevel = nullptr;
  std::uint32_t* mSuccessorAbove = nullptr;
};

/// One element of a count sweep, over the triplets with predecessors: the
/// sum of the roads of the predecessors exactly one level below.
///
/// As the level sweep, in place: a sweep may read a predecessor already
/// raised in the same sweep. Every count starts at or below its final value
/// and a sum of such counts is too, so counts only rise, and the atomic max
/// stores them safely against the concurrent reads; the fixed point is
/// unique. A sum beyond 32 bits saturates and raises \a overflow: roads are
/// never truncated silently.
class RoadCountUpdate
{
 public:
  using Index = StageOutput<TripletT>::Index;

  RoadCountUpdate(const Index* chained, const Index* neighbourOffsets, const Index* from, const std::uint32_t* level, std::uint32_t* roads, Index* changed,
                  Index* overflow)
    : mChained(chained), mNeighbourOffsets(neighbourOffsets), mFrom(from), mLevel(level), mRoads(roads), mChanged(changed), mOverflow(overflow)
  {
  }

  GPUdi() void operator()(Index element) const
  {
    const Index self = mChained[element];
    const std::uint32_t below = mLevel[self] - 1U;
    std::uint64_t sum = 0U;
    for (Index edge = mNeighbourOffsets[self]; edge < mNeighbourOffsets[self + 1U]; ++edge) {
      const Index predecessor = mFrom[edge];
      if (mLevel[predecessor] == below) {
        sum += load(predecessor);
      }
    }
    const auto count = saturatedRoads(sum, mOverflow);
    if (count > load(self)) {
      o2::gpu::CAMath::AtomicMax(mRoads + self, count);
      o2::gpu::CAMath::AtomicMax(mChanged, Index{1U});
    }
  }

 private:
  /// A count another element may be raising (`CAMath::AtomicLoad`).
  [[nodiscard]] GPUdi() std::uint32_t load(Index triplet) const { return o2::gpu::CAMath::AtomicLoad(mRoads + triplet); }

  const Index* mChained = nullptr;
  const Index* mNeighbourOffsets = nullptr;
  const Index* mFrom = nullptr;
  const std::uint32_t* mLevel = nullptr;
  std::uint32_t* mRoads = nullptr;
  Index* mChanged = nullptr;
  Index* mOverflow = nullptr;
};

/// The levels of \a triplets from the whole graph of \a neighbours, swept
/// over the triplets with predecessors until one changes nothing;
/// \a forEach(elements, update) runs update(element) for every element,
/// serially or in parallel. The sweeps that changed a level.
template <typename ForEach>
StageOutput<TripletT>::Index fixedPointLevels(ForEach forEach, const StageOutput<NeighbourT>& neighbours, StageOutput<TripletT>& triplets)
{
  using Index = StageOutput<TripletT>::Index;
  const auto chained = chainedTriplets(neighbours);
  Index changed = 0U;
  const LevelUpdate update{chained.data(), neighbours.offsets().data(), neighbours.view().from, triplets.view().level, &changed};
  Index sweeps = 0U;
  while (true) {
    changed = 0U;
    forEach(static_cast<Index>(chained.size()), update);
    if (changed == 0U) {
      return sweeps;
    }
    ++sweeps;
  }
}

/// \a count's successor flags and road counts by sweeps to their fixed point
/// over the whole graph of \a neighbours (its triplets with predecessors),
/// then its starts by \a startRule and their numbering as the
/// tracker's (`RoadCount::runPrepared`). \a linkSources maps a link to its
/// source layer, \a core the slice's rows per layer. \a forEach as
/// `fixedPointLevels`'s.
template <typename ForEach>
[[nodiscard]] StageStatus fixedPointRoads(RoadCount& count, ForEach forEach, const LayerID* linkSources, const MeasurementRange* core,
                                          const RoadStartRule& startRule, const StageOutput<TrackletT>& tracklets, const StageOutput<TripletT>& triplets,
                                          const StageOutput<NeighbourT>& neighbours, RoadCount::Index chunkSize = 256U)
{
  using Index = RoadCount::Index;
  const auto status = guardOutOfMemory([&] {
    const Index triplet = triplets.size();
    count.prepareOrdered(forEach, triplet);
    const auto counts = count.ordered(linkSources, core);
    forEach(triplet, [&](Index t) { counts.roads[t] = 0U; });
    const auto view = triplets.view();
    const Index* const offsets = neighbours.offsets().data();
    const Index* const from = neighbours.view().from;
    forEach(triplet, SuccessorMark{offsets, from, view.level, counts.successorAbove});
    forEach(triplet, RoadBase{tracklets.view(), view, linkSources, core, counts.roads});
    Index changed = 0U;
    Index overflow = 0U;
    const auto chained = chainedTriplets(neighbours);
    const RoadCountUpdate update{chained.data(), offsets, from, view.level, counts.roads, &changed, &overflow};
    while (true) {
      changed = 0U;
      forEach(static_cast<Index>(chained.size()), update);
      if (changed == 0U) {
        break;
      }
    }
    return overflow != 0U ? StageStatus::kIndexOverflow : StageStatus::kComplete;
  });
  if (status != StageStatus::kComplete) {
    return status;
  }
  return count.runPrepared(forEach, startRule, tracklets, triplets, chunkSize);
}

} // namespace o2::gtrk

#endif // GTRK_TEST_FIXEDPOINTREFERENCES_H_
