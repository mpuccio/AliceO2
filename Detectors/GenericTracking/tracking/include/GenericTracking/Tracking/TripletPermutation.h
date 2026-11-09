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

/// \file TripletPermutation.h
/// \brief The first tracklets of the triplet stage.

#ifndef GTRK_TRACKING_TRIPLETPERMUTATION_H_
#define GTRK_TRACKING_TRIPLETPERMUTATION_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/HostAlgorithms.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LayerRowMap.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/TrackletPaths.h"

#include <memory>
#include <span>

namespace o2::gtrk
{

/// Whether tracklet \a first is one of the triplet stage's first tracklets:
/// its target a source of the slice with tracklets of its own (\a offsets,
/// the tracklets' per source), on a tracklet path of the iteration's
/// minClusters (\a paths). The one predicate of
/// TripletPermutation::build and of the device's permutation (SliceGlue.cuh).
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool firstTrackletEligible(const IterationParameters<Det>& parameters, const LayerRowMap<Det>& sources, TrackletConstView tracklets, const std::uint32_t* offsets, TrackletPathBounds paths, std::uint32_t first)
{
  const auto link = tracklets.link[first];
  const auto layer = parameters.linkTargets[link];
  if (!sources.contains(layer, tracklets.target[first])) {
    return false;
  }
  const auto middle = sources.element(layer, tracklets.target[first]);
  const auto source = sources.element(parameters.linkSources[link], tracklets.source[first]);
  return offsets[middle] != offsets[middle + 1U] && paths.onPath(source, middle);
}

/// The first tracklets consumed by TripletStage, in arena order: the
/// tracklets whose target is a source of the slice with tracklets of its own,
/// on a tracklet path of the iteration's minClusters (\a paths). Every tracklet's own source already passed the
/// tracklet stage's gate, so this is a compaction, and the order is a pure
/// function of the inputs.
template <DetectorLayout Det>
class TripletPermutation
{
 public:
  using Index = typename StageOutput<TrackletT>::Index;
  using Sources = LayerRowMap<Det>;
  using Tracklets = StageOutput<TrackletT>;
  static_assert(sizeof(Index) == 4U);

  TripletPermutation() = default;
  ~TripletPermutation() = default;
  TripletPermutation(const TripletPermutation&) = delete;
  TripletPermutation& operator=(const TripletPermutation&) = delete;
  TripletPermutation(TripletPermutation&&) = delete;
  TripletPermutation& operator=(TripletPermutation&&) = delete;

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool.swap(pool);
    clearStorage();
  }

  /// \a forEach(count, body) runs body(i) for i in [0, count), serially or
  /// in parallel, over chunks of \a chunkSize tracklets (HostAlgorithms.h):
  /// the same order either way.
  template <typename ForEach>
  [[nodiscard]] StageStatus build(const IterationParameters<Det>& parameters, const Sources& sources, const Tracklets& tracklets, TrackletPathBounds paths, const ForEach& forEach, Index chunkSize)
  {
    return guardOutOfMemory([&] {
      const auto view = tracklets.view();
      const auto offsets = tracklets.offsets();
      // The compaction reads its predicate twice (count, place): one byte per
      // tracklet, evaluated once, instead of the gathers twice.
      resizeBoundedVector(mEligible, tracklets.size());
      forEach(tracklets.size(), [&](Index first) { mEligible[first] = firstTrackletEligible(parameters, sources, view, offsets.data(), paths, first) ? 1U : 0U; });
      host_algorithms::compactIndices(forEach, tracklets.size(), chunkSize, [&](Index first) { return mEligible[first] != 0U; }, mIndices);
      return StageStatus::kComplete;
    });
  }

  void clearStorage() noexcept { deepVectorClear(mMemoryPool.get(), mIndices, mEligible); }

  [[nodiscard]] std::span<const Index> indices() const { return mIndices; }
  [[nodiscard]] const Index* data() const { return mIndices.data(); }
  [[nodiscard]] Index size() const { return static_cast<Index>(mIndices.size()); }

 private:
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  BoundedVector<Index> mIndices;
  BoundedVector<std::uint8_t> mEligible;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRIPLETPERMUTATION_H_
