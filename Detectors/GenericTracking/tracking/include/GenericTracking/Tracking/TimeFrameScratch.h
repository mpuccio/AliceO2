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

/// \file TimeFrameScratch.h
/// \brief The host tracker's stage storage, reused across its slices and iterations.
///
/// Storage only: the host stages run on CpuBackend.h; the device tracker
/// keeps its own (DeviceTracker.cuh). For the tracklet stage it
/// holds the numbering of the slice's source rows and the stage's canonical
/// output; then the triplets and their neighbours.
/// The used-cluster mask lives for the timeframe and is deliberately not here.

#ifndef GTRK_TRACKING_TIMEFRAMESCRATCH_H_
#define GTRK_TRACKING_TIMEFRAMESCRATCH_H_

#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Tracking/LayerRowMap.h"
#include "GenericTracking/Tracking/MeasurementRange.h"
#include "GenericTracking/Tracking/Neighbour.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/Triplet.h"
#include "GenericTracking/Tracking/TrackletPaths.h"
#include "GenericTracking/Tracking/TripletEndIndex.h"
#include "GenericTracking/Tracking/TripletPaths.h"
#include "GenericTracking/Tracking/TripletPermutation.h"

#include <array>
#include <cstddef>
#include <memory>
#include <span>

namespace o2::gtrk
{

template <DetectorLayout Det>
class TimeFrameScratch
{
 public:
  using Sources = LayerRowMap<Det>;
  using Tracklets = StageOutput<TrackletT>;
  using Triplets = StageOutput<TripletT>;
  using Neighbours = StageOutput<NeighbourT>;
  using Permutation = TripletPermutation<Det>;
  static constexpr std::size_t kLayerCount = Sources::kLayerCount;

  /// The layers a tracklet can start on: those with an outgoing link.
  static constexpr std::array<bool, kLayerCount> kSourceLayers = [] {
    std::array<bool, kLayerCount> presented{};
    for (std::size_t layer = 0; layer < kLayerCount; ++layer) {
      presented[layer] = outLinks<Det>.count(static_cast<LayerID>(layer)) != 0U;
    }
    return presented;
  }();

  TimeFrameScratch() = default;
  ~TimeFrameScratch() = default;
  TimeFrameScratch(const TimeFrameScratch&) = delete;
  TimeFrameScratch& operator=(const TimeFrameScratch&) = delete;
  TimeFrameScratch(TimeFrameScratch&&) = delete;
  TimeFrameScratch& operator=(TimeFrameScratch&&) = delete;

  void setMemoryPool(const std::shared_ptr<BoundedMemoryResource>& pool)
  {
    mTracklets.setMemoryPool(pool);
    mTriplets.setMemoryPool(pool);
    mTrackletPaths.setMemoryPool(pool);
    mPermutation.setMemoryPool(pool);
    mEnding.setMemoryPool(pool);
    mTripletPaths.setMemoryPool(pool);
    mNeighbours.setMemoryPool(pool);
  }

  /// Build the tracklet source map: every source layer's halo rows. Chains
  /// the output depends on reach into the halo on every layer; the tracklet
  /// stage's gate skips the rows beyond them.
  [[nodiscard]] static Sources makeSources(std::span<const MeasurementRange, kLayerCount> halo) { return Sources{halo, kSourceLayers}; }

  /// Take the slice's numbered source rows, `makeSources` of its halo, and
  /// size the tracklet output for them.
  void beginSlice(const Sources& sources)
  {
    mSources = sources;
    mTracklets.prepare(mSources.size());
  }

  /// Release every allocation back to the pool.
  void clearStorage() noexcept
  {
    mTracklets.clearStorage();
    mTriplets.clearStorage();
    mTrackletPaths.clearStorage();
    mPermutation.clearStorage();
    mEnding.clearStorage();
    mTripletPaths.clearStorage();
    mNeighbours.clearStorage();
    mSources = Sources{};
  }

  [[nodiscard]] const Sources& sources() const { return mSources; }
  [[nodiscard]] Tracklets& tracklets() { return mTracklets; }
  [[nodiscard]] const Tracklets& tracklets() const { return mTracklets; }
  /// The longest tracklet paths through each measurement, which the
  /// permutation and the triplet stage read.
  [[nodiscard]] TrackletPaths<Det>& trackletPaths() { return mTrackletPaths; }
  /// The first tracklets the triplet stage works on, in its order.
  [[nodiscard]] Permutation& permutation() { return mPermutation; }
  [[nodiscard]] Triplets& triplets() { return mTriplets; }
  [[nodiscard]] const Triplets& triplets() const { return mTriplets; }
  /// The triplets ending with each tracklet: the neighbour stage's candidates.
  [[nodiscard]] TripletEndIndex& ending() { return mEnding; }
  /// The longest triplet paths through each tracklet, which the neighbour
  /// stage reads.
  [[nodiscard]] TripletPaths<Det>& tripletPaths() { return mTripletPaths; }
  [[nodiscard]] Neighbours& neighbours() { return mNeighbours; }
  [[nodiscard]] const Neighbours& neighbours() const { return mNeighbours; }

 private:
  Sources mSources;
  Tracklets mTracklets;
  TrackletPaths<Det> mTrackletPaths;
  Permutation mPermutation;
  Triplets mTriplets;
  TripletEndIndex mEnding;
  TripletPaths<Det> mTripletPaths;
  Neighbours mNeighbours;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TIMEFRAMESCRATCH_H_
