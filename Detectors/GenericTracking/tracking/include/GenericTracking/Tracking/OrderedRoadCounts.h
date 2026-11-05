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

/// \file OrderedRoadCounts.h
/// \brief The road counts the neighbour stage writes in layer order.

#ifndef GTRK_TRACKING_ORDEREDROADCOUNTS_H_
#define GTRK_TRACKING_ORDEREDROADCOUNTS_H_

#include "GenericTracking/Tracking/MeasurementRange.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/Triplet.h"

#include "GPUCommonMath.h"

#include <cstdint>
#include <limits>

namespace o2::gtrk
{

/// Whether \a triplet's innermost row is the slice's own: in the core of its
/// first tracklet's source layer.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool roadOwned(TrackletConstView tracklets, TripletConstView triplets, const LayerID* linkSources,
                                                        const MeasurementRange* core, std::uint32_t triplet)
{
  const auto first = triplets.first[triplet];
  const auto row = tracklets.source[first];
  const auto range = core[linkSources[tracklets.link[first]]];
  return row >= range.begin && row < range.end;
}

/// A road \a sum as a 32-bit count: saturated, and \a overflow raised,
/// beyond it -- roads are never truncated silently.
[[nodiscard]] GPUdi() std::uint32_t saturatedRoads(std::uint64_t sum, std::uint32_t* overflow)
{
  constexpr std::uint64_t kLargest = std::numeric_limits<std::uint32_t>::max();
  if (sum > kLargest) {
    o2::gpu::CAMath::AtomicMax(overflow, std::uint32_t{1U});
    return static_cast<std::uint32_t>(kLargest);
  }
  return static_cast<std::uint32_t>(sum);
}

/// Count roads while generating neighbours in layer order. Predecessor counts
/// are final before this element runs. Flags and overflow are initialized once
/// before the pass; replaying it overwrites counts and repeats idempotent marks.
/// A default instance (no `roads`) counts nothing.
struct OrderedRoadCounts {
  std::uint32_t* roads = nullptr;
  std::uint32_t* successorAbove = nullptr;
  std::uint32_t* overflow = nullptr;
  const LayerID* linkSources = nullptr;
  const MeasurementRange* core = nullptr;

  GPUdi() void base(std::uint32_t element, TrackletConstView tracklets, TripletConstView triplets) const
  {
    roads[element] = roadOwned(tracklets, triplets, linkSources, core, element) ? 1U : 0U;
  }

  GPUdi() void add(std::uint32_t predecessor, std::uint64_t& sum) const
  {
    sum += roads[predecessor];
    o2::gpu::CAMath::AtomicMax(successorAbove + predecessor, std::uint32_t{1U});
  }

  GPUdi() void finish(std::uint32_t element, std::uint64_t sum) const { roads[element] = saturatedRoads(sum, overflow); }
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_ORDEREDROADCOUNTS_H_
