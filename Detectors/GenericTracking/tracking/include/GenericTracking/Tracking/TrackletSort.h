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

/// \file TrackletSort.h
/// \brief Each source's tracklets in (link, lambda, target) order.
///
/// The triplet stage searches the tracklets leaving a
/// middle cluster by link and then by lambda: a binary search for the start of
/// its lambda window and a scan to its end. This puts every source's range of
/// the tracklet arena in that order, in place: on the host by
/// `std::ranges::sort` over a copy of the range's records
/// (`sortTrackletSegmentOnHost`), on a device by the library's segmented sort
/// on a stable pass by target, then by (link, lambda) (`trackletLinkLambdaKey`,
/// SliceGlue.cuh).

#ifndef GTRK_TRACKING_TRACKLETSORT_H_
#define GTRK_TRACKING_TRACKLETSORT_H_

#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonDef.h"
#include "GPUCommonMath.h"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace o2::gtrk
{

/// Whether a tracklet with key (\a link, \a lambda, \a target) sorts before one
/// with the other key.
///
/// Within one source the target fixes the target layer and so the link, and a
/// query visits each row once, so no two tracklets of a source share a target.
/// The key is therefore unique within a range: the sorted order does not
/// depend on the algorithm or the backend, and equal lambdas -- -0 and +0
/// among them -- are ordered by the target.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool trackletBefore(std::uint16_t link, Scalar lambda, std::uint32_t target, std::uint16_t otherLink,
                                                             Scalar otherLambda, std::uint32_t otherTarget)
{
  if (link != otherLink) {
    return link < otherLink;
  }
  if (lambda < otherLambda) {
    return true;
  }
  if (otherLambda < lambda) {
    return false;
  }
  return target < otherTarget;
}

/// The counters build's tally of one source's range of \a length tracklets.
GPUhdi() GPUAlwaysInline() void countSortSegment(std::uint32_t length)
{
  if (length < 2U) {
    return;
  }
  countWork(WorkCounter::kSortSegments);
  countWork(WorkCounter::kSortTracklets, length);
  if (length > 8U) {
    countWork(WorkCounter::kSortSegmentsOver8);
  }
  if (length > 32U) {
    countWork(WorkCounter::kSortSegmentsOver32);
    countWork(WorkCounter::kSortTrackletsOver32, length);
  }
  if (length > 128U) {
    countWork(WorkCounter::kSortSegmentsOver128);
  }
  if (length > 512U) {
    countWork(WorkCounter::kSortSegmentsOver512);
  }
}

/// \a lambda as an unsigned key in `trackletBefore`'s order of lambdas: the
/// float's bits made monotonic (every bit flipped for a negative, the sign bit
/// set otherwise). Zero and the denormals of both signs are first made +0: `<`
/// does not tell -0 from +0, and the trackers compare with denormals as zero --
/// the host's processes start with MXCSR's DAZ and FTZ set, the devices compile
/// with -ftz=true (CUDA) and -fgpu-flush-denormals-to-zero (HIP). A lambda is
/// never NaN (an atan2 of finite chords).
[[nodiscard]] GPUdi() std::uint32_t trackletLambdaKey(Scalar lambda)
{
  std::uint32_t bits = o2::gpu::CAMath::Float2UIntReint(lambda);
  constexpr std::uint32_t kSign = 0x80000000U;
  constexpr std::uint32_t kSmallestNormal = 0x00800000U;
  if ((bits & ~kSign) < kSmallestNormal) {
    bits = 0U;
  }
  return (bits & kSign) != 0U ? ~bits : (bits | kSign);
}

/// The link and lambda of a tracklet as one key in `trackletBefore`'s order,
/// the target aside: the device sorts on it after a stable pass by target.
[[nodiscard]] GPUdi() std::uint64_t trackletLinkLambdaKey(std::uint16_t link, Scalar lambda)
{
  return (std::uint64_t{link} << 32U) | std::uint64_t{trackletLambdaKey(lambda)};
}

/// Sort the tracklets in [\a begin, \a end) of \a tracklets in place, on the
/// host: `std::ranges::sort` over the
/// segment's records, copied into a buffer of the thread and written back.
/// Sorting the rows in place across the five columns is slower, by how much
/// depending on where the columns lie relative to each other;
/// the buffer holds the largest segment.
///
/// The key is unique within a source, so this is the order of
/// `trackletBefore`, for any length and any ties, in O(n log n); the library
/// picks its own way for short ranges and long ones.
inline void sortTrackletSegmentOnHost(const TrackletView& tracklets, std::uint32_t begin, std::uint32_t end)
{
  if constexpr (kWorkCounters) {
    countSortSegment(end - begin);
  }
  // Fewer than two rows are in order: no call into the library for them, the
  // most common case.
  if (end - begin < 2U) {
    return;
  }
  thread_local std::vector<Tracklet> records;
  const std::uint32_t count = end - begin;
  records.resize(count);
  for (std::uint32_t i = 0; i < count; ++i) {
    const auto at = begin + i;
    records[i] = Tracklet{
      .source = tracklets.source[at], .target = tracklets.target[at], .link = tracklets.link[at], .lambda = tracklets.lambda[at], .phi = tracklets.phi[at]};
  }
  std::ranges::sort(records, [](const Tracklet& a, const Tracklet& b) { return trackletBefore(a.link, a.lambda, a.target, b.link, b.lambda, b.target); });
  for (std::uint32_t i = 0; i < count; ++i) {
    const auto at = begin + i;
    tracklets.source[at] = records[i].source;
    tracklets.target[at] = records[i].target;
    tracklets.link[at] = records[i].link;
    tracklets.lambda[at] = records[i].lambda;
    tracklets.phi[at] = records[i].phi;
  }
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKLETSORT_H_
