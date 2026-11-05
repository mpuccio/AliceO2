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

/// \file TripletEndIndex.h
/// \brief The triplets ending with each tracklet: the neighbour stage's
///        candidates.

#ifndef GTRK_TRACKING_TRIPLETENDINDEX_H_
#define GTRK_TRACKING_TRIPLETENDINDEX_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/HostAlgorithms.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/Triplet.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <span>

namespace o2::gtrk
{

/// Triplets `T1 = (a, b, c)` and `T2 = (b, c, d)` are neighbours when
/// `T1.second == T2.first`, so the neighbour stage looks up, for each `T2`,
/// the triplets whose second tracklet is `T2.first`.
///
/// A counting sort of the triplet indices by their second tracklet: offsets
/// per tracklet, and within each tracklet the triplet indices in arena order.
/// The key (second tracklet, triplet index) is unique, so any sort -- a
/// device radix sort included -- gives the same index.
///
/// On the host the sort runs in two levels without atomics:
/// the arena in chunks of
/// `kChunkTriplets`, the tracklets in buckets of `kBucketTracklets`. Each chunk
/// counts its triplets per bucket; a serial scan of the chunk x bucket table,
/// small, gives each bucket its start and each chunk its cursor in each
/// bucket (the table is chunks x buckets, so the chunks are at most
/// \a maxChunks, the caller's tasks, and at least `kChunkTriplets` long, so the
/// table stays small); each chunk places its triplets at its cursors, so a bucket holds
/// its triplets in chunk order, which is arena order; then each bucket sorts
/// its own by their second tracklet, a counting sort over its tracklets alone,
/// into the offsets and the index.
class TripletEndIndex
{
 public:
  using Index = typename StageOutput<TripletT>::Index;
  static_assert(sizeof(Index) == 4U);
  static constexpr Index kChunkTriplets = 4096U;   ///< triplets per chunk of the first level
  static constexpr Index kBucketTracklets = 4096U; ///< tracklets per bucket: a bucket's counts fit a thread's cache

  TripletEndIndex() = default;
  ~TripletEndIndex() = default;
  TripletEndIndex(const TripletEndIndex&) = delete;
  TripletEndIndex& operator=(const TripletEndIndex&) = delete;
  TripletEndIndex(TripletEndIndex&&) = delete;
  TripletEndIndex& operator=(TripletEndIndex&&) = delete;

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool.swap(pool);
    clearStorage();
  }

  /// Index \a triplets, whose tracklet indices are below \a trackletCount, the
  /// chunks and the buckets through \a forEach(count, body): body(i) once for
  /// every i in [0, count), one task per i; at most \a maxChunks chunks,
  /// \a maxChunks > 0. The index is the same for any chunking.
  template <typename ForEach>
  [[nodiscard]] StageStatus build(Index trackletCount, const StageOutput<TripletT>& triplets, const ForEach& forEach, Index maxChunks)
  {
    return guardOutOfMemory([&] {
      const auto view = triplets.view();
      const Index count = triplets.size();
      const Index buckets = 1U + (trackletCount / kBucketTracklets);
      const Index perChunk = std::max(kChunkTriplets, count / maxChunks + 1U);
      const host_algorithms::Chunking chunking{count, perChunk};
      const Index chunks = chunking.count();
      // Kept from slice to slice and written whole before they are read.
      resizeBoundedVector(mOffsets, static_cast<std::size_t>(trackletCount) + 1U, mMemoryPool.get());
      resizeBoundedVector(mTriplets, count, mMemoryPool.get());
      resizeBoundedVector(mByBucket, count, mMemoryPool.get());
      resizeBoundedVector(mCursors, static_cast<std::size_t>(chunks) * buckets, mMemoryPool.get());
      resizeBoundedVector(mBucketStarts, static_cast<std::size_t>(buckets) + 1U, mMemoryPool.get());
      // 1. Each chunk's triplets per bucket.
      forEach(chunks, [&](Index chunk) {
        Index* row = mCursors.data() + (static_cast<std::size_t>(chunk) * buckets);
        std::fill_n(row, buckets, Index{0U});
        const Index end = chunking.end(chunk);
        for (Index triplet = chunking.begin(chunk); triplet < end; ++triplet) {
          ++row[view.second[triplet] / kBucketTracklets];
        }
      });
      // 2. Where each bucket begins, and each chunk's cursor in it: the
      // exclusive sums in (bucket, chunk) order.
      Index total = 0U;
      for (Index bucket = 0U; bucket < buckets; ++bucket) {
        mBucketStarts[bucket] = total;
        for (Index chunk = 0U; chunk < chunks; ++chunk) {
          Index& cell = mCursors[(static_cast<std::size_t>(chunk) * buckets) + bucket];
          const Index counted = cell;
          cell = total;
          total += counted;
        }
      }
      mBucketStarts[buckets] = total;
      // 3. Each chunk places its triplets at its cursors, in arena order: a
      // bucket's triplets end in chunk order, so in arena order.
      forEach(chunks, [&](Index chunk) {
        Index* cursor = mCursors.data() + (static_cast<std::size_t>(chunk) * buckets);
        const Index end = chunking.end(chunk);
        for (Index triplet = chunking.begin(chunk); triplet < end; ++triplet) {
          mByBucket[cursor[view.second[triplet] / kBucketTracklets]++] = triplet;
        }
      });
      // 4. Each bucket by second tracklet: the counts, then each tracklet's
      // offset from the bucket's start, then the placement behind a cursor per
      // tracklet, the bucket's triplets taken in arena order.
      forEach(buckets, [&](Index bucket) {
        const Index firstTracklet = bucket * kBucketTracklets;
        const Index endTracklet = std::min(trackletCount, firstTracklet + kBucketTracklets);
        const Index begin = mBucketStarts[bucket];
        const Index end = mBucketStarts[bucket + 1U];
        std::array<Index, kBucketTracklets + 1U> local{};
        for (Index at = begin; at < end; ++at) {
          ++local[(view.second[mByBucket[at]] - firstTracklet) + 1U];
        }
        local[0] = begin;
        for (Index tracklet = firstTracklet; tracklet < endTracklet; ++tracklet) {
          local[(tracklet - firstTracklet) + 1U] += local[tracklet - firstTracklet];
          mOffsets[tracklet] = local[tracklet - firstTracklet];
        }
        for (Index at = begin; at < end; ++at) {
          const Index triplet = mByBucket[at];
          mTriplets[local[view.second[triplet] - firstTracklet]++] = triplet;
        }
      });
      mOffsets[trackletCount] = count;
      return StageStatus::kComplete;
    });
  }

  void clearStorage() noexcept { deepVectorClear(mMemoryPool.get(), mOffsets, mTriplets, mByBucket, mCursors, mBucketStarts); }

  /// The raw columns the device stage reads: `trackletCount + 1` offsets and
  /// the triplet indices they delimit.
  [[nodiscard]] const Index* offsets() const { return mOffsets.data(); }
  [[nodiscard]] const Index* triplets() const { return mTriplets.data(); }

 private:
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  BoundedVector<Index> mOffsets;
  BoundedVector<Index> mTriplets;
  BoundedVector<Index> mByBucket;     ///< the triplets by bucket, in arena order within one (the first level)
  BoundedVector<Index> mCursors;      ///< chunks x buckets: the counts, then the cursors
  BoundedVector<Index> mBucketStarts; ///< buckets + 1
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRIPLETENDINDEX_H_
