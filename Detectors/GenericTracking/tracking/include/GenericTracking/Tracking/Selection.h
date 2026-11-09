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

/// \file Selection.h
/// \brief An iteration's global greedy selection, in a form the host and a
///        device run alike.
///
/// The records of every slice are concatenated, deduplicated by their
/// cluster set and ranked by the quality order -- sorts of fixed-width keys,
/// radix sorts on a device -- and selected by rounds of bidding, which end
/// with exactly the result of the sequential greedy walk. The records and the
/// used mask stay where they are; only the flag "is any record still live"
/// crosses to the controller per round.
///
/// The used mask holds, per measurement row, free, shared (taken only at
/// positions among a track's `k` innermost) or exclusive. Two records
/// conflict when they have a common row at a position beyond the first `k` in
/// either of them.

#ifndef GTRK_TRACKING_SELECTION_H_
#define GTRK_TRACKING_SELECTION_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/Columns.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Base/HostAlgorithms.h"
#include "GenericTracking/Tracking/Profile.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/StageOutput.h"

#include "GPUCommonMath.h"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>

namespace o2::gtrk
{

namespace selection
{

/// Whether a record's cluster at \a position is shareable.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool shareable(std::uint32_t position, std::uint32_t sharedInnermost) { return position < sharedInnermost; }

/// Whether a record conflicts with the tracks already in \a mask.
template <std::size_t K>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool conflicts(const ClusterList<K>& clusters, std::uint32_t count, const std::uint8_t* mask, std::uint32_t sharedInnermost)
{
  for (std::uint32_t position = 0; position < count; ++position) {
    const std::uint8_t state = mask[clusters.index[position]];
    if (state == kExclusive || (state == kShared && !shareable(position, sharedInnermost))) {
      return true;
    }
  }
  return false;
}

/// The cluster set, sorted: the deduplication key.
template <std::size_t K>
[[nodiscard]] GPUhdi() GPUAlwaysInline() ClusterList<K> sortedKey(const ClusterList<K>& clusters, std::uint32_t count)
{
  ClusterList<K> key = clusters;
  insertionSort(key.index, count, [](std::uint32_t a, std::uint32_t b) { return a < b; });
  return key;
}

template <std::size_t K>
[[nodiscard]] GPUhdi() GPUAlwaysInline() int compareKeys(const ClusterList<K>& lhs, const ClusterList<K>& rhs)
{
  for (std::size_t i = 0; i < K; ++i) {
    if (lhs.index[i] != rhs.index[i]) {
      return lhs.index[i] < rhs.index[i] ? -1 : 1;
    }
  }
  return 0;
}

/// The width of `qualityKey` in bits.
inline constexpr int kQualityKeyBits = 48;

/// The quality order's key: clusters descending, kinks ascending, chi2
/// ascending, over `kQualityKeyBits` bits; the cluster set breaks ties. The fits refuse a
/// negative or NaN chi2 (`updateUsable`); clamped here as well, so that a
/// record that slipped through orders first rather than by its sign bit.
[[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint64_t qualityKey(std::uint8_t count, std::uint8_t kinks, Scalar chi2)
{
  const Scalar quality = chi2 > 0.F ? chi2 : 0.F;
  return (std::uint64_t{255U - count} << 40U) | (std::uint64_t{kinks} << 32U) | std::uint64_t{std::bit_cast<std::uint32_t>(quality)};
}

/// Whether a duplicate is the record it duplicates -- the same count, kinks,
/// chi2 bits and clusters in their order: a record found by more than one
/// slice is the same record in each.
template <std::size_t K>
[[nodiscard]] GPUdi() bool sameRecord(const ClusterList<K>& clusters, std::uint8_t count, std::uint8_t kinks, Scalar chi2, const ClusterList<K>& otherClusters, std::uint8_t otherCount, std::uint8_t otherKinks, Scalar otherChi2)
{
  return count == otherCount && kinks == otherKinks && std::bit_cast<std::uint32_t>(chi2) == std::bit_cast<std::uint32_t>(otherChi2) && compareKeys(clusters, otherClusters) == 0;
}

/// One element per rank: the rank arrays the rounds read -- the ranked
/// record's clusters and count, its priority (the best the highest), live and
/// not accepted.
template <std::size_t K>
struct RankFill {
  const ClusterList<K>* clusters = nullptr;
  const std::uint8_t* count = nullptr;
  const std::uint32_t* ranked = nullptr;
  std::uint32_t unique = 0U;
  ClusterList<K>* rankClusters = nullptr;
  std::uint8_t* rankCount = nullptr;
  std::uint32_t* priority = nullptr;
  std::uint8_t* live = nullptr;
  std::uint8_t* accepted = nullptr;

  GPUhdi() void operator()(std::uint32_t rank) const
  {
    const std::uint32_t record = ranked[rank];
    rankClusters[rank] = clusters[record];
    rankCount[rank] = count[record];
    priority[rank] = unique - rank;
    live[rank] = 1U;
    accepted[rank] = 0U;
  }
};

/// One element per live record: clear its clusters' bid slots.
template <std::size_t K>
struct BidReset {
  const ClusterList<K>* clusters = nullptr;
  const std::uint8_t* count = nullptr;
  const std::uint8_t* live = nullptr;
  std::uint32_t* slotAll = nullptr;
  std::uint32_t* slotExclusive = nullptr;

  GPUhdi() void operator()(std::uint32_t record) const
  {
    if (live[record] == 0U) {
      return;
    }
    for (std::uint32_t position = 0; position < count[record]; ++position) {
      slotAll[clusters[record].index[position]] = 0U;
      slotExclusive[clusters[record].index[position]] = 0U;
    }
  }
};

/// One element per live record: bid its priority on every cluster, and at an
/// exclusive position also on the exclusive users' slot.
template <std::size_t K>
struct Bid {
  const ClusterList<K>* clusters = nullptr;
  const std::uint8_t* count = nullptr;
  const std::uint8_t* live = nullptr;
  const std::uint32_t* priority = nullptr;
  std::uint32_t* slotAll = nullptr;
  std::uint32_t* slotExclusive = nullptr;
  std::uint32_t sharedInnermost = 0U;

  GPUdi() void operator()(std::uint32_t record) const
  {
    if (live[record] == 0U) {
      return;
    }
    for (std::uint32_t position = 0; position < count[record]; ++position) {
      const std::uint32_t row = clusters[record].index[position];
      o2::gpu::CAMath::AtomicMax(slotAll + row, priority[record]);
      if (!shareable(position, sharedInnermost)) {
        o2::gpu::CAMath::AtomicMax(slotExclusive + row, priority[record]);
      }
    }
  }
};

/// One element per live record: accepted if it wins every cluster -- holds
/// the maximum at an exclusive position, outbids every exclusive user at a
/// shareable one -- and then writes the mask. Two records accepted in one
/// round never conflict, so their writes never disagree.
template <std::size_t K>
struct Decide {
  const ClusterList<K>* clusters = nullptr;
  const std::uint8_t* count = nullptr;
  std::uint8_t* live = nullptr;
  std::uint8_t* accepted = nullptr;
  const std::uint32_t* priority = nullptr;
  const std::uint32_t* slotAll = nullptr;
  const std::uint32_t* slotExclusive = nullptr;
  std::uint8_t* mask = nullptr;
  std::uint32_t sharedInnermost = 0U;

  GPUhdi() void operator()(std::uint32_t record) const
  {
    if (live[record] == 0U) {
      return;
    }
    const std::uint32_t own = priority[record];
    for (std::uint32_t position = 0; position < count[record]; ++position) {
      const std::uint32_t row = clusters[record].index[position];
      const bool wins = shareable(position, sharedInnermost) ? slotExclusive[row] < own : slotAll[row] == own;
      if (!wins) {
        return;
      }
    }
    accepted[record] = 1U;
    live[record] = 0U;
    // Two accepted records sharing an innermost row store the same byte.
    for (std::uint32_t position = 0; position < count[record]; ++position) {
      o2::gpu::CAMath::AtomicStore(mask + clusters[record].index[position], static_cast<std::uint8_t>(shareable(position, sharedInnermost) ? kShared : kExclusive));
    }
  }
};

/// One element per live record: out if it conflicts with the mask; \a stillLive
/// raised otherwise.
template <std::size_t K>
struct Prune {
  const ClusterList<K>* clusters = nullptr;
  const std::uint8_t* count = nullptr;
  std::uint8_t* live = nullptr;
  const std::uint8_t* mask = nullptr;
  std::uint32_t sharedInnermost = 0U;
  std::uint32_t* stillLive = nullptr;

  GPUdi() void operator()(std::uint32_t record) const
  {
    if (live[record] == 0U) {
      return;
    }
    if (conflicts(clusters[record], count[record], mask, sharedInnermost)) {
      live[record] = 0U;
      return;
    }
    o2::gpu::CAMath::AtomicMax(stillLive, std::uint32_t{1U});
  }
};

} // namespace selection

/// One iteration's global selection over the records of every slice.
template <std::size_t K>
class GlobalSelection
{
 public:
  using Records = StageOutput<SelectionRecordOf<K>::template RecordT>;
  using Index = std::uint32_t;

  GlobalSelection() = default;
  ~GlobalSelection() = default;
  GlobalSelection(const GlobalSelection&) = delete;
  GlobalSelection& operator=(const GlobalSelection&) = delete;
  GlobalSelection(GlobalSelection&&) = delete;
  GlobalSelection& operator=(GlobalSelection&&) = delete;

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool.swap(pool);
    clearStorage();
  }

  /// Start an iteration: no records.
  void clear() { mRecords.resize(0U); }

  /// Append one slice's records; \a forEach(count, body) runs body(i) for i
  /// in [0, count), serially or in parallel (the copy).
  template <typename ForEach>
  [[nodiscard]] StageStatus append(const Records& records, const ForEach& forEach)
  {
    return guardOutOfMemory([&] {
      const std::size_t begin = mRecords.size();
      const std::size_t count = records.size();
      mRecords.resize(begin + count);
      copyRows(records.view(), begin, count, forEach);
      return StageStatus::kComplete;
    });
  }

  /// Append every record of \a other: an iteration's records made by the
  /// previous iteration's slices (one graph for both rounds). Records whose rows an earlier track holds are out before the
  /// selection's first round, as any.
  template <typename ForEach>
  [[nodiscard]] StageStatus append(const GlobalSelection& other, const ForEach& forEach)
  {
    return guardOutOfMemory([&] {
      const std::size_t begin = mRecords.size();
      const std::size_t count = other.mRecords.size();
      mRecords.resize(begin + count);
      copyRows(other.mRecords.view(), begin, count, forEach);
      return StageStatus::kComplete;
    });
  }

  /// Select into \a mask, one entry per measurement row, which holds the
  /// earlier iterations' tracks. A record's `k` innermost clusters are
  /// shareable, \a sharedInnermost. \a forEach(elements, body) runs
  /// body(element) for every element, serially or in parallel, and
  /// \a sort(indices, less) sorts a span of indices by a strict total order,
  /// on the same threads (`CpuStage::sortIndices`).
  template <typename ForEach, typename Sort>
  [[nodiscard]] StageStatus select(ForEach forEach, Sort sort, std::uint32_t sharedInnermost, std::span<std::uint8_t> mask)
  {
    using profile::Phase;
    return guardOutOfMemory([&] {
      const auto records = mRecords.view();
      const auto total = static_cast<Index>(mRecords.size());
      mDuplicates = 0U;
      mMismatches = 0U;

      // Deduplication: the records sorted by their cluster set, the first of
      // each run kept; a record found by more than one slice is the same
      // record in each.
      {
        const profile::Zone zone{Phase::kSelectionKeys};
        // The arrays are kept from iteration to iteration and filled in
        // parallel, not reallocated.
        resizeBoundedVector(mKeys, total, mMemoryPool.get());
        forEach(total, [&](Index record) { mKeys[record] = selection::sortedKey(records.clusters[record], records.count[record]); });
      }
      {
        const profile::Zone zone{Phase::kSelectionOrder};
        resizeBoundedVector(mOrder, total, mMemoryPool.get());
        forEach(total, [&](Index record) { mOrder[record] = record; });
        // A total order -- the tie-break is the record index -- so the sorted
        // result is unique and a parallel unstable sort gives the serial one's.
        sort(std::span<Index>{mOrder}, [&](Index lhs, Index rhs) {
          const int order = selection::compareKeys(mKeys[lhs], mKeys[rhs]);
          if (order != 0) {
            return order < 0;
          }
          return lhs < rhs;
        });
      }
      std::optional<profile::Zone> zone{std::in_place, Phase::kSelectionDeduplicate};
      mUnique.resize(0U);
      for (Index i = 0; i < total; ++i) {
        const Index record = mOrder[i];
        if (i != 0U && selection::compareKeys(mKeys[record], mKeys[mOrder[i - 1U]]) == 0) {
          const Index kept = mUnique.back();
          ++mDuplicates;
          const bool same = selection::sameRecord(records.clusters[record], records.count[record], records.kinks[record], records.chi2[record], records.clusters[kept], records.count[kept], records.kinks[kept], records.chi2[kept]);
          mMismatches += same ? 0U : 1U;
          continue;
        }
        mUnique.push_back(record);
      }

      // The quality order: clusters descending, kinks ascending, chi2
      // ascending (`selection::qualityKey`, the device's sort key), the
      // cluster set ascending; the rank gives each record its unique priority.
      // Also a total order: the deduplication above left one record per
      // cluster set, so the final key comparison never returns 0.
      zone.reset();
      zone.emplace(Phase::kSelectionRank);
      sort(std::span<Index>{mUnique}, [&](Index lhs, Index rhs) {
        const auto lhsKey = selection::qualityKey(records.count[lhs], records.kinks[lhs], records.chi2[lhs]);
        const auto rhsKey = selection::qualityKey(records.count[rhs], records.kinks[rhs], records.chi2[rhs]);
        if (lhsKey != rhsKey) {
          return lhsKey < rhsKey;
        }
        return selection::compareKeys(mKeys[lhs], mKeys[rhs]) < 0;
      });
      const auto unique = static_cast<Index>(mUnique.size());
      zone.reset();
      zone.emplace(Phase::kSelectionFill);
      resizeBoundedVector(mClusters, unique, mMemoryPool.get());
      resizeBoundedVector(mCount, unique, mMemoryPool.get());
      resizeBoundedVector(mPriority, unique, mMemoryPool.get());
      resizeBoundedVector(mLive, unique, mMemoryPool.get());
      resizeBoundedVector(mAccepted, unique, mMemoryPool.get());
      forEach(unique, selection::RankFill<K>{
                        .clusters = records.clusters,
                        .count = records.count,
                        .ranked = mUnique.data(),
                        .unique = unique,
                        .rankClusters = mClusters.data(),
                        .rankCount = mCount.data(),
                        .priority = mPriority.data(),
                        .live = mLive.data(),
                        .accepted = mAccepted.data(),
                      });

      // The rounds; a record conflicting with an earlier iteration's track is
      // out before the first.
      zone.reset();
      zone.emplace(Phase::kSelectionSlots);
      // Kept from call to call and zeroed in parallel, not reallocated.
      mSlotAll.resize(mask.size());
      mSlotExclusive.resize(mask.size());
      forEach(static_cast<Index>(mask.size()), [&](Index row) {
        mSlotAll[row] = 0U;
        mSlotExclusive[row] = 0U;
      });
      zone.reset();
      zone.emplace(Phase::kSelectionRounds);
      Index stillLive = 0U;
      const selection::Prune<K> prune{mClusters.data(), mCount.data(), mLive.data(), mask.data(), sharedInnermost, &stillLive};
      const selection::BidReset<K> reset{mClusters.data(), mCount.data(), mLive.data(), mSlotAll.data(), mSlotExclusive.data()};
      const selection::Bid<K> bid{mClusters.data(), mCount.data(), mLive.data(), mPriority.data(), mSlotAll.data(), mSlotExclusive.data(), sharedInnermost};
      const selection::Decide<K> decide{mClusters.data(), mCount.data(), mLive.data(), mAccepted.data(), mPriority.data(), mSlotAll.data(), mSlotExclusive.data(), mask.data(), sharedInnermost};
      forEach(unique, prune);
      mRounds = 0U;
      while (stillLive != 0U) {
        forEach(unique, reset);
        forEach(unique, bid);
        forEach(unique, decide);
        stillLive = 0U;
        forEach(unique, prune);
        ++mRounds;
      }

      zone.reset();
      zone.emplace(Phase::kSelectionCompact);
      mSelected.resize(0U);
      for (Index rank = 0; rank < unique; ++rank) {
        if (mAccepted[rank] != 0U) {
          mSelected.push_back(mUnique[rank]);
        }
      }
      return StageStatus::kComplete;
    });
  }

  void clearStorage() noexcept
  {
    clearRecordStorage();
    deepVectorClear(mMemoryPool.get(), mSlotAll, mSlotExclusive);
  }

  /// What grows with an iteration's records back to the pool; the slots, one
  /// per row of the timeframe whatever the iteration, stay.
  void clearRecordStorage() noexcept
  {
    mRecords.setMemoryResource(mMemoryPool.get());
    deepVectorClear(mMemoryPool.get(), mKeys, mOrder, mUnique, mClusters, mCount, mPriority, mLive, mAccepted, mSelected);
  }

  /// The records appended this iteration.
  [[nodiscard]] typename BoundedColumns<SelectionRecordOf<K>::template RecordT>::ConstView records() const { return mRecords.view(); }
  [[nodiscard]] Index recordCount() const { return static_cast<Index>(mRecords.size()); }
  /// The accepted records, as indices into records(), best first.
  [[nodiscard]] std::span<const Index> selected() const { return mSelected; }
  /// The deduplicated records, best first: the order the selection follows.
  [[nodiscard]] std::span<const Index> ranked() const { return mUnique; }
  [[nodiscard]] Index duplicates() const { return mDuplicates; }
  /// Duplicates that differ from the record they duplicate: a record found by
  /// more than one slice is the same record in each, so this is zero, or the
  /// output is wrong.
  [[nodiscard]] Index mismatches() const { return mMismatches; }
  [[nodiscard]] Index rounds() const { return mRounds; }

 private:
  /// The \a count rows of \a from to rows \a begin on, every column of the
  /// record, whatever it holds, the rows in \a forEach's chunks.
  template <typename View, typename ForEach>
  void copyRows(const View& from, std::size_t begin, std::size_t count, const ForEach& forEach)
  {
    const auto to = mRecords.view();
    forEach(static_cast<Index>(count), [&](Index row) { from.apply([&](const auto*... source) { to.apply([&](auto*... target) { ((target[begin + row] = source[row]), ...); }); }); });
  }

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  BoundedColumns<SelectionRecordOf<K>::template RecordT> mRecords;
  BoundedVector<ClusterList<K>> mKeys;
  BoundedVector<Index> mOrder;
  BoundedVector<Index> mUnique;
  BoundedVector<ClusterList<K>> mClusters;
  BoundedVector<std::uint8_t> mCount;
  BoundedVector<Index> mPriority;
  BoundedVector<std::uint8_t> mLive;
  BoundedVector<std::uint8_t> mAccepted;
  BoundedVector<Index> mSlotAll;
  BoundedVector<Index> mSlotExclusive;
  BoundedVector<Index> mSelected;
  Index mDuplicates = 0U;
  Index mMismatches = 0U;
  Index mRounds = 0U;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SELECTION_H_
