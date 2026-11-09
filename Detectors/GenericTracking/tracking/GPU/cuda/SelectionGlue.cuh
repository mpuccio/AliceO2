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

/// \file SelectionGlue.cuh
/// \brief An iteration's global selection on a device: GlobalSelection::select
///        over device storage.
///
/// Only storage, sorts, scans, key packing and orchestration are here; the
/// selection's logic is Selection.h's -- sortedKey, compareKeys, sameRecord,
/// BidReset, Bid, Decide, Prune -- called as it is:
///
///   keys        each record's sorted cluster set, word-major, kNoCluster
///               remapped to the row count so a word fits bit_width(rows)
///   dedup order K stable radix passes from the last word to the first,
///               starting from the identity: (key, record index) order
///   runs        heads where the key changes; each position's head by an
///               inclusive max-scan; duplicates and mismatches counted
///   quality     one stable radix pass on `selection::qualityKey` over the
///               unique records, in key order: ranked()
///   rounds      Prune, then BidReset, Bid, Decide, Prune while any record
///               is live, one flag read per pass
///   selected    the accepted ranks, as record indices, best first
///
/// Outputs go to the base end in the caller's scope, temporaries to the
/// scratch end. Included by the .cu that instantiates a tracker, and by the
/// test kernels, in the backend's namespace.

#ifndef GTRK_TRACKING_GPU_SELECTIONGLUE_CUH_
#define GTRK_TRACKING_GPU_SELECTIONGLUE_CUH_

#include "DeviceAlgorithms.cuh"
#include "Platform.cuh"

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/Selection.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::selection_glue
{

using Index = std::uint32_t;

/// What the selection read back and where it left its result.
struct Selected {
  DeviceSpan<Index> ranked;   ///< the unique records, best first, as record indices
  DeviceSpan<Index> selected; ///< the accepted records, best first, as record indices
  Index unique = 0U;
  Index duplicates = 0U;
  Index mismatches = 0U;
  Index rounds = 0U;
  Index tracks = 0U; ///< entries of `selected`
  std::uint32_t readbacks = 0U;
  DeviceSpan<Index> slots; ///< the rounds' `slotAll`, one per row, for the extension's bids to reuse
};

// Unnamed, so that two sources of one library may include this header.
namespace detail
{
namespace
{

/// Word j of record r's sorted cluster set at keys[j * records + r].
template <std::size_t K>
struct WriteKeys {
  const ClusterList<K>* clusters = nullptr;
  const std::uint8_t* count = nullptr;
  Index records = 0U;
  Index noCluster = 0U; ///< the row count: above every row
  Index* keys = nullptr;

  GPUdi() void operator()(Index record) const
  {
    const auto key = selection::sortedKey(clusters[record], count[record]);
    for (std::size_t word = 0; word < K; ++word) {
      const Index value = key.index[word];
      keys[(word * records) + record] = value == ClusterList<K>::kNoCluster ? noCluster : value;
    }
  }
};

struct Iota {
  Index* out = nullptr;
  GPUdi() void operator()(Index i) const { out[i] = i; }
};

/// A run's head where the key differs from the previous position's; the
/// position itself there, 0 elsewhere, for the max-scan.
template <std::size_t K>
struct MarkHeads {
  const Index* keys = nullptr; ///< word-major
  const Index* permutation = nullptr;
  Index records = 0U;
  std::uint8_t* head = nullptr;
  Index* headPosition = nullptr;

  GPUdi() void operator()(Index position) const
  {
    bool differs = position == 0U;
    if (!differs) {
      const Index record = permutation[position];
      const Index previous = permutation[position - 1U];
      for (std::size_t word = 0; word < K && !differs; ++word) {
        differs = keys[(word * records) + record] != keys[(word * records) + previous];
      }
    }
    head[position] = differs ? 1U : 0U;
    headPosition[position] = differs ? position : 0U;
  }
};

struct MaxOf {
  GPUdi() Index operator()(Index lhs, Index rhs) const { return lhs < rhs ? rhs : lhs; }
};

/// counts[0] += 1 for every duplicate, counts[1] += 1 for every one that is
/// not the record it duplicates: GlobalSelection::select's counts.
template <std::size_t K>
struct CountDuplicates {
  const ClusterList<K>* clusters = nullptr;
  const std::uint8_t* count = nullptr;
  const std::uint8_t* kinks = nullptr;
  const Scalar* chi2 = nullptr;
  const Index* permutation = nullptr;
  const std::uint8_t* head = nullptr;
  const Index* runHead = nullptr;
  Index* counts = nullptr;

  GPUdi() void operator()(Index position) const
  {
    if (head[position] != 0U) {
      return;
    }
    const Index record = permutation[position];
    const Index kept = permutation[runHead[position]];
    o2::gpu::CAMath::AtomicAdd(counts, Index{1U});
    if (!selection::sameRecord(clusters[record], count[record], kinks[record], chi2[record], clusters[kept], count[kept], kinks[kept], chi2[kept])) {
      o2::gpu::CAMath::AtomicAdd(counts + 1, Index{1U});
    }
  }
};

struct IsHead {
  const std::uint8_t* head = nullptr;
  GPUdi() bool operator()(Index position) const { return head[position] != 0U; }
};

/// The quality order's keys (`selection::qualityKey`); the order among
/// equal keys, the cluster set's, comes from the sort's stability.
struct QualityKeys {
  const std::uint8_t* count = nullptr;
  const std::uint8_t* kinks = nullptr;
  const Scalar* chi2 = nullptr;
  const Index* unique = nullptr;
  std::uint64_t* keys = nullptr;

  GPUdi() void operator()(Index i) const
  {
    const Index record = unique[i];
    keys[i] = selection::qualityKey(count[record], kinks[record], chi2[record]);
  }
};

struct IsAccepted {
  const std::uint8_t* accepted = nullptr;
  GPUdi() bool operator()(Index rank) const { return accepted[rank] != 0U; }
};

} // namespace
} // namespace detail

/// GlobalSelection::select over \a records (\a recordCount of them) into \a
/// mask (\a rows entries): the same ranked and selected records, duplicates,
/// mismatches and rounds. A record's \a sharedInnermost innermost clusters
/// are shareable.
template <std::size_t K, typename RecordView>
[[nodiscard]] StageStatus select(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, const RecordView& records, Index recordCount, std::uint8_t* mask, Index rows, std::uint32_t sharedInnermost, Selected& result)
{
  result = {};
  algorithms::detail::requireItems(recordCount);
  if (recordCount == 0U) {
    if (!result.ranked.allocate(arena, 0U) || !result.selected.allocate(arena, 0U)) {
      return StageStatus::kOutOfMemory;
    }
    return StageStatus::kComplete;
  }
  const int keyEndBit = algorithms::keyBits(rows);

  // The unique records in key order, their duplicates counted.
  DeviceSpan<Index> unique;
  {
    auto scratch = arena.scope(ArenaEnd::kScratch);
    DeviceSpan<Index> keys;
    DeviceSpan<Index> word;
    DeviceSpan<Index> sortedWord;
    DeviceSpan<Index> permutation;
    DeviceSpan<Index> nextPermutation;
    if (!keys.allocate(arena, K * std::size_t{recordCount}, ArenaEnd::kScratch) || !word.allocate(arena, recordCount, ArenaEnd::kScratch) || !sortedWord.allocate(arena, recordCount, ArenaEnd::kScratch) || !permutation.allocate(arena, recordCount, ArenaEnd::kScratch) || !nextPermutation.allocate(arena, recordCount, ArenaEnd::kScratch)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(runtime, stream, recordCount, detail::WriteKeys<K>{.clusters = records.clusters, .count = records.count, .records = recordCount, .noCluster = rows, .keys = keys.data()});
    algorithms::forEach(runtime, stream, recordCount, detail::Iota{.out = permutation.data()});
    Index* current = permutation.data();
    Index* next = nextPermutation.data();
    for (std::size_t w = K; w-- > 0U;) {
      algorithms::forEach(runtime, stream, recordCount, algorithms::Gather<Index>{.from = std::as_const(keys).data() + (w * recordCount), .indices = current, .to = word.data()});
      const StageStatus sorted = algorithms::sortPairs(runtime, arena, stream, std::as_const(word).data(), sortedWord.data(), current, next, recordCount, keyEndBit);
      if (sorted != StageStatus::kComplete) {
        return sorted;
      }
      std::swap(current, next);
    }

    DeviceSpan<std::uint8_t> head;
    DeviceSpan<Index> headPosition;
    DeviceSpan<Index> runHead;
    DeviceSpan<Index> counts;
    if (!head.allocate(arena, recordCount, ArenaEnd::kScratch) || !headPosition.allocate(arena, recordCount, ArenaEnd::kScratch) || !runHead.allocate(arena, recordCount, ArenaEnd::kScratch) || !counts.allocate(arena, 2U, ArenaEnd::kScratch)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(runtime, stream, recordCount,
                        detail::MarkHeads<K>{
                          .keys = std::as_const(keys).data(),
                          .permutation = current,
                          .records = recordCount,
                          .head = head.data(),
                          .headPosition = headPosition.data(),
                        });
    const StageStatus scanned = algorithms::inclusiveScan(runtime, arena, stream, std::as_const(headPosition).data(), runHead.data(), recordCount, detail::MaxOf{});
    if (scanned != StageStatus::kComplete) {
      return scanned;
    }
    runtime.fill(counts.data(), 0U, 2U * sizeof(Index), stream);
    algorithms::forEach(runtime, stream, recordCount,
                        detail::CountDuplicates<K>{
                          .clusters = records.clusters,
                          .count = records.count,
                          .kinks = records.kinks,
                          .chi2 = records.chi2,
                          .permutation = current,
                          .head = std::as_const(head).data(),
                          .runHead = std::as_const(runHead).data(),
                          .counts = counts.data(),
                        });
    // Read back with the heads' count: selectIndices synchronises.
    std::array<Index, 2> counted{};
    download(counts, counted.data(), 2U, stream);
    ++result.readbacks;

    // The heads' positions on the scratch end, then their records: `unique`
    // on the base end.
    DeviceSpan<Index> positions;
    Index heads = 0U;
    const StageStatus selected = algorithms::selectIndices(runtime, arena, stream, recordCount, detail::IsHead{.head = std::as_const(head).data()}, positions, heads, ArenaEnd::kScratch);
    ++result.readbacks;
    if (selected != StageStatus::kComplete) {
      return selected;
    }
    result.duplicates = counted[0];
    result.mismatches = counted[1];
    result.unique = heads;
    if (!unique.allocate(arena, heads)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(runtime, stream, heads, algorithms::Gather<Index>{.from = current, .indices = std::as_const(positions).data(), .to = unique.data()});
  }

  // The quality order.
  const Index uniqueCount = result.unique;
  if (!result.ranked.allocate(arena, uniqueCount)) {
    return StageStatus::kOutOfMemory;
  }
  {
    auto scratch = arena.scope(ArenaEnd::kScratch);
    DeviceSpan<std::uint64_t> keys;
    DeviceSpan<std::uint64_t> sortedKeys;
    if (!keys.allocate(arena, uniqueCount, ArenaEnd::kScratch) || !sortedKeys.allocate(arena, uniqueCount, ArenaEnd::kScratch)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(runtime, stream, uniqueCount,
                        detail::QualityKeys{
                          .count = records.count,
                          .kinks = records.kinks,
                          .chi2 = records.chi2,
                          .unique = std::as_const(unique).data(),
                          .keys = keys.data(),
                        });
    const StageStatus sorted = algorithms::sortPairs(runtime, arena, stream, std::as_const(keys).data(), sortedKeys.data(), std::as_const(unique).data(), result.ranked.data(), uniqueCount, selection::kQualityKeyBits);
    if (sorted != StageStatus::kComplete) {
      return sorted;
    }
  }

  // The rank arrays and the rounds.
  DeviceSpan<ClusterList<K>> rankClusters;
  DeviceSpan<std::uint8_t> rankCount;
  DeviceSpan<Index> priority;
  DeviceSpan<std::uint8_t> live;
  DeviceSpan<std::uint8_t> accepted;
  DeviceSpan<Index> slotAll;
  DeviceSpan<Index> slotExclusive;
  if (!rankClusters.allocate(arena, uniqueCount) || !rankCount.allocate(arena, uniqueCount) || !priority.allocate(arena, uniqueCount) || !live.allocate(arena, uniqueCount) || !accepted.allocate(arena, uniqueCount) || !slotAll.allocate(arena, rows) || !slotExclusive.allocate(arena, rows)) {
    return StageStatus::kOutOfMemory;
  }
  algorithms::forEach(runtime, stream, uniqueCount,
                      selection::RankFill<K>{
                        .clusters = records.clusters,
                        .count = records.count,
                        .ranked = std::as_const(result.ranked).data(),
                        .unique = uniqueCount,
                        .rankClusters = rankClusters.data(),
                        .rankCount = rankCount.data(),
                        .priority = priority.data(),
                        .live = live.data(),
                        .accepted = accepted.data(),
                      });
  runtime.fill(slotAll.data(), 0U, std::size_t{rows} * sizeof(Index), stream);
  runtime.fill(slotExclusive.data(), 0U, std::size_t{rows} * sizeof(Index), stream);
  const selection::BidReset<K> reset{rankClusters.data(), rankCount.data(), live.data(), slotAll.data(), slotExclusive.data()};
  const selection::Bid<K> bid{rankClusters.data(), rankCount.data(), live.data(), priority.data(), slotAll.data(), slotExclusive.data(), sharedInnermost};
  const selection::Decide<K> decide{rankClusters.data(), rankCount.data(), live.data(), accepted.data(), priority.data(), slotAll.data(), slotExclusive.data(), mask, sharedInnermost};
  // The first pass prunes the records an earlier iteration's tracks exclude;
  // each later one is a round. A pass that leaves a record live raises the
  // flag, so the raising passes are the host's rounds; each round accepts at
  // least the best live record, so there are at most `uniqueCount`.
  bool first = true;
  const auto rounds = algorithms::sweepToFixedPoint<1U>(runtime, arena, stream, uniqueCount, [&](Index* flags) {
    if (!first) {
      algorithms::forEach(runtime, stream, uniqueCount, reset);
      algorithms::forEach(runtime, stream, uniqueCount, bid);
      algorithms::forEach(runtime, stream, uniqueCount, decide);
    }
    first = false;
    algorithms::forEach(runtime, stream, uniqueCount, selection::Prune<K>{rankClusters.data(), rankCount.data(), live.data(), mask, sharedInnermost, flags});
  });
  if (rounds.status != StageStatus::kComplete) {
    return rounds.status;
  }
  result.rounds = rounds.changingSweeps;
  result.slots = slotAll;
  result.readbacks += rounds.changingSweeps + 1U;

  // The accepted ranks, best first, on the scratch end, then as record
  // indices.
  auto scratch = arena.scope(ArenaEnd::kScratch);
  DeviceSpan<Index> acceptedRanks;
  Index tracks = 0U;
  const StageStatus compacted = algorithms::selectIndices(runtime, arena, stream, uniqueCount, detail::IsAccepted{.accepted = std::as_const(accepted).data()}, acceptedRanks, tracks, ArenaEnd::kScratch);
  ++result.readbacks;
  if (compacted != StageStatus::kComplete) {
    return compacted;
  }
  if (!result.selected.allocate(arena, tracks)) {
    return StageStatus::kOutOfMemory;
  }
  algorithms::forEach(runtime, stream, tracks, algorithms::Gather<Index>{.from = std::as_const(result.ranked).data(), .indices = std::as_const(acceptedRanks).data(), .to = result.selected.data()});
  result.tracks = tracks;
  return StageStatus::kComplete;
}

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::selection_glue

#endif // GTRK_TRACKING_GPU_SELECTIONGLUE_CUH_
