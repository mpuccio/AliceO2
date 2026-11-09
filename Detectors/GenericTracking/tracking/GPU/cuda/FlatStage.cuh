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

/// \file FlatStage.cuh
/// \brief A stage run over a flat candidate list: the walk per element, the
///        test per candidate.
///
/// An element's candidates are what its walk admits, in the walk's order; a
/// stage's records are the candidates its test accepts, in that order, with
/// the element's count. `runFlat` runs them apart:
///
/// count     one thread an element: the walk, counting; the exclusive scan
///           of the counts is each element's pair offset in the list.
/// chunk     a range [a, b) of the list sized to the arena's free bytes;
///           the elements whose candidates lie in it are [sourceBegin,
///           sourceEnd) of the listed elements (one may straddle both ends).
/// fill      one thread an element of the chunk (a group of the adapter's
///           `kFillLanes` where it has more than one): the walk again, writing the
///           candidates in [a, b) as pairs at their offsets, and the
///           element's context the test reads.
/// test      `selectIndices` over the pairs with the adapter's test: the
///           accepted pairs' indices, stable, so in list order.
/// emit      one thread an accepted pair: its record into the chunk's
///           accepted columns; the elements' counts from the runs of equal
///           sources in the selected list.
/// assemble  the counts' scan is the offsets; one chunk's columns are the
///           output, several are gathered in chunk order.
///
/// The records and their order are the serial stage's: the walk is the
/// stage's own traversal, the test its own test, the stable select keeps the
/// walk's order, and chunks are consecutive ranges of it. No ordinal limit
/// on the candidates: a source may hold more than the output admits of
/// records (the output's limit per element, `kMaxPerElement`, is checked on
/// the accepted counts as the staging checks it).
///
/// The adapter (one per stage, copied to the kernels by value):
///   std::uint32_t count(element) const                       -- the walk, counting
///   void fill(element, position, Sink& sink) const           -- the walk, `sink.put(row, key)` per candidate; the
///                                                               element's context at `position - sourceBegin`
///   and, with `kFillLanes` above 1, fill(element, position, Sink&, group) -- the walk by the group's lanes in step,
///                                                               a step's candidates by `sink.putStep` in lane order
///   std::uint32_t sourceOf(const CandidatePair&) const       -- the listed position the pair's key holds
///   bool test(k, const CandidatePair&) const                 -- k the pair's index in the chunk, for what the test keeps for emit
///   Record emit(k, const CandidatePair&, element) const
///   bool beginChunk(DeviceArena&, sourceBegin, sourceCount, pairs) -- host: the chunk's contexts and per-pair storage, in the
///                                                               current scratch scope
///   std::size_t contextBytes(sources) const                  -- host: what beginChunk takes for that many sources
///   std::size_t pairBytes() const                            -- host: and per pair
///
/// \a listed may be nullptr: the elements are then 0 to listedCount - 1.
///
/// Included only by the .cu sources of a backend library or its test library,
/// as GpuStage.cuh.

#ifndef GTRK_TRACKING_GPU_FLATSTAGE_CUH_
#define GTRK_TRACKING_GPU_FLATSTAGE_CUH_

#include "DeviceAlgorithms.cuh"
#include "Platform.cuh"
#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStageOutput.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/StageOutput.h"

#include "GPUCommonAlgorithm.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::flat
{

/// One candidate: its row and the adapter's key (the listed position of its
/// element, with whatever the adapter packs beside it).
struct CandidatePair {
  std::uint32_t row = 0U;
  std::uint32_t key = 0U;
};

/// What the fill hands an element's walk: `put` every candidate in the walk's
/// order; those in the chunk are written, the others only counted past.
struct Sink {
  CandidatePair* pairs = nullptr;
  std::uint32_t next = 0U;       ///< the candidate's position in the list
  std::uint32_t chunkBegin = 0U; ///< the chunk [chunkBegin, chunkEnd)
  std::uint32_t chunkEnd = 0U;
  GPUdi() void put(std::uint32_t row, std::uint32_t key)
  {
    const std::uint32_t at = next++;
    if (at >= chunkBegin && at < chunkEnd) {
      pairs[at - chunkBegin] = CandidatePair{.row = row, .key = key};
    }
  }
  /// A group's step: \a mask the lanes with a candidate;
  /// this lane's, if \a has, after those of the lanes before it, and `next`
  /// past them all.
  GPUdi() void putStep(std::uint64_t mask, std::uint32_t lane, bool has, std::uint32_t row, std::uint32_t key)
  {
    if (has) {
      const std::uint32_t at = next + static_cast<std::uint32_t>(__popcll(mask & ((std::uint64_t{1} << lane) - 1U)));
      if (at >= chunkBegin && at < chunkEnd) {
        pairs[at - chunkBegin] = CandidatePair{.row = row, .key = key};
      }
    }
    next += static_cast<std::uint32_t>(__popcll(mask));
  }
};

/// A test hook: chunks of at most this many pairs (0: the arena's free bytes
/// decide); set only by the tests, directly or through
/// `DeviceTracker::setFlatChunkPairs`.
struct FlatOptions {
  std::size_t maxChunkPairs = 0U;
};

namespace detail
{

/// The element at a listed position: the list's, or the position itself.
GPUdi() std::uint32_t elementAt(const std::uint32_t* listed, std::uint32_t position) { return listed != nullptr ? listed[position] : position; }

/// The lanes an adapter's fill walks an element with: its `kFillLanes`, or one.
template <typename Adapter>
constexpr std::uint32_t fillLanes()
{
  if constexpr (requires { Adapter::kFillLanes; }) {
    return Adapter::kFillLanes;
  } else {
    return 1U;
  }
}

/// The candidates of each listed element, and the 0 after the last that
/// makes the scan's last offset the total.
template <typename Adapter>
struct CountElements {
  Adapter adapter;
  const std::uint32_t* listed;
  std::uint32_t listedCount;
  std::uint32_t* counts;
  GPUdi() void operator()(std::uint32_t i) const { counts[i] = i < listedCount ? adapter.count(elementAt(listed, i)) : 0U; }
};

/// A run's words on the device, one block: zeroed by one fill, the sum and
/// the largest read back by one copy.
struct FlatWords {
  Counter sum;            ///< the list's length in 64 bits (ReduceCounts)
  std::uint32_t largest;  ///< its longest element (ReduceCounts)
  std::uint32_t range[2]; ///< a chunk's elements (chunkRangeKernel)
  std::uint32_t overflow; ///< an element past the output's limit (AddCounts)
};

/// The chunk's elements: [0] the first listed position whose candidates reach
/// past \a a, [1] the first whose candidates start at or past \a b. A
/// template, as every kernel of a header: one instance a library.
template <typename = void>
GPUg() void chunkRangeKernel(const std::uint32_t* offsets, std::uint32_t listedCount, std::uint32_t a, std::uint32_t b, std::uint32_t* range)
{
  if (blockIdx.x != 0U || threadIdx.x != 0U) {
    return;
  }
  range[0] = o2::gpu::CAAlgo::partitionPoint(0U, listedCount, [&](std::uint32_t i) { return offsets[i + 1U] <= a; });
  range[1] = o2::gpu::CAAlgo::partitionPoint(0U, listedCount, [&](std::uint32_t i) { return offsets[i] < b; });
}

template <typename Adapter>
struct FillElements {
  Adapter adapter;
  const std::uint32_t* listed;
  const std::uint32_t* offsets;
  std::uint32_t sourceBegin;
  CandidatePair* pairs;
  std::uint32_t chunkBegin;
  std::uint32_t chunkEnd;
  GPUdi() void operator()(std::uint32_t j) const
  {
    const std::uint32_t position = sourceBegin + j;
    // A source without candidates: no pair of it, so no test or record reads its context.
    if (offsets[position] == offsets[position + 1U]) {
      return;
    }
    Sink sink{.pairs = pairs, .next = offsets[position], .chunkBegin = chunkBegin, .chunkEnd = chunkEnd};
    adapter.fill(elementAt(listed, position), position, sink);
  }
  /// The same by a group of the adapter's `kFillLanes` lanes.
  template <typename Group>
  GPUdi() void together(std::uint32_t j, const Group& group) const
  {
    const std::uint32_t position = sourceBegin + j;
    if (offsets[position] == offsets[position + 1U]) {
      return;
    }
    Sink sink{.pairs = pairs, .next = offsets[position], .chunkBegin = chunkBegin, .chunkEnd = chunkEnd};
    adapter.fill(elementAt(listed, position), position, sink, group);
  }
};

template <typename Adapter>
struct TestPairs {
  Adapter adapter;
  const CandidatePair* pairs;
  GPUdi() bool operator()(std::uint32_t k) const { return adapter.test(k, pairs[k]); }
};

/// One accepted pair's record into the chunk's columns, and the runs of equal
/// sources: where a source's run starts and ends in the selected list.
template <template <template <class> class> class S, typename Adapter>
struct EmitPairs {
  Adapter adapter;
  const std::uint32_t* listed;
  const CandidatePair* pairs;
  const std::uint32_t* selected;
  std::uint32_t selectedCount;
  std::uint32_t sourceBegin;
  typename DeviceColumns<S>::View accepted;
  std::uint32_t* runStart; ///< per chunk element, the first selected index of its run
  std::uint32_t* runEnd;   ///< one past the last
  GPUdi() void operator()(std::uint32_t k) const
  {
    const CandidatePair pair = pairs[selected[k]];
    const std::uint32_t source = adapter.sourceOf(pair);
    const auto record = adapter.emit(selected[k], pair, elementAt(listed, source));
    accepted.apply([&](auto*... column) { record.apply([&](const auto&... value) { ((column[k] = value), ...); }); });
    if (k == 0U || adapter.sourceOf(pairs[selected[k - 1U]]) != source) {
      runStart[source - sourceBegin] = k;
    }
    if (k + 1U == selectedCount || adapter.sourceOf(pairs[selected[k + 1U]]) != source) {
      runEnd[source - sourceBegin] = k + 1U;
    }
  }
};

/// The counts summed in 64 bits and their largest, into \a sum and
/// \a largest: whether the 32-bit offsets hold, and the list's longest element.
struct ReduceCounts {
  const std::uint32_t* counts;
  std::uint32_t count;
  Counter* sum;
  std::uint32_t* largest;
  GPUdi() void operator()(std::uint32_t i) const
  {
    // One thread of the grid-stride loop adds every `stride`-th count.
    const std::uint32_t stride = blockDim.x * gridDim.x;
    Counter mine = 0U;
    std::uint32_t top = 0U;
    for (std::uint32_t k = i; k < count; k += stride) {
      mine += counts[k];
      top = top > counts[k] ? top : counts[k];
    }
    if (mine != 0U) {
      atomicAdd(sum, mine);
    }
    if (top != 0U) {
      atomicMax(largest, top);
    }
  }
};

/// The chunk's elements' accepted counts added to the output's, and the
/// largest seen (a flag past the output's limit per element).
struct AddCounts {
  const std::uint32_t* listed;
  std::uint32_t sourceBegin;
  const std::uint32_t* runStart;
  const std::uint32_t* runEnd;
  std::uint32_t* counts;
  std::uint32_t* overflow;
  std::uint32_t limit;
  GPUdi() void operator()(std::uint32_t j) const
  {
    const std::uint32_t added = runEnd[j] - runStart[j];
    const std::uint32_t element = elementAt(listed, sourceBegin + j);
    const std::uint32_t total = counts[element] + added;
    counts[element] = total;
    if (total > limit) {
      *overflow = 1U;
    }
  }
};

} // namespace detail

/// The stage \a adapter describes over the \a listedCount listed elements,
/// into \a output (prepared for its element count by the caller; the others'
/// counts 0): the records of the serial stage in its order. The list and
/// the chunks' working storage on the scratch end within this call; the
/// output and the chunks' accepted records on the base end, in the caller's
/// scope: a run of several chunks holds its records twice until that scope
/// closes. Statuses as a `GpuStage` run's.
template <template <template <class> class> class S, typename Adapter>
[[nodiscard]] StageStatus runFlat(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, const std::uint32_t* listed, std::uint32_t listedCount, Adapter adapter, DeviceStageOutput<S>& output, FlatRunReport& report, const FlatOptions& options = {})
{
  using Index = typename DeviceStageOutput<S>::Index;
  using Record = typename DeviceStageOutput<S>::Record;
  report = {.freeBytes = arena.capacity() - arena.used()};
  const Index elements = output.elementCount();
  // CUB's item counts are ints: past them the run ends with kIndexOverflow,
  // as a `GpuStage` run's.
  if (std::size_t{elements} + 1U > static_cast<std::size_t>(std::numeric_limits<int>::max()) || std::size_t{listedCount} + 1U > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    return StageStatus::kIndexOverflow;
  }
  auto* counts = output.counts().data();
  runtime.fill(counts, 0U, (std::size_t{elements} + 1U) * sizeof(Index), stream);

  auto scratch = arena.scope(ArenaEnd::kScratch);
  // The list's offsets: the candidates of each listed element, counted and scanned.
  DeviceSpan<std::uint32_t> candidateCounts;
  DeviceSpan<std::uint32_t> offsets;
  DeviceSpan<detail::FlatWords> words;
  if (!candidateCounts.allocate(arena, std::size_t{listedCount} + 1U, ArenaEnd::kScratch) || !offsets.allocate(arena, std::size_t{listedCount} + 1U, ArenaEnd::kScratch) || !words.allocate(arena, 1U, ArenaEnd::kScratch)) {
    return StageStatus::kOutOfMemory;
  }
  detail::FlatWords* const device = words.data();
  runtime.fill(device, 0U, sizeof(detail::FlatWords), stream);
  algorithms::forEach(runtime, stream, listedCount + 1U, detail::CountElements<Adapter>{adapter, listed, listedCount, candidateCounts.data()});
  StageStatus status = algorithms::exclusiveScan(runtime, arena, stream, candidateCounts.data(), offsets.data(), std::size_t{listedCount} + 1U);
  if (status != StageStatus::kComplete) {
    return status;
  }
  // The list's length in 64 bits -- the 32-bit offsets must hold it -- and
  // its longest element; one readback.
  // A few thousand threads, each over its stride of the counts: two atomics each.
  algorithms::forEach(runtime, stream, std::min<std::uint32_t>(listedCount, 2048U), detail::ReduceCounts{candidateCounts.data(), listedCount, &device->sum, &device->largest});
  detail::FlatWords host{};
  download(words, &host, 1U, stream);
  stream.synchronize();
  ++report.readbacks;
  report.largest = host.largest;
  if (host.sum > std::numeric_limits<std::uint32_t>::max()) {
    return StageStatus::kIndexOverflow;
  }
  const auto total = static_cast<std::uint32_t>(host.sum);
  report.candidates = total;

  // The chunks: a pair, the select's index, flag and output for it (9 bytes),
  // the library's temporary storage (~8), the adapter's per-pair storage, and
  // at most one accepted record a pair on the base end; the contexts of the
  // chunk's elements and their run bounds beside.
  const std::size_t perPair = sizeof(CandidatePair) + 4U + 1U + 4U + 8U + adapter.pairBytes() + sizeof(Record);
  DeviceColumnChunks<S> chunks;
  DeviceColumns<S> accepted;
  std::uint32_t a = 0U;
  while (a < total || (total == 0U && report.chunks == 0U)) {
    auto chunkScope = arena.scope(ArenaEnd::kScratch);
    ++report.chunks;
    // The chunk: the rest of the list, halved until what it needs -- with the
    // contexts of the elements it reaches -- fits the arena's free bytes;
    // one pair at least, so the run progresses or fails for memory.
    const std::size_t free = arena.capacity() - arena.used();
    std::size_t chunkPairs = std::size_t{total} - a;
    if (options.maxChunkPairs > 0U) {
      chunkPairs = std::min(chunkPairs, options.maxChunkPairs);
    }
    chunkPairs = std::min<std::size_t>(chunkPairs, static_cast<std::size_t>(std::numeric_limits<int>::max()) - 1U);
    std::uint32_t b = a;
    std::uint32_t range[2]{0U, listedCount};
    for (;;) {
      b = static_cast<std::uint32_t>(std::size_t{a} + chunkPairs);
      // The chunk's elements: every listed one when the chunk is the whole list.
      range[0] = 0U;
      range[1] = listedCount;
      if (listedCount > 0U && !(a == 0U && b == total)) {
        detail::chunkRangeKernel<><<<1, 1, 0U, nativeStream(runtime, stream)>>>(offsets.data(), listedCount, a, b, device->range);
        checkLaunch();
        runtime.copyToHost(range, device->range, sizeof(range), stream);
        stream.synchronize();
        ++report.readbacks;
      }
      const std::uint32_t reached = range[1] > range[0] ? range[1] - range[0] : 0U;
      const std::size_t need = (perPair * std::max<std::size_t>(chunkPairs, 1U)) + adapter.contextBytes(reached) + (2U * sizeof(std::uint32_t) * reached);
      if (need <= free || total == 0U) {
        break;
      }
      if (chunkPairs <= 1U) {
        return StageStatus::kOutOfMemory;
      }
      chunkPairs /= 2U;
    }
    const std::uint32_t sourceBegin = range[0];
    const std::uint32_t sourceCount = range[1] > range[0] ? range[1] - range[0] : 0U;

    // The run bounds of the chunk's elements, starts then ends, in one block.
    DeviceSpan<CandidatePair> pairs;
    DeviceSpan<std::uint32_t> runBounds;
    const std::size_t boundCount = std::max<std::size_t>(sourceCount, 1U);
    if (!pairs.allocate(arena, std::max<std::size_t>(b - a, 1U), ArenaEnd::kScratch) || !runBounds.allocate(arena, 2U * boundCount, ArenaEnd::kScratch) || !adapter.beginChunk(arena, sourceBegin, sourceCount, b - a)) {
      return StageStatus::kOutOfMemory;
    }
    std::uint32_t* const runStart = runBounds.data();
    std::uint32_t* const runEnd = runBounds.data() + boundCount;
    runtime.fill(runBounds.data(), 0U, 2U * boundCount * sizeof(std::uint32_t), stream);
    const detail::FillElements<Adapter> fill{adapter, listed, offsets.data(), sourceBegin, pairs.data(), a, b};
    if constexpr (detail::fillLanes<Adapter>() > 1U) {
      algorithms::forEachGroup<detail::fillLanes<Adapter>()>(runtime, stream, sourceCount, fill);
    } else {
      algorithms::forEach(runtime, stream, sourceCount, fill);
    }

    DeviceSpan<std::uint32_t> selected;
    std::uint32_t selectedCount = 0U;
    if (b > a) {
      status = algorithms::selectIndices(runtime, arena, stream, b - a, detail::TestPairs<Adapter>{adapter, pairs.data()}, selected, selectedCount, ArenaEnd::kScratch);
      if (status != StageStatus::kComplete) {
        return status;
      }
      ++report.readbacks;
    }

    if (!accepted.allocate(arena, selectedCount, ArenaEnd::kBase)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(runtime, stream, selectedCount, detail::EmitPairs<S, Adapter>{adapter, listed, pairs.data(), selectedCount > 0U ? selected.data() : nullptr, selectedCount, sourceBegin, accepted.view(), runStart, runEnd});
    algorithms::forEach(runtime, stream, sourceCount, detail::AddCounts{listed, sourceBegin, runStart, runEnd, counts, &device->overflow, static_cast<std::uint32_t>(DeviceStageOutput<S>::kMaxPerElement)});
    chunks.adopt(accepted);
    a = b;
    if (total == 0U) {
      break;
    }
  }

  // The offsets; the limit flag; the records from the chunks' rows.
  status = algorithms::exclusiveScan(runtime, arena, stream, counts, output.offsets().data(), std::size_t{elements} + 1U);
  if (status != StageStatus::kComplete) {
    return status;
  }
  std::uint32_t overflow = 0U;
  runtime.copyToHost(&overflow, &device->overflow, sizeof(overflow), stream);
  stream.synchronize();
  ++report.readbacks;
  if (overflow != 0U) {
    return StageStatus::kOrdinalOverflow;
  }
  if (chunks.rows() > std::numeric_limits<Index>::max()) {
    return StageStatus::kIndexOverflow;
  }
  const auto records = static_cast<Index>(chunks.rows());
  if (report.chunks == 1U) {
    output.adoptOutput(accepted, records);
    return StageStatus::kComplete;
  }
  if (!output.allocateOutput(arena, std::max<Index>(records, 1U))) {
    return StageStatus::kOutOfMemory;
  }
  output.setSize(records);
  chunks.gather(output.output(), 0U, stream);
  return StageStatus::kComplete;
}

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::flat

#endif // GTRK_TRACKING_GPU_FLATSTAGE_CUH_
