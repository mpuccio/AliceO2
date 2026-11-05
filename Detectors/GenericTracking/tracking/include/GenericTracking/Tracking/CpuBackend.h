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

/// \file CpuBackend.h
/// \brief Running one stage's generate, scan and place on the host threads.
///
/// Elements are cut into chunks of contiguous ordinals. Chunks run in
/// parallel; inside a chunk the elements run in order, appending their records
/// to the running arena slot's buffer, so a chunk's records sit contiguous in
/// whichever buffer took them. After a two-level scan --
/// chunk totals, then elements inside each chunk -- each chunk is one
/// contiguous copy to `offset[first element]`. The output is the serial
/// reference's, whatever the thread count, chunk size or scheduling.
///
/// The tracker gives every slot's buffer the staging per thread
/// CapacityEstimator predicts; a miss grows that buffer for the run and tells
/// the estimator. The one hard limit is the bounded pool: kOutOfMemory, which
/// ends the run.

#ifndef GTRK_TRACKING_CPUBACKEND_H_
#define GTRK_TRACKING_CPUBACKEND_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/Columns.h"
#include "GenericTracking/Base/HostAlgorithms.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/Profile.h"
#include "GenericTracking/Tracking/StageOutput.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <utility>
#include <vector>

namespace o2::gtrk
{

template <template <template <class> class> class S>
class CpuStage
{
 public:
  using Output = StageOutput<S>;
  using Index = typename Output::Index;
  using Record = typename Output::Record;
  using Columns = BoundedColumns<S>;

  /// One task-arena slot's staging: its columns and rows in use.
  struct ThreadBuffer {
    Columns columns;
    std::size_t used = 0;
    bool grew = false;

    void append(const Record& record)
    {
      if (used == columns.size()) {
        grew = true;
        columns.resize(std::max<std::size_t>(kInitialRows, 2U * columns.size()));
      }
      columns.write(used++, record);
    }
  };

  /// What an element's generator emits into: the running arena slot's buffer.
  class Emitter
  {
   public:
    void emit(const Record& record)
    {
      mBuffer->append(record);
      ++mCount;
    }
    [[nodiscard]] Index count() const { return mCount; }

   private:
    friend class CpuStage;
    explicit Emitter(ThreadBuffer& buffer) : mBuffer(&buffer) {}

    ThreadBuffer* mBuffer = nullptr;
    Index mCount = 0;
  };

  /// On \a arena, which outlives the stage: the tracker's one arena, shared
  /// by its stages and loops so that TBB's workers stay with it from one loop
  /// to the next.
  explicit CpuStage(detail::TaskArena& arena) : mArena{arena}, mSlots{arena.maxConcurrency()} {}
  ~CpuStage() = default;
  CpuStage(const CpuStage&) = delete;
  CpuStage& operator=(const CpuStage&) = delete;
  CpuStage(CpuStage&&) = delete;
  CpuStage& operator=(CpuStage&&) = delete;

  /// Bind one buffer per task-arena slot and the chunk records to \a pool.
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    mBuffers.clear();
    mBuffers.reserve(static_cast<std::size_t>(mArena.maxConcurrency()));
    for (int slot = 0; slot < mArena.maxConcurrency(); ++slot) {
      mBuffers.push_back(std::make_unique<ThreadBuffer>());
    }
    releaseStorage();
  }

  /// Each thread's buffer at \a perThread records from its next run on, kept
  /// where it is larger: the timeframe's staging.
  /// The run allocates it,
  /// inside the stage's guard, so a pool too small is the stage's status.
  void prepareStagingPerThread(std::size_t perThread) { mPreparedRows = perThread; }

  /// The most records one buffer held at the end of a run since
  /// `resetLargestUse`, and whether a buffer grew in a run since.
  [[nodiscard]] std::size_t largestUse() const noexcept { return mLargestUse; }
  [[nodiscard]] bool grewSinceReset() const noexcept { return mGrew; }
  void resetLargestUse() noexcept
  {
    mLargestUse = 0U;
    mGrew = false;
  }

  /// The staging and the chunk lists back to the pool;
  /// the buffer objects stay,
  /// and the prepared size, which the next run allocates again
  /// (the tracker releases every iteration).
  void releaseStorage() noexcept
  {
    for (const auto& value : mBuffers) {
      value->columns.setMemoryResource(mMemoryPool.get());
      value->used = 0U;
      value->grew = false;
    }
    deepVectorClear(mMemoryPool.get(), mChunks, mLayerChunks, mGroupChunks);
  }

  [[nodiscard]] int maxConcurrency() const noexcept { return mArena.maxConcurrency(); }

  /// The profile phases `run` reports its reset, generate, scan and place
  /// as; none by default.
  void setPhases(const profile::StagePhases& phases) { mPhases = phases; }

  /// Allocated records across the reusable task-arena staging buffers.
  [[nodiscard]] std::size_t stagingCapacity() const
  {
    std::size_t result = 0U;
    for (const auto& buffer : mBuffers) {
      result += buffer->columns.size();
    }
    return result;
  }

  class Staged;

  /// Run \a generate(element, emitter) for every element of \a output, in
  /// chunks of \a elementsPerChunk, and place the result canonically. A
  /// generator with `filterStaged(staged)` sees every record staged before
  /// the place and drops some (`Staged`):
  /// the output holds the kept ones in canonical order, each element's count
  /// the kept count; the overflow checks and the staging's use are the
  /// generator's raw records'.
  template <typename Generate>
  StageStatus run(Output& output, Index elementsPerChunk, Generate generate)
  {
    {
      // Generate writes every element's count and the scan every offset:
      // only the logical state needs clearing.
      const profile::Zone zone{mPhases.reset};
      output.clearResult();
    }
    const Index elements = output.elementCount();
    const host_algorithms::Chunking chunking{elements, elementsPerChunk};
    const Index chunkSize = chunking.size;
    std::size_t chunkCount = chunking.count();
    // A layered generator (`Tracker::LayeredNeighbours`): its elements by
    // layer, each layer's own chunks.
    constexpr bool layered = requires {
      generate.layerBoundaries;
      generate.layerGroups;
    };
    return guardOutOfMemory([&] {
      auto counts = output.counts();
      if constexpr (layered) {
        mLayerChunks.resize(generate.layerBoundaries.size());
        chunkCount = 0U;
        for (std::size_t layer = 0U; layer < generate.layerGroups.order.size(); ++layer) {
          mLayerChunks[layer] = static_cast<Index>(chunkCount);
          const Index count = generate.layerBoundaries[layer + 1U] - generate.layerBoundaries[layer];
          chunkCount += count == 0U ? 0U : 1U + ((count - 1U) / chunkSize);
        }
        mLayerChunks.back() = static_cast<Index>(chunkCount);
      }

      // generate
      {
        const profile::Zone zone{mPhases.generate};
        // Each slot's staging in its own task: a buffer the iteration's
        // release emptied is allocated and filled here.
        detail::parallelFor(mArena, 0, mBuffers.size(), 1, [&](std::size_t first, std::size_t last) {
          for (std::size_t at = first; at < last; ++at) {
            auto& columns = mBuffers[at]->columns;
            columns.resizeRewritten(std::max(columns.size(), mPreparedRows));
          }
        });
        mChunks.resize(chunkCount);
        if constexpr (layered) {
          for (std::size_t layer = 0U; layer < generate.layerGroups.order.size(); ++layer) {
            for (auto chunk = mLayerChunks[layer]; chunk < mLayerChunks[layer + 1U]; ++chunk) {
              const auto first = generate.layerBoundaries[layer] + (chunk - mLayerChunks[layer]) * chunkSize;
              mChunks[chunk].firstElement = first;
              mChunks[chunk].endElement = static_cast<Index>(std::min<std::size_t>(generate.layerBoundaries[layer + 1U], std::size_t{first} + chunkSize));
            }
          }
        }
        for (const auto& buffer : mBuffers) {
          auto& local = *buffer;
          local.used = 0;
          local.grew = false;
        }
        const profile::Region region{mSlots, mArena.maxConcurrency(), elements};
        const auto generateChunk = [&](std::size_t chunk) {
          const auto chunkBegin = region.begin();
          const int slot = detail::currentArenaThreadIndex();
          assert(slot >= 0 && slot < mArena.maxConcurrency());
          auto& buffer = *mBuffers[static_cast<std::size_t>(slot)];
          auto& record = mChunks[chunk];
          if constexpr (!layered) {
            record.firstElement = chunking.begin(static_cast<Index>(chunk));
            record.endElement = chunking.end(static_cast<Index>(chunk));
          }
          record = {.firstElement = record.firstElement, .endElement = record.endElement, .buffer = &buffer, .begin = buffer.used};
          for (Index element = record.firstElement; element < record.endElement; ++element) {
            Emitter emitter{buffer};
            generate(element, emitter);
            counts[element] = emitter.count();
          }
          record.count = buffer.used - record.begin;
          region.end(slot, chunkBegin);
        };
        if constexpr (layered) {
          // The layers of one depth of the link DAG have no link between them
          // and every predecessor at a smaller depth (LayerLinks.h): one pass
          // per depth over its layers' chunks, listed depth by depth, with
          // every predecessor level final before the pass that reads it. The
          // chunk records keep arena order for the one scan/place below.
          const auto& groups = generate.layerGroups;
          mGroupChunks.resize(chunkCount);
          std::size_t listed = 0U;
          for (std::size_t group = 0U; group < groups.count; ++group) {
            const auto begin = listed;
            for (auto at = groups.groups[group]; at < groups.groups[group + 1U]; ++at) {
              const auto layer = groups.order[at];
              for (auto chunk = mLayerChunks[layer]; chunk < mLayerChunks[static_cast<std::size_t>(layer) + 1U]; ++chunk) {
                mGroupChunks[listed++] = chunk;
              }
            }
            if (listed != begin) {
              detail::parallelFor(mArena, begin, listed, 1U, [&](std::size_t first, std::size_t last) {
                for (std::size_t at = first; at < last; ++at) {
                  generateChunk(mGroupChunks[at]);
                }
              });
            }
          }
        } else {
          detail::parallelFor(mArena, 0, chunkCount, 1, [&](std::size_t first, std::size_t last) {
            for (std::size_t chunk = first; chunk < last; ++chunk) {
              generateChunk(chunk);
            }
          });
        }
      }

      // What the buffers held, for the staging the next timeframe is given
      // (`largestUse`).
      for (const auto& buffer : mBuffers) {
        mLargestUse = std::max(mLargestUse, buffer->used);
        mGrew = mGrew || buffer->grew;
      }

      // scan: chunk totals serially, then the elements inside each chunk, the
      // chunks in parallel
      std::uint64_t total = 0;
      bool ordinalOverflow = false;
      const auto scan = [&] {
        const profile::Zone zone{mPhases.scan};
        total = 0;
        for (auto& record : mChunks) {
          record.offset = total;
          total += record.count;
        }
        if (total > std::numeric_limits<Index>::max()) {
          return StageStatus::kIndexOverflow;
        }
        auto offsets = output.offsets();
        // A region, so the tree reads the loop as the parallel work it is.
        const profile::Region region{mSlots, mArena.maxConcurrency(), elements};
        detail::parallelFor(mArena, 0, chunkCount, 1, [&](std::size_t first, std::size_t last) {
          for (std::size_t chunk = first; chunk < last; ++chunk) {
            const auto chunkBegin = region.begin();
            auto& record = mChunks[chunk];
            auto running = static_cast<Index>(record.offset);
            bool overflow = false;
            for (Index element = record.firstElement; element < record.endElement; ++element) {
              offsets[element] = running;
              running += counts[element];
              overflow = overflow || counts[element] > Output::kMaxPerElement;
            }
            record.ordinalOverflow = overflow;
            region.end(detail::currentArenaThreadIndex(), chunkBegin);
          }
        });
        ordinalOverflow = std::ranges::any_of(mChunks, [](const Chunk& record) { return record.ordinalOverflow; });
        offsets[elements] = static_cast<Index>(total);
        return StageStatus::kComplete;
      };
      if (const auto scanned = scan(); scanned != StageStatus::kComplete) {
        return scanned;
      }
      profile::count(mPhases.generate, 0U, total);
      if (ordinalOverflow) {
        return StageStatus::kOrdinalOverflow;
      }
      if constexpr (requires(const Staged& staged) { generate.filterStaged(staged); }) {
        generate.filterStaged(Staged{*this, counts});
        static_cast<void>(scan());
      }

      // place: one contiguous copy per chunk
      {
        const profile::Zone zone{mPhases.place};
        output.resizeOutput(static_cast<Index>(total));
        const profile::Region region{mSlots, mArena.maxConcurrency(), total};
        detail::parallelFor(mArena, 0, chunkCount, 1, [&](std::size_t first, std::size_t last) {
          for (std::size_t chunk = first; chunk < last; ++chunk) {
            const auto chunkBegin = region.begin();
            const auto& record = mChunks[chunk];
            record.buffer->columns.copyTo(record.begin, output.output(), record.offset, record.count);
            region.end(detail::currentArenaThreadIndex(), chunkBegin);
          }
        });
      }
      return StageStatus::kComplete;
    });
  }

  /// Run \a body(element) for every element in [0, \a elements), in chunks of
  /// \a elementsPerChunk. For work on disjoint parts of an output, such as
  /// sorting each element's range, which a device runs one thread per element.
  template <typename Body>
  void forEachElement(Index elements, Index elementsPerChunk, Body body)
  {
    const host_algorithms::Chunking chunking{elements, elementsPerChunk};
    const profile::Region region{mSlots, mArena.maxConcurrency(), elements};
    detail::parallelFor(mArena, 0, chunking.count(), 1, [&](std::size_t first, std::size_t last) {
      for (std::size_t chunk = first; chunk < last; ++chunk) {
        const auto chunkBegin = region.begin();
        const auto end = chunking.end(static_cast<Index>(chunk));
        for (auto element = chunking.begin(static_cast<Index>(chunk)); element < end; ++element) {
          body(element);
        }
        region.end(detail::currentArenaThreadIndex(), chunkBegin);
      }
    });
  }

  /// Sort \a indices by \a less, a strict total order on them, on this
  /// stage's task arena (`detail::parallelSort`): the serial sort's result,
  /// on the threads the arena allows and no others.
  template <typename Less>
  void sortIndices(std::span<std::uint32_t> indices, Less less)
  {
    detail::parallelSort(mArena, indices, less);
  }

  /// A filtering generator's view of the records staged by its run
  /// (`run`, `filterStaged`), after the raw scan: each element's records, and
  /// their compaction.
  class Staged
  {
   public:
    /// body(element, columns, at) for each record the elements in
    /// [\a begin, \a end) staged: row \a at of the staging's \a columns (a
    /// const view), so that a body reads only the columns it needs. The
    /// chunks in parallel, each element's records in order on one thread.
    template <typename Body>
    void forEach(Index begin, Index end, const Body& body) const
    {
      if (begin >= end) {
        return;
      }
      auto& stage = *mStage;
      const auto& chunks = stage.mChunks;
      // The chunks holding the elements: from the last that starts at or before begin.
      const auto first = static_cast<std::size_t>((std::ranges::upper_bound(chunks, begin, {}, &Chunk::firstElement) - chunks.begin()) - 1);
      const auto last = static_cast<std::size_t>(std::ranges::lower_bound(chunks, end, {}, &Chunk::firstElement) - chunks.begin());
      // A region, so the tree reads the loop as the parallel work it is.
      const profile::Region region{stage.mSlots, stage.mArena.maxConcurrency(), end - begin};
      detail::parallelFor(stage.mArena, first, last, 1, [&](std::size_t from, std::size_t to) {
        for (std::size_t at = from; at < to; ++at) {
          const auto chunkBegin = region.begin();
          const auto& chunk = chunks[at];
          const auto columns = std::as_const(chunk.buffer->columns).view();
          std::size_t record = chunk.begin;
          for (Index element = chunk.firstElement; element < chunk.endElement; ++element) {
            const std::size_t next = record + mCounts[element];
            if (element >= begin && element < end) {
              for (; record < next; ++record) {
                body(element, columns, record);
              }
            }
            record = next;
          }
          region.end(detail::currentArenaThreadIndex(), chunkBegin);
        }
      });
    }

    /// The staged records kept where \a keep(element, columns, at) holds
    /// (as `forEach`), in order, every chunk compacted in place in parallel,
    /// the elements' counts lowered to their kept records'.
    template <typename Keep>
    void compact(const Keep& keep) const
    {
      auto& stage = *mStage;
      const profile::Region region{stage.mSlots, stage.mArena.maxConcurrency(), mCounts.size()};
      detail::parallelFor(stage.mArena, 0, stage.mChunks.size(), 1, [&](std::size_t first, std::size_t last) {
        for (std::size_t at = first; at < last; ++at) {
          const auto chunkBegin = region.begin();
          auto& chunk = stage.mChunks[at];
          auto& columns = chunk.buffer->columns;
          const auto view = std::as_const(columns).view();
          std::size_t read = chunk.begin;
          std::size_t write = chunk.begin;
          for (Index element = chunk.firstElement; element < chunk.endElement; ++element) {
            Index kept = 0U;
            for (Index i = 0U; i < mCounts[element]; ++i, ++read) {
              if (keep(element, view, read)) {
                if (write != read) {
                  columns.write(write, columns.read(read));
                }
                ++write;
                ++kept;
              }
            }
            mCounts[element] = kept;
          }
          chunk.count = write - chunk.begin;
          region.end(detail::currentArenaThreadIndex(), chunkBegin);
        }
      });
    }

   private:
    friend class CpuStage;
    Staged(CpuStage& stage, std::span<Index> counts) : mStage(&stage), mCounts(counts) {}

    CpuStage* mStage = nullptr;
    std::span<Index> mCounts;
  };

 private:
  static constexpr std::size_t kInitialRows = 1024;

  /// Where one chunk's records sit before placement.
  struct Chunk {
    Index firstElement = 0;
    Index endElement = 0;
    ThreadBuffer* buffer = nullptr;
    std::size_t begin = 0;
    std::size_t count = 0;
    std::uint64_t offset = 0;
    bool ordinalOverflow = false; ///< an element of the chunk emitted more than an ordinal numbers
  };

  // Pool must outlive allocator-backed members, the thread buffers included.
  detail::TaskArena& mArena;
  profile::RegionSlots mSlots; ///< one per arena slot, for the chunk timings
  profile::StagePhases mPhases;
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  std::vector<std::unique_ptr<ThreadBuffer>> mBuffers;
  BoundedVector<Chunk> mChunks;
  BoundedVector<Index> mLayerChunks;
  BoundedVector<Index> mGroupChunks; ///< the layered generation's chunks, depth by depth
  std::size_t mPreparedRows = kInitialRows;
  std::size_t mLargestUse = 0U; ///< `largestUse`
  bool mGrew = false;           ///< `grewSinceReset`
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_CPUBACKEND_H_
