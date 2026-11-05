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

/// \file StageOutput.h
/// \brief One stage's output in canonical order: generate, scan, place.
///
/// A stage runs one element after another in the serial reference; each
/// element emits a variable number of records of `S`. The canonical output is
/// what that serial loop produces, in that order: element 0's records, then element 1's, and so
/// on. Any backend reproduces it with three steps:
///
///  - generate: every element emits its records and its count;
///  - scan: an exclusive scan of the counts gives `offset`, and element `e`
///    owns `[offset[e], offset[e+1])` of the output;
///  - place: record `k` of element `e` goes to `offset[e] + k`.
///
/// This owns the counts, the offsets and the placed output; each backend
/// stages its records itself (CpuBackend.h per thread, GpuStage.cuh by
/// element and ordinal, the tests' serial reference test/SerialReferences.h as a
/// device does). Nothing here knows what an element is; a stage maps element
/// ordinals to its own inputs.

#ifndef GTRK_TRACKING_STAGEOUTPUT_H_
#define GTRK_TRACKING_STAGEOUTPUT_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/Columns.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <span>
#include <utility>

namespace o2::gtrk
{

/// How a stage run ended. Only kComplete leaves a usable output; everything
/// else is discarded whole and tells the controller which rung of the overflow
/// ladder to take.
enum class StageStatus : std::uint8_t {
  kComplete,
  kStagingOverflow, ///< staging was too small; the backend's staged count is the exact size a retry needs
  kOrdinalOverflow, ///< an element emitted more records than an ordinal can number
  kIndexOverflow,   ///< the total does not fit the output index type
  kOutOfMemory,     ///< the bounded pool refused an allocation: the run ends with this status
};

/// Run \a body, turning a refusal of the bounded pool into kOutOfMemory.
template <typename Body>
StageStatus guardOutOfMemory(Body body)
{
  try {
    return body();
  } catch (const std::bad_alloc&) {
    return StageStatus::kOutOfMemory;
  }
}

template <template <template <class> class> class S>
class StageOutput
{
 public:
  using Index = std::uint32_t;
  using Ordinal = std::uint16_t;
  using Record = S<memlayout::value>;
  using Columns = BoundedColumns<S>;

  /// The most records one element may emit. The limit is the staged ordinal's
  /// width and holds on every backend, so that they fail alike.
  static constexpr std::size_t kMaxPerElement = std::size_t{std::numeric_limits<Ordinal>::max()} + 1U;

  StageOutput() = default;
  ~StageOutput() = default;
  StageOutput(const StageOutput&) = delete;
  StageOutput& operator=(const StageOutput&) = delete;
  StageOutput(StageOutput&&) = delete;
  StageOutput& operator=(StageOutput&&) = delete;

  /// Reseat every allocation on \a pool. Call before the first prepare().
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    releaseStorage();
  }

  /// Size the counts and offsets for \a elements, and clear the result. The
  /// counts and offsets keep whatever they held: every run writes them whole
  /// before reading them.
  void prepare(Index elements)
  {
    resizeBoundedVector(mCount, elements);
    resizeBoundedVector(mOffset, std::size_t{elements} + 1U);
    clearResult();
  }

  /// No records; the counts and offsets left as they are, for a run that
  /// writes them all.
  void clearResult() noexcept { mSize = 0; }

  /// Release every allocation back to the pool.
  void clearStorage() noexcept
  {
    releaseStorage();
    mSize = 0;
  }

  // ---- building blocks for the backends ----

  /// Size the placed output to \a total records, which completes it once
  /// every record is placed. The columns only grow: rows past \a total are
  /// kept, so that a run no smaller than an earlier one zeroes nothing it then
  /// overwrites; grown, they are rebuilt at
  /// \a total, as the place writes them all. The output is `size()` rows.
  void resizeOutput(Index total)
  {
    if (mOutput.size() < total) {
      mOutput.resizeRewritten(total);
    }
    mSize = total;
  }

  [[nodiscard]] std::span<Index> counts() { return mCount; }
  [[nodiscard]] std::span<Index> offsets() { return mOffset; }
  [[nodiscard]] Columns& output() { return mOutput; }

  // ---- the result ----

  [[nodiscard]] Index elementCount() const { return static_cast<Index>(mCount.size()); }
  [[nodiscard]] std::span<const Index> counts() const { return mCount; }
  [[nodiscard]] std::span<const Index> offsets() const { return mOffset; }
  [[nodiscard]] Index size() const { return mSize; }
  [[nodiscard]] typename Columns::View view() { return mOutput.view(); }
  [[nodiscard]] typename Columns::ConstView view() const { return mOutput.view(); }

 private:
  void releaseStorage() noexcept
  {
    auto* resource = mMemoryPool.get();
    mOutput.setMemoryResource(resource);
    deepVectorClear(resource, mCount, mOffset);
  }

  // Pool must outlive allocator-backed members.
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  Columns mOutput;
  BoundedVector<Index> mCount;
  BoundedVector<Index> mOffset;
  Index mSize = 0;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_STAGEOUTPUT_H_
