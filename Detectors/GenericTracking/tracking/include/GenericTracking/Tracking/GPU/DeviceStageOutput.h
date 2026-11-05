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

/// \file DeviceStageOutput.h
/// \brief One stage's canonical output in device memory: the device side of
///        StageOutput<S>.
///
/// `count` and `offset` have `n + 1` entries each - the last count is zero so
/// that one exclusive scan gives `offset[n]`, the total - and the placed
/// columns hold the stage's records in serial order. Everything lives on the
/// base end of the arena, in the caller's scope; the working storage of a run
/// is its stage's (GpuStage.cuh, FlatStage.cuh), on the scratch end; a flat run
/// of several chunks also leaves its chunks' records on the base end
/// (FlatStage.cuh).

#ifndef GTRK_TRACKING_GPU_DEVICESTAGEOUTPUT_H_
#define GTRK_TRACKING_GPU_DEVICESTAGEOUTPUT_H_

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/StageOutput.h"

#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace o2::gtrk
{

/// What a stage's run over a flat candidate list reports (`flat::runFlat`).
struct FlatRunReport {
  std::uint64_t candidates = 0U; ///< the list's candidates, every element's
  std::uint32_t largest = 0U;    ///< the most candidates of one element
  std::uint32_t chunks = 0U;     ///< the chunks the list was run in
  std::uint32_t readbacks = 0U;
  std::size_t freeBytes = 0U; ///< the arena's free bytes when the run starts
};

template <template <template <class> class> class S>
class DeviceStageOutput
{
 public:
  using Host = StageOutput<S>;
  using Index = typename Host::Index;
  using Ordinal = typename Host::Ordinal;
  using Record = typename Host::Record;
  using Columns = DeviceColumns<S>;
  static constexpr std::size_t kMaxPerElement = Host::kMaxPerElement;

  /// Counts and offsets for \a elements elements. False if the arena refuses.
  [[nodiscard]] bool prepare(DeviceArena& arena, Index elements)
  {
    mElements = elements;
    mSize = 0U;
    mOutputReady = false;
    return mCounts.allocate(arena, std::size_t{elements} + 1U) && mOffsets.allocate(arena, std::size_t{elements} + 1U);
  }

  /// The placed columns, once a run knows the total. False if the arena refuses.
  [[nodiscard]] bool allocateOutput(DeviceArena& arena, Index total)
  {
    mOutputReady = mOutput.allocate(arena, total);
    return mOutputReady;
  }

  /// A completed run's record count.
  void setSize(Index total) noexcept { mSize = total; }

  /// \a columns, already holding the \a total placed records, as the output
  /// (a flat run's one chunk is its output).
  void adoptOutput(const Columns& columns, Index total) noexcept
  {
    mOutput = columns;
    mOutputReady = true;
    mSize = total;
  }

  [[nodiscard]] Index elementCount() const noexcept { return mElements; }
  [[nodiscard]] Index size() const noexcept { return mSize; }
  [[nodiscard]] DeviceSpan<Index>& counts() noexcept { return mCounts; }
  [[nodiscard]] const DeviceSpan<Index>& counts() const noexcept { return mCounts; }
  [[nodiscard]] DeviceSpan<Index>& offsets() noexcept { return mOffsets; }
  [[nodiscard]] const DeviceSpan<Index>& offsets() const noexcept { return mOffsets; }
  [[nodiscard]] Columns& output()
  {
    requireOutput();
    return mOutput;
  }
  [[nodiscard]] const Columns& output() const
  {
    requireOutput();
    return mOutput;
  }
  [[nodiscard]] typename Columns::View view() { return output().view(); }
  [[nodiscard]] typename Columns::ConstView view() const { return output().view(); }

 private:
  void requireOutput() const
  {
    if (!mOutputReady) {
      throw std::logic_error("the device stage output has no placed records: no run has completed");
    }
  }

  Index mElements = 0U;
  Index mSize = 0U;
  bool mOutputReady = false;
  DeviceSpan<Index> mCounts;
  DeviceSpan<Index> mOffsets;
  Columns mOutput;
};

/// Copy a completed device stage into \a host - counts, offsets and the placed
/// columns - so that it compares with a host StageOutput byte for byte. The
/// copies are ordered on \a stream; the caller synchronises before reading.
template <template <template <class> class> class S>
void download(const DeviceStageOutput<S>& device, StageOutput<S>& host, DeviceStream& stream)
{
  const auto elements = device.elementCount();
  host.prepare(elements);
  download(device.counts(), host.counts().data(), elements, stream);
  download(device.offsets(), host.offsets().data(), std::size_t{elements} + 1U, stream);
  host.resizeOutput(device.size());
  download(device.output(), host.output(), device.size(), stream);
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_DEVICESTAGEOUTPUT_H_
