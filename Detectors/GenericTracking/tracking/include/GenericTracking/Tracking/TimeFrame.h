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

/// \file TimeFrame.h
/// \brief Read-only grid-ordered measurements and block-aligned slices.

#ifndef GTRK_TRACKING_TIMEFRAME_H_
#define GTRK_TRACKING_TIMEFRAME_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/MeasurementRange.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>

namespace o2::gtrk
{

/// The measurement columns and grid offsets belong to the caller -- the DPL
/// input, a replay's storage or the device tracker's allocations
/// (DeviceTracker.cuh): the frame borrows both, with the geometry that
/// interprets their order.
template <DetectorLayout Det>
class TimeFrame
{
 public:
  using Index = std::uint32_t;
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);

  /// Core blocks partition the timeframe. Halo blocks may overlap any number
  /// of slices. The row ranges use global indices into the shared SoA columns.
  class TimeSlice
  {
   public:
    TimeSlice(const TimeFrame& frame, BlockRange coreRange, std::int64_t haloPS)
    {
      assert(coreRange.begin >= 0 && coreRange.begin <= coreRange.end && coreRange.end <= frame.grid().geometry().blockCount() && haloPS >= 0);
      const auto blocks = frame.grid().geometry().slice(coreRange, haloPS);
      coreBlocks = blocks.core;
      haloBlocks = blocks.halo;
      for (std::size_t layer = 0; layer < kLayerCount; ++layer) {
        const auto id = static_cast<LayerID>(layer);
        core[layer] = frame.grid().rows(id, coreBlocks);
        halo[layer] = frame.grid().rows(id, haloBlocks);
      }
    }

    BlockRange coreBlocks{};
    BlockRange haloBlocks{};
    std::array<MeasurementRange, kLayerCount> core{};
    std::array<MeasurementRange, kLayerCount> halo{};
  };

  /// Requires equal column lengths, `bucketCount() + 1` offsets, monotonic
  /// offsets and nondecreasing box starts within every bucket, as the
  /// producer writes them.
  TimeFrame(MeasurementView measurements, const GridGeometry<Det>& geometry, std::span<const Index> offsets)
    : TimeFrame(measurements, geometry, offsets, offsets.empty() ? 0U : offsets.back())
  {
  }

  /// The same frame over columns and offsets the host cannot read -- a
  /// device copy -- given the row count, which the other constructor reads
  /// from the last offset. Nothing here dereferences them; slices are cut
  /// from the host frame.
  TimeFrame(MeasurementView measurements, const GridGeometry<Det>& geometry, std::span<const Index> offsets, Index rows)
    : mGrid(measurements, geometry, offsets), mSize(rows)
  {
  }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() Index size() const { return mSize; }
  [[nodiscard]] GPUhdi() GPUAlwaysInline() auto measurement(Index index) const
  {
    assert(index < size());
    return measurements()[index];
  }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() const MeasurementGrid<Det>& grid() const { return mGrid; }
  /// The columns the frame borrows.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() MeasurementView measurements() const { return mGrid.measurements(); }

 private:
  MeasurementGrid<Det> mGrid;
  Index mSize = 0;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TIMEFRAME_H_
