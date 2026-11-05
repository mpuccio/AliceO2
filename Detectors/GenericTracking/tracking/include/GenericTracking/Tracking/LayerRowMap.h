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

/// \file LayerRowMap.h
/// \brief A dense numbering of the rows of chosen layers in a slice.
///
/// A stage that runs one element per measurement -- the tracklet stage, one
/// per source -- sees a slice as one contiguous row range per layer
/// (`TimeFrame<Det>::TimeSlice`), and presents only some layers. This numbers
/// the presented rows `0 .. size()-1`, layer after layer in row order, which
/// is the stage's element order and so its canonical output order. The
/// mapping is arithmetic both ways: a row to its element through the layer's
/// offset, an element to its layer by a binary search over the layer offsets.

#ifndef GTRK_TRACKING_LAYERROWMAP_H_
#define GTRK_TRACKING_LAYERROWMAP_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/MeasurementRange.h"

#include "GPUCommonAlgorithm.h"
#include "GPUCommonDef.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>

namespace o2::gtrk
{

template <DetectorLayout Det>
struct LayerRowMap {
  using Index = std::uint32_t;
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);

  struct Row {
    LayerID layer = kInvalidLayerID;
    Index row = 0;
  };

  LayerRowMap() = default;

  /// \a rows per layer; only the layers \a presented marks contribute.
  LayerRowMap(std::span<const MeasurementRange, kLayerCount> rows, const std::array<bool, kLayerCount>& presented)
  {
    for (std::size_t layer = 0; layer < kLayerCount; ++layer) {
      ranges[layer] = presented[layer] ? rows[layer] : MeasurementRange{rows[layer].begin, rows[layer].begin};
      offsets[layer + 1U] = offsets[layer] + ranges[layer].size();
    }
  }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() Index size() const { return offsets[kLayerCount]; }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool contains(LayerID layer, Index row) const
  {
    return static_cast<std::size_t>(layer) < kLayerCount && row >= ranges[layer].begin && row < ranges[layer].end;
  }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() Index element(LayerID layer, Index row) const
  {
    assert(contains(layer, row));
    return offsets[layer] + (row - ranges[layer].begin);
  }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() Row row(Index element) const
  {
    assert(element < size());
    const std::size_t low =
      o2::gpu::CAAlgo::partitionPoint(std::size_t{1}, kLayerCount - 1U, [&](std::size_t layer) { return offsets[layer] <= element; }) - 1U;
    return {.layer = static_cast<LayerID>(low), .row = ranges[low].begin + (element - offsets[low])};
  }

  std::array<MeasurementRange, kLayerCount> ranges{};
  std::array<Index, kLayerCount + 1U> offsets{};
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_LAYERROWMAP_H_
