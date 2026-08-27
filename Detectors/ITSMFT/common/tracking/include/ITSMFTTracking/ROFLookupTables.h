// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef TRACKINGITSU_INCLUDE_ROFOVERLAPTABLE_H_
#define TRACKINGITSU_INCLUDE_ROFOVERLAPTABLE_H_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>
#include <ranges>

#ifndef GPUCA_GPUCODE
#include <format>
#include "Framework/Logger.h"
#endif

#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/RangeReference.h"
#include "DataFormatsITS/TimeEstBC.h"
#include "DataFormatsITS/Vertex.h"
#include "GPUCommonMath.h"
#include "GPUCommonDef.h"
#include "ITSMFTTracking/ROFViews.h"

namespace o2::its
{

using LayerTiming = o2::itsmft::tracking::ROFTimingLayer;

// Base class for lookup to define layers
template <int32_t NLayers>
class LayerTimingBase
{
 protected:
  LayerTiming mLayers[NLayers];

 public:
  using T = LayerTiming::BCType;
  LayerTimingBase() = default;

  GPUh() void defineLayer(int32_t layer, T nROFsTF, T rofLength, T rofDelay, T rofBias, T rofTE)
  {
    assert(layer >= 0 && layer < NLayers);
    mLayers[layer] = {nROFsTF, rofLength, rofDelay, rofBias, rofTE};
  }

  GPUh() void defineLayer(int32_t layer, const LayerTiming& timing)
  {
    assert(layer >= 0 && layer < NLayers);
    mLayers[layer] = timing;
  }

  GPUhdi() const LayerTiming& getLayer(int32_t layer) const
  {
    assert(layer >= 0 && layer < NLayers);
    return mLayers[layer];
  }

  GPUhdi() constexpr int32_t getEntries() noexcept { return NLayers; }

#ifndef GPUCA_GPUCODE
  GPUh() void print() const
  {
    LOGP(info, "Imposed time structure:");
    for (int32_t iL{0}; iL < NLayers; ++iL) {
      LOGP(info, "\tLayer:{} {}", iL, mLayers[iL].asString());
    }
  }
#endif
};

// GPU friendly view of the table below
template <int32_t NLayers, typename TableEntry, typename TableIndex>
using ROFOverlapTableView = o2::itsmft::tracking::ROFOverlapView<TableEntry, TableIndex>;

// Precalculated lookup table to find overlapping ROFs in another layer given a ROF index in the current layer
template <int32_t NLayers>
class ROFOverlapTable : public LayerTimingBase<NLayers>
{
 public:
  using T = LayerTimingBase<NLayers>::T;
  using TableEntry = dataformats::RangeReference<T, T>;
  using TableIndex = dataformats::RangeReference<T, T>;

  using View = ROFOverlapTableView<NLayers, TableEntry, TableIndex>;
  ROFOverlapTable() = default;

  GPUh() void init()
  {
    std::vector<TableEntry> table[NLayers][NLayers];
    for (int32_t i{0}; i < NLayers; ++i) {
      for (int32_t j{0}; j < NLayers; ++j) {
        if (i != j) { // we do not need self-lookup
          buildMapping(i, j, table[i][j]);
        }
      }
    }
    flatten(table);
  }

  GPUh() View getView() const
  {
    View view;
    view.mFlatTable = mFlatTable.data();
    view.mIndices = mIndices;
    view.mLayers = this->mLayers;
    view.mLayerCount = NLayers;
    return view;
  }

  GPUh() View getDeviceView(const TableEntry* deviceFlatTablePtr, const TableIndex* deviceIndicesPtr, const LayerTiming* deviceLayerTimingPtr) const
  {
    View view;
    view.mFlatTable = deviceFlatTablePtr;
    view.mIndices = deviceIndicesPtr;
    view.mLayers = deviceLayerTimingPtr;
    view.mLayerCount = NLayers;
    return view;
  }

  GPUh() size_t getFlatTableSize() const noexcept { return mFlatTable.size(); }
  static GPUh() constexpr size_t getIndicesSize() { return static_cast<size_t>(NLayers * NLayers); }

 private:
  GPUh() void buildMapping(int32_t from, int32_t to, std::vector<TableEntry>& table)
  {
    const auto& layerFrom = this->mLayers[from];
    const auto& layerTo = this->mLayers[to];
    table.resize(layerFrom.mNROFsTF);

    for (int32_t iROF{0}; iROF < layerFrom.mNROFsTF; ++iROF) {
      int64_t fromStart = o2::gpu::CAMath::Max((int64_t)layerFrom.getROFStartInBC(iROF) - (int64_t)layerFrom.mROFAddTimeErr, int64_t(0));
      int64_t fromEnd = (int64_t)layerFrom.getROFEndInBC(iROF) + layerFrom.mROFAddTimeErr;

      int32_t firstROFTo = o2::gpu::CAMath::Max(0, (int32_t)((fromStart - (int64_t)layerTo.mROFAddTimeErr - (int64_t)layerTo.mROFDelay - (int64_t)layerTo.mROFBias) / (int64_t)layerTo.mROFLength));
      auto lastROFTo = (int32_t)((fromEnd + (int64_t)layerTo.mROFAddTimeErr - (int64_t)layerTo.mROFDelay - (int64_t)layerTo.mROFBias - 1) / (int64_t)layerTo.mROFLength);
      firstROFTo = o2::gpu::CAMath::Max(0, firstROFTo);
      lastROFTo = o2::gpu::CAMath::Min((int32_t)layerTo.mNROFsTF - 1, lastROFTo);

      while (firstROFTo <= lastROFTo) {
        int64_t toStart = o2::gpu::CAMath::Max((int64_t)layerTo.getROFStartInBC(firstROFTo) - (int64_t)layerTo.mROFAddTimeErr, int64_t(0));
        int64_t toEnd = (int64_t)layerTo.getROFEndInBC(firstROFTo) + layerTo.mROFAddTimeErr;
        if (toEnd > fromStart && toStart < fromEnd) {
          break;
        }
        ++firstROFTo;
      }
      while (lastROFTo >= firstROFTo) {
        int64_t toStart = o2::gpu::CAMath::Max((int64_t)layerTo.getROFStartInBC(lastROFTo) - (int64_t)layerTo.mROFAddTimeErr, int64_t(0));
        int64_t toEnd = (int64_t)layerTo.getROFEndInBC(lastROFTo) + layerTo.mROFAddTimeErr;
        if (toEnd > fromStart && toStart < fromEnd) {
          break;
        }
        --lastROFTo;
      }
      int32_t count = (firstROFTo <= lastROFTo) ? (lastROFTo - firstROFTo + 1) : 0;
      table[iROF] = {static_cast<T>(firstROFTo), static_cast<T>(count)};
    }
  }

  GPUh() void flatten(const std::vector<TableEntry> table[NLayers][NLayers])
  {
    size_t total{0};
    for (int32_t i{0}; i < NLayers; ++i) {
      for (int32_t j{0}; j < NLayers; ++j) {
        if (i != j) { // we do not need self-lookup
          total += table[i][j].size();
        }
      }
    }

    mFlatTable.reserve(total);

    for (int32_t i{0}; i < NLayers; ++i) {
      for (int32_t j{0}; j < NLayers; ++j) {
        size_t idx = (i * NLayers) + j;
        if (i != j) {
          mIndices[idx].setFirstEntry(static_cast<T>(mFlatTable.size()));
          mIndices[idx].setEntries(static_cast<T>(table[i][j].size()));
          mFlatTable.insert(mFlatTable.end(), table[i][j].begin(), table[i][j].end());
        } else {
          mIndices[idx] = {0, 0};
        }
      }
    }
  }

  TableIndex mIndices[NLayers * NLayers];
  std::vector<TableEntry> mFlatTable;
};

// GPU friendly view of the table below
template <int32_t NLayers, typename TableEntry, typename TableIndex>
using ROFVertexLookupTableView = o2::itsmft::tracking::ROFVertexLookupView<TableEntry, TableIndex>;

// Precalculated lookup table to find vertices compatible with ROFs
// Given a layer and ROF index, returns the range of vertices that overlap in time.
// The vertex time is defined as symmetrical [t0-e,t0+e]
// It needs to be guaranteed that the input vertices are sorted by their lower-bound!
// additionally compatibliyty has to be queried per vertex!
template <int32_t NLayers>
class ROFVertexLookupTable : public LayerTimingBase<NLayers>
{
 public:
  using T = LayerTimingBase<NLayers>::T;
  using BCType = LayerTiming::BCType;
  using TableEntry = dataformats::RangeReference<T, T>;
  using TableIndex = dataformats::RangeReference<T, T>;
  using View = ROFVertexLookupTableView<NLayers, TableEntry, TableIndex>;

  ROFVertexLookupTable() = default;

  GPUh() size_t getFlatTableSize() const noexcept { return mFlatTable.size(); }
  static GPUh() constexpr size_t getIndicesSize() { return NLayers; }

  // Build the lookup table given a sorted array of vertices
  // vertices must be sorted by timestamp, then by error (secondary)
  GPUh() void init(const Vertex* vertices, size_t nVertices)
  {
    if (nVertices > std::numeric_limits<T>::max()) {
      LOGF(fatal, "too many vertices %zu, max supported is %u", nVertices, std::numeric_limits<T>::max());
    }

    std::vector<TableEntry> table[NLayers];
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      buildMapping(layer, vertices, nVertices, table[layer]);
    }
    flatten(table);
  }

  // Pre-allocated needed memory, then use update(...)
  GPUh() void init()
  {
    size_t total{0};
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      total += this->mLayers[layer].mNROFsTF;
    }
    mFlatTable.resize(total, {0, 0});
    size_t offset = 0;
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      size_t nROFs = this->mLayers[layer].mNROFsTF;
      mIndices[layer].setFirstEntry(static_cast<T>(offset));
      mIndices[layer].setEntries(static_cast<T>(nROFs));
      offset += nROFs;
    }
  }

  // Recalculate lookup table with new vertices
  GPUh() void update(const Vertex* vertices, size_t nVertices)
  {
    size_t offset = 0;
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      const auto& idx = mIndices[layer];
      size_t nROFs = idx.getEntries();
      for (size_t iROF = 0; iROF < nROFs; ++iROF) {
        updateROFMapping(layer, iROF, vertices, nVertices, offset + iROF);
      }
      offset += nROFs;
    }
  }

  GPUh() View getView() const
  {
    View view;
    view.mFlatTable = mFlatTable.data();
    view.mIndices = mIndices;
    view.mLayers = this->mLayers;
    view.mLayerCount = NLayers;
    return view;
  }

  GPUh() View getDeviceView(const TableEntry* deviceFlatTablePtr, const TableIndex* deviceIndicesPtr, const LayerTiming* deviceLayerTimingPtr) const
  {
    View view;
    view.mFlatTable = deviceFlatTablePtr;
    view.mIndices = deviceIndicesPtr;
    view.mLayers = deviceLayerTimingPtr;
    view.mLayerCount = NLayers;
    return view;
  }

 private:
  // Build the mapping for one layer
  GPUh() void buildMapping(int32_t layer, const Vertex* vertices, size_t nVertices, std::vector<TableEntry>& table)
  {
    const auto& layerDef = this->mLayers[layer];
    table.resize(layerDef.mNROFsTF);
    size_t vertexSearchStart = 0;
    for (int32_t iROF{0}; iROF < layerDef.mNROFsTF; ++iROF) {
      int64_t rofLower = o2::gpu::CAMath::Max((int64_t)layerDef.getROFStartInBC(iROF) - (int64_t)layerDef.mROFAddTimeErr, int64_t(0));
      int64_t rofUpper = (int64_t)layerDef.getROFEndInBC(iROF) + layerDef.mROFAddTimeErr;
      size_t lastVertex = binarySearchFirst(vertices, nVertices, vertexSearchStart, rofUpper);
      size_t firstVertex = vertexSearchStart;
      while (firstVertex < lastVertex) {
        auto vUpper = (int64_t)vertices[firstVertex].getTimeStamp().upper();
        if (vUpper > rofLower) {
          break;
        }
        ++firstVertex;
      }
      size_t count = (lastVertex > firstVertex) ? (lastVertex - firstVertex) : 0;
      table[iROF] = {static_cast<T>(firstVertex), static_cast<T>(count)};
      vertexSearchStart = firstVertex;
    }
  }

  // Update a single ROF's vertex mapping
  GPUh() void updateROFMapping(int32_t layer, size_t iROF, const Vertex* vertices, size_t nVertices, size_t flatTableIdx)
  {
    const auto& layerDef = this->mLayers[layer];
    int64_t rofLower = o2::gpu::CAMath::Max((int64_t)layerDef.getROFStartInBC(iROF) - (int64_t)layerDef.mROFAddTimeErr, int64_t(0));
    int64_t rofUpper = (int64_t)layerDef.getROFEndInBC(iROF) + layerDef.mROFAddTimeErr;
    size_t lastVertex = binarySearchFirst(vertices, nVertices, 0, rofUpper);
    size_t firstVertex = 0;
    while (firstVertex < lastVertex) {
      int64_t vUpper = (int64_t)vertices[firstVertex].getTimeStamp().getTimeStamp() +
                       (int64_t)vertices[firstVertex].getTimeStamp().getTimeStampError();
      if (vUpper > rofLower) {
        break;
      }
      ++firstVertex;
    }
    size_t count = (lastVertex > firstVertex) ? (lastVertex - firstVertex) : 0;
    mFlatTable[flatTableIdx].setFirstEntry(static_cast<T>(firstVertex));
    mFlatTable[flatTableIdx].setEntries(static_cast<T>(count));
  }

  // Binary search for first vertex where lowerBC >= targetBC
  GPUh() size_t binarySearchFirst(const Vertex* vertices, size_t nVertices, size_t searchStart, BCType targetBC) const
  {
    size_t left = searchStart;
    size_t right = nVertices;
    while (left < right) {
      size_t mid = left + ((right - left) / 2);
      int64_t lower = (int64_t)vertices[mid].getTimeStamp().lower();
      if (lower < targetBC) {
        left = mid + 1;
      } else {
        right = mid;
      }
    }
    return left;
  }

  // Compress the temporary table into a single flat table
  GPUh() void flatten(const std::vector<TableEntry> table[NLayers])
  {
    // Count total entries
    size_t total{0};
    for (int32_t i{0}; i < NLayers; ++i) {
      total += table[i].size();
    }

    mFlatTable.reserve(total);

    // Build flat table and indices
    for (int32_t i{0}; i < NLayers; ++i) {
      mIndices[i].setFirstEntry(static_cast<T>(mFlatTable.size()));
      mIndices[i].setEntries(static_cast<T>(table[i].size()));
      mFlatTable.insert(mFlatTable.end(), table[i].begin(), table[i].end());
    }
  }

  TableIndex mIndices[NLayers];
  std::vector<TableEntry> mFlatTable;
};

// GPU-friendly view of the ROF mask table
template <int32_t NLayers, typename TableEntry, typename TableIndex>
using ROFMaskTableView = o2::itsmft::tracking::ROFMaskView<TableEntry, TableIndex>;

// Per-ROF per-layer boolean mask (uint8_t for GPU compatibility).
template <int32_t NLayers>
class ROFMaskTable : public LayerTimingBase<NLayers>
{
 public:
  using T = LayerTimingBase<NLayers>::T;
  using BCRange = dataformats::RangeReference<T, T>;
  using TableIndex = uint32_t;
  using TableEntry = uint8_t;
  using View = ROFMaskTableView<NLayers, TableEntry, TableIndex>;

  ROFMaskTable() = default;
  GPUh() explicit ROFMaskTable(const LayerTimingBase<NLayers>& timingBase) : LayerTimingBase<NLayers>(timingBase) { init(); }

  GPUh() void init()
  {
    int32_t totalROFs = 0;
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      mLayerROFOffsets[layer] = totalROFs;
      totalROFs += this->getLayer(layer).mNROFsTF;
    }
    mLayerROFOffsets[NLayers] = totalROFs; // sentinel
    mFlatMask.resize(totalROFs, 0u);
  }

  GPUh() size_t getFlatMaskSize() const noexcept { return mFlatMask.size(); }

  GPUh() void setROFEnabled(int32_t layer, int32_t rofId, uint8_t state = 1) noexcept
  {
    assert(layer >= 0 && layer < NLayers);
    assert(rofId >= 0 && rofId < mLayerROFOffsets[layer + 1] - mLayerROFOffsets[layer]);
    mFlatMask[mLayerROFOffsets[layer] + rofId] = state;
  }

  GPUh() void setROFsEnabled(int32_t layer, int32_t firstRof, int32_t nRofs, uint8_t state = 1) noexcept
  {
    assert(layer >= 0 && layer < NLayers);
    assert(firstRof >= 0);
    assert(firstRof + nRofs <= mLayerROFOffsets[layer + 1] - mLayerROFOffsets[layer]);
    std::memset(mFlatMask.data() + mLayerROFOffsets[layer] + firstRof, state, nRofs);
  }

  // Enable all ROFs in all layers that are time-compatible with the given BC range
  GPUh() void selectROF(const BCRange& t)
  {
    const int32_t bcStart = t.getFirstEntry();
    const int32_t bcEnd = t.getEntriesBound();
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      const auto& lay = this->getLayer(layer);
      const int32_t offset = mLayerROFOffsets[layer];
      for (int32_t rofId{0}; rofId < lay.mNROFsTF; ++rofId) {
        if (static_cast<int32_t>(lay.getROFStartInBC(rofId)) < bcEnd &&
            static_cast<int32_t>(lay.getROFEndInBC(rofId)) > bcStart) {
          mFlatMask[offset + rofId] = 1u;
        }
      }
    }
  }

  // Reset mask to 0, then enable all ROFs compatible with any of the given BC ranges
  GPUh() void selectROFs(const std::vector<BCRange>& ts)
  {
    resetMask();
    for (const auto& t : ts) {
      selectROF(t);
    }
  }

  GPUh() void resetMask(uint8_t s = 0u)
  {
    std::memset(mFlatMask.data(), s, mFlatMask.size());
  }

  GPUh() void invertMask()
  {
    std::ranges::transform(mFlatMask, mFlatMask.begin(), [](uint8_t x) { return 1 - x; });
  }

  GPUh() void swap(ROFMaskTable& other) noexcept
  {
    std::swap(mFlatMask, other.mFlatMask);
    std::swap(mLayerROFOffsets, other.mLayerROFOffsets);
  }

  GPUh() View getView() const
  {
    View view;
    view.mFlatMask = mFlatMask.data();
    view.mLayerROFOffsets = mLayerROFOffsets;
    view.mLayerCount = NLayers;
    return view;
  }

  GPUh() View getDeviceView(const TableEntry* deviceFlatMaskPtr, const TableIndex* deviceOffsetPtr) const
  {
    View view;
    view.mFlatMask = deviceFlatMaskPtr;
    view.mLayerROFOffsets = deviceOffsetPtr;
    view.mLayerCount = NLayers;
    return view;
  }

 private:
  TableIndex mLayerROFOffsets[NLayers + 1] = {0};
  std::vector<TableEntry> mFlatMask;
};

} // namespace o2::its

#endif
