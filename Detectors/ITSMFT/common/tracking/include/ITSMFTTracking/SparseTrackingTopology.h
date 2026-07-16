// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SPARSETRACKINGTOPOLOGY_H_
#define ALICEO2_ITSMFT_TRACKING_SPARSETRACKINGTOPOLOGY_H_

#include <cstdint>
#include <type_traits>

#ifndef GPUCA_GPUCODE
#include <algorithm>
#include <vector>
#endif

#include "GPUCommonDef.h"
#include "ITSMFTTracking/SurfaceId.h"
#include "ITSMFTTracking/SurfaceMask.h"

namespace o2::itsmft::tracking
{

struct SurfaceTransition {
  SurfaceId from{};
  SurfaceId to{};
  SurfaceMask skippedSurfaces{};
  uint16_t policyTag{0};
  uint16_t flags{0};
};

struct SurfaceCellTopology {
  TransitionId firstTransition{};
  TransitionId secondTransition{};
  SurfaceMask hitSurfaces{};
};

struct TopologyRange {
  uint32_t firstEntry{0};
  uint32_t entries{0};

  GPUhdi() constexpr uint32_t getFirstEntry() const noexcept { return firstEntry; }
  GPUhdi() constexpr uint32_t getEntries() const noexcept { return entries; }
  GPUhdi() constexpr uint32_t getEntriesBound() const noexcept { return firstEntry + entries; }
};

struct SparseTrackingTopologyView {
  const SurfaceTransition* transitions{nullptr};
  const SurfaceCellTopology* cells{nullptr};
  const uint32_t* cellsByFirstTransitionOffsets{nullptr};
  const CellTopologyId* cellsByFirstTransition{nullptr};
  SurfaceMask seedingSurfaces{};
  uint32_t nTransitions{0};
  uint32_t nCells{0};

  GPUhdi() const SurfaceTransition& getTransition(TransitionId id) const { return transitions[id.value()]; }
  GPUhdi() const SurfaceCellTopology& getCell(CellTopologyId id) const { return cells[id.value()]; }
  GPUhdi() TopologyRange getCellsStartingWithTransition(TransitionId transition) const
  {
    const uint32_t id = transition.value();
    return TopologyRange{cellsByFirstTransitionOffsets[id], cellsByFirstTransitionOffsets[id + 1] - cellsByFirstTransitionOffsets[id]};
  }
};

static_assert(std::is_standard_layout_v<SurfaceTransition> && std::is_trivially_copyable_v<SurfaceTransition>);
static_assert(std::is_standard_layout_v<SurfaceCellTopology> && std::is_trivially_copyable_v<SurfaceCellTopology>);
static_assert(std::is_standard_layout_v<TopologyRange> && std::is_trivially_copyable_v<TopologyRange>);
static_assert(std::is_standard_layout_v<SparseTrackingTopologyView> && std::is_trivially_copyable_v<SparseTrackingTopologyView>);
static_assert(sizeof(SurfaceTransition) == 12);
static_assert(sizeof(SurfaceCellTopology) == 8);

#ifndef GPUCA_GPUCODE

enum class TopologyBuildError : uint8_t {
  None,
  InvalidSurfaceCount,
  InvalidSurface,
  SelfTransition,
  DuplicateTransition,
  TooManyTransitions,
  InvalidTransition,
  DisconnectedTransitions,
  RepeatedSurface,
  DuplicateCell,
  TooManyCells,
  AlreadyFinalized,
  NotFinalized
};

class SparseTrackingTopology
{
 public:
  explicit SparseTrackingTopology(uint32_t surfaceCount, SurfaceMask seedingSurfaces = {})
    : mSurfaceCount{surfaceCount}, mSeedingSurfaces{seedingSurfaces}
  {
    if (surfaceCount > MaxLayoutSurfaces) {
      mError = TopologyBuildError::InvalidSurfaceCount;
    } else if (!seedingSurfaces.isSubsetOf(activeSurfaceMask())) {
      mError = TopologyBuildError::InvalidSurface;
    }
  }

  TransitionId addTransition(SurfaceTransition transition)
  {
    if (!canModify()) {
      return TransitionId::invalid();
    }
    if (!isSurfaceInLayout(transition.from) || !isSurfaceInLayout(transition.to) ||
        !transition.skippedSurfaces.isSubsetOf(activeSurfaceMask()) ||
        transition.skippedSurfaces.has(transition.from) || transition.skippedSurfaces.has(transition.to)) {
      mError = TopologyBuildError::InvalidSurface;
      return TransitionId::invalid();
    }
    if (transition.from == transition.to) {
      mError = TopologyBuildError::SelfTransition;
      return TransitionId::invalid();
    }
    const auto duplicate = std::find_if(mTransitions.begin(), mTransitions.end(), [&](const auto& existing) {
      return existing.from == transition.from && existing.to == transition.to;
    });
    if (duplicate != mTransitions.end()) {
      mError = TopologyBuildError::DuplicateTransition;
      return TransitionId::invalid();
    }
    if (mTransitions.size() >= MaxLayoutTransitions) {
      mError = TopologyBuildError::TooManyTransitions;
      return TransitionId::invalid();
    }
    const auto id = TransitionId{static_cast<uint16_t>(mTransitions.size())};
    mTransitions.push_back(transition);
    return id;
  }

  CellTopologyId addCell(TransitionId first, TransitionId second)
  {
    if (!canModify()) {
      return CellTopologyId::invalid();
    }
    if (!isTransitionInLayout(first) || !isTransitionInLayout(second)) {
      mError = TopologyBuildError::InvalidTransition;
      return CellTopologyId::invalid();
    }
    const auto& firstTransition = mTransitions[first.value()];
    const auto& secondTransition = mTransitions[second.value()];
    if (firstTransition.to != secondTransition.from) {
      mError = TopologyBuildError::DisconnectedTransitions;
      return CellTopologyId::invalid();
    }
    if (firstTransition.from == secondTransition.to) {
      mError = TopologyBuildError::RepeatedSurface;
      return CellTopologyId::invalid();
    }
    const auto duplicate = std::find_if(mCells.begin(), mCells.end(), [&](const auto& existing) {
      return existing.firstTransition == first && existing.secondTransition == second;
    });
    if (duplicate != mCells.end()) {
      mError = TopologyBuildError::DuplicateCell;
      return CellTopologyId::invalid();
    }
    if (mCells.size() >= MaxLayoutCellTopologies) {
      mError = TopologyBuildError::TooManyCells;
      return CellTopologyId::invalid();
    }

    SurfaceMask hits;
    hits.set(firstTransition.from);
    hits.set(firstTransition.to);
    hits.set(secondTransition.to);
    const auto id = CellTopologyId{static_cast<uint16_t>(mCells.size())};
    mCells.push_back(SurfaceCellTopology{first, second, hits});
    return id;
  }

  bool finalize()
  {
    if (mFinalized) {
      return true;
    }
    if (mError != TopologyBuildError::None) {
      return false;
    }
    mCellsByFirstTransitionOffsets.assign(mTransitions.size() + 1, 0);
    for (const auto& cell : mCells) {
      ++mCellsByFirstTransitionOffsets[cell.firstTransition.value() + 1];
    }
    for (size_t i = 1; i < mCellsByFirstTransitionOffsets.size(); ++i) {
      mCellsByFirstTransitionOffsets[i] += mCellsByFirstTransitionOffsets[i - 1];
    }
    mCellsByFirstTransition.resize(mCells.size());
    auto cursor = mCellsByFirstTransitionOffsets;
    for (uint32_t cell = 0; cell < mCells.size(); ++cell) {
      const auto transition = mCells[cell].firstTransition.value();
      mCellsByFirstTransition[cursor[transition]++] = CellTopologyId{static_cast<uint16_t>(cell)};
    }
    mFinalized = true;
    return true;
  }

  SparseTrackingTopologyView getView() const noexcept
  {
    if (!mFinalized) {
      return {};
    }
    return getDeviceView(mTransitions.data(), mCells.data(), mCellsByFirstTransitionOffsets.data(), mCellsByFirstTransition.data());
  }

  SparseTrackingTopologyView getDeviceView(const SurfaceTransition* transitions,
                                           const SurfaceCellTopology* cells,
                                           const uint32_t* offsets,
                                           const CellTopologyId* cellsByFirstTransition) const noexcept
  {
    if (!mFinalized) {
      return {};
    }
    return SparseTrackingTopologyView{transitions, cells, offsets, cellsByFirstTransition, mSeedingSurfaces,
                                      static_cast<uint32_t>(mTransitions.size()), static_cast<uint32_t>(mCells.size())};
  }

  uint32_t getSurfaceCount() const noexcept { return mSurfaceCount; }
  TopologyBuildError getError() const noexcept { return mError; }
  bool isFinalized() const noexcept { return mFinalized; }
  const auto& getTransitions() const noexcept { return mTransitions; }
  const auto& getCells() const noexcept { return mCells; }
  const auto& getCellsByFirstTransitionOffsets() const noexcept { return mCellsByFirstTransitionOffsets; }
  const auto& getCellsByFirstTransition() const noexcept { return mCellsByFirstTransition; }

 private:
  bool canModify()
  {
    if (mFinalized) {
      mError = TopologyBuildError::AlreadyFinalized;
      return false;
    }
    return mError == TopologyBuildError::None;
  }
  bool isSurfaceInLayout(SurfaceId id) const noexcept { return id.isValid() && id.value() < mSurfaceCount; }
  bool isTransitionInLayout(TransitionId id) const noexcept { return id.isValid() && id.value() < mTransitions.size(); }
  SurfaceMask activeSurfaceMask() const noexcept
  {
    return SurfaceMask{mSurfaceCount == MaxLayoutSurfaces ? uint32_t{0xffffffff} : ((uint32_t{1} << mSurfaceCount) - 1)};
  }

  uint32_t mSurfaceCount{0};
  SurfaceMask mSeedingSurfaces{};
  TopologyBuildError mError{TopologyBuildError::None};
  bool mFinalized{false};
  std::vector<SurfaceTransition> mTransitions;
  std::vector<SurfaceCellTopology> mCells;
  std::vector<uint32_t> mCellsByFirstTransitionOffsets;
  std::vector<CellTopologyId> mCellsByFirstTransition;
};

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SPARSETRACKINGTOPOLOGY_H_ */
