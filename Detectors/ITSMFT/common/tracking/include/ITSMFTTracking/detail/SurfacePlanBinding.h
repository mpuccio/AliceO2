// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEPLANBINDING_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEPLANBINDING_H_

// Private plan binding used by the common tracker. It owns validated ordered
// surface positions and source-qualified compact transition/cell schedules;
// it owns no graph, workspace, event data, or detector-specific policy.

#ifndef GPUCA_GPUCODE

#include <algorithm>
#include <cstdint>
#include <memory>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <vector>

#include <gsl/span>

#include "ITSMFTTracking/SurfaceGraph.h"

namespace o2::itsmft::tracking
{

enum class SurfacePlanBindingError : uint8_t {
  None,
  InvalidSurfaceMask,
  InvalidLegacySurfaceOrder,
  InvalidTopology,
  InvalidSurface,
  CrossBoundaryTransition,
  CrossBoundaryCell,
  IncompleteTransitionMapping,
  IncompleteCellMapping,
  DuplicateScratchSlot
};

class SurfacePlanBinding
{
 public:
  struct BuildResult {
    std::unique_ptr<SurfacePlanBinding> binding{};
    SurfacePlanBindingError error{SurfacePlanBindingError::None};
    bool ok() const noexcept { return static_cast<bool>(binding); }
  };

  // `orderedSurfaces` maps each plan position to a global SurfaceId.
  static BuildResult build(const SurfaceGraphView& globalLayout,
                           SurfaceMask ownedSurfaces,
                           gsl::span<const SurfaceId> orderedSurfaces)
  {
    auto result = std::make_unique<SurfacePlanBinding>();
    if (globalLayout.surfaces == nullptr || globalLayout.nSurfaces == 0 ||
        !ownedSurfaces.isSubsetOf(SurfaceMask{globalLayout.nSurfaces == MaxLayoutSurfaces ? uint32_t{0xffffffff} : (uint32_t{1} << globalLayout.nSurfaces) - 1}) ||
        ownedSurfaces.count() != static_cast<int>(orderedSurfaces.size())) {
      return {{}, SurfacePlanBindingError::InvalidSurfaceMask};
    }

    result->mOwnedSurfaces = ownedSurfaces;
    // Retain the validated positional order; it is the runtime traversal
    // authority. The inverse map serves sparse global-ID lookups.
    result->mOrderedSurfaces.assign(orderedSurfaces.begin(), orderedSurfaces.end());
    result->mOwnedSurfaceIndexBySurface.assign(globalLayout.nSurfaces, -1);
    for (uint16_t position = 0; position < orderedSurfaces.size(); ++position) {
      const auto surface = orderedSurfaces[position];
      if (!surface.isValid() || surface.value() >= globalLayout.nSurfaces || !ownedSurfaces.has(surface) ||
          result->mOwnedSurfaceIndexBySurface[surface.value()] >= 0) {
        return {{}, SurfacePlanBindingError::InvalidLegacySurfaceOrder};
      }
      result->mOwnedSurfaceIndexBySurface[surface.value()] = static_cast<int16_t>(position);
    }
    for (uint16_t surface = 0; surface < globalLayout.nSurfaces; ++surface) {
      if (ownedSurfaces.has(SurfaceId{surface}) && result->mOwnedSurfaceIndexBySurface[surface] < 0) {
        return {{}, SurfacePlanBindingError::InvalidLegacySurfaceOrder};
      }
    }

    if ((globalLayout.nTransitions != 0 && globalLayout.transitions == nullptr) ||
        (globalLayout.nCells != 0 && globalLayout.cells == nullptr)) {
      return {{}, SurfacePlanBindingError::InvalidTopology};
    }
    const auto& topology = globalLayout;
    std::vector<uint32_t> indegree(topology.nSurfaces, 0);
    std::vector<uint32_t> rank(topology.nSurfaces, 0);
    for (uint32_t id = 0; id < topology.nTransitions; ++id) {
      const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(id)});
      if (!transition.from.isValid() || !transition.to.isValid() ||
          transition.from.value() >= topology.nSurfaces || transition.to.value() >= topology.nSurfaces) {
        return {{}, SurfacePlanBindingError::InvalidTopology};
      }
      ++indegree[transition.to.value()];
    }
    const auto laterSurface = [](SurfaceId lhs, SurfaceId rhs) { return rhs < lhs; };
    std::priority_queue<SurfaceId, std::vector<SurfaceId>, decltype(laterSurface)> ready{laterSurface};
    for (uint16_t surface = 0; surface < topology.nSurfaces; ++surface) {
      if (indegree[surface] == 0) {
        ready.push(SurfaceId{surface});
      }
    }
    uint32_t visited = 0;
    while (!ready.empty()) {
      const auto surface = ready.top();
      ready.pop();
      ++visited;
      for (uint32_t id = 0; id < topology.nTransitions; ++id) {
        const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(id)});
        if (transition.from != surface) {
          continue;
        }
        rank[transition.to.value()] = std::max(rank[transition.to.value()], rank[surface.value()] + 1);
        if (--indegree[transition.to.value()] == 0) {
          ready.push(transition.to);
        }
      }
    }
    if (visited != topology.nSurfaces) {
      return {{}, SurfacePlanBindingError::InvalidTopology};
    }
    std::vector<uint16_t> componentParent(topology.nSurfaces);
    std::iota(componentParent.begin(), componentParent.end(), uint16_t{0});
    const auto componentRoot = [&componentParent](uint16_t surface) {
      while (componentParent[surface] != surface) {
        componentParent[surface] = componentParent[componentParent[surface]];
        surface = componentParent[surface];
      }
      return surface;
    };
    for (uint32_t id = 0; id < topology.nTransitions; ++id) {
      const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(id)});
      const auto fromRoot = componentRoot(transition.from.value());
      const auto toRoot = componentRoot(transition.to.value());
      if (fromRoot != toRoot) {
        componentParent[toRoot] = fromRoot;
      }
    }
    std::vector<uint32_t> componentOrder(topology.nSurfaces, std::numeric_limits<uint32_t>::max());
    for (uint32_t position = 0; position < orderedSurfaces.size(); ++position) {
      const auto root = componentRoot(orderedSurfaces[position].value());
      componentOrder[root] = std::min(componentOrder[root], position);
    }
    for (uint32_t id = 0; id < topology.nCells; ++id) {
      const auto& cell = topology.getCell(CellTopologyId{static_cast<uint16_t>(id)});
      if (!cell.firstTransition.isValid() || !cell.secondTransition.isValid() ||
          cell.firstTransition.value() >= topology.nTransitions || cell.secondTransition.value() >= topology.nTransitions) {
        return {{}, SurfacePlanBindingError::InvalidTopology};
      }
    }
    result->mScratchTransitionSlot.assign(topology.nTransitions, -1);
    result->mScratchCellSlot.assign(topology.nCells, -1);

    // Validate every global transition first. A boundary edge is invalid
    // regardless of whether it would otherwise be filtered from this
    // binding's own scope.
    for (uint32_t id = 0; id < topology.nTransitions; ++id) {
      const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(id)});
      const bool fromOwned = ownedSurfaces.has(transition.from);
      const bool toOwned = ownedSurfaces.has(transition.to);
      if (fromOwned != toOwned) {
        return {{}, SurfacePlanBindingError::CrossBoundaryTransition};
      }
      if (!fromOwned) {
        continue;
      }
      if (!transition.skippedSurfaces.isSubsetOf(ownedSurfaces)) {
        return {{}, SurfacePlanBindingError::CrossBoundaryTransition};
      }
    }

    // Filter the immutable global-ID order only by ownership; disconnected
    // SurfaceKind components remain part of this one binding and are batched
    // by TrackerTraits from their endpoint descriptors.
    for (uint32_t rawId = 0; rawId < topology.nTransitions; ++rawId) {
      const auto id = TransitionId{static_cast<uint16_t>(rawId)};
      const auto& transition = topology.getTransition(id);
      if (!ownedSurfaces.has(transition.from)) {
        continue;
      }
      if (result->mScratchTransitionSlot[id.value()] >= 0) {
        return {{}, SurfacePlanBindingError::DuplicateScratchSlot};
      }
      result->mScratchTransitionSlot[id.value()] = static_cast<int16_t>(result->mGlobalTransitions.size());
      result->mGlobalTransitions.push_back(id);
    }
    for (uint32_t id = 0; id < topology.nTransitions; ++id) {
      const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(id)});
      if (ownedSurfaces.has(transition.from) && result->mScratchTransitionSlot[id] < 0) {
        return {{}, SurfacePlanBindingError::IncompleteTransitionMapping};
      }
    }

    for (uint32_t id = 0; id < topology.nCells; ++id) {
      const auto cellId = CellTopologyId{static_cast<uint16_t>(id)};
      const auto& cell = topology.getCell(cellId);
      const bool firstOwned = cell.firstTransition.isValid() && cell.firstTransition.value() < topology.nTransitions &&
                              result->mScratchTransitionSlot[cell.firstTransition.value()] >= 0;
      const bool secondOwned = cell.secondTransition.isValid() && cell.secondTransition.value() < topology.nTransitions &&
                               result->mScratchTransitionSlot[cell.secondTransition.value()] >= 0;
      if (firstOwned != secondOwned || (firstOwned && !cell.hitSurfaces.isSubsetOf(ownedSurfaces))) {
        return {{}, SurfacePlanBindingError::CrossBoundaryCell};
      }
    }
    for (uint32_t rawId = 0; rawId < topology.nCells; ++rawId) {
      const auto id = CellTopologyId{static_cast<uint16_t>(rawId)};
      const auto& cell = topology.getCell(id);
      if (result->mScratchTransitionSlot[cell.firstTransition.value()] < 0) {
        continue;
      }
      if (result->mScratchCellSlot[id.value()] >= 0) {
        return {{}, SurfacePlanBindingError::DuplicateScratchSlot};
      }
      result->mScratchCellSlot[id.value()] = static_cast<int16_t>(result->mGlobalCells.size());
      result->mGlobalCells.push_back(id);
    }
    for (uint32_t id = 0; id < topology.nCells; ++id) {
      const auto& cell = topology.getCell(CellTopologyId{static_cast<uint16_t>(id)});
      if (result->mScratchTransitionSlot[cell.firstTransition.value()] >= 0 && result->mScratchCellSlot[id] < 0) {
        return {{}, SurfacePlanBindingError::IncompleteCellMapping};
      }
    }
    for (uint32_t rawId = 0; rawId < topology.nCells; ++rawId) {
      const auto id = CellTopologyId{static_cast<uint16_t>(rawId)};
      const auto& cell = topology.getCell(id);
      if (result->mScratchCellSlot[id.value()] >= 0 &&
          topology.seedingSurfaces.has(topology.getTransition(cell.secondTransition).to)) {
        result->mGlobalRoadStartCells.push_back(id);
      }
    }
    std::sort(result->mGlobalRoadStartCells.begin(), result->mGlobalRoadStartCells.end(), [&](CellTopologyId lhs, CellTopologyId rhs) {
      const auto lhsTarget = topology.getTransition(topology.getCell(lhs).secondTransition).to;
      const auto rhsTarget = topology.getTransition(topology.getCell(rhs).secondTransition).to;
      const auto lhsComponent = componentOrder[componentRoot(lhsTarget.value())];
      const auto rhsComponent = componentOrder[componentRoot(rhsTarget.value())];
      if (lhsComponent != rhsComponent) {
        return lhsComponent < rhsComponent;
      }
      return rank[lhsTarget.value()] != rank[rhsTarget.value()] ? rank[lhsTarget.value()] < rank[rhsTarget.value()] : lhs < rhs;
    });
    result->mRoadStartComponentOffsets.push_back(0);
    uint32_t previousComponent = std::numeric_limits<uint32_t>::max();
    for (uint32_t offset = 0; offset < result->mGlobalRoadStartCells.size(); ++offset) {
      const auto cell = result->mGlobalRoadStartCells[offset];
      const auto target = topology.getTransition(topology.getCell(cell).secondTransition).to;
      const auto component = componentOrder[componentRoot(target.value())];
      if (component != previousComponent && offset != 0) {
        result->mRoadStartComponentOffsets.push_back(offset);
      }
      previousComponent = component;
    }
    result->mRoadStartComponentOffsets.push_back(static_cast<uint32_t>(result->mGlobalRoadStartCells.size()));
    // Same source cells as mGlobalCells above, but in rank-sorted
    // (neighbour-schedule) order rather than ascending CellTopologyId order.
    // Every id here is already confirmed owned (mScratchCellSlot valid) by
    // the mGlobalCells loop above, so this is a pure reordering of the same
    // set, never a superset/subset of it.
    std::vector<CellTopologyId> scheduledCells = result->mGlobalCells;
    std::sort(scheduledCells.begin(), scheduledCells.end(), [&](CellTopologyId lhs, CellTopologyId rhs) {
      const auto lhsTarget = topology.getTransition(topology.getCell(lhs).secondTransition).to;
      const auto rhsTarget = topology.getTransition(topology.getCell(rhs).secondTransition).to;
      const auto lhsComponent = componentOrder[componentRoot(lhsTarget.value())];
      const auto rhsComponent = componentOrder[componentRoot(rhsTarget.value())];
      if (lhsComponent != rhsComponent) {
        return lhsComponent < rhsComponent;
      }
      return rank[lhsTarget.value()] != rank[rhsTarget.value()] ? rank[lhsTarget.value()] < rank[rhsTarget.value()] : lhs < rhs;
    });
    for (const auto id : scheduledCells) {
      if (result->mScratchCellSlot[id.value()] >= 0) {
        result->mGlobalScheduledCells.push_back(id);
      }
    }
    return {std::move(result), SurfacePlanBindingError::None};
  }

  SurfaceMask getOwnedSurfaces() const noexcept { return mOwnedSurfaces; }
  gsl::span<const SurfaceId> getOrderedSurfaces() const noexcept { return mOrderedSurfaces; }
  std::optional<uint16_t> getOwnedSurfaceIndex(SurfaceId id) const noexcept { return getSlot(mOwnedSurfaceIndexBySurface, id); }
  std::optional<uint16_t> getScratchTransitionSlot(TransitionId id) const noexcept { return getSlot(mScratchTransitionSlot, id); }
  std::optional<uint16_t> getScratchCellSlot(CellTopologyId id) const noexcept { return getSlot(mScratchCellSlot, id); }
  gsl::span<const TransitionId> getGlobalTransitions() const noexcept { return mGlobalTransitions; }
  gsl::span<const CellTopologyId> getGlobalCells() const noexcept { return mGlobalCells; }
  gsl::span<const CellTopologyId> getGlobalRoadStartCells() const noexcept { return mGlobalRoadStartCells; }
  gsl::span<const uint32_t> getRoadStartComponentOffsets() const noexcept { return mRoadStartComponentOffsets; }
  gsl::span<const CellTopologyId> getGlobalScheduledCells() const noexcept { return mGlobalScheduledCells; }

 private:
  template <typename Id>
  static std::optional<uint16_t> getSlot(const std::vector<int16_t>& map, Id id) noexcept
  {
    if (!id.isValid() || id.value() >= map.size() || map[id.value()] < 0) {
      return std::nullopt;
    }
    return static_cast<uint16_t>(map[id.value()]);
  }

  SurfaceMask mOwnedSurfaces{};
  std::vector<SurfaceId> mOrderedSurfaces;
  std::vector<int16_t> mOwnedSurfaceIndexBySurface;
  std::vector<int16_t> mScratchTransitionSlot;
  std::vector<int16_t> mScratchCellSlot;
  std::vector<TransitionId> mGlobalTransitions;
  std::vector<CellTopologyId> mGlobalCells;
  std::vector<CellTopologyId> mGlobalRoadStartCells;
  std::vector<uint32_t> mRoadStartComponentOffsets;
  std::vector<CellTopologyId> mGlobalScheduledCells;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_SURFACEPLANBINDING_H_
