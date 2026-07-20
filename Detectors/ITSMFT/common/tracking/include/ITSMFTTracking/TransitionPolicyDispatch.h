// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYDISPATCH_H_
#define ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYDISPATCH_H_

#include <cstdint>

#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/TransitionPolicy.h"
#include "ITSMFTTracking/TransitionPolicyState.h"

// Host-only: this is the outer-loop dispatch boundary itself, never a
// candidate/hot-loop body, so it has no device-side counterpart (D007/
// Architecture.md 10.1). GPU execution instead launches one specialized
// kernel per active (stage, family), selected from this same grouping.
#ifndef GPUCA_GPUCODE

#include <algorithm>
#include <array>
#include <queue>
#include <vector>

#include <gsl/gsl>

namespace o2::itsmft::tracking
{

enum class TransitionPolicyScheduleError : uint8_t {
  None,
  MissingSurfaceData,
  MissingTopologyData,
  InvalidTransitionSurface,
  InvalidCellTransition,
  CyclicTopology
};

/// Host-side, one-shot grouping of an already-validated DetectorLayoutView's
/// transitions and cell topologies by TransitionPolicyTag. Built once outside
/// candidate/hot loops; `dispatchTransitionPolicies` then issues at most one
/// call per active tag into a template-specialized policy implementation.
class TransitionPolicyGrouping
{
 public:
  explicit TransitionPolicyGrouping(const DetectorLayoutView& layout)
  {
    if (layout.nSurfaces != 0 && layout.surfaces == nullptr) {
      mScheduleError = TransitionPolicyScheduleError::MissingSurfaceData;
      return;
    }
    const auto& topology = layout.topology;
    if ((topology.nTransitions != 0 && topology.transitions == nullptr) ||
        (topology.nCells != 0 && topology.cells == nullptr)) {
      mScheduleError = TransitionPolicyScheduleError::MissingTopologyData;
      return;
    }

    std::vector<uint32_t> indegree(layout.nSurfaces, 0);
    std::vector<uint32_t> rank(layout.nSurfaces, 0);
    for (uint32_t i = 0; i < topology.nTransitions; ++i) {
      const auto id = TransitionId{static_cast<uint16_t>(i)};
      const auto& transition = topology.getTransition(id);
      if (!transition.from.isValid() || !transition.to.isValid() ||
          transition.from.value() >= layout.nSurfaces || transition.to.value() >= layout.nSurfaces) {
        mScheduleError = TransitionPolicyScheduleError::InvalidTransitionSurface;
        clear();
        return;
      }
      ++indegree[transition.to.value()];
      auto* group = groupFor(transition.policyTag);
      if (group != nullptr) {
        group->transitions.push_back(id);
      }
    }

    const auto laterSurface = [](SurfaceId lhs, SurfaceId rhs) { return rhs < lhs; };
    std::priority_queue<SurfaceId, std::vector<SurfaceId>, decltype(laterSurface)> ready{laterSurface};
    for (uint16_t surface = 0; surface < layout.nSurfaces; ++surface) {
      if (indegree[surface] == 0) {
        ready.push(SurfaceId{surface});
      }
    }
    uint32_t visited = 0;
    while (!ready.empty()) {
      const auto surface = ready.top();
      ready.pop();
      ++visited;
      for (uint32_t i = 0; i < topology.nTransitions; ++i) {
        const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(i)});
        if (transition.from != surface) {
          continue;
        }
        rank[transition.to.value()] = std::max(rank[transition.to.value()], rank[surface.value()] + 1);
        if (--indegree[transition.to.value()] == 0) {
          ready.push(transition.to);
        }
      }
    }
    if (visited != layout.nSurfaces) {
      mScheduleError = TransitionPolicyScheduleError::CyclicTopology;
      clear();
      return;
    }

    for (uint32_t i = 0; i < topology.nCells; ++i) {
      const auto id = CellTopologyId{static_cast<uint16_t>(i)};
      const auto& cell = topology.getCell(id);
      if (!cell.firstTransition.isValid() || !cell.secondTransition.isValid() ||
          cell.firstTransition.value() >= topology.nTransitions || cell.secondTransition.value() >= topology.nTransitions) {
        mScheduleError = TransitionPolicyScheduleError::InvalidCellTransition;
        clear();
        return;
      }
      auto* group = groupFor(topology.getTransition(cell.firstTransition).policyTag);
      if (group != nullptr) {
        group->cells.push_back(id);
        group->scheduledCells.push_back(id);
        // Road-start eligibility (Architecture.md Sec 10, D003): a cell may
        // start a road iff its traversal endpoint -- the `to` SurfaceId of
        // its *second* transition, never a numeric/highest-bit reading of
        // hitSurfaces -- is one of the layout's seeding surfaces. Both
        // `endpoint` (SurfaceId) and `topology.seedingSurfaces` (SurfaceMask)
        // live in the global-SurfaceId space; this is a SurfaceMask::has()
        // test, never a legacy vector index. Appended in the same ascending
        // CellTopologyId order as `cells` above, so `roadStartCells` is
        // exactly that loop's eligible subsequence -- not sorted, not
        // rank-ordered like `scheduledCells` (which stays neighbour-schedule
        // specific and must not be reused here).
        const auto endpoint = topology.getTransition(cell.secondTransition).to;
        if (topology.seedingSurfaces.has(endpoint)) {
          group->roadStartCells.push_back(id);
        }
      }
    }
    for (auto& group : mGroups) {
      std::sort(group.scheduledCells.begin(), group.scheduledCells.end(), [&](CellTopologyId lhs, CellTopologyId rhs) {
        const auto lhsTarget = topology.getTransition(topology.getCell(lhs).secondTransition).to;
        const auto rhsTarget = topology.getTransition(topology.getCell(rhs).secondTransition).to;
        return rank[lhsTarget.value()] != rank[rhsTarget.value()] ? rank[lhsTarget.value()] < rank[rhsTarget.value()] : lhs < rhs;
      });
    }
  }

  bool valid() const noexcept { return mScheduleError == TransitionPolicyScheduleError::None; }
  TransitionPolicyScheduleError getScheduleError() const noexcept { return mScheduleError; }

  bool hasTag(TransitionPolicyTag tag) const noexcept
  {
    const auto* group = findGroup(tag);
    return group != nullptr && !group->transitions.empty();
  }

  gsl::span<const TransitionId> transitionsForTag(TransitionPolicyTag tag) const noexcept
  {
    const auto* group = findGroup(tag);
    return group != nullptr ? gsl::span<const TransitionId>(group->transitions) : gsl::span<const TransitionId>();
  }

  gsl::span<const CellTopologyId> cellsForTag(TransitionPolicyTag tag) const noexcept
  {
    const auto* group = findGroup(tag);
    return group != nullptr ? gsl::span<const CellTopologyId>(group->cells) : gsl::span<const CellTopologyId>();
  }

  gsl::span<const CellTopologyId> scheduledCellsForTag(TransitionPolicyTag tag) const noexcept
  {
    const auto* group = findGroup(tag);
    return group != nullptr ? gsl::span<const CellTopologyId>(group->scheduledCells) : gsl::span<const CellTopologyId>();
  }

  /// Deterministic, ascending-CellTopologyId subsequence of `cellsForTag`
  /// whose traversal endpoint (topology.getTransition(cell.secondTransition).to)
  /// is a layout seeding surface. Built once per grouping construction
  /// (TrackerTraits::initialiseTimeFrame(), Architecture.md Sec 10.1); never
  /// rank-sorted, unlike `scheduledCellsForTag`. An empty seeding mask, or a
  /// seeding surface that terminates no cell, both yield a valid empty/short
  /// span -- not a schedule error.
  gsl::span<const CellTopologyId> roadStartCellsForTag(TransitionPolicyTag tag) const noexcept
  {
    const auto* group = findGroup(tag);
    return group != nullptr ? gsl::span<const CellTopologyId>(group->roadStartCells) : gsl::span<const CellTopologyId>();
  }

 private:
  struct Group {
    std::vector<TransitionId> transitions;
    std::vector<CellTopologyId> cells;
    std::vector<CellTopologyId> scheduledCells;
    std::vector<CellTopologyId> roadStartCells;
  };

  void clear() noexcept
  {
    for (auto& group : mGroups) {
      group.transitions.clear();
      group.cells.clear();
      group.scheduledCells.clear();
      group.roadStartCells.clear();
    }
  }

  static constexpr int indexOf(TransitionPolicyTag tag) noexcept
  {
    switch (tag) {
      case TransitionPolicyTag::CylinderCylinder:
        return 0;
      case TransitionPolicyTag::DiskDisk:
        return 1;
      case TransitionPolicyTag::Invalid:
        return -1;
    }
    return -1;
  }

  Group* groupFor(TransitionPolicyTag tag) noexcept
  {
    const auto idx = indexOf(tag);
    return idx >= 0 ? &mGroups[idx] : nullptr;
  }
  const Group* findGroup(TransitionPolicyTag tag) const noexcept
  {
    const auto idx = indexOf(tag);
    return idx >= 0 ? &mGroups[idx] : nullptr;
  }

  std::array<Group, 2> mGroups{};
  TransitionPolicyScheduleError mScheduleError{TransitionPolicyScheduleError::None};
};

/// Outer-loop dispatch: for each TransitionPolicyTag with active work, calls
/// `visitor(TransitionPolicyTraits<Tag>{}, transitionIds, cellIds)` exactly
/// once. `visitor` must be a generic callable (templated on the traits type);
/// its body selects a template-specialized policy implementation from
/// `Traits::Family` / `Traits::Params` / `Traits::SeedState` and must not
/// branch on detector ID. A disconnected layout combining both Stage-A
/// families therefore triggers exactly two calls; a single-family layout
/// triggers exactly one. This function itself must never be invoked from
/// inside a candidate, cluster, neighbour, road, or refit loop.
template <typename Visitor>
void dispatchTransitionPolicies(const TransitionPolicyGrouping& grouping, Visitor&& visitor)
{
  if (grouping.hasTag(TransitionPolicyTag::CylinderCylinder)) {
    visitor(TransitionPolicyTraits<TransitionPolicyTag::CylinderCylinder>{},
            grouping.transitionsForTag(TransitionPolicyTag::CylinderCylinder),
            grouping.cellsForTag(TransitionPolicyTag::CylinderCylinder));
  }
  if (grouping.hasTag(TransitionPolicyTag::DiskDisk)) {
    visitor(TransitionPolicyTraits<TransitionPolicyTag::DiskDisk>{},
            grouping.transitionsForTag(TransitionPolicyTag::DiskDisk),
            grouping.cellsForTag(TransitionPolicyTag::DiskDisk));
  }
}

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYDISPATCH_H_ */
