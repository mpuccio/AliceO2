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

#include <array>
#include <vector>

#include <gsl/gsl>

namespace o2::itsmft::tracking
{

/// Host-side, one-shot grouping of an already-validated DetectorLayoutView's
/// transitions and cell topologies by TransitionPolicyTag. Built once outside
/// candidate/hot loops; `dispatchTransitionPolicies` then issues at most one
/// call per active tag into a template-specialized policy implementation.
class TransitionPolicyGrouping
{
 public:
  explicit TransitionPolicyGrouping(const DetectorLayoutView& layout)
  {
    if (layout.surfaces == nullptr) {
      return;
    }
    const auto& topology = layout.topology;
    for (uint32_t i = 0; i < topology.nTransitions; ++i) {
      const auto id = TransitionId{static_cast<uint16_t>(i)};
      auto* group = groupFor(topology.getTransition(id).policyTag);
      if (group != nullptr) {
        group->transitions.push_back(id);
      }
    }
    for (uint32_t i = 0; i < topology.nCells; ++i) {
      const auto id = CellTopologyId{static_cast<uint16_t>(i)};
      const auto& cell = topology.getCell(id);
      auto* group = groupFor(topology.getTransition(cell.firstTransition).policyTag);
      if (group != nullptr) {
        group->cells.push_back(id);
      }
    }
  }

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

 private:
  struct Group {
    std::vector<TransitionId> transitions;
    std::vector<CellTopologyId> cells;
  };

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
