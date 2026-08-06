// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEPLANBINDING_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEPLANBINDING_H_

// M6f (Detectors/ITSMFT/common/tracking/doc/design/0002-m6-generic-workspace-
// migration.md Sec 3.2, 7, 9; GenericTrackingEngineMigration.md M6): the one
// detector-neutral plan binding used by the common tracker. It replaced the
// temporary traversal bridge when M6d/M6e wiring completed; no second binding
// implementation remains in common production code. Confined to detail/ --
// private hot-loop-dispatch implementation, never a public/adapter-facing
// contract, like every other header under this directory (see
// TransitionPolicy.h's own file-level doc):
// TransitionPolicyTag must never be nameable outside this boundary.
//
// The two behavioral changes from the former traversal binding::build() that
// the design note specifies:
// 1. no detector-identity allow-list (the old `detector != ITS && detector
//    != MFT` gate is gone -- ownership/topology validation below already
//    fully determines correctness without it);
// 2. `expectedKind`/`expectedPolicy` are caller-supplied parameters instead
//    of being derived internally from a `detector` switch. They are no
//    longer guaranteed mutually consistent by construction, so build() adds
//    one new preflight check for that -- IncompatibleExpectedPolicyKind
//    below -- using the same isSurfaceKindCompatible() the rest of this
//    library already uses as its single shared policy/surface-kind
//    compatibility rule (TransitionPolicy.h).
//
// A third `detector` use the design note's "only two lines" framing did not
// separately enumerate -- the retired traversal binding::build() also checked
// every legacy-ordered surface's SurfaceDescriptor::detectorId against the
// caller-supplied `detector` (SurfaceDetectorMismatch) -- is deliberately
// dropped, not replaced. An earlier revision of this type reintroduced an
// equivalent check (every owned surface sharing one consistent detectorId)
// without an external parameter; that was itself a hidden constraint this
// type must not carry: SurfacePlanBinding is generic over its own
// SurfaceId set and must not assume "one binding, one detector" at all, so
// a future participant whose owned surfaces legitimately span more than one
// detectorId (e.g. a merged/aggregate plan) is not artificially rejected.
// detectorId-based ownership bookkeeping, if a caller ever needs it, belongs
// in that caller (which already knows its own semantics for "owner"), never
// in this detector-neutral type.

#ifndef GPUCA_GPUCODE

#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

#include <gsl/span>

#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/detail/TransitionPolicyDispatch.h"

namespace o2::itsmft::tracking
{

enum class SurfacePlanBindingError : uint8_t {
  None,
  InvalidSource,
  IncompatibleExpectedPolicyKind,
  InvalidSurfaceMask,
  InvalidLegacySurfaceOrder,
  InvalidTopology,
  InvalidPolicySurface,
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

  // `orderedSurfaces` maps position i to a global SurfaceId -- the only
  // global-to-positional projection here, exactly as
  // the retired traversal binding::build()'s own `legacySurfaceOrder` parameter.
  // `expectedKind`/`expectedPolicy` are the caller's own already-derived
  // family selection (e.g. from its own plan's surface kinds, mirroring
  // M4b's stateFamilyOf(SurfaceKind)/decision 8's
  // transitionPolicyTagForSurfaceKind()) -- never derived here from a
  // detector identity, since this type accepts none.
  static BuildResult build(const SurfaceGraphView& globalLayout,
                           ClusterSourceId source,
                           SurfaceMask ownedSurfaces,
                           gsl::span<const SurfaceId> orderedSurfaces,
                           SurfaceKind expectedKind,
                           TransitionPolicyTag expectedPolicy)
  {
    auto result = std::make_unique<SurfacePlanBinding>();
    if (!isSurfaceKindCompatible(expectedPolicy, expectedKind)) {
      return {{}, SurfacePlanBindingError::IncompatibleExpectedPolicyKind};
    }
    if (!source.isValid()) {
      return {{}, SurfacePlanBindingError::InvalidSource};
    }
    if (globalLayout.surfaces == nullptr || globalLayout.nSurfaces == 0 ||
        !ownedSurfaces.isSubsetOf(SurfaceMask{globalLayout.nSurfaces == MaxLayoutSurfaces ? uint32_t{0xffffffff} : (uint32_t{1} << globalLayout.nSurfaces) - 1}) ||
        ownedSurfaces.count() != static_cast<int>(orderedSurfaces.size())) {
      return {{}, SurfacePlanBindingError::InvalidSurfaceMask};
    }

    result->mSource = source;
    result->mOwnedSurfaces = ownedSurfaces;
    // Retain the validated positional order itself.  The inverse map below
    // is useful for sparse global ids, but it cannot answer the hot-loop
    // question "which surface is position i?" without reconstructing an
    // order from a numeric SurfaceId.  The plan's ordered positions are the
    // runtime traversal authority; the vector is immutable after build().
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

    const TransitionPolicyGrouping grouping{globalLayout};
    if (!grouping.valid()) {
      return {{}, SurfacePlanBindingError::InvalidTopology};
    }
    const auto& topology = globalLayout;
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
      if (globalLayout.getSurface(transition.from).kind != expectedKind ||
          globalLayout.getSurface(transition.to).kind != expectedKind) {
        return {{}, SurfacePlanBindingError::InvalidPolicySurface};
      }
    }

    // Existing grouping supplies the authoritative global-ID order. Filter
    // it only by ownership; never rebuild/rebase a detector-local topology.
    for (const auto id : grouping.transitionsForTag(expectedPolicy)) {
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
    for (const auto id : grouping.cellsForTag(expectedPolicy)) {
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
    for (const auto id : grouping.roadStartCellsForTag(expectedPolicy)) {
      if (result->mScratchCellSlot[id.value()] >= 0) {
        result->mGlobalRoadStartCells.push_back(id);
      }
    }
    // Same source grouping as mGlobalCells above, but in
    // TransitionPolicyGrouping::scheduledCellsForTag()'s rank-sorted
    // (neighbour-schedule) order rather than ascending CellTopologyId order.
    // Every id here is already confirmed owned (mScratchCellSlot valid) by
    // the mGlobalCells loop above, so this is a pure reordering of the same
    // set, never a superset/subset of it.
    for (const auto id : grouping.scheduledCellsForTag(expectedPolicy)) {
      if (result->mScratchCellSlot[id.value()] >= 0) {
        result->mGlobalScheduledCells.push_back(id);
      }
    }
    return {std::move(result), SurfacePlanBindingError::None};
  }

  ClusterSourceId getSource() const noexcept { return mSource; }
  SurfaceMask getOwnedSurfaces() const noexcept { return mOwnedSurfaces; }
  gsl::span<const SurfaceId> getOrderedSurfaces() const noexcept { return mOrderedSurfaces; }
  std::optional<uint16_t> getOwnedSurfaceIndex(SurfaceId id) const noexcept { return getSlot(mOwnedSurfaceIndexBySurface, id); }
  std::optional<uint16_t> getScratchTransitionSlot(TransitionId id) const noexcept { return getSlot(mScratchTransitionSlot, id); }
  std::optional<uint16_t> getScratchCellSlot(CellTopologyId id) const noexcept { return getSlot(mScratchCellSlot, id); }
  gsl::span<const TransitionId> getGlobalTransitions() const noexcept { return mGlobalTransitions; }
  gsl::span<const CellTopologyId> getGlobalCells() const noexcept { return mGlobalCells; }
  gsl::span<const CellTopologyId> getGlobalRoadStartCells() const noexcept { return mGlobalRoadStartCells; }
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

  ClusterSourceId mSource{};
  SurfaceMask mOwnedSurfaces{};
  std::vector<SurfaceId> mOrderedSurfaces;
  std::vector<int16_t> mOwnedSurfaceIndexBySurface;
  std::vector<int16_t> mScratchTransitionSlot;
  std::vector<int16_t> mScratchCellSlot;
  std::vector<TransitionId> mGlobalTransitions;
  std::vector<CellTopologyId> mGlobalCells;
  std::vector<CellTopologyId> mGlobalRoadStartCells;
  std::vector<CellTopologyId> mGlobalScheduledCells;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_SURFACEPLANBINDING_H_
