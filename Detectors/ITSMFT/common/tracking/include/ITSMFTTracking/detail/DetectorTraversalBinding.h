// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORTRAVERSALBINDING_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORTRAVERSALBINDING_H_

// Host-only bridge for the first combined ITS+MFT execution slice.  It does
// not own or project a topology: it retains global topology identifiers and
// provides only the compact legacy-scratch storage slots for one detector.

#ifndef GPUCA_GPUCODE

#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

#include <gsl/span>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/detail/TransitionPolicyDispatch.h"

namespace o2::itsmft::tracking
{

enum class DetectorTraversalBindingError : uint8_t {
  None,
  UnsupportedDetector,
  InvalidSource,
  InvalidSurfaceMask,
  InvalidLegacySurfaceOrder,
  SurfaceDetectorMismatch,
  InvalidTopology,
  InvalidPolicySurface,
  CrossBoundaryTransition,
  CrossBoundaryCell,
  IncompleteTransitionMapping,
  IncompleteCellMapping,
  DuplicateScratchSlot
};

class DetectorTraversalBinding
{
 public:
  struct BuildResult {
    std::unique_ptr<DetectorTraversalBinding> binding{};
    DetectorTraversalBindingError error{DetectorTraversalBindingError::None};
    bool ok() const noexcept { return static_cast<bool>(binding); }
  };

  // `legacySurfaceOrder` maps legacy layer i to a global SurfaceId. It is the
  // only global-to-legacy-layer projection in the combined path; no local
  // detector catalog or topology is accepted here.
  static BuildResult build(const DetectorLayoutView& globalLayout,
                           o2::detectors::DetID::ID detector,
                           ClusterSourceId source,
                           SurfaceMask ownedSurfaces,
                           gsl::span<const SurfaceId> legacySurfaceOrder)
  {
    auto result = std::make_unique<DetectorTraversalBinding>();
    if ((detector != o2::detectors::DetID::ITS && detector != o2::detectors::DetID::MFT)) {
      return {{}, DetectorTraversalBindingError::UnsupportedDetector};
    }
    if (!source.isValid()) {
      return {{}, DetectorTraversalBindingError::InvalidSource};
    }
    if (globalLayout.surfaces == nullptr || globalLayout.nSurfaces == 0 ||
        !ownedSurfaces.isSubsetOf(SurfaceMask{globalLayout.nSurfaces == MaxLayoutSurfaces ? uint32_t{0xffffffff} : (uint32_t{1} << globalLayout.nSurfaces) - 1}) ||
        ownedSurfaces.count() != static_cast<int>(legacySurfaceOrder.size())) {
      return {{}, DetectorTraversalBindingError::InvalidSurfaceMask};
    }

    result->mDetector = detector;
    result->mSource = source;
    result->mOwnedSurfaces = ownedSurfaces;
    result->mLegacyLayerBySurface.assign(globalLayout.nSurfaces, -1);
    for (uint16_t layer = 0; layer < legacySurfaceOrder.size(); ++layer) {
      const auto surface = legacySurfaceOrder[layer];
      if (!surface.isValid() || surface.value() >= globalLayout.nSurfaces || !ownedSurfaces.has(surface) ||
          result->mLegacyLayerBySurface[surface.value()] >= 0) {
        return {{}, DetectorTraversalBindingError::InvalidLegacySurfaceOrder};
      }
      if (globalLayout.getSurface(surface).detectorId != static_cast<uint8_t>(detector)) {
        return {{}, DetectorTraversalBindingError::SurfaceDetectorMismatch};
      }
      result->mLegacyLayerBySurface[surface.value()] = static_cast<int16_t>(layer);
    }
    for (uint16_t surface = 0; surface < globalLayout.nSurfaces; ++surface) {
      if (ownedSurfaces.has(SurfaceId{surface}) && result->mLegacyLayerBySurface[surface] < 0) {
        return {{}, DetectorTraversalBindingError::InvalidLegacySurfaceOrder};
      }
    }

    const auto expectedKind = detector == o2::detectors::DetID::ITS ? SurfaceKind::Cylinder : SurfaceKind::Disk;
    const auto expectedPolicy = detector == o2::detectors::DetID::ITS ? TransitionPolicyTag::CylinderCylinder : TransitionPolicyTag::DiskDisk;
    const TransitionPolicyGrouping grouping{globalLayout};
    if (!grouping.valid()) {
      return {{}, DetectorTraversalBindingError::InvalidTopology};
    }
    const auto& topology = globalLayout.topology;
    result->mScratchTransitionSlot.assign(topology.nTransitions, -1);
    result->mScratchCellSlot.assign(topology.nCells, -1);

    // Validate every global transition first. A boundary edge is invalid for
    // either binding, even if it would otherwise be filtered from one side.
    for (uint32_t id = 0; id < topology.nTransitions; ++id) {
      const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(id)});
      const bool fromOwned = ownedSurfaces.has(transition.from);
      const bool toOwned = ownedSurfaces.has(transition.to);
      if (fromOwned != toOwned) {
        return {{}, DetectorTraversalBindingError::CrossBoundaryTransition};
      }
      if (!fromOwned) {
        continue;
      }
      if (!transition.skippedSurfaces.isSubsetOf(ownedSurfaces)) {
        return {{}, DetectorTraversalBindingError::CrossBoundaryTransition};
      }
      if (globalLayout.getSurface(transition.from).kind != expectedKind ||
          globalLayout.getSurface(transition.to).kind != expectedKind) {
        return {{}, DetectorTraversalBindingError::InvalidPolicySurface};
      }
    }

    // Existing grouping supplies the authoritative global-ID order. Filter it
    // only by ownership; never rebuild/rebase a detector-local topology.
    for (const auto id : grouping.transitionsForTag(expectedPolicy)) {
      const auto& transition = topology.getTransition(id);
      if (!ownedSurfaces.has(transition.from)) {
        continue;
      }
      if (result->mScratchTransitionSlot[id.value()] >= 0) {
        return {{}, DetectorTraversalBindingError::DuplicateScratchSlot};
      }
      result->mScratchTransitionSlot[id.value()] = static_cast<int16_t>(result->mGlobalTransitions.size());
      result->mGlobalTransitions.push_back(id);
    }
    for (uint32_t id = 0; id < topology.nTransitions; ++id) {
      const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(id)});
      if (ownedSurfaces.has(transition.from) && result->mScratchTransitionSlot[id] < 0) {
        return {{}, DetectorTraversalBindingError::IncompleteTransitionMapping};
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
        return {{}, DetectorTraversalBindingError::CrossBoundaryCell};
      }
    }
    for (const auto id : grouping.cellsForTag(expectedPolicy)) {
      const auto& cell = topology.getCell(id);
      if (result->mScratchTransitionSlot[cell.firstTransition.value()] < 0) {
        continue;
      }
      if (result->mScratchCellSlot[id.value()] >= 0) {
        return {{}, DetectorTraversalBindingError::DuplicateScratchSlot};
      }
      result->mScratchCellSlot[id.value()] = static_cast<int16_t>(result->mGlobalCells.size());
      result->mGlobalCells.push_back(id);
    }
    for (uint32_t id = 0; id < topology.nCells; ++id) {
      const auto& cell = topology.getCell(CellTopologyId{static_cast<uint16_t>(id)});
      if (result->mScratchTransitionSlot[cell.firstTransition.value()] >= 0 && result->mScratchCellSlot[id] < 0) {
        return {{}, DetectorTraversalBindingError::IncompleteCellMapping};
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
    // findCellsNeighboursForPolicy() (TrackerTraits.cxx) is the only
    // consumer; every id here is already confirmed owned (mScratchCellSlot
    // valid) by the mGlobalCells loop above, so this is a pure reordering of
    // the same set, never a superset/subset of it.
    for (const auto id : grouping.scheduledCellsForTag(expectedPolicy)) {
      if (result->mScratchCellSlot[id.value()] >= 0) {
        result->mGlobalScheduledCells.push_back(id);
      }
    }
    return {std::move(result), DetectorTraversalBindingError::None};
  }

  o2::detectors::DetID::ID getDetector() const noexcept { return mDetector; }
  ClusterSourceId getSource() const noexcept { return mSource; }
  SurfaceMask getOwnedSurfaces() const noexcept { return mOwnedSurfaces; }
  std::optional<uint16_t> getLegacyLayer(SurfaceId id) const noexcept { return getSlot(mLegacyLayerBySurface, id); }
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

  o2::detectors::DetID::ID mDetector{o2::detectors::DetID::ITS};
  ClusterSourceId mSource{};
  SurfaceMask mOwnedSurfaces{};
  std::vector<int16_t> mLegacyLayerBySurface;
  std::vector<int16_t> mScratchTransitionSlot;
  std::vector<int16_t> mScratchCellSlot;
  std::vector<TransitionId> mGlobalTransitions;
  std::vector<CellTopologyId> mGlobalCells;
  std::vector<CellTopologyId> mGlobalRoadStartCells;
  std::vector<CellTopologyId> mGlobalScheduledCells;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif // ALICEO2_ITSMFT_TRACKING_DETECTORTRAVERSALBINDING_H_
