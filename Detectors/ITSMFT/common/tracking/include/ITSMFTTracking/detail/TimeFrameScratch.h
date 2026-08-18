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
///
/// \file TimeFrameScratch.h
/// \brief Runtime-plan-owned, detector-neutral CA workspace.
///
/// Host storage follows the runtime surface graph; device capacities remain
/// fixed. TimeFrame owns the workspace, while adapters own raw ROFs.
#ifndef ALICEO2_ITSMFT_TRACKING_TimeFrameScratch_H_
#define ALICEO2_ITSMFT_TRACKING_TimeFrameScratch_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <vector>

#include <gsl/gsl>

#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/GenericTrack.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TraversalTopology.h"
#include "ITSMFTTracking/TrackingPrimitives.h"
#include "ITSMFTTracking/detail/TrackingKernelParameters.h"
#include "ITStracking/BoundedAllocator.h"
#include <optional>
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::itsmft::tracking
{

// Per-iteration derived traversal state. The TimeFrame workspace owns this
// storage; Tracker builds a short-lived view for a traversal call.
struct TraversalWorkspace {
  TrackingKernelParameters kernelParameters{};
  AttachHitConfigView attachHitConfig{};
  std::vector<NominalSurfaceMaterial> layerMaterial;
  std::vector<gsl::span<const GlobalMeasurement>> layerGlobalMeasurements;
  std::vector<float> diskLayerReferenceZ;
  gsl::span<const float> diskLayerReferenceZView{};
  o2::its::bounded_vector<TrackingCandidate> acceptedTracks;
  // Per-pass traversal plan. The graph remains immutable configuration; the
  // tracker derives this selected topology and its compact scratch mapping.
  SurfaceMask activeSurfaces{};
  std::vector<LayerId> orderedSurfaces;
  std::vector<int16_t> surfaceSlotById;
  std::vector<int16_t> edgeSlotById;
  std::vector<int16_t> cellSlotById;
  std::vector<EdgeId> edges;
  std::vector<CellPathId> cells;
  std::vector<CellPathId> roadStartCells;
  std::vector<uint32_t> roadStartComponentOffsets;
  std::vector<CellPathId> scheduledCells;
  SurfaceCatalogView topologyCatalog{};
  TraversalTopology topology;
  bool valid{false};

  std::optional<uint16_t> getSurfaceSlot(LayerId id) const noexcept;
  std::optional<uint16_t> getEdgeSlot(EdgeId id) const noexcept;
  std::optional<uint16_t> getCellSlot(CellPathId id) const noexcept;
  TraversalTopologyView getTopologyView() const noexcept { return topology.getView(topologyCatalog); }

  void reset(std::pmr::memory_resource* resource) noexcept
  {
    kernelParameters = {};
    attachHitConfig = {};
    layerMaterial.clear();
    layerGlobalMeasurements.clear();
    diskLayerReferenceZ.clear();
    diskLayerReferenceZView = {};
    o2::its::deepVectorClear(acceptedTracks, resource);
    activeSurfaces = {};
    orderedSurfaces = {};
    surfaceSlotById.clear();
    edgeSlotById.clear();
    cellSlotById.clear();
    edges.clear();
    cells.clear();
    roadStartCells.clear();
    roadStartComponentOffsets.clear();
    scheduledCells.clear();
    topologyCatalog = {};
    topology = {};
    valid = false;
  }
};

/// Detector-neutral CA state rebuilt for each tracking iteration.
class TimeFrameScratch
{
 private:
  // Pool must outlive allocator-backed members.
  std::shared_ptr<o2::its::BoundedMemoryResource> mMemoryPool;

 public:
  TimeFrameScratch() = default;
  ~TimeFrameScratch() = default;
  TimeFrameScratch(const TimeFrameScratch&) = delete;
  TimeFrameScratch& operator=(const TimeFrameScratch&) = delete;
  TimeFrameScratch(TimeFrameScratch&&) = delete;
  TimeFrameScratch& operator=(TimeFrameScratch&&) = delete;

  /// Size surface, edge and cell storage; setMemoryPool() comes first.
  void adoptPlan(std::size_t nOwnedSurfaces, std::size_t nEdges, std::size_t nCells);
  void configureTraversalWorkspaces(std::size_t nIterations);
  TraversalWorkspace& getTraversalWorkspace(std::size_t iteration) { return mTraversalWorkspaces.at(iteration); }
  const TraversalWorkspace& getTraversalWorkspace(std::size_t iteration) const { return mTraversalWorkspaces.at(iteration); }
  std::size_t getNTraversalWorkspaces() const noexcept { return mTraversalWorkspaces.size(); }

  std::size_t getNOwnedSurfaces() const noexcept { return mNOwnedSurfaces; }
  std::size_t getNEdges() const noexcept { return mNEdges; }
  std::size_t getNCells() const noexcept { return mNCells; }

  /// Clear iteration state without changing plan sizes.
  void reset();

  /// Reseat allocator-backed containers.
  void setMemoryPool(std::shared_ptr<o2::its::BoundedMemoryResource> pool);
  auto& getMemoryPool() const noexcept { return mMemoryPool; }
  float getEdgePhiCut(int edgeId) const { return mEdgePhiCuts[edgeId]; }
  float getEdgeMSAngle(int edgeId) const { return mEdgeMSAngles[edgeId]; }
  auto& getEdgePhiCuts() { return mEdgePhiCuts; }
  auto& getEdgeMSAngles() { return mEdgeMSAngles; }
  float getPositionResolution(int layer) const { return mPositionResolution[layer]; }
  auto& getPositionResolutions() { return mPositionResolution; }

  auto& getTrackletsLabel(int layer) { return mTrackletLabels[layer]; }
  auto& getCellsLabel(int layer) { return mCellLabels[layer]; }

  void initialise(TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers, int iteration,
                  const IndexTableUtilsCore& indexTableConfig, TraversalTopologyView topology,
                  gsl::span<const EdgeId> edgeIds, gsl::span<const CellPathId> cellIds,
                  gsl::span<const LayerId> orderedSurfaces,
                  gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements);
  void initialise(TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers, int iteration,
                  gsl::span<const IndexTableUtilsCore> indexTableConfigs, TraversalTopologyView topology,
                  gsl::span<const EdgeId> edgeIds, gsl::span<const CellPathId> cellIds,
                  gsl::span<const LayerId> orderedSurfaces,
                  gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements);

  auto& getTracklets() { return mTracklets; }
  auto& getTrackletsLookupTable() { return mTrackletsLookupTable; }

  auto& getCells() { return mCells; }
  const auto& getCells() const { return mCells; }

  auto& getCellsLookupTable() { return mCellsLookupTable; }
  auto& getCellsNeighbours() { return mCellsNeighbours; }
  auto& getCellsNeighboursTopology() { return mCellsNeighboursTopology; }
  auto& getCellsNeighboursLUT() { return mCellsNeighboursLUT; }
  size_t getNumberOfCells() const;
  size_t getNumberOfTracklets() const;
  size_t getNumberOfNeighbours() const;

  // ---- Per-iteration surface and CA construction state ----
  o2::its::bounded_vector<float> mPositionResolution;
  std::vector<o2::its::bounded_vector<Tracklet>> mTracklets;
  std::vector<o2::its::bounded_vector<int>> mTrackletsLookupTable;
  std::vector<o2::its::bounded_vector<o2::MCCompLabel>> mTrackletLabels;
  o2::its::bounded_vector<float> mEdgePhiCuts;
  o2::its::bounded_vector<float> mEdgeMSAngles;
  std::vector<o2::its::bounded_vector<CellSeed>> mCells;
  std::vector<o2::its::bounded_vector<int>> mCellsLookupTable;
  std::vector<o2::its::bounded_vector<int>> mCellsNeighbours;
  std::vector<o2::its::bounded_vector<int>> mCellsNeighboursTopology;
  std::vector<o2::its::bounded_vector<int>> mCellsNeighboursLUT;
  std::vector<o2::its::bounded_vector<o2::MCCompLabel>> mCellLabels;
  std::vector<TraversalWorkspace> mTraversalWorkspaces;

 private:
  std::size_t mNOwnedSurfaces{0};
  std::size_t mNEdges{0};
  std::size_t mNCells{0};
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TimeFrameScratch_H_ */
