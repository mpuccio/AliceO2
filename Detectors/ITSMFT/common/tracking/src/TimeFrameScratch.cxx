// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/detail/TimeFrameScratch.h"

#include <algorithm>
#include <limits>
#include <new>
#include <stdexcept>

#include "ITSMFTTracking/IndexTableConfiguration.h"

namespace o2::itsmft::tracking
{

using o2::its::clearResizeBoundedVector;
using o2::its::deepVectorClear;

namespace
{
template <typename Id>
std::optional<uint16_t> traversalSlot(const std::vector<int16_t>& slots, Id id) noexcept
{
  if (!id.isValid() || id.value() >= slots.size() || slots[id.value()] < 0) {
    return std::nullopt;
  }
  return static_cast<uint16_t>(slots[id.value()]);
}
} // namespace

std::optional<uint16_t> TraversalWorkspace::getSurfaceSlot(LayerId id) const noexcept { return traversalSlot(surfaceSlotById, id); }
std::optional<uint16_t> TraversalWorkspace::getEdgeSlot(EdgeId id) const noexcept { return traversalSlot(edgeSlotById, id); }
std::optional<uint16_t> TraversalWorkspace::getCellSlot(CellPathId id) const noexcept { return traversalSlot(cellSlotById, id); }

void TimeFrameScratch::adoptPlan(std::size_t nOwnedSurfaces, std::size_t nEdges, std::size_t nCells)
{
  mNOwnedSurfaces = nOwnedSurfaces;
  mNEdges = nEdges;
  mNCells = nCells;
  clearResizeBoundedVector(mPositionResolution, nOwnedSurfaces, mMemoryPool.get());
  clearResizeBoundedVector(mTracklets, nEdges, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletsLookupTable, nEdges, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletLabels, nEdges, mMemoryPool.get());
  clearResizeBoundedVector(mEdgePhiCuts, nEdges, mMemoryPool.get());
  clearResizeBoundedVector(mEdgeMSAngles, nEdges, mMemoryPool.get());
  clearResizeBoundedVector(mCells, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsLookupTable, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighbours, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursTopology, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursLUT, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellLabels, nCells, mMemoryPool.get());
}

void TimeFrameScratch::configureTraversalWorkspaces(std::size_t nIterations)
{
  mTraversalWorkspaces.resize(nIterations);
  for (auto& workspace : mTraversalWorkspaces) {
    workspace.reset(mMemoryPool.get());
  }
}

void TimeFrameScratch::reset()
{
  deepVectorClear(mTracklets);
  deepVectorClear(mTrackletsLookupTable);
  deepVectorClear(mTrackletLabels);
  deepVectorClear(mCells);
  deepVectorClear(mCellsLookupTable);
  deepVectorClear(mCellsNeighbours);
  deepVectorClear(mCellsNeighboursTopology);
  deepVectorClear(mCellsNeighboursLUT);
  deepVectorClear(mCellLabels);
  deepVectorClear(mEdgePhiCuts);
  deepVectorClear(mEdgeMSAngles);
  deepVectorClear(mPositionResolution);
  for (auto& workspace : mTraversalWorkspaces) {
    workspace.reset(mMemoryPool.get());
  }
}

void TimeFrameScratch::setMemoryPool(std::shared_ptr<o2::its::BoundedMemoryResource> pool)
{
  mMemoryPool = std::move(pool);
  auto initVector = [&]<typename T>(bounded_vector<T>& vector) { deepVectorClear(vector, mMemoryPool.get()); };
  auto initContainers = [&]<typename Container>(Container& container) {
    for (auto& vector : container) {
      initVector(vector);
    }
  };
  initVector(mPositionResolution);
  initVector(mEdgePhiCuts);
  initVector(mEdgeMSAngles);
  initContainers(mTracklets);
  initContainers(mTrackletsLookupTable);
  initContainers(mTrackletLabels);
  initContainers(mCells);
  initContainers(mCellsLookupTable);
  initContainers(mCellsNeighbours);
  initContainers(mCellsNeighboursTopology);
  initContainers(mCellsNeighboursLUT);
  initContainers(mCellLabels);
  for (auto& workspace : mTraversalWorkspaces) {
    workspace.reset(mMemoryPool.get());
  }
}

std::size_t TimeFrameScratch::getNumberOfCells() const
{
  std::size_t result = 0;
  for (const auto& cells : mCells) {
    result += cells.size();
  }
  return result;
}

std::size_t TimeFrameScratch::getNumberOfTracklets() const
{
  std::size_t result = 0;
  for (const auto& tracklets : mTracklets) {
    result += tracklets.size();
  }
  return result;
}

std::size_t TimeFrameScratch::getNumberOfNeighbours() const
{
  std::size_t result = 0;
  for (const auto& neighbours : mCellsNeighbours) {
    result += neighbours.size();
  }
  return result;
}

void TimeFrameScratch::initialise(TimeFrame& frame, const TrackingParameters& parameters, int maxLayers, int iteration,
                                  const IndexTableUtilsCore& indexTableConfig, TraversalTopologyView topology,
                                  gsl::span<const EdgeId> edgeIds, gsl::span<const CellPathId> cellIds,
                                  gsl::span<const LayerId> orderedSurfaces,
                                  gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements)
{
  std::vector<IndexTableUtilsCore> configs(mNOwnedSurfaces, indexTableConfig);
  initialise(frame, parameters, maxLayers, iteration, configs, topology, edgeIds, cellIds,
             orderedSurfaces, layerMeasurements);
}

void TimeFrameScratch::initialise(TimeFrame& frame, const TrackingParameters& parameters, int maxLayers, int iteration,
                                  gsl::span<const IndexTableUtilsCore> indexTableConfigs, TraversalTopologyView topology,
                                  gsl::span<const EdgeId> edgeIds, gsl::span<const CellPathId> cellIds,
                                  gsl::span<const LayerId> orderedSurfaces,
                                  gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements)
{
  (void)iteration;
  if (orderedSurfaces.size() != mNOwnedSurfaces || indexTableConfigs.size() != mNOwnedSurfaces ||
      edgeIds.size() != mNEdges || cellIds.size() != mNCells ||
      edgeIds.size() > topology.nEdges || cellIds.size() > topology.nPaths ||
      (topology.nEdges != 0 && (topology.edges == nullptr || topology.pathsByFirstEdgeOffsets == nullptr)) ||
      (topology.nPaths != 0 && (topology.paths == nullptr || topology.pathsByFirstEdge == nullptr))) {
    throw std::logic_error{"TimeFrameScratch::initialise(): plan/sparse-topology extent mismatch"};
  }
  const auto surfaceSlot = [&](LayerId surface) {
    const auto it = std::find(orderedSurfaces.begin(), orderedSurfaces.end(), surface);
    if (it == orderedSurfaces.end()) {
      throw std::logic_error{"TimeFrameScratch::initialise(): sparse topology surface is not bound"};
    }
    return static_cast<std::size_t>(std::distance(orderedSurfaces.begin(), it));
  };
  for (const auto edgeId : edgeIds) {
    if (!edgeId.isValid() || edgeId.value() >= topology.nEdges) {
      throw std::logic_error{"TimeFrameScratch::initialise(): invalid sparse edge binding"};
    }
    const auto& edge = topology.getEdge(edgeId);
    (void)surfaceSlot(edge.from);
    (void)surfaceSlot(edge.to);
  }
  for (const auto cellId : cellIds) {
    if (!cellId.isValid() || cellId.value() >= topology.nPaths) {
      throw std::logic_error{"TimeFrameScratch::initialise(): invalid sparse cell binding"};
    }
    const auto& path = topology.getPath(cellId);
    if (!path.first.isValid() || !path.second.isValid() ||
        path.first.value() >= topology.nEdges || path.second.value() >= topology.nEdges) {
      throw std::logic_error{"TimeFrameScratch::initialise(): sparse cell references an invalid edge"};
    }
  }

  if (parameters.PassFlags[IterationStep::FirstPass]) {
    frame.mIndexTableUtils.assign(indexTableConfigs.begin(), indexTableConfigs.end());
    clearResizeBoundedVector(mPositionResolution, maxLayers, mMemoryPool.get());
    clearResizeBoundedVector(frame.mBogusClusters, maxLayers, frame.mMemoryPool.get());
    for (int layer = 0; layer < std::min(maxLayers, static_cast<int>(frame.mClusters.size())); ++layer) {
      clearResizeBoundedVector(frame.mClusters[layer], layerMeasurements[layer].size(), frame.mMemoryPool.get());
      clearResizeBoundedVector(frame.mUsedClusters[layer], layerMeasurements[layer].size(), frame.mMemoryPool.get());
    }
    for (std::size_t layer = 0; layer < mNOwnedSurfaces; ++layer) {
      std::size_t stride = 0;
      if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(frame.mIndexTableUtils[layer].getNrowBins()),
                                        static_cast<std::size_t>(frame.mIndexTableUtils[layer].getNcolBins()), stride) ||
          stride == std::numeric_limits<std::size_t>::max()) {
        throw std::bad_alloc{};
      }
      ++stride;
      std::size_t tableSize = 0;
      if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(frame.getNrof(static_cast<int>(layer))), stride, tableSize)) {
        throw std::bad_alloc{};
      }
      clearResizeBoundedVector(frame.mIndexTables[layer], tableSize, frame.mMemoryPool.get());
    }
    std::fill(frame.mMinR.begin(), frame.mMinR.end(), std::numeric_limits<float>::max());
    std::fill(frame.mMaxR.begin(), frame.mMaxR.end(), std::numeric_limits<float>::lowest());
    std::fill(frame.mMinZ.begin(), frame.mMinZ.end(), std::numeric_limits<float>::max());
    std::fill(frame.mMaxZ.begin(), frame.mMaxZ.end(), std::numeric_limits<float>::lowest());
  }

  clearResizeBoundedVector(mCells, mNCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsLookupTable, mNCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighbours, mNCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursTopology, mNCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursLUT, mNCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellLabels, mNCells, mMemoryPool.get());
  clearResizeBoundedVector(mTracklets, mNEdges, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletLabels, mNEdges, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletsLookupTable, mNEdges, mMemoryPool.get());
  clearResizeBoundedVector(mEdgePhiCuts, mNEdges, mMemoryPool.get());
  clearResizeBoundedVector(mEdgeMSAngles, mNEdges, mMemoryPool.get());

  if (parameters.PassFlags[IterationStep::RebuildClusterLUT]) {
    frame.prepareClusters(maxLayers, layerMeasurements);
  }
  for (std::size_t layer = 0; layer < mNOwnedSurfaces; ++layer) {
    mPositionResolution[layer] = o2::gpu::CAMath::Sqrt(
      0.5f * (parameters.SystError2Col[layer] + parameters.SystError2Row[layer]) +
      parameters.LayerResolution[layer] * parameters.LayerResolution[layer]);
  }
  for (int edge = 0; edge < static_cast<int>(mTracklets.size()); ++edge) {
    const auto fromSlot = surfaceSlot(topology.getEdge(edgeIds[edge]).from);
    deepVectorClear(mTracklets[edge]);
    deepVectorClear(mTrackletLabels[edge]);
    deepVectorClear(mTrackletsLookupTable[edge]);
    mTrackletsLookupTable[edge].resize(frame.mClusters[fromSlot].size() + 1, 0);
  }
  for (int cell = 0; cell < static_cast<int>(mCells.size()); ++cell) {
    deepVectorClear(mCells[cell]);
    deepVectorClear(mCellsLookupTable[cell]);
    deepVectorClear(mCellsNeighbours[cell]);
    deepVectorClear(mCellsNeighboursTopology[cell]);
    deepVectorClear(mCellsNeighboursLUT[cell]);
    deepVectorClear(mCellLabels[cell]);
  }
}

} // namespace o2::itsmft::tracking
