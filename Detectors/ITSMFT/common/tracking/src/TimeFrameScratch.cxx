// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/detail/TimeFrameScratch.h"

#include <stdexcept>

namespace o2::itsmft::tracking
{

using o2::its::clearResizeBoundedVector;
using o2::its::deepVectorClear;

void TimeFrameScratch::configureStorage(std::size_t nEdges, std::size_t nCells)
{
  mNEdges = nEdges;
  mNCells = nCells;
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
  for (auto& iteration : mIterations) {
    iteration.reset(mMemoryPool.get());
  }
}

void TimeFrameScratch::clearStorage() noexcept
{
  mTracklets.clear();
  mTrackletsLookupTable.clear();
  mTrackletLabels.clear();
  mCells.clear();
  mCellsLookupTable.clear();
  mCellsNeighbours.clear();
  mCellsNeighboursTopology.clear();
  mCellsNeighboursLUT.clear();
  mCellLabels.clear();
  deepVectorClear(mEdgePhiCuts);
  deepVectorClear(mEdgeMSAngles);
  mIterations.clear();
  mNEdges = 0;
  mNCells = 0;
}

void TimeFrameScratch::setMemoryPool(std::shared_ptr<o2::its::BoundedMemoryResource> pool)
{
  mMemoryPool = std::move(pool);
  auto initVector = [&]<typename T>(o2::its::bounded_vector<T>& vector) { deepVectorClear(vector, mMemoryPool.get()); };
  auto initContainers = [&]<typename Container>(Container& container) {
    for (auto& vector : container) {
      initVector(vector);
    }
  };
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
  for (auto& iteration : mIterations) {
    iteration.reset(mMemoryPool.get());
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

void TimeFrameScratch::beginIteration(std::size_t nEdges, std::size_t nCells,
                                      gsl::span<const std::size_t> trackletLookupSizes)
{
  if (nEdges > mNEdges || nCells > mNCells || trackletLookupSizes.size() != nEdges) {
    throw std::logic_error{"TimeFrameScratch::beginIteration(): requested storage exceeds configured capacity"};
  }

  clearResizeBoundedVector(mCells, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsLookupTable, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighbours, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursTopology, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursLUT, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellLabels, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mTracklets, nEdges, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletLabels, nEdges, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletsLookupTable, nEdges, mMemoryPool.get());
  clearResizeBoundedVector(mEdgePhiCuts, nEdges, mMemoryPool.get());
  clearResizeBoundedVector(mEdgeMSAngles, nEdges, mMemoryPool.get());

  for (std::size_t edge = 0; edge < nEdges; ++edge) {
    deepVectorClear(mTracklets[edge]);
    deepVectorClear(mTrackletLabels[edge]);
    deepVectorClear(mTrackletsLookupTable[edge]);
    mTrackletsLookupTable[edge].resize(trackletLookupSizes[edge] + 1, 0);
  }
  for (std::size_t cell = 0; cell < nCells; ++cell) {
    deepVectorClear(mCells[cell]);
    deepVectorClear(mCellsLookupTable[cell]);
    deepVectorClear(mCellsNeighbours[cell]);
    deepVectorClear(mCellsNeighboursTopology[cell]);
    deepVectorClear(mCellsNeighboursLUT[cell]);
    deepVectorClear(mCellLabels[cell]);
  }
}

} // namespace o2::itsmft::tracking
