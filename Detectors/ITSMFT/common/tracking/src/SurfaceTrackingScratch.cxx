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
/// \file SurfaceTrackingScratch.cxx
/// \brief M6c: see SurfaceTrackingScratch.h for the full design rationale.
///

#include "ITSMFTTracking/SurfaceTrackingScratch.h"

#include <algorithm>
#include <limits>

namespace o2::itsmft::tracking
{

using o2::its::bounded_vector;
using o2::its::clearResizeBoundedVector;
using o2::its::deepVectorClear;

void SurfaceTrackingScratch::adoptPlan(std::size_t nOwnedSurfaces, std::size_t nTransitions, std::size_t nCells)
{
  mNOwnedSurfaces = nOwnedSurfaces;
  mNTransitions = nTransitions;
  mNCells = nCells;

  // Group A -- one slot per owned surface.
  clearResizeBoundedVector(mClusters, nOwnedSurfaces, getMaybeFrameworkHostResource());
  clearResizeBoundedVector(mUnsortedClusters, nOwnedSurfaces, getMaybeFrameworkHostResource());
  clearResizeBoundedVector(mTrackingFrameInfo, nOwnedSurfaces, getMaybeFrameworkHostResource());
  clearResizeBoundedVector(mClusterExternalIndices, nOwnedSurfaces, mMemoryPool.get());
  clearResizeBoundedVector(mClusterSize, nOwnedSurfaces, mMemoryPool.get());
  clearResizeBoundedVector(mROFramesClusters, nOwnedSurfaces, getMaybeFrameworkHostResource());
  mClusterLabels.assign(nOwnedSurfaces, nullptr);
  clearResizeBoundedVector(mIndexTables, nOwnedSurfaces, getMaybeFrameworkHostResource());
  clearResizeBoundedVector(mUsedClusters, nOwnedSurfaces, getMaybeFrameworkHostResource());
  clearResizeBoundedVector(mNClustersPerROF, nOwnedSurfaces, mMemoryPool.get());
  mMinR.assign(nOwnedSurfaces, std::numeric_limits<float>::max());
  mMaxR.assign(nOwnedSurfaces, std::numeric_limits<float>::lowest());
  clearResizeBoundedVector(mBogusClusters, nOwnedSurfaces, mMemoryPool.get());
  clearResizeBoundedVector(mPositionResolution, nOwnedSurfaces, mMemoryPool.get());

  // Group B -- sparse transition/cell counts, already runtime today.
  clearResizeBoundedVector(mTracklets, nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletsLookupTable, nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletLabels, nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTransitionPhiCuts, nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTransitionMSAngles, nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mCells, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsLookupTable, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighbours, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursTopology, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursLUT, nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellLabels, nCells, mMemoryPool.get());
}

void SurfaceTrackingScratch::reset()
{
  // Group B.
  deepVectorClear(mTracklets);
  deepVectorClear(mTrackletsLookupTable);
  deepVectorClear(mCells);
  deepVectorClear(mCellsLookupTable);
  deepVectorClear(mCellsNeighbours);
  deepVectorClear(mCellsNeighboursTopology);
  deepVectorClear(mCellsNeighboursLUT);
  deepVectorClear(mTransitionPhiCuts);
  deepVectorClear(mTransitionMSAngles);

  // Group A (allocator-backed, always cleared regardless of framework allocator).
  deepVectorClear(mClusterExternalIndices);
  deepVectorClear(mNClustersPerROF);
  deepVectorClear(mBogusClusters);
  deepVectorClear(mPositionResolution);
  deepVectorClear(mClusterSize);
  deepVectorClear(mNTrackletsPerCluster);
  deepVectorClear(mNTrackletsPerClusterSum);

  // Group E.
  deepVectorClear(mPValphaX);

  // Group D.
  deepVectorClear(mTrackletsIndexROF);
  deepVectorClear(mTrackletClusters);
  deepVectorClear(mLines);
  mTotalTracklets = {0, 0};
  mTotalLines = 0;
  mNTotalLowPtVertices = 0;

  // If we use the external host allocator, the assumption is that we don't
  // clear that memory ourselves -- mirrors
  // LegacyTrackerScratch<NLayers>::resetScratch() exactly.
  if (!hasFrameworkAllocator()) {
    deepVectorClear(mClusters);
    deepVectorClear(mUsedClusters);
    deepVectorClear(mUnsortedClusters);
    deepVectorClear(mIndexTables);
    deepVectorClear(mTrackingFrameInfo);
    deepVectorClear(mROFramesClusters);
  }

  // Only needed to clear if we have MC info -- mirrors resetScratch() exactly.
  if (hasMCinformation()) {
    deepVectorClear(mLinesLabels);
    deepVectorClear(mTrackletLabels);
    deepVectorClear(mCellLabels);
  }

  // mClusterLabels holds non-owning pointers into caller-supplied MC label
  // containers, not owned storage -- reset to nullptr, not deepVectorClear'd,
  // and the vector is not resized (mirrors LegacyTrackerScratch<NLayers>'s
  // own std::array::fill(nullptr), which cannot resize).
  std::fill(mClusterLabels.begin(), mClusterLabels.end(), nullptr);
}

void SurfaceTrackingScratch::setMemoryPool(std::shared_ptr<o2::its::BoundedMemoryResource> pool)
{
  mMemoryPool = std::move(pool);

  auto initVector = [&]<typename T>(bounded_vector<T>& vec, bool useExternal = false) {
    std::pmr::memory_resource* mr = useExternal ? mExtMemoryPool.get() : mMemoryPool.get();
    deepVectorClear(vec, mr);
  };
  auto initContainers = [&]<typename Container>(Container& container, bool useExternal = false) {
    for (auto& v : container) {
      initVector(v, useExternal);
    }
  };

  // Host-only, mirrors LegacyTrackerScratch<NLayers>::setMemoryPool().
  initContainers(mClusterExternalIndices);
  initContainers(mNTrackletsPerCluster);
  initContainers(mNTrackletsPerClusterSum);
  initContainers(mNClustersPerROF);
  initVector(mTransitionPhiCuts);
  initVector(mTransitionMSAngles);
  initVector(mPositionResolution);
  initContainers(mClusterSize);
  initVector(mPValphaX);
  initVector(mBogusClusters);
  initContainers(mTrackletsIndexROF);
  initContainers(mTracklets);
  initContainers(mCells);
  initContainers(mCellsNeighbours);
  initContainers(mCellsLookupTable);
  // MC info (we don't know if we have MC).
  initContainers(mLinesLabels);
  initContainers(mTrackletLabels);
  initContainers(mCellLabels);
  // May use an externally provided allocator.
  initContainers(mClusters, hasFrameworkAllocator());
  initContainers(mUsedClusters, hasFrameworkAllocator());
  initContainers(mUnsortedClusters, hasFrameworkAllocator());
  initContainers(mIndexTables, hasFrameworkAllocator());
  initContainers(mTrackingFrameInfo, hasFrameworkAllocator());
  initContainers(mROFramesClusters, hasFrameworkAllocator());
}

void SurfaceTrackingScratch::setFrameworkAllocator(o2::its::ExternalAllocator* ext)
{
  mExternalAllocator = ext;
  mExtMemoryPool = std::make_shared<o2::its::BoundedMemoryResource>(mExternalAllocator);
}

namespace
{
template <typename T>
bool flatAllocatorMatches(const bounded_vector<T>& a, const bounded_vector<T>& b) noexcept
{
  return a.get_allocator().resource() == b.get_allocator().resource();
}

template <typename T, std::size_t N>
bool flatArrayAllocatorMatches(const std::array<bounded_vector<T>, N>& a, const std::array<bounded_vector<T>, N>& b) noexcept
{
  for (std::size_t i = 0; i < N; ++i) {
    if (!flatAllocatorMatches(a[i], b[i])) {
      return false;
    }
  }
  return true;
}
} // namespace

bool SurfaceTrackingScratch::allocatorsMatch(const SurfaceTrackingScratch& staged) const noexcept
{
  // Every flat bounded_vector<T> member swap() exchanges directly (the ones
  // whose allocators must compare equal for a well-defined
  // bounded_vector::swap()). Per-owned-surface/per-transition/per-cell
  // *outer* std::vector<bounded_vector<T>> members carry no such precondition
  // (see the header doc) and are intentionally not checked here.
  return flatAllocatorMatches(mBogusClusters, staged.mBogusClusters) &&
         flatAllocatorMatches(mPositionResolution, staged.mPositionResolution) &&
         flatAllocatorMatches(mTransitionPhiCuts, staged.mTransitionPhiCuts) &&
         flatAllocatorMatches(mTransitionMSAngles, staged.mTransitionMSAngles) &&
         flatAllocatorMatches(mPValphaX, staged.mPValphaX) &&
         flatArrayAllocatorMatches(mNTrackletsPerCluster, staged.mNTrackletsPerCluster) &&
         flatArrayAllocatorMatches(mNTrackletsPerClusterSum, staged.mNTrackletsPerClusterSum) &&
         flatArrayAllocatorMatches(mTrackletsIndexROF, staged.mTrackletsIndexROF);
}

void SurfaceTrackingScratch::swap(SurfaceTrackingScratch& other) noexcept
{
  static_assert(noexcept(std::declval<bounded_vector<float>&>().swap(std::declval<bounded_vector<float>&>())));
  static_assert(noexcept(std::declval<bounded_vector<int>&>().swap(std::declval<bounded_vector<int>&>())));
  static_assert(noexcept(std::declval<bounded_vector<std::array<float, 2>>&>().swap(std::declval<bounded_vector<std::array<float, 2>>&>())));

  std::swap(mNOwnedSurfaces, other.mNOwnedSurfaces);
  std::swap(mNTransitions, other.mNTransitions);
  std::swap(mNCells, other.mNCells);

  // Outer std::vector<bounded_vector<T>> containers: always safe, see the
  // header doc.
  mClusters.swap(other.mClusters);
  mUnsortedClusters.swap(other.mUnsortedClusters);
  mTrackingFrameInfo.swap(other.mTrackingFrameInfo);
  mClusterExternalIndices.swap(other.mClusterExternalIndices);
  mClusterSize.swap(other.mClusterSize);
  mROFramesClusters.swap(other.mROFramesClusters);
  mClusterLabels.swap(other.mClusterLabels);
  mIndexTables.swap(other.mIndexTables);
  mUsedClusters.swap(other.mUsedClusters);
  mNClustersPerROF.swap(other.mNClustersPerROF);
  mMinR.swap(other.mMinR);
  mMaxR.swap(other.mMaxR);
  mTracklets.swap(other.mTracklets);
  mTrackletsLookupTable.swap(other.mTrackletsLookupTable);
  mTrackletLabels.swap(other.mTrackletLabels);
  mCells.swap(other.mCells);
  mCellsLookupTable.swap(other.mCellsLookupTable);
  mCellsNeighbours.swap(other.mCellsNeighbours);
  mCellsNeighboursTopology.swap(other.mCellsNeighboursTopology);
  mCellsNeighboursLUT.swap(other.mCellsNeighboursLUT);
  mCellLabels.swap(other.mCellLabels);
  mLines.swap(other.mLines);
  mTrackletClusters.swap(other.mTrackletClusters);
  mNTrackletsPerROF.swap(other.mNTrackletsPerROF);
  mLinesLabels.swap(other.mLinesLabels);

  // Flat bounded_vector<T> containers: precondition allocatorsMatch(other).
  mBogusClusters.swap(other.mBogusClusters);
  mPositionResolution.swap(other.mPositionResolution);
  mTransitionPhiCuts.swap(other.mTransitionPhiCuts);
  mTransitionMSAngles.swap(other.mTransitionMSAngles);
  mPValphaX.swap(other.mPValphaX);
  for (std::size_t i = 0; i < mNTrackletsPerCluster.size(); ++i) {
    mNTrackletsPerCluster[i].swap(other.mNTrackletsPerCluster[i]);
    mNTrackletsPerClusterSum[i].swap(other.mNTrackletsPerClusterSum[i]);
    mTrackletsIndexROF[i].swap(other.mTrackletsIndexROF[i]);
  }

  std::swap(mTotalTracklets, other.mTotalTracklets);
  std::swap(mTotalLines, other.mTotalLines);
  std::swap(mNTotalLowPtVertices, other.mNTotalLowPtVertices);
  std::swap(mIsStaggered, other.mIsStaggered);

  // mMemoryPool/mExtMemoryPool/mExternalAllocator deliberately not swapped --
  // see the header doc.
}

} // namespace o2::itsmft::tracking
