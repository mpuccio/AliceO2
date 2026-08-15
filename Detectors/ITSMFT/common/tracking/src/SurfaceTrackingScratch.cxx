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
/// \brief SurfaceTrackingScratch implementation.
///

#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <new>
#include <stdexcept>
#include <string>

#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/SurfaceMask.h"
#include "ITStracking/MathUtils.h"

namespace o2::itsmft::tracking
{

// TimeFrame.h already defines bounded_vector<T> in this namespace.
using o2::its::clearResizeBoundedVector;
using o2::its::deepVectorClear;

namespace constants
{
using namespace o2::its::constants;
} // namespace constants

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

std::optional<uint16_t> TraversalWorkspace::getSurfaceSlot(SurfaceId id) const noexcept
{
  return traversalSlot(surfaceSlotById, id);
}

std::optional<uint16_t> TraversalWorkspace::getEdgeSlot(EdgeId id) const noexcept
{
  return traversalSlot(edgeSlotById, id);
}

std::optional<uint16_t> TraversalWorkspace::getCellSlot(CellPathId id) const noexcept
{
  return traversalSlot(cellSlotById, id);
}

void SurfaceTrackingScratch::adoptPlan(std::size_t nOwnedSurfaces, std::size_t nEdges, std::size_t nCells)
{
  mNOwnedSurfaces = nOwnedSurfaces;
  mNEdges = nEdges;
  mNCells = nCells;

  // Group A: one slot per owned surface.
  clearResizeBoundedVector(mClusters, nOwnedSurfaces, getMaybeFrameworkHostResource());
  clearResizeBoundedVector(mUnsortedClusters, nOwnedSurfaces, getMaybeFrameworkHostResource());
  clearResizeBoundedVector(mTrackingFrameInfo, nOwnedSurfaces, getMaybeFrameworkHostResource());
  clearResizeBoundedVector(mClusterExternalIndices, nOwnedSurfaces, mMemoryPool.get());
  clearResizeBoundedVector(mClusterSize, nOwnedSurfaces, mMemoryPool.get());
  clearResizeBoundedVector(mROFramesClusters, nOwnedSurfaces, getMaybeFrameworkHostResource());
  mClusterLabels.assign(nOwnedSurfaces, nullptr);
  mIndexTableUtils.assign(nOwnedSurfaces, IndexTableUtilsCore{});
  clearResizeBoundedVector(mIndexTables, nOwnedSurfaces, getMaybeFrameworkHostResource());
  clearResizeBoundedVector(mUsedClusters, nOwnedSurfaces, getMaybeFrameworkHostResource());
  clearResizeBoundedVector(mNClustersPerROF, nOwnedSurfaces, mMemoryPool.get());
  mMinR.assign(nOwnedSurfaces, std::numeric_limits<float>::max());
  mMaxR.assign(nOwnedSurfaces, std::numeric_limits<float>::lowest());
  mMinZ.assign(nOwnedSurfaces, std::numeric_limits<float>::max());
  mMaxZ.assign(nOwnedSurfaces, std::numeric_limits<float>::lowest());
  clearResizeBoundedVector(mBogusClusters, nOwnedSurfaces, mMemoryPool.get());
  clearResizeBoundedVector(mPositionResolution, nOwnedSurfaces, mMemoryPool.get());

  // Group B: runtime edge and cell counts.
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

void SurfaceTrackingScratch::configureTraversalWorkspaces(std::size_t nIterations)
{
  mTraversalWorkspaces.resize(nIterations);
  for (auto& workspace : mTraversalWorkspaces) {
    workspace.reset(mMemoryPool.get());
  }
}

void SurfaceTrackingScratch::reset()
{
  // Drop the non-owning event view at each event boundary.
  mROFViews = {};
  mROFViewsBySurface.clear();
  mROFLocalLayerBySurface.clear();
  mSourceBySurface.clear();
  mUseUPC = false;

  // Group B.
  deepVectorClear(mTracklets);
  deepVectorClear(mTrackletsLookupTable);
  deepVectorClear(mCells);
  deepVectorClear(mCellsLookupTable);
  deepVectorClear(mCellsNeighbours);
  deepVectorClear(mCellsNeighboursTopology);
  deepVectorClear(mCellsNeighboursLUT);
  deepVectorClear(mEdgePhiCuts);
  deepVectorClear(mEdgeMSAngles);
  for (auto& workspace : mTraversalWorkspaces) {
    workspace.reset(mMemoryPool.get());
  }

  // Group A: allocator-backed storage, always cleared.
  deepVectorClear(mClusterExternalIndices);
  deepVectorClear(mNClustersPerROF);
  deepVectorClear(mBogusClusters);
  deepVectorClear(mPositionResolution);
  deepVectorClear(mClusterSize);
  deepVectorClear(mNTrackletsPerCluster);
  deepVectorClear(mNTrackletsPerClusterSum);

  // Group D.
  deepVectorClear(mTrackletsIndexROF);
  deepVectorClear(mTrackletClusters);
  deepVectorClear(mLines);
  mTotalTracklets = {0, 0};
  mTotalLines = 0;
  mNTotalLowPtVertices = 0;

  // Framework-owned host memory is released by its owner.
  if (!hasFrameworkAllocator()) {
    deepVectorClear(mClusters);
    deepVectorClear(mUsedClusters);
    deepVectorClear(mUnsortedClusters);
    deepVectorClear(mIndexTables);
    deepVectorClear(mTrackingFrameInfo);
    deepVectorClear(mROFramesClusters);
  }

  // Clear MC labels only when MC data is present.
  if (hasMCinformation()) {
    deepVectorClear(mLinesLabels);
    deepVectorClear(mTrackletLabels);
    deepVectorClear(mCellLabels);
  }

  // These are non-owning pointers; reset them without freeing or resizing.
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

  // Host-only allocator bindings.
  initContainers(mClusterExternalIndices);
  initContainers(mNTrackletsPerCluster);
  initContainers(mNTrackletsPerClusterSum);
  initContainers(mNClustersPerROF);
  initVector(mEdgePhiCuts);
  initVector(mEdgeMSAngles);
  initVector(mPositionResolution);
  initContainers(mClusterSize);
  initVector(mBogusClusters);
  initContainers(mTrackletsIndexROF);
  initContainers(mTracklets);
  initContainers(mCells);
  initContainers(mCellsNeighbours);
  initContainers(mCellsLookupTable);
  // MC data (if present).
  initContainers(mLinesLabels);
  initContainers(mTrackletLabels);
  initContainers(mCellLabels);
  for (auto& workspace : mTraversalWorkspaces) {
    workspace.reset(mMemoryPool.get());
  }
  // These may use an external allocator.
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
  // Flat bounded_vector<T> swaps require equal allocators. Outer
  // std::vector<bounded_vector<T>> swaps do not, so they are not checked.
  return flatAllocatorMatches(mBogusClusters, staged.mBogusClusters) &&
         flatAllocatorMatches(mPositionResolution, staged.mPositionResolution) &&
         flatAllocatorMatches(mEdgePhiCuts, staged.mEdgePhiCuts) &&
         flatAllocatorMatches(mEdgeMSAngles, staged.mEdgeMSAngles) &&
         flatArrayAllocatorMatches(mNTrackletsPerCluster, staged.mNTrackletsPerCluster) &&
         flatArrayAllocatorMatches(mNTrackletsPerClusterSum, staged.mNTrackletsPerClusterSum) &&
         flatArrayAllocatorMatches(mTrackletsIndexROF, staged.mTrackletsIndexROF);
}

void SurfaceTrackingScratch::swapLoadedEvent(SurfaceTrackingScratch& other) noexcept
{
  // Swap only data loaded by loadNormalizedSource(); plan sizes, allocators,
  // and edge/cell capacity stay with the live workspace.
  mUnsortedClusters.swap(other.mUnsortedClusters);
  mTrackingFrameInfo.swap(other.mTrackingFrameInfo);
  mClusterExternalIndices.swap(other.mClusterExternalIndices);
  mClusterSize.swap(other.mClusterSize);
  mROFramesClusters.swap(other.mROFramesClusters);
  mClusterLabels.swap(other.mClusterLabels);
  for (int i = 0; i < 2; ++i) {
    mNTrackletsPerCluster[i].swap(other.mNTrackletsPerCluster[i]);
    mNTrackletsPerClusterSum[i].swap(other.mNTrackletsPerClusterSum[i]);
  }
  std::swap(mROFViews, other.mROFViews);
  mROFViewsBySurface.swap(other.mROFViewsBySurface);
  mROFLocalLayerBySurface.swap(other.mROFLocalLayerBySurface);
  mSourceBySurface.swap(other.mSourceBySurface);
  std::swap(mUseUPC, other.mUseUPC);
}

void SurfaceTrackingScratch::swap(SurfaceTrackingScratch& other) noexcept
{
  static_assert(noexcept(std::declval<bounded_vector<float>&>().swap(std::declval<bounded_vector<float>&>())));
  static_assert(noexcept(std::declval<bounded_vector<int>&>().swap(std::declval<bounded_vector<int>&>())));

  std::swap(mNOwnedSurfaces, other.mNOwnedSurfaces);
  std::swap(mNEdges, other.mNEdges);
  std::swap(mNCells, other.mNCells);

  // Outer vectors are always safe to swap; see the header documentation.
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
  mMinZ.swap(other.mMinZ);
  mMaxZ.swap(other.mMaxZ);
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

  // Flat bounded_vector<T> containers require allocatorsMatch(other).
  mBogusClusters.swap(other.mBogusClusters);
  mPositionResolution.swap(other.mPositionResolution);
  mEdgePhiCuts.swap(other.mEdgePhiCuts);
  mEdgeMSAngles.swap(other.mEdgeMSAngles);
  for (std::size_t i = 0; i < mNTrackletsPerCluster.size(); ++i) {
    mNTrackletsPerCluster[i].swap(other.mNTrackletsPerCluster[i]);
    mNTrackletsPerClusterSum[i].swap(other.mNTrackletsPerClusterSum[i]);
    mTrackletsIndexROF[i].swap(other.mTrackletsIndexROF[i]);
  }

  std::swap(mTotalTracklets, other.mTotalTracklets);
  std::swap(mTotalLines, other.mTotalLines);
  std::swap(mNTotalLowPtVertices, other.mNTotalLowPtVertices);
  std::swap(mIsStaggered, other.mIsStaggered);
  std::swap(mROFViews, other.mROFViews);
  mROFViewsBySurface.swap(other.mROFViewsBySurface);
  mROFLocalLayerBySurface.swap(other.mROFLocalLayerBySurface);
  mSourceBySurface.swap(other.mSourceBySurface);

  // Keep memory pools and the external allocator with their owner.
}

// Accessors and runtime-plan operations.

gsl::span<const int> SurfaceTrackingScratch::getROFrameClusters(int layerId) const
{
  return {&mROFramesClusters[layerId][0], static_cast<gsl::span<const int>::size_type>(mROFramesClusters[layerId].size())};
}

gsl::span<o2::its::Cluster> SurfaceTrackingScratch::getClustersOnLayer(int rofId, int layerId)
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<o2::its::Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

gsl::span<const o2::its::Cluster> SurfaceTrackingScratch::getClustersOnLayer(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<const o2::its::Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

gsl::span<uint8_t> SurfaceTrackingScratch::getUsedClustersROF(int rofId, int layerId)
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUsedClusters[layerId][startIdx], static_cast<gsl::span<uint8_t>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

gsl::span<const uint8_t> SurfaceTrackingScratch::getUsedClustersROF(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUsedClusters[layerId][startIdx], static_cast<gsl::span<const uint8_t>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

gsl::span<const o2::its::Cluster> SurfaceTrackingScratch::getClustersPerROFrange(int rofMin, int range, int layerId) const
{
  if (rofMin < 0 || rofMin >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofMin]};
  int endIdx{mROFramesClusters[layerId][o2::gpu::CAMath::Min(rofMin + range, getNrof(layerId))]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<o2::its::Cluster>::size_type>(endIdx - startIdx)};
}

gsl::span<const int> SurfaceTrackingScratch::getROFramesClustersPerROFrange(int rofMin, int range, int layerId) const
{
  int chkdRange{o2::gpu::CAMath::Min(range, getNrof(layerId) - rofMin)};
  return {&mROFramesClusters[layerId][rofMin], static_cast<gsl::span<int>::size_type>(chkdRange)};
}

gsl::span<const int> SurfaceTrackingScratch::getNClustersROFrange(int rofMin, int range, int layerId) const
{
  int chkdRange{o2::gpu::CAMath::Min(range, getNrof(layerId) - rofMin)};
  return {&mNClustersPerROF[layerId][rofMin], static_cast<gsl::span<int>::size_type>(chkdRange)};
}

int SurfaceTrackingScratch::getTotalClustersPerROFrange(int rofMin, int range, int layerId) const
{
  int startIdx{rofMin};
  int endIdx{o2::gpu::CAMath::Min(rofMin + range, getNrof(layerId))};
  return mROFramesClusters[layerId][endIdx] - mROFramesClusters[layerId][startIdx];
}

int SurfaceTrackingScratch::getClusterROF(int iLayer, int iCluster) const
{
  return static_cast<int>(std::lower_bound(mROFramesClusters[iLayer].begin(), mROFramesClusters[iLayer].end(), iCluster + 1) - mROFramesClusters[iLayer].begin() - 1);
}

gsl::span<const o2::its::Cluster> SurfaceTrackingScratch::getUnsortedClustersOnLayer(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUnsortedClusters[layerId][startIdx], static_cast<gsl::span<o2::its::Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

gsl::span<int> SurfaceTrackingScratch::getIndexTable(int rofId, int layer)
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const auto& utils = mIndexTableUtils[layer];
  const int tableSize = utils.getNrowBins() * utils.getNcolBins() + 1;
  return {&mIndexTables[layer][rofId * tableSize], static_cast<gsl::span<int>::size_type>(tableSize)};
}

gsl::span<unsigned char> SurfaceTrackingScratch::getUsedClusters(const int layer)
{
  return {&mUsedClusters[layer][0], static_cast<gsl::span<unsigned char>::size_type>(mUsedClusters[layer].size())};
}

gsl::span<int> SurfaceTrackingScratch::getNTrackletsCluster(int rofId, int combId)
{
  if (rofId < 0 || rofId >= getNrof(1)) {
    return {};
  }
  auto startIdx{mROFramesClusters[1][rofId]};
  return {&mNTrackletsPerCluster[combId][startIdx], static_cast<gsl::span<int>::size_type>(mROFramesClusters[1][rofId + 1] - startIdx)};
}

gsl::span<int> SurfaceTrackingScratch::getExclusiveNTrackletsCluster(int rofId, int combId)
{
  if (rofId < 0 || rofId >= getNrof(1)) {
    return {};
  }
  auto clusStartIdx{mROFramesClusters[1][rofId]};
  return {&mNTrackletsPerClusterSum[combId][clusStartIdx], static_cast<gsl::span<int>::size_type>(mROFramesClusters[1][rofId + 1] - clusStartIdx)};
}

gsl::span<o2::its::Tracklet> SurfaceTrackingScratch::getFoundTracklets(int rofId, int combId)
{
  if (rofId < 0 || rofId >= getNrof(1) || mTracklets[combId].empty()) {
    return {};
  }
  auto startIdx{mNTrackletsPerROF[combId][rofId]};
  return {&mTracklets[combId][startIdx], static_cast<gsl::span<o2::its::Tracklet>::size_type>(mNTrackletsPerROF[combId][rofId + 1] - startIdx)};
}

gsl::span<const o2::its::Tracklet> SurfaceTrackingScratch::getFoundTracklets(int rofId, int combId) const
{
  if (rofId < 0 || rofId >= getNrof(1)) {
    return {};
  }
  auto startIdx{mNTrackletsPerROF[combId][rofId]};
  return {&mTracklets[combId][startIdx], static_cast<gsl::span<o2::its::Tracklet>::size_type>(mNTrackletsPerROF[combId][rofId + 1] - startIdx)};
}

gsl::span<const MCCompLabel> SurfaceTrackingScratch::getLabelsFoundTracklets(int rofId, int combId) const
{
  if (rofId < 0 || rofId >= getNrof(1) || !hasMCinformation()) {
    return {};
  }
  auto startIdx{mNTrackletsPerROF[combId][rofId]};
  return {&mTrackletLabels[combId][startIdx], static_cast<gsl::span<o2::its::Tracklet>::size_type>(mNTrackletsPerROF[combId][rofId + 1] - startIdx)};
}

int SurfaceTrackingScratch::getTotalClusters() const
{
  size_t totalClusters{0};
  for (const auto& clusters : mUnsortedClusters) {
    totalClusters += clusters.size();
  }
  return static_cast<int>(totalClusters);
}

size_t SurfaceTrackingScratch::getNumberOfClusters() const
{
  size_t nClusters{0};
  for (const auto& layer : mClusters) {
    nClusters += layer.size();
  }
  return nClusters;
}

size_t SurfaceTrackingScratch::getNumberOfCells() const
{
  size_t nCells{0};
  for (const auto& layer : mCells) {
    nCells += layer.size();
  }
  return nCells;
}

size_t SurfaceTrackingScratch::getNumberOfTracklets() const
{
  size_t nTracklets{0};
  for (const auto& layer : mTracklets) {
    nTracklets += layer.size();
  }
  return nTracklets;
}

size_t SurfaceTrackingScratch::getNumberOfNeighbours() const
{
  size_t neigh{0};
  for (const auto& l : mCellsNeighbours) {
    neigh += l.size();
  }
  return neigh;
}

size_t SurfaceTrackingScratch::getNumberOfUsedClusters() const
{
  size_t nClusters = 0;
  for (const auto& layer : mUsedClusters) {
    nClusters += std::count(layer.begin(), layer.end(), true);
  }
  return nClusters;
}

const o2::its::TrackingFrameInfo& SurfaceTrackingScratch::getClusterTrackingFrameInfo(int layerId, const o2::its::Cluster& cl) const
{
  return mTrackingFrameInfo[layerId][cl.clusterId];
}

gsl::span<const Vertex> SurfaceTrackingScratch::getPrimaryVertices(const TimeFrame& frame, int layer, int rofId) const
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const auto& view = getROFViews(layer).vertexLookup;
  const auto& entry = view.getVertices(getROFLocalLayer(layer), rofId);
  const auto& vertices = frame.getPrimaryVertices();
  return {&vertices[entry.getFirstEntry()], static_cast<gsl::span<const Vertex>::size_type>(entry.getEntries())};
}

void SurfaceTrackingScratch::computeTrackletsPerROFScans()
{
  for (int iLayer = 0; iLayer < 2; ++iLayer) {
    for (unsigned int iRof{0}; iRof < static_cast<unsigned int>(getNrof(1)); ++iRof) {
      if (isROFEnabled(1, iRof)) {
        mTotalTracklets[iLayer] += mNTrackletsPerROF[iLayer][iRof];
      }
    }
    std::exclusive_scan(mNTrackletsPerROF[iLayer].begin(), mNTrackletsPerROF[iLayer].end(), mNTrackletsPerROF[iLayer].begin(), 0);
    std::exclusive_scan(mNTrackletsPerCluster[iLayer].begin(), mNTrackletsPerCluster[iLayer].end(), mNTrackletsPerClusterSum[iLayer].begin(), 0);
  }
}
void SurfaceTrackingScratch::prepareClusters(const TimeFrame& frame, const TrackingParameters& trkParam, const int maxLayers,
                                             gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements)
{
  struct ClusterHelper {
    float rowCoord;
    float r;
    int bin;
    int ind;
  };
  const std::array<float, 2> beamXY{frame.getBeamX(), frame.getBeamY()};
  const int stopLayer = std::min(maxLayers, static_cast<int>(mNOwnedSurfaces));
  for (int iLayer{0}; iLayer < stopLayer; ++iLayer) {
    const auto& utils = mIndexTableUtils[iLayer];
    const int colBinsCount = utils.getNcolBins();
    std::size_t numBins = 0;
    if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(utils.getNrowBins()),
                                      static_cast<std::size_t>(colBinsCount), numBins) ||
        numBins == std::numeric_limits<std::size_t>::max()) {
      throw std::bad_alloc{};
    }
    const std::size_t stride{numBins + 1};
    bounded_vector<ClusterHelper> cHelper(mMemoryPool.get());
    bounded_vector<int> clsPerBin(numBins, 0, mMemoryPool.get());
    bounded_vector<int> lutPerBin(numBins, 0, mMemoryPool.get());
    for (int rof{0}; rof < getNrof(iLayer); ++rof) {
      if (!isROFEnabled(iLayer, rof)) {
        continue;
      }
      const auto& unsortedClusters{getUnsortedClustersOnLayer(rof, iLayer)};
      const int clustersNum{static_cast<int>(unsortedClusters.size())};
      auto* tableBase = mIndexTables[iLayer].data() + rof * stride;

      cHelper.resize(clustersNum);

      const bool usePhiRBinning = utils.getCoordType() == o2::itsmft::IndexTableCoordType::PhiR;
      for (int iCluster{0}; iCluster < clustersNum; ++iCluster) {
        const o2::its::Cluster& c = unsortedClusters[iCluster];
        const auto& measurement = layerMeasurements[iLayer][c.clusterId];
        ClusterHelper& h = cHelper[iCluster];

        const float x = measurement.position.x - beamXY[0];
        const float y = measurement.position.y - beamXY[1];
        const float z = measurement.position.z;

        const float rowCoord = o2::its::math_utils::computePhi(x, y);
        const float colCoord = usePhiRBinning ? o2::its::math_utils::hypot(x, y) : z;
        int colBin{utils.getColBinIndex(iLayer, colCoord)};
        if (colBin < 0 || colBin >= colBinsCount) {
          colBin = std::clamp(colBin, 0, colBinsCount - 1);
          mBogusClusters[iLayer]++;
        }
        int bin = utils.getBinIndex(colBin, utils.getRowBinIndex(rowCoord));
        h.rowCoord = rowCoord;
        h.r = o2::its::math_utils::hypot(x, y);
        mMinR[iLayer] = o2::gpu::GPUCommonMath::Min(h.r, mMinR[iLayer]);
        mMaxR[iLayer] = o2::gpu::GPUCommonMath::Max(h.r, mMaxR[iLayer]);
        mMinZ[iLayer] = o2::gpu::GPUCommonMath::Min(z, mMinZ[iLayer]);
        mMaxZ[iLayer] = o2::gpu::GPUCommonMath::Max(z, mMaxZ[iLayer]);
        h.bin = bin;
        h.ind = clsPerBin[bin]++;
      }
      std::exclusive_scan(clsPerBin.begin(), clsPerBin.end(), lutPerBin.begin(), 0);

      auto clusters2beSorted{getClustersOnLayer(rof, iLayer)};
      for (int iCluster{0}; iCluster < clustersNum; ++iCluster) {
        const ClusterHelper& h = cHelper[iCluster];
        o2::its::Cluster& c = clusters2beSorted[lutPerBin[h.bin] + h.ind];

        c = unsortedClusters[iCluster];
        const auto& measurement = layerMeasurements[iLayer][c.clusterId];
        c.phi = h.rowCoord;
        c.radius = h.r;
        c.indexTableBinIndex = h.bin;
      }
      std::copy_n(lutPerBin.data(), clsPerBin.size(), tableBase);
      std::fill_n(tableBase + clsPerBin.size(), stride - clsPerBin.size(), clustersNum);

      std::fill(clsPerBin.begin(), clsPerBin.end(), 0);
      cHelper.clear();
    }
  }
}
void SurfaceTrackingScratch::initialise(const TimeFrame& frame, const TrackingParameters& trkParam, const int maxLayers, const int iteration,
                                        const IndexTableUtilsCore& indexTableConfig, TraversalTopologyView topology,
                                        gsl::span<const EdgeId> edgeIds, gsl::span<const CellPathId> cellIds,
                                        gsl::span<const SurfaceId> orderedSurfaces,
                                        gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements)
{
  std::vector<IndexTableUtilsCore> configs(mNOwnedSurfaces, indexTableConfig);
  initialise(frame, trkParam, maxLayers, iteration, configs, topology, edgeIds, cellIds,
             orderedSurfaces, layerMeasurements);
}

void SurfaceTrackingScratch::initialise(const TimeFrame& frame, const TrackingParameters& trkParam, const int maxLayers, const int iteration,
                                        gsl::span<const IndexTableUtilsCore> indexTableConfigs, TraversalTopologyView topology,
                                        gsl::span<const EdgeId> edgeIds, gsl::span<const CellPathId> cellIds,
                                        gsl::span<const SurfaceId> orderedSurfaces,
                                        gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements)
{
  (void)iteration;
  if (orderedSurfaces.size() != mNOwnedSurfaces || indexTableConfigs.size() != mNOwnedSurfaces ||
      edgeIds.size() != mNEdges || cellIds.size() != mNCells ||
      edgeIds.size() > topology.nEdges || cellIds.size() > topology.nPaths ||
      (topology.nEdges != 0 && (topology.edges == nullptr || topology.pathsByFirstEdgeOffsets == nullptr)) ||
      (topology.nPaths != 0 && (topology.paths == nullptr || topology.pathsByFirstEdge == nullptr))) {
    throw std::logic_error{"SurfaceTrackingScratch::initialise(): plan/sparse-topology extent mismatch"};
  }
  const auto surfaceSlot = [&](SurfaceId surface) {
    const auto it = std::find(orderedSurfaces.begin(), orderedSurfaces.end(), surface);
    if (it == orderedSurfaces.end()) {
      throw std::logic_error{"SurfaceTrackingScratch::initialise(): sparse topology surface is not bound"};
    }
    return static_cast<std::size_t>(std::distance(orderedSurfaces.begin(), it));
  };

  for (const auto edgeId : edgeIds) {
    if (!edgeId.isValid() || edgeId.value() >= topology.nEdges) {
      throw std::logic_error{"SurfaceTrackingScratch::initialise(): invalid sparse edge binding"};
    }
    const auto& edge = topology.getEdge(edgeId);
    (void)surfaceSlot(edge.from);
    (void)surfaceSlot(edge.to);
  }
  for (const auto cellId : cellIds) {
    if (!cellId.isValid() || cellId.value() >= topology.nPaths) {
      throw std::logic_error{"SurfaceTrackingScratch::initialise(): invalid sparse cell binding"};
    }
    const auto& path = topology.getPath(cellId);
    if (!path.first.isValid() || !path.second.isValid() ||
        path.first.value() >= topology.nEdges || path.second.value() >= topology.nEdges) {
      throw std::logic_error{"SurfaceTrackingScratch::initialise(): sparse cell references an invalid edge"};
    }
  }

  if (trkParam.PassFlags[IterationStep::FirstPass]) {
    deepVectorClear(mLines);
    deepVectorClear(mLinesLabels);
    clearResizeBoundedVector(mLinesLabels, getNrof(1), mMemoryPool.get());
    mIndexTableUtils.assign(indexTableConfigs.begin(), indexTableConfigs.end());
    clearResizeBoundedVector(mPositionResolution, maxLayers, mMemoryPool.get());
    clearResizeBoundedVector(mBogusClusters, maxLayers, mMemoryPool.get());
    deepVectorClear(mTrackletClusters);
    for (unsigned int iLayer{0}; iLayer < std::min(static_cast<int>(mClusters.size()), maxLayers); ++iLayer) {
      clearResizeBoundedVector(mClusters[iLayer], mUnsortedClusters[iLayer].size(), getMaybeFrameworkHostResource(maxLayers != static_cast<int>(mNOwnedSurfaces)));
      clearResizeBoundedVector(mUsedClusters[iLayer], mUnsortedClusters[iLayer].size(), getMaybeFrameworkHostResource(maxLayers != static_cast<int>(mNOwnedSurfaces)));
      mPositionResolution[iLayer] = o2::gpu::CAMath::Sqrt((0.5f * (trkParam.SystError2Col[iLayer] + trkParam.SystError2Row[iLayer])) + (trkParam.LayerResolution[iLayer] * trkParam.LayerResolution[iLayer]));
    }
    clearResizeBoundedVector(mLines, getNrof(1), mMemoryPool.get());
    clearResizeBoundedVector(mTrackletClusters, getNrof(1), mMemoryPool.get());

    for (std::size_t iLayer{0}; iLayer < mNOwnedSurfaces; ++iLayer) {
      std::size_t indexTableStride = 0;
      if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(mIndexTableUtils[iLayer].getNrowBins()),
                                        static_cast<std::size_t>(mIndexTableUtils[iLayer].getNcolBins()), indexTableStride) ||
          indexTableStride == std::numeric_limits<std::size_t>::max()) {
        throw std::bad_alloc{};
      }
      ++indexTableStride;
      std::size_t indexTableSize = 0;
      if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(getNrof(static_cast<int>(iLayer))), indexTableStride, indexTableSize)) {
        throw std::bad_alloc{};
      }
      clearResizeBoundedVector(mIndexTables[iLayer], indexTableSize, getMaybeFrameworkHostResource());
    }

    std::fill(mMinR.begin(), mMinR.end(), std::numeric_limits<float>::max());
    std::fill(mMaxR.begin(), mMaxR.end(), std::numeric_limits<float>::lowest());
    std::fill(mMinZ.begin(), mMinZ.end(), std::numeric_limits<float>::max());
    std::fill(mMaxZ.begin(), mMaxZ.end(), std::numeric_limits<float>::lowest());
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
  mNTrackletsPerROF.resize(2);
  for (auto& v : mNTrackletsPerROF) {
    v = bounded_vector<int>(getNrof(1) + 1, 0, mMemoryPool.get());
  }
  if (trkParam.PassFlags[IterationStep::RebuildClusterLUT]) {
    prepareClusters(frame, trkParam, maxLayers, layerMeasurements);
  }
  mTotalTracklets = {0, 0};
  if (maxLayers < static_cast<int>(mNOwnedSurfaces)) { // Vertexer only, but in both iterations
    for (int iLayer{0}; iLayer < maxLayers; ++iLayer) {
      deepVectorClear(mUsedClusters[iLayer]);
      clearResizeBoundedVector(mUsedClusters[iLayer], mUnsortedClusters[iLayer].size(), mMemoryPool.get());
    }
  }

  for (std::size_t iLayer{0}; iLayer < mNOwnedSurfaces; ++iLayer) {
    mPositionResolution[iLayer] = o2::gpu::CAMath::Sqrt((0.5f * (trkParam.SystError2Col[iLayer] + trkParam.SystError2Row[iLayer])) + (trkParam.LayerResolution[iLayer] * trkParam.LayerResolution[iLayer]));
  }

  for (int edgeId{0}; edgeId < static_cast<int>(mTracklets.size()); ++edgeId) {
    const auto& edge = topology.getEdge(edgeIds[edgeId]);
    const auto fromSlot = surfaceSlot(edge.from);
    deepVectorClear(mTracklets[edgeId]);
    deepVectorClear(mTrackletLabels[edgeId]);
    deepVectorClear(mTrackletsLookupTable[edgeId]);
    mTrackletsLookupTable[edgeId].resize(mClusters[fromSlot].size() + 1, 0);
  }

  for (int cellId{0}; cellId < static_cast<int>(mCells.size()); ++cellId) {
    deepVectorClear(mCells[cellId]);
    deepVectorClear(mCellsLookupTable[cellId]);
    deepVectorClear(mCellsNeighbours[cellId]);
    deepVectorClear(mCellsNeighboursTopology[cellId]);
    deepVectorClear(mCellsNeighboursLUT[cellId]);
    deepVectorClear(mCellLabels[cellId]);
  }
}

} // namespace o2::itsmft::tracking
