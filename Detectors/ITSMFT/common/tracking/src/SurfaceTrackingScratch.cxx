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
#include <numeric>
#include <new>
#include <stdexcept>
#include <string>

#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/LayerMask.h"
#include "ITStracking/MathUtils.h"

namespace o2::itsmft::tracking
{

// bounded_vector<T> is already brought into this namespace by TimeFrame.h's
// own `template <typename T> using bounded_vector = o2::its::bounded_vector<T>;`
// -- a using-declaration for the same name here would conflict with it.
using o2::its::clearResizeBoundedVector;
using o2::its::deepVectorClear;

namespace constants
{
using namespace o2::its::constants;
} // namespace constants

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
  // the former fixed-layer scratch<NLayers>::resetScratch() exactly.
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
  // and the vector is not resized (mirrors the former fixed-layer scratch<NLayers>'s
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

  // Host-only, mirrors the former fixed-layer scratch<NLayers>::setMemoryPool().
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

// ---------------------------------------------------------------------------
// M6d: the remaining accessor surface TrackerTraits<NLayers>/Tracker<NLayers>/
// SurfacePlanTrackingParticipant needs, ported mechanically from the former
// fixed-layer scratch with
// with every NLayers-bound array index replaced by the equivalent runtime
// vector index -- no algorithm/formula change anywhere in this section.
// ---------------------------------------------------------------------------

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

int SurfaceTrackingScratch::getClusterROF(int iLayer, int iCluster)
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
  const int tableSize = mIndexTableUtils.getNrowBins() * mIndexTableUtils.getNcolBins() + 1;
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

template <int NLayers>
gsl::span<const Vertex> SurfaceTrackingScratch::getPrimaryVertices(const TimeFrame& frame, int layer, int rofId) const
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const auto& view = getROFVertexLookupTableView<NLayers>();
  const auto& entry = view.getVertices(layer, rofId);
  const auto& vertices = frame.getPrimaryVertices();
  return {&vertices[entry.getFirstEntry()], static_cast<gsl::span<const Vertex>::size_type>(entry.getEntries())};
}
template gsl::span<const Vertex> SurfaceTrackingScratch::getPrimaryVertices<ITSNLayers>(const TimeFrame&, int, int) const;
template gsl::span<const Vertex> SurfaceTrackingScratch::getPrimaryVertices<SurfaceTrackingScratch::MFTNLayers>(const TimeFrame&, int, int) const;

template <int NLayers>
void SurfaceTrackingScratch::computeTrackletsPerROFScans()
{
  const auto& maskView = getROFMaskView<NLayers>();
  for (int iLayer = 0; iLayer < 2; ++iLayer) {
    for (unsigned int iRof{0}; iRof < static_cast<unsigned int>(getNrof(1)); ++iRof) {
      if (maskView.isROFEnabled(1, iRof)) {
        mTotalTracklets[iLayer] += mNTrackletsPerROF[iLayer][iRof];
      }
    }
    std::exclusive_scan(mNTrackletsPerROF[iLayer].begin(), mNTrackletsPerROF[iLayer].end(), mNTrackletsPerROF[iLayer].begin(), 0);
    std::exclusive_scan(mNTrackletsPerCluster[iLayer].begin(), mNTrackletsPerCluster[iLayer].end(), mNTrackletsPerClusterSum[iLayer].begin(), 0);
  }
}
template void SurfaceTrackingScratch::computeTrackletsPerROFScans<ITSNLayers>();
template void SurfaceTrackingScratch::computeTrackletsPerROFScans<SurfaceTrackingScratch::MFTNLayers>();

template <int NLayers>
void SurfaceTrackingScratch::initTrackerTopologies(gsl::span<const TrackingParameters> trkParams, int maxLayers)
{
  checkSupportedNLayers<NLayers>();
  if constexpr (NLayers == ITSNLayers) {
    mTrackerTopologiesITS.resize(trkParams.size());
    for (size_t iteration = 0; iteration < trkParams.size(); ++iteration) {
      const int iterationMaxLayers = std::min(maxLayers, trkParams[iteration].NLayers);
      mTrackerTopologiesITS[iteration].init(iterationMaxLayers, trkParams[iteration].MaxHoles, LayerMask{trkParams[iteration].HoleLayerMask});
    }
  } else {
    mTrackerTopologiesMFT.resize(trkParams.size());
    for (size_t iteration = 0; iteration < trkParams.size(); ++iteration) {
      const int iterationMaxLayers = std::min(maxLayers, trkParams[iteration].NLayers);
      mTrackerTopologiesMFT[iteration].init(iterationMaxLayers, trkParams[iteration].MaxHoles, LayerMask{trkParams[iteration].HoleLayerMask});
    }
  }
}
template void SurfaceTrackingScratch::initTrackerTopologies<ITSNLayers>(gsl::span<const TrackingParameters>, int);
template void SurfaceTrackingScratch::initTrackerTopologies<SurfaceTrackingScratch::MFTNLayers>(gsl::span<const TrackingParameters>, int);

template <int NLayers>
void SurfaceTrackingScratch::initDefaultTrackingTopology(const TrackingParameters& trkParam, int maxLayers)
{
  checkSupportedNLayers<NLayers>();
  if constexpr (NLayers == ITSNLayers) {
    mDefaultTrackingTopologyITS.init(maxLayers, trkParam.MaxHoles, LayerMask{trkParam.HoleLayerMask});
  } else {
    mDefaultTrackingTopologyMFT.init(maxLayers, trkParam.MaxHoles, LayerMask{trkParam.HoleLayerMask});
  }
}
template void SurfaceTrackingScratch::initDefaultTrackingTopology<ITSNLayers>(const TrackingParameters&, int);
template void SurfaceTrackingScratch::initDefaultTrackingTopology<SurfaceTrackingScratch::MFTNLayers>(const TrackingParameters&, int);

template <int NLayers>
void SurfaceTrackingScratch::initVertexingTopology(const TrackingParameters& trkParam)
{
  checkSupportedNLayers<NLayers>();
  if constexpr (NLayers == ITSNLayers) {
    mVertexingTopologyITS.init(3, trkParam.MaxHoles, LayerMask{trkParam.HoleLayerMask});
  } else {
    mVertexingTopologyMFT.init(3, trkParam.MaxHoles, LayerMask{trkParam.HoleLayerMask});
  }
}
template void SurfaceTrackingScratch::initVertexingTopology<ITSNLayers>(const TrackingParameters&);
template void SurfaceTrackingScratch::initVertexingTopology<SurfaceTrackingScratch::MFTNLayers>(const TrackingParameters&);

template <int NLayers>
void SurfaceTrackingScratch::prepareClusters(const TimeFrame& frame, const TrackingParameters& trkParam, const int maxLayers,
                                             gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements)
{
  struct ClusterHelper {
    float rowCoord;
    float r;
    int bin;
    int ind;
  };
  const int colBinsCount = mIndexTableUtils.getNcolBins();
  std::size_t numBins = 0;
  if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(mIndexTableUtils.getNrowBins()),
                                    static_cast<std::size_t>(colBinsCount), numBins) ||
      numBins == std::numeric_limits<std::size_t>::max()) {
    throw std::bad_alloc{};
  }
  const std::size_t stride{numBins + 1};
  bounded_vector<ClusterHelper> cHelper(mMemoryPool.get());
  bounded_vector<int> clsPerBin(numBins, 0, mMemoryPool.get());
  bounded_vector<int> lutPerBin(numBins, 0, mMemoryPool.get());
  const std::array<float, 2> beamXY{frame.getBeamX(), frame.getBeamY()};
  const auto& maskView = getROFMaskView<NLayers>();
  for (int iLayer{0}, stopLayer = std::min(trkParam.NLayers, maxLayers); iLayer < stopLayer; ++iLayer) {
    for (int rof{0}; rof < getNrof(iLayer); ++rof) {
      if (!maskView.isROFEnabled(iLayer, rof)) {
        continue;
      }
      const auto& unsortedClusters{getUnsortedClustersOnLayer(rof, iLayer)};
      const int clustersNum{static_cast<int>(unsortedClusters.size())};
      auto* tableBase = mIndexTables[iLayer].data() + rof * stride;

      cHelper.resize(clustersNum);

      const bool useXYBinning = mIndexTableUtils.getCoordType() == o2::itsmft::IndexTableCoordType::XY;
      for (int iCluster{0}; iCluster < clustersNum; ++iCluster) {
        const o2::its::Cluster& c = unsortedClusters[iCluster];
        const auto& measurement = layerMeasurements[iLayer][c.clusterId];
        ClusterHelper& h = cHelper[iCluster];

        const float x = measurement.global.x - (useXYBinning ? 0.f : beamXY[0]);
        const float y = measurement.global.y - (useXYBinning ? 0.f : beamXY[1]);
        const float z = measurement.global.z;

        const float rowCoord = useXYBinning ? measurement.global.y : o2::its::math_utils::computePhi(x, y);
        const float colCoord = useXYBinning ? measurement.global.x : z;
        int colBin{mIndexTableUtils.getColBinIndex(iLayer, colCoord)};
        if (colBin < 0 || colBin >= colBinsCount) {
          colBin = std::clamp(colBin, 0, colBinsCount - 1);
          mBogusClusters[iLayer]++;
        }
        int bin = mIndexTableUtils.getBinIndex(colBin, mIndexTableUtils.getRowBinIndex(rowCoord));
        h.rowCoord = rowCoord;
        h.r = o2::its::math_utils::hypot(x, y);
        mMinR[iLayer] = o2::gpu::GPUCommonMath::Min(h.r, mMinR[iLayer]);
        mMaxR[iLayer] = o2::gpu::GPUCommonMath::Max(h.r, mMaxR[iLayer]);
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
        c.phi = useXYBinning ? o2::its::math_utils::computePhi(measurement.global.x, measurement.global.y) : h.rowCoord;
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
template void SurfaceTrackingScratch::prepareClusters<ITSNLayers>(const TimeFrame&, const TrackingParameters&, int, gsl::span<const gsl::span<const SurfaceMeasurement>>);
template void SurfaceTrackingScratch::prepareClusters<SurfaceTrackingScratch::MFTNLayers>(const TimeFrame&, const TrackingParameters&, int, gsl::span<const gsl::span<const SurfaceMeasurement>>);

template <int NLayers>
void SurfaceTrackingScratch::initialise(const TimeFrame& frame, const TrackingParameters& trkParam, const int maxLayers, const int iteration,
                                        const IndexTableUtilsN& indexTableConfig,
                                        gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements)
{
  checkSupportedNLayers<NLayers>();
  // M6e1: restored to the full three-way ternary the former fixed-layer scratch<NLayers>::
  // initialise() itself uses -- iteration is always a concrete value on every
  // production call path (TrackerTraits<NLayers>::initialiseTimeFrame()
  // never passes UnusedIndex), so the false branch below is still never
  // actually taken; this is fidelity restoration, not a behavior change (see
  // initDefaultTrackingTopology()'s own doc for why it is no longer dead code
  // even though this branch remains unreachable).
  auto& trackingTopologyView = getTrackingTopologyViewMutable<NLayers>();
  if constexpr (NLayers == ITSNLayers) {
    trackingTopologyView = iteration != constants::UnusedIndex ? mTrackerTopologiesITS[iteration].getView() : (maxLayers == 3 ? mVertexingTopologyITS.getView() : mDefaultTrackingTopologyITS.getView());
  } else {
    trackingTopologyView = iteration != constants::UnusedIndex ? mTrackerTopologiesMFT[iteration].getView() : (maxLayers == 3 ? mVertexingTopologyMFT.getView() : mDefaultTrackingTopologyMFT.getView());
  }

  if (trkParam.PassFlags[IterationStep::FirstPass]) {
    deepVectorClear(mLines);
    deepVectorClear(mLinesLabels);
    clearResizeBoundedVector(mLinesLabels, getNrof(1), mMemoryPool.get());
    mIndexTableUtils = indexTableConfig;
    clearResizeBoundedVector(mPositionResolution, trkParam.NLayers, mMemoryPool.get());
    clearResizeBoundedVector(mBogusClusters, trkParam.NLayers, mMemoryPool.get());
    deepVectorClear(mTrackletClusters);
    for (unsigned int iLayer{0}; iLayer < std::min(static_cast<int>(mClusters.size()), maxLayers); ++iLayer) {
      clearResizeBoundedVector(mClusters[iLayer], mUnsortedClusters[iLayer].size(), getMaybeFrameworkHostResource(maxLayers != static_cast<int>(mNOwnedSurfaces)));
      clearResizeBoundedVector(mUsedClusters[iLayer], mUnsortedClusters[iLayer].size(), getMaybeFrameworkHostResource(maxLayers != static_cast<int>(mNOwnedSurfaces)));
      mPositionResolution[iLayer] = o2::gpu::CAMath::Sqrt((0.5f * (trkParam.SystError2Col[iLayer] + trkParam.SystError2Row[iLayer])) + (trkParam.LayerResolution[iLayer] * trkParam.LayerResolution[iLayer]));
    }
    clearResizeBoundedVector(mLines, getNrof(1), mMemoryPool.get());
    clearResizeBoundedVector(mTrackletClusters, getNrof(1), mMemoryPool.get());

    std::size_t indexTableStride = 0;
    if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(mIndexTableUtils.getNrowBins()),
                                      static_cast<std::size_t>(mIndexTableUtils.getNcolBins()), indexTableStride) ||
        indexTableStride == std::numeric_limits<std::size_t>::max()) {
      throw std::bad_alloc{};
    }
    ++indexTableStride;
    for (std::size_t iLayer{0}; iLayer < mNOwnedSurfaces; ++iLayer) {
      std::size_t indexTableSize = 0;
      if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(getNrof(static_cast<int>(iLayer))), indexTableStride, indexTableSize)) {
        throw std::bad_alloc{};
      }
      clearResizeBoundedVector(mIndexTables[iLayer], indexTableSize, getMaybeFrameworkHostResource());
    }

    std::fill(mMinR.begin(), mMinR.end(), std::numeric_limits<float>::max());
    std::fill(mMaxR.begin(), mMaxR.end(), std::numeric_limits<float>::lowest());
  }
  clearResizeBoundedVector(mCells, trackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsLookupTable, trackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighbours, trackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursTopology, trackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellsNeighboursLUT, trackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mCellLabels, trackingTopologyView.nCells, mMemoryPool.get());
  clearResizeBoundedVector(mTracklets, trackingTopologyView.nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletLabels, trackingTopologyView.nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTrackletsLookupTable, trackingTopologyView.nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTransitionPhiCuts, trackingTopologyView.nTransitions, mMemoryPool.get());
  clearResizeBoundedVector(mTransitionMSAngles, trackingTopologyView.nTransitions, mMemoryPool.get());
  mNTrackletsPerROF.resize(2);
  for (auto& v : mNTrackletsPerROF) {
    v = bounded_vector<int>(getNrof(1) + 1, 0, mMemoryPool.get());
  }
  if (trkParam.PassFlags[IterationStep::RebuildClusterLUT]) {
    prepareClusters<NLayers>(frame, trkParam, maxLayers, layerMeasurements);
  }
  mTotalTracklets = {0, 0};
  if (maxLayers < trkParam.NLayers) { // Vertexer only, but in both iterations
    for (int iLayer{0}; iLayer < maxLayers; ++iLayer) {
      deepVectorClear(mUsedClusters[iLayer]);
      clearResizeBoundedVector(mUsedClusters[iLayer], mUnsortedClusters[iLayer].size(), mMemoryPool.get());
    }
  }

  for (std::size_t iLayer{0}; iLayer < mNOwnedSurfaces; ++iLayer) {
    mPositionResolution[iLayer] = o2::gpu::CAMath::Sqrt((0.5f * (trkParam.SystError2Col[iLayer] + trkParam.SystError2Row[iLayer])) + (trkParam.LayerResolution[iLayer] * trkParam.LayerResolution[iLayer]));
  }

  for (int transitionId{0}; transitionId < static_cast<int>(mTracklets.size()); ++transitionId) {
    const auto& transition = trackingTopologyView.getTransition(transitionId);
    deepVectorClear(mTracklets[transitionId]);
    deepVectorClear(mTrackletLabels[transitionId]);
    deepVectorClear(mTrackletsLookupTable[transitionId]);
    mTrackletsLookupTable[transitionId].resize(mClusters[transition.fromLayer].size() + 1, 0);
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
template void SurfaceTrackingScratch::initialise<ITSNLayers>(const TimeFrame&, const TrackingParameters&, int, int, const IndexTableUtilsN&, gsl::span<const gsl::span<const SurfaceMeasurement>>);
template void SurfaceTrackingScratch::initialise<SurfaceTrackingScratch::MFTNLayers>(const TimeFrame&, const TrackingParameters&, int, int, const IndexTableUtilsN&, gsl::span<const gsl::span<const SurfaceMeasurement>>);

#ifndef GPUCA_GPUCODE
namespace
{
class NormalizedBackfillAllocatorMismatch final : public std::logic_error
{
 public:
  explicit NormalizedBackfillAllocatorMismatch(int layer)
    : std::logic_error("SurfaceTrackingScratch::loadNormalizedSource(): staged/live memory-resource mismatch on layer " + std::to_string(layer))
  {
  }
};
} // namespace

LoadSourcesResult SurfaceTrackingScratch::loadNormalizedSource(
  TimeFrame& frame,
  const ClusterDecoder& decoder,
  const o2::InteractionRecord& origin,
  const ROFTimingConfig& timing,
  gsl::span<const itsmft::CompClusterExt> clusters,
  gsl::span<const unsigned char> patterns,
  gsl::span<const o2::itsmft::ROFRecord> rofs,
  const itsmft::TopologyDictionary* dictionary,
  const dataformats::MCTruthContainer<MCCompLabel>* labels,
  o2::detectors::DetID::ID detId,
  gsl::span<const SurfaceId> orderedSurfaces,
  SurfaceCatalogView catalogView,
  bool applySysErrors)
{
  // M6e2: this scratch is now shared by ITS too (previously MFT-only) --
  // matches the former fixed-layer scratch<NLayers>'s own ITS-or-MFT preflight.
  constexpr ClusterSourceId kSourceId{0};
  if (detId != o2::detectors::DetID::MFT && detId != o2::detectors::DetID::ITS) {
    return {MultiSourceLoadError::UnsupportedDetector, kSourceId};
  }
  if (catalogView.surfaces == nullptr || catalogView.nSurfaces == 0) {
    return {MultiSourceLoadError::SurfaceCatalogNotConfigured, kSourceId};
  }
  if (orderedSurfaces.size() != mNOwnedSurfaces) {
    return {MultiSourceLoadError::InvalidLayerMapping, kSourceId};
  }
  SurfaceMask mappedSurfaceSeen{};
  for (const auto& surfaceId : orderedSurfaces) {
    if (!surfaceId.isValid() || surfaceId.value() >= catalogView.nSurfaces) {
      return {MultiSourceLoadError::InvalidLayerMapping, kSourceId};
    }
    if (mappedSurfaceSeen.has(surfaceId)) {
      return {MultiSourceLoadError::InvalidLayerMapping, kSourceId};
    }
    mappedSurfaceSeen.set(surfaceId);
  }
  for (const auto& surfaceId : orderedSurfaces) {
    if (catalogView.getSurface(surfaceId).detectorId != static_cast<uint8_t>(detId)) {
      return {MultiSourceLoadError::DetectorSurfaceMismatch, kSourceId};
    }
  }

  const gsl::span<const SurfaceId> layerToSurface = orderedSurfaces;
  const std::size_t nOwnedSurfaces = orderedSurfaces.size();

  MultiSourceFrame staged;
  ClusterSourceInput src;
  src.id = kSourceId;
  src.detector = detId;
  src.clusters = clusters;
  src.patterns = patterns;
  src.rofs = rofs;
  src.dictionary = dictionary;
  src.labels = labels;
  src.layerToSurface = layerToSurface;
  src.timing = timing;
  src.decoder = &decoder;
  src.applySysErrors = applySysErrors;

  const auto result = loadSources(staged, catalogView, gsl::span<const ClusterSourceInput>(&src, 1), origin);
  if (!result.ok()) {
    return result;
  }

  const bool isMFT = (detId == o2::detectors::DetID::MFT);
  auto* pool = mMemoryPool.get();
  const auto nROFs = static_cast<size_t>(rofs.size());

  std::vector<bounded_vector<o2::its::Cluster>> stagedUnsortedClusters;
  std::vector<bounded_vector<o2::its::TrackingFrameInfo>> stagedTrackingFrameInfo;
  std::vector<bounded_vector<int>> stagedClusterExternalIndices;
  std::vector<bounded_vector<uint8_t>> stagedClusterSize;
  std::vector<bounded_vector<int>> stagedROFramesClusters;
  std::vector<const dataformats::MCTruthContainer<MCCompLabel>*> stagedClusterLabels(nOwnedSurfaces, nullptr);
  stagedUnsortedClusters.reserve(nOwnedSurfaces);
  stagedTrackingFrameInfo.reserve(nOwnedSurfaces);
  stagedClusterExternalIndices.reserve(nOwnedSurfaces);
  stagedClusterSize.reserve(nOwnedSurfaces);
  stagedROFramesClusters.reserve(nOwnedSurfaces);

  for (std::size_t layer = 0; layer < nOwnedSurfaces; ++layer) {
    const auto measurements = staged.getSurfaceMeasurements(layerToSurface[layer]);

    auto* mr = getMaybeFrameworkHostResource();
    auto* unsortedClustersMr = mr != nullptr ? mr : mUnsortedClusters[layer].get_allocator().resource();
    auto* trackingFrameInfoMr = mr != nullptr ? mr : mTrackingFrameInfo[layer].get_allocator().resource();
    auto* clusterExternalIndicesMr = pool != nullptr ? pool : mClusterExternalIndices[layer].get_allocator().resource();
    auto* clusterSizeMr = pool != nullptr ? pool : mClusterSize[layer].get_allocator().resource();
    auto* rofFramesClustersMr = mr != nullptr ? mr : mROFramesClusters[layer].get_allocator().resource();

    stagedUnsortedClusters.emplace_back(std::pmr::polymorphic_allocator<o2::its::Cluster>{unsortedClustersMr});
    stagedTrackingFrameInfo.emplace_back(std::pmr::polymorphic_allocator<o2::its::TrackingFrameInfo>{trackingFrameInfoMr});
    stagedClusterExternalIndices.emplace_back(std::pmr::polymorphic_allocator<int>{clusterExternalIndicesMr});
    stagedClusterSize.emplace_back(measurements.size(), uint8_t{0}, std::pmr::polymorphic_allocator<uint8_t>{clusterSizeMr});
    stagedROFramesClusters.emplace_back(nROFs + 1, 0, std::pmr::polymorphic_allocator<int>{rofFramesClustersMr});

    size_t mi{0};
    for (const auto& m : measurements) {
      o2::its::TrackingFrameInfo tfInfo;
      if (isMFT) {
        // Recreate the established synthetic legacy MFT representation
        // (TrackingFrameInfoAdapters.h::makeTrackingFrameInfo<MFT>) from the
        // normalized global position and row/column covariance -- ported
        // byte-for-byte from the former fixed-layer scratch<NLayers>::
        // loadNormalizedSource() (the former fixed-layer scratch.cxx).
        tfInfo = o2::its::TrackingFrameInfo{
          m.global.x, m.global.y, m.global.z,
          m.global.x, 0.f,
          std::array<float, 2>{m.global.y, m.global.z},
          std::array<float, 3>{m.covariance.uu, m.covariance.uv, m.covariance.vv}};
      } else {
        // ITS: as above, ported byte-for-byte from
        // the former fixed-layer scratch<NLayers>::loadNormalizedSource().
        tfInfo = o2::its::TrackingFrameInfo{
          m.global.x, m.global.y, m.global.z,
          m.frame.q, m.frame.frameAngle,
          std::array<float, 2>{m.frame.u, m.frame.v},
          std::array<float, 3>{m.covariance.uu, m.covariance.uv, m.covariance.vv}};
      }
      stagedTrackingFrameInfo[layer].push_back(tfInfo);
      stagedUnsortedClusters[layer].emplace_back(m.global.x, m.global.y, m.global.z, static_cast<int>(stagedUnsortedClusters[layer].size()));
      stagedClusterExternalIndices[layer].push_back(static_cast<int>(m.cluster.index));
      stagedClusterSize[layer][mi] = static_cast<uint8_t>(std::clamp(m.shape.nPixels, 0u, 255u));
      ++mi;
    }

    size_t mj{0};
    for (size_t r = 0; r < nROFs; ++r) {
      while (mj < measurements.size() && measurements[mj].sourceROF == static_cast<uint32_t>(r)) {
        ++mj;
      }
      stagedROFramesClusters[layer][r + 1] = static_cast<int>(mj);
    }

    stagedClusterLabels[layer] = labels;
  }

  const size_t nClustersLayer1 = stagedUnsortedClusters[1].size();
  auto* nTrackletsPerClusterMr = pool != nullptr ? pool : mNTrackletsPerCluster[0].get_allocator().resource();
  auto* nTrackletsPerClusterSumMr = pool != nullptr ? pool : mNTrackletsPerClusterSum[0].get_allocator().resource();
  std::array<bounded_vector<int>, 2> stagedNTrackletsPerCluster{
    bounded_vector<int>(nClustersLayer1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterMr}),
    bounded_vector<int>(nClustersLayer1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterMr})};
  std::array<bounded_vector<int>, 2> stagedNTrackletsPerClusterSum{
    bounded_vector<int>(nClustersLayer1 + 1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterSumMr}),
    bounded_vector<int>(nClustersLayer1 + 1, 0, std::pmr::polymorphic_allocator<int>{nTrackletsPerClusterSumMr})};

  for (std::size_t layer = 0; layer < nOwnedSurfaces; ++layer) {
    if (mUnsortedClusters[layer].get_allocator().resource() != stagedUnsortedClusters[layer].get_allocator().resource() ||
        mTrackingFrameInfo[layer].get_allocator().resource() != stagedTrackingFrameInfo[layer].get_allocator().resource() ||
        mClusterExternalIndices[layer].get_allocator().resource() != stagedClusterExternalIndices[layer].get_allocator().resource() ||
        mClusterSize[layer].get_allocator().resource() != stagedClusterSize[layer].get_allocator().resource() ||
        mROFramesClusters[layer].get_allocator().resource() != stagedROFramesClusters[layer].get_allocator().resource()) {
      throw NormalizedBackfillAllocatorMismatch(static_cast<int>(layer));
    }
  }
  for (int i = 0; i < 2; ++i) {
    if (mNTrackletsPerCluster[i].get_allocator().resource() != stagedNTrackletsPerCluster[i].get_allocator().resource() ||
        mNTrackletsPerClusterSum[i].get_allocator().resource() != stagedNTrackletsPerClusterSum[i].get_allocator().resource()) {
      throw NormalizedBackfillAllocatorMismatch(-1);
    }
  }

  static_assert(noexcept(std::declval<bounded_vector<o2::its::Cluster>&>().swap(std::declval<bounded_vector<o2::its::Cluster>&>())));
  static_assert(noexcept(std::declval<bounded_vector<int>&>().swap(std::declval<bounded_vector<int>&>())));

  frame.commitNormalizedFrame(std::move(staged));
  for (std::size_t layer = 0; layer < nOwnedSurfaces; ++layer) {
    mUnsortedClusters[layer].swap(stagedUnsortedClusters[layer]);
    mTrackingFrameInfo[layer].swap(stagedTrackingFrameInfo[layer]);
    mClusterExternalIndices[layer].swap(stagedClusterExternalIndices[layer]);
    mClusterSize[layer].swap(stagedClusterSize[layer]);
    mROFramesClusters[layer].swap(stagedROFramesClusters[layer]);
    mClusterLabels[layer] = stagedClusterLabels[layer];
  }
  for (int i = 0; i < 2; ++i) {
    mNTrackletsPerCluster[i].swap(stagedNTrackletsPerCluster[i]);
    mNTrackletsPerClusterSum[i].swap(stagedNTrackletsPerClusterSum[i]);
  }

  return result;
}
#endif

} // namespace o2::itsmft::tracking
