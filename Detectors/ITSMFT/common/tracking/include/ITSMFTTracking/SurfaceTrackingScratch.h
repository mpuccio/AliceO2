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
/// \file SurfaceTrackingScratch.h
/// \brief Runtime-plan-owned, detector-neutral CA workspace.
///
/// The workspace sizes all host-side event storage from adoptPlan(): one slot
/// per ordered owned surface and one compact slot per bound sparse transition
/// or cell. The fixed-capacity device value types remain independent of this
/// host workspace.
///
/// Timing, ROF assignment, and masks are supplied as one non-owning
/// RuntimeROFViews event context. Fixed-capacity table builders stay with the
/// application adapter that owns their lifetime; this class has no detector
/// table or topology object to select.
///
/// SparseTrackingTopologyView and SurfacePlanBinding ids are passed to
/// initialise() explicitly. Traversal order is therefore the binding's
/// ordered-surface/transition/cell order, never a numeric SurfaceId order.
///
/// This type never owns a plan or binding. It borrows the caller's
/// SparseTrackingTopologyView and runtime ROF context for the current event;
/// the event owner/adapter owns raw ROFs and the event lifecycle.
#ifndef ALICEO2_ITSMFT_TRACKING_SURFACETRACKINGSCRATCH_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACETRACKINGSCRATCH_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <vector>

#include <gsl/gsl>

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITS/Vertex.h"

#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/ROFViews.h"
#include "ITSMFTTracking/SparseTrackingTopology.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/ClusterLines.h"
#include "ITStracking/ExternalAllocator.h"
#include "ITStracking/Tracklet.h"
#ifndef GPUCA_GPUCODE
#include <optional>
#include "ITSMFTTracking/SurfaceCatalogView.h"
#endif
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include "DetectorsCommonDataFormats/DetID.h"

namespace o2::itsmft
{
class CompClusterExt;
class TopologyDictionary;
class ROFRecord;
} // namespace o2::itsmft

namespace o2::itsmft::tracking
{

class MultiSourceTimeFrameLoader;
class ClusterDecoder;

/// Non-templated, detector-neutral CA working state. Detector-specific table
/// builders and raw ROF ownership stay at the application boundary; this
/// class retains only the runtime views needed by the current event.
class SurfaceTrackingScratch
{
 private:
  // ---- Group E: memory/allocator/device plumbing ----
  // Declared first so pool owners outlive every allocator-backed member below.
  // C++ destroys non-static data members in reverse declaration order, so
  // C++ destroys non-static data members in reverse declaration order, so
  // declaring these two pool owners first guarantees they are destroyed
  // *last*, after every pmr/bounded_vector member declared below that could
  // still hold memory allocated through them. An earlier revision of this
  // header declared these last on the theory that the underlying
  // BoundedMemoryResource is always externally owned; that is not
  // guaranteed (a caller may hand this class the sole owning shared_ptr, as
  // this type's own focused tests do), and declaring them last caused a
  // real reproducible segfault in ~SurfaceTrackingScratch() when the pool's
  // last reference was released before the vectors allocated through it
  // were destroyed. This ordering is part of the allocator ownership contract.
  std::shared_ptr<o2::its::BoundedMemoryResource> mExtMemoryPool;
  std::shared_ptr<o2::its::BoundedMemoryResource> mMemoryPool;
  o2::its::ExternalAllocator* mExternalAllocator{nullptr};
  bool mIsStaggered{false};
  // Flagged in the M6 design note as possibly dead (no production
  // read/write site found by that audit); carried over structurally rather
  // than dropped, since this milestone does not re-run that verification.
  o2::its::bounded_vector<std::array<float, 2>> mPValphaX;

 public:
  // M6d: LoadTargetImpl-equivalent staging code (MultiSourceTimeFrameLoader.cxx)
  // needs direct access to mExternalAllocator/mExtMemoryPool to preserve
  // allocator identity across stage()/commit(), exactly mirroring
  // the workspace's allocator ownership contract.
  friend class MultiSourceTimeFrameLoader;

  SurfaceTrackingScratch() = default;
  ~SurfaceTrackingScratch() = default;
  SurfaceTrackingScratch(const SurfaceTrackingScratch&) = delete;
  SurfaceTrackingScratch& operator=(const SurfaceTrackingScratch&) = delete;
  SurfaceTrackingScratch(SurfaceTrackingScratch&&) = delete;
  SurfaceTrackingScratch& operator=(SurfaceTrackingScratch&&) = delete;

  /// Sizes every Group A container to `nOwnedSurfaces` (one slot per owned
  /// surface, and every Group B container to `nTransitions`/`nCells`
  /// (already-runtime sparse-topology counts). Precondition: setMemoryPool() has
  /// already been called -- this never allocates through a null resource
  /// silently, it inherits whichever resource the owner already configured.
  /// Group D is not sized here (never plan-sized -- see the file doc); Group
  /// E is unaffected (memory/allocator plumbing is set up independently via
  /// setMemoryPool()/setFrameworkAllocator(), not by plan adoption).
  void adoptPlan(std::size_t nOwnedSurfaces, std::size_t nTransitions, std::size_t nCells);

  std::size_t getNOwnedSurfaces() const noexcept { return mNOwnedSurfaces; }
  std::size_t getNTransitions() const noexcept { return mNTransitions; }
  std::size_t getNCells() const noexcept { return mNCells; }

  /// Clears scratch-owned working state in place -- mirrors
  /// the previous workspace reset contract (same deep-clear-vs-framework-
  /// allocator-skip, same MC-info-conditional
  /// clearing, same non-owning mClusterLabels re-nulling instead of
  /// deep-clear). Never touches a TimeFrame, even implicitly, and never
  /// changes the adopted plan sizing (getNOwnedSurfaces()/getNTransitions()/
  /// getNCells() are unaffected -- only each container's *contents* are
  /// cleared; the operation never shrinks the plan-sized outer
  /// arrays either). Matches resetScratch()'s own name too now (M6d wires
  /// this in where production code calls it).
  void reset();
  void resetScratch() { reset(); }

  /// memory management -- Group E: reseat every
  /// allocator-backed container onto the new resource via a deep clear, so
  /// every subsequent allocation happens through the caller's pool.
  void setMemoryPool(std::shared_ptr<o2::its::BoundedMemoryResource> pool);
  auto& getMemoryPool() const noexcept { return mMemoryPool; }
  void setFrameworkAllocator(o2::its::ExternalAllocator* ext);
  auto getFrameworkAllocator() const noexcept { return mExternalAllocator; }
  bool hasFrameworkAllocator() const noexcept { return mExternalAllocator != nullptr; }
  std::pmr::memory_resource* getMaybeFrameworkHostResource(bool forceHost = false)
  {
    return (hasFrameworkAllocator() && !forceHost) ? mExtMemoryPool.get() : mMemoryPool.get();
  }

  void setIsStaggered(bool b) noexcept { mIsStaggered = b; }
  bool isStaggered() const noexcept { return mIsStaggered; }

  bool hasMCinformation() const noexcept { return !mClusterLabels.empty() && mClusterLabels[0] != nullptr; }

  /// Staging/swap support for the atomic loader migration -- uses the same
  /// stage-then-commit discipline as the normalized loader, generalized to
  /// every plan-sized container here.
  ///
  /// allocatorsMatch() is the precondition a caller must check before
  /// swap(): every *flat* bounded_vector member here (the ones swap()
  /// exchanges via bounded_vector::swap(), which is only well-defined when
  /// both sides' allocators compare equal -- std::pmr::polymorphic_allocator
  /// neither propagates on swap nor is ever always_equal) must already share
  /// its counterpart's memory-resource pointer. The per-owned-surface/
  /// per-transition/per-cell *outer* containers (std::vector<bounded_vector<T>>)
  /// carry no such precondition: swapping the outer std::vector only
  /// exchanges its own internal pointers/size/capacity and never touches an
  /// inner bounded_vector's allocator, so those are always safe to swap
  /// regardless of allocator identity.
  bool allocatorsMatch(const SurfaceTrackingScratch& staged) const noexcept;

  /// Precondition: allocatorsMatch(other) (checked by the caller before
  /// commit, mirroring loadNormalizedSource()'s own invariant gate -- not
  /// re-checked here so this can stay noexcept). Never swaps mMemoryPool/
  /// mExtMemoryPool/mExternalAllocator: allocator identity is owner-bound,
  /// not staged data, exactly as loadNormalizedSource() never swaps
  /// mMemoryPool itself, only the vectors allocated through it.
  void swap(SurfaceTrackingScratch& other) noexcept;

  // ---- read-in data: every loop uses the runtime ordered-surface span.
  // detId preflight accepts both common-CA detectors.
#ifndef GPUCA_GPUCODE
  LoadSourcesResult loadNormalizedSource(TimeFrame& frame,
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
                                         bool applySysErrors = true);
#endif

  int getTotalClusters() const;
  bool empty() const { return getTotalClusters() == 0; }
  int getSortedIndex(int rofId, int layer, int idx) const { return mROFramesClusters[layer][rofId] + idx; }
  int getSortedStartIndex(const int rofId, const int layer) const { return mROFramesClusters[layer][rofId]; }
  int getNrof(int layer) const { return static_cast<int>(mROFramesClusters[layer].size()) - 1; }

  auto& getMinRs() { return mMinR; }
  auto& getMaxRs() { return mMaxR; }
  float getMinR(int layer) const { return mMinR[layer]; }
  float getMaxR(int layer) const { return mMaxR[layer]; }
  float getTransitionPhiCut(int transitionId) const { return mTransitionPhiCuts[transitionId]; }
  float getTransitionMSAngle(int transitionId) const { return mTransitionMSAngles[transitionId]; }
  auto& getTransitionPhiCuts() { return mTransitionPhiCuts; }
  auto& getTransitionMSAngles() { return mTransitionMSAngles; }
  float getPositionResolution(int layer) const { return mPositionResolution[layer]; }
  auto& getPositionResolutions() { return mPositionResolution; }

  gsl::span<o2::its::Cluster> getClustersOnLayer(int rofId, int layerId);
  gsl::span<const o2::its::Cluster> getClustersOnLayer(int rofId, int layerId) const;
  gsl::span<const o2::its::Cluster> getClustersPerROFrange(int rofMin, int range, int layerId) const;
  gsl::span<const o2::its::Cluster> getUnsortedClustersOnLayer(int rofId, int layerId) const;
  gsl::span<uint8_t> getUsedClustersROF(int rofId, int layerId);
  gsl::span<const uint8_t> getUsedClustersROF(int rofId, int layerId) const;
  gsl::span<const int> getROFramesClustersPerROFrange(int rofMin, int range, int layerId) const;
  gsl::span<const int> getROFrameClusters(int layerId) const;
  gsl::span<const int> getNClustersROFrange(int rofMin, int range, int layerId) const;
  gsl::span<int> getIndexTable(int rofId, int layerId);
  const auto& getTrackingFrameInfoOnLayer(int layerId) const { return mTrackingFrameInfo[layerId]; }

  // Navigation and event timing views. Fixed-capacity detector tables are
  // owned by their application adapters; common scratch retains one
  // non-owning runtime view for the current event.
  const auto& getIndexTableUtils() const { return mIndexTableUtils; }
  void setROFViews(RuntimeROFViews views) noexcept
  {
    mROFViews = views;
    mUseUPC = false;
  }
  const RuntimeROFViews& getROFViews() const noexcept { return mROFViews; }
  const RuntimeROFOverlapView& getROFOverlapView() const noexcept { return mROFViews.overlap; }
  const RuntimeROFVertexLookupView& getROFVertexLookupView() const noexcept { return mROFViews.vertexLookup; }
  const RuntimeROFMaskView& getROFMaskView() const noexcept { return mUseUPC ? mROFViews.upcMask : mROFViews.mask; }
  void useUPCMask() noexcept { mUseUPC = true; }
  gsl::span<const Vertex> getPrimaryVertices(const TimeFrame& frame, int layer, int rofId) const;

  const o2::its::TrackingFrameInfo& getClusterTrackingFrameInfo(int layerId, const o2::its::Cluster& cl) const;
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const o2::its::Cluster& cl) const { return getClusterLabels(layerId, cl.clusterId); }
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const int clId) const { return mClusterLabels[(mIsStaggered ? layerId : 0)]->getLabels(mClusterExternalIndices[layerId][clId]); }
  int getClusterExternalIndex(int layerId, const int clId) const { return mClusterExternalIndices[layerId][clId]; }
  int getClusterSize(int layer, int clusterId) const { return mClusterSize[layer][clusterId]; }
  void setClusterSize(int layer, o2::its::bounded_vector<uint8_t>& v) { mClusterSize[layer] = std::move(v); }

  auto& getTrackletsLabel(int layer) { return mTrackletLabels[layer]; }
  auto& getCellsLabel(int layer) { return mCellLabels[layer]; }

  void initialise(const TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers, int iteration,
                  const IndexTableUtilsCore& indexTableConfig, SparseTrackingTopologyView topology,
                  gsl::span<const TransitionId> transitionIds, gsl::span<const CellTopologyId> cellIds,
                  gsl::span<const SurfaceId> orderedSurfaces,
                  gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements);

  bool isClusterUsed(int layer, int clusterId) const { return mUsedClusters[layer][clusterId]; }
  void markUsedCluster(int layer, int clusterId) { mUsedClusters[layer][clusterId] = true; }
  gsl::span<unsigned char> getUsedClusters(const int layer);

  auto& getTracklets() { return mTracklets; }
  auto& getTrackletsLookupTable() { return mTrackletsLookupTable; }

  auto& getClusters() { return mClusters; }
  auto& getUnsortedClusters() { return mUnsortedClusters; }
  const auto& getUnsortedClusters() const { return mUnsortedClusters; }
  int getClusterROF(int iLayer, int iCluster) const;
  auto& getCells() { return mCells; }

  auto& getCellsLookupTable() { return mCellsLookupTable; }
  auto& getCellsNeighbours() { return mCellsNeighbours; }
  auto& getCellsNeighboursTopology() { return mCellsNeighboursTopology; }
  auto& getCellsNeighboursLUT() { return mCellsNeighboursLUT; }
  auto& getLinesLabel(const int rofId) { return mLinesLabels[rofId]; }

  size_t getNumberOfClusters() const;
  size_t getNumberOfCells() const;
  size_t getNumberOfTracklets() const;
  size_t getNumberOfNeighbours() const;
  size_t getNumberOfUsedClusters() const;

  int hasBogusClusters() const { return std::accumulate(mBogusClusters.begin(), mBogusClusters.end(), 0); }

  template <typename... T>
  void addClusterToLayer(int layer, T&&... args)
  {
    mUnsortedClusters[layer].emplace_back(std::forward<T>(args)...);
  }
  template <typename... T>
  void addTrackingFrameInfoToLayer(int layer, T&&... args)
  {
    mTrackingFrameInfo[layer].emplace_back(std::forward<T>(args)...);
  }
  void addClusterExternalIndexToLayer(int layer, const int idx) { mClusterExternalIndices[layer].push_back(idx); }

  // ---- Group A: per-(owned-surface) cluster/index-table cache ----
  std::vector<o2::its::bounded_vector<o2::its::Cluster>> mClusters;
  std::vector<o2::its::bounded_vector<o2::its::Cluster>> mUnsortedClusters;
  std::vector<o2::its::bounded_vector<o2::its::TrackingFrameInfo>> mTrackingFrameInfo;
  std::vector<o2::its::bounded_vector<int>> mClusterExternalIndices;
  std::vector<o2::its::bounded_vector<uint8_t>> mClusterSize;
  std::vector<o2::its::bounded_vector<int>> mROFramesClusters;
  std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*> mClusterLabels;
  std::vector<o2::its::bounded_vector<int>> mIndexTables;
  std::vector<o2::its::bounded_vector<uint8_t>> mUsedClusters;
  std::vector<o2::its::bounded_vector<int>> mNClustersPerROF;
  std::vector<float> mMinR;
  std::vector<float> mMaxR;
  o2::its::bounded_vector<int> mBogusClusters;
  o2::its::bounded_vector<float> mPositionResolution;
  // Not per-owned-surface: fixed at the two tracklet combinations (positions
  // 0-1, 1-2), sized from cluster count at load time.
  std::array<o2::its::bounded_vector<int>, 2> mNTrackletsPerCluster;
  std::array<o2::its::bounded_vector<int>, 2> mNTrackletsPerClusterSum;

  // ---- Group B: plan-sized CA construction/result transients ----
  std::vector<o2::its::bounded_vector<o2::its::Tracklet>> mTracklets;
  std::vector<o2::its::bounded_vector<int>> mTrackletsLookupTable;
  std::vector<o2::its::bounded_vector<o2::MCCompLabel>> mTrackletLabels;
  o2::its::bounded_vector<float> mTransitionPhiCuts;
  o2::its::bounded_vector<float> mTransitionMSAngles;
  std::vector<o2::its::bounded_vector<CellSeed>> mCells;
  std::vector<o2::its::bounded_vector<int>> mCellsLookupTable;
  std::vector<o2::its::bounded_vector<int>> mCellsNeighbours;
  std::vector<o2::its::bounded_vector<int>> mCellsNeighboursTopology;
  std::vector<o2::its::bounded_vector<int>> mCellsNeighboursLUT;
  std::vector<o2::its::bounded_vector<o2::MCCompLabel>> mCellLabels;
  // The shared navigation auxiliary.
  IndexTableUtilsCore mIndexTableUtils;

  // ---- Group D: vertexer working scratch ----
  // Never plan-sized (bound by ROF count and the fixed pair count of 2).
  // Not resized by adoptPlan(): ROF count is not known at plan-adoption time.
  std::vector<o2::its::bounded_vector<o2::its::Line>> mLines;
  std::vector<o2::its::bounded_vector<o2::its::ClusterLines>> mTrackletClusters;
  std::vector<o2::its::bounded_vector<int>> mNTrackletsPerROF;
  std::array<o2::its::bounded_vector<int>, 2> mTrackletsIndexROF;
  std::vector<o2::its::bounded_vector<o2::MCCompLabel>> mLinesLabels;
  std::array<uint32_t, 2> mTotalTracklets{0, 0};
  uint32_t mTotalLines{0};
  unsigned int mNTotalLowPtVertices{0};
  void computeTrackletsPerROFScans();
  int& getNTrackletsROF(int rofId, int combId) { return mNTrackletsPerROF[combId][rofId]; }
  auto& getLines(int rofId) { return mLines[rofId]; }
  int getNLinesTotal() const noexcept { return mTotalLines; }
  void setNLinesTotal(uint32_t a) noexcept { mTotalLines = a; }
  auto& getTrackletClusters(int rofId) { return mTrackletClusters[rofId]; }
  gsl::span<const o2::its::Tracklet> getFoundTracklets(int rofId, int combId) const;
  gsl::span<o2::its::Tracklet> getFoundTracklets(int rofId, int combId);
  gsl::span<const MCCompLabel> getLabelsFoundTracklets(int rofId, int combId) const;
  gsl::span<int> getNTrackletsCluster(int rofId, int combId);
  gsl::span<int> getExclusiveNTrackletsCluster(int rofId, int combId);
  uint32_t getTotalTrackletsTF(const int iLayer) { return mTotalTracklets[iLayer]; }
  int getTotalClustersPerROFrange(int rofMin, int range, int layerId) const;

 private:
  void prepareClusters(const TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers,
                       gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements);

  RuntimeROFViews mROFViews{};
  bool mUseUPC{false};

  std::size_t mNOwnedSurfaces{0};
  std::size_t mNTransitions{0};
  std::size_t mNCells{0};
};

// M6d: the SurfaceTrackingScratch overload of resetTimeFrameEvent() --
// Tracker::clustersToTracks()'s recoverable-failure path calls this
// unqualified (ADL). Same
// reset-scratch-then-wipe-frame sequencing, same "not the future combined-
// owner policy" caveat as the original.
inline void resetTimeFrameEvent(TimeFrame& frame, SurfaceTrackingScratch& scratch) noexcept
{
  scratch.reset();
  frame.wipe();
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACETRACKINGSCRATCH_H_ */
