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
/// \brief M6c: the non-templated, detector-neutral successor to the generic
/// portions of LegacyTrackerScratch<NLayers>
/// (doc/design/0002-m6-generic-workspace-migration.md Sec 4, 9).
///
/// SurfaceTrackingScratch owns exactly the M6a audit's Groups A, B, D, and E
/// (LegacyTrackerScratch's per-owned-surface cluster/index-table cache,
/// plan-sized CA construction/result transients, vertexer working scratch,
/// and memory/allocator/device plumbing) with every std::array<T, NLayers>
/// bound replaced by a runtime count supplied at plan-adoption time
/// (adoptPlan()): one slot per owned surface for Group A, and the sparse
/// transition/cell counts a SurfacePlanBinding already exposes for Group B.
/// Group D (vertexer scratch) was never actually NLayers-sized -- it carries
/// over unchanged in shape, still resized later from ROF count, not here.
///
/// Deliberately excluded, and left for a later milestone to decide (not
/// silently dropped -- see the M6c handoff report for the explicit reasoning):
/// - the legacy per-owned-surface index-table *binning configuration*
///   (IndexTableUtils<NLayers>) -- itself an NLayers-templated type that
///   would need its own generic redesign, out of this milestone's authorized
///   scope;
/// - the CA topology *view* objects (TrackingTopology<NLayers> and its
///   vertexing/default/per-iteration instances) -- Group B's sizing is taken
///   directly from the caller-supplied transition/cell counts instead, which
///   already supersedes needing an owned topology-view object for sizing;
/// - Group C (mTracks/mTracksLabel legacy result staging, detector-output
///   typed) -- adapter-private, out of scope per the design note;
/// - raw ROFs, ROF-timing/overlap/vertex-lookup tables, ROF masks, and every
///   per-collider-detector compatibility sidecar -- explicitly required to
///   stay outside this type.
///
/// Like LegacyTrackerScratch, this type never owns or references a TimeFrame
/// or any binding/plan object -- adoptPlan() takes plain runtime counts
/// (nOwnedSurfaces, nTransitions, nCells), the exact three numbers a
/// detector-neutral surface-plan binding's own owned-surface count and
/// global-transition/-cell span sizes already expose, without this header
/// including that binding's own (detail/-confined) header at all: this type
/// must carry no detector-identity, per-detector-layer-count, per-source-index,
/// output/data-processing-pipeline-layer, or detail/-confined
/// policy-key/state-family knowledge of any kind, and this keeps that true
/// structurally, not just by convention (see the M6c handoff report's
/// dependency-scan test for the
/// exact, checkable form of this claim).
///
/// M6c only: this type is additive and unused by production. No participant
/// wires it in; LegacyTrackerScratch stays production-live and unmodified.

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACETRACKINGSCRATCH_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACETRACKINGSCRATCH_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "ITSMFTTracking/Cell.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/ClusterLines.h"
#include "ITStracking/ExternalAllocator.h"
#include "ITStracking/Tracklet.h"
#include "DataFormatsITS/TrackITS.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

namespace o2::itsmft::tracking
{

/// Non-templated, detector-neutral CA working state. See the file-level doc
/// above for exactly what this does and does not own.
class SurfaceTrackingScratch
{
 private:
  // ---- Group E: memory/allocator/device plumbing ----
  // Declared first, exactly like LegacyTrackerScratch<NLayers>'s equivalent
  // members and for the identical reason (see that class's own doc comment):
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
  // were destroyed. Fixed by matching LegacyTrackerScratch<NLayers>'s
  // ordering exactly, not by asserting an ownership guarantee this type
  // cannot actually enforce.
  std::shared_ptr<o2::its::BoundedMemoryResource> mExtMemoryPool;
  std::shared_ptr<o2::its::BoundedMemoryResource> mMemoryPool;
  o2::its::ExternalAllocator* mExternalAllocator{nullptr};
  bool mIsStaggered{false};
  // Flagged in the M6 design note as possibly dead (no production
  // read/write site found by that audit); carried over structurally rather
  // than dropped, since this milestone does not re-run that verification.
  o2::its::bounded_vector<std::array<float, 2>> mPValphaX;

 public:
  SurfaceTrackingScratch() = default;
  ~SurfaceTrackingScratch() = default;
  SurfaceTrackingScratch(const SurfaceTrackingScratch&) = delete;
  SurfaceTrackingScratch& operator=(const SurfaceTrackingScratch&) = delete;
  SurfaceTrackingScratch(SurfaceTrackingScratch&&) = delete;
  SurfaceTrackingScratch& operator=(SurfaceTrackingScratch&&) = delete;

  /// Sizes every Group A container to `nOwnedSurfaces` (one slot per owned
  /// surface, replacing every LegacyTrackerScratch<NLayers>
  /// std::array<T, NLayers>) and every Group B container to `nTransitions`/
  /// `nCells` (already-runtime sparse-topology counts). Precondition, exactly
  /// as LegacyTrackerScratch<NLayers>::initialise(): setMemoryPool() has
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
  /// LegacyTrackerScratch<NLayers>::resetScratch() member-for-member (same
  /// deep-clear-vs-framework-allocator-skip, same MC-info-conditional
  /// clearing, same non-owning mClusterLabels re-nulling instead of
  /// deep-clear). Never touches a TimeFrame, even implicitly, and never
  /// changes the adopted plan sizing (getNOwnedSurfaces()/getNTransitions()/
  /// getNCells() are unaffected -- only each container's *contents* are
  /// cleared, exactly as the original never shrinks its NLayers-wide outer
  /// arrays either).
  void reset();

  /// memory management -- Group E. Doc mirrors
  /// LegacyTrackerScratch<NLayers>::setMemoryPool(): reseats every
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

  /// Staging/swap support for M6d's future atomic loader migration --
  /// mirrors the stage-then-commit discipline
  /// LegacyTrackerScratch<NLayers>::loadNormalizedSource() already
  /// implements for Group A, generalized to every plan-sized container here.
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

  // ---- Group A: legacy per-(owned-surface) cluster/index-table cache ----
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
  // Not per-owned-surface: fixed at the two tracklet combinations (layers
  // 0-1, 1-2), sized from cluster count at load time, exactly as in
  // LegacyTrackerScratch<NLayers>.
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

  // ---- Group D: vertexer working scratch ----
  // Never actually NLayers-sized in LegacyTrackerScratch<NLayers> either
  // (bound by ROF count and the fixed layer-pair index of 2); carries over
  // unchanged in shape. Not resized by adoptPlan() -- ROF count is not known
  // at plan-adoption time, exactly as it was never known at NLayers-scratch
  // construction time either.
  std::vector<o2::its::bounded_vector<o2::its::Line>> mLines;
  std::vector<o2::its::bounded_vector<o2::its::ClusterLines>> mTrackletClusters;
  std::vector<o2::its::bounded_vector<int>> mNTrackletsPerROF;
  std::array<o2::its::bounded_vector<int>, 2> mTrackletsIndexROF;
  std::vector<o2::its::bounded_vector<o2::MCCompLabel>> mLinesLabels;
  std::array<uint32_t, 2> mTotalTracklets{0, 0};
  uint32_t mTotalLines{0};
  unsigned int mNTotalLowPtVertices{0};

 private:
  std::size_t mNOwnedSurfaces{0};
  std::size_t mNTransitions{0};
  std::size_t mNCells{0};
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACETRACKINGSCRATCH_H_ */
