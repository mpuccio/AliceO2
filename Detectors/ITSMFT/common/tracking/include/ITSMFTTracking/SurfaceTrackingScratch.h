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
/// \brief M6c/M6d: the non-templated, detector-neutral successor to the
/// generic portions of the former fixed-layer scratch<NLayers>
/// (doc/design/0002-m6-generic-workspace-migration.md Sec 4, 9).
///
/// SurfaceTrackingScratch owns the M6a audit's Groups A, B, D, and E
/// (the former fixed-layer scratch's per-owned-surface cluster/index-table cache,
/// plan-sized CA construction/result transients, vertexer working scratch,
/// and memory/allocator/device plumbing) with every std::array<T, NLayers>
/// bound replaced by a runtime count supplied at plan-adoption time
/// (adoptPlan()): one slot per owned surface for Group A, and the sparse
/// transition/cell counts a SurfacePlanBinding already exposes for Group B.
/// Group D (vertexer scratch) was never actually NLayers-sized -- it carries
/// over unchanged in shape, still resized later from ROF count, not here.
///
/// M6d (production wiring for MFT, see the design note's own M6c-section
/// addendum) added the now-retired Group C output staging and a small set
/// of auxiliary NLayers-templated types (ROFOverlapTable, ROFVertexLookupTable,
/// ROFMaskTable, TrackingTopology) this scratch owns. IndexTableUtils is a
/// true exception: since M6e2 it is alias-erased to one shared, non-templated
/// IndexTableUtilsCore (IndexTableUtils.h), so a single member serves both
/// detectors and needs no duplication.
///
/// M6e2 (this scratch becomes shared by ITS(7) too, not MFT(10)-only) found
/// that ROFOverlapTable<N>/ROFVertexLookupTable<N>/ROFMaskTable<N>/
/// TrackingTopology<N> are genuinely N-sized (fixed-capacity internal
/// storage, e.g. TrackingTopology<N>::MaxTransitions = N*(N-1)/2), so a
/// single MFT(10)-shaped instance of any of them is wrong for ITS(7) --
/// not just a compile-time mismatch but a real correctness bug (wrong
/// transition/cell counts, wrong ROF-layer array width). This scratch
/// therefore stores these four groups exactly like Group C: two
/// separately-typed members (one ITS(7)-shaped, one MFT(10)-shaped, always
/// exactly one populated per instance -- each participant owns its own
/// instance, for exactly one detector), selected at compile time via a
/// template accessor `<NLayers>`. Every other member of this class remains
/// genuinely detector/layer-count-agnostic, sized only from adoptPlan()'s
/// runtime counts.
///
/// Like the former fixed-layer scratch, this type never owns or references a plan/
/// binding object -- adoptPlan() takes plain runtime counts (nOwnedSurfaces,
/// nTransitions, nCells), the exact three numbers a detector-neutral
/// surface-plan binding's own owned-surface count and global-transition/
/// -cell span sizes already expose, without this header including that
/// binding's own (detail/-confined) header at all. It does now reference
/// TimeFrame (initialise()/loadNormalizedSource()/updateROFVertexLookupTable()/
/// getPrimaryVertices() all take one, mirroring the former fixed-layer scratch's own
/// signatures exactly) -- M6c's original "never touches TimeFrame" framing
/// applied to that milestone's unwired scope only; M6d's whole point is
/// wiring this type into the same TimeFrame-cooperating role
/// the former fixed-layer scratch already has for ITS.
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
#include "ITSMFTTracking/MFTCATrack.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingTopology.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/ClusterLines.h"
#include "ITStracking/ExternalAllocator.h"
#include "ITStracking/ROFLookupTables.h"
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

/// Non-templated, detector-neutral CA working state, plus (M6d/M6e2) a set of
/// still-NLayers-templated auxiliary types this scratch owns dual-typed
/// (ITS(7)+MFT(10)). See the file-level doc above for exactly what this does
/// and does not own.
class SurfaceTrackingScratch
{
 private:
  // ---- Group E: memory/allocator/device plumbing ----
  // Declared first, exactly like the former fixed-layer scratch<NLayers>'s equivalent
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
  // were destroyed. Fixed by matching the former fixed-layer scratch<NLayers>'s
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
  // M6d: LoadTargetImpl-equivalent staging code (MultiSourceTimeFrameLoader.cxx)
  // needs direct access to mExternalAllocator/mExtMemoryPool to preserve
  // allocator identity across stage()/commit(), exactly mirroring
  // the former fixed-layer scratch<NLayers>'s own sole-friend contract.
  friend class MultiSourceTimeFrameLoader;

  static constexpr int MFTNLayers = o2::mft::constants::mft::LayersNumber;
  // IndexTableUtils is alias-erased (IndexTableUtils.h) to one shared,
  // non-templated IndexTableUtilsCore regardless of N -- a single member
  // (mIndexTableUtils below) already serves both detectors, no duplication
  // needed. CellSeed and TrackSeed are likewise not detector-typed and are
  // used directly by TrackerTraits/DetectorTraits, not stored as scratch
  // members.
  using IndexTableUtilsN = o2::itsmft::IndexTableUtils<MFTNLayers>;

  SurfaceTrackingScratch() = default;
  ~SurfaceTrackingScratch() = default;
  SurfaceTrackingScratch(const SurfaceTrackingScratch&) = delete;
  SurfaceTrackingScratch& operator=(const SurfaceTrackingScratch&) = delete;
  SurfaceTrackingScratch(SurfaceTrackingScratch&&) = delete;
  SurfaceTrackingScratch& operator=(SurfaceTrackingScratch&&) = delete;

  /// Sizes every Group A container to `nOwnedSurfaces` (one slot per owned
  /// surface, replacing every the former fixed-layer scratch<NLayers>
  /// std::array<T, NLayers>) and every Group B container to `nTransitions`/
  /// `nCells` (already-runtime sparse-topology counts). Precondition, exactly
  /// as the former fixed-layer scratch<NLayers>::initialise(): setMemoryPool() has
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
  /// the former fixed-layer scratch<NLayers>::resetScratch() member-for-member (same
  /// deep-clear-vs-framework-allocator-skip, same MC-info-conditional
  /// clearing, same non-owning mClusterLabels re-nulling instead of
  /// deep-clear). Never touches a TimeFrame, even implicitly, and never
  /// changes the adopted plan sizing (getNOwnedSurfaces()/getNTransitions()/
  /// getNCells() are unaffected -- only each container's *contents* are
  /// cleared, exactly as the original never shrinks its NLayers-wide outer
  /// arrays either). Matches resetScratch()'s own name too now (M6d wires
  /// this in where production code calls it).
  void reset();
  void resetScratch() { reset(); }

  /// memory management -- Group E. Doc mirrors
  /// the former fixed-layer scratch<NLayers>::setMemoryPool(): reseats every
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

  /// Staging/swap support for M6d's atomic loader migration -- mirrors the
  /// stage-then-commit discipline the former fixed-layer scratch<NLayers>::
  /// loadNormalizedSource() already implements for Group A, generalized to
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

  // ---- read-in data (M6d): mirrors
  // the former fixed-layer scratch<NLayers>::loadNormalizedSource() exactly, except
  // every NLayers-bound loop becomes a runtime orderedSurfaces.size() (==
  // getNOwnedSurfaces()) loop. detId preflight accepts both common-CA
  // detectors. See the former fixed-layer scratch.cxx for
  // the byte-for-byte original this was ported from.
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

  // navigation tables. getIndexTableUtils() is a single shared member (see
  // the type-alias comment above). Everything else here is dual-typed
  // (M6e2, see file doc) and selected at compile time via an explicit
  // `<NLayers>` template accessor -- mirrors Group C's getTracks<NLayers>()
  // pattern exactly, for the same reason (ROFOverlapTable<N>/
  // ROFVertexLookupTable<N>/ROFMaskTable<N>/TrackingTopology<N> are
  // genuinely N-sized production types, not detector-neutral). Every
  // Group-C-style call sites inside the shared Tracker<NLayers>/
  // TrackerTraits<NLayers> bodies go through the scratchXxx<NLayers>() free
  // dispatcher functions below instead of calling these accessors directly.
  const auto& getIndexTableUtils() const { return mIndexTableUtils; }

  template <int NLayers>
  static constexpr void checkSupportedNLayers() noexcept
  {
    static_assert(NLayers == ITSNLayers || NLayers == MFTNLayers, "SurfaceTrackingScratch's per-detector accessors support ITS (7) and MFT (10) only");
  }

  template <int NLayers>
  const auto& getROFOverlapTable() const noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      return mROFOverlapTableITS;
    } else {
      return mROFOverlapTableMFT;
    }
  }
  template <int NLayers>
  const auto& getROFOverlapTableView() const noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      return mROFOverlapTableViewITS;
    } else {
      return mROFOverlapTableViewMFT;
    }
  }
  template <int NLayers>
  const auto& getTrackerTopologies() const noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      return mTrackerTopologiesITS;
    } else {
      return mTrackerTopologiesMFT;
    }
  }
  template <int NLayers>
  const auto& getTrackingTopologyView() const noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      return mTrackingTopologyViewITS;
    } else {
      return mTrackingTopologyViewMFT;
    }
  }
  template <int NLayers>
  void setROFOverlapTable(o2::its::ROFOverlapTable<NLayers> table) noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      mROFOverlapTableITS = std::move(table);
      mROFOverlapTableViewITS = mROFOverlapTableITS.getView();
    } else {
      mROFOverlapTableMFT = std::move(table);
      mROFOverlapTableViewMFT = mROFOverlapTableMFT.getView();
    }
  }
  template <int NLayers>
  const auto& getROFVertexLookupTable() const noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      return mROFVertexLookupTableITS;
    } else {
      return mROFVertexLookupTableMFT;
    }
  }
  template <int NLayers>
  const auto& getROFVertexLookupTableView() const noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      return mROFVertexLookupTableViewITS;
    } else {
      return mROFVertexLookupTableViewMFT;
    }
  }
  template <int NLayers>
  void setROFVertexLookupTable(o2::its::ROFVertexLookupTable<NLayers> table) noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      mROFVertexLookupTableITS = std::move(table);
      mROFVertexLookupTableViewITS = mROFVertexLookupTableITS.getView();
    } else {
      mROFVertexLookupTableMFT = std::move(table);
      mROFVertexLookupTableViewMFT = mROFVertexLookupTableMFT.getView();
    }
  }
  template <int NLayers>
  gsl::span<const Vertex> getPrimaryVertices(const TimeFrame& frame, int layer, int rofId) const;
  template <int NLayers>
  void updateROFVertexLookupTable(const TimeFrame& frame) noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      mROFVertexLookupTableITS.update(frame.getPrimaryVertices().data(), frame.getPrimaryVertices().size());
    } else {
      mROFVertexLookupTableMFT.update(frame.getPrimaryVertices().data(), frame.getPrimaryVertices().size());
    }
  }
  template <int NLayers>
  void setMultiplicityCutMask(o2::its::ROFMaskTable<NLayers> cutMask) noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      mMultiplicityCutMaskITS = std::move(cutMask);
      mROFMaskViewITS = mROFMaskITS->getView();
    } else {
      mMultiplicityCutMaskMFT = std::move(cutMask);
      mROFMaskViewMFT = mROFMaskMFT->getView();
    }
  }
  template <int NLayers>
  void useMultiplictyMask() noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      mROFMaskITS = &mMultiplicityCutMaskITS;
      mROFMaskViewITS = mROFMaskITS->getView();
    } else {
      mROFMaskMFT = &mMultiplicityCutMaskMFT;
      mROFMaskViewMFT = mROFMaskMFT->getView();
    }
  }
  template <int NLayers>
  void setUPCCutMask(o2::its::ROFMaskTable<NLayers> cutMask) noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      mUPCCutMaskITS = std::move(cutMask);
    } else {
      mUPCCutMaskMFT = std::move(cutMask);
    }
  }
  template <int NLayers>
  void useUPCMask() noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      mROFMaskITS = &mUPCCutMaskITS;
      mROFMaskViewITS = mROFMaskITS->getView();
    } else {
      mROFMaskMFT = &mUPCCutMaskMFT;
      mROFMaskViewMFT = mROFMaskMFT->getView();
    }
  }
  template <int NLayers>
  const auto& getROFMaskView() const noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      return mROFMaskViewITS;
    } else {
      return mROFMaskViewMFT;
    }
  }

  const o2::its::TrackingFrameInfo& getClusterTrackingFrameInfo(int layerId, const o2::its::Cluster& cl) const;
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const o2::its::Cluster& cl) const { return getClusterLabels(layerId, cl.clusterId); }
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const int clId) const { return mClusterLabels[(mIsStaggered ? layerId : 0)]->getLabels(mClusterExternalIndices[layerId][clId]); }
  int getClusterExternalIndex(int layerId, const int clId) const { return mClusterExternalIndices[layerId][clId]; }
  int getClusterSize(int layer, int clusterId) const { return mClusterSize[layer][clusterId]; }
  void setClusterSize(int layer, o2::its::bounded_vector<uint8_t>& v) { mClusterSize[layer] = std::move(v); }

  auto& getTrackletsLabel(int layer) { return mTrackletLabels[layer]; }
  auto& getCellsLabel(int layer) { return mCellLabels[layer]; }

  /// M6d: TrackerTraits<NLayers>::initialiseTimeFrame()'s
  /// step 3-5 needs this typed the same way `mScratch->initialise()`'s
  /// existing call site already builds it -- see SurfaceTrackingScratch::IndexTableUtilsN.
  /// M6e2: reads/writes the dual ROF-mask/topology members above, so this is
  /// now a template on NLayers too; `layerMeasurements` takes a runtime-sized
  /// span-of-spans (not a fixed-width LayerMeasurementSpans<N> array) since
  /// nothing in its body needs a compile-time width, only a runtime layer
  /// index bounded by `maxLayers` -- both LayerMeasurementSpans<7> and <10>
  /// std::arrays convert to it implicitly via gsl::span's array constructor.
  template <int NLayers>
  void initialise(const TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers, int iteration,
                  const IndexTableUtilsN& indexTableConfig,
                  gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements);

  bool isClusterUsed(int layer, int clusterId) const { return mUsedClusters[layer][clusterId]; }
  void markUsedCluster(int layer, int clusterId) { mUsedClusters[layer][clusterId] = true; }
  gsl::span<unsigned char> getUsedClusters(const int layer);

  auto& getTracklets() { return mTracklets; }
  auto& getTrackletsLookupTable() { return mTrackletsLookupTable; }

  auto& getClusters() { return mClusters; }
  auto& getUnsortedClusters() { return mUnsortedClusters; }
  const auto& getUnsortedClusters() const { return mUnsortedClusters; }
  int getClusterROF(int iLayer, int iCluster);
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
  // the former fixed-layer scratch<NLayers>.
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
  // The one truly shared navigation auxiliary (see the type-alias comment
  // above).
  IndexTableUtilsN mIndexTableUtils;
  // M6d/M6e2: the NLayers-templated topology auxiliary (§4.1's deferred
  // item), dual-typed since M6e2 -- see file doc.
  std::vector<o2::itsmft::tracking::TrackingTopology<ITSNLayers>> mTrackerTopologiesITS;
  typename o2::itsmft::tracking::TrackingTopology<ITSNLayers>::View mTrackingTopologyViewITS;
  std::vector<o2::itsmft::tracking::TrackingTopology<MFTNLayers>> mTrackerTopologiesMFT;
  typename o2::itsmft::tracking::TrackingTopology<MFTNLayers>::View mTrackingTopologyViewMFT;
  // M6e1: default-constructed unless initDefaultTrackingTopology()/
  // initVertexingTopology() is called; see those methods' own doc for why
  // the former is no longer dead code (initVertexingTopology() still is).
  o2::itsmft::tracking::TrackingTopology<ITSNLayers> mDefaultTrackingTopologyITS;
  o2::itsmft::tracking::TrackingTopology<ITSNLayers> mVertexingTopologyITS;
  o2::itsmft::tracking::TrackingTopology<MFTNLayers> mDefaultTrackingTopologyMFT;
  o2::itsmft::tracking::TrackingTopology<MFTNLayers> mVertexingTopologyMFT;

  // ---- Group D: vertexer working scratch ----
  // Never actually NLayers-sized in the former fixed-layer scratch<NLayers> either
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
  template <int NLayers>
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

  /// initTrackerTopologies() (M6d): production caller is
  /// SurfacePlanTrackingParticipant<...>::configureRofTables() (combined MFT),
  /// unchanged body, generic via the fixed scratch type. Mirrors
  /// the former fixed-layer scratch<NLayers>::initTrackerTopologies() exactly.
  /// M6e2: writes the dual mTrackerTopologiesITS/MFT member, so this is now
  /// a template on NLayers too (see the file doc's Group-C-style rationale).
  template <int NLayers>
  void initTrackerTopologies(gsl::span<const TrackingParameters> trkParams, int maxLayers = NLayers);
  /// M6e1 correction: M6d's own claim that initVertexingTopology()/
  /// initDefaultTrackingTopology() have "zero production callers" was scoped
  /// only to the combined-participant-path files that milestone's own audit
  /// read (TrackerTraits.cxx/CATracker.cxx/SurfacePlanTrackingParticipant.cxx)
  /// -- never TrackingInterface.cxx, which was out of M6d's scope entirely.
  /// ITSMFTTrackingInterface<NLayers>::configureTrackingTopology() calls
  /// initDefaultTrackingTopology() unconditionally once per event, for both
  /// ITS and MFT -- so this scratch needs it too, now that it backs the
  /// standalone MFT interface. Mirrors
  /// the former fixed-layer scratch<NLayers>::initDefaultTrackingTopology() exactly.
  template <int NLayers>
  void initDefaultTrackingTopology(const TrackingParameters& trkParam, int maxLayers);
  /// initVertexingTopology() still has zero production callers even after
  /// M6e1 (grep-confirmed across TrackingInterface.cxx too) -- ported anyway
  /// for structural parity with the former fixed-layer scratch's own three-method
  /// topology-init group, since leaving only its sibling ported would be a
  /// more confusing asymmetry than one extra dead one-line mirror. Mirrors
  /// the former fixed-layer scratch<NLayers>::initVertexingTopology() exactly.
  template <int NLayers>
  void initVertexingTopology(const TrackingParameters& trkParam);

 private:
  template <int NLayers>
  void prepareClusters(const TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers,
                       gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements);

  // initialise() assigns mTrackingTopologyView{ITS,MFT} itself (not just
  // reads it), so it needs a mutable counterpart to the public const
  // getTrackingTopologyView<NLayers>() accessor above.
  template <int NLayers>
  auto& getTrackingTopologyViewMutable() noexcept
  {
    checkSupportedNLayers<NLayers>();
    if constexpr (NLayers == ITSNLayers) {
      return mTrackingTopologyViewITS;
    } else {
      return mTrackingTopologyViewMFT;
    }
  }

  o2::its::ROFOverlapTable<ITSNLayers> mROFOverlapTableITS;
  typename o2::its::ROFOverlapTable<ITSNLayers>::View mROFOverlapTableViewITS;
  o2::its::ROFOverlapTable<MFTNLayers> mROFOverlapTableMFT;
  typename o2::its::ROFOverlapTable<MFTNLayers>::View mROFOverlapTableViewMFT;
  o2::its::ROFVertexLookupTable<ITSNLayers> mROFVertexLookupTableITS;
  typename o2::its::ROFVertexLookupTable<ITSNLayers>::View mROFVertexLookupTableViewITS;
  o2::its::ROFVertexLookupTable<MFTNLayers> mROFVertexLookupTableMFT;
  typename o2::its::ROFVertexLookupTable<MFTNLayers>::View mROFVertexLookupTableViewMFT;
  o2::its::ROFMaskTable<ITSNLayers> mMultiplicityCutMaskITS;
  o2::its::ROFMaskTable<ITSNLayers> mUPCCutMaskITS;
  o2::its::ROFMaskTable<ITSNLayers>* mROFMaskITS = &mMultiplicityCutMaskITS;
  typename o2::its::ROFMaskTable<ITSNLayers>::View mROFMaskViewITS;
  o2::its::ROFMaskTable<MFTNLayers> mMultiplicityCutMaskMFT;
  o2::its::ROFMaskTable<MFTNLayers> mUPCCutMaskMFT;
  o2::its::ROFMaskTable<MFTNLayers>* mROFMaskMFT = &mMultiplicityCutMaskMFT;
  typename o2::its::ROFMaskTable<MFTNLayers>::View mROFMaskViewMFT;

  std::size_t mNOwnedSurfaces{0};
  std::size_t mNTransitions{0};
  std::size_t mNCells{0};
};

// M6d: the SurfaceTrackingScratch overload of resetTimeFrameEvent() --
// Tracker<NLayers>::clustersToTracks()'s recoverable-failure path calls this
// unqualified (ADL). Same
// reset-scratch-then-wipe-frame sequencing, same "not the future combined-
// owner policy" caveat as the original.
inline void resetTimeFrameEvent(TimeFrame& frame, SurfaceTrackingScratch& scratch) noexcept
{
  scratch.reset();
  frame.wipe();
}

template <int NLayers>
inline auto& scratchROFOverlapTable(SurfaceTrackingScratch& scratch) noexcept
{
  return scratch.template getROFOverlapTable<NLayers>();
}

template <int NLayers>
inline auto& scratchROFOverlapTableView(SurfaceTrackingScratch& scratch) noexcept
{
  return scratch.template getROFOverlapTableView<NLayers>();
}

template <int NLayers>
inline const auto& scratchROFOverlapTableView(const SurfaceTrackingScratch& scratch) noexcept
{
  return scratch.template getROFOverlapTableView<NLayers>();
}

template <int NLayers>
inline auto& scratchROFVertexLookupTableView(SurfaceTrackingScratch& scratch) noexcept
{
  return scratch.template getROFVertexLookupTableView<NLayers>();
}

template <int NLayers>
inline auto& scratchROFMaskView(SurfaceTrackingScratch& scratch) noexcept
{
  return scratch.template getROFMaskView<NLayers>();
}

template <int NLayers>
inline auto& scratchTrackerTopologies(SurfaceTrackingScratch& scratch) noexcept
{
  return scratch.template getTrackerTopologies<NLayers>();
}

template <int NLayers>
inline auto& scratchTrackingTopologyView(SurfaceTrackingScratch& scratch) noexcept
{
  return scratch.template getTrackingTopologyView<NLayers>();
}

template <int NLayers>
inline void scratchUseUPCMask(SurfaceTrackingScratch& scratch) noexcept
{
  scratch.template useUPCMask<NLayers>();
}

template <int NLayers>
inline gsl::span<const Vertex> scratchGetPrimaryVertices(SurfaceTrackingScratch& scratch, const TimeFrame& frame, int layer, int rofId)
{
  return scratch.template getPrimaryVertices<NLayers>(frame, layer, rofId);
}

template <int NLayers>
inline void scratchInitialise(SurfaceTrackingScratch& scratch, const TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers, int iteration,
                              const typename SurfaceTrackingScratch::IndexTableUtilsN& indexTableConfig,
                              gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements)
{
  scratch.template initialise<NLayers>(frame, trkParam, maxLayers, iteration, indexTableConfig, layerMeasurements);
}

template <int NLayers>
inline void scratchInitTrackerTopologies(SurfaceTrackingScratch& scratch, gsl::span<const TrackingParameters> trkParams, int maxLayers = NLayers)
{
  scratch.template initTrackerTopologies<NLayers>(trkParams, maxLayers);
}

template <int NLayers>
inline void scratchInitDefaultTrackingTopology(SurfaceTrackingScratch& scratch, const TrackingParameters& trkParam, int maxLayers)
{
  scratch.template initDefaultTrackingTopology<NLayers>(trkParam, maxLayers);
}

template <int NLayers>
inline void scratchInitVertexingTopology(SurfaceTrackingScratch& scratch, const TrackingParameters& trkParam)
{
  scratch.template initVertexingTopology<NLayers>(trkParam);
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACETRACKINGSCRATCH_H_ */
