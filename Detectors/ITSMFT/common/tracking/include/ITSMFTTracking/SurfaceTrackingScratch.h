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
/// generic portions of LegacyTrackerScratch<NLayers>
/// (doc/design/0002-m6-generic-workspace-migration.md Sec 4, 9).
///
/// SurfaceTrackingScratch owns the M6a audit's Groups A, B, D, and E
/// (LegacyTrackerScratch's per-owned-surface cluster/index-table cache,
/// plan-sized CA construction/result transients, vertexer working scratch,
/// and memory/allocator/device plumbing) with every std::array<T, NLayers>
/// bound replaced by a runtime count supplied at plan-adoption time
/// (adoptPlan()): one slot per owned surface for Group A, and the sparse
/// transition/cell counts a SurfacePlanBinding already exposes for Group B.
/// Group D (vertexer scratch) was never actually NLayers-sized -- it carries
/// over unchanged in shape, still resized later from ROF count, not here.
///
/// M6d (production wiring for MFT, see the design note's own M6c-section
/// addendum) adds Group C (mTracks/mTracksLabel -- still the sole production
/// detector-typed output staging path until M6e retires it) and a small set
/// of auxiliary NLayers-templated types (IndexTableUtils, ROFOverlapTable,
/// ROFVertexLookupTable, ROFMaskTable, TrackingTopology) this scratch owns
/// hardcoded to o2::mft::constants::mft::LayersNumber (=10). This is a
/// deliberate, narrow exception to "must not know 7/10 layers": these
/// specific auxiliary types are themselves still NLayers-templated
/// production types this milestone does not redesign (that redesign is its
/// own future scope, exactly like M6c's own IndexTableUtils/topology-view
/// deferral note already flagged), and SurfaceTrackingScratch is, for now,
/// only ever used by MFT (LegacyCATrackingParticipant<MFTNLayers,
/// SurfaceTrackingScratch, SurfacePlanBinding>) -- every other member of
/// this class remains genuinely detector/layer-count-agnostic, sized only
/// from adoptPlan()'s runtime counts.
///
/// Like LegacyTrackerScratch, this type never owns or references a plan/
/// binding object -- adoptPlan() takes plain runtime counts (nOwnedSurfaces,
/// nTransitions, nCells), the exact three numbers a detector-neutral
/// surface-plan binding's own owned-surface count and global-transition/
/// -cell span sizes already expose, without this header including that
/// binding's own (detail/-confined) header at all. It does now reference
/// TimeFrame (initialise()/loadNormalizedSource()/updateROFVertexLookupTable()/
/// getPrimaryVertices() all take one, mirroring LegacyTrackerScratch's own
/// signatures exactly) -- M6c's original "never touches TimeFrame" framing
/// applied to that milestone's unwired scope only; M6d's whole point is
/// wiring this type into the same TimeFrame-cooperating role
/// LegacyTrackerScratch already has for ITS.
#ifndef ALICEO2_ITSMFT_TRACKING_SURFACETRACKINGSCRATCH_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACETRACKINGSCRATCH_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
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

/// One slot per legacy-layer-equivalent (owned-surface) position -- mirrors
/// the original per-detector legacy scratch type's own
/// LayerMeasurementSpans<NLayers> shape, hardcoded to MFT's own NLayers (see
/// file doc for why this one type stays NLayers-templated at MFT's own
/// value).
using LayerMeasurementSpansMFT = std::array<gsl::span<const SurfaceMeasurement>, o2::mft::constants::mft::LayersNumber>;

/// Non-templated, detector-neutral CA working state, plus (M6d) MFT's own
/// hardcoded auxiliary types. See the file-level doc above for exactly what
/// this does and does not own.
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
  // M6d: LoadTargetImpl-equivalent staging code (MultiSourceTimeFrameLoader.cxx)
  // needs direct access to mExternalAllocator/mExtMemoryPool to preserve
  // allocator identity across stage()/commit(), exactly mirroring
  // LegacyTrackerScratch<NLayers>'s own sole-friend contract.
  friend class MultiSourceTimeFrameLoader;

  static constexpr int MFTNLayers = o2::mft::constants::mft::LayersNumber;
  using IndexTableUtilsN = o2::itsmft::IndexTableUtils<MFTNLayers>;
  using ROFOverlapTableN = o2::its::ROFOverlapTable<MFTNLayers>;
  using ROFVertexLookupTableN = o2::its::ROFVertexLookupTable<MFTNLayers>;
  using ROFMaskTableN = o2::its::ROFMaskTable<MFTNLayers>;
  using TrackingTopologyN = o2::itsmft::tracking::TrackingTopology<MFTNLayers>;
  using CellSeedN = CellSeed;
  using TrackSeedN = TrackSeedTpl<MFTNLayers>;

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
  /// arrays either). Matches resetScratch()'s own name too now (M6d wires
  /// this in where production code calls it).
  void reset();
  void resetScratch() { reset(); }

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

  /// Staging/swap support for M6d's atomic loader migration -- mirrors the
  /// stage-then-commit discipline LegacyTrackerScratch<NLayers>::
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
  // LegacyTrackerScratch<NLayers>::loadNormalizedSource() exactly, except
  // every NLayers-bound loop becomes a runtime orderedSurfaces.size() (==
  // getNOwnedSurfaces()) loop. detId preflight is narrowed to MFT only,
  // since this scratch type is MFT-only. See LegacyTrackerScratch.cxx for
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

  // navigation tables (M6d): mirror LegacyTrackerScratch's own accessor
  // names exactly (the compile-time seam TrackerTraits<NLayers, ScratchT,
  // BindingT> relies on) -- hardcoded to MFTNLayers, see file doc.
  const auto& getIndexTableUtils() const { return mIndexTableUtils; }
  const auto& getROFOverlapTable() const { return mROFOverlapTable; }
  const auto& getROFOverlapTableView() const { return mROFOverlapTableView; }
  const auto& getTrackerTopologies() const { return mTrackerTopologies; }
  const auto& getTrackingTopologyView() const { return mTrackingTopologyView; }
  void setROFOverlapTable(ROFOverlapTableN table)
  {
    mROFOverlapTable = std::move(table);
    mROFOverlapTableView = mROFOverlapTable.getView();
  }
  const auto& getROFVertexLookupTable() const { return mROFVertexLookupTable; }
  const auto& getROFVertexLookupTableView() const { return mROFVertexLookupTableView; }
  void setROFVertexLookupTable(ROFVertexLookupTableN table)
  {
    mROFVertexLookupTable = std::move(table);
    mROFVertexLookupTableView = mROFVertexLookupTable.getView();
  }
  gsl::span<const Vertex> getPrimaryVertices(const TimeFrame& frame, int layer, int rofId) const;
  void updateROFVertexLookupTable(const TimeFrame& frame) { mROFVertexLookupTable.update(frame.getPrimaryVertices().data(), frame.getPrimaryVertices().size()); }
  void setMultiplicityCutMask(ROFMaskTableN cutMask)
  {
    mMultiplicityCutMask = std::move(cutMask);
    mROFMaskView = mROFMask->getView();
  }
  void useMultiplictyMask() noexcept
  {
    mROFMask = &mMultiplicityCutMask;
    mROFMaskView = mROFMask->getView();
  }
  void setUPCCutMask(ROFMaskTableN cutMask) { mUPCCutMask = std::move(cutMask); }
  void useUPCMask() noexcept
  {
    mROFMask = &mUPCCutMask;
    mROFMaskView = mROFMask->getView();
  }
  const auto& getROFMaskView() const { return mROFMaskView; }

  const o2::its::TrackingFrameInfo& getClusterTrackingFrameInfo(int layerId, const o2::its::Cluster& cl) const;
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const o2::its::Cluster& cl) const { return getClusterLabels(layerId, cl.clusterId); }
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const int clId) const { return mClusterLabels[(mIsStaggered ? layerId : 0)]->getLabels(mClusterExternalIndices[layerId][clId]); }
  int getClusterExternalIndex(int layerId, const int clId) const { return mClusterExternalIndices[layerId][clId]; }
  int getClusterSize(int layer, int clusterId) const { return mClusterSize[layer][clusterId]; }
  void setClusterSize(int layer, o2::its::bounded_vector<uint8_t>& v) { mClusterSize[layer] = std::move(v); }

  auto& getTrackletsLabel(int layer) { return mTrackletLabels[layer]; }
  auto& getCellsLabel(int layer) { return mCellLabels[layer]; }

  /// M6d: TrackerTraits<NLayers, ScratchT, BindingT>::initialiseTimeFrame()'s
  /// step 3-5 needs this typed the same way `mScratch->initialise()`'s
  /// existing call site already builds it -- see ScratchN::IndexTableUtilsN.
  void initialise(const TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers, int iteration,
                  const IndexTableUtilsN& indexTableConfig,
                  const LayerMeasurementSpansMFT& layerMeasurements);

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
  // ---- Group C (M6d): legacy per-detector result staging. Still the sole
  // production source of detector-typed output today (M6a Group C
  // classification -- CommonTrack is already populated in parallel via
  // AcceptedTrackShadowPublisher, but the legacy-typed copy here is not
  // retired until M6e). CATrackType<MFTNLayers> resolves to MFTCATrack via
  // MFTCATrack.h's own CATrackTypeHelper<MFTNLayers> specialization. ----
  auto& getTracks() { return mTracks; }
  const auto& getTracks() const { return mTracks; }
  auto& getTracksLabel() { return mTracksLabel; }
  const auto& getTracksLabel() const { return mTracksLabel; }
  auto& getLinesLabel(const int rofId) { return mLinesLabels[rofId]; }

  size_t getNumberOfClusters() const;
  size_t getNumberOfCells() const;
  size_t getNumberOfTracklets() const;
  size_t getNumberOfNeighbours() const;
  size_t getNumberOfTracks() const { return mTracks.size(); }
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
  // M6d: the three NLayers-templated navigation/topology auxiliaries
  // (§4.1's deferred item), hardcoded to MFTNLayers -- see file doc.
  IndexTableUtilsN mIndexTableUtils;
  std::vector<TrackingTopologyN> mTrackerTopologies;
  typename TrackingTopologyN::View mTrackingTopologyView;
  // M6e1: default-constructed unless initDefaultTrackingTopology()/
  // initVertexingTopology() is called; see those methods' own doc for why
  // the former is no longer dead code (initVertexingTopology() still is).
  TrackingTopologyN mDefaultTrackingTopology;
  TrackingTopologyN mVertexingTopology;

  // ---- Group C (M6d): legacy per-detector result staging ----
  o2::its::bounded_vector<CATrackType<MFTNLayers>> mTracks;
  o2::its::bounded_vector<o2::MCCompLabel> mTracksLabel;

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
  /// LegacyCATrackingParticipant<...>::configureRofTables() (combined MFT),
  /// unchanged body, generic via ScratchN. Mirrors
  /// LegacyTrackerScratch<NLayers>::initTrackerTopologies() exactly.
  void initTrackerTopologies(gsl::span<const TrackingParameters> trkParams, int maxLayers = MFTNLayers);
  /// M6e1 correction: M6d's own claim that initVertexingTopology()/
  /// initDefaultTrackingTopology() have "zero production callers" was scoped
  /// only to the combined-participant-path files that milestone's own audit
  /// read (TrackerTraits.cxx/CATracker.cxx/LegacyCATrackingParticipant.cxx)
  /// -- never TrackingInterface.cxx, which was out of M6d's scope entirely.
  /// ITSMFTTrackingInterface<NLayers>::configureTrackingTopology() (the
  /// standalone-path owner M6e1 migrates) calls
  /// initDefaultTrackingTopology() unconditionally once per event, for both
  /// ITS and MFT -- so this scratch type needs it too, now that it backs the
  /// standalone MFT interface. Mirrors
  /// LegacyTrackerScratch<NLayers>::initDefaultTrackingTopology() exactly.
  void initDefaultTrackingTopology(const TrackingParameters& trkParam, int maxLayers);
  /// initVertexingTopology() still has zero production callers even after
  /// M6e1 (grep-confirmed across TrackingInterface.cxx too) -- ported anyway
  /// for structural parity with LegacyTrackerScratch's own three-method
  /// topology-init group, since leaving only its sibling ported would be a
  /// more confusing asymmetry than one extra dead one-line mirror. Mirrors
  /// LegacyTrackerScratch<NLayers>::initVertexingTopology() exactly.
  void initVertexingTopology(const TrackingParameters& trkParam);

 private:
  void prepareClusters(const TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers,
                       const LayerMeasurementSpansMFT& layerMeasurements);

  ROFOverlapTableN mROFOverlapTable;
  typename ROFOverlapTableN::View mROFOverlapTableView;
  ROFVertexLookupTableN mROFVertexLookupTable;
  typename ROFVertexLookupTableN::View mROFVertexLookupTableView;
  ROFMaskTableN mMultiplicityCutMask;
  ROFMaskTableN mUPCCutMask;
  ROFMaskTableN* mROFMask = &mMultiplicityCutMask;
  typename ROFMaskTableN::View mROFMaskView;

  std::size_t mNOwnedSurfaces{0};
  std::size_t mNTransitions{0};
  std::size_t mNCells{0};
};

// M6d: the SurfaceTrackingScratch overload of the original per-detector
// legacy scratch type's own resetTimeFrameEvent() free function --
// Tracker<NLayers, ScratchT, BindingT>::clustersToTracks()'s recoverable-
// failure path calls this unqualified (ADL), so an overload matching
// whichever ScratchT is in play must exist in this namespace. Same
// reset-scratch-then-wipe-frame sequencing, same "not the future combined-
// owner policy" caveat as the original.
inline void resetTimeFrameEvent(TimeFrame& frame, SurfaceTrackingScratch& scratch) noexcept
{
  scratch.reset();
  frame.wipe();
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACETRACKINGSCRATCH_H_ */
