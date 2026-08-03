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
/// \file LegacyTrackerScratch.h
/// \brief Gate 4 B3.1: temporary, per-detector legacy CA scratch/result state.
///
/// LegacyTrackerScratch<NLayers> is everything the old, merged
/// TimeFrame<NLayers> owned that is NLayers-templated or otherwise
/// single-detector-scoped: legacy per-layer cluster/index-table/topology
/// arrays, CA hot-loop scratch (tracklets/cells/neighbours), the legacy
/// per-detector result containers (mTracks/mTracksLabel), and vertexer
/// working scratch. It holds no permanent event data of its own -- vertices,
/// beam state, Bz, CommonTrack/TrackClusterReference storage, and the
/// normalized per-SurfaceId measurements all live on the one non-templated
/// TimeFrame (ITSMFTTracking/TimeFrame.h) instead.
///
/// Lifetime contract: LegacyTrackerScratch<NLayers> never owns or stores a
/// reference to the TimeFrame it works alongside -- every method that needs
/// TimeFrame-owned data (beam position, the normalized frame, vertices)
/// takes it as an explicit parameter. This is deliberate: it means neither
/// type can "contain or reference" the other, so there is no ambiguity about
/// which one outlives which -- ownership is entirely in the hands of
/// whichever class constructs both (today, ITSMFTTrackingInterface<NLayers>;
/// see TrackingInterface.h, which declares its TimeFrame member before its
/// LegacyTrackerScratch<NLayers> member so C++'s reverse-declaration-order
/// destruction tears the scratch down first).
///
/// Reset: resetScratch() clears only scratch-owned state below and never
/// touches a TimeFrame, even implicitly -- see resetTimeFrameEvent() (this
/// header) for the one owner-level operation that sequences a scratch reset
/// with a TimeFrame wipe.
///
/// mDetId is not stored here either: LegacyTrackerScratch<NLayers>'s only
/// real runtime use of a detector identity (mapping a raw cluster to its
/// legacy layer during the pre-normalized loadROFrameData() path) is fully
/// determined by this class's own compile-time NLayers
/// (detIdFromNLayers<NLayers>()), so callers pass it explicitly rather than
/// it being independently stored/settable state.

#ifndef ALICEO2_ITSMFT_TRACKING_LEGACYTRACKERSCRATCH_H_
#define ALICEO2_ITSMFT_TRACKING_LEGACYTRACKERSCRATCH_H_

#include <array>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>
#include <gsl/gsl>

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITS/Vertex.h"

#include "ITStracking/BoundedAllocator.h"
#include "ITSMFTTracking/Cell.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/ClusterLines.h"
#include "ITStracking/ExternalAllocator.h"
#include "ITStracking/ROFLookupTables.h"
#include "ITStracking/Tracklet.h"

#include "ITSMFTTracking/MFTCATrack.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/LayerMask.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingTopology.h"
#ifndef GPUCA_GPUCODE
#include <optional>
#include "ITSMFTTracking/SurfaceCatalogView.h"
#endif
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include "DetectorsCommonDataFormats/DetID.h"

namespace o2
{
namespace itsmft
{
class CompClusterExt;
class TopologyDictionary;
class ROFRecord;
} // namespace itsmft

namespace itsmft::tracking
{

class MultiSourceTimeFrameLoader;

// Re-use ITS CA tracking data structures; only LegacyTrackerScratch and
// index-table I/O are detector-aware.
using Cluster = o2::its::Cluster;
using TrackingFrameInfo = o2::its::TrackingFrameInfo;
using Tracklet = o2::its::Tracklet;
using LayerMask = o2::itsmft::tracking::LayerMask;
using Line = o2::its::Line;
using ClusterLines = o2::its::ClusterLines;
using TrackITSExt = o2::its::TrackITSExt;
using ExternalAllocator = o2::its::ExternalAllocator;

// Non-owning per-layer normalized measurements used while building the
// transient sorted locator/navigation cache. LegacyTrackerScratch
// deliberately retains no copy: TrackerTraits owns the validated lifetime
// contract.
template <int NLayers>
using LayerMeasurementSpans = std::array<gsl::span<const SurfaceMeasurement>, NLayers>;

namespace constants
{
using namespace o2::its::constants;
using o2::itsmft::tracking::ITSNLayers;
using o2::itsmft::tracking::nLayersForDet;
} // namespace constants

template <int NLayers>
struct LegacyTrackerScratch {
  using IndexTableUtilsN = o2::itsmft::IndexTableUtils<NLayers>;
  using ROFOverlapTableN = o2::its::ROFOverlapTable<NLayers>;
  using ROFVertexLookupTableN = o2::its::ROFVertexLookupTable<NLayers>;
  using ROFMaskTableN = o2::its::ROFMaskTable<NLayers>;
  using TrackingTopologyN = o2::itsmft::tracking::TrackingTopology<NLayers>;
  using CellSeedN = o2::itsmft::tracking::CellSeedN<NLayers>;
  using TrackSeedN = o2::itsmft::tracking::TrackSeedN<NLayers>;

  LegacyTrackerScratch() = default;
  virtual ~LegacyTrackerScratch() = default;

  // The dormant combined-owner loader stages a second scratch with these
  // exact allocator owners, then exchanges only the normalized-load legacy
  // backfill in its final no-throw commit. It is deliberately the sole
  // friend: neither TimeFrame nor a tracker gains access to scratch state.
  friend class MultiSourceTimeFrameLoader;

 protected:
  // Both host memory-pool owners below must be declared -- and therefore
  // destroyed -- before every pmr/bounded_vector member that may allocate
  // through them (mClusters, mUnsortedClusters, mClusterExternalIndices,
  // mROFramesClusters, mClusterSize, and every other bounded_vector-typed
  // member declared later in this class). C++ destroys non-static data
  // members in reverse declaration order, so declaring these resource
  // owners first -- ahead of every other data member in this class --
  // guarantees they are destroyed last, after every vector that could
  // still hold memory allocated from them has already released it back.
  // Declaring either shared_ptr after any allocator-backed member would let
  // this scratch release the last external reference to that
  // BoundedMemoryResource while such a member still holds a buffer
  // allocated from it, leaving that member's later destructor call into
  // now-freed memory.
  std::shared_ptr<BoundedMemoryResource> mExtMemoryPool; // host memory pool managed by the framework
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  ExternalAllocator* mExternalAllocator{nullptr};

 public:
  // read-in data
  void loadROFrameData(gsl::span<const o2::itsmft::ROFRecord> rofs,
                       gsl::span<const itsmft::CompClusterExt> clusters,
                       gsl::span<const unsigned char>::iterator& pattIt,
                       const itsmft::TopologyDictionary* dict,
                       int layer,
                       const dataformats::MCTruthContainer<MCCompLabel>* mcLabels = nullptr,
                       o2::detectors::DetID::ID detId = o2::detectors::DetID::ITS);
  void resetROFrameData(int iLayer);
  void prepareROFrameData(gsl::span<const itsmft::CompClusterExt> clusters, int layer, o2::detectors::DetID::ID detId);

#ifndef GPUCA_GPUCODE
  // Owner-level current-bridge load operation (Gate 4 B3.1): stages `frame`'s
  // normalized-frame update and this scratch's own legacy backfill, then
  // commits both or neither. `frame` is the caller's shared, non-templated
  // TimeFrame -- the plan's one owner (ITSMFTTrackingInterface) supplies the
  // layer-to-surface mapping and the borrowed SurfaceCatalogView explicitly,
  // both derived once from its own immutable DetectorLayoutSet.
  //
  // Exactly one source is ever submitted here, and loadSources() requires
  // dense, zero-based source IDs, so the source ID is always
  // ClusterSourceId{0} -- not a caller-supplied parameter -- and every
  // ClusterRef produced by this call carries that source.
  //
  // `applySysErrors` defaults to true to match loadROFrameData()'s own
  // existing behavior above, which calls loadClusterTrackingFrameInfo<DetId>()
  // with its default applySysErrors=true.
  //
  // Preflight (in order, before any allocation/decoding/mutation) rejects:
  // NLayers not matching detId (UnsupportedDetector); an empty/unconfigured
  // `catalogView` (SurfaceCatalogNotConfigured); orderedSurfaces.size() !=
  // NLayers (InvalidLayerMapping); any mapped SurfaceId invalid, out of
  // range, or duplicated (InvalidLayerMapping); any mapped
  // SurfaceDescriptor.detectorId != detId (DetectorSurfaceMismatch). Every
  // preflight failure returns source ClusterSourceId{0}, consistent with
  // this single-source API.
  //
  // Preflight, decoding, and the complete legacy backfill are staged before
  // either representation is committed. A returned failing (non-ok()) result
  // or an exception during staging therefore leaves both `frame` (its
  // normalized owner and CommonTrack/TrackClusterReference storage) and this
  // scratch's every legacy compatibility structure unchanged. The final
  // commit uses only checked equal-allocator swaps, a no-throw
  // frame.commitNormalizedFrame() call, and no-throw owner moves.
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
  int getNrof(int layer) const { return mROFramesClusters[layer].size() - 1; }

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

  gsl::span<Cluster> getClustersOnLayer(int rofId, int layerId);
  gsl::span<const Cluster> getClustersOnLayer(int rofId, int layerId) const;
  gsl::span<const Cluster> getClustersPerROFrange(int rofMin, int range, int layerId) const;
  gsl::span<const Cluster> getUnsortedClustersOnLayer(int rofId, int layerId) const;
  gsl::span<uint8_t> getUsedClustersROF(int rofId, int layerId);
  gsl::span<const uint8_t> getUsedClustersROF(int rofId, int layerId) const;
  gsl::span<const int> getROFramesClustersPerROFrange(int rofMin, int range, int layerId) const;
  gsl::span<const int> getROFrameClusters(int layerId) const;
  gsl::span<const int> getNClustersROFrange(int rofMin, int range, int layerId) const;
  gsl::span<int> getIndexTable(int rofId, int layerId);
  const auto& getTrackingFrameInfoOnLayer(int layerId) const { return mTrackingFrameInfo[layerId]; }

  // navigation tables
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
  // `frame`'s primary vertices, restricted to the (layer, rofId) window this
  // scratch's own ROF-vertex-lookup table defines. Cross-cutting by
  // construction (the window is scratch-owned, the vertices themselves are
  // frame-owned), so `frame` is an explicit parameter rather than either
  // owner storing a reference to the other.
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

  const TrackingFrameInfo& getClusterTrackingFrameInfo(int layerId, const Cluster& cl) const;
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const Cluster& cl) const { return getClusterLabels(layerId, cl.clusterId); }
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const int clId) const { return mClusterLabels[((mIsStaggered) ? layerId : 0)]->getLabels(mClusterExternalIndices[layerId][clId]); }
  int getClusterExternalIndex(int layerId, const int clId) const { return mClusterExternalIndices[layerId][clId]; }
  int getClusterSize(int layer, int clusterId) const { return mClusterSize[layer][clusterId]; }
  void setClusterSize(int layer, bounded_vector<uint8_t>& v) { mClusterSize[layer] = std::move(v); }

  auto& getTrackletsLabel(int layer) { return mTrackletLabels[layer]; }
  auto& getCellsLabel(int layer) { return mCellLabels[layer]; }

  bool hasMCinformation() const { return mClusterLabels[0] != nullptr; }
  void initVertexingTopology(const TrackingParameters& trkParam);
  void initDefaultTrackingTopology(const TrackingParameters& trkParam, const int maxLayers = NLayers);
  void initTrackerTopologies(gsl::span<const TrackingParameters> trkParams, const int maxLayers = NLayers);
  // `indexTableConfig` must already be validated by the caller (see
  // ITSMFTTracking/IndexTableConfiguration.h::bindIndexTableConfiguration);
  // LegacyTrackerScratch never inspects a transition policy or DetID to
  // derive it, and never fails to configure/allocate on this account -- it
  // simply commits an already-valid value. `frame`'s beam position feeds
  // prepareClusters()'s XY binning below. The only call site
  // (TrackerTraits::initialiseTimeFrame) always supplies every argument
  // explicitly, so no parameter here carries a default.
  void initialise(const TimeFrame& frame, const TrackingParameters& trkParam, const int maxLayers, const int iteration,
                  const IndexTableUtilsN& indexTableConfig,
                  const LayerMeasurementSpans<NLayers>& layerMeasurements);

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
  auto& getTracks() { return mTracks; }
  const auto& getTracks() const { return mTracks; }
  auto& getTracksLabel() { return mTracksLabel; }
  const auto& getTracksLabel() const { return mTracksLabel; }
  auto& getLinesLabel(const int rofId) { return mLinesLabels[rofId]; }

  size_t getNumberOfClusters() const;
  virtual size_t getNumberOfCells() const;
  virtual size_t getNumberOfTracklets() const;
  virtual size_t getNumberOfNeighbours() const;
  size_t getNumberOfTracks() const;
  size_t getNumberOfUsedClusters() const;

  /// memory management (own copy of the caller's pool -- see setMemoryPool()
  /// doc below for why two independent shared_ptr owners of the same
  /// underlying resource is safe here, unlike wipe()/reset() coordination).
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  auto& getMemoryPool() const noexcept { return mMemoryPool; }
  bool checkMemory(unsigned long max) { return getArtefactsMemory() < max; }
  unsigned long getArtefactsMemory() const;
  void printArtefactsMemory() const;

  auto getFrameworkAllocator() { return mExternalAllocator; };
  void setFrameworkAllocator(ExternalAllocator* ext);
  bool hasFrameworkAllocator() const noexcept { return mExternalAllocator != nullptr; }
  std::pmr::memory_resource* getMaybeFrameworkHostResource(bool forceHost = false) { return (hasFrameworkAllocator() && !forceHost) ? mExtMemoryPool.get() : mMemoryPool.get(); }

  /// staggering
  void setIsStaggered(bool b) noexcept { mIsStaggered = b; }

  // Vertexer
  void computeTrackletsPerROFScans();
  void computeTracletsPerClusterScans();
  int& getNTrackletsROF(int rofId, int combId) { return mNTrackletsPerROF[combId][rofId]; }
  auto& getLines(int rofId) { return mLines[rofId]; }
  int getNLinesTotal() const noexcept { return mTotalLines; }
  void setNLinesTotal(uint32_t a) noexcept { mTotalLines = a; }
  auto& getTrackletClusters(int rofId) { return mTrackletClusters[rofId]; }
  gsl::span<const Tracklet> getFoundTracklets(int rofId, int combId) const;
  gsl::span<Tracklet> getFoundTracklets(int rofId, int combId);
  gsl::span<const MCCompLabel> getLabelsFoundTracklets(int rofId, int combId) const;
  gsl::span<int> getNTrackletsCluster(int rofId, int combId);
  gsl::span<int> getExclusiveNTrackletsCluster(int rofId, int combId);
  uint32_t getTotalTrackletsTF(const int iLayer) { return mTotalTracklets[iLayer]; }
  int getTotalClustersPerROFrange(int rofMin, int range, int layerId) const;
  // \Vertexer

  int hasBogusClusters() const { return std::accumulate(mBogusClusters.begin(), mBogusClusters.end(), 0); }

  template <typename... T>
  void addClusterToLayer(int layer, T&&... args);
  template <typename... T>
  void addTrackingFrameInfoToLayer(int layer, T&&... args);
  void addTrackingFrameInfoToLayer(int layer, const TrackingFrameInfo& tfInfo) { mTrackingFrameInfo[layer].push_back(tfInfo); }
  void addClusterExternalIndexToLayer(int layer, const int idx) { mClusterExternalIndices[layer].push_back(idx); }

  std::array<bounded_vector<Cluster>, NLayers> mClusters;
  std::array<bounded_vector<TrackingFrameInfo>, NLayers> mTrackingFrameInfo;
  std::array<bounded_vector<int>, NLayers> mClusterExternalIndices;
  std::array<bounded_vector<int>, NLayers> mROFramesClusters;
  std::array<const dataformats::MCTruthContainer<MCCompLabel>*, NLayers> mClusterLabels{nullptr};
  std::array<bounded_vector<int>, 2> mNTrackletsPerCluster;
  std::array<bounded_vector<int>, 2> mNTrackletsPerClusterSum;
  std::array<bounded_vector<int>, NLayers> mNClustersPerROF;
  std::array<bounded_vector<int>, NLayers> mIndexTables;
  std::vector<bounded_vector<int>> mTrackletsLookupTable;
  std::array<bounded_vector<uint8_t>, NLayers> mUsedClusters;

  std::array<bounded_vector<Cluster>, NLayers> mUnsortedClusters;
  std::vector<bounded_vector<Tracklet>> mTracklets;
  std::vector<bounded_vector<CellSeedN>> mCells;
  bounded_vector<CATrackType<NLayers>> mTracks;
  bounded_vector<MCCompLabel> mTracksLabel;
  std::vector<bounded_vector<int>> mCellsNeighbours;
  std::vector<bounded_vector<int>> mCellsNeighboursTopology;
  std::vector<bounded_vector<int>> mCellsLookupTable;

  virtual void resetScratch();

  // interface
  virtual bool isGPU() const noexcept { return false; }
  virtual const char* getName() const noexcept { return "CPU"; }

 protected:
  void prepareClusters(const TimeFrame& frame, const TrackingParameters& trkParam, const int maxLayers,
                       const LayerMeasurementSpans<NLayers>& layerMeasurements);
  // mExtMemoryPool/mMemoryPool/mExternalAllocator are declared at the very
  // top of this class (right after the constructor/destructor) instead of
  // here, so they are destroyed last, after every pmr/bounded_vector member
  // below -- see that declaration's own doc comment for why.
  unsigned int mNTotalLowPtVertices = 0;
  std::array<float, NLayers> mMinR;
  std::array<float, NLayers> mMaxR;
  bounded_vector<float> mTransitionPhiCuts;
  bounded_vector<float> mTransitionMSAngles;
  bounded_vector<float> mPositionResolution;
  std::array<bounded_vector<uint8_t>, NLayers> mClusterSize;

  bounded_vector<std::array<float, 2>> mPValphaX; /// PV x and alpha for track propagation
  std::vector<bounded_vector<MCCompLabel>> mTrackletLabels;
  std::vector<bounded_vector<MCCompLabel>> mCellLabels;
  std::vector<bounded_vector<int>> mCellsNeighboursLUT;
  bounded_vector<int> mBogusClusters; /// keep track of clusters with wild coordinates

  // Vertexer working scratch -- intermediate state consumed while finding
  // vertices, distinct from the *result* (TimeFrame::getPrimaryVertices()),
  // which is event data owned by the shared TimeFrame instead.
  std::vector<bounded_vector<int>> mNTrackletsPerROF;
  std::vector<bounded_vector<Line>> mLines;
  std::vector<bounded_vector<ClusterLines>> mTrackletClusters;
  std::array<bounded_vector<int>, 2> mTrackletsIndexROF;
  std::vector<bounded_vector<MCCompLabel>> mLinesLabels;
  std::array<uint32_t, 2> mTotalTracklets = {0, 0};
  uint32_t mTotalLines = 0;
  // \Vertexer

  // lookup tables
  IndexTableUtilsN mIndexTableUtils;
  ROFOverlapTableN mROFOverlapTable;
  ROFOverlapTableN::View mROFOverlapTableView;
  TrackingTopologyN mVertexingTopology;
  TrackingTopologyN mDefaultTrackingTopology;
  std::vector<TrackingTopologyN> mTrackerTopologies;
  typename TrackingTopologyN::View mTrackingTopologyView;
  ROFVertexLookupTableN mROFVertexLookupTable;
  ROFVertexLookupTableN::View mROFVertexLookupTableView;
  ROFMaskTableN mMultiplicityCutMask;
  ROFMaskTableN mUPCCutMask;
  ROFMaskTableN* mROFMask = &mMultiplicityCutMask;
  ROFMaskTableN::View mROFMaskView;

  bool mIsStaggered{false};
};

// Internal helper for the owner-level reset operation: clears `scratch` then
// wipes `frame`, in that order, exactly once each. Neither step is ever
// reordered or split across two separate caller-visible calls -- see
// ITSMFTTrackingInterface<NLayers>::resetEvent(), the current single-
// detector-bridge wrapper around this, and Tracker<NLayers>'s own
// recoverable-failure handling, which uses this same helper directly.
//
// This is NOT the future combined-owner policy: a future owner with several
// participating LegacyTrackerScratch instances (one per detector) sharing
// one TimeFrame must reset every participating scratch first and only then
// wipe the one shared TimeFrame once -- resetting one detector's scratch
// must never, by itself, wipe another detector's still-valid event data, and
// the shared TimeFrame must never be wiped while any participating scratch
// might still reference (and thus dangle against) the frame content that
// wipe would discard. Whether a given recoverable failure warrants a whole-
// event wipe at all, versus resetting only the failing detector's own
// scratch, is that future combined owner's decision to make -- never encoded
// here, and never encoded in Tracker<NLayers> itself (see CATracker.h/.cxx).
template <int NLayers>
inline void resetTimeFrameEvent(TimeFrame& frame, LegacyTrackerScratch<NLayers>& scratch) noexcept
{
  scratch.resetScratch();
  frame.wipe();
}

template <int NLayers>
inline gsl::span<const int> LegacyTrackerScratch<NLayers>::getROFrameClusters(int layerId) const
{
  return {&mROFramesClusters[layerId][0], static_cast<gsl::span<const int>::size_type>(mROFramesClusters[layerId].size())};
}

template <int NLayers>
inline gsl::span<Cluster> LegacyTrackerScratch<NLayers>::getClustersOnLayer(int rofId, int layerId)
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const Cluster> LegacyTrackerScratch<NLayers>::getClustersOnLayer(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<const Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<uint8_t> LegacyTrackerScratch<NLayers>::getUsedClustersROF(int rofId, int layerId)
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUsedClusters[layerId][startIdx], static_cast<gsl::span<uint8_t>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const uint8_t> LegacyTrackerScratch<NLayers>::getUsedClustersROF(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUsedClusters[layerId][startIdx], static_cast<gsl::span<const uint8_t>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const Cluster> LegacyTrackerScratch<NLayers>::getClustersPerROFrange(int rofMin, int range, int layerId) const
{
  if (rofMin < 0 || rofMin >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofMin]}; // First cluster of rofMin
  int endIdx{mROFramesClusters[layerId][o2::gpu::CAMath::Min(rofMin + range, getNrof(layerId))]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<Cluster>::size_type>(endIdx - startIdx)};
}

template <int NLayers>
inline gsl::span<const int> LegacyTrackerScratch<NLayers>::getROFramesClustersPerROFrange(int rofMin, int range, int layerId) const
{
  int chkdRange{o2::gpu::CAMath::Min(range, getNrof(layerId) - rofMin)};
  return {&mROFramesClusters[layerId][rofMin], static_cast<gsl::span<int>::size_type>(chkdRange)};
}

template <int NLayers>
inline gsl::span<const int> LegacyTrackerScratch<NLayers>::getNClustersROFrange(int rofMin, int range, int layerId) const
{
  int chkdRange{o2::gpu::CAMath::Min(range, getNrof(layerId) - rofMin)};
  return {&mNClustersPerROF[layerId][rofMin], static_cast<gsl::span<int>::size_type>(chkdRange)};
}

template <int NLayers>
inline int LegacyTrackerScratch<NLayers>::getTotalClustersPerROFrange(int rofMin, int range, int layerId) const
{
  int startIdx{rofMin}; // First cluster of rofMin
  int endIdx{o2::gpu::CAMath::Min(rofMin + range, getNrof(layerId))};
  return mROFramesClusters[layerId][endIdx] - mROFramesClusters[layerId][startIdx];
}

template <int NLayers>
inline int LegacyTrackerScratch<NLayers>::getClusterROF(int iLayer, int iCluster)
{
  return std::lower_bound(mROFramesClusters[iLayer].begin(), mROFramesClusters[iLayer].end(), iCluster + 1) - mROFramesClusters[iLayer].begin() - 1;
}

template <int NLayers>
inline gsl::span<const Cluster> LegacyTrackerScratch<NLayers>::getUnsortedClustersOnLayer(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUnsortedClusters[layerId][startIdx], static_cast<gsl::span<Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<int> LegacyTrackerScratch<NLayers>::getIndexTable(int rofId, int layer)
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const int tableSize = mIndexTableUtils.getNrowBins() * mIndexTableUtils.getNcolBins() + 1;
  return {&mIndexTables[layer][rofId * tableSize], static_cast<gsl::span<int>::size_type>(tableSize)};
}

template <int NLayers>
template <typename... T>
void LegacyTrackerScratch<NLayers>::addClusterToLayer(int layer, T&&... values)
{
  mUnsortedClusters[layer].emplace_back(std::forward<T>(values)...);
}

template <int NLayers>
template <typename... T>
void LegacyTrackerScratch<NLayers>::addTrackingFrameInfoToLayer(int layer, T&&... values)
{
  mTrackingFrameInfo[layer].emplace_back(std::forward<T>(values)...);
}

template <int NLayers>
inline gsl::span<uint8_t> LegacyTrackerScratch<NLayers>::getUsedClusters(const int layer)
{
  return {&mUsedClusters[layer][0], static_cast<gsl::span<uint8_t>::size_type>(mUsedClusters[layer].size())};
}

template <int NLayers>
inline gsl::span<int> LegacyTrackerScratch<NLayers>::getNTrackletsCluster(int rofId, int combId)
{
  if (rofId < 0 || rofId >= getNrof(1)) {
    return {};
  }
  auto startIdx{mROFramesClusters[1][rofId]};
  return {&mNTrackletsPerCluster[combId][startIdx], static_cast<gsl::span<int>::size_type>(mROFramesClusters[1][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<int> LegacyTrackerScratch<NLayers>::getExclusiveNTrackletsCluster(int rofId, int combId)
{
  if (rofId < 0 || rofId >= getNrof(1)) {
    return {};
  }
  auto clusStartIdx{mROFramesClusters[1][rofId]};

  return {&mNTrackletsPerClusterSum[combId][clusStartIdx], static_cast<gsl::span<int>::size_type>(mROFramesClusters[1][rofId + 1] - clusStartIdx)};
}

template <int NLayers>
inline gsl::span<Tracklet> LegacyTrackerScratch<NLayers>::getFoundTracklets(int rofId, int combId)
{
  if (rofId < 0 || rofId >= getNrof(1) || mTracklets[combId].empty()) {
    return {};
  }
  auto startIdx{mNTrackletsPerROF[combId][rofId]};
  return {&mTracklets[combId][startIdx], static_cast<gsl::span<Tracklet>::size_type>(mNTrackletsPerROF[combId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const Tracklet> LegacyTrackerScratch<NLayers>::getFoundTracklets(int rofId, int combId) const
{
  if (rofId < 0 || rofId >= getNrof(1)) {
    return {};
  }
  auto startIdx{mNTrackletsPerROF[combId][rofId]};
  return {&mTracklets[combId][startIdx], static_cast<gsl::span<Tracklet>::size_type>(mNTrackletsPerROF[combId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const MCCompLabel> LegacyTrackerScratch<NLayers>::getLabelsFoundTracklets(int rofId, int combId) const
{
  if (rofId < 0 || rofId >= getNrof(1) || !hasMCinformation()) {
    return {};
  }
  auto startIdx{mNTrackletsPerROF[combId][rofId]};
  return {&mTrackletLabels[combId][startIdx], static_cast<gsl::span<Tracklet>::size_type>(mNTrackletsPerROF[combId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline int LegacyTrackerScratch<NLayers>::getTotalClusters() const
{
  size_t totalClusters{0};
  for (const auto& clusters : mUnsortedClusters) {
    totalClusters += clusters.size();
  }
  return int(totalClusters);
}

template <int NLayers>
inline size_t LegacyTrackerScratch<NLayers>::getNumberOfClusters() const
{
  size_t nClusters{0};
  for (const auto& layer : mClusters) {
    nClusters += layer.size();
  }
  return nClusters;
}

template <int NLayers>
inline size_t LegacyTrackerScratch<NLayers>::getNumberOfCells() const
{
  size_t nCells{0};
  for (const auto& layer : mCells) {
    nCells += layer.size();
  }
  return nCells;
}

template <int NLayers>
inline size_t LegacyTrackerScratch<NLayers>::getNumberOfTracklets() const
{
  size_t nTracklets{0};
  for (const auto& layer : mTracklets) {
    nTracklets += layer.size();
  }
  return nTracklets;
}

template <int NLayers>
inline size_t LegacyTrackerScratch<NLayers>::getNumberOfNeighbours() const
{
  size_t neigh{0};
  for (const auto& l : mCellsNeighbours) {
    neigh += l.size();
  }
  return neigh;
}

template <int NLayers>
inline size_t LegacyTrackerScratch<NLayers>::getNumberOfTracks() const
{
  return mTracks.size();
}

template <int NLayers>
inline size_t LegacyTrackerScratch<NLayers>::getNumberOfUsedClusters() const
{
  size_t nClusters = 0;
  for (const auto& layer : mUsedClusters) {
    nClusters += std::count(layer.begin(), layer.end(), true);
  }
  return nClusters;
}

template <int NLayers>
inline const TrackingFrameInfo& LegacyTrackerScratch<NLayers>::getClusterTrackingFrameInfo(int layerId, const Cluster& cl) const
{
  return mTrackingFrameInfo[layerId][cl.clusterId];
}

template <int NLayers>
inline gsl::span<const Vertex> LegacyTrackerScratch<NLayers>::getPrimaryVertices(const TimeFrame& frame, int layer, int rofId) const
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const auto& entry = mROFVertexLookupTableView.getVertices(layer, rofId);
  const auto& vertices = frame.getPrimaryVertices();
  return {&vertices[entry.getFirstEntry()], static_cast<gsl::span<const Vertex>::size_type>(entry.getEntries())};
}

using LegacyTrackerScratchITS = LegacyTrackerScratch<ITSNLayers>;
/// MFT CA scratch: NLayers = half-disk CA layers (see MFTFwdTrackHelpers.h).
using LegacyTrackerScratchMFT = LegacyTrackerScratch<o2::mft::constants::mft::LayersNumber>;

} // namespace itsmft::tracking
} // namespace o2

#endif /* ALICEO2_ITSMFT_TRACKING_LEGACYTRACKERSCRATCH_H_ */
