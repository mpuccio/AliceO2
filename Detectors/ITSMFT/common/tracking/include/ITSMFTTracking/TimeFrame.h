// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
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

#ifndef ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_
#define ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_

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
#include "ITSMFTTracking/CommonTrack.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/LayerMask.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceMeasurementIndex.h"
#include "ITSMFTTracking/TrackingTopology.h"
#ifndef GPUCA_GPUCODE
#include <optional>
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"
#endif
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsCommonDataFormats/DetID.h"

namespace o2
{
namespace gpu
{
class GPUChainITS;
}

namespace itsmft
{
class CompClusterExt;
class TopologyDictionary;
class ROFRecord;
} // namespace itsmft

namespace itsmft::tracking
{

// Re-use ITS CA tracking data structures; only TimeFrame and index-table I/O are detector-aware.
using Cluster = o2::its::Cluster;
using TrackingFrameInfo = o2::its::TrackingFrameInfo;
using Tracklet = o2::its::Tracklet;
using LayerMask = o2::itsmft::tracking::LayerMask;
using BoundedMemoryResource = o2::its::BoundedMemoryResource;
template <typename T>
using bounded_vector = o2::its::bounded_vector<T>;
using ExternalAllocator = o2::its::ExternalAllocator;
using Line = o2::its::Line;
using ClusterLines = o2::its::ClusterLines;
using Vertex = o2::its::Vertex;
using VertexLabel = o2::its::VertexLabel;
using TrackITSExt = o2::its::TrackITSExt;

// Non-owning per-layer normalized measurements used while building the
// transient sorted locator/navigation cache. TimeFrame deliberately retains
// no copy: TrackerTraits owns the validated lifetime contract.
template <int NLayers>
using LayerMeasurementSpans = std::array<gsl::span<const SurfaceMeasurement>, NLayers>;

namespace constants
{
using namespace o2::its::constants;
using o2::itsmft::tracking::ITSNLayers;
using o2::itsmft::tracking::nLayersForDet;
} // namespace constants

template <int NLayers>
struct TimeFrame {
  using IndexTableUtilsN = o2::itsmft::IndexTableUtils<NLayers>;
  using ROFOverlapTableN = o2::its::ROFOverlapTable<NLayers>;
  using ROFVertexLookupTableN = o2::its::ROFVertexLookupTable<NLayers>;
  using ROFMaskTableN = o2::its::ROFMaskTable<NLayers>;
  using TrackingTopologyN = o2::itsmft::tracking::TrackingTopology<NLayers>;
  using CellSeedN = o2::itsmft::tracking::CellSeedN<NLayers>;
  using TrackSeedN = o2::itsmft::tracking::TrackSeedN<NLayers>;

  TimeFrame() = default;
  virtual ~TimeFrame() = default;

  const Vertex& getPrimaryVertex(const int ivtx) const { return mPrimaryVertices[ivtx]; }
  auto& getPrimaryVertices() { return mPrimaryVertices; };
  auto getPrimaryVerticesNum() { return mPrimaryVertices.size(); };
  const auto& getPrimaryVertices() const { return mPrimaryVertices; };
  auto& getPrimaryVerticesLabels() { return mPrimaryVerticesLabels; };
  gsl::span<const Vertex> getPrimaryVertices(int layer, int rofId) const;
  void addPrimaryVertex(const Vertex& vertex);
  void addPrimaryVertexLabel(const VertexLabel& label) { mPrimaryVerticesLabels.push_back(label); }

  // read-in data
  void loadROFrameData(gsl::span<const o2::itsmft::ROFRecord> rofs,
                       gsl::span<const itsmft::CompClusterExt> clusters,
                       gsl::span<const unsigned char>::iterator& pattIt,
                       const itsmft::TopologyDictionary* dict,
                       int layer,
                       const dataformats::MCTruthContainer<MCCompLabel>* mcLabels = nullptr,
                       o2::detectors::DetID::ID detId = o2::detectors::DetID::ITS);
  o2::detectors::DetID::ID getDetId() const { return mDetId; }
  void setDetId(o2::detectors::DetID::ID detId) { mDetId = detId; }
  void resetROFrameData(int iLayer);
  void prepareROFrameData(gsl::span<const itsmft::CompClusterExt> clusters, int layer);

  // Non-owning, read-only access to the normalized owner/view associated
  // with this TimeFrame by the most recent successful loadNormalizedSource()
  // call. Empty/default until that first succeeds, and after wipe() (see
  // below): wipe() unconditionally clears this owner in place. The
  // `const MultiSourceFrame&` returned by getNormalizedFrame() is a
  // reference to that same long-lived member object -- it remains valid
  // and safe to dereference across a wipe() call, it simply then observes
  // the owner's newly cleared (empty) state. What wipe() does invalidate is
  // any MultiSourceFrameView or gsl::span (getSurfaceMeasurements(),
  // getSourceIntervals(), getLabels(), getView()) obtained *before* the
  // wipe() call: those hold pointers into the owner's internal buffers,
  // which clear() may reallocate/free, so they must be re-obtained
  // afterwards rather than reused.
  const MultiSourceFrame& getNormalizedFrame() const noexcept { return mNormalizedFrame; }
  MultiSourceFrameView getNormalizedFrameView() const noexcept { return mNormalizedFrame.getView(); }

#ifndef GPUCA_GPUCODE
  // Gate 2 compatibility boundary: loads one single-detector cluster stream
  // -- the same raw ROFRecord/CompClusterExt/pattern/dictionary/label inputs
  // loadROFrameData() above already consumes -- through the normalized
  // MultiSourceFrame owner (loadSources(), single decode per cluster via the
  // caller-supplied `decoder`), then backfills this TimeFrame's existing
  // legacy compatibility structures (unsorted clusters, TrackingFrameInfo,
  // external indices, cluster sizes, ROF boundaries, label lookup) from the
  // committed normalized measurements. Purely additive: loadROFrameData()
  // and its production callers are unchanged, and this entry point does not
  // decode any compact cluster a second time.
  //
  // Unlike the Gate 1 signature, this no longer accepts an externally
  // supplied DetectorLayoutView or layer-to-surface mapping: both are
  // derived, from a single current DetectorLayoutSet snapshot obtained once,
  // via ensureDetectorLayouts()/getDetectorLayouts() below -- the canonical
  // SurfaceCatalogView from that set's owned catalog, and the layer-to-
  // surface mapping from its configurationKey.orderedSurfaces. No tracking
  // iteration (DetectorLayout) is selected or required, so a canonical
  // catalog configured with zero tracking iterations loads normally. Because
  // this now touches host-only layout ownership, the declaration lives
  // inside this GPUCA_GPUCODE guard.
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
  // NLayers not matching detId (UnsupportedDetector); no catalog owner ever
  // stored (SurfaceCatalogNotConfigured); a stored owner that is not current
  // (SurfaceCatalogStale); configurationKey.catalogRequest.detector != detId
  // (DetectorSurfaceMismatch); orderedSurfaces.size() != NLayers
  // (InvalidLayerMapping); any mapped SurfaceId invalid, out of range, or
  // duplicated (InvalidLayerMapping); any mapped SurfaceDescriptor.detectorId
  // != detId (DetectorSurfaceMismatch). Every preflight failure returns
  // source ClusterSourceId{0}, consistent with this single-source API.
  //
  // Preflight, decoding, and the complete legacy backfill are staged before
  // either representation is committed. A returned failing (non-ok()) result
  // or an exception during staging therefore leaves both the normalized owner
  // and every legacy compatibility accessor unchanged. The final commit uses
  // only checked equal-allocator swaps and no-throw owner moves.
  LoadSourcesResult loadNormalizedSource(const ClusterDecoder& decoder,
                                         const o2::InteractionRecord& origin,
                                         const ROFTimingConfig& timing,
                                         gsl::span<const itsmft::CompClusterExt> clusters,
                                         gsl::span<const unsigned char> patterns,
                                         gsl::span<const o2::itsmft::ROFRecord> rofs,
                                         const itsmft::TopologyDictionary* dictionary,
                                         const dataformats::MCTruthContainer<MCCompLabel>* labels,
                                         o2::detectors::DetID::ID detId,
                                         bool applySysErrors = true);

  // Builds (or reuses, if already current) this TimeFrame's DetectorLayoutSet
  // from `provider`'s resolved surface catalog. Per-surface nominal material
  // is not a parameter here: it lives directly on each SurfaceDescriptor in
  // the resolved catalog (SurfaceDescriptor::material) and is therefore
  // supplied by `provider`, not threaded through this call.
  DetectorLayoutSetBuildResult ensureDetectorLayouts(const DetectorSurfaceCatalogProvider* provider,
                                                     const DetectorSurfaceCatalogRequest& catalogRequest,
                                                     gsl::span<const SurfaceId> orderedSurfaces,
                                                     TransitionPolicyTag policyTag,
                                                     gsl::span<const TrackingParameters> trackingParameters);

  // Invalidates the current DetectorLayoutSet, forcing a rebuild on the next
  // ensureDetectorLayouts() call. A change to nominal surface material is a
  // change to the surface description (SurfaceDescriptor::material), so it
  // is invalidated through this same geometry path -- there is no separate
  // material invalidation entry point.
  void invalidateDetectorLayouts() noexcept;
  DetectorGeometryEpoch getRequiredDetectorGeometryEpoch() const noexcept { return mRequiredDetectorGeometryEpoch; }
  bool detectorLayoutsCurrent() const noexcept;
  bool hasStoredDetectorLayouts() const noexcept { return mDetectorLayouts.has_value(); }
  const std::vector<SurfaceDescriptor>* getSurfaceCatalog() const noexcept
  {
    const auto* layouts = getDetectorLayouts();
    return layouts ? &layouts->getSurfaceCatalog() : nullptr;
  }
  gsl::span<const SurfaceDescriptor> getSurfaceCatalogView() const noexcept
  {
    const auto* catalog = getSurfaceCatalog();
    return catalog ? gsl::span<const SurfaceDescriptor>{catalog->data(), catalog->size()} : gsl::span<const SurfaceDescriptor>{};
  }
  const DetectorLayoutSet* getDetectorLayouts() const noexcept
  {
    return detectorLayoutsCurrent() ? &*mDetectorLayouts : nullptr;
  }
  const DetectorLayout* getDetectorLayout(size_t iteration) const noexcept
  {
    const auto* layouts = getDetectorLayouts();
    return layouts ? layouts->getLayout(iteration) : nullptr;
  }
  DetectorLayoutView getDetectorLayoutView(size_t iteration) const noexcept
  {
    const auto* layouts = getDetectorLayouts();
    return layouts ? layouts->getLayoutView(iteration) : DetectorLayoutView{};
  }
#endif

  int getTotalClusters() const;
  bool empty() const { return getTotalClusters() == 0; }
  int getSortedIndex(int rofId, int layer, int idx) const { return mROFramesClusters[layer][rofId] + idx; }
  int getSortedStartIndex(const int rofId, const int layer) const { return mROFramesClusters[layer][rofId]; }
  int getNrof(int layer) const { return mROFramesClusters[layer].size() - 1; }

  void resetBeamXY(const float x, const float y, const float w = 0);
  void setBeamPosition(const float x, const float y, const float s2, const float base = 50.f, const float systematic = 0.f)
  {
    isBeamPositionOverridden = true;
    resetBeamXY(x, y, s2 / o2::gpu::CAMath::Sqrt((base * base) + systematic));
  }

  float getBeamX() const { return mBeamPos[0]; }
  float getBeamY() const { return mBeamPos[1]; }
  std::array<float, 2>& getBeamXY() { return mBeamPos; }

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
  void updateROFVertexLookupTable() { mROFVertexLookupTable.update(mPrimaryVertices.data(), mPrimaryVertices.size()); }
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
  // TimeFrame never inspects a TransitionPolicyTag or DetID to derive it, and
  // never fails to configure/allocate on this account -- it simply commits
  // an already-valid value. The only call site (TrackerTraits::
  // initialiseTimeFrame) always supplies every argument explicitly, so no
  // parameter here carries a default.
  void initialise(const TrackingParameters& trkParam, const int maxLayers, const int iteration,
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

  // Detector-neutral common-CA result storage (Gate 4 CommonTrack
  // foundation; ITSMFTTracking/CommonTrack.h). CommonTrack itself carries no
  // NLayers dependency; this is a temporary bridge holding it inside
  // TimeFrame<NLayers> alongside the legacy per-NLayers artefacts above.
  // Unpopulated by this slice -- no production or test call site writes
  // through these accessors from CA seeds yet -- but already event data:
  // both containers are wiped together with every other per-event CA
  // artefact by wipe() (see wipe()'s own doc) and are only meaningful
  // together with the TimeFrame's normalized frame (getNormalizedFrame())
  // that was current when a CommonTrack was built (see CommonTrack.h's own
  // lifetime doc).
  auto& getCommonTracks() { return mCommonTracks; }
  const auto& getCommonTracks() const { return mCommonTracks; }
  // Flat, TimeFrame-owned array of SurfaceMeasurementIndex; a CommonTrack's
  // [firstClusterRef, clusterRefEnd) range (CommonTrack.h) is a half-open
  // range of *positions* into this array, in traversal order (inner to
  // outer). Each element is, in turn, a canonical position into the
  // flattened SurfaceMeasurement array owned by this TimeFrame's normalized
  // frame (getNormalizedFrame()/getNormalizedFrameView()) -- resolved via
  // MultiSourceFrame::getMeasurement()/MultiSourceFrameView::getMeasurement().
  auto& getTrackClusterIndices() { return mTrackClusterIndices; }
  const auto& getTrackClusterIndices() const { return mTrackClusterIndices; }

  size_t getNumberOfClusters() const;
  virtual size_t getNumberOfCells() const;
  virtual size_t getNumberOfTracklets() const;
  virtual size_t getNumberOfNeighbours() const;
  size_t getNumberOfTracks() const;
  size_t getNumberOfUsedClusters() const;

  /// memory management
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  auto& getMemoryPool() const noexcept { return mMemoryPool; }
  bool checkMemory(unsigned long max) { return getArtefactsMemory() < max; }
  unsigned long getArtefactsMemory() const;
  void printArtefactsMemory() const;

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

  void setBz(float bz) { mBz = bz; }
  float getBz() const { return mBz; }

  /// State if memory will be externally managed by the GPU framework
  ExternalAllocator* mExternalAllocator{nullptr};
  // Both host memory-pool owners below must be declared -- and therefore
  // destroyed -- before every pmr/bounded_vector member that may allocate
  // through them (mClusters, mUnsortedClusters, mClusterExternalIndices,
  // mROFramesClusters, mClusterSize, and every other bounded_vector-typed
  // member declared later in this class, public or protected). C++ destroys
  // non-static data members in reverse declaration order, so declaring
  // these resource owners first guarantees they are destroyed last, after
  // every vector that could still hold memory allocated from them has
  // already released it back. Declaring either shared_ptr after any
  // allocator-backed member would let TimeFrame release the last external
  // reference to that BoundedMemoryResource while such a member still holds
  // a buffer allocated from it, leaving that member's later destructor call
  // into now-freed memory.
  std::shared_ptr<BoundedMemoryResource> mExtMemoryPool; // host memory pool managed by the framework
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  auto getFrameworkAllocator() { return mExternalAllocator; };
  void setFrameworkAllocator(ExternalAllocator* ext);
  bool hasFrameworkAllocator() const noexcept { return mExternalAllocator != nullptr; }
  std::pmr::memory_resource* getMaybeFrameworkHostResource(bool forceHost = false) { return (hasFrameworkAllocator() && !forceHost) ? mExtMemoryPool.get() : mMemoryPool.get(); }

  // Propagator
  const o2::base::PropagatorImpl<float>* getDevicePropagator() const { return mPropagatorDevice; }
  virtual void setDevicePropagator(const o2::base::PropagatorImpl<float>* /*unused*/) {};

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

  // Gate 4 CommonTrack foundation (ITSMFTTracking/CommonTrack.h): detector-
  // neutral common-CA result storage. Neither element type depends on
  // NLayers; see getCommonTracks()/getTrackClusterIndices() above for the
  // ownership/lifetime contract.
  bounded_vector<CommonTrack> mCommonTracks;
  bounded_vector<SurfaceMeasurementIndex> mTrackClusterIndices;

  const o2::base::PropagatorImpl<float>* mPropagatorDevice = nullptr; // Needed only for GPU
  o2::detectors::DetID::ID mDetId{o2::detectors::DetID::ITS};

  virtual void wipe();

  // interface
  virtual bool isGPU() const noexcept { return false; }
  virtual const char* getName() const noexcept { return "CPU"; }

 protected:
  void prepareClusters(const TrackingParameters& trkParam, const int maxLayers,
                       const LayerMeasurementSpans<NLayers>& layerMeasurements);
  float mBz = 5.;
  unsigned int mNTotalLowPtVertices = 0;
  int mBeamPosWeight = 0;
  std::array<float, 2> mBeamPos = {0.f, 0.f};
  bool isBeamPositionOverridden = false;
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

  // Vertexer
  bounded_vector<Vertex> mPrimaryVertices;
  bounded_vector<VertexLabel> mPrimaryVerticesLabels;
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

  // Normalized owner associated by loadNormalizedSource(); host-only, never
  // GPU-managed or dictionary-serialized (see getNormalizedFrame()). Does
  // not itself hold pmr/bounded-vector allocations (MultiSourceFrame's own
  // members are plain std::vector<T> with the default allocator), so it has
  // no ordering dependency on mMemoryPool/mExtMemoryPool above.
  MultiSourceFrame mNormalizedFrame;

#ifndef GPUCA_GPUCODE
  // Semantic configuration, unlike event artefacts, survives wipe().
  std::optional<DetectorLayoutSet> mDetectorLayouts;
  std::optional<DetectorLayoutConfigurationKey> mRequiredDetectorLayoutConfiguration;
  DetectorGeometryEpoch mRequiredDetectorGeometryEpoch{InitialDetectorGeometryEpoch};
#endif
};

template <int NLayers>
gsl::span<const Vertex> TimeFrame<NLayers>::getPrimaryVertices(int layer, int rofId) const
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const auto& entry = mROFVertexLookupTableView.getVertices(layer, rofId);
  return {&mPrimaryVertices[entry.getFirstEntry()], static_cast<gsl::span<const Vertex>::size_type>(entry.getEntries())};
}

template <int NLayers>
inline void TimeFrame<NLayers>::resetBeamXY(const float x, const float y, const float w)
{
  mBeamPos[0] = x;
  mBeamPos[1] = y;
  mBeamPosWeight = w;
}

template <int NLayers>
inline gsl::span<const int> TimeFrame<NLayers>::getROFrameClusters(int layerId) const
{
  return {&mROFramesClusters[layerId][0], static_cast<gsl::span<const int>::size_type>(mROFramesClusters[layerId].size())};
}

template <int NLayers>
inline gsl::span<Cluster> TimeFrame<NLayers>::getClustersOnLayer(int rofId, int layerId)
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const Cluster> TimeFrame<NLayers>::getClustersOnLayer(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<const Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<uint8_t> TimeFrame<NLayers>::getUsedClustersROF(int rofId, int layerId)
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUsedClusters[layerId][startIdx], static_cast<gsl::span<uint8_t>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const uint8_t> TimeFrame<NLayers>::getUsedClustersROF(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUsedClusters[layerId][startIdx], static_cast<gsl::span<const uint8_t>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const Cluster> TimeFrame<NLayers>::getClustersPerROFrange(int rofMin, int range, int layerId) const
{
  if (rofMin < 0 || rofMin >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofMin]}; // First cluster of rofMin
  int endIdx{mROFramesClusters[layerId][o2::gpu::CAMath::Min(rofMin + range, getNrof(layerId))]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<Cluster>::size_type>(endIdx - startIdx)};
}

template <int NLayers>
inline gsl::span<const int> TimeFrame<NLayers>::getROFramesClustersPerROFrange(int rofMin, int range, int layerId) const
{
  int chkdRange{o2::gpu::CAMath::Min(range, getNrof(layerId) - rofMin)};
  return {&mROFramesClusters[layerId][rofMin], static_cast<gsl::span<int>::size_type>(chkdRange)};
}

template <int NLayers>
inline gsl::span<const int> TimeFrame<NLayers>::getNClustersROFrange(int rofMin, int range, int layerId) const
{
  int chkdRange{o2::gpu::CAMath::Min(range, getNrof(layerId) - rofMin)};
  return {&mNClustersPerROF[layerId][rofMin], static_cast<gsl::span<int>::size_type>(chkdRange)};
}

template <int NLayers>
inline int TimeFrame<NLayers>::getTotalClustersPerROFrange(int rofMin, int range, int layerId) const
{
  int startIdx{rofMin}; // First cluster of rofMin
  int endIdx{o2::gpu::CAMath::Min(rofMin + range, getNrof(layerId))};
  return mROFramesClusters[layerId][endIdx] - mROFramesClusters[layerId][startIdx];
}

template <int NLayers>
inline int TimeFrame<NLayers>::getClusterROF(int iLayer, int iCluster)
{
  return std::lower_bound(mROFramesClusters[iLayer].begin(), mROFramesClusters[iLayer].end(), iCluster + 1) - mROFramesClusters[iLayer].begin() - 1;
}

template <int NLayers>
inline gsl::span<const Cluster> TimeFrame<NLayers>::getUnsortedClustersOnLayer(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUnsortedClusters[layerId][startIdx], static_cast<gsl::span<Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<int> TimeFrame<NLayers>::getIndexTable(int rofId, int layer)
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const int tableSize = mIndexTableUtils.getNrowBins() * mIndexTableUtils.getNcolBins() + 1;
  return {&mIndexTables[layer][rofId * tableSize], static_cast<gsl::span<int>::size_type>(tableSize)};
}

template <int NLayers>
template <typename... T>
void TimeFrame<NLayers>::addClusterToLayer(int layer, T&&... values)
{
  mUnsortedClusters[layer].emplace_back(std::forward<T>(values)...);
}

template <int NLayers>
template <typename... T>
void TimeFrame<NLayers>::addTrackingFrameInfoToLayer(int layer, T&&... values)
{
  mTrackingFrameInfo[layer].emplace_back(std::forward<T>(values)...);
}

template <int NLayers>
inline gsl::span<uint8_t> TimeFrame<NLayers>::getUsedClusters(const int layer)
{
  return {&mUsedClusters[layer][0], static_cast<gsl::span<uint8_t>::size_type>(mUsedClusters[layer].size())};
}

template <int NLayers>
inline gsl::span<int> TimeFrame<NLayers>::getNTrackletsCluster(int rofId, int combId)
{
  if (rofId < 0 || rofId >= getNrof(1)) {
    return {};
  }
  auto startIdx{mROFramesClusters[1][rofId]};
  return {&mNTrackletsPerCluster[combId][startIdx], static_cast<gsl::span<int>::size_type>(mROFramesClusters[1][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<int> TimeFrame<NLayers>::getExclusiveNTrackletsCluster(int rofId, int combId)
{
  if (rofId < 0 || rofId >= getNrof(1)) {
    return {};
  }
  auto clusStartIdx{mROFramesClusters[1][rofId]};

  return {&mNTrackletsPerClusterSum[combId][clusStartIdx], static_cast<gsl::span<int>::size_type>(mROFramesClusters[1][rofId + 1] - clusStartIdx)};
}

template <int NLayers>
inline gsl::span<Tracklet> TimeFrame<NLayers>::getFoundTracklets(int rofId, int combId)
{
  if (rofId < 0 || rofId >= getNrof(1) || mTracklets[combId].empty()) {
    return {};
  }
  auto startIdx{mNTrackletsPerROF[combId][rofId]};
  return {&mTracklets[combId][startIdx], static_cast<gsl::span<Tracklet>::size_type>(mNTrackletsPerROF[combId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const Tracklet> TimeFrame<NLayers>::getFoundTracklets(int rofId, int combId) const
{
  if (rofId < 0 || rofId >= getNrof(1)) {
    return {};
  }
  auto startIdx{mNTrackletsPerROF[combId][rofId]};
  return {&mTracklets[combId][startIdx], static_cast<gsl::span<Tracklet>::size_type>(mNTrackletsPerROF[combId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const MCCompLabel> TimeFrame<NLayers>::getLabelsFoundTracklets(int rofId, int combId) const
{
  if (rofId < 0 || rofId >= getNrof(1) || !hasMCinformation()) {
    return {};
  }
  auto startIdx{mNTrackletsPerROF[combId][rofId]};
  return {&mTrackletLabels[combId][startIdx], static_cast<gsl::span<Tracklet>::size_type>(mNTrackletsPerROF[combId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline int TimeFrame<NLayers>::getTotalClusters() const
{
  size_t totalClusters{0};
  for (const auto& clusters : mUnsortedClusters) {
    totalClusters += clusters.size();
  }
  return int(totalClusters);
}

template <int NLayers>
inline size_t TimeFrame<NLayers>::getNumberOfClusters() const
{
  size_t nClusters{0};
  for (const auto& layer : mClusters) {
    nClusters += layer.size();
  }
  return nClusters;
}

template <int NLayers>
inline size_t TimeFrame<NLayers>::getNumberOfCells() const
{
  size_t nCells{0};
  for (const auto& layer : mCells) {
    nCells += layer.size();
  }
  return nCells;
}

template <int NLayers>
inline size_t TimeFrame<NLayers>::getNumberOfTracklets() const
{
  size_t nTracklets{0};
  for (const auto& layer : mTracklets) {
    nTracklets += layer.size();
  }
  return nTracklets;
}

template <int NLayers>
inline size_t TimeFrame<NLayers>::getNumberOfNeighbours() const
{
  size_t neigh{0};
  for (const auto& l : mCellsNeighbours) {
    neigh += l.size();
  }
  return neigh;
}

template <int NLayers>
inline size_t TimeFrame<NLayers>::getNumberOfTracks() const
{
  return mTracks.size();
}

template <int NLayers>
inline size_t TimeFrame<NLayers>::getNumberOfUsedClusters() const
{
  size_t nClusters = 0;
  for (const auto& layer : mUsedClusters) {
    nClusters += std::count(layer.begin(), layer.end(), true);
  }
  return nClusters;
}

template <int NLayers>
inline const TrackingFrameInfo& TimeFrame<NLayers>::getClusterTrackingFrameInfo(int layerId, const Cluster& cl) const
{
  return mTrackingFrameInfo[layerId][cl.clusterId];
}

using TimeFrameITS = TimeFrame<ITSNLayers>;
/// MFT CA TimeFrame: NLayers = half-disk CA layers (see MFTFwdTrackHelpers.h).
using TimeFrameMFT = TimeFrame<o2::mft::constants::mft::LayersNumber>;

} // namespace itsmft::tracking
} // namespace o2

#endif /* ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_ */
