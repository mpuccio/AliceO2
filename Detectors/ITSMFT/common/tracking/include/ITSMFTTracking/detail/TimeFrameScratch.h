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
/// \file TimeFrameScratch.h
/// \brief Runtime-plan-owned, detector-neutral CA workspace.
///
/// Host storage follows the runtime surface graph; device capacities remain
/// fixed. TimeFrame owns the workspace, while adapters own raw ROFs.
#ifndef ALICEO2_ITSMFT_TRACKING_TimeFrameScratch_H_
#define ALICEO2_ITSMFT_TRACKING_TimeFrameScratch_H_

#include <algorithm>
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
#include "ITSMFTTracking/GenericTrack.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/ROFViews.h"
#include "ITSMFTTracking/TraversalTopology.h"
#include "ITSMFTTracking/detail/TrackingKernelParameters.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/ClusterLines.h"
#include "ITStracking/ExternalAllocator.h"
#include "ITStracking/Tracklet.h"
#ifndef GPUCA_GPUCODE
#include <optional>
#include "ITSMFTTracking/SurfaceDescriptor.h"
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

class ClusterDecoder;
struct ClusterSourceInput;
struct LoadSourcesResult;

// Per-iteration derived traversal state. The TimeFrame workspace owns this
// storage; Tracker builds a short-lived view for a traversal call.
struct TraversalWorkspace {
  TrackingKernelParameters kernelParameters{};
  AttachHitConfigView attachHitConfig{};
  std::vector<NominalSurfaceMaterial> layerMaterial;
  std::vector<gsl::span<const SurfaceMeasurement>> layerMeasurements;
  std::vector<gsl::span<const GlobalMeasurement>> layerGlobalMeasurements;
  std::vector<float> diskLayerReferenceZ;
  gsl::span<const float> diskLayerReferenceZView{};
  o2::its::bounded_vector<TrackingCandidate> acceptedTracks;
  // Per-pass traversal plan. The graph remains immutable configuration; the
  // tracker derives this selected topology and its compact scratch mapping.
  SurfaceMask activeSurfaces{};
  std::vector<LayerId> orderedSurfaces;
  std::vector<int16_t> surfaceSlotById;
  std::vector<int16_t> edgeSlotById;
  std::vector<int16_t> cellSlotById;
  std::vector<EdgeId> edges;
  std::vector<CellPathId> cells;
  std::vector<CellPathId> roadStartCells;
  std::vector<uint32_t> roadStartComponentOffsets;
  std::vector<CellPathId> scheduledCells;
  SurfaceCatalogView topologyCatalog{};
  TraversalTopology topology;
  bool valid{false};

  std::optional<uint16_t> getSurfaceSlot(LayerId id) const noexcept;
  std::optional<uint16_t> getEdgeSlot(EdgeId id) const noexcept;
  std::optional<uint16_t> getCellSlot(CellPathId id) const noexcept;
  TraversalTopologyView getTopologyView() const noexcept { return topology.getView(topologyCatalog); }

  void reset(std::pmr::memory_resource* resource) noexcept
  {
    kernelParameters = {};
    attachHitConfig = {};
    layerMaterial.clear();
    layerMeasurements.clear();
    layerGlobalMeasurements.clear();
    diskLayerReferenceZ.clear();
    diskLayerReferenceZView = {};
    o2::its::deepVectorClear(acceptedTracks, resource);
    activeSurfaces = {};
    orderedSurfaces = {};
    surfaceSlotById.clear();
    edgeSlotById.clear();
    cellSlotById.clear();
    edges.clear();
    cells.clear();
    roadStartCells.clear();
    roadStartComponentOffsets.clear();
    scheduledCells.clear();
    topologyCatalog = {};
    topology = {};
    valid = false;
  }
};

/// Detector-neutral CA working state with current-event views.
class TimeFrameScratch
{
 private:
  // Pools must outlive allocator-backed members.
  std::shared_ptr<o2::its::BoundedMemoryResource> mExtMemoryPool;
  std::shared_ptr<o2::its::BoundedMemoryResource> mMemoryPool;
  o2::its::ExternalAllocator* mExternalAllocator{nullptr};
  bool mIsStaggered{false};

 public:
  TimeFrameScratch() = default;
  ~TimeFrameScratch() = default;
  TimeFrameScratch(const TimeFrameScratch&) = delete;
  TimeFrameScratch& operator=(const TimeFrameScratch&) = delete;
  TimeFrameScratch(TimeFrameScratch&&) = delete;
  TimeFrameScratch& operator=(TimeFrameScratch&&) = delete;

  /// Size surface, edge and cell storage; setMemoryPool() comes first.
  void adoptPlan(std::size_t nOwnedSurfaces, std::size_t nEdges, std::size_t nCells);
  void configureTraversalWorkspaces(std::size_t nIterations);
  TraversalWorkspace& getTraversalWorkspace(std::size_t iteration) { return mTraversalWorkspaces.at(iteration); }
  const TraversalWorkspace& getTraversalWorkspace(std::size_t iteration) const { return mTraversalWorkspaces.at(iteration); }
  std::size_t getNTraversalWorkspaces() const noexcept { return mTraversalWorkspaces.size(); }

  std::size_t getNOwnedSurfaces() const noexcept { return mNOwnedSurfaces; }
  std::size_t getNEdges() const noexcept { return mNEdges; }
  std::size_t getNCells() const noexcept { return mNCells; }

  /// Clear event state without changing plan sizes.
  void reset();

  /// Reseat allocator-backed containers.
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

  /// Atomic-load staging requires matching flat-container allocators; nested
  /// vector swaps preserve their own allocator identity.
  bool allocatorsMatch(const TimeFrameScratch& staged) const noexcept;

  /// Precondition: allocatorsMatch(other); owner-bound allocators are retained.
  void swap(TimeFrameScratch& other) noexcept;

  // ---- Read-in data: loops use the runtime ordered-surface span.
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
                                         gsl::span<const LayerId> orderedSurfaces,
                                         SurfaceCatalogView catalogView,
                                         bool applySysErrors = true);
  LoadSourcesResult backfillNormalizedSources(const TimeFrame& measurements,
                                              gsl::span<const ClusterSourceInput> sources,
                                              gsl::span<const LayerId> orderedSurfaces,
                                              SurfaceCatalogView catalog);
#endif

  int getTotalClusters() const;
  bool empty() const { return getTotalClusters() == 0; }
  int getSortedIndex(int rofId, int layer, int idx) const { return mROFramesClusters[layer][rofId] + idx; }
  int getSortedStartIndex(const int rofId, const int layer) const { return mROFramesClusters[layer][rofId]; }
  int getNrof(int layer) const { return static_cast<int>(mROFramesClusters[layer].size()) - 1; }
#ifndef GPUCA_GPUCODE
  std::optional<ClusterSourceId> getSurfaceSource(int layer) const noexcept
  {
    return layer >= 0 && static_cast<std::size_t>(layer) < mSourceBySurface.size() && mSourceBySurface[layer].isValid()
             ? std::optional<ClusterSourceId>{mSourceBySurface[layer]}
             : std::nullopt;
  }
  bool setSurfaceSources(gsl::span<const ClusterSourceId> sources)
  {
    if (sources.size() != mNOwnedSurfaces || std::any_of(sources.begin(), sources.end(), [](ClusterSourceId source) { return !source.isValid(); })) {
      return false;
    }
    mSourceBySurface.assign(sources.begin(), sources.end());
    return true;
  }
#endif

  auto& getMinRs() { return mMinR; }
  auto& getMaxRs() { return mMaxR; }
  float getMinR(int layer) const { return mMinR[layer]; }
  float getMaxR(int layer) const { return mMaxR[layer]; }
  float getMinZ(int layer) const { return mMinZ[layer]; }
  float getMaxZ(int layer) const { return mMaxZ[layer]; }
  float getEdgePhiCut(int edgeId) const { return mEdgePhiCuts[edgeId]; }
  float getEdgeMSAngle(int edgeId) const { return mEdgeMSAngles[edgeId]; }
  auto& getEdgePhiCuts() { return mEdgePhiCuts; }
  auto& getEdgeMSAngles() { return mEdgeMSAngles; }
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

  // Navigation and event timing views. Adapters own detector tables; scratch
  // retains one non-owning runtime view for the current event.
  const auto& getIndexTableUtils() const { return mIndexTableUtils.front(); }
  const auto& getIndexTableUtils(int layer) const { return mIndexTableUtils[layer]; }
  void setROFViews(RuntimeROFViews views) noexcept
  {
    mROFViews = views;
    mROFViewsBySurface.assign(mNOwnedSurfaces, views);
    mROFLocalLayerBySurface.resize(mNOwnedSurfaces);
    for (uint16_t layer = 0; layer < mNOwnedSurfaces; ++layer) {
      mROFLocalLayerBySurface[layer] = layer;
    }
    mUseUPC = false;
  }
  const RuntimeROFViews& getROFViews() const noexcept { return mROFViews; }
  const RuntimeROFViews& getROFViews(int layer) const noexcept { return mROFViewsBySurface.empty() ? mROFViews : mROFViewsBySurface[layer]; }
  int getROFLocalLayer(int layer) const noexcept { return mROFLocalLayerBySurface.empty() ? layer : mROFLocalLayerBySurface[layer]; }
  const ROFTimingLayer& getROFTiming(int layer) const noexcept { return getROFViews(layer).overlap.getLayer(getROFLocalLayer(layer)); }
  const RuntimeROFTableEntry& getROFOverlap(int fromLayer, int toLayer, int rof) const noexcept
  {
    return getROFViews(fromLayer).overlap.getOverlap(getROFLocalLayer(fromLayer), getROFLocalLayer(toLayer), rof);
  }
  bool isROFEnabled(int layer, int rof) const noexcept
  {
    const auto& views = getROFViews(layer);
    const auto& mask = mUseUPC ? views.upcMask : views.mask;
    return mask.isROFEnabled(getROFLocalLayer(layer), rof);
  }
  bool isVertexCompatible(int layer, int rof, const Vertex& vertex) const noexcept
  {
    return getROFViews(layer).vertexLookup.isVertexCompatible(getROFLocalLayer(layer), rof, vertex);
  }
  o2::its::TimeEstBC getROFTimeStamp(int fromLayer, int fromROF, int toLayer, int toROF) const noexcept
  {
    return getROFViews(fromLayer).overlap.getTimeStamp(getROFLocalLayer(fromLayer), fromROF,
                                                       getROFLocalLayer(toLayer), toROF);
  }
  int getMaxVerticesPerROF() const noexcept
  {
    int result = 0;
    if (mROFViewsBySurface.empty()) {
      return mROFViews.vertexLookup.getMaxVerticesPerROF();
    }
    for (const auto& views : mROFViewsBySurface) {
      result = std::max(result, views.vertexLookup.getMaxVerticesPerROF());
    }
    return result;
  }
  const RuntimeROFOverlapView& getROFOverlapView() const noexcept { return mROFViews.overlap; }
  const RuntimeROFVertexLookupView& getROFVertexLookupView() const noexcept { return mROFViews.vertexLookup; }
  const RuntimeROFMaskView& getROFMaskView() const noexcept { return mUseUPC ? mROFViews.upcMask : mROFViews.mask; }
  void useUPCMask() noexcept { mUseUPC = true; }
  gsl::span<const Vertex> getPrimaryVertices(const TimeFrame& frame, int layer, int rofId) const;

  const o2::its::TrackingFrameInfo& getClusterTrackingFrameInfo(int layerId, const o2::its::Cluster& cl) const;
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const o2::its::Cluster& cl) const { return getClusterLabels(layerId, cl.clusterId); }
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const int clId) const { return mClusterLabels[(mIsStaggered ? layerId : 0)]->getLabels(mClusterExternalIndices[layerId][clId]); }
  bool hasClusterExternalIndex(int layerId, int clId) const noexcept
  {
    return layerId >= 0 && layerId < static_cast<int>(mClusterExternalIndices.size()) && clId >= 0 &&
           clId < static_cast<int>(mClusterExternalIndices[layerId].size());
  }
  int getClusterExternalIndex(int layerId, const int clId) const { return mClusterExternalIndices[layerId][clId]; }
  int getClusterSize(int layer, int clusterId) const { return mClusterSize[layer][clusterId]; }
  void setClusterSize(int layer, o2::its::bounded_vector<uint8_t>& v) { mClusterSize[layer] = std::move(v); }

  auto& getTrackletsLabel(int layer) { return mTrackletLabels[layer]; }
  auto& getCellsLabel(int layer) { return mCellLabels[layer]; }

  void initialise(const TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers, int iteration,
                  const IndexTableUtilsCore& indexTableConfig, TraversalTopologyView topology,
                  gsl::span<const EdgeId> edgeIds, gsl::span<const CellPathId> cellIds,
                  gsl::span<const LayerId> orderedSurfaces,
                  gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements);
  void initialise(const TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers, int iteration,
                  gsl::span<const IndexTableUtilsCore> indexTableConfigs, TraversalTopologyView topology,
                  gsl::span<const EdgeId> edgeIds, gsl::span<const CellPathId> cellIds,
                  gsl::span<const LayerId> orderedSurfaces,
                  gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements);

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
  const auto& getCells() const { return mCells; }

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
  std::vector<float> mMinZ;
  std::vector<float> mMaxZ;
  o2::its::bounded_vector<int> mBogusClusters;
  o2::its::bounded_vector<float> mPositionResolution;
  // Two fixed tracklet combinations (positions 0-1 and 1-2), sized at load.
  std::array<o2::its::bounded_vector<int>, 2> mNTrackletsPerCluster;
  std::array<o2::its::bounded_vector<int>, 2> mNTrackletsPerClusterSum;

  // ---- Group B: plan-sized CA construction/result transients ----
  std::vector<o2::its::bounded_vector<o2::its::Tracklet>> mTracklets;
  std::vector<o2::its::bounded_vector<int>> mTrackletsLookupTable;
  std::vector<o2::its::bounded_vector<o2::MCCompLabel>> mTrackletLabels;
  o2::its::bounded_vector<float> mEdgePhiCuts;
  o2::its::bounded_vector<float> mEdgeMSAngles;
  std::vector<o2::its::bounded_vector<CellSeed>> mCells;
  std::vector<o2::its::bounded_vector<int>> mCellsLookupTable;
  std::vector<o2::its::bounded_vector<int>> mCellsNeighbours;
  std::vector<o2::its::bounded_vector<int>> mCellsNeighboursTopology;
  std::vector<o2::its::bounded_vector<int>> mCellsNeighboursLUT;
  std::vector<o2::its::bounded_vector<o2::MCCompLabel>> mCellLabels;
  // Shared navigation auxiliary.
  std::vector<IndexTableUtilsCore> mIndexTableUtils;
  std::vector<TraversalWorkspace> mTraversalWorkspaces;

  // ---- Group D: vertexer working scratch ----
  // Sized by ROF count and the two fixed pairs, not by adoptPlan().
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
  friend struct TimeFrame;

  void swapLoadedEvent(TimeFrameScratch& other) noexcept;

  void prepareClusters(const TimeFrame& frame, const TrackingParameters& trkParam, int maxLayers,
                       gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements);

  RuntimeROFViews mROFViews{};
  std::vector<RuntimeROFViews> mROFViewsBySurface;
  std::vector<uint16_t> mROFLocalLayerBySurface;
  std::vector<ClusterSourceId> mSourceBySurface;
  bool mUseUPC{false};

  std::size_t mNOwnedSurfaces{0};
  std::size_t mNEdges{0};
  std::size_t mNCells{0};
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TimeFrameScratch_H_ */
