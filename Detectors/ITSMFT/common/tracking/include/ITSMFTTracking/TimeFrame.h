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
/// \file TimeFrame.h
/// \brief Passive common event owner.
///
/// TimeFrame owns configuration, event measurements and navigation, generic
/// results, one tracking workspace, and allocator/capacity state. The
/// application owns raw ROFs, publication state, and workflow state.

#ifndef ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_
#define ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_

#include <memory>
#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include <gsl/gsl>

#include "DataFormatsITS/Vertex.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "ITSMFTTracking/GenericTrack.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/GlobalMeasurement.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceLayout.h"
#include "ITSMFTTracking/TrackingPrimitives.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/ROFViews.h"
#include "ITStracking/BoundedAllocator.h"

namespace o2::itsmft::tracking
{

class TimeFrameScratch;

using BoundedMemoryResource = o2::its::BoundedMemoryResource;
using Vertex = o2::its::Vertex;
using VertexLabel = o2::its::VertexLabel;
template <typename T>
using bounded_vector = o2::its::bounded_vector<T>;

struct TrackingWorkspaceCapacity {
  std::size_t ownedSurfaces = 0;
  std::size_t edges = 0;
  std::size_t cells = 0;
};

struct TimeFrame {
  TimeFrame() = default;
  TimeFrame(const TimeFrame&) = delete;
  TimeFrame& operator=(const TimeFrame&) = delete;
  virtual ~TimeFrame();

  const Vertex& getPrimaryVertex(const int ivtx) const { return mPrimaryVertices[ivtx]; }
  auto& getPrimaryVertices() { return mPrimaryVertices; };
  auto getPrimaryVerticesNum() { return mPrimaryVertices.size(); };
  const auto& getPrimaryVertices() const { return mPrimaryVertices; };
  auto& getPrimaryVerticesLabels() { return mPrimaryVerticesLabels; };
  void addPrimaryVertex(const Vertex& vertex);
  void addPrimaryVertexLabel(const VertexLabel& label) { mPrimaryVerticesLabels.push_back(label); }

  void resetBeamXY(const float x, const float y, const float w = 0);
  void setBeamPosition(const float x, const float y, const float s2, const float base = 50.f, const float systematic = 0.f)
  {
    isBeamPositionOverridden = true;
    resetBeamXY(x, y, s2 / o2::gpu::CAMath::Sqrt((base * base) + systematic));
    mBeamPositionVariance = s2;
  }

  float getBeamX() const { return mBeamPos[0]; }
  float getBeamY() const { return mBeamPos[1]; }
  float getBeamPositionVariance() const { return mBeamPositionVariance; }
  std::array<float, 2>& getBeamXY() { return mBeamPos; }

  void setBz(float bz) { mBz = bz; }
  float getBz() const { return mBz; }

  gsl::span<const GlobalMeasurement> getGlobalMeasurements(LayerId surface) const;
  gsl::span<GlobalMeasurement> getGlobalMeasurements(LayerId surface);
  const SurfaceMeasurement* getSurfaceMeasurement(LayerId layer, uint32_t clusterId) const noexcept;
  gsl::span<const o2::MCCompLabel> getLabels(LayerId layer, uint32_t clusterId) const;
  uint32_t getNMeasurementSurfaces() const noexcept { return static_cast<uint32_t>(mLayerGlobalMeasurements.size()); }
  std::size_t getTotalMeasurements() const noexcept;

  int getTotalClusters() const { return static_cast<int>(getTotalMeasurements()); }
  bool empty() const { return getTotalMeasurements() == 0; }
  int getSortedIndex(int rofId, int layer, int idx) const { return mROFramesClusters[layer][rofId] + idx; }
  int getSortedStartIndex(int rofId, int layer) const { return mROFramesClusters[layer][rofId]; }
  int getNrof(int layer) const
  {
    return mROFramesClusters[layer].empty() ? 0 : static_cast<int>(mROFramesClusters[layer].size()) - 1;
  }
  gsl::span<GlobalMeasurement> getClustersOnLayer(int rofId, int layer);
  gsl::span<const GlobalMeasurement> getClustersOnLayer(int rofId, int layer) const;
  auto& getClusters() noexcept { return mLayerGlobalMeasurements; }
  const auto& getClusters() const noexcept { return mLayerGlobalMeasurements; }
  gsl::span<const GlobalMeasurement> getClustersPerROFrange(int rofMin, int range, int layer) const;
  gsl::span<const int> getROFramesClustersPerROFrange(int rofMin, int range, int layer) const;
  gsl::span<const int> getROFrameClusters(int layer) const;
  gsl::span<int> getIndexTable(int rofId, int layer);
  int getClusterROF(int layer, int cluster) const;
  int getTotalClustersPerROFrange(int rofMin, int range, int layer) const;

  bool isClusterUsed(int layer, uint32_t clusterId) const;
  void markUsedCluster(int layer, uint32_t clusterId);
  gsl::span<unsigned char> getUsedClusters(int layer);
  std::size_t getNumberOfClusters() const;
  std::size_t getNumberOfUsedClusters() const;

  float getMinR(int layer) const { return mMinR[layer]; }
  float getMaxR(int layer) const { return mMaxR[layer]; }
  float getMinZ(int layer) const { return mMinZ[layer]; }
  float getMaxZ(int layer) const { return mMaxZ[layer]; }
  const auto& getIndexTableUtils() const { return mIndexTableUtils.front(); }
  const auto& getIndexTableUtils(int layer) const { return mIndexTableUtils[layer]; }

  void setROFViews(RuntimeROFViews views) noexcept;
  const RuntimeROFViews& getROFViews() const noexcept { return mROFViews; }
  const RuntimeROFViews& getROFViews(int layer) const noexcept { return mROFViewsBySurface.empty() ? mROFViews : mROFViewsBySurface[layer]; }
  int getROFLocalLayer(int layer) const noexcept { return mROFLocalLayerBySurface.empty() ? layer : mROFLocalLayerBySurface[layer]; }
  const ROFTimingLayer& getROFTiming(int layer) const noexcept { return getROFViews(layer).overlap.getLayer(getROFLocalLayer(layer)); }
  const RuntimeROFTableEntry& getROFOverlap(int fromLayer, int toLayer, int rof) const noexcept;
  bool isROFEnabled(int layer, int rof) const noexcept;
  bool isVertexCompatible(int layer, int rof, const Vertex& vertex) const noexcept;
  o2::its::TimeEstBC getROFTimeStamp(int fromLayer, int fromROF, int toLayer, int toROF) const noexcept;
  int getMaxVerticesPerROF() const noexcept;
  const RuntimeROFOverlapView& getROFOverlapView() const noexcept { return mROFViews.overlap; }
  const RuntimeROFVertexLookupView& getROFVertexLookupView() const noexcept { return mROFViews.vertexLookup; }
  const RuntimeROFMaskView& getROFMaskView() const noexcept { return mUseUPC ? mROFViews.upcMask : mROFViews.mask; }
  void useUPCMask() noexcept { mUseUPC = true; }
  gsl::span<const Vertex> getPrimaryVertices(int layer, int rofId) const;

  bool hasMCinformation() const noexcept;
  gsl::span<const MCCompLabel> getClusterLabels(int layer, int cluster) const;

  // Loader staging; committed only after complete decoding and validation.
  void assignLoadedMeasurements(std::vector<std::vector<GlobalMeasurement>>&& perSurfaceGlobalMeasurements,
                                std::vector<std::vector<SurfaceMeasurement>>&& perSurfaceMeasurements,
                                std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>&& perSurfaceLabels,
                                bool hasMCInformation);
  void assignLoadedEventNavigation(std::vector<std::vector<int>>&& rofBoundaries,
                                   RuntimeROFViews defaultViews,
                                   std::vector<RuntimeROFViews>&& viewsBySurface,
                                   std::vector<uint16_t>&& localLayerBySurface);
  void commitMeasurements(TimeFrame& staged) noexcept;

  // Atomically replace the event after preflight.
  bool commitLoadedEvent(TimeFrame& staged) noexcept;

  // Clear event state while preserving configuration and allocator identity.
  void resetTimeFrame() noexcept;
  std::size_t getEventResetCount() const noexcept { return mEventResetCount; }

  TimeFrameScratch& getWorkspace();
  const TimeFrameScratch& getWorkspace() const;

  bool commitConfiguration(std::vector<SurfaceLayout>&& layouts,
                           std::vector<TrackingParameters>&& parameters,
                           std::vector<TrackingWorkspaceCapacity>&& capacities,
                           std::shared_ptr<BoundedMemoryResource> memoryPool);
  bool isConfigured() const noexcept { return mConfigurationValid; }
  std::size_t getNIterations() const noexcept { return mLayouts.size(); }
  const SurfaceLayout& getLayout(std::size_t iteration) const { return mLayouts.at(iteration); }
  const std::vector<TrackingParameters>& getTrackingParameters() const noexcept;
  const TrackingParameters* getTrackingParameters(std::size_t iteration) const noexcept;
  const TrackingWorkspaceCapacity* getWorkspaceCapacity(std::size_t iteration) const noexcept;

  // Results are valid only with this event's normalized measurements.
  auto& getGenericTracks() { return mGenericTracks; }
  const auto& getGenericTracks() const { return mGenericTracks; }
  // Flat inner-to-outer references; IDs are stable pre-sort positions in the
  // TimeFrame-owned per-surface arrays.
  auto& getTrackClusterIndices() { return mTrackClusterIndices; }
  const auto& getTrackClusterIndices() const { return mTrackClusterIndices; }

  /// memory management
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  auto& getMemoryPool() const noexcept { return mMemoryPool; }

 private:
  // Must outlive containers allocated from it (reverse destruction order).
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;

  // Event and cross-iteration tracking state.
  std::vector<std::vector<int>> mROFramesClusters;
  std::vector<bounded_vector<int>> mIndexTables;
  std::vector<std::vector<uint8_t>> mLayerUsedClusters;
  std::vector<IndexTableUtilsCore> mIndexTableUtils;
  std::vector<float> mMinR;
  std::vector<float> mMaxR;
  std::vector<float> mMinZ;
  std::vector<float> mMaxZ;

  RuntimeROFViews mROFViews{};
  std::vector<RuntimeROFViews> mROFViewsBySurface;
  std::vector<uint16_t> mROFLocalLayerBySurface;
  bool mUseUPC{false};

  float mBz = 5.;
  unsigned int mNTotalLowPtVertices = 0;
  int mBeamPosWeight = 0;
  std::array<float, 2> mBeamPos = {0.f, 0.f};
  float mBeamPositionVariance = 0.f;
  bool isBeamPositionOverridden = false;

  bounded_vector<Vertex> mPrimaryVertices;
  bounded_vector<VertexLabel> mPrimaryVerticesLabels;

  bounded_vector<GenericTrack> mGenericTracks;
  bounded_vector<TrackClusterReference> mTrackClusterIndices;

  std::vector<std::vector<GlobalMeasurement>> mLayerGlobalMeasurements;
  std::vector<std::vector<SurfaceMeasurement>> mLayerSurfaceMeasurements;
  std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>> mLayerClusterLabels;
  bool mHasMCInformation{false};

  bool mConfigurationValid = false;
  std::vector<SurfaceLayout> mLayouts;
  std::vector<TrackingParameters> mTrackingParameters;
  std::vector<TrackingWorkspaceCapacity> mWorkspaceCapacities;
  struct WorkspaceDeleter {
    void operator()(TimeFrameScratch* workspace) const noexcept;
  };
  std::unique_ptr<TimeFrameScratch, WorkspaceDeleter> mWorkspace;
  std::size_t mEventResetCount{0};

  void swapMeasurements(TimeFrame& other) noexcept;
  void publishConfiguration(TimeFrame& staged) noexcept;
  void configureEventStorage(std::size_t nOwnedSurfaces);
  void prepareClusters(int maxLayers);
  friend class TimeFrameScratch;
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_ */
