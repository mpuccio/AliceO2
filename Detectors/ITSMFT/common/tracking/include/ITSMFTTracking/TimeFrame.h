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
/// \file TimeFrame.h
/// \brief Passive common TimeFrame owner.
///
/// TimeFrame owns the invariant detector layout, measurements and navigation,
/// generic results, tracking scratch, and allocator state. The
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
#include "ITSMFTTracking/CapacityEstimator.h"
#include "ITSMFTTracking/GlobalMeasurement.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/DetectorConfiguration.h"
#include "ITSMFTTracking/TrackingPrimitives.h"
#include "ITSMFTTracking/IndexTableConfigurationSet.h"
#include "ITSMFTTracking/ROFViews.h"
#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITSMFTTracking/BoundedAllocator.h"

namespace o2::itsmft::tracking
{

using Vertex = o2::its::Vertex;
using VertexLabel = o2::its::VertexLabel;

// A layer of the timeframe as the CA steps read it, on the host or on the
// device: clusters sorted by ROF and index-table bin, and, by cluster id,
// surface measurements and used flags.
struct LayerView {
  const GlobalMeasurement* clusters{};
  size_t nClusters{};
  const SurfaceMeasurement* measurements{};
  size_t nMeasurements{};
  const uint8_t* used{};
  size_t nUsed{};
  const int* rofClusters{}; // first cluster of every ROF, plus the total
  int nROFs{};
  const uint8_t* rofEnabled{};
  const int* indexTables{}; // one table of tableSize entries per ROF
  size_t tableSize{};

  GPUhdi() bool isUsed(uint32_t clusterId) const { return clusterId < nUsed && used[clusterId]; }
  GPUhdi() const SurfaceMeasurement* measurement(uint32_t clusterId) const { return clusterId < nMeasurements ? measurements + clusterId : nullptr; }
};
struct FrameView {
  std::array<LayerView, MaxLayoutSurfaces> layers{};
  size_t nLayers{};
  const Vertex* vertices{};

  GPUhdi() const GlobalMeasurement* cluster(int layer, int index) const
  {
    return layer >= 0 && size_t(layer) < nLayers && index >= 0 && size_t(index) < layers[layer].nClusters ? layers[layer].clusters + index : nullptr;
  }
};

struct TimeFrame {
  TimeFrame() = default;
  TimeFrame(const TimeFrame&) = delete;
  TimeFrame& operator=(const TimeFrame&) = delete;
  virtual ~TimeFrame() = default;

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
  void addMeasurement(LayerId surface, GlobalMeasurement global,
                      const SurfaceMeasurement& measurement);
  void addMeasurement(LayerId surface, GlobalMeasurement global,
                      const SurfaceMeasurement& measurement,
                      gsl::span<const o2::MCCompLabel> labels);
  void setHasMCInformation(bool value) noexcept { mHasMCInformation = value; }
  const SurfaceMeasurement* getSurfaceMeasurement(LayerId layer, uint32_t clusterId) const noexcept;
  // A layer's surface measurements, indexed by cluster id.
  gsl::span<const SurfaceMeasurement> getSurfaceMeasurements(LayerId layer) const noexcept
  {
    return layer.isValid() && layer.value() < mLayerSurfaceMeasurements.size() ? gsl::make_span(mLayerSurfaceMeasurements[layer.value()]) : gsl::span<const SurfaceMeasurement>{};
  }
  gsl::span<const o2::MCCompLabel> getLabels(LayerId layer, uint32_t clusterId) const;
  uint32_t getNMeasurementSurfaces() const noexcept { return static_cast<uint32_t>(mLayerGlobalMeasurements.size()); }
  std::size_t getTotalMeasurements() const noexcept;

  int getTotalClusters() const { return static_cast<int>(getTotalMeasurements()); }
  bool empty() const { return getTotalMeasurements() == 0; }
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
  // Index-table bin of every cluster of a layer, in sorted cluster order,
  // and whether prepareClusters built each ROF's table (unbuilt ones stay
  // zero). Together they reproduce the index tables exactly.
  gsl::span<const int> getClusterBins(int layer) const { return mClusterBins[layer]; }
  gsl::span<const uint8_t> getIndexTableBuilt(int layer) const { return mIndexTableBuilt[layer]; }
  int getClusterROF(int layer, int cluster) const;
  int getTotalClustersPerROFrange(int rofMin, int range, int layer) const;

  // The frame as the CA steps read it, with the given sorted clusters of
  // every layer (the frame's own, as the tracker prepares them).
  FrameView makeView(gsl::span<const gsl::span<const GlobalMeasurement>> clusters) const;
  bool isClusterUsed(int layer, uint32_t clusterId) const;
  void markUsedCluster(int layer, uint32_t clusterId);
  gsl::span<unsigned char> getUsedClusters(int layer);
  std::size_t getNumberOfClusters() const;
  std::size_t getNumberOfUsedClusters() const;

  float getMinR(int layer) const { return mMinR[layer]; }
  float getMaxR(int layer) const { return mMaxR[layer]; }
  float getMinZ(int layer) const { return mMinZ[layer]; }
  float getMaxZ(int layer) const { return mMaxZ[layer]; }
  const auto& getIndexTableUtils() const { return mIndexTableUtils[0]; }
  const auto& getIndexTableUtils(int layer) const { return mIndexTableUtils[layer]; }

  void setROFViews(RuntimeROFViews views) noexcept;
  void setROFViews(std::size_t position, RuntimeROFViews views, uint16_t localLayer);
  void setROFClusters(std::size_t position, gsl::span<const int> boundaries);
  const RuntimeROFViews& getROFViews() const noexcept { return mROFViews; }
  const RuntimeROFViews& getROFViews(int layer) const noexcept { return mROFViewsBySurface.empty() ? mROFViews : mROFViewsBySurface[layer]; }
  int getROFLocalLayer(int layer) const noexcept { return mROFLocalLayerBySurface.empty() ? layer : mROFLocalLayerBySurface[layer]; }
  const ROFTimingLayer& getROFTiming(int layer) const noexcept { return getROFViews(layer).overlap.getLayer(getROFLocalLayer(layer)); }
  bool isROFEnabled(int layer, int rof) const noexcept;
  int getMaxVerticesPerROF() const noexcept;
  const RuntimeROFOverlapView& getROFOverlapView() const noexcept { return mROFViews.overlap; }
  const RuntimeROFMaskView& getROFMaskView() const noexcept { return mUseUPC ? mROFViews.upcMask : mROFViews.mask; }
  void useUPCMask() noexcept { mUseUPC = true; }
  gsl::span<const Vertex> getPrimaryVertices(int layer, int rofId) const;

  bool hasMCinformation() const noexcept;
  gsl::span<const MCCompLabel> getClusterLabels(int layer, int cluster) const;

  // Clear TimeFrame data while preserving configuration and allocator identity.
  void resetTimeFrame() noexcept;

  TimeFrameScratch& getScratch();
  const TimeFrameScratch& getScratch() const;
  CapacityEstimator& getCapacityEstimator() noexcept { return mCapacityEstimator; }
  const CapacityEstimator& getCapacityEstimator() const noexcept { return mCapacityEstimator; }

  bool configure(DetectorConfiguration&& layout, std::size_t maxEdges, std::size_t maxCells,
                 std::shared_ptr<BoundedMemoryResource> memoryPool);
  bool isConfigured() const noexcept { return mConfigurationValid; }
  const DetectorConfiguration& getDetectorConfiguration() const noexcept { return mDetectorConfiguration; }

  // Results are valid only with this TimeFrame's measurements.
  auto& getGenericTracks() { return mGenericTracks; }
  const auto& getGenericTracks() const { return mGenericTracks; }
  auto& getTrackLabels() { return mTrackLabels; }
  const auto& getTrackLabels() const { return mTrackLabels; }
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

  // TimeFrame and cross-iteration tracking state.
  std::vector<std::vector<int>> mROFramesClusters;
  std::vector<bounded_vector<int>> mIndexTables;
  std::vector<bounded_vector<int>> mClusterBins;
  std::vector<bounded_vector<uint8_t>> mIndexTableBuilt;
  std::vector<std::vector<uint8_t>> mLayerUsedClusters;
  IndexTableConfigurationSet mIndexTableUtils;
  std::vector<float> mMinR;
  std::vector<float> mMaxR;
  std::vector<float> mMinZ;
  std::vector<float> mMaxZ;

  RuntimeROFViews mROFViews{};
  std::vector<RuntimeROFViews> mROFViewsBySurface;
  std::vector<uint16_t> mROFLocalLayerBySurface;
  bool mUseUPC{false};

  float mBz = 5.;
  int mBeamPosWeight = 0;
  std::array<float, 2> mBeamPos = {0.f, 0.f};
  float mBeamPositionVariance = 0.f;
  bool isBeamPositionOverridden = false;

  bounded_vector<Vertex> mPrimaryVertices;
  bounded_vector<VertexLabel> mPrimaryVerticesLabels;

  bounded_vector<GenericTrack> mGenericTracks;
  bounded_vector<MCCompLabel> mTrackLabels;
  bounded_vector<TrackClusterReference> mTrackClusterIndices;

  std::vector<std::vector<GlobalMeasurement>> mLayerGlobalMeasurements;
  std::vector<std::vector<SurfaceMeasurement>> mLayerSurfaceMeasurements;
  std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>> mLayerClusterLabels;
  bool mHasMCInformation{false};

  bool mConfigurationValid = false;
  DetectorConfiguration mDetectorConfiguration;
  TimeFrameScratch mScratch;
  CapacityEstimator mCapacityEstimator;
  // Without host tables (a backend that builds them from the cluster bins)
  // only the bins and built flags are prepared.
  void prepareIndexTables(const IndexTableConfigurationSet& indexTableConfigs, bool hostTables = true);
  void prepareClusters(int maxLayers);
  friend class Tracker;
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_ */
