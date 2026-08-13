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
/// TimeFrame owns configuration, event measurements, generic results, one
/// tracking workspace, and allocator/capacity state. The application owns
/// raw ROFs, timing tables, publication state, sidecars, and workflow state.

#ifndef ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_
#define ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_

#include <memory>
#include <array>
#include <cstddef>
#include <vector>

#include <gsl/gsl>

#include "DataFormatsITS/Vertex.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "ITSMFTTracking/CommonTrack.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/MeasurementView.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/Cluster.h"

namespace o2::itsmft::tracking
{

class SurfaceTrackingScratch;

using BoundedMemoryResource = o2::its::BoundedMemoryResource;
using Vertex = o2::its::Vertex;
using VertexLabel = o2::its::VertexLabel;
template <typename T>
using bounded_vector = o2::its::bounded_vector<T>;

struct TrackingWorkspaceCapacity {
  std::size_t ownedSurfaces = 0;
  std::size_t transitions = 0;
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
  }

  float getBeamX() const { return mBeamPos[0]; }
  float getBeamY() const { return mBeamPos[1]; }
  std::array<float, 2>& getBeamXY() { return mBeamPos; }

  void setBz(float bz) { mBz = bz; }
  float getBz() const { return mBz; }

  MeasurementView getMeasurementView() const noexcept;
  gsl::span<const SurfaceMeasurement> getSurfaceMeasurements(SurfaceId surface) const;
  gsl::span<const GlobalMeasurement> getGlobalMeasurements(SurfaceId surface) const;
  const GlobalMeasurement* getGlobalMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept;
  const SurfaceMeasurement* getSurfaceMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept;
  gsl::span<const ROFIntervalBC> getSourceIntervals(ClusterSourceId source) const;
  gsl::span<const o2::MCCompLabel> getLabels(ClusterRef cluster) const;
  uint32_t getNMeasurementSurfaces() const noexcept { return static_cast<uint32_t>(mPerSurfaceMeasurements.size()); }
  std::size_t getTotalMeasurements() const noexcept;

  // Loader-only staging: commit swaps staged measurements after decoding and
  // workspace backfill succeed.
  void assignLoadedMeasurements(std::vector<std::vector<GlobalMeasurement>>&& perSurfaceGlobalMeasurements,
                                std::vector<std::vector<SurfaceMeasurement>>&& perSurfaceMeasurements,
                                std::vector<ROFIntervalBC>&& rofIntervals,
                                std::vector<uint32_t>&& sourceROFOffsets,
                                std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>&& labelSources);
  void commitMeasurements(TimeFrame& staged) noexcept;

  // Commit a preflighted multi-source event, resetting the previous event once
  // before installing measurements and the staged workspace.
  bool commitLoadedEvent(TimeFrame& staged,
                         std::unique_ptr<SurfaceTrackingScratch>&& stagedWorkspace) noexcept;

  // Reset event state while preserving configuration and allocator identity.
  // External owners call this once per event.
  void resetEvent() noexcept;
  std::size_t getEventResetCount() const noexcept { return mEventResetCount; }

  SurfaceTrackingScratch& getWorkspace();
  const SurfaceTrackingScratch& getWorkspace() const;

  bool commitConfiguration(std::vector<SurfaceGraph>&& graphs,
                           std::vector<TrackingParameters>&& parameters,
                           std::vector<std::unique_ptr<SurfacePlanBinding>>&& bindings,
                           std::vector<TrackingWorkspaceCapacity>&& capacities,
                           std::shared_ptr<BoundedMemoryResource> memoryPool);
  bool isConfigured() const noexcept { return mConfigurationValid; }
  std::size_t getNIterations() const noexcept { return mGraphs.size(); }
  const std::vector<SurfaceGraph>& getGraphs() const noexcept { return mGraphs; }
  const std::vector<TrackingParameters>& getTrackingParameters() const noexcept;
  const SurfaceGraph& getGraph(std::size_t iteration) const { return mGraphs.at(iteration); }
  const TrackingParameters* getTrackingParameters(std::size_t iteration) const noexcept;
  const SurfacePlanBinding* getBinding(std::size_t iteration) const noexcept;
  const TrackingWorkspaceCapacity* getWorkspaceCapacity(std::size_t iteration) const noexcept;

  // Detector-neutral common-CA results (see CommonTrack.h). Valid only with
  // the current event measurements; CommonTrack has no NLayers dependency.
  auto& getCommonTracks() { return mCommonTracks; }
  const auto& getCommonTracks() const { return mCommonTracks; }
  // Flat TrackClusterReference array. A CommonTrack's
  // [firstClusterRef, clusterRefEnd) range contains positions in inner-to-
  // outer traversal order. Each entry pairs a SurfaceId with an index local
  // to that surface's measurements, resolved by getGlobalMeasurement(); it is
  // not a global measurement position.
  auto& getTrackClusterIndices() { return mTrackClusterIndices; }
  const auto& getTrackClusterIndices() const { return mTrackClusterIndices; }

  /// memory management
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  auto& getMemoryPool() const noexcept { return mMemoryPool; }

  // Declared before vectors that allocate from it so reverse destruction order
  // releases those allocations before the resource. TimeFrame containers do
  // not use the framework/GPU allocator; SurfaceTrackingScratch owns that
  // compatibility path for its per-layer containers.
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;

 private:
  float mBz = 5.;
  unsigned int mNTotalLowPtVertices = 0;
  int mBeamPosWeight = 0;
  std::array<float, 2> mBeamPos = {0.f, 0.f};
  bool isBeamPositionOverridden = false;

  bounded_vector<Vertex> mPrimaryVertices;
  bounded_vector<VertexLabel> mPrimaryVerticesLabels;

  // Common-CA results, valid only with normalized content from the same event.
  bounded_vector<CommonTrack> mCommonTracks;
  bounded_vector<TrackClusterReference> mTrackClusterIndices;

  std::vector<std::vector<GlobalMeasurement>> mPerSurfaceGlobalMeasurements;
  std::vector<std::vector<SurfaceMeasurement>> mPerSurfaceMeasurements;
  std::vector<SurfaceMeasurementSpan> mMeasurementSpans;
  std::vector<ROFIntervalBC> mROFIntervals;
  std::vector<uint32_t> mSourceROFOffsets;
  std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*> mLabelSources;

  bool mConfigurationValid = false;
  std::vector<SurfaceGraph> mGraphs;
  std::vector<TrackingParameters> mTrackingParameters;
  std::vector<std::unique_ptr<SurfacePlanBinding>> mBindings;
  std::vector<TrackingWorkspaceCapacity> mWorkspaceCapacities;
  struct WorkspaceDeleter {
    void operator()(SurfaceTrackingScratch* workspace) const noexcept;
  };
  std::unique_ptr<SurfaceTrackingScratch, WorkspaceDeleter> mWorkspace;
  std::size_t mEventResetCount{0};

  void clearEventData() noexcept;
  void clearMeasurements() noexcept;
  void rebuildMeasurementSpans();
  void swapMeasurements(TimeFrame& other) noexcept;
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_ */
