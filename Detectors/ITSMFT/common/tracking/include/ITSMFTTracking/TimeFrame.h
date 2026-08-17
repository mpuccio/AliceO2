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
#include "ITSMFTTracking/GenericTrack.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/GlobalMeasurement.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceLayout.h"
#include "ITSMFTTracking/SurfaceTiming.h"
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
  }

  float getBeamX() const { return mBeamPos[0]; }
  float getBeamY() const { return mBeamPos[1]; }
  std::array<float, 2>& getBeamXY() { return mBeamPos; }

  void setBz(float bz) { mBz = bz; }
  float getBz() const { return mBz; }

  gsl::span<const SurfaceMeasurement> getSurfaceMeasurements(SurfaceId surface) const;
  gsl::span<const GlobalMeasurement> getGlobalMeasurements(SurfaceId surface) const;
  const GlobalMeasurement* getGlobalMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept;
  const SurfaceMeasurement* getSurfaceMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept;
  gsl::span<const ROFIntervalBC> getSourceIntervals(ClusterSourceId source) const;
  gsl::span<const o2::MCCompLabel> getLabels(ClusterRef cluster) const;
  uint32_t getNMeasurementSurfaces() const noexcept { return static_cast<uint32_t>(mPerSurfaceMeasurements.size()); }
  std::size_t getTotalMeasurements() const noexcept;

  // Loader staging; committed only after decoding and workspace backfill.
  void assignLoadedMeasurements(std::vector<std::vector<GlobalMeasurement>>&& perSurfaceGlobalMeasurements,
                                std::vector<std::vector<SurfaceMeasurement>>&& perSurfaceMeasurements,
                                std::vector<ROFIntervalBC>&& rofIntervals,
                                std::vector<uint32_t>&& sourceROFOffsets,
                                std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>&& labelSources);
  void commitMeasurements(TimeFrame& staged) noexcept;

  // Atomically replace the event after preflight.
  bool commitLoadedEvent(TimeFrame& staged,
                         std::unique_ptr<SurfaceTrackingScratch>&& stagedWorkspace) noexcept;

  // Clear event state while preserving configuration and allocator identity.
  void resetTimeFrame() noexcept;
  std::size_t getEventResetCount() const noexcept { return mEventResetCount; }

  SurfaceTrackingScratch& getWorkspace();
  const SurfaceTrackingScratch& getWorkspace() const;

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
  // Flat inner-to-outer references; indices are local to their surface.
  auto& getTrackClusterIndices() { return mTrackClusterIndices; }
  const auto& getTrackClusterIndices() const { return mTrackClusterIndices; }

  /// memory management
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  auto& getMemoryPool() const noexcept { return mMemoryPool; }

  // Must outlive containers allocated from it (reverse destruction order).
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;

 private:
  float mBz = 5.;
  unsigned int mNTotalLowPtVertices = 0;
  int mBeamPosWeight = 0;
  std::array<float, 2> mBeamPos = {0.f, 0.f};
  bool isBeamPositionOverridden = false;

  bounded_vector<Vertex> mPrimaryVertices;
  bounded_vector<VertexLabel> mPrimaryVerticesLabels;

  bounded_vector<GenericTrack> mGenericTracks;
  bounded_vector<TrackClusterReference> mTrackClusterIndices;

  std::vector<std::vector<GlobalMeasurement>> mPerSurfaceGlobalMeasurements;
  std::vector<std::vector<SurfaceMeasurement>> mPerSurfaceMeasurements;
  std::vector<ROFIntervalBC> mROFIntervals;
  std::vector<uint32_t> mSourceROFOffsets;
  std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*> mLabelSources;

  bool mConfigurationValid = false;
  std::vector<SurfaceLayout> mLayouts;
  std::vector<TrackingParameters> mTrackingParameters;
  std::vector<TrackingWorkspaceCapacity> mWorkspaceCapacities;
  struct WorkspaceDeleter {
    void operator()(SurfaceTrackingScratch* workspace) const noexcept;
  };
  std::unique_ptr<SurfaceTrackingScratch, WorkspaceDeleter> mWorkspace;
  std::size_t mEventResetCount{0};

  void swapMeasurements(TimeFrame& other) noexcept;
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_ */
