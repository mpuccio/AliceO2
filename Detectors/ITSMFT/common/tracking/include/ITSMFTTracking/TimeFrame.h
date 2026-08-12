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
/// TimeFrame owns configuration, event measurements, generic results,
/// one tracking workspace, and its allocator/capacity state. Raw ROFs,
/// timing-table storage, publication state, typed sidecars, and workflow
/// lifecycle remain with the application boundary.

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

  // Loader-only staging hooks. A staged TimeFrame carries measurement data
  // only; commit swaps it after all decoding and workspace backfill succeed.
  void assignLoadedMeasurements(std::vector<std::vector<GlobalMeasurement>>&& perSurfaceGlobalMeasurements,
                                std::vector<std::vector<SurfaceMeasurement>>&& perSurfaceMeasurements,
                                std::vector<ROFIntervalBC>&& rofIntervals,
                                std::vector<uint32_t>&& sourceROFOffsets,
                                std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>&& labelSources);
  void commitMeasurements(TimeFrame& staged) noexcept;

  // Commits a complete multi-source event after loader preflight. The
  // operation resets the prior event exactly once before installing all
  // measurement data and the staged global workspace.
  bool commitLoadedEvent(TimeFrame& staged,
                         std::unique_ptr<SurfaceTrackingScratch>&& stagedWorkspace) noexcept;

  // One generic event-state reset. It preserves static configuration and
  // allocator/capacity identity while clearing the tracking workspace and
  // event result. External owners call this once for a whole event.
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

  // Detector-neutral common-CA result storage (ITSMFTTracking/CommonTrack.h).
  // CommonTrack itself carries no NLayers dependency. Only meaningful
  // together with the current event measurements -- see
  // CommonTrack.h's own lifetime doc.
  auto& getCommonTracks() { return mCommonTracks; }
  const auto& getCommonTracks() const { return mCommonTracks; }
  // Flat array of TrackClusterReference (CommonTrack.h); a CommonTrack's
  // [firstClusterRef, clusterRefEnd) range is a half-open range of
  // *positions* into this array, in traversal order (inner to outer). Each
  // element pairs a SurfaceId with a SurfaceMeasurementIndex local to that
  // surface's own measurement array -- resolved via
  // getGlobalMeasurement(reference.surface, reference.index)
  // -- never a global/flattened measurement position.
  auto& getTrackClusterIndices() { return mTrackClusterIndices; }
  const auto& getTrackClusterIndices() const { return mTrackClusterIndices; }

  /// memory management
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  auto& getMemoryPool() const noexcept { return mMemoryPool; }

  // Must be declared -- and therefore destroyed -- before every
  // pmr/bounded_vector member that may allocate through it (mPrimaryVertices,
  // mPrimaryVerticesLabels, mCommonTracks, mTrackClusterIndices). C++
  // destroys non-static data members in reverse declaration order, so
  // declaring this resource owner first guarantees it is destroyed last,
  // after every vector that could still hold memory allocated from it has
  // already released it back. TimeFrame's own containers never use the
  // framework/GPU allocator: that mechanism (mExtMemoryPool,
  // hasFrameworkAllocator(), getMaybeFrameworkHostResource()) is used
  // exclusively by SurfaceTrackingScratch's per-layer compatibility
  // containers and lives there instead.
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;

 private:
  float mBz = 5.;
  unsigned int mNTotalLowPtVertices = 0;
  int mBeamPosWeight = 0;
  std::array<float, 2> mBeamPos = {0.f, 0.f};
  bool isBeamPositionOverridden = false;

  bounded_vector<Vertex> mPrimaryVertices;
  bounded_vector<VertexLabel> mPrimaryVerticesLabels;

  // Detector-neutral common-CA result storage. Results are valid only with
  // the normalized frame content from the same event.
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
