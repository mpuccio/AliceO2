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
/// \brief Gate 4 B3.1: the permanent, non-templated common event owner.
///
/// TimeFrame owns exactly the detector- and NLayers-independent event state:
/// normalized per-SurfaceId measurements (MultiSourceFrame), the detector-
/// neutral CommonTrack/TrackClusterReference result storage, primary
/// vertices, beam state, Bz, and the memory-pool ownership those event
/// vectors allocate from. It stores no NLayers-templated type, no per-layer
/// array, and no single-detector identity (mDetId is gone -- callers that
/// need "which detector" pass it explicitly, since it was always fully
/// determined by the caller's own compile-time NLayers). Generic CA workspace
/// is owned here as source-qualified private implementation state. Raw ROFs,
/// timing-table storage, publication state, and typed sidecars remain outside
/// this entity.

#ifndef ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_
#define ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_

#include <memory>
#include <array>
#include <cstddef>
#include <vector>

#include <gsl/gsl>

#include "DataFormatsITS/Vertex.h"
#include "ITSMFTTracking/CommonTrack.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
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
  using BindingSet = std::vector<std::unique_ptr<SurfacePlanBinding>>;

  TimeFrame() = default;
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

  // Non-owning, read-only access to the normalized owner/view associated
  // with this TimeFrame by the most recent successful commitNormalizedFrame()
  // call (see SurfaceTrackingScratch::loadNormalizedSource(), the
  // owner-level load operation this is the TimeFrame-side half of). Empty/
  // default until that first succeeds, and after resetEvent(): resetEvent()
  // unconditionally clears this owner in place. The `const MultiSourceFrame&`
  // returned by getNormalizedFrame() is a reference to that same long-lived
  // member object -- it remains valid and safe to dereference across a
  // resetEvent() call, it simply then observes the owner's newly cleared
  // (empty) state. What resetEvent() does invalidate is any MultiSourceFrameView or
  // gsl::span (getSurfaceMeasurements(), getSourceIntervals(), getLabels(),
  // getView()) obtained *before* the resetEvent() call: those hold pointers into
  // the owner's internal buffers, which clear() may reallocate/free, so they
  // must be re-obtained afterwards rather than reused.
  const MultiSourceFrame& getNormalizedFrame() const noexcept { return mNormalizedFrame; }
  MultiSourceFrameView getNormalizedFrameView() const noexcept { return mNormalizedFrame.getView(); }

  // Internal commit primitive for the owner-level load operation
  // (SurfaceTrackingScratch::loadNormalizedSource(), which stages both
  // this TimeFrame's normalized update and its own legacy backfill before
  // calling this) -- not general API, and never throws: `staged` must
  // already be the fully-built replacement frame. Swaps it in and clears
  // mCommonTracks/mTrackClusterIndices in the same commit, since a
  // CommonTrack/TrackClusterReference built against the previous normalized
  // frame is meaningless once this replaces it (see CommonTrack.h's own
  // lifetime doc).
  void commitNormalizedFrame(MultiSourceFrame&& staged) noexcept;

  // Commits a complete multi-source event after loader preflight. The
  // operation resets the prior event exactly once before installing all
  // normalized data and staged source workspaces.
  bool commitLoadedEvent(MultiSourceFrame&& staged, gsl::span<const ClusterSourceId> sources,
                         std::vector<std::unique_ptr<SurfaceTrackingScratch>>&& stagedWorkspaces) noexcept;

  // One generic event-state reset. It preserves static configuration and
  // allocator/capacity identity while clearing every source workspace and
  // event result. External owners call this once for a whole event.
  void resetEvent() noexcept;
  std::size_t getEventResetCount() const noexcept { return mEventResetCount; }

  SurfaceTrackingScratch& getWorkspace(ClusterSourceId source);
  const SurfaceTrackingScratch& getWorkspace(ClusterSourceId source) const;

  bool commitConfiguration(std::vector<SurfaceGraph>&& graphs,
                           std::vector<std::vector<TrackingParameters>>&& parameters,
                           std::vector<BindingSet>&& bindings,
                           std::vector<std::vector<TrackingWorkspaceCapacity>>&& capacities,
                           std::shared_ptr<BoundedMemoryResource> memoryPool);
  bool isConfigured() const noexcept { return mConfigurationValid; }
  std::size_t getNIterations() const noexcept { return mGraphs.size(); }
  std::size_t getNConfiguredSources() const noexcept { return mBindings.empty() ? 0 : mBindings.front().size(); }
  const std::vector<SurfaceGraph>& getGraphs() const noexcept { return mGraphs; }
  const std::vector<TrackingParameters>& getTrackingParameters() const noexcept;
  const std::vector<TrackingParameters>& getTrackingParameters(ClusterSourceId source) const noexcept;
  const SurfaceGraph& getGraph(std::size_t iteration) const { return mGraphs.at(iteration); }
  const TrackingParameters* getTrackingParameters(std::size_t iteration, ClusterSourceId source) const noexcept;
  const SurfacePlanBinding* getBinding(std::size_t iteration, ClusterSourceId source) const noexcept;
  const TrackingWorkspaceCapacity* getWorkspaceCapacity(std::size_t iteration, ClusterSourceId source) const noexcept;

  // Detector-neutral common-CA result storage (ITSMFTTracking/CommonTrack.h).
  // CommonTrack itself carries no NLayers dependency. Only meaningful
  // together with getNormalizedFrame()'s current content -- see
  // CommonTrack.h's own lifetime doc.
  auto& getCommonTracks() { return mCommonTracks; }
  const auto& getCommonTracks() const { return mCommonTracks; }
  // Flat array of TrackClusterReference (CommonTrack.h); a CommonTrack's
  // [firstClusterRef, clusterRefEnd) range is a half-open range of
  // *positions* into this array, in traversal order (inner to outer). Each
  // element pairs a SurfaceId with a SurfaceMeasurementIndex local to that
  // surface's own measurement array -- resolved via
  // getNormalizedFrame().getMeasurement(reference.surface, reference.index)
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

  // Detector-neutral common-CA result storage (Gate 4 CommonTrack
  // foundation; ITSMFTTracking/CommonTrack.h). Unpopulated by B3.1 -- no
  // production or test call site writes through these accessors from CA
  // seeds yet (see CommonTrack.h/SurfaceTrackingScratch.h for the shadow-
  // population contract a later slice adds).
  bounded_vector<CommonTrack> mCommonTracks;
  bounded_vector<TrackClusterReference> mTrackClusterIndices;

  // Normalized owner associated by the owner-level load operation
  // (SurfaceTrackingScratch::loadNormalizedSource()); host-only,
  // never GPU-managed or dictionary-serialized (see getNormalizedFrame()).
  // Does not itself hold pmr/bounded-vector allocations (MultiSourceFrame's
  // own members are plain std::vector<T> with the default allocator), so it
  // has no ordering dependency on mMemoryPool above.
  MultiSourceFrame mNormalizedFrame;

  bool mConfigurationValid = false;
  std::vector<SurfaceGraph> mGraphs;
  struct SourceParameters {
    ClusterSourceId source{};
    std::vector<TrackingParameters> values;
  };
  std::vector<SourceParameters> mTrackingParameters;
  std::vector<BindingSet> mBindings;
  std::vector<std::vector<TrackingWorkspaceCapacity>> mWorkspaceCapacities;
  struct WorkspaceDeleter {
    void operator()(SurfaceTrackingScratch* workspace) const noexcept;
  };
  struct WorkspaceEntry {
    ClusterSourceId source{};
    std::unique_ptr<SurfaceTrackingScratch, WorkspaceDeleter> workspace;
  };
  std::vector<WorkspaceEntry> mWorkspaces;
  std::size_t mEventResetCount{0};

  void clearEventData() noexcept;
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TIMEFRAME_H_ */
