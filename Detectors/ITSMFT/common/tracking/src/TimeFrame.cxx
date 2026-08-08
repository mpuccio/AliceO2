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
/// \file TimeFrame.cxx
/// \brief
///

#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"

#include <algorithm>
#include <stdexcept>

namespace o2::itsmft::tracking
{

using o2::its::deepVectorClear;

TimeFrame::~TimeFrame() = default;

void TimeFrame::WorkspaceDeleter::operator()(SurfaceTrackingScratch* workspace) const noexcept
{
  delete workspace;
}

void TimeFrame::addPrimaryVertex(const Vertex& vert)
{
  mPrimaryVertices.emplace_back(vert);
  if (!isBeamPositionOverridden) {
    const float w = vert.getNContributors();
    mBeamPos[0] = (mBeamPos[0] * mBeamPosWeight + vert.getX() * w) / (mBeamPosWeight + w);
    mBeamPos[1] = (mBeamPos[1] * mBeamPosWeight + vert.getY() * w) / (mBeamPosWeight + w);
    mBeamPosWeight += w;
  }
}

void TimeFrame::resetBeamXY(const float x, const float y, const float w)
{
  mBeamPos[0] = x;
  mBeamPos[1] = y;
  mBeamPosWeight = w;
}

void TimeFrame::commitNormalizedFrame(MultiSourceFrame&& staged) noexcept
{
  resetEvent();
  mNormalizedFrame = std::move(staged);
}

bool TimeFrame::commitLoadedEvent(MultiSourceFrame&& staged,
                                  std::unique_ptr<SurfaceTrackingScratch>&& stagedWorkspace) noexcept
{
  if (!mConfigurationValid || !mWorkspace || !stagedWorkspace) {
    return false;
  }
  if (stagedWorkspace->getNOwnedSurfaces() != mWorkspace->getNOwnedSurfaces() ||
      stagedWorkspace->getMemoryPool() != mMemoryPool || !mWorkspace->allocatorsMatch(*stagedWorkspace)) {
    return false;
  }

  resetEvent();
  mNormalizedFrame = std::move(staged);
  mWorkspace->swapLoadedEvent(*stagedWorkspace);
  return true;
}

bool TimeFrame::commitConfiguration(std::vector<SurfaceGraph>&& graphs,
                                    std::vector<TrackingParameters>&& parameters,
                                    std::vector<std::unique_ptr<SurfacePlanBinding>>&& bindings,
                                    std::vector<TrackingWorkspaceCapacity>&& capacities,
                                    std::shared_ptr<BoundedMemoryResource> memoryPool)
{
  if (!memoryPool || graphs.empty() || graphs.size() != parameters.size() ||
      graphs.size() != bindings.size() || graphs.size() != capacities.size()) {
    return false;
  }
  for (std::size_t iteration = 0; iteration < graphs.size(); ++iteration) {
    if (!bindings[iteration] || capacities[iteration].ownedSurfaces == 0) {
      return false;
    }
  }

  TrackingWorkspaceCapacity capacity{};
  for (const auto& iterationCapacity : capacities) {
    capacity.ownedSurfaces = std::max(capacity.ownedSurfaces, iterationCapacity.ownedSurfaces);
    capacity.transitions = std::max(capacity.transitions, iterationCapacity.transitions);
    capacity.cells = std::max(capacity.cells, iterationCapacity.cells);
  }
  std::unique_ptr<SurfaceTrackingScratch, WorkspaceDeleter> workspace;
  try {
    workspace.reset(new SurfaceTrackingScratch);
    workspace->setMemoryPool(memoryPool);
    workspace->adoptPlan(capacity.ownedSurfaces, capacity.transitions, capacity.cells);
  } catch (const std::exception&) {
    return false;
  }

  // All validation is complete before any static or event state is replaced.
  // The prepared workspaces make the final replacement coherent as well.
  setMemoryPool(memoryPool);
  mGraphs = std::move(graphs);
  mTrackingParameters = std::move(parameters);
  mBindings = std::move(bindings);
  mWorkspaceCapacities = std::move(capacities);
  mWorkspace = std::move(workspace);
  mConfigurationValid = true;
  return true;
}

SurfaceTrackingScratch& TimeFrame::getWorkspace()
{
  if (!mWorkspace) {
    throw std::logic_error{"TimeFrame workspace is not configured"};
  }
  return *mWorkspace;
}

const SurfaceTrackingScratch& TimeFrame::getWorkspace() const
{
  if (!mWorkspace) {
    throw std::logic_error{"TimeFrame workspace is not configured"};
  }
  return *mWorkspace;
}

void TimeFrame::resetEvent() noexcept
{
  if (mWorkspace) {
    mWorkspace->reset();
  }
  clearEventData();
  ++mEventResetCount;
}

const std::vector<TrackingParameters>& TimeFrame::getTrackingParameters() const noexcept
{
  return mTrackingParameters;
}

const TrackingParameters* TimeFrame::getTrackingParameters(std::size_t iteration) const noexcept
{
  return iteration < mTrackingParameters.size() ? &mTrackingParameters[iteration] : nullptr;
}

const SurfacePlanBinding* TimeFrame::getBinding(std::size_t iteration) const noexcept
{
  if (!mConfigurationValid || iteration >= mBindings.size()) {
    return nullptr;
  }
  return mBindings[iteration].get();
}

const TrackingWorkspaceCapacity* TimeFrame::getWorkspaceCapacity(std::size_t iteration) const noexcept
{
  if (!mConfigurationValid || iteration >= mWorkspaceCapacities.size()) {
    return nullptr;
  }
  return &mWorkspaceCapacities[iteration];
}

void TimeFrame::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
{
  mMemoryPool = pool;

  auto initVector = [&]<typename T>(bounded_vector<T>& vec) {
    deepVectorClear(vec, mMemoryPool.get());
  };

  initVector(mPrimaryVertices);
  initVector(mPrimaryVerticesLabels);
  initVector(mCommonTracks);
  initVector(mTrackClusterIndices);
}

void TimeFrame::clearEventData() noexcept
{
  deepVectorClear(mPrimaryVertices);
  deepVectorClear(mPrimaryVerticesLabels);
  // Common-CA result storage is per-event data and is invalidated together: a
  // CommonTrack's cluster-reference range
  // and every SurfaceMeasurementIndex it reaches are only meaningful
  // alongside the normalized frame current when it was built (CommonTrack.h),
  // so both must be cleared here, unconditionally.
  deepVectorClear(mCommonTracks);
  deepVectorClear(mTrackClusterIndices);
  // Event-owned normalized data, cleared unconditionally (MultiSourceFrame is
  // host-only and never framework/GPU-managed). This clears mNormalizedFrame
  // in place: a getNormalizedFrame() reference obtained before this call
  // still refers to that same live member and remains safe to use, now
  // observing its cleared (empty) state. Any MultiSourceFrameView or
  // gsl::span obtained before this call (from getNormalizedFrameView(),
  // getSurfaceMeasurements(), getSourceIntervals(), getLabels()) is
  // invalidated by it, since clear() may reallocate/free the buffers those
  // point into.
  mNormalizedFrame.clear();
}

} // namespace o2::itsmft::tracking
