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

bool TimeFrame::commitLoadedEvent(MultiSourceFrame&& staged, gsl::span<const ClusterSourceId> sources,
                                  std::vector<std::unique_ptr<SurfaceTrackingScratch>>&& stagedWorkspaces) noexcept
{
  if (!mConfigurationValid || sources.size() != mWorkspaces.size() || stagedWorkspaces.size() != mWorkspaces.size()) {
    return false;
  }
  for (std::size_t i = 0; i < mWorkspaces.size(); ++i) {
    const auto workspaceIt = std::find_if(mWorkspaces.begin(), mWorkspaces.end(), [source = sources[i]](const auto& entry) {
      return entry.source == source && entry.workspace != nullptr;
    });
    if (workspaceIt == mWorkspaces.end() || !stagedWorkspaces[i] ||
        stagedWorkspaces[i]->getNOwnedSurfaces() != workspaceIt->workspace->getNOwnedSurfaces() ||
        stagedWorkspaces[i]->getMemoryPool() != mMemoryPool || !workspaceIt->workspace->allocatorsMatch(*stagedWorkspaces[i])) {
      return false;
    }
  }

  resetEvent();
  mNormalizedFrame = std::move(staged);
  for (std::size_t i = 0; i < mWorkspaces.size(); ++i) {
    const auto workspaceIt = std::find_if(mWorkspaces.begin(), mWorkspaces.end(), [source = sources[i]](const auto& entry) {
      return entry.source == source;
    });
    workspaceIt->workspace->swapLoadedEvent(*stagedWorkspaces[i]);
  }
  return true;
}

bool TimeFrame::commitConfiguration(std::vector<SurfaceGraph>&& graphs,
                                    std::vector<std::vector<TrackingParameters>>&& parameters,
                                    std::vector<BindingSet>&& bindings,
                                    std::vector<std::vector<TrackingWorkspaceCapacity>>&& capacities,
                                    std::shared_ptr<BoundedMemoryResource> memoryPool)
{
  if (!memoryPool || graphs.empty() || graphs.size() != parameters.size() ||
      graphs.size() != bindings.size() || graphs.size() != capacities.size()) {
    return false;
  }
  for (std::size_t iteration = 0; iteration < graphs.size(); ++iteration) {
    if (bindings[iteration].empty() || bindings[iteration].size() != capacities[iteration].size() ||
        bindings[iteration].size() != parameters[iteration].size()) {
      return false;
    }
    for (std::size_t binding = 0; binding < bindings[iteration].size(); ++binding) {
      if (!bindings[iteration][binding] || capacities[iteration][binding].ownedSurfaces == 0) {
        return false;
      }
    }
  }

  std::vector<SourceParameters> sourceParameters;
  sourceParameters.reserve(bindings.front().size());
  for (std::size_t binding = 0; binding < bindings.front().size(); ++binding) {
    sourceParameters.push_back(SourceParameters{bindings.front()[binding]->getSource(), {}});
    sourceParameters.back().values.reserve(graphs.size());
  }
  for (std::size_t iteration = 0; iteration < graphs.size(); ++iteration) {
    for (std::size_t binding = 0; binding < bindings[iteration].size(); ++binding) {
      const auto source = bindings[iteration][binding]->getSource();
      const auto it = std::find_if(sourceParameters.begin(), sourceParameters.end(), [source](const auto& entry) {
        return entry.source == source;
      });
      if (it == sourceParameters.end()) {
        return false;
      }
      it->values.push_back(parameters[iteration][binding]);
    }
  }
  for (const auto& entry : sourceParameters) {
    if (entry.values.size() != graphs.size()) {
      return false;
    }
  }

  std::vector<WorkspaceEntry> workspaces;
  workspaces.reserve(sourceParameters.size());
  try {
    for (const auto& sourceParametersEntry : sourceParameters) {
      TrackingWorkspaceCapacity capacity{};
      for (std::size_t iteration = 0; iteration < bindings.size(); ++iteration) {
        for (std::size_t binding = 0; binding < bindings[iteration].size(); ++binding) {
          if (bindings[iteration][binding]->getSource() == sourceParametersEntry.source) {
            capacity.ownedSurfaces = std::max(capacity.ownedSurfaces, capacities[iteration][binding].ownedSurfaces);
            capacity.transitions = std::max(capacity.transitions, capacities[iteration][binding].transitions);
            capacity.cells = std::max(capacity.cells, capacities[iteration][binding].cells);
          }
        }
      }
      std::unique_ptr<SurfaceTrackingScratch, WorkspaceDeleter> workspace{new SurfaceTrackingScratch};
      workspace->setMemoryPool(memoryPool);
      workspace->adoptPlan(capacity.ownedSurfaces, capacity.transitions, capacity.cells);
      workspaces.push_back(WorkspaceEntry{sourceParametersEntry.source, std::move(workspace)});
    }
  } catch (const std::exception&) {
    return false;
  }

  // All validation is complete before any static or event state is replaced.
  // The prepared workspaces make the final replacement coherent as well.
  setMemoryPool(memoryPool);
  mGraphs = std::move(graphs);
  mTrackingParameters = std::move(sourceParameters);
  mBindings = std::move(bindings);
  mWorkspaceCapacities = std::move(capacities);
  mWorkspaces = std::move(workspaces);
  mConfigurationValid = true;
  return true;
}

SurfaceTrackingScratch& TimeFrame::getWorkspace(ClusterSourceId source)
{
  const auto it = std::find_if(mWorkspaces.begin(), mWorkspaces.end(), [source](const auto& entry) {
    return entry.source == source && entry.workspace != nullptr;
  });
  if (it == mWorkspaces.end()) {
    throw std::logic_error{"TimeFrame workspace is not configured for the requested source"};
  }
  return *it->workspace;
}

const SurfaceTrackingScratch& TimeFrame::getWorkspace(ClusterSourceId source) const
{
  const auto it = std::find_if(mWorkspaces.begin(), mWorkspaces.end(), [source](const auto& entry) {
    return entry.source == source && entry.workspace != nullptr;
  });
  if (it == mWorkspaces.end()) {
    throw std::logic_error{"TimeFrame workspace is not configured for the requested source"};
  }
  return *it->workspace;
}

void TimeFrame::resetEvent() noexcept
{
  for (auto& entry : mWorkspaces) {
    if (entry.workspace != nullptr) {
      entry.workspace->reset();
    }
  }
  clearEventData();
  ++mEventResetCount;
}

const std::vector<TrackingParameters>& TimeFrame::getTrackingParameters() const noexcept
{
  static const std::vector<TrackingParameters> empty;
  return mTrackingParameters.empty() ? empty : mTrackingParameters.front().values;
}

const std::vector<TrackingParameters>& TimeFrame::getTrackingParameters(ClusterSourceId source) const noexcept
{
  static const std::vector<TrackingParameters> empty;
  const auto it = std::find_if(mTrackingParameters.begin(), mTrackingParameters.end(), [source](const auto& entry) {
    return entry.source == source;
  });
  return it == mTrackingParameters.end() ? empty : it->values;
}

const TrackingParameters* TimeFrame::getTrackingParameters(std::size_t iteration, ClusterSourceId source) const noexcept
{
  const auto& parameters = getTrackingParameters(source);
  return iteration < parameters.size() ? &parameters[iteration] : nullptr;
}

const SurfacePlanBinding* TimeFrame::getBinding(std::size_t iteration, ClusterSourceId source) const noexcept
{
  if (!mConfigurationValid || iteration >= mBindings.size()) {
    return nullptr;
  }
  const auto& bindings = mBindings[iteration];
  const auto it = std::find_if(bindings.begin(), bindings.end(), [source](const auto& binding) {
    return binding && binding->getSource() == source;
  });
  return it == bindings.end() ? nullptr : it->get();
}

const TrackingWorkspaceCapacity* TimeFrame::getWorkspaceCapacity(std::size_t iteration, ClusterSourceId source) const noexcept
{
  if (!mConfigurationValid || iteration >= mBindings.size() || iteration >= mWorkspaceCapacities.size()) {
    return nullptr;
  }
  const auto& bindings = mBindings[iteration];
  const auto& capacities = mWorkspaceCapacities[iteration];
  for (std::size_t i = 0; i < bindings.size() && i < capacities.size(); ++i) {
    if (bindings[i] && bindings[i]->getSource() == source) {
      return &capacities[i];
    }
  }
  return nullptr;
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
