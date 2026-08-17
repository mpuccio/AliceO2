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
#include "ITSMFTTracking/detail/TimeFrameScratch.h"

#include <algorithm>
#include <stdexcept>

namespace o2::itsmft::tracking
{

using o2::its::deepVectorClear;

TimeFrame::~TimeFrame() = default;

void TimeFrame::WorkspaceDeleter::operator()(TimeFrameScratch* workspace) const noexcept
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

void TimeFrame::swapMeasurements(TimeFrame& other) noexcept
{
  mLayerGlobalMeasurements.swap(other.mLayerGlobalMeasurements);
  mLayerSurfaceMeasurements.swap(other.mLayerSurfaceMeasurements);
  mLabelSources.swap(other.mLabelSources);
}

void TimeFrame::assignLoadedMeasurements(std::vector<std::vector<GlobalMeasurement>>&& globals,
                                         std::vector<std::vector<SurfaceMeasurement>>&& measurements,
                                         std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>&& labels)
{
  mLayerGlobalMeasurements = std::move(globals);
  mLayerSurfaceMeasurements = std::move(measurements);
  mLabelSources = std::move(labels);
}

void TimeFrame::commitMeasurements(TimeFrame& staged) noexcept
{
  resetTimeFrame();
  swapMeasurements(staged);
}

gsl::span<const SurfaceMeasurement> TimeFrame::getSurfaceMeasurements(LayerId surface) const
{
  return surface.isValid() && surface.value() < mLayerSurfaceMeasurements.size() ? gsl::make_span(mLayerSurfaceMeasurements[surface.value()]) : gsl::span<const SurfaceMeasurement>{};
}

gsl::span<const GlobalMeasurement> TimeFrame::getGlobalMeasurements(LayerId surface) const
{
  return surface.isValid() && surface.value() < mLayerGlobalMeasurements.size() ? gsl::make_span(mLayerGlobalMeasurements[surface.value()]) : gsl::span<const GlobalMeasurement>{};
}

const GlobalMeasurement* TimeFrame::getGlobalMeasurement(LayerId surface, MeasurementIndex index) const noexcept
{
  if (!surface.isValid() || surface.value() >= mLayerGlobalMeasurements.size() || !index.isValid()) {
    return nullptr;
  }
  const auto& measurements = mLayerGlobalMeasurements[surface.value()];
  return index.value() < measurements.size() ? &measurements[index.value()] : nullptr;
}

const SurfaceMeasurement* TimeFrame::getSurfaceMeasurement(LayerId surface, MeasurementIndex index) const noexcept
{
  if (!surface.isValid() || surface.value() >= mLayerSurfaceMeasurements.size() || !index.isValid()) {
    return nullptr;
  }
  const auto& measurements = mLayerSurfaceMeasurements[surface.value()];
  return index.value() < measurements.size() ? &measurements[index.value()] : nullptr;
}

gsl::span<const o2::MCCompLabel> TimeFrame::getLabels(ClusterRef cluster) const
{
  if (!cluster.source.isValid() || cluster.source.value() >= mLabelSources.size() || mLabelSources[cluster.source.value()] == nullptr) {
    return {};
  }
  return mLabelSources[cluster.source.value()]->getLabels(cluster.index);
}

std::size_t TimeFrame::getTotalMeasurements() const noexcept
{
  std::size_t total = 0;
  for (const auto& measurements : mLayerSurfaceMeasurements) {
    total += measurements.size();
  }
  return total;
}

bool TimeFrame::commitLoadedEvent(TimeFrame& staged,
                                  std::unique_ptr<TimeFrameScratch>&& stagedWorkspace) noexcept
{
  if (!mConfigurationValid || !mWorkspace || !stagedWorkspace) {
    return false;
  }
  if (stagedWorkspace->getNOwnedSurfaces() != mWorkspace->getNOwnedSurfaces() ||
      stagedWorkspace->getMemoryPool() != mMemoryPool || !mWorkspace->allocatorsMatch(*stagedWorkspace)) {
    return false;
  }

  resetTimeFrame();
  swapMeasurements(staged);
  mWorkspace->swapLoadedEvent(*stagedWorkspace);
  return true;
}

bool TimeFrame::commitConfiguration(std::vector<SurfaceLayout>&& layouts,
                                    std::vector<TrackingParameters>&& parameters,
                                    std::vector<TrackingWorkspaceCapacity>&& capacities,
                                    std::shared_ptr<BoundedMemoryResource> memoryPool)
{
  if (!memoryPool || layouts.empty() || layouts.size() != parameters.size() ||
      layouts.size() != capacities.size()) {
    return false;
  }
  for (std::size_t iteration = 0; iteration < layouts.size(); ++iteration) {
    if (!layouts[iteration].valid() || capacities[iteration].ownedSurfaces == 0) {
      return false;
    }
  }

  TrackingWorkspaceCapacity capacity{};
  for (const auto& iterationCapacity : capacities) {
    capacity.ownedSurfaces = std::max(capacity.ownedSurfaces, iterationCapacity.ownedSurfaces);
    capacity.edges = std::max(capacity.edges, iterationCapacity.edges);
    capacity.cells = std::max(capacity.cells, iterationCapacity.cells);
  }
  std::unique_ptr<TimeFrameScratch, WorkspaceDeleter> workspace;
  try {
    workspace.reset(new TimeFrameScratch);
    workspace->setMemoryPool(memoryPool);
    workspace->adoptPlan(capacity.ownedSurfaces, capacity.edges, capacity.cells);
    workspace->configureTraversalWorkspaces(parameters.size());
  } catch (const std::exception&) {
    return false;
  }

  setMemoryPool(memoryPool);
  mLayouts = std::move(layouts);
  mTrackingParameters = std::move(parameters);
  mWorkspaceCapacities = std::move(capacities);
  mWorkspace = std::move(workspace);
  mConfigurationValid = true;
  return true;
}

TimeFrameScratch& TimeFrame::getWorkspace()
{
  if (!mWorkspace) {
    throw std::logic_error{"TimeFrame workspace is not configured"};
  }
  return *mWorkspace;
}

const TimeFrameScratch& TimeFrame::getWorkspace() const
{
  if (!mWorkspace) {
    throw std::logic_error{"TimeFrame workspace is not configured"};
  }
  return *mWorkspace;
}

void TimeFrame::resetTimeFrame() noexcept
{
  if (mWorkspace) {
    mWorkspace->reset();
  }
  deepVectorClear(mPrimaryVertices);
  deepVectorClear(mPrimaryVerticesLabels);
  // Common tracks and their cluster references are valid only for the current
  // normalized event, so clear both together.
  deepVectorClear(mGenericTracks);
  deepVectorClear(mTrackClusterIndices);
  mLayerGlobalMeasurements.clear();
  mLayerSurfaceMeasurements.clear();
  mLabelSources.clear();
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
  initVector(mGenericTracks);
  initVector(mTrackClusterIndices);
}

} // namespace o2::itsmft::tracking
