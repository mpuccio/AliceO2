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

void TimeFrame::rebuildMeasurementSpans()
{
  mMeasurementSpans.resize(mPerSurfaceMeasurements.size());
  for (std::size_t surface = 0; surface < mMeasurementSpans.size(); ++surface) {
    const auto& local = mPerSurfaceMeasurements[surface];
    const auto& global = mPerSurfaceGlobalMeasurements[surface];
    mMeasurementSpans[surface] = {global.empty() ? nullptr : global.data(), local.empty() ? nullptr : local.data(), static_cast<uint32_t>(local.size())};
  }
}

void TimeFrame::swapMeasurements(TimeFrame& other) noexcept
{
  mPerSurfaceGlobalMeasurements.swap(other.mPerSurfaceGlobalMeasurements);
  mPerSurfaceMeasurements.swap(other.mPerSurfaceMeasurements);
  mMeasurementSpans.swap(other.mMeasurementSpans);
  mROFIntervals.swap(other.mROFIntervals);
  mSourceROFOffsets.swap(other.mSourceROFOffsets);
  mLabelSources.swap(other.mLabelSources);
}

void TimeFrame::assignLoadedMeasurements(std::vector<std::vector<GlobalMeasurement>>&& globals,
                                         std::vector<std::vector<SurfaceMeasurement>>&& measurements,
                                         std::vector<ROFIntervalBC>&& intervals,
                                         std::vector<uint32_t>&& offsets,
                                         std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>&& labels)
{
  mPerSurfaceGlobalMeasurements = std::move(globals);
  mPerSurfaceMeasurements = std::move(measurements);
  mROFIntervals = std::move(intervals);
  mSourceROFOffsets = std::move(offsets);
  mLabelSources = std::move(labels);
  rebuildMeasurementSpans();
}

void TimeFrame::commitMeasurements(TimeFrame& staged) noexcept
{
  resetEvent();
  swapMeasurements(staged);
}

MeasurementView TimeFrame::getMeasurementView() const noexcept
{
  return {mMeasurementSpans.empty() ? nullptr : mMeasurementSpans.data(), static_cast<uint32_t>(mMeasurementSpans.size()),
          mROFIntervals.empty() ? nullptr : mROFIntervals.data(), mSourceROFOffsets.empty() ? nullptr : mSourceROFOffsets.data(),
          mSourceROFOffsets.empty() ? 0u : static_cast<uint32_t>(mSourceROFOffsets.size() - 1)};
}

gsl::span<const SurfaceMeasurement> TimeFrame::getSurfaceMeasurements(SurfaceId surface) const
{
  return surface.isValid() && surface.value() < mPerSurfaceMeasurements.size() ? gsl::make_span(mPerSurfaceMeasurements[surface.value()]) : gsl::span<const SurfaceMeasurement>{};
}

gsl::span<const GlobalMeasurement> TimeFrame::getGlobalMeasurements(SurfaceId surface) const
{
  return surface.isValid() && surface.value() < mPerSurfaceGlobalMeasurements.size() ? gsl::make_span(mPerSurfaceGlobalMeasurements[surface.value()]) : gsl::span<const GlobalMeasurement>{};
}

const GlobalMeasurement* TimeFrame::getGlobalMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept
{
  return getMeasurementView().getGlobalMeasurement(surface, index);
}

const SurfaceMeasurement* TimeFrame::getSurfaceMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept
{
  return getMeasurementView().getSurfaceMeasurement(surface, index);
}

gsl::span<const ROFIntervalBC> TimeFrame::getSourceIntervals(ClusterSourceId source) const
{
  if (!source.isValid() || source.value() + 1 >= mSourceROFOffsets.size()) {
    return {};
  }
  return gsl::make_span(mROFIntervals).subspan(mSourceROFOffsets[source.value()], mSourceROFOffsets[source.value() + 1] - mSourceROFOffsets[source.value()]);
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
  for (const auto& measurements : mPerSurfaceMeasurements) {
    total += measurements.size();
  }
  return total;
}

bool TimeFrame::commitLoadedEvent(TimeFrame& staged,
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
  swapMeasurements(staged);
  mWorkspace->swapLoadedEvent(*stagedWorkspace);
  return true;
}

bool TimeFrame::commitConfiguration(std::vector<SurfaceGraph>&& graphs,
                                    std::vector<TrackingParameters>&& parameters,
                                    std::vector<TrackingWorkspaceCapacity>&& capacities,
                                    std::shared_ptr<BoundedMemoryResource> memoryPool)
{
  if (!memoryPool || graphs.empty() || graphs.size() != parameters.size() ||
      graphs.size() != capacities.size()) {
    return false;
  }
  for (std::size_t iteration = 0; iteration < graphs.size(); ++iteration) {
    if (!graphs[iteration].valid() || capacities[iteration].ownedSurfaces == 0) {
      return false;
    }
  }

  TrackingWorkspaceCapacity capacity{};
  for (const auto& iterationCapacity : capacities) {
    capacity.ownedSurfaces = std::max(capacity.ownedSurfaces, iterationCapacity.ownedSurfaces);
    capacity.edges = std::max(capacity.edges, iterationCapacity.edges);
    capacity.cells = std::max(capacity.cells, iterationCapacity.cells);
  }
  std::unique_ptr<SurfaceTrackingScratch, WorkspaceDeleter> workspace;
  try {
    workspace.reset(new SurfaceTrackingScratch);
    workspace->setMemoryPool(memoryPool);
    workspace->adoptPlan(capacity.ownedSurfaces, capacity.edges, capacity.cells);
    workspace->configureTraversalWorkspaces(parameters.size());
  } catch (const std::exception&) {
    return false;
  }

  // Replace state only after validation and workspace preparation succeed.
  setMemoryPool(memoryPool);
  mGraphs = std::move(graphs);
  mTrackingParameters = std::move(parameters);
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

void TimeFrame::clearEventData() noexcept
{
  deepVectorClear(mPrimaryVertices);
  deepVectorClear(mPrimaryVerticesLabels);
  // Common tracks and their cluster references are valid only for the current
  // normalized event, so clear both together.
  deepVectorClear(mGenericTracks);
  deepVectorClear(mTrackClusterIndices);
  clearMeasurements();
}

void TimeFrame::clearMeasurements() noexcept
{
  mPerSurfaceGlobalMeasurements.clear();
  mPerSurfaceMeasurements.clear();
  mMeasurementSpans.clear();
  mROFIntervals.clear();
  mSourceROFOffsets.clear();
  mLabelSources.clear();
}

} // namespace o2::itsmft::tracking
