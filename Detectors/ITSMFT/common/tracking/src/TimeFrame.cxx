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
#include <limits>
#include <new>
#include <numeric>
#include <stdexcept>
#include <type_traits>

#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITStracking/MathUtils.h"

namespace o2::itsmft::tracking
{

using o2::its::deepVectorClear;

namespace
{

template <typename T>
void replaceBoundedVector(bounded_vector<T>& destination, bounded_vector<T>& source) noexcept
{
  static_assert(std::is_nothrow_move_constructible_v<bounded_vector<T>>);
  destination.~bounded_vector<T>();
  new (&destination) bounded_vector<T>(std::move(source));
}

} // namespace

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
  mLayerUsedClusters.swap(other.mLayerUsedClusters);
  mLayerClusterLabels.swap(other.mLayerClusterLabels);
  std::swap(mHasMCInformation, other.mHasMCInformation);
  mROFramesClusters.swap(other.mROFramesClusters);
  std::swap(mROFViews, other.mROFViews);
  mROFViewsBySurface.swap(other.mROFViewsBySurface);
  mROFLocalLayerBySurface.swap(other.mROFLocalLayerBySurface);
  std::swap(mUseUPC, other.mUseUPC);
}

void TimeFrame::assignLoadedEventNavigation(std::vector<std::vector<int>>&& rofBoundaries,
                                            RuntimeROFViews defaultViews,
                                            std::vector<RuntimeROFViews>&& viewsBySurface,
                                            std::vector<uint16_t>&& localLayerBySurface)
{
  mROFramesClusters = std::move(rofBoundaries);
  mROFViews = defaultViews;
  mROFViewsBySurface = std::move(viewsBySurface);
  mROFLocalLayerBySurface = std::move(localLayerBySurface);
  mUseUPC = false;
}

void TimeFrame::assignLoadedMeasurements(std::vector<std::vector<GlobalMeasurement>>&& globals,
                                         std::vector<std::vector<SurfaceMeasurement>>&& measurements,
                                         std::vector<o2::dataformats::MCTruthContainer<o2::MCCompLabel>>&& labels,
                                         bool hasMCInformation)
{
  mLayerGlobalMeasurements = std::move(globals);
  mLayerSurfaceMeasurements = std::move(measurements);
  mLayerUsedClusters.clear();
  mLayerUsedClusters.reserve(mLayerSurfaceMeasurements.size());
  for (const auto& layerMeasurements : mLayerSurfaceMeasurements) {
    mLayerUsedClusters.emplace_back(layerMeasurements.size(), uint8_t{0});
  }
  mLayerClusterLabels = std::move(labels);
  mHasMCInformation = hasMCInformation;
}

void TimeFrame::commitMeasurements(TimeFrame& staged) noexcept
{
  resetTimeFrame();
  swapMeasurements(staged);
}

gsl::span<const GlobalMeasurement> TimeFrame::getGlobalMeasurements(LayerId surface) const
{
  return surface.isValid() && surface.value() < mLayerGlobalMeasurements.size() ? gsl::make_span(mLayerGlobalMeasurements[surface.value()]) : gsl::span<const GlobalMeasurement>{};
}

gsl::span<GlobalMeasurement> TimeFrame::getGlobalMeasurements(LayerId surface)
{
  return surface.isValid() && surface.value() < mLayerGlobalMeasurements.size() ? gsl::make_span(mLayerGlobalMeasurements[surface.value()]) : gsl::span<GlobalMeasurement>{};
}

const SurfaceMeasurement* TimeFrame::getSurfaceMeasurement(LayerId layer, uint32_t clusterId) const noexcept
{
  if (!layer.isValid() || layer.value() >= mLayerSurfaceMeasurements.size()) {
    return nullptr;
  }
  const auto& measurements = mLayerSurfaceMeasurements[layer.value()];
  return clusterId < measurements.size() ? &measurements[clusterId] : nullptr;
}

gsl::span<const o2::MCCompLabel> TimeFrame::getLabels(LayerId layer, uint32_t clusterId) const
{
  if (!layer.isValid() || layer.value() >= mLayerClusterLabels.size()) {
    return {};
  }
  return mLayerClusterLabels[layer.value()].getLabels(clusterId);
}

std::size_t TimeFrame::getTotalMeasurements() const noexcept
{
  std::size_t total = 0;
  for (const auto& measurements : mLayerGlobalMeasurements) {
    total += measurements.size();
  }
  return total;
}

gsl::span<GlobalMeasurement> TimeFrame::getClustersOnLayer(int rofId, int layer)
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const int first = mROFramesClusters[layer][rofId];
  return {mLayerGlobalMeasurements[layer].data() + first,
          static_cast<gsl::span<GlobalMeasurement>::size_type>(mROFramesClusters[layer][rofId + 1] - first)};
}

gsl::span<const GlobalMeasurement> TimeFrame::getClustersOnLayer(int rofId, int layer) const
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const int first = mROFramesClusters[layer][rofId];
  return {mLayerGlobalMeasurements[layer].data() + first,
          static_cast<gsl::span<const GlobalMeasurement>::size_type>(mROFramesClusters[layer][rofId + 1] - first)};
}

gsl::span<const GlobalMeasurement> TimeFrame::getClustersPerROFrange(int rofMin, int range, int layer) const
{
  if (rofMin < 0 || rofMin >= getNrof(layer)) {
    return {};
  }
  const int first = mROFramesClusters[layer][rofMin];
  const int last = mROFramesClusters[layer][o2::gpu::CAMath::Min(rofMin + range, getNrof(layer))];
  return {mLayerGlobalMeasurements[layer].data() + first, static_cast<gsl::span<const GlobalMeasurement>::size_type>(last - first)};
}

gsl::span<const int> TimeFrame::getROFramesClustersPerROFrange(int rofMin, int range, int layer) const
{
  const int checkedRange = o2::gpu::CAMath::Min(range, getNrof(layer) - rofMin);
  return {mROFramesClusters[layer].data() + rofMin, static_cast<gsl::span<const int>::size_type>(checkedRange)};
}

gsl::span<const int> TimeFrame::getROFrameClusters(int layer) const
{
  return gsl::make_span(mROFramesClusters[layer]);
}

gsl::span<int> TimeFrame::getIndexTable(int rofId, int layer)
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const int tableSize = mIndexTableUtils[layer].getNrowBins() * mIndexTableUtils[layer].getNcolBins() + 1;
  return {mIndexTables[layer].data() + rofId * tableSize, static_cast<gsl::span<int>::size_type>(tableSize)};
}

int TimeFrame::getClusterROF(int layer, int cluster) const
{
  return static_cast<int>(std::lower_bound(mROFramesClusters[layer].begin(), mROFramesClusters[layer].end(), cluster + 1) -
                          mROFramesClusters[layer].begin() - 1);
}

int TimeFrame::getTotalClustersPerROFrange(int rofMin, int range, int layer) const
{
  const int last = o2::gpu::CAMath::Min(rofMin + range, getNrof(layer));
  return mROFramesClusters[layer][last] - mROFramesClusters[layer][rofMin];
}

gsl::span<unsigned char> TimeFrame::getUsedClusters(int layer)
{
  return layer >= 0 && static_cast<std::size_t>(layer) < mLayerUsedClusters.size() ? gsl::make_span(mLayerUsedClusters[layer]) : gsl::span<unsigned char>{};
}

bool TimeFrame::isClusterUsed(int layer, uint32_t clusterId) const
{
  return layer >= 0 && static_cast<std::size_t>(layer) < mLayerUsedClusters.size() && clusterId < mLayerUsedClusters[layer].size() && mLayerUsedClusters[layer][clusterId] != 0;
}

void TimeFrame::markUsedCluster(int layer, uint32_t clusterId)
{
  if (layer >= 0 && static_cast<std::size_t>(layer) < mLayerUsedClusters.size() && clusterId < mLayerUsedClusters[layer].size()) {
    mLayerUsedClusters[layer][clusterId] = 1;
  }
}

std::size_t TimeFrame::getNumberOfClusters() const
{
  return std::accumulate(mLayerGlobalMeasurements.begin(), mLayerGlobalMeasurements.end(), std::size_t{0},
                         [](std::size_t total, const auto& layer) { return total + layer.size(); });
}

std::size_t TimeFrame::getNumberOfUsedClusters() const
{
  return std::accumulate(mLayerUsedClusters.begin(), mLayerUsedClusters.end(), std::size_t{0}, [](std::size_t total, const auto& layer) {
    return total + static_cast<std::size_t>(std::count(layer.begin(), layer.end(), uint8_t{1}));
  });
}

void TimeFrame::setROFViews(RuntimeROFViews views) noexcept
{
  mROFViews = views;
  mROFViewsBySurface.assign(mWorkspaceCapacities.empty() ? 0 : mWorkspaceCapacities.front().ownedSurfaces, views);
  mROFLocalLayerBySurface.resize(mROFViewsBySurface.size());
  std::iota(mROFLocalLayerBySurface.begin(), mROFLocalLayerBySurface.end(), uint16_t{0});
  mUseUPC = false;
}

const RuntimeROFTableEntry& TimeFrame::getROFOverlap(int fromLayer, int toLayer, int rof) const noexcept
{
  return getROFViews(fromLayer).overlap.getOverlap(getROFLocalLayer(fromLayer), getROFLocalLayer(toLayer), rof);
}

bool TimeFrame::isROFEnabled(int layer, int rof) const noexcept
{
  const auto& views = getROFViews(layer);
  return (mUseUPC ? views.upcMask : views.mask).isROFEnabled(getROFLocalLayer(layer), rof);
}

bool TimeFrame::isVertexCompatible(int layer, int rof, const Vertex& vertex) const noexcept
{
  return getROFViews(layer).vertexLookup.isVertexCompatible(getROFLocalLayer(layer), rof, vertex);
}

o2::its::TimeEstBC TimeFrame::getROFTimeStamp(int fromLayer, int fromROF, int toLayer, int toROF) const noexcept
{
  return getROFViews(fromLayer).overlap.getTimeStamp(getROFLocalLayer(fromLayer), fromROF,
                                                     getROFLocalLayer(toLayer), toROF);
}

int TimeFrame::getMaxVerticesPerROF() const noexcept
{
  if (mROFViewsBySurface.empty()) {
    return mROFViews.vertexLookup.getMaxVerticesPerROF();
  }
  int result = 0;
  for (const auto& views : mROFViewsBySurface) {
    result = std::max(result, views.vertexLookup.getMaxVerticesPerROF());
  }
  return result;
}

gsl::span<const Vertex> TimeFrame::getPrimaryVertices(int layer, int rofId) const
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const auto& entry = getROFViews(layer).vertexLookup.getVertices(getROFLocalLayer(layer), rofId);
  return {mPrimaryVertices.data() + entry.getFirstEntry(),
          static_cast<gsl::span<const Vertex>::size_type>(entry.getEntries())};
}

bool TimeFrame::hasMCinformation() const noexcept
{
  return mHasMCInformation;
}

gsl::span<const MCCompLabel> TimeFrame::getClusterLabels(int layer, int cluster) const
{
  if (layer < 0 || static_cast<std::size_t>(layer) >= mLayerGlobalMeasurements.size() || cluster < 0 || static_cast<std::size_t>(cluster) >= mLayerGlobalMeasurements[layer].size()) {
    return {};
  }
  return getLabels(LayerId{static_cast<uint16_t>(layer)}, mLayerGlobalMeasurements[layer][cluster].clusterId);
}

bool TimeFrame::commitLoadedEvent(TimeFrame& staged) noexcept
{
  if (!mConfigurationValid || !mWorkspace || staged.mROFramesClusters.size() != mROFramesClusters.size()) {
    return false;
  }

  resetTimeFrame();
  swapMeasurements(staged);
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

  TimeFrame staged;
  try {
    staged.mWorkspace.reset(new TimeFrameScratch);
    staged.mWorkspace->setMemoryPool(memoryPool);
    staged.mWorkspace->adoptPlan(capacity.ownedSurfaces, capacity.edges, capacity.cells);
    staged.mWorkspace->configureTraversalWorkspaces(parameters.size());

    staged.mROFramesClusters = mROFramesClusters;
    staged.setMemoryPool(memoryPool);
    staged.configureEventStorage(capacity.ownedSurfaces);
  } catch (const std::bad_alloc&) {
    return false;
  }

  staged.mLayouts = std::move(layouts);
  staged.mTrackingParameters = std::move(parameters);
  staged.mWorkspaceCapacities = std::move(capacities);
  staged.mConfigurationValid = true;
  publishConfiguration(staged);
  return true;
}

void TimeFrame::publishConfiguration(TimeFrame& staged) noexcept
{
  static_assert(std::is_nothrow_move_assignable_v<decltype(mROFramesClusters)>);
  static_assert(std::is_nothrow_move_assignable_v<decltype(mIndexTables)>);
  static_assert(std::is_nothrow_move_assignable_v<decltype(mIndexTableUtils)>);
  static_assert(std::is_nothrow_move_assignable_v<decltype(mLayouts)>);
  static_assert(std::is_nothrow_move_assignable_v<decltype(mTrackingParameters)>);
  static_assert(std::is_nothrow_move_assignable_v<decltype(mWorkspaceCapacities)>);

  // Destroy every object backed by the old pool while the frame still owns
  // that pool. PMR vectors cannot be swapped across unequal resources, but
  // their allocator-preserving move construction is noexcept.
  mWorkspace.reset();
  replaceBoundedVector(mPrimaryVertices, staged.mPrimaryVertices);
  replaceBoundedVector(mPrimaryVerticesLabels, staged.mPrimaryVerticesLabels);
  replaceBoundedVector(mGenericTracks, staged.mGenericTracks);
  replaceBoundedVector(mTrackClusterIndices, staged.mTrackClusterIndices);

  mROFramesClusters = std::move(staged.mROFramesClusters);
  mIndexTables = std::move(staged.mIndexTables);
  mIndexTableUtils = std::move(staged.mIndexTableUtils);
  mMinR = std::move(staged.mMinR);
  mMaxR = std::move(staged.mMaxR);
  mMinZ = std::move(staged.mMinZ);
  mMaxZ = std::move(staged.mMaxZ);
  mLayouts = std::move(staged.mLayouts);
  mTrackingParameters = std::move(staged.mTrackingParameters);
  mWorkspaceCapacities = std::move(staged.mWorkspaceCapacities);
  mWorkspace = std::move(staged.mWorkspace);
  mConfigurationValid = staged.mConfigurationValid;

  // All objects using the old resource are gone. Publish the new owner last;
  // staged keeps the old resource alive until its moved-from members die.
  mMemoryPool.swap(staged.mMemoryPool);
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
  mLayerUsedClusters.clear();
  mLayerClusterLabels.clear();
  mHasMCInformation = false;
  mROFViews = {};
  mROFViewsBySurface.clear();
  mROFLocalLayerBySurface.clear();
  mUseUPC = false;
  for (auto& boundaries : mROFramesClusters) {
    boundaries.clear();
  }
  deepVectorClear(mIndexTables);
  std::fill(mMinR.begin(), mMinR.end(), std::numeric_limits<float>::max());
  std::fill(mMaxR.begin(), mMaxR.end(), std::numeric_limits<float>::lowest());
  std::fill(mMinZ.begin(), mMinZ.end(), std::numeric_limits<float>::max());
  std::fill(mMaxZ.begin(), mMaxZ.end(), std::numeric_limits<float>::lowest());
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
  for (auto& table : mIndexTables) {
    initVector(table);
  }
}

void TimeFrame::configureEventStorage(std::size_t nOwnedSurfaces)
{
  mROFramesClusters.resize(nOwnedSurfaces);
  o2::its::clearResizeBoundedVector(mIndexTables, nOwnedSurfaces, mMemoryPool.get());
  mIndexTableUtils.assign(nOwnedSurfaces, IndexTableUtilsCore{});
  mMinR.assign(nOwnedSurfaces, std::numeric_limits<float>::max());
  mMaxR.assign(nOwnedSurfaces, std::numeric_limits<float>::lowest());
  mMinZ.assign(nOwnedSurfaces, std::numeric_limits<float>::max());
  mMaxZ.assign(nOwnedSurfaces, std::numeric_limits<float>::lowest());
}

void TimeFrame::prepareClusters(int maxLayers)
{
  struct SortingHelper {
    int bin;
    int indexWithinBin;
    int measurementIndex;
  };

  const int stopLayer = std::min(maxLayers, static_cast<int>(mLayerGlobalMeasurements.size()));
  for (int layer = 0; layer < stopLayer; ++layer) {
    const auto& utils = mIndexTableUtils[layer];
    const int colBinsCount = utils.getNcolBins();
    std::size_t numBins = 0;
    if (!checkedIndexTableSizeProduct(static_cast<std::size_t>(utils.getNrowBins()),
                                      static_cast<std::size_t>(colBinsCount), numBins) ||
        numBins == std::numeric_limits<std::size_t>::max()) {
      throw std::bad_alloc{};
    }
    const std::size_t stride = numBins + 1;
    bounded_vector<SortingHelper> helpers(mMemoryPool.get());
    bounded_vector<GlobalMeasurement> sortedMeasurements(mMemoryPool.get());
    bounded_vector<int> counts(numBins, 0, mMemoryPool.get());
    bounded_vector<int> offsets(numBins, 0, mMemoryPool.get());

    for (int rof = 0; rof < getNrof(layer); ++rof) {
      if (!isROFEnabled(layer, rof)) {
        continue;
      }
      const int first = mROFramesClusters[layer][rof];
      const int last = mROFramesClusters[layer][rof + 1];
      const int count = last - first;
      auto* tableBase = mIndexTables[layer].data() + rof * stride;
      helpers.resize(count);
      sortedMeasurements.resize(count);
      const bool usePhiRBinning = utils.getCoordType() == o2::itsmft::IndexTableCoordType::PhiR;

      for (int local = 0; local < count; ++local) {
        const int measurementIndex = first + local;
        const auto& measurement = mLayerGlobalMeasurements[layer][measurementIndex];
        auto& helper = helpers[local];
        int colBin = utils.getColBinIndex(layer, usePhiRBinning ? measurement.radius : measurement.z);
        if (colBin < 0 || colBin >= colBinsCount) {
          colBin = std::clamp(colBin, 0, colBinsCount - 1);
        }
        helper.bin = utils.getBinIndex(colBin, utils.getRowBinIndex(measurement.phi));
        helper.indexWithinBin = counts[helper.bin]++;
        helper.measurementIndex = measurementIndex;
        mMinR[layer] = o2::gpu::GPUCommonMath::Min(measurement.radius, mMinR[layer]);
        mMaxR[layer] = o2::gpu::GPUCommonMath::Max(measurement.radius, mMaxR[layer]);
        mMinZ[layer] = o2::gpu::GPUCommonMath::Min(measurement.z, mMinZ[layer]);
        mMaxZ[layer] = o2::gpu::GPUCommonMath::Max(measurement.z, mMaxZ[layer]);
      }
      std::exclusive_scan(counts.begin(), counts.end(), offsets.begin(), 0);

      for (const auto& helper : helpers) {
        sortedMeasurements[offsets[helper.bin] + helper.indexWithinBin] = mLayerGlobalMeasurements[layer][helper.measurementIndex];
      }
      std::copy(sortedMeasurements.begin(), sortedMeasurements.end(), mLayerGlobalMeasurements[layer].begin() + first);
      std::copy_n(offsets.data(), counts.size(), tableBase);
      std::fill_n(tableBase + counts.size(), stride - counts.size(), count);
      std::fill(counts.begin(), counts.end(), 0);
      helpers.clear();
      sortedMeasurements.clear();
    }
  }
}

} // namespace o2::itsmft::tracking
