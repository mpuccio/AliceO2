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
#include "ITSMFTTrackingGPU/TrackerTraitsGPU.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <oneapi/tbb.h>
#include "Framework/Logger.h"
#include "ITSMFTTracking/CapacityEstimator.h"

namespace o2::itsmft::tracking
{
namespace gpu
{
inline namespace ITSMFT_TRACKING_GPU_BACKEND
{
namespace
{
bool artefactLabels(const IterationContext& context)
{
  return context.frame.hasMCinformation() && context.configuration.parameters.CreateArtefactLabels;
}
template <typename T>
std::vector<T> copyOf(const bounded_vector<T>& values)
{
  return {values.begin(), values.end()};
}
} // namespace

const char* TrackerTraitsGPU::getName() const noexcept
{
#ifdef ITSMFT_TRACKING_HIP
  return "HIP tracklets + cells + neighbours + roads + refit";
#else
  return "CUDA tracklets + cells + neighbours + roads + refit";
#endif
}

void TrackerTraitsGPU::computeLayerTracklets(IterationContext& context, int iteration, int iVertex)
{
  // As the legacy tracker: the timeframe is uploaded before its first
  // tracklet stage (once the clusters are sorted), each iteration's ROF masks
  // before its first vertex.
  if (!mTimeFrameGPU.isLoaded()) {
    mTimeFrameGPU.loadTimeFrame(context.frame, context.configuration.topology.nLayers);
  }
  if (iVertex <= 0) {
    mTimeFrameGPU.loadROFEnabled(context.frame);
  }
  auto& scratch = context.scratch;
  const auto nEdges = scratch.getTracklets().size();
  for (size_t edge = 0; edge < nEdges; ++edge) {
    scratch.getTracklets()[edge].clear();
    scratch.getTrackletsLabel(edge).clear();
    std::fill(scratch.getTrackletsLookupTable()[edge].begin(), scratch.getTrackletsLookupTable()[edge].end(), 0);
  }
  mTimeFrameGPU.initialiseTracklets(nEdges);
  const bool download = validateTracklets() || validateCells() || artefactLabels(context);
  const auto edgeIds = context.configuration.edgeIds();
  auto& estimator = context.frame.getCapacityEstimator();
  executeInArena([&] {
    tbb::parallel_for(0, static_cast<int>(edgeIds.size()), [&](const int index) {
      const auto edgeId = edgeIds[index];
      const int edge = edgeId.value();
      const auto search = prepareTrackletSearch(context, iVertex, edgeId);
      auto& results = mTimeFrameGPU.tracklets(edge);
      const auto key = CapacityEstimator::makeKey(SlabSite::Tracklets, iteration, iVertex + 1, edgeId);
      runOnSlab(estimator, key, mTimeFrameGPU.view().layers[search.edgeCache.fromLayer].nClusters, [&](const int capacity) {
        results.reserve(capacity);
        return computeTracklets(mTimeFrameGPU, search, edge);
      });
      if (download) {
        mTimeFrameGPU.download(results.items, scratch.getTracklets()[edge]);
        mTimeFrameGPU.download(results.lookup, scratch.getTrackletsLookupTable()[edge]);
      }
    });
  });
  createTrackletLabels(context);
  if (!validateTracklets()) {
    return;
  }
  std::vector<std::vector<Tracklet>> tracklets(nEdges);
  std::vector<std::vector<int>> lookups(nEdges);
  for (const auto edgeId : edgeIds) {
    tracklets[edgeId.value()] = copyOf(scratch.getTracklets()[edgeId.value()]);
    lookups[edgeId.value()] = copyOf(scratch.getTrackletsLookupTable()[edgeId.value()]);
  }
  TrackerTraits::computeLayerTracklets(context, iteration, iVertex);
  const auto close = [](float a, float b) {
    return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= 2.e-6f * std::max({1.f, std::abs(a), std::abs(b)});
  };
  const auto same = [&](const Tracklet& a, const Tracklet& b) {
    return a == b && a.getTimeStamp().getTimeStamp() == b.getTimeStamp().getTimeStamp() &&
           a.getTimeStamp().getTimeStampError() == b.getTimeStamp().getTimeStampError() && close(a.tanLambda, b.tanLambda) && close(a.phi, b.phi);
  };
  for (const auto edgeId : edgeIds) {
    const auto edge = edgeId.value();
    auto& cpu = scratch.getTracklets()[edge];
    const auto& gpu = tracklets[edge];
    const auto& cpuLookup = scratch.getTrackletsLookupTable()[edge];
    if (!std::equal(cpu.begin(), cpu.end(), gpu.begin(), gpu.end(), same) ||
        !std::equal(lookups[edge].begin(), lookups[edge].end(), cpuLookup.begin(), cpuLookup.end())) {
      throw std::runtime_error{"CPU/GPU tracklet mismatch on edge " + std::to_string(edge)};
    }
    LOGP(info, "CPU/GPU tracklet validation: iteration={} vertex={} edge={} candidates={} matched", iteration, iVertex, edge, cpu.size());
    std::copy(gpu.begin(), gpu.end(), cpu.begin()); // later device stages read the device's
  }
}

void TrackerTraitsGPU::adoptHostTracklets(IterationContext& context)
{
  if (!mTimeFrameGPU.isLoaded()) {
    mTimeFrameGPU.loadTimeFrame(context.frame, context.configuration.topology.nLayers);
  }
  auto& scratch = context.scratch;
  mTimeFrameGPU.initialiseTracklets(scratch.getTracklets().size());
  for (size_t edge = 0; edge < scratch.getTracklets().size(); ++edge) {
    auto& results = mTimeFrameGPU.tracklets(edge);
    const auto& tracklets = scratch.getTracklets()[edge];
    const auto& lookup = scratch.getTrackletsLookupTable()[edge];
    mTimeFrameGPU.upload(results.items, std::span<const Tracklet>{tracklets.data(), tracklets.size()});
    mTimeFrameGPU.upload(results.lookup, std::span<const int>{lookup.data(), lookup.size()});
  }
}

void TrackerTraitsGPU::computeLayerCells(IterationContext& context, int iteration)
{
  auto& scratch = context.scratch;
  const auto nPaths = scratch.getCells().size();
  for (size_t path = 0; path < nPaths; ++path) {
    deepVectorClear(scratch.getCells()[path]);
    deepVectorClear(scratch.getCellsLookupTable()[path]);
    if (artefactLabels(context)) {
      deepVectorClear(scratch.getCellsLabel(path));
    }
  }
  mTimeFrameGPU.initialiseCells(nPaths);
  const bool download = validateCells() || validateNeighbours() || validateRoads() || artefactLabels(context);
  const auto cellIds = context.configuration.cellIds();
  auto& estimator = context.frame.getCapacityEstimator();
  executeInArena([&] {
    tbb::parallel_for(0, static_cast<int>(cellIds.size()), [&](const int index) {
      const auto cellId = cellIds[index];
      const int path = cellId.value();
      const int firstEdge = context.topology.getPath(cellId).first.value();
      const int secondEdge = context.topology.getPath(cellId).second.value();
      const auto nTracklets = mTimeFrameGPU.tracklets(firstEdge).items.size;
      if (!nTracklets || !mTimeFrameGPU.tracklets(secondEdge).items.size) {
        return;
      }
      const auto search = prepareCellSearch(context, cellId);
      auto& results = mTimeFrameGPU.cells(path);
      const auto key = CapacityEstimator::makeKey(SlabSite::Cells, iteration, 0, cellId);
      runOnSlab(estimator, key, static_cast<double>(nTracklets), [&](const int capacity) {
        results.reserve(capacity);
        return computeCells(mTimeFrameGPU, search, firstEdge, secondEdge, path);
      });
      if (download) {
        mTimeFrameGPU.download(results.items, scratch.getCells()[path]);
        mTimeFrameGPU.download(results.lookup, scratch.getCellsLookupTable()[path]);
      }
    });
  });
  if (!validateCells()) {
    finishCells(context);
    return;
  }
  std::vector<std::vector<Triplet>> cells(nPaths);
  std::vector<std::vector<int>> lookups(nPaths);
  for (const auto cellId : cellIds) {
    cells[cellId.value()] = copyOf(scratch.getCells()[cellId.value()]);
    lookups[cellId.value()] = copyOf(scratch.getCellsLookupTable()[cellId.value()]);
  }
  TrackerTraits::computeLayerCells(context, iteration); // also finishes the stage
  const auto same = [](const Triplet& a, const Triplet& b) {
    return a.getFirstTrackletIndex() == b.getFirstTrackletIndex() && a.getSecondTrackletIndex() == b.getSecondTrackletIndex() &&
           a.getClusters() == b.getClusters() && a.getHitLayerMask().value() == b.getHitLayerMask().value() && a.getLevel() == b.getLevel() &&
           a.getTimeStamp().getTimeStamp() == b.getTimeStamp().getTimeStamp() &&
           a.getTimeStamp().getTimeStampError() == b.getTimeStamp().getTimeStampError() && cellFactorsEquivalent(b.tripletFactor(), a.tripletFactor());
  };
  for (const auto cellId : cellIds) {
    const auto path = cellId.value();
    auto& cpu = scratch.getCells()[path];
    const auto& gpu = cells[path];
    const auto& cpuLookup = scratch.getCellsLookupTable()[path];
    const bool noLookups = cpuLookup.empty() && std::all_of(lookups[path].begin(), lookups[path].end(), [](int v) { return v == 0; });
    if (!std::equal(cpu.begin(), cpu.end(), gpu.begin(), gpu.end(), same) ||
        !(noLookups || std::equal(lookups[path].begin(), lookups[path].end(), cpuLookup.begin(), cpuLookup.end()))) {
      throw std::runtime_error{"CPU/GPU cell mismatch on path " + std::to_string(path)};
    }
    LOGP(info, "CPU/GPU cell validation: iteration={} path={} cells={} matched", iteration, path, cpu.size());
    std::copy(gpu.begin(), gpu.end(), cpu.begin()); // later device stages read the device's
  }
}

void TrackerTraitsGPU::findCellsNeighbours(IterationContext& context, int iteration)
{
  auto& scratch = context.scratch;
  const auto nPaths = scratch.getCells().size();
  for (size_t path = 0; path < nPaths; ++path) {
    deepVectorClear(scratch.getCellsNeighbours()[path]);
    deepVectorClear(scratch.getCellsNeighboursLUT()[path]);
  }
  mTimeFrameGPU.initialiseNeighbours(nPaths);
  auto& estimator = context.frame.getCapacityEstimator();
  const auto nCells = [&](CellPathId path) { return mTimeFrameGPU.cells(path.value()).items.size; };
  forEachNeighbourTarget(context, nCells, [&](CellPathId target, std::span<NeighbourSearch> sources, size_t nSourceCells) {
    auto& results = mTimeFrameGPU.neighbours(target.value());
    const auto key = CapacityEstimator::makeKey(SlabSite::Neighbours, iteration, 0, target);
    runOnSlab(estimator, key, static_cast<double>(nSourceCells), [&](const int capacity) {
      results.reserve(capacity);
      return computeNeighbours(mTimeFrameGPU, target.value(), sources);
    });
  });
  if (!validateNeighbours()) {
    if (validateRoads()) {
      downloadRoadGraph(context); // the host builds the road jobs to validate them
    }
    return;
  }
  // The CPU stage starts from the cells as the cell stage left them.
  TrackerTraits::findCellsNeighbours(context, iteration);
  std::vector<std::vector<int>> levels(nPaths), lookups(nPaths);
  std::vector<std::vector<CellNeighbour>> neighbours(nPaths);
  for (const auto cellId : context.configuration.cellIds()) {
    const int path = cellId.value();
    for (const auto& cell : scratch.getCells()[path]) {
      levels[path].push_back(cell.getLevel());
    }
    lookups[path] = copyOf(scratch.getCellsNeighboursLUT()[path]);
    neighbours[path] = copyOf(scratch.getCellsNeighbours()[path]);
    deepVectorClear(scratch.getCellsNeighbours()[path]);
    deepVectorClear(scratch.getCellsNeighboursLUT()[path]);
  }
  downloadRoadGraph(context);
  const auto same = [](const CellNeighbour& a, const CellNeighbour& b) { return a.cell == b.cell && a.cellPath == b.cellPath && a.nextCell == b.nextCell; };
  const auto sameLevel = [](const Triplet& cell, int level) { return cell.getLevel() == level; };
  size_t nNeighbours = 0;
  for (const auto cellId : context.configuration.cellIds()) {
    const int path = cellId.value();
    const auto& cells = scratch.getCells()[path];
    const auto& gpuLookup = scratch.getCellsNeighboursLUT()[path];
    const auto& gpuNeighbours = scratch.getCellsNeighbours()[path];
    if (!std::equal(cells.begin(), cells.end(), levels[path].begin(), levels[path].end(), sameLevel) ||
        !std::equal(lookups[path].begin(), lookups[path].end(), gpuLookup.begin(), gpuLookup.end()) ||
        !std::equal(neighbours[path].begin(), neighbours[path].end(), gpuNeighbours.begin(), gpuNeighbours.end(), same)) {
      throw std::runtime_error{"CPU/GPU neighbour mismatch on path " + std::to_string(path)};
    }
    nNeighbours += gpuNeighbours.size();
  }
  LOGP(info, "GPU neighbour validation: {} neighbours", nNeighbours);
}

void TrackerTraitsGPU::downloadRoadGraph(IterationContext& context)
{
  auto& scratch = context.scratch;
  for (const auto cellId : context.configuration.cellIds()) {
    const int path = cellId.value();
    mTimeFrameGPU.download(mTimeFrameGPU.cells(path).items, scratch.getCells()[path]);
    const auto& neighbours = mTimeFrameGPU.neighbours(path);
    if (neighbours.items.size) {
      mTimeFrameGPU.download(neighbours.items, scratch.getCellsNeighbours()[path]);
      mTimeFrameGPU.download(neighbours.lookup, scratch.getCellsNeighboursLUT()[path]);
    }
  }
}

void TrackerTraitsGPU::findRoads(IterationContext& context, int iteration)
{
  // As the legacy tracker: every start path's roads are extended level by
  // level on the device, and each start level's selected seeds are refitted
  // and sorted there; only the accepted tracks come back.
  bindRoadGraph(mTimeFrameGPU);
  auto& scratch = context.scratch;
  const auto& pool = scratch.getMemoryPool();
  const auto catalog = context.topology.getSurfaceCatalogView();
  const std::span<const SurfaceDescriptor> surfaces{catalog.surfaces, catalog.nSurfaces};
  const float maxChi2 = context.configuration.kernelParameters.maxChi2ClusterAttachment;
  auto firstClusters = makeFirstClusters(context);
  forEachRoadStartLevel(context, iteration, [&](std::span<const CellPathId> starts, int startLevel) {
    const auto selector = makeRoadSeedSelector(context, startLevel);
    mTimeFrameGPU.trackSeeds.size = 0;
    for (const auto start : starts) {
      if (!mTimeFrameGPU.cells(start.value()).items.size) {
        continue;
      }
      RoadExtension extension{start.value(), false, startLevel, context.configuration.topology.nLayers, context.bz, maxChi2, surfaces};
      // With validation, every extension is checked against the host's.
      bounded_vector<RoadSeedEmission> seeds(pool.get());
      const auto extend = [&](const auto& hostSeeds) {
        const size_t n = extendRoads(mTimeFrameGPU, extension);
        if (!validateRoads()) {
          return n;
        }
        bounded_vector<RoadJob> jobs(pool.get());
        bounded_vector<RoadTarget> targets(pool.get());
        buildRoadJobs(context, start.value(), extension.currentLevel, hostSeeds, jobs, targets);
        bounded_vector<RoadResult> results(pool.get());
        mTimeFrameGPU.download(mTimeFrameGPU.roads[mTimeFrameGPU.lastRoads], results);
        size_t cursor = 0;
        for (const auto& job : jobs) {
          forEachRoad(job.seed, job.initialize ? &job.initial : nullptr, targets.data(), job.firstTarget, job.endTarget, context.bz, maxChi2,
                      [&](const RoadSeedEmission& expected) {
                        if (cursor >= results.size() || results[cursor].job != job.source || !roadSeedsEquivalent(expected, results[cursor].emission)) {
                          throw std::runtime_error{"GPU road validation: mismatch at extension " + std::to_string(cursor)};
                        }
                        ++cursor;
                      });
        }
        if (cursor != results.size()) {
          throw std::runtime_error{"GPU road validation: extension count mismatch"};
        }
        LOGP(info, "GPU road validation: {} extensions", cursor);
        seeds.clear();
        for (const auto& result : results) {
          seeds.push_back(result.emission);
        }
        return n;
      };
      size_t n = extend(scratch.getCells()[start.value()]);
      for (int level = startLevel - 1; level >= 2 && n; --level) {
        extension.reusePreviousSeeds = true;
        extension.currentLevel = level;
        const auto previous = seeds;
        n = extend(previous);
      }
      if (n) {
        selectRoadSeeds(mTimeFrameGPU, selector);
      }
    }
    if (!mTimeFrameGPU.trackSeeds.size) {
      return;
    }
    const auto& minPt = context.configuration.parameters.MinPt;
    refitTrackSeeds(mTimeFrameGPU, surfaces, makeRefitParameters(context), {minPt.data(), std::min<size_t>(minPt.size(), MaxLayoutSurfaces)});
    if (validateRefit()) {
      bounded_vector<TrackSeed> seeds(pool.get());
      bounded_vector<RefitResult> results(pool.get());
      mTimeFrameGPU.download(mTimeFrameGPU.trackSeeds, seeds);
      mTimeFrameGPU.download(mTimeFrameGPU.refitResults, results);
      const auto& jobs = prepareRefitJobs(context, seeds);
      const auto parameters = makeRefitParameters(context);
      for (size_t i = 0; i < jobs.size(); ++i) {
        if (!refitsEquivalent(evaluateRefit(jobs[i], catalog, parameters), results[i])) {
          throw std::runtime_error{"GPU refit mismatch at seed " + std::to_string(i)};
        }
      }
      LOGP(info, "Validated {} GPU refits against CPU", jobs.size());
    }
    mTimeFrameGPU.download(mTimeFrameGPU.tracks, mTracks);
    acceptSortedTracks(context, iteration, mTracks, firstClusters);
  });
}
} // namespace ITSMFT_TRACKING_GPU_BACKEND
} // namespace gpu

#ifdef ITSMFT_TRACKING_HIP
std::unique_ptr<TrackerTraits> createTrackerTraitsHIP()
#else
std::unique_ptr<TrackerTraits> createTrackerTraitsCUDA()
#endif
{
  return std::make_unique<gpu::TrackerTraitsGPU>();
}
} // namespace o2::itsmft::tracking
