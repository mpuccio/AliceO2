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
///
/// \file TrackerTraits.cxx
/// \brief
///

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_sort.h>

#include "CommonConstants/MathConstants.h"
#include "Framework/Logger.h"
#include "GPUCommonMath.h"
#include "ITSMFTTracking/TrackSeed.h"
#include "ITSMFTTracking/BoundedAllocator.h"
#include "ITSMFTTracking/Triplet.h"
#include "ITSMFTTracking/CapacityEstimator.h"
#include "ITSMFTTracking/SlabBumpAllocator.h"
#include "ITSMFTTracking/Constants.h"
#include "ITSMFTTracking/MathUtils.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/RefitDriver.h"
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/LayerMask.h"
#include "ITSMFTTracking/TripletFitting.h"
#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/CandidateFinding.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::itsmft::tracking
{

namespace math_utils = o2::its::math_utils;
using o2::its::TimeEstBC;

namespace
{

void reserveGenericTrackPublication(TimeFrame& frame, std::size_t candidateCount, std::size_t maxReferencesPerTrack)
{
  auto& tracks = frame.getGenericTracks();
  auto& references = frame.getTrackClusterIndices();
  if (candidateCount > tracks.max_size() - tracks.size() ||
      (maxReferencesPerTrack != 0 && candidateCount > (references.max_size() - references.size()) / maxReferencesPerTrack)) {
    throw std::length_error{"GenericTrack publication exceeds the output container capacity"};
  }
  tracks.reserve(tracks.size() + candidateCount);
  references.reserve(references.size() + candidateCount * maxReferencesPerTrack);
}

bool appendGenericTrack(TimeFrame& frame,
                        const TrackingCandidate& candidate,
                        gsl::span<const gsl::span<const GlobalMeasurement>> layerMeasurements)
{
  GenericTrack track = candidate.track;
  track.hitLayers = {};
  std::vector<TrackClusterReference> resolvedReferences;
  resolvedReferences.reserve(layerMeasurements.size());
  for (std::size_t position = 0; position < layerMeasurements.size(); ++position) {
    const int localIndex = candidate.getClusterIndex(static_cast<int>(position));
    if (localIndex == o2::its::constants::UnusedIndex) {
      continue;
    }
    if (localIndex < 0 || static_cast<std::size_t>(localIndex) >= layerMeasurements[position].size()) {
      return false;
    }
    const auto& measurement = layerMeasurements[position][localIndex];
    const TrackClusterReference reference{LayerId{static_cast<uint16_t>(position)}, 0, measurement.clusterId};
    if (!reference.isValid()) {
      return false;
    }
    resolvedReferences.push_back(reference);
    track.hitLayers.set(static_cast<int>(position));
  }
  if (!track.innerState.hasRecognizedKind() || !track.outerState.hasRecognizedKind() ||
      !std::isfinite(track.timestamp.getTimeStamp()) ||
      !std::isfinite(track.timestamp.getTimeStampError()) || track.timestamp.getTimeStampError() <= 0.f || resolvedReferences.empty()) {
    return false;
  }

  auto& tracks = frame.getGenericTracks();
  auto& references = frame.getTrackClusterIndices();
  const auto oldTrackSize = tracks.size();
  const auto oldReferenceSize = references.size();
  if (oldTrackSize > std::numeric_limits<uint32_t>::max() || oldReferenceSize > std::numeric_limits<uint32_t>::max() ||
      resolvedReferences.size() > std::numeric_limits<uint32_t>::max() - oldReferenceSize) {
    return false;
  }

  try {
    for (const auto& reference : resolvedReferences) {
      references.push_back(reference);
    }
    track.firstClusterRef = static_cast<uint32_t>(oldReferenceSize);
    track.clusterRefEnd = static_cast<uint32_t>(references.size());
    tracks.push_back(track);
  } catch (...) {
    references.resize(oldReferenceSize);
    tracks.resize(oldTrackSize);
    throw;
  }
  return true;
}

} // namespace

void TrackerTraits::runTraversal(IterationContext& view)
{
  if (view.iteration < 0) {
    throw std::invalid_argument{"CA traversal: iteration out of range (iteration " + std::to_string(view.iteration) + ")"};
  }
  int maxNvertices{-1};
  if (view.configuration.parameters.PerPrimaryVertexProcessing) {
    maxNvertices = view.frame.getMaxVerticesPerROF();
  }
  const auto timed = [&](Stage stage, auto&& run) {
    const auto start = std::chrono::steady_clock::now();
    run();
    mStageMs[stage] += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
  };
  int iVertex = std::min(maxNvertices, 0);
  do {
    beginPass(view);
    try {
      timed(Tracklets, [&] { computeLayerTracklets(view, view.iteration, iVertex); });
      timed(Cells, [&] { computeLayerCells(view, view.iteration); });
      timed(Neighbours, [&] { findCellsNeighbours(view, view.iteration); });
      timed(Roads, [&] { findRoads(view, view.iteration); });
    } catch (...) {
      endPass(view);
      throw;
    }
    endPass(view);
  } while (++iVertex < maxNvertices);
}

void TrackerTraits::computeLayerTracklets(IterationContext& context, const int iteration, int iVertex)
{
  auto& scratch = context.scratch;
  for (size_t edge = 0; edge < scratch.getTracklets().size(); ++edge) {
    scratch.getTracklets()[edge].clear();
    scratch.getTrackletsLabel(edge).clear();
    std::fill(scratch.getTrackletsLookupTable()[edge].begin(), scratch.getTrackletsLookupTable()[edge].end(), 0);
  }
  const auto frame = context.frame.makeView(context.layerGlobalMeasurements);
  const auto edgeIds = context.configuration.edgeIds();
  const int maxConcurrency = std::max(1, mTaskArena->max_concurrency());
  const int nConcurrentSinks = std::min(static_cast<int>(edgeIds.size()), maxConcurrency);
  auto& estimator = context.frame.getCapacityEstimator();
  mTaskArena->execute([&] {
    tbb::parallel_for(0, static_cast<int>(edgeIds.size()), [&](const int edgeIndex) {
      const auto edgeId = edgeIds[edgeIndex];
      const auto edge = prepareTrackletSearch(context, iVertex, edgeId);
      const auto& from = frame.layers[edge.edgeCache.fromLayer];
      const int nSources = static_cast<int>(from.nClusters);
      const auto key = CapacityEstimator::makeKey(SlabSite::Tracklets, iteration, iVertex + 1, edgeId);
      const auto capacity = estimator.capacity(key, nSources);
      UnorderedSlabSink<Tracklet> sink{{.capacity = capacity, .nThreads = maxConcurrency, .nConcurrentSinks = nConcurrentSinks}, scratch.getMemoryPool().get()};
      tbb::parallel_for(0, nSources, [&](const int source) {
        auto& handle = sink.local();
        forEachTracklet(frame, edge, source, [&](const Tracklet& tracklet) { handle.emplace(tracklet); });
      });
      const auto stats = sink.stats();
      auto& tracklets = scratch.getTracklets()[edgeId.value()];
      sink.finalizeUnordered(tracklets);
      estimator.update(key, nSources, stats.requested, stats.capacity, stats.emitted, stats.spilled, stats.overflowed, stats.memoryLimited);
      // The same pair can be found from several vertices, always with the
      // same content, so equal tracklets may be sorted in any order.
      tbb::parallel_sort(tracklets.begin(), tracklets.end());
      tracklets.erase(std::unique(tracklets.begin(), tracklets.end()), tracklets.end());
      auto& lookup = scratch.getTrackletsLookupTable()[edgeId.value()];
      for (const auto& tracklet : tracklets) {
        ++lookup[tracklet.firstClusterIndex + 1];
      }
      std::inclusive_scan(lookup.begin(), lookup.end(), lookup.begin());
    });
  });
  createTrackletLabels(context);
}

void TrackerTraits::createTrackletLabels(IterationContext& context)
{
  auto& scratch = context.scratch;
  auto* mFrame = &context.frame;
  if (!mFrame->hasMCinformation() || !context.configuration.parameters.CreateArtefactLabels) {
    return;
  }
  const auto edgeIds = context.configuration.edgeIds();
  const auto& topology = context.topology;
  executeInArena([&] {
    tbb::parallel_for(0, static_cast<int>(edgeIds.size()), [&](const int edgeIndex) {
      const auto edgeId = edgeIds[edgeIndex];
      const auto& edge = topology.getEdge(edgeId);
      const int fromLayer = edge.from.value();
      const int toLayer = edge.to.value();
      for (auto& trk : scratch.getTracklets()[edgeId.value()]) {
        MCCompLabel label;
        const auto currentId = mFrame->getClusters()[fromLayer][trk.firstClusterIndex].clusterId;
        const auto nextId = mFrame->getClusters()[toLayer][trk.secondClusterIndex].clusterId;
        for (const auto& lab1 : mFrame->getLabels(LayerId{static_cast<uint16_t>(fromLayer)}, currentId)) {
          for (const auto& lab2 : mFrame->getLabels(LayerId{static_cast<uint16_t>(toLayer)}, nextId)) {
            if (lab1 == lab2 && lab1.isValid()) {
              label = lab1;
              break;
            }
          }
          if (label.isValid()) {
            break;
          }
        }
        scratch.getTrackletsLabel(edgeId.value()).emplace_back(label);
      }
    });
  });
}

TrackletSearch TrackerTraits::prepareTrackletSearch(IterationContext& context, int iVertex, EdgeId edgeId) const
{
  const auto& frame = context.frame;
  const auto& parameters = context.configuration.parameters;
  const auto& edge = context.topology.getEdge(edgeId);
  const int fromLayer = edge.from.value();
  const int toLayer = edge.to.value();
  const int fromLocal = frame.getROFLocalLayer(fromLayer);
  const int toLocal = frame.getROFLocalLayer(toLayer);
  const auto& overlapView = frame.getROFViews(fromLayer).overlap;
  const auto& overlap = overlapView.mIndices[(fromLocal * overlapView.mLayerCount) + toLocal];
  TrackletSearch search;
  search.overlapRanges = overlapView.mFlatTable + overlap.getFirstEntry();
  search.nOverlapRanges = overlap.getEntries();
  if (!parameters.UseDiamond) {
    const auto& vertexView = frame.getROFViews(fromLayer).vertexLookup;
    const auto& vertices = vertexView.mIndices[fromLocal];
    search.vertexRanges = vertexView.mFlatTable + vertices.getFirstEntry();
    search.nVertexRanges = vertices.getEntries();
  }
  search.fromLayerTiming = overlapView.getLayer(fromLocal);
  search.toLayerTiming = overlapView.getLayer(toLocal);
  search.edgeCache = {fromLayer, toLayer, context.detectorConfiguration.getRepresentativeRadius(edge.from),
                      context.detectorConfiguration.getRepresentativeRadius(edge.to), frame.getMinR(toLayer), frame.getMaxR(toLayer),
                      frame.getMinZ(toLayer), frame.getMaxZ(toLayer), context.detectorConfiguration.positionResolutions[fromLayer],
                      context.scratch.getEdgeMSAngle(edgeId.value()), context.scratch.getEdgePhiCut(edgeId.value())};
  search.indexTableUtils = frame.getIndexTableUtils(toLayer);
  search.nSigmaCut = context.configuration.kernelParameters.nSigmaCut;
  search.beamPositionVariance = frame.getBeamPositionVariance();
  search.useDiamond = parameters.UseDiamond;
  search.diamondBase = Vertex(parameters.Diamond, parameters.DiamondCov, 1, 1.f);
  search.selectUPC = parameters.PassFlags[IterationStep::SelectUPCVertices];
  search.iVertex = iVertex;
  search.kind = context.topology.getSurface(edge.from).kind;
  return search;
}

void TrackerTraits::computeLayerCells(IterationContext& context, const int iteration)
{
  auto& scratch = context.scratch;
  for (size_t path = 0; path < scratch.getCells().size(); ++path) {
    deepVectorClear(scratch.getCells()[path]);
    deepVectorClear(scratch.getCellsLookupTable()[path]);
    if (context.frame.hasMCinformation() && context.configuration.parameters.CreateArtefactLabels) {
      deepVectorClear(scratch.getCellsLabel(path));
    }
  }
  const auto frame = context.frame.makeView(context.layerGlobalMeasurements);
  const int maxConcurrency = std::max(1, mTaskArena->max_concurrency());
  auto& estimator = context.frame.getCapacityEstimator();
  const auto cellIds = context.configuration.cellIds();
  const int nConcurrentSinks = std::min(static_cast<int>(cellIds.size()), maxConcurrency);
  // Paths write disjoint outputs, so they run concurrently, as the legacy
  // cell topologies do.
  mTaskArena->execute([&] {
    tbb::parallel_for(size_t{0}, static_cast<size_t>(cellIds.size()), [&](size_t index) {
      const auto cellId = cellIds[index];
      const auto& path = context.topology.getPath(cellId);
      const auto& first = scratch.getTracklets()[path.first.value()];
      const auto& second = scratch.getTracklets()[path.second.value()];
      if (first.empty() || second.empty()) {
        return;
      }
      auto search = prepareCellSearch(context, cellId);
      search.first = first.data();
      search.second = second.data();
      search.secondLookup = scratch.getTrackletsLookupTable()[path.second.value()].data();
      const int nFirst = static_cast<int>(first.size());
      const auto key = CapacityEstimator::makeKey(SlabSite::Cells, iteration, 0, cellId);
      GroupedSlabSink<Triplet> sink{{.capacity = estimator.capacity(key, nFirst), .nThreads = maxConcurrency, .nConcurrentSinks = nConcurrentSinks}, scratch.getMemoryPool().get()};
      tbb::parallel_for(0, nFirst, [&](const int tracklet) {
        auto& handle = sink.local();
        handle.beginProducer(tracklet);
        forEachCell(frame, search, tracklet, [&](const Triplet& cell) { handle.emplace(cell); });
      });
      const auto stats = sink.stats();
      sink.finalizeGrouped(first.size(), scratch.getCellsLookupTable()[cellId.value()], scratch.getCells()[cellId.value()]);
      estimator.update(key, nFirst, stats.requested, stats.capacity, stats.emitted, stats.spilled, stats.overflowed, stats.memoryLimited);
    });
  });
  finishCells(context);
}

CellSearch TrackerTraits::prepareCellSearch(IterationContext& context, CellPathId cellId) const
{
  const auto& topology = context.topology;
  const auto& path = topology.getPath(cellId);
  const auto& firstEdge = topology.getEdge(path.first);
  const auto& kernel = context.configuration.kernelParameters;
  CellSearch search;
  search.layers = {firstEdge.from.value(), firstEdge.to.value(), topology.getEdge(path.second).to.value()};
  search.parameters = {kernel.nSigmaCut * context.scratch.getEdgeMSAngle(path.second.value()),
                       std::abs(o2::constants::math::B2C * context.bz) / kernel.trackletMinPt,
                       topology.getSurface(firstEdge.to).kind == SurfaceKind::Disk};
  return search;
}

void TrackerTraits::finishCells(IterationContext& context)
{
  auto& scratch = context.scratch;
  if (context.frame.hasMCinformation() && context.configuration.parameters.CreateArtefactLabels) {
    for (const auto cellId : context.configuration.cellIds()) {
      const auto& path = context.topology.getPath(cellId);
      auto& labels = scratch.getCellsLabel(cellId.value());
      const auto& cells = scratch.getCells()[cellId.value()];
      labels.reserve(cells.size());
      for (const auto& cell : cells) {
        MCCompLabel currentLab{scratch.getTrackletsLabel(path.first.value())[cell.getFirstTrackletIndex()]};
        MCCompLabel nextLab{scratch.getTrackletsLabel(path.second.value())[cell.getSecondTrackletIndex()]};
        labels.emplace_back(currentLab == nextLab ? currentLab : MCCompLabel());
      }
    }
  }
  const auto scratchEdgeCount = scratch.getTracklets().size();
  for (size_t edgeId = 0; edgeId < scratchEdgeCount; ++edgeId) {
    deepVectorClear(scratch.getTracklets()[edgeId]);
    deepVectorClear(scratch.getTrackletsLabel(edgeId));
  }
}

void TrackerTraits::findCellsNeighbours(IterationContext& context, const int iteration)
{
  auto& scratch = context.scratch;
  for (size_t path = 0; path < scratch.getCellsNeighbours().size(); ++path) {
    deepVectorClear(scratch.getCellsNeighbours()[path]);
    deepVectorClear(scratch.getCellsNeighboursLUT()[path]);
  }
  const auto frame = context.frame.makeView(context.layerGlobalMeasurements);
  const int maxConcurrency = std::max(1, mTaskArena->max_concurrency());
  const int nConcurrentSinks = std::min(static_cast<int>(context.configuration.topology.scheduledPaths.size()), maxConcurrency);
  auto& estimator = context.frame.getCapacityEstimator();
  const auto nCells = [&](CellPathId path) { return scratch.getCells()[path.value()].size(); };
  forEachNeighbourTarget(context, nCells, [&](CellPathId target, std::span<NeighbourSearch> sources, size_t nSourceCells) {
    auto& cells = scratch.getCells()[target.value()];
    const auto& lookup = scratch.getCellsLookupTable()[target.value()];
    const auto key = CapacityEstimator::makeKey(SlabSite::Neighbours, iteration, 0, target);
    UnorderedSlabSink<CellNeighbour> sink{{.capacity = estimator.capacity(key, nSourceCells), .nThreads = maxConcurrency, .nConcurrentSinks = nConcurrentSinks},
                                          scratch.getMemoryPool().get()};
    for (auto& search : sources) {
      const auto& sourceCells = scratch.getCells()[search.sourcePath];
      search.sources = sourceCells.data();
      search.targets = cells.data();
      search.targetLookup = lookup.data();
      search.lookupSize = lookup.size();
      tbb::parallel_for(0, static_cast<int>(sourceCells.size()), [&](int source) {
        forEachNeighbour(frame, search, source, [&](const CellNeighbour& neighbour) { sink.local().emplace(neighbour); });
      });
    }
    const auto stats = sink.stats();
    auto& neighbours = scratch.getCellsNeighbours()[target.value()];
    sink.finalizeUnordered(neighbours);
    estimator.update(key, nSourceCells, stats.requested, stats.capacity, stats.emitted, stats.spilled, stats.overflowed, stats.memoryLimited);
    if (neighbours.empty()) {
      return;
    }
    // Every (target cell, source path, source cell) key is unique.
    tbb::parallel_sort(neighbours.begin(), neighbours.end(), [](const auto& a, const auto& b) {
      return std::tie(a.nextCell, a.cellPath, a.cell) < std::tie(b.nextCell, b.cellPath, b.cell);
    });
    auto& neighbourLookup = scratch.getCellsNeighboursLUT()[target.value()];
    neighbourLookup.assign(cells.size() + 1, 0);
    for (const auto& neighbour : neighbours) {
      ++neighbourLookup[neighbour.nextCell + 1];
      auto& cell = cells[neighbour.nextCell];
      cell.setLevel(std::max(cell.getLevel(), scratch.getCells()[neighbour.cellPath][neighbour.cell].getLevel() + 1));
    }
    std::inclusive_scan(neighbourLookup.begin(), neighbourLookup.end(), neighbourLookup.begin());
  });
  for (auto& lookup : scratch.getCellsLookupTable()) {
    deepVectorClear(lookup);
  }
}

bool TrackerTraits::buildTrackSeed(IterationContext& context, int, const Triplet& cell, TrackSeed& output) const
{
  SeedInput input;
  const auto frame = context.frame.makeView(context.layerGlobalMeasurements);
  throwRoadError(makeSeedInput(frame, context.topology.getSurfaceCatalogView(), resolveCell(frame, cell), &input));
  return initializeTrackSeed(input, context.bz, context.configuration.kernelParameters.maxChi2ClusterAttachment, output);
}

template <typename InputSeed>
void TrackerTraits::buildRoadJobs(IterationContext& context, int startPath, int currentLevel, const bounded_vector<InputSeed>& seeds,
                                  bounded_vector<RoadJob>& jobs, bounded_vector<RoadTarget>& targets) const
{
  constexpr bool startCells = std::is_same_v<InputSeed, Triplet>;
  auto& scratch = context.scratch;
  const auto& pool = scratch.getMemoryPool();
  const auto frame = context.frame.makeView(context.layerGlobalMeasurements);
  bounded_vector<RoadGraphPath> graph(scratch.getCells().size(), pool.get());
  for (size_t path = 0; path < graph.size(); ++path) {
    const auto& cells = scratch.getCells()[path];
    const auto& lookup = scratch.getCellsNeighboursLUT()[path];
    const auto& neighbours = scratch.getCellsNeighbours()[path];
    graph[path] = {cells.data(), cells.size(), lookup.data(), lookup.size(), neighbours.data(), neighbours.size()};
  }
  const RoadStep step{graph.data(), graph.size(), context.configuration.topology.nLayers, currentLevel, context.topology.getSurfaceCatalogView()};
  std::atomic<unsigned> errors{0};
  // The neighbour range of every seed; only seeds with neighbours get a job.
  bounded_vector<std::array<int, 3>> ranges(seeds.size(), pool.get()); // path, begin, end
  tbb::parallel_for(size_t{0}, seeds.size(), [&](size_t i) {
    auto& [path, begin, end] = ranges[i];
    unsigned error;
    if constexpr (startCells) {
      path = startPath;
      error = roadRange(frame, step, path, static_cast<int>(i), seeds[i].getLevel(), true, begin, end);
    } else {
      path = seeds[i].cellPathId;
      error = roadRange(frame, step, path, seeds[i].cellId, seeds[i].seed.getLevel(), false, begin, end);
    }
    if (error) {
      errors |= error;
    }
  });
  throwRoadError(errors);
  jobs.clear();
  size_t nTargets = 0;
  for (size_t i = 0; i < seeds.size(); ++i) {
    if (const size_t n = ranges[i][2] - ranges[i][1]) {
      jobs.push_back({{}, nTargets, nTargets + n, {}, startCells, i});
      nTargets += n;
    }
  }
  targets.assign(nTargets, RoadTarget{});
  tbb::parallel_for(size_t{0}, jobs.size(), [&](size_t index) {
    auto& job = jobs[index];
    const auto& [path, begin, end] = ranges[job.source];
    unsigned error = 0;
    if constexpr (startCells) {
      error = makeSeedInput(frame, step.catalog, resolveCell(frame, seeds[job.source]), &job.initial);
    } else {
      job.seed = seeds[job.source].seed;
    }
    for (int neighbour = begin; neighbour < end; ++neighbour) {
      error |= makeRoadTarget(frame, step, graph[path].neighbours[neighbour], targets[job.firstTarget + neighbour - begin]);
    }
    if (error) {
      errors |= error;
    }
  });
  throwRoadError(errors);
}
template void TrackerTraits::buildRoadJobs<Triplet>(IterationContext&, int, int, const bounded_vector<Triplet>&,
                                                    bounded_vector<RoadJob>&, bounded_vector<RoadTarget>&) const;
template void TrackerTraits::buildRoadJobs<RoadSeedEmission>(IterationContext&, int, int, const bounded_vector<RoadSeedEmission>&,
                                                             bounded_vector<RoadJob>&, bounded_vector<RoadTarget>&) const;

RoadSeedSelector TrackerTraits::makeRoadSeedSelector(IterationContext& context, int startLevel) const
{
  const auto& trkParam = context.configuration.parameters;
  // Filter roads by absolute q/pT in parameters[4]'s units, identically for
  // both families. Non-finite values fail the finite-bound comparison.
  constexpr float maxAbsQOverPt = 1.e3f;
  // Missing layers may be allowed, but do not count toward MinTrackLength.
  return {~context.topology.seedingLayers, context.frame.getDetectorConfiguration().getHoleLayers(), trkParam.MaxHoles,
          trkParam.MinTrackLength, maxAbsQOverPt, trkParam.MaxChi2NDF * ((startLevel + 2) * 2 - 5)};
}

void TrackerTraits::findRoads(IterationContext& context, const int iteration)
{
  auto& scratch = context.scratch;
  const auto& pool = scratch.getMemoryPool();
  const float maxChi2 = context.configuration.kernelParameters.maxChi2ClusterAttachment;
  auto& estimator = context.frame.getCapacityEstimator();
  auto firstClusters = makeFirstClusters(context);
  forEachRoadStartLevel(context, iteration, [&](std::span<const CellPathId> starts, int startLevel) {
    const auto selector = makeRoadSeedSelector(context, startLevel);
    bounded_vector<TrackSeed> trackSeeds(pool.get());
    for (const auto start : starts) {
      if (scratch.getCells()[start.value()].empty()) {
        continue;
      }
      bounded_vector<RoadSeedEmission> current(pool.get()), next(pool.get());
      const auto extend = [&](const auto& seeds, int level) {
        bounded_vector<RoadJob> jobs(pool.get());
        bounded_vector<RoadTarget> targets(pool.get());
        const auto key = CapacityEstimator::makeKey(SlabSite::Roads, iteration, CapacityEstimator::makeVariant(startLevel, level), start);
        mTaskArena->execute([&] {
          buildRoadJobs(context, start.value(), level, seeds, jobs, targets);
          GroupedSlabSink<RoadSeedEmission> sink{{.capacity = estimator.capacity(key, seeds.size()), .nThreads = std::max(1, mTaskArena->max_concurrency())}, pool.get()};
          tbb::parallel_for(size_t{0}, jobs.size(), [&](size_t index) {
            const auto& job = jobs[index];
            auto& handle = sink.local();
            handle.beginProducer(static_cast<int>(job.source));
            forEachRoad(job.seed, job.initialize ? &job.initial : nullptr, targets.data(), job.firstTarget, job.endTarget, context.bz, maxChi2,
                        [&](const RoadSeedEmission& emission) { handle.emplace(emission); });
          });
          const auto stats = sink.stats();
          bounded_vector<int> lookup(pool.get());
          sink.finalizeGrouped(seeds.size(), lookup, next);
          estimator.update(key, seeds.size(), stats.requested, stats.capacity, stats.emitted, stats.spilled, stats.overflowed, stats.memoryLimited);
        });
      };
      extend(scratch.getCells()[start.value()], startLevel);
      for (int level = startLevel - 1; level >= 2 && !next.empty(); --level) {
        current.swap(next);
        deepVectorClear(next);
        extend(current, level);
      }
      deepVectorClear(current);
      for (const auto& road : next) {
        if (selector(road.seed)) {
          trackSeeds.push_back(road.seed);
        }
      }
    }
    if (trackSeeds.empty()) {
      return;
    }
    bounded_vector<TrackingCandidate> tracks(pool.get());
    refitTracks(context, trackSeeds, tracks);
    // As o2::its::track::isBetter (longer, then lower chi2); ties keep seed
    // order, as the device sort does. The candidates are large: sort compact
    // keys, then move every candidate once.
    struct SortKey {
      int nClusters;
      float chi2;
      int index;
    };
    bounded_vector<SortKey> keys(tracks.size(), pool.get());
    for (size_t i = 0; i < tracks.size(); ++i) {
      keys[i] = {tracks[i].getNumberOfClusters(), tracks[i].track.chi2, static_cast<int>(i)};
    }
    std::sort(keys.begin(), keys.end(), [](const SortKey& a, const SortKey& b) {
      if (a.nClusters != b.nClusters) {
        return a.nClusters > b.nClusters;
      }
      if (a.chi2 < b.chi2 || b.chi2 < a.chi2) {
        return a.chi2 < b.chi2;
      }
      return a.index < b.index;
    });
    bounded_vector<TrackingCandidate> sorted(pool.get());
    sorted.reserve(tracks.size());
    for (const auto& key : keys) {
      sorted.push_back(std::move(tracks[key.index]));
    }
    acceptSortedTracks(context, iteration, sorted, firstClusters);
  });
}

bounded_vector<bounded_vector<int>> TrackerTraits::makeFirstClusters(IterationContext& context) const
{
  const auto& pool = context.scratch.getMemoryPool();
  const int activeSurfaceCount = context.configuration.topology.nLayers;
  bounded_vector<bounded_vector<int>> firstClusters(activeSurfaceCount, bounded_vector<int>(pool.get()), pool.get());
  firstClusters.resize(activeSurfaceCount);
  return firstClusters;
}

void TrackerTraits::forEachRoadStartLevel(IterationContext& context, int iteration,
                                          const std::function<void(std::span<const CellPathId>, int)>& run) const
{
  const auto& trkParam = context.configuration.parameters;
  // Road starts are the binding's seeding-eligible sparse-plan subsequence.
  // CellPathId values use compact slots; LayerId directly indexes layout-owned
  // layer data.
  const gsl::span<const CellPathId> roadStartCells = context.configuration.topology.roadStartPaths;
  const int cellsPerRoad = context.topology.seedingLayers.count() - 2;
  const auto& componentOffsets = context.configuration.topology.roadStartComponentOffsets;
  if (componentOffsets.empty() || componentOffsets.front() != 0 || componentOffsets.back() != roadStartCells.size()) {
    throw std::invalid_argument{"CA traversal: sparse topology mismatch (iteration " + std::to_string(iteration) + ")"};
  }
  for (size_t component = 0; component + 1 < componentOffsets.size(); ++component) {
    const auto componentRoadStarts = roadStartCells.subspan(componentOffsets[component],
                                                            componentOffsets[component + 1] - componentOffsets[component]);
    for (int startLevel{cellsPerRoad}; startLevel >= trkParam.CellMinimumLevel(); --startLevel) {
      run({componentRoadStarts.data(), componentRoadStarts.size()}, startLevel);
    }
  }
}

void TrackerTraits::acceptSortedTracks(IterationContext& context, int iteration, bounded_vector<TrackingCandidate>& tracks,
                                       bounded_vector<bounded_vector<int>>& firstClusters)
{
  acceptTracks(context, iteration, tracks, firstClusters);
  usedClustersChanged(context);
}

void TrackerTraits::bindHostBuffers(const std::shared_ptr<BoundedMemoryResource>& pool)
{
  if (mHostBuffersPool == pool) {
    return;
  }
  // Rebuild in place: pmr allocators do not propagate on assignment. The old
  // pool is still held here, so the old buffers are released into it.
  deepVectorClear(mRefitJobs, pool.get());
  mHostBuffersPool = pool;
}

RefitParameters TrackerTraits::makeRefitParameters(IterationContext& context) const
{
  const auto& p = context.configuration.parameters;
  return {context.bz, p.MaxChi2ClusterAttachment, p.MaxChi2NDF, p.ShiftRefToCluster, p.RepeatRefitOut};
}

const bounded_vector<RefitJob>& TrackerTraits::prepareRefitJobs(IterationContext& context, std::span<const TrackSeed> seeds)
{
  const auto& p = context.configuration.parameters;
  bindHostBuffers(context.scratch.getMemoryPool());
  auto& jobs = mRefitJobs;
  jobs.resize(seeds.size());
  const auto frame = context.frame.makeView(context.layerGlobalMeasurements);
  mTaskArena->execute([&] {
    tbb::parallel_for(size_t{0}, seeds.size(), [&](size_t i) { prepareRefitJob(frame, seeds[i], p.MinPt.data(), p.MinPt.size(), jobs[i]); });
  });
  return jobs;
}

void TrackerTraits::refitTracks(IterationContext& context, std::span<const TrackSeed> seeds, bounded_vector<TrackingCandidate>& tracks)
{
  const auto catalog = context.topology.getSurfaceCatalogView();
  const auto parameters = makeRefitParameters(context);
  const auto& jobs = prepareRefitJobs(context, seeds);
  bounded_vector<RefitResult> output(jobs.size(), RefitResult{}, context.scratch.getMemoryPool().get());
  mTaskArena->execute([&] {
    tbb::parallel_for(size_t{0}, jobs.size(), [&](size_t i) { output[i] = evaluateRefit(jobs[i], catalog, parameters); });
  });
  tracks.clear();
  tracks.reserve(std::count_if(output.begin(), output.end(), [](const auto& result) { return result.accepted; }));
  for (size_t i = 0; i < output.size(); ++i) {
    if (!output[i].accepted) {
      continue;
    }
    TrackingCandidate track;
    track.seed = seeds[i];
    track.track.innerState = output[i].inner;
    track.track.outerState = output[i].outer;
    track.track.chi2 = output[i].chi2;
    tracks.push_back(std::move(track));
  }
}

void TrackerTraits::acceptTracks(IterationContext& context, int iteration,
                                 bounded_vector<TrackingCandidate>& tracks,
                                 bounded_vector<bounded_vector<int>>& firstClusters)
{
  auto* scratch = &context.scratch;
  auto* mFrame = &context.frame;
  const auto& trkParam = context.configuration.parameters;
  const auto& mLayerGlobalMeasurements = context.layerGlobalMeasurements;
  const int activeSurfaceCount = context.configuration.topology.nLayers;
  reserveGenericTrackPublication(*mFrame, tracks.size(), static_cast<std::size_t>(activeSurfaceCount));
  for (auto& track : tracks) {
    int nShared = 0;
    bool isFirstShared{false};
    int firstLayer{-1}, firstCluster{-1};
    for (int iLayer{0}; iLayer < activeSurfaceCount; ++iLayer) {
      if (track.getClusterIndex(iLayer) == o2::its::constants::UnusedIndex) {
        continue;
      }
      const auto clusterId = mLayerGlobalMeasurements[iLayer][track.getClusterIndex(iLayer)].clusterId;
      bool isShared = mFrame->isClusterUsed(iLayer, clusterId);
      nShared += int(isShared);
      if (firstLayer < 0) {
        firstCluster = track.getClusterIndex(iLayer);
        isFirstShared = isShared && trkParam.AllowSharingFirstCluster && std::find(firstClusters[iLayer].begin(), firstClusters[iLayer].end(), firstCluster) != firstClusters[iLayer].end();
        firstLayer = iLayer;
      }
    }

    /// do not account for the first cluster in the shared clusters number if it is allowed
    if (nShared - int(isFirstShared && trkParam.AllowSharingFirstCluster) > trkParam.SharedMaxClusters) {
      continue;
    }

    bool firstCls{true}, nominalCompatible{true};
    TimeEstBC nominalTS, expandedTS;
    float smallestROFHalf = std::numeric_limits<float>::max();
    for (int iLayer{0}; iLayer < activeSurfaceCount; ++iLayer) {
      if (track.getClusterIndex(iLayer) == o2::its::constants::UnusedIndex) {
        continue;
      }
      smallestROFHalf = std::min(smallestROFHalf, mFrame->getROFTiming(iLayer).mROFLength * 0.5f);
      const auto clusterId = mLayerGlobalMeasurements[iLayer][track.getClusterIndex(iLayer)].clusterId;
      mFrame->markUsedCluster(iLayer, clusterId);
      const int currentROF = mLayerGlobalMeasurements[iLayer][track.getClusterIndex(iLayer)].rof;
      const auto nominalROFTS = mFrame->getROFTiming(iLayer).getROFTimeBounds(currentROF);
      const auto expandedROFTS = mFrame->getROFTiming(iLayer).getROFTimeBounds(currentROF, true);
      if (firstCls) {
        firstCls = false;
        nominalTS = nominalROFTS;
        expandedTS = expandedROFTS;
      } else {
        if (nominalCompatible) {
          if (nominalTS.isCompatible(nominalROFTS)) {
            nominalTS += nominalROFTS;
          } else {
            nominalCompatible = false;
          }
        }
        if (!expandedTS.isCompatible(expandedROFTS)) {
          LOGP(fatal, "TS {}+/-{} are incompatible with {}+/-{}, this should not happen!", expandedROFTS.getTimeStamp(), expandedROFTS.getTimeStampError(), expandedTS.getTimeStamp(), expandedTS.getTimeStampError());
        }
        expandedTS += expandedROFTS;
      }
    }
    track.track.timestamp = (nominalCompatible ? nominalTS : expandedTS).makeSymmetrical();
    // Match the legacy track timestamp, including its uncertainty clamp.
    track.track.timestamp.setTimeStampError(std::min(track.track.timestamp.getTimeStampError(), smallestROFHalf));
    if (!appendGenericTrack(*mFrame, track, mLayerGlobalMeasurements)) {
      LOGP(fatal, "GenericTrack publication failed for an accepted CA track");
    }

    if (trkParam.AllowSharingFirstCluster) {
      firstClusters[firstLayer].push_back(firstCluster);
    }
  }
}

void TrackerTraits::setNThreads(int n, std::shared_ptr<tbb::task_arena>& arena)
{
#if defined(OPTIMISATION_OUTPUT)
  mTaskArena = std::make_shared<tbb::task_arena>(1);
#else
  if (arena == nullptr) {
    mTaskArena = std::make_shared<tbb::task_arena>(std::abs(n));
    LOGP(info, "Setting tracker with {} threads.", n);
  } else {
    mTaskArena = arena;
  }
#endif
}

} // namespace o2::itsmft::tracking
