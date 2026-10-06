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
/// \file TrackerTraits.h
/// \brief Shared CA tracker traits: same ITS-style tracklet/cell/road logic; MFT uses x-y LUT and forward refit
///

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKERTRAITS_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKERTRAITS_H_

#include <array>
#include <optional>
#include <functional>
#include <utility>
#include <vector>

#include <gsl/span>
#include <oneapi/tbb.h>

#include "ITSMFTTracking/TrackingKernels.h"
#include "ITSMFTTracking/TrackSeed.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/GenericTrack.h"
#include "ITSMFTTracking/IterationConfiguration.h"
#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/detail/TrackingKernelParameters.h"
#include "ITSMFTTracking/BoundedAllocator.h"
#include "ITSMFTTracking/ExternalAllocator.h"

namespace o2::itsmft::tracking
{

struct TrackerTestAccess;

struct IterationContext {
  int iteration{-1};
  TimeFrame& frame;
  TimeFrameScratch& scratch;
  TraversalTopologyView topology{};
  const DetectorConfiguration& detectorConfiguration;
  const IterationConfiguration& configuration;
  // Borrows the caller's span sequence and the frame's measurements. Both must
  // outlive this synchronous traversal; loading/resetting the frame invalidates it.
  gsl::span<const gsl::span<const GlobalMeasurement>> layerGlobalMeasurements;
  float bz{0.f};

  IterationContext(int iterationValue, TimeFrame& frameValue, TimeFrameScratch& scratchValue,
                   TraversalTopologyView topologyValue, const IterationConfiguration& configurationValue,
                   gsl::span<const gsl::span<const GlobalMeasurement>> layerGlobalMeasurementsValue,
                   float bzValue)
    : iteration{iterationValue}, frame{frameValue}, scratch{scratchValue}, topology{topologyValue}, detectorConfiguration{frameValue.getDetectorConfiguration()}, configuration{configurationValue}, layerGlobalMeasurements{layerGlobalMeasurementsValue}, bz{bzValue}
  {
  }
};

// Backend implementation of a traversal supplied explicitly by Tracker.
class TrackerTraits
{
 public:
  virtual ~TrackerTraits() = default;
  // The production caller supplies all event and iteration state explicitly.
  void runTraversal(IterationContext& view);

  // Wall time of each stage, in ms, accumulated by runTraversal since the
  // last reset. GPU stages end synchronised (they read back counts), so this
  // includes their device time.
  enum Stage { Tracklets,
               Cells,
               Neighbours,
               Roads,
               NStages };
  const std::array<double, NStages>& getStageTimes() const noexcept { return mStageMs; }
  void resetStageTimes() noexcept { mStageMs.fill(0.); }

  void setValidateRefit(bool value) { mValidateRefit = value; }
  void setValidateRoads(bool value) { mValidateRoads = value; }
  void setValidateNeighbours(bool value) { mValidateNeighbours = value; }
  void setValidateCells(bool value) { mValidateCells = value; }
  void setValidateTracklets(bool value) { mValidateTracklets = value; }
  virtual const char* getName() const noexcept { return "CPU"; }
  virtual bool isGPU() const noexcept { return false; }
  // Called once per timeframe, before its first iteration. No-op on CPU;
  // GPU traits use it to invalidate any device-resident cache that persists
  // data across iterations/vertices within one timeframe.
  virtual void beginTimeframe() {}
  // Whether the stages read the host index tables. A backend that builds
  // them itself from the cluster bins spares the host from filling them.
  virtual bool usesHostIndexTables() const noexcept { return true; }
  // Makes a GPU backend take its device memory from the GPU reconstruction
  // framework, as the legacy TimeFrame::setFrameworkAllocator. No-op on CPU.
  virtual void setFrameworkAllocator(ExternalAllocator*) {}
  void setNThreads(int n, std::shared_ptr<tbb::task_arena>& arena);
  int getNThreads() { return mTaskArena->max_concurrency(); }
  // Runs f in this tracker's task arena, so parallel work inside it respects
  // the configured thread count.
  template <typename F>
  void executeInArena(F&& f)
  {
    if (mTaskArena) {
      mTaskArena->execute(std::forward<F>(f));
    } else {
      std::forward<F>(f)();
    }
  }

 protected:
  // GPU traits override the numerical stages. Host scheduling and shared-hit
  // arbitration remain common to both execution paths.
  virtual void computeLayerTracklets(IterationContext& context, int iteration, int iVertex);
  virtual void computeLayerCells(IterationContext& context, int iteration);
  virtual void findCellsNeighbours(IterationContext& context, int iteration);
  virtual void findRoads(IterationContext& context, int iteration);

  // One edge's tracklet search, pointing at the host frame's ROF tables.
  TrackletSearch prepareTrackletSearch(IterationContext& context, int iVertex, EdgeId edgeId) const;
  // MC labels of the scratch tracklets, if artefact labels are enabled.
  void createTrackletLabels(IterationContext& context);
  bool validateTracklets() const noexcept { return mValidateTracklets; }

  // One path's cell search, without its tracklets.
  CellSearch prepareCellSearch(IterationContext& context, CellPathId cellId) const;
  // After the cell stage: MC labels of the scratch cells, if artefact labels
  // are enabled, then the host tracklets are released.
  void finishCells(IterationContext& context);
  bool validateCells() const noexcept { return mValidateCells; }

  bool validateNeighbours() const noexcept { return mValidateNeighbours; }
  // Runs every target path of the neighbour search with its source paths
  // (cell pointers unset), in order of the targets' outer layer so that their
  // sources' levels are final; targets on one outer layer run concurrently.
  template <typename NCells, typename Run>
  void forEachNeighbourTarget(IterationContext& context, NCells nCells, Run run)
  {
    const auto& topology = context.topology;
    const auto& scheduled = context.configuration.topology.scheduledPaths;
    const float maxChi2 = context.configuration.kernelParameters.maxChi2ClusterAttachment;
    const auto variance = [&](CellPathId path) {
      const float angle = context.scratch.getEdgeMSAngle(topology.getPath(path).second.value());
      return angle * angle;
    };
    const auto outerLayer = [&](CellPathId path) { return topology.getEdge(topology.getPath(path).second).to; };
    executeInArena([&] {
      for (size_t begin = 0, end = 0; begin < scheduled.size(); begin = end) {
        while (end < scheduled.size() && outerLayer(scheduled[end]) == outerLayer(scheduled[begin])) {
          ++end;
        }
        tbb::parallel_for(begin, end, [&](size_t index) {
          const auto target = scheduled[index];
          if (!nCells(target)) {
            return;
          }
          std::vector<NeighbourSearch> sources;
          size_t nSourceCells = 0;
          for (const auto source : scheduled) {
            if (topology.getPath(source).second == topology.getPath(target).first && nCells(source)) {
              sources.push_back({nullptr, source.value(), nullptr, nullptr, 0, {variance(source), variance(target)}, maxChi2});
              nSourceCells += nCells(source);
            }
          }
          if (!sources.empty()) {
            run(target, std::span<NeighbourSearch>{sources}, nSourceCells);
          }
        });
      }
    });
  }
  bool validateRoads() const noexcept { return mValidateRoads; }

  void refitTracks(IterationContext& context, std::span<const TrackSeed> seeds, bounded_vector<TrackingCandidate>& tracks);
  // The host's refit of every seed, as the device refit must reproduce it.
  RefitParameters makeRefitParameters(IterationContext& context) const;
  const bounded_vector<RefitJob>& prepareRefitJobs(IterationContext& context, std::span<const TrackSeed> seeds);
  bool validateRefit() const noexcept { return mValidateRefit; }

  // Called around every traversal pass, the four stages on one selection of
  // vertices; endPass is called also when a stage throws.
  virtual void beginPass(IterationContext&) {}
  virtual void endPass(IterationContext&) {}
  // Called after accepted tracks mark their clusters used.
  virtual void usedClustersChanged(IterationContext&) {}
  // Test support: makes tracklets written into the host scratch visible to a
  // backend that keeps them on the device.
  virtual void adoptHostTracklets(IterationContext&) {}

  // Shared road-stage steps. forEachRoadStartLevel runs every component's
  // start levels in order; acceptRoadSeeds refits a level's seeds, sorts and
  // accepts the tracks.
  bounded_vector<bounded_vector<int>> makeFirstClusters(IterationContext& context) const;
  void forEachRoadStartLevel(IterationContext& context, int iteration,
                             const std::function<void(std::span<const CellPathId>, int)>& run) const;
  RoadSeedSelector makeRoadSeedSelector(IterationContext& context, int startLevel) const;
  // Accepts tracks sorted longest first, then by chi2 (ties in seed order).
  void acceptSortedTracks(IterationContext& context, int iteration, bounded_vector<TrackingCandidate>& tracks,
                          bounded_vector<bounded_vector<int>>& firstClusters);
  // The jobs and targets of one road extension of startPath's cells or of
  // the previous extension's seeds.
  template <typename InputSeed>
  void buildRoadJobs(IterationContext& context, int startPath, int currentLevel, const bounded_vector<InputSeed>& seeds,
                     bounded_vector<RoadJob>& jobs, bounded_vector<RoadTarget>& targets) const;

 private:
  friend struct TrackerTestAccess;

  void acceptTracks(IterationContext& context, int iteration,
                    bounded_vector<TrackingCandidate>& tracks,
                    bounded_vector<bounded_vector<int>>& firstClusters);

  bool buildTrackSeed(IterationContext& context, int cellPathId,
                      const Triplet& cell, TrackSeed& output) const;


  bool mValidateTracklets{false};
  bool mValidateCells{false};
  bool mValidateNeighbours{false};
  bool mValidateRefit{false};
  bool mValidateRoads{false};
  std::array<double, NStages> mStageMs{};
  std::shared_ptr<tbb::task_arena> mTaskArena;
  // Host staging reused across calls: allocating these afresh re-faults tens
  // of MB every time. The pool is held so the buffers (declared after it,
  // destroyed first) never outlive their allocator.
  void bindHostBuffers(const std::shared_ptr<BoundedMemoryResource>& pool);
  std::shared_ptr<BoundedMemoryResource> mHostBuffersPool;
  bounded_vector<RefitJob> mRefitJobs;
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKERTRAITS_H_ */
