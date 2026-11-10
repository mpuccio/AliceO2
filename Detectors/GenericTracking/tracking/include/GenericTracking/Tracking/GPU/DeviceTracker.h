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

/// \file DeviceTracker.h
/// \brief A timeframe's tracking on a device: every slice of every
///        iteration, each iteration's selection, refit and follow, and
///        the tracks.
///
/// A host interface without a vendor type, implemented once per backend
/// library and detector (GPU/cuda/DeviceTracker.cuh). It owns a
/// device arena, a stream, the iterations' parameters, the device copy of
/// the current timeframe with its plan and used mask, and the tracks. It
/// runs a TimeframePlan, the host tracker's own (SlicePlan.h), so both run
/// on the same numbers; `run` drives it alone, a comparison call by call.

#ifndef GTRK_TRACKING_GPU_DEVICETRACKER_H_
#define GTRK_TRACKING_GPU_DEVICETRACKER_H_

#include "GenericTracking/Base/Columns.h"
#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceStageOutput.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/SlicePlan.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/TrackRefit.h"
#include "GenericTracking/Tracking/VertexStage.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/Neighbour.h"
#include "GenericTracking/Tracking/Triplet.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <span>
#include <vector>

namespace o2::gtrk
{

/// What a slice leaves on the device, downloaded for a comparison with the
/// host tracker's. The stage outputs are the
/// host types, so the existing comparisons apply.
template <DetectorLayout Det>
struct DeviceSliceDownload {
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;

  explicit DeviceSliceDownload(const std::shared_ptr<BoundedMemoryResource>& pool)
  {
    tracklets.setMemoryPool(pool);
    triplets.setMemoryPool(pool);
    neighbours.setMemoryPool(pool);
    records.setMemoryPool(pool);
  }

  StageOutput<TrackletT> tracklets;
  StageOutput<TripletT> triplets; ///< with their final levels
  StageOutput<NeighbourT> neighbours;
  std::vector<std::uint32_t> chained;
  std::vector<std::uint32_t> successorAbove;
  std::vector<std::uint32_t> roads;
  std::vector<std::uint32_t> starts;
  std::vector<std::uint32_t> roadOffsets;
  StageOutput<SelectionRecordOf<kMaxClusters>::template RecordT> records;
};

/// What an iteration leaves on the device, downloaded for a comparison: the
/// records of every slice, gathered in slice order;
/// the selection's order and choice, the follow's proposals, the
/// iteration's tracks and the used mask after it.
template <DetectorLayout Det>
struct DeviceIterationDownload {
 private:
  // Own the allocator until after the columns have been destroyed.
  std::shared_ptr<BoundedMemoryResource> mPool;

 public:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using Records = BoundedColumns<SelectionRecordOf<kMaxClusters>::template RecordT>;
  using Tracks = typename TrackRefit<Det>::Tracks;
  using Proposal = typename Follow<Det>::Proposal;

  explicit DeviceIterationDownload(const std::shared_ptr<BoundedMemoryResource>& pool) : mPool(pool)
  {
    records.setMemoryResource(mPool.get());
    tracks.setMemoryResource(mPool.get());
  }

  Records records;
  std::vector<std::uint32_t> ranked;   ///< the unique records, best first, as indices into `records`
  std::vector<std::uint32_t> selected; ///< the accepted ones, best first
  std::vector<Proposal> proposals;     ///< one per track after the follow's refit; none without the follow
  Tracks tracks;                       ///< the iteration's, refitted and extended
  std::vector<std::uint8_t> used;      ///< the mask after the iteration
};

/// The tracker on a device, for a detector and its material (the tests' UniformMaterial,
/// o2::base::MatLayerCylSet). Its calls follow a lifecycle of nested scopes:
///
///   attach            the device memory to work in, the tracker's only one
///                     until detach: nothing of an earlier block stays
///   configureRun      the run: parameters, geometry, material, propagator
///   beginTimeframe    a timeframe: measurements, offsets, the plan, used mask
///   beginIteration    an iteration: its slices' records kept as they come
///   runSlice          a slice of the plan, as many as the iteration has
///   endIteration      the iteration's records gathered and selected, the
///                     tracks refitted, extended and kept; the next may begin
///   downloadTracks    every iteration's tracks so far
///   endTimeframe      the timeframe closed: the host frame and plan no
///                     longer read
///   detach            the queued work finished and the block let go: the
///                     run and the timeframe gone with it
///
/// The tracking device attaches and configures every timeframe: its block may
/// be another one each time (the GPU reconstruction's pool after the TPC).
///
/// A timeframe as the tracking device runs it, the vertexer's pass first if
/// the run has one: `findTimeframeVertices`, then `trackTimeframe`.
///
/// A call out of order throws std::logic_error. A failure inside a call
/// closes what that call opened; a failed slice fails its iteration, whose
/// later slices are refused and whose endIteration returns the failure.
template <DetectorLayout Det, typename Material>
class DeviceTracker
{
 public:
  /// What the last slice did.
  struct SliceReport {
    std::uint32_t tracklets = 0U;
    std::uint32_t permutation = 0U; ///< first tracklets the triplet stage worked
    std::uint32_t triplets = 0U;
    std::uint32_t neighbours = 0U;
    std::uint32_t chained = 0U;
    std::uint32_t roads = 0U;
    std::uint32_t records = 0U;
    std::uint32_t readbacks = 0U;
    bool stagingPredictionMissed = false; ///< by any of the stages
    /// The flat runs: the tracklet stage's
    /// where the backend runs it flat (its candidates past `cheap`), the
    /// triplet stage's (its window candidates over the first tracklets).
    std::uint32_t trackletSources = 0U; ///< the sources with a query the flat tracklet run listed
    FlatRunReport trackletRun;
    FlatRunReport tripletRun;
  };

  /// What the last iteration's end did: its records, its selection, refit
  /// and follow.
  struct IterationReport {
    std::size_t records = 0U;
    std::uint32_t duplicates = 0U;
    std::uint32_t mismatches = 0U; ///< duplicates that differ: must be none -- a record found by more than one slice is the same record in each
    std::uint32_t rounds = 0U;
    std::uint32_t tracks = 0U;
    std::uint32_t readbacks = 0U;    ///< of the selection and the length cut
    std::size_t arenaPeakBytes = 0U; ///< the whole arena, both ends, at its largest from beginIteration to endIteration
  };
  using Tracks = typename TrackRefit<Det>::Tracks;

  DeviceTracker() = default;
  virtual ~DeviceTracker() = default;
  DeviceTracker(const DeviceTracker&) = delete;
  DeviceTracker& operator=(const DeviceTracker&) = delete;
  DeviceTracker(DeviceTracker&&) = delete;
  DeviceTracker& operator=(DeviceTracker&&) = delete;

  /// Work in \a block: what the tracker held goes (detach), configureRun is
  /// due before the next timeframe.
  virtual void attach(DeviceBlock block) = 0;
  /// Finish the queued work, then close the run and every scope inside it and
  /// let go of the block: it may be reused once this returns.
  virtual void detach() noexcept = 0;

  /// Rebuild the run: every open scope closed, the arena reset, then every
  /// iteration's parameters, the geometry and the material uploaded and the
  /// propagator kept. kOutOfMemory, and no run configured, if the arena
  /// refuses. A FlatObject material's host image is kept from one call to
  /// the next and made again only for another object: the same object is
  /// taken to be the same material.
  [[nodiscard]] virtual StageStatus configureRun(std::span<const IterationParameters<Det>> iterations, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator) = 0;

  /// Upload \a frame -- its measurement columns and grid offsets -- the
  /// device part of \a plan and a used mask with every row free, replacing
  /// the previous timeframe, its tracks and any open iteration. \a plan must
  /// be of \a frame and of the configured run, std::invalid_argument before
  /// anything changes otherwise; both must outlive the timeframe.
  [[nodiscard]] virtual StageStatus beginTimeframe(const TimeFrame<Det>& frame, const TimeframePlan<Det>& plan) = 0;

  /// The staging every staged stage starts its next run with (the flat runs
  /// have none); each then keeps what its runs needed. 2^20
  /// records by default.
  virtual void setStagingCapacity(std::size_t records) = 0;
  /// The flat runs' chunks held to \a pairs at most (0: the arena's free
  /// bytes decide): a test hook, so that chunks
  /// cut inside sources and first tracklets on a fixture; never a tracker's.
  virtual void setFlatChunkPairs(std::size_t pairs) = 0;

  /// Open iteration \a iteration of the configured ones; the previous
  /// iteration's gathered records go.
  virtual void beginIteration(std::size_t iteration) = 0;

  /// Run slice \a slice of the open iteration's plan, tracklets to records,
  /// as the host tracker runs it. With \a download, the slice's outputs come
  /// down into it.
  [[nodiscard]] virtual StageStatus runSlice(std::size_t slice, DeviceSliceDownload<Det>* download) = 0;

  /// Close the open iteration: its records gathered, selected, the selected
  /// ones refitted and extended, the tracks kept for the timeframe and the
  /// used mask committed; with \a download brought down. The first failure
  /// of its slices or of its end, if any: a failed iteration leaves the used
  /// mask and the tracks as they were before it.
  [[nodiscard]] virtual StageStatus endIteration(DeviceIterationDownload<Det>* download) = 0;

  /// Every iteration's tracks so far, in iteration order, each iteration's
  /// best first, into \a tracks; not with an iteration open.
  virtual void downloadTracks(Tracks& tracks) = 0;

  /// The vertexing pass: iteration \a pass
  /// of the configured run -- the pass, configured after the tracking's --
  /// slice by slice as its plan has them, tracklets and triplets, then the
  /// vertexer on them with \a settings for the frames the slice's core owns
  /// (`OwnedFrames`), \a beamSpot giving the vertices' transverse spread.
  /// The vertices, sorted by
  /// time and each frame's by z, and their contributors into \a result, as
  /// `VertexingPass::find` gives them; its counts too, and optionally their
  /// claimed triplets in \a claims; no beam estimate. Not with an
  /// iteration open; the used mask and the tracks stay as they are.
  [[nodiscard]] virtual StageStatus findVertices(std::size_t pass, const VertexerSettings& settings, const VertexConstraint& beamSpot, VertexerResult& result, std::vector<TripletClaim>* claims = nullptr) = 0;

  /// Replace the open timeframe's plan with \a plan, of the same frame and
  /// the same run (std::invalid_argument before anything changes
  /// otherwise): its device part uploaded -- the slices, the vertex windows,
  /// the seeds -- the measurements, the used mask and the tracks kept. The
  /// tracking's plan with the vertexer's list, after `findVertices`. Not with
  /// an iteration open; \a plan must outlive the timeframe.
  [[nodiscard]] virtual StageStatus replacePlan(const TimeframePlan<Det>& plan) = 0;

  /// Close the open timeframe, and any iteration in it: from here until the
  /// next beginTimeframe no call reads its host frame or plan. The stages
  /// keep the staging they learnt.
  virtual void endTimeframe() noexcept = 0;

  /// What a timeframe's calls did, each part on the host's clock.
  struct TimeframeReport {
    double uploadSeconds = 0.;     ///< the frame and its first plan (beginTimeframe)
    double vertexSeconds = 0.;     ///< the vertexer's pass (findVertices)
    double planUploadSeconds = 0.; ///< the tracking's plan with the vertexer's list (replacePlan)
    std::vector<double> iterationSeconds;
    std::vector<IterationReport> iterations;
    std::vector<std::vector<SliceReport>> slices; ///< per iteration, per slice run
    double downloadSeconds = 0.;
  };

  /// A timeframe's vertices: \a frame and \a plan -- the run's every
  /// iteration, the pass last, with an empty vertex list -- uploaded, then
  /// the pass, iteration \a pass, into \a result (findVertices). The
  /// timeframe stays open for trackTimeframe; \a frame and \a plan must
  /// outlive it.
  [[nodiscard]] StageStatus findTimeframeVertices(const TimeFrame<Det>& frame, const TimeframePlan<Det>& plan, std::size_t pass, const VertexerSettings& settings, const VertexConstraint& beamSpot, VertexerResult& result, TimeframeReport& report, std::vector<TripletClaim>* claims = nullptr)
  {
    report = {};
    auto start = Clock::now();
    StageStatus status = beginTimeframe(frame, plan);
    report.uploadSeconds = since(start);
    if (status != StageStatus::kComplete) {
      return status;
    }
    start = Clock::now();
    status = findVertices(pass, settings, beamSpot, result, claims);
    report.vertexSeconds = since(start);
    return status;
  }

  /// A timeframe's tracks: its iterations [0, \a iterations) as \a plan
  /// has them, then every track into \a tracks. With \a replaced, \a plan
  /// -- the tracking's, with the vertexer's list -- replaces the open
  /// timeframe's (findTimeframeVertices'); else \a frame and \a plan begin
  /// the timeframe. Stops at the first failure and returns it, as
  /// Tracker::run does; tracks the host cannot hold are kOutOfMemory. \a
  /// frame and \a plan must outlive the timeframe.
  [[nodiscard]] StageStatus trackTimeframe(const TimeFrame<Det>& frame, const TimeframePlan<Det>& plan, std::size_t iterations, bool replaced, Tracks& tracks, TimeframeReport& report)
  {
    if (!replaced) {
      report = {};
    }
    auto start = Clock::now();
    StageStatus status = replaced ? replacePlan(plan) : beginTimeframe(frame, plan);
    (replaced ? report.planUploadSeconds : report.uploadSeconds) = since(start);
    for (std::size_t iteration = 0; iteration < iterations && status == StageStatus::kComplete; ++iteration) {
      start = Clock::now();
      beginIteration(iteration);
      report.slices.emplace_back();
      StageStatus slices = StageStatus::kComplete;
      for (std::size_t slice = 0; slice < plan.slices.size() && slices == StageStatus::kComplete; ++slice) {
        slices = runSlice(slice, nullptr);
        if (slices == StageStatus::kComplete) {
          report.slices.back().push_back(lastSlice());
        }
      }
      status = endIteration(nullptr);
      report.iterationSeconds.push_back(since(start));
      report.iterations.push_back(lastIteration());
    }
    if (status != StageStatus::kComplete) {
      return status;
    }
    start = Clock::now();
    try {
      downloadTracks(tracks);
    } catch (const std::bad_alloc&) {
      return StageStatus::kOutOfMemory;
    }
    report.downloadSeconds = since(start);
    return StageStatus::kComplete;
  }

  [[nodiscard]] virtual const SliceReport& lastSlice() const = 0;
  [[nodiscard]] virtual const IterationReport& lastIteration() const = 0;
  /// The arena's size, what each end holds now, and the largest part of it
  /// used so far.
  [[nodiscard]] virtual std::size_t arenaBytes() const = 0;
  [[nodiscard]] virtual std::size_t arenaUsedBytes(ArenaEnd end) const = 0;
  [[nodiscard]] virtual std::size_t arenaPeakBytes() const = 0;

 private:
  using Clock = std::chrono::steady_clock;
  [[nodiscard]] static double since(Clock::time_point from) { return std::chrono::duration<double>(Clock::now() - from).count(); }
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_DEVICETRACKER_H_
