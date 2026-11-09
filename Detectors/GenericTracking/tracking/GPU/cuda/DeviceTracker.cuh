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

/// \file DeviceTracker.cuh
/// \brief DeviceTracker<Det, Material> on one backend: a timeframe's slices,
///        iterations and tracks.
///
/// What lives where in the arena:
///
///   run        the iterations' parameters, the geometry, the material, until
///              configureRun again or detach; the arena itself, over the
///              attached block
///   timeframe  base end: the measurement columns, the grid offsets, the used
///              mask, the TimeFrame over them and the plan's device part,
///              until the next beginTimeframe
///   tracks     scratch end, at its bottom: every iteration's tracks, one
///              chunk per iteration, until the next beginTimeframe
///   chunks     scratch end: the open iteration's records, one chunk per
///              slice, until endIteration has gathered them
///   slice      base end: the slice's outputs, until the next runSlice or
///              endIteration
///   iteration  base end, inside endIteration: the gathered records, the
///              selection's arrays, the working mask, the iteration's tracks
///              and proposals
///
/// Staging, scan and sort storage are on the scratch end, released when the
/// pass that took them returns. Included by the .cu that instantiates it for
/// a detector, in the backend's namespace, like GpuStage.cuh.

#ifndef GTRK_TRACKING_GPU_DEVICETRACKER_CUH_
#define GTRK_TRACKING_GPU_DEVICETRACKER_CUH_

#include "DeviceAlgorithms.cuh"
#include "GpuStage.cuh"
#include "RoadGlue.cuh"
#include "Platform.cuh"
#include "PathGlue.cuh"
#include "SelectionGlue.cuh"
#include "SliceGlue.cuh"
#include "VertexGlue.cuh"

#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStageOutput.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/GPU/DeviceTracker.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/NeighbourStage.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/RoadCount.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/SlicePlan.h"
#include "GenericTracking/Tracking/TrackRefit.h"
#include "GenericTracking/Tracking/TrackletStage.h"
#include "FlatStage.cuh"
#include "GenericTracking/Tracking/TripletStage.h"

#include "FlatObject.h"

#include <algorithm>
#include <array>
#include <string>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>
#include <type_traits>
#include <utility>
#include <span>
#include <stdexcept>
#include <vector>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
{

namespace tracker_detail
{

/// The lanes that follow one track together, on every iteration: a step's
/// window holds many
/// candidates, each an update, and the kink step confirms every row it admits
/// by the next step's search -- long windows a group shares. 16 on both
/// backends, for warps of 32 and waves of 64 alike.
inline constexpr std::uint32_t kFollowLanes = 16U;

/// \a stage's follow of \a count tracks: `kFollowLanes` lanes in step a
/// track (`FollowStage::together`). An iteration that decides the end again
/// (`followEnd`) has no kink allowance (the settings refuse both) and is
/// followed one thread a track: its stage has no group kernel.
template <typename Stage>
void follow(DeviceRuntime& runtime, DeviceStream& stream, std::uint32_t count, const Stage& stage)
{
  if constexpr (Stage::kDecidesEnd) {
    algorithms::forEach(runtime, stream, count, stage);
  } else {
    algorithms::forEachGroup<kFollowLanes>(runtime, stream, count, stage);
  }
}

/// The device's launches of the follow's passes (`follow::runSequence`) on
/// \a stream: a pass one thread a track, the proposals by `follow`, the bids
/// zeroed by a fill.
struct FollowLaunch {
  DeviceRuntime& runtime;
  DeviceStream& stream;

  template <typename Body>
  void pass(profile::Phase /*phase*/, std::uint32_t count, const Body& body) const
  {
    algorithms::forEach(runtime, stream, count, body);
  }

  template <typename Stage>
  void propose(profile::Phase /*phase*/, std::uint32_t count, const Stage& stage) const
  {
    tracker_detail::follow(runtime, stream, count, stage);
  }

  void clearBids(std::uint32_t* bids, std::size_t rows) const { runtime.fill(bids, 0U, rows * sizeof(std::uint32_t), stream); }
};

/// Whether a source has a query (`TrackletStage::hasQuery`), for the
/// tracklet launch's list.
template <typename Stage>
struct HasQuery {
  Stage stage;
  GPUdi() bool operator()(std::uint32_t source) const { return stage.hasQuery(source); }
};

/// The tracks holding at least \a minimum clusters: the length after the
/// follow.
template <DetectorLayout Det>
struct KeptAfterFollow {
  const IterationParameters<Det>* parameters = nullptr;
  const std::uint8_t* count = nullptr;
  const BoundTrackState* state = nullptr;
  GPUdi() bool operator()(std::uint32_t track) const { return keptAfterFollow(*parameters, count[track], state[track]); }
};

/// The kept tracks' rows in \a mask as \a after has them: the mask before the
/// selection plus what the kept tracks hold, as the host's `dropShortTracks`
/// rebuilds it. Two kept tracks sharing a row store the same byte.
template <std::size_t K>
struct RestoreRows {
  const ClusterList<K>* clusters = nullptr;
  const std::uint8_t* count = nullptr;
  const std::uint8_t* after = nullptr;
  std::uint8_t* mask = nullptr;
  GPUdi() void operator()(std::uint32_t track) const
  {
    for (std::uint32_t position = 0; position < count[track]; ++position) {
      const auto row = clusters[track].index[position];
      mask[row] = after[row];
    }
  }
};

/// The tracklet stage as a flat run:
/// the walk per source is `TrackletStage::forEachCandidate`, the contexts
/// its source's (`SourceContext`, the vertex cache, `BetweenSource`) and each
/// query's (`QueryContext`, the vertex intervals), the test `testCandidate`,
/// the record `record`. A pair's key is the source's slot in its chunk (24
/// bits) and its query's position among the source's out-links (8).
template <typename Stage>
struct FlatTracklets {
  using SourceContext = typename Stage::SourceContext;
  using QueryContext = typename Stage::QueryContext;
  using CompatibleVertex = typename Stage::CompatibleVertex;
  using Between = typename Stage::Between;
  static constexpr std::uint32_t kVertexCapacity = Stage::kVertexCapacity;
  static constexpr std::uint32_t kQueryBits = 8U;
  /// The lanes that walk one source together in the fill: one thread a
  /// source, or a group of lanes (`WarpGroup`).
  static constexpr std::uint32_t kFillLanes = gpu::kFillLanes;
  struct Source {
    SourceContext context;
    CompatibleVertex vertices[kVertexCapacity];
    Between between;
  };
  struct Query {
    QueryContext context;
    SecondRange intervals[kVertexCapacity];
    bool hasIntervals = false;
  };
  /// A fill's bookkeeping over one source's walk (`note`).
  struct FillState {
    bool stored = false;
    std::uint32_t query = 0U;
    std::size_t lastLink = 0U;
  };
  Stage stage;
  std::uint32_t queriesPerSource = 1U; ///< the most out-links of a layer
  Source* sources = nullptr;           ///< per chunk element
  Query* queries = nullptr;            ///< per chunk element, `queriesPerSource` each
  std::uint32_t sourceBegin = 0U;

  [[nodiscard]] GPUdi() std::uint32_t count(std::uint32_t element) const { return stage.candidateCount(element); }
  template <typename Sink>
  GPUdi() void fill(std::uint32_t element, std::uint32_t position, Sink& sink) const
  {
    const std::uint32_t slot = position - sourceBegin;
    FillState state;
    const auto one = [&](const SourceContext& source, const CompatibleVertex* cachedVertices, const Between& between, const QueryContext& context, const SecondRange* intervals, const MeasurementView&, std::uint32_t targetRow) { sink.put(targetRow, note(state, slot, true, source, cachedVertices, between, context, intervals)); };
    stage.forEachCandidate(element, one);
  }
  /// `fill` by the lanes of \a group in step (`forEachCandidate`): a step's
  /// candidates in lane order at the sink's next positions, so the list is
  /// the serial walk's. The bookkeeping on the step's ballot, alike on every
  /// lane; lane 0 stores.
  template <typename Sink, typename Group>
  GPUdi() void fill(std::uint32_t element, std::uint32_t position, Sink& sink, const Group& group) const
  {
    const std::uint32_t slot = position - sourceBegin;
    const std::uint32_t lane = group.lane();
    FillState state;
    const auto step = [&](const SourceContext& source, const CompatibleVertex* cachedVertices, const Between& between, const QueryContext& context, const SecondRange* intervals, const MeasurementView&, std::uint32_t targetRow, bool has) {
      const std::uint64_t mask = group.ballot(has);
      if (mask == 0U) {
        return;
      }
      sink.putStep(mask, lane, has, targetRow, note(state, slot, lane == 0U, source, cachedVertices, between, context, intervals));
    };
    stage.forEachCandidate(element, step, group);
  }
  /// A candidate's key, its source and query noted in \a state: the source's
  /// context stored at its first candidate, a query's at its first; a
  /// query's number counts the queries with candidates before it. \a store
  /// whether this thread writes the contexts.
  [[nodiscard]] GPUdi() std::uint32_t note(FillState& state, std::uint32_t slot, bool store, const SourceContext& source, const CompatibleVertex* cachedVertices, const Between& between, const QueryContext& context, const SecondRange* intervals) const
  {
    if (!state.stored) {
      if (store) {
        sources[slot].context = source;
        for (std::uint32_t i = 0; i < source.cachedCount; ++i) {
          sources[slot].vertices[i] = cachedVertices[i];
        }
        sources[slot].between = between;
        keep(slot, 0U, context, intervals);
      }
      state.stored = true;
      state.lastLink = context.linkIndex;
      state.query = 0U;
    } else if (context.linkIndex != state.lastLink) {
      state.lastLink = context.linkIndex;
      ++state.query;
      if (store) {
        keep(slot, state.query, context, intervals);
      }
    }
    return (slot << kQueryBits) | state.query;
  }
  GPUdi() void keep(std::uint32_t slot, std::uint32_t query, const QueryContext& context, const SecondRange* intervals) const
  {
    Query& out = queries[(std::size_t{slot} * queriesPerSource) + query];
    out.context = context;
    out.hasIntervals = intervals != nullptr;
    if (intervals != nullptr) {
      for (std::uint32_t i = 0; i < kVertexCapacity; ++i) {
        out.intervals[i] = intervals[i];
      }
    }
  }
  [[nodiscard]] GPUdi() std::uint32_t sourceOf(const flat::CandidatePair& pair) const { return (pair.key >> kQueryBits) + sourceBegin; }
  [[nodiscard]] GPUdi() bool test(std::uint32_t, const flat::CandidatePair& pair) const
  {
    const Source& source = sources[pair.key >> kQueryBits];
    const Query& query = queries[(std::size_t{pair.key >> kQueryBits} * queriesPerSource) + (pair.key & ((1U << kQueryBits) - 1U))];
    const MeasurementView columns = stage.columns();
    const PairPoint target = pairPointAt(columns, pair.row);
    return stage.testCandidate(source.context, source.vertices, source.between, query.context, query.hasIntervals ? query.intervals : nullptr, columns, pair.row, target);
  }
  [[nodiscard]] GPUdi() Tracklet emit(std::uint32_t, const flat::CandidatePair& pair, std::uint32_t) const
  {
    const Source& source = sources[pair.key >> kQueryBits];
    const Query& query = queries[(std::size_t{pair.key >> kQueryBits} * queriesPerSource) + (pair.key & ((1U << kQueryBits) - 1U))];
    const PairPoint target = pairPointAt(stage.columns(), pair.row);
    return stage.record(source.context.at.row, pair.row, query.context.linkIndex, source.context.point, target);
  }
  [[nodiscard]] bool beginChunk(DeviceArena& arena, std::uint32_t begin, std::uint32_t count, std::uint32_t)
  {
    if (count > (std::uint32_t{1U} << (32U - kQueryBits))) {
      return false; // more sources in the chunk than the key numbers: the run fails for memory, which halves the chunk
    }
    sourceBegin = begin;
    DeviceSpan<Source> sourceStorage;
    DeviceSpan<Query> queryStorage;
    if (!sourceStorage.allocate(arena, std::max<std::size_t>(count, 1U), ArenaEnd::kScratch) || !queryStorage.allocate(arena, std::max<std::size_t>(std::size_t{count} * queriesPerSource, 1U), ArenaEnd::kScratch)) {
      return false;
    }
    sources = sourceStorage.data();
    queries = queryStorage.data();
    return true;
  }
  [[nodiscard]] std::size_t contextBytes(std::uint32_t count) const { return std::size_t{count} * (sizeof(Source) + (std::size_t{queriesPerSource} * sizeof(Query))); }
  [[nodiscard]] std::size_t pairBytes() const { return 0U; }
};

/// The triplet stage as a flat run:
/// the walk per first tracklet is `TripletStage::prepareFirst` and
/// `forEachWindowCandidate`, the context its `FirstContext`, the test
/// `testCandidate` with its fit kept per pair for `record`; with \a Shape,
/// the shape-testing stage.
template <DetectorLayout Det, bool Shape = false>
struct FlatTriplets {
  using Stage = TripletStage<Det, Shape>;
  using Context = typename Stage::FirstContext;
  Stage stage;
  Context* contexts = nullptr;      ///< per chunk element
  TripletFitResult* fits = nullptr; ///< per chunk pair, where the test passed
  std::uint8_t* kinked = nullptr;
  std::uint32_t sourceBegin = 0U;

  [[nodiscard]] GPUdi() std::uint32_t count(std::uint32_t element) const { return stage.candidateCount(element); }
  template <typename Sink>
  GPUdi() void fill(std::uint32_t element, std::uint32_t position, Sink& sink) const
  {
    Context context;
    if (!stage.prepareFirst(element, context)) {
      return; // no candidates: the context is never read
    }
    contexts[position - sourceBegin] = context;
    const auto one = [&](std::uint32_t candidateIndex) { sink.put(candidateIndex, position); };
    stage.forEachWindowCandidate(context, one);
  }
  [[nodiscard]] GPUdi() std::uint32_t sourceOf(const flat::CandidatePair& pair) const { return pair.key; }
  [[nodiscard]] GPUdi() bool test(std::uint32_t k, const flat::CandidatePair& pair) const
  {
    TripletFitResult fit;
    bool isKinked = false;
    if (!stage.testCandidate(contexts[pair.key - sourceBegin], pair.row, fit, isKinked)) {
      return false;
    }
    fits[k] = fit;
    kinked[k] = isKinked ? std::uint8_t{1} : std::uint8_t{0};
    return true;
  }
  [[nodiscard]] GPUdi() Triplet emit(std::uint32_t k, const flat::CandidatePair& pair, std::uint32_t) const { return stage.record(contexts[pair.key - sourceBegin], pair.row, fits[k], kinked[k] != 0U); }
  [[nodiscard]] bool beginChunk(DeviceArena& arena, std::uint32_t begin, std::uint32_t sources, std::uint32_t pairs)
  {
    sourceBegin = begin;
    DeviceSpan<Context> contextStorage;
    DeviceSpan<TripletFitResult> fitStorage;
    DeviceSpan<std::uint8_t> kinkedStorage;
    if (!contextStorage.allocate(arena, std::max<std::size_t>(sources, 1U), ArenaEnd::kScratch) || !fitStorage.allocate(arena, std::max<std::size_t>(pairs, 1U), ArenaEnd::kScratch) || !kinkedStorage.allocate(arena, std::max<std::size_t>(pairs, 1U), ArenaEnd::kScratch)) {
      return false;
    }
    contexts = contextStorage.data();
    fits = fitStorage.data();
    kinked = kinkedStorage.data();
    return true;
  }
  [[nodiscard]] std::size_t contextBytes(std::uint32_t sources) const { return std::size_t{sources} * sizeof(Context); }
  [[nodiscard]] std::size_t pairBytes() const { return sizeof(TripletFitResult) + 1U; }
};

} // namespace tracker_detail

/// Where a material lives on the device: an O2 FlatObject -- the
/// MatLayerCylSet -- with its buffer, anything else as a value. Decided by
/// the base class, not by trivial copyability: O2 compiles a FlatObject's
/// copy constructor and destructor only outside GPU translation units, so in
/// this one a MatLayerCylSet looks trivially copyable, and copied as bytes it
/// keeps its host buffer pointers.
template <typename Material, bool Plain = !std::is_base_of_v<o2::gpu::FlatObject, Material>>
struct DeviceMaterial {
  static_assert(std::is_trivially_copyable_v<Material>, "a material that is not a FlatObject uploads as bytes");
  DeviceObject<Material> storage;
  [[nodiscard]] bool upload(DeviceArena& arena, const Material& material, DeviceStream& stream)
  {
    if (!storage.allocate(arena)) {
      return false;
    }
    uploadObject(material, storage, stream);
    return true;
  }
  void release() noexcept { storage = {}; }
  [[nodiscard]] const Material* data() const { return storage.data(); }
};

template <typename Material>
struct DeviceMaterial<Material, false> {
  DeviceFlatObject<Material> storage;
  [[nodiscard]] bool upload(DeviceArena& arena, const Material& material, DeviceStream& stream) { return storage.upload(arena, material, stream); }
  /// The device copy let go, the host image kept for the next run's upload.
  void release() noexcept { storage.release(); }
  [[nodiscard]] const Material* data() const { return storage.data(); }
};

template <DetectorLayout Det, typename Material>
class DeviceTrackerImpl final : public DeviceTracker<Det, Material>
{
 public:
  using Report = typename DeviceTracker<Det, Material>::SliceReport;
  using Frame = TimeFrame<Det>;
  using Index = std::uint32_t;
  static constexpr std::size_t kLayerCount = LayerRowMap<Det>::kLayerCount;
  static constexpr std::size_t kInitialStaging = std::size_t{1} << 20U;
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using Records = DeviceStageOutput<SelectionRecordOf<kMaxClusters>::template RecordT>;
  using IterationReport = typename DeviceTracker<Det, Material>::IterationReport;
  using Tracks = typename DeviceTracker<Det, Material>::Tracks;
  using TrackRecord = TrackOf<kMaxClusters>;
  using Proposal = FollowProposal<kMaxClusters>;
  static_assert(std::is_trivially_copyable_v<TrackingGeometry<Det>>, "the geometry uploads as bytes");

  /// The arena on \a block.
  DeviceTrackerImpl(DeviceRuntime& runtime, DeviceBlock block)
    : mRuntime(&runtime),
      mArena(runtime, block),
      mStream(runtime.createStream()),
      mTrackletStage(runtime, mArena, mStream),
      mNeighbourStage(runtime, mArena, mStream),
      mRecordStage(runtime, mArena, mStream)
  {
    setStagingCapacity(kInitialStaging);
  }

  ~DeviceTrackerImpl() final
  {
    // Work still queued reads the arena: finish it before the arena goes.
    mStream.synchronize();
    closeFrom(Scope::kRun);
  }

  void attach(DeviceBlock block) final
  {
    detach();
    mArena.rebind(block);
  }

  void detach() noexcept final
  {
    // Work still queued reads the block: finish it before the block goes.
    mStream.synchronize();
    closeFrom(Scope::kRun);
    mArena.rebind({});
  }
  DeviceTrackerImpl(const DeviceTrackerImpl&) = delete;
  DeviceTrackerImpl& operator=(const DeviceTrackerImpl&) = delete;
  DeviceTrackerImpl(DeviceTrackerImpl&&) = delete;
  DeviceTrackerImpl& operator=(DeviceTrackerImpl&&) = delete;

  [[nodiscard]] StageStatus configureRun(std::span<const IterationParameters<Det>> iterations, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator) final
  {
    mStream.synchronize();
    closeFrom(Scope::kRun);
    if (!mParameters.allocate(mArena, iterations.size()) || !mGeometry.allocate(mArena) || !mMaterial.upload(mArena, material, mStream)) {
      closeFrom(Scope::kRun);
      return StageStatus::kOutOfMemory;
    }
    upload(iterations.data(), mParameters, iterations.size(), mStream);
    uploadObject(geometry, mGeometry, mStream);
    // The sensors by place, once, where an iteration follows.
    mSensorCells = {};
    if (std::ranges::any_of(iterations, [](const auto& parameters) { return parameters.follow || parameters.followEnd; })) {
      // Built on the host once per geometry: configureRun runs every timeframe.
      if (!mHostCells || std::memcmp(mHostCellsGeometry.data(), static_cast<const void*>(&geometry), sizeof(geometry)) != 0) {
        mHostCells.emplace(geometry);
        std::memcpy(mHostCellsGeometry.data(), static_cast<const void*>(&geometry), sizeof(geometry));
      }
      const SensorCells<Det>& cells = *mHostCells;
      if (!mCellLayers.allocate(mArena, cells.layers().size()) || !mCellOffsets.allocate(mArena, cells.offsets().size()) || !mCellSensors.allocate(mArena, std::max<std::size_t>(cells.sensors().size(), 1U))) {
        closeFrom(Scope::kRun);
        return StageStatus::kOutOfMemory;
      }
      upload(cells.layers().data(), mCellLayers, cells.layers().size(), mStream);
      upload(cells.offsets().data(), mCellOffsets, cells.offsets().size(), mStream);
      upload(cells.sensors().data(), mCellSensors, cells.sensors().size(), mStream);
      mSensorCells = {.layers = mCellLayers.data(), .offsets = mCellOffsets.data(), .sensors = mCellSensors.data()};
    }
    mStream.synchronize(); // the host objects may go once this returns
    mHostParameters.assign(iterations.begin(), iterations.end());
    mPropagator = propagator;
    mIterationCount = iterations.size();
    mConfigured = true;
    return StageStatus::kComplete;
  }

  [[nodiscard]] StageStatus beginTimeframe(const Frame& frame, const TimeframePlan<Det>& plan) final
  {
    if (!mConfigured) {
      throw std::logic_error("beginTimeframe before configureRun");
    }
    plan.check(frame, mHostParameters); // before anything changes
    mStream.synchronize();
    closeFrom(Scope::kTimeframe);
    mTimeframeScope.emplace(mArena.scope());
    const std::size_t rows = frame.size();
    const auto offsets = frame.grid().offsets();
    if (!mColumns.allocate(mArena, rows) || !mOffsets.allocate(mArena, offsets.size()) || !mUsed.allocate(mArena, rows) || !mFrame.allocate(mArena) || !uploadPlan(plan)) {
      closeFrom(Scope::kTimeframe);
      return StageStatus::kOutOfMemory;
    }
    upload<MeasurementT>(frame.measurements(), rows, mColumns, mStream);
    upload(offsets.data(), mOffsets, offsets.size(), mStream);
    mRuntime->fill(mUsed.data(), std::uint8_t{kFree}, rows, mStream);
    // The frame is a view: built on the host over the device columns, and
    // copied as the value it is.
    const Frame device{std::as_const(mColumns).view(), frame.grid().geometry(), std::span<const Index>{mOffsets.data(), offsets.size()}, static_cast<Index>(rows)};
    uploadObject(device, mFrame, mStream);
    mStream.synchronize();
    // The tracks' chunks, at the bottom of the scratch end for the timeframe.
    mTracksScope.emplace(mArena.scope(ArenaEnd::kScratch));
    mRows = rows;
    mHostFrame = &frame;
    mPlan = &plan;
    return StageStatus::kComplete;
  }

  void setFlatChunkPairs(std::size_t pairs) final { mFlatOptions.maxChunkPairs = pairs; }

  void setStagingCapacity(std::size_t records) final
  {
    mTrackletStage.prepareStaging(records);
    mNeighbourStage.prepareStaging(records);
    mRecordStage.prepareStaging(records);
  }

  void beginIteration(std::size_t iteration) final
  {
    requireTimeframe();
    if (mIterationOpen) {
      throw std::logic_error("beginIteration with an iteration open");
    }
    if (iteration >= mIterationCount) {
      throw std::invalid_argument("beginIteration of an iteration the run does not have");
    }
    // Every call before returns with the stream drained.
    // The iteration's own arena peak, in a window endIteration closes.
    mOuterArenaPeak = mArena.peak();
    mArena.resetPeak();
    mChunksScope.emplace(mArena.scope(ArenaEnd::kScratch));
    mIteration = iteration;
    mIterationStatus = StageStatus::kComplete;
    mIterationOpen = true;
  }

  [[nodiscard]] StageStatus runSlice(std::size_t slice, DeviceSliceDownload<Det>* download) final
  {
    if (!mIterationOpen) {
      throw std::logic_error("runSlice outside an open iteration");
    }
    if (mIterationStatus != StageStatus::kComplete) {
      throw std::logic_error("runSlice in an iteration a slice has failed");
    }
    if (slice >= mPlan->slices.size()) {
      throw std::logic_error("runSlice of a slice the plan does not have");
    }
    openSlice();
    // An iteration that takes the previous one's graph has its records from
    // that one's second fit (one graph for both rounds): no slice of its own,
    // as on the host.
    if (mHostParameters[mIteration].sharedGraph) {
      return StageStatus::kComplete;
    }
    const StageStatus status = runStages(mIteration, slice, download);
    mStream.synchronize(); // the host objects the slice reads are read until here
    if (status != StageStatus::kComplete) {
      mIterationStatus = status;
      closeFrom(Scope::kSlice);
    }
    return status;
  }

  [[nodiscard]] StageStatus endIteration(DeviceIterationDownload<Det>* download) final
  {
    if (!mIterationOpen) {
      throw std::logic_error("endIteration without an open iteration");
    }
    closeFrom(Scope::kSlice); // runSlice returns with the stream drained
    StageStatus status = mIterationStatus;
    const bool shared = mHostParameters[mIteration].sharedGraph;
    const bool sharesNext = mIteration + 1U < mIterationCount && mHostParameters[mIteration + 1U].sharedGraph;
    mIterationReport = {.records = mChunks.rows()};
    if (status == StageStatus::kComplete) {
      // The records of every slice, in slice order, into one array on the
      // base end: the selection's input. A
      // shared iteration's are its predecessor's second fit, kept for it in a
      // scope of their own below the iteration's (`mSharedScope`).
      if (sharesNext) {
        mSharedScope.emplace(mArena.scope());
        mSharedGathered = {};
        if (!mSharedGathered.allocate(mArena, mSharedChunks.rows())) {
          status = StageStatus::kOutOfMemory;
        } else {
          mSharedChunks.gather(mSharedGathered, 0U, mStream);
        }
      }
      mIterationScope.emplace(mArena.scope());
      mGathered = {};
      const std::size_t recordRows = shared ? mSharedGathered.size() : mChunks.rows();
      mIterationReport.records = recordRows;
      if (status == StageStatus::kComplete && !mGathered.allocate(mArena, recordRows)) {
        status = StageStatus::kOutOfMemory;
      }
      if (status == StageStatus::kComplete) {
        if (shared) {
          copyColumns(mSharedGathered, 0U, mGathered, 0U, recordRows, mStream);
        } else {
          mChunks.gather(mGathered, 0U, mStream);
        }
        if (download != nullptr) {
          gtrk::download(mGathered, download->records, recordRows, mStream);
        }
      }
      // The gather is queued before anything that reuses the chunks' space.
      closeFrom(Scope::kChunks);
      if (status == StageStatus::kComplete) {
        status = endStages(download);
      }
      mStream.synchronize(); // the downloads land, and nothing queued reads the iteration scope
    }
    closeFrom(Scope::kChunks);
    closeFrom(Scope::kIteration);
    // The shared records go with the iteration that took them, or with a
    // failure of the one that made them.
    if (shared || status != StageStatus::kComplete) {
      mSharedGathered = {};
      mSharedScope.reset();
    }
    mIterationReport.arenaPeakBytes = mArena.peak();
    mArena.raisePeak(mOuterArenaPeak);
    return status;
  }

  void downloadTracks(Tracks& tracks) final
  {
    requireNoIteration("downloadTracks");
    tracks.resize(mTrackChunks.rows());
    mTrackChunks.download(tracks, 0U, mStream);
    mStream.synchronize();
  }

  [[nodiscard]] StageStatus findVertices(std::size_t pass, const VertexerSettings& settings, const VertexConstraint& beamSpot, VertexerResult& result) final
  {
    requireNoIteration("findVertices");
    if (pass >= mIterationCount) {
      throw std::invalid_argument("findVertices of an iteration the run does not have");
    }
    // Slice by slice, each its own frames' vertices; the lists in time order.
    result = {};
    const auto& slices = mPlan->slices;
    StageStatus status = StageStatus::kComplete;
    for (std::size_t slice = 0; slice < slices.size() && status == StageStatus::kComplete; ++slice) {
      openSlice();
      DeviceStageOutput<TrackletT> tracklets;
      DeviceSpan<Index> permutation;
      DeviceStageOutput<TripletT> triplets;
      status = tripletStages(pass, pass, slice, tracklets, permutation, triplets);
      VertexerResult found;
      if (status == StageStatus::kComplete) {
        status = vertex_glue::findVertices<Det>(*mRuntime, mArena, mStream, std::as_const(mColumns).view(), static_cast<Index>(mRows), constView(tracklets), constView(triplets), triplets.size(), mParameters.data() + pass, settings, beamSpot, ownedFrames(mHostFrame->grid().geometry(), slices[slice].slice.coreBlocks), found);
      }
      // A complete vertexer has synchronised after its last readback.
      if (status != StageStatus::kComplete) {
        mStream.synchronize();
      }
      closeFrom(Scope::kSlice);
      appendVertices(result, found);
    }
    return status;
  }

  [[nodiscard]] StageStatus replacePlan(const TimeframePlan<Det>& plan) final
  {
    requireNoIteration("replacePlan");
    plan.check(*mHostFrame, mHostParameters); // before anything changes
    mStream.synchronize();
    closeFrom(Scope::kSlice);
    // The new device part goes on top of the timeframe's base end, the old
    // one's few kB staying until the timeframe ends.
    if (!uploadPlan(plan)) {
      return StageStatus::kOutOfMemory;
    }
    mPlan = &plan;
    return StageStatus::kComplete;
  }

  void endTimeframe() noexcept final
  {
    // Work still queued reads the timeframe's storage.
    mStream.synchronize();
    closeFrom(Scope::kTimeframe);
  }

  [[nodiscard]] const IterationReport& lastIteration() const final { return mIterationReport; }

  [[nodiscard]] const Report& lastSlice() const final { return mReport; }
  [[nodiscard]] std::size_t arenaBytes() const final { return mArena.capacity(); }
  [[nodiscard]] std::size_t arenaUsedBytes(ArenaEnd end) const final { return mArena.used(end); }
  [[nodiscard]] std::size_t arenaPeakBytes() const final { return mArena.peak(); }

 private:
  struct LayeredNeighbours {
    NeighbourStage<Det> body;
    const Index* layerBoundaries;
    std::array<LayerID, kLayerCount> layerOrder = ::o2::gtrk::layerDepthGroups<Det>.order;
  };

  /// The plan's device part, in the timeframe scope: every slice's source map and core ranges, once for every
  /// iteration -- the core ranges with an iteration's start-layer mask per
  /// iteration, where one has a mask --, each iteration's vertex windows, the
  /// beam line and then the vertices, the follow's bound.
  /// False, with the scope to be closed by the caller, if the arena refuses.
  [[nodiscard]] bool uploadPlan(const TimeframePlan<Det>& plan)
  {
    std::vector<LayerRowMap<Det>> sources;
    std::vector<MeasurementRange> cores;
    std::vector<MeasurementRange> maskedCores; // per iteration, with its start-layer mask (RoadCount.h, `roadCores`)
    const std::vector<VertexWindowParameters<Det>>& windows = plan.vertexWindows;
    for (const auto& slice : plan.slices) {
      sources.push_back(slice.sources);
      cores.insert(cores.end(), slice.slice.core.begin(), slice.slice.core.end());
    }
    const bool startMasks = std::ranges::any_of(mHostParameters, [](const IterationParameters<Det>& one) { return one.roadStartLayers != 0U; });
    if (startMasks) {
      for (std::size_t iteration = 0U; iteration < plan.parameters.size(); ++iteration) {
        for (const auto& slice : plan.slices) {
          // A graph shared with the next iteration has its start layers too (Tracker.h).
          const auto masked = roadCores(std::span<const MeasurementRange, kLayerCount>{slice.slice.core}, mHostParameters[iteration].roadStartLayers);
          maskedCores.insert(maskedCores.end(), masked.begin(), masked.end());
        }
      }
    }
    std::vector<VertexSeed> seeds{plan.beamLine};
    seeds.insert(seeds.end(), plan.vertices.begin(), plan.vertices.end());
    if (!mSliceSources.allocate(mArena, sources.size()) || !mSliceCores.allocate(mArena, cores.size()) || (startMasks && !mSliceRoadCores.allocate(mArena, maskedCores.size())) || !mVertexWindows.allocate(mArena, windows.size()) || !mSeeds.allocate(mArena, seeds.size()) || !mLayerMaxWidthPS.allocate(mArena, kLayerCount) || !mTargetBounds.allocate(mArena, kLayerCount)) {
      return false;
    }
    upload(sources.data(), mSliceSources, sources.size(), mStream);
    upload(cores.data(), mSliceCores, cores.size(), mStream);
    if (startMasks) {
      upload(maskedCores.data(), mSliceRoadCores, maskedCores.size(), mStream);
    }
    mRoadCoresMasked = startMasks;
    upload(windows.data(), mVertexWindows, windows.size(), mStream);
    upload(seeds.data(), mSeeds, seeds.size(), mStream);
    upload(plan.layerMaxWidthPS.data(), mLayerMaxWidthPS, kLayerCount, mStream);
    upload(plan.layerTargetBounds.data(), mTargetBounds, kLayerCount, mStream);
    mStream.synchronize(); // the host vectors go when this returns
    return true;
  }

  template <template <template <class> class> class S>
  [[nodiscard]] static typename DeviceColumns<S>::ConstView constView(const DeviceStageOutput<S>& output)
  {
    return output.view();
  }

  /// A stage's run over a flat candidate list
  /// into \a output, its report into \a run, its readbacks counted.
  template <template <template <class> class> class S, typename Adapter>
  StageStatus runFlat(const Index* listed, Index listedCount, const Adapter& adapter, DeviceStageOutput<S>& output, FlatRunReport& run)
  {
    const StageStatus status = flat::runFlat(*mRuntime, mArena, mStream, listed, listedCount, adapter, output, run, mFlatOptions);
    mReport.readbacks += run.readbacks;
    return status;
  }

  /// A stage's run, counted: its readback, and a second one after a missed
  /// prediction.
  template <template <template <class> class> class S, typename Generate>
  StageStatus runStage(GpuStage<S>& gpuStage, DeviceStageOutput<S>& output, const Generate& generate)
  {
    return counted(gpuStage, gpuStage.run(output, generate));
  }

  /// The same over the \a listedCount elements \a listed names only, when set
  /// (`GpuStage::run`).
  template <template <template <class> class> class S, typename Generate>
  StageStatus runStage(GpuStage<S>& gpuStage, DeviceStageOutput<S>& output, const Generate& generate, const Index* listed, Index listedCount)
  {
    return counted(gpuStage, gpuStage.run(output, generate, listed, listedCount));
  }

  /// \a status of \a gpuStage's last run, its readbacks counted.
  template <template <template <class> class> class S>
  StageStatus counted(const GpuStage<S>& gpuStage, StageStatus status)
  {
    const bool missed = gpuStage.stagingPredictionMissed();
    mReport.readbacks += missed ? 2U : 1U;
    mReport.stagingPredictionMissed = mReport.stagingPredictionMissed || missed;
    return status;
  }

  /// A slice's tracklets and triplets, the first stages every slice runs
  /// and all the vertexing pass runs, into \a tracklets and \a triplets in
  /// the slice scope, with the triplet stage's first-tracklet order in
  /// \a permutation. The slice is \a iteration's; the parameters are
  /// \a graph's -- the successor's where it takes this iteration's graph,
  /// this iteration's otherwise.
  StageStatus tripletStages(std::size_t iteration, std::size_t graph, std::size_t sliceIndex, DeviceStageOutput<TrackletT>& tracklets, DeviceSpan<Index>& permutation, DeviceStageOutput<TripletT>& triplets)
  {
    const IterationParameters<Det>* parameters = mParameters.data() + graph;
    const IterationParameters<Det>& hostParameters = mHostParameters[graph];
    const auto& slicePlan = mPlan->slices[sliceIndex];
    // The slices' sources and cores, uploaded once for every iteration.
    const LayerRowMap<Det>* sources = mSliceSources.data() + sliceIndex;
    const MeasurementRange* core = mSliceCores.data() + (sliceIndex * kLayerCount);
    const Index sourceCount = slicePlan.sources.size();

    // The tracklet stage as the host tracker builds it from the plan, bound to
    // the device copies.
    const auto inputs = slicePlan.inputs(hostParameters, mPlan->layerMaxWidthPS);
    const bool beamLine = hostParameters.vertexSource == VertexSource::kBeamLine;
    if (!tracklets.prepare(mArena, sourceCount)) {
      return StageStatus::kOutOfMemory;
    }
    // The beam-line stage for a beam-line iteration, its displaced one for a
    // displaced iteration;
    // each testing the shapes where the iteration does.
    StageStatus status = withTrackletVertices(hostParameters.vertexSource, hostParameters.impactMax > 0.F, [&](auto beamLineOnly, auto displaced) {
      return withShapeTest(hostParameters.shapeTest, [&](auto shape) {
        const TrackletStageFor<Det, decltype(beamLineOnly)::value, decltype(displaced)::value, decltype(shape)::value> stage{hostParameters, slicePlan.sources, *mHostFrame, inputs.targetBlocks, mPlan->beamLine, mPlan->vertices, inputs.maxCandidateWidthPS, mPlan->vertexWindows[iteration], nullptr, inputs.gate};
        const auto bound = stage.bind({
          .parameters = parameters,
          .sources = sources,
          .frame = mFrame.data(),
          .used = mUsed.data(),
          .vertexWindows = mVertexWindows.data() + iteration,
          .vertices = beamLine ? mSeeds.data() : mSeeds.data() + 1,
          .gateCore = core,
        });
        // A thread for the sources with a query only: a dead source's thread holds its warp's lanes idle while the
        // others walk; it emits nothing, so its count is the 0 the stage sets.
        auto scratch = mArena.scope(ArenaEnd::kScratch);
        DeviceSpan<Index> listed;
        Index listedCount = 0U;
        const StageStatus selected = algorithms::selectIndices(*mRuntime, mArena, mStream, sourceCount, tracker_detail::HasQuery<decltype(bound)>{bound}, listed, listedCount, ArenaEnd::kScratch);
        ++mReport.readbacks; // the list's length
        if (selected != StageStatus::kComplete) {
          return selected;
        }
        // The flat run: the sources' candidates listed, tested one a thread.
        std::uint32_t queriesPerSource = 1U;
        for (std::size_t layer = 0; layer < kLayerCount; ++layer) {
          queriesPerSource = std::max(queriesPerSource, hostParameters.outLinkOffsets[layer + 1U] - hostParameters.outLinkOffsets[layer]);
        }
        if (queriesPerSource > (1U << tracker_detail::FlatTracklets<decltype(bound)>::kQueryBits)) {
          throw std::logic_error("a layer has more out-links than a flat pair's key numbers");
        }
        mReport.trackletSources = listedCount;
        return runFlat(listedCount > 0U ? listed.data() : nullptr, listedCount, tracker_detail::FlatTracklets<decltype(bound)>{.stage = bound, .queriesPerSource = queriesPerSource}, tracklets, mReport.trackletRun);
      });
    });
    if (status != StageStatus::kComplete) {
      return status;
    }
    // Each source's tracklets in (link, lambda, target) order, which the
    // triplet stage searches.
    status = glue::sortTracklets(*mRuntime, mArena, mStream, tracklets);
    if (status != StageStatus::kComplete) {
      return status;
    }
    mReport.tracklets = tracklets.size();

    // The longest tracklet paths through each measurement, which the
    // permutation and the triplet stage read.
    glue::DevicePathBounds trackletPaths;
    status = glue::trackletPaths<Det>(*mRuntime, mArena, mStream, parameters, slicePlan.sources, constView(tracklets), tracklets.offsets().data(), trackletPaths);
    if (status != StageStatus::kComplete) {
      return status;
    }
    const auto paths = trackletPaths.view(hostParameters.minClusters);

    // The triplet stage's inputs: the first-tracklet order; the targets'
    // variance bounds are the plan's, over the whole frame, uploaded with it.
    Index firstTracklets = 0U;
    status = glue::permutation<Det>(*mRuntime, mArena, mStream, constView(tracklets), tracklets.size(), parameters, sources, permutation, firstTracklets, tracklets.offsets().data(), paths);
    ++mReport.readbacks;
    if (status != StageStatus::kComplete) {
      return status;
    }
    mReport.permutation = firstTracklets;

    if (!triplets.prepare(mArena, firstTracklets)) {
      return StageStatus::kOutOfMemory;
    }
    // The flat run: the windows'
    // candidates listed per first tracklet, tested one a thread; the
    // shape-testing stage where the iteration tests the shapes.
    status = withShapeTest(hostParameters.shapeTest, [&](auto shape) {
      const TripletStage<Det, decltype(shape)::value> tripletStage{parameters, sources, mFrame.data(), constView(tracklets), tracklets.offsets().data(), mTargetBounds.data(), permutation.data(), paths};
      return runFlat(nullptr, firstTracklets, tracker_detail::FlatTriplets<Det, decltype(shape)::value>{.stage = tripletStage}, triplets, mReport.tripletRun);
    });
    if (status != StageStatus::kComplete) {
      return status;
    }
    mReport.triplets = triplets.size();
    return StageStatus::kComplete;
  }

  StageStatus runStages(std::size_t iteration, std::size_t sliceIndex, DeviceSliceDownload<Det>* download)
  {
    // The graph of a slice whose successor takes it is built with the
    // successor's parameters (its kink allowance), as on the host; only the
    // nominal round's fit below takes this iteration's.
    const bool sharesNext = iteration + 1U < mIterationCount && mHostParameters[iteration + 1U].sharedGraph;
    const std::size_t graph = sharesNext ? iteration + 1U : iteration;
    const IterationParameters<Det>* parameters = mParameters.data() + graph;
    const IterationParameters<Det>& hostParameters = mHostParameters[graph];
    // The unmasked cores by slice; the masked ones by iteration and slice (`uploadPlan`).
    const MeasurementRange* core = mRoadCoresMasked ? mSliceRoadCores.data() + (((iteration * mPlan->slices.size()) + sliceIndex) * kLayerCount) : mSliceCores.data() + (sliceIndex * kLayerCount);
    DeviceStageOutput<TrackletT> tracklets;
    DeviceSpan<Index> permutation;
    DeviceStageOutput<TripletT> triplets;
    StageStatus status = tripletStages(iteration, graph, sliceIndex, tracklets, permutation, triplets);
    if (status != StageStatus::kComplete) {
      return status;
    }

    // The longest triplet paths through each tracklet, from the triplets
    // alone, and the triplets kept only where they are on one:
    // what the neighbour stage reads.
    glue::DevicePathBounds tripletPaths;
    status = glue::tripletPaths<Det>(*mRuntime, mArena, mStream, parameters, mPlan->slices[sliceIndex].sources, constView(tracklets), tracklets.size(), tracklets.offsets().data(), constView(triplets), triplets.offsets().data(), permutation.data(), mReport.permutation, tripletPaths);
    if (status != StageStatus::kComplete) {
      return status;
    }
    const auto paths = tripletPaths.view(hostParameters.minClusters);
    {
      DeviceStageOutput<TripletT> kept;
      status = glue::onPathTriplets(*mRuntime, mArena, mStream, paths, triplets, kept, tripletPaths);
      ++mReport.readbacks;
      if (status != StageStatus::kComplete) {
        return status;
      }
      triplets = kept;
    }
    mReport.triplets = triplets.size();

    DeviceSpan<Index> endOffsets;
    DeviceSpan<Index> endTriplets;
    status = glue::endIndex(*mRuntime, mArena, mStream, constView(triplets), triplets.size(), tracklets.size(), endOffsets, endTriplets);
    if (status != StageStatus::kComplete) {
      return status;
    }

    DeviceStageOutput<NeighbourT> neighbours;
    if (!neighbours.prepare(mArena, triplets.size())) {
      return StageStatus::kOutOfMemory;
    }
    road_glue::OrderedCounts counts;
    status = counts.prepare(*mRuntime, mArena, mStream, triplets.size());
    if (status != StageStatus::kComplete) {
      return status;
    }
    const NeighbourStage<Det> neighbourStage{parameters, mFrame.data(), constView(tracklets), constView(triplets), endOffsets.data(), endTriplets.data(), paths, triplets.view().level, counts.view(road_glue::linkSourcesOf(parameters), core)};
    status = runStage(mNeighbourStage, neighbours, LayeredNeighbours{neighbourStage, tripletPaths.layerTriplets.data()});
    if (status != StageStatus::kComplete) {
      return status;
    }
    mReport.neighbours = neighbours.size();

    // The chained count and the overflow flag in one readback; with a
    // download, also the chained triplets' list, which the tests compare
    // with the host's.
    Index chainedCount = 0U;
    status = counts.summarize(*mRuntime, mStream, constView(triplets), triplets.size(), chainedCount);
    mReport.readbacks += triplets.size() != 0U ? 1U : 0U;
    DeviceSpan<Index> chained;
    Index chainedListed = 0U;
    if (status == StageStatus::kComplete && download != nullptr) {
      status = glue::chained(*mRuntime, mArena, mStream, neighbours.offsets().data(), triplets.size(), chained, chainedListed);
      ++mReport.readbacks;
    }
    mReport.chained = chainedCount;
    if (status != StageStatus::kComplete) {
      return status;
    }

    // Counts and successor flags were finalized with the neighbours; overflow
    // was checked together with the chained summary above.
    const auto& successorAbove = counts.successorAbove;
    const auto& roads = counts.roads;
    DeviceSpan<Index> starts;
    Index startCount = 0U;
    status = road_glue::roadStarts(*mRuntime, mArena, mStream, constView(tracklets), constView(triplets), triplets.size(), successorAbove.data(), roads.data(), road_glue::roadStartRuleOf(parameters, hostParameters), starts, startCount);
    ++mReport.readbacks;
    if (status != StageStatus::kComplete) {
      return status;
    }
    DeviceSpan<Index> roadOffsets;
    Index roadCount = 0U;
    status = road_glue::roadOffsets(*mRuntime, mArena, mStream, starts.data(), startCount, roads.data(), roadOffsets, roadCount);
    ++mReport.readbacks;
    if (status != StageStatus::kComplete) {
      return status;
    }
    mReport.roads = roadCount;

    // The road fit, and its records kept for the iteration.
    Records records;
    if (!records.prepare(mArena, roadCount)) {
      return StageStatus::kOutOfMemory;
    }
    const RoadWalk walk{starts.data(), roadOffsets.data(), startCount, neighbours.offsets().data(), neighbours.view().from, triplets.view().level, roads.data()};
    // An iteration with a kink allowance takes it on the fit's last step
    // (against the road-end fakes): its own instantiation.
    // An iteration with the beam constraint fits on the stage that also
    // tests the beam line, as the host.
    // With a successor that takes this graph (one graph for both rounds; the
    // host: Tracker::runRoads), the nominal round's
    // roads are fitted with this iteration's settings into its records and
    // every road with the successor's into the successor's (`mSharedChunks`),
    // which its endIteration takes as its records.
    const auto fitRoads = [&]<bool NominalOnly>(std::size_t owner, Records& into) {
      const IterationParameters<Det>* ownerParameters = mParameters.data() + owner;
      const IterationParameters<Det>& ownerHost = mHostParameters[owner];
      const auto run = [&]<bool BeamLine, bool LastStepAllowance, bool InnerTrim>() {
        const RoadFitStage<Det, Material, BeamLine, NominalOnly, LastStepAllowance, InnerTrim> fit{ownerParameters, mFrame.data(), mGeometry.data(), mMaterial.data(), mPropagator, constView(tracklets), constView(triplets), walk};
        return runStage(mRecordStage, into, fit);
      };
      // An iteration with innerTrim fits a refused road again without its
      // innermost clusters, as the host.
      return withBools([&]<bool BeamLine, bool LastStepAllowance, bool InnerTrim>(std::bool_constant<BeamLine>, std::bool_constant<LastStepAllowance>, std::bool_constant<InnerTrim>) { return run.template operator()<BeamLine, LastStepAllowance, InnerTrim>(); }, ownerHost.beamConstraint, takesLastStepAllowance(ownerHost), ownerHost.innerTrim != 0U);
    };
    status = sharesNext ? fitRoads.template operator()<true>(iteration, records) : fitRoads.template operator()<false>(iteration, records);
    if (status != StageStatus::kComplete) {
      return status;
    }
    mReport.records = records.size();
    if (!mChunks.append(mArena, ArenaEnd::kScratch, records.output(), 0U, records.size(), mStream)) {
      return StageStatus::kOutOfMemory;
    }
    if (sharesNext) {
      Records shared;
      if (!shared.prepare(mArena, roadCount)) {
        return StageStatus::kOutOfMemory;
      }
      status = fitRoads.template operator()<false>(iteration + 1U, shared);
      if (status != StageStatus::kComplete) {
        return status;
      }
      if (!mSharedChunks.append(mArena, ArenaEnd::kScratch, shared.output(), 0U, shared.size(), mStream)) {
        return StageStatus::kOutOfMemory;
      }
    }

    if (download != nullptr) {
      gtrk::download(successorAbove, download->successorAbove, triplets.size(), mStream);
      gtrk::download(roads, download->roads, triplets.size(), mStream);
      gtrk::download(starts, download->starts, startCount, mStream);
      gtrk::download(roadOffsets, download->roadOffsets, std::size_t{startCount} + 1U, mStream);
      gtrk::download(records, download->records, mStream);
      gtrk::download(tracklets, download->tracklets, mStream);
      gtrk::download(triplets, download->triplets, mStream);
      gtrk::download(neighbours, download->neighbours, mStream);
      gtrk::download(chained, download->chained, chainedListed, mStream);
    }
    return StageStatus::kComplete; // runSlice synchronises
  }

  /// The end of an iteration after its gather, in the iteration scope: the
  /// selection, the refit and the
  /// follow on a working copy of the used mask, the tracks appended to
  /// the timeframe's, and only then the working mask committed. A failure
  /// before the commit leaves the used mask and the tracks as they were.
  [[nodiscard]] StageStatus endStages(DeviceIterationDownload<Det>* download)
  {
    const Index rows = static_cast<Index>(mRows);
    const auto records = std::as_const(mGathered).view();
    const auto recordCount = static_cast<Index>(mGathered.size());
    const IterationParameters<Det>& hostParameters = mHostParameters[mIteration];
    const IterationParameters<Det>* parameters = mParameters.data() + mIteration;
    auto& report = mIterationReport;

    DeviceSpan<std::uint8_t> working;
    if (!working.allocate(mArena, rows)) {
      return StageStatus::kOutOfMemory;
    }
    mRuntime->copyOnDevice(working.data(), std::as_const(mUsed).data(), rows, mStream);

    selection_glue::Selected selected;
    StageStatus status = selection_glue::select<kMaxClusters>(*mRuntime, mArena, mStream, records, recordCount, working.data(), rows, hostParameters.sharedInnermost, selected);
    report.readbacks = selected.readbacks;
    if (status != StageStatus::kComplete) {
      return status;
    }
    report.duplicates = selected.duplicates;
    report.mismatches = selected.mismatches;
    report.rounds = selected.rounds;
    report.tracks = selected.tracks;

    // The refit, one element per selected record, into the iteration's tracks.
    const Index trackCount = selected.tracks;
    DeviceColumns<TrackRecord::template RecordT> tracks;
    DeviceSpan<Proposal> proposals;
    DeviceSpan<follow::KinkSteps> kinkSteps;
    DeviceSpan<ClusterList<kMaxClusters>> followRows;
    DeviceSpan<std::uint8_t> followCounts;
    DeviceSpan<BoundTrackState> endStarts;
    const Index followCount = hostParameters.follow || hostParameters.followEnd ? trackCount : 0U;
    const Index kinkStepCount = Follow<Det>::takesKinkSteps(hostParameters) ? followCount : 0U;
    if (!tracks.allocate(mArena, trackCount) || !proposals.allocate(mArena, followCount) || !kinkSteps.allocate(mArena, kinkStepCount) || !followRows.allocate(mArena, followCount) || !followCounts.allocate(mArena, followCount) || !endStarts.allocate(mArena, hostParameters.followEnd ? trackCount : 0U)) {
      return StageStatus::kOutOfMemory;
    }
    const TrackRefitElement<Det, Material> refit{mFrame.data(), mGeometry.data(), mMaterial.data(), mPropagator, records.clusters, records.count, records.outwardState, records.chi2, std::as_const(selected.selected).data(), tracks.view(), 0U, static_cast<std::uint8_t>(mIteration), kinkBudget(hostParameters)};
    algorithms::forEach(*mRuntime, mStream, trackCount, refit);

    // The follow, the host's sequence
    // (`follow::runSequence`), its bids in the rounds' slots. The host's
    // priority, trackCount - (trackOffset + track), is the iteration's count
    // less the track.
    if ((hostParameters.follow || hostParameters.followEnd) && trackCount != 0U) {
      runFollow(mIteration, tracks, trackCount, false, proposals.data(), kinkStepCount != 0U ? kinkSteps.data() : nullptr, followRows, followCounts, selected.slots, rows, working, hostParameters.followEnd ? endStarts.data() : nullptr);
    }

    // A length or a pT floor after the follow (the host:
    // Tracker::dropShortTracks): the
    // iteration's tracks below them dropped, the kept ones and their proposals compacted, and the mask
    // rebuilt as the one before the selection plus the kept tracks' rows as
    // the follow left them.
    Index keptCount = trackCount;
    if (dropsAfterFollow(hostParameters) && trackCount != 0U) {
      auto scratch = mArena.scope(ArenaEnd::kScratch);
      DeviceSpan<Index> kept;
      status = algorithms::selectIndices(*mRuntime, mArena, mStream, trackCount, tracker_detail::KeptAfterFollow<Det>{parameters, std::as_const(tracks).view().count, std::as_const(tracks).view().state}, kept, keptCount, ArenaEnd::kScratch);
      ++report.readbacks;
      if (status != StageStatus::kComplete) {
        return status;
      }
      if (keptCount != trackCount) {
        DeviceColumns<TrackRecord::template RecordT> keptTracks;
        DeviceSpan<Proposal> keptProposals;
        DeviceSpan<std::uint8_t> cut;
        if (!keptTracks.allocate(mArena, keptCount) || !keptProposals.allocate(mArena, followCount != 0U ? keptCount : 0U) || !cut.allocate(mArena, rows)) {
          return StageStatus::kOutOfMemory;
        }
        device_storage_detail::forEachColumnPair(std::as_const(tracks).view(), keptTracks.view(), [&](auto* to, const auto* from) { algorithms::forEach(*mRuntime, mStream, keptCount, algorithms::Gather<std::remove_pointer_t<decltype(to)>>{from, kept.data(), to}); });
        if (followCount != 0U) {
          algorithms::forEach(*mRuntime, mStream, keptCount, algorithms::Gather<Proposal>{std::as_const(proposals).data(), kept.data(), keptProposals.data()});
        }
        mRuntime->copyOnDevice(cut.data(), std::as_const(mUsed).data(), rows, mStream);
        algorithms::forEach(*mRuntime, mStream, keptCount, tracker_detail::RestoreRows<kMaxClusters>{std::as_const(keptTracks).view().clusters, std::as_const(keptTracks).view().count, std::as_const(working).data(), cut.data()});
        tracks = keptTracks;
        proposals = keptProposals;
        working = cut;
      }
    }

    // The successor's recovery pass (the host: Tracker::runIteration's
    // followEarlier): every track so far followed
    // again with the successor's settings before it makes its seeds, the rows
    // it takes its own. Here, with this iteration's tracks final and before
    // the commit, nothing runs in between: the tracks so far gathered with
    // this iteration's, followed, and put back as one chunk; the old chunks'
    // storage stays until the timeframe ends.
    bool committed = false;
    if (mIteration + 1U < mIterationCount && mHostParameters[mIteration + 1U].followEarlier && mHostParameters[mIteration + 1U].follow && mTrackChunks.rows() + keptCount != 0U) {
      status = followEarlier(mIteration + 1U, tracks, keptCount, working, selected.slots, rows);
      if (status != StageStatus::kComplete) {
        return status;
      }
      committed = true;
    }

    if (download != nullptr) {
      gtrk::download(selected.ranked, download->ranked, selected.unique, mStream);
      gtrk::download(selected.selected, download->selected, trackCount, mStream);
      gtrk::download(proposals, download->proposals, proposals.size(), mStream);
      gtrk::download(tracks, download->tracks, keptCount, mStream);
      gtrk::download(working, download->used, rows, mStream);
    }

    // The last fallible step, then the commit.
    if (!committed && !mTrackChunks.append(mArena, ArenaEnd::kScratch, tracks, 0U, keptCount, mStream)) {
      return StageStatus::kOutOfMemory;
    }
    mRuntime->copyOnDevice(mUsed.data(), std::as_const(working).data(), rows, mStream);
    return StageStatus::kComplete; // endIteration synchronises
  }

  /// The recovery pass of iteration \a successor over every track so far --
  /// the chunks' and this iteration's \a tracks (\a trackCount of them) --
  /// with \a working as the used mask, its bids in \a slots (one per row,
  /// \a rows): the host's `Follow::run` with `earlier`. The tracks, updated,
  /// replace the chunks as one.
  [[nodiscard]] StageStatus followEarlier(std::size_t successor, const DeviceColumns<TrackRecord::template RecordT>& tracks, Index trackCount, DeviceSpan<std::uint8_t>& working, DeviceSpan<Index>& slots, Index rows)
  {
    const IterationParameters<Det>& host = mHostParameters[successor];
    const auto earlier = static_cast<Index>(mTrackChunks.rows());
    const Index total = earlier + trackCount;
    DeviceColumns<TrackRecord::template RecordT> all;
    DeviceSpan<Proposal> proposals;
    DeviceSpan<follow::KinkSteps> kinkSteps;
    DeviceSpan<ClusterList<kMaxClusters>> followRows;
    DeviceSpan<std::uint8_t> followCounts;
    const Index kinkStepCount = Follow<Det>::takesKinkSteps(host) ? total : 0U;
    if (!all.allocate(mArena, total) || !proposals.allocate(mArena, total) || !kinkSteps.allocate(mArena, kinkStepCount) || !followRows.allocate(mArena, total) || !followCounts.allocate(mArena, total)) {
      return StageStatus::kOutOfMemory;
    }
    mTrackChunks.gather(all, 0U, mStream);
    copyColumns(tracks, 0U, all, earlier, trackCount, mStream);
    runFollow(successor, all, total, true, proposals.data(), kinkStepCount != 0U ? kinkSteps.data() : nullptr, followRows, followCounts, slots, rows, working, nullptr);
    mTrackChunks.clear();
    if (!mTrackChunks.append(mArena, ArenaEnd::kScratch, all, 0U, total, mStream)) {
      return StageStatus::kOutOfMemory;
    }
    return StageStatus::kComplete;
  }

  /// The follow of iteration \a iteration over \a count tracks of \a tracks
  /// (`follow::runSequence`, as the host's `Follow::run`): \a earlier the
  /// earlier iterations' tracks (`followEarlier`), the proposals into
  /// \a settled, the kink steps into \a steps where the iteration keeps
  /// them, the combined rows into \a followRows and \a followCounts, the
  /// bids in \a slots (one per row, \a rows of them), the marks in
  /// \a working, the end starts in \a endStarts where the iteration decides
  /// the end again.
  void runFollow(std::size_t iteration, DeviceColumns<TrackRecord::template RecordT>& tracks, Index count, bool earlier, Proposal* settled, follow::KinkSteps* steps, DeviceSpan<ClusterList<kMaxClusters>>& followRows, DeviceSpan<std::uint8_t>& followCounts, DeviceSpan<Index>& slots, Index rows, DeviceSpan<std::uint8_t>& working, BoundTrackState* endStarts)
  {
    const IterationParameters<Det>& host = mHostParameters[iteration];
    follow::runSequence(tracker_detail::FollowLaunch{.runtime = *mRuntime, .stream = mStream}, follow::Inputs<Det, Material>{
                                                                                                 .parameters = mParameters.data() + iteration,
                                                                                                 .choices = &host,
                                                                                                 .frame = mFrame.data(),
                                                                                                 .geometry = mGeometry.data(),
                                                                                                 .material = mMaterial.data(),
                                                                                                 .propagator = mPropagator,
                                                                                                 .tracks = tracks.view(),
                                                                                                 .constTracks = std::as_const(tracks).view(),
                                                                                                 .trackOffset = 0U,
                                                                                                 .trackCount = count,
                                                                                                 .inward = host.sharedInnermost == 0U,
                                                                                                 .earlier = earlier,
                                                                                                 .used = working.data(),
                                                                                                 .rowCount = rows,
                                                                                                 .maxCandidateWidthPS = std::as_const(mLayerMaxWidthPS).data(),
                                                                                                 .cells = mSensorCells,
                                                                                                 .proposals = settled,
                                                                                                 .bids = slots.data(),
                                                                                                 .rows = followRows.data(),
                                                                                                 .counts = followCounts.data(),
                                                                                                 .steps = steps,
                                                                                                 .endStarts = endStarts,
                                                                                               });
  }

  /// `requireTimeframe`, and no iteration open: \a call names the caller.
  void requireNoIteration(const char* call) const
  {
    requireTimeframe();
    if (mIterationOpen) {
      throw std::logic_error(std::string{call} + " with an iteration open");
    }
  }

  /// A slice's scope on the arena, the previous one closed and the report
  /// cleared.
  void openSlice()
  {
    // Every call before returns with the stream drained.
    closeFrom(Scope::kSlice);
    mSliceScope.emplace(mArena.scope());
    mReport = {};
  }

  void requireTimeframe() const
  {
    if (!mTimeframeScope) {
      throw std::logic_error("no timeframe has been uploaded to the device");
    }
  }

  /// The tracker's scopes: the run is the
  /// arena itself; the timeframe holds, on the base end, the slice or the
  /// iteration's -- never both at once -- and, on the scratch end, the
  /// tracks and above them the open iteration's record chunks.
  enum class Scope : std::uint8_t { kRun,
                                    kTimeframe,
                                    kChunks,
                                    kIteration,
                                    kSlice };

  /// Close \a scope and every scope inside it, innermost first: the only
  /// place scopes are closed. The slice goes with each; the chunks and the
  /// iteration's scope go with their own call or with the timeframe, the
  /// tracks with the timeframe.
  void closeFrom(Scope scope) noexcept
  {
    mSliceScope.reset();
    switch (scope) {
      case Scope::kSlice:
        return;
      case Scope::kChunks:
        closeChunks();
        return;
      case Scope::kIteration:
        mIterationScope.reset();
        mGathered = {};
        return;
      case Scope::kTimeframe:
      case Scope::kRun:
        break;
    }
    mIterationScope.reset();
    mGathered = {};
    mSharedGathered = {};
    mSharedScope.reset();
    closeChunks();
    mTrackChunks.clear();
    mTracksScope.reset();
    mTimeframeScope.reset();
    mRows = 0U;
    mHostFrame = nullptr;
    mPlan = nullptr;
    if (scope == Scope::kTimeframe) {
      return;
    }
    mArena.reset();
    mConfigured = false;
    mIterationCount = 0U;
    mParameters = {};
    mHostParameters.clear();
    mGeometry = {};
    mMaterial.release();
    mCellLayers = {};
    mCellOffsets = {};
    mCellSensors = {};
    mSensorCells = {};
  }

  void closeChunks() noexcept
  {
    mChunks.clear();
    mSharedChunks.clear();
    mChunksScope.reset();
    mIterationOpen = false;
  }

  DeviceRuntime* mRuntime = nullptr;
  DeviceArena mArena;
  DeviceStream mStream;
  GpuStage<TrackletT> mTrackletStage;
  GpuStage<NeighbourT> mNeighbourStage;
  GpuStage<SelectionRecordOf<kMaxClusters>::template RecordT> mRecordStage;
  DeviceColumnChunks<SelectionRecordOf<kMaxClusters>::template RecordT> mChunks;
  DeviceColumnChunks<SelectionRecordOf<kMaxClusters>::template RecordT> mSharedChunks; ///< the successor's records, where it takes this graph
  DeviceColumns<SelectionRecordOf<kMaxClusters>::template RecordT> mGathered;
  DeviceColumns<SelectionRecordOf<kMaxClusters>::template RecordT> mSharedGathered; ///< `mSharedChunks` gathered, in `mSharedScope`
  DeviceColumnChunks<TrackRecord::template RecordT> mTrackChunks;
  IterationReport mIterationReport;
  std::size_t mOuterArenaPeak = 0U; ///< the arena's peak before the open iteration's window
  flat::FlatOptions mFlatOptions;   ///< the flat runs' test hook (`setFlatChunkPairs`)
  Report mReport;

  DeviceSpan<IterationParameters<Det>> mParameters;
  std::vector<IterationParameters<Det>> mHostParameters; ///< the same, for the values the host reads
  DeviceObject<TrackingGeometry<Det>> mGeometry;
  DeviceMaterial<Material> mMaterial;
  DeviceSpan<SensorCellLayer> mCellLayers; ///< the sensors by place, for the outward follow
  DeviceSpan<std::uint32_t> mCellOffsets;
  DeviceSpan<std::uint32_t> mCellSensors;
  SensorCellsView mSensorCells;
  /// The sensor cells on the host and the geometry they were built from.
  std::optional<SensorCells<Det>> mHostCells;
  std::array<std::byte, sizeof(TrackingGeometry<Det>)> mHostCellsGeometry{};
  Propagator mPropagator;
  std::size_t mIterationCount = 0U;
  bool mConfigured = false;

  std::optional<ArenaScope> mTimeframeScope;
  DeviceColumns<MeasurementT> mColumns;
  DeviceSpan<Index> mOffsets;
  DeviceSpan<std::uint8_t> mUsed;
  DeviceObject<Frame> mFrame;
  std::size_t mRows = 0U;
  const Frame* mHostFrame = nullptr;
  const TimeframePlan<Det>* mPlan = nullptr;
  DeviceSpan<LayerRowMap<Det>> mSliceSources;
  DeviceSpan<MeasurementRange> mSliceCores;
  DeviceSpan<MeasurementRange> mSliceRoadCores; ///< mSliceCores with each iteration's start-layer mask, where one has a mask
  bool mRoadCoresMasked = false;
  DeviceSpan<VertexWindowParameters<Det>> mVertexWindows;
  DeviceSpan<VertexSeed> mSeeds; ///< the beam line, then the vertices
  DeviceSpan<std::uint32_t> mLayerMaxWidthPS;
  DeviceSpan<TripletTargetBounds> mTargetBounds; ///< the plan's, over the whole frame
  std::optional<ArenaScope> mTracksScope;

  std::optional<ArenaScope> mChunksScope;
  std::optional<ArenaScope> mIterationScope;
  std::optional<ArenaScope> mSharedScope; ///< the shared records', below the iteration's that takes them
  std::optional<ArenaScope> mSliceScope;
  std::size_t mIteration = 0U;
  StageStatus mIterationStatus = StageStatus::kComplete;
  bool mIterationOpen = false;
};

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE

#endif // GTRK_TRACKING_GPU_DEVICETRACKER_CUH_
