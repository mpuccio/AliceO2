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

/// \file TrackingInterface.h
/// \brief One detector's tracking on a timeframe of measurements, for the
///        device that holds the run's conditions.
///
/// It owns what the tracking is made of: the geometry, the material's budgets
/// and scattering table, the tracking geometry, the tracker and the own
/// vertexer. It borrows the run's conditions (`TrackingConditions`): the
/// field, the measurement depths and the material LUT stay the enclosing
/// device's, which requests them; the interface makes no request of its own.
///
/// Per timeframe, in this order:
/// - `configureFor`: `configure` on the first timeframe and for a new shipped
///   beam line;
/// - `checkInput`: the columns, the grid offsets, the variances, the seeds;
/// - `vertices`: the input seeds, or the own vertexer's;
/// - `track`: plan and run;
/// - `endTimeframe`, through the guard `endTimeframeOnExit` taken at the
///   top: on success, on a failure, on an exception. It ends what began, and
///   the results are invalid after.
///
/// A timeframe that does not complete: out of memory (`memoryLimitGB`) is
/// dropped with `dropTFUponFailure`, the result then reported as dropped;
/// without it, and for every other status, the run fails.

#ifndef GTRK_WORKFLOW_TRACKINGINTERFACE_H_
#define GTRK_WORKFLOW_TRACKINGINTERFACE_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/GPU/DeviceContext.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceTracker.h"
#include "GenericTracking/Tracking/LayerMaterial.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/Tracker.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/VertexingPass.h"
#include "GenericTracking/Workflow/BeamLine.h"
#include "GenericTracking/Workflow/GeometryProvider.h"
#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Workflow/MaterialLUT.h"
#include "GenericTracking/Workflow/RecoParam.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <fairlogger/Logger.h>
#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk::workflow
{

/// Whether every one of \a columns has \a rows rows.
template <class Columns>
[[nodiscard]] bool measurementColumnsHaveSize(const Columns& columns, std::size_t rows)
{
  return columns.apply([rows](const auto&... column) { return ((column.size() == rows) && ...); });
}

/// Whether a measurement's variances can enter the fit: finite and positive,
/// the spatial block positive definite. A zero variance divides the update's
/// whitening and the U-D factorisation by zero, and nothing downstream
/// reports the NaN that makes.
[[nodiscard]] inline bool measurementVariancesValid(Scalar varLoc0, Scalar covLoc01, Scalar varLoc1) { return isFiniteScalar(varLoc0) && isFiniteScalar(covLoc01) && isFiniteScalar(varLoc1) && varLoc0 > 0.F && varLoc1 > 0.F && covLoc01 * covLoc01 < varLoc0 * varLoc1; }

/// What the tracking borrows from the device that holds the run's conditions.
struct TrackingConditions {
  double bzKiloGauss = 0.;                            ///< the run's nominal field
  std::vector<double> depthsCM;                       ///< each layer's measurement depth; empty for 0 on every layer
  const o2::base::MatLayerCylSet* material = nullptr; ///< the run's material LUT (MaterialLUT.h)
};

/// The view of \a conditions (Conditions.h's, any detector's, with the
/// material LUT: `Conditions<Det, true>`), once in for the timeframe.
template <class Conditions>
  requires requires(const Conditions& conditions) { conditions.materialLut(); }
[[nodiscard]] TrackingConditions trackingConditions(const Conditions& conditions)
{
  return {.bzKiloGauss = nominalBzKiloGauss(), .depthsCM = conditions.measurementDepthsCM(), .material = conditions.materialLut()};
}

/// Where the tracking runs: the host, or a GPU backend (the tracking
/// device's options `--gtrk-gpu-backend`, `--gtrk-gpu-device`, `--gtrk-gpu-arena-mb`).
struct TrackingDevice {
  std::string backend;     ///< empty: the host; else "cuda" or "hip"
  int device = 0;          ///< the backend's device index
  std::size_t arenaMB = 0; ///< the device tracker's arena in O2's GPU pool [MB]
};

template <DetectorLayout Det>
class TrackingInterface
{
 public:
  using Material = o2::base::MatLayerCylSet;
  using Controller = Tracker<Det, Material>;
  using Vertexing = VertexingPass<Det, Material>;
  using Frame = TimeFrame<Det>;
  using TimeframeResult = typename Controller::TimeframeResult;
  using Device = DeviceTracker<Det, Material>;

  /// \a withVertices: the native seeds (DataSpecs.h's `kVertices`) are an input.
  explicit TrackingInterface(bool withVertices) : mWithVertices(withVertices) {}

  /// Where the tracking runs, \a device: on a GPU backend, the GPU -- its
  /// reconstruction and the device tracker's arena (DeviceContext) -- before
  /// any timeframe; nothing on the host.
  void initDevice(TrackingDevice device)
  {
    mDeviceOptions = std::move(device);
    if (mDeviceOptions.backend.empty()) {
      return;
    }
    if (mDeviceOptions.arenaMB == 0U) {
      LOGP(fatal, "--gtrk-gpu-backend {}: --gtrk-gpu-arena-mb, the device tracker's arena, is required", mDeviceOptions.backend);
      return;
    }
    const std::array<std::size_t, 1> arena{mDeviceOptions.arenaMB << 20U};
    mContext = std::make_unique<DeviceContext>(mDeviceOptions.backend, mDeviceOptions.device, arena);
    mBlock = mContext->block(0U);
    const auto& properties = mContext->runtime().properties();
    LOGP(info, "the tracking on {} device {}: {} ({}), arena {} MB", mContext->runtime().backend(), mDeviceOptions.device, properties.name, properties.architecture, mDeviceOptions.arenaMB);
  }

  /// The tracking on \a reconstruction's GPU, which its owner initialised
  /// (the GPU reconstruction workflow's): a borrowing context, and each
  /// timeframe's device memory from `useBlock`.
  void initDevice(o2::gpu::GPUReconstruction& reconstruction)
  {
    mContext = std::make_unique<DeviceContext>(reconstruction);
    mDeviceOptions = {.backend = std::string{mContext->runtime().backend()}};
    const auto& properties = mContext->runtime().properties();
    LOGP(info, "the tracking on {} inside the GPU reconstruction: {} ({}), in the pool's free memory after the TPC", mContext->runtime().backend(), properties.name, properties.architecture);
  }

  /// The next timeframe's device memory, on a borrowed reconstruction: set
  /// before its vertices and tracking, valid until its end.
  void useBlock(DeviceBlock block) { mBlock = block; }

  /// The GPU's context; a GPU backend's only.
  [[nodiscard]] const DeviceContext& context() const { return *mContext; }

  /// The GPU released, the device tracker before its context: at the end
  /// of the processing, not at the process's exit.
  void releaseDevice() noexcept
  {
    mDevice.reset();
    mContext.reset();
  }

  [[nodiscard]] bool configured() const { return mTracker != nullptr; }

  /// `configure` from \a conditions (Conditions.h's, with the material LUT)
  /// if the trackers are not built for \a shipped: on the first timeframe,
  /// and for a new beam line. Whether it did.
  template <class Conditions>
  bool configureFor(const Conditions& conditions, const MeasurementGeometry<Det>& shipped)
  {
    if (mTracker && shipped.beamLine == mBeamLine) {
      return false;
    }
    configure(trackingConditions(conditions), shipped);
    return true;
  }

  /// Place the geometry before deriving a standalone grid from its layer
  /// envelopes. A workflow can leave this to configure, using its shipped grid.
  void initGeometry(std::span<const double> depthsCM)
  {
    if (!mGeometry) {
      mGeometry = makeGeometry<Det>(depthsCM);
    }
  }

  /// On the first call the geometry, the material's budgets and scattering
  /// table and the tracking geometry, from \a conditions; on every call the
  /// trackers, for \a shipped's beam line and grid.
  void configure(const TrackingConditions& conditions, const MeasurementGeometry<Det>& shipped)
  {
    const auto& param = recoParam<Det>();
    LOGP(info, "the field: nominal Bz {} kG, from the run's GRPMagField", conditions.bzKiloGauss);
    mBzKiloGauss = conditions.bzKiloGauss;
    const bool first = mMaterial == nullptr;
    initGeometry(conditions.depthsCM);
    // The beam line and the grid the measurements are relative to, the producer's.
    mBeamLine = shipped.beamLine;
    LOGP(info, "the beam line, the trackers' constraint: {}", describe(mBeamLine));
    mGrid.emplace(shipped.grid);
    mFrames.clear();
    for (std::size_t layer = 0; layer < GridGeometry<Det>::kLayerCount; ++layer) {
      mFrames.push_back(shipped.grid.frames(static_cast<LayerID>(layer)));
    }
    if (first) {
      configureMaterial(conditions.material);
    }
    const auto beamSpot = beamSpotConstraint(mBeamLine, param);
    const auto chunkSize = static_cast<typename Controller::Index>(std::max(param.trackletChunkSize, 1));
    // The tracker takes each iteration's configured source.
    std::vector<IterationSettings> settings;
    for (int iteration = 0; iteration < iterationCount(param); ++iteration) {
      settings.push_back(iterationSettings(param, iteration));
      if (settings.back().vertexSource == VertexSource::kVertexList && !mWithVertices && !param.vertexer) {
        LOGP(fatal, "iteration {} seeds from the vertex list (GTRKReco<DET>.iterationVertexSource), and this device has none: an input of seeds or GTRKReco<DET>.vertexer", iteration);
        return;
      }
    }
    // GTRKReco<DET>.vertexer: the own vertexer's list, in place of an input one.
    if (param.vertexer && mWithVertices) {
      LOGP(fatal, "GTRKReco<DET>.vertexer finds the vertex list, and an input supplies one: take one or the other");
      return;
    }
    mVertexing.reset();
    // One pool for the vertexer and the tracker: memoryLimitGB bounds what
    // both hold at once (the vertexer's pass stays while the tracking runs).
    mPool = trackerMemoryPool(param);
    if (param.vertexer) {
      mVertexing = std::make_unique<Vertexing>(vertexingPassSettings(param), vertexerSettings(param), beamSpot, mScattering, mPool, chunkSize, std::max(param.threads, 0));
      LOGP(info, "the vertexer: a pass at {} GeV on layers 0-2, window {} um, claim {} um, at least {} contributors ({} in a frame without a vertex)", param.vertexerMinPtGeV, param.vertexerWindowUM, param.vertexerClaimUM, param.vertexerMinContributors, param.vertexerAloneMinContributors);
    }
    mTrackerSource = iterationSettings(param, 0).vertexSource;
    mTrackerLabel = sourceLabel(settings);
    mTracker = std::make_unique<Controller>(settings, beamSpot, mScattering, mPool, chunkSize, std::max(param.threads, 0));
    if (param.memoryLimitGB > 0.) {
      LOGP(info, "the tracker's memory, its vertexer's included, limited to {:.3f} GB (GTRKReco<DET>.memoryLimitGB)", param.memoryLimitGB);
    }
    // The road fit and the refit: the sensors, the material LUT, the field.
    if (!mTracking) {
      mTracking = std::make_unique<TrackingGeometry<Det>>(*mGeometry);
    }
    mTracker->setTrackModel(*mTracking, *mMaterial, propagator());
    if (mContext) {
      configureDevice(settings);
    }
  }

  /// The input's checks: every column \a columns as long as the others, the
  /// grid's prefix \a offsets over its \a bucketCount buckets ending at the
  /// row count, every measurement's variances fit for the fit, \a seeds
  /// time-sorted with finite z. Nullptr if it may be tracked, else why not.
  template <class Columns>
  [[nodiscard]] static const char* checkInput(const Columns& columns, std::span<const std::uint32_t> offsets, std::size_t bucketCount, std::span<const VertexSeed> seeds)
  {
    if (!std::ranges::is_sorted(seeds, {}, &VertexSeed::timePS) || !std::ranges::all_of(seeds, [](const VertexSeed& seed) { return std::isfinite(seed.z); })) {
      return "vertex seeds are not time-sorted or finite";
    }
    const auto n = columns.loc0.size();
    const bool validOffsets = n <= std::numeric_limits<std::uint32_t>::max() && offsets.size() == bucketCount + 1U && offsets.front() == 0U && offsets.back() == n && std::ranges::is_sorted(offsets);
    if (!measurementColumnsHaveSize(columns, n) || !validOffsets) {
      return "measurement column lengths or grid offsets are invalid";
    }
    for (std::size_t row = 0; row < n; ++row) {
      if (!measurementVariancesValid(columns.varLoc0[row], columns.covLoc01[row], columns.varLoc1[row])) {
        return "a measurement has a variance the fit cannot take: not finite, not positive, or a singular spatial block";
      }
    }
    return nullptr;
  }

  /// The timeframe's vertex list: \a input, or the own vertexer's, which
  /// keeps its triplets' claims with \a keepClaims. Valid until endTimeframe.
  struct Vertices {
    std::span<const VertexSeed> seeds;
    const VertexerResult* found = nullptr;            ///< the own vertexer's list and contributors, host or device; none without it
    const typename Vertexing::Result* pass = nullptr; ///< the own vertexer's on the host, with its claims; none without it or on a GPU
    std::span<const TripletClaim> claims;             ///< the own vertexer's claimed triplets when requested
    double seconds = 0.;
    double planSeconds = 0.; ///< the initial device plan, before upload and vertex finding
    bool dropped = false;    ///< out of memory, dropped (dropTFUponFailure): no tracking
  };
  [[nodiscard]] Vertices vertices(const Frame& frame, std::span<const VertexSeed> input, bool keepClaims)
  {
    if (!mVertexing) {
      return {.seeds = input};
    }
    const auto& param = recoParam<Det>();
    const auto start = std::chrono::steady_clock::now();
    mVertexed = true;
    if (mDevice) {
      // The run's every iteration, the pass last, with an empty list.
      mDevicePlan = planTimeframe<Det>(frame, mDeviceIterations, mDeviceWindows, mTracker->beamLine(), {}, frame.grid().geometry().blocksPerSlice(param.sliceWidthPS), param.sliceHaloPS);
      const double planSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      if (!beginDeviceTimeframe()) {
        return {.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(), .planSeconds = planSeconds, .dropped = true};
      }
      const auto status = mDevice->findTimeframeVertices(frame, mDevicePlan, mDeviceIterations.size() - 1U, mVertexing->settings(), beamSpotConstraint(mBeamLine, param), mDeviceFound, mDeviceReport, keepClaims ? &mDeviceClaims : nullptr);
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      if (status != StageStatus::kComplete) {
        return {.seconds = seconds, .planSeconds = planSeconds, .dropped = dropOrFail(status, "the vertexing pass")};
      }
      return {.seeds = mDeviceFound.vertices, .found = &mDeviceFound, .claims = keepClaims ? std::span<const TripletClaim>{mDeviceClaims} : std::span<const TripletClaim>{}, .seconds = seconds, .planSeconds = planSeconds};
    }
    mPass.emplace(mVertexing->find(frame, frame.grid().geometry().blocksPerSlice(param.sliceWidthPS), param.sliceHaloPS, keepClaims));
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (mPass->status != StageStatus::kComplete) {
      return {.pass = &*mPass, .seconds = seconds, .dropped = dropOrFail(mPass->status, "the vertexing pass")};
    }
    return {.seeds = mPass->found.vertices, .found = &mPass->found, .pass = &*mPass, .claims = mPass->claims, .seconds = seconds};
  }

  /// The tracker's run on \a frame from \a seeds. Valid until endTimeframe.
  struct Tracked {
    TimeframeResult result;
    double seconds = 0.;
    double planSeconds = 0.; ///< included in seconds: the tracking plan, before the run
    bool dropped = false;    ///< out of memory, dropped (dropTFUponFailure): no tracks
  };
  [[nodiscard]] Tracked track(const Frame& frame, std::span<const VertexSeed> seeds)
  {
    const auto& param = recoParam<Det>();
    const auto blocksPerSlice = frame.grid().geometry().blocksPerSlice(param.sliceWidthPS);
    const auto start = std::chrono::steady_clock::now();
    if (mDevice) {
      return trackOnDevice(frame, seeds, blocksPerSlice, start);
    }
    mTracked = true;
    const auto plan = mTracker->planTimeframe(frame, seeds, blocksPerSlice, param.sliceHaloPS);
    const double planSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    Tracked tracked{.result = mTracker->run(frame, plan), .planSeconds = planSeconds};
    tracked.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (tracked.result.status != StageStatus::kComplete) {
      tracked.dropped = dropOrFail(tracked.result.status, "the tracking");
    }
    return tracked;
  }

  /// The end of the timeframe: what began, the vertexer's and the tracker's
  /// storage, back to the pool, the staging learnt. Its results are invalid
  /// after.
  void endTimeframe() noexcept
  {
    mLastArenaPeak = 0U;
    if (mDeviceTimeframe) {
      mDevice->endTimeframe();
      mLastArenaPeak = mDevice->arenaPeakBytes();
      mDevice->detach();
      mDeviceTracks.setMemoryResource(mDevicePool.get());
      mDevicePlan = {};
      mDeviceTrackingPlan = {};
      mDeviceFound = {};
      mDeviceClaims.clear();
      mDeviceTimeframe = false;
      mVertexed = false;
    }
    if (mTracked) {
      mTracker->endTimeframe();
      mTracked = false;
    }
    if (mVertexed && !mDevice) {
      mVertexing->endTimeframe();
    }
    mVertexed = false;
    mPass.reset();
  }

  /// Calls endTimeframe when it goes out of scope.
  class TimeframeGuard
  {
   public:
    explicit TimeframeGuard(TrackingInterface& tracking) : mTracking(&tracking) {}
    TimeframeGuard(const TimeframeGuard&) = delete;
    TimeframeGuard& operator=(const TimeframeGuard&) = delete;
    TimeframeGuard(TimeframeGuard&&) = delete;
    TimeframeGuard& operator=(TimeframeGuard&&) = delete;
    ~TimeframeGuard() { mTracking->endTimeframe(); }

   private:
    TrackingInterface* mTracking;
  };

  /// The guard of one timeframe, taken at its top.
  [[nodiscard]] TimeframeGuard endTimeframeOnExit() { return TimeframeGuard{*this}; }

  /// The timeframes dropped so far (dropTFUponFailure).
  [[nodiscard]] std::size_t droppedTimeframes() const { return mDropped; }

  [[nodiscard]] const std::shared_ptr<BoundedMemoryResource>& pool() const { return mPool; }
  [[nodiscard]] Vertexing* vertexing() { return mVertexing.get(); }
  [[nodiscard]] const typename Device::TimeframeReport& deviceReport() const { return mDeviceReport; }
  [[nodiscard]] const DeviceProperties& deviceProperties() const { return mContext->runtime().properties(); }
  [[nodiscard]] std::size_t deviceArenaPeakBytes() const { return mDevice->arenaPeakBytes(); }
  /// The device arena's peak in the last timeframe that ended, 0 before one.
  [[nodiscard]] std::size_t lastArenaPeakBytes() const { return mLastArenaPeak; }

  [[nodiscard]] Controller& tracker() { return *mTracker; }
  [[nodiscard]] const Controller& tracker() const { return *mTracker; }
  [[nodiscard]] const Geometry<Det>& geometry() const { return *mGeometry; }
  [[nodiscard]] const TrackingGeometry<Det>& tracking() const { return *mTracking; }
  [[nodiscard]] const Material& material() const { return *mMaterial; }
  [[nodiscard]] const ScatteringTable<Det>& scattering() const { return mScattering; }
  [[nodiscard]] const BeamLine& beamLine() const { return mBeamLine; }
  [[nodiscard]] const GridGeometry<Det>& grid() const
  {
    if (!mGrid) {
      throw std::logic_error("GenericTracking: the grid is asked for before the producer shipped it");
    }
    return *mGrid;
  }
  [[nodiscard]] const std::vector<FrameTiming>& frames() const { return mFrames; }
  /// The tracker's first iteration's source, and its sources as the log lines name them.
  [[nodiscard]] VertexSource trackerSource() const { return mTrackerSource; }
  [[nodiscard]] const std::string& trackerLabel() const { return mTrackerLabel; }
  /// The track model's propagator: the run's field, the mass hypothesis.
  [[nodiscard]] Propagator propagator() const { return {.bzKiloGauss = static_cast<Scalar>(mBzKiloGauss), .massGeV = static_cast<Scalar>(recoParam<Det>().massHypothesisGeV)}; }

 private:
  /// The layers' and the passive scatterers' budgets from \a material and
  /// the scattering table: once per run.
  void configureMaterial(const Material* material)
  {
    mMaterial = material;
    if (mMaterial == nullptr) {
      LOGP(fatal, "no material LUT from the CCDB ({})", materialLutPath<Det>());
      return;
    }
    const auto budget = layerMaterial<Det>(*mMaterial);
    const auto passiveBudget = passiveMaterial<Det>(*mMaterial);
    mScattering = scatteringTable<Det>(budget, passiveBudget);
    std::string layers;
    for (std::size_t layer = 0; layer < budget.size(); ++layer) {
      layers += fmt::format(" {}: {:.3f} %", layer, 100. * budget[layer]);
    }
    LOGP(info, "layer budgets, x/X0 across each layer's band (LUT):{}", layers);
    std::size_t passiveIndex = 0U;
    for (const double x2X0 : passiveBudget) {
      const auto& passive = mScattering.passives[passiveIndex];
      LOGP(info, "passive scatterer {}: {} at {:.3f} cm, x/X0 {:.3e} (LUT)", passiveIndex++, passive.barrel ? "cylinder" : "disc", passive.at.position, x2X0);
    }
  }

  /// The device tracker, made once, and the run it is configured with every
  /// timeframe (beginDeviceTimeframe): the host tracker's iterations, \a
  /// settings their windows, with the pass last if the own vertexer runs. On
  /// the first timeframe and on every new beam line.
  void configureDevice(const std::vector<IterationSettings>& settings)
  {
    if (!mDevice) {
      mDevice = mContext->runtime().makeTracker<Det>(DeviceBlock{});
      mDevicePool = trackerMemoryPool(recoParam<Det>());
      mDeviceTracks.setMemoryResource(mDevicePool.get());
    }
    mDeviceIterations.clear();
    mDeviceWindows.clear();
    for (std::size_t iteration = 0; iteration < mTracker->iterationCount(); ++iteration) {
      mDeviceIterations.push_back(mTracker->iteration(iteration));
      mDeviceWindows.push_back(settings[iteration].windows);
    }
    if (mVertexing) {
      mDeviceIterations.push_back(mVertexing->passParameters());
      mDeviceWindows.push_back(vertexingPassSettings(recoParam<Det>()).windows);
    }
  }

  /// The timeframe's device work begun, once per timeframe: the block
  /// attached and the run configured on it. False if the timeframe is dropped:
  /// the block does not hold the run (dropTFUponFailure).
  [[nodiscard]] bool beginDeviceTimeframe()
  {
    if (mDeviceTimeframe) {
      return true;
    }
    mDeviceTimeframe = true; // endTimeframe detaches
    mDevice->attach(mBlock);
    const auto status = mDevice->configureRun(mDeviceIterations, *mTracking, *mMaterial, propagator());
    return status == StageStatus::kComplete || !dropOrFail(status, "the run's data on the device");
  }

  /// The tracking on the device: with the own vertexer, the tracking's plan
  /// with its list replacing the timeframe's; else the timeframe begun with
  /// the tracking's plan from \a seeds.
  [[nodiscard]] Tracked trackOnDevice(const Frame& frame, std::span<const VertexSeed> seeds, std::int32_t blocksPerSlice, std::chrono::steady_clock::time_point start)
  {
    const auto& param = recoParam<Det>();
    const bool replaced = mVertexing != nullptr;
    const std::size_t iterations = mTracker->iterationCount();
    mDeviceTrackingPlan = replaced ? planTimeframe<Det>(frame, mDeviceIterations, mDeviceWindows, mTracker->beamLine(), seeds, blocksPerSlice, param.sliceHaloPS) : mTracker->planTimeframe(frame, seeds, blocksPerSlice, param.sliceHaloPS);
    Tracked tracked;
    tracked.planSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (!beginDeviceTimeframe()) {
      tracked.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      tracked.dropped = true;
      return tracked;
    }
    tracked.result.status = mDevice->trackTimeframe(frame, mDeviceTrackingPlan, iterations, replaced, mDeviceTracks, mDeviceReport);
    tracked.result.tracks = &mDeviceTracks;
    tracked.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (tracked.result.status != StageStatus::kComplete) {
      tracked.dropped = dropOrFail(tracked.result.status, "the tracking");
    }
    return tracked;
  }

  /// A timeframe whose \a what ended with \a status: dropped, if it ran out
  /// of memory and dropTFUponFailure is on; else the run fails.
  bool dropOrFail(StageStatus status, const char* what)
  {
    if (status == StageStatus::kOutOfMemory && recoParam<Det>().dropTFUponFailure) {
      ++mDropped;
      LOGP(error, "{} ran out of memory (GTRKReco<DET>.memoryLimitGB, or the device's --gtrk-gpu-arena-mb): the timeframe dropped, {} so far (GTRKReco<DET>.dropTFUponFailure)", what, mDropped);
      return true;
    }
    LOGP(fatal, "{} failed with status {}", what, static_cast<int>(status));
    return false;
  }

  bool mWithVertices = false;
  double mBzKiloGauss = 0.;
  std::unique_ptr<const Geometry<Det>> mGeometry;
  const Material* mMaterial = nullptr; ///< the condition's (MaterialLUT.h)
  ScatteringTable<Det> mScattering;
  std::unique_ptr<TrackingGeometry<Det>> mTracking;
  BeamLine mBeamLine;                     ///< the producer's, the trackers were built for
  std::optional<GridGeometry<Det>> mGrid; ///< the grid's layout, the producer's (MeasurementGeometry)
  std::vector<FrameTiming> mFrames;       ///< each layer's readout frames, the grid's
  std::shared_ptr<BoundedMemoryResource> mPool;
  std::unique_ptr<Controller> mTracker;  ///< each iteration from its configured source (GTRKReco<DET>.iterationVertexSource)
  std::unique_ptr<Vertexing> mVertexing; ///< the own vertexer before the tracking (GTRKReco<DET>.vertexer); none without
  std::optional<typename Vertexing::Result> mPass;
  VertexSource mTrackerSource = VertexSource::kBeamLine;
  std::string mTrackerLabel;
  TrackingDevice mDeviceOptions;
  std::unique_ptr<DeviceContext> mContext;                 ///< the GPU, on a GPU backend: its reconstruction, the arena's block and the runtime
  DeviceBlock mBlock;                                      ///< the device memory each timeframe attaches
  std::unique_ptr<Device> mDevice;                         ///< the device tracker, on the context; released before it
  std::vector<IterationParameters<Det>> mDeviceIterations; ///< as configured on the device: the host tracker's, then the pass
  std::vector<SearchWindowParams> mDeviceWindows;          ///< the same iterations' windows
  /// The timeframe's, which the device reads until endTimeframe: the plan
  /// of the vertexer's pass, the tracking's plan, the vertexer's list, the
  /// tracks downloaded into the pool, the calls' report.
  TimeframePlan<Det> mDevicePlan;
  TimeframePlan<Det> mDeviceTrackingPlan;
  std::vector<TripletClaim> mDeviceClaims; ///< optional MC contributor identities
  VertexerResult mDeviceFound;
  std::shared_ptr<BoundedMemoryResource> mDevicePool;
  typename TrackRefit<Det>::Tracks mDeviceTracks;
  typename Device::TimeframeReport mDeviceReport;
  bool mDeviceTimeframe = false; ///< the device holds a timeframe
  bool mVertexed = false;        ///< the vertexer ran this timeframe
  bool mTracked = false;         ///< the tracker ran this timeframe
  std::size_t mDropped = 0U;
  std::size_t mLastArenaPeak = 0U;
};

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_TRACKINGINTERFACE_H_
