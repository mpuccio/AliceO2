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
/// \file TrackingInterface.cxx
/// \brief
///

#include "ITSMFTTracking/DetectorTrackingOperationAdapterSupport.h"
#include "ITSMFTTracking/TrackingInterface.h"

#include <array>
#include <algorithm>
#include <limits>
#include <type_traits>

#include "ITStracking/Constants.h"

#include <oneapi/tbb/task_arena.h>

#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/IRFrame.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/TimeFrameLoadFailure.h"

#include <format>

namespace o2::itsmft::tracking
{

namespace
{
template <o2::detectors::DetID::ID detId>
consteval const char* detName()
{
  return detId == o2::detectors::DetID::MFT ? "MFT" : "ITS";
}

// Both kITSStaticSurfaceCatalog and kMFTStaticSurfaceCatalog are dense and
// local to their own detector (StaticDetectorCatalogs.h / ITSSurfaceSpec.h /
// MFTSurfaceSpec.h): surface i's id is always SurfaceId{i}. The owner's
// legacy-layer traversal order is therefore always the trivial identity
// sequence, for either detector -- no external/runtime value to supply.
template <int NLayers>
constexpr std::array<SurfaceId, NLayers> identitySurfaceOrder()
{
  std::array<SurfaceId, NLayers> order{};
  for (int i = 0; i < NLayers; ++i) {
    order[i] = SurfaceId{static_cast<uint16_t>(i)};
  }
  return order;
}

} // namespace

#ifndef GPUCA_GPUCODE
template <int NLayers>
ITSMFTTrackingInterface<NLayers>::ITSMFTTrackingInterface(bool useMC,
                                                          o2::itsmft::TrackingMode::Type mode,
                                                          bool overrideBeamEst)
  : ITSMFTTrackingInterface(useMC, mode, overrideBeamEst, nullptr)
{
}

template <int NLayers>
ITSMFTTrackingInterface<NLayers>::ITSMFTTrackingInterface(bool useMC,
                                                          o2::itsmft::TrackingMode::Type mode,
                                                          bool overrideBeamEst,
                                                          std::unique_ptr<ClusterDecoder> clusterDecoder)
  : mUseMC(useMC), mOverrideBeamEstimation(overrideBeamEst), mTrackingMode(mode), mClusterDecoder(clusterDecoder != nullptr ? std::move(clusterDecoder) : std::make_unique<GeometryClusterDecoder<DetId>>())
{
}
#else
template <int NLayers>
ITSMFTTrackingInterface<NLayers>::ITSMFTTrackingInterface(bool useMC,
                                                          o2::itsmft::TrackingMode::Type mode,
                                                          bool overrideBeamEst)
  : mUseMC(useMC), mOverrideBeamEstimation(overrideBeamEst), mTrackingMode(mode)
{
}
#endif

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::initialise()
{
  const auto parameters = resolveTrackingParameters();
  initialiseTracker(parameters);
}

template <int NLayers>
std::vector<o2::itsmft::TrackingParameters> ITSMFTTrackingInterface<NLayers>::resolveTrackingParameters()
{
  auto mode = mTrackingMode;
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    // The common ITS interface must not consult TrackerParamRef<ITS>::get()
    // (o2::its::TrackerParamConfig, the frozen legacy "ITSCATrackerParam"
    // namespace -- see TrackingConfigParam.h's doc comment and this class's
    // own nThreads isolation in initialiseTracker() below for the same
    // rationale) for its tracking mode. The common workflow
    // constructor/setTrackingMode() is the sole owner of mTrackingMode; only
    // its Unset default is resolved here.
    if (mode == o2::itsmft::TrackingMode::Unset) {
      mode = o2::itsmft::TrackingMode::Sync;
    }
  } else {
    const auto& tc = o2::itsmft::tracking::TrackerParamRef<DetId>::get();
    if (auto parmode = static_cast<o2::itsmft::TrackingMode::Type>(tc.trackingMode);
        mode == o2::itsmft::TrackingMode::Unset ||
        (parmode != o2::itsmft::TrackingMode::Unset && mode != parmode)) {
      if (parmode != o2::itsmft::TrackingMode::Unset) {
        LOGP(info, "{} CA tracking mode overwritten by configurable params from {} to {}",
             detName<DetId>(), o2::itsmft::TrackingMode::toString(mode), o2::itsmft::TrackingMode::toString(parmode));
        mode = parmode;
      }
    }
    if (mode == o2::itsmft::TrackingMode::Unset) {
      mode = o2::itsmft::TrackingMode::Sync;
    }
  }
  mTrackingMode = mode;
  auto trackParams = o2::itsmft::TrackingMode::getTrackingParameters(DetId, mode);
  LOGP(info, "{} CA tracker initialized in {} mode with {} iteration(s)",
       detName<DetId>(), o2::itsmft::TrackingMode::toString(mode), trackParams.size());
  for (size_t i = 0; i < trackParams.size(); ++i) {
    LOGP(info, "  iteration {}: {}", i, trackParams[i].asString());
  }
  return trackParams;
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::initialiseTracker(const std::vector<o2::itsmft::TrackingParameters>& trackParams)
{
  if (trackParams.empty()) {
    return;
  }
  mTrackerTraits = std::make_unique<TrackerTraits>();
  std::shared_ptr<tbb::task_arena> taskArena;

  int nThreads;
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    // The common ITS interface must consume its own dedicated
    // ITSCommonCATrackerParam.nThreads, never TrackerParamRef<ITS>::get()
    // (o2::its::TrackerParamConfig, the frozen legacy "ITSCATrackerParam"
    // namespace -- see TrackingConfigParam.h's doc comment on why
    // ITSCommonCATrackerParam is kept a deliberately distinct type). A
    // non-positive value is a construction-time misconfiguration, not
    // per-TF data, so it is rejected here, once, before a tbb::task_arena
    // is ever built from it -- TrackerTraits::setNThreads silently
    // std::abs()-es its argument, which is MFT's pre-existing, untouched
    // behavior via TrackerParamRef<MFT> below, not something ITS should
    // rely on for its own dedicated knob.
    const auto& itsCommonParam = o2::itsmft::ITSCommonCATrackerParam::Instance();
    if (itsCommonParam.nThreads <= 0) {
      LOGP(fatal, "ITS CA tracker requires ITSCommonCATrackerParam.nThreads > 0, got {}", itsCommonParam.nThreads);
    }
    nThreads = itsCommonParam.nThreads;
  } else {
    nThreads = o2::itsmft::tracking::TrackerParamRef<DetId>::get().nThreads;
  }
  mTrackerTraits->setNThreads(nThreads, taskArena);
  mTracker = std::make_unique<Tracker>(mTrackerTraits.get());
  static constexpr auto kOrderedSurfaces = identitySurfaceOrder<NLayers>();
  TrackerInitialization configuration;
  std::shared_ptr<BoundedMemoryResourceN> memoryPool;
  size_t maxMemory = std::numeric_limits<size_t>::max();
  if (trackParams[0].MaxMemory != maxMemory) {
    maxMemory = trackParams[0].MaxMemory;
  }
  memoryPool = std::make_shared<BoundedMemoryResourceN>(maxMemory);
  configuration.memoryPool = memoryPool;
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    configuration.catalog = SurfaceCatalogView{kITSStaticSurfaceCatalog.data(), static_cast<uint32_t>(kITSStaticSurfaceCatalog.size())};
  } else {
    configuration.catalog = SurfaceCatalogView{kMFTStaticSurfaceCatalog.data(), static_cast<uint32_t>(kMFTStaticSurfaceCatalog.size())};
  }
  configuration.iterations.reserve(trackParams.size());
  for (const auto& params : trackParams) {
    TrackerIterationConfiguration iteration;
    iteration.parameters = {params};
    SurfaceGraphSubgraph subgraph;
    subgraph.orderedSurfaces.assign(kOrderedSurfaces.begin(), kOrderedSurfaces.end());
    subgraph.maxHoles = params.MaxHoles;
    subgraph.holeSurfaces = positionalSurfaceMask(params.HoleLayerMask, kOrderedSurfaces, NLayers);
    subgraph.seedingSurfaces = positionalSurfaceMask(params.StartLayerMask, kOrderedSurfaces, NLayers);
    iteration.graphSubgraphs.push_back(subgraph);
    SurfaceMask owned;
    for (uint16_t i = 0; i < NLayers; ++i) {
      owned.set(SurfaceId{i});
    }
    constexpr SurfaceKind expectedKind = (DetId == o2::detectors::DetID::ITS) ? SurfaceKind::Cylinder : SurfaceKind::Disk;
    constexpr TransitionPolicyTag expectedPolicy = (DetId == o2::detectors::DetID::ITS) ? TransitionPolicyTag::CylinderCylinder : TransitionPolicyTag::DiskDisk;
    iteration.bindings.push_back(SurfacePlanBinding::Declaration{ClusterSourceId{0}, owned,
                                                                 std::vector<SurfaceId>{kOrderedSurfaces.begin(), kOrderedSurfaces.end()},
                                                                 expectedKind, expectedPolicy});
    configuration.iterations.push_back(std::move(iteration));
  }
  mTracker->setSource(ClusterSourceId{0});
  const auto result = mTracker->initialize(mFrame, configuration);
  if (!result.ok()) {
    LOGP(fatal, "{} CA tracker failed to initialize static configuration (error={} iteration={} graph={} binding={})",
         detName<DetId>(), static_cast<int>(result.error), result.failedIteration,
         static_cast<int>(result.graphError), static_cast<int>(result.bindingError));
  }
  mTrackerTraits->setMemoryPool(mFrame.getMemoryPool());
}

template <int NLayers>
float ITSMFTTrackingInterface<NLayers>::processTimeFrame(gsl::span<const o2::itsmft::ROFRecord> rofs,
                                                         gsl::span<const o2::itsmft::CompClusterExt> clusters,
                                                         gsl::span<const unsigned char> patterns,
                                                         const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels,
                                                         gsl::span<const o2::dataformats::IRFrame> irFrames,
                                                         RuntimeROFViews rofViews)
{
  if (mFrame.getTrackingParameters().empty()) {
    LOGP(info, "{} CA tracking mode is off, skipping TimeFrame processing", detName<DetId>());
    return 0.f;
  }
  if (rofViews.overlap.mLayerCount <= 0 && mBoundROFViews != nullptr) {
    rofViews = *mBoundROFViews;
  }

  // mDict == nullptr is a structural configuration failure (the dictionary
  // was never configured on this interface), not per-TF data -- it belongs
  // inside the same failure boundary as every other loading-phase check
  // below, not as a bare LOGP(fatal,...) ahead of it, so that it wipes
  // TimeFrame state and propagates through the identical typed contract
  // rather than a separate, uncaught fatal path.
  try {
    if (mDict == nullptr) {
      throw TimeFrameLoadException{TimeFrameLoadFailureReason::DictionaryNotConfigured,
                                   std::format("{} CA tracker cluster dictionary is not available", detName<DetId>())};
    }
    loadTimeFrame(rofs, clusters, patterns, labels, irFrames, rofViews);
  } catch (const RecoverableLoadFailure& err) {
    // Typed, per-TF malformed loading input (see TimeFrameLoadFailure.h /
    // isRecoverableLoadError()).
    LOGP(error, "{} CA loading recoverably failed: {}", detName<DetId>(), err.what());
    resetEvent();
    if (mFrame.getTrackingParameters()[0].DropTFUponFailure) {
      return kDroppedTimeFrameResult;
    }
    throw;
  } catch (const BoundedMemoryResourceN::MemoryLimitExceeded& err) {
    // Recoverable, per-TF resource failure during loading -- same
    // classification/gating as the identical catch clause in
    // Tracker::clustersToTracks() (Tracker.cxx).
    LOGP(error, "{} CA loading exceeded memory limit: {}", detName<DetId>(), err.what());
    resetEvent();
    if (mFrame.getTrackingParameters()[0].DropTFUponFailure) {
      return kDroppedTimeFrameResult;
    }
    throw;
  } catch (const std::bad_alloc& err) {
    LOGP(error, "{} CA loading allocation failed: {}", detName<DetId>(), err.what());
    resetEvent();
    if (mFrame.getTrackingParameters()[0].DropTFUponFailure) {
      return kDroppedTimeFrameResult;
    }
    throw;
  } catch (const TimeFrameLoadException& err) {
    // Structural loading-boundary failure: never gated by DropTFUponFailure.
    LOGP(error, "{} CA loading hit a structural failure: {}", detName<DetId>(), err.what());
    resetEvent();
    throw;
  } catch (const std::exception& err) {
    // Unclassified: not a recognized recoverable-resource or recoverable
    // loading-data failure, so treated as structural regardless of
    // DropTFUponFailure, and rethrown by its original type (not wrapped).
    LOGP(error, "{} CA loading failed with an unclassified exception; treating as structural: {}", detName<DetId>(), err.what());
    resetEvent();
    throw;
  }
  onTimeFrameLoaded();

  const float elapsedMs = runTracking();
  onTrackingFinished(elapsedMs);

  return elapsedMs;
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::loadTimeFrame(gsl::span<const o2::itsmft::ROFRecord> rofs,
                                                     gsl::span<const o2::itsmft::CompClusterExt> clusters,
                                                     gsl::span<const unsigned char> patterns,
                                                     const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels,
                                                     gsl::span<const o2::dataformats::IRFrame> irFrames,
                                                     RuntimeROFViews rofViews)
{
  validateROFInput(rofs);

  mFrame.setBz(o2::base::Propagator::Instance()->getNominalBz());
  configureBeamPosition();

  // origin: the TF-relative BC anchor every source ROF is converted against
  // (Architecture.md 7.2). The first ROF's own real interaction record, or
  // the explicit default InteractionRecord{} when there are no ROFs at all.
  const o2::InteractionRecord origin = rofs.empty() ? o2::InteractionRecord{} : rofs.front().getBCData();

  const auto* binding = mFrame.getBinding(0, ClusterSourceId{0});
  if (binding == nullptr) {
    throw TimeFrameLoadException{TimeFrameLoadFailureReason::DictionaryNotConfigured,
                                 "CA tracker has no configured source binding"};
  }
  const auto orderedSurfaces = binding->getOrderedSurfaces();
  ClusterSourceInput source;
  source.id = ClusterSourceId{0};
  source.detector = DetId;
  source.clusters = clusters;
  source.patterns = patterns;
  source.rofs = rofs;
  source.dictionary = mDict;
  source.labels = labels;
  source.layerToSurface = orderedSurfaces;
  if (rofViews.overlap.mLayerCount <= 0) {
    rofViews = mFrame.getWorkspace(ClusterSourceId{0}).getROFViews();
  }
  if (rofViews.overlap.mLayerCount <= 0) {
    throw TimeFrameLoadException{TimeFrameLoadFailureReason::NonUniformROFTiming,
                                 "CA tracker received no adapter-owned runtime ROF timing view"};
  }
  const auto& clock = rofViews.overlap.getLayer(0);
  source.timing = ROFTimingConfig{clock.mROFLength, clock.mROFDelay, clock.mROFBias, clock.mROFAddTimeErr};
  source.decoder = mClusterDecoder.get();
  source.rofViews = rofViews;
  const auto result = MultiSourceTimeFrameLoader::load(
    mFrame, gsl::span<const ClusterSourceInput>{&source, 1}, mFrame.getGraph(0).getSurfaceCatalog(), origin);
  if (!result.ok()) {
    if (isRecoverableLoadError(result.error, result.timingDetail)) {
      throw RecoverableLoadFailure{result};
    }
    throw TimeFrameLoadException{result};
  }
  const auto& scratch = mFrame.getWorkspace(ClusterSourceId{0});

  LOGP(info, "{} CA loaded {} clusters from {} ROFs into TimeFrame ({} pattern bytes, MC={})",
       detName<DetId>(), scratch.getTotalClusters(), rofs.size(), patterns.size(), labels != nullptr);

  for (std::size_t iLayer = 0; iLayer < orderedSurfaces.size(); ++iLayer) {
    LOGP(info, "  layer {}: {} ROF slots", iLayer, scratch.getNrof(iLayer));
  }
}

template <int NLayers>
float ITSMFTTrackingInterface<NLayers>::runTracking()
{
  if (!mTracker || !mFrame.isConfigured() || mFrame.getTrackingParameters().empty()) {
    return 0.f;
  }
  mTracker->setBz(mFrame.getBz());
  // Gate 4 C2 Slice 2: clustersToTracks() now returns a typed TrackingResult
  // instead of a float+sentinel; this is the sole place that maps it back to
  // this interface's own external float/kDroppedTimeFrameResult contract
  // (still consumed by CATrackerSpec.cxx/its MFT counterpart via
  // isDroppedTimeFrame()). A structural or non-dropped-recoverable failure
  // still propagates as a thrown exception straight out of clustersToTracks()
  // -- unchanged, never caught here -- so the only outcome this mapping ever
  // observes is Success or RecoverableDropped.
  const auto result = mTracker->clustersToTracks(*this);
  if (result.outcome == TrackingOutcome::RecoverableDropped) {
    LOGP(warn, "{} CA tracking failed for this TF", detName<DetId>());
    return kDroppedTimeFrameResult;
  }
  const auto nTracks = mFrame.getCommonTracks().size();
  LOGP(info, "{} CA tracking produced {} tracks in {:.2f} ms", detName<DetId>(), nTracks, result.elapsedMs);
  return result.elapsedMs;
}

template <int NLayers>
bool ITSMFTTrackingInterface<NLayers>::refitSeed(const TrackSeed& seed,
                                                 const TrackingParameters& params,
                                                 float bz,
                                                 SurfaceTrackingScratch& scratch,
                                                 gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                                                 SurfaceCatalogView surfaceCatalog,
                                                 ClusterSourceId expectedSource,
                                                 TrackingCandidate& candidate)
{
  return detail::refitDetectorSeed<DetId>(seed, params, bz, scratch, layerMeasurements, surfaceCatalog, expectedSource, candidate);
}

template <int NLayers>
bool ITSMFTTrackingInterface<NLayers>::completeAccepted(gsl::span<const TrackingCandidate> candidates,
                                                        const TrackingParameters& params,
                                                        const SurfaceTrackingScratch& scratch,
                                                        bool final)
{
  return mDetectorPublicationAdapter != nullptr && mDetectorPublicationAdapter->completeAccepted(candidates, params, scratch, final);
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::resetAdapterState() noexcept
{
  if (mDetectorPublicationAdapter != nullptr) {
    mDetectorPublicationAdapter->reset();
  }
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::configureBeamPosition()
{
  if (mFrame.getTrackingParameters().empty()) {
    return;
  }
  const auto& p = mFrame.getTrackingParameters()[0];
  detail::configureAdapterBeamPosition<DetId>(mFrame, p, mMeanVertex, mOverrideBeamEstimation);
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::validateROFInput(gsl::span<const o2::itsmft::ROFRecord> rofs) const
{
  const auto views = getScratch().getROFViews();
  if (views.overlap.mLayerCount <= 0) {
    return;
  }
  const auto expectedROFsTF = views.overlap.getLayer(0).mNROFsTF;
  if (rofs.size() != expectedROFsTF) {
    LOGP(warn, "{} CA ROF count differs from continuous timing expectation: received {} expected {}",
         detName<DetId>(), rofs.size(), expectedROFsTF);
  }
}

template class ITSMFTTrackingInterface<7>;
template class ITSMFTTrackingInterface<10>;

} // namespace o2::itsmft::tracking
