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
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "ITSMFTTracking/ROFTimingUniformity.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/TimeFrameLoadFailure.h"
#include "MFTTracking/MFTTrackingParam.h"

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

bool rofOverlapsIRFrames(const o2::itsmft::ROFRecord& rof, int rofLengthInBC, gsl::span<const o2::dataformats::IRFrame> irFrames)
{
  o2::InteractionRecord rofStart{rof.getBCData()};
  o2::InteractionRecord rofEnd = rofStart + rofLengthInBC - 1;
  o2::dataformats::IRFrame ref(rofStart, rofEnd);
  for (const auto& ir : irFrames) {
    if (ir.info > 0) {
      const auto overlap = ref.getOverlap(ir);
      if (overlap.isValid()) {
        return true;
      }
    }
  }
  return false;
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
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    mDetectorPublicationAdapter.adoptITSSharedClusterCompatibility(
      &static_cast<ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar);
  }
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    mDetectorPublicationAdapter.adoptMFTPublicationCompatibility(
      &static_cast<MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar);
  }

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
                                                         gsl::span<const o2::dataformats::IRFrame> irFrames)
{
  mPublicationClock.reset();
  if (mFrame.getTrackingParameters().empty()) {
    LOGP(info, "{} CA tracking mode is off, skipping TimeFrame processing", detName<DetId>());
    return 0.f;
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
    loadTimeFrame(rofs, clusters, patterns, labels, irFrames);
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
                                                     gsl::span<const o2::dataformats::IRFrame> irFrames)
{
  // Throws TimeFrameLoadException (NonUniformROFTiming) before touching
  // the frame workspace if per-layer DPLAlpideParam values disagree; otherwise
  // configures its ROF overlap/vertex-lookup views as a side
  // effect and returns the single source-level ROFTimingConfig loadNormalizedSource() needs below.
  const ROFTimingConfig timing = configureROFLookupTables();
  validateROFInput(rofs);

  mFrame.setBz(o2::base::Propagator::Instance()->getNominalBz());
  configureBeamPosition();

  configureROFMask(rofs, irFrames);

  // origin: the TF-relative BC anchor every source ROF is converted against
  // (Architecture.md 7.2). The first ROF's own real interaction record, or
  // the explicit default InteractionRecord{} when there are no ROFs at all.
  const o2::InteractionRecord origin = rofs.empty() ? o2::InteractionRecord{} : rofs.front().getBCData();

  // Transactional normalized-source loading spans both mFrame and its
  // frame-owned workspace: on any failure below, mFrame's normalized frame and every workspace
  // compatibility container are left exactly as
  // they were before this call.
  const auto* binding = mFrame.getBinding(0, ClusterSourceId{0});
  if (binding == nullptr) {
    throw TimeFrameLoadException{TimeFrameLoadFailureReason::DictionaryNotConfigured,
                                 "CA tracker has no configured source binding"};
  }
  const auto orderedSurfaces = binding->getOrderedSurfaces();
  auto& scratch = mFrame.getWorkspace(ClusterSourceId{0});
  const auto result = scratch.loadNormalizedSource(mFrame, *mClusterDecoder, origin, timing, clusters, patterns, rofs, mDict, labels, DetId,
                                                     orderedSurfaces, mFrame.getGraph(0).getSurfaceCatalog());
  if (!result.ok()) {
    if (isRecoverableLoadError(result.error, result.timingDetail)) {
      throw RecoverableLoadFailure{result};
    }
    throw TimeFrameLoadException{result};
  }
  // Copy the finalized clock layer while scratch is live. The export is a
  // value boundary for workflows, never a retained overlap-table reference.
  mPublicationClock.emplace(scratch.getROFOverlapView().getClockLayer());
  // A successful normalized-frame replacement invalidates every CommonTrack
  // index, so clear detector-local compatibility entries in the same owner-
  // level operation. Failed loads return above without changing either.
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    static_cast<ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
  }
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    static_cast<MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
  }

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
  return mDetectorPublicationAdapter.completeAccepted(candidates, params, scratch, final);
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::resetAdapterState() noexcept
{
  mDetectorPublicationAdapter.reset();
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    static_cast<ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
  }
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    static_cast<MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
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
ROFTimingConfig ITSMFTTrackingInterface<NLayers>::configureROFLookupTables()
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    const bool continuous = o2::base::GRPGeomHelper::instance().getGRPECS()->isDetContinuousReadOut(DetId);
    LOGP(info, "{} CA tracker RO: continuous = {}", detName<DetId>(), continuous);
    mMFTTriggered = !continuous;
    const auto& alpParams = o2::itsmft::DPLAlpideParam<DetId>::Instance();
    if (mMFTTriggered) {
      mMFTROFrameLengthInBC = std::max(1, static_cast<int>(alpParams.roFrameLengthTrig / (o2::constants::lhc::LHCBunchSpacingNS * 1e3)));
    } else {
      mMFTROFrameLengthInBC = alpParams.roFrameLengthInBC;
    }
  }

  const auto& par = o2::itsmft::DPLAlpideParam<DetId>::Instance();
  const int nOrbitsPerTF = o2::base::GRPGeomHelper::getNHBFPerTF();

  // Per-layer LayerTiming, and the timing-validation this whole boundary
  // exists for, are both completed before anything below constructs or
  // commits ROFOverlapTable/ROFVertexLookupTable and frame-workspace views -- so a
  // rejected configuration never leaves TimeFrame partially updated (see
  // loadTimeFrame()'s mutation-inventory contract).
  //
  // A non-positive configured ROF length would divide by zero computing
  // nROFsPerOrbit below (o2::detectors::ID::DPLAlpideParam's
  // roFrameLayerLengthInBC/roFrameLengthInBC are plain runtime-configurable
  // ints, so this is reachable through misconfiguration, not just a
  // theoretical type-level concern); checked per layer, before that
  // division, and reported through the same typed structural failure as
  // every other timing problem in this boundary.
  std::array<o2::its::LayerTiming, NLayers> layerTimings{};
  for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
    const auto rofLengthInBC = par.getROFLengthInBC(iLayer);
    if (rofLengthInBC <= 0) {
      throw TimeFrameLoadException{
        TimeFrameLoadFailureReason::NonUniformROFTiming,
        std::format("{} CA per-layer ROF timing configuration has a non-positive ROF length ({}) on layer {}",
                    detName<DetId>(), rofLengthInBC, iLayer)};
    }
    const unsigned int nROFsPerOrbit = o2::constants::lhc::LHCMaxBunches / static_cast<unsigned int>(rofLengthInBC);
    layerTimings[iLayer] = o2::its::LayerTiming{
      .mNROFsTF = nROFsPerOrbit * static_cast<unsigned int>(nOrbitsPerTF),
      .mROFLength = static_cast<uint32_t>(rofLengthInBC),
      .mROFDelay = static_cast<uint32_t>(par.getROFDelayInBC(iLayer)),
      .mROFBias = static_cast<uint32_t>(par.getROFBiasInBC(iLayer)),
      .mROFAddTimeErr = mFrame.getTrackingParameters().empty()
                          ? static_cast<uint32_t>(o2::itsmft::tracking::TrackerParamRef<DetId>::get().addTimeError[iLayer])
                          : mFrame.getTrackingParameters()[0].AddTimeError[iLayer]};
    // A zero ROF count (rofLengthInBC > LHCMaxBunches makes nROFsPerOrbit
    // integer-divide to 0, or a misconfigured nOrbitsPerTF <= 0) leaves no
    // real ROF to anchor a diamond-vertex TF interval envelope on
    // (TrackerTraits::computeLayerTrackletsForPolicy indexes ROF 0 and ROF
    // mNROFsTF-1); reject here, at the same structural loading boundary as
    // every other timing misconfiguration, rather than let it reach that
    // indexing later.
    if (layerTimings[iLayer].mNROFsTF == 0) {
      throw TimeFrameLoadException{
        TimeFrameLoadFailureReason::ZeroROFCount,
        std::format("{} CA per-layer ROF timing configuration yields zero ROFs per TimeFrame on layer {} (rofLengthInBC={}, nOrbitsPerTF={})",
                    detName<DetId>(), iLayer, rofLengthInBC, nOrbitsPerTF)};
    }
  }

  // TimeFrame::loadNormalizedSource() takes one source-level ROFTimingConfig,
  // but DPLAlpideParam supports genuine per-layer staggering (see
  // ROFTimingUniformity.h). Both ITS and MFT default every staggering
  // override to zero, so production configurations are uniform out of the
  // box; a divergent configuration -- including one that only disagrees on
  // mNROFsTF -- is rejected here, structurally, through the same typed
  // failure as every other timing problem in this boundary, never a bare
  // LOGP(fatal, ...) and never silently collapsed. A separate per-layer
  // mNROFsTF equality check is unnecessary once this passes: mNROFsTF is a
  // pure function of mROFLength and the single shared nOrbitsPerTF above, so
  // a uniform mROFLength across all layers already guarantees a uniform
  // mNROFsTF.
  const auto uniformTiming = deriveUniformROFTimingConfig(layerTimings);
  if (!uniformTiming.uniform) {
    throw TimeFrameLoadException{
      TimeFrameLoadFailureReason::NonUniformROFTiming,
      std::format("{} CA per-layer ROF timing configuration is not uniform; source-level normalized loading requires equal length/delay/bias/addTimeErr across all {} layers",
                  detName<DetId>(), NLayers)};
  }

  ROFOverlapTableN rofTable;
  ROFVertexLookupTableN vtxTable;
  for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
    rofTable.defineLayer(iLayer, layerTimings[iLayer]);
    vtxTable.defineLayer(iLayer, layerTimings[iLayer]);
  }
  rofTable.init();
  vtxTable.init();
  mROFOverlapTable = std::move(rofTable);
  mROFVertexLookupTable = std::move(vtxTable);
  getScratch().setROFViews(RuntimeROFViews{mROFOverlapTable.getView(), mROFVertexLookupTable.getView(), {}, {}});

  return uniformTiming.config;
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::configureROFMask(gsl::span<const o2::itsmft::ROFRecord> rofs,
                                                        gsl::span<const o2::dataformats::IRFrame> irFrames)
{
  ROFMaskTableN mask{mROFOverlapTable};
  mask.resetMask();

  if constexpr (DetId == o2::detectors::DetID::MFT) {
    const auto& trackingParam = o2::mft::MFTTrackingParam::Instance();
    const bool useIrFilter = trackingParam.irFramesOnly && !irFrames.empty();

    if (trackingParam.isMultCutRequested()) {
      LOGP(info, "{} CA multiplicity filter: Min nClusters = {} ; Max nClusters = {}",
           detName<DetId>(), trackingParam.cutMultClusLow, trackingParam.cutMultClusHigh);
    }
    if (useIrFilter) {
      LOGP(info, "{} CA IRFrame filter enabled with {} ITS IR frames", detName<DetId>(), irFrames.size());
    }

    const auto nROFs = getScratch().getROFOverlapView().getLayer(0).mNROFsTF;
    for (int iRof = 0; iRof < static_cast<int>(nROFs); ++iRof) {
      bool accept = true;
      if (iRof < static_cast<int>(rofs.size())) {
        if (useIrFilter) {
          accept = rofOverlapsIRFrames(rofs[iRof], mMFTROFrameLengthInBC, irFrames);
        }
        if (accept && trackingParam.isMultCutRequested()) {
          accept = trackingParam.isPassingMultCut(rofs[iRof].getNEntries());
        }
      }
      if (accept) {
        for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
          mask.setROFEnabled(iLayer, iRof, 1);
        }
      }
    }
  } else {
    for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
      mask.setROFsEnabled(iLayer, 0, getScratch().getROFOverlapView().getLayer(iLayer).mNROFsTF, 1);
    }
  }

  mMultiplicityMask = std::move(mask);
  getScratch().setROFViews(RuntimeROFViews{mROFOverlapTable.getView(), mROFVertexLookupTable.getView(), mMultiplicityMask.getView(), mUPCMask.getView()});
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::validateROFInput(gsl::span<const o2::itsmft::ROFRecord> rofs) const
{
  const auto expectedROFsTF = getScratch().getROFOverlapView().getLayer(0).mNROFsTF;
  if (rofs.size() != expectedROFsTF) {
    LOGP(warn, "{} CA ROF count differs from continuous timing expectation: received {} expected {}",
         detName<DetId>(), rofs.size(), expectedROFsTF);
  }
}

template class ITSMFTTrackingInterface<7>;
template class ITSMFTTrackingInterface<10>;

} // namespace o2::itsmft::tracking
