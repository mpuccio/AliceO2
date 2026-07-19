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

#include "ITSMFTTracking/DetectorTraits.h"
#include "ITSMFTTracking/TrackingInterface.h"

#include <limits>

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
                                                          std::unique_ptr<DetectorSurfaceCatalogProvider> catalogProvider,
                                                          std::unique_ptr<ClusterDecoder> clusterDecoder)
  : mUseMC(useMC), mOverrideBeamEstimation(overrideBeamEst), mTrackingMode(mode), mDetectorSurfaceCatalogProvider(std::move(catalogProvider)),
    mClusterDecoder(clusterDecoder != nullptr ? std::move(clusterDecoder) : std::make_unique<GeometryClusterDecoder<DetId>>())
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
  resolveTrackingParameters();
  initialiseMemoryPool();
  initialiseTracker();
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::resolveTrackingParameters()
{
  const auto& tc = o2::itsmft::tracking::TrackerParamRef<DetId>::get();
  auto mode = mTrackingMode;
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
  mTrackingMode = mode;
  mTrackParams = o2::itsmft::TrackingMode::getTrackingParameters(DetId, mode);
  LOGP(info, "{} CA tracker initialized in {} mode with {} iteration(s)",
       detName<DetId>(), o2::itsmft::TrackingMode::toString(mode), mTrackParams.size());
  for (size_t i = 0; i < mTrackParams.size(); ++i) {
    LOGP(info, "  iteration {}: {}", i, mTrackParams[i].asString());
  }
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::initialiseMemoryPool()
{
  size_t maxMemory = std::numeric_limits<size_t>::max();
  if (!mTrackParams.empty() && mTrackParams[0].MaxMemory != maxMemory) {
    maxMemory = mTrackParams[0].MaxMemory;
  }
  mMemoryPool = std::make_shared<BoundedMemoryResourceN>(maxMemory);
  mTimeFrame.setMemoryPool(mMemoryPool);
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::initialiseTracker()
{
  if (mTrackParams.empty()) {
    return;
  }
  const auto& tc = o2::itsmft::tracking::TrackerParamRef<DetId>::get();
  mTrackerTraits = std::make_unique<TrackerTraitsN>();
  mTrackerTraits->setMemoryPool(mMemoryPool);
  std::shared_ptr<tbb::task_arena> taskArena;
  mTrackerTraits->setNThreads(tc.nThreads, taskArena);
  mTracker = std::make_unique<TrackerN>(mTrackerTraits.get());
}

template <int NLayers>
float ITSMFTTrackingInterface<NLayers>::processTimeFrame(gsl::span<const o2::itsmft::ROFRecord> rofs,
                                                         gsl::span<const o2::itsmft::CompClusterExt> clusters,
                                                         gsl::span<const unsigned char> patterns,
                                                         const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels,
                                                         gsl::span<const o2::dataformats::IRFrame> irFrames)
{
  if (mTrackParams.empty()) {
    LOGP(info, "{} CA tracking mode is off, skipping TimeFrame processing", detName<DetId>());
    return 0.f;
  }
  if (mDict == nullptr) {
    LOGP(fatal, "{} CA tracker cluster dictionary is not available", detName<DetId>());
  }

  loadTimeFrame(rofs, clusters, patterns, labels, irFrames);
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
  // mTimeFrame if per-layer DPLAlpideParam values disagree; otherwise
  // configures mTimeFrame's ROF overlap/vertex-lookup tables as a side
  // effect and returns the single source-level ROFTimingConfig loadNormalizedSource() needs below.
  const ROFTimingConfig timing = configureROFLookupTables();
  validateROFInput(rofs);

  mTimeFrame.setBz(o2::base::Propagator::Instance()->getNominalBz());
  configureBeamPosition();

  configureROFMask(rofs, irFrames);

  // origin: the TF-relative BC anchor every source ROF is converted against
  // (Architecture.md 7.2). The first ROF's own real interaction record, or
  // the explicit default InteractionRecord{} when there are no ROFs at all.
  const o2::InteractionRecord origin = rofs.empty() ? o2::InteractionRecord{} : rofs.front().getBCData();

  // Transactional (TimeFrame::loadNormalizedSource()/loadSources()): on any
  // failure below, mTimeFrame's normalized frame and every legacy
  // compatibility container are left exactly as they were before this call.
  const auto result = mTimeFrame.loadNormalizedSource(*mClusterDecoder, origin, timing, clusters, patterns, rofs, mDict, labels, DetId);
  if (!result.ok()) {
    if (isRecoverableLoadError(result.error, result.timingDetail)) {
      throw RecoverableLoadFailure{result};
    }
    throw TimeFrameLoadException{result};
  }

  configureTrackingTopology();

  LOGP(info, "{} CA loaded {} clusters from {} ROFs into TimeFrame ({} pattern bytes, MC={})",
       detName<DetId>(), mTimeFrame.getTotalClusters(), rofs.size(), patterns.size(), labels != nullptr);

  for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
    LOGP(info, "  layer {}: {} ROF slots", iLayer, mTimeFrame.getNrof(iLayer));
  }
}

template <int NLayers>
float ITSMFTTrackingInterface<NLayers>::runTracking()
{
  if (!mTracker || mTrackParams.empty()) {
    return 0.f;
  }
  mTracker->adoptTimeFrame(mTimeFrame);
  mTracker->setParameters(mTrackParams);
  mTracker->setMemoryPool(mMemoryPool);
  mTracker->setBz(mTimeFrame.getBz());
  const float elapsedMs = mTracker->clustersToTracks();
  if (elapsedMs < 0.f) {
    LOGP(warn, "{} CA tracking failed for this TF", detName<DetId>());
    return elapsedMs;
  }
  LOGP(info, "{} CA tracking produced {} tracks in {:.2f} ms",
       detName<DetId>(), mTimeFrame.getNumberOfTracks(), elapsedMs);
  return elapsedMs;
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::configureTrackingTopology()
{
  if (mTrackParams.empty()) {
    return;
  }
  mTimeFrame.initDefaultTrackingTopology(mTrackParams[0], NLayers);
  mTimeFrame.initTrackerTopologies(mTrackParams);
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::configureBeamPosition()
{
  if (mTrackParams.empty()) {
    return;
  }
  const auto& p = mTrackParams[0];
  TrackingLoadPolicyN<NLayers>::configureBeamPosition(mTimeFrame, p, mMeanVertex, mOverrideBeamEstimation);
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

  std::array<o2::its::LayerTiming, NLayers> layerTimings{};
  ROFOverlapTableN rofTable;
  ROFVertexLookupTableN vtxTable;
  for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
    const unsigned int nROFsPerOrbit = o2::constants::lhc::LHCMaxBunches / par.getROFLengthInBC(iLayer);
    layerTimings[iLayer] = o2::its::LayerTiming{
      .mNROFsTF = nROFsPerOrbit * static_cast<unsigned int>(nOrbitsPerTF),
      .mROFLength = static_cast<uint32_t>(par.getROFLengthInBC(iLayer)),
      .mROFDelay = static_cast<uint32_t>(par.getROFDelayInBC(iLayer)),
      .mROFBias = static_cast<uint32_t>(par.getROFBiasInBC(iLayer)),
      .mROFAddTimeErr = mTrackParams.empty()
                         ? static_cast<uint32_t>(o2::itsmft::tracking::TrackerParamRef<DetId>::get().addTimeError[iLayer])
                         : mTrackParams[0].AddTimeError[iLayer]};
    rofTable.defineLayer(iLayer, layerTimings[iLayer]);
    vtxTable.defineLayer(iLayer, layerTimings[iLayer]);
  }

  const auto nROFsLayer0 = rofTable.getLayer(0).mNROFsTF;
  for (int iLayer = 1; iLayer < NLayers; ++iLayer) {
    if (rofTable.getLayer(iLayer).mNROFsTF != nROFsLayer0) {
      LOGP(fatal,
           "{} CA single CLUSTERSROF input requires identical mNROFsTF on all {} layers (layer 0: {}, layer {}: {})",
           detName<DetId>(), NLayers, nROFsLayer0, iLayer, rofTable.getLayer(iLayer).mNROFsTF);
    }
  }

  // TimeFrame::loadNormalizedSource() takes one source-level ROFTimingConfig,
  // but DPLAlpideParam supports genuine per-layer staggering (see
  // ROFTimingUniformity.h). Both ITS and MFT default every staggering
  // override to zero, so production configurations are uniform out of the
  // box; a divergent configuration is rejected here, structurally, before
  // any mTimeFrame mutation below -- never silently collapsed.
  const auto uniformTiming = deriveUniformROFTimingConfig(layerTimings);
  if (!uniformTiming.uniform) {
    throw TimeFrameLoadException{
      TimeFrameLoadFailureReason::NonUniformROFTiming,
      std::format("{} CA per-layer ROF timing configuration is not uniform; source-level normalized loading requires equal length/delay/bias/addTimeErr across all {} layers",
                  detName<DetId>(), NLayers)};
  }

  rofTable.init();
  mTimeFrame.setROFOverlapTable(std::move(rofTable));
  vtxTable.init();
  mTimeFrame.setROFVertexLookupTable(std::move(vtxTable));

  return uniformTiming.config;
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::configureROFMask(gsl::span<const o2::itsmft::ROFRecord> rofs,
                                                          gsl::span<const o2::dataformats::IRFrame> irFrames)
{
  ROFMaskTableN mask{mTimeFrame.getROFOverlapTable()};
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

    const auto nROFs = mTimeFrame.getROFOverlapTableView().getLayer(0).mNROFsTF;
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
      mask.setROFsEnabled(iLayer, 0, mTimeFrame.getROFOverlapTableView().getLayer(iLayer).mNROFsTF, 1);
    }
  }

  mTimeFrame.setMultiplicityCutMask(std::move(mask));
}

template <int NLayers>
void ITSMFTTrackingInterface<NLayers>::validateROFInput(gsl::span<const o2::itsmft::ROFRecord> rofs) const
{
  const auto expectedROFsTF = mTimeFrame.getROFOverlapTableView().getLayer(0).mNROFsTF;
  if (rofs.size() != expectedROFsTF) {
    LOGP(warn, "{} CA ROF count differs from continuous timing expectation: received {} expected {}",
         detName<DetId>(), rofs.size(), expectedROFsTF);
  }
}

template class ITSMFTTrackingInterface<7>;
template class ITSMFTTrackingInterface<10>;

} // namespace o2::itsmft::tracking
