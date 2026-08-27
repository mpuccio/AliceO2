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

/// @file   CATrackerSpec.cxx

#include "MFTWorkflow/CATrackerSpec.h"

#include <array>
#include <algorithm>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
#include <stdexcept>

#include <gsl/span>

#include "CommonDataFormat/IRFrame.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DetectorsBase/GeometryManager.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Logger.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/GenericTrackOutputAdapter.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/SurfaceTiming.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "DetectorsBase/Propagator.h"
#include <oneapi/tbb/task_arena.h>
#include "CommonConstants/LHCConstants.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTTracking/Constants.h"
#include "MFTTracking/MFTTrackingParam.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::framework;

namespace o2::mft
{

namespace
{
using namespace o2::itsmft::tracking;

template <int NLayers>
constexpr std::array<LayerId, NLayers> detectorLocalToLayoutLayers()
{
  std::array<LayerId, NLayers> order{};
  for (int i = 0; i < NLayers; ++i) {
    order[i] = LayerId{static_cast<uint16_t>(i)};
  }
  return order;
}

inline constexpr auto kLayerToLayout = detectorLocalToLayoutLayers<MFTNLayers>();

bool rofOverlapsIRFrames(const o2::itsmft::ROFRecord& rof, int rofLengthInBC,
                         gsl::span<const o2::dataformats::IRFrame> irFrames)
{
  o2::InteractionRecord start{rof.getBCData()};
  const o2::InteractionRecord end = start + rofLengthInBC - 1;
  const o2::dataformats::IRFrame reference{start, end};
  for (const auto& ir : irFrames) {
    if (ir.info > 0 && reference.getOverlap(ir).isValid()) {
      return true;
    }
  }
  return false;
}

} // namespace

CATrackerDPL::CATrackerDPL(std::shared_ptr<o2::base::GRPGeomRequest> gr, bool useMC,
                           o2::itsmft::TrackingMode::Type trMode)
  : mGGCCDBRequest(std::move(gr)), mUseMC(useMC), mTrackingMode(trMode)
{
  mClusterDecoder = std::make_unique<o2::itsmft::tracking::MFTGeometryClusterDecoder>();
}

void CATrackerDPL::configureROFViews(gsl::span<const o2::itsmft::ROFRecord> rofs,
                                     gsl::span<const o2::dataformats::IRFrame> irFrames)
{
  const auto& detector = mTracker->getDetectorConfiguration();
  const auto& alpParams = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::MFT>::Instance();
  const bool continuous = o2::base::GRPGeomHelper::instance().getGRPECS()->isDetContinuousReadOut(o2::detectors::DetID::MFT);
  mMFTROFrameLengthInBC = continuous ? alpParams.roFrameLengthInBC : std::max(1, static_cast<int>(alpParams.roFrameLengthTrig / (o2::constants::lhc::LHCBunchSpacingNS * 1e3)));
  const int nOrbitsPerTF = o2::base::GRPGeomHelper::getNHBFPerTF();
  std::array<o2::its::LayerTiming, o2::itsmft::tracking::MFTNLayers> layerTimings{};
  for (int layer = 0; layer < o2::itsmft::tracking::MFTNLayers; ++layer) {
    const auto length = alpParams.getROFLengthInBC(layer);
    if (length <= 0) {
      throw o2::itsmft::tracking::TimeFrameLoadException{
        o2::itsmft::tracking::TimeFrameLoadFailureReason::NonUniformROFTiming,
        "MFT CA per-layer ROF timing has a non-positive ROF length"};
    }
    const auto nROFsPerOrbit = o2::constants::lhc::LHCMaxBunches / static_cast<unsigned int>(length);
    layerTimings[layer] = o2::its::LayerTiming{
      .mNROFsTF = nROFsPerOrbit * static_cast<unsigned int>(nOrbitsPerTF),
      .mROFLength = static_cast<uint32_t>(length),
      .mROFDelay = static_cast<uint32_t>(alpParams.getROFDelayInBC(layer)),
      .mROFBias = static_cast<uint32_t>(alpParams.getROFBiasInBC(layer)),
      .mROFAddTimeErr = detector.addTimeError[layer]};
    if (layerTimings[layer].mNROFsTF == 0) {
      throw o2::itsmft::tracking::TimeFrameLoadException{
        o2::itsmft::tracking::TimeFrameLoadFailureReason::ZeroROFCount,
        "MFT CA per-layer ROF timing yields zero ROFs per TimeFrame"};
    }
  }
  const auto uniform = o2::itsmft::tracking::deriveUniformROFTimingConfig(layerTimings);
  if (!uniform.uniform) {
    throw o2::itsmft::tracking::TimeFrameLoadException{
      o2::itsmft::tracking::TimeFrameLoadFailureReason::NonUniformROFTiming,
      "MFT CA per-layer ROF timing configuration is not uniform"};
  }

  o2::its::ROFOverlapTable<o2::itsmft::tracking::MFTNLayers> overlap;
  o2::its::ROFVertexLookupTable<o2::itsmft::tracking::MFTNLayers> vertex;
  for (int layer = 0; layer < o2::itsmft::tracking::MFTNLayers; ++layer) {
    overlap.defineLayer(layer, layerTimings[layer]);
    vertex.defineLayer(layer, layerTimings[layer]);
  }
  overlap.init();
  vertex.init();
  o2::its::ROFMaskTable<o2::itsmft::tracking::MFTNLayers> mask{overlap};
  mask.resetMask();
  const auto& trackingParam = o2::mft::MFTTrackingParam::Instance();
  const bool useIrFilter = trackingParam.irFramesOnly && !irFrames.empty();
  const auto nROFs = layerTimings[0].mNROFsTF;
  for (int rof = 0; rof < static_cast<int>(nROFs); ++rof) {
    bool accept = rof >= static_cast<int>(rofs.size());
    if (!accept) {
      accept = (!useIrFilter || rofOverlapsIRFrames(rofs[rof], mMFTROFrameLengthInBC, irFrames)) &&
               (!trackingParam.isMultCutRequested() || trackingParam.isPassingMultCut(rofs[rof].getNEntries()));
    }
    if (accept) {
      for (int layer = 0; layer < o2::itsmft::tracking::MFTNLayers; ++layer) {
        mask.setROFEnabled(layer, rof, 1);
      }
    }
  }
  mROFOverlapTable = std::move(overlap);
  mROFVertexLookupTable = std::move(vertex);
  mMultiplicityMask = std::move(mask);
  mFrame.setROFViews({mROFOverlapTable.getView(), mROFVertexLookupTable.getView(), mMultiplicityMask.getView(), mUPCMask.getView()});
}

void CATrackerDPL::invalidatePublication() noexcept
{
  mPublicationClock.reset();
  mExternalIndicesBySurface.clear();
  mClusterSizesBySurface.clear();
  mFrame.setROFViews({});
}

void CATrackerDPL::initialiseTracking()
{
  auto mode = mTrackingMode;
  const auto& trackerParams = o2::itsmft::tracking::TrackerParamRef<o2::detectors::DetID::MFT>::get();
  const auto configuredMode = static_cast<o2::itsmft::TrackingMode::Type>(trackerParams.trackingMode);
  if (mode == o2::itsmft::TrackingMode::Unset ||
      (configuredMode != o2::itsmft::TrackingMode::Unset && mode != configuredMode)) {
    if (configuredMode != o2::itsmft::TrackingMode::Unset) {
      LOGP(info, "MFT CA tracking mode overwritten by configurable params from {} to {}",
           o2::itsmft::TrackingMode::toString(mode), o2::itsmft::TrackingMode::toString(configuredMode));
      mode = configuredMode;
    }
  }
  if (mode == o2::itsmft::TrackingMode::Unset) {
    mode = o2::itsmft::TrackingMode::Sync;
  }
  mTrackingMode = mode;
  const auto parameters = o2::itsmft::TrackingMode::getTrackingParameters(o2::detectors::DetID::MFT, mode);
  LOGP(info, "MFT CA tracker initialized in {} mode with {} iteration(s)",
       o2::itsmft::TrackingMode::toString(mode), parameters.size());
  if (parameters.empty()) {
    return;
  }

  mTrackerTraits = std::make_unique<o2::itsmft::tracking::TrackerTraits>();
  std::shared_ptr<tbb::task_arena> taskArena;
  mTrackerTraits->setNThreads(trackerParams.nThreads, taskArena);

  o2::itsmft::tracking::TrackerInitialization configuration;
  configuration.catalog = o2::itsmft::tracking::SurfaceCatalogView{o2::itsmft::tracking::kMFTStaticSurfaceCatalog.data(),
                                                                   static_cast<uint32_t>(o2::itsmft::tracking::kMFTStaticSurfaceCatalog.size())};
  configuration.layout = o2::itsmft::tracking::makeDetectorLayout(
    o2::itsmft::tracking::LayerMask{trackerParams.holeLayerMask});
  configuration.memoryPool = std::make_shared<o2::itsmft::tracking::BoundedMemoryResource>(parameters.front().MaxMemory);
  configuration.parameters = parameters;

  mTracker = std::make_unique<o2::itsmft::tracking::Tracker>();
  const auto result = mTracker->initialize(mFrame, configuration);
  if (!result.ok()) {
    LOGP(fatal, "MFT CA tracker failed to initialize static configuration (error={} iteration={} layout={})",
         static_cast<int>(result.error), result.failedIteration, static_cast<int>(result.layoutError));
  }
}

o2::itsmft::tracking::TrackingOutcome CATrackerDPL::processTimeFrame(
  gsl::span<const o2::itsmft::ROFRecord> rofs,
  gsl::span<const o2::itsmft::CompClusterExt> clusters,
  gsl::span<const unsigned char> patterns,
  const o2::dataformats::MCTruthContainer<MCCompLabel>* labels,
  gsl::span<const o2::dataformats::IRFrame> irFrames)
{
  if (!isActive()) {
    LOGP(info, "MFT CA tracking mode is off, skipping TimeFrame processing");
    return o2::itsmft::tracking::TrackingOutcome::Success;
  }
  const auto& params = mTracker->getIterationConfigurations().front().parameters;
  mFrame.setBz(o2::base::Propagator::Instance()->getNominalBz());
  const auto views = mFrame.getROFViews();
  if (views.overlap.mLayerCount > 0 && rofs.size() != views.overlap.getLayer(0).mNROFsTF) {
    LOGP(warn, "MFT CA ROF count differs from continuous timing expectation: received {} expected {}",
         rofs.size(), views.overlap.getLayer(0).mNROFsTF);
  }
  try {
    if (mDictionary == nullptr) {
      throw o2::itsmft::tracking::TimeFrameLoadException{
        o2::itsmft::tracking::TimeFrameLoadFailureReason::DictionaryNotConfigured,
        "MFT CA tracker cluster dictionary is not available"};
    }
    const auto rofViews = mFrame.getROFViews();
    if (rofViews.overlap.mLayerCount <= 0) {
      throw o2::itsmft::tracking::TimeFrameLoadException{
        o2::itsmft::tracking::TimeFrameLoadFailureReason::NonUniformROFTiming,
        "MFT CA tracker received no adapter-owned runtime ROF timing view"};
    }
    const auto& clock = rofViews.overlap.getLayer(0);
    o2::itsmft::tracking::ClusterSourceInput source;
    source.id = o2::itsmft::tracking::ClusterSourceId{0};
    source.detector = o2::detectors::DetID::MFT;
    source.clusters = clusters;
    source.patterns = patterns;
    source.rofs = rofs;
    source.dictionary = mDictionary;
    source.labels = labels;
    source.layerToSurface = kLayerToLayout;
    source.timing = o2::itsmft::tracking::ROFTimingConfig{clock.mROFLength, clock.mROFDelay, clock.mROFBias, clock.mROFAddTimeErr};
    source.decoder = mClusterDecoder.get();
    source.rofViews = rofViews;
    const auto origin = rofs.empty() ? o2::InteractionRecord{} : rofs.front().getBCData();
    const auto loaded = o2::itsmft::tracking::loadTimeFrameSources(
      mFrame, gsl::span<const o2::itsmft::tracking::ClusterSourceInput>{&source, 1}, mFrame.getLayout().getSurfaceCatalog(), origin,
      &mExternalIndicesBySurface, &mClusterSizesBySurface);
    if (!loaded.ok()) {
      if (o2::itsmft::tracking::isRecoverableLoadError(loaded.error, loaded.timingDetail)) {
        throw o2::itsmft::tracking::RecoverableLoadFailure{loaded};
      }
      throw o2::itsmft::tracking::TimeFrameLoadException{loaded};
    }
  } catch (const o2::itsmft::tracking::RecoverableLoadFailure& failure) {
    LOGP(error, "MFT CA loading recoverably failed: {}", failure.what());
    resetTimeFrame();
    if (params.DropTFUponFailure) {
      return o2::itsmft::tracking::TrackingOutcome::RecoverableDropped;
    }
    throw;
  } catch (const o2::itsmft::tracking::BoundedMemoryResource::MemoryLimitExceeded& failure) {
    LOGP(error, "MFT CA loading exceeded memory limit: {}", failure.what());
    resetTimeFrame();
    if (params.DropTFUponFailure) {
      return o2::itsmft::tracking::TrackingOutcome::RecoverableDropped;
    }
    throw;
  } catch (const std::bad_alloc& failure) {
    LOGP(error, "MFT CA loading allocation failed: {}", failure.what());
    resetTimeFrame();
    if (params.DropTFUponFailure) {
      return o2::itsmft::tracking::TrackingOutcome::RecoverableDropped;
    }
    throw;
  } catch (const o2::itsmft::tracking::TimeFrameLoadException& failure) {
    LOGP(error, "MFT CA loading hit a structural failure: {}", failure.what());
    resetTimeFrame();
    throw;
  } catch (const std::exception& failure) {
    LOGP(error, "MFT CA loading failed with an unclassified exception: {}", failure.what());
    resetTimeFrame();
    throw;
  }

  const auto result = mTracker->run(mFrame, *mTrackerTraits);
  if (result.outcome == o2::itsmft::tracking::TrackingOutcome::RecoverableDropped) {
    LOGP(warn, "MFT CA tracking failed for this TF");
  } else {
    LOGP(info, "MFT CA tracking produced {} tracks in {:.2f} ms", mFrame.getGenericTracks().size(), result.elapsedMs);
  }
  return result.outcome;
}

void CATrackerDPL::resetTimeFrame() noexcept
{
  mExternalIndicesBySurface.clear();
  mClusterSizesBySurface.clear();
  mFrame.resetTimeFrame();
}

CATrackerPublicationAction decideCATrackerPublicationAction(bool trackerActive, o2::itsmft::tracking::TrackingOutcome outcome) noexcept
{
  if (!trackerActive) {
    return CATrackerPublicationAction::PublishInactiveEmpty;
  }
  if (outcome == o2::itsmft::tracking::TrackingOutcome::RecoverableDropped) {
    return CATrackerPublicationAction::SkipDroppedTimeFrame;
  }
  return CATrackerPublicationAction::PublishActiveResult;
}

void CATrackerDPL::init(InitContext&)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
}

void CATrackerDPL::run(ProcessingContext& pc)
{
  updateTimeDependentParams(pc);

  auto rofsinput = pc.inputs().get<const std::vector<o2::itsmft::ROFRecord>>("ROframes");

  if (decideCATrackerPublicationAction(isActive(), o2::itsmft::tracking::TrackingOutcome::Success) == CATrackerPublicationAction::PublishInactiveEmpty) {
    // Existing production behavior, preserved exactly: publish the input
    // ROFs verbatim (their firstEntry/nEntries are not rewritten here) plus
    // empty track/cluster-index/seed-pattern outputs, when the tracker is
    // not configured to run.
    pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{"MFT", "MFTTrackROF", 0},
                                                          rofsinput.begin(), rofsinput.end());
    pc.outputs().make<std::vector<o2::mft::TrackMFT>>(Output{"MFT", "TRACKS", 0});
    pc.outputs().make<std::vector<int>>(Output{"MFT", "TRACKCLSID", 0});
    pc.outputs().make<std::vector<uint16_t>>(Output{"MFT", "TRACKSEEDPAT", 0});
    return;
  }

  auto compClusters = pc.inputs().get<const std::vector<o2::itsmft::CompClusterExt>>("compClusters");
  gsl::span<const unsigned char> patterns = pc.inputs().get<gsl::span<unsigned char>>("patterns");

  const dataformats::MCTruthContainer<MCCompLabel>* labels = nullptr;
  if (mUseMC && pc.inputs().getPos("labels") >= 0) {
    labels = pc.inputs().get<const dataformats::MCTruthContainer<MCCompLabel>*>("labels").release();
  }

  gsl::span<const o2::dataformats::IRFrame> irFrames;
  if (pc.inputs().getPos("IRFramesITS") >= 0) {
    irFrames = pc.inputs().get<gsl::span<o2::dataformats::IRFrame>>("IRFramesITS");
  }

  LOGP(info, "MFT CA input pulled {} compressed clusters in {} RO frames ({} pattern bytes)",
       compClusters.size(), rofsinput.size(), patterns.size());

  try {
    configureROFViews(gsl::span<const o2::itsmft::ROFRecord>(rofsinput.data(), rofsinput.size()), irFrames);
    const auto trackingResult = processTimeFrame(gsl::span<const o2::itsmft::ROFRecord>(rofsinput.data(), rofsinput.size()),
                                                 gsl::span<const o2::itsmft::CompClusterExt>(compClusters.data(), compClusters.size()),
                                                 patterns, labels, irFrames);

    if (decideCATrackerPublicationAction(isActive(), trackingResult) == CATrackerPublicationAction::SkipDroppedTimeFrame) {
      LOGP(error, "MFT CA tracking dropped this TimeFrame ({} ROFs, {} clusters); publishing nothing and continuing with the next TimeFrame",
           rofsinput.size(), compClusters.size());
      invalidatePublication();
      return;
    }

    {
      mPublicationClock.emplace(mROFOverlapTable.getView().getClockLayer());
      const o2::itsmft::tracking::GenericTrackPublicationContext context{
        o2::detectors::DetID::MFT, o2::itsmft::tracking::ClusterSourceId{0},
        gsl::span<const o2::itsmft::ROFRecord>{rofsinput.data(), rofsinput.size()}, *mPublicationClock,
        kLayerToLayout,
        &mExternalIndicesBySurface, &mClusterSizesBySurface};
      o2::itsmft::tracking::GenericTrackOutputAdapterError error = o2::itsmft::tracking::GenericTrackOutputAdapterError::None;
      const auto staged = o2::itsmft::tracking::stageMFTGenericTrackOutput(mFrame, context, mUseMC, error);
      if (!staged) {
        throw std::runtime_error{"MFT GenericTrack output staging failed"};
      }

      auto& trackROFs = pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{"MFT", "MFTTrackROF", 0},
                                                                              staged->trackROFs.begin(), staged->trackROFs.end());
      auto& allTracksMFT = pc.outputs().make<std::vector<o2::mft::TrackMFT>>(Output{"MFT", "TRACKS", 0});
      allTracksMFT.assign(staged->tracks.begin(), staged->tracks.end());
      auto& allClusIdx = pc.outputs().make<std::vector<int>>(Output{"MFT", "TRACKCLSID", 0});
      allClusIdx.assign(staged->clusterIndices.begin(), staged->clusterIndices.end());
      auto& allSeedPatterns = pc.outputs().make<std::vector<uint16_t>>(Output{"MFT", "TRACKSEEDPAT", 0});
      allSeedPatterns.assign(staged->seedPatterns.begin(), staged->seedPatterns.end());
      LOGP(info, "MFT CA pushed {} tracks in {} ROFs", allTracksMFT.size(), trackROFs.size());
      if (mUseMC) {
        pc.outputs().snapshot(Output{"MFT", "TRACKSMCTR", 0}, staged->labels);
        LOGP(info, "MFT CA pushed {} track MC labels", staged->labels.size());
      }
    }
  } catch (...) {
    resetTimeFrame();
    invalidatePublication();
    throw;
  }

  resetTimeFrame();
  invalidatePublication();
}

void CATrackerDPL::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  if (!mTrackingInitialised) {
    mTrackingInitialised = true;
    initialiseTracking();
  }
  static bool initOnceDone = false;
  if (!initOnceDone) {
    initOnceDone = true;
    if (pc.inputs().getPos("mftTGeo") >= 0) {
      pc.inputs().get<o2::mft::GeometryTGeo*>("mftTGeo");
    }
    pc.inputs().get<o2::itsmft::TopologyDictionary*>("cldict");
    o2::mft::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                                                o2::math_utils::TransformType::T2GRot,
                                                                                o2::math_utils::TransformType::T2G,
                                                                                o2::math_utils::TransformType::L2G));
  }
}

void CATrackerDPL::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("MFT", "CLUSDICT", 0)) {
    LOG(info) << "MFT CA input cluster dictionary updated";
    mDictionary = static_cast<const o2::itsmft::TopologyDictionary*>(obj);
    return;
  }
  if (matcher == ConcreteDataMatcher("MFT", "GEOMTGEO", 0)) {
    LOG(info) << "MFT CA input GeometryTGeo loaded from CCDB";
    o2::mft::GeometryTGeo::adopt(static_cast<o2::mft::GeometryTGeo*>(obj));
    o2::mft::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                                                o2::math_utils::TransformType::T2GRot,
                                                                                o2::math_utils::TransformType::T2G,
                                                                                o2::math_utils::TransformType::L2G));
    // The catalog has static process lifetime; geometry adoption remains
    // necessary for raw cluster decoding.
    return;
  }
}

DataProcessorSpec getCATrackerSpec(bool useMC, bool useGeom, bool useIRFrames, o2::itsmft::TrackingMode::Type trMode)
{
  std::vector<InputSpec> inputs;
  inputs.emplace_back("compClusters", "MFT", "COMPCLUSTERS", 0, Lifetime::Timeframe);
  inputs.emplace_back("patterns", "MFT", "PATTERNS", 0, Lifetime::Timeframe);
  inputs.emplace_back("ROframes", "MFT", "CLUSTERSROF", 0, Lifetime::Timeframe);
  inputs.emplace_back("cldict", "MFT", "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec("MFT/Calib/ClusterDictionary"));

  if (useMC) {
    inputs.emplace_back("labels", "MFT", "CLUSTERSMCTR", 0, Lifetime::Timeframe);
  }

  const auto& trackingParam = MFTTrackingParam::Instance();
  if (useIRFrames || trackingParam.irFramesOnly) {
    inputs.emplace_back("IRFramesITS", "ITS", "IRFRAMES", 0, Lifetime::Timeframe);
  }

  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,
                                                              true,
                                                              false,
                                                              true,
                                                              true,
                                                              useGeom ? o2::base::GRPGeomRequest::Aligned : o2::base::GRPGeomRequest::None,
                                                              inputs,
                                                              true);
  if (!useGeom) {
    ggRequest->addInput({"mftTGeo", "MFT", "GEOMTGEO", 0, Lifetime::Condition, framework::ccdbParamSpec("MFT/Config/Geometry")}, inputs);
  }

  std::vector<OutputSpec> outputs;
  outputs.emplace_back("MFT", "TRACKS", 0, Lifetime::Timeframe);
  outputs.emplace_back("MFT", "MFTTrackROF", 0, Lifetime::Timeframe);
  outputs.emplace_back("MFT", "TRACKCLSID", 0, Lifetime::Timeframe);
  outputs.emplace_back("MFT", "TRACKSEEDPAT", 0, Lifetime::Timeframe);
  if (useMC) {
    outputs.emplace_back("MFT", "TRACKSMCTR", 0, Lifetime::Timeframe);
  }

  return DataProcessorSpec{
    "mft-ca-tracker",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<CATrackerDPL>(ggRequest, useMC, trMode)},
    Options{}};
}

} // namespace o2::mft
