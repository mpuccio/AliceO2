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

#include "ITSCAWorkflow/CATrackerSpec.h"

#include <stdexcept>
#include <array>
#include <limits>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

#include <gsl/span>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsBase/GeometryManager.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Logger.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/GenericTrackOutputAdapter.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/SurfaceTiming.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITSMFTTracking/BoundedAllocator.h"
#include "CommonConstants/LHCConstants.h"
#include "DetectorsBase/Propagator.h"
#include <oneapi/tbb/task_arena.h>
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::framework;

namespace o2::its::ca
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

inline constexpr auto kLayerToLayout = detectorLocalToLayoutLayers<ITSNLayers>();

template <int NLayers>
bool completePublication(DetectorPublicationAdapter<NLayers>& publication,
                         const TimeFrame& frame,
                         const Tracker& tracker,
                         ClusterSourceId source,
                         const TrackingResult& result)
{
  const auto configurations = tracker.getIterationConfigurations();
  std::size_t firstTrack = 0;
  for (std::size_t iteration = 0; iteration < configurations.size(); ++iteration) {
    if (iteration >= result.acceptedTrackCounts.size() ||
        result.acceptedTrackCounts[iteration] > frame.getGenericTracks().size() - firstTrack) {
      return false;
    }
    std::vector<uint32_t> trackIndices(result.acceptedTrackCounts[iteration]);
    std::iota(trackIndices.begin(), trackIndices.end(), static_cast<uint32_t>(firstTrack));
    if (!publication.completeAccepted(trackIndices, configurations[iteration].parameters, frame, iteration + 1 == configurations.size())) {
      return false;
    }
    firstTrack += result.acceptedTrackCounts[iteration];
  }
  return firstTrack == frame.getGenericTracks().size();
}

} // namespace

CATrackerDPL::CATrackerDPL(std::shared_ptr<o2::base::GRPGeomRequest> gr, bool useMC,
                           o2::itsmft::TrackingMode::Type trMode)
  : mGGCCDBRequest(std::move(gr)), mUseMC(useMC), mTrackingMode(trMode)
{
  mClusterDecoder = std::make_unique<o2::itsmft::tracking::ITSGeometryClusterDecoder>();
  mPublicationAdapter.adoptITSSharedClusterCompatibility(&mCompatibility);
}

void CATrackerDPL::configureROFViews(gsl::span<const o2::itsmft::ROFRecord> rofs)
{
  const auto& detector = mTracker->getDetectorConfiguration();
  const auto& alpParams = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance();
  const int nOrbitsPerTF = o2::base::GRPGeomHelper::getNHBFPerTF();
  std::array<o2::its::LayerTiming, o2::itsmft::tracking::ITSNLayers> layerTimings{};
  for (int layer = 0; layer < o2::itsmft::tracking::ITSNLayers; ++layer) {
    const auto length = alpParams.getROFLengthInBC(layer);
    if (length <= 0) {
      throw o2::itsmft::tracking::TimeFrameLoadException{
        o2::itsmft::tracking::TimeFrameLoadFailureReason::NonUniformROFTiming,
        "ITS CA per-layer ROF timing has a non-positive ROF length"};
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
        "ITS CA per-layer ROF timing yields zero ROFs per TimeFrame"};
    }
  }
  const auto uniform = o2::itsmft::tracking::deriveUniformROFTimingConfig(layerTimings);
  if (!uniform.uniform) {
    throw o2::itsmft::tracking::TimeFrameLoadException{
      o2::itsmft::tracking::TimeFrameLoadFailureReason::NonUniformROFTiming,
      "ITS CA per-layer ROF timing configuration is not uniform"};
  }

  o2::its::ROFOverlapTable<o2::itsmft::tracking::ITSNLayers> overlap;
  o2::its::ROFVertexLookupTable<o2::itsmft::tracking::ITSNLayers> vertex;
  for (int layer = 0; layer < o2::itsmft::tracking::ITSNLayers; ++layer) {
    overlap.defineLayer(layer, layerTimings[layer]);
    vertex.defineLayer(layer, layerTimings[layer]);
  }
  overlap.init();
  vertex.init();
  o2::its::ROFMaskTable<o2::itsmft::tracking::ITSNLayers> mask{overlap};
  mask.resetMask();
  for (int layer = 0; layer < o2::itsmft::tracking::ITSNLayers; ++layer) {
    mask.setROFsEnabled(layer, 0, static_cast<int>(layerTimings[layer].mNROFsTF), 1);
  }
  mROFOverlapTable = std::move(overlap);
  mROFVertexLookupTable = std::move(vertex);
  mMultiplicityMask = std::move(mask);
  mPublicationAdapter.reset();
  mFrame.setROFViews({mROFOverlapTable.getView(), mROFVertexLookupTable.getView(), mMultiplicityMask.getView(), mUPCMask.getView()});
  (void)rofs;
}

void CATrackerDPL::invalidatePublication() noexcept
{
  mPublicationAdapter.reset();
  mPublicationClock.reset();
  mExternalIndicesBySurface.clear();
  mClusterSizesBySurface.clear();
  mFrame.setROFViews({});
}

void CATrackerDPL::initialiseTracking()
{
  auto mode = mTrackingMode;
  if (mode == o2::itsmft::TrackingMode::Unset) {
    mode = o2::itsmft::TrackingMode::Sync;
  }
  mTrackingMode = mode;
  const auto parameters = o2::itsmft::TrackingMode::getTrackingParameters(o2::detectors::DetID::ITS, mode);
  LOGP(info, "ITS CA tracker initialized in {} mode with {} iteration(s)",
       o2::itsmft::TrackingMode::toString(mode), parameters.size());
  if (parameters.empty()) {
    return;
  }

  mTrackerTraits = std::make_unique<o2::itsmft::tracking::TrackerTraits>();
  std::shared_ptr<tbb::task_arena> taskArena;
  const auto& commonParams = o2::itsmft::ITSCommonCATrackerParam::Instance();
  if (commonParams.nThreads <= 0) {
    LOGP(fatal, "ITS CA tracker requires ITSCommonCATrackerParam.nThreads > 0, got {}", commonParams.nThreads);
  }
  mTrackerTraits->setNThreads(commonParams.nThreads, taskArena);

  o2::itsmft::tracking::TrackerInitialization configuration;
  configuration.catalog = o2::itsmft::tracking::SurfaceCatalogView{o2::itsmft::tracking::kITSStaticSurfaceCatalog.data(),
                                                                   static_cast<uint32_t>(o2::itsmft::tracking::kITSStaticSurfaceCatalog.size())};
  configuration.layout = o2::itsmft::tracking::makeDetectorLayout(
    o2::itsmft::tracking::LayerMask{commonParams.holeLayerMask});
  const auto maxMemory = parameters.front().MaxMemory;
  configuration.memoryPool = std::make_shared<o2::itsmft::tracking::BoundedMemoryResource>(maxMemory);
  configuration.parameters = parameters;

  mTracker = std::make_unique<o2::itsmft::tracking::Tracker>();
  const auto result = mTracker->initialize(mFrame, configuration);
  if (!result.ok()) {
    LOGP(fatal, "ITS CA tracker failed to initialize static configuration (error={} iteration={} layout={})",
         static_cast<int>(result.error), result.failedIteration, static_cast<int>(result.layoutError));
  }
}

o2::itsmft::tracking::TrackingOutcome CATrackerDPL::processTimeFrame(
  gsl::span<const o2::itsmft::ROFRecord> rofs,
  gsl::span<const o2::itsmft::CompClusterExt> clusters,
  gsl::span<const unsigned char> patterns,
  const o2::dataformats::MCTruthContainer<MCCompLabel>* labels)
{
  if (!isActive()) {
    LOGP(info, "ITS CA tracking mode is off, skipping TimeFrame processing");
    return o2::itsmft::tracking::TrackingOutcome::Success;
  }
  const auto& params = mTracker->getIterationConfigurations().front().parameters;
  mFrame.setBz(o2::base::Propagator::Instance()->getNominalBz());
  const auto views = mFrame.getROFViews();
  if (views.overlap.mLayerCount > 0 && rofs.size() != views.overlap.getLayer(0).mNROFsTF) {
    LOGP(warn, "ITS CA ROF count differs from continuous timing expectation: received {} expected {}",
         rofs.size(), views.overlap.getLayer(0).mNROFsTF);
  }
  try {
    if (mDictionary == nullptr) {
      throw o2::itsmft::tracking::TimeFrameLoadException{
        o2::itsmft::tracking::TimeFrameLoadFailureReason::DictionaryNotConfigured,
        "ITS CA tracker cluster dictionary is not available"};
    }
    const auto rofViews = mFrame.getROFViews();
    if (rofViews.overlap.mLayerCount <= 0) {
      throw o2::itsmft::tracking::TimeFrameLoadException{
        o2::itsmft::tracking::TimeFrameLoadFailureReason::NonUniformROFTiming,
        "ITS CA tracker received no adapter-owned runtime ROF timing view"};
    }
    const auto& clock = rofViews.overlap.getLayer(0);
    o2::itsmft::tracking::ClusterSourceInput source;
    source.id = o2::itsmft::tracking::ClusterSourceId{0};
    source.detector = o2::detectors::DetID::ITS;
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
  } catch (const o2::itsmft::tracking::RecoverableLoadFailure& error) {
    LOGP(error, "ITS CA loading recoverably failed: {}", error.what());
    resetTimeFrame();
    if (params.DropTFUponFailure) {
      return o2::itsmft::tracking::TrackingOutcome::RecoverableDropped;
    }
    throw;
  } catch (const o2::itsmft::tracking::BoundedMemoryResource::MemoryLimitExceeded& error) {
    LOGP(error, "ITS CA loading exceeded memory limit: {}", error.what());
    resetTimeFrame();
    if (params.DropTFUponFailure) {
      return o2::itsmft::tracking::TrackingOutcome::RecoverableDropped;
    }
    throw;
  } catch (const std::bad_alloc& error) {
    LOGP(error, "ITS CA loading allocation failed: {}", error.what());
    resetTimeFrame();
    if (params.DropTFUponFailure) {
      return o2::itsmft::tracking::TrackingOutcome::RecoverableDropped;
    }
    throw;
  } catch (const o2::itsmft::tracking::TimeFrameLoadException& error) {
    LOGP(error, "ITS CA loading hit a structural failure: {}", error.what());
    resetTimeFrame();
    throw;
  } catch (const std::exception& error) {
    LOGP(error, "ITS CA loading failed with an unclassified exception: {}", error.what());
    resetTimeFrame();
    throw;
  }

  o2::itsmft::tracking::TrackingResult result;
  try {
    result = mTracker->run(mFrame, *mTrackerTraits);
    if (result.outcome == o2::itsmft::tracking::TrackingOutcome::RecoverableDropped) {
      mPublicationAdapter.reset();
    } else if (!completePublication(mPublicationAdapter, mFrame, *mTracker, o2::itsmft::tracking::ClusterSourceId{0}, result)) {
      mPublicationAdapter.reset();
      throw std::runtime_error{"failed to seal ITS tracking compatibility"};
    }
  } catch (...) {
    mPublicationAdapter.reset();
    throw;
  }
  if (result.outcome == o2::itsmft::tracking::TrackingOutcome::RecoverableDropped) {
    LOGP(warn, "ITS CA tracking failed for this TF");
  } else {
    LOGP(info, "ITS CA tracking produced {} tracks in {:.2f} ms", mFrame.getGenericTracks().size(), result.elapsedMs);
  }
  return result.outcome;
}

void CATrackerDPL::resetTimeFrame() noexcept
{
  mPublicationAdapter.reset();
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
    pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{"ITS", "ITSTrackROF", 0},
                                                          rofsinput.begin(), rofsinput.end());
    pc.outputs().make<std::vector<o2::its::TrackITS>>(Output{"ITS", "TRACKS", 0});
    pc.outputs().make<std::vector<int>>(Output{"ITS", "TRACKCLSID", 0});
    return;
  }

  auto compClusters = pc.inputs().get<const std::vector<o2::itsmft::CompClusterExt>>("compClusters");
  gsl::span<const unsigned char> patterns = pc.inputs().get<gsl::span<unsigned char>>("patterns");

  const dataformats::MCTruthContainer<MCCompLabel>* labels = nullptr;
  if (mUseMC && pc.inputs().getPos("labels") >= 0) {
    labels = pc.inputs().get<const dataformats::MCTruthContainer<MCCompLabel>*>("labels").release();
  }

  LOGP(info, "ITS CA input pulled {} compressed clusters in {} RO frames ({} pattern bytes)",
       compClusters.size(), rofsinput.size(), patterns.size());

  try {
    configureROFViews(gsl::span<const o2::itsmft::ROFRecord>(rofsinput.data(), rofsinput.size()));
    const auto trackingResult = processTimeFrame(gsl::span<const o2::itsmft::ROFRecord>(rofsinput.data(), rofsinput.size()),
                                                 gsl::span<const o2::itsmft::CompClusterExt>(compClusters.data(), compClusters.size()),
                                                 patterns, labels);

    if (decideCATrackerPublicationAction(isActive(), trackingResult) == CATrackerPublicationAction::SkipDroppedTimeFrame) {
      LOGP(error, "ITS CA tracking dropped this TimeFrame ({} ROFs, {} clusters); publishing nothing and continuing with the next TimeFrame",
           rofsinput.size(), compClusters.size());
      invalidatePublication();
      return;
    }

    {
      mPublicationClock.emplace(mROFOverlapTable.getView().getClockLayer());
      const o2::itsmft::tracking::GenericTrackPublicationContext context{
        o2::detectors::DetID::ITS, o2::itsmft::tracking::ClusterSourceId{0},
        gsl::span<const o2::itsmft::ROFRecord>{rofsinput.data(), rofsinput.size()}, *mPublicationClock,
        kLayerToLayout,
        &mExternalIndicesBySurface, &mClusterSizesBySurface};
      o2::itsmft::tracking::GenericTrackOutputAdapterError error = o2::itsmft::tracking::GenericTrackOutputAdapterError::None;
      const auto staged = o2::itsmft::tracking::stageITSGenericTrackOutput(mFrame, context, mCompatibility, mUseMC, error);
      if (!staged) {
        throw std::runtime_error{"ITS GenericTrack output staging failed"};
      }

      auto& trackROFs = pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{"ITS", "ITSTrackROF", 0},
                                                                              staged->trackROFs.begin(), staged->trackROFs.end());
      auto& allTracksITS = pc.outputs().make<std::vector<o2::its::TrackITS>>(Output{"ITS", "TRACKS", 0});
      allTracksITS.assign(staged->tracks.begin(), staged->tracks.end());
      auto& allClusIdx = pc.outputs().make<std::vector<int>>(Output{"ITS", "TRACKCLSID", 0});
      allClusIdx.assign(staged->clusterIndices.begin(), staged->clusterIndices.end());
      LOGP(info, "ITS CA pushed {} tracks in {} ROFs", allTracksITS.size(), trackROFs.size());
      if (mUseMC) {
        pc.outputs().snapshot(Output{"ITS", "TRACKSMCTR", 0}, staged->labels);
        LOGP(info, "ITS CA pushed {} track MC labels", staged->labels.size());
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
    if (pc.inputs().getPos("itsTGeo") >= 0) {
      pc.inputs().get<o2::its::GeometryTGeo*>("itsTGeo");
    }
    pc.inputs().get<o2::itsmft::TopologyDictionary*>("itscldict");
    o2::its::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                                                o2::math_utils::TransformType::T2GRot,
                                                                                o2::math_utils::TransformType::T2G));
  }
}

void CATrackerDPL::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
    LOG(info) << "ITS CA input cluster dictionary updated";
    mDictionary = static_cast<const o2::itsmft::TopologyDictionary*>(obj);
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "GEOMTGEO", 0)) {
    LOG(info) << "ITS CA input GeometryTGeo loaded from CCDB";
    o2::its::GeometryTGeo::adopt(static_cast<o2::its::GeometryTGeo*>(obj));
    o2::its::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                                                o2::math_utils::TransformType::T2GRot,
                                                                                o2::math_utils::TransformType::T2G));
    // The catalog has static process lifetime; geometry adoption remains
    // necessary for raw cluster decoding.
    return;
  }
}

DataProcessorSpec getCATrackerSpec(bool useMC, bool useGeom, o2::itsmft::TrackingMode::Type trMode)
{
  std::vector<InputSpec> inputs;
  inputs.emplace_back("compClusters", "ITS", "COMPCLUSTERS", 0, Lifetime::Timeframe);
  inputs.emplace_back("patterns", "ITS", "PATTERNS", 0, Lifetime::Timeframe);
  inputs.emplace_back("ROframes", "ITS", "CLUSTERSROF", 0, Lifetime::Timeframe);
  inputs.emplace_back("itscldict", "ITS", "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec("ITS/Calib/ClusterDictionary"));

  if (useMC) {
    inputs.emplace_back("labels", "ITS", "CLUSTERSMCTR", 0, Lifetime::Timeframe);
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
    ggRequest->addInput({"itsTGeo", "ITS", "GEOMTGEO", 0, Lifetime::Condition, framework::ccdbParamSpec("ITS/Config/Geometry")}, inputs);
  }

  std::vector<OutputSpec> outputs;
  outputs.emplace_back("ITS", "TRACKS", 0, Lifetime::Timeframe);
  outputs.emplace_back("ITS", "TRACKCLSID", 0, Lifetime::Timeframe);
  outputs.emplace_back("ITS", "ITSTrackROF", 0, Lifetime::Timeframe);
  if (useMC) {
    outputs.emplace_back("ITS", "TRACKSMCTR", 0, Lifetime::Timeframe);
  }

  return DataProcessorSpec{
    "its-ca-tracker",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<CATrackerDPL>(ggRequest, useMC, trMode)},
    Options{}};
}

} // namespace o2::its::ca
