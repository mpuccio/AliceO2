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
#include "ITSMFTTracking/CommonTrackOutputAdapter.h"
#include "ITSMFTTracking/ROFTimingUniformity.h"
#include "ITSMFTTracking/TimeFrameLoadFailure.h"
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

void CATrackerDPL::configureROFViews(gsl::span<const o2::itsmft::ROFRecord> rofs,
                                     gsl::span<const o2::dataformats::IRFrame> irFrames)
{
  const auto& params = mTracking.getTrackingParameters().at(0);
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
      .mROFAddTimeErr = static_cast<uint32_t>(params.AddTimeError[layer])};
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
  mPublicationAdapter.reset();
  mTracking.getScratch().setROFViews({mROFOverlapTable.getView(), mROFVertexLookupTable.getView(), mMultiplicityMask.getView(), mUPCMask.getView()});
}

void CATrackerDPL::invalidatePublication() noexcept
{
  mPublicationAdapter.reset();
  mPublicationClock.reset();
  mTracking.getScratch().setROFViews({});
}

static_assert(o2::itsmft::tracking::ITSMFTTrackingInterfaceMFT::DetId == o2::detectors::DetID::MFT);

CATrackerPublicationAction decideCATrackerPublicationAction(bool trackerActive, float trackingResult) noexcept
{
  if (!trackerActive) {
    return CATrackerPublicationAction::PublishInactiveEmpty;
  }
  if (o2::itsmft::tracking::isDroppedTimeFrame(trackingResult)) {
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

  if (decideCATrackerPublicationAction(mTracking.isActive(), 0.f) == CATrackerPublicationAction::PublishInactiveEmpty) {
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

  // A structural or unclassified failure inside processTimeFrame() throws
  // uncaught out of this function -- no output is created on that path, and
  // DPL treats the escaping exception as fatal for this device. Only a
  // recoverable, dropped-and-wiped TimeFrame returns here as a sentinel
  // value; see Tracker.h/Tracker.cxx for the classification.
  const auto resetCount = mTracking.getTimeFrame().getEventResetCount();
  try {
    configureROFViews(gsl::span<const o2::itsmft::ROFRecord>(rofsinput.data(), rofsinput.size()), irFrames);
    const float trackingResult = mTracking.processTimeFrame(gsl::span<const o2::itsmft::ROFRecord>(rofsinput.data(), rofsinput.size()),
                                                            gsl::span<const o2::itsmft::CompClusterExt>(compClusters.data(), compClusters.size()),
                                                            patterns,
                                                            labels,
                                                            irFrames,
                                                            mTracking.getScratch().getROFViews());

    if (decideCATrackerPublicationAction(mTracking.isActive(), trackingResult) == CATrackerPublicationAction::SkipDroppedTimeFrame) {
      LOGP(error, "MFT CA tracking dropped this TimeFrame ({} ROFs, {} clusters); publishing nothing and continuing with the next TimeFrame",
           rofsinput.size(), compClusters.size());
      invalidatePublication();
      return;
    }

    {
      mPublicationClock.emplace(mROFOverlapTable.getView().getClockLayer());
      const auto exportContext = mTracking.getCommonTrackPublicationExport(*mPublicationClock);
      const auto* compatibility = &mCompatibility;
      if (!exportContext || compatibility == nullptr) {
        throw std::runtime_error{"MFT CommonTrack output publication context is unavailable"};
      }
      const o2::itsmft::tracking::CommonTrackPublicationContext context{
        exportContext->detector, exportContext->source,
        gsl::span<const o2::itsmft::ROFRecord>{rofsinput.data(), rofsinput.size()}, exportContext->clock, exportContext->orderedSurfaces};
      o2::itsmft::tracking::CommonTrackOutputAdapterError error = o2::itsmft::tracking::CommonTrackOutputAdapterError::None;
      const auto staged = o2::itsmft::tracking::stageMFTCommonTrackOutput(mTracking.getTimeFrame(), context, *compatibility, mUseMC, error);
      if (!staged) {
        throw std::runtime_error{"MFT CommonTrack output staging failed"};
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
    if (mTracking.getTimeFrame().getEventResetCount() == resetCount) {
      mTracking.resetEvent();
    }
    invalidatePublication();
    throw;
  }

  mTracking.resetEvent();
  invalidatePublication();
}

void CATrackerDPL::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  if (!mTrackingInitialised) {
    mTrackingInitialised = true;
    mTracking.initialise();
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
    mTracking.setClusterDictionary(static_cast<const o2::itsmft::TopologyDictionary*>(obj));
    return;
  }
  if (matcher == ConcreteDataMatcher("MFT", "GEOMTGEO", 0)) {
    LOG(info) << "MFT CA input GeometryTGeo loaded from CCDB";
    o2::mft::GeometryTGeo::adopt(static_cast<o2::mft::GeometryTGeo*>(obj));
    o2::mft::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                                                o2::math_utils::TransformType::T2GRot,
                                                                                o2::math_utils::TransformType::T2G,
                                                                                o2::math_utils::TransformType::L2G));
    // Gate 4 B2 Slice 2: the tracking catalog is a static, process-lifetime
    // table (StaticDetectorCatalogs.h), immune to alignment/geometry updates
    // by design -- there is nothing left to invalidate here. GeometryTGeo
    // adoption above stays: raw cluster decoding still needs it.
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
