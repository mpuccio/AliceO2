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

#include <utility>
#include <vector>

#include <gsl/span>

#include "CommonDataFormat/IRFrame.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DetectorsBase/GeometryManager.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Logger.h"
#include "ITSMFTTracking/CATracker.h"
#include "ITSMFTTracking/MFTCATrack.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTTracking/Constants.h"
#include "MFTTracking/MFTTrackingParam.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::framework;

namespace o2::mft
{

static_assert(o2::itsmft::tracking::ITSMFTTrackingInterfaceMFT::DetId == o2::detectors::DetID::MFT);

namespace
{
template <typename TracksVec, typename ClusterIdxVec, typename ROFVec, typename LabelsVec, typename SeedPatternVec>
void fillMFTOutputs(const o2::itsmft::tracking::TimeFrameMFT& tf,
                    gsl::span<const o2::itsmft::ROFRecord> inputROFs,
                    TracksVec& tracks,
                    ClusterIdxVec& clusterIndices,
                    ROFVec& trackROFs,
                    LabelsVec& trackLabels,
                    SeedPatternVec& seedPatterns,
                    bool useMC)
{
  trackROFs.assign(inputROFs.begin(), inputROFs.end());
  for (auto& rof : trackROFs) {
    rof.setFirstEntry(0);
    rof.setNEntries(0);
  }

  const auto& tracksIn = tf.getTracks();
  tracks.reserve(tracksIn.size());
  seedPatterns.reserve(tracksIn.size());
  if (useMC) {
    trackLabels.reserve(tracksIn.size());
  }

  const auto& clockLayer = tf.getROFOverlapTableView().getClockLayer();
  std::vector<int> rofEntries(trackROFs.size() + 1, 0);

  for (size_t iTrk = 0; iTrk < tracksIn.size(); ++iTrk) {
    const auto& src = tracksIn[iTrk];
    auto dst = src.getTrack();
    dst.setExternalClusterIndexOffset(clusterIndices.size());
    int nPoints = 0;
    for (int layer = o2::itsmft::tracking::MFTCATrack::MaxClusters; layer--;) {
      if (!src.hasHitOnLayer(layer)) {
        continue;
      }
      const int extIdx = src.getClusterIndex(layer);
      if (extIdx < 0) {
        continue;
      }
      // src's per-layer cluster size was already captured correctly by
      // Tracker<NLayers>::rectifyClusterIndices() (CATracker.cxx) while the
      // layer-local cluster identity was still available; `extIdx` here is
      // already the external/global identity, and mClusterSize is a
      // layer-local vector, so it must never be re-queried with `extIdx`.
      dst.setClusterSize(layer, src.getClusterSize(layer));
      clusterIndices.push_back(extIdx);
      ++nPoints;
    }
    dst.setNumberOfPoints(nPoints);
    tracks.push_back(dst);
    seedPatterns.push_back(src.getSeedPattern());

    if (useMC && iTrk < tf.getTracksLabel().size()) {
      trackLabels.push_back(tf.getTracksLabel()[iTrk]);
    }

    const auto rof = clockLayer.getROF(src.getTimeStamp());
    if (rof >= 0 && rof < static_cast<int>(trackROFs.size())) {
      ++rofEntries[rof];
    }
  }

  std::exclusive_scan(rofEntries.begin(), rofEntries.end(), rofEntries.begin(), 0);
  for (size_t iROF = 0; iROF < trackROFs.size(); ++iROF) {
    trackROFs[iROF].setFirstEntry(rofEntries[iROF]);
    trackROFs[iROF].setNEntries(rofEntries[iROF + 1] - rofEntries[iROF]);
  }
}
} // namespace

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
  // value; see CATracker.h/CATracker.cxx for the classification.
  const float trackingResult = mTracking.processTimeFrame(gsl::span<const o2::itsmft::ROFRecord>(rofsinput.data(), rofsinput.size()),
                                                          gsl::span<const o2::itsmft::CompClusterExt>(compClusters.data(), compClusters.size()),
                                                          patterns,
                                                          labels,
                                                          irFrames);

  if (decideCATrackerPublicationAction(mTracking.isActive(), trackingResult) == CATrackerPublicationAction::SkipDroppedTimeFrame) {
    LOGP(error, "MFT CA tracking dropped this TimeFrame ({} ROFs, {} clusters); publishing nothing and continuing with the next TimeFrame",
         rofsinput.size(), compClusters.size());
    return;
  }

  auto& trackROFs = pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{"MFT", "MFTTrackROF", 0},
                                                                          rofsinput.begin(), rofsinput.end());
  auto& allTracksMFT = pc.outputs().make<std::vector<o2::mft::TrackMFT>>(Output{"MFT", "TRACKS", 0});
  auto& allClusIdx = pc.outputs().make<std::vector<int>>(Output{"MFT", "TRACKCLSID", 0});
  auto& allSeedPatterns = pc.outputs().make<std::vector<uint16_t>>(Output{"MFT", "TRACKSEEDPAT", 0});
  std::vector<o2::MCCompLabel> allTrackLabels;

  fillMFTOutputs(mTracking.getTimeFrame(),
                 gsl::span<const o2::itsmft::ROFRecord>(rofsinput.data(), rofsinput.size()),
                 allTracksMFT,
                 allClusIdx,
                 trackROFs,
                 allTrackLabels,
                 allSeedPatterns,
                 mUseMC);

  LOGP(info, "MFT CA pushed {} tracks in {} ROFs", allTracksMFT.size(), trackROFs.size());

  if (mUseMC) {
    pc.outputs().snapshot(Output{"MFT", "TRACKSMCTR", 0}, allTrackLabels);
    LOGP(info, "MFT CA pushed {} track MC labels", allTrackLabels.size());
  }

  mTracking.clearTimeFrame();
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
