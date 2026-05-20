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

#include <vector>

#include <gsl/span>

#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/IRFrame.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DeviceSpec.h"
#include "Framework/Logger.h"
#include "MFTBase/GeometryTGeo.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::framework;

namespace o2::mft
{

void CATrackerDPL::init(InitContext&)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
}

void CATrackerDPL::run(ProcessingContext& pc)
{
  updateTimeDependentParams(pc);

  auto compClusters = pc.inputs().get<const std::vector<o2::itsmft::CompClusterExt>>("compClusters");
  gsl::span<const unsigned char> patterns = pc.inputs().get<gsl::span<unsigned char>>("patterns");
  auto rofs = pc.inputs().get<const std::vector<o2::itsmft::ROFRecord>>("ROframes");
  initialiseROFTable(gsl::span<const o2::itsmft::ROFRecord>(rofs.data(), rofs.size()));

  const dataformats::MCTruthContainer<MCCompLabel>* labels = nullptr;
  if (mUseMC) {
    labels = pc.inputs().get<const dataformats::MCTruthContainer<MCCompLabel>*>("labels").release();
  }

  LOGP(info, "MFT CA input pulled {} compressed clusters, {} pattern bytes, {} RO frames, dictionary={}, MC labels={}",
       compClusters.size(), patterns.size(), rofs.size(), mDict != nullptr, labels != nullptr);

  size_t nClusterRefs = 0;
  for (const auto& rof : rofs) {
    nClusterRefs += rof.getNEntries();
  }
  LOGP(info, "MFT CA input ROF cluster references: {}", nClusterRefs);

  if (pc.inputs().getPos("IRFramesITS") >= 0) {
    auto irFrames = pc.inputs().get<gsl::span<o2::dataformats::IRFrame>>("IRFramesITS");
    LOGP(info, "MFT CA input pulled {} ITS IR frames", irFrames.size());
  }

  // Next step: replace this inspection point with a loader into a generalized
  // ITS TimeFrame interface that accepts MFT layer geometry and cluster data.
}

void CATrackerDPL::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  static bool initOnceDone = false;
  if (!initOnceDone) {
    initOnceDone = true;
    if (pc.inputs().getPos("mftTGeo") >= 0) {
      pc.inputs().get<o2::mft::GeometryTGeo*>("mftTGeo");
    }
    pc.inputs().get<o2::itsmft::TopologyDictionary*>("cldict");
    o2::mft::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                                                 o2::math_utils::TransformType::T2GRot,
                                                                                 o2::math_utils::TransformType::T2G));
  }
}

void CATrackerDPL::initialiseROFTable(gsl::span<const o2::itsmft::ROFRecord> rofs)
{
  if (!o2::base::GRPGeomHelper::instance().getGRPECS()->isDetContinuousReadOut(o2::detectors::DetID::MFT)) {
    LOGP(fatal, "MFT CA tracker currently supports only continuous readout");
  }

  const auto& par = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::MFT>::Instance();
  const int nOrbitsPerTF = o2::base::GRPGeomHelper::getNHBFPerTF();
  const int timingSourceLayer = 0; // MFT CA disks share timing for now; keep per-CA-layer definition for future staggering.
  const unsigned int nROFsPerOrbit = o2::constants::lhc::LHCMaxBunches / par.getROFLengthInBC(timingSourceLayer);
  const auto nROFsTF = nROFsPerOrbit * nOrbitsPerTF;

  ROFOverlapTable rofTable;
  for (int caLayer = 0; caLayer < NCALayers; ++caLayer) {
    const o2::its::LayerTiming timing{
      .mNROFsTF = nROFsTF,
      .mROFLength = static_cast<uint32_t>(par.getROFLengthInBC(timingSourceLayer)),
      .mROFDelay = static_cast<uint32_t>(par.getROFDelayInBC(timingSourceLayer)),
      .mROFBias = static_cast<uint32_t>(par.getROFBiasInBC(timingSourceLayer)),
      .mROFAddTimeErr = 0};
    rofTable.defineLayer(caLayer, timing);
  }
  rofTable.init();
  mROFTable = std::move(rofTable);
  mROFTableView = mROFTable.getView();

  if (rofs.size() != nROFsTF) {
    LOGP(warn, "MFT CA ROF count differs from continuous timing expectation: received {} expected {}", rofs.size(), nROFsTF);
  }
  LOGP(info, "MFT CA ROF lookup table initialised for {} CA layers, {} ROFs/TF, ROF length {} BC, bias {} BC",
       NCALayers, nROFsTF, par.getROFLengthInBC(timingSourceLayer), par.getROFBiasInBC(timingSourceLayer));
}

void CATrackerDPL::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("MFT", "CLUSDICT", 0)) {
    LOG(info) << "MFT CA input cluster dictionary updated";
    mDict = static_cast<const o2::itsmft::TopologyDictionary*>(obj);
    return;
  }
  if (matcher == ConcreteDataMatcher("MFT", "GEOMTGEO", 0)) {
    LOG(info) << "MFT CA input GeometryTGeo loaded from CCDB";
    o2::mft::GeometryTGeo::adopt(static_cast<o2::mft::GeometryTGeo*>(obj));
    return;
  }
}

DataProcessorSpec getCATrackerSpec(bool useMC, bool useGeom, bool useIRFrames)
{
  std::vector<InputSpec> inputs;
  inputs.emplace_back("compClusters", "MFT", "COMPCLUSTERS", 0, Lifetime::Timeframe);
  inputs.emplace_back("patterns", "MFT", "PATTERNS", 0, Lifetime::Timeframe);
  inputs.emplace_back("ROframes", "MFT", "CLUSTERSROF", 0, Lifetime::Timeframe);
  inputs.emplace_back("cldict", "MFT", "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec("MFT/Calib/ClusterDictionary"));

  if (useMC) {
    inputs.emplace_back("labels", "MFT", "CLUSTERSMCTR", 0, Lifetime::Timeframe);
  }
  if (useIRFrames) {
    inputs.emplace_back("IRFramesITS", "ITS", "IRFRAMES", 0, Lifetime::Timeframe);
  }

  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                                                                        // orbitResetTime
                                                              true,                                                                         // GRPECS=true
                                                              false,                                                                        // GRPLHCIF
                                                              true,                                                                         // GRPMagField
                                                              false,                                                                        // askMatLUT
                                                              useGeom ? o2::base::GRPGeomRequest::Aligned : o2::base::GRPGeomRequest::None, // geometry
                                                              inputs,
                                                              true);
  if (!useGeom) {
    ggRequest->addInput({"mftTGeo", "MFT", "GEOMTGEO", 0, Lifetime::Condition, framework::ccdbParamSpec("MFT/Config/Geometry")}, inputs);
  }

  return DataProcessorSpec{
    "mft-ca-tracker",
    inputs,
    {},
    AlgorithmSpec{adaptFromTask<CATrackerDPL>(ggRequest, useMC)},
    Options{}};
}

} // namespace o2::mft
