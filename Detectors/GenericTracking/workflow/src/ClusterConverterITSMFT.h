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

/// \file ClusterConverterITSMFT.h
/// \brief The cluster converter of the ITS and the MFT, one format
///        (`o2::itsmft::CompClusterExt`, Clusters.h): ITS.cxx's and MFT.cxx's.

#ifndef GTRK_WORKFLOW_CLUSTERCONVERTERITSMFT_H_
#define GTRK_WORKFLOW_CLUSTERCONVERTERITSMFT_H_

#include "ITSMFTInput.h"

#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Clustering/Cluster.h"
#include "GenericTracking/Workflow/Conditions.h"
#include "GenericTracking/Workflow/DataSpecs.h"

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Logger.h"
#include "Framework/Task.h"
#include "Framework/TimingInfo.h"
#include "SimulationDataFormat/ConstMCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include <gsl/span>

#include <span>
#include <string>
#include <vector>

namespace o2::gtrk::workflow
{

/// O2's compact clusters of the ITS or the MFT, as the O2 cluster reader
/// publishes them, turned into GTRK/<DET>CLUSTERS for the measurement producer
/// and GTRK/<DET>CLUSMCTR for the MC truth. The order is O2's, so a
/// measurement's cluster index addresses streams concatenated in layer order.
template <DetectorLayout Det>
class ClusterConverterTask : public o2::framework::Task
{
 public:
  ClusterConverterTask(bool useMC, bool staggered) : mUseMC(useMC), mStreams(staggered ? Detector<Det>::nLayers : 1) {}

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final
  {
    if (matcher == o2::framework::ConcreteDataMatcher(Conditions<Det>::origin(), "ALPIDEPARAM", 0)) {
      return; // the CCDB fetch updates the singleton
    }
    if (matcher == o2::framework::ConcreteDataMatcher(Conditions<Det>::origin(), "CLUSDICT", 0)) {
      mDictionary = static_cast<const o2::itsmft::TopologyDictionary*>(obj);
      LOGP(info, "{} topology dictionary from the CCDB: {} patterns", Det::name, mDictionary->getSize());
    }
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    pc.inputs().get<o2::itsmft::TopologyDictionary*>("dictionary"); // delivers a new one through finaliseCCDB
    using AlpideParam = o2::itsmft::DPLAlpideParam<DetIDOf<Det>::value>;
    pc.inputs().get<AlpideParam*>("alpideParam");
    const auto input = readITSMFTInput(pc, mUseMC, mStreams);
    const auto& param = AlpideParam::Instance();
    const auto orbit = pc.services().get<o2::framework::TimingInfo>().firstTForbit;

    std::vector<FrameTiming> timings;
    timings.reserve(Detector<Det>::nLayers);
    for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
      timings.push_back({.lengthBC = param.getROFLengthInBC(layer), .firstStartBC = param.getROFDelayInBC(layer) + param.getROFBiasInBC(layer)});
      if (mStreams == 1 && timings.back() != timings.front()) {
        throw std::runtime_error("staggered ITS/MFT timing needs per-layer cluster input (--enable-its-staggering or --enable-mft-staggering)");
      }
    }
    const auto clusters = convertITSMFTInput(input, *mDictionary, timings, orbit);
    std::size_t rofCount = 0U;
    for (const auto& stream : input.streams()) {
      rofCount += stream.rofs.size();
    }
    LOGP(info, "{}: {} clusters in {} readout frames from {} streams for the timeframe at orbit {}", Det::name, clusters.size(), rofCount, mStreams, orbit);
    pc.outputs().snapshot(o2::framework::Output{gDataOrigin, description<Det>(kClusters), 0}, clusters);

    if (mUseMC) {
      auto& flat = pc.outputs().make<o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel>>(o2::framework::Output{gDataOrigin, description<Det>(kClustersMCTruth), 0});
      concatenateLabels(input).flatten_to(flat);
    }
  }

 private:
  bool mUseMC;
  int mStreams;
  const o2::itsmft::TopologyDictionary* mDictionary = nullptr; ///< the CCDB cache's, from the run's ITS/MFT Calib/ClusterDictionary
};

/// The converter's device for the ITS or the MFT.
template <DetectorLayout Det>
o2::framework::DataProcessorSpec getITSMFTClusterConverterSpec(bool useMC, bool staggered)
{
  using namespace o2::framework;
  const auto origin = Conditions<Det>::origin();
  const int streams = staggered ? Detector<Det>::nLayers : 1;
  std::vector<InputSpec> inputs;
  addITSMFTClusterInputs(inputs, origin, useMC, streams);
  const auto detector = origin.template as<std::string>();
  inputs.emplace_back("dictionary", origin, "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec(detector + "/Calib/ClusterDictionary"));
  inputs.emplace_back("alpideParam", origin, "ALPIDEPARAM", 0, Lifetime::Condition, ccdbParamSpec(detector + "/Config/AlpideParam"));
  std::vector<OutputSpec> outputs{OutputSpec{gDataOrigin, description<Det>(kClusters), 0, Lifetime::Timeframe}};
  if (useMC) {
    outputs.emplace_back(gDataOrigin, description<Det>(kClustersMCTruth), 0, Lifetime::Timeframe);
  }
  return DataProcessorSpec{
    .name = deviceName<Det>("cluster-converter"),
    .inputs = inputs,
    .outputs = outputs,
    .algorithm = AlgorithmSpec{adaptFromTask<ClusterConverterTask<Det>>(useMC, staggered)},
    .options = Options{},
  };
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_CLUSTERCONVERTERITSMFT_H_
