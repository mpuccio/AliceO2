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
/// measurement's cluster index is O2's cluster index.
template <DetectorLayout Det>
class ClusterConverterTask : public o2::framework::Task
{
 public:
  explicit ClusterConverterTask(bool useMC) : mUseMC(useMC) {}

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final
  {
    if (matcher == o2::framework::ConcreteDataMatcher(Conditions<Det>::origin(), "CLUSDICT", 0)) {
      mDictionary = static_cast<const o2::itsmft::TopologyDictionary*>(obj);
      LOGP(info, "{} topology dictionary from the CCDB: {} patterns", Det::name, mDictionary->getSize());
    }
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    pc.inputs().get<o2::itsmft::TopologyDictionary*>("dictionary"); // delivers a new one through finaliseCCDB
    const auto compact = pc.inputs().get<gsl::span<o2::itsmft::CompClusterExt>>("compclusters");
    const auto patterns = pc.inputs().get<gsl::span<unsigned char>>("patterns");
    const auto rofs = pc.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("rofs");
    const auto& timing = pc.services().get<o2::framework::TimingInfo>();

    std::vector<Cluster> clusters;
    if (!rofs.empty()) {
      const auto local = alice::localClusters(std::span{compact.data(), compact.size()}, std::span{patterns.data(), patterns.size()}, *mDictionary);
      const std::span<const o2::itsmft::ROFRecord> frames{rofs.data(), rofs.size()};
      clusters = alice::toClusters(local, frames, alice::readoutFrames(frames), timing.firstTForbit);
    }
    LOGP(info, "{}: {} clusters in {} readout frames for the timeframe at orbit {}", Det::name, clusters.size(), rofs.size(), timing.firstTForbit);
    pc.outputs().snapshot(o2::framework::Output{gDataOrigin, description<Det>(kClusters), 0}, clusters);

    if (mUseMC) {
      const auto labels = pc.inputs().get<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>("labels");
      if (labels->getIndexedSize() != clusters.size() && !clusters.empty()) {
        LOGP(fatal, "{}: {} labelled clusters for {} clusters", Det::name, labels->getIndexedSize(), clusters.size());
      }
      auto& flat = pc.outputs().make<o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel>>(o2::framework::Output{gDataOrigin, description<Det>(kClustersMCTruth), 0});
      labels->flatten_to(flat);
    }
  }

 private:
  bool mUseMC;
  const o2::itsmft::TopologyDictionary* mDictionary = nullptr; ///< the CCDB cache's, from the run's ITS/MFT Calib/ClusterDictionary
};

/// The converter's device for the ITS or the MFT.
template <DetectorLayout Det>
o2::framework::DataProcessorSpec getITSMFTClusterConverterSpec(bool useMC)
{
  using namespace o2::framework;
  const auto origin = Conditions<Det>::origin();
  std::vector<InputSpec> inputs{InputSpec{"compclusters", origin, "COMPCLUSTERS", 0, Lifetime::Timeframe}, InputSpec{"patterns", origin, "PATTERNS", 0, Lifetime::Timeframe}, InputSpec{"rofs", origin, "CLUSTERSROF", 0, Lifetime::Timeframe}, InputSpec{"dictionary", origin, "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec(origin.template as<std::string>() + "/Calib/ClusterDictionary")}};
  std::vector<OutputSpec> outputs{OutputSpec{gDataOrigin, description<Det>(kClusters), 0, Lifetime::Timeframe}};
  if (useMC) {
    inputs.emplace_back("labels", origin, "CLUSTERSMCTR", 0, Lifetime::Timeframe);
    outputs.emplace_back(gDataOrigin, description<Det>(kClustersMCTruth), 0, Lifetime::Timeframe);
  }
  return DataProcessorSpec{
    .name = deviceName<Det>("cluster-converter"),
    .inputs = inputs,
    .outputs = outputs,
    .algorithm = AlgorithmSpec{adaptFromTask<ClusterConverterTask<Det>>(useMC)},
    .options = Options{},
  };
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_CLUSTERCONVERTERITSMFT_H_
