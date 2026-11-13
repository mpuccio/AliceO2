// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// SPDX-License-Identifier: GPL-3.0-only

// Run the native IOTOF digitizer and clusterer on an existing production's
// hits and collision context. No replacement/smearing model is used here.
#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsIOTOF/Cluster.h"
#include "DetectorsBase/GeometryManager.h"
#include "IOTOFSimulation/Digitizer.h"
#include "IOTOFSimulation/DPLDigitizerParam.h"
#include "IOTOFReconstruction/Clusterer.h"
#include "SimulationDataFormat/DigitizationContext.h"
#include "SimulationDataFormat/IOMCTruthContainerView.h"
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>

int main(int argc, char** argv)
{
  try {
    namespace bpo = boost::program_options;
    std::string config, geometry, contextPath, hitPath, output;
    unsigned seed = 42, firstOrbit = 0;
    bpo::options_description options("Native TF3 digitization and clustering of frozen simulation inputs");
    options.add_options()("help,h", "show options")                                                          //
      ("production-config", bpo::value(&config)->required(), "simulation configuration INI")                 //
      ("geometry", bpo::value(&geometry)->required(), "same aligned geometry")                               //
      ("context", bpo::value(&contextPath)->required(), "same collision context")                            //
      ("tof-hits", bpo::value(&hitPath)->required(), "source-0 TF3 hit file")                                //
      ("output-dir", bpo::value(&output)->required(), "new output directory; refuses existing ROOT outputs") //
      ("seed", bpo::value(&seed)->default_value(42), "digitization seed")                                    //
      ("first-orbit", bpo::value(&firstOrbit)->default_value(0), "timeframe first orbit");
    bpo::variables_map args;
    bpo::store(bpo::parse_command_line(argc, argv, options), args);
    if (args.contains("help")) {
      std::cout << options << '\n';
      return 0;
    }
    bpo::notify(args);
    const auto directory = std::filesystem::absolute(output);
    std::filesystem::create_directories(directory);
    if (std::filesystem::exists(directory / "tf3digits.root") || std::filesystem::exists(directory / "tf3clusters.root")) {
      throw std::runtime_error("output already exists; use a new directory");
    }
    o2::conf::ConfigurableParam::updateFromFile(config, "IOTOFBase,ChipSpecificsParam,TF3DigitizerParam");
    o2::base::GeometryManager::loadGeometry(geometry, false, true);
    auto* geom = o2::iotof::GeometryTGeo::Instance();
    geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
    gRandom->SetSeed(seed);
    std::unique_ptr<o2::steer::DigitizationContext> context(o2::steer::DigitizationContext::loadFromFile(contextPath));
    std::unique_ptr<TFile> file(TFile::Open(hitPath.c_str()));
    auto* tree = file ? file->Get<TTree>("o2sim") : nullptr;
    std::vector<o2::itsmft::Hit> hits;
    auto* hitPointer = &hits;
    if (!context || !tree || tree->SetBranchAddress("TF3Hit", &hitPointer) < 0) {
      throw std::runtime_error("missing context or TF3 hit tree");
    }
    std::vector<o2::iotof::Digit> digits;
    std::vector<o2::itsmft::ROFRecord> rofs;
    o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
    o2::iotof::Digitizer digitizer;
    digitizer.setGeometry(geom);
    digitizer.setContinuous(o2::iotof::DPLDigitizerParam::Instance().continuous);
    if (!digitizer.isContinuous()) {
      throw std::runtime_error("this single-timeframe driver requires continuous TF3 readout");
    }
    digitizer.setDigits(&digits);
    digitizer.setMCLabels(&labels);
    digitizer.setROFRecords(&rofs);
    digitizer.setROFRecordIR(o2::InteractionRecord(0, firstOrbit));
    digitizer.init();
    const auto& times = context->getEventRecords();
    const auto& parts = context->getEventParts();
    std::size_t events = 0, hitCount = 0;
    std::set<int> sourceEvents;
    for (std::size_t collision = 0; collision < times.size(); ++collision) {
      digitizer.setEventTime(times[collision]);
      for (const auto& part : parts[collision]) {
        if (part.sourceID != 0) {
          throw std::runtime_error("provide a source-0-only context for this driver");
        }
        if (!sourceEvents.insert(part.entryID).second) {
          throw std::runtime_error("reused source events require occurrence-aware study identities");
        }
        if (part.entryID < 0 || part.entryID >= tree->GetEntries() || tree->GetEntry(part.entryID) <= 0) {
          throw std::runtime_error("cannot read a context-referenced TF3 hit event");
        }
        digitizer.process(&hits, part.entryID, part.sourceID);
        ++events;
        hitCount += hits.size();
      }
    }
    digitizer.fillOutputContainer();
    o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel> flat;
    labels.flatten_to(flat);
    o2::dataformats::ConstMCTruthContainerView<o2::MCCompLabel> view(flat);
    std::vector<o2::iotof::Cluster> clusters;
    std::vector<unsigned char> patterns;
    std::vector<o2::itsmft::ROFRecord> clusterROFs;
    o2::dataformats::MCTruthContainer<o2::MCCompLabel> clusterLabels;
    o2::iotof::Clusterer clusterer;
    std::filesystem::current_path(directory); // native topology-map output
    clusterer.process(digits, rofs, clusters, patterns, clusterROFs, &view, &clusterLabels);
    if (labels.getIndexedSize() != digits.size() || clusterLabels.getIndexedSize() != clusters.size()) {
      throw std::runtime_error("native MC label rows do not match measurements");
    }
    {
      TFile result((directory / "tf3digits.root").c_str(), "RECREATE");
      TTree outputTree("o2sim", "Native IOTOF digits");
      o2::dataformats::IOMCTruthContainerView ioLabels(flat);
      std::vector<o2::itsmft::MC2ROFRecord> mc2rof;
      outputTree.Branch("TF3Digit", &digits);
      outputTree.Branch("TF3DigitROF", &rofs);
      outputTree.Branch("TF3DigitMCTruth", &ioLabels);
      outputTree.Branch("TF3DigitMC2ROF", &mc2rof);
      outputTree.Fill();
      outputTree.Write();
    }
    {
      TFile result((directory / "tf3clusters.root").c_str(), "RECREATE");
      TTree outputTree("o2sim", "Native IOTOF clusters with BC/fine timing");
      outputTree.Branch("TF3Cluster", &clusters);
      outputTree.Branch("TF3ClusterPatt", &patterns);
      outputTree.Branch("TF3ClusterROF", &clusterROFs);
      outputTree.Branch("TF3ClusterMCTruth", &clusterLabels);
      outputTree.Fill();
      outputTree.Write();
    }
    o2::conf::ConfigurableParam::writeINI((directory / "tof-config.ini").string());
    nlohmann::ordered_json report{{"production_config", config}, {"geometry", geometry}, {"context", contextPath}, {"tof_hits", hitPath}, {"first_orbit", firstOrbit}, {"events", events}, {"transport_hits", hitCount}, {"digits", digits.size()}, {"digit_label_rows", labels.getIndexedSize()}, {"clusters", clusters.size()}, {"cluster_label_rows", clusterLabels.getIndexedSize()}, {"digit_rofs", rofs.size()}, {"cluster_rofs", clusterROFs.size()}, {"seed", seed}, {"response_maps", "disabled, matching the native workflow default; flat response from production config"}};
    std::ofstream(directory / "preparation.json") << report.dump(2) << '\n';
    std::cout << report.dump(2) << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
