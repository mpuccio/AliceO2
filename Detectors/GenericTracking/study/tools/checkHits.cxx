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

/// \file checkHits.cxx
/// \brief Check that the MC hits of the ITS, MFT or ALICE 3 TRK lie on the
///        silicon of their chips as the tracking places them.
///
/// Loads the geometry, places every chip as GenericTracking/Alice/Placement.h
/// does and takes each MC hit into its chip's frame: both ends of every hit
/// must be on its chip's silicon. It also reports how far into the bounds the
/// hits reach, whether the bounds hold silicon no hit uses. Exits 1 if an end
/// lies off its chip's silicon by more than its float precision; the reach
/// only informs.

#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Study/Alice/Hits.h"
#include "GenericTracking/Study/Statistics.h"
#include <boost/program_options.hpp>
#include "GenericTracking/Tools/SensorFrame.h"

#include "DetectorsBase/GeometryManager.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <format>
#include <iostream>
#include <limits>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace bpo = boost::program_options;

namespace o2::gtrk::study
{

namespace
{

/// Both ends of every MC hit in its chip's frame -- where Geant4 entered and
/// left the silicon -- how many are off it by more than their own precision,
/// and by how much at most [cm], per layer. O2 stores a hit's ends as floats,
/// so an end on the silicon's surface can come out beyond it by the float step
/// of its largest coordinate. (Not the hit's middle: on a curved sensor the
/// chord of a grazing track sags inside the inner surface.) And per layer the
/// ends' reach into the bounds, |loc0| over the half width and |loc1| over
/// the half length in the unrolled frame: the 99.9th percentile and the
/// largest.
template <DetectorLayout Det>
int checkHits(const Geometry<Det>& geometry, study::Chips chips, const std::string& file, const std::string& discFile)
{
  study::HitFile hits(file, chips, discFile);
  std::vector<long> total(Detector<Det>::nLayers, 0);
  std::vector<long> off(Detector<Det>::nLayers, 0);
  std::vector<double> worst(Detector<Det>::nLayers, 0.);
  std::vector<std::vector<float>> across(Detector<Det>::nLayers);
  std::vector<std::vector<float>> along(Detector<Det>::nLayers);
  for (int event = 0; event < hits.nEvents(); ++event) {
    for (const auto& segment : hits.event(event)) {
      const auto& sensor = geometry.sensor(segment.chip);
      const auto& spec = geometry.specOf(segment.chip);
      for (const auto& end : {segment.start, segment.end}) {
        double frame[3];
        tools::toLocal(sensor, end.data(), frame);
        double unrolled[3];
        tools::frameToUnrolled(spec, frame, unrolled);
        const double beyond = std::max({std::abs(unrolled[0]) - spec.halfX, std::abs(unrolled[1]) - spec.halfY, spec.surfaceMin() - unrolled[2], unrolled[2] - spec.surfaceMax(), 0.});
        const double precision = static_cast<double>(std::numeric_limits<float>::epsilon()) * std::max({std::abs(end[0]), std::abs(end[1]), std::abs(end[2])});
        ++total[sensor.layer];
        across[sensor.layer].push_back(static_cast<float>(std::abs(unrolled[0]) / spec.halfX));
        along[sensor.layer].push_back(static_cast<float>(std::abs(unrolled[1]) / spec.halfY));
        if (beyond > precision) {
          ++off[sensor.layer];
          worst[sensor.layer] = std::max(worst[sensor.layer], beyond);
        }
      }
    }
  }
  std::cout << std::format("\nMC hits of {}{}, both ends of each in its chip's frame:\n", file, discFile.empty() ? "" : " and " + discFile);
  const auto reach = [](std::vector<float>& values) -> std::pair<float, float> {
    if (values.empty()) {
      return {0.F, 0.F};
    }
    return {study::quantile(std::span{values}, 0.999), *std::ranges::max_element(values)};
  };
  int status = 0;
  for (std::size_t layer = 0; layer < total.size(); ++layer) {
    const auto [acrossQ, acrossMax] = reach(across[layer]);
    const auto [alongQ, alongMax] = reach(along[layer]);
    std::cout << std::format(
      "layer {}: {} hit ends, {} off the silicon beyond their float precision, by at most {:.3f} um; reach into the bounds, "
      "99.9 % / largest: |loc0| / half width {:.4f} / {:.4f}, |loc1| / half length {:.4f} / {:.4f}\n",
      layer, total[layer], off[layer], worst[layer] * 1e4, acrossQ, acrossMax, alongQ, alongMax);
    status = off[layer] > 0 ? 1 : status;
  }
  return status;
}

template <DetectorLayout Det>
int runFor(study::Chips chips, const std::string& file, const std::string& discFile)
{
  const Geometry<Det> geometry(alice::placeChecked<Det>());
  return checkHits(geometry, chips, file, discFile);
}

int run(int argc, char** argv)
{
  std::string geometry;
  std::string detector;
  std::string hits;
  std::string ft3Hits;
  bpo::options_description description("Check that the MC hits of the ITS, MFT or ALICE 3 TRK lie on their chips' silicon");
  description.add_options()                                                                                    //
    ("help,h", "this message")                                                                                 //
    ("geometry,g", bpo::value(&geometry)->default_value("o2sim_geometry-aligned.root"), "o2sim geometry file") //
    ("detector,d", bpo::value(&detector)->default_value("ITS"), "ITS, MFT or TRKFT3")                          //
    ("hits", bpo::value(&hits)->required(), "the detector's MC hits, e.g. sgn_HitsTRK.root")                   //
    ("ft3-hits", bpo::value(&ft3Hits)->default_value("sgn_HitsFT3.root"), "TRKFT3: FT3's MC hits, checked too");
  bpo::variables_map vm;
  try {
    bpo::store(bpo::parse_command_line(argc, argv, description), vm);
    if (vm.contains("help")) {
      std::cout << description << '\n';
      return 0;
    }
    bpo::notify(vm);
  } catch (const bpo::error& error) {
    std::cerr << error.what() << "\n\n"
              << description << '\n';
    return 1;
  }
  o2::base::GeometryManager::loadGeometry(geometry);
  if (detector == "ITS") {
    return runFor<detectors::ITS>(study::Chips::kITS, hits, "");
  }
  if (detector == "MFT") {
    return runFor<detectors::MFT>(study::Chips::kMFT, hits, "");
  }
  if (detector == "TRKFT3") {
    return runFor<detectors::TRKFT3>(study::Chips::kTRK, hits, ft3Hits);
  }
  std::cerr << "unknown detector '" << detector << "'\n";
  return 1;
}

} // namespace

} // namespace o2::gtrk::study

int main(int argc, char** argv)
{
  try {
    return o2::gtrk::study::run(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
