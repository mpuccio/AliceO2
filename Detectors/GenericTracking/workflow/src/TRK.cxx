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

/// \file TRK.cxx
/// \brief ALICE 3's own, TRK's barrel and FT3's discs: its chips placed, its
///        cluster format, its dictionaries' depths, frames and field (TRK.h),
///        its cluster converter and its measurement producer.

#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Alice/TRK.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Clustering/Cluster.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Workflow/ClusterConverterSpec.h"
#include "GenericTracking/Workflow/ConditionsTRK.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/MeasurementProducerSpec.h"

#include "MeasurementProducer.h"

#include "DataFormatsTRKFT3/Cluster.h"
#include "DataFormatsTRKFT3/ROFRecord.h"
#include "DetectorsRaw/HBFUtils.h"
#include "FT3Base/GeometryTGeo.h"
#include "Field/ALICE3MagneticField.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Logger.h"
#include "Framework/Task.h"
#include "Framework/TimingInfo.h"
#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"
#include "SimulationDataFormat/ConstMCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TRKBase/AlmiraParam.h"
#include "TRKBase/GeometryTGeo.h"
#include "TRKBase/Specs.h"
#include "TRKReconstruction/TopologyDictionary.h"

#include <gsl/span>

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <initializer_list>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace o2::gtrk::alice
{

namespace
{

// TRKFT3.h restates O2's TRK sensors; here they are held to O2's, to the bit.
namespace specs = o2::trk::constants;
static_assert(detectors::trkft3::kVDRadius == specs::VD::petal::layer::radii);
static_assert(detectors::trkft3::kVDPetalGap == specs::VD::petal::layer::gaps);
static_assert(exactlyEqual(detectors::trkft3::vdWidth(0), specs::VD::petal::layer::width[0]));
static_assert(exactlyEqual(detectors::trkft3::vdWidth(1), specs::VD::petal::layer::width[1]));
static_assert(exactlyEqual(detectors::trkft3::vdWidth(2), specs::VD::petal::layer::width[2]));
static_assert(exactlyEqual(detectors::trkft3::kVDLength, specs::VD::petal::layer::length));
static_assert(exactlyEqual(detectors::trkft3::kVDThickness, specs::VD::silicon::thickness));
static_assert(exactlyEqual(detectors::trkft3::kVDPitch, specs::VD::petal::layer::pitchX) && exactlyEqual(detectors::trkft3::kVDPitch, specs::VD::petal::layer::pitchZ));
static_assert(exactlyEqual(detectors::trkft3::kChipWidth, specs::moduleMLOT::chip::width));
static_assert(exactlyEqual(detectors::trkft3::kChipPassiveEdgeReadOut, specs::moduleMLOT::chip::passiveEdgeReadOut));
static_assert(exactlyEqual(detectors::trkft3::kChipLength, specs::moduleMLOT::chip::length));
static_assert(exactlyEqual(detectors::trkft3::kChipThickness, specs::moduleMLOT::silicon::thickness));
static_assert(exactlyEqual(detectors::trkft3::kChipPitch, specs::moduleMLOT::chip::pitchX) && exactlyEqual(detectors::trkft3::kChipPitch, specs::moduleMLOT::chip::pitchZ));

} // namespace

template <>
std::vector<Sensor> placeChips<detectors::TRKFT3>()
{
  // The VD's curved sensors (sub-detector 0) and the ML/OT chips, whose
  // sensor volume is the active matrix, so their origin needs no shift; then
  // FT3's chips, the same chip with the same volume, their layers after the
  // barrel's.
  auto& trk = *o2::trk::GeometryTGeo::Instance();
  auto& ft3 = *o2::ft3::GeometryTGeo::Instance();
  trk.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  ft3.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  std::vector<Sensor> sensors;
  sensors.reserve(static_cast<std::size_t>(trk.getNumberOfChips()) + static_cast<std::size_t>(ft3.getNumberOfChips()));
  for (int chip = 0; chip < trk.getNumberOfChips(); ++chip) {
    const auto& l2g = trk.getMatrixL2G(chip);
    const int layer = trk.getLayerTRK(chip);
    sensors.push_back(trk.getSubDetID(chip) == 0 ? curvedSensor(l2g, layer) : flatSensor(l2g, layer, 0.));
  }
  const int firstDiscLayer = Detector<detectors::TRKFT3>::nBarrelLayers;
  for (int chip = 0; chip < ft3.getNumberOfChips(); ++chip) {
    sensors.push_back(flatSensor(ft3.getMatrixL2G(chip), firstDiscLayer + ft3.getLayer(chip), 0.));
  }
  return sensors;
}

namespace
{

/// The local clusters of TRK's or FT3's \a clusters, each chip numbered
/// \a firstSensor plus O2's chip ID; \a subDetectors, the sub-detectors their
/// clusters may name (a cluster naming another is refused).
template <typename ClusterType, typename Dictionary>
std::vector<LocalCluster> alice3LocalClusters(std::span<const ClusterType> clusters, std::span<const unsigned char> patterns, const Dictionary& dictionary, std::uint32_t firstSensor, std::initializer_list<int> subDetectors)
{
  std::vector<LocalCluster> result;
  result.reserve(clusters.size());
  std::size_t offset = 0;
  for (const auto& cluster : clusters) {
    if (std::ranges::find(subDetectors, static_cast<int>(cluster.subDetID)) == subDetectors.end()) {
      throw std::runtime_error(std::format("a cluster on chip {} names sub-detector {}", cluster.chipID, cluster.subDetID));
    }
    // The pattern: its row span, its column span, then its pixels' bits row
    // by row, the first pixel in the highest bit (O2's Clusterer<TRK>).
    if (offset + 2U > patterns.size()) {
      throw std::runtime_error("TRK cluster patterns ran out");
    }
    const auto* pattern = patterns.data() + offset;
    const auto rowSpan = pattern[0];
    const auto colSpan = pattern[1];
    const auto bytes = 2U + ((static_cast<std::size_t>(rowSpan) * colSpan + 7U) / 8U);
    if (offset + bytes > patterns.size()) {
      throw std::runtime_error("TRK cluster patterns ran out");
    }
    offset += bytes;
    int fired = 0;
    for (std::size_t byte = 2; byte < bytes; ++byte) {
      fired += std::popcount(pattern[byte]);
    }
    LocalCluster local{.chip = firstSensor + cluster.chipID, .nPixels = static_cast<std::uint32_t>(fired), .rowSpan = rowSpan, .colSpan = colSpan};
    const int layer = Dictionary::layerOf(cluster);
    if (layer < 0 || layer >= dictionary.getNumberOfLayers()) {
      throw std::runtime_error(std::format("a cluster on chip {} is on no layer of the topology dictionary", cluster.chipID));
    }
    const int id = dictionary.findID(layer, pattern);
    const auto position = dictionary.getFlatCoordinates(cluster, pattern, id);
    const auto& entry = dictionary.getEntry(layer, id);
    local.x = position.X();
    local.z = position.Z();
    local.sigma2X = entry.xWidth * entry.xWidth;
    local.sigma2Z = entry.zWidth * entry.zWidth;
    result.push_back(local);
  }
  if (offset != patterns.size()) {
    throw std::runtime_error(std::format("{} bytes of TRK cluster patterns left over", patterns.size() - offset));
  }
  return result;
}

} // namespace

std::vector<LocalCluster> trkLocalClusters(std::span<const o2::trkft3::TRKCluster> clusters, std::span<const unsigned char> patterns, const TrkDictionary& dictionary) { return alice3LocalClusters(clusters, patterns, dictionary, 0U, {0, 1}); }

std::vector<LocalCluster> ft3LocalClusters(std::span<const o2::trkft3::FT3Cluster> clusters, std::span<const unsigned char> patterns, const Ft3Dictionary& dictionary) { return alice3LocalClusters(clusters, patterns, dictionary, static_cast<std::uint32_t>(detectors::TRKFT3::kFirstDiscSensor), {2}); }

std::vector<Cluster> toClusters(std::span<const LocalCluster> local, std::span<const o2::trkft3::ROFRecord> rofs, ReadoutFrames frames, std::uint32_t tfFirstOrbit)
{
  // The TRK's sensor volume is its active matrix: our origin is O2's (TRK.cxx's placeChips).
  return detail::toClusters(local, rofs, frames, tfFirstOrbit, 0.);
}

ReadoutFrames readoutFrames(std::span<const o2::trkft3::ROFRecord> rofs) { return detail::readoutFrames(rofs); }

std::vector<double> measurementDepthsCM(const TrkDictionary& trk, const Ft3Dictionary& ft3)
{
  std::vector<double> depths;
  depths.reserve(static_cast<std::size_t>(trk.getNumberOfLayers()) + static_cast<std::size_t>(ft3.getNumberOfLayers()));
  for (int layer = 0; layer < trk.getNumberOfLayers(); ++layer) {
    depths.push_back(-static_cast<double>(trk.getLayer(layer).depth));
  }
  for (int layer = 0; layer < ft3.getNumberOfLayers(); ++layer) {
    depths.push_back(-static_cast<double>(ft3.getLayer(layer).depth));
  }
  return depths;
}

std::vector<FrameTiming> trkFrameTiming()
{
  const auto& param = o2::trk::AlmiraParam::Instance();
  constexpr int kLayers = Detector<detectors::TRKFT3>::nLayers;
  // Our layers are the parameter's entries: TRK's barrel layers, then FT3's
  // discs, which take entries 11-22.
  static_assert(kLayers == static_cast<int>(o2::trk::AlmiraParam::kNLayers));
  std::vector<FrameTiming> frames;
  frames.reserve(kLayers);
  for (int layer = 0; layer < kLayers; ++layer) {
    frames.push_back({.lengthBC = param.getROFLengthInBC(layer), .firstStartBC = param.getROFBiasInBC(layer) + param.getROFDelayInBC(layer)});
  }
  return frames;
}

double alice3BzKiloGauss()
{
  const char* macro = std::getenv("ALICE3_MAGFIELD_MACRO");
  if (macro == nullptr || !std::filesystem::exists(macro)) {
    throw std::runtime_error(std::string{"ALICE3_MAGFIELD_MACRO names no file ("} + (macro == nullptr ? "unset" : macro) + "): the ALICE 3 field is unknown, and O2 would take a toy one");
  }
  o2::field::ALICE3MagneticField field;
  const double origin[3] = {0., 0., 0.};
  double b[3] = {0., 0., 0.};
  field.Field(origin, b);
  return b[2];
}

} // namespace o2::gtrk::alice

namespace o2::gtrk::workflow
{

using namespace o2::framework;

namespace
{

/// ALICE 3's TRK and FT3 clusters, one stream per barrel layer and one per
/// disc as O2's cluster readers publish them (TRK/COMPCLUSTERS/<layer>,
/// FT3/COMPCLUSTERS/<disc> and the rest), turned into GTRK/TRKCLUSTERS and
/// GTRK/TRKCLUSMCTR: each layer's
/// stream converted with its own readout frames, each cluster with its
/// position and errors from the topology dictionaries, CCDB conditions,
/// then the layers
/// concatenated in order, barrel first, so a measurement's cluster index is
/// its position in that concatenation and the labels follow it.
class TrkClusterConverterTask : public Task
{
 public:
  static constexpr int kLayers = Detector<detectors::TRKFT3>::nLayers;
  static constexpr int kBarrelLayers = Detector<detectors::TRKFT3>::nBarrelLayers;

  /// Our layer \a layer's stream: TRK's barrel layer, or FT3's disc after them.
  [[nodiscard]] static o2::header::DataOrigin originOf(int layer) { return layer < kBarrelLayers ? o2::header::DataOrigin{"TRK"} : o2::header::DataOrigin{"FT3"}; }
  [[nodiscard]] static o2::header::DataHeader::SubSpecificationType subspecOf(int layer) { return static_cast<o2::header::DataHeader::SubSpecificationType>(layer < kBarrelLayers ? layer : layer - kBarrelLayers); }

  TrkClusterConverterTask(bool useMC, TrkDictionaryConditions dictionaries) : mUseMC(useMC), mDictionaries(dictionaries) {}

  void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj) final { mDictionaries.finaliseCCDB(matcher, obj); }

  void run(ProcessingContext& pc) final
  {
    TrkDictionaryConditions::checkUpdates(pc);
    const auto& timing = pc.services().get<o2::framework::TimingInfo>();
    std::vector<Cluster> clusters;
    o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
    for (int layer = 0; layer < kLayers; ++layer) {
      const auto patterns = pc.inputs().get<gsl::span<unsigned char>>(std::format("patterns_{}", layer));
      const auto rofs = pc.inputs().get<gsl::span<o2::trkft3::ROFRecord>>(std::format("rofs_{}", layer));
      if (rofs.empty()) {
        continue;
      }
      // The timeframe's first orbit as DPL announces it, held to the one the
      // frames imply, as for the ITS.
      const auto implied = o2::raw::HBFUtils::Instance().getFirstIRofTF(rofs.front().getBCData()).orbit;
      if (implied != timing.firstTForbit) {
        LOGP(fatal, "TRKFT3 layer {}: the readout frames belong to the timeframe starting at orbit {}, DPL says {}", layer, implied, timing.firstTForbit);
      }
      const std::span<const unsigned char> patternBytes{patterns.data(), patterns.size()};
      const auto local = [&]() {
        if (layer < kBarrelLayers) {
          const auto compact = pc.inputs().get<gsl::span<o2::trkft3::TRKCluster>>(std::format("compclusters_{}", layer));
          return alice::trkLocalClusters(std::span{compact.data(), compact.size()}, patternBytes, *mDictionaries.trk());
        }
        const auto compact = pc.inputs().get<gsl::span<o2::trkft3::FT3Cluster>>(std::format("compclusters_{}", layer));
        return alice::ft3LocalClusters(std::span{compact.data(), compact.size()}, patternBytes, *mDictionaries.ft3());
      }();
      const std::span<const o2::trkft3::ROFRecord> frames{rofs.data(), rofs.size()};
      const auto converted = alice::toClusters(local, frames, alice::readoutFrames(frames), timing.firstTForbit);
      clusters.insert(clusters.end(), converted.begin(), converted.end());
      if (mUseMC) {
        const auto layerLabels = pc.inputs().get<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>(std::format("labels_{}", layer));
        if (layerLabels->getIndexedSize() != converted.size()) {
          LOGP(fatal, "TRKFT3 layer {}: {} labelled clusters for {} clusters", layer, layerLabels->getIndexedSize(), converted.size());
        }
        labels.mergeAtBack(*layerLabels);
      }
    }
    LOGP(info, "TRKFT3: {} clusters on {} layers for the timeframe at orbit {}", clusters.size(), kLayers, timing.firstTForbit);
    pc.outputs().snapshot(Output{gDataOrigin, description<detectors::TRKFT3>(kClusters), 0}, clusters);
    if (mUseMC) {
      auto& flat = pc.outputs().make<o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel>>(Output{gDataOrigin, description<detectors::TRKFT3>(kClustersMCTruth), 0});
      labels.flatten_to(flat);
    }
  }

 private:
  bool mUseMC;
  TrkDictionaryConditions mDictionaries;
};

} // namespace

template <>
DataProcessorSpec getClusterConverterSpec<detectors::TRKFT3>(bool useMC)
{
  std::vector<InputSpec> inputs;
  std::vector<OutputSpec> outputs{OutputSpec{gDataOrigin, description<detectors::TRKFT3>(kClusters), 0, Lifetime::Timeframe}};
  for (int layer = 0; layer < TrkClusterConverterTask::kLayers; ++layer) {
    const auto origin = TrkClusterConverterTask::originOf(layer);
    const auto subspec = TrkClusterConverterTask::subspecOf(layer);
    inputs.emplace_back(std::format("compclusters_{}", layer), origin, "COMPCLUSTERS", subspec, Lifetime::Timeframe);
    inputs.emplace_back(std::format("patterns_{}", layer), origin, "PATTERNS", subspec, Lifetime::Timeframe);
    inputs.emplace_back(std::format("rofs_{}", layer), origin, "CLUSTERSROF", subspec, Lifetime::Timeframe);
    if (useMC) {
      inputs.emplace_back(std::format("labels_{}", layer), origin, "CLUSTERSMCTR", subspec, Lifetime::Timeframe);
    }
  }
  if (useMC) {
    outputs.emplace_back(gDataOrigin, description<detectors::TRKFT3>(kClustersMCTruth), 0, Lifetime::Timeframe);
  }
  const TrkDictionaryConditions dictionaries(inputs);
  return DataProcessorSpec{
    .name = deviceName<detectors::TRKFT3>("cluster-converter"),
    .inputs = inputs,
    .outputs = outputs,
    .algorithm = AlgorithmSpec{adaptFromTask<TrkClusterConverterTask>(useMC, dictionaries)},
    .options = Options{},
  };
}

template DataProcessorSpec getMeasurementProducerSpec<detectors::TRKFT3>();

} // namespace o2::gtrk::workflow
