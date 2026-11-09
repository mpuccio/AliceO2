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

/// \file clusterClosure.cxx
/// \brief Convert one timeframe of O2's ITS or MFT clusters and hold the
///        result to O2.
///
/// Reads a cluster file, the topology dictionary, the geometry and the
/// digitiser's configuration (for HBFUtils; the readout frames' length is read
/// off the frames themselves), converts the clusters as the converter device
/// will (GenericTracking/Alice/Clusters.h), and checks, cluster by cluster:
///
/// - the global position through our double geometry against O2's L2G matrix
///   applied to the same local position in double, and through our float
///   tracking surfaces against the float position O2's trackers load
///   (`getMatrixL2G(chip) * locXYZ`);
/// - the layer against `GeometryTGeo::getLayer`;
/// - that the position is inside the active matrix, which is what the
///   in-plane origin shift is for;
/// - that the time is inside the timeframe, and the labels as many as the
///   clusters;
/// - that the pattern's spans hold its pixels: n pixels over r rows and c
///   columns, r and c at least 1 and n <= r c (not max(r, c) <= n: O2's
///   clusterer makes a few clusters with rows between their pixels).
///
/// It prints the counts per layer, the readout frames and the error spread,
/// and exits non-zero if a check fails.
///
/// For TRKFT3 (-d TRKFT3, --cluster-dictionaries) every barrel layer's and
/// disc's stream is converted at the dictionaries' measurement depths and held
/// to O2 the same way, with the frames against the digitiser's configuration;
/// then each signal cluster against its particle's MC hit on the chip: the
/// residuals at its layer's measurement depth by sub-detector, axis and span,
/// with their tails (the clusters without such a hit are counted, not failed).
///
/// The closure of GTRK's converter and geometry against O2 on real files:
/// run it whenever the geometry, the placement or the converter changes.

#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Study/Alice/Hits.h"
#include "GenericTracking/Study/Statistics.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Alice/TRK.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tools/SensorFrame.h"

#include "CommonConstants/LHCConstants.h"
#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsTRKFT3/Cluster.h"
#include "DataFormatsTRKFT3/ROFRecord.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsRaw/HBFUtils.h"
#include "ITSBase/GeometryTGeo.h"
#include "MFTBase/GeometryTGeo.h"
#include "FT3Base/GeometryTGeo.h"
#include "TRKBase/GeometryTGeo.h"
#include "TRKBase/SegmentationChip.h"
#include "TRKReconstruction/TopologyDictionary.h"
#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include <TFile.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <format>
#include <functional>
#include <map>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace bpo = boost::program_options;

namespace o2::gtrk::study
{
namespace
{

constexpr double kDoubleTolerance = 1e-6; // [cm]
constexpr double kFloatTolerance = 1e-4;  // [cm]

/// The centre of \a cluster's time box [ns from the timeframe's start]: its
/// readout frame's middle, which the timeframe must hold.
[[nodiscard]] double boxCentreNS(const Cluster& cluster) { return 1e-3 * (static_cast<double>(cluster.timePS) + (0.5 * static_cast<double>(cluster.widthPS))); }

struct Options {
  std::string detector;
  std::string clusters;
  std::string dictionary;
  std::string geometry;
  std::string config;
  std::string hits;
  std::string ft3Clusters;
  std::string ft3Hits;
  std::string clusterDictionaries;
};

/// One layer's stream of an ALICE 3 cluster file (TRK's or FT3's).
template <typename ClusterType>
struct Alice3Stream {
  std::vector<ClusterType> compact;
  std::vector<unsigned char> patterns;
  std::vector<o2::trkft3::ROFRecord> rofs;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
};

/// Stream \a stream of \a tree, whose branches start with \a prefix (TRK, FT3).
template <typename ClusterType>
Alice3Stream<ClusterType> readStream(TTree& tree, const std::string& prefix, int stream)
{
  Alice3Stream<ClusterType> result;
  auto* compactPtr = &result.compact;
  auto* patternsPtr = &result.patterns;
  auto* rofsPtr = &result.rofs;
  auto* labelsPtr = &result.labels;
  tree.SetBranchAddress(std::format("{}ClusterComp_{}", prefix, stream).c_str(), &compactPtr);
  tree.SetBranchAddress(std::format("{}ClusterPatt_{}", prefix, stream).c_str(), &patternsPtr);
  tree.SetBranchAddress(std::format("{}ClustersROF_{}", prefix, stream).c_str(), &rofsPtr);
  tree.SetBranchAddress(std::format("{}ClusterMCTruth_{}", prefix, stream).c_str(), &labelsPtr);
  tree.GetEntry(0);
  tree.ResetBranchAddresses();
  return result;
}

/// The o2sim tree of the cluster file \a path, held open by \a file.
TTree& clusterTree(const std::string& path, std::unique_ptr<TFile>& file)
{
  file.reset(TFile::Open(path.c_str()));
  if (!file || file->IsZombie()) {
    throw std::runtime_error("cannot open " + path);
  }
  auto* tree = file->Get<TTree>("o2sim");
  if (tree == nullptr || tree->GetEntries() < 1) {
    throw std::runtime_error(path + ": no o2sim tree, or an empty one");
  }
  return *tree;
}

/// One timeframe of one detector's cluster file.
struct ClusterFile {
  std::vector<o2::itsmft::CompClusterExt> clusters;
  std::vector<unsigned char> patterns;
  std::vector<o2::itsmft::ROFRecord> rofs;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
  bool hasLabels = false;
};

ClusterFile readClusters(const std::string& path, const std::string& prefix)
{
  std::unique_ptr<TFile> file;
  auto* tree = &clusterTree(path, file);
  ClusterFile result;
  auto* clusters = &result.clusters;
  auto* patterns = &result.patterns;
  auto* rofs = &result.rofs;
  auto* labels = &result.labels;
  tree->SetBranchAddress((prefix + "ClusterComp").c_str(), &clusters);
  tree->SetBranchAddress((prefix + "ClusterPatt").c_str(), &patterns);
  tree->SetBranchAddress((prefix + "ClustersROF").c_str(), &rofs);
  result.hasLabels = tree->GetBranch((prefix + "ClusterMCTruth").c_str()) != nullptr;
  if (result.hasLabels) {
    tree->SetBranchAddress((prefix + "ClusterMCTruth").c_str(), &labels);
  }
  tree->GetEntry(0);
  return result;
}

struct Deviation {
  double worst = 0.;
  std::size_t at = 0;

  void add(double value, std::size_t index)
  {
    if (value > worst) {
      worst = value;
      at = index;
    }
  }
};

template <DetectorLayout Det, typename GeometryTGeo>
int runFor(const Options& options, GeometryTGeo& o2Geometry, const std::string& prefix)
{
  const auto sensors = alice::placeChecked<Det>();
  const auto geometry = std::make_unique<const Geometry<Det>>(sensors);
  const auto tracking = std::make_unique<const TrackingGeometry<Det>>(*geometry);
  o2Geometry.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));

  const std::unique_ptr<o2::itsmft::TopologyDictionary> dictionary{o2::itsmft::TopologyDictionary::loadFrom(options.dictionary)};
  if (!dictionary) {
    throw std::runtime_error("cannot read the topology dictionary " + options.dictionary);
  }
  const auto file = readClusters(options.clusters, prefix);
  const auto local = alice::localClusters(file.clusters, file.patterns, *dictionary);
  const auto frames = alice::readoutFrames(file.rofs);
  const auto& hbf = o2::raw::HBFUtils::Instance();
  const auto tfFirstOrbit = hbf.getFirstSampledTFIR().orbit;
  const auto clusters = alice::toClusters(local, file.rofs, frames, tfFirstOrbit);
  const double durationNS = static_cast<double>(hbf.getNOrbitsPerTF()) * o2::constants::lhc::LHCOrbitNS;

  Deviation exact;
  Deviation narrow;
  std::size_t wrongLayer = 0;
  std::size_t outside = 0;
  std::size_t outsideTime = 0;
  std::size_t badSpans = 0;
  std::vector<std::size_t> perLayer(static_cast<std::size_t>(Detector<Det>::nLayers), 0U);
  std::vector<float> rowSpans;
  std::vector<float> colSpans;
  rowSpans.reserve(clusters.size());
  colSpans.reserve(clusters.size());
  std::vector<float> sigmaX;
  std::vector<float> sigmaY;
  sigmaX.reserve(clusters.size());
  sigmaY.reserve(clusters.size());
  for (std::size_t i = 0; i < clusters.size(); ++i) {
    const auto& cluster = clusters[i];
    const auto chip = static_cast<int>(cluster.chipID);
    const auto& l2g = o2Geometry.getMatrixL2G(chip);

    // Double: O2's matrix and our sensor, on the same local position.
    const auto o2Double = l2g(o2::math_utils::Point3D<double>{local[i].x, 0., local[i].z});
    const double ours[3] = {static_cast<double>(cluster.xLocal), static_cast<double>(cluster.yLocal), 0.};
    double global[3];
    geometry->toGlobal(cluster.chipID, ours, global);
    exact.add(std::max({std::abs(global[0] - o2Double.X()), std::abs(global[1] - o2Double.Y()), std::abs(global[2] - o2Double.Z())}), i);

    // Float: what O2's trackers load, and what our tracking surfaces give.
    const auto o2Float = l2g * o2::math_utils::Point3D<float>{local[i].x, 0.F, local[i].z};
    const auto& surface = tracking->surface(cluster.chipID);
    Scalar point[3];
    surfaceToGlobal(surface, cluster.xLocal, cluster.yLocal, point);
    narrow.add(std::max({std::abs(static_cast<double>(point[0] - o2Float.X())), std::abs(static_cast<double>(point[1] - o2Float.Y())), std::abs(static_cast<double>(point[2] - o2Float.Z()))}), i);

    const auto layer = geometry->sensor(cluster.chipID).layer;
    wrongLayer += static_cast<int>(layer) != o2Geometry.getLayer(chip) ? 1U : 0U;
    // The active area: the float tracking surface's rectangular bounds on
    // the unrolled coordinates, the ones the tracking reads.
    outside += surface.bounds.isInside(cluster.xLocal, cluster.yLocal) ? 0U : 1U;
    outsideTime += boxCentreNS(cluster) < durationNS ? 0U : 1U;
    const auto rows = static_cast<std::uint32_t>(cluster.rowSpan);
    const auto columns = static_cast<std::uint32_t>(cluster.colSpan);
    if (rows == 0U || columns == 0U || local[i].nPixels > rows * columns) {
      if (badSpans < 10U) {
        const auto id = file.clusters[i].getPatternID();
        std::cout << std::format("  spans {} x {} for {} pixels: cluster {}, pattern ID {}{}\n", rows, columns, local[i].nPixels, i, id, id != o2::itsmft::CompCluster::InvalidPatternID && dictionary->isGroup(id) ? ", a group" : "");
      }
      ++badSpans;
    }
    rowSpans.push_back(static_cast<float>(rows));
    colSpans.push_back(static_cast<float>(columns));
    ++perLayer[layer];
    sigmaX.push_back(cluster.sigmaXLocal);
    sigmaY.push_back(cluster.sigmaYLocal);
  }

  std::size_t filled = 0;
  for (const auto& rof : file.rofs) {
    filled += rof.getNEntries() > 0 ? 1U : 0U;
  }
  // The error spread as quantiles: a few huge clusters own the mean.
  const auto quantiles = [](std::vector<float> values) {
    const auto at = [&values](double q) { return values.empty() ? 0. : 1e4 * static_cast<double>(study::quantile(std::span{values}, q)); };
    return std::format("{:.2f} / {:.2f} / {:.2f} / {:.2f}", at(0.5), at(0.9), at(0.99), at(1.));
  };
  const auto spanQuantiles = [](std::vector<float> values) {
    const auto at = [&values](double q) { return values.empty() ? 0.F : study::quantile(std::span{values}, q); };
    return std::format("{} / {} / {} / {}", at(0.5), at(0.9), at(0.99), at(1.));
  };
  std::cout << std::format("{}: {} clusters in {} readout frames of {} BC ({} with clusters), first TF orbit {}\n", Det::name, clusters.size(), file.rofs.size(), frames.lengthBC, filled, tfFirstOrbit) << std::format("  per layer:");
  for (std::size_t layer = 0; layer < perLayer.size(); ++layer) {
    std::cout << std::format(" {}:{}", layer, perLayer[layer]);
  }
  std::cout << std::format("\n  sigma across, median / 90 % / 99 % / max: {} um\n", quantiles(sigmaX)) << std::format("  sigma along,  median / 90 % / 99 % / max: {} um\n", quantiles(sigmaY)) << std::format("  double position: worst {:.3g} cm (cluster {}), tolerance {:.0e}\n", exact.worst, exact.at, kDoubleTolerance) << std::format("  float position:  worst {:.3g} cm (cluster {}), tolerance {:.0e}\n", narrow.worst, narrow.at, kFloatTolerance) << std::format("  wrong layer {}, outside the active matrix {}, outside the timeframe {}\n", wrongLayer, outside, outsideTime) << std::format("  rows, median / 90 % / 99 % / max: {}; columns: {}; spans that cannot hold their pixels {}\n", spanQuantiles(rowSpans), spanQuantiles(colSpans), badSpans);

  // Which topologies the errors come from: the most used, with the
  // dictionary's errors and whether each is a group.
  std::map<int, std::size_t> use;
  std::size_t invalid = 0;
  std::size_t groups = 0;
  for (const auto& cluster : file.clusters) {
    const auto id = cluster.getPatternID();
    if (id == o2::itsmft::CompCluster::InvalidPatternID) {
      ++invalid;
      continue;
    }
    groups += dictionary->isGroup(id) ? 1U : 0U;
    ++use[id];
  }
  std::vector<std::pair<std::size_t, int>> top;
  top.reserve(use.size());
  for (const auto& [id, count] : use) {
    top.emplace_back(count, id);
  }
  std::ranges::sort(top, std::greater{});
  std::cout << std::format("  topologies: {} used, {} clusters in groups, {} with an invalid ID; the most used:\n", use.size(), groups, invalid);
  for (std::size_t i = 0; i < std::min<std::size_t>(top.size(), 6U); ++i) {
    const auto id = top[i].second;
    std::cout << std::format("    id {:5} x{:7}: {:2} pixels, error {:6.2f} x {:6.2f} um{}\n", id, top[i].first, dictionary->getNpixels(id), 1e4 * static_cast<double>(dictionary->getErrX(id)), 1e4 * static_cast<double>(dictionary->getErrZ(id)), dictionary->isGroup(id) ? ", group" : "");
  }

  bool ok = exact.worst <= kDoubleTolerance && narrow.worst <= kFloatTolerance && wrongLayer == 0U && outside == 0U && outsideTime == 0U && badSpans == 0U;
  if (file.hasLabels) {
    const bool labelled = file.labels.getIndexedSize() == clusters.size();
    std::cout << std::format("  MC labels for {} of {} clusters\n", file.labels.getIndexedSize(), clusters.size());
    ok = ok && labelled;
  }
  std::cout << (ok ? "closure: pass\n" : "closure: FAIL\n");
  return ok ? 0 : 1;
}

/// The spread of residuals [cm]: their count, mean, rms, the core's width,
/// half the 15.9 - 84.1 % range, which a few delta electrons do not move, and
/// the fractions further than 3 and 5 core widths from the median: the tails
/// a fit taking the core as its error meets.
struct Spread {
  std::vector<double> values;

  [[nodiscard]] std::string describe() const
  {
    if (values.empty()) {
      return "      -";
    }
    auto sorted = values;
    std::ranges::sort(sorted);
    double sum = 0.;
    double sum2 = 0.;
    for (const double value : sorted) {
      sum += value;
      sum2 += value * value;
    }
    const auto n = static_cast<double>(sorted.size());
    const double mean = sum / n;
    const auto at = [&sorted](double q) { return study::quantile(std::span{sorted}, q); };
    const double median = at(0.5);
    const double width = 0.5 * (at(0.841) - at(0.159));
    const auto beyond = [&](double k) { return 100. * static_cast<double>(std::ranges::count_if(sorted, [&](double value) { return std::abs(value - median) > k * width; })) / n; };
    // The normalised median absolute deviation beside it: the width the pulls measure (TrackPulls).
    std::vector<double> deviations(sorted.size());
    std::ranges::transform(sorted, deviations.begin(), [&](double value) { return std::abs(value - median); });
    const double mad = 1.4826 * study::median(deviations);
    return std::format("{:8} {:6.2f} {:6.2f} {:6.2f} {:6.2f} | >3 {:5.2f} % >5 {:5.2f} %", sorted.size(), 1e4 * mean, 1e4 * std::sqrt(std::max(sum2 / n - (mean * mean), 0.)), 1e4 * width, 1e4 * mad, beyond(3.), beyond(5.));
  }
};

/// The TRK and FT3: the stream of
/// every barrel layer and every disc converted on its own, held to O2's centre
/// of gravity and L2G matrix, to O2's layer, the active area, the digitiser's
/// frames and the labels; then each signal cluster against its MC hit, the
/// residuals of its centre of gravity per sub-detector (VD, ML/OT, FT3), axis
/// and span, at its layer's measurement depth. Every cluster is converted as
/// the converter device does, by its topology dictionary's entry
/// (--cluster-dictionaries),
/// and held to O2's position of the same entry at the layer's depth, our
/// geometry placed at the dictionaries' depths; the residuals are the
/// dictionary's, taken at those depths.
int runTRK(const Options& options)
{
  if (options.clusterDictionaries.empty()) {
    throw std::runtime_error("TRKFT3 needs --cluster-dictionaries");
  }
  using Det = detectors::TRKFT3;
  using Segmentation = o2::trk::SegmentationChip;
  constexpr int kLayers = Detector<Det>::nLayers;
  constexpr int kBarrelLayers = Detector<Det>::nBarrelLayers;
  constexpr auto kFirstDisc = static_cast<std::uint32_t>(Det::kFirstDiscSensor);
  auto& o2Geometry = *o2::trk::GeometryTGeo::Instance();
  o2Geometry.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  auto& ft3Geometry = *o2::ft3::GeometryTGeo::Instance();
  ft3Geometry.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  // The dictionaries' CCDB snapshots below the directory, as the devices take them.
  const auto load = [&]<typename Dictionary>(std::unique_ptr<Dictionary>& dictionary, const char* path) { dictionary.reset(Dictionary::loadFrom(options.clusterDictionaries + "/" + path + "/snapshot.root")); };
  std::unique_ptr<alice::TrkDictionary> trkDictionary;
  std::unique_ptr<alice::Ft3Dictionary> ft3Dictionary;
  load(trkDictionary, "TRK/Calib/ClusterDictionary");
  load(ft3Dictionary, "FT3/Calib/ClusterDictionary");
  const auto depths = alice::measurementDepthsCM(*trkDictionary, *ft3Dictionary);
  const auto geometry = std::make_unique<const Geometry<Det>>(alice::placeChecked<Det>(), depths);
  // The digitiser's frames, from its configuration.
  o2::conf::ConfigurableParam::updateFromFile(options.config, "TRKAlmiraParam");
  const auto frames = alice::trkFrameTiming();
  const auto& hbf = o2::raw::HBFUtils::Instance();
  const auto tfFirstOrbit = hbf.getFirstSampledTFIR().orbit;
  const double durationNS = static_cast<double>(hbf.getNOrbitsPerTF()) * o2::constants::lhc::LHCOrbitNS;

  std::unique_ptr<TFile> trkFile;
  std::unique_ptr<TFile> ft3File;
  auto& trkTree = clusterTree(options.clusters, trkFile);
  auto& ft3Tree = clusterTree(options.ft3Clusters, ft3File);
  study::HitFile hits(options.hits, study::Chips::kTRK, options.ft3Hits);
  // Each event's hits by (track, chip), in file order, built on first use:
  // a signal cluster looks up its particle's crossings of its chip without
  // scanning its event (Pb-Pb events are large).
  using HitIndex = std::unordered_map<std::uint64_t, std::vector<const study::HitSegment*>>;
  const auto hitKey = [](int track, std::uint32_t chip) { return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(track)) << 32U) | chip; };
  std::vector<std::optional<HitIndex>> hitIndex(static_cast<std::size_t>(hits.nEvents()));
  const auto segmentsOf = [&](int event, int track, std::uint32_t chip) -> const std::vector<const study::HitSegment*>* {
    if (event < 0 || event >= hits.nEvents()) {
      return nullptr;
    }
    auto& index = hitIndex[static_cast<std::size_t>(event)];
    if (!index) {
      index.emplace();
      for (const auto& segment : hits.event(event)) {
        (*index)[hitKey(segment.track, segment.chip)].push_back(&segment);
      }
    }
    const auto found = index->find(hitKey(track, chip));
    return found == index->end() ? nullptr : &found->second;
  };

  Deviation exact;
  Deviation narrow;
  std::size_t total = 0;
  std::size_t wrongLayer = 0;
  std::size_t outside = 0;
  std::size_t outsideTime = 0;
  std::size_t badFrames = 0;
  std::size_t unlabelled = 0;
  std::size_t signal = 0;
  std::size_t noHit = 0;
  // Every signal cluster with its particle's hit on its chip: the cluster's
  // unrolled position, the hit's two ends unrolled (across, along, depth),
  // its sub-detector and its spans.
  struct Match {
    int layer = 0;
    int subDet = 0;
    int rowSpan = 0;
    int colSpan = 0;
    std::array<double, 2> cluster{};
    std::array<std::array<double, 3>, 2> ends{};
  };
  std::vector<Match> matches;
  std::cout << "TRKFT3 layer | clusters | frames: records' length, configured length@first start\n";
  // One layer's stream: TRK's barrel layer, or FT3's disc after them.
  const auto checkLayer = [&]<typename ClusterType>(int layer, const Alice3Stream<ClusterType>& stream) {
    constexpr bool kDisc = std::is_same_v<ClusterType, o2::trkft3::FT3Cluster>;
    const auto& [compact, patterns, rofs, labels] = stream;
    const auto local = [&]() {
      if constexpr (kDisc) {
        return alice::ft3LocalClusters(compact, patterns, *ft3Dictionary);
      } else {
        return alice::trkLocalClusters(compact, patterns, *trkDictionary);
      }
    }();
    const auto spacing = alice::readoutFrames(std::span<const o2::trkft3::ROFRecord>{rofs});
    badFrames += spacing.lengthBC != frames[static_cast<std::size_t>(layer)].lengthBC ? 1U : 0U;
    const auto clusters = alice::toClusters(local, std::span<const o2::trkft3::ROFRecord>{rofs}, spacing, tfFirstOrbit);
    unlabelled += labels.getIndexedSize() == clusters.size() ? 0U : 1U;
    std::cout << std::format("{:9} | {:8} | {:3} BC, {}@{}\n", layer, clusters.size(), spacing.lengthBC, frames[static_cast<std::size_t>(layer)].lengthBC, frames[static_cast<std::size_t>(layer)].firstStartBC);

    std::size_t offset = 0;
    for (std::size_t i = 0; i < clusters.size(); ++i, ++total) {
      const auto& cluster = clusters[i];
      const auto& o2Cluster = compact[i];
      const auto* pattern = patterns.data() + offset;
      const std::size_t bytes = 2U + ((static_cast<std::size_t>(pattern[0]) * pattern[1] + 7U) / 8U);
      offset += bytes;
      const auto chip = static_cast<int>(kDisc ? cluster.chipID - kFirstDisc : cluster.chipID); // O2's
      const auto& l2g = kDisc ? ft3Geometry.getMatrixL2G(chip) : o2Geometry.getMatrixL2G(chip);

      // O2's position: its dictionary entry's at the layer's depth, curved on
      // the VD, through its matrix.
      const auto o2Local = [&]() {
        if constexpr (kDisc) {
          return ft3Dictionary->getClusterCoordinates(o2Cluster, pattern);
        } else {
          return trkDictionary->getClusterCoordinates(o2Cluster, pattern);
        }
      }();
      const auto o2Double = l2g(o2::math_utils::Point3D<double>{o2Local.X(), o2Local.Y(), o2Local.Z()});
      const auto& spec = geometry->specOf(cluster.chipID);
      double frame[3];
      spec.unrolledToFrame(static_cast<double>(cluster.xLocal), static_cast<double>(cluster.yLocal), geometry->measurementDepth(geometry->sensor(cluster.chipID).layer), frame);
      double global[3];
      geometry->toGlobal(cluster.chipID, frame, global);
      exact.add(std::max({std::abs(global[0] - o2Double.X()), std::abs(global[1] - o2Double.Y()), std::abs(global[2] - o2Double.Z())}), total);
      // Float: what O2's tracker loads, and what our producer ships.
      const auto o2Float = l2g * o2Local;
      const auto point = gridPoint<Det>(*geometry, cluster.chipID, static_cast<double>(cluster.xLocal), static_cast<double>(cluster.yLocal));
      narrow.add(std::max({std::abs(static_cast<double>(point.x - o2Float.X())), std::abs(static_cast<double>(point.y - o2Float.Y())), std::abs(static_cast<double>(point.z - o2Float.Z()))}), total);

      const auto ourLayer = static_cast<int>(geometry->sensor(cluster.chipID).layer);
      const int o2Layer = kDisc ? kBarrelLayers + ft3Geometry.getLayer(chip) : o2Geometry.getLayerTRK(chip);
      wrongLayer += ourLayer != layer || ourLayer != o2Layer ? 1U : 0U;
      // The active area: the double sensor spec's rectangular bounds on the
      // unrolled coordinates (a curved VD sensor's too).
      outside += tools::contains(spec, static_cast<double>(cluster.xLocal), static_cast<double>(cluster.yLocal)) ? 0U : 1U;
      outsideTime += boxCentreNS(cluster) < durationNS ? 0U : 1U;

      // The residual to the MC hit of its first signal label, at the middle
      // of the hit's path through the silicon, in the unrolled frame.
      const o2::MCCompLabel* label = nullptr;
      if (i < labels.getIndexedSize()) {
        for (const auto& candidate : labels.getLabels(static_cast<std::uint32_t>(i))) {
          if (candidate.isValid() && !candidate.isNoise() && candidate.getSourceID() == 0) {
            label = &candidate;
            break;
          }
        }
      }
      if (label == nullptr) {
        continue;
      }
      ++signal;
      // The particle's crossing of the chip: the segment whose middle is
      // nearest the cluster, where several (a looper's).
      double best = std::numeric_limits<double>::max();
      bool found = false;
      Match match{
        .layer = layer,
        .subDet = o2Cluster.subDetID,
        .rowSpan = pattern[0],
        .colSpan = pattern[1],
        .cluster = {static_cast<double>(cluster.xLocal), static_cast<double>(cluster.yLocal)},
      };
      static const std::vector<const study::HitSegment*> kNone;
      const auto* crossings = segmentsOf(label->getEventID(), label->getTrackID(), cluster.chipID);
      for (const auto* segment : crossings != nullptr ? *crossings : kNone) {
        std::array<std::array<double, 3>, 2> ends{};
        for (std::size_t end = 0; end < 2; ++end) {
          double hitFrame[3];
          tools::toLocal(geometry->sensor(cluster.chipID), end == 0 ? segment->start.data() : segment->end.data(), hitFrame);
          tools::frameToUnrolled(spec, hitFrame, ends[end].data());
        }
        const double distance = std::hypot(match.cluster[0] - (0.5 * (ends[0][0] + ends[1][0])), match.cluster[1] - (0.5 * (ends[0][1] + ends[1][1])));
        if (distance < best) {
          best = distance;
          match.ends = ends;
          found = true;
        }
      }
      if (!found) {
        ++noHit;
        continue;
      }
      matches.push_back(match);
    }
  };
  for (int layer = 0; layer < kLayers; ++layer) {
    if (layer < kBarrelLayers) {
      checkLayer(layer, readStream<o2::trkft3::TRKCluster>(trkTree, "TRK", layer));
    } else {
      checkLayer(layer, readStream<o2::trkft3::FT3Cluster>(ft3Tree, "FT3", layer - kBarrelLayers));
    }
  }

  std::cout << std::format(
                 "TRKFT3: {} clusters; double position: worst {:.3g} cm (cluster {}), tolerance {:.0e}; float position: worst {:.3g} cm "
                 "(cluster {}), tolerance {:.0e}\n",
                 total, exact.worst, exact.at, kDoubleTolerance, narrow.worst, narrow.at, kFloatTolerance)
            << std::format(
                 "  wrong layer {}, outside the active area {}, outside the timeframe {}, layers with other frames than configured {}, "
                 "layers whose labels do not match {}\n",
                 wrongLayer, outside, outsideTime, badFrames, unlabelled)
            << std::format("  signal clusters {}, without their particle's hit on the chip {}\n", signal, noHit);

  // The residuals, cluster minus the particle's crossing of the plane at
  // its layer's measurement depth -- the dictionaries', the converter's
  // (the hit's ends interpolated) [um], by sub-detector, axis and span, with
  // the tails a fit taking the core as its error meets. Every signal cluster
  // counts, at the matrix's edge and the pieces of a split crossing too
  // (where O2's CheckTopologyDictionary.C leaves them out).
  constexpr std::size_t kSubDets = 3;
  const std::array<const char*, kSubDets> subDets{"VD", "ML/OT", "FT3"};
  const std::array<std::array<double, 2>, kSubDets> pitches{{{static_cast<double>(Segmentation::PitchRowVD), static_cast<double>(Segmentation::PitchColVD)}, {static_cast<double>(Segmentation::PitchRowMLOT), static_cast<double>(Segmentation::PitchColMLOT)}, {static_cast<double>(Segmentation::PitchRowMLOT), static_cast<double>(Segmentation::PitchColMLOT)}}};
  const auto residualAt = [](const Match& match, double depth) {
    const auto& [a, b] = match.ends;
    // On the segment: a particle that starts or stops in the silicon
    // crosses no deeper plane than its own end.
    const double t = std::abs(b[2] - a[2]) > 1e-9 ? std::clamp((depth - a[2]) / (b[2] - a[2]), 0., 1.) : 0.5;
    return std::array<double, 2>{match.cluster[0] - (a[0] + (t * (b[0] - a[0]))), match.cluster[1] - (a[1] + (t * (b[1] - a[1])))};
  };
  std::cout << "  each layer's measurement depth, the dictionaries' (our unrolled z) [um]:";
  for (int layer = 0; layer < kLayers; ++layer) {
    std::cout << std::format(" {}:{:.1f}", layer, 1e4 * geometry->measurementDepth(static_cast<LayerID>(layer)));
  }
  std::cout << '\n';

  const std::array<const char*, 2> axes{"across (rows)", "along (columns)"};
  constexpr std::size_t kSpans = 7; // 1-5, 6+; 0 every span
  std::array<std::array<std::array<Spread, kSpans>, 2>, kSubDets> residuals{};
  for (const auto& match : matches) {
    const auto subDet = static_cast<std::size_t>(match.subDet);
    const auto residual = residualAt(match, geometry->measurementDepth(static_cast<LayerID>(match.layer)));
    const std::array<int, 2> spans{match.rowSpan, match.colSpan};
    for (std::size_t axis = 0; axis < 2; ++axis) {
      residuals[subDet][axis][0].values.push_back(residual[axis]);
      residuals[subDet][axis][static_cast<std::size_t>(std::min(spans[axis], 6))].values.push_back(residual[axis]);
    }
  }
  std::cout << "  residuals at each layer's measurement depth [um]: clusters, mean, rms, core (half the 68 % range), 1.4826 MAD\n";
  for (std::size_t subDet = 0; subDet < kSubDets; ++subDet) {
    for (std::size_t axis = 0; axis < 2; ++axis) {
      const double pitch = 1e4 * pitches[subDet][axis];
      std::cout << std::format("  {} {}: pitch {:.1f}, pitch/sqrt(12) {:.2f}, pitch/2 {:.2f}\n", subDets[subDet], axes[axis], pitch, pitch / std::sqrt(12.), 0.5 * pitch) << std::format("    every span {}\n", residuals[subDet][axis][0].describe());
      for (std::size_t span = 1; span < kSpans; ++span) {
        std::cout << std::format("    span {}{} {}\n", span, span + 1 == kSpans ? "+" : " ", residuals[subDet][axis][span].describe());
      }
    }
  }

  const bool ok = exact.worst <= kDoubleTolerance && narrow.worst <= kFloatTolerance && wrongLayer == 0U && outside == 0U && outsideTime == 0U && badFrames == 0U && unlabelled == 0U;
  std::cout << (ok ? "closure: pass\n" : "closure: FAIL\n");
  return ok ? 0 : 1;
}

int run(int argc, char** argv)
{
  Options options;
  bpo::options_description description("Convert one timeframe of O2's ITS, MFT or TRKFT3 clusters and hold the result to O2");
  description.add_options()                                                                                                                                   //
    ("help,h", "this message")                                                                                                                                //
    ("detector,d", bpo::value(&options.detector)->default_value("ITS"), "ITS, MFT or TRKFT3")                                                                 //
    ("clusters,c", bpo::value(&options.clusters)->default_value(""), "cluster file; o2clus_its.root, mftclusters.root or o2clus_trk.root if empty")           //
    ("dictionary", bpo::value(&options.dictionary)->default_value(""), "ITS, MFT: topology dictionary (a CCDB snapshot)")                                     //
    ("hits", bpo::value(&options.hits)->default_value("sgn_HitsTRK.root"), "TRKFT3: the MC hits the residuals are taken to")                                  //
    ("ft3-clusters", bpo::value(&options.ft3Clusters)->default_value("o2clus_ft3.root"), "TRKFT3: FT3's cluster file, the discs'")                            //
    ("ft3-hits", bpo::value(&options.ft3Hits)->default_value("sgn_HitsFT3.root"), "TRKFT3: FT3's MC hits")                                                    //
    ("cluster-dictionaries", bpo::value(&options.clusterDictionaries)->default_value(""), "TRKFT3: the topology dictionaries' snapshot directory (required)") //
    ("geometry,g", bpo::value(&options.geometry)->default_value("o2sim_geometry-aligned.root"), "o2sim geometry")                                             //
    ("config", bpo::value(&options.config)->default_value("o2simdigitizerworkflow_configuration.ini"), "the digitiser's configuration, for HBFUtils (and the TRK's frames)");
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

  o2::conf::ConfigurableParam::updateFromFile(options.config, "HBFUtils");
  o2::base::GeometryManager::loadGeometry(options.geometry);
  if (options.detector == "TRKFT3") {
    if (options.clusters.empty()) {
      options.clusters = "o2clus_trk.root";
    }
    return runTRK(options);
  }
  if (options.dictionary.empty()) {
    throw std::runtime_error("the ITS and the MFT need --dictionary");
  }
  if (options.detector == "ITS") {
    if (options.clusters.empty()) {
      options.clusters = "o2clus_its.root";
    }
    return runFor<detectors::ITS>(options, *o2::its::GeometryTGeo::Instance(), "ITS");
  }
  if (options.detector == "MFT") {
    if (options.clusters.empty()) {
      options.clusters = "mftclusters.root";
    }
    return runFor<detectors::MFT>(options, *o2::mft::GeometryTGeo::Instance(), "MFT");
  }
  std::cerr << "unknown detector '" << options.detector << "'\n";
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
