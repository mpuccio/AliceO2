// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// SPDX-License-Identifier: GPL-3.0-only

// Match real reconstructed GTRK states to native TOF clusters and use MC labels
// only for evaluation. An optional transport-hit proxy mode is kept for diagnostics.
#include "GenericTracking/TOF/Clusters.h"
#include "GenericTracking/TOF/Transport.h"
#include "GenericTracking/Study/TrackingInput.h"
#include "GenericTracking/Workflow/TrackingInterface.h"
#include "GenericTracking/Alice/TRK.h"
#include "DataFormatsIOTOF/Cluster.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsRaw/HBFUtils.h"
#include "ITSMFTSimulation/Hit.h"
#include "IOTOFBase/GeometryTGeo.h"
#include "IOTOFBase/IOTOFBaseParam.h"
#include "SimulationDataFormat/DigitizationContext.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include <TList.h>
#include <TStreamerInfo.h>
#include <TFile.h>
#include <TTree.h>
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <unordered_map>

namespace
{
using namespace o2::gtrk;
using Det = detectors::TRKFT3;
using JSON = nlohmann::ordered_json;
struct TGeoMaterial {
  o2::base::MatBudget getMatBudget(float x0, float y0, float z0, float x1, float y1, float z1) const
  {
    return o2::base::GeometryManager::meanMaterialBudget(x0, y0, z0, x1, y1, z1);
  }
};
struct Selected {
  tof::Track track;
  std::uint64_t particle;
  study::ParticleTruth truth;
};
struct Counts {
  std::size_t eligible = 0, matched = 0, correct = 0, wrong = 0, withoutTruthHit = 0, spatial = 0, coarseTime = 0, ambiguous = 0;
  JSON json() const
  {
    return {{"eligible", eligible}, {"ambiguous_matches", ambiguous}, {"matched_eligible", matched}, {"correct", correct}, {"wrong_eligible", wrong}, {"matched_without_truth_hit", withoutTruthHit}, {"with_spatial_truth_candidate", spatial}, {"with_coarse_time_truth_candidate", coarseTime}, {"total_matches", matched + withoutTruthHit}, {"purity", matched + withoutTruthHit ? JSON(static_cast<double>(correct) / (matched + withoutTruthHit)) : JSON(nullptr)}, {"correct_efficiency", eligible ? JSON(static_cast<double>(correct) / eligible) : JSON(nullptr)}};
  }
};

int run(int argc, char** argv)
{
  namespace bpo = boost::program_options;
  std::string inputFile, geometryFile, materialFile, hitFile, contextFile, clusterFile, output;
  int maxTracks = 20, threads = 5;
  unsigned seed = 42, firstOrbit = 0;
  double minPt = 0.5, maxEta = 1., preselectionNS = 100.;
  bpo::options_description options("TOF matching study with native clusters or a transport-hit proxy");
  options.add_options()("help,h", "show options")("input", bpo::value(&inputFile)->required(), "frozen GTRK tracking input, with .ini and .truth")("geometry", bpo::value(&geometryFile)->required(), "same aligned ALICE3 geometry")("tracking-material", bpo::value(&materialFile)->required(), "frozen tracker LUT; matching uses TGeo")("tof-clusters", bpo::value(&clusterFile), "native TF3 clusters with MC labels; replaces hit model")("tof-hits", bpo::value(&hitFile), "source 0 TF3 transport-hit file")("context", bpo::value(&contextFile), "same timeframe's collision context")("output", bpo::value(&output)->required(), "output prefix for JSON and per-track CSV")("max-tracks", bpo::value(&maxTracks)->default_value(20), "fixed-seed sample of unique good primary tracks")("min-pt", bpo::value(&minPt)->default_value(0.5), "minimum generated pT [GeV]")("max-abs-eta", bpo::value(&maxEta)->default_value(1.), "maximum absolute generated eta")("seed", bpo::value(&seed)->default_value(42), "measurement smearing and track sampling seed")("first-orbit", bpo::value(&firstOrbit)->default_value(0), "timeframe anchor orbit")("threads", bpo::value(&threads)->default_value(5), "tracking CPU threads")("preselection-ns", bpo::value(&preselectionNS)->default_value(100.), "hit-time margin around the track readout interval; truth losses counted");
  bpo::variables_map values;
  bpo::store(bpo::parse_command_line(argc, argv, options), values);
  if (values.contains("help")) {
    std::cout << options << '\n';
    return 0;
  }
  bpo::notify(values);
  if (clusterFile.empty() && (hitFile.empty() || contextFile.empty())) {
    throw std::invalid_argument("provide --tof-clusters or both --tof-hits and --context");
  }
  if (!clusterFile.empty() && (!hitFile.empty() || !contextFile.empty())) {
    throw std::invalid_argument("cluster and transport-hit inputs are mutually exclusive");
  }
  if (maxTracks <= 0 || minPt < 0. || maxEta <= 0. || preselectionNS < 0.) {
    throw std::invalid_argument("invalid sample size or kinematic/time selection");
  }
  if (!std::filesystem::path(output).parent_path().empty()) {
    std::filesystem::create_directories(std::filesystem::path(output).parent_path());
  }
  o2::conf::ConfigurableParam::updateFromFile(inputFile + ".ini", "GTRKRecoTRK,HBFUtils");
  o2::conf::ConfigurableParam::updateFromString("GTRKRecoTRK.threads=" + std::to_string(threads));
  auto input = study::readTrackingInput(inputFile);
  const auto truth = study::readTrackingTruth(inputFile + ".truth", input);
  o2::base::GeometryManager::loadGeometry(geometryFile, false, true);
  const auto bz = alice::alice3BzKiloGauss();
  o2::base::Propagator::Instance(true)->setNominalBz(static_cast<float>(bz));
  auto material = std::unique_ptr<o2::base::MatLayerCylSet>(o2::base::MatLayerCylSet::loadFromFile(materialFile));
  if (!material) {
    throw std::runtime_error("cannot load tracking material");
  }
  const auto frames = study::layerFrames<Det>(input);
  const auto& param = recoParam<Det>();
  workflow::TrackingInterface<Det> interface(!param.vertexer);
  interface.initGeometry(input.depths);
  const auto grid = workflow::makeGridGeometry<Det>(param, input.beamLine, frames, o2::raw::HBFUtils::Instance().getNOrbitsPerTF());
  interface.configure({bz, input.depths, material.get()}, {input.beamLine, grid});
  if (const auto* invalid = interface.checkInput(input.columns, input.offsets, grid.bucketCount(), input.vertices)) {
    throw std::runtime_error(invalid);
  }
  const TimeFrame<Det> frame{measurementView(input.columns), grid, input.offsets};
  const auto end = interface.endTimeframeOnExit();
  const auto vertexed = interface.vertices(frame, input.vertices, false);
  const auto tracked = interface.track(frame, vertexed.seeds);
  if (tracked.dropped || !tracked.result.tracks) {
    throw std::runtime_error("tracking failed");
  }
  const auto& records = *tracked.result.tracks;
  std::vector<Selected> selected;
  std::set<std::uint64_t> seen;
  const tof::Timestamp anchor{static_cast<std::uint64_t>(firstOrbit) * o2::constants::lhc::LHCMaxBunches, 0.F};
  for (std::size_t t = 0; t < records.size(); ++t) {
    const auto record = records.read(t);
    if (record.status != PropagationStatus::kOk) {
      continue;
    }
    const auto id = workflow::goodParticle(truth.rows, {record.clusters.index, record.count});
    const auto found = truth.particles.find(id);
    if (found == truth.particles.end() || found->second.origin != study::ParticleOrigin::kPhysicalPrimary ||
        found->second.pt < minPt || std::abs(found->second.eta) >= maxEta || ((id >> 50) & 255U) != 0 || !seen.insert(id).second) {
      continue;
    }
    selected.push_back({{.state = record.state, .outwardState = record.outwardState, .surface = interface.tracking().planeOf(record.state), .outwardSurface = interface.tracking().planeOf(record.outwardState), .boxStart = tof::shifted(anchor, record.timePS * 1e-3), .boxWidthNS = record.widthPS * 1e-3, .valid = true},
                        id,
                        found->second});
  }
  const auto population = selected.size();
  std::mt19937_64 random(seed);
  std::shuffle(selected.begin(), selected.end(), random);
  if (selected.size() > static_cast<std::size_t>(maxTracks)) {
    selected.resize(maxTracks);
  }
  std::cerr << "Tracked " << records.size() << "; population " << population << "; selected " << selected.size() << '\n';

  auto* geometry = o2::iotof::GeometryTGeo::Instance();
  geometry->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  std::vector<tof::Hit> hits;
  std::vector<std::vector<std::uint64_t>> hitTruth; // indexed by original cluster/measurement identity
  std::unordered_map<std::uint64_t, std::array<std::size_t, 2>> availability;
  std::size_t rawCount = 0, outsidePixels = 0, inputClusters = 0, ambiguousClusters = 0;
  const auto carries = [](const auto& labels, std::uint64_t identity) {
    return std::find(labels.begin(), labels.end(), identity) != labels.end();
  };
  if (clusterFile.empty()) {
    std::unique_ptr<o2::steer::DigitizationContext> context(o2::steer::DigitizationContext::loadFromFile(contextFile));
    if (!context) {
      throw std::runtime_error("cannot load collision context");
    }
    std::unordered_map<int, tof::Timestamp> eventTimes;
    const auto& parts = context->getEventParts();
    const auto& collisions = context->getEventRecords();
    for (std::size_t c = 0; c < parts.size(); ++c) {
      for (const auto& part : parts[c]) {
        if (part.sourceID != 0) {
          continue;
        }
        if (!eventTimes.emplace(part.entryID, tof::shifted({static_cast<std::uint64_t>(collisions[c].toLong()), 0.F}, collisions[c].getTimeOffsetWrtBC())).second) {
          throw std::runtime_error("reused source events require an occurrence-aware study; refusing to mix collisions");
        }
      }
    }
    std::vector<o2::iotof::Cluster> templates;
    const auto chips = geometry->getITOFNumberOfChips() + geometry->getOTOFNumberOfChips();
    for (int chip = 0; chip < chips; ++chip) {
      templates.emplace_back(0, 0, 1, 1, 1, 0, chip, 0.F, 0);
    }
    const auto sensors = tof::convertClusters(templates, *geometry, 0.020F);
    if (sensors.size() != templates.size()) {
      throw std::runtime_error("could not place every TOF sensor");
    }
    std::unique_ptr<TFile> file(TFile::Open(hitFile.c_str()));
    auto* tree = file ? file->Get<TTree>("o2sim") : nullptr;
    std::vector<o2::itsmft::Hit> raw;
    auto* rawPointer = &raw;
    if (!tree || tree->SetBranchAddress("TF3Hit", &rawPointer) < 0) {
      throw std::runtime_error("missing TF3 hits");
    }
    std::normal_distribution<double> smear(0., 0.020);
    const auto& specs = o2::iotof::ChipSpecificsParam::Instance();
    for (Long64_t event = 0; event < tree->GetEntries(); ++event) {
      tree->GetEntry(event);
      const auto collision = eventTimes.find(static_cast<int>(event));
      if (collision == eventTimes.end()) {
        continue;
      }
      for (const auto& point : raw) {
        ++rawCount;
        if (point.GetDetectorID() < 0 || point.GetDetectorID() >= chips || point.GetTrackID() < 0) {
          throw std::runtime_error("hit identity outside the supported source/barrel geometry");
        }
        auto measurement = sensors[point.GetDetectorID()];
        const Scalar global[3]{point.GetX(), point.GetY(), point.GetZ()};
        Scalar local[3];
        surfaceToLocal(measurement.surface, global, local);
        const auto row = std::lround((measurement.loc0 - local[0]) / specs.PitchRow);
        const auto column = std::lround((local[1] - measurement.loc1) / specs.PitchCol);
        if (row < 0 || row >= specs.NRows || column < 0 || column >= specs.NCols) {
          ++outsidePixels;
          continue;
        }
        measurement.loc0 -= static_cast<float>(row) * specs.PitchRow;
        measurement.loc1 += static_cast<float>(column) * specs.PitchCol;
        measurement.time = tof::shifted(collision->second, point.GetTime() * 1e9 + smear(random));
        measurement.index = static_cast<std::int32_t>(hits.size());
        const auto identity = o2::MCCompLabel(point.GetTrackID(), static_cast<int>(event), 0).getTrackEventSourceID();
        ++availability[identity][measurement.layer == tof::Layer::Inner ? 0 : 1];
        hits.push_back(measurement);
        hitTruth.push_back({identity});
      }
    }
  } else {
    std::unique_ptr<TFile> file(TFile::Open(clusterFile.c_str()));
    auto* tree = file ? file->Get<TTree>("o2sim") : nullptr;
    const auto* schema = file ? dynamic_cast<TStreamerInfo*>(file->GetStreamerInfoList()->FindObject("o2::iotof::Cluster")) : nullptr;
    if (!tree || tree->GetEntries() != 1 || !schema || schema->GetClassVersion() < 3) {
      throw std::runtime_error("provide one native cluster timeframe with BC/fine timing");
    }
    std::vector<o2::iotof::Cluster> clusters;
    auto* clusterPointer = &clusters;
    o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
    auto* labelPointer = &labels;
    if (tree->SetBranchAddress("TF3Cluster", &clusterPointer) < 0 ||
        tree->SetBranchAddress("TF3ClusterMCTruth", &labelPointer) < 0 || tree->GetEntry(0) <= 0) {
      throw std::runtime_error("missing native TOF clusters or MC labels");
    }
    inputClusters = clusters.size();
    if (labels.getIndexedSize() != inputClusters) {
      throw std::runtime_error("TOF cluster MC label rows do not match cluster rows");
    }
    hitTruth.resize(inputClusters);
    for (std::size_t cluster = 0; cluster < inputClusters; ++cluster) {
      auto& identities = hitTruth[cluster];
      for (const auto& label : labels.getLabels(cluster)) {
        if (label.isValid()) {
          identities.push_back(label.getTrackEventSourceID());
        }
      }
      std::sort(identities.begin(), identities.end());
      identities.erase(std::unique(identities.begin(), identities.end()), identities.end());
    }
    hits = tof::convertClusters(clusters, *geometry, 0.020F);
    for (const auto& hit : hits) {
      ambiguousClusters += hitTruth[hit.index].size() > 1;
      for (const auto identity : hitTruth[hit.index]) {
        ++availability[identity][hit.layer == tof::Layer::Inner ? 0 : 1];
      }
    }
  }
  std::vector<tof::Track> tracks;
  for (const auto& item : selected) {
    tracks.push_back(item.track);
  }
  tof::Parameters parameters;
  const auto acceptsTime = [&](const tof::Track& track, const tof::Hit& hit) {
    const double time = tof::differenceNS(hit.time, track.boxStart);
    return time >= -preselectionNS && time <= track.boxWidthNS + preselectionNS;
  };
  std::size_t preselectionLost = 0;
  for (const auto& item : selected) {
    for (std::size_t h = 0; h < hits.size(); ++h) {
      preselectionLost += carries(hitTruth[hits[h].index], item.particle) && !acceptsTime(item.track, hits[h]);
    }
  }
  std::size_t exactCalls = 0;
  std::vector<std::array<bool, 2>> trueSpatial(tracks.size()), trueTime(tracks.size());
  std::vector<std::array<double, 2>> minimumChi2(tracks.size(), {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()});
  const auto start = std::chrono::steady_clock::now();
  const auto matches = tof::match(tracks, hits, [&](const tof::Track& track, const tof::Hit& hit) {
    if (!acceptsTime(track, hit)) {
      return tof::Prediction{};
    }
    ++exactCalls;
    return tof::predict(track, hit, interface.propagator(), TGeoMaterial{}, parameters);
  });
  const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  // Truth is consulted only after matching. Evaluate correct-hit compatibility
  // independently of cluster ownership, so arbitration cannot appear as a
  // propagation failure merely because a hit was already claimed.
  std::vector<std::array<bool, 2>> geometric(tracks.size()), inside(tracks.size()), covarianceOK(tracks.size());
  for (std::size_t t = 0; t < tracks.size(); ++t) {
    const auto& track = tracks[t];
    for (std::size_t h = 0; h < hits.size(); ++h) {
      const auto& hit = hits[h];
      if (!carries(hitTruth[hits[h].index], selected[t].particle) || !acceptsTime(track, hit)) {
        continue;
      }
      const bool outer = hit.layer == tof::Layer::Outer;
      const auto layer = outer ? 1U : 0U;
      const auto spatial = tof::transport(outer ? track.outwardState : track.state, outer ? track.outwardSurface : track.surface,
                                          hit.surface, interface.propagator(), TGeoMaterial{}, PropagationDirection::kForward, parameters, false);
      geometric[t][layer] = geometric[t][layer] || spatial.ok;
      inside[t][layer] = inside[t][layer] || (spatial.ok && spatial.inside);
      const auto withCovariance = tof::transport(outer ? track.outwardState : track.state, outer ? track.outwardSurface : track.surface,
                                                 hit.surface, interface.propagator(), TGeoMaterial{}, PropagationDirection::kForward, parameters, true);
      covarianceOK[t][layer] = covarianceOK[t][layer] || (withCovariance.ok && withCovariance.inside);
      const auto prediction = tof::predict(track, hit, interface.propagator(), TGeoMaterial{}, parameters);
      if (std::isfinite(prediction.spatialChi2)) {
        minimumChi2[t][layer] = std::min(minimumChi2[t][layer], static_cast<double>(prediction.spatialChi2));
      }
      if (prediction.spatialChi2 >= 0.F && prediction.spatialChi2 <= parameters.maxSpatialChi2) {
        trueSpatial[t][layer] = true;
        for (std::size_t mass = 0; mass < tof::kHypotheses; ++mass) {
          const auto atReference = tof::differenceNS(hit.time, track.boxStart) - prediction.flightNS[mass];
          trueTime[t][layer] = trueTime[t][layer] || (prediction.valid[mass] && atReference >= -parameters.timeMarginNS && atReference < track.boxWidthNS + parameters.timeMarginNS);
        }
      }
    }
  }
  std::array<Counts, 2> counts;
  std::array<std::array<std::size_t, 6>, 2> losses{};
  std::size_t bothEligible = 0, bothCorrect = 0, guidedEligible = 0, guidedCorrect = 0;
  std::ofstream csv(output + ".csv");
  csv << "particle,pt_gen,eta_gen,pdg,inner_available,outer_available,inner_match,outer_match,inner_correct,outer_correct,inner_candidates,outer_candidates,truth_inner_min_chi2,truth_outer_min_chi2,truth_inner_time_candidate,truth_outer_time_candidate,pt_reco,inner_transport,inner_inside,inner_covariance,outer_transport,outer_inside,outer_covariance\n";
  for (std::size_t t = 0; t < tracks.size(); ++t) {
    const auto& item = selected[t];
    const auto& m = matches[t];
    const auto available = availability[item.particle];
    const std::array<int, 2> assigned{m.inner, m.outer};
    std::array<bool, 2> correct{};
    for (std::size_t layer = 0; layer < 2; ++layer) {
      const bool matched = assigned[layer] >= 0;
      correct[layer] = matched && carries(hitTruth[assigned[layer]], item.particle);
      counts[layer].ambiguous += matched && hitTruth[assigned[layer]].size() > 1;
      if (available[layer]) {
        ++counts[layer].eligible;
        counts[layer].spatial += trueSpatial[t][layer];
        counts[layer].coarseTime += trueTime[t][layer];
        if (!correct[layer]) {
          const auto reason = !geometric[t][layer] ? 0U : !inside[t][layer]       ? 1U
                                                        : !covarianceOK[t][layer] ? 2U
                                                        : !trueSpatial[t][layer]  ? 3U
                                                        : !trueTime[t][layer]     ? 4U
                                                                                  : 5U;
          ++losses[layer][reason];
        }
        counts[layer].matched += matched;
        counts[layer].correct += correct[layer];
        counts[layer].wrong += matched && !correct[layer];
      } else {
        counts[layer].withoutTruthHit += matched;
      }
    }
    bothEligible += available[0] && available[1];
    bothCorrect += available[0] && available[1] && correct[0] && correct[1];
    guidedEligible += correct[1] && available[0];
    guidedCorrect += correct[1] && available[0] && correct[0];
    csv << item.particle << ',' << item.truth.pt << ',' << item.truth.eta << ',' << item.truth.pdg << ',' << available[0] << ',' << available[1] << ','
        << m.inner << ',' << m.outer << ',' << correct[0] << ',' << correct[1] << ',' << m.innerCandidates << ',' << m.outerCandidates << ',' << minimumChi2[t][0] << ',' << minimumChi2[t][1] << ',' << trueTime[t][0] << ',' << trueTime[t][1] << ',' << std::sin(item.track.state.theta()) / std::abs(item.track.state.qOverP()) << ',' << geometric[t][0] << ',' << inside[t][0] << ',' << covarianceOK[t][0] << ',' << geometric[t][1] << ',' << inside[t][1] << ',' << covarianceOK[t][1] << '\n';
  }
  const auto lossReport = [&](std::size_t layer) {
    return JSON{{"propagation_or_stopped", losses[layer][0]}, {"mean_outside_sensor", losses[layer][1]}, {"covariance_transport", losses[layer][2]}, {"spatial_chi2", losses[layer][3]}, {"flight_or_coarse_time", losses[layer][4]}, {"guided_time_or_assignment", losses[layer][5]}};
  };
  JSON report{{"kind", clusterFile.empty() ? "hit-based matching diagnostic; not digitized-cluster efficiency" : "native TOF cluster matching; sampled reconstructed-track competition"},
              {"input", inputFile},
              {"geometry", geometryFile},
              {"tof_hits", hitFile},
              {"tof_clusters", clusterFile},
              {"input_clusters", inputClusters},
              {"unusable_clusters", inputClusters ? inputClusters - hits.size() : 0},
              {"ambiguous_usable_clusters", ambiguousClusters},
              {"context", contextFile},
              {"tracking_material", materialFile},
              {"matching_material", "TGeo mean material on each bounded step"},
              {"competition", "sampled reconstructed tracks only; complete measurement clutter"},
              {"timing_input", clusterFile.empty() ? "20 ps Gaussian added to transport hits" : "native cluster BC/fine time, no added smearing"},
              {"timing_sigma_ps", 20.0},
              {"tracking_seconds", tracked.seconds},
              {"nominal_bz_kG", bz},
              {"seed", seed},
              {"primary_min_pt_GeV", minPt},
              {"primary_max_abs_eta", maxEta},
              {"reconstructed_tracks", records.size()},
              {"selected_population", population},
              {"sampled_tracks", selected.size()},
              {"raw_hits", rawCount},
              {"outside_pixels", outsidePixels},
              {"usable_hits", hits.size()},
              {"usable_measurements", hits.size()},
              {"preselection_margin_ns", preselectionNS},
              {"true_hits_lost_by_preselection", preselectionLost},
              {"inner", counts[0].json()},
              {"outer", counts[1].json()},
              {"inner_losses", lossReport(0)},
              {"outer_losses", lossReport(1)},
              {"both_eligible", bothEligible},
              {"both_correct", bothCorrect},
              {"guided_inner_eligible", guidedEligible},
              {"guided_inner_correct", guidedCorrect},
              {"matching_seconds", seconds},
              {"exact_predictions", exactCalls}};
  std::ofstream(output + ".json") << report.dump(2) << '\n';
  std::cout << report.dump(2) << '\n';
  return preselectionLost == 0 ? 0 : 2;
}
} // namespace

int main(int argc, char** argv)
{
  try {
    return run(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
