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

/// \file TrackBench.h
/// \brief o2-gtrk-track-bench's run of one detector, `track<Det>`, and what it
///        uses: one translation unit per detector instantiates it
///        (trackBench<Det>.cxx); trackBench.cxx parses the options and dispatches.

#ifndef GTRK_STUDY_TRACKBENCH_H_
#define GTRK_STUDY_TRACKBENCH_H_

#include "GenericTracking/Study/Alice/Hits.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Study/Digest.h"
#include "GenericTracking/Study/ProfileReport.h"
#include "GenericTracking/Study/Statistics.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/Profile.h"
#include "GenericTracking/Study/StageComparison.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/Tracker.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/WorkCounters.h"
#include "GenericTracking/Study/Validation.h"
#include "GenericTracking/Workflow/GeometryProvider.h"
#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Workflow/MaterialLUT.h"
#include "GenericTracking/Workflow/RecoParam.h"
#include "GenericTracking/Workflow/TrackingInterface.h"
#include "GenericTracking/Study/TrackingInput.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DetectorsRaw/HBFUtils.h"
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/NameConf.h"
#include "DetectorsBase/MatLayerCylSet.h"

#include <TGeoManager.h>

#include "FitInputs.h"
#include <boost/program_options.hpp>
#include "GenericTracking/Study/StageFunnel.h"
#include "CompleteAnatomy.h"
#include "FollowAnatomy.h"
#include "SeedingStrategy.h"
#include "TrackBreakdown.h"
#include "VertexStudy.h"

#include <sys/resource.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <bit>
#include <type_traits>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace o2::gtrk::study
{

namespace bpo = boost::program_options;
template <DetectorLayout Det>
using Controller = Tracker<Det, o2::base::MatLayerCylSet>;

// Read once per run: its order follows the command line, not the padding.
struct Options { // NOLINT(clang-analyzer-optin.performance.Padding)
  std::string detector;
  std::string input;
  std::string truth;
  std::string config;
  std::string keyValues;
  std::string ccdbUrl; ///< --ccdb-url: the conditions' CCDB (a sample's snapshots: file://<run dir>/ccdb)
  int threads = -1;
  unsigned repeat = 3U;
  bool tree = true;
  bool iterationTrees = false;
  bool digestStages = false;
  bool funnel = false;
  int vertexStudy = -1;
  bool verticesOnly = false; ///< --vertices-only: the vertexing pass alone, its vertices against the truth
  double beamWidthCM = 0.;  ///< --beam-width: the beam's transverse width the windows take, instead of the dump's
  std::string vertexList;   ///< --vertex-list: the vertexer's vertices against the truth, as CSV; empty for none
  std::string vertexReplay; ///< --vertex-replay: the replay's vertices and collisions, as CSV with this prefix; empty for none
  std::string vertexPrior;  ///< --vertex-prior: an earlier pass's replay vertices the prior-* variants start from
  double findableMinPtGeV = 0.1;
  std::uint64_t findableAllLayers = 0U; ///< the score's "on every layer" mask; 0, every layer
  std::string seedCoverage;             ///< --seed-coverage: layer masks, comma separated
  std::string seedingStrategy;          ///< --seeding-strategy: the truth-removal study's JSON lines, one per pass
  std::string seedingSteps;             ///< --seeding-steps: its pT steps [GeV], comma separated
  bool seedingNoRemoval = false;        ///< --seeding-no-removal: every pass on the whole frame
  std::string seedingShadow;            ///< --seeding-shadow: the shadow passes, comma separated
  bool seedingFit = false;              ///< --seeding-fit: seeded by the road fit's pure records
  std::string seedingParticles;         ///< --seeding-particles: per findable primary, a CSV line
  std::string seedingPopulation;        ///< --seeding-population: the particles the study scores
  std::string seedingReal;              ///< --seeding-real: the real run's iterations scored as the study's passes, JSON lines
  std::string saveState;                ///< --save-state: the checkpoint after --state-iterations iterations
  std::size_t stateIterations = 0U;     ///< --state-iterations: how many iterations --save-state runs
  std::string resumeState;              ///< --resume-state: the run goes on from this checkpoint
  std::string holdClusters;             ///< --hold-clusters: clusters, one per line, taken before the first iteration
  std::string fitInputs;
  bool noProfile = false;
  bool dram = false;
  std::string trace;
  std::string json;
  std::string tracksOut;
  std::string clonesOut;   ///< --clones-out
  std::string rowsOut;     ///< --rows-out
  std::string fakesOut;    ///< --fakes-out
  std::string followSteps; ///< --follow-steps
  std::string cloneJoins;  ///< --clone-joins
  std::string backend;
  int gpuDevice = 0;
  std::size_t gpuArenaMB = 0U;
  bool hostReference = false;
  bool rebin = false;
  bool followStops = false;       ///< --follow-stops
  bool completeAnatomy = false;   ///< --complete-anatomy
  std::string stepAnatomy;        ///< --step-anatomy: the CSV of the replayed steps
  std::string truthWalk;          ///< --truth-walk: the CSV of the truth walks' steps
  std::string fateOut;            ///< --fate-out: the CSV of every class particle's fate
  std::string breakdown;          ///< --breakdown: where the last run's tracks sit, a JSON line
  std::string breakdownParticles; ///< --breakdown-particles: the same per particle of the class, a CSV
  std::string inwardHistory;      ///< --inward-history: the CSV of the inward refit's states
  std::string endDecisions;       ///< --end-decisions: the CSV of every follow walk's end
  std::string endEdits;           ///< --end-edits: with it, the CSV of each track's truth-free end edits
  int anatomyPdg = 0;             ///< --anatomy-pdg
  bool anatomyFindable = false;   ///< --anatomy-findable
  std::string anatomyBeamLine;    ///< --anatomy-beam-line
  std::string anatomyFollowLoss;  ///< --anatomy-follow-loss
  int truthWalkIteration = -1;    ///< --truth-walk-iteration: the iteration of the truth walk (-1: the first holding the class's mask)
  std::string hits;               ///< --hits: the MC hits (sgn_Hits<det>.root) the steps are held to
  std::string ft3Hits;            ///< --ft3-hits: FT3's MC hits beside TRK's
  std::string kine;               ///< --kine: the MC kinematics (sgn_Kine.root), the anatomy's clutter test
};

/// The lines of one section of an INI file, `[name]` excluded.
inline std::vector<std::string> iniSection(const std::string& path, const std::string& name)
{
  std::ifstream in{path};
  std::vector<std::string> lines;
  bool inside = false;
  for (std::string line; std::getline(in, line);) {
    if (!line.empty() && line.front() == '[') {
      inside = line == "[" + name + "]";
      continue;
    }
    if (inside && !line.empty()) {
      lines.push_back(line);
    }
  }
  return lines;
}

/// Load the sink's parameters and check that every line of them reads back as
/// written: a value that does not survive the round trip would make this a
/// different configuration from the sink's. A key the file lacks -- added to
/// RecoSettings after the dump -- is at its default, as the sink ran it;
/// a key RecoSettings does not have is an
/// error. \a section is the detector's parameter, GTRKReco<DET>.
inline void loadConfiguration(const std::string& path, const std::string& section)
{
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("no configuration " + path + "; the sink writes it beside the input");
  }
  if (iniSection(path, section).empty()) {
    throw std::runtime_error(std::format("{} has no [{}]: written before each detector had its own parameter; re-dump it", path, section));
  }
  o2::conf::ConfigurableParam::updateFromFile(path, section + ",HBFUtils");
  const auto echo = (std::filesystem::temp_directory_path() / std::format("gtrk-track-bench-{}.ini", getpid())).string();
  study::writeParameters(echo);
  for (const std::string& name : {section, std::string{"HBFUtils"}}) {
    const auto read = iniSection(echo, name);
    for (const auto& line : iniSection(path, name)) {
      if (std::ranges::find(read, line) == read.end()) {
        std::filesystem::remove(echo);
        throw std::runtime_error(std::format("[{}] of {} does not read back as written: {}", name, path, line));
      }
    }
  }
  std::filesystem::remove(echo);
}

/// One slice's stage outputs as digests: where two runs first part, when
/// their track digests differ.
template <DetectorLayout Det>
void printStageDigests(const typename Controller<Det>::SliceObservation& slice)
{
  const auto stage = [](const auto& output) {
    Digest digest;
    digest.column(output.offsets());
    digest.columns(output.view(), output.size());
    return digest.result();
  };
  const auto& result = *slice.result;
  Digest roads;
  roads.column(slice.roads->roads());
  roads.column(slice.roads->starts());
  roads.column(slice.roads->roadOffsets());
  std::cout << std::format(
    "  stage digests, iteration {} slice [{}, {}): tracklets {:016x} triplets {:016x} neighbours {:016x} roads {:016x} records {:016x}\n", slice.iteration,
    slice.slice->coreBlocks.begin, slice.slice->coreBlocks.end, stage(*result.tracklets), stage(*result.triplets), stage(*result.neighbours), roads.result(),
    stage(*slice.records));
}

/// Every accepted triplet's and neighbour's fit inputs of one slice, with the
/// scattering variances the stages fitted them with: the iteration's table,
/// a kinked triplet's middle layer scaled by the kink allowance
/// (`graphKinkScale2`) in its own fit and in its neighbours'.
template <DetectorLayout Det>
void collectFitInputs(const typename Controller<Det>::SliceObservation& slice, const IterationParameters<Det>& parameters, study::FitInputs& inputs)
{
  const auto& frame = *slice.frame;
  const auto& result = *slice.result;
  const auto tracklets = result.tracklets->view();
  const auto triplets = result.triplets->view();
  const auto neighbours = result.neighbours->view();
  const auto neighbourOffsets = result.neighbours->offsets();
  const auto point = [&](auto row) { return pairPoint(frame.measurement(row)); };
  const auto iteration = static_cast<std::uint32_t>(slice.iteration);
  const Scalar tripletChi2Max = parameters.pairCut.nSigma * parameters.pairCut.nSigma;
  const auto kinkScale = [&](std::size_t triplet) { return triplets.kinked[triplet] != 0U ? parameters.graphKinkScale2 : 1.F; };
  for (std::size_t triplet = 0; triplet < result.triplets->size(); ++triplet) {
    const auto first = triplets.first[triplet];
    const auto second = triplets.second[triplet];
    const LayerID middleLayer = parameters.linkTargets[tracklets.link[first]];
    const auto a = point(tracklets.source[first]);
    const auto b = point(tracklets.target[first]);
    const auto c = point(tracklets.target[second]);
    inputs.triplets.push_back({.first = a,
                               .middle = b,
                               .third = c,
                               .variance = tripletScattering(parameters.scattering, parameters.pairCut, middleLayer, a, b, c, detail::referenceSinTheta(a, c),
                                                             kinkScale(triplet)),
                               .chi2Max = tripletChi2Max,
                               .curvature = triplets.curvature[triplet],
                               .varCurvature = triplets.varCurvature[triplet],
                               .iteration = iteration});
  }
  for (std::size_t element = 0; element + 1U < neighbourOffsets.size(); ++element) {
    const auto shared = triplets.first[element];
    const auto last = triplets.second[element];
    const std::size_t link = tracklets.link[shared];
    const auto b = point(tracklets.source[shared]);
    const auto c = point(tracklets.target[shared]);
    const auto d = point(tracklets.target[last]);
    for (auto neighbour = neighbourOffsets[element]; neighbour < neighbourOffsets[element + 1U]; ++neighbour) {
      const auto from = neighbours.from[neighbour];
      const auto a = point(tracklets.source[triplets.first[from]]);
      const auto scattering =
        adjacentScattering(parameters.scattering, parameters.pairCut, parameters.linkSources[link], parameters.linkTargets[link], a, b, c, d,
                           detail::referenceSinTheta(a, c), detail::referenceSinTheta(b, d), kinkScale(from), kinkScale(element));
      inputs.neighbours.push_back({.a = a,
                                   .b = b,
                                   .c = c,
                                   .d = d,
                                   .varianceSecond = scattering.second,
                                   .varianceThird = scattering.third,
                                   .covariance = scattering.covariance,
                                   .chi2Max = parameters.neighbourChi2Max,
                                   .iteration = iteration});
    }
  }
}

/// One iteration's work counts, where the build counts (GTRK_WORK_COUNTERS).
inline void printWork(const WorkCounts& work, std::size_t iteration)
{
  using enum WorkCounter;
  const auto at = [&](WorkCounter counter) { return work[static_cast<std::size_t>(counter)]; };
  const auto perFit = [&](WorkCounter nanoseconds, WorkCounter fits) {
    return at(fits) == 0U ? 0. : static_cast<double>(at(nanoseconds)) / static_cast<double>(at(fits));
  };
  std::cout << std::format(
    "  work, iteration {}, triplets: {} first tracklets ({} without a middle), {} in the dip-angle window; rejected {} by time, {} by "
    "curvature; {} fits, {} failed, {} over chi2, {} over chi2 at their own momentum; {} accepted ({} with the kink allowance); {:.1f} ns/fit\n",
    iteration, at(kTripletFirsts), at(kTripletNoMiddle), at(kTripletWindow), at(kTripletTime), at(kTripletCurvature), at(kTripletFits), at(kTripletFitFailed),
    at(kTripletChi2), at(kTripletOwnMomentum), at(kTripletFits) - at(kTripletFitFailed) - at(kTripletChi2) - at(kTripletOwnMomentum), at(kTripletKinked),
    perFit(kTripletFitNS, kTripletFits));
  std::cout << std::format(
    "  work, iteration {}, neighbours: {} triplets, {} candidates ({} both kinked); rejected {} by curvature, {} by time; {} fits, {} failed, {} over "
    "chi2; {} accepted; {:.1f} ns/fit\n",
    iteration, at(kNeighbourTriplets), at(kNeighbourCandidates), at(kNeighbourBothKinked), at(kNeighbourCurvature), at(kNeighbourTime), at(kNeighbourFits),
    at(kNeighbourFitFailed), at(kNeighbourChi2), at(kNeighbourFits) - at(kNeighbourFitFailed) - at(kNeighbourChi2), perFit(kNeighbourFitNS, kNeighbourFits));
  std::cout << std::format(
    "  work, iteration {}, roads: {} without a common time, {} fitted, {} of {} filter steps taken; rejected {} by chi2, {} by propagation, {} "
    "below the pT floor, {} by the beam line; {} accepted ({} of them without their innermost clusters, of {} such fits)\n",
    iteration, at(kRoadTime), at(kRoadFits), at(kRoadSteps), at(kRoadStepsPossible), at(kRoadChi2), at(kRoadPropagation), at(kRoadFloor), at(kRoadBeam),
    at(kRoadAccepted), at(kRoadTrimAccepted), at(kRoadTrimFits));
  std::cout << std::format(
    "  work, iteration {}, tracklets: {} sources, {} queries, {} not made as every pair would cross another layer; {} candidates held by earlier "
    "tracks, {} crossing another layer, {} refused at the bound on a kink between the layers, {} pair tests, {} past the transverse test, {} past "
    "the longitudinal one, {} accepted; inside their query's window: {} candidates, {} refused at the bound, {} pair tests, {} longitudinal tests; {} accepted "
    "outside it\n",
    iteration, at(kTrackletSources), at(kTrackletQueries), at(kTrackletHoleQueries), at(kTrackletUsed), at(kTrackletHoles), at(kPairBetweenBound),
    at(kPairTests), at(kPairTransverse), at(kPairLongitudinal), at(kTrackletAccepted), at(kTrackletInWindow), at(kPairBetweenBoundIn), at(kPairTestsIn),
    at(kPairLongitudinalTestsIn), at(kTrackletAcceptedOutside));
  std::cout << std::format(
    "  work, iteration {}, grid (tracklets): {} queries over {} blocks; {} ranges over {} cells, {} of them searched; {} "
    "rows, {} in the time bound, {} past the time cut\n",
    iteration, at(kGridQueries), at(kGridBlocks), at(kGridRanges), at(kGridCells), at(kGridSearchedCells), at(kGridRows), at(kGridInBound), at(kGridVisits));
  std::cout << std::format(
    "  work, iteration {}, grid (follow): {} queries over {} blocks; {} ranges over {} cells, {} of them searched; {} "
    "rows, {} in the time bound, {} past the time cut; {} updates; {} rows proposed, {} outbid, {} with a failed refit\n",
    iteration, at(kFollowGridQueries), at(kFollowGridBlocks), at(kFollowGridRanges), at(kFollowGridCells), at(kFollowGridSearchedCells), at(kFollowGridRows),
    at(kFollowGridInBound), at(kFollowGridVisits), at(kFollowUpdates), at(kFollowProposed), at(kFollowOutbid), at(kFollowRefitRejected));
  std::cout << std::format(
    "  work, iteration {}, follow: {} layers without a row tested, {} of them holes; {} lookaheads; {} kink passes, {} taken, {} "
    "updates, {} witness updates, {} rows dropped with a witness; {} last rows unconfirmed, {} kept on their density, {} refused at "
    "the layout's end, {} against the beam line; {} road ends fitted again, {} replaced\n",
    iteration, at(kFollowHoleTests), at(kFollowHoles), at(kFollowLookaheads), at(kFollowKinkSteps), at(kFollowKinksTaken), at(kFollowKinkUpdates),
    at(kFollowWitnessUpdates), at(kFollowWitnessTrimmed), at(kFollowUnconfirmed), at(kFollowDensityKept), at(kFollowEndRefused), at(kFollowBeamEndRefused),
    at(kFollowEndFits), at(kFollowEndReplaced));
  if (at(kTrackletShape) != 0U || at(kTripletShape) != 0U) {
    std::cout << std::format("  work, iteration {}, shapes: {} tracklet candidates and {} thirds refused\n", iteration, at(kTrackletShape),
                             at(kTripletShape));
  }
  if (at(kFollowTrimFits) != 0U) {
    std::cout << std::format("  work, iteration {}, follow trim: {} tracks refitted without their last outward row, {} trimmed, {} updates\n",
                             iteration, at(kFollowTrimFits), at(kFollowTrimmed), at(kFollowTrimUpdates));
  }
  if (at(kFollowEarlierProposed) != 0U) {
    std::cout << std::format("  work, iteration {}, earlier iterations' follow: {} rows proposed, {} outbid, {} with a failed refit\n", iteration,
                             at(kFollowEarlierProposed), at(kFollowEarlierOutbid), at(kFollowEarlierRefitRejected));
  }
  std::cout << std::format("  work, iteration {}, propagations: {} plane crossings, {} at |n . d| <= 0.1\n", iteration, at(kPropagations),
                           at(kPropagationsGrazing));
  std::cout << std::format(
    "  work, iteration {}, tracklet sort: {} sources with 2+ tracklets ({} tracklets, {:.1f} each); {} with more than 8, {} more than 32 ({} "
    "tracklets), {} more than 128, {} more than 512\n",
    iteration, at(kSortSegments), at(kSortTracklets),
    at(kSortSegments) == 0U ? 0. : static_cast<double>(at(kSortTracklets)) / static_cast<double>(at(kSortSegments)), at(kSortSegmentsOver8),
    at(kSortSegmentsOver32), at(kSortTrackletsOver32), at(kSortSegmentsOver128), at(kSortSegmentsOver512));
}

/// One iteration's tracklet counts per link, the links with a query.
template <DetectorLayout Det>
void printLinkWork(const LinkWorkCounts& work, std::size_t iteration)
{
  for (std::size_t link = 0; link < work.size(); ++link) {
    const auto& counts = work[link];
    if (counts[static_cast<std::size_t>(LinkCounter::kQueries)] == 0U) {
      continue;
    }
    std::cout << std::format("  work, iteration {}, link {} ({} -> {}): {} queries, {} pair candidates, {} tracklets\n", iteration, link,
                             Det::layerLinks[link].from, Det::layerLinks[link].to, counts[static_cast<std::size_t>(LinkCounter::kQueries)],
                             counts[static_cast<std::size_t>(LinkCounter::kCandidates)], counts[static_cast<std::size_t>(LinkCounter::kTracklets)]);
  }
}

/// \a values as a JSON array.
template <std::size_t N>
std::string jsonArray(const std::array<std::uint64_t, N>& values)
{
  std::string out = "[";
  for (std::size_t index = 0; index < N; ++index) {
    out += std::format("{}{}", index == 0U ? "" : ",", values[index]);
  }
  return out + "]";
}

template <std::size_t N>
std::string jsonArray(const std::array<double, N>& values)
{
  std::string out = "[";
  for (std::size_t index = 0; index < N; ++index) {
    out += std::format("{}{:.6g}", index == 0U ? "" : ",", values[index]);
  }
  return out + "]";
}

inline std::string jsonScore(const study::FindableScore& r)
{
  std::string species;
  for (const auto& [pdg, counts] : r.bySpecies) {
    species += std::format(R"({}"{}":{})", species.empty() ? "" : ",", pdg, jsonArray(counts));
  }
  return std::format(R"({{"findable":{},"found":{},"fakeOnly":{},"cloned":{},"cloneTracks":{},"withFake":{},"findableByPt":{},"foundByPt":{},)"
                     R"("findableByEta":{},"foundByEta":{},"findableByPhi":{},"foundByPhi":{},"complete":{},"completeByPt":{},"layersFound":{},)"
                     R"("layersOf":{},"layerFraction":{:.6g},"layerFractionByPt":{},"bySpecies":{{{}}}}})",
                     r.findable, r.found, r.fakeOnly, r.cloned, r.cloneTracks, r.withFake, jsonArray(r.findableByPt), jsonArray(r.foundByPt),
                     jsonArray(r.findableByEta), jsonArray(r.foundByEta), jsonArray(r.findableByPhi), jsonArray(r.foundByPhi), r.complete,
                     jsonArray(r.completeByPt), r.layersFound, r.layersOf, r.layerFraction, jsonArray(r.layerFractionByPt), species);
}

/// One score per particle origin, in ParticleOrigin's order.
template <std::size_t N>
std::string jsonScores(const std::array<study::FindableScore, N>& scores)
{
  std::string out = "[";
  for (std::size_t origin = 0; origin < scores.size(); ++origin) {
    out += (origin == 0U ? "" : ",") + jsonScore(scores[origin]);
  }
  return out + "]";
}

/// The decay daughters: per DecayFamily, per production bin, [on minLayers
/// or more, found of them, findable, found of them].
inline std::string jsonDecays(const study::TrackScore& score, const study::ComparisonSettings& settings)
{
  std::string edges = "[";
  for (std::size_t index = 0; index < settings.productionEdges.size(); ++index) {
    edges += std::format("{}{}", index == 0U ? "" : ",", settings.productionEdges[index]);
  }
  std::string families = "[";
  for (std::size_t family = 0; family < score.decays.size(); ++family) {
    families += family == 0U ? "[" : ",[";
    for (std::size_t bin = 0; bin < score.decays[family].size(); ++bin) {
      const auto& counts = score.decays[family][bin];
      families += std::format("{}[{},{},{},{}]", bin == 0U ? "" : ",", counts.onMinLayers, counts.foundOnMinLayers, counts.findable, counts.found);
    }
    families += "]";
  }
  return std::format(R"({{"productionEdges":{}],"alongZ":{},"byFamily":{}]}})", edges, settings.productionAlongZ, families);
}

/// \a score as JSON: counts only; the bins are TrackComparison.h's, the
/// species |PDG codes|.
inline std::string jsonScore(const study::TrackScore& score, const study::ComparisonSettings& settings, double seconds)
{
  return std::format(R"({{"seconds":{:.6f},"tracks":{},"fakes":{},"tracksByLength":{},"fakesByLength":{},"allLayers":{},"minLayers":{},"o2":{},)"
                     R"("originAllLayers":{},"originMinLayers":{},"originFindable":{},"primaryRegion":{},"decays":{}}})",
                     seconds, score.tracks, score.fakes, jsonArray(score.tracksByLength), jsonArray(score.fakesByLength), jsonScore(score.allLayers),
                     jsonScore(score.minLayers), jsonScore(score.o2), jsonScores(score.originAllLayers), jsonScores(score.originMinLayers),
                     jsonScores(score.originFindable), jsonScores(score.primaryRegion), jsonDecays(score, settings));
}

/// The seed sets' coverage of the findable: by class, the fraction with a measurement on every layer of a set,
/// and on all but one, for each set of \a masks and for any of them.
inline void printSeedCoverage(const study::TrackingTruth& truth, const study::ComparisonSettings& comparison, std::span<const std::uint64_t> masks)
{
  using study::ParticleOrigin;
  constexpr std::size_t kClasses = 5U;
  constexpr std::array<std::string_view, kClasses> kNames{"primaries", "barrel", "disc", "pT 0.1-0.2", "decays"};
  // Layer counts up to this many, the last slot holding the longer ones.
  constexpr std::size_t kLengthSlots = 12U;
  struct Counts {
    std::uint64_t findable = 0U;
    std::vector<std::uint64_t> all;
    std::vector<std::uint64_t> allButOne;
    std::array<std::uint64_t, kLengthSlots + 1U> layers{};
  };
  std::array<Counts, kClasses> counts;
  for (auto& count : counts) {
    count.all.assign(masks.size() + 1U, 0U);
    count.allButOne.assign(masks.size() + 1U, 0U);
  }
  for (const auto& [id, layers] : study::particleLayers(truth.rows, comparison)) {
    const auto found = truth.particles.find(id);
    if (found == truth.particles.end() || !study::isFindable(found->second.pt, layers, comparison)) {
      continue;
    }
    const auto& particle = found->second;
    std::array<bool, kClasses> in{};
    if (particle.origin == ParticleOrigin::kPhysicalPrimary) {
      const bool disc = (layers.layers & ~comparison.barrelLayers) != 0U;
      in = {true, !disc, disc, particle.pt < 0.2F, false};
    } else if (particle.origin == ParticleOrigin::kDecay) {
      in[4] = true;
    }
    bool anyAll = false;
    bool anyAllButOne = false;
    std::vector<bool> all(masks.size());
    std::vector<bool> allButOne(masks.size());
    for (std::size_t set = 0; set < masks.size(); ++set) {
      const int missing = std::popcount(masks[set] & ~layers.layers);
      all[set] = missing == 0;
      allButOne[set] = missing <= 1;
      anyAll = anyAll || all[set];
      anyAllButOne = anyAllButOne || allButOne[set];
    }
    for (std::size_t c = 0; c < kClasses; ++c) {
      if (!in[c]) {
        continue;
      }
      ++counts[c].findable;
      ++counts[c].layers[std::min<std::size_t>(static_cast<std::size_t>(std::popcount(layers.layers)), kLengthSlots)];
      for (std::size_t set = 0; set < masks.size(); ++set) {
        counts[c].all[set] += all[set] ? 1U : 0U;
        counts[c].allButOne[set] += allButOne[set] ? 1U : 0U;
      }
      counts[c].all.back() += anyAll ? 1U : 0U;
      counts[c].allButOne.back() += anyAllButOne ? 1U : 0U;
    }
  }
  using study::percent;
  std::cout << "seed coverage of the findable, % on every layer of a set (on all but one):\n";
  for (std::size_t c = 0; c < kClasses; ++c) {
    std::string line = std::format("  {:<11} {:8d}:", kNames[c], counts[c].findable);
    for (std::size_t set = 0; set <= masks.size(); ++set) {
      line += std::format(" {} {:.2f} ({:.2f})", set < masks.size() ? std::format("set {}", masks[set]) : std::string{"any"},
                          percent(counts[c].all[set], counts[c].findable), percent(counts[c].allButOne[set], counts[c].findable));
    }
    std::cout << line << '\n';
  }
  std::cout << "  the same by the layers they have a measurement on, % with at least n:\n";
  for (std::size_t c = 0; c < kClasses; ++c) {
    std::string line = std::format("  {:<11}", kNames[c]);
    std::uint64_t atLeast = counts[c].findable;
    for (std::size_t n = 0; n <= kLengthSlots; ++n) {
      if (n >= 4U) {
        line += std::format(" {}{}: {:.2f}", n, n == kLengthSlots ? "+" : "", percent(atLeast, counts[c].findable));
      }
      atLeast -= counts[c].layers[n];
    }
    std::cout << line << '\n';
  }
}

/// The process's peak resident memory so far [bytes].
inline std::uint64_t maxRssBytes()
{
  rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
  constexpr std::uint64_t kBytesPerKiB = 1024U;
  return static_cast<std::uint64_t>(usage.ru_maxrss) * kBytesPerKiB; // NOLINT(cppcoreguidelines-pro-type-union-access): glibc's rusage
}

/// \a tracks into \a path as text, one per line: the iteration, the clusters
/// as cluster numbers, the chi2 -- to diff two runs' tracks.
template <DetectorLayout Det>
void writeTracks(const std::string& path, const typename TrackRefit<Det>::Tracks& tracks, const TimeFrame<Det>& frame)
{
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot write " + path);
  }
  const auto columns = frame.measurements();
  const auto view = tracks.view();
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    out << static_cast<int>(view.iteration[track]);
    for (std::size_t at = 0; at < view.count[track]; ++at) {
      out << ' ' << columns.clusterIndex[view.clusters[track].index[at]];
    }
    out << std::format(" chi2 {:.6g}\n", static_cast<double>(view.chi2[track]));
  }
}

/// The checkpoint file of --save-state and --resume-state: a header, then
/// each part raw -- the first iterations' settings, which --resume-state checks against its own, their summaries, the used mask and
/// the tracks. The header's fingerprints hold what else the iterations read:
/// the frame (every measurement column, so also the grid's row order, which
/// the configuration's every iteration sizes) and the conditions (the beam
/// line, the dump's vertices, the slicing). Read back by the same build only
/// (the header's sizes; the caller's file names carry the build).
struct CheckpointHeader {
  std::array<char, 8> magic{'G', 'T', 'R', 'K', 'C', 'K', 'P', '2'};
  std::uint64_t frame = 0U;
  std::uint64_t conditions = 0U;
  std::uint64_t iterations = 0U;
  std::uint64_t rows = 0U;
  std::uint64_t tracks = 0U;
  std::uint64_t settingsBytes = 0U;
  std::uint64_t summaryBytes = 0U;
  std::uint64_t trackBytes = 0U;
};

/// The checkpoint's frame fingerprint: every measurement column of \a frame.
template <DetectorLayout Det>
std::uint64_t checkpointFrame(const TimeFrame<Det>& frame)
{
  Digest digest;
  digest.columns(frame.measurements(), frame.size());
  return digest.result();
}

/// Its conditions' fingerprint: \a beamLine, \a vertices, \a blocksPerSlice and \a haloPS.
inline std::uint64_t checkpointConditions(const workflow::BeamLine& beamLine, std::span<const VertexSeed> vertices, std::int32_t blocksPerSlice,
                                          std::int64_t haloPS)
{
  Digest digest;
  for (const double value : {beamLine.x, beamLine.y, beamLine.z, beamLine.sigmaX, beamLine.sigmaY, beamLine.sigmaZ}) {
    digest.add(value);
  }
  digest.add(vertices.size());
  for (const auto& vertex : vertices) {
    digest.add(vertex.timePS);
    digest.add(vertex.widthPS);
    digest.add(vertex.z);
    digest.add(vertex.varZ);
    digest.add(vertex.varTransverse);
  }
  digest.add(blocksPerSlice);
  digest.add(haloPS);
  return digest.result();
}

template <typename Checkpoint>
void writeCheckpoint(const std::string& path, const Checkpoint& state, std::span<const IterationSettings> settings, std::uint64_t frame,
                     std::uint64_t conditions)
{
  using Summary = typename decltype(state.summaries)::value_type;
  using Record = typename decltype(state.tracks)::value_type;
  static_assert(std::is_trivially_copyable_v<IterationSettings> && std::is_trivially_copyable_v<Summary> && std::is_trivially_copyable_v<Record>);
  std::ofstream out(path, std::ios::binary);
  if (!out) {
    throw std::runtime_error("cannot write " + path);
  }
  const CheckpointHeader header{.frame = frame,
                                .conditions = conditions,
                                .iterations = state.summaries.size(),
                                .rows = state.used.size(),
                                .tracks = state.tracks.size(),
                                .settingsBytes = sizeof(IterationSettings),
                                .summaryBytes = sizeof(Summary),
                                .trackBytes = sizeof(Record)};
  const auto write = [&](const auto* data, std::size_t count) {
    out.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(count * sizeof(*data))); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
  };
  write(&header, 1U);
  write(settings.data(), state.summaries.size());
  write(state.summaries.data(), state.summaries.size());
  write(state.used.data(), state.used.size());
  write(state.tracks.data(), state.tracks.size());
  if (!out) {
    throw std::runtime_error("cannot write " + path);
  }
}

/// The checkpoint in \a path, whose iterations must be the first of \a settings,
/// on the frame and conditions of the fingerprints \a frame and \a conditions.
template <typename Checkpoint>
Checkpoint readCheckpoint(const std::string& path, std::span<const IterationSettings> settings, std::uint64_t frame, std::uint64_t conditions)
{
  using Summary = typename decltype(Checkpoint::summaries)::value_type;
  using Record = typename decltype(Checkpoint::tracks)::value_type;
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("cannot read " + path);
  }
  const auto read = [&](auto* data, std::size_t count) {
    in.read(reinterpret_cast<char*>(data), static_cast<std::streamsize>(count * sizeof(*data))); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    if (!in) {
      throw std::runtime_error(path + " is cut short");
    }
  };
  CheckpointHeader header;
  const auto magic = header.magic;
  read(&header, 1U);
  if (header.magic != magic || header.settingsBytes != sizeof(IterationSettings) || header.summaryBytes != sizeof(Summary) ||
      header.trackBytes != sizeof(Record)) {
    throw std::runtime_error(path + " is not a checkpoint of this build");
  }
  const auto bytes = sizeof(header) + (header.iterations * (header.settingsBytes + header.summaryBytes)) + header.rows + (header.tracks * header.trackBytes);
  if (std::filesystem::file_size(path) != bytes) {
    throw std::runtime_error(std::format("{} has {} bytes, its header promises {}", path, std::filesystem::file_size(path), bytes));
  }
  if (header.frame != frame || header.conditions != conditions) {
    throw std::runtime_error(path + " was saved on another frame or with other conditions (beam line, vertices, slicing)");
  }
  if (header.iterations > settings.size()) {
    throw std::runtime_error(std::format("{} holds {} iterations, the configuration has {}", path, header.iterations, settings.size()));
  }
  std::vector<IterationSettings> saved(header.iterations);
  read(saved.data(), saved.size());
  for (std::size_t iteration = 0; iteration < saved.size(); ++iteration) {
    if (!(saved[iteration] == settings[iteration])) {
      throw std::runtime_error(std::format("{}: iteration {} was run with other settings than the configuration's", path, iteration));
    }
  }
  Checkpoint state;
  state.summaries.resize(header.iterations);
  state.used.resize(header.rows);
  state.tracks.resize(header.tracks);
  read(state.summaries.data(), state.summaries.size());
  read(state.used.data(), state.used.size());
  read(state.tracks.data(), state.tracks.size());
  return state;
}

/// The follow's rows of \a iteration's own tracks, inward and outward, into
/// \a followed by track (for --fakes-out).
template <typename Observation>
void recordFollowed(const Observation& iteration, std::vector<std::array<std::uint8_t, 2>>& followed)
{
  if (iteration.iteration == 0U) {
    followed.clear();
  }
  followed.resize(iteration.tracks->size());
  for (std::size_t i = 0; i < iteration.proposals.size(); ++i) {
    followed[iteration.trackOffset + i] = {iteration.proposals[i].inwardCount, iteration.proposals[i].outwardCount};
  }
}

/// Every track of the last run against \a truth, for the study of the
/// follow's fakes: one `T` line per
/// track -- iteration, clusters, its own iteration's follow's rows inward and
/// outward (\a followed, from the iteration observer) and the earlier
/// iterations' follows' (`followEarlier`: what the track added beyond its
/// own follow's), good (1) or not, its majority particle's rows on it, then
/// of the rows not carrying that particle those in the road, its follow's
/// inward, outward and the earlier follows', and how many of them carry no
/// particle at all;
/// then the majority particle (0 for noise only), its origin (0 primary, 1
/// decay, 2 material, 3 none), pT, whether it is a findable physical
/// primary, its rows in the frame, of them held by a good track of its own,
/// held by another track, and the layer of the track's first wrong row;
/// then the track's chi2 and, from \a stagesOf (the funnel, with --funnel),
/// that particle's longest chain of its own in the track's iteration at each
/// stage, tracklets to the road fit, its road loss, its fit loss and its
/// adjacent pairs not kept by each PairLoss, then the best chi2 of its
/// longest records and the chi2 of the track's own record (\a recordChi2 of
/// its road's rows) (-1 each without), then the track's eta at its fit's
/// state, its majority particle's (-99 for noise only) and the bits of the
/// layers its rows sit on.
/// A fake has a `W` line per wrong row: the track, the row's position on it
/// and the road's first and end positions, its layer, its part (0 road, 1
/// follow inward, 2 outward, 3 earlier follows), whether it carries no
/// particle, its first particle, that one's origin and pT and whether it has
/// a good track, then the majority particle's rows on that layer: how many,
/// free, held by a good track of its own, held by another track; then the
/// wrong row's errors in r-phi and z [cm], and for the nearest of those own
/// rows its errors and its distance to the wrong row, transverse and in z
/// [cm] (-1 each without one); then the PDG codes of the majority particle,
/// of the wrong row's particle and of that one's mother, and the wrong row's
/// particle's production radius [cm] (for the road fakes).
template <DetectorLayout Det, typename StagesOf, typename RecordChi2>
void writeFakes(const std::string& path, const typename TrackRefit<Det>::Tracks& tracks, const study::TrackingTruth& truth,
                const study::ComparisonSettings& settings, const TimeFrame<Det>& frame, const std::vector<std::array<std::uint8_t, 2>>& followed,
                StagesOf stagesOf, RecordChi2 recordChi2)
{
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot write " + path);
  }
  const auto view = tracks.view();
  const auto rows = trackRows(tracks);
  std::vector<std::uint64_t> goodOf(rows.size(), study::kNoTrackParticle);
  // Per row, the track holding it.
  std::vector<std::uint32_t> holder;
  holder.assign(truth.rows.size(), std::numeric_limits<std::uint32_t>::max());
  for (std::size_t track = 0; track < rows.size(); ++track) {
    goodOf[track] = study::goodParticle(truth.rows, rows[track]);
    for (const auto row : rows[track]) {
      holder[row] = static_cast<std::uint32_t>(track);
    }
  }
  auto rowsOf = study::rowsOfParticles(truth.rows);
  const auto layers = study::particleLayers(truth.rows, settings);
  const std::unordered_set<std::uint64_t> withGood(goodOf.begin(), goodOf.end());
  for (std::size_t track = 0; track < rows.size(); ++track) {
    const auto& own = rows[track];
    const auto count = static_cast<std::uint32_t>(own.size());
    // The rows in order: the earlier follows' inward, the own follow's
    // inward, the road, the own follow's outward, the earlier follows'
    // outward.
    const std::uint32_t followIn = track < followed.size() ? followed[track][0] : 0U;
    const std::uint32_t followOut = track < followed.size() ? followed[track][1] : 0U;
    const std::uint32_t earlierIn = view.addedInward[track] >= followIn ? view.addedInward[track] - followIn : 0U;
    const std::uint32_t earlierOut = view.addedOutward[track] >= followOut ? view.addedOutward[track] - followOut : 0U;
    const auto majority = study::majorityParticle(truth.rows, own);
    std::array<std::uint32_t, 4> wrong{};
    std::uint32_t noise = 0U;
    int firstWrongLayer = -1;
    std::vector<std::pair<std::uint32_t, std::size_t>> wrongRows; // position, part
    for (std::uint32_t i = 0; i < count; ++i) {
      const auto& at = truth.rows[own[i]];
      if (majority.particle != study::kNoTrackParticle && at.carries(majority.particle)) {
        continue;
      }
      std::size_t part = 0U; // the road
      if (i < earlierIn || i >= count - earlierOut) {
        part = 3U;
      } else if (i < earlierIn + followIn) {
        part = 1U;
      } else if (i >= count - earlierOut - followOut) {
        part = 2U;
      }
      ++wrong[part];
      wrongRows.emplace_back(i, part);
      noise += at.count == 0U ? 1U : 0U;
      if (firstWrongLayer < 0) {
        firstWrongLayer = at.layer;
      }
    }
    int origin = 3;
    float pt = 0.F;
    int findable = 0;
    std::uint32_t total = 0U;
    std::uint32_t heldGood = 0U;
    std::uint32_t heldOther = 0U;
    if (majority.particle != study::kNoTrackParticle) {
      const auto found = truth.particles.find(majority.particle);
      if (found != truth.particles.end()) {
        origin = static_cast<int>(found->second.origin);
        pt = found->second.pt;
        const auto chain = layers.find(majority.particle);
        findable = found->second.origin == study::ParticleOrigin::kPhysicalPrimary && chain != layers.end() &&
                       study::isFindable(found->second.pt, chain->second, settings)
                     ? 1
                     : 0;
      }
      for (const auto row : rowsOf[majority.particle]) {
        ++total;
        const auto by = holder[row];
        if (by == std::numeric_limits<std::uint32_t>::max() || by == track) {
          continue;
        }
        (goodOf[by] == majority.particle ? heldGood : heldOther) += 1U;
      }
    }
    out << std::format("T {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {:.4f} {} {} {} {} {}", static_cast<int>(view.iteration[track]), count, followIn,
                       followOut, earlierIn + earlierOut, goodOf[track] != study::kNoTrackParticle ? 1 : 0, majority.carried, wrong[0], wrong[1], wrong[2],
                       wrong[3], noise, majority.particle == study::kNoTrackParticle ? 0U : majority.particle, origin, findable, static_cast<double>(pt),
                       total, heldGood, heldOther, firstWrongLayer, track);
    out << std::format(" {:.6g}", static_cast<double>(view.chi2[track]));
    const auto stages = majority.particle == study::kNoTrackParticle ? std::nullopt : stagesOf(majority.particle, view.iteration[track]);
    for (std::size_t stage = 0; stage < 5U; ++stage) {
      out << ' ' << (stages ? static_cast<int>(stages->longest[stage]) : -1);
    }
    out << ' ' << (stages ? static_cast<int>(stages->roadLoss) : -1) << ' ' << (stages ? static_cast<int>(stages->fitLoss) : -1);
    for (std::size_t reason = 0; reason < study::kPairLosses; ++reason) {
      out << ' ' << (stages ? static_cast<std::int64_t>(stages->lostPairs[reason]) : -1);
    }
    const std::uint32_t roadBegin = earlierIn + followIn;
    const std::uint32_t roadEnd = count - earlierOut - followOut;
    const auto ownRecord = recordChi2(view.iteration[track], std::span<const std::uint32_t>{own}.subspan(roadBegin, roadEnd - roadBegin));
    out << std::format(" {:.6g} {:.6g}", stages ? static_cast<double>(stages->recordChi2) : -1., ownRecord ? static_cast<double>(*ownRecord) : -1.);
    const auto majorityFound = truth.particles.find(majority.particle);
    std::uint64_t layerBits = 0U;
    for (const auto row : own) {
      layerBits |= std::uint64_t{1} << truth.rows[row].layer;
    }
    out << std::format(" {:.4f} {:.4f} {}\n", -std::log(std::tan(static_cast<double>(view.state[track].theta()) / 2.)),
                       majorityFound != truth.particles.end() ? static_cast<double>(majorityFound->second.eta) : -99., layerBits);
    if (goodOf[track] != study::kNoTrackParticle) {
      continue;
    }
    for (const auto& [position, part] : wrongRows) {
      const auto& at = truth.rows[own[position]];
      const std::uint64_t particle = at.count == 0U ? 0U : at.particles[0];
      int wrongOrigin = 3;
      float wrongPt = 0.F;
      std::array<std::int32_t, 3> pdgs{}; // the majority's, the wrong particle's, its mother's
      float wrongRadius = -1.F;
      if (const auto found = truth.particles.find(particle); at.count != 0U && found != truth.particles.end()) {
        wrongOrigin = static_cast<int>(found->second.origin);
        wrongPt = found->second.pt;
        pdgs[1] = found->second.pdg;
        pdgs[2] = found->second.motherPdg;
        wrongRadius = std::hypot(found->second.vx, found->second.vy);
      }
      if (const auto found = truth.particles.find(majority.particle); found != truth.particles.end()) {
        pdgs[0] = found->second.pdg;
      }
      std::array<std::uint32_t, 4> onLayer{}; // all, free, held by a good track of its own, held by another
      const auto wrongAt = frame.measurement(own[position]);
      std::array<double, 4> nearest{-1., -1., -1., -1.}; // its sigma r-phi and z, distance transverse and in z
      double nearestDistance = std::numeric_limits<double>::max();
      if (majority.particle != study::kNoTrackParticle) {
        for (const auto row : rowsOf[majority.particle]) {
          if (truth.rows[row].layer != at.layer) {
            continue;
          }
          const auto ownAt = frame.measurement(row);
          const double transverse =
            std::hypot(static_cast<double>(ownAt.x) - static_cast<double>(wrongAt.x), static_cast<double>(ownAt.y) - static_cast<double>(wrongAt.y));
          const double dz = std::abs(static_cast<double>(ownAt.z) - static_cast<double>(wrongAt.z));
          if (std::hypot(transverse, dz) < nearestDistance) {
            nearestDistance = std::hypot(transverse, dz);
            nearest = {std::sqrt(static_cast<double>(ownAt.varRPhi)), std::sqrt(static_cast<double>(ownAt.varZ)), transverse, dz};
          }
          ++onLayer[0];
          const auto by = holder[row];
          if (by == std::numeric_limits<std::uint32_t>::max()) {
            ++onLayer[1];
          } else {
            ++onLayer[goodOf[by] == majority.particle ? 2 : 3];
          }
        }
      }
      out << std::format("W {} {} {} {} {} {} {} {} {} {:.4f} {} {} {} {} {} {:.5g} {:.5g} {:.5g} {:.5g} {:.5g} {:.5g} {} {} {} {:.4g}\n", track, position,
                         roadBegin, roadEnd, static_cast<int>(at.layer), part, at.count == 0U ? 1 : 0, particle, wrongOrigin, static_cast<double>(wrongPt),
                         at.count != 0U && withGood.contains(particle) ? 1 : 0, onLayer[0], onLayer[1], onLayer[2], onLayer[3],
                         std::sqrt(static_cast<double>(wrongAt.varRPhi)), std::sqrt(static_cast<double>(wrongAt.varZ)), nearest[0], nearest[1], nearest[2],
                         nearest[3], pdgs[0], pdgs[1], pdgs[2], static_cast<double>(wrongRadius));
    }
  }
}

/// The findable physical primaries and their good tracks against \a truth,
/// for the clone study: one `P` line per
/// particle (id, pT, eta, phi, PDG, vz, layer bits, chain, on a disc, good
/// tracks), a `T` line per good track (iteration, rows) and, for a particle
/// with two good tracks or more, an `R` line per row carrying it (row,
/// layer, r, z, labels).
template <DetectorLayout Det>
void writeClones(const std::string& path, const typename TrackRefit<Det>::Tracks& tracks, const study::TrackingTruth& truth,
                 const study::ComparisonSettings& settings)
{
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot write " + path);
  }
  const auto rows = trackRows(tracks);
  const auto iterations = tracks.view().iteration;
  std::unordered_map<std::uint64_t, std::vector<std::size_t>> goodTracks;
  for (std::size_t track = 0; track < rows.size(); ++track) {
    const auto particle = study::goodParticle(truth.rows, rows[track]);
    if (particle != study::kNoTrackParticle) {
      goodTracks[particle].push_back(track);
    }
  }
  auto rowsOf = study::rowsOfParticles(truth.rows);
  for (const auto& [id, layers] : study::particleLayers(truth.rows, settings)) {
    const auto found = truth.particles.find(id);
    if (found == truth.particles.end()) {
      continue;
    }
    const auto& particle = found->second;
    if (particle.origin != study::ParticleOrigin::kPhysicalPrimary || !study::isFindable(particle.pt, layers, settings)) {
      continue;
    }
    const auto& good = goodTracks[id];
    const bool disc = (layers.layers & ~settings.barrelLayers) != 0U;
    out << std::format("P {} {:.5g} {:.5g} {:.5g} {} {:.5g} {} {} {} {}\n", id, particle.pt, particle.eta, particle.phi, particle.pdg, particle.vz,
                       layers.layers, layers.chain, disc ? 1 : 0, good.size());
    for (const auto track : good) {
      out << "T " << static_cast<int>(iterations[track]);
      for (const auto row : rows[track]) {
        out << ' ' << row;
      }
      out << '\n';
    }
    if (good.size() < 2U) {
      continue;
    }
    for (const auto row : rowsOf[id]) {
      const auto& r = truth.rows[row];
      out << std::format("R {} {} {:.4f} {:.4f} {}\n", row, static_cast<int>(r.layer), r.r, r.z, static_cast<int>(r.count));
    }
  }
}

/// Every row of the frame against \a truth, for occupancy studies:
/// one line per row -- its layer, global x, y and z [cm], its
/// particles' count, then its first particle's origin (0 primary, 1 decay,
/// 2 material; 3 for noise), pT [GeV], eta, PDG code and production vertex
/// x, y, z [cm] (zeros for noise), whether a good track of that particle
/// exists in \a tracks, and whether a track of \a tracks holds the row.
template <DetectorLayout Det>
void writeRowsOut(const std::string& path, const typename TrackRefit<Det>::Tracks& tracks, const study::TrackingTruth& truth, const TimeFrame<Det>& frame)
{
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot write " + path);
  }
  std::unordered_set<std::uint64_t> withGood;
  std::vector<bool> held(truth.rows.size());
  for (const auto& rows : trackRows(tracks)) {
    withGood.insert(study::goodParticle(truth.rows, rows));
    for (const auto row : rows) {
      held[row] = true;
    }
  }
  for (std::uint32_t row = 0; row < truth.rows.size(); ++row) {
    const auto& at = truth.rows[row];
    const auto measurement = frame.measurement(row);
    const auto cluster = static_cast<std::uint64_t>(frame.measurements().clusterIndex[row]);
    out << std::format("{} {:.4f} {:.4f} {:.4f} {}", static_cast<int>(at.layer), static_cast<double>(measurement.x), static_cast<double>(measurement.y),
                       static_cast<double>(measurement.z), static_cast<int>(at.count));
    const auto found = at.count == 0U ? truth.particles.end() : truth.particles.find(at.particles[0]);
    if (found == truth.particles.end()) {
      out << std::format(" 3 0 0 0 0 0 0 0 {} {}\n", held[row] ? 1 : 0, cluster);
      continue;
    }
    const auto& particle = found->second;
    out << std::format(" {} {:.4f} {:.4f} {} {:.4g} {:.4g} {:.4g} {} {} {}\n", static_cast<int>(particle.origin), static_cast<double>(particle.pt),
                       static_cast<double>(particle.eta), particle.pdg, static_cast<double>(particle.vx), static_cast<double>(particle.vy),
                       static_cast<double>(particle.vz), withGood.contains(at.particles[0]) ? 1 : 0, held[row] ? 1 : 0, cluster);
  }
}

/// The clone join: each findable physical primary with two good
/// tracks or more, its tracks in the order of their outermost rows' distance
/// from the origin; each one's outward state updated with the next one's rows
/// in order, along the truth. At each step the true row's chi2, and the
/// competition: every other row on that layer that the step's own predicted
/// window and the inner track's time box admit, whose update passes the gate.
/// One `J` line per join (particle, pT, eta, the two tracks' iterations,
/// their counts) and an `S` line per step (layer, the true row's chi2, 1 if
/// the window holds it, the wrong rows passing, the smallest wrong chi2, the
/// gate) into \a path; every row a track holds at the end is counted, free
/// or not.
template <DetectorLayout Det, typename Material, typename ParametersOf>
void writeCloneJoins(const std::string& path, const typename TrackRefit<Det>::Tracks& tracks, const study::TrackingTruth& truth,
                     const study::ComparisonSettings& settings, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material,
                     const Propagator& propagator, ParametersOf parametersOf)
{
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot write " + path);
  }
  const auto& measurements = frame.measurements();
  // The widest row's box per layer, which bounds the time query.
  std::array<std::uint32_t, static_cast<std::size_t>(Detector<Det>::nLayers)> widest{};
  for (std::uint32_t row = 0; row < truth.rows.size(); ++row) {
    const auto box = timeBoxAt(measurements, row);
    auto& width = widest[truth.rows[row].layer];
    width = std::max(width, box.upper() - box.lower());
  }
  const auto view = tracks.view();
  const auto rows = trackRows(tracks);
  std::unordered_map<std::uint64_t, std::vector<std::size_t>> goodTracks;
  for (std::size_t track = 0; track < rows.size(); ++track) {
    const auto particle = study::goodParticle(truth.rows, rows[track]);
    if (particle != study::kNoTrackParticle && view.status[track] == PropagationStatus::kOk) {
      goodTracks[particle].push_back(track);
    }
  }
  const auto distance = [&](std::uint32_t row) {
    const auto& at = truth.rows[row];
    return (at.r * at.r) + (at.z * at.z);
  };
  for (const auto& [id, layers] : study::particleLayers(truth.rows, settings)) {
    const auto found = truth.particles.find(id);
    const auto good = goodTracks.find(id);
    if (found == truth.particles.end() || good == goodTracks.end() || good->second.size() < 2U) {
      continue;
    }
    const auto& particle = found->second;
    if (particle.origin != study::ParticleOrigin::kPhysicalPrimary || !study::isFindable(particle.pt, layers, settings)) {
      continue;
    }
    auto ordered = good->second;
    const auto outermost = [&](std::size_t track) {
      float most = 0.F;
      for (const auto row : rows[track]) {
        most = std::max(most, distance(row));
      }
      return most;
    };
    std::ranges::sort(ordered, [&](std::size_t a, std::size_t b) { return outermost(a) < outermost(b); });
    for (std::size_t at = 0; at + 1U < ordered.size(); ++at) {
      const auto inner = ordered[at];
      const auto outer = ordered[at + 1U];
      std::vector<std::uint32_t> next(rows[outer].begin(), rows[outer].end());
      if (next.empty()) {
        continue;
      }
      std::ranges::sort(next, [&](std::uint32_t a, std::uint32_t b) { return distance(a) < distance(b); });
      out << std::format("J {} {:.5g} {:.5g} {} {} {} {} {}\n", id, particle.pt, particle.eta, static_cast<int>(view.iteration[inner]),
                         static_cast<int>(view.iteration[outer]), rows[inner].size(), rows[outer].size(), distance(next.front()) > outermost(inner) ? 1 : 0);
      const auto& parameters = parametersOf(view.iteration[inner]);
      BoundTrackState state = view.outwardState[inner];
      CovarianceFactors factors = factorCovariance(state.cov);
      auto box = trackTimeBox(measurements, view.clusters[inner], view.count[inner]);
      for (const auto row : next) {
        const LayerID target = truth.rows[row].layer;
        const auto prediction = predictLayerCrossing(propagator, parameters, geometry.planeOf(state), state, target, PropagationDirection::kForward);
        if (prediction.status != PropagationStatus::kOk) {
          out << std::format("S {} nopred\n", static_cast<int>(target));
          break;
        }
        CandidateWindow window{.earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(widest[target]) + 1,
                               .latestPS = static_cast<std::int64_t>(box.upper()) - 1,
                               .afterPS = box.lower(),
                               .phiCentre = prediction.phi,
                               .phiHalfWidth = std::numbers::pi_v<Scalar>,
                               .secondCentre = prediction.second,
                               .secondHalfWidth = 1.e4F};
        follow::predictedWindow<false>(propagator, parameters, geometry, material, state, factors, target, PropagationDirection::kForward, prediction, window,
                                       1.F);
        bool inWindow = false;
        std::uint32_t wrongPassing = 0U;
        Scalar wrongBest = std::numeric_limits<Scalar>::infinity();
        frame.grid().visitCandidates(target, window, frame.grid().geometry().allBlocks(), [&](std::uint32_t candidate) {
          const auto point = frame.measurement(candidate);
          if (!spatiallyCompatible(window, point.phi, secondCoordinate<Det>(target, point.r, point.z))) {
            return;
          }
          if (candidate == row) {
            inWindow = true;
            return;
          }
          BoundTrackState updated;
          CovarianceFactors updatedFactors;
          Scalar chi2 = 0.F;
          if (follow::updateCandidate<false>(propagator, geometry, material, measurementAt(frame, candidate), PropagationDirection::kForward, state, factors,
                                             updated, updatedFactors, chi2, 1.F) &&
              chi2 <= parameters.clusterChi2Max) {
            ++wrongPassing;
            wrongBest = std::min(wrongBest, chi2);
          }
        });
        BoundTrackState updated;
        CovarianceFactors updatedFactors;
        Scalar chi2 = 0.F;
        if (!follow::updateCandidate<false>(propagator, geometry, material, measurementAt(frame, row), PropagationDirection::kForward, state, factors, updated,
                                            updatedFactors, chi2, 1.F)) {
          out << std::format("S {} noupdate {} {}\n", static_cast<int>(target), inWindow ? 1 : 0, wrongPassing);
          break;
        }
        out << std::format("S {} {:.4g} {} {} {:.4g} {:.4g}\n", static_cast<int>(target), chi2, inWindow ? 1 : 0, wrongPassing, wrongBest,
                           parameters.clusterChi2Max);
        state = updated;
        factors = updatedFactors;
        box += timeBoxAt(measurements, row);
      }
    }
  }
}

/// A block's width [ps], for printing.
template <DetectorLayout Det>
double blockWidthPS(const GridGeometry<Det>& geometry)
{
  return static_cast<double>(geometry.blockWidth()) / static_cast<double>(geometry.partsPerPS());
}

/// Each row's layer, from its sensor in \a detector.
template <DetectorLayout Det>
std::vector<LayerID> rowLayers(const study::TrackingInput& input, const Geometry<Det>& detector)
{
  std::vector<LayerID> layers(input.columns.surface.size());
  for (std::size_t row = 0; row < layers.size(); ++row) {
    layers[row] = detector.sensor(input.columns.surface[row]).layer;
  }
  return layers;
}

/// Whether \a input's rows are binned by \a grid: as many buckets, and every
/// row in the bucket \a grid puts it in. A dump does not say which grid it
/// was binned with; its rows do.
template <DetectorLayout Det>
bool binnedBy(const study::TrackingInput& input, std::span<const LayerID> layers, const GridGeometry<Det>& grid)
{
  const auto& columns = input.columns;
  if (input.offsets.size() != std::size_t{grid.bucketCount()} + 1U || input.offsets.back() != layers.size()) {
    return false;
  }
  for (std::uint32_t bucket = 0; bucket < grid.bucketCount(); ++bucket) {
    for (auto row = input.offsets[bucket]; row < input.offsets[bucket + 1U]; ++row) {
      if (grid.bucketAt(layers[row], columns.timePS[row], columns.phi[row], secondCoordinate<Det>(layers[row], columns.r[row], columns.z[row])) != bucket) {
        return false;
      }
    }
  }
  return true;
}

/// \a input's rows, with their \a layers, as the producer would have binned
/// them by \a to: in cluster order, each into its bucket of \a to by its
/// layer, time, azimuth and second coordinate, by the producer's counting
/// sort, which keeps cluster order within a bucket. \a truth's rows follow theirs. Returns the seconds the binning
/// took.
template <DetectorLayout Det>
double rebin(study::TrackingInput& input, study::TrackingTruth* truth, std::span<const LayerID> layers, const GridGeometry<Det>& to)
{
  const auto start = std::chrono::steady_clock::now();
  auto& columns = input.columns;
  const std::size_t rows = columns.clusterIndex.size();
  std::vector<std::uint32_t> order(rows);
  std::iota(order.begin(), order.end(), 0U);
  std::ranges::sort(order, [&](std::uint32_t a, std::uint32_t b) { return columns.clusterIndex[a] < columns.clusterIndex[b]; });

  std::vector<std::uint32_t> offsets(std::size_t{to.bucketCount()} + 1U);
  MeasurementGridBuilder builder{std::span<std::uint32_t>{offsets}};
  std::vector<std::uint32_t> buckets(rows);
  for (const auto row : order) {
    const auto layer = layers[row];
    buckets[row] = to.bucketAt(layer, columns.timePS[row], columns.phi[row], secondCoordinate<Det>(layer, columns.r[row], columns.z[row]));
    builder.count(buckets[row]);
  }
  builder.prefixSum();
  std::vector<std::uint32_t> source(rows); // the old row of each new one
  for (const auto row : order) {
    source[builder.place(buckets[row])] = row;
  }
  builder.finish();

  const auto permute = [&](auto& column) {
    auto old = column;
    for (std::size_t row = 0; row < rows; ++row) {
      column[row] = old[source[row]];
    }
  };
  columns.apply([&](auto&... column) { (permute(column), ...); });
  if (truth != nullptr) {
    permute(truth->rows);
  }
  input.offsets = std::move(offsets);
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

/// The device's vertex list against the host's on the same frame:
/// both sorted by time and each frame's by
/// z, so paired by place -- how many bit for bit, how many with the same box
/// and contributors, and the largest z difference. A list of another length
/// pairs as far as both go, and says so.
inline std::string describeVertexAgreement(const VertexerResult& device, const VertexerResult& host)
{
  const std::size_t paired = std::min(device.vertices.size(), host.vertices.size());
  std::size_t identical = 0U;
  std::size_t sameClaim = 0U;
  double maxDzUM = 0.;
  const auto bits = [](Scalar value) { return std::bit_cast<std::uint32_t>(value); };
  for (std::size_t v = 0; v < paired; ++v) {
    const auto& d = device.vertices[v];
    const auto& h = host.vertices[v];
    const bool claim = d.timePS == h.timePS && d.widthPS == h.widthPS && device.contributors[v] == host.contributors[v];
    sameClaim += claim ? 1U : 0U;
    identical += claim && bits(d.z) == bits(h.z) && bits(d.varZ) == bits(h.varZ) && bits(d.varTransverse) == bits(h.varTransverse) ? 1U : 0U;
    maxDzUM = std::max(maxDzUM, std::abs(static_cast<double>(d.z) - static_cast<double>(h.z)) * 1.e4);
  }
  return std::format(
    "vertexer, device against host: {} vertices against {}{}; of {} paired, {} bit for bit, {} with the same box and contributors, |dz| "
    "at most {:.2f} um; candidates {} against {}",
    device.vertices.size(), host.vertices.size(), device.vertices.size() == host.vertices.size() ? "" : " (NOT the same count)", paired, identical, sameClaim,
    maxDzUM, device.candidates, host.candidates);
}

/// --seeding-steps \a text's pT steps [GeV]: one or more, comma separated,
/// each a whole number token and a positive momentum.
[[nodiscard]] inline std::vector<double> seedingSteps(const std::string& text)
{
  std::vector<double> steps;
  for (const auto part : text | std::views::split(',')) {
    const std::string token{part.begin(), part.end()};
    std::size_t parsed = 0U;
    double step = 0.;
    try {
      step = std::stod(token, &parsed);
    } catch (const std::logic_error&) {
      parsed = 0U;
    }
    if (parsed == 0U || parsed != token.size() || !std::isfinite(step) || step <= 0.) {
      throw std::runtime_error("--seeding-steps " + text + ": one or more steps, each a positive momentum [GeV]");
    }
    steps.push_back(step);
  }
  if (steps.empty()) {
    throw std::runtime_error("--seeding-steps " + text + ": one or more steps, each a positive momentum [GeV]");
  }
  return steps;
}

/// The truth-removal study: every
/// configured iteration's graph alone, in order, each on the rows the earlier
/// passes' seeded particles left; the report on stdout, a JSON line per pass
/// to --seeding-strategy. Each pass starts from a released pool -- the
/// timeframe ended and begun again, the held rows held again -- so that its
/// peak is its own. With --resume-state the checkpoint's iterations are the
/// real run's: its good tracks seed before
/// the first pass, every pass starts from its used mask, and the passes are
/// the configured iterations after it. 0 when every pass completed, 2 when
/// one failed (a pass over the pool's limit), which ends the study.
template <DetectorLayout Det>
int seedingStrategy(const Options& options, Controller<Det>& tracker, const TimeFrame<Det>& frame, const study::TrackingTruth& truth,
                    const study::ComparisonSettings& comparison, std::span<const IterationSettings> settings, const workflow::BeamLine& beamLine,
                    std::span<const VertexSeed> vertices, std::int32_t blocksPerSlice, BoundedMemoryResource& pool, const TrackingGeometry<Det>& tracking,
                    const o2::base::MatLayerCylSet& material, const Propagator& propagator)
{
  using Clock = std::chrono::steady_clock;
  const auto& param = recoParam<Det>();
  if (options.backend != "cpu" || param.vertexer) {
    throw std::runtime_error("--seeding-strategy runs the host tracker on the dump's vertices: not with a GPU backend or GTRKReco<DET>.vertexer");
  }
  auto steps = seedingSteps(options.seedingSteps);
  std::vector<std::uint8_t> shadow(settings.size(), std::uint8_t{0});
  for (const auto part : options.seedingShadow | std::views::split(',')) {
    if (part.empty()) {
      continue;
    }
    const auto pass = std::stoul(std::string{part.begin(), part.end()});
    if (pass + 1U >= settings.size()) {
      throw std::runtime_error("--seeding-shadow " + options.seedingShadow + ": a shadow pass needs a pass after it");
    }
    shadow[pass] = 1U;
  }
  study::SeedingPopulation population = study::SeedingPopulation::kPrimaries;
  if (options.seedingPopulation == "decay") {
    population = study::SeedingPopulation::kDecays;
  } else if (options.seedingPopulation == "v0") {
    population = study::SeedingPopulation::kV0s;
  } else if (options.seedingPopulation != "primary") {
    throw std::runtime_error("--seeding-population " + options.seedingPopulation + ": primary, decay or v0");
  }
  std::optional<typename Controller<Det>::Checkpoint> resumeFrom;
  if (!options.resumeState.empty()) {
    resumeFrom = readCheckpoint<typename Controller<Det>::Checkpoint>(options.resumeState, settings, checkpointFrame(frame),
                                                                      checkpointConditions(beamLine, vertices, blocksPerSlice, param.sliceHaloPS));
  }
  const std::size_t first = resumeFrom ? resumeFrom->summaries.size() : 0U;
  std::ofstream json(options.seedingStrategy);
  if (!json) {
    throw std::runtime_error("cannot write " + options.seedingStrategy);
  }
  tracker.setProfiling(false);
  const auto plan = tracker.planTimeframe(frame, vertices, blocksPerSlice, param.sliceHaloPS);
  if (tracker.beginTimeframe(frame) != StageStatus::kComplete) {
    throw std::runtime_error("--seeding-strategy: the timeframe's used mask does not fit the pool");
  }
  study::SeedingStrategy<Det> study(truth, comparison, std::move(steps), !options.seedingNoRemoval, options.seedingFit, population);
  std::cout << std::format("seeding strategy: {} passes, {}, scoring the findable {}\n", settings.size() - first,
                           options.seedingNoRemoval ? "each on the whole frame" : "each on the rows the earlier passes' seeded particles left",
                           study::kSeedingPopulationNames[static_cast<std::size_t>(population)]);
  if (resumeFrom) {
    tracker.restoreUsed(resumeFrom->used);
    study.startAfter(first, resumeFrom->tracks, resumeFrom->used, std::cout);
  }
  int status = 0;
  std::vector<std::uint32_t> held;
  for (std::size_t pass = first; pass < settings.size(); ++pass) {
    if (pass != first) {
      tracker.endTimeframe();
      if (tracker.beginTimeframe(frame) != StageStatus::kComplete) {
        throw std::runtime_error("--seeding-strategy: the timeframe's used mask does not fit the pool");
      }
      if (resumeFrom) {
        tracker.restoreUsed(resumeFrom->used);
      }
      tracker.holdRows(held);
    }
    study.beginPass(shadow[pass] != 0U);
    pool.resetPeakMemory();
    const auto workBefore = workCounts();
    double seconds = 0.;
    double fitSeconds = 0.;
    auto start = Clock::now();
    const auto ran = tracker.runGraph(pass, frame, plan, options.seedingFit, [&](const auto& result, const auto& roads, const auto* records, double fitted) {
      seconds += std::chrono::duration<double>(Clock::now() - start).count();
      fitSeconds += fitted;
      study.observe(result, roads, records);
      start = Clock::now();
    });
    if (ran != StageStatus::kComplete) {
      std::cout << std::format("pass {}: the graph failed with status {} (infeasible); the study ends here\n", pass, static_cast<int>(ran));
      status = 2;
      break;
    }
    const auto workAfter = workCounts();
    const auto counted = [&](WorkCounter counter) {
      const auto at = static_cast<std::size_t>(counter);
      return workAfter[at] - workBefore[at];
    };
    // The road fit stage's tests on one road (`roadFitRefusal`).
    const auto& parameters = tracker.iteration(pass);
    const typename study::SeedingStrategy<Det>::Refit refit = [&](std::span<const std::uint32_t> road) {
      return roadFitRefusal(parameters, propagator, tracking, material, frame, road);
    };
    const auto rows = study.endPass(settings[pass],
                                    {.seconds = seconds,
                                     .poolPeakBytes = pool.getPeakMemory(),
                                     .counted = kWorkCounters,
                                     .tripletFits = counted(WorkCounter::kTripletFits),
                                     .neighbourFits = counted(WorkCounter::kNeighbourFits),
                                     .fitSeconds = fitSeconds},
                                    std::cout, &json, options.seedingFit ? refit : typename study::SeedingStrategy<Det>::Refit{});
    if (!json) {
      throw std::runtime_error("cannot write " + options.seedingStrategy);
    }
    tracker.holdRows(rows);
    held.insert(held.end(), rows.begin(), rows.end());
  }
  tracker.endTimeframe();
  json.close();
  if (!json) {
    throw std::runtime_error("cannot write " + options.seedingStrategy);
  }
  if (!options.seedingParticles.empty()) {
    std::ofstream particles(options.seedingParticles);
    study.writeParticles(particles);
    particles.close();
    if (!particles) {
      throw std::runtime_error("cannot write " + options.seedingParticles);
    }
  }
  return status;
}

/// The truth studies that watch every iteration: the follow's rows per track for --fakes-out and the follow's
/// replay for --follow-steps and --follow-stops.
template <DetectorLayout Det>
struct IterationStudies {
  std::vector<std::array<std::uint8_t, 2>> followed;
  std::optional<FollowAnatomy<Det, o2::base::MatLayerCylSet>> follow;
  bool fakes = false;

  IterationStudies(const Options& options, const study::TrackingTruth* truth) : fakes(!options.fakesOut.empty())
  {
    if ((fakes || !options.followSteps.empty() || options.followStops) && truth == nullptr) {
      throw std::runtime_error("--fakes-out, --follow-steps and --follow-stops need the MC truth");
    }
    if (!options.followSteps.empty() || options.followStops) {
      follow.emplace(options.followSteps, *truth);
    }
  }
  [[nodiscard]] bool active() const { return fakes || follow.has_value(); }

  template <typename Observation>
  void observe(const Observation& iteration, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry,
               const o2::base::MatLayerCylSet& material, const Propagator& propagator)
  {
    if (follow) {
      follow->observe(iteration, parameters, frame, geometry, material, propagator);
    }
    if (fakes) {
      recordFollowed(iteration, followed);
    }
  }
};

/// The last run's truth studies the options ask for: the fakes,
/// the clones, the rows and the clone joins.
template <DetectorLayout Det, typename Material, typename ParametersOf, typename StagesOf, typename RecordChi2>
void writeTruthStudies(const Options& options, const typename TrackRefit<Det>::Tracks& tracks, const study::TrackingTruth& truth,
                       const study::ComparisonSettings& comparison, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& tracking,
                       const Material& material, const Propagator& propagator, const std::vector<std::array<std::uint8_t, 2>>& followed,
                       ParametersOf parametersOf, StagesOf stagesOf, RecordChi2 recordChi2)
{
  if (!options.fakesOut.empty()) {
    writeFakes<Det>(options.fakesOut, tracks, truth, comparison, frame, followed, stagesOf, recordChi2);
  }
  if (!options.clonesOut.empty()) {
    writeClones<Det>(options.clonesOut, tracks, truth, comparison);
  }
  if (!options.rowsOut.empty()) {
    writeRowsOut<Det>(options.rowsOut, tracks, truth, frame);
  }
  if (!options.cloneJoins.empty()) {
    writeCloneJoins<Det>(options.cloneJoins, tracks, truth, comparison, frame, tracking, material, propagator, parametersOf);
  }
}

/// Track the dumped input \a options.repeat times with the parameters already
/// loaded; 0 if every run gave the same digest.
template <DetectorLayout Det>
int track(const Options& options) // NOLINT(readability-function-size): the bench's driver, one mode per option
{
  using Frame = TimeFrame<Det>;
  using Controller = study::Controller<Det>;
  const auto& param = recoParam<Det>();

  const auto start = std::chrono::steady_clock::now();
  auto input = study::readTrackingInput(options.input);
  const auto readSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  std::cout << std::format("{} measurements, {} buckets, {} vertices read from {} in {:.3f} s\n", input.columns.clusterIndex.size(), input.offsets.size() - 1U,
                           input.vertices.size(), options.input, readSeconds)
            << std::format("the beam line the rows are about: {}\n", workflow::describe(input.beamLine));
  if (options.beamWidthCM > 0.) {
    input.beamLine.sigmaX = input.beamLine.sigmaY = options.beamWidthCM;
    std::cout << std::format("the beam's transverse width taken as {} cm (--beam-width)\n", options.beamWidthCM);
  }
  // The MC truth to score the tracks by.
  const auto truthPath = options.truth.empty() ? options.input + ".truth" : options.truth;
  std::optional<study::TrackingTruth> truth;
  if (!options.truth.empty() || std::filesystem::exists(truthPath)) {
    truth = study::readTrackingTruth(truthPath, input);
    std::cout << std::format("{} truth rows, {} particles read from {}\n", truth->rows.size(), truth->particles.size(), truthPath);
  }
  // The geometry first: a placed layout's extents, which the grid and the
  // windows below are sized from, are the loaded geometry's (Base/LayerEnvelopes.h);
  // its measurement depths the dump's, where its rows were placed.
  // The conditions, as a standalone O2 program reads them: from the CCDB at
  // --ccdb-url (trackBench.cxx), at the paths the sink takes them from.
  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  if (ccdb.get<TGeoManager>("GLO/Config/GeometryAligned") == nullptr) {
    throw std::runtime_error("no GLO/Config/GeometryAligned at " + options.ccdbUrl);
  }
  // The run's orbits per timeframe, as the measurement producer takes them
  // (Conditions::orbitsPerTimeframe): its GRPECS's; ALICE 3's simulation
  // writes none, and the dump's HBFUtils is its.
  int orbitsPerTimeframe = 0;
  if constexpr (std::same_as<Det, detectors::TRKFT3>) {
    orbitsPerTimeframe = o2::raw::HBFUtils::Instance().getNOrbitsPerTF();
  } else {
    const auto* ecs = ccdb.get<o2::parameters::GRPECSObject>("GLO/Config/GRPECS");
    if (ecs == nullptr) {
      throw std::runtime_error("no GLO/Config/GRPECS at " + options.ccdbUrl);
    }
    orbitsPerTimeframe = ecs->getNHBFPerTF();
  }
  // The grid this build derives from the dump's parameters and the trial's
  // keys below; the rows are rebinned into it unless they already are.
  // Each row's layer is its sensor's, so
  // a dump binned by a build with another link table reads all the same.
  const auto frames = study::layerFrames<Det>(input);
  // A tuning trial's configuration.
  if (!options.keyValues.empty()) {
    o2::conf::ConfigurableParam::updateFromString(options.keyValues);
  }
  if (options.threads >= 0) {
    o2::conf::ConfigurableParam::setValue(std::string{RecoParamOf<Det>::type::Instance().getName()} + ".threads", std::to_string(options.threads));
  }
  workflow::checkRecoParam<Det>(param);
  workflow::TrackingInterface<Det> interface{!param.vertexer};
  interface.initGeometry(input.depths);
  const auto geometry = workflow::makeGridGeometry<Det>(param, input.beamLine, frames, orbitsPerTimeframe);
  const auto* material = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb.get<o2::base::MatLayerCylSet>(workflow::materialLutPath<Det>()));
  const bool device = options.backend != "cpu";
  if (device) {
    interface.initDevice({.backend = options.backend, .device = options.gpuDevice, .arenaMB = options.gpuArenaMB});
  }
  interface.configure({workflow::nominalBzKiloGauss(), input.depths, material}, {input.beamLine, geometry});
  auto& tracker = interface.tracker();
  const auto& pool = interface.pool();
  auto* vertexing = interface.vertexing();
  const auto& tracking = interface.tracking();
  const auto& scattering = interface.scattering();
  const auto propagator = interface.propagator();
  double rebinSeconds = 0.;
  const auto layers = rowLayers<Det>(input, interface.geometry());
  const bool rebinned = options.rebin || !binnedBy<Det>(input, layers, geometry);
  if (rebinned) {
    rebinSeconds = rebin<Det>(input, truth ? &*truth : nullptr, layers, geometry);
  }
  std::cout << std::format("grid{}: {} blocks of {:.0f} ps; phi x second cells by layer:", rebinned ? std::format(", rebinned in {:.3f} s", rebinSeconds) : "",
                           geometry.blockCount(), blockWidthPS(geometry));
  for (std::size_t layer = 0; layer < GridGeometry<Det>::kLayerCount; ++layer) {
    const auto& axes = geometry.axes(static_cast<LayerID>(layer));
    std::cout << std::format(" {}x{}", axes.nPhi, axes.nSecond);
  }
  std::cout << std::format("; {} buckets, offsets {:.1f} MB\n", geometry.bucketCount(),
                           static_cast<double>(input.offsets.size() * sizeof(std::uint32_t)) / 1.e6);
  const Frame frame{measurementView(input.columns), geometry, std::span<const typename Frame::Index>{input.offsets}};
  if (const auto* why = interface.checkInput(input.columns, input.offsets, geometry.bucketCount(), input.vertices)) {
    throw std::runtime_error(why);
  }
  std::string phiCells;
  std::string secondCells;
  for (std::size_t layer = 0; layer < GridGeometry<Det>::kLayerCount; ++layer) {
    const auto& axes = geometry.axes(static_cast<LayerID>(layer));
    phiCells += std::format("{}{}", layer == 0U ? "" : ",", axes.nPhi);
    secondCells += std::format("{}{}", layer == 0U ? "" : ",", axes.nSecond);
  }
  const auto gridJson = std::format(
    R"("grid":{{"blocks":{},"blockWidthPS":{:.3f},"phi":[{}],"second":[{}],"buckets":{},"offsetBytes":{},"rebinSeconds":{:.6f}}},)", geometry.blockCount(),
    blockWidthPS(geometry), phiCells, secondCells, geometry.bucketCount(), input.offsets.size() * sizeof(std::uint32_t), rebinSeconds);
  auto comparison = study::comparisonSettings<Det>(4, options.findableMinPtGeV);
  comparison.allLayersMask = options.findableAllLayers;
  if (!options.seedCoverage.empty()) {
    if (!truth) {
      throw std::runtime_error("--seed-coverage needs the MC truth");
    }
    std::vector<std::uint64_t> masks;
    for (const auto part : options.seedCoverage | std::views::split(',')) {
      masks.push_back(std::stoull(std::string{part.begin(), part.end()}));
    }
    printSeedCoverage(*truth, comparison, masks);
  }

  std::vector<IterationSettings> settings;
  for (int iteration = 0; iteration < workflow::iterationCount(param); ++iteration) {
    settings.push_back(workflow::iterationSettings(param, iteration));
    if (settings.back().vertexSource == VertexSource::kVertexList && input.vertices.empty() && !param.vertexer) {
      throw std::runtime_error(std::format("iteration {} seeds from the vertex list, and {} holds none", iteration, options.input));
    }
  }
  const auto chunkSize = static_cast<typename Controller::Index>(std::max(param.trackletChunkSize, 1));
  tracker.setProfiling(!options.noProfile);
  if (options.dram && !profile::countDram()) {
    throw std::runtime_error("--dram: the memory controllers' counters (amd_umc) cannot be opened");
  }
  if (device) {
    if (options.funnel || options.digestStages || !options.fitInputs.empty() || !options.trace.empty()) {
      throw std::runtime_error("--funnel, --digest-stages, --fit-inputs and --trace observe the host tracker: not with a GPU backend");
    }
  } else if (options.hostReference) {
    throw std::runtime_error("--host-reference compares a GPU backend's tracks with the host's: it needs --backend cuda or hip");
  }
  // --fit-inputs: the last run's, with each iteration's parameters the
  // tracker's own.
  study::FitInputs fitInputs;
  bool captureFits = false;
  // --funnel: the last run's stage funnel of the findable.
  if (options.funnel && !truth) {
    throw std::runtime_error("--funnel needs the MC truth");
  }
  std::optional<study::StageFunnel<Det, Controller>> funnel;
  // --vertex-study: the vertexer on one iteration's triplets against the
  // truth, and its replay (VertexStudy.h).
  if (options.vertexStudy >= 0 && (!truth || device || std::cmp_greater_equal(options.vertexStudy, settings.size()))) {
    throw std::runtime_error("--vertex-study needs the MC truth, the host tracker and one of the configured iterations");
  }
  std::optional<study::VertexStudy<Det, Controller>> vertexStudy;
  if (options.verticesOnly && (!truth || device || !param.vertexer)) {
    throw std::runtime_error("--vertices-only needs the MC truth, the host and GTRKReco<DET>.vertexer");
  }
  if (options.digestStages || !options.fitInputs.empty() || options.funnel || options.vertexStudy >= 0) {
    tracker.setSliceObserver([&](const Controller& observed, const typename Controller::SliceObservation& slice) {
      if (options.digestStages) {
        printStageDigests<Det>(slice);
      }
      if (captureFits) {
        collectFitInputs<Det>(slice, observed.iteration(slice.iteration), fitInputs);
      }
      if (funnel) {
        funnel->observeSlice(observed, slice);
      }
      if (vertexStudy) {
        vertexStudy->observeSlice(observed, slice);
      }
    });
  }
  // Each iteration's work counts: the counts at its end less those at the
  // previous one's.
  std::vector<WorkCounts> work;
  WorkCounts workBefore{};
  std::vector<LinkWorkCounts> linkWork;
  LinkWorkCounts linkWorkBefore;
  if ((options.completeAnatomy && !truth) || (options.completeAnatomy && options.findableAllLayers == 0U)) {
    throw std::runtime_error("--complete-anatomy needs the MC truth and --findable-all-layers");
  }
  if (!options.endEdits.empty() && options.endDecisions.empty()) {
    throw std::runtime_error("--end-edits needs --end-decisions");
  }
  std::optional<study::CompleteAnatomy<Det>> completeAnatomy;
  if (options.completeAnatomy) {
    std::vector<std::uint64_t> roadLayers;
    roadLayers.reserve(settings.size());
    for (const auto& iteration : settings) {
      roadLayers.push_back(iteration.layerMask);
    }
    completeAnatomy.emplace(*truth, options.findableAllLayers, options.findableMinPtGeV, std::move(roadLayers),
                            typename study::CompleteAnatomy<Det>::ClassCuts{.pdg = options.anatomyPdg,
                                                                            .findable = options.anatomyFindable ? &comparison : nullptr});
    std::vector<double> floors;
    floors.reserve(settings.size());
    for (const auto& iteration : settings) {
      floors.push_back(iteration.windows.minPtGeV);
    }
    completeAnatomy->setIterationFloors(std::move(floors));
    completeAnatomy->setBeamLine(vertexSeed(workflow::beamSpotConstraint(input.beamLine, param)));
    if (!options.anatomyBeamLine.empty()) {
      completeAnatomy->keepBeamLine();
    }
    if (!options.anatomyFollowLoss.empty()) {
      completeAnatomy->keepFollowLoss();
    }
    if (!options.stepAnatomy.empty() || !options.hits.empty()) {
      completeAnatomy->open(options.stepAnatomy, options.hits, options.hits.empty() ? std::string{} : options.ft3Hits);
    }
    if (!options.truthWalk.empty()) {
      if (options.hits.empty()) {
        throw std::runtime_error("--truth-walk needs --hits");
      }
      completeAnatomy->openTruthWalk(options.truthWalk);
    }
    if (options.truthWalkIteration >= 0) {
      completeAnatomy->setWalkIteration(static_cast<std::size_t>(options.truthWalkIteration));
    }
    if (!options.inwardHistory.empty()) {
      if (options.hits.empty()) {
        throw std::runtime_error("--inward-history needs --hits");
      }
      completeAnatomy->openInwardHistory(options.inwardHistory);
    }
    if (!options.endDecisions.empty()) {
      completeAnatomy->openEndDecisions(options.endDecisions);
    }
    if (!options.endEdits.empty()) {
      completeAnatomy->openEndEdits(options.endEdits);
    }
    if (!options.kine.empty()) {
      if (options.hits.empty()) {
        throw std::runtime_error("--kine needs --hits");
      }
      completeAnatomy->openKine(options.kine);
    }
  }
  IterationStudies<Det> studies(options, truth ? &*truth : nullptr);
  // The real run's iterations as the truth-removal
  // study's passes, scored by their good tracks and the used mask.
  // One study per run: a repetition scores its own iterations from the start,
  // its JSON written anew, not on from the previous run's study; the last
  // run's is the one left.
  std::optional<study::SeedingStrategy<Det>> seedingReal;
  std::ofstream seedingRealJson;
  std::size_t seedingRealNext = 0U; ///< the iteration after the last one scored
  const auto beginSeedingReal = [&]() {
    seedingReal.emplace(*truth, comparison, seedingSteps(options.seedingSteps), false, false);
    seedingRealJson.close();
    seedingRealJson.open(options.seedingReal);
    if (!seedingRealJson) {
      throw std::runtime_error("cannot write " + options.seedingReal);
    }
  };
  if (!options.seedingReal.empty()) {
    if (!truth) {
      throw std::runtime_error("--seeding-real needs the MC truth");
    }
    beginSeedingReal();
  }
  if (kWorkCounters || options.funnel || completeAnatomy || studies.active() || seedingReal) {
    tracker.setIterationObserver([&](const Controller& observed, const typename Controller::IterationObservation& iteration) {
      if (seedingReal) {
        // A run's iterations come in order: an earlier one begins the next run.
        if (iteration.iteration < seedingRealNext) {
          beginSeedingReal();
        }
        seedingRealNext = iteration.iteration + 1U;
        const auto& summary = *iteration.summary;
        seedingReal->observeReal(iteration.tracks->view(), iteration.trackOffset, iteration.tracks->size(), iteration.used, summary.roads);
        double seconds = 0.;
        for (const auto phase :
             {profile::Phase::kSlices, profile::Phase::kRoads, profile::Phase::kSelection, profile::Phase::kRefit, profile::Phase::kFollow}) {
          seconds += summary.seconds(phase);
        }
        static_cast<void>(seedingReal->endPass(settings[iteration.iteration], {.seconds = seconds}, std::cout, &seedingRealJson));
      }
      studies.observe(iteration, observed.iteration(iteration.iteration), frame, tracking, *material, propagator);
      if (completeAnatomy) {
        completeAnatomy->observe(iteration, observed.iteration(iteration.iteration), frame, tracking, *material, propagator);
      }
      if constexpr (kWorkCounters) {
        const auto now = workCounts();
        WorkCounts& counts = work.emplace_back();
        for (std::size_t counter = 0; counter < kWorkCounterCount; ++counter) {
          counts[counter] = now[counter] - workBefore[counter];
        }
        workBefore = now;
        auto links = linkWorkCounts();
        LinkWorkCounts& linkCounts = linkWork.emplace_back(links);
        for (std::size_t link = 0; link < std::min(links.size(), linkWorkBefore.size()); ++link) {
          for (std::size_t counter = 0; counter < kLinkCounterCount; ++counter) {
            linkCounts[link][counter] -= linkWorkBefore[link][counter];
          }
        }
        linkWorkBefore = std::move(links);
      }
      if (funnel) {
        funnel->observeIteration(observed, iteration);
      }
    });
  }
  const int threads = detail::TaskArena{std::max(param.threads, 0)}.maxConcurrency();
  const auto label = workflow::sourceLabel(settings);
  std::cout << std::format("{}, {} iterations, {} arena threads, chunks of {}{}\n", label, workflow::iterationCount(param), threads, chunkSize,
                           kWorkCounters ? "; work counters on: the times are not the default build's" : "");

  std::ofstream json;
  if (!options.json.empty()) {
    json.open(options.json, std::ios::app);
    if (!json) {
      throw std::runtime_error("cannot append to " + options.json);
    }
  }
  profile::Trace trace{static_cast<std::size_t>(threads)};
  // Against the truth, under the sink's rules and in its lines;
  // outside the time.
  const auto scoreAndPrint = [&](const study::TrackRows& rows) -> std::optional<std::pair<study::TrackScore, double>> {
    if (!truth) {
      return std::nullopt;
    }
    const auto scoreStart = std::chrono::steady_clock::now();
    auto score = study::scoreTracks(truth->rows, truth->particles, rows, comparison);
    const double scoreSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - scoreStart).count();
    std::cout << std::format("score against the truth, in {:.3f} s:\n", scoreSeconds);
    for (const auto& line : study::describeScore(score, comparison, "GTRK")) {
      std::cout << line << '\n';
    }
    return std::pair{std::move(score), scoreSeconds};
  };

  std::vector<double> walls;
  std::vector<std::uint64_t> digests;
  std::vector<VertexSeed> foundVertices;
  const auto blocksPerSlice = geometry.blocksPerSlice(param.sliceWidthPS);
  if (!options.seedingStrategy.empty()) {
    if (!truth) {
      throw std::runtime_error("--seeding-strategy needs the MC truth");
    }
    return seedingStrategy<Det>(options, tracker, frame, *truth, comparison, settings, input.beamLine, input.vertices, blocksPerSlice, *pool, tracking,
                                *material, propagator);
  }
  const bool deviceVertexing = vertexing && device;
  bool referenceSame = true;
  // A checkpoint: saved after the first iterations, or
  // the run going on from one -- on the host, and from one without the
  // observers, which would miss its iterations.
  if (!options.saveState.empty() && (options.stateIterations == 0U || options.stateIterations > settings.size() || !options.resumeState.empty())) {
    throw std::runtime_error("--save-state needs --state-iterations between 1 and the configuration's iterations, and no --resume-state");
  }
  const auto checkpointFingerprints = [&] {
    return std::pair{checkpointFrame(frame), checkpointConditions(input.beamLine, input.vertices, blocksPerSlice, param.sliceHaloPS)};
  };
  std::optional<typename Controller::Checkpoint> resumeFrom;
  if (!options.resumeState.empty()) {
    if (device || param.vertexer || kWorkCounters || options.funnel || completeAnatomy || studies.active() || seedingReal || options.digestStages ||
        !options.fitInputs.empty() || options.vertexStudy >= 0) {
      throw std::runtime_error("--resume-state runs the host tracker on the dump's vertices without observers");
    }
    const auto [frameFingerprint, conditions] = checkpointFingerprints();
    resumeFrom = readCheckpoint<typename Controller::Checkpoint>(options.resumeState, settings, frameFingerprint, conditions);
    std::cout << std::format("resuming from {}: {} iterations, {} tracks\n", options.resumeState, resumeFrom->summaries.size(), resumeFrom->tracks.size());
  }
  // Clusters taken before the run, by their cluster
  // number -- rows move with the grid -- as a checkpoint of no iterations
  // whose used mask holds their rows, so the observers miss nothing.
  if (!options.holdClusters.empty()) {
    if (device || param.vertexer || resumeFrom) {
      throw std::runtime_error("--hold-clusters runs the host tracker, without GTRKReco<DET>.vertexer or --resume-state");
    }
    std::ifstream in(options.holdClusters);
    if (!in) {
      throw std::runtime_error("cannot read " + options.holdClusters);
    }
    const auto columns = frame.measurements();
    std::unordered_map<std::uint64_t, std::uint32_t> rowOf;
    rowOf.reserve(frame.size());
    for (std::uint32_t row = 0; row < frame.size(); ++row) {
      if (!rowOf.emplace(columns.clusterIndex[row], row).second) {
        throw std::runtime_error(std::format("--hold-clusters: cluster {} on two rows", static_cast<std::uint64_t>(columns.clusterIndex[row])));
      }
    }
    resumeFrom.emplace();
    resumeFrom->used.assign(frame.size(), std::uint8_t{kFree});
    std::size_t held = 0U;
    for (std::uint64_t cluster = 0; in >> cluster;) {
      const auto found = rowOf.find(cluster);
      if (found == rowOf.end()) {
        throw std::runtime_error(std::format("--hold-clusters: no cluster {} in the frame", cluster));
      }
      held += resumeFrom->used[found->second] == kFree ? 1U : 0U;
      resumeFrom->used[found->second] = kExclusive;
    }
    std::cout << std::format("holding {} clusters of {} from {}\n", held, frame.size(), options.holdClusters);
  }
  if (!options.saveState.empty() && (device || param.vertexer)) {
    throw std::runtime_error("--save-state runs the host tracker on the dump's vertices");
  }
  for (unsigned repetition = 0; repetition < std::max(options.repeat, 1U); ++repetition) {
    const bool last = repetition + 1U == std::max(options.repeat, 1U);
    const auto timeframe = interface.endTimeframeOnExit();
    trace.clear();
    tracker.setTrace(!options.trace.empty() && last ? &trace : nullptr);
    captureFits = !options.fitInputs.empty() && last;
    funnel.reset();
    vertexStudy.reset();
    if (options.vertexStudy >= 0 && last) {
      const auto iteration = static_cast<std::size_t>(options.vertexStudy);
      vertexStudy.emplace(frame, truth->rows, truth->particles, iteration, settings[iteration].windows.minPtGeV, comparison,
                          workflow::vertexerSettings(param),
                          typename study::VertexStudy<Det, Controller>::Settings{.windowUM = param.vertexerWindowUM,
                                                                                 .claimUM = param.vertexerClaimUM,
                                                                                 .minContributors = static_cast<std::uint32_t>(param.vertexerMinContributors),
                                                                                 .significance = param.vertexerSignificance,
                                                                                 .maxChi2RZ = param.vertexerMaxChi2RZ,
                                                                                 .maxDcaUM = param.vertexerMaxDcaCM * 1.e4},
                          threads);
      vertexStudy->setBeamLine(input.beamLine.x, input.beamLine.y);
      if (!options.vertexPrior.empty()) {
        vertexStudy->setPriorVertices(options.vertexPrior);
      }
    }
    if (options.funnel && last) {
      std::vector<double> minPtGeV;
      minPtGeV.reserve(settings.size());
      for (const auto& iteration : settings) {
        minPtGeV.push_back(iteration.windows.minPtGeV);
      }
      funnel.emplace(frame, truth->rows, truth->particles, comparison, vertexSeed(workflow::beamSpotConstraint(input.beamLine, param)), std::move(minPtGeV));
      funnel->setRoadFit(tracking, *material, propagator);
      funnel->keepReportInputs(!options.fakesOut.empty());
    }
    // --vertices-only keeps the first run's claims, which its score takes.
    const auto vertexed = interface.vertices(frame, input.vertices, options.verticesOnly && repetition == 0U);
    if (vertexed.dropped) {
      throw std::runtime_error("the vertexing pass ran out of memory");
    }
    const auto vertices = vertexed.seeds;
    const double vertexSeconds = !device ? vertexed.seconds : 0.;
    if (vertexing && !device) {
      const auto& found = *vertexed.found;
      std::cout << std::format("vertexer: {} vertices from {} candidates of {} triplets in {:.4f} s\n", found.vertices.size(), found.candidates,
                               found.triplets, vertexSeconds);
      std::cout << std::format("digest, vertices: {:016x}\n", vertexDigest(found.vertices, found.contributors));
    }
    if (options.verticesOnly) {
      // The pass alone: every run timed, the first scored against the MC
      // collisions (VertexQA.h), as the sink's vertex QA line.
      const auto& found = *vertexed.found;
      const auto vertexListDigest = vertexDigest(found.vertices, found.contributors);
      walls.push_back(vertexSeconds);
      digests.push_back(vertexListDigest);
      std::optional<VertexScore> score;
      if (repetition == 0U) {
        score = scoreVertices<Det>(truth->rows, truth->particles, vertexed.pass->claims, found);
        std::cout << describe(*score) << '\n';
      }
      if (json.is_open()) {
        json << std::format(R"({{"input":"{}","run":{},"threads":{},"planSeconds":0,"runSeconds":{:.9f},"maxRssBytes":{},"digest":"{:016x}","vertices":{})",
                            options.input, repetition, threads, vertexSeconds, maxRssBytes(), vertexListDigest, found.vertices.size());
        if (score) {
          json << std::format(R"(,"score":{{"vertices":{{"vertices":{},"collisions":{},"found":{},"findable":{},"covered":{},"fakes":{},"splits":{},"merged":{}}}}})",
                              score->vertices, score->collisions, score->found, score->findable, score->covered, score->fakes, score->splits, score->merged);
        }
        json << "}\n";
        json.flush();
      }
      tracker.endTimeframe();
      continue;
    }
    profile::Profile planProfile;
    std::optional<TimeframePlan<Det>> plan;
    double planSeconds = 0.;
    if (!device) {
      const profile::Activation activation{planProfile, pool.get()};
      const profile::Zone zone{profile::Phase::kPlan};
      plan.emplace(tracker.planTimeframe(frame, vertices, blocksPerSlice, param.sliceHaloPS));
    }
    planSeconds = planProfile[profile::Phase::kPlan].seconds + vertexSeconds;
    if (device) {
      const auto run = interface.track(frame, vertices);
      if (run.dropped) {
        throw std::runtime_error("the tracking ran out of device memory");
      }
      const auto& deviceTracks = *run.result.tracks;
      const auto& deviceRun = interface.deviceReport();
      planSeconds = deviceVertexing ? vertexed.planSeconds : run.planSeconds;
      const double deviceSeconds = vertexed.seconds + run.seconds - planSeconds;
      const double deviceVertexSeconds = deviceRun.vertexSeconds + run.planSeconds + deviceRun.planUploadSeconds;
      if (deviceVertexing) {
        std::cout << std::format("vertexer on {}: {} vertices from {} candidates of {} triplets in {:.4f} s, the plan with them uploaded\n", options.backend,
                                 vertexed.found->vertices.size(), vertexed.found->candidates, vertexed.found->triplets, deviceVertexSeconds);
        std::cout << std::format("digest, vertices: {:016x}\n", vertexDigest(vertexed.found->vertices, vertexed.found->contributors));
        if (last) {
          // The host's list on the same frame, against the device's.
          auto host = vertexing->find(frame, blocksPerSlice, param.sliceHaloPS);
          if (host.status != StageStatus::kComplete) {
            throw std::runtime_error(std::format("the host's vertexing pass failed with status {}", static_cast<int>(host.status)));
          }
          foundVertices = host.found.vertices;
          std::cout << describeVertexAgreement(*vertexed.found, host.found) << '\n';
        }
      }
      const auto deviceDigest = tracksDigest(deviceTracks);
      walls.push_back(planSeconds + deviceSeconds);
      digests.push_back(deviceDigest);
      // The tracks in cluster numbers, as the host's lines carry them: equal
      // across layouts, and to the host's in the deterministic build.
      const auto deviceClusterDigest = tracksClusterDigest<Det>(deviceTracks, frame);
      if (!options.tracksOut.empty() && last) {
        writeTracks<Det>(options.tracksOut, deviceTracks, frame);
      }
      std::cout << std::format("\n=== run {}{}: tracks, {}, {}: {} in {:.4f} s (plan {:.4f} s, upload {:.4f} s, iterations", repetition,
                               repetition == 0U ? " (cold)" : "", label, options.backend, deviceTracks.size(), planSeconds + deviceSeconds, planSeconds,
                               deviceRun.uploadSeconds);
      for (const double seconds : deviceRun.iterationSeconds) {
        std::cout << std::format(" {:.4f}", seconds);
      }
      std::cout << std::format(" s, download {:.4f} s); arena peak {:.1f} MB\n", deviceRun.downloadSeconds,
                               static_cast<double>(interface.deviceArenaPeakBytes()) / 1.e6);
      std::cout << std::format("digest, tracks: {:016x}\n", deviceDigest);
      std::cout << std::format("digest, tracks in clusters: {:016x}\n", deviceClusterDigest);
      for (std::size_t iteration = 0; iteration < deviceRun.iterations.size(); ++iteration) {
        const auto& report = deviceRun.iterations[iteration];
        std::cout << std::format(
          "  iteration {}: {} records, {} duplicates, {} mismatches, {} rounds, {} tracks, {} readbacks; arena peak {:.1f} MB\n", iteration,
          report.records, report.duplicates, report.mismatches, report.rounds, report.tracks, report.readbacks, static_cast<double>(report.arenaPeakBytes) / 1.e6);
      }
      const auto score = scoreAndPrint(trackRows(deviceTracks));
      if (json.is_open()) {
        json << std::format(R"({{"input":"{}",{}"run":{},"backend":"{}","device":"{}","planSeconds":{:.9f},"runSeconds":{:.9f},"uploadSeconds":{:.9f},)"
                            R"("downloadSeconds":{:.9f},"arenaPeakBytes":{},"maxRssBytes":{},"digest":"{:016x}","tracksClusterDigest":"{:016x}","tracks":{},)"
                            R"("iterations":[)",
                            options.input, gridJson, repetition, options.backend, interface.deviceProperties().name, planSeconds, deviceSeconds,
                            deviceRun.uploadSeconds, deviceRun.downloadSeconds, interface.deviceArenaPeakBytes(), maxRssBytes(), deviceDigest,
                            deviceClusterDigest, deviceTracks.size());
        for (std::size_t iteration = 0; iteration < deviceRun.iterations.size(); ++iteration) {
          const auto& report = deviceRun.iterations[iteration];
          json << (iteration == 0U ? "" : ",")
               << std::format(R"({{"seconds":{:.9f},"records":{},"rounds":{},"tracks":{},)"
                              R"("arenaPeakBytes":{},"roads":{},"slices":[)",
                              deviceRun.iterationSeconds[iteration], report.records, report.rounds, report.tracks, report.arenaPeakBytes,
                              std::accumulate(deviceRun.slices[iteration].begin(), deviceRun.slices[iteration].end(), std::uint64_t{0U},
                                              [](auto total, const auto& slice) { return total + slice.roads; }));
          for (std::size_t index = 0U; index < deviceRun.slices[iteration].size(); ++index) {
            const auto& slice = deviceRun.slices[iteration][index];
            json << (index == 0U ? "" : ",")
                 << std::format(R"({{"tracklets":{},"permutation":{},"triplets":{},"neighbours":{},"roads":{},"records":{},)"
                                R"("chained":{},"readbacks":{},)"
                                R"("candidatePairs":{},"candidateSources":{},"candidateMax":{},"candidateChunks":{},)"
                                R"("freeBytesAtTracklets":{},"tripletCandidates":{},"tripletCandidateMax":{},"tripletChunks":{},)"
                                R"("freeBytesAtTriplets":{}}})",
                                slice.tracklets, slice.permutation, slice.triplets, slice.neighbours, slice.roads, slice.records,
                                slice.chained, slice.readbacks, slice.trackletRun.candidates, slice.trackletSources,
                                slice.trackletRun.largest, slice.trackletRun.chunks, slice.trackletRun.freeBytes, slice.tripletRun.candidates,
                                slice.tripletRun.largest, slice.tripletRun.chunks, slice.tripletRun.freeBytes);
          }
          json << "]}";
        }
        json << ']';
        if (score) {
          json << R"(,"score":)" << jsonScore(score->first, comparison, score->second);
        }
        json << "}\n";
        json.flush();
      }
      if (options.hostReference && last) {
        const std::span<const VertexSeed> vertices = vertexing ? std::span<const VertexSeed>{foundVertices} : std::span<const VertexSeed>{input.vertices};
        const auto plan = tracker.planTimeframe(frame, vertices, blocksPerSlice, param.sliceHaloPS);
        const auto result = tracker.run(frame, plan);
        if (result.status != StageStatus::kComplete) {
          throw std::runtime_error(std::format("the host reference failed with status {}", static_cast<int>(result.status)));
        }
        Digest hostDigest;
        hostDigest.columns(result.tracks->view(), result.tracks->size());
        const auto compared = compareColumns<TrackOf<maxChainLayers<Det>>::template RecordT>(result.tracks->view(), result.tracks->size(),
                                                                                             std::as_const(deviceTracks).view(), deviceTracks.size());
        referenceSame = compared.sameElements && compared.recordsDiffering == 0U;
        std::cout << std::format("\nhost reference: {} tracks, digest, tracks: {:016x}; device {} tracks: ", result.tracks->size(), hostDigest.result(),
                                 deviceTracks.size());
        if (!compared.sameElements) {
          std::cout << "the track counts differ\n";
        } else if (referenceSame) {
          std::cout << "identical\n";
        } else {
          std::cout << std::format("{} tracks differ; by column (memlayout order):", compared.recordsDiffering);
          for (std::size_t column = 0; column < compared.differingByColumn.size(); ++column) {
            if (compared.differingByColumn[column] != 0U) {
              std::cout << std::format(" {}: {} (max |diff| {:.3g})", column, compared.differingByColumn[column], compared.maxAbsDifference[column]);
            }
          }
          std::cout << '\n';
        }
        const auto hostScore = scoreAndPrint(trackRows(*result.tracks));
        static_cast<void>(hostScore);
        tracker.endTimeframe();
      }
      continue;
    }
    work.clear();
    workBefore = workCounts();
    linkWork.clear();
    linkWorkBefore = linkWorkCounts();
    const auto runStart = std::chrono::steady_clock::now();
    const auto result = [&] {
      if (resumeFrom) {
        return tracker.resume(frame, *plan, *resumeFrom);
      }
      if (!options.saveState.empty()) {
        return tracker.run(frame, *plan, options.stateIterations);
      }
      return tracker.run(frame, *plan);
    }();
    const double runSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - runStart).count();
    if (result.status != StageStatus::kComplete) {
      throw std::runtime_error(std::format("run {} failed with status {}", repetition, static_cast<int>(result.status)));
    }
    if (!options.saveState.empty() && last) {
      const auto [frameFingerprint, conditions] = checkpointFingerprints();
      writeCheckpoint(options.saveState, tracker.checkpoint(result), settings, frameFingerprint, conditions);
      std::cout << std::format("checkpoint after {} iterations into {}\n", options.stateIterations, options.saveState);
    }
    const auto digest = trackingDigest(result, tracker.used());
    walls.push_back(planSeconds + runSeconds);
    digests.push_back(digest);

    std::cout << std::format("\n=== run {}{}: tracks, {}: {} in {:.3f} s (plan {:.3f} s, run {:.3f} s); pool peak {:.1f} MB\n", repetition,
                             repetition == 0U ? " (cold)" : "", label, result.tracks->size(), planSeconds + runSeconds, planSeconds, runSeconds,
                             static_cast<double>(result.profile[profile::Phase::kRun].peakBytes) / 1.e6);
    std::cout << std::format("digest, {}: {:016x}\n", label, digest);
    // The stages' largest outputs of any slice:
    // what the pool's peak is made of.
    std::cout << std::format("largest slice: {} tracklets, {} triplets, {} neighbours; staging {} records\n", tracker.trackletHighWater(),
                             tracker.tripletHighWater(), tracker.neighbourHighWater(), tracker.stagingHighWater());
    // The tracks alone: the digest also counts the selection records' bytes,
    // which a record's layout changes without changing a track.
    std::cout << std::format("digest, tracks: {:016x}\n", tracksDigest(*result.tracks));
    // In cluster numbers: unchanged by another grid layout, which renumbers
    // the rows.
    std::cout << std::format("digest, clusters: {:016x}\n", clusterDigest(result, frame, tracker.used()));
    // The tracks alone in cluster numbers, sorted by their clusters: unchanged
    // by the slicing too, which changes the summaries' counts and the tracks'
    // order.
    std::cout << std::format("digest, tracks in clusters: {:016x}\n", tracksClusterDigest<Det>(*result.tracks, frame));
    if (!options.tracksOut.empty() && last) {
      writeTracks<Det>(options.tracksOut, *result.tracks, frame);
    }
    if (truth && last) {
      writeTruthStudies<Det>(
        options, *result.tracks, *truth, comparison, frame, tracking, *material, propagator, studies.followed,
        [&](std::size_t iteration) -> const IterationParameters<Det>& { return tracker.iteration(iteration); },
        [&](std::uint64_t id, std::size_t iteration) { return funnel ? funnel->stages(id, iteration) : std::nullopt; },
        [&](std::size_t iteration, std::span<const std::uint32_t> rows) { return funnel ? funnel->recordChi2(iteration, rows) : std::nullopt; });
    }
    for (std::size_t iteration = 0; iteration < result.iterations.size(); ++iteration) {
      for (const auto& line : study::describeIteration(result.iterations[iteration], iteration)) {
        std::cout << line << '\n';
      }
      if (iteration < work.size()) {
        printWork(work[iteration], iteration);
      }
      if (iteration < linkWork.size()) {
        printLinkWork<Det>(linkWork[iteration], iteration);
      }
    }
    const auto hostRows = trackRows(*result.tracks);
    const auto score = scoreAndPrint(hostRows);
    if (truth && last && (!options.breakdown.empty() || !options.breakdownParticles.empty())) {
      const auto view = result.tracks->view();
      const study::TrackBreakdown breakdown(truth->rows, truth->particles, hostRows, comparison,
                                            std::span<const std::uint8_t>{view.iteration, result.tracks->size()});
      if (!options.breakdown.empty()) {
        breakdown.write(options.breakdown);
      }
      if (!options.breakdownParticles.empty()) {
        breakdown.writeParticles(options.breakdownParticles);
      }
    }
    if (studies.follow && options.followStops) {
      studies.follow->printStops();
    }
    if (completeAnatomy) {
      completeAnatomy->print();
      if (!options.fateOut.empty()) {
        completeAnatomy->writeFates(options.fateOut);
      }
      // The seeding funnel: a pure road from the stage
      // funnel where it ran.
      completeAnatomy->printSeeding([&](std::uint64_t id, std::size_t iteration) -> std::optional<std::pair<std::uint32_t, int>> {
        if (!funnel) {
          return std::nullopt;
        }
        const auto stages = funnel->stages(id, iteration);
        if (!stages) {
          return std::nullopt;
        }
        return std::pair{stages->longest[static_cast<std::size_t>(study::FunnelStage::kRoads)], static_cast<int>(stages->fitLoss)};
      });
      completeAnatomy->printRecordLosses();
      if (!options.anatomyBeamLine.empty()) {
        completeAnatomy->writeBeamLine(options.anatomyBeamLine);
      }
      if (!options.anatomyFollowLoss.empty()) {
        completeAnatomy->writeFollowLoss(options.anatomyFollowLoss);
      }
    }
    if (funnel) {
      for (const auto& line : funnel->report(trackRows(*result.tracks))) {
        std::cout << line << '\n';
      }
      std::vector<std::size_t> trackIterations(result.tracks->size());
      for (std::size_t track = 0; track < trackIterations.size(); ++track) {
        trackIterations[track] = result.tracks->view().iteration[track];
      }
      for (const auto& line : funnel->lengthReport(trackRows(*result.tracks), trackIterations, tracker)) {
        std::cout << line << '\n';
      }
    }
    if (vertexStudy) {
      for (const auto& line : vertexStudy->report()) {
        std::cout << line << '\n';
      }
      if (!options.vertexList.empty()) {
        vertexStudy->writeStageVertices(options.vertexList);
        std::cout << std::format("the vertexer's vertices written to {}\n", options.vertexList);
      }
      if (!options.vertexReplay.empty()) {
        vertexStudy->writeReplay(options.vertexReplay);
        std::cout << std::format("the replay's lists written to {}-{{vertices,collisions,dropped}}.csv and {}-owned.txt\n", options.vertexReplay,
                                 options.vertexReplay);
      }
      std::cout << std::format("peak RSS {:.2f} GB\n", static_cast<double>(maxRssBytes()) / 1.e9);
    }
    if (options.tree) {
      study::report(std::cout, result.profile, profile::Phase::kRun, threads);
    }
    if (options.iterationTrees) {
      for (std::size_t iteration = 0; iteration < result.iterations.size(); ++iteration) {
        std::cout << std::format("--- iteration {}\n", iteration);
        study::report(std::cout, result.iterations[iteration].profile, profile::Phase::kIteration, threads);
      }
    }
    if (json.is_open()) {
      json << std::format(
        R"({{"input":"{}",{}"run":{},"threads":{},"chunkSize":{},"planSeconds":{:.9f},"runSeconds":{:.9f},"maxRssBytes":{},"poolPeakBytes":{},)"
        R"("digest":"{:016x}","clusterDigest":"{:016x}","tracksClusterDigest":"{:016x}","tracks":{},)",
        options.input, gridJson, repetition, threads, chunkSize, planSeconds, runSeconds, maxRssBytes(), result.profile[profile::Phase::kRun].peakBytes, digest,
        clusterDigest(result, frame, tracker.used()), tracksClusterDigest<Det>(*result.tracks, frame), result.tracks->size());
      json << R"("iterations":[)";
      for (std::size_t iteration = 0; iteration < result.iterations.size(); ++iteration) {
        const auto& summary = result.iterations[iteration];
        json << (iteration == 0U ? "" : ",")
             << std::format(R"({{"roads":{},"records":{},"rounds":{},"tracks":{},"followed":{},"followedOutward":{},"followedInward":{},"dropped":{},)"
                            R"("followedEarlier":{},"followedEarlierOutward":{},"followedEarlierInward":{},)",
                            summary.roads, summary.records, summary.rounds, summary.tracks, summary.tracksFollowed, summary.clustersFollowedOutward,
                            summary.clustersFollowedInward, summary.tracksDropped, summary.tracksFollowedEarlier, summary.clustersFollowedEarlierOutward,
                            summary.clustersFollowedEarlierInward);
        // The counters build's work counts, in WorkCounter's order.
        if (iteration < work.size()) {
          json << R"("work":)" << jsonArray(work[iteration]) << ',';
        }
        json << R"("profile":)";
        study::writeJson(json, summary.profile);
        json << '}';
      }
      json << R"(],"profile":)";
      auto profiled = result.profile;
      profiled += planProfile;
      study::writeJson(json, profiled);
      if (score) {
        json << R"(,"score":)" << jsonScore(score->first, comparison, score->second);
      }
      json << "}\n";
      json.flush();
    }
    if (last && !options.trace.empty()) {
      std::ofstream out{options.trace};
      trace.write(out);
      std::cout << std::format("trace written to {}\n", options.trace);
    }
    // The timeframe's memory back to the pool: the next
    // run pays first touch again, as a real run's timeframes do.
    tracker.endTimeframe();
  }
  if (!options.fitInputs.empty()) {
    study::writeFitInputs(options.fitInputs, fitInputs);
    std::cout << std::format("fit inputs of {} triplets and {} neighbours written to {}\n", fitInputs.triplets.size(), fitInputs.neighbours.size(),
                             options.fitInputs);
  }

  if (seedingReal && !options.seedingParticles.empty()) {
    std::ofstream particles(options.seedingParticles);
    seedingReal->writeParticles(particles);
    particles.close();
    if (!particles) {
      throw std::runtime_error("cannot write " + options.seedingParticles);
    }
  }
  const bool same = std::ranges::all_of(digests, [&](std::uint64_t digest) { return digest == digests.front(); });
  std::vector<double> warm{walls.begin() + (walls.size() > 1U ? 1 : 0), walls.end()};
  std::ranges::sort(warm);
  std::cout << std::format("\n{} runs {}: cold {:.4f} s; warm min {:.4f} s, median {:.4f} s, max {:.4f} s; digests {}\n", walls.size(),
                           device ? std::format("on {}", options.backend) : std::format("at {} threads", threads), walls.front(), warm.front(),
                           warm[warm.size() / 2U], warm.back(), same ? "all equal" : "DIFFER");
  return same && referenceSame ? 0 : 2;
}

extern template int track<detectors::ITS>(const Options& options);
extern template int track<detectors::MFT>(const Options& options);
extern template int track<detectors::TRKFT3>(const Options& options);

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_TRACKBENCH_H_
