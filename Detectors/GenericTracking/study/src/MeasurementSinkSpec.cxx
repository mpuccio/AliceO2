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

#include "GenericTracking/Study/MCVertexReaderSpec.h"
#include "GenericTracking/Study/Digest.h"
#include "GenericTracking/Study/MeasurementSinkSpec.h"

#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Study/Alice/Hits.h"
#include "GenericTracking/Study/Alice/Tracks.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/LayerMaterial.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/SlicePlan.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/Tracker.h"
#include "GenericTracking/Study/Validation.h"
#include "GenericTracking/Tracking/VertexingPass.h"
#include "GenericTracking/Study/VertexQA.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/BeamLine.h"
#include "GenericTracking/Workflow/Conditions.h"
#include "GenericTracking/Workflow/ConditionsTRK.h"
#include "GenericTracking/Workflow/TrackingInterface.h"
#include "GenericTracking/Workflow/TrackerSpec.h"
#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Study/TrackPulls.h"
#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Workflow/RecoParam.h"
#include "GenericTracking/Study/StageFunnel.h"
#include "GenericTracking/Study/ProfileReport.h"
#include "GenericTracking/Study/Statistics.h"
#include "GenericTracking/Study/TrackingInput.h"
#include "GenericTracking/Study/SpeciesOracle.h"
#include "GenericTracking/Study/TrackingQA.h"

#include "Framework/DataProcessorSpec.h"
#include "Framework/EndOfStreamContext.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/Logger.h"
#include "Framework/Task.h"
#include "Framework/TimingInfo.h"
#include "SimulationDataFormat/ConstMCTruthContainer.h"
#include "SimulationDataFormat/DigitizationContext.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/O2DatabasePDG.h"
#include "Steer/MCKinematicsReader.h"

#include "CommonUtils/StringUtils.h"
#include "CommonUtils/TreeStreamRedirector.h"

#include <fmt/format.h>
#include <gsl/span>

#include <algorithm>
#include <array>
#include <initializer_list>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <optional>
#include <memory>
#include <ranges>
#include <span>
#include <stdexcept>
#include <tuple>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2::framework;

namespace o2::gtrk::study
{

namespace
{

// The sink is compiled once per detector, which GTRK_SINK_DETECTOR names
// (CMakeLists.txt): everything below is bound to it through this one alias.
#ifndef GTRK_SINK_DETECTOR
#error "GTRK_SINK_DETECTOR names the sink's detector"
#endif
using Det = detectors::GTRK_SINK_DETECTOR;
using Frame = TimeFrame<Det>;
using Tracklets = StageOutput<TrackletT>;
using Controller = workflow::TrackingInterface<Det>::Controller;
/// The stage funnel of a tracker's own run (StageFunnel.h), read by the text
/// reports per iteration.
using Funnel = StageFunnel<Det, Controller>;

/// One tracklet's bytes across its columns.
constexpr std::size_t kTrackletBytes = Tracklet{}.apply([](const auto&... field) { return (sizeof(field) + ... + 0U); });
constexpr std::size_t kCountBins = 20U;
constexpr std::array<Scalar, 13> kPtEdges{0.F, 0.1F, 0.2F, 0.3F, 0.5F, 1.F, 1.25F, 1.5F, 2.F, 3.F, 5.F, 10.F, std::numeric_limits<Scalar>::infinity()};
constexpr std::array<Scalar, 17> kCotThetaEdges{
  -std::numeric_limits<Scalar>::infinity(), -3.F, -2.F, -1.5F, -1.F, -0.75F, -0.5F, -0.25F, 0.F, 0.25F, 0.5F, 0.75F, 1.F, 1.5F, 2.F, 3.F,
  std::numeric_limits<Scalar>::infinity()};
constexpr std::size_t kPtBins = kPtEdges.size() - 1U;
constexpr std::size_t kCotThetaBins = kCotThetaEdges.size() - 1U;
using TruthBreakdown = std::array<std::array<std::array<std::uint64_t, kCotThetaBins>, kPtBins>, nLinks<Det>>;

struct TrackletRunSummary {
  std::uint64_t sources = 0U;
  std::uint64_t tracklets = 0U;
  std::uint64_t haloSourceTracklets = 0U;
  std::uint64_t haloTriplets = 0U;
  std::array<std::uint64_t, kCountBins> countHistogram{};
  std::uint32_t maxPerSource = 0U;
  std::size_t trackletHighWater = 0U;
  std::size_t stagingHighWater = 0U;
  std::size_t stagingPredictionMisses = 0U;
  std::size_t poolHighWater = 0U;
  std::array<std::uint64_t, nLinks<Det>> truePairs{};
  std::array<std::uint64_t, nLinks<Det>> truePairsKept{};
  TruthBreakdown truthBreakdown{};
  TruthBreakdown truthBreakdownKept{};
  // Triplets: a true triplet is two MC-true pairs (a, b) and (b, c) of one
  // particle, a row being every particle's it carries among its labels.
  std::uint64_t firstTracklets = 0U;
  std::uint64_t triplets = 0U;
  std::uint64_t trueTriplets = 0U;
  std::uint64_t trueTripletsPairsKept = 0U; ///< both pairs kept as tracklets
  std::uint64_t trueTripletsKept = 0U;
  std::uint64_t trueTripletsKinked = 0U; ///< of them, kept only with the kink allowance
  std::uint64_t tripletsKinked = 0U;     ///< all core-owned triplets kept only with it
  std::array<std::uint64_t, kPtBins> trueTripletsPairsKeptByPt{};
  std::array<std::uint64_t, kPtBins> trueTripletsKeptByPt{};
  double curvaturePullSum = 0.;
  double curvaturePullSum2 = 0.;
  double curvatureResidualSum2 = 0.; ///< relative to the MC curvature
  // The same pull by the particle's pT: the triplet's scattering is taken at
  // pT_min, so above it the pull's rms falls below 1 by design.
  std::array<double, kPtBins> curvaturePullSumByPt{};
  std::array<double, kPtBins> curvaturePullSum2ByPt{};
  std::size_t tripletHighWater = 0U;
  // What the triplets that are not MC-true are.
  std::uint64_t notTrueCombinatorial = 0U; ///< a, b, c not one particle
  std::uint64_t notTrueSecondary = 0U;
  std::uint64_t notTrueSoftPrimary = 0U;  ///< a primary below the iteration's pT_min
  std::uint64_t notTrueOtherPrimary = 0U; ///< a primary above it, off its first pass or not consecutive
  // Neighbours, each owned by the slice whose core holds its first
  // measurement. A true neighbour is two true triplets (a, b, c), (b, c, d).
  std::uint64_t neighbours = 0U;
  std::uint64_t trueNeighbours = 0U;
  std::uint64_t trueNeighboursTripletsKept = 0U; ///< both true triplets kept
  std::uint64_t trueNeighboursKept = 0U;
  std::array<std::uint64_t, kPtBins> trueNeighboursTripletsKeptByPt{};
  std::array<std::uint64_t, kPtBins> trueNeighboursKeptByPt{};
  std::array<std::uint64_t, kCotThetaBins> trueNeighboursTripletsKeptByCot{};
  std::array<std::uint64_t, kCotThetaBins> trueNeighboursKeptByCot{};
  std::uint64_t neighboursNotTrueCombinatorial = 0U;
  std::uint64_t neighboursNotTrueSecondary = 0U;
  std::uint64_t neighboursNotTrueSoftPrimary = 0U;
  std::uint64_t neighboursNotTrueOtherPrimary = 0U;
  std::array<std::uint64_t, kCountBins> neighbourHistogram{}; ///< predecessors per core-owned triplet
  std::uint32_t maxNeighbours = 0U;
  // Per MC track, the highest level of its kept true triplets against the
  // longest chain of them the true neighbours allow, by that length.
  std::array<std::uint64_t, 12> tracksByChain{};
  std::array<std::uint64_t, 12> tracksReachingChain{};
  std::array<std::uint64_t, 12> tracksBelowChain{};
  std::array<std::uint64_t, 12> tracksAboveChain{};
  std::uint64_t tripletsWorked = 0U; ///< core and halo, summed over the slices
  std::uint64_t chained = 0U;        ///< of those, swept: with a predecessor
  std::uint64_t slices = 0U;
  std::size_t neighbourHighWater = 0U;
};

/// A true triplet as the owning slice found it.
struct TrueTripletState {
  bool kept = false;
  std::uint32_t level = 0U;
};

struct TruePair {
  std::uint64_t key = 0U;
  std::uint64_t particle = 0U; ///< whose pair it is (MCCompLabel::getTrackEventSourceID)
  std::uint16_t link = 0U;
  Scalar pt = 0.F;       ///< of the MC track [GeV]
  Scalar cotTheta = 0.F; ///< pz / pT of the MC track at its origin
  Scalar originR = 0.F;  ///< transverse distance of the MC track's origin from the beam [cm]
};

/// The longest MC-true chain along a pair: one particle's pairs joined end
/// to start, its length in clusters up to the pair's target, and the pair
/// before it on it.
struct TrueChain {
  static constexpr std::size_t kNoPair = std::numeric_limits<std::size_t>::max();
  std::uint32_t length = 2U;
  std::size_t previous = kNoPair;
};

/// Each of \a pairs' (in key order) longest chain.
[[nodiscard]] inline std::vector<TrueChain> trueChains(std::span<const TruePair> pairs)
{
  std::vector<TrueChain> chains(pairs.size());
  for (bool changed = true; changed;) {
    changed = false;
    for (std::size_t pair = 0; pair < pairs.size(); ++pair) {
      const auto target = pairs[pair].key & 0xffffffffU;
      const auto next = std::ranges::lower_bound(pairs, target << 32U, {}, &TruePair::key);
      for (auto at = next; at != pairs.end() && (at->key >> 32U) == target; ++at) {
        const auto successor = static_cast<std::size_t>(at - pairs.begin());
        if (at->particle == pairs[pair].particle && chains[pair].length + 1U > chains[successor].length) {
          chains[successor] = {.length = chains[pair].length + 1U, .previous = pair};
          changed = true;
        }
      }
    }
  }
  return chains;
}

template <std::size_t N>
[[nodiscard]] std::size_t binOf(Scalar value, const std::array<Scalar, N>& edges)
{
  const auto upper = std::ranges::upper_bound(edges, value);
  if (upper == edges.begin()) {
    return 0U;
  }
  return std::min<std::size_t>(static_cast<std::size_t>(upper - edges.begin() - 1), N - 2U);
}

/// One MC-true kept triplet's curvature pull and relative residual, overall
/// and in its particle's pT bin.
void addCurvaturePull(TrackletRunSummary& summary, Scalar pt, double pull, double relativeResidual)
{
  summary.curvaturePullSum += pull;
  summary.curvaturePullSum2 += pull * pull;
  summary.curvaturePullSumByPt[binOf(pt, kPtEdges)] += pull;
  summary.curvaturePullSum2ByPt[binOf(pt, kPtEdges)] += pull * pull;
  summary.curvatureResidualSum2 += relativeResidual * relativeResidual;
}

/// The MC-true kept triplets' curvature pulls, overall and by pT.
void logTripletCurvature(const TrackletRunSummary& summary)
{
  const auto count = static_cast<double>(summary.trueTripletsKept);
  if (count > 0.) {
    const double mean = summary.curvaturePullSum / count;
    LOGP(info, "  curvature against MC 1/R(p): pull mean {:.3f}, rms {:.3f}; relative residual rms {:.4f}", mean,
         std::sqrt(std::max(summary.curvaturePullSum2 / count - (mean * mean), 0.)), std::sqrt(summary.curvatureResidualSum2 / count));
    std::string pullByPt;
    for (std::size_t bin = 1; bin < kPtBins; ++bin) {
      const auto inBin = static_cast<double>(summary.trueTripletsKeptByPt[bin]);
      if (inBin > 0.) {
        const double binMean = summary.curvaturePullSumByPt[bin] / inBin;
        pullByPt += fmt::format(" {}-{}: {:.2f} ({:.2f})", kPtEdges[bin], kPtEdges[bin + 1U], binMean,
                                std::sqrt(std::max(summary.curvaturePullSum2ByPt[bin] / inBin - (binMean * binMean), 0.)));
      }
    }
    LOGP(info, "    by pT [GeV], pull mean (rms), the scattering at pT_min:{}", pullByPt);
  }
}

[[nodiscard]] std::size_t countBin(std::uint32_t count)
{
  if (count < 16U) {
    return count;
  }
  if (count < 32U) {
    return 16U;
  }
  if (count < 64U) {
    return 17U;
  }
  if (count < 128U) {
    return 18U;
  }
  return 19U;
}

/// One true crossing of a layer by a particle, from the MC hits (Hits.h): the
/// chip segment's midpoint, the momentum entering it and the particle's mass.
struct TrueCrossing {
  std::uint32_t layer = 0U;
  Frame::Index row = 0U; ///< the cluster the crossing made
  std::array<double, 3> position{};
  std::array<double, 3> momentum{};
  double mass = 0.;
  int pdg = 0;
};

/// A particle's crossings of its rows' chips, in layer order; empty without
/// the hits, or where one row has no crossing.
using TrueCrossings = std::function<std::vector<TrueCrossing>(std::uint64_t particle, std::span<const Frame::Index> rows)>;

/// What the MC hits say of one crossing (for the tracks still found
/// shorter): how far the particle's cluster sits from where it
/// crossed, how hard it had scattered on the way from the previous crossing,
/// and its pT there.
struct TrueCrossingReport {
  std::uint32_t layer = 0U;
  /// The space angle between the momentum entering this layer and the helix
  /// continued from the previous crossing, over the previous layer's Highland
  /// theta0 at the true momentum: Rayleigh with unit scale for the model's
  /// Gaussian core, so 5 is a hard scatter (P = 4e-6). NaN on the first
  /// crossing.
  double kink = std::numeric_limits<double>::quiet_NaN();
  /// The same angle over theta0 from the LUT's x/X0 along the straight path
  /// between the two crossings: what the material actually crossed allows,
  /// against the model's per-layer attribution. NaN without a LUT.
  double kinkPath = std::numeric_limits<double>::quiet_NaN();
  double pathX2X0 = 0.;         ///< the LUT's x/X0 along that path
  double pathNormalX2X0 = 0.;   ///< the same at normal incidence: over the path factor at the previous crossing
  std::uint32_t scatterer = 0U; ///< the layer the kink is attributed to: the previous crossing's
  double offsetCM = 0.;         ///< the cluster's distance from the crossing, the larger of r-phi and z
  double pt = 0.;               ///< entering the layer [GeV]
  double z = 0.;                ///< of the crossing [cm]
  int pdg = 0;
};

/// \a crossings' reports, one per layer (a second crossing of the same layer,
/// in a stave overlap, is skipped).
std::vector<TrueCrossingReport> trueAnatomy(const Frame& frame, const ScatteringTable<Det>& scattering, double bzKiloGauss,
                                            std::span<const TrueCrossing> crossings, const o2::base::MatLayerCylSet* lut = nullptr)
{
  std::vector<TrueCrossingReport> result;
  const TrueCrossing* previous = nullptr;
  for (const auto& crossing : crossings) {
    if (previous != nullptr && previous->layer == crossing.layer) {
      continue;
    }
    const auto m = frame.measurement(crossing.row);
    const auto& position = crossing.position;
    const auto& p = crossing.momentum;
    TrueCrossingReport report{.layer = crossing.layer, .pt = std::hypot(p[0], p[1]), .z = position[2], .pdg = crossing.pdg};
    {
      const double r = std::hypot(position[0], position[1]);
      const double dx = static_cast<double>(m.x) - position[0];
      const double dy = static_cast<double>(m.y) - position[1];
      const double across = r > 0. ? ((position[0] * dy) - (position[1] * dx)) / r : std::hypot(dx, dy);
      report.offsetCM = std::max(std::abs(across), std::abs(static_cast<double>(m.z) - position[2]));
    }
    if (previous != nullptr) {
      const auto& pa = previous->momentum;
      const double ptA = std::hypot(pa[0], pa[1]);
      const double pA = std::hypot(ptA, pa[2]);
      const double pB = std::hypot(report.pt, p[2]);
      if (ptA > 0. && pA > 0. && pB > 0.) {
        // The helix from the previous crossing turns by the angle its chord
        // subtends, towards the side the chord lies on; the dip stays.
        const double cx = position[0] - previous->position[0];
        const double cy = position[1] - previous->position[1];
        const double chord = std::hypot(cx, cy);
        const double radius = curvatureRadius(ptA, bzKiloGauss);
        const double alpha = radius > 0. && std::isfinite(radius) ? 2. * std::asin(std::min(chord / (2. * radius), 1.)) : 0.;
        const double sign = (pa[0] * cy) - (pa[1] * cx) >= 0. ? 1. : -1.;
        const double phi = std::atan2(pa[1], pa[0]) + (sign * alpha);
        const double lambda = std::atan2(pa[2], ptA);
        const std::array<double, 3> expected{std::cos(lambda) * std::cos(phi), std::cos(lambda) * std::sin(phi), std::sin(lambda)};
        const std::array<double, 3> actual{p[0] / pB, p[1] / pB, p[2] / pB};
        const double dot = (expected[0] * actual[0]) + (expected[1] * actual[1]) + (expected[2] * actual[2]);
        const std::array<double, 3> cross{(expected[1] * actual[2]) - (expected[2] * actual[1]), (expected[2] * actual[0]) - (expected[0] * actual[2]),
                                          (expected[0] * actual[1]) - (expected[1] * actual[0])};
        const double angle = std::atan2(std::hypot(cross[0], cross[1], cross[2]), dot);
        // Highland's theta0 for the previous layer at the true momentum and
        // incidence: k x (path factor) / (beta p)^2.
        const auto id = static_cast<LayerID>(previous->layer);
        const double rPrevious = std::hypot(previous->position[0], previous->position[1]);
        double normal = std::abs(pa[2]);
        if (Detector<Det>::isBarrelLayer(id)) {
          normal = rPrevious > 0. ? std::abs((pa[0] * previous->position[0]) + (pa[1] * previous->position[1])) / rPrevious : ptA;
        }
        const double pathFactor = normal > 0. ? pA / normal : 1.;
        const double beta = pA / std::hypot(pA, previous->mass);
        // NOLINTNEXTLINE(clang-analyzer-security.ArrayBound): a hit's layer is one of Det's
        const double theta0 = std::sqrt(static_cast<double>(scattering.layers[previous->layer].k) * pathFactor) / (beta * pA);
        report.kink = theta0 > 0. ? angle / theta0 : 0.;
        report.scatterer = previous->layer;
        if (lut != nullptr) {
          // The LUT along the straight path between the two crossings; its
          // x/X0 is the path's own, so no path factor.
          const auto budget =
            lut->getMatBudget(static_cast<float>(previous->position[0]), static_cast<float>(previous->position[1]), static_cast<float>(previous->position[2]),
                              static_cast<float>(position[0]), static_cast<float>(position[1]), static_cast<float>(position[2]));
          report.pathX2X0 = static_cast<double>(budget.meanX2X0);
          report.pathNormalX2X0 = report.pathX2X0 / pathFactor;
          const double theta0Path = std::sqrt(detail::highlandFactor(report.pathX2X0)) / (beta * pA);
          report.kinkPath = theta0Path > 0. ? angle / theta0Path : 0.;
        }
      }
    }
    result.push_back(report);
    previous = &crossing;
  }
  return result;
}

/// The bins the kinks are reported in: a Rayleigh of unit scale puts 39.3 %
/// below 1, 86.5 % below 2, 98.9 % below 3; 5 and beyond is a hard scatter.
constexpr std::array<double, 6> kKinkEdges{0., 1., 2., 3., 5., 10.};
constexpr double kHardKink = 5.;
constexpr double kTailKink = 3.;
[[nodiscard]] std::size_t kinkBin(double kink)
{
  std::size_t bin = 0U;
  for (std::size_t i = 1; i < kKinkEdges.size(); ++i) {
    if (kink >= kKinkEdges[i]) {
      bin = i;
    }
  }
  return bin;
}
[[nodiscard]] std::string kinkBinsLine(std::span<const std::uint64_t> counts)
{
  std::string line;
  for (std::size_t bin = 0; bin < kKinkEdges.size(); ++bin) {
    const auto upper = bin + 1U < kKinkEdges.size() ? fmt::format("{:.0f}", kKinkEdges[bin + 1U]) : std::string{"inf"};
    line += fmt::format(" [{:.0f},{}) {}", kKinkEdges[bin], upper, counts[bin]);
  }
  return line;
}

/// A cluster this far from its crossing [cm] is off it: six pixel sigmas.
constexpr double kOffsetCM = 0.003;

/// Every primary's MC-true pairs: consecutive layers of its first outward
/// pass, a row being every particle's it carries among its labels (\a rows,
/// the comparison's truth). A pair carried by two particles is the first's.
/// \a mcOf: a particle's MC track, by id.
template <typename McOf>
[[nodiscard]] std::vector<TruePair> mcTruePairs(const Frame& frame, std::span<const TruthRow> rows, const McOf& mcOf, const Geometry<Det>& geometry,
                                                const RecoSettings& param, double bzKiloGauss)
{
  struct SignalRow {
    std::uint64_t track = 0U;
    Frame::Index row = 0U;
  };

  std::vector<SignalRow> signalRows;
  signalRows.reserve(frame.size() / 4U);
  for (Frame::Index row = 0U; row < frame.size(); ++row) {
    for (std::uint8_t i = 0; i < rows[row].count; ++i) {
      signalRows.push_back({.track = rows[row].particles[i], .row = row});
    }
  }
  std::ranges::stable_sort(signalRows, {}, &SignalRow::track);

  // Down to the lowest iteration's pT_min: every particle some iteration is for.
  double minPtGeV = workflow::iterationSettings(param, 0).windows.minPtGeV;
  for (int iteration = 1; iteration < workflow::iterationCount(param); ++iteration) {
    minPtGeV = std::min(minPtGeV, workflow::iterationSettings(param, iteration).windows.minPtGeV);
  }
  std::vector<Frame::Index> firstPass;
  std::vector<TruePair> result;
  for (std::size_t begin = 0; begin < signalRows.size();) {
    std::size_t end = begin + 1U;
    while (end < signalRows.size() && signalRows[end].track == signalRows[begin].track) {
      ++end;
    }
    const auto* mc = mcOf(signalRows[begin].track);
    if (mc == nullptr || !mc->isPrimary() || mc->GetPt() < minPtGeV) {
      begin = end;
      continue;
    }

    firstPass.clear();
    for (std::size_t at = begin; at < end; ++at) {
      const auto row = signalRows[at].row;
      const auto measurement = frame.measurement(row);
      const auto layer = geometry.sensor(measurement.surface).layer;
      if (onFirstOutwardPass<Det>(*mc, layer, measurement.phi, secondCoordinate<Det>(layer, measurement.r, measurement.z), bzKiloGauss)) {
        firstPass.push_back(row);
      }
    }
    const auto distance = [&](Frame::Index row) {
      const auto measurement = frame.measurement(row);
      return std::hypot(static_cast<double>(measurement.x) - mc->Vx(), static_cast<double>(measurement.y) - mc->Vy(),
                        static_cast<double>(measurement.z) - mc->Vz());
    };
    std::ranges::sort(firstPass, {}, distance);

    for (std::size_t from = 0; from < firstPass.size();) {
      const auto layerOf = [&](std::size_t at) { return geometry.sensor(frame.measurement(firstPass[at]).surface).layer; };
      const auto fromLayer = layerOf(from);
      std::size_t fromEnd = from + 1U;
      while (fromEnd < firstPass.size() && layerOf(fromEnd) == fromLayer) {
        ++fromEnd;
      }
      if (fromEnd == firstPass.size()) {
        break;
      }
      const auto toLayer = layerOf(fromEnd);
      std::size_t toEnd = fromEnd + 1U;
      while (toEnd < firstPass.size() && layerOf(toEnd) == toLayer) {
        ++toEnd;
      }
      const auto link = linkOf<Det>(fromLayer, toLayer);
      if (link != kNoLink) {
        for (std::size_t source = from; source < fromEnd; ++source) {
          for (std::size_t target = fromEnd; target < toEnd; ++target) {
            result.push_back({.key = pairKey(firstPass[source], firstPass[target]),
                              .particle = signalRows[begin].track,
                              .link = static_cast<std::uint16_t>(link),
                              .pt = static_cast<Scalar>(mc->GetPt()),
                              .cotTheta = static_cast<Scalar>(mc->GetStartVertexMomentumZ() / mc->GetPt()),
                              .originR = static_cast<Scalar>(std::hypot(mc->Vx(), mc->Vy()))});
          }
        }
      }
      from = fromEnd;
    }
    begin = end;
  }

  std::ranges::stable_sort(result, {}, &TruePair::key);
  const auto duplicates = std::ranges::unique(result, {}, &TruePair::key);
  result.erase(duplicates.begin(), duplicates.end());
  return result;
}

std::string histogramText(const std::array<std::uint64_t, kCountBins>& histogram)
{
  std::string result;
  for (std::size_t bin = 0; bin < histogram.size(); ++bin) {
    if (!result.empty()) {
      result += ' ';
    }
    if (bin < 16U) {
      result += std::to_string(bin);
    } else if (bin == 16U) {
      result += "16-31";
    } else if (bin == 17U) {
      result += "32-63";
    } else if (bin == 18U) {
      result += "64-127";
    } else {
      result += "128+";
    }
    result += ':';
    result += std::to_string(histogram[bin]);
  }
  return result;
}

/// The species the reports split by: e, mu, pi, K, p and the others.
inline constexpr std::size_t kSpeciesCount = 6U;

/// A PDG code's species index (below kSpeciesCount).
[[nodiscard]] std::size_t speciesOf(int pdg)
{
  switch (std::abs(pdg)) {
    case 11:
      return 0U;
    case 13:
      return 1U;
    case 211:
      return 2U;
    case 321:
      return 3U;
    case 2212:
      return 4U;
    default:
      return 5U;
  }
}

/// The true pairs the tracklet stage lost, against the beam line: the first
/// requirement of a tracklet each fails, by kind of link and pT; for the
/// curvature and longitudinal cuts, how far beyond them (the pull), to see a
/// tail's shape. Pairs on a link the iteration does not take (one that skips a
/// layer, in an iteration without holes) cannot be found and are counted
/// apart.
/// The cut is the stage funnel's (pairCutLoss, StageFunnel.h).
void reportLostPairs(const Frame& frame, const IterationParameters<Det>& parameters, const VertexSeed& beamLine, std::span<const TruePair> truePairs,
                     const std::vector<bool>& pairKept, std::span<const TruthRow> rows, const TrueCrossings& crossings = {},
                     const ScatteringTable<Det>* scattering = nullptr, double bzKiloGauss = 0., double minPtGeV = 0.,
                     const o2::base::MatLayerCylSet* material = nullptr)
{
  if (!truePairs.empty()) {
    // The pairs the two spatial cuts lost, for the MC hits' verdict below.
    struct LostPair {
      std::size_t pair;
      std::size_t pullBin;
      PairLoss reason;
    };
    std::vector<LostPair> lostSpatial;
    constexpr auto kReasons = kPairLosses;
    std::array<std::array<std::uint64_t, kReasons>, 3> byKind{}; // barrel-barrel, with a disc, disc-disc
    std::array<std::array<std::uint64_t, kReasons>, kPtBins> byPt{};
    std::array<std::uint64_t, 3> lostByKind{};
    std::array<std::uint64_t, 3> pairsByKind{};
    std::array<std::uint64_t, kPtBins> pairsByPt{};
    std::array<std::uint64_t, kPtBins> keptByPt{};
    std::uint64_t notTaken = 0U;
    // Where the particle came from, transversely from the beam: a primary of
    // a displaced decay is no track from the beam line.
    constexpr std::array<Scalar, 6> kOriginEdges{0.F, 0.005F, 0.01F, 0.03F, 0.1F, std::numeric_limits<Scalar>::infinity()};
    std::array<std::uint64_t, kOriginEdges.size() - 1U> originAll{};
    std::array<std::uint64_t, kOriginEdges.size() - 1U> originCurvature{};
    constexpr std::array<double, 6> kPullEdges{3., 4., 5., 7., 10., 20.};
    std::array<std::uint64_t, kPullEdges.size()> longitudinalPulls{};
    std::array<std::uint64_t, kPullEdges.size()> curvaturePulls{};
    // The pairs lost to time: the target's box start less the source's, in
    // source boxes, by the link's inner layer.
    constexpr std::size_t kTickSlots = 7U; // <= -3, -2, -1, 0, 1, 2, >= 3
    std::array<std::array<std::uint64_t, kTickSlots>, Frame::kLayerCount> timeTicksByLayer{};
    for (std::size_t pair = 0; pair < truePairs.size(); ++pair) {
      const auto& truth = truePairs[pair];
      const auto from = Det::layerLinks[truth.link].from;
      const auto to = Det::layerLinks[truth.link].to;
      std::size_t kind = 1U;
      if (Detector<Det>::isBarrelLayer(from) && Detector<Det>::isBarrelLayer(to)) {
        kind = 0U;
      } else if (!Detector<Det>::isBarrelLayer(from) && !Detector<Det>::isBarrelLayer(to)) {
        kind = 2U;
      }
      if (!linkTaken(parameters, truth.link)) {
        ++notTaken;
        continue;
      }
      const auto ptBin = binOf(truth.pt, kPtEdges);
      ++pairsByKind[kind];
      ++pairsByPt[ptBin];
      ++originAll[binOf(truth.originR, kOriginEdges)];
      if (pairKept[pair]) {
        ++keptByPt[ptBin];
        continue;
      }
      ++lostByKind[kind];
      const auto sourceIndex = static_cast<Frame::Index>(truth.key >> 32U);
      const auto targetIndex = static_cast<Frame::Index>(truth.key & 0xffffffffU);
      const auto verdict = pairCutLoss(frame, parameters, beamLine, truth.link, sourceIndex, targetIndex);
      const auto reason = verdict.loss;
      const auto& terms = verdict.terms;
      const auto pullBin = [&](double pull, auto& histogram) {
        std::size_t at = 0U;
        for (std::size_t bin = kPullEdges.size(); bin-- > 0U;) {
          if (pull >= kPullEdges[bin]) {
            at = bin;
            break;
          }
        }
        ++histogram[at];
        return at;
      };
      if (reason == PairLoss::kTime) {
        // How many of the source's boxes apart the two start.
        const auto sourceRow = frame.measurement(sourceIndex);
        const auto targetRow = frame.measurement(targetIndex);
        const auto ticks =
          std::lround((static_cast<double>(targetRow.timePS) - static_cast<double>(sourceRow.timePS)) / static_cast<double>(sourceRow.widthPS));
        ++timeTicksByLayer[from][static_cast<std::size_t>(std::clamp<long>(ticks, -3L, 3L) + 3L)];
      } else if (reason == PairLoss::kCurvature) {
        ++originCurvature[binOf(truth.originR, kOriginEdges)];
        const auto pullAt = pullBin(
          (static_cast<double>(terms.curvature) - static_cast<double>(terms.bound)) / std::sqrt(static_cast<double>(terms.varCurvature)), curvaturePulls);
        lostSpatial.push_back({.pair = pair, .pullBin = pullAt, .reason = reason});
      } else if (reason == PairLoss::kLongitudinal) {
        const auto pullAt = pullBin(std::sqrt(static_cast<double>(terms.collinearity2) / static_cast<double>(terms.varCollinearity)), longitudinalPulls);
        lostSpatial.push_back({.pair = pair, .pullBin = pullAt, .reason = reason});
      }
      ++byKind[kind][static_cast<std::size_t>(reason)];
      ++byPt[ptBin][static_cast<std::size_t>(reason)];
    }
    {
      std::string ticksLine;
      for (std::size_t layer = 0; layer < timeTicksByLayer.size(); ++layer) {
        const auto& slots = timeTicksByLayer[layer];
        if (std::accumulate(slots.begin(), slots.end(), std::uint64_t{0}) == 0U) {
          continue;
        }
        ticksLine +=
          fmt::format(" {}->: <=-3 {} -2 {} -1 {} 0 {} +1 {} +2 {} >=+3 {};", layer, slots[0], slots[1], slots[2], slots[3], slots[4], slots[5], slots[6]);
      }
      LOGP(info, "    true pairs lost to time, target less source in source boxes, by the link's inner layer:{}", ticksLine);
    }
    constexpr std::array<const char*, 3> kKindNames{"barrel-barrel", "with a disc", "disc-disc"};
    const auto taken = std::accumulate(pairsByKind.begin(), pairsByKind.end(), std::uint64_t{0});
    const auto lost = std::accumulate(lostByKind.begin(), lostByKind.end(), std::uint64_t{0});
    LOGP(info, "  true pairs on the links the iteration takes: {}/{} found ({:.3f}%); {} more on links it does not take", taken - lost, taken,
         percent(taken - lost, taken), notTaken);
    for (std::size_t kind = 0; kind < kKindNames.size(); ++kind) {
      std::string reasons;
      for (const auto reason : kPairCutLosses) {
        reasons += fmt::format(" {} {}", kPairLossNames[static_cast<std::size_t>(reason)], byKind[kind][static_cast<std::size_t>(reason)]);
      }
      LOGP(info, "  true pairs lost, {:<13}: {}/{} ({:.3f}%):{}", kKindNames[kind], lostByKind[kind], pairsByKind[kind],
           percent(lostByKind[kind], pairsByKind[kind]), reasons);
    }
    for (std::size_t bin = 1; bin < kPtBins; ++bin) {
      if (pairsByPt[bin] == 0U) {
        continue;
      }
      std::string reasons;
      for (const auto reason : kPairCutLosses) {
        reasons += fmt::format(" {} {}", kPairLossNames[static_cast<std::size_t>(reason)], byPt[bin][static_cast<std::size_t>(reason)]);
      }
      LOGP(info, "    pT [{:.2f},{:.2f}): {}/{} found ({:.2f}%), {} lost:{}", kPtEdges[bin], kPtEdges[bin + 1U], keptByPt[bin], pairsByPt[bin],
           percent(keptByPt[bin], pairsByPt[bin]), pairsByPt[bin] - keptByPt[bin], reasons);
    }
    std::string longitudinal;
    std::string curvature;
    for (std::size_t bin = 0; bin < kPullEdges.size(); ++bin) {
      const auto upper = bin + 1U < kPullEdges.size() ? fmt::format("{:.0f}", kPullEdges[bin + 1U]) : std::string{"inf"};
      longitudinal += fmt::format(" [{:.0f},{}) {}", kPullEdges[bin], upper, longitudinalPulls[bin]);
      curvature += fmt::format(" [{:.0f},{}) {}", kPullEdges[bin], upper, curvaturePulls[bin]);
    }
    std::string origins;
    for (std::size_t bin = 0; bin < originAll.size(); ++bin) {
      const auto upper = bin + 1U < originAll.size() ? fmt::format("{:.0f}", 1.e4 * static_cast<double>(kOriginEdges[bin + 1U])) : std::string{"inf"};
      origins += fmt::format(" [{:.0f},{}) {}/{}", 1.e4 * static_cast<double>(kOriginEdges[bin]), upper, originCurvature[bin], originAll[bin]);
    }
    LOGP(info, "    origin radius [um] of the particle, pairs lost to the curvature cut / all:{}", origins);
    LOGP(info, "    longitudinal pulls of those lost to it:{}", longitudinal);
    LOGP(info, "    curvature pulls beyond the bound of those lost to it:{}", curvature);
    // The MC hits' verdict on each lost pair: the particle scattered hard
    // somewhere up to the target (a kink of 5 theta0 or more, theta0 of the
    // material actually crossed where a LUT is loaded), a cluster of
    // the pair sits off its crossing, the particle had slowed below the
    // floor by the source -- or none of these, and the cut itself lost it.
    if (crossings && scattering != nullptr && !lostSpatial.empty()) {
      std::unordered_map<std::uint64_t, std::vector<Frame::Index>> rowsOf;
      for (const auto& spatial : lostSpatial) {
        rowsOf[truePairs[spatial.pair].particle];
      }
      for (Frame::Index row = 0U; row < frame.size() && row < rows.size(); ++row) {
        for (std::uint8_t i = 0; i < rows[row].count; ++i) {
          if (const auto found = rowsOf.find(rows[row].particles[i]); found != rowsOf.end()) {
            found->second.push_back(row);
          }
        }
      }
      std::unordered_map<std::uint64_t, std::vector<TrueCrossingReport>> reportsOf;
      for (auto& [id, particleRows] : rowsOf) {
        reportsOf[id] = trueAnatomy(frame, *scattering, bzKiloGauss, crossings(id, particleRows), material);
      }
      enum Verdict : std::uint8_t { kVerdictKink,
                                    kVerdictOff,
                                    kVerdictSoft,
                                    kVerdictDrop,
                                    kVerdictDisplaced,
                                    kVerdictNone,
                                    kVerdictNoHits,
                                    kVerdicts };
      constexpr std::array<const char*, kVerdicts> kVerdictNames{
        "kink >= 5 theta0", "cluster off", "pT below the floor", "pT drops > 10 % on the link", "origin beyond 50 um", "none", "no hits"};
      constexpr std::array<const char*, kSpeciesCount> kSpeciesNames{"e", "mu", "pi", "K", "p", "other"};
      std::array<std::array<std::uint64_t, kSpeciesCount>, 2> noneBySpecies{};
      std::array<std::array<std::array<std::uint64_t, kVerdicts>, kPullEdges.size()>, 2> byPull{}; // curvature, longitudinal
      std::array<std::array<std::uint64_t, nLinks<Det>>, 2> noneByLink{};
      std::array<std::array<std::uint64_t, Frame::kLayerCount>, 2> kinkAtScatterer{}; // the kinked ones: the layer that scattered them hardest
      for (const auto& spatial : lostSpatial) {
        const auto& truth = truePairs[spatial.pair];
        const auto from = static_cast<std::uint32_t>(Det::layerLinks[truth.link].from);
        const auto to = static_cast<std::uint32_t>(Det::layerLinks[truth.link].to);
        const auto which = spatial.reason == PairLoss::kCurvature ? 0U : 1U;
        Verdict verdict = kVerdictNoHits;
        const auto found = reportsOf.find(truth.particle);
        if (found != reportsOf.end() && !found->second.empty()) {
          verdict = kVerdictNone;
          double hardest = 0.;
          std::uint32_t hardestAt = 0U;
          bool off = false;
          auto ptAtSource = static_cast<double>(truth.pt);
          double ptAtTarget = ptAtSource;
          int pdg = 0;
          for (const auto& report : found->second) {
            pdg = report.pdg;
            if (report.layer > to) {
              break;
            }
            const double kink = std::isfinite(report.kinkPath) ? report.kinkPath : report.kink;
            if (std::isfinite(kink) && kink > hardest) {
              hardest = kink;
              hardestAt = report.scatterer;
            }
            if ((report.layer == from || report.layer == to) && report.offsetCM > kOffsetCM) {
              off = true;
            }
            if (report.layer == from) {
              ptAtSource = report.pt;
            }
            if (report.layer == to) {
              ptAtTarget = report.pt;
            }
          }
          if (hardest >= kHardKink) {
            verdict = kVerdictKink;
            ++kinkAtScatterer[which][hardestAt];
          } else if (off) {
            verdict = kVerdictOff;
          } else if (ptAtSource < minPtGeV) {
            verdict = kVerdictSoft;
          } else if (ptAtTarget < 0.9 * ptAtSource) {
            verdict = kVerdictDrop;
          } else if (truth.originR > 0.005F) {
            verdict = kVerdictDisplaced;
          } else {
            ++noneByLink[which][truth.link];
            ++noneBySpecies[which][speciesOf(pdg)];
          }
        }
        ++byPull[which][spatial.pullBin][verdict];
      }
      constexpr std::array<const char*, 2> kCutNames{"curvature", "longitudinal"};
      for (std::size_t which = 0; which < 2U; ++which) {
        std::string line;
        for (std::size_t bin = 0; bin < kPullEdges.size(); ++bin) {
          const auto upper = bin + 1U < kPullEdges.size() ? fmt::format("{:.0f}", kPullEdges[bin + 1U]) : std::string{"inf"};
          line += fmt::format(" [{:.0f},{}):", kPullEdges[bin], upper);
          for (std::size_t verdict = 0; verdict < kVerdicts; ++verdict) {
            line += fmt::format(" {} {}", kVerdictNames[verdict], byPull[which][bin][verdict]);
          }
          line += ";";
        }
        LOGP(info, "    the MC hits on the pairs the {} cut lost, by pull:{}", kCutNames[which], line);
        std::string links;
        for (std::size_t link = 0; link < nLinks<Det>; ++link) {
          if (noneByLink[which][link] != 0U) {
            links += fmt::format(" {}->{}: {}", Det::layerLinks[link].from, Det::layerLinks[link].to, noneByLink[which][link]);
          }
        }
        std::string scatterers;
        for (std::size_t layer = 0; layer < Frame::kLayerCount; ++layer) {
          if (kinkAtScatterer[which][layer] != 0U) {
            scatterers += fmt::format(" {}: {}", layer, kinkAtScatterer[which][layer]);
          }
        }
        std::string species;
        for (std::size_t i = 0; i < kSpeciesNames.size(); ++i) {
          species += fmt::format(" {} {}", kSpeciesNames[i], noneBySpecies[which][i]);
        }
        LOGP(info, "      none of these, by the link:{}, by species:{}; the kinked, by the layer that scattered them hardest:{}", links, species, scatterers);
      }
    }
  }
}

/// One iteration of a tracker's own run against the MC truth, slice by slice
/// as the run's observers see them: its
/// tracklets, triplets and neighbours against the MC-true pairs, triplets and
/// neighbours above the iteration's pT_min on the clusters it found free; the
/// lost pairs' anatomy. A
/// cluster is every particle's it carries among its labels (TruthRow), the
/// rule of the comparison and the stage funnel.
class StageTruth
{
 public:
  /// At the iteration's first slice, with the used mask it read.
  StageTruth(const Frame& frame, const Controller& tracker, std::size_t iteration, const RecoSettings& param, const Propagator& propagator,
             std::span<const TruePair> allPairs, std::span<const TruthRow> rows, const ParticleTable& particles, const ScatteringTable<Det>& scattering,
             const o2::base::MatLayerCylSet* material, const TrueCrossings& crossings)
    : mFrame{&frame},
      mParam{&param},
      mIteration{iteration},
      mSettings{workflow::iterationSettings(param, static_cast<int>(iteration))},
      mParameters{tracker.iteration(iteration)},
      mBeamLine{tracker.beamLine()},
      mBzKiloGauss{propagator.bzKiloGauss},
      mRows{rows},
      mParticles{&particles},
      mScattering{&scattering},
      mMaterial{material},
      mCrossings{&crossings},
      mUsed(tracker.used().begin(), tracker.used().end()),
      mInitialPredictionMisses{tracker.stagingPredictionMisses()}
  {
    auto& summary = mSummary;
    for (const auto& pair : allPairs) {
      if (static_cast<double>(pair.pt) >= mSettings.windows.minPtGeV && mUsed[pair.key >> 32U] == 0U && mUsed[pair.key & 0xffffffffU] == 0U) {
        mTruePairs.push_back(pair);
      }
    }
    const auto& truePairs = mTruePairs;
    for (const auto& pair : truePairs) {
      ++summary.truePairs[pair.link];
      ++summary.truthBreakdown[pair.link][binOf(pair.pt, kPtEdges)][binOf(pair.cotTheta, kCotThetaEdges)];
    }
    // Every MC-true triplet, as the indices of its two pairs, one particle's,
    // in truePairs.
    for (std::size_t first = 0; first < truePairs.size(); ++first) {
      const auto middle = static_cast<std::uint32_t>(truePairs[first].key & 0xffffffffU);
      const auto next = std::ranges::lower_bound(truePairs, pairKey(middle, 0U), {}, &TruePair::key);
      for (auto second = next; second != truePairs.end() && (second->key >> 32U) == middle; ++second) {
        if (second->particle == truePairs[first].particle) {
          mTrueTriplets.emplace_back(first, static_cast<std::size_t>(second - truePairs.begin()));
        }
      }
    }
    summary.trueTriplets = mTrueTriplets.size();
    mPairKept.assign(truePairs.size(), false);
    // Every MC-true neighbour, as the indices of its two triplets in
    // trueTriplets; both lists are in lexicographic order.
    for (std::size_t first = 0; first < mTrueTriplets.size(); ++first) {
      const auto shared = mTrueTriplets[first].second;
      const auto next = std::ranges::lower_bound(mTrueTriplets, std::pair<std::size_t, std::size_t>{shared, 0U});
      for (auto second = next; second != mTrueTriplets.end() && second->first == shared; ++second) {
        mTrueChains.emplace_back(first, static_cast<std::size_t>(second - mTrueTriplets.begin()));
      }
    }
    summary.trueNeighbours = mTrueChains.size();
    mTrueTripletStates.assign(mTrueTriplets.size(), TrueTripletState{});
    mChainKept.assign(mTrueChains.size(), false);
    if (!mRows.empty()) {
      LOGP(info, "iteration {} in the run: {} MC-true pairs above {} GeV on the clusters it found free", mIteration, truePairs.size(),
           mSettings.windows.minPtGeV);
    }
  }

  void observeSlice(const Controller& tracker, const Controller::SliceObservation& observation)
  {
    const auto& parameters = mParameters;
    const auto& truePairs = mTruePairs;
    const auto& trueTriplets = mTrueTriplets;
    const auto& trueChains = mTrueChains;
    auto& summary = mSummary;
    auto& pairKept = mPairKept;
    auto& trueTripletStates = mTrueTripletStates;
    auto& chainKept = mChainKept;
    const auto& slice = *observation.slice;
    const auto& run = *observation.result;
    const auto truePairIndex = [&](Frame::Index source, Frame::Index target) {
      const auto key = pairKey(source, target);
      const auto found = std::ranges::lower_bound(truePairs, key, {}, &TruePair::key);
      return found != truePairs.end() && found->key == key ? static_cast<std::size_t>(found - truePairs.begin()) : truePairs.size();
    };
    const auto trueTripletIndex = [&](std::size_t firstPair, std::size_t secondPair) {
      const std::pair<std::size_t, std::size_t> key{firstPair, secondPair};
      const auto found = std::ranges::lower_bound(trueTriplets, key);
      return found != trueTriplets.end() && *found == key ? static_cast<std::size_t>(found - trueTriplets.begin()) : trueTriplets.size();
    };
    // Core sources only: the halo rows are worked within the slice gate so the
    // core's output does not depend on the slicing, and counting them would
    // count them once per slice that reaches them.
    const auto& sources = tracker.sources();
    const auto counts = run.tracklets->counts();
    for (Tracklets::Index element = 0U; element < run.tracklets->elementCount(); ++element) {
      const auto at = sources.row(element);
      const auto core = slice.core[at.layer];
      if (at.row < core.begin || at.row >= core.end) {
        continue;
      }
      ++summary.sources;
      ++summary.countHistogram[countBin(counts[element])];
      summary.maxPerSource = std::max(summary.maxPerSource, counts[element]);
    }
    const auto tracklets = run.tracklets->view();
    for (Tracklets::Index i = 0U; i < run.tracklets->size(); ++i) {
      const auto link = static_cast<std::size_t>(tracklets[i].link);
      const LayerID sourceLayer = parameters.linkSources[link];
      const auto core = slice.core[sourceLayer];
      const bool coreSource = tracklets[i].source >= core.begin && tracklets[i].source < core.end;
      if (!coreSource) {
        ++summary.haloSourceTracklets;
        continue;
      }
      ++summary.tracklets;
      const auto key = pairKey(tracklets[i].source, tracklets[i].target);
      const auto found = std::ranges::lower_bound(truePairs, key, {}, &TruePair::key);
      if (found != truePairs.end() && found->key == key) {
        ++summary.truePairsKept[link];
        ++summary.truthBreakdownKept[link][binOf(found->pt, kPtEdges)][binOf(found->cotTheta, kCotThetaEdges)];
        pairKept[static_cast<std::size_t>(found - truePairs.begin())] = true;
      }
    }

    // Core-owned triplets only -- first measurement in the core -- which
    // partitions them over the slices.
    const auto& triplets = *run.triplets;
    summary.firstTracklets += triplets.elementCount();
    const auto tripletView = triplets.view();
    for (Tracklets::Index i = 0U; i < triplets.size(); ++i) {
      const auto first = tripletView.first[i];
      const auto second = tripletView.second[i];
      const auto firstCore = slice.core[parameters.linkSources[tracklets[first].link]];
      if (tracklets[first].source < firstCore.begin || tracklets[first].source >= firstCore.end) {
        ++summary.haloTriplets;
        continue;
      }
      ++summary.triplets;
      const bool kinkedTriplet = tripletView.kinked[i] != 0U;
      summary.tripletsKinked += kinkedTriplet ? 1U : 0U;
      const auto firstPair = std::ranges::lower_bound(truePairs, pairKey(tracklets[first].source, tracklets[first].target), {}, &TruePair::key);
      const auto secondPair = std::ranges::lower_bound(truePairs, pairKey(tracklets[second].source, tracklets[second].target), {}, &TruePair::key);
      const bool trueFirst = firstPair != truePairs.end() && firstPair->key == pairKey(tracklets[first].source, tracklets[first].target);
      const bool trueSecond = secondPair != truePairs.end() && secondPair->key == pairKey(tracklets[second].source, tracklets[second].target);
      if (!trueFirst || !trueSecond) {
        classify({tracklets[first].source, tracklets[first].target, tracklets[second].target}, summary.notTrueCombinatorial, summary.notTrueSecondary,
                 summary.notTrueSoftPrimary, summary.notTrueOtherPrimary);
        continue;
      }
      ++summary.trueTripletsKept;
      summary.trueTripletsKinked += kinkedTriplet ? 1U : 0U;
      ++summary.trueTripletsKeptByPt[binOf(firstPair->pt, kPtEdges)];
      const auto trueIndex =
        trueTripletIndex(static_cast<std::size_t>(firstPair - truePairs.begin()), static_cast<std::size_t>(secondPair - truePairs.begin()));
      if (trueIndex < trueTriplets.size()) {
        trueTripletStates[trueIndex] = {.kept = true, .level = tripletView.level[i]};
      }
      // The fit's curvature is the 3D one, 1/R at the momentum p.
      const double momentum = static_cast<double>(firstPair->pt) * std::hypot(1., static_cast<double>(firstPair->cotTheta));
      const double truth = 1. / curvatureRadius(momentum, mBzKiloGauss);
      const double fitted = std::abs(static_cast<double>(tripletView.curvature[i]));
      addCurvaturePull(summary, firstPair->pt, (fitted - truth) / std::sqrt(static_cast<double>(tripletView.varCurvature[i])), (fitted - truth) / truth);
    }

    // Neighbours: each edge is owned by the slice whose core holds its first
    // measurement a; the predecessors per triplet by the triplet's owner.
    const auto& neighbours = *run.neighbours;
    const auto* const from = neighbours.view().from;
    const auto neighbourOffsets = neighbours.offsets();
    summary.tripletsWorked += triplets.size();
    // The triplets with predecessors.
    for (std::size_t triplet = 0; triplet + 1U < neighbourOffsets.size(); ++triplet) {
      summary.chained += neighbourOffsets[triplet + 1U] != neighbourOffsets[triplet] ? 1U : 0U;
    }
    ++summary.slices;
    const auto inCore = [&](Tracklets::Index tracklet) {
      const auto core = slice.core[parameters.linkSources[tracklets[tracklet].link]];
      return tracklets[tracklet].source >= core.begin && tracklets[tracklet].source < core.end;
    };
    for (Tracklets::Index to = 0U; to < triplets.size(); ++to) {
      const auto shared = tripletView.first[to];
      if (inCore(shared)) {
        const auto count = neighbourOffsets[to + 1U] - neighbourOffsets[to];
        ++summary.neighbourHistogram[countBin(count)];
        summary.maxNeighbours = std::max(summary.maxNeighbours, count);
      }
      for (auto edge = neighbourOffsets[to]; edge < neighbourOffsets[to + 1U]; ++edge) {
        const auto first = tripletView.first[from[edge]];
        if (!inCore(first)) {
          continue;
        }
        ++summary.neighbours;
        const auto a = tracklets[first].source;
        const auto b = tracklets[shared].source;
        const auto c = tracklets[shared].target;
        const auto d = tracklets[tripletView.second[to]].target;
        const auto ab = truePairIndex(a, b);
        const auto bc = truePairIndex(b, c);
        const auto cd = truePairIndex(c, d);
        const auto firstTriplet = ab < truePairs.size() && bc < truePairs.size() ? trueTripletIndex(ab, bc) : trueTriplets.size();
        const auto secondTriplet = bc < truePairs.size() && cd < truePairs.size() ? trueTripletIndex(bc, cd) : trueTriplets.size();
        if (firstTriplet < trueTriplets.size() && secondTriplet < trueTriplets.size()) {
          const std::pair<std::size_t, std::size_t> key{firstTriplet, secondTriplet};
          const auto chain = std::ranges::lower_bound(trueChains, key);
          if (chain != trueChains.end() && *chain == key) {
            chainKept[static_cast<std::size_t>(chain - trueChains.begin())] = true;
            continue;
          }
        }
        classify({a, b, c, d}, summary.neighboursNotTrueCombinatorial, summary.neighboursNotTrueSecondary, summary.neighboursNotTrueSoftPrimary,
                 summary.neighboursNotTrueOtherPrimary);
      }
    }
  }

  /// The iteration's end: its report, and its summary.
  [[nodiscard]] TrackletRunSummary finish(const Controller& tracker, const Controller::IterationObservation& observed)
  {
    const auto& frame = *mFrame;
    const auto& param = *mParam;
    const auto& parameters = mParameters;
    const auto& settings = mSettings;
    const auto& beamLine = mBeamLine;
    const auto& truePairs = mTruePairs;
    const auto& trueTriplets = mTrueTriplets;
    const auto& trueChains = mTrueChains;
    const auto& pairKept = mPairKept;
    const auto& trueTripletStates = mTrueTripletStates;
    const auto& chainKept = mChainKept;
    const auto& rows = mRows;
    const auto& crossings = *mCrossings;
    const auto& scattering = *mScattering;
    const auto* material = mMaterial;
    const auto vertexSource = parameters.vertexSource;
    const auto initialPredictionMisses = mInitialPredictionMisses;
    const std::int64_t requiredHalo = sliceReach(frameBounds(frame), frame.grid().geometry()).haloPS();
    auto& summary = mSummary;
    if (vertexSource == VertexSource::kBeamLine) {
      reportLostPairs(frame, parameters, beamLine, truePairs, pairKept, rows, crossings, &scattering, mBzKiloGauss, settings.windows.minPtGeV,
                      material);
    }

    // The true neighbours whose triplets were both kept, and of them the ones
    // found.
    for (std::size_t chain = 0; chain < trueChains.size(); ++chain) {
      const auto& [firstTriplet, secondTriplet] = trueChains[chain];
      if (!trueTripletStates[firstTriplet].kept || !trueTripletStates[secondTriplet].kept) {
        continue;
      }
      const auto& ab = truePairs[trueTriplets[firstTriplet].first];
      const auto ptBin = binOf(ab.pt, kPtEdges);
      const auto cotBin = binOf(ab.cotTheta, kCotThetaEdges);
      ++summary.trueNeighboursTripletsKept;
      ++summary.trueNeighboursTripletsKeptByPt[ptBin];
      ++summary.trueNeighboursTripletsKeptByCot[cotBin];
      if (chainKept[chain]) {
        ++summary.trueNeighboursKept;
        ++summary.trueNeighboursKeptByPt[ptBin];
        ++summary.trueNeighboursKeptByCot[cotBin];
      }
    }

    // Per MC track: the longest chain its kept true triplets allow, through
    // the true neighbours between them, against the highest level they got.
    {
      std::vector<std::uint32_t> longest(trueTriplets.size(), 0U);
      for (std::size_t triplet = 0; triplet < trueTriplets.size(); ++triplet) {
        longest[triplet] = trueTripletStates[triplet].kept ? 1U : 0U;
      }
      for (bool changed = true; changed;) {
        changed = false;
        for (const auto& [first, second] : trueChains) {
          if (longest[first] != 0U && longest[second] != 0U && longest[first] + 1U > longest[second]) {
            longest[second] = longest[first] + 1U;
            changed = true;
          }
        }
      }
      std::vector<std::pair<std::uint64_t, std::pair<std::uint32_t, std::uint32_t>>> tracks; // track, (longest, level)
      for (std::size_t triplet = 0; triplet < trueTriplets.size(); ++triplet) {
        if (trueTripletStates[triplet].kept) {
          tracks.push_back({truePairs[trueTriplets[triplet].first].particle, {longest[triplet], trueTripletStates[triplet].level}});
        }
      }
      std::ranges::sort(tracks);
      for (std::size_t begin = 0; begin < tracks.size();) {
        std::size_t end = begin;
        std::uint32_t chain = 0U;
        std::uint32_t level = 0U;
        for (; end < tracks.size() && tracks[end].first == tracks[begin].first; ++end) {
          chain = std::max(chain, tracks[end].second.first);
          level = std::max(level, tracks[end].second.second);
        }
        const auto bin = std::min<std::size_t>(chain, summary.tracksByChain.size() - 1U);
        ++summary.tracksByChain[bin];
        summary.tracksReachingChain[bin] += level == chain ? 1U : 0U;
        summary.tracksBelowChain[bin] += level < chain ? 1U : 0U;
        summary.tracksAboveChain[bin] += level > chain ? 1U : 0U;
        begin = end;
      }
    }

    for (const auto& [first, second] : trueTriplets) {
      if (pairKept[first] && pairKept[second]) {
        ++summary.trueTripletsPairsKept;
        ++summary.trueTripletsPairsKeptByPt[binOf(truePairs[first].pt, kPtEdges)];
      }
    }
    summary.tripletHighWater = tracker.tripletHighWater();
    summary.neighbourHighWater = tracker.neighbourHighWater();
    summary.trackletHighWater = tracker.trackletHighWater();
    summary.stagingHighWater = tracker.stagingHighWater();
    summary.stagingPredictionMisses = tracker.stagingPredictionMisses() - initialPredictionMisses;
    summary.poolHighWater = tracker.memoryHighWater();

    const auto label = fmt::format("iteration {}, {}", mIteration, vertexSource == VertexSource::kBeamLine ? "beam line" : "vertex list");
    const double trackletSeconds = observed.summary->seconds(profile::Phase::kTracklets);
    const double timePerSourceNS = summary.sources == 0U ? 0. : 1.e9 * trackletSeconds / static_cast<double>(summary.sources);
    LOGP(info, "tracklets, {}: {} sources, {} accepted, {:.3f} s, {:.1f} ns/source", label, summary.sources, summary.tracklets, trackletSeconds,
         timePerSourceNS);
    const auto allTracklets = summary.tracklets + summary.haloSourceTracklets;
    const double haloTrackletFraction = allTracklets == 0U ? 0. : static_cast<double>(summary.haloSourceTracklets) / static_cast<double>(allTracklets);
    LOGP(info, "  halo-source tracklets, within the slice gate: {}/{} ({:.3f}%)", summary.haloSourceTracklets, allTracklets, 100. * haloTrackletFraction);
    LOGP(info, "  halo configured/required/used {}/{}/{} ps", param.sliceHaloPS, requiredHalo, observed.summary->haloPS);
    LOGP(info,
         "  arena/staging high-water {} / {} records ({:.1f} / {:.1f} MB), "
         "pool {:.1f} MB, {} staging prediction misses",
         summary.trackletHighWater, summary.stagingHighWater, static_cast<double>(summary.trackletHighWater * kTrackletBytes) / 1.e6,
         static_cast<double>(summary.stagingHighWater * kTrackletBytes) / 1.e6, static_cast<double>(summary.poolHighWater) / 1.e6,
         summary.stagingPredictionMisses);
    LOGP(info, "  accepted per source, max {}: {}", summary.maxPerSource, histogramText(summary.countHistogram));
    std::uint64_t trueTotal = 0U;
    std::uint64_t trueKept = 0U;
    for (std::size_t link = 0; link < nLinks<Det>; ++link) {
      trueTotal += summary.truePairs[link];
      trueKept += summary.truePairsKept[link];
    }
    if (trueTotal != 0U) {
      LOGP(info, "  MC-true primary pairs: {}/{} kept ({:.3f}%)", trueKept, trueTotal, percent(trueKept, trueTotal));
      for (std::size_t link = 0; link < nLinks<Det>; ++link) {
        LOGP(info, "    {:>2} -> {:<2}: {}/{} ({:.3f}%)", Det::layerLinks[link].from, Det::layerLinks[link].to, summary.truePairsKept[link],
             summary.truePairs[link], percent(summary.truePairsKept[link], summary.truePairs[link]));
      }
    }

    const double tripletSeconds = observed.summary->seconds(profile::Phase::kTriplets);
    const double timePerFirstNS = summary.firstTracklets == 0U ? 0. : 1.e9 * tripletSeconds / static_cast<double>(summary.firstTracklets);
    LOGP(info, "triplets, {}: {} from {} first tracklets, {:.3f} s, {:.1f} ns/first tracklet", label, summary.triplets, summary.firstTracklets,
         tripletSeconds, timePerFirstNS);
    LOGP(info, "  arena high-water {} records ({:.1f} MB); halo-owned triplets, within the slice gate: {} ({:.3f}%)", summary.tripletHighWater,
         static_cast<double>(summary.tripletHighWater * sizeof(Triplet)) / 1.e6, summary.haloTriplets,
         percent(summary.haloTriplets, summary.triplets + summary.haloTriplets));
    if (summary.trueTriplets != 0U) {
      LOGP(info, "  MC-true triplets: {}/{} kept ({:.3f}%); with both pairs kept as tracklets {}/{} ({:.3f}%)", summary.trueTripletsKept, summary.trueTriplets,
           percent(summary.trueTripletsKept, summary.trueTriplets), summary.trueTripletsKept, summary.trueTripletsPairsKept,
           percent(summary.trueTripletsKept, summary.trueTripletsPairsKept));
      if (summary.tripletsKinked != 0U) {
        LOGP(info, "  kept only with the kink allowance: {} triplets ({:.3f}%), {} of them MC-true ({:.3f}% of the kinked; {:.3f}% of the true kept)",
             summary.tripletsKinked, percent(summary.tripletsKinked, summary.triplets), summary.trueTripletsKinked,
             percent(summary.trueTripletsKinked, summary.tripletsKinked), percent(summary.trueTripletsKinked, summary.trueTripletsKept));
      }
      const auto notTrue = summary.triplets - summary.trueTripletsKept;
      LOGP(info, "  not MC-true: {}/{} ({:.3f}%)", notTrue, summary.triplets, percent(notTrue, summary.triplets));
      LOGP(info,
           "    one particle: secondary {} ({:.2f}%), primary below pT_min {} ({:.2f}%), other primary {} ({:.2f}%); "
           "combinatorial {} ({:.2f}%)",
           summary.notTrueSecondary, percent(summary.notTrueSecondary, summary.triplets), summary.notTrueSoftPrimary,
           percent(summary.notTrueSoftPrimary, summary.triplets), summary.notTrueOtherPrimary, percent(summary.notTrueOtherPrimary, summary.triplets),
           summary.notTrueCombinatorial, percent(summary.notTrueCombinatorial, summary.triplets));
      std::string byPt;
      for (std::size_t bin = 1; bin < kPtBins; ++bin) {
        byPt +=
          fmt::format(" {}-{}: {:.2f}%", kPtEdges[bin], kPtEdges[bin + 1U], percent(summary.trueTripletsKeptByPt[bin], summary.trueTripletsPairsKeptByPt[bin]));
      }
      LOGP(info, "  kept with both pairs kept, by pT [GeV]:{}", byPt);
      logTripletCurvature(summary);
    }

    const double neighbourSeconds = observed.summary->seconds(profile::Phase::kNeighbours);
    const double neighbourPerTripletNS = summary.tripletsWorked == 0U ? 0. : 1.e9 * neighbourSeconds / static_cast<double>(summary.tripletsWorked);
    LOGP(info, "neighbours, {}: {} from {} triplets worked, {:.3f} s, {:.1f} ns/triplet", label, summary.neighbours, summary.tripletsWorked,
         neighbourSeconds, neighbourPerTripletNS);
    LOGP(info, "  arena high-water {} records ({:.1f} MB); swept {:.2f}% of the triplets", summary.neighbourHighWater,
         static_cast<double>(summary.neighbourHighWater * sizeof(Neighbour)) / 1.e6, percent(summary.chained, summary.tripletsWorked));
    LOGP(info, "  predecessors per triplet, max {}: {}", summary.maxNeighbours, histogramText(summary.neighbourHistogram));
    if (summary.trueNeighbours != 0U) {
      LOGP(info, "  MC-true neighbours: {}/{} kept; with both triplets kept {}/{} ({:.3f}%)", summary.trueNeighboursKept, summary.trueNeighbours,
           summary.trueNeighboursKept, summary.trueNeighboursTripletsKept, percent(summary.trueNeighboursKept, summary.trueNeighboursTripletsKept));
      std::string byPt;
      for (std::size_t bin = 1; bin < kPtBins; ++bin) {
        byPt += fmt::format(" [{:.2f},{:.2f}) {:.2f}%", kPtEdges[bin], kPtEdges[bin + 1U],
                            percent(summary.trueNeighboursKeptByPt[bin], summary.trueNeighboursTripletsKeptByPt[bin]));
      }
      LOGP(info, "    kept with both triplets kept, by pT [GeV]:{}", byPt);
      std::string byCot;
      for (std::size_t bin = 0; bin < kCotThetaBins; ++bin) {
        if (summary.trueNeighboursTripletsKeptByCot[bin] != 0U) {
          byCot += fmt::format(" [{:.2f},{:.2f}) {:.2f}%", kCotThetaEdges[bin], kCotThetaEdges[bin + 1U],
                               percent(summary.trueNeighboursKeptByCot[bin], summary.trueNeighboursTripletsKeptByCot[bin]));
        }
      }
      LOGP(info, "    by cot(theta):{}", byCot);
      const auto notTrue = summary.neighbours - summary.trueNeighboursKept;
      LOGP(info,
           "  not MC-true: {}/{} ({:.3f}%); one particle: secondary {:.2f}%, primary below pT_min {:.2f}%, other primary {:.2f}%; "
           "combinatorial {:.2f}%",
           notTrue, summary.neighbours, percent(notTrue, summary.neighbours), percent(summary.neighboursNotTrueSecondary, summary.neighbours),
           percent(summary.neighboursNotTrueSoftPrimary, summary.neighbours), percent(summary.neighboursNotTrueOtherPrimary, summary.neighbours),
           percent(summary.neighboursNotTrueCombinatorial, summary.neighbours));
      std::string byChain;
      for (std::size_t chain = 1; chain < summary.tracksByChain.size(); ++chain) {
        if (summary.tracksByChain[chain] != 0U) {
          byChain += fmt::format(" {}: {} ({:.2f}% reach it, {:.2f}% below, {:.2f}% above)", chain, summary.tracksByChain[chain],
                                 percent(summary.tracksReachingChain[chain], summary.tracksByChain[chain]),
                                 percent(summary.tracksBelowChain[chain], summary.tracksByChain[chain]),
                                 percent(summary.tracksAboveChain[chain], summary.tracksByChain[chain]));
        }
      }
      LOGP(info, "  MC tracks by their longest chain of kept true triplets, and the level they reach:{}", byChain);
    }
    return summary;
  }

 private:
  /// Counts a stage's element that is no MC-true one by what made it: no
  /// particle on every row of \a of, a secondary's, a primary's below the
  /// iteration's pT_min, another primary's.
  void classify(std::initializer_list<Frame::Index> of, std::uint64_t& combinatorial, std::uint64_t& secondary, std::uint64_t& softPrimary,
                std::uint64_t& otherPrimary) const
  {
    if (mRows.empty()) {
      return;
    }
    const auto id = goodParticle(mRows, std::span<const Frame::Index>{of.begin(), of.size()});
    const auto truth = id == kNoTrackParticle ? mParticles->end() : mParticles->find(id);
    if (truth == mParticles->end()) {
      ++combinatorial;
    } else if (!truth->second.primary) {
      ++secondary;
    } else if (static_cast<double>(truth->second.pt) < mSettings.windows.minPtGeV) {
      ++softPrimary;
    } else {
      ++otherPrimary;
    }
  }

  const Frame* mFrame;
  const RecoSettings* mParam;
  std::size_t mIteration;
  IterationSettings mSettings;
  IterationParameters<Det> mParameters;
  VertexSeed mBeamLine;
  double mBzKiloGauss;
  std::span<const TruthRow> mRows;
  const ParticleTable* mParticles;
  const ScatteringTable<Det>* mScattering;
  const o2::base::MatLayerCylSet* mMaterial;
  const TrueCrossings* mCrossings;
  std::vector<std::uint8_t> mUsed; ///< the mask the iteration read
  std::size_t mInitialPredictionMisses;
  std::vector<TruePair> mTruePairs; ///< above pT_min, on free clusters
  std::vector<std::pair<std::size_t, std::size_t>> mTrueTriplets;
  std::vector<std::pair<std::size_t, std::size_t>> mTrueChains;
  std::vector<bool> mPairKept;
  std::vector<TrueTripletState> mTrueTripletStates;
  std::vector<bool> mChainKept;
  TrackletRunSummary mSummary;
};

/// The MC-true chains whose fit fails: where
/// -- the pass, the propagation status, the step as layer to layer -- and
/// what -- pT, eta, species, origin, the chain's layers -- with the chain
/// reference's pT against the true one, and every chain printed up to a cap.
struct TrueRoadFitFailure {
  std::uint64_t id = 0U; ///< the particle
  std::vector<std::uint32_t> rows;
  TrackFitResult fit;
};

void reportTrueRoadFitFailures(const Frame& frame, std::span<const TrueRoadFitFailure> failures, const ParticleTable& particles,
                               const TrackingGeometry<Det>& tracking, const Propagator& propagator, const TrueCrossings& crossings,
                               const ScatteringTable<Det>& scattering, const o2::base::MatLayerCylSet& material)
{
  if (failures.empty()) {
    return;
  }
  constexpr std::array<const char*, 3> kPasses{"none", "outward", "inward"};
  constexpr std::array<const char*, 8> kStatuses{"ok", "no intersection", "wrong direction", "out of bounds",
                                                 "grazing", "not finite", "too many turns", "stopped"};
  constexpr std::size_t kExamples = 12U;
  std::map<std::pair<std::uint8_t, std::uint8_t>, std::uint64_t> byPassStatus;
  std::map<std::tuple<std::uint8_t, int, int>, std::uint64_t> bySteps;
  std::map<int, std::uint64_t> bySpecies;
  std::map<std::size_t, std::uint64_t> byCount;
  std::array<std::uint64_t, 2> byPrimary{};
  // From the MC hits, in the chain's order: the particle loops -- it turns
  // more than half a turn between two crossings (|dz| pT / |pz| over its
  // radius at the first), a second pass --; it turns back otherwise -- its
  // momentum points inwards at a barrel crossing, or its pz changes sign --;
  // it loses more than a fifth of its momentum (brems); it scatters hard
  // (a kink of 5 theta0 or more); none of these; no hits.
  enum Kind : std::uint8_t { kLoops,
                             kTurnsBack,
                             kLoses,
                             kHard,
                             kNone,
                             kNoHits,
                             kKinds };
  constexpr std::array<const char*, kKinds> kKindNames{
    "loops (> half a turn between crossings)", "turns back", "loses > 20 % of p", "hard scatter (>= 5 theta0)", "none of these", "no hits"};
  std::array<std::uint64_t, kKinds> byKind{};
  std::size_t examples = 0U;
  for (const auto& failure : failures) {
    const auto& rows = failure.rows;
    const auto& fit = failure.fit;
    const int n = static_cast<int>(rows.size());
    const auto layerOf = [&](int index) { return static_cast<int>(tracking.surface(frame.measurement(rows[static_cast<std::size_t>(index)]).surface).layer); };
    const int to = fit.failedMeasurement;
    const int from = fit.failedPass == FitPass::kInward ? to + 1 : to - 1;
    const bool stepInChain = std::min(from, to) >= 0 && std::max(from, to) < n;
    ++byPassStatus[{static_cast<std::uint8_t>(fit.failedPass), static_cast<std::uint8_t>(fit.status)}];
    ++bySteps[{static_cast<std::uint8_t>(fit.failedPass), stepInChain ? layerOf(from) : -1, to >= 0 && to < n ? layerOf(to) : -1}];
    ++byCount[rows.size()];
    const auto found = particles.find(failure.id);
    const ParticleTruth truth = found == particles.end() ? ParticleTruth{} : found->second;
    ++bySpecies[truth.pdg];
    ++byPrimary[truth.primary ? 0U : 1U];
    Kind kind = kNoHits;
    std::string hits;
    if (crossings) {
      auto sorted = crossings(failure.id, std::span<const Frame::Index>{rows});
      if (sorted.size() == rows.size()) {
        std::vector<TrueCrossing> ordered;
        ordered.reserve(rows.size());
        for (const auto row : rows) {
          ordered.push_back(*std::ranges::find(sorted, row, &TrueCrossing::row));
        }
        const auto reports = trueAnatomy(frame, scattering, static_cast<double>(propagator.bzKiloGauss), ordered, &material);
        bool back = false;
        double maxTurn = 0.;
        const double pzFirst = ordered.front().momentum[2];
        for (std::size_t i = 0; i < ordered.size(); ++i) {
          const auto& crossing = ordered[i];
          const auto& x = crossing.position;
          const auto& p = crossing.momentum;
          back = back || (Detector<Det>::isBarrelLayer(static_cast<LayerID>(crossing.layer)) && (x[0] * p[0]) + (x[1] * p[1]) < 0.) || pzFirst * p[2] < 0.;
          if (i > 0U) {
            const auto& before = ordered[i - 1U];
            const double pt = std::hypot(before.momentum[0], before.momentum[1]);
            const double turn =
              std::abs(x[2] - before.position[2]) * pt / std::abs(before.momentum[2]) / curvatureRadius(pt, static_cast<double>(propagator.bzKiloGauss));
            maxTurn = std::isfinite(turn) ? std::max(maxTurn, turn) : maxTurn;
          }
          hits += fmt::format(" {}(p {:.3f}, pT {:.3f})", crossing.layer, std::hypot(p[0], p[1], p[2]), std::hypot(p[0], p[1]));
        }
        const auto& first = ordered.front().momentum;
        const auto& last = ordered.back().momentum;
        const double ratio = std::hypot(last[0], last[1], last[2]) / std::hypot(first[0], first[1], first[2]);
        double hardest = 0.;
        for (const auto& report : reports) {
          hardest = std::isfinite(report.kink) ? std::max(hardest, report.kink) : hardest;
        }
        hits += fmt::format("; p last/first {:.3f}, hardest kink {:.1f} theta0, largest turn {:.2f} rad", ratio, hardest, maxTurn);
        if (maxTurn > std::numbers::pi) {
          kind = kLoops;
        } else if (back) {
          kind = kTurnsBack;
        } else if (ratio < 0.8) {
          kind = kLoses;
        } else if (hardest >= kHardKink) {
          kind = kHard;
        } else {
          kind = kNone;
        }
      }
    }
    ++byKind[kind];
    if (examples < kExamples) {
      ++examples;
      const auto seed = makeChainSeed(propagator, tracking, RowChain<Det>{.frame = &frame, .rows = rows.data()}, n);
      const auto qOverP = static_cast<double>(std::abs(seed.state.qOverP()));
      const double seedPt = qOverP > 0. ? std::abs(std::sin(static_cast<double>(seed.state.theta()))) / qOverP : 0.;
      std::string layers;
      for (int i = 0; i < n; ++i) {
        const auto m = frame.measurement(rows[static_cast<std::size_t>(i)]);
        layers += fmt::format(" {}(r {:.2f}, z {:.2f}, phi {:.3f})", layerOf(i), m.r, m.z, std::atan2(m.y, m.x));
      }
      LOGP(info,
           "    true chain {:#x}: {} {} at index {} of {}; pdg {} (mother {}, process {}), {} at ({:.2f}, {:.2f}, {:.2f}); MC pT {:.3f}, eta {:.3f}; seed pT "
           "{:.3f}, theta {:.4f}; {}; layers{}; hits{}",
           failure.id, kPasses[static_cast<std::size_t>(fit.failedPass)], kStatuses[static_cast<std::size_t>(fit.status)], to, n, truth.pdg, truth.motherPdg,
           static_cast<int>(truth.process), truth.primary ? "primary" : "secondary", truth.vx, truth.vy, truth.vz, truth.pt, truth.eta, seedPt,
           seed.state.theta(), kKindNames[kind], layers, hits);
    }
  }
  LOGP(info, "  true chains whose fit fails: {}; primary {}, secondary {}", failures.size(), byPrimary[0], byPrimary[1]);
  for (const auto& [key, count] : byPassStatus) {
    LOGP(info, "    {}, {}: {}", kPasses[key.first], kStatuses[key.second], count);
  }
  std::string steps;
  for (const auto& [key, count] : bySteps) {
    steps += fmt::format(" {} {}->{}: {};", kPasses[std::get<0>(key)], std::get<1>(key), std::get<2>(key), count);
  }
  LOGP(info, "    by the failing step (pass, layer to layer):{}", steps);
  std::string species;
  for (const auto& [pdg, count] : bySpecies) {
    species += fmt::format(" {}: {}", pdg, count);
  }
  LOGP(info, "    by species:{}", species);
  std::string counts;
  for (const auto& [count, chains] : byCount) {
    counts += fmt::format(" {}: {}", count, chains);
  }
  LOGP(info, "    by clusters:{}", counts);
  std::string kinds;
  for (std::size_t kind = 0; kind < kKinds; ++kind) {
    kinds += fmt::format("{} {}: {}", kind == 0U ? "" : ";", kKindNames[kind], byKind[kind]);
  }
  LOGP(info, "    by the MC hits:{}", kinds);
}

/// The road fit on the MC-true roads: every particle the
/// first iteration could find, by its longest true chain -- first-pass
/// clusters joined by true pairs, the innermost `maxChainLayers` of it --
/// fitted as the road fit and the refit fit it. Against the cut: how many fail
/// by length, next to what a calibrated chi2 loses at nSigma. Where the chi2
/// comes from: the outward chi2 per degree of freedom, and per layer the
/// median of each spatial update's increment against the median of its chi2
/// distribution -- the square of the scale its errors are off by, robust to a
/// tail -- and the fraction beyond the nSigma quantile. An increment counts
/// once its pass has fixed the helix (three earlier clusters, two dof); both
/// passes, so the inner layers are seen too.
void reportTrueRoadChi2(const Frame& frame, const RecoSettings& param, const Controller& tracker, std::span<const TruePair> truePairs,
                        const ParticleTable& particles, const TrackingGeometry<Det>& tracking, const o2::base::MatLayerCylSet& material,
                        const Propagator& propagator, const TrueCrossings& crossings, const ScatteringTable<Det>& scattering)
{
  constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  constexpr double kMedianChi2Two = 1.3862943611198906; // 2 ln 2
  const auto settings = workflow::iterationSettings(param, 0);
  const auto& parameters = tracker.iteration(0);
  const auto minClusters = static_cast<std::uint32_t>(std::max(settings.minClusters, 4));
  const double nSigma = settings.windows.nSigma;
  const double tailTwo = chiSquareQuantile(2, nSigma);

  // The longest true chain along each pair, one particle's, and each
  // particle's longest: its last pair and pT.
  const auto chains = trueChains(truePairs);
  std::unordered_map<std::uint64_t, std::pair<std::size_t, Scalar>> ends;
  for (std::size_t pair = 0; pair < truePairs.size(); ++pair) {
    auto [end, inserted] = ends.try_emplace(truePairs[pair].particle, pair, truePairs[pair].pt);
    if (!inserted && chains[pair].length > chains[end->second.first].length) {
      end->second.first = pair;
    }
  }

  struct ByLength {
    std::uint64_t chains = 0U;
    std::uint64_t failed = 0U;
    std::uint64_t fitFailed = 0U;
  };
  std::array<ByLength, kMaxClusters + 1U> byLength{};
  std::array<std::vector<float>, static_cast<std::size_t>(Detector<Det>::nLayers)> spatial;
  std::array<std::uint64_t, static_cast<std::size_t>(Detector<Det>::nLayers)> worstOfFailed{};
  double chi2Sum = 0.;
  double ndfSum = 0.;
  std::uint64_t kinked = 0U;                         ///< chains whose outward pass took a kink
  std::vector<float> chi2PerNdf;                     ///< every chain's outward chi2 / ndf
  std::array<std::uint64_t, 4> failedByBeyondGate{}; ///< failing chains by their outward steps beyond the kink gate: 0, 1, 2, 3+
  std::array<std::uint64_t, 4> failedByWorstShare{}; ///< failing chains by the worst step's share of their chi2: <25, 25-50, 50-75, >75 %
  std::array<std::uint64_t, 3> diffuseByKind{};      ///< failing chains without a step beyond the gate, by their layers: barrel only, discs only, both
  // Why the kink did not save a failing chain:
  // the budget spent before the miss; the miss on the last step, where no kink
  // is taken; a miss beyond the gate elsewhere
  // although kinks were left; no step beyond the gate. For each, how many the
  // inward pass -- which kinks one step later, at the last cluster too -- passes.
  enum Why : std::uint8_t { kWhySpent,
                            kWhyLastStep,
                            kWhyOtherMiss,
                            kWhyDiffuse,
                            kWhys };
  constexpr std::array<const char*, kWhys> kWhyNames{"the budget spent", "the miss on the last step", "a miss elsewhere", "no step beyond the gate"};
  std::array<std::uint64_t, kWhys> failedByWhy{};
  std::array<std::uint64_t, kWhys> inwardPassesByWhy{};
  std::array<std::uint64_t, kWhys> fullLengthByWhy{}; ///< of them, chains over every layer
  // The chains whose fit fails, for their anatomy below.
  std::vector<TrueRoadFitFailure> failures;
  std::vector<std::uint32_t> rows;
  for (const auto& [id, end] : ends) {
    if (static_cast<double>(end.second) < settings.windows.minPtGeV || chains[end.first].length < minClusters) {
      continue;
    }
    rows.clear();
    auto pair = end.first;
    rows.push_back(static_cast<std::uint32_t>(truePairs[pair].key & 0xffffffffU));
    for (; chains[pair].previous != TrueChain::kNoPair; pair = chains[pair].previous) {
      rows.push_back(static_cast<std::uint32_t>(truePairs[pair].key >> 32U));
    }
    rows.push_back(static_cast<std::uint32_t>(truePairs[pair].key >> 32U));
    std::ranges::reverse(rows);
    rows.resize(std::min(rows.size(), kMaxClusters));
    const int n = static_cast<int>(rows.size());
    std::array<float, kMaxClusters> outwardSpatial{};
    std::array<float, kMaxClusters> inwardSpatial{};
    const auto observe = [&](PropagationDirection direction, int index, Scalar increment) {
      (direction == PropagationDirection::kForward ? outwardSpatial : inwardSpatial)[static_cast<std::size_t>(index)] = increment;
    };
    const auto fit =
      fitMeasurements(propagator, tracking, material, RowChain<Det>{.frame = &frame, .rows = rows.data()}, n, kinkBudget(parameters), observe);
    auto& length = byLength[static_cast<std::size_t>(n)];
    ++length.chains;
    if (fit.status != PropagationStatus::kOk) {
      ++length.fitFailed;
      failures.push_back({.id = id, .rows = rows, .fit = fit});
      continue;
    }
    const int outwardNdf = (2 * n) - kNBoundParams - (2 * fit.outwardKinks);
    chi2Sum += static_cast<double>(fit.outwardChi2);
    ndfSum += static_cast<double>(outwardNdf);
    const auto layerOf = [&](int index) {
      return static_cast<std::size_t>(tracking.surface(frame.measurement(rows[static_cast<std::size_t>(index)]).surface).layer);
    };
    for (int index = 0; index < n; ++index) {
      const auto at = static_cast<std::size_t>(index);
      if (index >= 3) {
        spatial[layerOf(index)].push_back(outwardSpatial[at]);
      }
      if (index <= n - 4) {
        spatial[layerOf(index)].push_back(inwardSpatial[at]);
      }
    }
    kinked += fit.outwardKinks != 0U ? 1U : 0U;
    chi2PerNdf.push_back(fit.outwardChi2 / static_cast<float>(outwardNdf));
    if (!(fit.outwardChi2 <= roadChi2Cut(parameters, static_cast<std::uint32_t>(n), fit.outwardKinks))) {
      ++length.failed;
      int worst = 0;
      for (int index = 1; index < n; ++index) {
        if (outwardSpatial[static_cast<std::size_t>(index)] > outwardSpatial[static_cast<std::size_t>(worst)]) {
          worst = index;
        }
      }
      ++worstOfFailed[layerOf(worst)];
      std::size_t beyond = 0U;
      for (int index = 1; index < n; ++index) {
        beyond += outwardSpatial[static_cast<std::size_t>(index)] > parameters.kinkChi2Min ? 1U : 0U;
      }
      ++failedByBeyondGate[std::min<std::size_t>(beyond, 3U)];
      {
        Why why = kWhyDiffuse;
        if (beyond != 0U) {
          const bool lastBeyond = outwardSpatial[static_cast<std::size_t>(n - 1)] > parameters.kinkChi2Min;
          if (fit.outwardKinks >= parameters.maxKinks) {
            why = kWhySpent;
          } else if (lastBeyond && beyond == 1U) {
            why = kWhyLastStep;
          } else {
            why = kWhyOtherMiss;
          }
        }
        ++failedByWhy[why];
        inwardPassesByWhy[why] += fit.chi2 <= roadChi2Cut(parameters, static_cast<std::uint32_t>(n), fit.kinks) ? 1U : 0U;
        fullLengthByWhy[why] += static_cast<std::size_t>(n) >= static_cast<std::size_t>(Detector<Det>::nLayers) ? 1U : 0U;
      }
      if (beyond == 0U) {
        bool barrel = false;
        bool disc = false;
        for (int index = 0; index < n; ++index) {
          (Detector<Det>::isBarrelLayer(static_cast<LayerID>(layerOf(index))) ? barrel : disc) = true;
        }
        std::size_t kind = 0U;
        if (barrel && disc) {
          kind = 2U;
        } else if (disc) {
          kind = 1U;
        }
        ++diffuseByKind[kind];
      }
      const auto worstShare = static_cast<double>(outwardSpatial[static_cast<std::size_t>(worst)]) / std::max(static_cast<double>(fit.outwardChi2), 1e-30);
      ++failedByWorstShare[std::min<std::size_t>(static_cast<std::size_t>(worstShare * 4.), 3U)];
    }
  }

  ByLength all;
  std::string lengths;
  for (std::size_t n = 0; n < byLength.size(); ++n) {
    if (byLength[n].chains != 0U) {
      all.chains += byLength[n].chains;
      all.failed += byLength[n].failed;
      all.fitFailed += byLength[n].fitFailed;
      lengths += fmt::format(" {}: {} {:.2f}%", n, byLength[n].chains, percent(byLength[n].failed, byLength[n].chains));
    }
  }
  LOGP(info,
       "true roads, iteration 0 (pT_min {} GeV, {} clusters or more): {} chains; the road cut fails {} ({:.2f}%), a calibrated chi2 {:.2f}%; fit failed {}; "
       "a kink taken by {} ({:.2f}%, at most {})",
       settings.windows.minPtGeV, minClusters, all.chains, all.failed, percent(all.failed, all.chains), 100. * std::erfc(nSigma / std::sqrt(2.)), all.fitFailed,
       kinked, percent(kinked, all.chains), static_cast<int>(parameters.maxKinks));
  LOGP(info, "  by clusters, chains and failing:{}", lengths);
  reportTrueRoadFitFailures(frame, failures, particles, tracking, propagator, crossings, scattering, material);
  if (!chi2PerNdf.empty()) {
    std::ranges::sort(chi2PerNdf);
    const auto quantile = [&](double q) {
      return static_cast<double>(chi2PerNdf[std::min(chi2PerNdf.size() - 1U, static_cast<std::size_t>(q * static_cast<double>(chi2PerNdf.size())))]);
    };
    LOGP(info, "  outward chi2/ndf quantiles: 50% {:.2f}, 90% {:.2f}, 99% {:.2f}, 99.9% {:.2f}, max {:.1f}", quantile(0.5), quantile(0.9), quantile(0.99),
         quantile(0.999), static_cast<double>(chi2PerNdf.back()));
  }
  LOGP(info,
       "  failing chains by outward steps beyond the kink gate: 0 {}, 1 {}, 2 {}, 3+ {}; by the worst step's share of the chi2: <25% {}, 25-50% {}, "
       "50-75% {}, >75% {}",
       failedByBeyondGate[0], failedByBeyondGate[1], failedByBeyondGate[2], failedByBeyondGate[3], failedByWorstShare[0], failedByWorstShare[1],
       failedByWorstShare[2], failedByWorstShare[3]);
  {
    std::string whys;
    for (std::size_t why = 0; why < kWhys; ++why) {
      whys += fmt::format("{} {} {} (all-layer {}; the inward pass passes {})", why == 0U ? "" : ";", kWhyNames[why], failedByWhy[why], fullLengthByWhy[why],
                          inwardPassesByWhy[why]);
    }
    LOGP(info, "  failing chains by why the kink did not save them:{}", whys);
  }
  LOGP(info, "  failing chains without a step beyond the gate: barrel only {}, discs only {}, both {}", diffuseByKind[0], diffuseByKind[1], diffuseByKind[2]);
  LOGP(info, "  chi2/ndf {:.3f}", chi2Sum / std::max(ndfSum, 1.));
  for (std::size_t layer = 0; layer < spatial.size(); ++layer) {
    if (spatial[layer].empty()) {
      continue;
    }
    const auto beyondTwo = static_cast<std::uint64_t>(std::ranges::count_if(spatial[layer], [&](float v) { return static_cast<double>(v) > tailTwo; }));
    const auto count = spatial[layer].size();
    std::vector<double> updates(spatial[layer].begin(), spatial[layer].end()); // exact: the float values' median
    LOGP(info, "  layer {:>2} ({}): {:>7} updates, scale^2 {:.2f}, beyond {:.2f}%; worst of the failing {}", layer,
         Detector<Det>::isBarrelLayer(static_cast<LayerID>(layer)) ? "barrel" : "disc", count, median(updates) / kMedianChi2Two,
         percent(beyondTwo, count), worstOfFailed[layer]);
  }
}

/// Where the fakes go wrong, under the comparison's rule (TrackComparison.h:
/// the label is the particle most clusters carry, one cluster without it
/// makes a fake), by the track's layers -- every layer, the innermost four,
/// the outermost four, other -- and, for the fakes with one wrong cluster, by
/// what the label's own cluster on that layer could have made of it: none
/// there (the particle ended before the layer, or skipped it), or there, and
/// the track with it in place of the wrong one refitted as the road fit:
/// whether it passes the road cut, and whether it would have outranked the
/// fake in the selection (fewer kinks, then the lower outward chi2). A
/// passing alternative that outranks the fake never reached the selection,
/// so the label's stage record says where its chain was lost (the fake's
/// iteration's, in the tracker's own run: \a funnel; none for a label not
/// findable by the funnel's rule); one
/// that ranks below it was beaten by the wrong cluster on chi2. Also the
/// refit's chi2 per degree of freedom of the good tracks against the fakes',
/// for what a cut on it would reach.
void reportFakeAnatomy(const Frame& frame, const RecoSettings& param, std::span<const TruthRow> rows, const ParticleTable& particles, const Funnel* funnel,
                       const Controller& tracker, const TrackRefit<Det>::Tracks& tracks, const TrackingGeometry<Det>& tracking,
                       const o2::base::MatLayerCylSet& material, const Propagator& propagator, const char* label)
{
  constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  constexpr std::size_t kClasses = 4U;
  constexpr std::array<const char*, kClasses> kClassNames{"every layer", "innermost four", "outermost four", "other layers"};
  const auto nLayers = static_cast<std::uint32_t>(Detector<Det>::nLayers);
  const std::uint32_t innerFour = nLayers >= 4U ? 0xFU : 0U;
  std::uint32_t outerFour = 0U;
  for (std::uint32_t layer = nLayers >= 4U ? nLayers - 4U : 0U; layer < std::min(nLayers, 32U); ++layer) {
    outerFour |= 1U << layer;
  }
  const auto classOf = [&](std::uint32_t pattern) -> std::size_t {
    if (std::cmp_equal(std::popcount(pattern), nLayers)) {
      return 0U;
    }
    if (pattern == innerFour) {
      return 1U;
    }
    if (pattern == outerFour) {
      return 2U;
    }
    return 3U;
  };
  const auto primary = [&](std::uint64_t id) {
    const auto found = particles.find(id);
    return found != particles.end() && found->second.primary;
  };

  // Every track: its label and class; the fakes kept for the second look.
  struct Fake {
    std::size_t track = 0U;
    std::uint64_t label = kNoTrackParticle;
    std::uint32_t most = 0U;
    std::size_t cls = 0U;
  };
  std::vector<Fake> fakes;
  std::array<std::array<std::vector<float>, 2>, kClasses> chi2PerNdf; // [class][good, fake]
  std::array<std::array<std::uint64_t, 2>, kClasses> tracksByClass{};
  std::array<std::array<std::uint64_t, 2>, kClasses> kinkedByClass{};
  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> ownRows; // the fakes' labels -> every row carrying them
  const auto view = tracks.view();
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    const auto count = view.count[track];
    const std::span<const std::uint32_t> trackRows{view.clusters[track].index, count};
    const auto [id, most] = majorityParticle(rows, trackRows);
    std::uint32_t pattern = 0U;
    for (const auto row : trackRows) {
      pattern |= 1U << rows[row].layer;
    }
    const auto cls = classOf(pattern);
    const bool good = id != kNoTrackParticle && most == count;
    const std::size_t which = good ? 0U : 1U;
    ++tracksByClass[cls][which];
    kinkedByClass[cls][which] += view.kinks[track] != 0U ? 1U : 0U;
    if (view.ndf[track] > 0) {
      chi2PerNdf[cls][which].push_back(view.chi2[track] / static_cast<float>(view.ndf[track]));
    }
    if (!good && id != kNoTrackParticle) {
      fakes.push_back({.track = track, .label = id, .most = most, .cls = cls});
      ownRows.try_emplace(id);
    }
  }
  for (Frame::Index row = 0U; row < frame.size(); ++row) {
    const auto& truth = rows[row];
    for (std::uint8_t p = 0; p < truth.count; ++p) {
      const auto found = ownRows.find(truth.particles[p]);
      if (found != ownRows.end()) {
        found->second.push_back(row);
      }
    }
  }

  // The fakes with one wrong cluster, and the alternative with the label's own.
  struct ClassTally {
    std::uint64_t oneWrong = 0U;
    std::array<std::uint64_t, 3> wrongAt{}; ///< innermost, outermost, inside
    std::uint64_t wrongSecondary = 0U;      ///< the wrong cluster carries no primary
    std::uint64_t holeEnded = 0U;           ///< the label has no cluster on the layer, nor beyond it
    std::uint64_t holeSkipped = 0U;         ///< none on the layer, some beyond
    std::uint64_t altFitFailed = 0U;
    std::uint64_t altFailsCut = 0U;
    std::uint64_t altAbove = 0U;                ///< passes and outranks the fake
    std::uint64_t altBelow = 0U;                ///< passes and ranks below it
    std::array<std::uint64_t, 7> aboveLostAt{}; ///< tracklets, triplets, neighbours, roads, road cut, not lost, no record
    std::array<std::uint64_t, 3> noRecord{};    ///< of the last: a primary below the iteration's pT floor, a secondary, other
    std::vector<float> belowGap;                ///< ranks below: the alternative's outward chi2 less the fake's
  };
  std::array<ClassTally, kClasses> tally{};
  std::array<std::uint32_t, kMaxClusters> trackRows{};
  std::array<std::uint32_t, kMaxClusters> altRows{};
  for (const auto& fake : fakes) {
    const auto count = view.count[fake.track];
    if (count - fake.most != 1U) {
      continue;
    }
    auto& t = tally[fake.cls];
    ++t.oneWrong;
    std::uint32_t wrong = 0U;
    for (std::uint32_t i = 0; i < count; ++i) {
      trackRows[i] = view.clusters[fake.track].index[i];
      if (!rows[trackRows[i]].carries(fake.label)) {
        wrong = i;
      }
    }
    std::size_t at = 2U;
    if (wrong == 0U) {
      at = 0U;
    } else if (wrong == count - 1U) {
      at = 1U;
    }
    ++t.wrongAt[at];
    const auto& wrongTruth = rows[trackRows[wrong]];
    t.wrongSecondary += std::none_of(wrongTruth.particles.begin(), wrongTruth.particles.begin() + wrongTruth.count, primary) ? 1U : 0U;
    const auto layer = wrongTruth.layer;
    const auto& own = ownRows[fake.label];
    std::uint32_t labelLayers = 0U;
    for (const auto row : own) {
      labelLayers |= 1U << rows[row].layer;
    }
    if ((labelLayers & (1U << layer)) == 0U) {
      ++((labelLayers >> layer) == 0U ? t.holeEnded : t.holeSkipped);
      continue;
    }
    const auto& parameters = tracker.iteration(view.iteration[fake.track]);
    const int n = static_cast<int>(count);
    const auto fitRows = [&](const std::uint32_t* chain) {
      return fitMeasurements(propagator, tracking, material, RowChain<Det>{.frame = &frame, .rows = chain}, n, kinkBudget(parameters));
    };
    const auto fakeFit = fitRows(trackRows.data());
    const auto passes = [&](const TrackFitResult& fit) {
      return fit.status == PropagationStatus::kOk && fit.outwardChi2 <= roadChi2Cut(parameters, count, fit.outwardKinks);
    };
    // The best alternative: a passing one with the lowest chi2, else the lowest.
    std::optional<TrackFitResult> best;
    for (const auto row : own) {
      if (rows[row].layer != layer) {
        continue;
      }
      altRows = trackRows;
      altRows[wrong] = row;
      const auto fit = fitRows(altRows.data());
      if (fit.status != PropagationStatus::kOk) {
        if (!best) {
          best = fit;
        }
        continue;
      }
      const bool better = !best || best->status != PropagationStatus::kOk || (passes(fit) && !passes(*best)) ||
                          (passes(fit) == passes(*best) && fit.outwardChi2 < best->outwardChi2);
      if (better) {
        best = fit;
      }
    }
    if (!best || best->status != PropagationStatus::kOk) {
      ++t.altFitFailed;
      continue;
    }
    if (!passes(*best)) {
      ++t.altFailsCut;
      continue;
    }
    const bool above = fakeFit.status != PropagationStatus::kOk || best->outwardKinks < fakeFit.outwardKinks ||
                       (best->outwardKinks == fakeFit.outwardKinks && best->outwardChi2 < fakeFit.outwardChi2);
    if (!above) {
      ++t.altBelow;
      if (fakeFit.status == PropagationStatus::kOk) {
        t.belowGap.push_back(best->outwardChi2 - fakeFit.outwardChi2);
      }
      continue;
    }
    ++t.altAbove;
    std::size_t lostAt = 6U;
    if (const auto stage = funnel != nullptr ? funnel->stages(fake.label, view.iteration[fake.track]) : std::nullopt; stage.has_value()) {
      lostAt = static_cast<std::size_t>(firstLoss(stage->longest, count, false));
    }
    ++t.aboveLostAt[lostAt];
    if (lostAt == 6U) {
      const auto found = particles.find(fake.label);
      std::size_t kind = 2U;
      if (found != particles.end() && !found->second.primary) {
        kind = 1U;
      } else if (found != particles.end() &&
                 static_cast<double>(found->second.pt) < workflow::iterationSettings(param, static_cast<int>(view.iteration[fake.track])).windows.minPtGeV) {
        kind = 0U;
      }
      ++t.noRecord[kind];
    }
  }

  const auto quantile = [](std::vector<float>& values, double q) {
    if (values.empty()) {
      return 0.;
    }
    const auto at = std::min(values.size() - 1U, static_cast<std::size_t>(q * static_cast<double>(values.size())));
    std::ranges::nth_element(values, values.begin() + static_cast<std::ptrdiff_t>(at));
    return static_cast<double>(values[at]);
  };
  LOGP(info, "fake anatomy, {}: the comparison's rule, by the track's layers; the one-wrong-cluster fakes by the label's own cluster on that layer", label);
  for (std::size_t cls = 0; cls < kClasses; ++cls) {
    const auto good = tracksByClass[cls][0];
    const auto bad = tracksByClass[cls][1];
    if (good + bad == 0U) {
      continue;
    }
    auto& goodChi2 = chi2PerNdf[cls][0];
    auto& fakeChi2 = chi2PerNdf[cls][1];
    const double good90 = quantile(goodChi2, 0.9);
    const double good99 = quantile(goodChi2, 0.99);
    const auto fakesAbove = [&](double bound) {
      return percent(static_cast<std::uint64_t>(std::ranges::count_if(fakeChi2, [&](float v) { return static_cast<double>(v) > bound; })), fakeChi2.size());
    };
    LOGP(info,
         "  {:<14} {:>7} tracks, {:>6} fake ({:.2f} %); kinked {:.1f} % of the good, {:.1f} % of the fakes; refit chi2/ndf 50/90/99 %: good "
         "{:.2f}/{:.2f}/{:.2f}, fakes {:.2f}/{:.2f}/{:.2f}; fakes above the good's 90 % {:.1f} %, above its 99 % {:.1f} %",
         kClassNames[cls], good + bad, bad, percent(bad, good + bad), percent(kinkedByClass[cls][0], good), percent(kinkedByClass[cls][1], bad),
         quantile(goodChi2, 0.5), good90, good99, quantile(fakeChi2, 0.5), quantile(fakeChi2, 0.9), quantile(fakeChi2, 0.99), fakesAbove(good90),
         fakesAbove(good99));
    const auto& t = tally[cls];
    if (t.oneWrong == 0U) {
      continue;
    }
    const auto swaps = t.altFitFailed + t.altFailsCut + t.altAbove + t.altBelow;
    auto& gap = tally[cls].belowGap;
    LOGP(info,
         "    one wrong cluster {} ({:.1f} % of the fakes): innermost {}, outermost {}, inside {}; it carries no primary {:.1f} %; the label has none on "
         "that layer {} (ended before it {}, skipped it {}), its own there {}: with it the track fails to fit {}, fails the road cut {}, passes and "
         "outranks the fake {} (its chain lost at tracklets {}, triplets {}, neighbours {}, roads {}, the road cut {}; not lost {}; no record {}: a "
         "primary below the pT floor {}, a secondary {}, other {}), passes but ranks below the fake {} (the chi2 gap, 10/50/90 %: {:.2f}/{:.2f}/{:.2f})",
         t.oneWrong, percent(t.oneWrong, bad), t.wrongAt[0], t.wrongAt[1], t.wrongAt[2], percent(t.wrongSecondary, t.oneWrong), t.holeEnded + t.holeSkipped,
         t.holeEnded, t.holeSkipped, swaps, t.altFitFailed, t.altFailsCut, t.altAbove, t.aboveLostAt[0], t.aboveLostAt[1], t.aboveLostAt[2], t.aboveLostAt[3],
         t.aboveLostAt[4], t.aboveLostAt[5], t.aboveLostAt[6], t.noRecord[0], t.noRecord[1], t.noRecord[2], t.altBelow, quantile(gap, 0.1), quantile(gap, 0.5),
         quantile(gap, 0.9));
  }
}

/// Every track's time box, its rows' intersection: none may be empty; narrower than a row's box only where boxes are
/// widened or staggered. And the frames a track's rows come from, told apart
/// by their boxes' starts: more than one only with a margin.
void reportTrackTimeBoxes(const Frame& frame, const typename TrackRefit<Det>::Tracks& tracks)
{
  const auto view = tracks.view();
  std::uint64_t empty = 0U;
  std::uint64_t narrowed = 0U;
  std::array<std::uint64_t, 3> byFrames{}; // 1, 2, 3 or more
  double widthSumNS = 0.;
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    const auto& clusters = view.clusters[track];
    const auto count = view.count[track];
    const auto box = trackTimeBox(frame.measurements(), clusters, count);
    empty += box.width == 0U ? 1U : 0U;
    narrowed += box.width < timeBoxAt(frame.measurements(), clusters.index[0]).width ? 1U : 0U;
    widthSumNS += 1.e-3 * static_cast<double>(box.width);
    std::array<std::uint32_t, std::extent_v<std::remove_cvref_t<decltype(clusters.index)>>> starts{};
    const auto rows = std::min<std::size_t>(count, starts.size());
    for (std::size_t i = 0; i < rows; ++i) {
      starts[i] = timeBoxAt(frame.measurements(), clusters.index[i]).start;
    }
    const auto used = std::span{starts}.first(rows);
    std::ranges::sort(used);
    const auto frames = static_cast<std::size_t>(std::ranges::distance(used.begin(), std::ranges::unique(used).begin()));
    ++byFrames[std::min<std::size_t>(frames, 3U) - 1U];
  }
  LOGP(info, "track time boxes: {} tracks, {} empty, {} narrower than their first row's, mean width {:.1f} ns; by frames: 1: {} 2: {} 3+: {}", tracks.size(),
       empty, narrowed, tracks.size() > 0U ? widthSumNS / static_cast<double>(tracks.size()) : 0., byFrames[0], byFrames[1], byFrames[2]);
}

/// The pT edges [GeV] of measureTracks' and reportIncomplete's tables.
constexpr std::array<Scalar, 11> kTrackPtEdges{0.F, 0.5F, 0.75F, 1.F, 1.25F, 1.5F, 2.F, 3.F, 5.F, 10.F, std::numeric_limits<Scalar>::infinity()};
constexpr std::size_t kTrackPtBins = kTrackPtEdges.size() - 1U;

/// A particle lost at the tracklets: its lost true pairs of \a survived into
/// \a pairs by the cut, or one more in \a none if it has none recorded.
void addLostPairs(const Funnel::IterationStages& survived, std::array<std::uint64_t, kPairLosses>& pairs, std::uint64_t& none)
{
  bool any = false;
  for (std::size_t reason = 0; reason < kPairLosses; ++reason) {
    pairs[reason] += survived.lostPairs[reason];
    any |= survived.lostPairs[reason] != 0U;
  }
  none += any ? 0U : 1U;
}

/// Found incomplete: a particle whose true
/// chain spans every layer, found with a good track short of one. By the
/// iteration of its longest good track; whether the particle's pT is below
/// that iteration's floor or its floor for a full-length road (the stub
/// passed, the full road could not); which layers the track lacks; whether
/// a lacking cluster is held by another track, good or fake, or is free;
/// and, if every lacking cluster is free, the first stage where the
/// particle's chain fell short in that iteration of the tracker's own run (\a
/// funnel; none without one). \a particles:
/// measureTracks's, by id, with `pt` and `longest`; \a goodTracks: each
/// particle's good tracks.
template <typename Particles>
void reportIncomplete(const Frame& frame, const RecoSettings& param, double bzKiloGauss, const TrackingGeometry<Det>& tracking, const typename TrackRefit<Det>::Tracks& tracks,
                      std::span<const TruthRow> truthRows, const Funnel* funnel, const Particles& particles,
                      const std::unordered_map<std::uint64_t, std::vector<std::size_t>>& goodTracks, std::size_t iterations, const char* set,
                      const TrueCrossings& crossings = {}, const ScatteringTable<Det>* scattering = nullptr, const o2::base::MatLayerCylSet* material = nullptr)
{
  const auto view = tracks.view();
  constexpr auto nLayers = static_cast<std::uint32_t>(Detector<Det>::nLayers);
  std::vector<std::uint32_t> holder(frame.size(), std::numeric_limits<std::uint32_t>::max());
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    for (std::uint32_t i = 0; i < view.count[track]; ++i) {
      holder[view.clusters[track].index[i]] = static_cast<std::uint32_t>(track);
    }
  }
  std::vector<bool> isGood(tracks.size(), false);
  for (const auto& [id, list] : goodTracks) {
    for (const auto track : list) {
      isGood[track] = true;
    }
  }
  std::unordered_map<std::uint64_t, std::vector<Frame::Index>> rowsOf; // the all-layer particles' rows
  for (const auto& [id, particle] : particles) {
    if (particle.longest >= nLayers && goodTracks.contains(id)) {
      rowsOf[id];
    }
  }
  for (Frame::Index row = 0U; row < frame.size(); ++row) {
    for (std::uint8_t i = 0; i < truthRows[row].count; ++i) {
      if (const auto found = rowsOf.find(truthRows[row].particles[i]); found != rowsOf.end()) {
        found->second.push_back(row);
      }
    }
  }
  const auto layerOfRow = [&](Frame::Index row) { return static_cast<std::uint32_t>(tracking.surface(frame.measurement(row).surface).layer); };
  enum Why : std::uint8_t {
    kWhyBelowFloor,
    kWhyTaken,
    kWhyFreeTracklets,
    kWhyFreeTriplets,
    kWhyFreeNeighbours,
    kWhyFreeRoads,
    kWhyFreeFit,
    kWhyFreeSelection,
    kWhyFreeLater,
    kWhys
  };
  constexpr std::array<const char*, kWhys> kWhyNames{"below the iteration's floors", "a lacking cluster taken by another track",
                                                     "free, chain lost at tracklets", "free, at triplets",
                                                     "free, at neighbours", "free, at roads",
                                                     "free, at the road fit", "free, in the selection",
                                                     "free, not in the funnel"};
  std::array<std::uint64_t, kWhys> whys{};
  std::array<std::array<std::uint64_t, 3>, kWhys> whysByPt{}; // below 0.3 GeV, 0.3-1, above 1
  constexpr auto kLayerSlots = static_cast<std::size_t>(Detector<Det>::nLayers);
  constexpr std::size_t kNotALink = nLinks<Det>;                               // the link that broke is none of the table's
  std::array<std::array<std::uint64_t, nLinks<Det> + 1U>, kWhys> brokenLink{}; // by why, the link that broke
  std::array<std::uint64_t, kPairLosses> trackletLostPairs{};
  std::uint64_t trackletLostNone = 0U;
  std::vector<std::uint64_t> byIteration(iterations, 0U);
  std::uint64_t allLayerFound = 0U;
  std::uint64_t incomplete = 0U;
  std::uint64_t lacksOuter = 0U;  // only layers beyond the track's outermost
  std::uint64_t lacksInner = 0U;  // only layers inside its innermost
  std::uint64_t lacksMiddle = 0U; // a hole, or both ends
  std::uint64_t takenByGood = 0U;
  std::uint64_t takenByFake = 0U;
  std::array<std::uint64_t, kTrackPtBins> incompleteByPt{};
  std::array<std::uint64_t, kTrackPtBins> allLayerByPt{};
  // The MC hits along every all-layer particle found: the kink at each
  // crossing over its layer's theta0, as a check of the scattering model
  // itself, and the clusters off their crossing; for the incomplete, the
  // hardest kink, whether at the link that broke, by why.
  const bool withHits = static_cast<bool>(crossings) && scattering != nullptr;
  std::uint64_t particlesWithHits = 0U;
  std::array<std::array<std::uint64_t, kKinkEdges.size()>, kLayerSlots> kinkByScatterer{};
  // The same by |z| of the crossing the kink shows at: the material between
  // the barrels depends on z (the services), the layer's band average does not.
  constexpr std::array<double, 4> kZEdges{0., 10., 20., std::numeric_limits<double>::infinity()};
  const auto zBinOf = [&](double z) {
    std::size_t zBin = 0U;
    while (zBin + 2U < kZEdges.size() && std::abs(z) >= kZEdges[zBin + 1U]) {
      ++zBin;
    }
    return zBin;
  };
  std::array<std::array<std::array<std::uint64_t, kKinkEdges.size()>, kZEdges.size() - 1U>, kLayerSlots> kinkByScattererZ{};
  // Over the LUT's own theta0 along the path, and that path's mean x/X0.
  std::array<std::array<std::uint64_t, kKinkEdges.size()>, kLayerSlots> kinkPathByScatterer{};
  std::array<double, kLayerSlots> pathX2X0Sum{};
  std::array<std::uint64_t, kLayerSlots> pathCount{};
  std::array<std::array<double, kZEdges.size() - 1U>, kLayerSlots> pathNormalSumZ{}; // normal-incidence x/X0 by |z| of the outer crossing
  std::array<std::array<std::uint64_t, kZEdges.size() - 1U>, kLayerSlots> pathCountZ{};
  std::array<std::uint64_t, kLayerSlots> offByLayer{};
  std::array<std::array<std::uint64_t, kKinkEdges.size()>, kWhys> incompleteKink{};
  std::array<std::array<std::uint64_t, kLayerSlots>, kWhys> incompleteKinkLayer{}; // the hardest kink's scatterer, for the hard ones
  std::array<std::uint64_t, kWhys> incompleteKinkAtBroken{};
  std::array<std::uint64_t, kWhys> incompleteOff{};
  std::array<std::uint64_t, kWhys> incompleteNoHits{};
  for (const auto& [id, rows] : rowsOf) {
    const auto& list = goodTracks.at(id);
    const auto best = *std::ranges::max_element(list, {}, [&](std::size_t track) { return view.count[track]; });
    const auto& particle = particles.at(id);
    const auto ptBin = binOf(particle.pt, kTrackPtEdges);
    ++allLayerFound;
    ++allLayerByPt[ptBin];
    std::vector<TrueCrossingReport> reports;
    double hardest = 0.;
    std::uint32_t hardestAt = 0U;
    bool anyOff = false;
    if (withHits) {
      reports = trueAnatomy(frame, *scattering, bzKiloGauss, crossings(id, rows), material);
      particlesWithHits += reports.empty() ? 0U : 1U;
      for (const auto& report : reports) {
        if (std::isfinite(report.kinkPath)) {
          const auto slot = std::min<std::size_t>(report.scatterer, kLayerSlots - 1U);
          ++kinkPathByScatterer[slot][kinkBin(report.kinkPath)];
          pathX2X0Sum[slot] += report.pathX2X0;
          ++pathCount[slot];
          const auto zBin = zBinOf(report.z);
          pathNormalSumZ[slot][zBin] += report.pathNormalX2X0;
          ++pathCountZ[slot][zBin];
        }
        // The hardest kink by the theta0 of the material actually crossed, with a LUT.
        if (const double kink = std::isfinite(report.kinkPath) ? report.kinkPath : report.kink; std::isfinite(kink) && kink > hardest) {
          hardest = kink;
          hardestAt = report.scatterer;
        }
        if (std::isfinite(report.kink)) {
          ++kinkByScatterer[std::min<std::size_t>(report.scatterer, kLayerSlots - 1U)][kinkBin(report.kink)];
          ++kinkByScattererZ[std::min<std::size_t>(report.scatterer, kLayerSlots - 1U)][zBinOf(report.z)][kinkBin(report.kink)];
        }
        if (report.offsetCM > kOffsetCM) {
          ++offByLayer[std::min<std::size_t>(report.layer, kLayerSlots - 1U)];
          anyOff = true;
        }
      }
    }
    if (view.count[best] >= nLayers) {
      continue;
    }
    ++incomplete;
    ++incompleteByPt[ptBin];
    const auto iteration = std::min<std::size_t>(view.iteration[best], iterations - 1U);
    ++byIteration[iteration];
    std::uint32_t has = 0U;
    for (std::uint32_t i = 0; i < view.count[best]; ++i) {
      has |= 1U << layerOfRow(view.clusters[best].index[i]);
    }
    const auto innermost = static_cast<std::uint32_t>(std::countr_zero(has));
    const auto outermost = static_cast<std::uint32_t>(std::bit_width(has)) - 1U;
    bool outer = true;
    bool inner = true;
    bool anyTaken = false;
    for (const auto row : rows) {
      const auto layer = layerOfRow(row);
      if ((has & (1U << layer)) != 0U) {
        continue;
      }
      outer &= layer > outermost;
      inner &= layer < innermost;
      const auto other = holder[row];
      if (other != std::numeric_limits<std::uint32_t>::max() && other != best) {
        anyTaken = true;
        ++(isGood[other] ? takenByGood : takenByFake);
      }
    }
    if (outer) {
      ++lacksOuter;
    } else if (inner) {
      ++lacksInner;
    } else {
      ++lacksMiddle;
    }
    const auto settings = workflow::iterationSettings(param, static_cast<int>(iteration));
    const double fullFloor = nLayers < settings.minPtGeVByClusters.size() ? settings.minPtGeVByClusters[nLayers] : 0.;
    const auto survived = funnel != nullptr ? funnel->stages(id, iteration) : std::nullopt;
    Why why = kWhyFreeLater;
    if (static_cast<double>(particle.pt) < std::max(settings.windows.minPtGeV, fullFloor)) {
      why = kWhyBelowFloor;
    } else if (anyTaken) {
      why = kWhyTaken;
    } else if (survived.has_value()) {
      switch (firstLoss(survived->longest, nLayers, false)) {
        case FunnelStage::kTracklets:
          why = kWhyFreeTracklets;
          addLostPairs(*survived, trackletLostPairs, trackletLostNone);
          break;
        case FunnelStage::kTriplets:
          why = kWhyFreeTriplets;
          break;
        case FunnelStage::kNeighbours:
          why = kWhyFreeNeighbours;
          break;
        case FunnelStage::kRoads:
          why = kWhyFreeRoads;
          break;
        case FunnelStage::kFit:
          why = kWhyFreeFit;
          break;
        default:
          why = kWhyFreeSelection;
          break;
      }
    }
    ++whys[why];
    std::size_t ptSlot = 2U;
    if (particle.pt < 0.3F) {
      ptSlot = 0U;
    } else if (particle.pt < 1.F) {
      ptSlot = 1U;
    }
    ++whysByPt[why][ptSlot];
    // The link that broke: into the first lacking layer beyond the track, or
    // out of the last lacking layer inside it -- from its inner layer to the
    // next (the particle has a row on every layer), as the table's link.
    std::uint32_t brokenInner = outermost;
    if (!outer) {
      brokenInner = 0U;
      for (const auto row : rows) {
        const auto layer = layerOfRow(row);
        if ((has & (1U << layer)) == 0U && layer < innermost) {
          brokenInner = std::max(brokenInner, layer);
        }
      }
    }
    const auto broken = std::ranges::find_if(Det::layerLinks, [&](const auto& link) {
      return std::cmp_equal(link.from, brokenInner) && std::cmp_equal(link.to, brokenInner + 1U);
    });
    ++brokenLink[why][static_cast<std::size_t>(broken - Det::layerLinks.begin())];
    if (withHits) {
      if (reports.empty()) {
        ++incompleteNoHits[why];
      } else {
        ++incompleteKink[why][kinkBin(hardest)];
        if (hardest >= kHardKink) {
          ++incompleteKinkLayer[why][std::min<std::size_t>(hardestAt, kLayerSlots - 1U)];
          incompleteKinkAtBroken[why] += hardestAt == brokenInner ? 1U : 0U;
        }
        incompleteOff[why] += anyOff ? 1U : 0U;
      }
    }
  }
  std::string whyLine;
  for (std::size_t why = 0; why < kWhys; ++why) {
    whyLine += fmt::format("{} {} {} ({}/{}/{} below 0.3, 0.3-1, above 1 GeV)", why == 0U ? "" : ",", kWhyNames[why], whys[why], whysByPt[why][0],
                           whysByPt[why][1], whysByPt[why][2]);
  }
  std::string iterationLine;
  for (std::size_t iteration = 0; iteration < iterations; ++iteration) {
    iterationLine += fmt::format(" {}: {}", iteration, byIteration[iteration]);
  }
  std::string ptLine;
  for (std::size_t bin = 0; bin < kTrackPtBins; ++bin) {
    if (allLayerByPt[bin] != 0U) {
      ptLine += fmt::format(" [{:.2f},{:.2f}) {:.1f}", kTrackPtEdges[bin], kTrackPtEdges[bin + 1U], percent(incompleteByPt[bin], allLayerByPt[bin]));
    }
  }
  LOGP(info,
       "  found incomplete, {}: {} of {} all-layer particles found ({:.2f}%); by the iteration of the longest good track:{}; lacking layers beyond it "
       "{}, inside it {}, a hole or both ends {}; lacking clusters held by a good track {}, by a fake {}",
       set, incomplete, allLayerFound, percent(incomplete, allLayerFound), iterationLine, lacksOuter, lacksInner, lacksMiddle, takenByGood, takenByFake);
  LOGP(info, "    why:{}", whyLine);
  std::string cuts;
  for (std::size_t reason = 0; reason < kPairLosses; ++reason) {
    cuts += fmt::format(" {} {}", kPairLossNames[reason], trackletLostPairs[reason]);
  }
  LOGP(info, "    the tracklet losses' adjacent-layer pairs of their own their iteration did not keep, by the first test each fails:{}; none recorded {}", cuts,
       trackletLostNone);
  for (std::size_t why = 0; why < kWhys; ++why) {
    if (whys[why] == 0U) {
      continue;
    }
    std::string links;
    for (std::size_t link = 0; link < nLinks<Det>; ++link) {
      if (brokenLink[why][link] != 0U) {
        links += fmt::format(" {}->{}: {}", Det::layerLinks[link].from, Det::layerLinks[link].to, brokenLink[why][link]);
      }
    }
    if (brokenLink[why][kNotALink] != 0U) {
      links += fmt::format(" no such link: {}", brokenLink[why][kNotALink]);
    }
    LOGP(info, "    {}, by the link that broke:{}", kWhyNames[why], links);
  }
  LOGP(info, "    incomplete % by pT [GeV]:{}", ptLine);
  if (withHits) {
    std::string spectrum;
    std::string off;
    for (std::size_t layer = 0; layer < kLayerSlots; ++layer) {
      if (std::accumulate(kinkByScatterer[layer].begin(), kinkByScatterer[layer].end(), std::uint64_t{0}) != 0U) {
        spectrum += fmt::format(" {}:{};", layer, kinkBinsLine(kinkByScatterer[layer]));
      }
      if (offByLayer[layer] != 0U) {
        off += fmt::format(" {}: {}", layer, offByLayer[layer]);
      }
    }
    LOGP(info,
         "    the MC hits on the {} all-layer particles found with hits -- the kink over the layer's theta0, by the layer that scattered:{} clusters off "
         "their crossing by more than {:.0f} um, by layer:{}",
         particlesWithHits, spectrum, 1e4 * kOffsetCM, off);
    std::string byZ;
    for (std::size_t layer = 0; layer < kLayerSlots; ++layer) {
      std::string bins;
      for (std::size_t zBin = 0; zBin + 1U < kZEdges.size(); ++zBin) {
        const auto& counts = kinkByScattererZ[layer][zBin];
        const auto all = std::accumulate(counts.begin(), counts.end(), std::uint64_t{0});
        if (all == 0U) {
          continue;
        }
        std::uint64_t tail = 0U;
        std::uint64_t hard = 0U;
        for (std::size_t bin = 0; bin < kKinkEdges.size(); ++bin) {
          tail += kKinkEdges[bin] >= kTailKink ? counts[bin] : 0U;
          hard += kKinkEdges[bin] >= kHardKink ? counts[bin] : 0U;
        }
        bins += fmt::format(" |z| {:.0f}-{}: {} ({:.1f} % >= 3, {:.2f} % >= 5)", kZEdges[zBin],
                            std::isfinite(kZEdges[zBin + 1U]) ? fmt::format("{:.0f}", kZEdges[zBin + 1U]) : std::string{"inf"}, all, percent(tail, all),
                            percent(hard, all));
      }
      if (!bins.empty()) {
        byZ += fmt::format(" layer {}:{};", layer, bins);
      }
    }
    LOGP(info, "    the kink over theta0 by the layer that scattered and |z| [cm] of the crossing it shows at:{}", byZ);
    if (material != nullptr) {
      std::string path;
      for (std::size_t layer = 0; layer < kLayerSlots; ++layer) {
        if (pathCount[layer] != 0U) {
          path += fmt::format(" {} (mean x/X0 {:.3f} %):{};", layer, 100. * pathX2X0Sum[layer] / static_cast<double>(pathCount[layer]),
                              kinkBinsLine(kinkPathByScatterer[layer]));
        }
      }
      LOGP(info, "    the same kink over theta0 from the LUT's x/X0 along the path between the two crossings, by the layer left:{}", path);
      std::string normal;
      for (std::size_t layer = 0; layer < kLayerSlots; ++layer) {
        std::string bins;
        for (std::size_t zBin = 0; zBin + 1U < kZEdges.size(); ++zBin) {
          if (pathCountZ[layer][zBin] != 0U) {
            bins += fmt::format(" |z| {:.0f}-{} {:.3f} %", kZEdges[zBin],
                                std::isfinite(kZEdges[zBin + 1U]) ? fmt::format("{:.0f}", kZEdges[zBin + 1U]) : std::string{"inf"},
                                100. * pathNormalSumZ[layer][zBin] / static_cast<double>(pathCountZ[layer][zBin]));
          }
        }
        if (!bins.empty()) {
          normal += fmt::format(" {}:{};", layer, bins);
        }
      }
      LOGP(info, "    the LUT's x/X0 along the path at normal incidence, by the layer left and |z| [cm] of the crossing reached:{}", normal);
    }
    for (std::size_t why = 0; why < kWhys; ++why) {
      if (whys[why] == 0U) {
        continue;
      }
      std::string layers;
      for (std::size_t layer = 0; layer < kLayerSlots; ++layer) {
        if (incompleteKinkLayer[why][layer] != 0U) {
          layers += fmt::format(" {}: {}", layer, incompleteKinkLayer[why][layer]);
        }
      }
      LOGP(info,
           "    {}: hardest kink over theta0{}; the kinks >= 5 by the layer that scattered:{}, at the link that broke {}; a cluster off its crossing {}; no "
           "hits {}",
           kWhyNames[why], kinkBinsLine(incompleteKink[why]), layers, incompleteKinkAtBroken[why], incompleteOff[why], incompleteNoHits[why]);
    }
  }
}

/// The lost true pairs by the cut, and the particles with none recorded.
[[nodiscard]] std::string describeLostPairs(const std::array<std::uint64_t, kPairLosses>& pairs, std::uint64_t none)
{
  std::string text;
  for (std::size_t reason = 0; reason < kPairLosses; ++reason) {
    if (pairs[reason] != 0U) {
      text += fmt::format(" {} {}", kPairLossNames[reason], pairs[reason]);
    }
  }
  return fmt::format("{}; particles with none {}", text.empty() ? " -" : text, none);
}

/// The findable particles by their origin's distance in z from the beam spot,
/// in its sigmaZ: the beam line's pair cut reaches nSigma of them.
class OriginZTally
{
 public:
  explicit OriginZTally(const VertexConstraint& beamSpot) : mBeamSpot(beamSpot) {}

  void addFindable(Scalar originZ, bool found)
  {
    const auto bin = binOfZ(originZ);
    ++mFindable[bin];
    mFound[bin] += found ? 1U : 0U;
  }

  void addTrackletLoss(Scalar originZ) { ++mTrackletLosses[binOfZ(originZ)]; }

  void report(const char* set) const
  {
    std::string line;
    for (std::size_t bin = 0; bin < kBins; ++bin) {
      line += fmt::format(" [{:g},{:g}) {}/{} ({:.2f}%), {} lost at the tracklets;", kEdges[bin], kEdges[bin + 1U], mFound[bin], mFindable[bin],
                          percent(mFound[bin], mFindable[bin]), mTrackletLosses[bin]);
    }
    LOGP(info, "  {}, found of the findable by the origin's |z - z_beam| in the beam spot's sigmaZ ({:.2f} cm):{}", set, mBeamSpot.sigmaZ, line);
  }

 private:
  static constexpr std::array<Scalar, 6> kEdges{0.F, 1.F, 2.F, 3.F, 4.F, std::numeric_limits<Scalar>::infinity()};
  static constexpr std::size_t kBins = kEdges.size() - 1U;

  [[nodiscard]] std::size_t binOfZ(Scalar originZ) const
  {
    return binOf(static_cast<Scalar>(std::abs(static_cast<double>(originZ) - mBeamSpot.z) / mBeamSpot.sigmaZ), kEdges);
  }

  VertexConstraint mBeamSpot;
  std::array<std::uint64_t, kBins> mFindable{};
  std::array<std::uint64_t, kBins> mFound{};
  std::array<std::uint64_t, kBins> mTrackletLosses{};
};

/// The tracker's run on the timeframe: its digests, the per-iteration line
/// and the track time boxes; with MC, the tracks against the truth. A particle is **findable** by the
/// comparison's rule (comparisonSettings, particleLayers, isFindable: its pT
/// and the layers its clusters cover), in two sets, the physical primaries
/// and the rest. A track is **good** if every cluster carries one particle,
/// a **fake** otherwise (ITS's convention); a findable particle is **found**
/// if a good track is its; each further good track of it is a **clone**.
/// \a tracked, the comparisons, gets the tracks of a complete run first. The
/// run is \a interface's; its timeframe ends with the caller's guard.
void measureTracks(const Frame& frame, const RecoSettings& param, const VertexConstraint& beamSpot, std::span<const VertexSeed> vertices,
                   std::span<const TruthRow> truthRows, const ParticleTable& truth, const Funnel* funnel, workflow::TrackingInterface<Det>& interface,
                   const std::function<void(const TrackRefit<Det>::Tracks&)>& tracked = {}, const TrueCrossings& crossings = {})
{
  const auto* label = interface.trackerLabel().c_str();
  const auto& tracking = interface.tracking();
  const auto& material = interface.material();
  const auto* scattering = &interface.scattering();
  // Findable by the comparison's rule (isFindable): the funnel's, the QA's.
  const auto first = workflow::iterationSettings(param, 0);
  const auto findableRule = comparisonSettings<Det>();
  const auto minLayers = static_cast<std::uint32_t>(findableRule.minLayers);
  const auto iterations = static_cast<std::size_t>(workflow::iterationCount(param));
  auto& tracker = interface.tracker();
  const auto run = interface.track(frame, vertices);
  if (run.dropped) {
    return;
  }
  const auto& result = run.result;
  const double seconds = run.seconds;
  if (tracked) {
    tracked(*result.tracks);
  }
  const auto& tracks = *result.tracks;
  const auto view = tracks.view();
  LOGP(info, "tracks, {}: {} in {:.3f} s", label, tracks.size(), seconds);
  // The output's fingerprint: what a speed change must leave alone.
  LOGP(info, "digest, {}: {:016x}", label, trackingDigest(result, tracker.used()));
  LOGP(info, "digest, {}, clusters: {:016x}", label, clusterDigest(result, frame, tracker.used()));
  using profile::Phase;
  for (std::size_t iteration = 0; iteration < result.iterations.size(); ++iteration) {
    const auto& summary = result.iterations[iteration];
    for (const auto& line : study::describeIteration(summary, iteration)) {
      LOGP(info, "{}", line);
    }
    // The serial passes between the parallel stages, inside the slice and road
    // times above, from the profile. The share is of those two.
    LOGP(info, "    slices by stage: {:.3f} s tracklets, {:.3f} s triplets, {:.3f} s neighbours", summary.seconds(Phase::kTracklets),
         summary.seconds(Phase::kTriplets), summary.seconds(Phase::kNeighbours));
    if (summary.roads == 0U) {
      continue; // an empty timeframe times nothing; the shares would be noise
    }
    const auto sum = [&](std::initializer_list<Phase> phases) {
      double total = 0.;
      for (const auto phase : phases) {
        total += summary.seconds(phase);
      }
      return total;
    };
    const double resets = sum({Phase::kTrackletReset, Phase::kTripletReset, Phase::kNeighbourReset, Phase::kRoadFitReset});
    const double scans = sum({Phase::kTrackletScan, Phase::kTripletScan, Phase::kNeighbourScan, Phase::kRoadFitScan});
    const double serial = resets + scans + sum({Phase::kPermutation, Phase::kEndIndex, Phase::kRoadStarts});
    const double staged = summary.seconds(Phase::kSlices) + summary.seconds(Phase::kRoads);
    LOGP(info,
         "    serial host passes: {:.3f} s of {:.3f} s slices+roads ({:.1f} %); {:.3f} s stage resets, {:.3f} s stage scans, "
         "{:.3f} s permutation, {:.3f} s end index, {:.3f} s road scans",
         serial, staged, staged > 0. ? 100. * serial / staged : 0., resets, scans, summary.seconds(Phase::kPermutation), summary.seconds(Phase::kEndIndex),
         summary.seconds(Phase::kRoadStarts));
  }
  reportTrackTimeBoxes(frame, tracks);
  if (truthRows.empty()) {
    return;
  }

  // The MC particles: every particle on a row, with its longest chain
  // (particleLayers) and the kinds of layers it has; findable by the
  // comparison's rule, in two sets: the physical primaries and the rest.
  struct Particle {
    Scalar pt = 0.F;
    Scalar cotTheta = 0.F;
    Scalar originZ = 0.F;
    std::uint32_t longest = 0U; ///< its longest chain (particleLayers)
    bool barrel = false;        ///< has a barrel layer
    bool disc = false;          ///< has a disc layer
    bool findable = false;
    std::size_t set = 0U;      ///< 0 a physical primary, 1 another particle
    std::uint32_t good = 0U;   ///< good tracks of it
    std::uint8_t foundBy = 0U; ///< the iteration of its first good track
  };
  std::unordered_map<std::uint64_t, Particle> particles;
  for (const auto& [id, layers] : particleLayers(truthRows, findableRule)) {
    const auto found = truth.find(id);
    if (found == truth.end()) {
      continue;
    }
    const auto& t = found->second;
    Particle particle{.pt = t.pt,
                      .cotTheta = static_cast<Scalar>(std::sinh(static_cast<double>(t.eta))),
                      .originZ = t.vz,
                      .findable = isFindable(t.pt, layers, findableRule),
                      .set = t.origin == ParticleOrigin::kPhysicalPrimary ? 0U : 1U};
    particle.longest = layers.chain;
    for (auto bits = layers.layers; bits != 0U; bits &= bits - 1U) {
      (Detector<Det>::isBarrelLayer(static_cast<LayerID>(std::countr_zero(bits))) ? particle.barrel : particle.disc) = true;
    }
    particles.emplace(id, particle);
  }

  // The tracks.
  constexpr std::size_t kCountSlots = 18U;
  std::array<std::uint64_t, kCountSlots> tracksByCount{};
  std::array<std::uint64_t, kCountSlots> fakesByCount{};
  std::array<std::uint64_t, kCountSlots> extensionsByCount{};
  std::array<std::uint64_t, 5> fakePurity{}; // majority fraction: <50, 50-75, 75-90, 90-100, noise-only
  std::uint64_t good = 0U;
  std::uint64_t goodNotFindable = 0U;
  std::uint64_t extendedTracks = 0U;
  std::uint64_t addedClusters = 0U;
  std::uint64_t ownAddedClusters = 0U;
  std::vector<std::uint64_t> tracksByIteration(iterations, 0U);
  std::vector<std::uint64_t> fakesByIteration(iterations, 0U);
  std::uint64_t kinkedTracks = 0U;
  std::uint64_t kinkedFakes = 0U;
  std::array<std::uint64_t, static_cast<std::size_t>(Detector<Det>::nLayers)> kinksByLayer{};
  // Who holds a particle's rows among the fakes: one whose majority is the
  // particle itself -- it was found with a wrong cluster -- or another.
  std::unordered_map<std::uint64_t, bool> heldByOwnFake;
  std::unordered_map<std::uint64_t, bool> heldByOtherFake;
  // Every good track of a particle, in output order: iteration, then best first.
  std::unordered_map<std::uint64_t, std::vector<std::size_t>> goodTracks;
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    const auto count = view.count[track];
    const auto iteration = std::min<std::size_t>(view.iteration[track], iterations - 1U);
    ++tracksByIteration[iteration];
    ++tracksByCount[std::min<std::size_t>(count, kCountSlots - 1U)];
    const std::span<const std::uint32_t> trackRows{view.clusters[track].index, count};
    const std::uint32_t addedInward = view.addedInward[track];
    const std::uint32_t addedOutward = view.addedOutward[track];
    const std::uint32_t added = addedInward + addedOutward;
    ++extensionsByCount[std::min<std::size_t>(added, kCountSlots - 1U)];
    if (added != 0U) {
      ++extendedTracks;
      addedClusters += added;
      const std::uint32_t originalEnd = count - addedOutward;
      const auto owner = majorityParticle(truthRows, trackRows.subspan(addedInward, originalEnd - addedInward)).particle;
      for (std::uint32_t i = 0; i < count; ++i) {
        const bool addedHere = i < addedInward || i >= originalEnd;
        ownAddedClusters += addedHere && owner != kNoTrackParticle && truthRows[trackRows[i]].carries(owner) ? 1U : 0U;
      }
    }
    const auto goodId = goodParticle(truthRows, trackRows);
    const bool one = goodId != kNoTrackParticle;
    if (view.kinks[track] != 0U) {
      ++kinkedTracks;
      kinkedFakes += one ? 0U : 1U;
      const auto kinkRow = view.clusters[track].index[view.kinkCluster[track]];
      ++kinksByLayer[static_cast<std::size_t>(tracking.surface(frame.measurement(kinkRow).surface).layer)];
    }
    if (one) {
      ++good;
      goodTracks[goodId].push_back(track);
      const auto particle = particles.find(goodId);
      if (particle != particles.end() && particle->second.findable) {
        if (particle->second.good == 0U) {
          particle->second.foundBy = static_cast<std::uint8_t>(iteration);
        }
        ++particle->second.good;
      } else {
        ++goodNotFindable;
      }
      continue;
    }
    ++fakesByCount[std::min<std::size_t>(count, kCountSlots - 1U)];
    ++fakesByIteration[iteration];
    const auto [majorityId, majority] = majorityParticle(truthRows, trackRows);
    for (const auto row : trackRows) {
      for (std::uint8_t i = 0; i < truthRows[row].count; ++i) {
        const auto id = truthRows[row].particles[i];
        (id == majorityId ? heldByOwnFake : heldByOtherFake)[id] = true;
      }
    }
    const double purity = static_cast<double>(majority) / static_cast<double>(count);
    std::size_t bin = 3U;
    if (majority == 0U) {
      bin = 4U;
    } else if (purity < 0.5) {
      bin = 0U;
    } else if (purity < 0.75) {
      bin = 1U;
    } else if (purity < 0.9) {
      bin = 2U;
    }
    ++fakePurity[bin];
  }

  // The particles, per set: findable or not, found, cloned; by pT,
  // cot(theta), kind of layers and longest run; the missed by the first stage
  // of the first iteration without a chain of the minimum of their own.
  enum Loss : std::uint8_t { kTracklets,
                             kTriplets,
                             kNeighbours,
                             kRoads,
                             kFit,
                             kSelectionOwnFake,
                             kSelectionOtherFake,
                             kSelectionElse,
                             kLosses };
  struct SetTally {
    explicit SetTally(const VertexConstraint& beamSpot) : byOriginZ{beamSpot} {}
    std::uint64_t withCluster = 0U; ///< above the findable floor with a cluster
    std::uint64_t findable = 0U;
    std::uint64_t found = 0U;
    std::uint64_t clones = 0U;
    std::uint64_t findableAboveFirst = 0U; ///< the funnel's denominator
    std::uint64_t missedBelowFirst = 0U;
    std::array<std::uint64_t, 3> findableByKind{};
    std::array<std::uint64_t, 3> foundByKind{};
    std::array<std::uint64_t, kTrackPtBins> findableByPt{};
    std::array<std::uint64_t, kTrackPtBins> foundByPt{};
    std::array<std::uint64_t, kCotThetaBins> findableByCot{};
    std::array<std::uint64_t, kCotThetaBins> foundByCot{};
    std::array<std::uint64_t, kCountSlots> findableByChain{};
    std::array<std::uint64_t, kCountSlots> foundByChain{};
    std::array<std::uint64_t, kCountSlots> notFindableByChain{};
    std::array<std::uint64_t, kLosses> losses{};
    std::uint64_t lostAtRoadsRaised = 0U;
    // Of those lost at the tracklet stage: their adjacent-layer pairs of their
    // own the first iteration did not keep, by the first test each fails (the
    // funnel's); and the particles with none recorded (from the vertex list).
    std::array<std::uint64_t, kPairLosses> trackletLossPairs{};
    std::uint64_t trackletLossNone = 0U;
    OriginZTally byOriginZ;
  };
  constexpr std::array<const char*, 2> kSetNames{"physical primaries", "other findable"};
  std::array<SetTally, 2> sets{SetTally{beamSpot}, SetTally{beamSpot}};
  std::vector<std::uint64_t> foundByIteration(iterations, 0U);
  for (const auto& [id, particle] : particles) {
    if (static_cast<double>(particle.pt) < findableRule.minPtGeV) {
      continue;
    }
    auto& set = sets[particle.set];
    ++set.withCluster;
    const auto chain = std::min<std::size_t>(particle.longest, kCountSlots - 1U);
    if (!particle.findable) {
      ++set.notFindableByChain[chain];
      continue;
    }
    std::size_t kind = 0U; // barrel only
    if (particle.barrel && particle.disc) {
      kind = 2U;
    } else if (particle.disc) {
      kind = 1U;
    }
    const bool hit = particle.good != 0U;
    ++set.findable;
    set.findableAboveFirst += static_cast<double>(particle.pt) >= first.windows.minPtGeV ? 1U : 0U;
    ++set.findableByKind[kind];
    ++set.findableByPt[binOf(particle.pt, kTrackPtEdges)];
    ++set.findableByCot[binOf(particle.cotTheta, kCotThetaEdges)];
    ++set.findableByChain[chain];
    set.byOriginZ.addFindable(particle.originZ, hit);
    if (hit) {
      ++set.found;
      ++set.foundByKind[kind];
      ++set.foundByPt[binOf(particle.pt, kTrackPtEdges)];
      ++foundByIteration[particle.foundBy];
      ++set.foundByCot[binOf(particle.cotTheta, kCotThetaEdges)];
      ++set.foundByChain[chain];
      set.clones += particle.good - 1U;
      continue;
    }
    // The funnel's stages are the first iteration's: a softer particle is
    // counted apart.
    if (static_cast<double>(particle.pt) < first.windows.minPtGeV) {
      ++set.missedBelowFirst;
      continue;
    }
    // Every findable particle is the funnel's: the same rule.
    const auto survived = funnel != nullptr ? funnel->stages(id, 0U) : std::nullopt;
    if (!survived.has_value()) {
      continue;
    }
    Loss loss = kSelectionElse;
    switch (firstLoss(survived->longest, minLayers, false)) {
      case FunnelStage::kTracklets:
        loss = kTracklets;
        addLostPairs(*survived, set.trackletLossPairs, set.trackletLossNone);
        set.byOriginZ.addTrackletLoss(particle.originZ);
        break;
      case FunnelStage::kTriplets:
        loss = kTriplets;
        break;
      case FunnelStage::kNeighbours:
        loss = kNeighbours;
        break;
      case FunnelStage::kRoads:
        loss = kRoads;
        set.lostAtRoadsRaised +=
          survived->roadLoss == RoadLoss::kRaised || survived->roadLoss == RoadLoss::kRaisedBelow || survived->roadLoss == RoadLoss::kRaisedUnlinked ? 1U : 0U;
        break;
      case FunnelStage::kFit:
        loss = kFit;
        break;
      default:
        if (heldByOwnFake.contains(id)) {
          loss = kSelectionOwnFake;
        } else if (heldByOtherFake.contains(id)) {
          loss = kSelectionOtherFake;
        }
        break;
    }
    ++set.losses[loss];
  }
  for (std::size_t index = 0; index < sets.size(); ++index) {
    std::unordered_map<std::uint64_t, Particle> ofSet;
    for (const auto& [id, particle] : particles) {
      if (particle.findable && particle.set == index) {
        ofSet.emplace(id, particle);
      }
    }
    reportIncomplete(frame, param, interface.propagator().bzKiloGauss, tracking, tracks, truthRows, funnel, ofSet, goodTracks, iterations, kSetNames[index], crossings, scattering, &material);
  }
  const std::uint64_t fakes = tracks.size() - good;
  LOGP(info, "  tracks: {} fakes of {} ({:.2f}%)", fakes, tracks.size(), percent(fakes, tracks.size()));
  constexpr std::array<const char*, 3> kKinds{"barrel layers", "disc layers", "both"};
  for (std::size_t index = 0; index < sets.size(); ++index) {
    const auto& set = sets[index];
    LOGP(info, "  MC, {}: {} above {} GeV with a cluster; findable (a chain of {} layers): {}", kSetNames[index], set.withCluster, findableRule.minPtGeV,
         minLayers, set.findable);
    LOGP(info, "  efficiency, {}: {:.2f}% of the findable, {:.2f}% of all above {} GeV with a cluster; clones {:.2f}%", kSetNames[index],
         percent(set.found, set.findable), percent(set.found, set.withCluster), findableRule.minPtGeV, percent(set.clones, set.found));
    for (std::size_t kind = 0; kind < kKinds.size(); ++kind) {
      LOGP(info, "    {:<13} findable {:>7}, found {:.2f}%", kKinds[kind], set.findableByKind[kind], percent(set.foundByKind[kind], set.findableByKind[kind]));
    }
    std::string byPt;
    for (std::size_t bin = 0; bin < kTrackPtBins; ++bin) {
      if (set.findableByPt[bin] != 0U) {
        byPt += fmt::format(" [{:.2f},{:.2f}) {} {:.2f}%", kTrackPtEdges[bin], kTrackPtEdges[bin + 1U], set.findableByPt[bin],
                            percent(set.foundByPt[bin], set.findableByPt[bin]));
      }
    }
    LOGP(info, "    by pT [GeV]:{}", byPt);
    std::string byCot;
    for (std::size_t bin = 0; bin < kCotThetaBins; ++bin) {
      if (set.findableByCot[bin] != 0U) {
        byCot += fmt::format(" [{:.2f},{:.2f}) {:.2f}%", kCotThetaEdges[bin], kCotThetaEdges[bin + 1U], percent(set.foundByCot[bin], set.findableByCot[bin]));
      }
    }
    LOGP(info, "    by cot(theta):{}", byCot);
    std::string byChain;
    std::string notFindable;
    for (std::size_t chain = 1; chain < kCountSlots; ++chain) {
      if (set.findableByChain[chain] != 0U) {
        byChain += fmt::format(" {}: {} {:.2f}%", chain, set.findableByChain[chain], percent(set.foundByChain[chain], set.findableByChain[chain]));
      }
      if (set.notFindableByChain[chain] != 0U) {
        notFindable += fmt::format(" {}: {}", chain, set.notFindableByChain[chain]);
      }
    }
    LOGP(info, "    findable by their longest chain: count, found:{}", byChain);
    LOGP(info, "    not findable, by their longest chain:{}", notFindable);
  }
  std::string byCount;
  for (std::size_t count = 1; count < kCountSlots; ++count) {
    if (tracksByCount[count] != 0U) {
      byCount += fmt::format(" {}: {} ({:.2f}% fake)", count, tracksByCount[count], percent(fakesByCount[count], tracksByCount[count]));
    }
  }
  LOGP(info, "  tracks by clusters:{}", byCount);
  std::string kinkLayers;
  for (std::size_t layer = 0; layer < kinksByLayer.size(); ++layer) {
    if (kinksByLayer[layer] != 0U) {
      kinkLayers += fmt::format(" {}: {}", layer, kinksByLayer[layer]);
    }
  }
  LOGP(info, "  kinks: {} tracks ({:.2f}%), {:.2f}% of them fake; by the layer of the kink's cluster:{}", kinkedTracks, percent(kinkedTracks, tracks.size()),
       percent(kinkedFakes, kinkedTracks), kinkLayers);
  std::string extensions;
  for (std::size_t count = 1; count < kCountSlots; ++count) {
    if (extensionsByCount[count] != 0U) {
      extensions += fmt::format(" {}: {}", count, extensionsByCount[count]);
    }
  }
  LOGP(info, "  follow: {} tracks, {} clusters, {:.2f} per followed track; tracks by clusters added:{}", extendedTracks, addedClusters,
       extendedTracks == 0U ? 0. : static_cast<double>(addedClusters) / static_cast<double>(extendedTracks), extensions);
  LOGP(info, "  follow MC: {}/{} added clusters ({:.2f}%) belong to the original track's majority particle", ownAddedClusters, addedClusters,
       percent(ownAddedClusters, addedClusters));
  LOGP(info, "  fakes by their majority particle's share: <50% {}, 50-75% {}, 75-90% {}, 90-100% {}, noise only {}", fakePurity[0], fakePurity[1],
       fakePurity[2], fakePurity[3], fakePurity[4]);
  LOGP(info, "  good tracks of a particle not findable (below {} GeV, or no chain of {} layers): {}", findableRule.minPtGeV, minLayers, goodNotFindable);
  for (std::size_t iteration = 0; iteration < result.iterations.size(); ++iteration) {
    const auto settings = workflow::iterationSettings(param, static_cast<int>(iteration));
    std::string shape = fmt::format(", {} clusters or more", std::max(settings.minClusters, 4));
    if (!settings.skipLinks) {
      shape += ", no holes";
    }
    if (settings.layerMask != 0U) {
      shape += ", links within layers";
      for (auto bits = settings.layerMask; bits != 0U; bits &= bits - 1U) {
        shape += fmt::format(" {}", std::countr_zero(bits));
      }
    }
    if (settings.roadEndLayers != 0U) {
      shape += ", roads ending on layers";
      for (auto bits = settings.roadEndLayers; bits != 0U; bits &= bits - 1U) {
        shape += fmt::format(" {}", std::countr_zero(bits));
      }
    }
    std::string floors;
    for (std::size_t count = 0; count < settings.minPtGeVByClusters.size(); ++count) {
      if (settings.minPtGeVByClusters[count] > 0.) {
        floors += fmt::format(" {}: {:.3g}", count, settings.minPtGeVByClusters[count]);
      }
    }
    if (!floors.empty()) {
      shape += ", pT floor by clusters" + floors;
    }
    LOGP(info, "  iteration {} (pT_min {} GeV{}): {} tracks, {:.2f}% fake; first to find {} particles", iteration, settings.windows.minPtGeV, shape,
         tracksByIteration[iteration], percent(fakesByIteration[iteration], tracksByIteration[iteration]), foundByIteration[iteration]);
  }
  // The clones: each further good track B of a particle against its first, A,
  // by the iterations that found them and by where B lies -- inside A's
  // innermost cluster, beyond its outermost (by distance from the origin), on
  // layers A has too, or interleaved with A.
  std::array<std::pair<Frame::Index, Frame::Index>, Frame::kLayerCount> layerRows{};
  for (std::size_t layer = 0; layer < Frame::kLayerCount; ++layer) {
    const auto rows = frame.grid().layerRows(static_cast<LayerID>(layer));
    layerRows[layer] = {rows.begin, rows.end};
  }
  const auto layerOf = [&](Frame::Index row) {
    for (std::size_t layer = 0; layer < Frame::kLayerCount; ++layer) {
      if (row >= layerRows[layer].first && row < layerRows[layer].second) {
        return layer;
      }
    }
    return Frame::kLayerCount;
  };
  const auto distance = [&](Frame::Index row) {
    const auto m = frame.measurement(row);
    return std::hypot(static_cast<double>(m.x), static_cast<double>(m.y), static_cast<double>(m.z));
  };
  enum Relation : std::uint8_t { kInner,
                                 kOuter,
                                 kSameLayers,
                                 kInterleaved,
                                 kRelations };
  // A split -- the clone inside or beyond the first -- by the first stage of
  // the first iteration whose longest surviving chain of the particle is
  // shorter than the two tracks together: there the particle's chain broke.
  enum Split : std::uint8_t {
    kSplitTracklets,
    kSplitTriplets,
    kSplitNeighbours,
    kSplitRoads,
    kSplitFit,
    kSplitElse,
    kSplitBelowFirst,
    kSplitNotFindable,
    kSplits
  };
  std::array<std::uint64_t, kSplits> splits{};
  std::array<std::uint64_t, kPairLosses> splitLostPairs{};
  std::uint64_t splitLostNone = 0U;
  std::vector<std::vector<std::array<std::uint64_t, kRelations>>> cloneRelations(iterations, std::vector<std::array<std::uint64_t, kRelations>>(iterations));
  std::array<std::uint64_t, kCountSlots> cloneClusters{};
  for (const auto& [id, list] : goodTracks) {
    if (list.size() < 2U) {
      continue;
    }
    const auto a = list.front();
    double innerA = std::numeric_limits<double>::max();
    double outerA = 0.;
    std::array<bool, Frame::kLayerCount> layersA{};
    for (std::uint32_t i = 0; i < view.count[a]; ++i) {
      const auto row = view.clusters[a].index[i];
      innerA = std::min(innerA, distance(row));
      outerA = std::max(outerA, distance(row));
      layersA[std::min(layerOf(row), Frame::kLayerCount - 1U)] = true;
    }
    for (std::size_t other = 1; other < list.size(); ++other) {
      const auto b = list[other];
      double innerB = std::numeric_limits<double>::max();
      double outerB = 0.;
      bool sameLayers = true;
      for (std::uint32_t i = 0; i < view.count[b]; ++i) {
        const auto row = view.clusters[b].index[i];
        innerB = std::min(innerB, distance(row));
        outerB = std::max(outerB, distance(row));
        sameLayers &= layersA[std::min(layerOf(row), Frame::kLayerCount - 1U)];
      }
      Relation relation = kInterleaved;
      if (outerB < innerA) {
        relation = kInner;
      } else if (innerB > outerA) {
        relation = kOuter;
      } else if (sameLayers) {
        relation = kSameLayers;
      }
      const auto iterationA = std::min<std::size_t>(view.iteration[a], iterations - 1U);
      const auto iterationB = std::min<std::size_t>(view.iteration[b], iterations - 1U);
      ++cloneRelations[iterationA][iterationB][relation];
      if (relation == kInner || relation == kOuter) {
        const auto need = static_cast<std::uint32_t>(view.count[a]) + view.count[b];
        const auto stage = funnel != nullptr ? funnel->stages(id, 0U) : std::nullopt;
        const auto splitParticle = particles.find(id);
        Split split = kSplitElse;
        if (splitParticle == particles.end() || !splitParticle->second.findable || !stage.has_value()) {
          split = kSplitNotFindable;
        } else if (static_cast<double>(splitParticle->second.pt) < first.windows.minPtGeV) {
          // The first iteration's windows admit a softer particle's short
          // fragment; the iteration it belongs to then finds the rest.
          split = kSplitBelowFirst;
        } else {
          switch (firstLoss(stage->longest, need, false)) {
            case FunnelStage::kTracklets:
              split = kSplitTracklets;
              addLostPairs(*stage, splitLostPairs, splitLostNone);
              break;
            case FunnelStage::kTriplets:
              split = kSplitTriplets;
              break;
            case FunnelStage::kNeighbours:
              split = kSplitNeighbours;
              break;
            case FunnelStage::kRoads:
              split = kSplitRoads;
              break;
            case FunnelStage::kFit:
              split = kSplitFit;
              break;
            default:
              break;
          }
        }
        ++splits[split];
      }
      ++cloneClusters[std::min<std::size_t>(view.count[b], kCountSlots - 1U)];
    }
  }
  for (std::size_t iterationA = 0; iterationA < result.iterations.size(); ++iterationA) {
    for (std::size_t iterationB = iterationA; iterationB < result.iterations.size(); ++iterationB) {
      const auto& relations = cloneRelations[iterationA][iterationB];
      const auto total = relations[kInner] + relations[kOuter] + relations[kSameLayers] + relations[kInterleaved];
      if (total != 0U) {
        LOGP(info, "  clones, first found by iteration {}, clone by {}: {} -- inside the first {}, beyond it {}, on its layers {}, interleaved {}", iterationA,
             iterationB, total, relations[kInner], relations[kOuter], relations[kSameLayers], relations[kInterleaved]);
      }
    }
  }
  std::string cloneCounts;
  for (std::size_t count = 1; count < kCountSlots; ++count) {
    if (cloneClusters[count] != 0U) {
      cloneCounts += fmt::format(" {}: {}", count, cloneClusters[count]);
    }
  }
  LOGP(info, "  clones by their clusters:{}", cloneCounts);
  LOGP(info,
       "  split clones (inside or beyond the first), by the first stage of the first iteration whose longest chain of the particle is "
       "shorter than both tracks together: tracklets {}, triplets {}, neighbours {}, roads {}, road fit {}, none {}; "
       "below the first iteration's pT_min {}; not findable {}",
       splits[kSplitTracklets], splits[kSplitTriplets], splits[kSplitNeighbours], splits[kSplitRoads], splits[kSplitFit], splits[kSplitElse],
       splits[kSplitBelowFirst], splits[kSplitNotFindable]);
  LOGP(info, "    the tracklet splits' adjacent-layer pairs of their own the first iteration did not keep, by the first test:{}",
       describeLostPairs(splitLostPairs, splitLostNone));
  for (std::size_t index = 0; index < sets.size(); ++index) {
    const auto& set = sets[index];
    const auto& losses = set.losses;
    const auto missed = set.findable - set.found - set.missedBelowFirst;
    LOGP(info,
         "  missed {} above the first iteration's pT_min {} ({:.2f}% of their findable): the first stage without a chain of {} clusters "
         "entirely theirs -- tracklets {}, triplets {}, neighbours {}, roads {} ({} raised by another's predecessor), fit {}; lost in the "
         "selection to a fake of their own {}, to another's fake {}, otherwise {}; missed below it {}",
         kSetNames[index], missed, percent(missed, set.findableAboveFirst), minLayers, losses[kTracklets], losses[kTriplets], losses[kNeighbours],
         losses[kRoads], set.lostAtRoadsRaised, losses[kFit], losses[kSelectionOwnFake], losses[kSelectionOtherFake], losses[kSelectionElse],
         set.missedBelowFirst);
    LOGP(info, "    the tracklet losses' adjacent-layer pairs of their own the first iteration did not keep, by the first test:{}",
         describeLostPairs(set.trackletLossPairs, set.trackletLossNone));
    set.byOriginZ.report(kSetNames[index]);
  }
}

void writeTruthBreakdown(std::ostream& output, std::size_t timeframe, const char* vertexSource, const TrackletRunSummary& summary)
{
  for (std::size_t link = 0; link < nLinks<Det>; ++link) {
    for (std::size_t pt = 0; pt < kPtBins; ++pt) {
      for (std::size_t cot = 0; cot < kCotThetaBins; ++cot) {
        const auto total = summary.truthBreakdown[link][pt][cot];
        if (total == 0U) {
          continue;
        }
        output << timeframe << ',' << vertexSource << ',' << static_cast<unsigned int>(Det::layerLinks[link].from) << ','
               << static_cast<unsigned int>(Det::layerLinks[link].to) << ',' << kPtEdges[pt] << ',' << kPtEdges[pt + 1U] << ',' << kCotThetaEdges[cot] << ','
               << kCotThetaEdges[cot + 1U] << ',' << summary.truthBreakdownKept[link][pt][cot] << ',' << total << '\n';
      }
    }
  }
}

/// O2's MC hits for \a D, or none if \a path is not there: the pulls are an
/// addition to the comparison, not a condition of it.
template <DetectorLayout D>
std::unique_ptr<study::HitFile> openHits(const std::string& path, const std::string& discPath)
{
  try {
    return std::make_unique<study::HitFile>(path, study::chipsOf<D>(), discPath);
  } catch (const std::runtime_error& error) {
    LOGP(warning, "no pulls: {}", error.what());
    return nullptr;
  }
}

/// Binds the producer's message after checking its layout and positions
/// against the sink's geometry, and tracks it with the study's observers.
class MeasurementSinkTask : public Task
{
 public:
  MeasurementSinkTask(bool useMC, bool withVertices, TrackerProducts products, std::shared_ptr<workflow::Conditions<Det, true>> conditions)
    : mUseMC(useMC),
      mWithVertices(withVertices),
      mTrackerProducts(products),
      mInterface(withVertices),
      mConditions(std::move(conditions))
  {
  }

  void init(InitContext& ic) final
  {
    // O2's tracks of the same clusters, to put through the same truth rules:
    // the ITS's or the MFT's file; no other detector has one.
    if (const auto o2Tracks = ic.options().get<std::string>("o2-tracks"); !o2Tracks.empty()) {
      mO2Tracks = std::make_unique<study::TrackFile>(o2Tracks, study::chipsOf<Det>());
    }
    // The MC hits: the pulls with --o2-tracks, the true kinks behind the
    // lost pairs and the tracks found incomplete, and the QA's outward legs.
    if (mUseMC) {
      mHits =
        openHits<Det>(ic.options().get<std::string>("mc-hits"), std::is_same_v<Det, detectors::TRKFT3> ? ic.options().get<std::string>("mc-ft3-hits") : "");
    }
    if (mUseMC) {
      const auto contextFile = o2::utils::Str::concat_string(o2::utils::Str::rectifyDirectory(ic.options().get<std::string>("input-dir")),
                                                             ic.options().get<std::string>("collision-context"));
      mContext.reset(o2::steer::DigitizationContext::loadFromFile(contextFile));
      if (!mContext) {
        LOGP(fatal, "cannot read collision context {}", contextFile);
        return;
      }
      mKinematics = std::make_unique<o2::steer::MCKinematicsReader>(mContext.get());
      if (!mKinematics->isInitialized()) {
        LOGP(fatal, "cannot open MC kinematics referenced by {}", contextFile);
        return;
      }
      const auto truthOutput = ic.options().get<std::string>("tracklet-truth-output");
      if (!truthOutput.empty()) {
        mTruthOutput.open(truthOutput, std::ios::out | std::ios::trunc);
        if (!mTruthOutput) {
          LOGP(fatal, "cannot create tracklet truth output {}", truthOutput);
          return;
        }
        mTruthOutput << "timeframe,vertex_source,from,to,pt_min_gev,"
                        "pt_max_gev,cot_theta_min,cot_theta_max,"
                        "kept,true_pairs\n";
      }
    }

    mTruthReports = !ic.options().get<bool>("no-truth-reports");
    mSpeciesOracle = ic.options().get<bool>("species-oracle");
    // The MC truth QA's file.
    if (const auto qaOutput = ic.options().get<std::string>("qa-output"); !qaOutput.empty()) {
      if (mUseMC) {
        mQAOutput = std::make_unique<o2::utils::TreeStreamRedirector>(qaOutput.c_str(), "recreate");
      } else {
        LOGP(warning, "--qa-output {}: without MC there is nothing to fill, no file written", qaOutput);
      }
    }

    mDumpInput = ic.options().get<std::string>("dump-tracking-input");
    // The dump's parameters now, before the conditions: the CCDB's
    // AlpideParam registers a second instance, whose section the INI writer
    // refuses as a duplicate. The configuration does not change after init.
    if (!mDumpInput.empty()) {
      writeParameters(mDumpInput + ".ini");
    }
    mDumpTimeframe = static_cast<std::size_t>(std::max(ic.options().get<int>("dump-timeframe"), 0));
    mDumpOnly = ic.options().get<bool>("dump-only");
    mConditions->init(ic);
  }

  void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj) final
  {
    mConditions->finaliseCCDB(matcher, obj);
  }

  /// The study's observers on the tracker, again whenever the interface
  /// builds one (TrackingInterface.h's `configureFor`).
  void observeTracker()
  {
    auto& tracker = mInterface.tracker();
    // The iteration lines and the stage truth's times read the profile.
    tracker.setProfiling(true);
    // The host tracker's stages: into its stage truth, every iteration
    // against the clusters it found free; with MC into its stage funnel, read
    // by the text reports and the QA's seeded.
    tracker.setIterationObserver([this](const Controller& observed, const Controller::IterationObservation& iteration) {
      if (mObservers.funnel.has_value()) {
        mObservers.funnel->observeIteration(observed, iteration);
      }
      if (mObservers.stage.has_value()) {
        mObservers.summaries.push_back(mObservers.stage->finish(observed, iteration));
        mObservers.stage.reset();
      }
    });
    tracker.setSliceObserver([this](const Controller& observed, const Controller::SliceObservation& slice) {
      if (!mObservers.stage.has_value() && mObservers.openStage) {
        mObservers.openStage(observed, slice.iteration, mObservers.stage);
      }
      if (mObservers.stage.has_value()) {
        mObservers.stage->observeSlice(observed, slice);
      }
      if (mObservers.funnel.has_value()) {
        mObservers.funnel->observeSlice(observed, slice);
      }
    });
  }

  void run(ProcessingContext& pc) final
  {
    mConditions->checkUpdates(pc);
    LOGP(info, "{}: the timeframe at orbit {}", Det::name, pc.services().get<o2::framework::TimingInfo>().firstTForbit);
    // The tracking's timeframe ends on every way out of this one.
    const auto ended = mInterface.endTimeframeOnExit();
    // The trackers are built for the timeframe's beam line and grid, again when they change.
    const auto input = workflow::readTimeframe<Det>(pc, mWithVertices);
    if (mInterface.configureFor(*mConditions, input.shipped)) {
      observeTracker();
    }
    const auto& columns = input.columns;
    const auto offsets = input.offsets;

    const auto view = measurementView(columns);

    const auto n = columns.loc0.size();
    static_assert(std::is_same_v<decltype(view[0].timePS), const std::uint32_t&>, "the view must not be writable");

    const auto& param = recoParam<Det>();
    const auto& geometry = mInterface.grid();
    if (const char* invalid = workflow::TrackingInterface<Det>::checkInput(columns, offsets, static_cast<std::size_t>(geometry.bucketCount()), input.seeds)) {
      LOGP(fatal, "{}", invalid);
      return;
    }

    // Every row's bucket is recomputed from the shipped position alone, as the
    // tracker will see it, and the shipped position -- Cartesian, radius,
    // azimuth and projected variances -- is compared bit for bit with what the
    // geometry gives for the row's own sensor and local measurement.
    const auto sameBits = [](Scalar a, Scalar b) { return std::bit_cast<std::uint32_t>(a) == std::bit_cast<std::uint32_t>(b); };
    const auto& sensors = mInterface.geometry();
    const auto& beamLine = mInterface.beamLine();
    bool bucketsValid = true;
    bool positionsExact = true;
    for (std::uint32_t bucket = 0; bucket < geometry.bucketCount(); ++bucket) {
      for (auto row = offsets[bucket]; row < offsets[bucket + 1U]; ++row) {
        const auto measurement = view[row];
        if (row != offsets[bucket]) {
          bucketsValid &= view[row - 1U].timePS <= measurement.timePS;
        }
        if (measurement.surface >= static_cast<SensorID>(Det::nSensors)) {
          bucketsValid = false;
          continue;
        }
        const auto layer = sensors.sensor(measurement.surface).layer;
        bucketsValid &= geometry.bucketAt(layer, measurement.timePS, measurement.phi, secondCoordinate<Det>(layer, measurement.r, measurement.z)) == bucket;
        const auto point = gridPoint<Det>(sensors, measurement.surface, static_cast<double>(measurement.loc0), static_cast<double>(measurement.loc1),
                                          beamLine.x, beamLine.y);
        const auto variances = pointVariances<Det>(sensors, measurement.surface, static_cast<double>(measurement.loc0), point, measurement.varLoc0,
                                                   measurement.covLoc01, measurement.varLoc1);
        positionsExact &= sameBits(point.x, measurement.x) && sameBits(point.y, measurement.y) && sameBits(point.z, measurement.z) &&
                          sameBits(point.r, measurement.r) && sameBits(point.phi, measurement.phi) && sameBits(variances.z, measurement.varZ) &&
                          sameBits(variances.r, measurement.varR) && sameBits(variances.rPhi, measurement.varRPhi);
      }
    }
    if (!bucketsValid || !positionsExact) {
      LOGP(fatal, "measurements do not match their grid buckets, bucket time order or geometry");
      return;
    }

    const Frame frame{view, geometry, offsets};
    // The vertex list: the input's, or the own vertexer's, which keeps its
    // triplets' claims for the QA with MC.
    const auto found = mInterface.vertices(frame, input.seeds, mUseMC);
    if (found.dropped) {
      ++mTimeframe;
      return;
    }
    const auto* pass = found.pass;
    const std::span<const VertexSeed> seeds = found.seeds;
    if (pass != nullptr) {
      LOGP(info, "vertexer: {} vertices from {} candidates of {} triplets in {:.4f} s", seeds.size(), pass->found.candidates, pass->found.triplets, found.seconds);
    }
    LOGP(info, "bound {} measurements over {} columns, {} layers and {} buckets", n, kNMeasurementColumns, Frame::kLayerCount, geometry.bucketCount());
    for (std::size_t layer = 0; layer < Frame::kLayerCount; ++layer) {
      const auto rows = frame.grid().layerRows(static_cast<LayerID>(layer));
      LOGP(info, "  layer {:>2}  {:>9} measurements  [{:>9}, {:>9})", layer, rows.size(), rows.begin, rows.end);
    }

    const auto blocksPerSlice = geometry.blocksPerSlice(param.sliceWidthPS);
    LOGP(info, "slicing {} blocks in groups of {} with halo {} ps", geometry.blockCount(), blocksPerSlice, param.sliceHaloPS);

    std::vector<TruePair> truePairs;
    const bool dumping = !mDumpInput.empty() && n != 0U && mTimeframe == mDumpTimeframe;
    std::vector<TruthRow> comparisonRows;
    ParticleTable comparisonParticles;
    if (mUseMC) {
      const auto labels = pc.inputs().get<o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel>>("labels");
      // Every row's particles: a cluster is every particle's it carries among
      // its labels, for every report of the sink and the dump's truth.
      comparisonTruth(frame, labels, comparisonRows, comparisonParticles);
      const auto mcOf = [&](std::uint64_t id) -> const o2::MCTrack* {
        const auto found = mParticleLabels.find(id);
        return found == mParticleLabels.end() ? nullptr : mKinematics->getTrack(found->second.source, found->second.event, found->second.track);
      };
      if (mTruthReports) {
        truePairs = mcTruePairs(frame, comparisonRows, mcOf, mInterface.geometry(), param, mInterface.propagator().bzKiloGauss);
        LOGP(info, "  MC-true primary pairs  {}", truePairs.size());
      }
      if (pass) {
        const auto score = scoreVertices<Det>(comparisonRows, comparisonParticles, pass->claims, pass->found);
        LOGP(info, "{}", describe(score));
        LOGP(info, "{}", describePulls(score));
      }
    }
    if (mTrackerProducts == TrackerProducts::kDigests) {
      logTrackerDigests(pc, frame);
    }
    if (dumping) {
      // The tracker's input for o2-gtrk-track-bench, once per process: the
      // configured timeframe (--dump-timeframe), if it has measurements, and
      // with MC the truth to score its tracks by.
      std::vector<double> depths;
      depths.reserve(Detector<Det>::nLayers);
      for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
        depths.push_back(mInterface.geometry().measurementDepth(static_cast<LayerID>(layer)));
      }
      writeTrackingInput(mDumpInput, view, n, offsets, seeds, mInterface.beamLine(), mInterface.frames(), depths);
      LOGP(info, "tracking input written to {} and {}.ini", mDumpInput, mDumpInput);
      if (mUseMC) {
        writeTrackingTruth(mDumpInput + ".truth", std::span<const std::uint32_t>{columns.clusterIndex.data(), n}, comparisonRows, comparisonParticles);
        LOGP(info, "tracking truth written to {}.truth ({} particles)", mDumpInput, comparisonParticles.size());
      }
      mDumpInput.clear();
      if (mDumpOnly) {
        ++mTimeframe;
        return;
      }
    }
    const auto beamSpot = workflow::beamSpotConstraint(mInterface.beamLine());
    const Propagator propagator = mInterface.propagator();
    // The MC hits behind the truth's particles, for the lost pairs' and the
    // incomplete tracks' anatomy.
    const TrueCrossings crossingsOf = mUseMC ? trueCrossings(frame) : TrueCrossings{};
    // The tracker's observers for this timeframe: every iteration's stages
    // against the truth, from its first slice, and with MC the stage funnel.
    mObservers.summaries.clear();
    if (mTruthReports) {
      mObservers.openStage = [&](const Controller& tracker, std::size_t iteration, std::optional<StageTruth>& stage) {
        stage.emplace(frame, tracker, iteration, param, propagator, truePairs, comparisonRows, comparisonParticles, mInterface.scattering(),
                      &mInterface.material(), crossingsOf);
      };
    }
    if (!comparisonRows.empty() && mTruthReports) {
      reportTrueRoadChi2(frame, param, mInterface.tracker(), truePairs, comparisonParticles, mInterface.tracking(), mInterface.material(), propagator, crossingsOf, mInterface.scattering());
    }
    // O2's tracks of this timeframe as rows, once, for the comparison and
    // the QA, with MC.
    const auto theirs = mO2Tracks && !comparisonRows.empty() ? o2TrackRows(frame) : TrackRows{};
    // --qa-output: the tracker's tracks into the QA, after what else reads
    // them.
    using Tracked = std::function<void(const TrackRefit<Det>::Tracks&)>;
    const auto withQA = [&](Tracked tracked) -> Tracked {
      if (!mQAOutput || comparisonRows.empty()) {
        return tracked;
      }
      return [&, tracked = std::move(tracked)](const TrackRefit<Det>::Tracks& tracks) {
        if (tracked) {
          tracked(tracks);
        }
        const auto& funnel = mObservers.funnel;
        streamQA(frame, comparisonRows, comparisonParticles, funnel.has_value() ? &*funnel : nullptr, tracks, mO2Tracks ? &theirs : nullptr);
      };
    };
    const auto comparisons = [&]() -> Tracked {
      if (!comparisonRows.empty()) {
        return [&](const TrackRefit<Det>::Tracks& tracks) {
          if (mO2Tracks) {
            compareWithO2(frame, comparisonRows, comparisonParticles, tracks, theirs);
          }
          if (!mTruthReports) {
            return;
          }
          const auto& funnel = mObservers.funnel;
          reportFakeAnatomy(frame, param, comparisonRows, comparisonParticles, funnel.has_value() ? &*funnel : nullptr, mInterface.tracker(), tracks, mInterface.tracking(),
                            mInterface.material(), propagator, mInterface.trackerLabel().c_str());
        };
      }
      if (mO2Tracks && !mUseMC) {
        // Real data: no truth, O2's tracks of the same clusters are the yardstick.
        return [&](const TrackRefit<Det>::Tracks& tracks) { compareWithO2ByClusters(frame, tracks); };
      }
      return {};
    };
    // The tracker's stage funnel of this timeframe, with MC.
    const auto makeFunnel = [&](std::optional<Funnel>& funnel) {
      if (comparisonRows.empty()) {
        return;
      }
      std::vector<double> minPtGeV(static_cast<std::size_t>(workflow::iterationCount(param)));
      for (std::size_t iteration = 0; iteration < minPtGeV.size(); ++iteration) {
        minPtGeV[iteration] = workflow::iterationSettings(param, static_cast<int>(iteration)).windows.minPtGeV;
      }
      funnel.emplace(frame, comparisonRows, comparisonParticles, comparisonSettings<Det>(), vertexSeed(beamSpot), std::move(minPtGeV));
      // The refused roads refitted for their reason, as the bench's funnel.
      funnel->setRoadFit(mInterface.tracking(), mInterface.material(), propagator);
    };
    const auto funnelOf = [](const std::optional<Funnel>& funnel) { return funnel.has_value() ? &*funnel : nullptr; };
    // --compare-tracks: the device's products against this tracker's, before
    // what else reads the tracks.
    const auto withDevice = [&](Tracked tracked) -> Tracked {
      if (mTrackerProducts != TrackerProducts::kCompare) {
        return tracked;
      }
      return [&, tracked = std::move(tracked)](const TrackRefit<Det>::Tracks& tracks) {
        compareWithDevice(pc, tracks, pass != nullptr ? &pass->found : nullptr);
        if (tracked) {
          tracked(tracks);
        }
      };
    };
    makeFunnel(mObservers.funnel);
    measureTracks(frame, param, beamSpot, seeds, comparisonRows, comparisonParticles,
                  funnelOf(mObservers.funnel), mInterface, withDevice(withQA(comparisons())), crossingsOf);
    if (mTruthOutput) {
      // The first iteration's.
      if (!mObservers.summaries.empty()) {
        writeTruthBreakdown(mTruthOutput, mTimeframe, mInterface.trackerSource() == VertexSource::kBeamLine ? "beam_line" : "vertex_list", mObservers.summaries[0]);
      }
      mTruthOutput.flush();
    }
    mObservers.funnel.reset();
    mObservers.stage.reset();
    mObservers.openStage = nullptr;
    ++mTimeframe;
  }

  void endOfStream(EndOfStreamContext& /*ec*/) final
  {
    if (!mQAOutput) {
      return;
    }
    const std::string name = mQAOutput->GetFile()->GetName();
    if (mInterface.configured()) {
      // The MC events: the signal source's, as the truth is.
      streamRun(*mQAOutput, workflow::iterationCount(recoParam<Det>()), etaReach(PlacedEnvelopes<Det>::all()), comparisonSettings<Det>(),
                mKinematics ? mKinematics->getNEvents(0) : 0U);
    }
    mQAOutput->Close();
    mQAOutput.reset();
    LOGP(info, "the QA written to {}", name);
  }

 private:
  /// --tracker-digests: the tracking device's products of this timeframe
  /// (TrackerSpec.h) on \a frame, digested as the bench digests a device's
  /// tracks and its vertexer's list. One line per timeframe.
  void logTrackerDigests(ProcessingContext& pc, const Frame& frame)
  {
    using Record = workflow::TrackRecord<Det>;
    const auto records = pc.inputs().get<gsl::span<Record>>("deviceTracks");
    const auto vertices = pc.inputs().get<gsl::span<VertexSeed>>("deviceVertices");
    const auto contributors = pc.inputs().get<gsl::span<std::uint32_t>>("deviceContributors");
    BoundedMemoryResource resource;
    const auto tracks = recordsTracks<Det>(std::span<const Record>{records.data(), records.size()}, &resource);
    LOGP(info, "tracker-digests, {}, orbit {}: {} tracks, {} vertices; digest, tracks: {:016x}; digest, tracks in clusters: {:016x}; digest, vertices: {:016x}",
         Det::name, pc.services().get<o2::framework::TimingInfo>().firstTForbit, tracks.size(), vertices.size(), tracksDigest(tracks),
         tracksClusterDigest<Det>(tracks, frame),
         vertexDigest(std::span<const VertexSeed>{vertices.data(), vertices.size()}, std::span<const std::uint32_t>{contributors.data(), contributors.size()}));
  }

  /// --compare-tracks: the tracking device's products of this timeframe
  /// (TrackerSpec.h) against this tracker's \a tracks and the own vertexer's
  /// \a found (none without it): every record field by field, failed refits
  /// included, every vertex and contributor count. One line per timeframe;
  /// the run fails on the first difference.
  void compareWithDevice(ProcessingContext& pc, const typename TrackRefit<Det>::Tracks& tracks, const VertexerResult* found)
  {
    using Record = workflow::TrackRecord<Det>;
    const auto theirs = pc.inputs().get<gsl::span<Record>>("deviceTracks");
    const auto theirVertices = pc.inputs().get<gsl::span<VertexSeed>>("deviceVertices");
    const auto theirContributors = pc.inputs().get<gsl::span<std::uint32_t>>("deviceContributors");
    const std::span<const VertexSeed> vertices = found != nullptr ? std::span<const VertexSeed>{found->vertices} : std::span<const VertexSeed>{};
    const std::span<const std::uint32_t> contributors = found != nullptr ? std::span<const std::uint32_t>{found->contributors} : std::span<const std::uint32_t>{};
    const auto orbit = pc.services().get<o2::framework::TimingInfo>().firstTForbit;
    std::string difference;
    if (theirs.size() != tracks.size()) {
      difference = fmt::format("the device has {} tracks", theirs.size());
    } else {
      for (std::size_t track = 0; track < tracks.size(); ++track) {
        if (const auto field = differingField(tracks.read(track), theirs[track]); field >= 0) {
          difference = fmt::format("track {}, field {} (memlayout order)", track, field);
          break;
        }
      }
    }
    const auto sameBytes = []<class T>(std::span<const T> mine, gsl::span<const T> other) {
      return mine.size() == other.size() && (mine.empty() || std::memcmp(mine.data(), other.data(), mine.size_bytes()) == 0);
    };
    if (difference.empty() && !sameBytes(vertices, theirVertices)) {
      difference = fmt::format("the vertices (the device {})", theirVertices.size());
    }
    if (difference.empty() && !sameBytes(contributors, theirContributors)) {
      difference = "the vertices' contributors";
    }
    if (!difference.empty()) {
      LOGP(fatal, "compare-tracks, {}, orbit {}: {} tracks, {} vertices: differ, {}", Det::name, orbit, tracks.size(), vertices.size(), difference);
      return;
    }
    LOGP(info, "compare-tracks, {}, orbit {}: {} tracks, {} vertices: equal", Det::name, orbit, tracks.size(), vertices.size());
  }

  /// The first field of \a mine whose bytes differ from \a theirs, in
  /// memlayout's column order; -1 if none. Field by field: the record's
  /// padding between fields is not compared, and no field has padding of
  /// its own (BoundTrackState.h).
  template <class Record>
  [[nodiscard]] static int differingField(const Record& mine, const Record& theirs)
  {
    constexpr auto kFields = static_cast<std::size_t>(Record::Column::ColumnCount);
    using Fields = std::array<std::pair<const void*, std::size_t>, kFields>;
    const auto fieldsOf = [](const Record& record) {
      Fields fields{};
      record.apply([&](const auto&... field) {
        std::size_t index = 0;
        ((fields[index++] = {static_cast<const void*>(&field), sizeof(field)}), ...);
        return 0;
      });
      return fields;
    };
    const auto left = fieldsOf(mine);
    const auto right = fieldsOf(theirs);
    for (std::size_t field = 0; field < kFields; ++field) {
      if (std::memcmp(left[field].first, right[field].first, left[field].second) != 0) {
        return static_cast<int>(field);
      }
    }
    return -1;
  }

  /// The tracker's tracks of this timeframe into the QA's trees;
  /// \a theirs O2's (`o2TrackRows`), none
  /// without --o2-tracks.
  void streamQA(const Frame& frame, std::span<const TruthRow> rows, const ParticleTable& particles, const Funnel* funnel,
                const typename TrackRefit<Det>::Tracks& tracks, const TrackRows* theirs)
  {
    const auto view = tracks.view();
    const TrackRows ours = trackRows(tracks);
    std::vector<TrackKinematics> kinematics(tracks.size());
    std::vector<QAResolution> resolution;
    std::vector<QAOracle> oracle;
    const Propagator propagator = mInterface.propagator();
    for (std::size_t track = 0; track < tracks.size(); ++track) {
      const auto& state = view.state[track];
      const auto theta = static_cast<double>(state.theta());
      kinematics[track] = {.pt = static_cast<float>(std::sin(theta) / std::abs(static_cast<double>(state.qOverP()))),
                           .eta = static_cast<float>(-std::log(std::tan(theta / 2.))),
                           .phi = state.phi(),
                           .iteration = view.iteration[track],
                           .chi2 = view.chi2[track],
                           .ndf = view.ndf[track]};
      // A good track against its particle's crossing, with the hits.
      const auto particle = goodParticle(rows, ours[track]);
      const auto innermost = std::ranges::find_if(ours[track], [&](std::uint32_t row) { return frame.measurement(row).surface == state.surface; });
      const auto found = particles.find(particle);
      if (!mHits || particle == kNoTrackParticle || innermost == ours[track].end() || found == particles.end()) {
        continue;
      }
      const auto truth = stateTruth(frame, state, particle, *innermost);
      if (!truth.has_value()) {
        continue;
      }
      QAResolution entry{.pt = found->second.pt,
                         .eta = found->second.eta,
                         .phi = wrapPhi(static_cast<double>(found->second.phi)),
                         .charge = static_cast<std::int8_t>(std::lround(mParticleLabels[particle].charge)),
                         .pdg = found->second.pdg,
                         .origin = static_cast<std::uint8_t>(found->second.origin),
                         .iteration = view.iteration[track],
                         .kinks = view.kinks[track],
                         .clusters = view.count[track],
                         .layer = static_cast<std::uint8_t>(rows[*innermost].layer),
                         .fitPt = static_cast<float>(truth->fittedPt),
                         .ptResidual = static_cast<float>((truth->fittedPt / truth->pt) - 1.),
                         .ptSigma = static_cast<float>(truth->ptSigma)};
      for (std::size_t i = 0; i < 5U; ++i) {
        entry.residual[i] = static_cast<float>(truth->residual[i]);
        entry.sigma[i] = static_cast<float>(truth->sigma[i]);
      }
      vertexApproach(frame, state, particle, *innermost, propagator, found->second, entry);
      if (mSpeciesOracle) {
        speciesOracle(frame, view, track, particle, *innermost, propagator, entry, oracle);
      }
      resolution.push_back(entry);
    }
    // Each particle's decay geometry and the funnel's verdict,
    // against the beam line in the
    // run's nominal field.
    std::unordered_set<std::uint64_t> found;
    for (const auto& track : ours) {
      found.insert(goodParticle(rows, track));
    }
    const double bz = mInterface.propagator().bzKiloGauss;
    const double beamX = mInterface.beamLine().x;
    const double beamY = mInterface.beamLine().y;
    // Each signal particle's outward leg in its MC hits (TrackingQA.h), its
    // event's segments ordered by track once.
    std::unordered_map<std::uint64_t, std::uint8_t> outward;
    std::unordered_map<std::uint64_t, std::vector<float>> legMomenta; // --species-oracle
    if (mHits) {
      std::unordered_map<int, std::vector<std::uint64_t>> byEvent;
      for (const auto& [id, truth] : particles) {
        if (const auto label = mParticleLabels.find(id); label != mParticleLabels.end() && label->second.source == 0) {
          byEvent[label->second.event].push_back(id);
        }
      }
      const auto trackOf = [](const study::HitSegment* segment) { return segment->track; };
      std::vector<const study::HitSegment*> segments;
      std::vector<HitCrossing> crossings;
      for (const auto& [event, ids] : byEvent) {
        segments.clear();
        for (const auto& segment : mHits->event(event)) {
          segments.push_back(&segment);
        }
        std::ranges::sort(segments, {}, trackOf);
        for (const auto id : ids) {
          crossings.clear();
          for (const auto* segment : std::ranges::equal_range(segments, mParticleLabels[id].track, {}, trackOf)) {
            crossings.push_back({.time = segment->time,
                                 .position = {0.5 * (segment->start[0] + segment->end[0]), 0.5 * (segment->start[1] + segment->end[1]),
                                              0.5 * (segment->start[2] + segment->end[2])},
                                 .layer = mInterface.geometry().sensor(segment->chip).layer,
                                 .momentum = segment->momentum});
          }
          const auto& truth = particles.at(id);
          const auto leg = outwardLeg(crossings, {static_cast<double>(truth.vx), static_cast<double>(truth.vy), static_cast<double>(truth.vz)});
          outward[id] = outwardLayers(leg);
          if (mSpeciesOracle && (std::abs(truth.pdg) == 11 || (std::abs(truth.pdg) == 211 && id % 16U == 0U))) {
            legMomenta[id] = momentaAlongLeg(leg);
          }
        }
      }
    }
    ParticleExtrasTable extras;
    for (const auto& [id, truth] : particles) {
      const auto label = mParticleLabels.find(id);
      if (label == mParticleLabels.end()) {
        continue;
      }
      const auto& helix = label->second.helix;
      auto& extra = extras[id];
      if (const auto leg = outward.find(id); leg != outward.end()) {
        extra.outwardLayers = leg->second;
      }
      extra.vr = static_cast<float>(std::hypot(helix.position[0] - beamX, helix.position[1] - beamY));
      const auto impact = helixImpact(helix, bz, beamX, beamY);
      extra.d0 = static_cast<float>(impact.d0);
      extra.z0 = static_cast<float>(impact.z - mInterface.beamLine().z);
      if (const auto& mother = label->second.mother; mother.has_value()) {
        const auto direction = helixDirectionAt(*mother, bz, helix.position[0], helix.position[1]);
        const double p = std::hypot(helix.momentum[0], helix.momentum[1], helix.momentum[2]);
        const double cosine = ((direction[0] * helix.momentum[0]) + (direction[1] * helix.momentum[1]) + (direction[2] * helix.momentum[2])) / p;
        extra.openingAngle = static_cast<float>(std::acos(std::clamp(cosine, -1., 1.)));
        extra.motherMiss = static_cast<float>(helixMiss(*mother, bz, helix.position[0], helix.position[1]));
      }
      const auto verdict = funnel != nullptr ? funnel->particleVerdict(id) : std::nullopt;
      if (found.contains(id)) {
        extra.lostAt = static_cast<std::uint8_t>(FunnelStage::kFound);
      } else if (verdict.has_value()) {
        extra.lostAt = static_cast<std::uint8_t>(verdict->lostAt);
        // Why, where it was lost.
        if (verdict->lostAt == FunnelStage::kRoads && verdict->roadLoss.has_value()) {
          extra.roadLoss = static_cast<std::uint8_t>(*verdict->roadLoss);
          extra.roadLosses = verdict->roadLosses;
          if (verdict->raisedFit.has_value()) {
            extra.raisedFit = static_cast<std::uint8_t>(*verdict->raisedFit);
          }
        }
        if (verdict->lostAt == FunnelStage::kFit && verdict->fitLoss.has_value()) {
          extra.fitLoss = static_cast<std::uint8_t>(*verdict->fitLoss);
          extra.fitRoadLayers = verdict->fitRoadLayers;
          extra.fitLosses = verdict->fitLosses;
        }
      }
      if (verdict.has_value() && verdict->pairLoss.has_value()) {
        extra.pairLoss = static_cast<std::uint8_t>(*verdict->pairLoss);
      }
    }
    stream(*mQAOutput, trackingQA(rows, particles, ours, kinematics, funnel != nullptr ? funnel->seededBy() : SeededBy{}, extras, theirs,
                                  workflow::iterationCount(recoParam<Det>()), comparisonSettings<Det>()));
    streamResolution(*mQAOutput, resolution);
    streamOracle(*mQAOutput, oracle);
    if (!legMomenta.empty()) {
      const auto settings = comparisonSettings<Det>();
      const auto chains = particleLayers(rows, settings);
      // Each followed particle's fakes (StageFunnel::report's label rule: the
      // particle most of a fake's measurements carry).
      std::unordered_map<std::uint64_t, std::pair<std::uint8_t, std::uint8_t>> fakes; // fewest wrong, its layer
      for (const auto& track : ours) {
        if (goodParticle(rows, track) != kNoTrackParticle) {
          continue;
        }
        const auto [label, carried] = majorityParticle(rows, track);
        const std::size_t most = carried;
        if (most == 0U || most == track.size() || !legMomenta.contains(label)) {
          continue;
        }
        const auto wrong = static_cast<std::uint8_t>(track.size() - most);
        std::uint8_t layer = 0xFFU;
        if (wrong == 1U) {
          const auto row = *std::ranges::find_if(track, [&](std::uint32_t r) { return !rows[r].carries(label); });
          layer = static_cast<std::uint8_t>(rows[row].layer);
        }
        auto [entry, inserted] = fakes.try_emplace(label, wrong, layer);
        if (!inserted && wrong < entry->second.first) {
          entry->second = {wrong, layer};
        }
      }
      std::vector<QAMomenta> momenta;
      for (auto& [id, pAtLayer] : legMomenta) {
        const auto& truth = particles.at(id);
        const auto chain = chains.find(id);
        if (chain == chains.end() || !isFindable(truth.pt, chain->second, settings)) {
          continue;
        }
        const auto& extra = extras.at(id);
        momenta.push_back({.pdg = truth.pdg,
                           .origin = static_cast<std::uint8_t>(truth.origin),
                           .pt = truth.pt,
                           .eta = truth.eta,
                           .p0 = static_cast<float>(static_cast<double>(truth.pt) * std::cosh(static_cast<double>(truth.eta))),
                           .found = found.contains(id),
                           .lostAt = extra.lostAt,
                           .pairLoss = extra.pairLoss,
                           .layers = chain->second.layers,
                           .fakeWrong = fakes.contains(id) ? fakes.at(id).first : std::uint8_t{0},
                           .fakeWrongLayer = fakes.contains(id) ? fakes.at(id).second : std::uint8_t{0xFF},
                           .pAtLayer = std::move(pAtLayer)});
      }
      streamMomenta(*mQAOutput, momenta);
    }
  }

  /// |p| entering each layer's first crossing on the outward \a leg
  /// (`outwardLeg`); NaN for a layer without one [GeV].
  [[nodiscard]] static std::vector<float> momentaAlongLeg(std::span<const HitCrossing> leg)
  {
    std::vector<float> pAtLayer(static_cast<std::size_t>(Detector<Det>::nLayers), std::numeric_limits<float>::quiet_NaN());
    for (const auto& crossing : leg) {
      auto& p = pAtLayer.at(crossing.layer);
      if (std::isnan(p)) {
        p = static_cast<float>(std::hypot(crossing.momentum[0], crossing.momentum[1], crossing.momentum[2]));
      }
    }
    return pAtLayer;
  }

  /// Our beam-line tracks and O2's for this timeframe, under one set of truth
  /// rules and cluster set by cluster set (TrackComparison.h). O2's tracks
  /// name O2's cluster indices, which are our measurements' clusterIndex.
  /// Every row's signal particles, all of them, and those particles' truth:
  /// what O2's label rule counts (TrackComparison.h).
  void comparisonTruth(const Frame& frame, const o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel>& labels, std::vector<TruthRow>& rows,
                       ParticleTable& particles)
  {
    rows.assign(frame.size(), {});
    mParticleLabels.clear();
    for (Frame::Index row = 0U; row < frame.size(); ++row) {
      const auto measurement = frame.measurement(row);
      auto& truth = rows[row];
      truth.layer = mInterface.geometry().sensor(measurement.surface).layer;
      truth.r = measurement.r;
      truth.z = measurement.z;
      const auto all = labels.getLabels(measurement.clusterIndex);
      static_cast<workflow::RowParticles&>(truth) = workflow::rowParticles(std::span<const o2::MCCompLabel>{all.data(), all.size()});
      // The particle table: each particle the rows carry, once.
      for (const auto& label : all) {
        if (!carriesParticle(label)) {
          continue;
        }
        const auto id = label.getTrackEventSourceID(); // source, event and track in one word
        if (!truth.carries(id)) {
          continue; // past the rule's kMaxParticles
        }
        if (!particles.contains(id)) {
          const auto* mc = mKinematics->getTrack(label);
          const o2::MCTrack* mother = nullptr;
          if (mc == nullptr) {
            particles[id] = ParticleTruth{};
          } else {
            const auto& event = mKinematics->getTracks(label.getSourceID(), label.getEventID());
            mother = o2::mcutils::MCTrackNavigator::getMother(*mc, event);
            particles[id] = ParticleTruth{.pt = static_cast<float>(mc->GetPt()),
                                          .primary = mc->isPrimary(),
                                          .eta = static_cast<float>(mc->GetEta()),
                                          .phi = static_cast<float>(mc->GetPhi()),
                                          .origin = particleOrigin(*mc, event),
                                          .pdg = mc->GetPdgCode(),
                                          .motherPdg = mother == nullptr ? 0 : mother->GetPdgCode(),
                                          .process = static_cast<std::uint8_t>(mc->getProcess()),
                                          .vx = static_cast<float>(mc->Vx()),
                                          .vy = static_cast<float>(mc->Vy()),
                                          .vz = static_cast<float>(mc->Vz())};
          }
          const auto* pdg = mc == nullptr ? nullptr : o2::O2DatabasePDG::Instance()->GetParticle(mc->GetPdgCode());
          auto& entry = mParticleLabels[id];
          entry = {.event = label.getEventID(),
                   .track = label.getTrackID(),
                   .source = label.getSourceID(),
                   .charge = pdg == nullptr ? 0. : pdg->Charge() / 3.,
                   .pdg = mc == nullptr ? 0 : mc->GetPdgCode()};
          if (mc != nullptr) {
            const auto helixOf = [](const o2::MCTrack& track, double charge) {
              return TruthHelix{.position = {track.Vx(), track.Vy(), track.Vz()}, .momentum = {track.Px(), track.Py(), track.Pz()}, .charge = charge};
            };
            entry.helix = helixOf(*mc, entry.charge);
            if (mother != nullptr) {
              const auto* motherPdg = o2::O2DatabasePDG::Instance()->GetParticle(mother->GetPdgCode());
              entry.mother = helixOf(*mother, motherPdg == nullptr ? 0. : motherPdg->Charge() / 3.);
            }
          }
        }
      }
    }
  }

  /// Without MC truth: GTRK's tracks against
  /// O2's of the same clusters by the clusters they share. An O2 track is
  /// found when a GTRK track holds four or more of its clusters, and
  /// reproduced when one holds exactly its clusters; a GTRK track is O2's
  /// when it shares four or more with an O2 track, and GTRK's own otherwise.
  /// Both by the track's cluster count.
  void compareWithO2ByClusters(const Frame& frame, const typename TrackRefit<Det>::Tracks& tracks)
  {
    constexpr std::uint32_t kSharedToMatch = 4U;
    const auto rowOfCluster = rowsOfClusters(frame);
    // Which GTRK tracks hold each row.
    const auto view = tracks.view();
    std::vector<std::vector<std::uint32_t>> tracksOfRow(frame.size());
    std::vector<std::uint32_t> ourCount(tracks.size(), 0U);
    for (std::size_t track = 0; track < tracks.size(); ++track) {
      ourCount[track] = view.count[track];
      for (std::uint32_t i = 0; i < view.count[track]; ++i) {
        tracksOfRow[view.clusters[track].index[i]].push_back(static_cast<std::uint32_t>(track));
      }
    }
    const auto theirs = mO2Tracks->read(static_cast<long long>(mTimeframe));
    constexpr auto kMaxClusters = static_cast<std::size_t>(Detector<Det>::nLayers) + 1U;
    std::array<std::size_t, kMaxClusters> o2Total{};
    std::array<std::size_t, kMaxClusters> o2Found{};
    std::array<std::size_t, kMaxClusters> o2Reproduced{};
    std::vector<std::uint32_t> bestShared(tracks.size(), 0U); // per GTRK track, the most clusters shared with any O2 track
    std::vector<std::uint32_t> sharedWith;                    // scratch: per GTRK track, the clusters shared with the O2 track at hand
    sharedWith.assign(tracks.size(), 0U);
    std::vector<std::uint32_t> touched;
    std::size_t unknownClusters = 0U;
    // Of O2's full-length tracks that GTRK finds shorter: the layers GTRK's
    // track lacks, and GTRK's track length.
    std::array<std::size_t, kMaxClusters> lackingLayer{};
    std::array<std::size_t, kMaxClusters> shorterLength{};
    std::size_t fullFoundShorter = 0U;
    std::size_t shorterLackingTaken = 0U; // of those, with a lacking cluster held by another GTRK track
    std::size_t shorterLackingFree = 0U;  // with every lacking cluster free
    // O2's full-length tracks by their pT and |tan lambda| at the innermost
    // cluster (the ITS only): found full-length, found shorter, missed.
    constexpr std::array<double, 7> kO2PtEdges{0., 0.1, 0.2, 0.3, 0.5, 1., std::numeric_limits<double>::infinity()};
    constexpr std::array<double, 6> kO2TglEdges{0., 0.3, 0.6, 0.9, 1.2, std::numeric_limits<double>::infinity()};
    std::array<std::array<std::size_t, 3>, kO2PtEdges.size() - 1U> fullByPt{};
    std::array<std::array<std::size_t, 3>, kO2TglEdges.size() - 1U> fullByTgl{};
    std::vector<o2::track::TrackParCov> states;
    if constexpr (std::same_as<Det, detectors::ITS>) {
      states = mO2Tracks->innerStates(static_cast<long long>(mTimeframe));
    }
    const auto binOf = [](const auto& edges, double value) {
      std::size_t bin = 0U;
      while (bin + 2U < edges.size() && value >= edges[bin + 1U]) {
        ++bin;
      }
      return bin;
    };
    const auto layerOfRow = [&](std::uint32_t row) { return static_cast<std::size_t>(mInterface.geometry().sensor(frame.measurement(row).surface).layer); };
    for (std::size_t theirIndex = 0; theirIndex < theirs.size(); ++theirIndex) {
      const auto& track = theirs[theirIndex];
      touched.clear();
      std::uint32_t best = 0U;
      std::uint32_t bestTrack = std::numeric_limits<std::uint32_t>::max();
      for (const auto cluster : track) {
        if (cluster >= rowOfCluster.size() || rowOfCluster[cluster] == std::numeric_limits<std::uint32_t>::max()) {
          ++unknownClusters;
          continue;
        }
        for (const auto ours : tracksOfRow[rowOfCluster[cluster]]) {
          if (sharedWith[ours]++ == 0U) {
            touched.push_back(ours);
          }
          if (sharedWith[ours] > best) {
            best = sharedWith[ours];
            bestTrack = ours;
          }
        }
      }
      const auto n = std::min(track.size(), kMaxClusters - 1U);
      ++o2Total[n];
      if (track.size() == static_cast<std::size_t>(Detector<Det>::nLayers) && theirIndex < states.size()) {
        const auto& state = states[theirIndex];
        std::size_t outcome = 0U; // full, shorter, missed
        if (best < kSharedToMatch) {
          outcome = 2U;
        } else if (ourCount[bestTrack] < track.size()) {
          outcome = 1U;
        }
        const auto q2pt = std::abs(static_cast<double>(state.getParam(4)));
        const auto pt = q2pt > 0. ? 1. / q2pt : std::numeric_limits<double>::infinity();
        ++fullByPt[binOf(kO2PtEdges, pt)][outcome];
        ++fullByTgl[binOf(kO2TglEdges, std::abs(static_cast<double>(state.getParam(3))))][outcome];
      }
      if (best >= kSharedToMatch) {
        ++o2Found[n];
        if (track.size() == static_cast<std::size_t>(Detector<Det>::nLayers) && ourCount[bestTrack] < track.size()) {
          ++fullFoundShorter;
          ++shorterLength[std::min<std::size_t>(ourCount[bestTrack], kMaxClusters - 1U)];
          bool taken = false;
          for (const auto cluster : track) {
            if (cluster >= rowOfCluster.size() || rowOfCluster[cluster] == std::numeric_limits<std::uint32_t>::max()) {
              continue; // counted in unknownClusters above
            }
            const auto row = rowOfCluster[cluster];
            const auto& ours = tracksOfRow[row];
            if (std::ranges::find(ours, bestTrack) == ours.end()) {
              ++lackingLayer[std::min(layerOfRow(row), kMaxClusters - 1U)];
              taken |= !ours.empty();
            }
          }
          ++(taken ? shorterLackingTaken : shorterLackingFree);
        }
      }
      for (const auto ours : touched) {
        if (sharedWith[ours] == track.size() && sharedWith[ours] == ourCount[ours]) {
          ++o2Reproduced[n];
          break;
        }
      }
      for (const auto ours : touched) {
        bestShared[ours] = std::max(bestShared[ours], sharedWith[ours]);
        sharedWith[ours] = 0U;
      }
    }
    std::array<std::size_t, kMaxClusters> ourTotal{};
    std::array<std::size_t, kMaxClusters> ourMatched{};
    std::map<std::uint32_t, std::size_t> shortPatterns; // GTRK's tracks short of a layer, by their layer mask
    for (std::size_t track = 0; track < tracks.size(); ++track) {
      const auto n = std::min<std::size_t>(ourCount[track], kMaxClusters - 1U);
      ++ourTotal[n];
      if (bestShared[track] >= kSharedToMatch) {
        ++ourMatched[n];
      }
      if (ourCount[track] < static_cast<std::uint32_t>(Detector<Det>::nLayers)) {
        std::uint32_t mask = 0U;
        for (std::uint32_t i = 0; i < view.count[track]; ++i) {
          mask |= 1U << layerOfRow(view.clusters[track].index[i]);
        }
        ++shortPatterns[mask];
      }
    }
    const auto sum = [](const auto& counts) { return std::accumulate(counts.begin(), counts.end(), std::size_t{0}); };
    const auto byClusters = [&](const auto& part, const auto& whole) {
      std::string line;
      for (std::size_t n = kSharedToMatch; n < kMaxClusters; ++n) {
        if (whole[n] != 0U) {
          line += std::format(" {}: {}/{} ({:.1f})", n, part[n], whole[n], percent(part[n], whole[n]));
        }
      }
      return line;
    };
    LOGP(info, "  against O2's tracks of these clusters, by shared clusters (no truth; {} or more shared is a match):", kSharedToMatch);
    LOGP(info, "    O2: {} tracks, {} found by GTRK ({:.2f} %), {} reproduced exactly ({:.2f} %); found by clusters (%):{}", theirs.size(), sum(o2Found),
         percent(sum(o2Found), theirs.size()), sum(o2Reproduced), percent(sum(o2Reproduced), theirs.size()), byClusters(o2Found, o2Total));
    LOGP(info, "    GTRK: {} tracks, {} with an O2 track ({:.2f} %), {} GTRK's own ({:.2f} %); with O2 by clusters (%):{}", tracks.size(), sum(ourMatched),
         percent(sum(ourMatched), tracks.size()), tracks.size() - sum(ourMatched), percent(tracks.size() - sum(ourMatched), tracks.size()),
         byClusters(ourMatched, ourTotal));
    if (fullFoundShorter != 0U) {
      std::string lacking;
      for (std::size_t layer = 0; layer < static_cast<std::size_t>(Detector<Det>::nLayers); ++layer) {
        lacking += std::format(" {}: {:.1f}", layer, percent(lackingLayer[layer], fullFoundShorter));
      }
      std::string lengths;
      for (std::size_t n = kSharedToMatch; n < kMaxClusters; ++n) {
        if (shorterLength[n] != 0U) {
          lengths += std::format(" {}: {:.1f}", n, percent(shorterLength[n], fullFoundShorter));
        }
      }
      LOGP(info,
           "    O2's full-length tracks GTRK finds shorter: {}; the layer GTRK's track lacks (% of them, several per track):{}; GTRK's length (%):{}; "
           "a lacking cluster held by another GTRK track {} ({:.1f} %), all free {}",
           fullFoundShorter, lacking, lengths, shorterLackingTaken, percent(shorterLackingTaken, fullFoundShorter), shorterLackingFree);
    }
    if (!states.empty()) {
      const auto bins = [&](const auto& edges, const auto& counts) {
        std::string line;
        for (std::size_t bin = 0; bin + 1U < edges.size(); ++bin) {
          const auto total = counts[bin][0] + counts[bin][1] + counts[bin][2];
          if (total != 0U) {
            line += std::format(" [{:.1f},{:.1f}) {}: {:.1f}/{:.1f}/{:.1f}", edges[bin], edges[bin + 1U], total, percent(counts[bin][0], total),
                                percent(counts[bin][1], total), percent(counts[bin][2], total));
          }
        }
        return line;
      };
      LOGP(info, "    O2's full-length tracks by pT [GeV] (count: % full / % shorter / % missed):{}", bins(kO2PtEdges, fullByPt));
      LOGP(info, "    O2's full-length tracks by |tan lambda| (count: % full / % shorter / % missed):{}", bins(kO2TglEdges, fullByTgl));
    }
    {
      std::vector<std::pair<std::size_t, std::uint32_t>> top;
      top.reserve(shortPatterns.size());
      for (const auto& [mask, count] : shortPatterns) {
        top.emplace_back(count, mask);
      }
      std::ranges::sort(top, std::ranges::greater{});
      std::string line;
      std::size_t shortTotal = 0U;
      for (const auto& [count, mask] : top) {
        shortTotal += count;
      }
      for (std::size_t i = 0; i < std::min<std::size_t>(top.size(), 6U); ++i) {
        std::string layers;
        for (std::size_t layer = 0; layer < static_cast<std::size_t>(Detector<Det>::nLayers); ++layer) {
          if ((top[i].second & (1U << layer)) != 0U) {
            layers += std::to_string(layer);
          }
        }
        line += std::format(" {}: {} ({:.1f})", layers, top[i].first, percent(top[i].first, shortTotal));
      }
      LOGP(info, "    GTRK's tracks short of a layer, by layers (%):{}", line);
    }
    if (unknownClusters != 0U) {
      LOGP(warning, "    {} clusters of O2's tracks are not among the {} measurements", unknownClusters, frame.size());
    }
  }

  /// Each cluster index's measurement row, the maximum for none: one
  /// measurement per cluster, so a cluster index at or past the row count
  /// has none either.
  [[nodiscard]] static std::vector<std::uint32_t> rowsOfClusters(const Frame& frame)
  {
    std::vector<std::uint32_t> rowOfCluster(frame.size(), std::numeric_limits<std::uint32_t>::max());
    for (Frame::Index row = 0U; row < frame.size(); ++row) {
      const auto cluster = frame.measurement(row).clusterIndex;
      if (cluster < rowOfCluster.size()) {
        rowOfCluster[cluster] = row;
      }
    }
    return rowOfCluster;
  }

  /// O2's tracks of this timeframe as our measurement rows: O2's tracks name
  /// O2's cluster indices, which are our measurements' clusterIndex.
  [[nodiscard]] TrackRows o2TrackRows(const Frame& frame) const
  {
    const auto rowOfCluster = rowsOfClusters(frame);
    const auto clusters = mO2Tracks->read(static_cast<long long>(mTimeframe));
    TrackRows theirs;
    theirs.reserve(clusters.size());
    for (const auto& track : clusters) {
      auto& rowsOfTrack = theirs.emplace_back();
      for (const auto cluster : track) {
        if (cluster >= rowOfCluster.size() || rowOfCluster[cluster] == std::numeric_limits<std::uint32_t>::max()) {
          LOGP(fatal, "an O2 track holds cluster {}, which is not among the {} measurements", cluster, frame.size());
        }
        rowsOfTrack.push_back(rowOfCluster[cluster]);
      }
    }
    return theirs;
  }

  /// Our tracks against O2's, \a theirs (`o2TrackRows`), with MC.
  void compareWithO2(const Frame& frame, std::span<const TruthRow> rows, const ParticleTable& particles, const typename TrackRefit<Det>::Tracks& tracks,
                     const TrackRows& theirs)
  {
    const TrackRows ours = trackRows(tracks);
    for (const auto& line : compareTrackSets(rows, particles, ours, theirs, comparisonSettings<Det>(), "GTRK", "O2")) {
      LOGP(info, "{}", line);
    }
    // The rule held to O2's own: its track labels' fake flags.
    const auto flags = mO2Tracks->fakeFlags(static_cast<long long>(mTimeframe));
    const auto fakes = std::ranges::count(flags, true);
    LOGP(info, "  O2's own track labels: {} of {} fake ({:.2f} %)", fakes, flags.size(), percent(static_cast<std::uint64_t>(fakes), flags.size()));
    if (mHits) {
      measurePulls(frame, rows, tracks, theirs);
    }
  }

  /// The species the q/p pulls are split by: electrons, muons, pions, kaons,
  /// protons, and the rest.
  static constexpr std::array<const char*, kSpeciesCount> kSpecies{"electrons", "muons", "pions", "kaons", "protons", "others"};

  /// The crossing of \a chip by \a particle nearest to \a near, from the MC
  /// hits; none if the particle is not the signal's or crossed no such chip.
  [[nodiscard]] const study::HitSegment* crossingOf(std::uint64_t particle, std::uint32_t chip, const std::array<double, 3>& near)
  {
    const auto found = mParticleLabels.find(particle);
    if (found == mParticleLabels.end() || found->second.source != 0) {
      return nullptr;
    }
    const study::HitSegment* best = nullptr;
    double bestDistance = std::numeric_limits<double>::infinity();
    for (const auto& segment : mHits->event(found->second.event)) {
      if (segment.track != found->second.track || segment.chip != chip) {
        continue;
      }
      double distance = 0.;
      for (std::size_t i = 0; i < 3U; ++i) {
        const double d = (0.5 * (segment.start[i] + segment.end[i])) - near[i];
        distance += d * d;
      }
      if (distance < bestDistance) {
        bestDistance = distance;
        best = &segment;
      }
    }
    return best;
  }

  /// The MC hits' crossings of a particle's rows' chips, in layer order, for
  /// trueAnatomy; an empty function without the hits, and an empty list where
  /// a row has no crossing of the particle's.
  [[nodiscard]] TrueCrossings trueCrossings(const Frame& frame)
  {
    if (!mHits) {
      return {};
    }
    return [this, &frame](std::uint64_t particle, std::span<const Frame::Index> rows) {
      std::vector<TrueCrossing> result;
      const auto found = mParticleLabels.find(particle);
      if (found == mParticleLabels.end()) {
        return result;
      }
      const auto* pdg = o2::O2DatabasePDG::Instance()->GetParticle(found->second.pdg);
      const double mass = pdg == nullptr ? 0.13957 : pdg->Mass();
      for (const auto row : rows) {
        const auto m = frame.measurement(row);
        const auto* segment = crossingOf(particle, m.surface, {static_cast<double>(m.x), static_cast<double>(m.y), static_cast<double>(m.z)});
        if (segment == nullptr) {
          return std::vector<TrueCrossing>{};
        }
        result.push_back(
          {.layer = static_cast<std::uint32_t>(mInterface.geometry().sensor(m.surface).layer),
           .row = row,
           .position = {0.5 * (segment->start[0] + segment->end[0]), 0.5 * (segment->start[1] + segment->end[1]), 0.5 * (segment->start[2] + segment->end[2])},
           .momentum = segment->momentum,
           .mass = mass,
           .pdg = found->second.pdg});
      }
      std::ranges::stable_sort(result, {}, &TrueCrossing::layer);
      return result;
    };
  }

  /// A fitted state against its particle's MC crossing of the state's
  /// sensor, the one nearest the measurement at \a row.
  struct StateTruth {
    std::array<double, 5> residual{}; ///< fitted less true: loc0, loc1, phi, theta, q/p
    std::array<double, 5> sigma{};    ///< the fit's
    double pt = 0.;                   ///< the truth's [GeV]
    double momentum = 0.;             ///< the truth's [GeV]
    double fittedPt = 0.;             ///< [GeV]
    double ptSigma = 0.;              ///< of fittedPt, J C J^T over theta and q/p [GeV]
  };

  /// The truth on the plane the state is bound to (TrackingGeometry::planeOf:
  /// a flat sensor's own at its measurement depth, a curved sensor's tangent
  /// plane at the state's angle): where the MC segment's normal coordinate is
  /// zero; its direction and momentum there. None without the segment, for a
  /// neutral particle or a segment along the plane.
  [[nodiscard]] std::optional<StateTruth> stateTruth(const Frame& frame, const BoundTrackState& state, std::uint64_t particle, Frame::Index row)
  {
    const auto m = frame.measurement(row);
    const auto* segment = crossingOf(particle, state.surface, {static_cast<double>(m.x), static_cast<double>(m.y), static_cast<double>(m.z)});
    const double charge = mParticleLabels[particle].charge;
    if (segment == nullptr || exactlyEqual(charge, 0.) || !mInterface.configured()) {
      return std::nullopt;
    }
    const TrackingSurface plane = mInterface.tracking().planeOf(state); // a copy: a curved sensor's is made on the spot
    const auto local = [&](const std::array<double, 3>& global) {
      std::array<double, 3> l{};
      for (std::size_t i = 0; i < 3U; ++i) {
        const double d = global[i] - static_cast<double>(plane.origin[i]);
        l[0] += d * static_cast<double>(plane.axis0[i]);
        l[1] += d * static_cast<double>(plane.axis1[i]);
        l[2] += d * static_cast<double>(plane.normal[i]);
      }
      return l;
    };
    const auto start = local(segment->start);
    const auto end = local(segment->end);
    if (std::abs(start[2] - end[2]) < 1e-9) {
      return std::nullopt;
    }
    const double t = start[2] / (start[2] - end[2]);
    const auto& p = segment->momentum;
    StateTruth truth;
    truth.pt = std::hypot(p[0], p[1]);
    const double momentum = std::hypot(truth.pt, p[2]);
    truth.momentum = momentum;
    truth.residual = {static_cast<double>(state.par[kLoc0]) - (start[0] + (t * (end[0] - start[0]))),
                      static_cast<double>(state.par[kLoc1]) - (start[1] + (t * (end[1] - start[1]))),
                      std::remainder(static_cast<double>(state.par[kPhi]) - std::atan2(p[1], p[0]), 2. * std::numbers::pi),
                      static_cast<double>(state.par[kTheta]) - std::atan2(truth.pt, p[2]), static_cast<double>(state.par[kQOverP]) - (charge / momentum)};
    for (int i = 0; i < 5; ++i) {
      truth.sigma[static_cast<std::size_t>(i)] = std::sqrt(std::max(static_cast<double>(state.covAt(i, i)), 0.));
    }
    // pT = sin(theta) / |q/p|: d/dtheta = cos(theta) / |q/p|, d/d(q/p) = -pT / (q/p).
    const auto theta = static_cast<double>(state.par[kTheta]);
    const auto qOverP = static_cast<double>(state.par[kQOverP]);
    truth.fittedPt = std::sin(theta) / std::abs(qOverP);
    const double dTheta = std::cos(theta) / std::abs(qOverP);
    const double dQOverP = -truth.fittedPt / qOverP;
    const double variance = (dTheta * dTheta * static_cast<double>(state.covAt(kTheta, kTheta))) +
                            (2. * dTheta * dQOverP * static_cast<double>(state.covAt(kTheta, kQOverP))) +
                            (dQOverP * dQOverP * static_cast<double>(state.covAt(kQOverP, kQOverP)));
    truth.ptSigma = std::sqrt(std::max(variance, 0.));
    return truth;
  }

  /// The oracle refits of good track \a track of \a particle
  /// (SpeciesOracle.h) into \a out, each against the particle's
  /// crossing of \a row's sensor, \a entry its resolution entry; none for a
  /// species without an oracle mass.
  void speciesOracle(const Frame& frame, const typename TrackRefit<Det>::Tracks::ConstView& view, std::size_t track, std::uint64_t particle, Frame::Index row,
                     const Propagator& propagator, const QAResolution& entry, std::vector<QAOracle>& out)
  {
    const auto mass = oracleMass(entry.pdg);
    if (!mass.has_value()) {
      return;
    }
    const RowChain<Det> chain{.frame = &frame, .rows = view.clusters[track].index};
    const int count = view.count[track];
    // The momentum lost inside the fitted chain: the crossing of the
    // outermost measurement's sensor.
    float pOuter = std::numeric_limits<float>::quiet_NaN();
    const auto outer = frame.measurement(view.clusters[track].index[count - 1]);
    if (const auto* segment = crossingOf(particle, outer.surface, {static_cast<double>(outer.x), static_cast<double>(outer.y), static_cast<double>(outer.z)});
        segment != nullptr) {
      pOuter = static_cast<float>(std::hypot(segment->momentum[0], segment->momentum[1], segment->momentum[2]));
    }
    for (std::size_t v = 0; v < kOracleVariants; ++v) {
      const auto fit = oracleFit(static_cast<OracleVariant>(v), *mass, propagator, mInterface.tracking(), mInterface.material(), chain, count, view.outwardState[track]);
      QAOracle o{.pdg = entry.pdg,
                 .charge = entry.charge,
                 .origin = entry.origin,
                 .eta = entry.eta,
                 .layer = entry.layer,
                 .kinks = entry.kinks,
                 .clusters = entry.clusters,
                 .productionOffset = entry.productionOffset,
                 .pOuter = pOuter,
                 .variant = static_cast<std::uint8_t>(v),
                 .status = static_cast<std::uint8_t>(fit.status),
                 .chi2 = fit.chi2,
                 .ndf = fit.ndf};
      const auto truth = fit.status == PropagationStatus::kOk ? stateTruth(frame, fit.state, particle, row) : std::nullopt;
      if (truth.has_value()) {
        o.momentum = static_cast<float>(truth->momentum);
        o.ptResidual = static_cast<float>((truth->fittedPt / truth->pt) - 1.);
        for (std::size_t i = 0; i < 5U; ++i) {
          o.residual[i] = static_cast<float>(truth->residual[i]);
          o.sigma[i] = static_cast<float>(truth->sigma[i]);
        }
      } else if (fit.status == PropagationStatus::kOk) {
        o.status = std::numeric_limits<std::uint8_t>::max();
      }
      out.push_back(o);
    }
  }

  /// \a state's closest approach to its particle's MC event vertex into
  /// \a entry: carried to the plane
  /// through the vertex along the beam, normal to its transverse direction
  /// at the perigee (`beamLinePlane`), loc0 there d0 and loc1 less the
  /// vertex's z z0; their sigmas J C J^T, J the analytic Jacobian; and the
  /// scattering between the vertex and \a row, the innermost measurement, as
  /// the beam-line update takes it (`sourceScattering`) but without the
  /// measurement's own layer, whose stave the state is on: theta_k^2 s_k^2
  /// for d0, over sin^2(theta) for z0. Left NaN where the propagation fails.
  /// And \a truth's production point's distance from the vertex.
  void vertexApproach(const Frame& frame, const BoundTrackState& state, std::uint64_t particle, Frame::Index row, const Propagator& propagator,
                      const ParticleTruth& truth, QAResolution& entry) const
  {
    const auto& label = mParticleLabels.at(particle);
    if (!mKinematics || label.source < 0 || label.event < 0) {
      return;
    }
    const auto& header = mKinematics->getMCEventHeader(label.source, label.event);
    entry.productionOffset = static_cast<float>(
      std::hypot(static_cast<double>(truth.vx) - header.GetX(), static_cast<double>(truth.vy) - header.GetY(), static_cast<double>(truth.vz) - header.GetZ()));
    const auto& surface = mInterface.tracking().planeOf(state);
    const auto plane = beamLinePlane(static_cast<Scalar>(header.GetX()), static_cast<Scalar>(header.GetY()), perigeePhi(propagator, state, surface));
    SolenoidJacobian jacobian;
    BoundTrackState predicted;
    const auto status = jacobianAnalytic(propagator, state, surface, plane, PropagationDirection::kBackward, jacobian, &predicted);
    if (status != PropagationStatus::kOk && status != PropagationStatus::kOutOfBounds) {
      return;
    }
    const auto variance = [&](const Scalar* h) {
      double v = 0.;
      for (int i = 0; i < kNBoundParams; ++i) {
        for (int j = 0; j < kNBoundParams; ++j) {
          v += static_cast<double>(h[i]) * static_cast<double>(state.covAt(i, j)) * static_cast<double>(h[j]);
        }
      }
      return v;
    };
    const double sinTheta = std::sin(static_cast<double>(state.theta()));
    const auto layer = mInterface.tracking().surface(state.surface).layer;
    const auto source = pairPointAt(frame.measurements(), row);
    PairCut ownMomentum;
    const Scalar betaPT = betaMomentum(1.F / std::abs(state.qOverP()), propagator.massGeV) * std::abs(std::sin(state.theta()));
    ownMomentum.invBetaMomentum2 = 1.F / (betaPT * betaPT);
    const VertexSeed vertex{.z = static_cast<Scalar>(header.GetZ())};
    // sourceScattering's sum without the source's own layer.
    const Scalar dz = source.z - vertex.z;
    const Scalar slope = dz / source.r;
    const Scalar absDz = std::abs(dz);
    const Scalar invAbsDz = absDz > 0.F ? 1.F / absDz : 0.F;
    const Scalar discFactor = source.r * source.r * source.r * invAbsDz * invAbsDz * invAbsDz;
    Scalar inside = 0.F;
    Scalar term = 0.F;
    for (std::size_t other = 0; other < ScatteringTable<Det>::kLayerCount; ++other) {
      const auto id = static_cast<LayerID>(other);
      if (id != layer &&
          detail::crossedBeforeSource(mInterface.scattering().layers[other], Detector<Det>::isBarrelLayer(id), false, vertex, source, dz, slope, absDz, discFactor, term)) {
        inside += term;
      }
    }
    if constexpr (ScatteringTable<Det>::kPassiveCount > 0) {
      for (const auto& passive : mInterface.scattering().passives) {
        if (detail::crossedBeforeSource(passive.at, passive.barrel, false, vertex, source, dz, slope, absDz, discFactor, term)) {
          inside += term;
        }
      }
    }
    inside *= std::sqrt(1.F + (slope * slope)) * ownMomentum.invBetaMomentum2;
    entry.dca = {predicted.loc0(), static_cast<float>(static_cast<double>(predicted.loc1()) - header.GetZ())};
    entry.dcaSigma = {static_cast<float>(std::sqrt(std::max(variance(jacobian.loc0), 0.))),
                      static_cast<float>(std::sqrt(std::max(variance(jacobian.loc1), 0.)))};
    entry.dcaScattering = {std::sqrt(inside), static_cast<float>(static_cast<double>(std::sqrt(inside)) / std::abs(sinTheta))};
  }

  /// The fitted state at the innermost measurement of every good track, ours
  /// and O2's, against the MC hit on that measurement's chip, each in its
  /// own parameters and covariance, and both as common residuals
  /// (TrackPulls.h). Ours are bound to the sensor plane; O2's are in its
  /// tracking frame, which only the ITS's tracks are.
  void measurePulls(const Frame& frame, std::span<const TruthRow> rows, const typename TrackRefit<Det>::Tracks& tracks, const TrackRows& theirs)
  {
    const std::vector<PullParameter> common{{.name = "across", .unit = "um", .scale = 1e4},
                                            {.name = "along", .unit = "um", .scale = 1e4},
                                            {.name = "phi", .unit = "mrad", .scale = 1e3},
                                            {.name = "dip", .unit = "mrad", .scale = 1e3},
                                            {.name = "pT", .unit = "%", .scale = 1e2}};
    PullTable ours{"GTRK at the innermost measurement, bound parameters",
                   {{.name = "loc0", .unit = "um", .scale = 1e4},
                    {.name = "loc1", .unit = "um", .scale = 1e4},
                    {.name = "phi", .unit = "mrad", .scale = 1e3},
                    {.name = "theta", .unit = "mrad", .scale = 1e3},
                    {.name = "q/p", .unit = "1/GeV", .scale = 1.}}};
    PullTable oursCommon{"GTRK, common residuals", common};
    // The q/p pull and the pT bias split by species, whose mass the fit's
    // energy loss does not know, and by kink.
    const std::vector<PullParameter> qOverPAndPt{{.name = "q/p", .unit = "1/GeV", .scale = 1.}, {.name = "pT", .unit = "%", .scale = 1e2}};
    const auto split = [&](const std::string& who, std::span<const char* const> classes) {
      std::vector<PullTable> tables;
      tables.reserve(classes.size());
      for (const auto* name : classes) {
        tables.emplace_back(fmt::format("{}, {}", who, name), qOverPAndPt);
      }
      return tables;
    };
    auto oursBySpecies = split("GTRK q/p by species", kSpecies);
    constexpr std::array<const char*, 2> kKinkClasses{"without a kink", "with a kink"};
    auto oursByKink = split("GTRK q/p", kKinkClasses);
    const auto view = tracks.view();
    std::vector<std::uint32_t> trackRows;
    for (std::size_t track = 0; track < tracks.size(); ++track) {
      const auto& clusters = view.clusters[track];
      trackRows.assign(std::begin(clusters.index), std::begin(clusters.index) + view.count[track]);
      const auto particle = goodParticle(rows, trackRows);
      const auto& state = view.state[track];
      const auto innermost = std::ranges::find_if(trackRows, [&](std::uint32_t row) { return frame.measurement(row).surface == state.surface; });
      if (particle == kNoTrackParticle || innermost == trackRows.end()) {
        continue;
      }
      const auto truth = stateTruth(frame, state, particle, *innermost);
      if (!truth.has_value()) {
        continue;
      }
      const auto& residual = truth->residual;
      const auto& sigma = truth->sigma;
      const double pt = truth->pt;
      ours.add(pt, residual, sigma);
      const double fittedPt = truth->fittedPt;
      oursCommon.add(pt, std::array<double, 5>{residual[0], residual[1], residual[2], -residual[3], (fittedPt / pt) - 1.}, std::array<double, 5>{});
      const std::array<double, 2> qOverP{residual[4], (fittedPt / pt) - 1.};
      const std::array<double, 2> qOverPSigma{sigma[4], 0.};
      oursBySpecies[speciesOf(mParticleLabels[particle].pdg)].add(pt, qOverP, qOverPSigma);
      oursByKink[view.kinks[track] == 0U ? 0U : 1U].add(pt, qOverP, qOverPSigma);
    }
    std::vector<std::string> lines = ours.report();
    const auto append = [&](const PullTable& table) {
      const auto more = table.report();
      lines.insert(lines.end(), more.begin(), more.end());
    };
    append(oursCommon);
    for (const auto* tables : {&oursBySpecies, &oursByKink}) {
      for (const auto& table : *tables) {
        append(table);
      }
    }

    if constexpr (std::is_same_v<Det, detectors::ITS>) {
      const auto point = [&](Frame::Index row) {
        const auto m = frame.measurement(row);
        return std::array<double, 3>{static_cast<double>(m.x), static_cast<double>(m.y), static_cast<double>(m.z)};
      };
      const auto states = mO2Tracks->innerStates(static_cast<long long>(mTimeframe));
      PullTable o2{"O2 at the innermost cluster, tracking-frame parameters",
                   {{.name = "y", .unit = "um", .scale = 1e4},
                    {.name = "z", .unit = "um", .scale = 1e4},
                    {.name = "snp", .unit = "1e-3", .scale = 1e3},
                    {.name = "tgl", .unit = "1e-3", .scale = 1e3},
                    {.name = "q/pT", .unit = "1/GeV", .scale = 1.}}};
      PullTable o2Common{"O2, common residuals", common};
      auto o2BySpecies = split("O2 q/pT by species", kSpecies);
      for (std::size_t track = 0; track < std::min(theirs.size(), states.size()); ++track) {
        const auto particle = goodParticle(rows, theirs[track]);
        const auto innermost = std::ranges::min_element(theirs[track], {}, [&](std::uint32_t row) { return rows[row].layer; });
        if (particle == kNoTrackParticle || innermost == theirs[track].end()) {
          continue;
        }
        const auto chip = frame.measurement(*innermost).surface;
        const auto* segment = crossingOf(particle, chip, point(*innermost));
        const double charge = mParticleLabels[particle].charge;
        if (segment == nullptr || exactlyEqual(charge, 0.)) {
          continue;
        }
        const auto& state = states[track];
        const double c = std::cos(static_cast<double>(state.getAlpha()));
        const double s = std::sin(static_cast<double>(state.getAlpha()));
        const auto rotate = [&](const std::array<double, 3>& v) { return std::array<double, 3>{(c * v[0]) + (s * v[1]), (c * v[1]) - (s * v[0]), v[2]}; };
        const auto a = rotate(segment->start);
        const auto b = rotate(segment->end);
        const auto p = rotate(segment->momentum);
        if (std::abs(b[0] - a[0]) < 1e-9) {
          continue;
        }
        const double t = (static_cast<double>(state.getX()) - a[0]) / (b[0] - a[0]);
        const double pt = std::hypot(p[0], p[1]);
        const std::array<double, 5> truth{a[1] + (t * (b[1] - a[1])), a[2] + (t * (b[2] - a[2])), p[1] / pt, p[2] / pt, charge / pt};
        std::array<double, 5> residual{};
        std::array<double, 5> sigma{};
        for (std::size_t i = 0; i < 5U; ++i) {
          residual[i] = static_cast<double>(state.getParam(i)) - truth[i];
          sigma[i] = std::sqrt(std::max(static_cast<double>(state.getDiagError2(i)), 0.));
        }
        o2.add(pt, residual, sigma);
        o2Common.add(pt,
                     std::array<double, 5>{residual[0], residual[1], std::asin(std::clamp(static_cast<double>(state.getParam(2)), -1., 1.)) - std::asin(truth[2]),
                                           std::atan(static_cast<double>(state.getParam(3))) - std::atan(truth[3]), (std::abs(truth[4]) / std::abs(static_cast<double>(state.getParam(4)))) - 1.},
                     std::array<double, 5>{});
        o2BySpecies[speciesOf(mParticleLabels[particle].pdg)].add(pt, std::array<double, 2>{residual[4], (std::abs(truth[4]) / std::abs(static_cast<double>(state.getParam(4)))) - 1.},
                                                                  std::array<double, 2>{sigma[4], 0.});
      }
      append(o2);
      append(o2Common);
      for (const auto& table : o2BySpecies) {
        append(table);
      }
    }
    LOGP(info, "pulls against the MC hits");
    for (const auto& line : lines) {
      LOGP(info, "{}", line);
    }
  }

  std::unique_ptr<study::TrackFile> mO2Tracks; ///< O2's tracks of these clusters; none unless --o2-tracks
  std::unique_ptr<study::HitFile> mHits;       ///< O2's MC hits, with MC; none if the file is not there
  /// Every particle the comparison saw, as its MC label and charge.
  struct ParticleLabel {
    int event = -1;
    int track = -1;
    int source = -1;
    double charge = 0.;
    int pdg = 0;
    TruthHelix helix{}; ///< at production
    // NOLINTNEXTLINE(readability-redundant-member-init): GCC's -Wmissing-field-initializers wants it
    std::optional<TruthHelix> mother{}; ///< the mother's at its own production
  };
  std::unordered_map<std::uint64_t, ParticleLabel> mParticleLabels;
  bool mUseMC = true;
  bool mSpeciesOracle = false;                               ///< --species-oracle: the `oracle` QA tree
  bool mTruthReports = true;                                 ///< the stages' truth, the true roads' chi2 and the fakes' anatomy in the log; on unless --no-truth-reports, data runs too
  bool mWithVertices = true;                                 ///< the vertex list is an input
  TrackerProducts mTrackerProducts = TrackerProducts::kNone; ///< the tracking device's products, an input: compared, or digested
  workflow::TrackingInterface<Det> mInterface;               ///< the tracking: geometry, material, tracker, vertexer (TrackingInterface.h)
  std::unique_ptr<o2::steer::DigitizationContext> mContext;
  std::unique_ptr<o2::steer::MCKinematicsReader> mKinematics;
  std::shared_ptr<workflow::Conditions<Det, true>> mConditions; ///< the run's field, geometry, MeanVertex and material LUT (Conditions.h)
  std::ofstream mTruthOutput;
  std::unique_ptr<o2::utils::TreeStreamRedirector> mQAOutput; ///< --qa-output's file; none without it
  /// The tracker's observers over one timeframe: the open iteration's stage
  /// truth, made at its first slice, and each finished iteration's summary;
  /// with MC its stage funnel.
  struct RunObservers {
    std::optional<StageTruth> stage;
    std::vector<TrackletRunSummary> summaries;
    std::function<void(const Controller&, std::size_t, std::optional<StageTruth>&)> openStage; ///< set per timeframe
    std::optional<Funnel> funnel;
  };
  RunObservers mObservers;
  std::string mDumpInput;          ///< where the dumped timeframe's tracking input goes; empty: nowhere
  std::size_t mDumpTimeframe = 0U; ///< which timeframe, counting from 0
  bool mDumpOnly = false;          ///< the dumped timeframe is not tracked
  std::size_t mTimeframe = 0U;
};

/// The detector's MC hits in an O2 production's timeframe directory.
[[nodiscard]] constexpr const char* hitFileName()
{
  if constexpr (std::is_same_v<Det, detectors::MFT>) {
    return "sgn_HitsMFT.root";
  } else if constexpr (std::is_same_v<Det, detectors::TRKFT3>) {
    return "sgn_HitsTRK.root";
  } else {
    return "sgn_HitsITS.root";
  }
}

} // namespace

template <>
DataProcessorSpec getMeasurementSinkSpec<Det>(bool useMC, bool withVertices, TrackerProducts products)
{
  auto inputs = workflow::measurementInputs<Det>(withVertices);
  if (useMC) {
    inputs.emplace_back("labels", workflow::gDataOrigin, workflow::description<Det>(workflow::kClustersMCTruth), 0, Lifetime::Timeframe);
  }
  if (products != TrackerProducts::kNone) {
    inputs.emplace_back("deviceTracks", workflow::gDataOrigin, workflow::description<Det>(workflow::kTracks), 0, Lifetime::Timeframe);
    inputs.emplace_back("deviceVertices", workflow::gDataOrigin, workflow::description<Det>(workflow::kVertexSeeds), 0, Lifetime::Timeframe);
    inputs.emplace_back("deviceContributors", workflow::gDataOrigin, workflow::description<Det>(workflow::kVertexContributors), 0, Lifetime::Timeframe);
  }
  auto conditions = std::make_shared<workflow::Conditions<Det, true>>(inputs);

  DataProcessorSpec spec{
    .name = workflow::deviceName<Det>("measurement-sink"),
    .inputs = inputs,
    .outputs = Outputs{},
    .algorithm = AlgorithmSpec{adaptFromTask<MeasurementSinkTask>(useMC, withVertices, products, conditions)},
    .options = Options{
      {"collision-context", VariantType::String, "collisioncontext.root", {"collision context and MC-event mapping"}},
      {"input-dir", VariantType::String, "none", {"input directory"}},
      {"tracklet-truth-output", VariantType::String, "gtrk_tracklet_truth.csv", {"per-link pT and cot(theta) truth breakdown; empty disables"}},
      {"dump-tracking-input", VariantType::String, "", {"write one timeframe's tracking input here, for o2-gtrk-track-bench"}},
      {"dump-timeframe", VariantType::Int, 0, {"which timeframe --dump-tracking-input writes, counting from 0"}},
      {"dump-only", VariantType::Bool, false, {"with --dump-tracking-input, the dumped timeframe is not tracked"}},
      {"qa-output", VariantType::String, "", {"the MC truth QA as ROOT trees; empty disables"}},
      {"species-oracle",
       VariantType::Bool,
       false,
       {"with --qa-output, refit every good track of e, mu, pi, K, p with its true mass: the oracle tree"}},
      {"no-truth-reports", VariantType::Bool, false, {"with MC, no truth reports in the log (stages, true roads' chi2, fake anatomy): for the QA alone"}},
      {"o2-tracks", VariantType::String, "", {"O2's tracks of the same clusters (o2trac_its.root, mfttracks.root), compared under one set of truth rules"}},
      {"mc-hits", VariantType::String, hitFileName(), {"O2's MC hits, for the pulls, the true kinks and the QA's outward legs"}}}};
  if constexpr (std::is_same_v<Det, detectors::TRKFT3>) {
    spec.options.push_back({"mc-ft3-hits", VariantType::String, "sgn_HitsFT3.root", {"O2's MC hits of FT3's discs, with --mc-hits TRK's"}});
  }
  return spec;
}

} // namespace o2::gtrk::study
