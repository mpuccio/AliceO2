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

/// \file StageFunnel.h
/// \brief Where the findable particles are lost, stage by stage, in the
///        tracker's own run.
///
/// A particle above the pT floor with measurements on a chain of `minLayers`
/// layers is findable (TrackComparison.h). The funnel observes every
/// slice of every iteration and records, per findable particle, the longest
/// chain in clusters that is entirely its own -- every measurement carries it
/// among its labels -- at each stage: kept tracklets joined end to end, kept
/// triplets joined by a shared tracklet, kept triplets joined by kept
/// neighbours, the roads the walk gives, and the roads that passed the fit
/// and the pT floor (the selection's records). The longest of any iteration
/// counts. A particle without a good track is lost at the first stage whose
/// chain is shorter than `minLayers`; with a record of its own, at the
/// selection.
///
/// For the tracklets, the reason: each adjacent-layer pair of the particle's
/// own measurements that no iteration kept is tested again with the
/// iteration's parameters, in the order the stage tests -- the link not taken,
/// a measurement held by an earlier iteration's track, the source interval,
/// time, the turning point, the curvature (the pair's circle through the beam
/// line below the iteration's pT_min), the longitudinal test -- and what passes
/// them all was lost elsewhere: the grid's window, or the slicing. An
/// iteration from the vertex list records no reasons.
///
/// A particle is seeded in an iteration if one of that iteration's roads,
/// before the road fit, is entirely its own.
///
/// Beyond the stages, for the bench: the length funnel (`lengthReport`, the
/// found against their ceiling) and every record's chi2 (`recordChi2`, for
/// `--fakes-out`), whose inputs are kept only when asked (`keepReportInputs`).

#ifndef GTRK_STUDY_STAGEFUNNEL_H_
#define GTRK_STUDY_STAGEFUNNEL_H_

#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LinkHoles.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackletStage.h"
#include "GenericTracking/Study/Statistics.h"
#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/Validation.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace o2::gtrk::study
{

/// The stages of the funnel, in order; kFound last.
enum class FunnelStage : std::uint8_t { kTracklets,
                                        kTriplets,
                                        kNeighbours,
                                        kRoads,
                                        kFit,
                                        kSelection,
                                        kFound,
                                        kCount };
inline constexpr std::size_t kFunnelStages = static_cast<std::size_t>(FunnelStage::kCount);
inline constexpr std::array<const char*, kFunnelStages> kFunnelStageNames{"tracklets", "triplets", "neighbours", "roads", "road fit", "selection", "found"};

/// Why an adjacent-layer pair of a particle was not kept, in the order the
/// tracklet stage tests.
/// The curvature test is split by the particle's pT: below the iteration's
/// pT_min the loss is the iteration's choice; above it, the pair's circle
/// through the beam line is not the particle's.
enum class PairLoss : std::uint8_t { kLinkNotTaken,
                                     kUsed,
                                     kSourceInterval,
                                     kTime,
                                     kTurning,
                                     kCurvatureBelowPtMin,
                                     kCurvature,
                                     kLongitudinal,
                                     kWindow,
                                     kHole,
                                     kElsewhere,
                                     kCount };
inline constexpr std::size_t kPairLosses = static_cast<std::size_t>(PairLoss::kCount);
inline constexpr std::array<const char*, kPairLosses> kPairLossNames{"link not taken", "used", "source interval", "time", "turning point", "curvature below pT_min", "curvature above pT_min", "longitudinal", "outside the query's cells", "crosses another layer", "window or slicing"};

/// The losses the pair cut itself gives (pairCutLoss), in its order.
inline constexpr std::array<PairLoss, 8> kPairCutLosses{PairLoss::kSourceInterval, PairLoss::kTime, PairLoss::kTurning, PairLoss::kCurvature, PairLoss::kLongitudinal, PairLoss::kWindow, PairLoss::kHole, PairLoss::kElsewhere};

/// The pair cut on a pair of rows on \a link from \a vertex: the first test the
/// tracklet stage applies that it fails -- the source interval, time, the
/// turning point, the curvature, the longitudinal test -- or kElsewhere if it
/// passes them all; and the spatial terms as far as they were evaluated. In a
/// displaced iteration the curvature is the displaced test's transverse part
/// (its own turning point included) and the longitudinal its z; no terms.
/// Includes the query's cells and direct-pair rule; the stage's used,
/// shape and bound tests are not replayed.
struct PairCutVerdict {
  PairLoss loss = PairLoss::kElsewhere;
  PairTerms terms{};
  bool emptyDirectRange = false; ///< no direct part of the query, rather than a pair crossing a hole
};

template <DetectorLayout Det>
[[nodiscard]] PairCutVerdict pairCutLoss(const TimeFrame<Det>& frame, const IterationParameters<Det>& parameters, const VertexSeed& vertex, std::size_t link, std::uint32_t sourceRow, std::uint32_t targetRow)
{
  PairCutVerdict verdict;
  const auto sourceMeasurement = frame.measurement(sourceRow);
  const auto targetMeasurement = frame.measurement(targetRow);
  const auto source = pairPoint(sourceMeasurement);
  const auto target = pairPoint(targetMeasurement);
  const auto from = Det::layerLinks[link].from;
  const Scalar second = secondCoordinate<Det>(from, source.r, source.z);
  if (second < parameters.sourceMin[link] || second > parameters.sourceMax[link]) {
    verdict.loss = PairLoss::kSourceInterval;
    return verdict;
  }
  if (!timeBox(sourceMeasurement).isCompatible(timeBox(targetMeasurement))) {
    verdict.loss = PairLoss::kTime;
    return verdict;
  }
  const auto to = Det::layerLinks[link].to;
  // Only this target's membership matters: its width gives the same spatial
  // cells and suffices for the temporal extent of a compatible pair.
  auto window = candidateWindow(parameters.windows[link], timeBox(sourceMeasurement), sourceMeasurement.phi, source.r, source.z, timeBox(targetMeasurement).width);
  const auto holes = parameters.holeLayers[link];
  if (holes != 0U) {
    const auto& span = parameters.layerSpans[to];
    const auto direct = directPairRange(holes, parameters.barrelLayers, parameters.layerSpans, to, source.r, source.z, span.low, span.high);
    const Scalar low = std::max(direct.low, window.secondCentre - window.secondHalfWidth);
    const Scalar high = std::min(direct.high, window.secondCentre + window.secondHalfWidth);
    if (!(low <= high)) {
      verdict.loss = PairLoss::kHole;
      verdict.emptyDirectRange = true;
      return verdict;
    }
    window.secondCentre = 0.5F * (low + high);
    window.secondHalfWidth = 0.5F * (high - low);
  }
  const auto& grid = frame.grid();
  const auto cells = grid.cellSpan(to, window, grid.geometry().allBlocks());
  const Scalar phi = targetMeasurement.phi < 0.F ? targetMeasurement.phi + (2 * std::numbers::pi_v<Scalar>) : targetMeasurement.phi;
  const auto phiBin = grid.geometry().phiBinOf(to, phi);
  const auto secondBin = grid.geometry().secondBinOf(to, secondCoordinate<Det>(to, target.r, target.z));
  const bool inPhi = (phiBin >= cells.firstPhi && phiBin < cells.phiEnd) || phiBin < cells.wrapEnd;
  const bool inSecond = secondBin >= cells.firstSecond && secondBin <= cells.lastSecond;
  if (cells.blocks.empty() || !inPhi || !inSecond) {
    verdict.loss = PairLoss::kWindow;
    return verdict;
  }
  const auto directPair = [&] {
    if (verdict.loss == PairLoss::kElsewhere && holes != 0U && segmentLeavesHole(holes, parameters.barrelLayers, parameters.layerSpans, source.r, source.z, target.r, target.z)) {
      verdict.loss = PairLoss::kHole;
    }
    return verdict;
  };
  const Scalar scatter = sourceScattering<Det>(parameters.scattering, parameters.pairCut, vertex, from, source);
  const Scalar between = parameters.passiveBetween[link] ? betweenScattering<Det>(parameters.scattering, parameters.pairCut, vertex, source, target, parameters.betweenPassives[link]) : 0.F;
  if (parameters.impactMax > 0.F) {
    // A displaced iteration's test (TrackletStage.h, `displacedPairOf`): its
    // transverse part as the curvature, its z as the longitudinal test.
    const auto result = displacedPairOf(parameters, vertex, source, scatter, target, between);
    if (!result.transverse) {
      verdict.loss = PairLoss::kCurvature;
    } else if (!result.longitudinal) {
      verdict.loss = PairLoss::kLongitudinal;
    }
    return directPair();
  }
  verdict.terms = transversePairTerms(parameters.pairCut, vertex, source, scatter, target, between);
  if (!verdict.terms.beforeTurn) {
    verdict.loss = PairLoss::kTurning;
  } else if (!transverseCompatible(parameters.pairCut, verdict.terms)) {
    verdict.loss = PairLoss::kCurvature;
  } else {
    addLongitudinalPairTerms(verdict.terms, parameters.pairCut, vertex, source, scatter, target, between);
    if (!longitudinalCompatible(parameters.pairCut, verdict.terms)) {
      verdict.loss = PairLoss::kLongitudinal;
    }
  }
  return directPair();
}

/// Why a particle with a chain of kept triplets and neighbours of its own got
/// no road of its own: its outermost triplet below the iteration's start
/// level, raised above its own chain's level by a predecessor not its own (the
/// walk descends one level at a time) -- inside its chain, the walk leaving
/// it for the other's; below it, no triplet of its own before its innermost,
/// so the walk descends through the whole chain and on into the other's; or
/// with a triplet of its own before it but no neighbour to it (the neighbour
/// cut, or ordered generation linking only the highest level that has one) --
/// continued outward by a successor not its own (so not a start), or ending
/// off the iteration's end layers.
enum class RoadLoss : std::uint8_t { kNoChain,
                                     kStartLevel,
                                     kRaised,
                                     kRaisedBelow,
                                     kRaisedUnlinked,
                                     kContinued,
                                     kEndLayer,
                                     kOther,
                                     kCount };
inline constexpr std::size_t kRoadLosses = static_cast<std::size_t>(RoadLoss::kCount);
inline constexpr std::array<const char*, kRoadLosses> kRoadLossNames{"no chain this iteration", "below the start level", "raised inside its chain", "raised below its chain", "raised, own predecessor not linked", "continued by another's successor", "off the end layers", "other"};

/// Why a particle's longest road of its own was no record: its boxes share
/// no time, the fit's chi2, its propagation, the iteration's pT floor for the
/// road's length, the beam line; other if the road fit stage's tests take it
/// (`roadFitRefusal`).
enum class FitLoss : std::uint8_t { kNoRoad,
                                    kChi2,
                                    kPropagation,
                                    kPtFloor,
                                    kBeamLine,
                                    kOther,
                                    kTime,
                                    kCount };
inline constexpr std::size_t kFitLosses = static_cast<std::size_t>(FitLoss::kCount);
inline constexpr std::array<const char*, kFitLosses> kFitLossNames{"no road this iteration", "chi2", "propagation", "pT floor by length", "beam line", "other", "time"};

/// What a missed particle's fake -- a track it labels, with a measurement not
/// its own -- got wrong: one measurement, the innermost; one elsewhere; more.
enum class FakeKind : std::uint8_t { kNone,
                                     kInnermost,
                                     kOtherOne,
                                     kMore,
                                     kCount };
inline constexpr std::size_t kFakeKinds = static_cast<std::size_t>(FakeKind::kCount);
inline constexpr std::array<const char*, kFakeKinds> kFakeKindNames{"no fake of its own", "one wrong, the innermost", "one wrong elsewhere", "more wrong"};

/// The first stage a particle's chains fall short of \a minimum at, from the
/// longest chain of each stage; kSelection if all reach it, kFound if found.
[[nodiscard]] inline FunnelStage firstLoss(const std::array<std::uint32_t, 5>& longest, std::uint32_t minimum, bool found)
{
  if (found) {
    return FunnelStage::kFound;
  }
  for (std::size_t stage = 0; stage < longest.size(); ++stage) {
    if (longest[stage] < minimum) {
      return static_cast<FunnelStage>(stage);
    }
  }
  return FunnelStage::kSelection;
}

/// The longest chain in nodes ending at each node of \a edges that one ends
/// at, each edge from a node to one on a later layer; \a layerOf orders them.
template <typename LayerOf>
[[nodiscard]] std::unordered_map<std::uint32_t, std::uint32_t> chainLengths(std::vector<std::pair<std::uint32_t, std::uint32_t>> edges, LayerOf layerOf)
{
  std::ranges::sort(edges, {}, [&](const auto& edge) { return layerOf(edge.first); });
  std::unordered_map<std::uint32_t, std::uint32_t> length;
  for (const auto& [from, to] : edges) {
    const auto before = length.find(from);
    const std::uint32_t reached = (before == length.end() ? 1U : before->second) + 1U;
    auto& at = length[to];
    at = std::max(at, reached);
  }
  return length;
}

/// The longest chain in nodes of \a edges; 0 without any.
template <typename LayerOf>
[[nodiscard]] std::uint32_t longestChain(std::vector<std::pair<std::uint32_t, std::uint32_t>> edges, LayerOf layerOf)
{
  std::uint32_t result = 0U;
  for (const auto& [node, length] : chainLengths(std::move(edges), layerOf)) {
    result = std::max(result, length);
  }
  return result;
}

/// A pair of rows as one key, the source's in the upper half.
[[nodiscard]] inline std::uint64_t pairKey(std::uint32_t source, std::uint32_t target) { return (static_cast<std::uint64_t>(source) << 32U) | static_cast<std::uint64_t>(target); }

template <DetectorLayout Det, typename Controller>
class StageFunnel
{
 public:
  using Frame = TimeFrame<Det>;
  using SliceObservation = typename Controller::SliceObservation;
  using IterationObservation = typename Controller::IterationObservation;

  /// \a minPtGeV: each iteration's pT_min. Throws std::invalid_argument for
  /// more iterations than a seeded mask holds.
  StageFunnel(const Frame& frame, std::span<const TruthRow> rows, const ParticleTable& truth, const ComparisonSettings& settings, const VertexSeed& beamLine, std::vector<double> minPtGeV)
    : mRows{rows},
      mSettings{settings},
      mBeamLine{beamLine},
      mFrame{&frame},
      mMinPtGeV{std::move(minPtGeV)}
  {
    if (mMinPtGeV.size() > static_cast<std::size_t>(std::numeric_limits<std::uint64_t>::digits)) {
      throw std::invalid_argument(std::format("stage funnel: {} iterations, a seeded mask holds {}", mMinPtGeV.size(), std::numeric_limits<std::uint64_t>::digits));
    }
    for (const auto& [id, bits] : particleLayers(rows, settings)) {
      const auto found = truth.find(id);
      if (found == truth.end() || !isFindable(found->second.pt, bits, settings)) {
        continue;
      }
      mIndex[id] = static_cast<std::uint32_t>(mParticles.size());
      mParticles.push_back({
        .id = id,
        .pt = found->second.pt,
        .species = std::abs(found->second.pdg),
        .group = groupOf(found->second),
        .bin = productionBin(found->second, settings),
        .rows = {},
        .longestBy = {},
        .reasons = {},
        .roadLosses = {},
        .fitLosses = {},
        .recordChi2 = {},
        .pairsKept = 0U,
        .lostPairs = {},
        .fitRoadLayers = {},
        .raisedFits = {},
      });
    }
    // Each row's findable particles, and each particle's rows.
    mRowOffsets.assign(rows.size() + 1U, 0U);
    for (std::size_t row = 0; row < rows.size(); ++row) {
      mRowOffsets[row + 1U] = mRowOffsets[row];
      for (std::uint8_t i = 0; i < rows[row].count; ++i) {
        const auto found = mIndex.find(rows[row].particles[i]);
        if (found != mIndex.end()) {
          mRowParticles.push_back(found->second);
          mParticles[found->second].rows.push_back(static_cast<std::uint32_t>(row));
          ++mRowOffsets[row + 1U];
        }
      }
    }
  }

  void observeSlice(const Controller& tracker, const SliceObservation& slice)
  {
    beginIteration();
    const auto& parameters = tracker.iteration(slice.iteration);
    const auto& result = *slice.result;
    const auto tracklets = result.tracklets->view();
    const auto triplets = result.triplets->view();
    for (std::size_t i = 0; i < result.tracklets->size(); ++i) {
      const std::array<std::uint32_t, 2> ends{tracklets.source[i], tracklets.target[i]};
      const bool owned = forEachOwner(ends, [&](std::uint32_t particle) { mPairs[particle].emplace_back(ends[0], ends[1]); });
      if (owned) {
        mKept.insert(pairKey(ends[0], ends[1]));
      }
    }
    for (std::size_t i = 0; i < result.triplets->size(); ++i) {
      const auto first = triplets.first[i];
      const std::array<std::uint32_t, 3> rows{tracklets.source[first], tracklets.target[first], tracklets.target[triplets.second[i]]};
      // The road counts of the slice that owns it: its first measurement in the core.
      const auto core = slice.slice->core[parameters.linkSources[tracklets.link[first]]];
      const bool owned = rows[0] >= core.begin && rows[0] < core.end;
      forEachOwner(rows, [&](std::uint32_t particle) {
        mTriplets[particle].emplace_back(pairKey(rows[0], rows[1]), pairKey(rows[1], rows[2]));
        if (owned) {
          mTripletInfo[particle][{pairKey(rows[0], rows[1]), rows[2]}] = {
            .level = triplets.level[i],
            .successorAbove = slice.roads->successorAbove()[i] != 0U,
            .endLayer = parameters.linkTargets[tracklets.link[triplets.second[i]]],
          };
        }
      });
    }
    const auto from = result.neighbours->view().from;
    const auto offsets = result.neighbours->offsets();
    for (std::size_t to = 0; to + 1U < offsets.size(); ++to) {
      const auto shared = triplets.first[to];
      const auto b = tracklets.source[shared];
      const auto c = tracklets.target[shared];
      const auto d = tracklets.target[triplets.second[to]];
      for (auto edge = offsets[to]; edge < offsets[to + 1U]; ++edge) {
        const std::array<std::uint32_t, 4> rows{tracklets.source[triplets.first[from[edge]]], b, c, d};
        forEachOwner(rows, [&](std::uint32_t particle) { mNeighbours[particle].push_back(rows); });
      }
    }
    constexpr std::size_t kMax = maxChainLayers<Det>;
    const auto walk = slice.roads->walk(*result.triplets, *result.neighbours);
    for (std::uint32_t road = 0; road < slice.roads->roadCount(); ++road) {
      std::array<std::uint32_t, kMax - 2U> chain{};
      const auto length = walk(road, chain.data());
      std::array<std::uint32_t, kMax> rows{};
      rows[0] = tracklets.source[triplets.first[chain[0]]];
      rows[1] = tracklets.target[triplets.first[chain[0]]];
      for (std::uint32_t i = 0; i < length; ++i) {
        rows[i + 2U] = tracklets.target[triplets.second[chain[i]]];
      }
      const auto count = length + 2U;
      forEachOwner(std::span<const std::uint32_t>{rows.data(), count}, [&](std::uint32_t particle) {
        raise(particle, slice.iteration, FunnelStage::kRoads, count);
        mParticles[particle].seededBy |= std::uint64_t{1} << slice.iteration;
        auto& own = mOwnRoad[particle];
        if (count > own.size()) {
          own.assign(rows.begin(), rows.begin() + count);
        }
      });
    }
    const auto records = slice.records->view();
    for (std::size_t record = 0; record < slice.records->size(); ++record) {
      const auto count = static_cast<std::uint32_t>(records.count[record]);
      const std::span<const std::uint32_t> rows{records.clusters[record].index, count};
      const auto chi2 = static_cast<float>(records.chi2[record]);
      if (mKeepRecordChi2s) {
        const auto [kept, inserted] = mRecordChi2.try_emplace(recordKey(slice.iteration, rows), chi2);
        if (!inserted) {
          kept->second = std::min(kept->second, chi2);
        }
      }
      forEachOwner(rows, [&](std::uint32_t particle) {
        auto& p = mParticles[particle];
        const auto before = slice.iteration < p.longestBy.size() ? p.longestBy[slice.iteration][static_cast<std::size_t>(FunnelStage::kFit)] : 0U;
        raise(particle, slice.iteration, FunnelStage::kFit, count);
        p.recordChi2.resize(std::max(p.recordChi2.size(), slice.iteration + 1U), -1.F);
        auto& best = p.recordChi2[slice.iteration];
        if (count > before || best < 0.F) {
          best = chi2;
        } else if (count == before) {
          best = std::min(best, chi2);
        }
        auto& own = mOwnRecord[particle];
        own = std::max(own, count);
      });
    }
  }

  /// The end of an iteration: its chains, the reasons for its lost pairs, and
  /// the used mask the next one reads.
  void observeIteration(const Controller& tracker, const IterationObservation& iteration)
  {
    beginIteration();
    const auto layerOf = [&](std::uint32_t row) { return mRows[row].layer; };
    if (mKeepOwnKept) {
      for (const auto& [particle, edges] : mPairs) {
        auto& kept = mOwnKept[particle];
        for (const auto& [a, b] : edges) {
          kept.insert(pairKey(a, b));
        }
      }
    }
    for (auto& [particle, edges] : mPairs) {
      raise(particle, iteration.iteration, FunnelStage::kTracklets, longestChain(std::move(edges), layerOf));
    }
    for (auto& [particle, edges] : mTriplets) {
      // A triplet is an edge between its two tracklets; one more measurement each.
      std::vector<std::pair<std::uint32_t, std::uint32_t>> compact;
      std::unordered_map<std::uint64_t, std::uint32_t> ids;
      std::vector<LayerID> layers;
      const auto id = [&](std::uint64_t key) {
        const auto [at, inserted] = ids.try_emplace(key, static_cast<std::uint32_t>(ids.size()));
        if (inserted) {
          layers.push_back(mRows[static_cast<std::uint32_t>(key >> 32U)].layer);
        }
        return at->second;
      };
      for (const auto& [first, second] : edges) {
        compact.emplace_back(id(first), id(second));
      }
      raise(particle, iteration.iteration, FunnelStage::kTriplets, longestChain(std::move(compact), [&](std::uint32_t node) { return layers[node]; }) + 1U);
    }
    for (auto& [particle, chains] : mNeighbours) {
      std::vector<std::pair<std::uint32_t, std::uint32_t>> compact;
      std::map<std::pair<std::uint64_t, std::uint32_t>, std::uint32_t> ids;
      std::vector<LayerID> layers;
      const auto id = [&](std::uint32_t a, std::uint32_t b, std::uint32_t c) {
        const auto [at, inserted] = ids.try_emplace({pairKey(a, b), c}, static_cast<std::uint32_t>(ids.size()));
        if (inserted) {
          layers.push_back(mRows[a].layer);
        }
        return at->second;
      };
      for (const auto& rows : chains) {
        compact.emplace_back(id(rows[0], rows[1], rows[2]), id(rows[1], rows[2], rows[3]));
      }
      // The chain's outermost triplet: the longest chain's end.
      std::uint32_t longest = 0U;
      std::uint32_t end = 0U;
      for (const auto& [node, length] : chainLengths(std::move(compact), [&](std::uint32_t node) { return layers[node]; })) {
        if (length > longest) {
          longest = length;
          end = node;
        }
      }
      raise(particle, iteration.iteration, FunnelStage::kNeighbours, longest + 2U);
      if (longest + 2U >= static_cast<std::uint32_t>(mSettings.minLayers)) {
        const auto key = std::ranges::find_if(ids, [&](const auto& entry) { return entry.second == end; })->first;
        roadLoss(particle, iteration.iteration, tracker.iteration(iteration.iteration), longest, key, chains);
      }
    }
    lostPairReasons(tracker.iteration(iteration.iteration), iteration.iteration);
    fitLosses(tracker.iteration(iteration.iteration), iteration.iteration);
    raisedFitLosses(tracker.iteration(iteration.iteration), iteration.iteration);
    mIterations = std::max(mIterations, iteration.iteration + 1U);
    mUsed.assign(iteration.used.begin(), iteration.used.end());
    mPairs.clear();
    mTriplets.clear();
    mTripletInfo.clear();
    mRaisedChain.clear();
    mNeighbours.clear();
    mKept.clear();
    mOwnRoad.clear();
    mOwnRecord.clear();
  }

  /// What only the bench's reports read, kept from here on (none by
  /// default): each particle's own kept pairs, for `lengthReport`, and with
  /// \a recordChi2s every record's chi2, for `recordChi2`.
  void keepReportInputs(bool recordChi2s)
  {
    mKeepOwnKept = true;
    mKeepRecordChi2s = recordChi2s;
  }

  /// The road fit to refit the roads that were no record with: the run's.
  void setRoadFit(const TrackingGeometry<Det>& tracking, const o2::base::MatLayerCylSet& material, const Propagator& propagator)
  {
    mTracking = &tracking;
    mMaterial = &material;
    mPropagator = propagator;
  }

  /// What one iteration did with a findable particle: its longest chain of its
  /// own at each stage (FunnelStage, tracklets to the road fit), its
  /// adjacent-layer pairs the iteration did not keep by the first test each
  /// fails, if its chain of neighbours reached the minimum without a road,
  /// why, and if its road was no record, why.
  struct IterationStages {
    std::array<std::uint32_t, 5> longest{};
    std::array<std::uint32_t, kPairLosses> lostPairs{};
    RoadLoss roadLoss = RoadLoss::kNoChain;
    FitLoss fitLoss = FitLoss::kNoRoad;
    float recordChi2 = -1.F; ///< the best chi2 of its longest records; -1 without one
  };

  /// Particle \a id's in \a iteration; none for a particle not findable.
  [[nodiscard]] std::optional<IterationStages> stages(std::uint64_t id, std::size_t iteration) const
  {
    const auto found = mIndex.find(id);
    if (found == mIndex.end()) {
      return std::nullopt;
    }
    const auto& particle = mParticles[found->second];
    IterationStages result;
    if (iteration < particle.longestBy.size()) {
      result.longest = particle.longestBy[iteration];
    }
    if (iteration < particle.reasons.size()) {
      result.lostPairs = particle.reasons[iteration];
    }
    if (iteration < particle.roadLosses.size()) {
      result.roadLoss = particle.roadLosses[iteration];
    }
    if (iteration < particle.fitLosses.size()) {
      result.fitLoss = particle.fitLosses[iteration];
    }
    if (iteration < particle.recordChi2.size()) {
      result.recordChi2 = particle.recordChi2[iteration];
    }
    return result;
  }

  /// Particle \a id's verdict over the whole run: the first stage without a chain of its own as long as the minimum,
  /// as if not found; of its adjacent-layer pairs no iteration kept, the test
  /// most fail first over the iterations whose pT_min it clears; why its chain
  /// of neighbours gave no road, and why its longest road of its own was no
  /// record, with that road's layers -- each in the last iteration that saw
  /// one, and every reason any iteration gave, a bit per RoadLoss and FitLoss;
  /// for a chain raised from below, the fit of its own part.
  /// None where nothing was lost; none at all for a particle not findable.
  struct ParticleVerdict {
    FunnelStage lostAt = FunnelStage::kSelection;
    std::optional<PairLoss> pairLoss;
    std::optional<RoadLoss> roadLoss;
    std::optional<FitLoss> fitLoss;
    std::uint64_t fitRoadLayers = 0U; ///< the refused road's layers, a bit each; 0 without a fit loss
    std::uint8_t roadLosses = 0U;     ///< a bit per RoadLoss of any iteration, kNoChain's never
    std::uint8_t fitLosses = 0U;      ///< a bit per FitLoss of any iteration, kNoRoad's never
    std::optional<FitLoss> raisedFit; ///< the fit of a raised chain's own part, kOther if it takes it (`raisedFitLosses`)
  };
  [[nodiscard]] std::optional<ParticleVerdict> particleVerdict(std::uint64_t id) const
  {
    const auto found = mIndex.find(id);
    if (found == mIndex.end()) {
      return std::nullopt;
    }
    const auto& particle = mParticles[found->second];
    std::array<std::uint32_t, kPairLosses> lost{};
    for (const auto& pair : particle.lostPairs) {
      if (((particle.pairsKept >> pair.from) & 1U) != 0U || (pair.iteration < mMinPtGeV.size() && static_cast<double>(particle.pt) < mMinPtGeV[pair.iteration])) {
        continue;
      }
      ++lost[static_cast<std::size_t>(pair.loss)];
    }
    ParticleVerdict result{
      .lostAt = firstLoss(particle.longest, static_cast<std::uint32_t>(mSettings.minLayers), false),
      .pairLoss = std::nullopt,
      .roadLoss = std::nullopt,
      .fitLoss = std::nullopt,
      .fitRoadLayers = 0U,
      .roadLosses = 0U,
      .fitLosses = 0U,
      .raisedFit = std::nullopt,
    };
    if (const auto* const most = std::ranges::max_element(lost); *most > 0U) {
      result.pairLoss = static_cast<PairLoss>(most - lost.begin());
    }
    for (const auto loss : particle.roadLosses) {
      if (loss != RoadLoss::kNoChain) {
        result.roadLoss = loss;
        result.roadLosses |= static_cast<std::uint8_t>(1U << static_cast<unsigned>(loss));
      }
    }
    for (std::size_t iteration = 0; iteration < particle.fitLosses.size(); ++iteration) {
      if (const auto loss = particle.fitLosses[iteration]; loss != FitLoss::kNoRoad) {
        result.fitLoss = loss;
        result.fitRoadLayers = particle.fitRoadLayers[iteration];
        result.fitLosses |= static_cast<std::uint8_t>(1U << static_cast<unsigned>(loss));
      }
    }
    for (const auto loss : particle.raisedFits) {
      if (loss != FitLoss::kNoRoad) {
        result.raisedFit = loss;
      }
    }
    return result;
  }

  /// The chi2 of \a iteration's record on exactly \a rows, in any order
  /// (the best of its duplicates); none if no slice gave one, or without
  /// `keepReportInputs(true)` (for the road fakes).
  [[nodiscard]] std::optional<float> recordChi2(std::size_t iteration, std::span<const std::uint32_t> rows) const
  {
    const auto found = mRecordChi2.find(recordKey(iteration, rows));
    if (found == mRecordChi2.end()) {
      return std::nullopt;
    }
    return found->second;
  }

  /// Each findable particle's iterations with a road of its own, a bit each,
  /// by particle id.
  [[nodiscard]] std::unordered_map<std::uint64_t, std::uint64_t> seededBy() const
  {
    std::unordered_map<std::uint64_t, std::uint64_t> result;
    for (const auto& particle : mParticles) {
      result[particle.id] = particle.seededBy;
    }
    return result;
  }

  /// The report, against the run's \a tracks.
  [[nodiscard]] std::vector<std::string> report(const TrackRows& tracks) const
  {
    std::vector<bool> found(mParticles.size(), false);
    std::vector<FakeKind> fakes(mParticles.size(), FakeKind::kNone);
    for (const auto& track : tracks) {
      const auto particle = goodParticle(mRows, track);
      const auto at = mIndex.find(particle);
      if (at != mIndex.end()) {
        found[at->second] = true;
        continue;
      }
      // A fake: its label is the findable particle most of its measurements carry.
      std::uint32_t label = 0U;
      std::size_t most = 0U;
      for (const auto row : track) {
        for (auto candidate = mRowOffsets[row]; candidate < mRowOffsets[row + 1U]; ++candidate) {
          const auto id = mParticles[mRowParticles[candidate]].id;
          const auto count = static_cast<std::size_t>(std::ranges::count_if(track, [&](std::uint32_t other) { return mRows[other].carries(id); }));
          if (count > most) {
            most = count;
            label = mRowParticles[candidate];
          }
        }
      }
      if (most == 0U || most == track.size()) {
        continue;
      }
      const auto id = mParticles[label].id;
      const auto innermost = *std::ranges::min_element(track, {}, [&](std::uint32_t row) { return mRows[row].layer; });
      auto kind = FakeKind::kMore;
      if (most + 1U == track.size()) {
        kind = mRows[innermost].carries(id) ? FakeKind::kOtherOne : FakeKind::kInnermost;
      }
      auto& fake = fakes[label];
      fake = fake == FakeKind::kNone ? kind : std::min(fake, kind);
    }
    const auto bins = mSettings.productionEdges.size() + 1U;
    struct Row {
      std::uint64_t findable = 0U;
      std::array<std::uint64_t, kFunnelStages> byStage{};
      std::vector<std::array<std::uint64_t, kPairLosses>> reasons;
      std::vector<std::array<std::uint64_t, kRoadLosses>> roadLosses;
      std::vector<std::array<std::uint64_t, kFitLosses>> fitLosses;
      std::array<std::uint64_t, kFakeKinds> fakes{}; ///< of the missed
      /// Lost at the road fit, by species: all of them, and by the chi2 in the
      /// last iteration that gave them a road.
      std::map<std::int32_t, std::array<std::uint64_t, 2>> fitBySpecies;
    };
    // Physical primaries, material, other decays: one row each; V0s and
    // hyperons, decays in flight: one per production bin, then their sum.
    std::vector<std::pair<std::string, Row>> rows;
    const auto rowOf = [&](const Particle& particle) -> std::size_t {
      switch (particle.group) {
        case kPrimary:
          return 0U;
        case kV0:
          return 1U + particle.bin;
        case kInFlight:
          return 2U + bins + particle.bin;
        case kOtherDecay:
          return 3U + (2U * bins);
        default:
          return 4U + (2U * bins);
      }
    };
    rows.resize(5U + (2U * bins));
    std::string edges = "0";
    for (const auto edge : mSettings.productionEdges) {
      edges += std::format(", {:.1f}", edge);
    }
    const auto binName = [&](std::size_t bin) {
      const auto low = bin == 0U ? 0. : mSettings.productionEdges[bin - 1U];
      return bin < mSettings.productionEdges.size() ? std::format("[{:.1f},{:.1f})", low, mSettings.productionEdges[bin]) : std::format("[{:.1f},inf)", low);
    };
    rows[0].first = "physical primaries";
    for (std::size_t bin = 0; bin < bins; ++bin) {
      rows[1U + bin].first = "V0s, hyperons " + binName(bin);
      rows[2U + bins + bin].first = "in flight " + binName(bin);
    }
    rows[1U + bins].first = "V0s, hyperons, all";
    rows[2U + (2U * bins)].first = "in flight, all";
    rows[3U + (2U * bins)].first = "other decays";
    rows[4U + (2U * bins)].first = "material";
    const auto add = [&](Row& row, const Particle& particle, FunnelStage stage, FakeKind fake) {
      ++row.findable;
      ++row.byStage[static_cast<std::size_t>(stage)];
      if (stage != FunnelStage::kFound) {
        ++row.fakes[static_cast<std::size_t>(fake)];
      }
      if (stage == FunnelStage::kFit) {
        auto last = FitLoss::kNoRoad;
        for (const auto loss : particle.fitLosses) {
          last = loss == FitLoss::kNoRoad ? last : loss;
        }
        auto& species = row.fitBySpecies[particle.species];
        ++species[0];
        species[1] += last == FitLoss::kChi2 ? 1U : 0U;
        row.fitLosses.resize(std::max(row.fitLosses.size(), mIterations));
        for (std::size_t iteration = 0; iteration < mIterations; ++iteration) {
          const auto loss = iteration < particle.fitLosses.size() ? particle.fitLosses[iteration] : FitLoss::kNoRoad;
          ++row.fitLosses[iteration][static_cast<std::size_t>(loss)];
        }
      }
      if (stage == FunnelStage::kRoads) {
        row.roadLosses.resize(std::max(row.roadLosses.size(), mIterations));
        for (std::size_t iteration = 0; iteration < mIterations; ++iteration) {
          const auto loss = iteration < particle.roadLosses.size() ? particle.roadLosses[iteration] : RoadLoss::kNoChain;
          ++row.roadLosses[iteration][static_cast<std::size_t>(loss)];
        }
      }
      if (stage == FunnelStage::kTracklets) {
        row.reasons.resize(std::max(row.reasons.size(), particle.reasons.size()));
        for (std::size_t iteration = 0; iteration < particle.reasons.size(); ++iteration) {
          for (std::size_t reason = 0; reason < kPairLosses; ++reason) {
            row.reasons[iteration][reason] += particle.reasons[iteration][reason];
          }
        }
      }
    };
    for (std::size_t index = 0; index < mParticles.size(); ++index) {
      const auto& particle = mParticles[index];
      const auto stage = firstLoss(particle.longest, static_cast<std::uint32_t>(mSettings.minLayers), found[index]);
      add(rows[rowOf(particle)].second, particle, stage, fakes[index]);
      if (particle.group == kV0) {
        add(rows[1U + bins].second, particle, stage, fakes[index]);
      } else if (particle.group == kInFlight) {
        add(rows[2U + (2U * bins)].second, particle, stage, fakes[index]);
      }
    }
    std::vector<std::string> lines{
      std::format("stage funnel: the findable above {} GeV on a chain of {} layers, production {} [cm] in bins from {}; "
                  "found %, then the missed by the first stage without a chain of {} of their own (% of the findable)",
                  mSettings.minPtGeV, mSettings.minLayers, mSettings.productionAlongZ ? "|z|" : "radius", edges, mSettings.minLayers)};
    std::string seeded = std::format("  seeded, a road of its own, of {} findable, by iteration:", mParticles.size());
    for (std::size_t iteration = 0; iteration < mIterations; ++iteration) {
      seeded += std::format(" {} {}", iteration, std::ranges::count_if(mParticles, [&](const Particle& p) { return ((p.seededBy >> iteration) & 1U) != 0U; }));
    }
    lines.push_back(seeded + std::format(", any {}", std::ranges::count_if(mParticles, [](const Particle& p) { return p.seededBy != 0U; })));
    for (const auto& [name, row] : rows) {
      if (row.findable == 0U) {
        continue;
      }
      std::string line = std::format("  {:<28} {:>6} findable, found {:5.1f} %; lost at", name, row.findable, percent(row.byStage[static_cast<std::size_t>(FunnelStage::kFound)], row.findable));
      for (std::size_t stage = 0; stage + 1U < kFunnelStages; ++stage) {
        line += std::format(" {} {:.1f}", kFunnelStageNames[stage], percent(row.byStage[stage], row.findable));
      }
      lines.push_back(line);
      for (std::size_t iteration = 0; iteration < row.reasons.size(); ++iteration) {
        std::uint64_t total = 0U;
        for (const auto count : row.reasons[iteration]) {
          total += count;
        }
        std::string reasons = std::format("    lost at tracklets, iteration {}: {} adjacent pairs not kept:", iteration, total);
        for (std::size_t reason = 0; reason < kPairLosses; ++reason) {
          reasons += std::format(" {} {:.1f}", kPairLossNames[reason], percent(row.reasons[iteration][reason], total));
        }
        lines.push_back(reasons);
      }
      for (std::size_t iteration = 0; iteration < row.roadLosses.size(); ++iteration) {
        const auto lost = row.byStage[static_cast<std::size_t>(FunnelStage::kRoads)];
        std::string reasons = std::format("    lost at roads, iteration {}, of {}:", iteration, lost);
        for (std::size_t reason = 0; reason < kRoadLosses; ++reason) {
          reasons += std::format(" {} {:.1f}", kRoadLossNames[reason], percent(row.roadLosses[iteration][reason], lost));
        }
        lines.push_back(reasons);
      }
      for (std::size_t iteration = 0; iteration < row.fitLosses.size(); ++iteration) {
        const auto lost = row.byStage[static_cast<std::size_t>(FunnelStage::kFit)];
        std::string reasons = std::format("    lost at the road fit, iteration {}, of {}:", iteration, lost);
        for (std::size_t reason = 0; reason < kFitLosses; ++reason) {
          reasons += std::format(" {} {:.1f}", kFitLossNames[reason], percent(row.fitLosses[iteration][reason], lost));
        }
        lines.push_back(reasons);
      }
      const auto missed = row.findable - row.byStage[static_cast<std::size_t>(FunnelStage::kFound)];
      std::string fakeLine = std::format("    the {} missed, by their fake:", missed);
      for (std::size_t kind = 0; kind < kFakeKinds; ++kind) {
        fakeLine += std::format(" {} {:.1f}", kFakeKindNames[kind], percent(row.fakes[kind], missed));
      }
      lines.push_back(fakeLine);
      if (!row.fitBySpecies.empty()) {
        std::string species = "    lost at the road fit by species, all (by chi2):";
        for (const auto& [pdg, counts] : row.fitBySpecies) {
          species += std::format(" {} {} ({})", speciesName(pdg), counts[0], counts[1]);
        }
        lines.push_back(species);
      }
    }
    return lines;
  }

  /// Why the particles found are found shorter than they could be:
  /// for each findable particle with a
  /// good track among \a tracks, its ceiling -- the longest chain its own
  /// measurements make over the layout's links, any link, whatever an
  /// iteration takes -- against its longest good track, and, if shorter, the
  /// first stage whose longest chain of its own, over the iterations, falls
  /// below the ceiling; the selection if its records reached it.
  ///
  /// A particle short at the tracklets with a direct pair not kept: each
  /// pair of its longest chain of direct pairs that no iteration kept, by the
  /// first test it fails in the iteration of its track, \a iterations[i]
  /// being track i's.
  [[nodiscard]] std::vector<std::string> lengthReport(const TrackRows& tracks, std::span<const std::size_t> iterations, const Controller& tracker) const
  {
    std::vector<std::uint32_t> found(mParticles.size(), 0U);
    std::vector<std::size_t> foundIn(mParticles.size(), 0U);
    for (std::size_t track = 0; track < tracks.size(); ++track) {
      const auto at = mIndex.find(goodParticle(mRows, tracks[track]));
      if (at != mIndex.end() && tracks[track].size() > found[at->second]) {
        found[at->second] = static_cast<std::uint32_t>(tracks[track].size());
        foundIn[at->second] = iterations[track];
      }
    }
    constexpr std::size_t kLayerBits = std::numeric_limits<std::uint64_t>::digits;
    const auto links = [] {
      std::array<std::uint64_t, kLayerBits> result{};
      for (const auto& link : Det::layerLinks) {
        result[link.from] |= std::uint64_t{1} << link.to;
      }
      return result;
    }();
    const std::uint64_t barrels = barrelLayerBits<Det>();
    std::array<LayerSpan, kLayerBits> spans{};
    for (std::size_t layer = 0; layer < Detector<Det>::nLayers; ++layer) {
      spans[layer] = layerSpan<Det>(static_cast<LayerID>(layer));
    }
    std::array<std::uint32_t, kLayerBits> order{};
    for (std::size_t position = 0; position < layerOrder<Det>.size(); ++position) {
      order[layerOrder<Det>[position]] = static_cast<std::uint32_t>(position);
    }
    struct Counts {
      std::uint64_t found = 0U;
      std::uint64_t complete = 0U;
      std::uint64_t belowLayers = 0U; ///< the ceiling below its distinct layers
      std::array<std::uint64_t, kFunnelStages - 1U> shortAt{};
      std::array<std::uint64_t, kFunnelStages - 1U> missingAt{};      ///< ceiling - found, summed
      std::array<std::uint64_t, kFunnelStages - 1U> stageShortfall{}; ///< ceiling - that stage's chain, summed
      /// Short at the tracklets: their chain as long as the direct pairs'
      /// longest, so only a pair over a layer crossed without a row joins the
      /// rest; or shorter, a direct pair not kept.
      std::array<std::uint64_t, 2> trackletsBy{};
      std::uint64_t overMiss = 0U;                           ///< of them, the direct pairs' longest below the ceiling
      std::array<std::uint64_t, kPairLosses> missingPairs{}; ///< the direct pairs not kept, by reason
      std::array<std::uint64_t, kPairLosses> firstMissing{}; ///< the innermost of each particle's
      std::map<std::pair<LayerID, LayerID>, std::uint64_t> missingLinks;
    };
    std::array<Counts, 2> counts; // physical primaries, all
    for (std::size_t index = 0; index < mParticles.size(); ++index) {
      if (found[index] == 0U) {
        continue;
      }
      const auto& particle = mParticles[index];
      std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
      std::uint64_t layers = 0U;
      for (const auto a : particle.rows) {
        layers |= std::uint64_t{1} << mRows[a].layer;
        for (const auto b : particle.rows) {
          if (((links[mRows[a].layer] >> mRows[b].layer) & 1U) != 0U) {
            edges.emplace_back(a, b);
          }
        }
      }
      // The same over the direct pairs alone, which cross no other layer
      // inside its span: what an iteration without skips can join.
      std::vector<std::pair<std::uint32_t, std::uint32_t>> direct;
      for (const auto& [a, b] : edges) {
        const auto holes = linkHoleLayers<Det>(mRows[a].layer, mRows[b].layer);
        if (!segmentLeavesHole(holes, barrels, spans.data(), mRows[a].r, mRows[a].z, mRows[b].r, mRows[b].z)) {
          direct.emplace_back(a, b);
        }
      }
      const auto layerOrderOf = [&](std::uint32_t row) { return order[mRows[row].layer]; };
      const auto directPath = longestPath(direct, layerOrderOf);
      const auto directCeiling = static_cast<std::uint32_t>(directPath.size());
      const auto ceiling = std::max(longestChain(std::move(edges), layerOrderOf), found[index]);
      std::size_t stage = kFunnelStages - 2U; // the selection
      for (std::size_t s = 0; s < particle.longest.size(); ++s) {
        if (particle.longest[s] < ceiling) {
          stage = s;
          break;
        }
      }
      for (std::size_t group = 0; group < counts.size(); ++group) {
        if (group == 0U && particle.group != kPrimary) {
          continue;
        }
        auto& count = counts[group];
        ++count.found;
        count.belowLayers += ceiling < static_cast<std::uint32_t>(std::popcount(layers)) ? 1U : 0U;
        if (found[index] >= ceiling) {
          ++count.complete;
          continue;
        }
        ++count.shortAt[stage];
        if (stage == 0U) {
          ++count.trackletsBy[particle.longest[0] >= directCeiling ? 0U : 1U];
          count.overMiss += directCeiling < ceiling ? 1U : 0U;
          if (particle.longest[0] < directCeiling) {
            const auto& parameters = tracker.iteration(foundIn[index]);
            const auto kept = mOwnKept.find(static_cast<std::uint32_t>(index));
            bool first = true;
            for (std::size_t i = 1; i < directPath.size(); ++i) {
              const auto a = directPath[i - 1U];
              const auto b = directPath[i];
              if (kept != mOwnKept.end() && kept->second.contains(pairKey(a, b))) {
                continue;
              }
              const auto link = linkOf<Det>(mRows[a].layer, mRows[b].layer);
              auto loss = link == kNoLink || !linkTaken(parameters, link) ? PairLoss::kLinkNotTaken : pairCutLoss(*mFrame, parameters, mBeamLine, link, a, b).loss;
              if (loss == PairLoss::kCurvature && foundIn[index] < mMinPtGeV.size() && static_cast<double>(particle.pt) < mMinPtGeV[foundIn[index]]) {
                loss = PairLoss::kCurvatureBelowPtMin;
              }
              ++count.missingPairs[static_cast<std::size_t>(loss)];
              ++count.missingLinks[{mRows[a].layer, mRows[b].layer}];
              if (first) {
                ++count.firstMissing[static_cast<std::size_t>(loss)];
                first = false;
              }
            }
          }
        }
        count.missingAt[stage] += ceiling - found[index];
        count.stageShortfall[stage] += stage < particle.longest.size() ? ceiling - particle.longest[stage] : 0U;
      }
    }
    std::vector<std::string> lines{
      "length funnel: the found against their ceiling, the longest chain of their own over the layout's links; the shorter by the first stage "
      "whose longest chain of their own falls below it: particles (layers the found track misses, layers that stage's chain misses)"};
    constexpr std::array<const char*, 2> kNames{"physical primaries", "all"};
    for (std::size_t group = 0; group < counts.size(); ++group) {
      const auto& count = counts[group];
      std::string line = std::format("  {:<18} {} found, complete {:.1f} %, ceiling below their layers {:.1f} %; shorter at", kNames[group], count.found, percent(count.complete, count.found), percent(count.belowLayers, count.found));
      for (std::size_t stage = 0; stage < count.shortAt.size(); ++stage) {
        line += std::format(" {} {} ({}, {})", kFunnelStageNames[stage], count.shortAt[stage], count.missingAt[stage], count.stageShortfall[stage]);
      }
      lines.push_back(line);
      lines.push_back(
        std::format("    short at the tracklets: {} with a chain as long as their direct pairs' longest (only a pair over a layer crossed "
                    "without a row joins the rest), {} with a direct pair not kept; {} whose direct pairs' longest is below the ceiling",
                    count.trackletsBy[0], count.trackletsBy[1], count.overMiss));
      std::string reasons = "    their direct pairs not kept, by the first test failing in the iteration of their track, all (the innermost):";
      for (std::size_t reason = 0; reason < kPairLosses; ++reason) {
        reasons += std::format(" {} {} ({})", kPairLossNames[reason], count.missingPairs[reason], count.firstMissing[reason]);
      }
      lines.push_back(reasons);
      std::vector<std::pair<std::uint64_t, std::pair<LayerID, LayerID>>> byLink;
      byLink.reserve(count.missingLinks.size());
      for (const auto& [link, pairs] : count.missingLinks) {
        byLink.emplace_back(pairs, link);
      }
      std::ranges::sort(byLink, std::greater{});
      std::string top = "    by link (layers), the most:";
      for (std::size_t i = 0; i < std::min<std::size_t>(byLink.size(), 12U); ++i) {
        top += std::format(" {}-{}: {}", byLink[i].second.first, byLink[i].second.second, byLink[i].first);
      }
      lines.push_back(top);
    }
    return lines;
  }

  /// The longest chain of \a edges, as its nodes in order; \a orderOf orders
  /// the nodes' layers.
  template <typename OrderOf>
  [[nodiscard]] static std::vector<std::uint32_t> longestPath(std::vector<std::pair<std::uint32_t, std::uint32_t>> edges, OrderOf orderOf)
  {
    std::ranges::sort(edges, {}, [&](const auto& edge) { return orderOf(edge.first); });
    std::unordered_map<std::uint32_t, std::pair<std::uint32_t, std::uint32_t>> best; ///< node: length, predecessor
    std::uint32_t end = 0U;
    std::uint32_t longest = 0U;
    for (const auto& [from, to] : edges) {
      const auto before = best.find(from);
      const std::uint32_t reached = (before == best.end() ? 1U : before->second.first) + 1U;
      auto& at = best[to];
      if (reached > at.first) {
        at = {reached, from};
        if (reached > longest) {
          longest = reached;
          end = to;
        }
      }
    }
    std::vector<std::uint32_t> path;
    if (longest == 0U) {
      return path;
    }
    for (std::uint32_t node = end;;) {
      path.push_back(node);
      const auto at = best.find(node);
      if (at == best.end()) {
        break;
      }
      node = at->second.second;
    }
    std::ranges::reverse(path);
    return path;
  }

 private:
  enum Group : std::uint8_t { kPrimary,
                              kV0,
                              kInFlight,
                              kOtherDecay,
                              kMaterial };

  /// An adjacent-layer pair of a particle's own measurements an iteration did
  /// not keep: the iteration, the pair's inner layer, the first test failing.
  struct LostPair {
    std::size_t iteration = 0U;
    LayerID from = 0;
    PairLoss loss = PairLoss::kElsewhere;
  };

  struct Particle {
    std::uint64_t id = 0U;
    float pt = 0.F;           ///< [GeV]
    std::int32_t species = 0; ///< |PDG code|
    Group group = kMaterial;
    std::size_t bin = 0U;
    std::vector<std::uint32_t> rows;
    std::array<std::uint32_t, 5> longest{};                      ///< per stage, tracklets to the fit
    std::vector<std::array<std::uint32_t, 5>> longestBy;         ///< the same per iteration
    std::uint64_t seededBy = 0U;                                 ///< a bit per iteration with a road of its own
    std::vector<std::array<std::uint32_t, kPairLosses>> reasons; ///< per iteration
    std::vector<RoadLoss> roadLosses;                            ///< per iteration
    std::vector<FitLoss> fitLosses;                              ///< per iteration
    std::vector<float> recordChi2;                               ///< per iteration, the best chi2 of its longest records
    /// The inner layers of its adjacent-layer pairs of its own some iteration
    /// kept, a bit each.
    std::uint64_t pairsKept = 0U;
    std::vector<LostPair> lostPairs;          ///< every iteration's pairs of its own not kept there
    std::vector<std::uint64_t> fitRoadLayers; ///< per iteration, the layers of the road `fitLosses` fitted, a bit each
    std::vector<FitLoss> raisedFits;          ///< per iteration, the fit of a raised chain's own part (`raisedFitLosses`)
  };

  /// Why each particle's longest road of its own in \a iteration, if it
  /// reaches the minimum and gave no record as long, was no record.
  void fitLosses(const IterationParameters<Det>& parameters, std::size_t iteration)
  {
    if (mTracking == nullptr || !mPropagator) {
      return;
    }
    const auto minimum = static_cast<std::size_t>(mSettings.minLayers);
    for (const auto& [particle, rows] : mOwnRoad) {
      if (rows.size() < minimum) {
        continue;
      }
      const auto record = mOwnRecord.find(particle);
      if (record != mOwnRecord.end() && record->second >= minimum) {
        continue;
      }
      auto& losses = mParticles[particle].fitLosses;
      losses.resize(std::max(losses.size(), iteration + 1U), FitLoss::kNoRoad);
      auto& layers = mParticles[particle].fitRoadLayers;
      layers.resize(std::max(layers.size(), iteration + 1U), 0U);
      for (const auto row : rows) {
        layers[iteration] |= std::uint64_t{1} << mRows[row].layer;
      }
      losses[iteration] = fitLoss(parameters, *mPropagator, rows);
    }
  }

  /// The fit of the own part of each chain the walk raised from below it
  /// (`raisedHow`) with no record of its own: kOther if the fit, the floor
  /// and the beam line all take it -- a road without the other's triplets
  /// below would have been a record.
  void raisedFitLosses(const IterationParameters<Det>& parameters, std::size_t iteration)
  {
    if (mTracking == nullptr || !mPropagator) {
      return;
    }
    const auto minimum = static_cast<std::size_t>(mSettings.minLayers);
    for (const auto& [particle, rows] : mRaisedChain) {
      const auto record = mOwnRecord.find(particle);
      if (rows.size() < minimum || (record != mOwnRecord.end() && record->second >= minimum)) {
        continue;
      }
      auto& fits = mParticles[particle].raisedFits;
      fits.resize(std::max(fits.size(), iteration + 1U), FitLoss::kNoRoad);
      fits[iteration] = fitLoss(parameters, *mPropagator, rows);
    }
  }

  /// Why the road of \a rows is no record (`roadFitRefusal`): its boxes
  /// share no time, the fit's chi2 or propagation, the pT floor for its
  /// length, the beam line; kOther if none.
  [[nodiscard]] FitLoss fitLoss(const IterationParameters<Det>& parameters, const Propagator& propagator, const std::vector<std::uint32_t>& rows) const
  {
    switch (roadFitRefusal(parameters, propagator, *mTracking, *mMaterial, *mFrame, std::span<const std::uint32_t>{rows})) {
      case FitRefusal::kTime:
        return FitLoss::kTime;
      case FitRefusal::kChi2:
        return FitLoss::kChi2;
      case FitRefusal::kPropagation:
        return FitLoss::kPropagation;
      case FitRefusal::kPtFloor:
        return FitLoss::kPtFloor;
      case FitRefusal::kBeamLine:
        return FitLoss::kBeamLine;
      case FitRefusal::kTaken:
        break;
    }
    return FitLoss::kOther;
  }

  /// What the road counting saw of a triplet.
  struct TripletInfo {
    std::uint32_t level = 0U;
    bool successorAbove = false;
    LayerID endLayer = 0;
  };

  /// Why \a particle's chain of \a length triplets ending at \a end is no road
  /// of its own in \a iteration; \a chains, its neighbours.
  void roadLoss(std::uint32_t particle, std::size_t iteration, const IterationParameters<Det>& parameters, std::uint32_t length, const std::pair<std::uint64_t, std::uint32_t>& end, const std::vector<std::array<std::uint32_t, 4>>& chains)
  {
    auto& losses = mParticles[particle].roadLosses;
    losses.resize(std::max(losses.size(), iteration + 1U), RoadLoss::kNoChain);
    const auto& infos = mTripletInfo[particle];
    const auto found = infos.find(end);
    auto loss = RoadLoss::kOther;
    if (found != infos.end()) {
      const auto& info = found->second;
      if (info.level > length) {
        std::vector<std::uint32_t> own;
        loss = raisedHow(infos, end, chains, mTriplets[particle], own);
        if (loss != RoadLoss::kRaised) {
          mRaisedChain[particle] = std::move(own);
        }
      } else if (info.level < parameters.minRoadStartLevel) {
        loss = RoadLoss::kStartLevel;
      } else if (parameters.roadEndLayers != 0U && ((parameters.roadEndLayers >> info.endLayer) & 1U) == 0U) {
        loss = RoadLoss::kEndLayer;
      } else if (parameters.roadEndLayers == 0U && info.successorAbove) {
        loss = RoadLoss::kContinued;
      }
    }
    losses[iteration] = loss;
  }

  /// How a chain whose end \a end is raised loses its road: the walk from
  /// it, descending through triplets of its own exactly one level below
  /// (\a chains, its neighbours), reaches one without a neighbour of its own
  /// below -- with no triplet of its own before it (\a triplets, as pairs of
  /// tracklets) or with one -- or never does: it leaves the chain. In the
  /// first two, \a own is the walked part of the chain, its rows innermost
  /// first: the road the walk would give without the other's triplets below.
  [[nodiscard]] static RoadLoss raisedHow(const std::map<std::pair<std::uint64_t, std::uint32_t>, TripletInfo>& infos, const std::pair<std::uint64_t, std::uint32_t>& end, const std::vector<std::array<std::uint32_t, 4>>& chains, const std::vector<std::pair<std::uint64_t, std::uint64_t>>& triplets, std::vector<std::uint32_t>& own)
  {
    using Key = std::pair<std::uint64_t, std::uint32_t>;
    std::vector<Key> stack{end};
    std::map<Key, Key> above{{end, end}}; // each triplet walked, and the one it was walked from
    while (!stack.empty()) {
      const auto current = stack.back();
      stack.pop_back();
      const auto level = infos.at(current).level;
      bool predecessor = false;
      for (const auto& rows : chains) {
        if (pairKey(rows[1], rows[2]) != current.first || rows[3] != current.second) {
          continue;
        }
        const std::pair<std::uint64_t, std::uint32_t> below{pairKey(rows[0], rows[1]), rows[2]};
        const auto found = infos.find(below);
        if (found == infos.end()) {
          continue;
        }
        predecessor = true;
        if (found->second.level + 1U == level && above.try_emplace(below, current).second) {
          stack.push_back(below);
        }
      }
      if (!predecessor) {
        own = {static_cast<std::uint32_t>(current.first >> 32U), static_cast<std::uint32_t>(current.first), current.second};
        for (auto at = current; at != end;) {
          at = above.at(at);
          own.push_back(at.second);
        }
        const bool before = std::ranges::any_of(triplets, [&](const auto& triplet) { return triplet.second == current.first; });
        return before ? RoadLoss::kRaisedUnlinked : RoadLoss::kRaisedBelow;
      }
    }
    return RoadLoss::kRaised;
  }

  [[nodiscard]] static Group groupOf(const ParticleTruth& truth)
  {
    switch (truth.origin) {
      case ParticleOrigin::kPhysicalPrimary:
        return kPrimary;
      case ParticleOrigin::kDecay: {
        const auto family = decayFamily(truth.motherPdg);
        if (family == DecayFamily::kOther) {
          return kOtherDecay;
        }
        return inFlight(family) ? kInFlight : kV0;
      }
      default:
        return kMaterial;
    }
  }

  /// The first iteration reads an empty mask.
  void beginIteration()
  {
    if (mUsed.empty()) {
      mUsed.assign(mRows.size(), std::uint8_t{0});
    }
  }

  /// Calls \a body with every findable particle all of \a rows carry; whether
  /// there was one.
  template <typename Body>
  bool forEachOwner(std::span<const std::uint32_t> rows, Body body) const
  {
    bool any = false;
    forEachCommonParticle(mRows, rows, [&](std::uint64_t id) {
      if (const auto at = mIndex.find(id); at != mIndex.end()) {
        body(at->second);
        any = true;
      }
    });
    return any;
  }

  /// A record's key: \a iteration and its rows, sorted, hashed (FNV-1a).
  [[nodiscard]] static std::uint64_t recordKey(std::size_t iteration, std::span<const std::uint32_t> rows)
  {
    std::vector<std::uint32_t> sorted(rows.begin(), rows.end());
    std::ranges::sort(sorted);
    std::uint64_t hash = 14695981039346656037ULL;
    const auto mix = [&](std::uint64_t value) {
      hash ^= value;
      hash *= 1099511628211ULL;
    };
    mix(iteration);
    for (const auto row : sorted) {
      mix(row);
    }
    return hash;
  }

  void raise(std::uint32_t particle, std::size_t iteration, FunnelStage stage, std::uint32_t length)
  {
    auto& p = mParticles[particle];
    auto& longest = p.longest[static_cast<std::size_t>(stage)];
    longest = std::max(longest, length);
    p.longestBy.resize(std::max(p.longestBy.size(), iteration + 1U));
    auto& here = p.longestBy[iteration][static_cast<std::size_t>(stage)];
    here = std::max(here, length);
  }

  /// Every findable particle's adjacent-layer pairs of its own measurements
  /// that this iteration did not keep, by the first test each fails.
  void lostPairReasons(const IterationParameters<Det>& parameters, std::size_t iteration)
  {
    // The cut's vertex is the beam line's here; an iteration from the vertex
    // list keeps none.
    if (parameters.vertexSource != VertexSource::kBeamLine) {
      return;
    }
    // Each particle its own: in parallel, the same as one after the other.
    detail::parallelFor(mArena, 0U, mParticles.size(), kParticlesPerTask, [&](std::size_t begin, std::size_t end) {
      for (std::size_t index = begin; index < end; ++index) {
        lostPairReasons(parameters, iteration, mParticles[index]);
      }
    });
  }

  /// lostPairReasons for one particle.
  void lostPairReasons(const IterationParameters<Det>& parameters, std::size_t iteration, Particle& particle) const
  {
    auto& reasons = particle.reasons;
    reasons.resize(std::max(reasons.size(), iteration + 1U));
    for (std::size_t first = 0; first < particle.rows.size(); ++first) {
      const auto fromLayer = mRows[particle.rows[first]].layer;
      const auto toLayer = static_cast<LayerID>(fromLayer + 1U);
      const auto link = linkOf<Det>(fromLayer, toLayer);
      // The layer pair once, from its first measurement on the inner layer.
      if (std::any_of(particle.rows.begin(), particle.rows.begin() + static_cast<std::ptrdiff_t>(first), [&](std::uint32_t row) { return mRows[row].layer == fromLayer; })) {
        continue;
      }
      bool any = false;
      bool kept = false;
      auto furthest = PairLoss::kLinkNotTaken;
      for (const auto source : particle.rows) {
        if (mRows[source].layer != fromLayer) {
          continue;
        }
        for (const auto target : particle.rows) {
          if (mRows[target].layer != toLayer) {
            continue;
          }
          any = true;
          if (mKept.contains(pairKey(source, target))) {
            kept = true;
            particle.pairsKept |= std::uint64_t{1} << fromLayer;
            break;
          }
          auto loss = link == kNoLink ? PairLoss::kLinkNotTaken : pairLoss(parameters, link, source, target);
          if (loss == PairLoss::kCurvature && iteration < mMinPtGeV.size() && static_cast<double>(particle.pt) < mMinPtGeV[iteration]) {
            loss = PairLoss::kCurvatureBelowPtMin;
          }
          furthest = std::max(furthest, loss);
        }
      }
      if (any && !kept) {
        ++reasons[iteration][static_cast<std::size_t>(furthest)];
        particle.lostPairs.push_back({.iteration = iteration, .from = fromLayer, .loss = furthest});
      }
    }
  }

  [[nodiscard]] PairLoss pairLoss(const IterationParameters<Det>& parameters, std::size_t link, std::uint32_t sourceRow, std::uint32_t targetRow) const
  {
    if (!linkTaken(parameters, link)) {
      return PairLoss::kLinkNotTaken;
    }
    if (mUsed[sourceRow] != 0U || mUsed[targetRow] != 0U) {
      return PairLoss::kUsed;
    }
    return pairCutLoss(*mFrame, parameters, mBeamLine, link, sourceRow, targetRow).loss;
  }

  static constexpr std::size_t kParticlesPerTask = 256U;

  std::span<const TruthRow> mRows;
  ComparisonSettings mSettings;
  VertexSeed mBeamLine;
  const Frame* mFrame = nullptr;
  const TrackingGeometry<Det>* mTracking = nullptr;
  const o2::base::MatLayerCylSet* mMaterial = nullptr;
  std::optional<Propagator> mPropagator;
  std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> mOwnRoad;     ///< each particle's longest road of its own, this iteration
  std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> mRaisedChain; ///< the own part of a chain raised from below it, this iteration
  std::unordered_map<std::uint32_t, std::uint32_t> mOwnRecord;                ///< and its longest record
  std::unordered_map<std::uint64_t, float> mRecordChi2;                       ///< every record's chi2 by recordKey
  std::vector<double> mMinPtGeV;
  std::vector<Particle> mParticles;
  detail::TaskArena mArena; ///< lostPairReasons' threads: every core
  std::unordered_map<std::uint64_t, std::uint32_t> mIndex;
  std::vector<std::uint32_t> mRowOffsets;
  std::vector<std::uint32_t> mRowParticles;
  // The open iteration's, over its slices.
  std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>> mPairs;
  std::unordered_map<std::uint32_t, std::vector<std::pair<std::uint64_t, std::uint64_t>>> mTriplets;
  std::unordered_map<std::uint32_t, std::map<std::pair<std::uint64_t, std::uint32_t>, TripletInfo>> mTripletInfo;
  std::unordered_map<std::uint32_t, std::vector<std::array<std::uint32_t, 4>>> mNeighbours;
  std::unordered_set<std::uint64_t> mKept;
  std::unordered_map<std::uint32_t, std::unordered_set<std::uint64_t>> mOwnKept; ///< each particle's own kept pairs, over the iterations
  std::vector<std::uint8_t> mUsed;                                               ///< the mask the open iteration read
  bool mKeepOwnKept = false;                                                     ///< `keepReportInputs`
  bool mKeepRecordChi2s = false;
  std::size_t mIterations = 0U; ///< observed
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_STAGEFUNNEL_H_
