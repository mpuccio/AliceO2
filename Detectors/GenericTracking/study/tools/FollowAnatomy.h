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

/// \file FollowAnatomy.h
/// \brief The follow's steps on correct roads, against the truth: a bench
///        observer.
///
/// After each iteration with the follow, every one of its tracks whose rows
/// before the follow are all one particle's -- a correct road -- is followed
/// again, outward and inward, from the states the follow started from and
/// with the used mask as the follow read it (the tracker's copy). The replay
/// must propose what the tracking proposed: the rows the track kept are a
/// prefix of the replay's in each direction. Each step is recorded: the
/// targets in their order, the free rows in each window with their chi2, the
/// kink pass and the row taken.
///
/// Per track and direction (0 outward, 1 inward) an `F` line: iteration,
/// track, the particle, its pT, the rows before the follow, the direction,
/// the rows the replay proposed, kept, of them wrong, whether the kept rows
/// are a prefix of the replay's (1), the first wrong step (-1 for none), its
/// class, how the direction ended, and the wrong rows proposed. Per step an
/// `S` line: iteration, track, direction, step, the layer it took a row on
/// (-1 for none), whether it took one, that row's chi2, whether it is the
/// particle's, the rows passing the gate on that layer, the particle's own
/// row there (its chi2, -1 where the step did not see one), the class of the
/// step (`kClassNames`), whether a kink step took it, the rows the lookahead ranked
/// (`lookaheadConfirmation`; 0 where one passed, or in a greedy pass),
/// whether the row taken was confirmed there (1, 0; -1 not ranked) and
/// whether the particle's own row on that layer was (1, 0; -1 not ranked or
/// none passing the gate). After each `S`
/// line a `Q` line per layer the step targeted that holds a row of the
/// particle: iteration, track, direction, step, the layer, whether
/// that row is free, then for the row predicted from the step's own state
/// (the best of the particle's rows there, by chi2) its chi2, its two pulls,
/// the prediction's two widths and the row's two errors [cm] (NaN where the
/// prediction fails).
///
/// `--follow-stops` prints, at the end, why the directions of the correct
/// roads ended (`printStops`).

#ifndef GTRK_STUDY_FOLLOWANATOMY_H_
#define GTRK_STUDY_FOLLOWANATOMY_H_

#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/SlicePlan.h"
#include "GenericTracking/Tracking/TrackFitter.h"
#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Study/TrackingInput.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

/// What a step that took a row did, against its particle's own rows:
/// - own: it took the particle's row;
/// - swap: the particle's row on that layer passed the gate and lost;
/// - refused: the particle's free row on that layer was seen, above the gate;
/// - not in window: the particle has a free row on that layer, not seen;
/// - refused nearer: the particle's free row was on a nearer target the step
///   searched first and did not take;
/// - held: the particle's row on that layer is another track's;
/// - filled hole: no row of the particle on that layer, but free ones further
///   along;
/// - past the end: no FREE row of the particle further along -- none, or only
///   rows other tracks hold.
/// "Further along" is a larger r^2 + z^2 than the direction's last row so far
/// outward, a smaller one inward: a proxy for the rows ahead, wrong for a
/// particle that turns back.
enum class StepClass : std::uint8_t { kOwn,
                                      kSwap,
                                      kRefused,
                                      kNotInWindow,
                                      kRefusedNearer,
                                      kHeld,
                                      kFilledHole,
                                      kPastEnd,
                                      kCount };
inline constexpr const char* kClassNames[] = {"own", "swap", "refused", "not-in-window", "refused-nearer", "held", "filled-hole", "past-end"};

/// `FollowEnd`'s names: no target, the holes, no row (gaps), the capacity;
/// with dashes, a word each, for the `F` lines split at whitespace.
inline constexpr const char* kEndNames[] = {"no-target", "holes", "no-row", "capacity"};

template <DetectorLayout Det, typename Material>
class FollowAnatomy
{
 public:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using Proposal = FollowProposal<kMaxClusters>;

  /// Steps into \a path, none where it is empty (`--follow-stops` alone).
  FollowAnatomy(const std::string& path, const study::TrackingTruth& truth) : mPath(path), mTruth(&truth)
  {
    if (!mPath.empty()) {
      mOut.open(mPath);
      if (!mOut) {
        throw std::runtime_error("cannot write " + path);
      }
    }
    mRowsOf = study::rowsOfParticles(truth.rows);
  }

  /// One iteration of the tracking, seen by its observer.
  template <typename Observation>
  void observe(const Observation& iteration, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    // One run's steps: a repetition writes its own from the start.
    if (iteration.iteration == 0U) {
      if (!mPath.empty()) {
        mOut.close();
        mOut.open(mPath);
      }
      mStops = {};
      mReplayed = 0U;
      mAgree = 0U;
    }
    if (iteration.followStart.empty()) {
      return;
    }
    if (mCells.empty()) {
      mCells = SensorCells<Det>(geometry);
      mWidths = frameBounds(frame).widthPS;
    }
    if (parameters.followEnd) {
      throw std::runtime_error("--follow-steps and --follow-stops cannot replay an iteration with iterationFollowEnd");
    }
    // The follow as the iteration ran it: the beam line, the greedy choice
    // and the end decisions on the density and the beam line.
    withBools(
      [&]<bool BeamLine, bool Greedy, bool EndDensity, bool BeamEnd>(std::bool_constant<BeamLine>, std::bool_constant<Greedy>, std::bool_constant<EndDensity>, std::bool_constant<BeamEnd>) {
        if constexpr (BeamLine || !BeamEnd) {
          replay<BeamLine, Greedy, EndDensity, BeamEnd>(iteration, parameters, frame, geometry, material, propagator);
        }
      },
      parameters.beamConstraint, parameters.followGreedy, Follow<Det>::endDensity(parameters), parameters.followBeamEnd);
  }

  /// Why the directions of the correct roads ended (`--follow-stops`): per
  /// direction and end, the follows whose particle still had a free row
  /// further along, and those without; and how many replays agreed with the
  /// tracking.
  void printStops() const
  {
    std::cout << std::format("follow stops of the correct roads: {} replayed, {} agree with the tracking\n", mReplayed, mAgree);
    for (std::size_t direction = 0; direction < 2U; ++direction) {
      std::string line = std::format("  {}: end (particle's free row further along / none):", direction == 0U ? "outward" : "inward");
      for (std::size_t end = 0; end < mStops[direction].size(); ++end) {
        line += std::format(" {} {}/{}", kEndNames[end], mStops[direction][end][1], mStops[direction][end][0]);
      }
      std::cout << line << '\n';
    }
  }

 private:
  struct Candidate {
    LayerID layer = 0;
    std::uint32_t row = 0;
    Scalar chi2 = 0.F;
  };
  struct Step {
    std::size_t direction = 0U; ///< 0 outward, 1 inward
    BoundTrackState state;
    CovarianceFactors factors;
    std::vector<LayerID> targets;
    std::vector<Candidate> candidates;
    std::uint32_t taken = ClusterList<kMaxClusters>::kNoCluster;
    Scalar chi2 = 0.F;
    bool kinked = false;
    /// The lookahead's confirmations (`lookaheadConfirmation`): each passing row it
    /// ranked, and whether a row of the next step confirmed it.
    std::vector<std::pair<std::uint32_t, bool>> confirmations;
  };
  struct Recorder {
    std::vector<Step> steps;
    std::size_t direction = 0U;
    std::array<FollowEnd, 2> end{FollowEnd::kCapacity, FollowEnd::kCapacity};
  };
  /// The stage's observer, writing into a recorder.
  struct Observer : follow::NoObserver {
    Recorder* recorder = nullptr;
    void begin(PropagationDirection direction) const { recorder->direction = direction == PropagationDirection::kForward ? 0U : 1U; }
    void step(LayerID /*layer*/, const BoundTrackState& state, const CovarianceFactors& factors, std::uint32_t /*rows*/) const
    {
      auto& step = recorder->steps.emplace_back();
      step.direction = recorder->direction;
      step.state = state;
      step.factors = factors;
    }
    void target(LayerID layer, const CandidateWindow& /*window*/) const { recorder->steps.back().targets.push_back(layer); }
    /// The rows the update was made for, its chi2 (infinity where it failed).
    void candidate(LayerID layer, std::uint32_t row, StepOutcome outcome, Scalar chi2) const
    {
      if (outcome == StepOutcome::kUsed || outcome == StepOutcome::kOutsideWindow || outcome == StepOutcome::kHole) {
        return;
      }
      const Scalar seen = outcome == StepOutcome::kNoUpdate ? std::numeric_limits<Scalar>::infinity() : chi2;
      recorder->steps.back().candidates.push_back({.layer = layer, .row = row, .chi2 = seen});
    }
    void lookaheadConfirmation(std::uint32_t row, Scalar /*chi2*/, std::uint32_t next, Scalar /*nextChi2*/) const { recorder->steps.back().confirmations.emplace_back(row, next != ClusterList<kMaxClusters>::kNoCluster); }
    void taken(std::uint32_t row, Scalar chi2, bool kinked) const
    {
      recorder->steps.back().taken = row;
      recorder->steps.back().chi2 = chi2;
      recorder->steps.back().kinked = kinked;
    }
    void end(FollowEnd why) const { recorder->end[recorder->direction] = why; }
  };

  /// What a track's checks read.
  struct Context {
    std::size_t iteration = 0U;
    std::uint32_t track = 0U;
    std::uint64_t particle = 0U;
    const IterationParameters<Det>* parameters = nullptr;
    const TimeFrame<Det>* frame = nullptr;
    const TrackingGeometry<Det>* geometry = nullptr;
    const Material* material = nullptr;
    const Propagator* propagator = nullptr;
    const std::vector<std::uint8_t>* used = nullptr;
  };

  /// The iteration's follow again, with \a BeamLine, \a Greedy,
  /// \a EndDensity and \a BeamEnd as it ran.
  template <bool BeamLine, bool Greedy, bool EndDensity, bool BeamEnd, typename Observation>
  void replay(const Observation& iteration, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    const auto view = std::as_const(*iteration.tracks).view();
    const auto trackCount = static_cast<std::uint32_t>(std::min(iteration.followStart.size(), iteration.proposals.size()));
    // The mask as the follow read it: the tracker's copy, taken with the
    // starting states.
    std::vector<std::uint8_t> used(iteration.followUsed.begin(), iteration.followUsed.end());
    const FollowStage<Det, Material, BeamLine, false, Greedy, EndDensity, BeamEnd> stage{&parameters, &frame, &geometry, &material, propagator, view, used.data(), mWidths.data(), mCells.view(), nullptr, nullptr, parameters.sharedInnermost == 0U, iteration.trackOffset};
    for (std::uint32_t i = 0; i < trackCount; ++i) {
      const auto track = iteration.trackOffset + i;
      if (view.status[track] != PropagationStatus::kOk) {
        continue;
      }
      // The road: the track's rows less the follow's, inward ones first.
      const auto& kept = iteration.proposals[i];
      const std::uint32_t before = iteration.followStartCount[i];
      ClusterList<kMaxClusters> road;
      std::ranges::fill(road.index, ClusterList<kMaxClusters>::kNoCluster);
      for (std::uint32_t at = 0; at < before; ++at) {
        road.index[at] = view.clusters[track].index[kept.inwardCount + at];
      }
      const auto particle = study::goodParticle(mTruth->rows, std::span<const std::uint32_t>{road.index, before});
      if (particle == study::kNoTrackParticle) {
        continue;
      }
      Recorder recorder;
      Proposal proposal;
      std::ranges::fill(proposal.outward, ClusterList<kMaxClusters>::kNoCluster);
      std::ranges::fill(proposal.inward, ClusterList<kMaxClusters>::kNoCluster);
      follow::KinkSteps steps;
      stage.follow(iteration.followStartInner[i], iteration.followStart[i], road, before, proposal, steps, Observer{{}, &recorder});
      const Context context{
        .iteration = iteration.iteration,
        .track = track,
        .particle = particle,
        .parameters = &parameters,
        .frame = &frame,
        .geometry = &geometry,
        .material = &material,
        .propagator = &propagator,
        .used = &used,
      };
      const bool agrees = write(context, road, before, kept, proposal, recorder);
      ++mReplayed;
      mAgree += agrees ? 1U : 0U;
    }
  }

  /// The prediction of \a row from \a state and \a factors in \a direction:
  /// the update's chi2, the two pulls on the row's own plane, the
  /// prediction's two widths and the row's two errors; NaN where it fails.
  [[nodiscard]] std::array<Scalar, 7> predict(const Context& context, const BoundTrackState& state, const CovarianceFactors& factors, std::uint32_t row, PropagationDirection direction) const
  {
    constexpr Scalar kNaN = std::numeric_limits<Scalar>::quiet_NaN();
    std::array<Scalar, 7> result{};
    result.fill(kNaN);
    const auto measurement = measurementAt(*context.frame, row);
    BoundTrackState updated;
    CovarianceFactors updatedFactors;
    Scalar chi2 = 0.F;
    if (follow::updateCandidate<false>(*context.propagator, *context.geometry, *context.material, measurement, direction, state, factors, updated, updatedFactors, chi2, 1.F)) {
      result[0] = chi2;
    }
    BoundTrackState predicted;
    CovarianceFactors predictedFactors;
    if (follow::transportCandidate<false>(*context.propagator, *context.geometry, *context.material, context.geometry->planeOf(measurement), direction, state, factors, predicted, predictedFactors, 1.F)) {
      covarianceFromFactors(predictedFactors, predicted.cov);
      result[1] = (measurement.loc0 - predicted.loc0()) / std::sqrt(predicted.variance(kLoc0) + measurement.varLoc0);
      result[2] = (measurement.loc1 - predicted.loc1()) / std::sqrt(predicted.variance(kLoc1) + measurement.varLoc1);
      result[3] = std::sqrt(predicted.variance(kLoc0));
      result[4] = std::sqrt(predicted.variance(kLoc1));
      result[5] = std::sqrt(measurement.varLoc0);
      result[6] = std::sqrt(measurement.varLoc1);
    }
    return result;
  }

  [[nodiscard]] Scalar distance2(std::uint32_t row) const
  {
    const auto& at = mTruth->rows[row];
    return (at.r * at.r) + (at.z * at.z);
  }

  /// The `S` and `F` lines of a track's replay and its stops; whether the
  /// rows the tracking kept are a prefix of the replay's in both directions.
  bool write(const Context& context, const ClusterList<kMaxClusters>& road, std::uint32_t before, const Proposal& kept, const Proposal& proposal, const Recorder& recorder)
  {
    const auto& truthRows = mTruth->rows;
    const auto owns = [&](std::uint32_t row) { return truthRows[row].carries(context.particle); };
    const auto found = mTruth->particles.find(context.particle);
    bool agrees = true;
    for (std::size_t direction = 0; direction < 2U; ++direction) {
      const bool outward = direction == 0U;
      const std::uint32_t* proposed = outward ? proposal.outward : proposal.inward;
      const std::uint32_t proposedCount = outward ? proposal.outwardCount : proposal.inwardCount;
      const std::uint32_t* keptRows = outward ? kept.outward : kept.inward;
      const std::uint32_t keptCount = outward ? kept.outwardCount : kept.inwardCount;
      bool prefix = keptCount <= proposedCount;
      std::uint32_t keptWrong = 0U;
      for (std::uint32_t i = 0; i < keptCount; ++i) {
        prefix = prefix && proposed[i] == keptRows[i];
        keptWrong += owns(keptRows[i]) ? 0U : 1U;
      }
      agrees = agrees && prefix;
      // The reach so far: the road's outermost (innermost) row, then each
      // taken one.
      Scalar reached = outward ? 0.F : std::numeric_limits<Scalar>::max();
      for (std::uint32_t i = 0; i < before; ++i) {
        reached = outward ? std::max(reached, distance2(road.index[i])) : std::min(reached, distance2(road.index[i]));
      }
      int firstWrong = -1;
      auto firstClass = StepClass::kOwn;
      std::size_t k = 0U;
      for (const auto& step : recorder.steps) {
        if (step.direction != direction) {
          continue;
        }
        int takenLayer = -1;
        bool takenOwn = false;
        std::uint32_t passing = 0U;
        Scalar ownChi2 = -1.F;
        std::uint32_t ownRow = ClusterList<kMaxClusters>::kNoCluster;
        auto stepClass = StepClass::kOwn;
        if (step.taken != ClusterList<kMaxClusters>::kNoCluster) {
          takenLayer = truthRows[step.taken].layer;
          takenOwn = owns(step.taken);
          for (const auto& candidate : step.candidates) {
            if (candidate.layer != takenLayer) {
              continue;
            }
            passing += candidate.chi2 <= context.parameters->clusterChi2Max ? 1U : 0U;
            if (owns(candidate.row) && (ownChi2 < 0.F || candidate.chi2 < ownChi2)) {
              ownChi2 = candidate.chi2;
              ownRow = candidate.row;
            }
          }
          if (!takenOwn) {
            stepClass = classify(step, takenLayer, ownChi2, outward, reached, *context.used, *context.parameters, context.particle);
            if (firstWrong < 0) {
              firstWrong = static_cast<int>(k);
              firstClass = stepClass;
            }
          }
          reached = outward ? std::max(reached, distance2(step.taken)) : std::min(reached, distance2(step.taken));
        }
        if (!mPath.empty()) {
          const auto ranked = std::ranges::find(step.confirmations, step.taken, &std::pair<std::uint32_t, bool>::first);
          const int confirmed = ranked == step.confirmations.end() ? -1 : (ranked->second ? 1 : 0);
          const auto ownRanked = std::ranges::find(step.confirmations, ownRow, &std::pair<std::uint32_t, bool>::first);
          const int ownConfirmed = ownRow == ClusterList<kMaxClusters>::kNoCluster || ownRanked == step.confirmations.end() ? -1 : (ownRanked->second ? 1 : 0);
          mOut << std::format("S {} {} {} {} {} {} {:.4f} {} {} {:.4f} {} {} {} {} {}\n", context.iteration, context.track, direction, k, takenLayer, step.taken == ClusterList<kMaxClusters>::kNoCluster ? 0 : 1, static_cast<double>(step.chi2), takenOwn ? 1 : 0, passing, static_cast<double>(ownChi2), kClassNames[static_cast<std::size_t>(stepClass)], step.kinked ? 1 : 0, step.confirmations.size(), confirmed, ownConfirmed);
          writePredictions(context, step, k);
        }
        ++k;
      }
      std::uint32_t proposedWrong = 0U;
      for (std::uint32_t i = 0; i < proposedCount; ++i) {
        proposedWrong += owns(proposed[i]) ? 0U : 1U;
      }
      const auto end = recorder.end[direction];
      if (k != 0U) {
        const bool beyond = std::ranges::any_of(mRowsOf[context.particle], [&](std::uint32_t row) { return (*context.used)[row] == kFree && (outward ? distance2(row) > reached : distance2(row) < reached); });
        ++mStops[direction][static_cast<std::size_t>(end)][beyond ? 1U : 0U];
      }
      if (!mPath.empty()) {
        mOut << std::format("F {} {} {} {:.4f} {} {} {} {} {} {} {} {} {} {}\n", context.iteration, context.track, context.particle, found != mTruth->particles.end() ? static_cast<double>(found->second.pt) : 0., before, direction, proposedCount, keptCount, keptWrong, prefix ? 1 : 0, firstWrong, kClassNames[static_cast<std::size_t>(firstClass)], kEndNames[static_cast<std::size_t>(end)], proposedWrong);
      }
    }
    return agrees;
  }

  /// The `Q` lines of step \a k: per layer it targeted, the particle's row
  /// there predicted from the step's own state (`predict`), the best by chi2.
  void writePredictions(const Context& context, const Step& step, std::size_t k)
  {
    const auto direction = step.direction == 0U ? PropagationDirection::kForward : PropagationDirection::kBackward;
    for (const auto target : step.targets) {
      std::array<Scalar, 7> best{};
      std::uint32_t bestRow = ClusterList<kMaxClusters>::kNoCluster;
      for (const auto row : mRowsOf[context.particle]) {
        if (mTruth->rows[row].layer != target) {
          continue;
        }
        const auto prediction = predict(context, step.state, step.factors, row, direction);
        if (bestRow == ClusterList<kMaxClusters>::kNoCluster || prediction[0] < best[0]) {
          best = prediction;
          bestRow = row;
        }
      }
      if (bestRow == ClusterList<kMaxClusters>::kNoCluster) {
        continue;
      }
      mOut << std::format("Q {} {} {} {} {} {} {:.4g} {:.4g} {:.4g} {:.4g} {:.4g} {:.4g} {:.4g}\n", context.iteration, context.track, step.direction, k, static_cast<int>(target), (*context.used)[bestRow] == kFree ? 1 : 0, static_cast<double>(best[0]), static_cast<double>(best[1]), static_cast<double>(best[2]), static_cast<double>(best[3]), static_cast<double>(best[4]), static_cast<double>(best[5]), static_cast<double>(best[6]));
    }
  }

  /// The class of a step that took a row not of the particle (`StepClass`).
  [[nodiscard]] StepClass classify(const Step& step, int takenLayer, Scalar ownChi2, bool outward, Scalar reached, const std::vector<std::uint8_t>& used, const IterationParameters<Det>& parameters, std::uint64_t particle) const
  {
    if (!(ownChi2 < 0.F)) {
      return ownChi2 <= parameters.clusterChi2Max ? StepClass::kSwap : StepClass::kRefused;
    }
    // The targets searched up to the taken one, in their order.
    const auto upTo = std::ranges::find(step.targets, static_cast<LayerID>(takenLayer));
    bool ownFreeOnTaken = false;
    bool ownHeldOnTaken = false;
    bool ownFreeNearer = false;
    bool ownFreeBeyond = false;
    for (const auto row : mRowsOf.at(particle)) {
      if (outward ? !(distance2(row) > reached) : !(distance2(row) < reached)) {
        continue;
      }
      const LayerID layer = mTruth->rows[row].layer;
      const bool free = used[row] == kFree;
      ownFreeBeyond = ownFreeBeyond || free;
      if (layer == takenLayer) {
        (free ? ownFreeOnTaken : ownHeldOnTaken) = true;
      } else if (free && std::find(step.targets.begin(), upTo, layer) != upTo) {
        ownFreeNearer = true;
      }
    }
    if (ownFreeNearer) {
      return StepClass::kRefusedNearer;
    }
    if (ownFreeOnTaken) {
      return StepClass::kNotInWindow;
    }
    if (ownHeldOnTaken) {
      return StepClass::kHeld;
    }
    return ownFreeBeyond ? StepClass::kFilledHole : StepClass::kPastEnd;
  }

  std::string mPath;
  std::ofstream mOut;
  const study::TrackingTruth* mTruth = nullptr;
  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> mRowsOf;
  SensorCells<Det> mCells;
  std::array<std::uint32_t, static_cast<std::size_t>(Detector<Det>::nLayers)> mWidths{};
  /// Per direction and end, the follows whose particle had no free row
  /// further along, and those that had one.
  std::array<std::array<std::array<std::uint64_t, 2>, 4>, 2> mStops{};
  std::uint64_t mReplayed = 0U;
  std::uint64_t mAgree = 0U;
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_FOLLOWANATOMY_H_
