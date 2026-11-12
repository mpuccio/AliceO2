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

/// \file CompleteAnatomy.h
/// \brief The bench's anatomy of the findable on every layer of a mask:
///        found whole, short and why, not found and why, against the
///        decisions the tracker actually made, and the particles' true kinks
///        from the MC hits.
///
/// The class is the physical primaries above the score's floor with a
/// measurement on every layer of the mask (`--findable-all-layers`),
/// narrowed by `ClassCuts` (`--anatomy-*`). A particle's rows are taken in
/// the order it crossed them: with the MC hits (`--hits`), by their time of
/// flight, its path order; without them, or for a particle
/// with a row the hits do not have, by the distance from its production
/// vertex, which grows along a first outward pass but falls on a looper's
/// way back.
///
/// Its parts: each class particle's fate after every iteration (complete,
/// short and why, not found and why, with the graph audit of a road never
/// made; `--fate-out`); the seeding funnel, the record losses and the
/// dropped records' follow replayed (`printSeeding`); the beam-line test
/// and the follow's losses against the MC (`--anatomy-beam-line`,
/// `--anatomy-follow-loss`); the path order and the true kinks, the truth
/// walk (`--truth-walk`), the inward refit's states (`--inward-history`)
/// and the steps against the MC (`--step-anatomy`), with the hits; the
/// fakes' transitions, the follow's end evidence and decisions
/// (`--end-decisions`, `--end-edits`), the own steps, the hole calls and the
/// road ends; the clutter test of the follow's steps with the MC
/// kinematics (`--kine`) and who blocks the not found.
///
/// For the tracks an iteration made, its follow is done again, outward and
/// inward: every track's proposal by the follow's own code
/// (`FollowStage::follow`) from the refit of its road, on the used mask as
/// the follow found it, with an observer on every decision; then the
/// bidding and the refit's acceptance, as the follow makes them (not the
/// refit's prefix second chance, its trim of an ambiguous last row or the
/// crossing rule). What that gives is checked against what the
/// iteration kept (`fidelity`), and a short track's stop is the decision
/// that made it: its particle's next row lost in the bidding, dropped with
/// a kink step's witness, dropped with the refit (failed or above its
/// chi2), or, at the proposal step, never queried, held, outside the
/// window, crossing a layer, without an update, above the chi2 gate, or
/// passing with another row taken (by the lookahead or not); where the step
/// had the kink allowance, what its kink step did with the row instead: no
/// budget left, refused under the scaled gate, passing it but unconfirmed,
/// or another confirmed row taken. With followEarlier the earlier tracks'
/// follow before the iteration's seeds is replayed the same way, from the
/// tracks and the mask as the previous iteration left them.

#ifndef GTRK_STUDY_COMPLETEANATOMY_H_
#define GTRK_STUDY_COMPLETEANATOMY_H_

#include "GenericTracking/Study/Alice/Hits.h"
#include "GenericTracking/Tracking/BeamConstraint.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/NeighbourCut.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/TrackFitter.h"
#include "GenericTracking/Tracking/TripletCut.h"
#include "GenericTracking/Study/Validation.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/StageFunnel.h"
#include "GenericTracking/Study/Statistics.h"
#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Study/TrackingInput.h"

#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/O2DatabasePDG.h"
#include "Steer/MCKinematicsReader.h"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <numbers>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

template <DetectorLayout Det>
// NOLINTNEXTLINE(clang-analyzer-optin.performance.Padding): a study's members, ordered by what they hold
class CompleteAnatomy
{
 public:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  static constexpr std::size_t kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);

  /// Why a short track is short: its particle's next row inward, the first
  /// layer of the mask it misses.
  enum Short : std::uint8_t {
    kGapInside,       ///< a layer missing between the track's own
    kLostBid,         ///< proposed, lost in the bidding on that row or one before
    kWitnessLost,     ///< proposed and won, dropped with a kink step's witness lost in the bidding
    kRefitFailed,     ///< proposed and won, dropped with the extended refit: the fit failed
    kRefitChi2,       ///< proposed and won, dropped with the extended refit: above its chi2
    kNotQueried,      ///< the step's grid query never visited it
    kUsed,            ///< held by a track when the step ran
    kWindow,          ///< visited, outside the step's window
    kHole,            ///< crossing another layer
    kNoUpdate,        ///< no usable update
    kChi2,            ///< above the chi2 gate
    kOtherTaken,      ///< passing, another row taken
    kLookaheadOther,  ///< passing, the lookahead took another row
    kKinkBudget,      ///< refused by the step, its kink step not made: the budget spent
    kKinkRefused,     ///< refused by the step and by its kink step's scaled gate (or window)
    kKinkUnconfirmed, ///< refused by the step, passing its kink step's, nothing confirming it
    kKinkOther,       ///< refused by the step, passing its kink step's, confirmed, another row taken
    kNoStep,          ///< the proposal never stepped from the layer outside it
    kShortReasons
  };
  static constexpr std::array<const char*, kShortReasons> kShortNames{"gap inside", "lost the bid", "dropped with a kink step's witness", "refit failed", "refit chi2", "not queried", "held", "outside the window", "crossing a layer", "no update", "chi2", "passing, another row taken", "passing, the lookahead took another", "kink budget spent", "kink step refused it", "kink step unconfirmed", "kink step took another", "no step from the layer outside"};
  /// Who holds a particle's own row at the end (the short tracks that lost
  /// the bid): an earlier iteration's track, or this one's.
  enum Holder : std::uint8_t { kNoOwnRow = 0U,
                               kRowFree = 1U,
                               kHeldByFake = 2U,
                               kHeldByOther = 3U,
                               kHeldByOwn = 4U,
                               kThisIteration = 3U,
                               kHolderKinds = 8U };
  static constexpr std::array<const char*, kHolderKinds> kHolderNames{"no own row", "free", "an earlier fake", "an earlier track of another particle", "an earlier track of its own", "a fake of this iteration", "a track of another particle of this iteration", "a track of its own of this iteration"};

  /// Why a class particle has no good track (`missedReason`); the values
  /// are the `--fate-out` codes.
  enum Missed : std::uint8_t { kLostToFake,
                               kTakenByOther,
                               kRoadChi2,
                               kRoadPropagation,
                               kFloor,
                               kBeamLine,
                               kNoRoad,
                               kTooFewRows,
                               kRoadTime,
                               kMissedReasons };
  static constexpr std::array<const char*, kMissedReasons> kMissedNames{"road rows held by a fake", "held by another particle's track", "road fit chi2", "road fit propagation", "pT floor", "beam line", "passing: no road made", "fewer road rows than the pass's clusters", "road rows share no time"};

  /// The first predicate of the graph that refuses a "no road made"
  /// particle's true outer rows, the triplets inner first, then the
  /// neighbours; a triplet's with the kink allowance
  /// where the iteration has one and the nominal cuts refuse.
  enum Graph : std::uint8_t {
    kGraphLambda,
    kGraphTime,
    kGraphCurvature,
    kGraphFitFailed,
    kGraphChi2,
    kGraphOwnMomentum,
    kGraphBothKinked,
    kGraphNeighbourCurvature,
    kGraphNeighbourTime,
    kGraphNeighbourFitFailed,
    kGraphNeighbourChi2,
    kGraphLinkRange,    ///< a row outside the z range its link admits sources from
    kGraphLinkNotTaken, ///< the layout's link between two consecutive true rows is not one the iteration takes
    kGraphStart,        ///< a start-layer mask: the chain's innermost row on a layer roads may not start on
    kGraphPairTime,     ///< the tracklet stage's own tests, as `pairCutLoss` orders them
    kGraphPairTurning,
    kGraphPairCurvature,
    kGraphPairLongitudinal,
    kGraphPairWindow, ///< the target outside the cells the source's query visits
    kGraphPairHole,   ///< no direct pair: the segment crosses another layer's silicon (LinkHoles.h)
    kGraphPasses,
    kGraphReasons
  };
  static constexpr std::array<const char*, kGraphReasons> kGraphNames{"triplet: dip window", "triplet: time", "triplet: curvature", "triplet: fit failed", "triplet: fit chi2", "triplet: own-momentum chi2", "neighbour: both triplets kinked", "neighbour: curvature agreement", "neighbour: time", "neighbour: fit failed", "neighbour: fit chi2", "link: source outside its z range", "link: not taken by the iteration", "start: the chain's innermost row on a layer roads may not start on", "pair: time", "pair: turning point", "pair: curvature", "pair: longitudinal", "pair: outside the query's cells", "pair: crosses another layer", "every predicate passes"};

  /// The class restricted to a species and optionally the QA's findable.
  struct ClassCuts {
    std::int32_t pdg = 0; ///< |PDG code| of the class's species; 0, every species
    /// Only the findable by these settings (`isFindable`, the QA's rule:
    /// a chain of `minLayers`); none, no such cut.
    const study::ComparisonSettings* findable = nullptr;
  };

  /// \a roadLayers: each iteration's `layerMask`, 0 for every layer (kept
  /// with every bit set, so that its readers test bits only).
  CompleteAnatomy(const study::TrackingTruth& truth, std::uint64_t mask, double minPtGeV, std::vector<std::uint64_t> roadLayers, ClassCuts cuts = {})
    : mTruth(&truth),
      mMask(mask),
      mMinPtGeV(minPtGeV),
      mRoadLayers([&] {
        std::ranges::replace(roadLayers, std::uint64_t{0}, ~std::uint64_t{0});
        return std::move(roadLayers);
      }()),
      mRowsOf(study::rowsOfParticles(truth.rows)),
      mWalkIteration([&] {
        for (std::size_t i = 0; i < mRoadLayers.size(); ++i) {
          if ((mRoadLayers[i] & mask) == mask) {
            return i;
          }
        }
        return std::size_t{0};
      }())
  {
    const auto chains = cuts.findable != nullptr ? study::particleLayers(truth.rows, *cuts.findable) : std::unordered_map<std::uint64_t, study::ParticleLayers>{};
    for (auto& [particle, rows] : mRowsOf) {
      const auto found = truth.particles.find(particle);
      std::uint64_t layers = 0U;
      for (const auto row : rows) {
        layers |= std::uint64_t{1} << truth.rows[row].layer;
      }
      // The class's reach in pT: physical primaries in |eta| < 1, every one and
      // those on every layer of the mask, in steps of 50 MeV to 1 GeV.
      if (found != truth.particles.end() && found->second.origin == study::ParticleOrigin::kPhysicalPrimary && std::abs(found->second.eta) < 1.F) {
        const auto bin = std::min<std::size_t>(static_cast<std::size_t>(found->second.pt / 0.05F), mReach.size() - 1U);
        ++mReach[bin][0];
        if ((layers & mMask) == mMask) {
          ++mReach[bin][1];
          mLowestPt = std::min(mLowestPt, static_cast<double>(found->second.pt));
        }
      }
      if (found != truth.particles.end() && found->second.origin == study::ParticleOrigin::kPhysicalPrimary && static_cast<double>(found->second.pt) >= minPtGeV && (layers & mMask) == mMask && (cuts.pdg == 0 || std::abs(found->second.pdg) == cuts.pdg) && (cuts.findable == nullptr || (chains.contains(particle) && study::isFindable(found->second.pt, chains.at(particle), *cuts.findable)))) {
        mIndex[particle] = static_cast<std::uint32_t>(mFindable.size());
        mFindable.push_back(particle);
      }
      sortByTraversal(particle, rows);
    }
    mFate.assign(mFindable.size(), kNotFound);
    mReason.assign(mFindable.size(), 0U);
    mGraph.assign(mFindable.size(), kGraphReasons);
    mStop.assign(mFindable.size(), kShortReasons);
    mStopDetail.assign(mFindable.size(), 0U);
    mStopLayer.assign(mFindable.size(), kInvalidLayerID);
    mStopTaken.assign(mFindable.size(), kNoTaken);
    mOutStop.assign(mFindable.size(), kShortReasons);
    mOutStopLayer.assign(mFindable.size(), kInvalidLayerID);
    mOutHolder.assign(mFindable.size(), kHolderKinds);
    mOutRefitWrong.assign(mFindable.size(), 0U);
    mOutRefitRows.assign(mFindable.size(), 0U);
    mOutRefitRatio.assign(mFindable.size(), 0);
    mRefitFailures.clear();
    mRefitReplay.clear();
    mOutRefitKink.assign(mFindable.size(), 3U);
    mOutStopTaken.assign(mFindable.size(), kNoTaken);
    mInNext.assign(mFindable.size(), 255U);
    mOutNext.assign(mFindable.size(), 255U);
    mInDropKink.assign(mFindable.size(), 255U);
    mOutDropKink.assign(mFindable.size(), 255U);
    mOutStepKink.assign(mFindable.size(), 255U);
    mInStepKink.assign(mFindable.size(), 255U);
    mLayersHeld.assign(mFindable.size(), 0U);
    mMissedKink.assign(mFindable.size(), 3U);
    mInChi2.assign(mFindable.size(), 0.F);
    mOutChi2.assign(mFindable.size(), 0.F);
    mInNextChi2.assign(mFindable.size(), 0.F);
    mOutNextChi2.assign(mFindable.size(), 0.F);
    mRecordLoss.assign(mFindable.size(), std::string{});
    mFurthest.assign(mFindable.size(), kSeedNotEligible);
  }

  /// The CSV of the inward refit's states: each state of
  /// the inward pass -- the seed, then each filtered state -- of the refits
  /// failing inward and of a control of good low-momentum refits, against
  /// the MC truth at the row.
  void openInwardHistory(const std::string& csv)
  {
    mInwardCsv = std::make_unique<std::ofstream>(csv);
    if (!*mInwardCsv) {
      throw std::runtime_error("cannot write " + csv);
    }
    *mInwardCsv << "kind,particle,pt,pdg,n,status,failStep,k,row,layer,leg,p,pTrue,phi,phiTrue,theta,thetaTrue,sigQOverPRel,sigPhi,chi2,loc0,loc1,loc0True,loc1True,"
                   "qOverP,charge,stepRoot\n";
  }

  /// The CSV of every follow walk's end, written after each iteration
  /// (`endDecisions`).
  void openEndDecisions(const std::string& csv)
  {
    mEndCsv = std::make_unique<std::ofstream>(csv);
    if (!*mEndCsv) {
      throw std::runtime_error("cannot write " + csv);
    }
    *mEndCsv << "track,iteration,direction,particle,pt,end,crossed,links,unconfirmedDrop,walked,bidKept,kept,lastRow,lastLayer,lastChi2,lastKink,lastKind,lastBeyond,"
                "droppedRow,droppedLayer,droppedChi2,droppedKind,droppedBeyond,survives,refit,agrees,trackGood,trackCount,trackLayers,ownLayers,"
                "endWindows,endSeen,endRefused,lastStepChi2,lastSqrtDet,lastRhoWindow,lastRhoLayer,lastReach,lastOverlapPS,lastBoxPS,lastRowPS,lastRest,"
                "lastBeamWith,lastBeamWithout,lastFloatDet,lastFactorDet,"
                "droppedStepChi2,droppedSqrtDet,droppedRhoWindow,droppedRhoLayer,droppedReach,droppedOverlapPS,droppedBoxPS,droppedRowPS,droppedRest,"
                "droppedBeamWith,droppedBeamWithout,droppedFloatDet,droppedFactorDet,"
                "endBeamChi2,endBeamPasses,refRefused,refStep,refRow,refKind,refLast,refAlt,refAltKind,refNew,refNewKind,lastRefChi2,lastRefPasses,refInvalid,"
                "lastSearch,edit,editCount,editFit,editKept,editLayers\n";
  }

  /// The CSV of each track's end edits (`endDecisions`): the
  /// track as it is, with its rows, and each truth-free edit of its ends.
  void openEndEdits(const std::string& csv)
  {
    mEditCsv = std::make_unique<std::ofstream>(csv);
    if (!*mEditCsv) {
      throw std::runtime_error("cannot write " + csv);
    }
    *mEditCsv << "track,iteration,in,out,count,fit,kept,good,majority,layers,rows\n";
  }

  /// The MC kinematics (\a path, a `<prefix>_Kine.root` such as
  /// `sgn_Kine.root`), each MC track's mother, for the clutter test's
  /// lineage: every event read now through O2's reader, as
  /// the test looks the mothers up from several threads.
  void openKine(const std::string& path)
  {
    constexpr std::string_view kSuffix = "_Kine.root";
    if (!path.ends_with(kSuffix) || !std::filesystem::exists(path)) {
      throw std::runtime_error("--kine: no file " + path + " (a <prefix>_Kine.root)");
    }
    o2::steer::MCKinematicsReader reader(path.substr(0, path.size() - kSuffix.size()), o2::steer::MCKinematicsReader::Mode::kMCKine);
    mMothers.resize(reader.getNEvents(0));
    if (mMothers.empty()) {
      throw std::runtime_error("--kine: no MC events in " + path);
    }
    for (std::size_t event = 0; event < mMothers.size(); ++event) {
      const auto& tracks = reader.getTracks(0, static_cast<int>(event));
      mMothers[event].reserve(tracks.size());
      for (const auto& track : tracks) {
        mMothers[event].push_back(track.getMotherTrackId());
      }
      // Only the mothers are kept: the reader's copy of the event goes.
      reader.releaseTracksForSourceAndEvent(0, static_cast<int>(event));
    }
  }

  /// The truth walk on iteration \a iteration instead of the first holding
  /// the class's mask (for the inward steps into the VD of the iterations
  /// starting outside it).
  void setWalkIteration(std::size_t iteration) { mWalkIteration = iteration; }

  /// The CSV of the truth walks' steps:
  /// each step of every findable's walk along its own rows, inward from
  /// the true road's refit and outward from its outward state, whatever the
  /// step's chi2, against the cluster and the MC crossing.
  void openTruthWalk(const std::string& csv)
  {
    mWalkCsv = std::make_unique<std::ofstream>(csv);
    if (!*mWalkCsv) {
      throw std::runtime_error("cannot write " + csv);
    }
    *mWalkCsv << "direction,iteration,particle,pt,eta,pdg,fromLayer,layer,status,r0,r1,i00,i11,p00,p11,chi2,c0,c1,kinkLocal,kinkUpstream,onLayer\n";
  }

  /// The fate of every particle of the class after the last iteration
  /// observed: one CSV line each, for a comparison between runs.
  void writeFates(const std::string& csv) const
  {
    std::ofstream out(csv);
    if (!out) {
      throw std::runtime_error("cannot write " + csv);
    }
    out << "particle,pt,eta,pdg,fate,layers,stop,stopLayer,outStop,outStopLayer,reason,crossed,rowsBack,overlaps,seedStage\n";
    for (std::size_t index = 0; index < mFindable.size(); ++index) {
      const auto found = mTruth->particles.find(mFindable[index]);
      const float pt = found != mTruth->particles.end() ? found->second.pt : 0.F;
      const float eta = found != mTruth->particles.end() ? found->second.eta : 0.F;
      const int pdg = found != mTruth->particles.end() ? found->second.pdg : 0;
      // The layers it crossed, its rows on later legs and its
      // overlaps (-1 without a timed path), its last iteration's SeedStage.
      std::uint64_t crossed = 0U;
      int back = -1;
      int overlaps = -1;
      const auto& own = mRowsOf.at(mFindable[index]);
      for (const auto row : own) {
        crossed |= std::uint64_t{1} << mTruth->rows[row].layer;
      }
      if (!mLegs.empty() && mTimed[index] != 0U) {
        back = 0;
        overlaps = 0;
        for (std::size_t k = 0; k < own.size(); ++k) {
          back += mLegs[index][k] != 0U ? 1 : 0;
          for (std::size_t j = 0; j < k; ++j) {
            overlaps += mLegs[index][j] == mLegs[index][k] && mTruth->rows[own[j]].layer == mTruth->rows[own[k]].layer ? 1 : 0;
          }
        }
      }
      const int seedStage = mSeedStage.empty() ? -1 : static_cast<int>(mSeedStage.back()[index]);
      out << std::format("{},{:.4f},{:.3f},{},{},{},{},{},{},{},{},{},{},{},{}\n", mFindable[index], static_cast<double>(pt), static_cast<double>(eta), pdg, static_cast<int>(mFate[index]), mLayersHeld[index], mFate[index] == kShortFate ? static_cast<int>(mStop[index]) : -1, mFate[index] == kShortFate && mStopLayer[index] < kLayers ? static_cast<int>(mStopLayer[index]) : -1, mFate[index] == kShortFate && mStop[index] == kGapInside ? static_cast<int>(mOutStop[index]) : -1, mFate[index] == kShortFate && mOutStopLayer[index] < kLayers ? static_cast<int>(mOutStopLayer[index]) : -1, mFate[index] == kNotFound ? static_cast<int>(mReason[index]) : -1, crossed, back, overlaps, seedStage);
    }
  }

  /// The CSV of the replayed steps, if \a csv is given; the MC hits, if
  /// \a hits is, which the steps and the true kinks are held to.
  void open(const std::string& csv, const std::string& hits, const std::string& discHits = {})
  {
    if (!csv.empty()) {
      mCsv = std::make_unique<std::ofstream>(csv);
      if (!*mCsv) {
        throw std::runtime_error("cannot write " + csv);
      }
      *mCsv << "kind,particle,iteration,layer,pt,pdg,labels,rowsOnLayer,chi2,r0,r1,v00,v01,v11,p00,p01,p11,hit,c0,c1,q0,q1,dphi,dtheta,pphi,ptheta,hitP,stateP,"
               "startHit,sdphi,sdtheta,sphi,stheta,staveKink,stavePhiKink,staveTheta0,staveX2x0,ref,chi2NoBeam,qn0,qn1,pn00,pn11,chi2Joint,qj0,qj1,pj00,pj11\n";
    }
    if (!hits.empty()) {
      if constexpr (std::is_same_v<Det, detectors::TRKFT3>) {
        mHits = std::make_unique<study::HitFile>(hits, study::Chips::kTRK, discHits); // FT3's hits beside TRK's (the disc crossings)
      } else if constexpr (std::is_same_v<Det, detectors::ITS>) {
        mHits = std::make_unique<study::HitFile>(hits, study::Chips::kITS);
      } else if constexpr (std::is_same_v<Det, detectors::MFT>) {
        mHits = std::make_unique<study::HitFile>(hits, study::Chips::kMFT);
      } else {
        throw std::runtime_error("--hits: no MC hits for this detector");
      }
      // Every event read now: the truth walks look the segments up from
      // several threads, which must find them and never read.
      for (int event = 0; event < mHits->nEvents(); ++event) {
        static_cast<void>(mHits->event(event));
      }
      // Each MC track's segments, for `segmentOf` (the path order looks up
      // every row of every particle).
      for (int event = 0; event < mHits->nEvents(); ++event) {
        for (const auto& segment : mHits->event(event)) {
          mSegmentsOf[segmentKey(event, segment.track)].push_back(&segment);
        }
      }
    }
  }

  template <typename Observation, typename Material>
  void observe(const Observation& iteration, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    Phases phases;
    const auto now = [] { return std::chrono::steady_clock::now(); };
    const auto seconds = [](auto from, auto to) { return std::chrono::duration<double>(to - from).count(); };
    auto clock = now();
    if (mHits && !mPathOrdered) {
      pathOrder(frame);
      mPathOrdered = true;
    }
    if (mHits && !mKinksDone) {
      trueKinks(frame, material, propagator);
      mKinksDone = true;
      phases.kinks = seconds(clock, now());
      clock = now();
    }
    if (mCells.empty()) {
      mCells = SensorCells<Det>(geometry);
    }
    if (mMaxWidthPS.empty()) {
      mMaxWidthPS.assign(kLayers, 0U);
      const auto columns = frame.measurements();
      for (std::uint32_t row = 0; row < mTruth->rows.size(); ++row) {
        const auto layer = mTruth->rows[row].layer;
        mMaxWidthPS[layer] = std::max(mMaxWidthPS[layer], columns.widthPS[row]);
        auto& bounds = mBounds[layer];
        bounds.varZ = std::max(bounds.varZ, columns.varZ[row]);
        bounds.varR = std::max(bounds.varR, columns.varR[row]);
        bounds.varRPhi = std::max(bounds.varRPhi, columns.varRPhi[row]);
      }
    }
    // The truth walk once, on the first iteration whose road layers hold
    // the class's mask -- its road is then the class's whole chain, the
    // barrel-to-the-rest leg included -- else on the first.
    if (mHits && !mWalkDone && iteration.iteration == mWalkIteration) {
      truthWalk(iteration.iteration, parameters, frame, geometry, material, propagator);
      mWalkDone = true;
      phases.walk = seconds(clock, now());
      clock = now();
    }
    const auto view = iteration.tracks->view();
    const auto tracks = iteration.tracks->size();
    std::vector<std::uint32_t> holder(mTruth->rows.size(), kNoTrack);
    std::vector<std::uint64_t> particles(tracks);
    std::unordered_map<std::uint64_t, std::uint32_t> goodTrack;
    for (std::size_t track = 0; track < tracks; ++track) {
      const std::span<const std::uint32_t> rows{view.clusters[track].index, view.count[track]};
      particles[track] = study::goodParticle(mTruth->rows, rows);
      for (const auto row : rows) {
        holder[row] = static_cast<std::uint32_t>(track);
      }
      if (particles[track] != study::kNoTrackParticle) {
        goodTrack.emplace(particles[track], static_cast<std::uint32_t>(track));
      }
    }
    // The iteration's follow of its own tracks; with followEarlier the
    // earlier tracks' before it, from the previous observation's snapshot.
    // The mask as the follow read it: the tracker's copy, taken with the
    // starting states. An iteration without a follow has no copy
    // and no proposals: the mask after it, which nothing replays.
    std::vector<std::uint8_t> ownUsed = iteration.followUsed.empty() ? std::vector<std::uint8_t>(iteration.used.begin(), iteration.used.end()) : std::vector<std::uint8_t>(iteration.followUsed.begin(), iteration.followUsed.end());
    mFollowUsedNow = ownUsed;
    const auto replayed = replayExtension(iteration, iteration.trackOffset, static_cast<std::uint32_t>(tracks), iteration.proposals, iteration.kinkSteps, ownUsed, false, parameters, frame, geometry, material, propagator);
    std::vector<Replay> recovered;
    if (!iteration.earlierProposals.empty() && iteration.earlierProposals.size() == mSnapshot.size()) {
      recovered = replayExtension(iteration, 0U, iteration.trackOffset, iteration.earlierProposals, iteration.earlierKinkSteps, mSnapshotUsed, true, parameters, frame, geometry, material, propagator);
    }
    phases.replay = seconds(clock, now());
    clock = now();
    const auto replayOf = [&](std::uint32_t track) -> const Replay* {
      if (track >= iteration.trackOffset) {
        return track - iteration.trackOffset < replayed.size() ? &replayed[track - iteration.trackOffset] : nullptr;
      }
      return track < recovered.size() ? &recovered[track] : nullptr;
    };

    // The fate is decided anew each iteration; a short track's diagnosis is
    // made once, when the iteration that made it is observed (the only one
    // replaying it without followEarlier), and kept: the fields below are
    // reset per track where it is replayed, not here (else a track made by an
    // earlier pass reads "not replayed").
    std::fill(mFate.begin(), mFate.end(), kNotFound);
    mLastGate = parameters.clusterChi2Max;
    mLayersFound = {};
    // The graph audit's histograms are the iteration's, like the fates.
    mAuditFloor = iteration.iteration < mFloors.size() ? mFloors[iteration.iteration] : 0.;
    if (mTrackOffsets.size() <= iteration.iteration) {
      // An iteration not observed made no tracks: the last
      // offset, which keeps the list sorted for the holders' lookup.
      mTrackOffsets.resize(iteration.iteration + 1U, mTrackOffsets.empty() ? 0U : mTrackOffsets.back());
    }
    mTrackOffsets[iteration.iteration] = iteration.trackOffset;
    mFakeHolders.clear();
    mCurvatureTriplets.clear();
    mLinkRangeRefusals.clear();
    mPairRefusals.clear();
    mRoadFitRefusals.clear();
    for (std::uint32_t index = 0; index < mFindable.size(); ++index) {
      const auto particle = mFindable[index];
      const auto& own = mRowsOf.at(particle);
      const auto good = goodTrack.find(particle);
      mLayersHeld[index] = 0U;
      if (good == goodTrack.end()) {
        mReason[index] = missedReason(own, holder, particles, iteration.iteration, parameters, frame, geometry, material, propagator);
        if (mReason[index] == kLostToFake) {
          // The fake holding the first of its road rows: the iteration that
          // made it (whose fakes take the rows).
          for (const auto row : roadRows(own, iteration.iteration)) {
            if (holder[row] != kNoTrack && particles[holder[row]] == study::kNoTrackParticle) {
              const auto made = static_cast<std::size_t>(std::upper_bound(mTrackOffsets.begin(), mTrackOffsets.end(), holder[row]) - mTrackOffsets.begin());
              const double pt = truePt(mTruth->particles.find(particle));
              ++mFakeHolders[std::format("pT {} the floor, a fake of iteration {}", pt >= mAuditFloor ? "at or above" : "below", made == 0U ? 0U : made - 1U)];
              break;
            }
          }
        }
        if (mReason[index] == kNoRoad || mReason[index] == kRoadTime) {
          mGraph[index] = graphAudit(roadRows(own, iteration.iteration), parameters, frame, particle);
        } else if (mReason[index] == kRoadChi2 || mReason[index] == kBeamLine) {
          // The road refused whole: its largest true kink over the road's
          // layers, for the floor split.
          const auto road = roadRows(own, iteration.iteration);
          const auto kink = largestKink(particle, mTruth->rows[road.front()].layer, mTruth->rows[road.back()].layer);
          mMissedKink[index] = kinkBin(kink);
          const auto found = mTruth->particles.find(particle);
          const double floor = iteration.iteration < mFloors.size() ? mFloors[iteration.iteration] : 0.;
          if (mReason[index] == kRoadChi2 && found != mTruth->particles.end() && static_cast<double>(found->second.pt) >= floor) {
            describeRefusedRoad(found->second, kink, road, parameters, frame, geometry, material, propagator);
          }
        }
        continue;
      }
      const auto track = good->second;
      std::uint64_t layers = 0U;
      LayerID innermost = std::numeric_limits<LayerID>::max();
      for (std::uint8_t i = 0; i < view.count[track]; ++i) {
        const auto layer = mTruth->rows[view.clusters[track].index[i]].layer;
        layers |= std::uint64_t{1} << layer;
        innermost = std::min(innermost, layer);
      }
      mLayersHeld[index] = layers;
      ++mLayersFound[static_cast<std::size_t>(std::popcount(layers & mMask))];
      const auto missing = mMask & ~layers;
      if (missing == 0U) {
        mFate[index] = kComplete;
        // The control of the refit's increments: a complete track this
        // iteration followed outward, its refit kept.
        if (const Replay* replay = replayOf(track); replay != nullptr && !replay->outward.empty() && replay->outward.front().kept != 0U && replay->outward.front().refit == kRefitKept) {
          addIncrements(mKeptIncrements, replay->outward.front());
        }
        // The inward states' control: a good refit of a
        // low-momentum track, the follow's rows in it.
        if (const Replay* replay = replayOf(track); mHits && replay != nullptr && replay->refit == kRefitKept && replay->refitCount != 0U && pathPtBin(particle) < 3U) {
          inwardStudy(kInwardControl, particle, replay->refitRows, parameters, frame, geometry, material, propagator);
        }
        continue;
      }
      mFate[index] = kShortFate;
      const Replay* replay = replayOf(track);
      if (replay == nullptr) {
        // Diagnosed when its iteration was observed; that diagnosis stands.
        continue;
      }
      mStopTaken[index] = kNoTaken;
      mOutStop[index] = kShortReasons;
      mOutStopLayer[index] = kInvalidLayerID;
      mOutHolder[index] = kHolderKinds;
      mOutStopTaken[index] = kNoTaken;
      mInNext[index] = 255U;
      mOutNext[index] = 255U;
      mInDropKink[index] = 255U;
      mOutDropKink[index] = 255U;
      mOutStepKink[index] = 255U;
      mInStepKink[index] = 255U;
      if ((missing >> innermost) != 0U) {
        mStop[index] = kGapInside;
        mStopLayer[index] = kInvalidLayerID;
        // Outward: the first layer of the mask it misses
        // outward of its innermost, decided at the outward walk's step from
        // the track's layer nearest below it.
        if (!replay->outward.empty()) {
          const auto gap = static_cast<LayerID>(std::countr_zero(missing & (~std::uint64_t{0} << innermost)));
          const auto below = layers & ((std::uint64_t{1} << gap) - 1U);
          const auto outside = static_cast<LayerID>(63 - std::countl_zero(below));
          const auto ownRow = ownRowOn(particle, gap);
          mOutStop[index] = stopReason(replay->outward.front(), ownRow, outside, particle, mStopDetail[index]);
          mOutStepKink[index] = kinkBin(largestKink(particle, outside, gap));
          mOutStopLayer[index] = gap;
          if (mOutStop[index] == kWindow || mOutStop[index] == kNotQueried) {
            const auto& out = replay->outward.front();
            const auto step = std::ranges::find(out.steps, outside, &Step::layer);
            if (step != out.steps.end()) {
              std::array<float, 4> offsets{};
              const auto kind = mOutStop[index] == kWindow ? 0U : 1U;
              ++mOutWindowAudit[kind][windowAudit(*step, ownRow, gap, parameters, frame, geometry, material, propagator, PropagationDirection::kForward, &offsets)];
              mOutWindowOffsets[kind].push_back(offsets);
            }
          }
          if (mOutStop[index] == kRefitFailed) {
            const auto& out = replay->outward.front();
            ++mRefitFailures[std::format("{} / pass {} / layer {}", static_cast<int>(out.refitStatus), static_cast<int>(out.refitFailedPass), static_cast<int>(out.refitFailedLayer))];
            ++mRefitFailures[std::format("trimmed: {} of {} outward rows kept, chi2 {}", out.refitTrimmedKept, out.kept, out.refitTrimmedChi2 ? "under" : "over")];
            std::uint32_t wrong = 0U;
            for (std::uint32_t i = 0; i < out.kept; ++i) {
              wrong += mTruth->rows[out.raw[i]].carries(particle) ? 0U : 1U;
            }
            const char* rows = wrong == 0U ? "own rows" : (wrong == 1U ? "1 wrong row" : "2+ wrong rows");
            ++mRefitReplay[std::format("{} | 1. all", rows)];
            if (mHits && out.refitFailedPass == FitPass::kInward) {
              inwardStudy(wrong == 0U ? kInwardOwn : kInwardWrong, particle, replay->refitRows, parameters, frame, geometry, material, propagator);
            }
            const auto& step = out.failedStep;
            if (!step.replayed) {
              ++mRefitReplay[std::format("{} | 2. not replayed", rows)];
            } else {
              ++mRefitReplay[std::format("{} | 2. pass {}: status {}", rows, static_cast<int>(out.refitFailedPass), static_cast<int>(step.state))];
              const auto missBin = [](float miss) {
                constexpr std::array<float, 4> kEdges{0.01F, 0.1F, 0.5F, 2.F};
                constexpr std::array<const char*, 5> kNames{"< 0.01", "0.01-0.1", "0.1-0.5", "0.5-2", ">= 2"};
                return kNames[static_cast<std::size_t>(std::ranges::upper_bound(kEdges, miss) - kEdges.begin())];
              };
              ++mRefitReplay[std::format("{} | 4. pass {}: failing step {} (1: the first from the pass's start), seed's miss of its first row [cm] {}", rows, static_cast<int>(out.refitFailedPass), std::min(step.step, 4), missBin(step.startMiss))];
              if (step.barrel) {
                const auto found = mTruth->particles.find(particle);
                const double trueReach = found != mTruth->particles.end() ? 2. / std::abs(propagator.curvature(1.F / found->second.pt, 1.F)) : 0.;
                const auto bin = [](double value) {
                  constexpr std::array<double, 7> kEdges{-20., -5., -1., 0., 1., 5., 20.};
                  constexpr std::array<const char*, 8> kNames{"< -20", "-20..-5", "-5..-1", "-1..0", "0..1", "1..5", "5..20", ">= 20"};
                  return kNames[static_cast<std::size_t>(std::ranges::upper_bound(kEdges, value) - kEdges.begin())];
                };
                ++mRefitReplay[std::format("{} | 3a. pass {}: far side beyond the plane [cm] {}", rows, static_cast<int>(out.refitFailedPass), bin(static_cast<double>(step.farMargin)))];
                ++mRefitReplay[std::format("{} | 3b. pass {}: plane beyond the near side {}", rows, static_cast<int>(out.refitFailedPass), bin(static_cast<double>(step.nearMargin)))];
                ++mRefitReplay[std::format("{} | 3c. true 2R - plane distance {}", rows, bin(trueReach - static_cast<double>(step.planeDistance)))];
              }
            }
          }
          if (mOutStop[index] == kRefitChi2) {
            // The rows the refit took outward: how many not the particle's,
            // and its chi2 against the cut.
            const auto& out = replay->outward.front();
            std::uint32_t wrong = 0U;
            for (std::uint32_t i = 0; i < out.kept; ++i) {
              wrong += mTruth->rows[out.raw[i]].carries(particle) ? 0U : 1U;
            }
            mOutRefitWrong[index] = static_cast<std::uint8_t>(std::min(wrong, 2U) + (std::ranges::any_of(out.steps, &Step::kink) ? 3U : 0U));
            mOutRefitRows[index] = static_cast<std::uint8_t>(out.kept);
            // The particle's largest true kink from the road's first layer to the last outward row's.
            const auto last = out.kept == 0U ? gap : mTruth->rows[out.raw[out.kept - 1U]].layer;
            const auto kink = largestKink(particle, mTruth->rows[out.road.index[0]].layer, last);
            mOutRefitKink[index] = kinkBin(kink);
            if (wrong == 0U) {
              addIncrements(mRefusedIncrements, out);
            }
            mOutRefitRatio[index] = out.refitChi2Max > 0 ? out.refitChi2 / out.refitChi2Max : 0;
          }
          // Who has the own row there at the end: none, a fake, another
          // particle's good track, another track of its own; of this
          // iteration or an earlier one.
          if (ownRow == kNoTrack) {
            mOutHolder[index] = kNoOwnRow;
          } else if (holder[ownRow] == kNoTrack) {
            mOutHolder[index] = kRowFree;
          } else {
            const auto by = holder[ownRow];
            const std::uint8_t kind = particles[by] == study::kNoTrackParticle ? kHeldByFake : (particles[by] == particle ? kHeldByOwn : kHeldByOther);
            mOutHolder[index] = static_cast<std::uint8_t>(kind + (by >= iteration.trackOffset ? std::uint8_t{kThisIteration} : std::uint8_t{0U}));
          }
          mOutStopTaken[index] = takenKind(replay->outward.front(), ownRow, outside, particle);
          mOutChi2[index] = candidateChi2(replay->outward.front(), ownRow, outside);
          // A direct unconfirmed drop: the step took the own row and nothing later
          // (the bid, a witness, the refit) accounts for its loss.
          if (mOutStopTaken[index] == kOwnDropped && (mOutStop[index] == kOtherTaken || mOutStop[index] == kLookaheadOther)) {
            std::tie(mOutNext[index], mOutNextChi2[index]) = afterDrop(replay->outward.front(), ownRow, particle);
            // The particle's true deflection over the step that refused its
            // next row: between its crossings of the dropped row's layer
            // and the next.
            mOutDropKink[index] = kinkBin(largestKink(particle, gap, static_cast<LayerID>(gap + 1U)));
            if (mOutNext[index] == static_cast<std::uint8_t>(StepOutcome::kPasses) && mDebugDrops < 5U) {
              ++mDebugDrops;
              dumpWalk(replay->outward.front(), particle, ownRow);
            }
          }
        }
        continue;
      }
      // The first layer inward of the track's innermost that it misses.
      auto target = static_cast<LayerID>(innermost - 1U);
      while (((missing >> target) & 1U) == 0U) {
        --target;
      }
      std::uint32_t row = kNoTrack;
      for (const auto candidate : own) {
        if (mTruth->rows[candidate].layer == target) {
          row = candidate;
          break;
        }
      }
      mStop[index] = stopReason(*replay, row, static_cast<LayerID>(target + 1U), particle, mStopDetail[index]);
      mStopLayer[index] = target;
      mStopTaken[index] = takenKind(*replay, row, static_cast<LayerID>(target + 1U), particle);
      mInChi2[index] = candidateChi2(*replay, row, static_cast<LayerID>(target + 1U));
      if (mStopTaken[index] == kOwnDropped && (mStop[index] == kOtherTaken || mStop[index] == kLookaheadOther)) {
        std::tie(mInNext[index], mInNextChi2[index]) = afterDrop(*replay, row, particle);
        mInDropKink[index] = target == 0U ? std::uint8_t{3U} : kinkBin(largestKink(particle, static_cast<LayerID>(target - 1U), target));
        // And over the step before, the one that took the dropped row (an
        // upstream scatter the position update absorbed, not the direction).
        mInStepKink[index] = kinkBin(largestKink(particle, target, static_cast<LayerID>(target + 1U)));
      }
      if (mStop[index] == kWindow || mStop[index] == kNotQueried) {
        const auto step = std::ranges::find(replay->steps, static_cast<LayerID>(target + 1U), &Step::layer);
        if (step != replay->steps.end()) {
          ++mWindowAudit[mStop[index] == kWindow ? 0U : 1U][windowAudit(*step, row, target, parameters, frame, geometry, material, propagator, PropagationDirection::kBackward)];
        }
      }
    }
    pathReport(iteration, goodTrack, replayOf);
    seedingFunnel(iteration, parameters, frame, geometry, material, propagator, holder, particles, goodTrack);
    phases.diagnosis = seconds(clock, now());
    clock = now();
    progress(iteration.iteration);
    // The control: each step a complete track of this iteration took onto its
    // particle's own row, from the state the proposal had.
    if (mCsv) {
      for (std::uint32_t index = 0; index < mFindable.size(); ++index) {
        const auto good = goodTrack.find(mFindable[index]);
        if (good == goodTrack.end() || replayOf(good->second) == nullptr) {
          continue;
        }
        const auto& replay = *replayOf(good->second);
        const char* kind = mFate[index] == kComplete ? "control" : nullptr;
        for (std::size_t stepIndex = 0; stepIndex < replay.steps.size(); ++stepIndex) {
          const auto& step = replay.steps[stepIndex];
          const auto own = ownRowOn(mFindable[index], static_cast<LayerID>(step.layer - 1U));
          if (own == kNoTrack) {
            continue;
          }
          const auto seen = std::ranges::find(step.candidates, own, &Candidate::row);
          const bool refused = mFate[index] == kShortFate && step.choice != own;
          if (kind == nullptr && !refused) {
            continue;
          }
          const char* label = kind;
          if (refused) {
            label = seen == step.candidates.end() ? "notQueried" : kOutcomeNames[static_cast<std::size_t>(seen->outcome)];
          }
          record(label, mFindable[index], iteration.iteration, replay, stepIndex, own, seen == step.candidates.end() ? -1.F : seen->chi2, parameters, frame, geometry, material, propagator);
        }
      }
    }
    std::array<double, 5> rest{};
    fakeTransitions(iteration, replayed, recovered);
    if (!mMothers.empty()) {
      clutterStudy(iteration, replayed, parameters, frame, geometry, material, propagator);
      decisionStudy(iteration, replayed, parameters, frame, geometry, material, propagator, goodTrack);
    }
    blockerStudy(iteration, holder, particles);
    rest[0] = seconds(clock, now());
    clock = now();
    endEvidence(iteration, replayed);
    rest[1] = seconds(clock, now());
    clock = now();
    lookaheadPairs(iteration, replayed);
    rest[2] = seconds(clock, now());
    clock = now();
    endPolicies(iteration, replayed, parameters);
    if (mEndCsv) {
      endDecisions(iteration, replayed, parameters, frame, geometry, material, propagator);
    }
    if (mHits) {
      holeCalls(iteration, replayed, parameters, geometry, propagator);
      ownSteps(iteration, replayed, parameters);
    }
    rest[3] = seconds(clock, now());
    clock = now();
    if (parameters.kinkScale2 > 1.F) {
      roadEndAudit(iteration, parameters, frame, geometry, material, propagator);
    }
    rest[4] = seconds(clock, now());
    std::cout << std::format(
      "  anatomy after iteration {} [s]: true kinks {:.1f}, truth walk {:.1f}, replay {:.1f}, diagnosis {:.1f}, fakes {:.1f}, evidence {:.1f}, pairs {:.1f}, "
      "ends {:.1f}, road ends {:.1f}\n",
      iteration.iteration, phases.kinks, phases.walk, phases.replay, phases.diagnosis, rest[0], rest[1], rest[2], rest[3], rest[4]);
    // What the next iteration's follow of the earlier tracks starts from
    // (followEarlier): the tracks and the mask as this one leaves them.
    mSnapshot.resize(tracks);
    for (std::size_t track = 0; track < tracks; ++track) {
      mSnapshot[track] = {.clusters = view.clusters[track], .count = view.count[track], .state = view.state[track], .outwardState = view.outwardState[track]};
    }
    mSnapshotUsed.assign(iteration.used.begin(), iteration.used.end());
  }

  /// How far an iteration took a particle of the class not found before it,
  /// in the funnel's order: eligible -- its floor, and
  /// own rows free on as many of the iteration's road layers as its seeds
  /// hold -- a record of its own (the road fit and the floor passed), that
  /// record selected, and then kept good, kept as a fake (its follow took a
  /// wrong row), or dropped (the length cut after the follow). A pure road
  /// before the fit is the funnel's (`printSeeding`).
  enum SeedStage : std::uint8_t { kSeedFoundBefore,
                                  kSeedNotEligible,
                                  kSeedEligible,
                                  kSeedRecord,
                                  kSeedSelected,
                                  kSeedDropped,
                                  kSeedKeptFake,
                                  kSeedKept,
                                  kSeedStages };
  static constexpr std::array<const char*, kSeedStages> kSeedStageNames{"found before", "not eligible", "eligible", "a record of its own, not selected", "selected", "selected, dropped", "selected, kept as a fake", "selected, kept good"};

  /// `SeedStage` of every class particle in \a iteration, and for the ones
  /// dropped the follow replayed from their record (`replayDropped`).
  template <typename Observation, typename Material>
  void seedingFunnel(const Observation& iteration, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, const std::vector<std::uint32_t>& holder, const std::vector<std::uint64_t>& particles, const std::unordered_map<std::uint64_t, std::uint32_t>& goodTrack)
  {
    const auto it = iteration.iteration;
    if (mSeedStage.size() <= it) {
      mSeedStage.resize(it + 1U);
      mMinClusters.resize(it + 1U, 0U);
      mDropped.resize(it + 1U);
    }
    mMinClusters[it] = parameters.minClusters;
    auto& stage = mSeedStage[it];
    stage.assign(mFindable.size(), kSeedNotEligible);
    for (std::uint32_t track = 0; track < iteration.trackOffset && track < particles.size(); ++track) {
      if (const auto found = mIndex.find(particles[track]); found != mIndex.end()) {
        stage[found->second] = kSeedFoundBefore;
      }
    }
    const double floor = it < mFloors.size() ? mFloors[it] : 0.;
    const std::uint64_t roadLayers = it < mRoadLayers.size() ? mRoadLayers[it] : ~std::uint64_t{0};
    // Eligible on the mask the previous iteration left: with followEarlier
    // the earlier tracks' follow takes rows before this iteration's seeds,
    // which read free here (the observation carries no mask at seeding).
    for (std::uint32_t index = 0; index < mFindable.size(); ++index) {
      if (stage[index] == kSeedFoundBefore) {
        continue;
      }
      const auto found = mTruth->particles.find(mFindable[index]);
      if (found == mTruth->particles.end() || static_cast<double>(found->second.pt) < floor) {
        continue;
      }
      std::uint64_t free = 0U;
      for (const auto row : mRowsOf.at(mFindable[index])) {
        const auto bit = std::uint64_t{1} << mTruth->rows[row].layer;
        free |= (roadLayers & bit) != 0U && (mSnapshotUsed.empty() || mSnapshotUsed[row] == kFree) ? bit : 0U;
      }
      stage[index] = std::cmp_greater_equal(std::popcount(free), parameters.minClusters) ? kSeedEligible : kSeedNotEligible;
    }
    // The records of their own, and the selected ones.
    const auto& selection = *iteration.selection;
    const auto records = selection.records();
    std::vector<std::uint8_t> chosen(static_cast<std::size_t>(selection.recordCount()), 0U);
    for (const auto record : selection.selected()) {
      chosen[static_cast<std::size_t>(record)] = 1U;
    }
    std::vector<std::uint32_t> selectedRecord(mFindable.size(), kNoTrack);
    std::vector<std::uint32_t> bestRecord(mFindable.size(), kNoTrack); // its own, not selected: the best by the quality key
    for (std::size_t record = 0; record < chosen.size(); ++record) {
      const std::span<const std::uint32_t> rows{records.clusters[record].index, records.count[record]};
      const auto found = mIndex.find(study::goodParticle(mTruth->rows, rows));
      if (found == mIndex.end() || stage[found->second] == kSeedFoundBefore) {
        continue;
      }
      auto& one = stage[found->second];
      one = std::max<std::uint8_t>(one, chosen[record] != 0U ? kSeedSelected : kSeedRecord);
      if (chosen[record] != 0U) {
        selectedRecord[found->second] = static_cast<std::uint32_t>(record);
      } else if (auto& best = bestRecord[found->second]; best == kNoTrack || selection::qualityKey(records.count[record], records.kinks[record], records.chi2[record]) < selection::qualityKey(records.count[best], records.kinks[best], records.chi2[best])) {
        best = static_cast<std::uint32_t>(record);
      }
    }
    recordLosses(iteration, stage, bestRecord, chosen);
    beamLineStudy(iteration, stage, parameters, frame, geometry, material, propagator);
    std::vector<std::uint32_t> dropped;
    for (std::uint32_t index = 0; index < mFindable.size(); ++index) {
      const auto good = goodTrack.find(mFindable[index]);
      if (good != goodTrack.end() && good->second >= iteration.trackOffset) {
        stage[index] = kSeedKept;
        continue;
      }
      if (selectedRecord[index] == kNoTrack) {
        continue;
      }
      const auto first = records.clusters[selectedRecord[index]].index[0];
      const auto by = holder[first];
      if (by != kNoTrack && by >= iteration.trackOffset) {
        stage[index] = particles[by] == mFindable[index] ? kSeedKept : kSeedKeptFake;
      } else {
        stage[index] = kSeedDropped;
        dropped.push_back(index);
      }
      if (mKeepFollowLoss && (stage[index] == kSeedDropped || stage[index] == kSeedKeptFake)) {
        mFollowLossRows.push_back(followLoss(iteration, index, stage[index], selectedRecord[index], stage[index] == kSeedKeptFake ? by : kNoTrack, parameters, frame));
      }
    }
    for (std::uint32_t index = 0; index < mFindable.size(); ++index) {
      if (stage[index] != kSeedFoundBefore) {
        mFurthest[index] = std::max(mFurthest[index], stage[index]);
      }
    }
    replayDropped(iteration, dropped, selectedRecord, parameters, frame, geometry, material, propagator);
    // The drops' verdicts, from the replay.
    for (auto i = mFollowLossRows.size(); mKeepFollowLoss && i-- > 0U && mFollowLossRows[i].iteration == iteration.iteration;) {
      if (mFollowLossRows[i].stage == kSeedDropped) {
        const auto why = mLastDropWhy.find(mFollowLossRows[i].index);
        mFollowLossRows[i].why = why != mLastDropWhy.end() ? why->second : std::string{"not replayed"};
      }
    }
  }

  /// Why a class particle's own records lost the selection:
  /// its best own record (the quality key's) against the selected records
  /// sharing a row with it, the best of them by the same key -- another
  /// particle's good record, a fake most of whose rows are this particle's,
  /// or another's fake; and how it ranks: longer, as long with fewer kinks,
  /// as long and kinks with a lower chi2, or ranked below the loser (the
  /// selection is not a plain greedy order) -- or no selected record shares a
  /// row. Kept per particle from the last iteration where it had a record
  /// and none selected; `printRecordLosses` tallies the missed.
  template <typename Observation>
  void recordLosses(const Observation& iteration, const std::vector<std::uint8_t>& stage, const std::vector<std::uint32_t>& bestRecord, const std::vector<std::uint8_t>& chosen)
  {
    const auto records = iteration.selection->records();
    // The selected records by row, once.
    std::unordered_map<std::uint32_t, std::uint32_t> selectedByRow;
    for (std::size_t record = 0; record < chosen.size(); ++record) {
      if (chosen[record] == 0U) {
        continue;
      }
      for (std::uint8_t i = 0; i < records.count[record]; ++i) {
        selectedByRow[records.clusters[record].index[i]] = static_cast<std::uint32_t>(record);
      }
    }
    for (std::uint32_t index = 0; index < mFindable.size(); ++index) {
      if (stage[index] != kSeedRecord || bestRecord[index] == kNoTrack) {
        continue;
      }
      const auto own = bestRecord[index];
      const auto ownKey = selection::qualityKey(records.count[own], records.kinks[own], records.chi2[own]);
      // The selection prunes a record holding a row an earlier iteration's
      // track took before its first round (`selection::Prune`); read on the
      // mask the previous iteration left, which misses the rows a
      // followEarlier follow takes before the seeds.
      if (!mSnapshotUsed.empty() && std::ranges::any_of(std::span<const std::uint32_t>{records.clusters[own].index, records.count[own]}, [&](std::uint32_t row) { return mSnapshotUsed[row] != kFree; })) {
        mRecordLoss[index] = "a row taken by an earlier iteration";
        continue;
      }
      std::uint32_t winner = kNoTrack;
      std::size_t shared = 0U;
      for (std::uint8_t i = 0; i < records.count[own]; ++i) {
        const auto at = selectedByRow.find(records.clusters[own].index[i]);
        if (at == selectedByRow.end()) {
          continue;
        }
        ++shared;
        if (winner == kNoTrack || selection::qualityKey(records.count[at->second], records.kinks[at->second], records.chi2[at->second]) < selection::qualityKey(records.count[winner], records.kinks[winner], records.chi2[winner])) {
          winner = at->second;
        }
      }
      if (winner == kNoTrack) {
        mRecordLoss[index] = "no selected record shares a row";
        continue;
      }
      const std::span<const std::uint32_t> rows{records.clusters[winner].index, records.count[winner]};
      const auto good = study::goodParticle(mTruth->rows, rows);
      const char* kind = good != study::kNoTrackParticle ? "another's good record" : study::majorityParticle(mTruth->rows, rows).particle == mFindable[index] ? "a fake mostly its own"
                                                                                                                                                              : "another's fake";
      const auto winnerKey = selection::qualityKey(records.count[winner], records.kinks[winner], records.chi2[winner]);
      const char* why = winnerKey > ownKey ? "ranked below it" : records.count[winner] > records.count[own] ? "longer"
                                                               : records.kinks[winner] < records.kinks[own] ? "as long, fewer kinks"
                                                               : winnerKey < ownKey                         ? "as long, lower chi2"
                                                                                                            : "tied, the cluster set";
      mRecordLoss[index] = std::format("{} | {} | {} of {} rows shared", kind, why, shared, records.count[own]);
    }
  }

  /// One class particle's own road in a pass with the beam constraint: its own free rows on the pass's
  /// road layers, the first crossing of each, fitted as the road fit does;
  /// for an accepted road above its floor, the beam-line test's inward pass
  /// and update (`clearsBeamLine`) taken apart, and the update again with
  /// the inward state's direction, momentum or both replaced by the MC's at
  /// the road's first row (the fitted covariance kept).
  struct BeamLineRow {
    std::uint32_t index = 0U; ///< into mFindable
    std::uint32_t iteration = 0U;
    std::uint8_t stage = 0U; ///< its SeedStage in the pass
    std::uint8_t count = 0U;
    LayerID first = 0;
    LayerID last = 0;
    std::uint8_t fit = 0U; ///< 0 accepted, 1 the chi2, 2 the propagation, 3 the floor
    float roadChi2 = 0.F;
    std::int8_t inward = -1; ///< the inward pass's PropagationStatus; -1 not run
    std::uint8_t inwardKinks = 0U;
    float chi2 = -1.F; ///< the nominal update's; -1 where it failed
    float residual = 0.F;
    float variance = 0.F;
    float scattering = 0.F;
    float allowanceChi2 = -1.F; ///< with the kink allowance, where the nominal one is above the gate
    float gate = 0.F;
    float dxy = 0.F;             ///< the production vertex's distance from the beam line [cm]
    float pFirst = -1.F;         ///< MC |p| entering the road's first row; -1 without its hit
    float pLast = -1.F;          ///< and its last
    float mcChi2 = -1.F;         ///< the update with the MC direction and momentum
    float mcAnglesChi2 = -1.F;   ///< with the MC direction alone
    float mcMomentumChi2 = -1.F; ///< with the MC momentum alone
    float dPhi = 0.F;            ///< the inward state less the MC [rad]
    float dTheta = 0.F;
    float dQOverPRel = 0.F; ///< (q/p - q/p_MC) / |q/p_MC|
  };

  /// `BeamLineRow` of every class particle eligible in \a iteration (stage
  /// \a stage), if the beam-line rows are kept (`keepBeamLine`). Its own
  /// rows read free on the mask the previous iteration left (with
  /// followEarlier the earlier tracks' follow takes rows before the seeds).
  /// The road fit stage refuses a road whose rows share no time before its
  /// fit; this study has no such test: its question is the fit's and the
  /// beam line's.
  template <typename Observation, typename Material>
  void beamLineStudy(const Observation& iteration, const std::vector<std::uint8_t>& stage, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    if (!mKeepBeamLine || !parameters.beamConstraint) {
      return;
    }
    for (std::uint32_t index = 0; index < mFindable.size(); ++index) {
      if (stage[index] < kSeedEligible) {
        continue;
      }
      const auto particle = mFindable[index];
      std::vector<std::uint32_t> road;
      for (const auto row : roadRows(mRowsOf.at(particle), iteration.iteration)) {
        if (mSnapshotUsed.empty() || mSnapshotUsed[row] == kFree) {
          road.push_back(row);
        }
      }
      if (road.size() < parameters.minClusters) {
        continue;
      }
      const auto& truth = mTruth->particles.at(particle);
      BeamLineRow entry{
        .index = index,
        .iteration = static_cast<std::uint32_t>(iteration.iteration),
        .stage = stage[index],
        .count = static_cast<std::uint8_t>(road.size()),
        .first = mTruth->rows[road.front()].layer,
        .last = mTruth->rows[road.back()].layer,
        .gate = parameters.beamChi2Max,
        .dxy = std::hypot(truth.vx - parameters.beamX, truth.vy - parameters.beamY),
      };
      const auto* firstHit = segmentOfRow(particle, frame, road.front());
      const auto* lastHit = segmentOfRow(particle, frame, road.back());
      entry.pFirst = firstHit != nullptr ? static_cast<float>(std::hypot(firstHit->momentum[0], firstHit->momentum[1], firstHit->momentum[2])) : -1.F;
      entry.pLast = lastHit != nullptr ? static_cast<float>(std::hypot(lastHit->momentum[0], lastHit->momentum[1], lastHit->momentum[2])) : -1.F;
      const auto count = static_cast<std::uint32_t>(road.size());
      const RowChain<Det> chain{.frame = &frame, .rows = road.data()};
      const auto fit = fitRoad(parameters, propagator, geometry, material, chain, static_cast<int>(count));
      entry.roadChi2 = fit.chi2;
      entry.fit = fit.status == RoadFitStatus::kRejected ? 1U : (fit.status == RoadFitStatus::kPropagationFailed ? 2U : 0U);
      if (entry.fit == 0U && !clearsPtFloor(parameters, count, fit)) {
        entry.fit = 3U;
      }
      if (entry.fit != 0U) {
        mBeamLineRows.push_back(entry);
        continue;
      }
      const auto inward = fitInward(propagator, geometry, material, chain, static_cast<int>(count), fit.state, kinkBudget(parameters));
      entry.inward = static_cast<std::int8_t>(inward.status);
      entry.inwardKinks = static_cast<std::uint8_t>(inward.kinks);
      if (inward.status != PropagationStatus::kOk) {
        mBeamLineRows.push_back(entry);
        continue;
      }
      const auto& plane = geometry.planeOf(inward.state);
      const LayerID layer = geometry.surface(inward.state.surface).layer;
      const auto point = pairPointAt(frame.measurements(), road.front());
      const auto updated = [&](BoundTrackState state, Scalar scale) { return updateWithBeamLine(propagator, parameters, plane, layer, point, state, scale); };
      const auto nominal = updated(inward.state, 1.F);
      if (nominal.ok) {
        entry.chi2 = nominal.chi2;
        entry.residual = nominal.residual;
        entry.variance = nominal.variance;
        entry.scattering = nominal.scattering;
        if (!(nominal.chi2 <= parameters.beamChi2Max) && parameters.kinkScale2 > 1.F) {
          const auto allowance = updated(inward.state, parameters.kinkScale2);
          entry.allowanceChi2 = allowance.ok ? allowance.chi2 : -1.F;
        }
      }
      if (firstHit != nullptr) {
        const auto& p = firstHit->momentum;
        const double pMc = std::hypot(p[0], p[1], p[2]);
        const auto phiMc = static_cast<Scalar>(std::atan2(p[1], p[0]));
        const auto thetaMc = static_cast<Scalar>(std::atan2(std::hypot(p[0], p[1]), p[2]));
        const auto qOverPMc = static_cast<Scalar>(std::copysign(1. / pMc, static_cast<double>(inward.state.qOverP())));
        entry.dPhi = static_cast<float>(std::remainder(static_cast<double>(inward.state.phi() - phiMc), 2. * std::numbers::pi));
        entry.dTheta = inward.state.theta() - thetaMc;
        entry.dQOverPRel = static_cast<float>((static_cast<double>(inward.state.qOverP()) - static_cast<double>(qOverPMc)) * pMc);
        const auto replaced = [&](bool angles, bool momentum) {
          BoundTrackState state = inward.state;
          if (angles) {
            state.par[kPhi] = phiMc;
            state.par[kTheta] = thetaMc;
          }
          if (momentum) {
            state.par[kQOverP] = qOverPMc;
          }
          const auto update = updated(state, 1.F);
          return update.ok ? update.chi2 : -1.F;
        };
        entry.mcChi2 = replaced(true, true);
        entry.mcAnglesChi2 = replaced(true, false);
        entry.mcMomentumChi2 = replaced(false, true);
      }
      mBeamLineRows.push_back(entry);
    }
  }

  /// Keep `beamLineStudy`'s rows, for `writeBeamLine`.
  void keepBeamLine() { mKeepBeamLine = true; }

  /// `beamLineStudy`'s rows as a CSV, with each particle's truth and its fate
  /// at the end.
  void writeBeamLine(const std::string& csv) const
  {
    std::ofstream out(csv);
    if (!out) {
      throw std::runtime_error("cannot write " + csv);
    }
    out << "particle,pdg,motherPdg,origin,pt,eta,fate,furthest,iteration,stage,count,first,last,fit,roadChi2,inward,inwardKinks,chi2,residual,variance,"
           "scattering,allowanceChi2,gate,dxy,pFirst,pLast,mcChi2,mcAnglesChi2,mcMomentumChi2,dPhi,dTheta,dQOverPRel\n";
    for (const auto& row : mBeamLineRows) {
      const auto particle = mFindable[row.index];
      const auto& truth = mTruth->particles.at(particle);
      out << std::format(
        "{},{},{},{},{:.4f},{:.3f},{},{},{},{},{},{},{},{},{:.4g},{},{},{:.4g},{:.4g},{:.4g},{:.4g},{:.4g},{:.4g},{:.4g},{:.5g},{:.5g},{:.4g},{:.4g},{:.4g},{:."
        "4g},{:.4g},{:.4g}\n",
        particle, truth.pdg, truth.motherPdg, static_cast<int>(truth.origin), truth.pt, truth.eta, static_cast<int>(mFate[row.index]), static_cast<int>(mFurthest[row.index]), row.iteration, static_cast<int>(row.stage), static_cast<int>(row.count), static_cast<int>(row.first), static_cast<int>(row.last), static_cast<int>(row.fit), row.roadChi2, static_cast<int>(row.inward), static_cast<int>(row.inwardKinks), row.chi2, row.residual, row.variance, row.scattering, row.allowanceChi2, row.gate, row.dxy, row.pFirst, row.pLast, row.mcChi2, row.mcAnglesChi2, row.mcMomentumChi2, row.dPhi, row.dTheta, row.dQOverPRel);
    }
  }

  /// A class particle's selected record that ended dropped or kept as a fake
  /// in an iteration, against the MC: the particle's leg of the record's last row in path
  /// order, its |p| at each own row (the hit's entrance), where it dropped
  /// most; for a fake, the first wrong row after the record in the track's
  /// chain and the momentum lost over the step that took it.
  struct FollowLossRow {
    std::uint32_t index = 0U; ///< into mFindable
    std::uint32_t iteration = 0U;
    std::uint8_t stage = 0U; ///< kSeedDropped or kSeedKeptFake
    std::uint8_t recordCount = 0U;
    std::uint8_t minAfter = 0U;       ///< the iteration's length after the follow
    float dxy = 0.F;                  ///< the production vertex's distance from the beam line [cm]
    std::int8_t legLayers = -1;       ///< the distinct layers of the leg's own rows; -1 without path order
    std::int8_t beyond = -1;          ///< the leg's own rows after the record's last
    float pLegFirst = -1.F;           ///< MC |p| at the leg's first own row
    float pRecordLast = -1.F;         ///< at the record's last row
    float pLegLast = -1.F;            ///< at the leg's last own row
    float rLegLast = -1.F;            ///< that row's radius [cm]
    float maxDrop = 0.F;              ///< the largest drop of |p| between consecutive own rows of the leg
    std::int8_t maxDropLayer = -1;    ///< the layer of the row after it
    std::int8_t maxDropWhere = -1;    ///< 0 up to the record's first row, 1 inside the record, 2 beyond it
    std::int8_t trackCount = -1;      ///< a fake: its rows
    std::int8_t wrong = -1;           ///< and how many not the particle's
    std::int8_t wrongLayer = -1;      ///< the first wrong row after the record's last in the chain; -1 none
    std::int8_t wrongBefore = -1;     ///< 1 if a wrong row precedes the record's first in the chain
    std::int8_t ownOnWrongLayer = -1; ///< 1 if the leg has an own row on that layer after the last own row held before it
    float pBeforeWrong = -1.F;        ///< MC |p| at the last own row the track holds before the wrong row
    float pAtWrong = -1.F;            ///< at the leg's own row on the wrong row's layer
    std::int8_t afterWrong = -1;      ///< the track's rows after that wrong row in the chain (0: it is the last)
    std::int8_t ownAfterWrong = -1;   ///< of them, the particle's own
    std::string why;                  ///< dropped: `replayDropped`'s verdict
  };

  /// `FollowLossRow` of \a index in \a iteration.
  template <typename Observation>
  [[nodiscard]] FollowLossRow followLoss(const Observation& iteration, std::uint32_t index, std::uint8_t stage, std::uint32_t record, std::uint32_t track, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame) const
  {
    const auto particle = mFindable[index];
    const auto records = iteration.selection->records();
    const auto count = records.count[record];
    FollowLossRow row;
    row.index = index;
    row.iteration = static_cast<std::uint32_t>(iteration.iteration);
    row.stage = stage;
    row.recordCount = count;
    row.minAfter = static_cast<std::uint8_t>(parameters.minClustersAfterFollow);
    const auto& truth = mTruth->particles.at(particle);
    row.dxy = std::hypot(truth.vx - parameters.beamX, truth.vy - parameters.beamY);
    const auto& own = mRowsOf.at(particle);
    const auto momentum = [&](std::uint32_t at) {
      const auto* hit = segmentOfRow(particle, frame, at);
      return hit != nullptr ? static_cast<float>(std::hypot(hit->momentum[0], hit->momentum[1], hit->momentum[2])) : -1.F;
    };
    const std::span<const std::uint32_t> recordRows{records.clusters[record].index, count};
    std::size_t first = own.size();
    std::size_t last = own.size();
    for (std::size_t k = 0; k < own.size(); ++k) {
      if (std::ranges::find(recordRows, own[k]) != recordRows.end()) {
        first = first == own.size() ? k : first;
        last = k;
      }
    }
    const bool ordered = !mLegs.empty() && mTimed[index] != 0U && last < own.size();
    std::vector<std::size_t> leg; // own row positions on the record's leg, in path order
    if (ordered) {
      const auto& legs = mLegs[index];
      std::uint64_t layers = 0U;
      for (std::size_t k = 0; k < own.size(); ++k) {
        if (legs[k] == legs[last]) {
          leg.push_back(k);
          layers |= std::uint64_t{1} << mTruth->rows[own[k]].layer;
        }
      }
      row.legLayers = static_cast<std::int8_t>(std::popcount(layers));
      row.beyond = static_cast<std::int8_t>(std::ranges::count_if(leg, [&](std::size_t k) { return k > last; }));
      row.pLegFirst = momentum(own[leg.front()]);
      row.pRecordLast = momentum(own[last]);
      row.pLegLast = momentum(own[leg.back()]);
      const auto point = frame.measurement(own[leg.back()]);
      row.rLegLast = std::hypot(point.x, point.y);
      for (std::size_t i = 1; i < leg.size(); ++i) {
        const auto before = momentum(own[leg[i - 1U]]);
        const auto after = momentum(own[leg[i]]);
        if (before > 0.F && after > 0.F && 1.F - (after / before) > row.maxDrop) {
          row.maxDrop = 1.F - (after / before);
          row.maxDropLayer = static_cast<std::int8_t>(mTruth->rows[own[leg[i]]].layer);
          row.maxDropWhere = static_cast<std::int8_t>(leg[i] <= first ? 0 : (leg[i] <= last ? 1 : 2));
        }
      }
    }
    if (stage != kSeedKeptFake || track == kNoTrack) {
      return row;
    }
    const auto view = iteration.tracks->view();
    const std::span<const std::uint32_t> rows{view.clusters[track].index, view.count[track]};
    row.trackCount = static_cast<std::int8_t>(rows.size());
    row.wrong = static_cast<std::int8_t>(std::ranges::count_if(rows, [&](std::uint32_t r) { return !mTruth->rows[r].carries(particle); }));
    std::size_t recordFirst = rows.size();
    std::size_t recordLast = rows.size();
    for (std::size_t i = 0; i < rows.size(); ++i) {
      if (std::ranges::find(recordRows, rows[i]) != recordRows.end()) {
        recordFirst = recordFirst == rows.size() ? i : recordFirst;
        recordLast = i;
      }
    }
    if (recordLast == rows.size()) {
      return row;
    }
    row.wrongBefore = std::any_of(rows.begin(), rows.begin() + static_cast<std::ptrdiff_t>(recordFirst), [&](std::uint32_t r) { return !mTruth->rows[r].carries(particle); }) ? 1 : 0;
    std::size_t held = recordLast; // the last own row before the first wrong one
    for (std::size_t i = recordLast + 1U; i < rows.size(); ++i) {
      if (mTruth->rows[rows[i]].carries(particle)) {
        held = i;
        continue;
      }
      const auto layer = mTruth->rows[rows[i]].layer;
      row.wrongLayer = static_cast<std::int8_t>(layer);
      row.afterWrong = static_cast<std::int8_t>(rows.size() - i - 1U);
      row.ownAfterWrong = static_cast<std::int8_t>(std::count_if(rows.begin() + static_cast<std::ptrdiff_t>(i) + 1, rows.end(), [&](std::uint32_t r) { return mTruth->rows[r].carries(particle); }));
      row.pBeforeWrong = momentum(rows[held]);
      if (ordered) {
        const auto heldAt = static_cast<std::size_t>(std::ranges::find(own, rows[held]) - own.begin());
        row.ownOnWrongLayer = 0;
        for (const auto k : leg) {
          if (k > heldAt && mTruth->rows[own[k]].layer == layer) {
            row.ownOnWrongLayer = 1;
            row.pAtWrong = momentum(own[k]);
            break;
          }
        }
      }
      break;
    }
    return row;
  }

  /// Keep `followLoss`'s rows, for `writeFollowLoss`.
  void keepFollowLoss() { mKeepFollowLoss = true; }

  /// `followLoss`'s rows as a CSV, with each particle's truth and its fate at
  /// the end.
  void writeFollowLoss(const std::string& csv) const
  {
    std::ofstream out(csv);
    if (!out) {
      throw std::runtime_error("cannot write " + csv);
    }
    out << "particle,pdg,motherPdg,pt,eta,dxy,fate,furthest,iteration,stage,recordCount,minAfter,legLayers,beyond,pLegFirst,pRecordLast,pLegLast,rLegLast,"
           "maxDrop,maxDropLayer,maxDropWhere,trackCount,wrong,wrongLayer,wrongBefore,ownOnWrongLayer,pBeforeWrong,pAtWrong,afterWrong,ownAfterWrong,why\n";
    for (const auto& row : mFollowLossRows) {
      const auto particle = mFindable[row.index];
      const auto& truth = mTruth->particles.at(particle);
      out << std::format("{},{},{},{:.4f},{:.3f},{:.5f},{},{},{},{},{},{},{},{},{:.5g},{:.5g},{:.5g},{:.4g},{:.4g},{},{},{},{},{},{},{},{:.5g},{:.5g},{},{},\"{}\"\n", particle, truth.pdg, truth.motherPdg, truth.pt, truth.eta, row.dxy, static_cast<int>(mFate[row.index]), static_cast<int>(mFurthest[row.index]), row.iteration, static_cast<int>(row.stage), static_cast<int>(row.recordCount), static_cast<int>(row.minAfter), static_cast<int>(row.legLayers), static_cast<int>(row.beyond), row.pLegFirst, row.pRecordLast, row.pLegLast, row.rLegLast, row.maxDrop, static_cast<int>(row.maxDropLayer), static_cast<int>(row.maxDropWhere), static_cast<int>(row.trackCount), static_cast<int>(row.wrong), static_cast<int>(row.wrongLayer), static_cast<int>(row.wrongBefore), static_cast<int>(row.ownOnWrongLayer), row.pBeforeWrong, row.pAtWrong, static_cast<int>(row.afterWrong), static_cast<int>(row.ownAfterWrong), row.why);
    }
  }

  /// The missed class particles by their furthest stage over the
  /// iterations, and for the ones whose furthest was an own record never
  /// selected, why it lost (`recordLosses`), by pT.
  void printRecordLosses() const
  {
    if (mSeedStage.empty() || mFurthest.empty()) {
      return;
    }
    std::map<std::string, std::size_t> furthest; // NOLINT(misc-const-correctness): written in a dependent call the check does not see
    std::map<std::string, std::size_t> losses;   // NOLINT(misc-const-correctness): written in a dependent call the check does not see
    std::size_t missed = 0U;
    for (std::uint32_t index = 0; index < mFindable.size(); ++index) {
      if (mFate[index] != kNotFound) {
        continue;
      }
      ++missed;
      const auto pt = kPathPtNames[pathPtBin(mFindable[index])];
      ++furthest[std::format("{} | {}", pt, kSeedStageNames[mFurthest[index]])];
      if (mFurthest[index] == kSeedRecord) {
        ++losses[std::format("{} | {}", pt, mRecordLoss[index].empty() ? "no reason kept" : mRecordLoss[index])];
      }
    }
    std::string text = std::format("  the class's {} missed by their furthest stage over the iterations, pT [GeV] | stage: particles", missed);
    for (const auto& [key, n] : furthest) {
      text += std::format("\n    {}: {}", key, n);
    }
    text += "\n  of them, the furthest an own record never selected: what beat it (the best selected record sharing a row) | how it ranks | rows shared";
    for (const auto& [key, n] : losses) {
      text += std::format("\n    {}: {}", key, n);
    }
    std::cout << text << '\n';
  }

  /// The follow of each dropped record again: from the
  /// record's rows refitted, the outward walk and the inward walk on the
  /// mask the follow read, and the extended refit of what they took -- no
  /// bidding -- and why the track stayed below the length after the follow:
  /// the walks short (and the particle's next row beyond them on its leg,
  /// as the step from its last row decided it, or none: only the way back
  /// beyond), or the refit dropping what they took. With each, the
  /// particle's second rows on a layer and its own rows free beyond the
  /// record on its leg: what the follow could have added.
  template <typename Observation, typename Material>
  void replayDropped(const Observation& iteration, const std::vector<std::uint32_t>& dropped, const std::vector<std::uint32_t>& selectedRecord, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    auto& lines = mDropped[iteration.iteration];
    lines.clear();
    mDroppedWalks.clear();
    mLastDropWhy.clear();
    if (dropped.empty() || !parameters.follow || iteration.followUsed.empty()) {
      return;
    }
    if (!mMothers.empty()) {
      mDroppedWalks.resize(dropped.size());
    }
    const auto records = iteration.selection->records();
    std::vector<std::string> keys(dropped.size());
    // The stage as the iteration instantiated it (see replayExtension).
    withFollowStage(iteration, parameters, frame, geometry, material, propagator, iteration.followUsed.data(), iteration.trackOffset, [&](const auto& stage) {
      o2::gtrk::detail::parallelFor(*mArena, 0U, dropped.size(), 64, [&](std::size_t a, std::size_t b) {
        for (std::size_t i = a; i < b; ++i) {
          const auto index = dropped[i];
          const auto particle = mFindable[index];
          const auto record = selectedRecord[index];
          ClusterList<kMaxClusters> road{};
          const auto count = records.count[record];
          for (std::size_t k = 0; k < kMaxClusters; ++k) {
            road.index[k] = k < count ? records.clusters[record].index[k] : kNoRow;
          }
          // The particle's rows: second rows, and own rows free beyond the
          // record's outermost on its leg.
          const auto& own = mRowsOf.at(particle);
          std::string seconds = "no second row";
          std::size_t beyondFree = 0U;
          std::size_t last = own.size();
          if (!mLegs.empty() && mTimed[index] != 0U) {
            const auto& legs = mLegs[index];
            bool overlap = false;
            bool back = false;
            std::array<std::uint8_t, kLayers> legOnLayer{};
            legOnLayer.fill(255U);
            for (std::size_t k = 0; k < own.size(); ++k) {
              const auto layer = mTruth->rows[own[k]].layer;
              if (legOnLayer[layer] != 255U) {
                overlap = overlap || legOnLayer[layer] == legs[k];
                back = back || legOnLayer[layer] != legs[k];
              }
              legOnLayer[layer] = legs[k];
              if (std::find(road.index, road.index + count, own[k]) != road.index + count) {
                last = k;
              }
            }
            seconds = overlap ? (back ? "an overlap and the way back" : "an overlap") : (back ? "the way back" : "no second row");
            for (std::size_t k = last + 1U; last < own.size() && k < own.size(); ++k) {
              beyondFree += legs[k] == legs[last] && iteration.followUsed[own[k]] == kFree ? 1U : 0U;
            }
          }
          // The starts as the iteration's refit made them (TrackRefit.h):
          // the inward pass from the record's outward state, which is the
          // track's outward state.
          const auto& outwardState = records.outwardState[record];
          const auto fit = fitInward(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = road.index}, static_cast<int>(count), outwardState, kinkBudget(parameters));
          std::string why;
          std::string extra = "first wrong row: n/a | cut before it: n/a | own rows taken from later legs: n/a";
          if (fit.status != PropagationStatus::kOk) {
            why = "the record's refit failed";
          } else {
            Walk outward;
            Walk inward;
            const auto walkedProposal = recordFollow(stage, fit.state, outwardState, road, count, inward, outward);
            const auto* outRows = walkedProposal.outward;
            const auto* inRows = walkedProposal.inward;
            const auto walked = walkedProposal.outwardCount;
            const auto inWalked = walkedProposal.inwardCount;
            if (!mDroppedWalks.empty()) {
              mDroppedWalks[i] = {.outward = outward, .inward = inward, .road = road, .roadCount = static_cast<std::uint8_t>(count), .walked = true};
            }
            // Each direction cut back to its kink steps' witnesses, as
            // before the refit (`KeepWitnesses`); no bidding.
            const bool kinkSteps = parameters.kinkScale2 > 1.F && parameters.maxKinks != 0U;
            const auto keptOut = kinkSteps ? follow::keptWithWitnesses(walked, outward.witnessed) : walked;
            const auto keptIn = kinkSteps ? follow::keptWithWitnesses(inWalked, inward.witnessed) : inWalked;
            std::vector<std::uint32_t> combined;
            for (auto k = static_cast<std::size_t>(keptIn); k-- > 0U;) {
              combined.push_back(inRows[k]);
            }
            combined.insert(combined.end(), road.index, road.index + count);
            combined.insert(combined.end(), outRows, outRows + keptOut);
            const auto need = parameters.minClustersAfterFollow > count ? parameters.minClustersAfterFollow - count : 0U;
            const auto added = static_cast<std::uint32_t>(keptOut + keptIn);
            std::size_t wrong = 0U;
            for (std::size_t k = 0; k < added; ++k) {
              wrong += mTruth->rows[k < keptOut ? outRows[k] : inRows[k - keptOut]].carries(particle) ? 0U : 1U;
            }
            const std::string took = std::format("took +{} ({} not its own) of the +{} needed", added, wrong, need);
            // The outward walk's first wrong row: its layer and the particle's
            // own row at that step; the walk cut before it, refitted (a
            // truth-aided bound, not a policy); its own rows from later legs.
            std::string firstWrong = "first wrong row: none outward";
            std::string cut = "cut before it: n/a";
            std::size_t wayBack = 0U;
            if (!mLegs.empty() && mTimed[index] != 0U && last < own.size()) {
              for (std::size_t k = 0; k < keptOut; ++k) {
                const auto here = std::ranges::find(own, outRows[k]);
                wayBack += here != own.end() && mLegs[index][static_cast<std::size_t>(here - own.begin())] > mLegs[index][last] ? 1U : 0U;
              }
            }
            for (std::size_t w = 0; w < keptOut; ++w) {
              if (mTruth->rows[outRows[w]].carries(particle)) {
                continue;
              }
              const auto layer = mTruth->rows[outRows[w]].layer;
              const auto step = std::ranges::find(outward.steps, outRows[w], &Step::choice);
              std::string own0 = "no own row on its layer";
              bool onLayer = false;
              for (const auto row : own) {
                if (mTruth->rows[row].layer != layer || step == outward.steps.end()) {
                  continue;
                }
                onLayer = true;
                const auto& seen = step->kink ? step->scaled : step->candidates;
                const auto at = std::ranges::find(seen, row, &Candidate::row);
                const auto taken = std::ranges::find(seen, outRows[w], &Candidate::row);
                if (at == seen.end()) {
                  own0 = "own row not visited";
                } else if (at->outcome == StepOutcome::kPasses) {
                  own0 = taken != seen.end() && taken->chi2 > at->chi2 ? "own row passing, taken higher" : "own row passing, taken lower";
                  break;
                } else {
                  own0 = std::format("own row {}", kOutcomeNames[static_cast<std::size_t>(at->outcome)]);
                }
              }
              if (!onLayer && !mLegs.empty() && mTimed[index] != 0U) {
                // Own rows of the leg beyond the last own row the walk held.
                std::size_t at = last;
                for (std::size_t k = 0; k < w; ++k) {
                  const auto here = std::ranges::find(own, outRows[k]);
                  at = here != own.end() ? std::max(at, static_cast<std::size_t>(here - own.begin())) : at;
                }
                bool sameLeg = false;
                bool later = false;
                for (std::size_t k = at + 1U; at < own.size() && k < own.size(); ++k) {
                  sameLeg = sameLeg || mLegs[index][k] == mLegs[index][at];
                  later = later || mLegs[index][k] != mLegs[index][at];
                }
                own0 = sameLeg ? "no own row on its layer, own rows of its leg beyond" : (later ? "past its last on its leg, the way back beyond" : "past its last row");
              }
              firstWrong = std::format("first wrong row: layer {}, {}{}", static_cast<int>(layer), own0, step != outward.steps.end() && step->kink ? " (kink step)" : "");
              std::vector<std::uint32_t> trimmed(road.index, road.index + count);
              // The prefix's own witnesses kept too: a kink step whose
              // witness the cut removes goes with it.
              const auto prefix = kinkSteps ? follow::keptWithWitnesses(static_cast<std::uint8_t>(w), outward.witnessed) : static_cast<std::uint8_t>(w);
              trimmed.insert(trimmed.end(), outRows, outRows + prefix);
              if (trimmed.size() < parameters.minClustersAfterFollow) {
                cut = std::format("cut before it: {} rows, short of the length", trimmed.size());
              } else {
                const auto refit = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = trimmed.data()}, static_cast<int>(trimmed.size()), kinkBudget(parameters));
                // Kept by the refit and the length: and its rows still free
                // after the iteration's bids and marks?
                const bool held = std::any_of(outRows, outRows + prefix, [&](std::uint32_t row) { return iteration.used[row] != kFree; });
                cut = refit.status != PropagationStatus::kOk ? std::format("cut before it: the refit fails in the {} pass", refit.failedPass == FitPass::kInward ? "inward" : "outward") : (follow::passesRoadChi2(parameters, refit) ? std::string{held ? "cut before it: kept, a row held after the iteration" : "cut before it: kept"} : std::string{"cut before it: chi2 over the cut"});
              }
              break;
            }
            extra = std::format("{} | {} | own rows taken from later legs: {}", firstWrong, cut, wayBack);
            if (added < need) {
              // The particle's next row beyond the outward walk on its leg.
              std::size_t at = last;
              for (std::size_t k = 0; k < walked; ++k) {
                const auto here = std::ranges::find(own, outRows[k]);
                at = here != own.end() ? std::max(at == own.size() ? 0U : at, static_cast<std::size_t>(here - own.begin())) : at;
              }
              std::uint32_t next = kNoTrack;
              bool laterLeg = false;
              if (at < own.size() && !mLegs.empty() && mTimed[index] != 0U) {
                const auto& legs = mLegs[index];
                for (std::size_t k = at + 1U; k < own.size(); ++k) {
                  if (legs[k] == legs[at]) {
                    next = own[k];
                    break;
                  }
                  laterLeg = true;
                }
              }
              if (next == kNoTrack) {
                why = std::format("the walks short, {}: no own row beyond on its leg, {}", took, laterLeg ? "the way back beyond" : "none beyond");
              } else {
                outward.bidKept = walked;
                outward.kept = walked;
                outward.refit = kRefitKept;
                std::uint8_t detail = 0U;
                const auto reason = stopReason(outward, next, mTruth->rows[own[at]].layer, particle, detail);
                const bool droppedUnconfirmed = outward.dropped == next;
                why = std::format("the walks short, {}: its next row (layer {}) {}", took, static_cast<int>(mTruth->rows[next].layer), droppedUnconfirmed ? "taken, dropped after (the unconfirmed drop)" : kShortNames[reason]);
              }
            } else {
              const auto refit = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = combined.data()}, static_cast<int>(combined.size()), kinkBudget(parameters));
              if (refit.status != PropagationStatus::kOk) {
                why = std::format("the walks long enough, {}: the refit failed, status {} in the {} pass", took, static_cast<int>(refit.status), refit.failedPass == FitPass::kInward ? "inward" : "outward");
              } else if (!(follow::passesRoadChi2(parameters, refit))) {
                why = std::format("the walks long enough, {}: the refit's chi2 over the cut", took);
              } else {
                why = std::format("the walks long enough, {}: the refit keeps it (the bid, or the iteration's follow differs)", took);
              }
            }
          }
          keys[i] = std::format("{} | {} | {} | second rows: {} | own rows free beyond on its leg: {}", kPathPtNames[pathPtBin(particle)], why, extra, seconds, std::min<std::size_t>(beyondFree, 3U) == 3U ? std::string{"3+"} : std::to_string(beyondFree));
        }
      });
    });
    for (std::size_t i = 0; i < keys.size(); ++i) {
      ++lines[keys[i]];
      if (mKeepFollowLoss) {
        mLastDropWhy[dropped[i]] = keys[i];
      }
    }
  }

  /// The seeding funnel of every iteration, the class's
  /// particles not found before it, by pT; with \a longestRoad (the stage
  /// funnel's longest road of a particle's own in an iteration and why its
  /// road made no record, a `study::FitLoss`, if the funnel ran) the
  /// eligible without a record split by a pure road and that reason; and
  /// the dropped records' replays. Particle counts: one selected record a
  /// particle.
  void printSeeding(const std::function<std::optional<std::pair<std::uint32_t, int>>(std::uint64_t, std::size_t)>& longestRoad) const
  {
    for (std::size_t it = 0; it < mSeedStage.size(); ++it) {
      const auto& stage = mSeedStage[it];
      if (stage.empty()) {
        continue;
      }
      std::array<std::array<std::size_t, kSeedStages + 2U>, kPathPtNames.size()> counts{};
      std::array<std::map<std::string, std::size_t>, kPathPtNames.size()> pureWhy{};
      for (std::uint32_t index = 0; index < stage.size(); ++index) {
        const auto bin = pathPtBin(mFindable[index]);
        ++counts[bin][stage[index]];
        if (stage[index] == kSeedEligible) {
          const auto road = longestRoad ? longestRoad(mFindable[index], it) : std::nullopt;
          const bool pure = road && road->first >= mMinClusters[it];
          counts[bin][kSeedStages + (pure ? 1U : 0U)] += road ? 1U : 0U;
          if (pure) {
            ++pureWhy[bin][study::kFitLossNames[static_cast<std::size_t>(road->second)]];
          }
        }
      }
      // "selected" is only a step on the way: every selected record ends
      // dropped, kept as a fake or kept good, so its column is not printed.
      std::string text = std::format("  seeding funnel of iteration {}, by pT [GeV]:", it);
      for (std::size_t s = 0; s < kSeedStages; ++s) {
        text += s == kSeedSelected ? std::string{} : std::format(" {}{}", kSeedStageNames[s], s + 1U < kSeedStages ? " /" : "");
      }
      text += " | eligible without a record: no pure road / a pure road, no record (the stage funnel's; by its reason)";
      for (std::size_t b = 0; b < counts.size(); ++b) {
        text += std::format("\n    {}:", kPathPtNames[b]);
        for (std::size_t s = 0; s < kSeedStages; ++s) {
          text += s == kSeedSelected ? std::string{} : std::format(" {}", counts[b][s]);
        }
        text += std::format(" | {} / {} (", counts[b][kSeedStages], counts[b][kSeedStages + 1U]);
        for (const auto& [why, n] : pureWhy[b]) {
          text += std::format("{} {}; ", why, n);
        }
        text += ")";
      }
      if (it < mDropped.size() && !mDropped[it].empty()) {
        text += std::format("\n    the dropped of iteration {}, the follow replayed from the record (pT | why | the particle's rows):", it);
        for (const auto& [key, count] : mDropped[it]) {
          text += std::format("\n      {}: {}", key, count);
        }
      }
      std::cout << text << '\n';
    }
    // The not found at the end, by their last iteration's stage.
    if (!mSeedStage.empty()) {
      std::array<std::array<std::size_t, kSeedStages>, kMissedReasons> cross{};
      const auto& last = mSeedStage.back();
      for (std::uint32_t index = 0; index < last.size(); ++index) {
        if (mFate[index] == kNotFound) {
          ++cross[mReason[index]][last[index]];
        }
      }
      std::string text = "  the not found by the anatomy's reason, against the last iteration's seeding stage:";
      for (std::size_t reason = 0; reason < kMissedReasons; ++reason) {
        std::size_t total = 0U;
        for (const auto n : cross[reason]) {
          total += n;
        }
        if (total == 0U) {
          continue;
        }
        text += std::format("\n    {} ({}):", kMissedNames[reason], total);
        for (std::size_t s = 0; s < kSeedStages; ++s) {
          text += cross[reason][s] != 0U ? std::format(" {} {};", kSeedStageNames[s], cross[reason][s]) : std::string{};
        }
      }
      std::cout << text << '\n';
    }
  }

  /// The pT bins of the path tables [GeV]: the 0.1 and 0.15 rungs apart.
  static constexpr std::array<double, 4> kPathPtEdges{0.15, 0.2, 0.3, 0.5};
  static constexpr std::array<const char*, 5> kPathPtNames{"0.1-0.15", "0.15-0.2", "0.2-0.3", "0.3-0.5", ">= 0.5"};
  [[nodiscard]] std::size_t pathPtBin(std::uint64_t particle) const
  {
    const double pt = truePt(mTruth->particles.find(particle));
    return static_cast<std::size_t>(std::ranges::upper_bound(kPathPtEdges, pt) - kPathPtEdges.begin());
  }

  /// How a follow walk ended (`FollowEnd`), as the printed tables name it.
  static constexpr std::array<const char*, 4> kFollowEndNames{"no target", "holes", "no row", "capacity"};
  /// A replayed step's outcome on a row, as the printed tables name it: the
  /// six `StepOutcome`s, then `kNotVisited` and `kNoOwnNext`.
  static constexpr std::array<const char*, 8> kStepOutcomeNames{"held", "outside the window", "crossing a layer", "no update", "chi2", "passing", "not visited", "no own next row"};
  /// The names of the bins of `kWalkPtEdges` [GeV].
  static constexpr std::array<const char*, 5> kWalkPtNames{"< 0.2", "0.2-0.3", "0.3-0.6", "0.6-1", ">= 1"};

  /// Where a track's outermost own row (in path order) leaves the particle:
  /// at its last row; the next an overlap (its layer again, its leg); rows of
  /// its leg beyond (the walk could go on); only rows of later legs beyond
  /// (it ends at the turning point, the way back left).
  enum PathEnd : std::uint8_t { kEndLast,
                                kEndOverlap,
                                kEndLegBeyond,
                                kEndTurning,
                                kPathEnds };
  static constexpr std::array<const char*, kPathEnds> kPathEndNames{"at its last row", "an overlap next", "rows of its leg beyond", "at the turning point, only later legs beyond"};
  /// \a index's (a findable's) path end after the own rows \a held, a bit
  /// per path position; kPathEnds without a timed path or a held row.
  [[nodiscard]] PathEnd pathEnd(std::uint32_t index, const std::vector<bool>& held) const
  {
    if (mLegs.empty() || mTimed[index] == 0U) {
      return kPathEnds;
    }
    const auto& legs = mLegs[index];
    const auto& rows = mRowsOf.at(mFindable[index]);
    std::size_t last = rows.size();
    for (std::size_t k = 0; k < rows.size(); ++k) {
      last = held[k] ? k : last;
    }
    if (last == rows.size()) {
      return kPathEnds;
    }
    if (last + 1U == rows.size()) {
      return kEndLast;
    }
    if (legs[last + 1U] == legs[last] && mTruth->rows[rows[last + 1U]].layer == mTruth->rows[rows[last]].layer) {
      return kEndOverlap;
    }
    for (std::size_t k = last + 1U; k < rows.size(); ++k) {
      if (legs[k] == legs[last]) {
        return kEndLegBeyond;
      }
    }
    return kEndTurning;
  }

  /// The class in path order, after each iteration (the
  /// last one's printed): the particles' rows by leg, their second rows on a
  /// layer split into overlaps (the same leg) and returns (a later one); the
  /// rows a track could reach -- the distinct layers of the first outward
  /// leg, the tracker's one row per layer, an upper bound (no holes, links), and
  /// with the way back's distinct layers too; where each good track ends
  /// against its particle's rows (`PathEnd`); for the iteration's own
  /// tracks, by how their outward walk ended.
  template <typename Observation, typename ReplayOf>
  void pathReport(const Observation& iteration, const std::unordered_map<std::uint64_t, std::uint32_t>& goodTrack, const ReplayOf& replayOf)
  {
    if (mLegs.empty()) {
      return;
    }
    const auto view = iteration.tracks->view();
    struct Bin {
      std::size_t particles = 0U, timed = 0U, withReturn = 0U, withOverlap = 0U, overlapOnly = 0U, returnOnly = 0U, both = 0U;
      std::size_t rows0 = 0U, rows1 = 0U, rows2 = 0U, ceiling0 = 0U, ceilingBack = 0U, reach7 = 0U, reach7Back = 0U;
      std::size_t found = 0U, held = 0U, heldReturn = 0U, heldLater = 0U, ceilingFound = 0U;
      std::array<std::size_t, kPathEnds + 1U> ends{};
    };
    std::array<Bin, kPathPtNames.size()> bins{};
    std::array<std::array<std::size_t, kPathEnds + 1U>, 4> walkEnds{}; // [FollowEnd][PathEnd], the iteration's own good tracks
    std::vector<bool> held;
    for (std::uint32_t index = 0; index < mFindable.size(); ++index) {
      const auto particle = mFindable[index];
      auto& bin = bins[pathPtBin(particle)];
      ++bin.particles;
      if (mTimed[index] == 0U) {
        continue;
      }
      ++bin.timed;
      const auto& rows = mRowsOf.at(particle);
      const auto& legs = mLegs[index];
      std::uint64_t layers0 = 0U;
      std::uint64_t layers1 = 0U;
      bool overlap = false;
      bool back = false;
      std::array<std::uint8_t, kLayers> legOnLayer{};
      legOnLayer.fill(255U);
      for (std::size_t k = 0; k < rows.size(); ++k) {
        const auto layer = mTruth->rows[rows[k]].layer;
        const auto bit = std::uint64_t{1} << layer;
        bin.rows0 += legs[k] == 0U ? 1U : 0U;
        bin.rows1 += legs[k] == 1U ? 1U : 0U;
        bin.rows2 += legs[k] >= 2U ? 1U : 0U;
        layers0 |= legs[k] == 0U ? bit : 0U;
        layers1 |= legs[k] == 1U ? bit : 0U;
        if (legOnLayer[layer] != 255U) {
          overlap = overlap || legOnLayer[layer] == legs[k];
          back = back || legOnLayer[layer] != legs[k];
        }
        legOnLayer[layer] = legs[k];
      }
      const bool returns = std::ranges::any_of(legs, [](std::uint8_t leg) { return leg != 0U; });
      bin.withReturn += returns ? 1U : 0U;
      bin.withOverlap += overlap ? 1U : 0U;
      bin.overlapOnly += overlap && !back ? 1U : 0U;
      bin.returnOnly += back && !overlap ? 1U : 0U;
      bin.both += overlap && back ? 1U : 0U;
      const auto ceiling0 = static_cast<std::size_t>(std::popcount(layers0));
      const auto ceilingBack = ceiling0 + static_cast<std::size_t>(std::popcount(layers1));
      bin.ceiling0 += ceiling0;
      bin.ceilingBack += ceilingBack;
      bin.reach7 += ceiling0 >= 7U ? 1U : 0U;
      bin.reach7Back += ceilingBack >= 7U ? 1U : 0U;
      const auto good = goodTrack.find(particle);
      if (good == goodTrack.end()) {
        continue;
      }
      const auto track = good->second;
      held.assign(rows.size(), false);
      std::size_t count = 0U;
      for (std::uint8_t i = 0; i < view.count[track]; ++i) {
        const auto at = std::ranges::find(rows, view.clusters[track].index[i]);
        if (at != rows.end()) {
          held[static_cast<std::size_t>(at - rows.begin())] = true;
          ++count;
          const auto leg = legs[static_cast<std::size_t>(at - rows.begin())];
          bin.heldReturn += leg == 1U ? 1U : 0U;
          bin.heldLater += leg >= 2U ? 1U : 0U;
        }
      }
      ++bin.found;
      bin.held += count;
      bin.ceilingFound += ceiling0;
      const auto end = pathEnd(index, held);
      ++bin.ends[end];
      if (track >= iteration.trackOffset) {
        if (const auto* replay = replayOf(track); replay != nullptr && !replay->outward.empty() && replay->outward.front().agrees) {
          ++walkEnds[static_cast<std::size_t>(replay->outward.front().end)][end];
        }
      }
    }
    const auto mean = [](std::size_t sum, std::size_t n) { return n == 0U ? 0. : static_cast<double>(sum) / static_cast<double>(n); };
    std::string text = std::format(
      "  path order, after iteration {}: rows by the MC hits' time; reordered against the distance from the vertex: class {} of {} (first "
      "row on a "
      "layer changed {}), others {} of {}; kept in the distance's order (a row without its segment): class {}, others {}",
      iteration.iteration, mPathCounts[1][1], mPathCounts[1][0] + mPathCounts[1][1] + mPathCounts[1][2], mPathCounts[1][3], mPathCounts[0][1], mPathCounts[0][0] + mPathCounts[0][1] + mPathCounts[0][2], mPathCounts[1][2], mPathCounts[0][2]);
    text +=
      "\n    by pT [GeV]: particles (timed) | a second row on a layer: overlap only (its leg) / a later leg only / both | mean rows out (leg 0) / the first "
      "way "
      "back (leg 1) / later legs | distinct layers out (today's ceiling, an upper bound) / with the first way back | % reaching 7 out / with the first way "
      "back";
    for (std::size_t b = 0; b < bins.size(); ++b) {
      const auto& one = bins[b];
      text += std::format("\n      {}: {} ({}) | {} / {} / {} | {:.2f} / {:.2f} / {:.2f} | {:.2f} / {:.2f} | {:.1f} / {:.1f}", kPathPtNames[b], one.particles, one.timed, one.overlapOnly, one.returnOnly, one.both, mean(one.rows0, one.timed), mean(one.rows1, one.timed), mean(one.rows2, one.timed), mean(one.ceiling0, one.timed), mean(one.ceilingBack, one.timed), study::percent(one.reach7, one.timed), study::percent(one.reach7Back, one.timed));
    }
    text += "\n    the good tracks by pT: found | mean own rows held / the ceiling out | rows held from the first way back / later legs | where the track ends:";
    for (std::size_t end = 0; end < kPathEnds; ++end) {
      text += std::format(" {} /", kPathEndNames[end]);
    }
    for (std::size_t b = 0; b < bins.size(); ++b) {
      const auto& one = bins[b];
      text += std::format("\n      {}: {} | {:.2f} / {:.2f} | {} / {} |", kPathPtNames[b], one.found, mean(one.held, one.found), mean(one.ceilingFound, one.found), one.heldReturn, one.heldLater);
      for (std::size_t end = 0; end < kPathEnds; ++end) {
        text += std::format(" {}", one.ends[end]);
      }
    }
    text += std::format("\n    the iteration's own good tracks, by how the outward walk ended (replayed, agreeing) and where the track ends:");
    for (std::size_t end = 0; end < kFollowEndNames.size(); ++end) {
      text += std::format("\n      {}:", kFollowEndNames[end]);
      for (std::size_t at = 0; at < kPathEnds; ++at) {
        text += std::format(" {} {};", kPathEndNames[at], walkEnds[end][at]);
      }
    }
    mPathLines = text;
  }

  /// \a particle's largest true kink, either plane, over theta0, between
  /// crossings of consecutive layers from \a first to \a last; -1 without one.
  /// The kinks are grouped by particle (`mKinkBegin`): the truth walk asks
  /// for them at every step, too often for a scan of them all.
  [[nodiscard]] float largestKink(std::uint64_t particle, LayerID first, LayerID last) const
  {
    const auto at = mIndex.find(particle);
    if (at == mIndex.end() || at->second + 1U >= mKinkBegin.size()) {
      return -1.F;
    }
    float largest = -1.F;
    for (auto i = mKinkBegin[at->second]; i < mKinkBegin[at->second + 1U]; ++i) {
      const auto& kink = mKinks[i];
      if (kink.layer >= first && kink.to <= last) {
        largest = std::max({largest, std::abs(kink.theta), std::abs(kink.phi)});
      }
    }
    return largest;
  }

  /// \a particle's true kink over the crossing ending at \a layer (\a into)
  /// or starting there, either plane over theta0; -1 without one (the
  /// crossing before C1 is ML 4 -> C1, not 10 -> 11).
  [[nodiscard]] float kinkAt(std::uint64_t particle, LayerID layer, bool into) const
  {
    const auto at = mIndex.find(particle);
    if (at == mIndex.end() || at->second + 1U >= mKinkBegin.size()) {
      return -1.F;
    }
    for (auto i = mKinkBegin[at->second]; i < mKinkBegin[at->second + 1U]; ++i) {
      const auto& kink = mKinks[i];
      if ((into ? kink.to : kink.layer) == layer) {
        return std::max(std::abs(kink.theta), std::abs(kink.phi));
      }
    }
    return -1.F;
  }

  void print() const
  {
    std::array<std::size_t, 3> fates{};
    for (const auto fate : mFate) {
      ++fates[fate];
    }
    const auto n = mFindable.size();
    const auto pc = [&](std::size_t count) { return study::percent(count, n); };
    std::cout << std::format(
      "complete anatomy: {} physical primaries on every layer of mask {}: complete {} ({:.2f} %), short {} ({:.2f} %), not found {} "
      "({:.2f} %); the replayed follow agrees with the iteration's on {} of {} tracks\n",
      n, mMask, fates[kComplete], pc(fates[kComplete]), fates[kShortFate], pc(fates[kShortFate]), fates[kNotFound], pc(fates[kNotFound]), mAgree, mReplayed);
    printFidelity();
    {
      std::string text = "  good tracks of the findable by the layers of the mask they hold:";
      for (std::size_t layers = 0; layers < mLayersFound.size(); ++layers) {
        if (mLayersFound[layers] != 0U) {
          text += std::format(" {}: {} ({:.2f} %);", layers, mLayersFound[layers], pc(mLayersFound[layers]));
        }
      }
      std::cout << text << '\n';
    }
    std::array<std::size_t, kShortReasons + 1U> stops{};
    std::array<std::size_t, kMissedReasons> missed{};
    for (std::size_t index = 0; index < n; ++index) {
      if (mFate[index] == kShortFate) {
        ++stops[mStop[index]];
      } else if (mFate[index] == kNotFound) {
        ++missed[mReason[index]];
      }
    }
    std::string line = "  short, by the decision that stopped it:";
    for (std::size_t reason = 0; reason < kShortReasons; ++reason) {
      if (stops[reason] != 0U) {
        line += std::format(" {} {} ({:.2f} %);", kShortNames[reason], stops[reason], pc(stops[reason]));
      }
    }
    std::cout << line << '\n';
    // The blocks of the follow's replay, on a run whose iterations followed
    // (none on the ITS: they would be empty or zero).
    const bool followed = mReplayed + mOutReplayed != 0U;
    if (followed) {
      std::array<std::size_t, 8> kinkRefused{};
      std::array<std::array<std::size_t, 8>, 2> unconfirmed{};
      for (std::size_t index = 0; index < n; ++index) {
        if (mFate[index] == kShortFate && mStop[index] == kKinkRefused) {
          ++kinkRefused[mStopDetail[index] % 8U];
        } else if (mFate[index] == kShortFate && mStop[index] == kKinkUnconfirmed) {
          ++unconfirmed[mStopDetail[index] / 8U][mStopDetail[index] % 8U];
        }
      }
      line = "  kink step refused it, under the scaled gate:";
      for (std::size_t outcome = 0; outcome < kStepOutcomeNames.size(); ++outcome) {
        if (kinkRefused[outcome] != 0U) {
          line += std::format(" {} {};", kStepOutcomeNames[outcome], kinkRefused[outcome]);
        }
      }
      std::cout << line << '\n';
      for (std::size_t beam = 0; beam < 2U; ++beam) {
        line = beam == 0U ? "  kink step unconfirmed, the particle's next row under the nominal gate:" : "  kink step unconfirmed, the beam line:";
        for (std::size_t outcome = 0; outcome < kStepOutcomeNames.size(); ++outcome) {
          if (unconfirmed[beam][outcome] != 0U) {
            line += std::format(" {} {};", kStepOutcomeNames[outcome], unconfirmed[beam][outcome]);
          }
        }
        std::cout << line << '\n';
      }
    }
    line = "  short, by the layer it missed:";
    std::array<std::size_t, kLayers> stopLayers{};
    for (std::size_t index = 0; index < n; ++index) {
      if (mFate[index] == kShortFate && mStopLayer[index] < kLayers) {
        ++stopLayers[mStopLayer[index]];
      }
    }
    for (std::size_t layer = 0; layer < kLayers; ++layer) {
      if (stopLayers[layer] != 0U) {
        line += std::format(" {}:{}", layer, stopLayers[layer]);
      }
    }
    std::cout << line << '\n';
    // What the passing step took instead of the own row,
    // inward; then the short of an outer layer, by the outward step.
    const auto takenLine = [&](const char* title, const std::vector<std::uint8_t>& taken, const std::vector<std::uint8_t>& stop) {
      std::array<std::array<std::size_t, 2>, kTakenKinds> counts{};
      for (std::size_t index = 0; index < n; ++index) {
        if (mFate[index] == kShortFate && taken[index] % kLowerChi2 != kNoTaken && (stop[index] == kOtherTaken || stop[index] == kLookaheadOther)) {
          ++counts[taken[index] % kLowerChi2][taken[index] / kLowerChi2];
        }
      }
      std::string text = title;
      for (std::size_t kind = 1; kind < kTakenKinds; ++kind) {
        if (counts[kind][0] + counts[kind][1] != 0U) {
          text += std::format(" {} {} ({} with a lower chi2);", kTakenNames[kind], counts[kind][0] + counts[kind][1], counts[kind][1]);
        }
      }
      std::cout << text << '\n';
    };
    if (followed) {
      takenLine("  short, its own row passing, the row taken instead:", mStopTaken, mStop);
    }
    {
      std::array<std::size_t, kShortReasons + 1U> outStops{};
      std::array<std::size_t, kLayers> outLayers{};
      std::size_t outShort = 0U;
      for (std::size_t index = 0; index < n; ++index) {
        if (mFate[index] == kShortFate && mStop[index] == kGapInside) {
          ++outShort;
          ++outStops[mOutStop[index]];
          if (mOutStopLayer[index] < kLayers) {
            ++outLayers[mOutStopLayer[index]];
          }
        }
      }
      line = std::format("  short of an outer layer ({}, {:.2f} %; outward walks replayed {}, {} agree), by the decision at the first it missed:", outShort, pc(outShort), mOutReplayed, mOutAgree);
      for (std::size_t reason = 0; reason <= kShortReasons; ++reason) {
        if (outStops[reason] != 0U) {
          line += std::format(" {} {} ({:.2f} %);", reason < kShortReasons ? kShortNames[reason] : "not replayed", outStops[reason], pc(outStops[reason]));
        }
      }
      std::cout << line << '\n';
      line = "  short of an outer layer, by the first it missed:";
      for (std::size_t layer = 0; layer < kLayers; ++layer) {
        if (outLayers[layer] != 0U) {
          line += std::format(" {}:{}", layer, outLayers[layer]);
        }
      }
      std::cout << line << '\n';
      if (followed) {
        takenLine("  short of an outer layer, its own row passing, the row taken instead:", mOutStopTaken, mOutStop);
        line =
          "  short of an outer layer, refit failed, by status (1 no intersection, 4 grazing, 5 not finite, 6 too many turns, 7 stopped) / pass (1 outward, 2 "
          "inward) / layer:";
        for (const auto& [key, count] : mRefitFailures) {
          line += std::format(" {}: {};", key, count);
        }
        std::cout << line << '\n';
        std::cout << "  short of an outer layer, refit failed, replayed by the outward rows not the particle's (2: the failing step's propagation status; 3: the "
                     "filtered state's circle against the failing plane, barrel planes):\n";
        for (const auto& [key, count] : mRefitReplay) {
          std::cout << std::format("    {}: {}\n", key, count);
        }
        {
          // The refit-chi2 stops: the outward rows the refit took, the wrong
          // ones among them, and the chi2 over the cut.
          std::array<std::array<std::size_t, 6>, 6> bins{};
          std::array<std::size_t, kMaxClusters + 1U> rows{};
          constexpr std::array<double, 5> kEdges{1.5, 2., 3., 5., 10.};
          for (std::size_t index = 0; index < n; ++index) {
            if (mFate[index] == kShortFate && mStop[index] == kGapInside && mOutStop[index] == kRefitChi2) {
              const auto ratio = static_cast<double>(mOutRefitRatio[index]);
              const auto bin = static_cast<std::size_t>(std::ranges::upper_bound(kEdges, ratio) - kEdges.begin());
              ++bins[mOutRefitWrong[index]][bin];
              ++rows[mOutRefitRows[index]];
            }
          }
          line =
            "  short of an outer layer, refit chi2: outward rows not the particle's 0 / 1 / 2+ (k: with a kink step), by chi2 over the cut (1-1.5 / 1.5-2 / 2-3 "
            "/ 3-5 / 5-10 / >10):";
          for (std::size_t wrong = 0; wrong < 6U; ++wrong) {
            line += std::format(" [{}{}] {} {} {} {} {} {};", wrong % 3U, wrong >= 3U ? "k" : "", bins[wrong][0], bins[wrong][1], bins[wrong][2], bins[wrong][3], bins[wrong][4], bins[wrong][5]);
          }
          std::array<std::array<std::size_t, 6>, 4> kinkBins{};
          for (std::size_t index = 0; index < n; ++index) {
            if (mFate[index] == kShortFate && mStop[index] == kGapInside && mOutStop[index] == kRefitChi2 && mOutRefitWrong[index] % 3U == 0U) {
              const auto ratio = static_cast<double>(mOutRefitRatio[index]);
              ++kinkBins[mOutRefitKink[index]][static_cast<std::size_t>(std::ranges::upper_bound(kEdges, ratio) - kEdges.begin())];
            }
          }
          constexpr std::array<const char*, 4> kKinkNames{"<= 3 theta0", "3-5", "> 5", "no hits"};
          line += " | all own, by the largest true kink over the rows:";
          for (std::size_t kink = 0; kink < 4U; ++kink) {
            line += std::format(" [{}] {} {} {} {} {} {};", kKinkNames[kink], kinkBins[kink][0], kinkBins[kink][1], kinkBins[kink][2], kinkBins[kink][3], kinkBins[kink][4], kinkBins[kink][5]);
          }
          line += " outward rows taken:";
          for (std::size_t count = 0; count <= kMaxClusters; ++count) {
            if (rows[count] != 0U) {
              line += std::format(" {}:{}", count, rows[count]);
            }
          }
          std::cout << line << '\n';
        }
        for (const auto& [title, group] : {std::pair<const char*, const Increments*>{"refused (refit chi2, all own rows)", &mRefusedIncrements}, std::pair<const char*, const Increments*>{"control (complete, refit kept)", &mKeptIncrements}}) {
          line = std::format("  refit inward increments by layer, {} {}: layer mean / % > 9.21 / % the largest:", title, group->tracks);
          for (std::size_t layer = 0; layer < kLayers; ++layer) {
            if (group->count[layer] != 0U) {
              const auto c = static_cast<double>(group->count[layer]);
              line += std::format(" {}: {:.2f} / {:.1f} / {:.1f};", layer, group->sum[layer] / c, study::percent(group->above[layer], group->count[layer]), study::percent(group->largest[layer], group->tracks));
            }
          }
          std::cout << line << '\n';
        }
        // Who holds the own row on the first layer missed, by the decision there.
        for (const std::uint8_t reason : {static_cast<std::uint8_t>(kLostBid), static_cast<std::uint8_t>(kShortReasons)}) {
          std::array<std::size_t, kHolderKinds> holders{};
          for (std::size_t index = 0; index < n; ++index) {
            if (mFate[index] == kShortFate && mStop[index] == kGapInside && mOutStop[index] == reason && mOutHolder[index] < kHolderKinds) {
              ++holders[mOutHolder[index]];
            }
          }
          line = std::format("  short of an outer layer, {}: the own row on that layer at the end:", reason < kShortReasons ? kShortNames[reason] : "not replayed");
          for (std::size_t kind = 0; kind < kHolderKinds; ++kind) {
            if (holders[kind] != 0U) {
              line += std::format(" {} {};", kHolderNames[kind], holders[kind]);
            }
          }
          std::cout << line << '\n';
        }
      }
    }
    {
      // The unconfirmed drops' next step, and the chi2 of the own
      // rows refused by the gate, against it.
      constexpr std::size_t kNextOutcomes = kNotVisited + 1U; // the step outcomes and "not visited"
      constexpr std::array<double, 6> kChi2Edges{0., 15., 25., 50., 100., 1.e300};
      const auto chi2Bins = [&](const std::vector<Scalar>& chi2, auto&& select) {
        std::array<std::size_t, kChi2Edges.size() - 1U> bins{};
        for (std::size_t index = 0; index < n; ++index) {
          if (select(index)) {
            for (std::size_t bin = 0; bin + 1U < kChi2Edges.size(); ++bin) {
              const auto value = static_cast<double>(chi2[index]);
              bins[bin] += value >= kChi2Edges[bin] && value < kChi2Edges[bin + 1U] ? 1U : 0U;
            }
          }
        }
        return std::format("gate-15 {} / 15-25 {} / 25-50 {} / 50-100 {} / >100 {}", bins[0], bins[1], bins[2], bins[3], bins[4]);
      };
      if (followed) {
        for (std::size_t outward = 0; outward < 2U; ++outward) {
          const auto& next = outward == 0U ? mInNext : mOutNext;
          const auto& nextChi2 = outward == 0U ? mInNextChi2 : mOutNextChi2;
          std::array<std::size_t, kNextOutcomes> counts{};
          for (std::size_t index = 0; index < n; ++index) {
            if (mFate[index] == kShortFate && next[index] < kNextOutcomes) {
              ++counts[next[index]];
            }
          }
          std::string text = std::format("  {}, unconfirmed drops: the next step with the particle's own row there:", outward == 0U ? "inward" : "outward");
          for (std::size_t outcome = 0; outcome < kNextOutcomes; ++outcome) {
            if (counts[outcome] != 0U) {
              text += std::format(" {} {};", kStepOutcomeNames[outcome], counts[outcome]);
            }
          }
          std::cout << text << '\n';
          std::cout << std::format("    its chi2 where the gate ({:.2f}) refused it: {}\n", static_cast<double>(mLastGate), chi2Bins(nextChi2, [&](std::size_t index) { return mFate[index] == kShortFate && next[index] == static_cast<std::uint8_t>(StepOutcome::kChi2); }));
          const auto& own = outward == 0U ? mInChi2 : mOutChi2;
          const auto& stop = outward == 0U ? mStop : mOutStop;
          std::cout << std::format("    the own row refused by the gate at the deciding step, its chi2: {}\n", chi2Bins(own, [&](std::size_t index) { return mFate[index] == kShortFate && stop[index] == kChi2; }));
          // The particle's true kink over the step that refused its next row:
          // by that step's outcome on the own row.
          const auto& dropKink = outward == 0U ? mInDropKink : mOutDropKink;
          std::string kinkText = "    the particle's true kink over that next step, over theta0 (<= 3 / 3-5 / > 5 / no hits), by the step's outcome on the own row:";
          for (std::size_t outcome = 0; outcome < kNextOutcomes; ++outcome) {
            std::array<std::size_t, 4> bins{};
            for (std::size_t index = 0; index < n; ++index) {
              if (mFate[index] == kShortFate && next[index] == outcome && dropKink[index] < 4U) {
                ++bins[dropKink[index]];
              }
            }
            if (bins[0] + bins[1] + bins[2] + bins[3] != 0U) {
              kinkText += std::format(" {} {}/{}/{}/{};", kStepOutcomeNames[outcome], bins[0], bins[1], bins[2], bins[3]);
            }
          }
          std::cout << kinkText << '\n';
          // The same over the step before, the one that took the dropped row
          // (outward: the step from the layer below the gap, mOutStepKink).
          const auto& stepKink = outward == 0U ? mInStepKink : mOutStepKink;
          kinkText =
            "    the particle's true kink over the step before, the one that took the dropped row (<= 3 / 3-5 / > 5 / no hits), by the next step's outcome on "
            "the own row:";
          for (std::size_t outcome = 0; outcome < kNextOutcomes; ++outcome) {
            std::array<std::size_t, 4> bins{};
            for (std::size_t index = 0; index < n; ++index) {
              if (mFate[index] == kShortFate && next[index] == outcome && stepKink[index] < 4U) {
                ++bins[stepKink[index]];
              }
            }
            if (bins[0] + bins[1] + bins[2] + bins[3] != 0U) {
              kinkText += std::format(" {} {}/{}/{}/{};", kStepOutcomeNames[outcome], bins[0], bins[1], bins[2], bins[3]);
            }
          }
          std::cout << kinkText << '\n';
        }
      }
      // The outward kink step's refusal of the own row, under the scaled
      // gate (the inward one is printed with the short reasons), and the
      // outward decisions by the particle's pT.
      {
        if (followed) {
          std::array<std::size_t, kNextOutcomes> refused{};
          for (std::size_t index = 0; index < n; ++index) {
            if (mFate[index] == kShortFate && mStop[index] == kGapInside && mOutStop[index] == kKinkRefused && mStopDetail[index] % 8U < kNextOutcomes) {
              ++refused[mStopDetail[index] % 8U];
            }
          }
          std::string text = "  short of an outer layer, kink step refused it, under the scaled gate:";
          for (std::size_t outcome = 0; outcome < kNextOutcomes; ++outcome) {
            if (refused[outcome] != 0U) {
              text += std::format(" {} {};", kStepOutcomeNames[outcome], refused[outcome]);
            }
          }
          std::cout << text << '\n';
          // The particle's true kink over the step from the layer below the
          // gap, by the decision there.
          text =
            "  short of an outer layer, by the decision at the first it missed, the particle's true kink over the step from the layer below (<= 3 / 3-5 / > 5 "
            "theta0 / no hits):";
          for (std::size_t reason = 0; reason <= kShortReasons; ++reason) {
            std::array<std::size_t, 4> bins{};
            for (std::size_t index = 0; index < n; ++index) {
              if (mFate[index] == kShortFate && mStop[index] == kGapInside && mOutStop[index] == reason && mOutStepKink[index] < 4U) {
                ++bins[mOutStepKink[index]];
              }
            }
            if (bins[0] + bins[1] + bins[2] + bins[3] != 0U) {
              text += std::format(" {} {}/{}/{}/{};", reason < kShortReasons ? kShortNames[reason] : "not replayed", bins[0], bins[1], bins[2], bins[3]);
            }
          }
          std::cout << text << '\n';
        }
        for (std::size_t bin = 0; bin <= kWalkPtEdges.size(); ++bin) {
          std::array<std::size_t, kShortReasons + 1U> stopsByPt{};
          std::array<std::size_t, kLayers> layersByPt{};
          std::size_t total = 0U;
          std::size_t dropped = 0U;
          for (std::size_t index = 0; index < n; ++index) {
            if (mFate[index] != kShortFate || mStop[index] != kGapInside) {
              continue;
            }
            const double pt = truePt(mTruth->particles.find(mFindable[index]));
            if (std::cmp_not_equal(std::ranges::upper_bound(kWalkPtEdges, pt) - kWalkPtEdges.begin(), bin)) {
              continue;
            }
            ++total;
            ++stopsByPt[mOutStop[index]];
            dropped += mOutStopTaken[index] == kOwnDropped && (mOutStop[index] == kOtherTaken || mOutStop[index] == kLookaheadOther) ? 1U : 0U;
            if (mOutStopLayer[index] < kLayers) {
              ++layersByPt[mOutStopLayer[index]];
            }
          }
          if (total == 0U) {
            continue;
          }
          std::string byPt = std::format("  short of an outer layer, pT {} GeV: {} (the own row taken and dropped after, by the unconfirmed drop: {}); by the decision:", kWalkPtNames[bin], total, dropped);
          for (std::size_t reason = 0; reason <= kShortReasons; ++reason) {
            if (stopsByPt[reason] != 0U) {
              byPt += std::format(" {} {};", reason < kShortReasons ? kShortNames[reason] : "not replayed", stopsByPt[reason]);
            }
          }
          byPt += " | by the first layer missed:";
          for (std::size_t layer = 0; layer < kLayers; ++layer) {
            if (layersByPt[layer] != 0U) {
              byPt += std::format(" {}:{}", layer, layersByPt[layer]);
            }
          }
          std::cout << byPt << '\n';
        }
      }
    }
    if (followed) {
      for (std::size_t kind = 0; kind < 2U; ++kind) {
        for (std::size_t direction = 0; direction < 2U; ++direction) {
          const auto& audit = direction == 0U ? mWindowAudit[kind] : mOutWindowAudit[kind];
          std::cout << std::format("  {}{}, the row against the windows (pair window / covariance window at nSigma / at the gate's sqrt(chi2 max)): ", direction == 0U ? "" : "outward, ", kind == 0U ? "outside the window" : "not queried");
          for (std::size_t flags = 0; flags < 8U; ++flags) {
            if (audit[flags] != 0U) {
              std::cout << std::format("{}{}{} {}; ", (flags & 1U) != 0U ? "P" : "p", (flags & 2U) != 0U ? "N" : "n", (flags & 4U) != 0U ? "G" : "g", audit[flags]);
            }
          }
          std::cout << "(capital: inside)\n";
        }
        // How far outside the gate's window the outward rows lie, per axis.
        auto offsets = mOutWindowOffsets[kind];
        if (!offsets.empty()) {
          std::string offsetLine = std::format("  outward, {}, the row's offset from the window's centre over the gate window's half-width, 50/90/99 % (1/50/99 % signed):", kind == 0U ? "outside the window" : "not queried");
          for (std::size_t axis = 0; axis < 4U; ++axis) {
            std::ranges::sort(offsets, {}, [axis](const std::array<float, 4>& offset) { return offset[axis]; });
            const auto quantile = [&](double q) { return offsets[std::min(offsets.size() - 1U, static_cast<std::size_t>(q * static_cast<double>(offsets.size())))][axis]; };
            constexpr std::array<const char*, 4> kAxes{"azimuth", "second", "azimuth times the charge's sign", "azimuth half-width [um]"};
            const bool signedAxis = axis == 2U;
            offsetLine += std::format(" {} {:.2f} {:.2f} {:.2f};", kAxes[axis], quantile(signedAxis ? 0.01 : 0.5), quantile(signedAxis ? 0.5 : 0.9), quantile(0.99));
          }
          std::cout << offsetLine << '\n';
        }
      }
    }
    line = "  not found (the last iteration's cuts):";
    for (std::size_t reason = 0; reason < kMissedReasons; ++reason) {
      if (missed[reason] != 0U) {
        line += std::format(" {} {} ({:.2f} %);", kMissedNames[reason], missed[reason], pc(missed[reason]));
      }
    }
    std::cout << line << '\n';
    // The same above and below the last iteration's pT floor: a
    // road refused by the triplet's curvature below the floor is the floor.
    {
      const double floor = mFloors.empty() ? 0. : mFloors.back();
      for (const bool above : {true, false}) {
        std::array<std::size_t, kMissedReasons> part{};
        std::size_t total = 0U;
        std::size_t population = 0U;
        for (std::size_t index = 0; index < n; ++index) {
          const double pt = truePt(mTruth->particles.find(mFindable[index]));
          if ((pt >= floor) != above) {
            continue;
          }
          ++population;
          if (mFate[index] == kNotFound) {
            ++total;
            ++part[mReason[index]];
          }
        }
        std::string text = std::format("  not found, pT {} the last floor {:.2f} GeV ({} of {}):", above ? "at or above" : "below", floor, total, population);
        for (std::size_t reason = 0; reason < kMissedReasons; ++reason) {
          if (part[reason] != 0U) {
            text += std::format(" {} {};", kMissedNames[reason], part[reason]);
          }
        }
        // Of those with no road (made, or its rows sharing no time), the
        // first graph predicate refusing the true rows.
        std::array<std::size_t, kGraphReasons + 1U> graph{};
        for (std::size_t index = 0; index < n; ++index) {
          const double pt = truePt(mTruth->particles.find(mFindable[index]));
          if ((pt >= floor) == above && mFate[index] == kNotFound && (mReason[index] == kNoRoad || mReason[index] == kRoadTime)) {
            ++graph[std::min<std::size_t>(mGraph[index], kGraphReasons)];
          }
        }
        text += " | no road made or its rows sharing no time, by the first predicate refusing the true rows:";
        for (std::size_t reason = 0; reason <= kGraphReasons; ++reason) {
          if (graph[reason] != 0U) {
            text += std::format(" {} {};", reason < kGraphReasons ? kGraphNames[reason] : "not audited", graph[reason]);
          }
        }
        // The roads refused whole, by their largest true kink over the
        // road's layers.
        for (const auto reason : {kRoadChi2, kBeamLine}) {
          std::array<std::size_t, 4> bins{};
          for (std::size_t index = 0; index < n; ++index) {
            const double pt = truePt(mTruth->particles.find(mFindable[index]));
            if ((pt >= floor) == above && mFate[index] == kNotFound && mReason[index] == reason) {
              ++bins[mMissedKink[index]];
            }
          }
          if (bins[0] + bins[1] + bins[2] + bins[3] != 0U) {
            text += std::format(" | {} by the road's true kink (<= 3 / 3-5 / > 5 theta0 / no hits): {}/{}/{}/{}", kMissedNames[reason], bins[0], bins[1], bins[2], bins[3]);
          }
        }
        std::cout << text << '\n';
      }
    }
    {
      std::array<std::array<std::size_t, 3>, kGraphReasons + 1U> graph{};
      std::vector<float> worstRoad(mFindable.size(), -1.F);
      // Between the road's layers: the last iteration's.
      const std::uint64_t road = mRoadLayers.empty() ? ~std::uint64_t{0} : mRoadLayers.back();
      for (const auto& kink : mKinks) {
        if (((road >> kink.layer) & 1U) != 0U && ((road >> kink.to) & 1U) != 0U) {
          worstRoad[kink.particle] = std::max({worstRoad[kink.particle], std::abs(kink.theta), std::abs(kink.phi)});
        }
      }
      for (std::size_t index = 0; index < n; ++index) {
        if (mFate[index] == kNotFound && (mReason[index] == kNoRoad || mReason[index] == kRoadTime)) {
          std::size_t kind = worstRoad[index] > 3.F ? 1U : 0U;
          if (worstRoad[index] < 0.F) {
            kind = 2U;
          }
          ++graph[mGraph[index]][kind];
        }
      }
      line = "  no road made or its rows sharing no time, the first graph predicate refusing the true rows (road kink <= 3 theta0 / > 3 / no hits):";
      for (std::size_t reason = 0; reason < kGraphReasons; ++reason) {
        if (graph[reason][0] + graph[reason][1] + graph[reason][2] != 0U) {
          line += std::format(" {} {}/{}/{};", kGraphNames[reason], graph[reason][0], graph[reason][1], graph[reason][2]);
        }
      }
      std::cout << line << '\n';
      if (!mCurvatureTriplets.empty()) {
        std::string text = "  the curvature predicate's refusals of the true rows, by the triplet's layers:";
        for (const auto& [layers, count] : mCurvatureTriplets) {
          text += std::format(" {} {};", layers, count);
        }
        std::cout << text << '\n';
      }
      if (!mLinkRangeRefusals.empty()) {
        // By the link and the side, the distance outside summed into 1-cm
        // bins would be too fine: the map's keys hold it to 0.1 cm, so
        // print the links' totals and the farthest.
        std::map<std::string, std::pair<std::size_t, double>> byLink;
        for (const auto& [key, count] : mLinkRangeRefusals) {
          const auto at = key.find(" by ");
          const auto distance = std::stod(key.substr(at + 4U));
          auto& entry = byLink[key.substr(0, at)];
          entry.first += count;
          entry.second = std::max(entry.second, distance);
        }
        std::string text = "  rows outside their link's source z range, by the link (count, the farthest outside [cm]):";
        for (const auto& [link, entry] : byLink) {
          text += std::format(" {} {} ({:.1f});", link, entry.first, entry.second);
        }
        std::cout << text << '\n';
      }
      const auto sorted = [](const std::map<std::string, std::size_t>& counts) {
        std::vector<std::pair<std::string, std::size_t>> entries(counts.begin(), counts.end());
        std::ranges::stable_sort(entries, std::ranges::greater{}, &std::pair<std::string, std::size_t>::second);
        return entries;
      };
      if (!mFakeHolders.empty()) {
        std::string text = "  road rows held by a fake, by the iteration that made the fake holding the first:";
        for (const auto& [key, count] : sorted(mFakeHolders)) {
          text += std::format(" {} {};", key, count);
        }
        std::cout << text << '\n';
      }
      if (!mPairRefusals.empty()) {
        std::string text = "  the pair tests' refusals of the true rows, by the link and the true kink over it:";
        for (const auto& [key, count] : sorted(mPairRefusals)) {
          text += std::format(" {} {};", key, count);
        }
        std::cout << text << '\n';
      }
      if (!mRoadFitRefusals.empty()) {
        std::cout << "  the road fit's chi2 refusals of true roads at or above the floor, by side, species, rows, the road's true kink, chi2 over the "
                     "kink-free cut (a lower bound), and the fit of the barrel rows and the other rows alone:\n";
        for (const auto& [key, count] : sorted(mRoadFitRefusals)) {
          std::cout << std::format("    {} : {}\n", key, count);
        }
      }
      std::cout << std::format("  the graph audit held to the stages' own cuts: {} triplet and neighbour decisions, {} decided otherwise\n", mAuditDecisions, mAuditDisagreements);
    }
    if (mWalkDone) {
      // The beam spot the class was produced in: its production vertices.
      double sx = 0.;
      double sy = 0.;
      double mx = 0.;
      double my = 0.;
      std::size_t vertices = 0U;
      for (const auto particle : mFindable) {
        const auto found = mTruth->particles.find(particle);
        if (found == mTruth->particles.end()) {
          continue;
        }
        const auto vx = static_cast<double>(found->second.vx);
        const auto vy = static_cast<double>(found->second.vy);
        mx += vx;
        my += vy;
        sx += vx * vx;
        sy += vy * vy;
        ++vertices;
      }
      if (vertices != 0U) {
        const auto c = static_cast<double>(vertices);
        mx /= c;
        my /= c;
        std::vector<double> ax;
        std::vector<double> ay;
        std::size_t beyond = 0U;
        for (const auto particle : mFindable) {
          const auto found = mTruth->particles.find(particle);
          if (found != mTruth->particles.end()) {
            ax.push_back(std::abs(static_cast<double>(found->second.vx)));
            ay.push_back(std::abs(static_cast<double>(found->second.vy)));
            beyond += std::hypot(static_cast<double>(found->second.vx), static_cast<double>(found->second.vy)) > 0.01 ? 1U : 0U;
          }
        }
        // The median |v| of a Gaussian is 0.6745 sigma.
        std::cout << std::format(
          "  the class's production vertices: mean ({:.2f}, {:.2f}) um, rms ({:.2f}, {:.2f}) um, core sigma from the median |v| ({:.2f}, "
          "{:.2f}) um, beyond 100 um from the axis {:.2f} %\n",
          1.e4 * mx, 1.e4 * my, 1.e4 * std::sqrt(std::max(0., (sx / c) - (mx * mx))), 1.e4 * std::sqrt(std::max(0., (sy / c) - (my * my))), 1.e4 * study::median(ax) / 0.6745, 1.e4 * study::median(ay) / 0.6745, study::percent(beyond, vertices));
      }
      constexpr std::array<const char*, 3> kConstructions{"follow", "clusters only", "joint"};
      std::cout << "  truth walk (every findable, its own rows inward whatever their chi2): prediction against the MC crossing, pull rms loc0 / loc1 (% beyond 3)\n";
      constexpr std::array<const char*, 3> kSelections{"all", "quiet", "quiet, prompt"};
      for (std::size_t selection = 0; selection < 3U; ++selection) {
        for (std::size_t layer = 0; layer < kLayers; ++layer) {
          if (mWalk[selection][0][layer].n == 0U) {
            continue;
          }
          std::string text = std::format("    {} target {}:", kSelections[selection], layer);
          for (std::size_t construction = 0; construction < 3U; ++construction) {
            const auto& sum = mWalk[selection][construction][layer];
            if (sum.n == 0U) {
              continue;
            }
            const auto count = static_cast<double>(sum.n);
            text += std::format(" | {} {}: {:.2f} ({:.2f}) / {:.2f} ({:.2f})", kConstructions[construction], sum.n, std::sqrt(sum.s0 / count), study::percent(sum.out0, sum.n), std::sqrt(sum.s1 / count), study::percent(sum.out1, sum.n));
          }
          std::cout << text << '\n';
        }
      }
      // The same walks against the cluster, unselected, the
      // gate's quantity, by target layer and pT: the robust core and the
      // tails.
      std::cout << std::format(
        "  truth walk against the cluster, unselected (by target layer and pT [GeV]): n | median pull loc0 / loc1 | core width "
        "(1.4826 MAD) loc0 / loc1 | % beyond 3 loc0 / loc1 | % chi2 above the gate ({:.2f}) | transport failed\n",
        static_cast<double>(mWalkGate));
      const auto robust = [](std::vector<float> values) {
        std::ranges::sort(values);
        const auto med = static_cast<double>(values[values.size() / 2U]);
        std::vector<float> deviations(values.size());
        for (std::size_t i = 0; i < values.size(); ++i) {
          deviations[i] = std::abs(values[i] - static_cast<float>(med));
        }
        std::ranges::sort(deviations);
        const double width = 1.4826 * static_cast<double>(deviations[deviations.size() / 2U]);
        std::size_t beyond = 0U;
        for (const auto value : values) {
          beyond += std::abs(value) > 3.F ? 1U : 0U;
        }
        return std::tuple{med, width, study::percent(beyond, values.size())};
      };
      for (std::size_t direction = 0; direction < 3U; ++direction) {
        for (std::size_t layer = 0; layer < kLayers; ++layer) {
          for (std::size_t bin = 0; bin <= kWalkPtEdges.size(); ++bin) {
            const auto& pulls = mWalkPulls[direction][layer][bin];
            if (pulls.pull0.empty() && pulls.fails == 0U) {
              continue;
            }
            if (pulls.pull0.empty()) {
              std::cout << std::format("    {} target {} pT {}: 0 | transport failed {}\n", kWalkDirectionNames[direction], layer, kWalkPtNames[bin], pulls.fails);
              continue;
            }
            const auto [med0, width0, beyond0] = robust(pulls.pull0);
            const auto [med1, width1, beyond1] = robust(pulls.pull1);
            std::cout << std::format("    {} target {} pT {}: {} | {:+.3f} / {:+.3f} | {:.3f} / {:.3f} | {:.2f} / {:.2f} | {:.2f} | {}\n", kWalkDirectionNames[direction], layer, kWalkPtNames[bin], pulls.pull0.size(), med0, med1, width0, width1, beyond0, beyond1, study::percent(pulls.overGate, pulls.pull0.size()), pulls.fails);
          }
        }
      }
    }
    for (const auto& text : mTransitionLines) {
      std::cout << text << '\n';
    }
    for (const auto& text : mProgressLines) {
      std::cout << text << '\n';
    }
    {
      std::string text = std::format("  physical primaries in |eta| < 1 on every layer of the mask, of all, by pT (the lowest {:.3f} GeV):", mLowestPt);
      for (std::size_t bin = 0; bin < mReach.size(); ++bin) {
        if (mReach[bin][0] == 0U) {
          continue;
        }
        const std::string range = bin + 1U == mReach.size() ? std::format(">= {:.2f}", 0.05 * static_cast<double>(bin)) : std::format("{:.2f}-{:.2f}", 0.05 * static_cast<double>(bin), 0.05 * static_cast<double>(bin + 1U));
        text += std::format("\n    {} GeV: {} of {} ({:.1f} %)", range, mReach[bin][1], mReach[bin][0], study::percent(mReach[bin][1], mReach[bin][0]));
      }
      std::cout << text << '\n';
    }
    printKinks();
    if (!mPathLines.empty()) {
      std::cout << mPathLines << '\n';
    }
    printInward();
    printClutter();
    printDecisions();
    printBlockers();
  }

  /// Each iteration's pT floor [GeV], for the class's fates by pT range
  /// after each iteration (`progress`).
  void setIterationFloors(std::vector<double> floors) { mFloors = std::move(floors); }

  /// The beam line the tracklet stage's pair cut pairs through, for the
  /// graph audit.
  void setBeamLine(const VertexSeed& beamLine) { mBeamLine = beamLine; }

 private:
  static constexpr std::uint32_t kNoTrack = std::numeric_limits<std::uint32_t>::max();
  static constexpr std::uint32_t kNoRow = ClusterList<kMaxClusters>::kNoCluster;
  enum Fate : std::uint8_t { kComplete,
                             kShortFate,
                             kNotFound };
  static constexpr std::array<const char*, 6> kOutcomeNames{"used", "window", "hole", "noUpdate", "chi2", "passesOtherTaken"};

  struct Candidate {
    std::uint32_t row = 0U;
    StepOutcome outcome = StepOutcome::kPasses;
    Scalar chi2 = 0.F;
    LayerID layer = 0;
  };
  /// A layer a step searched and its window.
  struct TargetSeen {
    LayerID layer = 0;
    CandidateWindow window;
  };
  /// A row a recovery search's confirmation visited for \a row: \a next,
  /// kNoRow for the beam line, and its outcome under the nominal gate.
  struct Witness {
    std::uint32_t row = 0U;
    std::uint32_t next = 0U;
    StepOutcome outcome = StepOutcome::kPasses;
    Scalar chi2 = 0.F;
  };
  /// One step of a replayed proposal: the layer it stepped from, its state
  /// and factors before it, the row it stepped from, what it visited, what
  /// it took; of the recovery pass's searches, the kink budget left where
  /// the kink step was considered (-1: not), whether it looked ahead, the
  /// rows their search visited under its gate and their witnesses, and
  /// whether it took a kink step.
  struct Step {
    LayerID layer = 0;
    BoundTrackState state;
    CovarianceFactors factors;
    std::uint32_t from = 0U;
    std::vector<Candidate> candidates;
    std::vector<TargetSeen> targets;
    std::uint32_t choice = 0U;
    int kinksLeft = -1;
    bool lookahead = false;
    bool kink = false;
    std::vector<Candidate> scaled;
    std::vector<Witness> witnesses;
    /// The lookahead's confirmations at a nominal step more than one row passed:
    /// each passing row, its chi2, the next step's row that
    /// confirmed it (kNoCluster: none) and that row's chi2.
    struct Lookahead {
      std::uint32_t row = 0U;
      Scalar chi2 = 0.F;
      std::uint32_t next = 0U;
      Scalar nextChi2 = 0.F;
    };
    std::vector<Lookahead> lookaheads;
    /// The nominal pass's targets without a row, each with whether
    /// `crossesSensor` called it crossed.
    std::vector<std::pair<LayerID, bool>> holes;
  };
  enum RefitOutcome : std::uint8_t { kRefitKept,
                                     kRefitFailedFit,
                                     kRefitAboveChi2 };
  /// A failed refit's step replayed: the pass run again up to the
  /// measurement it failed on, the filtered state's propagation there tried,
  /// and how far its transverse circle reaches beyond the failing plane along
  /// the plane's transverse normal, and how far its near side stays short of
  /// it [cm]: it crosses the plane only where both are positive.
  struct FailedStep {
    bool replayed = false;
    bool barrel = false; ///< the plane's normal transverse: the margins are measured
    PropagationStatus state = PropagationStatus::kOk;
    float farMargin = 0.F;     ///< the circle's far side beyond the plane
    float nearMargin = 0.F;    ///< the plane beyond the circle's near side
    float planeDistance = 0.F; ///< the plane's distance from the beam axis
    int step = 0;              ///< the failing step, counted from the pass's first measurement
    float startMiss = 0.F;     ///< the pass's seed against its first measurement, on its surface [cm]
  };
  /// One replayed walk of the follow: its raw proposal, its steps, the
  /// rows the bidding kept of it, of them the rows every kept kink step's
  /// witness is among, and what the refit did with them.
  struct Walk {
    ClusterList<kMaxClusters> road{};
    std::uint8_t roadCount = 0U;
    std::vector<std::uint32_t> raw;
    std::vector<Step> steps;
    std::uint64_t witnessed = 0U;
    std::uint32_t bidKept = 0U;
    std::uint32_t kept = 0U;
    RefitOutcome refit = kRefitKept;
    Scalar refitChi2 = 0;    ///< the extended refit's chi2
    Scalar refitChi2Max = 0; ///< and its cut, the road fit's at its ndf
    /// A failed refit: the step's status, its pass and the layer it propagated to.
    PropagationStatus refitStatus = PropagationStatus::kOk;
    FitPass refitFailedPass = FitPass::kNone;
    LayerID refitFailedLayer = kInvalidLayerID;
    /// The failed refit tried again without its outermost outward rows, one
    /// at a time: the outward rows the first fit that succeeds keeps (-1:
    /// none), and whether its chi2 is under the cut.
    int refitTrimmedKept = -1;
    bool refitTrimmedChi2 = false;
    FailedStep failedStep;
    /// The refit's inward pass, the chi2 it decides on: each measurement's
    /// increment and layer, in the chain's order.
    std::array<Scalar, kMaxClusters> refitIncrement{};
    std::array<LayerID, kMaxClusters> refitLayer{};
    std::uint8_t refitCount = 0U;
    /// The extended refit's rows as `replayExtension` fits them: the inward
    /// rows kept, innermost first, then the track's other rows (the road
    /// and the outward rows kept; for a snapshot's replay the track's rows
    /// less its inward proposal). Empty without a refit. The inward study's
    /// chain.
    std::vector<std::uint32_t> refitRows;
    /// How the walk ended, the row the unconfirmed drop removed (kNoRow: none), the last row
    /// refused after it -- at the layout's end (EndDensity) or against the
    /// beam line (BeamEnd) -- and whether its rows
    /// begin with the ones the iteration kept.
    FollowEnd end = FollowEnd::kCapacity;
    std::uint32_t dropped = kNoRow;
    std::uint32_t refused = kNoRow;
    bool agrees = false;
  };
  /// A track's replayed follow: its inward walk, and the outward walk
  /// replayed too: one, its kept rows the iteration's, including a
  /// snapshot's replay.
  struct Replay : Walk {
    std::vector<Walk> outward;
  };

  struct Recorder : follow::NoObserver {
    mutable Walk* replay = nullptr;
    std::uint32_t* from = nullptr;
    Walk* outward = nullptr;
    Walk* inward = nullptr;
    std::uint32_t outerFrom = kNoRow;
    std::uint32_t innerFrom = kNoRow;
    void begin(PropagationDirection direction) const
    {
      if (outward != nullptr) {
        replay = direction == PropagationDirection::kForward ? outward : inward;
        *from = direction == PropagationDirection::kForward ? outerFrom : innerFrom;
      }
    }
    void step(LayerID layer, const BoundTrackState& state, const CovarianceFactors& factors, std::uint32_t /*rows*/) const
    {
      Step one;
      one.layer = layer;
      one.state = state;
      one.factors = factors;
      one.from = *from;
      one.choice = kNoRow;
      replay->steps.push_back(std::move(one));
    }
    void target(LayerID layer, const CandidateWindow& window) const { replay->steps.back().targets.push_back({.layer = layer, .window = window}); }
    void hole(LayerID layer, bool crossed) const { replay->steps.back().holes.emplace_back(layer, crossed); }
    void candidate(LayerID layer, std::uint32_t row, StepOutcome outcome, Scalar chi2) const { replay->steps.back().candidates.push_back({.row = row, .outcome = outcome, .chi2 = chi2, .layer = layer}); }
    void taken(std::uint32_t row, Scalar /*chi2*/, bool kinked) const
    {
      replay->steps.back().choice = row;
      replay->steps.back().kink = kinked;
      *from = row;
    }
    void kink(std::uint32_t kinksLeft) const { replay->steps.back().kinksLeft = static_cast<int>(kinksLeft); }
    void lookahead(std::uint32_t /*passing*/) const { replay->steps.back().lookahead = true; }
    void scaled(LayerID /*layer*/, std::uint32_t row, StepOutcome outcome, Scalar chi2) const { replay->steps.back().scaled.push_back({.row = row, .outcome = outcome, .chi2 = chi2}); }
    void confirmation(std::uint32_t row, bool confirmed, std::uint32_t next, Scalar chi2) const { replay->steps.back().witnesses.push_back({.row = row, .next = next, .outcome = confirmed ? StepOutcome::kPasses : StepOutcome::kChi2, .chi2 = chi2}); }
    void lookaheadConfirmation(std::uint32_t row, Scalar chi2, std::uint32_t next, Scalar nextChi2) const { replay->steps.back().lookaheads.push_back({.row = row, .chi2 = chi2, .next = next, .nextChi2 = nextChi2}); }
    void unconfirmed(std::uint32_t row) const { replay->dropped = row; }
    void refusedEnd(std::uint32_t row) const { replay->refused = row; }
    void refusedBeamEnd(std::uint32_t row) const { replay->refused = row; }
    void end(FollowEnd why) const { replay->end = why; }
  };
  /// The iteration's follow stage, with the same dispatch for both replays.
  template <typename Observation, typename Material, typename Body>
  void withFollowStage(const Observation& iteration, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, const std::uint8_t* used, std::uint32_t first, const Body& body)
  {
    withBools(
      [&]<bool BeamLine, bool Greedy, bool EndDensity, bool BeamEnd>(std::bool_constant<BeamLine>, std::bool_constant<Greedy>, std::bool_constant<EndDensity>, std::bool_constant<BeamEnd>) {
        if constexpr (BeamLine || !BeamEnd) {
          const FollowStage<Det, Material, BeamLine, false, Greedy, EndDensity, BeamEnd> stage{&parameters, &frame, &geometry, &material, propagator, iteration.tracks->view(), used, mMaxWidthPS.data(), mCells.view(), nullptr, nullptr, parameters.sharedInnermost == 0U, first};
          body(stage);
        }
      },
      parameters.beamConstraint, parameters.followGreedy, Follow<Det>::endDensity(parameters), parameters.followBeamEnd);
  }

  /// Both walks through the stage's sequencing: outward first, its narrowed
  /// time box and walked count carried into the inward walk.
  template <typename Stage>
  static auto recordFollow(const Stage& stage, const BoundTrackState& inner, const BoundTrackState& outer, const ClusterList<kMaxClusters>& road, std::uint8_t count, Walk& inward, Walk& outward)
  {
    FollowProposal<kMaxClusters> proposal{};
    follow::KinkSteps steps{};
    std::uint32_t from = road.index[0];
    const Recorder recorder{{}, &inward, &from, &outward, &inward, road.index[count - 1U], road.index[0]};
    stage.follow(inner, outer, road, count, proposal, steps, recorder);
    outward.raw.assign(proposal.outward, proposal.outward + proposal.outwardCount);
    outward.witnessed = steps.outwardWitnessed;
    inward.raw.assign(proposal.inward, proposal.inward + proposal.inwardCount);
    inward.witnessed = steps.inwardWitnessed;
    return proposal;
  }

  /// A track as an iteration left it, for the next one's follow of it.
  struct Snapshot {
    ClusterList<kMaxClusters> clusters{};
    std::uint8_t count = 0U;
    BoundTrackState state;
    BoundTrackState outwardState;
  };

  /// A track's road: its \a rows less those its follow added (\a proposal's
  /// inward and outward rows), all of them without a proposal (no follow).
  /// Not `roadRows`, the truth's road.
  template <typename Proposal>
  [[nodiscard]] static std::vector<std::uint32_t> roadOf(std::span<const std::uint32_t> rows, const Proposal* proposal)
  {
    std::vector<std::uint32_t> road;
    for (const auto row : rows) {
      if (proposal == nullptr || (std::find(proposal->inward, proposal->inward + proposal->inwardCount, row) == proposal->inward + proposal->inwardCount && std::find(proposal->outward, proposal->outward + proposal->outwardCount, row) == proposal->outward + proposal->outwardCount)) {
        road.push_back(row);
      }
    }
    return road;
  }
  /// The true pT [GeV] of the particle \a found points to, 0 for none (the
  /// table's end).
  [[nodiscard]] double truePt(study::ParticleTable::const_iterator found) const { return found != mTruth->particles.end() ? static_cast<double>(found->second.pt) : 0.; }
  /// The production vertex in r-z [cm] of the particle \a found points to,
  /// the origin for none (the table's end): the reference of `pathDistance`.
  [[nodiscard]] std::array<double, 2> vertexRZ(study::ParticleTable::const_iterator found) const
  {
    if (found == mTruth->particles.end()) {
      return {0., 0.};
    }
    return {std::hypot(static_cast<double>(found->second.vx), static_cast<double>(found->second.vy)), static_cast<double>(found->second.vz)};
  }
  /// How far \a row is from a particle's production \a vertex (`vertexRZ`)
  /// in r-z [cm]: the order of "further along the particle's path" for the
  /// fake transitions, the end decisions and the hole calls.
  [[nodiscard]] double pathDistance(std::uint32_t row, const std::array<double, 2>& vertex) const { return std::hypot(static_cast<double>(mTruth->rows[row].r) - vertex[0], static_cast<double>(mTruth->rows[row].z) - vertex[1]); }

  /// Where the iteration's fakes came from: each fake track after the iteration, by its length, as an
  /// earlier track it was before -- a fake kept, a fake lengthened, a good
  /// track made fake by the rows its follow again added -- or as the
  /// iteration's own, its road holding a wrong row or its follow adding
  /// one. Of the rows made wrong by a follow, the first: the step that
  /// took it (a kink step, a step more than one row passed, a step one row
  /// passed, outward) and what that step did with the particle's own row on
  /// its layer.
  template <typename Observation>
  void fakeTransitions(const Observation& iteration, const std::vector<Replay>& replayed, const std::vector<Replay>& recovered)
  {
    const auto view = iteration.tracks->view();
    const auto tracks = iteration.tracks->size();
    enum Kind : std::uint8_t { kFakeKept,
                               kFakeLonger,
                               kMadeFake,
                               kRoadFake,
                               kFollowFake,
                               kKinds };
    enum Source : std::uint8_t { kKinkSource,
                                 kSeveralPassed,
                                 kOnePassed,
                                 kSources };
    constexpr std::array<const char*, kKinds> kKindNames{"earlier fake, kept", "earlier fake, lengthened", "earlier good, made fake", "own: road holds a wrong row", "own: its follow adds one"};
    constexpr std::array<const char*, kSources> kSourceNames{"kink step", "more than one row passed", "one row passed"};
    // The particle's own rows at the step (StepOutcome's names first, the
    // best outcome of its rows on the layer), and the refinements: no
    // own row there with own rows further along the particle's path (by
    // `sortByTraversal`'s distance from its vertex) free or all held, or past
    // its last one; its own row passing with the taken row's chi2 below,
    // equal to or above its own (the lookahead's choice); the walk's replay
    // not agreeing with the rows the iteration kept.
    enum OwnStatus : std::uint8_t { kOwnNotVisited = 6,
                                    kOwnGapFree,
                                    kOwnNotReplayed,
                                    kOwnPastLast,
                                    kOwnPassingTakenHigher,
                                    kOwnPassingEqual,
                                    kOwnGapHeld,
                                    kOwnDisagrees,
                                    kOwnPastLastReturn,
                                    kOwnStatuses };
    constexpr std::array<const char*, kOwnStatuses> kOwnNames{"held", "outside the window", "crossing a layer", "no update", "chi2", "passing, taken lower", "not visited", "no own row, free own rows beyond", "not replayed", "no own row, past its last", "passing, taken higher", "passing, equal", "no own row, own rows beyond all held", "replay disagrees", "no own row, past its last on this leg, own rows on later legs"};
    constexpr std::array<double, 4> kPtEdges{0.2, 0.3, 0.6, 1e9};
    constexpr std::array<const char*, 4> kPtNames{"< 0.2", "0.2-0.3", "0.3-0.6", ">= 0.6"};
    std::array<std::array<std::array<std::size_t, kOwnStatuses>, kPtEdges.size()>, 2> byPt{}; // [inward, outward][pT][own status]
    std::array<std::array<std::size_t, kMaxClusters + 1U>, kKinds> kinds{};
    std::array<std::array<std::array<std::size_t, kOwnStatuses>, kSources>, 2> sources{}; // [inward, outward][source][own status]
    // The same for the class's particles,
    // and their first wrong row's layer by the own row's status there.
    std::array<std::array<std::array<std::size_t, kOwnStatuses>, kSources>, 2> classSources{};
    std::array<std::array<std::size_t, kOwnStatuses>, kLayers> classByLayer{};
    std::array<std::size_t, 4> roadWrong{}; // own road fakes: 0..2 wrong rows, 3 or more
    std::array<std::size_t, kLayers> roadWrongLayer{};
    std::array<std::size_t, 5> roadWhose{};              // the majority particle: the class, a primary missing a layer, one below the floor, not a primary, none
    std::array<std::size_t, 3> roadOwnThere{};           // its own row on a wrong row's layer: none, held by a track, free
    std::array<std::size_t, kLayers> followWrongLayer{}; // follow fakes: the first wrong row's layer
    std::array<std::size_t, kLayers> followRoadInner{};  // follow fakes: the road's innermost layer
    std::array<std::size_t, kLayers> goodRoadInner{};    // the iteration's own good tracks: the road's innermost layer
    std::array<std::size_t, kLayers> roadFakeInner{};    // own road fakes: the road's innermost layer
    std::vector<bool> held(mTruth->rows.size(), false);
    for (std::uint32_t track = 0; track < tracks; ++track) {
      for (std::uint8_t i = 0; i < view.count[track]; ++i) {
        held[view.clusters[track].index[i]] = true;
      }
    }
    const auto holderOf = [&](std::uint32_t row) { return static_cast<bool>(held[row]); };
    const auto wrongRow = [&](std::uint64_t particle, std::span<const std::uint32_t> rows) -> std::uint32_t {
      for (const auto row : rows) {
        if (!mTruth->rows[row].carries(particle)) {
          return row;
        }
      }
      return kNoTrack;
    };
    // The step of \a replay that took \a row, and the particle's own rows there:
    // the inward walk's, or with \a outward the outward walk's.
    const auto attribute = [&](const Replay* replay, std::uint64_t particle, std::uint32_t row, bool outward) {
      const Walk* walk = replay;
      if (outward) {
        walk = replay != nullptr && !replay->outward.empty() ? &replay->outward.front() : nullptr;
      }
      const auto found = mTruth->particles.find(particle);
      const double pt = truePt(found);
      const auto bin = static_cast<std::size_t>(std::ranges::upper_bound(kPtEdges, pt) - kPtEdges.begin());
      auto& counts = sources[outward ? 1U : 0U];
      auto& binned = byPt[outward ? 1U : 0U][std::min(bin, kPtEdges.size() - 1U)];
      const bool inClass = mIndex.contains(particle);
      const auto count = [&](Source source, std::size_t status) {
        ++counts[source][status];
        ++binned[status];
        if (inClass) {
          ++classSources[outward ? 1U : 0U][source][status];
          ++classByLayer[mTruth->rows[row].layer][status];
        }
      };
      if (walk == nullptr) {
        count(kOnePassed, kOwnNotReplayed);
        return;
      }
      if (!walk->agrees) {
        count(kOnePassed, kOwnDisagrees);
        return;
      }
      const auto step = std::ranges::find(walk->steps, row, &Step::choice);
      if (step == walk->steps.end()) {
        count(kOnePassed, kOwnNotReplayed);
        return;
      }
      std::size_t passed = 0U;
      for (const auto& candidate : step->candidates) {
        passed += candidate.outcome == StepOutcome::kPasses ? 1U : 0U;
      }
      Source source = passed > 1U ? kSeveralPassed : kOnePassed;
      if (step->kink) {
        source = kKinkSource;
      }
      const auto layer = mTruth->rows[row].layer;
      const auto& seen = step->kink ? step->scaled : step->candidates;
      const auto taken = std::ranges::find(seen, row, &Candidate::row);
      // Every own row on the layer: the best outcome among them (passing
      // first, then by StepOutcome's order from the end).
      const auto& ownRows = mRowsOf.at(particle);
      bool onLayer = false;
      const Candidate* best = nullptr;
      for (const auto own : ownRows) {
        if (mTruth->rows[own].layer != layer) {
          continue;
        }
        onLayer = true;
        const auto at = std::ranges::find(seen, own, &Candidate::row);
        if (at != seen.end() && (best == nullptr || static_cast<int>(at->outcome) > static_cast<int>(best->outcome) || (at->outcome == best->outcome && at->chi2 < best->chi2))) {
          best = &*at;
        }
      }
      std::size_t status = kOwnNotVisited;
      // The step's row in the particle's path order: where
      // it is the particle's and its rows are timed, the rows beyond it are
      // the later ones (outward) or the earlier (inward), and those on its
      // own leg are what the walk could still take.
      const auto inClassIndex = mIndex.find(particle);
      const auto fromAt = std::ranges::find(ownRows, step->from);
      const bool pathKnown = inClassIndex != mIndex.end() && !mLegs.empty() && mTimed[inClassIndex->second] != 0U && fromAt != ownRows.end();
      if (!onLayer && pathKnown) {
        const auto& legs = mLegs[inClassIndex->second];
        const auto from = static_cast<std::size_t>(fromAt - ownRows.begin());
        bool sameLeg = false;
        bool otherLeg = false;
        bool free = false;
        for (std::size_t k = 0; k < ownRows.size(); ++k) {
          if (outward ? k <= from : k >= from) {
            continue;
          }
          if (legs[k] == legs[from]) {
            sameLeg = true;
            free = free || (ownRows[k] < mFollowUsedNow.size() && mFollowUsedNow[ownRows[k]] == kFree);
          } else {
            otherLeg = true;
          }
        }
        status = sameLeg ? (free ? std::size_t{kOwnGapFree} : std::size_t{kOwnGapHeld}) : (otherLeg ? std::size_t{kOwnPastLastReturn} : std::size_t{kOwnPastLast});
      } else if (!onLayer) {
        // Further along the particle's path than the taken row, in the
        // walk's direction: free at the follow, or all held.
        const auto vertex = vertexRZ(found);
        const auto distance = [&](std::uint32_t one) { return pathDistance(one, vertex); };
        const double here = distance(row);
        bool beyond = false;
        bool free = false;
        for (const auto own : ownRows) {
          if (outward ? distance(own) > here : distance(own) < here) {
            beyond = true;
            free = free || (own < mFollowUsedNow.size() && mFollowUsedNow[own] == kFree);
          }
        }
        status = !beyond ? std::size_t{kOwnPastLast} : (free ? std::size_t{kOwnGapFree} : std::size_t{kOwnGapHeld});
      } else if (best != nullptr) {
        status = static_cast<std::size_t>(best->outcome);
        if (best->outcome == StepOutcome::kPasses && taken != seen.end()) {
          status = taken->chi2 > best->chi2 ? std::size_t{kOwnPassingTakenHigher} : (taken->chi2 < best->chi2 ? static_cast<std::size_t>(StepOutcome::kPasses) : std::size_t{kOwnPassingEqual});
        }
      }
      count(source, status);
    };
    for (std::uint32_t track = 0; track < tracks; ++track) {
      const std::span<const std::uint32_t> rows{view.clusters[track].index, view.count[track]};
      if (study::goodParticle(mTruth->rows, rows) != study::kNoTrackParticle) {
        if (track >= iteration.trackOffset) {
          const auto element = track - iteration.trackOffset;
          const auto* proposal = element < iteration.proposals.size() ? &iteration.proposals[element] : nullptr;
          std::uint8_t innermost = std::numeric_limits<std::uint8_t>::max();
          for (const auto row : roadOf(rows, proposal)) {
            innermost = std::min(innermost, static_cast<std::uint8_t>(mTruth->rows[row].layer));
          }
          if (innermost < kLayers) {
            ++goodRoadInner[innermost];
          }
        }
        continue;
      }
      const auto length = static_cast<std::size_t>(view.count[track]);
      if (track < iteration.trackOffset) {
        // An earlier track: as the previous iteration left it.
        if (track >= mSnapshot.size()) {
          continue;
        }
        const std::span<const std::uint32_t> before{mSnapshot[track].clusters.index, mSnapshot[track].count};
        const auto particle = study::goodParticle(mTruth->rows, before);
        if (particle == study::kNoTrackParticle) {
          ++kinds[before.size() == rows.size() ? kFakeKept : kFakeLonger][length];
          continue;
        }
        ++kinds[kMadeFake][length];
        const auto* proposal = track < iteration.earlierProposals.size() ? &iteration.earlierProposals[track] : nullptr;
        if (proposal == nullptr) {
          continue;
        }
        const auto inward = wrongRow(particle, {proposal->inward, proposal->inwardCount});
        attribute(track < recovered.size() ? &recovered[track] : nullptr, particle, inward != kNoTrack ? inward : wrongRow(particle, {proposal->outward, proposal->outwardCount}), inward == kNoTrack);
        continue;
      }
      // The iteration's own: its road, the rows its follow did not add.
      const auto element = track - iteration.trackOffset;
      // Without a follow (no proposal) the track is its road.
      const auto* proposal = element < iteration.proposals.size() ? &iteration.proposals[element] : nullptr;
      const auto road = roadOf(rows, proposal);
      const auto particle = study::goodParticle(mTruth->rows, road);
      if (particle == study::kNoTrackParticle) {
        ++kinds[kRoadFake][length];
        // Its wrong rows, where they are, and whose track it mostly is.
        const auto majority = study::majorityParticle(mTruth->rows, rows);
        ++roadWrong[std::min<std::size_t>(rows.size() - majority.carried, 3U)];
        ++roadFakeInner[mTruth->rows[std::ranges::min(road, {}, [&](std::uint32_t row) { return mTruth->rows[row].layer; })].layer];
        for (const auto row : road) {
          if (!mTruth->rows[row].carries(majority.particle)) {
            ++roadWrongLayer[mTruth->rows[row].layer];
            const auto own = ownRowOn(majority.particle, mTruth->rows[row].layer);
            ++roadOwnThere[own == kNoTrack ? 0U : (holderOf(own) ? 1U : 2U)];
          }
        }
        const auto found = mTruth->particles.find(majority.particle);
        std::size_t whose = 4U;
        if (mIndex.contains(majority.particle)) {
          whose = 0U;
        } else if (found != mTruth->particles.end() && found->second.origin != study::ParticleOrigin::kPhysicalPrimary) {
          whose = 3U;
        } else if (found != mTruth->particles.end()) {
          whose = static_cast<double>(found->second.pt) < mMinPtGeV ? 2U : 1U;
        }
        ++roadWhose[whose];
        continue;
      }
      ++kinds[kFollowFake][length];
      const auto inward = wrongRow(particle, {proposal->inward, proposal->inwardCount});
      const auto first = inward != kNoTrack ? inward : wrongRow(particle, {proposal->outward, proposal->outwardCount});
      if (first != kNoTrack) {
        ++followWrongLayer[mTruth->rows[first].layer];
      }
      const auto innermost = std::ranges::min(road, {}, [&](std::uint32_t row) { return mTruth->rows[row].layer; });
      ++followRoadInner[mTruth->rows[innermost].layer];
      attribute(element < replayed.size() ? &replayed[element] : nullptr, particle, inward != kNoTrack ? inward : wrongRow(particle, {proposal->outward, proposal->outwardCount}), inward == kNoTrack);
    }
    std::string text = std::format("  fakes after iteration {}, by where they came from (by length):", iteration.iteration);
    for (std::size_t kind = 0; kind < kKinds; ++kind) {
      std::size_t total = 0U;
      std::string lengths;
      for (std::size_t length = 0; length <= kMaxClusters; ++length) {
        total += kinds[kind][length];
        if (kinds[kind][length] != 0U) {
          lengths += std::format(" {}:{}", length, kinds[kind][length]);
        }
      }
      if (total != 0U) {
        text += std::format("\n    {} {} ({} )", kKindNames[kind], total, lengths);
      }
    }
    if (std::ranges::any_of(kinds[kRoadFake], [](std::size_t count) { return count != 0U; })) {
      text += std::format("\n    own road fakes: wrong rows 1: {}, 2: {}, 3+: {}; their road's wrong rows by layer:", roadWrong[1], roadWrong[2], roadWrong[3]);
      for (std::size_t layer = 0; layer < kLayers; ++layer) {
        if (roadWrongLayer[layer] != 0U) {
          text += std::format(" {}:{}", layer, roadWrongLayer[layer]);
        }
      }
      text += std::format(
        "; mostly the track of: the class {}, a physical primary missing a layer of the mask {}, one below the pT floor {}, not a primary {}, "
        "unknown {}; its own row on a wrong row's layer: none {}, held {}, free {}",
        roadWhose[0], roadWhose[1], roadWhose[2], roadWhose[3], roadWhose[4], roadOwnThere[0], roadOwnThere[1], roadOwnThere[2]);
      text += "; the road's innermost layer:";
      for (std::size_t layer = 0; layer < kLayers; ++layer) {
        if (roadFakeInner[layer] != 0U) {
          text += std::format(" {}:{}", layer, roadFakeInner[layer]);
        }
      }
    }
    if (std::ranges::any_of(kinds[kFollowFake], [](std::size_t count) { return count != 0U; })) {
      text += "\n    own follow fakes: the first wrong row's layer:";
      for (std::size_t layer = 0; layer < kLayers; ++layer) {
        if (followWrongLayer[layer] != 0U) {
          text += std::format(" {}:{}", layer, followWrongLayer[layer]);
        }
      }
      text += "; the road's innermost layer:";
      for (std::size_t layer = 0; layer < kLayers; ++layer) {
        if (followRoadInner[layer] != 0U) {
          text += std::format(" {}:{}", layer, followRoadInner[layer]);
        }
      }
    }
    text += "\n    own good tracks, the road's innermost layer:";
    for (std::size_t layer = 0; layer < kLayers; ++layer) {
      if (goodRoadInner[layer] != 0U) {
        text += std::format(" {}:{}", layer, goodRoadInner[layer]);
      }
    }
    text += "\n    the first wrong row a follow added, by the walk and the step that took it | the particle's own row at that step:";
    for (std::size_t walk = 0; walk < 2U; ++walk) {
      for (std::size_t source = 0; source < kSources; ++source) {
        std::size_t total = 0U;
        std::string own;
        for (std::size_t status = 0; status < kOwnNames.size(); ++status) {
          total += sources[walk][source][status];
          if (sources[walk][source][status] != 0U) {
            own += std::format(" {} {};", kOwnNames[status], sources[walk][source][status]);
          }
        }
        if (total != 0U) {
          text += std::format("\n      {}, {} {} |{}", walk == 0U ? "inward" : "outward", kSourceNames[source], total, own);
        }
      }
    }
    text += "\n    the same by the particle's pT:";
    for (std::size_t walk = 0; walk < 2U; ++walk) {
      for (std::size_t bin = 0; bin < kPtEdges.size(); ++bin) {
        std::size_t total = 0U;
        std::string own;
        for (std::size_t status = 0; status < kOwnNames.size(); ++status) {
          total += byPt[walk][bin][status];
          if (byPt[walk][bin][status] != 0U) {
            own += std::format(" {} {};", kOwnNames[status], byPt[walk][bin][status]);
          }
        }
        if (total != 0U) {
          text += std::format("\n      {}, pT {} GeV {} |{}", walk == 0U ? "inward" : "outward", kPtNames[bin], total, own);
        }
      }
    }
    text += "\n    the same for the class's particles:";
    for (std::size_t walk = 0; walk < 2U; ++walk) {
      for (std::size_t source = 0; source < kSources; ++source) {
        std::size_t total = 0U;
        std::string own;
        for (std::size_t status = 0; status < kOwnNames.size(); ++status) {
          total += classSources[walk][source][status];
          if (classSources[walk][source][status] != 0U) {
            own += std::format(" {} {};", kOwnNames[status], classSources[walk][source][status]);
          }
        }
        if (total != 0U) {
          text += std::format("\n      {}, {} {} |{}", walk == 0U ? "inward" : "outward", kSourceNames[source], total, own);
        }
      }
    }
    text += "\n    the class's particles, by the first wrong row's layer | the particle's own row at that step:";
    for (std::size_t layer = 0; layer < kLayers; ++layer) {
      std::size_t total = 0U;
      std::string own;
      for (std::size_t status = 0; status < kOwnNames.size(); ++status) {
        total += classByLayer[layer][status];
        if (classByLayer[layer][status] != 0U) {
          own += std::format(" {} {};", kOwnNames[status], classByLayer[layer][status]);
        }
      }
      if (total != 0U) {
        text += std::format("\n      layer {} {} |{}", layer, total, own);
      }
    }
    mTransitionLines.push_back(text);
  }

  /// What the unconfirmed drop saw where it removed a row,
  /// over the walks of the iteration's own good roads: whether the dropped
  /// row was the road's particle's, and what the step that found nothing
  /// (the walk's last) had seen: a row of the particle or of another within
  /// the window but above the gate (a failed update does not count), a row
  /// passing the kink pass's scaled gate (unconfirmed, or the pass not made
  /// for want of a budget or an allowance), or nothing of either. The
  /// strongest kind is counted: scaled before gated, the particle's before
  /// another's.
  template <typename Observation>
  void endEvidence(const Observation& iteration, const std::vector<Replay>& replayed)
  {
    enum Seen : std::uint8_t { kSeenOwnScaled,
                               kSeenOtherScaled,
                               kSeenOwnGated,
                               kSeenOtherGated,
                               kSeenNone,
                               kSeenNoKinkPass,
                               kSeenKinds };
    constexpr std::array<const char*, kSeenKinds> kSeenNames{"the particle's row passing the scaled gate", "another row passing the scaled gate", "the particle's row within the window, above the gate", "another row within the window, above the gate", "nothing in the window", "nothing in the window, no kink pass"};
    constexpr std::array<double, 4> kPtEdges{0.3, 0.6, 1., 1e9};
    constexpr std::array<const char*, 4> kPtNames{"< 0.3", "0.3-0.6", "0.6-1", ">= 1"};
    // [walk][dropped the particle's / another's][pT][seen]
    std::array<std::array<std::array<std::array<std::size_t, kSeenKinds>, kPtEdges.size()>, 2>, 2> counts{};
    for (const auto& replay : replayed) {
      const std::span<const std::uint32_t> road{replay.road.index, replay.roadCount};
      const auto particle = study::goodParticle(mTruth->rows, road);
      if (particle == study::kNoTrackParticle) {
        continue;
      }
      const auto found = mTruth->particles.find(particle);
      const double pt = truePt(found);
      const auto bin = std::min(static_cast<std::size_t>(std::ranges::upper_bound(kPtEdges, pt) - kPtEdges.begin()), kPtEdges.size() - 1U);
      const auto one = [&](const Walk& walk, std::size_t direction) {
        // Agreeing walks only; a gated row is one the update was made for
        // (a failed update is no evidence); the scaled evidence exists only
        // where the kink pass ran (a spent budget shows as none).
        if (walk.dropped == kNoRow || walk.steps.empty() || !walk.agrees) {
          return;
        }
        const auto& last = walk.steps.back();
        bool ownScaled = false;
        bool otherScaled = false;
        bool ownGated = false;
        bool otherGated = false;
        for (const auto& candidate : last.scaled) {
          if (candidate.outcome == StepOutcome::kPasses) {
            (mTruth->rows[candidate.row].carries(particle) ? ownScaled : otherScaled) = true;
          }
        }
        for (const auto& candidate : last.candidates) {
          if (candidate.outcome == StepOutcome::kChi2) {
            (mTruth->rows[candidate.row].carries(particle) ? ownGated : otherGated) = true;
          }
        }
        const std::size_t seen = ownScaled ? kSeenOwnScaled : (otherScaled ? kSeenOtherScaled : (ownGated ? kSeenOwnGated : (otherGated ? kSeenOtherGated : (last.kinksLeft < 0 ? kSeenNoKinkPass : kSeenNone))));
        ++counts[direction][mTruth->rows[walk.dropped].carries(particle) ? 0U : 1U][bin][seen];
      };
      one(replay, 0U);
      if (!replay.outward.empty()) {
        one(replay.outward.front(), 1U);
      }
    }
    std::string text = std::format("  the unconfirmed drops after iteration {} (the walks of the iteration's own good roads), by the dropped row and what the refusing step saw:", iteration.iteration);
    for (std::size_t direction = 0; direction < 2U; ++direction) {
      for (std::size_t own = 0; own < 2U; ++own) {
        for (std::size_t bin = 0; bin < kPtEdges.size(); ++bin) {
          std::size_t total = 0U;
          std::string seen;
          for (std::size_t kind = 0; kind < kSeenKinds; ++kind) {
            total += counts[direction][own][bin][kind];
            if (counts[direction][own][bin][kind] != 0U) {
              seen += std::format(" {} {};", kSeenNames[kind], counts[direction][own][bin][kind]);
            }
          }
          if (total != 0U) {
            text += std::format("\n    {}, dropped {} row, pT {} GeV {} |{}", direction == 0U ? "inward" : "outward", own == 0U ? "the particle's" : "another's", kPtNames[bin], total, seen);
          }
        }
      }
    }
    mTransitionLines.push_back(text);
  }

  /// The lookahead against the step's own chi2
  /// ranking, paired, at every nominal step of the iteration's own tracks'
  /// agreeing walks where more than one row passed and the road's particle
  /// is among them: whether the lookahead's choice and the smallest chi2 are the
  /// particle's row (both, the lookahead corrects, the lookahead spoils, neither), by the rows
  /// passing and the particle's pT; where they differ, what confirmed each
  /// (the particle's next row, another row, none).
  template <typename Observation>
  void lookaheadPairs(const Observation& iteration, const std::vector<Replay>& replayed)
  {
    enum Pair : std::uint8_t { kBoth,
                               kCorrects,
                               kSpoils,
                               kNeither,
                               kPairs };
    constexpr std::array<const char*, kPairs> kPairNames{"both the particle's", "the lookahead corrects", "the lookahead spoils", "neither"};
    constexpr std::array<const char*, 3> kPassing{"2 passing", "3 passing", "4+ passing"};
    constexpr std::array<double, 4> kPtEdges{0.2, 0.3, 0.6, 1e9};
    constexpr std::array<const char*, 4> kPtNames{"< 0.2", "0.2-0.3", "0.3-0.6", ">= 0.6"};
    constexpr std::array<const char*, 3> kWitness{"own next", "other next", "none"};
    std::array<std::array<std::size_t, kPairs>, 3> byPassing{};
    std::array<std::array<std::size_t, kPairs>, 4> byPt{};
    // Where the lookahead and chi2 differ: [corrects, spoils][the lookahead's choice's witness][the chi2 choice's witness]
    std::array<std::array<std::array<std::size_t, 3>, 3>, 2> witnesses{};
    std::size_t noOwn = 0U;
    std::size_t steps = 0U;
    for (std::uint32_t element = 0; element < replayed.size(); ++element) {
      const auto& replay = replayed[element];
      if (replay.roadCount == 0U) {
        continue;
      }
      const auto particle = study::goodParticle(mTruth->rows, std::span<const std::uint32_t>{replay.road.index, replay.roadCount});
      if (particle == study::kNoTrackParticle) {
        continue;
      }
      const auto found = mTruth->particles.find(particle);
      const double pt = truePt(found);
      const auto ptBin = std::min(static_cast<std::size_t>(std::ranges::upper_bound(kPtEdges, pt) - kPtEdges.begin()), kPtEdges.size() - 1U);
      const auto witnessOf = [&](const Step::Lookahead& one) { return one.next == ClusterList<kMaxClusters>::kNoCluster ? 2U : (mTruth->rows[one.next].carries(particle) ? 0U : 1U); };
      const auto walks = [&](const Walk& walk) {
        if (!walk.agrees) {
          return;
        }
        for (const auto& step : walk.steps) {
          if (step.kink || step.lookaheads.size() < 2U || step.choice == kNoRow) {
            continue;
          }
          ++steps;
          const auto r2 = std::ranges::find(step.lookaheads, step.choice, &Step::Lookahead::row);
          const auto chi2 = std::ranges::min_element(step.lookaheads, [](const Step::Lookahead& a, const Step::Lookahead& b) { return a.chi2 < b.chi2 || (!(b.chi2 < a.chi2) && a.row < b.row); });
          const bool ownPassing = std::ranges::any_of(step.lookaheads, [&](const Step::Lookahead& one) { return mTruth->rows[one.row].carries(particle); });
          if (!ownPassing || r2 == step.lookaheads.end()) {
            noOwn += ownPassing ? 0U : 1U;
            continue;
          }
          const bool r2Own = mTruth->rows[r2->row].carries(particle);
          const bool chi2Own = mTruth->rows[chi2->row].carries(particle);
          const Pair pair = r2Own ? (chi2Own ? kBoth : kCorrects) : (chi2Own ? kSpoils : kNeither);
          ++byPassing[std::min<std::size_t>(step.lookaheads.size() - 2U, 2U)][pair];
          ++byPt[ptBin][pair];
          if (pair == kCorrects || pair == kSpoils) {
            ++witnesses[pair == kCorrects ? 0U : 1U][witnessOf(*r2)][witnessOf(*chi2)];
          }
        }
      };
      walks(replay);
      for (const auto& outward : replay.outward) {
        walks(outward);
      }
    }
    std::string text = std::format(
      "  the lookahead against the smallest chi2 after iteration {} (nominal steps of agreeing walks, more than one row passing: {}, "
      "{} without the particle's row among them):",
      iteration.iteration, steps, noOwn);
    const auto line = [&](const char* name, const std::array<std::size_t, kPairs>& counts) {
      std::string one = std::format("\n    {}:", name);
      for (std::size_t pair = 0; pair < kPairs; ++pair) {
        one += std::format(" {} {};", kPairNames[pair], counts[pair]);
      }
      return one;
    };
    for (std::size_t k = 0; k < kPassing.size(); ++k) {
      text += line(kPassing[k], byPassing[k]);
    }
    for (std::size_t k = 0; k < kPtNames.size(); ++k) {
      text += line(std::format("pT {} GeV", kPtNames[k]).c_str(), byPt[k]);
    }
    for (std::size_t kind = 0; kind < 2U; ++kind) {
      text += std::format("\n    where the lookahead {}, the lookahead's choice confirmed by / the smallest chi2's by:", kind == 0U ? "corrects" : "spoils");
      for (std::size_t a = 0; a < 3U; ++a) {
        for (std::size_t b = 0; b < 3U; ++b) {
          if (witnesses[kind][a][b] != 0U) {
            text += std::format(" {} / {} {};", kWitness[a], kWitness[b], witnesses[kind][a][b]);
          }
        }
      }
    }
    mTransitionLines.push_back(text);
  }

  /// End rules replayed on the iteration's
  /// own tracks, against what the iteration did (the unconfirmed drop included): (a) every
  /// direction's last added row dropped; (b) the last added row dropped
  /// where one row alone passed its step and the walk ended searching
  /// (a hole crossed or nothing found), not at the layout's end. Only a
  /// walk's last row, kept whole by the bid; agreeing walks. For each: the
  /// rows dropped (the particle's or not), fakes made good, good tracks
  /// losing a true row, tracks under the length cut (good, fake).
  template <typename Observation>
  void endPolicies(const Observation& iteration, const std::vector<Replay>& replayed, const IterationParameters<Det>& parameters)
  {
    const auto view = iteration.tracks->view();
    enum Rule : std::uint8_t { kTrimLast,
                               kLoneUnconfirmed,
                               kRules };
    constexpr std::array<const char*, kRules> kRuleNames{"(a) the last row dropped", "(b) a lone last row dropped where the walk ended searching"};
    struct Counts {
      std::size_t tracks = 0U, wrongDropped = 0U, ownDropped = 0U, fakeToGood = 0U, goodLosesRow = 0U, goodUnderCut = 0U, fakeUnderCut = 0U;
    };
    std::array<Counts, kRules> counts{};
    std::array<std::array<std::size_t, 4>, 2> ends{}; // the walks' ends, [last row the particle's, not][FollowEnd]
    for (std::uint32_t element = 0; element < replayed.size() && element < iteration.proposals.size(); ++element) {
      const auto track = iteration.trackOffset + element;
      const auto& proposal = iteration.proposals[element];
      const auto& replay = replayed[element];
      const std::span<const std::uint32_t> rows{view.clusters[track].index, view.count[track]};
      const auto before = study::goodParticle(mTruth->rows, rows);
      const auto roadParticle = replay.roadCount == 0U ? study::kNoTrackParticle : study::goodParticle(mTruth->rows, std::span<const std::uint32_t>{replay.road.index, replay.roadCount});
      for (std::size_t rule = 0; rule < kRules; ++rule) {
        std::vector<std::uint32_t> dropped;
        const auto consider = [&](const Walk& walk, const std::uint32_t* kept, std::uint8_t keptCount) {
          if (!walk.agrees || keptCount == 0U || keptCount != walk.raw.size()) {
            return;
          }
          const auto last = kept[keptCount - 1U];
          if (rule == 0U) {
            ++ends[roadParticle != study::kNoTrackParticle && mTruth->rows[last].carries(roadParticle) ? 0U : 1U][static_cast<std::size_t>(walk.end)];
          }
          const auto step = std::ranges::find(walk.steps, last, &Step::choice);
          std::size_t passed = 0U;
          if (step != walk.steps.end()) {
            for (const auto& candidate : step->candidates) {
              passed += candidate.outcome == StepOutcome::kPasses ? 1U : 0U;
            }
          }
          const bool lone = step != walk.steps.end() && !step->kink && passed == 1U;
          const bool searching = walk.end == FollowEnd::kHoles || walk.end == FollowEnd::kNoRow;
          if (rule == kTrimLast || (lone && searching)) {
            dropped.push_back(last);
          }
        };
        if (!replay.outward.empty()) {
          consider(replay.outward.front(), proposal.outward, proposal.outwardCount);
        }
        consider(replay, proposal.inward, proposal.inwardCount);
        if (dropped.empty()) {
          continue;
        }
        auto& c = counts[rule];
        ++c.tracks;
        std::vector<std::uint32_t> after;
        for (const auto row : rows) {
          if (std::ranges::find(dropped, row) == dropped.end()) {
            after.push_back(row);
          }
        }
        const auto particle = before != study::kNoTrackParticle ? before : roadParticle;
        for (const auto row : dropped) {
          const bool own = particle != study::kNoTrackParticle && mTruth->rows[row].carries(particle);
          c.ownDropped += own ? 1U : 0U;
          c.wrongDropped += own ? 0U : 1U;
        }
        const bool under = parameters.minClustersAfterFollow != 0U && after.size() < parameters.minClustersAfterFollow;
        const auto afterParticle = study::goodParticle(mTruth->rows, after);
        if (before == study::kNoTrackParticle) {
          c.fakeUnderCut += under ? 1U : 0U;
          c.fakeToGood += !under && afterParticle != study::kNoTrackParticle ? 1U : 0U;
        } else {
          c.goodUnderCut += under ? 1U : 0U;
          c.goodLosesRow += !under ? 1U : 0U;
        }
      }
    }
    std::string text = std::format("  end rules after iteration {} (against the iteration as it ran, the unconfirmed drop included):", iteration.iteration);
    for (std::size_t own = 0; own < 2U; ++own) {
      text += std::format("\n    walks whose kept last row is {}, by how they ended:", own == 0U ? "the road's particle's" : "another's");
      for (std::size_t end = 0; end < kFollowEndNames.size(); ++end) {
        text += std::format(" {} {};", kFollowEndNames[end], ends[own][end]);
      }
    }
    for (std::size_t rule = 0; rule < kRules; ++rule) {
      const auto& c = counts[rule];
      text += std::format(
        "\n    {}: {} tracks; rows dropped {} not the particle's, {} its own; fakes made good {}, fakes under the length cut {}; good "
        "tracks losing a true row {}, good tracks under the length cut {}",
        kRuleNames[rule], c.tracks, c.wrongDropped, c.ownDropped, c.fakeToGood, c.fakeUnderCut, c.goodLosesRow, c.goodUnderCut);
    }
    mTransitionLines.push_back(text);
  }

  /// Every walk of the iteration's own tracks at its end, one
  /// CSV line each: the track, the direction, the road's particle (0: the
  /// road is not one particle's) and its pT; how the walk ended
  /// (`FollowEnd`), whether its last step called a sensor crossed (the
  /// unconfirmed drop's condition) and how many links leave the layer it ended on that way
  /// (none: the layout's end; else a `noTarget` is a failed prediction);
  /// the unconfirmed drop's verdict -- the last row dropped, kept on evidence (the keep's chi2
  /// or inward's near miss), kept unconfirmed (no sensor crossed), or no row
  /// walked; the rows walked, kept by the bid, kept with the witnesses; the
  /// last row walked (with whether a kink step took it) and the row the unconfirmed drop
  /// dropped, each with its layer, chi2 and what it is to the road's particle: its
  /// own, noise, past the end (no own row on that layer; then its own rows
  /// further along its path in the walk's direction, by the distance from
  /// its vertex as `attribute` measures it: free at the follow, all held,
  /// none), or a swap (an own row there free at the follow, or held);
  /// whether the last row is on the track the iteration left, the refit's
  /// outcome, and whether the replay agrees with the iteration; then that
  /// track: good or not, its rows, its layers and those of its rows carrying
  /// the road's particle (hex) -- the track without its wrong rows; then the
  /// track edited, on both its lines: (a) a fake whose
  /// rows not carrying the road's particle are exactly the wrong last rows
  /// its agreeing walks left on it, without them; (b) a good track with the
  /// particle's own rows the unconfirmed drop removed from its agreeing walks put back where
  /// they were dropped (after the outward end, before the inward one, when
  /// that end is still the walk's last row or, the walk keeping none, the
  /// road's); "-" for neither. The edit's rows, its refit as `follow::Refit`
  /// makes it (`fitMeasurements`, the road's chi2 cut: ok, failedFit,
  /// aboveChi2), whether it then stays (`keptAfterFollow`: the length after
  /// the follow and the pT floor for its count; -1 without a refit) and its
  /// layers (hex). Between them: what the ending step saw --
  /// windows opened, rows in them, rows the gate refused -- and the evidence
  /// for the last row (when it ends the track that way) and the unconfirmed
  /// drop's row, truth-free, from the state before the row: the nominal chi2 and
  /// sqrt|V| (`oneStep`; a kink step's own chi2 is the scaled one), the
  /// density of the window's other rows beyond the gate and the layer's, the
  /// reach of its layer (`reachChance`), its time box's overlap with the
  /// walk's box before it, both widths [ps]; the chi2 predicted for it from the
  /// rest of the track (a fit without it; the dropped row: the track's fit)
  /// and, inward, the beam line's chi2 at the inner state with and without
  /// it; -1 where missing. Parallel over the tracks, written in their order.
  template <typename Observation, typename Material>
  void endDecisions(const Observation& iteration, const std::vector<Replay>& replayed, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    const auto view = iteration.tracks->view();
    const auto gate = static_cast<double>(parameters.clusterChi2Max);
    constexpr std::array<const char*, 4> kEnds{"noTarget", "holes", "noRow", "capacity"};
    constexpr std::array<const char*, 3> kRefits{"kept", "failedFit", "aboveChi2"};
    layerOccupancy(frame, geometry);
    std::vector<std::string> text(replayed.size());
    std::vector<std::string> edits(mEditCsv ? replayed.size() : 0U);
    o2::gtrk::detail::parallelFor(*mArena, 0U, replayed.size(), 256, [&](std::size_t a, std::size_t b) {
      for (auto element = static_cast<std::uint32_t>(a); element < b; ++element) {
        const auto track = iteration.trackOffset + element;
        const auto& replay = replayed[element];
        const std::span<const std::uint32_t> rows{view.clusters[track].index, view.count[track]};
        const auto particle = replay.roadCount == 0U ? study::kNoTrackParticle : study::goodParticle(mTruth->rows, std::span<const std::uint32_t>{replay.road.index, replay.roadCount});
        const auto found = particle == study::kNoTrackParticle ? mTruth->particles.end() : mTruth->particles.find(particle);
        const double pt = truePt(found);
        const bool trackGood = study::goodParticle(mTruth->rows, rows) != study::kNoTrackParticle;
        std::uint64_t trackLayers = 0U;
        std::uint64_t ownLayers = 0U;
        for (const auto row : rows) {
          trackLayers |= std::uint64_t{1} << mTruth->rows[row].layer;
          ownLayers |= particle != study::kNoTrackParticle && mTruth->rows[row].carries(particle) ? std::uint64_t{1} << mTruth->rows[row].layer : 0U;
        }
        const auto vertex = vertexRZ(found);
        const auto distance = [&](std::uint32_t one) { return pathDistance(one, vertex); };
        // What \a row is to the road's particle, and for a row past its end
        // its own rows beyond.
        const auto kindOf = [&](std::uint32_t row, bool outward) -> std::pair<const char*, const char*> {
          if (particle == study::kNoTrackParticle) {
            return {"-", "-"};
          }
          const auto& at = mTruth->rows[row];
          if (at.carries(particle)) {
            return {"own", "-"};
          }
          if (at.count == 0U) {
            return {"noise", "-"};
          }
          const auto& ownRows = mRowsOf.at(particle);
          bool onLayer = false;
          bool freeOnLayer = false;
          for (const auto own : ownRows) {
            if (mTruth->rows[own].layer == at.layer) {
              onLayer = true;
              freeOnLayer = freeOnLayer || (own < mFollowUsedNow.size() && mFollowUsedNow[own] == kFree);
            }
          }
          if (onLayer) {
            return {freeOnLayer ? "swapFree" : "swapHeld", "-"};
          }
          const double here = distance(row);
          bool beyond = false;
          bool free = false;
          for (const auto own : ownRows) {
            if (outward ? distance(own) > here : distance(own) < here) {
              beyond = true;
              free = free || (own < mFollowUsedNow.size() && mFollowUsedNow[own] == kFree);
            }
          }
          return {"past", !beyond ? "none" : (free ? "free" : "held")};
        };
        // The chi2 of \a row at the step that took it, and whether a kink step did.
        const auto chi2Of = [](const Walk& walk, std::uint32_t row) -> std::pair<double, int> {
          const auto step = std::ranges::find(walk.steps, row, &Step::choice);
          if (step == walk.steps.end()) {
            return {-1., -1};
          }
          const auto& seen = step->kink ? step->scaled : step->candidates;
          const auto at = std::ranges::find(seen, row, &Candidate::row);
          return {at == seen.end() ? -1. : static_cast<double>(at->chi2), step->kink ? 1 : 0};
        };
        // The fits of the track's rows: all of them, without the
        // outermost, without the innermost, with a row before the innermost.
        const auto fitOf = [&](std::vector<std::uint32_t> list) {
          if (list.empty()) {
            throw std::logic_error("CompleteAnatomy: a fit of no rows");
          }
          return fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = list.data()}, static_cast<int>(list.size()), kinkBudget(parameters));
        };
        // The beam line's chi2 at the inner state of \a fit, whose innermost row is \a innermost; -1 without.
        const auto beamOf = [&](const auto& fit, std::uint32_t innermost) {
          if (fit.status != PropagationStatus::kOk) {
            return -1.;
          }
          BoundTrackState state = fit.state;
          const auto& plane = geometry.planeOf(state);
          const LayerID layer = geometry.surface(state.surface).layer;
          const auto update = updateWithBeamLine(propagator, parameters, plane, layer, pairPointAt(frame.measurements(), innermost), state);
          return update.ok ? static_cast<double>(update.chi2) : -1.;
        };
        std::optional<decltype(fitOf({}))> fitAll;
        const auto all = [&]() -> const auto& {
          if (!fitAll) {
            fitAll = fitOf(std::vector<std::uint32_t>(rows.begin(), rows.end()));
          }
          return *fitAll;
        };
        const auto& proposal = iteration.proposals[element];
        // The evidence for \a row (the walk's last, or the unconfirmed drop's), from
        // the state before it: the step that took it, its window's density
        // of other rows beyond the gate and the layer's, the reach of its
        // layer, its time box against the walk's before it; then the chi2
        // predicted from the rest of the track (\a rest: a fit without it;
        // \a restState its end state toward the row) and inward the beam
        // line's chi2 with and without it.
        struct Evidence {
          double chi2 = -1., sqrtDet = -1., rhoWindow = -1., rhoLayer = -1., reach = -1., overlapPS = -1., boxPS = -1., rowPS = -1., rest = -1., beamWith = -1., beamWithout = -1., floatDet = -1., factorDet = -1.;
        };
        const auto evidenceOf = [&](const Walk& walk, std::uint32_t row, bool outward) {
          Evidence out;
          const auto direction = outward ? PropagationDirection::kForward : PropagationDirection::kBackward;
          const auto step = std::ranges::find(walk.steps, row, &Step::choice);
          if (step == walk.steps.end()) {
            return out;
          }
          const LayerID layer = mTruth->rows[row].layer;
          const auto taken = oneStep(row, step->state, step->factors, direction, frame, geometry, material, propagator);
          out.chi2 = taken.valid ? taken.chi2 : -1.;
          out.sqrtDet = taken.valid ? taken.sqrtDet : -1.;
          out.floatDet = taken.valid ? taken.floatDet : -1.;
          out.factorDet = taken.valid ? taken.factorDet : -1.;
          const auto target = std::ranges::find(step->targets, layer, &TargetSeen::layer);
          if (target != step->targets.end()) {
            const auto& window = target->window;
            out.rhoWindow = windowRho(parameters, *step, layer, window, gate, taken.valid ? taken.sqrtDet : 0., direction, frame, geometry, material, propagator, [row](std::uint32_t candidate) { return candidate == row; }, -1.);
            out.rhoLayer = layerDensity(layer, window);
            out.reach = reachChance(parameters, step->state, step->factors, layer, window, direction, geometry, propagator, row);
          }
          // The time box before the row: the road's, the outward rows kept
          // (inward), the walk's rows before it.
          auto box = trackTimeBox(frame.measurements(), replay.road, replay.roadCount);
          if (!outward) {
            for (std::uint8_t i = 0; i < proposal.outwardCount; ++i) {
              box += timeBoxAt(frame.measurements(), proposal.outward[i]);
            }
          }
          for (const auto before : walk.raw) {
            if (before == row) {
              break;
            }
            box += timeBoxAt(frame.measurements(), before);
          }
          const auto mine = timeBoxAt(frame.measurements(), row);
          out.boxPS = static_cast<double>(box.width);
          out.rowPS = static_cast<double>(mine.width);
          out.overlapPS = static_cast<double>(std::min(box.upper(), mine.upper())) - static_cast<double>(std::max(box.lower(), mine.lower()));
          return out;
        };
        std::vector<std::string> lines;
        std::uint32_t wrongLast = 0U;
        std::uint32_t appendOut = kNoRow;
        std::uint32_t prependIn = kNoRow;
        bool allAgree = true;
        // The truth-free edits each walk allows: its last row refused
        // where it ends the track, the unconfirmed drop's row put back where it attaches.
        std::array<bool, 2> refuse{}; // inward, outward
        std::array<std::uint32_t, 2> restore{kNoRow, kNoRow};
        const auto one = [&](const Walk& walk, bool outward) {
          bool crossed = false;
          std::uint32_t links = 0U;
          std::size_t windows = 0U;
          std::size_t seen = 0U;
          std::size_t refused = 0U;
          if (!walk.steps.empty()) {
            const auto& ending = walk.steps.back();
            for (const auto& [layer, called] : ending.holes) {
              crossed = crossed || called;
            }
            const auto current = ending.layer;
            const auto* offsets = outward ? parameters.allOutLinkOffsets : parameters.allInLinkOffsets;
            links = offsets[current + 1U] - offsets[current];
            windows = ending.targets.size();
            seen = ending.candidates.size();
            for (const auto& candidate : ending.candidates) {
              refused += candidate.outcome == StepOutcome::kChi2 ? 1U : 0U;
            }
          }
          const bool dropped = walk.dropped != kNoRow;
          const char* verdict = dropped ? "dropped" : (walk.raw.empty() ? "none" : (crossed ? "evidence" : "unconfirmed"));
          const std::uint32_t last = walk.raw.empty() ? kNoRow : walk.raw.back();
          const auto [lastChi2, lastKink] = last == kNoRow ? std::pair{-1., -1} : chi2Of(walk, last);
          const auto [lastKind, lastBeyond] = last == kNoRow ? std::pair{"-", "-"} : kindOf(last, outward);
          const double droppedChi2 = dropped ? chi2Of(walk, walk.dropped).first : -1.;
          const auto [droppedKind, droppedBeyond] = dropped ? kindOf(walk.dropped, outward) : std::pair{"-", "-"};
          const bool survives = last != kNoRow && std::ranges::find(rows, last) != rows.end();
          wrongLast += survives && particle != study::kNoTrackParticle && !mTruth->rows[last].carries(particle) ? 1U : 0U;
          // Where a dropped row would go: past the walk's last row, or the road's end where the walk kept none.
          bool attachable = false;
          if (dropped && !rows.empty() && replay.roadCount != 0U) {
            const std::uint32_t end = last != kNoRow ? last : (outward ? replay.road.index[replay.roadCount - 1U] : replay.road.index[0]);
            attachable = (outward ? rows.back() : rows.front()) == end;
          }
          if (attachable && walk.agrees && particle != study::kNoTrackParticle && mTruth->rows[walk.dropped].carries(particle)) {
            (outward ? appendOut : prependIn) = walk.dropped;
          }
          refuse[outward ? 1 : 0] = walk.agrees && survives && rows.size() >= 4U && (outward ? rows.back() : rows.front()) == last;
          restore[outward ? 1 : 0] = walk.agrees && attachable && rows.size() < kMaxClusters ? walk.dropped : kNoRow;
          // The evidence, for every agreeing walk; the rules apply to all tracks.
          Evidence lastEvidence;
          Evidence droppedEvidence;
          if (walk.agrees) {
            const auto direction = outward ? PropagationDirection::kForward : PropagationDirection::kBackward;
            if (survives && rows.size() >= 4U && (outward ? rows.back() : rows.front()) == last) {
              lastEvidence = evidenceOf(walk, last, outward);
              std::vector<std::uint32_t> rest(rows.begin(), rows.end());
              outward ? rest.pop_back() : static_cast<void>(rest.erase(rest.begin()));
              const auto fit = fitOf(rest);
              if (fit.status == PropagationStatus::kOk) {
                const auto& state = outward ? fit.outwardState : fit.state;
                const auto at = oneStep(last, state, factorCovariance(state.cov), direction, frame, geometry, material, propagator);
                lastEvidence.rest = at.valid ? at.chi2 : -1.;
                if (!outward) {
                  lastEvidence.beamWithout = beamOf(fit, rest.front());
                  lastEvidence.beamWith = beamOf(all(), rows.front());
                }
              }
            }
            if (dropped) {
              droppedEvidence = evidenceOf(walk, walk.dropped, outward);
              if (attachable) {
                const auto& fit = all();
                if (fit.status == PropagationStatus::kOk) {
                  const auto& state = outward ? fit.outwardState : fit.state;
                  const auto at = oneStep(walk.dropped, state, factorCovariance(state.cov), direction, frame, geometry, material, propagator);
                  droppedEvidence.rest = at.valid ? at.chi2 : -1.;
                  if (!outward) {
                    std::vector<std::uint32_t> with{walk.dropped};
                    with.insert(with.end(), rows.begin(), rows.end());
                    droppedEvidence.beamWithout = beamOf(fit, rows.front());
                    droppedEvidence.beamWith = beamOf(fitOf(with), walk.dropped);
                  }
                }
              }
            }
          }
          const auto evidence = [](const Evidence& e) { return std::format("{:.5g},{:.5g},{:.5g},{:.5g},{:.4g},{:.0f},{:.0f},{:.0f},{:.5g},{:.5g},{:.5g},{:.9g},{:.9g}", e.chi2, e.sqrtDet, e.rhoWindow, e.rhoLayer, e.reach, e.overlapPS, e.boxPS, e.rowPS, e.rest, e.beamWith, e.beamWithout, e.floatDet, e.factorDet); };
          // The beam line inward: its update's chi2 at the state after
          // the last row walked, as `confirmation` makes it for a kink row at
          // the innermost layer, and whether it passes beamChi2Max; and the
          // beam line as the reference of each nominal step (`oneStep` with
          // it): the rows the walk took that the reference's prediction puts
          // above the gate, the first of them -- its position among the rows
          // taken, whether it is the walk's last kept row (1) or the unconfirmed
          // drop's (2), the step's best other free row the reference passes (the
          // first target with one, the smallest chi2) -- and, where it refuses
          // none, a row of the ending step the nominal gate refused that the
          // reference passes; the reference's chi2 of the walk's last kept row
          // on its own, whatever came before (-1: a kink step's, -2: the
          // transport or the beam update failed); the rows taken it could not
          // judge so. Kink steps are not judged; a row it cannot judge is not
          // refused.
          double endBeam = -1.;
          int endBeamPasses = -1;
          int refRefused = -1;
          int refInvalid = -1;
          double lastRef = -1.;
          int refStep = -1;
          int refLast = -1;
          std::uint32_t refRow = kNoRow;
          std::uint32_t refAlt = kNoRow;
          std::uint32_t refNew = kNoRow;
          if (!outward && walk.agrees) {
            const auto at = last == kNoRow ? walk.steps.end() : std::ranges::find(walk.steps, last, &Step::choice);
            if (at != walk.steps.end() && std::next(at) != walk.steps.end()) {
              BoundTrackState after = std::next(at)->state;
              covarianceFromFactors(std::next(at)->factors, after.cov);
              const auto update = updateWithBeamLine(propagator, parameters, geometry.planeOf(after), geometry.surface(after.surface).layer, pairPointAt(frame.measurements(), last), after);
              if (update.ok) {
                endBeam = static_cast<double>(update.chi2);
                endBeamPasses = update.chi2 <= parameters.beamChi2Max ? 1 : 0;
              }
            }
            refRefused = 0;
            refInvalid = 0;
            int position = 0;
            for (std::size_t s = 0; s < walk.steps.size(); ++s) {
              const auto& step = walk.steps[s];
              const auto referenced = [&](std::uint32_t row) {
                const auto update = oneStep(row, step.state, step.factors, PropagationDirection::kBackward, frame, geometry, material, propagator, &parameters);
                return update.valid ? update.chi2 : -2.;
              };
              const auto bestOther = [&](std::uint32_t besides) {
                for (const auto& target : step.targets) {
                  std::uint32_t best = kNoRow;
                  double smallest = gate;
                  for (const auto& candidate : step.candidates) {
                    if (candidate.layer != target.layer || candidate.row == besides || (candidate.outcome != StepOutcome::kPasses && candidate.outcome != StepOutcome::kChi2)) {
                      continue;
                    }
                    const double chi2 = referenced(candidate.row);
                    if (chi2 >= 0. && chi2 <= smallest) {
                      best = candidate.row;
                      smallest = chi2;
                    }
                  }
                  if (best != kNoRow) {
                    return best;
                  }
                }
                return kNoRow;
              };
              if (step.choice == kNoRow) {
                if (s + 1U == walk.steps.size() && refStep < 0) {
                  refNew = bestOther(kNoRow);
                }
                continue;
              }
              const double chi2 = step.kink ? -1. : referenced(step.choice);
              lastRef = step.choice == last ? chi2 : lastRef;
              refInvalid += chi2 < -1.5 ? 1 : 0;
              if (chi2 > gate) {
                ++refRefused;
                if (refStep < 0) {
                  refStep = position;
                  refRow = step.choice;
                  refLast = step.choice == last ? 1 : (step.choice == walk.dropped ? 2 : 0);
                  refAlt = bestOther(step.choice);
                }
              }
              ++position;
            }
          }
          const auto rowOf = [](std::uint32_t row) { return row == kNoRow ? std::int64_t{-1} : static_cast<std::int64_t>(row); };
          const auto kindName = [&](std::uint32_t row) { return row == kNoRow ? "-" : kindOf(row, outward).first; };
          // The step that took the last row: each
          // layer it searched, in its order, with what it found -- r the
          // row, h a hole (a sensor crossed without a row), m a miss (no
          // sensor crossed: no hole) -- "/"-separated; "-" without one.
          std::string search = "-";
          if (walk.agrees && last != kNoRow) {
            const auto at = std::ranges::find(walk.steps, last, &Step::choice);
            if (at != walk.steps.end()) {
              search.clear();
              for (const auto& target : at->targets) {
                const auto hole = std::ranges::find(at->holes, target.layer, &std::pair<LayerID, bool>::first);
                const char code = target.layer == mTruth->rows[last].layer ? 'r' : (hole == at->holes.end() ? '?' : (hole->second ? 'h' : 'm'));
                search += std::format("{}{}{}", search.empty() ? "" : "/", static_cast<int>(target.layer), code);
              }
            }
          }
          const auto beam = std::format("{:.5g},{},{},{},{},{},{},{},{},{},{},{:.5g},{},{},{}", endBeam, endBeamPasses, refRefused, refStep, rowOf(refRow), kindName(refRow), refLast, rowOf(refAlt), kindName(refAlt), rowOf(refNew), kindName(refNew), lastRef, lastRef >= 0. ? (lastRef <= gate ? 1 : 0) : -1, refInvalid, search);
          lines.push_back(std::format("{},{},{},{},{:.6g},{},{},{},{},{},{},{},{},{},{:.5g},{},{},{},{},{},{:.5g},{},{},{},{},{},{},{},{:x},{:x},{},{},{},{},{}", track, iteration.iteration, outward ? "out" : "in", particle == study::kNoTrackParticle ? 0U : particle, pt, kEnds[static_cast<std::size_t>(walk.end)], crossed ? 1 : 0, links, verdict, walk.raw.size(), walk.bidKept, walk.kept, last == kNoRow ? -1 : static_cast<std::int64_t>(last), last == kNoRow ? -1 : static_cast<int>(mTruth->rows[last].layer), lastChi2, lastKink, lastKind, lastBeyond, dropped ? static_cast<std::int64_t>(walk.dropped) : -1, dropped ? static_cast<int>(mTruth->rows[walk.dropped].layer) : -1, droppedChi2, droppedKind, droppedBeyond, survives ? 1 : 0, kRefits[replay.refit], walk.agrees ? 1 : 0, trackGood ? 1 : 0, rows.size(), trackLayers, ownLayers, windows, seen, refused, evidence(lastEvidence), evidence(droppedEvidence)) + "," + beam);
          allAgree = allAgree && walk.agrees;
        };
        if (!replay.outward.empty()) {
          one(replay.outward.front(), true);
        }
        one(replay, false);
        // The edit: (a) the fake without its wrong last rows, (b) the good track with its own dropped rows back.
        std::uint32_t wrongRows = 0U;
        for (const auto row : rows) {
          wrongRows += particle != study::kNoTrackParticle && !mTruth->rows[row].carries(particle) ? 1U : 0U;
        }
        std::vector<std::uint32_t> edited;
        const char* edit = "-";
        if (particle != study::kNoTrackParticle && allAgree && !trackGood && wrongRows != 0U && wrongRows == wrongLast && study::majorityParticle(mTruth->rows, rows).particle == particle) {
          edit = "a";
          for (const auto row : rows) {
            if (mTruth->rows[row].carries(particle)) {
              edited.push_back(row);
            }
          }
        } else if (trackGood && (appendOut != kNoRow || prependIn != kNoRow) && rows.size() + 2U <= kMaxClusters) {
          edit = "b";
          if (prependIn != kNoRow) {
            edited.push_back(prependIn);
          }
          edited.insert(edited.end(), rows.begin(), rows.end());
          if (appendOut != kNoRow) {
            edited.push_back(appendOut);
          }
        }
        std::string tail = ",-,-1,-,-1,0";
        if (*edit != '-') {
          const auto fit = fitOf(edited);
          const bool ok = fit.status == PropagationStatus::kOk;
          const bool under = follow::passesRoadChi2(parameters, fit);
          std::uint64_t layers = 0U;
          for (const auto row : edited) {
            layers |= std::uint64_t{1} << mTruth->rows[row].layer;
          }
          tail = std::format(",{},{},{},{},{:x}", edit, edited.size(), !ok ? "failedFit" : (under ? "ok" : "aboveChi2"), under ? (keptAfterFollow(parameters, static_cast<std::uint32_t>(edited.size()), fit.state) ? 1 : 0) : -1, layers);
        }
        for (const auto& line : lines) {
          text[element] += line;
          text[element] += tail;
          text[element] += '\n';
        }
        // The truth-free end edits: the track as it is, then every
        // combination of the edits its walks allow -- inward and outward
        // each none (0), the last row refused (R), the unconfirmed drop's row back (K)
        // -- refitted and floored as (a) and (b): the rows, the refit, kept
        // after the floors, the good particle (0: a fake), the majority's.
        if (mEditCsv) {
          const auto idOf = [](std::uint64_t id) { return id == study::kNoTrackParticle ? std::uint64_t{0} : id; };
          auto& out = edits[element];
          out += std::format("{},{},0,0,{},ok,1,{},{},{:x},", track, iteration.iteration, rows.size(), idOf(study::goodParticle(mTruth->rows, rows)), idOf(study::majorityParticle(mTruth->rows, rows).particle), trackLayers);
          for (std::size_t i = 0; i < rows.size(); ++i) {
            out += std::format("{}{}", i == 0U ? "" : ";", rows[i]);
          }
          out += '\n';
          for (const char in : {'0', 'R', 'K'}) {
            for (const char outer : {'0', 'R', 'K'}) {
              const auto allowed = [&](char change, std::size_t side) { return change == '0' || (change == 'R' ? refuse[side] : restore[side] != kNoRow); };
              if ((in == '0' && outer == '0') || !allowed(in, 0U) || !allowed(outer, 1U) || (in == 'K' && outer == 'K' && rows.size() + 2U > kMaxClusters)) {
                continue;
              }
              std::vector<std::uint32_t> list(rows.begin(), rows.end());
              if (outer == 'R') {
                list.pop_back();
              } else if (outer == 'K') {
                list.push_back(restore[1]);
              }
              if (in == 'R') {
                list.erase(list.begin());
              } else if (in == 'K') {
                list.insert(list.begin(), restore[0]);
              }
              const auto fit = fitOf(list);
              const bool ok = fit.status == PropagationStatus::kOk;
              const bool under = follow::passesRoadChi2(parameters, fit);
              std::uint64_t layers = 0U;
              for (const auto row : list) {
                layers |= std::uint64_t{1} << mTruth->rows[row].layer;
              }
              out += std::format("{},{},{},{},{},{},{},{},{},{:x},\n", track, iteration.iteration, in, outer, list.size(), !ok ? "failedFit" : (under ? "ok" : "aboveChi2"), under && keptAfterFollow(parameters, static_cast<std::uint32_t>(list.size()), fit.state) ? 1 : 0, idOf(study::goodParticle(mTruth->rows, list)), idOf(study::majorityParticle(mTruth->rows, list).particle), layers);
            }
          }
        }
      }
    });
    for (const auto& one : text) {
      *mEndCsv << one;
    }
    for (const auto& one : edits) {
      *mEditCsv << one;
    }
  }

  /// The particle's own row at every step of the follow: on the agreeing walks of the correct roads of the class's
  /// particles, each step's target layers the particle has a row on, the
  /// best outcome of its rows there -- passing (taken or not) or refused by
  /// the chi2 gate, with the chi2; outside the window; held; no update; not
  /// visited -- by direction, link, pT range and the particle's true kink
  /// over the step (the crossing pair from the inner layer to the outer,
  /// over theta0). Without a kink (<= 3 theta0) the update's chi2 is a
  /// chi-square of two degrees of freedom where the error model holds: the
  /// gate's tail is exp(-gate / 2), 0.27 % at 11.83. Also the outward unconfirmed drops
  /// on every agreeing walk: whether the dropped row is the road particle's,
  /// and whether the last step saw a row near (a row in its window the gate
  /// refused, inward's evidence) and whose -- the gain and the cost of
  /// inward's evidence rule outward.
  template <typename Observation>
  void ownSteps(const Observation& iteration, const std::vector<Replay>& replayed, const IterationParameters<Det>& parameters)
  {
    enum Outcome : std::uint8_t { kPass,
                                  kRefused,
                                  kOutside,
                                  kHeld,
                                  kNoUpdateOwn,
                                  kUnvisited,
                                  kOutcomes };
    constexpr std::array<double, 10> kChi2Edges{2., 4., 6., 9.21, 0., 15., 25., 50., 100., std::numeric_limits<double>::infinity()};
    struct Cell {
      std::array<std::size_t, kOutcomes> outcome{};
      std::array<std::size_t, kChi2Edges.size()> chi2{};
    };
    // [direction][link from][link to][pT range][kink bin]
    std::map<std::tuple<std::size_t, int, int, std::size_t, std::size_t>, Cell> cells;
    const auto gate = static_cast<double>(parameters.clusterChi2Max);
    auto edges = kChi2Edges;
    edges[4] = gate;
    std::ranges::sort(edges);
    // Unconfirmed drops outward: [dropped row the particle's][near row: none, the particle's, another's]
    std::array<std::array<std::size_t, 3>, 2> unconfirmedDrops{};
    // and by the dropped row's chi2 at the step that took it: [own][near][bin] (< 1 / 1-2 / 2-4 / 4-6 / 6-9.21 / >= 9.21)
    constexpr std::array<double, 6> kDropEdges{1., 2., 4., 6., 9.21, std::numeric_limits<double>::infinity()};
    std::array<std::array<std::array<std::size_t, kDropEdges.size()>, 3>, 2> r1Chi2{};
    for (std::uint32_t element = 0; element < replayed.size(); ++element) {
      const auto& replay = replayed[element];
      if (replay.roadCount == 0U) {
        continue;
      }
      const auto particle = study::goodParticle(mTruth->rows, std::span<const std::uint32_t>{replay.road.index, replay.roadCount});
      if (!replay.outward.empty() && replay.outward.front().agrees && replay.outward.front().dropped != kNoRow && particle != study::kNoTrackParticle) {
        const auto& walk = replay.outward.front();
        const bool ownDropped = mTruth->rows[walk.dropped].carries(particle);
        std::size_t near = 0U;
        if (!walk.steps.empty()) {
          for (const auto& candidate : walk.steps.back().candidates) {
            if (candidate.outcome == StepOutcome::kChi2) {
              near = std::max<std::size_t>(near, mTruth->rows[candidate.row].carries(particle) ? 1U : 2U);
            }
          }
        }
        // the particle's near row wins over another's
        if (near == 2U) {
          for (const auto& candidate : walk.steps.back().candidates) {
            near = candidate.outcome == StepOutcome::kChi2 && mTruth->rows[candidate.row].carries(particle) ? 1U : near;
          }
        }
        ++unconfirmedDrops[ownDropped ? 0U : 1U][near];
        for (const auto& step : walk.steps) {
          if (step.choice != walk.dropped) {
            continue;
          }
          const auto taken = std::ranges::find(step.candidates, walk.dropped, &Candidate::row);
          if (taken != step.candidates.end()) {
            const auto bin = static_cast<std::size_t>(std::ranges::upper_bound(kDropEdges, static_cast<double>(taken->chi2)) - kDropEdges.begin());
            ++r1Chi2[ownDropped ? 0U : 1U][near][std::min(bin, kDropEdges.size() - 1U)];
          }
        }
      }
      if (particle == study::kNoTrackParticle || !mIndex.contains(particle)) {
        continue;
      }
      const auto truth = mTruth->particles.find(particle);
      const double pt = truePt(truth);
      const std::size_t ptRange = pt < 0.6 ? 0U : (pt < 1. ? 1U : 2U);
      const auto tally = [&](const Walk& walk, std::size_t direction) {
        if (!walk.agrees) {
          return;
        }
        for (const auto& step : walk.steps) {
          for (const auto& target : step.targets) {
            const auto layer = target.layer;
            if (ownRowOn(particle, layer) == kNoTrack) {
              continue;
            }
            Outcome best = kUnvisited;
            double chi2 = std::numeric_limits<double>::infinity();
            for (const auto& candidate : step.candidates) {
              if (candidate.layer != layer || !mTruth->rows[candidate.row].carries(particle)) {
                continue;
              }
              const auto kind = candidate.outcome == StepOutcome::kPasses ? kPass : candidate.outcome == StepOutcome::kChi2          ? kRefused
                                                                                  : candidate.outcome == StepOutcome::kOutsideWindow ? kOutside
                                                                                  : candidate.outcome == StepOutcome::kUsed          ? kHeld
                                                                                                                                     : kNoUpdateOwn;
              if (kind < best || (kind == best && static_cast<double>(candidate.chi2) < chi2)) {
                best = kind;
                chi2 = kind <= kRefused ? static_cast<double>(candidate.chi2) : chi2;
              }
            }
            const auto outer = direction == 0U ? layer : step.layer;
            const auto bin = kinkBin(kinkAt(particle, outer, true));
            auto& cell = cells[{direction, static_cast<int>(step.layer), static_cast<int>(layer), ptRange, bin}];
            ++cell.outcome[best];
            if (best <= kRefused) {
              const auto at = static_cast<std::size_t>(std::ranges::upper_bound(edges, chi2) - edges.begin());
              ++cell.chi2[std::min(at, edges.size() - 1U)];
            }
          }
        }
      };
      if (!replay.outward.empty()) {
        tally(replay.outward.front(), 0U);
      }
      tally(replay, 1U);
    }
    // The chi2 bins' labels from the sorted edges, the gate's by name: where
    // it falls among the fixed edges is printed, not assumed.
    const auto edgeName = [&](double edge) { return edge == gate ? std::string{"gate"} : std::format("{}", edge); };
    std::string bins = "<" + edgeName(edges[0]);
    for (std::size_t k = 1; k < edges.size(); ++k) {
      bins += std::isinf(edges[k]) ? "/>" + edgeName(edges[k - 1]) : "/" + edgeName(edges[k - 1]) + "-" + edgeName(edges[k]);
    }
    std::string text = std::format("  own steps after iteration {} (gate {:.2f}): direction link pT kink: pass/refused/outside/held/noupdate/unvisited | chi2 {}", iteration.iteration, gate, bins);
    for (const auto& [key, cell] : cells) {
      const auto& [direction, from, to, ptRange, bin] = key;
      text += std::format("\n    O {} {}-{} {} {}: {}/{}/{}/{}/{}/{} |", direction, from, to, ptRange, bin, cell.outcome[0], cell.outcome[1], cell.outcome[2], cell.outcome[3], cell.outcome[4], cell.outcome[5]);
      for (const auto n : cell.chi2) {
        text += std::format(" {}", n);
      }
    }
    text += std::format("\n    outward unconfirmed drops (dropped row the particle's / another's; near row none / the particle's / another's): {}/{}/{} {}/{}/{}", unconfirmedDrops[0][0], unconfirmedDrops[0][1], unconfirmedDrops[0][2], unconfirmedDrops[1][0], unconfirmedDrops[1][1], unconfirmedDrops[1][2]);
    for (std::size_t own = 0; own < 2U; ++own) {
      for (std::size_t near = 0; near < 3U; ++near) {
        const auto& c = r1Chi2[own][near];
        text += std::format("\n    outward unconfirmed drops by the dropped row's chi2 (<1/1-2/2-4/4-6/6-9.21/>=9.21), dropped {} near {}: {}/{}/{}/{}/{}/{}", own == 0U ? "own" : "wrong", near == 0U ? "none" : (near == 1U ? "own" : "other"), c[0], c[1], c[2], c[3], c[4], c[5]);
      }
    }
    mTransitionLines.push_back(text);
  }

  /// \a particle's MC hit segments (`--hits`: the signal source's only);
  /// none (nullptr) for another source or a particle without hits.
  [[nodiscard]] const std::vector<const study::HitSegment*>* segmentsOf(std::uint64_t particle) const
  {
    const auto label = workflow::labelOfParticle(particle);
    if (label.getSourceID() != 0) {
      return nullptr;
    }
    const auto found = mSegmentsOf.find(segmentKey(label.getEventID(), label.getTrackID()));
    return found != mSegmentsOf.end() ? &found->second : nullptr;
  }

  /// The follow's hole calls against the MC hits: on
  /// the agreeing walks of the iteration's correct roads (every road row the
  /// road's particle's), each nominal target without a row, by direction and
  /// layer: whether `crossesSensor` called it crossed, whether the particle
  /// has an MC hit on that layer, and whether it has a row there. Crossed
  /// without a hit is a gap called a crossing; a hit and no row the
  /// detector's own miss; a row the particle's row refused. Also the calls
  /// that ended a walk (its end `holes`, the last step's last call), with
  /// whether the particle has a row further along than the walk's last --
  /// by the distance from its production vertex, outward farther, inward
  /// nearer: what the end loses.
  template <typename Observation>
  void holeCalls(const Observation& iteration, const std::vector<Replay>& replayed, const IterationParameters<Det>& parameters, const TrackingGeometry<Det>& geometry, const Propagator& propagator)
  {
    // The crossed ending calls' `reachChance` (the step's state drawn from
    // its covariance, as crossesSensor carried): [direction][hit][further][bin]
    constexpr std::array<double, 5> kChanceEdges{0.1, 0.5, 0.9, 0.99, 2.};
    std::array<std::array<std::array<std::array<std::size_t, kChanceEdges.size()>, 2>, 2>, 2> chances{};
    // [direction][ended the walk][layer][crossed][hit][row]
    using Cube = std::array<std::array<std::array<std::size_t, 2>, 2>, 2>;
    std::array<std::array<std::array<Cube, kLayers>, 2>, 2> calls{};
    // The ending calls: [direction][layer][crossed][hit][a row further along]
    std::array<std::array<Cube, kLayers>, 2> further{};
    for (std::uint32_t element = 0; element < replayed.size(); ++element) {
      const auto& replay = replayed[element];
      if (replay.roadCount == 0U) {
        continue;
      }
      const auto particle = study::goodParticle(mTruth->rows, std::span<const std::uint32_t>{replay.road.index, replay.roadCount});
      if (particle == study::kNoTrackParticle) {
        continue;
      }
      std::uint64_t hitLayers = 0U;
      const auto* segments = segmentsOf(particle);
      if (segments == nullptr) {
        continue;
      }
      for (const auto* segment : *segments) {
        hitLayers |= std::uint64_t{1} << geometry.surface(static_cast<SensorID>(segment->chip)).layer;
      }
      std::uint64_t rowLayers = 0U;
      if (const auto rows = mRowsOf.find(particle); rows != mRowsOf.end()) {
        for (const auto row : rows->second) {
          rowLayers |= std::uint64_t{1} << mTruth->rows[row].layer;
        }
      }
      const auto vertex = vertexRZ(mTruth->particles.find(particle));
      const auto distance = [&](std::uint32_t row) { return pathDistance(row, vertex); };
      const auto tally = [&](const Walk& walk, std::size_t direction) {
        if (!walk.agrees) {
          return;
        }
        // The walk's last row: its last step's taken row, else the road's end on that side.
        double last = direction == 0U ? 0. : std::numeric_limits<double>::infinity();
        for (std::uint8_t i = 0; i < replay.roadCount; ++i) {
          last = direction == 0U ? std::max(last, distance(replay.road.index[i])) : std::min(last, distance(replay.road.index[i]));
        }
        for (const auto& step : walk.steps) {
          if (step.choice != kNoRow) {
            last = distance(step.choice);
          }
        }
        bool beyond = false;
        if (const auto rows = mRowsOf.find(particle); rows != mRowsOf.end()) {
          for (const auto row : rows->second) {
            beyond = beyond || (direction == 0U ? distance(row) > last + 0.5 : distance(row) < last - 0.5);
          }
        }
        for (std::size_t s = 0; s < walk.steps.size(); ++s) {
          const auto& holes = walk.steps[s].holes;
          for (std::size_t h = 0; h < holes.size(); ++h) {
            const auto [layer, crossed] = holes[h];
            const bool ended = walk.end == FollowEnd::kHoles && s + 1U == walk.steps.size() && h + 1U == holes.size();
            ++calls[direction][ended ? 1U : 0U][layer][crossed ? 1U : 0U][(hitLayers >> layer) & 1U][(rowLayers >> layer) & 1U];
            if (ended) {
              ++further[direction][layer][crossed ? 1U : 0U][(hitLayers >> layer) & 1U][beyond ? 1U : 0U];
              const auto& step = walk.steps[s];
              const auto target = std::ranges::find(step.targets, layer, &TargetSeen::layer);
              if (crossed && target != step.targets.end()) {
                const double chance = reachChance(parameters, step.state, step.factors, layer, target->window, direction == 0U ? PropagationDirection::kForward : PropagationDirection::kBackward, geometry, propagator, element, 32);
                const auto bin = static_cast<std::size_t>(std::ranges::upper_bound(kChanceEdges, chance) - kChanceEdges.begin());
                ++chances[direction][(hitLayers >> layer) & 1U][beyond ? 1U : 0U][std::min(bin, kChanceEdges.size() - 1U)];
              }
            }
          }
        }
      };
      if (!replay.outward.empty()) {
        tally(replay.outward.front(), 0U);
      }
      tally(replay, 1U);
    }
    std::string text = std::format(
      "  hole calls after iteration {} on the correct roads' walks (crossed / not, by the particle's MC hit and row on "
      "the layer: hit+row hit-row nohit+row nohit-row):",
      iteration.iteration);
    for (std::size_t direction = 0; direction < 2U; ++direction) {
      for (std::size_t ended = 0; ended < 2U; ++ended) {
        text += std::format("\n    {}{}:", direction == 0U ? "outward" : "inward", ended == 1U ? ", the call that ended the walk" : "");
        for (std::size_t layer = 0; layer < kLayers; ++layer) {
          const auto& c = calls[direction][ended][layer];
          std::size_t all = 0U;
          for (const auto& a : c) {
            for (const auto& b : a) {
              all += b[0] + b[1];
            }
          }
          if (all == 0U) {
            continue;
          }
          text += std::format(" L{} crossed {}/{}/{}/{} not {}/{}/{}/{};", layer, c[1][1][1], c[1][1][0], c[1][0][1], c[1][0][0], c[0][1][1], c[0][1][0], c[0][0][1], c[0][0][0]);
        }
      }
      text += std::format("\n    {}, the ending calls by a row further along (crossed: hit+further hit-further nohit+further nohit-further):", direction == 0U ? "outward" : "inward");
      for (std::size_t layer = 0; layer < kLayers; ++layer) {
        const auto& c = further[direction][layer];
        if (c[1][1][1] + c[1][1][0] + c[1][0][1] + c[1][0][0] != 0U) {
          text += std::format(" F{} {}/{}/{}/{};", layer, c[1][1][1], c[1][1][0], c[1][0][1], c[1][0][0]);
        }
      }
      text += std::format("\n    {}, the crossed ending calls' reach chance (< 0.1 / 0.1-0.5 / 0.5-0.9 / 0.9-0.99 / > 0.99):", direction == 0U ? "outward" : "inward");
      constexpr std::array<const char*, 4> kKinds{"hit, further", "hit, none", "no hit, further", "no hit, none"};
      for (std::size_t hit = 0; hit < 2U; ++hit) {
        for (std::size_t beyond = 0; beyond < 2U; ++beyond) {
          const auto& c = chances[direction][1U - hit][1U - beyond];
          text += std::format(" {}: {}/{}/{}/{}/{};", kKinds[(2U * hit) + beyond], c[0], c[1], c[2], c[3], c[4]);
        }
      }
    }
    mTransitionLines.push_back(text);
  }

  /// Which of the kink allowance's tests each of the iteration's own roads
  /// needs, and how its outermost row fits the rest (the road-end fakes): the
  /// roads of fakes with one wrong row, the
  /// outermost, against the good roads, those of a particle with a true kink
  /// above 3 theta0 just inside that row apart. The allowance's tests, each
  /// against its nominal twin: a kinked triplet at the road's outer end, one
  /// elsewhere, the road fit's last-step kink (the nominal fit refuses, the
  /// one with the last step allowed accepts), the beam gate's tail. The
  /// outermost row's chi2: the road less that row fitted (nominal), carried
  /// to the row and updated with it.
  template <typename Observation, typename Material>
  void roadEndAudit(const Observation& iteration, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    const auto view = iteration.tracks->view();
    const auto tracks = iteration.tracks->size();
    const auto columns = frame.measurements();
    auto nominal = std::make_unique<IterationParameters<Det>>(parameters);
    nominal->kinkScale2 = 1.F;
    const KinkBudget nominalBudget{.max = static_cast<int>(parameters.maxKinks), .gate = parameters.kinkChi2Min};
    enum Population : std::uint8_t { kEndFake,
                                     kOtherFake,
                                     kGoodEndKink,
                                     kGoodQuiet,
                                     kGoodInnerKink,
                                     kGoodOther,
                                     kPopulations };
    constexpr std::array<const char*, kPopulations> kPopulationNames{"fakes, one wrong row, the outermost", "other road fakes", "good, a true kink > 3 theta0 inside the outermost row", "good, no true kink > 3 theta0 in the road", "good, a true kink > 3 theta0 inside the road, not at its end", "good, not the class's or without MC hits"};
    // The allowance's tests a road needs, a bit each.
    enum Need : std::uint8_t { kOuterTriplet = 1U,
                               kInnerTriplet = 2U,
                               kLastStep = 4U,
                               kBeamTail = 8U };
    constexpr std::array<double, 8> kChi2Edges{0., 2., 4., 6., 9., 13.8, 25., 100.};
    std::array<std::array<std::size_t, 16>, kPopulations> needs{};
    std::array<std::array<std::size_t, kChi2Edges.size() + 1U>, kPopulations> chi2s{};
    std::array<std::array<std::size_t, kChi2Edges.size() + 1U>, kPopulations> scaledChi2s{};
    std::array<std::size_t, kPopulations> scaledPass{};
    std::array<std::size_t, kPopulations> counts{};
    std::array<std::size_t, kPopulations> unfitted{};
    for (std::uint32_t track = iteration.trackOffset; track < tracks; ++track) {
      const auto element = track - iteration.trackOffset;
      const auto* proposal = element < iteration.proposals.size() ? &iteration.proposals[element] : nullptr;
      auto road = roadOf(std::span<const std::uint32_t>{view.clusters[track].index, view.count[track]}, proposal);
      if (road.size() < 4U) {
        continue;
      }
      std::ranges::sort(road, [&](std::uint32_t a, std::uint32_t b) { return mTruth->rows[a].layer < mTruth->rows[b].layer; });
      const auto outer = road.back();
      const auto outerLayer = mTruth->rows[outer].layer;
      // The population.
      const auto good = study::goodParticle(mTruth->rows, road);
      Population population = kOtherFake;
      if (good == study::kNoTrackParticle) {
        const auto majority = study::majorityParticle(mTruth->rows, road);
        std::size_t wrong = 0U;
        for (const auto row : road) {
          wrong += mTruth->rows[row].carries(majority.particle) ? 0U : 1U;
        }
        population = wrong == 1U && !mTruth->rows[outer].carries(majority.particle) ? kEndFake : kOtherFake;
      } else {
        const auto at = mIndex.find(good);
        population = kGoodOther;
        if (at != mIndex.end() && mKinksDone && at->second + 1U < mKinkBegin.size()) {
          population = kGoodQuiet;
          for (auto i = mKinkBegin[at->second]; i < mKinkBegin[at->second + 1U]; ++i) {
            const auto& kink = mKinks[i];
            if (kink.layer < mTruth->rows[road.front()].layer || kink.to > outerLayer) {
              continue;
            }
            const float size = std::max(std::abs(kink.theta), std::abs(kink.phi));
            if (size > 3.F && kink.to == outerLayer) {
              population = kGoodEndKink;
              break;
            }
            if (size > 3.F) {
              population = kGoodInnerKink;
            }
          }
        }
      }
      ++counts[population];
      // The triplets, nominal or kinked.
      unsigned need = 0U;
      for (std::size_t i = 0; i + 2U < road.size(); ++i) {
        const auto a = timedPointAt(columns, road[i]);
        const auto b = timedPointAt(columns, road[i + 1U]);
        const auto c = timedPointAt(columns, road[i + 2U]);
        TripletFitResult fit;
        KinkTerms terms;
        if (tripletAudit(parameters, mTruth->rows[road[i + 1U]].layer, a, b, c, 1.F, fit, terms) != kGraphPasses) {
          need |= i + 3U == road.size() ? kOuterTriplet : kInnerTriplet;
        }
      }
      // The road fit, its last step nominal and allowed.
      const RowChain<Det> chain{.frame = &frame, .rows = road.data()};
      const auto count = static_cast<int>(road.size());
      auto plainParameters = parameters;
      plainParameters.kinkScale2 = 1.F;
      plainParameters.maxKinks = nominalBudget.max;
      plainParameters.kinkChi2Min = nominalBudget.gate;
      const auto plain = fitRoad(plainParameters, propagator, geometry, material, chain, count);
      const auto allowed = fitRoad(parameters, propagator, geometry, material, chain, count);
      if (plain.status != RoadFitStatus::kAccepted && allowed.status == RoadFitStatus::kAccepted) {
        need |= kLastStep;
      }
      // The beam gate, nominal and with its tail.
      const auto& fitted = plain.status == RoadFitStatus::kAccepted ? plain : allowed;
      Scalar beamChi2 = 0.F;
      if (parameters.beamConstraint && fitted.status == RoadFitStatus::kAccepted && !clearsBeamLine(propagator, *nominal, frame, geometry, material, road.data(), static_cast<std::uint32_t>(road.size()), fitted.state, beamChi2)) {
        need |= kBeamTail;
      }
      ++needs[population][need];
      // The outermost row against the rest.
      const auto inner = fitOutward(propagator, geometry, material, chain, count - 1, parameters.roadChi2MaxByNdf, nominalBudget);
      BoundTrackState output;
      CovarianceFactors outputFactors;
      Scalar chi2 = 0.F;
      if (inner.status != RoadFitStatus::kAccepted || !follow::updateCandidate<false>(propagator, geometry, material, measurementAt(frame, outer), PropagationDirection::kForward, inner.state, factorCovariance(inner.state.cov), output, outputFactors, chi2, 1.F)) {
        ++unfitted[population];
        continue;
      }
      const auto binOf = [&](Scalar value) {
        std::size_t bin = 0U;
        while (bin < kChi2Edges.size() && static_cast<double>(value) >= kChi2Edges[bin]) {
          ++bin;
        }
        return bin;
      };
      ++chi2s[population][binOf(chi2)];
      // The same step with its scattering scaled by the allowance, as the
      // follow's kink step takes it, against the follow's gate.
      Scalar scaled = 0.F;
      if (follow::updateCandidate<true>(propagator, geometry, material, measurementAt(frame, outer), PropagationDirection::kForward, inner.state, factorCovariance(inner.state.cov), output, outputFactors, scaled, std::sqrt(parameters.kinkScale2))) {
        ++scaledChi2s[population][binOf(scaled)];
        scaledPass[population] += scaled <= parameters.clusterChi2Max ? 1U : 0U;
      }
    }
    std::string text = std::format("  the road ends after iteration {} (its own roads; the allowance's tests each needs, against its nominal twin):", iteration.iteration);
    constexpr std::array<const char*, 4> kNeedNames{"outer triplet", "inner triplet", "last-step kink", "beam tail"};
    for (std::size_t population = 0; population < kPopulations; ++population) {
      if (counts[population] == 0U) {
        continue;
      }
      text += std::format("\n    {} {}: needs", kPopulationNames[population], counts[population]);
      for (unsigned need = 0; need < 16U; ++need) {
        if (needs[population][need] == 0U) {
          continue;
        }
        std::string name;
        for (std::size_t bit = 0; bit < kNeedNames.size(); ++bit) {
          if (((need >> bit) & 1U) != 0U) {
            name += name.empty() ? kNeedNames[bit] : std::string{" + "} + kNeedNames[bit];
          }
        }
        text += std::format(" {} {};", name.empty() ? "none" : name, needs[population][need]);
      }
      text += " | the outermost row's chi2 on the rest:";
      for (std::size_t bin = 0; bin <= kChi2Edges.size(); ++bin) {
        if (chi2s[population][bin] != 0U) {
          text += bin == kChi2Edges.size() ? std::format(" >= {} {};", kChi2Edges.back(), chi2s[population][bin]) : std::format(" < {} {};", kChi2Edges[bin], chi2s[population][bin]);
        }
      }
      text += std::format(" the rest not fitted {} | scaled by the allowance:", unfitted[population]);
      for (std::size_t bin = 0; bin <= kChi2Edges.size(); ++bin) {
        if (scaledChi2s[population][bin] != 0U) {
          text += bin == kChi2Edges.size() ? std::format(" >= {} {};", kChi2Edges.back(), scaledChi2s[population][bin]) : std::format(" < {} {};", kChi2Edges[bin], scaledChi2s[population][bin]);
        }
      }
      text += std::format(" within the step gate ({:.1f}) {}", parameters.clusterChi2Max, scaledPass[population]);
    }
    mTransitionLines.push_back(text);
  }

  /// The class's fates after iteration \a iteration, cumulative: complete,
  /// short, not found and found at all, over the class and by the pT ranges
  /// the iterations' floors bound (from each floor up to the next higher).
  void progress(std::size_t iteration)
  {
    std::vector<double> edges = mFloors;
    std::ranges::sort(edges);
    edges.erase(std::ranges::unique(edges).begin(), edges.end());
    const std::size_t ranges = edges.size() + 1U;           // the last: the whole class
    std::vector<std::array<std::size_t, 3>> counts(ranges); // by Fate
    for (std::size_t index = 0; index < mFindable.size(); ++index) {
      const auto found = mTruth->particles.find(mFindable[index]);
      const double pt = truePt(found);
      ++counts.back()[mFate[index]];
      for (std::size_t range = edges.size(); range-- > 0U;) {
        if (pt >= edges[range]) {
          ++counts[range][mFate[index]];
          break;
        }
      }
    }
    std::string text = std::format("  after iteration {} (cumulative; complete / short / not found / found, % of the class in each range):", iteration);
    for (std::size_t range = ranges; range-- > 0U;) {
      const auto& c = counts[range];
      const std::size_t n = c[kComplete] + c[kShortFate] + c[kNotFound];
      if (n == 0U) {
        continue;
      }
      const auto pc = [&](std::size_t k) { return study::percent(k, n); };
      const std::string name = range == ranges - 1U ? std::string{"all"} : (range + 1U < edges.size() ? std::format("{}-{} GeV", edges[range], edges[range + 1U]) : std::format(">= {} GeV", edges[range]));
      text += std::format("\n    {} ({}): {:.2f} / {:.2f} / {:.2f} / {:.2f}", name, n, pc(c[kComplete]), pc(c[kShortFate]), pc(c[kNotFound]), pc(c[kComplete] + c[kShortFate]));
    }
    mProgressLines.push_back(text);
  }

  /// The inward pass of a refit of \a rows from \a seed, as
  /// `detail::inwardPass` runs it: its status, the chain
  /// index and the step it failed at, the state the failing step started
  /// from, and with \a record each state: the seed, then each filtered one
  /// with its row and increment, the covariance written. It is the backward
  /// loop of `detail::runPass` (`TrackFitter.h`), recording; kept in step
  /// with it, a difference counted as a mismatch (`mInwardMismatch`).
  struct InwardStep {
    std::uint32_t row = 0U;
    BoundTrackState state;
    Scalar chi2 = 0.F;
  };
  struct InwardRun {
    PropagationStatus status = PropagationStatus::kOk;
    int failedIndex = -1;
    int failedStep = 0;
    BoundTrackState before;
    std::vector<InwardStep> history;
  };
  template <typename Material>
  [[nodiscard]] static InwardRun runInward(const BoundTrackState& seed, const std::vector<std::uint32_t>& rows, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, bool record)
  {
    InwardRun out;
    const int n = static_cast<int>(rows.size());
    const RowChain<Det> chain{.frame = &frame, .rows = rows.data()};
    const auto budget = o2::gtrk::detail::kinkBudgetFor(kinkBudget(parameters), n);
    BoundTrackState state = seed;
    CovarianceFactors factors;
    if (record) {
      out.history.push_back({.row = rows[static_cast<std::size_t>(n - 1)], .state = seed, .chi2 = 0.F});
    }
    const auto first = o2::gtrk::detail::startOnMeasurement(state, factors, chain[n - 1]);
    if (!updateUsable(first, factors)) {
      out.status = PropagationStatus::kNotFinite;
      out.failedIndex = n - 1;
      return out;
    }
    Scalar chi2 = first.chi2;
    int kinks = 0;
    for (int step = 1; step < n; ++step) {
      const int index = n - 1 - step;
      const auto& stepProp = stepPropagator(propagator, geometry, chain, index + 1, index);
      out.before = state;
      covarianceFromFactors(factors, out.before.cov);
      const Scalar before = chi2;
      const auto status = o2::gtrk::detail::kinkStep(stepProp, geometry, material, chain[index], PropagationDirection::kBackward, state, factors, chi2, budget, step < n - 1, kinks);
      if (status != PropagationStatus::kOk) {
        out.status = status;
        out.failedIndex = index;
        out.failedStep = step;
        return out;
      }
      if (record) {
        InwardStep one{.row = rows[static_cast<std::size_t>(index)], .state = state, .chi2 = chi2 - before};
        covarianceFromFactors(factors, one.state.cov);
        out.history.push_back(one);
      }
    }
    return out;
  }

  /// Both crossings of a barrel plane's line by \a state's transverse circle,
  /// from the geometry here rather than by the solve's
  /// quadratic: each one's 3D path in \a direction's sense (positive: along
  /// it), turn [rad], offset along the line from the start's foot, local
  /// coordinates, whether inside the sensor, incidence and distance from the
  /// plane's centre; which the solve tries first (the smaller offset).
  struct Root {
    double path = 0.;
    double turn = 0.;
    double mu = 0.;
    double loc0 = 0.;
    double loc1 = 0.;
    double cosInc = 0.;
    double offCentre = 0.;
    bool inside = false;
  };
  struct Roots {
    bool barrel = false;
    bool crosses = false;
    std::array<Root, 2> root{};
    std::size_t first = 0U;
  };
  [[nodiscard]] static Roots rootsOf(const BoundTrackState& state, const TrackingSurface& from, const TrackingSurface& to, PropagationDirection direction, const Propagator& propagator)
  {
    Roots out;
    const auto nx = static_cast<double>(to.normal[0]);
    const auto ny = static_cast<double>(to.normal[1]);
    const double transverse = std::hypot(nx, ny);
    out.barrel = transverse > 0.9;
    if (!out.barrel) {
      return out;
    }
    Scalar global[3]{};
    surfaceToGlobal(from, state.loc0(), state.loc1(), global);
    const auto x0 = static_cast<double>(global[0]);
    const auto y0 = static_cast<double>(global[1]);
    const auto z0 = static_cast<double>(global[2]);
    const auto phi = static_cast<double>(state.phi());
    const double sinTheta = std::sin(static_cast<double>(state.theta()));
    const double cosTheta = std::cos(static_cast<double>(state.theta()));
    const double kappa = propagator.curvature(static_cast<double>(state.qOverP()), sinTheta);
    if (!(std::abs(kappa) > 0.)) {
      return out;
    }
    const double dx = std::cos(phi);
    const double dy = std::sin(phi);
    const double cx = x0 - (dy / kappa);
    const double cy = y0 + (dx / kappa);
    const double radius = 1. / std::abs(kappa);
    const double ux = nx / transverse;
    const double uy = ny / transverse;
    const double tx = -uy;
    const double ty = ux;
    const auto ox = static_cast<double>(to.origin[0]);
    const auto oy = static_cast<double>(to.origin[1]);
    const auto oz = static_cast<double>(to.origin[2]);
    const double d = ((cx - ox) * ux) + ((cy - oy) * uy);
    if (std::abs(d) > radius) {
      return out;
    }
    out.crosses = true;
    const double half = std::sqrt((radius * radius) - (d * d));
    const double fx = cx - (d * ux);
    const double fy = cy - (d * uy);
    const auto sense = static_cast<double>(static_cast<int>(direction));
    for (std::size_t i = 0; i < 2U; ++i) {
      const double s = i == 0U ? -1. : 1.;
      const double px = fx + (s * half * tx);
      const double py = fy + (s * half * ty);
      const double ax = px - x0;
      const double ay = py - y0;
      const double along = (ax * dx) + (ay * dy);
      const double across = (-ax * dy) + (ay * dx);
      const double turn = std::atan2(kappa * along, 1. - (kappa * across));
      const double ell = std::abs(turn) < 1e-6 ? along : turn / kappa;
      auto& root = out.root[i];
      root.turn = turn;
      root.path = sense * ell / sinTheta;
      root.mu = (ax * tx) + (ay * ty);
      const double pz = z0 + (ell / sinTheta * cosTheta);
      const double rel[3] = {px - ox, py - oy, pz - oz};
      root.loc0 = (rel[0] * static_cast<double>(to.axis0[0])) + (rel[1] * static_cast<double>(to.axis0[1])) + (rel[2] * static_cast<double>(to.axis0[2]));
      root.loc1 = (rel[0] * static_cast<double>(to.axis1[0])) + (rel[1] * static_cast<double>(to.axis1[1])) + (rel[2] * static_cast<double>(to.axis1[2]));
      root.inside = to.bounds.isInside(static_cast<Scalar>(root.loc0), static_cast<Scalar>(root.loc1));
      const double phi1 = phi + turn;
      root.cosInc = (static_cast<double>(to.normal[0]) * sinTheta * std::cos(phi1)) + (static_cast<double>(to.normal[1]) * sinTheta * std::sin(phi1)) + (static_cast<double>(to.normal[2]) * cosTheta);
      root.offCentre = std::hypot(rel[0], rel[1], rel[2]);
    }
    out.first = std::abs(out.root[0].mu) <= std::abs(out.root[1].mu) ? 0U : 1U;
    return out;
  }

  enum InwardKind : std::uint8_t { kInwardOwn,
                                   kInwardWrong,
                                   kInwardControl,
                                   kInwardKinds };
  static constexpr std::array<const char*, kInwardKinds> kInwardKindNames{"failing, own rows", "failing, wrong rows", "control"};
  /// The seeds the failing inward passes are tried again from: as the fit has
  /// it, the outward end state with its own
  /// covariance and that times 100, and the diffuse seed with the MC truth's
  /// momentum and direction at the outermost row, its q/p alone, its
  /// direction alone.
  static constexpr std::array<const char*, 6> kSeedNames{"as the fit (diffuse)", "the outward covariance", "the outward covariance x 100", "MC momentum and direction, diffuse", "MC q/p only, diffuse", "MC direction only, diffuse"};

  /// A particle's charge [e] from its PDG code, by O2's particle table as
  /// the sink takes it (`Charge()` in thirds); 0 for a code it does not know.
  [[nodiscard]] static int chargeOf(int pdg)
  {
    const auto* particle = o2::O2DatabasePDG::Instance()->GetParticle(pdg);
    return particle == nullptr ? 0 : static_cast<int>(std::lround(particle->Charge() / 3.));
  }
  /// The sign relating the state's q/p to the charge: q/p = sign x charge / p
  /// (Propagator's curvature convention, checked by the control's agreement,
  /// `mChargeAgree`).
  static constexpr double kChargeSign = 1.;

  /// One refit's inward pass against the truth: its states
  /// to the CSV and the tables by step, the pass from other seeds, both roots
  /// of a failing step on a barrel plane.
  template <typename Material>
  void inwardStudy(InwardKind kind, std::uint64_t particle, const std::vector<std::uint32_t>& rows, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    const int n = static_cast<int>(rows.size());
    if (n < 3) {
      return;
    }
    const RowChain<Det> chain{.frame = &frame, .rows = rows.data()};
    const auto fit = fitMeasurements(propagator, geometry, material, chain, n, kinkBudget(parameters));
    if (fit.failedPass == FitPass::kOutward) {
      ++mInwardMismatch[kind];
      return;
    }
    const auto chainSeed = makeChainSeed(propagator, geometry, chain, n);
    BoundTrackState seed = fit.outwardState;
    setDiffuseCovariance(seed, chainSeed.qOverPSigma);
    const auto run = runInward(seed, rows, parameters, frame, geometry, material, propagator, true);
    if ((run.status == PropagationStatus::kOk) != (kind == kInwardControl)) {
      ++mInwardMismatch[kind];
      return;
    }
    ++mInwardCount[kind];
    const auto found = mTruth->particles.find(particle);
    const auto pt = static_cast<float>(truePt(found));
    const int pdg = found != mTruth->particles.end() ? found->second.pdg : 0;
    const auto index = mIndex.find(particle);
    const auto& own = mRowsOf.at(particle);
    const auto legOf = [&](std::uint32_t row) -> int {
      if (index == mIndex.end() || mLegs.empty() || mTimed[index->second] == 0U) {
        return -1;
      }
      const auto at = std::ranges::find(own, row);
      return at == own.end() ? -1 : static_cast<int>(mLegs[index->second][static_cast<std::size_t>(at - own.begin())]);
    };
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    const int charge = chargeOf(pdg);
    for (std::size_t k = 0; k < run.history.size(); ++k) {
      const auto& one = run.history[k];
      const auto& state = one.state;
      const auto* segment = mTruth->rows[one.row].carries(particle) ? segmentOfRow(particle, frame, one.row) : nullptr;
      const double p = 1. / std::abs(static_cast<double>(state.qOverP()));
      const double sigmaQOverP = std::sqrt(std::max(0., static_cast<double>(state.variance(kQOverP))));
      const double sigmaPhi = std::sqrt(std::max(0., static_cast<double>(state.variance(kPhi))));
      double pTrue = nan;
      double phiTrue = nan;
      double thetaTrue = nan;
      std::array<double, 2> locTrue{nan, nan};
      if (segment != nullptr) {
        pTrue = std::hypot(segment->momentum[0], segment->momentum[1], segment->momentum[2]);
        phiTrue = std::atan2(segment->momentum[1], segment->momentum[0]);
        thetaTrue = std::atan2(std::hypot(segment->momentum[0], segment->momentum[1]), segment->momentum[2]);
        locTrue = crossingOn(*segment, geometry.planeOf(state));
        auto& by = mInwardSteps[kind][std::min<std::size_t>(k, kInwardStepSlots - 1U)];
        // Signed: the state's q/p against the truth's charge over its
        // momentum; the charge's sign agreement.
        const double truthQOverP = charge != 0 ? static_cast<double>(charge) * kChargeSign / pTrue : nan;
        if (charge != 0 && k == 0U) {
          ++mChargeAgree[kind][(static_cast<double>(state.qOverP()) * truthQOverP) > 0. ? 0U : 1U];
        }
        by.ratio.push_back(p / pTrue);
        by.pullQOverP.push_back((static_cast<double>(state.qOverP()) - truthQOverP) / sigmaQOverP);
        const double dphi = std::remainder(static_cast<double>(state.phi()) - phiTrue, 2. * std::numbers::pi);
        by.dphi.push_back(1e3 * dphi);
        by.pullPhi.push_back(dphi / sigmaPhi);
        by.relSigma.push_back(sigmaQOverP * p);
      }
      // Both crossings of this step's plane from the state it started from:
      // whether the one the solve takes is the particle's.
      int stepRoot = -1;
      if (k >= 1U && segment != nullptr) {
        BoundTrackState from = run.history[k - 1U].state;
        if (k == 1U) {
          const auto start = measurementAt(frame, run.history[0].row);
          from.par[kLoc0] = start.loc0;
          from.par[kLoc1] = start.loc1;
        }
        const auto& to = geometry.planeOf(state);
        const auto roots = rootsOf(from, geometry.planeOf(from), to, PropagationDirection::kBackward, propagator);
        stepRoot = !roots.barrel ? 0 : 3;
        if (roots.barrel && roots.crosses) {
          const auto& first = roots.root[roots.first];
          const auto& second = roots.root[1U - roots.first];
          const Root* taken = first.path >= 0. ? &first : (second.path >= 0. && second.offCentre < first.offCentre ? &second : nullptr);
          if (taken != nullptr) {
            const Root& other = taken == &first ? second : first;
            stepRoot = std::hypot(taken->loc0 - locTrue[0], taken->loc1 - locTrue[1]) <= std::hypot(other.loc0 - locTrue[0], other.loc1 - locTrue[1]) ? 1 : 2;
            // Whether the sensor's bounds tell the two apart, and the turns.
            const auto turnBin = [](double turn) { return std::abs(turn) < std::numbers::pi / 2. ? "< 90" : "90-180"; };
            ++mStepRootDetail[std::format("{} | {} | taken {} the sensor, turn {} deg, path {} | the other {} the sensor, turn {} deg, path {}", kInwardKindNames[kind], stepRoot == 1 ? "the solve's the particle's" : "the other the particle's", taken->inside ? "inside" : "outside", turnBin(taken->turn), taken->path >= 0. ? "ahead" : "behind", other.inside ? "inside" : "outside", turnBin(other.turn), other.path >= 0. ? "ahead" : "behind")];
            // A rule among the crossings ahead:
            // the one nearer the plane's centre, against the particle's.
            if (other.path >= 0.) {
              const bool nearerIsTaken = taken->offCentre <= other.offCentre;
              ++mStepRootDetail[std::format("{} | both ahead: the crossing nearer the plane's centre is {}", kInwardKindNames[kind], nearerIsTaken == (stepRoot == 1) ? "the particle's" : "not the particle's")];
              const double margin = std::abs(taken->offCentre - other.offCentre);
              ++mStepRootDetail[std::format("{} | both ahead, {}: the two centres' distances differ by {}", kInwardKindNames[kind], stepRoot == 1 ? "the solve's the particle's" : "the other the particle's", margin < 1. ? "< 1 cm" : (margin < 5. ? "1-5 cm" : (margin < 20. ? "5-20 cm" : ">= 20 cm")))];
            }
          }
        }
        ++mStepRoots[kind][std::min<std::size_t>(k, kInwardStepSlots - 1U)][static_cast<std::size_t>(stepRoot)];
      }
      if (mInwardCsv) {
        *mInwardCsv << std::format("{},{},{:.4f},{},{},{},{},{},{},{},{},{:.5g},{:.5g},{:.6g},{:.6g},{:.6g},{:.6g},{:.4g},{:.4g},{:.4g},{:.5g},{:.5g},{:.5g},{:.5g},{:.6g},{},{}\n", static_cast<int>(kind), particle, static_cast<double>(pt), pdg, n, static_cast<int>(run.status), run.failedStep, k, one.row, static_cast<int>(mTruth->rows[one.row].layer), legOf(one.row), p, pTrue, static_cast<double>(state.phi()), phiTrue, static_cast<double>(state.theta()), thetaTrue, sigmaQOverP * p, sigmaPhi, static_cast<double>(one.chi2), static_cast<double>(state.loc0()), static_cast<double>(state.loc1()), locTrue[0], locTrue[1], static_cast<double>(state.qOverP()), charge, stepRoot);
      }
    }
    ruleStudy(kind, particle, rows, parameters, frame, geometry, material, propagator);
    if (kind == kInwardControl) {
      return;
    }
    // Both roots at the failing step, on a barrel plane.
    {
      const auto& to = geometry.planeOf(chain[run.failedIndex]);
      const auto roots = rootsOf(run.before, geometry.planeOf(run.before), to, PropagationDirection::kBackward, propagator);
      const std::string key = std::format("{} | status {}", kInwardKindNames[kind], static_cast<int>(run.status));
      if (!roots.barrel) {
        ++mRoots[key + " | a disc plane"];
      } else if (!roots.crosses) {
        ++mRoots[key + " | the circle misses the line"];
      } else {
        const auto& first = roots.root[roots.first];
        const auto& second = roots.root[1U - roots.first];
        std::string what;
        if (first.path >= 0.) {
          what = "the first crossing ahead: the solve takes it";
        } else if (!(second.path >= 0.)) {
          what = "both crossings behind";
        } else if (second.offCentre < first.offCentre) {
          what = "the second ahead and nearer the centre: the solve takes it";
        } else {
          what = std::format("the second ahead, farther from the centre (the guard), {} the sensor, turn {}", second.inside ? "inside" : "outside", std::abs(second.turn) < std::numbers::pi / 2. ? "< 90 deg" : "90-180 deg");
        }
        ++mRoots[key + " | " + what];
        // Which crossing the particle's own row is at (own rows only).
        if (const auto* segment = mTruth->rows[rows[static_cast<std::size_t>(run.failedIndex)]].carries(particle) ? segmentOfRow(particle, frame, rows[static_cast<std::size_t>(run.failedIndex)]) : nullptr; segment != nullptr) {
          const auto truth = crossingOn(*segment, to);
          const double toFirst = std::hypot(first.loc0 - truth[0], first.loc1 - truth[1]);
          const double toSecond = std::hypot(second.loc0 - truth[0], second.loc1 - truth[1]);
          ++mRoots[key + std::format(" | the MC crossing nearer the {} crossing", toFirst <= toSecond ? "first" : "second")];
        }
      }
    }
    if (kind != kInwardOwn) {
      return;
    }
    // The same pass from other seeds.
    const auto* outer = segmentOfRow(particle, frame, rows.back());
    for (std::size_t variant = 0; variant < kSeedNames.size(); ++variant) {
      BoundTrackState other = seed;
      if (variant == 1U || variant == 2U) {
        other = fit.outwardState;
        if (variant == 2U) {
          for (auto& value : other.cov) {
            value *= 100.F;
          }
        }
      } else if (variant >= 3U) {
        if (outer == nullptr) {
          ++mSeedOutcomes[variant][3];
          continue;
        }
        const double pTrue = std::hypot(outer->momentum[0], outer->momentum[1], outer->momentum[2]);
        if (variant != 5U) {
          other.par[kQOverP] = static_cast<Scalar>(charge != 0 ? static_cast<double>(charge) * kChargeSign / pTrue : (seed.qOverP() < 0.F ? -1. : 1.) / pTrue);
        }
        if (variant != 4U) {
          other.par[kPhi] = static_cast<Scalar>(std::atan2(outer->momentum[1], outer->momentum[0]));
          other.par[kTheta] = static_cast<Scalar>(std::atan2(std::hypot(outer->momentum[0], outer->momentum[1]), outer->momentum[2]));
        }
      }
      const auto again = runInward(other, rows, parameters, frame, geometry, material, propagator, false);
      // [passes, fails at the same step, fails later, fails earlier]
      std::size_t outcome = 0U;
      if (again.status != PropagationStatus::kOk) {
        outcome = again.failedStep == run.failedStep ? 1U : (again.failedStep > run.failedStep ? 2U : 3U);
      }
      ++mSeedOutcomes[variant][outcome];
    }
  }

  /// The crossing rules and the seeds compared: the refit of the same rows as `follow::Refit` takes it, under
  /// each `CrossingRule`, from the diffuse seed as production and from the
  /// outward end state with its own covariance.
  static constexpr std::array<CrossingRule, 3> kCrossingRules{CrossingRule::kFirstAhead, CrossingRule::kNearerCentre, CrossingRule::kInsideFirst};
  static constexpr std::array<const char*, 3> kCrossingRuleNames{"the first ahead (default)", "nearer the centre", "inside first"};
  static constexpr std::array<const char*, 2> kRuleSeedNames{"diffuse seed", "outward covariance"};
  enum RuleOutcome : std::uint8_t { kRuleOutward,
                                    kRuleInward,
                                    kRuleChi2,
                                    kRuleAccepted,
                                    kRuleOutcomes };
  struct RuleTally {
    std::array<std::size_t, kRuleOutcomes> outcome{};
    /// Against the default from the same seed: accepted by both, the state
    /// or chi2 not bit for bit; accepted by the default only; by this only.
    std::array<std::size_t, 3> against{};
    /// The accepted fits' state on the innermost row against the MC truth.
    std::vector<double> pullQOverP, pullPhi, pullTheta, pullLoc0, pullLoc1;
  };

  /// One chain's refits under every rule and seed (`kCrossingRules`), each
  /// against the default's from the same seed.
  template <typename Material>
  void ruleStudy(InwardKind kind, std::uint64_t particle, const std::vector<std::uint32_t>& rows, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    const int n = static_cast<int>(rows.size());
    const RowChain<Det> chain{.frame = &frame, .rows = rows.data()};
    const auto found = mTruth->particles.find(particle);
    const int charge = found != mTruth->particles.end() ? chargeOf(found->second.pdg) : 0;
    const auto* inner = mTruth->rows[rows.front()].carries(particle) ? segmentOfRow(particle, frame, rows.front()) : nullptr;
    struct Outcome {
      RuleOutcome outcome = kRuleOutward;
      BoundTrackState state;
      Scalar chi2 = 0.F;
    };
    const auto refit = [&]<CrossingRule Rule>(bool outwardCovariance) {
      Outcome out;
      const auto one = fitMeasurements<Rule>(propagator, geometry, material, chain, n, kinkBudget(parameters));
      if (one.failedPass == FitPass::kOutward) {
        return out;
      }
      out.state = one.state;
      out.chi2 = one.chi2;
      int ndf = one.ndf;
      bool ok = one.status == PropagationStatus::kOk;
      if (outwardCovariance) {
        std::int16_t failed = -1;
        std::int16_t kinkAt = -1;
        int kinks = 0;
        out.state = BoundTrackState{};
        ok = o2::gtrk::detail::runPass<Rule>(propagator, geometry, material, chain, n, PropagationDirection::kBackward, one.outwardState, out.state, out.chi2, failed, kinkBudget(parameters), kinks, kinkAt, NoFitObserver{}) == PropagationStatus::kOk;
        ndf = (2 * n) - kNBoundParams - (2 * kinks);
      }
      if (!ok) {
        out.outcome = kRuleInward;
        return out;
      }
      out.outcome = out.chi2 <= parameters.roadChi2MaxByNdf[ndf] ? kRuleAccepted : kRuleChi2;
      return out;
    };
    for (std::size_t seed = 0; seed < kRuleSeedNames.size(); ++seed) {
      std::array<Outcome, kCrossingRules.size()> outcomes { refit.template operator()<kCrossingRules[0]>(seed == 1U),
                                                            refit.template operator()<kCrossingRules[1]>(seed == 1U),
                                                            refit.template operator()<kCrossingRules[2]>(seed == 1U) };
      for (std::size_t rule = 0; rule < kCrossingRules.size(); ++rule) {
        auto& tally = mRuleTallies[kind][rule][seed];
        const auto& one = outcomes[rule];
        ++tally.outcome[one.outcome];
        const bool accepted = one.outcome == kRuleAccepted;
        const bool byDefault = outcomes[0].outcome == kRuleAccepted;
        if (accepted && byDefault) {
          // Bit identity is intentional; these records have no padding.
          // NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
          tally.against[0] += std::memcmp(&one.state, &outcomes[0].state, sizeof(BoundTrackState)) != 0 || std::bit_cast<std::uint32_t>(one.chi2) != std::bit_cast<std::uint32_t>(outcomes[0].chi2) ? 1U : 0U;
        } else if (byDefault) {
          ++tally.against[1];
        } else if (accepted) {
          ++tally.against[2];
        }
        if (!accepted || inner == nullptr || one.state.surface != geometry.planeOf(chain[0]).id) {
          continue;
        }
        const auto& state = one.state;
        const double pTrue = std::hypot(inner->momentum[0], inner->momentum[1], inner->momentum[2]);
        const double sigma[kNBoundParams] = {std::sqrt(static_cast<double>(state.variance(kLoc0))), std::sqrt(static_cast<double>(state.variance(kLoc1))), std::sqrt(static_cast<double>(state.variance(kPhi))), std::sqrt(static_cast<double>(state.variance(kTheta))), std::sqrt(static_cast<double>(state.variance(kQOverP)))};
        const auto locTrue = crossingOn(*inner, geometry.planeOf(state));
        tally.pullLoc0.push_back((static_cast<double>(state.loc0()) - locTrue[0]) / sigma[kLoc0]);
        tally.pullLoc1.push_back((static_cast<double>(state.loc1()) - locTrue[1]) / sigma[kLoc1]);
        tally.pullPhi.push_back(std::remainder(static_cast<double>(state.phi()) - std::atan2(inner->momentum[1], inner->momentum[0]), 2. * std::numbers::pi) / sigma[kPhi]);
        tally.pullTheta.push_back((static_cast<double>(state.theta()) - std::atan2(std::hypot(inner->momentum[0], inner->momentum[1]), inner->momentum[2])) / sigma[kTheta]);
        if (charge != 0) {
          tally.pullQOverP.push_back((static_cast<double>(state.qOverP()) - (static_cast<double>(charge) * kChargeSign / pTrue)) / sigma[kQOverP]);
        }
      }
    }
  }

  /// The clutter test on the follow: at
  /// each step of a follow whose core -- the road and the rows the walks
  /// took before -- is one particle's, every time-compatible row in each
  /// target's window refitted from the step's state: its chi2 against the
  /// prediction and the predicted residual covariance V. The non-own rows
  /// beyond the gate in the window give the local clutter density rho; the
  /// uniform expectation of strangers within chi2 t is rho pi t sqrt|V|,
  /// summed over the targets where the window leaves room for rho, against
  /// the strangers seen on the same targets by their lineage, and the
  /// held-out own row's chi2 against chi2(2). The follows are the
  /// iteration's kept tracks' and its dropped records' (`replayDropped`).
  static constexpr std::array<double, 5> kClutterEdges{1., 2., 4., 8., 0.}; ///< the last is the gate
  enum Lineage : std::uint8_t { kLineOwnHeldOut,
                                kLineOwnOverlap,
                                kLineOwnOtherLeg,
                                kLineNoise,
                                kLineDelta,
                                kLineDecay,
                                kLineDaughterOther,
                                kLineGrandDaughter,
                                kLineMother,
                                kLineSibling,
                                kLineUnknown,
                                kLineUnrelated,
                                kLineages };
  static constexpr std::array<const char*, kLineages> kLineageNames{"own, held out", "own, an overlap on its leg", "own, another leg", "no signal label", "daughter: delta ray", "daughter: decay", "daughter: other", "granddaughter", "its mother", "sibling", "lineage unknown", "unrelated"};
  /// Why a target has no held-out own row: the particle's row on the layer
  /// and the core's leg exists but is not in the window (outside it, or its
  /// transport failed); the leg has rows beyond the layer but none on it (a
  /// gap); its rows there are on another leg; it turned (rows on later legs
  /// only); or it has no row after the core's last at all.
  enum NoOwn : std::uint8_t { kNoOwnOutside,
                              kNoOwnGap,
                              kNoOwnOtherLeg,
                              kNoOwnTurned,
                              kNoOwnEnded,
                              kNoOwnReasons };
  static constexpr std::array<const char*, kNoOwnReasons> kNoOwnNames{"own row outside the window", "a gap on its leg", "own row on another leg", "turned: rows on later legs only", "past its last row"};
  struct ClutterTally {
    std::size_t targets = 0U;
    std::size_t empty = 0U;    ///< targets with no row in the window
    std::size_t measured = 0U; ///< targets whose window leaves room beyond the gate for rho
    std::array<double, kClutterEdges.size()> expected{};
    /// The same within the gate from the layer's mean density: its rows in
    /// the window's time range over its sensors' area (rho, read from the
    /// grid's counts).
    double expectedLayer = 0.;
    std::array<std::array<std::size_t, kClutterEdges.size()>, kLineages> seen{}; ///< on the measured targets, cumulative, chi2 <= edge
    std::size_t ownInBox = 0U;                                                   ///< measured targets with the held-out own row in the window
    std::array<std::size_t, kClutterEdges.size()> ownSeen{};                     ///< its chi2 <= edge
    std::size_t anyStrangerInGate = 0U;                                          ///< measured targets with a non-own row within the gate
    double anyPredicted = 0.;                                                    ///< sum of 1 - exp(-E(gate)) over them
    double sqrtDetV = 0.;                                                        ///< over the measured targets, for the mean
    /// The targets without the held-out own row, by why (`NoOwn`): all, the
    /// ones a stranger was taken at, an own row of another leg taken.
    std::array<std::array<std::size_t, 3>, kNoOwnReasons> noOwn{};
    /// The own row there but not seen: not visited (beyond the window's cells
    /// or time), visited outside the box, in the box but its transport failed.
    std::array<std::size_t, 3> outsideWhy{};
    /// The swaps, the held-out own row within the gate and a stranger (not
    /// the particle's) taken on its layer: the own row available to the step
    /// (it passed in production) or not (held, refused), by whether the step
    /// looked ahead: [available][plain, lookahead, kink]; of the available,
    /// the stranger's chi2 the lower, and the one-step evidence (-chi2/2 -
    /// log sqrt|V|) preferring the own row.
    std::array<std::array<std::size_t, 3>, 2> swaps{};
    std::array<std::size_t, 3> swapsTakenLower{};
    std::array<std::size_t, 3> swapsEvidenceOwn{};
    std::size_t swapsVolumeDecides = 0U; ///< available: the evidence prefers the own row while its chi2 is the larger
    std::size_t ownOtherTaken = 0U;      ///< the held-out row in the gate, an own row of an overlap or another leg taken
    /// The lookahead against the greedy choice at the steps that looked ahead: of the
    /// rows that passed (the production outcomes: free, within the gate),
    /// the smallest chi2 against what the lookahead took, each the particle's own on
    /// its leg or not: [greedy own][lookahead own].
    std::array<std::array<std::size_t, 2>, 2> versusGreedy{};
    /// The recomputed chi2 against production's, for its rows that passed
    /// or were refused by the gate: |relative difference| < 1e-4 / 1e-3 /
    /// 1e-2 / more; and the gate decided differently.
    std::array<std::size_t, 5> chi2Agreement{};
  };
  /// \a particle's mother: none for a primary (kNoTrackParticle), nullopt
  /// where it cannot be looked up (another source, an event or track beyond
  /// the kinematics file).
  [[nodiscard]] std::optional<std::uint64_t> motherOf(std::uint64_t particle) const
  {
    const auto label = workflow::labelOfParticle(particle);
    const int track = label.getTrackID();
    const int event = label.getEventID();
    if (label.getSourceID() != 0 || event < 0 || static_cast<std::size_t>(event) >= mMothers.size()) {
      return std::nullopt;
    }
    const auto& mothers = mMothers[static_cast<std::size_t>(event)];
    if (track < 0 || static_cast<std::size_t>(track) >= mothers.size()) {
      return std::nullopt;
    }
    if (mothers[static_cast<std::size_t>(track)] < 0) {
      return study::kNoTrackParticle;
    }
    // The mother is of the same event and source.
    return o2::MCCompLabel(mothers[static_cast<std::size_t>(track)], event, label.getSourceID()).getTrackEventSourceID();
  }
  /// A stranger's row against \a particle: every label of the row, the
  /// closest relation of any (a row is unrelated only if all its labels are).
  [[nodiscard]] Lineage lineageOf(std::uint64_t particle, std::uint32_t row) const
  {
    const auto& truth = mTruth->rows[row];
    if (truth.count == 0U) {
      return kLineNoise;
    }
    const auto mine = motherOf(particle);
    auto best = kLineUnrelated;
    for (std::uint8_t k = 0; k < truth.count; ++k) {
      const auto other = truth.particles[k];
      const auto theirs = motherOf(other);
      Lineage line = kLineUnrelated;
      if (!theirs || !mine) {
        line = kLineUnknown;
      } else if (*theirs == particle) {
        const auto found = mTruth->particles.find(other);
        const int process = found != mTruth->particles.end() ? static_cast<int>(found->second.process) : -1;
        line = process == 9 ? kLineDelta : (process == 4 ? kLineDecay : kLineDaughterOther); // TMCProcess kPDeltaRay, kPDecay
      } else if (const auto grand = *theirs != study::kNoTrackParticle ? motherOf(*theirs) : std::optional<std::uint64_t>{}; grand && *grand == particle) {
        line = kLineGrandDaughter;
      } else if (other == *mine) {
        line = kLineMother;
      } else if (*mine != study::kNoTrackParticle && *theirs == *mine) {
        line = kLineSibling;
      }
      best = std::min(best, line);
    }
    return best;
  }
  /// \a layer's group for the clutter tables: TRKFT3's VD, ML or OT barrel
  /// or its discs, else barrel or discs.
  [[nodiscard]] static std::string layerGroupOf(LayerID layer, const IterationParameters<Det>& parameters)
  {
    const bool barrel = ((parameters.barrelLayers >> layer) & 1U) != 0U;
    if constexpr (std::is_same_v<Det, detectors::TRKFT3>) {
      // TRKFT3.h's barrel: the VD layers (kVDRadius), then the five middle
      // (ML) and three outer (OT) layers.
      constexpr std::size_t kVD = detectors::trkft3::kVDRadius.size();
      return !barrel ? "discs" : (layer < kVD ? "VD" : (layer < kVD + 5U ? "ML" : "OT"));
    } else {
      return barrel ? "barrel" : "discs";
    }
  }
  /// The positions of a particle's rows (`own`, in its path order) beyond
  /// \a lastAt in the walk's direction: later along the path outward,
  /// earlier inward.
  [[nodiscard]] static std::vector<std::size_t> positionsBeyond(std::size_t lastAt, std::size_t size, PropagationDirection direction)
  {
    std::vector<std::size_t> out;
    if (lastAt >= size) {
      return out;
    }
    if (direction == PropagationDirection::kForward) {
      for (std::size_t k = lastAt + 1U; k < size; ++k) {
        out.push_back(k);
      }
    } else {
      for (std::size_t k = lastAt; k-- > 0U;) {
        out.push_back(k);
      }
    }
    return out;
  }
  /// A followed road for `clutterStudy`: its walks, the outward one first as
  /// the follow runs them (either may be absent), and its rows.
  struct ClutterWalks {
    const Walk* outward = nullptr;
    const Walk* inward = nullptr;
    ClusterList<kMaxClusters> road{};
    std::uint8_t roadCount = 0U;
    const char* population = "";
  };
  /// The dropped records' walks of the iteration being observed, kept by
  /// `replayDropped` for `clutterStudy`.
  struct DroppedWalks {
    Walk outward;
    Walk inward;
    ClusterList<kMaxClusters> road{};
    std::uint8_t roadCount = 0U;
    bool walked = false;
  };

  /// Each layer's rows by their time boxes' ends, and its sensors' area
  /// [cm2], for `layerDensity`; made once per timeframe.
  void layerOccupancy(const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry)
  {
    if (!mLayerLower.empty()) {
      return;
    }
    mLayerLower.assign(kLayers, {});
    mLayerUpper.assign(kLayers, {});
    mLayerArea.assign(kLayers, 0.);
    const auto columns = frame.measurements();
    std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> boxes(kLayers);
    for (std::uint32_t row = 0; row < mTruth->rows.size(); ++row) {
      const auto box = timeBoxAt(columns, row);
      boxes[static_cast<std::size_t>(mTruth->rows[row].layer)].emplace_back(box.lower(), box.upper());
    }
    // By start, each with its end: the grid's test (`visitIfCompatible`).
    for (std::size_t layer = 0; layer < kLayers; ++layer) {
      std::ranges::sort(boxes[layer]);
      for (const auto& [lower, upper] : boxes[layer]) {
        mLayerLower[layer].push_back(lower);
        mLayerUpper[layer].push_back(upper);
      }
    }
    for (const auto& surface : geometry.surfaces()) {
      if (static_cast<std::size_t>(surface.layer) < kLayers) {
        mLayerArea[surface.layer] += 4. * static_cast<double>(surface.bounds.halfX) * static_cast<double>(surface.bounds.halfY);
      }
    }
  }
  /// \a layer's rows the grid would let through \a window's time test --
  /// a start within [earliestPS, latestPS] and an end after afterPS
  /// (`visitIfCompatible`) -- over the layer's
  /// sensors' area [1/cm2].
  [[nodiscard]] double layerDensity(LayerID layer, const CandidateWindow& window) const
  {
    const auto& lower = mLayerLower[layer];
    const auto& upper = mLayerUpper[layer];
    if (lower.empty() || !(mLayerArea[layer] > 0.)) {
      return 0.;
    }
    const auto earliest = static_cast<std::uint32_t>(std::clamp<std::int64_t>(window.earliestPS, 0, std::numeric_limits<std::uint32_t>::max()));
    const auto latest = static_cast<std::uint32_t>(std::clamp<std::int64_t>(window.latestPS, 0, std::numeric_limits<std::uint32_t>::max()));
    const auto first = static_cast<std::size_t>(std::ranges::lower_bound(lower, earliest) - lower.begin());
    const auto last = static_cast<std::size_t>(std::ranges::upper_bound(lower, latest) - lower.begin());
    std::size_t compatible = 0U;
    for (std::size_t k = first; k < last; ++k) {
      compatible += upper[k] > window.afterPS ? 1U : 0U;
    }
    return static_cast<double>(compatible) / mLayerArea[layer];
  }

  /// The innovation of \a measurement against the prediction \a predicted
  /// (its covariance filled), in double: V = H P H^T + R on the two local
  /// coordinates, its determinant and r^T V^-1 r (NaN where the determinant
  /// is not positive).
  struct Innovation {
    double chi2 = 0.;
    double det = 0.;
  };
  template <typename Measurement>
  [[nodiscard]] static Innovation innovation(const BoundTrackState& predicted, const Measurement& measurement)
  {
    const double v00 = static_cast<double>(predicted.variance(kLoc0)) + static_cast<double>(measurement.varLoc0);
    const double v11 = static_cast<double>(predicted.variance(kLoc1)) + static_cast<double>(measurement.varLoc1);
    const double v01 = static_cast<double>(predicted.covAt(kLoc0, kLoc1)) + static_cast<double>(measurement.covLoc01);
    const double det = (v00 * v11) - (v01 * v01);
    if (!(det > 0.)) {
      return {.chi2 = std::numeric_limits<double>::quiet_NaN(), .det = det};
    }
    const double r0 = static_cast<double>(measurement.loc0) - static_cast<double>(predicted.loc0());
    const double r1 = static_cast<double>(measurement.loc1) - static_cast<double>(predicted.loc1());
    return {.chi2 = ((r0 * r0 * v11) - (2. * r0 * r1 * v01) + (r1 * r1 * v00)) / det, .det = det};
  }

  /// The area [cm2] of \a window on \a layer: the phi half width (at most
  /// pi) at the barrel layer's radius or the window's second centre, times
  /// the second coordinate's full width.
  [[nodiscard]] static double windowArea(const IterationParameters<Det>& parameters, LayerID layer, const CandidateWindow& window)
  {
    const bool barrel = ((parameters.barrelLayers >> layer) & 1U) != 0U;
    const double radius = barrel ? static_cast<double>(parameters.layerSpans[layer].position) : static_cast<double>(window.secondCentre);
    return 4. * std::min(static_cast<double>(window.phiHalfWidth), std::numbers::pi) * radius * static_cast<double>(window.secondHalfWidth);
  }

  /// The clutter test: the follow's steps from
  /// a clean core, every time-compatible row in each target's window
  /// refitted from the step's state; the strangers within the gate against
  /// a uniform-density expectation, the held-out own row, the swaps -- per
  /// population, iteration, pT, walk and layers (`printClutter`).
  template <typename Observation, typename Material>
  void clutterStudy(const Observation& iteration, const std::vector<Replay>& replayed, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    const auto gate = static_cast<double>(parameters.clusterChi2Max);
    auto edges = kClutterEdges;
    edges.back() = gate;
    layerOccupancy(frame, geometry);
    std::vector<ClutterWalks> walks;
    walks.reserve(replayed.size() + mDroppedWalks.size());
    for (const auto& replay : replayed) {
      walks.push_back({
        .outward = replay.outward.empty() ? nullptr : &replay.outward.front(),
        .inward = &replay,
        .road = replay.road,
        .roadCount = replay.roadCount,
        .population = "kept",
      });
    }
    for (const auto& dropped : mDroppedWalks) {
      if (dropped.walked) {
        walks.push_back({.outward = &dropped.outward, .inward = &dropped.inward, .road = dropped.road, .roadCount = dropped.roadCount, .population = "dropped"});
      }
    }
    std::mutex merge;
    o2::gtrk::detail::parallelFor(*mArena, 0U, walks.size(), 64, [&](std::size_t a, std::size_t b) {
      std::map<std::string, ClutterTally> local;
      for (std::size_t element = a; element < b; ++element) {
        const auto& followed = walks[element];
        if (followed.roadCount == 0U) {
          continue;
        }
        const std::span<const std::uint32_t> road{followed.road.index, followed.roadCount};
        if (mTruth->rows[road[0]].count == 0U) {
          continue;
        }
        const auto particle = mTruth->rows[road[0]].particles[0];
        if (!std::ranges::all_of(road, [&](std::uint32_t row) { return mTruth->rows[row].carries(particle); })) {
          continue;
        }
        const auto index = mIndex.find(particle);
        if (index == mIndex.end()) {
          continue; // the class's particles: their rows by leg
        }
        const auto& own = mRowsOf.at(particle);
        const bool timed = !mLegs.empty() && mTimed[index->second] != 0U;
        const auto positionOf = [&](std::uint32_t row) -> std::size_t { return static_cast<std::size_t>(std::ranges::find(own, row) - own.begin()); };
        const auto legAt = [&](std::size_t k) -> int { return !timed || k >= own.size() ? -1 : static_cast<int>(mLegs[index->second][k]); };
        const auto legOf = [&](std::uint32_t row) -> int { return legAt(positionOf(row)); };
        const std::string ptName = kPathPtNames[pathPtBin(particle)];
        bool clean = true;
        std::uint32_t last = road.back();
        std::size_t core = road.size();
        const auto walkOf = [&](const Walk& walk, PropagationDirection direction) {
          last = direction == PropagationDirection::kForward ? road.back() : road.front();
          for (const auto& step : walk.steps) {
            if (!clean) {
              return;
            }
            const int coreLeg = legOf(last);
            const std::size_t lastAt = positionOf(last);
            for (const auto& target : step.targets) {
              auto& tally = local[std::format("{} | {} | {} | {} | {} | core {}", followed.population, iteration.iteration, ptName, direction == PropagationDirection::kForward ? "outward" : "inward", layerGroupOf(target.layer, parameters), std::min<std::size_t>(core, 9U))];
              ++tally.targets;
              const double boxArea = windowArea(parameters, target.layer, target.window);
              struct Seen {
                std::uint32_t row;
                double chi2;
                double sqrtDet;
                StepOutcome outcome;
              };
              std::vector<Seen> inBox;
              for (const auto& candidate : step.candidates) {
                if (candidate.layer != target.layer) {
                  continue;
                }
                const auto measurement = frame.measurement(candidate.row);
                if (!spatiallyCompatible(target.window, measurement.phi, secondCoordinate<Det>(target.layer, measurement.r, measurement.z))) {
                  continue;
                }
                const auto at = measurementAt(frame, candidate.row);
                BoundTrackState predicted;
                CovarianceFactors factors;
                if (!follow::transportCandidate<false>(propagator, geometry, material, geometry.planeOf(at), direction, step.state, step.factors, predicted, factors, 1.F)) {
                  continue;
                }
                covarianceFromFactors(factors, predicted.cov);
                const auto [chi2, det] = innovation(predicted, at);
                if (!(det > 0.)) {
                  continue;
                }
                inBox.push_back({.row = candidate.row, .chi2 = chi2, .sqrtDet = std::sqrt(det), .outcome = candidate.outcome});
                if (candidate.outcome == StepOutcome::kPasses || candidate.outcome == StepOutcome::kChi2) {
                  const auto production = static_cast<double>(candidate.chi2);
                  const double relative = std::abs(chi2 - production) / std::max(production, 1e-6);
                  ++tally.chi2Agreement[relative < 1e-4 ? 0U : (relative < 1e-3 ? 1U : (relative < 1e-2 ? 2U : 3U))];
                  tally.chi2Agreement[4] += (chi2 <= gate) != (candidate.outcome == StepOutcome::kPasses) ? 1U : 0U;
                }
              }
              // The held-out own row: the particle's on this layer and the
              // core's leg, the nearest to the prediction if there are two.
              std::uint32_t heldOut = kNoRow;
              double heldChi2 = std::numeric_limits<double>::infinity();
              for (const auto& one : inBox) {
                if (mTruth->rows[one.row].carries(particle) && legOf(one.row) == coreLeg && one.chi2 < heldChi2) {
                  heldOut = one.row;
                  heldChi2 = one.chi2;
                }
              }
              const auto lineage = [&](std::uint32_t row) {
                if (row == heldOut) {
                  return kLineOwnHeldOut;
                }
                if (mTruth->rows[row].carries(particle)) {
                  return legOf(row) == coreLeg ? kLineOwnOverlap : kLineOwnOtherLeg;
                }
                return lineageOf(particle, row);
              };
              const bool tookHere = step.choice != kNoRow && mTruth->rows[step.choice].layer == target.layer;
              if (inBox.empty()) {
                ++tally.empty;
              }
              if (heldOut == kNoRow) {
                // Why: the particle's rows on this layer and beyond the core's
                // last, by leg, in its path order.
                bool sameLegHere = false;
                bool otherLegHere = false;
                bool sameLegBeyond = false;
                bool laterLegs = false;
                const auto beyond = positionsBeyond(lastAt, own.size(), direction);
                for (const auto k : beyond) {
                  const bool here = mTruth->rows[own[k]].layer == target.layer;
                  const bool sameLeg = legAt(k) == coreLeg;
                  sameLegHere = sameLegHere || (here && sameLeg);
                  otherLegHere = otherLegHere || (here && !sameLeg);
                  sameLegBeyond = sameLegBeyond || sameLeg;
                  laterLegs = laterLegs || !sameLeg;
                }
                const NoOwn why = sameLegHere ? kNoOwnOutside : (sameLegBeyond ? kNoOwnGap : (otherLegHere ? kNoOwnOtherLeg : (laterLegs ? kNoOwnTurned : kNoOwnEnded)));
                auto& counts = tally.noOwn[why];
                ++counts[0];
                if (why == kNoOwnOutside) {
                  std::uint32_t mine = kNoRow;
                  for (const auto k : beyond) {
                    if (mTruth->rows[own[k]].layer == target.layer && legAt(k) == coreLeg) {
                      mine = own[k];
                      break;
                    }
                  }
                  const auto visited = std::ranges::find(step.candidates, mine, &Candidate::row);
                  std::size_t miss = 0U;
                  if (visited != step.candidates.end()) {
                    const auto measurement = frame.measurement(mine);
                    miss = spatiallyCompatible(target.window, measurement.phi, secondCoordinate<Det>(target.layer, measurement.r, measurement.z)) ? 2U : 1U;
                  }
                  ++tally.outsideWhy[miss];
                }
                if (tookHere) {
                  counts[1] += mTruth->rows[step.choice].carries(particle) ? 0U : 1U;
                  counts[2] += mTruth->rows[step.choice].carries(particle) ? 1U : 0U;
                }
              }
              if (inBox.empty()) {
                continue;
              }
              double sqrtDet = 0.;
              for (const auto& one : inBox) {
                sqrtDet += one.sqrtDet;
              }
              sqrtDet /= static_cast<double>(inBox.size());
              const double ellipse = std::numbers::pi * gate * sqrtDet;
              if (boxArea > 1.25 * ellipse) {
                // The measured targets: the strangers seen and expected on the same ones.
                ++tally.measured;
                tally.sqrtDetV += sqrtDet;
                std::size_t outside = 0U;
                bool stranger = false;
                for (const auto& one : inBox) {
                  const auto line = lineage(one.row);
                  const bool ownRow = line <= kLineOwnOtherLeg;
                  if (one.chi2 > gate) {
                    outside += ownRow ? 0U : 1U;
                    continue;
                  }
                  stranger = stranger || !ownRow;
                  for (std::size_t e = 0; e < edges.size(); ++e) {
                    tally.seen[line][e] += one.chi2 <= edges[e] ? 1U : 0U;
                  }
                }
                const double rho = static_cast<double>(outside) / (boxArea - ellipse);
                for (std::size_t e = 0; e < edges.size(); ++e) {
                  tally.expected[e] += rho * std::numbers::pi * edges[e] * sqrtDet;
                }
                tally.expectedLayer += layerDensity(target.layer, target.window) * ellipse;
                tally.anyStrangerInGate += stranger ? 1U : 0U;
                tally.anyPredicted += 1. - std::exp(-rho * ellipse);
                if (heldOut != kNoRow) {
                  ++tally.ownInBox;
                  for (std::size_t e = 0; e < edges.size(); ++e) {
                    tally.ownSeen[e] += heldChi2 <= edges[e] ? 1U : 0U;
                  }
                }
              }
              const auto takenSeen = std::ranges::find(inBox, step.choice, &Seen::row);
              const auto ownSeen = std::ranges::find(inBox, heldOut, &Seen::row);
              if (heldOut != kNoRow && heldChi2 <= gate && tookHere && step.choice != heldOut && takenSeen != inBox.end()) {
                if (mTruth->rows[step.choice].carries(particle)) {
                  ++tally.ownOtherTaken;
                } else {
                  const std::size_t how = step.kink ? 2U : (step.lookahead ? 1U : 0U);
                  const bool available = ownSeen->outcome == StepOutcome::kPasses;
                  ++tally.swaps[available ? 1U : 0U][how];
                  if (available) {
                    tally.swapsTakenLower[how] += takenSeen->chi2 < ownSeen->chi2 ? 1U : 0U;
                    const double evidenceOwn = (-0.5 * ownSeen->chi2) - std::log(ownSeen->sqrtDet);
                    const double evidenceTaken = (-0.5 * takenSeen->chi2) - std::log(takenSeen->sqrtDet);
                    tally.swapsEvidenceOwn[how] += evidenceOwn > evidenceTaken ? 1U : 0U;
                    tally.swapsVolumeDecides += evidenceOwn > evidenceTaken && ownSeen->chi2 > takenSeen->chi2 ? 1U : 0U;
                  }
                }
              }
              if (step.lookahead && tookHere) {
                std::uint32_t greedy = kNoRow;
                Scalar best = std::numeric_limits<Scalar>::infinity();
                for (const auto& candidate : step.candidates) {
                  if (candidate.layer == target.layer && candidate.outcome == StepOutcome::kPasses && (candidate.chi2 < best || (!(candidate.chi2 > best) && candidate.row < greedy))) {
                    best = candidate.chi2;
                    greedy = candidate.row;
                  }
                }
                const auto ownOnLeg = [&](std::uint32_t row) { return mTruth->rows[row].carries(particle) && legOf(row) == coreLeg ? 1U : 0U; };
                if (greedy != kNoRow) {
                  ++tally.versusGreedy[ownOnLeg(greedy)][ownOnLeg(step.choice)];
                }
              }
            }
            if (step.choice != kNoRow) {
              clean = mTruth->rows[step.choice].carries(particle);
              last = step.choice;
              ++core;
            }
          }
        };
        if (followed.outward != nullptr) {
          walkOf(*followed.outward, PropagationDirection::kForward);
        }
        if (followed.inward != nullptr) {
          walkOf(*followed.inward, PropagationDirection::kBackward);
        }
      }
      const std::scoped_lock lock(merge);
      for (const auto& [key, one] : local) {
        auto& sum = mClutter[key];
        sum.targets += one.targets;
        sum.empty += one.empty;
        sum.measured += one.measured;
        sum.ownInBox += one.ownInBox;
        sum.anyStrangerInGate += one.anyStrangerInGate;
        sum.anyPredicted += one.anyPredicted;
        sum.sqrtDetV += one.sqrtDetV;
        sum.expectedLayer += one.expectedLayer;
        for (std::size_t w = 0; w < 3U; ++w) {
          sum.outsideWhy[w] += one.outsideWhy[w];
        }
        for (std::size_t e = 0; e < edges.size(); ++e) {
          sum.expected[e] += one.expected[e];
          sum.ownSeen[e] += one.ownSeen[e];
          for (std::size_t l = 0; l < kLineages; ++l) {
            sum.seen[l][e] += one.seen[l][e];
          }
        }
        for (std::size_t r = 0; r < kNoOwnReasons; ++r) {
          for (std::size_t c = 0; c < 3U; ++c) {
            sum.noOwn[r][c] += one.noOwn[r][c];
          }
        }
        for (std::size_t h = 0; h < 3U; ++h) {
          sum.swaps[0][h] += one.swaps[0][h];
          sum.swaps[1][h] += one.swaps[1][h];
          sum.swapsTakenLower[h] += one.swapsTakenLower[h];
          sum.swapsEvidenceOwn[h] += one.swapsEvidenceOwn[h];
        }
        sum.swapsVolumeDecides += one.swapsVolumeDecides;
        sum.ownOtherTaken += one.ownOtherTaken;
        for (std::size_t g = 0; g < 2U; ++g) {
          for (std::size_t r = 0; r < 2U; ++r) {
            sum.versusGreedy[g][r] += one.versusGreedy[g][r];
          }
        }
        for (std::size_t c = 0; c < sum.chi2Agreement.size(); ++c) {
          sum.chi2Agreement[c] += one.chi2Agreement[c];
        }
      }
    });
  }

  /// `clutterStudy`'s tables.
  void printClutter() const
  {
    if (mClutter.empty()) {
      return;
    }
    std::string text = std::format(
      "  the clutter test: the follow's steps from a clean core, every time-compatible row in each target's window refitted "
      "from the step's state; per population (the iteration's kept tracks, its dropped records) | iteration | pT | walk | layers:\n"
      "    targets (empty, measured: room for rho beyond the gate) | on the measured: mean sqrt|V| [cm2] | the held-out own row in the window | "
      "strangers within chi2 <= 1 / 2 / 4 / 8 / gate, every non-own row (the unrelated among them) against the uniform expectation rho pi t sqrt|V| (rho "
      "from the window's non-own rows beyond the gate), and within the gate from the layer's mean density | targets with a stranger in the gate: seen / "
      "Poisson | the held-out row's chi2 <= 1 / 2 / 4 / 8 / "
      "gate as a fraction (chi2(2): {:.3f} / {:.3f} / {:.3f} / {:.3f} / gate) | the other lineages within the gate",
      1. - std::exp(-0.5), 1. - std::exp(-1.), 1. - std::exp(-2.), 1. - std::exp(-4.));
    for (const auto& [key, t] : mClutter) {
      text += std::format("\n    {}: {} ({}, {}) | {:.4f} | {} |", key, t.targets, t.empty, t.measured, t.measured != 0U ? t.sqrtDetV / static_cast<double>(t.measured) : 0., t.ownInBox);
      for (std::size_t e = 0; e < kClutterEdges.size(); ++e) {
        std::size_t strangers = 0U;
        for (std::size_t l = kLineNoise; l < kLineages; ++l) {
          strangers += t.seen[l][e];
        }
        text += std::format(" {} ({}) / {:.1f}{}", strangers, t.seen[kLineUnrelated][e], t.expected[e], e + 1U < kClutterEdges.size() ? ";" : " |");
      }
      text += std::format(" layer {:.1f} | {} / {:.1f} |", t.expectedLayer, t.anyStrangerInGate, t.anyPredicted);
      for (std::size_t e = 0; e < kClutterEdges.size(); ++e) {
        text += std::format(" {:.3f}", t.ownInBox != 0U ? static_cast<double>(t.ownSeen[e]) / static_cast<double>(t.ownInBox) : 0.);
      }
      text += " |";
      for (std::size_t l = kLineOwnOverlap; l < kLineUnrelated; ++l) {
        text += t.seen[l].back() != 0U ? std::format(" {} {};", kLineageNames[l], t.seen[l].back()) : std::string{};
      }
    }
    text +=
      "\n    the same keys: the targets without the held-out own row, by why: targets / a stranger taken / an own row of another leg taken | the swaps "
      "(the held-out row in the gate, a stranger taken), the own row not available to the step (held, refused) plain / lookahead / kink, and "
      "available: plain / lookahead / kink, each (the stranger's chi2 the lower, the one-step evidence -chi2/2 - log sqrt|V| preferring the own "
      "row); the own preferred while its chi2 is the larger | an own row of an overlap or another leg taken instead | the lookahead against the greedy "
      "smallest chi2 among the passing rows: both own / greedy own, lookahead not (it costs) / greedy not, lookahead own (it saves) / neither | the "
      "recomputed chi2 against production's: |relative difference| < 1e-4 / 1e-3 / 1e-2 / more, the gate decided otherwise";
    for (const auto& [key, t] : mClutter) {
      text += std::format("\n    {}:", key);
      for (std::size_t r = 0; r < kNoOwnReasons; ++r) {
        text += t.noOwn[r][0] != 0U ? std::format(" {} {} / {} / {};", kNoOwnNames[r], t.noOwn[r][0], t.noOwn[r][1], t.noOwn[r][2]) : std::string{};
      }
      text += std::format(" (seen not: not visited {}, outside the box {}, transport {})", t.outsideWhy[0], t.outsideWhy[1], t.outsideWhy[2]);
      text += std::format(" | not available {} / {} / {}; available", t.swaps[0][0], t.swaps[0][1], t.swaps[0][2]);
      for (std::size_t h = 0; h < 3U; ++h) {
        text += std::format(" {} ({}, {}){}", t.swaps[1][h], t.swapsTakenLower[h], t.swapsEvidenceOwn[h], h < 2U ? " /" : ";");
      }
      text += std::format(" {} | {} | {} / {} / {} / {} | {} / {} / {} / {}, {}", t.swapsVolumeDecides, t.ownOtherTaken, t.versusGreedy[1][1], t.versusGreedy[1][0], t.versusGreedy[0][1], t.versusGreedy[0][0], t.chi2Agreement[0], t.chi2Agreement[1], t.chi2Agreement[2], t.chi2Agreement[3], t.chi2Agreement[4]);
    }
    std::cout << text << '\n';
  }

  /// The first wrong rows' alternatives on the follow: at the first
  /// wrong row of each followed road whose road is one particle's -- the
  /// iteration's dropped records, and its kept tracks with one wrong row --
  /// the three alternatives of that step: the row taken (T), the
  /// particle's own row on that layer and leg (O), and stopping there (D,
  /// the prefix). For T and O, from the step's state: the chi2, `|V|`, the
  /// window's and the layer's clutter density, and the one-step evidence
  /// against stopping, log(g / rho) + log(p / (1 - p)), with p the chance of
  /// reaching a live sensor of the layer, sampled from the state's
  /// covariance. Each alternative's chain refitted as `follow::Refit`, and the
  /// two particles' fates; each chain's evidence (the innovations, the seed prior's determinant, flat in the
  /// diffuse limit), with the seed's covariance as made and 10 and 100 times
  /// wider: whether the alternatives' order holds; and the prefix's evidence
  /// fitted alone against the same rows' running sum inside the longer chain.
  enum Decision : std::uint8_t { kDecideTaken,
                                 kDecideOwn,
                                 kDecideStop,
                                 kDecisions };
  static constexpr std::array<const char*, kDecisions> kDecisionNames{"taken", "own", "stop"};
  enum OwnThere : std::uint8_t { kOwnAvailable,
                                 kOwnHeld,
                                 kOwnRefused,
                                 kOwnNotSeen,
                                 kOwnNone,
                                 kOwnTheres };
  static constexpr std::array<const char*, kOwnTheres> kOwnThereNames{"passing", "held by a track", "refused (gate, update)", "not in the window", "none on the layer and leg"};
  struct DecisionTally {
    std::size_t events = 0U;
    std::array<std::size_t, kOwnTheres> own{};
    std::array<std::size_t, kLineages> taken{};
    /// The one-step evidence's choice among T, O (where seen) and stopping:
    /// with the window's rho, with the layer's.
    std::array<std::size_t, kDecisions> stepLocal{};
    std::array<std::size_t, kDecisions> stepLayer{};
    std::array<std::size_t, 4> reach{}; ///< p < 0.25 / < 0.5 / < 0.9 / >= 0.9
    double reachSum = 0.;
    /// The refits: T accepted, O accepted (where there is one), the prefix
    /// accepted and as long as the length cut; refused or failed counted apart.
    std::size_t refitTaken = 0U;
    std::size_t refitOwn = 0U;
    std::size_t ownChains = 0U;
    std::size_t refitPrefix = 0U;
    std::size_t prefixLong = 0U;
    /// The particles at the iteration's end: the road's found (by any track),
    /// the taken row's found; the taken row's being a class particle.
    std::size_t particleFound = 0U;
    std::size_t strangerFound = 0U;
    std::size_t strangerParticles = 0U;
    /// The chains' evidence (S, the clutter of the last row
    /// and the reach included): the choice with the seed's covariance as
    /// made, the same choice at 10 and 100 times; it against the one-step
    /// choice; the prefix alone against its running sum (|difference| < 0.1
    /// / < 1 / more).
    std::array<std::size_t, kDecisions> chainChoice{};
    std::array<std::size_t, 2> chainStable{};
    std::size_t chainAgreesStep = 0U;
    std::size_t chains = 0U;
    std::array<std::size_t, 3> prefixInherited{};
  };

  /// The evidence of \a rows: the outward filter
  /// from the chain's seed (`makeChainSeed`), its covariance times \a scale,
  /// each step's innovation density (-chi2/2 - log sqrt|V| - log 2 pi), plus
  /// half the seed's log-determinant and 5/2 log 2 pi; no kinks. \a running,
  /// where given, gets the sum after each row (without the seed's terms).
  /// NaN where a step fails.
  template <typename Material>
  [[nodiscard]] double chainEvidence(const std::vector<std::uint32_t>& rows, double scale, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, std::vector<double>* running = nullptr) const
  {
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    const int n = static_cast<int>(rows.size());
    if (n < 3) {
      return nan;
    }
    const RowChain<Det> chain{.frame = &frame, .rows = rows.data()};
    BoundTrackState state = makeChainSeed(propagator, geometry, chain, n).state;
    for (auto& value : state.cov) {
      value = static_cast<Scalar>(static_cast<double>(value) * scale);
    }
    // log|P0| by Cholesky, in double.
    std::array<double, static_cast<std::size_t>(kNBoundParams) * kNBoundParams> l{}; // NOLINT(misc-const-correctness): written through at(), which the check does not follow
    const auto at = [](int row, int column) { return (static_cast<std::size_t>(row) * kNBoundParams) + static_cast<std::size_t>(column); };
    double logDet = 0.;
    for (int i = 0; i < kNBoundParams; ++i) {
      for (int j = 0; j <= i; ++j) {
        auto sum = static_cast<double>(state.cov[covIndex(i, j)]);
        for (int k = 0; k < j; ++k) {
          sum -= l[at(i, k)] * l[at(j, k)];
        }
        if (i == j) {
          if (!(sum > 0.)) {
            return nan;
          }
          l[at(i, i)] = std::sqrt(sum);
          logDet += std::log(sum);
        } else {
          l[at(i, j)] = sum / l[at(j, j)];
        }
      }
    }
    auto factors = factorCovariance(state.cov);
    double evidence = 0.;
    const double log2pi = std::log(2. * std::numbers::pi);
    for (int k = 0; k < n; ++k) {
      const auto measurement = chain[k];
      if (k > 0) {
        BoundTrackState next;
        CovarianceFactors nextFactors;
        if (!follow::transportCandidate<false>(propagator, geometry, material, geometry.planeOf(measurement), PropagationDirection::kForward, state, factors, next, nextFactors, 1.F)) {
          return nan;
        }
        state = next;
        factors = nextFactors;
      }
      BoundTrackState predicted = state;
      covarianceFromFactors(factors, predicted.cov);
      const auto [chi2, det] = innovation(predicted, measurement);
      if (!(det > 0.)) {
        return nan;
      }
      evidence += (-0.5 * chi2) - (0.5 * std::log(det)) - log2pi;
      const auto update = updateWithMeasurement(state, factors, measurement);
      if (!updateUsable(update, factors)) {
        return nan;
      }
      if (running != nullptr) {
        running->push_back(evidence);
      }
    }
    return evidence + (0.5 * logDet) + (2.5 * log2pi);
  }

  /// The chance that \a state, with the covariance \a factors, crosses a
  /// live sensor of \a layer in \a direction: \a samples draws of the state,
  /// each carried to the sensors listed within \a window about its centre
  /// (as `FollowStage::crossesSensor`); clipped away from 0 and 1 by half a
  /// draw.
  [[nodiscard]] double reachChance(const IterationParameters<Det>& parameters, const BoundTrackState& state, const CovarianceFactors& factors, LayerID layer, const CandidateWindow& window, PropagationDirection direction, const TrackingGeometry<Det>& geometry, const Propagator& propagator, std::uint32_t seed, int samples = 64) const
  {
    BoundTrackState mean = state;
    covarianceFromFactors(factors, mean.cov);
    std::array<double, static_cast<std::size_t>(kNBoundParams) * kNBoundParams> l{}; // NOLINT(misc-const-correctness): written through at(), which the check does not follow
    const auto at = [](int row, int column) { return (static_cast<std::size_t>(row) * kNBoundParams) + static_cast<std::size_t>(column); };
    for (int i = 0; i < kNBoundParams; ++i) {
      for (int j = 0; j <= i; ++j) {
        auto sum = static_cast<double>(mean.cov[covIndex(i, j)]);
        for (int k = 0; k < j; ++k) {
          sum -= l[at(i, k)] * l[at(j, k)];
        }
        if (i == j) {
          l[at(i, i)] = sum > 0. ? std::sqrt(sum) : 0.;
        } else {
          const double d = l[at(j, j)];
          l[at(i, j)] = d > 0. ? sum / d : 0.;
        }
      }
    }
    Scalar crossing[3]{};
    follow::crossingPoint(parameters, layer, {.phi = window.phiCentre, .second = window.secondCentre}, crossing);
    const auto& from = geometry.planeOf(state);
    std::mt19937 rng(seed);
    std::normal_distribution<double> gauss(0., 1.);
    int crossed = 0;
    for (int s = 0; s < samples; ++s) {
      std::array<double, kNBoundParams> normal{};
      for (auto& x : normal) {
        x = gauss(rng);
      }
      BoundTrackState sample = state;
      for (int i = 0; i < kNBoundParams; ++i) {
        double shift = 0.;
        for (int k = 0; k <= i; ++k) {
          shift += l[at(i, k)] * normal[static_cast<std::size_t>(k)];
        }
        sample.par[i] = static_cast<Scalar>(static_cast<double>(state.par[i]) + shift);
      }
      const auto phiReach = window.phiHalfWidth;
      const auto secondReach = window.secondHalfWidth;
      const bool hit = mCells.view().anySensor(layer, window.phiCentre, phiReach, window.secondCentre, secondReach, [&](std::uint32_t id) {
        const TrackingSurface& surface = geometry.surface(static_cast<SensorID>(id));
        TrackingSurface plane = surface;
        if (surface.radius > 0.F) {
          // The tangent plane at the stage's nominal layer crossing.
          Scalar along0 = 0.F;
          Scalar along1 = 0.F;
          for (int i = 0; i < 3; ++i) {
            along0 += (crossing[i] - surface.origin[i]) * surface.axis0[i];
            along1 += (crossing[i] - surface.origin[i]) * surface.axis1[i];
          }
          plane = tangentPlane(surface, o2::gpu::CAMath::ATan2(along1, along0));
        }
        BoundTrackState at = sample;
        const auto result = propagate(propagator, at, from, plane, direction);
        // The half-turn veto at the sample's own curvature.
        const ExactScalar curvature = propagator.curvature(static_cast<ExactScalar>(sample.qOverP()), 1.);
        return result.status == PropagationStatus::kOk && !follow::beyondHalfTurn(curvature * static_cast<ExactScalar>(result.pathLength));
      });
      crossed += hit ? 1 : 0;
    }
    const double half = 0.5 / static_cast<double>(samples);
    return std::clamp(static_cast<double>(crossed) / static_cast<double>(samples), half, 1. - half);
  }

  /// \a row's update from \a state (its \a factors) as one step of the
  /// follow makes it, carried to the row's own sensor (`transportCandidate`,
  /// nominal scattering): the chi2 and sqrt|V|, V = H P H^T + R the
  /// residual's predicted covariance. With
  /// \a beamLine the prediction is first updated with the beam line
  /// (`updateWithBeamLine`, the scattering inside the row's layer): the row
  /// against the beam line as the walk's reference. Without
  /// it also |V| as the follow's EndDensity has it, the product of the
  /// update's float innovations, and from the carried factors in double
  /// (the variant's quantity against its oracle).
  struct OneStep {
    bool valid = false;
    double chi2 = 0.;
    double sqrtDet = 0.;
    double floatDet = -1.;
    double factorDet = -1.;
  };
  /// The window's clutter beyond the nominal gate; callers choose which
  /// rows to exclude and the value for an area too small to measure it.
  template <typename Material, typename Skip>
  [[nodiscard]] double windowRho(const IterationParameters<Det>& parameters, const Step& step, LayerID layer, const CandidateWindow& window, double gate, double sqrtDet, PropagationDirection direction, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, const Skip& skip, double missing) const
  {
    const double boxArea = windowArea(parameters, layer, window);
    std::size_t outside = 0U;
    for (const auto& candidate : step.candidates) {
      if (candidate.layer != layer || skip(candidate.row)) {
        continue;
      }
      const auto one = oneStep(candidate.row, step.state, step.factors, direction, frame, geometry, material, propagator);
      const auto measurement = frame.measurement(candidate.row);
      outside += one.valid && one.chi2 > gate && spatiallyCompatible(window, measurement.phi, secondCoordinate<Det>(layer, measurement.r, measurement.z)) ? 1U : 0U;
    }
    const double ellipse = std::numbers::pi * gate * sqrtDet;
    return boxArea > 1.25 * ellipse ? std::max(static_cast<double>(outside), 0.5) / (boxArea - ellipse) : missing;
  }

  template <typename Material>
  [[nodiscard]] OneStep oneStep(std::uint32_t row, const BoundTrackState& state, const CovarianceFactors& factors, PropagationDirection direction, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, const IterationParameters<Det>* beamLine = nullptr) const
  {
    const auto measurement = measurementAt(frame, row);
    BoundTrackState predicted;
    CovarianceFactors carried;
    if (!follow::transportCandidate<false>(propagator, geometry, material, geometry.planeOf(measurement), direction, state, factors, predicted, carried, 1.F)) {
      return {};
    }
    double floatDet = -1.;
    double factorDet = -1.;
    if (beamLine == nullptr) {
      BoundTrackState updated = predicted;
      CovarianceFactors updatedFactors = carried;
      Scalar determinant = 0.F;
      static_cast<void>(updateWithMeasurement<true>(updated, updatedFactors, measurement, &determinant));
      floatDet = static_cast<double>(determinant);
      const auto covariance = [&](int i, int j) {
        double sum = 0.;
        for (int k = std::max(i, j); k < kNBoundParams; ++k) {
          const double ui = k == i ? 1. : static_cast<double>(carried.upper[o2::gtrk::detail::strictUpperIndex(i, k)]);
          const double uj = k == j ? 1. : static_cast<double>(carried.upper[o2::gtrk::detail::strictUpperIndex(j, k)]);
          sum += ui * static_cast<double>(carried.diagonal[k]) * uj;
        }
        return sum;
      };
      const double w00 = covariance(kLoc0, kLoc0) + static_cast<double>(measurement.varLoc0);
      const double w11 = covariance(kLoc1, kLoc1) + static_cast<double>(measurement.varLoc1);
      const double w01 = covariance(kLoc0, kLoc1) + static_cast<double>(measurement.covLoc01);
      factorDet = (w00 * w11) - (w01 * w01);
    }
    covarianceFromFactors(carried, predicted.cov);
    if (beamLine != nullptr && !updateWithBeamLine(propagator, *beamLine, geometry.planeOf(measurement), mTruth->rows[row].layer, pairPointAt(frame.measurements(), row), predicted).ok) {
      return {};
    }
    const auto [chi2, det] = innovation(predicted, measurement);
    if (!(det > 0.)) {
      return {};
    }
    return {.valid = true, .chi2 = chi2, .sqrtDet = std::sqrt(det), .floatDet = floatDet, .factorDet = factorDet};
  }

  template <typename Observation, typename Material>
  void decisionStudy(const Observation& iteration, const std::vector<Replay>& replayed, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, const std::unordered_map<std::uint64_t, std::uint32_t>& goodTrack)
  {
    const auto gate = static_cast<double>(parameters.clusterChi2Max);
    const auto view = iteration.tracks->view();
    std::vector<ClutterWalks> walks;
    for (std::size_t element = 0; element < replayed.size(); ++element) {
      // The kept tracks with one wrong row, their road one particle's.
      const auto track = iteration.trackOffset + static_cast<std::uint32_t>(element);
      const auto& replay = replayed[element];
      if (replay.roadCount == 0U || track >= iteration.tracks->size() || mTruth->rows[replay.road.index[0]].count == 0U) {
        continue;
      }
      const auto particle = mTruth->rows[replay.road.index[0]].particles[0];
      std::size_t wrong = 0U;
      for (std::uint8_t i = 0; i < view.count[track]; ++i) {
        wrong += mTruth->rows[view.clusters[track].index[i]].carries(particle) ? 0U : 1U;
      }
      if (wrong == 1U) {
        walks.push_back({
          .outward = replay.outward.empty() ? nullptr : &replay.outward.front(),
          .inward = &replay,
          .road = replay.road,
          .roadCount = replay.roadCount,
          .population = "kept, one wrong row",
        });
      }
    }
    for (const auto& dropped : mDroppedWalks) {
      if (dropped.walked) {
        walks.push_back({.outward = &dropped.outward, .inward = &dropped.inward, .road = dropped.road, .roadCount = dropped.roadCount, .population = "dropped"});
      }
    }
    std::mutex merge;
    o2::gtrk::detail::parallelFor(*mArena, 0U, walks.size(), 16, [&](std::size_t a, std::size_t b) {
      std::map<std::string, DecisionTally> local;
      for (std::size_t element = a; element < b; ++element) {
        const auto& followed = walks[element];
        const std::span<const std::uint32_t> road{followed.road.index, followed.roadCount};
        const auto particle = mTruth->rows[road[0]].particles[0];
        if (!std::ranges::all_of(road, [&](std::uint32_t row) { return mTruth->rows[row].carries(particle); })) {
          continue;
        }
        const auto index = mIndex.find(particle);
        if (index == mIndex.end()) {
          continue;
        }
        const auto& own = mRowsOf.at(particle);
        const bool timed = !mLegs.empty() && mTimed[index->second] != 0U;
        const auto positionOf = [&](std::uint32_t row) -> std::size_t { return static_cast<std::size_t>(std::ranges::find(own, row) - own.begin()); };
        const auto legAt = [&](std::size_t k) -> int { return !timed || k >= own.size() ? -1 : static_cast<int>(mLegs[index->second][k]); };
        // The first wrong row: the outward walk first, as the follow runs.
        const Step* at = nullptr;
        PropagationDirection direction = PropagationDirection::kForward;
        std::vector<std::uint32_t> outwardTaken;
        std::vector<std::uint32_t> inwardTaken;
        const auto scan = [&](const Walk* walk, PropagationDirection dir, std::vector<std::uint32_t>& taken) {
          if (walk == nullptr || at != nullptr) {
            return;
          }
          for (const auto& step : walk->steps) {
            if (step.choice == kNoRow) {
              continue;
            }
            if (!mTruth->rows[step.choice].carries(particle)) {
              at = &step;
              direction = dir;
              return;
            }
            taken.push_back(step.choice);
          }
        };
        scan(followed.outward, PropagationDirection::kForward, outwardTaken);
        scan(followed.inward, PropagationDirection::kBackward, inwardTaken);
        if (at == nullptr) {
          continue;
        }
        const auto& step = *at;
        const LayerID layer = mTruth->rows[step.choice].layer;
        const auto target = std::ranges::find(step.targets, layer, &TargetSeen::layer);
        if (target == step.targets.end()) {
          continue;
        }
        const std::uint32_t last = direction == PropagationDirection::kForward ? (outwardTaken.empty() ? road.back() : outwardTaken.back()) : (inwardTaken.empty() ? road.front() : inwardTaken.back());
        const std::size_t lastAt = positionOf(last);
        const int coreLeg = legAt(lastAt);
        // The particle's own row on the layer and leg, beyond the core's last.
        std::uint32_t mine = kNoRow;
        for (const auto k : positionsBeyond(lastAt, own.size(), direction)) {
          if (mTruth->rows[own[k]].layer == layer && legAt(k) == coreLeg) {
            mine = own[k];
            break;
          }
        }
        auto& tally = local[std::format("{} | {} | {} | {} | {}", followed.population, iteration.iteration, kPathPtNames[pathPtBin(particle)], direction == PropagationDirection::kForward ? "outward" : "inward", layerGroupOf(layer, parameters))];
        ++tally.events;
        ++tally.taken[lineageOf(particle, step.choice)];
        OwnThere there = kOwnNone;
        if (mine != kNoRow) {
          const auto seen = std::ranges::find(step.candidates, mine, &Candidate::row);
          there = seen == step.candidates.end() ? kOwnNotSeen : (seen->outcome == StepOutcome::kPasses ? kOwnAvailable : (seen->outcome == StepOutcome::kUsed ? kOwnHeld : kOwnRefused));
        }
        ++tally.own[there];
        // T and O from the step's state.
        struct Alternative {
          bool valid = false;
          double chi2 = 0.;
          double sqrtDet = 0.;
        };
        const auto alternative = [&](std::uint32_t row) {
          const auto one = row == kNoRow ? OneStep{} : oneStep(row, step.state, step.factors, direction, frame, geometry, material, propagator);
          return Alternative{.valid = one.valid, .chi2 = one.chi2, .sqrtDet = one.sqrtDet};
        };
        const auto taken = alternative(step.choice);
        const auto ownAlt = there == kOwnNone || there == kOwnNotSeen ? Alternative{} : alternative(mine);
        // rho: the window's rows beyond the gate (the clutter test's), and the layer's.
        const auto& window = target->window;
        const double rhoLocal = windowRho(parameters, step, layer, window, gate, taken.valid ? taken.sqrtDet : 0., direction, frame, geometry, material, propagator, [&](std::uint32_t row) { return mTruth->rows[row].carries(particle); }, std::numeric_limits<double>::quiet_NaN());
        const double rhoLayer = std::max(layerDensity(layer, window), 1e-12);
        const double p = reachChance(parameters, step.state, step.factors, layer, window, direction, geometry, propagator, step.choice);
        tally.reachSum += p;
        ++tally.reach[p < 0.25 ? 0U : (p < 0.5 ? 1U : (p < 0.9 ? 2U : 3U))];
        const double odds = std::log(p / (1. - p));
        const auto logG = [](const Alternative& one) { return (-0.5 * one.chi2) - std::log(one.sqrtDet) - std::log(2. * std::numbers::pi); };
        const auto choose = [&](double rho) {
          std::array<double, kDecisions> score{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), 0.};
          if (taken.valid) {
            score[kDecideTaken] = logG(taken) - std::log(rho) + odds;
          }
          if (ownAlt.valid) {
            score[kDecideOwn] = logG(ownAlt) - std::log(rho) + odds;
          }
          return static_cast<Decision>(std::ranges::max_element(score) - score.begin());
        };
        Decision stepChoice = kDecisions;
        if (std::isfinite(rhoLocal)) {
          stepChoice = choose(rhoLocal);
          ++tally.stepLocal[stepChoice];
        }
        ++tally.stepLayer[choose(rhoLayer)];
        // The chains, in the fit's order: the inward rows innermost first, the
        // road, the outward rows.
        const auto chainWith = [&](std::uint32_t extra) {
          std::vector<std::uint32_t> rows;
          std::vector<std::uint32_t> inward = inwardTaken;
          if (direction == PropagationDirection::kBackward && extra != kNoRow) {
            inward.push_back(extra);
          }
          for (auto k = inward.size(); k-- > 0U;) {
            rows.push_back(inward[k]);
          }
          rows.insert(rows.end(), road.begin(), road.end());
          rows.insert(rows.end(), outwardTaken.begin(), outwardTaken.end());
          if (direction == PropagationDirection::kForward && extra != kNoRow) {
            rows.push_back(extra);
          }
          return rows;
        };
        const auto accepted = [&](const std::vector<std::uint32_t>& rows) {
          const auto fit = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = rows.data()}, static_cast<int>(rows.size()), kinkBudget(parameters));
          return follow::passesRoadChi2(parameters, fit);
        };
        const auto withTaken = chainWith(step.choice);
        const auto prefix = chainWith(kNoRow);
        tally.refitTaken += accepted(withTaken) ? 1U : 0U;
        std::vector<std::uint32_t> withOwn;
        if (ownAlt.valid) {
          withOwn = chainWith(mine);
          ++tally.ownChains;
          tally.refitOwn += accepted(withOwn) ? 1U : 0U;
        }
        const bool prefixOk = accepted(prefix);
        tally.refitPrefix += prefixOk ? 1U : 0U;
        tally.prefixLong += prefixOk && prefix.size() >= parameters.minClustersAfterFollow ? 1U : 0U;
        tally.particleFound += goodTrack.contains(particle) ? 1U : 0U;
        if (mTruth->rows[step.choice].count != 0U) {
          const auto stranger = mTruth->rows[step.choice].particles[0];
          tally.strangerFound += goodTrack.contains(stranger) ? 1U : 0U;
          tally.strangerParticles += mIndex.contains(stranger) ? 1U : 0U;
        }
        // The chains' evidence, the last row's clutter and the reach in.
        if (std::isfinite(rhoLocal)) {
          std::array<std::array<double, kDecisions>, 3> score{};
          bool usable = true;
          // The running sums and the two chains' evidences at the seed's own
          // covariance (s = 0), for the prefix's inherited evidence below.
          std::vector<double> running;
          std::vector<double> prefixRunning;
          double withTakenEvidence = std::numeric_limits<double>::quiet_NaN();
          double alone = std::numeric_limits<double>::quiet_NaN();
          for (std::size_t s = 0; s < 3U; ++s) {
            const double scale = s == 0U ? 1. : (s == 1U ? 10. : 100.);
            const double takenEvidence = chainEvidence(withTaken, scale, frame, geometry, material, propagator, s == 0U ? &running : nullptr);
            const double stopEvidence = chainEvidence(prefix, scale, frame, geometry, material, propagator, s == 0U ? &prefixRunning : nullptr);
            if (s == 0U) {
              withTakenEvidence = takenEvidence;
              alone = stopEvidence;
            }
            score[s][kDecideTaken] = takenEvidence - std::log(rhoLocal) + std::log(p);
            score[s][kDecideOwn] = ownAlt.valid ? chainEvidence(withOwn, scale, frame, geometry, material, propagator) - std::log(rhoLocal) + std::log(p) : -std::numeric_limits<double>::infinity();
            score[s][kDecideStop] = stopEvidence + std::log(1. - p);
            usable = usable && std::isfinite(score[s][kDecideTaken]) && std::isfinite(score[s][kDecideStop]) && (!ownAlt.valid || std::isfinite(score[s][kDecideOwn]));
          }
          if (usable) {
            ++tally.chains;
            const auto best = [&](std::size_t s) { return static_cast<Decision>(std::ranges::max_element(score[s]) - score[s].begin()); };
            ++tally.chainChoice[best(0)];
            tally.chainStable[0] += best(1) == best(0) ? 1U : 0U;
            tally.chainStable[1] += best(2) == best(0) ? 1U : 0U;
            tally.chainAgreesStep += best(0) == stepChoice ? 1U : 0U;
            // The seed's own dependence on the chain: the prefix's evidence
            // fitted alone, with its own seed's terms, against its rows'
            // running sum inside the chain with T plus that chain's seed terms.
            if (direction == PropagationDirection::kForward && !prefixRunning.empty() && running.size() > prefixRunning.size() && std::isfinite(alone) && std::isfinite(withTakenEvidence)) {
              const double inherited = running[prefixRunning.size() - 1U] + (withTakenEvidence - running.back());
              const double gap = std::abs(alone - inherited);
              ++tally.prefixInherited[gap < 0.1 ? 0U : (gap < 1. ? 1U : 2U)];
            }
          }
        }
      }
      const std::scoped_lock lock(merge);
      for (const auto& [key, one] : local) {
        auto& sum = mDecisions[key];
        sum.events += one.events;
        sum.reachSum += one.reachSum;
        sum.refitTaken += one.refitTaken;
        sum.refitOwn += one.refitOwn;
        sum.ownChains += one.ownChains;
        sum.refitPrefix += one.refitPrefix;
        sum.prefixLong += one.prefixLong;
        sum.particleFound += one.particleFound;
        sum.strangerFound += one.strangerFound;
        sum.strangerParticles += one.strangerParticles;
        sum.chainAgreesStep += one.chainAgreesStep;
        sum.chains += one.chains;
        for (std::size_t i = 0; i < kOwnTheres; ++i) {
          sum.own[i] += one.own[i];
        }
        for (std::size_t i = 0; i < kLineages; ++i) {
          sum.taken[i] += one.taken[i];
        }
        for (std::size_t i = 0; i < kDecisions; ++i) {
          sum.stepLocal[i] += one.stepLocal[i];
          sum.stepLayer[i] += one.stepLayer[i];
          sum.chainChoice[i] += one.chainChoice[i];
        }
        for (std::size_t i = 0; i < 4U; ++i) {
          sum.reach[i] += one.reach[i];
        }
        for (std::size_t i = 0; i < 2U; ++i) {
          sum.chainStable[i] += one.chainStable[i];
        }
        for (std::size_t i = 0; i < 3U; ++i) {
          sum.prefixInherited[i] += one.prefixInherited[i];
        }
      }
    });
  }

  /// `decisionStudy`'s tables.
  void printDecisions() const
  {
    if (mDecisions.empty()) {
      return;
    }
    std::string text =
      "  the first wrong rows' alternatives: per population | iteration | pT | walk | layers of the wrong row:\n"
      "    events | the own row on that layer and leg: passing / held by a track / refused / not in the window / none | the taken row's lineage | "
      "the reach chance p: < 0.25 / < 0.5 / < 0.9 / >= 0.9 (mean) | the one-step evidence log(g/rho) + log(p/(1-p)) against stopping chooses taken / own / "
      "stop, with the window's rho; with the layer's | refits accepted: with the taken row, with the own row (of the chains that have one), the prefix "
      "(and as long as the length cut) | at the iteration's end: the particle found, the taken row's particle found (of them class particles) | "
      "the chains' evidence chooses taken / own / stop (of the chains computed), the same choice with the seed's covariance x 10 / x 100, agreeing with "
      "the one-step choice; the prefix alone against its rows' running sum in the chain with the taken row: |difference| < 0.1 / < 1 / more";
    for (const auto& [key, t] : mDecisions) {
      text += std::format("\n    {}: {} |", key, t.events);
      for (std::size_t i = 0; i < kOwnTheres; ++i) {
        text += std::format(" {}{}", t.own[i], i + 1U < kOwnTheres ? " /" : " |");
      }
      for (std::size_t i = 0; i < kLineages; ++i) {
        text += t.taken[i] != 0U ? std::format(" {} {};", kLineageNames[i], t.taken[i]) : std::string{};
      }
      text += std::format(
        " | {} / {} / {} / {} ({:.2f}) | {} / {} / {}; {} / {} / {} | {} / {} (of {}) / {} ({}) | {} / {} ({}) | {} / {} / {} (of {}), {} / {}, {} | "
        "{} / {} / {}",
        t.reach[0], t.reach[1], t.reach[2], t.reach[3], t.events != 0U ? t.reachSum / static_cast<double>(t.events) : 0., t.stepLocal[0], t.stepLocal[1], t.stepLocal[2], t.stepLayer[0], t.stepLayer[1], t.stepLayer[2], t.refitTaken, t.refitOwn, t.ownChains, t.refitPrefix, t.prefixLong, t.particleFound, t.strangerFound, t.strangerParticles, t.chainChoice[0], t.chainChoice[1], t.chainChoice[2], t.chains, t.chainStable[0], t.chainStable[1], t.chainAgreesStep, t.prefixInherited[0], t.prefixInherited[1], t.prefixInherited[2]);
    }
    std::cout << text << '\n';
  }

  /// The blockers: the class's particles not found after the iteration, by their seeding stage in it, and the rows
  /// of theirs held by tracks: how many, by how many distinct tracks (more
  /// than one: removing the first holder exposes another), fakes or another
  /// particle's good tracks, and the holders' iterations.
  struct BlockerTally {
    std::size_t particles = 0U;
    std::size_t blocked = 0U;             ///< with a row of theirs held
    std::array<std::size_t, 3> holders{}; ///< 1 / 2 / 3 or more distinct holders
    std::size_t heldRows = 0U;
    std::size_t byFake = 0U; ///< held rows whose holder is a fake
    std::size_t byGood = 0U; ///< held rows whose holder is another particle's good track
    std::size_t byThis = 0U; ///< held rows whose holder is the iteration's own
    /// Held rows whose holder is the particle's own track made fake (its
    /// majority particle, below the purity a good track needs), and the
    /// particles whose only holder is that.
    std::size_t byOwnFake = 0U;
    std::size_t onlyOwnFake = 0U;
  };
  template <typename Observation>
  void blockerStudy(const Observation& iteration, const std::vector<std::uint32_t>& holder, const std::vector<std::uint64_t>& particles)
  {
    const auto view = iteration.tracks->view();
    auto& tallies = mBlockers;
    tallies.clear(); // the last iteration's stand: the end
    const auto* stage = iteration.iteration < mSeedStage.size() && !mSeedStage[iteration.iteration].empty() ? &mSeedStage[iteration.iteration] : nullptr;
    for (std::uint32_t index = 0; index < mFindable.size(); ++index) {
      if (mFate[index] != kNotFound) {
        continue;
      }
      const auto particle = mFindable[index];
      auto& tally = tallies[std::format("{} | {}", kPathPtNames[pathPtBin(particle)], stage != nullptr ? kSeedStageNames[(*stage)[index]] : "no stage")];
      ++tally.particles;
      std::vector<std::uint32_t> distinct;
      bool otherHolder = false;
      for (const auto row : mRowsOf.at(particle)) {
        const auto track = holder[row];
        if (track == kNoTrack) {
          continue;
        }
        ++tally.heldRows;
        const std::span<const std::uint32_t> rows{view.clusters[track].index, view.count[track]};
        const bool ownFake = particles[track] == study::kNoTrackParticle && study::majorityParticle(mTruth->rows, rows).particle == particle;
        tally.byOwnFake += ownFake ? 1U : 0U;
        otherHolder = otherHolder || !ownFake;
        tally.byFake += particles[track] == study::kNoTrackParticle ? 1U : 0U;
        tally.byGood += particles[track] != study::kNoTrackParticle ? 1U : 0U;
        tally.byThis += track >= iteration.trackOffset ? 1U : 0U;
        if (std::ranges::find(distinct, track) == distinct.end()) {
          distinct.push_back(track);
        }
      }
      if (!distinct.empty()) {
        ++tally.blocked;
        tally.onlyOwnFake += otherHolder ? 0U : 1U;
        ++tally.holders[std::min<std::size_t>(distinct.size(), 3U) - 1U];
      }
    }
  }
  void printBlockers() const
  {
    if (mBlockers.empty()) {
      return;
    }
    std::string text =
      "  the class's particles not found at the end, by pT and their seeding stage in the last iteration: "
      "particles | with a row held by a track | of them by 1 / 2 / 3+ distinct holders | held rows: by a fake / by another particle's good "
      "track / by the last iteration's own tracks; by its own track made fake (the majority its, the purity short), and the particles held by that alone";
    for (const auto& [key, t] : mBlockers) {
      text += std::format("\n    {}: {} | {} | {} / {} / {} | {}: {} / {} / {}; {}, {}", key, t.particles, t.blocked, t.holders[0], t.holders[1], t.holders[2], t.heldRows, t.byFake, t.byGood, t.byThis, t.byOwnFake, t.onlyOwnFake);
    }
    std::cout << text << '\n';
  }

  /// `inwardStudy`'s tables.
  void printInward() const
  {
    if (mInwardCount[0] + mInwardCount[1] + mInwardCount[2] == 0U) {
      return;
    }
    const auto median = [](std::vector<double> values) {
      std::erase_if(values, [](double v) { return !std::isfinite(v); });
      return values.empty() ? std::numeric_limits<double>::quiet_NaN() : study::median(values);
    };
    const auto width = [&](const std::vector<double>& values) {
      const double centre = median(values);
      std::vector<double> deviations;
      deviations.reserve(values.size());
      for (const auto v : values) {
        deviations.push_back(std::abs(v - centre));
      }
      return 1.4826 * median(deviations);
    };
    std::string text = std::format(
      "  the refit's inward pass against the truth: passes replayed {} / {} / {} ({}, {}, {}); not reproduced {} / {} / {}\n"
      "    by state (0: the seed, the outward end state; k: after k updates inward), per kind: n | median p / p_true | % off by > 20 % | q/p pull median / "
      "width | median |dphi| [mrad] | phi pull width | median sigma(q/p) / (q/p)",
      mInwardCount[0], mInwardCount[1], mInwardCount[2], kInwardKindNames[0], kInwardKindNames[1], kInwardKindNames[2], mInwardMismatch[0], mInwardMismatch[1], mInwardMismatch[2]);
    for (std::size_t kind = 0; kind < kInwardKinds; ++kind) {
      for (std::size_t k = 0; k < kInwardStepSlots; ++k) {
        const auto& by = mInwardSteps[kind][k];
        if (by.ratio.empty()) {
          continue;
        }
        std::size_t off = 0U;
        std::vector<double> absPhi;
        for (std::size_t i = 0; i < by.ratio.size(); ++i) {
          off += std::abs(by.ratio[i] - 1.) > 0.2 ? 1U : 0U;
          absPhi.push_back(std::abs(by.dphi[i]));
        }
        text += std::format("\n    {} | {}{}: {} | {:.3f} | {:.1f} | {:+.2f} / {:.2f} | {:.2f} | {:.2f} | {:.3f}", kInwardKindNames[kind], k, k + 1U == kInwardStepSlots ? "+" : "", by.ratio.size(), median(by.ratio), 100. * static_cast<double>(off) / static_cast<double>(by.ratio.size()), median(by.pullQOverP), width(by.pullQOverP), median(absPhi), width(by.pullPhi), median(by.relSigma));
      }
    }
    text += std::format("\n    the seed's q/p sign against the MC charge (convention: q/p = {} charge / p): agree / disagree:", kChargeSign > 0. ? "+" : "-");
    for (std::size_t kind = 0; kind < kInwardKinds; ++kind) {
      text += std::format(" {} {} / {};", kInwardKindNames[kind], mChargeAgree[kind][0], mChargeAgree[kind][1]);
    }
    text += "\n    each successful inward step onto a barrel plane: the crossing the solve takes against the MC crossing (nearer / the other nearer / none taken)";
    for (std::size_t kind = 0; kind < kInwardKinds; ++kind) {
      for (std::size_t k = 1; k < kInwardStepSlots; ++k) {
        const auto& r = mStepRoots[kind][k];
        if (r[1] + r[2] + r[3] == 0U) {
          continue;
        }
        text += std::format("\n      {} | {}{}: {} / {} / {}", kInwardKindNames[kind], k, k + 1U == kInwardStepSlots ? "+" : "", r[1], r[2], r[3]);
      }
    }
    text += "\n    the same steps by the crossings' place (inside the sensor's bounds, turn, ahead or behind along the inward pass):";
    for (const auto& [key, count] : mStepRootDetail) {
      text += std::format("\n      {}: {}", key, count);
    }
    text += "\n    the failing own-row passes from other seeds (completion only, not acceptance): passes / fails at the same step / later / earlier (or no MC row)";
    for (std::size_t variant = 0; variant < kSeedNames.size(); ++variant) {
      text += std::format("\n      {}: {} / {} / {} / {}", kSeedNames[variant], mSeedOutcomes[variant][0], mSeedOutcomes[variant][1], mSeedOutcomes[variant][2], mSeedOutcomes[variant][3]);
    }
    text += "\n    the failing step's two crossings of the plane's line, from the filtered state it started from (status 1 no intersection, 2 wrong direction):";
    for (const auto& [key, count] : mRoots) {
      text += std::format("\n      {}: {}", key, count);
    }
    text +=
      "\n    the same rows refitted under each crossing rule: failed outward / failed inward / chi2 over the cut / accepted | against the "
      "default from the same seed: accepted by both, not bit for bit / the default's only / this one's only | the accepted's pulls on the innermost "
      "row, median / width: q/p, phi, theta, loc0, loc1";
    for (std::size_t kind = 0; kind < kInwardKinds; ++kind) {
      for (std::size_t seed = 0; seed < kRuleSeedNames.size(); ++seed) {
        for (std::size_t rule = 0; rule < kCrossingRules.size(); ++rule) {
          const auto& t = mRuleTallies[kind][rule][seed];
          text += std::format("\n      {} | {} | {}: {} / {} / {} / {} | {} / {} / {} |", kInwardKindNames[kind], kRuleSeedNames[seed], kCrossingRuleNames[rule], t.outcome[kRuleOutward], t.outcome[kRuleInward], t.outcome[kRuleChi2], t.outcome[kRuleAccepted], t.against[0], t.against[1], t.against[2]);
          for (const auto* pulls : {&t.pullQOverP, &t.pullPhi, &t.pullTheta, &t.pullLoc0, &t.pullLoc1}) {
            text += std::format(" {:+.2f} / {:.2f};", median(*pulls), width(*pulls));
          }
        }
      }
    }
    std::cout << text << '\n';
  }

  /// `FailedStep` of the refit \a fit of \a rows, failed.
  template <typename Material>
  [[nodiscard]] static FailedStep replayFailure(const std::vector<std::uint32_t>& rows, const TrackFitResult& fit, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    FailedStep out;
    const int n = static_cast<int>(rows.size());
    const int failed = fit.failedMeasurement;
    if ((fit.failedPass != FitPass::kOutward && fit.failedPass != FitPass::kInward) || failed < 0 || failed >= n) {
      return out;
    }
    const RowChain<Det> chain{.frame = &frame, .rows = rows.data()};
    // The failing plane's transverse normal and distance; a circle's farthest
    // point along that normal, its centre's projection plus its radius.
    const auto& target = geometry.planeOf(chain[failed]);
    const auto nx = static_cast<double>(target.normal[0]);
    const auto ny = static_cast<double>(target.normal[1]);
    const double transverse = std::hypot(nx, ny);
    out.barrel = transverse > 0.9;
    const double ux = transverse > 0. ? nx / transverse : 0.;
    const double uy = transverse > 0. ? ny / transverse : 0.;
    const double distance = std::abs((static_cast<double>(target.origin[0]) * ux) + (static_cast<double>(target.origin[1]) * uy));
    const double sign = ((static_cast<double>(target.origin[0]) * ux) + (static_cast<double>(target.origin[1]) * uy)) < 0. ? -1. : 1.;
    out.planeDistance = static_cast<float>(distance);
    const auto margins = [&](const BoundTrackState& one) {
      Scalar global[3]{};
      surfaceToGlobal(geometry.planeOf(one), one.loc0(), one.loc1(), global);
      const double curvature = propagator.curvature(one.qOverP(), std::sin(one.theta()));
      const double cx = static_cast<double>(global[0]) - (std::sin(static_cast<double>(one.phi())) / curvature);
      const double cy = static_cast<double>(global[1]) + (std::cos(static_cast<double>(one.phi())) / curvature);
      const double centre = sign * ((cx * ux) + (cy * uy));
      const double radius = 1. / std::abs(curvature);
      out.farMargin = static_cast<float>(centre + radius - distance);
      out.nearMargin = static_cast<float>(distance - (centre - radius));
    };
    const auto chainSeed = makeChainSeed(propagator, geometry, chain, n);
    const bool forward = fit.failedPass == FitPass::kOutward;
    BoundTrackState seed = chainSeed.state;
    if (!forward) {
      // `inwardPass`'s seed: the outward pass's end state.
      seed = fit.outwardState;
      setDiffuseCovariance(seed, chainSeed.qOverPSigma);
    }
    // `runPass` up to the failing measurement.
    const auto direction = forward ? PropagationDirection::kForward : PropagationDirection::kBackward;
    const int first = forward ? 0 : n - 1;
    const int increment = forward ? 1 : -1;
    const auto budget = o2::gtrk::detail::kinkBudgetFor(kinkBudget(parameters), n);
    BoundTrackState state = seed;
    CovarianceFactors factors;
    {
      const auto measurement = chain[first];
      out.startMiss = std::hypot(seed.loc0() - measurement.loc0, seed.loc1() - measurement.loc1);
    }
    o2::gtrk::detail::startOnMeasurement(state, factors, chain[first]);
    Scalar chi2 = 0.F;
    int kinks = 0;
    for (int step = 1; step < n; ++step) {
      const int index = first + (step * increment);
      const auto& stepProp = stepPropagator(propagator, geometry, chain, index - increment, index);
      if (index == failed) {
        out.step = step;
        BoundTrackState filtered = state;
        out.state = propagate(stepProp, filtered, geometry.planeOf(state), geometry.planeOf(chain[index]), direction).status;
        margins(state);
        out.replayed = true;
        return out;
      }
      if (o2::gtrk::detail::kinkStep(stepProp, geometry, material, chain[index], direction, state, factors, chi2, budget, step < n - 1, kinks) != PropagationStatus::kOk) {
        return out;
      }
    }
    return out;
  }

  /// One follow of the iteration again (see the file's comment): of the
  /// tracks from \a first to \a last, with their \a proposals and
  /// \a kinkSteps as it made them, on the mask \a used as it found it; each
  /// track from the refit of its road, or with \a fromSnapshot from the
  /// previous iteration's track whole (followEarlier).
  template <typename Observation, typename Proposals, typename Material>
  [[nodiscard]] std::vector<Replay> replayExtension(const Observation& iteration, std::uint32_t first, std::uint32_t last, const Proposals& proposals, std::span<const follow::KinkSteps> kinkSteps, std::vector<std::uint8_t>& used, bool fromSnapshot, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    const bool kinks = !kinkSteps.empty();
    const auto view = iteration.tracks->view();
    std::vector<Replay> replays(last - first);
    if (!parameters.follow || proposals.size() != replays.size()) {
      return replays;
    }
    // The stage as the iteration instantiated it (Follow's dispatch): the
    // choice among passing rows and the end decisions on the last row's
    // density change which rows a walk takes. FollowEnd stays off: its end
    // states are the follow's own.
    withFollowStage(iteration, parameters, frame, geometry, material, propagator, used.data(), first, [&](const auto& stage) {
      // Each track's walks on its own thread: the stage is read-only, the
      // replay is the track's.
      o2::gtrk::detail::parallelFor(*mArena, first, last, 256, [&](std::size_t a, std::size_t b) {
        for (auto track = static_cast<std::uint32_t>(a); track < b; ++track) {
          const auto& proposal = proposals[track - first];
          auto& replay = replays[track - first];
          // The road: the track's rows less the ones its follow added, in
          // the road's order, and its refit as the iteration had it -- or the
          // earlier track as it was, with its state.
          ClusterList<kMaxClusters> road{};
          std::uint8_t count = 0U;
          BoundTrackState start;
          BoundTrackState outer;
          if (fromSnapshot) {
            road = mSnapshot[track].clusters;
            count = mSnapshot[track].count;
            start = mSnapshot[track].state;
            outer = mSnapshot[track].outwardState;
          } else {
            for (std::uint8_t i = 0; i < view.count[track]; ++i) {
              const auto row = view.clusters[track].index[i];
              const bool added = std::find(proposal.inward, proposal.inward + proposal.inwardCount, row) != proposal.inward + proposal.inwardCount || std::find(proposal.outward, proposal.outward + proposal.outwardCount, row) != proposal.outward + proposal.outwardCount;
              if (!added) {
                road.index[count++] = row;
              }
            }
            for (auto i = static_cast<std::size_t>(count); i < kMaxClusters; ++i) {
              road.index[i] = kNoRow;
            }
            // The states the follow started from where the tracker kept them,
            // else the road's refit again.
            const auto element = track - first;
            if (element < iteration.followStart.size() && element < iteration.followStartInner.size()) {
              start = iteration.followStartInner[element];
              outer = iteration.followStart[element];
            } else {
              const auto fit = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = road.index}, static_cast<int>(count), kinkBudget(parameters));
              if (fit.status != PropagationStatus::kOk) {
                continue;
              }
              start = fit.state;
              outer = fit.outwardState;
            }
          }
          Walk outward;
          const auto walked = recordFollow(stage, start, outer, road, count, replay, outward);
          replay.road = road;
          replay.roadCount = count;
          replay.agrees = walked.inwardCount >= proposal.inwardCount && std::equal(proposal.inward, proposal.inward + proposal.inwardCount, walked.inward);
          outward.road = road;
          outward.roadCount = count;
          outward.bidKept = proposal.outwardCount;
          outward.kept = proposal.outwardCount;
          outward.agrees = walked.outwardCount >= proposal.outwardCount && std::equal(proposal.outward, proposal.outward + proposal.outwardCount, walked.outward);
          replay.outward.push_back(std::move(outward));
        }
      });
    });
    for (const auto& replay : replays) {
      if (!replay.outward.empty()) {
        ++mOutReplayed;
        mOutAgree += replay.outward.front().agrees ? 1U : 0U;
      }
    }
    // The bidding: each row to the proposing track first in output order;
    // each track keeps its proposal up to its first row lost (`Bid`,
    // `Resolve`), then up to the last kink step whose witness it keeps
    // (`KeepWitnesses`). Outward proposals bid with what the iteration kept
    // of them.
    std::vector<std::uint32_t> slots(mTruth->rows.size(), 0U);
    for (std::uint32_t track = first; track < last; ++track) {
      const std::uint32_t priority = last - track;
      for (const auto row : replays[track - first].raw) {
        slots[row] = std::max(slots[row], priority);
      }
      // The outward proposal as walked (the replayed walk, which the
      // iteration's agrees with), else what the iteration kept of it.
      const auto& proposal = proposals[track - first];
      if (!replays[track - first].outward.empty()) {
        for (const auto row : replays[track - first].outward.front().raw) {
          slots[row] = std::max(slots[row], priority);
        }
      } else {
        for (std::uint8_t i = 0; i < proposal.outwardCount; ++i) {
          slots[proposal.outward[i]] = std::max(slots[proposal.outward[i]], priority);
        }
      }
    }
    // Each track's bid outcome and refit on its own thread; the fidelity
    // counted after.
    std::vector<std::uint8_t> faithful(last - first, 0U);
    o2::gtrk::detail::parallelFor(*mArena, first, last, 256, [&](std::size_t a, std::size_t b) {
      for (auto track = static_cast<std::uint32_t>(a); track < b; ++track) {
        auto& replay = replays[track - first];
        const std::uint32_t priority = last - track;
        replay.bidKept = static_cast<std::uint32_t>(replay.raw.size());
        for (std::uint32_t i = 0; i < replay.raw.size(); ++i) {
          if (slots[replay.raw[i]] != priority) {
            replay.bidKept = i;
            break;
          }
        }
        replay.kept = kinks ? follow::keptWithWitnesses(static_cast<std::uint8_t>(replay.bidKept), replay.witnessed) : replay.bidKept;
        // The outward walk's bid the same way.
        const auto& proposal = proposals[track - first];
        for (auto& outward : replay.outward) {
          outward.bidKept = static_cast<std::uint32_t>(outward.raw.size());
          for (std::uint32_t i = 0; i < outward.raw.size(); ++i) {
            if (slots[outward.raw[i]] != priority) {
              outward.bidKept = i;
              break;
            }
          }
          outward.kept = kinks ? follow::keptWithWitnesses(static_cast<std::uint8_t>(outward.bidKept), outward.witnessed) : outward.bidKept;
        }
        const std::uint32_t outwardKept = replay.outward.empty() ? proposal.outwardCount : replay.outward.front().kept;
        // The extended refit's acceptance (`follow::Refit`): the inward rows
        // kept, the road, the outward rows kept.
        if (replay.kept + outwardKept != 0U) {
          std::vector<std::uint32_t> combined;
          for (auto i = replay.kept; i-- > 0U;) {
            combined.push_back(replay.raw[i]);
          }
          if (replay.outward.empty()) {
            for (std::uint8_t i = 0; i < view.count[track]; ++i) {
              const auto row = view.clusters[track].index[i];
              if (std::find(proposal.inward, proposal.inward + proposal.inwardCount, row) == proposal.inward + proposal.inwardCount) {
                combined.push_back(row);
              }
            }
          } else {
            combined.insert(combined.end(), replay.road.index, replay.road.index + replay.roadCount);
            const auto& outward = replay.outward.front().raw;
            combined.insert(combined.end(), outward.begin(), outward.begin() + outwardKept);
          }
          replay.refitRows = combined;
          replay.refitCount = static_cast<std::uint8_t>(std::min(combined.size(), kMaxClusters));
          for (std::size_t i = 0; i < replay.refitCount; ++i) {
            replay.refitLayer[i] = mTruth->rows[combined[i]].layer;
            replay.refitIncrement[i] = 0;
          }
          const auto increments = [&replay](PropagationDirection direction, int index, Scalar chi2) {
            if (direction == PropagationDirection::kBackward && index >= 0 && index < replay.refitCount) {
              replay.refitIncrement[static_cast<std::size_t>(index)] = chi2;
            }
          };
          const auto fit = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = combined.data()}, static_cast<int>(combined.size()), kinkBudget(parameters), increments);
          if (fit.status != PropagationStatus::kOk) {
            replay.refit = kRefitFailedFit;
            replay.refitStatus = fit.status;
            replay.refitFailedPass = fit.failedPass;
            replay.refitFailedLayer = fit.failedMeasurement >= 0 && static_cast<std::size_t>(fit.failedMeasurement) < combined.size() ? mTruth->rows[combined[static_cast<std::size_t>(fit.failedMeasurement)]].layer : kInvalidLayerID;
            replay.failedStep = replayFailure(combined, fit, parameters, frame, geometry, material, propagator);
            for (int kept = static_cast<int>(outwardKept) - 1; kept >= 0 && !replay.outward.empty(); --kept) {
              const auto size = combined.size() - (outwardKept - static_cast<std::size_t>(kept));
              const auto trimmed = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = combined.data()}, static_cast<int>(size), kinkBudget(parameters));
              if (trimmed.status == PropagationStatus::kOk) {
                replay.refitTrimmedKept = kept;
                replay.refitTrimmedChi2 = follow::passesRoadChi2(parameters, trimmed);
                break;
              }
            }
          } else if (!follow::passesRoadChi2(parameters, fit)) {
            replay.refit = kRefitAboveChi2;
          }
          if (fit.status == PropagationStatus::kOk) {
            replay.refitChi2 = fit.chi2;
            replay.refitChi2Max = parameters.roadChi2MaxByNdf[fit.ndf];
          }
        }
        for (auto& outward : replay.outward) {
          outward.refit = replay.refit;
          outward.refitChi2 = replay.refitChi2;
          outward.refitChi2Max = replay.refitChi2Max;
          outward.refitStatus = replay.refitStatus;
          outward.refitFailedPass = replay.refitFailedPass;
          outward.refitFailedLayer = replay.refitFailedLayer;
          outward.refitTrimmedKept = replay.refitTrimmedKept;
          outward.refitTrimmedChi2 = replay.refitTrimmedChi2;
          outward.failedStep = replay.failedStep;
          outward.refitIncrement = replay.refitIncrement;
          outward.refitLayer = replay.refitLayer;
          outward.refitCount = replay.refitCount;
        }
        // Fidelity: the rows the replay says the iteration kept are the ones it
        // kept, and its kink steps' witnesses are the ones it recorded.
        const std::uint32_t replayKept = replay.refit == kRefitKept ? replay.kept : 0U;
        bool agrees = replayKept == proposal.inwardCount && (!kinks || replay.witnessed == kinkSteps[track - first].inwardWitnessed);
        for (std::uint32_t i = 0; agrees && i < replayKept; ++i) {
          agrees = replay.raw[i] == proposal.inward[i];
        }
        // A disagreement where the replay refuses the refit and the
        // iteration kept rows inward: what the stage's prefix second chance
        // (`iterationFollowPrefix`), not replayed, would explain.
        faithful[track - first] = agrees ? kFaithful : (replay.refit != kRefitKept && proposal.inwardCount != 0U ? kRefusedKept : kUnfaithful);
      }
    });
    if (mFidelity.size() <= iteration.iteration) {
      mFidelity.resize(iteration.iteration + 1U);
    }
    auto& fidelity = mFidelity[iteration.iteration][fromSnapshot ? 1U : 0U];
    for (const auto one : faithful) {
      ++mReplayed;
      ++fidelity.replayed;
      mAgree += one == kFaithful ? 1U : 0U;
      fidelity.agree += one == kFaithful ? 1U : 0U;
      fidelity.refusedKept += one == kRefusedKept ? 1U : 0U;
    }
    return replays;
  }

  /// The replay's fidelity per iteration, the iteration's own tracks and the
  /// earlier tracks it extends (`followEarlier`) apart.
  void printFidelity() const
  {
    std::string text;
    for (std::size_t it = 0; it < mFidelity.size(); ++it) {
      for (std::size_t kind = 0; kind < 2U; ++kind) {
        const auto& f = mFidelity[it][kind];
        if (f.replayed != 0U) {
          text += std::format("\n    iteration {}, {}: {} of {} agree; of the {} others, the replay refuses the refit where the iteration kept rows inward {}", it, kind == 0U ? "its own tracks" : "the earlier tracks", f.agree, f.replayed, f.replayed - f.agree, f.refusedKept);
        }
      }
    }
    if (!text.empty()) {
      std::cout << "  the replayed follow against the iteration's, per iteration (not replayed: the refit's prefix second chance, its trim of an "
                   "ambiguous last row, the crossing rule):"
                << text << '\n';
    }
  }

  /// The decision that left \a row, the particle's next row inward, out of
  /// the track: in the bidding or the refit if the replay proposed it, else
  /// at the step from \a outside, the layer just outside it -- or, where the
  /// step had the kink allowance, at its kink step.
  [[nodiscard]] Short stopReason(const Walk& replay, std::uint32_t row, LayerID outside, std::uint64_t particle, std::uint8_t& detail) const
  {
    if (row == kNoTrack) {
      return kNoStep;
    }
    const auto at = std::ranges::find(replay.raw, row);
    if (at != replay.raw.end()) {
      const auto position = static_cast<std::uint32_t>(at - replay.raw.begin());
      if (position >= replay.bidKept) {
        return kLostBid;
      }
      if (position >= replay.kept) {
        return kWitnessLost;
      }
      return replay.refit == kRefitFailedFit ? kRefitFailed : kRefitChi2;
    }
    const auto step = std::ranges::find(replay.steps, outside, &Step::layer);
    if (step == replay.steps.end()) {
      return kNoStep;
    }
    const auto seen = std::ranges::find(step->candidates, row, &Candidate::row);
    if (step->kinksLeft >= 0 && (seen == step->candidates.end() || seen->outcome != StepOutcome::kUsed)) {
      // No row passed the step: its kink step decided.
      if (step->kinksLeft == 0) {
        return kKinkBudget;
      }
      const auto scaled = std::ranges::find(step->scaled, row, &Candidate::row);
      if (scaled == step->scaled.end() || scaled->outcome != StepOutcome::kPasses) {
        detail = static_cast<std::uint8_t>(scaled == step->scaled.end() ? kNotVisited : static_cast<std::size_t>(scaled->outcome));
        return kKinkRefused;
      }
      bool confirmed = false;
      for (const auto& witness : step->witnesses) {
        confirmed = confirmed || (witness.row == row && witness.outcome == StepOutcome::kPasses);
      }
      if (confirmed) {
        return kKinkOther;
      }
      // What became of the particle's own next row, or of the beam line. The
      // next row is taken on the layer numbered two below \a outside: the
      // class's barrel layers numbered consecutively from the beam outward
      // and an inward walk (an outward one's witness is on the layer above;
      // a skipped layer or a disc is not supported).
      const auto next = outside >= 2U ? ownRowOn(particle, static_cast<LayerID>(outside - 2U)) : kNoTrack;
      std::size_t beam = 0U;
      std::size_t outcome = next == kNoTrack ? kNoOwnNext : kNotVisited;
      for (const auto& witness : step->witnesses) {
        if (witness.row == row && witness.next == kNoRow) {
          beam = 1U;
          outcome = static_cast<std::size_t>(witness.outcome);
        } else if (witness.row == row && next != kNoTrack && witness.next == next) {
          outcome = static_cast<std::size_t>(witness.outcome);
        }
      }
      detail = static_cast<std::uint8_t>((8U * beam) + outcome);
      return kKinkUnconfirmed;
    }
    if (seen == step->candidates.end()) {
      return kNotQueried;
    }
    switch (seen->outcome) {
      case StepOutcome::kUsed:
        return kUsed;
      case StepOutcome::kOutsideWindow:
        return kWindow;
      case StepOutcome::kHole:
        return kHole;
      case StepOutcome::kNoUpdate:
        return kNoUpdate;
      case StepOutcome::kChi2:
        return kChi2;
      case StepOutcome::kPasses:
        return step->lookahead ? kLookaheadOther : kOtherTaken;
    }
    return kNoStep;
  }

  /// What the step from \a outside took where \a row, the particle's own,
  /// passed it too (kOtherTaken, kLookaheadOther): the
  /// taken row's `TakenKind`, plus `kLowerChi2` where its chi2 was below the
  /// own row's. kNoTaken where the step does not apply.
  enum TakenKind : std::uint8_t { kNoTaken,
                                  kOwnOtherLayer,
                                  kOtherSameLayer,
                                  kOtherOtherLayer,
                                  kNoiseSameLayer,
                                  kNoiseOtherLayer,
                                  kOwnDropped,
                                  kOwnRefusedLast,
                                  kNothingTaken,
                                  kTakenKinds };
  static constexpr std::uint8_t kLowerChi2 = 16U;
  static_assert(kTakenKinds <= kLowerChi2);
  static constexpr std::array<const char*, kTakenKinds> kTakenNames{"", "its own row on another layer", "another particle's on its layer", "another particle's on another layer", "noise on its layer", "noise on another layer", "the own row itself, dropped after (the unconfirmed drop: the next step crossed a sensor without a row)", "the own row itself, refused after as the walk's last (EndDensity, BeamEnd)", "nothing"};
  [[nodiscard]] std::uint8_t takenKind(const Walk& replay, std::uint32_t row, LayerID outside, std::uint64_t particle) const
  {
    if (row == kNoTrack) {
      return kNoTaken;
    }
    const auto step = std::ranges::find(replay.steps, outside, &Step::layer);
    if (step == replay.steps.end()) {
      return kNoTaken;
    }
    if (step->choice == row) {
      return replay.refused == row ? kOwnRefusedLast : kOwnDropped;
    }
    if (step->choice == kNoRow) {
      return kNothingTaken;
    }
    const auto own = std::ranges::find(step->candidates, row, &Candidate::row);
    const auto taken = std::ranges::find(step->candidates, step->choice, &Candidate::row);
    if (own == step->candidates.end() || own->outcome != StepOutcome::kPasses) {
      return kNoTaken;
    }
    const auto& truth = mTruth->rows[step->choice];
    const bool sameLayer = truth.layer == mTruth->rows[row].layer;
    std::uint8_t kind = kOtherOtherLayer;
    if (truth.carries(particle)) {
      kind = kOwnOtherLayer;
    } else if (truth.count == 0U) {
      kind = sameLayer ? kNoiseSameLayer : kNoiseOtherLayer;
    } else if (sameLayer) {
      kind = kOtherSameLayer;
    }
    const bool lower = taken != step->candidates.end() && taken->chi2 < own->chi2;
    return static_cast<std::uint8_t>(kind + (lower ? kLowerChi2 : 0U));
  }

  /// A true kink over theta0 binned as the refit-chi2 report bins it: <= 3,
  /// 3-5, > 5, no hits (-1).
  [[nodiscard]] static std::uint8_t kinkBin(float kink) { return static_cast<std::uint8_t>(kink < 0.F ? 3U : (kink <= 3.F ? 0U : (kink <= 5.F ? 1U : 2U))); }

  /// A replayed walk's steps, for one track: each step's layer,
  /// the row it stepped from, its choice, and every candidate with its
  /// layer, outcome, chi2 and whether it carries \a particle; \a own the
  /// row the diagnosis asked about.
  void dumpWalk(const Walk& walk, std::uint64_t particle, std::uint32_t own) const
  {
    std::cout << std::format("  walk dump: particle {} own row {} (layer {}), dropped {}, end {}, raw", particle, own, static_cast<int>(mTruth->rows[own].layer), walk.dropped, static_cast<int>(walk.end));
    for (const auto row : walk.raw) {
      std::cout << ' ' << row;
    }
    std::cout << '\n';
    for (const auto& step : walk.steps) {
      std::cout << std::format("    step from layer {} row {}: choice {}{}{}, candidates:", static_cast<int>(step.layer), step.from, step.choice, step.kink ? " (kink)" : "", step.lookahead ? " (lookahead)" : "");
      for (const auto& candidate : step.candidates) {
        std::cout << std::format(" [{} L{} {} {:.2f}{}]", candidate.row, static_cast<int>(mTruth->rows[candidate.row].layer), kOutcomeNames[static_cast<std::size_t>(candidate.outcome)], static_cast<double>(candidate.chi2), mTruth->rows[candidate.row].carries(particle) ? " own" : "");
      }
      std::cout << '\n';
    }
  }

  /// For a row the unconfirmed drop removed (kOwnDropped): what the step from it did with the
  /// particle's own row there, a `StepOutcome`, or 6 where none of the
  /// particle's rows was among its candidates; and that row's chi2.
  [[nodiscard]] std::pair<std::uint8_t, Scalar> afterDrop(const Walk& replay, std::uint32_t row, std::uint64_t particle) const
  {
    const auto step = std::ranges::find(replay.steps, row, &Step::from);
    if (step == replay.steps.end()) {
      return {6U, 0.F};
    }
    for (const auto& candidate : step->candidates) {
      if (candidate.row != row && mTruth->rows[candidate.row].carries(particle)) {
        return {static_cast<std::uint8_t>(candidate.outcome), candidate.chi2};
      }
    }
    return {6U, 0.F};
  }
  /// The chi2 of \a row among the candidates of the step from \a outside; 0
  /// where it was not one.
  [[nodiscard]] Scalar candidateChi2(const Walk& replay, std::uint32_t row, LayerID outside) const
  {
    const auto step = std::ranges::find(replay.steps, outside, &Step::layer);
    if (step == replay.steps.end()) {
      return 0.F;
    }
    const auto seen = std::ranges::find(step->candidates, row, &Candidate::row);
    return seen == step->candidates.end() ? 0.F : seen->chi2;
  }

  /// Where \a row lies against the step's windows: inside the pair window
  /// (1), inside the covariance window at the pair cut's nSigma (2), inside
  /// it at the update gate's sqrt(chi2 max) (4) -- the covariance windows
  /// unclipped by the pair window. The window's sigma is the gate's,
  /// sqrt(clusterChi2Max), so a scaled copy of the parameters sets it.
  /// The step walks in \a direction. With \a offsets, the row's distance
  /// from the window's centre in azimuth and in the second coordinate, each
  /// over the gate window's half-width (above 1: outside it on that axis);
  /// the azimuth's signed by the track's charge; the azimuth half-width
  /// along r phi [um].
  template <typename Material>
  [[nodiscard]] std::size_t windowAudit(const Step& step, std::uint32_t row, LayerID target, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, PropagationDirection direction, std::array<float, 4>* offsets = nullptr) const
  {
    const auto nominal = predictLayerCrossing(propagator, parameters, geometry.planeOf(step.state), step.state, target, direction);
    if (nominal.status != PropagationStatus::kOk) {
      return 0U;
    }
    // The link from the inner of the two layers to the outer.
    const auto link = direction == PropagationDirection::kBackward ? linkOf<Det>(target, step.layer) : linkOf<Det>(step.layer, target);
    const auto point = frame.measurement(row);
    const Scalar second = secondCoordinate<Det>(target, point.r, point.z);
    // The window's centre as the follow has it: about the beam line.
    Scalar phiCentre = 0.F;
    Scalar secondCentre = 0.F;
    follow::windowCentre(parameters, target, nominal, phiCentre, secondCentre);
    const auto inside = [&](Scalar phiHalfWidth, Scalar secondHalfWidth) {
      const CandidateWindow window{
        .earliestPS = 0,
        .latestPS = 0,
        .afterPS = 0U,
        .phiCentre = phiCentre,
        .phiHalfWidth = phiHalfWidth,
        .secondCentre = secondCentre,
        .secondHalfWidth = secondHalfWidth,
      };
      return spatiallyCompatible(window, point.phi, second);
    };
    const auto covariance = [&](Scalar sigma, std::array<float, 4>* ratios = nullptr) {
      auto scaled = std::make_unique<IterationParameters<Det>>(parameters);
      scaled->clusterChi2Max = sigma * sigma;
      CandidateWindow window{
        .earliestPS = 0,
        .latestPS = 0,
        .afterPS = 0U,
        .phiCentre = phiCentre,
        .phiHalfWidth = std::numeric_limits<Scalar>::infinity(),
        .secondCentre = secondCentre,
        .secondHalfWidth = std::numeric_limits<Scalar>::infinity(),
      };
      follow::predictedWindow<false>(propagator, *scaled, geometry, material, step.state, step.factors, target, direction, nominal, window, 1.F);
      if (ratios != nullptr) {
        const Scalar dPhi = std::remainder(point.phi - phiCentre, 2.F * std::numbers::pi_v<Scalar>);
        (*ratios)[0] = static_cast<float>(std::abs(dPhi) / window.phiHalfWidth);
        (*ratios)[1] = static_cast<float>(std::abs(second - secondCentre) / window.secondHalfWidth);
        (*ratios)[2] = static_cast<float>((step.state.qOverP() < 0.F ? -dPhi : dPhi) / window.phiHalfWidth);
        (*ratios)[3] = static_cast<float>(window.phiHalfWidth * point.r * 1.e4F);
      }
      return inside(window.phiHalfWidth, window.secondHalfWidth);
    };
    std::size_t flags = 0U;
    flags |= link != kNoLink && inside(parameters.windows[link].phiHalfWidth, parameters.windows[link].secondHalfWidth) ? 1U : 0U;
    flags |= covariance(parameters.pairCut.nSigma) ? 2U : 0U;
    flags |= covariance(std::sqrt(parameters.clusterChi2Max), offsets) ? 4U : 0U;
    return flags;
  }

  /// A true road the road fit refuses: the fit repeated for its
  /// chi2 against the cut (a lower bound: the fit stops at the first step
  /// over it), and the road's barrel rows and its other rows fitted alone,
  /// each where they make a road; binned by the particle's side, species and
  /// true kink.
  template <typename Material>
  void describeRefusedRoad(const study::ParticleTruth& particle, float kink, const std::vector<std::uint32_t>& road, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    const auto fitOf = [&](const std::vector<std::uint32_t>& rows) {
      if (rows.size() < 4U) {
        return std::string{"n/a"};
      }
      const auto fit = fitRoad(parameters, propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = rows.data()}, static_cast<int>(rows.size()));
      return std::string{fit.status == RoadFitStatus::kAccepted ? "passes" : (fit.status == RoadFitStatus::kRejected ? "chi2" : "propagation")};
    };
    const auto count = static_cast<std::uint32_t>(road.size());
    const auto fit = fitRoad(parameters, propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = road.data()}, static_cast<int>(count));
    // Against the kink-free cut: a lower bound of the ratio to the cut the
    // fit applied, whose kink count a refused fit need not report.
    const double ratio = static_cast<double>(fit.chi2) / static_cast<double>(roadChi2Cut(parameters, count, 0U));
    std::vector<std::uint32_t> barrel;
    std::vector<std::uint32_t> other;
    for (const auto row : road) {
      (((parameters.barrelLayers >> mTruth->rows[row].layer) & 1U) != 0U ? barrel : other).push_back(row);
    }
    const auto pdg = std::abs(particle.pdg);
    const char* species = pdg == 11 ? "e" : (pdg == 13 ? "mu" : (pdg == 211 ? "pi" : (pdg == 321 ? "K" : (pdg == 2212 ? "p" : "other"))));
    ++mRoadFitRefusals[std::format("{} {} {} rows kink {} chi2/cut {} | barrel rows alone {} | other rows alone {}", particle.eta < 0.F ? "-z" : "+z", species, count, kink < 0.F ? "no hits" : (kink <= 3.F ? "<= 3" : (kink <= 5.F ? "3-5" : "> 5")), ratio < 1.5 ? "1-1.5" : (ratio < 2. ? "1.5-2" : (ratio < 4. ? "2-4" : "> 4")), fitOf(barrel), fitOf(other))];
  }

  /// A particle without a good track: its rows on the iteration's road
  /// layers held by a fake or another particle's track, fewer of them than
  /// the iteration's chains reach (`minClusters`), or the road fit stage's
  /// refusal of them (`roadFitRefusal`, Validation.h, in the stage's order).
  template <typename Material>
  [[nodiscard]] std::uint8_t missedReason(const std::vector<std::uint32_t>& own, const std::vector<std::uint32_t>& holder, const std::vector<std::uint64_t>& particles, std::size_t iteration, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator) const
  {
    const auto road = roadRows(own, iteration);
    bool fake = false;
    bool other = false;
    for (const auto row : road) {
      if (holder[row] != kNoTrack) {
        (particles[holder[row]] == study::kNoTrackParticle ? fake : other) = true;
      }
    }
    if (fake) {
      return kLostToFake;
    }
    if (other) {
      return kTakenByOther;
    }
    if (road.empty() || road.size() < parameters.minClusters) {
      return kTooFewRows;
    }
    switch (roadFitRefusal(parameters, propagator, geometry, material, frame, std::span<const std::uint32_t>{road})) {
      case FitRefusal::kTime:
        return kRoadTime;
      case FitRefusal::kChi2:
        return kRoadChi2;
      case FitRefusal::kPropagation:
        return kRoadPropagation;
      case FitRefusal::kPtFloor:
        return kFloor;
      case FitRefusal::kBeamLine:
        return kBeamLine;
      case FitRefusal::kTaken:
        break;
    }
    return kNoRoad;
  }

  /// The particle's rows \a own on iteration \a iteration's road layers,
  /// the first it crossed on each.
  [[nodiscard]] std::vector<std::uint32_t> roadRows(const std::vector<std::uint32_t>& own, std::size_t iteration) const
  {
    std::vector<std::uint32_t> road;
    std::uint64_t taken = 0U;
    const std::uint64_t roadLayers = iteration < mRoadLayers.size() ? mRoadLayers[iteration] : ~std::uint64_t{0};
    for (const auto row : own) {
      const auto layer = mTruth->rows[row].layer;
      if (((roadLayers >> layer) & 1U) == 0U || ((taken >> layer) & 1U) != 0U) {
        continue;
      }
      taken |= std::uint64_t{1} << layer;
      road.push_back(row);
    }
    return road;
  }

  /// The triplet cut's predicates after the time (TripletCut.h), the middle
  /// layer's scattering times \a scale; kGraphPasses with \a fit and \a terms.
  [[nodiscard]] static Graph tripletAudit(const IterationParameters<Det>& parameters, LayerID middle, const TimedPoint& a, const TimedPoint& b, const TimedPoint& c, Scalar scale, TripletFitResult& fit, KinkTerms& terms)
  {
    if (!tripletCurvatureCompatible(parameters, middle, a.point, b.point, c.point, scale)) {
      return kGraphCurvature;
    }
    if (!fitTriplet(a.point, b.point, c.point, middle, parameters.scattering, parameters.pairCut, fit, terms, scale)) {
      return kGraphFitFailed;
    }
    const Scalar nSigma = parameters.pairCut.nSigma;
    if (!(fit.chi2 <= nSigma * nSigma)) {
      return kGraphChi2;
    }
    if (!o2::gtrk::detail::tripletOwnMomentumCompatible(parameters, middle, a.point, b.point, c.point, scale, fit, terms)) {
      return kGraphOwnMomentum;
    }
    return kGraphPasses;
  }

  /// The graph's predicates on \a road's consecutive triplets and their
  /// neighbours, in the stages' order, as the iteration takes them: a
  /// triplet the nominal cuts refuse tried again with the kink allowance
  /// (its refusal then recorded), two kinked triplets never neighbours, a
  /// kinked middle point's scattering scaled in the adjacent fit; the dip
  /// window with each layer's largest variances in the frame, which bounds
  /// the slice's. Each decision is held to the stage's own cut
  /// (tripletCompatible, neighbourCompatible), disagreements counted.
  [[nodiscard]] Graph graphAudit(const std::vector<std::uint32_t>& road, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, std::uint64_t particle = 0U) const
  {
    if (parameters.roadStartLayers == 0U) {
      return auditChain(road, parameters, frame, particle, true);
    }
    // A start-layer mask (`iterationRoadStartLayers`, RoadCount.h
    // `roadCores`): a road whose innermost measurement lies on a layer roads
    // may not start on is none, and its suffix with it, since the roads are
    // the maximal chains. The chain the stages make of the true rows is the
    // suffix from the first row every predicate passes from; its innermost
    // row's layer decides.
    const auto need = std::max<std::size_t>(parameters.minClusters, 3U);
    for (std::size_t k = 0; k + need <= road.size(); ++k) {
      const std::vector<std::uint32_t> chain(road.begin() + static_cast<std::ptrdiff_t>(k), road.end());
      if (auditChain(chain, parameters, frame, particle, false) == kGraphPasses) {
        return ((parameters.roadStartLayers >> mTruth->rows[road[k]].layer) & 1U) != 0U ? kGraphPasses : kGraphStart;
      }
    }
    return auditChain(road, parameters, frame, particle, true);
  }

  /// The first predicate refusing the chain \a road of true rows of
  /// \a particle, in the stages' order; \a record whether the refusal
  /// histograms count it.
  [[nodiscard]] Graph auditChain(const std::vector<std::uint32_t>& road, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, std::uint64_t particle, bool record) const
  {
    const auto columns = frame.measurements();
    const auto layerOf = [&](std::uint32_t row) { return mTruth->rows[row].layer; };
    // A pair refusal, by the link, the particle's true kink over it and
    // whether it was produced within the beam line's z reach.
    const auto pairRefused = [&](const char* what, std::size_t link) {
      if (!record) {
        return;
      }
      const auto kink = largestKink(particle, Det::layerLinks[link].from, Det::layerLinks[link].to);
      const auto found = mTruth->particles.find(particle);
      const double vz = found != mTruth->particles.end() ? std::abs(static_cast<double>(found->second.vz)) : 0.;
      const double reach = static_cast<double>(parameters.pairCut.nSigma) * std::sqrt(static_cast<double>(mBeamLine.varZ));
      const double pt = truePt(found);
      ++mPairRefusals[std::format("pT {} the floor: {} {}->{} kink {} |vz| {} the beam line's {:.0f} cm", pt >= mAuditFloor ? "at or above" : "below", what, static_cast<int>(Det::layerLinks[link].from), static_cast<int>(Det::layerLinks[link].to), kink < 0.F ? "no hits" : (kink <= 3.F ? "<= 3" : "> 3"), vz <= reach ? "within" : "beyond", reach)];
    };
    const auto dip = [&](const TimedPoint& from, const TimedPoint& to) { return std::atan2(static_cast<double>(to.point.z - from.point.z), std::hypot(static_cast<double>(to.point.x - from.point.x), static_cast<double>(to.point.y - from.point.y))); };
    std::vector<TripletFitResult> fits;
    std::vector<KinkTerms> terms;
    std::vector<bool> kinked;
    const bool allowance = parameters.graphKinkScale2 > 1.F;
    // The tracklet stage's own tests on every consecutive pair first, in the
    // stage's order: the link's source range, the time, the cells
    // the query visits (the window rounded out to the grid's cells, its second
    // range clipped to the direct part where the iteration has no skip
    // links: `prepareQuery`), the pair cut from the beam line
    // (`pairCutLoss`; this schedule's iterations pair through it), and the
    // direct-pair rule on the pair itself (`segmentLeavesHole`).
    for (std::size_t i = 0; i + 1U < road.size(); ++i) {
      const auto link = linkOf<Det>(layerOf(road[i]), layerOf(road[i + 1U]));
      if (link == kNoLink) {
        return kGraphLambda;
      }
      if (!linkTaken(parameters, link)) {
        return kGraphLinkNotTaken;
      }
      const auto source = frame.measurement(road[i]);
      const auto pair = study::pairCutLoss(frame, parameters, mBeamLine, link, road[i], road[i + 1U]);
      const auto verdict = pair.loss;
      if (verdict == study::PairLoss::kSourceInterval && !record) {
        return kGraphLinkRange;
      }
      if (verdict == study::PairLoss::kSourceInterval) {
        const auto z = static_cast<double>(secondCoordinate<Det>(Det::layerLinks[link].from, source.r, source.z));
        const auto low = static_cast<double>(parameters.sourceMin[link]);
        const auto high = static_cast<double>(parameters.sourceMax[link]);
        ++mLinkRangeRefusals[std::format("{}->{} {} [{:.1f}, {:.1f}] by {:.1f}", static_cast<int>(Det::layerLinks[link].from), static_cast<int>(Det::layerLinks[link].to), z < low ? "below" : "above", low, high, z < low ? low - z : z - high)];
        return kGraphLinkRange;
      }
      if (verdict == study::PairLoss::kTime) {
        pairRefused("time", link);
        return kGraphPairTime;
      }
      switch (verdict) {
        case study::PairLoss::kWindow:
          pairRefused("outside the cells", link);
          return kGraphPairWindow;
        case study::PairLoss::kHole:
          pairRefused(pair.emptyDirectRange ? "no direct part" : "crosses another layer", link);
          return kGraphPairHole;
        case study::PairLoss::kTurning:
          pairRefused("turning point", link);
          return kGraphPairTurning;
        case study::PairLoss::kCurvature:
          pairRefused("curvature", link);
          return kGraphPairCurvature;
        case study::PairLoss::kLongitudinal:
          pairRefused("longitudinal", link);
          return kGraphPairLongitudinal;
        default:
          break;
      }
    }
    for (std::size_t i = 0; i + 2U < road.size(); ++i) {
      const auto a = timedPointAt(columns, road[i]);
      const auto b = timedPointAt(columns, road[i + 1U]);
      const auto c = timedPointAt(columns, road[i + 2U]);
      const auto middle = layerOf(road[i + 1U]);
      const auto firstLink = linkOf<Det>(layerOf(road[i]), middle);
      const auto second = linkOf<Det>(middle, layerOf(road[i + 2U]));
      if (firstLink == kNoLink || second == kNoLink) {
        return kGraphLambda;
      }
      // A link admits sources from a z range (`LayerLink::sourceMin/Max`,
      // traced from the MC transitions): a row outside it makes no pair,
      // whatever the predicates say.
      // The source coordinate is the layer's second one: z on a barrel
      // layer, r on a disc (`TrackletStage`, `secondCoordinate`).
      const auto inRange = [&](std::size_t link, const TimedPoint& source) {
        const auto z = static_cast<double>(secondCoordinate<Det>(Det::layerLinks[link].from, source.point.r, source.point.z));
        return z >= Det::layerLinks[link].sourceMin && z <= Det::layerLinks[link].sourceMax;
      };
      if (!inRange(firstLink, a) || !inRange(second, b)) {
        const auto link = inRange(firstLink, a) ? second : firstLink;
        const auto& source = inRange(firstLink, a) ? b : a;
        const auto z = static_cast<double>(secondCoordinate<Det>(Det::layerLinks[link].from, source.point.r, source.point.z));
        if (record) {
          ++mLinkRangeRefusals[std::format("{}->{} {} [{:.1f}, {:.1f}] by {:.1f}", static_cast<int>(Det::layerLinks[link].from), static_cast<int>(Det::layerLinks[link].to), z < Det::layerLinks[link].sourceMin ? "below" : "above", Det::layerLinks[link].sourceMin, Det::layerLinks[link].sourceMax, z < Det::layerLinks[link].sourceMin ? Det::layerLinks[link].sourceMin - z : z - Det::layerLinks[link].sourceMax)];
        }
        return kGraphLinkRange;
      }
      // GCC 14 finds a path with `middle` one past the last layer once the
      // half width's helpers are inlined; `second` holds a link
      // out of `middle`, so there is none.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
      const auto half = static_cast<double>(tripletLambdaHalfWidth(parameters, second, middle, a.point, b.point, mBounds[layerOf(road[i + 2U])]));
#pragma GCC diagnostic pop
      if (std::abs(dip(b, c) - dip(a, b)) > half) {
        return kGraphLambda;
      }
      if (!tripletTimeCompatible(a, c)) {
        return kGraphTime;
      }
      TripletFitResult fit;
      KinkTerms kink;
      auto reason = tripletAudit(parameters, middle, a, b, c, 1.F, fit, kink);
      bool withAllowance = false;
      if (reason != kGraphPasses && allowance) {
        reason = tripletAudit(parameters, middle, a, b, c, parameters.graphKinkScale2, fit, kink);
        withAllowance = reason == kGraphPasses;
      }
      TripletFitResult stageFit;
      KinkTerms stageTerms;
      bool stageKinked = false;
      const bool stage = tripletCompatible(parameters, middle, a, b, c, stageFit, stageTerms, stageKinked);
      mAuditDecisions += record ? 1U : 0U;
      mAuditDisagreements += record && (stage != (reason == kGraphPasses) || (stage && stageKinked != withAllowance)) ? 1U : 0U;
      if (reason == kGraphCurvature) {
        // The disc pass's seeds: which triplet, and the predicate's
        // numbers for the first few.
        if (record) {
          ++mCurvatureTriplets[std::format("{}-{}-{}", static_cast<int>(layerOf(road[i])), static_cast<int>(middle), static_cast<int>(layerOf(road[i + 2U])))];
        }
        if (record && mDebugTriplets < 16U) {
          ++mDebugTriplets;
          const auto chords = tripletChords(a.point, b.point, c.point);
          const double cross = (static_cast<double>(b.point.x - a.point.x) * static_cast<double>(c.point.y - b.point.y)) - (static_cast<double>(b.point.y - a.point.y) * static_cast<double>(c.point.x - b.point.x));
          const double l01 = chords.length01;
          const double l12 = chords.length12;
          const double l02 = chords.length02;
          const double curvature = 2. * std::abs(cross) / (l01 * l12 * l02);
          const double dz = static_cast<double>(c.point.z) - static_cast<double>(a.point.z);
          const double cosine2 = (l02 * l02) / ((l02 * l02) + (dz * dz));
          const double scattering = static_cast<double>(tripletScattering(parameters.scattering, parameters.pairCut, middle, a.point, b.point, c.point, static_cast<Scalar>(std::sqrt(cosine2)), 1.F)) / cosine2;
          const auto ux = static_cast<Scalar>(static_cast<double>(c.point.x - a.point.x) / l02);
          const auto uy = static_cast<Scalar>(static_cast<double>(c.point.y - a.point.y) / l02);
          const double inFirst = 1. / l01;
          const double inSecond = 1. / l12;
          const double inMiddle = inFirst + inSecond;
          const double m0 = static_cast<double>(o2::gtrk::detail::varianceAcross(a.point, ux, uy)) * inFirst * inFirst;
          const double m1 = static_cast<double>(o2::gtrk::detail::varianceAcross(b.point, ux, uy)) * inMiddle * inMiddle;
          const double m2 = static_cast<double>(o2::gtrk::detail::varianceAcross(c.point, ux, uy)) * inSecond * inSecond;
          const double tolerance = static_cast<double>(parameters.pairCut.nSigma) * (2. / l02) * std::sqrt(scattering + m0 + m1 + m2);
          const auto found = mTruth->particles.find(particle);
          std::cout << std::format(
            "  curvature refused: pT {:.3f} eta {:.2f}, layers {}-{}-{}, (r, z) a ({:.2f}, {:.1f}) b ({:.2f}, {:.1f}) c ({:.2f}, {:.1f}), chords {:.2f} {:.2f} "
            "{:.2f}, curvature {:.5f} max {:.5f}, tolerance {:.5f} = nSigma {:.1f} x 2/L x sqrt(scattering {:.2e} + meas {:.2e} {:.2e} {:.2e}), cos2 {:.3f}\n",
            truePt(found), found != mTruth->particles.end() ? static_cast<double>(found->second.eta) : 0., static_cast<int>(layerOf(road[i])), static_cast<int>(middle), static_cast<int>(layerOf(road[i + 2U])), static_cast<double>(a.point.r), static_cast<double>(a.point.z), static_cast<double>(b.point.r), static_cast<double>(b.point.z), static_cast<double>(c.point.r), static_cast<double>(c.point.z), l01, l12, l02, curvature, static_cast<double>(parameters.pairCut.maxCurvature), tolerance, static_cast<double>(parameters.pairCut.nSigma), scattering, m0, m1, m2, cosine2);
        }
      }
      if (reason != kGraphPasses) {
        return reason;
      }
      fits.push_back(fit);
      terms.push_back(kink);
      kinked.push_back(withAllowance);
    }
    for (std::size_t i = 0; i + 3U < road.size(); ++i) {
      if (kinked[i] && kinked[i + 1U]) {
        return kGraphBothKinked;
      }
      if (!neighbourCurvatureCompatible(parameters, fits[i].curvature, fits[i].varCurvature, fits[i + 1U].curvature, fits[i + 1U].varCurvature)) {
        return kGraphNeighbourCurvature;
      }
      const auto a = timedPointAt(columns, road[i]);
      const auto b = timedPointAt(columns, road[i + 1U]);
      const auto c = timedPointAt(columns, road[i + 2U]);
      const auto d = timedPointAt(columns, road[i + 3U]);
      const Scalar secondScale = kinked[i] ? parameters.graphKinkScale2 : 1.F;
      const Scalar thirdScale = kinked[i + 1U] ? parameters.graphKinkScale2 : 1.F;
      auto reason = kGraphPasses;
      TripletFitResult adjacent;
      if (!neighbourTimeCompatible(a, d)) {
        reason = kGraphNeighbourTime;
      } else if (!fitAdjacentTriplets(terms[i], terms[i + 1U], a.point, b.point, c.point, d.point, layerOf(road[i + 1U]), layerOf(road[i + 2U]), parameters.scattering, parameters.pairCut, adjacent, secondScale, thirdScale)) {
        reason = kGraphNeighbourFitFailed;
      } else if (!(adjacent.chi2 <= parameters.neighbourChi2Max)) {
        reason = kGraphNeighbourChi2;
      }
      // The stage's own cut (NeighbourCut.h), its kink terms made from the points.
      TripletFitResult stageFit;
      KinkTerms stageTerms;
      bool stageMade = false;
      bool stageValid = false;
      const bool stage = neighbourCompatible(parameters, layerOf(road[i + 1U]), layerOf(road[i + 2U]), a, b, c, d, stageTerms, stageMade, stageValid, stageFit, secondScale, thirdScale);
      mAuditDecisions += record ? 1U : 0U;
      mAuditDisagreements += record && stage != (reason == kGraphPasses) ? 1U : 0U;
      if (reason != kGraphPasses) {
        return reason;
      }
    }
    return kGraphPasses;
  }

  /// \a rows of \a particle in the order it crossed them: by the distance
  /// from its production vertex -- `pathOrder`'s fallback, which falls on a
  /// looper's way back.
  void sortByTraversal(std::uint64_t particle, std::vector<std::uint32_t>& rows) const
  {
    const auto found = mTruth->particles.find(particle);
    const auto vertex = vertexRZ(found);
    const auto distance = [&](std::uint32_t row) { return pathDistance(row, vertex); };
    std::ranges::sort(rows, [&](std::uint32_t a, std::uint32_t b) { return distance(a) < distance(b) || (!(distance(b) < distance(a)) && a < b); });
  }

  /// The ensemble without gate selection: every
  /// findable particle along its own rows inward from its true road, the
  /// rows taken whatever their chi2; at each target inner layer three
  /// predictions held to the MC crossing: the follow's construction (the
  /// road's refit updated once with the beam line, then filtered through the
  /// rows), the clusters-only fit, and the joint reference (the
  /// clusters-only prediction updated at the target with the beam line).
  template <typename Material>
  void truthWalk(std::size_t iteration, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    mWalkGate = parameters.clusterChi2Max;
    std::vector<float> worst(mFindable.size(), 0.F);
    for (const auto& kink : mKinks) {
      worst[kink.particle] = std::max({worst[kink.particle], std::abs(kink.theta), std::abs(kink.phi)});
    }
    // Each particle's walks on its own thread, the sums, the pulls and the
    // CSV lines gathered per thread and merged after (the hits' events are
    // all read before, `open`); the CSV lines per range of particles, written
    // in the particles' order.
    struct Local {
      std::array<std::array<std::array<WalkSum, kLayers>, 3>, 3> walk{};
      std::array<std::array<std::array<WalkPulls, kWalkPtEdges.size() + 1U>, kLayers>, 3> pulls{};
      std::vector<std::pair<std::size_t, std::string>> csv;
    };
    std::vector<Local> locals(static_cast<std::size_t>(std::max(1, mArena->maxConcurrency())) + 1U);
    o2::gtrk::detail::parallelFor(*mArena, 0U, mFindable.size(), 64, [&](std::size_t a, std::size_t b) {
      auto& mine = locals[static_cast<std::size_t>(std::max(0, o2::gtrk::detail::currentArenaThreadIndex()))];
      std::string& csv = mine.csv.emplace_back(a, std::string{}).second;
      for (auto index = static_cast<std::uint32_t>(a); index < b; ++index) {
        const auto particle = mFindable[index];
        const auto& own = mRowsOf.at(particle);
        const auto road = roadRows(own, iteration);
        if (road.size() < 4U) {
          continue;
        }
        const auto firstRoadLayer = mTruth->rows[road.front()].layer;
        // One row per layer inside the road, innermost first.
        std::array<std::uint32_t, kLayers> inner{};
        inner.fill(kNoTrack);
        for (const auto row : own) {
          const auto layer = mTruth->rows[row].layer;
          if (layer < firstRoadLayer && ((mMask >> layer) & 1U) != 0U && inner[layer] == kNoTrack) {
            inner[layer] = row;
          }
        }
        // The follow's construction, input along the rows.
        const auto roadFit = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = road.data()}, static_cast<int>(road.size()), kinkBudget(parameters));
        if (roadFit.status != PropagationStatus::kOk) {
          continue;
        }
        BoundTrackState input = roadFit.state;
        static_cast<void>(updateWithBeamLine(propagator, parameters, geometry.planeOf(input), geometry.surface(input.surface).layer, pairPointAt(frame.measurements(), road.front()), input));
        CovarianceFactors inputFactors = factorCovariance(input.cov);
        const bool quiet = worst[index] < 2.F;
        // Produced within 50 um of the axis, 2.5 of the beam's widths.
        const auto truthOf = mTruth->particles.find(particle);
        const bool prompt = truthOf != mTruth->particles.end() && std::hypot(truthOf->second.vx, truthOf->second.vy) < 0.005F;
        const double pt = truePt(truthOf);
        const double eta = truthOf != mTruth->particles.end() ? static_cast<double>(truthOf->second.eta) : 0.;
        const int pdg = truthOf != mTruth->particles.end() ? truthOf->second.pdg : 0;
        const auto ptBin = static_cast<std::size_t>(std::ranges::upper_bound(kWalkPtEdges, pt) - kWalkPtEdges.begin());
        // The step's prediction at the
        // particle's own row \a row from \a state in \a dir, as the follow's
        // step makes it (`transportCandidate`, the nominal noise), against the
        // cluster -- the residuals, the innovation's diagonal, the update's
        // chi2, the gate's quantity -- and against the MC crossing; the true
        // kink over the step and over the one before. One CSV line and the
        // summary's pulls; whether the transport worked.
        const auto walkStep = [&](std::size_t direction, PropagationDirection dir, LayerID fromLayer, std::uint32_t row, const BoundTrackState& state, const CovarianceFactors& stateFactors) {
          const auto measurement = measurementAt(frame, row);
          const auto plane = geometry.planeOf(measurement);
          const auto layer = mTruth->rows[row].layer;
          const float local = largestKink(particle, std::min(fromLayer, layer), std::max(fromLayer, layer));
          const float upstream = kinkAt(particle, fromLayer, direction != 0U);
          std::size_t onLayer = 0U;
          for (const auto at : own) {
            onLayer += mTruth->rows[at].layer == layer ? 1U : 0U;
          }
          auto& pulls = mine.pulls[direction][layer][ptBin];
          BoundTrackState predicted;
          CovarianceFactors predictedFactors;
          if (!follow::transportCandidate<false>(propagator, geometry, material, plane, dir, state, stateFactors, predicted, predictedFactors, 1.F)) {
            ++pulls.fails;
            if (mWalkCsv) {
              csv += std::format("{},{},{},{:.4f},{:.3f},{},{},{},1,0,0,0,0,0,0,0,0,0,{:.3f},{:.3f},{}\n", direction, iteration, particle, pt, eta, pdg, static_cast<int>(fromLayer), static_cast<int>(layer), static_cast<double>(local), static_cast<double>(upstream), onLayer);
            }
            return false;
          }
          covarianceFromFactors(predictedFactors, predicted.cov);
          const double r0 = static_cast<double>(measurement.loc0) - static_cast<double>(predicted.loc0());
          const double r1 = static_cast<double>(measurement.loc1) - static_cast<double>(predicted.loc1());
          const auto p00 = static_cast<double>(predicted.covAt(kLoc0, kLoc0));
          const auto p11 = static_cast<double>(predicted.covAt(kLoc1, kLoc1));
          const double i00 = p00 + static_cast<double>(measurement.varLoc0);
          const double i11 = p11 + static_cast<double>(measurement.varLoc1);
          // The update on the factors, as the follow's step makes it
          // (`updateCandidate`), and its usability.
          BoundTrackState updated = predicted;
          CovarianceFactors updatedFactors = predictedFactors;
          const auto update = updateWithMeasurement(updated, updatedFactors, measurement);
          const double chi2 = updateUsable(update, updatedFactors) ? static_cast<double>(update.chi2) : -1.;
          const double pull0 = i00 > 0. ? r0 / std::sqrt(i00) : 0.;
          const double pull1 = i11 > 0. ? r1 / std::sqrt(i11) : 0.;
          pulls.pull0.push_back(static_cast<float>(pull0));
          pulls.pull1.push_back(static_cast<float>(pull1));
          pulls.overGate += chi2 > static_cast<double>(mWalkGate) ? 1U : 0U;
          double c0 = -99.;
          double c1 = -99.;
          if (const auto* segment = segmentOfRow(particle, frame, row); segment != nullptr) {
            const auto crossing = crossingOn(*segment, plane);
            c0 = p00 > 0. ? (static_cast<double>(predicted.loc0()) - crossing[0]) / std::sqrt(p00) : 0.;
            c1 = p11 > 0. ? (static_cast<double>(predicted.loc1()) - crossing[1]) / std::sqrt(p11) : 0.;
          }
          if (mWalkCsv) {
            csv += std::format("{},{},{},{:.4f},{:.3f},{},{},{},0,{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{:.4f},{:.3f},{:.3f},{:.3f},{:.3f},{}\n", direction, iteration, particle, pt, eta, pdg, static_cast<int>(fromLayer), static_cast<int>(layer), r0, r1, i00, i11, p00, p11, chi2, c0, c1, static_cast<double>(local), static_cast<double>(upstream), onLayer);
          }
          return true;
        };
        // Outward from the road's outward state along the particle's own rows
        // beyond the road, in crossing order, one per layer of the mask, each
        // taken whatever its chi2.
        {
          BoundTrackState state = roadFit.outwardState;
          CovarianceFactors factors = factorCovariance(state.cov);
          std::uint64_t seen = 0U;
          for (const auto row : road) {
            seen |= std::uint64_t{1} << mTruth->rows[row].layer;
          }
          bool beyond = false;
          for (const auto row : own) {
            if (!beyond) {
              beyond = row == road.back();
              continue;
            }
            const auto layer = mTruth->rows[row].layer;
            if (((mMask >> layer) & 1U) == 0U || ((seen >> layer) & 1U) != 0U) {
              continue;
            }
            seen |= std::uint64_t{1} << layer;
            const LayerID from = geometry.surface(state.surface).layer;
            if (!walkStep(1U, PropagationDirection::kForward, from, row, state, factors)) {
              break;
            }
            BoundTrackState output;
            CovarianceFactors outputFactors;
            Scalar chi2 = 0.F;
            if (!follow::updateCandidate<false>(propagator, geometry, material, measurementAt(frame, row), PropagationDirection::kForward, state, factors, output, outputFactors, chi2, 1.F)) {
              break;
            }
            state = output;
            factors = outputFactors;
            covarianceFromFactors(factors, state.cov);
          }
        }
        // The disc class's road-fit refusals: the road's barrel rows
        // fitted alone, then the step to each of its other rows (the discs) in
        // crossing order, each taken whatever its chi2 (direction 2).
        {
          std::vector<std::uint32_t> barrel;
          std::vector<std::uint32_t> other;
          for (const auto row : road) {
            (((parameters.barrelLayers >> mTruth->rows[row].layer) & 1U) != 0U ? barrel : other).push_back(row);
          }
          if (barrel.size() >= 4U && !other.empty()) {
            const auto barrelFit = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = barrel.data()}, static_cast<int>(barrel.size()), kinkBudget(parameters));
            if (barrelFit.status == PropagationStatus::kOk) {
              BoundTrackState state = barrelFit.outwardState;
              CovarianceFactors factors = factorCovariance(state.cov);
              for (const auto row : other) {
                const LayerID from = geometry.surface(state.surface).layer;
                if (!walkStep(2U, PropagationDirection::kForward, from, row, state, factors)) {
                  break;
                }
                BoundTrackState output;
                CovarianceFactors outputFactors;
                Scalar chi2 = 0.F;
                if (!follow::updateCandidate<false>(propagator, geometry, material, measurementAt(frame, row), PropagationDirection::kForward, state, factors, output, outputFactors, chi2, 1.F)) {
                  break;
                }
                state = output;
                factors = outputFactors;
                covarianceFromFactors(factors, state.cov);
              }
            }
          }
        }
        // Inward layer by layer by number: the class's barrel layers
        // consecutive from the beam outward (a skipped layer or a disc is
        // not supported).
        for (auto target = static_cast<int>(firstRoadLayer) - 1; target >= 0; --target) {
          const auto row = inner[static_cast<std::size_t>(target)];
          if (row == kNoTrack) {
            break;
          }
          const auto measurement = measurementAt(frame, row);
          const auto plane = geometry.planeOf(measurement);
          const auto* segment = segmentOfRow(particle, frame, row);
          static_cast<void>(walkStep(0U, PropagationDirection::kBackward, geometry.surface(input.surface).layer, row, input, inputFactors));
          // The clusters-only fit over the rows outside the target.
          std::vector<std::uint32_t> chain;
          for (auto layer = static_cast<std::size_t>(target) + 1U; layer < static_cast<std::size_t>(firstRoadLayer); ++layer) {
            chain.push_back(inner[layer]);
          }
          chain.insert(chain.end(), road.begin(), road.end());
          const auto clusters = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = chain.data()}, static_cast<int>(chain.size()), kinkBudget(parameters));
          std::array<BoundTrackState, 3> predicted{};
          std::array<bool, 3> ok{};
          CovarianceFactors factors;
          ok[0] = follow::transportCandidate<false>(propagator, geometry, material, plane, PropagationDirection::kBackward, input, inputFactors, predicted[0], factors, 1.F);
          if (ok[0]) {
            covarianceFromFactors(factors, predicted[0].cov);
          }
          if (clusters.status == PropagationStatus::kOk) {
            ok[1] = follow::transportCandidate<false>(propagator, geometry, material, plane, PropagationDirection::kBackward, clusters.state, factorCovariance(clusters.state.cov), predicted[1], factors, 1.F);
            if (ok[1]) {
              covarianceFromFactors(factors, predicted[1].cov);
              predicted[2] = predicted[1];
              ok[2] = updateWithBeamLine(propagator, parameters, plane, static_cast<LayerID>(target), pairPointAt(frame.measurements(), row), predicted[2]).ok;
            }
          }
          if (segment != nullptr) {
            const auto crossing = crossingOn(*segment, plane);
            for (std::size_t construction = 0; construction < 3U; ++construction) {
              if (!ok[construction]) {
                continue;
              }
              const auto& state = predicted[construction];
              const double pull0 = (static_cast<double>(state.loc0()) - crossing[0]) / std::sqrt(static_cast<double>(state.covAt(kLoc0, kLoc0)));
              const double pull1 = (static_cast<double>(state.loc1()) - crossing[1]) / std::sqrt(static_cast<double>(state.covAt(kLoc1, kLoc1)));
              for (std::size_t selection = 0; selection < 3U; ++selection) {
                if ((selection >= 1U && !quiet) || (selection == 2U && !prompt)) {
                  continue;
                }
                auto& sum = mine.walk[selection][construction][static_cast<std::size_t>(target)];
                ++sum.n;
                sum.s0 += pull0 * pull0;
                sum.s1 += pull1 * pull1;
                sum.out0 += std::abs(pull0) > 3. ? 1U : 0U;
                sum.out1 += std::abs(pull1) > 3. ? 1U : 0U;
              }
            }
          }
          // The follow's construction takes the row, whatever its chi2.
          BoundTrackState output;
          CovarianceFactors outputFactors;
          Scalar chi2 = 0.F;
          if (!follow::updateCandidate<false>(propagator, geometry, material, measurement, PropagationDirection::kBackward, input, inputFactors, output, outputFactors, chi2, 1.F)) {
            break;
          }
          input = output;
          inputFactors = outputFactors;
          covarianceFromFactors(inputFactors, input.cov);
        }
      }
    });
    for (auto& each : locals) {
      for (std::size_t s = 0; s < 3U; ++s) {
        for (std::size_t c = 0; c < 3U; ++c) {
          for (std::size_t layer = 0; layer < kLayers; ++layer) {
            auto& sum = mWalk[s][c][layer];
            const auto& one = each.walk[s][c][layer];
            sum.n += one.n;
            sum.s0 += one.s0;
            sum.s1 += one.s1;
            sum.out0 += one.out0;
            sum.out1 += one.out1;
          }
        }
      }
      for (std::size_t direction = 0; direction < 3U; ++direction) {
        for (std::size_t layer = 0; layer < kLayers; ++layer) {
          for (std::size_t bin = 0; bin <= kWalkPtEdges.size(); ++bin) {
            auto& pulls = mWalkPulls[direction][layer][bin];
            auto& one = each.pulls[direction][layer][bin];
            pulls.pull0.insert(pulls.pull0.end(), one.pull0.begin(), one.pull0.end());
            pulls.pull1.insert(pulls.pull1.end(), one.pull1.begin(), one.pull1.end());
            pulls.overGate += one.overGate;
            pulls.fails += one.fails;
          }
        }
      }
    }
    if (mWalkCsv) {
      std::vector<std::pair<std::size_t, std::string>*> ranges;
      for (auto& each : locals) {
        for (auto& range : each.csv) {
          ranges.push_back(&range);
        }
      }
      std::ranges::sort(ranges, {}, [](const auto* range) { return range->first; });
      for (const auto* range : ranges) {
        *mWalkCsv << range->second;
      }
    }
  }

  /// The MC segment's crossing of \a plane in its local coordinates.
  [[nodiscard]] static std::array<double, 2> crossingOn(const study::HitSegment& segment, const TrackingSurface& plane)
  {
    double along[3]{};
    double offset = 0.;
    double slope = 0.;
    for (std::size_t i = 0; i < 3U; ++i) {
      along[i] = segment.end[i] - segment.start[i];
      offset += (segment.start[i] - static_cast<double>(plane.origin[i])) * static_cast<double>(plane.normal[i]);
      slope += along[i] * static_cast<double>(plane.normal[i]);
    }
    const double t = std::abs(slope) > 0. ? -offset / slope : 0.5;
    const Scalar global[3] = {static_cast<Scalar>(segment.start[0] + (t * along[0])), static_cast<Scalar>(segment.start[1] + (t * along[1])), static_cast<Scalar>(segment.start[2] + (t * along[2]))};
    Scalar local[3]{};
    surfaceToLocal(plane, global, local);
    return {static_cast<double>(local[0]), static_cast<double>(local[1])};
  }

  /// The particle's first row on \a layer in its crossing order; none for a
  /// particle without rows (a fake's majority of noise rows).
  [[nodiscard]] std::uint32_t ownRowOn(std::uint64_t particle, LayerID layer) const
  {
    const auto rows = mRowsOf.find(particle);
    if (rows == mRowsOf.end()) {
      return kNoTrack;
    }
    for (const auto row : rows->second) {
      if (mTruth->rows[row].layer == layer) {
        return row;
      }
    }
    return kNoTrack;
  }

  /// The MC segment of \a particle on \a chip nearest \a near; none without the hits.
  [[nodiscard]] const study::HitSegment* segmentOf(std::uint64_t particle, std::uint32_t chip, const Scalar near[3]) const
  {
    const auto* segments = mHits ? segmentsOf(particle) : nullptr;
    if (segments == nullptr) {
      return nullptr;
    }
    const study::HitSegment* best = nullptr;
    double bestDistance = std::numeric_limits<double>::infinity();
    for (const auto* segment : *segments) {
      if (segment->chip != chip) {
        continue;
      }
      double d2 = 0.;
      for (std::size_t i = 0; i < 3U; ++i) {
        const double d = (0.5 * (segment->start[i] + segment->end[i])) - static_cast<double>(near[i]);
        d2 += d * d;
      }
      if (d2 < bestDistance) {
        bestDistance = d2;
        best = segment;
      }
    }
    return best;
  }

  [[nodiscard]] static std::uint64_t segmentKey(int event, int track) { return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(event)) << 32U) | static_cast<std::uint32_t>(track); }

  /// The radial direction of \a segment: +1 outward, -1 inward; from its
  /// ends, or where they are at one radius (a crossing at the turning point)
  /// from its momentum.
  [[nodiscard]] static int radialSign(const study::HitSegment& segment)
  {
    const double dr = std::hypot(segment.end[0], segment.end[1]) - std::hypot(segment.start[0], segment.start[1]);
    if (std::abs(dr) > 1e-4) {
      return dr > 0. ? 1 : -1;
    }
    const double along = (0.5 * (segment.start[0] + segment.end[0]) * segment.momentum[0]) + (0.5 * (segment.start[1] + segment.end[1]) * segment.momentum[1]);
    return along >= 0. ? 1 : -1;
  }

  /// Every particle's rows in path order: by its MC
  /// segments' time of flight, where every row has its segment; the distance
  /// from the vertex the constructor sorted by falls on a looper's way back,
  /// which comes nearer the vertex in r while it goes on in z and in time. A
  /// particle with a row the hits do not have keeps the distance's order. For
  /// the class, each row's leg (`mLegs`): 0 out to the first turning point, 1
  /// back from it, 2 out again, ...: where the radius along the path turns.
  void pathOrder(const TimeFrame<Det>& frame)
  {
    std::vector<std::pair<const std::uint64_t, std::vector<std::uint32_t>>*> all;
    all.reserve(mRowsOf.size());
    for (auto& entry : mRowsOf) {
      all.push_back(&entry);
    }
    // Per particle: 0 timed, its order the distance's; 1 timed and
    // reordered; 2 kept (a row without its segment); plus 4 where the first
    // row on some layer changed.
    std::vector<std::uint8_t> outcome(all.size(), 0U);
    mLegs.assign(mFindable.size(), {});
    mTimed.assign(mFindable.size(), 0U);
    o2::gtrk::detail::parallelFor(*mArena, 0U, all.size(), 1024, [&](std::size_t a, std::size_t b) {
      std::vector<std::pair<double, std::uint32_t>> keyed;
      std::vector<int> signs;
      for (std::size_t i = a; i < b; ++i) {
        const auto particle = all[i]->first;
        auto& rows = all[i]->second;
        keyed.clear();
        signs.clear();
        for (const auto row : rows) {
          const auto* segment = segmentOfRow(particle, frame, row);
          if (segment == nullptr) {
            break;
          }
          keyed.emplace_back(segment->time, row);
          signs.push_back(radialSign(*segment));
        }
        if (keyed.size() != rows.size()) {
          outcome[i] = 2U;
          continue;
        }
        std::array<std::uint32_t, kLayers> firstBefore{};
        firstBefore.fill(kNoTrack);
        for (const auto row : rows) {
          auto& first = firstBefore[mTruth->rows[row].layer];
          first = first == kNoTrack ? row : first;
        }
        // Sorted by time, the row as the tie-break; the signs with them.
        std::vector<std::size_t> order(rows.size());
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::ranges::sort(order, [&](std::size_t x, std::size_t y) { return keyed[x] < keyed[y]; });
        bool same = true;
        for (std::size_t k = 0; k < rows.size(); ++k) {
          same = same && rows[k] == keyed[order[k]].second;
        }
        std::vector<int> sortedSigns(rows.size());
        for (std::size_t k = 0; k < rows.size(); ++k) {
          rows[k] = keyed[order[k]].second;
          sortedSigns[k] = signs[order[k]];
        }
        std::array<std::uint32_t, kLayers> firstAfter{};
        firstAfter.fill(kNoTrack);
        for (const auto row : rows) {
          auto& first = firstAfter[mTruth->rows[row].layer];
          first = first == kNoTrack ? row : first;
        }
        outcome[i] = static_cast<std::uint8_t>((same ? 0U : 1U) | (firstBefore != firstAfter ? 4U : 0U));
        const auto index = mIndex.find(particle);
        if (index == mIndex.end()) {
          continue;
        }
        auto& legs = mLegs[index->second];
        legs.resize(rows.size());
        std::uint8_t leg = 0U;
        for (std::size_t k = 0; k < rows.size(); ++k) {
          if (k != 0U && sortedSigns[k] != sortedSigns[k - 1U] && leg < 255U) {
            ++leg;
          }
          legs[k] = leg;
        }
        mTimed[index->second] = 1U;
      }
    });
    for (std::size_t i = 0; i < all.size(); ++i) {
      const bool inClass = mIndex.contains(all[i]->first);
      auto& counts = mPathCounts[inClass ? 1U : 0U];
      ++counts[(outcome[i] & 3U) == 2U ? 2U : (outcome[i] & 1U)];
      counts[3] += (outcome[i] & 4U) != 0U ? 1U : 0U;
    }
  }

  [[nodiscard]] const study::HitSegment* segmentOfRow(std::uint64_t particle, const TimeFrame<Det>& frame, std::uint32_t row) const
  {
    const auto point = frame.measurement(row);
    const Scalar near[3] = {point.x, point.y, point.z};
    return segmentOf(particle, measurementAt(frame, row).surface, near);
  }

  static float midOf(const study::HitSegment& one, std::size_t i) { return static_cast<float>(0.5 * (one.start[i] + one.end[i])); }

  /// Highland's theta0 at \a p for the LUT's x/X0 between two crossings, and that x/X0.
  template <typename Material>
  static std::pair<double, double> theta0Between(const Material& material, const Propagator& propagator, const study::HitSegment& a, const study::HitSegment& b, double p)
  {
    const auto budget = material.getMatBudget(midOf(a, 0), midOf(a, 1), midOf(a, 2), midOf(b, 0), midOf(b, 1), midOf(b, 2));
    const auto x2x0 = static_cast<double>(budget.meanX2X0);
    const auto mass = static_cast<double>(propagator.massGeV);
    const double betaP = p * p / std::sqrt((p * p) + (mass * mass));
    return {x2x0 > 0. ? 0.0136 / betaP * std::sqrt(x2x0) * (1. + (0.038 * std::log(x2x0))) : 0., x2x0};
  }

  /// The true kinks from crossing \a a to crossing \a b, theta and the
  /// bending plane's, in radians: the momentum's change less the field's
  /// turn over the transverse arc at the true pT, its sign the one leaving
  /// the smaller kink (the turn is tens of mrad, a kink a few).
  static std::pair<double, double> kinksBetween(const study::HitSegment& a, const study::HitSegment& b, double bzKiloGauss)
  {
    const auto& pa = a.momentum;
    const auto& pb = b.momentum;
    const double theta = std::atan2(std::hypot(pb[0], pb[1]), pb[2]) - std::atan2(std::hypot(pa[0], pa[1]), pa[2]);
    const double pt = std::hypot(pa[0], pa[1]);
    const double p = std::hypot(pa[0], pa[1], pa[2]);
    const double radius = 100. * pt / (0.0299792458 * std::abs(bzKiloGauss));
    const double chord = std::hypot(static_cast<double>(midOf(b, 0) - midOf(a, 0)), static_cast<double>(midOf(b, 1) - midOf(a, 1)));
    const double turn = 2. * std::asin(std::min(1., 0.5 * chord / radius));
    const double dphi = std::remainder(std::atan2(pb[1], pb[0]) - std::atan2(pa[1], pa[0]), 2. * std::numbers::pi);
    const double phi = (std::abs(dphi - turn) < std::abs(dphi + turn) ? dphi - turn : dphi + turn) * pt / p;
    return {theta, phi};
  }

  /// Every findable particle's kinks between its crossings of consecutive
  /// layers of the mask, its first row on each in crossing order; which pairs
  /// were evaluated, a bit each; and whether it crossed a layer twice.
  template <typename Material>
  void trueKinks(const TimeFrame<Det>& frame, const Material& material, const Propagator& propagator)
  {
    mEvaluated.assign(mFindable.size(), 0U);
    mCrossed.assign(mFindable.size(), 0U);
    mRepeated.assign(mFindable.size(), false);
    for (std::uint32_t index = 0; index < mFindable.size(); ++index) {
      const auto particle = mFindable[index];
      std::array<const study::HitSegment*, kLayers> crossing{};
      std::uint64_t seen = 0U;
      // Every layer the particle crossed, the class's mask or not: a road
      // holds rows beyond the mask (the VD under a disc class), and its
      // "quiet" has to cover them.
      for (const auto row : mRowsOf.at(particle)) {
        const auto layer = mTruth->rows[row].layer;
        if (((seen >> layer) & 1U) != 0U) {
          mRepeated[index] = true;
          continue;
        }
        seen |= std::uint64_t{1} << layer;
        crossing[layer] = segmentOfRow(particle, frame, row);
        mMissingSegment += crossing[layer] == nullptr ? 1U : 0U;
      }
      // Between consecutive crossings of the mask's layers, whatever lies
      // between them (the disc classes' ML 4 -> C1 is not between layers L
      // and L + 1).
      std::size_t previous = kLayers;
      mCrossed[index] = seen;
      for (std::size_t layer = 0; layer < kLayers; ++layer) {
        if (crossing[layer] == nullptr) {
          continue;
        }
        const auto from = previous;
        previous = layer;
        if (from == kLayers) {
          continue;
        }
        const auto* a = crossing[from];
        const auto* b = crossing[layer];
        const double p = std::hypot(a->momentum[0], a->momentum[1], a->momentum[2]);
        const auto [theta0, x2x0] = theta0Between(material, propagator, *a, *b, p);
        if (!(theta0 > 0.)) {
          ++mZeroMaterial;
          continue;
        }
        const auto [theta, phi] = kinksBetween(*a, *b, static_cast<double>(propagator.bzKiloGauss));
        mEvaluated[index] |= std::uint64_t{1} << from;
        mKinks.push_back({
          .particle = index,
          .layer = static_cast<LayerID>(from),
          .to = static_cast<LayerID>(layer),
          .p = static_cast<float>(p),
          .theta = static_cast<float>(theta / theta0),
          .phi = static_cast<float>(phi / theta0),
          .x2x0 = static_cast<float>(x2x0),
        });
      }
    }
    // The kinks are pushed in particle order: each particle's first, for
    // `largestKink`.
    mKinkBegin.assign(mFindable.size() + 1U, 0U);
    for (const auto& kink : mKinks) {
      ++mKinkBegin[kink.particle + 1U];
    }
    for (std::size_t i = 1; i < mKinkBegin.size(); ++i) {
      mKinkBegin[i] += mKinkBegin[i - 1U];
    }
  }

  void printKinks() const
  {
    if (mKinks.empty()) {
      return;
    }
    // Every pair evaluated: a kink from every layer crossed but the last.
    const auto wholeOf = [&](std::size_t index) {
      const auto crossed = mCrossed[index];
      const auto pairs = crossed == 0U ? 0U : crossed & ~(std::uint64_t{1} << (63 - std::countl_zero(crossed)));
      return mEvaluated[index] == pairs;
    };
    std::size_t whole = 0U;
    std::size_t repeated = 0U;
    for (std::size_t index = 0; index < mFindable.size(); ++index) {
      whole += wholeOf(index) ? 1U : 0U;
      repeated += mRepeated[index] ? 1U : 0U;
    }
    std::cout << std::format(
      "  true kinks (MC hits): {} of {} particles with every pair evaluated; rows without a segment {}, pairs without material {}; particles "
      "crossing a layer twice {}\n",
      whole, mFindable.size(), mMissingSegment, mZeroMaterial, repeated);
    std::cout << "  kinks over theta0 between consecutive crossings (Gaussian: 0.27 % beyond 3): theta rms within 3, % beyond 3, beyond 5 | the bending "
                 "plane's the same\n";
    const auto summary = [&](const auto& select) {
      std::array<std::size_t, 2> count{};
      std::array<double, 2> core2{};
      std::array<std::array<std::size_t, 2>, 2> beyond{};
      std::size_t all = 0U;
      for (const auto& kink : mKinks) {
        if (!select(kink)) {
          continue;
        }
        ++all;
        for (std::size_t plane = 0; plane < 2U; ++plane) {
          const auto a = std::abs(static_cast<double>(plane == 0U ? kink.theta : kink.phi));
          if (a < 3.) {
            ++count[plane];
            core2[plane] += a * a;
          }
          beyond[plane][0] += a > 3. ? 1U : 0U;
          beyond[plane][1] += a > 5. ? 1U : 0U;
        }
      }
      const auto share = [&](std::size_t k) { return all != 0U ? 100. * static_cast<double>(k) / static_cast<double>(all) : 0.; };
      const auto rms = [&](std::size_t plane) { return count[plane] != 0U ? std::sqrt(core2[plane] / static_cast<double>(count[plane])) : 0.; };
      return std::format("{}: {:.2f} {:.2f} {:.3f} | {:.2f} {:.2f} {:.3f}", all, rms(0), share(beyond[0][0]), share(beyond[0][1]), rms(1), share(beyond[1][0]), share(beyond[1][1]));
    };
    for (std::size_t layer = 0; layer + 1U < kLayers; ++layer) {
      if (std::ranges::any_of(mKinks, [&](const TrueKink& kink) { return kink.layer == layer; })) {
        std::cout << std::format("    crossings from {}, {}\n", layer, summary([&](const TrueKink& kink) { return kink.layer == layer; }));
      }
    }
    constexpr std::array<double, 6> kPEdges{1., 1.5, 2., 3., 5., 1.e9};
    for (std::size_t bin = 0; bin + 1U < kPEdges.size(); ++bin) {
      std::cout << std::format("    p [{}, {}) GeV, {}\n", kPEdges[bin], kPEdges[bin + 1U], summary([&](const TrueKink& kink) { return static_cast<double>(kink.p) >= kPEdges[bin] && static_cast<double>(kink.p) < kPEdges[bin + 1U]; }));
    }
    // The particles by their largest kink, either plane; those without every
    // pair evaluated apart.
    std::vector<float> worst(mFindable.size(), 0.F);
    for (const auto& kink : mKinks) {
      worst[kink.particle] = std::max({worst[kink.particle], std::abs(kink.theta), std::abs(kink.phi)});
    }
    std::cout << "  the findable by their largest kink, either plane (theta0): n (% of all) | complete, short, not found (% of them); short and not found by reason\n";
    constexpr std::array<double, 7> kWorstEdges{0., 2., 3., 4., 5., 10., 1.e30};
    const auto row = [&](const std::string& name, const auto& select) {
      std::array<std::size_t, 3> fates{};
      std::array<std::size_t, kShortReasons + 1U> stops{};
      std::array<std::size_t, kMissedReasons> missed{};
      std::size_t count = 0U;
      for (std::size_t index = 0; index < mFindable.size(); ++index) {
        if (!select(index)) {
          continue;
        }
        ++count;
        ++fates[mFate[index]];
        if (mFate[index] == kShortFate) {
          ++stops[mStop[index]];
        } else if (mFate[index] == kNotFound) {
          ++missed[mReason[index]];
        }
      }
      const auto share = [&](std::size_t k) { return count != 0U ? 100. * static_cast<double>(k) / static_cast<double>(count) : 0.; };
      std::string line = std::format("    {}: {} ({:.2f} %) | {:.2f} {:.2f} {:.2f} |", name, count, 100. * static_cast<double>(count) / static_cast<double>(std::max<std::size_t>(mFindable.size(), 1U)), share(fates[0]), share(fates[1]), share(fates[2]));
      for (std::size_t reason = 0; reason < kShortReasons; ++reason) {
        if (stops[reason] != 0U) {
          line += std::format(" {} {:.2f};", kShortNames[reason], share(stops[reason]));
        }
      }
      line += " |";
      for (std::size_t reason = 0; reason < kMissedReasons; ++reason) {
        if (missed[reason] != 0U) {
          line += std::format(" {} {:.2f};", kMissedNames[reason], share(missed[reason]));
        }
      }
      std::cout << line << '\n';
    };
    for (std::size_t bin = 0; bin + 1U < kWorstEdges.size(); ++bin) {
      row(std::format("[{}, {})", kWorstEdges[bin], kWorstEdges[bin + 1U]), [&](std::size_t index) { return wholeOf(index) && static_cast<double>(worst[index]) >= kWorstEdges[bin] && static_cast<double>(worst[index]) < kWorstEdges[bin + 1U]; });
    }
    row("not every pair evaluated", [&](std::size_t index) { return !wholeOf(index); });
  }

  /// One CSV line: the step \a step's state, carried to its particle's row
  /// \a row as the step's update does, against the cluster and the MC
  /// crossing; the start's own stave's true kink.
  template <typename Material>
  void record(const char* kind, std::uint64_t particle, std::size_t iteration, const Walk& replay, std::size_t stepIndex, std::uint32_t row, Scalar chi2, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator)
  {
    const auto& step = replay.steps[stepIndex];
    const auto measurement = measurementAt(frame, row);
    const auto plane = geometry.planeOf(measurement);
    BoundTrackState predicted;
    CovarianceFactors predictedFactors;
    if (!follow::transportCandidate<false>(propagator, geometry, material, plane, PropagationDirection::kBackward, step.state, step.factors, predicted, predictedFactors, 1.F)) {
      return;
    }
    covarianceFromFactors(predictedFactors, predicted.cov);
    // The references: the clusters-only inward filter
    // over the road and the rows the proposal took before this step, carried
    // to the row's plane; and that prediction updated there with the beam
    // line, whose variance holds the scattering inside the target only --
    // the joint linear-Gaussian prediction with the beam point.
    // The rows taken before this step are the earlier steps' choices: the
    // walk takes at most one per step, and `raw` is the walk's result, short
    // of the rows its end dropped (the unconfirmed drop, the end refusals) that this step's
    // state still carries.
    std::vector<std::uint32_t> chain;
    for (auto i = std::min(stepIndex, replay.steps.size()); i-- > 0U;) {
      if (replay.steps[i].choice != kNoRow) {
        chain.push_back(replay.steps[i].choice);
      }
    }
    for (std::uint8_t i = 0; i < replay.roadCount; ++i) {
      chain.push_back(replay.road.index[i]);
    }
    const auto clusters = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = chain.data()}, static_cast<int>(chain.size()), kinkBudget(parameters));
    BoundTrackState noBeam;
    BoundTrackState joint;
    bool references = false;
    if (clusters.status == PropagationStatus::kOk) {
      CovarianceFactors noBeamFactors;
      references = follow::transportCandidate<false>(propagator, geometry, material, plane, PropagationDirection::kBackward, clusters.state, factorCovariance(clusters.state.cov), noBeam, noBeamFactors, 1.F);
      if (references) {
        covarianceFromFactors(noBeamFactors, noBeam.cov);
        joint = noBeam;
        references = updateWithBeamLine(propagator, parameters, plane, mTruth->rows[row].layer, pairPointAt(frame.measurements(), row), joint).ok;
      }
    }
    const auto chi2Of = [&](BoundTrackState state) {
      auto factors = factorCovariance(state.cov);
      return updateWithMeasurement(state, factors, measurement).chi2;
    };
    const Scalar chi2NoBeam = references ? chi2Of(noBeam) : -1.F;
    const Scalar chi2Joint = references ? chi2Of(joint) : -1.F;
    BoundTrackState start = step.state;
    covarianceFromFactors(step.factors, start.cov);
    const auto& truthRow = mTruth->rows[row];
    std::size_t onLayer = 0U;
    for (const auto at : mRowsOf.at(particle)) {
      onLayer += mTruth->rows[at].layer == truthRow.layer ? 1U : 0U;
    }
    const auto found = mTruth->particles.find(particle);
    const auto* segment = segmentOfRow(particle, frame, row);
    double c0 = 0.;
    double c1 = 0.;
    double q0 = 0.;
    double q1 = 0.;
    double qn0 = 0.;
    double qn1 = 0.;
    double qj0 = 0.;
    double qj1 = 0.;
    double dphi = 0.;
    double dtheta = 0.;
    double hitP = 0.;
    if (segment != nullptr) {
      const auto local = crossingOn(*segment, plane);
      c0 = static_cast<double>(measurement.loc0) - local[0];
      c1 = static_cast<double>(measurement.loc1) - local[1];
      q0 = static_cast<double>(predicted.loc0()) - local[0];
      q1 = static_cast<double>(predicted.loc1()) - local[1];
      qn0 = static_cast<double>(noBeam.loc0()) - local[0];
      qn1 = static_cast<double>(noBeam.loc1()) - local[1];
      qj0 = static_cast<double>(joint.loc0()) - local[0];
      qj1 = static_cast<double>(joint.loc1()) - local[1];
      const auto& p = segment->momentum;
      hitP = std::hypot(p[0], p[1], p[2]);
      dphi = std::remainder(static_cast<double>(predicted.phi()) - std::atan2(p[1], p[0]), 2. * std::numbers::pi);
      dtheta = static_cast<double>(predicted.theta()) - std::atan2(std::hypot(p[0], p[1]), p[2]);
    }
    // The start against the truth there, and across the start's own stave:
    // from entering the start's layer to entering the next one out.
    const auto* startSegment = segmentOfRow(particle, frame, step.from);
    double sdphi = 0.;
    double sdtheta = 0.;
    double staveTheta = 0.;
    double stavePhi = 0.;
    double staveTheta0 = 0.;
    double staveX2x0 = 0.;
    if (startSegment != nullptr) {
      const auto& ps = startSegment->momentum;
      sdphi = std::remainder(static_cast<double>(start.phi()) - std::atan2(ps[1], ps[0]), 2. * std::numbers::pi);
      sdtheta = static_cast<double>(start.theta()) - std::atan2(std::hypot(ps[0], ps[1]), ps[2]);
      // The stave beyond the start: the layer numbered one above (the
      // class's barrel layers consecutive from the beam outward).
      const auto next = ownRowOn(particle, static_cast<LayerID>(mTruth->rows[step.from].layer + 1U));
      const auto* nextSegment = next == kNoTrack ? nullptr : segmentOfRow(particle, frame, next);
      if (nextSegment != nullptr) {
        std::tie(staveTheta, stavePhi) = kinksBetween(*startSegment, *nextSegment, static_cast<double>(propagator.bzKiloGauss));
        std::tie(staveTheta0, staveX2x0) = theta0Between(material, propagator, *startSegment, *nextSegment, std::hypot(ps[0], ps[1], ps[2]));
      }
    }
    *mCsv << std::format(
      "{},{},{},{},{:.4f},{},{},{},{:.3f},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},"
      "{:.6g},{:.6g},{:.4f},{:.4f},{},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{:.6g},{},{:.3f},{:.6g},{:.6g},{:.6g},{:.6g},{:.3f},{:.6g},{:.6g},"
      "{:.6g},{:.6g}\n",
      kind, particle, iteration, static_cast<int>(truthRow.layer), found != mTruth->particles.end() ? found->second.pt : 0.F, found != mTruth->particles.end() ? found->second.pdg : 0, static_cast<int>(truthRow.count), onLayer, chi2, measurement.loc0 - predicted.loc0(), measurement.loc1 - predicted.loc1(), measurement.varLoc0, measurement.covLoc01, measurement.varLoc1, predicted.covAt(kLoc0, kLoc0), predicted.covAt(kLoc0, kLoc1), predicted.covAt(kLoc1, kLoc1), segment != nullptr ? 1 : 0, c0, c1, q0, q1, dphi, dtheta, predicted.covAt(kPhi, kPhi), predicted.covAt(kTheta, kTheta), hitP, 1. / std::abs(static_cast<double>(predicted.qOverP())), startSegment != nullptr ? 1 : 0, sdphi, sdtheta, start.covAt(kPhi, kPhi), start.covAt(kTheta, kTheta), staveTheta, stavePhi, staveTheta0, staveX2x0, references ? 1 : 0, chi2NoBeam, qn0, qn1, noBeam.covAt(kLoc0, kLoc0), noBeam.covAt(kLoc1, kLoc1), chi2Joint, qj0, qj1, joint.covAt(kLoc0, kLoc0), joint.covAt(kLoc1, kLoc1));
  }

  struct TrueKink {
    std::uint32_t particle = 0U; ///< its index in mFindable
    LayerID layer = 0;           ///< the inner of the two crossings
    LayerID to = 0;              ///< the outer: the next layer of the mask the particle crossed (ML 4 -> C1 is one)
    float p = 0.F;               ///< [GeV], entering it
    float theta = 0.F;           ///< the theta kink over theta0
    float phi = 0.F;             ///< the bending plane's
    float x2x0 = 0.F;
  };

  const study::TrackingTruth* mTruth = nullptr;
  std::uint64_t mMask = 0U;
  double mMinPtGeV = 0.;
  std::vector<std::uint64_t> mRoadLayers;
  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> mRowsOf;
  std::unordered_map<std::uint64_t, std::uint32_t> mIndex;
  std::vector<std::uint64_t> mFindable;
  std::vector<std::uint8_t> mFate;
  std::vector<std::uint8_t> mReason;       ///< a Missed, for those not found
  std::vector<std::uint8_t> mStop;         ///< a Short, for the short, from their iteration
  std::vector<LayerID> mStopLayer;         ///< for the short with a stop decision, the layer it missed
  std::vector<std::uint8_t> mStopTaken;    ///< for the short whose own row passed: `takenKind`
  std::vector<std::uint8_t> mOutStop;      ///< for the short of an outer layer: the outward step's Short
  std::vector<LayerID> mOutStopLayer;      ///< for those: the first outer layer it missed
  std::vector<std::uint8_t> mOutStopTaken; ///< for those whose own row passed: `takenKind`
  std::vector<std::uint8_t> mInNext;       ///< for an unconfirmed drop: the next step's outcome on the own row (`afterDrop`), inward
  std::vector<std::uint8_t> mOutNext;      ///< the same outward
  std::vector<std::uint8_t> mInDropKink;   ///< for an unconfirmed drop: the particle's true kink over that next step (<= 3 / 3-5 / > 5 theta0 / no hits), inward
  std::vector<std::uint8_t> mOutDropKink;  ///< the same outward
  std::vector<std::uint8_t> mOutStepKink;  ///< for the short of an outer layer: the particle's true kink over the step from the layer below the gap (binned as above)
  std::vector<std::uint8_t> mInStepKink;   ///< for an inward unconfirmed drop: the true kink over the step that took the dropped row
  std::vector<std::uint64_t> mLayersHeld;  ///< the good track's layers, a bit each (0: no good track), for `writeFates`
  std::unique_ptr<std::ofstream> mWalkCsv; ///< --truth-walk: one line per step of the truth walks
  /// The truth walks' predictions against the cluster, unselected, for the
  /// robust summary: [inward, outward][target layer][pT bin].
  struct WalkPulls {
    std::vector<float> pull0;
    std::vector<float> pull1;
    std::size_t overGate = 0U; ///< the update's chi2 above the gate
    std::size_t fails = 0U;    ///< the transport failed
  };
  static constexpr std::array<double, 4> kWalkPtEdges{0.2, 0.3, 0.6, 1.};
  std::array<std::array<std::array<WalkPulls, kWalkPtEdges.size() + 1U>, kLayers>, 3> mWalkPulls{};
  static constexpr std::array<const char*, 3> kWalkDirectionNames{"inward", "outward", "barrel to the rest"};
  Scalar mWalkGate = 0.F;                            ///< the gate the truth walks were held to
  std::size_t mDebugDrops = 0U;                      ///< unconfirmed drops whose next step saw the own row passing, dumped
  std::vector<Scalar> mInChi2;                       ///< the own row's chi2 at the deciding step, inward
  std::vector<Scalar> mOutChi2;                      ///< the same outward
  std::vector<std::uint8_t> mOutRefitWrong;          ///< refit-chi2 stops: the outward rows kept not the particle's (0, 1, 2+), +3 with a kink step
  std::vector<std::uint8_t> mOutRefitRows;           ///< and how many outward rows the refit took
  std::vector<Scalar> mOutRefitRatio;                ///< and its chi2 over the cut
  std::map<std::string, std::size_t> mRefitFailures; ///< refit-failed stops: status / pass / layer
  std::map<std::string, std::size_t> mRefitReplay;   ///< and their failing step replayed (`FailedStep`)
  /// The extended refit's inward increments by layer: the sum,
  /// the count, those above 9.21 (2 dof, 1 %), and how often the layer's is
  /// the track's largest.
  struct Increments {
    std::array<double, kLayers> sum{};
    std::array<std::size_t, kLayers> count{};
    std::array<std::size_t, kLayers> above{};
    std::array<std::size_t, kLayers> largest{};
    std::size_t tracks = 0U;
  };
  Increments mRefusedIncrements; ///< refit-chi2 stops, the rows all the particle's own
  Increments mKeptIncrements;    ///< complete tracks the iteration followed outward, refit kept
  void addIncrements(Increments& into, const Walk& replay)
  {
    ++into.tracks;
    std::size_t top = 0U;
    for (std::size_t i = 0; i < replay.refitCount; ++i) {
      const auto layer = replay.refitLayer[i];
      const auto increment = static_cast<double>(replay.refitIncrement[i]);
      into.sum[layer] += increment;
      ++into.count[layer];
      into.above[layer] += increment > 9.21 ? 1U : 0U;
      top = replay.refitIncrement[i] > replay.refitIncrement[top] ? i : top;
    }
    if (replay.refitCount != 0U) {
      ++into.largest[replay.refitLayer[top]];
    }
  }
  std::vector<std::uint8_t> mOutRefitKink;                            ///< and the largest true kink over the rows (<= 3 / 3-5 / > 5 theta0 / no hits)
  std::vector<std::uint8_t> mOutHolder;                               ///< for those: who holds the own row on that layer at the end (`kHolderNames`)
  std::vector<Scalar> mInNextChi2;                                    ///< for an unconfirmed drop: the own next row's chi2, inward
  std::vector<Scalar> mOutNextChi2;                                   ///< the same outward
  std::size_t mOutReplayed = 0U;                                      ///< outward walks replayed
  Scalar mLastGate = 0.F;                                             ///< the last iteration's update gate, clusterChi2Max
  std::size_t mOutAgree = 0U;                                         ///< of them, agreeing with the iteration's kept rows
  std::array<std::size_t, 65> mLayersFound{};                         ///< good tracks of the findable, by the mask's layers they hold
  std::array<std::array<std::size_t, 8>, 2> mWindowAudit{};           ///< window, not queried; by windowAudit's flags
  std::array<std::array<std::size_t, 8>, 2> mOutWindowAudit{};        ///< the same of the outward stops
  std::array<std::vector<std::array<float, 4>>, 2> mOutWindowOffsets; ///< their rows' offsets (`windowAudit`)
  std::array<TripletTargetBounds, kLayers> mBounds{};                 ///< each layer's largest variances in the frame
  mutable std::size_t mAuditDecisions = 0U;                           ///< the graph audit's triplet and neighbour decisions
  mutable std::map<std::string, std::size_t> mCurvatureTriplets;      ///< the graph audit's curvature refusals by the triplet's layers
  mutable std::map<std::string, std::size_t> mLinkRangeRefusals;      ///< the graph audit's link-range refusals by the link and the distance outside [cm], 1-cm bins
  VertexSeed mBeamLine;                                               ///< the beam line the pair cut is held to (`setBeamLine`)
  mutable std::map<std::string, std::size_t> mPairRefusals;           ///< the graph audit's pair refusals by the test, the link and the true kink
  double mAuditFloor = 0.;                                            ///< the iteration's pT floor, for the pair refusals' split
  std::vector<std::uint32_t> mTrackOffsets;                           ///< per iteration observed, its first track's index
  std::map<std::string, std::size_t> mFakeHolders;                    ///< the not-found whose road rows a fake holds, by the fake's iteration
  std::map<std::string, std::size_t> mRoadFitRefusals;                ///< the road fit's refusals of true roads, binned (`describeRefusedRoad`)
  std::vector<std::uint8_t> mMissedKink;                              ///< for a particle not found: the largest true kink over its road's layers (binned as `kinkBin`)

  mutable std::size_t mDebugTriplets = 0U;      ///< of them, dumped
  mutable std::size_t mAuditDisagreements = 0U; ///< of them, those the stage's own cut decided otherwise
  std::vector<std::uint8_t> mGraph;             ///< a Graph, for those not found for kNoRoad or kRoadTime
  struct WalkSum {
    std::size_t n = 0U;
    double s0 = 0.;
    double s1 = 0.;
    std::size_t out0 = 0U;
    std::size_t out1 = 0U;
  };
  /// truthWalk's pulls: [all, quiet, quiet and prompt][follow, clusters only, joint][target layer]
  std::array<std::array<std::array<WalkSum, kLayers>, 3>, 3> mWalk{};
  bool mWalkDone = false;
  std::size_t mWalkIteration = 0U; ///< the iteration the truth walk runs on (`observe`)
  std::vector<std::uint32_t> mMaxWidthPS;
  SensorCells<Det> mCells; ///< the sensors by place, which the follow's hole test reads
  std::size_t mReplayed = 0U;
  std::size_t mAgree = 0U;
  /// A replayed track's fidelity: as the iteration made it; a refused refit
  /// where the iteration kept rows; otherwise not.
  enum Fidelity : std::uint8_t { kUnfaithful,
                                 kFaithful,
                                 kRefusedKept };
  struct FidelityCount {
    std::size_t replayed = 0U;
    std::size_t agree = 0U;
    std::size_t refusedKept = 0U;
  };
  std::vector<std::array<FidelityCount, 2>> mFidelity;        ///< per iteration, its own tracks' replays and the earlier tracks' (`printFidelity`)
  std::vector<std::string> mTransitionLines;                  ///< fakeTransitions', one per iteration
  std::vector<std::uint8_t> mFollowUsedNow;                   ///< the mask the iteration's follow read
  std::vector<std::string> mProgressLines;                    ///< progress', one per iteration
  std::vector<double> mFloors;                                ///< each iteration's pT floor [GeV]
  std::array<std::array<std::size_t, 2>, 21> mReach{};        ///< physical primaries in |eta| < 1 by 50 MeV of pT: all, on every layer
  double mLowestPt = std::numeric_limits<double>::infinity(); ///< the lowest pT of one on every layer
  std::vector<Snapshot> mSnapshot;                            ///< the tracks as the last iteration observed left them
  std::vector<std::uint8_t> mSnapshotUsed;                    ///< and its mask
  /// A kink step's refusals of the particle's row by its outcome
  /// (StepOutcome, then not visited); its unconfirmed rows by what became
  /// of their witness, the particle's next row [0] or the beam line [1]
  /// (StepOutcome, not visited, no own next row).
  static constexpr std::size_t kNotVisited = 6U;
  static constexpr std::size_t kNoOwnNext = 7U;
  std::vector<std::uint8_t> mStopDetail; ///< for kKinkRefused and kKinkUnconfirmed: the outcome, plus 8 for the beam line
  std::unique_ptr<std::ofstream> mCsv;
  std::unique_ptr<study::HitFile> mHits;
  std::vector<std::vector<int>> mMothers;              ///< `openKine`'s: per MC event, each track's mother
  std::map<std::string, ClutterTally> mClutter;        ///< `clutterStudy`'s, by iteration, pT, walk and layers
  std::vector<DroppedWalks> mDroppedWalks;             ///< the observed iteration's dropped records' walks (`replayDropped`, for `clutterStudy`)
  std::map<std::string, DecisionTally> mDecisions;     ///< `decisionStudy`'s, by population, iteration, pT, walk and layers
  std::map<std::string, BlockerTally> mBlockers;       ///< `blockerStudy`'s, the last observed iteration's
  std::vector<std::vector<std::uint32_t>> mLayerLower; ///< per layer: its rows' time boxes' starts, sorted (`layerOccupancy`)
  std::vector<std::vector<std::uint32_t>> mLayerUpper; ///< and their ends
  std::vector<double> mLayerArea;                      ///< per layer: its sensors' area [cm2]
  /// Each MC track's segments, by `segmentKey`.
  std::unordered_map<std::uint64_t, std::vector<const study::HitSegment*>> mSegmentsOf;
  bool mPathOrdered = false;
  /// Per findable: each row's leg in path order (`pathOrder`), parallel to its rows; whether its rows are in path order.
  std::vector<std::vector<std::uint8_t>> mLegs;
  std::vector<std::uint8_t> mTimed;
  /// The path order's outcome, [not in the class, the class][timed as the distance had it, reordered, kept, first row on a layer changed].
  std::array<std::array<std::size_t, 4>, 2> mPathCounts{};
  std::string mPathLines;                                      ///< `pathReport`'s, the last iteration's
  std::vector<std::vector<std::uint8_t>> mSeedStage;           ///< per iteration, per findable: a SeedStage
  std::vector<std::uint8_t> mFurthest;                         ///< per findable: its furthest SeedStage over the iterations
  std::vector<std::string> mRecordLoss;                        ///< per findable: `recordLosses`' verdict, its last record never selected
  bool mKeepBeamLine = false;                                  ///< `keepBeamLine`
  std::vector<BeamLineRow> mBeamLineRows;                      ///< `beamLineStudy`, every iteration
  bool mKeepFollowLoss = false;                                ///< `keepFollowLoss`
  std::vector<FollowLossRow> mFollowLossRows;                  ///< `followLoss`, every iteration
  std::unordered_map<std::uint32_t, std::string> mLastDropWhy; ///< per findable: the last verdict of `replayDropped`
  std::vector<std::uint32_t> mMinClusters;                     ///< per iteration: its seeds' length
  std::vector<std::map<std::string, std::size_t>> mDropped;    ///< per iteration: `replayDropped`'s lines
  /// `inwardStudy`'s: per kind and state, the truth
  /// comparisons; the passes replayed and those the replay did not
  /// reproduce; the failing passes from other seeds; the failing steps'
  /// roots.
  static constexpr std::size_t kInwardStepSlots = 6U;
  struct InwardSlot {
    std::vector<double> ratio, pullQOverP, dphi, pullPhi, relSigma;
  };
  std::array<std::array<InwardSlot, kInwardStepSlots>, kInwardKinds> mInwardSteps{};
  std::array<std::size_t, kInwardKinds> mInwardCount{};
  std::array<std::size_t, kInwardKinds> mInwardMismatch{};
  std::array<std::array<std::size_t, 4>, kSeedNames.size()> mSeedOutcomes{};
  std::map<std::string, std::size_t> mRoots;
  /// `ruleStudy`'s, per kind, rule and seed.
  std::array<std::array<std::array<RuleTally, kRuleSeedNames.size()>, kCrossingRules.size()>, kInwardKinds> mRuleTallies{};
  std::unique_ptr<std::ofstream> mInwardCsv;
  std::unique_ptr<std::ofstream> mEndCsv;  ///< `endDecisions`' CSV
  std::unique_ptr<std::ofstream> mEditCsv; ///< and its end edits'
  std::array<std::array<std::size_t, 2>, kInwardKinds> mChargeAgree{};
  /// Per kind and step, the step's crossing [disc, the solve's the particle's, the other the particle's, none].
  std::array<std::array<std::array<std::size_t, 4>, kInwardStepSlots>, kInwardKinds> mStepRoots{};
  std::map<std::string, std::size_t> mStepRootDetail;
  bool mKinksDone = false;
  std::vector<TrueKink> mKinks;
  std::vector<std::uint32_t> mKinkBegin;                                                                  ///< particle index -> its first kink in mKinks (grouped by particle), one past the end last
  std::unique_ptr<o2::gtrk::detail::TaskArena> mArena = std::make_unique<o2::gtrk::detail::TaskArena>(0); ///< the replays' and the truth walks' threads
  /// Timings of the phases of `observe`, printed per iteration [s].
  struct Phases {
    double kinks = 0.;
    double walk = 0.;
    double replay = 0.;
    double diagnosis = 0.;
  };
  std::vector<std::uint64_t> mEvaluated;
  std::vector<std::uint64_t> mCrossed; ///< per findable, the layers with a crossing (one row per layer, the first)
  std::vector<bool> mRepeated;
  std::size_t mMissingSegment = 0U;
  std::size_t mZeroMaterial = 0U;
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_COMPLETEANATOMY_H_
