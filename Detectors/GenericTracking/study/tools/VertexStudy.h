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

/// \file VertexStudy.h
/// \brief The vertexer on one iteration's triplets against the MC
///        collisions, and its replay.
///
/// The study observes one iteration -- the vertexing pass, links within the
/// innermost layers -- and keeps each distinct triplet: its rows, its z0 from
/// a straight line in r-z weighted by the points' z variances, and its frame,
/// the middle measurement's time. It runs the vertexer (`Vertexer`,
/// VertexStage.h) on every slice's triplets and measures each triplet as the
/// vertexer does, without its cuts; it replays the vertexer's rounds on them
/// (`replayAll`), checked to be the vertexer bit for bit at its settings,
/// with the variants the vertexer does not run (the cuts in sigmas of the
/// triplet's own errors, an earlier pass's vertices first).
/// It prints the triplets' pulls in units of their own errors
/// (`describePulls`) and the beam point the claimed triplets agree on
/// (`estimateBeamPoint`), and writes the vertexer's vertices with their
/// makeup (`writeStageVertices`) and the replay's lists (`writeReplay`).
///
/// A triplet is true when one particle carries all three measurements; it
/// then belongs to that particle's collision when the particle is a primary.
/// A vertex belongs to the collision most of its candidates belong to;
/// background (fake or secondary triplets) as the majority makes it a fake. A
/// collision is found when a vertex belongs to it; a second vertex of the
/// same collision is a split; a collision whose own candidates went mostly to
/// another collision's vertex is merged. What the tracking loses with a
/// missed collision is counted as its findable primaries, the tracking QA's
/// (`isFindable`: physical primaries on a chain of four layers).

#ifndef GTRK_STUDY_VERTEXSTUDY_H_
#define GTRK_STUDY_VERTEXSTUDY_H_

#include "GenericTracking/Tracking/ShapeCut.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Study/Validation.h"
#include "GenericTracking/Tracking/VertexStage.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/VertexQA.h"
#include "GenericTracking/Study/Statistics.h"
#include "GenericTracking/Study/TrackComparison.h"

#include <algorithm>
#include <bit>
#include <array>
#include <cmath>
#include <cstdint>
#include <format>
#include <fstream>
#include <iostream>
#include <chrono>
#include <map>
#include <memory>
#include <limits>
#include <numbers>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

/// The beam point a timeframe's triplets agree on: the least-squares solution of their BeamTerms, in the frame the
/// measurements are in, so that it is the offset of the true beam line from
/// the one the rows were placed about. Online, with no calibrated MeanVertex
/// for the run yet, it is the beam position; offline it is the check on the
/// object, and the watch on a drift.
struct BeamEstimate {
  Scalar x = 0.F;      ///< [cm]
  Scalar y = 0.F;      ///< [cm]
  Scalar sigmaX = 0.F; ///< the solution's error from the residuals' spread [cm]
  Scalar sigmaY = 0.F;
  Scalar rmsCM = 0.F;     ///< the residuals' spread of the terms used [cm]
  std::size_t used = 0U;  ///< the terms the last pass kept
  std::size_t terms = 0U; ///< the terms offered
  bool valid = false;     ///< enough terms in enough directions
};

/// The beam point from \a terms, each counted once: a 2x2 least squares of
/// n . p = d, made robust by three passes -- every term within
/// \a maxDistanceCM of the axis (the rest debris, `kMaxBeamDistanceCM`), then
/// the residuals within 1 mm of the first solution, then within three times
/// their spread (at least 50 um) of the second -- so that the fake triplets,
/// whose distances are broad, and a beam line off by a real offset of a
/// millimetre both leave the answer alone. Invalid with fewer than 10 terms
/// or directions that do not span the plane.
[[nodiscard]] inline BeamEstimate estimateBeamPoint(std::span<const BeamTerm> terms, Scalar maxDistanceCM)
{
  BeamEstimate estimate;
  estimate.terms = terms.size();
  std::vector<std::uint8_t> keep(terms.size(), 0U);
  for (std::size_t i = 0; i < terms.size(); ++i) {
    keep[i] = terms[i].valid && std::abs(terms[i].d) < maxDistanceCM ? 1U : 0U;
  }
  Scalar x = 0.F;
  Scalar y = 0.F;
  Scalar rms = 0.F;
  std::size_t used = 0U;
  Scalar sxx = 0.F;
  Scalar sxy = 0.F;
  Scalar syy = 0.F;
  Scalar det = 0.F;
  for (int pass = 0; pass < 3; ++pass) {
    sxx = sxy = syy = 0.F;
    Scalar vx = 0.F;
    Scalar vy = 0.F;
    Scalar sw = 0.F;
    used = 0U;
    for (std::size_t i = 0; i < terms.size(); ++i) {
      if (keep[i] == 0U) {
        continue;
      }
      const auto& term = terms[i];
      sxx += term.nx * term.nx;
      sxy += term.nx * term.ny;
      syy += term.ny * term.ny;
      vx += term.nx * term.d;
      vy += term.ny * term.d;
      sw += 1.F;
      ++used;
    }
    det = (sxx * syy) - (sxy * sxy);
    // The directions must span the plane: the normal matrix's determinant at
    // least a tenth of its value for directions spread evenly, (trace/2)^2.
    const Scalar trace = sxx + syy;
    if (used < 10U || !(det > 0.F) || det < 0.1F * trace * trace / 4.F) {
      estimate.used = used;
      return estimate;
    }
    x = ((syy * vx) - (sxy * vy)) / det;
    y = ((sxx * vy) - (sxy * vx)) / det;
    Scalar sumSquares = 0.F;
    for (std::size_t i = 0; i < terms.size(); ++i) {
      if (keep[i] == 0U) {
        continue;
      }
      const Scalar residual = terms[i].d - (terms[i].nx * x) - (terms[i].ny * y);
      sumSquares += residual * residual;
    }
    rms = sw > 0.F ? std::sqrt(sumSquares / sw) : 0.F;
    if (pass == 2) {
      break;
    }
    const Scalar cut = pass == 0 ? 0.1F : std::max(3.F * rms, 0.005F);
    for (std::size_t i = 0; i < terms.size(); ++i) {
      if (keep[i] == 0U) {
        continue;
      }
      const Scalar residual = terms[i].d - (terms[i].nx * x) - (terms[i].ny * y);
      keep[i] = std::abs(residual) < cut ? 1U : 0U;
    }
  }
  estimate.x = x;
  estimate.y = y;
  estimate.rmsCM = rms;
  estimate.used = used;
  // The covariance of the solution: rms^2 times the inverse normal matrix.
  estimate.sigmaX = rms * std::sqrt(syy / det);
  estimate.sigmaY = rms * std::sqrt(sxx / det);
  estimate.valid = true;
  return estimate;
}

/// The beam estimate ignores a triplet whose circle passes farther than
/// this from the axis [cm]: debris, no track from a collision.
inline constexpr Scalar kMaxBeamDistanceCM = 0.5F;

template <DetectorLayout Det, typename Controller>
class VertexStudy
{
 public:
  /// The vertexer's settings as the study's tables print them and its own
  /// candidates are cut with (the makeup's reuse shares), from the configuration, in double and integer micrometres; the vertexer's
  /// (`VertexerSettings`) narrow them to `Scalar` and `int32_t`.
  struct Settings {
    std::int64_t windowUM = 1000; ///< w: the window a peak is counted in
    std::int64_t claimUM = 1000;  ///< the half-width around the peak's centre whose free candidates the vertex takes
    std::uint32_t minContributors = 3U;
    double significance = 0.; ///< k: the peak must exceed the sidebands' mean b by k sqrt(b)
    double maxChi2RZ = 0.;    ///< candidates above it do not count; 0 for no cut
    double maxDcaUM = 0.;     ///< the same for the transverse distance to the axis [um]
  };

  /// The study of \a iteration's triplets, its pass's pT floor \a
  /// passMinPtGeV, against \a rows and \a particles, the collisions' findable
  /// primaries by the tracking QA's \a comparison; the vertexer run on them
  /// with \a stage on \a threads, the tables at \a table (the same settings
  /// from the configuration).
  VertexStudy(const TimeFrame<Det>& frame, std::span<const study::TruthRow> rows, const study::ParticleTable& particles, std::size_t iteration,
              double passMinPtGeV, const study::ComparisonSettings& comparison, const VertexerSettings& stage, const Settings& table, int threads)
    : mFrame(&frame),
      mRows(rows),
      mParticles(&particles),
      mIteration(iteration),
      mPassMinPtGeV(passMinPtGeV),
      mComparison(comparison),
      mStage(stage, VertexConstraint{}),
      mTable(table),
      mStagePool(std::make_unique<BoundedMemoryResource>()),
      mArena(std::make_unique<o2::gtrk::detail::TaskArena>(threads))
  {
  }

  void observeSlice(const Controller& tracker, const typename Controller::SliceObservation& slice)
  {
    if (slice.iteration != mIteration) {
      return;
    }
    mPass = tracker.iteration(slice.iteration);
    mBeamVarTransverse = tracker.beamLine().varTransverse;
    const auto& result = *slice.result;
    const auto tracklets = result.tracklets->view();
    const auto triplets = result.triplets->view();
    ++mSlices;
    const auto columns = slice.frame->measurements();
    for (std::size_t t = 0; t < result.tracklets->size(); ++t) {
      ++mTrackletsByLink[tracklets.link[t]];
      mTrackletsShaped[tracklets.link[t]] += pairShaped(columns, tracklets.source[t], tracklets.target[t]) ? 1U : 0U;
    }
    // The vertexer itself on this slice's triplets.
    const auto start = std::chrono::steady_clock::now();
    const ScatteringErrors<Det> errors{.table = &mPass->scattering, .cut = mPass->pairCut, .varTransverse = mBeamVarTransverse};
    // Its work arrays on a pool of their own, for their peak.
    mStagePool->resetPeakMemory();
    const auto found = mStage.find(slice.frame->measurements(), tracklets, triplets, result.triplets->size(), mArena.get(), true, errors, mStagePool.get(),
                                   ownedFrames(slice.frame->grid().geometry(), slice.slice->coreBlocks));
    mStageSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    mStagePoolPeak = std::max(mStagePoolPeak, mStagePool->getPeakMemory());
    // Each triplet as the vertexer measures it, without its cuts, for the
    // replay (`replay`): its status, chi2, z, dca, weight and box.
    VertexerSettings open = mStage.settings();
    open.maxChi2RZ = 0.F;
    open.maxDcaCM = 0.F;
    const auto measurements = slice.frame->measurements();
    // The slice's triplets in its order and the frames it owns: the
    // vertexer deduplicates per slice, ties in triplet order.
    auto& sliceTriplets = mSliceTriplets.emplace_back(ownedFrames(slice.frame->grid().geometry(), slice.slice->coreBlocks), std::vector<Triple>{});
    sliceTriplets.second.reserve(result.triplets->size());
    // Each triplet's key and beam term, for the beam point of its claimed
    // candidates in the vertexer's order: by key, ties in triplet order.
    std::vector<std::uint64_t> keys(result.triplets->size(), 0U);
    std::vector<BeamTerm> terms(result.triplets->size());
    for (std::size_t t = 0; t < result.triplets->size(); ++t) {
      VertexCandidateFit fit;
      StageFit stageFit;
      BeamTerm term;
      stageFit.status =
        measureCandidate(measurements, tracklets, triplets, static_cast<std::uint32_t>(t), open, errors, fit, stageFit.weight, stageFit.box, &term);
      stageFit.chi2RZ = fit.chi2RZ;
      stageFit.zUM = stageFit.status == VertexCandidateStatus::kAccepted ? vertexKeyZUM(fit.key) : 0;
      stageFit.dcaCM = term.d;
      keys[t] = fit.key;
      terms[t] = term;
      const auto first = triplets.first[t];
      mStageFits.emplace_back(Triple{tracklets.source[first], tracklets.target[first], tracklets.target[triplets.second[t]]}, stageFit);
      sliceTriplets.second.push_back(mStageFits.back().first);
    }
    mStageCandidates += found.candidates;
    for (std::size_t t = 0; t < result.triplets->size(); ++t) {
      const auto first = triplets.first[t];
      const auto second = triplets.second[t];
      const auto vertex = found.tripletVertex[t];
      if (vertex == VertexerResult::kNotOwned) {
        continue; // another slice's frame
      }
      // A slice can hold one triple twice (two triplets on the same
      // rows); the vertexer keeps one of them, the other deduplicated
      // (-2): the claim is the kept one's.
      const auto claim = vertex < 0 ? std::int64_t{vertex} : vertex + static_cast<std::int64_t>(mStageVertices.size());
      const auto [at, inserted] = mStageClaims.try_emplace({tracklets.source[first], tracklets.target[first], tracklets.target[second]}, claim);
      if (!inserted) {
        ++mStageRepeats;
        at->second = std::max(at->second, claim);
      }
    }
    for (std::size_t v = 0; v < found.vertices.size(); ++v) {
      mStageVertices.push_back(found.vertices[v]);
      mStageContributors.push_back(found.contributors[v]);
    }
    mStageRounds = std::max(mStageRounds, found.rounds);
    std::vector<std::uint32_t> claimed;
    for (std::size_t t = 0; t < result.triplets->size(); ++t) {
      if (found.tripletVertex[t] >= 0) {
        claimed.push_back(static_cast<std::uint32_t>(t));
      }
    }
    std::ranges::stable_sort(claimed, {}, [&](std::uint32_t t) { return keys[t]; });
    std::vector<BeamTerm> claimedTerms;
    claimedTerms.reserve(claimed.size());
    for (const auto t : claimed) {
      claimedTerms.push_back(terms[t]);
    }
    const auto beam = estimateBeamPoint(claimedTerms, kMaxBeamDistanceCM);
    if (beam.terms > mStageBeam.terms) {
      mStageBeam = beam; // the slice with the most to say; one slice in the usual run
    }
  }

  /// The beam line the rows are about, for the estimate's absolute position.
  void setBeamLine(double x, double y)
  {
    mBeamX = x;
    mBeamY = y;
  }

  /// An earlier pass's vertices: the `stage` rows of a
  /// replay's <prefix>-vertices.csv, listed ones, by their frame's time.
  /// With them the replay adds the variants `prior-*`, which start each
  /// frame from them: their claims taken, their triplets within the claim
  /// owning their rows, before the first round.
  void setPriorVertices(const std::string& path)
  {
    std::ifstream in(path);
    if (!in) {
      throw std::runtime_error("cannot read " + path);
    }
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
      std::vector<std::string> fields;
      for (const auto part : line | std::views::split(',')) {
        fields.emplace_back(part.begin(), part.end());
      }
      // variant, id, listed, framePS, round, windowStartUM, centreUM, zUM, ...
      if (fields.size() < 8U || fields[0] != "stage" || std::stoll(fields[2]) < 0) {
        continue;
      }
      mPriorVertices[static_cast<std::uint32_t>(std::stoul(fields[3]))].push_back(static_cast<Scalar>(std::stod(fields[7])));
    }
  }

  /// Everything, as lines: the pass's counts, the candidates against the
  /// truth, the triplets' pulls, the vertexer's vertices against the
  /// collisions, the replay, and the vertexer's fakes and missed collisions.
  [[nodiscard]] std::vector<std::string> report()
  {
    std::vector<std::string> lines;
    // Every slice's triplets, then the distinct ones: a triple twice in a
    // slice or in two slices' halos is one candidate.
    const auto observed = mStageFits.size();
    std::ranges::sort(mStageFits, {}, [](const auto& entry) { return entry.first; });
    const auto repeated = std::ranges::unique(mStageFits, {}, [](const auto& entry) { return entry.first; });
    mStageFits.erase(repeated.begin(), repeated.end());
    std::string links;
    for (const auto& [link, count] : mTrackletsByLink) {
      links += std::format(" {}->{} {}", Det::layerLinks[link].from, Det::layerLinks[link].to, count);
      links += std::format(" (shape {})", mTrackletsShaped[link]);
    }
    lines.push_back(std::format("vertex study: iteration {}, {} slices; tracklets by link:{}; triplets {} ({} distinct)", mIteration, mSlices, links, observed,
                                mStageFits.size()));
    buildCollisions();
    buildCandidates();
    describeCandidates(lines);
    describePulls(lines);
    const auto stage = stageScore();
    lines.emplace_back("the vertexer (VertexStage.h) at the table's settings:");
    lines.push_back(describe(stage, mTable) + std::format(" | rounds {}", mStageRounds));
    lines.push_back(std::format("the vertexer took {:.4f} s for {} candidates, {} threads; its work arrays' peak {:.1f} MB", mStageSeconds, mStageCandidates,
                                mArena->maxConcurrency(), static_cast<double>(mStagePoolPeak) * 1.e-6));
    const auto& beam = mStageBeam;
    lines.push_back(beam.valid
                      ? std::format("the beam point from the triplets' circles: ({:+.1f}, {:+.1f}) um off the rows' beam line, +- ({:.1f}, "
                                    "{:.1f}) um, residuals' rms {:.0f} um, {} of {} terms; the rows about ({:.4f}, {:.4f}) cm, so the beam at ({:.4f}, "
                                    "{:.4f}) cm",
                                    static_cast<double>(beam.x) * 1.e4, static_cast<double>(beam.y) * 1.e4, static_cast<double>(beam.sigmaX) * 1.e4,
                                    static_cast<double>(beam.sigmaY) * 1.e4, static_cast<double>(beam.rmsCM) * 1.e4, beam.used, beam.terms, mBeamX, mBeamY,
                                    mBeamX + static_cast<double>(beam.x), mBeamY + static_cast<double>(beam.y))
                      : std::format("the beam point from the triplets' circles: no estimate ({} of {} terms usable)", beam.used, beam.terms));
    // The vertexer's list (its debris rule included) by multiplicity; its
    // diagnostics after the replay.
    const auto efficiency = efficiencyByMultiplicity(mTable, stage, "the vertexer's ");
    lines.push_back(efficiency.table);
    replayAll(lines);
    lines.emplace_back("the vertexer's fakes and missed collisions:");
    lines.insert(lines.end(), efficiency.diagnostics.begin(), efficiency.diagnostics.end());
    return lines;
  }

  /// The replay's vertices, collisions and drops, one CSV
  /// each, \a prefix-vertices.csv, -collisions.csv and -dropped.csv, every
  /// variant's rows with its name, for tools/vertexReplay.py, and the first
  /// variant's owned rows as their clusters, \a prefix-owned.txt. The drops
  /// are the true triplets an ownership dropped whose collision is not the
  /// owner's majority: what a real neighbour loses.
  void writeReplay(const std::string& prefix) const
  {
    std::ofstream vertices(prefix + "-vertices.csv");
    std::ofstream collisions(prefix + "-collisions.csv");
    std::ofstream drops(prefix + "-dropped.csv");
    if (!vertices || !collisions || !drops) {
      throw std::runtime_error(std::format("cannot write the replay's lists {}-*.csv", prefix));
    }
    vertices << "variant,id,listed,framePS,round,windowStartUM,centreUM,zUM,window,sidebands,hidden,contributors,selected,dropped,reuse,"
                "collision,trueZUM,findable,reference,startPS,widthPS,denseShare,softShare,noParticleShare\n";
    collisions << "variant,collision,trueZUM,framePS,inner,findable,found,vertexZUM,contributors,round,accepted,kept,core,chi2,dca,"
                  "degenerate,deduplicated,free,own,other,fake,dropped,dense,denseFindable,denseDzUM,denseContributors,vetoedWindow,vetoedSidebands\n";
    // The vertexer's owned rows (the first run's) as their clusters, one per
    // line -- rows move with the grid, which a pass's floor rebins -- for a
    // later pass to hold (trackBench --hold-clusters).
    std::ofstream owned(prefix + "-owned.txt");
    if (!owned) {
      throw std::runtime_error(std::format("cannot write {}-owned.txt", prefix));
    }
    const auto columns = mFrame->measurements();
    for (const auto row : mReplays.front().owned) {
      owned << columns.clusterIndex[row] << '\n';
    }
    drops << "variant,collision,findable,zUM,owner,ownerCollision,ownerZUM,source,sourceCollision,sourceState,sourceZUM,row,rowCollision\n";
    constexpr std::array<std::string_view, kReplayStates> kStateNames{"chi2", "dca", "degenerate", "deduplicated", "kept",
                                                                      "own", "other", "fake", "dropped"};
    for (const auto& run : mReplays) {
      for (std::size_t v = 0; v < run.vertices.size(); ++v) {
        const auto& r = run.vertices[v];
        const auto* collision = r.collision < 0 ? nullptr : &mCollisions[static_cast<std::size_t>(r.collision)];
        vertices << std::format("{},{},{},{},{},{},{:.1f},{:.1f},{},{},{},{},{},{},{:.4f},{},{:.1f},{},{},{},{},{:.3f},{:.3f},{:.3f}\n", run.name, v, r.listed,
                                mFrameTimes[static_cast<std::size_t>(r.frame)], r.round, r.windowStartUM, r.centreUM, r.zUM, r.window, r.sidebands, r.hidden,
                                r.contributors, r.selected, r.dropped, r.reuse, r.collision, collision == nullptr ? 0. : collision->z * 1.e4,
                                collision == nullptr ? 0U : collision->findable, r.reference, r.startPS, r.widthPS, r.denseShare, r.softShare,
                                r.noParticleShare);
      }
      for (std::size_t c = 0; c < mCollisions.size(); ++c) {
        const auto& collision = mCollisions[c];
        const auto& s = run.collisions[c];
        const auto best = run.score.best[c];
        const auto* vertex = best < 0 ? nullptr : &run.vertices[static_cast<std::size_t>(run.listedToId[static_cast<std::size_t>(best)])];
        const auto* dense = s.dense < 0 ? nullptr : &mCollisions[static_cast<std::size_t>(s.dense)];
        collisions << std::format("{},{},{:.1f},{},{},{},{:d},{:.1f},{},{},{},{},{}", run.name, c, collision.z * 1.e4,
                                  collision.frame < 0 ? 0U : mFrameTimes[static_cast<std::size_t>(collision.frame)], collision.inner, collision.findable,
                                  vertex != nullptr, vertex == nullptr ? 0. : vertex->zUM, vertex == nullptr ? 0U : vertex->contributors,
                                  vertex == nullptr ? 0U : vertex->round, s.accepted, s.kept, s.core);
        for (const auto count : s.states) {
          collisions << ',' << count;
        }
        collisions << std::format(",{},{},{:.1f},{},{},{}\n", s.dense, dense == nullptr ? 0U : dense->findable,
                                  dense == nullptr ? 0. : std::abs(collision.z - dense->z) * 1.e4, s.denseContributors, s.vetoedWindow, s.vetoedSidebands);
      }
      for (const auto& drop : run.drops) {
        const auto& candidate = mCandidates[drop.candidate];
        const auto& owner = run.vertices[static_cast<std::size_t>(drop.owner)];
        if (owner.collision == candidate.collision) {
          continue;
        }
        const auto& source = mCandidates[drop.source];
        drops << std::format("{},{},{},{},{},{},{:.1f},{},{},{},{},{},{}\n", run.name, candidate.collision,
                             mCollisions[static_cast<std::size_t>(candidate.collision)].findable, candidate.fit.zUM, drop.owner, owner.collision, owner.zUM,
                             drop.source, source.collision, kStateNames[static_cast<std::size_t>(run.base[drop.source])], source.fit.zUM, drop.row,
                             rowCollision(drop.row));
      }
    }
  }

  /// The vertexer's vertices into \a path, one CSV line each with its
  /// majority collision (-1 a fake) and that collision's true z and frame,
  /// for the cuts on the list (tools/vertexMakeup.py).
  void writeStageVertices(const std::string& path) const
  {
    std::ofstream out(path);
    if (!out) {
      throw std::runtime_error(std::format("cannot write the vertex list {}", path));
    }
    const auto score = stageScore();
    const auto made = makeup(score.vertices, score.claimedBy);
    out << "frame,timePS,widthPS,zUM,contributors,collision,trueZUM,collisionFrame,findable,inner,strongest,trueShare,pairShare,noParticle,dense,own,"
           "other,dcaUM,chi2RZ,radiusM,shared,reuseClaim,reuseWindow\n";
    for (std::size_t v = 0; v < score.vertices.size(); ++v) {
      const auto& vertex = score.vertices[v];
      const auto& seed = mStageVertices[v];
      const auto* collision = vertex.collision < 0 ? nullptr : &mCollisions[static_cast<std::size_t>(vertex.collision)];
      const auto& m = made[v];
      out << std::format("{},{},{},{:.1f},{},{},{:.1f},{},{},{},{:d},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.2f},{:.3f},{:.3f},{:.1f},{:.4f},{:.4f}\n",
                         vertex.frame, seed.timePS, seed.widthPS, vertex.zUM, vertex.contributors, vertex.collision,
                         collision == nullptr ? 0. : collision->z * 1.e4, collision == nullptr ? -1 : collision->frame,
                         collision == nullptr ? 0U : collision->findable, collision == nullptr ? 0U : collision->inner, m.strongest, m.shares[0], m.shares[1],
                         m.shares[2], m.shares[3], m.shares[4], m.shares[5], m.dcaUM, m.chi2RZ, m.radiusM, m.shared, m.reuseClaim, m.reuseWindow);
    }
  }

 private:
  using Triple = std::array<std::uint32_t, 3>;

  struct Collision {
    double z = 0.;                ///< the median of its measured primaries' production z [cm]
    std::int64_t frame = -1;      ///< the frame most of its layer-1 primary measurements are in
    std::uint32_t inner = 0U;     ///< primaries above the pass's pT floor on layers 0-2, the vertexing pass's
    std::uint32_t primaries = 0U; ///< primaries with any measurement
    /// Physical primaries findable as the tracking QA has them (`isFindable`:
    /// a chain of four layers at its floor).
    std::uint32_t findable = 0U;
  };

  /// A triplet as the vertexer's map measures it (`measureCandidate`)
  /// without its cuts: kAccepted or kDegenerate, the r-z chi2, z at the beam
  /// line [um], the circle's signed distance to the axis [cm], the weight in
  /// a vertex's mean and the time box.
  struct StageFit {
    VertexCandidateStatus status = VertexCandidateStatus::kDegenerate;
    Scalar chi2RZ = 0.F;
    std::int64_t zUM = 0;
    Scalar dcaCM = 0.F;
    std::uint64_t weight = 0U;
    TimeEstPS box;
  };

  struct Candidate {
    Triple triple{};
    std::int64_t frame = 0;
    std::int64_t zUM = 0;
    double z0 = 0.;
    double sigmaZ0 = 0.;
    std::int64_t collision = -1;     ///< index into mCollisions, -1 for background
    bool trueTriplet = false;        ///< one particle on all three measurements
    bool straddles = false;          ///< not all three measurements in one frame
    std::uint8_t sector = 0U;        ///< the middle measurement's azimuth in sixteenths
    std::uint32_t sharedFirst = 1U;  ///< triplets with this one's first tracklet, itself included
    std::uint32_t sharedSecond = 1U; ///< the same for its second tracklet
    std::uint64_t particle = 0U;     ///< the particle on all three, 0 for none
    double chi2RZ = 0.;              ///< the straight line's chi2 in r-z on the points' z errors (1 dof)
    double dcaXY = 0.;               ///< the circle through the three points: its distance of closest approach to the beam axis [cm]
    double radius = 0.;              ///< and its radius [cm], infinite for a straight line: the transverse momentum's order
    StageFit fit{};                  ///< the vertexer's measure of it, for the replay
  };

  void buildCollisions()
  {
    // Per particle: the layers its measurements are on, and its layer-1 frames.
    std::unordered_map<std::uint64_t, std::uint64_t> layers;
    std::unordered_map<std::uint64_t, std::map<std::uint32_t, std::uint32_t>> frames; // by box start
    const auto columns = mFrame->measurements();
    for (std::size_t row = 0; row < mRows.size(); ++row) {
      const auto& truth = mRows[row];
      for (std::size_t i = 0; i < truth.count; ++i) {
        layers[truth.particles[i]] |= std::uint64_t{1} << truth.layer;
        if (truth.layer == 1U) {
          ++frames[truth.particles[i]][columns.timePS[row]];
        }
      }
    }
    std::unordered_map<std::uint64_t, std::map<std::uint32_t, std::uint32_t>> collisionFrames;
    // The collision's z: the median of its measured primaries' production z.
    // The generator's primaries include the daughters of the decays it made,
    // some of them displaced, so no single primary will do.
    std::vector<std::vector<double>> productionZ;
    // Every generator primary's collision, measured or not (`VertexQA.h`'s
    // index takes those with a measured primary).
    for (const auto& [particle, info] : *mParticles) {
      if (!info.primary) {
        continue;
      }
      const auto key = collisionOf(particle);
      auto [at, inserted] = mCollisionIndex.try_emplace(key, mCollisions.size());
      if (inserted) {
        mCollisions.emplace_back();
        productionZ.emplace_back();
      }
      auto& collision = mCollisions[at->second];
      const auto found = layers.find(particle);
      if (found == layers.end()) {
        continue;
      }
      ++collision.primaries;
      productionZ[at->second].push_back(static_cast<double>(info.vz));
      const auto mask = found->second;
      if (static_cast<double>(info.pt) >= mPassMinPtGeV && (mask & 7U) == 7U) {
        ++collision.inner;
      }
      if (const auto f = frames.find(particle); f != frames.end()) {
        for (const auto& [time, count] : f->second) {
          collisionFrames[key][time] += count;
        }
      }
    }
    for (std::size_t c = 0; c < mCollisions.size(); ++c) {
      auto& zs = productionZ[c];
      if (!zs.empty()) {
        mCollisions[c].z = study::median(zs);
      }
    }
    for (const auto& [key, histogram] : collisionFrames) {
      const auto best = std::ranges::max_element(histogram, {}, [](const auto& entry) { return entry.second; });
      mCollisions[mCollisionIndex[key]].frame = frameIndex(best->first);
    }
    for (const auto& [particle, chain] : study::particleLayers(mRows, mComparison)) {
      const auto info = mParticles->find(particle);
      if (info == mParticles->end() || info->second.origin != study::ParticleOrigin::kPhysicalPrimary ||
          !study::isFindable(info->second.pt, chain, mComparison)) {
        continue;
      }
      if (const auto at = mCollisionIndex.find(collisionOf(particle)); at != mCollisionIndex.end()) {
        ++mCollisions[at->second].findable;
      }
    }
  }

  std::int64_t frameIndex(std::uint32_t timePS)
  {
    const auto [at, inserted] = mFrameIndex.try_emplace(timePS, static_cast<std::int64_t>(mFrameIndex.size()));
    return at->second;
  }

  void buildCandidates()
  {
    mCandidates.clear();
    const auto columns = mFrame->measurements();
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> firstShared;
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> secondShared;
    for (const auto& [triple, fit] : mStageFits) {
      ++firstShared[{triple[0], triple[1]}];
      ++secondShared[{triple[1], triple[2]}];
    }
    for (const auto& [triple, fit] : mStageFits) {
      // The weighted straight line z = z0 + slope r.
      double sw = 0.;
      double swr = 0.;
      double swr2 = 0.;
      double swz = 0.;
      double swrz = 0.;
      for (const auto row : triple) {
        const double w = 1. / std::max(static_cast<double>(columns.varZ[row]), 1.e-12);
        const double r = columns.r[row];
        const double z = columns.z[row];
        sw += w;
        swr += w * r;
        swr2 += w * r * r;
        swz += w * z;
        swrz += w * r * z;
      }
      const double det = (sw * swr2) - (swr * swr);
      const double z0 = ((swr2 * swz) - (swr * swrz)) / det;
      const double slope = ((sw * swrz) - (swr * swz)) / det;
      double chi2 = 0.;
      for (const auto row : triple) {
        const double residual = static_cast<double>(columns.z[row]) - z0 - (slope * static_cast<double>(columns.r[row]));
        chi2 += residual * residual / std::max(static_cast<double>(columns.varZ[row]), 1.e-12);
      }
      // The circle through the three points in x-y, and its distance to the axis.
      const double ax = columns.x[triple[0]], ay = columns.y[triple[0]];
      const double bx = columns.x[triple[1]], by = columns.y[triple[1]];
      const double cx = columns.x[triple[2]], cy = columns.y[triple[2]];
      const double d = 2. * ((ax * (by - cy)) + (bx * (cy - ay)) + (cx * (ay - by)));
      double dca = 0.;
      double radius = std::numeric_limits<double>::infinity();
      if (std::abs(d) > 1.e-12) {
        const double a2 = (ax * ax) + (ay * ay), b2 = (bx * bx) + (by * by), c2 = (cx * cx) + (cy * cy);
        const double ux = ((a2 * (by - cy)) + (b2 * (cy - ay)) + (c2 * (ay - by))) / d;
        const double uy = ((a2 * (cx - bx)) + (b2 * (ax - cx)) + (c2 * (bx - ax))) / d;
        dca = std::abs(std::hypot(ux, uy) - std::hypot(ax - ux, ay - uy));
        radius = std::hypot(ax - ux, ay - uy);
      } else {
        dca = std::abs((ax * cy) - (ay * cx)) / std::hypot(cx - ax, cy - ay);
      }
      const auto frame = frameIndex(columns.timePS[triple[1]]);
      Candidate candidate{.frame = frame,
                          .zUM = std::llround(z0 * 1.e4),
                          .z0 = z0,
                          .sigmaZ0 = std::sqrt(swr2 / det),
                          .straddles = frameIndex(columns.timePS[triple[0]]) != frame || frameIndex(columns.timePS[triple[2]]) != frame};
      candidate.chi2RZ = chi2;
      candidate.dcaXY = dca;
      candidate.radius = radius;
      candidate.sharedFirst = firstShared.at({triple[0], triple[1]});
      candidate.sharedSecond = secondShared.at({triple[1], triple[2]});
      candidate.sector = static_cast<std::uint8_t>(std::min(15., std::floor(static_cast<double>(columns.phi[triple[1]]) * 16. / (2. * std::numbers::pi))));
      if (const auto particle = commonParticle(mRows, triple)) {
        candidate.trueTriplet = true;
        candidate.particle = *particle;
        const auto info = mParticles->find(*particle);
        if (info != mParticles->end() && info->second.primary) {
          candidate.collision = static_cast<std::int64_t>(mCollisionIndex.at(collisionOf(*particle)));
        }
      }
      candidate.triple = triple;
      candidate.fit = fit;
      mCandidates.push_back(candidate);
    }
    std::ranges::sort(mCandidates, {}, [](const Candidate& c) { return std::pair{c.frame, c.zUM}; });
    // Every frame is indexed by now: their times by index, and the
    // candidates by their triple, for the replay.
    mFrameTimes.assign(mFrameIndex.size(), 0U);
    for (const auto& [time, frame] : mFrameIndex) {
      mFrameTimes[static_cast<std::size_t>(frame)] = time;
    }
    mByTriple.resize(mCandidates.size());
    for (std::size_t i = 0; i < mByTriple.size(); ++i) {
      mByTriple[i] = i;
    }
    std::ranges::sort(mByTriple, {}, [&](std::size_t i) { return mCandidates[i].triple; });
  }

  /// The candidate on \a triple (`mByTriple`).
  [[nodiscard]] std::size_t candidateOf(const Triple& triple) const
  {
    return *std::ranges::lower_bound(mByTriple, triple, {}, [&](std::size_t i) { return mCandidates[i].triple; });
  }

  /// The collision of \a row's first label: -2 for no particle, -1 for one
  /// whose collision is not in the index.
  [[nodiscard]] std::int64_t rowCollision(std::uint32_t row) const
  {
    const auto& truth = mRows[row];
    if (truth.count == 0U) {
      return -2;
    }
    const auto at = mCollisionIndex.find(collisionOf(truth.particles[0]));
    return at == mCollisionIndex.end() ? -1 : static_cast<std::int64_t>(at->second);
  }

  /// The element at quantile \a q of \a sorted, the one at min(n - 1, q n):
  /// every quantile the study prints. \a sorted must not be empty.
  [[nodiscard]] static double quantileAt(const std::vector<double>& sorted, double q)
  {
    return sorted[std::min(sorted.size() - 1U, static_cast<std::size_t>(q * static_cast<double>(sorted.size())))];
  }

  /// 68/95/99 % and the largest of \a values times \a scale, "n/a" for none.
  static std::string quantiles(std::vector<double> values, double scale)
  {
    if (values.empty()) {
      return "n/a";
    }
    std::ranges::sort(values);
    return std::format("{:.0f} {:.0f} {:.0f} {:.0f}", quantileAt(values, 0.68) * scale, quantileAt(values, 0.95) * scale, quantileAt(values, 0.99) * scale,
                       values.back() * scale);
  }

  void describeCandidates(std::vector<std::string>& lines) const
  {
    std::size_t primary = 0;
    std::size_t secondary = 0;
    std::size_t straddling = 0;
    std::size_t primaryStraddling = 0;
    std::vector<double> residuals;
    std::vector<double> pulls;
    for (const auto& c : mCandidates) {
      straddling += c.straddles ? 1U : 0U;
      if (c.collision >= 0) {
        ++primary;
        primaryStraddling += c.straddles ? 1U : 0U;
        const double dz = c.z0 - mCollisions[static_cast<std::size_t>(c.collision)].z;
        residuals.push_back(std::abs(dz));
        pulls.push_back(std::abs(dz) / c.sigmaZ0);
      } else if (c.trueTriplet) {
        ++secondary;
      }
    }
    const auto n = static_cast<double>(std::max<std::size_t>(mCandidates.size(), 1U));
    lines.push_back(std::format("candidates {}: primary {} ({:.1f} %), secondary {} ({:.1f} %), fake {} ({:.1f} %); straddling frames {} ({} primary)",
                                mCandidates.size(), primary, 100. * static_cast<double>(primary) / n, secondary, 100. * static_cast<double>(secondary) / n,
                                mCandidates.size() - primary - secondary, 100. * static_cast<double>(mCandidates.size() - primary - secondary) / n, straddling,
                                primaryStraddling));
    lines.push_back(std::format("primary candidates' |z0 - z_true| 68/95/99/max [um]: {}; |pull| on the points' errors: {}", quantiles(residuals, 1.e4),
                                quantiles(pulls, 1.)));
    {
      std::vector<double> chi2Primary, chi2Fake, dcaPrimary, dcaFake;
      for (const auto& c : mCandidates) {
        if (c.collision >= 0) {
          chi2Primary.push_back(c.chi2RZ);
          dcaPrimary.push_back(c.dcaXY);
        } else if (!c.trueTriplet) {
          chi2Fake.push_back(c.chi2RZ);
          dcaFake.push_back(c.dcaXY);
        }
      }
      const auto spread = [](std::vector<double> values, double scale) {
        if (values.empty()) {
          return std::string{"n/a"};
        }
        std::ranges::sort(values);
        const auto at = [&](double q) { return quantileAt(values, q) * scale; };
        return std::format("10% {:.3g} 50% {:.3g} 90% {:.3g} 99% {:.3g}", at(0.1), at(0.5), at(0.9), at(0.99));
      };
      lines.push_back(std::format("r-z chi2: primary {} | fake {}", spread(chi2Primary, 1.), spread(chi2Fake, 1.)));
      lines.push_back(std::format("dca xy [um]: primary {} | fake {}", spread(dcaPrimary, 1.e4), spread(dcaFake, 1.e4)));
    }
    std::size_t frames = 0;
    std::size_t busiest = 0;
    for (std::size_t begin = 0; begin < mCandidates.size();) {
      auto end = begin;
      while (end < mCandidates.size() && mCandidates[end].frame == mCandidates[begin].frame) {
        ++end;
      }
      ++frames;
      busiest = std::max(busiest, end - begin);
      begin = end;
    }
    std::uint32_t withInner = 0;
    std::uint32_t findable = 0;
    std::uint32_t findableWithoutInner = 0;
    std::map<std::int64_t, std::uint32_t> perFrame;
    for (const auto& collision : mCollisions) {
      withInner += collision.inner > 0U ? 1U : 0U;
      findable += collision.findable;
      findableWithoutInner += collision.inner == 0U ? collision.findable : 0U;
      if (collision.primaries > 0U) {
        ++perFrame[collision.frame];
      }
    }
    std::array<std::uint32_t, 5> sharing{};
    for (const auto& [frame, count] : perFrame) {
      ++sharing[std::min<std::size_t>(count, sharing.size() - 1U)];
    }
    lines.push_back(std::format(
      "frames with candidates {}, the busiest {} candidates; collisions {} ({} with an inner primary above {:.2f} GeV), "
      "findable primaries {} ({} in collisions without one); frames by collisions 1/2/3/4+: {} {} {} {}",
      frames, busiest, mCollisions.size(), withInner, mPassMinPtGeV, findable, findableWithoutInner, sharing[1], sharing[2], sharing[3], sharing[4]));
  }

  /// Each triplet's two tests in units of its own error, at
  /// its own momentum -- pT from the circle's radius (|B2C Bz| R), 1/(beta p)^2
  /// for the pion hypothesis of the pass (the triplet's own beta, not the
  /// floor's), none for a straight line.
  /// The transverse one: the circle's signed distance to the axis over the
  /// variance its own Jacobian gives -- by finite differences, the points'
  /// r-phi errors, and each barrel scatterer inside the third point (the
  /// layers and the passives the line crosses) as an azimuthal kick that
  /// moves the points beyond it by kick (r_i - r_s) across the track -- and
  /// the beam's transverse spread.
  /// The r-z one: the dip angle's kink at the middle point over the points'
  /// z errors, each chord with its own derivative, and the kink's scattering
  /// (`tripletScattering`: the middle point's layer and the passives
  /// between) rescaled from the floor's 1/(beta p)^2 to the triplet's.
  struct TripletPulls {
    double dca = 0.;
    double rz = 0.;
    double ptGeV = 0.; ///< from the radius, infinite for a straight line
  };

  /// The signed distance of the circle through three transverse points to
  /// the axis [cm], as the vertexer's map has it; a line's for collinear ones.
  static double circleDca(double ax, double ay, double bx, double by, double cx, double cy)
  {
    const double ux0 = ax - bx;
    const double uy0 = ay - by;
    const double vx0 = cx - bx;
    const double vy0 = cy - by;
    const double twiceArea = 2. * ((ux0 * vy0) - (uy0 * vx0));
    if (std::abs(twiceArea) > 1.e-18) {
      const double a2 = (ux0 * ux0) + (uy0 * uy0);
      const double c2 = (vx0 * vx0) + (vy0 * vy0);
      const double ux = ((vy0 * a2) - (uy0 * c2)) / twiceArea;
      const double uy = ((ux0 * c2) - (vx0 * a2)) / twiceArea;
      const double radius = std::hypot(ux, uy);
      const double centre = std::hypot(bx + ux, by + uy);
      const double power = (bx * bx) + (by * by) + (2. * ((bx * ux) + (by * uy)));
      return power / (centre + radius);
    }
    const double dx = cx - ax;
    const double dy = cy - ay;
    return ((ax * cy) - (ay * cx)) / std::hypot(dx, dy);
  }

  [[nodiscard]] TripletPulls tripletPulls(const Candidate& c) const
  {
    const auto columns = mFrame->measurements();
    const auto a = pairPointAt(columns, c.triple[0]);
    const auto b = pairPointAt(columns, c.triple[1]);
    const auto d = pairPointAt(columns, c.triple[2]);
    const auto& table = mPass->scattering;
    const auto& cut = mPass->pairCut;
    TripletPulls out;
    const bool curved = std::isfinite(c.radius) && c.radius > 0.;
    out.ptGeV = curved ? static_cast<double>(cut.momentumCurvature) * c.radius : std::numeric_limits<double>::infinity();
    const double slope = static_cast<double>(d.z - a.z) / static_cast<double>(d.r - a.r);
    const double g = std::sqrt(1. + (slope * slope));
    const double mass2 = static_cast<double>(cut.mass2);
    // 1/(beta p)^2 at the triplet's pT, and the floor's, as the pair cut has it.
    const double ownInvPt2 = curved ? (1. + (mass2 / (out.ptGeV * out.ptGeV))) / (out.ptGeV * out.ptGeV) : 0.;
    const double scale = cut.invBetaMomentum2 > 0.F ? ownInvPt2 / static_cast<double>(cut.invBetaMomentum2) : 0.;
    // Transverse.
    const std::array<std::array<double, 2>, 3> points{{{a.x, a.y}, {b.x, b.y}, {d.x, d.y}}};
    const std::array<double, 3> radii{a.r, b.r, d.r};
    const auto dcaOf = [&](const std::array<std::array<double, 2>, 3>& p) { return circleDca(p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1]); };
    const double dca0 = dcaOf(points);
    constexpr double kStep = 1.e-5; // [cm]
    double varPoints = 0.;
    const std::array<double, 3> varRPhi{a.varRPhi, b.varRPhi, d.varRPhi};
    for (std::size_t i = 0; i < 3U; ++i) {
      auto moved = points;
      moved[i][0] += -kStep * points[i][1] / radii[i];
      moved[i][1] += kStep * points[i][0] / radii[i];
      const double derivative = (dcaOf(moved) - dca0) / kStep;
      varPoints += derivative * derivative * varRPhi[i];
    }
    // Across the track: perpendicular to the chord from the first point to the third.
    const double chordX = points[2][0] - points[0][0];
    const double chordY = points[2][1] - points[0][1];
    const double chord = std::hypot(chordX, chordY);
    const double nx = -chordY / chord;
    const double ny = chordX / chord;
    double varScattering = 0.;
    const auto kick = [&](double rs, double k) {
      if (!(rs < radii[2]) || !(k > 0.)) {
        return;
      }
      auto moved = points;
      for (std::size_t i = 0; i < 3U; ++i) {
        if (radii[i] > rs) {
          moved[i][0] += kStep * (radii[i] - rs) * nx;
          moved[i][1] += kStep * (radii[i] - rs) * ny;
        }
      }
      const double derivative = (dcaOf(moved) - dca0) / kStep;
      // theta0^2 = k g / (beta p)^2 at p = pT g; the azimuthal kick theta0 g.
      const double theta2 = k * g * (1. + (mass2 / (out.ptGeV * out.ptGeV * g * g))) / (out.ptGeV * out.ptGeV * g * g);
      varScattering += derivative * derivative * theta2 * g * g;
    };
    if (curved) {
      const double z0 = c.z0;
      const auto crosses = [&](const ScatteringLayer& layer) {
        const double zAt = z0 + (slope * static_cast<double>(layer.position));
        return zAt >= static_cast<double>(layer.low) && zAt <= static_cast<double>(layer.high);
      };
      for (std::size_t layer = 0; layer < ScatteringTable<Det>::kLayerCount; ++layer) {
        if (Detector<Det>::isBarrelLayer(static_cast<LayerID>(layer)) && crosses(table.layers[layer])) {
          kick(static_cast<double>(table.layers[layer].position), static_cast<double>(table.layers[layer].k));
        }
      }
      if constexpr (ScatteringTable<Det>::kPassiveCount > 0) {
        for (std::size_t i = 0; i < ScatteringTable<Det>::kPassiveCount; ++i) {
          const auto& passive = table.passives[i];
          if (passive.barrel && crosses(passive.at)) {
            kick(static_cast<double>(passive.at.position), static_cast<double>(passive.at.k));
          }
        }
      }
    }
    out.dca = static_cast<double>(c.fit.dcaCM) / std::sqrt(varPoints + varScattering + static_cast<double>(mBeamVarTransverse));
    // r-z: the kink of the dip angle at the middle point, each chord with its own derivative.
    const double d1 = static_cast<double>(b.r - a.r);
    const double d2 = static_cast<double>(d.r - b.r);
    const double t1 = static_cast<double>(b.z - a.z) / d1;
    const double t2 = static_cast<double>(d.z - b.z) / d2;
    const double a1 = 1. / (1. + (t1 * t1));
    const double a2 = 1. / (1. + (t2 * t2));
    const double varKink = (a1 * a1 * static_cast<double>(a.varZ) / (d1 * d1)) +
                           ((a1 / d1) + (a2 / d2)) * ((a1 / d1) + (a2 / d2)) * static_cast<double>(b.varZ) +
                           (a2 * a2 * static_cast<double>(d.varZ) / (d2 * d2));
    const double varKinkScattering = static_cast<double>(tripletScattering<Det>(table, cut, mRows[c.triple[1]].layer, a, b, d)) * scale;
    out.rz = (std::atan(t2) - std::atan(t1)) / std::sqrt(varKink + varKinkScattering);
    return out;
  }

  /// The pulls of the true primary triplets by their true pT and |eta|, the
  /// true secondaries' and the combinatorial ones' (|pull| at 50/68/95/99 %,
  /// the transverse one's mean), the radius's pT against the truth, and how
  /// many of each class cuts in sigmas keep against the absolute cuts.
  void describePulls(std::vector<std::string>& lines) const
  {
    if (!mPass) {
      return; // no slice of the iteration
    }
    struct Class {
      std::string name;
      std::vector<double> dca{};
      std::vector<double> rz{};
      std::vector<double> ptRatio{};
      double dcaSum = 0.;
      std::array<std::uint64_t, 5> both{}, onlyDca{}, onlyRz{};
      std::uint64_t absolute = 0U;
      std::uint64_t absoluteChi2 = 0U;
      std::uint64_t absoluteDca = 0U;
    };
    constexpr std::array<double, 5> kSigmas{2., 3., 4., 5., 6.};
    constexpr std::array<double, 6> kPtEdges{0., 0.3, 0.5, 1., 2., 1.e30};
    std::vector<Class> classes;
    for (std::size_t bin = 0; bin + 1U < kPtEdges.size(); ++bin) {
      classes.push_back({.name = bin + 2U < kPtEdges.size() ? std::format("true primary, pT {}-{} GeV", kPtEdges[bin], kPtEdges[bin + 1U])
                                                            : std::format("true primary, pT > {} GeV", kPtEdges[bin])});
    }
    classes.push_back({.name = "true primary, |eta| < 1"});
    classes.push_back({.name = "true primary, |eta| > 1"});
    classes.push_back({.name = "true primary, all"});
    classes.push_back({.name = "true secondary"});
    classes.push_back({.name = "combinatorial"});
    const auto& settings = mStage.settings();
    for (const auto& c : mCandidates) {
      if (c.fit.status != VertexCandidateStatus::kAccepted) {
        continue;
      }
      const auto pulls = tripletPulls(c);
      std::vector<std::size_t> into;
      double ptTrue = 0.;
      if (c.trueTriplet && c.collision >= 0) {
        const auto& info = mParticles->at(c.particle);
        ptTrue = static_cast<double>(info.pt);
        std::size_t bin = 0;
        while (bin + 2U < kPtEdges.size() && ptTrue >= kPtEdges[bin + 1U]) {
          ++bin;
        }
        into = {bin, std::abs(info.eta) < 1.F ? kPtEdges.size() - 1U : kPtEdges.size(), kPtEdges.size() + 1U};
      } else {
        into = {c.trueTriplet ? kPtEdges.size() + 2U : kPtEdges.size() + 3U};
      }
      const bool chi2In = !(settings.maxChi2RZ > 0.F && c.fit.chi2RZ > settings.maxChi2RZ);
      const bool dcaAbsolute = !(settings.maxDcaCM > 0.F && std::abs(c.fit.dcaCM) > settings.maxDcaCM);
      const bool absolute = chi2In && dcaAbsolute;
      for (const auto k : into) {
        auto& out = classes[k];
        out.dca.push_back(std::abs(pulls.dca));
        out.rz.push_back(std::abs(pulls.rz));
        out.dcaSum += pulls.dca;
        if (ptTrue > 0.) {
          out.ptRatio.push_back(pulls.ptGeV / ptTrue);
        }
        out.absolute += absolute ? 1U : 0U;
        out.absoluteChi2 += chi2In ? 1U : 0U;
        out.absoluteDca += dcaAbsolute ? 1U : 0U;
        for (std::size_t s = 0; s < kSigmas.size(); ++s) {
          const bool dcaIn = std::abs(pulls.dca) < kSigmas[s];
          const bool rzIn = std::abs(pulls.rz) < kSigmas[s];
          out.both[s] += dcaIn && rzIn ? 1U : 0U;
          out.onlyDca[s] += dcaIn ? 1U : 0U;
          out.onlyRz[s] += rzIn ? 1U : 0U;
        }
      }
    }
    for (auto& out : classes) {
      std::ranges::sort(out.dca);
      std::ranges::sort(out.rz);
      std::ranges::sort(out.ptRatio);
    }
    const auto at = [](const std::vector<double>& sorted, double q) { return sorted.empty() ? std::numeric_limits<double>::quiet_NaN() : quantileAt(sorted, q); };
    lines.push_back(
      std::format("triplet pulls: the triplets past the pass's cuts; |pull| 50/68/95/99 %, transverse and r-z; the transverse "
                  "pull's mean; the radius's pT over the true one 16/50/84 %; kept by both pulls < 2/3/4/5/6 sigma; by the absolute cuts "
                  "(chi2 {}, dca {} um)",
                  settings.maxChi2RZ, settings.maxDcaCM * 1.e4F));
    for (auto& out : classes) {
      const auto n = out.dca.size();
      if (n == 0U) {
        continue;
      }
      const auto share = [&](std::uint64_t kept) { return std::format("{:.3f}", static_cast<double>(kept) / static_cast<double>(n)); };
      std::string kept;
      for (std::size_t s = 0; s < kSigmas.size(); ++s) {
        kept += std::format(" {}", share(out.both[s]));
      }
      std::string dcaOnly;
      std::string rzOnly;
      for (std::size_t s = 0; s < kSigmas.size(); ++s) {
        dcaOnly += std::format(" {}", share(out.onlyDca[s]));
        rzOnly += std::format(" {}", share(out.onlyRz[s]));
      }
      const double mean = out.dcaSum / static_cast<double>(n);
      lines.push_back(
        std::format("  {}: {} | dca {:.2f} {:.2f} {:.2f} {:.2f} mean {:+.2f} | rz {:.2f} {:.2f} {:.2f} {:.2f} | pT ratio {:.2f} {:.2f} {:.2f} | "
                    "kept both{} | dca only{} | rz only{} | absolute {} (chi2 {}, dca {})",
                    out.name, n, at(out.dca, 0.5), at(out.dca, 0.68), at(out.dca, 0.95), at(out.dca, 0.99), mean, at(out.rz, 0.5), at(out.rz, 0.68),
                    at(out.rz, 0.95), at(out.rz, 0.99), at(out.ptRatio, 0.16), at(out.ptRatio, 0.5), at(out.ptRatio, 0.84), kept, dcaOnly, rzOnly,
                    share(out.absolute), share(out.absoluteChi2), share(out.absoluteDca)));
    }
  }

  struct Vertex {
    std::int64_t frame = 0;
    double zUM = 0.;
    std::uint32_t contributors = 0U;
    std::int64_t collision = -1; ///< the majority; -1 background
  };

  struct Score {
    std::vector<Vertex> vertices;
    std::vector<std::int64_t> claimedBy; ///< per candidate, its vertex, -1 free, -2 cut or not kept
    std::vector<std::int64_t> best;      ///< per collision, its vertex with the most contributors, or -1
    std::uint32_t fakes = 0U;
    std::uint32_t splits = 0U;
    std::uint32_t merged = 0U;
  };

  /// The vertexer's vertices and claims scored: each vertex's frame and
  /// majority from the candidates it claimed.
  [[nodiscard]] Score stageScore() const
  {
    std::vector<std::int64_t> claimedBy;
    auto vertices = stageVertices(claimedBy);
    return scoreClaims(std::move(vertices), std::move(claimedBy));
  }

  /// The vertexer's vertices as `stageScore` scores them, and which vertex
  /// claimed each candidate.
  std::vector<Vertex> stageVertices(std::vector<std::int64_t>& claimedBy) const
  {
    claimedBy.assign(mCandidates.size(), -2);
    std::vector<Vertex> vertices(mStageVertices.size());
    for (std::size_t i = 0; i < mCandidates.size(); ++i) {
      const auto found = mStageClaims.find(mCandidates[i].triple);
      claimedBy[i] = found == mStageClaims.end() ? -2 : found->second;
      if (claimedBy[i] >= 0) {
        auto& vertex = vertices[static_cast<std::size_t>(claimedBy[i])];
        vertex.frame = mCandidates[i].frame;
        ++vertex.contributors;
      }
    }
    for (std::size_t v = 0; v < vertices.size(); ++v) {
      vertices[v].zUM = static_cast<double>(mStageVertices[v].z) * 1.e4;
    }
    return vertices;
  }

  [[nodiscard]] Score scoreClaims(std::vector<Vertex> vertices, std::vector<std::int64_t> claimedBy) const
  {
    Score result;
    result.vertices = std::move(vertices);
    result.claimedBy = std::move(claimedBy);
    std::vector<std::map<std::int64_t, std::uint32_t>> labels(result.vertices.size());
    std::vector<std::map<std::int64_t, std::uint32_t>> wentTo(mCollisions.size());
    std::vector<std::uint32_t> contributors;
    for (const auto& vertex : result.vertices) {
      contributors.push_back(vertex.contributors);
    }
    for (std::size_t i = 0; i < mCandidates.size(); ++i) {
      const auto collision = mCandidates[i].collision;
      const auto vertex = result.claimedBy[i];
      if (vertex >= 0) {
        ++labels[static_cast<std::size_t>(vertex)][collision];
      }
      if (collision >= 0) {
        ++wentTo[static_cast<std::size_t>(collision)][vertex];
      }
    }
    auto assigned = scoreVertexClaims(labels, contributors, wentTo);
    for (std::size_t v = 0; v < result.vertices.size(); ++v) {
      result.vertices[v].collision = assigned.collision[v];
    }
    result.best = std::move(assigned.best);
    result.fakes = assigned.fakes;
    result.splits = assigned.splits;
    result.merged = assigned.merged;
    return result;
  }

  /// \a result as one row at \a settings: the settings; the vertices, the
  /// collisions found, the efficiency for those with an inner primary, the
  /// findable primaries covered; fakes, splits, merged; the most vertices in
  /// one frame, in three consecutive frames with candidates, frames over 8;
  /// the found vertices' |dz| 68/95/99 % and largest [um]; the median
  /// contributors.
  [[nodiscard]] std::string describe(const Score& result, const Settings& settings) const
  {
    std::uint32_t found = 0U;
    std::uint32_t withInner = 0U;
    std::uint32_t foundWithInner = 0U;
    std::uint64_t covered = 0U;
    std::uint64_t findable = 0U;
    std::vector<double> residuals;
    std::vector<double> contributors;
    for (std::size_t c = 0; c < mCollisions.size(); ++c) {
      const auto& collision = mCollisions[c];
      findable += collision.findable;
      withInner += collision.inner > 0U ? 1U : 0U;
      if (result.best[c] < 0) {
        continue;
      }
      const auto& vertex = result.vertices[static_cast<std::size_t>(result.best[c])];
      ++found;
      foundWithInner += collision.inner > 0U ? 1U : 0U;
      covered += collision.findable;
      residuals.push_back(std::abs(vertex.zUM - (collision.z * 1.e4)));
      contributors.push_back(vertex.contributors);
    }
    std::ranges::sort(contributors);
    // Vertices per frame, the frames in time order.
    std::map<std::uint32_t, std::uint32_t> perFrame;
    for (const auto& [time, frame] : mFrameIndex) {
      perFrame[time] = 0U;
    }
    for (const auto& vertex : result.vertices) {
      ++perFrame[mFrameTimes[static_cast<std::size_t>(vertex.frame)]];
    }
    std::vector<std::uint32_t> counts;
    for (const auto& [time, count] : perFrame) {
      if (count > 0U) {
        counts.push_back(count);
      }
    }
    std::uint32_t most = 0U;
    std::uint32_t mostThree = 0U;
    std::uint32_t over = 0U;
    for (std::size_t f = 0; f < counts.size(); ++f) {
      most = std::max(most, counts[f]);
      over += counts[f] > 8U ? 1U : 0U;
      mostThree = std::max(mostThree, counts[f] + (f > 0U ? counts[f - 1U] : 0U) + (f + 1U < counts.size() ? counts[f + 1U] : 0U));
    }
    return std::format(
      "{:4.0f} {:4.0f} {:5d} {:5d} {:4d} {:2.0f} | {:8d} {:6d}/{:<5d} {:6.1f} % {:8d}/{:<8d} | {:5d} {:6d} {:6d} | {:4d} {:4d} {:4d} | {:>24} | {:.0f}",
      settings.maxChi2RZ, settings.maxDcaUM, settings.windowUM, settings.claimUM, settings.minContributors, settings.significance, result.vertices.size(),
      found, mCollisions.size(),
      100. * static_cast<double>(foundWithInner) / std::max(1., static_cast<double>(withInner)), covered, findable, result.fakes, result.splits, result.merged,
      most, mostThree, over, quantiles(residuals, 1.), contributors.empty() ? 0. : contributors[contributors.size() / 2U]);
  }

  /// The efficiency table, its lines, and the diagnostics of the fakes and
  /// the missed collisions.
  struct Efficiency {
    std::string table;
    std::vector<std::string> diagnostics;
  };

  /// Each frame's strongest of \a vertices: the first of the most
  /// contributors.
  [[nodiscard]] static std::map<std::int64_t, std::size_t> strongestByFrame(const std::vector<Vertex>& vertices)
  {
    std::map<std::int64_t, std::size_t> strongest;
    for (std::size_t v = 0; v < vertices.size(); ++v) {
      const auto [at, inserted] = strongest.try_emplace(vertices[v].frame, v);
      if (!inserted && vertices[at->second].contributors < vertices[v].contributors) {
        at->second = v;
      }
    }
    return strongest;
  }

  /// \a result's efficiency by the collisions' inner primaries and by their
  /// findable ones at \a settings, \a who naming whose list it is, and the
  /// diagnostics of its fakes and missed collisions.
  [[nodiscard]] Efficiency efficiencyByMultiplicity(const Settings& settings, const Score& result, std::string_view who) const
  {
    Efficiency out;
    constexpr std::array<std::uint32_t, 8> kEdges{0U, 1U, 3U, 5U, 10U, 20U, 50U, 100U};
    std::array<std::uint32_t, kEdges.size()> all{};
    std::array<std::uint32_t, kEdges.size()> found{};
    std::array<std::uint64_t, kEdges.size()> findable{};
    std::array<std::uint64_t, kEdges.size()> covered{};
    for (std::size_t c = 0; c < mCollisions.size(); ++c) {
      const auto& collision = mCollisions[c];
      std::size_t bin = 0;
      while (bin + 1U < kEdges.size() && collision.inner >= kEdges[bin + 1U]) {
        ++bin;
      }
      ++all[bin];
      findable[bin] += collision.findable;
      if (result.best[c] >= 0) {
        ++found[bin];
        covered[bin] += collision.findable;
      }
    }
    {
      // Where the fakes are: their contributors and their distance to their
      // frame's strongest vertex; the same for the real vertices that are not
      // their frame's strongest (pile-up).
      const auto strongest = strongestByFrame(result.vertices);
      std::vector<double> fakeContributors;
      std::vector<double> fakeDistance;
      std::vector<double> fakeRatio;
      std::vector<double> pileContributors;
      std::vector<double> pileDistance;
      std::vector<double> pileRatio;
      for (std::size_t v = 0; v < result.vertices.size(); ++v) {
        const auto& vertex = result.vertices[v];
        const auto& top = result.vertices[strongest.at(vertex.frame)];
        if (strongest.at(vertex.frame) == v) {
          continue;
        }
        const double distance = std::abs(vertex.zUM - top.zUM);
        const double ratio = static_cast<double>(vertex.contributors) / static_cast<double>(top.contributors);
        if (vertex.collision < 0) {
          fakeContributors.push_back(vertex.contributors);
          fakeDistance.push_back(distance);
          fakeRatio.push_back(ratio);
        } else if (result.best[static_cast<std::size_t>(vertex.collision)] == static_cast<std::int64_t>(v)) {
          pileContributors.push_back(vertex.contributors);
          pileDistance.push_back(distance);
          pileRatio.push_back(ratio);
        }
      }
      const auto lowQuantiles = [](std::vector<double> values, double scale) {
        if (values.empty()) {
          return std::string{"n/a"};
        }
        std::ranges::sort(values);
        const auto at = [&](double q) { return quantileAt(values, q) * scale; };
        return std::format("min {:.3g} 5% {:.3g} 50% {:.3g} 95% {:.3g} max {:.3g}", values.front() * scale, at(0.05), at(0.5), at(0.95), values.back() * scale);
      };
      {
        std::vector<std::uint16_t> sectors(result.vertices.size(), 0U);
        std::vector<std::uint32_t> trueCount(result.vertices.size(), 0U);
        std::vector<std::map<std::uint64_t, std::uint32_t>> byParticle(result.vertices.size());
        for (std::size_t i = 0; i < mCandidates.size(); ++i) {
          const auto v = result.claimedBy[i];
          if (v < 0) {
            continue;
          }
          const auto& c = mCandidates[i];
          sectors[static_cast<std::size_t>(v)] |= static_cast<std::uint16_t>(1U << c.sector);
          if (c.trueTriplet) {
            ++trueCount[static_cast<std::size_t>(v)];
            ++byParticle[static_cast<std::size_t>(v)][c.particle];
          }
        }
        std::vector<double> fakeSectors, realSectors, fakeTrue, fakeTop, realTop;
        for (std::size_t v = 0; v < result.vertices.size(); ++v) {
          const auto& vertex = result.vertices[v];
          const double n = vertex.contributors;
          std::uint32_t top = 0U;
          for (const auto& [particle, count] : byParticle[v]) {
            top = std::max(top, count);
          }
          if (vertex.collision < 0) {
            fakeSectors.push_back(std::popcount(sectors[v]));
            fakeTrue.push_back(trueCount[v] / n);
            fakeTop.push_back(top / n);
          } else if (result.best[static_cast<std::size_t>(vertex.collision)] == static_cast<std::int64_t>(v)) {
            realSectors.push_back(std::popcount(sectors[v]));
            realTop.push_back(top / n);
          }
        }
        out.diagnostics.push_back(std::format("  fakes: azimuth sectors of 16 hit {}; true-triplet share {}; one particle's share {}", lowQuantiles(fakeSectors, 1.),
                                              lowQuantiles(fakeTrue, 1.), lowQuantiles(fakeTop, 1.)));
        out.diagnostics.push_back(
          std::format("  real: azimuth sectors of 16 hit {}; one particle's share {}", lowQuantiles(realSectors, 1.), lowQuantiles(realTop, 1.)));
      }
      out.diagnostics.push_back(std::format("  fakes not their frame's strongest ({}): contributors {}; |dz| to the strongest [mm] {}; contributors / strongest's {}",
                                            fakeContributors.size(), lowQuantiles(fakeContributors, 1.), lowQuantiles(fakeDistance, 1.e-3),
                                            lowQuantiles(fakeRatio, 1.)));
      out.diagnostics.push_back(std::format("  real pile-up vertices ({}): contributors {}; |dz| to the strongest [mm] {}; contributors / strongest's {}",
                                            pileContributors.size(), lowQuantiles(pileContributors, 1.), lowQuantiles(pileDistance, 1.e-3),
                                            lowQuantiles(pileRatio, 1.)));
    }
    // The missed collisions with 20 inner primaries or more, and why.
    for (std::size_t c = 0; c < mCollisions.size(); ++c) {
      const auto& collision = mCollisions[c];
      if (result.best[c] >= 0 || collision.inner < 20U) {
        continue;
      }
      std::map<std::int64_t, std::uint32_t> claims;
      std::uint32_t own = 0U;
      for (std::size_t i = 0; i < mCandidates.size(); ++i) {
        if (mCandidates[i].collision == static_cast<std::int64_t>(c)) {
          ++own;
          ++claims[result.claimedBy[i]];
        }
      }
      std::string where;
      for (const auto& [claim, count] : claims) {
        where += claim >= 0 ? std::format(" vertex of collision {} at {:.4f} cm: {}", result.vertices[static_cast<std::size_t>(claim)].collision,
                                          result.vertices[static_cast<std::size_t>(claim)].zUM * 1.e-4, count)
                            : std::format(" {}: {}", claim == -1 ? "free" : "cut", count);
      }
      out.diagnostics.push_back(std::format("  missed collision {}: z {:.4f} cm, inner {}, findable {}, candidates {};{}", c, collision.z, collision.inner,
                                            collision.findable, own, where));
    }
    std::string& line = out.table;
    line = std::format(
      "{}efficiency by inner primaries (chi2 {}, dca {} um, w {} um, claim {} um, minC {}, k {}): bin found/all (findable primaries covered/all):", who,
      settings.maxChi2RZ, settings.maxDcaUM, settings.windowUM, settings.claimUM, settings.minContributors, settings.significance);
    for (std::size_t bin = 0; bin < kEdges.size(); ++bin) {
      line += std::format(" [{}{}] {}/{} ({}/{})", kEdges[bin], bin + 1U < kEdges.size() ? std::format("-{}", kEdges[bin + 1U] - 1U) : std::string{"+"},
                          found[bin], all[bin], covered[bin], findable[bin]);
    }
    // The same by the collision's findable primaries, which do not depend on
    // the pass's pT floor: passes at different floors compare bin by bin.
    constexpr std::array<std::uint32_t, 10> kFindableEdges{0U, 10U, 20U, 50U, 100U, 200U, 500U, 1000U, 2500U, 5000U};
    std::array<std::uint32_t, kFindableEdges.size()> allByFindable{};
    std::array<std::uint32_t, kFindableEdges.size()> foundByFindable{};
    for (std::size_t c = 0; c < mCollisions.size(); ++c) {
      std::size_t bin = 0;
      while (bin + 1U < kFindableEdges.size() && mCollisions[c].findable >= kFindableEdges[bin + 1U]) {
        ++bin;
      }
      ++allByFindable[bin];
      foundByFindable[bin] += result.best[c] >= 0 ? 1U : 0U;
    }
    line += std::format("\n{}efficiency by findable primaries: bin found/all:", who);
    for (std::size_t bin = 0; bin < kFindableEdges.size(); ++bin) {
      line += std::format(" [{}{}] {}/{}", kFindableEdges[bin],
                          bin + 1U < kFindableEdges.size() ? std::format("-{}", kFindableEdges[bin + 1U] - 1U) : std::string{"+"}, foundByFindable[bin],
                          allByFindable[bin]);
    }
    return out;
  }

  /// What a vertex is made of, over its claimed candidates:
  /// whether it is its frame's strongest, the share of true triplets, of
  /// triplets with two measurements of one particle, and of the measurements
  /// by origin -- no particle (QED, noise), the frame's strongest vertex's
  /// collision (dense; 0 for the strongest itself), the vertex's own, another
  /// -- and the medians of the candidates' transverse dca [um], r-z chi2,
  /// circle radius [m] and most shared tracklet.
  struct Makeup {
    bool strongest = false;
    /// Of its candidates, the share with a measurement that a candidate past
    /// the table's cuts within the frame's strongest vertex's claim (window)
    /// of its z also has: what that vertex would own if it owned its
    /// candidates' measurements; 0 for the strongest itself.
    double reuseClaim = 0.;
    double reuseWindow = 0.;
    std::array<double, 6> shares{}; ///< true, pair, no particle, dense, own, other
    double dcaUM = 0.;
    double chi2RZ = 0.;
    double radiusM = 0.;
    double shared = 0.;
  };

  [[nodiscard]] std::vector<Makeup> makeup(const std::vector<Vertex>& vertices, const std::vector<std::int64_t>& claims) const
  {
    const auto strongest = strongestByFrame(vertices);
    std::vector<std::vector<std::size_t>> claimed(vertices.size());
    for (std::size_t i = 0; i < mCandidates.size(); ++i) {
      if (claims[i] >= 0) {
        claimed[static_cast<std::size_t>(claims[i])].push_back(i);
      }
    }
    // Per frame, the measurements of the candidates past the table's cuts
    // within the strongest vertex's claim and window.
    const auto& cuts = mTable;
    std::map<std::int64_t, std::pair<std::unordered_set<std::uint32_t>, std::unordered_set<std::uint32_t>>> owned;
    for (const auto& [frame, top] : strongest) {
      owned[frame];
    }
    for (const auto& c : mCandidates) {
      const auto frame = strongest.find(c.frame);
      if (frame == strongest.end() || (cuts.maxChi2RZ > 0. && c.chi2RZ > cuts.maxChi2RZ) || (cuts.maxDcaUM > 0. && c.dcaXY * 1.e4 > cuts.maxDcaUM)) {
        continue;
      }
      const double distance = std::abs(static_cast<double>(c.zUM) - vertices[frame->second].zUM);
      auto& [inClaim, inWindow] = owned[c.frame];
      if (distance <= static_cast<double>(cuts.claimUM)) {
        inClaim.insert(c.triple.begin(), c.triple.end());
      }
      if (distance <= static_cast<double>(cuts.windowUM)) {
        inWindow.insert(c.triple.begin(), c.triple.end());
      }
    }
    std::vector<Makeup> out(vertices.size());
    for (std::size_t v = 0; v < vertices.size(); ++v) {
      const auto& vertex = vertices[v];
      const auto top = strongest.at(vertex.frame);
      const auto dense = top == v ? std::int64_t{-3} : vertices[top].collision;
      auto& m = out[v];
      m.strongest = top == v;
      std::vector<double> dca, chi2, radius, shared;
      for (const auto i : claimed[v]) {
        const auto& c = mCandidates[i];
        m.shares[0] += c.trueTriplet ? 1. : 0.;
        bool pair = false;
        for (const auto& [first, second] : {std::pair{0U, 1U}, std::pair{1U, 2U}, std::pair{0U, 2U}}) {
          const auto& row = mRows[c.triple[first]];
          for (std::size_t k = 0; k < row.count && !pair; ++k) {
            pair = mRows[c.triple[second]].carries(row.particles[k]);
          }
        }
        m.shares[1] += pair ? 1. : 0.;
        if (!m.strongest) {
          const auto& [inClaim, inWindow] = owned.at(vertex.frame);
          m.reuseClaim += std::ranges::any_of(c.triple, [&](std::uint32_t row) { return inClaim.contains(row); }) ? 1. : 0.;
          m.reuseWindow += std::ranges::any_of(c.triple, [&](std::uint32_t row) { return inWindow.contains(row); }) ? 1. : 0.;
        }
        for (const auto row : c.triple) {
          const auto origin = rowCollision(row);
          m.shares[origin == -2 ? 2U : (origin >= 0 && origin == dense)            ? 3U
                                     : (origin >= 0 && origin == vertex.collision) ? 4U
                                                                                   : 5U] += 1. / 3.;
        }
        dca.push_back(c.dcaXY * 1.e4);
        chi2.push_back(c.chi2RZ);
        radius.push_back(std::isfinite(c.radius) ? c.radius * 1.e-2 : 1.e3);
        shared.push_back(std::max(c.sharedFirst, c.sharedSecond));
      }
      const auto n = std::max<double>(1., static_cast<double>(claimed[v].size()));
      for (auto& share : m.shares) {
        share /= n;
      }
      m.reuseClaim /= n;
      m.reuseWindow /= n;
      m.dcaUM = study::median(dca);
      m.chi2RZ = study::median(chi2);
      m.radiusM = study::median(radius);
      m.shared = study::median(shared);
    }
    return out;
  }

  // ---- The replay ----
  //
  // The vertexer's rounds (`Vertexer::find`, `findInFrame`) on the
  // timeframe's distinct triplets as it measures them (`StageFit`): its cuts,
  // the deduplication per slice, per frame the rounds -- the peak test, the
  // centre, the claim, the weighted selection, the ownership -- and the list
  // in the vertexer's order. Every variant owns, as the vertexer does, every
  // measured triplet within the claim; at the vertexer's cuts
  // the replay is the vertexer (checked in `replayAll`).

  /// Where a candidate ended, per collision in `ReplayCollision::states`.
  enum class ReplayState : std::uint8_t { kChi2,
                                          kDca,
                                          kDegenerate,
                                          kDeduplicated,
                                          kFree,
                                          kOwn,
                                          kOther,
                                          kFake,
                                          kDropped,
                                          kCount };
  static constexpr std::size_t kReplayStates = static_cast<std::size_t>(ReplayState::kCount);

  /// The candidates' cuts: the vertexer's absolute ones, the r-z pull in
  /// place of the chi2 (the dca absolute), both pulls, or both pulls with
  /// the absolute dca as a cap (`tripletPulls`).
  enum class CutMode : std::uint8_t { kAbsolute,
                                      kRzPull,
                                      kPulls,
                                      kPullsCapped };

  struct ReplayVariant {
    std::string name;
    CutMode cuts = CutMode::kAbsolute;
    double sigma = 0.;  ///< the pulls' cut
    bool prior = false; ///< starting from `setPriorVertices`'
  };

  struct ReplayVertex {
    std::int64_t frame = 0;
    std::uint32_t round = 0U;
    std::int64_t windowStartUM = 0; ///< the winning window's first candidate
    double centreUM = 0.;           ///< the window's free candidates' mean, the claim's centre
    double zUM = 0.;                ///< the fitted z, the selection's weighted mean
    std::uint32_t window = 0U;      ///< W: the window's free candidates
    std::uint32_t sidebands = 0U;   ///< S: the two sidebands' free candidates
    std::uint32_t hidden = 0U;      ///< kept candidates in the sidebands not free: claimed or dropped before
    std::uint32_t contributors = 0U;
    std::uint32_t selected = 0U; ///< of them, in the fit
    std::uint32_t dropped = 0U;  ///< free candidates its ownership dropped
    /// Of its claimed candidates, the share with a measurement of a triplet
    /// past the cuts within an earlier vertex's claim in its frame.
    double reuse = 0.;
    std::int64_t listed = -1;    ///< its index in the list
    std::int64_t collision = -1; ///< its claimed candidates' majority, -1 background
    std::int64_t reference = -1; ///< the list index of its frame's strongest listed vertex, the first of the most contributors
    std::uint32_t startPS = 0U;
    std::uint32_t widthPS = 0U;
    VertexSeed seed{}; ///< as the vertexer lists it
    /// Of its claimed candidates' measurements, the shares whose first label
    /// is the frame's strongest other vertex's collision, a particle below the
    /// pass's pT floor, and no particle.
    double denseShare = 0.;
    double softShare = 0.;
    double noParticleShare = 0.;
  };

  /// A true triplet of a collision dropped by a vertex's ownership: the row
  /// it shared with an owned triplet, that triplet (the source), the owner.
  struct ReplayDrop {
    std::size_t candidate = 0U;
    std::int64_t owner = 0; ///< the replay vertex
    std::size_t source = 0U;
    std::uint32_t row = 0U;
  };

  /// A window that held the contributor floor when its frame's rounds ended:
  /// the sideband test refused it.
  struct ReplayVeto {
    std::int64_t frame = 0;
    std::int64_t zUM = 0; ///< its start
    std::uint32_t window = 0U;
    std::uint32_t sidebands = 0U;
  };

  struct ReplayCollision {
    std::uint32_t accepted = 0U; ///< its true triplets past the cuts
    std::uint32_t kept = 0U;     ///< of them, kept by the deduplication
    std::uint32_t core = 0U;     ///< of those, within the window of its true z
    std::array<std::uint32_t, kReplayStates> states{};
    std::int64_t dense = -1; ///< the collision of its frame's strongest listed vertex, -1 a fake or none
    std::uint32_t denseContributors = 0U;
    /// The refused window over its true z in its frame with the most
    /// candidates, its count and its sidebands'; 0 for none.
    std::uint32_t vetoedWindow = 0U;
    std::uint32_t vetoedSidebands = 0U;
  };

  struct ReplayRun {
    std::string name;
    std::vector<ReplayVertex> vertices{}; ///< every vertex found, in the order found
    std::vector<std::int64_t> listedToId{};
    Score score{};                   ///< its claims per candidate: the listed vertex that claimed it, -1 free, -2 not kept
    std::vector<ReplayState> base{}; ///< per candidate: its state after the cuts and the deduplication
    std::vector<ReplayCollision> collisions{};
    std::uint64_t dropped = 0U;
    std::vector<ReplayDrop> drops{};
    std::vector<ReplayVeto> vetoes{};
    std::vector<std::uint32_t> owned{}; ///< every row a vertex owned, in the order owned
  };

  /// The candidates' states after the cuts and the deduplication: kChi2,
  /// kDca, kDegenerate, kDeduplicated, or kFree for the kept. The
  /// deduplication as the vertexer's: per slice, over the slice's accepted
  /// triplets in increasing chi2, ties in the slice's triplet order; a
  /// candidate's state the one of the slice that owns its frame.
  [[nodiscard]] std::vector<ReplayState> replayCandidates(const ReplayVariant& variant) const
  {
    const auto& s = mStage.settings();
    std::vector<ReplayState> state(mCandidates.size(), ReplayState::kFree);
    std::vector<std::uint8_t> passes(mCandidates.size(), 0U);
    for (std::size_t i = 0; i < mCandidates.size(); ++i) {
      const auto& c = mCandidates[i];
      const bool chi2Out = variant.cuts == CutMode::kAbsolute ? s.maxChi2RZ > 0.F && c.fit.chi2RZ > s.maxChi2RZ : !(std::abs(mPulls[i].rz) < variant.sigma);
      const bool dcaAbsoluteOut = s.maxDcaCM > 0.F && std::abs(c.fit.dcaCM) > s.maxDcaCM;
      const bool dcaPullOut = !(std::abs(mPulls[i].dca) < variant.sigma);
      const bool dcaOut = variant.cuts == CutMode::kPulls ? dcaPullOut : variant.cuts == CutMode::kPullsCapped ? dcaPullOut || dcaAbsoluteOut
                                                                                                               : dcaAbsoluteOut;
      if (c.fit.status != VertexCandidateStatus::kAccepted) {
        state[i] = ReplayState::kDegenerate;
      } else if (chi2Out) {
        state[i] = ReplayState::kChi2;
      } else if (dcaOut) {
        state[i] = ReplayState::kDca;
      } else {
        passes[i] = 1U;
      }
    }
    const auto columns = mFrame->measurements();
    std::vector<std::uint8_t> taken(mFrame->size(), 0U);
    std::vector<std::size_t> accepted;
    std::vector<std::uint8_t> inSlice(mCandidates.size(), 0U);
    std::vector<std::uint8_t> inOwner(mCandidates.size(), 0U);
    for (const auto& [owned, triples] : mSliceTriplets) {
      accepted.clear();
      for (const auto& triple : triples) {
        if (owned.owns(columns.timePS[triple[1]])) {
          inOwner[candidateOf(triple)] = 1U; // the ownership's sources are the owning slice's
        }
        // A triple twice in the slice: its second triplet is the
        // deduplication's, the candidate's state the first's.
        if (const auto i = candidateOf(triple); passes[i] != 0U && inSlice[i] == 0U) {
          inSlice[i] = 1U;
          accepted.push_back(i);
        }
      }
      std::ranges::stable_sort(accepted, {}, [&](std::size_t i) { return mCandidates[i].fit.chi2RZ; });
      for (const auto i : accepted) {
        const auto& rows = mCandidates[i].triple;
        const bool kept = std::ranges::none_of(rows, [&](std::uint32_t row) { return taken[row] != 0U; });
        if (kept) {
          for (const auto row : rows) {
            taken[row] = 1U;
          }
        }
        if (owned.owns(columns.timePS[rows[1]])) {
          state[i] = kept ? ReplayState::kFree : ReplayState::kDeduplicated;
          passes[i] = 2U;
        }
      }
      for (const auto i : accepted) {
        inSlice[i] = 0U;
        for (const auto row : mCandidates[i].triple) {
          taken[row] = 0U;
        }
      }
    }
    if (std::ranges::count(passes, std::uint8_t{1U}) != 0) {
      throw std::logic_error("the replay: a candidate past the cuts in no slice's frames");
    }
    if (std::ranges::count(inOwner, std::uint8_t{0U}) != 0) {
      throw std::logic_error("the replay: a triplet not in the slice that owns its frame");
    }
    return state;
  }

  [[nodiscard]] ReplayRun replay(const ReplayVariant& variant, const std::vector<ReplayState>& base) const
  {
    const auto& s = mStage.settings();
    const std::int64_t window = s.windowUM;
    const std::int64_t claim = s.claimUM;
    const Scalar k2 = 2.F * s.significance * s.significance;
    ReplayRun run{.name = variant.name};
    const auto n = mCandidates.size();
    std::vector<std::int64_t> claimedBy(n, -1); // the replay vertex
    std::vector<std::uint8_t> dropped(n, 0U);
    std::vector<std::int64_t> rowOwner(mFrame->size(), -1);
    std::vector<std::size_t> rowSource(mFrame->size(), 0U);
    std::vector<std::uint8_t> rowSeen(mFrame->size(), 0U);  // for the reuse share
    std::vector<std::uint8_t> rowPrior(mFrame->size(), 0U); // owned by an earlier pass's vertex
    std::vector<std::uint32_t> touched;
    const auto zOf = [&](std::size_t i) { return mCandidates[i].fit.zUM; };
    for (std::size_t begin = 0; begin < n;) {
      auto end = begin;
      while (end < n && mCandidates[end].frame == mCandidates[begin].frame) {
        ++end;
      }
      std::vector<std::size_t> kept;
      std::vector<std::size_t> sources;  // the ownership's triplets, by z
      std::vector<std::size_t> accepted; // the reuse share's, by z
      for (auto i = begin; i < end; ++i) {
        const bool pastCuts = base[i] == ReplayState::kFree || base[i] == ReplayState::kDeduplicated;
        if (base[i] == ReplayState::kFree) {
          kept.push_back(i);
        }
        if (pastCuts) {
          accepted.push_back(i);
        }
        if (base[i] != ReplayState::kDegenerate) {
          sources.push_back(i);
        }
      }
      std::ranges::stable_sort(kept, {}, zOf);
      std::ranges::stable_sort(sources, {}, zOf);
      std::ranges::stable_sort(accepted, {}, zOf);
      const auto m = kept.size();
      std::vector<std::int64_t> z(m);
      for (std::size_t k = 0; k < m; ++k) {
        z[k] = zOf(kept[k]);
      }
      const auto at = [&](std::int64_t value) { return static_cast<std::size_t>(std::ranges::lower_bound(z, value) - z.begin()); };
      const auto isFree = [&](std::size_t k) { return claimedBy[kept[k]] < 0 && dropped[kept[k]] == 0U; };
      // The triplets of \a list within \a half of \a centre.
      // As the vertexer's ownership (`withinZ`).
      const auto within = [&](const std::vector<std::size_t>& list, Scalar centre, std::int64_t half, const auto& body) {
        const auto reach = static_cast<std::int64_t>(centre);
        for (auto it = std::ranges::lower_bound(list, reach - half - 1, {}, zOf); it != list.end() && zOf(*it) <= reach + half + 1; ++it) {
          if (withinZ(zOf(*it), centre, static_cast<Scalar>(half))) {
            body(*it);
          }
        }
      };
      // The earlier pass's vertices first -- their claims
      // taken, their sources' rows owned, the free candidates on those rows
      // dropped -- as rounds before the first.
      if (variant.prior) {
        if (const auto found = mPriorVertices.find(mFrameTimes[static_cast<std::size_t>(mCandidates[begin].frame)]); found != mPriorVertices.end()) {
          for (const auto prior : found->second) {
            for (auto k = at(static_cast<std::int64_t>(prior) - claim - 1); k < m && static_cast<Scalar>(z[k]) <= prior + static_cast<Scalar>(claim); ++k) {
              if (isFree(k) && withinZ(z[k], prior, static_cast<Scalar>(claim))) {
                dropped[kept[k]] = 1U;
              }
            }
            within(sources, prior, claim, [&](std::size_t i) {
              for (const auto row : mCandidates[i].triple) {
                if (rowPrior[row] == 0U) {
                  rowPrior[row] = 1U;
                  touched.push_back(row);
                }
              }
            });
          }
          for (std::size_t k = 0; k < m; ++k) {
            if (isFree(k) && std::ranges::any_of(mCandidates[kept[k]].triple, [&](std::uint32_t row) { return rowPrior[row] != 0U; })) {
              dropped[kept[k]] = 1U;
              ++run.dropped;
            }
          }
        }
      }
      std::vector<std::uint32_t> count(m + 1U, 0U);
      std::uint32_t rounds = 0U;
      while (true) {
        ++rounds;
        for (std::size_t k = 0; k < m; ++k) {
          count[k + 1U] = count[k] + (isFree(k) ? 1U : 0U);
        }
        const auto countIn = [&](std::int64_t low, std::int64_t high) { return count[at(high)] - count[at(low)]; };
        std::int64_t bestExcess = 0;
        std::size_t best = m;
        for (std::size_t k = 0; k < m; ++k) {
          if (!isFree(k)) {
            continue;
          }
          const auto excess = peakExcess(z[k], window, claim, s.minContributors, k2, countIn);
          if (excess > bestExcess) {
            bestExcess = excess;
            best = k;
          }
        }
        if (best == m) {
          for (std::size_t k = 0; k < m; ++k) {
            if (isFree(k) && countIn(z[k], z[k] + window) >= s.minContributors) {
              run.vetoes.push_back({.frame = mCandidates[begin].frame,
                                    .zUM = z[k],
                                    .window = countIn(z[k], z[k] + window),
                                    .sidebands = countIn(z[k] - claim - window, z[k] - claim) + countIn(z[k] + window + claim, z[k] + (2 * window) + claim)});
            }
          }
          break;
        }
        const auto id = static_cast<std::int64_t>(run.vertices.size());
        ReplayVertex vertex{.frame = mCandidates[begin].frame, .round = rounds, .windowStartUM = z[best]};
        const auto windowEnd = at(z[best] + window);
        std::int64_t sum = 0;
        for (auto k = best; k < windowEnd; ++k) {
          sum += isFree(k) ? z[k] : 0;
        }
        const auto members = static_cast<std::int64_t>(count[windowEnd] - count[best]);
        const Scalar centre = static_cast<Scalar>(sum) / static_cast<Scalar>(members);
        vertex.centreUM = centre;
        vertex.window = countIn(z[best], z[best] + window);
        const std::array<std::pair<std::int64_t, std::int64_t>, 2> sides{
          {{z[best] - claim - window, z[best] - claim}, {z[best] + window + claim, z[best] + (2 * window) + claim}}};
        for (const auto& [low, high] : sides) {
          vertex.sidebands += countIn(low, high);
          vertex.hidden += static_cast<std::uint32_t>(at(high) - at(low)) - countIn(low, high);
        }
        // The claim and the weighted selection, as `findInFrame`.
        std::vector<std::int64_t> claimedZ;
        std::vector<std::uint64_t> claimedWeight;
        std::vector<std::size_t> claimedIndex;
        std::uint32_t earliest = std::numeric_limits<std::uint32_t>::max();
        std::uint32_t latest = 0U;
        for (auto k = at(static_cast<std::int64_t>(centre) - claim - 1); k < m && static_cast<Scalar>(z[k]) <= centre + static_cast<Scalar>(claim); ++k) {
          const Scalar distance = static_cast<Scalar>(z[k]) > centre ? static_cast<Scalar>(z[k]) - centre : centre - static_cast<Scalar>(z[k]);
          if (!isFree(k) || distance > static_cast<Scalar>(claim)) {
            continue;
          }
          const auto& fit = mCandidates[kept[k]].fit;
          claimedBy[kept[k]] = id;
          ++vertex.contributors;
          earliest = std::min(earliest, fit.box.lower());
          latest = std::max(latest, fit.box.upper());
          claimedZ.push_back(z[k]);
          claimedWeight.push_back(fit.weight);
          claimedIndex.push_back(kept[k]);
        }
        std::vector<std::uint8_t> selected(claimedZ.size(), 0U);
        std::uint64_t weightSum = 0U;
        std::int64_t weightedZSum = 0;
        std::int64_t selectedCount = 0;
        const auto select = [&](std::size_t k) {
          selected[k] = 1U;
          weightSum += claimedWeight[k];
          weightedZSum += static_cast<std::int64_t>(claimedWeight[k]) * claimedZ[k];
          ++selectedCount;
        };
        for (std::size_t k = 0; k < claimedZ.size(); ++k) {
          if (withinZ(claimedZ[k], centre, static_cast<Scalar>(window))) {
            select(k);
          }
        }
        for (bool grew = true; grew;) {
          grew = false;
          const auto moments = zMoments(weightSum, weightedZSum);
          for (std::size_t k = 0; k < claimedZ.size(); ++k) {
            if (selected[k] == 0U && withinError(claimedZ[k], claimedWeight[k], moments, s.clipNSigma)) {
              select(k);
              grew = true;
            }
          }
        }
        const auto settled = zMoments(weightSum, weightedZSum);
        std::uint64_t scatterSum = 0U;
        for (std::size_t k = 0; k < claimedZ.size(); ++k) {
          if (selected[k] != 0U) {
            scatterSum += candidateScatter(claimedZ[k], claimedWeight[k], settled.mean);
          }
        }
        const auto seed = vertexSeedOf(
          {.weightSum = weightSum, .weightedZSum = weightedZSum, .scatterSum = scatterSum, .count = selectedCount, .earliestPS = earliest, .latestPS = latest},
          VertexConstraint{});
        vertex.zUM = static_cast<double>(seed.z) * 1.e4;
        vertex.selected = static_cast<std::uint32_t>(selectedCount);
        vertex.startPS = seed.timePS;
        vertex.widthPS = seed.widthPS;
        vertex.seed = seed;
        // The reuse share against the earlier vertices of the frame, then
        // this one's triplets past the cuts within its claim.
        std::uint32_t reused = 0U;
        for (const auto i : claimedIndex) {
          reused += std::ranges::any_of(mCandidates[i].triple, [&](std::uint32_t row) { return rowSeen[row] != 0U; }) ? 1U : 0U;
        }
        vertex.reuse = vertex.contributors == 0U ? 0. : static_cast<double>(reused) / static_cast<double>(vertex.contributors);
        within(accepted, centre, claim, [&](std::size_t i) {
          for (const auto row : mCandidates[i].triple) {
            rowSeen[row] = 1U;
            touched.push_back(row);
          }
        });
        // The ownership, and the free candidates it drops.
        within(sources, centre, claim, [&](std::size_t i) {
          for (const auto row : mCandidates[i].triple) {
            if (rowOwner[row] < 0) {
              rowOwner[row] = id;
              rowSource[row] = i;
              touched.push_back(row);
              run.owned.push_back(row);
            }
          }
        });
        for (std::size_t k = 0; k < m; ++k) {
          const auto& rows = mCandidates[kept[k]].triple;
          const auto owned = std::ranges::find_if(rows, [&](std::uint32_t row) { return rowOwner[row] >= 0; });
          if (isFree(k) && owned != rows.end()) {
            dropped[kept[k]] = 1U;
            ++vertex.dropped;
            if (mCandidates[kept[k]].collision >= 0) {
              run.drops.push_back({.candidate = kept[k], .owner = id, .source = rowSource[*owned], .row = *owned});
            }
          }
        }
        run.dropped += vertex.dropped;
        run.vertices.push_back(vertex);
      }
      for (const auto row : touched) {
        rowOwner[row] = -1;
        rowSeen[row] = 0U;
        rowPrior[row] = 0U;
      }
      touched.clear();
      begin = end;
    }
    // Each vertex's majority collision.
    std::vector<std::map<std::int64_t, std::uint32_t>> labels(run.vertices.size());
    for (std::size_t i = 0; i < n; ++i) {
      if (claimedBy[i] >= 0) {
        ++labels[static_cast<std::size_t>(claimedBy[i])][mCandidates[i].collision];
      }
    }
    for (std::size_t v = 0; v < run.vertices.size(); ++v) {
      run.vertices[v].collision = std::ranges::max_element(labels[v], {}, [](const auto& entry) { return entry.second; })->first;
    }
    // The list in the vertexer's order -- frames in time, a frame's vertices
    // by z, ties in the order found.
    std::vector<std::size_t> order(run.vertices.size());
    for (std::size_t v = 0; v < order.size(); ++v) {
      order[v] = v;
    }
    std::ranges::stable_sort(order, {},
                             [&](std::size_t v) { return std::pair{mFrameTimes[static_cast<std::size_t>(run.vertices[v].frame)], run.vertices[v].seed.z}; });
    for (std::size_t l = 0; l < order.size(); ++l) {
      run.vertices[order[l]].listed = static_cast<std::int64_t>(l);
    }
    std::vector<Vertex> vertices;
    std::map<std::int64_t, std::int64_t> strongest; // per frame, its first listed vertex of the most contributors
    for (const auto v : order) {
      const auto& vertex = run.vertices[v];
      if (vertex.listed < 0) {
        continue;
      }
      run.listedToId.push_back(static_cast<std::int64_t>(v));
      vertices.push_back({.frame = vertex.frame, .zUM = vertex.zUM, .contributors = vertex.contributors, .collision = vertex.collision});
      const auto [top, inserted] = strongest.try_emplace(vertex.frame, vertex.listed);
      if (!inserted && vertices[static_cast<std::size_t>(top->second)].contributors < vertex.contributors) {
        top->second = vertex.listed;
      }
    }
    for (auto& vertex : run.vertices) {
      const auto top = strongest.find(vertex.frame);
      vertex.reference = top == strongest.end() ? -1 : top->second;
    }
    for (std::size_t i = 0; i < n; ++i) {
      if (claimedBy[i] < 0) {
        continue;
      }
      auto& vertex = run.vertices[static_cast<std::size_t>(claimedBy[i])];
      const auto dense = vertex.reference < 0 || vertex.reference == vertex.listed
                           ? std::int64_t{-3}
                           : run.vertices[static_cast<std::size_t>(run.listedToId[static_cast<std::size_t>(vertex.reference)])].collision;
      for (const auto row : mCandidates[i].triple) {
        const auto& truth = mRows[row];
        if (truth.count == 0U) {
          vertex.noParticleShare += 1.;
          continue;
        }
        const auto particle = truth.particles[0];
        const auto at = mCollisionIndex.find(collisionOf(particle));
        vertex.denseShare += at != mCollisionIndex.end() && static_cast<std::int64_t>(at->second) == dense ? 1. : 0.;
        const auto info = mParticles->find(particle);
        vertex.softShare += info != mParticles->end() && static_cast<double>(info->second.pt) < mPassMinPtGeV ? 1. : 0.;
      }
    }
    for (auto& vertex : run.vertices) {
      const double rows = 3. * std::max(1., static_cast<double>(vertex.contributors));
      vertex.denseShare /= rows;
      vertex.softShare /= rows;
      vertex.noParticleShare /= rows;
    }
    std::vector<std::int64_t> claims(n, -2);
    for (std::size_t i = 0; i < n; ++i) {
      if (base[i] == ReplayState::kFree) {
        claims[i] = claimedBy[i] < 0 ? -1 : run.vertices[static_cast<std::size_t>(claimedBy[i])].listed;
      }
    }
    run.score = scoreClaims(std::move(vertices), std::move(claims));
    // Per collision: its true triplets' states, and its frame's strongest.
    run.collisions.assign(mCollisions.size(), {});
    for (std::size_t i = 0; i < n; ++i) {
      const auto& c = mCandidates[i];
      if (c.collision < 0) {
        continue;
      }
      auto& out = run.collisions[static_cast<std::size_t>(c.collision)];
      auto state = base[i];
      if (state == ReplayState::kFree || state == ReplayState::kDeduplicated) {
        ++out.accepted;
      }
      if (state == ReplayState::kFree) {
        ++out.kept;
        out.core +=
          std::abs(static_cast<double>(c.fit.zUM) - (mCollisions[static_cast<std::size_t>(c.collision)].z * 1.e4)) <= static_cast<double>(window) ? 1U : 0U;
        if (dropped[i] != 0U) {
          state = ReplayState::kDropped;
        } else if (claimedBy[i] >= 0) {
          const auto& vertex = run.vertices[static_cast<std::size_t>(claimedBy[i])];
          state = vertex.collision == c.collision ? ReplayState::kOwn : vertex.collision >= 0 ? ReplayState::kOther
                                                                                              : ReplayState::kFake;
        }
      }
      ++out.states[static_cast<std::size_t>(state)];
    }
    for (std::size_t c = 0; c < mCollisions.size(); ++c) {
      const auto top = strongest.find(mCollisions[c].frame);
      if (top != strongest.end()) {
        const auto& vertex = run.vertices[static_cast<std::size_t>(run.listedToId[static_cast<std::size_t>(top->second)])];
        run.collisions[c].dense = vertex.collision;
        run.collisions[c].denseContributors = vertex.contributors;
      }
    }
    for (const auto& veto : run.vetoes) {
      for (std::size_t c = 0; c < mCollisions.size(); ++c) {
        const auto trueZ = mCollisions[c].z * 1.e4;
        auto& out = run.collisions[c];
        if (mCollisions[c].frame == veto.frame && static_cast<double>(veto.zUM) <= trueZ && trueZ < static_cast<double>(veto.zUM + window) &&
            veto.window > out.vetoedWindow) {
          out.vetoedWindow = veto.window;
          out.vetoedSidebands = veto.sidebands;
        }
      }
    }
    return run;
  }

  /// Whether the slope intervals of \a from and \a to admit their line: a
  /// tracklet's line, for each observed tracklet in the
  /// per-link shape counts (`observeSlice`).
  template <typename Columns>
  [[nodiscard]] static bool pairShaped(const Columns& columns, std::uint32_t from, std::uint32_t to)
  {
    const auto admits = [&](std::uint32_t at) {
      return shapeAdmits(columns, at, columns.x[from], columns.y[from], columns.z[from], columns.x[to], columns.y[to], columns.z[to]);
    };
    return admits(from) && admits(to);
  }

  /// The replay's variants, its check against the vertexer, and a line each.
  void replayAll(std::vector<std::string>& lines)
  {
    mPulls.clear();
    mPulls.reserve(mCandidates.size());
    for (const auto& c : mCandidates) {
      mPulls.push_back(c.fit.status == VertexCandidateStatus::kAccepted ? tripletPulls(c) : TripletPulls{});
    }
    mReplays.clear();
    // The vertexer (first, the gate's), and its cuts in sigmas of the
    // candidate's own errors.
    std::vector<ReplayVariant> variants{
      {.name = "stage"},
      {.name = "rz3", .cuts = CutMode::kRzPull, .sigma = 3.},
      {.name = "pulls3", .cuts = CutMode::kPulls, .sigma = 3.},
      {.name = "pulls3-cap", .cuts = CutMode::kPullsCapped, .sigma = 3.},
    };
    if (!mPriorVertices.empty()) {
      variants.push_back({.name = "prior-stage", .prior = true});
      variants.push_back({.name = "prior-rz3", .cuts = CutMode::kRzPull, .sigma = 3., .prior = true});
      variants.push_back({.name = "prior-pulls3", .cuts = CutMode::kPulls, .sigma = 3., .prior = true});
    }
    for (const auto& variant : variants) {
      auto base = replayCandidates(variant);
      mReplays.push_back(replay(variant, base));
      mReplays.back().base = std::move(base);
    }
    // The gate: the replay at the vertexer's settings is the vertexer --
    // its list in order (box, z and variance bitwise, contributors) and
    // every candidate's claim, free (-1) and not kept (-2) included.
    std::vector<std::int64_t> stageClaims;
    static_cast<void>(stageVertices(stageClaims));
    const auto& first = mReplays.front();
    const auto& firstClaims = first.score.claimedBy;
    std::size_t sameVertices = 0U;
    const auto bits = [](Scalar value) { return std::bit_cast<std::uint32_t>(value); };
    for (std::size_t l = 0; l < std::min(first.listedToId.size(), mStageVertices.size()); ++l) {
      const auto& replayed = first.vertices[static_cast<std::size_t>(first.listedToId[l])];
      const auto& staged = mStageVertices[l];
      sameVertices += replayed.seed.timePS == staged.timePS && replayed.seed.widthPS == staged.widthPS && bits(replayed.seed.z) == bits(staged.z) &&
                          bits(replayed.seed.varZ) == bits(staged.varZ) && replayed.contributors == mStageContributors[l]
                        ? 1U
                        : 0U;
    }
    std::size_t sameClaims = 0U;
    for (std::size_t i = 0; i < mCandidates.size(); ++i) {
      sameClaims += firstClaims[i] == stageClaims[i] ? 1U : 0U;
    }
    const bool identical = sameVertices == first.listedToId.size() && sameVertices == mStageVertices.size() && sameClaims == mCandidates.size();
    const auto check = std::format(
      "replay of the vertexer: {}; vertices identical in order (box, z, variance, contributors) {} of the replay's {} and "
      "the vertexer's {}; candidates with the same claim {} of {}; triples twice in a slice {}",
      identical ? "IDENTICAL" : "DIFFERENT", sameVertices, first.listedToId.size(), mStageVertices.size(), sameClaims, mCandidates.size(), mStageRepeats);
    if (!identical) {
      std::cerr << check << '\n';
      throw std::runtime_error("the replay differs from the vertexer");
    }
    lines.push_back(check);
    constexpr std::array<std::uint32_t, 13> kEdges{0U, 1U, 2U, 3U, 4U, 10U, 20U, 50U, 100U, 200U, 500U, 1000U, 2500U};
    for (const auto& run : mReplays) {
      std::uint32_t found = 0U;
      std::array<std::uint32_t, kEdges.size()> all{};
      std::array<std::uint32_t, kEdges.size()> foundIn{};
      for (std::size_t c = 0; c < mCollisions.size(); ++c) {
        std::size_t bin = 0;
        while (bin + 1U < kEdges.size() && mCollisions[c].findable >= kEdges[bin + 1U]) {
          ++bin;
        }
        ++all[bin];
        const bool hit = run.score.best[c] >= 0;
        found += hit ? 1U : 0U;
        foundIn[bin] += hit ? 1U : 0U;
      }
      // The triplets the pass made (not refused by their shapes), those past
      // the cuts, and of those the share with one collision's particle on
      // all three measurements.
      std::size_t made = 0U;
      std::size_t pastCuts = 0U;
      std::size_t pastCutsTrue = 0U;
      for (std::size_t i = 0; i < mCandidates.size(); ++i) {
        const auto state = run.base[i];
        ++made;
        if (state == ReplayState::kFree || state == ReplayState::kDeduplicated) {
          ++pastCuts;
          pastCutsTrue += mCandidates[i].trueTriplet && mCandidates[i].collision >= 0 ? 1U : 0U;
        }
      }
      std::string line = std::format(
        "replay {}: vertices {} of {} found, collisions {}/{}, fakes {} splits {} merged {}, candidates dropped {}; triplets {}, past the cuts {}, "
        "true {}; by findable (QA):",
        run.name, run.listedToId.size(), run.vertices.size(), found, mCollisions.size(), run.score.fakes, run.score.splits, run.score.merged, run.dropped, made,
        pastCuts, pastCutsTrue);
      for (std::size_t bin = 0; bin < kEdges.size(); ++bin) {
        line += std::format(" [{}{}] {}/{}", kEdges[bin], bin + 1U < kEdges.size() ? std::format("-{}", kEdges[bin + 1U] - 1U) : std::string{"+"}, foundIn[bin],
                            all[bin]);
      }
      lines.push_back(line);
    }
  }

  const TimeFrame<Det>* mFrame;
  std::span<const study::TruthRow> mRows;
  const study::ParticleTable* mParticles;
  std::size_t mIteration;
  double mPassMinPtGeV;
  study::ComparisonSettings mComparison; ///< the tracking QA's findability, the collisions' findable primaries (`isFindable`)
  Vertexer mStage;
  Settings mTable; ///< the vertexer's settings from the configuration, for the tables and the makeup
  std::unique_ptr<BoundedMemoryResource> mStagePool;
  std::unique_ptr<o2::gtrk::detail::TaskArena> mArena;
  std::size_t mSlices = 0;
  std::map<std::uint16_t, std::size_t> mTrackletsByLink;
  std::map<std::uint16_t, std::size_t> mTrackletsShaped; ///< of them, those whose slope intervals admit their line
  std::vector<Collision> mCollisions;
  std::unordered_map<std::uint64_t, std::size_t> mCollisionIndex;
  std::map<std::uint32_t, std::int64_t> mFrameIndex; ///< by box start
  std::vector<Candidate> mCandidates;
  std::vector<std::size_t> mByTriple;            ///< the candidates' indices by their triple
  std::optional<IterationParameters<Det>> mPass; ///< the observed pass's, for the triplets' pulls
  Scalar mBeamVarTransverse = 0.F;               ///< the beam's transverse spread squared, the same
  std::map<Triple, std::int64_t> mStageClaims;
  std::vector<VertexSeed> mStageVertices;
  std::vector<std::uint32_t> mStageContributors;
  std::size_t mStageRepeats = 0U;  ///< triplets whose triple the slice already had
  std::size_t mStagePoolPeak = 0U; ///< the vertexer's work arrays' peak over the slices [bytes]
  std::vector<std::pair<Triple, StageFit>> mStageFits;
  std::vector<std::pair<OwnedFrames, std::vector<Triple>>> mSliceTriplets; ///< per slice, its triplets and the frames it owns
  std::map<std::uint32_t, std::vector<Scalar>> mPriorVertices;             ///< an earlier pass's vertices' z [um] by frame time
  std::vector<std::uint32_t> mFrameTimes;                                  ///< by frame index, its box start
  std::vector<ReplayRun> mReplays;
  std::vector<TripletPulls> mPulls; ///< per candidate, its pulls (`tripletPulls`)
  std::uint32_t mStageRounds = 0U;
  double mStageSeconds = 0.;
  std::size_t mStageCandidates = 0U;
  BeamEstimate mStageBeam;
  double mBeamX = 0.; ///< the rows' beam line, from the dump
  double mBeamY = 0.;
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_VERTEXSTUDY_H_
