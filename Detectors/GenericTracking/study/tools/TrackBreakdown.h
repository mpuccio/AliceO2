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

/// \file TrackBreakdown.h
/// \brief Where a run's tracks sit: the tracks by
///        their number of measurements, and the findable particles of one
///        origin found complete, incomplete or not, with where an incomplete
///        one's missing measurements went.
///
/// Under the score's rules (TrackComparison.h): a track is good when every
/// row carries one particle, a fake otherwise; a particle's good track with
/// the most layers is its best, the others are clones; a particle is found
/// complete when its best track has a row on every layer the particle has
/// one on. Per track length: the tracks, the fakes, the clones, and the
/// best tracks of a complete particle.
///
/// The class is the score's findable particles of one population (the
/// physical primaries, the decay daughters or the V0s'; a chain of
/// `minLayers`, the score's pT floor). An incomplete particle's missing
/// layers are the layers it has a row on that its best track has none on.
/// Each missing layer's rows are looked up in the tracks: in another good
/// track of the same particle (a split), in a fake labelled with the
/// particle (its majority), in another fake, in another particle's good
/// track (the row carries both), or in none -- free after the run. With
/// several rows on the layer the first of that order holds. The layer is
/// also placed against the best track: before its first row, after its last
/// or between them, by the distance from the beam axis at the particle's
/// production z, hypot(r, z - vz) (no MC hits: the order of a first outward
/// pass; for a primary from the beam spot the distance from its vertex but
/// for the transverse offset).
///
/// Per particle of the class (`writeParticles`, a CSV): its truth, its layers,
/// its status, its best track's layers and iteration, and its missing layers
/// by fate and place, a bit per layer -- for splits by layer pattern, pT or
/// the pass that found it (tools/incompleteParticles.py). A particle not
/// found has every layer's fate, no place: whether a later pass could still
/// take its rows. Its inbound layers:
/// those whose rows all lie past the turning point -- each nearer the axis
/// than a row before it along the path, the rows ordered by |z - vz| (a
/// primary's longitudinal distance grows along its helix). Rows within the
/// margin in |z - vz| are not ordered against each other (pz near 0, rows
/// of one plane). A heuristic: a hard scatter can turn a track too.

#ifndef GTRK_STUDY_TRACKBREAKDOWN_H_
#define GTRK_STUDY_TRACKBREAKDOWN_H_

#include "GenericTracking/Study/TrackComparison.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <format>
#include <fstream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::gtrk::study
{

class TrackBreakdown
{
 public:
  /// Where a missing layer's row is after the run.
  enum Fate : std::uint8_t { kSplit,
                             kOwnFake,
                             kFake,
                             kOtherGood,
                             kFree,
                             kFates };
  /// Where a missing layer sits against the best track's rows.
  enum Place : std::uint8_t { kBefore,
                              kBetween,
                              kAfter,
                              kPlaces };
  static constexpr std::size_t kLengths = study::kScoreLengthSlots;
  static constexpr std::size_t kMaxLayers = 64;

  /// A particle's status.
  enum Status : std::uint8_t { kComplete,
                               kIncomplete,
                               kFakeOnly,
                               kNotFound };
  /// One particle of the class: missing = layers & ~best; the missing layers'
  /// fates and places a bit per layer (between: missing & ~before & ~after).
  struct Particle {
    std::uint64_t id = 0U;
    study::ParticleTruth truth;
    std::uint64_t layers = 0U;
    Status status = kNotFound;
    std::uint64_t best = 0U;         ///< the best track's layers
    std::uint8_t bestIteration = 0U; ///< the iteration that made it
    std::uint32_t goodTracks = 0U;
    std::array<std::uint64_t, kFates> fates{};
    std::uint64_t before = 0U;
    std::uint64_t after = 0U;
    std::uint64_t inbound = 0U;
  };

  /// \a tracks of one run scored against \a rows and \a truth under \a comparison,
  /// the class's particles those of \a population; \a iterations, the
  /// iteration of each track, for the particles' records.
  TrackBreakdown(std::span<const study::TruthRow> rows, const study::ParticleTable& truth, const study::TrackRows& tracks,
                 const study::ComparisonSettings& comparison, study::Population population, std::span<const std::uint8_t> iterations)
  {
    const auto length = [](std::size_t size) { return std::min(size, kLengths - 1U); };
    std::vector<std::uint64_t> good(tracks.size(), study::kNoTrackParticle);
    std::vector<std::uint64_t> majority(tracks.size(), study::kNoTrackParticle); // a fake's majority particle
    std::vector<std::uint64_t> layersOf(tracks.size(), 0U);
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> tracksOf;
    std::vector<std::vector<std::uint32_t>> tracksOfRow(rows.size());
    for (std::uint32_t t = 0; t < tracks.size(); ++t) {
      const auto& track = tracks[t];
      for (const auto row : track) {
        layersOf[t] |= std::uint64_t{1} << rows[row].layer;
        tracksOfRow[row].push_back(t);
      }
      good[t] = study::goodParticle(rows, track);
      ++mTracks[length(track.size())];
      if (good[t] == study::kNoTrackParticle) {
        ++mFakes[length(track.size())];
        majority[t] = study::majorityParticle(rows, track).particle;
      } else {
        tracksOf[good[t]].push_back(t);
      }
    }
    // Each particle's best track: the most layers, the first such.
    std::unordered_map<std::uint64_t, std::uint32_t> best;
    for (const auto& [particle, owned] : tracksOf) {
      const auto top = *std::ranges::max_element(owned, {}, [&](std::uint32_t t) { return std::popcount(layersOf[t]); });
      best[particle] = top;
      for (const auto t : owned) {
        if (t != top) {
          ++mClones[length(tracks[t].size())];
        }
      }
    }
    // Every particle's rows, for the missing layers.
    auto rowsOf = study::rowsOfParticles(rows);
    const auto layers = study::particleLayers(rows, comparison);
    for (const auto& [particle, top] : best) {
      if (layersOf[top] == layers.at(particle).layers) {
        ++mCompleteTracks[length(tracks[top].size())];
      }
    }
    for (const auto& [id, mine] : layers) {
      const auto found = truth.find(id);
      if (found == truth.end() || !study::inPopulation(found->second, population) ||
          !study::isFindable(found->second.pt, mine, comparison)) {
        continue;
      }
      ++mFindable;
      const auto layerCount = static_cast<std::size_t>(std::popcount(mine.layers));
      const auto at = best.find(id);
      auto& record =
        mParticles.emplace_back(Particle{.id = id, .truth = found->second, .layers = mine.layers, .inbound = inbound(rows, rowsOf[id], found->second)});
      if (at == best.end()) {
        const bool withFake = std::ranges::any_of(rowsOf[id], [&](std::uint32_t row) {
          return std::ranges::any_of(tracksOfRow[row], [&](std::uint32_t t) {
            return good[t] == study::kNoTrackParticle && majority[t] == id;
          });
        });
        ++(withFake ? mFakeOnly : mNotFound);
        record.status = withFake ? kFakeOnly : kNotFound;
      } else {
        record.best = layersOf[at->second];
        record.bestIteration = iterations[at->second];
        record.goodTracks = static_cast<std::uint32_t>(tracksOf[id].size());
        if ((mine.layers & ~record.best) == 0U) {
          ++mComplete[std::min(layerCount, kMaxLayers - 1U)];
          record.status = kComplete;
          continue;
        }
        record.status = kIncomplete;
        ++mIncomplete[std::min(layerCount, kMaxLayers - 1U)];
        ++mMissingCount[std::min(static_cast<std::size_t>(std::popcount(mine.layers & ~record.best)), kLengths - 1U)];
      }
      // The missing layers -- every layer of a particle not found -- and,
      // against a best track, their place along the particle's path, by the
      // distance from the beam axis at its production z.
      const auto missing = mine.layers & ~record.best;
      const bool placed = at != best.end();
      const auto distance = [&](std::uint32_t row) {
        return std::hypot(static_cast<double>(rows[row].r), static_cast<double>(rows[row].z) - static_cast<double>(found->second.vz));
      };
      double first = std::numeric_limits<double>::infinity();
      double last = 0.;
      if (placed) {
        for (const auto row : tracks[at->second]) {
          first = std::min(first, distance(row));
          last = std::max(last, distance(row));
        }
      }
      std::array<std::uint8_t, kMaxLayers> fate{};
      std::array<double, kMaxLayers> where{};
      fate.fill(kFree);
      where.fill(std::numeric_limits<double>::infinity());
      for (const auto row : rowsOf[id]) {
        const auto layer = rows[row].layer;
        if (((missing >> layer) & 1U) == 0U) {
          continue;
        }
        where[layer] = std::min(where[layer], distance(row));
        for (const auto t : tracksOfRow[row]) {
          const auto kind = good[t] == id                        ? kSplit
                            : good[t] != study::kNoTrackParticle ? kOtherGood
                            : majority[t] == id                  ? kOwnFake
                                                                 : kFake;
          fate[layer] = std::min<std::uint8_t>(fate[layer], kind);
        }
      }
      for (std::size_t layer = 0; layer < kMaxLayers; ++layer) {
        if (((missing >> layer) & 1U) == 0U) {
          continue;
        }
        const auto bit = std::uint64_t{1} << layer;
        record.fates[fate[layer]] |= bit;
        if (placed) {
          const auto place = where[layer] < first ? kBefore : (where[layer] > last ? kAfter : kBetween);
          ++mMissing[layer][fate[layer]][place];
          record.before |= place == kBefore ? bit : 0U;
          record.after |= place == kAfter ? bit : 0U;
        }
      }
    }
  }

  /// The breakdown as one JSON line appended to \a path.
  void write(const std::string& path) const
  {
    std::ofstream out(path, std::ios::app);
    if (!out) {
      throw std::runtime_error("cannot append to " + path);
    }
    const auto list = [](const auto& values) {
      std::string text = "[";
      for (std::size_t i = 0; i < values.size(); ++i) {
        text += std::format("{}{}", i == 0U ? "" : ",", values[i]);
      }
      return text + "]";
    };
    std::string missing = "[";
    for (std::size_t layer = 0; layer < kMaxLayers; ++layer) {
      missing += layer == 0U ? "[" : ",[";
      for (std::size_t fate = 0; fate < kFates; ++fate) {
        missing += std::format("{}{}", fate == 0U ? "" : ",", list(mMissing[layer][fate]));
      }
      missing += "]";
    }
    missing += "]";
    out << std::format(R"({{"tracks":{},"fakes":{},"clones":{},"completeTracks":{},"findable":{},"notFound":{},"fakeOnly":{},"complete":{},"incomplete":{},)"
                       R"("missingCount":{},"missing":{}}})",
                       list(mTracks), list(mFakes), list(mClones), list(mCompleteTracks), mFindable, mNotFound, mFakeOnly, list(mComplete), list(mIncomplete),
                       list(mMissingCount), missing)
        << '\n';
  }

  /// The layers of \a mine (rows of \a truth's particle) whose rows are all
  /// past its turning point: nearer the axis than a row at least
  /// kTurnMarginCm before it in |z - vz|, by kTurnMarginCm.
  static std::uint64_t inbound(std::span<const study::TruthRow> rows, std::vector<std::uint32_t> mine, const study::ParticleTruth& truth)
  {
    constexpr float kTurnMarginCm = 0.5F;
    const auto along = [&](std::uint32_t row) { return std::abs(rows[row].z - truth.vz); };
    std::ranges::sort(mine, {}, along);
    std::uint64_t in = 0U;
    std::uint64_t out = 0U;
    float reached = 0.F;
    std::size_t earlier = 0U;
    for (const auto row : mine) {
      for (; earlier < mine.size() && along(mine[earlier]) < along(row) - kTurnMarginCm; ++earlier) {
        reached = std::max(reached, rows[mine[earlier]].r);
      }
      (rows[row].r < reached - kTurnMarginCm ? in : out) |= std::uint64_t{1} << rows[row].layer;
    }
    return in & ~out;
  }

  /// The class's particles as CSV rows appended to \a path, the header first
  /// when the file is new or empty (masks in hex).
  void writeParticles(const std::string& path) const
  {
    std::ifstream before(path);
    const bool fresh = !before.good() || before.peek() == std::ifstream::traits_type::eof();
    std::ofstream out(path, std::ios::app);
    if (!out) {
      throw std::runtime_error("cannot append to " + path);
    }
    if (fresh) {
      out << "id,pdg,pt,eta,phi,vx,vy,vz,layers,status,best,iteration,goodTracks,split,ownFake,fake,otherGood,free,before,after,inbound\n";
    }
    for (const auto& p : mParticles) {
      out << std::format("{},{},{:.9g},{:.9g},{:.9g},{:.9g},{:.9g},{:.9g},{:x},{},{:x},{},{},{:x},{:x},{:x},{:x},{:x},{:x},{:x},{:x}\n", p.id, p.truth.pdg,
                         p.truth.pt, p.truth.eta, p.truth.phi, p.truth.vx, p.truth.vy, p.truth.vz, p.layers, static_cast<int>(p.status), p.best,
                         static_cast<int>(p.bestIteration), p.goodTracks, p.fates[kSplit], p.fates[kOwnFake], p.fates[kFake], p.fates[kOtherGood],
                         p.fates[kFree], p.before, p.after, p.inbound);
    }
  }

 private:
  std::vector<Particle> mParticles;
  /// By track length (the last slot the longer ones too).
  std::array<std::uint64_t, kLengths> mTracks{};
  std::array<std::uint64_t, kLengths> mFakes{};
  std::array<std::uint64_t, kLengths> mClones{};
  std::array<std::uint64_t, kLengths> mCompleteTracks{};
  /// The class: not found, found only with a fake, and found complete or
  /// incomplete by the particle's layer count.
  std::uint64_t mFindable = 0U;
  std::uint64_t mNotFound = 0U;
  std::uint64_t mFakeOnly = 0U;
  std::array<std::uint64_t, kMaxLayers> mComplete{};
  std::array<std::uint64_t, kMaxLayers> mIncomplete{};
  /// The incomplete by their number of missing layers.
  std::array<std::uint64_t, kLengths> mMissingCount{};
  /// The missing layers by layer, fate and place.
  std::array<std::array<std::array<std::uint64_t, kPlaces>, kFates>, kMaxLayers> mMissing{};
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_TRACKBREAKDOWN_H_
