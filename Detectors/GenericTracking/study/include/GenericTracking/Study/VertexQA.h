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

/// \file VertexQA.h
/// \brief The own vertexer's vertices against the MC collisions, under the
///        vertex study's rules (the bench's vertex study, VertexStudy.h).
///
/// A collision is a generator event, its z the median of its measured
/// primaries' production z (the generator's primaries include displaced
/// decay daughters, so no single one will do). A triplet of the pass belongs
/// to a collision when one primary of it carries all three measurements. A
/// vertex belongs to the collision most of its claimed triplets belong to,
/// and is a fake when background -- no such primary -- is the most. A
/// collision is found when a vertex belongs to it, the one with the most
/// contributors; any other is a split. A collision without one whose
/// triplets went mostly to another collision's vertex is merged. What a
/// missed collision costs the tracking is its findable primaries: 0.1 GeV or
/// more on a chain of four layers, the sink's findable (TrackComparison.h).

#ifndef GTRK_STUDY_VERTEXQA_H_
#define GTRK_STUDY_VERTEXQA_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/VertexStage.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/Statistics.h"
#include "GenericTracking/Study/TrackComparison.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <map>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::gtrk::study
{

/// One timeframe's vertices against its collisions.
struct VertexScore {
  std::size_t vertices = 0U;
  std::size_t collisions = 0U; ///< with a measured primary
  std::size_t found = 0U;
  std::uint64_t findable = 0U; ///< the collisions' findable primaries
  std::uint64_t covered = 0U;  ///< those of the found collisions
  std::size_t fakes = 0U;
  std::size_t splits = 0U;
  std::size_t merged = 0U;
  std::vector<double> residualsUM; ///< |z - z_true| of each found collision's vertex [um]
  std::vector<double> pulls;       ///< (z - z_true) / its sigma z, the same vertices
  /// A collision no vertex is the best of, and where its triplets went.
  struct Missed {
    double z = 0.; ///< the median production z of its primaries [cm]
    std::uint32_t findable = 0U;
    std::uint32_t triplets = 0U; ///< the pass's triplets of one of its primaries, in the slices owning their frames
    std::uint32_t cut = 0U;      ///< of those, none (kNone: cut, degenerate or deduplicated)
    std::uint32_t free = 0U;     ///< a candidate no vertex claimed (kFree)
    std::uint32_t other = 0U;    ///< claimed by a vertex another collision (or none) is the majority of
  };
  std::vector<Missed> missed;
};

/// The shared majority, best-vertex, split and merge rule. Ordered maps
/// make the first collision/vertex win a tie; equal contributors keep the
/// first vertex. Negative labels denote background or unclaimed triplets.
struct VertexClaimsScore {
  std::vector<std::int64_t> collision;
  std::vector<std::int64_t> best;
  std::uint32_t fakes = 0U;
  std::uint32_t splits = 0U;
  std::uint32_t merged = 0U;
};

[[nodiscard]] inline VertexClaimsScore scoreVertexClaims(std::span<const std::map<std::int64_t, std::uint32_t>> labels, std::span<const std::uint32_t> contributors, std::span<const std::map<std::int64_t, std::uint32_t>> wentTo)
{
  VertexClaimsScore result;
  result.best.assign(wentTo.size(), -1);
  for (std::size_t v = 0; v < labels.size(); ++v) {
    const auto collision = labels[v].empty() ? -1 : std::ranges::max_element(labels[v], {}, [](const auto& entry) { return entry.second; })->first;
    result.collision.push_back(collision);
    if (collision < 0) {
      ++result.fakes;
      continue;
    }
    auto& top = result.best[static_cast<std::size_t>(collision)];
    result.splits += top >= 0 ? 1U : 0U;
    if (top < 0 || contributors[static_cast<std::size_t>(top)] < contributors[v]) {
      top = static_cast<std::int64_t>(v);
    }
  }
  for (std::size_t c = 0; c < wentTo.size(); ++c) {
    if (result.best[c] < 0 && !wentTo[c].empty()) {
      const auto top = std::ranges::max_element(wentTo[c], {}, [](const auto& entry) { return entry.second; })->first;
      result.merged += top >= 0 ? 1U : 0U;
    }
  }
  return result;
}

/// \a found's vertices against the collisions of \a rows (one per frame row)
/// and \a particles, from the pass's \a claims (`VertexingPass::Result`).
template <DetectorLayout Det>
[[nodiscard]] VertexScore scoreVertices(std::span<const TruthRow> rows, const ParticleTable& particles, std::span<const TripletClaim> claims, const VertexerResult& found)
{
  // Per particle, its layers and longest chain: the findable's rule.
  const auto settings = comparisonSettings<Det>();
  const auto layers = particleLayers(rows, settings);
  struct Collision {
    std::vector<double> productionZ;
    std::uint32_t findable = 0U;
  };
  // The collisions with a measured primary, by `collisionOf` (the vertex
  // study's index takes every generator primary's).
  std::unordered_map<std::uint64_t, std::size_t> index;
  std::vector<Collision> collisions;
  for (const auto& [particle, info] : particles) {
    const auto measured = layers.find(particle);
    if (!info.primary || measured == layers.end()) {
      continue;
    }
    const auto [at, inserted] = index.try_emplace(collisionOf(particle), collisions.size());
    if (inserted) {
      collisions.emplace_back();
    }
    auto& collision = collisions[at->second];
    collision.productionZ.push_back(static_cast<double>(info.vz));
    collision.findable += isFindable(info.pt, measured->second, settings) ? 1U : 0U;
  }
  // Each triplet's collision, -1 for background.
  const auto tripletCollision = [&](const TripletClaim& claim) -> std::int64_t {
    const auto particle = commonParticle(rows, claim.rows);
    if (!particle) {
      return -1;
    }
    const auto info = particles.find(*particle);
    return info != particles.end() && info->second.primary ? static_cast<std::int64_t>(index.at(collisionOf(*particle))) : -1;
  };
  VertexScore score;
  score.vertices = found.vertices.size();
  score.collisions = collisions.size();
  std::vector<std::map<std::int64_t, std::uint32_t>> labels(found.vertices.size());
  std::vector<std::map<std::int64_t, std::uint32_t>> wentTo(collisions.size()); // per collision: the vertices its triplets went to, -1 free, -2 cut
  for (const auto& claim : claims) {
    const auto vertex = claim.vertex;
    const auto collision = tripletCollision(claim);
    if (vertex >= 0) {
      ++labels[static_cast<std::size_t>(vertex)][collision];
    }
    if (collision >= 0) {
      ++wentTo[static_cast<std::size_t>(collision)][vertex];
    }
  }
  const auto assigned = scoreVertexClaims(labels, found.contributors, wentTo);
  const auto& best = assigned.best;
  score.fakes = assigned.fakes;
  score.splits = assigned.splits;
  score.merged = assigned.merged;
  for (std::size_t c = 0; c < collisions.size(); ++c) {
    auto& collision = collisions[c];
    score.findable += collision.findable;
    if (best[c] < 0) {
      VertexScore::Missed miss{.z = median(collision.productionZ), .findable = collision.findable};
      for (const auto& [vertex, count] : wentTo[c]) {
        if (vertex == VertexerResult::kNotOwned) {
          continue; // a halo's copy, counted in the slice that owns its frame
        }
        miss.triplets += count;
        (vertex == VertexerResult::kNone ? miss.cut : vertex == VertexerResult::kFree ? miss.free
                                                                                      : miss.other) += count;
      }
      score.missed.push_back(miss);
      continue;
    }
    ++score.found;
    score.covered += collision.findable;
    const double z = median(collision.productionZ); // a collision has a primary: never empty
    const auto& vertex = found.vertices[static_cast<std::size_t>(best[c])];
    score.residualsUM.push_back(std::abs((static_cast<double>(vertex.z) - z) * 1.e4));
    if (vertex.varZ > 0.F) {
      score.pulls.push_back((static_cast<double>(vertex.z) - z) / std::sqrt(static_cast<double>(vertex.varZ)));
    }
  }
  return score;
}

/// \a score as the sink logs it.
[[nodiscard]] inline std::string describe(const VertexScore& score)
{
  auto residuals = score.residualsUM;
  std::ranges::sort(residuals);
  const auto at = [&](double q) { return residuals.empty() ? 0. : residuals[std::min(residuals.size() - 1U, static_cast<std::size_t>(q * static_cast<double>(residuals.size())))]; };
  return std::format(
    "vertex QA: {} vertices; collisions found {}/{}, their findable primaries {}/{} ({:.2f} %); fakes {}, splits {}, merged {}; "
    "|dz| 68/95/99/max [um] {:.0f} {:.0f} {:.0f} {:.0f}",
    score.vertices, score.found, score.collisions, score.covered, score.findable, 100. * static_cast<double>(score.covered) / std::max(1., static_cast<double>(score.findable)), score.fakes, score.splits, score.merged, at(0.68), at(0.95), at(0.99), residuals.empty() ? 0. : residuals.back());
}

/// \a score's missed collisions by their findable primaries: per bin the
/// collisions, their findable primaries and where their triplets went.
[[nodiscard]] inline std::string describeMissed(const VertexScore& score)
{
  constexpr std::array<std::uint32_t, 7> kLow{0U, 1U, 5U, 10U, 20U, 50U, 100U};
  std::string line = "vertex QA, the missed collisions by findable primaries (collisions, findable; their triplets: none / free / another's):";
  for (std::size_t bin = 0; bin < kLow.size(); ++bin) {
    const auto high = bin + 1U < kLow.size() ? kLow[bin + 1U] : std::numeric_limits<std::uint32_t>::max();
    std::uint64_t collisions = 0U;
    std::uint64_t findable = 0U;
    std::uint64_t cut = 0U;
    std::uint64_t free = 0U;
    std::uint64_t other = 0U;
    for (const auto& miss : score.missed) {
      if (miss.findable >= kLow[bin] && miss.findable < high) {
        ++collisions;
        findable += miss.findable;
        cut += miss.cut;
        free += miss.free;
        other += miss.other;
      }
    }
    line += std::format(" [{}{}] {}, {}; {} / {} / {}", kLow[bin], high == std::numeric_limits<std::uint32_t>::max() ? "+" : std::format("-{}", high - 1U), collisions, findable, cut, free, other);
  }
  return line;
}

/// \a score's pulls, (z - z_true) / sigma z of each found collision's vertex,
/// their mean and rms and every one of them, for sums over timeframes.
[[nodiscard]] inline std::string describePulls(const VertexScore& score)
{
  double sum = 0.;
  double sum2 = 0.;
  for (const double pull : score.pulls) {
    sum += pull;
    sum2 += pull * pull;
  }
  const auto n = static_cast<double>(std::max<std::size_t>(score.pulls.size(), 1U));
  std::string line = std::format("vertex pulls: {} vertices, mean {:+.3f}, rms {:.3f}:", score.pulls.size(), sum / n, std::sqrt(sum2 / n));
  for (const double pull : score.pulls) {
    line += std::format(" {:.3f}", pull);
  }
  return line;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_VERTEXQA_H_
