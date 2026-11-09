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

/// \file TrackLabel.h
/// \brief The one rule that gives a track its MC particle: what a measurement's
///        labels carry, and which particle a track's measurements carry.
///
/// The study scores tracks with it and the O2-format adaptors label them with
/// it, so the two cannot drift apart.

#ifndef GTRK_WORKFLOW_TRACKLABEL_H_
#define GTRK_WORKFLOW_TRACKLABEL_H_

#include "SimulationDataFormat/MCCompLabel.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>

namespace o2::gtrk::workflow
{

/// Whether \a label names a particle the truth counts: O2's valid, which is
/// set and not noise, the QED background counting as noise.
[[nodiscard]] inline bool carriesParticle(const o2::MCCompLabel& label) { return label.isValid(); }

/// The particles one measurement's cluster carries, each identified by its
/// source, event and track in one word (`MCCompLabel::getTrackEventSourceID`).
struct RowParticles {
  static constexpr std::size_t kMaxParticles = 4;
  std::array<std::uint64_t, kMaxParticles> particles{}; ///< the first `count` are valid
  std::uint8_t count = 0;
  /// Whether `particles[0]` is the cluster's first label and O2's QA counts it
  /// (`CheckTracksCA.C`: valid, correct, from the signal): O2's findable
  /// reads only that one.
  bool firstLabelCounts = false;

  [[nodiscard]] bool carries(std::uint64_t particle) const { return std::find(particles.begin(), particles.begin() + count, particle) != particles.begin() + count; }
};

/// The particles of a cluster with the labels \a labels, in their order: those
/// that carry one (`carriesParticle`), each once, the first kMaxParticles.
[[nodiscard]] inline RowParticles rowParticles(std::span<const o2::MCCompLabel> labels)
{
  RowParticles row;
  row.firstLabelCounts = !labels.empty() && labels[0].isValid() && labels[0].isCorrect() && labels[0].getSourceID() == 0;
  for (const auto& label : labels) {
    if (!carriesParticle(label)) {
      continue;
    }
    const auto id = label.getTrackEventSourceID();
    if (row.carries(id) || row.count == RowParticles::kMaxParticles) {
      continue;
    }
    row.particles[row.count++] = id;
  }
  return row;
}

/// The label of \a particle, a `getTrackEventSourceID()` word, with the
/// fake flag \a fake: the inverse of the word's packing.
[[nodiscard]] inline o2::MCCompLabel labelOfParticle(std::uint64_t particle, bool fake = false)
{
  using L = o2::MCCompLabel;
  return {static_cast<int>(particle & L::maskTrackID), static_cast<int>((particle >> L::nbitsTrackID) & L::maskEvID), static_cast<int>((particle >> (L::nbitsTrackID + L::nbitsEvID)) & L::maskSrcID), fake};
}

/// No particle: a fake's good particle, or noise only.
inline constexpr std::uint64_t kNoTrackParticle = std::numeric_limits<std::uint64_t>::max();

/// The particle every row of \a track carries, the first such among the
/// first row's; kNoTrackParticle for a fake. \a rows are indexed by the
/// track's rows and hold RowParticles' members.
template <typename Rows>
[[nodiscard]] std::uint64_t goodParticle(const Rows& rows, std::span<const std::uint32_t> track)
{
  if (track.empty()) {
    return kNoTrackParticle;
  }
  const auto& first = rows[track.front()];
  for (std::uint8_t i = 0; i < first.count; ++i) {
    const auto candidate = first.particles[i];
    if (std::ranges::all_of(track, [&](std::uint32_t row) { return rows[row].carries(candidate); })) {
      return candidate;
    }
  }
  return kNoTrackParticle;
}

/// The particle most rows of \a track carry, the first such in track order,
/// and on how many: a fake's label; kNoTrackParticle and 0 on noise only.
struct MajorityParticle {
  std::uint64_t particle = kNoTrackParticle;
  std::uint32_t carried = 0U;
};

template <typename Rows>
[[nodiscard]] MajorityParticle majorityParticle(const Rows& rows, std::span<const std::uint32_t> track)
{
  MajorityParticle result;
  for (const auto row : track) {
    for (std::uint8_t i = 0; i < rows[row].count; ++i) {
      const auto candidate = rows[row].particles[i];
      const auto carried = static_cast<std::uint32_t>(std::ranges::count_if(track, [&](std::uint32_t other) { return rows[other].carries(candidate); }));
      if (carried > result.carried) {
        result = {.particle = candidate, .carried = carried};
      }
    }
  }
  return result;
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_TRACKLABEL_H_
