// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef GTRK_WORKFLOW_VERTEXLABEL_H_
#define GTRK_WORKFLOW_VERTEXLABEL_H_

#include "GenericTracking/Tracking/VertexStage.h"
#include "GenericTracking/Workflow/TrackLabel.h"

#include <array>
#include <map>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace o2::gtrk::workflow
{

struct VertexTruth {
  o2::MCCompLabel label;
  float purity = 0.F;
};

/// One vote per claimed, cluster-disjoint triplet, from the particle common
/// to its three measurements. Vote for its (source, event), ignoring track ID.
/// Mixed/noise/QED triplets cast no vote but remain in the purity denominator.
/// The largest vote wins (smallest source/event on ties); without a strict
/// majority its label is fake. No valid votes gives an unset label, purity 0.
/// RowLabels maps measurement rows to their original clusters' MC labels.
/// Claims must cover all contributors, before any adapter vertex reordering.
template <class RowLabels>
[[nodiscard]] std::vector<VertexTruth> vertexTruth(std::span<const TripletClaim> claims, std::span<const std::uint32_t> contributors, const RowLabels& rowLabels)
{
  using Collision = std::pair<int, int>; // source, event
  std::vector<std::map<Collision, std::uint32_t>> votes(contributors.size());
  std::vector<std::uint32_t> counts(contributors.size(), 0U);
  for (const auto& claim : claims) {
    if (claim.vertex < 0) {
      continue;
    }
    const auto v = static_cast<std::size_t>(claim.vertex);
    if (v >= contributors.size()) {
      throw std::runtime_error("vertex MC claim refers to a missing vertex");
    }
    ++counts[v];
    std::array<RowParticles, 3> rows;
    for (std::size_t i = 0; i < rows.size(); ++i) {
      rows[i] = rowParticles(rowLabels(claim.rows[i]));
    }
    constexpr std::array<std::uint32_t, 3> order{0U, 1U, 2U};
    if (const auto particle = goodParticle(rows, order); particle != kNoTrackParticle) {
      const auto label = labelOfParticle(particle);
      ++votes[v][{label.getSourceID(), label.getEventID()}];
    }
  }
  std::vector<VertexTruth> result(contributors.size());
  for (std::size_t v = 0; v < result.size(); ++v) {
    if (counts[v] != contributors[v]) {
      throw std::runtime_error("vertex MC claims do not match its contributor count");
    }
    std::uint32_t best = 0U;
    for (const auto& [collision, count] : votes[v]) {
      if (count > best) {
        best = count;
        result[v].label = {o2::MCCompLabel::maxTrackID(), collision.second, collision.first, count <= contributors[v] / 2U};
      }
    }
    result[v].purity = contributors[v] == 0U ? 0.F : static_cast<float>(best) / static_cast<float>(contributors[v]);
  }
  return result;
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_VERTEXLABEL_H_
