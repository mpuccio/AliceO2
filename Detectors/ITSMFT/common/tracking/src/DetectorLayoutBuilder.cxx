// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file DetectorLayoutBuilder.cxx
/// \brief Explicit-order construction of a DetectorLayout from a surface catalog
///

#include "ITSMFTTracking/DetectorLayoutBuilder.h"

namespace o2::itsmft::tracking
{

namespace
{
struct BuiltTransition {
  TransitionId id{};
  SurfaceId from{};
  SurfaceId to{};
  SurfaceMask skipped{};
};
} // namespace

DetectorLayoutBuildResult DetectorLayoutBuilder::build() const
{
  DetectorLayoutBuildResult result;

  // A catalog larger than the 32-surface limit cannot be represented by
  // SurfaceMask at all: bits for ids >= MaxLayoutSurfaces are silently
  // dropped by SurfaceMask::set/has rather than reported. Reject it here,
  // before any subgraph validation below relies on mask arithmetic being
  // faithful, and report it through the same TopologyBuildError value
  // SparseTrackingTopology itself would raise for the same condition.
  if (mCatalog.nSurfaces > MaxLayoutSurfaces) {
    result.error = DetectorLayoutBuildError::TopologyRejected;
    result.topologyError = TopologyBuildError::InvalidSurfaceCount;
    return result;
  }

  // The builder's own subgraph-membership checks below key off `id.value() <
  // mCatalog.nSurfaces`, which is only meaningful if catalog entries are
  // dense (catalog[i].id == i). DetectorLayout::validate() checks this too,
  // but only after the topology is fully built; check it up front so a
  // non-dense catalog fails fast with the same diagnostic it would
  // eventually get.
  for (uint32_t i = 0; i < mCatalog.nSurfaces; ++i) {
    if (mCatalog.surfaces[i].id != SurfaceId{static_cast<uint16_t>(i)}) {
      result.error = DetectorLayoutBuildError::LayoutRejected;
      result.layoutError = DetectorLayoutError::NonDenseSurfaceIds;
      return result;
    }
  }

  SurfaceMask seenAcrossSubgraphs{};
  SurfaceMask combinedSeeding{};
  for (const auto& subgraph : mSubgraphs) {
    if (subgraph.orderedSurfaces.empty()) {
      result.error = DetectorLayoutBuildError::EmptySubgraph;
      return result;
    }
    if (subgraph.maxHoles < 0) {
      result.error = DetectorLayoutBuildError::NegativeMaxHoles;
      return result;
    }

    // Every surface in a subgraph must share one SurfaceKind (ADR 0007
    // decision 8: the subgraph's family is derived from its own catalog
    // entries, never asserted by the caller as a separate policy tag).
    // Checked explicitly rather than left to fall out of addTransition(): a
    // singleton subgraph never calls addTransition at all, so this would
    // otherwise go unreported. `expectedKind` is taken from the first
    // surface once its id has been validated below -- never read before
    // that, since an out-of-range first id must fail
    // InvalidSubgraphSurfaceId, not an out-of-bounds catalog access.
    SurfaceMask subgraphSurfaces{};
    std::optional<SurfaceKind> expectedKind{};
    for (const auto& id : subgraph.orderedSurfaces) {
      if (!id.isValid() || id.value() >= mCatalog.nSurfaces) {
        result.error = DetectorLayoutBuildError::InvalidSubgraphSurfaceId;
        return result;
      }
      if (subgraphSurfaces.has(id)) {
        result.error = DetectorLayoutBuildError::DuplicateSurfaceInSubgraph;
        return result;
      }
      if (seenAcrossSubgraphs.has(id)) {
        result.error = DetectorLayoutBuildError::SurfaceDuplicatedAcrossSubgraphs;
        return result;
      }
      const auto kind = mCatalog.surfaces[id.value()].kind;
      if (!expectedKind.has_value()) {
        expectedKind = kind;
      } else if (kind != *expectedKind) {
        result.error = DetectorLayoutBuildError::LayoutRejected;
        result.layoutError = DetectorLayoutError::PolicySurfaceKindMismatch;
        return result;
      }
      subgraphSurfaces.set(id);
    }
    seenAcrossSubgraphs |= subgraphSurfaces;

    if (!subgraph.holeSurfaces.isSubsetOf(subgraphSurfaces)) {
      result.error = DetectorLayoutBuildError::HoleSurfacesOutsideSubgraph;
      return result;
    }
    if (!subgraph.seedingSurfaces.isSubsetOf(subgraphSurfaces)) {
      result.error = DetectorLayoutBuildError::SeedingSurfacesOutsideSubgraph;
      return result;
    }
    combinedSeeding |= subgraph.seedingSurfaces;
  }

  SparseTrackingTopology topology{mCatalog.nSurfaces, combinedSeeding};
  if (topology.getError() != TopologyBuildError::None) {
    result.error = DetectorLayoutBuildError::TopologyRejected;
    result.topologyError = topology.getError();
    return result;
  }

  for (const auto& subgraph : mSubgraphs) {
    const auto& ordered = subgraph.orderedSurfaces;
    std::vector<BuiltTransition> transitions;

    // Transitions and their skipped-surface masks are derived purely from
    // positions in `ordered`, so a permutation of global ids that does not
    // increase along the traversal is handled identically to a monotonic
    // one (see the non-monotonic-chain test).
    for (size_t posFrom = 0; posFrom < ordered.size(); ++posFrom) {
      for (size_t posTo = posFrom + 1; posTo < ordered.size(); ++posTo) {
        SurfaceMask skipped{};
        for (size_t k = posFrom + 1; k < posTo; ++k) {
          skipped.set(ordered[k]);
        }
        if (skipped.count() > subgraph.maxHoles || !skipped.isSubsetOf(subgraph.holeSurfaces)) {
          continue;
        }
        const auto id = topology.addTransition(SurfaceTransition{ordered[posFrom], ordered[posTo], skipped, 0});
        if (!id.isValid()) {
          result.error = DetectorLayoutBuildError::TopologyRejected;
          result.topologyError = topology.getError();
          return result;
        }
        transitions.push_back(BuiltTransition{id, ordered[posFrom], ordered[posTo], skipped});
      }
    }

    for (const auto& first : transitions) {
      for (const auto& second : transitions) {
        if (first.to != second.from) {
          continue;
        }
        const auto combinedSkipped = first.skipped | second.skipped;
        if (combinedSkipped.count() > subgraph.maxHoles) {
          continue;
        }
        const auto cellId = topology.addCell(first.id, second.id);
        if (!cellId.isValid()) {
          result.error = DetectorLayoutBuildError::TopologyRejected;
          result.topologyError = topology.getError();
          return result;
        }
      }
    }
  }

  if (!topology.finalize()) {
    result.error = DetectorLayoutBuildError::TopologyRejected;
    result.topologyError = topology.getError();
    return result;
  }

  DetectorLayout layout{gsl::span<const SurfaceDescriptor>{mCatalog.surfaces, mCatalog.nSurfaces}, std::move(topology)};
  if (!layout.valid()) {
    result.error = DetectorLayoutBuildError::LayoutRejected;
    result.layoutError = layout.getError();
    return result;
  }

  result.layout.emplace(std::move(layout));
  return result;
}

} // namespace o2::itsmft::tracking
