// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SEEDANCHOR_H_
#define ALICEO2_ITSMFT_TRACKING_SEEDANCHOR_H_

#include <cstdint>
#include <type_traits>

// M4b: SeedAnchor is an unrelated legacy concept that was previously bundled
// into the state-family classification header by historical accident (both
// used to live in the single pre-M4 TransitionPolicy.h). It has nothing to
// do with the legacy hot-loop-dispatch tag's containment or with
// state-family/SurfaceKind classification: it never names a surface kind, a
// state family, or that tag.
//
// M5a classification (reviewed, not moved): SeedAnchor is retained in this
// public, non-detail header only because it is the parameter type of a
// public generic leaf operation's own declared signature --
// barrel::/forward::buildSeed(SeedAnchor, ...), Sec 4 of
// doc/design/0001-descriptor-driven-operation-boundary.md. It is not itself
// a true generic operation *concept* the way SurfaceDescriptor/SurfaceKind
// are: it is a temporary Stage-B refit-primitive-slice selector, needed only
// because buildSeed's `Inner` case exists solely to reproduce the frozen
// legacy ITS refit's reverse anchor convention (see the enum doc below); its
// own anchor-taking overload has no live production caller today (the
// production path always uses the anchor-less overload, which internally
// defaults to `Outer`). It stays here, unmoved, because no correctness issue
// requires relocating it now and the M5a design note above classifies the
// operation it parametrizes as an irreducible family-local leaf, not
// duplicated orchestration. If a later containment slice determines
// buildSeed's `Inner` case itself is legacy-only (i.e. no descriptor-driven
// contract ever calls it with `Inner`), relocate SeedAnchor to detail/ at
// that point instead of leaving it here as a false generic-concept signal.
namespace o2::itsmft::tracking
{

// Selects which of the three hits in a {inner, middle, outer} candidate
// supplies the anchor/reference frame for seed construction (Stage-B
// refit-primitive slice, design report Sec 5). Explicit values are locked
// (never renumbered) because they may be threaded through device-facing
// call sites in a later slice. Outer is the current accepted
// buildCellSeed/buildSeed anchor (referenceCoordinate/alpha/covariance
// come from the outer measurement's own tracking frame). Inner is the
// frozen ITS `reverse=true` anchor
// (o2::its::track::buildTrackSeed/seedTrackForRefit,
// ITStracking/TrackHelpers.h): referenceCoordinate/alpha/covariance come
// from the inner measurement's own tracking frame instead, with the
// legacy sign flip applied to snp/q2pt/tgl so the local direction
// convention stays consistent with the swapped anchor. This is a plain
// selector, not a reverse-traversal flag: it never encodes propagation
// direction, material-correction direction, or fit-leg order by itself.
enum class SeedAnchor : uint8_t {
  Inner = 0,
  Outer = 1
};

static_assert(std::is_standard_layout_v<SeedAnchor> && std::is_trivially_copyable_v<SeedAnchor>);
static_assert(std::is_same_v<std::underlying_type_t<SeedAnchor>, uint8_t>);
static_assert(sizeof(SeedAnchor) == sizeof(uint8_t));
static_assert(static_cast<uint8_t>(SeedAnchor::Inner) == 0);
static_assert(static_cast<uint8_t>(SeedAnchor::Outer) == 1);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SEEDANCHOR_H_ */
