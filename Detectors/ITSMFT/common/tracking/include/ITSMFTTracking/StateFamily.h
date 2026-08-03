// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_STATEFAMILY_H_
#define ALICEO2_ITSMFT_TRACKING_STATEFAMILY_H_

#include <cstdint>
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"

// ADR 0007 decisions 7-8 / M4: StateFamily and SurfaceKind are the permanent
// public generic-core concepts; the legacy hot-loop-dispatch tag those
// decisions confine behind the detail/ boundary is a temporary
// implementation detail, not a future public/core abstraction. This header
// lets any generic/adapter-facing code classify a surface's family without
// naming that legacy tag at all.
namespace o2::itsmft::tracking
{

enum class StateFamily : uint8_t {
  Invalid,
  Barrel,
  Forward
};

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

/// The state family every surface of `kind` belongs to. This is the public,
/// tag-free bridge between SurfaceKind and StateFamily (ADR 0007 decision 8):
/// generic/adapter code classifies a surface's family directly from its own
/// SurfaceDescriptor::kind, never by naming the legacy hot-loop-dispatch tag.
GPUhdi() constexpr StateFamily stateFamilyOf(SurfaceKind kind) noexcept
{
  switch (kind) {
    case SurfaceKind::Cylinder:
      return StateFamily::Barrel;
    case SurfaceKind::Disk:
      return StateFamily::Forward;
  }
  return StateFamily::Invalid;
}

static_assert(std::is_same_v<std::underlying_type_t<StateFamily>, uint8_t>);
static_assert(sizeof(StateFamily) == sizeof(uint8_t));
static_assert(stateFamilyOf(SurfaceKind::Cylinder) == StateFamily::Barrel);
static_assert(stateFamilyOf(SurfaceKind::Disk) == StateFamily::Forward);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_STATEFAMILY_H_ */
