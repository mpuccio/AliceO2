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

// ADR 0007 decisions 7-8 / M4-M4b: StateFamily and SurfaceKind are the
// permanent public generic-core concepts; the legacy hot-loop-dispatch tag
// those decisions confine behind the detail/ boundary is a temporary
// implementation detail, not a future public/core abstraction. This header
// lets state-representation types (SurfaceKinematicState,
// SurfaceLinearizationReference) and the barrel::/forward:: operations that
// consume them classify a surface's family without naming that legacy tag at
// all.
//
// M4b correction: StateFamily is a state-representation-local concept only.
// It must never be used as a substitute transition-policy/dispatch key by
// public topology, traversal, scheduling, binding-count, or adapter-facing
// observability APIs -- those must use SurfaceDescriptor/SurfaceKind
// directly, or (for anything still dispatched by the legacy hot-loop tag)
// stay confined behind the same detail/ boundary as that tag. See
// detail::stateFamilyFromNLayers() (detail/TransitionPolicyState.h) for the
// one legacy NLayers-to-family dispatch bridge this correction moved out of
// the public Cell.h.
namespace o2::itsmft::tracking
{

enum class StateFamily : uint8_t {
  Invalid,
  Barrel,
  Forward
};

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
