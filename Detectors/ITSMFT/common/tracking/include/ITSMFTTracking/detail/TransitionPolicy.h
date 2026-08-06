// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_DETAIL_TRANSITIONPOLICY_H_
#define ALICEO2_ITSMFT_TRACKING_DETAIL_TRANSITIONPOLICY_H_

// M4 (GenericTrackingEngineMigration.md; ADR 0007 decisions 7-8): this header
// -- and every other header under ITSMFTTracking/detail/ -- is a temporary
// legacy hot-loop-dispatch implementation detail, not a public/adapter-facing
// API. TransitionPolicyTag must never be nameable outside this detail
// boundary; the permanent public concepts (SurfaceKind, StateFamily) live in
// ITSMFTTracking/SurfaceDescriptor.h / ITSMFTTracking/StateFamily.h.

#include <cstdint>
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/StateFamily.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"

namespace o2::itsmft::tracking
{

enum class TransitionPolicyTag : uint16_t {
  Invalid,
  CylinderCylinder,
  DiskDisk
};

GPUhdi() constexpr bool isKnownTransitionPolicyTag(TransitionPolicyTag tag) noexcept
{
  switch (tag) {
    case TransitionPolicyTag::Invalid:
    case TransitionPolicyTag::CylinderCylinder:
    case TransitionPolicyTag::DiskDisk:
      return true;
  }
  return false;
}

GPUhdi() constexpr bool isStageATransitionPolicyTagEnabled(TransitionPolicyTag tag) noexcept
{
  return tag == TransitionPolicyTag::CylinderCylinder || tag == TransitionPolicyTag::DiskDisk;
}

GPUhdi() constexpr StateFamily stateFamilyOf(TransitionPolicyTag tag) noexcept
{
  switch (tag) {
    case TransitionPolicyTag::CylinderCylinder:
      return StateFamily::Barrel;
    case TransitionPolicyTag::DiskDisk:
      return StateFamily::Forward;
    case TransitionPolicyTag::Invalid:
      return StateFamily::Invalid;
  }
  return StateFamily::Invalid;
}

/// True if a transition carrying `tag` may legally connect surfaces of
/// `kind`. This is the single shared definition of the policy/surface-kind
/// compatibility rule; SurfaceGraph construction-time validation and the
/// policy/state traits both call this instead of each encoding their own copy.
GPUhdi() constexpr bool isSurfaceKindCompatible(TransitionPolicyTag tag, SurfaceKind kind) noexcept
{
  switch (tag) {
    case TransitionPolicyTag::CylinderCylinder:
      return kind == SurfaceKind::Cylinder;
    case TransitionPolicyTag::DiskDisk:
      return kind == SurfaceKind::Disk;
    case TransitionPolicyTag::Invalid:
      return false;
  }
  return false;
}

/// The unique TransitionPolicyTag whose isSurfaceKindCompatible() accepts
/// `kind` (current same-family-only topology, ADR 0007 decision 8): the
/// inverse of isSurfaceKindCompatible, used to derive a transition's tag from
/// its endpoint SurfaceDescriptor::kind now that SurfaceTransition no longer
/// stores the tag itself (M4). Every current production transition's kind
/// and tag were already consistent by construction (SurfaceGraphBuilder
/// only ever built same-kind subgraphs), so this derivation reproduces
/// exactly the value the removed stored field always held.
GPUhdi() constexpr TransitionPolicyTag transitionPolicyTagForSurfaceKind(SurfaceKind kind) noexcept
{
  switch (kind) {
    case SurfaceKind::Cylinder:
      return TransitionPolicyTag::CylinderCylinder;
    case SurfaceKind::Disk:
      return TransitionPolicyTag::DiskDisk;
  }
  return TransitionPolicyTag::Invalid;
}

static_assert(std::is_same_v<std::underlying_type_t<TransitionPolicyTag>, uint16_t>);
static_assert(sizeof(TransitionPolicyTag) == sizeof(uint16_t));
static_assert(transitionPolicyTagForSurfaceKind(SurfaceKind::Cylinder) == TransitionPolicyTag::CylinderCylinder);
static_assert(transitionPolicyTagForSurfaceKind(SurfaceKind::Disk) == TransitionPolicyTag::DiskDisk);
static_assert(isSurfaceKindCompatible(transitionPolicyTagForSurfaceKind(SurfaceKind::Cylinder), SurfaceKind::Cylinder));
static_assert(isSurfaceKindCompatible(transitionPolicyTagForSurfaceKind(SurfaceKind::Disk), SurfaceKind::Disk));

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_DETAIL_TRANSITIONPOLICY_H_ */
