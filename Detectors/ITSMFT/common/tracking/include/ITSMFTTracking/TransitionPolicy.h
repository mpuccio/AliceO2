// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICY_H_
#define ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICY_H_

#include <cstdint>
#include <type_traits>

#include "GPUCommonDef.h"

namespace o2::itsmft::tracking
{

enum class StateFamily : uint8_t {
  Invalid,
  Barrel,
  Forward
};

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

static_assert(std::is_same_v<std::underlying_type_t<StateFamily>, uint8_t>);
static_assert(std::is_same_v<std::underlying_type_t<TransitionPolicyTag>, uint16_t>);
static_assert(sizeof(StateFamily) == sizeof(uint8_t));
static_assert(sizeof(TransitionPolicyTag) == sizeof(uint16_t));

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICY_H_ */
