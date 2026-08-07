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

// StateFamily describes the state representation associated with a
// SurfaceDescriptor::kind. It is local to state-representation operations;
// topology and traversal code use SurfaceKind directly.
namespace o2::itsmft::tracking
{

enum class StateFamily : uint8_t {
  Invalid,
  Barrel,
  Forward
};

/// Maps a SurfaceKind to its state representation. It is a representation tag,
/// never topology or dispatch policy; derive it from SurfaceKind.
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
