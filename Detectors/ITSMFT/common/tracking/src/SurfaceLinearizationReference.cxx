// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/SurfaceKinematicState.h"

#include <cmath>

namespace o2::itsmft::tracking
{

bool makeLinearizationReference(const SurfaceKinematicState& state, SurfaceLinearizationReference& out) noexcept
{
  if (!state.hasRecognizedKind()) {
    return false;
  }
  SurfaceLinearizationReference scratch{};
  for (uint8_t i = 0; i < 5; ++i) {
    if (!std::isfinite(state.parameters[i])) {
      return false;
    }
    scratch.parameters[i] = state.parameters[i];
  }
  if (!std::isfinite(state.referenceCoordinate) || !std::isfinite(state.alpha)) {
    return false;
  }
  scratch.referenceCoordinate = state.referenceCoordinate;
  scratch.alpha = state.alpha;
  scratch.kind = state.kind;
  out = scratch;
  return true;
}

} // namespace o2::itsmft::tracking
