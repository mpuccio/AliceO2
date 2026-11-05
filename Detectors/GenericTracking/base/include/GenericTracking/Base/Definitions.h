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

/// \file Definitions.h
/// \brief Index types and GPU annotations shared by every GTRK header.
///
/// This is the single place GPUCommonDef.h is pulled in from, so that the
/// host/device annotations arrive identically everywhere.
///
/// Nothing in this file may depend on ROOT or on anything that cannot be
/// compiled by a device compiler. The STL is limited to <cstdint> and <limits>,
/// both of which are pure compile-time facilities here.

#ifndef GTRK_BASE_DEFINITIONS_H_
#define GTRK_BASE_DEFINITIONS_H_

#include "GPUCommonDef.h"

#include <cstdint>
#include <limits>

// The function qualifiers. GTRK writes them in O2's vocabulary (GPUCommonDef.h),
// so that a header reads the same here and in O2:
//
//   GPUhdi() GPUAlwaysInline()
//              host and device, always inlined -- every helper a kernel can
//              reach: the cuts, terms, fits, steps, accessors and math a stage
//              calls per candidate. Forced, so that a caller growing past
//              GCC's inlining limits cannot push one out of line and add a
//              call per candidate.
//   GPUhdnii() host and device, `inline` on the host and never inlined on the
//              device -- the heavy leaf every stage carries a state through,
//              which the device otherwise copies into every search and window
//              of every variant.
//   GPUhdi()   host and device, `inline`, the compiler decides -- the stage
//              bodies only: a stage's per-element `operator()`, the walks and
//              drivers around the candidates, and the few large computations
//              a compiler keeps out of line; test/testQualifiers.cxx lists them.
//   GPUdi()    device only (`inline` on the host) -- what calls one of O2's
//              device-only functions, `CAMath`'s atomics or the material
//              budget's `getMatBudget`; the device
//              compiler refuses it as `GPUhdi()`, so this is never a choice.
//   GPUg()     a kernel.
//
// Host-only functions carry plain C++ (`inline`, `constexpr`) and no qualifier.
// An attribute comes first, the qualifier next, the other specifiers after it:
// `[[nodiscard]] GPUhdi() GPUAlwaysInline() constexpr bool f(...)`.

namespace o2::gtrk
{

/// Global index of a sensor (chip) within a Layout.
///
/// Dense and stable: sensor IDs run 0 ... nSensors-1 in a fixed geometric
/// order, which is what lets every downstream stage index chip-keyed arrays
/// directly instead of going through a map. It is the detector's own chip
/// index (O2's chip ID for the ITS, the MFT and the TRK), so a cluster's
/// sensorID needs no decoding.
using SensorID = std::uint32_t;

/// Index of a layer within a Layout. Layers are few; the narrow type keeps the
/// per-cluster records small, which matters once they are uploaded per ROF.
using LayerID = std::uint16_t;

/// Sentinel for "no sensor". Chosen as the maximum representable value rather
/// than -1 so that the unsigned comparison `id < nSensors` rejects it without a
/// separate check -- the form that survives being written as a device-side
/// bounds test.
inline constexpr SensorID kInvalidSensorID = std::numeric_limits<SensorID>::max();

/// Sentinel for "no layer", same reasoning as kInvalidSensorID.
inline constexpr LayerID kInvalidLayerID = std::numeric_limits<LayerID>::max();

/// `a == b` for floating point, for values that are equal by construction
/// rather than by arithmetic: a switched-off field, an absent neighbour's
/// zero, a constant a test pins down. The same result as `==` for every input,
/// NaN and signed zeros included, written without `==` because GCC's
/// -Wfloat-equal, unlike Clang's, also flags comparisons with a literal.
template <typename T>
GPUhdi() constexpr bool exactlyEqual(T a, T b)
{
  return a <= b && a >= b;
}

} // namespace o2::gtrk

#endif // GTRK_BASE_DEFINITIONS_H_
