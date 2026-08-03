// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACELINEARIZATIONREFERENCE_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACELINEARIZATIONREFERENCE_H_

#include <cstddef>
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/StateFamily.h"

namespace o2::itsmft::tracking
{

// Transient, covariance-free trajectory used only as a Kalman
// linearization point (legacy name "linRef") during an ITS-shaped refit
// leg -- the parameter-only companion to SurfaceKinematicState that
// rotate/propagate/reference-shift operations transport alongside the
// fitted state. Searched first: no existing public O2 parameter-only type
// satisfies this contract without pulling in a legacy dependency --
// o2::track::TrackParametrization<float> (TrackPar), the direct base class
// of the excluded TrackParCovF, has no StateFamily tag (so it cannot
// represent both the Barrel and Forward parameter conventions in one
// dependency-light POD without misusing its fields), and pulls in the wide
// legacy header footprint (BaseCluster.h, TrackUtils.h, PID.h, ...) this
// slice is required to avoid.
//
// Deliberately excludes covariance, chi2, cluster bookkeeping, and
// PID/absCharge: a linearization reference always accompanies exactly one
// SurfaceKinematicState (the fitted state being linearized), which already
// carries the particle hypothesis. Operations that need absCharge/pid
// (e.g. curvature) read them from that paired state, never from a
// duplicated copy here. This also means a SurfaceLinearizationReference
// has no meaning on its own; it is only ever valid in the context of the
// SurfaceKinematicState it was constructed from or is paired with.
//
// Barrel: (Y, Z, Snp, Tgl, Q2Pt), referenceCoordinate is local X, alpha is
// frame angle -- identical parameter/frame convention to
// SurfaceKinematicState's own Barrel interpretation.
// Forward: (X, Y, Phi, Tanl, InvQPt), referenceCoordinate is global Z,
// alpha is unused (zero) -- identical to SurfaceKinematicState's Forward
// interpretation.
//
// This is a current in-memory layout commitment (standard-layout,
// trivially-copyable, explicit size/alignment/offsets below), not yet a
// durable serialized or device ABI promise -- matching the same
// current-layout-only caveat MaterialPhysics.h already documents for its
// own PODs.
struct SurfaceLinearizationReference {
  float parameters[5]{};
  float referenceCoordinate{0.f};
  float alpha{0.f};
  StateFamily family{StateFamily::Invalid};

  GPUhdi() constexpr bool hasRecognizedFamily() const noexcept { return family == StateFamily::Barrel || family == StateFamily::Forward; }
};

static_assert(std::is_standard_layout_v<SurfaceLinearizationReference>);
static_assert(std::is_trivially_copyable_v<SurfaceLinearizationReference>);
static_assert(sizeof(SurfaceLinearizationReference) == 32);
static_assert(alignof(SurfaceLinearizationReference) == 4);
static_assert(offsetof(SurfaceLinearizationReference, parameters) == 0);
static_assert(offsetof(SurfaceLinearizationReference, referenceCoordinate) == 20);
static_assert(offsetof(SurfaceLinearizationReference, alpha) == 24);
static_assert(offsetof(SurfaceLinearizationReference, family) == 28);

#ifndef GPUCA_GPUCODE

// Host-only for this slice (no GPU/device-readiness claim, matching
// buildSeed/buildCellSeed/attachHit's own host-only contract): validating
// factory/transactional-mutation primitive. Copies only
// parameters/referenceCoordinate/alpha/family out of `state` into `out`.
// This same operation serves both roles required by the refit orchestration
// this slice prepares for: constructing the initial linearization reference
// from a freshly built seed, and re-anchoring it at a later fit leg's
// starting state (the legacy `linRef = track.paramOut;` idiom). Fails if
// `state.hasRecognizedFamily()` is false or any copied field is
// non-finite; `out` is left completely unchanged on failure.
bool makeLinearizationReference(const SurfaceKinematicState& state, SurfaceLinearizationReference& out) noexcept;

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_SURFACELINEARIZATIONREFERENCE_H_
