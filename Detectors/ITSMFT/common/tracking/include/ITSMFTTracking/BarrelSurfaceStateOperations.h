// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_BARRELSURFACESTATEOPERATIONS_H_
#define ALICEO2_ITSMFT_TRACKING_BARRELSURFACESTATEOPERATIONS_H_

#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"

namespace o2::itsmft::tracking::barrel
{

// Rotates the local barrel frame to targetAlpha. On success alpha is
// canonicalized to [-pi, pi], while referenceCoordinate and the local-X plane
// are transformed with the frame. The operation fails if the rotated state
// would not have a forward-going local direction.
bool rotate(SurfaceKinematicState& state, float targetAlpha, OperationFailureReason& reason) noexcept;

// Propagates in the current local barrel frame to targetX in a uniform Bz.
// On success referenceCoordinate is exactly targetX and alpha is unchanged.
bool propagate(SurfaceKinematicState& state, float targetX, float bz, OperationFailureReason& reason) noexcept;

// Barrel measurements use frame.u/frame.v as local Y/Z, respectively. chi2
// is unchanged on failure.
bool predictedChi2(const SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
                   OperationFailureReason& reason) noexcept;

// Applies the local Y/Z two-dimensional Kalman update and returns its chi2
// increment. State and chi2 are unchanged on failure.
bool update(SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
            OperationFailureReason& reason) noexcept;

// PID/absCharge-aware material correction using the shared scalar kernel
// (MaterialPhysics.h): derives physical momentum from the pre-material
// state, evaluates calculateMaterialPhysics() with the state's own PID and
// absCharge, and projects the resulting energy-loss/scattering/straggling
// contributions onto the packed covariance.
//
// Diagnostics contract (result is always the scalar MaterialOperationResult,
// see MaterialPhysics.h for field semantics):
//  - Preflight failure (state validation or physical-momentum derivation
//    fails before the scalar kernel runs): momentumBeforeGeV, every other
//    numeric diagnostic, and flags are all zero/None; only failure is set.
//  - Scalar-kernel failure: the scalar kernel's own result is returned
//    unchanged (momentumBeforeGeV echoes the derived input momentum; all
//    other diagnostics are that kernel's deterministic zero/None values).
//  - Projection failure (the scalar kernel succeeds but the projected
//    scratch state fails post-projection validation): momentumBeforeGeV is
//    preserved from the scalar result; every other numeric diagnostic and
//    flags are zero/None; only failure is set.
//  - Success: the scalar kernel's successful result is returned unchanged.
// In every case above, the complete 92-byte state is left byte-for-byte
// unchanged unless the operation returns a successful result (result.ok()).
material::MaterialOperationResult correctForMaterial(SurfaceKinematicState& state, material::IntegratedMaterialBudget materialBudget,
                                                     material::MaterialTraversalDirection direction) noexcept;

// Same-family compatibility chi2 between two barrel states expressed at the
// same reference X and alpha. Residuals are the direct (unwrapped)
// differences of (Y, Z, Snp, Tgl, Q2Pt). Neither input is mutated and chi2 is
// unchanged on failure.
bool stateChi2(const SurfaceKinematicState& reference, const SurfaceKinematicState& candidate, float& chi2,
               OperationFailureReason& reason) noexcept;

#ifndef GPUCA_GPUCODE

// Slice A (Stage-B design report Sec 8/11), migrated in the Stage-B
// normalized-measurement slice to consume SurfaceMeasurement directly:
// builds the initial *outer*-anchored SurfaceKinematicState seed for a
// cylindrical three-hit Cell candidate, transcribing the exact
// o2::its::track::buildTrackSeed initialization (ITStracking/TrackHelpers.h)
// directly onto SurfaceKinematicState. Production must never call
// buildTrackSeed or construct the legacy barrel track-parametrization-with-
// error type it returns; this operation reproduces the identical parameter
// meaning, covariance initialization, alpha/reference-coordinate convention
// and field/sign convention using only float arithmetic and the packed
// SurfaceKinematicState layout.
//
// Input order matches this operation family's existing {inner, middle,
// outer} contract, never inferred from numeric layer/SurfaceId/radius/z:
// `measurementInner`/`measurementMiddle` are read only for their
// `measurement.global` (x, y, z) position; `measurementOuter` additionally
// supplies the seed's reference frame (`measurement.frame.q`/
// `measurement.frame.frameAngle`) and its own measured (u, v) position
// (`measurement.frame.u`/`measurement.frame.v`) and covariance
// (`measurement.covariance`). The outer measurement's global position is not
// a parameter -- buildTrackSeed never reads it either, only the outer hit's
// already-frame-expressed tracking-frame fields.
//
// Output anchor/reference frame: referenceCoordinate == measurementOuter.
// frame.q and alpha == measurementOuter.frame.frameAngle -- the *outer*
// measurement's own tracking frame. This is deliberately not the Cell's
// eventual inner-anchored frame (anchor contract, design report Sec 5): this
// operation produces only the initial outer seed that the existing
// outer->middle->inner buildCellSeed sequence subsequently
// rotates/propagates/updates inward, hit by hit, to reach the Cell's actual
// innermost-surface anchor. Calling buildSeed alone does not produce a
// complete Cell state.
//
// Compatibility hypothesis (absCharge, pid): supplied by the caller and
// never defaulted here. SurfaceKinematicState's own absCharge=0 default
// would silently select neutral-material behavior in later material
// operations, so this operation always sets state.absCharge/state.pid from
// the caller-supplied values on success. Unchecked caller precondition
// (matching every other operation in this file family): absCharge/pid are
// caller-fixed constants, not per-candidate data, and are not independently
// validated here. flags is always set to 0.
//
// Failure vocabulary: OperationFailureReason::NonFiniteInput if any raw
// input is not finite; OperationFailureReason::NonFiniteOutput if the fully
// constructed candidate state is not finite (the retained formula has no
// other explicit rejection -- its own Almost0/VeryBig degenerate-geometry
// fallbacks are preserved verbatim and only surface here if they ultimately
// produce a non-finite parameter or covariance entry).
//
// Transactional: constructed entirely in local scratch; outState is
// committed only on complete success. On any failure outState is left
// exactly as passed in, byte-for-byte.
bool buildSeed(const SurfaceMeasurement& measurementInner, const SurfaceMeasurement& measurementMiddle,
               const SurfaceMeasurement& measurementOuter, float bz,
               uint8_t absCharge, o2::track::PID pid,
               SurfaceKinematicState& outState, OperationFailureReason& reason) noexcept;

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking::barrel

#endif // ALICEO2_ITSMFT_TRACKING_BARRELSURFACESTATEOPERATIONS_H_
