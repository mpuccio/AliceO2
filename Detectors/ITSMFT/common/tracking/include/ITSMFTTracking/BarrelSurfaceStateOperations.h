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
#include "ITSMFTTracking/StateFamily.h"

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

// PID/absCharge-aware material correction through the shared scalar kernel;
// the projected covariance is committed transactionally. Diagnostics follow
// MaterialOperationResult: preflight/projection failures zero all fields
// except the prescribed momentumBeforeGeV and failure; scalar results pass
// through unchanged. The 92-byte state is unchanged unless result.ok().
material::MaterialOperationResult correctForMaterial(SurfaceKinematicState& state, material::IntegratedMaterialBudget materialBudget,
                                                     material::MaterialTraversalDirection direction) noexcept;

// Same-family compatibility chi2 between two barrel states expressed at the
// same reference X and alpha. Residuals are the direct (unwrapped)
// differences of (Y, Z, Snp, Tgl, Q2Pt). Neither input is mutated and chi2 is
// unchanged on failure.
bool stateChi2(const SurfaceKinematicState& reference, const SurfaceKinematicState& candidate, float& chi2,
               OperationFailureReason& reason) noexcept;

#ifndef GPUCA_GPUCODE

// Builds an outer-anchored cylindrical three-hit seed. Inputs are ordered
// {inner, middle, outer}; the outer frame supplies referenceCoordinate,
// alpha, local position, and covariance. absCharge/pid are caller-fixed and
// copied on success; flags is zero. Degenerate formula fallbacks are retained.
// Construction is scratch-then-commit: outState is byte-for-byte unchanged
// on NonFiniteInput or NonFiniteOutput.
bool buildSeed(const SurfaceMeasurement& measurementInner, const SurfaceMeasurement& measurementMiddle,
               const SurfaceMeasurement& measurementOuter, float bz,
               uint8_t absCharge, o2::track::PID pid,
               SurfaceKinematicState& outState, OperationFailureReason& reason) noexcept;

// Rotates a fitted barrel state around a paired linearization reference.
// Pairing requires exact family/referenceCoordinate/alpha equality; state and
// reference parameters may differ. Preconditions retain the established
// |snp| and post-rotation checks. Both objects are scratch-then-commit and
// unchanged byte-for-byte on failure.
bool rotate(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef, float targetAlpha, float bz,
            OperationFailureReason& reason) noexcept;

// Propagates a fitted barrel state using a paired linearization reference.
// The exact pairing contract is that of rotate. A zero step updates both
// reference coordinates to targetX without changing other parameters.
// UnreachableTarget/PropagationFailure leave both objects unchanged.
bool propagate(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef, float targetX, float bz,
               OperationFailureReason& reason) noexcept;

// Sets only the barrel linearization reference's Y/Z parameters from the
// measurement's local u/v coordinates. Reference frame and other parameters
// are untouched. Requires a barrel reference and finite u/v; failure leaves
// the reference unchanged.
bool shiftReferenceToMeasurement(SurfaceLinearizationReference& linRef, const SurfaceMeasurement& measurement,
                                 OperationFailureReason& reason) noexcept;

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking::barrel

#endif // ALICEO2_ITSMFT_TRACKING_BARRELSURFACESTATEOPERATIONS_H_
