// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_BARRELSURFACESTATEOPERATIONS_H_
#define ALICEO2_ITSMFT_TRACKING_BARRELSURFACESTATEOPERATIONS_H_

#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/GlobalMeasurement.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"

namespace o2::itsmft::tracking::barrel
{

// Rotates the local barrel frame to targetAlpha, canonicalizing alpha to
// [-pi, pi] and transforming referenceCoordinate and the local-X plane.
// Fails if the rotated state has no forward-going local direction.
bool rotate(SurfaceKinematicState& state, float targetAlpha, OperationFailureReason& reason) noexcept;

// Propagates in the current local barrel frame to targetX in uniform Bz.
// On success referenceCoordinate equals targetX and alpha is unchanged.
bool propagate(SurfaceKinematicState& state, float targetX, float bz, OperationFailureReason& reason) noexcept;

// Uses measurement frame.u/frame.v as local Y/Z. chi2 is unchanged on failure.
bool predictedChi2(const SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
                   OperationFailureReason& reason) noexcept;

// Applies the local Y/Z two-dimensional Kalman update and returns its chi2
// increment. State and chi2 remain unchanged on failure.
bool update(SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
            OperationFailureReason& reason) noexcept;

// Applies PID/absCharge-aware material correction through the shared scalar
// kernel and commits projected covariance transactionally. Diagnostics follow
// MaterialOperationResult: preflight/projection failures zero all fields
// except momentumBeforeGeV and failure; scalar results pass through unchanged.
// The state is unchanged unless result.ok().
material::MaterialOperationResult correctForMaterial(SurfaceKinematicState& state, material::IntegratedMaterialBudget materialBudget,
                                                     material::MaterialTraversalDirection direction) noexcept;

// Computes same-cylinder compatibility chi2 for states at the same reference X
// and alpha. Residuals are direct (unwrapped) differences of (Y, Z, Snp, Tgl,
// Q2Pt). Inputs are not mutated; chi2 is unchanged on failure.
bool stateChi2(const SurfaceKinematicState& reference, const SurfaceKinematicState& candidate, float& chi2,
               OperationFailureReason& reason) noexcept;

#ifndef GPUCA_GPUCODE

// Builds an outer-anchored cylindrical three-hit seed. Inner and middle
// positions are global; the outer native frame supplies anchor, local
// position, and covariance. Copies caller-fixed absCharge/pid on success and
// sets flags to zero. Degenerate formula fallbacks are retained. Construction
// is scratch-then-commit, leaving outState byte-for-byte unchanged on
// NonFiniteInput/Output.
bool buildSeed(const GlobalPoint3F& globalInner, const GlobalPoint3F& globalMiddle,
               const SurfaceMeasurement& measurementOuter, float bz,
               uint8_t absCharge, o2::track::PID pid,
               SurfaceKinematicState& outState, OperationFailureReason& reason) noexcept;

// Rotates a fitted barrel state around a paired linearization reference.
// Pairing requires exact kind/referenceCoordinate/alpha equality, while
// state and reference parameters may differ. Both objects are
// scratch-then-commit and byte-for-byte unchanged on failure.
bool rotate(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef, float targetAlpha, float bz,
            OperationFailureReason& reason) noexcept;

// Propagates a fitted barrel state using a paired linearization reference;
// pairing follows rotate. A zero step updates both reference coordinates to
// targetX. UnreachableTarget/PropagationFailure leave both objects unchanged.
bool propagate(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef, float targetX, float bz,
               OperationFailureReason& reason) noexcept;

// Sets only the barrel linearization reference Y/Z parameters from local u/v.
// The frame and other parameters are untouched. Requires a barrel reference
// and finite u/v; failure leaves the reference unchanged.
bool shiftReferenceToMeasurement(SurfaceLinearizationReference& linRef, const SurfaceMeasurement& measurement,
                                 OperationFailureReason& reason) noexcept;

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking::barrel

#endif // ALICEO2_ITSMFT_TRACKING_BARRELSURFACESTATEOPERATIONS_H_
