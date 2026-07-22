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

} // namespace o2::itsmft::tracking::barrel

#endif // ALICEO2_ITSMFT_TRACKING_BARRELSURFACESTATEOPERATIONS_H_
