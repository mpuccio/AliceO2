// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_FORWARDSURFACESTATEOPERATIONS_H_
#define ALICEO2_ITSMFT_TRACKING_FORWARDSURFACESTATEOPERATIONS_H_

#include <cstdint>

#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"

namespace o2::itsmft::tracking::forward
{

// Disk measurements use global x/y as u/v, including the full uu/uv/vv
// covariance. chi2 is unchanged on failure.
bool predictedChi2(const SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
                   OperationFailureReason& reason) noexcept;

// Returns the chi2 increment through chi2 and applies the 2D Kalman update.
// Both outputs are unchanged on failure.
bool update(SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
            OperationFailureReason& reason) noexcept;

// Variance core used by the forward Highland operation.
constexpr float highlandTheta2(float inverseMomentum, float xOverX0) noexcept
{
  const float theta = 0.0136f * inverseMomentum;
  return theta * theta * xOverX0;
}

// Forward Highland MCS correction; this overload is mass/PID/absCharge
// agnostic and remains separate from the scalar material kernel below.
bool correctForMaterial(SurfaceKinematicState& state, float xOverX0, OperationFailureReason& reason) noexcept;

// PID/absCharge-aware material correction through the shared scalar kernel;
// the projected covariance is committed transactionally. Diagnostics follow
// MaterialOperationResult and the 92-byte state is unchanged unless result.ok().
material::MaterialOperationResult correctForMaterial(SurfaceKinematicState& state, material::IntegratedMaterialBudget materialBudget,
                                                     material::MaterialTraversalDirection direction) noexcept;

// Same-family compatibility chi2 between two forward states expressed at the
// same reference Z. Residuals are the direct (unwrapped) differences of
// (X, Y, Phi, Tanl, InvQPt), preserving the legacy same-family behavior of
// detail::mftFwdStateChi2. Neither input is mutated and chi2 is unchanged on
// failure.
bool stateChi2(const SurfaceKinematicState& reference, const SurfaceKinematicState& candidate, float& chi2,
               OperationFailureReason& reason) noexcept;

#ifndef GPUCA_GPUCODE

// Builds an outer-anchored forward/disk three-hit seed from ordered
// {inner, middle, outer} measurements. The outer frame supplies z and
// covariance; trackletMinPt uses the established `(trackletMinPt > 0.f) ?
// 1.f/trackletMinPt : 0.f` fallback. Finite-input checks precede the strict
// boundary requiring inner z to exceed outer z by 1e-6f and all relevant
// separations to exceed 1e-6f, followed by the NonFiniteOutput check.
// Construction is scratch-then-commit and leaves outState unchanged on failure.
bool buildSeed(const SurfaceMeasurement& measurementInner, const SurfaceMeasurement& measurementMiddle,
               const SurfaceMeasurement& measurementOuter, float bz, float trackletMinPt,
               uint8_t absCharge, o2::track::PID pid,
               SurfaceKinematicState& outState, OperationFailureReason& reason) noexcept;

// Sets a forward reference's X/Y parameters from measurement u/v (q=z),
// without changing referenceCoordinate, alpha, Phi, Tanl, or InvQPt.
// Requires a forward reference and finite u/v; failure leaves it unchanged.
bool shiftReferenceToMeasurement(SurfaceLinearizationReference& linRef, const SurfaceMeasurement& measurement,
                                 OperationFailureReason& reason) noexcept;

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking::forward

#endif // ALICEO2_ITSMFT_TRACKING_FORWARDSURFACESTATEOPERATIONS_H_
