// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_DETAIL_SURFACESTATEOPERATIONS_H_
#define ALICEO2_ITSMFT_TRACKING_DETAIL_SURFACESTATEOPERATIONS_H_

#include "ITSMFTTracking/GlobalMeasurement.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"

// Coordinate-family leaves used only by Propagator, candidate construction,
// and their numerical tests. Production callers use Propagator's
// descriptor/state-driven API rather than selecting a family themselves.
namespace o2::itsmft::tracking::detail
{
namespace barrel
{
bool rotate(SurfaceKinematicState& state, float targetAlpha, OperationFailureReason& reason) noexcept;
bool propagate(SurfaceKinematicState& state, float targetX, float bz, OperationFailureReason& reason) noexcept;
bool predictedChi2(const SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
                   OperationFailureReason& reason) noexcept;
bool update(SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
            OperationFailureReason& reason) noexcept;
material::MaterialOperationResult correctForMaterial(SurfaceKinematicState& state, material::IntegratedMaterialBudget materialBudget,
                                                     material::MaterialTraversalDirection direction) noexcept;
bool stateChi2(const SurfaceKinematicState& reference, const SurfaceKinematicState& candidate, float& chi2,
               OperationFailureReason& reason) noexcept;

#ifndef GPUCA_GPUCODE
bool buildSeed(const GlobalPoint3F& globalInner, const GlobalPoint3F& globalMiddle,
               const SurfaceMeasurement& measurementOuter, float bz,
               uint8_t absCharge, o2::track::PID pid,
               SurfaceKinematicState& outState, OperationFailureReason& reason) noexcept;
bool rotate(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef, float targetAlpha, float bz,
            OperationFailureReason& reason) noexcept;
bool propagate(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef, float targetX, float bz,
               OperationFailureReason& reason) noexcept;
bool shiftReferenceToMeasurement(SurfaceLinearizationReference& linRef, const SurfaceMeasurement& measurement,
                                 OperationFailureReason& reason) noexcept;
#endif
} // namespace barrel

namespace forward
{
bool propagate(SurfaceKinematicState& state, float targetZ, float bz, OperationFailureReason& reason) noexcept;
bool propagate(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef,
               float targetZ, float bz, OperationFailureReason& reason) noexcept;
bool predictedChi2(const SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
                   OperationFailureReason& reason) noexcept;
bool update(SurfaceKinematicState& state, const SurfaceMeasurement& measurement, float& chi2,
            OperationFailureReason& reason) noexcept;
constexpr float highlandTheta2(float inverseMomentum, float xOverX0) noexcept
{
  const float theta = 0.0136f * inverseMomentum;
  return theta * theta * xOverX0;
}
bool correctForMaterial(SurfaceKinematicState& state, float xOverX0, OperationFailureReason& reason) noexcept;
material::MaterialOperationResult correctForMaterial(SurfaceKinematicState& state, material::IntegratedMaterialBudget materialBudget,
                                                     material::MaterialTraversalDirection direction) noexcept;
bool stateChi2(const SurfaceKinematicState& reference, const SurfaceKinematicState& candidate, float& chi2,
               OperationFailureReason& reason) noexcept;

#ifndef GPUCA_GPUCODE
bool buildSeed(const SurfaceMeasurement& measurementInner, const SurfaceMeasurement& measurementMiddle,
               const SurfaceMeasurement& measurementOuter, float bz, float trackletMinPt,
               uint8_t absCharge, o2::track::PID pid,
               SurfaceKinematicState& outState, OperationFailureReason& reason) noexcept;
bool shiftReferenceToMeasurement(SurfaceLinearizationReference& linRef, const SurfaceMeasurement& measurement,
                                 OperationFailureReason& reason) noexcept;
#endif
} // namespace forward
} // namespace o2::itsmft::tracking::detail

#endif
