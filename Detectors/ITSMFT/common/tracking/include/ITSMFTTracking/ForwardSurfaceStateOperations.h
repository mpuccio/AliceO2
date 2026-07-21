// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_FORWARDSURFACESTATEOPERATIONS_H_
#define ALICEO2_ITSMFT_TRACKING_FORWARDSURFACESTATEOPERATIONS_H_

#include <cstdint>

#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"

namespace o2::itsmft::tracking::forward
{

// Values intentionally match the four legacy MFT propagation choices. The
// template argument binds the model outside a tracking hot loop.
enum class PropagationModel : uint8_t {
  Linear,
  Quadratic,
  Helix,
  Optimized // Helix parameters with the quadratic covariance Jacobian.
};

// All operations below are host-only and float-native. Parameters,
// covariance, local arithmetic, and committed results remain float; no legacy
// forward track object is constructed. Mutating operations use scratch then
// commit and therefore leave state byte-for-byte unchanged on failure.
template <PropagationModel Model>
bool propagate(SurfaceKinematicState& state, float targetZ, float bz, OperationFailureReason& reason) noexcept;

template <>
bool propagate<PropagationModel::Linear>(SurfaceKinematicState&, float, float, OperationFailureReason&) noexcept;
template <>
bool propagate<PropagationModel::Quadratic>(SurfaceKinematicState&, float, float, OperationFailureReason&) noexcept;
template <>
bool propagate<PropagationModel::Helix>(SurfaceKinematicState&, float, float, OperationFailureReason&) noexcept;
template <>
bool propagate<PropagationModel::Optimized>(SurfaceKinematicState&, float, float, OperationFailureReason&) noexcept;

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

// Preserves the legacy forward, mass/PID/absCharge-agnostic MCS behavior.
bool correctForMaterial(SurfaceKinematicState& state, float xOverX0, OperationFailureReason& reason) noexcept;

} // namespace o2::itsmft::tracking::forward

#endif // ALICEO2_ITSMFT_TRACKING_FORWARDSURFACESTATEOPERATIONS_H_
