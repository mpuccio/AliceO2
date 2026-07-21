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

namespace o2::itsmft::tracking
{

// Values intentionally match the four legacy MFT propagation choices. The
// template argument binds the model outside a tracking hot loop.
enum class ForwardPropagationModel : uint8_t {
  Linear,
  Quadratic,
  Helix,
  Optimized // Helix parameters with the quadratic covariance Jacobian.
};

// This is the subset needed by the isolated forward primitive slice. Values
// may be appended when the complete Stage-B operation contract is accepted.
enum class OperationFailureReason : uint8_t {
  SourceFamilyMismatch,
  NonFiniteInput,
  NonFiniteOutput,
  InvalidCovariance,
  UnreachableTarget,
  PropagationFailure,
  MaterialFailure,
  PredictedChi2Failure,
  UpdateFailure
};

// All operations below are host-only and float-native. Parameters,
// covariance, local arithmetic, and committed results remain float; no legacy
// forward track object is constructed. Mutating operations use scratch then
// commit and therefore leave state byte-for-byte unchanged on failure.
template <ForwardPropagationModel Model>
bool propagateToDisk(SurfaceKinematicState& state, float targetZ, float bz, OperationFailureReason& reason) noexcept;

template <>
bool propagateToDisk<ForwardPropagationModel::Linear>(SurfaceKinematicState&, float, float, OperationFailureReason&) noexcept;
template <>
bool propagateToDisk<ForwardPropagationModel::Quadratic>(SurfaceKinematicState&, float, float, OperationFailureReason&) noexcept;
template <>
bool propagateToDisk<ForwardPropagationModel::Helix>(SurfaceKinematicState&, float, float, OperationFailureReason&) noexcept;
template <>
bool propagateToDisk<ForwardPropagationModel::Optimized>(SurfaceKinematicState&, float, float, OperationFailureReason&) noexcept;

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

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_FORWARDSURFACESTATEOPERATIONS_H_
