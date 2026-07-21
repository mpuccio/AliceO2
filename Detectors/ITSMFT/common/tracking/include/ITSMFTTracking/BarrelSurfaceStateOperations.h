// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_BARRELSURFACESTATEOPERATIONS_H_
#define ALICEO2_ITSMFT_TRACKING_BARRELSURFACESTATEOPERATIONS_H_

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

} // namespace o2::itsmft::tracking::barrel

#endif // ALICEO2_ITSMFT_TRACKING_BARRELSURFACESTATEOPERATIONS_H_
