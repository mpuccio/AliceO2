// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file BoundTrackState.h
/// \brief Five parameters and their covariance, bound to one sensor plane.
///
/// This is data and nothing else.  It has no opinion about how it moves; that
/// belongs to the propagator, a separate type.  Free functions propagate and update it; its
/// members are accessors and `decorrelate`, so the record can be copied into
/// device memory.

#ifndef GTRK_TRACKING_BOUNDTRACKSTATE_H_
#define GTRK_TRACKING_BOUNDTRACKSTATE_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Tracking/Definitions.h"

namespace o2::gtrk
{

/// A track state on a sensor plane: the five geometric parameters. A track's
/// time is its rows' box, kept beside it.
///
/// The plane is the sensor `surface` names; on a curved sensor it is the
/// sensor's tangent plane at the sensor-local angle `surfacePhi`
/// (TrackingGeometry::planeOf). That angle
/// names the frame, it is not a parameter: an update leaves it alone, and it
/// has no covariance. On a flat sensor it is 0.
struct BoundTrackState {
  Scalar par[kNBoundParams]{};         ///< loc0, loc1, phi, theta, q/p
  Scalar cov[kBoundCovSize]{};         ///< packed lower triangle; see covIndex()
  SensorID surface = kInvalidSensorID; ///< the sensor these parameters are bound to
  Scalar surfacePhi = 0.F;             ///< a curved sensor's tangent plane: its angle in the sensor [rad]

  GPUhdi() GPUAlwaysInline() Scalar loc0() const { return par[kLoc0]; }
  GPUhdi() GPUAlwaysInline() Scalar loc1() const { return par[kLoc1]; }
  GPUhdi() GPUAlwaysInline() Scalar phi() const { return par[kPhi]; }
  GPUhdi() GPUAlwaysInline() Scalar theta() const { return par[kTheta]; }
  GPUhdi() GPUAlwaysInline() Scalar qOverP() const { return par[kQOverP]; }

  /// Covariance element (i, j).  Symmetric, so the caller never has to order
  /// the indices.
  GPUhdi() GPUAlwaysInline() Scalar covAt(int i, int j) const { return cov[covIndex(i, j)]; }
  GPUhdi() GPUAlwaysInline() Scalar& covAt(int i, int j) { return cov[covIndex(i, j)]; }

  GPUhdi() GPUAlwaysInline() Scalar variance(int i) const { return cov[varIndex(i)]; }
  GPUhdi() GPUAlwaysInline() Scalar& variance(int i) { return cov[varIndex(i)]; }

  /// Drop the covariances of parameters \a a and \a b with every other
  /// parameter; their variances and the rest stay.
  GPUhdi() GPUAlwaysInline() void decorrelate(int a, int b)
  {
    for (int i = 0; i < kNBoundParams; ++i) {
      if (i != a) {
        covAt(i, a) = 0.F;
      }
      if (i != b) {
        covAt(i, b) = 0.F;
      }
    }
  }
};

// Its bytes are its value: a record copies and compares it as bytes
// (Study/StageComparison.h).
static_assert(sizeof(BoundTrackState) == (sizeof(Scalar) * (kNBoundParams + kBoundCovSize + 1)) + sizeof(SensorID), "a state has no padding");

/// A state's covariance as P = U D U^T: U unit upper triangular, of which the
/// strict upper triangle is stored (MathUtils.h, strictUpperIndex), and D
/// diagonal. A filter pass carries its covariance in this form from one
/// update to the next: the step's prediction is formed on the factors
/// (transportFactors) and the update updates them (updateWithMeasurement), so
/// no predicted P is ever rounded to float. Rounding it loses the variance a
/// precise measurement leaves once a diffuse direction is carried over a long
/// lever arm.
struct CovarianceFactors {
  Scalar upper[kBoundUpperSize]{};
  Scalar diagonal[kNBoundParams]{};
};

/// Whether every factor in D is positive; written so that NaN fails, without
/// an algorithm device code lacks.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool positiveDiagonal(const CovarianceFactors& factors)
{
  bool positive = true;
  for (const Scalar d : factors.diagonal) {
    positive = positive && d > 0.F;
  }
  return positive;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_BOUNDTRACKSTATE_H_
