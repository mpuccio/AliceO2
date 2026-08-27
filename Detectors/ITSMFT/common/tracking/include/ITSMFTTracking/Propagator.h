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

#ifndef ALICEO2_ITSMFT_TRACKING_PROPAGATOR_H_
#define ALICEO2_ITSMFT_TRACKING_PROPAGATOR_H_

#include "GPUCommonDef.h"

#ifndef GPUCA_GPUCODE

#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceTrackState.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"

// Descriptor-driven propagation using the material and kind resolved from
// SurfaceDescriptor and SurfaceCatalogView.
namespace o2::itsmft::tracking
{

class Propagator
{
 public:
  // Advance a state in its current surface convention to measurement, apply
  // the supplied nominal material, validate the predicted residual, and
  // update it. Candidate construction uses this state-driven operation after
  // building a family-local three-hit seed; it does not select a coordinate
  // family itself. State and chi2 are committed only after complete success.
  static bool attachMeasurement(SurfaceTrackState& state, const SurfaceMeasurement& measurement,
                                NominalSurfaceMaterial material, float bz,
                                material::MaterialTraversalDirection direction,
                                bool chi2GateEnabled, float maxChi2, float& chi2,
                                OperationFailureReason& reason) noexcept;

  // Compatibility chi2 for two states in the same surface convention. The
  // coordinate convention is selected from the states, never by the caller.
  static bool stateChi2(const SurfaceTrackState& reference, const SurfaceTrackState& candidate,
                        float& chi2, OperationFailureReason& reason) noexcept;

  // Propagate in the state’s current surface convention to its target
  // reference coordinate. Disk transport uses helix propagation for
  // |bz| > 0.01f and linear transport otherwise. Both objects are unchanged
  // on failure when a linearization reference is supplied.
  static bool propagateToReference(SurfaceTrackState& state, float targetReferenceCoordinate, float bz,
                                   OperationFailureReason& reason) noexcept;
  static bool propagateToReference(SurfaceTrackState& state, SurfaceTrackParameters& linRef,
                                   float targetReferenceCoordinate, float bz,
                                   OperationFailureReason& reason) noexcept;

  // Re-express state in targetKind's parameter/reference convention at the
  // current reference surface. A matching kind is a no-op.
  //
  // This is a first-order Jacobian reparametrization between Barrel's
  // radius/angle reference and Forward's z reference. The new reference's
  // fitted variance is discarded; the newly freed coordinate receives a
  // loose, uninformative ceiling variance.
  //
  // Preserves absCharge, PID, and all fields outside the parameter
  // convention. On failure, state is unchanged.
  static bool convertKind(SurfaceTrackState& state, SurfaceKind targetKind,
                          OperationFailureReason& reason) noexcept;

  // Propagate to a measurement, converting the state to the target surface
  // kind when needed, then applying material, the chi2 gate, and the update.
  // State, reference, and chi2 are committed only after complete success.
  //
  // The incoming chi2 must be finite and non-negative. maxChi2 is validated
  // the same way when the gate is enabled.
  static bool propagateToMeasurement(SurfaceTrackState& state, SurfaceTrackParameters& linRef,
                                     const SurfaceDescriptor& targetSurface, const SurfaceMeasurement& targetMeasurement,
                                     float bz, material::MaterialTraversalDirection direction,
                                     bool chi2GateEnabled, float maxChi2, float& chi2,
                                     bool shiftReferenceToMeasurement, OperationFailureReason& reason) noexcept;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_PROPAGATOR_H_ */
