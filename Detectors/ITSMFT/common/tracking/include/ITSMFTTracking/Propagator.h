// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_PROPAGATOR_H_
#define ALICEO2_ITSMFT_TRACKING_PROPAGATOR_H_

#include "GPUCommonDef.h"

#ifndef GPUCA_GPUCODE

#include <gsl/span>

#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/StateFamily.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"

// Generic descriptor-driven propagation of SurfaceKinematicState using the
// nominal material stored in SurfaceDescriptor and resolved via
// SurfaceCatalogView.
//
// Every routing decision in this class inspects the actual
// SurfaceDescriptor::kind of the surface at hand (equivalently, the
// StateFamily it implies via stateFamilyOf(), StateFamily.h) at the call
// site. Nothing here names the confined legacy hot-loop-dispatch tag
// (detail/TrackingKernelParameters.h) or a persisted family/kind pair, and
// nothing here names ITS, MFT, a detector ID, a fixed layer count, a source
// position, a workflow, DPL, a writer, or an output type (ADR 0007 decisions
// 7, 8, 10).
namespace o2::itsmft::tracking
{

class Propagator
{
 public:
  // Propagates a forward state with the accepted model used by disk cell and
  // refit operations: full helix transport for |bz| > 0.01f, and linear
  // transport otherwise. The model selection is centralized here so forward
  // production callers do not grow independent propagation shortcuts.
  static bool propagateForward(SurfaceKinematicState& state, float targetZ, float bz,
                               OperationFailureReason& reason) noexcept;

  // State-family conversion.
  // Re-expresses `state` (and, if supplied, its paired `linRef`) in
  // `targetFamily`'s parameter/reference convention, evaluated at the
  // state's *current* reference surface (task requirement 4: a real
  // coordinate/state/covariance transformation, never a family-flag
  // mutation and never a shortcut that pretends the state is already on
  // the target surface). `state.family == targetFamily` is accepted as a
  // no-op success.
  //
  // The conversion is a first-order (linearized-Jacobian) reparametrization
  // between two different "one coordinate is the fixed reference, the rest
  // are fitted parameters" gauges (Barrel fixes radius+frame angle; Forward
  // fixes only z). Converting necessarily moves which coordinate is fixed,
  // so the coordinate that *becomes* the new fixed reference has its own
  // fitted variance discarded (exactly as every family already discards its
  // own reference coordinate's variance -- SurfaceKinematicState has no
  // variance slot for referenceCoordinate at all), while the coordinate
  // that is newly *freed* into a fitted parameter is assigned a loose,
  // uninformative ceiling variance (reusing this library's own existing
  // "reset to an uninformative diagonal" ceiling constants, see
  // Propagator.cxx) rather than a fabricated precise value. This is a
  // documented, honest engineering choice, not a claim of information-
  // preserving round-trip conversion.
  //
  // Preserves absCharge, PID, and every field of `state` outside the
  // family-specific parameter convention. Transactional: both `state` and
  // `linRef` (if supplied) are left completely unchanged, byte-for-byte, on
  // any failure.
  static bool convertFamily(SurfaceKinematicState& state, SurfaceLinearizationReference* linRef,
                            StateFamily targetFamily, OperationFailureReason& reason) noexcept;

  // Propagation to a measurement.
  // Contract (task spec items 1-7):
  //  1. `targetSurface`/`targetMeasurement` are the explicit target
  //     surface/measurement context; `state`/`linRef` are the explicit
  //     current/source-surface context.
  //  2. Compatibility is `state.family == stateFamilyOf(targetSurface.kind)`.
  //  3. On mismatch, `state`/`linRef` are converted (convertFamily above,
  //     targeting stateFamilyOf(targetSurface.kind)) before any rotation or
  //     propagation is attempted.
  //  4. See convertFamily's own doc: a real transform, never a flag flip.
  //  5. charge/PID/reference/linearization state survive intact -- assigned
  //     only by the barrel::/forward:: primitives' own already-documented
  //     transactional contracts (BarrelSurfaceStateOperations.h /
  //     ForwardSurfaceStateOperations.h), which this function composes but
  //     does not reimplement.
  //  6. Material is `targetSurface.material`, applied through the existing
  //     PID/absCharge-aware material::calculateMaterialPhysics() kernel via
  //     barrel::/forward::correctForMaterial -- never a second material
  //     formula (see the class doc above).
  //  7. Conversion + rotate/propagate + material + chi2 gate + update +
  //     optional reference shift all run in local scratch; `state`/
  //     `linRef`/`chi2` are committed only together, on complete success.
  //     `reason` is always written on failure.
  //
  // Chi2 hardening on entry is shared by all native refit leaves: the incoming
  // `chi2` accumulator
  // must be finite and non-negative; `maxChi2` is validated the same way,
  // but only when `chi2GateEnabled` is true.
  static bool propagateToMeasurement(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef,
                                     const SurfaceDescriptor& targetSurface, const SurfaceMeasurement& targetMeasurement,
                                     float bz, material::MaterialTraversalDirection direction,
                                     bool chi2GateEnabled, float maxChi2, float& chi2,
                                     bool shiftReferenceToMeasurement, OperationFailureReason& reason) noexcept;

  // Shared cylinder/disk leg orchestration.
  // Shared descriptor-driven leg orchestration with the identical per-slot
  // contract
  // (a slot with an invalid `cluster` is a hole, skipped before its
  // `surface`/`surfaceCatalog` association is even inspected; a present
  // slot's `measurement.surface` is validated against `surfaceCatalog`
  // exactly as that function documents, failing with
  // OperationFailureReason::InvalidSurfaceCatalogAssociation; the same chi2
  // hardening; the same scratch-then-commit transactionality; an empty or
  // all-hole leg is unconditional success) -- but selects
  // Propagator::propagateToMeasurement per slot from the slot's own
  // resolved SurfaceDescriptor instead of a compile-time family helper, so
  // one call serves a leg mixing Barrel and Forward slots without any
  // family branch of its own (ADR 0007 decision 10's shared orchestration).
  static bool driveRefitLeg(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef,
                            float& chi2, uint32_t& acceptedHitCount,
                            gsl::span<const SurfaceMeasurement> orderedSlots, SurfaceCatalogView surfaceCatalog,
                            float bz, material::MaterialTraversalDirection direction,
                            bool shiftReferenceToMeasurement, float maxChi2, OperationFailureReason& reason) noexcept;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_PROPAGATOR_H_ */
