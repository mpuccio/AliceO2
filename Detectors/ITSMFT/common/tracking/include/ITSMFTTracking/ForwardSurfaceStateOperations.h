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
#include "ITSMFTTracking/SurfaceLinearizationReference.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"
#include "ITSMFTTracking/TransitionPolicy.h"

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
// This is a retained characterization boundary: it is not redirected to the
// PID/absCharge-aware physics below, even where both are available.
bool correctForMaterial(SurfaceKinematicState& state, float xOverX0, OperationFailureReason& reason) noexcept;

// PID/absCharge-aware material correction using the shared scalar kernel
// (MaterialPhysics.h): derives physical momentum from the pre-material
// state, evaluates calculateMaterialPhysics() with the state's own PID and
// absCharge, and projects the resulting energy-loss/scattering/straggling
// contributions onto the packed covariance. This is an independent overload
// from the legacy correctForMaterial() above, which remains unchanged.
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

// Same-family compatibility chi2 between two forward states expressed at the
// same reference Z. Residuals are the direct (unwrapped) differences of
// (X, Y, Phi, Tanl, InvQPt), preserving the legacy same-family behavior of
// detail::mftFwdStateChi2. Neither input is mutated and chi2 is unchanged on
// failure.
bool stateChi2(const SurfaceKinematicState& reference, const SurfaceKinematicState& candidate, float& chi2,
               OperationFailureReason& reason) noexcept;

#ifndef GPUCA_GPUCODE

// Slice A (Stage-B design report Sec 8/11), migrated in the Stage-B
// normalized-measurement slice to consume SurfaceMeasurement directly:
// builds the initial *outer*-anchored SurfaceKinematicState seed for a
// forward/disk three-hit Cell candidate, transcribing the current
// closed-form direction/covariance initialization at the start of
// buildCellSeed<DiskDisk> (TransitionPolicyOperations.cxx) /
// detail::mftFwdFitCellClusters (MFTFwdTrackHelpers.h) directly onto
// SurfaceKinematicState in float arithmetic only. The legacy initializer
// stores its seed in a ROOT::Math::SVector<double,5>/
// SMatrix<double,5,5,MatRepSym<double,5>> pair; this operation reproduces
// the identical formula and strict-boundary/fallback behavior without ever
// holding a double-precision parameter or covariance value. Production must
// never construct the legacy MFT forward track-parametrization-with-error
// type this reproduces.
//
// Input order matches this operation family's existing {inner, middle,
// outer} contract, never inferred from numeric layer/SurfaceId/radius/z:
// measurementInner/measurementMiddle/measurementOuter supply
// `measurement.global` (x, y, z); measurementOuter additionally supplies the
// outer measurement's measured (u, v) covariance (`measurement.covariance`)
// used to seed the diagonal. trackletMinPt is the same configured
// minimum-pT scale buildCellSeed<DiskDisk> already reads from
// DiskDiskPolicyParams; its established `(trackletMinPt > 0.f) ?
// 1.f/trackletMinPt : 0.f` fallback is preserved verbatim, not
// re-validated here.
//
// Strict boundary (preserved exactly, and reported through the dedicated
// OperationFailureReason::SeedGeometryDegenerate rather than a generic
// non-finite reason -- these are hard rejections the legacy initializer
// already applies before any direction estimate is computed, not NaN/Inf
// artifacts of arithmetic that ran anyway): rejects unless
// measurementInner.global.z is strictly greater than
// measurementOuter.global.z by more than 1e-6f, and unless every one of the
// inner-middle and inner-outer transverse/longitudinal separations exceeds
// its established 1e-6f minimum.
//
// Output anchor/reference frame: Forward family, alpha == 0 (forward has no
// rotation frame), referenceCoordinate == measurementOuter.frame.q -- the
// *outer* measurement's own reference z (the accepted disk adapter contract
// guarantees frame.q == global.z). This is deliberately not the Cell's
// eventual inner-anchored frame (anchor contract, design report Sec 5): this
// operation produces only the initial outer seed that the existing
// outer->middle->inner buildCellSeed sequence subsequently attaches inward,
// hit by hit, to reach the Cell's actual innermost-surface anchor. Calling
// buildSeed alone does not produce a complete Cell state.
//
// Compatibility hypothesis (absCharge, pid): supplied by the caller and
// never defaulted here -- see barrel::buildSeed's identical rationale.
// Unchecked caller precondition: absCharge/pid/trackletMinPt are
// caller-fixed configuration, not independently validated here. flags is
// always set to 0.
//
// Failure precedence: (1) OperationFailureReason::NonFiniteInput if any raw
// input is not finite; (2) OperationFailureReason::SeedGeometryDegenerate if
// the strict-boundary geometry checks above reject; (3)
// OperationFailureReason::NonFiniteOutput if the fully constructed candidate
// state is still not finite despite passing (1) and (2).
//
// Transactional: constructed entirely in local scratch; outState is
// committed only on complete success. On any failure outState is left
// exactly as passed in, byte-for-byte.
bool buildSeed(const SurfaceMeasurement& measurementInner, const SurfaceMeasurement& measurementMiddle,
               const SurfaceMeasurement& measurementOuter, float bz, float trackletMinPt,
               uint8_t absCharge, o2::track::PID pid,
               SurfaceKinematicState& outState, OperationFailureReason& reason) noexcept;

// Stage-B refit-primitive slice: builds the initial anchor-selected
// SurfaceKinematicState seed for a forward/disk three-hit candidate.
// `SeedAnchor::Outer` reproduces `buildSeed` above exactly (byte-for-byte;
// that overload now delegates here). `SeedAnchor::Inner` anchors at
// measurementInner's own frame/reference/covariance instead
// (`referenceCoordinate == measurementInner.frame.q`, `parameters[0]/[1] ==
// measurementInner.global.x/y`, covariance diagonal seeded from
// measurementInner.covariance) using the same {inner, middle, outer}
// physical ordering and physically identical direction estimate. Unlike
// barrel::buildSeed's Inner anchor, no sign flip is applied to phi/tanl/
// invQPt: there is no legacy MFT "reverse" seed formula to reproduce, and
// this operation's closed-form phi/tanl estimate (secant lines through the
// same three global positions) is anchor-symmetric by construction --
// changing which measurement supplies the reference frame does not change
// the estimated direction of the same physical trajectory. The strict
// z-ordering/degenerate-separation boundary (SeedGeometryDegenerate) is
// unchanged: it is a physical-ordering check on the three hits, not an
// anchor-dependent one. An unrecognized `SeedAnchor` value fails with
// `OperationFailureReason::InvalidSeedAnchor` (not NonFiniteInput -- the
// raw measurement/bz/absCharge/pid inputs may be perfectly well-formed;
// it is the anchor selector itself that is invalid) and leaves `outState`
// unchanged.
// Same failure-precedence, absCharge/pid contract, and scratch-then-commit
// transactionality as `buildSeed` above.
bool buildSeed(SeedAnchor anchor, const SurfaceMeasurement& measurementInner, const SurfaceMeasurement& measurementMiddle,
               const SurfaceMeasurement& measurementOuter, float bz, float trackletMinPt,
               uint8_t absCharge, o2::track::PID pid,
               SurfaceKinematicState& outState, OperationFailureReason& reason) noexcept;

// Stage-B refit-primitive slice (linRef-aware propagate): the Disk
// counterpart to barrel::propagate(state, linRef, ...) above, adopting the
// same linearization structure -- propagate the reference trajectory
// non-linearly with the accepted forward Model, evaluate that Model's own
// Jacobian at the reference's pre-propagation parameters (rather than at
// `state`'s own, which is what the plain `propagate<Model>` above does),
// then transport `state` as reference_after + Jacobian*(state_before -
// reference_before) for parameters and as Jacobian*Cov*Jacobian^T
// (congruence) for covariance. There is no legacy MFT linRef oracle to
// transcribe (the legacy MFT forward track-parametrization-with-error type
// has no linearization-reference refit counterpart), so this reuses the
// accepted forward float-native propagation mathematics (the same
// position-update and Jacobian formulas `propagate<Model>` already
// implements) rather than an independent derivation.
//
// `linRef` carries no absCharge/pid of its own (see
// SurfaceLinearizationReference.h); its own propagation and Jacobian
// evaluation always use `state.absCharge`/`state.pid` implicitly through
// the shared per-Model helper, matching barrel's identical substitution.
//
// Same per-Model preconditions as `propagate<Model>` above (Linear/
// Quadratic require only tanl!=0 for a nonzero step; Helix additionally
// requires bz!=0 and invQPt!=0), evaluated against the *reference's*
// pre-propagation parameters. Optimized combines a Helix position update
// for both `linRef` and `state` with a Quadratic-Jacobian covariance
// transport, exactly mirroring the plain `propagate<Optimized>` composition
// above.
//
// Fitted-state/linRef pairing precondition, checked before the above:
// exact `state.family == linRef.family`
// (OperationFailureReason::SourceFamilyMismatch) and exact
// `state.referenceCoordinate == linRef.referenceCoordinate`
// (OperationFailureReason::ReferenceCoordinateMismatch) -- both bit
// comparisons, not tolerance checks, matching barrel::rotate/propagate's
// identical contract. No alpha check: Forward's alpha is always 0/unused.
// Parameters may differ between `state` and `linRef` -- that is the
// entire purpose of a linearization reference -- but the anchor may not.
//
// Full scratch-then-commit transactionality: both `state` and `linRef` are
// left completely unchanged (byte-for-byte) on any failure.
template <PropagationModel Model>
bool propagate(SurfaceKinematicState& state, SurfaceLinearizationReference& linRef, float targetZ, float bz,
               OperationFailureReason& reason) noexcept;

template <>
bool propagate<PropagationModel::Linear>(SurfaceKinematicState&, SurfaceLinearizationReference&, float, float, OperationFailureReason&) noexcept;
template <>
bool propagate<PropagationModel::Quadratic>(SurfaceKinematicState&, SurfaceLinearizationReference&, float, float, OperationFailureReason&) noexcept;
template <>
bool propagate<PropagationModel::Helix>(SurfaceKinematicState&, SurfaceLinearizationReference&, float, float, OperationFailureReason&) noexcept;
template <>
bool propagate<PropagationModel::Optimized>(SurfaceKinematicState&, SurfaceLinearizationReference&, float, float, OperationFailureReason&) noexcept;

// Explicit reference-shift operation required by the future
// ShiftRefToCluster leg option (design report Sec 8/11), the Disk
// counterpart to barrel::shiftReferenceToMeasurement above. The disk
// measurement frame maps normal/measured coordinates as q=z, u=x, v=y
// (Architecture.md Sec 6.5); `parameters[0]/[1]` (X, Y) are the disk
// family's direct position parameters, the same structural role Y/Z play
// for barrel -- so this sets `linRef.parameters[0]/[1]` to
// `measurement.frame.u`/`measurement.frame.v` (equivalently
// `measurement.global.x`/`measurement.global.y` under the disk adapter
// contract), never touching `linRef.referenceCoordinate`, Phi, Tanl, or
// InvQPt. `linRef.alpha` is likewise untouched (always 0 for Forward,
// unused). Fails (leaving `linRef` completely unchanged) if `linRef.family
// != StateFamily::Forward` or if `measurement.frame.u`/`measurement.frame.v`
// is not finite.
bool shiftReferenceToMeasurement(SurfaceLinearizationReference& linRef, const SurfaceMeasurement& measurement,
                                 OperationFailureReason& reason) noexcept;

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking::forward

#endif // ALICEO2_ITSMFT_TRACKING_FORWARDSURFACESTATEOPERATIONS_H_
