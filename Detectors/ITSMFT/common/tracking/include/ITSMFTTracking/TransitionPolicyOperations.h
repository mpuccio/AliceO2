// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYOPERATIONS_H_
#define ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYOPERATIONS_H_

#include <array>

#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/TransitionPolicyState.h"

#ifndef GPUCA_GPUCODE
#include <gsl/span>

#include "DetectorsBase/Propagator.h"

namespace o2::dataformats
{
template <typename Stamp>
class Vertex;
}

namespace o2::its
{
class TimeEstBC;
using Vertex = o2::dataformats::Vertex<TimeEstBC>;
struct TrackingFrameInfo;
struct Cluster;
} // namespace o2::its

namespace o2::itsmft
{
template <int NLayers>
class IndexTableUtils;
}
#endif

namespace o2::itsmft::tracking
{

/// D007 policy-boundary operation (Architecture.md 10): true if `nextCell`
/// may extend a road built on `currentCell`. cellsAreCompatible and attachHit
/// are now migrated to this boundary; projectSearchWindow, buildCellSeed, and
/// finalRefit remain TrackerTraits-internal until later slices. Dispatch is a
/// compile-time tag switch -- one specialization per Tag, selected by the
/// caller once outside any candidate/neighbour loop -- never a runtime detector
/// branch inside this operation or its callers' hot loops. Instantiating the
/// primary template for an unsupported tag is a compile error rather than a
/// silent fallback (mirrors TransitionPolicyTraits).
///
/// Host-only for this CPU migration slice: the disk (DiskDisk) specialization
/// is defined in TransitionPolicyOperations.cxx against the existing MFT
/// forward-state helpers (TrackParCovFwd::propagateToZhelix/Linear and
/// friends), which are themselves host-only. This operation must not be
/// declared or assumed device-compatible until a device-capable forward
/// state/propagator exists; do not add GPU qualifiers here without first
/// making that dependency chain device-compatible. The public header
/// deliberately does not include MFTFwdTrackHelpers.h, so the MFT-specific
/// helper stays behind this implementation boundary rather than leaking MFT
/// constants/TimeFrame/MFTCATrack dependencies into a common policy header.
template <TransitionPolicyTag Tag>
bool cellsAreCompatible(const CellSeedTpl<typename TransitionPolicyTraits<Tag>::SeedState>& currentCell,
                        const CellSeedTpl<typename TransitionPolicyTraits<Tag>::SeedState>& nextCell,
                        float bz,
                        const typename TransitionPolicyTraits<Tag>::Params& params);

/// Barrel formula: rotate/propagate `nextCell` into `currentCell`'s frame and
/// accept if the predicted chi2 stays within the bound. `nextCell` is never
/// mutated -- the propagated state is a local copy -- so callers may pass a
/// stored cell directly without pre-copying it for this call. Defined out of
/// line in TransitionPolicyOperations.cxx.
template <>
bool cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(
  const CellSeedTpl<o2::track::TrackParCovF>& currentCell,
  const CellSeedTpl<o2::track::TrackParCovF>& nextCell,
  float bz,
  const CylinderCylinderPolicyParams& params);

/// Disk formula: delegates to the existing MFT forward-state helper, which
/// already performs its own internal copy/propagation; neither input cell is
/// mutated. Defined out of line in TransitionPolicyOperations.cxx, the only
/// translation unit permitted to include MFTFwdTrackHelpers.h on behalf of
/// this policy operation.
template <>
bool cellsAreCompatible<TransitionPolicyTag::DiskDisk>(
  const CellSeedTpl<o2::track::TrackParCovFwd>& currentCell,
  const CellSeedTpl<o2::track::TrackParCovFwd>& nextCell,
  float bz,
  const DiskDiskPolicyParams& params);

#ifndef GPUCA_GPUCODE

/// Stage-A compatibility state for tracklet projection. `fromLayer` and
/// `toLayer` are legacy layout-local indices into TimeFrame, IndexTableUtils,
/// TrackingParameters, and MFT helper arrays. They are not global SurfaceIds;
/// replacing them belongs to the later layout-driven indexing migration.
template <TransitionPolicyTag Tag>
struct TrackletProjectionState;

template <>
struct TrackletProjectionState<TransitionPolicyTag::CylinderCylinder> {
  int fromLayer;
  int toLayer;
  float meanDeltaR;
  float targetMinR;
  float targetMaxR;
  float sourcePositionResolution;
  float transitionMSAngle;
  float transitionPhiCut;
};

template <>
struct TrackletProjectionState<TransitionPolicyTag::DiskDisk> {
  int fromLayer;
  int toLayer;
  float fromZ;
  float toZ;
  float meanDeltaZ;
  float sourceReferenceRadius;
  float transitionMSAngle;
  float transitionBendingAngle;
};

template <TransitionPolicyTag Tag>
struct TrackletSearchWindow;

template <>
struct TrackletSearchWindow<TransitionPolicyTag::CylinderCylinder> {
  int4 bins;
  float tanLambda;
  float sigmaZ;
  float phiCut;
  float nSigmaCut;

  bool acceptCandidate(const o2::its::Cluster& source,
                       const o2::its::Cluster& target,
                       float& tanLambdaOut) const;
};

template <>
struct TrackletSearchWindow<TransitionPolicyTag::DiskDisk> {
  int4 bins;
  float xProj;
  float yProj;
  float sigmaX;
  float sigmaY;
  float meanDeltaZ;
  float nSigmaCut;

  bool acceptCandidate(const o2::its::Cluster& source,
                       const o2::its::Cluster& target,
                       float& tanLambdaOut) const;
};

/// Projects one source cluster/vertex pair and builds the family-specific LUT
/// window. Implemented and explicitly instantiated in the host translation
/// unit: DiskDisk reuses MFTFwdTrackHelpers, which must not leak into this
/// declarations-only header. On rejection `out` is left unchanged.
template <TransitionPolicyTag Tag, int NLayers>
bool projectSearchWindow(const o2::its::Cluster& source,
                         const o2::its::TrackingFrameInfo& sourceHit,
                         const o2::its::Vertex& vertex,
                         const TrackletProjectionState<Tag>& transitionState,
                         float bz,
                         const o2::itsmft::IndexTableUtils<NLayers>& indexUtils,
                         const typename TransitionPolicyTraits<Tag>::Params& params,
                         TrackletSearchWindow<Tag>& out);

/// D007 attach-hit policy operation. The typed family state and parameter
/// block are selected once by the caller's outer policy dispatch. `xOverX0`
/// is the already-selected material budget for the hit surface; no legacy
/// TrackingParameters object or detector identity crosses this boundary.
/// `o2::its::TrackingFrameInfo` is a temporary Gate 3 compatibility boundary,
/// not a detector-neutral measurement contract; production normalized loading
/// must eventually supply SurfaceMeasurement directly.
///
/// Host-only: CylinderCylinder calls the host Propagator singleton and
/// DiskDisk calls the host forward-state propagation/update chain. The two
/// specializations deliberately preserve their existing failure mutation
/// contracts: the cylinder state is modified as each successful step is
/// applied, while the disk state and chi2 are committed only after the full
/// attachment succeeds.
template <TransitionPolicyTag Tag>
bool attachHit(typename TransitionPolicyTraits<Tag>::SeedState& state,
               const o2::its::TrackingFrameInfo& hit,
               float xOverX0,
               o2::base::PropagatorF::MatCorrType corrType,
               float bz,
               float& chi2,
               const typename TransitionPolicyTraits<Tag>::Params& params);

template <>
bool attachHit<TransitionPolicyTag::CylinderCylinder>(
  o2::track::TrackParCovF& state,
  const o2::its::TrackingFrameInfo& hit,
  float xOverX0,
  o2::base::PropagatorF::MatCorrType corrType,
  float bz,
  float& chi2,
  const CylinderCylinderPolicyParams& params);

template <>
bool attachHit<TransitionPolicyTag::DiskDisk>(
  o2::track::TrackParCovFwd& state,
  const o2::its::TrackingFrameInfo& hit,
  float xOverX0,
  o2::base::PropagatorF::MatCorrType corrType,
  float bz,
  float& chi2,
  const DiskDiskPolicyParams& params);

/// D007 policy-boundary operation (Architecture.md Sec 10): seeds and fits
/// the three-cluster CA cell state for one candidate cell. `xOverX0` is
/// always ordered {inner, middle, outer}; each specialization reads only the
/// slots its own attachment order consumes (documented per specialization
/// below) -- the caller resolves these floats from the same per-iteration
/// AttachHitPolicyConfigView.layerxX0 binding used by attachHit, outside any
/// candidate/neighbour loop. The MFT geometric road pre-cut (legacy
/// CellRoadRCut / detail::validateMFTCellClusters) is deliberately not part
/// of this operation: it depends on nominal per-surface layer position, not
/// on the measurements passed here, and remains a TrackerTraits-owned
/// orchestration guard evaluated before this operation is called -- exactly
/// as the tracklet timestamp/deltaTanLambda checks already are for
/// cellsAreCompatible/attachHit. This operation does not add any input
/// validation beyond what the legacy inline implementations perform; valid
/// policy parameters and material arrays are preconditions established by
/// the existing one-time binding/initialization boundary, not rechecked
/// here.
///
/// Output contract: `outState`/`chi2` are committed only if this operation
/// returns true. On any rejection (fit-precondition, rotation, propagation,
/// material, update, or chi2-cut failure) `outState` and `chi2` are left
/// exactly as passed in -- the fit runs into local scratch that is discarded
/// on failure, never partially exposed to the caller. This is a stricter,
/// explicitly reusable contract than the current inline
/// TrackerTraits::computeLayerCells, whose equivalent locals happen not to
/// be read by any caller after a failed fit today.
///
/// Host-only for this slice, like attachHit: CylinderCylinder calls
/// o2::its::track::buildTrackSeed and TrackParCovF::rotate/propagateTo/
/// correctForMaterial/update directly -- no o2::base::Propagator singleton,
/// no MatCorrType; correctForMaterial is applied unconditionally, exactly as
/// the current inline barrel branch does. This must not be unified with
/// attachHit<CylinderCylinder>'s Propagator-based material path. DiskDisk
/// calls the existing host-only MFT forward-state helpers (ROOT::Math
/// based). Defining this out of line in a host translation unit is not
/// itself evidence of device readiness: a device-callable implementation
/// (in particular for CylinderCylinder, whose underlying TrackParCovF
/// primitives already have a GPU cell-building kernel counterpart in
/// ITS/tracking/GPU) requires separate verification before this operation is
/// declared or compiled for device code. Both specializations therefore stay
/// behind this GPUCA_GPUCODE guard, unlike cellsAreCompatible, which is left
/// unchanged by this slice.
///
/// o2::its::Cluster and o2::its::TrackingFrameInfo are forward-declared
/// above and used here only as the documented, time-boxed Gate 3
/// compatibility boundary (mirrors attachHit's TrackingFrameInfo use); this
/// header does not include ITStracking/Cluster.h.
template <TransitionPolicyTag Tag>
bool buildCellSeed(const o2::its::Cluster& clusterInner,
                   const o2::its::Cluster& clusterMiddle,
                   const o2::its::Cluster& clusterOuter,
                   const o2::its::TrackingFrameInfo& hitInner,
                   const o2::its::TrackingFrameInfo& hitMiddle,
                   const o2::its::TrackingFrameInfo& hitOuter,
                   const std::array<float, 3>& xOverX0,
                   float bz,
                   typename TransitionPolicyTraits<Tag>::SeedState& outState,
                   float& chi2,
                   const typename TransitionPolicyTraits<Tag>::Params& params);

/// Barrel formula: o2::its::track::buildTrackSeed(clusterInner,
/// clusterMiddle, hitOuter, bz) followed by middle-then-inner
/// rotate/propagateTo/correctForMaterial/update (xOverX0[1] then xOverX0[0];
/// xOverX0[2] is unused -- the outer surface contributes only through
/// hitOuter inside buildTrackSeed, never a separate attach step). The chi2
/// cut (params.maxChi2ClusterAttachment) is enforced only on the final
/// (inner) step, matching the legacy `!iC` condition. Defined out of line in
/// TransitionPolicyOperations.cxx.
template <>
bool buildCellSeed<TransitionPolicyTag::CylinderCylinder>(
  const o2::its::Cluster& clusterInner,
  const o2::its::Cluster& clusterMiddle,
  const o2::its::Cluster& clusterOuter,
  const o2::its::TrackingFrameInfo& hitInner,
  const o2::its::TrackingFrameInfo& hitMiddle,
  const o2::its::TrackingFrameInfo& hitOuter,
  const std::array<float, 3>& xOverX0,
  float bz,
  o2::track::TrackParCovF& outState,
  float& chi2,
  const CylinderCylinderPolicyParams& params);

/// Disk formula: the same closed-form outward direction estimate and
/// outer/middle/inner MFT Kalman attachment (detail::mftFwdAttachCluster) as
/// detail::mftFwdFitCellClusters, reading its three clusters/hits directly
/// instead of through a TimeFrame. xOverX0 is consumed in outer/middle/inner
/// order (xOverX0[2], then [1], then [0]); the chi2 cut is enforced only on
/// the final (inner) step. Defined out of line in
/// TransitionPolicyOperations.cxx, the only translation unit permitted to
/// include MFTFwdTrackHelpers.h on behalf of this policy operation.
template <>
bool buildCellSeed<TransitionPolicyTag::DiskDisk>(
  const o2::its::Cluster& clusterInner,
  const o2::its::Cluster& clusterMiddle,
  const o2::its::Cluster& clusterOuter,
  const o2::its::TrackingFrameInfo& hitInner,
  const o2::its::TrackingFrameInfo& hitMiddle,
  const o2::its::TrackingFrameInfo& hitOuter,
  const std::array<float, 3>& xOverX0,
  float bz,
  o2::track::TrackParCovFwd& outState,
  float& chi2,
  const DiskDiskPolicyParams& params);

/// Gate 3 transition-preparation slice (Architecture.md Sec 10/10.1):
/// per-layer multiple-scattering angle, one specialization per Tag. Called
/// once per iteration/layer by the caller's outer policy dispatch, never per
/// candidate. Each specialization receives every physical fact it needs
/// explicitly; in particular the DiskDisk specialization does not call
/// mftLayerZ()/LayerZCoordinate() or otherwise infer MFT geometry -- the
/// caller supplies `referenceCoordinate` via the explicit, time-boxed
/// compatibility binding below (bindLegacyMFTReferenceCoordinates()).
/// Instantiating the primary template for an unsupported tag is a compile
/// error rather than a silent fallback (mirrors TransitionPolicyTraits).
template <TransitionPolicyTag Tag>
struct LayerScatteringInputs;

/// Barrel: mass/momentum radiation-length formula (math_utils::MSangle),
/// unchanged from the frozen ITS expression.
template <>
struct LayerScatteringInputs<TransitionPolicyTag::CylinderCylinder> {
  float layerxX0;
};

/// Disk: adds the incidence-angle (pseudo-tan(lambda)) correction that
/// mftLayerMSAngle() already applies. `layerRadius` is
/// TrackingParameters::LayerRadii[layer] reused as MFT's nominal transverse
/// reference radius at `referenceCoordinate` (z) -- the same legacy field
/// overload as today, not something this slice changes.
template <>
struct LayerScatteringInputs<TransitionPolicyTag::DiskDisk> {
  float layerxX0;
  float layerRadius;
  float referenceCoordinate;
};

template <TransitionPolicyTag Tag>
float layerMultipleScatteringAngle(const LayerScatteringInputs<Tag>& inputs, float trackletMinPt);

template <>
float layerMultipleScatteringAngle<TransitionPolicyTag::CylinderCylinder>(
  const LayerScatteringInputs<TransitionPolicyTag::CylinderCylinder>& inputs, float trackletMinPt);

template <>
float layerMultipleScatteringAngle<TransitionPolicyTag::DiskDisk>(
  const LayerScatteringInputs<TransitionPolicyTag::DiskDisk>& inputs, float trackletMinPt);

/// Explicit, time-boxed Gate 3 compatibility binding: supplies the legacy MFT
/// half-disk z coordinates (o2::mft::constants::mft::LayerZCoordinate()) used
/// by legacy mftLayerMSAngle() and required for
/// layerMultipleScatteringAngle<DiskDisk> to reproduce the accepted 91-track /
/// hash 826dc653cd936a472929c600c97c140b MFT common-CA baseline bit-for-bit.
/// Deliberately NOT SurfaceDescriptor::referenceCoordinate: the real-geometry
/// catalog values differ from these legacy constants, and substituting one
/// for the other would silently change the accepted physics baseline rather
/// than migrate it -- that substitution is a separate, later, independently
/// replay-validated decision. Isolated here so the generic DiskDisk operation
/// above never references MFT constants directly; a future disk layout
/// supplies its own LayerScatteringInputs<DiskDisk>::referenceCoordinate
/// without changing layerMultipleScatteringAngle's signature or body.
struct DiskDiskReferenceCoordinateView {
  gsl::span<const float> perLayerReferenceZ; // legacy layout-local layer index
  bool isValid(size_t expectedLayers) const noexcept { return perLayerReferenceZ.size() >= expectedLayers; }
};

/// Bound from static-storage-duration legacy constants (see .cxx) -- the
/// returned view never dangles and needs no per-iteration staging.
DiskDiskReferenceCoordinateView bindLegacyMFTReferenceCoordinates() noexcept;

/// Narrowly isolated legacy-arithmetic compatibility step (integration review
/// finding): the two legacy branches this slice replaces compare
/// `0.5 * oneOverR` (CylinderCylinder: `0.5` is a double literal, so the
/// multiply is evaluated in double after promoting `oneOverR`) against
/// `0.5f * oneOverR` (DiskDisk: stays in float) before the same `>= 1.f / r2`
/// clamp. This is not a policy/physics difference -- it is an accidental
/// literal difference in the frozen legacy code -- but it is preserved
/// bit-for-bit rather than canonicalized, since no evidence yet justifies
/// unifying the precision and doing so would require its own boundary tests
/// and replay evidence. Preserving it here, isolated, avoids duplicating the
/// (genuinely shared) remainder of the transition formula just to keep this
/// one literal. The caller owns the loop-carried `oneOverR` and must call
/// this once per transition, in increasing legacy transitionId order,
/// threading the returned value into the next call -- `oneOverR` is not
/// reset per transition.
template <TransitionPolicyTag Tag>
float clampTransitionCurvature(float oneOverR, float r2) noexcept;

template <>
float clampTransitionCurvature<TransitionPolicyTag::CylinderCylinder>(float oneOverR, float r2) noexcept;

template <>
float clampTransitionCurvature<TransitionPolicyTag::DiskDisk>(float oneOverR, float r2) noexcept;

/// Tag-independent: the arithmetic after the family-specific curvature clamp
/// (clampTransitionCurvature) is verified byte-identical between
/// CylinderCylinder and DiskDisk in the frozen legacy code (both the common
/// TimeFrame's former if-constexpr branches and, for the barrel case, the
/// separate frozen ITS-only TimeFrame::initialise()). It is therefore not a
/// policy specialization point without new evidence of divergence; adding a
/// Tag parameter here would duplicate a formula that does not currently
/// differ. Performs the integrated multiple-scattering sum over
/// `[fromLayer, toLayer)` (half-open, legacy layout-local layer indices) and
/// the transverse bending/phi-acceptance half-width, clamped to [., PI/2].
/// Never fails (pure, total float arithmetic, matching legacy: degenerate
/// inputs such as `r2==0` propagate through to a non-finite result exactly as
/// today, they are not rejected here).
struct TransitionScatteringBendingPrep {
  float msAngle; // integrated multiple-scattering angle across the transition
  float phiCut;  // combined phi-acceptance / bending half-width, radians; called
                 // "transitionPhiCut" by CylinderCylinder and
                 // "transitionBendingAngle" by DiskDisk downstream -- same
                 // quantity, two domain names
};

TransitionScatteringBendingPrep prepareTransitionScatteringAndBending(
  gsl::span<const float> perLayerMSAngle, int fromLayer, int toLayer,
  float r1, float r2, float clampedOneOverR, float res1, float res2) noexcept;

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYOPERATIONS_H_ */
