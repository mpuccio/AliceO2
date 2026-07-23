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
#include <cmath>
#include <limits>

#include <gsl/span>

#include "DetectorsBase/Propagator.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"

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
/// block are selected once by the caller's outer policy dispatch. `material`
/// is the already-selected authoritative nominal material for the hit
/// surface (SurfaceDescriptor::material, resolved once per iteration into
/// TrackerTraits::mLayerMaterial); no legacy TrackingParameters object or
/// detector identity crosses this boundary. `o2::its::TrackingFrameInfo` is a
/// temporary Gate 3 compatibility boundary, not a detector-neutral
/// measurement contract; production normalized loading must eventually
/// supply SurfaceMeasurement directly.
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
               const NominalSurfaceMaterial& material,
               o2::base::PropagatorF::MatCorrType corrType,
               float bz,
               float& chi2,
               const typename TransitionPolicyTraits<Tag>::Params& params);

/// Barrel formula: calls the same legacy correction as before, now reading
/// its two arguments directly off `material` (xOverX0, arealDensityGPerCm2)
/// instead of recomputing the silicon-equivalent areal density inline --
/// bit-for-bit identical to the previous `xOverX0 * Radl * Rho`, since the
/// catalog population uses the exact same expression (GeometrySurfaceCatalogProvider.cxx).
template <>
bool attachHit<TransitionPolicyTag::CylinderCylinder>(
  o2::track::TrackParCovF& state,
  const o2::its::TrackingFrameInfo& hit,
  const NominalSurfaceMaterial& material,
  o2::base::PropagatorF::MatCorrType corrType,
  float bz,
  float& chi2,
  const CylinderCylinderPolicyParams& params);

/// Disk formula: passes only `material.xOverX0` to the unchanged legacy MFT
/// addMCSEffect path (detail::mftFwdAttachCluster) -- MCS-only, exactly as
/// before. `material.arealDensityGPerCm2` is deliberately not read: this
/// slice does not activate MFT energy loss.
template <>
bool attachHit<TransitionPolicyTag::DiskDisk>(
  o2::track::TrackParCovFwd& state,
  const o2::its::TrackingFrameInfo& hit,
  const NominalSurfaceMaterial& material,
  o2::base::PropagatorF::MatCorrType corrType,
  float bz,
  float& chi2,
  const DiskDiskPolicyParams& params);

/// D007 policy-boundary operation (Architecture.md Sec 10): seeds and fits
/// the three-cluster CA cell state for one candidate cell. `material` is
/// always ordered {inner, middle, outer}; each specialization reads only the
/// slots its own attachment order consumes (documented per specialization
/// below) -- the caller resolves these from the same per-iteration
/// AttachHitPolicyConfigView.layerMaterial binding used by attachHit, outside
/// any candidate/neighbour loop. The MFT geometric road pre-cut
/// (CellRoadRCut / passesCellRoadPrecut<DiskDisk>, below) is deliberately not
/// part of this operation: it depends on nominal per-surface layer position, not
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
                   const std::array<NominalSurfaceMaterial, 3>& material,
                   float bz,
                   typename TransitionPolicyTraits<Tag>::SeedState& outState,
                   float& chi2,
                   const typename TransitionPolicyTraits<Tag>::Params& params);

/// Barrel formula: o2::its::track::buildTrackSeed(clusterInner,
/// clusterMiddle, hitOuter, bz) followed by middle-then-inner
/// rotate/propagateTo/correctForMaterial/update (material[1] then
/// material[0]; material[2]/outer is unused -- the outer surface contributes
/// only through hitOuter inside buildTrackSeed, never a separate attach
/// step). correctForMaterial reads material.xOverX0/arealDensityGPerCm2
/// directly, bit-for-bit identical to the previous inline
/// `xOverX0 * Radl * Rho` (see attachHit<CylinderCylinder>'s doc). The chi2
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
  const std::array<NominalSurfaceMaterial, 3>& material,
  float bz,
  o2::track::TrackParCovF& outState,
  float& chi2,
  const CylinderCylinderPolicyParams& params);

/// Disk formula: the same closed-form outward direction estimate and
/// outer/middle/inner MFT Kalman attachment (detail::mftFwdAttachCluster) as
/// detail::mftFwdFitCellClusters, reading its three clusters/hits directly
/// instead of through a TimeFrame. `material` is consumed in outer/middle/inner
/// order (material[2], then [1], then [0]), passing only each slot's
/// `.xOverX0` to the unchanged MCS-only legacy path (arealDensityGPerCm2 is
/// not read here, matching attachHit<DiskDisk>); the chi2 cut is enforced
/// only on the final (inner) step. Defined out of line in
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
  const std::array<NominalSurfaceMaterial, 3>& material,
  float bz,
  o2::track::TrackParCovFwd& outState,
  float& chi2,
  const DiskDiskPolicyParams& params);

/// Stage-B Slice B (design report Sec 8/11): native SurfaceKinematicState
/// overload of the same-named operation above, distinguished purely by
/// argument type (outState here is SurfaceKinematicState&, never
/// TrackParCovF/TrackParCovFwd) -- ordinary overload resolution on an
/// explicit Tag call site, not a TransitionPolicyTraits<Tag>::SeedState
/// change. Additive and unwired: no production caller exists yet: this
/// slice only composes barrel::/forward:: primitives (buildSeed, rotate/
/// propagate, the PID/absCharge-aware correctForMaterial overload,
/// predictedChi2, update) that are already characterized against their own
/// legacy oracles; it does not itself construct or convert through
/// TrackParCovF/TrackParCovFwd.
///
/// Both specializations: composed from barrel::/forward::buildSeed for the
/// initial outer anchor, then step inward hit-by-hit exactly as the legacy
/// operation above does (same material-slot order, same chi2-cut step
/// placement), applying MaterialTraversalDirection::OppositeMomentum via the
/// PID/absCharge-aware correctForMaterial(state, IntegratedMaterialBudget,
/// direction) overload -- never the legacy MCS-only
/// forward::correctForMaterial(state, xOverX0, reason) overload. `absCharge`/
/// `pid` are caller-supplied compatibility values threaded straight into
/// buildSeed (which sets state.absCharge/state.pid/state.flags=0 on the
/// constructed seed; not re-set by this operation). DiskDisk additionally
/// activates PID-aware energy loss/straggling by reading
/// material[*].arealDensityGPerCm2 (unlike the legacy DiskDisk operation
/// above, which is MCS-only) -- an intentional physics difference, not a
/// legacy-equivalence claim.
///
/// Measurement input: the private compatibility projection in
/// TransitionPolicyOperations.cxx converts each o2::its::TrackingFrameInfo
/// hit into the minimum SurfaceMeasurement fields the native update
/// operation reads (barrel: local Y/Z and full yy/yz/zz; forward: global X/Y
/// and the diagonal-only uu/vv, matching the established
/// mftFwdAttachCluster contract) -- this is not a claim that
/// TrackingFrameInfo is a normalized SurfaceMeasurement.
///
/// Transactional and byte-identical-on-failure exactly like the legacy
/// operation above: outState/chi2/reason are the only outputs, and
/// outState/chi2 are committed only on complete success. reason is always
/// written on failure and otherwise left unspecified on success (mirrors
/// every other native operation's own reason contract in this file family).
///
/// Host-only for this slice (requirement: no GPU/device claim): kept behind
/// this GPUCA_GPUCODE guard like the legacy buildCellSeed/attachHit above,
/// even though its own underlying barrel::/forward:: primitives (other than
/// buildSeed) are individually declared unconditionally -- no device
/// readiness has been established for this composed operation as a whole.
template <TransitionPolicyTag Tag>
bool buildCellSeed(const o2::its::Cluster& clusterInner,
                   const o2::its::Cluster& clusterMiddle,
                   const o2::its::Cluster& clusterOuter,
                   const o2::its::TrackingFrameInfo& hitInner,
                   const o2::its::TrackingFrameInfo& hitMiddle,
                   const o2::its::TrackingFrameInfo& hitOuter,
                   const std::array<NominalSurfaceMaterial, 3>& material,
                   float bz,
                   uint8_t absCharge,
                   o2::track::PID pid,
                   SurfaceKinematicState& outState,
                   float& chi2,
                   const typename TransitionPolicyTraits<Tag>::Params& params,
                   OperationFailureReason& reason) noexcept;

/// Barrel formula: barrel::buildSeed(clusterInner, clusterMiddle, hitOuter,
/// bz, absCharge, pid, ...) for the outer anchor, then middle-then-inner
/// rotate/propagate/correctForMaterial/predictedChi2/update using material[1]
/// then material[0] (material[2]/outer unused, matching the legacy operation
/// above exactly). The chi2 cut (params.maxChi2ClusterAttachment) is enforced
/// only on the final (inner) step, with no `< 0.f` rejection -- matching the
/// legacy buildCellSeed<CylinderCylinder> exactly (unlike native attachHit
/// below, which does reject `< 0.f`). Finishes anchored at the inner hit's
/// frame/reference coordinate. Defined out of line in
/// TransitionPolicyOperations.cxx.
template <>
bool buildCellSeed<TransitionPolicyTag::CylinderCylinder>(
  const o2::its::Cluster& clusterInner,
  const o2::its::Cluster& clusterMiddle,
  const o2::its::Cluster& clusterOuter,
  const o2::its::TrackingFrameInfo& hitInner,
  const o2::its::TrackingFrameInfo& hitMiddle,
  const o2::its::TrackingFrameInfo& hitOuter,
  const std::array<NominalSurfaceMaterial, 3>& material,
  float bz,
  uint8_t absCharge,
  o2::track::PID pid,
  SurfaceKinematicState& outState,
  float& chi2,
  const CylinderCylinderPolicyParams& params,
  OperationFailureReason& reason) noexcept;

/// Disk formula: forward::buildSeed(clusterInner, clusterMiddle,
/// clusterOuter, hitOuter, bz, params.trackletMinPt, absCharge, pid, ...) for
/// the outer anchor, then outer/middle/inner propagation with the accepted
/// forward model (threshold dispatch to forward::propagate<Helix>/<Linear>
/// reproducing detail::mftFwdPropagateToZ's own |bz|>0.01f selection exactly
/// -- PropagationModel::Optimized is not used by this slice) and
/// correctForMaterial using material[2], material[1], material[0] in that
/// order -- all three slots consumed (unlike the legacy DiskDisk operation
/// above, which is MCS-only and reads only .xOverX0): this activates
/// PID-aware energy loss and straggling and is intentionally not required to
/// reproduce the legacy MCS-only numerical result. Chi2 cut enforced only on
/// the final (inner) step, no `< 0.f` rejection (matches the legacy
/// buildCellSeed<DiskDisk> exactly). Finishes anchored at the inner hit's z.
/// Defined out of line in TransitionPolicyOperations.cxx.
template <>
bool buildCellSeed<TransitionPolicyTag::DiskDisk>(
  const o2::its::Cluster& clusterInner,
  const o2::its::Cluster& clusterMiddle,
  const o2::its::Cluster& clusterOuter,
  const o2::its::TrackingFrameInfo& hitInner,
  const o2::its::TrackingFrameInfo& hitMiddle,
  const o2::its::TrackingFrameInfo& hitOuter,
  const std::array<NominalSurfaceMaterial, 3>& material,
  float bz,
  uint8_t absCharge,
  o2::track::PID pid,
  SurfaceKinematicState& outState,
  float& chi2,
  const DiskDiskPolicyParams& params,
  OperationFailureReason& reason) noexcept;

/// Stage-B Slice B (design report Sec 8/11): native SurfaceKinematicState
/// overload of the legacy attachHit<Tag> above, distinguished by argument
/// type (state here is SurfaceKinematicState&, never TrackParCovF/
/// TrackParCovFwd) and by omitting the legacy corrType parameter entirely
/// (requirement: no corrType on the new native attachHit API -- the later
/// activation slice rejects non-NONE CylinderCylinder material modes before
/// candidate processing, not this operation). Additive and unwired.
///
/// In-place semantics like the legacy operation: `state`/`chi2` are both
/// input and output, updated only on complete success; on any failure both
/// are left exactly as passed in (byte-for-byte), and `reason` is always
/// written. Composed entirely from barrel::/forward:: primitives: rotate
/// (barrel only)/propagate (accepted forward model for DiskDisk, matching
/// buildCellSeed above), the PID/absCharge-aware correctForMaterial overload
/// with MaterialTraversalDirection::OppositeMomentum, predictedChi2, and
/// update, using the same private TrackingFrameInfo->SurfaceMeasurement
/// projection as native buildCellSeed.
///
/// Host-only for this slice, kept behind this GPUCA_GPUCODE guard like the
/// legacy attachHit/buildCellSeed above.
template <TransitionPolicyTag Tag>
bool attachHit(SurfaceKinematicState& state,
               const o2::its::TrackingFrameInfo& hit,
               const NominalSurfaceMaterial& material,
               float bz,
               float& chi2,
               const typename TransitionPolicyTraits<Tag>::Params& params,
               OperationFailureReason& reason) noexcept;

/// Barrel formula: rotate to hit.alphaTrackingFrame, propagate to
/// hit.xTrackingFrame, PID/absCharge-aware correctForMaterial(material,
/// OppositeMomentum), predictedChi2, reject if predictedChi2 >
/// params.maxChi2ClusterAttachment or predictedChi2 < 0.f (matching the
/// legacy attachHit<CylinderCylinder> exactly, including its `< 0.f` guard),
/// then update. Defined out of line in TransitionPolicyOperations.cxx.
template <>
bool attachHit<TransitionPolicyTag::CylinderCylinder>(
  SurfaceKinematicState& state,
  const o2::its::TrackingFrameInfo& hit,
  const NominalSurfaceMaterial& material,
  float bz,
  float& chi2,
  const CylinderCylinderPolicyParams& params,
  OperationFailureReason& reason) noexcept;

/// Disk formula: propagate to hit.zCoordinate with the accepted forward
/// model (see native buildCellSeed<DiskDisk> doc), PID/absCharge-aware
/// correctForMaterial(material, OppositeMomentum) -- both .xOverX0 and
/// .arealDensityGPerCm2 are read, activating PID-aware energy loss/
/// straggling unlike the legacy MCS-only attachHit<DiskDisk> above --
/// predictedChi2, reject only if predictedChi2 > params.maxChi2ClusterAttachment
/// (no `< 0.f` guard, matching the legacy attachHit<DiskDisk> exactly), then
/// update. Defined out of line in TransitionPolicyOperations.cxx.
template <>
bool attachHit<TransitionPolicyTag::DiskDisk>(
  SurfaceKinematicState& state,
  const o2::its::TrackingFrameInfo& hit,
  const NominalSurfaceMaterial& material,
  float bz,
  float& chi2,
  const DiskDiskPolicyParams& params,
  OperationFailureReason& reason) noexcept;

/// Stage-B Slice B (design report Sec 8/11): native SurfaceKinematicState
/// overload of the legacy cellsAreCompatible<Tag> at the top of this file,
/// distinguished by argument type/arity (current/next here are
/// SurfaceKinematicState, never CellSeedTpl<...>, and this overload takes
/// two extra explicit cluster-index parameters). Additive and unwired.
///
/// General contract (both families): `current` is the inner Cell state,
/// `next` is the outer Cell state; `next` is copied into local scratch,
/// transported inward to current's reference surface/frame, and evaluated
/// against untouched `current` -- neither input is ever mutated (enforced by
/// the type system: both are `const&`). Threshold comparison is inclusive
/// (`<=`), exactly as today.
///
/// `currentSecondClusterIndex`/`nextFirstClusterIndex`: the temporary,
/// explicitly threaded raw legacy cluster-index continuity input
/// (requirement: this native state-only API cannot silently lose
/// detail::mftFwdCellsAreCompatible's `currentCell.getSecondClusterIndex()
/// == nextCell.getFirstClusterIndex()` check). CylinderCylinder ignores both
/// (mirrors passesCellRoadPrecut<CylinderCylinder>'s own unused-parameter
/// pattern above); DiskDisk checks them first, exactly matching legacy
/// precedence. This raw-index input is temporary and will be replaced by a
/// topology surface plus SurfaceMeasurementIndex in the measurement-index
/// slice; continuity is never inferred from state fields.
///
/// No `reason` out-param, unlike native buildCellSeed/attachHit above: this
/// mirrors the legacy cellsAreCompatible<Tag>'s own bool-only contract
/// exactly (it has none today either); rotation/propagation/family-mismatch
/// causes are already covered by the focused barrel::/forward:: primitive
/// tests and simply surface here as `false`.
///
/// Host-only for this slice (requirement: no GPU/device claim). This
/// deliberately differs from the legacy cellsAreCompatible<Tag> above, which
/// is declared unconditionally (outside any GPUCA_GPUCODE guard): that
/// declaration's host-only status is carried entirely by its .cxx
/// definition, with prior evidence of an eventual device-portable
/// dependency chain. This new operation has no such evidence yet, so it is
/// declared behind this guard instead of inheriting the unconditional
/// declaration.
template <TransitionPolicyTag Tag>
bool cellsAreCompatible(const SurfaceKinematicState& current,
                        const SurfaceKinematicState& next,
                        int currentSecondClusterIndex,
                        int nextFirstClusterIndex,
                        float bz,
                        const typename TransitionPolicyTraits<Tag>::Params& params) noexcept;

/// Barrel formula: copy `next` into scratch, rotate to current.alpha,
/// propagate to current.referenceCoordinate, then barrel::stateChi2(current,
/// scratch, ...) <= params.maxChi2ClusterAttachment -- matching the legacy
/// cellsAreCompatible<CylinderCylinder> formula exactly. Cluster-index
/// parameters are ignored (CylinderCylinder has no continuity check, legacy
/// or native). Defined out of line in TransitionPolicyOperations.cxx.
template <>
bool cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(
  const SurfaceKinematicState& current,
  const SurfaceKinematicState& next,
  int currentSecondClusterIndex,
  int nextFirstClusterIndex,
  float bz,
  const CylinderCylinderPolicyParams& params) noexcept;

/// Disk formula: checked first, `currentSecondClusterIndex ==
/// nextFirstClusterIndex` (see the primary template's doc); then copy `next`
/// into scratch, propagate to current.referenceCoordinate with the accepted
/// forward model, then forward::stateChi2(current, scratch, ...) <=
/// params.maxChi2ClusterAttachment -- matching
/// detail::mftFwdCellsAreCompatible's formula and precedence exactly.
/// Defined out of line in TransitionPolicyOperations.cxx.
template <>
bool cellsAreCompatible<TransitionPolicyTag::DiskDisk>(
  const SurfaceKinematicState& current,
  const SurfaceKinematicState& next,
  int currentSecondClusterIndex,
  int nextFirstClusterIndex,
  float bz,
  const DiskDiskPolicyParams& params) noexcept;

/// D007 policy-boundary operation (Architecture.md Sec 10, Gate 3 cell-road
/// pre-cut slice): true if the three-cluster candidate {inner, middle, outer}
/// should be considered by buildCellSeed<Tag> at all. This is the last
/// detector-branch removed from TrackerTraits::computeLayerCellsForPolicy's
/// candidate loop -- the caller now issues one unconditional call for both
/// families, exactly like cellsAreCompatible/attachHit/buildCellSeed.
///
/// Operates on `GlobalPoint3F` (SurfaceMeasurement.h), the existing
/// detector-neutral 3-float global-position primitive, rather than
/// `o2::its::Cluster`: this header therefore does not include
/// ITStracking/Cluster.h, MFTFwdTrackHelpers.h, or any MFT/detector geometry
/// header. The caller converts its three Cluster positions to GlobalPoint3F
/// immediately before the call (TrackerTraits.cxx).
///
/// Defined inline, in this header, unlike cellsAreCompatible/attachHit/
/// buildCellSeed: this predicate calls nothing host-only external (no
/// Propagator singleton, no TrackParCovFwd Kalman machinery), so there is no
/// reason to hide an implementation behind a separate translation unit. It
/// runs once per CA candidate, before the expensive fit it guards, so being
/// visible for ordinary same-translation-unit inlining in the caller
/// (TrackerTraits.cxx already includes this header) matters -- this is
/// described as inlineable/same-TU, not as a guarantee that any particular
/// compiler/flags combination will inline or eliminate every construction.
///
/// `layerInner`/`layerMiddle`/`layerOuter` are legacy layout-local layer
/// indices (the same index space as `TrackletProjectionState<Tag>` and
/// `LayerScatteringInputs<Tag>::referenceCoordinate` already document).
/// `perLayerReferenceZ` is the caller's once-per-iteration cached span from
/// `DiskDiskReferenceCoordinateView`/`bindLegacyMFTReferenceCoordinates()` --
/// this operation does not call `mftLayerZ()`, `MFT::LayerZCoordinate()`, or
/// any detector constant lookup itself. Unchecked caller precondition
/// (matching every other operation in this file): `perLayerReferenceZ` must
/// have a valid element at each of `layerInner`, `layerMiddle`, `layerOuter`
/// for the `DiskDisk` specialization; the `CylinderCylinder` specialization
/// never reads any of its layer/span/cluster-position arguments.
///
/// Pure predicate: no argument is mutated, on acceptance or rejection.
/// noexcept: the arithmetic is total for valid inputs (see the DiskDisk
/// specialization doc for the exact, non-blanket NaN/Inf behavior); this
/// slice adds no new finite/range/physics validation beyond the legacy
/// formula's own behavior.
template <TransitionPolicyTag Tag>
bool passesCellRoadPrecut(const GlobalPoint3F& pointInner,
                          const GlobalPoint3F& pointMiddle,
                          const GlobalPoint3F& pointOuter,
                          int layerInner, int layerMiddle, int layerOuter,
                          gsl::span<const float> perLayerReferenceZ,
                          const typename TransitionPolicyTraits<Tag>::Params& params) noexcept;

/// CylinderCylinder has no geometric road pre-cut; compile-time no-op so the
/// caller's candidate loop stays a single unconditional call for both
/// families (Architecture.md 10.1: dispatch is a compile-time tag switch to
/// statically compiled functions, never a runtime branch in the hot loop).
/// Reads none of its arguments -- including an empty/default `perLayerReferenceZ`,
/// which is always safe here.
template <>
inline bool passesCellRoadPrecut<TransitionPolicyTag::CylinderCylinder>(
  const GlobalPoint3F&, const GlobalPoint3F&, const GlobalPoint3F&,
  int, int, int, gsl::span<const float>, const CylinderCylinderPolicyParams&) noexcept
{
  return true;
}

/// Disk formula: the legacy detail::validateMFTCellClusters() three-check
/// road pre-cut (CellRoadRCut / ROADclsRCut), re-expressed on GlobalPoint3F
/// and an explicit per-layer reference z instead of o2::its::Cluster and an
/// internal mftLayerZ() lookup. Cluster ordering is strictly {inner, middle,
/// outer}; every distance-function argument permutation and the asymmetric
/// conical-scale argument order are preserved exactly, as is strict `<` for
/// all three checks and short-circuit evaluation of the combining `&&`.
/// `params.cellRoadRCut` is the unsquared cut (matching
/// DiskDiskPolicyParams::cellRoadRCut's existing contract); it is squared
/// exactly once, here.
template <>
inline bool passesCellRoadPrecut<TransitionPolicyTag::DiskDisk>(
  const GlobalPoint3F& pointInner, const GlobalPoint3F& pointMiddle, const GlobalPoint3F& pointOuter,
  int layerInner, int layerMiddle, int layerOuter,
  gsl::span<const float> perLayerReferenceZ,
  const DiskDiskPolicyParams& params) noexcept
{
  // Squared transverse distance from `point` to the seed line from->to
  // (legacy detail::mftDistanceToSeedSquared / MFT getDistanceToSeed, moved
  // here verbatim -- this arithmetic never referenced MFT constants, only
  // measured coordinates). |dzSeed| < 1e-9f (near-degenerate seed line,
  // parallel to the measurement plane) returns FLT_MAX, guaranteeing
  // rejection against any finite cut.
  auto distanceToSeedLineSquared = [](const GlobalPoint3F& from, const GlobalPoint3F& to, const GlobalPoint3F& point) -> float {
    const float dxSeed = to.x - from.x;
    const float dySeed = to.y - from.y;
    const float dzSeed = to.z - from.z;
    if (std::abs(dzSeed) < 1e-9f) {
      return std::numeric_limits<float>::max();
    }
    const float invdzSeed = (point.z - from.z) / dzSeed;
    const float xSeed = from.x + dxSeed * invdzSeed;
    const float ySeed = from.y + dySeed * invdzSeed;
    const float dx = point.x - xSeed;
    const float dy = point.y - ySeed;
    return dx * dx + dy * dy;
  };

  // Conical road scale (legacy detail::mftConicalRoadR2Scale), taking zFrom/
  // zTo explicitly instead of calling mftLayerZ(layerFrom/layerTo)
  // internally. |zFrom| < 1e-6f (near-zero reference z) falls back to 1.f
  // (no conical correction). Preserved exactly as `1.f + (zTo - zFrom) /
  // zFrom` -- do NOT simplify to zTo/zFrom, which would change float
  // rounding.
  auto conicalRoadR2Scale = [](float zFrom, float zTo) -> float {
    if (std::abs(zFrom) < 1e-6f) {
      return 1.f;
    }
    const float dCone = 1.f + (zTo - zFrom) / zFrom;
    return dCone * dCone;
  };

  const float zInner = perLayerReferenceZ[layerInner];
  const float zMiddle = perLayerReferenceZ[layerMiddle];
  const float zOuter = perLayerReferenceZ[layerOuter];
  const float r2Cut = params.cellRoadRCut * params.cellRoadRCut;

  return distanceToSeedLineSquared(pointInner, pointOuter, pointMiddle) < r2Cut * conicalRoadR2Scale(zInner, zMiddle) &&
         distanceToSeedLineSquared(pointInner, pointMiddle, pointOuter) < r2Cut * conicalRoadR2Scale(zInner, zOuter) &&
         distanceToSeedLineSquared(pointMiddle, pointOuter, pointInner) < r2Cut * conicalRoadR2Scale(zMiddle, zInner);
}

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

/// Unchecked caller preconditions (noexcept reflects the arithmetic being
/// total for valid inputs, not that these are validated here -- matching
/// every other operation in this file, which likewise do not re-validate
/// caller-supplied indices/ranges):
///  - `perLayerMSAngle` must have a valid element at every index in
///    `[fromLayer, toLayer)`; it is indexed without a bounds check.
///  - `fromLayer` and `toLayer` must be validated legacy layout-local layer
///    indices for the active topology (the same index space
///    TrackletProjectionState<Tag>/TrackerTraits already document and
///    validate elsewhere via TrackingTopology/DetectorLayoutView), not
///    arbitrary integers or global SurfaceIds.
///  - `fromLayer <= toLayer`; the accumulation loop is a no-op (`ms2 == 0`)
///    for `fromLayer == toLayer` and undefined for `fromLayer > toLayer`
///    only in the sense of iterating zero times, but callers must not rely
///    on that as a substitute for passing a valid transition's own
///    fromLayer/toLayer pair.
/// This slice adds no new runtime validation for these and does not change
/// any legacy behavior; violating them is a caller bug, exactly as it would
/// have been in the inline code this operation replaces.
TransitionScatteringBendingPrep prepareTransitionScatteringAndBending(
  gsl::span<const float> perLayerMSAngle, int fromLayer, int toLayer,
  float r1, float r2, float clampedOneOverR, float res1, float res2) noexcept;

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRANSITIONPOLICYOPERATIONS_H_ */
