// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_NATIVECYLINDERCYLINDERREFITDRIVER_H_
#define ALICEO2_ITSMFT_TRACKING_NATIVECYLINDERCYLINDERREFITDRIVER_H_

#include "GPUCommonDef.h"

#ifndef GPUCA_GPUCODE

#include <array>
#include <cmath>

#include <gsl/span>

#include "CommonConstants/MathConstants.h"
#include "DataFormatsITS/TrackITS.h"
#include "ITSMFTTracking/RefitLegAssembly.h"
#include "ITSMFTTracking/SurfaceCatalogView.h"
#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"
#include "ITSMFTTracking/TransitionPolicyOperations.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

namespace o2::itsmft::tracking
{

/// Gate 3 Slice B prerequisite: native equivalent of the frozen ITS
/// `resetTrackCovariance(o2::track::TrackParCov&)` (ITSMFTTracking/
/// TrackHelpers.h) for a Barrel-family `SurfaceKinematicState`. Applied by
/// the frozen `refitTrack` driver once per leg (immediately before each of
/// its three `fitTrack` calls) so every leg's Kalman filter starts from a
/// loose, uninformative covariance while keeping the previous leg's fitted
/// parameters -- not merely a one-time seed-preparation step. Exact port:
/// `resetCovariance(0.f)` zeroes every element then sets the diagonal to
/// `{kCY2max, kCZ2max, kCSnp2max, kCTgl2max, kC1Pt2max}` (o2::track::
/// TrackParametrizationWithError<float>::resetCovariance, DataFormats/
/// Reconstruction/src/TrackParametrizationWithError.cxx); the trailing
/// `setCov(getQ2Pt()*getQ2Pt()*getCov()[kSigQ2Pt2], kSigQ2Pt2)` then
/// rescales only the freshly-reset Q2Pt diagonal entry by the state's own
/// q/pT squared. Reads/writes only `state.covariance`; `state.parameters`
/// (including q/pT) is never modified.
///
/// Unwired: no production call site uses this in this slice.
GPUhdi() void resetCylinderCylinderCovarianceForRefit(SurfaceKinematicState& state) noexcept
{
  for (auto& element : state.covariance) {
    element = 0.f;
  }
  state.covariance[packedCovarianceIndex(0, 0)] = o2::track::kCY2max;
  state.covariance[packedCovarianceIndex(1, 1)] = o2::track::kCZ2max;
  state.covariance[packedCovarianceIndex(2, 2)] = o2::track::kCSnp2max;
  state.covariance[packedCovarianceIndex(3, 3)] = o2::track::kCTgl2max;
  const float q2pt = state.parameters[4];
  state.covariance[packedCovarianceIndex(4, 4)] = q2pt * q2pt * o2::track::kC1Pt2max;
}

/// Native equivalent of `o2::track::TrackParametrization<float>::getPt()`
/// (DataFormats/Reconstruction/include/ReconstructionDataFormats/
/// TrackParametrization.h), for the frozen `refitTrackSeed`'s trailing
/// `if (minPt > 0.f && track.getPt() < minPt) return false;` MinPt check.
/// Exact port of `getPtInv()`'s `absCharge > 1` scaling and `MinPTInv`
/// floor (avoids a division by an arbitrarily small q/pT).
GPUhdi() float getPtFromQOverPt(float q2pt, uint8_t absCharge) noexcept
{
  float ptInv = std::abs(q2pt);
  if (ptInv < o2::track::MinPTInv) {
    ptInv = o2::track::MinPTInv;
  }
  if (absCharge > 1) {
    ptInv /= static_cast<float>(absCharge);
  }
  return 1.f / ptInv;
}

/// Gate 3 Slice B: unwired, native `SurfaceKinematicState`-based refit
/// driver reproducing the frozen ITS `refitTrack`/`refitTrackSeed` leg
/// sequencing (ITSMFTTracking/TrackHelpers.h) on `driveRefitLeg<
/// TransitionPolicyTag::CylinderCylinder>` instead of the frozen per-hit
/// `fitTrack` loop. No production call site uses this in this slice; it is
/// not called from DetectorTraits::refitSeed.
///
/// Leg sequencing, exactly matching `refitTrack`:
///   Leg A (inward-index, `[0, NLayers)`, step +1): seeded from `seed.
///     state()` as-is (see the stated scope limitation below), direction
///     AlongMomentum, maxQoverPt = VeryBig. Its result becomes `outParamOut`
///     unless `repeatRefitOut` succeeds below.
///   Leg B (outward-index, `[NLayers-1, -1)`, step -1): re-seeded from leg
///     A's result (fresh `SurfaceLinearizationReference`, fresh covariance
///     via `resetCylinderCylinderCovarianceForRefit`), direction
///     OppositeMomentum, maxQoverPt = 50.f. Its result and chi2 are
///     *always* the reported `outParamIn`/`outChi2`, regardless of
///     `repeatRefitOut` -- matching `refitTrack`'s own `track.paramIn`/
///     `track.setChi2(saveChi2)` restoration after an optional leg C.
///   Leg C (optional, only if `repeatRefitOut`): re-seeded from leg B's
///     result, same direction/maxQoverPt as leg A. If it fails, the *whole*
///     driver fails (matching `refitTrack`'s `if (!fitSuccess) return
///     false;` inside its `repeatRefitOut` block -- leg B's already-fitted
///     result is not returned as a fallback). If it succeeds, its result
///     becomes `outParamOut` (replacing leg A's).
/// Each leg's own final acceptance (`|Q2Pt| < maxQoverPt && chi2 <
/// maxChi2NDF*(acceptedHitCount*2-5)`) is re-applied here after
/// `driveRefitLeg` itself succeeds, since `driveRefitLeg`'s own contract
/// covers only the per-hit chi2 gate, not this whole-leg formula; failure
/// here fails the whole driver with `OperationFailureReason::
/// LegAcceptanceFailure`, matching `fitTrack`'s own trailing return
/// condition exactly (same VeryBig/50.f thresholds per leg).
///
/// After leg B, the frozen `refitTrackSeed`'s MinPt check is reproduced
/// exactly: `params.MinPt[NLayers - seed.getHitLayerMask().count()]`
/// (`LayerMask::count()` is the same population count `TrackITSInternal::
/// getNClusters()` tracks) compared against `getPtFromQOverPt` on leg B's
/// result; failure reports `OperationFailureReason::MinPtFailure`.
///
/// `shiftReferenceToMeasurement`/`maxChi2ClusterAttachment` are passed to
/// every `driveRefitLeg` call unchanged, matching `ctx.shiftRefToCluster`/
/// `ctx.maxChi2ClusterAttachment`'s own per-call-identical usage in the
/// frozen driver.
///
/// Stated scope limitation (Gate 3 Slice B, not implemented here): the
/// frozen `seedTrackForRefit`'s conditional mid-track geometric reseed
/// (`ncl < reseedIfShorter && ncl > 2`, re-deriving the initial
/// parametrization via `buildTrackSeed`/`selectReseedMidLayer` from raw
/// `Cluster`/`TrackingFrameInfo`) is *not* reproduced by this driver. This
/// driver always starts leg A from `seed.state()` unchanged -- callers
/// (and this slice's own tests) must not rely on this driver for seeds
/// whose cluster count would trigger that legacy reseed path. Reproducing
/// it natively is deferred to a follow-up slice.
///
/// Host-only (no GPU/device-readiness claim), like driveRefitLeg/refitHit
/// above, whose contracts this driver inherits unchanged.
template <int NLayers>
bool nativeRefitTrackCylinderCylinder(
  const TrackSeedN<NLayers>& seed,
  const LayerMeasurementSpans<NLayers>& layerMeasurements,
  SurfaceCatalogView surfaceCatalog,
  float bz,
  bool shiftReferenceToMeasurement,
  float maxChi2ClusterAttachment,
  float maxChi2NDF,
  bool repeatRefitOut,
  gsl::span<const float> minPt,
  SurfaceKinematicState& outParamIn,
  SurfaceKinematicState& outParamOut,
  float& outChi2,
  OperationFailureReason& reason) noexcept
{
  constexpr auto Tag = TransitionPolicyTag::CylinderCylinder;

  auto legAcceptable = [](const SurfaceKinematicState& state, float chi2, uint32_t acceptedHitCount,
                          float maxQoverPt, float maxChi2NDFValue) noexcept -> bool {
    if (!(std::abs(state.parameters[4]) < maxQoverPt)) {
      return false;
    }
    return chi2 < maxChi2NDFValue * static_cast<float>(static_cast<int>(acceptedHitCount) * 2 - 5);
  };

  // --- Leg A: inward-index, seeded from the seed's own state as-is. ---
  SurfaceKinematicState stateA = seed.state();
  SurfaceLinearizationReference linRefA{};
  if (!makeLinearizationReference(stateA, linRefA)) {
    reason = OperationFailureReason::SourceFamilyMismatch;
    return false;
  }
  resetCylinderCylinderCovarianceForRefit(stateA);
  float chi2A = 0.f;
  uint32_t acceptedA = 0;
  std::array<SurfaceMeasurement, NLayers> slotsBufferA{};
  const auto slotsA = assembleRefitLegSlots<NLayers>(seed, layerMeasurements, 0, NLayers, 1, slotsBufferA);
  if (!driveRefitLeg<Tag>(stateA, linRefA, chi2A, acceptedA, slotsA, surfaceCatalog, bz,
                          material::MaterialTraversalDirection::AlongMomentum, shiftReferenceToMeasurement,
                          maxChi2ClusterAttachment, reason)) {
    return false;
  }
  if (!legAcceptable(stateA, chi2A, acceptedA, o2::constants::math::VeryBig, maxChi2NDF)) {
    reason = OperationFailureReason::LegAcceptanceFailure;
    return false;
  }

  // --- Leg B: outward-index, re-seeded from leg A's result. Its result and
  // chi2 are always the reported outParamIn/outChi2. ---
  SurfaceKinematicState stateB = stateA;
  SurfaceLinearizationReference linRefB{};
  if (!makeLinearizationReference(stateB, linRefB)) {
    reason = OperationFailureReason::SourceFamilyMismatch;
    return false;
  }
  resetCylinderCylinderCovarianceForRefit(stateB);
  float chi2B = 0.f;
  uint32_t acceptedB = 0;
  std::array<SurfaceMeasurement, NLayers> slotsBufferB{};
  const auto slotsB = assembleRefitLegSlots<NLayers>(seed, layerMeasurements, NLayers - 1, -1, -1, slotsBufferB);
  if (!driveRefitLeg<Tag>(stateB, linRefB, chi2B, acceptedB, slotsB, surfaceCatalog, bz,
                          material::MaterialTraversalDirection::OppositeMomentum, shiftReferenceToMeasurement,
                          maxChi2ClusterAttachment, reason)) {
    return false;
  }
  if (!legAcceptable(stateB, chi2B, acceptedB, 50.f, maxChi2NDF)) {
    reason = OperationFailureReason::LegAcceptanceFailure;
    return false;
  }

  // --- MinPt check (frozen refitTrackSeed), keyed on the seed's own
  // attached-cluster count -- unaffected by which leg produced stateB. ---
  const int nClAttached = seed.getHitLayerMask().count();
  const int minPtSlot = NLayers - nClAttached;
  if (minPtSlot >= 0 && minPtSlot < static_cast<int>(minPt.size())) {
    const float minPtThreshold = minPt[minPtSlot];
    if (minPtThreshold > 0.f && getPtFromQOverPt(stateB.parameters[4], stateB.absCharge) < minPtThreshold) {
      reason = OperationFailureReason::MinPtFailure;
      return false;
    }
  }

  // --- Optional leg C: re-seeded from leg B's result, inward-index again.
  // Failing leg C fails the whole driver (matching refitTrack exactly);
  // succeeding leg C's result replaces leg A's as outParamOut. ---
  SurfaceKinematicState stateOut = stateA;
  if (repeatRefitOut) {
    SurfaceKinematicState stateC = stateB;
    SurfaceLinearizationReference linRefC{};
    if (!makeLinearizationReference(stateC, linRefC)) {
      reason = OperationFailureReason::SourceFamilyMismatch;
      return false;
    }
    resetCylinderCylinderCovarianceForRefit(stateC);
    float chi2C = 0.f;
    uint32_t acceptedC = 0;
    std::array<SurfaceMeasurement, NLayers> slotsBufferC{};
    const auto slotsC = assembleRefitLegSlots<NLayers>(seed, layerMeasurements, 0, NLayers, 1, slotsBufferC);
    if (!driveRefitLeg<Tag>(stateC, linRefC, chi2C, acceptedC, slotsC, surfaceCatalog, bz,
                            material::MaterialTraversalDirection::AlongMomentum, shiftReferenceToMeasurement,
                            maxChi2ClusterAttachment, reason)) {
      return false;
    }
    if (!legAcceptable(stateC, chi2C, acceptedC, o2::constants::math::VeryBig, maxChi2NDF)) {
      reason = OperationFailureReason::LegAcceptanceFailure;
      return false;
    }
    stateOut = stateC;
  }

  outParamIn = stateB;
  outParamOut = stateOut;
  outChi2 = chi2B;
  return true;
}

/// Gate 3 Slice B: converts a successful `nativeRefitTrackCylinderCylinder`
/// result into a `TrackITSExt`, field for field matching `makeTrackITSExt`
/// (ITStracking/TrackITSInternal.h) exactly: `paramIn`/`paramOut` via the
/// existing, already-generic `legacy::exportBarrelTrackParCov` (reused here
/// a second time -- the same function `DetectorTraits.cxx`'s frozen ITS
/// branch already calls once at entry), `chi2` passthrough, `timeStamp`
/// from `seed.getTimeStamp().makeSymmetrical()`, and per-layer external
/// cluster index/pattern passthrough from `seed.getCluster(layer)` --
/// refit never touches hit assignment, only kinematics/chi2. Returns false
/// (leaving `outTrack` untouched) if either export fails, which can only
/// happen if `paramIn`/`paramOut` are not `StateFamily::Barrel` -- a
/// caller/precondition violation, since every state this driver produces is
/// Barrel by construction.
///
/// Unwired: no production call site uses this in this slice.
template <int NLayers>
bool exportNativeRefitToTrackITSExt(
  const TrackSeedN<NLayers>& seed,
  const SurfaceKinematicState& paramIn,
  const SurfaceKinematicState& paramOut,
  float chi2,
  o2::its::TrackITSExt& outTrack) noexcept
{
  o2::track::TrackParCovF legacyParamIn{};
  o2::track::TrackParCovF legacyParamOut{};
  if (!legacy::exportBarrelTrackParCov(paramIn, legacyParamIn) ||
      !legacy::exportBarrelTrackParCov(paramOut, legacyParamOut)) {
    return false;
  }

  o2::its::TrackITSExt scratch{};
  scratch.getParamIn() = legacyParamIn;
  scratch.getParamOut() = legacyParamOut;
  scratch.setChi2(chi2);
  scratch.getTimeStamp() = seed.getTimeStamp().makeSymmetrical();
  for (int layer = 0; layer < NLayers; ++layer) {
    const int clsIdx = seed.getCluster(layer);
    if (clsIdx != o2::its::constants::UnusedIndex) {
      scratch.setExternalClusterIndex(layer, clsIdx, true);
    }
  }
  outTrack = scratch;
  return true;
}

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_NATIVECYLINDERCYLINDERREFITDRIVER_H_ */
