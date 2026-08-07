// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_NATIVEREFITDRIVER_H_
#define ALICEO2_ITSMFT_TRACKING_NATIVEREFITDRIVER_H_

#include "GPUCommonDef.h"

#ifndef GPUCA_GPUCODE

#include <cmath>
#include <vector>

#include <gsl/span>

#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/RefitLegAssembly.h"
#include "ITSMFTTracking/SurfaceCatalogView.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

// M5d: the shared (cylinder- and disk-common) whole-seed refit driver, built
// entirely on Propagator (Propagator.h). This is the "one tracklet/cell/
// road/refit flow" ADR 0007 decision 10 and the GenericTrackingEngineMigration
// M5 plan describe for the refit stage: fitTrackSeedLegs below
// contains no family/Tag branch of its own -- every family difference is
// already confined inside Propagator::propagateToMeasurement's own
// descriptor-driven dispatch. It is the single native refit entry point used by
// the adapter refit operation.
namespace o2::itsmft::tracking
{

// Reproduces the "loose, uninformative diagonal" covariance reset every leg
// of a from-scratch refit starts from (fresh Kalman filter, previous leg's
// fitted parameters retained), family-dispatched by state.family rather than
// by Tag. Barrel reuses the exact ceiling constants
// the established barrel covariance formula (bit-for-bit identical, just resolved at
// runtime instead of compiled once for CylinderCylinder). Forward has no
// legacy analogue to port (the frozen MFT TrackFitter/TrackLTF Kalman engine
// this milestone removes from production never exposed a comparable
// "reset to uninformative" primitive on its own parametrization) -- these are
// new native ceiling constants, chosen by the same rationale as their Barrel
// counterparts: a diagonal loose enough to cover this parameter's full
// physical range. X/Y position and InvQPt/Tanl reuse the Barrel position/
// curvature/slope ceilings directly (same physical quantity, same units);
// Phi is a full angle (unlike Barrel's bounded Snp=sin(local angle)), so its
// own ceiling is (pi)^2 -- "SigmaPhi <= pi", the same "cover the full range"
// rationale kCSnp2max's own comment states for Snp's [-1, 1] range.
GPUhdi() void resetCovarianceForRefit(SurfaceKinematicState& state) noexcept
{
  for (auto& element : state.covariance) {
    element = 0.f;
  }
  if (state.family == StateFamily::Barrel) {
    state.covariance[packedCovarianceIndex(0, 0)] = o2::track::kCY2max;
    state.covariance[packedCovarianceIndex(1, 1)] = o2::track::kCZ2max;
    state.covariance[packedCovarianceIndex(2, 2)] = o2::track::kCSnp2max;
    state.covariance[packedCovarianceIndex(3, 3)] = o2::track::kCTgl2max;
    const float q2pt = state.parameters[4];
    state.covariance[packedCovarianceIndex(4, 4)] = q2pt * q2pt * o2::track::kC1Pt2max;
  } else {
    constexpr float kCPhi2maxForward = o2::constants::math::PI * o2::constants::math::PI;
    state.covariance[packedCovarianceIndex(0, 0)] = o2::track::kCY2max;
    state.covariance[packedCovarianceIndex(1, 1)] = o2::track::kCY2max;
    state.covariance[packedCovarianceIndex(2, 2)] = kCPhi2maxForward;
    state.covariance[packedCovarianceIndex(3, 3)] = o2::track::kCTgl2max;
    const float invQPt = state.parameters[4];
    state.covariance[packedCovarianceIndex(4, 4)] = invQPt * invQPt * o2::track::kC1Pt2max;
  }
}

// Same physical formula as the native minimum-pT helper: parameters[4] is the
// raw signed q/pT for both families
// (SurfaceKinematicState.h / ForwardStateView::getQ2Pt() doc), so no
// family branch is needed here either.
GPUhdi() float ptFromQOverPt(float q2pt, uint8_t absCharge) noexcept
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

// Shared (family-blind) whole-seed refit: reproduces the three-leg
// (inward-index A, outward-index B, optional inward-index C) sequencing and
// per-leg acceptance-gate structure nativeRefitTrackCylinderCylinder
// documents in full in the native operation -- same leg
// direction/maxQoverPt/acceptance formula, same MinPt check keyed on the
// seed's own attached-cluster count -- but built on Propagator::driveRefitLeg
// (descriptor-driven) instead of driveRefitLeg<Tag>, so the identical
// function body serves a Barrel seed (ITS) or a Forward seed (MFT) with no
// specialization. `reseedIfShorter` is not carried over: this is a new
// production entry point, not a byte-for-byte port, and no caller of this
// function threads a nonzero value through it (see the adapter refit helper).
//
// Transactional exactly like nativeRefitTrackCylinderCylinder: outParamIn/
// outParamOut/outChi2 are committed only on complete success.
inline bool fitTrackSeedLegs(
  const TrackSeed& seed,
  gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
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
  resetCovarianceForRefit(stateA);
  float chi2A = 0.f;
  uint32_t acceptedA = 0;
  const int activeSurfaceCount = static_cast<int>(layerMeasurements.size());
  std::vector<SurfaceMeasurement> slotsBufferA(static_cast<std::size_t>(activeSurfaceCount));
  const auto slotsA = assembleRefitLegSlots(seed, layerMeasurements, 0, activeSurfaceCount, 1, slotsBufferA);
  if (!Propagator::driveRefitLeg(stateA, linRefA, chi2A, acceptedA, slotsA, surfaceCatalog, bz,
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
  resetCovarianceForRefit(stateB);
  float chi2B = 0.f;
  uint32_t acceptedB = 0;
  std::vector<SurfaceMeasurement> slotsBufferB(static_cast<std::size_t>(activeSurfaceCount));
  const auto slotsB = assembleRefitLegSlots(seed, layerMeasurements, activeSurfaceCount - 1, -1, -1, slotsBufferB);
  if (!Propagator::driveRefitLeg(stateB, linRefB, chi2B, acceptedB, slotsB, surfaceCatalog, bz,
                                 material::MaterialTraversalDirection::OppositeMomentum, shiftReferenceToMeasurement,
                                 maxChi2ClusterAttachment, reason)) {
    return false;
  }
  if (!legAcceptable(stateB, chi2B, acceptedB, 50.f, maxChi2NDF)) {
    reason = OperationFailureReason::LegAcceptanceFailure;
    return false;
  }

  // --- MinPt check, keyed on the seed's own attached-cluster count. ---
  const int nClAttached = seed.getHitLayerMask().count();
  const int minPtSlot = activeSurfaceCount - nClAttached;
  if (minPtSlot >= 0 && minPtSlot < static_cast<int>(minPt.size())) {
    const float minPtThreshold = minPt[minPtSlot];
    if (minPtThreshold > 0.f && ptFromQOverPt(stateB.parameters[4], stateB.absCharge) < minPtThreshold) {
      reason = OperationFailureReason::MinPtFailure;
      return false;
    }
  }

  // --- Optional leg C: re-seeded from leg B's result, inward-index again. ---
  SurfaceKinematicState stateOut = stateA;
  if (repeatRefitOut) {
    SurfaceKinematicState stateC = stateB;
    SurfaceLinearizationReference linRefC{};
    if (!makeLinearizationReference(stateC, linRefC)) {
      reason = OperationFailureReason::SourceFamilyMismatch;
      return false;
    }
    resetCovarianceForRefit(stateC);
    float chi2C = 0.f;
    uint32_t acceptedC = 0;
    std::vector<SurfaceMeasurement> slotsBufferC(static_cast<std::size_t>(activeSurfaceCount));
    const auto slotsC = assembleRefitLegSlots(seed, layerMeasurements, 0, activeSurfaceCount, 1, slotsBufferC);
    if (!Propagator::driveRefitLeg(stateC, linRefC, chi2C, acceptedC, slotsC, surfaceCatalog, bz,
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

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_NATIVEREFITDRIVER_H_ */
