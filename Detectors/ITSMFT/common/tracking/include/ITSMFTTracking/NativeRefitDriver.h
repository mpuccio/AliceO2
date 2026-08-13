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
#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/GlobalMeasurement.h"
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceStateOperationResult.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

// Shared cylinder/disk refit entry point built on descriptor-driven
// Propagator operations; cylinder/disk dispatch remains inside propagation.
namespace o2::itsmft::tracking
{

/// Builds one traversal-ordered leg of local measurement slots for native
/// refit from a `TrackSeed`'s already-attached, layer-indexed cluster
/// bookkeeping. `[start, end)` stepping by `step` is the caller-supplied
/// layer-index range; `step == +1` walks inward and `step == -1` walks
/// outward. Holes have `present == false`,
/// and valid cluster indices are looked up in the corresponding layer span.
/// Slots are written in traversal order, so decreasing legs remain reversed
/// rather than being silently reordered by source layer index.
inline gsl::span<const RefitMeasurementSlot> assembleRefitLegSlots(
  const TrackSeed& seed,
  gsl::span<const gsl::span<const GlobalMeasurement>> layerGlobals,
  gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
  int start, int end, int step,
  gsl::span<RefitMeasurementSlot> out) noexcept
{
  int position = 0;
  for (int surfacePosition = start; surfacePosition != end && position < static_cast<int>(out.size()); surfacePosition += step) {
    const int clsIdx = seed.getCluster(surfacePosition);
    if (clsIdx == o2::its::constants::UnusedIndex) {
      out[position++] = {};
      continue;
    }
    out[position++] = RefitMeasurementSlot{layerMeasurements[surfacePosition][clsIdx],
                                           layerGlobals[surfacePosition][clsIdx].surface, true};
  }
  return gsl::span<const RefitMeasurementSlot>(out.data(), position);
}

// Reset each refit leg to a loose diagonal. Barrel uses the existing ceilings;
// Forward uses the same position/slope/curvature units and (pi)^2 for Phi.
GPUhdi() void resetCovarianceForRefit(SurfaceKinematicState& state) noexcept
{
  for (auto& element : state.covariance) {
    element = 0.f;
  }
  if (state.kind == SurfaceKind::Cylinder) {
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
// raw signed q/pT for both coordinate conventions
// (SurfaceKinematicState.h / ForwardStateView::getQ2Pt() doc), so no
// SurfaceKind branch is needed here either.
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

// Shared three-leg refit: inward A, outward B, and optional inward C. The
// acceptance formula and MinPt check use the seed's attached-cluster count.
// Outputs are committed only on complete success.
inline bool fitTrackSeedLegs(
  const TrackSeed& seed,
  gsl::span<const gsl::span<const GlobalMeasurement>> layerGlobals,
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
    reason = OperationFailureReason::SourceSurfaceKindMismatch;
    return false;
  }
  resetCovarianceForRefit(stateA);
  float chi2A = 0.f;
  uint32_t acceptedA = 0;
  const int activeSurfaceCount = static_cast<int>(layerMeasurements.size());
  std::vector<RefitMeasurementSlot> slotsBufferA(static_cast<std::size_t>(activeSurfaceCount));
  const auto slotsA = assembleRefitLegSlots(seed, layerGlobals, layerMeasurements, 0, activeSurfaceCount, 1, slotsBufferA);
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
    reason = OperationFailureReason::SourceSurfaceKindMismatch;
    return false;
  }
  resetCovarianceForRefit(stateB);
  float chi2B = 0.f;
  uint32_t acceptedB = 0;
  std::vector<RefitMeasurementSlot> slotsBufferB(static_cast<std::size_t>(activeSurfaceCount));
  const auto slotsB = assembleRefitLegSlots(seed, layerGlobals, layerMeasurements, activeSurfaceCount - 1, -1, -1, slotsBufferB);
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
      reason = OperationFailureReason::SourceSurfaceKindMismatch;
      return false;
    }
    resetCovarianceForRefit(stateC);
    float chi2C = 0.f;
    uint32_t acceptedC = 0;
    std::vector<RefitMeasurementSlot> slotsBufferC(static_cast<std::size_t>(activeSurfaceCount));
    const auto slotsC = assembleRefitLegSlots(seed, layerGlobals, layerMeasurements, 0, activeSurfaceCount, 1, slotsBufferC);
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
