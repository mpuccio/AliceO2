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

#ifndef ALICEO2_ITSMFT_TRACKING_REFITDRIVER_H_
#define ALICEO2_ITSMFT_TRACKING_REFITDRIVER_H_

#include "GPUCommonDef.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <gsl/span>

#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/TrackSeed.h"
#include "ITSMFTTracking/GlobalMeasurement.h"
#ifndef GPUCA_GPUCODE
#include "ITSMFTTracking/TimeFrame.h"
#endif
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

// Descriptor-driven refit built on Propagator operations.
namespace o2::itsmft::tracking
{

template <typename T>
struct RefitSpan {
  T* data;
  size_t count;
  GPUdi() size_t size() const { return count; }
  GPUdi() T& operator[](size_t i) const { return data[i]; }
  GPUdi() T& front() const { return data[0]; }
  GPUdi() T& back() const { return data[count - 1]; }
  GPUdi() T* begin() const { return data; }
  GPUdi() T* end() const { return data + count; }
};

namespace detail
{

constexpr float MinCircleFitBz = 0.01f; // kG

struct CircleFitPoint {
  float x, y;
  float xx, xy, yy;
};

// Preserve cancellation in a*b-c*d with two fused multiply-add operations.
GPUdi() float circleDifferenceOfProducts(float a, float b, float c, float d)
{
  const float cd = c * d;
  return std::fma(a, b, -cd) + std::fma(-c, d, cd);
}

struct CircleFloatDifference {
  float hi, lo;
};

// Return the rounded difference and its residual; do not reassociate these sums.
GPUdi() CircleFloatDifference circleTwoDiff(float a, float b)
{
  const float hi = a - b;
  const float bv = a - hi;
  return {hi, (a - (hi + bv)) + (bv - b)};
}

// Fit y = a + b*x + c*(x*x + y*y) in a frame centered on the chord.
// Compensate coordinate differences and the chord determinant to preserve
// the small sagitta in float. Cache invariant transforms for the four
// covariance-reweighting iterations; all fit arithmetic is single precision.
GPUdi() float estimateCircleQOverPt(RefitSpan<const CircleFitPoint> points, float bz) noexcept
{
  const float invalid = std::numeric_limits<float>::quiet_NaN();
  if (points.size() < 3 || points.size() > MaxLayoutSurfaces || std::abs(bz) < MinCircleFitBz) {
    return invalid;
  }
  const float x0 = points.front().x, y0 = points.front().y;
  const auto dx = circleTwoDiff(points.back().x, x0);
  const auto dy = circleTwoDiff(points.back().y, y0);
  float lengthSquared = std::fma(dx.hi, dx.hi, dy.hi * dy.hi);
  lengthSquared += 2.f * std::fma(dx.hi, dx.lo, dy.hi * dy.lo);
  const float length = std::sqrt(lengthSquared);
  if (!(length > 0.f) || !std::isfinite(length)) {
    return invalid;
  }
  const float cs = dx.hi / length, sn = dy.hi / length;
  const float invLengthSquared = 1.f / lengthSquared;
  struct CachedPoint {
    float x, y, r2, xx, xy, yy;
  };
  std::array<CachedPoint, MaxLayoutSurfaces> cache;
  for (std::size_t i = 0; i < points.size(); ++i) {
    const auto& in = points[i];
    const auto px = circleTwoDiff(in.x, x0);
    const auto py = circleTwoDiff(in.y, y0);
    // Retain subtraction residuals before dividing the small determinant.
    float cross = circleDifferenceOfProducts(dx.hi, py.hi, dy.hi, px.hi);
    float dot = std::fma(dx.hi, px.hi, dy.hi * py.hi);

    float correction = std::fma(dx.hi, py.lo, dx.lo * py.hi);
    correction = std::fma(-dy.hi, px.lo, correction);
    correction = std::fma(-dy.lo, px.hi, correction);
    correction += circleDifferenceOfProducts(dx.lo, py.lo, dy.lo, px.lo);
    cross += correction;
    dot += std::fma(dx.hi, px.lo, std::fma(dx.lo, px.hi, std::fma(dy.hi, py.lo, dy.lo * py.hi)));

    const float x = dot * invLengthSquared - .5f;
    const float y = cross * invLengthSquared;
    const float xx = in.xx, xy = in.xy, yy = in.yy;
    cache[i] = {x, y, std::fma(x, x, y * y),
                std::fma(cs * cs, xx, std::fma(2.f * cs * sn, xy, sn * sn * yy)) * invLengthSquared,
                std::fma(-cs * sn, xx, std::fma(std::fma(cs, cs, -sn * sn), xy, cs * sn * yy)) * invLengthSquared,
                std::fma(sn * sn, xx, std::fma(-2.f * cs * sn, xy, cs * cs * yy)) * invLengthSquared};
  }
  std::array<float, 3> fit{};
  for (int iteration = 0; iteration < 4; ++iteration) {
    float matrix[3][4]{};
    for (const auto& point : RefitSpan<const CachedPoint>{cache.data(), points.size()}) {
      const float nx = std::fma(-2.f * fit[2], point.x, -fit[1]);
      const float ny = std::fma(-2.f * fit[2], point.y, 1.f);
      const float variance = std::fma(nx * nx, point.xx, std::fma(2.f * nx * ny, point.xy, ny * ny * point.yy));
      if (!(variance > 0.f) || !std::isfinite(variance)) {
        return invalid;
      }

      const float weight = 1.f / variance, basis[4] = {1.f, point.x, point.r2, point.y};
      for (int i = 0; i < 3; ++i) {
        const float weighted = weight * basis[i];
        for (int j = i; j < 4; ++j) {
          matrix[i][j] = std::fma(weighted, basis[j], matrix[i][j]);
        }
      }
    }

    matrix[1][0] = matrix[0][1];
    matrix[2][0] = matrix[0][2];
    matrix[2][1] = matrix[1][2];
    // Solve the three normal equations with partial pivoting.
    for (int i = 0; i < 3; ++i) {
      int pivot = i;
      for (int j = i + 1; j < 3; ++j) {
        if (std::abs(matrix[j][i]) > std::abs(matrix[pivot][i])) {
          pivot = j;
        }
      }
      for (int k = i; k < 4; ++k) {
        std::swap(matrix[i][k], matrix[pivot][k]);
      }
      const float diagonal = matrix[i][i];
      if (std::abs(diagonal) < 1.e-15f) {
        return invalid;
      }
      for (int k = i; k < 4; ++k) {
        matrix[i][k] /= diagonal;
      }
      for (int j = 0; j < 3; ++j) {
        if (j == i) {
          continue;
        }
        const float factor = matrix[j][i];
        for (int k = i; k < 4; ++k) {
          matrix[j][k] = std::fma(-factor, matrix[i][k], matrix[j][k]);
        }
      }
    }
    for (int i = 0; i < 3; ++i) {
      fit[i] = matrix[i][3];
    }
  }
  const float discriminant = std::fma(-4.f * fit[0], fit[2], std::fma(fit[1], fit[1], 1.f));
  return discriminant > 0.f ? 2.f * fit[2] / (length * std::sqrt(discriminant) * bz * o2::constants::math::B2C) : invalid;
}

struct RefitMeasurementSlot {
  SurfaceMeasurement measurement{};
  LayerId surface{};
  bool present{false};
};

// Holes are skipped; present slots must resolve to a descriptor. Commit state,
// reference, chi2 and count only after the full leg succeeds.
GPUdi() bool driveRefitLeg(SurfaceTrackState& state, SurfaceTrackParameters& linRef,
                           float& chi2, uint32_t& acceptedHitCount,
                           RefitSpan<const RefitMeasurementSlot> orderedSlots, SurfaceCatalogView surfaceCatalog,
                           float bz, material::MaterialTraversalDirection direction,
                           bool shiftReferenceToMeasurement, float maxChi2) noexcept
{
  if (chi2 < 0.f) {
    return false;
  }

  SurfaceTrackState scratchState = state;
  SurfaceTrackParameters scratchLinRef = linRef;
  float scratchChi2 = chi2;
  uint32_t scratchAcceptedHitCount = 0;
  constexpr uint32_t kChi2GateMinAcceptedHits = 3;
  for (const auto& slot : orderedSlots) {
    if (!slot.present) {
      continue;
    }
    if (!slot.surface.isValid() || !(surfaceCatalog.nSurfaces == 0 || surfaceCatalog.surfaces != nullptr) ||
        !(slot.surface.value() < surfaceCatalog.nSurfaces)) {
      return false;
    }
    const auto& covariance = slot.measurement.covariance;
    if (!std::isfinite(covariance.uu) || !std::isfinite(covariance.uv) || !std::isfinite(covariance.vv) ||
        covariance.uu < 0.f || covariance.vv < 0.f) {
      return false;
    }
    const SurfaceDescriptor& descriptor = surfaceCatalog.surfaces[slot.surface.value()];
    if (!Propagator::propagateToMeasurement(scratchState, scratchLinRef, descriptor, slot.measurement, bz, direction,
                                            scratchAcceptedHitCount >= kChi2GateMinAcceptedHits, maxChi2, scratchChi2,
                                            shiftReferenceToMeasurement)) {
      return false;
    }
    ++scratchAcceptedHitCount;
  }
  state = scratchState;
  linRef = scratchLinRef;
  chi2 = scratchChi2;
  acceptedHitCount = scratchAcceptedHitCount;
  return true;
}

} // namespace detail

// Common first-pass prior for the two position coordinates, direction and q/pT.
GPUdi() void resetCovarianceForRefit(SurfaceTrackState& state) noexcept
{
  for (auto& element : state.covariance) {
    element = 0.f;
  }
  for (int i = 0; i < 4; ++i) {
    state.covariance[packedCovarianceIndex(i, i)] = 1.f;
  }
  // This is the variance, not the standard deviation.
  state.covariance[packedCovarianceIndex(4, 4)] = std::clamp(std::abs(state.parameters[4]), 1.f, 10.f);
}

// Start a subsequent leg with five times the previous parameter uncertainties.
GPUdi() void inflateDiagonalCovarianceForRefit(SurfaceTrackState& state) noexcept
{
  constexpr float varianceInflation = 25.f;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < i; ++j) {
      state.covariance[packedCovarianceIndex(i, j)] = 0.f;
    }
    state.covariance[packedCovarianceIndex(i, i)] *= varianceInflation;
  }
}

// parameters[4] is signed q/pT for both coordinate conventions.
GPUdi() float ptFromQOverPt(float q2pt, uint8_t absCharge) noexcept
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

// Refit inward, outward, then optionally inward again; commit on success.
GPUdi() bool fitPreparedTrackStateLegs(
  const SurfaceTrackState& seedState,
  RefitSpan<const detail::RefitMeasurementSlot> slots,
  RefitSpan<const detail::CircleFitPoint> points,
  SurfaceCatalogView surfaceCatalog,
  float bz,
  bool shiftReferenceToMeasurement,
  float maxChi2ClusterAttachment,
  float maxChi2NDF,
  bool repeatRefitOut,
  float minPtThreshold,
  SurfaceTrackState& outParamIn,
  SurfaceTrackState& outParamOut,
  float& outChi2) noexcept
{
  if (slots.size() == 0 || slots.size() > MaxLayoutSurfaces) {
    return false;
  }
  std::array<detail::RefitMeasurementSlot, MaxLayoutSurfaces> reverseSlots{};
  for (size_t i = 0; i < slots.size(); ++i) {
    reverseSlots[i] = slots[slots.size() - 1 - i];
  }
  auto legAcceptable = [](const SurfaceTrackState& state, float chi2, uint32_t acceptedHitCount,
                          float maxQoverPt, float maxChi2NDFValue) noexcept -> bool {
    if (!(std::abs(state.parameters[4]) < maxQoverPt)) {
      return false;
    }
    return chi2 < maxChi2NDFValue * static_cast<float>(static_cast<int>(acceptedHitCount) * 2 - 5);
  };

  // Leg A: inward.
  SurfaceTrackState stateA = seedState;
  if (!std::isfinite(bz)) {
    return false;
  }
  // There is no curvature constraint with the field off; keep the CA seed.
  if (std::abs(bz) >= detail::MinCircleFitBz) {
    const float qOverPt = detail::estimateCircleQOverPt(points, bz);
    if (!std::isfinite(qOverPt)) {
      return false;
    }
    stateA.parameters[4] = qOverPt;
  }
  SurfaceTrackParameters linRefA{stateA};
  resetCovarianceForRefit(stateA);
  float chi2A = 0.f;
  uint32_t acceptedA = 0;
  const auto slotsA = slots;
  if (!detail::driveRefitLeg(stateA, linRefA, chi2A, acceptedA, slotsA, surfaceCatalog, bz,
                             material::MaterialTraversalDirection::AlongMomentum, shiftReferenceToMeasurement,
                             maxChi2ClusterAttachment)) {
    return false;
  }
  if (!legAcceptable(stateA, chi2A, acceptedA, o2::constants::math::VeryBig, maxChi2NDF)) {
    return false;
  }

  // Leg B: outward; this is the reported inner result.
  SurfaceTrackState stateB = stateA;
  SurfaceTrackParameters linRefB{stateB};
  inflateDiagonalCovarianceForRefit(stateB);
  float chi2B = 0.f;
  uint32_t acceptedB = 0;
  const RefitSpan<const detail::RefitMeasurementSlot> slotsB{reverseSlots.data(), slots.size()};
  if (!detail::driveRefitLeg(stateB, linRefB, chi2B, acceptedB, slotsB, surfaceCatalog, bz,
                             material::MaterialTraversalDirection::OppositeMomentum, shiftReferenceToMeasurement,
                             maxChi2ClusterAttachment)) {
    return false;
  }
  if (!legAcceptable(stateB, chi2B, acceptedB, 50.f, maxChi2NDF)) {
    return false;
  }

  if (minPtThreshold > 0.f && ptFromQOverPt(stateB.parameters[4], stateB.absCharge) < minPtThreshold) {
    return false;
  }

  // Optional leg C: inward again.
  SurfaceTrackState stateOut = stateA;
  if (repeatRefitOut) {
    SurfaceTrackState stateC = stateB;
    SurfaceTrackParameters linRefC{stateC};
    inflateDiagonalCovarianceForRefit(stateC);
    float chi2C = 0.f;
    uint32_t acceptedC = 0;
    const auto slotsC = slots;
    if (!detail::driveRefitLeg(stateC, linRefC, chi2C, acceptedC, slotsC, surfaceCatalog, bz,
                               material::MaterialTraversalDirection::AlongMomentum, shiftReferenceToMeasurement,
                               maxChi2ClusterAttachment)) {
      return false;
    }
    if (!legAcceptable(stateC, chi2C, acceptedC, o2::constants::math::VeryBig, maxChi2NDF)) {
      return false;
    }
    stateOut = stateC;
  }

  outParamIn = stateB;
  outParamOut = stateOut;
  outChi2 = chi2B;
  return true;
}

// A seed's refit measurements: one slot per layer (holes stay empty) and the
// circle-fit points of its hits. Returns the number of slots, 0 if one of its
// measurements is missing.
GPUhdi() size_t prepareRefitSlots(const FrameView& frame, const TrackSeed& seed, detail::RefitMeasurementSlot* slots,
                                  detail::CircleFitPoint* points, size_t& nPoints)
{
  nPoints = 0;
  for (size_t layer = 0; layer < frame.nLayers; ++layer) {
    const int index = seed.getCluster(static_cast<int>(layer));
    slots[layer] = {};
    if (index == o2::its::constants::UnusedIndex) {
      continue;
    }
    const auto* cluster = frame.cluster(static_cast<int>(layer), index);
    const auto* measurement = cluster ? frame.layers[layer].measurement(cluster->clusterId) : nullptr;
    if (!measurement) {
      nPoints = 0;
      return 0;
    }
    slots[layer] = {*measurement, LayerId{static_cast<uint16_t>(layer)}, true};
    points[nPoints++] = {cluster->x, cluster->y, cluster->covariance.xx, cluster->covariance.xy, cluster->covariance.yy};
  }
  return frame.nLayers;
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_REFITDRIVER_H_ */
