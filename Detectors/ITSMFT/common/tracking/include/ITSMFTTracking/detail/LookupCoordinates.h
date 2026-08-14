// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_LOOKUPCOORDINATES_H_
#define ALICEO2_ITSMFT_TRACKING_LOOKUPCOORDINATES_H_

#include <algorithm>
#include <cmath>

#include "CommonConstants/MathConstants.h"
#include "GPUCommonMath.h"
#include "GPUCommonDef.h"
#include "ITSMFTTracking/GlobalMeasurement.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"

namespace o2::itsmft::tracking
{

struct LookupCoordinates {
  float phi{0.f};
  float transverse{0.f};
  float covariance[3]{}; // (transverse, transverse-phi, phi)
};

struct LookupWindow {
  float phiMin{0.f};
  float phiMax{0.f};
  float transverseMin{0.f};
  float transverseMax{0.f};
  bool wrapsPhi{false};
};

GPUhdi() float normalizeLookupPhi(float phi) noexcept
{
  phi -= o2::constants::math::TwoPI * o2::gpu::GPUCommonMath::Floor(phi * (1.f / o2::constants::math::TwoPI));
  return phi;
}

GPUhdi() bool makeLookupCoordinates(const SurfaceDescriptor& surface,
                                    const GlobalMeasurement& measurement,
                                    LookupCoordinates& out) noexcept
{
  const float x = measurement.position.x;
  const float y = measurement.position.y;
  const float r2 = x * x + y * y;
  if (!(surface.kind == SurfaceKind::Cylinder || surface.kind == SurfaceKind::Disk) ||
      !o2::gpu::GPUCommonMath::Finite(x) || !o2::gpu::GPUCommonMath::Finite(y) || !(r2 > 0.f)) {
    return false;
  }
  const float r = o2::gpu::GPUCommonMath::Sqrt(r2);
  const float inverseR = 1.f / r;
  const float inverseR2 = inverseR * inverseR;
  const float phiX = -y * inverseR2;
  const float phiY = x * inverseR2;
  const auto& c = measurement.covariance;
  const float varPhi = phiX * phiX * c.xx + 2.f * phiX * phiY * c.xy + phiY * phiY * c.yy;
  if (!o2::gpu::GPUCommonMath::Finite(varPhi)) {
    return false;
  }
  if (surface.kind == SurfaceKind::Cylinder) {
    const float covZPhi = c.xz * phiX + c.yz * phiY;
    if (!o2::gpu::GPUCommonMath::Finite(measurement.position.z) ||
        !o2::gpu::GPUCommonMath::Finite(covZPhi) || !o2::gpu::GPUCommonMath::Finite(c.zz)) {
      return false;
    }
    out = {normalizeLookupPhi(o2::gpu::GPUCommonMath::ATan2(y, x)), measurement.position.z, {c.zz, covZPhi, varPhi}};
    return true;
  }
  const float radialX = x * inverseR;
  const float radialY = y * inverseR;
  const float varR = radialX * radialX * c.xx + 2.f * radialX * radialY * c.xy + radialY * radialY * c.yy;
  const float covRPhi = radialX * phiX * c.xx + (radialX * phiY + radialY * phiX) * c.xy + radialY * phiY * c.yy;
  if (!o2::gpu::GPUCommonMath::Finite(varR) || !o2::gpu::GPUCommonMath::Finite(covRPhi)) {
    return false;
  }
  out = {normalizeLookupPhi(o2::gpu::GPUCommonMath::ATan2(y, x)), r, {varR, covRPhi, varPhi}};
  return true;
}

inline bool makeLookupWindow(const LookupCoordinates& coordinates, SurfaceChartRange range,
                             float nSigmaCut, LookupWindow& out) noexcept
{
  const float determinant = coordinates.covariance[0] * coordinates.covariance[2] -
                            coordinates.covariance[1] * coordinates.covariance[1];
  if (!range.isValid() || !std::isfinite(nSigmaCut) || !(nSigmaCut >= 0.f) ||
      !std::isfinite(coordinates.phi) || !std::isfinite(coordinates.transverse) ||
      !(coordinates.covariance[0] >= 0.f) || !(coordinates.covariance[2] >= 0.f) ||
      !(determinant >= 0.f)) {
    return false;
  }
  const float phiHalfWidth = nSigmaCut * std::sqrt(coordinates.covariance[2]);
  const float transverseHalfWidth = nSigmaCut * std::sqrt(coordinates.covariance[0]);
  if (!std::isfinite(phiHalfWidth) || !std::isfinite(transverseHalfWidth)) {
    return false;
  }
  out.transverseMin = std::max(range.min, coordinates.transverse - transverseHalfWidth);
  out.transverseMax = std::min(range.max, coordinates.transverse + transverseHalfWidth);
  if (!(out.transverseMin <= out.transverseMax)) {
    return false;
  }
  if (phiHalfWidth >= o2::constants::math::PI) {
    out.phiMin = 0.f;
    out.phiMax = o2::constants::math::TwoPI;
    out.wrapsPhi = false;
  } else {
    out.phiMin = normalizeLookupPhi(coordinates.phi - phiHalfWidth);
    out.phiMax = normalizeLookupPhi(coordinates.phi + phiHalfWidth);
    out.wrapsPhi = out.phiMin > out.phiMax;
  }
  return true;
}

} // namespace o2::itsmft::tracking

#endif
