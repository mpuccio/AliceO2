// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_DETAIL_TRACKINGKERNELPARAMETERS_H_
#define ALICEO2_ITSMFT_TRACKING_DETAIL_TRACKINGKERNELPARAMETERS_H_

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "GPUCommonDef.h"
#include "GPUCommonMath.h"

namespace o2::itsmft::tracking
{

/// True for finite single-precision values, on host and device alike.
GPUhdi() bool isFiniteParam(float x) noexcept
{
  return (o2::gpu::GPUCommonMath::Float2UIntReint(x) & 0x7f800000u) != 0x7f800000u;
}

/// Compact device-facing tracking configuration. Lengths are in cm, momentum in GeV/c,
/// angles and their resolutions in radians, and chi-square quantities are
/// dimensionless.
struct TrackingKernelParameters {
  float trackletMinPt{0.3f};
  float nSigmaCut{5.f};
  float maxChi2ClusterAttachment{60.f};
  float maxChi2NDF{30.f};
  float pvResolution{1.e-2f};

  GPUhdi() bool isValid() const noexcept
  {
    if (!isFiniteParam(trackletMinPt) || trackletMinPt <= 0.f ||
        !isFiniteParam(nSigmaCut) || nSigmaCut <= 0.f ||
        !isFiniteParam(maxChi2ClusterAttachment) || maxChi2ClusterAttachment <= 0.f ||
        !isFiniteParam(maxChi2NDF) || maxChi2NDF <= 0.f) {
      return false;
    }
    return isFiniteParam(pvResolution) && pvResolution >= 0.f;
  }
};

static_assert(std::is_standard_layout_v<TrackingKernelParameters>);
static_assert(std::is_trivially_copyable_v<TrackingKernelParameters>);
static_assert(sizeof(TrackingKernelParameters) == 20);
static_assert(alignof(TrackingKernelParameters) == alignof(float));
static_assert(offsetof(TrackingKernelParameters, trackletMinPt) == 0);
static_assert(offsetof(TrackingKernelParameters, nSigmaCut) == 4);
static_assert(offsetof(TrackingKernelParameters, maxChi2ClusterAttachment) == 8);
static_assert(offsetof(TrackingKernelParameters, maxChi2NDF) == 12);
static_assert(offsetof(TrackingKernelParameters, pvResolution) == 16);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_DETAIL_TRACKINGKERNELPARAMETERS_H_ */
