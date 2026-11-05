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

/// \file Trig.h
/// \brief The transcendental functions of shared tracking code, overloaded for
///        float and double, for host and device alike.
///
/// Outside O2's deterministic mode they are `CAMath`'s: `Sin`, `ATan2`, ... in
/// float and `Sind`, `ATan2d`, ... in double. `sincosOf` is
/// `o2::math_utils::detail::sincos`, one call on hosts and devices, and the
/// spelling O2's own track code uses.
///
/// O2's deterministic mode is different. It computes the float functions in
/// double and rounds, which almost always hides the host's and the device's
/// libraries disagreeing in the last place of the double; the double
/// functions it leaves native, and there nothing hides it (PortableMath.h). So in that
/// mode every function here is PortableMath's,
/// in both precisions: the same arithmetic on the host and on every device,
/// the float ones in double and rounded, as O2 does.

#ifndef GTRK_TRACKING_TRIG_H_
#define GTRK_TRACKING_TRIG_H_

#include "GenericTracking/Base/Definitions.h"

#include "GenericTracking/Base/PortableMath.h"

#include "GPUCommonMath.h"
#include "MathUtils/detail/trigonometric.h"

// The documented libdevice float implementation, called directly so CUDA's
// --use-fast-math does not substitute the lower-accuracy intrinsic.
#if defined(__CUDACC__) && !defined(__HIPCC__)
extern "C" __device__ void __nv_sincosf(float, float*, float*);
#endif

namespace o2::gtrk
{

GPUhdi() GPUAlwaysInline() float sinOf(float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::sin(static_cast<double>(x)));
#else
  return o2::gpu::CAMath::Sin(x);
#endif
}

GPUhdi() GPUAlwaysInline() float cosOf(float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::cos(static_cast<double>(x)));
#else
  return o2::gpu::CAMath::Cos(x);
#endif
}

GPUhdi() GPUAlwaysInline() float asinOf(float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::asin(static_cast<double>(x)));
#else
  return o2::gpu::CAMath::ASin(x);
#endif
}

GPUhdi() GPUAlwaysInline() float atanOf(float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::atan(static_cast<double>(x)));
#else
  return o2::gpu::CAMath::ATan(x);
#endif
}

GPUhdi() GPUAlwaysInline() float logOf(float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::log(static_cast<double>(x)));
#else
  return o2::gpu::CAMath::Log(x);
#endif
}

GPUhdi() GPUAlwaysInline() float atan2Of(float y, float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::atan2(static_cast<double>(y), static_cast<double>(x)));
#else
  return o2::gpu::CAMath::ATan2(y, x);
#endif
}

GPUhdi() GPUAlwaysInline() void sincosOf(float x, float& s, float& c)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  s = sinOf(x);
  c = cosOf(x);
#else
  math_utils::detail::sincos<float>(x, s, c);
#endif
}

#ifdef GPUCA_DETERMINISTIC_MODE
GPUhdi() GPUAlwaysInline() double sinOf(double x) { return portable::sin(x); }
GPUhdi() GPUAlwaysInline() double cosOf(double x) { return portable::cos(x); }
GPUhdi() GPUAlwaysInline() double asinOf(double x) { return portable::asin(x); }
GPUhdi() GPUAlwaysInline() double atanOf(double x) { return portable::atan(x); }
GPUhdi() GPUAlwaysInline() double atan2Of(double y, double x) { return portable::atan2(y, x); }
#else
GPUhdi() GPUAlwaysInline() double sinOf(double x) { return o2::gpu::CAMath::Sind(x); }
GPUhdi() GPUAlwaysInline() double cosOf(double x) { return o2::gpu::CAMath::Cosd(x); }
GPUhdi() GPUAlwaysInline() double asinOf(double x) { return o2::gpu::CAMath::ASind(x); }
GPUhdi() GPUAlwaysInline() double atanOf(double x) { return o2::gpu::CAMath::ATand(x); }
GPUhdi() GPUAlwaysInline() double atan2Of(double y, double x) { return o2::gpu::CAMath::ATan2d(y, x); }
#endif

GPUhdi() GPUAlwaysInline() void sincosOf(double x, double& s, double& c)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  s = portable::sin(x);
  c = portable::cos(x);
#else
  math_utils::detail::sincos<double>(x, s, c);
#endif
}

/// sin and cos of a direction for the double plane solve. Under CUDA's
/// -use_fast_math (outside the deterministic mode) `sincosf` becomes the
/// intrinsic, too coarse here; libdevice's precise one is called instead.
GPUhdi() GPUAlwaysInline() void directionSincosOf(float angle, double& sine, double& cosine)
{
  float s = 0.F;
  float c = 0.F;
#if defined(__CUDA_ARCH__) && !defined(GPUCA_DETERMINISTIC_MODE)
  __nv_sincosf(angle, &s, &c);
#else
  sincosOf(angle, s, c);
#endif
  sine = s;
  cosine = c;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRIG_H_
