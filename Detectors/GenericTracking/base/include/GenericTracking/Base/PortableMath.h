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

/// \file PortableMath.h
/// \brief sin, cos, atan, atan2, asin and log in double, the same bits on the
///        host and on every device.
///
/// The deterministic build asks the host and the devices for bit-identical
/// results, and their math libraries do not give them: none of glibc's,
/// CUDA's or ROCm's double sin, cos or atan2 is correctly rounded, and each
/// rounds differently in the last place. CUDA's atan2 differs from glibc's
/// on many of the propagator's turn angles, by up to 2 ulp. O2's deterministic
/// mode survives this in float, since a double result rounded to float hides
/// an ulp of double almost always; the propagator's arithmetic is double, and
/// the ulp reaches the tracks' times.
///
/// So the deterministic build computes them here, from fdlibm's algorithms as
/// FreeBSD's msun has them: only +, -, *, / and sqrt, which IEEE 754 rounds
/// correctly everywhere, and integer operations on the bits. With contraction
/// off -- -ffp-contract=off, nvcc's --fmad=false, as the deterministic build
/// sets -- every step is the same operation on the same operands on any
/// target. The accuracy is fdlibm's: under 1 ulp, atan2 under 2.
///
/// Arguments of sin and cos beyond 2^20 pi/2 (1.6e6) fall back to the native
/// functions: fdlibm reduces them with a 1584-bit table of 2/pi, and no
/// tracking angle comes near. A NaN argument is returned as it is, and an
/// invalid one gives one quiet NaN, so even those are the same bits.

#ifndef GTRK_BASE_PORTABLEMATH_H_
#define GTRK_BASE_PORTABLEMATH_H_

#include "GenericTracking/Base/Definitions.h"

#include "GPUCommonMath.h"

#ifndef GPUCA_GPUCODE_DEVICE
#include <bit>
#endif
#include <cstdint>
#include <numbers>

namespace o2::gtrk::portable
{

namespace detail
{

GPUhdi() GPUAlwaysInline() std::uint64_t bits(double x)
{
#ifdef GPUCA_GPUCODE_DEVICE
  return static_cast<std::uint64_t>(__double_as_longlong(x));
#else
  return std::bit_cast<std::uint64_t>(x);
#endif
}

GPUhdi() GPUAlwaysInline() double fromBits(std::uint64_t b)
{
#ifdef GPUCA_GPUCODE_DEVICE
  return __longlong_as_double(static_cast<long long>(b));
#else
  return std::bit_cast<double>(b);
#endif
}

/// The high word, fdlibm's GET_HIGH_WORD.
GPUhdi() GPUAlwaysInline() std::uint32_t high(double x) { return static_cast<std::uint32_t>(bits(x) >> 32U); }
GPUhdi() GPUAlwaysInline() std::uint32_t low(double x) { return static_cast<std::uint32_t>(bits(x)); }
GPUhdi() GPUAlwaysInline() double withHigh(double x, std::uint32_t word) { return fromBits((std::uint64_t{word} << 32U) | low(x)); }
GPUhdi() GPUAlwaysInline() double withLowZero(double x) { return fromBits(bits(x) & 0xFFFFFFFF00000000ULL); }
GPUhdi() GPUAlwaysInline() double absOf(double x) { return fromBits(bits(x) & 0x7FFFFFFFFFFFFFFFULL); }
GPUhdi() GPUAlwaysInline() bool negative(double x) { return (bits(x) >> 63U) != 0U; }
/// By the bits, not by comparison: -Wfloat-equal, and no fast-math
/// assumption about NaNs can remove them.
GPUhdi() GPUAlwaysInline() bool isNaN(double x) { return (bits(x) & 0x7FFFFFFFFFFFFFFFULL) > 0x7FF0000000000000ULL; }
GPUhdi() GPUAlwaysInline() bool isZero(double x) { return (bits(x) & 0x7FFFFFFFFFFFFFFFULL) == 0U; }
GPUhdi() GPUAlwaysInline() int biasedExponent(double x) { return static_cast<int>((high(x) >> 20U) & 0x7FFU); }

GPUhdi() GPUAlwaysInline() double quietNaN() { return fromBits(0x7FF8000000000000ULL); }
GPUhdi() GPUAlwaysInline() double infinity() { return fromBits(0x7FF0000000000000ULL); }

/// pi/2, pi/4 and pi/2's tail (pi/2 - kPio2), as fdlibm splits them.
constexpr double kPio2 = 1.57079632679489655800e+00;
constexpr double kPio4 = 7.85398163397448278999e-01;
constexpr double kPio2Lo = 6.12323399573676603587e-17;

/// sin on [-pi/4, pi/4] of x + y, y the tail of a reduced argument;
/// \a hasTail false when y is zero. fdlibm's __kernel_sin.
GPUhdi() GPUAlwaysInline() double kernelSin(double x, double y, bool hasTail)
{
  constexpr double kS1 = -1.66666666666666324348e-01;
  constexpr double kS2 = 8.33333333332248946124e-03;
  constexpr double kS3 = -1.98412698298579493134e-04;
  constexpr double kS4 = 2.75573137070700676789e-06;
  constexpr double kS5 = -2.50507602534068634195e-08;
  constexpr double kS6 = 1.58969099521155010221e-10;
  const double z = x * x;
  const double w = z * z;
  const double r = kS2 + z * (kS3 + z * kS4) + z * w * (kS5 + z * kS6);
  const double v = z * x;
  if (!hasTail) {
    return x + v * (kS1 + z * r);
  }
  return x - ((z * (0.5 * y - v * r) - y) - v * kS1);
}

/// cos on [-pi/4, pi/4] of x + y. fdlibm's __kernel_cos, as msun has it.
GPUhdi() GPUAlwaysInline() double kernelCos(double x, double y)
{
  constexpr double kC1 = 4.16666666666666019037e-02;
  constexpr double kC2 = -1.38888888888741095749e-03;
  constexpr double kC3 = 2.48015872894767294178e-05;
  constexpr double kC4 = -2.75573143513906633035e-07;
  constexpr double kC5 = 2.08757232129817482790e-09;
  constexpr double kC6 = -1.13596475577881948265e-11;
  const double z = x * x;
  double w = z * z;
  const double r = z * (kC1 + z * (kC2 + z * kC3)) + w * w * (kC4 + z * (kC5 + z * kC6));
  const double hz = 0.5 * z;
  w = 1.0 - hz;
  return w + (((1.0 - w) - hz) + (z * r - x * y));
}

/// x - n pi/2 as y0 + y1, and n, for pi/4 < |x| < 2^20 pi/2: fdlibm's
/// __ieee754_rem_pio2, the medium-size path, pi/2 in three parts of 33 bits
/// so that n times each is exact.
GPUhdi() GPUAlwaysInline() int remPio2(double x, double& y0, double& y1)
{
  constexpr double kInvPio2 = 6.36619772367581382433e-01;
  constexpr double kPio2Part1 = 1.57079632673412561417e+00;
  constexpr double kPio2Tail1 = 6.07710050650619224932e-11;
  constexpr double kPio2Part2 = 6.07710050630396597660e-11;
  constexpr double kPio2Tail2 = 2.02226624879595063154e-21;
  constexpr double kPio2Part3 = 2.02226624871116645580e-21;
  constexpr double kPio2Tail3 = 8.47842766036889956997e-32;
  // 1.5 * 2^52: adding and subtracting it rounds to the nearest integer.
  constexpr double kToInteger = 6755399441055744.0;
  const double fn = (x * kInvPio2 + kToInteger) - kToInteger;
  const int n = static_cast<int>(fn);
  double r = x - fn * kPio2Part1;
  double w = fn * kPio2Tail1;
  const int j = biasedExponent(x);
  y0 = r - w;
  if (j - biasedExponent(y0) > 16) { // cancellation: a second round, 118 bits
    double t = r;
    w = fn * kPio2Part2;
    r = t - w;
    w = fn * kPio2Tail2 - ((t - r) - w);
    y0 = r - w;
    if (j - biasedExponent(y0) > 49) { // and a third, 151 bits
      t = r;
      w = fn * kPio2Part3;
      r = t - w;
      w = fn * kPio2Tail3 - ((t - r) - w);
      y0 = r - w;
    }
  }
  y1 = (r - y0) - w;
  return n;
}

/// atan(c) for c = 1/2, 1, 3/2 and infinity, as a high and a low part. A
/// switch, not a local array: HIP, compiling with real calls and relocatable
/// device code, reads wrong values from such an array.
GPUhdi() GPUAlwaysInline() double atanHi(int id)
{
  switch (id) {
    case 0:
      return 4.63647609000806093515e-01;
    case 1:
      return kPio4;
    case 2:
      return 9.82793723247329054082e-01;
    default:
      return kPio2;
  }
}

GPUhdi() GPUAlwaysInline() double atanLo(int id)
{
  switch (id) {
    case 0:
      return 2.26987774529616870924e-17;
    case 1:
      return 3.06161699786838301793e-17;
    case 2:
      return 1.39033110312309984516e-17;
    default:
      return kPio2Lo;
  }
}

/// |x| at or above 2^20 pi/2, where remPio2 no longer applies.
constexpr std::uint32_t kLargeArgument = 0x413921FBU;

} // namespace detail

GPUhdi() GPUAlwaysInline() double sin(double x)
{
  const std::uint32_t ix = detail::high(x) & 0x7FFFFFFFU;
  if (ix <= 0x3FE921FBU) {  // |x| <~ pi/4
    if (ix < 0x3E500000U) { // |x| < 2^-26: sin x is x
      return x;
    }
    return detail::kernelSin(x, 0., false);
  }
  if (ix >= 0x7FF00000U) {
    return detail::isNaN(x) ? x : detail::quietNaN();
  }
  if (ix >= detail::kLargeArgument) {
    return o2::gpu::CAMath::Sind(x); // beyond remPio2's reach
  }
  double y0 = 0.;
  double y1 = 0.;
  switch (detail::remPio2(x, y0, y1) & 3) {
    case 0:
      return detail::kernelSin(y0, y1, true);
    case 1:
      return detail::kernelCos(y0, y1);
    case 2:
      return -detail::kernelSin(y0, y1, true);
    default:
      return -detail::kernelCos(y0, y1);
  }
}

GPUhdi() GPUAlwaysInline() double cos(double x)
{
  const std::uint32_t ix = detail::high(x) & 0x7FFFFFFFU;
  if (ix <= 0x3FE921FBU) {
    if (ix < 0x3E46A09EU) { // |x| < 2^-27 sqrt(2): cos x is 1
      return 1.;
    }
    return detail::kernelCos(x, 0.);
  }
  if (ix >= 0x7FF00000U) {
    return detail::isNaN(x) ? x : detail::quietNaN();
  }
  if (ix >= detail::kLargeArgument) {
    return o2::gpu::CAMath::Cosd(x); // beyond remPio2's reach
  }
  double y0 = 0.;
  double y1 = 0.;
  switch (detail::remPio2(x, y0, y1) & 3) {
    case 0:
      return detail::kernelCos(y0, y1);
    case 1:
      return -detail::kernelSin(y0, y1, true);
    case 2:
      return -detail::kernelCos(y0, y1);
    default:
      return detail::kernelSin(y0, y1, true);
  }
}

/// fdlibm's s_atan.c: atan(x) = atan(c) + atan((x - c) / (1 + x c)) for c
/// one of 0, 1/2, 1, 3/2 and infinity, then an odd polynomial.
GPUhdi() GPUAlwaysInline() double atan(double x)
{
  constexpr double kAt0 = 3.33333333333329318027e-01;
  constexpr double kAt1 = -1.99999999998764832476e-01;
  constexpr double kAt2 = 1.42857142725034663711e-01;
  constexpr double kAt3 = -1.11111104054623557880e-01;
  constexpr double kAt4 = 9.09088713343650656196e-02;
  constexpr double kAt5 = -7.69187620504482999495e-02;
  constexpr double kAt6 = 6.66107313738753120669e-02;
  constexpr double kAt7 = -5.83357013379057348645e-02;
  constexpr double kAt8 = 4.97687799461593236017e-02;
  constexpr double kAt9 = -3.65315727442169155270e-02;
  constexpr double kAt10 = 1.62858201153657823623e-02;

  const bool sign = detail::negative(x);
  const std::uint32_t ix = detail::high(x) & 0x7FFFFFFFU;
  if (ix >= 0x44100000U) { // |x| >= 2^66, or not finite
    if (detail::isNaN(x)) {
      return x;
    }
    return sign ? -detail::atanHi(3) - detail::atanLo(3) : detail::atanHi(3) + detail::atanLo(3);
  }
  int id = -1;
  if (ix < 0x3FDC0000U) {   // |x| < 0.4375
    if (ix < 0x3E400000U) { // |x| < 2^-27: atan x is x
      return x;
    }
  } else {
    x = detail::absOf(x);
    if (ix < 0x3FF30000U) {   // |x| < 1.1875
      if (ix < 0x3FE60000U) { // 7/16 <= |x| < 11/16
        id = 0;
        x = (2.0 * x - 1.0) / (2.0 + x);
      } else { // 11/16 <= |x| < 19/16
        id = 1;
        x = (x - 1.0) / (x + 1.0);
      }
    } else if (ix < 0x40038000U) { // |x| < 2.4375
      id = 2;
      x = (x - 1.5) / (1.0 + 1.5 * x);
    } else { // 2.4375 <= |x| < 2^66
      id = 3;
      x = -1.0 / x;
    }
  }
  const double z = x * x;
  const double w = z * z;
  // The odd and the even terms of the polynomial, apart.
  const double s1 = z * (kAt0 + w * (kAt2 + w * (kAt4 + w * (kAt6 + w * (kAt8 + w * kAt10)))));
  const double s2 = w * (kAt1 + w * (kAt3 + w * (kAt5 + w * (kAt7 + w * kAt9))));
  if (id < 0) {
    return x - x * (s1 + s2);
  }
  const double result = detail::atanHi(id) - ((x * (s1 + s2) - detail::atanLo(id)) - x);
  return sign ? -result : result;
}

/// fdlibm's e_atan2.c: atan(|y / x|), placed in its quadrant.
GPUhdi() GPUAlwaysInline() double atan2(double y, double x)
{
  constexpr double kPi = std::numbers::pi;
  constexpr double kPiLo = 2.0 * detail::kPio2Lo;
  constexpr double kPiOver2 = detail::kPio2;
  constexpr double kPiOver4 = detail::kPio4;

  if (detail::isNaN(x)) {
    return x;
  }
  if (detail::isNaN(y)) {
    return y;
  }
  if (detail::bits(x) == 0x3FF0000000000000ULL) { // x = 1
    return atan(y);
  }
  const bool yNegative = detail::negative(y);
  const bool xNegative = detail::negative(x);
  const auto ix = static_cast<std::int32_t>(detail::high(x) & 0x7FFFFFFFU);
  const auto iy = static_cast<std::int32_t>(detail::high(y) & 0x7FFFFFFFU);

  if (detail::isZero(y)) {
    if (!xNegative) {
      return y; // atan2(+-0, +anything) = +-0
    }
    return yNegative ? -kPi : kPi; // atan2(+-0, -anything) = +-pi
  }
  if (detail::isZero(x)) {
    return yNegative ? -kPiOver2 : kPiOver2;
  }
  if (ix == 0x7FF00000) { // x infinite
    if (iy == 0x7FF00000) {
      const double angle = xNegative ? 3.0 * kPiOver4 : kPiOver4;
      return yNegative ? -angle : angle;
    }
    const double angle = xNegative ? kPi : 0.;
    return yNegative ? -angle : angle;
  }
  if (iy == 0x7FF00000) {
    return yNegative ? -kPiOver2 : kPiOver2;
  }

  double z = 0.;
  const std::int32_t k = (iy - ix) >> 20;
  bool leftHalf = xNegative;
  if (k > 60) { // |y / x| > 2^60
    z = kPiOver2 + 0.5 * kPiLo;
    leftHalf = false;
  } else if (xNegative && k < -60) { // 0 > |y| / x > -2^-60
    z = 0.;
  } else {
    z = atan(detail::absOf(y / x));
  }
  if (!leftHalf) {
    return yNegative ? -z : z;
  }
  return yNegative ? (z - kPiLo) - kPi : kPi - (z - kPiLo);
}

/// fdlibm's e_asin.c: a rational approximation below 1/2, and above it
/// asin x = pi/2 - 2 asin(sqrt((1 - x) / 2)).
GPUhdi() GPUAlwaysInline() double asin(double x)
{
  constexpr double kPio2Hi = detail::kPio2;
  constexpr double kPio2Lo = detail::kPio2Lo;
  constexpr double kPio4Hi = detail::kPio4;
  constexpr double kPS0 = 1.66666666666666657415e-01;
  constexpr double kPS1 = -3.25565818622400915405e-01;
  constexpr double kPS2 = 2.01212532134862925881e-01;
  constexpr double kPS3 = -4.00555345006794114027e-02;
  constexpr double kPS4 = 7.91534994289814532176e-04;
  constexpr double kPS5 = 3.47933107596021167570e-05;
  constexpr double kQS1 = -2.40339491173441421878e+00;
  constexpr double kQS2 = 2.02094576023350569471e+00;
  constexpr double kQS3 = -6.88283971605453293030e-01;
  constexpr double kQS4 = 7.70381505559019352791e-02;

  const std::uint32_t ix = detail::high(x) & 0x7FFFFFFFU;
  if (ix >= 0x3FF00000U) { // |x| >= 1, or not a number
    if (detail::isNaN(x)) {
      return x;
    }
    if (((ix - 0x3FF00000U) | detail::low(x)) == 0U) { // asin(+-1) = +-pi/2
      return x * kPio2Hi + x * kPio2Lo;
    }
    return detail::quietNaN();
  }
  if (ix < 0x3FE00000U) {   // |x| < 1/2
    if (ix < 0x3E500000U) { // |x| < 2^-26: asin x is x
      return x;
    }
    const double t = x * x;
    const double p = t * (kPS0 + t * (kPS1 + t * (kPS2 + t * (kPS3 + t * (kPS4 + t * kPS5)))));
    const double q = 1.0 + t * (kQS1 + t * (kQS2 + t * (kQS3 + t * kQS4)));
    return x + x * (p / q);
  }
  // 1/2 <= |x| < 1
  const double w = 1.0 - detail::absOf(x);
  double t = w * 0.5;
  double p = t * (kPS0 + t * (kPS1 + t * (kPS2 + t * (kPS3 + t * (kPS4 + t * kPS5)))));
  double q = 1.0 + t * (kQS1 + t * (kQS2 + t * (kQS3 + t * kQS4)));
  const double s = o2::gpu::CAMath::Sqrtd(t);
  if (ix >= 0x3FEF3333U) { // |x| > 0.975
    t = kPio2Hi - (2.0 * (s + s * (p / q)) - kPio2Lo);
  } else {
    // s split as its high word and the rest, c, for the extra bits.
    const double high = detail::withLowZero(s);
    const double c = (t - high * high) / (s + high);
    const double r = p / q;
    p = 2.0 * s * r - (kPio2Lo - 2.0 * c);
    q = kPio4Hi - 2.0 * high;
    t = kPio4Hi - (p - q);
  }
  return detail::negative(x) ? -t : t;
}

/// fdlibm's e_log.c: x = 2^k (1 + f), sqrt(2)/2 < 1 + f < sqrt(2), and
/// log(1 + f) from s = f / (2 + f).
GPUhdi() GPUAlwaysInline() double log(double x)
{
  constexpr double kLn2Hi = 6.93147180369123816490e-01; // NOLINT(modernize-use-std-numbers): ln 2 split, not ln 2
  constexpr double kLn2Lo = 1.90821492927058770002e-10;
  constexpr double kTwo54 = 1.80143985094819840000e+16;
  constexpr double kLg1 = 6.666666666666735130e-01;
  constexpr double kLg2 = 3.999999999940941908e-01;
  constexpr double kLg3 = 2.857142874366239149e-01;
  constexpr double kLg4 = 2.222219843214978396e-01;
  constexpr double kLg5 = 1.818357216161805012e-01;
  constexpr double kLg6 = 1.531383769920937332e-01;
  constexpr double kLg7 = 1.479819860511658591e-01;

  if (detail::isNaN(x)) {
    return x;
  }
  auto hx = static_cast<std::int32_t>(detail::high(x));
  std::int32_t k = 0;
  if (hx < 0x00100000) { // x < 2^-1022: zero, negative or subnormal
    if (((static_cast<std::uint32_t>(hx) & 0x7FFFFFFFU) | detail::low(x)) == 0U) {
      return -detail::infinity();
    }
    if (hx < 0) {
      return detail::quietNaN();
    }
    k -= 54;
    x *= kTwo54;
    hx = static_cast<std::int32_t>(detail::high(x));
  }
  if (hx >= 0x7FF00000) {
    return x; // +infinity
  }
  k += (hx >> 20) - 1023;
  hx &= 0x000FFFFF;
  const std::int32_t i = (hx + 0x95F64) & 0x100000;
  x = detail::withHigh(x, static_cast<std::uint32_t>(hx | (i ^ 0x3FF00000))); // x or x/2 into [sqrt(2)/2, sqrt(2))
  k += i >> 20;
  const double f = x - 1.0;
  const auto dk = static_cast<double>(k);
  if ((0x000FFFFF & (2 + hx)) < 3) { // -2^-20 <= f < 2^-20
    if (detail::isZero(f)) {
      return k == 0 ? 0. : dk * kLn2Hi + dk * kLn2Lo;
    }
    const double r = f * f * (0.5 - 0.33333333333333333 * f);
    return k == 0 ? f - r : dk * kLn2Hi - ((r - dk * kLn2Lo) - f);
  }
  const double s = f / (2.0 + f);
  const double z = s * s;
  const double w = z * z;
  const double t1 = w * (kLg2 + w * (kLg4 + w * kLg6));
  const double t2 = z * (kLg1 + w * (kLg3 + w * (kLg5 + w * kLg7)));
  const double r = t2 + t1;
  if (((hx - 0x6147A) | (0x6B851 - hx)) > 0) {
    const double hfsq = 0.5 * f * f;
    return k == 0 ? f - (hfsq - s * (hfsq + r)) : dk * kLn2Hi - ((hfsq - (s * (hfsq + r) + dk * kLn2Lo)) - f);
  }
  return k == 0 ? f - s * (f - r) : dk * kLn2Hi - ((s * (f - r) - dk * kLn2Lo) - f);
}

} // namespace o2::gtrk::portable

#endif // GTRK_BASE_PORTABLEMATH_H_
