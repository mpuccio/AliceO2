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

/// \file testPortableMath.cxx
/// \brief The deterministic build's sin, cos, atan, atan2, asin and log
///        (PortableMath.h) against long double: fdlibm's accuracy, and the
///        special values. That they give the device the same bits is
///        testDeviceRuntime.cxx's.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK PortableMath
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Base/PortableMath.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <vector>

namespace portable = o2::gtrk::portable;

namespace
{

/// |value - exact| in units of the last place of \a value.
double ulpError(double value, long double exact)
{
  const double ulp = std::nextafter(std::abs(value), std::numeric_limits<double>::infinity()) - std::abs(value);
  return static_cast<double>(std::abs(static_cast<long double>(value) - exact) / static_cast<long double>(ulp));
}

std::vector<double> uniform(double low, double high, std::size_t count, std::uint64_t seed)
{
  std::mt19937_64 generator{seed};
  std::uniform_real_distribution<double> distribution{low, high};
  std::vector<double> values(count);
  for (auto& value : values) {
    value = distribution(generator);
  }
  return values;
}

/// Arguments within a few ulp of multiples of pi/2, where the reduction
/// cancels most: the 2nd and 3rd rounds of remPio2.
std::vector<double> nearMultiplesOfHalfPi()
{
  std::vector<double> values;
  for (int n = -2000; n <= 2000; ++n) {
    double x = n * (std::numbers::pi / 2.);
    for (int step = 0; step < 4; ++step) {
      values.push_back(x);
      x = std::nextafter(x, std::numeric_limits<double>::infinity());
    }
  }
  values.push_back(6381956970095103.0 * 0x1p-50); // 5.67e0, a hard case of Kahan's
  values.push_back(1.e6);
  values.push_back(-1.6e6);
  return values;
}

bool sameBits(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

} // namespace

BOOST_AUTO_TEST_CASE(SinAndCosWithinAnUlp)
{
  auto values = uniform(-50., 50., 1'000'000U, 1U);
  const auto small = uniform(-1.e-3, 1.e-3, 100'000U, 2U);
  const auto hard = nearMultiplesOfHalfPi();
  values.insert(values.end(), small.begin(), small.end());
  values.insert(values.end(), hard.begin(), hard.end());
  double worstSin = 0.;
  double worstCos = 0.;
  for (const double x : values) {
    worstSin = std::max(worstSin, ulpError(portable::sin(x), std::sin(static_cast<long double>(x))));
    worstCos = std::max(worstCos, ulpError(portable::cos(x), std::cos(static_cast<long double>(x))));
  }
  BOOST_TEST_MESSAGE("sin " << worstSin << " ulp, cos " << worstCos << " ulp");
  BOOST_TEST(worstSin < 1.);
  BOOST_TEST(worstCos < 1.);
}

BOOST_AUTO_TEST_CASE(AtanAndAsinWithinAnUlp)
{
  auto values = uniform(-4., 4., 1'000'000U, 3U);
  const auto wide = uniform(-1.e8, 1.e8, 100'000U, 4U);
  values.insert(values.end(), wide.begin(), wide.end());
  double worstAtan = 0.;
  for (const double x : values) {
    worstAtan = std::max(worstAtan, ulpError(portable::atan(x), std::atan(static_cast<long double>(x))));
  }
  double worstAsin = 0.;
  for (const double x : uniform(-1., 1., 1'000'000U, 5U)) {
    worstAsin = std::max(worstAsin, ulpError(portable::asin(x), std::asin(static_cast<long double>(x))));
  }
  for (const double x : {0.49999999999999994, 0.5, 0.97499999999999998, 0.975, 0.99999999999999989}) {
    worstAsin = std::max(worstAsin, ulpError(portable::asin(x), std::asin(static_cast<long double>(x))));
    worstAsin = std::max(worstAsin, ulpError(portable::asin(-x), std::asin(static_cast<long double>(-x))));
  }
  BOOST_TEST_MESSAGE("atan " << worstAtan << " ulp, asin " << worstAsin << " ulp");
  BOOST_TEST(worstAtan < 1.);
  BOOST_TEST(worstAsin < 1.);
}

BOOST_AUTO_TEST_CASE(Atan2WithinTwoUlpInEveryQuadrant)
{
  const auto ys = uniform(-3., 3., 1'000'000U, 6U);
  const auto xs = uniform(-3., 3., 1'000'000U, 7U);
  double worst = 0.;
  for (std::size_t i = 0; i < ys.size(); ++i) {
    worst = std::max(worst, ulpError(portable::atan2(ys[i], xs[i]), std::atan2(static_cast<long double>(ys[i]), static_cast<long double>(xs[i]))));
  }
  // The propagator's turn: sine and cosine of up to about 100 degrees.
  const auto sines = uniform(-1., 1., 1'000'000U, 8U);
  const auto cosines = uniform(-0.2, 1., 1'000'000U, 9U);
  for (std::size_t i = 0; i < sines.size(); ++i) {
    worst =
      std::max(worst, ulpError(portable::atan2(sines[i], cosines[i]), std::atan2(static_cast<long double>(sines[i]), static_cast<long double>(cosines[i]))));
  }
  BOOST_TEST_MESSAGE("atan2 " << worst << " ulp");
  BOOST_TEST(worst < 2.);
}

BOOST_AUTO_TEST_CASE(LogWithinAnUlp)
{
  auto values = uniform(0., 4., 1'000'000U, 10U);
  const auto nearOne = uniform(1. - 1.e-6, 1. + 1.e-6, 100'000U, 11U);
  values.insert(values.end(), nearOne.begin(), nearOne.end());
  // No subnormals: the test links O2's libO2GPUDataTypes.so, which sets the
  // host's DAZ at load (crtfastmath.o), and under DAZ one reads as zero.
  for (const double x : {1.e-300, std::numeric_limits<double>::min(), 1.e-200, 1.e300, std::numeric_limits<double>::max(), 2., 0.5}) {
    values.push_back(x);
  }
  double worst = 0.;
  for (const double x : values) {
    if (x > 0.) {
      worst = std::max(worst, ulpError(portable::log(x), std::log(static_cast<long double>(x))));
    }
  }
  BOOST_TEST_MESSAGE("log " << worst << " ulp");
  BOOST_TEST(worst < 1.);
}

BOOST_AUTO_TEST_CASE(SpecialValues)
{
  BOOST_TEST(sameBits(portable::sin(0.), 0.));
  BOOST_TEST(sameBits(portable::sin(-0.), -0.));
  BOOST_TEST(sameBits(portable::cos(0.), 1.));
  BOOST_TEST(sameBits(portable::cos(-0.), 1.));
  BOOST_TEST(std::isnan(portable::sin(kInf)));
  BOOST_TEST(std::isnan(portable::cos(-kInf)));
  BOOST_TEST(std::isnan(portable::sin(kNaN)));

  BOOST_TEST(sameBits(portable::atan(-0.), -0.));
  BOOST_TEST(sameBits(portable::atan(kInf), std::numbers::pi / 2.));
  BOOST_TEST(sameBits(portable::atan(-kInf), -std::numbers::pi / 2.));
  BOOST_TEST(std::isnan(portable::atan(kNaN)));

  BOOST_TEST(sameBits(portable::atan2(0., 1.), 0.));
  BOOST_TEST(sameBits(portable::atan2(-0., 1.), -0.));
  BOOST_TEST(sameBits(portable::atan2(0., 0.), 0.));
  BOOST_TEST(sameBits(portable::atan2(-0., 0.), -0.));
  BOOST_TEST(sameBits(portable::atan2(0., -0.), std::numbers::pi));
  BOOST_TEST(sameBits(portable::atan2(-0., -0.), -std::numbers::pi));
  BOOST_TEST(sameBits(portable::atan2(0., -1.), std::numbers::pi));
  BOOST_TEST(sameBits(portable::atan2(-0., -1.), -std::numbers::pi));
  BOOST_TEST(sameBits(portable::atan2(1., 0.), std::numbers::pi / 2.));
  BOOST_TEST(sameBits(portable::atan2(-1., -0.), -std::numbers::pi / 2.));
  BOOST_TEST(sameBits(portable::atan2(kInf, kInf), std::numbers::pi / 4.));
  BOOST_TEST(sameBits(portable::atan2(-kInf, -kInf), -3. * std::numbers::pi / 4.));
  BOOST_TEST(sameBits(portable::atan2(1., kInf), 0.));
  BOOST_TEST(sameBits(portable::atan2(-1., kInf), -0.));
  BOOST_TEST(sameBits(portable::atan2(1., -kInf), std::numbers::pi));
  BOOST_TEST(sameBits(portable::atan2(-kInf, 1.), -std::numbers::pi / 2.));
  BOOST_TEST(sameBits(portable::atan2(1.e300, 1.e-300), std::numbers::pi / 2.));
  BOOST_TEST(sameBits(portable::atan2(1.e-300, -1.e300), std::numbers::pi));
  BOOST_TEST(sameBits(portable::atan2(0.5, 1.), portable::atan(0.5)));
  BOOST_TEST(std::isnan(portable::atan2(kNaN, 1.)));
  BOOST_TEST(std::isnan(portable::atan2(1., kNaN)));

  BOOST_TEST(sameBits(portable::asin(1.), std::numbers::pi / 2.));
  BOOST_TEST(sameBits(portable::asin(-1.), -std::numbers::pi / 2.));
  BOOST_TEST(sameBits(portable::asin(-0.), -0.));
  BOOST_TEST(std::isnan(portable::asin(1.0000000000000002)));
  BOOST_TEST(std::isnan(portable::asin(-2.)));

  BOOST_TEST(sameBits(portable::log(1.), 0.));
  BOOST_TEST(sameBits(portable::log(0.), -kInf));
  BOOST_TEST(sameBits(portable::log(-0.), -kInf));
  BOOST_TEST(sameBits(portable::log(kInf), kInf));
  BOOST_TEST(std::isnan(portable::log(-1.)));
  BOOST_TEST(std::isnan(portable::log(-kInf)));
  BOOST_TEST(std::isnan(portable::log(kNaN)));
}

/// Invalid arguments give one NaN, the same bits everywhere, not whichever
/// NaN the target's arithmetic makes.
BOOST_AUTO_TEST_CASE(InvalidArgumentsGiveTheQuietNaN)
{
  const auto quiet = std::bit_cast<double>(0x7FF8000000000000ULL);
  BOOST_TEST(sameBits(portable::sin(kInf), quiet));
  BOOST_TEST(sameBits(portable::cos(-kInf), quiet));
  BOOST_TEST(sameBits(portable::asin(2.), quiet));
  BOOST_TEST(sameBits(portable::log(-1.), quiet));
}
