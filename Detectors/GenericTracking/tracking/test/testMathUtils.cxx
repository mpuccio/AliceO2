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

/// \file testMathUtils.cxx
/// \brief Exercises scalar helpers shared by the tracking code.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK MathUtils
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/MathUtils.h"

#include "CommonConstants/MathConstants.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <type_traits>

namespace o2::gtrk
{

BOOST_AUTO_TEST_CASE(TransverseArcHasTheStraightAndCurvedLimits)
{
  static_assert(std::is_same_v<decltype(transverseArcAtRadius(1.F, 2.F)), float>);
  static_assert(std::is_same_v<decltype(transverseArcAtRadius(1., 2.)), double>);

  BOOST_TEST(transverseArcAtRadius(0., 4.) == 4.);
  BOOST_TEST(transverseArcAtRadius(-1., 4.) == 4.);

  constexpr double curvature = 0.2;
  constexpr double radius = 4.;
  const double expected = 2. * std::asin(0.5 * curvature * radius) / curvature;
  BOOST_TEST(transverseArcAtRadius(curvature, radius) == expected, boost::test_tools::tolerance(1e-14));
}

/// The quotient from one asin with its derivative is, bit for bit, the
/// quotient alone, on both sides of the
/// series' switch and at the ends.
BOOST_AUTO_TEST_CASE(AsinOverArgumentWithDerivativeIsTheQuotient)
{
  int mismatches = 0;
  for (int i = -20000; i <= 20000; ++i) {
    const float x = static_cast<float>(i) / 20000.F;
    const auto both = detail::asinOverArgumentWithDerivative(x);
    mismatches += std::bit_cast<std::uint32_t>(both.value) != std::bit_cast<std::uint32_t>(detail::asinOverArgument(x)) ? 1 : 0;
  }
  for (const float x : {0.05F, -0.05F, std::nextafter(0.05F, 0.F), 1.F, -1.F}) {
    const auto both = detail::asinOverArgumentWithDerivative(x);
    mismatches += std::bit_cast<std::uint32_t>(both.value) != std::bit_cast<std::uint32_t>(detail::asinOverArgument(x)) ? 1 : 0;
  }
  BOOST_TEST(mismatches == 0);
}

BOOST_AUTO_TEST_CASE(ChiSquareQuantilesMatchKnownValues)
{
  // Three degrees of freedom at the 68.27 percent and 99.73 percent
  // two-sided Gaussian confidence levels.
  BOOST_TEST(chiSquareQuantile(3, 1.) == 3.5267, boost::test_tools::tolerance(1.e-4));
  BOOST_TEST(chiSquareQuantile(3, 3.) == 14.1564, boost::test_tools::tolerance(1.e-4));
  // One degree of freedom is n^2 exactly; twenty was checked by numerical
  // integration of the density.
  BOOST_TEST(chiSquareQuantile(1, 3.) == 9., boost::test_tools::tolerance(1.e-6));
  BOOST_TEST(chiSquareQuantile(20, 3.) == 42.0802, boost::test_tools::tolerance(1.e-4));
}

} // namespace o2::gtrk
