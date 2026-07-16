// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT SurfaceTiming
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <limits>

#include "CommonDataFormat/InteractionRecord.h"
#include "ITSMFTTracking/SurfaceTiming.h"

using namespace o2::itsmft::tracking;

BOOST_AUTO_TEST_CASE(ROFIntervalBCViewsAreDeviceFriendly)
{
  static_assert(std::is_standard_layout_v<ROFIntervalBC>);
  static_assert(std::is_trivially_copyable_v<ROFIntervalBC>);
  static_assert(std::is_standard_layout_v<ROFTimingConfig>);
  static_assert(std::is_trivially_copyable_v<ROFTimingConfig>);
  BOOST_CHECK(std::is_standard_layout_v<ROFIntervalBC>);
  BOOST_CHECK(std::is_trivially_copyable_v<ROFIntervalBC>);
}

BOOST_AUTO_TEST_CASE(BeginEndFollowLayerTimingSignConvention)
{
  // Mirrors o2::its::LayerTiming::getROFStartInBC/getROFEndInBC: start is the
  // ROF's own BC plus delay plus bias; end is start plus the readout length.
  const o2::InteractionRecord origin{100, 0};
  const o2::InteractionRecord rofIR{150, 0}; // 50 BC after origin
  const ROFTimingConfig cfg{/*rofLength*/ 40, /*rofDelay*/ 5, /*rofBias*/ -2, /*rofAddTimeErr*/ 3};

  const auto built = computeROFIntervalBC(rofIR, origin, cfg, 7);
  BOOST_REQUIRE(built.ok());
  BOOST_CHECK_EQUAL(built.interval.begin, 53); // 50 + 5 - 2
  BOOST_CHECK_EQUAL(built.interval.end, 93);   // begin + rofLength
  BOOST_CHECK_EQUAL(built.interval.sourceROF, 7u);
  BOOST_CHECK(built.interval.isValid());
}

BOOST_AUTO_TEST_CASE(NegativeTFRelativeBCIsLegal)
{
  const o2::InteractionRecord origin{200, 0};
  const o2::InteractionRecord rofIR{50, 0}; // 150 BC before origin
  const ROFTimingConfig cfg{40, 0, 0, 0};

  const auto built = computeROFIntervalBC(rofIR, origin, cfg, 0);
  BOOST_REQUIRE(built.ok());
  BOOST_CHECK_EQUAL(built.interval.begin, -150);
  BOOST_CHECK_EQUAL(built.interval.end, -110);
  BOOST_CHECK(built.interval.isValid());
}

BOOST_AUTO_TEST_CASE(InvalidROFLengthIsRejected)
{
  const o2::InteractionRecord origin{0, 0};
  const ROFTimingConfig cfg{0, 0, 0, 0};
  const auto built = computeROFIntervalBC(origin, origin, cfg, 0);
  BOOST_CHECK(!built.ok());
  BOOST_CHECK(built.error == TimingBuildError::InvalidROFLength);
}

BOOST_AUTO_TEST_CASE(OverflowIsDetectedAndChecked)
{
  const o2::InteractionRecord origin{0, 0};
  const o2::InteractionRecord rofIR{0, 0};
  ROFTimingConfig cfg{1, std::numeric_limits<TFBC>::max(), 1, 0};
  const auto built = computeROFIntervalBC(rofIR, origin, cfg, 0);
  BOOST_CHECK(!built.ok());
  BOOST_CHECK(built.error == TimingBuildError::Overflow);
}

BOOST_AUTO_TEST_CASE(WidenAppliesSymmetricMarginOnDemandOnly)
{
  const ROFIntervalBC interval{10, 20, 3, 0};
  const auto widened = widen(interval, 5);
  BOOST_CHECK_EQUAL(widened.begin, 5);
  BOOST_CHECK_EQUAL(widened.end, 25);
  // The base interval itself is never mutated by widen().
  BOOST_CHECK_EQUAL(interval.begin, 10);
  BOOST_CHECK_EQUAL(interval.end, 20);
}

BOOST_AUTO_TEST_CASE(IntersectionIsHalfOpenAndIgnoresROFOrdinal)
{
  const ROFIntervalBC a{0, 10, 5, 0};
  const ROFIntervalBC touching{10, 20, 5, 0}; // same sourceROF as `a`, adjacent
  const ROFIntervalBC overlapping{9, 15, 99, 0}; // different sourceROF, overlaps
  const ROFIntervalBC disjoint{100, 110, 5, 0};

  BOOST_CHECK(!intersects(a, touching));    // half-open: touching does not intersect
  BOOST_CHECK(intersects(a, overlapping));  // overlap decided by time, not ROF equality
  BOOST_CHECK(!intersects(a, disjoint));
  BOOST_CHECK(intersects(a, a));
}

BOOST_AUTO_TEST_CASE(DefaultIntervalIsInvalidSentinel)
{
  constexpr ROFIntervalBC interval{};
  BOOST_CHECK_EQUAL(interval.sourceROF, std::numeric_limits<uint32_t>::max());
  BOOST_CHECK(interval.isValid()); // begin == end == 0
  BOOST_CHECK_EQUAL(interval.length(), 0);
}
