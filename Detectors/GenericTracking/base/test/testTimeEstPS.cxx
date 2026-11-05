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

/// \file testTimeEstPS.cxx
/// \brief The time box in picoseconds and its edges rounded from BC
///        (GenericTracking/Base/TimeEstPS.h).

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK time box
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/TimeEstPS.h"

#include "CommonConstants/LHCConstants.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace o2::gtrk
{

namespace
{

/// The BC of a 48-orbit timeframe, beyond any the boxes hold.
constexpr std::int64_t kBCs = std::int64_t{48} * o2::constants::lhc::LHCMaxBunches;
constexpr std::int64_t kFrameBC = 198;

} // namespace

BOOST_AUTO_TEST_CASE(BunchCrossingsArePicosecondsOfO2sSpacingRoundedDown)
{
  for (std::int64_t bc = 0; bc <= kBCs; bc += 997) {
    const double exactPS = static_cast<double>(bc) * o2::constants::lhc::LHCBunchSpacingNS * 1000.;
    const double below = exactPS - static_cast<double>(psOfBC(bc));
    BOOST_TEST(below > -1e-3);
    BOOST_TEST(below < 1. + 1e-3);
  }
  BOOST_TEST(psOfBC(0) == 0);
  BOOST_TEST(psOfBC(kFrameBC) == 4'940'255);
  BOOST_TEST(psOfBC(2 * kFrameBC) == 9'880'510);
}

BOOST_AUTO_TEST_CASE(PicosecondsGiveTheirBunchCrossingBack)
{
  for (std::int64_t bc = 0; bc <= kBCs; ++bc) {
    if (bcOfPS(psOfBC(bc)) != bc) {
      BOOST_FAIL("bcOfPS(psOfBC(" << bc << ")) is " << bcOfPS(psOfBC(bc)));
    }
  }
}

BOOST_AUTO_TEST_CASE(FramesAreEquallyWideAndOverlapOnlyWithAMargin)
{
  const std::int64_t frames = (std::int64_t{32} * o2::constants::lhc::LHCMaxBunches) / kFrameBC;
  for (std::int64_t k = 0; k + 2 < frames; ++k) {
    const auto start = 64 + (k * kFrameBC); // a frame bias as data has one
    const auto box = frameBox(start, start + kFrameBC, 0);
    const auto next = frameBox(start + kFrameBC, start + (2 * kFrameBC), 0);
    BOOST_REQUIRE(box.width == static_cast<std::uint32_t>(psOfBC(kFrameBC)));
    BOOST_REQUIRE(box.upper() <= next.lower());
    BOOST_REQUIRE(next.lower() - box.upper() <= 1U);
    BOOST_REQUIRE(!box.isCompatible(next));
    const auto wide = frameBox(start, start + kFrameBC, 30);
    BOOST_REQUIRE(wide.width == static_cast<std::uint32_t>(psOfBC(kFrameBC + 60)));
    BOOST_REQUIRE(wide.isCompatible(frameBox(start + kFrameBC, start + (2 * kFrameBC), 30)));
    BOOST_REQUIRE(!wide.isCompatible(frameBox(start + (2 * kFrameBC), start + (3 * kFrameBC), 30)));
  }
}

/// Two boxes overlap exactly when their BC intervals do, whatever their
/// lengths -- staggered frames, margins -- and wherever they start.
BOOST_AUTO_TEST_CASE(BoxesOverlapExactlyWhenTheirBunchCrossingsDo)
{
  constexpr std::array<std::int64_t, 8> kLengths{1, 2, 60, 197, 198, 228, 258, 891};
  std::int64_t checked = 0;
  for (std::int64_t a = 0; a < 1200; a += 3) {
    for (const auto la : kLengths) {
      for (std::int64_t b = a - 900 < 0 ? 0 : a - 900; b < a + 900; ++b) {
        for (const auto lb : kLengths) {
          const bool bc = a < b + lb && b < a + la;
          if (frameBox(a, a + la, 0).isCompatible(frameBox(b, b + lb, 0)) != bc) {
            BOOST_FAIL("[" << a << ", " << a + la << ") and [" << b << ", " << b + lb << ") BC");
          }
          ++checked;
        }
      }
    }
  }
  BOOST_TEST(checked > 1'000'000);
}

/// A box is one of a layer's frames exactly when `frameBox` made it from
/// one of them with no margin.
BOOST_AUTO_TEST_CASE(ABoxIsAFrameOnlyAsFrameBoxMakesIt)
{
  constexpr FrameTiming kFrames{.lengthBC = 22, .firstStartBC = 64};
  for (std::int64_t frame = 0; frame < 5184; ++frame) {
    const auto start = 64 + (frame * 22);
    BOOST_REQUIRE(kFrames.isFrame(frameBox(start, start + 22, 0)));
  }
  BOOST_TEST(!kFrames.isFrame(frameBox(64 + 22, 64 + 44, 1))); // widened
  BOOST_TEST(!kFrames.isFrame(frameBox(64 + 11, 64 + 33, 0))); // between frame starts
  BOOST_TEST(!kFrames.isFrame(frameBox(64 - 22, 64, 0)));      // before the first frame
  BOOST_TEST(!kFrames.isFrame(frameBox(64, 64 + 81, 0)));      // another length
  BOOST_TEST(!kFrames.isFrame({.start = 1U + static_cast<std::uint32_t>(psOfBC(64)), .width = static_cast<std::uint32_t>(psOfBC(22))}));
  BOOST_TEST(!FrameTiming{}.isFrame(frameBox(0, 22, 0))); // no frames
}

BOOST_AUTO_TEST_CASE(TheFirstFrameIsClampedAtTheTimeframesStart)
{
  const auto box = frameBox(0, kFrameBC, 30);
  BOOST_TEST(box.lower() == 0U);
  BOOST_TEST(box.width == static_cast<std::uint32_t>(psOfBC(kFrameBC + 30)));
}

BOOST_AUTO_TEST_CASE(TheIntersectionNarrowsAndTheEndSaturates)
{
  TimeEstPS box{.start = 0U, .width = 10U};
  BOOST_TEST(box.isCompatible(TimeEstPS{.start = 5U, .width = 10U}));
  box += TimeEstPS{.start = 5U, .width = 10U};
  box += TimeEstPS{.start = 9U, .width = 10U};
  BOOST_TEST(box.lower() == 9U);
  BOOST_TEST(box.upper() == 10U);
  BOOST_TEST(!box.isCompatible(TimeEstPS{.start = 10U, .width = 5U}));

  const TimeEstPS late{.start = std::numeric_limits<std::uint32_t>::max() - 3U, .width = 10U};
  BOOST_TEST(late.upper() == std::numeric_limits<std::uint32_t>::max());
}

} // namespace o2::gtrk
