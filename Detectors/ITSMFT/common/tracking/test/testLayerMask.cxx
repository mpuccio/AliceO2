// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE ITSMFT LayerMask
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "ITSMFTTracking/SurfaceMask.h"

/// Baseline characterization of o2::itsmft::tracking::LayerMask (Wave 0).
/// These tests record the *current* behavior of the shared 16-bit mask
/// ahead of the surface-mask migration described in Architecture.md (D002).
/// They do not encode the future wider-mask design.

using o2::itsmft::tracking::LayerMask;

BOOST_AUTO_TEST_CASE(layermask_default_and_basic_construction)
{
  LayerMask empty{};
  BOOST_CHECK(empty.empty());
  BOOST_CHECK_EQUAL(empty.count(), 0);
  BOOST_CHECK_EQUAL(empty.value(), 0);

  LayerMask single{static_cast<uint16_t>(1 << 3)};
  BOOST_CHECK(!single.empty());
  BOOST_CHECK(single.has(3));
  BOOST_CHECK(!single.has(2));
  BOOST_CHECK_EQUAL(single.count(), 1);

  LayerMask triple{0, 2, 4};
  BOOST_CHECK_EQUAL(triple.value(), 0b10101);
  BOOST_CHECK(triple.has(0));
  BOOST_CHECK(triple.has(2));
  BOOST_CHECK(triple.has(4));
  BOOST_CHECK(!triple.has(1));
  BOOST_CHECK(!triple.has(3));
}

BOOST_AUTO_TEST_CASE(layermask_set_and_conversion)
{
  LayerMask mask{};
  mask.set(1);
  mask.set(5);
  BOOST_CHECK_EQUAL(static_cast<uint16_t>(mask), 0b100010);
  BOOST_CHECK_EQUAL(mask.count(), 2);
}

BOOST_AUTO_TEST_CASE(layermask_bitwise_operators)
{
  const LayerMask a{0x77}; // 0111 0111
  const LayerMask b{0x0f}; // 0000 1111

  BOOST_CHECK_EQUAL((a & b).value(), 0x07);
  BOOST_CHECK_EQUAL((a | b).value(), 0x7f);
  BOOST_CHECK_EQUAL((~LayerMask{0x00ff}).value(), uint32_t{0xffffff00});

  LayerMask c{a};
  c &= b;
  BOOST_CHECK_EQUAL(c.value(), 0x07);

  LayerMask d{a};
  d |= b;
  BOOST_CHECK_EQUAL(d.value(), 0x7f);
}

BOOST_AUTO_TEST_CASE(layermask_holes_and_length)
{
  const LayerMask layer3Hole{0x77}; // layers 0,1,2,4,5,6 (7-layer ITS-like mask with a hole at 3)
  BOOST_CHECK_EQUAL(layer3Hole.count(), 6);
  BOOST_CHECK_EQUAL(layer3Hole.length(), 7);
  BOOST_CHECK_EQUAL(layer3Hole.holeMask().value(), 0x08);
  BOOST_CHECK(layer3Hole.isAllowed(1, 0x08));
  BOOST_CHECK(!layer3Hole.isAllowed(0, 0x08));
  BOOST_CHECK(!layer3Hole.isAllowed(1, 0x10)); // hole present but at a non-allowed position

  const LayerMask missingLeadingLayer0{0x7e}; // layers 1..6
  BOOST_CHECK_EQUAL(missingLeadingLayer0.count(), 6);
  BOOST_CHECK_EQUAL(missingLeadingLayer0.length(), 6);
  BOOST_CHECK_EQUAL(missingLeadingLayer0.holeMask().value(), 0x00);

  const LayerMask missingTrailingLayer6{0x3f}; // layers 0..5
  BOOST_CHECK_EQUAL(missingTrailingLayer6.count(), 6);
  BOOST_CHECK_EQUAL(missingTrailingLayer6.length(), 6);
  BOOST_CHECK_EQUAL(missingTrailingLayer6.holeMask().value(), 0x00);
}

BOOST_AUTO_TEST_CASE(layermask_holes_beyond_seven_layers_mft_like)
{
  // 10-layer MFT-like mask with a hole at layer 4: layers 0,1,2,3,5,6,7,8,9
  const LayerMask mftHole{static_cast<uint16_t>(0x3ff & ~(1 << 4))};
  BOOST_CHECK_EQUAL(mftHole.count(), 9);
  BOOST_CHECK_EQUAL(mftHole.length(), 10);
  BOOST_CHECK_EQUAL(mftHole.holeMask().value(), 1 << 4);
  BOOST_CHECK(mftHole.isAllowed(1, 1 << 4));
  BOOST_CHECK(!mftHole.isAllowed(0, 1 << 4));
}

BOOST_AUTO_TEST_CASE(layermask_first_last_empty_sentinel)
{
  const LayerMask empty{};
  BOOST_CHECK_EQUAL(empty.first(), o2::its::constants::UnusedIndex);
  BOOST_CHECK_EQUAL(empty.last(), o2::its::constants::UnusedIndex);
  BOOST_CHECK_EQUAL(empty.length(), 0);

  const LayerMask single{static_cast<uint16_t>(1 << 9)}; // topmost MFT layer
  BOOST_CHECK_EQUAL(single.first(), 9);
  BOOST_CHECK_EQUAL(single.last(), 9);
  BOOST_CHECK_EQUAL(single.length(), 1);
}

BOOST_AUTO_TEST_CASE(layermask_slot)
{
  const LayerMask mask{0, 2, 5}; // layers 0,2,5 populated
  BOOST_CHECK_EQUAL(mask.slot(0), 0);
  BOOST_CHECK_EQUAL(mask.slot(2), 1);
  BOOST_CHECK_EQUAL(mask.slot(5), 2);
  BOOST_CHECK_EQUAL(mask.slot(1), o2::its::constants::UnusedIndex);
  BOOST_CHECK_EQUAL(mask.slot(4), o2::its::constants::UnusedIndex);
}

BOOST_AUTO_TEST_CASE(layermask_span)
{
  BOOST_CHECK_EQUAL(LayerMask::span(0, 6).value(), 0x7f);   // full 7-layer ITS span
  BOOST_CHECK_EQUAL(LayerMask::span(0, 9).value(), 0x3ff);  // full 10-layer MFT span
  BOOST_CHECK_EQUAL(LayerMask::span(2, 4).value(), 0b11100);
  BOOST_CHECK_EQUAL(LayerMask::span(3, 3).value(), 1 << 3); // single-layer span
  BOOST_CHECK_EQUAL(LayerMask::span(4, 2).value(), 0);      // inverted range -> empty
}

BOOST_AUTO_TEST_CASE(layermask_skipped)
{
  BOOST_CHECK_EQUAL(LayerMask::skipped(0, 1).value(), 0);           // adjacent layers, nothing skipped
  BOOST_CHECK_EQUAL(LayerMask::skipped(0, 0).value(), 0);           // degenerate range
  BOOST_CHECK_EQUAL(LayerMask::skipped(1, 3).value(), 1 << 2);      // one layer skipped
  BOOST_CHECK_EQUAL(LayerMask::skipped(0, 6).value(), 0b0111110);   // skip everything strictly between 0 and 6
}

BOOST_AUTO_TEST_CASE(layermask_is_subset_of)
{
  const LayerMask mask{0x0f};
  BOOST_CHECK(mask.isSubsetOf(0xff));
  BOOST_CHECK(mask.isSubsetOf(0x0f));
  BOOST_CHECK(!mask.isSubsetOf(0x07));
}

BOOST_AUTO_TEST_CASE(layermask_as_string_is_stable_width)
{
  const LayerMask mask{0x7f};
  BOOST_CHECK_EQUAL(mask.asString(), "00000000000000000000000001111111");
}

BOOST_AUTO_TEST_CASE(layermask_covers_the_complete_surface_rank_domain)
{
  static_assert(sizeof(LayerMask) == sizeof(uint32_t));

  LayerMask mask{};
  mask.set(15);
  BOOST_CHECK(mask.has(15));
  BOOST_CHECK_EQUAL(mask.value(), uint32_t{1} << 15);

  LayerMask combined{};
  combined.set(16);
  BOOST_CHECK(combined.has(16));
  BOOST_CHECK_EQUAL(combined.value(), uint32_t{1} << 16);
}
