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

/// \file testTimeFrameScratch.cxx
/// \brief The slice's source rows, numbered, and the scratch that holds them.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TimeFrameScratch
#include <boost/test/unit_test.hpp>

#include "ToyBarrel3.h"
#include "GenericTracking/Tracking/LayerRowMap.h"
#include "GenericTracking/Tracking/TimeFrameScratch.h"

#include <array>
#include <memory>

namespace o2::gtrk
{

namespace
{

using Det = detectors::ToyBarrel3;
using Scratch = TimeFrameScratch<Det>;
using Map = LayerRowMap<Det>;

constexpr std::array<MeasurementRange, Map::kLayerCount> kRows{{{.begin = 10, .end = 13}, {.begin = 20, .end = 22}, {.begin = 30, .end = 34}}};
constexpr std::array<MeasurementRange, Map::kLayerCount> kHalo{{{.begin = 8, .end = 15}, {.begin = 18, .end = 24}, {.begin = 28, .end = 36}}};

} // namespace

BOOST_AUTO_TEST_CASE(PresentedRowsHaveOneDenseOrder)
{
  const Map map{kRows, {true, true, true}};
  BOOST_TEST(map.size() == 9U);
  BOOST_TEST(map.offsets[1] == 3U);
  BOOST_TEST(map.offsets[2] == 5U);
  BOOST_TEST(map.element(0, 12U) == 2U);
  BOOST_TEST(map.element(1, 20U) == 3U);
  BOOST_TEST(map.element(2, 33U) == 8U);
  for (Map::Index element = 0; element < map.size(); ++element) {
    const auto at = map.row(element);
    BOOST_TEST(map.element(at.layer, at.row) == element);
  }
  BOOST_TEST(!map.contains(0, 13U));
  BOOST_TEST(!map.contains(3, 10U));
}

BOOST_AUTO_TEST_CASE(AnEmptyOrUnpresentedLayerIsSkipped)
{
  // The middle layer contributes nothing, so its offset repeats.
  const Map map{kRows, {true, false, true}};
  BOOST_TEST(map.size() == 7U);
  BOOST_TEST(!map.contains(1, 20U));
  BOOST_TEST(map.row(2U).layer == 0U);
  BOOST_TEST(map.row(3U).layer == 2U);
  BOOST_TEST(map.row(3U).row == 30U);
  BOOST_TEST(map.row(6U).row == 33U);
}

BOOST_AUTO_TEST_CASE(TrackletSourcesAreTheLayersWithAnOutgoingLink)
{
  // ToyBarrel3 links 0 -> 1 -> 2; the outermost layer starts nothing.
  BOOST_TEST(Scratch::kSourceLayers[0]);
  BOOST_TEST(Scratch::kSourceLayers[1]);
  BOOST_TEST(!Scratch::kSourceLayers[2]);

  // Every source layer contributes its halo rows; the tracklet stage's gate
  // decides which it works.
  Scratch scratch;
  scratch.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  scratch.beginSlice(Scratch::makeSources(kHalo));
  BOOST_TEST(scratch.sources().size() == 13U);
  BOOST_TEST(scratch.tracklets().elementCount() == 13U);
  BOOST_TEST(scratch.tracklets().offsets().size() == 14U);
  BOOST_TEST(scratch.sources().contains(0, 9U));
  BOOST_TEST(scratch.sources().contains(1, 18U));
  BOOST_TEST(!scratch.sources().contains(2, 30U));
}

} // namespace o2::gtrk
