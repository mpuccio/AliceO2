// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Pure deriveUniformROFTimingConfig() tests: no DPLAlpideParam/CCDB/GRPGeom
// singleton dependency, synthetic o2::its::LayerTiming arrays only.

#define BOOST_TEST_MODULE ITSMFT ROFTimingUniformity
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>

#include "ITSMFTTracking/ROFTimingUniformity.h"

using namespace o2::itsmft::tracking;
using o2::its::LayerTiming;

BOOST_AUTO_TEST_CASE(EmptyInputIsNotUniform)
{
  const auto result = deriveUniformROFTimingConfig({});
  BOOST_CHECK(!result.uniform);
}

BOOST_AUTO_TEST_CASE(SingleLayerIsTriviallyUniform)
{
  const std::array<LayerTiming, 1> layers{LayerTiming{.mNROFsTF = 10, .mROFLength = 40, .mROFDelay = 5, .mROFBias = 2, .mROFAddTimeErr = 1}};
  const auto result = deriveUniformROFTimingConfig(layers);
  BOOST_REQUIRE(result.uniform);
  BOOST_CHECK_EQUAL(result.config.rofLength, 40);
  BOOST_CHECK_EQUAL(result.config.rofDelay, 5);
  BOOST_CHECK_EQUAL(result.config.rofBias, 2);
  BOOST_CHECK_EQUAL(result.config.rofAddTimeErr, 1);
}

BOOST_AUTO_TEST_CASE(MatchedPerLayerValuesAreUniform)
{
  // Mirrors real ITS/MFT production defaults: every layer resolves to the
  // same shared global value because DPLAlpideParam's per-layer staggering
  // overrides all default to zero.
  std::array<LayerTiming, 7> layers{};
  for (auto& lt : layers) {
    lt = LayerTiming{.mNROFsTF = 100, .mROFLength = 594, .mROFDelay = 0, .mROFBias = 64, .mROFAddTimeErr = 0};
  }
  const auto result = deriveUniformROFTimingConfig(layers);
  BOOST_REQUIRE(result.uniform);
  BOOST_CHECK_EQUAL(result.config.rofLength, 594);
  BOOST_CHECK_EQUAL(result.config.rofBias, 64);
}

BOOST_AUTO_TEST_CASE(DivergentROFLengthIsRejected)
{
  std::array<LayerTiming, 3> layers{
    LayerTiming{.mROFLength = 40, .mROFDelay = 0, .mROFBias = 0, .mROFAddTimeErr = 0},
    LayerTiming{.mROFLength = 40, .mROFDelay = 0, .mROFBias = 0, .mROFAddTimeErr = 0},
    LayerTiming{.mROFLength = 44, .mROFDelay = 0, .mROFBias = 0, .mROFAddTimeErr = 0}}; // staggered length
  BOOST_CHECK(!deriveUniformROFTimingConfig(layers).uniform);
}

BOOST_AUTO_TEST_CASE(DivergentROFDelayIsRejected)
{
  std::array<LayerTiming, 2> layers{
    LayerTiming{.mROFLength = 40, .mROFDelay = 0, .mROFBias = 0, .mROFAddTimeErr = 0},
    LayerTiming{.mROFLength = 40, .mROFDelay = 3, .mROFBias = 0, .mROFAddTimeErr = 0}};
  BOOST_CHECK(!deriveUniformROFTimingConfig(layers).uniform);
}

BOOST_AUTO_TEST_CASE(DivergentROFBiasIsRejected)
{
  std::array<LayerTiming, 2> layers{
    LayerTiming{.mROFLength = 40, .mROFDelay = 0, .mROFBias = 64, .mROFAddTimeErr = 0},
    LayerTiming{.mROFLength = 40, .mROFDelay = 0, .mROFBias = 60, .mROFAddTimeErr = 0}};
  BOOST_CHECK(!deriveUniformROFTimingConfig(layers).uniform);
}

BOOST_AUTO_TEST_CASE(DivergentAddTimeErrIsRejected)
{
  std::array<LayerTiming, 2> layers{
    LayerTiming{.mROFLength = 40, .mROFDelay = 0, .mROFBias = 0, .mROFAddTimeErr = 0},
    LayerTiming{.mROFLength = 40, .mROFDelay = 0, .mROFBias = 0, .mROFAddTimeErr = 5}};
  BOOST_CHECK(!deriveUniformROFTimingConfig(layers).uniform);
}
