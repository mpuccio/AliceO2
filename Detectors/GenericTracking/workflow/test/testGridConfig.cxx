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

/// \file testGridConfig.cxx
/// \brief The grid's time blocks from the layers' readout frames.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK grid configuration
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Workflow/RecoParam.h"
#include "ToyBarrel3.h"
#include "GenericTracking/Tracking/GridSizing.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "CommonConstants/LHCConstants.h"

#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>

namespace o2::gtrk::workflow
{

namespace
{

constexpr std::int64_t kPartsPerPS = kPSPerBCDenominator;
constexpr std::int64_t kPartsPerBC = kPSPerBCNumerator;

} // namespace

/// With frames of two lengths the blocks are whole frames of the longer,
/// from its first start; each layer's frames travel with the parameters.
BOOST_AUTO_TEST_CASE(TheBlocksFollowTheLongestFrame)
{
  const RecoSettings param{};
  BOOST_REQUIRE(param.blockWidthPS == 10'000'000);
  std::vector<FrameTiming> frames(3U, {.lengthBC = 22, .firstStartBC = 0});
  frames.resize(11U, {.lengthBC = 81, .firstStartBC = 4});
  const auto params = gridParams(param, frames);
  // 10 us is 4.95 frames of 81 BC: five a block.
  BOOST_TEST(params.partsPerPS == kPartsPerPS);
  BOOST_TEST(params.blockWidth == kPartsPerBC * 5 * 81);
  BOOST_TEST(params.blockOrigin == (4 * kPartsPerBC) - kPartsPerPS);
  BOOST_TEST(std::ranges::equal(frames, std::span{params.frames}.first(frames.size())));
}

/// One length on every layer: blocks of two readout frames, their edges a
/// picosecond before the frames' starts. Every frame's box in a 32-orbit
/// timeframe is binned in its block, and the zero-margin candidate range of a
/// block's second frame lies in that block.
BOOST_AUTO_TEST_CASE(OneLengthKeepsItsBlocks)
{
  const RecoSettings param{};
  const std::vector<FrameTiming> frames(7U, {.lengthBC = 198, .firstStartBC = 64});
  const auto params = gridParams(param, frames);
  BOOST_TEST(params.blockWidth == kPartsPerBC * 2 * 198);
  BOOST_TEST(params.blockOrigin == (64 * kPartsPerBC) - kPartsPerPS);

  constexpr std::int64_t kFrameBC = 198;
  constexpr std::int64_t kFirstBC = 64;
  constexpr std::int64_t kDurationBC = std::int64_t{32} * o2::constants::lhc::LHCMaxBunches;
  const IterationSettings settings{};
  const auto windows = gridLayerWindows<detectors::ToyBarrel3>(std::span{&settings, 1U}, VertexConstraint{});
  const GridGeometry<detectors::ToyBarrel3> geometry{psOfBC(kDurationBC), params, windows};
  const auto frameCount = (kDurationBC - kFirstBC) / kFrameBC;
  for (std::int64_t frame = 0; frame < frameCount; ++frame) {
    const auto start = kFirstBC + (frame * kFrameBC);
    const auto block = static_cast<std::int32_t>(frame / 2);
    const auto box = frameBox(start, start + kFrameBC, 0);
    BOOST_REQUIRE_EQUAL(geometry.blockOfPS(box.lower()), block);
    if (frame % 2 == 1) {
      const LayerPairWindow pair{.phiHalfWidth = 0.1F, .secondHalfWidth = 1.F};
      const auto window = candidateWindow(pair, box, 0.F, 0.F, 0.F, box.width);
      BOOST_REQUIRE_EQUAL(geometry.blockOfPS(window.earliestPS), block);
      BOOST_REQUIRE_EQUAL(geometry.blockOfPS(window.latestPS), block);
    }
  }
  BOOST_TEST(geometry.blockStart(1) == params.blockOrigin + params.blockWidth);
}

} // namespace o2::gtrk::workflow
