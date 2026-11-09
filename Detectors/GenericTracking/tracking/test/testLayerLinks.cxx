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

/// \file testLayerLinks.cxx
/// \brief What follows from a detector's link table, and the grid it sizes.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK LayerLinks
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Base/LayerLinks.h"
#include "ToyBarrel3.h"
#include "ToyBarrel4Disc8.h"
#include "GenericTracking/Tracking/GridSizing.h"
#include "GenericTracking/Tracking/LinkHoles.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"

#include "GridTestSupport.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace o2::gtrk
{

namespace
{

using detectors::ToyBarrel3;
using detectors::ToyBarrel4Disc8;

/// Whether `directPairRange` leaves anything.
bool directPairPossible(std::uint64_t holes, std::uint64_t barrels, const LayerSpan* spans, LayerID target, Scalar ra, Scalar za, Scalar low, Scalar high)
{
  const auto range = directPairRange(holes, barrels, spans, target, ra, za, low, high);
  return range.low <= range.high;
}

/// ToyBarrel3's layers with links whose order is not the layer numbering:
/// 0 -> 2 -> 1, and 0 -> 1 directly past 2.
struct Reordered {
  static constexpr const char* name = "Reordered";
  static constexpr auto barrelSpecs = ToyBarrel3::barrelSpecs;
  static constexpr auto discSpecs = ToyBarrel3::discSpecs;
  static constexpr auto barrelLayers = ToyBarrel3::barrelLayers;
  static constexpr auto discLayers = ToyBarrel3::discLayers;
  static constexpr std::array<LayerLink, 3> layerLinks{{
    {.from = 2, .to = 1, .sourceMin = -1., .sourceMax = 1., .skipsLayer = false},
    {.from = 0, .to = 2, .sourceMin = -1., .sourceMax = 1., .skipsLayer = false},
    {.from = 0, .to = 1, .sourceMin = -0.5, .sourceMax = 0.5, .skipsLayer = true},
  }};
  static constexpr int nSensors = ToyBarrel3::nSensors;
};

/// ToyBarrel3's inner layers with a two-module outer one, |z| 2 cm at r 8:
/// 0 -> 2 passes beyond layer 1's |z| 1 cm for steep tracks and crosses it
/// for the others; from |z| below 0.5 on layer 0 every pair crosses it.
struct PastTheEnd {
  static constexpr const char* name = "PastTheEnd";
  static constexpr std::array<BarrelLayerSpec, 3> barrelSpecs{{
    ToyBarrel3::barrelSpecs[0],
    ToyBarrel3::barrelSpecs[1],
    {.r = 8., .nStaves = 4, .nModulesPerStave = 2, .phi0 = 0., .tilt = 0., .moduleGap = 0., .sensor = ToyBarrel3::sensor},
  }};
  static constexpr auto discSpecs = ToyBarrel3::discSpecs;
  static constexpr auto barrelLayers = placedLayers(barrelSpecs);
  static constexpr auto discLayers = placedLayers(discSpecs);
  static constexpr std::array<LayerLink, 3> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false},
    {.from = 0, .to = 2, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false},
    {.from = 1, .to = 2, .sourceMin = -1.5, .sourceMax = 1.5, .skipsLayer = false},
  }};
  static constexpr int nSensors = countSensors(barrelLayers, discLayers);
};

[[maybe_unused]] const bool kReorderedPublished = publishToy<Reordered>();
[[maybe_unused]] const bool kPastTheEndPublished = publishToy<PastTheEnd>();

/// Every layer's span and the barrel bits, as an iteration holds them.
template <DetectorLayout Det>
auto spansOf()
{
  std::array<LayerSpan, static_cast<std::size_t>(Detector<Det>::nLayers)> spans{};
  for (std::size_t layer = 0; layer < spans.size(); ++layer) {
    spans[layer] = layerSpan<Det>(static_cast<LayerID>(layer));
  }
  return spans;
}

/// For every link of \a Det with layers a pair can cross: points over the
/// source's and the target's extents, the per-source check against every
/// pair. A direct pair from a source the check refuses fails; returns the
/// sources refused and the sources checked.
template <DetectorLayout Det>
std::pair<std::size_t, std::size_t> checkEveryDirectPairIsPossible()
{
  constexpr int kSteps = 25;
  const auto spans = spansOf<Det>();
  const auto barrels = barrelLayerBits<Det>();
  const auto along = [](Scalar low, Scalar high, int i, int steps) { return low + ((high - low) * static_cast<Scalar>(i) / static_cast<Scalar>(steps - 1)); };
  // A point of a layer's extent, as (r, z).
  const auto point = [&](LayerID layer, int i, int j, int positionSteps) {
    const auto& span = spans[layer];
    const Scalar position = along(span.positionMin, span.positionMax, i, positionSteps);
    const Scalar second = along(span.low, span.high, j, kSteps);
    return ((barrels >> layer) & 1U) != 0U ? std::pair{position, second} : std::pair{second, position};
  };
  std::size_t refused = 0;
  std::size_t checked = 0;
  for (const auto& link : Det::layerLinks) {
    const auto holes = linkHoleLayers<Det>(link.from, link.to);
    if (holes == 0U) {
      continue;
    }
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < kSteps; ++j) {
        const auto [ra, za] = point(link.from, i, j, 3);
        const bool possible = directPairPossible(holes, barrels, spans.data(), link.to, ra, za, spans[link.to].low, spans[link.to].high);
        ++checked;
        refused += possible ? 0U : 1U;
        for (int k = 0; k < 3 && !possible; ++k) {
          for (int l = 0; l < kSteps; ++l) {
            const auto [rb, zb] = point(link.to, k, l, 3);
            BOOST_TEST_CONTEXT(Det::name << " link " << link.from << " -> " << link.to << " from (" << ra << ", " << za << ") to (" << rb << ", " << zb << ")") { BOOST_TEST(segmentLeavesHole(holes, barrels, spans.data(), ra, za, rb, zb)); }
          }
        }
      }
    }
  }
  return {refused, checked};
}

/// Where each layer sits in `layerOrder`.
template <DetectorLayout Det>
constexpr auto positions()
{
  std::array<std::size_t, static_cast<std::size_t>(Detector<Det>::nLayers)> result{};
  for (std::size_t i = 0; i < result.size(); ++i) {
    result[layerOrder<Det>[i]] = i;
  }
  return result;
}

template <DetectorLayout Det>
constexpr bool everyLinkPointsForward()
{
  constexpr auto position = positions<Det>();
  return std::ranges::all_of(Det::layerLinks, [&](const LayerLink& link) { return position[link.from] < position[link.to]; });
}

/// Every link appears once among the out-links of its source and once among
/// the in-links of its target.
template <DetectorLayout Det>
constexpr bool indexedOnceEachWay()
{
  for (std::size_t i = 0; i < nLinks<Det>; ++i) {
    const auto& link = Det::layerLinks[i];
    std::size_t out = 0;
    for (auto k = outLinks<Det>.begin(link.from); k < outLinks<Det>.end(link.from); ++k) {
      out += outLinks<Det>.links[k] == i ? 1U : 0U;
    }
    std::size_t in = 0;
    for (auto k = inLinks<Det>.begin(link.to); k < inLinks<Det>.end(link.to); ++k) {
      in += inLinks<Det>.links[k] == i ? 1U : 0U;
    }
    if (out != 1 || in != 1) {
      return false;
    }
  }
  return outLinks<Det>.offsets.back() == nLinks<Det> && inLinks<Det>.offsets.back() == nLinks<Det>;
}

} // namespace

BOOST_AUTO_TEST_CASE(TheOrderFollowsTheLinksNotTheNumbering)
{
  static_assert(layerOrder<Reordered> == std::array<LayerID, 3>{0, 2, 1});
  static_assert(everyLinkPointsForward<Reordered>());
  static_assert(indexedOnceEachWay<Reordered>());
  // 0 -> 2 -> 1 is the longest chain; 0 -> 1 skips 2 and ends one.
  static_assert(maxChainLayers<Reordered> == 3);

  // Out-links keep table order within a layer.
  static_assert(outLinks<Reordered>.count(0) == 2);
  static_assert(Reordered::layerLinks[outLinks<Reordered>.links[outLinks<Reordered>.begin(0)]].to == 2);
  static_assert(inLinks<Reordered>.count(1) == 2);
  static_assert(outLinks<Reordered>.count(1) == 0);
  BOOST_TEST(true);
}

BOOST_AUTO_TEST_CASE(TheBarrelAndDiscTableIsAcyclicAndIndexed)
{
  // Compiling layerOrder is the acyclicity check; these pin what it gives.
  static_assert(everyLinkPointsForward<ToyBarrel4Disc8>());
  static_assert(indexedOnceEachWay<ToyBarrel4Disc8>());

  BOOST_TEST(maxChainLayers<ToyBarrel4Disc8> >= 4U); // at least the whole barrel
}

/// A crossing link must not size cells: its window is the target's whole
/// second axis, and taken as a bound it would leave one cell there. With the
/// tracklet queries alone (no extension), every layer a same-kind link enters
/// keeps both axes binned, and a layer only crossing links enter is not sized
/// at all: one cell per axis.
///
/// At the 1 GeV of the first iteration, with every other default but a beam
/// spot narrow against the toy's layers, whose windows it would otherwise
/// spread over a whole barrel layer.
BOOST_AUTO_TEST_CASE(CrossingLinksDoNotCollapseTheGrid)
{
  SearchWindowParams windows;
  windows.minPtGeV = 1.;
  const IterationSettings tracklets{.windows = windows, .follow = false};
  const auto layers = gridLayerWindows<ToyBarrel4Disc8>(std::span<const IterationSettings>{&tracklets, 1U}, VertexConstraint{.sigmaZ = 1.});
  const GridGeometry<ToyBarrel4Disc8> grid{100'000, test::withToyFrames(GridParams{.blockWidth = 100'000}), layers};
  for (int layer = 0; layer < Detector<ToyBarrel4Disc8>::nLayers; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    bool sameKindIn = false;
    for (const auto& link : ToyBarrel4Disc8::layerLinks) {
      sameKindIn = sameKindIn || (link.to == id && Detector<ToyBarrel4Disc8>::isBarrelLayer(link.from) == Detector<ToyBarrel4Disc8>::isBarrelLayer(id));
    }
    const auto& axes = grid.axes(id);
    BOOST_TEST_CONTEXT("layer " << layer)
    {
      if (sameKindIn) {
        BOOST_TEST(axes.nPhi > 1);
        BOOST_TEST(axes.nSecond > 1);
      } else {
        BOOST_TEST(axes.cellCount() == 1);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(APairLeavesAHoleWhereItsSegmentCrossesALayer)
{
  BOOST_TEST(linkHoleLayers<PastTheEnd>(0, 1) == 0U);
  BOOST_TEST(linkHoleLayers<PastTheEnd>(0, 2) == 0b010U);
  BOOST_TEST(linkHoleLayers<PastTheEnd>(1, 2) == 0U);
  const auto spans = spansOf<PastTheEnd>();
  const auto barrels = barrelLayerBits<PastTheEnd>();
  BOOST_TEST(barrels == 0b111U);
  // Layer 1 at r 4 over |z| 1: from z 0.9 at r 2 to z 1.9 at r 8 the segment
  // is at z 1.23 there; at z 0 throughout it crosses.
  BOOST_TEST(!segmentLeavesHole(0b010U, barrels, spans.data(), 2.F, 0.9F, 8.F, 1.9F));
  BOOST_TEST(segmentLeavesHole(0b010U, barrels, spans.data(), 2.F, 0.F, 8.F, 0.F));
  // From z 0.9 the segments reaching z beyond 1.2 at r 8 pass; the check
  // keeps the source. From z 0 layer 1's shadow is |z| 3, past layer 2's 2.
  BOOST_TEST(directPairPossible(0b010U, barrels, spans.data(), 2, 2.F, 0.9F, spans[2].low, spans[2].high));
  BOOST_TEST(!directPairPossible(0b010U, barrels, spans.data(), 2, 2.F, 0.F, spans[2].low, spans[2].high));
}

/// The per-source check refuses only sources with no direct pair, on a
/// barrel toy and on barrels with discs; and it does refuse some.
BOOST_AUTO_TEST_CASE(ThePerSourceCheckNeverRefusesADirectPair)
{
  const auto [toyRefused, toyChecked] = checkEveryDirectPairIsPossible<PastTheEnd>();
  BOOST_TEST(toyRefused > 0U);
  BOOST_TEST(toyRefused < toyChecked);
  const auto [mixedRefused, mixedChecked] = checkEveryDirectPairIsPossible<ToyBarrel4Disc8>();
  BOOST_TEST(mixedRefused > 0U);
  BOOST_TEST(mixedRefused < mixedChecked);
}

/// Without skips a link whose pairs can cross a layer is queried only at that
/// layer's ends: it does not size a layer a link without holes enters. With
/// skips every source queries it, and it does.
BOOST_AUTO_TEST_CASE(EdgeOnlyLinksDoNotSizeTheGrid)
{
  SearchWindowParams windows;
  windows.minPtGeV = 0.2;
  windows.bzKiloGauss = 20.;
  const IterationSettings without{.windows = windows, .minClusters = 2, .skipLinks = false, .follow = false};
  const IterationSettings with{.windows = windows, .minClusters = 2, .skipLinks = true, .follow = false};
  const auto noSkips = gridLayerWindows<PastTheEnd>(std::span<const IterationSettings>{&without, 1U}, VertexConstraint{});
  const auto skips = gridLayerWindows<PastTheEnd>(std::span<const IterationSettings>{&with, 1U}, VertexConstraint{});
  const auto adjacent = layerPairWindow<PastTheEnd>(1, 2, windows, VertexConstraint{});
  const auto skipping = layerPairWindow<PastTheEnd>(0, 2, windows, VertexConstraint{});
  BOOST_REQUIRE(skipping.phiHalfWidth > adjacent.phiHalfWidth);
  BOOST_TEST(noSkips[2].phiHalfWidth == static_cast<double>(adjacent.phiHalfWidth));
  BOOST_TEST(skips[2].phiHalfWidth == static_cast<double>(skipping.phiHalfWidth));
}

BOOST_AUTO_TEST_CASE(SegmentsCrossBarrelsAndDiscsWithinTheirExtent)
{
  BOOST_TEST(segmentCrosses(true, {.position = 4.F, .low = -1.F, .high = 1.F}, 2.F, 0.F, 8.F, 0.F));
  BOOST_TEST(!segmentCrosses(true, {.position = 4.F, .low = -1.F, .high = 1.F}, 2.F, 0.9F, 8.F, 3.9F));
  BOOST_TEST(segmentCrosses(false, {.position = 10.F, .low = 3.F, .high = 20.F}, 2.F, 5.F, 12.F, 15.F));
  BOOST_TEST(!segmentCrosses(false, {.position = 10.F, .low = 3.F, .high = 20.F}, 1.F, 5.F, 2.5F, 15.F));
}

} // namespace o2::gtrk
