// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 3 Slice A: focused coverage for assembleRefitLegSlots
// (RefitLegAssembly.h) -- the ordered-slot assembly prerequisite for the
// native refit driver (Slice B, unwired). Exercises both traversal
// directions explicitly: a forward (increasing-index) leg and a reverse
// (decreasing-index) leg must produce different, correctly-ordered slot
// sequences from the same seed/measurement fixture -- this is the exact
// bug class the design review flagged (writing slots by legacy layer index
// instead of by traversal position would silently return the same
// ascending-layer order for both directions).

#define BOOST_TEST_MODULE ITSMFTRefitLegAssembly
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstring>
#include <vector>

#include "ITSMFTTracking/RefitLegAssembly.h"

using namespace o2::itsmft::tracking;

namespace
{

constexpr int NLayers = ITSNLayers;

// Distinguishable, non-hole measurement fingerprinted by legacy layer and
// cluster index, so a slot's origin can be verified unambiguously.
SurfaceMeasurement makeMeasurement(int legacyLayer, int clusterIndex)
{
  SurfaceMeasurement m{};
  m.surface = SurfaceId{static_cast<uint16_t>(100 * legacyLayer + clusterIndex)};
  m.cluster = ClusterRef{ClusterSourceId{0}, static_cast<uint32_t>(100 * legacyLayer + clusterIndex)};
  return m;
}

bool bitEqual(const SurfaceMeasurement& lhs, const SurfaceMeasurement& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(SurfaceMeasurement)) == 0;
}

// Fixture: clusters attached at even layers (0,2,4,6), holes at odd layers
// (1,3,5). Each even layer's measurement span holds two distinguishable
// entries; the seed picks a different slot per layer (0 or 1) so a
// layer-index/cluster-index transposition bug would also be caught.
struct Fixture {
  std::array<std::vector<SurfaceMeasurement>, NLayers> storage;
  LayerMeasurementSpans<NLayers> layerMeasurements{};
  TrackSeedN<NLayers> seed{};

  Fixture()
  {
    for (int layer = 0; layer < NLayers; ++layer) {
      storage[layer] = {makeMeasurement(layer, 0), makeMeasurement(layer, 1)};
      layerMeasurements[layer] = gsl::span<const SurfaceMeasurement>(storage[layer]);
    }
    seed.getClusters()[0] = 0;
    seed.getClusters()[2] = 1;
    seed.getClusters()[4] = 0;
    seed.getClusters()[6] = 1;
    // layers 1, 3, 5 keep their default UnusedIndex -> holes.
  }
};

} // namespace

BOOST_AUTO_TEST_CASE(ForwardLegProducesIncreasingLayerOrder)
{
  Fixture fx;
  std::array<SurfaceMeasurement, NLayers> out{};
  const auto slots = assembleRefitLegSlots<NLayers>(fx.seed, fx.layerMeasurements, 0, NLayers, 1, out);

  BOOST_REQUIRE_EQUAL(slots.size(), NLayers);
  BOOST_CHECK(bitEqual(slots[0], fx.storage[0][0])); // layer 0, cluster 0
  BOOST_CHECK(!slots[1].cluster.isValid());          // layer 1, hole
  BOOST_CHECK(bitEqual(slots[2], fx.storage[2][1])); // layer 2, cluster 1
  BOOST_CHECK(!slots[3].cluster.isValid());          // layer 3, hole
  BOOST_CHECK(bitEqual(slots[4], fx.storage[4][0])); // layer 4, cluster 0
  BOOST_CHECK(!slots[5].cluster.isValid());          // layer 5, hole
  BOOST_CHECK(bitEqual(slots[6], fx.storage[6][1])); // layer 6, cluster 1
}

BOOST_AUTO_TEST_CASE(ReverseLegProducesDecreasingLayerOrder)
{
  Fixture fx;
  std::array<SurfaceMeasurement, NLayers> out{};
  const auto slots = assembleRefitLegSlots<NLayers>(fx.seed, fx.layerMeasurements, NLayers - 1, -1, -1, out);

  BOOST_REQUIRE_EQUAL(slots.size(), NLayers);
  // Exact mirror of the forward-leg order: position 0 is legacy layer 6,
  // position 6 is legacy layer 0. A buggy implementation that writes
  // out[legacyLayer] instead of out[position++] would return the SAME
  // ascending-layer sequence as the forward test above; these assertions
  // fail against that bug.
  BOOST_CHECK(bitEqual(slots[0], fx.storage[6][1])); // layer 6, cluster 1
  BOOST_CHECK(!slots[1].cluster.isValid());          // layer 5, hole
  BOOST_CHECK(bitEqual(slots[2], fx.storage[4][0])); // layer 4, cluster 0
  BOOST_CHECK(!slots[3].cluster.isValid());          // layer 3, hole
  BOOST_CHECK(bitEqual(slots[4], fx.storage[2][1])); // layer 2, cluster 1
  BOOST_CHECK(!slots[5].cluster.isValid());          // layer 1, hole
  BOOST_CHECK(bitEqual(slots[6], fx.storage[0][0])); // layer 0, cluster 0
}

BOOST_AUTO_TEST_CASE(ForwardAndReverseOrdersAreNotEqual)
{
  // Directly asserts the two directions disagree in general (rules out a
  // regression where both directions accidentally return the same order,
  // e.g. because the bug produces an order that is its own reverse for a
  // symmetric fixture -- this fixture is deliberately asymmetric: the
  // cluster slot chosen per layer, 0/1/0/1, is not palindromic).
  Fixture fx;
  std::array<SurfaceMeasurement, NLayers> outFwd{};
  std::array<SurfaceMeasurement, NLayers> outRev{};
  const auto fwd = assembleRefitLegSlots<NLayers>(fx.seed, fx.layerMeasurements, 0, NLayers, 1, outFwd);
  const auto rev = assembleRefitLegSlots<NLayers>(fx.seed, fx.layerMeasurements, NLayers - 1, -1, -1, outRev);

  BOOST_REQUIRE_EQUAL(fwd.size(), rev.size());
  bool anyDifferent = false;
  for (size_t i = 0; i < fwd.size(); ++i) {
    if (!bitEqual(fwd[i], rev[i])) {
      anyDifferent = true;
      break;
    }
  }
  BOOST_CHECK(anyDifferent);
}

BOOST_AUTO_TEST_CASE(AllHoleLegProducesAllInvalidSlots)
{
  Fixture fx;
  TrackSeedN<NLayers> allHoleSeed{}; // every layer left at UnusedIndex
  std::array<SurfaceMeasurement, NLayers> out{};
  const auto slots = assembleRefitLegSlots<NLayers>(allHoleSeed, fx.layerMeasurements, 0, NLayers, 1, out);

  BOOST_REQUIRE_EQUAL(slots.size(), NLayers);
  for (const auto& slot : slots) {
    BOOST_CHECK(!slot.cluster.isValid());
  }
}

BOOST_AUTO_TEST_CASE(PartialRangeReturnsOnlyThePopulatedPrefix)
{
  // Not exercised by the frozen refitTrack driver (its three legs always
  // span the complete [0,NLayers)/[NLayers-1,-1) range), but the [0,
  // position) contract must hold generally: a caller must never assume the
  // returned span has NLayers elements.
  Fixture fx;
  std::array<SurfaceMeasurement, NLayers> out{};
  const auto slots = assembleRefitLegSlots<NLayers>(fx.seed, fx.layerMeasurements, 2, 5, 1, out);

  BOOST_REQUIRE_EQUAL(slots.size(), 3);
  BOOST_CHECK(bitEqual(slots[0], fx.storage[2][1])); // layer 2
  BOOST_CHECK(!slots[1].cluster.isValid());          // layer 3, hole
  BOOST_CHECK(bitEqual(slots[2], fx.storage[4][0])); // layer 4
}
