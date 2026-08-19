// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Focused coverage for assembleRefitLegSlots in the native refit driver.
// Exercises both traversal
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
#include <memory>
#include <vector>

#include "ITSMFTTracking/RefitDriver.h"

using namespace o2::itsmft::tracking;

namespace
{

constexpr int NLayers = ITSNLayers;

// Distinguishable, non-hole measurement fingerprinted by legacy layer and
// cluster index, so a slot's origin can be verified unambiguously.
SurfaceMeasurement makeMeasurement(int legacyLayer, int clusterIndex)
{
  SurfaceMeasurement m{};
  m.frame = {static_cast<float>(legacyLayer), static_cast<float>(clusterIndex),
             static_cast<float>(100 * legacyLayer + clusterIndex), 0.f};
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
  std::array<std::vector<GlobalMeasurement>, NLayers> globalStorage;
  std::vector<gsl::span<const GlobalMeasurement>> layerGlobals = std::vector<gsl::span<const GlobalMeasurement>>(NLayers);
  std::array<LayerId, NLayers> orderedSurfaces{};
  TimeFrame frame;
  TrackSeed seed{};

  Fixture()
  {
    for (int layer = 0; layer < NLayers; ++layer) {
      storage[layer] = {makeMeasurement(layer, 0), makeMeasurement(layer, 1)};
      globalStorage[layer].resize(2);
      globalStorage[layer][0].clusterId = 0;
      globalStorage[layer][1].clusterId = 1;
      layerGlobals[layer] = gsl::span<const GlobalMeasurement>(globalStorage[layer]);
      orderedSurfaces[layer] = LayerId{static_cast<uint16_t>(layer)};
    }
    std::array<SurfaceDescriptor, NLayers> surfaces{};
    for (int layer = 0; layer < NLayers; ++layer) {
      surfaces[layer].id = orderedSurfaces[layer];
      surfaces[layer].kind = SurfaceKind::Cylinder;
    }
    BOOST_REQUIRE(frame.configure(SurfaceLayout{surfaces, makeSurfaceLayoutChain(orderedSurfaces)}, 0, 0,
                                  std::make_shared<BoundedMemoryResource>()));
    for (int layer = 0; layer < NLayers; ++layer) {
      for (std::size_t cluster = 0; cluster < globalStorage[layer].size(); ++cluster) {
        frame.addMeasurement(orderedSurfaces[layer], globalStorage[layer][cluster], storage[layer][cluster]);
      }
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
  std::array<detail::RefitMeasurementSlot, NLayers> out{};
  bool valid = false;
  const auto slots = detail::assembleRefitLegSlots(fx.seed, fx.frame, fx.layerGlobals, fx.orderedSurfaces, 0, NLayers, 1, out, valid);

  BOOST_REQUIRE(valid);
  BOOST_REQUIRE_EQUAL(slots.size(), NLayers);
  BOOST_CHECK(bitEqual(slots[0].measurement, fx.storage[0][0]));
  BOOST_CHECK(!slots[1].present);
  BOOST_CHECK(bitEqual(slots[2].measurement, fx.storage[2][1]));
  BOOST_CHECK(!slots[3].present);
  BOOST_CHECK(bitEqual(slots[4].measurement, fx.storage[4][0]));
  BOOST_CHECK(!slots[5].present);
  BOOST_CHECK(bitEqual(slots[6].measurement, fx.storage[6][1]));
}

BOOST_AUTO_TEST_CASE(ReverseLegProducesDecreasingLayerOrder)
{
  Fixture fx;
  std::array<detail::RefitMeasurementSlot, NLayers> out{};
  bool valid = false;
  const auto slots = detail::assembleRefitLegSlots(fx.seed, fx.frame, fx.layerGlobals, fx.orderedSurfaces, NLayers - 1, -1, -1, out, valid);

  BOOST_REQUIRE(valid);
  BOOST_REQUIRE_EQUAL(slots.size(), NLayers);
  // Exact mirror of the forward-leg order: position 0 is legacy layer 6,
  // position 6 is legacy layer 0. A buggy implementation that writes
  // out[legacyLayer] instead of out[position++] would return the SAME
  // ascending-layer sequence as the forward test above; these assertions
  // fail against that bug.
  BOOST_CHECK(bitEqual(slots[0].measurement, fx.storage[6][1]));
  BOOST_CHECK(!slots[1].present);
  BOOST_CHECK(bitEqual(slots[2].measurement, fx.storage[4][0]));
  BOOST_CHECK(!slots[3].present);
  BOOST_CHECK(bitEqual(slots[4].measurement, fx.storage[2][1]));
  BOOST_CHECK(!slots[5].present);
  BOOST_CHECK(bitEqual(slots[6].measurement, fx.storage[0][0]));
}

BOOST_AUTO_TEST_CASE(ForwardAndReverseOrdersAreNotEqual)
{
  // Directly asserts the two directions disagree in general (rules out a
  // regression where both directions accidentally return the same order,
  // e.g. because the bug produces an order that is its own reverse for a
  // symmetric fixture -- this fixture is deliberately asymmetric: the
  // cluster slot chosen per layer, 0/1/0/1, is not palindromic).
  Fixture fx;
  std::array<detail::RefitMeasurementSlot, NLayers> outFwd{};
  std::array<detail::RefitMeasurementSlot, NLayers> outRev{};
  bool validFwd = false;
  bool validRev = false;
  const auto fwd = detail::assembleRefitLegSlots(fx.seed, fx.frame, fx.layerGlobals, fx.orderedSurfaces, 0, NLayers, 1, outFwd, validFwd);
  const auto rev = detail::assembleRefitLegSlots(fx.seed, fx.frame, fx.layerGlobals, fx.orderedSurfaces, NLayers - 1, -1, -1, outRev, validRev);

  BOOST_REQUIRE(validFwd);
  BOOST_REQUIRE(validRev);
  BOOST_REQUIRE_EQUAL(fwd.size(), rev.size());
  bool anyDifferent = false;
  for (size_t i = 0; i < fwd.size(); ++i) {
    if (fwd[i].present != rev[i].present || !bitEqual(fwd[i].measurement, rev[i].measurement)) {
      anyDifferent = true;
      break;
    }
  }
  BOOST_CHECK(anyDifferent);
}

BOOST_AUTO_TEST_CASE(AllHoleLegProducesAllInvalidSlots)
{
  Fixture fx;
  TrackSeed allHoleSeed{}; // every layer left at UnusedIndex
  std::array<detail::RefitMeasurementSlot, NLayers> out{};
  bool valid = false;
  const auto slots = detail::assembleRefitLegSlots(allHoleSeed, fx.frame, fx.layerGlobals, fx.orderedSurfaces, 0, NLayers, 1, out, valid);

  BOOST_REQUIRE(valid);
  BOOST_REQUIRE_EQUAL(slots.size(), NLayers);
  for (const auto& slot : slots) {
    BOOST_CHECK(!slot.present);
  }
}

BOOST_AUTO_TEST_CASE(PartialRangeReturnsOnlyThePopulatedPrefix)
{
  // Not exercised by the frozen refitTrack driver (its three legs always
  // span the complete [0,NLayers)/[NLayers-1,-1) range), but the [0,
  // position) contract must hold generally: a caller must never assume the
  // returned span has NLayers elements.
  Fixture fx;
  std::array<detail::RefitMeasurementSlot, NLayers> out{};
  bool valid = false;
  const auto slots = detail::assembleRefitLegSlots(fx.seed, fx.frame, fx.layerGlobals, fx.orderedSurfaces, 2, 5, 1, out, valid);

  BOOST_REQUIRE(valid);
  BOOST_REQUIRE_EQUAL(slots.size(), 3);
  BOOST_CHECK(bitEqual(slots[0].measurement, fx.storage[2][1]));
  BOOST_CHECK(!slots[1].present);
  BOOST_CHECK(bitEqual(slots[2].measurement, fx.storage[4][0]));
}
