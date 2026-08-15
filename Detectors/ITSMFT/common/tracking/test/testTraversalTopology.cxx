// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.

#define BOOST_TEST_MODULE ITSMFT TraversalTopology
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <type_traits>
#include <vector>

#include "ITSMFTTracking/TraversalTopology.h"

namespace
{
using namespace o2::itsmft::tracking;

std::vector<SurfaceDescriptor> catalog(uint16_t count)
{
  std::vector<SurfaceDescriptor> result;
  result.reserve(count);
  for (uint16_t id = 0; id < count; ++id) {
    result.push_back(SurfaceDescriptor{SurfaceId{id}, id, 0, SurfaceKind::Cylinder});
  }
  return result;
}

SurfaceLayout makeLayout(std::vector<SurfaceId> ordered,
                         std::vector<uint16_t> componentOffsets = {0},
                         int maxHoles = 0,
                         SurfaceMask holeSurfaces = {},
                         SurfaceMask seedingSurfaces = {})
{
  const auto surfaces = catalog(static_cast<uint16_t>(std::max<uint16_t>(4, ordered.size())));
  SurfaceLayoutDefinition definition;
  definition.orderedSurfaces = std::move(ordered);
  definition.componentOffsets = std::move(componentOffsets);
  definition.maxHoles = maxHoles;
  definition.holeSurfaces = holeSurfaces;
  definition.seedingSurfaces = seedingSurfaces;
  return SurfaceLayout{surfaces, std::move(definition)};
}

std::vector<SurfaceId> chain(std::initializer_list<uint16_t> ids)
{
  std::vector<SurfaceId> result;
  result.reserve(ids.size());
  for (const auto id : ids) {
    result.emplace_back(id);
  }
  return result;
}

SurfaceMask mask(std::initializer_list<uint16_t> ids)
{
  SurfaceMask result;
  for (const auto id : ids) {
    result.set(SurfaceId{id});
  }
  return result;
}

const Edge* findEdge(const TraversalTopology& topology, SurfaceId from, SurfaceId to)
{
  const auto edge = std::find_if(topology.edges.begin(), topology.edges.end(), [&](const auto& candidate) {
    return candidate.from == from && candidate.to == to;
  });
  return edge == topology.edges.end() ? nullptr : &*edge;
}
} // namespace

BOOST_AUTO_TEST_CASE(CellPathContainsOnlyTwoEdgeIds)
{
  static_assert(std::is_standard_layout_v<CellPath>);
  static_assert(std::is_trivially_copyable_v<CellPath>);
  static_assert(std::is_same_v<decltype(CellPath::first), EdgeId>);
  static_assert(std::is_same_v<decltype(CellPath::second), EdgeId>);
  static_assert(sizeof(CellPath) == sizeof(EdgeId) + sizeof(EdgeId));
  BOOST_CHECK_EQUAL(sizeof(CellPath), 4u);
}

BOOST_AUTO_TEST_CASE(ComponentBoundariesRejectCrossComponentEdges)
{
  const auto layout = makeLayout(chain({0, 1, 2, 3}), {0, 2});
  const auto result = deriveTraversalTopology(layout);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result.topology->edges.size(), 2u);
  BOOST_CHECK(findEdge(*result.topology, SurfaceId{1}, SurfaceId{2}) == nullptr);
}

BOOST_AUTO_TEST_CASE(AllActiveChainDerivesEdgesAndCellPaths)
{
  const auto layout = makeLayout(chain({0, 1, 2, 3}));
  const auto result = deriveTraversalTopology(layout);
  BOOST_REQUIRE(result.ok());
  const auto& topology = *result.topology;
  BOOST_CHECK_EQUAL(topology.orderedSurfaces.size(), 4u);
  BOOST_CHECK_EQUAL(topology.activeSurfaceList.size(), 4u);
  BOOST_CHECK_EQUAL(topology.edges.size(), 3u);
  BOOST_CHECK_EQUAL(topology.paths.size(), 2u);
  BOOST_CHECK(topology.edges[0].from == SurfaceId{0});
  BOOST_CHECK(topology.edges[0].to == SurfaceId{1});
  BOOST_CHECK(topology.edges[1].from == SurfaceId{1});
  BOOST_CHECK(topology.edges[1].to == SurfaceId{2});
  BOOST_CHECK(topology.edges[2].from == SurfaceId{2});
  BOOST_CHECK(topology.edges[2].to == SurfaceId{3});
  BOOST_CHECK(topology.paths[0].first == EdgeId{0});
  BOOST_CHECK(topology.paths[0].second == EdgeId{1});
  BOOST_CHECK(topology.paths[1].first == EdgeId{1});
  BOOST_CHECK(topology.paths[1].second == EdgeId{2});
}

BOOST_AUTO_TEST_CASE(DisabledMiddleSurfaceRetainsAdmittedBridge)
{
  const auto layout = makeLayout(chain({0, 1, 2, 3}), {0}, 1, mask({1}));
  const auto result = deriveTraversalTopology(layout, mask({1}));
  BOOST_REQUIRE(result.ok());
  const auto& topology = *result.topology;
  BOOST_CHECK_EQUAL(topology.activeSurfaceList.size(), 3u);
  BOOST_CHECK_EQUAL(topology.edges.size(), 2u);
  BOOST_CHECK_EQUAL(topology.paths.size(), 1u);
  const auto* bridge = findEdge(topology, SurfaceId{0}, SurfaceId{2});
  BOOST_REQUIRE(bridge != nullptr);
  BOOST_CHECK(bridge->skippedSurfaces == mask({1}));
  BOOST_CHECK(topology.paths[0].first == EdgeId{0});
  BOOST_CHECK(topology.paths[0].second == EdgeId{1});
}

BOOST_AUTO_TEST_CASE(DisabledEndpointOmitsItsEdges)
{
  const auto layout = makeLayout(chain({0, 1, 2, 3}), {0}, 1, mask({1}));
  const auto result = deriveTraversalTopology(layout, mask({0}));
  BOOST_REQUIRE(result.ok());
  for (const auto& edge : result.topology->edges) {
    BOOST_CHECK(edge.from != SurfaceId{0});
    BOOST_CHECK(edge.to != SurfaceId{0});
  }
  BOOST_CHECK(findEdge(*result.topology, SurfaceId{1}, SurfaceId{2}) != nullptr);
}

BOOST_AUTO_TEST_CASE(InvalidDerivationIsTransactional)
{
  const auto layout = makeLayout(chain({0, 1, 2, 3}));
  const auto result = deriveTraversalTopology(layout, mask({7}));
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(!result.topology.has_value());
  BOOST_CHECK(result.error == TraversalTopologyError::DisabledSurfaceOutsideLayout);

  SurfaceLayout invalid;
  const auto invalidResult = deriveTraversalTopology(invalid);
  BOOST_CHECK(!invalidResult.ok());
  BOOST_CHECK(!invalidResult.topology.has_value());
  BOOST_CHECK(invalidResult.error == TraversalTopologyError::InvalidLayout);
}
