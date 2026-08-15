// Copyright 2019-2026 CERN and copyright holders of ALICE O2.

#define BOOST_TEST_MODULE ITSMFT SurfaceLayout
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <vector>

#include "ITSMFTTracking/TraversalTopology.h"

namespace
{
using namespace o2::itsmft::tracking;

std::vector<SurfaceDescriptor> catalog(uint16_t count, SurfaceKind kind = SurfaceKind::Cylinder)
{
  std::vector<SurfaceDescriptor> result;
  for (uint16_t id = 0; id < count; ++id) {
    result.emplace_back(SurfaceId{id}, id, 0, kind);
  }
  return result;
}

std::vector<SurfaceId> order(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> result;
  for (uint16_t id = 0; id < count; ++id) {
    result.emplace_back(static_cast<uint16_t>(first + id));
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
} // namespace

BOOST_AUTO_TEST_CASE(SurfaceMaskCoversThirtyTwoGlobalSurfaces)
{
  SurfaceMask surfaces;
  surfaces.set(SurfaceId{0});
  surfaces.set(SurfaceId{16});
  surfaces.set(SurfaceId{31});
  BOOST_CHECK(surfaces.has(SurfaceId{0}));
  BOOST_CHECK(surfaces.has(SurfaceId{16}));
  BOOST_CHECK(surfaces.has(SurfaceId{31}));
  BOOST_CHECK_EQUAL(surfaces.count(), 3);
}

BOOST_AUTO_TEST_CASE(LayoutValidatesLimitsAndNonContiguousCatalogIds)
{
  const auto surfaces = catalog(33);
  const auto layout = SurfaceLayout{surfaces, makeSurfaceLayoutChain(order(0, 33))};
  BOOST_CHECK(layout.getError() == SurfaceLayoutError::TooManySurfaces);

  const std::vector<SurfaceDescriptor> sparse{{SurfaceId{0}, 0, 0, SurfaceKind::Cylinder},
                                               {SurfaceId{2}, 1, 0, SurfaceKind::Cylinder}};
  const std::vector<SurfaceId> sparseOrder{SurfaceId{0}, SurfaceId{2}};
  const auto valid = SurfaceLayout{sparse, makeSurfaceLayoutChain(sparseOrder)};
  BOOST_CHECK(valid.valid());
  BOOST_CHECK(valid.getSurfaceCatalog().getSurface(SurfaceId{2}).id == SurfaceId{2});
}

BOOST_AUTO_TEST_CASE(ComponentBoundariesAndKindIndependentCatalogs)
{
  const auto mixed = std::vector<SurfaceDescriptor>{{SurfaceId{0}, 0, 0, SurfaceKind::Cylinder},
                                                     {SurfaceId{1}, 1, 0, SurfaceKind::Cylinder},
                                                     {SurfaceId{2}, 0, 8, SurfaceKind::Disk},
                                                     {SurfaceId{3}, 1, 8, SurfaceKind::Disk}};
  SurfaceLayoutDefinition definition;
  definition.orderedSurfaces = order(0, 4);
  definition.componentOffsets = {0, 2};
  const auto layout = SurfaceLayout{mixed, std::move(definition)};
  BOOST_REQUIRE(layout.valid());
  BOOST_CHECK(layout.sameComponent(0, 1));
  BOOST_CHECK(!layout.sameComponent(1, 2));

  const auto topology = deriveTraversalTopology(layout);
  BOOST_REQUIRE(topology.ok());
  BOOST_CHECK_EQUAL(topology.topology->edges.size(), 2u);
  BOOST_CHECK(std::all_of(topology.topology->edges.begin(), topology.topology->edges.end(), [](const Edge& edge) {
    return edge.from.value() / 2 == edge.to.value() / 2;
  }));
}

BOOST_AUTO_TEST_CASE(HoleAndSeedPoliciesProduceSparseTopology)
{
  SurfaceLayoutDefinition definition;
  definition.orderedSurfaces = order(0, 4);
  definition.maxHoles = 1;
  definition.holeSurfaces = mask({1});
  definition.seedingSurfaces = mask({2});
  const std::vector<SurfaceDescriptor> surfaces = catalog(4);
  const auto layout = SurfaceLayout{surfaces, std::move(definition)};
  const auto result = deriveTraversalTopology(layout, mask({1}));
  BOOST_REQUIRE(result.ok());
  const auto& topology = *result.topology;
  BOOST_CHECK_EQUAL(topology.activeSurfaceList.size(), 3u);
  BOOST_CHECK_EQUAL(topology.edges.size(), 2u);
  BOOST_CHECK_EQUAL(topology.paths.size(), 1u);
  BOOST_CHECK(topology.edges[0].from == SurfaceId{0});
  BOOST_CHECK(topology.edges[0].to == SurfaceId{2});
  BOOST_REQUIRE_EQUAL(topology.roadStartPaths.size(), 1u);
  BOOST_CHECK(topology.getView(layout.getSurfaceCatalog()).getPath(topology.roadStartPaths.front()).first == EdgeId{0});
}

BOOST_AUTO_TEST_CASE(InvalidLayoutAndDisabledSurfaceDerivationIsTransactional)
{
  const auto surfaces = catalog(4);
  const auto layout = SurfaceLayout{surfaces, makeSurfaceLayoutChain(order(0, 4))};
  const auto invalidDisabled = deriveTraversalTopology(layout, mask({7}));
  BOOST_CHECK(!invalidDisabled.ok());
  BOOST_CHECK(!invalidDisabled.topology.has_value());
  BOOST_CHECK(invalidDisabled.error == TraversalTopologyError::DisabledSurfaceOutsideLayout);

  const auto invalidLayout = deriveTraversalTopology(SurfaceLayout{});
  BOOST_CHECK(!invalidLayout.ok());
  BOOST_CHECK(!invalidLayout.topology.has_value());
}
