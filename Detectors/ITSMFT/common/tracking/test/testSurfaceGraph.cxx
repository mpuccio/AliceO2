// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#define BOOST_TEST_MODULE ITSMFT SurfaceGraph
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <type_traits>
#include <vector>

#include "ITSMFTTracking/SurfaceGraphBuilder.h"

using namespace o2::itsmft::tracking;

namespace
{
SurfaceDescriptor surface(uint16_t id, SurfaceKind kind = SurfaceKind::Cylinder)
{
  return SurfaceDescriptor{SurfaceId{id}, id, 0, kind};
}

std::vector<SurfaceDescriptor> catalog()
{
  return {surface(0), surface(1), surface(2, SurfaceKind::Disk), surface(3, SurfaceKind::Disk)};
}

SurfaceGraph buildGraph(const std::vector<SurfaceDescriptor>& surfaces, std::vector<SurfaceId> order)
{
  SurfaceGraph graph{gsl::span<const SurfaceDescriptor>{surfaces}};
  graph.setOrderedSurfaces(std::move(order));
  BOOST_REQUIRE(graph.addEdge(Edge{SurfaceId{3}, SurfaceId{2}, {}, 0}).isValid());
  BOOST_REQUIRE(graph.addEdge(Edge{SurfaceId{1}, SurfaceId{0}, {}, 0}).isValid());
  BOOST_REQUIRE(graph.finalize());
  return graph;
}
} // namespace

static_assert(std::is_standard_layout_v<SurfaceGraphView>);
static_assert(std::is_trivially_copyable_v<SurfaceGraphView>);
static_assert(sizeof(SurfaceGraphView) == 96);
static_assert(alignof(SurfaceGraphView) == 8);

BOOST_AUTO_TEST_CASE(DefaultViewIsEmpty)
{
  const SurfaceGraphView view{};
  BOOST_CHECK(view.surfaces == nullptr);
  BOOST_CHECK_EQUAL(view.nSurfaces, 0u);
  BOOST_CHECK(view.orderedSurfaces == nullptr);
}

BOOST_AUTO_TEST_CASE(ViewCarriesOneOwnerCatalogOrderAndTopology)
{
  const auto surfaces = catalog();
  auto graph = buildGraph(surfaces, {SurfaceId{3}, SurfaceId{1}, SurfaceId{2}, SurfaceId{0}});
  const auto view = graph.getView();

  BOOST_REQUIRE(graph.valid());
  BOOST_CHECK(view.surfaces == graph.getSurfaceCatalog().surfaces);
  BOOST_CHECK(view.orderedSurfaces == graph.getOrderedSurfaces().data());
  BOOST_CHECK_EQUAL(view.nSurfaces, 4u);
  BOOST_CHECK_EQUAL(view.nOrderedSurfaces, 4u);
  BOOST_CHECK_EQUAL(view.cylinderSurfaces.value(), 0x3u);
  BOOST_CHECK_EQUAL(view.diskSurfaces.value(), 0xcu);
  BOOST_CHECK(view.orderedSurfaces[0] == SurfaceId{3});
  BOOST_CHECK(view.orderedSurfaces[1] == SurfaceId{1});
  BOOST_CHECK_EQUAL(view.nEdges, 2u);
  BOOST_CHECK_EQUAL(view.nCells, 0u);
  BOOST_CHECK(view.getEdge(EdgeId{0}).from == SurfaceId{3});
  BOOST_CHECK(view.getEdge(EdgeId{0}).to == SurfaceId{2});
}

BOOST_AUTO_TEST_CASE(ViewIsDeviceSafePOD)
{
  const auto surfaces = catalog();
  auto graph = buildGraph(surfaces, {SurfaceId{0}, SurfaceId{1}, SurfaceId{2}, SurfaceId{3}});
  const auto view = graph.getView();
  const auto copy = view;
  BOOST_CHECK(copy.surfaces == view.surfaces);
  BOOST_CHECK(copy.edges == view.edges);
  BOOST_CHECK(copy.cells == view.cells);
}

BOOST_AUTO_TEST_CASE(InvalidGraphReturnsEmptyView)
{
  SurfaceGraph topology{3};
  BOOST_REQUIRE(topology.finalize());
  const std::vector<SurfaceDescriptor> surfaces{surface(0), surface(1)};
  SurfaceGraph graph{gsl::span<const SurfaceDescriptor>{surfaces}, std::move(topology)};
  BOOST_CHECK(!graph.valid());
  const auto view = graph.getView();
  BOOST_CHECK(view.surfaces == nullptr);
  BOOST_CHECK_EQUAL(view.nSurfaces, 0u);
}

BOOST_AUTO_TEST_CASE(NonContiguousGlobalIdsUseOneGraphLookup)
{
  const std::vector<SurfaceDescriptor> surfaces{surface(1), surface(4), surface(7)};
  const std::vector<SurfaceId> ordered{SurfaceId{7}, SurfaceId{1}, SurfaceId{4}};
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())},
                              makeSurfaceChain(ordered, 1, SurfaceMask{uint32_t{1} << 1}, SurfaceMask{uint32_t{1} << 7})};
  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());
  const auto view = result.graph->getView();
  BOOST_CHECK_EQUAL(view.nSurfaces, 3u);
  BOOST_CHECK(view.getSurface(SurfaceId{7}).id == SurfaceId{7});
  BOOST_CHECK(view.getSurface(SurfaceId{1}).id == SurfaceId{1});
  BOOST_CHECK(view.orderedSurfaces[0] == SurfaceId{7});
  BOOST_CHECK(view.seedingSurfaces.has(SurfaceId{7}));
  BOOST_CHECK(view.getEdge(EdgeId{0}).from == SurfaceId{7});
  BOOST_CHECK(view.getEdge(EdgeId{0}).to == SurfaceId{1});
}
