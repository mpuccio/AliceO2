// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT SurfaceGraphBuilder
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include "ITSMFTTracking/SurfaceGraphBuilder.h"

namespace
{
using namespace o2::itsmft::tracking;

SurfaceDescriptor surface(uint16_t id, SurfaceKind kind = SurfaceKind::Cylinder)
{
  return SurfaceDescriptor{SurfaceId{id}, id, 0, kind};
}

std::vector<SurfaceDescriptor> denseCatalog(uint16_t count, SurfaceKind kind = SurfaceKind::Cylinder)
{
  std::vector<SurfaceDescriptor> catalog;
  catalog.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    catalog.push_back(surface(i, kind));
  }
  return catalog;
}

// SurfaceGraphBuilder borrows (SurfaceCatalogView), so every call site
// below binds its catalog to a named local first -- the view must not
// outlive the vector it points into, but only needs to survive the single
// build() call (SurfaceGraph's own borrowed-span constructor never
// retains it past construction).
SurfaceCatalogView asView(const std::vector<SurfaceDescriptor>& catalog)
{
  return SurfaceCatalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
}

std::vector<SurfaceId> orderedIds(std::initializer_list<uint16_t> ids)
{
  std::vector<SurfaceId> out;
  out.reserve(ids.size());
  for (auto id : ids) {
    out.push_back(SurfaceId{id});
  }
  return out;
}

SurfaceMask maskOf(std::initializer_list<uint16_t> ids)
{
  SurfaceMask mask;
  for (auto id : ids) {
    mask.set(SurfaceId{id});
  }
  return mask;
}

void checkCsrConsistency(const SurfaceGraphView& view)
{
  uint32_t total = 0;
  for (uint32_t t = 0; t < view.nTransitions; ++t) {
    const auto range = view.getCellsStartingWithTransition(TransitionId{static_cast<uint16_t>(t)});
    total += range.getEntries();
    for (uint32_t k = 0; k < range.getEntries(); ++k) {
      const auto cellId = view.cellsByFirstTransition[range.getFirstEntry() + k];
      BOOST_CHECK_EQUAL(view.getCell(cellId).firstTransition.value(), t);
    }
  }
  BOOST_CHECK_EQUAL(total, view.nCells);
}

void checkSparseTopology(const SurfaceGraphView& view)
{
  BOOST_REQUIRE_GT(view.nTransitions, 0u);
  BOOST_REQUIRE_GT(view.nCells, 0u);
  for (uint32_t t = 0; t < view.nTransitions; ++t) {
    const auto& transition = view.getTransition(TransitionId{static_cast<uint16_t>(t)});
    BOOST_CHECK(transition.from.isValid());
    BOOST_CHECK(transition.to.isValid());
    BOOST_CHECK(transition.from != transition.to);
  }
  for (uint32_t c = 0; c < view.nCells; ++c) {
    const auto& cell = view.getCell(CellTopologyId{static_cast<uint16_t>(c)});
    BOOST_CHECK(cell.firstTransition.value() < view.nTransitions);
    BOOST_CHECK(cell.secondTransition.value() < view.nTransitions);
    BOOST_CHECK_EQUAL(cell.hitSurfaces.count(), 3);
  }

  checkCsrConsistency(view);
}
} // namespace

BOOST_AUTO_TEST_CASE(TraversalFollowsSuppliedOrderNotNumericSurfaceId)
{
  // Catalog ids are dense (0..4) but the physical traversal order is the
  // deliberately non-monotonic chain 3 -> 1 -> 4 -> 0; surface 2 is present
  // in the catalog but not activated by this subgraph.
  const auto catalog = denseCatalog(5);
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({3, 1, 4, 0}), 0, SurfaceMask{}, SurfaceMask{}});

  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());
  const auto view = result.graph->getView();

  BOOST_CHECK_EQUAL(view.nTransitions, 3u);
  BOOST_CHECK(view.getTransition(TransitionId{0}).from == SurfaceId{3});
  BOOST_CHECK(view.getTransition(TransitionId{0}).to == SurfaceId{1});
  BOOST_CHECK(view.getTransition(TransitionId{1}).from == SurfaceId{1});
  BOOST_CHECK(view.getTransition(TransitionId{1}).to == SurfaceId{4});
  BOOST_CHECK(view.getTransition(TransitionId{2}).from == SurfaceId{4});
  BOOST_CHECK(view.getTransition(TransitionId{2}).to == SurfaceId{0});

  // Numeric-SurfaceId-adjacent pairs that are not adjacent in the supplied
  // traversal must not appear.
  for (uint32_t t = 0; t < view.nTransitions; ++t) {
    const auto& transition = view.getTransition(TransitionId{static_cast<uint16_t>(t)});
    BOOST_CHECK(!(transition.from == SurfaceId{0} && transition.to == SurfaceId{1}));
    BOOST_CHECK(!(transition.from == SurfaceId{1} && transition.to == SurfaceId{2}));
    BOOST_CHECK(!(transition.from == SurfaceId{2} && transition.to == SurfaceId{3}));
  }

  BOOST_CHECK_EQUAL(view.nCells, 2u);
  BOOST_CHECK_EQUAL(view.getCell(CellTopologyId{0}).hitSurfaces.value(), maskOf({3, 1, 4}).value());
  BOOST_CHECK_EQUAL(view.getCell(CellTopologyId{1}).hitSurfaces.value(), maskOf({1, 4, 0}).value());

  checkCsrConsistency(view);
}

BOOST_AUTO_TEST_CASE(SparseTopologySevenNoHoles)
{
  const auto catalog = denseCatalog(7);
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 1, 2, 3, 4, 5, 6}), 0, SurfaceMask{}, SurfaceMask{}});
  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());

  checkSparseTopology(result.graph->getView());
}

BOOST_AUTO_TEST_CASE(SparseTopologySevenSingleAllowedHole)
{
  const auto catalog = denseCatalog(7);
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 1, 2, 3, 4, 5, 6}), 1, maskOf({3}), SurfaceMask{}});
  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());

  const auto view = result.graph->getView();
  checkSparseTopology(view);
  BOOST_CHECK_EQUAL(view.nTransitions, 7u);
  BOOST_CHECK_EQUAL(view.nCells, 7u);
}

BOOST_AUTO_TEST_CASE(SparseTopologyTenNoHoles)
{
  const auto catalog = denseCatalog(10, SurfaceKind::Disk);
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), 0, SurfaceMask{}, SurfaceMask{}});
  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());
  checkSparseTopology(result.graph->getView());
}

BOOST_AUTO_TEST_CASE(SparseTopologyTenSingleAllowedHole)
{
  const auto catalog = denseCatalog(10, SurfaceKind::Disk);
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), 1, maskOf({5}), SurfaceMask{}});
  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());

  const auto view = result.graph->getView();
  checkSparseTopology(view);
  BOOST_CHECK_EQUAL(view.nTransitions, 10u);
  BOOST_CHECK_EQUAL(view.nCells, 10u);
}

BOOST_AUTO_TEST_CASE(SingleCallDisconnectedCylinderAndDiskLayout)
{
  // One builder, one build() call, two disjoint subgraphs: a 7-cylinder
  // ITS-like stack (ids 0-6) and a 10-disk MFT-like stack (ids 7-16) over a
  // single 17-surface catalog.
  std::vector<SurfaceDescriptor> catalog = denseCatalog(7, SurfaceKind::Cylinder);
  for (uint16_t id = 7; id < 17; ++id) {
    catalog.push_back(surface(id, SurfaceKind::Disk));
  }

  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 1, 2, 3, 4, 5, 6}), 0, SurfaceMask{}, maskOf({0})});
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({7, 8, 9, 10, 11, 12, 13, 14, 15, 16}), 0, SurfaceMask{}, maskOf({7})});

  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(result.graph->valid());

  const auto masks = computeSurfaceKindMasks(catalog);
  const auto layoutView = result.graph->getView();
  BOOST_CHECK_EQUAL(layoutView.nSurfaces, 17u);
  BOOST_CHECK_EQUAL(layoutView.cylinderSurfaces.value(), 0x7fu); // bits 0-6
  BOOST_CHECK_EQUAL(layoutView.diskSurfaces.value(), 0x1ff80u);  // bits 7-16

  const auto view = layoutView;
  BOOST_CHECK_EQUAL(view.nTransitions, 15u); // 6 + 9
  BOOST_CHECK_EQUAL(view.nCells, 13u);       // 5 + 8

  // Disconnected: the last ITS transition (5 -> 6) has no successors, and no
  // cell crosses the cylinder/disk boundary.
  BOOST_CHECK_EQUAL(view.getCellsStartingWithTransition(TransitionId{5}).getEntries(), 0u);
  for (uint32_t c = 0; c < view.nCells; ++c) {
    const auto& cell = view.getCell(CellTopologyId{static_cast<uint16_t>(c)});
    const bool allCylinder = cell.hitSurfaces.isSubsetOf(layoutView.cylinderSurfaces);
    const bool allDisk = cell.hitSurfaces.isSubsetOf(layoutView.diskSurfaces);
    BOOST_CHECK(allCylinder || allDisk);
  }

  checkCsrConsistency(view);
}

BOOST_AUTO_TEST_CASE(CatalogAboveThirtyTwoSurfacesIsRejected)
{
  const auto catalog = denseCatalog(33);
  SurfaceGraphBuilder builder{asView(catalog)};
  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == SurfaceGraphBuildError::TopologyRejected);
  BOOST_CHECK(result.topologyError == SurfaceGraphTopologyError::InvalidSurfaceCount);
}

BOOST_AUTO_TEST_CASE(NonContiguousCatalogSurfaceIdsAreSupported)
{
  std::vector<SurfaceDescriptor> catalog;
  catalog.push_back(surface(0));
  catalog.push_back(surface(2));
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 2}), 0, SurfaceMask{}, maskOf({2})});

  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result.graph->getView().nSurfaces, 2u);
  BOOST_CHECK(result.graph->getView().getSurface(SurfaceId{2}).id == SurfaceId{2});
}

BOOST_AUTO_TEST_CASE(InvalidSubgraphSurfaceIdsAreRejected)
{
  {
    // Out-of-range id.
    const auto catalog = denseCatalog(3);
    SurfaceGraphBuilder builder{asView(catalog)};
    builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 5}), 0, SurfaceMask{}, SurfaceMask{}});
    const auto result = builder.build();
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == SurfaceGraphBuildError::InvalidSubgraphSurfaceId);
  }
  {
    // The invalid-sentinel id.
    const auto catalog = denseCatalog(3);
    SurfaceGraphBuilder builder{asView(catalog)};
    std::vector<SurfaceId> ordered{SurfaceId{0}, SurfaceId::invalid()};
    builder.addSubgraph(SurfaceGraphSubgraph{std::move(ordered), 0, SurfaceMask{}, SurfaceMask{}});
    const auto result = builder.build();
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == SurfaceGraphBuildError::InvalidSubgraphSurfaceId);
  }
}

BOOST_AUTO_TEST_CASE(DuplicateSurfaceWithinSubgraphIsRejected)
{
  const auto catalog = denseCatalog(3);
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 1, 0}), 0, SurfaceMask{}, SurfaceMask{}});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == SurfaceGraphBuildError::DuplicateSurfaceInSubgraph);
}

BOOST_AUTO_TEST_CASE(SurfaceDuplicatedAcrossSubgraphsIsRejected)
{
  const auto catalog = denseCatalog(4);
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 1}), 0, SurfaceMask{}, SurfaceMask{}});
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({1, 2}), 0, SurfaceMask{}, SurfaceMask{}});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == SurfaceGraphBuildError::SurfaceDuplicatedAcrossSubgraphs);
}

BOOST_AUTO_TEST_CASE(SubgraphMixingSurfaceKindsIsRejected)
{
  // M4 (ADR 0007 decision 8): a subgraph no longer asserts an external policy
  // tag the builder validates against the catalog -- its expected SurfaceKind
  // is derived from its own first surface, and every other surface in the
  // subgraph must match it. Catalog surface 0 is Cylinder, surface 1 is Disk;
  // a subgraph spanning both is rejected on that basis alone.
  std::vector<SurfaceDescriptor> catalog{surface(0, SurfaceKind::Cylinder), surface(1, SurfaceKind::Disk)};
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 1}), 0, SurfaceMask{}, SurfaceMask{}});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == SurfaceGraphBuildError::GraphRejected);
  BOOST_CHECK(result.graphError == SurfaceGraphError::PolicySurfaceKindMismatch);
}

BOOST_AUTO_TEST_CASE(SingletonSubgraphsOfEitherKindAreAccepted)
{
  // A singleton subgraph never calls addTransition, so the (former) tag
  // validation this used to exercise up front no longer has anything to
  // check -- a single surface's own kind is trivially internally consistent
  // with itself, for either SurfaceKind.
  {
    const auto catalog = denseCatalog(1, SurfaceKind::Disk);
    SurfaceGraphBuilder builder{asView(catalog)};
    builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0}), 0, SurfaceMask{}, SurfaceMask{}});
    BOOST_CHECK(builder.build().ok());
  }
  {
    const auto catalog = denseCatalog(1, SurfaceKind::Cylinder);
    SurfaceGraphBuilder builder{asView(catalog)};
    builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0}), 0, SurfaceMask{}, SurfaceMask{}});
    BOOST_CHECK(builder.build().ok());
  }
}

BOOST_AUTO_TEST_CASE(NegativeMaxHolesIsRejected)
{
  // Explicit contract: a negative maxHoles is rejected outright rather than
  // silently normalized to zero.
  const auto catalog = denseCatalog(2);
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 1}), -1, SurfaceMask{}, SurfaceMask{}});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == SurfaceGraphBuildError::NegativeMaxHoles);
}

BOOST_AUTO_TEST_CASE(HoleSurfacesOutsideSubgraphAreRejected)
{
  const auto catalog = denseCatalog(3);
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 1}), 1, maskOf({2}), SurfaceMask{}});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == SurfaceGraphBuildError::HoleSurfacesOutsideSubgraph);
}

BOOST_AUTO_TEST_CASE(SeedingSurfacesOutsideSubgraphAreRejected)
{
  const auto catalog = denseCatalog(3);
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0, 1}), 0, SurfaceMask{}, maskOf({2})});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == SurfaceGraphBuildError::SeedingSurfacesOutsideSubgraph);
}

BOOST_AUTO_TEST_CASE(EmptyCatalogWithNoSubgraphsIsATrivialValidLayout)
{
  SurfaceGraphBuilder builder{SurfaceCatalogView{}};
  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result.graph->getView().nSurfaces, 0u);
  BOOST_CHECK_EQUAL(result.graph->getView().nTransitions, 0u);
}

BOOST_AUTO_TEST_CASE(EmptySubgraphIsRejected)
{
  const auto catalog = denseCatalog(2);
  SurfaceGraphBuilder builder{asView(catalog)};
  builder.addSubgraph(SurfaceGraphSubgraph{{}, 0, SurfaceMask{}, SurfaceMask{}});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == SurfaceGraphBuildError::EmptySubgraph);
}

BOOST_AUTO_TEST_CASE(EmptyCatalogWithNonEmptySubgraphIsRejected)
{
  SurfaceGraphBuilder builder{SurfaceCatalogView{}};
  builder.addSubgraph(SurfaceGraphSubgraph{orderedIds({0}), 0, SurfaceMask{}, SurfaceMask{}});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == SurfaceGraphBuildError::InvalidSubgraphSurfaceId);
}
