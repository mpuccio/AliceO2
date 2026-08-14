// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details.

#define BOOST_TEST_MODULE ITSMFT TraversalWorkspacePlan
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <vector>

#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/Tracker.h"

#include "TraversalTestSupport.h"

using namespace o2::itsmft::tracking;

namespace
{
SurfaceGraph makeGraph()
{
  std::array<SurfaceDescriptor, 4> catalog{};
  std::vector<SurfaceId> ordered;
  for (uint16_t id = 0; id < catalog.size(); ++id) {
    catalog[id] = SurfaceDescriptor{SurfaceId{id}, id, 0, SurfaceKind::Cylinder};
    ordered.push_back(SurfaceId{id});
  }
  const auto definition = makeSurfaceChain(ordered, 1, SurfaceMask{uint32_t{1} << 1}, SurfaceMask{uint32_t{1} << 3});
  const auto result = SurfaceGraphBuilder{SurfaceCatalogView{catalog.data(), static_cast<uint32_t>(catalog.size())}, definition}.build();
  BOOST_REQUIRE(result.ok());
  return *result.graph;
}
} // namespace

BOOST_AUTO_TEST_CASE(PassWorkspaceDerivesSelectedTopologyFromTheStaticGraph)
{
  const auto graph = makeGraph();
  Tracker tracker;
  TraversalWorkspace first;
  TraversalWorkspace second;

  TrackerTestAccess::preparePlan(tracker, first, graph.getView());
  TrackerTestAccess::preparePlan(tracker, second, graph.getView(), LayerMask{uint32_t{1} << 1});

  BOOST_CHECK_EQUAL(first.orderedSurfaces.size(), 4u);
  BOOST_CHECK_EQUAL(first.edges.size(), 4u);
  BOOST_CHECK_EQUAL(first.cells.size(), 3u);
  BOOST_CHECK_EQUAL(second.orderedSurfaces.size(), 4u);
  BOOST_CHECK_EQUAL(second.activeSurfaces.count(), 3);
  BOOST_CHECK_EQUAL(second.edges.size(), 2u);
  BOOST_CHECK_EQUAL(second.cells.size(), 1u);
  BOOST_CHECK(second.getSurfaceSlot(SurfaceId{1}).has_value());
  BOOST_CHECK(second.getEdgeSlot(EdgeId{1}).has_value()); // 0 -> 2, skipping disabled surface 1
  BOOST_CHECK_EQUAL(second.scheduledCells.size(), second.cells.size());
}

BOOST_AUTO_TEST_CASE(FailedPlanPreparationLeavesTheWorkspaceInvalidAndEmpty)
{
  Tracker tracker;
  TraversalWorkspace workspace;
  const SurfaceGraphView invalid{};
  BOOST_CHECK_THROW(TrackerTestAccess::preparePlan(tracker, workspace, invalid), TraversalException);
  BOOST_CHECK(!workspace.valid);
  BOOST_CHECK(workspace.orderedSurfaces.empty());
  BOOST_CHECK(workspace.edges.empty());
  BOOST_CHECK(workspace.cells.empty());
}
