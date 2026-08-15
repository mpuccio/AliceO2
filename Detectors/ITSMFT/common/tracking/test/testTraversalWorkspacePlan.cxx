// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details.

#define BOOST_TEST_MODULE ITSMFT TraversalWorkspacePlan
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "ITSMFTTracking/Tracker.h"

#include "TraversalTestSupport.h"

using namespace o2::itsmft::tracking;

namespace
{
SurfaceLayout makeLayout()
{
  std::array<SurfaceDescriptor, 4> catalog{};
  std::vector<SurfaceId> ordered;
  for (uint16_t id = 0; id < catalog.size(); ++id) {
    catalog[id] = SurfaceDescriptor{SurfaceId{id}, id, 0, SurfaceKind::Cylinder};
    ordered.push_back(SurfaceId{id});
  }
  SurfaceLayoutDefinition definition;
  definition.orderedSurfaces = std::move(ordered);
  definition.maxHoles = 1;
  definition.holeSurfaces = SurfaceMask{uint32_t{1} << 1};
  definition.seedingSurfaces = SurfaceMask{uint32_t{1} << 3};
  return SurfaceLayout{gsl::span<const SurfaceDescriptor>{catalog.data(), catalog.size()}, std::move(definition)};
}
} // namespace

BOOST_AUTO_TEST_CASE(PassWorkspaceDerivesSelectedTopologyFromTheStaticLayout)
{
  const auto layout = makeLayout();
  Tracker tracker;
  TraversalWorkspace first;
  TraversalWorkspace second;

  TrackerTestAccess::preparePlan(tracker, first, layout);
  TrackerTestAccess::preparePlan(tracker, second, layout, SurfaceMask{uint32_t{1} << 1});

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
  const SurfaceLayout invalid{};
  BOOST_CHECK_THROW(TrackerTestAccess::preparePlan(tracker, workspace, invalid), TraversalException);
  BOOST_CHECK(!workspace.valid);
  BOOST_CHECK(workspace.orderedSurfaces.empty());
  BOOST_CHECK(workspace.edges.empty());
  BOOST_CHECK(workspace.cells.empty());
  BOOST_CHECK(workspace.topology.edges.empty());
  BOOST_CHECK(workspace.topology.paths.empty());
  BOOST_CHECK_EQUAL(workspace.getTopologyView().nEdges, 0u);
  BOOST_CHECK_EQUAL(workspace.getTopologyView().nPaths, 0u);
}

BOOST_AUTO_TEST_CASE(KernelViewBorrowsWorkspaceTopology)
{
  TraversalWorkspace workspace;
  workspace.topology.orderedSurfaces = {SurfaceId{0}, SurfaceId{1}};
  workspace.topology.activeSurfaceList = workspace.topology.orderedSurfaces;
  workspace.topology.surfacePositionById.assign(MaxLayoutSurfaces, -1);
  workspace.topology.surfacePositionById[0] = 0;
  workspace.topology.surfacePositionById[1] = 1;
  workspace.topology.activeSurfaces = SurfaceMask{uint32_t{0b11}};
  workspace.topology.edges.push_back(Edge{SurfaceId{0}, SurfaceId{1}});
  workspace.topology.paths.push_back(CellPath{EdgeId{0}, EdgeId{0}});
  workspace.topology.pathsByFirstEdgeOffsets = {0, 1};
  workspace.topology.pathsByFirstEdge.push_back(CellPathId{0});
  workspace.topology.scheduledPaths.push_back(CellPathId{0});
  workspace.topology.roadStartPaths.push_back(CellPathId{0});
  workspace.topology.roadStartComponentOffsets = {0, 1};

  const auto view = workspace.getTopologyView();
  BOOST_CHECK(view.orderedSurfaces == workspace.topology.orderedSurfaces.data());
  BOOST_CHECK(view.edges == workspace.topology.edges.data());
  BOOST_CHECK(view.paths == workspace.topology.paths.data());
  BOOST_CHECK(view.pathsByFirstEdgeOffsets == workspace.topology.pathsByFirstEdgeOffsets.data());
  BOOST_CHECK(view.pathsByFirstEdge == workspace.topology.pathsByFirstEdge.data());
  BOOST_CHECK(view.scheduledPaths == workspace.topology.scheduledPaths.data());
  BOOST_CHECK(view.roadStartPaths == workspace.topology.roadStartPaths.data());
  BOOST_CHECK_EQUAL(view.nEdges, 1u);
  BOOST_CHECK_EQUAL(view.nPaths, 1u);
}

BOOST_AUTO_TEST_CASE(ResetClearsAndInvalidatesOwnedTopology)
{
  TraversalWorkspace workspace;
  workspace.valid = true;
  workspace.topology.edges.push_back(Edge{SurfaceId{0}, SurfaceId{1}});
  workspace.topology.paths.push_back(CellPath{EdgeId{0}, EdgeId{0}});
  workspace.topology.activeSurfaces = SurfaceMask{uint32_t{0b11}};

  workspace.reset(nullptr);

  BOOST_CHECK(!workspace.valid);
  BOOST_CHECK(workspace.topology.edges.empty());
  BOOST_CHECK(workspace.topology.paths.empty());
  BOOST_CHECK(workspace.topology.activeSurfaces.empty());
  const auto view = workspace.getTopologyView();
  BOOST_CHECK_EQUAL(view.nEdges, 0u);
  BOOST_CHECK_EQUAL(view.nPaths, 0u);
  BOOST_CHECK(view.activeSurfaces.empty());
}

BOOST_AUTO_TEST_CASE(ExecutionSourcesUseCurrentTopology)
{
  namespace fs = std::filesystem;
  const auto trackingRoot = fs::path{__FILE__}.parent_path().parent_path();
  const std::array<fs::path, 7> sources{
    trackingRoot / "include/ITSMFTTracking/TrackerTraits.h",
    trackingRoot / "src/TrackerTraits.cxx",
    trackingRoot / "include/ITSMFTTracking/detail/SurfaceTrackingScratch.h",
    trackingRoot / "src/SurfaceTrackingScratch.cxx",
    trackingRoot / "include/ITSMFTTracking/TraversalTopology.h",
    trackingRoot / "src/TraversalTopology.cxx",
    trackingRoot / "src/Tracker.cxx"};
  for (const auto& path : sources) {
    std::ifstream input{path};
    const std::string source{std::istreambuf_iterator<char>{input}, {}};
    BOOST_REQUIRE_MESSAGE(!source.empty(), "cannot read execution source " << path.string());
    BOOST_CHECK_MESSAGE(source.find("template <typename LayoutView>") == std::string::npos,
                        "workspace view retains ignored LayoutView constructor: " << path.string());
    BOOST_CHECK_MESSAGE(source.find("SurfaceMask skippedSurfaces") == std::string::npos,
                        "lean Edge contract retains skipped-surface storage: " << path.string());
    BOOST_CHECK_MESSAGE(source.find("uint16_t flags") == std::string::npos,
                        "lean Edge contract retains unused flags storage: " << path.string());
  }
}
