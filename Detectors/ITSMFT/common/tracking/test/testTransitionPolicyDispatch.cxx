// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT TransitionPolicyDispatch
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <map>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/TransitionPolicyDispatch.h"

using namespace o2::itsmft::tracking;

namespace
{

SurfaceTransition adjacent(uint16_t from, uint16_t to, TransitionPolicyTag policyTag)
{
  return SurfaceTransition{SurfaceId{from}, SurfaceId{to}, SurfaceMask{}, policyTag, 0};
}

SurfaceDescriptor surface(uint16_t id, SurfaceKind kind)
{
  return SurfaceDescriptor{SurfaceId{id}, id, 0, kind};
}

/// A chain of `nSurfaces` surfaces of one kind, adjacent-only transitions and
/// cells (no holes) -- ITS-like when Cylinder/CylinderCylinder, MFT-like when
/// Disk/DiskDisk. Mirrors the fixture style already accepted in
/// testSurfaceLayout.cxx.
DetectorLayout buildChainLayout(uint16_t nSurfaces, SurfaceKind kind, TransitionPolicyTag tag)
{
  SparseTrackingTopology topology{nSurfaces};
  std::vector<TransitionId> transitions;
  for (uint16_t s = 0; s + 1 < nSurfaces; ++s) {
    transitions.push_back(topology.addTransition(adjacent(s, s + 1, tag)));
  }
  for (size_t t = 0; t + 1 < transitions.size(); ++t) {
    topology.addCell(transitions[t], transitions[t + 1]);
  }
  BOOST_REQUIRE(topology.finalize());

  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t s = 0; s < nSurfaces; ++s) {
    surfaces.push_back(surface(s, kind));
  }
  return DetectorLayout{std::move(surfaces), std::move(topology)};
}

/// A disconnected combined layout: an ITS-like cylinder chain over the first
/// `nCylinders` surfaces and an MFT-like disk chain over the next `nDisks`
/// surfaces, in one shared global surface-id space (Architecture.md §8, the
/// first combined milestone).
DetectorLayout buildCombinedDisconnectedLayout(uint16_t nCylinders, uint16_t nDisks)
{
  const uint16_t total = nCylinders + nDisks;
  SparseTrackingTopology topology{total};
  std::vector<TransitionId> cylinderTransitions;
  std::vector<TransitionId> diskTransitions;
  for (uint16_t s = 0; s + 1 < nCylinders; ++s) {
    cylinderTransitions.push_back(topology.addTransition(adjacent(s, s + 1, TransitionPolicyTag::CylinderCylinder)));
  }
  for (uint16_t s = nCylinders; s + 1 < total; ++s) {
    diskTransitions.push_back(topology.addTransition(adjacent(s, s + 1, TransitionPolicyTag::DiskDisk)));
  }
  for (size_t t = 0; t + 1 < cylinderTransitions.size(); ++t) {
    topology.addCell(cylinderTransitions[t], cylinderTransitions[t + 1]);
  }
  for (size_t t = 0; t + 1 < diskTransitions.size(); ++t) {
    topology.addCell(diskTransitions[t], diskTransitions[t + 1]);
  }
  BOOST_REQUIRE(topology.finalize());

  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t s = 0; s < nCylinders; ++s) {
    surfaces.push_back(surface(s, SurfaceKind::Cylinder));
  }
  for (uint16_t s = nCylinders; s < total; ++s) {
    surfaces.push_back(surface(s, SurfaceKind::Disk));
  }
  return DetectorLayout{std::move(surfaces), std::move(topology)};
}

SurfaceMask maskOf(uint16_t id)
{
  SurfaceMask mask;
  mask.set(SurfaceId{id});
  return mask;
}

DetectorLayout buildIdentityLayout(uint16_t nSurfaces, SurfaceKind kind, TransitionPolicyTag tag, int maxHoles, uint16_t hole)
{
  std::vector<SurfaceDescriptor> surfaces;
  std::vector<SurfaceId> order;
  for (uint16_t id = 0; id < nSurfaces; ++id) {
    surfaces.push_back(surface(id, kind));
    order.push_back(SurfaceId{id});
  }
  DetectorLayoutBuilder builder{std::move(surfaces)};
  const auto result = builder.addSubgraph(DetectorLayoutSubgraph{std::move(order), maxHoles, maxHoles ? maskOf(hole) : SurfaceMask{}, SurfaceMask{}, tag}).build();
  BOOST_REQUIRE(result.ok());
  return std::move(*result.layout);
}

void checkIdentitySchedule(uint16_t nSurfaces, SurfaceKind kind, TransitionPolicyTag tag, int maxHoles, uint16_t hole)
{
  const auto layout = buildIdentityLayout(nSurfaces, kind, tag, maxHoles, hole);
  TransitionPolicyGrouping grouping{layout.getView()};
  BOOST_REQUIRE(grouping.valid());
  const auto scheduled = grouping.scheduledCellsForTag(tag);
  BOOST_REQUIRE_EQUAL(scheduled.size(), layout.getTopology().getCells().size());
  for (uint16_t id = 0; id < scheduled.size(); ++id) {
    BOOST_CHECK(scheduled[id] == CellTopologyId{id});
  }
}

/// Records which policy tags were dispatched and how much work each carried.
/// Its call operator is a template so the compiler selects the specialized
/// policy behaviour at compile time (D007); it never inspects a detector ID,
/// so the same visitor works unmodified for ITS-like, MFT-like, and combined
/// layouts.
struct RecordingVisitor {
  std::map<TransitionPolicyTag, size_t> transitionCounts;
  std::map<TransitionPolicyTag, size_t> cellCounts;
  std::map<TransitionPolicyTag, StateFamily> families;
  int calls{0};

  template <typename Traits>
  void operator()(Traits, gsl::span<const TransitionId> transitions, gsl::span<const CellTopologyId> cells)
  {
    ++calls;
    transitionCounts[Traits::Tag] = transitions.size();
    cellCounts[Traits::Tag] = cells.size();
    families[Traits::Tag] = Traits::Family;
  }
};

} // namespace

BOOST_AUTO_TEST_CASE(GroupingIsEmptyForAnEmptyLayoutView)
{
  TransitionPolicyGrouping grouping{DetectorLayoutView{}};
  BOOST_CHECK(grouping.valid());
  BOOST_CHECK(!grouping.hasTag(TransitionPolicyTag::CylinderCylinder));
  BOOST_CHECK(!grouping.hasTag(TransitionPolicyTag::DiskDisk));
  BOOST_CHECK_EQUAL(grouping.transitionsForTag(TransitionPolicyTag::CylinderCylinder).size(), 0u);
  BOOST_CHECK_EQUAL(grouping.cellsForTag(TransitionPolicyTag::DiskDisk).size(), 0u);
}

BOOST_AUTO_TEST_CASE(IdentityLayoutsHaveExactLegacyCellSchedule)
{
  checkIdentitySchedule(7, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder, 0, 0);
  checkIdentitySchedule(7, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder, 1, 3);
  checkIdentitySchedule(10, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk, 0, 0);
  checkIdentitySchedule(10, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk, 1, 5);
}

BOOST_AUTO_TEST_CASE(NonMonotonicSurfaceIdsFollowGraphRank)
{
  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t id = 0; id < 5; ++id) {
    surfaces.push_back(surface(id, SurfaceKind::Cylinder));
  }
  DetectorLayoutBuilder builder{std::move(surfaces)};
  auto result = builder.addSubgraph(DetectorLayoutSubgraph{{SurfaceId{3}, SurfaceId{1}, SurfaceId{4}, SurfaceId{0}}, 0, {}, {}, TransitionPolicyTag::CylinderCylinder}).build();
  BOOST_REQUIRE(result.ok());
  TransitionPolicyGrouping grouping{result.layout->getView()};
  BOOST_REQUIRE(grouping.valid());
  const auto scheduled = grouping.scheduledCellsForTag(TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE_EQUAL(scheduled.size(), 2u);
  BOOST_CHECK(scheduled[0] == CellTopologyId{0}); // target surface 4, graph rank 2
  BOOST_CHECK(scheduled[1] == CellTopologyId{1}); // target surface 0, graph rank 3
}

BOOST_AUTO_TEST_CASE(DisconnectedComponentsAreScheduledDeterministically)
{
  const auto layout = buildCombinedDisconnectedLayout(7, 10);
  for (int repeat = 0; repeat < 8; ++repeat) {
    TransitionPolicyGrouping grouping{layout.getView()};
    BOOST_REQUIRE(grouping.valid());
    const auto cylinder = grouping.scheduledCellsForTag(TransitionPolicyTag::CylinderCylinder);
    const auto disk = grouping.scheduledCellsForTag(TransitionPolicyTag::DiskDisk);
    BOOST_REQUIRE_EQUAL(cylinder.size(), 5u);
    BOOST_REQUIRE_EQUAL(disk.size(), 8u);
    for (uint16_t i = 0; i < cylinder.size(); ++i) {
      BOOST_CHECK(cylinder[i] == CellTopologyId{i});
    }
    for (uint16_t i = 0; i < disk.size(); ++i) {
      BOOST_CHECK(disk[i] == CellTopologyId{static_cast<uint16_t>(i + 5)});
    }
  }
}

BOOST_AUTO_TEST_CASE(CyclicTopologyIsRejectedExplicitly)
{
  SparseTrackingTopology topology{3};
  BOOST_REQUIRE(topology.addTransition(adjacent(0, 1, TransitionPolicyTag::CylinderCylinder)).isValid());
  BOOST_REQUIRE(topology.addTransition(adjacent(1, 2, TransitionPolicyTag::CylinderCylinder)).isValid());
  BOOST_REQUIRE(topology.addTransition(adjacent(2, 0, TransitionPolicyTag::CylinderCylinder)).isValid());
  BOOST_REQUIRE(topology.finalize());
  DetectorLayout layout{{surface(0, SurfaceKind::Cylinder), surface(1, SurfaceKind::Cylinder), surface(2, SurfaceKind::Cylinder)}, std::move(topology)};
  BOOST_REQUIRE(layout.valid());

  TransitionPolicyGrouping grouping{layout.getView()};
  BOOST_CHECK(!grouping.valid());
  BOOST_CHECK(grouping.getScheduleError() == TransitionPolicyScheduleError::CyclicTopology);
  BOOST_CHECK(grouping.transitionsForTag(TransitionPolicyTag::CylinderCylinder).empty());
  BOOST_CHECK(grouping.scheduledCellsForTag(TransitionPolicyTag::CylinderCylinder).empty());
}

BOOST_AUTO_TEST_CASE(ItsLikeLayoutOnlyGroupsCylinderCylinderWork)
{
  const auto layout = buildChainLayout(7, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(layout.valid());
  TransitionPolicyGrouping grouping{layout.getView()};

  BOOST_CHECK(grouping.hasTag(TransitionPolicyTag::CylinderCylinder));
  BOOST_CHECK(!grouping.hasTag(TransitionPolicyTag::DiskDisk));
  BOOST_CHECK_EQUAL(grouping.transitionsForTag(TransitionPolicyTag::CylinderCylinder).size(), 6u);
  BOOST_CHECK_EQUAL(grouping.cellsForTag(TransitionPolicyTag::CylinderCylinder).size(), 5u);

  RecordingVisitor visitor;
  dispatchTransitionPolicies(grouping, visitor);
  BOOST_CHECK_EQUAL(visitor.calls, 1);
  BOOST_CHECK(visitor.families.at(TransitionPolicyTag::CylinderCylinder) == StateFamily::Barrel);
  BOOST_CHECK_EQUAL(visitor.transitionCounts.at(TransitionPolicyTag::CylinderCylinder), 6u);
  BOOST_CHECK_EQUAL(visitor.cellCounts.at(TransitionPolicyTag::CylinderCylinder), 5u);
  BOOST_CHECK(visitor.families.find(TransitionPolicyTag::DiskDisk) == visitor.families.end());
}

BOOST_AUTO_TEST_CASE(MftLikeLayoutOnlyGroupsDiskDiskWork)
{
  const auto layout = buildChainLayout(10, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk);
  BOOST_REQUIRE(layout.valid());
  TransitionPolicyGrouping grouping{layout.getView()};

  BOOST_CHECK(!grouping.hasTag(TransitionPolicyTag::CylinderCylinder));
  BOOST_CHECK(grouping.hasTag(TransitionPolicyTag::DiskDisk));
  BOOST_CHECK_EQUAL(grouping.transitionsForTag(TransitionPolicyTag::DiskDisk).size(), 9u);
  BOOST_CHECK_EQUAL(grouping.cellsForTag(TransitionPolicyTag::DiskDisk).size(), 8u);

  RecordingVisitor visitor;
  dispatchTransitionPolicies(grouping, visitor);
  BOOST_CHECK_EQUAL(visitor.calls, 1);
  BOOST_CHECK(visitor.families.at(TransitionPolicyTag::DiskDisk) == StateFamily::Forward);
  BOOST_CHECK_EQUAL(visitor.transitionCounts.at(TransitionPolicyTag::DiskDisk), 9u);
  BOOST_CHECK_EQUAL(visitor.cellCounts.at(TransitionPolicyTag::DiskDisk), 8u);
  BOOST_CHECK(visitor.families.find(TransitionPolicyTag::CylinderCylinder) == visitor.families.end());
}

BOOST_AUTO_TEST_CASE(CombinedDisconnectedLayoutDispatchesBothFamiliesWithoutDetectorBranching)
{
  // Same chain lengths as the two single-detector cases above, so the
  // combined per-family counts must equal the corresponding standalone
  // layout's counts -- proof that grouping/dispatch depends only on the
  // topology's policy tags, not on which detector produced the layout.
  const auto combined = buildCombinedDisconnectedLayout(7, 10);
  BOOST_REQUIRE(combined.valid());
  TransitionPolicyGrouping grouping{combined.getView()};

  BOOST_CHECK(grouping.hasTag(TransitionPolicyTag::CylinderCylinder));
  BOOST_CHECK(grouping.hasTag(TransitionPolicyTag::DiskDisk));

  RecordingVisitor visitor;
  dispatchTransitionPolicies(grouping, visitor);
  BOOST_CHECK_EQUAL(visitor.calls, 2);
  BOOST_CHECK_EQUAL(visitor.transitionCounts.at(TransitionPolicyTag::CylinderCylinder), 6u);
  BOOST_CHECK_EQUAL(visitor.cellCounts.at(TransitionPolicyTag::CylinderCylinder), 5u);
  BOOST_CHECK_EQUAL(visitor.transitionCounts.at(TransitionPolicyTag::DiskDisk), 9u);
  BOOST_CHECK_EQUAL(visitor.cellCounts.at(TransitionPolicyTag::DiskDisk), 8u);
  BOOST_CHECK(visitor.families.at(TransitionPolicyTag::CylinderCylinder) == StateFamily::Barrel);
  BOOST_CHECK(visitor.families.at(TransitionPolicyTag::DiskDisk) == StateFamily::Forward);
}

BOOST_AUTO_TEST_CASE(GroupedTransitionsOnlyConnectSurfacesOfTheExpectedKind)
{
  // Policy/surface compatibility: every transition the grouping assigns to a
  // tag must connect surfaces of that tag's expected kind. DetectorLayout
  // construction already enforces this (Gate 1); this proves the dispatch
  // boundary's own grouping preserves the guarantee end-to-end.
  const auto combined = buildCombinedDisconnectedLayout(7, 10);
  BOOST_REQUIRE(combined.valid());
  const auto view = combined.getView();
  TransitionPolicyGrouping grouping{view};

  for (const auto tag : {TransitionPolicyTag::CylinderCylinder, TransitionPolicyTag::DiskDisk}) {
    const auto transitionIds = grouping.transitionsForTag(tag);
    BOOST_REQUIRE_GT(transitionIds.size(), 0u);
    for (const auto id : transitionIds) {
      const auto& transition = view.topology.getTransition(id);
      const auto fromKind = view.getSurface(transition.from).kind;
      const auto toKind = view.getSurface(transition.to).kind;
      BOOST_CHECK(isSurfaceKindCompatible(tag, fromKind));
      BOOST_CHECK(isSurfaceKindCompatible(tag, toKind));
    }
  }
}

BOOST_AUTO_TEST_CASE(DispatchIsANoOpForALayoutWithNoActiveWork)
{
  SparseTrackingTopology topology{2};
  BOOST_REQUIRE(topology.finalize());
  DetectorLayout layout{{surface(0, SurfaceKind::Cylinder), surface(1, SurfaceKind::Cylinder)}, std::move(topology)};
  BOOST_REQUIRE(layout.valid());
  TransitionPolicyGrouping grouping{layout.getView()};

  RecordingVisitor visitor;
  dispatchTransitionPolicies(grouping, visitor);
  BOOST_CHECK_EQUAL(visitor.calls, 0);
}
