// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT TransitionPolicyDispatch
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <algorithm>
#include <array>
#include <map>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/detail/TransitionPolicyDispatch.h"

using namespace o2::itsmft::tracking;

namespace
{

SurfaceTransition adjacent(uint16_t from, uint16_t to)
{
  return SurfaceTransition{SurfaceId{from}, SurfaceId{to}, SurfaceMask{}, 0};
}

SurfaceDescriptor surface(uint16_t id, SurfaceKind kind)
{
  return SurfaceDescriptor{SurfaceId{id}, id, 0, kind};
}

SurfaceMask maskOf(uint16_t id)
{
  SurfaceMask mask;
  mask.set(SurfaceId{id});
  return mask;
}

/// SurfaceGraph no longer owns a surface copy (Slice 3, shared ownership):
/// it borrows a caller-supplied catalog span only for construction/view
/// assembly. Test fixtures that build one SurfaceGraph in isolation (no
/// std::vector<SurfaceGraph> involved) keep their own catalog alongside it so
/// `.getView()` keeps working as a zero-argument call at every existing call
/// site below.
struct BuiltLayout {
  SurfaceGraph layout;
  std::vector<SurfaceDescriptor> surfaces;

  bool valid() const noexcept { return layout.valid(); }
  SurfaceGraphView getView() const noexcept
  {
    return layout.getView();
  }
};

/// A chain of `nSurfaces` surfaces of one kind, adjacent-only transitions and
/// cells (no holes) -- ITS-like when Cylinder/CylinderCylinder, MFT-like when
/// Disk/DiskDisk. Mirrors the fixture style already accepted in
/// testSurfaceLayout.cxx. `seedingSurfaces` defaults to empty so existing
/// (pre-road-start) callers are unaffected.
BuiltLayout buildChainLayout(uint16_t nSurfaces, SurfaceKind kind, TransitionPolicyTag tag, SurfaceMask seedingSurfaces = {})
{
  SurfaceGraph topology{nSurfaces, seedingSurfaces};
  std::vector<TransitionId> transitions;
  for (uint16_t s = 0; s + 1 < nSurfaces; ++s) {
    transitions.push_back(topology.addTransition(adjacent(s, s + 1)));
  }
  for (size_t t = 0; t + 1 < transitions.size(); ++t) {
    topology.addCell(transitions[t], transitions[t + 1]);
  }
  BOOST_REQUIRE(topology.finalize());

  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t s = 0; s < nSurfaces; ++s) {
    surfaces.push_back(surface(s, kind));
  }
  return BuiltLayout{SurfaceGraph{surfaces, std::move(topology)}, std::move(surfaces)};
}

/// Two disconnected chains sharing one TransitionPolicyTag: surfaces
/// [0, nA) and [nA, nA+nB), each an adjacent-only chain of `kind`. Used to
/// prove road-start selection needs no per-component special-casing (the
/// grouping's rank computation already treats disconnected components
/// independently -- TransitionPolicyDispatch.h's Kahn's-algorithm loop).
BuiltLayout buildTwoChainsSameTag(uint16_t nA, uint16_t nB, SurfaceKind kind, TransitionPolicyTag tag, SurfaceMask seedingSurfaces = {})
{
  const uint16_t total = static_cast<uint16_t>(nA + nB);
  SurfaceGraph topology{total, seedingSurfaces};
  std::vector<TransitionId> transitionsA, transitionsB;
  for (uint16_t s = 0; s + 1 < nA; ++s) {
    transitionsA.push_back(topology.addTransition(adjacent(s, s + 1)));
  }
  for (uint16_t s = nA; s + 1 < total; ++s) {
    transitionsB.push_back(topology.addTransition(adjacent(s, s + 1)));
  }
  for (size_t t = 0; t + 1 < transitionsA.size(); ++t) {
    topology.addCell(transitionsA[t], transitionsA[t + 1]);
  }
  for (size_t t = 0; t + 1 < transitionsB.size(); ++t) {
    topology.addCell(transitionsB[t], transitionsB[t + 1]);
  }
  BOOST_REQUIRE(topology.finalize());

  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t s = 0; s < total; ++s) {
    surfaces.push_back(surface(s, kind));
  }
  return BuiltLayout{SurfaceGraph{surfaces, std::move(topology)}, std::move(surfaces)};
}

/// A disconnected combined layout: an ITS-like cylinder chain over the first
/// `nCylinders` surfaces and an MFT-like disk chain over the next `nDisks`
/// surfaces, in one shared global surface-id space (Architecture.md §8, the
/// first combined milestone).
BuiltLayout buildCombinedDisconnectedLayout(uint16_t nCylinders, uint16_t nDisks, SurfaceMask seedingSurfaces = {})
{
  const uint16_t total = nCylinders + nDisks;
  SurfaceGraph topology{total, seedingSurfaces};
  std::vector<TransitionId> cylinderTransitions;
  std::vector<TransitionId> diskTransitions;
  for (uint16_t s = 0; s + 1 < nCylinders; ++s) {
    cylinderTransitions.push_back(topology.addTransition(adjacent(s, s + 1)));
  }
  for (uint16_t s = nCylinders; s + 1 < total; ++s) {
    diskTransitions.push_back(topology.addTransition(adjacent(s, s + 1)));
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
  return BuiltLayout{SurfaceGraph{surfaces, std::move(topology)}, std::move(surfaces)};
}

BuiltLayout buildIdentityLayout(uint16_t nSurfaces, SurfaceKind kind, TransitionPolicyTag tag, int maxHoles, uint16_t hole, SurfaceMask seedingSurfaces = {})
{
  std::vector<SurfaceDescriptor> surfaces;
  std::vector<SurfaceId> order;
  for (uint16_t id = 0; id < nSurfaces; ++id) {
    surfaces.push_back(surface(id, kind));
    order.push_back(SurfaceId{id});
  }
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}};
  auto result = builder.addSubgraph(SurfaceGraphSubgraph{std::move(order), maxHoles, maxHoles ? maskOf(hole) : SurfaceMask{}, seedingSurfaces}).build();
  BOOST_REQUIRE(result.ok());
  return BuiltLayout{std::move(*result.graph), std::move(surfaces)};
}

void checkIdentitySchedule(uint16_t nSurfaces, SurfaceKind kind, TransitionPolicyTag tag, int maxHoles, uint16_t hole)
{
  const auto layout = buildIdentityLayout(nSurfaces, kind, tag, maxHoles, hole);
  TransitionPolicyGrouping grouping{layout.getView()};
  BOOST_REQUIRE(grouping.valid());
  const auto scheduled = grouping.scheduledCellsForTag(tag);
  BOOST_REQUIRE_EQUAL(scheduled.size(), layout.layout.getCells().size());
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
  TransitionPolicyGrouping grouping{SurfaceGraphView{}};
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
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}};
  auto result = builder.addSubgraph(SurfaceGraphSubgraph{{SurfaceId{3}, SurfaceId{1}, SurfaceId{4}, SurfaceId{0}}, 0, {}, {}}).build();
  BOOST_REQUIRE(result.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  TransitionPolicyGrouping grouping{result.graph->getView()};
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
  SurfaceGraph topology{3};
  BOOST_REQUIRE(topology.addTransition(adjacent(0, 1)).isValid());
  BOOST_REQUIRE(topology.addTransition(adjacent(1, 2)).isValid());
  BOOST_REQUIRE(topology.addTransition(adjacent(2, 0)).isValid());
  BOOST_REQUIRE(topology.finalize());
  const std::vector<SurfaceDescriptor> surfaces{surface(0, SurfaceKind::Cylinder), surface(1, SurfaceKind::Cylinder), surface(2, SurfaceKind::Cylinder)};
  SurfaceGraph layout{surfaces, std::move(topology)};
  BOOST_REQUIRE(layout.valid());

  const auto masks = computeSurfaceKindMasks(surfaces);
  TransitionPolicyGrouping grouping{layout.getView()};
  BOOST_CHECK(!grouping.valid());
  BOOST_CHECK(grouping.getScheduleError() == TransitionPolicyScheduleError::CyclicTopology);
  BOOST_CHECK(grouping.transitionsForTag(TransitionPolicyTag::CylinderCylinder).empty());
  BOOST_CHECK(grouping.scheduledCellsForTag(TransitionPolicyTag::CylinderCylinder).empty());
  BOOST_CHECK(grouping.cellsForTag(TransitionPolicyTag::CylinderCylinder).empty());
  BOOST_CHECK(grouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder).empty());
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
  // tag must connect surfaces of that tag's expected kind. SurfaceGraph
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
      const auto& transition = view.getTransition(id);
      const auto fromKind = view.getSurface(transition.from).kind;
      const auto toKind = view.getSurface(transition.to).kind;
      BOOST_CHECK(isSurfaceKindCompatible(tag, fromKind));
      BOOST_CHECK(isSurfaceKindCompatible(tag, toKind));
    }
  }
}

BOOST_AUTO_TEST_CASE(PolicyTagDerivationMatchesEndpointSurfaceKindForEveryItsAndMftTransition)
{
  // M4 (GenericTrackingEngineMigration.md; ADR 0007 decision 8): SurfaceTransition
  // no longer stores a policy tag -- TransitionPolicyGrouping derives it from
  // each transition's endpoint SurfaceDescriptor::kind instead
  // (transitionPolicyTagForSurfaceKind(), detail/TransitionPolicy.h). This
  // proves that derivation is byte-for-byte equivalent to the removed stored
  // field for every transition in an ITS-shaped (7 Cylinder layers) plus
  // MFT-shaped (10 Disk layers) combined layout -- the current production
  // topology shape -- and that the grouping built from it places every
  // transition/cell in exactly the tag bucket the endpoint kind implies.
  constexpr uint16_t kItsLikeLayers = 7;  // ITSNLayers
  constexpr uint16_t kMftLikeLayers = 10; // o2::mft::constants::mft::LayersNumber
  const auto combined = buildCombinedDisconnectedLayout(kItsLikeLayers, kMftLikeLayers);
  BOOST_REQUIRE(combined.valid());
  const auto view = combined.getView();
  TransitionPolicyGrouping grouping{view};
  BOOST_REQUIRE(grouping.valid());

  size_t checkedTransitions = 0;
  for (uint32_t id = 0; id < view.nTransitions; ++id) {
    const auto transitionId = TransitionId{static_cast<uint16_t>(id)};
    const auto& transition = view.getTransition(transitionId);
    const auto fromTag = transitionPolicyTagForSurfaceKind(view.getSurface(transition.from).kind);
    const auto toTag = transitionPolicyTagForSurfaceKind(view.getSurface(transition.to).kind);
    // Same-family-only topology (ADR 0007 decision 8): both endpoints must
    // derive the same tag -- SurfaceGraph::validate()'s MixedSurfaceTransition
    // check already guarantees fromKind == toKind for every valid layout.
    BOOST_CHECK(fromTag == toTag);
    BOOST_CHECK(fromTag != TransitionPolicyTag::Invalid);
    const auto ownGroup = grouping.transitionsForTag(fromTag);
    BOOST_CHECK(std::find(ownGroup.begin(), ownGroup.end(), transitionId) != ownGroup.end());
    ++checkedTransitions;
  }
  BOOST_CHECK_EQUAL(checkedTransitions, view.nTransitions);
  BOOST_CHECK_EQUAL(grouping.transitionsForTag(TransitionPolicyTag::CylinderCylinder).size(), kItsLikeLayers - 1u);
  BOOST_CHECK_EQUAL(grouping.transitionsForTag(TransitionPolicyTag::DiskDisk).size(), kMftLikeLayers - 1u);
}

BOOST_AUTO_TEST_CASE(DispatchIsANoOpForALayoutWithNoActiveWork)
{
  SurfaceGraph topology{2};
  BOOST_REQUIRE(topology.finalize());
  const std::vector<SurfaceDescriptor> surfaces{surface(0, SurfaceKind::Cylinder), surface(1, SurfaceKind::Cylinder)};
  SurfaceGraph layout{surfaces, std::move(topology)};
  BOOST_REQUIRE(layout.valid());
  const auto masks = computeSurfaceKindMasks(surfaces);
  TransitionPolicyGrouping grouping{layout.getView()};

  RecordingVisitor visitor;
  dispatchTransitionPolicies(grouping, visitor);
  BOOST_CHECK_EQUAL(visitor.calls, 0);
}

// --- Gate 3 road-start selection (Architecture.md Sec 10, D003) ---------

BOOST_AUTO_TEST_CASE(RoadStartCellsMatchIdentityItsLikeSeeding)
{
  // Full seeding mask (all 7 layers, the ITS StartLayerMask default in
  // Configuration.h): every cell in a monotonic identity chain is a road
  // start, i.e. roadStartCellsForTag must equal cellsForTag exactly.
  const auto full = buildChainLayout(7, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder, SurfaceMask{0x7Fu});
  TransitionPolicyGrouping fullGrouping{full.getView()};
  BOOST_REQUIRE(fullGrouping.valid());
  const auto allCells = fullGrouping.cellsForTag(TransitionPolicyTag::CylinderCylinder);
  const auto allStarts = fullGrouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE_EQUAL(allStarts.size(), allCells.size());
  for (size_t i = 0; i < allCells.size(); ++i) {
    BOOST_CHECK(allStarts[i] == allCells[i]);
  }

  // Restricted mask: only the cell(s) whose transition endpoint is surface 6
  // (the outermost layer) qualify.
  const auto restricted = buildChainLayout(7, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder, maskOf(6));
  TransitionPolicyGrouping restrictedGrouping{restricted.getView()};
  BOOST_REQUIRE(restrictedGrouping.valid());
  const auto restrictedCells = restrictedGrouping.cellsForTag(TransitionPolicyTag::CylinderCylinder);
  const auto restrictedStarts = restrictedGrouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  const auto view = restricted.getView();
  BOOST_REQUIRE_GT(restrictedStarts.size(), 0u);
  BOOST_CHECK(std::is_sorted(restrictedStarts.begin(), restrictedStarts.end()));
  for (const auto id : restrictedCells) {
    const auto& cell = view.getCell(id);
    const bool isEndpoint6 = view.getTransition(cell.secondTransition).to == SurfaceId{6};
    const bool selected = std::find(restrictedStarts.begin(), restrictedStarts.end(), id) != restrictedStarts.end();
    BOOST_CHECK_EQUAL(isEndpoint6, selected);
  }
}

BOOST_AUTO_TEST_CASE(RoadStartCellsMatchIdentityMftLikeSeeding)
{
  // Same proof as the ITS-like case above, over a 10-disk identity chain
  // (MFT's StartLayerMask default is also "all layers set").
  const auto full = buildChainLayout(10, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk, SurfaceMask{0x3FFu});
  TransitionPolicyGrouping fullGrouping{full.getView()};
  BOOST_REQUIRE(fullGrouping.valid());
  const auto allCells = fullGrouping.cellsForTag(TransitionPolicyTag::DiskDisk);
  const auto allStarts = fullGrouping.roadStartCellsForTag(TransitionPolicyTag::DiskDisk);
  BOOST_REQUIRE_EQUAL(allStarts.size(), allCells.size());
  for (size_t i = 0; i < allCells.size(); ++i) {
    BOOST_CHECK(allStarts[i] == allCells[i]);
  }

  const auto restricted = buildChainLayout(10, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk, maskOf(9));
  TransitionPolicyGrouping restrictedGrouping{restricted.getView()};
  BOOST_REQUIRE(restrictedGrouping.valid());
  const auto restrictedCells = restrictedGrouping.cellsForTag(TransitionPolicyTag::DiskDisk);
  const auto restrictedStarts = restrictedGrouping.roadStartCellsForTag(TransitionPolicyTag::DiskDisk);
  const auto view = restricted.getView();
  BOOST_REQUIRE_GT(restrictedStarts.size(), 0u);
  BOOST_CHECK(std::is_sorted(restrictedStarts.begin(), restrictedStarts.end()));
  for (const auto id : restrictedCells) {
    const auto& cell = view.getCell(id);
    const bool isEndpoint9 = view.getTransition(cell.secondTransition).to == SurfaceId{9};
    const bool selected = std::find(restrictedStarts.begin(), restrictedStarts.end(), id) != restrictedStarts.end();
    BOOST_CHECK_EQUAL(isEndpoint9, selected);
  }
}

BOOST_AUTO_TEST_CASE(RoadStartEndpointWinsOverNumericHighestHitSurfaceBit)
{
  // Non-monotonic SurfaceId order: order = {3, 1, 4, 0} at positions 0..3,
  // matching NonMonotonicSurfaceIdsFollowGraphRank above. Cells:
  //   cell0: transitions (3->1),(1->4), hitSurfaces {3,1,4}, endpoint = 4
  //   cell1: transitions (1->4),(4->0), hitSurfaces {1,4,0}, endpoint = 0
  // For cell1, the numerically greatest set bit of hitSurfaces is 4 (not 0):
  // a LayerMask::last()-style oracle would misidentify the endpoint. Seeding
  // surface 4 makes this concrete: the correct (transition-endpoint) answer
  // selects only cell0; a numeric-highest-bit reading of hitSurfaces would
  // wrongly also select cell1 (whose hitSurfaces contains 4).
  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t id = 0; id < 5; ++id) {
    surfaces.push_back(surface(id, SurfaceKind::Cylinder));
  }
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}};
  auto result = builder.addSubgraph(SurfaceGraphSubgraph{{SurfaceId{3}, SurfaceId{1}, SurfaceId{4}, SurfaceId{0}}, 0, {}, maskOf(4)}).build();
  BOOST_REQUIRE(result.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  TransitionPolicyGrouping grouping{result.graph->getView()};
  BOOST_REQUIRE(grouping.valid());
  const auto starts = grouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE_EQUAL(starts.size(), 1u);
  BOOST_CHECK(starts[0] == CellTopologyId{0});

  // Seeding surface 0 (cell1's true, transition-graph endpoint) selects only
  // cell1, even though surface 0 is not the numerically greatest bit of
  // either cell's hitSurfaces.
  std::vector<SurfaceDescriptor> surfaces2;
  for (uint16_t id = 0; id < 5; ++id) {
    surfaces2.push_back(surface(id, SurfaceKind::Cylinder));
  }
  SurfaceGraphBuilder builder2{SurfaceCatalogView{surfaces2.data(), static_cast<uint32_t>(surfaces2.size())}};
  auto result2 = builder2.addSubgraph(SurfaceGraphSubgraph{{SurfaceId{3}, SurfaceId{1}, SurfaceId{4}, SurfaceId{0}}, 0, {}, maskOf(0)}).build();
  BOOST_REQUIRE(result2.ok());
  const auto masks2 = computeSurfaceKindMasks(surfaces2);
  TransitionPolicyGrouping grouping2{result2.graph->getView()};
  BOOST_REQUIRE(grouping2.valid());
  const auto starts2 = grouping2.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE_EQUAL(starts2.size(), 1u);
  BOOST_CHECK(starts2[0] == CellTopologyId{1});
}

BOOST_AUTO_TEST_CASE(RoadStartCellsSeparateDisconnectedSamePolicyComponents)
{
  // Two disconnected CylinderCylinder chains sharing surfaces [0,4) and
  // [4,8). Seeding surfaces 3 (component A's endpoint) and 7 (component B's
  // endpoint) must each select exactly one cell, from the correct component,
  // in ascending CellTopologyId order (component A's cells are constructed
  // before component B's).
  const auto layout = buildTwoChainsSameTag(4, 4, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder, maskOf(3) | maskOf(7));
  TransitionPolicyGrouping grouping{layout.getView()};
  BOOST_REQUIRE(grouping.valid());
  const auto cells = grouping.cellsForTag(TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE_EQUAL(cells.size(), 4u); // 2 cells per 4-surface chain
  const auto starts = grouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE_EQUAL(starts.size(), 2u);
  BOOST_CHECK(std::is_sorted(starts.begin(), starts.end()));

  const auto view = layout.getView();
  for (const auto id : starts) {
    const auto& cell = view.getCell(id);
    const auto endpoint = view.getTransition(cell.secondTransition).to;
    BOOST_CHECK(endpoint == SurfaceId{3} || endpoint == SurfaceId{7});
  }
  // Ascending order means the component-A cell (lower CellTopologyId, built
  // first) precedes the component-B cell.
  BOOST_CHECK(starts[0].value() < starts[1].value());
}

BOOST_AUTO_TEST_CASE(RoadStartCellsHandleMultipleSeedingSurfaces)
{
  const auto layout = buildChainLayout(7, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder, maskOf(4) | maskOf(6));
  TransitionPolicyGrouping grouping{layout.getView()};
  BOOST_REQUIRE(grouping.valid());
  const auto starts = grouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE_EQUAL(starts.size(), 2u);
  BOOST_CHECK(std::is_sorted(starts.begin(), starts.end()));
  const auto view = layout.getView();
  for (const auto id : starts) {
    const auto& cell = view.getCell(id);
    const auto endpoint = view.getTransition(cell.secondTransition).to;
    BOOST_CHECK(endpoint == SurfaceId{4} || endpoint == SurfaceId{6});
  }
}

BOOST_AUTO_TEST_CASE(RoadStartCellsAreEmptyForAnEmptySeedingMask)
{
  const auto layout = buildChainLayout(7, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder, SurfaceMask{});
  TransitionPolicyGrouping grouping{layout.getView()};
  BOOST_REQUIRE(grouping.valid());
  BOOST_CHECK(!grouping.cellsForTag(TransitionPolicyTag::CylinderCylinder).empty());
  BOOST_CHECK(grouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder).empty());
}

BOOST_AUTO_TEST_CASE(RoadStartCellsAreEmptyWhenSeedingSurfaceTerminatesNoCell)
{
  // Surface 0 (the very first position of a monotonic chain) is never the
  // `to` endpoint of any transition, so it can never be a cell's
  // secondTransition.to -- a seeding surface that is valid but unreachable
  // as an endpoint.
  const auto layout = buildChainLayout(7, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder, maskOf(0));
  TransitionPolicyGrouping grouping{layout.getView()};
  BOOST_REQUIRE(grouping.valid());
  BOOST_CHECK(!grouping.cellsForTag(TransitionPolicyTag::CylinderCylinder).empty());
  BOOST_CHECK(grouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder).empty());
}

BOOST_AUTO_TEST_CASE(RoadStartCellsToleratesSkippedNonAdjacentTransitions)
{
  // maxHoles=1 with hole surface 3 permits the skip transition 2->4
  // (skipping 3) alongside the adjacent ones. Seeding surface 5 is reachable
  // both via the plain adjacent chain (3,4,5) and via the hole-skipping
  // chain (2,4,5) [first transition 2->4 skips 3, second transition 4->5].
  // The endpoint definition (secondTransition.to) must keep selecting
  // exactly the cells ending at the seeded surface regardless of whether
  // their first transition is adjacent or hole-skipping.
  const auto layout = buildIdentityLayout(7, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder, 1, 3, maskOf(5));
  TransitionPolicyGrouping grouping{layout.getView()};
  BOOST_REQUIRE(grouping.valid());
  const auto starts = grouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE_GT(starts.size(), 0u);
  BOOST_CHECK(std::is_sorted(starts.begin(), starts.end()));
  const auto view = layout.getView();
  bool sawSkippedTransition = false;
  for (const auto id : starts) {
    const auto& cell = view.getCell(id);
    BOOST_CHECK(view.getTransition(cell.secondTransition).to == SurfaceId{5});
    if (!view.getTransition(cell.firstTransition).skippedSurfaces.empty() ||
        !view.getTransition(cell.secondTransition).skippedSurfaces.empty()) {
      sawSkippedTransition = true;
    }
  }
  BOOST_CHECK(sawSkippedTransition);
}

BOOST_AUTO_TEST_CASE(RoadStartCellsSeparateCylinderAndDiskSpansInACombinedGrouping)
{
  // Grouping-only combined fixture (mixed-family production activation is
  // rejected elsewhere, in TrackerTraits::initialiseTimeFrame -- not here):
  // seeding surface 6 is the cylinder chain's endpoint, surface 16 is the
  // disk chain's endpoint (surfaces 7..16 for 10 disks). Each tag's
  // roadStartCellsForTag must contain only cells from its own family.
  const auto layout = buildCombinedDisconnectedLayout(7, 10, maskOf(6) | maskOf(16));
  TransitionPolicyGrouping grouping{layout.getView()};
  BOOST_REQUIRE(grouping.valid());

  const auto cylinderStarts = grouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  const auto diskStarts = grouping.roadStartCellsForTag(TransitionPolicyTag::DiskDisk);
  BOOST_REQUIRE_EQUAL(cylinderStarts.size(), 1u);
  BOOST_REQUIRE_EQUAL(diskStarts.size(), 1u);

  const auto view = layout.getView();
  const auto& cylinderCell = view.getCell(cylinderStarts[0]);
  const auto& diskCell = view.getCell(diskStarts[0]);
  BOOST_CHECK(view.getTransition(cylinderCell.secondTransition).to == SurfaceId{6});
  BOOST_CHECK(view.getTransition(diskCell.secondTransition).to == SurfaceId{16});
  // Cross-check per-policy ownership matches the existing cellsForTag split.
  const auto cylinderCells = grouping.cellsForTag(TransitionPolicyTag::CylinderCylinder);
  const auto diskCells = grouping.cellsForTag(TransitionPolicyTag::DiskDisk);
  BOOST_CHECK(std::find(cylinderCells.begin(), cylinderCells.end(), cylinderStarts[0]) != cylinderCells.end());
  BOOST_CHECK(std::find(diskCells.begin(), diskCells.end(), diskStarts[0]) != diskCells.end());
}

BOOST_AUTO_TEST_CASE(ConstructorFailureClearsRoadStartCellsAlongsideAllGroups)
{
  // A hand-built raw view (SurfaceGraphView/SurfaceGraphView are
  // trivially-copyable PODs, already used directly by
  // GroupingIsEmptyForAnEmptyLayoutView above) whose first cell is valid and
  // road-start-eligible, and whose second cell references an out-of-range
  // TransitionId. The grouping must reject the whole schedule
  // (InvalidCellTransition) and clear() must wipe every group's
  // roadStartCells -- not just the ones the offending cell would have
  // touched -- exactly as it already does for transitions/cells/scheduledCells.
  std::array<SurfaceDescriptor, 3> surfaces{
    surface(0, SurfaceKind::Cylinder), surface(1, SurfaceKind::Cylinder), surface(2, SurfaceKind::Cylinder)};
  std::array<SurfaceTransition, 2> transitions{
    adjacent(0, 1),
    adjacent(1, 2)};
  std::array<SurfaceCellTopology, 2> cells{
    SurfaceCellTopology{TransitionId{0}, TransitionId{1}, SurfaceMask{}},   // valid: endpoint surface 2, seeded
    SurfaceCellTopology{TransitionId{0}, TransitionId{99}, SurfaceMask{}}}; // invalid: out-of-range secondTransition
  std::array<uint32_t, 3> offsets{0, 2, 2};
  std::array<CellTopologyId, 2> byFirstTransition{CellTopologyId{0}, CellTopologyId{1}};

  SurfaceMask cylinderSurfaces;
  for (const auto& descriptor : surfaces) {
    cylinderSurfaces.set(descriptor.id);
  }
  SurfaceGraphView layoutView{surfaces.data(), static_cast<uint32_t>(surfaces.size()), nullptr, 0, cylinderSurfaces, {}, transitions.data(), cells.data(), offsets.data(), byFirstTransition.data(), maskOf(2), static_cast<uint32_t>(transitions.size()), static_cast<uint32_t>(cells.size())};

  TransitionPolicyGrouping grouping{layoutView};
  BOOST_CHECK(!grouping.valid());
  BOOST_CHECK(grouping.getScheduleError() == TransitionPolicyScheduleError::InvalidCellTransition);
  for (const auto tag : {TransitionPolicyTag::CylinderCylinder, TransitionPolicyTag::DiskDisk}) {
    BOOST_CHECK(grouping.transitionsForTag(tag).empty());
    BOOST_CHECK(grouping.cellsForTag(tag).empty());
    BOOST_CHECK(grouping.scheduledCellsForTag(tag).empty());
    BOOST_CHECK(grouping.roadStartCellsForTag(tag).empty());
  }
}
