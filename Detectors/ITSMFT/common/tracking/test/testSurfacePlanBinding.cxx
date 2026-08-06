// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// M6b (Detectors/ITSMFT/common/tracking/doc/design/0002-m6-generic-workspace-migration.md
// Sec 3.2, 7, 9): focused coverage for the detector-neutral SurfacePlanBinding
// slot and topology contract. Production-participant wiring is covered by the
// migration tests.

#define BOOST_TEST_MODULE ITSMFT SurfacePlanBinding
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"

namespace
{
using namespace o2::itsmft::tracking;

std::vector<SurfaceId> ordered(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

SurfaceMask maskOf(uint16_t id)
{
  SurfaceMask mask;
  mask.set(SurfaceId{id});
  return mask;
}

SurfaceDescriptor surfaceWithOwner(uint16_t id, SurfaceKind kind, uint8_t detectorId)
{
  return SurfaceDescriptor{SurfaceId{id}, id, detectorId, kind};
}

// --- Shared ITS+MFT combined-catalog fixture. Per-file-local fixtures are
// this test directory's existing convention. ---
struct CombinedLayout {
  SurfaceGraph layout;
  SurfaceGraphView view;

  CombinedLayout()
    : layout{build()}, view{layout.getView()}
  {
  }

 private:
  static std::pair<SurfaceMask, SurfaceMask> masks()
  {
    return computeSurfaceKindMasks(kITSMFTCombinedStaticSurfaceCatalog);
  }
  static SurfaceMask cylinderMask() { return masks().first; }
  static SurfaceMask diskMask() { return masks().second; }
  static SurfaceGraph build()
  {
    SurfaceMask itsHoles;
    itsHoles.set(SurfaceId{3});
    SurfaceMask mftHoles;
    mftHoles.set(SurfaceId{12});
    SurfaceMask seeds;
    seeds.set(SurfaceId{6});
    seeds.set(SurfaceId{16});
    const SurfaceCatalogView catalog{kITSMFTCombinedStaticSurfaceCatalog.data(), static_cast<uint32_t>(kITSMFTCombinedStaticSurfaceCatalog.size())};
    SurfaceGraphBuilder builder{catalog};
    builder.addSubgraph({ordered(0, ITSNLayers), 1, itsHoles, seeds & SurfaceMask{uint32_t{0x7f}}});
    builder.addSubgraph({ordered(ITSNLayers, MFTNLayers), 1, mftHoles, seeds & ~SurfaceMask{uint32_t{0x7f}}});
    auto built = builder.build();
    BOOST_REQUIRE(built.ok());
    return std::move(*built.graph);
  }
};

SurfaceMask itsMask() { return SurfaceMask{uint32_t{0x7f}}; }
SurfaceMask mftMask() { return SurfaceMask{uint32_t{0x1ff80}}; }

void checkBindingCoversOwnedTopology(const SurfacePlanBinding& binding, const SurfaceGraphView& global)
{
  BOOST_CHECK(binding.getSource().isValid());
  for (uint16_t s = 0; s < global.nSurfaces; ++s) {
    if (binding.getOwnedSurfaces().has(SurfaceId{s})) {
      BOOST_REQUIRE(binding.getOwnedSurfaceIndex(SurfaceId{s}));
    }
  }
  for (uint32_t t = 0; t < global.nTransitions; ++t) {
    const auto id = TransitionId{static_cast<uint16_t>(t)};
    const auto& transition = global.getTransition(id);
    if (binding.getOwnedSurfaces().has(transition.from)) {
      BOOST_REQUIRE(binding.getOwnedSurfaces().has(transition.to));
      BOOST_REQUIRE(binding.getScratchTransitionSlot(id));
    }
  }
  for (uint32_t c = 0; c < global.nCells; ++c) {
    const auto id = CellTopologyId{static_cast<uint16_t>(c)};
    const auto& cell = global.getCell(id);
    const bool ownedTransition = binding.getScratchTransitionSlot(cell.firstTransition).has_value() ||
                                 binding.getScratchTransitionSlot(cell.secondTransition).has_value();
    if (ownedTransition) {
      BOOST_REQUIRE(binding.getScratchCellSlot(id));
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(SurfacePlanBindingMapsCombinedItsAndMftPlans)
{
  CombinedLayout combined;

  const auto its = SurfacePlanBinding::build(combined.view, ClusterSourceId{0}, itsMask(), ordered(0, ITSNLayers),
                                             SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  const auto mft = SurfacePlanBinding::build(combined.view, ClusterSourceId{1}, mftMask(), ordered(ITSNLayers, MFTNLayers),
                                             SurfaceKind::Disk, TransitionPolicyTag::DiskDisk);
  BOOST_REQUIRE(its.ok());
  BOOST_REQUIRE(mft.ok());

  BOOST_CHECK(its.binding->getOwnedSurfaces() == itsMask());
  BOOST_CHECK(mft.binding->getOwnedSurfaces() == mftMask());
  checkBindingCoversOwnedTopology(*its.binding, combined.view);
  checkBindingCoversOwnedTopology(*mft.binding, combined.view);
}

BOOST_AUTO_TEST_CASE(SurfacePlanBindingBuildsForASyntheticNonItsMftDetector)
{
  // Requirement 2b: a synthetic detector this library has never heard of
  // (detectorId 250, no corresponding o2::detectors::DetID::ID value at all)
  // must build successfully given internally consistent expected kind/policy.
  constexpr uint8_t kSyntheticDetectorId = 250;
  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t id = 0; id < 4; ++id) {
    surfaces.push_back(surfaceWithOwner(id, SurfaceKind::Cylinder, kSyntheticDetectorId));
  }
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}};
  auto built = builder.addSubgraph({ordered(0, 4), 0, SurfaceMask{}, maskOf(3)}).build();
  BOOST_REQUIRE(built.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = built.graph->getView();

  SurfaceMask owned;
  for (uint16_t id = 0; id < 4; ++id) {
    owned.set(SurfaceId{id});
  }
  const auto result = SurfacePlanBinding::build(view, ClusterSourceId{7}, owned, ordered(0, 4),
                                                SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result.binding->getGlobalTransitions().size(), 3u);
  BOOST_CHECK_EQUAL(result.binding->getGlobalCells().size(), 2u);
  for (uint16_t id = 0; id < 4; ++id) {
    BOOST_REQUIRE(result.binding->getOwnedSurfaceIndex(SurfaceId{id}));
    BOOST_CHECK_EQUAL(*result.binding->getOwnedSurfaceIndex(SurfaceId{id}), id);
  }
}

BOOST_AUTO_TEST_CASE(SparsePlanPositionsAreTheOnlyRuntimeCountAndOrderAuthority)
{
  // The catalog remains dense because SurfaceGraph deliberately uses dense
  // global ids, while the application plan owns a sparse, non-identity order.
  // This is the shape that catches accidental numeric-SurfaceId traversal.
  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t id = 0; id < 8; ++id) {
    surfaces.push_back(surfaceWithOwner(id, SurfaceKind::Cylinder, 250));
  }
  const std::vector<SurfaceId> planOrder{SurfaceId{5}, SurfaceId{2}, SurfaceId{7}, SurfaceId{1}};
  SurfaceMask owned;
  for (const auto id : planOrder) {
    owned.set(id);
  }
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}};
  auto built = builder.addSubgraph({planOrder, 0, SurfaceMask{}, maskOf(5)}).build();
  BOOST_REQUIRE(built.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = built.graph->getView();
  const auto bindingResult = SurfacePlanBinding::build(view, ClusterSourceId{7}, owned, planOrder,
                                                       SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(bindingResult.ok());
  const auto& binding = *bindingResult.binding;

  BOOST_CHECK(std::equal(binding.getOrderedSurfaces().begin(), binding.getOrderedSurfaces().end(), planOrder.begin(), planOrder.end()));
  for (std::size_t position = 0; position < planOrder.size(); ++position) {
    BOOST_REQUIRE(binding.getOwnedSurfaceIndex(planOrder[position]));
    BOOST_CHECK_EQUAL(*binding.getOwnedSurfaceIndex(planOrder[position]), position);
  }

  SurfaceTrackingScratch scratch;
  scratch.adoptPlan(binding.getOrderedSurfaces().size(), binding.getGlobalTransitions().size(), binding.getGlobalCells().size());
  BOOST_CHECK_EQUAL(scratch.getNOwnedSurfaces(), planOrder.size());
  BOOST_CHECK_EQUAL(scratch.getNTransitions(), binding.getGlobalTransitions().size());
  BOOST_CHECK_EQUAL(scratch.getNCells(), binding.getGlobalCells().size());
  for (std::size_t slot = 0; slot < binding.getGlobalTransitions().size(); ++slot) {
    BOOST_REQUIRE(binding.getScratchTransitionSlot(binding.getGlobalTransitions()[slot]));
    BOOST_CHECK_EQUAL(*binding.getScratchTransitionSlot(binding.getGlobalTransitions()[slot]), slot);
  }
  for (std::size_t slot = 0; slot < binding.getGlobalCells().size(); ++slot) {
    BOOST_REQUIRE(binding.getScratchCellSlot(binding.getGlobalCells()[slot]));
    BOOST_CHECK_EQUAL(*binding.getScratchCellSlot(binding.getGlobalCells()[slot]), slot);
  }

  std::vector<std::vector<SurfaceMeasurement>> measurements(planOrder.size());
  std::vector<gsl::span<const SurfaceMeasurement>> measurementViews;
  measurementViews.reserve(planOrder.size());
  for (std::size_t position = 0; position < planOrder.size(); ++position) {
    SurfaceMeasurement measurement;
    measurement.surface = binding.getOrderedSurfaces()[position];
    measurement.cluster = ClusterRef{binding.getSource(), static_cast<uint32_t>(100 + position)};
    measurements[position].push_back(measurement);
    measurementViews.emplace_back(measurements[position]);
  }
  for (std::size_t position = 0; position < measurementViews.size(); ++position) {
    BOOST_CHECK_EQUAL(measurementViews[position][0].surface.value(), planOrder[position].value());
    BOOST_CHECK(measurementViews[position][0].cluster.source == binding.getSource());
  }

  TrackSeed seed;
  SurfaceMask activePositions;
  for (std::size_t position = 0; position < planOrder.size(); ++position) {
    seed.setCluster(static_cast<int>(position), static_cast<int>(100 + position));
    activePositions.set(SurfaceId{static_cast<uint16_t>(position)});
  }
  seed.setSurfaceMask(activePositions);
  BOOST_CHECK_EQUAL(seed.getActiveSurfaceCount(), planOrder.size());
  BOOST_CHECK(seed.getSurfaceMask() == SurfaceMask{uint32_t{0x0f}});
  for (std::size_t position = 0; position < planOrder.size(); ++position) {
    BOOST_CHECK(seed.hasCluster(static_cast<int>(position)));
    BOOST_CHECK_EQUAL(seed.getCluster(static_cast<int>(position)), 100 + position);
  }
  // The first plan position is global SurfaceId 5, but the TrackSeed bit is
  // position 0.  This assertion makes a global-id-as-layer-index regression
  // fail even though all four positions are active.
  BOOST_CHECK(!seed.getSurfaceMask().has(planOrder.front()));
}

BOOST_AUTO_TEST_CASE(RejectsInvalidSource)
{
  CombinedLayout combined;
  const auto result = SurfacePlanBinding::build(combined.view, ClusterSourceId{}, itsMask(), ordered(0, ITSNLayers),
                                                SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidSource);
}

BOOST_AUTO_TEST_CASE(RejectsIncompatibleExpectedPolicyKind)
{
  CombinedLayout combined;
  const auto result = SurfacePlanBinding::build(combined.view, ClusterSourceId{0}, itsMask(), ordered(0, ITSNLayers),
                                                SurfaceKind::Disk, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::IncompatibleExpectedPolicyKind);
}

BOOST_AUTO_TEST_CASE(RejectsSurfaceMaskSizeMismatch)
{
  CombinedLayout combined;
  const auto result = SurfacePlanBinding::build(combined.view, ClusterSourceId{0}, itsMask(), ordered(0, ITSNLayers - 1),
                                                SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidSurfaceMask);
}

BOOST_AUTO_TEST_CASE(RejectsSurfaceMaskNotASubsetOfLayout)
{
  std::vector<SurfaceDescriptor> surfaces{surfaceWithOwner(0, SurfaceKind::Cylinder, 250), surfaceWithOwner(1, SurfaceKind::Cylinder, 250)};
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}};
  auto built = builder.addSubgraph({ordered(0, 2), 0, SurfaceMask{}, SurfaceMask{}}).build();
  BOOST_REQUIRE(built.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = built.graph->getView();

  SurfaceMask outOfRange = maskOf(0) | maskOf(1) | maskOf(5); // surface 5 does not exist in this 2-surface layout
  const std::vector<SurfaceId> order{SurfaceId{0}, SurfaceId{1}, SurfaceId{5}};
  const auto result = SurfacePlanBinding::build(view, ClusterSourceId{0}, outOfRange, order,
                                                SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidSurfaceMask);
}

BOOST_AUTO_TEST_CASE(RejectsDuplicateSurfaceInOrder)
{
  CombinedLayout combined;
  auto order = ordered(0, ITSNLayers);
  order[1] = order[0]; // duplicate
  const auto result = SurfacePlanBinding::build(combined.view, ClusterSourceId{0}, itsMask(), order,
                                                SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidLegacySurfaceOrder);
}

BOOST_AUTO_TEST_CASE(RejectsUnownedSurfaceInOrder)
{
  CombinedLayout combined;
  SurfaceMask owned = maskOf(0);
  const std::vector<SurfaceId> order{SurfaceId{1}}; // count matches (1) but surface 1 is not owned
  const auto result = SurfacePlanBinding::build(combined.view, ClusterSourceId{0}, owned, order,
                                                SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidLegacySurfaceOrder);
}

BOOST_AUTO_TEST_CASE(SurfacePlanBindingBuildsAcrossMultipleDistinctDetectorIdentitiesInOneOwnedSet)
{
  // SurfacePlanBinding must be generic over its own owned SurfaceId set: it
  // must not assume "one binding, one detector". Three compatible (Cylinder)
  // surfaces spanning two distinct synthetic detectorIds (250, 251), one
  // valid source, consistent expected kind/policy, and an internally valid
  // chain topology -- this must build successfully, unlike
  // The former binding's own single-detector-scoped contract.
  std::vector<SurfaceDescriptor> surfaces{
    surfaceWithOwner(0, SurfaceKind::Cylinder, 250),
    surfaceWithOwner(1, SurfaceKind::Cylinder, 250),
    surfaceWithOwner(2, SurfaceKind::Cylinder, 251)}; // second detectorId, still owned by the same binding
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}};
  auto built = builder.addSubgraph({ordered(0, 3), 0, SurfaceMask{}, maskOf(2)}).build();
  BOOST_REQUIRE(built.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = built.graph->getView();

  SurfaceMask owned = maskOf(0) | maskOf(1) | maskOf(2);
  const auto result = SurfacePlanBinding::build(view, ClusterSourceId{0}, owned, ordered(0, 3),
                                                SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result.binding->getGlobalTransitions().size(), 2u);
  BOOST_CHECK_EQUAL(result.binding->getGlobalCells().size(), 1u);
  for (uint16_t id = 0; id < 3; ++id) {
    BOOST_REQUIRE(result.binding->getOwnedSurfaceIndex(SurfaceId{id}));
    BOOST_CHECK_EQUAL(*result.binding->getOwnedSurfaceIndex(SurfaceId{id}), id);
  }
}

BOOST_AUTO_TEST_CASE(RejectsPolicySurfaceKindMismatch)
{
  std::vector<SurfaceDescriptor> surfaces{surfaceWithOwner(0, SurfaceKind::Disk, 250), surfaceWithOwner(1, SurfaceKind::Disk, 250)};
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}};
  auto built = builder.addSubgraph({ordered(0, 2), 0, SurfaceMask{}, SurfaceMask{}}).build();
  BOOST_REQUIRE(built.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = built.graph->getView();

  SurfaceMask owned = maskOf(0) | maskOf(1);
  // expectedKind/expectedPolicy are self-consistent (Cylinder/CylinderCylinder)
  // but disagree with the actual Disk-kind surfaces in `view`.
  const auto result = SurfacePlanBinding::build(view, ClusterSourceId{0}, owned, ordered(0, 2),
                                                SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidPolicySurface);
}

BOOST_AUTO_TEST_CASE(RejectsCrossBoundaryTransition)
{
  std::vector<SurfaceDescriptor> surfaces{surfaceWithOwner(0, SurfaceKind::Cylinder, 250), surfaceWithOwner(1, SurfaceKind::Cylinder, 250)};
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}};
  auto built = builder.addSubgraph({ordered(0, 2), 0, SurfaceMask{}, SurfaceMask{}}).build();
  BOOST_REQUIRE(built.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = built.graph->getView();

  // Own only surface 0: the 0->1 transition has fromOwned=true, toOwned=false.
  SurfaceMask owned = maskOf(0);
  const std::vector<SurfaceId> order{SurfaceId{0}};
  const auto result = SurfacePlanBinding::build(view, ClusterSourceId{0}, owned, order,
                                                SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::CrossBoundaryTransition);
}

BOOST_AUTO_TEST_CASE(RejectsCrossBoundaryCell)
{
  // Hand-built raw view (SurfaceGraphView/SurfaceGraphView are
  // trivially-copyable PODs -- same technique testTransitionPolicyDispatch.cxx's
  // ConstructorFailureClearsRoadStartCellsAlongsideAllGroups already uses):
  // every transition individually satisfies the from/to ownership-parity
  // check (0,1,2 all owned), but the one cell's own hitSurfaces mask includes
  // surface 3 -- present in the catalog, untouched by either transition, and
  // deliberately not owned. Only reachable via direct corruption of the
  // topology's own per-cell metadata, exactly like the transition-id
  // corruption case this technique was first proven on.
  std::array<SurfaceDescriptor, 4> surfaces{
    surfaceWithOwner(0, SurfaceKind::Cylinder, 250), surfaceWithOwner(1, SurfaceKind::Cylinder, 250),
    surfaceWithOwner(2, SurfaceKind::Cylinder, 250), surfaceWithOwner(3, SurfaceKind::Cylinder, 250)};
  std::array<SurfaceTransition, 2> transitions{
    SurfaceTransition{SurfaceId{0}, SurfaceId{1}, SurfaceMask{}, 0},
    SurfaceTransition{SurfaceId{1}, SurfaceId{2}, SurfaceMask{}, 0}};
  SurfaceMask hitSurfaces = maskOf(0) | maskOf(1) | maskOf(2) | maskOf(3); // surface 3: deliberate corruption
  std::array<SurfaceCellTopology, 1> cells{SurfaceCellTopology{TransitionId{0}, TransitionId{1}, hitSurfaces}};
  std::array<uint32_t, 3> offsets{0, 1, 1};
  std::array<CellTopologyId, 1> byFirstTransition{CellTopologyId{0}};

  SurfaceMask cylinderSurfaces;
  for (const auto& descriptor : surfaces) {
    cylinderSurfaces.set(descriptor.id);
  }
  SurfaceGraphView layoutView{surfaces.data(), static_cast<uint32_t>(surfaces.size()), nullptr, 0, cylinderSurfaces, {}, transitions.data(), cells.data(), offsets.data(), byFirstTransition.data(), {}, static_cast<uint32_t>(transitions.size()), static_cast<uint32_t>(cells.size())};

  SurfaceMask owned = maskOf(0) | maskOf(1) | maskOf(2); // surface 3 not owned
  const std::vector<SurfaceId> order{SurfaceId{0}, SurfaceId{1}, SurfaceId{2}};
  const auto result = SurfacePlanBinding::build(layoutView, ClusterSourceId{0}, owned, order,
                                                SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::CrossBoundaryCell);
}

BOOST_AUTO_TEST_CASE(RejectsInvalidTopology)
{
  // Same corruption technique as RejectsCrossBoundaryCell above, but this
  // time the second cell references an out-of-range TransitionId, so
  // TransitionPolicyGrouping's own constructor rejects the whole schedule
  // (InvalidCellTransition) before SurfacePlanBinding::build() ever reaches
  // its own per-cell ownership checks.
  std::array<SurfaceDescriptor, 3> surfaces{
    surfaceWithOwner(0, SurfaceKind::Cylinder, 250), surfaceWithOwner(1, SurfaceKind::Cylinder, 250), surfaceWithOwner(2, SurfaceKind::Cylinder, 250)};
  std::array<SurfaceTransition, 2> transitions{
    SurfaceTransition{SurfaceId{0}, SurfaceId{1}, SurfaceMask{}, 0},
    SurfaceTransition{SurfaceId{1}, SurfaceId{2}, SurfaceMask{}, 0}};
  std::array<SurfaceCellTopology, 2> cells{
    SurfaceCellTopology{TransitionId{0}, TransitionId{1}, maskOf(0) | maskOf(1) | maskOf(2)}, // valid
    SurfaceCellTopology{TransitionId{0}, TransitionId{99}, SurfaceMask{}}};                   // out-of-range secondTransition
  std::array<uint32_t, 3> offsets{0, 2, 2};
  std::array<CellTopologyId, 2> byFirstTransition{CellTopologyId{0}, CellTopologyId{1}};

  SurfaceMask cylinderSurfaces;
  for (const auto& descriptor : surfaces) {
    cylinderSurfaces.set(descriptor.id);
  }
  SurfaceGraphView layoutView{surfaces.data(), static_cast<uint32_t>(surfaces.size()), nullptr, 0, cylinderSurfaces, {}, transitions.data(), cells.data(), offsets.data(), byFirstTransition.data(), {}, static_cast<uint32_t>(transitions.size()), static_cast<uint32_t>(cells.size())};

  SurfaceMask owned = maskOf(0) | maskOf(1) | maskOf(2);
  const auto result = SurfacePlanBinding::build(layoutView, ClusterSourceId{0}, owned, ordered(0, 3),
                                                SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidTopology);
}
