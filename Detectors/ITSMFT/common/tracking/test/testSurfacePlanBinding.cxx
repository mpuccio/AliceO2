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
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
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

SurfaceGraphDefinition appendDefinition(SurfaceGraphDefinition first, const SurfaceGraphDefinition& second)
{
  const auto offset = static_cast<uint16_t>(first.orderedSurfaces.size());
  first.orderedSurfaces.insert(first.orderedSurfaces.end(), second.orderedSurfaces.begin(), second.orderedSurfaces.end());
  for (const auto pair : second.basePairs) {
    first.basePairs.push_back(SurfaceAdjacencyPair{static_cast<uint16_t>(pair.fromIndex + offset),
                                                   static_cast<uint16_t>(pair.toIndex + offset)});
  }
  first.maxHoles = std::max(first.maxHoles, second.maxHoles);
  first.holeSurfaces |= second.holeSurfaces;
  first.seedingSurfaces |= second.seedingSurfaces;
  return first;
}

SurfaceDescriptor surfaceWithOwner(uint16_t id, SurfaceKind kind, uint8_t detectorId)
{
  return SurfaceDescriptor{SurfaceId{id}, id, detectorId, kind};
}

SurfaceLink adjacent(uint16_t from, uint16_t to)
{
  return SurfaceLink{SurfaceId{from}, SurfaceId{to}, SurfaceMask{}, 0};
}

struct BuiltLayout {
  SurfaceGraph layout;
  std::vector<SurfaceDescriptor> surfaces;

  SurfaceGraphView getView() const noexcept { return layout.getView(); }
};

BuiltLayout buildChainLayout(uint16_t nSurfaces, SurfaceKind kind, SurfaceMask seedingSurfaces = {})
{
  SurfaceGraph topology{nSurfaces, seedingSurfaces};
  std::vector<LinkId> links;
  for (uint16_t surface = 0; surface + 1 < nSurfaces; ++surface) {
    links.push_back(topology.addLink(adjacent(surface, surface + 1)));
  }
  for (size_t link = 0; link + 1 < links.size(); ++link) {
    topology.addCell(links[link], links[link + 1]);
  }
  BOOST_REQUIRE(topology.finalize());
  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t id = 0; id < nSurfaces; ++id) {
    surfaces.push_back(surfaceWithOwner(id, kind, 250));
  }
  return BuiltLayout{SurfaceGraph{surfaces, std::move(topology)}, std::move(surfaces)};
}

BuiltLayout buildTwoChainsLayout(uint16_t firstSurfaces, uint16_t secondSurfaces, SurfaceKind kind, SurfaceMask seedingSurfaces = {})
{
  const uint16_t totalSurfaces = static_cast<uint16_t>(firstSurfaces + secondSurfaces);
  SurfaceGraph topology{totalSurfaces, seedingSurfaces};
  std::vector<LinkId> firstLinks;
  std::vector<LinkId> secondLinks;
  for (uint16_t surface = 0; surface + 1 < firstSurfaces; ++surface) {
    firstLinks.push_back(topology.addLink(adjacent(surface, surface + 1)));
  }
  for (uint16_t surface = firstSurfaces; surface + 1 < totalSurfaces; ++surface) {
    secondLinks.push_back(topology.addLink(adjacent(surface, surface + 1)));
  }
  for (size_t link = 0; link + 1 < firstLinks.size(); ++link) {
    topology.addCell(firstLinks[link], firstLinks[link + 1]);
  }
  for (size_t link = 0; link + 1 < secondLinks.size(); ++link) {
    topology.addCell(secondLinks[link], secondLinks[link + 1]);
  }
  BOOST_REQUIRE(topology.finalize());
  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t id = 0; id < totalSurfaces; ++id) {
    surfaces.push_back(surfaceWithOwner(id, kind, 250));
  }
  return BuiltLayout{SurfaceGraph{surfaces, std::move(topology)}, std::move(surfaces)};
}

SurfaceMask allSurfaces(uint16_t nSurfaces)
{
  SurfaceMask result;
  for (uint16_t id = 0; id < nSurfaces; ++id) {
    result.set(SurfaceId{id});
  }
  return result;
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
    auto definition = makeSurfaceChain(ordered(0, ITSNLayers), 1, itsHoles, seeds & SurfaceMask{uint32_t{0x7f}});
    definition = appendDefinition(std::move(definition), makeSurfaceChain(ordered(ITSNLayers, MFTNLayers), 1, mftHoles,
                                                                          seeds & ~SurfaceMask{uint32_t{0x7f}}));
    SurfaceGraphBuilder builder{catalog, std::move(definition)};
    auto built = builder.build();
    BOOST_REQUIRE(built.ok());
    return std::move(*built.graph);
  }
};

SurfaceMask itsMask() { return SurfaceMask{uint32_t{0x7f}}; }
SurfaceMask mftMask() { return SurfaceMask{uint32_t{0x1ff80}}; }

void checkBindingCoversOwnedTopology(const SurfacePlanBinding& binding, const SurfaceGraphView& global)
{
  for (uint16_t s = 0; s < global.nSurfaces; ++s) {
    if (binding.getOwnedSurfaces().has(SurfaceId{s})) {
      BOOST_REQUIRE(binding.getOwnedSurfaceIndex(SurfaceId{s}));
    }
  }
  for (uint32_t t = 0; t < global.nLinks; ++t) {
    const auto id = LinkId{static_cast<uint16_t>(t)};
    const auto& link = global.getLink(id);
    if (binding.getOwnedSurfaces().has(link.from)) {
      BOOST_REQUIRE(binding.getOwnedSurfaces().has(link.to));
      BOOST_REQUIRE(binding.getScratchLinkSlot(id));
    }
  }
  for (uint32_t c = 0; c < global.nCells; ++c) {
    const auto id = CellTopologyId{static_cast<uint16_t>(c)};
    const auto& cell = global.getCell(id);
    const bool ownedLink = binding.getScratchLinkSlot(cell.firstLink).has_value() ||
                           binding.getScratchLinkSlot(cell.secondLink).has_value();
    if (ownedLink) {
      BOOST_REQUIRE(binding.getScratchCellSlot(id));
    }
  }
}

BOOST_AUTO_TEST_CASE(SurfacePlanBindingPreservesExactGlobalAndScheduledCellOrder)
{
  auto layout = buildChainLayout(7, SurfaceKind::Cylinder, maskOf(2) | maskOf(4));
  const auto order = ordered(0, 7);
  const auto result = SurfacePlanBinding::build(layout.getView(), allSurfaces(7), order);
  BOOST_REQUIRE(result.ok());

  const auto links = result.binding->getGlobalLinks();
  const auto cells = result.binding->getGlobalCells();
  const auto scheduled = result.binding->getGlobalScheduledCells();
  const auto roadStarts = result.binding->getGlobalRoadStartCells();
  BOOST_REQUIRE_EQUAL(links.size(), 6u);
  BOOST_REQUIRE_EQUAL(cells.size(), 5u);
  BOOST_REQUIRE_EQUAL(scheduled.size(), 5u);
  BOOST_REQUIRE_EQUAL(roadStarts.size(), 2u);
  for (uint16_t id = 0; id < links.size(); ++id) {
    BOOST_CHECK(links[id] == LinkId{id});
  }
  for (uint16_t id = 0; id < cells.size(); ++id) {
    BOOST_CHECK(cells[id] == CellTopologyId{id});
    BOOST_CHECK(scheduled[id] == CellTopologyId{id});
  }
  BOOST_CHECK(roadStarts[0] == CellTopologyId{0});
  BOOST_CHECK(roadStarts[1] == CellTopologyId{2});
}

BOOST_AUTO_TEST_CASE(SurfacePlanBindingPreservesRankOrderForNonMonotonicPlan)
{
  const std::vector<SurfaceId> order{SurfaceId{3}, SurfaceId{1}, SurfaceId{4}, SurfaceId{0}};
  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t id = 0; id < 5; ++id) {
    surfaces.push_back(surfaceWithOwner(id, SurfaceKind::Cylinder, 250));
  }
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}, makeSurfaceChain(order)};
  auto built = builder.build();
  BOOST_REQUIRE(built.ok());
  SurfaceMask owned;
  for (const auto surface : order) {
    owned.set(surface);
  }
  const auto result = SurfacePlanBinding::build(built.graph->getView(), owned, order);
  BOOST_REQUIRE(result.ok());
  const auto scheduled = result.binding->getGlobalScheduledCells();
  BOOST_REQUIRE_EQUAL(scheduled.size(), 2u);
  BOOST_CHECK(scheduled[0] == CellTopologyId{0});
  BOOST_CHECK(scheduled[1] == CellTopologyId{1});
}

BOOST_AUTO_TEST_CASE(SurfacePlanBindingSeparatesDisconnectedDetectorSchedules)
{
  auto layout = buildTwoChainsLayout(7, 10, SurfaceKind::Cylinder, maskOf(2) | maskOf(9));
  const auto firstResult = SurfacePlanBinding::build(layout.getView(), allSurfaces(7), ordered(0, 7));
  SurfaceMask second;
  for (uint16_t id = 7; id < 17; ++id) {
    second.set(SurfaceId{id});
  }
  const auto secondResult = SurfacePlanBinding::build(layout.getView(), second, ordered(7, 10));
  BOOST_REQUIRE(firstResult.ok());
  BOOST_REQUIRE(secondResult.ok());
  BOOST_REQUIRE_EQUAL(firstResult.binding->getGlobalScheduledCells().size(), 5u);
  BOOST_REQUIRE_EQUAL(secondResult.binding->getGlobalScheduledCells().size(), 8u);
  BOOST_CHECK(firstResult.binding->getGlobalScheduledCells()[0] == CellTopologyId{0});
  BOOST_CHECK(secondResult.binding->getGlobalScheduledCells()[0] == CellTopologyId{5});
  BOOST_REQUIRE_EQUAL(firstResult.binding->getGlobalRoadStartCells().size(), 1u);
  BOOST_REQUIRE_EQUAL(secondResult.binding->getGlobalRoadStartCells().size(), 1u);
  BOOST_CHECK(firstResult.binding->getGlobalRoadStartCells()[0] == CellTopologyId{0});
  BOOST_CHECK(secondResult.binding->getGlobalRoadStartCells()[0] == CellTopologyId{5});
}

BOOST_AUTO_TEST_CASE(SurfacePlanBindingRejectsCyclicTopology)
{
  SurfaceGraph topology{3};
  BOOST_REQUIRE(topology.addLink(adjacent(0, 1)).isValid());
  BOOST_REQUIRE(topology.addLink(adjacent(1, 2)).isValid());
  BOOST_REQUIRE(topology.addLink(adjacent(2, 0)).isValid());
  BOOST_REQUIRE(topology.finalize());
  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t id = 0; id < 3; ++id) {
    surfaces.push_back(surfaceWithOwner(id, SurfaceKind::Cylinder, 250));
  }
  SurfaceGraph layout{surfaces, std::move(topology)};
  const auto result = SurfacePlanBinding::build(layout.getView(), allSurfaces(3), ordered(0, 3));
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidTopology);
}

} // namespace

BOOST_AUTO_TEST_CASE(SurfacePlanBindingMapsCombinedItsAndMftPlans)
{
  CombinedLayout combined;

  const auto its = SurfacePlanBinding::build(combined.view, itsMask(), ordered(0, ITSNLayers));
  const auto mft = SurfacePlanBinding::build(combined.view, mftMask(), ordered(ITSNLayers, MFTNLayers));
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
  // must build successfully given internally consistent expected kind/configuration.
  constexpr uint8_t kSyntheticDetectorId = 250;
  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t id = 0; id < 4; ++id) {
    surfaces.push_back(surfaceWithOwner(id, SurfaceKind::Cylinder, kSyntheticDetectorId));
  }
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())},
                              makeSurfaceChain(ordered(0, 4), 0, SurfaceMask{}, maskOf(3))};
  auto built = builder.build();
  BOOST_REQUIRE(built.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = built.graph->getView();

  SurfaceMask owned;
  for (uint16_t id = 0; id < 4; ++id) {
    owned.set(SurfaceId{id});
  }
  const auto result = SurfacePlanBinding::build(view, owned, ordered(0, 4));
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result.binding->getGlobalLinks().size(), 3u);
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
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())},
                              makeSurfaceChain(planOrder, 0, SurfaceMask{}, maskOf(5))};
  auto built = builder.build();
  BOOST_REQUIRE(built.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = built.graph->getView();
  const auto bindingResult = SurfacePlanBinding::build(view, owned, planOrder);
  BOOST_REQUIRE(bindingResult.ok());
  const auto& binding = *bindingResult.binding;

  BOOST_CHECK(std::equal(binding.getOrderedSurfaces().begin(), binding.getOrderedSurfaces().end(), planOrder.begin(), planOrder.end()));
  for (std::size_t position = 0; position < planOrder.size(); ++position) {
    BOOST_REQUIRE(binding.getOwnedSurfaceIndex(planOrder[position]));
    BOOST_CHECK_EQUAL(*binding.getOwnedSurfaceIndex(planOrder[position]), position);
  }

  SurfaceTrackingScratch scratch;
  scratch.adoptPlan(binding.getOrderedSurfaces().size(), binding.getGlobalLinks().size(), binding.getGlobalCells().size());
  BOOST_CHECK_EQUAL(scratch.getNOwnedSurfaces(), planOrder.size());
  BOOST_CHECK_EQUAL(scratch.getNLinks(), binding.getGlobalLinks().size());
  BOOST_CHECK_EQUAL(scratch.getNCells(), binding.getGlobalCells().size());
  for (std::size_t slot = 0; slot < binding.getGlobalLinks().size(); ++slot) {
    BOOST_REQUIRE(binding.getScratchLinkSlot(binding.getGlobalLinks()[slot]));
    BOOST_CHECK_EQUAL(*binding.getScratchLinkSlot(binding.getGlobalLinks()[slot]), slot);
  }
  for (std::size_t slot = 0; slot < binding.getGlobalCells().size(); ++slot) {
    BOOST_REQUIRE(binding.getScratchCellSlot(binding.getGlobalCells()[slot]));
    BOOST_CHECK_EQUAL(*binding.getScratchCellSlot(binding.getGlobalCells()[slot]), slot);
  }

  std::vector<std::vector<GlobalMeasurement>> measurements(planOrder.size());
  std::vector<gsl::span<const GlobalMeasurement>> measurementViews;
  measurementViews.reserve(planOrder.size());
  for (std::size_t position = 0; position < planOrder.size(); ++position) {
    GlobalMeasurement measurement;
    measurement.surface = binding.getOrderedSurfaces()[position];
    measurement.cluster = ClusterRef{ClusterSourceId{7}, static_cast<uint32_t>(100 + position)};
    measurements[position].push_back(measurement);
    measurementViews.emplace_back(measurements[position]);
  }
  for (std::size_t position = 0; position < measurementViews.size(); ++position) {
    BOOST_CHECK_EQUAL(measurementViews[position][0].surface.value(), planOrder[position].value());
    BOOST_CHECK(measurementViews[position][0].cluster.source == ClusterSourceId{7});
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

BOOST_AUTO_TEST_CASE(RejectsSurfaceMaskSizeMismatch)
{
  CombinedLayout combined;
  const auto result = SurfacePlanBinding::build(combined.view, itsMask(), ordered(0, ITSNLayers - 1));
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidSurfaceMask);
}

BOOST_AUTO_TEST_CASE(RejectsSurfaceMaskNotASubsetOfLayout)
{
  std::vector<SurfaceDescriptor> surfaces{surfaceWithOwner(0, SurfaceKind::Cylinder, 250), surfaceWithOwner(1, SurfaceKind::Cylinder, 250)};
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}, makeSurfaceChain(ordered(0, 2))};
  auto built = builder.build();
  BOOST_REQUIRE(built.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = built.graph->getView();

  SurfaceMask outOfRange = maskOf(0) | maskOf(1) | maskOf(5); // surface 5 does not exist in this 2-surface layout
  const std::vector<SurfaceId> order{SurfaceId{0}, SurfaceId{1}, SurfaceId{5}};
  const auto result = SurfacePlanBinding::build(view, outOfRange, order);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidSurfaceMask);
}

BOOST_AUTO_TEST_CASE(RejectsDuplicateSurfaceInOrder)
{
  CombinedLayout combined;
  auto order = ordered(0, ITSNLayers);
  order[1] = order[0]; // duplicate
  const auto result = SurfacePlanBinding::build(combined.view, itsMask(), order);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidLegacySurfaceOrder);
}

BOOST_AUTO_TEST_CASE(RejectsUnownedSurfaceInOrder)
{
  CombinedLayout combined;
  SurfaceMask owned = maskOf(0);
  const std::vector<SurfaceId> order{SurfaceId{1}}; // count matches (1) but surface 1 is not owned
  const auto result = SurfacePlanBinding::build(combined.view, owned, order);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidLegacySurfaceOrder);
}

BOOST_AUTO_TEST_CASE(SurfacePlanBindingBuildsAcrossMultipleDistinctDetectorIdentitiesInOneOwnedSet)
{
  // SurfacePlanBinding must be generic over its own owned SurfaceId set: it
  // must not assume "one binding, one detector". Three compatible (Cylinder)
  // surfaces spanning two distinct synthetic detectorIds (250, 251), one
  // valid source, consistent expected kind/configuration, and an internally valid
  // chain topology -- this must build successfully, unlike
  // The former binding's own single-detector-scoped contract.
  std::vector<SurfaceDescriptor> surfaces{
    surfaceWithOwner(0, SurfaceKind::Cylinder, 250),
    surfaceWithOwner(1, SurfaceKind::Cylinder, 250),
    surfaceWithOwner(2, SurfaceKind::Cylinder, 251)}; // second detectorId, still owned by the same binding
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())},
                              makeSurfaceChain(ordered(0, 3), 0, SurfaceMask{}, maskOf(2))};
  auto built = builder.build();
  BOOST_REQUIRE(built.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = built.graph->getView();

  SurfaceMask owned = maskOf(0) | maskOf(1) | maskOf(2);
  const auto result = SurfacePlanBinding::build(view, owned, ordered(0, 3));
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result.binding->getGlobalLinks().size(), 2u);
  BOOST_CHECK_EQUAL(result.binding->getGlobalCells().size(), 1u);
  for (uint16_t id = 0; id < 3; ++id) {
    BOOST_REQUIRE(result.binding->getOwnedSurfaceIndex(SurfaceId{id}));
    BOOST_CHECK_EQUAL(*result.binding->getOwnedSurfaceIndex(SurfaceId{id}), id);
  }
}

BOOST_AUTO_TEST_CASE(RejectsCrossBoundaryLink)
{
  std::vector<SurfaceDescriptor> surfaces{surfaceWithOwner(0, SurfaceKind::Cylinder, 250), surfaceWithOwner(1, SurfaceKind::Cylinder, 250)};
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}, makeSurfaceChain(ordered(0, 2))};
  auto built = builder.build();
  BOOST_REQUIRE(built.ok());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = built.graph->getView();

  // Own only surface 0: the 0->1 link has fromOwned=true, toOwned=false.
  SurfaceMask owned = maskOf(0);
  const std::vector<SurfaceId> order{SurfaceId{0}};
  const auto result = SurfacePlanBinding::build(view, owned, order);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::CrossBoundaryLink);
}

BOOST_AUTO_TEST_CASE(RejectsCrossBoundaryCell)
{
  // Hand-built raw view (the graph-view types are trivially-copyable PODs):
  // every link individually satisfies the from/to ownership-parity
  // check (0,1,2 all owned), but the one cell's own hitSurfaces mask includes
  // surface 3 -- present in the catalog, untouched by either link, and
  // deliberately not owned. Only reachable via direct corruption of the
  // topology's own per-cell metadata, exactly like the link-id
  // corruption case this technique was first proven on.
  std::array<SurfaceDescriptor, 4> surfaces{
    surfaceWithOwner(0, SurfaceKind::Cylinder, 250), surfaceWithOwner(1, SurfaceKind::Cylinder, 250),
    surfaceWithOwner(2, SurfaceKind::Cylinder, 250), surfaceWithOwner(3, SurfaceKind::Cylinder, 250)};
  std::array<SurfaceLink, 2> links{
    SurfaceLink{SurfaceId{0}, SurfaceId{1}, SurfaceMask{}, 0},
    SurfaceLink{SurfaceId{1}, SurfaceId{2}, SurfaceMask{}, 0}};
  SurfaceMask hitSurfaces = maskOf(0) | maskOf(1) | maskOf(2) | maskOf(3); // surface 3: deliberate corruption
  std::array<SurfaceCellTopology, 1> cells{SurfaceCellTopology{LinkId{0}, LinkId{1}, hitSurfaces}};
  std::array<uint32_t, 3> offsets{0, 1, 1};
  std::array<CellTopologyId, 1> byFirstLink{CellTopologyId{0}};

  SurfaceMask cylinderSurfaces;
  for (const auto& descriptor : surfaces) {
    cylinderSurfaces.set(descriptor.id);
  }
  SurfaceGraphView layoutView{surfaces.data(), static_cast<uint32_t>(surfaces.size()), nullptr, 0, cylinderSurfaces, {}, links.data(), cells.data(), offsets.data(), byFirstLink.data(), {}, static_cast<uint32_t>(links.size()), static_cast<uint32_t>(cells.size())};

  SurfaceMask owned = maskOf(0) | maskOf(1) | maskOf(2); // surface 3 not owned
  const std::vector<SurfaceId> order{SurfaceId{0}, SurfaceId{1}, SurfaceId{2}};
  const auto result = SurfacePlanBinding::build(layoutView, owned, order);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::CrossBoundaryCell);
}

BOOST_AUTO_TEST_CASE(RejectsInvalidTopology)
{
  // Same corruption technique as RejectsCrossBoundaryCell above, but this
  // time the second cell references an out-of-range LinkId, so
  // SurfacePlanBinding's topology validation rejects the whole schedule
  // (InvalidCellLink) before SurfacePlanBinding::build() ever reaches
  // its own per-cell ownership checks.
  std::array<SurfaceDescriptor, 3> surfaces{
    surfaceWithOwner(0, SurfaceKind::Cylinder, 250), surfaceWithOwner(1, SurfaceKind::Cylinder, 250), surfaceWithOwner(2, SurfaceKind::Cylinder, 250)};
  std::array<SurfaceLink, 2> links{
    SurfaceLink{SurfaceId{0}, SurfaceId{1}, SurfaceMask{}, 0},
    SurfaceLink{SurfaceId{1}, SurfaceId{2}, SurfaceMask{}, 0}};
  std::array<SurfaceCellTopology, 2> cells{
    SurfaceCellTopology{LinkId{0}, LinkId{1}, maskOf(0) | maskOf(1) | maskOf(2)}, // valid
    SurfaceCellTopology{LinkId{0}, LinkId{99}, SurfaceMask{}}};                   // out-of-range secondLink
  std::array<uint32_t, 3> offsets{0, 2, 2};
  std::array<CellTopologyId, 2> byFirstLink{CellTopologyId{0}, CellTopologyId{1}};

  SurfaceMask cylinderSurfaces;
  for (const auto& descriptor : surfaces) {
    cylinderSurfaces.set(descriptor.id);
  }
  SurfaceGraphView layoutView{surfaces.data(), static_cast<uint32_t>(surfaces.size()), nullptr, 0, cylinderSurfaces, {}, links.data(), cells.data(), offsets.data(), byFirstLink.data(), {}, static_cast<uint32_t>(links.size()), static_cast<uint32_t>(cells.size())};

  SurfaceMask owned = maskOf(0) | maskOf(1) | maskOf(2);
  const auto result = SurfacePlanBinding::build(layoutView, owned, ordered(0, 3));
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == SurfacePlanBindingError::InvalidTopology);
}
