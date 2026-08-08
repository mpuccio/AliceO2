// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details.
// All rights not expressly granted are reserved.

#define BOOST_TEST_MODULE ITSMFT pair - list graph prototype
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cstring>
#include <type_traits>
#include <vector>

#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "PairListGraphPrototype.h"

namespace
{
using namespace o2::itsmft::tracking;
using namespace o2::itsmft::tracking::test;

static_assert(std::is_trivially_copyable_v<SurfaceTransition>);
static_assert(std::is_trivially_copyable_v<SurfaceCellTopology>);
static_assert(std::is_trivially_copyable_v<TransitionId>);
static_assert(std::is_trivially_copyable_v<CellTopologyId>);
static_assert(std::is_trivially_copyable_v<SurfaceMask>);
static_assert(std::is_trivially_copyable_v<uint32_t>);

template <typename T>
void checkBytes(gsl::span<const T> actual, gsl::span<const T> expected)
{
  BOOST_REQUIRE_EQUAL(actual.size(), expected.size());
  if (actual.empty()) {
    return;
  }
  if constexpr (std::is_same_v<T, SurfaceTransition>) {
    for (size_t i = 0; i < actual.size(); ++i) {
      BOOST_CHECK(actual[i].from == expected[i].from);
      BOOST_CHECK(actual[i].to == expected[i].to);
      BOOST_CHECK(actual[i].skippedSurfaces == expected[i].skippedSurfaces);
      BOOST_CHECK_EQUAL(actual[i].flags, expected[i].flags);
    }
    // SurfaceTransition has two tail padding bytes (its value fields total
    // ten bytes). Normalize those unspecified bytes before the raw comparison;
    // all value fields above are still checked and the type is asserted
    // trivially copyable at file scope.
    std::vector<SurfaceTransition> canonicalActual(actual.size());
    std::vector<SurfaceTransition> canonicalExpected(expected.size());
    for (size_t i = 0; i < actual.size(); ++i) {
      canonicalActual[i].from = actual[i].from;
      canonicalActual[i].to = actual[i].to;
      canonicalActual[i].skippedSurfaces = actual[i].skippedSurfaces;
      canonicalActual[i].flags = actual[i].flags;
      canonicalExpected[i].from = expected[i].from;
      canonicalExpected[i].to = expected[i].to;
      canonicalExpected[i].skippedSurfaces = expected[i].skippedSurfaces;
      canonicalExpected[i].flags = expected[i].flags;
    }
    BOOST_CHECK_EQUAL(std::memcmp(canonicalActual.data(), canonicalExpected.data(), canonicalActual.size() * sizeof(SurfaceTransition)), 0);
    return;
  }
  BOOST_CHECK_EQUAL(std::memcmp(actual.data(), expected.data(), actual.size_bytes()), 0);
}

template <typename T>
void checkBytes(const std::vector<T>& actual, const std::vector<T>& expected)
{
  checkBytes(gsl::span<const T>{actual}, gsl::span<const T>{expected});
}

SurfaceMask maskOf(std::initializer_list<uint16_t> ids)
{
  SurfaceMask result;
  for (const auto id : ids) {
    result.set(SurfaceId{id});
  }
  return result;
}

std::vector<SurfaceDescriptor> catalog(uint16_t count, SurfaceKind kind = SurfaceKind::Cylinder)
{
  std::vector<SurfaceDescriptor> result;
  for (uint16_t id = 0; id < count; ++id) {
    result.push_back(SurfaceDescriptor{SurfaceId{id}, id, 0, kind});
  }
  return result;
}

std::vector<SurfaceId> ids(std::initializer_list<uint16_t> values)
{
  std::vector<SurfaceId> result;
  for (const auto value : values) {
    result.emplace_back(value);
  }
  return result;
}

PairListComponentInput component(std::initializer_list<uint16_t> values)
{
  PairListComponentInput result;
  result.activeSurfaces = ids(values);
  for (uint16_t index = 0; index + 1 < result.activeSurfaces.size(); ++index) {
    result.basePairs.push_back(PairListBasePair{index, static_cast<uint16_t>(index + 1)});
  }
  return result;
}

PairListGraphInput input(std::initializer_list<PairListComponentInput> components,
                         PairListHolePolicy holePolicy = {}, SurfaceMask seedingMask = {})
{
  return PairListGraphInput{std::vector<PairListComponentInput>{components}, holePolicy, seedingMask};
}

SurfaceGraphBuildResult buildCurrent(const std::vector<SurfaceDescriptor>& surfaces,
                                     const PairListGraphInput& input)
{
  SurfaceGraphDefinition definition;
  for (const auto& item : input.components) {
    SurfaceMask active;
    for (const auto id : item.activeSurfaces) {
      active.set(id);
    }
    const auto offset = static_cast<uint16_t>(definition.orderedSurfaces.size());
    const auto componentDefinition = makeSurfaceChain(item.activeSurfaces, input.holePolicy.maxSkipped,
                                                      input.holePolicy.skippableSurfaceMask & active,
                                                      input.seedingMask & active);
    definition.orderedSurfaces.insert(definition.orderedSurfaces.end(), componentDefinition.orderedSurfaces.begin(), componentDefinition.orderedSurfaces.end());
    for (const auto pair : componentDefinition.basePairs) {
      definition.basePairs.push_back(SurfaceAdjacencyPair{static_cast<uint16_t>(pair.fromIndex + offset),
                                                          static_cast<uint16_t>(pair.toIndex + offset)});
    }
    definition.maxHoles = std::max(definition.maxHoles, componentDefinition.maxHoles);
    definition.holeSurfaces |= componentDefinition.holeSurfaces;
    definition.seedingSurfaces |= componentDefinition.seedingSurfaces;
  }
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}, std::move(definition)};
  return builder.build();
}

void compareTopology(const PairListGraph& prototype, const SurfaceGraphView& current)
{
  {
    BOOST_TEST_CONTEXT("transitions")
    {
      checkBytes(prototype.transitions, std::vector<SurfaceTransition>(current.transitions, current.transitions + current.nTransitions));
    }
  }
  std::vector<SurfaceMask> expectedWitnesses;
  for (uint32_t id = 0; id < current.nTransitions; ++id) {
    expectedWitnesses.push_back(current.getTransition(TransitionId{static_cast<uint16_t>(id)}).skippedSurfaces);
  }
  {
    BOOST_TEST_CONTEXT("witnesses") { checkBytes(prototype.skippedWitnesses, expectedWitnesses); }
    BOOST_TEST_CONTEXT("cells") { checkBytes(prototype.cells, std::vector<SurfaceCellTopology>(current.cells, current.cells + current.nCells)); }
    BOOST_TEST_CONTEXT("offsets")
    {
      checkBytes(prototype.cellsByFirstTransitionOffsets,
                 std::vector<uint32_t>(current.cellsByFirstTransitionOffsets, current.cellsByFirstTransitionOffsets + current.nTransitions + 1));
    }
    BOOST_TEST_CONTEXT("csr entries")
    {
      checkBytes(prototype.cellsByFirstTransition,
                 std::vector<CellTopologyId>(current.cellsByFirstTransition, current.cellsByFirstTransition + current.nCells));
    }
  }
}

void compareBinding(const PairListGraph& prototype, const SurfaceGraphView& view, const PairListGraphInput& input)
{
  SurfaceMask owned;
  std::vector<SurfaceId> orderedSurfaces;
  for (const auto& component : input.components) {
    for (const auto id : component.activeSurfaces) {
      owned.set(id);
      orderedSurfaces.push_back(id);
    }
  }
  const auto bindingResult = SurfacePlanBinding::build(view, owned, orderedSurfaces);
  BOOST_REQUIRE(bindingResult.ok());
  std::vector<TransitionId> expectedTransitions;
  std::vector<CellTopologyId> expectedCells;
  for (uint16_t id = 0; id < prototype.transitions.size(); ++id) {
    expectedTransitions.emplace_back(id);
  }
  for (uint16_t id = 0; id < prototype.cells.size(); ++id) {
    expectedCells.emplace_back(id);
  }
  checkBytes(expectedTransitions, std::vector<TransitionId>(bindingResult.binding->getGlobalTransitions().begin(),
                                                            bindingResult.binding->getGlobalTransitions().end()));
  checkBytes(expectedCells, std::vector<CellTopologyId>(bindingResult.binding->getGlobalCells().begin(),
                                                        bindingResult.binding->getGlobalCells().end()));
  checkBytes(prototype.scheduledCells, std::vector<CellTopologyId>(bindingResult.binding->getGlobalScheduledCells().begin(),
                                                                   bindingResult.binding->getGlobalScheduledCells().end()));
  checkBytes(prototype.roadStartCells, std::vector<CellTopologyId>(bindingResult.binding->getGlobalRoadStartCells().begin(),
                                                                   bindingResult.binding->getGlobalRoadStartCells().end()));
}

void checkCase(const std::vector<SurfaceDescriptor>& surfaces, const PairListGraphInput& input)
{
  const auto prototypeResult = derivePairListGraph(input);
  BOOST_REQUIRE(prototypeResult.ok());
  const auto currentResult = buildCurrent(surfaces, input);
  BOOST_REQUIRE(currentResult.ok());
  compareTopology(*prototypeResult.graph, currentResult.graph->getView());
  const auto view = currentResult.graph->getView();
  compareBinding(*prototypeResult.graph, view, input);
}
} // namespace

BOOST_AUTO_TEST_CASE(ITSLikeSevenCylinderChainNoHoles)
{
  checkCase(catalog(7), input({component({0, 1, 2, 3, 4, 5, 6})}, {}, maskOf({2, 4})));
}

BOOST_AUTO_TEST_CASE(ITSLikeSevenCylinderChainOneHole)
{
  checkCase(catalog(7), input({component({0, 1, 2, 3, 4, 5, 6})}, {1, maskOf({3})}, maskOf({0, 5})));
}

BOOST_AUTO_TEST_CASE(MFTLikeTenDiskChainNoHoles)
{
  checkCase(catalog(10, SurfaceKind::Disk), input({component({0, 1, 2, 3, 4, 5, 6, 7, 8, 9})}));
}

BOOST_AUTO_TEST_CASE(MFTLikeTenDiskChainOneHole)
{
  checkCase(catalog(10, SurfaceKind::Disk), input({component({0, 1, 2, 3, 4, 5, 6, 7, 8, 9})}, {1, maskOf({5})}, maskOf({1, 8})));
}

BOOST_AUTO_TEST_CASE(CombinedDisconnectedComponentsNoHoles)
{
  auto graphInput = input({component({0, 1, 2, 3, 4, 5, 6}),
                           component({7, 8, 9, 10, 11, 12, 13, 14, 15, 16})},
                          {}, maskOf({2, 9, 15}));
  auto surfaces = catalog(17);
  for (uint16_t id = 7; id < 17; ++id) {
    surfaces[id].kind = SurfaceKind::Disk;
  }
  const auto prototypeResult = derivePairListGraph(graphInput);
  BOOST_REQUIRE(prototypeResult.ok());
  const auto currentResult = buildCurrent(surfaces, graphInput);
  BOOST_REQUIRE(currentResult.ok());
  compareTopology(*prototypeResult.graph, currentResult.graph->getView());
  compareBinding(*prototypeResult.graph, currentResult.graph->getView(), graphInput);
}

BOOST_AUTO_TEST_CASE(IndependentAuthoritiesChangeOnlyTheirDerivedOutputs)
{
  const auto active = component({0, 1, 2, 3, 4});
  const PairListHolePolicy oneHole{1, maskOf({2})};
  const auto noSeed = input({active}, oneHole, {});
  const auto seed = input({active}, oneHole, maskOf({3}));
  const auto noSeedResult = derivePairListGraph(noSeed);
  const auto seedResult = derivePairListGraph(seed);
  BOOST_REQUIRE(noSeedResult.ok());
  BOOST_REQUIRE(seedResult.ok());
  checkBytes(noSeedResult.graph->transitions, seedResult.graph->transitions);
  checkBytes(noSeedResult.graph->cells, seedResult.graph->cells);
  checkBytes(noSeedResult.graph->cellsByFirstTransitionOffsets, seedResult.graph->cellsByFirstTransitionOffsets);
  checkBytes(noSeedResult.graph->cellsByFirstTransition, seedResult.graph->cellsByFirstTransition);
  checkBytes(noSeedResult.graph->scheduledCells, seedResult.graph->scheduledCells);
  BOOST_CHECK(noSeedResult.graph->roadStartCells != seedResult.graph->roadStartCells);

  const auto noHole = input({active}, {}, maskOf({3}));
  const auto holeResult = derivePairListGraph(noHole);
  const auto expandedResult = derivePairListGraph(seed);
  BOOST_REQUIRE(holeResult.ok());
  BOOST_REQUIRE(expandedResult.ok());
  checkBytes(noHole.components[0].basePairs, seed.components[0].basePairs);
  BOOST_CHECK_NE(holeResult.graph->transitions.size(), expandedResult.graph->transitions.size());
  auto catalogFive = catalog(5);
  const auto noHoleCurrent = buildCurrent(catalogFive, noHole);
  const auto expandedCurrent = buildCurrent(catalogFive, seed);
  BOOST_REQUIRE(noHoleCurrent.ok());
  BOOST_REQUIRE(expandedCurrent.ok());
  compareTopology(*holeResult.graph, noHoleCurrent.graph->getView());
  compareTopology(*expandedResult.graph, expandedCurrent.graph->getView());

  auto missingEdgeComponent = component({0, 1, 2, 3});
  missingEdgeComponent.basePairs.erase(missingEdgeComponent.basePairs.begin() + 1);
  const auto missingEdge = input({missingEdgeComponent}, oneHole, maskOf({3}));
  const auto missingEdgeResult = derivePairListGraph(missingEdge);
  BOOST_REQUIRE(missingEdgeResult.ok());
  BOOST_CHECK_EQUAL(missingEdgeResult.graph->transitions.size(), 2u);
  for (const auto transition : missingEdgeResult.graph->transitions) {
    BOOST_CHECK(!(transition.from == SurfaceId{0} && transition.to == SurfaceId{2}));
    BOOST_CHECK(!(transition.from == SurfaceId{1} && transition.to == SurfaceId{3}));
  }

  const auto changedActive = derivePairListGraph(input({component({0, 1, 2, 4, 5})}, oneHole, maskOf({2})));
  BOOST_REQUIRE(changedActive.ok());
  BOOST_CHECK(changedActive.graph->transitions[3].to != noSeedResult.graph->transitions[3].to);
}

BOOST_AUTO_TEST_CASE(NonMonotonicGlobalSurfaceIds)
{
  checkCase(catalog(5), input({component({3, 1, 4, 0})}));
}

BOOST_AUTO_TEST_CASE(EmptyTopology)
{
  const PairListGraphInput input;
  const auto result = derivePairListGraph(input);
  BOOST_REQUIRE(result.ok());
  const auto current = buildCurrent({}, input);
  BOOST_REQUIRE(current.ok());
  compareTopology(*result.graph, current.graph->getView());
  BOOST_CHECK(result.graph->scheduledCells.empty());
  BOOST_CHECK(result.graph->roadStartCells.empty());
}

BOOST_AUTO_TEST_CASE(MultipleAllowedHolesAndIndependentMasks)
{
  checkCase(catalog(7), input({component({0, 1, 2, 3, 4, 5, 6})}, {2, maskOf({2, 4})}, maskOf({0, 6})));
  auto rejected = component({0, 1, 2, 3});
  rejected.basePairs.pop_back();
  const auto result = derivePairListGraph(input({std::move(rejected)}, {1, maskOf({2})}, maskOf({0})));
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result.graph->transitions.size(), 2u);
  BOOST_CHECK(result.graph->transitions[0].from == SurfaceId{0});
  BOOST_CHECK(result.graph->transitions[0].to == SurfaceId{1});
}

BOOST_AUTO_TEST_CASE(ValidationRejectsIndependentAuthorityViolations)
{
  PairListComponentInput empty;
  BOOST_CHECK(derivePairListGraph(input({std::move(empty)})).error == PairListGraphError::EmptyComponent);
  PairListComponentInput invalidId;
  invalidId.activeSurfaces = {SurfaceId::invalid()};
  BOOST_CHECK(derivePairListGraph(input({std::move(invalidId)})).error == PairListGraphError::InvalidSurfaceId);
  PairListComponentInput outOfRangeId;
  outOfRangeId.activeSurfaces = {SurfaceId{MaxLayoutSurfaces}};
  BOOST_CHECK(derivePairListGraph(input({std::move(outOfRangeId)})).error == PairListGraphError::InvalidSurfaceId);
  BOOST_CHECK(derivePairListGraph(input({component({0, 1})}, {-1, {}})).error == PairListGraphError::NegativeHoleBudget);
  BOOST_CHECK(derivePairListGraph(input({component({0, 1})}, {0, maskOf({3})})).error == PairListGraphError::HoleMaskOutsideActive);
  BOOST_CHECK(derivePairListGraph(input({component({0, 1})}, {}, maskOf({3}))).error == PairListGraphError::SeedingMaskOutsideActive);
  auto duplicate = component({0, 1});
  duplicate.activeSurfaces.push_back(SurfaceId{1});
  BOOST_CHECK(derivePairListGraph(input({std::move(duplicate)})).error == PairListGraphError::DuplicateActiveSurface);
  auto malformed = component({0, 1, 2});
  malformed.basePairs = {{0, 2}};
  BOOST_CHECK(derivePairListGraph(input({std::move(malformed)})).error == PairListGraphError::InvalidPair);
  PairListComponentInput selfPair;
  selfPair.activeSurfaces = ids({0, 1});
  selfPair.basePairs = {{0, 0}};
  BOOST_CHECK(derivePairListGraph(input({std::move(selfPair)})).error == PairListGraphError::InvalidPair);
  auto duplicatePair = component({0, 1});
  duplicatePair.basePairs.push_back({0, 1});
  BOOST_CHECK(derivePairListGraph(input({std::move(duplicatePair)})).error == PairListGraphError::DuplicatePair);
  PairListComponentInput cross;
  cross.activeSurfaces = ids({0, 1});
  cross.basePairs = {{1, 0}};
  BOOST_CHECK(derivePairListGraph(input({std::move(cross)})).error == PairListGraphError::InvalidPair);
}
