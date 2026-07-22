// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT DetectorLayoutBuilder
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/TrackingTopology.h"

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

void checkCsrConsistency(const SparseTrackingTopologyView& view)
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

// Compares a built SparseTrackingTopologyView against the legacy templated
// TrackingTopology<NLayers>::View it should be exactly equivalent to, under
// the identity mapping SurfaceId{i} <-> layer i (i.e. a monotonic,
// hole-free-in-position orderedSurfaces = {0, 1, ..., NLayers-1}).
template <int NLayers>
void checkParity(const SparseTrackingTopologyView& built, const typename TrackingTopology<NLayers>::View& reference)
{
  BOOST_REQUIRE_EQUAL(built.nTransitions, static_cast<uint32_t>(reference.nTransitions));
  for (uint32_t t = 0; t < built.nTransitions; ++t) {
    const auto& builtTransition = built.getTransition(TransitionId{static_cast<uint16_t>(t)});
    const auto& refTransition = reference.getTransition(t);
    BOOST_CHECK_EQUAL(builtTransition.from.value(), refTransition.fromLayer);
    BOOST_CHECK_EQUAL(builtTransition.to.value(), refTransition.toLayer);
  }

  BOOST_REQUIRE_EQUAL(built.nCells, static_cast<uint32_t>(reference.nCells));
  for (uint32_t c = 0; c < built.nCells; ++c) {
    const auto& builtCell = built.getCell(CellTopologyId{static_cast<uint16_t>(c)});
    const auto& refCell = reference.getCell(c);
    BOOST_CHECK_EQUAL(builtCell.firstTransition.value(), refCell.firstTransition);
    BOOST_CHECK_EQUAL(builtCell.secondTransition.value(), refCell.secondTransition);
    BOOST_CHECK_EQUAL(builtCell.hitSurfaces.value(), static_cast<uint32_t>(refCell.hitLayerMask.value()));
  }

  checkCsrConsistency(built);
}
} // namespace

BOOST_AUTO_TEST_CASE(TraversalFollowsSuppliedOrderNotNumericSurfaceId)
{
  // Catalog ids are dense (0..4) but the physical traversal order is the
  // deliberately non-monotonic chain 3 -> 1 -> 4 -> 0; surface 2 is present
  // in the catalog but not activated by this subgraph.
  DetectorLayoutBuilder builder{denseCatalog(5)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({3, 1, 4, 0}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});

  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());
  const auto view = result.layout->getTopology().getView();

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

BOOST_AUTO_TEST_CASE(ExactParityWithTrackingTopologySevenNoHoles)
{
  DetectorLayoutBuilder builder{denseCatalog(7)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1, 2, 3, 4, 5, 6}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});
  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());

  TrackingTopology<7> reference;
  reference.init(7, 0, LayerMask{0});

  checkParity<7>(result.layout->getTopology().getView(), reference.getView());
}

BOOST_AUTO_TEST_CASE(ExactParityWithTrackingTopologySevenSingleAllowedHole)
{
  // Mirrors the production ITS iteration shape exercised in
  // testTrackingTopology.cxx: MaxHoles=1, HoleLayerMask=1<<3.
  DetectorLayoutBuilder builder{denseCatalog(7)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1, 2, 3, 4, 5, 6}), 1, maskOf({3}), SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});
  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());

  TrackingTopology<7> reference;
  reference.init(7, 1, LayerMask{static_cast<uint16_t>(1 << 3)});

  checkParity<7>(result.layout->getTopology().getView(), reference.getView());
}

BOOST_AUTO_TEST_CASE(ExactParityWithTrackingTopologyTenNoHoles)
{
  DetectorLayoutBuilder builder{denseCatalog(10, SurfaceKind::Disk)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::DiskDisk});
  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());

  TrackingTopology<10> reference;
  reference.init(10, 0, LayerMask{0});

  checkParity<10>(result.layout->getTopology().getView(), reference.getView());
}

BOOST_AUTO_TEST_CASE(ExactParityWithTrackingTopologyTenSingleAllowedHole)
{
  // Mirrors the production MFT iteration shape exercised in
  // testTrackingTopology.cxx: MaxHoles=1, HoleLayerMask=1<<5.
  DetectorLayoutBuilder builder{denseCatalog(10, SurfaceKind::Disk)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), 1, maskOf({5}), SurfaceMask{}, TransitionPolicyTag::DiskDisk});
  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());

  TrackingTopology<10> reference;
  reference.init(10, 1, LayerMask{static_cast<uint16_t>(1 << 5)});

  checkParity<10>(result.layout->getTopology().getView(), reference.getView());
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

  DetectorLayoutBuilder builder{catalog};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1, 2, 3, 4, 5, 6}), 0, SurfaceMask{}, maskOf({0}), TransitionPolicyTag::CylinderCylinder});
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({7, 8, 9, 10, 11, 12, 13, 14, 15, 16}), 0, SurfaceMask{}, maskOf({7}), TransitionPolicyTag::DiskDisk});

  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(result.layout->valid());

  const auto masks = computeSurfaceKindMasks(catalog);
  const std::vector<NominalSurfaceMaterialBudget> zeroMaterial(catalog.size());
  const auto layoutView = result.layout->getView(catalog, masks.first, masks.second, zeroMaterial);
  BOOST_CHECK_EQUAL(layoutView.nSurfaces, 17u);
  BOOST_CHECK_EQUAL(layoutView.cylinderSurfaces.value(), 0x7fu);           // bits 0-6
  BOOST_CHECK_EQUAL(layoutView.diskSurfaces.value(), 0x1ff80u);            // bits 7-16

  const auto view = layoutView.topology;
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
  DetectorLayoutBuilder builder{denseCatalog(33)};
  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::TopologyRejected);
  BOOST_CHECK(result.topologyError == TopologyBuildError::InvalidSurfaceCount);
}

BOOST_AUTO_TEST_CASE(NonDenseCatalogSurfaceIdsAreRejected)
{
  std::vector<SurfaceDescriptor> catalog;
  catalog.push_back(surface(0));
  catalog.push_back(surface(2)); // should be 1 to stay dense
  DetectorLayoutBuilder builder{std::move(catalog)};

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::LayoutRejected);
  BOOST_CHECK(result.layoutError == DetectorLayoutError::NonDenseSurfaceIds);
}

BOOST_AUTO_TEST_CASE(InvalidSubgraphSurfaceIdsAreRejected)
{
  {
    // Out-of-range id.
    DetectorLayoutBuilder builder{denseCatalog(3)};
    builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 5}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});
    const auto result = builder.build();
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == DetectorLayoutBuildError::InvalidSubgraphSurfaceId);
  }
  {
    // The invalid-sentinel id.
    DetectorLayoutBuilder builder{denseCatalog(3)};
    std::vector<SurfaceId> ordered{SurfaceId{0}, SurfaceId::invalid()};
    builder.addSubgraph(DetectorLayoutSubgraph{std::move(ordered), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});
    const auto result = builder.build();
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == DetectorLayoutBuildError::InvalidSubgraphSurfaceId);
  }
}

BOOST_AUTO_TEST_CASE(DuplicateSurfaceWithinSubgraphIsRejected)
{
  DetectorLayoutBuilder builder{denseCatalog(3)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1, 0}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::DuplicateSurfaceInSubgraph);
}

BOOST_AUTO_TEST_CASE(SurfaceDuplicatedAcrossSubgraphsIsRejected)
{
  DetectorLayoutBuilder builder{denseCatalog(4)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({1, 2}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::SurfaceDuplicatedAcrossSubgraphs);
}

BOOST_AUTO_TEST_CASE(InvalidPolicyTagIsRejected)
{
  DetectorLayoutBuilder builder{denseCatalog(2)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::Invalid});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::TopologyRejected);
  BOOST_CHECK(result.topologyError == TopologyBuildError::InvalidPolicyTag);
}

BOOST_AUTO_TEST_CASE(PolicySurfaceKindMismatchIsRejected)
{
  // Both surfaces are disks, but the subgraph is tagged for cylinders.
  DetectorLayoutBuilder builder{denseCatalog(2, SurfaceKind::Disk)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::LayoutRejected);
  BOOST_CHECK(result.layoutError == DetectorLayoutError::PolicySurfaceKindMismatch);
}

BOOST_AUTO_TEST_CASE(SingletonSubgraphInvalidPolicyTagIsRejected)
{
  // A singleton subgraph never calls addTransition, so this must be caught
  // by explicit up-front validation rather than falling out of transition
  // enumeration.
  DetectorLayoutBuilder builder{denseCatalog(1)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::Invalid});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::TopologyRejected);
  BOOST_CHECK(result.topologyError == TopologyBuildError::InvalidPolicyTag);
}

BOOST_AUTO_TEST_CASE(SingletonDiskSurfaceTaggedCylinderCylinderIsRejected)
{
  DetectorLayoutBuilder builder{denseCatalog(1, SurfaceKind::Disk)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::LayoutRejected);
  BOOST_CHECK(result.layoutError == DetectorLayoutError::PolicySurfaceKindMismatch);
}

BOOST_AUTO_TEST_CASE(SingletonCylinderSurfaceTaggedDiskDiskIsRejected)
{
  DetectorLayoutBuilder builder{denseCatalog(1, SurfaceKind::Cylinder)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::DiskDisk});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::LayoutRejected);
  BOOST_CHECK(result.layoutError == DetectorLayoutError::PolicySurfaceKindMismatch);
}

BOOST_AUTO_TEST_CASE(NegativeMaxHolesIsRejected)
{
  // Explicit contract: a negative maxHoles is rejected outright rather than
  // silently normalized to zero.
  DetectorLayoutBuilder builder{denseCatalog(2)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1}), -1, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::NegativeMaxHoles);
}

BOOST_AUTO_TEST_CASE(HoleSurfacesOutsideSubgraphAreRejected)
{
  DetectorLayoutBuilder builder{denseCatalog(3)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1}), 1, maskOf({2}), SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::HoleSurfacesOutsideSubgraph);
}

BOOST_AUTO_TEST_CASE(SeedingSurfacesOutsideSubgraphAreRejected)
{
  DetectorLayoutBuilder builder{denseCatalog(3)};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0, 1}), 0, SurfaceMask{}, maskOf({2}), TransitionPolicyTag::CylinderCylinder});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::SeedingSurfacesOutsideSubgraph);
}

BOOST_AUTO_TEST_CASE(EmptyCatalogWithNoSubgraphsIsATrivialValidLayout)
{
  DetectorLayoutBuilder builder{{}};
  const auto result = builder.build();
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result.layout->getView({}, {}, {}, {}).nSurfaces, 0u);
  BOOST_CHECK_EQUAL(result.layout->getTopology().getView().nTransitions, 0u);
}

BOOST_AUTO_TEST_CASE(EmptySubgraphIsRejected)
{
  DetectorLayoutBuilder builder{denseCatalog(2)};
  builder.addSubgraph(DetectorLayoutSubgraph{{}, 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::EmptySubgraph);
}

BOOST_AUTO_TEST_CASE(EmptyCatalogWithNonEmptySubgraphIsRejected)
{
  DetectorLayoutBuilder builder{{}};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds({0}), 0, SurfaceMask{}, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});

  const auto result = builder.build();
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutBuildError::InvalidSubgraphSurfaceId);
}
