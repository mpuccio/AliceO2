// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT SurfaceLayout
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/DetectorLayout.h"

namespace
{
using namespace o2::itsmft::tracking;

SurfaceTransition adjacent(uint16_t from, uint16_t to, TransitionPolicyTag policyTag = TransitionPolicyTag::CylinderCylinder)
{
  return SurfaceTransition{SurfaceId{from}, SurfaceId{to}, SurfaceMask{}, policyTag, 0};
}

SurfaceDescriptor surface(uint16_t id, SurfaceKind kind)
{
  return SurfaceDescriptor{SurfaceId{id}, id, 0, kind};
}
} // namespace

BOOST_AUTO_TEST_CASE(PolicyTagsMapToStateFamilies)
{
  BOOST_CHECK(stateFamilyOf(TransitionPolicyTag::Invalid) == StateFamily::Invalid);
  BOOST_CHECK(stateFamilyOf(TransitionPolicyTag::CylinderCylinder) == StateFamily::Barrel);
  BOOST_CHECK(stateFamilyOf(TransitionPolicyTag::DiskDisk) == StateFamily::Forward);
  BOOST_CHECK(stateFamilyOf(static_cast<TransitionPolicyTag>(3)) == StateFamily::Invalid);

  BOOST_CHECK(isKnownTransitionPolicyTag(TransitionPolicyTag::Invalid));
  BOOST_CHECK(isKnownTransitionPolicyTag(TransitionPolicyTag::CylinderCylinder));
  BOOST_CHECK(isKnownTransitionPolicyTag(TransitionPolicyTag::DiskDisk));
  BOOST_CHECK(!isKnownTransitionPolicyTag(static_cast<TransitionPolicyTag>(3)));
  BOOST_CHECK(!isStageATransitionPolicyTagEnabled(TransitionPolicyTag::Invalid));
  BOOST_CHECK(isStageATransitionPolicyTagEnabled(TransitionPolicyTag::CylinderCylinder));
  BOOST_CHECK(isStageATransitionPolicyTagEnabled(TransitionPolicyTag::DiskDisk));
}

BOOST_AUTO_TEST_CASE(PolicyTagAbiIsStable)
{
  BOOST_CHECK_EQUAL(sizeof(StateFamily), sizeof(uint8_t));
  BOOST_CHECK_EQUAL(sizeof(TransitionPolicyTag), sizeof(uint16_t));
  BOOST_CHECK_EQUAL(sizeof(SurfaceTransition), 12u);
  BOOST_CHECK_EQUAL(sizeof(SurfaceCellTopology), 8u);
}

BOOST_AUTO_TEST_CASE(SurfaceMaskCoversThirtyTwoGlobalSurfaces)
{
  SurfaceMask mask;
  mask.set(SurfaceId{0});
  mask.set(SurfaceId{16});
  mask.set(SurfaceId{31});

  BOOST_CHECK(mask.has(SurfaceId{0}));
  BOOST_CHECK(mask.has(SurfaceId{16}));
  BOOST_CHECK(mask.has(SurfaceId{31}));
  BOOST_CHECK_EQUAL(mask.count(), 3);
  BOOST_CHECK_EQUAL(mask.first(), 0);
  BOOST_CHECK_EQUAL(mask.last(), 31);

  mask.reset(SurfaceId{16});
  BOOST_CHECK(!mask.has(SurfaceId{16}));
  BOOST_CHECK_EQUAL(mask.value(), (uint32_t{1} | (uint32_t{1} << 31)));
}

BOOST_AUTO_TEST_CASE(ParallelTransitionsAreRejected)
{
  SparseTrackingTopology topology{3};
  BOOST_CHECK(topology.addTransition(adjacent(0, 1)).isValid());
  BOOST_CHECK(!topology.addTransition(adjacent(0, 1)).isValid());
  BOOST_CHECK(topology.getError() == TopologyBuildError::DuplicateTransition);
  BOOST_CHECK(!topology.finalize());
}

BOOST_AUTO_TEST_CASE(InvalidAndUnknownPolicyTagsAreRejected)
{
  SparseTrackingTopology invalid{2};
  BOOST_CHECK(!invalid.addTransition(adjacent(0, 1, TransitionPolicyTag::Invalid)).isValid());
  BOOST_CHECK(invalid.getError() == TopologyBuildError::InvalidPolicyTag);

  SparseTrackingTopology unknown{2};
  BOOST_CHECK(!unknown.addTransition(adjacent(0, 1, static_cast<TransitionPolicyTag>(0xffff))).isValid());
  BOOST_CHECK(unknown.getError() == TopologyBuildError::InvalidPolicyTag);
}

BOOST_AUTO_TEST_CASE(LayoutLimitsAndDenseIdsAreValidated)
{
  SparseTrackingTopology tooLarge{33};
  BOOST_CHECK(tooLarge.getError() == TopologyBuildError::InvalidSurfaceCount);

  SparseTrackingTopology topology{2};
  BOOST_REQUIRE(topology.finalize());
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.push_back(SurfaceDescriptor{SurfaceId{0}, 0, 0, SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{2}, 1, 0, SurfaceKind::Cylinder});
  DetectorLayout layout{surfaces, std::move(topology)};
  BOOST_CHECK(!layout.valid());
  BOOST_CHECK(layout.getError() == DetectorLayoutError::NonDenseSurfaceIds);
}

BOOST_AUTO_TEST_CASE(CellCannotReturnToItsFirstSurface)
{
  SparseTrackingTopology topology{2};
  const auto outward = topology.addTransition(adjacent(0, 1));
  const auto returning = topology.addTransition(adjacent(1, 0));
  BOOST_REQUIRE(outward.isValid());
  BOOST_REQUIRE(returning.isValid());

  BOOST_CHECK(!topology.addCell(outward, returning).isValid());
  BOOST_CHECK(topology.getError() == TopologyBuildError::RepeatedSurface);
}

BOOST_AUTO_TEST_CASE(CellCannotMixPolicies)
{
  SparseTrackingTopology topology{3};
  const auto barrel = topology.addTransition(adjacent(0, 1, TransitionPolicyTag::CylinderCylinder));
  const auto forward = topology.addTransition(adjacent(1, 2, TransitionPolicyTag::DiskDisk));
  BOOST_REQUIRE(barrel.isValid());
  BOOST_REQUIRE(forward.isValid());

  BOOST_CHECK(!topology.addCell(barrel, forward).isValid());
  BOOST_CHECK(topology.getError() == TopologyBuildError::MixedPolicyCell);
}

BOOST_AUTO_TEST_CASE(DisconnectedCombinedTopologyIsSparse)
{
  SparseTrackingTopology topology{17};
  std::vector<TransitionId> transitions;

  for (uint16_t surface = 0; surface < 6; ++surface) {
    transitions.push_back(topology.addTransition(adjacent(surface, surface + 1)));
  }
  for (uint16_t surface = 7; surface < 16; ++surface) {
    transitions.push_back(topology.addTransition(adjacent(surface, surface + 1, TransitionPolicyTag::DiskDisk)));
  }
  BOOST_REQUIRE_EQUAL(transitions.size(), 15u);

  for (uint16_t transition = 0; transition < 5; ++transition) {
    BOOST_REQUIRE(topology.addCell(transitions[transition], transitions[transition + 1]).isValid());
  }
  for (uint16_t transition = 6; transition < 14; ++transition) {
    BOOST_REQUIRE(topology.addCell(transitions[transition], transitions[transition + 1]).isValid());
  }

  BOOST_REQUIRE(topology.finalize());
  const auto view = topology.getView();
  BOOST_CHECK_EQUAL(view.nTransitions, 15u);
  BOOST_CHECK_EQUAL(view.nCells, 13u);
  BOOST_CHECK_EQUAL(view.getCell(CellTopologyId{0}).hitSurfaces.count(), 3);

  const auto firstSuccessors = view.getCellsStartingWithTransition(transitions[0]);
  BOOST_CHECK_EQUAL(firstSuccessors.getFirstEntry(), 0u);
  BOOST_CHECK_EQUAL(firstSuccessors.getEntries(), 1u);
  BOOST_CHECK(view.cellsByFirstTransition[firstSuccessors.getFirstEntry()] == CellTopologyId{0});

  // The final ITS transition has no successor and is not connected to MFT.
  BOOST_CHECK_EQUAL(view.getCellsStartingWithTransition(transitions[5]).getEntries(), 0u);
}

BOOST_AUTO_TEST_CASE(CylinderAndDiskLayoutsAcceptMatchingPolicies)
{
  SparseTrackingTopology cylinderTopology{2};
  BOOST_REQUIRE(cylinderTopology.addTransition(adjacent(0, 1, TransitionPolicyTag::CylinderCylinder)).isValid());
  BOOST_REQUIRE(cylinderTopology.finalize());
  const std::vector<SurfaceDescriptor> cylinderSurfaces{surface(0, SurfaceKind::Cylinder), surface(1, SurfaceKind::Cylinder)};
  DetectorLayout cylinderLayout{cylinderSurfaces, std::move(cylinderTopology)};
  BOOST_CHECK(cylinderLayout.valid());

  SparseTrackingTopology diskTopology{2};
  BOOST_REQUIRE(diskTopology.addTransition(adjacent(0, 1, TransitionPolicyTag::DiskDisk)).isValid());
  BOOST_REQUIRE(diskTopology.finalize());
  const std::vector<SurfaceDescriptor> diskSurfaces{surface(0, SurfaceKind::Disk), surface(1, SurfaceKind::Disk)};
  DetectorLayout diskLayout{diskSurfaces, std::move(diskTopology)};
  BOOST_CHECK(diskLayout.valid());
}

BOOST_AUTO_TEST_CASE(DisconnectedCombinedLayoutAcceptsBothPolicies)
{
  SparseTrackingTopology topology{4};
  BOOST_REQUIRE(topology.addTransition(adjacent(0, 1, TransitionPolicyTag::CylinderCylinder)).isValid());
  BOOST_REQUIRE(topology.addTransition(adjacent(2, 3, TransitionPolicyTag::DiskDisk)).isValid());
  BOOST_REQUIRE(topology.finalize());

  const std::vector<SurfaceDescriptor> surfaces{surface(0, SurfaceKind::Cylinder), surface(1, SurfaceKind::Cylinder),
                                                surface(2, SurfaceKind::Disk), surface(3, SurfaceKind::Disk)};
  DetectorLayout layout{surfaces, std::move(topology)};
  BOOST_CHECK(layout.valid());
}

BOOST_AUTO_TEST_CASE(LayoutRejectsPolicySurfaceKindMismatch)
{
  SparseTrackingTopology diskPolicy{2};
  BOOST_REQUIRE(diskPolicy.addTransition(adjacent(0, 1, TransitionPolicyTag::DiskDisk)).isValid());
  BOOST_REQUIRE(diskPolicy.finalize());

  const std::vector<SurfaceDescriptor> cylinderSurfaces{surface(0, SurfaceKind::Cylinder), surface(1, SurfaceKind::Cylinder)};
  DetectorLayout cylinders{cylinderSurfaces, std::move(diskPolicy)};
  BOOST_CHECK(!cylinders.valid());
  BOOST_CHECK(cylinders.getError() == DetectorLayoutError::PolicySurfaceKindMismatch);

  SparseTrackingTopology cylinderPolicy{2};
  BOOST_REQUIRE(cylinderPolicy.addTransition(adjacent(0, 1, TransitionPolicyTag::CylinderCylinder)).isValid());
  BOOST_REQUIRE(cylinderPolicy.finalize());

  const std::vector<SurfaceDescriptor> diskSurfaces{surface(0, SurfaceKind::Disk), surface(1, SurfaceKind::Disk)};
  DetectorLayout disks{diskSurfaces, std::move(cylinderPolicy)};
  BOOST_CHECK(!disks.valid());
  BOOST_CHECK(disks.getError() == DetectorLayoutError::PolicySurfaceKindMismatch);
}

BOOST_AUTO_TEST_CASE(LayoutRejectsCylinderDiskTransitions)
{
  SparseTrackingTopology cylinderToDisk{2};
  BOOST_REQUIRE(cylinderToDisk.addTransition(adjacent(0, 1, TransitionPolicyTag::CylinderCylinder)).isValid());
  BOOST_REQUIRE(cylinderToDisk.finalize());
  const std::vector<SurfaceDescriptor> outwardSurfaces{surface(0, SurfaceKind::Cylinder), surface(1, SurfaceKind::Disk)};
  DetectorLayout outward{outwardSurfaces, std::move(cylinderToDisk)};
  BOOST_CHECK(!outward.valid());
  BOOST_CHECK(outward.getError() == DetectorLayoutError::MixedSurfaceTransition);

  SparseTrackingTopology diskToCylinder{2};
  BOOST_REQUIRE(diskToCylinder.addTransition(adjacent(0, 1, TransitionPolicyTag::DiskDisk)).isValid());
  BOOST_REQUIRE(diskToCylinder.finalize());
  const std::vector<SurfaceDescriptor> inwardSurfaces{surface(0, SurfaceKind::Disk), surface(1, SurfaceKind::Cylinder)};
  DetectorLayout inward{inwardSurfaces, std::move(diskToCylinder)};
  BOOST_CHECK(!inward.valid());
  BOOST_CHECK(inward.getError() == DetectorLayoutError::MixedSurfaceTransition);
}

BOOST_AUTO_TEST_CASE(LayoutCachesKindMasksInTheGlobalIdSpace)
{
  SparseTrackingTopology topology{4};
  BOOST_REQUIRE(topology.finalize());

  std::vector<SurfaceDescriptor> surfaces;
  surfaces.push_back(SurfaceDescriptor{SurfaceId{0}, 0, 0, SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{1}, 1, 0, SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{2}, 0, 8, SurfaceKind::Disk});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{3}, 1, 8, SurfaceKind::Disk});

  // Cylinder/disk masks are now a full-catalogue property computed once by
  // the shared owner (DetectorLayoutSet in production; computeSurfaceKindMasks()
  // directly here, since this test builds one DetectorLayout in isolation),
  // not cached per DetectorLayout.
  const auto masks = computeSurfaceKindMasks(surfaces);
  DetectorLayout layout{surfaces, std::move(topology)};
  BOOST_REQUIRE(layout.valid());
  const auto view = layout.getView(surfaces, masks.first, masks.second);
  BOOST_CHECK_EQUAL(view.nSurfaces, 4u);
  BOOST_CHECK_EQUAL(view.cylinderSurfaces.value(), 0x3u);
  BOOST_CHECK_EQUAL(view.diskSurfaces.value(), 0xcu);
}
