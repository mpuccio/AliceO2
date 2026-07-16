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

SurfaceTransition adjacent(uint16_t from, uint16_t to)
{
  return SurfaceTransition{SurfaceId{from}, SurfaceId{to}, SurfaceMask{}, 0, 0};
}
} // namespace

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

BOOST_AUTO_TEST_CASE(LayoutLimitsAndDenseIdsAreValidated)
{
  SparseTrackingTopology tooLarge{33};
  BOOST_CHECK(tooLarge.getError() == TopologyBuildError::InvalidSurfaceCount);

  SparseTrackingTopology topology{2};
  BOOST_REQUIRE(topology.finalize());
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.push_back(SurfaceDescriptor{SurfaceId{0}, 0, 0, SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{2}, 1, 0, SurfaceKind::Cylinder});
  DetectorLayout layout{std::move(surfaces), std::move(topology)};
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

BOOST_AUTO_TEST_CASE(DisconnectedCombinedTopologyIsSparse)
{
  SparseTrackingTopology topology{17};
  std::vector<TransitionId> transitions;

  for (uint16_t surface = 0; surface < 6; ++surface) {
    transitions.push_back(topology.addTransition(adjacent(surface, surface + 1)));
  }
  for (uint16_t surface = 7; surface < 16; ++surface) {
    transitions.push_back(topology.addTransition(adjacent(surface, surface + 1)));
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

BOOST_AUTO_TEST_CASE(LayoutCachesKindMasksInTheGlobalIdSpace)
{
  SparseTrackingTopology topology{4};
  BOOST_REQUIRE(topology.finalize());

  std::vector<SurfaceDescriptor> surfaces;
  surfaces.push_back(SurfaceDescriptor{SurfaceId{0}, 0, 0, SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{1}, 1, 0, SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{2}, 0, 8, SurfaceKind::Disk});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{3}, 1, 8, SurfaceKind::Disk});

  DetectorLayout layout{std::move(surfaces), std::move(topology)};
  BOOST_REQUIRE(layout.valid());
  const auto view = layout.getView();
  BOOST_CHECK_EQUAL(view.nSurfaces, 4u);
  BOOST_CHECK_EQUAL(view.cylinderSurfaces.value(), 0x3u);
  BOOST_CHECK_EQUAL(view.diskSurfaces.value(), 0xcu);
}
