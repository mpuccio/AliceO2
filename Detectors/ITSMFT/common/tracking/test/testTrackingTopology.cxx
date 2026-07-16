// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE ITSMFT TrackingTopology
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "ITSMFTTracking/TrackingTopology.h"

/// Baseline characterization of o2::itsmft::tracking::TrackingTopology (Wave 0).
///
/// NLayers=7 mirrors the current ITS barrel configuration and NLayers=10
/// mirrors the current MFT disk configuration (see DetectorTraits.cxx /
/// Configuration.h: detIdFromNLayers<7>() -> ITS, detIdFromNLayers<10>() -> MFT).
/// These tests record what the shared topology currently builds; they do not
/// encode the future explicit sparse graph/surface-mask design described in
/// Architecture.md.

using o2::itsmft::tracking::LayerMask;
using o2::itsmft::tracking::TrackingTopology;

namespace
{
template <int NLayers>
void checkCellsByFirstTransitionConsistency(const typename TrackingTopology<NLayers>::View& view)
{
  int total = 0;
  for (int t = 0; t < view.nTransitions; ++t) {
    const auto range = view.getCellsStartingWithTransition(t);
    total += range.getEntries();
    for (int k = 0; k < range.getEntries(); ++k) {
      const auto cellId = view.cellsByFirstTransition[range.getFirstEntry() + k];
      BOOST_CHECK_EQUAL(view.getCell(cellId).firstTransition, t);
    }
  }
  BOOST_CHECK_EQUAL(total, view.nCells);
}
} // namespace

BOOST_AUTO_TEST_CASE(static_capacity_formulas)
{
  // NLayers=7 (ITS-like)
  BOOST_CHECK_EQUAL(TrackingTopology<7>::MaxTransitions, 21); // 7*6/2
  BOOST_CHECK_EQUAL(TrackingTopology<7>::MaxCells, 35);       // 7*6*5/6

  // NLayers=10 (MFT-like)
  BOOST_CHECK_EQUAL(TrackingTopology<10>::MaxTransitions, 45); // 10*9/2
  BOOST_CHECK_EQUAL(TrackingTopology<10>::MaxCells, 120);      // 10*9*8/6
}

BOOST_AUTO_TEST_CASE(its_like_topology_no_holes)
{
  TrackingTopology<7> topo;
  topo.init(7, 0, LayerMask{0});
  const auto view = topo.getView();
  view.print();

  // With no holes allowed, only adjacent-layer transitions/cells survive.
  BOOST_CHECK_EQUAL(view.nTransitions, 6); // NLayers - 1
  BOOST_CHECK_EQUAL(view.nCells, 5);       // NLayers - 2

  for (int i = 0; i < view.nTransitions; ++i) {
    const auto& t = view.getTransition(i);
    BOOST_CHECK_EQUAL(t.fromLayer, i);
    BOOST_CHECK_EQUAL(t.toLayer, i + 1);
  }
  for (int i = 0; i < view.nCells; ++i) {
    const auto& cell = view.getCell(i);
    BOOST_CHECK_EQUAL(cell.firstTransition, i);
    BOOST_CHECK_EQUAL(cell.secondTransition, i + 1);
    BOOST_CHECK(cell.hitLayerMask.holeMask().empty());
    BOOST_CHECK_EQUAL(cell.hitLayerMask.count(), 3);
  }

  checkCellsByFirstTransitionConsistency<7>(view);
}

BOOST_AUTO_TEST_CASE(mft_like_topology_no_holes)
{
  TrackingTopology<10> topo;
  topo.init(10, 0, LayerMask{0});
  const auto view = topo.getView();
  view.print();

  BOOST_CHECK_EQUAL(view.nTransitions, 9); // NLayers - 1
  BOOST_CHECK_EQUAL(view.nCells, 8);       // NLayers - 2

  for (int i = 0; i < view.nTransitions; ++i) {
    const auto& t = view.getTransition(i);
    BOOST_CHECK_EQUAL(t.fromLayer, i);
    BOOST_CHECK_EQUAL(t.toLayer, i + 1);
  }
  for (int i = 0; i < view.nCells; ++i) {
    const auto& cell = view.getCell(i);
    BOOST_CHECK(cell.hitLayerMask.holeMask().empty());
    BOOST_CHECK_EQUAL(cell.hitLayerMask.count(), 3);
  }

  checkCellsByFirstTransitionConsistency<10>(view);
}

BOOST_AUTO_TEST_CASE(its_like_topology_restricts_maxlayers_below_nlayers)
{
  // TrackingTopology<7> used with fewer active layers than the compile-time
  // maximum (e.g. a reduced/partial-layout iteration).
  TrackingTopology<7> topo;
  topo.init(4, 0, LayerMask{0});
  const auto view = topo.getView();

  BOOST_CHECK_EQUAL(view.nTransitions, 3);
  BOOST_CHECK_EQUAL(view.nCells, 2);
  for (int i = 0; i < view.nTransitions; ++i) {
    BOOST_CHECK_LT(view.getTransition(i).toLayer, 4);
  }
}

BOOST_AUTO_TEST_CASE(its_like_topology_single_allowed_hole)
{
  // maxHoles=1, only layer 3 may be skipped: reproduces the shape of a
  // production ITS iteration configured with MaxHoles=1 / HoleLayerMask=1<<3.
  TrackingTopology<7> topo;
  topo.init(7, 1, LayerMask{static_cast<uint16_t>(1 << 3)});
  const auto view = topo.getView();
  view.print();

  // Hand-enumerated: 6 adjacent transitions + 1 skip transition (2->4).
  BOOST_CHECK_EQUAL(view.nTransitions, 7);
  // 5 contiguous cells + 2 hole-containing cells (1,2,4) and (2,4,5).
  BOOST_CHECK_EQUAL(view.nCells, 7);

  bool hasHoleTransition = false;
  for (int i = 0; i < view.nTransitions; ++i) {
    const auto& t = view.getTransition(i);
    hasHoleTransition |= (t.fromLayer == 2 && t.toLayer == 4);
    BOOST_CHECK(LayerMask::skipped(t.fromLayer, t.toLayer).isAllowedHoleMask(1, 1 << 3));
  }
  BOOST_CHECK(hasHoleTransition);

  int holeCells = 0;
  for (int i = 0; i < view.nCells; ++i) {
    const auto& cell = view.getCell(i);
    BOOST_CHECK(cell.hitLayerMask.isAllowed(1, 1 << 3));
    if (!cell.hitLayerMask.holeMask().empty()) {
      ++holeCells;
      BOOST_CHECK_EQUAL(cell.hitLayerMask.holeMask().value(), 1 << 3);
    }
  }
  BOOST_CHECK_EQUAL(holeCells, 2);

  checkCellsByFirstTransitionConsistency<7>(view);
}

BOOST_AUTO_TEST_CASE(its_like_topology_rejects_wrong_hole_layer)
{
  TrackingTopology<7> topo;
  topo.init(7, 1, LayerMask{static_cast<uint16_t>(1 << 3)});
  const auto view = topo.getView();

  for (int i = 0; i < view.nTransitions; ++i) {
    const auto& t = view.getTransition(i);
    BOOST_CHECK(!(t.fromLayer == 0 && t.toLayer == 2));
    BOOST_CHECK(!(t.fromLayer == 1 && t.toLayer == 3));
    BOOST_CHECK(!(t.fromLayer == 3 && t.toLayer == 5));
    BOOST_CHECK(!(t.fromLayer == 4 && t.toLayer == 6));
  }
  for (int i = 0; i < view.nCells; ++i) {
    BOOST_CHECK(view.getCell(i).hitLayerMask.holeMask().isSubsetOf(1 << 3));
  }
}

BOOST_AUTO_TEST_CASE(mft_like_topology_single_allowed_hole)
{
  // maxHoles=1, only layer 5 may be skipped (mid-stack MFT-like hole).
  TrackingTopology<10> topo;
  topo.init(10, 1, LayerMask{static_cast<uint16_t>(1 << 5)});
  const auto view = topo.getView();
  view.print();

  // Same combinatorial pattern as the 7-layer case, shifted to a 10-layer
  // stack: 9 adjacent transitions + 1 skip transition (4->6);
  // 8 contiguous cells + 2 hole-containing cells (3,4,6) and (4,6,7).
  BOOST_CHECK_EQUAL(view.nTransitions, 10);
  BOOST_CHECK_EQUAL(view.nCells, 10);

  bool hasHoleTransition = false;
  for (int i = 0; i < view.nTransitions; ++i) {
    const auto& t = view.getTransition(i);
    hasHoleTransition |= (t.fromLayer == 4 && t.toLayer == 6);
    BOOST_CHECK(LayerMask::skipped(t.fromLayer, t.toLayer).isAllowedHoleMask(1, 1 << 5));
  }
  BOOST_CHECK(hasHoleTransition);

  int holeCells = 0;
  for (int i = 0; i < view.nCells; ++i) {
    const auto& cell = view.getCell(i);
    BOOST_CHECK(cell.hitLayerMask.isAllowed(1, 1 << 5));
    if (!cell.hitLayerMask.holeMask().empty()) {
      ++holeCells;
      BOOST_CHECK_EQUAL(cell.hitLayerMask.holeMask().value(), 1 << 5);
    }
  }
  BOOST_CHECK_EQUAL(holeCells, 2);

  checkCellsByFirstTransitionConsistency<10>(view);
}

BOOST_AUTO_TEST_CASE(topology_reinit_clears_previous_state)
{
  // init() must be idempotent/re-entrant: TimeFrame re-initializes the
  // default topology per-iteration (see TimeFrame::initDefaultTrackingTopology).
  TrackingTopology<7> topo;
  topo.init(7, 1, LayerMask{static_cast<uint16_t>(1 << 3)});
  BOOST_CHECK_EQUAL(topo.getView().nTransitions, 7);

  topo.init(7, 0, LayerMask{0});
  const auto view = topo.getView();
  BOOST_CHECK_EQUAL(view.nTransitions, 6);
  BOOST_CHECK_EQUAL(view.nCells, 5);
  for (int i = 0; i < view.nTransitions; ++i) {
    BOOST_CHECK_EQUAL(view.getTransition(i).toLayer, view.getTransition(i).fromLayer + 1);
  }
}
