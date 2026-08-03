// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B2 Slice 1: replacement combined-layout test proving the static,
// compile-time-concatenated ITS+MFT catalog (StaticDetectorCatalogs.h)
// reproduces the same one-global-id-space guarantees that the now-deleted
// testCombinedSurfaceCatalogBuilder.cxx's
// CombinedLayoutPreservesCurrentHoleSemanticsUnderGlobalIds established for
// the (also now-deleted) runtime CombinedSurfaceCatalogBuilder -- now from a
// static, borrowed source instead of a per-call runtime build. Gate 4 B2
// Slice 3 deleted CombinedSurfaceCatalogBuilder and its test once this
// replacement was already in place and passing (never independently
// beforehand): this is now the sole combined-catalog coverage.
//
// Proves, all at construction time, with no tracker/TimeFrame/workflow
// involved:
//  - dense global SurfaceIds 0..16 across the concatenated catalog;
//  - each surface's detector-qualified local identity (detectorId,
//    detectorSurfaceIndex, kind) is exactly what its own single-detector
//    static catalog declares -- unaffected by global rebasing;
//  - a SurfaceCatalogView constructed over the static storage borrows it
//    (pointer identity to the inline constexpr array's own .data()/element
//    addresses), proving no copy is required to view it;
//  - built under the shared global numbering, via the SAME (unmodified)
//    DetectorLayoutBuilder used today, each detector's own current
//    production-shaped hole/transition topology (MaxHoles=1,
//    HoleLayerMask=1<<3 for ITS; MaxHoles=1, HoleLayerMask=1<<5 for MFT --
//    the shapes testTrackingTopology.cxx documents as production and
//    testCombinedSurfaceCatalogBuilder.cxx already exercised under the
//    runtime-built catalog) is reproduced exactly against a standalone
//    TrackingTopology<NLayers> reference, using each detector's existing
//    ordered-surfaces convention (contiguous ids in traversal order);
//  - zero transitions or cells cross the ITS/MFT boundary;
//  - CSR (transition -> cell) consistency holds;
//  - no runtime component/subgraph concept is introduced: DetectorLayoutBuilder
//    and DetectorLayoutSubgraph are the same pre-existing types every
//    detector already uses standalone today; this test adds no new type,
//    only a new catalog *source*.
//
// Gate 4 B2 Slice 2 landed after this test: DetectorLayoutBuilder now
// borrows SurfaceCatalogView directly, so the topology test below borrows
// the static storage the same way the borrowing-proof test case already
// does -- no materialized vector anywhere in this file.

#define BOOST_TEST_MODULE ITSMFT CombinedStaticSurfaceCatalogTopology
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/SurfaceCatalogView.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITSMFTTracking/TrackingTopology.h"

namespace
{
using namespace o2::itsmft::tracking;

std::vector<SurfaceId> orderedIds(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> out;
  out.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    out.push_back(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return out;
}

void checkCsrConsistency(const SparseTrackingTopologyView& view)
{
  uint32_t total = 0;
  for (uint32_t t = 0; t < view.nTransitions; ++t) {
    const auto range = view.getCellsStartingWithTransition(TransitionId{static_cast<uint16_t>(t)});
    total += range.getEntries();
  }
  BOOST_CHECK_EQUAL(total, view.nCells);
}
} // namespace

BOOST_AUTO_TEST_CASE(CombinedStaticCatalogHasDenseGlobalIds)
{
  BOOST_REQUIRE_EQUAL(kITSMFTCombinedStaticSurfaceCatalog.size(), static_cast<std::size_t>(ITSNLayers + MFTNLayers));
  for (uint16_t i = 0; i < kITSMFTCombinedStaticSurfaceCatalog.size(); ++i) {
    BOOST_CHECK(kITSMFTCombinedStaticSurfaceCatalog[i].id == SurfaceId{i});
  }
}

BOOST_AUTO_TEST_CASE(CombinedStaticCatalogPreservesDetectorQualifiedLocalIdentity)
{
  for (uint16_t i = 0; i < ITSNLayers; ++i) {
    const auto& d = kITSMFTCombinedStaticSurfaceCatalog[i];
    BOOST_CHECK_EQUAL(d.detectorSurfaceIndex, i);
    BOOST_CHECK_EQUAL(d.detectorId, static_cast<uint8_t>(o2::detectors::DetID::ITS));
    BOOST_CHECK(d.kind == SurfaceKind::Cylinder);
    // Global rebasing touches only `id`; detector-local identity and
    // geometry must survive unchanged relative to the standalone
    // single-detector projection.
    BOOST_CHECK_EQUAL(d.detectorSurfaceIndex, kITSStaticSurfaceCatalog[i].detectorSurfaceIndex);
    BOOST_CHECK_EQUAL(d.detectorId, kITSStaticSurfaceCatalog[i].detectorId);
  }
  for (uint16_t i = 0; i < MFTNLayers; ++i) {
    const auto& d = kITSMFTCombinedStaticSurfaceCatalog[ITSNLayers + i];
    BOOST_CHECK_EQUAL(d.detectorSurfaceIndex, i);
    BOOST_CHECK_EQUAL(d.detectorId, static_cast<uint8_t>(o2::detectors::DetID::MFT));
    BOOST_CHECK(d.kind == SurfaceKind::Disk);
    BOOST_CHECK_EQUAL(d.detectorSurfaceIndex, kMFTStaticSurfaceCatalog[i].detectorSurfaceIndex);
    BOOST_CHECK_EQUAL(d.detectorId, kMFTStaticSurfaceCatalog[i].detectorId);
  }
}

BOOST_AUTO_TEST_CASE(SurfaceCatalogViewBorrowsTheStaticStorageWithoutCopying)
{
  const SurfaceCatalogView view{kITSMFTCombinedStaticSurfaceCatalog.data(),
                                static_cast<uint32_t>(kITSMFTCombinedStaticSurfaceCatalog.size())};
  // Pointer identity, not value equality: proves the view borrows the
  // inline constexpr array's own process-lifetime storage rather than a
  // copy of it.
  BOOST_CHECK(view.surfaces == kITSMFTCombinedStaticSurfaceCatalog.data());
  BOOST_CHECK_EQUAL(view.nSurfaces, kITSMFTCombinedStaticSurfaceCatalog.size());
  for (uint16_t i = 0; i < view.nSurfaces; ++i) {
    BOOST_CHECK(&view.getSurface(SurfaceId{i}) == &kITSMFTCombinedStaticSurfaceCatalog[i]);
  }
}

// The central regression this test exists to establish -- the static,
// compile-time-concatenated catalog's global-id offset, fed through the SAME
// (unmodified) DetectorLayoutBuilder used today, reproduces each detector's
// own current production-shaped hole/transition topology exactly, with zero
// transitions or cells crossing the ITS/MFT boundary.
BOOST_AUTO_TEST_CASE(CombinedStaticCatalogPreservesCurrentTopologyWithNoBoundaryCrossing)
{
  SurfaceMask itsHoleMask;
  itsHoleMask.set(SurfaceId{3});
  SurfaceMask mftHoleMask;
  mftHoleMask.set(SurfaceId{ITSNLayers + 5});

  // Gate 4 B2 Slice 2: DetectorLayoutBuilder now borrows SurfaceCatalogView
  // directly -- no materialized vector needed any more (the earlier Slice 1
  // version of this test built one as a deliberate, temporary exception;
  // that exception is gone now that the builder itself borrows).
  const SurfaceCatalogView catalogView{kITSMFTCombinedStaticSurfaceCatalog.data(),
                                       static_cast<uint32_t>(kITSMFTCombinedStaticSurfaceCatalog.size())};
  DetectorLayoutBuilder builder{catalogView};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds(0, ITSNLayers), 1, itsHoleMask, SurfaceMask{}});
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds(ITSNLayers, MFTNLayers), 1, mftHoleMask, SurfaceMask{}});
  const auto layoutResult = builder.build();
  BOOST_REQUIRE(layoutResult.ok());
  const auto view = layoutResult.layout->getTopology().getView();
  const auto masks = computeSurfaceKindMasks(kITSMFTCombinedStaticSurfaceCatalog);

  // Independent standalone references, exactly as production configures
  // each detector today (local 0-based ids), unrelated to the combined
  // catalog above -- the same shapes testTrackingTopology.cxx documents as
  // production (MaxHoles=1/HoleLayerMask=1<<3 for ITS, 1<<5 for MFT).
  TrackingTopology<ITSNLayers> itsReference;
  itsReference.init(ITSNLayers, 1, LayerMask{static_cast<uint16_t>(1 << 3)});
  TrackingTopology<MFTNLayers> mftReference;
  mftReference.init(MFTNLayers, 1, LayerMask{static_cast<uint16_t>(1 << 5)});

  BOOST_CHECK_EQUAL(view.nTransitions, static_cast<uint32_t>(itsReference.getView().nTransitions) +
                                         static_cast<uint32_t>(mftReference.getView().nTransitions));
  BOOST_CHECK_EQUAL(view.nCells, static_cast<uint32_t>(itsReference.getView().nCells) +
                                   static_cast<uint32_t>(mftReference.getView().nCells));

  // Every ITS transition/cell reproduces the standalone reference, offset by
  // nothing (ITS is first in the concatenation); every MFT one reproduces
  // its own standalone reference, offset by +ITSNLayers.
  const auto& itsRefView = itsReference.getView();
  for (uint32_t t = 0; t < itsRefView.nTransitions; ++t) {
    const auto& built = view.getTransition(TransitionId{static_cast<uint16_t>(t)});
    const auto& ref = itsRefView.getTransition(t);
    BOOST_CHECK_EQUAL(built.from.value(), ref.fromLayer);
    BOOST_CHECK_EQUAL(built.to.value(), ref.toLayer);
    BOOST_CHECK(masks.first.has(built.from) && masks.first.has(built.to));
  }
  const auto& mftRefView = mftReference.getView();
  for (uint32_t t = 0; t < mftRefView.nTransitions; ++t) {
    const auto& built = view.getTransition(TransitionId{static_cast<uint16_t>(itsRefView.nTransitions + t)});
    const auto& ref = mftRefView.getTransition(t);
    BOOST_CHECK_EQUAL(built.from.value(), ref.fromLayer + ITSNLayers);
    BOOST_CHECK_EQUAL(built.to.value(), ref.toLayer + ITSNLayers);
    BOOST_CHECK(masks.second.has(built.from) && masks.second.has(built.to));
  }

  // No declared edge crosses the ITS/MFT boundary: every transition and
  // every cell's hit-surface mask stays within one detector's global range.
  for (uint32_t t = 0; t < view.nTransitions; ++t) {
    const auto& transition = view.getTransition(TransitionId{static_cast<uint16_t>(t)});
    const bool bothCylinder = masks.first.has(transition.from) && masks.first.has(transition.to);
    const bool bothDisk = masks.second.has(transition.from) && masks.second.has(transition.to);
    BOOST_CHECK(bothCylinder || bothDisk);
  }
  for (uint32_t c = 0; c < view.nCells; ++c) {
    const auto& cell = view.getCell(CellTopologyId{static_cast<uint16_t>(c)});
    const bool allCylinder = cell.hitSurfaces.isSubsetOf(masks.first);
    const bool allDisk = cell.hitSurfaces.isSubsetOf(masks.second);
    BOOST_CHECK(allCylinder || allDisk);
  }

  checkCsrConsistency(view);
}
