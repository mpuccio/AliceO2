// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT CombinedSurfaceCatalogBuilder
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <vector>

#include "ITSMFTTracking/CombinedSurfaceCatalogBuilder.h"
#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITSMFTTracking/TrackingTopology.h"

namespace
{
using namespace o2::itsmft::tracking;

// Test-only stand-in for a real geometry-backed DetectorSurfaceCatalogProvider
// (ITSSurfaceCatalogProvider/MFTSurfaceCatalogProvider), so this file needs no
// initialized detector geometry singleton. It reproduces exactly the request
// contract every real provider enforces (Gate 2 GeometrySurfaceCatalogProvider
// behavior): detector must match, firstSurface must be 0, count must match.
class MockProvider final : public DetectorSurfaceCatalogProvider
{
 public:
  MockProvider(o2::detectors::DetID::ID detector, SurfaceKind kind, uint32_t count)
    : mDetector{detector}, mKind{kind}, mCount{count}
  {
  }

  void failNextCall(DetectorSurfaceCatalogError error) { mForcedError = error; }

  DetectorSurfaceCatalogResult buildCatalog(const DetectorSurfaceCatalogRequest& request) const final
  {
    if (mForcedError != DetectorSurfaceCatalogError::None) {
      return {{}, mForcedError};
    }
    if (request.detector != mDetector || request.firstSurface != SurfaceId{0} || request.detectorSurfaceCount != mCount) {
      return {{}, DetectorSurfaceCatalogError::InvalidRequest};
    }
    DetectorSurfaceCatalogResult result;
    result.catalog.reserve(mCount);
    for (uint16_t i = 0; i < mCount; ++i) {
      SurfaceDescriptor descriptor{SurfaceId{i}, i, static_cast<uint8_t>(mDetector), mKind, 0,
                                   static_cast<float>(10 + i), static_cast<float>(9 + i), static_cast<float>(11 + i)};
      descriptor.material.xOverX0 = 0.01f * (i + 1);
      descriptor.material.arealDensityGPerCm2 = 0.02f * (i + 1);
      result.catalog.push_back(descriptor);
    }
    return result;
  }

 private:
  o2::detectors::DetID::ID mDetector;
  SurfaceKind mKind;
  uint32_t mCount;
  DetectorSurfaceCatalogError mForcedError{DetectorSurfaceCatalogError::None};
};

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

BOOST_AUTO_TEST_CASE(CombinedCatalogIsDenseAndGloballyOffset)
{
  MockProvider itsProvider{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, ITSNLayers};
  MockProvider mftProvider{o2::detectors::DetID::MFT, SurfaceKind::Disk, MFTNLayers};

  const std::array<CombinedCatalogSlice, 2> slices{{{&itsProvider, o2::detectors::DetID::ITS, ITSNLayers},
                                                    {&mftProvider, o2::detectors::DetID::MFT, MFTNLayers}}};
  const auto result = buildCombinedSurfaceCatalog(slices);
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result.catalog.size(), ITSNLayers + MFTNLayers);

  // Dense global ids, contiguous per slice.
  for (uint16_t i = 0; i < result.catalog.size(); ++i) {
    BOOST_CHECK(result.catalog[i].id == SurfaceId{i});
  }

  // ITS slice: global ids [0,7), detector-local index unchanged, own detectorId/kind.
  for (uint16_t i = 0; i < ITSNLayers; ++i) {
    const auto& d = result.catalog[i];
    BOOST_CHECK_EQUAL(d.detectorSurfaceIndex, i);
    BOOST_CHECK_EQUAL(d.detectorId, static_cast<uint8_t>(o2::detectors::DetID::ITS));
    BOOST_CHECK(d.kind == SurfaceKind::Cylinder);
  }
  // MFT slice: global ids [7,17), detector-local index restarts at 0.
  for (uint16_t i = 0; i < MFTNLayers; ++i) {
    const auto& d = result.catalog[ITSNLayers + i];
    BOOST_CHECK_EQUAL(d.detectorSurfaceIndex, i);
    BOOST_CHECK_EQUAL(d.detectorId, static_cast<uint8_t>(o2::detectors::DetID::MFT));
    BOOST_CHECK(d.kind == SurfaceKind::Disk);
  }
}

BOOST_AUTO_TEST_CASE(NeitherProviderIsQueriedWithANonZeroOffset)
{
  // Precondition this whole design leans on: every provider still only ever
  // sees firstSurface=0, exactly as it is invoked standalone today. This is
  // enforced by MockProvider::buildCatalog() itself rejecting anything else;
  // a passing CombinedCatalogIsDenseAndGloballyOffset case above already
  // proves both slices succeeded, so both providers must have seen exactly
  // that request shape.
  MockProvider itsProvider{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, ITSNLayers};
  MockProvider mftProvider{o2::detectors::DetID::MFT, SurfaceKind::Disk, MFTNLayers};
  const std::array<CombinedCatalogSlice, 2> slices{{{&itsProvider, o2::detectors::DetID::ITS, ITSNLayers},
                                                    {&mftProvider, o2::detectors::DetID::MFT, MFTNLayers}}};
  BOOST_CHECK(buildCombinedSurfaceCatalog(slices).ok());
}

BOOST_AUTO_TEST_CASE(NullProviderSliceIsRejected)
{
  const std::array<CombinedCatalogSlice, 1> slices{{{nullptr, o2::detectors::DetID::ITS, ITSNLayers}}};
  const auto result = buildCombinedSurfaceCatalog(slices);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorSurfaceCatalogError::InvalidRequest);
  BOOST_CHECK_EQUAL(result.failedSlice, 0u);
  BOOST_CHECK(result.catalog.empty());
}

BOOST_AUTO_TEST_CASE(FailingProviderIsPropagatedTransactionally)
{
  MockProvider itsProvider{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, ITSNLayers};
  MockProvider mftProvider{o2::detectors::DetID::MFT, SurfaceKind::Disk, MFTNLayers};
  mftProvider.failNextCall(DetectorSurfaceCatalogError::GeometryNotInitialized);

  const std::array<CombinedCatalogSlice, 2> slices{{{&itsProvider, o2::detectors::DetID::ITS, ITSNLayers},
                                                    {&mftProvider, o2::detectors::DetID::MFT, MFTNLayers}}};
  const auto result = buildCombinedSurfaceCatalog(slices);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorSurfaceCatalogError::GeometryNotInitialized);
  BOOST_CHECK_EQUAL(result.failedSlice, 1u);
  // Transactional: the first (successful) slice's own success does not leak
  // into a partial result.
  BOOST_CHECK(result.catalog.empty());
}

BOOST_AUTO_TEST_CASE(SliceReturningWrongSizedCatalogIsRejected)
{
  // A provider that (incorrectly) returns a different count than declared.
  MockProvider wrongSized{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, ITSNLayers};
  const std::array<CombinedCatalogSlice, 1> slices{{{&wrongSized, o2::detectors::DetID::ITS, ITSNLayers - 1}}};
  const auto result = buildCombinedSurfaceCatalog(slices);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorSurfaceCatalogError::InvalidRequest);
  BOOST_CHECK_EQUAL(result.failedSlice, 0u);
}

BOOST_AUTO_TEST_CASE(SlicesExceedingCapacityAreRejectedBeforeQueryingFurtherProviders)
{
  MockProvider big1{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, 20};
  MockProvider big2{o2::detectors::DetID::MFT, SurfaceKind::Disk, 20};
  const std::array<CombinedCatalogSlice, 2> slices{{{&big1, o2::detectors::DetID::ITS, 20},
                                                    {&big2, o2::detectors::DetID::MFT, 20}}};
  const auto result = buildCombinedSurfaceCatalog(slices);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorSurfaceCatalogError::InvalidRequest);
  BOOST_CHECK_EQUAL(result.failedSlice, 1u); // 20 + 20 > 32, caught when the second slice is reached
}

// The central regression this slice exists to establish: the combined
// catalog's global-id offset, fed through the SAME (unmodified)
// DetectorLayoutBuilder used today, reproduces each detector's own current
// production hole/transition semantics exactly -- ITS's MaxHoles=1/
// HoleLayerMask=1<<3 and MFT's MaxHoles=1/HoleLayerMask=1<<5 (the shapes
// already exercised in testTrackingTopology.cxx and
// testDetectorLayoutBuilder.cxx) -- under the shared global numbering, with
// zero transitions or cells crossing the ITS/MFT boundary. No hole
// semantics, transition enumeration, or physics behavior changes; only the
// numeric id space each detector's edges are expressed in.
BOOST_AUTO_TEST_CASE(CombinedLayoutPreservesCurrentHoleSemanticsUnderGlobalIds)
{
  MockProvider itsProvider{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, ITSNLayers};
  MockProvider mftProvider{o2::detectors::DetID::MFT, SurfaceKind::Disk, MFTNLayers};
  const std::array<CombinedCatalogSlice, 2> slices{{{&itsProvider, o2::detectors::DetID::ITS, ITSNLayers},
                                                    {&mftProvider, o2::detectors::DetID::MFT, MFTNLayers}}};
  const auto combined = buildCombinedSurfaceCatalog(slices);
  BOOST_REQUIRE(combined.ok());
  BOOST_REQUIRE_EQUAL(combined.catalog.size(), ITSNLayers + MFTNLayers);

  SurfaceMask itsHoleMask;
  itsHoleMask.set(SurfaceId{3});
  SurfaceMask mftHoleMask;
  mftHoleMask.set(SurfaceId{ITSNLayers + 5});

  DetectorLayoutBuilder builder{SurfaceCatalogView{combined.catalog.data(), static_cast<uint32_t>(combined.catalog.size())}};
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds(0, ITSNLayers), 1, itsHoleMask, SurfaceMask{}, TransitionPolicyTag::CylinderCylinder});
  builder.addSubgraph(DetectorLayoutSubgraph{orderedIds(ITSNLayers, MFTNLayers), 1, mftHoleMask, SurfaceMask{}, TransitionPolicyTag::DiskDisk});
  const auto layoutResult = builder.build();
  BOOST_REQUIRE(layoutResult.ok());
  const auto view = layoutResult.layout->getTopology().getView();

  // Independent standalone references, exactly as production configures
  // each detector today (local 0-based ids), unrelated to the combined
  // catalog above.
  TrackingTopology<ITSNLayers> itsReference;
  itsReference.init(ITSNLayers, 1, LayerMask{static_cast<uint16_t>(1 << 3)});
  TrackingTopology<MFTNLayers> mftReference;
  mftReference.init(MFTNLayers, 1, LayerMask{static_cast<uint16_t>(1 << 5)});

  BOOST_CHECK_EQUAL(view.nTransitions, static_cast<uint32_t>(itsReference.getView().nTransitions) +
                                         static_cast<uint32_t>(mftReference.getView().nTransitions));
  BOOST_CHECK_EQUAL(view.nCells, static_cast<uint32_t>(itsReference.getView().nCells) +
                                   static_cast<uint32_t>(mftReference.getView().nCells));

  // Every ITS transition/cell reproduces the standalone reference, offset by
  // nothing (ITS is first in the allocation); every MFT one reproduces its
  // own standalone reference, offset by +ITSNLayers.
  const auto& itsRefView = itsReference.getView();
  for (uint32_t t = 0; t < itsRefView.nTransitions; ++t) {
    const auto& built = view.getTransition(TransitionId{static_cast<uint16_t>(t)});
    const auto& ref = itsRefView.getTransition(t);
    BOOST_CHECK_EQUAL(built.from.value(), ref.fromLayer);
    BOOST_CHECK_EQUAL(built.to.value(), ref.toLayer);
    BOOST_CHECK(built.policyTag == TransitionPolicyTag::CylinderCylinder);
  }
  const auto& mftRefView = mftReference.getView();
  for (uint32_t t = 0; t < mftRefView.nTransitions; ++t) {
    const auto& built = view.getTransition(TransitionId{static_cast<uint16_t>(itsRefView.nTransitions + t)});
    const auto& ref = mftRefView.getTransition(t);
    BOOST_CHECK_EQUAL(built.from.value(), ref.fromLayer + ITSNLayers);
    BOOST_CHECK_EQUAL(built.to.value(), ref.toLayer + ITSNLayers);
    BOOST_CHECK(built.policyTag == TransitionPolicyTag::DiskDisk);
  }

  // No declared edge crosses the ITS/MFT boundary: every transition and
  // every cell's hit-surface mask stays within one detector's global range.
  const auto masks = computeSurfaceKindMasks(combined.catalog);
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
