// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#define BOOST_TEST_MODULE ITSMFT DetectorTraversalBinding
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utility>
#include <vector>

#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/DetectorTraversalBinding.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/TrackingTopology.h"

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

struct CombinedLayout {
  DetectorLayout layout;
  DetectorLayoutView view;

  CombinedLayout()
    : layout{build()}, view{layout.getView(kITSMFTCombinedStaticSurfaceCatalog, cylinderMask(), diskMask())}
  {
  }

 private:
  static std::pair<SurfaceMask, SurfaceMask> masks()
  {
    return computeSurfaceKindMasks(kITSMFTCombinedStaticSurfaceCatalog);
  }
  static SurfaceMask cylinderMask() { return masks().first; }
  static SurfaceMask diskMask() { return masks().second; }
  static DetectorLayout build()
  {
    SurfaceMask itsHoles;
    itsHoles.set(SurfaceId{3});
    SurfaceMask mftHoles;
    mftHoles.set(SurfaceId{12});
    SurfaceMask seeds;
    seeds.set(SurfaceId{6});
    seeds.set(SurfaceId{16});
    const SurfaceCatalogView catalog{kITSMFTCombinedStaticSurfaceCatalog.data(), static_cast<uint32_t>(kITSMFTCombinedStaticSurfaceCatalog.size())};
    DetectorLayoutBuilder builder{catalog};
    builder.addSubgraph({ordered(0, ITSNLayers), 1, itsHoles, seeds & SurfaceMask{uint32_t{0x7f}}, TransitionPolicyTag::CylinderCylinder});
    builder.addSubgraph({ordered(ITSNLayers, MFTNLayers), 1, mftHoles, seeds & ~SurfaceMask{uint32_t{0x7f}}, TransitionPolicyTag::DiskDisk});
    auto built = builder.build();
    BOOST_REQUIRE(built.ok());
    return std::move(*built.layout);
  }
};

SurfaceMask itsMask() { return SurfaceMask{uint32_t{0x7f}}; }
SurfaceMask mftMask() { return SurfaceMask{uint32_t{0x1ff80}}; }

template <int NLayers>
void checkParity(const DetectorTraversalBinding& binding, const DetectorLayoutView& global, uint16_t surfaceOffset, uint16_t holeLayer)
{
  TrackingTopology<NLayers> local;
  local.init(NLayers, 1, LayerMask{static_cast<uint16_t>(uint16_t{1} << holeLayer)});
  const auto localView = local.getView();
  BOOST_REQUIRE_EQUAL(binding.getGlobalTransitions().size(), localView.nTransitions);
  BOOST_REQUIRE_EQUAL(binding.getGlobalCells().size(), localView.nCells);
  for (uint32_t slot = 0; slot < localView.nTransitions; ++slot) {
    const auto globalId = binding.getGlobalTransitions()[slot];
    BOOST_REQUIRE(binding.getScratchTransitionSlot(globalId));
    BOOST_CHECK_EQUAL(*binding.getScratchTransitionSlot(globalId), slot);
    const auto& transition = global.topology.getTransition(globalId);
    const auto& expected = localView.getTransition(slot);
    BOOST_CHECK_EQUAL(*binding.getLegacyLayer(transition.from), expected.fromLayer);
    BOOST_CHECK_EQUAL(*binding.getLegacyLayer(transition.to), expected.toLayer);
    BOOST_CHECK_EQUAL(transition.from.value(), expected.fromLayer + surfaceOffset);
    BOOST_CHECK_EQUAL(transition.to.value(), expected.toLayer + surfaceOffset);
  }
  for (uint32_t slot = 0; slot < localView.nCells; ++slot) {
    const auto globalId = binding.getGlobalCells()[slot];
    BOOST_REQUIRE(binding.getScratchCellSlot(globalId));
    BOOST_CHECK_EQUAL(*binding.getScratchCellSlot(globalId), slot);
    const auto& cell = global.topology.getCell(globalId);
    const auto& expected = localView.getCell(slot);
    BOOST_CHECK_EQUAL(*binding.getScratchTransitionSlot(cell.firstTransition), expected.firstTransition);
    BOOST_CHECK_EQUAL(*binding.getScratchTransitionSlot(cell.secondTransition), expected.secondTransition);
  }
}
} // namespace

BOOST_AUTO_TEST_CASE(GlobalBindingsAreCompleteDisjointAndRetainGlobalIds)
{
  CombinedLayout combined;
  const auto its = DetectorTraversalBinding::build(combined.view, o2::detectors::DetID::ITS, ClusterSourceId{0}, itsMask(), ordered(0, ITSNLayers));
  const auto mft = DetectorTraversalBinding::build(combined.view, o2::detectors::DetID::MFT, ClusterSourceId{1}, mftMask(), ordered(ITSNLayers, MFTNLayers));
  BOOST_REQUIRE(its.ok());
  BOOST_REQUIRE(mft.ok());
  const auto& itsBinding = *its.binding;
  const auto& mftBinding = *mft.binding;

  BOOST_CHECK_EQUAL(itsBinding.getSource().value(), 0);
  BOOST_CHECK_EQUAL(mftBinding.getSource().value(), 1);
  for (uint16_t surface = 0; surface < ITSNLayers; ++surface) {
    BOOST_REQUIRE(itsBinding.getLegacyLayer(SurfaceId{surface}));
    BOOST_CHECK_EQUAL(*itsBinding.getLegacyLayer(SurfaceId{surface}), surface);
    BOOST_CHECK(!mftBinding.getLegacyLayer(SurfaceId{surface}));
  }
  for (uint16_t layer = 0; layer < MFTNLayers; ++layer) {
    const SurfaceId surface{static_cast<uint16_t>(ITSNLayers + layer)};
    BOOST_REQUIRE(mftBinding.getLegacyLayer(surface));
    BOOST_CHECK_EQUAL(*mftBinding.getLegacyLayer(surface), layer);
    BOOST_CHECK(!itsBinding.getLegacyLayer(surface));
  }
  BOOST_CHECK(!mftBinding.getLegacyLayer(SurfaceId{0}));
  BOOST_CHECK(!mftBinding.getScratchTransitionSlot(TransitionId::invalid()));
  BOOST_CHECK(!itsBinding.getScratchCellSlot(CellTopologyId::invalid()));

  for (const auto id : itsBinding.getGlobalTransitions()) {
    const auto& transition = combined.view.topology.getTransition(id);
    BOOST_CHECK(transition.from.value() < ITSNLayers);
    BOOST_CHECK(transition.to.value() < ITSNLayers);
  }
  for (const auto id : mftBinding.getGlobalTransitions()) {
    const auto& transition = combined.view.topology.getTransition(id);
    BOOST_CHECK_GE(transition.from.value(), ITSNLayers);
    BOOST_CHECK_GE(transition.to.value(), ITSNLayers);
  }
  checkParity<ITSNLayers>(itsBinding, combined.view, 0, 3);
  checkParity<MFTNLayers>(mftBinding, combined.view, ITSNLayers, 5);
}

BOOST_AUTO_TEST_CASE(FilteredRoadStartsKeepTheirGlobalCellIdentity)
{
  CombinedLayout combined;
  const auto its = DetectorTraversalBinding::build(combined.view, o2::detectors::DetID::ITS, ClusterSourceId{0}, itsMask(), ordered(0, ITSNLayers));
  const auto mft = DetectorTraversalBinding::build(combined.view, o2::detectors::DetID::MFT, ClusterSourceId{1}, mftMask(), ordered(ITSNLayers, MFTNLayers));
  BOOST_REQUIRE(its.ok());
  BOOST_REQUIRE(mft.ok());
  for (const auto id : its.binding->getGlobalRoadStartCells()) {
    BOOST_CHECK(its.binding->getScratchCellSlot(id));
    BOOST_CHECK(!mft.binding->getScratchCellSlot(id));
  }
  for (const auto id : mft.binding->getGlobalRoadStartCells()) {
    BOOST_CHECK(mft.binding->getScratchCellSlot(id));
    BOOST_CHECK(!its.binding->getScratchCellSlot(id));
  }
}
