// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Stage-B cell-state activation: focused coverage for the composition-based
// common CellSeed/TrackSeedTpl<NLayers> representation (Cell.h). Proves the
// activation slice's structural claims -- composition rather than
// inheritance, one common CellSeed type regardless of NLayers, TrackSeed
// differing only by cluster-array width, complete Cell->TrackSeed transfer,
// and one common getQ2Pt() formula -- without imposing a durable byte-offset
// ABI lock (only SurfaceKinematicState.h itself carries that lock).

#define BOOST_TEST_MODULE ITSMFT CellRepresentation
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cstring>
#include <type_traits>

#include "ITSMFTTracking/Cell.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackFwd.h"

using namespace o2::itsmft::tracking;

namespace
{

SurfaceKinematicState makeDistinctState(StateFamily family)
{
  SurfaceKinematicState state{};
  for (int i = 0; i < 5; ++i) {
    state.parameters[i] = 1.f + static_cast<float>(i);
  }
  for (int i = 0; i < 15; ++i) {
    state.covariance[i] = 10.f + static_cast<float>(i);
  }
  state.referenceCoordinate = 42.f;
  state.alpha = 0.25f;
  state.family = family;
  state.flags = 0;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  return state;
}

} // namespace

// --- Composition, not inheritance -------------------------------------------

BOOST_AUTO_TEST_CASE(CellSeedComposesSurfaceKinematicStateAndDoesNotInheritLegacyTrackTypes)
{
  BOOST_CHECK(!(std::is_base_of_v<o2::track::TrackParCovF, CellSeed>));
  BOOST_CHECK(!(std::is_base_of_v<o2::track::TrackParCovFwd, CellSeed>));
  BOOST_CHECK(!(std::is_base_of_v<o2::track::TrackParCovF, TrackSeedTpl<ITSNLayers>>));
  BOOST_CHECK(!(std::is_base_of_v<o2::track::TrackParCovFwd, TrackSeedTpl<MFTNLayers>>));

  CellSeed cell{};
  const auto state = makeDistinctState(StateFamily::Barrel);
  cell.state() = state;
  BOOST_CHECK_EQUAL(std::memcmp(&cell.state(), &state, sizeof(SurfaceKinematicState)), 0);

  // state() is a real, mutable reference to the composed member, not a copy.
  cell.state().parameters[4] = 99.f;
  BOOST_CHECK_EQUAL(cell.state().parameters[4], 99.f);
}

// --- One common CellSeed type regardless of NLayers -------------------------

BOOST_AUTO_TEST_CASE(CellSeedNResolvesToTheSameTypeForEveryNLayers)
{
  static_assert(std::is_same_v<CellSeedN<ITSNLayers>, CellSeed>);
  static_assert(std::is_same_v<CellSeedN<MFTNLayers>, CellSeed>);
  static_assert(std::is_same_v<CellSeedN<ITSNLayers>, CellSeedN<MFTNLayers>>);
  BOOST_CHECK(true); // compile-time proof above; this case exists so it appears in test output.
}

// --- TrackSeed differs only by NLayers cluster-array width ------------------

BOOST_AUTO_TEST_CASE(TrackSeedNDiffersOnlyByClusterArrayWidth)
{
  static_assert(!std::is_same_v<TrackSeedN<ITSNLayers>, TrackSeedN<MFTNLayers>>);
  static_assert(std::is_same_v<TrackSeedN<ITSNLayers>, TrackSeedTpl<ITSNLayers>>);
  static_assert(std::is_same_v<TrackSeedN<MFTNLayers>, TrackSeedTpl<MFTNLayers>>);

  // Both wrap the identical SeedMetadataBase<NLayers> shape apart from the
  // NLayers-sized raw int cluster array; the only expected size delta is
  // exactly that many extra ints (no other family-dependent storage).
  constexpr size_t expectedDelta = static_cast<size_t>(MFTNLayers - ITSNLayers) * sizeof(int);
  const size_t actualDelta = sizeof(TrackSeedN<MFTNLayers>) - sizeof(TrackSeedN<ITSNLayers>);
  BOOST_CHECK_EQUAL(actualDelta, expectedDelta);
}

// --- Cell -> TrackSeed copies every state byte and every metadata field -----

BOOST_AUTO_TEST_CASE(TrackSeedConstructionFromCellCopiesStateAndMetadataCompletely)
{
  constexpr int innerLayer = 2;
  const auto state = makeDistinctState(StateFamily::Barrel);
  const o2::its::TimeEstBC time{static_cast<uint32_t>(123), static_cast<uint16_t>(4)};
  CellSeed cell{innerLayer, 10, 20, 30, 5, 6, state, 7.5f, time};

  TrackSeedN<ITSNLayers> seed{cell};

  BOOST_CHECK_EQUAL(std::memcmp(&seed.state(), &cell.state(), sizeof(SurfaceKinematicState)), 0);
  BOOST_CHECK_EQUAL(seed.getChi2(), cell.getChi2());
  BOOST_CHECK_EQUAL(seed.getLevel(), cell.getLevel());
  BOOST_CHECK_EQUAL(seed.getTimeStamp().getTimeStamp(), cell.getTimeStamp().getTimeStamp());
  BOOST_CHECK_EQUAL(seed.getTimeStamp().getTimeStampError(), cell.getTimeStamp().getTimeStampError());
  BOOST_CHECK_EQUAL(seed.getHitLayerMask().value(), cell.getHitLayerMask().value());
  BOOST_CHECK_EQUAL(seed.getFirstTrackletIndex(), cell.getFirstTrackletIndex());
  BOOST_CHECK_EQUAL(seed.getSecondTrackletIndex(), cell.getSecondTrackletIndex());
  BOOST_CHECK_EQUAL(seed.getCluster(innerLayer), cell.getFirstClusterIndex());
  BOOST_CHECK_EQUAL(seed.getCluster(innerLayer + 1), cell.getSecondClusterIndex());
  BOOST_CHECK_EQUAL(seed.getCluster(innerLayer + 2), cell.getThirdClusterIndex());
}

// --- getQ2Pt: one common formula, no family branch --------------------------

BOOST_AUTO_TEST_CASE(GetQ2PtSquaresSlotFourIdenticallyForBarrelAndForward)
{
  const auto barrelState = makeDistinctState(StateFamily::Barrel);
  const auto forwardState = makeDistinctState(StateFamily::Forward);
  const o2::its::TimeEstBC time{};

  CellSeed barrelCell{0, 0, 1, 2, 0, 1, barrelState, 0.f, time};
  CellSeed forwardCell{0, 0, 1, 2, 0, 1, forwardState, 0.f, time};

  const float expected = barrelState.parameters[4] * barrelState.parameters[4];
  BOOST_REQUIRE_EQUAL(barrelState.parameters[4], forwardState.parameters[4]);
  BOOST_CHECK_EQUAL(barrelCell.getQ2Pt(), expected);
  BOOST_CHECK_EQUAL(forwardCell.getQ2Pt(), expected);
  BOOST_CHECK_EQUAL(barrelCell.getQ2Pt(), forwardCell.getQ2Pt());
}

// --- Size/alignment characterization (no static ABI lock) ------------------

BOOST_AUTO_TEST_CASE(CellSeedAndTrackSeedSizeAlignmentCharacterization)
{
  // Reported, not locked: these values may legitimately change with future
  // metadata/state evolution. Only cross-checked against SurfaceKinematicState
  // (the one durably ABI-locked type, SurfaceKinematicState.h) so a
  // regression that shrinks CellSeed below its composed member would fail
  // loudly, without pinning an exact byte count here.
  BOOST_TEST_MESSAGE("sizeof(SurfaceKinematicState) = " << sizeof(SurfaceKinematicState));
  BOOST_TEST_MESSAGE("sizeof(CellSeed) = " << sizeof(CellSeed) << ", alignof(CellSeed) = " << alignof(CellSeed));
  BOOST_TEST_MESSAGE("sizeof(TrackSeedN<ITSNLayers>) = " << sizeof(TrackSeedN<ITSNLayers>)
                                                         << ", alignof = " << alignof(TrackSeedN<ITSNLayers>));
  BOOST_TEST_MESSAGE("sizeof(TrackSeedN<MFTNLayers>) = " << sizeof(TrackSeedN<MFTNLayers>)
                                                         << ", alignof = " << alignof(TrackSeedN<MFTNLayers>));

  BOOST_CHECK_GE(sizeof(CellSeed), sizeof(SurfaceKinematicState));
  BOOST_CHECK_GE(sizeof(TrackSeedN<ITSNLayers>), sizeof(SurfaceKinematicState));
  BOOST_CHECK_GE(sizeof(TrackSeedN<MFTNLayers>), sizeof(TrackSeedN<ITSNLayers>));
  BOOST_CHECK_EQUAL(alignof(CellSeed), alignof(SurfaceKinematicState));
}
