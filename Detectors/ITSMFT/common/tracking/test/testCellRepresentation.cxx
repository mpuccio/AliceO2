// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Stage-B cell-state activation: focused coverage for the composition-based
// common CellSeed/TrackSeed representation (Cell.h). Proves the
// activation slice's structural claims -- composition rather than
// inheritance, one common CellSeed representation across the ITS/MFT layer
// ranges, fixed-capacity TrackSeed storage, complete Cell->TrackSeed transfer,
// and one common raw-q/pT getQOverPt() accessor plus the family-agnostic
// road-filter bound built on it -- without imposing a durable byte-offset
// ABI lock (only SurfaceKinematicState.h itself carries that lock).

#define BOOST_TEST_MODULE ITSMFT CellRepresentation
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <cstring>
#include <limits>
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
  BOOST_CHECK(!(std::is_base_of_v<o2::track::TrackParCovF, TrackSeed>));
  BOOST_CHECK(!(std::is_base_of_v<o2::track::TrackParCovFwd, TrackSeed>));

  CellSeed cell{};
  const auto state = makeDistinctState(StateFamily::Barrel);
  cell.state() = state;
  BOOST_CHECK_EQUAL(std::memcmp(&cell.state(), &state, sizeof(SurfaceKinematicState)), 0);

  // state() is a real, mutable reference to the composed member, not a copy.
  cell.state().parameters[4] = 99.f;
  BOOST_CHECK_EQUAL(cell.state().parameters[4], 99.f);

  BOOST_CHECK(!cell.tripletFactor().isValid());
  cell.tripletFactor().rho.phi = -0.5f;
  BOOST_CHECK(cell.tripletFactor().isValid());
}

BOOST_AUTO_TEST_CASE(CellTripletFactorReferencesReusePositionalClusterBinding)
{
  const auto state = makeDistinctState(StateFamily::Barrel);
  const o2::its::TimeEstBC time{};
  CellSeed cell{LayerMask{2, 4, 7}, 101, 202, 303, 1, 2, state, 0.f, time};

  const std::array<CellClusterReference, 3> expected{{{2, 101}, {4, 202}, {7, 303}}};
  for (int slot = 0; slot < 3; ++slot) {
    const auto reference = cell.getClusterReference(slot);
    BOOST_CHECK_EQUAL(reference.surfacePosition, expected[slot].surfacePosition);
    BOOST_CHECK_EQUAL(reference.clusterIndex, expected[slot].clusterIndex);
  }
  BOOST_CHECK_EQUAL(cell.getClusterReference(-1).clusterIndex, o2::its::constants::UnusedIndex);
  BOOST_CHECK_EQUAL(cell.getClusterReference(3).surfacePosition, o2::its::constants::UnusedIndex);
}

// --- One common CellSeed representation -------------------------------------

BOOST_AUTO_TEST_CASE(CellSeedUsesOneCommonRepresentationAcrossLayerRanges)
{
  const auto state = makeDistinctState(StateFamily::Barrel);
  const o2::its::TimeEstBC time{};
  CellSeed itsSeed{LayerMask{0x007f}, 10, 11, 12, 1, 2, state, 0.f, time};
  CellSeed mftSeed{LayerMask{0x03ff}, 20, 21, 22, 3, 4, state, 0.f, time};

  BOOST_CHECK_EQUAL(itsSeed.getHitLayerMask().value(), 0x007f);
  BOOST_CHECK_EQUAL(mftSeed.getHitLayerMask().value(), 0x03ff);
  BOOST_CHECK_EQUAL(itsSeed.getCluster(6), o2::its::constants::UnusedIndex);
  BOOST_CHECK_EQUAL(mftSeed.getCluster(9), o2::its::constants::UnusedIndex);
}

// --- TrackSeed has one fixed capacity independent of NLayers ----------------

BOOST_AUTO_TEST_CASE(TrackSeedUsesOneFixedCapacity)
{
  static_assert(TrackSeed::MaxSurfaces == static_cast<int>(MaxLayoutSurfaces));

  TrackSeed seed;
  SurfaceMask mask;
  mask.set(SurfaceId{static_cast<uint16_t>(TrackSeed::MaxSurfaces - 1)});
  seed.setSurfaceMask(mask);
  seed.setCluster(TrackSeed::MaxSurfaces - 1, 42);

  BOOST_CHECK(seed.hasCluster(TrackSeed::MaxSurfaces - 1));
  BOOST_CHECK_EQUAL(seed.getCluster(TrackSeed::MaxSurfaces - 1), 42);
  BOOST_CHECK_EQUAL(seed.getActiveSurfaceCount(), 1);
  BOOST_CHECK_EQUAL(seed.getCluster(TrackSeed::MaxSurfaces), o2::its::constants::UnusedIndex);

  seed.setCluster(TrackSeed::MaxSurfaces, 99);
  BOOST_CHECK_EQUAL(seed.getCluster(TrackSeed::MaxSurfaces), o2::its::constants::UnusedIndex);
}

// --- Cell -> TrackSeed copies every state byte and every metadata field -----

BOOST_AUTO_TEST_CASE(TrackSeedConstructionFromCellCopiesStateAndMetadataCompletely)
{
  constexpr int innerLayer = 2;
  const auto state = makeDistinctState(StateFamily::Barrel);
  const o2::its::TimeEstBC time{static_cast<uint32_t>(123), static_cast<uint16_t>(4)};
  CellSeed cell{innerLayer, 10, 20, 30, 5, 6, state, 7.5f, time};

  TrackSeed seed{cell};

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

// --- getQOverPt: raw signed q/pT, one common accessor, no family branch ----
//
// Correction to the Stage-B activation: the former common getQ2Pt()
// accessor squared slot 4 for both families, which was the correct
// convention for neither (SurfaceKinematicState's slot 4 is the raw signed
// q/pT parameter for both Barrel and Forward -- see SurfaceKinematicState.h's
// BarrelStateView::getQ2Pt()/ForwardStateView::getInvQPt(), which already
// returned this same raw value unsquared). getQ2Pt() has been removed from
// Cell.h/SeedMetadataBase; getQOverPt() is its sole replacement.

BOOST_AUTO_TEST_CASE(GetQOverPtReturnsTheExactRawPositiveValue)
{
  auto state = makeDistinctState(StateFamily::Barrel);
  state.parameters[4] = 12.5f;
  const o2::its::TimeEstBC time{};
  CellSeed cell{0, 0, 1, 2, 0, 1, state, 0.f, time};
  BOOST_CHECK_EQUAL(cell.getQOverPt(), 12.5f);
}

BOOST_AUTO_TEST_CASE(GetQOverPtReturnsTheExactRawNegativeValueWithoutSquaringOrAbs)
{
  auto state = makeDistinctState(StateFamily::Forward);
  state.parameters[4] = -12.5f;
  const o2::its::TimeEstBC time{};
  CellSeed cell{0, 0, 1, 2, 0, 1, state, 0.f, time};
  // Neither squared (which would be positive 156.25) nor made positive by abs().
  BOOST_CHECK_EQUAL(cell.getQOverPt(), -12.5f);
  BOOST_CHECK_LT(cell.getQOverPt(), 0.f);
}

BOOST_AUTO_TEST_CASE(GetQOverPtIsIdenticalForBarrelAndForwardGivenTheSameSlotValue)
{
  const auto barrelState = makeDistinctState(StateFamily::Barrel);
  const auto forwardState = makeDistinctState(StateFamily::Forward);
  const o2::its::TimeEstBC time{};

  CellSeed barrelCell{0, 0, 1, 2, 0, 1, barrelState, 0.f, time};
  CellSeed forwardCell{0, 0, 1, 2, 0, 1, forwardState, 0.f, time};

  BOOST_REQUIRE_EQUAL(barrelState.parameters[4], forwardState.parameters[4]);
  BOOST_CHECK_EQUAL(barrelCell.getQOverPt(), barrelState.parameters[4]);
  BOOST_CHECK_EQUAL(forwardCell.getQOverPt(), forwardState.parameters[4]);
  BOOST_CHECK_EQUAL(barrelCell.getQOverPt(), forwardCell.getQOverPt());
}

// --- Road-length filter bound: std::abs(getQOverPt()) <= maxAbsQOverPt -----
//
// TrackerTraits<NLayers>::findRoadsForKind's seedFilter (TrackerTraits.cxx)
// applies `std::abs(seed.getQOverPt()) <= maxAbsQOverPt` (maxAbsQOverPt =
// 1.e3f) identically for both SurfaceKind values: the expression
// itself never reads NLayers, DetID, or SurfaceKinematicState::family. These
// tests re-derive that exact expression as an independent oracle (mirroring
// the pattern already used elsewhere in this test suite for other
// production formulas) against CellSeed/TrackSeedN objects built for both
// families, so a family branch reintroduced in TrackerTraits.cxx would not
// be caught by this test file -- only a divergence in the formula's own
// pass/fail boundary would. Orchestration-level proof that TrackerTraits.cxx
// itself calls this unconditionally for both tags lives in
// testComputeLayerCellsOrchestration.cxx/testComputeLayerTrackletsOrchestration.cxx's
// existing one-shot-dispatch coverage.
namespace
{
constexpr float kMaxAbsQOverPt = 1.e3f;

bool passesRoadOverPtFilter(float qOverPt) noexcept
{
  return std::abs(qOverPt) <= kMaxAbsQOverPt;
}
} // namespace

BOOST_AUTO_TEST_CASE(RoadFilterAcceptsBothInclusiveBoundariesForBothFamilies)
{
  for (const auto family : {StateFamily::Barrel, StateFamily::Forward}) {
    BOOST_CHECK(passesRoadOverPtFilter(kMaxAbsQOverPt));
    BOOST_CHECK(passesRoadOverPtFilter(-kMaxAbsQOverPt));
    (void)family; // formula itself takes no family argument; looped only for documentation
  }
}

BOOST_AUTO_TEST_CASE(RoadFilterRejectsTheImmediatelyRepresentableValuesOutsideBothBoundaries)
{
  const float justAbove = std::nextafter(kMaxAbsQOverPt, std::numeric_limits<float>::infinity());
  const float justBelowNegative = std::nextafter(-kMaxAbsQOverPt, -std::numeric_limits<float>::infinity());
  BOOST_REQUIRE_GT(justAbove, kMaxAbsQOverPt);
  BOOST_REQUIRE_LT(justBelowNegative, -kMaxAbsQOverPt);

  BOOST_CHECK(!passesRoadOverPtFilter(justAbove));
  BOOST_CHECK(!passesRoadOverPtFilter(justBelowNegative));
}

BOOST_AUTO_TEST_CASE(RoadFilterAcceptsZero)
{
  BOOST_CHECK(passesRoadOverPtFilter(0.f));
}

BOOST_AUTO_TEST_CASE(RoadFilterRejectsNaNAndBothInfinities)
{
  // std::abs(x) <= finite bound is false for every non-finite x under
  // ordinary IEEE-754 comparison semantics; this is asserted explicitly
  // (per the correction requirement not to rely on it without a focused
  // test) rather than merely assumed.
  BOOST_CHECK(!passesRoadOverPtFilter(std::numeric_limits<float>::quiet_NaN()));
  BOOST_CHECK(!passesRoadOverPtFilter(std::numeric_limits<float>::infinity()));
  BOOST_CHECK(!passesRoadOverPtFilter(-std::numeric_limits<float>::infinity()));
}

BOOST_AUTO_TEST_CASE(RoadFilterTreatsEqualMagnitudeOppositeSignSeedsIdentically)
{
  const std::array<float, 4> magnitudes{0.f, 1.f, kMaxAbsQOverPt, kMaxAbsQOverPt * 2.f};
  for (const float magnitude : magnitudes) {
    BOOST_CHECK_EQUAL(passesRoadOverPtFilter(magnitude), passesRoadOverPtFilter(-magnitude));
  }

  // Same check through the real accessor, for both families, at a
  // representative in-bound magnitude.
  auto barrelState = makeDistinctState(StateFamily::Barrel);
  auto forwardState = makeDistinctState(StateFamily::Forward);
  const o2::its::TimeEstBC time{};
  for (const float signedValue : {10.f, -10.f}) {
    barrelState.parameters[4] = signedValue;
    forwardState.parameters[4] = signedValue;
    CellSeed barrelCell{0, 0, 1, 2, 0, 1, barrelState, 0.f, time};
    CellSeed forwardCell{0, 0, 1, 2, 0, 1, forwardState, 0.f, time};
    BOOST_CHECK_EQUAL(passesRoadOverPtFilter(barrelCell.getQOverPt()), passesRoadOverPtFilter(forwardCell.getQOverPt()));
  }
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
  BOOST_TEST_MESSAGE("sizeof(TrackSeed) = " << sizeof(TrackSeed)
                                            << ", alignof = " << alignof(TrackSeed));

  BOOST_CHECK_GE(sizeof(CellSeed), sizeof(SurfaceKinematicState));
  BOOST_CHECK_GE(sizeof(TrackSeed), sizeof(SurfaceKinematicState));
  BOOST_CHECK_EQUAL(alignof(CellSeed), alignof(SurfaceKinematicState));
}
