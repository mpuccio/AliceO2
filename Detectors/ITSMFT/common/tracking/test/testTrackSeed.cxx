// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// Gate 4 M6f: focused coverage for the fixed-capacity, GPU-portable,
// non-templated TrackSeed. Container/value-type validation mirrors
// testCellRepresentation.cxx's style for the live Cell.h types.

#define BOOST_TEST_MODULE ITSMFT TrackSeed
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cstring>
#include <type_traits>

#include "ITSMFTTracking/Cell.h"

using namespace o2::itsmft::tracking;

namespace
{
SurfaceKinematicState makeDistinctState(SurfaceKind family)
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
  state.kind = family;
  state.flags = 0;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  return state;
}
} // namespace

// --- Capacity source: MaxLayoutSurfaces, no invented bound -------------------

BOOST_AUTO_TEST_CASE(CapacityIsExactlyMaxLayoutSurfaces)
{
  // The one existing generic compile-time capacity constant this library
  // already uses to bound every owned-surface position (LayerMask,
  // the retired topology owner -- IdTypes.h). TrackSeed
  // reuses it rather than inventing a new one.
  static_assert(TrackSeed::MaxSurfaces == static_cast<int>(MaxLayoutSurfaces));
  static_assert(TrackSeed::MaxSurfaces == 32);
  BOOST_CHECK(true); // compile-time proof above; case exists so it appears in test output.
}

// --- Not a template / not a detector-specific alias --------------------------

BOOST_AUTO_TEST_CASE(TrackSeedIsAConcreteNonTemplatedType)
{
  static_assert(std::is_final_v<TrackSeed>);
  static_assert(std::is_same_v<decltype(TrackSeed{}), TrackSeed>);
  BOOST_CHECK(true);
}

// --- Construction from CellSeed: complete metadata transfer, holes preserved -

BOOST_AUTO_TEST_CASE(ConstructionFromCellCopiesStateAndMetadataCompletely)
{
  constexpr int innerLayer = 2;
  const auto state = makeDistinctState(SurfaceKind::Cylinder);
  const o2::its::TimeEstBC time{static_cast<uint32_t>(123), static_cast<uint16_t>(4)};
  CellSeed cell{innerLayer, 10, 20, 30, 5, 6, time};

  TrackSeed seed{cell, state, 7.5f};

  BOOST_CHECK_EQUAL(std::memcmp(&seed.state(), &state, sizeof(SurfaceKinematicState)), 0);
  BOOST_CHECK_EQUAL(seed.getChi2(), 7.5f);
  BOOST_CHECK_EQUAL(seed.getLevel(), cell.getLevel());
  BOOST_CHECK_EQUAL(seed.getTimeStamp().getTimeStamp(), cell.getTimeStamp().getTimeStamp());
  BOOST_CHECK_EQUAL(seed.getTimeStamp().getTimeStampError(), cell.getTimeStamp().getTimeStampError());
  BOOST_CHECK_EQUAL(seed.getFirstTrackletIndex(), cell.getFirstTrackletIndex());
  BOOST_CHECK_EQUAL(seed.getSecondTrackletIndex(), cell.getSecondTrackletIndex());
  BOOST_CHECK_EQUAL(seed.getActiveLayerCount(), 3);
  BOOST_CHECK_EQUAL(seed.getCluster(innerLayer), cell.getFirstClusterIndex());
  BOOST_CHECK_EQUAL(seed.getCluster(innerLayer + 1), cell.getSecondClusterIndex());
  BOOST_CHECK_EQUAL(seed.getCluster(innerLayer + 2), cell.getThirdClusterIndex());
  BOOST_CHECK_EQUAL(seed.getFirstClusterIndex(), cell.getFirstClusterIndex());
  BOOST_CHECK_EQUAL(seed.getSecondClusterIndex(), cell.getSecondClusterIndex());
  BOOST_CHECK_EQUAL(seed.getThirdClusterIndex(), cell.getThirdClusterIndex());

  // Every position outside the cell's own 3 hit layers is a hole: not
  // active, and reads back UnusedIndex.
  for (int position = 0; position < 16; ++position) {
    if (position == innerLayer || position == innerLayer + 1 || position == innerLayer + 2) {
      BOOST_CHECK(seed.hasCluster(position));
    } else {
      BOOST_CHECK(!seed.hasCluster(position));
      BOOST_CHECK_EQUAL(seed.getCluster(position), o2::its::constants::UnusedIndex);
    }
  }
}

// --- Holes and all active slots, at full MaxSurfaces capacity ----------------

BOOST_AUTO_TEST_CASE(SupportsHolesAndAllActiveSlotsAtFullCapacity)
{
  TrackSeed seed{};
  LayerMask mask;
  for (int position = 0; position < TrackSeed::MaxSurfaces; ++position) {
    if (position % 2 == 0) {
      seed.setCluster(position, position * 10);
      mask.set(static_cast<uint16_t>(position));
    } else {
      seed.setCluster(position, o2::its::constants::UnusedIndex); // explicit hole
    }
  }
  seed.setHitLayerMask(mask);

  BOOST_CHECK_EQUAL(seed.getActiveLayerCount(), TrackSeed::MaxSurfaces / 2);
  for (int position = 0; position < TrackSeed::MaxSurfaces; ++position) {
    if (position % 2 == 0) {
      BOOST_CHECK(seed.hasCluster(position));
      BOOST_CHECK_EQUAL(seed.getCluster(position), position * 10);
    } else {
      BOOST_CHECK(!seed.hasCluster(position));
      BOOST_CHECK_EQUAL(seed.getCluster(position), o2::its::constants::UnusedIndex);
    }
  }

  // All MaxSurfaces positions active at once (no holes): the full-capacity
  // boundary itself, not just "some slots".
  LayerMask full;
  for (int position = 0; position < TrackSeed::MaxSurfaces; ++position) {
    seed.setCluster(position, position);
    full.set(static_cast<uint16_t>(position));
  }
  seed.setHitLayerMask(full);
  BOOST_CHECK_EQUAL(seed.getActiveLayerCount(), TrackSeed::MaxSurfaces);
  for (int position = 0; position < TrackSeed::MaxSurfaces; ++position) {
    BOOST_CHECK(seed.hasCluster(position));
    BOOST_CHECK_EQUAL(seed.getCluster(position), position);
  }
}

// --- Rejects/outcomes safely at the fixed-capacity boundary ------------------

BOOST_AUTO_TEST_CASE(GetClusterIsBoundsCheckedBeyondCapacity)
{
  TrackSeed seed{};
  seed.setCluster(0, 111);

  BOOST_CHECK_EQUAL(seed.getCluster(-1), o2::its::constants::UnusedIndex);
  BOOST_CHECK_EQUAL(seed.getCluster(TrackSeed::MaxSurfaces), o2::its::constants::UnusedIndex);
  BOOST_CHECK_EQUAL(seed.getCluster(TrackSeed::MaxSurfaces + 1000), o2::its::constants::UnusedIndex);
  BOOST_CHECK(!seed.hasCluster(-1));
  BOOST_CHECK(!seed.hasCluster(TrackSeed::MaxSurfaces));

  // Out-of-range writes are safely dropped, not undefined behavior, and
  // never corrupt an in-range slot.
  seed.setCluster(-1, 999);
  seed.setCluster(TrackSeed::MaxSurfaces, 999);
  seed.setCluster(TrackSeed::MaxSurfaces * 4, 999);
  BOOST_CHECK_EQUAL(seed.getCluster(0), 111);
}

// --- getClusterBySlot boundary: fewer active slots than requested ------------

BOOST_AUTO_TEST_CASE(SlotAccessorsReturnUnusedIndexBeyondActiveCount)
{
  const o2::its::TimeEstBC time{};
  const auto state = makeDistinctState(SurfaceKind::Disk);
  CellSeed cell{0, 100, 200, 300, 1, 2, time};
  TrackSeed seed{cell, state, 0.f}; // exactly 3 active slots

  BOOST_CHECK_EQUAL(seed.getFirstClusterIndex(), 100);
  BOOST_CHECK_EQUAL(seed.getSecondClusterIndex(), 200);
  BOOST_CHECK_EQUAL(seed.getThirdClusterIndex(), 300);
  BOOST_CHECK_EQUAL(seed.getActiveLayerCount(), 3);

  // An empty seed (no active slots at all) safely returns UnusedIndex for
  // every requested slot rather than reading past an empty mask.
  TrackSeed empty{};
  BOOST_CHECK_EQUAL(empty.getFirstClusterIndex(), o2::its::constants::UnusedIndex);
  BOOST_CHECK_EQUAL(empty.getActiveLayerCount(), 0);
}

// --- getQOverPt: same raw-signed-value convention as CellSeed -------------

BOOST_AUTO_TEST_CASE(GetQOverPtReturnsTheExactRawValueWithoutSquaringOrAbs)
{
  auto state = makeDistinctState(SurfaceKind::Disk);
  state.parameters[4] = -12.5f;
  const o2::its::TimeEstBC time{};
  CellSeed cell{0, 0, 1, 2, 0, 1, time};
  TrackSeed seed{cell, state, 0.f};
  BOOST_CHECK_EQUAL(seed.getQOverPt(), -12.5f);
  BOOST_CHECK_LT(seed.getQOverPt(), 0.f);
}

// --- GPU value-type property: the exact applicable one, not an invalid one --

BOOST_AUTO_TEST_CASE(TrackSeedIsTriviallyCopyableButNotStandardLayout)
{
  // TrackSeed is GPUhd()-annotated throughout, exactly like CellSeed, and
  // must remain copyable by value across the
  // host/device boundary: is_trivially_copyable is the property this
  // actually requires (byte-for-byte copyable, no user-provided copy/move/
  // destructor logic) and the one it satisfies -- already static_assert'd
  // in Cell.h itself, re-asserted here so a regression fails this focused
  // test suite too, not just a downstream build.
  static_assert(std::is_trivially_copyable_v<TrackSeed>);

  // Deliberately NOT standard-layout: its embedded o2::its::TimeEstBC
  // (o2::dataformats::TimeStampWithError<uint32_t, uint16_t>, deriving from
  // TimeStamp<uint32_t>) declares non-static data members in more than one
  // class of its own hierarchy. This is a property of TimeEstBC itself, not
  // anything TrackSeed adds -- CellSeed embeds the same member and is not
  // standard-layout either (it
  // carries a standard-layout static_assert in this codebase). Asserted
  // here explicitly, not merely omitted, so this is a stated fact rather
  // than an untested assumption.
  static_assert(!std::is_standard_layout_v<TrackSeed>);

  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE(TrackSeedSizeCharacterization)
{
  BOOST_TEST_MESSAGE("sizeof(TrackSeed) = " << sizeof(TrackSeed) << ", alignof(TrackSeed) = " << alignof(TrackSeed));
  BOOST_CHECK_GE(sizeof(TrackSeed), sizeof(SurfaceKinematicState));
}
