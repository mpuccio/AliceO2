// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFTSurfaceLinearizationReference
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cstring>
#include <limits>
#include <type_traits>

#include "ITSMFTTracking/SurfaceKinematicState.h"

namespace
{
using namespace o2::itsmft::tracking;

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

SurfaceKinematicState makeBarrelState()
{
  SurfaceKinematicState state{};
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.2f;
  state.parameters[3] = -0.35f;
  state.parameters[4] = 0.8f;
  state.referenceCoordinate = 4.f;
  state.alpha = 0.3f;
  state.kind = SurfaceKind::Cylinder;
  state.flags = 0x5a;
  state.absCharge = 1;
  state.pid = o2::track::PID::Kaon;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
    }
  }
  return state;
}

SurfaceKinematicState makeForwardState()
{
  SurfaceKinematicState state{};
  state.parameters[0] = 5.f;
  state.parameters[1] = -2.f;
  state.parameters[2] = 0.4f;
  state.parameters[3] = -1.2f;
  state.parameters[4] = 0.3f;
  state.referenceCoordinate = -45.f;
  state.alpha = 0.f;
  state.kind = SurfaceKind::Disk;
  state.flags = 0x11;
  state.absCharge = 1;
  state.pid = o2::track::PID::Muon;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.02f * (row + 1) : 0.0001f * (row + column + 1);
    }
  }
  return state;
}

SurfaceLinearizationReference sentinel()
{
  SurfaceLinearizationReference ref{};
  ref.kind = SurfaceKind::Cylinder;
  ref.referenceCoordinate = 999.f;
  ref.alpha = -1.f;
  for (uint8_t i = 0; i < 5; ++i) {
    ref.parameters[i] = 42.f + i;
  }
  return ref;
}

} // namespace

BOOST_AUTO_TEST_CASE(DefaultReferenceIsInvalid)
{
  SurfaceLinearizationReference ref{};
  BOOST_CHECK(ref.kind == SurfaceKind::Undefined);
  BOOST_CHECK(!ref.hasRecognizedKind());
  for (float value : ref.parameters) {
    BOOST_CHECK_EQUAL(value, 0.f);
  }
  BOOST_CHECK_EQUAL(ref.referenceCoordinate, 0.f);
  BOOST_CHECK_EQUAL(ref.alpha, 0.f);
}

BOOST_AUTO_TEST_CASE(LayoutIsCurrentInMemoryCommitment)
{
  static_assert(std::is_standard_layout_v<SurfaceLinearizationReference>);
  static_assert(std::is_trivially_copyable_v<SurfaceLinearizationReference>);
  static_assert(sizeof(SurfaceLinearizationReference) == 32);
  static_assert(alignof(SurfaceLinearizationReference) == 4);
}

BOOST_AUTO_TEST_CASE(FactoryMapsBarrelStateFieldsExactly)
{
  const auto state = makeBarrelState();
  SurfaceLinearizationReference ref{};
  BOOST_REQUIRE(makeLinearizationReference(state, ref));
  BOOST_CHECK(ref.kind == SurfaceKind::Cylinder);
  BOOST_CHECK_EQUAL(ref.referenceCoordinate, state.referenceCoordinate);
  BOOST_CHECK_EQUAL(ref.alpha, state.alpha);
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_EQUAL(ref.parameters[i], state.parameters[i]);
  }
}

BOOST_AUTO_TEST_CASE(FactoryMapsForwardStateFieldsExactly)
{
  const auto state = makeForwardState();
  SurfaceLinearizationReference ref{};
  BOOST_REQUIRE(makeLinearizationReference(state, ref));
  BOOST_CHECK(ref.kind == SurfaceKind::Disk);
  BOOST_CHECK_EQUAL(ref.referenceCoordinate, state.referenceCoordinate);
  BOOST_CHECK_EQUAL(ref.alpha, state.alpha);
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_EQUAL(ref.parameters[i], state.parameters[i]);
  }
}

BOOST_AUTO_TEST_CASE(FactoryRejectsUnrecognizedFamily)
{
  auto state = makeBarrelState();
  state.kind = SurfaceKind::Undefined;
  auto ref = sentinel();
  const auto before = ref;
  BOOST_CHECK(!makeLinearizationReference(state, ref));
  BOOST_CHECK(bitEqual(ref, before));
}

BOOST_AUTO_TEST_CASE(FactoryTransactionalOnNonFiniteParameter)
{
  auto state = makeBarrelState();
  state.parameters[2] = std::numeric_limits<float>::quiet_NaN();
  auto ref = sentinel();
  const auto before = ref;
  BOOST_CHECK(!makeLinearizationReference(state, ref));
  BOOST_CHECK(bitEqual(ref, before));

  auto stateBadRef = makeBarrelState();
  stateBadRef.referenceCoordinate = std::numeric_limits<float>::infinity();
  auto ref2 = sentinel();
  const auto before2 = ref2;
  BOOST_CHECK(!makeLinearizationReference(stateBadRef, ref2));
  BOOST_CHECK(bitEqual(ref2, before2));

  auto stateBadAlpha = makeBarrelState();
  stateBadAlpha.alpha = std::numeric_limits<float>::quiet_NaN();
  auto ref3 = sentinel();
  const auto before3 = ref3;
  BOOST_CHECK(!makeLinearizationReference(stateBadAlpha, ref3));
  BOOST_CHECK(bitEqual(ref3, before3));
}

// Covariance/PID/absCharge independence: two states with identical
// parameters/referenceCoordinate/alpha/family but deliberately different
// covariance, PID, absCharge, and flags must produce byte-identical
// SurfaceLinearizationReference outputs, since the reference has no fields
// for any of those.
BOOST_AUTO_TEST_CASE(FactoryIsIndependentOfCovariancePidCharge)
{
  auto stateA = makeBarrelState();
  auto stateB = makeBarrelState();
  stateB.absCharge = 2;
  stateB.pid = o2::track::PID::Proton;
  stateB.flags = 0x00;
  for (uint8_t i = 0; i < 15; ++i) {
    stateB.covariance[i] = stateA.covariance[i] * 7.f + 3.f;
  }
  BOOST_REQUIRE(!bitEqual(stateA, stateB));

  SurfaceLinearizationReference refA{};
  SurfaceLinearizationReference refB{};
  BOOST_REQUIRE(makeLinearizationReference(stateA, refA));
  BOOST_REQUIRE(makeLinearizationReference(stateB, refB));
  BOOST_CHECK(bitEqual(refA, refB));
}

// Re-anchoring: the same factory serves as the transactional mutation
// primitive used to re-anchor an existing reference at a later fit-leg's
// state (the legacy `linRef = track.paramOut;` idiom).
BOOST_AUTO_TEST_CASE(FactoryReanchorsExistingReference)
{
  const auto stateA = makeBarrelState();
  auto stateB = makeBarrelState();
  stateB.parameters[0] += 1.f;
  stateB.referenceCoordinate += 2.f;

  SurfaceLinearizationReference ref{};
  BOOST_REQUIRE(makeLinearizationReference(stateA, ref));
  BOOST_REQUIRE(makeLinearizationReference(stateB, ref));
  BOOST_CHECK_EQUAL(ref.parameters[0], stateB.parameters[0]);
  BOOST_CHECK_EQUAL(ref.referenceCoordinate, stateB.referenceCoordinate);
}
