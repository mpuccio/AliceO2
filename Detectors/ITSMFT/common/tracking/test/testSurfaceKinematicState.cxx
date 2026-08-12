// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFTSurfaceKinematicState
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <type_traits>

#include "ITSMFTTracking/detail/SurfaceKinematicStateLegacyAdapters.h"

namespace
{
using namespace o2::itsmft::tracking;
namespace legacy = o2::itsmft::tracking::legacy;

o2::track::TrackParCovF makeBarrelState()
{
  o2::track::TrackParCovF::params_t parameters{1.f, 2.f, 0.3f, -0.4f, 0.5f};
  o2::track::TrackParCovF::covMat_t covariance{};
  for (uint8_t i = 0; i < 15; ++i) {
    covariance[i] = 10.f + i;
  }
  return {3.f, -0.2f, parameters, covariance, 2, o2::track::PID::Kaon};
}

o2::track::TrackParCovFwd makeForwardState()
{
  o2::track::SMatrix5 parameters{};
  o2::track::SMatrix55Sym covariance{};
  for (uint8_t i = 0; i < 5; ++i) {
    parameters(i) = 0.123456789 + i * 10.25;
  }
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      covariance(row, column) = 1.e-6 * (1 + packedCovarianceIndex(row, column));
    }
  }
  return {42.123456789, parameters, covariance, 7.};
}

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}
} // namespace

BOOST_AUTO_TEST_CASE(AbiIsLocked)
{
  static_assert(std::is_standard_layout_v<SurfaceKinematicState>);
  static_assert(std::is_trivially_copyable_v<SurfaceKinematicState>);
  static_assert(std::is_standard_layout_v<BarrelStateView>);
  static_assert(std::is_trivially_copyable_v<BarrelStateView>);
  static_assert(std::is_standard_layout_v<ForwardStateView>);
  static_assert(std::is_trivially_copyable_v<ForwardStateView>);

  BOOST_CHECK_EQUAL(sizeof(SurfaceKinematicState), 92U);
  BOOST_CHECK_EQUAL(alignof(SurfaceKinematicState), 4U);
  BOOST_CHECK_EQUAL(offsetof(SurfaceKinematicState, parameters), 0U);
  BOOST_CHECK_EQUAL(offsetof(SurfaceKinematicState, covariance), 20U);
  BOOST_CHECK_EQUAL(offsetof(SurfaceKinematicState, referenceCoordinate), 80U);
  BOOST_CHECK_EQUAL(offsetof(SurfaceKinematicState, alpha), 84U);
  BOOST_CHECK_EQUAL(offsetof(SurfaceKinematicState, family), 88U);
  BOOST_CHECK_EQUAL(offsetof(SurfaceKinematicState, flags), 89U);
  BOOST_CHECK_EQUAL(offsetof(SurfaceKinematicState, absCharge), 90U);
  BOOST_CHECK_EQUAL(offsetof(SurfaceKinematicState, pid), 91U);
  BOOST_CHECK_EQUAL(sizeof(BarrelStateView), sizeof(SurfaceKinematicState*));
  BOOST_CHECK_EQUAL(alignof(BarrelStateView), alignof(SurfaceKinematicState*));
  BOOST_CHECK_EQUAL(sizeof(ForwardStateView), sizeof(SurfaceKinematicState*));
  BOOST_CHECK_EQUAL(alignof(ForwardStateView), alignof(SurfaceKinematicState*));
}

BOOST_AUTO_TEST_CASE(DefaultStateIsDeterministicallyInvalid)
{
  const SurfaceKinematicState state{};
  for (const auto value : state.parameters) {
    BOOST_CHECK_EQUAL(value, 0.f);
  }
  for (const auto value : state.covariance) {
    BOOST_CHECK_EQUAL(value, 0.f);
  }
  BOOST_CHECK_EQUAL(state.referenceCoordinate, 0.f);
  BOOST_CHECK_EQUAL(state.alpha, 0.f);
  BOOST_CHECK(state.kind == SurfaceKind::Undefined);
  BOOST_CHECK(!state.hasRecognizedKind());
  BOOST_CHECK_EQUAL(state.flags, 0U);
  BOOST_CHECK_EQUAL(state.absCharge, 0U);
  BOOST_CHECK_EQUAL(state.pid.getID(), o2::track::PID::Pion);
}

BOOST_AUTO_TEST_CASE(BarrelViewMapsBarrelCoordinates)
{
  SurfaceKinematicState state{};
  state.kind = SurfaceKind::Cylinder;
  state.referenceCoordinate = 7.f;
  state.alpha = -0.7f;
  for (uint8_t i = 0; i < 5; ++i) {
    state.parameters[i] = 2.f + i;
  }
  BarrelStateView view;
  BOOST_CHECK(!view.isValid());
  BOOST_REQUIRE(makeBarrelStateView(state, view));
  BOOST_CHECK(view.isValid());
  BOOST_CHECK_EQUAL(view.getY(), 2.f);
  BOOST_CHECK_EQUAL(view.getZ(), 3.f);
  BOOST_CHECK_EQUAL(view.getSnp(), 4.f);
  BOOST_CHECK_EQUAL(view.getTgl(), 5.f);
  BOOST_CHECK_EQUAL(view.getQ2Pt(), 6.f);
  BOOST_CHECK_EQUAL(view.getReferenceX(), 7.f);
  BOOST_CHECK_EQUAL(view.getAlpha(), -0.7f);
  view.setCovariance(4, 1, 12.f);
  BOOST_CHECK_EQUAL(state.covariance[packedCovarianceIndex(4, 1)], 12.f);
}

BOOST_AUTO_TEST_CASE(ForwardViewMapsForwardCoordinates)
{
  SurfaceKinematicState state{};
  state.kind = SurfaceKind::Disk;
  state.referenceCoordinate = -12.f;
  for (uint8_t i = 0; i < 5; ++i) {
    state.parameters[i] = 3.f + i;
  }
  ForwardStateView view;
  BOOST_CHECK(!view.isValid());
  BOOST_REQUIRE(makeForwardStateView(state, view));
  BOOST_CHECK(view.isValid());
  BOOST_CHECK_EQUAL(view.getX(), 3.f);
  BOOST_CHECK_EQUAL(view.getY(), 4.f);
  BOOST_CHECK_EQUAL(view.getPhi(), 5.f);
  BOOST_CHECK_EQUAL(view.getTanl(), 6.f);
  BOOST_CHECK_EQUAL(view.getInvQPt(), 7.f);
  BOOST_CHECK_EQUAL(view.getReferenceZ(), -12.f);
  view.setCovariance(3, 2, 19.f);
  BOOST_CHECK_EQUAL(state.covariance[packedCovarianceIndex(3, 2)], 19.f);
}

BOOST_AUTO_TEST_CASE(ViewFactoriesRejectFamilyMismatchWithoutMutatingView)
{
  SurfaceKinematicState barrel{};
  barrel.kind = SurfaceKind::Cylinder;
  BarrelStateView barrelView;
  BOOST_REQUIRE(makeBarrelStateView(barrel, barrelView));
  SurfaceKinematicState forward{};
  forward.kind = SurfaceKind::Disk;
  BOOST_CHECK(!makeBarrelStateView(forward, barrelView));
  barrelView.setY(11.f);
  BOOST_CHECK_EQUAL(barrel.parameters[0], 11.f);
  BOOST_CHECK_EQUAL(forward.parameters[0], 0.f);

  ForwardStateView forwardView;
  BOOST_REQUIRE(makeForwardStateView(forward, forwardView));
  BOOST_CHECK(!makeForwardStateView(barrel, forwardView));
  forwardView.setX(13.f);
  BOOST_CHECK_EQUAL(forward.parameters[0], 13.f);
  BOOST_CHECK_EQUAL(barrel.parameters[0], 11.f);
}

BOOST_AUTO_TEST_CASE(BarrelAdapterHasExactDeclaredParityAndChargePidRoundTrip)
{
  const auto barrel = makeBarrelState();
  SurfaceKinematicState state{};
  BOOST_REQUIRE(legacy::importBarrelTrackParCov(barrel, state));
  BOOST_CHECK(state.kind == SurfaceKind::Cylinder);
  BOOST_CHECK_EQUAL(state.absCharge, 2U);
  BOOST_CHECK_EQUAL(state.pid.getID(), o2::track::PID::Kaon);
  o2::track::TrackParCovF restored{};
  BOOST_REQUIRE(legacy::exportBarrelTrackParCov(state, restored));
  BOOST_CHECK_EQUAL(restored.getX(), barrel.getX());
  BOOST_CHECK_EQUAL(restored.getAlpha(), barrel.getAlpha());
  BOOST_CHECK_EQUAL(restored.getAbsCharge(), barrel.getAbsCharge());
  BOOST_CHECK_EQUAL(restored.getPID().getID(), barrel.getPID().getID());
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_EQUAL(restored.getParam(i), barrel.getParam(i));
  }
  for (uint8_t i = 0; i < 15; ++i) {
    BOOST_CHECK_EQUAL(restored.getCov()[i], barrel.getCov()[i]);
  }

  const auto sentinel = restored;
  state.kind = SurfaceKind::Disk;
  BOOST_CHECK(!legacy::exportBarrelTrackParCov(state, restored));
  BOOST_CHECK(bitEqual(restored, sentinel));
}

BOOST_AUTO_TEST_CASE(ForwardAdapterNarrowsWithExplicitToleranceAndRejectsNonFiniteWithoutMutation)
{
  const auto forward = makeForwardState();
  SurfaceKinematicState state{};
  BOOST_REQUIRE(legacy::importLegacyForwardTrackParCov(forward, state));
  BOOST_CHECK(state.kind == SurfaceKind::Disk);
  BOOST_CHECK_EQUAL(state.alpha, 0.f);
  BOOST_CHECK_EQUAL(state.absCharge, 1U);
  BOOST_CHECK_EQUAL(state.pid.getID(), o2::track::PID::Pion);

  BOOST_CHECK(state.hasRecognizedKind());
  BOOST_CHECK_CLOSE_FRACTION(state.referenceCoordinate, static_cast<float>(forward.getZ()), 1.e-7);
  BOOST_CHECK_EQUAL(state.alpha, 0.f);
  BOOST_CHECK(state.kind == SurfaceKind::Disk);
  BOOST_CHECK_EQUAL(state.absCharge, 1U);
  BOOST_CHECK_EQUAL(state.pid.getID(), o2::track::PID::Pion);
  const double expectedParameters[] = {forward.getX(), forward.getY(), forward.getPhi(), forward.getTanl(), forward.getInvQPt()};
  for (uint8_t i = 0; i < 5; ++i) {
    BOOST_CHECK_CLOSE_FRACTION(state.parameters[i], static_cast<float>(expectedParameters[i]), 1.e-7);
  }
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      BOOST_CHECK_CLOSE_FRACTION(state.covariance[packedCovarianceIndex(row, column)],
                                 static_cast<float>(forward.getCovariances()(row, column)), 1.e-7);
    }
  }

  const SurfaceKinematicState sentinel = state;
  auto invalid = forward;
  invalid.setX(std::numeric_limits<double>::infinity());
  BOOST_CHECK(!legacy::importLegacyForwardTrackParCov(invalid, state));
  BOOST_CHECK(bitEqual(state, sentinel));

  invalid = forward;
  invalid.setZ(std::numeric_limits<double>::max());
  BOOST_CHECK(!legacy::importLegacyForwardTrackParCov(invalid, state));
  BOOST_CHECK(bitEqual(state, sentinel));

  invalid = forward;
  invalid.setX(std::numeric_limits<double>::max());
  BOOST_CHECK(!legacy::importLegacyForwardTrackParCov(invalid, state));
  BOOST_CHECK(bitEqual(state, sentinel));

  invalid = forward;
  auto covariance = invalid.getCovariances();
  covariance(3, 1) = std::numeric_limits<double>::max();
  invalid.setCovariances(covariance);
  BOOST_CHECK(!legacy::importLegacyForwardTrackParCov(invalid, state));
  BOOST_CHECK(bitEqual(state, sentinel));
}
