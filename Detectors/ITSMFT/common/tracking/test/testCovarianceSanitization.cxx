// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// M5d covariance-validity correction (doc/decisions/0008-native-refit-activation.md,
// covariance-fault-localization investigation): focused, deterministic
// regression coverage for sanitizeCovariance() (SurfaceKinematicState.h) and
// its eight call sites (barrel rotate/propagate x2 overloads/update, forward
// propagate<Model> x2 overloads/update). Several fixtures below reproduce a
// real captured production failure verbatim (exact state/covariance/
// measurement values from a checksummed replay of the
// pp-20ev-run303000-seed20260716-daily20260717 fixture, candidate keys
// "13,6,6,5,4,9,5" (ITS) and "68,71,73,67,72,73,62,76,80,-1" (MFT)) rather
// than a synthetic approximation, per the covariance-fault-localization
// investigation's minimal-reproducer design.

#define BOOST_TEST_MODULE ITSMFTCovarianceSanitization
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/SurfaceLinearizationReference.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

namespace
{
using namespace o2::itsmft::tracking;

bool allDiagonalsNonNegative(const SurfaceKinematicState& state)
{
  for (uint8_t i = 0; i < 5; ++i) {
    if (state.covariance[packedCovarianceIndex(i, i)] < 0.f) {
      return false;
    }
  }
  return true;
}

bool closeTo(float a, float b, float absTol = 5.e-4f, float relTol = 2.e-3f)
{
  const float diff = std::fabs(a - b);
  return diff <= absTol || diff <= relTol * std::fabs(b);
}

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}
} // namespace

// --- 1. sanitizeCovariance() itself: the core policy, in isolation. -------

BOOST_AUTO_TEST_CASE(SanitizeCovarianceAbsNegativeDiagonal)
{
  SurfaceKinematicState state{};
  state.covariance[packedCovarianceIndex(0, 0)] = -0.25f;
  state.covariance[packedCovarianceIndex(1, 1)] = 0.5f;
  state.covariance[packedCovarianceIndex(2, 2)] = 0.5f;
  state.covariance[packedCovarianceIndex(3, 3)] = 0.5f;
  state.covariance[packedCovarianceIndex(4, 4)] = 0.5f;
  const float maxDiagonal[5] = {1.f, 1.f, 1.f, 1.f, 1.f};
  sanitizeCovariance(state, maxDiagonal);
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(0, 0)], 0.25f, 1e-4f);
  BOOST_CHECK(allDiagonalsNonNegative(state));
}

BOOST_AUTO_TEST_CASE(SanitizeCovarianceClampsOverRangeAndRescalesOffDiagonal)
{
  SurfaceKinematicState state{};
  state.covariance[packedCovarianceIndex(0, 0)] = 4.f; // 4x the max below.
  state.covariance[packedCovarianceIndex(1, 0)] = 1.f; // Shares row/column 0.
  state.covariance[packedCovarianceIndex(2, 0)] = 2.f;
  state.covariance[packedCovarianceIndex(1, 1)] = 0.3f;
  state.covariance[packedCovarianceIndex(2, 2)] = 0.3f;
  state.covariance[packedCovarianceIndex(3, 3)] = 0.3f;
  state.covariance[packedCovarianceIndex(4, 4)] = 0.3f;
  const float maxDiagonal[5] = {1.f, 1.f, 1.f, 1.f, 1.f};
  sanitizeCovariance(state, maxDiagonal);
  // scale = sqrt(max/old) = sqrt(1/4) = 0.5.
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(0, 0)], 1.f, 1e-4f);
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(1, 0)], 0.5f, 1e-4f);
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(2, 0)], 1.f, 1e-4f);
  // Untouched entries not sharing the clamped row/column.
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(1, 1)], 0.3f, 1e-4f);
}

BOOST_AUTO_TEST_CASE(SanitizeCovariancePreservesSymmetryByConstruction)
{
  // Packed lower-triangular storage: only one entry exists per (row,column)
  // pair, so "symmetry" is a representation invariant, not a check --
  // packedCovarianceIndex(i,j) == packedCovarianceIndex(j,i) is exercised
  // directly by every read/write sanitizeCovariance performs.
  SurfaceKinematicState state{};
  state.covariance[packedCovarianceIndex(3, 1)] = 5.f;
  const float maxDiagonal[5] = {1.f, 1.f, 1.f, 1.f, 1.f};
  sanitizeCovariance(state, maxDiagonal);
  BOOST_CHECK_EQUAL(packedCovarianceIndex(1, 3), packedCovarianceIndex(3, 1));
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(1, 3)], 5.f, 1e-4f);
}

// --- 2. ITS legB reproducer: barrel::update() on the exact captured real --
// prior state/covariance and measurement (candidate "13,6,6,5,4,9,5", hit 5)
// that produced OperationFailureReason::MaterialFailure /
// MaterialFailureReason::InvalidCovariance before this correction (posterior
// Q2Pt-Q2Pt diagonal = -0.032802999, real production value, captured
// verbatim from the checksummed 20-event replay).

BOOST_AUTO_TEST_CASE(ITSLegBReproducerNowSanitizesToValidCovariance)
{
  SurfaceKinematicState state{};
  state.family = StateFamily::Barrel;
  state.referenceCoordinate = 3.76323366f;
  state.alpha = -0.12901926f;
  state.parameters[0] = 0.642236829f;
  state.parameters[1] = -6.11814785f;
  state.parameters[2] = 0.167980343f;
  state.parameters[3] = -1.58871007f;
  state.parameters[4] = 1.2842629f;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  const float cov[15] = {
    0.0615117364f, -0.0162716303f, 0.00781002454f, -0.00648284703f, 0.00164899346f, 0.000680086901f,
    0.000152464694f, -0.000262976653f, -1.178005e-05f, 1.45164713e-05f,
    -0.22814776f, 0.0546577908f, 0.0237723477f, -0.000194984852f, 0.822642863f};
  for (int i = 0; i < 15; ++i) {
    state.covariance[i] = cov[i];
  }

  SurfaceMeasurement meas{};
  meas.frame.u = 0.633100867f;
  meas.frame.v = -6.10807085f;
  meas.covariance.uu = 1.18710993e-07f;
  meas.covariance.uv = 0.f;
  meas.covariance.vv = 3.60069805e-07f;

  float chi2 = 0.f;
  OperationFailureReason reason{};
  const bool ok = barrel::update(state, meas, chi2, reason);

  BOOST_REQUIRE(ok);
  BOOST_CHECK(allDiagonalsNonNegative(state));
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(4, 4)], 0.0328048468f, 5.f); // sign-flipped, matches production magnitude within float tolerance.
}

// --- 3. MFT reproducer: forward::update() on the exact captured real ------
// prior state/covariance and measurement (candidate
// "68,71,73,67,72,73,62,76,80,-1", legB, hit 3) that produced a
// Q2Pt-Q2Pt diagonal of -52.064167 (real production value) before this
// correction.

BOOST_AUTO_TEST_CASE(MFTReproducerNowSanitizesToValidCovariance)
{
  SurfaceKinematicState state{};
  state.family = StateFamily::Forward;
  state.referenceCoordinate = -67.6889038f;
  state.alpha = 0.f;
  state.parameters[0] = -3.40663648f;
  state.parameters[1] = -3.04799104f;
  state.parameters[2] = -2.40926218f;
  state.parameters[3] = -15.2632132f;
  state.parameters[4] = -0.0805783421f;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  const float cov[15] = {
    5.45968469e-05f, 7.92145147e-05f, 2.34508334e-05f, 0.000361069426f, 9.60996113e-05f, 0.00121101411f,
    0.00140110496f, 0.0018511567f, 0.00489055132f, 0.0514357649f,
    0.117691882f, 0.0217776336f, 0.452787817f, 1.21933973f, 168.588654f};
  for (int i = 0; i < 15; ++i) {
    state.covariance[i] = cov[i];
  }

  SurfaceMeasurement meas{};
  meas.global.x = -3.4059999f;
  meas.global.y = -3.04678011f;
  meas.covariance.uu = 4.4239976e-05f;
  meas.covariance.uv = 0.f;
  meas.covariance.vv = 0.000105412393f;

  float chi2 = 0.f;
  OperationFailureReason reason{};
  const bool ok = forward::update(state, meas, chi2, reason);

  BOOST_REQUIRE(ok);
  BOOST_CHECK(allDiagonalsNonNegative(state));
}

// --- 4. Large-step propagation invariant: barrel::propagate(state, linRef, --
// ...) on the exact captured real inputs that fed the ITS legB reproducer
// above (the immediately preceding hit) must itself leave the covariance
// invariant satisfied before the next update() ever runs, even though the
// raw off-diagonal transport for this large (~-15.5cm) step is what first
// makes the matrix's Y-Q2Pt correlation exceed 1 (the precondition
// update() then reveals as a negative diagonal -- see the
// covariance-fault-localization investigation, ADR 0008).

BOOST_AUTO_TEST_CASE(LargeStepPropagationPreservesInvariantBeforeUpdate)
{
  SurfaceKinematicState state{};
  state.family = StateFamily::Barrel;
  state.referenceCoordinate = 19.2192478f;
  state.alpha = -0.12901926f;
  state.parameters[0] = 3.03678966f;
  state.parameters[1] = -30.9622726f;
  state.parameters[2] = 0.138186395f;
  state.parameters[3] = -1.58871007f;
  state.parameters[4] = 1.2842629f;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  const float cov[15] = {
    1.98605832e-07f, -7.50241043e-08f, 2.4906555e-07f, -9.13163856e-09f, -3.06560999e-08f, 1.98362322e-05f,
    3.80933152e-09f, -3.28565477e-08f, -7.25654581e-06f, 1.45164713e-05f,
    -1.76052566e-07f, -8.04973183e-07f, 0.00468764221f, -0.000194984852f, 0.822642863f};
  for (int i = 0; i < 15; ++i) {
    state.covariance[i] = cov[i];
  }

  SurfaceLinearizationReference linRef{};
  linRef.family = StateFamily::Barrel;
  linRef.referenceCoordinate = 19.2192478f;
  linRef.alpha = -0.12901926f;
  linRef.parameters[0] = 3.03678894f;
  linRef.parameters[1] = -30.962265f;
  linRef.parameters[2] = 0.137899101f;
  linRef.parameters[3] = -1.58717895f;
  linRef.parameters[4] = 1.21498108f;

  const float targetX = 3.76323366f;
  const float bz = 5.00675011f;
  OperationFailureReason reason{};
  const bool ok = barrel::propagate(state, linRef, targetX, bz, reason);

  BOOST_REQUIRE(ok);
  BOOST_CHECK(allDiagonalsNonNegative(state));
  // Cross-check against the real captured production post-propagate values
  // (also independently hand-rederived in the investigation, float32 and
  // float64, agreeing to 7 significant digits -- confirming this large-dx
  // Jacobian transport is precision-independent, not a rounding artifact).
  BOOST_CHECK(closeTo(state.covariance[packedCovarianceIndex(0, 0)], 0.0615117364f));
  BOOST_CHECK(closeTo(state.covariance[packedCovarianceIndex(4, 0)], -0.22814776f));
  BOOST_CHECK(closeTo(state.covariance[packedCovarianceIndex(4, 4)], 0.822642863f));
}

// --- 5. Every rotate/propagate/update independently sanitizes, both -------
// families. Each case below uses a deliberate zero-step (rotate: delta==0;
// propagate: dx/dz==0) or an otherwise-trivial transport so the operation's
// own transform is a documented no-op/identity on the covariance, isolating
// the sanitization call itself as the only thing that can explain a clamped
// result -- rather than depending on a from-scratch derivation of each
// operation's own Jacobian to predict a non-trivial expected output.

SurfaceKinematicState makeOverRangeBarrelState()
{
  SurfaceKinematicState state{};
  state.family = StateFamily::Barrel;
  state.referenceCoordinate = 4.f;
  state.alpha = 0.3f;
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.2f;
  state.parameters[3] = -0.35f;
  state.parameters[4] = 0.05f; // Small |Q2Pt| so Q2Pt-Q2Pt max isn't reached trivially by other tests.
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  state.covariance[packedCovarianceIndex(0, 0)] = 50.f * o2::track::kCY2max; // Deliberately over range.
  state.covariance[packedCovarianceIndex(1, 1)] = 0.01f;
  state.covariance[packedCovarianceIndex(2, 2)] = 0.01f;
  state.covariance[packedCovarianceIndex(3, 3)] = 0.01f;
  state.covariance[packedCovarianceIndex(4, 4)] = 0.01f;
  return state;
}

BOOST_AUTO_TEST_CASE(BarrelRotateSanitizesOnZeroDeltaTrivialStep)
{
  SurfaceKinematicState state = makeOverRangeBarrelState();
  OperationFailureReason reason{};
  const bool ok = barrel::rotate(state, state.alpha, reason); // delta == 0: ratio == 1, transform is identity.
  BOOST_REQUIRE(ok);
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(0, 0)], o2::track::kCY2max, 1e-3f);
}

BOOST_AUTO_TEST_CASE(BarrelPropagateSanitizesOnZeroDxTrivialStep)
{
  SurfaceKinematicState state = makeOverRangeBarrelState();
  OperationFailureReason reason{};
  const bool ok = barrel::propagate(state, state.referenceCoordinate, 0.5f, reason); // dx == 0: early-return path.
  BOOST_REQUIRE(ok);
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(0, 0)], o2::track::kCY2max, 1e-3f);
}

BOOST_AUTO_TEST_CASE(BarrelUpdateSanitizesReproducer)
{
  // Same fixture and assertion as ITSLegBReproducerNowSanitizesToValidCovariance
  // above; kept as a separate, minimally-named case so "update sanitizes" is
  // independently visible in the test list without relying on the reproducer
  // test's name to convey it.
  SurfaceKinematicState state{};
  state.family = StateFamily::Barrel;
  state.referenceCoordinate = 3.76323366f;
  state.alpha = -0.12901926f;
  state.parameters[0] = 0.642236829f;
  state.parameters[1] = -6.11814785f;
  state.parameters[2] = 0.167980343f;
  state.parameters[3] = -1.58871007f;
  state.parameters[4] = 1.2842629f;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  const float cov[15] = {
    0.0615117364f, -0.0162716303f, 0.00781002454f, -0.00648284703f, 0.00164899346f, 0.000680086901f,
    0.000152464694f, -0.000262976653f, -1.178005e-05f, 1.45164713e-05f,
    -0.22814776f, 0.0546577908f, 0.0237723477f, -0.000194984852f, 0.822642863f};
  for (int i = 0; i < 15; ++i) {
    state.covariance[i] = cov[i];
  }
  SurfaceMeasurement meas{};
  meas.frame.u = 0.633100867f;
  meas.frame.v = -6.10807085f;
  meas.covariance.uu = 1.18710993e-07f;
  meas.covariance.vv = 3.60069805e-07f;
  float chi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(barrel::update(state, meas, chi2, reason));
  BOOST_CHECK(allDiagonalsNonNegative(state));
}

BOOST_AUTO_TEST_CASE(BarrelLinRefRotateSanitizesOnZeroDeltaTrivialStep)
{
  SurfaceKinematicState state = makeOverRangeBarrelState();
  SurfaceLinearizationReference linRef{};
  linRef.family = StateFamily::Barrel;
  linRef.referenceCoordinate = state.referenceCoordinate;
  linRef.alpha = state.alpha;
  for (int i = 0; i < 5; ++i) {
    linRef.parameters[i] = state.parameters[i];
  }
  OperationFailureReason reason{};
  const bool ok = barrel::rotate(state, linRef, state.alpha, 0.5f, reason);
  BOOST_REQUIRE(ok);
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(0, 0)], o2::track::kCY2max, 1e-3f);
}

BOOST_AUTO_TEST_CASE(BarrelLinRefPropagateSanitizesLargeStep)
{
  // Same fixture and assertion as LargeStepPropagationPreservesInvariantBeforeUpdate
  // above; kept as a separate, minimally-named case for the same reason as
  // BarrelUpdateSanitizesReproducer.
  SurfaceKinematicState state{};
  state.family = StateFamily::Barrel;
  state.referenceCoordinate = 19.2192478f;
  state.alpha = -0.12901926f;
  state.parameters[0] = 3.03678966f;
  state.parameters[1] = -30.9622726f;
  state.parameters[2] = 0.138186395f;
  state.parameters[3] = -1.58871007f;
  state.parameters[4] = 1.2842629f;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  const float cov[15] = {
    1.98605832e-07f, -7.50241043e-08f, 2.4906555e-07f, -9.13163856e-09f, -3.06560999e-08f, 1.98362322e-05f,
    3.80933152e-09f, -3.28565477e-08f, -7.25654581e-06f, 1.45164713e-05f,
    -1.76052566e-07f, -8.04973183e-07f, 0.00468764221f, -0.000194984852f, 0.822642863f};
  for (int i = 0; i < 15; ++i) {
    state.covariance[i] = cov[i];
  }
  SurfaceLinearizationReference linRef{};
  linRef.family = StateFamily::Barrel;
  linRef.referenceCoordinate = 19.2192478f;
  linRef.alpha = -0.12901926f;
  linRef.parameters[0] = 3.03678894f;
  linRef.parameters[1] = -30.962265f;
  linRef.parameters[2] = 0.137899101f;
  linRef.parameters[3] = -1.58717895f;
  linRef.parameters[4] = 1.21498108f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(barrel::propagate(state, linRef, 3.76323366f, 5.00675011f, reason));
  BOOST_CHECK(allDiagonalsNonNegative(state));
}

SurfaceKinematicState makeOverRangeForwardState()
{
  SurfaceKinematicState state{};
  state.family = StateFamily::Forward;
  state.referenceCoordinate = -40.f;
  state.alpha = 0.f;
  state.parameters[0] = 1.f;
  state.parameters[1] = -1.f;
  state.parameters[2] = 0.1f;
  state.parameters[3] = -2.f;
  state.parameters[4] = 0.05f;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  state.covariance[packedCovarianceIndex(0, 0)] = 50.f * o2::track::kCY2max; // Deliberately over range.
  state.covariance[packedCovarianceIndex(1, 1)] = 0.01f;
  state.covariance[packedCovarianceIndex(2, 2)] = 0.01f;
  state.covariance[packedCovarianceIndex(3, 3)] = 0.01f;
  state.covariance[packedCovarianceIndex(4, 4)] = 0.01f;
  return state;
}

BOOST_AUTO_TEST_CASE(ForwardPropagateSanitizesOnZeroDzTrivialStep)
{
  SurfaceKinematicState state = makeOverRangeForwardState();
  OperationFailureReason reason{};
  const bool ok = forward::propagate<forward::PropagationModel::Linear>(state, state.referenceCoordinate, 0.5f, reason);
  BOOST_REQUIRE(ok);
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(0, 0)], o2::track::kCY2max, 1e-3f);
}

BOOST_AUTO_TEST_CASE(ForwardLinRefPropagateSanitizesOnZeroDzTrivialStep)
{
  SurfaceKinematicState state = makeOverRangeForwardState();
  SurfaceLinearizationReference linRef{};
  linRef.family = StateFamily::Forward;
  linRef.referenceCoordinate = state.referenceCoordinate;
  for (int i = 0; i < 5; ++i) {
    linRef.parameters[i] = state.parameters[i];
  }
  OperationFailureReason reason{};
  const bool ok = forward::propagate<forward::PropagationModel::Linear>(state, linRef, state.referenceCoordinate, 0.5f, reason);
  BOOST_REQUIRE(ok);
  BOOST_CHECK_CLOSE(state.covariance[packedCovarianceIndex(0, 0)], o2::track::kCY2max, 1e-3f);
}

BOOST_AUTO_TEST_CASE(ForwardUpdateSanitizesReproducer)
{
  SurfaceKinematicState state{};
  state.family = StateFamily::Forward;
  state.referenceCoordinate = -67.6889038f;
  state.alpha = 0.f;
  state.parameters[0] = -3.40663648f;
  state.parameters[1] = -3.04799104f;
  state.parameters[2] = -2.40926218f;
  state.parameters[3] = -15.2632132f;
  state.parameters[4] = -0.0805783421f;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  const float cov[15] = {
    5.45968469e-05f, 7.92145147e-05f, 2.34508334e-05f, 0.000361069426f, 9.60996113e-05f, 0.00121101411f,
    0.00140110496f, 0.0018511567f, 0.00489055132f, 0.0514357649f,
    0.117691882f, 0.0217776336f, 0.452787817f, 1.21933973f, 168.588654f};
  for (int i = 0; i < 15; ++i) {
    state.covariance[i] = cov[i];
  }
  SurfaceMeasurement meas{};
  meas.global.x = -3.4059999f;
  meas.global.y = -3.04678011f;
  meas.covariance.uu = 4.4239976e-05f;
  meas.covariance.vv = 0.000105412393f;
  float chi2 = 0.f;
  OperationFailureReason reason{};
  BOOST_REQUIRE(forward::update(state, meas, chi2, reason));
  BOOST_CHECK(allDiagonalsNonNegative(state));
}

// --- 6. preflightValidate remains strict: a deliberately malformed --------
// *externally supplied* state (never touched by propagate/rotate/update) is
// still rejected by correctForMaterial's preflight, proving the fix does not
// weaken or bypass that check -- it only ensures the Propagator's own
// internal callers never hand it an invalid state in normal operation.

BOOST_AUTO_TEST_CASE(MalformedExternalBarrelCovarianceStillRejectedByPreflight)
{
  SurfaceKinematicState state{};
  state.family = StateFamily::Barrel;
  state.referenceCoordinate = 4.f;
  state.alpha = 0.3f;
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.2f;
  state.parameters[3] = -0.35f;
  state.parameters[4] = 0.8f;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  for (uint8_t i = 0; i < 5; ++i) {
    state.covariance[packedCovarianceIndex(i, i)] = 0.01f;
  }
  state.covariance[packedCovarianceIndex(4, 4)] = -0.01f; // Deliberately invalid, constructed directly.

  const material::IntegratedMaterialBudget budget{0.01f, 0.05f};
  const auto result = barrel::correctForMaterial(state, budget, material::MaterialTraversalDirection::AlongMomentum);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.failure == material::MaterialFailureReason::InvalidCovariance);
}

BOOST_AUTO_TEST_CASE(MalformedExternalForwardCovarianceStillRejectedByPreflight)
{
  SurfaceKinematicState state{};
  state.family = StateFamily::Forward;
  state.referenceCoordinate = -40.f;
  state.alpha = 0.f;
  state.parameters[0] = 1.f;
  state.parameters[1] = -1.f;
  state.parameters[2] = 0.1f;
  state.parameters[3] = -2.f;
  state.parameters[4] = 0.05f;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  for (uint8_t i = 0; i < 5; ++i) {
    state.covariance[packedCovarianceIndex(i, i)] = 0.01f;
  }
  state.covariance[packedCovarianceIndex(2, 2)] = -0.01f; // Deliberately invalid, constructed directly.

  const material::IntegratedMaterialBudget budget{0.01f, 0.05f};
  const auto result = forward::correctForMaterial(state, budget, material::MaterialTraversalDirection::AlongMomentum);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.failure == material::MaterialFailureReason::InvalidCovariance);
}

// --- 7. Operation failure remains transactional: a failing update/rotate/-
// propagate call must leave the input state byte-for-byte unchanged -- the
// new sanitization call must never run (and never partially mutate state)
// on a failure path.

BOOST_AUTO_TEST_CASE(FailingBarrelUpdateLeavesStateUnchanged)
{
  SurfaceKinematicState state{};
  state.family = StateFamily::Barrel;
  state.referenceCoordinate = 4.f;
  state.alpha = 0.3f;
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.2f;
  state.parameters[3] = -0.35f;
  state.parameters[4] = 0.8f;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  for (uint8_t i = 0; i < 5; ++i) {
    state.covariance[packedCovarianceIndex(i, i)] = 0.01f;
  }
  const SurfaceKinematicState original = state;

  SurfaceMeasurement meas{};
  meas.frame.u = std::numeric_limits<float>::quiet_NaN(); // Non-finite measurement: must fail before any mutation.
  meas.frame.v = 0.f;
  meas.covariance = {0.01f, 0.f, 0.01f};
  float chi2 = 0.f;
  OperationFailureReason reason{};
  const bool ok = barrel::update(state, meas, chi2, reason);

  BOOST_CHECK(!ok);
  BOOST_CHECK(bitEqual(state, original));
}

BOOST_AUTO_TEST_CASE(FailingBarrelRotateLeavesStateUnchanged)
{
  SurfaceKinematicState state{};
  state.family = StateFamily::Barrel;
  state.referenceCoordinate = 4.f;
  state.alpha = 0.3f;
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 1.5f; // |Snp| >= 1: rotate must reject before touching anything.
  state.parameters[3] = -0.35f;
  state.parameters[4] = 0.8f;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  for (uint8_t i = 0; i < 5; ++i) {
    state.covariance[packedCovarianceIndex(i, i)] = 0.01f;
  }
  const SurfaceKinematicState original = state;

  OperationFailureReason reason{};
  const bool ok = barrel::rotate(state, state.alpha + 3.0f, reason); // Large rotation: local direction inversion.

  BOOST_CHECK(!ok);
  BOOST_CHECK(bitEqual(state, original));
}
