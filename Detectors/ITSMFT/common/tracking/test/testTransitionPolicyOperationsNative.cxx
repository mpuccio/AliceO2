// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT TransitionPolicyOperationsNative
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <cmath>
#include <cstring>
#include <limits>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/BarrelSurfaceStateOperations.h"
#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/TransitionPolicyOperations.h"
#include "ITStracking/Cluster.h"

/// Stage-B Slice B focused coverage: the native SurfaceKinematicState
/// overloads of buildCellSeed<Tag>/cellsAreCompatible<Tag>/attachHit<Tag>
/// (TransitionPolicyOperations.h/.cxx). These operations are additive and
/// unwired -- there is no production caller and no TrackerTraits/TimeFrame
/// orchestration exercised here.
///
/// Oracle strategy (see the Stage-B design report's "Precision/oracles"
/// split, referenced from TransitionPolicyOperations.h):
///  - Formula-preserving evidence (seed construction, rotation/propagation,
///    predicted chi2/update, state compatibility) is checked by an
///    independent, hand-written re-transcription of each operation's own
///    documented call sequence, built directly on the already-oracle-tested
///    barrel::/forward:: primitives (BarrelSurfaceStateOperations.h,
///    ForwardSurfaceStateOperations.h -- each already characterized against
///    its own legacy oracle in testBarrelSurfaceStateOperations.cxx /
///    testForwardSurfaceStateOperations.cxx). A bit-identical match against
///    this independent replay is strong evidence that the production
///    orchestration (step order, material-slot selection, chi2-cut
///    placement, measurement projection) is correct, without re-deriving
///    the already-tested primitive formulas themselves.
///  - Intentional material-physics differences (PID/absCharge-aware barrel
///    covariance correction; newly active forward energy loss/straggling)
///    are validated by observing that charge/PID and areal-density inputs
///    measurably change the result, not by comparing to legacy MCS-only
///    output.
using namespace o2::itsmft::tracking;

namespace
{

template <typename T>
bool bitEqual(const T& lhs, const T& rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

// --- Barrel fixtures -------------------------------------------------------
// Reuses the same well-separated, non-collinear geometry already established
// in testBarrelSurfaceStateOperations.cxx's buildSeed fixtures, extended
// with two inward attach hits sharing the outer hit's frame angle (delta=0
// rotation) so the nominal-success fixtures below stay clear of rotation's
// own cos(phi)>=0 boundary.

// Small transverse offsets relative to radius (nearly straight, realistic
// high-pT trajectory) -- unlike testBarrelSurfaceStateOperations.cxx's own
// buildSeed-only fixture (which deliberately uses larger offsets that are
// fine for a finiteness/oracle-only check but yield an extreme, sub-10-MeV
// curvature unsuited to exercising material correction here).
o2::its::Cluster barrelClusterInner() { return o2::its::Cluster{1.8f, 0.0010f, -0.60f, 0}; }
o2::its::Cluster barrelClusterMiddle() { return o2::its::Cluster{3.0f, 0.0035f, 0.10f, 1}; }
o2::its::Cluster barrelClusterOuter() { return o2::its::Cluster{4.4f, 0.0050f, -0.8f, 2}; }

o2::its::TrackingFrameInfo barrelHitOuter()
{
  return o2::its::TrackingFrameInfo{0.f, 0.f, 0.f, 4.4f, 0.f, {0.0050f, -0.8f}, {0.02f, 0.001f, 0.03f}};
}
o2::its::TrackingFrameInfo barrelHitMiddle()
{
  return o2::its::TrackingFrameInfo{0.f, 0.f, 0.f, 3.0f, 0.f, {0.0035f, 0.10f}, {0.02f, 0.001f, 0.03f}};
}
o2::its::TrackingFrameInfo barrelHitInner()
{
  return o2::its::TrackingFrameInfo{0.f, 0.f, 0.f, 1.8f, 0.f, {0.0010f, -0.60f}, {0.02f, 0.001f, 0.03f}};
}

constexpr float BarrelBz = 0.5f;

NominalSurfaceMaterial barrelMaterialMiddle() { return NominalSurfaceMaterial{0.02f, 0.002f}; }
NominalSurfaceMaterial barrelMaterialInner() { return NominalSurfaceMaterial{0.01f, 0.001f}; }
// Deliberately invalid: buildCellSeed<CylinderCylinder> must never read the
// outer material slot (the outer surface only contributes through hitOuter
// inside barrel::buildSeed, matching the legacy operation exactly).
NominalSurfaceMaterial barrelMaterialOuterPoison()
{
  const float nan = std::numeric_limits<float>::quiet_NaN();
  return NominalSurfaceMaterial{nan, nan};
}

std::array<NominalSurfaceMaterial, 3> barrelMaterial()
{
  return {barrelMaterialInner(), barrelMaterialMiddle(), barrelMaterialOuterPoison()};
}

SurfaceMeasurement barrelMeasurementFromHit(const o2::its::TrackingFrameInfo& hit)
{
  SurfaceMeasurement measurement{};
  measurement.frame.u = hit.positionTrackingFrame[0];
  measurement.frame.v = hit.positionTrackingFrame[1];
  measurement.covariance.uu = hit.covarianceTrackingFrame[0];
  measurement.covariance.uv = hit.covarianceTrackingFrame[1];
  measurement.covariance.vv = hit.covarianceTrackingFrame[2];
  return measurement;
}

// Independent re-transcription of native buildCellSeed<CylinderCylinder>'s
// own documented call sequence, built directly on the public barrel::
// primitives rather than calling the production function under test.
// `nSteps` limits how many of the middle/inner attach steps run (1 = middle
// only, used to read out the inner step's own marginal predicted chi2 for
// boundary tests), mirroring the nSteps convention already established by
// legacyBarrelCellFit in testTransitionPolicyOperations.cxx.
bool replayBuildCellSeedBarrel(uint8_t absCharge, o2::track::PID pid, float maxChi2,
                               SurfaceKinematicState& outState, float& outChi2, float& lastStepChi2,
                               OperationFailureReason& reason, int nSteps = 2)
{
  SurfaceKinematicState state{};
  if (!barrel::buildSeed(barrelClusterInner(), barrelClusterMiddle(), barrelHitOuter(), BarrelBz,
                         absCharge, pid, state, reason)) {
    return false;
  }
  const std::array<o2::its::TrackingFrameInfo, 2> steps{barrelHitMiddle(), barrelHitInner()};
  const std::array<NominalSurfaceMaterial, 2> stepMaterial{barrelMaterialMiddle(), barrelMaterialInner()};
  float chi2 = 0.f;
  float stepChi2 = 0.f;
  for (int step = 0; step < nSteps; ++step) {
    const bool isLast = (nSteps == 2) && (step == 1);
    const auto& hit = steps[step];
    if (!barrel::rotate(state, hit.alphaTrackingFrame, reason)) {
      return false;
    }
    if (!barrel::propagate(state, hit.xTrackingFrame, BarrelBz, reason)) {
      return false;
    }
    const auto& stepMat = stepMaterial[step];
    const auto materialResult = barrel::correctForMaterial(
      state, material::IntegratedMaterialBudget{stepMat.xOverX0, stepMat.arealDensityGPerCm2},
      material::MaterialTraversalDirection::OppositeMomentum);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    const SurfaceMeasurement measurement = barrelMeasurementFromHit(hit);
    float predChi2 = 0.f;
    if (!barrel::predictedChi2(state, measurement, predChi2, reason)) {
      return false;
    }
    stepChi2 = predChi2;
    if (isLast && predChi2 > maxChi2) {
      reason = OperationFailureReason::PredictedChi2Failure;
      return false;
    }
    float updateChi2 = 0.f;
    if (!barrel::update(state, measurement, updateChi2, reason)) {
      return false;
    }
    chi2 += updateChi2;
  }
  outState = state;
  outChi2 = chi2;
  lastStepChi2 = stepChi2;
  return true;
}

// --- Disk fixtures -----------------------------------------------------
// Reuses testForwardSurfaceStateOperations.cxx's buildSeed geometry: inner z
// greater (less negative) than outer z, well-separated transverse positions.

o2::its::Cluster diskClusterInner() { return o2::its::Cluster{2.0f, 1.0f, -40.0f, 0}; }
o2::its::Cluster diskClusterMiddle() { return o2::its::Cluster{2.5f, 1.3f, -55.0f, 1}; }
o2::its::Cluster diskClusterOuter() { return o2::its::Cluster{3.0f, 1.6f, -70.0f, 2}; }

o2::its::TrackingFrameInfo diskHitOuter()
{
  return o2::its::TrackingFrameInfo{3.0f, 1.6f, -70.0f, 0.f, 0.f, {0.f, 0.f}, {0.05f, 0.f, 0.07f}};
}
o2::its::TrackingFrameInfo diskHitMiddle()
{
  return o2::its::TrackingFrameInfo{2.5f, 1.3f, -55.0f, 0.f, 0.f, {0.f, 0.f}, {0.05f, 0.f, 0.07f}};
}
o2::its::TrackingFrameInfo diskHitInner()
{
  return o2::its::TrackingFrameInfo{2.0f, 1.0f, -40.0f, 0.f, 0.f, {0.f, 0.f}, {0.05f, 0.f, 0.07f}};
}

constexpr float DiskBz = 5.f;
constexpr float DiskTrackletMinPt = 0.3f;

NominalSurfaceMaterial diskMaterialInner() { return NominalSurfaceMaterial{0.02f, 0.002f}; }
NominalSurfaceMaterial diskMaterialMiddle() { return NominalSurfaceMaterial{0.015f, 0.0015f}; }
NominalSurfaceMaterial diskMaterialOuter() { return NominalSurfaceMaterial{0.01f, 0.001f}; }

std::array<NominalSurfaceMaterial, 3> diskMaterial()
{
  return {diskMaterialInner(), diskMaterialMiddle(), diskMaterialOuter()};
}

SurfaceMeasurement diskMeasurementFromHit(const o2::its::TrackingFrameInfo& hit)
{
  SurfaceMeasurement measurement{};
  measurement.global.x = hit.xCoordinate;
  measurement.global.y = hit.yCoordinate;
  measurement.covariance.uu = hit.covarianceTrackingFrame[0];
  measurement.covariance.uv = 0.f;
  measurement.covariance.vv = hit.covarianceTrackingFrame[2];
  return measurement;
}

// "The accepted forward model": reproduces detail::mftFwdPropagateToZ's own
// |bz|>0.01f threshold dispatch to Helix/Linear (see
// TransitionPolicyOperations.cxx's forwardPropagateAcceptedModel, which is
// private to that translation unit -- this is an independent, deliberately
// duplicated one-line re-transcription, not a call into production code).
bool replayForwardPropagateAcceptedModel(SurfaceKinematicState& state, float targetZ, float bz,
                                         OperationFailureReason& reason)
{
  if (std::abs(bz) > 0.01f) {
    return forward::propagate<forward::PropagationModel::Helix>(state, targetZ, bz, reason);
  }
  return forward::propagate<forward::PropagationModel::Linear>(state, targetZ, bz, reason);
}

// Independent re-transcription of native buildCellSeed<DiskDisk>'s own
// documented call sequence (outer, middle, inner; material[2]/[1]/[0]).
bool replayBuildCellSeedDisk(uint8_t absCharge, o2::track::PID pid, float maxChi2,
                             SurfaceKinematicState& outState, float& outChi2, float& lastStepChi2,
                             OperationFailureReason& reason, int nSteps = 3)
{
  SurfaceKinematicState state{};
  if (!forward::buildSeed(diskClusterInner(), diskClusterMiddle(), diskClusterOuter(), diskHitOuter(),
                          DiskBz, DiskTrackletMinPt, absCharge, pid, state, reason)) {
    return false;
  }
  const std::array<o2::its::TrackingFrameInfo, 3> steps{diskHitOuter(), diskHitMiddle(), diskHitInner()};
  const std::array<NominalSurfaceMaterial, 3> stepMaterial{diskMaterialOuter(), diskMaterialMiddle(), diskMaterialInner()};
  float chi2 = 0.f;
  float stepChi2 = 0.f;
  for (int step = 0; step < nSteps; ++step) {
    const bool isLast = (nSteps == 3) && (step == 2);
    const auto& hit = steps[step];
    if (!replayForwardPropagateAcceptedModel(state, hit.zCoordinate, DiskBz, reason)) {
      return false;
    }
    const auto& stepMat = stepMaterial[step];
    const auto materialResult = forward::correctForMaterial(
      state, material::IntegratedMaterialBudget{stepMat.xOverX0, stepMat.arealDensityGPerCm2},
      material::MaterialTraversalDirection::OppositeMomentum);
    if (!materialResult.ok()) {
      reason = OperationFailureReason::MaterialFailure;
      return false;
    }
    const SurfaceMeasurement measurement = diskMeasurementFromHit(hit);
    float predChi2 = 0.f;
    if (!forward::predictedChi2(state, measurement, predChi2, reason)) {
      return false;
    }
    stepChi2 = predChi2;
    if (isLast && predChi2 > maxChi2) {
      reason = OperationFailureReason::PredictedChi2Failure;
      return false;
    }
    float updateChi2 = 0.f;
    if (!forward::update(state, measurement, updateChi2, reason)) {
      return false;
    }
    chi2 += updateChi2;
  }
  outState = state;
  outChi2 = chi2;
  lastStepChi2 = stepChi2;
  return true;
}

// --- attachHit fixtures -----------------------------------------------

SurfaceKinematicState barrelAttachState()
{
  SurfaceKinematicState state{};
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.2f;
  state.parameters[3] = -0.35f;
  state.parameters[4] = 0.8f;
  state.referenceCoordinate = 4.f;
  state.alpha = 0.3f;
  state.family = StateFamily::Barrel;
  state.absCharge = 1;
  state.pid = o2::track::PID::Kaon;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
    }
  }
  return state;
}

o2::its::TrackingFrameInfo barrelAttachHit()
{
  return o2::its::TrackingFrameInfo{0.f, 0.f, 0.f, 2.5f, 0.3f, {0.8f, -0.45f}, {0.04f, 0.012f, 0.09f}};
}

constexpr float BarrelAttachBz = 5.f;

NominalSurfaceMaterial barrelAttachMaterial() { return NominalSurfaceMaterial{0.01f, 0.001f}; }

SurfaceKinematicState diskAttachState()
{
  SurfaceKinematicState state{};
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.35f;
  state.parameters[3] = -2.5f;
  state.parameters[4] = 0.8f;
  state.referenceCoordinate = -45.f;
  state.family = StateFamily::Forward;
  state.absCharge = 2;
  state.pid = o2::track::PID::Pion;
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      state.covariance[packedCovarianceIndex(row, column)] = row == column ? 0.01f * (row + 1) : 0.0002f * (row + column + 1);
    }
  }
  return state;
}

o2::its::TrackingFrameInfo diskAttachHit()
{
  return o2::its::TrackingFrameInfo{0.8f, -0.45f, -50.f, 0.f, 0.f, {0.f, 0.f}, {0.04f, 0.f, 0.09f}};
}

constexpr float DiskAttachBz = 5.f;

NominalSurfaceMaterial diskAttachMaterial() { return NominalSurfaceMaterial{0.02f, 0.002f}; }

} // namespace

// ===========================================================================
// buildCellSeed<CylinderCylinder>
// ===========================================================================

BOOST_AUTO_TEST_CASE(BuildCellSeedBarrelSuccessMatchesIndependentReplayAndCharacterizesResult)
{
  SurfaceKinematicState replayState{};
  float replayChi2 = 0.f;
  float lastStepChi2 = 0.f;
  OperationFailureReason replayReason{};
  BOOST_REQUIRE(replayBuildCellSeedBarrel(1, o2::track::PID::Kaon, 1.e6f, replayState, replayChi2, lastStepChi2, replayReason));

  SurfaceKinematicState outState{};
  float chi2 = 0.f;
  OperationFailureReason reason{};
  CylinderCylinderPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::CylinderCylinder>(
    barrelClusterInner(), barrelClusterMiddle(), barrelClusterOuter(),
    barrelHitInner(), barrelHitMiddle(), barrelHitOuter(), barrelMaterial(), BarrelBz,
    1, o2::track::PID::Kaon, outState, chi2, params, reason));

  BOOST_CHECK(bitEqual(outState, replayState));
  BOOST_CHECK_EQUAL(chi2, replayChi2);
  BOOST_CHECK(std::isfinite(chi2));
  BOOST_CHECK_GT(chi2, 0.f);

  // Finishes anchored at the innermost hit's own frame.
  BOOST_CHECK_EQUAL(outState.referenceCoordinate, barrelHitInner().xTrackingFrame);
  BOOST_CHECK_CLOSE(outState.alpha, barrelHitInner().alphaTrackingFrame, 1.e-3f);

  BOOST_CHECK(outState.family == StateFamily::Barrel);
  BOOST_CHECK_EQUAL(outState.flags, uint8_t{0});
  BOOST_CHECK_EQUAL(outState.absCharge, uint8_t{1});
  BOOST_CHECK(outState.pid == o2::track::PID::Kaon);
}

BOOST_AUTO_TEST_CASE(BuildCellSeedBarrelExactChi2ThresholdBoundary)
{
  SurfaceKinematicState replayState{};
  float replayChi2 = 0.f;
  float lastStepChi2 = 0.f;
  OperationFailureReason replayReason{};
  BOOST_REQUIRE(replayBuildCellSeedBarrel(1, o2::track::PID::Pion, 1.e6f, replayState, replayChi2, lastStepChi2, replayReason));
  BOOST_REQUIRE_GT(lastStepChi2, 0.f);

  SurfaceKinematicState outState{};
  float chi2 = 0.f;
  OperationFailureReason reason{};
  CylinderCylinderPolicyParams accept;
  accept.maxChi2ClusterAttachment = lastStepChi2; // exact threshold: inclusive accept
  BOOST_CHECK(buildCellSeed<TransitionPolicyTag::CylinderCylinder>(
    barrelClusterInner(), barrelClusterMiddle(), barrelClusterOuter(),
    barrelHitInner(), barrelHitMiddle(), barrelHitOuter(), barrelMaterial(), BarrelBz,
    1, o2::track::PID::Pion, outState, chi2, accept, reason));

  const auto sentinel = barrelAttachState();
  auto before = sentinel;
  float chi2Before = -123.f;
  auto rejectedState = before;
  float rejectedChi2 = chi2Before;
  CylinderCylinderPolicyParams reject;
  reject.maxChi2ClusterAttachment = std::nextafter(lastStepChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::CylinderCylinder>(
    barrelClusterInner(), barrelClusterMiddle(), barrelClusterOuter(),
    barrelHitInner(), barrelHitMiddle(), barrelHitOuter(), barrelMaterial(), BarrelBz,
    1, o2::track::PID::Pion, rejectedState, rejectedChi2, reject, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(rejectedState, before));
  BOOST_CHECK_EQUAL(rejectedChi2, chi2Before);
}

BOOST_AUTO_TEST_CASE(BuildCellSeedBarrelFailuresAreTransactional)
{
  const auto sentinel = barrelAttachState();

  auto checkFailure = [&](const o2::its::Cluster& clusterInner, const o2::its::TrackingFrameInfo& hitMiddle,
                          const std::array<NominalSurfaceMaterial, 3>& material, OperationFailureReason expected) {
    auto state = sentinel;
    float chi2 = -321.f;
    const auto before = state;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    CylinderCylinderPolicyParams params;
    params.maxChi2ClusterAttachment = 1.e6f;
    BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::CylinderCylinder>(
      clusterInner, barrelClusterMiddle(), barrelClusterOuter(), barrelHitInner(), hitMiddle, barrelHitOuter(),
      material, BarrelBz, 1, o2::track::PID::Pion, state, chi2, params, reason));
    BOOST_CHECK(reason == expected);
    BOOST_CHECK(bitEqual(state, before));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  };

  // Seed-construction failure: non-finite cluster coordinate.
  {
    auto badInner = barrelClusterInner();
    badInner.xCoordinate = std::numeric_limits<float>::quiet_NaN();
    checkFailure(badInner, barrelHitMiddle(), barrelMaterial(), OperationFailureReason::NonFiniteInput);
  }

  // Rotation failure: hitMiddle's frame angle is far enough from the seed's
  // outer-anchored alpha that rotate's cos(phi)>=0 precondition fails
  // (mirrors testTransitionPolicyOperations.cxx's own rotation-failure
  // margin).
  {
    auto farHit = barrelHitMiddle();
    farHit.alphaTrackingFrame = barrelHitOuter().alphaTrackingFrame + 3.f;
    checkFailure(barrelClusterInner(), farHit, barrelMaterial(), OperationFailureReason::RotationFailure);
  }

  // Propagation failure: an unreachable target x far outside the seed's
  // curvature-bounded range.
  {
    auto farHit = barrelHitMiddle();
    farHit.xTrackingFrame = -50000.f;
    auto state = sentinel;
    float chi2 = -321.f;
    const auto before = state;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    CylinderCylinderPolicyParams params;
    params.maxChi2ClusterAttachment = 1.e6f;
    BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::CylinderCylinder>(
      barrelClusterInner(), barrelClusterMiddle(), barrelClusterOuter(), barrelHitInner(), farHit, barrelHitOuter(),
      barrelMaterial(), BarrelBz, 1, o2::track::PID::Pion, state, chi2, params, reason));
    BOOST_CHECK(reason == OperationFailureReason::UnreachableTarget || reason == OperationFailureReason::PropagationFailure);
    BOOST_CHECK(bitEqual(state, before));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  }

  // Material failure: non-finite middle material.
  {
    auto badMaterial = barrelMaterial();
    badMaterial[1] = NominalSurfaceMaterial{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
    checkFailure(barrelClusterInner(), barrelHitMiddle(), badMaterial, OperationFailureReason::MaterialFailure);
  }
}

BOOST_AUTO_TEST_CASE(BuildCellSeedBarrelIsByteDeterministic)
{
  SurfaceKinematicState first{};
  SurfaceKinematicState second{};
  float chi2First = 0.f;
  float chi2Second = 0.f;
  OperationFailureReason reason{};
  CylinderCylinderPolicyParams params;
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::CylinderCylinder>(
    barrelClusterInner(), barrelClusterMiddle(), barrelClusterOuter(), barrelHitInner(), barrelHitMiddle(),
    barrelHitOuter(), barrelMaterial(), BarrelBz, 1, o2::track::PID::Kaon, first, chi2First, params, reason));
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::CylinderCylinder>(
    barrelClusterInner(), barrelClusterMiddle(), barrelClusterOuter(), barrelHitInner(), barrelHitMiddle(),
    barrelHitOuter(), barrelMaterial(), BarrelBz, 1, o2::track::PID::Kaon, second, chi2Second, params, reason));
  BOOST_CHECK(bitEqual(first, second));
  BOOST_CHECK_EQUAL(chi2First, chi2Second);
}

// ===========================================================================
// buildCellSeed<DiskDisk>
// ===========================================================================

BOOST_AUTO_TEST_CASE(BuildCellSeedDiskSuccessMatchesIndependentReplayAndCharacterizesResult)
{
  SurfaceKinematicState replayState{};
  float replayChi2 = 0.f;
  float lastStepChi2 = 0.f;
  OperationFailureReason replayReason{};
  BOOST_REQUIRE(replayBuildCellSeedDisk(2, o2::track::PID::Pion, 1.e6f, replayState, replayChi2, lastStepChi2, replayReason));

  SurfaceKinematicState outState{};
  float chi2 = 0.f;
  OperationFailureReason reason{};
  DiskDiskPolicyParams params;
  params.trackletMinPt = DiskTrackletMinPt;
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::DiskDisk>(
    diskClusterInner(), diskClusterMiddle(), diskClusterOuter(),
    diskHitInner(), diskHitMiddle(), diskHitOuter(), diskMaterial(), DiskBz,
    2, o2::track::PID::Pion, outState, chi2, params, reason));

  BOOST_CHECK(bitEqual(outState, replayState));
  BOOST_CHECK_EQUAL(chi2, replayChi2);
  BOOST_CHECK(std::isfinite(chi2));
  BOOST_CHECK_GT(chi2, 0.f);

  BOOST_CHECK_EQUAL(outState.referenceCoordinate, diskHitInner().zCoordinate);
  BOOST_CHECK_EQUAL(outState.alpha, 0.f);
  BOOST_CHECK(outState.family == StateFamily::Forward);
  BOOST_CHECK_EQUAL(outState.flags, uint8_t{0});
  BOOST_CHECK_EQUAL(outState.absCharge, uint8_t{2});
  BOOST_CHECK(outState.pid == o2::track::PID::Pion);
}

BOOST_AUTO_TEST_CASE(BuildCellSeedDiskExactChi2ThresholdBoundary)
{
  SurfaceKinematicState replayState{};
  float replayChi2 = 0.f;
  float lastStepChi2 = 0.f;
  OperationFailureReason replayReason{};
  BOOST_REQUIRE(replayBuildCellSeedDisk(2, o2::track::PID::Pion, 1.e6f, replayState, replayChi2, lastStepChi2, replayReason));
  BOOST_REQUIRE_GT(lastStepChi2, 0.f);

  SurfaceKinematicState outState{};
  float chi2 = 0.f;
  OperationFailureReason reason{};
  DiskDiskPolicyParams accept;
  accept.trackletMinPt = DiskTrackletMinPt;
  accept.maxChi2ClusterAttachment = lastStepChi2;
  BOOST_CHECK(buildCellSeed<TransitionPolicyTag::DiskDisk>(
    diskClusterInner(), diskClusterMiddle(), diskClusterOuter(),
    diskHitInner(), diskHitMiddle(), diskHitOuter(), diskMaterial(), DiskBz,
    2, o2::track::PID::Pion, outState, chi2, accept, reason));

  const auto sentinel = diskAttachState();
  auto rejectedState = sentinel;
  float rejectedChi2 = -123.f;
  const auto before = rejectedState;
  const float chi2Before = rejectedChi2;
  DiskDiskPolicyParams reject;
  reject.trackletMinPt = DiskTrackletMinPt;
  reject.maxChi2ClusterAttachment = std::nextafter(lastStepChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::DiskDisk>(
    diskClusterInner(), diskClusterMiddle(), diskClusterOuter(),
    diskHitInner(), diskHitMiddle(), diskHitOuter(), diskMaterial(), DiskBz,
    2, o2::track::PID::Pion, rejectedState, rejectedChi2, reject, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(rejectedState, before));
  BOOST_CHECK_EQUAL(rejectedChi2, chi2Before);
}

BOOST_AUTO_TEST_CASE(BuildCellSeedDiskFailuresAreTransactional)
{
  const auto sentinel = diskAttachState();

  // Seed-construction failure: the strict inner/outer z-ordering boundary
  // (SeedGeometryDegenerate), transcribed verbatim from
  // detail::mftFwdFitCellClusters.
  {
    auto state = sentinel;
    float chi2 = -321.f;
    const auto before = state;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    DiskDiskPolicyParams params;
    params.trackletMinPt = DiskTrackletMinPt;
    params.maxChi2ClusterAttachment = 1.e6f;
    BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::DiskDisk>(
      diskClusterOuter() /* swapped: inner<=outer violates the strict boundary */, diskClusterMiddle(), diskClusterInner(),
      diskHitInner(), diskHitMiddle(), diskHitOuter(), diskMaterial(), DiskBz,
      2, o2::track::PID::Pion, state, chi2, params, reason));
    BOOST_CHECK(reason == OperationFailureReason::SeedGeometryDegenerate);
    BOOST_CHECK(bitEqual(state, before));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  }

  // Propagation-stage failure: a non-finite inner-hit target z is rejected
  // by forward::propagate's own finite-target precondition before any
  // commit.
  {
    auto badHitInner = diskHitInner();
    badHitInner.zCoordinate = std::numeric_limits<float>::infinity();
    auto state = sentinel;
    float chi2 = -321.f;
    const auto before = state;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    DiskDiskPolicyParams params;
    params.trackletMinPt = DiskTrackletMinPt;
    params.maxChi2ClusterAttachment = 1.e6f;
    BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::DiskDisk>(
      diskClusterInner(), diskClusterMiddle(), diskClusterOuter(),
      badHitInner, diskHitMiddle(), diskHitOuter(), diskMaterial(), DiskBz,
      2, o2::track::PID::Pion, state, chi2, params, reason));
    BOOST_CHECK(reason == OperationFailureReason::NonFiniteInput);
    BOOST_CHECK(bitEqual(state, before));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  }

  // Material failure: non-finite middle material.
  {
    auto badMaterial = diskMaterial();
    badMaterial[1] = NominalSurfaceMaterial{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
    auto state = sentinel;
    float chi2 = -321.f;
    const auto before = state;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    DiskDiskPolicyParams params;
    params.trackletMinPt = DiskTrackletMinPt;
    params.maxChi2ClusterAttachment = 1.e6f;
    BOOST_CHECK(!buildCellSeed<TransitionPolicyTag::DiskDisk>(
      diskClusterInner(), diskClusterMiddle(), diskClusterOuter(),
      diskHitInner(), diskHitMiddle(), diskHitOuter(), badMaterial, DiskBz,
      2, o2::track::PID::Pion, state, chi2, params, reason));
    BOOST_CHECK(reason == OperationFailureReason::MaterialFailure);
    BOOST_CHECK(bitEqual(state, before));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  }
}

BOOST_AUTO_TEST_CASE(BuildCellSeedDiskIsByteDeterministic)
{
  SurfaceKinematicState first{};
  SurfaceKinematicState second{};
  float chi2First = 0.f;
  float chi2Second = 0.f;
  OperationFailureReason reason{};
  DiskDiskPolicyParams params;
  params.trackletMinPt = DiskTrackletMinPt;
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::DiskDisk>(
    diskClusterInner(), diskClusterMiddle(), diskClusterOuter(), diskHitInner(), diskHitMiddle(), diskHitOuter(),
    diskMaterial(), DiskBz, 2, o2::track::PID::Pion, first, chi2First, params, reason));
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::DiskDisk>(
    diskClusterInner(), diskClusterMiddle(), diskClusterOuter(), diskHitInner(), diskHitMiddle(), diskHitOuter(),
    diskMaterial(), DiskBz, 2, o2::track::PID::Pion, second, chi2Second, params, reason));
  BOOST_CHECK(bitEqual(first, second));
  BOOST_CHECK_EQUAL(chi2First, chi2Second);
}

BOOST_AUTO_TEST_CASE(BuildCellSeedDiskActivatesEnergyLossUnlikeLegacyMcsOnlyPath)
{
  // Intentional physics difference (kickoff message): unlike the legacy,
  // MCS-only buildCellSeed<DiskDisk>, the native path reads
  // material.arealDensityGPerCm2 and therefore changes the fitted q/pT
  // (parameters[4]) when energy loss is active, not only the covariance.
  auto noEnergyLoss = diskMaterial();
  for (auto& m : noEnergyLoss) {
    m.arealDensityGPerCm2 = 0.f;
  }
  auto withEnergyLoss = diskMaterial(); // arealDensityGPerCm2 > 0 for all three slots

  SurfaceKinematicState stateNoLoss{};
  SurfaceKinematicState stateWithLoss{};
  float chi2NoLoss = 0.f;
  float chi2WithLoss = 0.f;
  OperationFailureReason reason{};
  DiskDiskPolicyParams params;
  params.trackletMinPt = DiskTrackletMinPt;
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::DiskDisk>(
    diskClusterInner(), diskClusterMiddle(), diskClusterOuter(), diskHitInner(), diskHitMiddle(), diskHitOuter(),
    noEnergyLoss, DiskBz, 2, o2::track::PID::Pion, stateNoLoss, chi2NoLoss, params, reason));
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::DiskDisk>(
    diskClusterInner(), diskClusterMiddle(), diskClusterOuter(), diskHitInner(), diskHitMiddle(), diskHitOuter(),
    withEnergyLoss, DiskBz, 2, o2::track::PID::Pion, stateWithLoss, chi2WithLoss, params, reason));

  BOOST_CHECK_NE(stateNoLoss.parameters[4], stateWithLoss.parameters[4]);
}

// ===========================================================================
// cellsAreCompatible<CylinderCylinder> / <DiskDisk>
// ===========================================================================

namespace
{
SurfaceKinematicState compatBarrelCurrent()
{
  SurfaceKinematicState state{};
  state.parameters[0] = 0.1f;
  state.parameters[1] = 1.f;
  state.parameters[2] = 0.05f;
  state.parameters[3] = 0.4f;
  state.parameters[4] = 0.2f;
  state.referenceCoordinate = 5.f;
  state.alpha = 0.3f;
  state.family = StateFamily::Barrel;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  state.covariance[packedCovarianceIndex(0, 0)] = 1.e-4f;
  state.covariance[packedCovarianceIndex(1, 1)] = 1.e-4f;
  state.covariance[packedCovarianceIndex(2, 2)] = 1.e-4f;
  state.covariance[packedCovarianceIndex(3, 3)] = 1.e-4f;
  state.covariance[packedCovarianceIndex(4, 4)] = 1.e-2f;
  return state;
}

SurfaceKinematicState compatBarrelNext()
{
  SurfaceKinematicState state{};
  state.parameters[0] = 0.2f;
  state.parameters[1] = 1.1f;
  state.parameters[2] = 0.06f;
  state.parameters[3] = 0.42f;
  state.parameters[4] = 0.21f;
  state.referenceCoordinate = 4.f;
  state.alpha = 0.3f;
  state.family = StateFamily::Barrel;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  state.covariance[packedCovarianceIndex(0, 0)] = 1.e-4f;
  state.covariance[packedCovarianceIndex(1, 1)] = 1.e-4f;
  state.covariance[packedCovarianceIndex(2, 2)] = 1.e-4f;
  state.covariance[packedCovarianceIndex(3, 3)] = 1.e-4f;
  state.covariance[packedCovarianceIndex(4, 4)] = 1.e-2f;
  return state;
}

constexpr float CompatBz = 0.5f;

SurfaceKinematicState compatDiskCurrent()
{
  SurfaceKinematicState state{};
  state.parameters[0] = 0.f;
  state.parameters[1] = 0.f;
  state.parameters[2] = 0.f;
  state.parameters[3] = -1.2f;
  state.parameters[4] = 0.1f;
  state.referenceCoordinate = 0.f;
  state.family = StateFamily::Forward;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  state.covariance[packedCovarianceIndex(0, 0)] = 1.e-2f;
  state.covariance[packedCovarianceIndex(1, 1)] = 1.e-2f;
  state.covariance[packedCovarianceIndex(2, 2)] = 1.e-3f;
  state.covariance[packedCovarianceIndex(3, 3)] = 1.e-3f;
  state.covariance[packedCovarianceIndex(4, 4)] = 1.e-2f;
  return state;
}

SurfaceKinematicState compatDiskNext()
{
  SurfaceKinematicState state{};
  state.parameters[0] = 0.05f;
  state.parameters[1] = -0.03f;
  state.parameters[2] = 0.02f;
  state.parameters[3] = -1.19f;
  state.parameters[4] = 0.11f;
  state.referenceCoordinate = -1.f;
  state.family = StateFamily::Forward;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  state.covariance[packedCovarianceIndex(0, 0)] = 1.e-2f;
  state.covariance[packedCovarianceIndex(1, 1)] = 1.e-2f;
  state.covariance[packedCovarianceIndex(2, 2)] = 1.e-3f;
  state.covariance[packedCovarianceIndex(3, 3)] = 1.e-3f;
  state.covariance[packedCovarianceIndex(4, 4)] = 1.e-2f;
  return state;
}

constexpr float CompatDiskBz = 0.5f;

} // namespace

BOOST_AUTO_TEST_CASE(CellsAreCompatibleBarrelAcceptsAtExactThresholdAndRejectsBelow)
{
  const auto current = compatBarrelCurrent();
  const auto next = compatBarrelNext();

  auto reference = next;
  OperationFailureReason reason{};
  BOOST_REQUIRE(barrel::rotate(reference, current.alpha, reason));
  BOOST_REQUIRE(barrel::propagate(reference, current.referenceCoordinate, CompatBz, reason));
  float refChi2 = 0.f;
  BOOST_REQUIRE(barrel::stateChi2(current, reference, refChi2, reason));
  BOOST_REQUIRE_GT(refChi2, 0.f);

  CylinderCylinderPolicyParams accept;
  accept.maxChi2ClusterAttachment = refChi2;
  BOOST_CHECK(cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(current, next, -1, -1, CompatBz, accept));

  CylinderCylinderPolicyParams reject;
  reject.maxChi2ClusterAttachment = std::nextafter(refChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(current, next, -1, -1, CompatBz, reject));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleBarrelRotationAndPropagationFailuresAreRejectedNotThrown)
{
  CylinderCylinderPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;

  // Rotation failure: far frame angle.
  {
    auto current = compatBarrelCurrent();
    current.alpha = 0.f;
    auto next = compatBarrelNext();
    next.alpha = 3.f;
    BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(current, next, -1, -1, CompatBz, permissive));
  }

  // Propagation failure: extreme q2pt on `next` makes the requested step
  // unreachable.
  {
    auto current = compatBarrelCurrent();
    auto next = compatBarrelNext();
    next.parameters[4] = 2000.f;
    BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(current, next, -1, -1, CompatBz, permissive));
  }
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleBarrelFamilyMismatchFailsClosed)
{
  auto current = compatBarrelCurrent();
  current.family = StateFamily::Forward;
  const auto next = compatBarrelNext();
  CylinderCylinderPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(current, next, -1, -1, CompatBz, permissive));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleBarrelIgnoresClusterIndicesAndNeverMutatesInputs)
{
  const auto current = compatBarrelCurrent();
  const auto next = compatBarrelNext();
  const auto currentBefore = current;
  const auto nextBefore = next;
  CylinderCylinderPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;

  const bool r1 = cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(current, next, 10, 10, CompatBz, permissive);
  const bool r2 = cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(current, next, 10, 999, CompatBz, permissive);
  BOOST_CHECK_EQUAL(r1, r2);
  BOOST_CHECK(bitEqual(current, currentBefore));
  BOOST_CHECK(bitEqual(next, nextBefore));

  const bool r3 = cellsAreCompatible<TransitionPolicyTag::CylinderCylinder>(current, next, 10, 10, CompatBz, permissive);
  BOOST_CHECK_EQUAL(r1, r3);
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleDiskAcceptsAtExactThresholdAndRejectsBelow)
{
  const auto current = compatDiskCurrent();
  const auto next = compatDiskNext();

  auto reference = next;
  OperationFailureReason reason{};
  BOOST_REQUIRE(forward::propagate<forward::PropagationModel::Linear>(reference, current.referenceCoordinate, 0.f, reason));
  float refChi2 = 0.f;
  BOOST_REQUIRE(forward::stateChi2(current, reference, refChi2, reason));
  BOOST_REQUIRE_GT(refChi2, 0.f);

  DiskDiskPolicyParams accept;
  accept.maxChi2ClusterAttachment = refChi2;
  BOOST_CHECK(cellsAreCompatible<TransitionPolicyTag::DiskDisk>(current, next, 42, 42, 0.f, accept));

  DiskDiskPolicyParams reject;
  reject.maxChi2ClusterAttachment = std::nextafter(refChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::DiskDisk>(current, next, 42, 42, 0.f, reject));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleDiskClusterIndexContinuityRejectedBeforeAnyPropagation)
{
  const auto current = compatDiskCurrent();
  const auto next = compatDiskNext();
  DiskDiskPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  // Otherwise-compatible states, but the temporary raw cluster-index
  // continuity input (see the header doc) rejects first.
  BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::DiskDisk>(current, next, 5, 6, 0.f, permissive));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleDiskPropagationFailureIsRejectedNotThrown)
{
  auto current = compatDiskCurrent();
  current.referenceCoordinate = -5.f; // dz != 0 relative to `next`
  auto next = compatDiskNext();
  next.parameters[3] = 0.f; // tanl == 0: propagateLinear's UnreachableTarget
  DiskDiskPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::DiskDisk>(current, next, 7, 7, 0.f, permissive));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleDiskFamilyMismatchFailsClosed)
{
  auto current = compatDiskCurrent();
  current.family = StateFamily::Barrel;
  const auto next = compatDiskNext();
  DiskDiskPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_CHECK(!cellsAreCompatible<TransitionPolicyTag::DiskDisk>(current, next, 3, 3, 0.f, permissive));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleDiskInputsUntouchedAndByteDeterministic)
{
  const auto current = compatDiskCurrent();
  const auto next = compatDiskNext();
  const auto currentBefore = current;
  const auto nextBefore = next;
  DiskDiskPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  const bool r1 = cellsAreCompatible<TransitionPolicyTag::DiskDisk>(current, next, 1, 1, 0.f, permissive);
  const bool r2 = cellsAreCompatible<TransitionPolicyTag::DiskDisk>(current, next, 1, 1, 0.f, permissive);
  BOOST_CHECK_EQUAL(r1, r2);
  BOOST_CHECK(bitEqual(current, currentBefore));
  BOOST_CHECK(bitEqual(next, nextBefore));
}

// ===========================================================================
// attachHit<CylinderCylinder>
// ===========================================================================

BOOST_AUTO_TEST_CASE(AttachHitBarrelSuccessAndExactChi2Threshold)
{
  const auto state0 = barrelAttachState();
  const auto hit = barrelAttachHit();
  const auto material = barrelAttachMaterial();

  auto probe = state0;
  float probeChi2 = 0.f;
  OperationFailureReason reason{};
  CylinderCylinderPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::CylinderCylinder>(probe, hit, material, BarrelAttachBz, probeChi2, permissive, reason));
  BOOST_REQUIRE_GT(probeChi2, 0.f);

  auto accepted = state0;
  float acceptedChi2 = 0.f;
  CylinderCylinderPolicyParams accept;
  accept.maxChi2ClusterAttachment = probeChi2;
  BOOST_CHECK(attachHit<TransitionPolicyTag::CylinderCylinder>(accepted, hit, material, BarrelAttachBz, acceptedChi2, accept, reason));
  BOOST_CHECK(bitEqual(accepted, probe));
  BOOST_CHECK_EQUAL(acceptedChi2, probeChi2);

  auto rejected = state0;
  float rejectedChi2 = -1.f;
  const auto before = rejected;
  const float chi2Before = rejectedChi2;
  CylinderCylinderPolicyParams reject;
  reject.maxChi2ClusterAttachment = std::nextafter(probeChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!attachHit<TransitionPolicyTag::CylinderCylinder>(rejected, hit, material, BarrelAttachBz, rejectedChi2, reject, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(rejected, before));
  BOOST_CHECK_EQUAL(rejectedChi2, chi2Before);
}

BOOST_AUTO_TEST_CASE(AttachHitBarrelEachFailureStagePreservesStateTransactionally)
{
  const auto state0 = barrelAttachState();
  CylinderCylinderPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;

  auto checkFailure = [&](const o2::its::TrackingFrameInfo& hit, const NominalSurfaceMaterial& material,
                          OperationFailureReason expectedOneOf1, OperationFailureReason expectedOneOf2) {
    auto state = state0;
    float chi2 = -1.f;
    const auto before = state;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    BOOST_CHECK(!attachHit<TransitionPolicyTag::CylinderCylinder>(state, hit, material, BarrelAttachBz, chi2, permissive, reason));
    BOOST_CHECK(reason == expectedOneOf1 || reason == expectedOneOf2);
    BOOST_CHECK(bitEqual(state, before));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  };

  // Rotation failure.
  {
    auto farHit = barrelAttachHit();
    farHit.alphaTrackingFrame = state0.alpha + 3.f;
    checkFailure(farHit, barrelAttachMaterial(), OperationFailureReason::RotationFailure, OperationFailureReason::RotationFailure);
  }

  // Propagation failure.
  {
    auto farHit = barrelAttachHit();
    farHit.xTrackingFrame = -50000.f;
    checkFailure(farHit, barrelAttachMaterial(), OperationFailureReason::UnreachableTarget, OperationFailureReason::PropagationFailure);
  }

  // Material failure.
  {
    const NominalSurfaceMaterial badMaterial{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
    checkFailure(barrelAttachHit(), badMaterial, OperationFailureReason::MaterialFailure, OperationFailureReason::MaterialFailure);
  }
}

BOOST_AUTO_TEST_CASE(AttachHitBarrelNegativeChi2IsRejectedMatchingLegacyInclusiveCut)
{
  // No-op rotate (hit shares state's alpha) and no-op propagate (hit's
  // target x equals state's own referenceCoordinate) plus a no-op material
  // budget ({0,0} is the documented unconditional no-op) keep the state
  // byte-identical to barrelAttachState() through predictedChi2, so its
  // known covariance can be reasoned about directly: a pathologically large
  // measurement cross-covariance (uv) makes the combined 2x2 determinant
  // negative, which residualInverse's own gate does not reject outright
  // (only exact-zero/non-finite determinants are), producing a negative
  // predicted chi2 -- the same `< 0.f` established rejection
  // attachHit<CylinderCylinder> already applies today. barrel::update shares
  // the identical residualInverse gate and therefore cannot independently
  // fail once predictedChi2 has already succeeded with the same inputs; the
  // two checks share one deterministic failure precedence (predictedChi2,
  // evaluated before update, is always the one that observes a bad
  // residualInverse first).
  auto state = barrelAttachState();
  auto hit = barrelAttachHit();
  hit.alphaTrackingFrame = state.alpha;
  hit.xTrackingFrame = state.referenceCoordinate;
  hit.covarianceTrackingFrame[1] = 50.f; // huge uv cross term
  const NominalSurfaceMaterial noopMaterial{0.f, 0.f};
  const auto before = state;
  float chi2 = -1.f;
  const float chi2Before = chi2;
  OperationFailureReason reason{};
  CylinderCylinderPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_CHECK(!attachHit<TransitionPolicyTag::CylinderCylinder>(state, hit, noopMaterial, BarrelAttachBz, chi2, permissive, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(state, before));
  BOOST_CHECK_EQUAL(chi2, chi2Before);
}

BOOST_AUTO_TEST_CASE(AttachHitBarrelIsChargeAwareUnlikeNeutralMaterialCorrection)
{
  // PID/absCharge-aware behavior: a neutral state (absCharge == 0) takes the
  // material kernel's documented unconditional no-op path, while a charged
  // state with identical kinematics picks up a Highland covariance
  // contribution -- the results must differ.
  auto neutral = barrelAttachState();
  neutral.absCharge = 0;
  auto charged = barrelAttachState();
  charged.absCharge = 1;
  const auto hit = barrelAttachHit();
  const auto material = barrelAttachMaterial();

  float neutralChi2 = 0.f;
  float chargedChi2 = 0.f;
  OperationFailureReason reason{};
  CylinderCylinderPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::CylinderCylinder>(neutral, hit, material, BarrelAttachBz, neutralChi2, permissive, reason));
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::CylinderCylinder>(charged, hit, material, BarrelAttachBz, chargedChi2, permissive, reason));
  BOOST_CHECK(!bitEqual(neutral, charged));
}

BOOST_AUTO_TEST_CASE(AttachHitBarrelIsByteDeterministic)
{
  auto first = barrelAttachState();
  auto second = barrelAttachState();
  float chi2First = 0.f;
  float chi2Second = 0.f;
  OperationFailureReason reason{};
  CylinderCylinderPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::CylinderCylinder>(first, barrelAttachHit(), barrelAttachMaterial(), BarrelAttachBz, chi2First, permissive, reason));
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::CylinderCylinder>(second, barrelAttachHit(), barrelAttachMaterial(), BarrelAttachBz, chi2Second, permissive, reason));
  BOOST_CHECK(bitEqual(first, second));
  BOOST_CHECK_EQUAL(chi2First, chi2Second);
}

// ===========================================================================
// attachHit<DiskDisk>
// ===========================================================================

BOOST_AUTO_TEST_CASE(AttachHitDiskSuccessAndExactChi2Threshold)
{
  const auto state0 = diskAttachState();
  const auto hit = diskAttachHit();
  const auto material = diskAttachMaterial();

  auto probe = state0;
  float probeChi2 = 0.f;
  OperationFailureReason reason{};
  DiskDiskPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(probe, hit, material, DiskAttachBz, probeChi2, permissive, reason));
  BOOST_REQUIRE_GT(probeChi2, 0.f);

  auto accepted = state0;
  float acceptedChi2 = 0.f;
  DiskDiskPolicyParams accept;
  accept.maxChi2ClusterAttachment = probeChi2;
  BOOST_CHECK(attachHit<TransitionPolicyTag::DiskDisk>(accepted, hit, material, DiskAttachBz, acceptedChi2, accept, reason));
  BOOST_CHECK(bitEqual(accepted, probe));
  BOOST_CHECK_EQUAL(acceptedChi2, probeChi2);

  auto rejected = state0;
  float rejectedChi2 = -1.f;
  const auto before = rejected;
  const float chi2Before = rejectedChi2;
  DiskDiskPolicyParams reject;
  reject.maxChi2ClusterAttachment = std::nextafter(probeChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!attachHit<TransitionPolicyTag::DiskDisk>(rejected, hit, material, DiskAttachBz, rejectedChi2, reject, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(rejected, before));
  BOOST_CHECK_EQUAL(rejectedChi2, chi2Before);
}

BOOST_AUTO_TEST_CASE(AttachHitDiskEachFailureStagePreservesStateTransactionally)
{
  const auto state0 = diskAttachState();
  DiskDiskPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;

  // Propagation failure: tanl == 0 at zero field rejects with
  // UnreachableTarget in forward::propagate<Linear>.
  {
    auto zeroTanl = state0;
    zeroTanl.parameters[3] = 0.f;
    auto hit = diskAttachHit();
    hit.zCoordinate = -60.f; // dz != 0
    auto state = zeroTanl;
    float chi2 = -1.f;
    const auto before = state;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    BOOST_CHECK(!attachHit<TransitionPolicyTag::DiskDisk>(state, hit, diskAttachMaterial(), 0.f, chi2, permissive, reason));
    BOOST_CHECK(reason == OperationFailureReason::UnreachableTarget);
    BOOST_CHECK(bitEqual(state, before));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  }

  // Material failure.
  {
    const NominalSurfaceMaterial badMaterial{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
    auto state = state0;
    float chi2 = -1.f;
    const auto before = state;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    BOOST_CHECK(!attachHit<TransitionPolicyTag::DiskDisk>(state, diskAttachHit(), badMaterial, DiskAttachBz, chi2, permissive, reason));
    BOOST_CHECK(reason == OperationFailureReason::MaterialFailure);
    BOOST_CHECK(bitEqual(state, before));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  }
}

BOOST_AUTO_TEST_CASE(AttachHitDiskIsChargeAwareUnlikeNeutralMaterialCorrection)
{
  auto neutral = diskAttachState();
  neutral.absCharge = 0;
  auto charged = diskAttachState();
  charged.absCharge = 1;
  const auto hit = diskAttachHit();
  const auto material = diskAttachMaterial();

  float neutralChi2 = 0.f;
  float chargedChi2 = 0.f;
  OperationFailureReason reason{};
  DiskDiskPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(neutral, hit, material, DiskAttachBz, neutralChi2, permissive, reason));
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(charged, hit, material, DiskAttachBz, chargedChi2, permissive, reason));
  BOOST_CHECK(!bitEqual(neutral, charged));
}

BOOST_AUTO_TEST_CASE(AttachHitDiskActivatesEnergyLossUnlikeLegacyMcsOnlyPath)
{
  auto noLossMaterial = diskAttachMaterial();
  noLossMaterial.arealDensityGPerCm2 = 0.f;
  auto withLossMaterial = diskAttachMaterial();

  auto stateNoLoss = diskAttachState();
  auto stateWithLoss = diskAttachState();
  const auto hit = diskAttachHit();
  float chi2NoLoss = 0.f;
  float chi2WithLoss = 0.f;
  OperationFailureReason reason{};
  DiskDiskPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(stateNoLoss, hit, noLossMaterial, DiskAttachBz, chi2NoLoss, permissive, reason));
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(stateWithLoss, hit, withLossMaterial, DiskAttachBz, chi2WithLoss, permissive, reason));
  BOOST_CHECK_NE(stateNoLoss.parameters[4], stateWithLoss.parameters[4]);
}

BOOST_AUTO_TEST_CASE(AttachHitDiskIsByteDeterministic)
{
  auto first = diskAttachState();
  auto second = diskAttachState();
  float chi2First = 0.f;
  float chi2Second = 0.f;
  OperationFailureReason reason{};
  DiskDiskPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(first, diskAttachHit(), diskAttachMaterial(), DiskAttachBz, chi2First, permissive, reason));
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(second, diskAttachHit(), diskAttachMaterial(), DiskAttachBz, chi2Second, permissive, reason));
  BOOST_CHECK(bitEqual(first, second));
  BOOST_CHECK_EQUAL(chi2First, chi2Second);
}

// ===========================================================================
// Compatibility-projection coverage (private TrackingFrameInfo ->
// SurfaceMeasurement boundary, exercised indirectly through attachHit).
// ===========================================================================

BOOST_AUTO_TEST_CASE(BarrelProjectionUsesFullCovarianceIncludingCrossTerm)
{
  const auto state0 = barrelAttachState();
  auto lowCrossTerm = barrelAttachHit();
  lowCrossTerm.covarianceTrackingFrame[1] = 0.f;
  auto highCrossTerm = barrelAttachHit();
  highCrossTerm.covarianceTrackingFrame[1] = 0.03f;

  auto stateLow = state0;
  auto stateHigh = state0;
  float chi2Low = 0.f;
  float chi2High = 0.f;
  OperationFailureReason reason{};
  CylinderCylinderPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::CylinderCylinder>(stateLow, lowCrossTerm, barrelAttachMaterial(), BarrelAttachBz, chi2Low, permissive, reason));
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::CylinderCylinder>(stateHigh, highCrossTerm, barrelAttachMaterial(), BarrelAttachBz, chi2High, permissive, reason));
  BOOST_CHECK_NE(chi2Low, chi2High);
}

BOOST_AUTO_TEST_CASE(ForwardProjectionIsDiagonalOnlyAndIgnoresUnreadTrackingFrameFields)
{
  const auto state0 = diskAttachState();
  const auto baseline = diskAttachHit();

  auto varyingCrossTerm = baseline;
  varyingCrossTerm.covarianceTrackingFrame[1] = 999.f; // forward never reads this slot

  auto varyingUnreadFields = baseline;
  varyingUnreadFields.xTrackingFrame = 12345.f;
  varyingUnreadFields.alphaTrackingFrame = 6.7f;
  varyingUnreadFields.positionTrackingFrame = {-999.f, 999.f};

  DiskDiskPolicyParams permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  OperationFailureReason reason{};

  auto stateBaseline = state0;
  float chi2Baseline = 0.f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(stateBaseline, baseline, diskAttachMaterial(), DiskAttachBz, chi2Baseline, permissive, reason));

  auto stateCrossTerm = state0;
  float chi2CrossTerm = 0.f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(stateCrossTerm, varyingCrossTerm, diskAttachMaterial(), DiskAttachBz, chi2CrossTerm, permissive, reason));
  BOOST_CHECK(bitEqual(stateBaseline, stateCrossTerm));
  BOOST_CHECK_EQUAL(chi2Baseline, chi2CrossTerm);

  auto stateUnreadFields = state0;
  float chi2UnreadFields = 0.f;
  BOOST_REQUIRE(attachHit<TransitionPolicyTag::DiskDisk>(stateUnreadFields, varyingUnreadFields, diskAttachMaterial(), DiskAttachBz, chi2UnreadFields, permissive, reason));
  BOOST_CHECK(bitEqual(stateBaseline, stateUnreadFields));
  BOOST_CHECK_EQUAL(chi2Baseline, chi2UnreadFields);
}
