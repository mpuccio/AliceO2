// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT CellFindingNative
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
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/detail/CandidateFinding.h"
#include "ITStracking/Cluster.h"

/// Focused coverage for the explicit cylinder/disk SurfaceKinematicState
/// cell-seed, compatibility, and hit-attachment leaves
/// (CandidateFinding.h/.cxx). The leaves are production-wired through
/// TrackerTraits; this test exercises their numerical contracts directly,
/// while the separate orchestration tests cover their TrackerTraits callers.
///
/// Oracle strategy (see the Stage-B design report's "Precision/oracles"
/// split, referenced from CandidateFinding.h):
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
// Deliberately invalid: buildCylinderCellSeed must never read the
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

// Test-local field-mapping helper (not a production API), combining a
// full-precision hit's tracking-frame fields (reference coordinate/frame
// angle/measured local coordinates/covariance -- everything predictedChi2/
// update/rotate/propagate read) into one SurfaceMeasurement. Also used
// directly (without an accompanying cluster) as the single measurement input
// for attachCylinderHit below.
SurfaceMeasurement barrelMeasurementFromHit(const o2::its::TrackingFrameInfo& hit)
{
  SurfaceMeasurement measurement{};
  measurement.frame.q = hit.xTrackingFrame;
  measurement.frame.frameAngle = hit.alphaTrackingFrame;
  measurement.frame.u = hit.positionTrackingFrame[0];
  measurement.frame.v = hit.positionTrackingFrame[1];
  measurement.covariance.uu = hit.covarianceTrackingFrame[0];
  measurement.covariance.uv = hit.covarianceTrackingFrame[1];
  measurement.covariance.vv = hit.covarianceTrackingFrame[2];
  return measurement;
}

// Combines a cluster's global position with a hit's tracking-frame fields:
// the single SurfaceMeasurement now standing in for the retired
// {Cluster, TrackingFrameInfo} pair at each of cell-seed leaves's three
// candidate positions.
SurfaceMeasurement barrelMeasurementFor(const o2::its::Cluster& cluster, const o2::its::TrackingFrameInfo& hit)
{
  auto measurement = barrelMeasurementFromHit(hit);
  measurement.global = {cluster.xCoordinate, cluster.yCoordinate, cluster.zCoordinate};
  return measurement;
}

SurfaceMeasurement barrelMeasurementInner() { return barrelMeasurementFor(barrelClusterInner(), barrelHitInner()); }
SurfaceMeasurement barrelMeasurementMiddle() { return barrelMeasurementFor(barrelClusterMiddle(), barrelHitMiddle()); }
SurfaceMeasurement barrelMeasurementOuter() { return barrelMeasurementFor(barrelClusterOuter(), barrelHitOuter()); }

// Independent re-transcription of native buildCylinderCellSeed's
// own documented call sequence, built directly on the public barrel::
// primitives rather than calling the production function under test.
// `nSteps` limits how many of the middle/inner attach steps run (1 = middle
// only, used to read out the inner step's own marginal predicted chi2 for
// boundary tests), mirroring the nSteps convention already established by
// the established nSteps convention in this test.
bool replayBuildCellSeedBarrel(uint8_t absCharge, o2::track::PID pid, float maxChi2,
                               SurfaceKinematicState& outState, float& outChi2, float& lastStepChi2,
                               OperationFailureReason& reason, int nSteps = 2)
{
  SurfaceKinematicState state{};
  if (!barrel::buildSeed(barrelMeasurementInner(), barrelMeasurementMiddle(), barrelMeasurementOuter(), BarrelBz,
                         absCharge, pid, state, reason)) {
    return false;
  }
  const std::array<SurfaceMeasurement, 2> steps{barrelMeasurementMiddle(), barrelMeasurementInner()};
  const std::array<NominalSurfaceMaterial, 2> stepMaterial{barrelMaterialMiddle(), barrelMaterialInner()};
  float chi2 = 0.f;
  float stepChi2 = 0.f;
  for (int step = 0; step < nSteps; ++step) {
    const bool isLast = (nSteps == 2) && (step == 1);
    const auto& measurement = steps[step];
    if (!barrel::rotate(state, measurement.frame.frameAngle, reason)) {
      return false;
    }
    if (!barrel::propagate(state, measurement.frame.q, BarrelBz, reason)) {
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

// Combines a cluster's global position (including z, read by
// forward::buildSeed's z-ordering/direction formula) with a hit's diagonal
// covariance and a reference z (Disk field mapping: reference z ->
// measurement.frame.q, contractually == global.z for the accepted disk
// adapter): the single SurfaceMeasurement now standing in for the retired
// {Cluster, TrackingFrameInfo} pair at each of cell-seed leaves's three
// candidate positions.
SurfaceMeasurement diskMeasurementFor(const o2::its::Cluster& cluster, const o2::its::TrackingFrameInfo& hit)
{
  auto measurement = diskMeasurementFromHit(hit);
  measurement.global = {cluster.xCoordinate, cluster.yCoordinate, cluster.zCoordinate};
  measurement.frame.q = cluster.zCoordinate;
  return measurement;
}

SurfaceMeasurement diskMeasurementInner() { return diskMeasurementFor(diskClusterInner(), diskHitInner()); }
SurfaceMeasurement diskMeasurementMiddle() { return diskMeasurementFor(diskClusterMiddle(), diskHitMiddle()); }
SurfaceMeasurement diskMeasurementOuter() { return diskMeasurementFor(diskClusterOuter(), diskHitOuter()); }

// Replays the accepted propagation boundary while this test independently
// checks the surrounding cell-seed operation order.
bool replayForwardPropagateAcceptedModel(SurfaceKinematicState& state, float targetZ, float bz,
                                         OperationFailureReason& reason)
{
  return Propagator::propagateForward(state, targetZ, bz, reason);
}

// Independent re-transcription of native buildDiskCellSeed's own
// documented call sequence (outer, middle, inner; material[2]/[1]/[0]).
bool replayBuildCellSeedDisk(uint8_t absCharge, o2::track::PID pid, float maxChi2,
                             SurfaceKinematicState& outState, float& outChi2, float& lastStepChi2,
                             OperationFailureReason& reason, int nSteps = 3)
{
  SurfaceKinematicState state{};
  if (!forward::buildSeed(diskMeasurementInner(), diskMeasurementMiddle(), diskMeasurementOuter(),
                          DiskBz, DiskTrackletMinPt, absCharge, pid, state, reason)) {
    return false;
  }
  const std::array<SurfaceMeasurement, 3> steps{diskMeasurementOuter(), diskMeasurementMiddle(), diskMeasurementInner()};
  const std::array<NominalSurfaceMaterial, 3> stepMaterial{diskMaterialOuter(), diskMaterialMiddle(), diskMaterialInner()};
  float chi2 = 0.f;
  float stepChi2 = 0.f;
  for (int step = 0; step < nSteps; ++step) {
    const bool isLast = (nSteps == 3) && (step == 2);
    const auto& measurement = steps[step];
    if (!replayForwardPropagateAcceptedModel(state, measurement.frame.q, DiskBz, reason)) {
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

// Test-local field-mapping helper: builds the SurfaceMeasurement
// attachDiskHit reads from a single legacy hit (Disk field mapping:
// global coordinates -> measurement.global, reference z -> measurement.
// frame.q [read as the propagate target, in place of the retired
// hit.zCoordinate], measured covariance -> measurement.covariance).
SurfaceMeasurement diskAttachMeasurementFrom(const o2::its::TrackingFrameInfo& hit)
{
  auto measurement = diskMeasurementFromHit(hit);
  measurement.frame.q = hit.zCoordinate;
  return measurement;
}

SurfaceMeasurement diskAttachMeasurement() { return diskAttachMeasurementFrom(diskAttachHit()); }

SurfaceDescriptor directionSurface(SurfaceKind kind, float referenceCoordinate)
{
  SurfaceDescriptor surface{};
  surface.kind = kind;
  surface.referenceCoordinate = referenceCoordinate;
  return surface;
}

SurfaceMeasurement diskDirectionMeasurement(float x, float y,
                                            float varianceX, float covarianceXY, float varianceY)
{
  SurfaceMeasurement measurement{};
  measurement.global = {x, y, 1234.f}; // disk observation takes z from the surface
  measurement.covariance = {varianceX, covarianceXY, varianceY};
  return measurement;
}

SurfaceMeasurement cylinderDirectionMeasurement(float q, float u, float z,
                                                float varianceU, float covarianceUZ, float varianceZ)
{
  SurfaceMeasurement measurement{};
  measurement.global = {999.f, 999.f, -1234.f}; // cylinder observation uses its tracking frame
  measurement.frame.q = q;
  measurement.frame.u = u;
  measurement.frame.v = z;
  measurement.covariance = {varianceU, covarianceUZ, varianceZ};
  return measurement;
}

SurfaceMeasurement rotateDiskMeasurement(const SurfaceMeasurement& measurement, double angle)
{
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  const double x = measurement.global.x;
  const double y = measurement.global.y;
  const double varianceX = measurement.covariance.uu;
  const double covarianceXY = measurement.covariance.uv;
  const double varianceY = measurement.covariance.vv;
  auto rotated = measurement;
  rotated.global.x = static_cast<float>(cosine * x - sine * y);
  rotated.global.y = static_cast<float>(sine * x + cosine * y);
  rotated.covariance.uu = static_cast<float>(cosine * cosine * varianceX -
                                             2. * cosine * sine * covarianceXY +
                                             sine * sine * varianceY);
  rotated.covariance.uv = static_cast<float>(cosine * sine * (varianceX - varianceY) +
                                             (cosine * cosine - sine * sine) * covarianceXY);
  rotated.covariance.vv = static_cast<float>(sine * sine * varianceX +
                                             2. * cosine * sine * covarianceXY +
                                             cosine * cosine * varianceY);
  return rotated;
}

} // namespace

// ===========================================================================
// cell direction compatibility
// ===========================================================================

BOOST_AUTO_TEST_CASE(CellDirectionCompatibilityMatchesFullAnalyticOracle)
{
  const std::array<DirectionObservation, 3> observations{{
    {2., 1., 0.04, 0.006, 0.09},
    {4., 4., 0.05, -0.004, 0.08},
    {7., 8., 0.03, 0.003, 0.07},
  }};
  CellDirectionCompatibility compatibility{};
  BOOST_REQUIRE(cellDirectionsAreCompatible(observations, 1.f, compatibility));
  BOOST_CHECK_CLOSE_FRACTION(compatibility.residual, 1., 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(compatibility.variance, 6.55, 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(compatibility.chi2, 1. / 6.55, 1.e-12);
}

BOOST_AUTO_TEST_CASE(DiskDirectionObservationProjectsFullGlobalXYCovariance)
{
  const auto surface = directionSurface(SurfaceKind::Disk, -50.f);
  const auto measurement = diskDirectionMeasurement(3.f, 4.f, 0.04f, 0.006f, 0.09f);
  DirectionObservation observation{};
  BOOST_REQUIRE(makeDiskDirectionObservation(surface, measurement, observation));
  BOOST_CHECK_CLOSE_FRACTION(observation.r, 5., 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(observation.z, -50., 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(observation.varianceR, 0.07776, 1.e-7);
  BOOST_CHECK_SMALL(observation.covarianceRZ, 1.e-15);
  BOOST_CHECK_SMALL(observation.varianceZ, 1.e-15);
}

BOOST_AUTO_TEST_CASE(CylinderDirectionObservationProjectsTrackingFrameCovariance)
{
  const auto surface = directionSurface(SurfaceKind::Cylinder, 3.5f);
  const auto measurement = cylinderDirectionMeasurement(3.f, 4.f, -2.2f, 0.04f, 0.006f, 0.09f);
  DirectionObservation observation{};
  BOOST_REQUIRE(makeCylinderDirectionObservation(surface, measurement, observation));
  BOOST_CHECK_CLOSE_FRACTION(observation.r, 5., 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(observation.z, -2.2, 1.e-6);
  BOOST_CHECK_CLOSE_FRACTION(observation.varianceR, 0.0256, 1.e-7);
  BOOST_CHECK_CLOSE_FRACTION(observation.covarianceRZ, 0.0048, 1.e-7);
  BOOST_CHECK_CLOSE_FRACTION(observation.varianceZ, 0.09, 1.e-6);
}

BOOST_AUTO_TEST_CASE(DiskDirectionCompatibilityIsInvariantUnderGlobalZRotation)
{
  const std::array<SurfaceDescriptor, 3> surfaces{{
    directionSurface(SurfaceKind::Disk, 1.f),
    directionSurface(SurfaceKind::Disk, 4.f),
    directionSurface(SurfaceKind::Disk, 8.f),
  }};
  const std::array<SurfaceMeasurement, 3> measurements{{
    diskDirectionMeasurement(1.2f, 1.6f, 0.04f, 0.006f, 0.09f),
    diskDirectionMeasurement(-2.4f, 3.2f, 0.05f, -0.004f, 0.08f),
    diskDirectionMeasurement(4.2f, -5.6f, 0.03f, 0.003f, 0.07f),
  }};
  std::array<DirectionObservation, 3> baseline{};
  std::array<DirectionObservation, 3> rotated{};
  for (std::size_t i = 0; i < measurements.size(); ++i) {
    BOOST_REQUIRE(makeDirectionObservation(surfaces[i], measurements[i], baseline[i]));
    const auto rotatedMeasurement = rotateDiskMeasurement(measurements[i], 0.73);
    BOOST_REQUIRE(makeDirectionObservation(surfaces[i], rotatedMeasurement, rotated[i]));
    BOOST_CHECK_CLOSE_FRACTION(rotated[i].r, baseline[i].r, 2.e-7);
    BOOST_CHECK_CLOSE_FRACTION(rotated[i].varianceR, baseline[i].varianceR, 5.e-7);
  }
  CellDirectionCompatibility baselineCompatibility{};
  CellDirectionCompatibility rotatedCompatibility{};
  BOOST_REQUIRE(cellDirectionsAreCompatible(baseline, 100.f, baselineCompatibility));
  BOOST_REQUIRE(cellDirectionsAreCompatible(rotated, 100.f, rotatedCompatibility));
  BOOST_CHECK_CLOSE_FRACTION(rotatedCompatibility.chi2, baselineCompatibility.chi2, 1.e-6);
}

BOOST_AUTO_TEST_CASE(CollinearCylinderAndDiskObservationsHaveZeroChi2)
{
  std::array<DirectionObservation, 3> cylinder{};
  std::array<DirectionObservation, 3> disk{};
  const std::array<float, 3> radii{1.f, 2.f, 4.f};
  const std::array<float, 3> z{10.f, 20.f, 40.f};
  for (std::size_t i = 0; i < radii.size(); ++i) {
    const auto cylinderSurface = directionSurface(SurfaceKind::Cylinder, radii[i]);
    const auto cylinderMeasurement = cylinderDirectionMeasurement(radii[i], 0.f, z[i], 0.02f, 0.001f, 0.03f);
    BOOST_REQUIRE(makeDirectionObservation(cylinderSurface, cylinderMeasurement, cylinder[i]));

    const auto diskSurface = directionSurface(SurfaceKind::Disk, z[i]);
    const auto diskMeasurement = diskDirectionMeasurement(radii[i], 0.f, 0.03f, 0.001f, 0.02f);
    BOOST_REQUIRE(makeDirectionObservation(diskSurface, diskMeasurement, disk[i]));
  }
  CellDirectionCompatibility cylinderCompatibility{};
  CellDirectionCompatibility diskCompatibility{};
  BOOST_REQUIRE(cellDirectionsAreCompatible(cylinder, 1.f, cylinderCompatibility));
  BOOST_REQUIRE(cellDirectionsAreCompatible(disk, 1.f, diskCompatibility));
  BOOST_CHECK_SMALL(cylinderCompatibility.chi2, 1.e-12);
  BOOST_CHECK_SMALL(diskCompatibility.chi2, 1.e-12);
}

BOOST_AUTO_TEST_CASE(CellDirectionCompatibilityUsesOnlyCommonNSigmaCut)
{
  const std::array<DirectionObservation, 3> observations{{
    {2., 1., 0.04, 0.006, 0.09},
    {4., 4., 0.05, -0.004, 0.08},
    {7., 8., 0.03, 0.003, 0.07},
  }};
  CellDirectionCompatibility pass{};
  CellDirectionCompatibility fail{};
  BOOST_CHECK(cellDirectionsAreCompatible(observations, 0.4f, pass));
  BOOST_CHECK(!cellDirectionsAreCompatible(observations, 0.39f, fail));
  BOOST_CHECK_EQUAL(pass.residual, fail.residual);
  BOOST_CHECK_EQUAL(pass.variance, fail.variance);
  BOOST_CHECK_EQUAL(pass.chi2, fail.chi2);
}

BOOST_AUTO_TEST_CASE(CellDirectionCompatibilityFailsClosedTransactionally)
{
  std::array<DirectionObservation, 3> observations{{
    {2., 1., 0.04, 0.006, 0.09},
    {4., 4., 0.05, -0.004, 0.08},
    {7., 8., 0.03, 0.003, 0.07},
  }};
  const CellDirectionCompatibility sentinel{11., 12., 13.};
  auto compatibility = sentinel;
  observations[1].covarianceRZ = 1.; // non-PSD
  BOOST_CHECK(!cellDirectionsAreCompatible(observations, 5.f, compatibility));
  BOOST_CHECK(bitEqual(compatibility, sentinel));

  observations[1] = {4., 4., 0.05, -0.004, 0.08};
  observations[2].z = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK(!cellDirectionsAreCompatible(observations, 5.f, compatibility));
  BOOST_CHECK(bitEqual(compatibility, sentinel));

  auto invalidDisk = diskDirectionMeasurement(3.f, 4.f, 0.04f, 1.f, 0.09f);
  DirectionObservation observationSentinel{1., 2., 3., 4., 5.};
  auto observation = observationSentinel;
  BOOST_CHECK(!makeDiskDirectionObservation(
    directionSurface(SurfaceKind::Disk, -50.f), invalidDisk, observation));
  BOOST_CHECK(bitEqual(observation, observationSentinel));

  auto invalidCylinder = cylinderDirectionMeasurement(3.f, 4.f, -2.2f, 0.04f, 1.f, 0.09f);
  BOOST_CHECK(!makeCylinderDirectionObservation(
    directionSurface(SurfaceKind::Cylinder, 3.5f), invalidCylinder, observation));
  BOOST_CHECK(bitEqual(observation, observationSentinel));

  invalidCylinder = cylinderDirectionMeasurement(0.f, 0.f, -2.2f, 0.04f, 0.006f, 0.09f);
  BOOST_CHECK(!makeCylinderDirectionObservation(
    directionSurface(SurfaceKind::Cylinder, 3.5f), invalidCylinder, observation));
  BOOST_CHECK(bitEqual(observation, observationSentinel));
}

// ===========================================================================
// buildCylinderCellSeed
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
  TrackingKernelParameters params;
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(buildCylinderCellSeed(
    barrelMeasurementInner(), barrelMeasurementMiddle(), barrelMeasurementOuter(), barrelMaterial(), BarrelBz,
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
  TrackingKernelParameters accept;
  accept.maxChi2ClusterAttachment = lastStepChi2; // exact threshold: inclusive accept
  BOOST_CHECK(buildCylinderCellSeed(
    barrelMeasurementInner(), barrelMeasurementMiddle(), barrelMeasurementOuter(), barrelMaterial(), BarrelBz,
    1, o2::track::PID::Pion, outState, chi2, accept, reason));

  const auto sentinel = barrelAttachState();
  auto before = sentinel;
  float chi2Before = -123.f;
  auto rejectedState = before;
  float rejectedChi2 = chi2Before;
  TrackingKernelParameters reject;
  reject.maxChi2ClusterAttachment = std::nextafter(lastStepChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!buildCylinderCellSeed(
    barrelMeasurementInner(), barrelMeasurementMiddle(), barrelMeasurementOuter(), barrelMaterial(), BarrelBz,
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
    TrackingKernelParameters params;
    params.maxChi2ClusterAttachment = 1.e6f;
    BOOST_CHECK(!buildCylinderCellSeed(
      barrelMeasurementFor(clusterInner, barrelHitInner()), barrelMeasurementFor(barrelClusterMiddle(), hitMiddle), barrelMeasurementOuter(),
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
  // (mirrors testCellFinding.cxx's own rotation-failure
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
    TrackingKernelParameters params;
    params.maxChi2ClusterAttachment = 1.e6f;
    BOOST_CHECK(!buildCylinderCellSeed(
      barrelMeasurementInner(), barrelMeasurementFor(barrelClusterMiddle(), farHit), barrelMeasurementOuter(),
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
  TrackingKernelParameters params;
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(buildCylinderCellSeed(
    barrelMeasurementInner(), barrelMeasurementMiddle(), barrelMeasurementOuter(), barrelMaterial(), BarrelBz, 1, o2::track::PID::Kaon, first, chi2First, params, reason));
  BOOST_REQUIRE(buildCylinderCellSeed(
    barrelMeasurementInner(), barrelMeasurementMiddle(), barrelMeasurementOuter(), barrelMaterial(), BarrelBz, 1, o2::track::PID::Kaon, second, chi2Second, params, reason));
  BOOST_CHECK(bitEqual(first, second));
  BOOST_CHECK_EQUAL(chi2First, chi2Second);
}

// ===========================================================================
// buildDiskCellSeed
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
  TrackingKernelParameters params;
  params.trackletMinPt = DiskTrackletMinPt;
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(buildDiskCellSeed(
    diskMeasurementInner(), diskMeasurementMiddle(), diskMeasurementOuter(), diskMaterial(), DiskBz,
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
  TrackingKernelParameters accept;
  accept.trackletMinPt = DiskTrackletMinPt;
  accept.maxChi2ClusterAttachment = lastStepChi2;
  BOOST_CHECK(buildDiskCellSeed(
    diskMeasurementInner(), diskMeasurementMiddle(), diskMeasurementOuter(), diskMaterial(), DiskBz,
    2, o2::track::PID::Pion, outState, chi2, accept, reason));

  const auto sentinel = diskAttachState();
  auto rejectedState = sentinel;
  float rejectedChi2 = -123.f;
  const auto before = rejectedState;
  const float chi2Before = rejectedChi2;
  TrackingKernelParameters reject;
  reject.trackletMinPt = DiskTrackletMinPt;
  reject.maxChi2ClusterAttachment = std::nextafter(lastStepChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!buildDiskCellSeed(
    diskMeasurementInner(), diskMeasurementMiddle(), diskMeasurementOuter(), diskMaterial(), DiskBz,
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
    TrackingKernelParameters params;
    params.trackletMinPt = DiskTrackletMinPt;
    params.maxChi2ClusterAttachment = 1.e6f;
    BOOST_CHECK(!buildDiskCellSeed(
      diskMeasurementOuter() /* swapped: inner<=outer violates the strict boundary */, diskMeasurementMiddle(), diskMeasurementInner(),
      diskMaterial(), DiskBz,
      2, o2::track::PID::Pion, state, chi2, params, reason));
    BOOST_CHECK(reason == OperationFailureReason::SeedGeometryDegenerate);
    BOOST_CHECK(bitEqual(state, before));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  }

  // Propagation-stage failure: a non-finite inner reference z (measurement.
  // frame.q) is rejected by forward::propagate's own finite-target
  // precondition before any commit.
  {
    auto badMeasurementInner = diskMeasurementInner();
    badMeasurementInner.frame.q = std::numeric_limits<float>::infinity();
    auto state = sentinel;
    float chi2 = -321.f;
    const auto before = state;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    TrackingKernelParameters params;
    params.trackletMinPt = DiskTrackletMinPt;
    params.maxChi2ClusterAttachment = 1.e6f;
    BOOST_CHECK(!buildDiskCellSeed(
      badMeasurementInner, diskMeasurementMiddle(), diskMeasurementOuter(), diskMaterial(), DiskBz,
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
    TrackingKernelParameters params;
    params.trackletMinPt = DiskTrackletMinPt;
    params.maxChi2ClusterAttachment = 1.e6f;
    BOOST_CHECK(!buildDiskCellSeed(
      diskMeasurementInner(), diskMeasurementMiddle(), diskMeasurementOuter(), badMaterial, DiskBz,
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
  TrackingKernelParameters params;
  params.trackletMinPt = DiskTrackletMinPt;
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(buildDiskCellSeed(
    diskMeasurementInner(), diskMeasurementMiddle(), diskMeasurementOuter(),
    diskMaterial(), DiskBz, 2, o2::track::PID::Pion, first, chi2First, params, reason));
  BOOST_REQUIRE(buildDiskCellSeed(
    diskMeasurementInner(), diskMeasurementMiddle(), diskMeasurementOuter(),
    diskMaterial(), DiskBz, 2, o2::track::PID::Pion, second, chi2Second, params, reason));
  BOOST_CHECK(bitEqual(first, second));
  BOOST_CHECK_EQUAL(chi2First, chi2Second);
}

BOOST_AUTO_TEST_CASE(BuildCellSeedDiskActivatesEnergyLossUnlikeLegacyMcsOnlyPath)
{
  // Intentional physics difference (kickoff message): unlike the legacy,
  // MCS-only buildDiskCellSeed, the native path reads
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
  TrackingKernelParameters params;
  params.trackletMinPt = DiskTrackletMinPt;
  params.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(buildDiskCellSeed(
    diskMeasurementInner(), diskMeasurementMiddle(), diskMeasurementOuter(),
    noEnergyLoss, DiskBz, 2, o2::track::PID::Pion, stateNoLoss, chi2NoLoss, params, reason));
  BOOST_REQUIRE(buildDiskCellSeed(
    diskMeasurementInner(), diskMeasurementMiddle(), diskMeasurementOuter(),
    withEnergyLoss, DiskBz, 2, o2::track::PID::Pion, stateWithLoss, chi2WithLoss, params, reason));

  BOOST_CHECK_NE(stateNoLoss.parameters[4], stateWithLoss.parameters[4]);
}

// ===========================================================================
// cell compatibility leaves
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

  TrackingKernelParameters accept;
  accept.maxChi2ClusterAttachment = refChi2;
  BOOST_CHECK(cellsCylinderAreCompatible(current, next, -1, -1, CompatBz, accept));

  TrackingKernelParameters reject;
  reject.maxChi2ClusterAttachment = std::nextafter(refChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!cellsCylinderAreCompatible(current, next, -1, -1, CompatBz, reject));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleBarrelRotationAndPropagationFailuresAreRejectedNotThrown)
{
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;

  // Rotation failure: far frame angle.
  {
    auto current = compatBarrelCurrent();
    current.alpha = 0.f;
    auto next = compatBarrelNext();
    next.alpha = 3.f;
    BOOST_CHECK(!cellsCylinderAreCompatible(current, next, -1, -1, CompatBz, permissive));
  }

  // Propagation failure: extreme q2pt on `next` makes the requested step
  // unreachable.
  {
    auto current = compatBarrelCurrent();
    auto next = compatBarrelNext();
    next.parameters[4] = 2000.f;
    BOOST_CHECK(!cellsCylinderAreCompatible(current, next, -1, -1, CompatBz, permissive));
  }
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleBarrelFamilyMismatchFailsClosed)
{
  auto current = compatBarrelCurrent();
  current.family = StateFamily::Forward;
  const auto next = compatBarrelNext();
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_CHECK(!cellsCylinderAreCompatible(current, next, -1, -1, CompatBz, permissive));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleBarrelIgnoresClusterIndicesAndNeverMutatesInputs)
{
  const auto current = compatBarrelCurrent();
  const auto next = compatBarrelNext();
  const auto currentBefore = current;
  const auto nextBefore = next;
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;

  const bool r1 = cellsCylinderAreCompatible(current, next, 10, 10, CompatBz, permissive);
  const bool r2 = cellsCylinderAreCompatible(current, next, 10, 999, CompatBz, permissive);
  BOOST_CHECK_EQUAL(r1, r2);
  BOOST_CHECK(bitEqual(current, currentBefore));
  BOOST_CHECK(bitEqual(next, nextBefore));

  const bool r3 = cellsCylinderAreCompatible(current, next, 10, 10, CompatBz, permissive);
  BOOST_CHECK_EQUAL(r1, r3);
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleDiskAcceptsAtExactThresholdAndRejectsBelow)
{
  const auto current = compatDiskCurrent();
  const auto next = compatDiskNext();

  auto reference = next;
  OperationFailureReason reason{};
  BOOST_REQUIRE(Propagator::propagateForward(reference, current.referenceCoordinate, 0.f, reason));
  float refChi2 = 0.f;
  BOOST_REQUIRE(forward::stateChi2(current, reference, refChi2, reason));
  BOOST_REQUIRE_GT(refChi2, 0.f);

  TrackingKernelParameters accept;
  accept.maxChi2ClusterAttachment = refChi2;
  BOOST_CHECK(cellsDiskAreCompatible(current, next, 42, 42, 0.f, accept));

  TrackingKernelParameters reject;
  reject.maxChi2ClusterAttachment = std::nextafter(refChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!cellsDiskAreCompatible(current, next, 42, 42, 0.f, reject));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleDiskClusterIndexContinuityRejectedBeforeAnyPropagation)
{
  const auto current = compatDiskCurrent();
  const auto next = compatDiskNext();
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  // Otherwise-compatible states, but the temporary raw cluster-index
  // continuity input (see the header doc) rejects first.
  BOOST_CHECK(!cellsDiskAreCompatible(current, next, 5, 6, 0.f, permissive));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleDiskPropagationFailureIsRejectedNotThrown)
{
  auto current = compatDiskCurrent();
  current.referenceCoordinate = -5.f; // dz != 0 relative to `next`
  auto next = compatDiskNext();
  next.parameters[3] = 0.f; // tanl == 0: propagateLinear's UnreachableTarget
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_CHECK(!cellsDiskAreCompatible(current, next, 7, 7, 0.f, permissive));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleDiskFamilyMismatchFailsClosed)
{
  auto current = compatDiskCurrent();
  current.family = StateFamily::Barrel;
  const auto next = compatDiskNext();
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_CHECK(!cellsDiskAreCompatible(current, next, 3, 3, 0.f, permissive));
}

BOOST_AUTO_TEST_CASE(CellsAreCompatibleDiskInputsUntouchedAndByteDeterministic)
{
  const auto current = compatDiskCurrent();
  const auto next = compatDiskNext();
  const auto currentBefore = current;
  const auto nextBefore = next;
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  const bool r1 = cellsDiskAreCompatible(current, next, 1, 1, 0.f, permissive);
  const bool r2 = cellsDiskAreCompatible(current, next, 1, 1, 0.f, permissive);
  BOOST_CHECK_EQUAL(r1, r2);
  BOOST_CHECK(bitEqual(current, currentBefore));
  BOOST_CHECK(bitEqual(next, nextBefore));
}

// ===========================================================================
// attachCylinderHit
// ===========================================================================

BOOST_AUTO_TEST_CASE(AttachHitBarrelSuccessAndExactChi2Threshold)
{
  const auto state0 = barrelAttachState();
  const auto hit = barrelMeasurementFromHit(barrelAttachHit());
  const auto material = barrelAttachMaterial();

  auto probe = state0;
  float probeChi2 = 0.f;
  OperationFailureReason reason{};
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachCylinderHit(probe, hit, material, BarrelAttachBz, probeChi2, permissive, reason));
  BOOST_REQUIRE_GT(probeChi2, 0.f);

  auto accepted = state0;
  float acceptedChi2 = 0.f;
  TrackingKernelParameters accept;
  accept.maxChi2ClusterAttachment = probeChi2;
  BOOST_CHECK(attachCylinderHit(accepted, hit, material, BarrelAttachBz, acceptedChi2, accept, reason));
  BOOST_CHECK(bitEqual(accepted, probe));
  BOOST_CHECK_EQUAL(acceptedChi2, probeChi2);

  auto rejected = state0;
  float rejectedChi2 = -1.f;
  const auto before = rejected;
  const float chi2Before = rejectedChi2;
  TrackingKernelParameters reject;
  reject.maxChi2ClusterAttachment = std::nextafter(probeChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!attachCylinderHit(rejected, hit, material, BarrelAttachBz, rejectedChi2, reject, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(rejected, before));
  BOOST_CHECK_EQUAL(rejectedChi2, chi2Before);
}

BOOST_AUTO_TEST_CASE(AttachHitBarrelEachFailureStagePreservesStateTransactionally)
{
  const auto state0 = barrelAttachState();
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;

  auto checkFailure = [&](const SurfaceMeasurement& measurement, const NominalSurfaceMaterial& material,
                          OperationFailureReason expectedOneOf1, OperationFailureReason expectedOneOf2) {
    auto state = state0;
    float chi2 = -1.f;
    const auto before = state;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    BOOST_CHECK(!attachCylinderHit(state, measurement, material, BarrelAttachBz, chi2, permissive, reason));
    BOOST_CHECK(reason == expectedOneOf1 || reason == expectedOneOf2);
    BOOST_CHECK(bitEqual(state, before));
    BOOST_CHECK_EQUAL(chi2, chi2Before);
  };

  // Rotation failure.
  {
    auto farHit = barrelAttachHit();
    farHit.alphaTrackingFrame = state0.alpha + 3.f;
    checkFailure(barrelMeasurementFromHit(farHit), barrelAttachMaterial(), OperationFailureReason::RotationFailure, OperationFailureReason::RotationFailure);
  }

  // Propagation failure.
  {
    auto farHit = barrelAttachHit();
    farHit.xTrackingFrame = -50000.f;
    checkFailure(barrelMeasurementFromHit(farHit), barrelAttachMaterial(), OperationFailureReason::UnreachableTarget, OperationFailureReason::PropagationFailure);
  }

  // Material failure.
  {
    const NominalSurfaceMaterial badMaterial{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
    checkFailure(barrelMeasurementFromHit(barrelAttachHit()), badMaterial, OperationFailureReason::MaterialFailure, OperationFailureReason::MaterialFailure);
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
  // attachCylinderHit already applies today. barrel::update shares
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
  const auto measurement = barrelMeasurementFromHit(hit);
  const NominalSurfaceMaterial noopMaterial{0.f, 0.f};
  const auto before = state;
  float chi2 = -1.f;
  const float chi2Before = chi2;
  OperationFailureReason reason{};
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_CHECK(!attachCylinderHit(state, measurement, noopMaterial, BarrelAttachBz, chi2, permissive, reason));
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
  const auto hit = barrelMeasurementFromHit(barrelAttachHit());
  const auto material = barrelAttachMaterial();

  float neutralChi2 = 0.f;
  float chargedChi2 = 0.f;
  OperationFailureReason reason{};
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachCylinderHit(neutral, hit, material, BarrelAttachBz, neutralChi2, permissive, reason));
  BOOST_REQUIRE(attachCylinderHit(charged, hit, material, BarrelAttachBz, chargedChi2, permissive, reason));
  BOOST_CHECK(!bitEqual(neutral, charged));
}

BOOST_AUTO_TEST_CASE(AttachHitBarrelIsByteDeterministic)
{
  auto first = barrelAttachState();
  auto second = barrelAttachState();
  float chi2First = 0.f;
  float chi2Second = 0.f;
  OperationFailureReason reason{};
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachCylinderHit(first, barrelMeasurementFromHit(barrelAttachHit()), barrelAttachMaterial(), BarrelAttachBz, chi2First, permissive, reason));
  BOOST_REQUIRE(attachCylinderHit(second, barrelMeasurementFromHit(barrelAttachHit()), barrelAttachMaterial(), BarrelAttachBz, chi2Second, permissive, reason));
  BOOST_CHECK(bitEqual(first, second));
  BOOST_CHECK_EQUAL(chi2First, chi2Second);
}

// ===========================================================================
// attachDiskHit
// ===========================================================================

BOOST_AUTO_TEST_CASE(AttachHitDiskSuccessAndExactChi2Threshold)
{
  const auto state0 = diskAttachState();
  const auto hit = diskAttachMeasurement();
  const auto material = diskAttachMaterial();

  auto probe = state0;
  float probeChi2 = 0.f;
  OperationFailureReason reason{};
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachDiskHit(probe, hit, material, DiskAttachBz, probeChi2, permissive, reason));
  BOOST_REQUIRE_GT(probeChi2, 0.f);

  auto accepted = state0;
  float acceptedChi2 = 0.f;
  TrackingKernelParameters accept;
  accept.maxChi2ClusterAttachment = probeChi2;
  BOOST_CHECK(attachDiskHit(accepted, hit, material, DiskAttachBz, acceptedChi2, accept, reason));
  BOOST_CHECK(bitEqual(accepted, probe));
  BOOST_CHECK_EQUAL(acceptedChi2, probeChi2);

  auto rejected = state0;
  float rejectedChi2 = -1.f;
  const auto before = rejected;
  const float chi2Before = rejectedChi2;
  TrackingKernelParameters reject;
  reject.maxChi2ClusterAttachment = std::nextafter(probeChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!attachDiskHit(rejected, hit, material, DiskAttachBz, rejectedChi2, reject, reason));
  BOOST_CHECK(reason == OperationFailureReason::PredictedChi2Failure);
  BOOST_CHECK(bitEqual(rejected, before));
  BOOST_CHECK_EQUAL(rejectedChi2, chi2Before);
}

BOOST_AUTO_TEST_CASE(AttachHitDiskEachFailureStagePreservesStateTransactionally)
{
  const auto state0 = diskAttachState();
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;

  // Propagation failure: tanl == 0 at zero field rejects with
  // UnreachableTarget.
  {
    auto zeroTanl = state0;
    zeroTanl.parameters[3] = 0.f;
    auto hit = diskAttachHit();
    hit.zCoordinate = -60.f; // dz != 0
    const auto measurement = diskAttachMeasurementFrom(hit);
    auto state = zeroTanl;
    float chi2 = -1.f;
    const auto before = state;
    const float chi2Before = chi2;
    OperationFailureReason reason{};
    BOOST_CHECK(!attachDiskHit(state, measurement, diskAttachMaterial(), 0.f, chi2, permissive, reason));
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
    BOOST_CHECK(!attachDiskHit(state, diskAttachMeasurement(), badMaterial, DiskAttachBz, chi2, permissive, reason));
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
  const auto hit = diskAttachMeasurement();
  const auto material = diskAttachMaterial();

  float neutralChi2 = 0.f;
  float chargedChi2 = 0.f;
  OperationFailureReason reason{};
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachDiskHit(neutral, hit, material, DiskAttachBz, neutralChi2, permissive, reason));
  BOOST_REQUIRE(attachDiskHit(charged, hit, material, DiskAttachBz, chargedChi2, permissive, reason));
  BOOST_CHECK(!bitEqual(neutral, charged));
}

BOOST_AUTO_TEST_CASE(AttachHitDiskActivatesEnergyLossUnlikeLegacyMcsOnlyPath)
{
  auto noLossMaterial = diskAttachMaterial();
  noLossMaterial.arealDensityGPerCm2 = 0.f;
  auto withLossMaterial = diskAttachMaterial();

  auto stateNoLoss = diskAttachState();
  auto stateWithLoss = diskAttachState();
  const auto hit = diskAttachMeasurement();
  float chi2NoLoss = 0.f;
  float chi2WithLoss = 0.f;
  OperationFailureReason reason{};
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachDiskHit(stateNoLoss, hit, noLossMaterial, DiskAttachBz, chi2NoLoss, permissive, reason));
  BOOST_REQUIRE(attachDiskHit(stateWithLoss, hit, withLossMaterial, DiskAttachBz, chi2WithLoss, permissive, reason));
  BOOST_CHECK_NE(stateNoLoss.parameters[4], stateWithLoss.parameters[4]);
}

BOOST_AUTO_TEST_CASE(AttachHitDiskIsByteDeterministic)
{
  auto first = diskAttachState();
  auto second = diskAttachState();
  float chi2First = 0.f;
  float chi2Second = 0.f;
  OperationFailureReason reason{};
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachDiskHit(first, diskAttachMeasurement(), diskAttachMaterial(), DiskAttachBz, chi2First, permissive, reason));
  BOOST_REQUIRE(attachDiskHit(second, diskAttachMeasurement(), diskAttachMaterial(), DiskAttachBz, chi2Second, permissive, reason));
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
  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  BOOST_REQUIRE(attachCylinderHit(stateLow, barrelMeasurementFromHit(lowCrossTerm), barrelAttachMaterial(), BarrelAttachBz, chi2Low, permissive, reason));
  BOOST_REQUIRE(attachCylinderHit(stateHigh, barrelMeasurementFromHit(highCrossTerm), barrelAttachMaterial(), BarrelAttachBz, chi2High, permissive, reason));
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

  TrackingKernelParameters permissive;
  permissive.maxChi2ClusterAttachment = 1.e6f;
  OperationFailureReason reason{};

  auto stateBaseline = state0;
  float chi2Baseline = 0.f;
  BOOST_REQUIRE(attachDiskHit(stateBaseline, diskAttachMeasurementFrom(baseline), diskAttachMaterial(), DiskAttachBz, chi2Baseline, permissive, reason));

  auto stateCrossTerm = state0;
  float chi2CrossTerm = 0.f;
  BOOST_REQUIRE(attachDiskHit(stateCrossTerm, diskAttachMeasurementFrom(varyingCrossTerm), diskAttachMaterial(), DiskAttachBz, chi2CrossTerm, permissive, reason));
  BOOST_CHECK(bitEqual(stateBaseline, stateCrossTerm));
  BOOST_CHECK_EQUAL(chi2Baseline, chi2CrossTerm);

  auto stateUnreadFields = state0;
  float chi2UnreadFields = 0.f;
  BOOST_REQUIRE(attachDiskHit(stateUnreadFields, diskAttachMeasurementFrom(varyingUnreadFields), diskAttachMaterial(), DiskAttachBz, chi2UnreadFields, permissive, reason));
  BOOST_CHECK(bitEqual(stateBaseline, stateUnreadFields));
  BOOST_CHECK_EQUAL(chi2Baseline, chi2UnreadFields);
}
