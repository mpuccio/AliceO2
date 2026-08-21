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

#include "CommonConstants/MathConstants.h"
#include "ITSMFTTracking/detail/SurfaceStateOperations.h"
#include "ITSMFTTracking/MaterialPhysics.h"
#include "ITSMFTTracking/Propagator.h"
#include "ITSMFTTracking/detail/CandidateFinding.h"
#include "ITSMFTTracking/detail/DirectionCompatibility.h"
#include "ITStracking/Cluster.h"

/// Focused coverage for the explicit cylinder/disk SurfaceKinematicState
/// compatibility and hit-attachment leaves
/// (CandidateFinding.h/.cxx). The leaves are production-wired through
/// TrackerTraits; this test exercises their numerical contracts directly,
/// while the separate orchestration tests cover their TrackerTraits callers.
///
/// Oracle strategy (see the Stage-B design report's "Precision/oracles"
/// split, referenced from CandidateFinding.h):
///  - Formula-preserving evidence (rotation/propagation, predicted
///    chi2/update, state compatibility) is checked by an
///    independent, hand-written re-transcription of each operation's own
///    documented call sequence, built directly on the already-oracle-tested
///    detail::barrel::/detail::forward:: primitives (BarrelSurfaceStateOperations.h,
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

bool attachMeasurement(SurfaceKinematicState& state, const SurfaceMeasurement& measurement,
                       NominalSurfaceMaterial material, float bz, float& chi2,
                       const TrackingKernelParameters& parameters,
                       OperationFailureReason& reason)
{
  return Propagator::attachMeasurement(state, measurement, material, bz,
                                       material::MaterialTraversalDirection::OppositeMomentum, true,
                                       parameters.maxChi2ClusterAttachment, chi2, reason);
}

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

SurfaceMeasurement diskMeasurementFromHit(const o2::its::TrackingFrameInfo& hit)
{
  SurfaceMeasurement measurement{};
  measurement.frame = {hit.zCoordinate, hit.xCoordinate, hit.yCoordinate, 0.f};
  measurement.covariance.uu = hit.covarianceTrackingFrame[0];
  measurement.covariance.uv = 0.f;
  measurement.covariance.vv = hit.covarianceTrackingFrame[2];
  return measurement;
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
  state.kind = SurfaceKind::Cylinder;
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
  state.kind = SurfaceKind::Disk;
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

GlobalMeasurement diskDirectionMeasurement(float x, float y,
                                           float varianceX, float covarianceXY, float varianceY)
{
  GlobalMeasurement measurement{};
  measurement.position = {x, y, 1234.f};
  measurement.radius = std::hypot(x, y);
  measurement.covariance = {varianceX, covarianceXY, 0.f, varianceY, 0.f, 0.f};
  return measurement;
}

GlobalMeasurement cylinderDirectionMeasurement(float q, float u, float z,
                                               float varianceU, float covarianceUZ, float varianceZ,
                                               float frameAngle = 0.f)
{
  const float sine = std::sin(frameAngle);
  const float cosine = std::cos(frameAngle);
  GlobalMeasurement measurement{};
  measurement.position = {q * cosine - u * sine, q * sine + u * cosine, z};
  measurement.radius = std::hypot(measurement.position.x, measurement.position.y);
  measurement.covariance = {sine * sine * varianceU,
                            -sine * cosine * varianceU,
                            -sine * covarianceUZ,
                            cosine * cosine * varianceU,
                            cosine * covarianceUZ,
                            varianceZ};
  return measurement;
}

TransverseDirectionObservation rotateTransverseObservation(
  const TransverseDirectionObservation& observation, double angle)
{
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  return {
    cosine * observation.x - sine * observation.y,
    sine * observation.x + cosine * observation.y,
    cosine * cosine * observation.varianceX -
      2. * cosine * sine * observation.covarianceXY +
      sine * sine * observation.varianceY,
    cosine * sine * (observation.varianceX - observation.varianceY) +
      (cosine * cosine - sine * sine) * observation.covarianceXY,
    sine * sine * observation.varianceX +
      2. * cosine * sine * observation.covarianceXY +
      cosine * cosine * observation.varianceY};
}

float transverseTrackletPhi(const TransverseDirectionObservation& first,
                            const TransverseDirectionObservation& second)
{
  return std::atan2(first.y - second.y, first.x - second.x);
}

GlobalMeasurement rotateDiskMeasurement(const GlobalMeasurement& measurement, double angle)
{
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  const double x = measurement.position.x;
  const double y = measurement.position.y;
  const double varianceX = measurement.covariance.xx;
  const double covarianceXY = measurement.covariance.xy;
  const double varianceY = measurement.covariance.yy;
  auto rotated = measurement;
  rotated.position.x = static_cast<float>(cosine * x - sine * y);
  rotated.position.y = static_cast<float>(sine * x + cosine * y);
  rotated.radius = std::hypot(rotated.position.x, rotated.position.y);
  rotated.covariance.xx = static_cast<float>(cosine * cosine * varianceX -
                                             2. * cosine * sine * covarianceXY +
                                             sine * sine * varianceY);
  rotated.covariance.xy = static_cast<float>(cosine * sine * (varianceX - varianceY) +
                                             (cosine * cosine - sine * sine) * covarianceXY);
  rotated.covariance.yy = static_cast<float>(sine * sine * varianceX +
                                             2. * cosine * sine * covarianceXY +
                                             cosine * cosine * varianceY);
  return rotated;
}

} // namespace

// ===========================================================================
// cell direction compatibility
// ===========================================================================

BOOST_AUTO_TEST_CASE(TransverseCompatibilityMatchesFullAnalyticCovarianceOracle)
{
  const std::array<TransverseDirectionObservation, 3> observations{{
    {0., 0., 0.04, 0.006, 0.09},
    {2., 1., 0.05, -0.004, 0.08},
    {5., 3., 0.03, 0.003, 0.07},
  }};
  const float firstPhi = transverseTrackletPhi(observations[0], observations[1]);
  const float secondPhi = transverseTrackletPhi(observations[1], observations[2]);
  TransverseDirectionCompatibility compatibility{};
  BOOST_REQUIRE(trackletDirectionsAreTransverselyCompatible(
    observations, firstPhi, secondPhi, {0.0025}, 0.f, 0.3f, 1.f, compatibility));

  const double expectedDeltaPhi = std::remainder(
    static_cast<double>(firstPhi) - secondPhi,
    2. * static_cast<double>(o2::constants::math::PI));
  BOOST_CHECK_CLOSE_FRACTION(compatibility.deltaPhi, expectedDeltaPhi, 1.e-12);
  BOOST_CHECK_SMALL(compatibility.maximumBending, 1.e-15);
  BOOST_CHECK_CLOSE_FRACTION(compatibility.variance, 0.06164035502958581, 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(compatibility.chi2,
                             expectedDeltaPhi * expectedDeltaPhi / 0.06164035502958581,
                             1.e-12);
}

BOOST_AUTO_TEST_CASE(TransverseCompatibilityUsesExactTrackletMinPtChordEnvelope)
{
  constexpr double radius = 100.;
  constexpr float bz = 5.f;
  const float boundaryPt = static_cast<float>(
    std::abs(o2::constants::math::B2C * bz) * radius);
  const auto point = [](double angle) {
    return TransverseDirectionObservation{
      radius * std::cos(angle), radius * std::sin(angle),
      1.e-8, 0., 1.e-8};
  };
  const std::array<TransverseDirectionObservation, 3> observations{{point(0.), point(0.1), point(0.25)}};
  const float firstPhi = transverseTrackletPhi(observations[0], observations[1]);
  const float secondPhi = transverseTrackletPhi(observations[1], observations[2]);

  TransverseDirectionCompatibility boundary{};
  BOOST_REQUIRE(trackletDirectionsAreTransverselyCompatible(
    observations, firstPhi, secondPhi, {}, bz, boundaryPt, 5.f, boundary));
  BOOST_CHECK_CLOSE_FRACTION(std::abs(boundary.deltaPhi), 0.125, 1.e-6);
  BOOST_CHECK_CLOSE_FRACTION(boundary.maximumBending, 0.125, 1.e-6);
  BOOST_CHECK_SMALL(boundary.chi2, 1.e-6);

  TransverseDirectionCompatibility aboveThreshold{};
  BOOST_CHECK(!trackletDirectionsAreTransverselyCompatible(
    observations, firstPhi, secondPhi, {}, bz, 2.f * boundaryPt, 5.f, aboveThreshold));
  BOOST_CHECK_LT(aboveThreshold.maximumBending, boundary.maximumBending);
}

BOOST_AUTO_TEST_CASE(TransverseCompatibilityNormalizesExcessWithOutgoingMS)
{
  const double tangent = std::tan(0.1);
  const std::array<TransverseDirectionObservation, 3> observations{{
    {0., 0., 1.e-10, 0., 1.e-10},
    {1., 0., 1.e-10, 0., 1.e-10},
    {2., tangent, 1.e-10, 0., 1.e-10},
  }};
  const float firstPhi = transverseTrackletPhi(observations[0], observations[1]);
  const float secondPhi = transverseTrackletPhi(observations[1], observations[2]);
  TransverseDirectionCompatibility withoutScattering{};
  TransverseDirectionCompatibility withScattering{};
  BOOST_CHECK(!trackletDirectionsAreTransverselyCompatible(
    observations, firstPhi, secondPhi, {}, 0.f, 0.3f, 5.f, withoutScattering));
  BOOST_REQUIRE(trackletDirectionsAreTransverselyCompatible(
    observations, firstPhi, secondPhi, {0.001}, 0.f, 0.3f, 5.f, withScattering));
  BOOST_CHECK_CLOSE_FRACTION(std::abs(withScattering.deltaPhi), 0.1, 2.e-6);
  BOOST_CHECK_GT(withScattering.variance, withoutScattering.variance);
  BOOST_CHECK_LT(withScattering.chi2, 25.);
}

BOOST_AUTO_TEST_CASE(TransverseCompatibilityIsInvariantUnderGlobalZRotation)
{
  const std::array<TransverseDirectionObservation, 3> observations{{
    {0.3, -0.7, 0.04, 0.006, 0.09},
    {2.2, 1.1, 0.05, -0.004, 0.08},
    {5.4, 3.6, 0.03, 0.003, 0.07},
  }};
  std::array<TransverseDirectionObservation, 3> rotated{};
  for (std::size_t i = 0; i < observations.size(); ++i) {
    rotated[i] = rotateTransverseObservation(observations[i], 0.73);
  }
  TransverseDirectionCompatibility baseline{};
  TransverseDirectionCompatibility transformed{};
  const float firstPhi = transverseTrackletPhi(observations[0], observations[1]);
  const float secondPhi = transverseTrackletPhi(observations[1], observations[2]);
  const float rotatedFirstPhi = transverseTrackletPhi(rotated[0], rotated[1]);
  const float rotatedSecondPhi = transverseTrackletPhi(rotated[1], rotated[2]);
  BOOST_REQUIRE(trackletDirectionsAreTransverselyCompatible(
    observations, firstPhi, secondPhi, {0.0025}, 0.5f, 0.3f, 100.f, baseline));
  BOOST_REQUIRE(trackletDirectionsAreTransverselyCompatible(
    rotated, rotatedFirstPhi, rotatedSecondPhi, {0.0025}, 0.5f, 0.3f, 100.f, transformed));
  BOOST_CHECK_CLOSE_FRACTION(transformed.deltaPhi, baseline.deltaPhi, 2.e-6);
  BOOST_CHECK_CLOSE_FRACTION(transformed.maximumBending, baseline.maximumBending, 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(transformed.variance, baseline.variance, 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(transformed.chi2, baseline.chi2, 3.e-6);
}

BOOST_AUTO_TEST_CASE(DiskTransverseObservationRetainsGlobalXYCovariance)
{
  const auto measurement = diskDirectionMeasurement(3.f, 4.f, 0.04f, 0.006f, 0.09f);
  TransverseDirectionObservation observation{};
  BOOST_REQUIRE(makeTransverseDirectionObservation(measurement, observation));
  BOOST_CHECK_CLOSE_FRACTION(observation.x, 3., 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(observation.y, 4., 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(observation.varianceX, 0.04, 1.e-7);
  BOOST_CHECK_CLOSE_FRACTION(observation.covarianceXY, 0.006, 1.e-7);
  BOOST_CHECK_CLOSE_FRACTION(observation.varianceY, 0.09, 1.e-7);
}

BOOST_AUTO_TEST_CASE(CylinderTransverseObservationRotatesTangentialCovariance)
{
  constexpr float angle = 0.7f;
  const auto measurement = cylinderDirectionMeasurement(
    3.f, 4.f, -2.2f, 0.04f, 0.006f, 0.09f, angle);
  TransverseDirectionObservation observation{};
  BOOST_REQUIRE(makeTransverseDirectionObservation(measurement, observation));
  const double sine = std::sin(angle);
  const double cosine = std::cos(angle);
  BOOST_CHECK_CLOSE_FRACTION(observation.x, 3. * cosine - 4. * sine, 1.e-7);
  BOOST_CHECK_CLOSE_FRACTION(observation.y, 3. * sine + 4. * cosine, 1.e-7);
  BOOST_CHECK_CLOSE_FRACTION(observation.varianceX, sine * sine * 0.04, 1.e-7);
  BOOST_CHECK_CLOSE_FRACTION(observation.covarianceXY, -sine * cosine * 0.04, 1.e-7);
  BOOST_CHECK_CLOSE_FRACTION(observation.varianceY, cosine * cosine * 0.04, 1.e-7);
}

BOOST_AUTO_TEST_CASE(TransverseCompatibilityFailsClosedTransactionally)
{
  std::array<TransverseDirectionObservation, 3> observations{{
    {0., 0., 0.04, 0.006, 0.09},
    {2., 1., 0.05, -0.004, 0.08},
    {5., 3., 0.03, 0.003, 0.07},
  }};
  const float firstPhi = transverseTrackletPhi(observations[0], observations[1]);
  const float secondPhi = transverseTrackletPhi(observations[1], observations[2]);
  const TransverseDirectionCompatibility sentinel{11., 12., 13., 14.};
  auto compatibility = sentinel;
  observations[1].covarianceXY = 1.;
  BOOST_CHECK(!trackletDirectionsAreTransverselyCompatible(
    observations, firstPhi, secondPhi, {}, 0.5f, 0.3f, 5.f, compatibility));
  BOOST_CHECK(bitEqual(compatibility, sentinel));

  observations[1] = {2., 1., 0.05, -0.004, 0.08};
  observations[2].x = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK(!trackletDirectionsAreTransverselyCompatible(
    observations, firstPhi, secondPhi, {}, 0.5f, 0.3f, 5.f, compatibility));
  BOOST_CHECK(bitEqual(compatibility, sentinel));

  observations[2] = observations[1];
  BOOST_CHECK(!trackletDirectionsAreTransverselyCompatible(
    observations, firstPhi, secondPhi, {}, 0.5f, 0.3f, 5.f, compatibility));
  BOOST_CHECK(bitEqual(compatibility, sentinel));

  observations[2] = {5., 3., 0.03, 0.003, 0.07};
  BOOST_CHECK(!trackletDirectionsAreTransverselyCompatible(
    observations, firstPhi, secondPhi, {-0.01}, 0.5f, 0.3f, 5.f, compatibility));
  BOOST_CHECK(bitEqual(compatibility, sentinel));
}

BOOST_AUTO_TEST_CASE(CellDirectionCompatibilityMatchesFullAnalyticOracle)
{
  const std::array<DirectionObservation, 3> observations{{
    {2., 1., 0.04, 0.006, 0.09},
    {4., 4., 0.05, -0.004, 0.08},
    {7., 8., 0.03, 0.003, 0.07},
  }};
  CellDirectionCompatibility compatibility{};
  BOOST_REQUIRE(cellDirectionsAreCompatible(observations, {}, 1.f, compatibility));
  BOOST_CHECK_CLOSE_FRACTION(compatibility.residual, 1., 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(compatibility.variance, 6.55, 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(compatibility.chi2, 1. / 6.55, 1.e-12);
}

BOOST_AUTO_TEST_CASE(CellDirectionCompatibilityTranslatesAngularProcessNoiseIntoKVariance)
{
  const std::array<DirectionObservation, 3> observations{{
    {2., 1., 0.04, 0.006, 0.09},
    {4., 4., 0.05, -0.004, 0.08},
    {7., 8., 0.03, 0.003, 0.07},
  }};
  CellDirectionCompatibility measurementOnly{};
  CellDirectionCompatibility withScattering{};
  BOOST_CHECK(!cellDirectionsAreCompatible(observations, {}, 0.35f, measurementOnly));
  BOOST_REQUIRE(cellDirectionsAreCompatible(observations, {0.01}, 0.35f, withScattering));

  // Segment dot product = (2,3).(3,4) = 18, so the process contribution is
  // 18^2 * 0.01 = 3.24.
  BOOST_CHECK_CLOSE_FRACTION(withScattering.residual, 1., 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(withScattering.variance, 6.55 + 3.24, 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(withScattering.chi2, 1. / (6.55 + 3.24), 1.e-12);
}

BOOST_AUTO_TEST_CASE(CellDirectionCompatibilityPropagatesCommonBeamUncertainty)
{
  std::array<GlobalMeasurement, 3> measurements{{
    diskDirectionMeasurement(2.f, 0.f, 0.04f, 0.f, 0.04f),
    diskDirectionMeasurement(2.f, 3.f, 0.04f, 0.f, 0.04f),
    diskDirectionMeasurement(-4.f, 3.f, 0.04f, 0.f, 0.04f),
  }};
  measurements[0].z = 1.f;
  measurements[1].z = 4.f;
  measurements[2].z = 8.f;
  std::array<DirectionObservation, 3> observations{};
  for (std::size_t i = 0; i < measurements.size(); ++i) {
    BOOST_REQUIRE(makeDirectionObservation(measurements[i], observations[i]));
  }

  CellDirectionCompatibility exactBeam{};
  CellDirectionCompatibility uncertainBeam{};
  BOOST_REQUIRE(cellDirectionsAreCompatible(observations, {}, 0.f, 100.f, exactBeam));
  BOOST_REQUIRE(cellDirectionsAreCompatible(observations, {}, 0.01f, 100.f, uncertainBeam));
  BOOST_CHECK_GT(uncertainBeam.variance, exactBeam.variance);
  BOOST_CHECK_LT(uncertainBeam.chi2, exactBeam.chi2);
}

BOOST_AUTO_TEST_CASE(DiskDirectionObservationProjectsFullGlobalXYCovariance)
{
  auto measurement = diskDirectionMeasurement(3.f, 4.f, 0.04f, 0.006f, 0.09f);
  measurement.position.z = -50.f;
  DirectionObservation observation{};
  BOOST_REQUIRE(makeDirectionObservation(measurement, observation));
  BOOST_CHECK_CLOSE_FRACTION(observation.r, 5., 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(observation.z, -50., 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(observation.varianceR, 0.07776, 1.e-7);
  BOOST_CHECK_SMALL(observation.covarianceRZ, 1.e-15);
  BOOST_CHECK_SMALL(observation.varianceZ, 1.e-15);
}

BOOST_AUTO_TEST_CASE(CylinderDirectionObservationProjectsTrackingFrameCovariance)
{
  const auto measurement = cylinderDirectionMeasurement(3.f, 4.f, -2.2f, 0.04f, 0.006f, 0.09f);
  DirectionObservation observation{};
  BOOST_REQUIRE(makeDirectionObservation(measurement, observation));
  BOOST_CHECK_CLOSE_FRACTION(observation.r, 5., 1.e-12);
  BOOST_CHECK_CLOSE_FRACTION(observation.z, -2.2, 1.e-6);
  BOOST_CHECK_CLOSE_FRACTION(observation.varianceR, 0.0256, 1.e-7);
  BOOST_CHECK_CLOSE_FRACTION(observation.covarianceRZ, 0.0048, 1.e-7);
  BOOST_CHECK_CLOSE_FRACTION(observation.varianceZ, 0.09, 1.e-6);
}

BOOST_AUTO_TEST_CASE(DirectionObservationRepairsOnlyProjectionRoundoff)
{
  GlobalMeasurement measurement{};
  measurement.position = {-1.15416992f, -1.91446579f, -0.0102339825f};
  measurement.radius = 2.23546124f;
  measurement.covariance = {1.70975554e-7f, -1.03069610e-7f, 0.f,
                            6.21336937e-8f, -0.f, 4.68642071e-7f};

  DirectionObservation observation{};
  BOOST_REQUIRE(makeDirectionObservation(measurement, observation));
  BOOST_CHECK_EQUAL(observation.varianceR, 0.);
  BOOST_CHECK_EQUAL(observation.covarianceRZ, 0.);
  BOOST_CHECK_CLOSE_FRACTION(observation.varianceZ, 4.68642071e-7, 1.e-7);
}

BOOST_AUTO_TEST_CASE(DiskDirectionCompatibilityIsInvariantUnderGlobalZRotation)
{
  const std::array<SurfaceDescriptor, 3> surfaces{{
    directionSurface(SurfaceKind::Disk, 1.f),
    directionSurface(SurfaceKind::Disk, 4.f),
    directionSurface(SurfaceKind::Disk, 8.f),
  }};
  std::array<GlobalMeasurement, 3> measurements{{
    diskDirectionMeasurement(1.2f, 1.6f, 0.04f, 0.006f, 0.09f),
    diskDirectionMeasurement(-2.4f, 3.2f, 0.05f, -0.004f, 0.08f),
    diskDirectionMeasurement(4.2f, -5.6f, 0.03f, 0.003f, 0.07f),
  }};
  for (std::size_t i = 0; i < measurements.size(); ++i) {
    measurements[i].position.z = surfaces[i].referenceCoordinate;
  }
  std::array<DirectionObservation, 3> baseline{};
  std::array<DirectionObservation, 3> rotated{};
  for (std::size_t i = 0; i < measurements.size(); ++i) {
    BOOST_REQUIRE(makeDirectionObservation(measurements[i], baseline[i]));
    const auto rotatedMeasurement = rotateDiskMeasurement(measurements[i], 0.73);
    BOOST_REQUIRE(makeDirectionObservation(rotatedMeasurement, rotated[i]));
    BOOST_CHECK_CLOSE_FRACTION(rotated[i].r, baseline[i].r, 2.e-7);
    BOOST_CHECK_CLOSE_FRACTION(rotated[i].varianceR, baseline[i].varianceR, 5.e-7);
  }
  CellDirectionCompatibility baselineCompatibility{};
  CellDirectionCompatibility rotatedCompatibility{};
  BOOST_REQUIRE(cellDirectionsAreCompatible(baseline, {0.0025}, 100.f, baselineCompatibility));
  BOOST_REQUIRE(cellDirectionsAreCompatible(rotated, {0.0025}, 100.f, rotatedCompatibility));
  BOOST_CHECK_CLOSE_FRACTION(rotatedCompatibility.chi2, baselineCompatibility.chi2, 5.e-6);
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
    BOOST_REQUIRE(makeDirectionObservation(cylinderMeasurement, cylinder[i]));

    const auto diskSurface = directionSurface(SurfaceKind::Disk, z[i]);
    auto diskMeasurement = diskDirectionMeasurement(radii[i], 0.f, 0.03f, 0.001f, 0.02f);
    diskMeasurement.position.z = diskSurface.referenceCoordinate;
    BOOST_REQUIRE(makeDirectionObservation(diskMeasurement, disk[i]));
  }
  CellDirectionCompatibility cylinderCompatibility{};
  CellDirectionCompatibility diskCompatibility{};
  BOOST_REQUIRE(cellDirectionsAreCompatible(cylinder, {}, 1.f, cylinderCompatibility));
  BOOST_REQUIRE(cellDirectionsAreCompatible(disk, {}, 1.f, diskCompatibility));
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
  BOOST_CHECK(cellDirectionsAreCompatible(observations, {}, 0.4f, pass));
  BOOST_CHECK(!cellDirectionsAreCompatible(observations, {}, 0.39f, fail));
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
  BOOST_CHECK(!cellDirectionsAreCompatible(observations, {}, 5.f, compatibility));
  BOOST_CHECK(bitEqual(compatibility, sentinel));

  observations[1] = {4., 4., 0.05, -0.004, 0.08};
  observations[2].z = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK(!cellDirectionsAreCompatible(observations, {}, 5.f, compatibility));
  BOOST_CHECK(bitEqual(compatibility, sentinel));

  observations[2] = {7., 8., 0.03, 0.003, 0.07};
  BOOST_CHECK(!cellDirectionsAreCompatible(observations, {-0.01}, 5.f, compatibility));
  BOOST_CHECK(bitEqual(compatibility, sentinel));
  BOOST_CHECK(!cellDirectionsAreCompatible(
    observations, {std::numeric_limits<double>::quiet_NaN()}, 5.f, compatibility));
  BOOST_CHECK(bitEqual(compatibility, sentinel));

  auto invalidDisk = diskDirectionMeasurement(3.f, 4.f, 0.04f, 1.f, 0.09f);
  DirectionObservation observationSentinel{1., 2., 3., 4., 5.};
  auto observation = observationSentinel;
  BOOST_CHECK(!makeDirectionObservation(invalidDisk, observation));
  BOOST_CHECK(bitEqual(observation, observationSentinel));

  auto invalidCylinder = cylinderDirectionMeasurement(3.f, 4.f, -2.2f, 0.04f, 1.f, 0.09f);
  BOOST_CHECK(!makeDirectionObservation(invalidCylinder, observation));
  BOOST_CHECK(bitEqual(observation, observationSentinel));

  invalidCylinder = cylinderDirectionMeasurement(0.f, 0.f, -2.2f, 0.04f, 0.006f, 0.09f);
  BOOST_CHECK(!makeDirectionObservation(invalidCylinder, observation));
  BOOST_CHECK(bitEqual(observation, observationSentinel));
}

// ===========================================================================
// Measurement attachment
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
  BOOST_REQUIRE(attachMeasurement(probe, hit, material, BarrelAttachBz, probeChi2, permissive, reason));
  BOOST_REQUIRE_GT(probeChi2, 0.f);

  auto accepted = state0;
  float acceptedChi2 = 0.f;
  TrackingKernelParameters accept;
  accept.maxChi2ClusterAttachment = probeChi2;
  BOOST_CHECK(attachMeasurement(accepted, hit, material, BarrelAttachBz, acceptedChi2, accept, reason));
  BOOST_CHECK(bitEqual(accepted, probe));
  BOOST_CHECK_EQUAL(acceptedChi2, probeChi2);

  auto rejected = state0;
  float rejectedChi2 = -1.f;
  const auto before = rejected;
  const float chi2Before = rejectedChi2;
  TrackingKernelParameters reject;
  reject.maxChi2ClusterAttachment = std::nextafter(probeChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!attachMeasurement(rejected, hit, material, BarrelAttachBz, rejectedChi2, reject, reason));
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
    BOOST_CHECK(!attachMeasurement(state, measurement, material, BarrelAttachBz, chi2, permissive, reason));
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
  // attachCylinderHit already applies today. detail::barrel::update shares
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
  BOOST_CHECK(!attachMeasurement(state, measurement, noopMaterial, BarrelAttachBz, chi2, permissive, reason));
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
  BOOST_REQUIRE(attachMeasurement(neutral, hit, material, BarrelAttachBz, neutralChi2, permissive, reason));
  BOOST_REQUIRE(attachMeasurement(charged, hit, material, BarrelAttachBz, chargedChi2, permissive, reason));
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
  BOOST_REQUIRE(attachMeasurement(first, barrelMeasurementFromHit(barrelAttachHit()), barrelAttachMaterial(), BarrelAttachBz, chi2First, permissive, reason));
  BOOST_REQUIRE(attachMeasurement(second, barrelMeasurementFromHit(barrelAttachHit()), barrelAttachMaterial(), BarrelAttachBz, chi2Second, permissive, reason));
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
  BOOST_REQUIRE(attachMeasurement(probe, hit, material, DiskAttachBz, probeChi2, permissive, reason));
  BOOST_REQUIRE_GT(probeChi2, 0.f);

  auto accepted = state0;
  float acceptedChi2 = 0.f;
  TrackingKernelParameters accept;
  accept.maxChi2ClusterAttachment = probeChi2;
  BOOST_CHECK(attachMeasurement(accepted, hit, material, DiskAttachBz, acceptedChi2, accept, reason));
  BOOST_CHECK(bitEqual(accepted, probe));
  BOOST_CHECK_EQUAL(acceptedChi2, probeChi2);

  auto rejected = state0;
  float rejectedChi2 = -1.f;
  const auto before = rejected;
  const float chi2Before = rejectedChi2;
  TrackingKernelParameters reject;
  reject.maxChi2ClusterAttachment = std::nextafter(probeChi2, -std::numeric_limits<float>::infinity());
  BOOST_CHECK(!attachMeasurement(rejected, hit, material, DiskAttachBz, rejectedChi2, reject, reason));
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
    BOOST_CHECK(!attachMeasurement(state, measurement, diskAttachMaterial(), 0.f, chi2, permissive, reason));
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
    BOOST_CHECK(!attachMeasurement(state, diskAttachMeasurement(), badMaterial, DiskAttachBz, chi2, permissive, reason));
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
  BOOST_REQUIRE(attachMeasurement(neutral, hit, material, DiskAttachBz, neutralChi2, permissive, reason));
  BOOST_REQUIRE(attachMeasurement(charged, hit, material, DiskAttachBz, chargedChi2, permissive, reason));
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
  BOOST_REQUIRE(attachMeasurement(stateNoLoss, hit, noLossMaterial, DiskAttachBz, chi2NoLoss, permissive, reason));
  BOOST_REQUIRE(attachMeasurement(stateWithLoss, hit, withLossMaterial, DiskAttachBz, chi2WithLoss, permissive, reason));
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
  BOOST_REQUIRE(attachMeasurement(first, diskAttachMeasurement(), diskAttachMaterial(), DiskAttachBz, chi2First, permissive, reason));
  BOOST_REQUIRE(attachMeasurement(second, diskAttachMeasurement(), diskAttachMaterial(), DiskAttachBz, chi2Second, permissive, reason));
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
  BOOST_REQUIRE(attachMeasurement(stateLow, barrelMeasurementFromHit(lowCrossTerm), barrelAttachMaterial(), BarrelAttachBz, chi2Low, permissive, reason));
  BOOST_REQUIRE(attachMeasurement(stateHigh, barrelMeasurementFromHit(highCrossTerm), barrelAttachMaterial(), BarrelAttachBz, chi2High, permissive, reason));
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
  BOOST_REQUIRE(attachMeasurement(stateBaseline, diskAttachMeasurementFrom(baseline), diskAttachMaterial(), DiskAttachBz, chi2Baseline, permissive, reason));

  auto stateCrossTerm = state0;
  float chi2CrossTerm = 0.f;
  BOOST_REQUIRE(attachMeasurement(stateCrossTerm, diskAttachMeasurementFrom(varyingCrossTerm), diskAttachMaterial(), DiskAttachBz, chi2CrossTerm, permissive, reason));
  BOOST_CHECK(bitEqual(stateBaseline, stateCrossTerm));
  BOOST_CHECK_EQUAL(chi2Baseline, chi2CrossTerm);

  auto stateUnreadFields = state0;
  float chi2UnreadFields = 0.f;
  BOOST_REQUIRE(attachMeasurement(stateUnreadFields, diskAttachMeasurementFrom(varyingUnreadFields), diskAttachMaterial(), DiskAttachBz, chi2UnreadFields, permissive, reason));
  BOOST_CHECK(bitEqual(stateBaseline, stateUnreadFields));
  BOOST_CHECK_EQUAL(chi2Baseline, chi2UnreadFields);
}
