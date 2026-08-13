// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT TrackletFinding
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <type_traits>
#include <utility>

#include <boost/test/unit_test.hpp>

#include <TGeoGlobalMagField.h>

#include "DataFormatsITS/Vertex.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "Field/MagneticField.h"
#include "GPUCommonMath.h"
#include "ITSMFTTracking/detail/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/detail/CandidateFinding.h"
#include "ITSMFTTracking/detail/TrackerTraversalPreparation.h"
#include "ITStracking/TrackHelpers.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

struct PropagatorFieldFixture {
  PropagatorFieldFixture()
  {
    if (!TGeoGlobalMagField::Instance()->GetField()) {
      TGeoGlobalMagField::Instance()->SetField(o2::field::MagneticField::createNominalField(5, true));
      TGeoGlobalMagField::Instance()->Lock();
    }
  }
};

BOOST_GLOBAL_FIXTURE(PropagatorFieldFixture);

/// Focused numerical-parity coverage for the first D007 surface-kind boundary
/// operation migrated off the legacy per-detector branch (Architecture.md
/// §10, cellsAreCompatible). These tests do not exercise TrackerTraits'
/// production traversal -- see the handoff note on scope.

namespace
{

constexpr float Bz = 0.5f;

o2::its::TrackingFrameInfo makeBarrelHit(float xTF, float alpha, float y, float z, float sigma2Y = 1.e-4f, float sigma2Z = 1.e-4f)
{
  return o2::its::TrackingFrameInfo{xTF, y, z, xTF, alpha, {y, z}, {sigma2Y, 0.f, sigma2Z}};
}

o2::its::TrackingFrameInfo makeDiskHit(float z, float x, float y, float sigma2X = 1.e-2f, float sigma2Y = 1.e-2f)
{
  return o2::its::TrackingFrameInfo{x, y, z, 0.f, 0.f, {x, y}, {sigma2X, 0.f, sigma2Y}};
}

o2::its::Vertex makeVertex(float x, float y, float z,
                           float sigma2X, float sigma2Y, float sigma2Z,
                           unsigned short contributors = 1)
{
  const float position[3]{x, y, z};
  const float covariance[6]{sigma2X, 0.f, sigma2Y, 0.f, 0.f, sigma2Z};
  return o2::its::Vertex{position, covariance, contributors, 1.f};
}

o2::its::Cluster makeGlobalCluster(float x, float y, float z, int id = 0)
{
  return o2::its::Cluster{x, y, z, id};
}

GlobalMeasurement makeMeasurement(float x, float y, float z, float uu = 1.e-4f, float vv = 1.e-4f, float uv = 0.f)
{
  GlobalMeasurement measurement{};
  measurement.position = {x, y, z};
  measurement.radius = std::hypot(x, y);
  measurement.covariance = {uu, uv, 0.f, vv, 0.f, 0.f};
  return measurement;
}

GlobalMeasurement makeMeasurement(const o2::its::Cluster& cluster, float uu = 1.e-4f, float vv = 1.e-4f, float uv = 0.f)
{
  return makeMeasurement(cluster.xCoordinate, cluster.yCoordinate, cluster.zCoordinate, uu, vv, uv);
}

TrackletProjectionCache makeCylinderProjectionCache(int fromLayer, int toLayer, float fromRadius, float toRadius,
                                                    float targetMinR, float targetMaxR, float sourcePositionResolution,
                                                    float transitionMSAngle, float transitionPhiCut)
{
  return {fromLayer, toLayer, fromRadius, toRadius, targetMinR, targetMaxR, sourcePositionResolution,
          0.f, 0.f, transitionMSAngle, transitionPhiCut, false};
}

TrackletProjectionCache makeDiskProjectionCache(int fromLayer, int toLayer, float fromRadius,
                                                float fromReferenceCoordinate, float toReferenceCoordinate,
                                                float transitionMSAngle, float transitionPhiCut)
{
  return {fromLayer, toLayer, fromRadius, 0.f, 0.f, 0.f, 0.f,
          fromReferenceCoordinate, toReferenceCoordinate, transitionMSAngle, transitionPhiCut, true};
}

template <typename Window>
void checkSearchWindowEqual(const Window& lhs, const Window& rhs)
{
  BOOST_CHECK_EQUAL(lhs.bins.x, rhs.bins.x);
  BOOST_CHECK_EQUAL(lhs.bins.y, rhs.bins.y);
  BOOST_CHECK_EQUAL(lhs.bins.z, rhs.bins.z);
  BOOST_CHECK_EQUAL(lhs.bins.w, rhs.bins.w);
  if constexpr (std::is_same_v<Window, CylinderTrackletSearchWindow>) {
    BOOST_CHECK_EQUAL(lhs.tanLambda, rhs.tanLambda);
    BOOST_CHECK_EQUAL(lhs.sigmaZ, rhs.sigmaZ);
    BOOST_CHECK_EQUAL(lhs.phiCut, rhs.phiCut);
    BOOST_CHECK_EQUAL(lhs.nSigmaCut, rhs.nSigmaCut);
  } else {
    BOOST_CHECK_EQUAL(lhs.xProj, rhs.xProj);
    BOOST_CHECK_EQUAL(lhs.yProj, rhs.yProj);
    BOOST_CHECK_EQUAL(lhs.sigmaX, rhs.sigmaX);
    BOOST_CHECK_EQUAL(lhs.sigmaY, rhs.sigmaY);
    BOOST_CHECK_EQUAL(lhs.nSigmaCut, rhs.nSigmaCut);
  }
}

NominalSurfaceMaterial toMaterial(float xOverX0)
{
  return NominalSurfaceMaterial{xOverX0, xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho};
}

std::array<NominalSurfaceMaterial, 3> toMaterial(const std::array<float, 3>& xOverX0)
{
  return {toMaterial(xOverX0[0]), toMaterial(xOverX0[1]), toMaterial(xOverX0[2])};
}

std::vector<NominalSurfaceMaterial> toMaterial(const std::vector<float>& xOverX0)
{
  std::vector<NominalSurfaceMaterial> material;
  material.reserve(xOverX0.size());
  for (const float x0 : xOverX0) {
    material.push_back(toMaterial(x0));
  }
  return material;
}

TrackingKernelParameters makeKernelParameters(const TrackingParameters& params, SurfaceKind kind)
{
  (void)kind;
  TrackingKernelParameters out;
  out.trackletMinPt = params.TrackletMinPt;
  out.nSigmaCut = params.NSigmaCut;
  out.maxChi2ClusterAttachment = params.MaxChi2ClusterAttachment;
  out.maxChi2NDF = params.MaxChi2NDF;
  out.pvResolution = params.PVres;
  return out;
}

} // namespace

BOOST_AUTO_TEST_CASE(BindingCopiesEveryFieldToTheCorrectSlot)
{
  // Distinct sentinel per field so a field-swap bug in the binding is caught.
  TrackingParameters legacy;
  legacy.TrackletMinPt = 1.11f;
  legacy.NSigmaCut = 3.33f;
  legacy.MaxChi2ClusterAttachment = 4.44f;
  legacy.MaxChi2NDF = 5.55f;
  legacy.PVres = 8.88f;
  legacy.LayerxX0 = {0.011f, 0.022f, 0.033f};
  legacy.CorrType = o2::base::PropagatorF::MatCorrType::USEMatCorrLUT;

  const auto barrel = makeKernelParameters(legacy, SurfaceKind::Cylinder);
  BOOST_CHECK_CLOSE(barrel.trackletMinPt, 1.11f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.nSigmaCut, 3.33f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.maxChi2ClusterAttachment, 4.44f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.maxChi2NDF, 5.55f, 1e-6);
  BOOST_CHECK_CLOSE(barrel.pvResolution, 8.88f, 1e-6);
  BOOST_CHECK(barrel.isValid());

  const auto disk = makeKernelParameters(legacy, SurfaceKind::Disk);
  BOOST_CHECK_CLOSE(disk.trackletMinPt, 1.11f, 1e-6);
  BOOST_CHECK_CLOSE(disk.nSigmaCut, 3.33f, 1e-6);
  BOOST_CHECK_CLOSE(disk.maxChi2ClusterAttachment, 4.44f, 1e-6);
  BOOST_CHECK_CLOSE(disk.maxChi2NDF, 5.55f, 1e-6);
  BOOST_CHECK(disk.isValid());

  const auto legacyMaterial = toMaterial(legacy.LayerxX0);
  const auto attach = bindAttachHitConfig(gsl::span<const NominalSurfaceMaterial>(legacyMaterial), legacy);
  BOOST_REQUIRE_EQUAL(attach.layerMaterial.size(), 3u);
  BOOST_CHECK_CLOSE(attach.layerMaterial[0].xOverX0, 0.011f, 1e-6);
  BOOST_CHECK_CLOSE(attach.layerMaterial[1].xOverX0, 0.022f, 1e-6);
  BOOST_CHECK_CLOSE(attach.layerMaterial[2].xOverX0, 0.033f, 1e-6);
  BOOST_CHECK(attach.corrType == o2::base::PropagatorF::MatCorrType::USEMatCorrLUT);
  BOOST_CHECK(attach.isValid(3));
  BOOST_CHECK(!attach.isValid(4));
}

BOOST_AUTO_TEST_CASE(BoundNonFiniteParametersAreDetectableThroughIsValid)
{
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const float inf = std::numeric_limits<float>::infinity();

  TrackingParameters legacy;
  legacy.TrackletMinPt = 1.11f;
  legacy.NSigmaCut = 3.33f;
  legacy.MaxChi2ClusterAttachment = 4.44f;
  legacy.MaxChi2NDF = 5.55f;

  auto barrelHealthy = legacy;
  BOOST_CHECK(makeKernelParameters(barrelHealthy, SurfaceKind::Cylinder).isValid());
  auto barrelNaN = legacy;
  barrelNaN.MaxChi2ClusterAttachment = nan;
  BOOST_CHECK(!makeKernelParameters(barrelNaN, SurfaceKind::Cylinder).isValid());
  auto barrelInf = legacy;
  barrelInf.NSigmaCut = inf;
  BOOST_CHECK(!makeKernelParameters(barrelInf, SurfaceKind::Cylinder).isValid());

  auto diskHealthy = legacy;
  BOOST_CHECK(makeKernelParameters(diskHealthy, SurfaceKind::Disk).isValid());
  auto diskNaN = legacy;
  diskNaN.MaxChi2ClusterAttachment = nan;
  BOOST_CHECK(!makeKernelParameters(diskNaN, SurfaceKind::Disk).isValid());
  auto diskInf = legacy;
  diskInf.NSigmaCut = inf;
  BOOST_CHECK(!makeKernelParameters(diskInf, SurfaceKind::Disk).isValid());

  auto materialNaN = legacy;
  materialNaN.LayerxX0[2] = nan;
  const auto materialNaNSpan = toMaterial(materialNaN.LayerxX0);
  BOOST_CHECK(!bindAttachHitConfig(gsl::span<const NominalSurfaceMaterial>(materialNaNSpan), materialNaN)
                 .isValid(materialNaN.LayerxX0.size()));
  auto materialInf = legacy;
  materialInf.LayerxX0[4] = inf;
  const auto materialInfSpan = toMaterial(materialInf.LayerxX0);
  BOOST_CHECK(!bindAttachHitConfig(gsl::span<const NominalSurfaceMaterial>(materialInfSpan), materialInf)
                 .isValid(materialInf.LayerxX0.size()));
  auto invalidCorrection = legacy;
  invalidCorrection.CorrType = static_cast<o2::base::PropagatorF::MatCorrType>(99);
  const auto invalidCorrectionMaterial = toMaterial(invalidCorrection.LayerxX0);
  BOOST_CHECK(!bindAttachHitConfig(gsl::span<const NominalSurfaceMaterial>(invalidCorrectionMaterial), invalidCorrection)
                 .isValid(invalidCorrection.LayerxX0.size()));
}

BOOST_AUTO_TEST_CASE(CylinderProjectSearchWindowMatchesInlineFormulaAndDirectPhiZBins)
{
  TrackingParameters legacy;
  legacy.PVres = 0.f; // valid: disables only the primary-vertex resolution term
  const auto params = makeKernelParameters(legacy, SurfaceKind::Cylinder);
  BOOST_REQUIRE(params.isValid());

  IndexTableUtilsCore indexUtils;
  indexUtils.setTrackingParameters(legacy);

  const auto source = makeGlobalCluster(2.f, 0.f, 0.5f);
  const auto sourceMeasurement = makeMeasurement(source);
  const auto vertex = makeVertex(0.f, 0.f, 0.f, 1.e-4f, 1.e-4f, 4.e-4f, 4);
  const auto state = makeCylinderProjectionCache(0, 3, 2.f, 4.f, 3.8f, 4.2f, 5.e-4f, 2.e-3f, 0.08f);

  CylinderTrackletSearchWindow window{};
  BOOST_REQUIRE((projectCylinderSearchWindow(
    sourceMeasurement, source, vertex, state, Bz, indexUtils, params, window)));

  const float inverseR0 = 1.f / source.radius;
  const float resolution = o2::gpu::CAMath::Sqrt(o2::its::math_utils::Sq(state.sourcePositionResolution) +
                                                 o2::its::math_utils::Sq(params.pvResolution) / float(vertex.getNContributors()));
  const float tanLambda = (source.zCoordinate - vertex.getZ()) * inverseR0;
  const float zAtTargetMinR = tanLambda * (state.targetMinR - source.radius) + source.zCoordinate;
  const float zAtTargetMaxR = tanLambda * (state.targetMaxR - source.radius) + source.zCoordinate;
  const float sqInvDeltaZ0 = 1.f / (o2::its::math_utils::Sq(source.zCoordinate - vertex.getZ()) + o2::its::constants::Tolerance);
  const float sigmaZ = o2::gpu::CAMath::Sqrt((o2::its::math_utils::Sq(resolution) * o2::its::math_utils::Sq(tanLambda) *
                                              ((o2::its::math_utils::Sq(inverseR0) + sqInvDeltaZ0) * o2::its::math_utils::Sq(state.toRadius - state.fromRadius) + 1.f)) +
                                             o2::its::math_utils::Sq((state.toRadius - state.fromRadius) * state.transitionMSAngle));
  const auto directBins = getBinsPhiZ(source.phi, state.toLayer, zAtTargetMinR, zAtTargetMaxR,
                                      sigmaZ * params.nSigmaCut, state.transitionPhiCut, indexUtils);

  BOOST_CHECK_EQUAL(window.bins.x, directBins.x);
  BOOST_CHECK_EQUAL(window.bins.y, directBins.y);
  BOOST_CHECK_EQUAL(window.bins.z, directBins.z);
  BOOST_CHECK_EQUAL(window.bins.w, directBins.w);
  BOOST_CHECK_EQUAL(window.tanLambda, tanLambda);
  BOOST_CHECK_EQUAL(window.sigmaZ, sigmaZ);

  legacy.PVres = 0.025f;
  const auto positivePVParams = makeKernelParameters(legacy, SurfaceKind::Cylinder);
  BOOST_REQUIRE(positivePVParams.isValid());
  CylinderTrackletSearchWindow positivePVWindow{};
  BOOST_REQUIRE((projectCylinderSearchWindow(
    sourceMeasurement, source, vertex, state, Bz, indexUtils, positivePVParams, positivePVWindow)));
  const float positivePVResolution = o2::gpu::CAMath::Sqrt(o2::its::math_utils::Sq(state.sourcePositionResolution) +
                                                           o2::its::math_utils::Sq(positivePVParams.pvResolution) / float(vertex.getNContributors()));
  const float positivePVSigmaZ = o2::gpu::CAMath::Sqrt((o2::its::math_utils::Sq(positivePVResolution) * o2::its::math_utils::Sq(tanLambda) *
                                                        ((o2::its::math_utils::Sq(inverseR0) + sqInvDeltaZ0) * o2::its::math_utils::Sq(state.toRadius - state.fromRadius) + 1.f)) +
                                                       o2::its::math_utils::Sq((state.toRadius - state.fromRadius) * state.transitionMSAngle));
  BOOST_CHECK_EQUAL(positivePVWindow.sigmaZ, positivePVSigmaZ);
  BOOST_CHECK_GT(positivePVWindow.sigmaZ, window.sigmaZ);

  const float targetRadius = 4.f;
  const float targetZ = tanLambda * (targetRadius - source.radius) + source.zCoordinate;
  const auto acceptedTarget = makeGlobalCluster(targetRadius, 0.f, targetZ);
  const auto acceptedTargetMeasurement = makeMeasurement(acceptedTarget);
  float acceptedTanLambda = -999.f;
  BOOST_CHECK(window.acceptCandidate(sourceMeasurement, source, acceptedTargetMeasurement, acceptedTarget, acceptedTanLambda));
  BOOST_CHECK_EQUAL(acceptedTanLambda, (source.zCoordinate - acceptedTarget.zCoordinate) / (source.radius - acceptedTarget.radius));

  const auto rejectedTarget = makeGlobalCluster(-targetRadius, 0.f, targetZ);
  const auto rejectedTargetMeasurement = makeMeasurement(rejectedTarget);
  float rejectedTanLambda = 123.f;
  BOOST_CHECK(!window.acceptCandidate(sourceMeasurement, source, rejectedTargetMeasurement, rejectedTarget, rejectedTanLambda));
  BOOST_CHECK_EQUAL(rejectedTanLambda, 123.f);
}

BOOST_AUTO_TEST_CASE(DiskProjectSearchWindowReusesHelpersAndDirectProjectedXYBins)
{
  TrackingParameters legacy;
  const auto params = makeKernelParameters(legacy, SurfaceKind::Disk);
  BOOST_REQUIRE(params.isValid());

  IndexTableUtilsCore indexUtils;
  std::array<float, 10> halfExtents{};
  halfExtents.fill(20.f);
  indexUtils.setIndexTableParams(IndexTableCoordType::XY, legacy.RowBins, legacy.ColBins, -20.f, 20.f, halfExtents);

  constexpr int fromLayer = 1;
  constexpr int toLayer = 4; // deliberately skipped/nonadjacent transition
  const float fromZ = detail::mftLayerZ(fromLayer);
  const float toZ = detail::mftLayerZ(toLayer);
  const auto source = makeGlobalCluster(1.2f, 0.7f, fromZ);
  const auto sourceMeasurement = makeMeasurement(source, 2.e-4f, 3.e-4f);
  const auto vertex = makeVertex(0.01f, -0.02f, 0.1f, 4.e-4f, 5.e-4f, 0.04f, 3);
  const auto state = makeDiskProjectionCache(fromLayer, toLayer, 2.f, fromZ, toZ, 3.e-3f, 0.04f);

  DiskTrackletSearchWindow window{};
  BOOST_REQUIRE((projectDiskSearchWindow(
    sourceMeasurement, source, vertex, state, Bz, indexUtils, params, window)));

  float expectedX = 0.f;
  float expectedY = 0.f;
  detail::mftTrackletProject(source.xCoordinate, source.yCoordinate, source.zCoordinate,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             fromLayer, toLayer, Bz, params.trackletMinPt, expectedX, expectedY);
  float expectedSigmaX = 0.f;
  float expectedSigmaY = 0.f;
  detail::mftTrackletSigmaXY(source.xCoordinate, source.yCoordinate,
                             vertex.getX(), vertex.getY(), vertex.getZ(),
                             sourceMeasurement.covariance.xx, sourceMeasurement.covariance.yy,
                             vertex.getSigmaX2(), vertex.getSigmaY2(), vertex.getSigmaZ2(),
                             fromLayer, toLayer, state.fromRadius, state.toReferenceCoordinate - state.fromReferenceCoordinate,
                             state.transitionMSAngle, state.transitionPhiCut,
                             expectedX, expectedY, expectedSigmaX, expectedSigmaY);

  const float zSpread = params.nSigmaCut * vertex.getSigmaZ();
  const float zVtxMin = vertex.getZ() - zSpread;
  const float zVtxMax = vertex.getZ() + zSpread;
  const float absZFrom = std::abs(fromZ);
  const float absZTo = std::abs(toZ);
  const float denomMin = zVtxMax + absZFrom;
  const float denomMax = absZFrom + zVtxMin;
  float radialRangeMin = (std::abs(denomMin) > 1.e-6f) ? source.radius * (zVtxMax + absZTo) / denomMin : source.radius;
  float radialRangeMax = (std::abs(denomMax) > 1.e-6f) ? source.radius * (absZTo + zVtxMin) / denomMax : source.radius;
  if (radialRangeMin > radialRangeMax) {
    std::swap(radialRangeMin, radialRangeMax);
  }
  const auto directBins = getBinsRectClusterAtProj(expectedX, expectedY, toLayer,
                                                   radialRangeMin, radialRangeMax,
                                                   expectedSigmaX * params.nSigmaCut,
                                                   expectedSigmaY * params.nSigmaCut,
                                                   indexUtils);

  BOOST_CHECK_EQUAL(window.bins.x, directBins.x);
  BOOST_CHECK_EQUAL(window.bins.y, directBins.y);
  BOOST_CHECK_EQUAL(window.bins.z, directBins.z);
  BOOST_CHECK_EQUAL(window.bins.w, directBins.w);
  BOOST_CHECK_EQUAL(window.xProj, expectedX);
  BOOST_CHECK_EQUAL(window.yProj, expectedY);
  BOOST_CHECK_EQUAL(window.sigmaX, expectedSigmaX);
  BOOST_CHECK_EQUAL(window.sigmaY, expectedSigmaY);

  const auto acceptedTarget = makeGlobalCluster(expectedX, expectedY, toZ);
  const auto acceptedTargetMeasurement = makeMeasurement(acceptedTarget);
  float acceptedTanLambda = -999.f;
  BOOST_CHECK(window.acceptCandidate(sourceMeasurement, acceptedTargetMeasurement, acceptedTanLambda));
  BOOST_CHECK_EQUAL(acceptedTanLambda,
                    (source.zCoordinate - acceptedTarget.zCoordinate) /
                      (source.radius - acceptedTarget.radius));

  const auto sameRadiusTarget = makeGlobalCluster(source.radius, 0.f, toZ);
  float rejectedTanLambda = 123.f;
  BOOST_CHECK(!window.acceptCandidate(sourceMeasurement, makeMeasurement(sameRadiusTarget), rejectedTanLambda));
  BOOST_CHECK_EQUAL(rejectedTanLambda, 123.f);
}

BOOST_AUTO_TEST_CASE(ProjectSearchWindowInvalidBinsLeaveEveryOutputFieldUnchanged)
{
  TrackingParameters legacy;

  IndexTableUtilsCore cylinderIndexUtils;
  cylinderIndexUtils.setTrackingParameters(legacy);
  const auto cylinderParams = makeKernelParameters(legacy, SurfaceKind::Cylinder);
  const auto cylinderSource = makeGlobalCluster(2.f, 0.f, 100.f);
  const auto cylinderMeasurement = makeMeasurement(cylinderSource);
  const auto cylinderVertex = makeVertex(0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
  const auto cylinderState = makeCylinderProjectionCache(0, 3, 2.f, 4.f, 3.8f, 4.2f, 5.e-4f, 2.e-3f, 0.08f);
  const CylinderTrackletSearchWindow cylinderSentinel{
    {101, 102, 103, 104}, 105.f, 106.f, 107.f, 108.f};
  auto cylinderOut = cylinderSentinel;
  BOOST_CHECK(!(projectCylinderSearchWindow(
    cylinderMeasurement, cylinderSource, cylinderVertex, cylinderState, Bz, cylinderIndexUtils, cylinderParams, cylinderOut)));
  checkSearchWindowEqual(cylinderOut, cylinderSentinel);

  IndexTableUtilsCore diskIndexUtils;
  std::array<float, 10> tinyHalfExtents{};
  tinyHalfExtents.fill(0.01f);
  diskIndexUtils.setIndexTableParams(IndexTableCoordType::XY, legacy.RowBins, legacy.ColBins, -0.01f, 0.01f, tinyHalfExtents);
  const auto diskParams = makeKernelParameters(legacy, SurfaceKind::Disk);
  constexpr int fromLayer = 0;
  constexpr int toLayer = 1;
  const float fromZ = detail::mftLayerZ(fromLayer);
  const float toZ = detail::mftLayerZ(toLayer);
  const auto diskSource = makeGlobalCluster(1.f, 0.5f, fromZ);
  const auto diskMeasurement = makeMeasurement(diskSource);
  const auto diskVertex = makeVertex(0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
  const auto diskState = makeDiskProjectionCache(fromLayer, toLayer, 2.f, fromZ, toZ, 3.e-3f, 0.04f);
  const DiskTrackletSearchWindow diskSentinel{
    {201, 202, 203, 204}, 205.f, 206.f, 207.f, 208.f, 210.f};
  auto diskOut = diskSentinel;
  BOOST_CHECK(!(projectDiskSearchWindow(
    diskMeasurement, diskSource, diskVertex, diskState, Bz, diskIndexUtils, diskParams, diskOut)));
  checkSearchWindowEqual(diskOut, diskSentinel);
}

BOOST_AUTO_TEST_CASE(CylinderCandidateUsesPeriodicPhiAndStrictSigmaAndPhiBoundaries)
{
  CylinderTrackletSearchWindow window{
    {}, 0.f, 1.f, 0.02f, 5.f};

  const float wrapEpsilon = 0.005f;
  const auto wrappedSource = makeGlobalCluster(std::cos(wrapEpsilon), -std::sin(wrapEpsilon), 0.f);
  const auto wrappedTarget = makeGlobalCluster(2.f * std::cos(wrapEpsilon), 2.f * std::sin(wrapEpsilon), 0.f);
  float tanLambda = -9.f;
  BOOST_REQUIRE(o2::its::math_utils::isPhiDifferenceBelow(wrappedSource.phi, wrappedTarget.phi, window.phiCut));
  BOOST_CHECK(window.acceptCandidate(makeMeasurement(wrappedSource), wrappedSource, makeMeasurement(wrappedTarget), wrappedTarget, tanLambda));

  const auto source = makeGlobalCluster(1.f, 0.f, 0.f);
  const auto sourceMeasurement = makeMeasurement(source);
  const auto exactSigmaTarget = makeGlobalCluster(2.f, 0.f, 5.f);
  tanLambda = 71.f;
  BOOST_CHECK(!window.acceptCandidate(sourceMeasurement, source, makeMeasurement(exactSigmaTarget), exactSigmaTarget, tanLambda));
  BOOST_CHECK_EQUAL(tanLambda, 71.f);
  const auto insideSigmaTarget = makeGlobalCluster(2.f, 0.f, std::nextafter(5.f, 0.f));
  BOOST_CHECK(window.acceptCandidate(sourceMeasurement, source, makeMeasurement(insideSigmaTarget), insideSigmaTarget, tanLambda));

  const auto phiTarget = makeGlobalCluster(2.f * std::cos(0.125f), 2.f * std::sin(0.125f), 0.f);
  window.phiCut = o2::gpu::CAMath::Abs(source.phi - phiTarget.phi);
  const bool directPhiDecision = o2::its::math_utils::isPhiDifferenceBelow(source.phi, phiTarget.phi, window.phiCut);
  tanLambda = 72.f;
  BOOST_CHECK_EQUAL(window.acceptCandidate(sourceMeasurement, source, makeMeasurement(phiTarget), phiTarget, tanLambda), directPhiDecision);
  BOOST_CHECK(!directPhiDecision);
  BOOST_CHECK_EQUAL(tanLambda, 72.f);
}

BOOST_AUTO_TEST_CASE(DiskProjectionLeafRetainsNearZeroDenominatorGuardAndRadialSwap)
{
  TrackingParameters legacy;
  const auto params = makeKernelParameters(legacy, SurfaceKind::Disk);
  constexpr int fromLayer = 0;
  constexpr int toLayer = 1;
  const float fromZ = detail::mftLayerZ(fromLayer);
  const float toZ = detail::mftLayerZ(toLayer);
  const auto source = makeGlobalCluster(1.f, 0.5f, fromZ);
  const auto sourceMeasurement = makeMeasurement(source);
  const auto state = makeDiskProjectionCache(fromLayer, toLayer, 2.f, fromZ, toZ, 3.e-3f, 0.04f);

  IndexTableUtilsCore indexUtils;
  std::array<float, 10> halfExtents{};
  halfExtents.fill(200.f);
  indexUtils.setIndexTableParams(IndexTableCoordType::XY, legacy.RowBins, legacy.ColBins, -200.f, 200.f, halfExtents);

  const auto straightVertex = makeVertex(0.1f, -0.2f, 0.3f, 4.e-4f, 5.e-4f, 0.04f);
  DiskTrackletSearchWindow straightWindow{};
  BOOST_REQUIRE((projectDiskSearchWindow(
    sourceMeasurement, source, straightVertex, state, 0.f, indexUtils, params, straightWindow)));
  float expectedX = 0.f;
  float expectedY = 0.f;
  detail::mftTrackletProject(source.xCoordinate, source.yCoordinate, source.zCoordinate,
                             straightVertex.getX(), straightVertex.getY(), straightVertex.getZ(),
                             fromLayer, toLayer, 0.f, params.trackletMinPt, expectedX, expectedY);
  BOOST_CHECK_EQUAL(straightWindow.xProj, expectedX);
  BOOST_CHECK_EQUAL(straightWindow.yProj, expectedY);

  const auto fallbackVertex = makeVertex(0.1f, -0.2f, fromZ, 4.e-4f, 5.e-4f, 0.f);
  DiskTrackletSearchWindow fallbackWindow{};
  BOOST_REQUIRE((projectDiskSearchWindow(
    sourceMeasurement, source, fallbackVertex, state, 0.f, indexUtils, params, fallbackWindow)));
  expectedX = expectedY = 0.f;
  detail::mftTrackletProject(source.xCoordinate, source.yCoordinate, source.zCoordinate,
                             fallbackVertex.getX(), fallbackVertex.getY(), fallbackVertex.getZ(),
                             fromLayer, toLayer, 0.f, params.trackletMinPt, expectedX, expectedY);
  BOOST_CHECK_EQUAL(expectedX, source.xCoordinate);
  BOOST_CHECK_EQUAL(expectedY, source.yCoordinate);
  BOOST_CHECK_EQUAL(fallbackWindow.xProj, expectedX);
  BOOST_CHECK_EQUAL(fallbackWindow.yProj, expectedY);

  const auto swapVertex = makeVertex(0.f, 0.f, fromZ + 0.1f, 0.f, 0.f, 0.01f);
  const float zSpread = params.nSigmaCut * swapVertex.getSigmaZ();
  const float zVtxMin = swapVertex.getZ() - zSpread;
  const float zVtxMax = swapVertex.getZ() + zSpread;
  const float absZFrom = std::abs(fromZ);
  const float absZTo = std::abs(toZ);
  const float rawRadialMin = source.radius * (zVtxMax + absZTo) / (zVtxMax + absZFrom);
  const float rawRadialMax = source.radius * (absZTo + zVtxMin) / (absZFrom + zVtxMin);
  BOOST_REQUIRE_GT(rawRadialMin, rawRadialMax);
  DiskTrackletSearchWindow swapWindow{};
  BOOST_REQUIRE((projectDiskSearchWindow(
    sourceMeasurement, source, swapVertex, state, 0.f, indexUtils, params, swapWindow)));
  expectedX = expectedY = 0.f;
  detail::mftTrackletProject(source.xCoordinate, source.yCoordinate, source.zCoordinate,
                             swapVertex.getX(), swapVertex.getY(), swapVertex.getZ(),
                             fromLayer, toLayer, 0.f, params.trackletMinPt, expectedX, expectedY);
  float expectedSigmaX = 0.f;
  float expectedSigmaY = 0.f;
  detail::mftTrackletSigmaXY(source.xCoordinate, source.yCoordinate,
                             swapVertex.getX(), swapVertex.getY(), swapVertex.getZ(),
                             sourceMeasurement.covariance.xx, sourceMeasurement.covariance.yy,
                             swapVertex.getSigmaX2(), swapVertex.getSigmaY2(), swapVertex.getSigmaZ2(),
                             fromLayer, toLayer, state.fromRadius, state.toReferenceCoordinate - state.fromReferenceCoordinate,
                             state.transitionMSAngle, state.transitionPhiCut,
                             expectedX, expectedY, expectedSigmaX, expectedSigmaY);
  const auto directBins = getBinsRectClusterAtProj(expectedX, expectedY, toLayer,
                                                   rawRadialMax, rawRadialMin,
                                                   expectedSigmaX * params.nSigmaCut,
                                                   expectedSigmaY * params.nSigmaCut,
                                                   indexUtils);
  BOOST_CHECK_EQUAL(swapWindow.bins.x, directBins.x);
  BOOST_CHECK_EQUAL(swapWindow.bins.y, directBins.y);
  BOOST_CHECK_EQUAL(swapWindow.bins.z, directBins.z);
  BOOST_CHECK_EQUAL(swapWindow.bins.w, directBins.w);
}

BOOST_AUTO_TEST_CASE(DiskCandidatePreservesInverseVarianceAndStrictBoundarySemantics)
{
  const auto source = makeGlobalCluster(1.f, 0.f, -45.f);
  const auto distantTarget = makeGlobalCluster(100.f, -80.f, -47.f);
  const auto sourceMeasurement = makeMeasurement(source);
  const auto distantTargetMeasurement = makeMeasurement(distantTarget);
  DiskTrackletSearchWindow zeroSigmaWindow{
    {}, 0.f, 0.f, 0.f, -1.f, 5.f};
  float tanLambda = -8.f;
  BOOST_CHECK(zeroSigmaWindow.acceptCandidate(sourceMeasurement, distantTargetMeasurement, tanLambda));
  BOOST_CHECK_CLOSE(tanLambda,
                    (source.zCoordinate - distantTarget.zCoordinate) /
                      (source.radius - distantTarget.radius),
                    1.e-4f);

  DiskTrackletSearchWindow chi2Window{
    {}, 0.f, 0.f, 1.f, 1.f, 5.f};
  const auto exactChi2Target = makeGlobalCluster(5.f, 0.f, -47.f);
  tanLambda = 81.f;
  BOOST_CHECK(!chi2Window.acceptCandidate(sourceMeasurement, makeMeasurement(exactChi2Target), tanLambda));
  BOOST_CHECK_EQUAL(tanLambda, 81.f);
  const auto insideChi2Target = makeGlobalCluster(std::nextafter(5.f, 0.f), 0.f, -47.f);
  BOOST_CHECK(chi2Window.acceptCandidate(sourceMeasurement, makeMeasurement(insideChi2Target), tanLambda));

  const auto sameRadiusTarget = makeGlobalCluster(1.f, 0.f, -47.f);
  tanLambda = 82.f;
  BOOST_CHECK(!chi2Window.acceptCandidate(sourceMeasurement, makeMeasurement(sameRadiusTarget), tanLambda));
  BOOST_CHECK_EQUAL(tanLambda, 82.f);
  const auto distinctRadiusTarget = makeGlobalCluster(0.99f, 0.f, -47.f);
  BOOST_CHECK(chi2Window.acceptCandidate(sourceMeasurement, makeMeasurement(distinctRadiusTarget), tanLambda));
  BOOST_CHECK_CLOSE(tanLambda, 200.f, 1.e-3f);
}

BOOST_AUTO_TEST_CASE(NormalizedMeasurementsRemainAuthoritativeOverPoisonedLocatorCoordinates)
{
  TrackingParameters cylinderParameters;
  cylinderParameters.PVres = 0.f;
  const auto cylinderKernelParameters = makeKernelParameters(cylinderParameters, SurfaceKind::Cylinder);
  IndexTableUtilsCore cylinderIndex;
  cylinderIndex.setTrackingParameters(cylinderParameters);
  const auto vertex = makeVertex(0.f, 0.f, 0.f, 1.e-4f, 1.e-4f, 4.e-4f, 4);
  const auto cylinderState = makeCylinderProjectionCache(0, 1, 2.f, 4.f, 3.8f, 4.2f, 5.e-4f, 2.e-3f, 0.08f);
  const auto sourceMeasurement = makeMeasurement(2.f, 0.f, 0.5f);
  const auto targetMeasurement = makeMeasurement(4.f, 0.f, 1.f);
  const auto source = makeGlobalCluster(2.f, 0.f, 0.5f);
  const auto target = makeGlobalCluster(4.f, 0.f, 1.f);

  CylinderTrackletSearchWindow baseline{};
  BOOST_REQUIRE((projectCylinderSearchWindow(
    sourceMeasurement, source, vertex, cylinderState, Bz, cylinderIndex, cylinderKernelParameters, baseline)));
  float baselineTanLambda = -1.f;
  BOOST_REQUIRE(baseline.acceptCandidate(sourceMeasurement, source, targetMeasurement, target, baselineTanLambda));
  const float baselinePhi = o2::gpu::GPUCommonMath::ATan2(sourceMeasurement.position.y - targetMeasurement.position.y,
                                                          sourceMeasurement.position.x - targetMeasurement.position.x);

  auto poisonedSource = source;
  auto poisonedTarget = target;
  poisonedSource.xCoordinate = -999.f;
  poisonedSource.yCoordinate = 888.f;
  poisonedSource.zCoordinate = -777.f;
  poisonedTarget.xCoordinate = 666.f;
  poisonedTarget.yCoordinate = -555.f;
  poisonedTarget.zCoordinate = 444.f;
  CylinderTrackletSearchWindow poisonedWindow{};
  BOOST_REQUIRE((projectCylinderSearchWindow(
    sourceMeasurement, poisonedSource, vertex, cylinderState, Bz, cylinderIndex, cylinderKernelParameters, poisonedWindow)));
  checkSearchWindowEqual(poisonedWindow, baseline);
  float poisonedTanLambda = -2.f;
  BOOST_REQUIRE(poisonedWindow.acceptCandidate(sourceMeasurement, poisonedSource, targetMeasurement, poisonedTarget, poisonedTanLambda));
  BOOST_CHECK_EQUAL(poisonedTanLambda, baselineTanLambda);
  BOOST_CHECK_EQUAL(o2::gpu::GPUCommonMath::ATan2(sourceMeasurement.position.y - targetMeasurement.position.y,
                                                  sourceMeasurement.position.x - targetMeasurement.position.x),
                    baselinePhi);

  auto poisonedNavigationCache = source;
  poisonedNavigationCache.radius = 4.f;
  CylinderTrackletSearchWindow cachePoisonedWindow{};
  BOOST_REQUIRE((projectCylinderSearchWindow(
    sourceMeasurement, poisonedNavigationCache, vertex, cylinderState, Bz, cylinderIndex, cylinderKernelParameters, cachePoisonedWindow)));
  BOOST_CHECK_NE(cachePoisonedWindow.tanLambda, baseline.tanLambda);

  TrackingParameters diskParameters;
  const auto diskKernelParameters = makeKernelParameters(diskParameters, SurfaceKind::Disk);
  IndexTableUtilsCore diskIndex;
  std::array<float, 10> halfExtents{};
  halfExtents.fill(20.f);
  diskIndex.setIndexTableParams(IndexTableCoordType::XY, diskParameters.RowBins, diskParameters.ColBins, -20.f, 20.f, halfExtents);
  const float fromZ = detail::mftLayerZ(0);
  const float toZ = detail::mftLayerZ(1);
  const auto diskMeasurement = makeMeasurement(1.f, 0.5f, fromZ, 2.e-4f, 3.e-4f, 7.f);
  auto diskLocator = makeGlobalCluster(1.f, 0.5f, fromZ);
  const auto diskState = makeDiskProjectionCache(0, 1, 2.f, fromZ, toZ, 3.e-3f, 0.04f);
  DiskTrackletSearchWindow diskBaseline{};
  BOOST_REQUIRE((projectDiskSearchWindow(
    diskMeasurement, diskLocator, vertex, diskState, Bz, diskIndex, diskKernelParameters, diskBaseline)));
  diskLocator.xCoordinate = 123.f;
  diskLocator.yCoordinate = -321.f;
  diskLocator.zCoordinate = 456.f;
  auto uvPoisoned = diskMeasurement;
  uvPoisoned.covariance.xy = -12345.f;
  DiskTrackletSearchWindow diskPoisoned{};
  BOOST_REQUIRE((projectDiskSearchWindow(
    uvPoisoned, diskLocator, vertex, diskState, Bz, diskIndex, diskKernelParameters, diskPoisoned)));
  checkSearchWindowEqual(diskPoisoned, diskBaseline);
}

/// Gate 3 transition-preparation slice coverage (relocated from
/// TimeFrame::initialise() into TrackerTraits::initialiseTimeFrame(); see
/// CandidateFinding.h family scattering leaves,
/// family curvature leaves, prepareTransitionScatteringAndBending, and
/// host LayerGeometryConfigView). These tests verify
/// exact legacy-formula parity, the family-specific arithmetic literal that
/// integration review required preserved (not canonicalized), and the
/// order-sensitive oneOverR ratchet -- independently of TrackerTraits'
/// production traversal (covered separately in
/// testComputeLayerTrackletsOrchestration.cxx).

namespace
{
/// Independent re-transcription of the shared arithmetic in the frozen
/// ITS-only TimeFrame::initialise() (ITS/tracking/src/TimeFrame.cxx:352-370),
/// which the (now-removed) common-CA non-MFT branch reproduced verbatim.
/// Deliberately re-derived here rather than calling
/// prepareTransitionScatteringAndBending, so a transcription mistake in
/// either the operation or this reference would show up as a mismatch.
TransitionScatteringBendingPrep referenceTransitionScatteringAndBending(
  gsl::span<const float> perLayerMSAngle, int fromLayer, int toLayer,
  float r1, float r2, float clampedOneOverR, float res1, float res2)
{
  float ms2 = 0.f;
  for (int layer = fromLayer; layer < toLayer; ++layer) {
    ms2 += o2::its::math_utils::Sq(perLayerMSAngle[layer]);
  }
  const float msAngle = o2::gpu::CAMath::Sqrt(ms2);
  const float cosTheta1half = o2::gpu::CAMath::Sqrt(1.f - o2::its::math_utils::Sq(0.5f * r1 * clampedOneOverR));
  const float cosTheta2half = o2::gpu::CAMath::Sqrt(1.f - o2::its::math_utils::Sq(0.5f * r2 * clampedOneOverR));
  const float x = (r2 * cosTheta1half) - (r1 * cosTheta2half);
  const float delta = o2::gpu::CAMath::Sqrt(1.f / (1.f - 0.25f * o2::its::math_utils::Sq(x * clampedOneOverR)) *
                                            (o2::its::math_utils::Sq((0.25f * r1 * r2 * o2::its::math_utils::Sq(clampedOneOverR) / cosTheta2half) + cosTheta1half) * o2::its::math_utils::Sq(res1) +
                                             o2::its::math_utils::Sq((0.25f * r1 * r2 * o2::its::math_utils::Sq(clampedOneOverR) / cosTheta1half) + cosTheta2half) * o2::its::math_utils::Sq(res2)));
  const float phiCut = o2::gpu::CAMath::Min(o2::gpu::CAMath::ASin(0.5f * x * clampedOneOverR) + 2.f * msAngle + delta, o2::constants::math::PI * 0.5f);
  return TransitionScatteringBendingPrep{msAngle, phiCut};
}
} // namespace

BOOST_AUTO_TEST_CASE(CylinderScatteringAngleMatchesFrozenITSFormula)
{
  // Bit-exact vs the frozen ITS expression (ITS/tracking/src/TimeFrame.cxx:347):
  // math_utils::MSangle(0.14f, trkParam.TrackletMinPt, trkParam.LayerxX0[iLayer]).
  const std::array<float, 4> xX0Values{0.f, -0.001f, 5.e-3f, 1.e-2f};
  const std::array<float, 3> trackletMinPtValues{0.1f, 0.3f, 2.5f};
  for (float xX0 : xX0Values) {
    for (float trackletMinPt : trackletMinPtValues) {
      const float reference = o2::its::math_utils::MSangle(0.14f, trackletMinPt, xX0);
      const float actual = cylinderLayerMultipleScatteringAngle(
        CylinderLayerScatteringInputs{xX0}, trackletMinPt);
      BOOST_CHECK_EQUAL(actual, reference);
    }
  }
  // xX0 <= 0 behavior, explicit: legacy MSangle maps this to zero, not a
  // rejection; the typed operation must not add validation beyond it.
  BOOST_CHECK_EQUAL(cylinderLayerMultipleScatteringAngle(
                      CylinderLayerScatteringInputs{0.f}, 0.3f),
                    0.f);
}

BOOST_AUTO_TEST_CASE(DiskScatteringAngleMatchesLegacyMftFormulaWithExplicitReferenceZ)
{
  // Bit-exact vs the legacy detail::mftLayerMSAngle(layer, params), except
  // the Disk operation receives referenceCoordinate/layerRadius
  // explicitly instead of calling mftLayerZ()/LayerZCoordinate() internally.
  // mftLayerZ() is used here only to construct the *expected* legacy value,
  // exactly as this operation's caller (TrackerTraits::initialiseTimeFrame(),
  // via bindLegacyMFTReferenceCoordinates()) is required to do.
  TrackingParameters legacy;
  resetDetectorDefaults(legacy, o2::detectors::DetID::MFT);
  for (int layer : {0, 3, o2::mft::constants::mft::LayersNumber - 1}) {
    const float referenceZ = detail::mftLayerZ(layer);
    const float radius = legacy.LayerRadii[layer];
    const float xX0 = legacy.LayerxX0[layer];

    const float reference = detail::mftLayerMSAngle(layer, legacy);
    const float actual = diskLayerMultipleScatteringAngle(
      DiskLayerScatteringInputs{xX0, radius, referenceZ}, legacy.TrackletMinPt);
    BOOST_CHECK_EQUAL(actual, reference);
  }

  // xX0 == 0 behavior, explicit: the legacy formula has no special case for
  // it (sqrt(0 * cscLambda) == 0), and this operation must not add one.
  const float referenceZ = detail::mftLayerZ(0);
  const float zeroX0Actual = diskLayerMultipleScatteringAngle(
    DiskLayerScatteringInputs{0.f, legacy.LayerRadii[0], referenceZ}, legacy.TrackletMinPt);
  BOOST_CHECK_EQUAL(zeroX0Actual, 0.f);
}

BOOST_AUTO_TEST_CASE(DiskScatteringAngleNearZeroReferenceRadiusFallback)
{
  // Legacy fallback: |rRef| <= 1e-6 => tanlRef = 0 (detail::mftLayerMSAngle),
  // rather than dividing by a near-zero radius.
  TrackingParameters legacy;
  resetDetectorDefaults(legacy, o2::detectors::DetID::MFT);
  legacy.LayerRadii[0] = 1.e-9f; // below the legacy 1e-6 fallback threshold
  const float referenceZ = detail::mftLayerZ(0);

  const float reference = detail::mftLayerMSAngle(0, legacy);
  const float actual = diskLayerMultipleScatteringAngle(
    DiskLayerScatteringInputs{legacy.LayerxX0[0], legacy.LayerRadii[0], referenceZ},
    legacy.TrackletMinPt);
  BOOST_CHECK_EQUAL(actual, reference);

  // Cross-check the fallback actually engages: tanlRef == 0 (rRef below the
  // 1e-6 threshold) makes absTanl == 0, which is *not* > 1e-6 either, so
  // cscLambda takes the near-parallel-incidence sentinel 1e6f, not 1 -- i.e.
  // this input is genuinely exercising the near-zero-radius branch, not
  // merely reproducing an unrelated formula.
  const float expectedWithSentinelCscLambda = 0.0136f * (1.f / legacy.TrackletMinPt) * std::sqrt(legacy.LayerxX0[0] * 1.e6f);
  BOOST_CHECK_EQUAL(reference, expectedWithSentinelCscLambda);
}

BOOST_AUTO_TEST_CASE(LayerGeometryConfigViewChecksSpanSizeOnlyNotNumericValues)
{
  TrackingParameters legacy;
  legacy.LayerRadii = {0.f, -1.f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f}; // degenerate/negative radii
  legacy.LayerxX0 = {5.e-3f, 5.e-3f, 5.e-3f, 1.e-2f, 1.e-2f, 1.e-2f, 1.e-2f};       // valid
  const auto legacyMaterial = toMaterial(legacy.LayerxX0);

  const auto attachHitConfig = bindAttachHitConfig(gsl::span<const NominalSurfaceMaterial>(legacyMaterial), legacy);
  BOOST_CHECK(attachHitConfig.isValid(7));

  const auto geometryConfig = bindLayerGeometryConfig(legacy, attachHitConfig);
  // Legacy TimeFrame::initialise() never rejected degenerate/zero/negative
  // radii; this slice must not silently start doing so.
  BOOST_CHECK(geometryConfig.isValid(7));
  BOOST_CHECK(!geometryConfig.isValid(8)); // span-size check still applies

  // Negative xOverX0 must be rejected -- by AttachHitConfigView, the
  // single established contract for that data. LayerGeometryConfigView
  // borrows the same (rejected) span rather than independently re-validating
  // it, so it must not be read as a numeric-validity signal on its own.
  auto corrupted = legacy;
  corrupted.LayerxX0[3] = -1.f;
  const auto corruptedMaterial = toMaterial(corrupted.LayerxX0);
  const auto corruptedAttachHitConfig = bindAttachHitConfig(gsl::span<const NominalSurfaceMaterial>(corruptedMaterial), corrupted);
  BOOST_CHECK(!corruptedAttachHitConfig.isValid(7));
  const auto corruptedGeometryConfig = bindLayerGeometryConfig(corrupted, corruptedAttachHitConfig);
  BOOST_CHECK(corruptedGeometryConfig.isValid(7)); // size-only: still reports valid
}

BOOST_AUTO_TEST_CASE(BindLayerGeometryConfigBorrowsAttachHitLayerMaterialSpan)
{
  TrackingParameters legacy;
  const auto legacyMaterial = toMaterial(legacy.LayerxX0);
  const auto attachHitConfig = bindAttachHitConfig(gsl::span<const NominalSurfaceMaterial>(legacyMaterial), legacy);
  const auto geometryConfig = bindLayerGeometryConfig(legacy, attachHitConfig);
  BOOST_CHECK_EQUAL(geometryConfig.layerMaterial.data(), attachHitConfig.layerMaterial.data());
  BOOST_CHECK_EQUAL(geometryConfig.layerMaterial.size(), attachHitConfig.layerMaterial.size());
  BOOST_CHECK_EQUAL(geometryConfig.layerRadii.data(), legacy.LayerRadii.data());
}

BOOST_AUTO_TEST_CASE(ClampTransitionCurvatureUsesOneCoordinateNeutralExpression)
{
  const std::array<std::pair<float, float>, 5> samples{{
    {0.001f, 50.f}, // clamp does not trigger
    {3.0f, 1.0f},   // clamp triggers
    {0.02f, 25.f},
    {0.5f, 0.9f},
    {0.0001f, 4.f},
  }};
  for (const auto& sample : samples) {
    const float oneOverR = sample.first;
    const float r2 = sample.second;

    const float actual = clampTransitionCurvature(oneOverR, r2);
    const float reference = (0.5f * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
    BOOST_CHECK_EQUAL(actual, reference);
  }
}

BOOST_AUTO_TEST_CASE(CurvatureClampIsTransitionLocal)
{
  constexpr float initialOneOverR = 3.f;
  const std::array<float, 3> outerRadii{1.f, 4.f, 0.5f};
  for (const auto outerRadius : outerRadii) {
    const auto forward = clampTransitionCurvature(initialOneOverR, outerRadius);
    const auto repeated = clampTransitionCurvature(initialOneOverR, outerRadius);
    BOOST_CHECK_EQUAL(forward, repeated);
  }
}

BOOST_AUTO_TEST_CASE(PrepareTransitionScatteringAndBendingMatchesFrozenFormulaForITSAndMFTShapedInputs)
{
  // ITS-shaped (cm-scale barrel radii from TrackingParameters defaults).
  {
    const std::array<float, 7> msAngles{1.e-3f, 1.1e-3f, 1.2e-3f, 2.e-3f, 2.1e-3f, 2.2e-3f, 2.3e-3f};
    constexpr int fromLayer = 0;
    constexpr int toLayer = 3; // half-open: sums layers 0,1,2 only
    constexpr float r1 = 2.33959f;
    constexpr float r2 = 19.6213f;
    constexpr float res1 = 5.e-4f;
    constexpr float res2 = 5.e-4f;
    const float oneOverR = clampTransitionCurvature(
      0.001f * 0.3f * std::abs(Bz) / 0.3f, r2);
    const gsl::span<const float> msSpan(msAngles.data(), msAngles.size());
    const auto actual = prepareTransitionScatteringAndBending(msSpan, fromLayer, toLayer, r1, r2, oneOverR, res1, res2);
    const auto reference = referenceTransitionScatteringAndBending(msSpan, fromLayer, toLayer, r1, r2, oneOverR, res1, res2);
    BOOST_CHECK_EQUAL(actual.msAngle, reference.msAngle);
    BOOST_CHECK_EQUAL(actual.phiCut, reference.phiCut);

    // Half-open range: layer index `toLayer` itself must not contribute.
    const auto includingToLayer = referenceTransitionScatteringAndBending(msSpan, fromLayer, toLayer + 1, r1, r2, oneOverR, res1, res2);
    BOOST_CHECK_NE(actual.msAngle, includingToLayer.msAngle);
  }

  // MFT-shaped, deliberately skipped/non-adjacent transition (fromLayer=1,
  // toLayer=4: sums layers 1,2,3, skipping layer 4 itself as the endpoint).
  {
    TrackingParameters mft;
    resetDetectorDefaults(mft, o2::detectors::DetID::MFT);
    std::array<float, o2::mft::constants::mft::LayersNumber> msAngles{};
    for (int layer = 0; layer < o2::mft::constants::mft::LayersNumber; ++layer) {
      msAngles[layer] = diskLayerMultipleScatteringAngle(
        DiskLayerScatteringInputs{mft.LayerxX0[layer], mft.LayerRadii[layer], detail::mftLayerZ(layer)},
        mft.TrackletMinPt);
    }
    constexpr int fromLayer = 1;
    constexpr int toLayer = 4;
    const float r1 = mft.LayerRadii[fromLayer];
    const float r2 = mft.LayerRadii[toLayer];
    constexpr float res1 = 5.e-4f;
    constexpr float res2 = 6.e-4f;
    const float oneOverR = clampTransitionCurvature(
      0.001f * 0.3f * std::abs(Bz) / mft.TrackletMinPt, r2);
    const gsl::span<const float> msSpan(msAngles.data(), msAngles.size());
    const auto actual = prepareTransitionScatteringAndBending(msSpan, fromLayer, toLayer, r1, r2, oneOverR, res1, res2);
    const auto reference = referenceTransitionScatteringAndBending(msSpan, fromLayer, toLayer, r1, r2, oneOverR, res1, res2);
    BOOST_CHECK_EQUAL(actual.msAngle, reference.msAngle);
    BOOST_CHECK_EQUAL(actual.phiCut, reference.phiCut);
  }
}

BOOST_AUTO_TEST_CASE(PrepareTransitionScatteringAndBendingZeroFieldAndDegenerateRadiusMatchLegacyFormula)
{
  const std::array<float, 3> msAngles{1.e-3f, 1.2e-3f, 1.4e-3f};
  const gsl::span<const float> msSpan(msAngles.data(), msAngles.size());

  // Zero field: oneOverR's initial value (before any clamp) is exactly 0,
  // matching the legacy `0.001f * 0.3f * std::abs(mBz) / trkParam.TrackletMinPt`.
  {
    const float zeroFieldOneOverR = 0.001f * 0.3f * std::abs(0.f) / 0.3f;
    BOOST_CHECK_EQUAL(zeroFieldOneOverR, 0.f);
    const float clamped = clampTransitionCurvature(zeroFieldOneOverR, 5.f);
    BOOST_CHECK_EQUAL(clamped, 0.f); // 0.5*0 >= 1/5 is false: clamp does not trigger
    const auto actual = prepareTransitionScatteringAndBending(msSpan, 0, 2, 2.f, 5.f, clamped, 5.e-4f, 5.e-4f);
    const auto reference = referenceTransitionScatteringAndBending(msSpan, 0, 2, 2.f, 5.f, clamped, 5.e-4f, 5.e-4f);
    BOOST_CHECK_EQUAL(actual.msAngle, reference.msAngle);
    BOOST_CHECK_EQUAL(actual.phiCut, reference.phiCut);
  }

  // Degenerate radius (r2 == 0): legacy does not reject this -- it flows
  // through to whatever the floating-point expression produces. This test
  // asserts parity with that expression, not any particular finiteness.
  {
    const float oneOverR = clampTransitionCurvature(0.01f, 0.f);
    const auto actual = prepareTransitionScatteringAndBending(msSpan, 0, 2, 2.f, 0.f, oneOverR, 5.e-4f, 5.e-4f);
    const auto reference = referenceTransitionScatteringAndBending(msSpan, 0, 2, 2.f, 0.f, oneOverR, 5.e-4f, 5.e-4f);
    // BOOST_CHECK_EQUAL on NaN is always false (NaN != NaN); compare the bit
    // pattern so a NaN-vs-NaN legacy-parity match is still recognized as a pass.
    BOOST_CHECK(std::memcmp(&actual.msAngle, &reference.msAngle, sizeof(float)) == 0);
    BOOST_CHECK(std::memcmp(&actual.phiCut, &reference.phiCut, sizeof(float)) == 0);
  }
}
