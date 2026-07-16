// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT SurfaceMeasurement
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <limits>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"

namespace
{
using namespace o2::itsmft::tracking;

constexpr ClusterShape Shape{9, 3, 4};
constexpr ClusterRef ITSCluster{ClusterSourceId{0}, 0};
constexpr DetectorSensorId ITSSensor{o2::detectors::DetID::ITS, 42};
} // namespace

BOOST_AUTO_TEST_CASE(ExternalIdentityIsSourceAndIndex)
{
  constexpr ClusterRef sourceZero{ClusterSourceId{0}, 0};
  constexpr ClusterRef sourceOne{ClusterSourceId{1}, 0};
  constexpr ClusterRef flaggedSourceZero{ClusterSourceId{0}, 0, 0xffff};

  BOOST_CHECK(sourceZero.isValid());
  BOOST_CHECK(sourceOne.isValid());
  BOOST_CHECK(sourceZero != sourceOne);
  BOOST_CHECK(sourceZero == flaggedSourceZero);
  BOOST_CHECK_EQUAL(sourceZero.flags, 0u);
  BOOST_CHECK_EQUAL(flaggedSourceZero.flags, 0xffffu);
  BOOST_CHECK_EQUAL(sourceZero.index, 0u);
  BOOST_CHECK_EQUAL(sourceOne.index, 0u);
}

BOOST_AUTO_TEST_CASE(FlagsHaveLockedABI)
{
  BOOST_CHECK_EQUAL(offsetof(ClusterRef, source), 0u);
  BOOST_CHECK_EQUAL(offsetof(ClusterRef, flags), 2u);
  BOOST_CHECK_EQUAL(offsetof(ClusterRef, index), 4u);
  BOOST_CHECK_EQUAL(offsetof(SurfaceMeasurement, surface), 68u);
  BOOST_CHECK_EQUAL(offsetof(SurfaceMeasurement, flags), 70u);
}

BOOST_AUTO_TEST_CASE(SensorIdentityIsDetectorQualified)
{
  constexpr DetectorSensorId its{o2::detectors::DetID::ITS, 7};
  constexpr DetectorSensorId mft{o2::detectors::DetID::MFT, 7};

  BOOST_CHECK(its.isValid());
  BOOST_CHECK(mft.isValid());
  BOOST_CHECK(its != mft);
  BOOST_CHECK_EQUAL(its.sensor, mft.sensor);
}

BOOST_AUTO_TEST_CASE(DefaultIdentifiersAreInvalid)
{
  constexpr ClusterSourceId source;
  constexpr ClusterRef cluster;
  constexpr DetectorSensorId sensor;
  constexpr SurfaceMeasurement measurement;

  BOOST_CHECK(!source.isValid());
  BOOST_CHECK_EQUAL(source.value(), std::numeric_limits<uint16_t>::max());
  BOOST_CHECK(!cluster.isValid());
  BOOST_CHECK(!sensor.isValid());
  BOOST_CHECK(!measurement.cluster.isValid());
  BOOST_CHECK(!measurement.sensor.isValid());
  BOOST_CHECK(!measurement.surface.isValid());
  BOOST_CHECK_EQUAL(cluster.flags, 0u);
  BOOST_CHECK_EQUAL(measurement.flags, 0u);
  BOOST_CHECK_EQUAL(measurement.sourceROF, std::numeric_limits<uint32_t>::max());
}

BOOST_AUTO_TEST_CASE(ITSFixturePreservesCylinderConventionAndMetadata)
{
  const o2::its::TrackingFrameInfo legacy{
    11.f, 12.f, 13.f,
    21.f, 0.25f,
    std::array<float, 2>{22.f, 23.f},
    std::array<float, 3>{0.1f, 0.02f, 0.3f}};
  const auto measurement = makeCylinderSurfaceMeasurement(legacy, ITSSensor, SurfaceId{3}, ITSCluster, 17, Shape);

  BOOST_CHECK_EQUAL(measurement.global.x, 11.f);
  BOOST_CHECK_EQUAL(measurement.global.y, 12.f);
  BOOST_CHECK_EQUAL(measurement.global.z, 13.f);
  BOOST_CHECK_EQUAL(measurement.frame.q, 21.f);
  BOOST_CHECK_EQUAL(measurement.frame.u, 22.f);
  BOOST_CHECK_EQUAL(measurement.frame.v, 23.f);
  BOOST_CHECK_EQUAL(measurement.frame.frameAngle, 0.25f);
  BOOST_CHECK_EQUAL(measurement.covariance.uu, 0.1f);
  BOOST_CHECK_EQUAL(measurement.covariance.uv, 0.02f);
  BOOST_CHECK_EQUAL(measurement.covariance.vv, 0.3f);
  BOOST_CHECK(measurement.sensor == ITSSensor);
  BOOST_CHECK(measurement.cluster == ITSCluster);
  BOOST_CHECK(measurement.surface == SurfaceId{3});
  BOOST_CHECK_EQUAL(measurement.flags, 0u);
  BOOST_CHECK_EQUAL(measurement.sourceROF, 17u);
  BOOST_CHECK_EQUAL(measurement.shape.nPixels, Shape.nPixels);
  BOOST_CHECK_EQUAL(measurement.shape.rowSpan, Shape.rowSpan);
  BOOST_CHECK_EQUAL(measurement.shape.columnSpan, Shape.columnSpan);
}

BOOST_AUTO_TEST_CASE(MFTFixtureUsesExplicitDiskCoordinatesAndXYCovariance)
{
  constexpr GlobalPoint3F global{31.f, 32.f, 33.f};
  constexpr SurfaceCovariance2F covarianceXY{0.4f, 0.05f, 0.6f};
  constexpr DetectorSensorId sensor{o2::detectors::DetID::MFT, 73};
  constexpr ClusterRef cluster{ClusterSourceId{1}, 0};
  const auto measurement = makeDiskSurfaceMeasurement(global, covarianceXY, sensor, SurfaceId{9}, cluster, 23, Shape);

  BOOST_CHECK_EQUAL(measurement.frame.q, global.z);
  BOOST_CHECK_EQUAL(measurement.frame.u, global.x);
  BOOST_CHECK_EQUAL(measurement.frame.v, global.y);
  BOOST_CHECK_EQUAL(measurement.frame.frameAngle, 0.f);
  BOOST_CHECK_EQUAL(measurement.covariance.uu, covarianceXY.uu);
  BOOST_CHECK_EQUAL(measurement.covariance.uv, covarianceXY.uv);
  BOOST_CHECK_EQUAL(measurement.covariance.vv, covarianceXY.vv);
  BOOST_CHECK(measurement.sensor == sensor);
  BOOST_CHECK(measurement.cluster == cluster);
  BOOST_CHECK(measurement.surface == SurfaceId{9});
  BOOST_CHECK_EQUAL(measurement.flags, 0u);
  BOOST_CHECK_EQUAL(measurement.sourceROF, 23u);
  BOOST_CHECK_EQUAL(measurement.shape.nPixels, Shape.nPixels);
  BOOST_CHECK_EQUAL(measurement.shape.rowSpan, Shape.rowSpan);
  BOOST_CHECK_EQUAL(measurement.shape.columnSpan, Shape.columnSpan);
}
