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
#include <cmath>
#include "ITSMFTTracking/ClusterDecoding.h"

namespace
{
using namespace o2::itsmft::tracking;

constexpr ClusterShape Shape{9, 3, 4};
} // namespace

BOOST_AUTO_TEST_CASE(CompactTypesHaveLockedABI)
{
  BOOST_CHECK_EQUAL(sizeof(GlobalMeasurement), 48u);
  BOOST_CHECK_EQUAL(offsetof(GlobalMeasurement, clusterId), 44u);
}

BOOST_AUTO_TEST_CASE(DefaultGlobalMeasurementHasInvalidClusterId)
{
  constexpr GlobalMeasurement measurement{};

  BOOST_CHECK(!measurement.hasValidClusterId());
}

BOOST_AUTO_TEST_CASE(ITSFixturePreservesCylinderConvention)
{
  const DecodedCluster decoded{
    {11.f, 12.f, 13.f},
    {21.f, 22.f, 23.f, 0.25f},
    {0.1f, 0.02f, 0.3f},
    Shape,
    3};
  const auto global = makeCylinderGlobalMeasurement(decoded, 7);
  const auto measurement = makeCylinderSurfaceMeasurement(decoded);

  BOOST_CHECK_EQUAL(global.position.x, 11.f);
  BOOST_CHECK_EQUAL(global.position.y, 12.f);
  BOOST_CHECK_EQUAL(global.position.z, 13.f);
  BOOST_CHECK_CLOSE(global.radius, std::hypot(11.f, 12.f), 1.e-5f);
  BOOST_CHECK_EQUAL(measurement.frame.q, 21.f);
  BOOST_CHECK_EQUAL(measurement.frame.u, 22.f);
  BOOST_CHECK_EQUAL(measurement.frame.v, 23.f);
  BOOST_CHECK_EQUAL(measurement.frame.frameAngle, 0.25f);
  BOOST_CHECK_EQUAL(measurement.covariance.uu, 0.1f);
  BOOST_CHECK_EQUAL(measurement.covariance.uv, 0.02f);
  BOOST_CHECK_EQUAL(measurement.covariance.vv, 0.3f);
  BOOST_CHECK_EQUAL(global.clusterId, 7u);
}

BOOST_AUTO_TEST_CASE(MFTFixtureUsesExplicitDiskCoordinatesAndXYCovariance)
{
  constexpr GlobalPoint3F global{31.f, 32.f, 33.f};
  constexpr SurfaceCovariance2F covarianceXY{0.4f, 0.05f, 0.6f};
  const DecodedCluster decoded{global, {}, covarianceXY, Shape, 9};
  const auto globalMeasurement = makeDiskGlobalMeasurement(decoded, 11);
  const auto measurement = makeDiskSurfaceMeasurement(decoded);

  BOOST_CHECK_EQUAL(measurement.frame.q, global.z);
  BOOST_CHECK_EQUAL(measurement.frame.u, global.x);
  BOOST_CHECK_EQUAL(measurement.frame.v, global.y);
  BOOST_CHECK_EQUAL(measurement.frame.frameAngle, 0.f);
  BOOST_CHECK_EQUAL(measurement.covariance.uu, covarianceXY.uu);
  BOOST_CHECK_EQUAL(measurement.covariance.uv, covarianceXY.uv);
  BOOST_CHECK_EQUAL(measurement.covariance.vv, covarianceXY.vv);
  BOOST_CHECK_EQUAL(globalMeasurement.clusterId, 11u);
}
