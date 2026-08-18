// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT SurfaceMeasurementAdapters
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <vector>

#include <gsl/gsl>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/ClusterDecoding.h"

namespace
{
using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

constexpr ClusterShape Shape{7, 2, 4};
constexpr DecodedCluster ITSDecoded{
  {11.f, 12.f, 13.f},
  {21.f, 22.f, 23.f, 0.25f},
  {0.1f, 0.02f, 0.3f},
  Shape,
  3};
constexpr DecodedCluster MFTDecoded{
  {31.f, 32.f, 33.f},
  {},
  {0.4f, 0.f, 0.6f},
  Shape,
  9};
} // namespace

BOOST_AUTO_TEST_CASE(ITSNormalizedMeasurementPreservesDecodedFacts)
{
  const auto global = makeCylinderGlobalMeasurement(ITSDecoded, 12345);
  const auto normalized = makeCylinderSurfaceMeasurement(ITSDecoded);

  BOOST_CHECK_EQUAL(global.position.x, ITSDecoded.global.x);
  BOOST_CHECK_EQUAL(global.position.y, ITSDecoded.global.y);
  BOOST_CHECK_EQUAL(global.position.z, ITSDecoded.global.z);
  BOOST_CHECK_EQUAL(normalized.frame.q, ITSDecoded.cylinderFrame.q);
  BOOST_CHECK_EQUAL(normalized.frame.u, ITSDecoded.cylinderFrame.u);
  BOOST_CHECK_EQUAL(normalized.frame.v, ITSDecoded.cylinderFrame.v);
  BOOST_CHECK_EQUAL(normalized.frame.frameAngle, ITSDecoded.cylinderFrame.frameAngle);
  BOOST_CHECK_EQUAL(normalized.covariance.uu, ITSDecoded.rowColumnCovariance.uu);
  BOOST_CHECK_EQUAL(normalized.covariance.uv, ITSDecoded.rowColumnCovariance.uv);
  BOOST_CHECK_EQUAL(normalized.covariance.vv, ITSDecoded.rowColumnCovariance.vv);
  BOOST_CHECK_EQUAL(global.clusterId, 12345u);
}

BOOST_AUTO_TEST_CASE(MFTNormalizedDiskMeasurementUsesDescriptorAxes)
{
  const auto global = makeDiskGlobalMeasurement(MFTDecoded, 67890);
  const auto normalized = makeDiskSurfaceMeasurement(MFTDecoded);

  BOOST_CHECK_EQUAL(normalized.frame.q, 33.f);
  BOOST_CHECK_EQUAL(normalized.frame.u, 31.f);
  BOOST_CHECK_EQUAL(normalized.frame.v, 32.f);
  BOOST_CHECK_EQUAL(normalized.frame.frameAngle, 0.f);
  BOOST_CHECK_EQUAL(normalized.covariance.uu, 0.4f); // ALPIDE row -> global x
  BOOST_CHECK_EQUAL(normalized.covariance.uv, 0.f);
  BOOST_CHECK_EQUAL(normalized.covariance.vv, 0.6f); // ALPIDE column -> global y
  BOOST_CHECK_EQUAL(global.clusterId, 67890u);
}

BOOST_AUTO_TEST_CASE(SystematicErrorAdjustedCovarianceIsPreservedByBothProjections)
{
  constexpr DecodedCluster decodedAfterSystematicErrors{
    {1.f, 2.f, 3.f},
    {4.f, 5.f, 6.f, 0.5f},
    {0.1f + 0.7f, 0.f, 0.2f + 0.9f},
    {1, 1, 1},
    2};
  const auto itsMeasurement = makeCylinderSurfaceMeasurement(decodedAfterSystematicErrors);
  const auto mftMeasurement = makeDiskSurfaceMeasurement(decodedAfterSystematicErrors);

  BOOST_CHECK_EQUAL(itsMeasurement.covariance.uu, 0.8f);
  BOOST_CHECK_EQUAL(itsMeasurement.covariance.vv, 1.1f);
  BOOST_CHECK_EQUAL(mftMeasurement.covariance.uu, 0.8f);
  BOOST_CHECK_EQUAL(mftMeasurement.covariance.vv, 1.1f);
}

BOOST_AUTO_TEST_CASE(ExplicitPatternIsConsumedExactlyOnceAndProvidesShape)
{
  const std::array<unsigned char, 6> patterns{2, 3, 0xf0, 1, 2, 0xc0};
  const gsl::span<const unsigned char> patternSpan{patterns};
  auto pattIt = patternSpan.begin();
  const TopologyDictionary dictionary;
  const CompClusterExt first{10, 20, CompCluster::InvalidPatternID, 1};
  const CompClusterExt second{30, 40, CompCluster::InvalidPatternID, 2};
  float varianceRow{0.f};
  float varianceColumn{0.f};
  ClusterShape firstShape{};
  ClusterShape secondShape{};

  o2::itsmft::ioutils::extractClusterData(first, pattIt, &dictionary, varianceRow, varianceColumn, nullptr, &firstShape);
  BOOST_CHECK_EQUAL(std::distance(patternSpan.begin(), pattIt), 3);
  o2::itsmft::ioutils::extractClusterData(second, pattIt, &dictionary, varianceRow, varianceColumn, nullptr, &secondShape);
  BOOST_CHECK_EQUAL(std::distance(patternSpan.begin(), pattIt), 6);

  BOOST_CHECK_EQUAL(firstShape.nPixels, 4u);
  BOOST_CHECK_EQUAL(firstShape.rowSpan, 2u);
  BOOST_CHECK_EQUAL(firstShape.columnSpan, 3u);
  BOOST_CHECK_EQUAL(secondShape.nPixels, 2u);
  BOOST_CHECK_EQUAL(secondShape.rowSpan, 1u);
  BOOST_CHECK_EQUAL(secondShape.columnSpan, 2u);
}

BOOST_AUTO_TEST_CASE(BoundedExplicitPatternConsumptionRejectsEveryTruncationBoundary)
{
  // 3x3 requires two payload bytes: encoded size is 2 + ceil(9/8) = 4.
  constexpr std::array<unsigned char, 4> encoded{3, 3, 0x80, 0x80};
  const CompClusterExt cluster{10, 20, CompCluster::InvalidPatternID, 1};
  const TopologyDictionary dictionary;

  for (size_t available = 0; available < encoded.size(); ++available) {
    BoundedPatternCursor patterns{gsl::span<const unsigned char>{encoded.data(), available}};
    const auto decoded = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, &dictionary);
    BOOST_CHECK(decoded.error == ClusterDecodeError::TruncatedExplicitPattern);
    BOOST_CHECK_EQUAL(patterns.consumed(), 0u);
  }

  BoundedPatternCursor patterns{encoded};
  const auto decoded = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, &dictionary);
  BOOST_REQUIRE(decoded.ok());
  BOOST_CHECK_EQUAL(patterns.consumed(), encoded.size());
  BOOST_CHECK(patterns.empty());
  BOOST_CHECK_EQUAL(decoded.shape.nPixels, 2u);
  BOOST_CHECK_EQUAL(decoded.shape.rowSpan, 3u);
  BOOST_CHECK_EQUAL(decoded.shape.columnSpan, 3u);
}

BOOST_AUTO_TEST_CASE(BoundedExplicitPatternRejectsMalformedSpansAndEmptyBitmap)
{
  const TopologyDictionary dictionary;
  const CompClusterExt cluster{10, 20, CompCluster::InvalidPatternID, 1};
  const std::array<std::vector<unsigned char>, 4> malformed{
    std::vector<unsigned char>{0, 1},
    std::vector<unsigned char>{1, 0},
    std::vector<unsigned char>{129, 1},
    std::vector<unsigned char>{1, 1, 0x00}};

  for (const auto& encoded : malformed) {
    BoundedPatternCursor patterns{encoded};
    const auto decoded = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, &dictionary);
    BOOST_CHECK(decoded.error == ClusterDecodeError::MalformedExplicitPattern);
    BOOST_CHECK_EQUAL(patterns.consumed(), 0u);
  }
}

BOOST_AUTO_TEST_CASE(BoundedDecodeReportsMissingDictionaryBeforeReadingPatterns)
{
  constexpr std::array<unsigned char, 3> encoded{1, 1, 0x80};
  BoundedPatternCursor patterns{encoded};
  const CompClusterExt cluster{10, 20, CompCluster::InvalidPatternID, 1};
  const auto decoded = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, nullptr);
  BOOST_CHECK(decoded.error == ClusterDecodeError::MissingDictionary);
  BOOST_CHECK_EQUAL(patterns.consumed(), 0u);
}

BOOST_AUTO_TEST_CASE(SensorAndLayerValidationKeepsBothBounds)
{
  using o2::itsmft::ioutils::detail::isLayerInDetector;
  using o2::itsmft::ioutils::detail::isSensorInGeometry;

  BOOST_CHECK(!isSensorInGeometry(-1, 100));
  BOOST_CHECK(isSensorInGeometry(0, 100));
  BOOST_CHECK(isSensorInGeometry(99, 100));
  BOOST_CHECK(!isSensorInGeometry(100, 100));
  BOOST_CHECK(!isLayerInDetector(-1, 7));
  BOOST_CHECK(isLayerInDetector(0, 7));
  BOOST_CHECK(isLayerInDetector(6, 7));
  BOOST_CHECK(!isLayerInDetector(7, 7));
  BOOST_CHECK(!isLayerInDetector(10, 10));
}
