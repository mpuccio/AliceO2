// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT detector surface catalog aggregation
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <limits>

#include "DetectorSurfaceCatalogAggregation.h"
#include "ITSMFTBase/SegmentationAlpide.h"

using namespace o2::itsmft::tracking::detail;

namespace
{
GeometryPoint identity(size_t, const GeometryPoint& point)
{
  return point;
}

SurfaceGeometryAggregationResult aggregate(const LocalActiveArea& area,
                                           SurfaceReferenceCoordinate reference,
                                           const LocalToGlobal& transform = identity)
{
  return aggregateSurfaceGeometry(1, 1, area, reference, [](size_t) { return 0; }, transform);
}
} // namespace

BOOST_AUTO_TEST_CASE(asymmetric_active_center_offset)
{
  using Segmentation = o2::itsmft::SegmentationAlpide;
  const LocalActiveArea area{-0.5 * Segmentation::SensorSizeRows + Segmentation::PassiveEdgeReadOut,
                             0.5 * Segmentation::SensorSizeRows - Segmentation::PassiveEdgeTop,
                             -0.5 * Segmentation::ActiveMatrixSizeCols,
                             0.5 * Segmentation::ActiveMatrixSizeCols};
  const auto result = aggregate(area, SurfaceReferenceCoordinate::MeanRadius,
                                [](size_t, const GeometryPoint& point) {
                                  return GeometryPoint{10. + point.x, point.y, point.z};
                                });
  BOOST_REQUIRE(result.ok());
  const double expectedCenterX = 10. + 0.5 * (area.xMin + area.xMax);
  BOOST_CHECK_CLOSE(result.surfaces[0].referenceCoordinate, expectedCenterX, 1.e-5);
  BOOST_CHECK(std::abs(0.5 * (area.xMin + area.xMax)) > 1.e-3);
}

BOOST_AUTO_TEST_CASE(edge_midpoint_is_exact_radial_minimum)
{
  const auto result = aggregate({5., 15., -3., 3.}, SurfaceReferenceCoordinate::MeanRadius,
                                [](size_t, const GeometryPoint& point) {
                                  return GeometryPoint{point.x, point.z, 0.};
                                });
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_CLOSE(result.surfaces[0].radialMin, 5., 1.e-10);
  BOOST_CHECK_CLOSE(result.surfaces[0].radialMax, std::sqrt(234.), 1.e-5);
}

BOOST_AUTO_TEST_CASE(line_degenerate_its_like_projection)
{
  const auto result = aggregate({-2., 2., -3., 3.}, SurfaceReferenceCoordinate::MeanRadius,
                                [](size_t, const GeometryPoint& point) {
                                  return GeometryPoint{5., point.x, point.z};
                                });
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_CLOSE(result.surfaces[0].radialMin, 5., 1.e-10);
  BOOST_CHECK_CLOSE(result.surfaces[0].radialMax, std::sqrt(29.), 1.e-5);
}

BOOST_AUTO_TEST_CASE(fully_degenerate_point_projection)
{
  const auto result = aggregate({-2., 2., -3., 3.}, SurfaceReferenceCoordinate::MeanRadius,
                                [](size_t, const GeometryPoint&) {
                                  return GeometryPoint{3., 4., 8.};
                                });
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_CLOSE(result.surfaces[0].radialMin, 5., 1.e-10);
  BOOST_CHECK_CLOSE(result.surfaces[0].radialMax, 5., 1.e-10);
  BOOST_CHECK_CLOSE(result.surfaces[0].referenceCoordinate, 5., 1.e-10);
}

BOOST_AUTO_TEST_CASE(rotated_and_tilted_projection_containing_origin)
{
  constexpr double inverseSqrt2 = 0.7071067811865475244;
  const auto result = aggregate({-2., 2., -3., 3.}, SurfaceReferenceCoordinate::MeanZ,
                                [](size_t, const GeometryPoint& point) {
                                  return GeometryPoint{inverseSqrt2 * point.x + 0.5 * point.z,
                                                       inverseSqrt2 * point.x - 0.5 * point.z,
                                                       17. + inverseSqrt2 * point.z};
                                });
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_SMALL(result.surfaces[0].radialMin, 1.e-7f);
  BOOST_CHECK_CLOSE(result.surfaces[0].radialMax, std::sqrt(8.5), 1.e-5);
  BOOST_CHECK_CLOSE(result.surfaces[0].referenceCoordinate, 17., 1.e-10);
}

BOOST_AUTO_TEST_CASE(out_of_range_surface)
{
  const auto result = aggregateSurfaceGeometry(1, 1, {-1., 1., -1., 1.},
                                               SurfaceReferenceCoordinate::MeanRadius,
                                               [](size_t) { return 1; }, identity);
  BOOST_CHECK(result.error == SurfaceGeometryAggregationError::SurfaceOutOfRange);
  BOOST_CHECK(result.surfaces.empty());
}

BOOST_AUTO_TEST_CASE(empty_bucket)
{
  const auto result = aggregateSurfaceGeometry(1, 2, {-1., 1., -1., 1.},
                                               SurfaceReferenceCoordinate::MeanRadius,
                                               [](size_t) { return 0; }, identity);
  BOOST_CHECK(result.error == SurfaceGeometryAggregationError::EmptySurface);
  BOOST_CHECK(result.surfaces.empty());
}

BOOST_AUTO_TEST_CASE(non_finite_transform_coordinate)
{
  const auto result = aggregate({-1., 1., -1., 1.}, SurfaceReferenceCoordinate::MeanRadius,
                                [](size_t, const GeometryPoint& point) {
                                  return GeometryPoint{point.x, std::numeric_limits<double>::quiet_NaN(), point.z};
                                });
  BOOST_CHECK(result.error == SurfaceGeometryAggregationError::NonFiniteCoordinate);
  BOOST_CHECK(result.surfaces.empty());
}

BOOST_AUTO_TEST_CASE(invalid_derived_bounds)
{
  const double largest = std::numeric_limits<double>::max();
  const auto result = aggregate({-1., 1., -1., 1.}, SurfaceReferenceCoordinate::MeanRadius,
                                [largest](size_t, const GeometryPoint& point) {
                                  return GeometryPoint{largest, largest, point.z};
                                });
  BOOST_CHECK(result.error == SurfaceGeometryAggregationError::InvalidDerivedGeometry);
  BOOST_CHECK(result.surfaces.empty());
}
