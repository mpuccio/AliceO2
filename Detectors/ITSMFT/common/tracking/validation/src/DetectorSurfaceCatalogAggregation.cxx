// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "DetectorSurfaceCatalogAggregation.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace o2::itsmft::tracking::detail
{
namespace
{
struct ProjectedPoint {
  double x;
  double y;
};

bool finite(const GeometryPoint& point)
{
  return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z);
}

double cross(const ProjectedPoint& a, const ProjectedPoint& b)
{
  return a.x * b.y - a.y * b.x;
}

double pointToSegmentDistance(const ProjectedPoint& point,
                              const ProjectedPoint& first,
                              const ProjectedPoint& second,
                              double lengthTolerance)
{
  const ProjectedPoint edge{second.x - first.x, second.y - first.y};
  const double lengthSquared = edge.x * edge.x + edge.y * edge.y;
  if (lengthSquared <= lengthTolerance * lengthTolerance) {
    return std::hypot(point.x - first.x, point.y - first.y);
  }
  const double projection = ((point.x - first.x) * edge.x + (point.y - first.y) * edge.y) / lengthSquared;
  const double clamped = std::clamp(projection, 0., 1.);
  return std::hypot(point.x - (first.x + clamped * edge.x),
                    point.y - (first.y + clamped * edge.y));
}

double minimumRadius(const std::array<ProjectedPoint, 4>& polygon)
{
  double scale = 1.;
  for (size_t i = 0; i < polygon.size(); ++i) {
    const auto& point = polygon[i];
    const auto& next = polygon[(i + 1) % polygon.size()];
    scale = std::max({scale, std::abs(point.x), std::abs(point.y),
                      std::hypot(next.x - point.x, next.y - point.y)});
  }

  // Scale tolerances to the coordinates and quantities being tested.
  constexpr double toleranceFactor = 64. * std::numeric_limits<double>::epsilon();
  const double lengthTolerance = toleranceFactor * scale;
  const double areaTolerance = toleranceFactor * scale * scale;

  double twiceArea = 0.;
  for (size_t i = 0; i < polygon.size(); ++i) {
    twiceArea += cross(polygon[i], polygon[(i + 1) % polygon.size()]);
  }

  if (std::abs(twiceArea) > areaTolerance) {
    bool nonNegative = true;
    bool nonPositive = true;
    for (size_t i = 0; i < polygon.size(); ++i) {
      const auto& first = polygon[i];
      const auto& second = polygon[(i + 1) % polygon.size()];
      const ProjectedPoint edge{second.x - first.x, second.y - first.y};
      const ProjectedPoint toOrigin{-first.x, -first.y};
      const double side = cross(edge, toOrigin);
      nonNegative = nonNegative && side >= -areaTolerance;
      nonPositive = nonPositive && side <= areaTolerance;
    }
    if (nonNegative || nonPositive) {
      return 0.;
    }
  }

  const ProjectedPoint origin{0., 0.};
  double result = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < polygon.size(); ++i) {
    result = std::min(result, pointToSegmentDistance(origin, polygon[i],
                                                     polygon[(i + 1) % polygon.size()],
                                                     lengthTolerance));
  }
  return result;
}

bool canNarrowToFloat(double value)
{
  return std::isfinite(value) && std::abs(value) <= std::numeric_limits<float>::max();
}

struct Accumulator {
  size_t count{0};
  double referenceSum{0.};
  double radialMin{std::numeric_limits<double>::infinity()};
  double radialMax{0.};
};
} // namespace

SurfaceGeometryAggregationResult aggregateSurfaceGeometry(size_t chipCount,
                                                          size_t surfaceCount,
                                                          const LocalActiveArea& activeArea,
                                                          SurfaceReferenceCoordinate referenceCoordinate,
                                                          const SurfaceForChip& surfaceForChip,
                                                          const LocalToGlobal& localToGlobal)
{
  if (!std::isfinite(activeArea.xMin) || !std::isfinite(activeArea.xMax) ||
      !std::isfinite(activeArea.zMin) || !std::isfinite(activeArea.zMax) ||
      activeArea.xMin > activeArea.xMax || activeArea.zMin > activeArea.zMax) {
    return {{}, SurfaceGeometryAggregationError::InvalidDerivedGeometry};
  }

  const double centerX = 0.5 * (activeArea.xMin + activeArea.xMax);
  if (!std::isfinite(centerX)) {
    return {{}, SurfaceGeometryAggregationError::InvalidDerivedGeometry};
  }
  const std::array<GeometryPoint, 4> localCorners{{
    {activeArea.xMin, 0., activeArea.zMin},
    {activeArea.xMax, 0., activeArea.zMin},
    {activeArea.xMax, 0., activeArea.zMax},
    {activeArea.xMin, 0., activeArea.zMax}}};
  const GeometryPoint localCenter{centerX, 0., 0.};
  std::vector<Accumulator> accumulators(surfaceCount);

  for (size_t chip = 0; chip < chipCount; ++chip) {
    const int surface = surfaceForChip(chip);
    if (surface < 0 || static_cast<size_t>(surface) >= surfaceCount) {
      return {{}, SurfaceGeometryAggregationError::SurfaceOutOfRange};
    }

    std::array<ProjectedPoint, 4> projected{};
    double chipRadialMax = 0.;
    for (size_t corner = 0; corner < localCorners.size(); ++corner) {
      const auto global = localToGlobal(chip, localCorners[corner]);
      if (!finite(global)) {
        return {{}, SurfaceGeometryAggregationError::NonFiniteCoordinate};
      }
      projected[corner] = {global.x, global.y};
      const double radius = std::hypot(global.x, global.y);
      if (!std::isfinite(radius)) {
        return {{}, SurfaceGeometryAggregationError::InvalidDerivedGeometry};
      }
      chipRadialMax = std::max(chipRadialMax, radius);
    }

    const auto globalCenter = localToGlobal(chip, localCenter);
    if (!finite(globalCenter)) {
      return {{}, SurfaceGeometryAggregationError::NonFiniteCoordinate};
    }
    const double reference = referenceCoordinate == SurfaceReferenceCoordinate::MeanRadius
                               ? std::hypot(globalCenter.x, globalCenter.y)
                               : globalCenter.z;
    const double chipRadialMin = minimumRadius(projected);
    if (!std::isfinite(reference) || !std::isfinite(chipRadialMin) ||
        !std::isfinite(chipRadialMax) || chipRadialMin < 0. ||
        chipRadialMin > chipRadialMax) {
      return {{}, SurfaceGeometryAggregationError::InvalidDerivedGeometry};
    }

    auto& accumulator = accumulators[surface];
    ++accumulator.count;
    accumulator.referenceSum += reference;
    accumulator.radialMin = std::min(accumulator.radialMin, chipRadialMin);
    accumulator.radialMax = std::max(accumulator.radialMax, chipRadialMax);
    if (!std::isfinite(accumulator.referenceSum)) {
      return {{}, SurfaceGeometryAggregationError::InvalidDerivedGeometry};
    }
  }

  SurfaceGeometryAggregationResult result;
  result.surfaces.reserve(surfaceCount);
  for (const auto& accumulator : accumulators) {
    if (accumulator.count == 0) {
      return {{}, SurfaceGeometryAggregationError::EmptySurface};
    }
    const double reference = accumulator.referenceSum / static_cast<double>(accumulator.count);
    if (!canNarrowToFloat(reference) || !canNarrowToFloat(accumulator.radialMin) ||
        !canNarrowToFloat(accumulator.radialMax) || accumulator.radialMin > accumulator.radialMax) {
      return {{}, SurfaceGeometryAggregationError::InvalidDerivedGeometry};
    }
    const AggregatedSurfaceGeometry narrowed{static_cast<float>(reference),
                                             static_cast<float>(accumulator.radialMin),
                                             static_cast<float>(accumulator.radialMax)};
    if (!std::isfinite(narrowed.referenceCoordinate) || !std::isfinite(narrowed.radialMin) ||
        !std::isfinite(narrowed.radialMax) || narrowed.radialMin > narrowed.radialMax) {
      return {{}, SurfaceGeometryAggregationError::InvalidDerivedGeometry};
    }
    result.surfaces.push_back(narrowed);
  }
  return result;
}

} // namespace o2::itsmft::tracking::detail
