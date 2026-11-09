// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "GenericTracking/Base/LayerEnvelopes.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace o2::gtrk
{

namespace
{

/// The eight corners of a sensor's active silicon, in the global frame.
std::array<std::array<double, 3>, 8> corners(const Sensor& sensor, const SensorSpec& spec)
{
  std::array<std::array<double, 3>, 8> result{};
  std::size_t n = 0;
  for (const double u : {-spec.halfX, spec.halfX}) {
    for (const double v : {-spec.halfY, spec.halfY}) {
      for (const double w : {-spec.halfThickness, spec.halfThickness}) {
        const double local[3] = {u, v, w};
        sensor.toGlobal(local, result[n++].data());
      }
    }
  }
  return result;
}

/// Distance from the beam line to the nearest point of a convex body, given
/// its corners: the distance from the origin to the convex hull of their
/// transverse projections, zero if the hull contains it.
double nearestRadius(const std::array<std::array<double, 3>, 8>& points)
{
  std::array<std::array<double, 2>, 8> p{};
  for (std::size_t i = 0; i < p.size(); ++i) {
    p[i] = {points[i][0], points[i][1]};
  }
  std::ranges::sort(p);
  const auto cross = [](const auto& o, const auto& a, const auto& b) { return ((a[0] - o[0]) * (b[1] - o[1])) - ((a[1] - o[1]) * (b[0] - o[0])); };
  // Andrew's monotone chain, counter-clockwise.
  std::array<std::array<double, 2>, 16> hull{};
  std::size_t k = 0;
  for (const auto& point : p) {
    while (k >= 2 && cross(hull[k - 2], hull[k - 1], point) <= 0.) {
      --k;
    }
    hull[k++] = point;
  }
  for (std::size_t i = p.size() - 1, lower = k + 1; i-- > 0;) {
    while (k >= lower && cross(hull[k - 2], hull[k - 1], p[i]) <= 0.) {
      --k;
    }
    hull[k++] = p[i];
  }
  const std::size_t n = k - 1; // the last point repeats the first
  bool inside = n >= 3;
  double nearest = std::numeric_limits<double>::infinity();
  const std::array<double, 2> origin{0., 0.};
  for (std::size_t i = 0; i < std::max<std::size_t>(n, 1); ++i) {
    const auto& a = hull[i];
    const auto& b = hull[(i + 1) % std::max<std::size_t>(n, 1)];
    inside = inside && cross(a, b, origin) >= 0.;
    const double dx = b[0] - a[0];
    const double dy = b[1] - a[1];
    const double length2 = (dx * dx) + (dy * dy);
    const double t = length2 > 0. ? std::clamp(-((a[0] * dx) + (a[1] * dy)) / length2, 0., 1.) : 0.;
    nearest = std::min(nearest, std::hypot(a[0] + (t * dx), a[1] + (t * dy)));
  }
  return inside ? 0. : nearest;
}

/// The hull in r and z of a curved sensor's silicon, exactly: an annular
/// sector about the sensor's axis, which must be the beam's, so that z is the
/// same along the arc. In the transverse plane the silicon is the points
/// c + rho d, with c where the axis crosses it, rho from the radius to the
/// outer surface and d within the arc's half-angle of its middle direction.
/// Its distance from the beam line is extremal on the arc's two ends or where
/// d points along +-c; along a direction it is least at rho = -c.d, clamped.
SensorExtent curvedExtentOf(const Sensor& sensor, const SensorSpec& spec)
{
  // The sensor's local z axis is column 2 of the local -> global rotation.
  if (!exactlyEqual(sensor.rot[2], 0.) || !exactlyEqual(sensor.rot[5], 0.) || !exactlyEqual(std::abs(sensor.rot[8]), 1.)) {
    throw std::runtime_error(std::format("a curved sensor of layer {} whose axis is not the beam's", sensor.layer));
  }
  const double rhoMin = spec.radius + spec.surfaceMin();
  const double rhoMax = spec.radius + spec.surfaceMax();
  const double halfAngle = spec.halfX / spec.radius;
  const std::array<double, 2> c{sensor.trans[0], sensor.trans[1]};
  // A transverse direction of the sensor's frame, in the global one.
  const auto direction = [&](double angle) {
    const double local[3] = {std::cos(angle), std::sin(angle), 0.};
    return std::array<double, 2>{(sensor.rot[0] * local[0]) + (sensor.rot[1] * local[1]), (sensor.rot[3] * local[0]) + (sensor.rot[4] * local[1])};
  };
  const auto middle = direction(0.);
  std::vector<std::array<double, 2>> directions{direction(-halfAngle), direction(halfAngle)};
  const double cNorm = std::hypot(c[0], c[1]);
  if (cNorm > 0.) {
    for (const double sign : {1., -1.}) {
      const std::array<double, 2> d{sign * c[0] / cNorm, sign * c[1] / cNorm};
      if ((d[0] * middle[0]) + (d[1] * middle[1]) >= std::cos(halfAngle)) {
        directions.push_back(d);
      }
    }
  }
  SensorExtent extent;
  for (const auto& d : directions) {
    const auto radiusAt = [&](double rho) { return std::hypot(c[0] + (rho * d[0]), c[1] + (rho * d[1])); };
    extent.rMin = std::min(extent.rMin, radiusAt(std::clamp(-((c[0] * d[0]) + (c[1] * d[1])), rhoMin, rhoMax)));
    extent.rMax = std::max({extent.rMax, radiusAt(rhoMin), radiusAt(rhoMax)});
  }
  for (const double y : {-spec.halfY, spec.halfY}) {
    const double z = (sensor.rot[8] * y) + sensor.trans[2];
    extent.zMin = std::min(extent.zMin, z);
    extent.zMax = std::max(extent.zMax, z);
  }
  return extent;
}

/// Where a sensor's silicon is centred, in the global frame: a flat sensor's
/// origin, a curved sensor's middle of the silicon at the middle of its arc.
std::array<double, 3> centreOf(const Sensor& sensor, const SensorSpec& spec)
{
  double local[3];
  spec.unrolledToFrame(0., 0., 0.5 * (spec.surfaceMin() + spec.surfaceMax()), local);
  std::array<double, 3> global{};
  sensor.toGlobal(local, global.data());
  return global;
}

// The value in units of the resolution, floored, ceiled or rounded, as the
// nearest double to that decimal: the integer divided by the units per cm.
double down(double value, double resolution) { return std::floor(value / resolution) / std::round(1. / resolution); }
double up(double value, double resolution) { return std::ceil(value / resolution) / std::round(1. / resolution); }
double nearest(double value, double resolution) { return std::round(value / resolution) / std::round(1. / resolution); }

} // namespace

SensorExtent extentOf(const Sensor& sensor, const SensorSpec& spec)
{
  if (spec.curved()) {
    return curvedExtentOf(sensor, spec);
  }
  const auto points = corners(sensor, spec);
  SensorExtent extent;
  extent.rMin = nearestRadius(points);
  for (const auto& point : points) {
    extent.rMax = std::max(extent.rMax, std::hypot(point[0], point[1]));
    extent.zMin = std::min(extent.zMin, point[2]);
    extent.zMax = std::max(extent.zMax, point[2]);
  }
  return extent;
}

std::vector<LayerEnvelope> measureEnvelopes(std::span<const Sensor> sensors, std::span<const SensorSpec> specs, int nBarrelLayers, double resolution)
{
  const int nLayers = static_cast<int>(specs.size());
  std::vector<SensorExtent> extents(static_cast<std::size_t>(nLayers));
  std::vector<double> position(static_cast<std::size_t>(nLayers), 0.);
  std::vector<std::size_t> count(static_cast<std::size_t>(nLayers), 0U);
  for (const auto& sensor : sensors) {
    if (std::cmp_greater_equal(sensor.layer, nLayers)) {
      throw std::runtime_error(std::format("sensor on layer {}, beyond the {} layers", sensor.layer, nLayers));
    }
    const auto& spec = specs[sensor.layer];
    const auto extent = extentOf(sensor, spec);
    auto& layer = extents[sensor.layer];
    layer.rMin = std::min(layer.rMin, extent.rMin);
    layer.rMax = std::max(layer.rMax, extent.rMax);
    layer.zMin = std::min(layer.zMin, extent.zMin);
    layer.zMax = std::max(layer.zMax, extent.zMax);
    const auto centre = centreOf(sensor, spec);
    position[sensor.layer] += std::cmp_less(sensor.layer, nBarrelLayers) ? std::hypot(centre[0], centre[1]) : centre[2];
    ++count[sensor.layer];
  }
  std::vector<LayerEnvelope> result;
  for (std::size_t layer = 0; layer < extents.size(); ++layer) {
    if (count[layer] == 0U) {
      throw std::runtime_error(std::format("no sensor on layer {}", layer));
    }
    const auto& e = extents[layer];
    const double centre = position[layer] / static_cast<double>(count[layer]);
    if (resolution > 0.) {
      result.push_back({
        .position = nearest(centre, resolution),
        .rMin = down(e.rMin, resolution),
        .rMax = up(e.rMax, resolution),
        .zMin = down(e.zMin, resolution),
        .zMax = up(e.zMax, resolution),
      });
    } else {
      result.push_back({.position = centre, .rMin = e.rMin, .rMax = e.rMax, .zMin = e.zMin, .zMax = e.zMax});
    }
  }
  return result;
}

std::string describeEnvelopes(std::span<const LayerEnvelope> envelopes, int nBarrelLayers)
{
  std::string result;
  for (std::size_t layer = 0; layer < envelopes.size(); ++layer) {
    const auto& e = envelopes[layer];
    const bool barrel = std::cmp_less(layer, nBarrelLayers);
    result += std::format("layer {}: {} {:.4f}, r [{:.4f}, {:.4f}], z [{:.4f}, {:.4f}]\n", layer, barrel ? "mean radius" : "mean z", e.position, e.rMin, e.rMax, e.zMin, e.zMax);
  }
  return result;
}

} // namespace o2::gtrk
