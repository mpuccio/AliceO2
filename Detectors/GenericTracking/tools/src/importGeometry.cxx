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

/// \file importGeometry.cxx
/// \brief Measure the ITS, MFT or ALICE 3 TRK layer table from an o2sim
///        geometry, check the compiled one against it, and check our sensors
///        against O2's.
///
/// Loads the geometry, places every chip as GenericTracking/Alice/Placement.h does and
/// prints each layer's chip count as the layout's `barrelLayers` and
/// `discLayers` members, then the layer envelopes the tracking would
/// measure from these placements (Base/LayerEnvelopes.h), for information:
/// where the chips sit is run-time data and nothing compiled holds it.
/// Then it checks the compiled counts against the
/// placements, which is what fails after an O2 geometry change until the table
/// is pasted again.
///
/// For the TRK it also holds our sensors to O2's: every chip's four corner
/// pixels and its middle one, placed by O2's
/// segmentation and L2G matrix and by our sensor, must coincide. (The MC hits
/// on their chips: o2-gtrk-check-hits.)

#include "GenericTracking/Tools/Material.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerEnvelopes.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"

#include "DetectorsBase/GeometryManager.h"
#include "MathUtils/Cartesian.h"
#include "FT3Base/GeometryTGeo.h"
#include "TRKBase/GeometryTGeo.h"
#include "TRKBase/SegmentationChip.h"

#include <boost/program_options.hpp>

#include <algorithm>
#include <chrono>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <exception>
#include <format>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace bpo = boost::program_options;

namespace o2::gtrk::tools
{

namespace
{

/// How a layer's sensor is named in the layout header.
template <DetectorLayout Det>
std::string sensorName(std::size_t layer)
{
  if constexpr (std::same_as<Det, detectors::TRKFT3>) {
    return layer < detectors::trkft3::kVDRadius.size() ? std::format("trkft3::vdSensor({})", layer) : std::string{"trkft3::kChip"};
  } else {
    return "alpide::kSensor";
  }
}

/// Our sensors against O2's TRK and FT3: pixel centres placed by O2's
/// segmentation and L2G matrix -- through `flatToCurved` on the VD -- and by
/// our sensor. Returns the largest distance [cm].
double checkTRKPixels(std::span<const Sensor> sensors)
{
  using Segmentation = o2::trk::SegmentationChip;
  auto& geometry = *o2::trk::GeometryTGeo::Instance();
  auto& ft3 = *o2::ft3::GeometryTGeo::Instance();
  constexpr auto kFirstDisc = static_cast<std::size_t>(detectors::TRKFT3::kFirstDiscSensor);
  double worst = 0.;
  for (std::size_t chip = 0; chip < sensors.size(); ++chip) {
    const bool disc = chip >= kFirstDisc;
    const int id = static_cast<int>(disc ? chip - kFirstDisc : chip); // O2's chip ID
    const int subDet = disc ? ft3.getSubDetID(id) : geometry.getSubDetID(id);
    const int layer = disc ? ft3.getLayer(id) : geometry.getLayer(id); // within the sub-detector
    const auto& spec = Geometry<detectors::TRKFT3>::layerSpec(sensors[chip].layer);
    const int nRows = subDet == 0 ? o2::trk::constants::VD::petal::layer::nRows[static_cast<std::size_t>(layer)] : o2::trk::constants::moduleMLOT::chip::nRows;
    const int nCols = subDet == 0 ? o2::trk::constants::VD::petal::layer::nCols : o2::trk::constants::moduleMLOT::chip::nCols;
    for (const auto [row, col] : std::array<std::array<int, 2>, 5>{{{0, 0}, {nRows - 1, 0}, {0, nCols - 1}, {nRows - 1, nCols - 1}, {nRows / 2, nCols / 2}}}) {
      float xRow = 0.F;
      float zCol = 0.F;
      Segmentation::detectorToLocalUnchecked(row, col, xRow, zCol, subDet, layer, -1);
      o2::math_utils::Point3D<double> o2Local{static_cast<double>(xRow), 0., static_cast<double>(zCol)};
      if (subDet == 0) {
        const auto curved = Segmentation::flatToCurved(layer, xRow, 0.F);
        o2Local.SetCoordinates(static_cast<double>(curved.X()), static_cast<double>(curved.Y()), static_cast<double>(zCol));
      }
      const auto expected = disc ? ft3.getMatrixL2G(id)(o2Local) : geometry.getMatrixL2G(id)(o2Local);
      double frame[3];
      spec.unrolledToFrame(static_cast<double>(xRow), static_cast<double>(zCol), 0., frame);
      double ours[3];
      sensors[chip].toGlobal(frame, ours);
      worst = std::max(worst, std::hypot(ours[0] - expected.X(), ours[1] - expected.Y(), ours[2] - expected.Z()));
    }
  }
  return worst;
}

/// One passive scatterer per gap between consecutive barrel layers, and one
/// between the beam line and the first,
/// measured from the loaded geometry: its band the gap between the layers'
/// chips, their sensors' metal stacks included (`tools::withStacks`), its
/// extent along z the outer layer's, its x/X0 the
/// gap's material at normal incidence, averaged over azimuth and z, and its
/// radius that material's x/X0-weighted mean radius. Read along radial
/// segments at random azimuths and z (a fixed seed): some of the material is
/// thin radial walls and structures that end along z, which a regular grid
/// hits by aliasing. The x/X0 across the whole gap at kSamples points, the
/// radius in steps of kStep at kStepSamples. Printed as the `passives` table
/// of the layout.
template <DetectorLayout Det>
void printPassives(const Geometry<Det>& geometry)
{
  constexpr int kSamples = 65536;
  constexpr double kStep = 0.01; // [cm]
  constexpr int kStepSamples = 2048;
  const auto budget = [](double r0, double r1, double phi, double z) { return static_cast<double>(o2::base::GeometryManager::meanMaterialBudget(o2::math_utils::Point3D<double>{r0 * std::cos(phi), r0 * std::sin(phi), z}, o2::math_utils::Point3D<double>{r1 * std::cos(phi), r1 * std::sin(phi), z}).meanX2X0); };
  std::mt19937 rng(20261007);
  std::uniform_real_distribution<double> uniformPhi(0., 2. * std::numbers::pi);
  const auto envelopes = PlacedEnvelopes<Det>::all();
  // The TRK's chips carry a metal stack beyond the silicon; the others' layers end at it.
  std::vector<LayerEnvelope> chips(envelopes.begin(), envelopes.end());
  if constexpr (std::same_as<Det, detectors::TRKFT3>) {
    const auto specs = Geometry<Det>::layerSpecs();
    chips = tools::withStacks(envelopes, geometry.sensors(), specs, tools::trkStackThickness(geometry.sensors(), specs, Detector<Det>::nBarrelLayers));
  }
  std::cout << std::format(
                 "\n// Measured by o2-gtrk-import-geometry --passives: one below each barrel layer, the first from the beam line, its\n"
                 "// x/X0 across the gap at {} random azimuths and z, its radius the x/X0-weighted mean of {} cm steps at {}.\n",
                 kSamples, kStep, kStepSamples)
            << std::format("  static constexpr std::array<PassiveSpec, {}> passives{{{{\n", Detector<Det>::nBarrelLayers);
  for (int layer = 0; layer < Detector<Det>::nBarrelLayers; ++layer) {
    const auto& outer = envelopes[static_cast<std::size_t>(layer)];
    const double low = layer == 0 ? 0. : chips[static_cast<std::size_t>(layer - 1)].rMax;
    const double high = chips[static_cast<std::size_t>(layer)].rMin;
    const double zMax = std::max(std::abs(outer.zMin), std::abs(outer.zMax));
    std::uniform_real_distribution<double> uniformZ(-zMax, zMax);
    double total = 0.;
    for (int sample = 0; sample < kSamples; ++sample) {
      total += budget(low, high, uniformPhi(rng), uniformZ(rng));
    }
    double sum = 0.;
    double weighted = 0.;
    const auto steps = static_cast<int>(std::ceil((high - low) / kStep));
    for (int sample = 0; sample < kStepSamples; ++sample) {
      const double phi = uniformPhi(rng);
      const double z = uniformZ(rng);
      for (int step = 0; step < steps; ++step) {
        const double r = low + (step * kStep);
        const double r1 = std::min(r + kStep, high);
        const double x = budget(r, r1, phi, z);
        sum += x;
        weighted += x * 0.5 * (r + r1);
      }
    }
    const double position = sum > 0. ? weighted / sum : 0.5 * (low + high);
    std::cout << std::format(
      "    {{.kind = PassiveSpec::Kind::kBarrel, .position = {:.4f}, .bandLow = {:.4f}, .bandHigh = {:.4f}, .low = {:.2f}, .high = {:.2f}, "
      ".x2X0 = {:.3e}}},\n",
      position, low, high, -zMax, zMax, total / kSamples);
  }
  std::cout << "  }};\n";
}

/// The material forward of the barrel, where the barrel's passive cylinders
/// do not reach (the material the disc passes' seeds cross): between consecutive z
/// boundaries of the layout on one side -- the barrel layers' ends and the
/// discs' faces -- and in radial bands between the barrel layers' radii, the
/// x/X0 along z at random (r, phi) in the band, its mean and its x/X0-weighted
/// mean z. Printed as a table, from which disc-kind passives are drawn.
template <DetectorLayout Det>
void printForwardMaterial(int side)
{
  constexpr int kSamples = 1024;
  constexpr double kStep = 0.05; // [cm]
  const auto budget = [](double r, double phi, double z0, double z1) { return static_cast<double>(o2::base::GeometryManager::meanMaterialBudget(o2::math_utils::Point3D<double>{r * std::cos(phi), r * std::sin(phi), z0}, o2::math_utils::Point3D<double>{r * std::cos(phi), r * std::sin(phi), z1}).meanX2X0); };
  std::mt19937 rng(20261022);
  std::uniform_real_distribution<double> uniformPhi(0., 2. * std::numbers::pi);
  const auto envelopes = PlacedEnvelopes<Det>::all();
  std::vector<double> boundaries;
  std::vector<double> radii{0.};
  for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
    const auto& envelope = envelopes[static_cast<std::size_t>(layer)];
    if (layer < Detector<Det>::nBarrelLayers) {
      boundaries.push_back(std::max(std::abs(envelope.zMin), std::abs(envelope.zMax)));
      radii.push_back(envelope.rMin);
      radii.push_back(envelope.rMax);
    } else if ((side > 0 && envelope.zMin > 0.) || (side < 0 && envelope.zMax < 0.)) {
      boundaries.push_back(std::abs(envelope.zMin));
      boundaries.push_back(std::abs(envelope.zMax));
    }
  }
  std::ranges::sort(boundaries);
  std::ranges::sort(radii);
  boundaries.erase(std::ranges::unique(boundaries).begin(), boundaries.end());
  radii.erase(std::ranges::unique(radii).begin(), radii.end());
  std::cout << std::format(
    "\nforward material on the {} side, x/X0 along z across each gap between the layout's z boundaries, in radial bands [cm]: "
    "mean x/X0 [%], its x/X0-weighted mean |z| ({} rays a cell)\n",
    side > 0 ? "+z" : "-z", kSamples);
  std::string header = "  |z| gap        ";
  for (std::size_t band = 0; band + 1U < radii.size(); ++band) {
    header += std::format(" | r {:5.1f}-{:5.1f}", radii[band], radii[band + 1U]);
  }
  std::cout << header << '\n';
  for (std::size_t gap = 0; gap + 1U < boundaries.size(); ++gap) {
    const double z0 = boundaries[gap];
    const double z1 = boundaries[gap + 1U];
    if (z1 - z0 < 0.2) {
      continue; // a layer's own thickness
    }
    std::string line = std::format("  {:6.1f}-{:6.1f}", z0, z1);
    for (std::size_t band = 0; band + 1U < radii.size(); ++band) {
      std::uniform_real_distribution<double> uniformR(radii[band], radii[band + 1U]);
      double total = 0.;
      double sum = 0.;
      double weighted = 0.;
      const auto steps = static_cast<int>(std::ceil((z1 - z0) / kStep));
      for (int sample = 0; sample < kSamples; ++sample) {
        const double r = uniformR(rng);
        const double phi = uniformPhi(rng);
        total += budget(r, phi, side * z0, side * z1);
        if (sample < kSamples / 8) {
          for (int step = 0; step < steps; ++step) {
            const double a = z0 + (step * kStep);
            const double b = std::min(a + kStep, z1);
            const double x = budget(r, phi, side * a, side * b);
            sum += x;
            weighted += x * 0.5 * (a + b);
          }
        }
      }
      const double mean = 100. * total / kSamples;
      line += mean > 0.005 ? std::format(" | {:5.2f} @ {:6.1f}", mean, sum > 0. ? weighted / sum : 0.5 * (z0 + z1)) : " |      -       ";
    }
    std::cout << line << '\n';
  }
}

/// The azimuth a sensor covers at one value of its layer's second
/// coordinate, an arc: from `start` [0, 2 pi) counter-clockwise over `length`.
struct Arc {
  double start = 0.;
  double length = 0.;
};

constexpr double kTwoPi = 2. * std::numbers::pi;

double wrap(double angle)
{
  const double wrapped = std::fmod(angle, kTwoPi);
  return wrapped < 0. ? wrapped + kTwoPi : wrapped;
}

/// The shorter arc between the azimuths of two points.
Arc arcBetween(const std::array<double, 2>& a, const std::array<double, 2>& b)
{
  const double from = std::atan2(a[1], a[0]);
  double delta = std::atan2(b[1], b[0]) - from;
  delta = std::remainder(delta, kTwoPi);
  return delta >= 0. ? Arc{.start = wrap(from), .length = delta} : Arc{.start = wrap(from + delta), .length = -delta};
}

/// A flat sensor's active area in its mid-plane, its corners in order around it, global.
std::array<std::array<double, 3>, 4> midPlane(const Sensor& sensor, const SensorSpec& spec)
{
  const std::array<std::array<double, 2>, 4> local{{{-spec.halfX, -spec.halfY}, {spec.halfX, -spec.halfY}, {spec.halfX, spec.halfY}, {-spec.halfX, spec.halfY}}};
  std::array<std::array<double, 3>, 4> result{};
  for (std::size_t i = 0; i < local.size(); ++i) {
    const double point[3] = {local[i][0], local[i][1], 0.};
    sensor.toGlobal(point, result[i].data());
  }
  return result;
}

/// The angular gap a union of arcs may leave and still count as closed: the
/// rounding of the azimuths of abutting sensors.
constexpr double kAzimuthRounding = 1e-9;

/// Whether \a arcs cover every azimuth.
bool coversCircle(std::vector<Arc>& arcs)
{
  std::vector<std::array<double, 2>> pieces;
  for (const auto& arc : arcs) {
    if (arc.length >= kTwoPi) {
      return true;
    }
    const double end = arc.start + arc.length;
    if (end > kTwoPi) {
      pieces.push_back({arc.start, kTwoPi});
      pieces.push_back({0., end - kTwoPi});
    } else {
      pieces.push_back({arc.start, end});
    }
  }
  std::ranges::sort(pieces);
  double reach = 0.;
  for (const auto& [low, high] : pieces) {
    if (low > reach + kAzimuthRounding) {
      return false;
    }
    reach = std::max(reach, high);
  }
  return reach >= kTwoPi - kAzimuthRounding;
}

/// A disc sensor's arcs on the circle of radius \a r about the beam line,
/// from its mid-plane projected along the beam: a convex quadrilateral.
void discArcs(const std::array<std::array<double, 3>, 4>& quad, double r, std::vector<Arc>& arcs)
{
  const auto inside = [&](double x, double y) {
    bool positive = false;
    bool negative = false;
    for (std::size_t i = 0; i < quad.size(); ++i) {
      const auto& a = quad[i];
      const auto& b = quad[(i + 1) % quad.size()];
      const double cross = ((b[0] - a[0]) * (y - a[1])) - ((b[1] - a[1]) * (x - a[0]));
      positive = positive || cross > 0.;
      negative = negative || cross < 0.;
    }
    return !(positive && negative);
  };
  std::vector<double> angles;
  for (std::size_t i = 0; i < quad.size(); ++i) {
    const auto& p = quad[i];
    const auto& q = quad[(i + 1) % quad.size()];
    const double dx = q[0] - p[0];
    const double dy = q[1] - p[1];
    const double a = (dx * dx) + (dy * dy);
    const double b = 2. * ((p[0] * dx) + (p[1] * dy));
    const double c = (p[0] * p[0]) + (p[1] * p[1]) - (r * r);
    const double discriminant = (b * b) - (4. * a * c);
    if (a <= 0. || discriminant < 0.) {
      continue;
    }
    for (const double sign : {-1., 1.}) {
      const double t = (-b + (sign * std::sqrt(discriminant))) / (2. * a);
      if (t >= 0. && t <= 1.) {
        angles.push_back(wrap(std::atan2(p[1] + (t * dy), p[0] + (t * dx))));
      }
    }
  }
  if (angles.empty()) {
    if (inside(r, 0.)) {
      arcs.push_back({.start = 0., .length = kTwoPi});
    }
    return;
  }
  std::ranges::sort(angles);
  for (std::size_t i = 0; i < angles.size(); ++i) {
    const double from = angles[i];
    const double length = i + 1 < angles.size() ? angles[i + 1] - from : angles[0] + kTwoPi - from;
    const double middle = from + (0.5 * length);
    if (length > 0. && inside(r * std::cos(middle), r * std::sin(middle))) {
      arcs.push_back({.start = from, .length = length});
    }
  }
}

/// Where segments p0-p1 and q0-q1 of the transverse plane cross, if they do.
std::optional<std::array<double, 2>> crossing(const std::array<double, 3>& p0, const std::array<double, 3>& p1, const std::array<double, 3>& q0, const std::array<double, 3>& q1)
{
  const double rx = p1[0] - p0[0];
  const double ry = p1[1] - p0[1];
  const double sx = q1[0] - q0[0];
  const double sy = q1[1] - q0[1];
  const double denominator = (rx * sy) - (ry * sx);
  if (exactlyEqual(denominator, 0.)) {
    return std::nullopt;
  }
  const double t = (((q0[0] - p0[0]) * sy) - ((q0[1] - p0[1]) * sx)) / denominator;
  const double u = (((q0[0] - p0[0]) * ry) - ((q0[1] - p0[1]) * rx)) / denominator;
  if (t < 0. || t > 1. || u < 0. || u > 1.) {
    return std::nullopt;
  }
  return std::array<double, 2>{p0[0] + (t * rx), p0[1] + (t * ry)};
}

/// A barrel sensor's arc at \a z: a flat one's cut with the plane at z, seen
/// from the beam line; a curved one's whole arc when z is on it.
void barrelArcs(const Sensor& sensor, const SensorSpec& spec, double z, std::vector<Arc>& arcs)
{
  if (spec.curved()) {
    std::array<std::array<double, 2>, 2> ends{};
    for (std::size_t i = 0; i < 2; ++i) {
      double frame[3];
      spec.unrolledToFrame(i == 0 ? -spec.halfX : spec.halfX, 0., spec.halfThickness, frame);
      double global[3];
      sensor.toGlobal(frame, global);
      ends[i] = {global[0], global[1]};
    }
    arcs.push_back(arcBetween(ends[0], ends[1]));
    return;
  }
  const auto quad = midPlane(sensor, spec);
  std::vector<std::array<double, 2>> cut;
  for (std::size_t i = 0; i < quad.size(); ++i) {
    const auto& p = quad[i];
    const auto& q = quad[(i + 1) % quad.size()];
    if ((p[2] - z) * (q[2] - z) > 0. || exactlyEqual(p[2], q[2])) {
      continue;
    }
    const double t = (z - p[2]) / (q[2] - p[2]);
    cut.push_back({p[0] + (t * (q[0] - p[0])), p[1] + (t * (q[1] - p[1]))});
  }
  if (cut.size() >= 2U) {
    arcs.push_back(arcBetween(cut.front(), cut.back()));
  }
}

/// Where each layer covers every azimuth: the intervals of its
/// second coordinate -- z on a barrel layer, r on a disc -- at which some
/// sensor's active silicon lies at every azimuth, in increasing order. A
/// barrel sensor counts with its cut by the plane at that z, seen from the
/// beam line; a disc sensor with its mid-plane projected along the beam.
/// Layers numbered barrel first, as for `measureEnvelopes`.
std::vector<std::vector<std::array<double, 2>>> measureCoverage(std::span<const Sensor> sensors, std::span<const SensorSpec> specs, int nBarrelLayers)
{
  const auto nLayers = specs.size();
  std::vector<std::vector<std::size_t>> byLayer(nLayers);
  for (std::size_t i = 0; i < sensors.size(); ++i) {
    byLayer[sensors[i].layer].push_back(i);
  }
  std::vector<std::vector<std::array<double, 2>>> result(nLayers);
  for (std::size_t layer = 0; layer < nLayers; ++layer) {
    const bool barrel = std::cmp_less(layer, nBarrelLayers);
    const auto& spec = specs[layer];
    // Each sensor's range of the second coordinate, and the values where the
    // coverage can change: a barrel sensor's ends in z; a disc sensor's
    // corners, its edges' nearest points and the crossings of its edges with
    // its neighbours'.
    struct Span {
      double low;
      double high;
      std::size_t sensor;
    };
    std::vector<Span> spans;
    std::vector<double> critical;
    std::vector<std::array<std::array<double, 3>, 4>> quads;
    for (const auto index : byLayer[layer]) {
      const auto extent = extentOf(sensors[index], spec);
      spans.push_back({.low = barrel ? extent.zMin : extent.rMin, .high = barrel ? extent.zMax : extent.rMax, .sensor = quads.size()});
      if (barrel) {
        critical.push_back(extent.zMin);
        critical.push_back(extent.zMax);
        if (!spec.curved()) {
          for (const auto& corner : midPlane(sensors[index], spec)) {
            critical.push_back(corner[2]);
          }
        }
        quads.push_back({});
        continue;
      }
      const auto quad = midPlane(sensors[index], spec);
      quads.push_back(quad);
      for (std::size_t i = 0; i < quad.size(); ++i) {
        const auto& p = quad[i];
        const auto& q = quad[(i + 1) % quad.size()];
        critical.push_back(std::hypot(p[0], p[1]));
        const double dx = q[0] - p[0];
        const double dy = q[1] - p[1];
        const double length2 = (dx * dx) + (dy * dy);
        const double t = length2 > 0. ? -((p[0] * dx) + (p[1] * dy)) / length2 : 0.;
        if (t > 0. && t < 1.) {
          critical.push_back(std::hypot(p[0] + (t * dx), p[1] + (t * dy)));
        }
      }
    }
    if (!barrel) {
      const auto box = [&](std::size_t i) {
        std::array<double, 4> b{std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
        for (const auto& corner : quads[i]) {
          b[0] = std::min(b[0], corner[0]);
          b[1] = std::max(b[1], corner[0]);
          b[2] = std::min(b[2], corner[1]);
          b[3] = std::max(b[3], corner[1]);
        }
        return b;
      };
      std::vector<std::array<double, 4>> boxes(quads.size());
      for (std::size_t i = 0; i < quads.size(); ++i) {
        boxes[i] = box(i);
      }
      for (std::size_t i = 0; i < quads.size(); ++i) {
        for (std::size_t j = i + 1; j < quads.size(); ++j) {
          if (boxes[i][0] > boxes[j][1] || boxes[j][0] > boxes[i][1] || boxes[i][2] > boxes[j][3] || boxes[j][2] > boxes[i][3]) {
            continue;
          }
          for (std::size_t a = 0; a < 4; ++a) {
            for (std::size_t b = 0; b < 4; ++b) {
              if (const auto point = crossing(quads[i][a], quads[i][(a + 1) % 4], quads[j][b], quads[j][(b + 1) % 4])) {
                critical.push_back(std::hypot((*point)[0], (*point)[1]));
              }
            }
          }
        }
      }
    }
    std::ranges::sort(critical);
    const auto duplicates = std::ranges::unique(critical);
    critical.erase(duplicates.begin(), duplicates.end());
    std::ranges::sort(spans, {}, &Span::low);
    // Between two neighbouring critical values the coverage is the same:
    // tested at the middle, swept in order with the sensors whose range
    // holds it.
    std::vector<Span> active;
    std::size_t next = 0;
    std::vector<Arc> arcs;
    auto& covered = result[layer];
    for (std::size_t i = 0; i + 1 < critical.size(); ++i) {
      const double value = 0.5 * (critical[i] + critical[i + 1]);
      while (next < spans.size() && spans[next].low <= value) {
        active.push_back(spans[next++]);
      }
      std::erase_if(active, [&](const Span& span) { return span.high < value; });
      arcs.clear();
      for (const auto& span : active) {
        if (barrel) {
          barrelArcs(sensors[byLayer[layer][span.sensor]], spec, value, arcs);
        } else {
          discArcs(quads[span.sensor], value, arcs);
        }
      }
      if (!coversCircle(arcs)) {
        continue;
      }
      if (!covered.empty() && exactlyEqual(covered.back()[1], critical[i])) {
        covered.back()[1] = critical[i + 1];
      } else {
        covered.push_back({critical[i], critical[i + 1]});
      }
    }
  }
  return result;
}

template <DetectorLayout Det>
int runFor(const std::string& file, bool passives)
{
  const auto sensors = alice::placeChips<Det>();
  constexpr auto specs = Geometry<Det>::layerSpecs();
  const auto measured = measureEnvelopes(sensors, specs, Detector<Det>::nBarrelLayers, kEnvelopeResolution);
  std::vector<int> counts(measured.size(), 0);
  for (const auto& sensor : sensors) {
    ++counts[sensor.layer];
  }
  std::cout << std::format("  // Chip counts measured by o2-gtrk-import-geometry from {}: {} chips.\n", file, sensors.size());
  const auto printLayers = [&](const char* member, std::size_t first, std::size_t end) {
    if (first == end) {
      return;
    }
    std::cout << std::format("  static constexpr std::array<PlacedLayerSpec, {}> {}{{{{\n", end - first, member);
    for (std::size_t layer = first; layer < end; ++layer) {
      std::cout << std::format("    {{.sensorCount = {}, .sensor = {}}},\n", counts[layer], sensorName<Det>(layer));
    }
    std::cout << "  }};\n";
  };
  constexpr auto kBarrel = static_cast<std::size_t>(Detector<Det>::nBarrelLayers);
  printLayers("barrelLayers", 0, kBarrel);
  printLayers("discLayers", kBarrel, measured.size());
  std::cout << "\n"
            << std::format("Envelopes of the active silicon, rounded outwards to {} um, which the tracking measures at load and nothing compiled holds:\n", kEnvelopeResolution * 1e4) << describeEnvelopes(measured, Detector<Det>::nBarrelLayers);
  // Where each layer covers every azimuth: z on a barrel layer, r on a disc.
  const auto start = std::chrono::steady_clock::now();
  const auto coverage = measureCoverage(sensors, specs, Detector<Det>::nBarrelLayers);
  std::cout << std::format("\nevery azimuth covered, measured in {:.3f} s:\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
  for (std::size_t layer = 0; layer < coverage.size(); ++layer) {
    std::string line = std::format("layer {}:", layer);
    for (const auto& [low, high] : coverage[layer]) {
      line += std::format(" [{:.4f}, {:.4f}]", low, high);
    }
    std::cout << line << '\n';
  }
  int status = 0;
  for (std::size_t layer = 0; layer < measured.size(); ++layer) {
    const int stated = Detector<Det>::sensorCount(static_cast<LayerID>(layer));
    if (stated != counts[layer]) {
      std::cout << std::format("\ncompiled {} layer {}: {} chips, the geometry has {}", Det::name, layer, stated, counts[layer]);
      status = 1;
    }
  }
  std::cout << std::format("\ncompiled {} table {}\n", Det::name, status == 0 ? "has every chip" : "is stale: paste the table above");
  if (status != 0) {
    return status;
  }
  if constexpr (std::same_as<Det, detectors::TRKFT3>) {
    const double worst = checkTRKPixels(sensors);
    std::cout << std::format("\npixel centres, ours against O2's segmentation and L2G: at most {:.3e} um apart\n", worst * 1e4);
    // O2's segmentation works in float: a few nm at these radii.
    status = worst < 1e-5 ? 0 : 1;
  }
  if (passives) {
    const Geometry<Det> geometry(sensors);
    printPassives(geometry);
    if constexpr (Detector<Det>::hasDiscs) {
      printForwardMaterial<Det>(+1);
      printForwardMaterial<Det>(-1);
    }
  }
  return status;
}

int run(int argc, char** argv)
{
  std::string file;
  std::string detector;
  bool passives = false;
  bpo::options_description description("Measure the ITS, MFT or ALICE 3 TRK layer table from an o2sim geometry");
  description.add_options()                                                                                //
    ("help,h", "this message")                                                                             //
    ("geometry,g", bpo::value(&file)->default_value("o2sim_geometry-aligned.root"), "o2sim geometry file") //
    ("detector,d", bpo::value(&detector)->default_value("ITS"), "ITS, MFT or TRKFT3")                      //
    ("passives", bpo::bool_switch(&passives), "measure and print the barrel's passive scatterers, one per gap between its layers");
  bpo::variables_map vm;
  try {
    bpo::store(bpo::parse_command_line(argc, argv, description), vm);
    bpo::notify(vm);
  } catch (const bpo::error& error) {
    std::cerr << error.what() << "\n\n"
              << description << '\n';
    return 1;
  }
  if (vm.contains("help")) {
    std::cout << description << '\n';
    return 0;
  }
  o2::base::GeometryManager::loadGeometry(file);
  if (detector == "ITS") {
    return runFor<detectors::ITS>(file, passives);
  }
  if (detector == "MFT") {
    return runFor<detectors::MFT>(file, passives);
  }
  if (detector == "TRKFT3") {
    return runFor<detectors::TRKFT3>(file, passives);
  }
  std::cerr << "unknown detector '" << detector << "'\n";
  return 1;
}

} // namespace

} // namespace o2::gtrk::tools

int main(int argc, char** argv)
{
  try {
    return o2::gtrk::tools::run(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
