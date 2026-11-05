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

/// \file ToyBarrel8Fixture.h
/// \brief Timeframes of helices on the deep ToyBarrel8, for the tests of the
///        neighbour stage, the level sweep and the roads.

#ifndef GTRK_TEST_TOYBARREL8FIXTURE_H_
#define GTRK_TEST_TOYBARREL8FIXTURE_H_

#include "ToyBarrel8.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/Tracker.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Base/Geometry.h"

#include "GridTestSupport.h"
#include "UniformMaterial.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numbers>
#include <span>
#include <vector>

namespace o2::gtrk::toybarrel8
{

using Det = detectors::ToyBarrel8;
using Frame = TimeFrame<Det>;
using TrackletOutput = StageOutput<TrackletT>;
using Triplets = StageOutput<TripletT>;
using Controller = Tracker<Det, UniformMaterial>;
constexpr std::size_t kLayers = 8U;
/// The helices of the fixtures by index: direction, dip slope and pT.
constexpr std::size_t kHelixTracks = 5U;

struct InputPoint {
  LayerID layer = kInvalidLayerID;
  Scalar x = 0.F;
  Scalar y = 0.F;
  Scalar z = 0.F;
  Scalar r = 0.F;
  Scalar phi = 0.F;
};

inline SearchWindowParams windowParameters()
{
  SearchWindowParams result;
  result.minPtGeV = 1.;
  result.bzKiloGauss = 20.;
  result.materialX2X0 = 0.;
  result.positionSigma = 0.003;
  result.nSigma = 3.;
  return result;
}

inline GridParams gridParameters()
{
  GridParams result;
  result.blockWidth = 20'000;
  return test::withToyFrames(result);
}

inline IterationSettings iterationSettings() { return {.windows = windowParameters(), .vertexSource = VertexSource::kBeamLine}; }

/// Where a track from the origin with direction \a phi0, dip slope \a slope
/// and \a ptGeV crosses \a layer.
inline InputPoint helixAt(double phi0, double slope, double ptGeV, LayerID layer)
{
  const double radius = layerDistance<Det>(layer);
  const double curvature = 1. / curvatureRadius(ptGeV, 20.);
  const double arc = 2. * std::asin(0.5 * curvature * radius) / curvature;
  const double direction = phi0 + (curvature * arc);
  const double x = (std::sin(direction) - std::sin(phi0)) / curvature;
  const double y = (std::cos(phi0) - std::cos(direction)) / curvature;
  double phi = std::atan2(y, x);
  if (phi < 0.) {
    phi += 2 * std::numbers::pi_v<double>;
  }
  return {.layer = layer,
          .x = static_cast<Scalar>(x),
          .y = static_cast<Scalar>(y),
          .z = static_cast<Scalar>(slope * arc),
          .r = static_cast<Scalar>(radius),
          .phi = static_cast<Scalar>(phi)};
}

inline InputPoint helixPoint(std::size_t track, LayerID layer)
{
  constexpr std::array<double, kHelixTracks> kPhi{{0.10, 1.35, 3.05, 4.75, 5.80}};
  constexpr std::array<double, kHelixTracks> kSlope{{0.04, -0.05, 0.03, -0.02, 0.01}};
  constexpr std::array<double, kHelixTracks> kPt{{1.2, 1.8, 3.0, 1.4, 2.0}};
  return helixAt(kPhi[track], kSlope[track], kPt[track], layer);
}

/// One measurement of a fixture; its index in the list is its cluster index.
struct Hit {
  InputPoint point;
  double timeNS = 10.;
  /// The sensor and local position, for hits made by sensorHelixHits; the
  /// other fixtures leave them unset, as the stages before the fit read only
  /// the global position.
  SensorID surface = kInvalidSensorID;
  Scalar loc0 = 0.F;
  Scalar loc1 = 0.F;
};

/// The detector's sensors, built once.
inline const TrackingGeometry<Det>& trackingGeometry()
{
  static const auto detector = toyGeometry<Det>();
  static const TrackingGeometry<Det> tracking{detector};
  return tracking;
}

/// The hits of a helix from the origin on the sensors it crosses, one per
/// layer: on each layer the first sensor whose plane the helix meets inside
/// the sensor's bounds, near where it crosses the layer's radius. The global
/// position and the local one are the same point, and the time is the flight
/// time from the origin at \a startNS, at the speed of light.
inline std::vector<Hit> sensorHelixHits(double phi0, double slope, double ptGeV, double startNS = 10.)
{
  const auto& tracking = trackingGeometry();
  const double curvature = 1. / curvatureRadius(ptGeV, 20.);
  const auto at = [&](double arc, double position[3]) {
    const double direction = phi0 + (curvature * arc);
    position[0] = (std::sin(direction) - std::sin(phi0)) / curvature;
    position[1] = (std::cos(phi0) - std::cos(direction)) / curvature;
    position[2] = slope * arc;
  };
  std::vector<Hit> hits;
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    const double radius = layerDistance<Det>(id);
    const double crossing = 2. * std::asin(0.5 * curvature * radius) / curvature;
    const SensorID begin = firstSensorOfLayer<Det>(id);
    const SensorID end = firstSensorOfLayer<Det>(static_cast<LayerID>(layer + 1));
    for (SensorID sensor = begin; sensor < end; ++sensor) {
      const auto& surface = tracking.surface(sensor);
      const auto height = [&](double arc) {
        double position[3]{};
        at(arc, position);
        double result = 0.;
        for (int i = 0; i < 3; ++i) {
          result += (position[i] - static_cast<double>(surface.origin[i])) * static_cast<double>(surface.normal[i]);
        }
        return result;
      };
      double low = 0.8 * crossing;
      double high = 1.3 * crossing;
      if ((height(low) < 0.) == (height(high) < 0.)) {
        continue;
      }
      for (int i = 0; i < 100; ++i) {
        const double middle = 0.5 * (low + high);
        ((height(middle) < 0.) == (height(low) < 0.) ? low : high) = middle;
      }
      const double arc = 0.5 * (low + high);
      double position[3]{};
      at(arc, position);
      const Scalar global[3]{static_cast<Scalar>(position[0]), static_cast<Scalar>(position[1]), static_cast<Scalar>(position[2])};
      Scalar local[3]{};
      surfaceToLocal(surface, global, local);
      if (!surface.bounds.isInside(local[0], local[1])) {
        continue;
      }
      double phi = std::atan2(position[1], position[0]);
      if (phi < 0.) {
        phi += 2 * std::numbers::pi_v<double>;
      }
      hits.push_back({.point = {.layer = id,
                                .x = global[0],
                                .y = global[1],
                                .z = global[2],
                                .r = static_cast<Scalar>(std::hypot(position[0], position[1])),
                                .phi = static_cast<Scalar>(phi)},
                      .timeNS = startNS + (arc * std::hypot(1., slope) / static_cast<double>(constants::physics::LightSpeedCm2NS)),
                      .surface = sensor,
                      .loc0 = local[0],
                      .loc1 = local[1]});
      break;
    }
  }
  return hits;
}

/// \a point moved along z and r-phi [cm].
inline InputPoint moved(InputPoint point, double dz, double dRPhi)
{
  const double phi = static_cast<double>(point.phi) + (dRPhi / static_cast<double>(point.r));
  point.x = static_cast<Scalar>(static_cast<double>(point.r) * std::cos(phi));
  point.y = static_cast<Scalar>(static_cast<double>(point.r) * std::sin(phi));
  point.phi = static_cast<Scalar>(phi);
  point.z += static_cast<Scalar>(dz);
  return point;
}

class Fixture
{
 public:
  static constexpr VertexConstraint kBeamConstraint{.z = 0., .sigmaZ = 0.02, .sigmaTransverse = 0.001};

  explicit Fixture(const std::vector<Hit>& hits, double durationNS = 20.)
    : Fixture(hits, durationNS, gridParameters(), [](const Hit& hit) { return test::toyTimeBox(hit.timeNS); })
  {
  }

  /// With the grid's \a params and each hit's box from \a boxOf.
  template <typename BoxOf>
  Fixture(const std::vector<Hit>& hits, double durationNS, const GridParams& params, const BoxOf& boxOf)
    : geometry(test::cappedGrid<Det>(durationNS, params, windowParameters(), kBeamConstraint, 16, 8)),
      offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U)
  {
    frame = test::frameOf(geometry, hits.size(), [&](std::size_t i) { return hits[i].point; }, [&](std::size_t i) { return boxOf(hits[i]); }, [&](std::size_t i) { const auto scale = static_cast<Scalar>(1U + (i % 3U)); return test::PointVariances{scale * 1.e-6F, scale * 2.e-6F, scale * 3.e-6F}; }, columns, offsets, [&](std::size_t i, std::uint32_t row) {
      columns.surface[row] = hits[i].surface;
      columns.loc0[row] = hits[i].loc0;
      columns.loc1[row] = hits[i].loc1;
      columns.varLoc0[row] = columns.varRPhi[row];
      columns.varLoc1[row] = columns.varZ[row]; });
  }

  [[nodiscard]] static constexpr std::size_t logicalIndex(std::size_t track, std::size_t layer) { return (kLayers * track) + layer; }

  [[nodiscard]] std::uint32_t row(std::size_t logical) const
  {
    return static_cast<std::uint32_t>(std::ranges::find(columns.clusterIndex, static_cast<std::uint32_t>(logical)) - columns.clusterIndex.begin());
  }

  GridGeometry<Det> geometry;
  std::vector<std::uint32_t> offsets;
  MeasurementStore columns;
  std::unique_ptr<Frame> frame;
};

/// The triplet whose measurements are the rows \a a, \a b, \a c.
inline std::uint32_t findTriplet(const TrackletOutput& tracklets, const Triplets& triplets, std::uint32_t a, std::uint32_t b, std::uint32_t c)
{
  const auto trackletView = tracklets.view();
  const auto view = triplets.view();
  for (std::uint32_t index = 0; index < triplets.size(); ++index) {
    const auto first = view.first[index];
    const auto second = view.second[index];
    if (trackletView.source[first] == a && trackletView.target[first] == b && trackletView.target[second] == c) {
      return index;
    }
  }
  return triplets.size();
}

namespace levels
{

/// The level fixture's cluster indices: track 0 through all eight layers,
/// a second hit next to it on layer 3, track 1 without its layer-3 hit, and
/// track 2 with hits on layers 5-7 only.
constexpr std::uint32_t kSecondOnLayer3 = 8U;
constexpr std::array<std::uint32_t, kLayers> kTrack1{{9U, 10U, 11U, 0U, 12U, 13U, 14U, 15U}}; // [3] unused: no hit
constexpr std::array<std::uint32_t, 3> kTrack2{{16U, 17U, 18U}};

inline std::vector<Hit> levelHits()
{
  std::vector<Hit> hits;
  hits.reserve(kTrack2.back() + 1U);
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    hits.push_back({.point = helixPoint(0U, static_cast<LayerID>(layer))});
  }
  hits.push_back({.point = moved(helixPoint(0U, 3), 0.0005, 0.)});
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    if (layer != 3U) {
      hits.push_back({.point = helixPoint(1U, static_cast<LayerID>(layer))});
    }
  }
  for (std::size_t layer = 5; layer < kLayers; ++layer) {
    hits.push_back({.point = helixPoint(2U, static_cast<LayerID>(layer))});
  }
  return hits;
}

} // namespace levels

/// The slicing fixture: tracks through all eight layers, spread over the
/// timeframe. Every other track's timestamps step one clock period after
/// layers 2 and 5 -- each pair, triplet and `a-d` span stays within the time
/// cut, and the chain crosses several 20 ns blocks. Every third track has a
/// second hit on layer 3, so some triplets have two predecessors.
constexpr double kSliceDurationNS = 1100.;
constexpr std::size_t kSliceTracks = 12U;

inline std::vector<Hit> sliceHits()
{
  constexpr std::array<double, kLayers> kSteps{{0., 0., 0., 25., 25., 25., 50., 50.}};
  std::vector<Hit> hits;
  for (std::size_t track = 0; track < kSliceTracks; ++track) {
    const double phi0 = std::fmod(0.3 + (0.9 * static_cast<double>(track)), 2 * std::numbers::pi_v<double>);
    const double slope = (track % 2U == 0U ? 0.03 : -0.02) + (0.002 * static_cast<double>(track));
    const double pt = 1.2 + (0.3 * static_cast<double>(track % 4U));
    const double start = 37. + (83. * static_cast<double>(track));
    for (std::size_t layer = 0; layer < kLayers; ++layer) {
      const double step = track % 2U == 0U ? kSteps[layer] : 0.;
      hits.push_back({.point = helixAt(phi0, slope, pt, static_cast<LayerID>(layer)), .timeNS = start + step});
    }
    if (track % 3U == 0U) {
      const double step = track % 2U == 0U ? kSteps[3] : 0.;
      hits.push_back({.point = moved(helixAt(phi0, slope, pt, 3), 0.0005, 0.), .timeNS = start + step});
    }
  }
  return hits;
}

} // namespace o2::gtrk::toybarrel8

#endif // GTRK_TEST_TOYBARREL8FIXTURE_H_
