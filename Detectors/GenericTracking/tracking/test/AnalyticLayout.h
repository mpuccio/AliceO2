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

/// \file AnalyticLayout.h
/// \brief The toy layouts' geometry: staves and disc rings described by a few
///        numbers, their sensors computed here and handed to Geometry as an
///        experiment's are.
///
/// A toy declares its analytic layers as `barrelSpecs` and `discSpecs`, and
/// derives the layout's `barrelLayers`/`discLayers` from them with
/// `placedLayers`. Chip IDs run barrel layers inner to outer, stave by
/// increasing azimuth, module by increasing z; then disc layers by increasing
/// z, ring inner to outer, module by increasing azimuth.

#ifndef GTRK_TEST_ANALYTICLAYOUT_H_
#define GTRK_TEST_ANALYTICLAYOUT_H_

#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/Sensor.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

namespace o2::gtrk
{

/// A ring of sensors at one radius within a disc. The sensors of a ring
/// overlap in azimuth, so alternate modules are displaced by +-gap/2 along the
/// beam, front and back; nModules is even for that alternation to close.
struct RingSpec {
  double r;       ///< radius of the sensor centres
  int nModules;   ///< sensors around the ring
  double phi0;    ///< azimuth of the first sensor
  double zOffset; ///< displacement of the whole ring from the disc's own z
  double gap;     ///< front-to-back separation of alternate modules
};

/// A barrel layer: nStaves staves around the beam, each carrying
/// nModulesPerStave sensors along it.
struct BarrelLayerSpec {
  double r;             ///< radius of the stave centres
  int nStaves;          ///< staves in azimuth
  int nModulesPerStave; ///< sensors along one stave
  double phi0;          ///< azimuth of the first stave
  double tilt;          ///< stave rotation about its own long axis
  double moduleGap;     ///< gap between consecutive sensors along a stave
  SensorSpec sensor;

  [[nodiscard]] constexpr int nSensors() const { return nStaves * nModulesPerStave; }

  /// Half the length spanned by the sensors of one stave.
  [[nodiscard]] constexpr double halfLength() const { return 0.5 * (static_cast<double>(nModulesPerStave) * (2. * sensor.halfY + moduleGap) - moduleGap); }
};

inline constexpr int kMaxRings = 5;

/// A disc layer: one or more concentric rings of one sensor shape at a fixed z.
struct DiscLayerSpec {
  double z; ///< position along the beam
  int nRings;
  RingSpec rings[kMaxRings];
  SensorSpec sensor;

  [[nodiscard]] constexpr int nSensors() const
  {
    int n = 0;
    for (int i = 0; i < nRings; ++i) {
      n += rings[i].nModules;
    }
    return n;
  }
};

/// The layout's layers of \a specs: each one's sensor count and shape.
template <typename Spec, std::size_t N>
[[nodiscard]] constexpr std::array<PlacedLayerSpec, N> placedLayers(const std::array<Spec, N>& specs)
{
  std::array<PlacedLayerSpec, N> layers{};
  for (std::size_t i = 0; i < N; ++i) {
    layers[i] = {.sensorCount = specs[i].nSensors(), .sensor = specs[i].sensor};
  }
  return layers;
}

template <typename Toy>
concept AnalyticToy = requires {
  Toy::barrelSpecs;
  Toy::discSpecs;
};

/// First chip ID of a layer.
template <AnalyticToy Toy>
[[nodiscard]] constexpr SensorID firstSensorOfLayer(LayerID layer)
{
  SensorID first = 0;
  LayerID seen = 0;
  for (const auto& spec : Toy::barrelSpecs) {
    if (seen++ == layer) {
      return first;
    }
    first += static_cast<SensorID>(spec.nSensors());
  }
  for (const auto& spec : Toy::discSpecs) {
    if (seen++ == layer) {
      return first;
    }
    first += static_cast<SensorID>(spec.nSensors());
  }
  return first;
}

/// The toy's sensors, indexed by chip ID.
///
/// A barrel sensor's local frame has y along the beam and z along the outward
/// normal, turned by the stave tilt:
///   z_global = ( cos a,  sin a, 0),  a = phi + tilt
///   y_global = (     0,      0, 1)
///   x_global = y x z = (-sin a, cos a, 0)
/// A disc sensor faces along the beam, its long axis pointing outwards:
///   z_global = (    0,     0, 1)
///   y_global = (cos p, sin p, 0)
///   x_global = y x z = (sin p, -cos p, 0)
template <AnalyticToy Toy>
[[nodiscard]] std::vector<Sensor> placeSensors()
{
  std::vector<Sensor> sensors;
  sensors.reserve(static_cast<std::size_t>(Toy::nSensors));
  LayerID layer = 0;
  for (const auto& spec : Toy::barrelSpecs) {
    const double dPhi = 2. * std::numbers::pi / static_cast<double>(spec.nStaves);
    const double pitchZ = 2. * spec.sensor.halfY + spec.moduleGap;
    const double z0 = -0.5 * static_cast<double>(spec.nModulesPerStave - 1) * pitchZ;
    for (int iStave = 0; iStave < spec.nStaves; ++iStave) {
      const double phi = spec.phi0 + dPhi * static_cast<double>(iStave);
      const double cosA = std::cos(phi + spec.tilt);
      const double sinA = std::sin(phi + spec.tilt);
      for (int iMod = 0; iMod < spec.nModulesPerStave; ++iMod) {
        sensors.push_back({.rot = {-sinA, 0., cosA, cosA, 0., sinA, 0., 1., 0.},
                           .trans = {spec.r * std::cos(phi), spec.r * std::sin(phi), z0 + pitchZ * static_cast<double>(iMod)},
                           .layer = layer});
      }
    }
    ++layer;
  }
  for (const auto& spec : Toy::discSpecs) {
    for (int iRing = 0; iRing < spec.nRings; ++iRing) {
      const RingSpec& ring = spec.rings[iRing];
      const double dPhi = 2. * std::numbers::pi / static_cast<double>(ring.nModules);
      for (int iMod = 0; iMod < ring.nModules; ++iMod) {
        const double phi = ring.phi0 + dPhi * static_cast<double>(iMod);
        const double cosP = std::cos(phi);
        const double sinP = std::sin(phi);
        sensors.push_back({.rot = {sinP, cosP, 0., -cosP, sinP, 0., 0., 0., 1.},
                           .trans = {ring.r * cosP, ring.r * sinP, spec.z + ring.zOffset + (iMod % 2 == 0 ? 0.5 * ring.gap : -0.5 * ring.gap)},
                           .layer = layer});
      }
    }
    ++layer;
  }
  return sensors;
}

/// The toy's geometry: its sensors placed, its envelopes measured and
/// published for the extent functions.
template <AnalyticToy Toy>
[[nodiscard]] Geometry<Toy> toyGeometry()
{
  const auto sensors = placeSensors<Toy>();
  return Geometry<Toy>{sensors};
}

/// Builds the toy's geometry once, which publishes its envelopes; a toy's
/// header calls it in a namespace-scope initialiser, so that the extent
/// functions answer from start-up on, as after an experiment's geometry is
/// loaded.
template <AnalyticToy Toy>
bool publishToy()
{
  static_cast<void>(toyGeometry<Toy>());
  return true;
}

} // namespace o2::gtrk

#endif // GTRK_TEST_ANALYTICLAYOUT_H_
