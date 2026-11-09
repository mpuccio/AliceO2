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

/// \file Geometry.h
/// \brief The runtime half of a detector description: where each sensor sits.
///
/// The sensors are the experiment's, handed in whole at construction (for the
/// ALICE detectors from O2's GeometryTGeo, Alice/Placement.h). They are
/// runtime rather than constexpr data because alignment corrections apply to
/// them.
///
/// The array extent is Det::nSensors, a compile-time constant, so the whole
/// geometry is one contiguous fixed-size block.
///
/// Everything here is double.  See Sensor.h for why, and for where the float
/// view the tracking wants gets built instead.

#ifndef GTRK_BASE_GEOMETRY_H_
#define GTRK_BASE_GEOMETRY_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/LayerEnvelopes.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/Sensor.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <span>
#include <stdexcept>

namespace o2::gtrk
{

/// Sensor placements of one detector, indexed by chip ID. The sensors need
/// not be grouped by layer.
template <DetectorLayout Det>
class Geometry
{
 public:
  /// Take the layout's sensors, indexed by chip ID. Checks the count and that
  /// every layer index exists, then measures each layer's envelope from the
  /// sensors (rounded outwards to kEnvelopeResolution) and publishes the
  /// table as this process's PlacedEnvelopes<Det>, which the extent functions
  /// of Tracking/SearchWindow.h read. \a depths: each layer's
  /// `measurementDepth`, or none for 0 on every layer. Throws
  /// std::invalid_argument on the counts, std::runtime_error on a layer
  /// without a sensor.
  explicit Geometry(std::span<const Sensor> sensors, std::span<const double> depths = {})
  {
    if (sensors.size() != static_cast<std::size_t>(Det::nSensors)) {
      throw std::invalid_argument("placed geometry: sensor count differs from the layout's");
    }
    if (!depths.empty() && depths.size() != static_cast<std::size_t>(Detector<Det>::nLayers)) {
      throw std::invalid_argument("placed geometry: measurement depths for another number of layers");
    }
    std::ranges::copy(depths, mDepths.begin());
    if (std::ranges::any_of(sensors, [](const Sensor& s) { return s.layer >= Detector<Det>::nLayers; })) {
      throw std::invalid_argument("placed geometry: layer index outside the layout");
    }
    std::ranges::copy(sensors, mSensors.begin());
    PlacedEnvelopes<Det>::publish(measureEnvelopes(sensors, layerSpecs(), Detector<Det>::nBarrelLayers, kEnvelopeResolution));
  }

  [[nodiscard]] const Sensor& sensor(SensorID id) const { return mSensors[id]; }
  [[nodiscard]] const std::array<Sensor, Det::nSensors>& sensors() const { return mSensors; }

  /// Where layer \a layer's measurements lie in its sensors' silicon: the
  /// unrolled z of the plane whose crossing a cluster's position follows
  /// [cm], about a flat sensor's mid-plane and outwards from a curved one's
  /// radius (SensorSpec::unrolledToFrame). The digitiser collects charge from
  /// part of the depth only; 0 unless
  /// given.
  [[nodiscard]] double measurementDepth(LayerID layer) const { return mDepths[layer]; }

  void toGlobal(SensorID id, const double l[3], double g[3]) const { mSensors[id].toGlobal(l, g); }

  /// The sensor shape of each layer, in layer order.
  [[nodiscard]] static constexpr std::array<SensorSpec, Detector<Det>::nLayers> layerSpecs()
  {
    std::array<SensorSpec, Detector<Det>::nLayers> specs{};
    for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
      specs[static_cast<std::size_t>(layer)] = layerSpec(static_cast<LayerID>(layer));
    }
    return specs;
  }

  /// The shape of a layer's sensors.
  [[nodiscard]] static constexpr const SensorSpec& layerSpec(LayerID layer) { return Detector<Det>::isBarrelLayer(layer) ? Det::barrelLayers[layer].sensor : Det::discLayers[Detector<Det>::discIndex(layer)].sensor; }

  /// The shape of the sensor a chip carries.
  [[nodiscard]] const SensorSpec& specOf(SensorID id) const { return layerSpec(mSensors[id].layer); }

 private:
  std::array<Sensor, Det::nSensors> mSensors{};
  std::array<double, Detector<Det>::nLayers> mDepths{}; ///< [cm]
};

} // namespace o2::gtrk

#endif // GTRK_BASE_GEOMETRY_H_
