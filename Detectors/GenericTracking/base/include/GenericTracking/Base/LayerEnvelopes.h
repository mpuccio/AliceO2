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

/// \file LayerEnvelopes.h
/// \brief Where a layout's layers are: measured from the sensors it was
///        loaded with, and held for the process.
///
/// A layout (the ITS, the MFT, the TRK) has no compile-time positions: where its
/// chips sit is the loaded geometry's business, aligned as that geometry is.
/// Constructing its Geometry measures each layer's envelope from the sensors
/// and publishes the table here; the extent functions of
/// Tracking/SearchWindow.h, and with them the grid, the search windows, the
/// scattering table and the material sampling, read it. Nothing on the device
/// does. The compiled layout keeps topology only: the layer count, each layer's
/// kind and sensor count, the links.

#ifndef GTRK_BASE_LAYERENVELOPES_H_
#define GTRK_BASE_LAYERENVELOPES_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/Sensor.h"

#include <array>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace o2::gtrk
{

/// The rounding a geometry measures its envelopes with [cm]: 1 um, the
/// minima down, the maxima up, the position to the nearest.
inline constexpr double kEnvelopeResolution = 1e-4;

/// The hull in r and z of one sensor's active silicon, with its effective
/// thickness. A curved sensor's axis must be the beam's (std::runtime_error).
struct SensorExtent {
  double rMin = std::numeric_limits<double>::infinity();
  double rMax = 0.;
  double zMin = std::numeric_limits<double>::infinity();
  double zMax = -std::numeric_limits<double>::infinity();
};

[[nodiscard]] SensorExtent extentOf(const Sensor& sensor, const SensorSpec& spec);

/// The envelope of each layer's active silicon: the hull in r and z of every
/// sensor's active matrix with its effective thickness, and as position the
/// mean radius of its sensors' centres on a barrel or their mean z on a disc.
/// \a specs holds each layer's sensor shape, layers numbered barrel first:
/// the first \a nBarrelLayers are barrel layers. A curved sensor's centre is
/// the middle of its silicon at the middle of its arc; its axis must be the
/// beam's. With a \a resolution, the minima are rounded down to it, the maxima
/// up and the position to the nearest; 0 leaves the values exact. Throws
/// std::runtime_error for a sensor beyond the layers, a layer without one or a
/// curved sensor whose axis is not the beam's.
[[nodiscard]] std::vector<LayerEnvelope> measureEnvelopes(std::span<const Sensor> sensors, std::span<const SensorSpec> specs, int nBarrelLayers,
                                                          double resolution = 0.);

/// One line per layer, for a log; layers numbered barrel first.
[[nodiscard]] std::string describeEnvelopes(std::span<const LayerEnvelope> envelopes, int nBarrelLayers);

/// The envelopes of the layout \a Det loaded in this process.
///
/// One table per layout per process, as O2's `GeometryTGeo::Instance()` holds
/// the one loaded detector for everything that asks: `Geometry<Det>` publishes
/// at construction, a later geometry of the same layout replaces the table,
/// and reading before any exists throws std::logic_error. Published once at
/// set-up, before any tracking thread runs, and read-only after.
template <DetectorLayout Det>
class PlacedEnvelopes
{
 public:
  static void publish(std::span<const LayerEnvelope> envelopes)
  {
    if (envelopes.size() != sTable.size()) {
      throw std::invalid_argument(std::string{Det::name} + ": " + std::to_string(envelopes.size()) + " envelopes for " + std::to_string(sTable.size()) +
                                  " layers");
    }
    for (std::size_t layer = 0; layer < sTable.size(); ++layer) {
      sTable[layer] = envelopes[layer];
    }
    sPublished = true;
  }

  [[nodiscard]] static const LayerEnvelope& get(LayerID layer)
  {
    if (!sPublished) {
      throw std::logic_error(std::string{Det::name} + ": no geometry loaded, so no layer envelopes; construct Geometry<" + Det::name + "> first");
    }
    if (layer >= sTable.size()) {
      throw std::out_of_range(std::string{Det::name} + ": layer " + std::to_string(layer) + " of " + std::to_string(sTable.size()));
    }
    return sTable[layer];
  }

  [[nodiscard]] static std::span<const LayerEnvelope> all()
  {
    if (!sPublished) {
      throw std::logic_error(std::string{Det::name} + ": no geometry loaded, so no layer envelopes; construct Geometry<" + Det::name + "> first");
    }
    return sTable;
  }

 private:
  inline static std::array<LayerEnvelope, Detector<Det>::nLayers> sTable{};
  inline static bool sPublished = false;
};

} // namespace o2::gtrk

#endif // GTRK_BASE_LAYERENVELOPES_H_
