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

/// \file ITS.h
/// \brief The ALICE ITS: seven barrel layers of ALPIDE chips, placed by O2.
///
/// A placed layout (PlacedLayerSpec): the chips' placements come from O2's
/// o2::its::GeometryTGeo at run time (GenericTracking/Alice/Placement.h), and what is
/// stated here is what the tracking needs at compile time -- each layer's chip
/// count and the chip. The tables are
/// generated, not typed: `o2-gtrk-import-geometry` prints the layers,
/// `o2-gtrk-layer-links` the links.

#ifndef GTRK_DETECTORS_ITS_H_
#define GTRK_DETECTORS_ITS_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/Alpide.h"

#include <array>

namespace o2::gtrk::detectors
{

struct ITS {
  static constexpr const char* name = "ITS";
  // clang-format off
  // Measured values that happen to lie near mathematical constants are not them.
  // NOLINTBEGIN(modernize-use-std-numbers)
  // Chip counts measured by o2-gtrk-import-geometry from an aligned o2sim_geometry-aligned.root: 24120 chips.
  static constexpr std::array<PlacedLayerSpec, 7> barrelLayers{{
    {.sensorCount = 108, .sensor = alpide::kSensor},
    {.sensorCount = 144, .sensor = alpide::kSensor},
    {.sensorCount = 180, .sensor = alpide::kSensor},
    {.sensorCount = 2688, .sensor = alpide::kSensor},
    {.sensorCount = 3360, .sensor = alpide::kSensor},
    {.sensorCount = 8232, .sensor = alpide::kSensor},
    {.sensorCount = 9408, .sensor = alpide::kSensor},
  }};

  // Traced by o2-gtrk-layer-links for ITS: pT >= 0.1 GeV at Bz = 5 kG,
  // beam spot z in [-18, 18] cm, |eta| <= 6, source intervals padded by 0.5 cm.
  // skipsLayer is traced from intervening layer extents, not alternative graph paths.
  // The comment on each link is the fraction of its source layer's crossings that take it.
  static constexpr std::array<LayerLink, 12> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = -14.06, .sourceMax = 14.06, .skipsLayer = false}, // 0.9833
    {.from = 0, .to = 3, .sourceMin = -14.06, .sourceMax = 14.06, .skipsLayer = false}, // 0.0126
    {.from = 0, .to = 5, .sourceMin = -13.89, .sourceMax = 13.89, .skipsLayer = false}, // 0.0042
    {.from = 1, .to = 2, .sourceMin = -14.06, .sourceMax = 14.06, .skipsLayer = false}, // 0.9858
    {.from = 1, .to = 3, .sourceMin = -14.06, .sourceMax = 14.06, .skipsLayer = false}, // 0.0133
    {.from = 1, .to = 5, .sourceMin = -12.69, .sourceMax = 12.69, .skipsLayer = false}, // 0.0009
    {.from = 2, .to = 3, .sourceMin = -14.06, .sourceMax = 14.06, .skipsLayer = false}, // 0.9957
    {.from = 2, .to = 5, .sourceMin = -14.05, .sourceMax = 14.05, .skipsLayer = false}, // 0.0043
    {.from = 3, .to = 4, .sourceMin = -37.84, .sourceMax = 37.84, .skipsLayer = false}, // 0.8916
    {.from = 3, .to = 5, .sourceMin = -42.71, .sourceMax = 42.71, .skipsLayer = false}, // 0.1084
    {.from = 4, .to = 5, .sourceMin = -42.71, .sourceMax = 42.71, .skipsLayer = false}, // 1.0000
    {.from = 5, .to = 6, .sourceMin = -67.33, .sourceMax = 67.33, .skipsLayer = false}, // 1.0000
  }};
  // NOLINTEND(modernize-use-std-numbers)
  // clang-format on
  static constexpr std::array<PlacedLayerSpec, 0> discLayers{};
  /// O2's beryllium pipe, inside the first layer.
  static constexpr std::array<PassiveSpec, 1> passives{kAliceBerylliumPipe};
  static constexpr int nSensors = countSensors(barrelLayers, discLayers);
};

} // namespace o2::gtrk::detectors

#endif // GTRK_DETECTORS_ITS_H_
