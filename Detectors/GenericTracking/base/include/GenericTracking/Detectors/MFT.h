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

/// \file MFT.h
/// \brief The ALICE MFT: ten disc faces of ALPIDE chips, placed by O2.
///
/// A placed layout (PlacedLayerSpec) like ITS.h, placed from
/// o2::mft::GeometryTGeo. A layer is one face of a half-disc pair: five discs,
/// two faces each, both halves at the same z. A face covers only part of its
/// disc -- the two faces of a disc together are what is hermetic -- so a track
/// usually crosses one face of a disc, not both, and the links skip faces.
///
/// Chip IDs are O2's, which run half, disc, ladder, chip and so are not
/// grouped by layer. The tables are generated: `o2-gtrk-import-geometry`
/// prints the layers, `o2-gtrk-layer-links --min-share 0.005` the links: the
/// ones carrying under 0.5 % of their source face's crossings are left out,
/// which saves tracking time for few findable tracks.

#ifndef GTRK_DETECTORS_MFT_H_
#define GTRK_DETECTORS_MFT_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/Alpide.h"

#include <array>

namespace o2::gtrk::detectors
{

struct MFT {
  static constexpr const char* name = "MFT";
  static constexpr std::array<PlacedLayerSpec, 0> barrelLayers{};
  // clang-format off
  // Measured values that happen to lie near mathematical constants are not them.
  // NOLINTBEGIN(modernize-use-std-numbers)
  // Chip counts measured by o2-gtrk-import-geometry from an aligned o2sim_geometry-aligned.root: 936 chips.
  static constexpr std::array<PlacedLayerSpec, 10> discLayers{{
    {.sensorCount = 66, .sensor = alpide::kSensor},
    {.sensorCount = 66, .sensor = alpide::kSensor},
    {.sensorCount = 66, .sensor = alpide::kSensor},
    {.sensorCount = 66, .sensor = alpide::kSensor},
    {.sensorCount = 82, .sensor = alpide::kSensor},
    {.sensorCount = 82, .sensor = alpide::kSensor},
    {.sensorCount = 118, .sensor = alpide::kSensor},
    {.sensorCount = 118, .sensor = alpide::kSensor},
    {.sensorCount = 136, .sensor = alpide::kSensor},
    {.sensorCount = 136, .sensor = alpide::kSensor},
  }};

  // Placed discs traced in 3D against the sensors of an aligned o2sim_geometry-aligned.root.
  // Traced by o2-gtrk-layer-links for MFT: pT >= 0.1 GeV at Bz = 5 kG,
  // beam spot z in [-18, 18] cm, |eta| <= 6, source intervals padded by 0.5 cm.
  // skipsLayer is traced from intervening layer extents, not alternative graph paths.
  // The comment on each link is the fraction of its source layer's crossings that take it.
  // Links below 0.005 of their source layer's crossings are left out (19 of 45).
  static constexpr std::array<LayerLink, 26> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = 1.80, .sourceMax = 12.01, .skipsLayer = false}, // 0.7276
    {.from = 0, .to = 2, .sourceMin = 1.80, .sourceMax = 12.12, .skipsLayer = true}, // 0.2596
    {.from = 0, .to = 4, .sourceMin = 3.37, .sourceMax = 11.96, .skipsLayer = true}, // 0.0072
    {.from = 1, .to = 2, .sourceMin = 1.80, .sourceMax = 11.96, .skipsLayer = false}, // 0.7236
    {.from = 1, .to = 3, .sourceMin = 1.80, .sourceMax = 12.13, .skipsLayer = true}, // 0.2591
    {.from = 1, .to = 4, .sourceMin = 3.36, .sourceMax = 12.32, .skipsLayer = true}, // 0.0132
    {.from = 2, .to = 3, .sourceMin = 1.80, .sourceMax = 12.01, .skipsLayer = false}, // 0.7208
    {.from = 2, .to = 4, .sourceMin = 1.80, .sourceMax = 12.26, .skipsLayer = true}, // 0.2622
    {.from = 2, .to = 5, .sourceMin = 3.37, .sourceMax = 12.42, .skipsLayer = true}, // 0.0131
    {.from = 3, .to = 4, .sourceMin = 1.80, .sourceMax = 12.64, .skipsLayer = false}, // 0.7446
    {.from = 3, .to = 5, .sourceMin = 1.80, .sourceMax = 11.83, .skipsLayer = true}, // 0.2489
    {.from = 3, .to = 6, .sourceMin = 3.34, .sourceMax = 12.63, .skipsLayer = true}, // 0.0051
    {.from = 4, .to = 5, .sourceMin = 1.80, .sourceMax = 13.67, .skipsLayer = false}, // 0.7600
    {.from = 4, .to = 6, .sourceMin = 1.87, .sourceMax = 12.88, .skipsLayer = true}, // 0.1696
    {.from = 4, .to = 7, .sourceMin = 1.83, .sourceMax = 13.35, .skipsLayer = true}, // 0.0577
    {.from = 4, .to = 8, .sourceMin = 1.82, .sourceMax = 13.35, .skipsLayer = true}, // 0.0112
    {.from = 5, .to = 6, .sourceMin = 1.95, .sourceMax = 13.35, .skipsLayer = false}, // 0.7651
    {.from = 5, .to = 7, .sourceMin = 1.88, .sourceMax = 12.86, .skipsLayer = false}, // 0.1962
    {.from = 5, .to = 8, .sourceMin = 1.84, .sourceMax = 13.47, .skipsLayer = false}, // 0.0180
    {.from = 5, .to = 9, .sourceMin = 1.80, .sourceMax = 13.41, .skipsLayer = false}, // 0.0207
    {.from = 6, .to = 7, .sourceMin = 2.90, .sourceMax = 16.32, .skipsLayer = false}, // 0.7482
    {.from = 6, .to = 8, .sourceMin = 2.91, .sourceMax = 15.72, .skipsLayer = true}, // 0.2030
    {.from = 6, .to = 9, .sourceMin = 2.94, .sourceMax = 16.03, .skipsLayer = true}, // 0.0489
    {.from = 7, .to = 8, .sourceMin = 2.90, .sourceMax = 16.53, .skipsLayer = false}, // 0.7957
    {.from = 7, .to = 9, .sourceMin = 2.90, .sourceMax = 16.02, .skipsLayer = false}, // 0.2043
    {.from = 8, .to = 9, .sourceMin = 3.30, .sourceMax = 17.14, .skipsLayer = false}, // 1.0000
  }};
  // NOLINTEND(modernize-use-std-numbers)
  // clang-format on
  /// O2's beryllium pipe, inside the first layer.
  static constexpr std::array<PassiveSpec, 1> passives{kAliceBerylliumPipe};
  static constexpr int nSensors = countSensors(barrelLayers, discLayers);
};

} // namespace o2::gtrk::detectors

#endif // GTRK_DETECTORS_MFT_H_
