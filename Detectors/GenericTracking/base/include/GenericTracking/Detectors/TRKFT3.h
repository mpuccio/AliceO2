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

/// \file TRKFT3.h
/// \brief The ALICE 3 tracker: the TRK barrel and the FT3 discs, placed by O2.
///
/// A placed layout (PlacedLayerSpec) like ITS.h, placed from O2's
/// o2::trk::GeometryTGeo. The barrel is
/// eleven layers: the vertex detector's (VD) three, each four curved sensors
/// that are quarter-cylinders about the beam; then the five middle (ML) and
/// three outer (OT) layers of flat chips. Then FT3's twelve discs of the same
/// flat chip, in FT3's order: six at -z from the inside out, then
/// their mirrors at +z. One layout for barrel and discs, because a forward
/// track crosses both.
///
/// The sensors are chip types, restated from O2's TRKBase/Specs.h, which
/// exists only with ENABLE_UPGRADES while the base library builds this
/// layout in every build; workflow/src/TRK.cxx holds them to O2's. Where
/// each chip sits is the loaded geometry's. Sensor IDs are O2's TRK chip IDs
/// for the barrel -- the twelve VD sensors first, petal by petal, then ML and
/// OT chip by chip -- and kFirstDiscSensor plus O2's FT3 chip ID for the discs. The
/// tables are generated: `o2-gtrk-import-geometry -d TRKFT3` prints the chip
/// counts and, with `--passives`, the passive scatterers; `o2-gtrk-layer-links
/// -d TRKFT3` the links.

#ifndef GTRK_DETECTORS_TRKFT3_H_
#define GTRK_DETECTORS_TRKFT3_H_

#include "GenericTracking/Base/LayerSpec.h"

#include <array>
#include <cstddef>
#include <numbers>

namespace o2::gtrk::detectors
{

namespace trkft3
{

// O2's units, so that each value below is O2's to the bit.
inline constexpr double kCm = 1.;
inline constexpr double kMm = 0.1;
inline constexpr double kUm = 1e-4;

/// A VD sensor: a quarter of a cylinder of radius kVDRadius less the gap
/// between petals, 50 cm long, 20 um of silicon, 10 um pixels. O2 measures its
/// flat coordinate on the inner surface (SensorSpec::radius).
inline constexpr std::array<double, 3> kVDRadius{0.5 * kCm, 1.2 * kCm, 2.5 * kCm};
inline constexpr std::array<double, 3> kVDPetalGap{1.63 * kMm, 1.2 * kMm, 1.2 * kMm};
inline constexpr double kVDLength = 50 * kCm;
inline constexpr double kVDThickness = 20 * kUm;
inline constexpr double kVDPitch = 10 * kUm;

/// The arc width of a VD sensor of layer \a layer [cm], as O2 computes it.
[[nodiscard]] constexpr double vdWidth(int layer)
{
  return (kVDRadius[static_cast<std::size_t>(layer)] * 2 * std::numbers::pi / 4) - kVDPetalGap[static_cast<std::size_t>(layer)];
}

[[nodiscard]] constexpr SensorSpec vdSensor(int layer)
{
  return {.shape = ShapeType::kCylinder,
          .halfX = 0.5 * vdWidth(layer),
          .halfY = 0.5 * kVDLength,
          .halfThickness = 0.5 * kVDThickness,
          .pitchX = kVDPitch,
          .pitchY = kVDPitch,
          .radius = kVDRadius[static_cast<std::size_t>(layer)]};
}

/// The ML, OT and FT3 chip: 2.5 x 2.9 cm, of which a 1.5 mm strip along one
/// long side is its readout and has no pixels; 100 um of silicon, 20 um
/// pixels. O2's sensor volume is the active matrix alone, in the barrel and in
/// FT3, so the chip's frame is centred on it.
inline constexpr double kChipWidth = 25 * kMm;
inline constexpr double kChipPassiveEdgeReadOut = 1.5 * kMm;
inline constexpr double kChipLength = 29 * kMm;
inline constexpr double kChipThickness = 100 * kUm;
inline constexpr double kChipPitch = 20 * kUm;
inline constexpr SensorSpec kChip{.shape = ShapeType::kRectangle,
                                  .halfX = 0.5 * (kChipWidth - kChipPassiveEdgeReadOut),
                                  .halfY = 0.5 * kChipLength,
                                  .halfThickness = 0.5 * kChipThickness,
                                  .pitchX = kChipPitch,
                                  .pitchY = kChipPitch};

} // namespace trkft3

struct TRKFT3 {
  static constexpr const char* name = "TRKFT3";
  // clang-format off
  // Chip counts measured by o2-gtrk-import-geometry from an ALICE 3 sgn_geometry-aligned.root: 50876 chips.
  static constexpr std::array<PlacedLayerSpec, 11> barrelLayers{{
    {.sensorCount = 4, .sensor = trkft3::vdSensor(0)},
    {.sensorCount = 4, .sensor = trkft3::vdSensor(1)},
    {.sensorCount = 4, .sensor = trkft3::vdSensor(2)},
    {.sensorCount = 880, .sensor = trkft3::kChip},
    {.sensorCount = 1232, .sensor = trkft3::kChip},
    {.sensorCount = 1584, .sensor = trkft3::kChip},
    {.sensorCount = 2288, .sensor = trkft3::kChip},
    {.sensorCount = 3344, .sensor = trkft3::kChip},
    {.sensorCount = 10560, .sensor = trkft3::kChip},
    {.sensorCount = 13376, .sensor = trkft3::kChip},
    {.sensorCount = 17600, .sensor = trkft3::kChip},
  }};
  // Chip counts measured by o2-gtrk-import-geometry from an ALICE 3 sgn_geometry-aligned.root: 13824 FT3 chips.
  static constexpr std::array<PlacedLayerSpec, 12> discLayers{{
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
  }};
  // Traced by o2-gtrk-layer-links -d TRKFT3 from an ALICE 3 geometry, the discs placed: pT >= 0.1 GeV at Bz = 20 kG
  // (the ALICE 3 field object's at the origin), beam spot z in [-18, 18] cm, |eta| <= 6, source intervals padded by 0.5 cm.
  // skipsLayer is traced from intervening layer extents, not alternative graph paths.
  // The comment on each link is the fraction of its source layer's crossings that take it.
  static constexpr std::array<LayerLink, 140> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = -21.42, .sourceMax = 21.42, .skipsLayer = false}, // 0.9664
    {.from = 0, .to = 3, .sourceMin = -21.75, .sourceMax = 21.75, .skipsLayer = false}, // 0.0041
    {.from = 0, .to = 12, .sourceMin = -22.64, .sourceMax = -18.67, .skipsLayer = false}, // 0.0023
    {.from = 0, .to = 13, .sourceMin = -23.75, .sourceMax = -16.92, .skipsLayer = false}, // 0.0067
    {.from = 0, .to = 14, .sourceMin = -22.73, .sourceMax = -18.67, .skipsLayer = true}, // 0.0001
    {.from = 0, .to = 15, .sourceMin = -23.69, .sourceMax = -18.67, .skipsLayer = true}, // 0.0009
    {.from = 0, .to = 16, .sourceMin = -24.96, .sourceMax = -15.17, .skipsLayer = false}, // 0.0046
    {.from = 0, .to = 18, .sourceMin = 18.67, .sourceMax = 22.64, .skipsLayer = false}, // 0.0023
    {.from = 0, .to = 19, .sourceMin = 16.92, .sourceMax = 23.75, .skipsLayer = false}, // 0.0067
    {.from = 0, .to = 20, .sourceMin = 18.67, .sourceMax = 22.73, .skipsLayer = true}, // 0.0001
    {.from = 0, .to = 21, .sourceMin = 18.67, .sourceMax = 23.69, .skipsLayer = true}, // 0.0009
    {.from = 0, .to = 22, .sourceMin = 15.17, .sourceMax = 24.96, .skipsLayer = false}, // 0.0046
    {.from = 1, .to = 2, .sourceMin = -21.86, .sourceMax = 21.86, .skipsLayer = false}, // 0.8459
    {.from = 1, .to = 3, .sourceMin = -25.50, .sourceMax = 25.50, .skipsLayer = false}, // 0.0796
    {.from = 1, .to = 12, .sourceMin = -25.50, .sourceMax = -11.50, .skipsLayer = false}, // 0.0069
    {.from = 1, .to = 13, .sourceMin = -25.50, .sourceMax = -8.38, .skipsLayer = false}, // 0.0183
    {.from = 1, .to = 14, .sourceMin = -25.50, .sourceMax = -13.05, .skipsLayer = true}, // 0.0004
    {.from = 1, .to = 15, .sourceMin = -25.50, .sourceMax = -9.93, .skipsLayer = true}, // 0.0021
    {.from = 1, .to = 16, .sourceMin = -25.50, .sourceMax = -6.82, .skipsLayer = false}, // 0.0095
    {.from = 1, .to = 18, .sourceMin = 11.50, .sourceMax = 25.50, .skipsLayer = false}, // 0.0069
    {.from = 1, .to = 19, .sourceMin = 8.38, .sourceMax = 25.50, .skipsLayer = false}, // 0.0183
    {.from = 1, .to = 20, .sourceMin = 13.05, .sourceMax = 25.50, .skipsLayer = true}, // 0.0004
    {.from = 1, .to = 21, .sourceMin = 9.93, .sourceMax = 25.50, .skipsLayer = true}, // 0.0021
    {.from = 1, .to = 22, .sourceMin = 6.82, .sourceMax = 25.50, .skipsLayer = false}, // 0.0095
    {.from = 2, .to = 3, .sourceMin = -25.50, .sourceMax = 25.50, .skipsLayer = false}, // 0.9432
    {.from = 2, .to = 12, .sourceMin = -25.50, .sourceMax = -10.19, .skipsLayer = false}, // 0.0031
    {.from = 2, .to = 13, .sourceMin = -25.50, .sourceMax = -10.19, .skipsLayer = false}, // 0.0165
    {.from = 2, .to = 14, .sourceMin = -25.50, .sourceMax = -18.11, .skipsLayer = true}, // 0.0001
    {.from = 2, .to = 15, .sourceMin = -25.50, .sourceMax = -10.19, .skipsLayer = true}, // 0.0016
    {.from = 2, .to = 16, .sourceMin = -25.50, .sourceMax = -10.19, .skipsLayer = false}, // 0.0071
    {.from = 2, .to = 18, .sourceMin = 10.19, .sourceMax = 25.50, .skipsLayer = false}, // 0.0031
    {.from = 2, .to = 19, .sourceMin = 10.19, .sourceMax = 25.50, .skipsLayer = false}, // 0.0165
    {.from = 2, .to = 20, .sourceMin = 18.11, .sourceMax = 25.50, .skipsLayer = true}, // 0.0001
    {.from = 2, .to = 21, .sourceMin = 10.19, .sourceMax = 25.50, .skipsLayer = true}, // 0.0016
    {.from = 2, .to = 22, .sourceMin = 10.19, .sourceMax = 25.50, .skipsLayer = false}, // 0.0071
    {.from = 3, .to = 4, .sourceMin = -54.45, .sourceMax = 54.45, .skipsLayer = false}, // 0.8386
    {.from = 3, .to = 11, .sourceMin = -60.88, .sourceMax = -45.34, .skipsLayer = false}, // 0.0100
    {.from = 3, .to = 12, .sourceMin = -64.49, .sourceMax = -45.34, .skipsLayer = false}, // 0.0484
    {.from = 3, .to = 13, .sourceMin = -64.49, .sourceMax = -45.34, .skipsLayer = true}, // 0.0189
    {.from = 3, .to = 14, .sourceMin = -64.49, .sourceMax = -45.46, .skipsLayer = true}, // 0.0005
    {.from = 3, .to = 15, .sourceMin = -64.49, .sourceMax = -45.37, .skipsLayer = true}, // 0.0009
    {.from = 3, .to = 16, .sourceMin = -64.49, .sourceMax = -45.54, .skipsLayer = true}, // 0.0020
    {.from = 3, .to = 17, .sourceMin = 45.34, .sourceMax = 60.88, .skipsLayer = false}, // 0.0100
    {.from = 3, .to = 18, .sourceMin = 45.34, .sourceMax = 64.49, .skipsLayer = false}, // 0.0484
    {.from = 3, .to = 19, .sourceMin = 45.34, .sourceMax = 64.49, .skipsLayer = true}, // 0.0189
    {.from = 3, .to = 20, .sourceMin = 45.46, .sourceMax = 64.49, .skipsLayer = true}, // 0.0005
    {.from = 3, .to = 21, .sourceMin = 45.37, .sourceMax = 64.49, .skipsLayer = true}, // 0.0009
    {.from = 3, .to = 22, .sourceMin = 45.54, .sourceMax = 64.49, .skipsLayer = true}, // 0.0020
    {.from = 4, .to = 5, .sourceMin = -53.12, .sourceMax = 53.12, .skipsLayer = false}, // 0.7777
    {.from = 4, .to = 11, .sourceMin = -64.49, .sourceMax = -42.63, .skipsLayer = false}, // 0.0820
    {.from = 4, .to = 12, .sourceMin = -64.49, .sourceMax = -42.63, .skipsLayer = true}, // 0.0279
    {.from = 4, .to = 13, .sourceMin = -64.49, .sourceMax = -42.67, .skipsLayer = true}, // 0.0012
    {.from = 4, .to = 14, .sourceMin = -64.49, .sourceMax = -43.27, .skipsLayer = true}, // 0.0000
    {.from = 4, .to = 15, .sourceMin = -64.49, .sourceMax = -43.40, .skipsLayer = true}, // 0.0000
    {.from = 4, .to = 16, .sourceMin = -64.48, .sourceMax = -48.06, .skipsLayer = true}, // 0.0000
    {.from = 4, .to = 17, .sourceMin = 42.63, .sourceMax = 64.49, .skipsLayer = false}, // 0.0820
    {.from = 4, .to = 18, .sourceMin = 42.63, .sourceMax = 64.49, .skipsLayer = true}, // 0.0279
    {.from = 4, .to = 19, .sourceMin = 42.67, .sourceMax = 64.49, .skipsLayer = true}, // 0.0012
    {.from = 4, .to = 20, .sourceMin = 43.27, .sourceMax = 64.49, .skipsLayer = true}, // 0.0000
    {.from = 4, .to = 21, .sourceMin = 43.40, .sourceMax = 64.49, .skipsLayer = true}, // 0.0000
    {.from = 4, .to = 22, .sourceMin = 48.06, .sourceMax = 64.48, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 6, .sourceMin = -46.99, .sourceMax = 46.99, .skipsLayer = false}, // 0.5203
    {.from = 5, .to = 8, .sourceMin = -48.04, .sourceMax = 48.04, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 11, .sourceMin = -64.49, .sourceMax = -30.29, .skipsLayer = false}, // 0.2372
    {.from = 5, .to = 12, .sourceMin = -64.49, .sourceMax = -30.32, .skipsLayer = true}, // 0.0026
    {.from = 5, .to = 13, .sourceMin = -64.49, .sourceMax = -31.27, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 14, .sourceMin = -64.48, .sourceMax = -32.10, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 15, .sourceMin = -63.91, .sourceMax = -34.79, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 16, .sourceMin = -58.92, .sourceMax = -32.99, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 17, .sourceMin = 30.29, .sourceMax = 64.49, .skipsLayer = false}, // 0.2372
    {.from = 5, .to = 18, .sourceMin = 30.32, .sourceMax = 64.49, .skipsLayer = true}, // 0.0026
    {.from = 5, .to = 19, .sourceMin = 31.27, .sourceMax = 64.49, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 20, .sourceMin = 32.10, .sourceMax = 64.48, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 21, .sourceMin = 34.79, .sourceMax = 63.91, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 22, .sourceMin = 32.99, .sourceMax = 58.92, .skipsLayer = true}, // 0.0000
    {.from = 6, .to = 7, .sourceMin = -47.81, .sourceMax = 47.81, .skipsLayer = false}, // 0.1654
    {.from = 6, .to = 8, .sourceMin = -64.49, .sourceMax = 64.49, .skipsLayer = false}, // 0.1386
    {.from = 6, .to = 11, .sourceMin = -64.49, .sourceMax = -25.80, .skipsLayer = false}, // 0.3430
    {.from = 6, .to = 12, .sourceMin = -64.49, .sourceMax = -26.95, .skipsLayer = true}, // 0.0011
    {.from = 6, .to = 13, .sourceMin = -64.45, .sourceMax = -37.92, .skipsLayer = true}, // 0.0000
    {.from = 6, .to = 14, .sourceMin = -64.49, .sourceMax = -30.59, .skipsLayer = false}, // 0.0038
    {.from = 6, .to = 15, .sourceMin = -64.48, .sourceMax = -33.05, .skipsLayer = true}, // 0.0001
    {.from = 6, .to = 16, .sourceMin = -64.46, .sourceMax = -36.83, .skipsLayer = true}, // 0.0000
    {.from = 6, .to = 17, .sourceMin = 25.80, .sourceMax = 64.49, .skipsLayer = false}, // 0.3430
    {.from = 6, .to = 18, .sourceMin = 26.95, .sourceMax = 64.49, .skipsLayer = true}, // 0.0011
    {.from = 6, .to = 19, .sourceMin = 37.92, .sourceMax = 64.45, .skipsLayer = true}, // 0.0000
    {.from = 6, .to = 20, .sourceMin = 30.59, .sourceMax = 64.49, .skipsLayer = false}, // 0.0038
    {.from = 6, .to = 21, .sourceMin = 33.05, .sourceMax = 64.48, .skipsLayer = true}, // 0.0001
    {.from = 6, .to = 22, .sourceMin = 36.83, .sourceMax = 64.46, .skipsLayer = true}, // 0.0000
    {.from = 7, .to = 8, .sourceMin = -64.49, .sourceMax = 64.49, .skipsLayer = false}, // 0.9844
    {.from = 7, .to = 11, .sourceMin = -64.49, .sourceMax = -47.69, .skipsLayer = false}, // 0.0068
    {.from = 7, .to = 14, .sourceMin = -64.48, .sourceMax = -59.59, .skipsLayer = false}, // 0.0010
    {.from = 7, .to = 17, .sourceMin = 47.69, .sourceMax = 64.49, .skipsLayer = false}, // 0.0068
    {.from = 7, .to = 20, .sourceMin = 59.59, .sourceMax = 64.48, .skipsLayer = false}, // 0.0010
    {.from = 8, .to = 9, .sourceMin = -101.31, .sourceMax = 101.31, .skipsLayer = false}, // 0.1297
    {.from = 8, .to = 14, .sourceMin = -128.88, .sourceMax = -70.15, .skipsLayer = false}, // 0.4250
    {.from = 8, .to = 15, .sourceMin = -128.88, .sourceMax = -77.47, .skipsLayer = true}, // 0.0101
    {.from = 8, .to = 16, .sourceMin = -128.88, .sourceMax = -87.08, .skipsLayer = true}, // 0.0001
    {.from = 8, .to = 20, .sourceMin = 70.15, .sourceMax = 128.88, .skipsLayer = false}, // 0.4250
    {.from = 8, .to = 21, .sourceMin = 77.47, .sourceMax = 128.88, .skipsLayer = true}, // 0.0101
    {.from = 8, .to = 22, .sourceMin = 87.08, .sourceMax = 128.88, .skipsLayer = true}, // 0.0001
    {.from = 9, .to = 10, .sourceMin = -101.25, .sourceMax = 101.25, .skipsLayer = false}, // 0.1171
    {.from = 9, .to = 14, .sourceMin = -128.88, .sourceMax = -91.29, .skipsLayer = false}, // 0.4406
    {.from = 9, .to = 15, .sourceMin = -128.88, .sourceMax = -110.80, .skipsLayer = true}, // 0.0009
    {.from = 9, .to = 20, .sourceMin = 91.29, .sourceMax = 128.88, .skipsLayer = false}, // 0.4406
    {.from = 9, .to = 21, .sourceMin = 110.80, .sourceMax = 128.88, .skipsLayer = true}, // 0.0009
    {.from = 11, .to = 8, .sourceMin = 23.87, .sourceMax = 35.17, .skipsLayer = false}, // 0.0952
    {.from = 11, .to = 12, .sourceMin = 9.50, .sourceMax = 33.31, .skipsLayer = false}, // 0.7973
    {.from = 11, .to = 13, .sourceMin = 9.50, .sourceMax = 30.66, .skipsLayer = true}, // 0.0066
    {.from = 11, .to = 14, .sourceMin = 9.56, .sourceMax = 35.14, .skipsLayer = false}, // 0.0969
    {.from = 11, .to = 15, .sourceMin = 9.53, .sourceMax = 34.13, .skipsLayer = true}, // 0.0039
    {.from = 11, .to = 16, .sourceMin = 9.96, .sourceMax = 32.62, .skipsLayer = true}, // 0.0001
    {.from = 12, .to = 8, .sourceMin = 33.24, .sourceMax = 35.16, .skipsLayer = false}, // 0.0000
    {.from = 12, .to = 13, .sourceMin = 9.50, .sourceMax = 33.83, .skipsLayer = false}, // 0.8259
    {.from = 12, .to = 14, .sourceMin = 10.57, .sourceMax = 35.17, .skipsLayer = false}, // 0.1645
    {.from = 12, .to = 15, .sourceMin = 9.57, .sourceMax = 35.16, .skipsLayer = true}, // 0.0088
    {.from = 12, .to = 16, .sourceMin = 9.55, .sourceMax = 34.99, .skipsLayer = true}, // 0.0008
    {.from = 13, .to = 14, .sourceMin = 12.06, .sourceMax = 35.17, .skipsLayer = false}, // 0.6962
    {.from = 13, .to = 15, .sourceMin = 9.90, .sourceMax = 35.16, .skipsLayer = false}, // 0.1845
    {.from = 13, .to = 16, .sourceMin = 9.51, .sourceMax = 35.05, .skipsLayer = false}, // 0.1193
    {.from = 14, .to = 15, .sourceMin = 15.18, .sourceMax = 68.75, .skipsLayer = false}, // 0.9588
    {.from = 14, .to = 16, .sourceMin = 15.18, .sourceMax = 63.29, .skipsLayer = true}, // 0.0412
    {.from = 15, .to = 16, .sourceMin = 15.18, .sourceMax = 68.29, .skipsLayer = false}, // 1.0000
    {.from = 17, .to = 8, .sourceMin = 23.87, .sourceMax = 35.17, .skipsLayer = false}, // 0.0952
    {.from = 17, .to = 18, .sourceMin = 9.50, .sourceMax = 33.31, .skipsLayer = false}, // 0.7973
    {.from = 17, .to = 19, .sourceMin = 9.50, .sourceMax = 30.66, .skipsLayer = true}, // 0.0066
    {.from = 17, .to = 20, .sourceMin = 9.56, .sourceMax = 35.14, .skipsLayer = false}, // 0.0969
    {.from = 17, .to = 21, .sourceMin = 9.53, .sourceMax = 34.13, .skipsLayer = true}, // 0.0039
    {.from = 17, .to = 22, .sourceMin = 9.96, .sourceMax = 32.62, .skipsLayer = true}, // 0.0001
    {.from = 18, .to = 8, .sourceMin = 33.24, .sourceMax = 35.16, .skipsLayer = false}, // 0.0000
    {.from = 18, .to = 19, .sourceMin = 9.50, .sourceMax = 33.83, .skipsLayer = false}, // 0.8259
    {.from = 18, .to = 20, .sourceMin = 10.57, .sourceMax = 35.17, .skipsLayer = false}, // 0.1645
    {.from = 18, .to = 21, .sourceMin = 9.57, .sourceMax = 35.16, .skipsLayer = true}, // 0.0088
    {.from = 18, .to = 22, .sourceMin = 9.55, .sourceMax = 34.99, .skipsLayer = true}, // 0.0008
    {.from = 19, .to = 20, .sourceMin = 12.06, .sourceMax = 35.17, .skipsLayer = false}, // 0.6962
    {.from = 19, .to = 21, .sourceMin = 9.90, .sourceMax = 35.16, .skipsLayer = false}, // 0.1845
    {.from = 19, .to = 22, .sourceMin = 9.51, .sourceMax = 35.05, .skipsLayer = false}, // 0.1193
    {.from = 20, .to = 21, .sourceMin = 15.18, .sourceMax = 68.75, .skipsLayer = false}, // 0.9588
    {.from = 20, .to = 22, .sourceMin = 15.18, .sourceMax = 63.29, .skipsLayer = true}, // 0.0412
    {.from = 21, .to = 22, .sourceMin = 15.18, .sourceMax = 68.29, .skipsLayer = false}, // 1.0000
  }};
  // Measured by o2-gtrk-import-geometry -d TRKFT3 --passives from O2's full ALICE 3 geometry: one per gap between
  // barrel layers, its x/X0 across the gap at 65536 random azimuths and z, its radius the x/X0-weighted mean of
  // 0.01 cm steps at 2048. Each owns its gap's material; a layer, its silicon envelope's. The gap from 20 to 30 cm
  // holds the inner TOF at 21 cm (2.0 % X0). A gap lies between the layers' chips, their sensors' metal stacks
  // included (tools::withStacks; the OT's first layer has its stack on the inner face): the stacks' material is the
  // layers', not the gaps'.
  static constexpr std::array<PassiveSpec, 11> passives{{
    {.kind = PassiveSpec::Kind::kBarrel, .position = 0.4899, .bandLow = 0.0000, .bandHigh = 0.4999, .low = -25.00, .high = 25.00, .x2X0 = 2.224e-03},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 0.8002, .bandLow = 0.5096, .bandHigh = 1.1999, .low = -25.00, .high = 25.00, .x2X0 = 1.228e-03},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 1.7708, .bandLow = 1.2096, .bandHigh = 2.5000, .low = -25.00, .high = 25.00, .x2X0 = 9.985e-04},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 3.7955, .bandLow = 2.5096, .bandHigh = 6.8657, .low = -63.98, .high = 63.98, .x2X0 = 1.132e-02},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 7.9710, .bandLow = 7.8977, .bandHigh = 8.8055, .low = -63.98, .high = 63.98, .x2X0 = 3.692e-04},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 10.0133, .bandLow = 9.8542, .bandHigh = 11.7623, .low = -63.98, .high = 63.98, .x2X0 = 4.094e-04},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 13.8492, .bandLow = 12.7722, .bandHigh = 18.7300, .low = -63.98, .high = 63.98, .x2X0 = 4.863e-04},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 21.2092, .bandLow = 20.1379, .bandHigh = 29.9999, .low = -63.98, .high = 63.98, .x2X0 = 2.027e-02},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 38.2067, .bandLow = 31.0792, .bandHigh = 44.8000, .low = -128.38, .high = 128.38, .x2X0 = 2.557e-02},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 52.1226, .bandLow = 45.3527, .bandHigh = 59.8000, .low = -128.38, .high = 128.38, .x2X0 = 4.899e-04},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 69.6851, .bandLow = 60.2898, .bandHigh = 79.8000, .low = -128.38, .high = 128.38, .x2X0 = 6.376e-04},
  }};
  // clang-format on
  static constexpr int nSensors = countSensors(barrelLayers, discLayers);
  /// The first disc sensor's ID: FT3's chip ID 0, after every barrel chip.
  static constexpr int kFirstDiscSensor = countSensors(barrelLayers, std::array<PlacedLayerSpec, 0>{});
};

} // namespace o2::gtrk::detectors

#endif // GTRK_DETECTORS_TRKFT3_H_
