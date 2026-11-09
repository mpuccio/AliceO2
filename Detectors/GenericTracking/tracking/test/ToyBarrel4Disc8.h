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

/// \file ToyBarrel4Disc8.h
/// \brief Four barrel layers and four discs on each side: the analytic layout
///        with both, for the tests.
///
/// The barrel's staves are tilted and carry several modules with gaps between
/// them, so some links skip a layer through the gaps; the inner three layers
/// end before the discs and the fourth reaches beyond them, so tracks run from
/// the barrel into the discs and from the discs back into the barrel. The discs
/// lie beyond the 3 sigma of the default 6 cm vertex spread, as the window
/// bounds require (a disc inside the vertex range bounds no kink). Each disc
/// has two rings, displaced front and back along the beam, reaching out to
/// just inside the fourth layer. Discs are ordered along the beam, from -z to
/// +z.

#ifndef GTRK_DETECTORS_TOYBARREL4DISC8_H_
#define GTRK_DETECTORS_TOYBARREL4DISC8_H_

#include "AnalyticLayout.h"

#include <array>

namespace o2::gtrk::detectors
{

namespace toy_barrel4_disc8
{

inline constexpr SensorSpec kBarrelSensor{
  .shape = ShapeType::kRectangle,
  .halfX = 0.5,
  .halfY = 2.0,
  .halfThickness = 0.005,
  .pitchX = 0.0025,
  .pitchY = 0.0025,
};

/// One sensor shape for both rings of a disc: a layout's layer has one.
inline constexpr SensorSpec kDiscSensor{.shape = ShapeType::kRectangle, .halfX = 1.8, .halfY = 4.1, .halfThickness = 0.005, .pitchX = 0.0025, .pitchY = 0.0025};

constexpr DiscLayerSpec disc(double z)
{
  return DiscLayerSpec{
    .z = z,
    .nRings = 2,
    .rings = {{.r = 6.6, .nModules = 24, .phi0 = 0., .zOffset = 0.2, .gap = 0.1}, {.r = 14.4, .nModules = 48, .phi0 = 0., .zOffset = -0.2, .gap = 0.1}},
    .sensor = kDiscSensor,
  };
}

} // namespace toy_barrel4_disc8

struct ToyBarrel4Disc8 {
  static constexpr const char* name = "ToyBarrel4Disc8";

  static constexpr std::array<BarrelLayerSpec, 4> barrelSpecs{{
    {.r = 3., .nStaves = 22, .nModulesPerStave = 5, .phi0 = 0., .tilt = 0.1, .moduleGap = 0.05, .sensor = toy_barrel4_disc8::kBarrelSensor},
    {.r = 6., .nStaves = 42, .nModulesPerStave = 5, .phi0 = 0., .tilt = 0.1, .moduleGap = 0.05, .sensor = toy_barrel4_disc8::kBarrelSensor},
    {.r = 10., .nStaves = 68, .nModulesPerStave = 7, .phi0 = 0., .tilt = 0.1, .moduleGap = 0.05, .sensor = toy_barrel4_disc8::kBarrelSensor},
    {.r = 20., .nStaves = 132, .nModulesPerStave = 19, .phi0 = 0., .tilt = 0.1, .moduleGap = 0.05, .sensor = toy_barrel4_disc8::kBarrelSensor},
  }};

  static constexpr std::array<DiscLayerSpec, 8> discSpecs{{
    toy_barrel4_disc8::disc(-34.),
    toy_barrel4_disc8::disc(-30.),
    toy_barrel4_disc8::disc(-26.),
    toy_barrel4_disc8::disc(-22.),
    toy_barrel4_disc8::disc(22.),
    toy_barrel4_disc8::disc(26.),
    toy_barrel4_disc8::disc(30.),
    toy_barrel4_disc8::disc(34.),
  }};

  // Traced as o2-gtrk-layer-links traces a layout, kept by hand:
  // pT >= 0.1 GeV at Bz = 20 kG, beam spot z in [-5, 5] cm, |eta| <= 6, source intervals padded by 0.5 cm.
  // The comment on each link is the fraction of its source layer's crossings that take it.
  static constexpr std::array<LayerLink, 25> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = -8.05, .sourceMax = 8.05, .skipsLayer = false},   // 0.6690
    {.from = 0, .to = 2, .sourceMin = -8.24, .sourceMax = 8.24, .skipsLayer = false},   // 0.0021
    {.from = 0, .to = 3, .sourceMin = -8.27, .sourceMax = -7.23, .skipsLayer = true},   // 0.0001
    {.from = 0, .to = 7, .sourceMin = -10.60, .sourceMax = -2.02, .skipsLayer = false}, // 0.1644
    {.from = 0, .to = 8, .sourceMin = 2.02, .sourceMax = 10.60, .skipsLayer = false},   // 0.1645
    {.from = 1, .to = 2, .sourceMin = -10.60, .sourceMax = 10.60, .skipsLayer = false}, // 0.8932
    {.from = 1, .to = 7, .sourceMin = -10.60, .sourceMax = -5.88, .skipsLayer = false}, // 0.0534
    {.from = 1, .to = 8, .sourceMin = 5.88, .sourceMax = 10.60, .skipsLayer = false},   // 0.0534
    {.from = 2, .to = 3, .sourceMin = -14.65, .sourceMax = 14.53, .skipsLayer = false}, // 0.8738
    {.from = 2, .to = 7, .sourceMin = -14.65, .sourceMax = -8.53, .skipsLayer = false}, // 0.0601
    {.from = 2, .to = 8, .sourceMin = 8.33, .sourceMax = 14.65, .skipsLayer = false},   // 0.0661
    {.from = 4, .to = 3, .sourceMin = 16.96, .sourceMax = 19.09, .skipsLayer = false},  // 1.0000
    {.from = 5, .to = 3, .sourceMin = 15.54, .sourceMax = 19.09, .skipsLayer = false},  // 0.0559
    {.from = 5, .to = 4, .sourceMin = 2.00, .sourceMax = 17.38, .skipsLayer = false},   // 0.9441
    {.from = 6, .to = 3, .sourceMin = 15.14, .sourceMax = 19.09, .skipsLayer = false},  // 0.0624
    {.from = 6, .to = 5, .sourceMin = 2.00, .sourceMax = 17.19, .skipsLayer = false},   // 0.9376
    {.from = 7, .to = 3, .sourceMin = 14.58, .sourceMax = 19.09, .skipsLayer = false},  // 0.0704
    {.from = 7, .to = 6, .sourceMin = 2.00, .sourceMax = 16.94, .skipsLayer = false},   // 0.9296
    {.from = 8, .to = 3, .sourceMin = 14.52, .sourceMax = 19.09, .skipsLayer = false},  // 0.0699
    {.from = 8, .to = 9, .sourceMin = 2.00, .sourceMax = 16.91, .skipsLayer = false},   // 0.9301
    {.from = 9, .to = 3, .sourceMin = 15.09, .sourceMax = 19.09, .skipsLayer = false},  // 0.0621
    {.from = 9, .to = 10, .sourceMin = 2.00, .sourceMax = 17.17, .skipsLayer = false},  // 0.9379
    {.from = 10, .to = 3, .sourceMin = 15.51, .sourceMax = 19.09, .skipsLayer = false}, // 0.0557
    {.from = 10, .to = 11, .sourceMin = 2.00, .sourceMax = 17.37, .skipsLayer = false}, // 0.9443
    {.from = 11, .to = 3, .sourceMin = 16.72, .sourceMax = 19.09, .skipsLayer = false}, // 1.0000
  }};

  static constexpr auto barrelLayers = placedLayers(barrelSpecs);
  static constexpr auto discLayers = placedLayers(discSpecs);
  static constexpr int nSensors = countSensors(barrelLayers, discLayers);
};

inline const bool kToyBarrel4Disc8Published = publishToy<ToyBarrel4Disc8>();

} // namespace o2::gtrk::detectors

#endif // GTRK_DETECTORS_TOYBARREL4DISC8_H_
