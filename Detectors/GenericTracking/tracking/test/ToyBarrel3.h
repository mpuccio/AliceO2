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

/// \file ToyBarrel3.h
/// \brief Three barrel layers, small enough to check by hand.
///
/// Not a stand-in for a real detector: a test fixture.  Twelve sensors in
/// total, so an exhaustive test over every chip is cheap, and the numbers are
/// round enough that a placement can be verified with a pencil.

#ifndef GTRK_DETECTORS_TOYBARREL3_H_
#define GTRK_DETECTORS_TOYBARREL3_H_

#include "AnalyticLayout.h"

#include <array>

namespace o2::gtrk::detectors
{

struct ToyBarrel3 {
  static constexpr const char* name = "ToyBarrel3";

  static constexpr SensorSpec sensor{
    .shape = ShapeType::kRectangle,
    .halfX = 0.5,
    .halfY = 1.0,
    .halfThickness = 0.005,
    .pitchX = 0.0025, // 25 um -> 400 rows
    .pitchY = 0.0025, //          800 columns
  };

  static constexpr std::array<BarrelLayerSpec, 3> barrelSpecs{{
    {.r = 2., .nStaves = 4, .nModulesPerStave = 1, .phi0 = 0., .tilt = 0., .moduleGap = 0., .sensor = sensor},
    {.r = 4., .nStaves = 4, .nModulesPerStave = 1, .phi0 = 0., .tilt = 0., .moduleGap = 0., .sensor = sensor},
    {.r = 8., .nStaves = 4, .nModulesPerStave = 1, .phi0 = 0., .tilt = 0., .moduleGap = 0., .sensor = sensor},
  }};

  static constexpr std::array<DiscLayerSpec, 0> discSpecs{};

  // Traced as o2-gtrk-layer-links traces a layout, kept by hand:
  // pT >= 0.1 GeV at Bz = 20 kG, beam spot z in [-1, 1] cm, |eta| <= 6, source intervals padded by 0.5 cm.
  // The comment on each link is the fraction of its source layer's crossings that take it.
  static constexpr std::array<LayerLink, 2> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = -1.50, .sourceMax = 1.50, .skipsLayer = false}, // 1.0000
    {.from = 1, .to = 2, .sourceMin = -1.50, .sourceMax = 1.50, .skipsLayer = false}, // 1.0000
  }};

  static constexpr auto barrelLayers = placedLayers(barrelSpecs);
  static constexpr auto discLayers = placedLayers(discSpecs);
  static constexpr int nSensors = countSensors(barrelLayers, discLayers);
};

inline const bool kToyBarrel3Published = publishToy<ToyBarrel3>();

} // namespace o2::gtrk::detectors

#endif // GTRK_DETECTORS_TOYBARREL3_H_
