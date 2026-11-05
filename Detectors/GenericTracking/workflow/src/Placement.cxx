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

#include "GenericTracking/Alice/Placement.h"

#include "GenericTracking/Detectors/Alpide.h"

#include "ITSMFTBase/SegmentationAlpide.h"

namespace o2::gtrk::alice
{

namespace
{

using Segmentation = o2::itsmft::SegmentationAlpide;

// Alpide.h restates O2's numbers so that the layouts need no O2 detector
// library; this is where they are held to them.
static_assert(detectors::alpide::kRows == Segmentation::NRows);
static_assert(detectors::alpide::kColumns == Segmentation::NCols);
static_assert(exactlyEqual(static_cast<float>(detectors::alpide::kPitchRow), Segmentation::PitchRow));
static_assert(exactlyEqual(static_cast<float>(detectors::alpide::kPitchColumn), Segmentation::PitchCol));
static_assert(exactlyEqual(static_cast<float>(detectors::alpide::kEffectiveThickness), Segmentation::SensorLayerThicknessEff));

/// The rotation and translation of one L2G matrix.
struct Matrix {
  double m[3][3];
  double t[3];
};

Matrix componentsOf(const o2::math_utils::Transform3D& l2g)
{
  Matrix matrix{};
  auto& m = matrix.m;
  auto& t = matrix.t;
  l2g.GetComponents(m[0][0], m[0][1], m[0][2], t[0], m[1][0], m[1][1], m[1][2], t[1], m[2][0], m[2][1], m[2][2], t[2]);
  return matrix;
}

} // namespace

const double kActiveCentreX = 0.5 * (static_cast<double>(Segmentation::PassiveEdgeReadOut) - static_cast<double>(Segmentation::PassiveEdgeTop));

Sensor flatSensor(const o2::math_utils::Transform3D& l2g, int layer, double activeCentreX)
{
  const auto [m, t] = componentsOf(l2g);
  Sensor sensor{};
  for (int i = 0; i < 3; ++i) {
    // Row-major local -> global, so column j is the global direction of our
    // local axis j.
    sensor.rot[(3 * i) + 0] = m[i][0];
    sensor.rot[(3 * i) + 1] = m[i][2];
    sensor.rot[(3 * i) + 2] = -m[i][1];
    sensor.trans[i] = t[i] + (activeCentreX * m[i][0]);
  }
  sensor.layer = static_cast<LayerID>(layer);
  return sensor;
}

Sensor curvedSensor(const o2::math_utils::Transform3D& l2g, int layer)
{
  const auto [m, t] = componentsOf(l2g);
  Sensor sensor{};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      sensor.rot[(3 * i) + j] = m[i][j];
    }
    sensor.trans[i] = t[i];
  }
  sensor.layer = static_cast<LayerID>(layer);
  return sensor;
}

} // namespace o2::gtrk::alice
