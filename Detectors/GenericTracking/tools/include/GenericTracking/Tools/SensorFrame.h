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

/// \file SensorFrame.h
/// \brief A point of the global frame in a sensor's own, and whether it is on
///        the sensor's silicon: the inverses the tracking does not need, for
///        the truth crossings and the geometry tools.

#ifndef GTRK_TOOLS_SENSORFRAME_H_
#define GTRK_TOOLS_SENSORFRAME_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/Sensor.h"

#include <cmath>

namespace o2::gtrk::tools
{

/// local = rot^T * (global - trans): a rotation's inverse is its transpose.
inline void toLocal(const Sensor& sensor, const double g[3], double l[3])
{
  const double d0 = g[0] - sensor.trans[0];
  const double d1 = g[1] - sensor.trans[1];
  const double d2 = g[2] - sensor.trans[2];
  l[0] = sensor.rot[0] * d0 + sensor.rot[3] * d1 + sensor.rot[6] * d2;
  l[1] = sensor.rot[1] * d0 + sensor.rot[4] * d1 + sensor.rot[7] * d2;
  l[2] = sensor.rot[2] * d0 + sensor.rot[5] * d1 + sensor.rot[8] * d2;
}

/// The inverse of SensorSpec::unrolledToFrame: a point of the sensor's frame
/// as its unrolled coordinates.
inline void frameToUnrolled(const SensorSpec& spec, const double frame[3], double unrolled[3])
{
  if (!spec.curved()) {
    unrolled[0] = frame[0];
    unrolled[1] = frame[1];
    unrolled[2] = frame[2];
    return;
  }
  unrolled[0] = spec.radius * std::atan2(frame[1], frame[0]);
  unrolled[1] = frame[2];
  unrolled[2] = std::hypot(frame[0], frame[1]) - spec.radius;
}

/// Whether a point's unrolled \a x, \a y (a flat sensor's frame; a curved
/// sensor's frame point goes through `frameToUnrolled` first) are on silicon.
[[nodiscard]] constexpr bool contains(const SensorSpec& spec, double x, double y) { return y >= -spec.halfY && y <= spec.halfY && x >= -spec.halfX && x <= spec.halfX; }

} // namespace o2::gtrk::tools

#endif // GTRK_TOOLS_SENSORFRAME_H_
