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

/// \file Sensor.h
/// \brief Placement of one sensor, and the transform between its frame and the
///        global one.
///
/// Plain doubles, no TGeoHMatrix: a Sensor is trivially copyable and an array
/// of them is one contiguous block.  The rotation is stored whole rather than
/// as angles because the callers apply it, they do not interpolate it.
///
/// This is the host-side geometry, and it is deliberately double.  It is built
/// once, so the precision costs nothing, and it keeps alignment and the
/// differencing of large coordinates clear of rounding.  Device code does not
/// use this type: the tracking builds its own float SoA at its boundary, which
/// is where narrow types start to pay for themselves.

#ifndef GTRK_BASE_SENSOR_H_
#define GTRK_BASE_SENSOR_H_

#include "GenericTracking/Base/Definitions.h"

namespace o2::gtrk
{

/// A sensor's position and orientation.
///
/// \c rot is row-major and maps local to global, so its columns are the local
/// axes expressed in global coordinates.
struct Sensor {
  double rot[9];   ///< row-major local -> global rotation
  double trans[3]; ///< sensor centre in the global frame
  LayerID layer;   ///< index into the barrel-then-disc layer numbering

  /// global = rot * local + trans
  void toGlobal(const double l[3], double g[3]) const
  {
    g[0] = rot[0] * l[0] + rot[1] * l[1] + rot[2] * l[2] + trans[0];
    g[1] = rot[3] * l[0] + rot[4] * l[1] + rot[5] * l[2] + trans[1];
    g[2] = rot[6] * l[0] + rot[7] * l[1] + rot[8] * l[2] + trans[2];
  }
};

} // namespace o2::gtrk

#endif // GTRK_BASE_SENSOR_H_
