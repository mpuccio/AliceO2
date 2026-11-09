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

/// \file Alpide.h
/// \brief The ALPIDE chip of the ALICE ITS and MFT, as a sensor shape.
///
/// Only the active matrix: 512 rows of 26.88 um across the chip, 1024
/// columns of 29.24 um along it, and the 28 um effective sensitive
/// thickness. The passive edges are not silicon a cluster can come from, and
/// the adapter centres each surface on the active matrix rather than on the
/// chip (GenericTracking/Alice/Placement.h), so the bounds stay symmetric.
///
/// Our local x is ALICE's local x (across the rows) and our y is ALICE's
/// local z (along the columns). The numbers are O2's
/// ITSMFTBase/SegmentationAlpide.h, restated here in double (O2's are float)
/// so that the layouts need no O2 detector library; Placement.cxx asserts
/// that they agree.

#ifndef GTRK_DETECTORS_ALPIDE_H_
#define GTRK_DETECTORS_ALPIDE_H_

#include "GenericTracking/Base/LayerSpec.h"

namespace o2::gtrk::detectors::alpide
{

inline constexpr int kRows = 512;
inline constexpr int kColumns = 1024;
inline constexpr double kPitchRow = 26.88e-4;        ///< across the chip [cm]
inline constexpr double kPitchColumn = 29.24e-4;     ///< along the chip [cm]
inline constexpr double kEffectiveThickness = 28e-4; ///< [cm]

inline constexpr SensorSpec kSensor{
  .shape = ShapeType::kRectangle,
  .halfX = 0.5 * kRows * kPitchRow,
  .halfY = 0.5 * kColumns * kPitchColumn,
  .halfThickness = 0.5 * kEffectiveThickness,
  .pitchX = kPitchRow,
  .pitchY = kPitchColumn,
};

} // namespace o2::gtrk::detectors::alpide

#endif // GTRK_DETECTORS_ALPIDE_H_
