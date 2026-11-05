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

/// \file Cluster.h
/// \brief A cluster as the converters write it from O2's: the measurement
///        producer's input (GTRK/<DET>CLUSTERS).

#ifndef GTRK_CLUSTERING_CLUSTER_H_
#define GTRK_CLUSTERING_CLUSTER_H_

#include "GenericTracking/Base/Definitions.h"

#include <cstdint>

namespace o2::gtrk
{

/// One reconstructed measurement on a sensor.
///
/// A plain record; the measurement producer turns them into the measurement
/// columns once, with the geometry.
///
/// A cluster is a sensor measurement and carries no global position: that is
/// the sensor ID and the local position through the geometry, and the
/// measurement producer computes it once, for every consumer downstream.
struct Cluster {
  SensorID chipID = kInvalidSensorID;

  float xLocal = 0.F; ///< binary centre of gravity in the sensor frame [cm]
  float yLocal = 0.F;

  /// The errors in the sensor frame [cm]: the topology dictionary's
  /// (Clusters.h, TRK.h).
  float sigmaXLocal = 0.F;
  float sigmaYLocal = 0.F;

  /// The readout frame's record as a time box `[timePS, timePS + widthPS)` in
  /// picoseconds from the timeframe's start, its edges rounded from BC
  /// (TimeEstPS.h). Written by every converter; the
  /// measurement producer moves it by the layer's frame start (FrameTiming.h:
  /// the AlpideParam's delay and bias, TRK's `TRKAlmiraParam`).
  std::uint32_t timePS = 0;
  std::uint32_t widthPS = 0;

  /// The pattern's rows (across, local x) and columns (along, local z); 0
  /// where the converter has no pattern. The producer
  /// turns them into the measurement's slope interval (ClusterShapeBand.h).
  std::uint8_t rowSpan = 0;
  std::uint8_t colSpan = 0;
};

} // namespace o2::gtrk

#endif // GTRK_CLUSTERING_CLUSTER_H_
