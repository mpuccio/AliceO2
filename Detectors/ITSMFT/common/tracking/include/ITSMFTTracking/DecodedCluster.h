// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_DECODEDCLUSTER_H_
#define ALICEO2_ITSMFT_TRACKING_DECODEDCLUSTER_H_

#include <cstdint>

#include "ITSMFTTracking/SurfaceMeasurement.h"

namespace o2::itsmft::tracking
{

// Host-side facts produced by compact-cluster and geometry decoding. This is
// deliberately independent of detector output and legacy tracking types.
struct DecodedCluster {
  GlobalPoint3F global{};
  // ITS geometry supplies its cylindrical tracking frame here. Disk
  // projection uses global coordinates directly.
  SurfaceFramePoint cylinderFrame{};
  // ALPIDE local row/column covariance. The detector projection determines
  // which normalized axes these values describe.
  SurfaceCovariance2F rowColumnCovariance{};
  ClusterShape shape{};
  uint32_t sensor{0};
  int layer{0};
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_DECODEDCLUSTER_H_ */
