// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENTADAPTERS_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENTADAPTERS_H_

#include "ITSMFTTracking/DecodedCluster.h"

namespace o2::itsmft::tracking
{

// Project decoded ITS facts into the accepted cylindrical convention.
inline SurfaceMeasurement makeCylinderSurfaceMeasurement(const DecodedCluster& decoded,
                                                         DetectorSensorId sensor,
                                                         SurfaceId surface,
                                                         ClusterRef cluster,
                                                         uint32_t sourceROF)
{
  return SurfaceMeasurement{
    decoded.global,
    decoded.cylinderFrame,
    decoded.rowColumnCovariance,
    sensor,
    cluster,
    decoded.shape,
    sourceROF,
    surface};
}

// Project decoded MFT facts into z-normal, global-x/global-y disk coordinates.
// ALPIDE row is established as global x and column as global y by the MFT
// geometry decoder. No legacy TrackingFrameInfo participates in this mapping.
inline SurfaceMeasurement makeDiskSurfaceMeasurement(const DecodedCluster& decoded,
                                                     DetectorSensorId sensor,
                                                     SurfaceId surface,
                                                     ClusterRef cluster,
                                                     uint32_t sourceROF)
{
  return SurfaceMeasurement{
    decoded.global,
    {decoded.global.z, decoded.global.x, decoded.global.y, 0.f},
    decoded.rowColumnCovariance,
    sensor,
    cluster,
    decoded.shape,
    sourceROF,
    surface};
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENTADAPTERS_H_ */
