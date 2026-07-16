// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENTADAPTERS_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENTADAPTERS_H_

#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITStracking/Cluster.h"

namespace o2::itsmft::tracking
{

// Host compatibility helper for the established ITS cylindrical-frame data.
inline SurfaceMeasurement makeCylinderSurfaceMeasurement(const o2::its::TrackingFrameInfo& info,
                                                         DetectorSensorId sensor,
                                                         SurfaceId surface,
                                                         ClusterRef cluster,
                                                         uint32_t sourceROF,
                                                         ClusterShape shape)
{
  return SurfaceMeasurement{
    {info.xCoordinate, info.yCoordinate, info.zCoordinate},
    {info.xTrackingFrame, info.positionTrackingFrame[0], info.positionTrackingFrame[1], info.alphaTrackingFrame},
    {info.covarianceTrackingFrame[0], info.covarianceTrackingFrame[1], info.covarianceTrackingFrame[2]},
    sensor,
    cluster,
    shape,
    sourceROF,
    surface};
}

// A disk measurement must be built from explicitly decoded global x/y
// covariance. The synthetic MFT TrackingFrameInfo is intentionally not an
// accepted input because its position and covariance axes are inconsistent.
inline SurfaceMeasurement makeDiskSurfaceMeasurement(GlobalPoint3F global,
                                                     SurfaceCovariance2F globalXYCovariance,
                                                     DetectorSensorId sensor,
                                                     SurfaceId surface,
                                                     ClusterRef cluster,
                                                     uint32_t sourceROF,
                                                     ClusterShape shape)
{
  return SurfaceMeasurement{
    global,
    {global.z, global.x, global.y, 0.f},
    globalXYCovariance,
    sensor,
    cluster,
    shape,
    sourceROF,
    surface};
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENTADAPTERS_H_ */
