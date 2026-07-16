// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKINGFRAMEINFOADAPTERS_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKINGFRAMEINFOADAPTERS_H_

#include <array>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITStracking/Cluster.h"

namespace o2::itsmft::tracking
{

// Temporary compatibility projection for the detector-facing legacy type.
// In particular, the MFT coordinates below are intentionally not disk-frame
// coordinates: existing production code consumes this exact synthetic layout.
template <o2::detectors::DetID::ID DetId>
inline o2::its::TrackingFrameInfo makeTrackingFrameInfo(const DecodedCluster& decoded)
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    return o2::its::TrackingFrameInfo{
      decoded.global.x, decoded.global.y, decoded.global.z,
      decoded.cylinderFrame.q, decoded.cylinderFrame.frameAngle,
      std::array<float, 2>{decoded.cylinderFrame.u, decoded.cylinderFrame.v},
      std::array<float, 3>{decoded.rowColumnCovariance.uu, decoded.rowColumnCovariance.uv, decoded.rowColumnCovariance.vv}};
  } else {
    return o2::its::TrackingFrameInfo{
      decoded.global.x, decoded.global.y, decoded.global.z,
      decoded.global.x, 0.f,
      std::array<float, 2>{decoded.global.y, decoded.global.z},
      std::array<float, 3>{decoded.rowColumnCovariance.uu, decoded.rowColumnCovariance.uv, decoded.rowColumnCovariance.vv}};
  }
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKINGFRAMEINFOADAPTERS_H_ */
