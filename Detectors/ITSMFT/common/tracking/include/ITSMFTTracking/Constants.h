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
///
/// \file Constants.h
/// \brief Detector-specific layer counts for shared ITSMFT CA tracking
///
/// MFT CA layer model (see also TimeFrameMFT and MFT CA workflow):
/// - MFTNLayers = 10 half-disk CA layers (same index as GeometryTGeo::getLayer and
///   MFTTracking/Constants.h LayersNumber). This is the single NLayers used by
///   TrackingParameters, TimeFrame, ROFOverlapTable, and cluster sorting.
/// - Five physical disks are addressed as disk = halfLayer / 2 (see mftDiskForHalfLayer).
///   Do not use disk count as CA NLayers.
/// - ROFOverlapTable<MFTNLayers> stores one LayerTiming per half-layer, filled from
///   MFTAlpideParam index L (roFrameLayer*InBC[L] with fallback to global defaults).
/// - Cluster input is one CLUSTERSROF vector per TF; loadROFrameData(..., layer=-1, ...)
///   distributes clusters to half-layers. All layers must then share the same mNROFsTF;
///   per-layer ROF length stagger is unsupported until the workflow provides per-layer ROFs.

#ifndef ALICEO2_ITSMFT_TRACKING_INCLUDE_CONSTANTS_H_
#define ALICEO2_ITSMFT_TRACKING_INCLUDE_CONSTANTS_H_

#include "DetectorsCommonDataFormats/DetID.h"

namespace o2::itsmft::tracking::constants
{

constexpr int ITSNLayers = 7;
constexpr int MFTNLayers = 10;
constexpr int MFTDisks = 5;
constexpr int MaxIter = 4; // same as o2::its::constants::MaxIter

constexpr int nLayersForDet(o2::detectors::DetID::ID detId)
{
  return detId == o2::detectors::DetID::MFT ? MFTNLayers : ITSNLayers;
}

/// Physical disk index (0..MFTDisks-1) for half-layer L in [0, MFTNLayers).
constexpr int mftDiskForHalfLayer(int halfLayer)
{
  return halfLayer / 2;
}

template <int NLayers>
constexpr o2::detectors::DetID::ID detIdFromNLayers()
{
  return NLayers == MFTNLayers ? o2::detectors::DetID::MFT : o2::detectors::DetID::ITS;
}

} // namespace o2::itsmft::tracking::constants

#endif /* ALICEO2_ITSMFT_TRACKING_INCLUDE_CONSTANTS_H_ */
