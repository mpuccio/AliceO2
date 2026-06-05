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

#ifndef ALICEO2_ITSMFT_TRACKING_INCLUDE_CONSTANTS_H_
#define ALICEO2_ITSMFT_TRACKING_INCLUDE_CONSTANTS_H_

#include "DetectorsCommonDataFormats/DetID.h"

namespace o2::itsmft::tracking::constants
{

constexpr int ITSNLayers = 7;
constexpr int MFTNLayers = 10;
constexpr int MaxIter = 4; // same as o2::its::constants::MaxIter

constexpr int nLayersForDet(o2::detectors::DetID::ID detId)
{
  return detId == o2::detectors::DetID::MFT ? MFTNLayers : ITSNLayers;
}

} // namespace o2::itsmft::tracking::constants

#endif /* ALICEO2_ITSMFT_TRACKING_INCLUDE_CONSTANTS_H_ */
