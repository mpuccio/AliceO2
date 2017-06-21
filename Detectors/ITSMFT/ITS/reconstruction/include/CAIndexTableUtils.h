// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAIndexTableUtils.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CAINDEXTABLEUTILS_H_
#define TRACKINGITSU_INCLUDE_CAINDEXTABLEUTILS_H_

#include "CAConstants.h"

namespace CAIndexTableUtils {
constexpr float getInverseZBinSize(const int);
constexpr int getZBinIndex(const float, const float);
constexpr int getPhiBinIndex(const float);
constexpr int getBinIndex(const int, const int);
}
;

constexpr float getInverseZCoordinate(const int layerIndex)
{

  return 0.5 * CAConstants::IndexTable::ZBins / CAConstants::ITS::LayersZCoordinate[layerIndex];
}

constexpr int CAIndexTableUtils::getZBinIndex(const float layerIndex, const float zCoordinate)
{

  return (zCoordinate + CAConstants::ITS::LayersZCoordinate[layerIndex])
      * CAConstants::IndexTable::InverseZBinSize[layerIndex];
}

constexpr int CAIndexTableUtils::getPhiBinIndex(const float currentPhi)
{

  return (currentPhi * CAConstants::IndexTable::InversePhiBinSize);
}

constexpr int CAIndexTableUtils::getBinIndex(const int zIndex, const int phiIndex)
{
  return std::min(phiIndex * CAConstants::IndexTable::PhiBins + zIndex,
      CAConstants::IndexTable::ZBins * CAConstants::IndexTable::PhiBins);
}

#endif /* TRACKINGITSU_INCLUDE_CAINDEXTABLEUTILS_H_ */
