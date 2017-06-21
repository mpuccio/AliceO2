// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAIndexTable.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "CAIndexTable.h"

#include <algorithm>

#include "CAIndexTableUtils.h"
#include "CAMathUtils.h"

CAIndexTable::CAIndexTable()
    : mLayerIndex { CAConstants::ITS::UnusedIndex }
{
}

CAIndexTable::CAIndexTable(const int layerIndex, const std::vector<CACluster>& clusters)
    : mLayerIndex { layerIndex }
{
  const int layerClustersNum = clusters.size();
  int previousBinIndex = 0;
  mTableBins[0] = 0;

  for (int iCluster = 1; iCluster < layerClustersNum; ++iCluster) {

    const int currentBinIndex = clusters[iCluster].indexTableBinIndex;

    if (currentBinIndex > previousBinIndex) {

      for (int iBin = previousBinIndex + 1; iBin <= currentBinIndex; ++iBin) {

        mTableBins[iBin] = iCluster;
      }

      previousBinIndex = currentBinIndex;
    }
  }

  for (int iBin = previousBinIndex + 1; iBin <= CAConstants::IndexTable::ZBins * CAConstants::IndexTable::PhiBins;
      iBin++) {

    mTableBins[iBin] = layerClustersNum;
  }

}

const std::vector<int> CAIndexTable::selectBins(const float zRangeMin, const float zRangeMax, const float phiRangeMin,
    const float phiRangeMax)
{
  std::vector<int> filteredBins;

  if (zRangeMax < -CAConstants::ITS::LayersZCoordinate[mLayerIndex]
      || zRangeMin > CAConstants::ITS::LayersZCoordinate[mLayerIndex] || zRangeMin > zRangeMax) {

    return filteredBins;
  }

  const int minZBinIndex = std::max(0, CAIndexTableUtils::getZBinIndex(mLayerIndex, zRangeMin));
  const int maxZBinIndex = std::min(CAConstants::IndexTable::ZBins - 1,
      CAIndexTableUtils::getZBinIndex(mLayerIndex, zRangeMax));
  const int zBinsNum = maxZBinIndex - minZBinIndex + 1;
  const int minPhiBinIndex = CAIndexTableUtils::getPhiBinIndex(CAMathUtils::getNormalizedPhiCoordinate(phiRangeMin));
  const int maxPhiBinIndex = CAIndexTableUtils::getPhiBinIndex(CAMathUtils::getNormalizedPhiCoordinate(phiRangeMax));

  int phiBinsNum = maxPhiBinIndex - minPhiBinIndex + 1;

  if (phiBinsNum < 0) {

    phiBinsNum += CAConstants::IndexTable::PhiBins;
  }

  filteredBins.reserve(phiBinsNum * zBinsNum);

  for (int iPhiBin = minPhiBinIndex, iPhiCount = 0; iPhiCount < phiBinsNum;
      iPhiBin = ++iPhiBin == CAConstants::IndexTable::PhiBins ? 0 : iPhiBin, iPhiCount++) {

    const int firstBinIndex = CAIndexTableUtils::getBinIndex(minZBinIndex, iPhiBin);
    const int maxBinIndex = firstBinIndex + zBinsNum;

    for (int iBinIndex = firstBinIndex; iBinIndex < maxBinIndex; ++iBinIndex) {

      if (mTableBins[iBinIndex] != mTableBins[iBinIndex + 1]) {

        filteredBins.emplace_back(iBinIndex);
      }
    }
  }

  return filteredBins;
}
