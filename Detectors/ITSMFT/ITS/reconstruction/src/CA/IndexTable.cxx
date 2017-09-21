// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file IndexTable.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "ITSReconstruction/CA/IndexTable.h"

#include <algorithm>

#include "ITSReconstruction/CA/MathUtils.h"

namespace o2
{
namespace ITS
{
namespace CA
{

IndexTable::IndexTable() : mLayerIndex{ Constants::ITS::UnusedIndex } {}

IndexTable::IndexTable(const int layerIndex, const std::vector<Cluster>& clusters) : mLayerIndex{ layerIndex }
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

  for (int iBin = previousBinIndex + 1; iBin <= Constants::IndexTable::ZBins * Constants::IndexTable::PhiBins;
       iBin++) {
    mTableBins[iBin] = layerClustersNum;
  }
}

const std::vector<int> IndexTable::selectBins(const float zRangeMin, const float zRangeMax, const float phiRangeMin,
                                                const float phiRangeMax)
{
  std::vector<int> filteredBins;

  if (zRangeMax < -Constants::ITS::LayersZCoordinate[mLayerIndex] ||
      zRangeMin > Constants::ITS::LayersZCoordinate[mLayerIndex] || zRangeMin > zRangeMax) {
    return filteredBins;
  }

  const int minZBinIndex = std::max(0, IndexTable::getZBinIndex(mLayerIndex, zRangeMin));
  const int maxZBinIndex =
    std::min(Constants::IndexTable::ZBins - 1, IndexTable::getZBinIndex(mLayerIndex, zRangeMax));
  const int zBinsNum = maxZBinIndex - minZBinIndex + 1;
  const int minPhiBinIndex = IndexTable::getPhiBinIndex(MathUtils::getNormalizedPhiCoordinate(phiRangeMin));
  const int maxPhiBinIndex = IndexTable::getPhiBinIndex(MathUtils::getNormalizedPhiCoordinate(phiRangeMax));

  int phiBinsNum = maxPhiBinIndex - minPhiBinIndex + 1;

  if (phiBinsNum < 0) {
    phiBinsNum += Constants::IndexTable::PhiBins;
  }

  filteredBins.reserve(phiBinsNum * zBinsNum);

  for (int iPhiBin = minPhiBinIndex, iPhiCount = 0; iPhiCount < phiBinsNum;
       iPhiBin = ++iPhiBin == Constants::IndexTable::PhiBins ? 0 : iPhiBin, iPhiCount++) {
    const int firstBinIndex = IndexTable::getBinIndex(minZBinIndex, iPhiBin);
    const int maxBinIndex = firstBinIndex + zBinsNum;

    for (int iBinIndex = firstBinIndex; iBinIndex < maxBinIndex; ++iBinIndex) {
      if (mTableBins[iBinIndex] != mTableBins[iBinIndex + 1]) {
        filteredBins.emplace_back(iBinIndex);
      }
    }
  }

  return filteredBins;
}

}
}
}
