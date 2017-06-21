// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAPrimaryVertexContext.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "CAPrimaryVertexContext.h"

#include <algorithm>
#include <cmath>

#include "CAConstants.h"
#include "CAEvent.h"
#include "CALayer.h"

CAPrimaryVertexContext::CAPrimaryVertexContext(const CAEvent& event, const int primaryVertexIndex)
    : primaryVertexIndex { primaryVertexIndex }
{
  for (int iLayer = 0; iLayer < CAConstants::ITS::LayersNumber; ++iLayer) {

    const CALayer& currentLayer = event.getLayer(iLayer);
    const int clustersNum = currentLayer.getClustersSize();
    clusters[iLayer].reserve(clustersNum);

    for (int iCluster = 0; iCluster < clustersNum; ++iCluster) {

      const CACluster& currentCluster = currentLayer.getCluster(iCluster);
      clusters[iLayer].emplace_back(iLayer, event.getPrimaryVertex(primaryVertexIndex), currentCluster);
    }

    std::sort(clusters[iLayer].begin(), clusters[iLayer].end(), [](CACluster& cluster1, CACluster& cluster2) {
      return cluster1.indexTableBinIndex < cluster2.indexTableBinIndex;
    });

    if (iLayer > 0) {

      indexTables[iLayer - 1] = CAIndexTable(iLayer, clusters[iLayer]);
    }

    if (iLayer < CAConstants::ITS::TrackletsPerRoad) {

      tracklets[iLayer].reserve(
          std::ceil(
              (CAConstants::Memory::TrackletsMemoryCoefficients[iLayer] * clustersNum)
                  * event.getLayer(iLayer + 1).getClustersSize()));
    }

    if (iLayer < CAConstants::ITS::CellsPerRoad) {

      cells[iLayer].reserve(
          std::ceil(
              ((CAConstants::Memory::CellsMemoryCoefficients[iLayer] * clustersNum)
                  * event.getLayer(iLayer + 1).getClustersSize()) * event.getLayer(iLayer + 2).getClustersSize()));
    }
  }
}

