// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file PrimaryVertexContext.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "ITSReconstruction/CA/PrimaryVertexContext.h"

#include <algorithm>
#include <cmath>

#include "ITSReconstruction/CA/Constants.h"
#include "ITSReconstruction/CA/Event.h"
#include "ITSReconstruction/CA/Layer.h"

namespace o2
{
namespace ITS
{
namespace CA
{

PrimaryVertexContext::PrimaryVertexContext(const Event& event, const int primaryVertexIndex)
  : primaryVertexIndex{ primaryVertexIndex }
{
  for (int iLayer = 0; iLayer < Constants::ITS::LayersNumber; ++iLayer) {
    const Layer& currentLayer = event.getLayer(iLayer);
    const int clustersNum = currentLayer.getClustersSize();
    clusters[iLayer].reserve(clustersNum);
    usedClusters[iLayer].resize(clustersNum,false);

    for (int iCluster = 0; iCluster < clustersNum; ++iCluster) {
      const Cluster& currentCluster = currentLayer.getCluster(iCluster);
      clusters[iLayer].emplace_back(iLayer, event.getPrimaryVertex(primaryVertexIndex), currentCluster);
    }

    ///TODO: Benchmark using separate vector of indices to reduce the memory per cluster at the cost of sorting and additional vector<int>
    std::sort(clusters[iLayer].begin(), clusters[iLayer].end(), [](Cluster& cluster1, Cluster& cluster2) {
      return cluster1.indexTableBinIndex < cluster2.indexTableBinIndex;
    });

    if (iLayer > 0) {
      indexTables[iLayer - 1] = IndexTable(iLayer, clusters[iLayer]);
    }

    if (iLayer < Constants::ITS::TrackletsPerRoad) {
      tracklets[iLayer].reserve(std::ceil((Constants::Memory::TrackletsMemoryCoefficients[iLayer] * clustersNum) *
                                          event.getLayer(iLayer + 1).getClustersSize()));
    }

    if (iLayer < Constants::ITS::CellsPerRoad) {
      cells[iLayer].reserve(std::ceil(((Constants::Memory::CellsMemoryCoefficients[iLayer] * clustersNum) *
                                       event.getLayer(iLayer + 1).getClustersSize()) *
                                       event.getLayer(iLayer + 2).getClustersSize()));
    }
  }
}

}
}
}
