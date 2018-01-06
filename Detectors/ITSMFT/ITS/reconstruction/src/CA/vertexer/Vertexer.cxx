// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Vertexer.cxx
/// \brief
/// \author matteo.concas@cern.ch
/// \author maximiliano.puccio@cern.ch

#include <algorithm>
#include "TString.h"
#include <iostream>
#include <cmath>

#include "ITSReconstruction/CA/Constants.h"
#include "ITSReconstruction/CA/Cluster.h"
#include "ITSReconstruction/CA/Event.h"
#include "ITSReconstruction/CA/Layer.h"
#include "ITSReconstruction/CA/vertexer/Vertexer.h"

namespace o2
{
namespace ITS
{
namespace CA
{

Vertexer::Vertexer(const Event& event) : mEvent{event}
{
  for(int iLayer { 0 }; iLayer < Constants::ITS::LayersNumberVertexer; ++iLayer) {
    const Layer& currentLayer { event.getLayer(iLayer) };
    const int clustersNum { currentLayer.getClustersSize() };
    mClusters[iLayer].clear();
    if(clustersNum > mClusters[iLayer].capacity()) {
      mClusters[iLayer].reserve(clustersNum);
    }
    for (int iCluster { 0 }; iCluster < clustersNum; ++iCluster) {
      mClusters[iLayer].emplace_back(currentLayer.getCluster(iCluster));
    }
  }
}

Vertexer::~Vertexer() {};

void Vertexer::initialize()
{
  for (int iLayer { 0 }; iLayer < Constants::ITS::LayersNumberVertexer; ++iLayer) {
    std::sort(mClusters[iLayer].begin(), mClusters[iLayer].end(), [](Cluster& cluster1, Cluster& cluster2) {
      return cluster1.indexTableBinIndex < cluster2.indexTableBinIndex;
    });
  }

  // Index Tables
  for (int iLayer { 0 }; iLayer < Constants::ITS::LayersNumberVertexer; ++iLayer) {
    const int clustersNum = static_cast<int>(mClusters[iLayer].size());
    int previousBinIndex { 0 };
     mIndexTables[iLayer][0] = 0;
    for (int iCluster { 0 }; iCluster < clustersNum; ++iCluster) {
      const int currentBinIndex { mClusters[iLayer][iCluster].indexTableBinIndex };
      if (currentBinIndex > previousBinIndex) {
        for (int iBin { previousBinIndex + 1 }; iBin <= currentBinIndex; ++iBin) {
          mIndexTables[iLayer][iBin] = iCluster;
        }
        previousBinIndex = currentBinIndex;
      }
    }
    // Fill remaining array with latest cluster index on the array of clusters
    for (int iBin { previousBinIndex + 1 }; iBin <= Constants::IndexTable::ZBins * Constants::IndexTable::PhiBins;
        iBin++) {
      mIndexTables[iLayer][iBin] = clustersNum;
    }
  }
  mVertexerInitialized = true;
}

void Vertexer::printIndexTables()
{
  for (int iTables { 0 }; iTables < Constants::ITS::LayersNumberVertexer; ++iTables) {
    Printf("Table %d", iTables);
    for (int iIndexPhi { 0 }; iIndexPhi < Constants::IndexTable::PhiBins; ++iIndexPhi)
    {
      for (int iIndexZeta { 0 }; iIndexZeta < Constants::IndexTable::ZBins; ++iIndexZeta){
        std::cout<<mIndexTables[iTables][iIndexZeta + Constants::IndexTable::ZBins * iIndexPhi]<<"\t";
      }
      Printf("");
    }
  }
}

void Vertexer::debugVertexerData()
{
  for (int iLayer { 0 }; iLayer < Constants::ITS::LayersNumberVertexer; ++iLayer) {
    for (int iCluster { 0 }; iCluster < mClusters[iLayer].size(); ++iCluster ) {
      Printf(" id: %d | it: %d", mClusters[iLayer][iCluster].clusterId, mClusters[iLayer][iCluster].indexTableBinIndex);
    }
  }
}

void Vertexer::computeTriplets()
{
  // A "tile" is one cell of an Index Table
  if (mVertexerInitialized) {
    // Loop over middle index table tiles
    for (int iMiddleTileIndex { 0 }; iMiddleTileIndex < Constants::IndexTable::ZBins; ++iMiddleTileIndex) {
       const int zTilesSpan = static_cast<int>(static_cast<float>((Constants::IndexTable::ZBins - 1)/2) - std::abs(static_cast<float>((Constants::IndexTable::ZBins - 1)/2) - iMiddleTileIndex));

      // Loop over middle index table clusters per tile
      for (int iMiddleClusterIndex = mIndexTables[1][iMiddleTileIndex]; iMiddleClusterIndex < mIndexTables[1][iMiddleTileIndex + 1]; ++iMiddleTileIndex) {

        for (int iTile { iMiddleTileIndex - zTilesSpan }; iTile <= iMiddleTileIndex + zTilesSpan; ++iTile) {
          const int offsetTile = Constants::IndexTable::ZBins * (Constants::IndexTable::PhiBins - 1) + iTile;
          for (int iUpperClusterIndex { mIndexTables[2][offsetTile]}; iUpperClusterIndex < mIndexTables[2][offsetTile + 1]; ++iUpperClusterIndex) {
            for (int iLowerClusterIndex { mIndexTables[2][iTile + Constants::IndexTable::ZBins]}; iLowerClusterIndex < mIndexTables[2][iTile + Constants::IndexTable::ZBins + 1]; ++iLowerClusterIndex) {
              mTriplets.emplace_back(std::array<Cluster, 3>{mClusters[2][iUpperClusterIndex], mClusters[1][iMiddleClusterIndex], mClusters[0][iLowerClusterIndex]});
            }
          }
          for (int iUpperClusterIndex { mIndexTables[2][iTile]}; iUpperClusterIndex < mIndexTables[2][iTile + 1]; ++iUpperClusterIndex) {
            for (int iLowerClusterIndex { mIndexTables[2][iTile]}; iLowerClusterIndex < mIndexTables[2][iTile + 1]; ++iLowerClusterIndex) {
              mTriplets.emplace_back(std::array<Cluster, 3>{mClusters[2][iUpperClusterIndex], mClusters[1][iMiddleClusterIndex], mClusters[0][iLowerClusterIndex]});
            }
          }
          for (int iUpperClusterIndex { mIndexTables[2][iTile + Constants::IndexTable::ZBins]}; iUpperClusterIndex < mIndexTables[2][iTile + Constants::IndexTable::ZBins + 1]; ++iUpperClusterIndex) {
            for (int iLowerClusterIndex { mIndexTables[2][offsetTile]}; iLowerClusterIndex < mIndexTables[2][offsetTile + 1]; ++iLowerClusterIndex) {
              mTriplets.emplace_back(std::array<Cluster, 3>{mClusters[2][iUpperClusterIndex], mClusters[1][iMiddleClusterIndex], mClusters[0][iLowerClusterIndex]});
            }
          }
        }
      }
    }

    for (int iMiddleTileIndex { Constants::IndexTable::ZBins }; iMiddleTileIndex < Constants::IndexTable::PhiBins * Constants::IndexTable::ZBins; ++iMiddleTileIndex) {
      const int zTilesSpan = static_cast<int>(static_cast<float>((Constants::IndexTable::ZBins - 1)/2) - std::abs(static_cast<float>((Constants::IndexTable::ZBins - 1)/2) - iMiddleTileIndex));

      // Loop over middle index table clusters per tile
      for (int iMiddleClusterIndex = mIndexTables[1][iMiddleTileIndex]; iMiddleClusterIndex < mIndexTables[1][iMiddleTileIndex + 1]; ++iMiddleTileIndex) {

        // Loop over index tables rows, taken three by three and opposite with the opposite tables
        for (int iTableRow { -1 }; iTableRow < 2; ++iTableRow ) {
          const int startUpperTile = iTableRow * Constants::IndexTable::ZBins + iMiddleTileIndex - zTilesSpan;
          const int finishUpperTile = iTableRow * Constants::IndexTable::ZBins + iMiddleTileIndex + zTilesSpan;
          const int startLowerTile = - iTableRow * Constants::IndexTable::ZBins + iMiddleTileIndex - zTilesSpan;
          const int finishLowerTile = - iTableRow * Constants::IndexTable::ZBins + iMiddleTileIndex + zTilesSpan;

          for (int iUpperTileIndex { startUpperTile }; iUpperTileIndex <= finishUpperTile; ++iUpperTileIndex) {
            for (int iUpperClusterIndex { mIndexTables[2][iUpperTileIndex] }; iUpperClusterIndex < mIndexTables[2][iUpperTileIndex + 1]; ++iUpperTileIndex) {
              for (int iLowerTileIndex { startLowerTile }; iLowerTileIndex <= finishLowerTile; ++iLowerTileIndex) {
                for (int iLowerClusterIndex { mIndexTables[0][iLowerTileIndex] }; iLowerClusterIndex < mIndexTables[0][iLowerTileIndex + 1]; ++iLowerTileIndex) {
                  mTriplets.emplace_back(std::array<Cluster, 3>{mClusters[2][iUpperClusterIndex], mClusters[1][iMiddleClusterIndex], mClusters[0][iLowerClusterIndex]});
                }
              }
            }
          }
        }
      }
    }
  }
}

}
}
}