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
#include <limits>
#include <iomanip>
#include <chrono>
#include "TH1F.h"

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

using Constants::IndexTable::ZBins;
using Constants::IndexTable::PhiBins;
using Constants::Math::TwoPi;

Vertexer::Vertexer(const Event& event) : 
mEvent{event},
mITSRadii{Constants::ITS::LayersRCoordinate()}
{
  mR1 = mITSRadii[1] - mITSRadii[0];
  mR2 = mITSRadii[2] - mITSRadii[1];
  std::array<float, Constants::ITS::LayersNumber> layersZCoordinate = Constants::ITS::LayersZCoordinate();
  
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

void Vertexer::initialize(const float zCut, const float phiCut)
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
    for (int iBin { previousBinIndex + 1 }; iBin <= ZBins * PhiBins; iBin++) {
      mIndexTables[iLayer][iBin] = clustersNum;
    }
  }
  
  mZCut = zCut;
  mPhiCut = phiCut;
  std::array<float, Constants::ITS::LayersNumber> layersZCoordinates = Constants::ITS::LayersZCoordinate();
  mZBinSize = layersZCoordinates[0]/ZBins;
  mPhiSpan = static_cast<int>( std::ceil( mPhiCut/(2 * TwoPi/PhiBins) ));
  mZSpan = static_cast<int>( std::ceil(mZCut/mZBinSize) );
  mVertexerInitialized = true;
}

void Vertexer::printIndexTables()
{
  for (int iTables { 0 }; iTables < Constants::ITS::LayersNumberVertexer; ++iTables) {
    Printf("Table %d", iTables);
    for (int iIndexPhi { 0 }; iIndexPhi < PhiBins; ++iIndexPhi) {
      for (int iIndexZeta { 0 }; iIndexZeta < ZBins; ++iIndexZeta) {
        std::cout<<mIndexTables[iTables][iIndexZeta + ZBins * iIndexPhi]<<"\t";
      }
      Printf("");
    }
    std::cout<<mIndexTables[iTables][ZBins * PhiBins]<<"\t";
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
  TH1F* h = new TH1F("h1", ";zetaproj;N", 5000, -100, 100);
  // 0 - nZbins * mPhiSpan
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  
  for ( int iMiddleTile { 0 }; iMiddleTile < ZBins * PhiBins; ++iMiddleTile ) {
    std::cout<<"middle tile :: "<< iMiddleTile <<std::endl;
    int lowestZInnerTile  { std::max(0, ZBins - static_cast<int>(std::ceil((ZBins - iMiddleTile % ZBins + 1)* (mR2 + mR1)/mR1))) };
    int highestZInnerTile { std::min(static_cast<int>(std::ceil((iMiddleTile % ZBins + 1) * (mR2 + mR1)/mR2 )), ZBins - 1) };
    int RowOffset = static_cast<int>(iMiddleTile/ZBins); // che riga sono? Apro in phi attorno a questa riga, sia in outer layer si in inner layer
    for ( int iMiddleCluster { mIndexTables[1][iMiddleTile] }; iMiddleCluster < mIndexTables[1][iMiddleTile + 1]; ++iMiddleCluster ) {
      for ( int iInnerRow { RowOffset - mPhiSpan }; iInnerRow <= RowOffset + mPhiSpan; ++iInnerRow ) {
        for ( int iInnerTile { lowestZInnerTile }; iInnerTile <= highestZInnerTile; ++iInnerTile) {
          int offsetInnerTile { (iInnerRow < 0) ? ZBins * (PhiBins + iInnerRow) + iInnerTile : iInnerRow * ZBins + iInnerTile };
          for ( int iInnerCluster { mIndexTables[0][offsetInnerTile] }; iInnerCluster < mIndexTables[0][offsetInnerTile + 1]; ++iInnerCluster ) {
            float zetaProjection { ((mClusters[1][iMiddleCluster].zCoordinate - mClusters[0][iInnerCluster].zCoordinate)/(mClusters[1][iMiddleCluster].rCoordinate - mClusters[0][iInnerCluster].rCoordinate)) * (mR2 + mR1) + mClusters[0][iInnerCluster].zCoordinate };
            int projectionTileZ { static_cast<int>(std::ceil( zetaProjection/mZBinSize )) }; // Phi ∈ [-mPhiSpan, +mPhiSpan] 
            //  std::cout<<" --- --- --- --- --- --- --- --- --- --- "<<std::endl;
            for ( int iOuterRow { RowOffset - mPhiSpan }; iOuterRow <= RowOffset + mPhiSpan; ++iOuterRow ) {
              // Il margine sinistro, se negativo, viene avvicinato a zero (bordo del detector), quello destro, nell'eventualita' sia negativo anche lui, non ha significato. 
              // Esco quindi dal loop. In entrambi i casi va quindi bene usare projectionTile + mZSpan
              // std::cout<<"\t";
              int lowestOuterTile { (projectionTileZ - mZSpan < 0) ? 0 : projectionTileZ - mZSpan };
              int highestOuterTile { (projectionTileZ + mZSpan >= ZBins) ? ZBins - 1 : projectionTileZ + mZSpan };
              for ( int iOuterTile { lowestOuterTile }; iOuterTile <= highestOuterTile; ++iOuterTile ) {
                // std::cout<<" "<<iOuterTile;
                // std::cout<<"\tiOuterTile "<<iOuterTile<<std::endl;
                int offsetOuterTile { (iOuterRow < 0) ? ZBins * (PhiBins + iOuterRow) + iOuterTile : ( (iOuterRow > PhiBins) ? iOuterRow - PhiBins : iOuterRow ) * ZBins + iOuterTile };
                // std::cout<<offsetOuterTile<<" ";
                for ( int iOuterCluster { mIndexTables[2][offsetOuterTile] }; iOuterCluster < mIndexTables[2][offsetOuterTile + 1]; ++iOuterCluster ) {
                  float zDiff { mClusters[0][iInnerCluster].zCoordinate - mClusters[2][iOuterCluster].zCoordinate };
                  mTriplets.emplace_back( std::array<int, 3>{ iInnerCluster, iMiddleCluster, iOuterCluster });
                  // if ( zDiff < mZCut ) { h->Fill( zDiff ); };
                } 
              }
              // std::cout<<std::endl;
            }
          }
        }
      }
    }
  }
  std::cout<< "mTriplets size is:" << mTriplets.size() << std::endl;
  end = std::chrono::system_clock::now();
  int elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end); 
    std::cout << "finished computation at " << std::ctime(&end_time) << "elapsed time: " << elapsed_milliseconds << "ms\n";
  h->Draw("SAME");
}

}
}
}