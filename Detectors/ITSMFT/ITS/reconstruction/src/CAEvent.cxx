// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAEvent.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "CAEvent.h"

#include <iostream>

CAEvent::CAEvent(const int eventId)
    : mEventId { eventId }
{
  for (int iLayer = 0; iLayer < CAConstants::ITS::LayersNumber; ++iLayer) {

    mLayers[iLayer] = CALayer(iLayer);
  }
}

void CAEvent::addPrimaryVertex(const float xCoordinate, const float yCoordinate, const float zCoordinate)
{
  mPrimaryVertices.emplace_back(std::array<float, 3> { { xCoordinate, yCoordinate, zCoordinate } });
}

void CAEvent::printPrimaryVertices() const
{
  const int verticesNum = mPrimaryVertices.size();

  for (int iVertex = 0; iVertex < verticesNum; ++iVertex) {

    const std::array<float, 3>& currentVertex = mPrimaryVertices[iVertex];

    std::cout << "-1\t" << currentVertex[0] << "\t" << currentVertex[1] << "\t" << currentVertex[2] << std::endl;
  }
}

void CAEvent::pushClusterToLayer(const int layerIndex, const int clusterId, const float xCoordinate,
    const float yCoordinate, const float zCoordinate, const float aplhaAngle, const int monteCarlo)
{
  mLayers[layerIndex].addCluster(clusterId, xCoordinate, yCoordinate, zCoordinate, aplhaAngle, monteCarlo);
}

const int CAEvent::getTotalClusters() const
{
  int totalClusters = 0;

  for (int iLayer = 0; iLayer < CAConstants::ITS::LayersNumber; ++iLayer) {

    totalClusters += mLayers[iLayer].getClustersSize();
  }

  return totalClusters;
}
