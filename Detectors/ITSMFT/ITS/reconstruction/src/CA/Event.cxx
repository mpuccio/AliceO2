// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CA/Event.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "ITSReconstruction/CA/Event.h"
#include "ITSMFTReconstruction/Cluster.h"

#include <iostream>

namespace o2
{
namespace ITS
{
namespace CA
{

Event::Event(const int eventId, const GeometryTGeo* geo) : mEventId{ eventId }, mGeom { geo }
{
  for (int iLayer = 0; iLayer < Constants::ITS::LayersNumber; ++iLayer) {
    mLayers[iLayer] = Layer(iLayer);
  }
}

void Event::loadClusters(const TClonesArray& clusters) {
  int numOfClusters = clusters.GetEntriesFast();
  if (numOfClusters == 0) {
    std::cout << "No clusters to load !" << std::endl;
    return;
  }

  for (int iCluster = 0; iCluster < numOfClusters; iCluster++) {
    ITSMFT::Cluster* cluster = (ITSMFT::Cluster*)clusters.UncheckedAt(iCluster);
    cluster->SetUniqueID(iCluster);
    const int layerId = mGeom->getLayer(cluster->getSensorID());
    float alphaRef = mGeom->getSensorRefAlpha( cluster->getSensorID());
    mLayers[layerId].addTrackingFrameInfo(cluster->getX(),alphaRef,std::array<float, 2>{cluster->getY(),cluster->getZ()},std::array<float, 3>{cluster->getSigmaY2(), cluster->getSigmaYZ(), cluster->getSigmaZ2()});
    ///FIXME: use Point3D interface for internal CA coordinate representation?
    auto xyz = cluster->getXYZGloRot(*mGeom);
    float globalXYZ[3]{0.f};
    xyz.GetCoordinates(globalXYZ);
    mLayers[layerId].addCluster(iCluster, globalXYZ);
  }
}

void Event::addPrimaryVertex(const float xCoordinate, const float yCoordinate, const float zCoordinate)
{
  mPrimaryVertices.emplace_back(std::array<float, 3>{ { xCoordinate, yCoordinate, zCoordinate } });
}

void Event::printPrimaryVertices() const
{
  const int verticesNum = mPrimaryVertices.size();

  for (int iVertex = 0; iVertex < verticesNum; ++iVertex) {
    const std::array<float, 3>& currentVertex = mPrimaryVertices[iVertex];

    std::cout << "-1\t" << currentVertex[0] << "\t" << currentVertex[1] << "\t" << currentVertex[2] << std::endl;
  }
}

const int Event::getTotalClusters(const int numberOfLayers) const
{
  int totalClusters = 0;

  for (int iLayer = 0; iLayer < numberOfLayers; ++iLayer) {
    totalClusters += mLayers[iLayer].getClustersSize();
  }

  return totalClusters;
}

}
}
}
