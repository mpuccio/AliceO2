// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CALayer.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "CALayer.h"

#include "CAConstants.h"

CALayer::CALayer()
    : mLayerIndex { CAConstants::ITS::UnusedIndex }
{
  //Nothing to do
}

CALayer::CALayer(const int layerIndex)
    : mLayerIndex { layerIndex }
{
  //Nothing to do
}

void CALayer::addCluster(const int clusterId, const float xCoordinate, const float yCoordinate, const float zCoordinate,
    const float alphaAngle, const int monteCarlo)
{
  mClusters.emplace_back(clusterId, xCoordinate, yCoordinate, zCoordinate, alphaAngle, monteCarlo);
}
