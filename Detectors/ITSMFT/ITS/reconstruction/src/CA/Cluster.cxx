// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CA/Cluster.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.chė
#include "ITSReconstruction/CA/Cluster.h"

#include "ITSReconstruction/CA/IndexTable.h"
#include "ITSReconstruction/CA/MathUtils.h"

namespace o2
{
namespace ITS
{
namespace CA
{

Cluster::Cluster(const int clusterId, const float xCoord, const float yCoord, const float zCoord, const int layerIndex)
  : clusterId{ clusterId },
    xCoordinate{ xCoord },
    yCoordinate{ yCoord },
    zCoordinate{ zCoord },
    phiCoordinate{ MathUtils::getNormalizedPhiCoordinate(MathUtils::calculatePhiCoordinate(xCoord, yCoord)) },
    rCoordinate{ MathUtils::calculateRCoordinate(xCoord, yCoord) },
    indexTableBinIndex{ IndexTable::getBinIndex(IndexTable::getZBinIndex(layerIndex, zCoord),
                                                IndexTable::getPhiBinIndex(phiCoordinate)) }
{
}

Cluster::Cluster(const int layerIndex, const std::array<float, 3>& primaryVertex, const Cluster& other)
  : clusterId{ other.clusterId },
    xCoordinate{ other.xCoordinate },
    yCoordinate{ other.yCoordinate },
    zCoordinate{ other.zCoordinate },
    phiCoordinate{ MathUtils::getNormalizedPhiCoordinate(
      MathUtils::calculatePhiCoordinate(xCoordinate - primaryVertex[0], yCoordinate - primaryVertex[1])) },
    rCoordinate{ MathUtils::calculateRCoordinate(xCoordinate - primaryVertex[0], yCoordinate - primaryVertex[1]) },
    indexTableBinIndex{ IndexTable::getBinIndex(IndexTable::getZBinIndex(layerIndex, zCoordinate),
                                                IndexTable::getPhiBinIndex(phiCoordinate)) }
{
}

TrackingFrameInfo::TrackingFrameInfo(float xTF, float alpha, std::array<float,2>&& posTF, std::array<float,3>&& covTF) :
xTrackingFrame{xTF},
alphaTrackingFrame{alpha},
positionTrackingFrame{posTF},
covarianceTrackingFrame{covTF} {
}

}
}
}
