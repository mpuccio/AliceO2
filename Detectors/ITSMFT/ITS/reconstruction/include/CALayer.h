// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CALayer.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CALAYER_H_
#define TRACKINGITSU_INCLUDE_CALAYER_H_

#include <vector>

#include "CACluster.h"

class CALayer final
{
  public:
    CALayer();
    CALayer(const int);

    int getLayerIndex() const;
    const std::vector<CACluster>& getClusters() const;
    const CACluster& getCluster(int) const;
    int getClustersSize() const;

    void addCluster(const int, const float, const float, const float, const float, const int);

  private:
    int mLayerIndex;
    std::vector<CACluster> mClusters;
};

inline int CALayer::getLayerIndex() const
{
  return mLayerIndex;
}

inline const std::vector<CACluster>& CALayer::getClusters() const
{
  return mClusters;
}

inline const CACluster& CALayer::getCluster(int clusterIndex) const
{
  return mClusters[clusterIndex];
}

inline int CALayer::getClustersSize() const
{
  return mClusters.size();
}

#endif /* TRACKINGITSU_INCLUDE_CALAYER_H_ */
