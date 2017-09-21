// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Layer.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_LAYER_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_LAYER_H_

#include <vector>
#include <utility>

#include "ITSReconstruction/CA/Cluster.h"

namespace o2
{
namespace ITS
{
namespace CA
{

class Layer final
{
 public:
  Layer();
  Layer(const int);

  int getLayerIndex() const;
  const std::vector<Cluster>& getClusters() const;
  const Cluster& getCluster(int) const;
  int getClustersSize() const;

  const std::vector<TrackingFrameInfo>& getTrackingFrameInfo() const;
  const TrackingFrameInfo& getTrackingFrameInfo(int) const;

  void addCluster(const int clusterId, const float cluster_xyz[3]);
  template<typename ... T> void addTrackingFrameInfo(T&& ... values);
 private:
  int mLayerIndex;
  std::vector<Cluster> mClusters;
  std::vector<TrackingFrameInfo> mTrackingFrameInfo;
};

inline int Layer::getLayerIndex() const { return mLayerIndex; }

inline const std::vector<Cluster>& Layer::getClusters() const { return mClusters; }

inline const Cluster& Layer::getCluster(int clusterIndex) const { return mClusters[clusterIndex]; }

inline int Layer::getClustersSize() const { return mClusters.size(); }

inline const std::vector<TrackingFrameInfo>& Layer::getTrackingFrameInfo() const { return mTrackingFrameInfo; }

inline const TrackingFrameInfo& Layer::getTrackingFrameInfo(int index) const { return mTrackingFrameInfo[index]; }

template<typename ... T> void Layer::addTrackingFrameInfo(T&& ... values) {
  mTrackingFrameInfo.emplace_back(std::forward<T>(values)...);
}

}
}
}

#endif /* O2_ITSMFT_RECONSTRUCTION_CA_LAYER_H_ */
