// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Cluster.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_CLUSTER_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_CLUSTER_H_

#include <array>

namespace o2
{
namespace ITS
{
namespace CA
{

struct Cluster final {
  ///TODO: Make interfaces uniform
  Cluster(const int id, const float x, const float y, const float z, const int layer);
  Cluster(const int layer, const std::array<float, 3>& vertex, const Cluster& );

  ///TODO: use Point3D
  int clusterId;
  float xCoordinate;
  float yCoordinate;
  float zCoordinate;
  float phiCoordinate;
  float rCoordinate;
  int indexTableBinIndex;
};

struct TrackingFrameInfo {
  TrackingFrameInfo(float xTF, float alpha, std::array<float,2>&& posTF, std::array<float,3>&& covTF);

  float xTrackingFrame;
  float alphaTrackingFrame;
  std::array<float,2> positionTrackingFrame;
  std::array<float,3> covarianceTrackingFrame;
};

}
}
}
#endif /* O2_ITSMFT_RECONSTRUCTION_CA_CLUSTER_H_ */
