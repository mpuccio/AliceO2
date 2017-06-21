// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAEvent.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CAEVENT_H_
#define TRACKINGITSU_INCLUDE_CAEVENT_H_

#include <array>
#include <vector>

#include "CAConstants.h"
#include "CALayer.h"

class CAEvent final
{

  public:
    explicit CAEvent(const int);

    int getEventId() const;
    const std::array<float, 3>& getPrimaryVertex(const int) const;
    const CALayer& getLayer(const int) const;
    int getPrimaryVerticesNum() const;

    void addPrimaryVertex(const float, const float, const float);
    void printPrimaryVertices() const;
    void pushClusterToLayer(const int, const int, const float, const float, const float, const float, const int);
    const int getTotalClusters() const;

  private:
    const int mEventId;
    std::vector<std::array<float, 3>> mPrimaryVertices;
    std::array<CALayer, CAConstants::ITS::LayersNumber> mLayers;
};

inline int CAEvent::getEventId() const
{
  return mEventId;
}

inline const std::array<float, 3>& CAEvent::getPrimaryVertex(const int vertexIndex) const
{
  return mPrimaryVertices[vertexIndex];
}

inline const CALayer& CAEvent::getLayer(const int layerIndex) const
{
  return mLayers[layerIndex];
}

inline int CAEvent::getPrimaryVerticesNum() const
{

  return mPrimaryVertices.size();
}

#endif /* TRACKINGITSU_INCLUDE_CAEVENT_H_ */
