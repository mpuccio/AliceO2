// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file Event.h
/// \brief 
///

#ifndef TRACKINGITSU_INCLUDE_EVENT_H_
#define TRACKINGITSU_INCLUDE_EVENT_H_

#include <array>
#include <vector>

#include "ITSReconstruction/CA/Constants.h"
#include "ITSReconstruction/CA/Layer.h"

namespace o2
{
namespace ITS
{
namespace CA
{

class Event
  final
  {
    public:
      explicit Event(const int);

      int getEventId() const;
      const float3& getPrimaryVertex(const int) const;
      const Layer& getLayer(const int) const;
      int getPrimaryVerticesNum() const;
      void addPrimaryVertex(const float, const float, const float);
      void printPrimaryVertices() const;
      void pushClusterToLayer(const int, const int, const float, const float, const float, const float, const int);
      int getTotalClusters() const;

    private:
      const int mEventId;
      std::vector<float3> mPrimaryVertices;
      std::array<Layer, Constants::ITS::LayersNumber> mLayers;
  };

  inline int Event::getEventId() const
  {
    return mEventId;
  }

  inline const float3& Event::getPrimaryVertex(const int vertexIndex) const
  {
    return mPrimaryVertices[vertexIndex];
  }

  inline const Layer& Event::getLayer(const int layerIndex) const
  {
    return mLayers[layerIndex];
  }

  inline int Event::getPrimaryVerticesNum() const
  {

    return mPrimaryVertices.size();
  }

}
}
}

#endif /* TRACKINGITSU_INCLUDE_EVENT_H_ */
