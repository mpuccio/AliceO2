// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Event.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_EVENT_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_EVENT_H_

#include <array>
#include <vector>

#include <TClonesArray.h>

#include "ITSBase/GeometryTGeo.h"
#include "ITSReconstruction/CA/Constants.h"
#include "ITSReconstruction/CA/Layer.h"

namespace o2
{
namespace ITS
{
namespace CA
{
 
class Event final
{
 public:
  explicit Event(const int id, const GeometryTGeo* geo);

  void loadClusters(const TClonesArray& clusters);

  int getEventId() const;
  const std::array<float, 3>& getPrimaryVertex(const int index) const;
  const Layer& getLayer(const int index) const;
  int getPrimaryVerticesNum() const;

  void addPrimaryVertex(const float x, const float y, const float z);
  void printPrimaryVertices() const;

  const int getTotalClusters(const int numberOfLayers = Constants::ITS::LayersNumber) const;

  const TrackingFrameInfo& getTrackingFrameInfo(int layer, int index) const;
  void setGeometry(GeometryTGeo* geom);

 private:
  const int mEventId;
  const GeometryTGeo* mGeom; /// interface to geometry
  std::vector<std::array<float, 3>> mPrimaryVertices;
  std::array<Layer, Constants::ITS::LayersNumber> mLayers;
};

inline int Event::getEventId() const { return mEventId; }

inline const std::array<float, 3>& Event::getPrimaryVertex(const int vertexIndex) const
{
  return mPrimaryVertices[vertexIndex];
}

inline const Layer& Event::getLayer(const int layerIndex) const { return mLayers[layerIndex]; }

inline int Event::getPrimaryVerticesNum() const { return mPrimaryVertices.size(); }

inline const TrackingFrameInfo& Event::getTrackingFrameInfo(int layer, int index) const { return mLayers[layer].getTrackingFrameInfo(index); }

}
}
}

#endif /* O2_ITSMFT_RECONSTRUCTION_CA_EVENT_H_ */
