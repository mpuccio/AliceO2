// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Tracker.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_TRACKER_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_TRACKER_H_

#include <iostream>
#include <vector>

#include "ITSReconstruction/CA/Cluster.h"
#include "ITSReconstruction/CA/Road.h"
#include "DetectorsBase/Track.h"

namespace o2
{
namespace ITS
{
namespace CA
{

class Event;
struct PrimaryVertexContext;

class Tracker final
{
 public:
  explicit Tracker(const Event&);

  Tracker(const Tracker&) = delete;
  Tracker& operator=(const Tracker&) = delete;

  std::vector<std::vector<Road>> clustersToTracks();

 private:
  void computeTracklets(PrimaryVertexContext&);
  void computeCells(PrimaryVertexContext&);
  void findCellsNeighbours(PrimaryVertexContext&);
  void findRoads(PrimaryVertexContext&);
  void findTracks(PrimaryVertexContext&);
  void traverseCellsTree(PrimaryVertexContext&, const int, const int);
  Base::Track::TrackParCov buildTrackSeed(const Cluster& c1, const Cluster& c2, const Cluster& c3, const TrackingFrameInfo& tf3);

  float            mBz;
  const Event&     mEvent;
};

}
}
}

#endif /* O2_ITSMFT_RECONSTRUCTION_CA_TRACKER_H_ */
