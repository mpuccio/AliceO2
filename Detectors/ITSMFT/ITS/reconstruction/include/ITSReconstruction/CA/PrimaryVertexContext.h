// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file PrimaryVertexContext.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_PRIMARYVERTEXCONTEXT_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_PRIMARYVERTEXCONTEXT_H_

#include <array>
#include <vector>

#include "ITSReconstruction/CA/Cell.h"
#include "ITSReconstruction/CA/Constants.h"
#include "ITSReconstruction/CA/IndexTable.h"
#include "ITSReconstruction/CA/Road.h"
#include "ITSReconstruction/CA/Tracklet.h"
#include "ITSReconstruction/CA/Track.h"

namespace o2
{
namespace ITS
{
namespace CA
{

class Event;

struct PrimaryVertexContext final {
  explicit PrimaryVertexContext(const Event&, const int);

  const int primaryVertexIndex;

  std::array<std::vector<Cluster>, Constants::ITS::LayersNumber> clusters;
  std::array<IndexTable, Constants::ITS::TrackletsPerRoad> indexTables;
  std::array<std::vector<bool>, Constants::ITS::LayersNumber> usedClusters;

  std::array<std::vector<Tracklet>, Constants::ITS::TrackletsPerRoad> tracklets;
  std::array<std::vector<int>, Constants::ITS::CellsPerRoad> trackletsLookupTable;

  std::array<std::vector<Cell>, Constants::ITS::CellsPerRoad> cells;
  std::array<std::vector<int>, Constants::ITS::CellsPerRoad - 1> cellsLookupTable;

  std::vector<Road> roads;
  std::vector<Track> tracks;
};

}
}
}

#endif /* O2_ITSMFT_RECONSTRUCTION_CA_PRIMARYVERTEXCONTEXT_H_ */
