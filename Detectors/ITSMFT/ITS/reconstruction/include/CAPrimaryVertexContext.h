// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAPrimaryVertexContext.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CATRACKERCONTEXT_H_
#define TRACKINGITSU_INCLUDE_CATRACKERCONTEXT_H_

#include <array>
#include <vector>

#include "CACell.h"
#include "CAConstants.h"
#include "CAIndexTable.h"
#include "CARoad.h"
#include "CATracklet.h"

class CAEvent;

struct CAPrimaryVertexContext final
{
    explicit CAPrimaryVertexContext(const CAEvent&, const int);

    const int primaryVertexIndex;

    std::array<std::vector<CACluster>, CAConstants::ITS::LayersNumber> clusters;
    std::array<CAIndexTable, CAConstants::ITS::TrackletsPerRoad> indexTables;

    std::array<std::vector<CATracklet>, CAConstants::ITS::TrackletsPerRoad> tracklets;
    std::array<std::vector<int>, CAConstants::ITS::CellsPerRoad> trackletsLookupTable;

    std::array<std::vector<CACell>, CAConstants::ITS::CellsPerRoad> cells;
    std::array<std::vector<int>, CAConstants::ITS::CellsPerRoad - 1> cellsLookupTable;

    std::vector<CARoad> roads;
};

#endif /* TRACKINGITSU_INCLUDE_CATRACKERCONTEXT_H_ */
