// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CARoad.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "CARoad.h"

CARoad::CARoad()
    : mCellIds { }, mRoadSize { }, mIsFakeRoad { }
{
  resetRoad();
}

CARoad::CARoad(int cellLayer, int cellId)
    : CARoad()
{
  addCell(cellLayer, cellId);
}

void CARoad::resetRoad()
{
  mCellIds.fill(CAConstants::ITS::UnusedIndex);
  mRoadSize = 0;
}

void CARoad::addCell(int cellLayer, int cellId)
{
  if (mCellIds[cellLayer] == CAConstants::ITS::UnusedIndex) {

    ++mRoadSize;
  }

  mCellIds[cellLayer] = cellId;
}
