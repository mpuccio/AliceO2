// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Cell.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "ITSReconstruction/CA/Cell.h"

namespace o2
{
namespace ITS
{
namespace CA
{

Cell::Cell(const int firstClusterIndex, const int secondClusterIndex, const int thirdClusterIndex,
               const int firstTrackletIndex, const int secondTrackletIndex,
               const std::array<float, 3>& normalVectorCoordinates, const float curvature)
  : mClusterIndices { firstClusterIndex, secondClusterIndex, thirdClusterIndex },
    mTrackletIndices { firstTrackletIndex, secondTrackletIndex },
    mNormalVectorCoordinates(normalVectorCoordinates),
    mCurvature{ curvature },
    mLevel{ 1 }
{
  // Nothing to do
}

void Cell::setLevel(const int level) { mLevel = level; }

bool Cell::combineCells(const Cell& otherCell, int otherCellId)
{
  if (mClusterIndices[1] == otherCell.getClusterIndex(2) &&
      mClusterIndices[0] == otherCell.getClusterIndex(1)) {
    mNeighbours.push_back(otherCellId);

    int otherCellLevel = otherCell.getLevel();

    if (otherCellLevel >= mLevel) {
      setLevel(otherCellLevel + 1);
    }

    return true;

  } else {
    return false;
  }
}

}
}
}
