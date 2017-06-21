// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CACell.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "CACell.h"

CACell::CACell(const int firstClusterIndex, const int secondClusterIndex, const int thirdClusterIndex,
    const int firstTrackletIndex, const int secondTrackletIndex, const std::array<float, 3>& normalVectorCoordinates,
    const float curvature)
    : mFirstClusterIndex { firstClusterIndex }, mSecondClusterIndex { secondClusterIndex }, mThirdClusterIndex {
        thirdClusterIndex }, mFirstTrackletIndex(firstTrackletIndex), mSecondTrackletIndex(secondTrackletIndex), mNormalVectorCoordinates(
        normalVectorCoordinates), mCurvature { curvature }, mLevel { 1 }
{
  // Nothing to do
}

void CACell::setLevel(const int level)
{
  mLevel = level;
}

bool CACell::combineCells(const CACell& otherCell, int otherCellId)
{
  if (mSecondClusterIndex == otherCell.getThirdClusterIndex()
      && mFirstClusterIndex == otherCell.getSecondClusterIndex()) {

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
