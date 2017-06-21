// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CACell.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CACELL_H_
#define TRACKINGITSU_INCLUDE_CACELL_H_

#include <array>
#include <vector>

class CACell final
{
  public:
    CACell(const int, const int, const int, const int, const int, const std::array<float, 3>&, const float);

    int getFirstClusterIndex() const;
    int getSecondClusterIndex() const;
    int getThirdClusterIndex() const;
    int getFirstTrackletIndex() const;
    int getSecondTrackletIndex() const;
    int getLevel() const;
    float getCurvature() const;
    int getNumberOfNeighbours() const;
    int getNeighbourCellId(const int) const;
    const std::array<float, 3>& getNormalVectorCoordinates() const;

    void setLevel(const int level);

    bool combineCells(const CACell&, int);

  private:
    const int mFirstClusterIndex;
    const int mSecondClusterIndex;
    const int mThirdClusterIndex;
    const int mFirstTrackletIndex;
    const int mSecondTrackletIndex;
    const std::array<float, 3> mNormalVectorCoordinates;
    const float mCurvature;
    int mLevel;
    std::vector<int> mNeighbours;
};

inline int CACell::getFirstClusterIndex() const
{
  return mFirstClusterIndex;
}

inline int CACell::getSecondClusterIndex() const
{
  return mSecondClusterIndex;
}

inline int CACell::getThirdClusterIndex() const
{
  return mThirdClusterIndex;
}

inline int CACell::getFirstTrackletIndex() const
{
  return mFirstTrackletIndex;
}

inline int CACell::getSecondTrackletIndex() const
{
  return mSecondTrackletIndex;
}

inline int CACell::getLevel() const
{
  return mLevel;
}

inline float CACell::getCurvature() const
{
  return mCurvature;
}

inline int CACell::getNumberOfNeighbours() const
{
  return mNeighbours.size();
}

inline int CACell::getNeighbourCellId(const int neighbourIndex) const
{
  return mNeighbours[neighbourIndex];
}

inline const std::array<float, 3>& CACell::getNormalVectorCoordinates() const
{
  return mNormalVectorCoordinates;
}

#endif /* TRACKINGITSU_INCLUDE_CACELL_H_ */
