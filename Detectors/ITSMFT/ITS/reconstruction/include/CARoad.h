// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CARoad.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CAROAD_H_
#define TRACKINGITSU_INCLUDE_CAROAD_H_

#include <array>

#include "CAConstants.h"

class CARoad final
{
  public:
    CARoad();
    CARoad(int, int);

    int getRoadSize() const;
    int getLabel() const;
    void setLabel(const int);
    bool isFakeRoad() const;
    void setFakeRoad(const bool);
    int &operator[](const int&);

    void resetRoad();
    void addCell(int, int);

  private:
    std::array<int, CAConstants::ITS::CellsPerRoad> mCellIds;
    int mRoadSize;
    int mLabel;
    bool mIsFakeRoad;
};

inline int CARoad::getRoadSize() const
{

  return mRoadSize;
}

inline int CARoad::getLabel() const
{

  return mLabel;
}

inline void CARoad::setLabel(const int label)
{

  mLabel = label;
}

inline int& CARoad::operator [](const int& i)
{
  return mCellIds[i];
}

inline bool CARoad::isFakeRoad() const
{
  return mIsFakeRoad;
}

inline void CARoad::setFakeRoad(const bool isFakeRoad)
{
  mIsFakeRoad = isFakeRoad;
}

#endif /* TRACKINGITSU_INCLUDE_CAROAD_H_ */
