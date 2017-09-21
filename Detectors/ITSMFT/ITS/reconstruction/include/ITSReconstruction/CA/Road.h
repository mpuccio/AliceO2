// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Road.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_ROAD_
#define O2_ITSMFT_RECONSTRUCTION_CA_ROAD_

#include <array>

#include "ITSReconstruction/CA/Constants.h"

namespace o2
{
namespace ITS
{
namespace CA
{

class Road final
{
 public:
  Road();
  Road(int, int);

  int getRoadSize() const;
  int& operator[](const int&);

  void resetRoad();
  void addCell(int, int);

 private:
  std::array<int, Constants::ITS::CellsPerRoad> mCellIds;
  int mRoadSize;
};

inline int Road::getRoadSize() const { return mRoadSize; }

inline int& Road::operator[](const int& i) { return mCellIds[i]; }

}
}
}

#endif /* O2_ITSMFT_RECONSTRUCTION_CA_ROAD */
