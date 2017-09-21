// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Cell.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_CELL_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_CELL_H_

#include <array>
#include <vector>

namespace o2
{
namespace ITS
{
namespace CA
{

class Cell final
{
 public:
  Cell(const int, const int, const int, const int, const int, const std::array<float, 3>&, const float);

  int getClusterIndex(const int i) const;
  int getTrackletIndex(const int i) const;

  int getLevel() const;
  float getCurvature() const;
  int getNumberOfNeighbours() const;
  int getNeighbourCellId(const int) const;
  const std::array<float, 3>& getNormalVectorCoordinates() const;

  void setLevel(const int level);

  bool combineCells(const Cell&, int);

 private:
  const std::array<int,3> mClusterIndices;               ///< Cluster indices from inside outward
  const std::array<int,2> mTrackletIndices;              ///< Tracklet indices from inside outward
  const std::array<float, 3> mNormalVectorCoordinates;
  const float mCurvature;
  int mLevel;
  std::vector<int> mNeighbours;
};

inline int Cell::getClusterIndex(const int i) const { return mClusterIndices[i]; }

inline int Cell::getTrackletIndex(const int i) const { return mTrackletIndices[i]; }

inline int Cell::getLevel() const { return mLevel; }

inline float Cell::getCurvature() const { return mCurvature; }

inline int Cell::getNumberOfNeighbours() const { return mNeighbours.size(); }

inline int Cell::getNeighbourCellId(const int neighbourIndex) const { return mNeighbours[neighbourIndex]; }

inline const std::array<float, 3>& Cell::getNormalVectorCoordinates() const { return mNormalVectorCoordinates; }

}
}
}

#endif /* O2_ITSMFT_RECONSTRUCTION_CA_CELL_H_ */
