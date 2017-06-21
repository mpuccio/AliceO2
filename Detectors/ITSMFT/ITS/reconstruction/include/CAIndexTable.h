// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAIndexTable.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CALOOKUPTABLE_H_
#define TRACKINGITSU_INCLUDE_CALOOKUPTABLE_H_

#include <array>
#include <vector>

#include <CACluster.h>
#include <CAConstants.h>

class CAIndexTable final
{
  public:
    CAIndexTable();
    CAIndexTable(const int, const std::vector<CACluster>&);

    const int getBin(const int) const;

    const std::vector<int> selectBins(const float, const float, const float, const float);

  private:
    int mLayerIndex;
    std::array<int, CAConstants::IndexTable::ZBins * CAConstants::IndexTable::PhiBins + 1> mTableBins;
};

inline const int CAIndexTable::getBin(const int binIndex) const
{

  return mTableBins[binIndex];
}

#endif /* TRACKINGITSU_INCLUDE_CALOOKUPTABLE_H_ */
