// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CA/IndexTable.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_INDEXTABLE_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_INDEXTABLE_H_

#include <array>
#include <vector>

#include "ITSReconstruction/CA/Cluster.h"
#include "ITSReconstruction/CA/Constants.h"

namespace o2
{
namespace ITS
{
namespace CA
{

class IndexTable final
{
 public:
  IndexTable();
  IndexTable(const int, const std::vector<Cluster>&);

  const int getBin(const int) const;

  const std::vector<int> selectBins(const float, const float, const float, const float);

  static constexpr int getZBinIndex(const float, const float);
  static constexpr int getPhiBinIndex(const float);
  static constexpr int getBinIndex(const int, const int);

 private:
  int mLayerIndex;
  std::array<int, Constants::IndexTable::ZBins * Constants::IndexTable::PhiBins + 1> mTableBins;
};

inline const int IndexTable::getBin(const int binIndex) const { return mTableBins[binIndex]; }

constexpr int IndexTable::getZBinIndex(const float layerIndex, const float zCoordinate)
{
  return (zCoordinate + Constants::ITS::LayersZCoordinate[layerIndex]) *
         Constants::IndexTable::InverseZBinSize[layerIndex];
}

constexpr int IndexTable::getPhiBinIndex(const float currentPhi)
{
  return (currentPhi * Constants::IndexTable::InversePhiBinSize);
}

constexpr int IndexTable::getBinIndex(const int zIndex, const int phiIndex)
{
  return std::min(phiIndex * Constants::IndexTable::PhiBins + zIndex,
                  Constants::IndexTable::ZBins * Constants::IndexTable::PhiBins);
}

}
}
}

#endif /* O2_ITSMFT_RECONSTRUCTION_CA_INDEXTABLE_H_ */
