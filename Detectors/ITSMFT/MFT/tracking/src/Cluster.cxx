// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file Cluster.cxx
///

#include "MFTTracking/Cluster.h"
#include "MFTTracking/IndexTableUtils.h"

#include "MathUtils/Utils.h"
#include "MathUtils/Cartesian2D.h"

namespace o2
{
namespace mft
{

Cluster::Cluster(const Float_t x, const Float_t y, const Float_t z, const Int_t index)
  : xCoordinate{x},
    yCoordinate{y},
    zCoordinate{z},
    phi{0.},
    rCoordinate{0.},
    clusterId{index},
    indexTableBin{0}
{
  auto clsPoint2D = Point2D<Float_t>(x, y);
  rCoordinate = clsPoint2D.R();
  phi = clsPoint2D.Phi();
  o2::utils::BringTo02PiGen(phi);
}

Cluster::Cluster(const Int_t layerIndex, const Cluster& other)
  : xCoordinate{other.xCoordinate},
    yCoordinate{other.yCoordinate},
    zCoordinate{other.zCoordinate},
    phi{0.},
    rCoordinate{0.},
    clusterId{other.clusterId},
    indexTableBin{index_table_utils::getBinIndex(index_table_utils::getRBinIndex(layerIndex, rCoordinate),
                                                 index_table_utils::getPhiBinIndex(phi))}
{
  auto clsPoint2D = Point2D<Float_t>(other.xCoordinate, other.yCoordinate);
  rCoordinate = clsPoint2D.R();
  phi = clsPoint2D.Phi();
  o2::utils::BringTo02PiGen(phi);
}

} // namespace mft
} // namespace o2
