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
/// \file TrackingUtils.cxx
/// \brief
///

#include "ITSReconstruction/CA/TrackingUtils.h"

#include <cmath>

#include "ITSReconstruction/CA/Constants.h"
#include "ITSReconstruction/CA/IndexTableUtils.h"
#include "ITSReconstruction/CA/MathUtils.h"

namespace o2
{
namespace ITS
{
namespace CA
{

GPU_DEVICE const int4 TrackingUtils::getBinsRect(const Cluster& currentCluster, const int layerIndex,
    const float directionZIntersection)
{
  const float zRangeMin = directionZIntersection - 2 * Constants::Thresholds::ZCoordinateCut;
  const float phiRangeMin = currentCluster.phiCoordinate - Constants::Thresholds::PhiCoordinateCut;
  const float zRangeMax = directionZIntersection + 2 * Constants::Thresholds::ZCoordinateCut;
  const float phiRangeMax = currentCluster.phiCoordinate + Constants::Thresholds::PhiCoordinateCut;

  if (zRangeMax < -Constants::ITS::LayersZCoordinate()[layerIndex + 1]
      || zRangeMin > Constants::ITS::LayersZCoordinate()[layerIndex + 1] || zRangeMin > zRangeMax) {

    return getEmptyBinsRect();
  }

  return int4 { MATH_MAX(0, IndexTableUtils::getZBinIndex(layerIndex + 1, zRangeMin)),
      IndexTableUtils::getPhiBinIndex(MathUtils::getNormalizedPhiCoordinate(phiRangeMin)), MATH_MIN(
          Constants::IndexTable::ZBins - 1, IndexTableUtils::getZBinIndex(layerIndex + 1, zRangeMax)),
      IndexTableUtils::getPhiBinIndex(MathUtils::getNormalizedPhiCoordinate(phiRangeMax)) };
}

}
}
}
