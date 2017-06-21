// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAMathUtils.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CAUTILS_H_
#define TRACKINGITSU_INCLUDE_CAUTILS_H_

#include <array>

#include "CAConstants.h"

namespace CAMathUtils {
float calculatePhiCoordinate(const float, const float);
float calculateRCoordinate(const float, const float);
constexpr float getNormalizedPhiCoordinate(const float);

constexpr std::array<float, 3> crossProduct(const std::array<float, 3>&, const std::array<float, 3>&);
}

constexpr float CAMathUtils::getNormalizedPhiCoordinate(const float phiCoordinate)
{
  return (phiCoordinate < 0) ? phiCoordinate + CAConstants::Math::TwoPi :
         (phiCoordinate > CAConstants::Math::TwoPi) ? phiCoordinate - CAConstants::Math::TwoPi : phiCoordinate;
}

constexpr std::array<float, 3> CAMathUtils::crossProduct(const std::array<float, 3>& firstVector,
    const std::array<float, 3>& secondVector)
{

  return std::array<float, 3> { { (firstVector[1] * secondVector[2]) - (firstVector[2] * secondVector[1]),
      (firstVector[2] * secondVector[0]) - (firstVector[0] * secondVector[2]), (firstVector[0] * secondVector[1])
          - (firstVector[1] * secondVector[0]) } };
}

#endif /* TRACKINGITSU_INCLUDE_CAUTILS_H_ */
