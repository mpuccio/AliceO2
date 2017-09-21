// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MathUtils.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_MATHUTILS_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_MATHUTILS_H_

#include <array>

#include "DetectorsBase/Constants.h"

namespace o2
{
namespace ITS
{
namespace CA
{

namespace MathUtils
{

float calculatePhiCoordinate(const float, const float);
float calculateRCoordinate(const float, const float);
float computeCurvature(float x1, float y1, float x2, float y2, float x3, float y3);
float computeCurvatureCentreX(float x1, float y1, float x2, float y2, float x3, float y3);
float computeTanDipAngle(float x1, float y1, float x2, float y2, float z1, float z2);

constexpr float getNormalizedPhiCoordinate(const float);
constexpr std::array<float, 3> crossProduct(const std::array<float, 3>&, const std::array<float, 3>&);

}

constexpr float MathUtils::getNormalizedPhiCoordinate(const float phiCoordinate)
{
  return (phiCoordinate < 0)
           ? phiCoordinate + Base::Constants::k2PI
           : (phiCoordinate > Base::Constants::k2PI) ? phiCoordinate - Base::Constants::k2PI : phiCoordinate;
}

constexpr std::array<float, 3> MathUtils::crossProduct(const std::array<float, 3>& firstVector,
                                                         const std::array<float, 3>& secondVector)
{
  return std::array<float, 3>{ { (firstVector[1] * secondVector[2]) - (firstVector[2] * secondVector[1]),
                                 (firstVector[2] * secondVector[0]) - (firstVector[0] * secondVector[2]),
                                 (firstVector[0] * secondVector[1]) - (firstVector[1] * secondVector[0]) } };
}

}
}
}

#endif /* O2_ITSMFT_RECONSTRUCTION_CA_MATHUTILS_H_ */
