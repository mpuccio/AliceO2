// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAMathUtils.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "CAMathUtils.h"

#include <cmath>

float CAMathUtils::calculatePhiCoordinate(const float xCoordinate, const float yCoordinate)
{
  return std::atan2(-yCoordinate, -xCoordinate) + CAConstants::Math::Pi;
}

float CAMathUtils::calculateRCoordinate(const float xCoordinate, const float yCoordinate)
{
  return std::sqrt(xCoordinate * xCoordinate + yCoordinate * yCoordinate);
}
