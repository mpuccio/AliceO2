// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Line.cxx
/// \brief
/// \author matteo.concas@cern.ch

#include <cmath>
#include "ITSReconstruction/CA/vertexer/Line.h"


namespace o2
{
namespace ITS
{
namespace CA
{

Line::Line( std::array<float, 3> firstPoint, std::array<float, 3> secondPoint ) :
mOriginPoint{ firstPoint }
{
  for ( int index { 0 }; index < 3; ++index) mCosinesDirector[index] = secondPoint[index] - firstPoint[index];
  float inverseNorm { 1.f/std::sqrt( mCosinesDirector[0]*mCosinesDirector[0] 
    + mCosinesDirector[1]*mCosinesDirector[1] + mCosinesDirector[2]*mCosinesDirector[2] ) };
  for ( int index { 0 }; index < 3; ++index) mCosinesDirector[index] *= inverseNorm; 
}

Line::~Line() {};
}
}
}