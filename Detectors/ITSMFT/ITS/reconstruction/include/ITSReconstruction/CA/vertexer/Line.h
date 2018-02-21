// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CA/vertexer/Line.h
/// \brief Definition of the ITS CA line

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_LINE_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_LINE_H_

#include <array>

namespace o2
{
namespace ITS
{
namespace CA
{
  
class Line
{
public:
  Line( std::array<float, 3> firstPoint, std::array<float, 3> secondPoint );
  virtual ~Line();
  Line& operator=(const Line&) = delete;

protected:
  std::array<float, 3> mOriginPoint;
  std::array<float, 3> mCosinesDirector;
};

}
}
}
#endif /* O2_ITSMFT_RECONSTRUCTION_CA_LINE_H_ */