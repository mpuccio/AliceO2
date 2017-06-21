// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CALabel.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CALABEL_H_
#define TRACKINGITSU_INCLUDE_CALABEL_H_

#include <ostream>

struct CALabel final
{
    CALabel(const int, const float, const float, const float, const int, const int);

    int monteCarloId;
    float transverseMomentum;
    float phiCoordinate;
    float pseudorapidity;
    int pdgCode;
    int numberOfClusters;

    friend std::ostream& operator<<(std::ostream&, const CALabel&);
};

#endif /* TRACKINGITSU_INCLUDE_CALABEL_H_ */
