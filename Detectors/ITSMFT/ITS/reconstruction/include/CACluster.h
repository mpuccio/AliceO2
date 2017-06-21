// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CACluster.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CACLUSTER_H_
#define TRACKINGITSU_INCLUDE_CACLUSTER_H_

#include <array>

struct CACluster final
{
    CACluster(const int, const float, const float, const float, const float, const int);
    CACluster(const int, const std::array<float, 3>&, const CACluster&);

    int clusterId;
    float xCoordinate;
    float yCoordinate;
    float zCoordinate;
    float alphaAngle;
    int monteCarloId;
    float phiCoordinate;
    float rCoordinate;
    int indexTableBinIndex;
};

#endif /* TRACKINGITSU_INCLUDE_CACLUSTER_H_ */
