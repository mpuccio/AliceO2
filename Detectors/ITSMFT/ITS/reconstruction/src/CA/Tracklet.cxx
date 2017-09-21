// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Tracklet.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "ITSReconstruction/CA/Tracklet.h"

o2::ITS::CA::Tracklet::Tracklet(const int firstClusterIndex, const int secondClusterIndex,
                                const float tanLambda, const float phiCoordinate)
  : firstClusterIndex{ firstClusterIndex },
    secondClusterIndex{ secondClusterIndex },
    tanLambda{ tanLambda },
    phiCoordinate{ phiCoordinate }
{
  // Nothing to do
}
