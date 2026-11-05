// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "GenericTracking/Tracking/GPU/MaterialPort.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <cstddef>

namespace o2::gtrk
{

std::size_t FlatObjectPort<o2::base::MatLayerCylSet>::bufferBytes(const o2::base::MatLayerCylSet& source) { return flatObjectBufferBytes(source); }

FlatObjectImage FlatObjectPort<o2::base::MatLayerCylSet>::image(const o2::base::MatLayerCylSet& source, char* futureBuffer)
{
  return flatObjectImage(source, futureBuffer);
}

} // namespace o2::gtrk
