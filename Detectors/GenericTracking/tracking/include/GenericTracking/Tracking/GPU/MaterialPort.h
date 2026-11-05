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

/// \file MaterialPort.h
/// \brief O2's material LUT ported to a device as a FlatObject.
///
/// Declared here, defined in host code (src/GPU/MaterialPort.cxx):
/// O2 compiles a FlatObject's host methods only outside GPU translation units.

#ifndef GTRK_TRACKING_GPU_MATERIALPORT_H_
#define GTRK_TRACKING_GPU_MATERIALPORT_H_

#include "GenericTracking/Tracking/GPU/DeviceStorage.h"

#include <cstddef>

namespace o2::base
{
class MatLayerCylSet;
}

namespace o2::gtrk
{

template <>
struct FlatObjectPort<o2::base::MatLayerCylSet> {
  static std::size_t bufferBytes(const o2::base::MatLayerCylSet& source);
  static FlatObjectImage image(const o2::base::MatLayerCylSet& source, char* futureBuffer);
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_MATERIALPORT_H_
