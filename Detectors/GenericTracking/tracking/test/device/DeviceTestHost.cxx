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

/// \file DeviceTestHost.cxx
/// \brief The host code of the test-kernel library: what O2 compiles only
///        outside GPU translation units.

#include "DeviceTestKernels.h"

#include "GenericTracking/Tracking/GPU/DeviceStorage.h"

namespace o2::gtrk
{

std::size_t FlatObjectPort<ToyFlatObject>::bufferBytes(const ToyFlatObject& source) { return flatObjectBufferBytes(source); }

FlatObjectImage FlatObjectPort<ToyFlatObject>::image(const ToyFlatObject& source, char* futureBuffer) { return flatObjectImage(source, futureBuffer); }

} // namespace o2::gtrk
