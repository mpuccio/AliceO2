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

/// \file DeviceTrackerITS.cu
/// \brief The backend's device tracker for the ITS.
///        One translation unit per layout,
///        so that the layouts compile in parallel.

#include "DeviceTracker.cuh"

#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Tracking/GPU/MaterialPort.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <cstddef>
#include <memory>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
{

std::unique_ptr<DeviceTracker<detectors::ITS, o2::base::MatLayerCylSet>> makeITSTracker(DeviceRuntime& runtime, DeviceBlock block)
{
  return std::make_unique<DeviceTrackerImpl<detectors::ITS, o2::base::MatLayerCylSet>>(runtime, block);
}

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
