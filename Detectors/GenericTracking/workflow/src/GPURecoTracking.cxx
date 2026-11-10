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

/// \file GPURecoTracking.cxx

#include "GenericTracking/Workflow/GPURecoTracking.h"

#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Workflow/TrackerSpec.h"

#include <stdexcept>
#include <string>

namespace o2::gtrk::workflow
{

std::unique_ptr<GPURecoTracking> makeGPURecoTracking(std::string_view detector, bool useMC)
{
  if (detector == "its") {
    return makeGPURecoTrackingOf<detectors::ITS>(useMC);
  }
  if (detector == "mft") {
    return makeGPURecoTrackingOf<detectors::MFT>(useMC);
  }
  throw std::invalid_argument("GenericTracking: no tracking inside the GPU reconstruction for " + std::string{detector});
}

} // namespace o2::gtrk::workflow
