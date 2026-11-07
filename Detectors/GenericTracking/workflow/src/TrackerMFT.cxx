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

/// \file TrackerMFT.cxx
/// \brief The MFT's tracking device (TrackerDevice.h): one unit per detector,
///        each instantiating the tracker, and its tracking inside the GPU
///        reconstruction (GPURecoTracking.h).

#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Workflow/Conditions.h"
#include "TrackerDevice.h"

namespace o2::gtrk
{

template o2::framework::DataProcessorSpec workflow::getTrackerSpec<detectors::MFT>(bool withVertices);
template std::unique_ptr<workflow::GPURecoTracking> workflow::makeGPURecoTrackingOf<detectors::MFT>();

} // namespace o2::gtrk
