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

/// @file   GPUWorkflowGenericTracking.cxx
/// @brief  GenericTracking's ITS and MFT tracking inside the GPU reconstruction, after the TPC's, in its pool

#include "GPUWorkflow/GPUWorkflowSpec.h"
#include "GPUO2Interface.h"
#include "GenericTracking/Workflow/GPURecoTracking.h"

namespace o2::gpu
{

void GPURecoWorkflowSpec::createGenericTracking()
{
  if (mSpecConfig.runGenericTrackingITS) {
    mGenericTracking.push_back(o2::gtrk::workflow::makeGPURecoTracking("its"));
  }
  if (mSpecConfig.runGenericTrackingMFT) {
    mGenericTracking.push_back(o2::gtrk::workflow::makeGPURecoTracking("mft"));
  }
}

void GPURecoWorkflowSpec::initFunctionGenericTracking(o2::framework::InitContext& ic)
{
  for (auto& tracking : mGenericTracking) {
    tracking->init(ic, *mGPUReco->GetReconstruction());
  }
}

void GPURecoWorkflowSpec::runGenericTracking(o2::framework::ProcessingContext& pc, bool reconstructed)
{
  for (auto& tracking : mGenericTracking) {
    tracking->run(pc, reconstructed);
  }
}

} // namespace o2::gpu
