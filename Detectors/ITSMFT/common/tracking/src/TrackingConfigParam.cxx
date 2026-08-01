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

#include "ITSMFTTracking/TrackingConfigParam.h"

namespace o2::itsmft
{
// Ensure MFT CA params are registered in the global param database.
// ITS production tracking uses o2::its::TrackerParamConfig in O2::ITStracking.
static auto& sMFTCATrackerParam = TrackerParamConfig<o2::detectors::DetID::MFT>::Instance();
} // namespace o2::itsmft

// Dedicated ITS common-CA opt-in configuration (workflow-onboarding Slice 1).
// Registered explicitly here, via the standard non-template O2ParamImpl
// macro -- does not touch, instantiate, or rename the still-dormant
// templated TrackerParamConfig<DetID::ITS> above.
O2ParamImpl(o2::itsmft::ITSCommonCATrackerParam);
O2ParamImpl(o2::itsmft::ITSMFTCommonCAOutputParam);
