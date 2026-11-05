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

#include "GenericTracking/Workflow/RecoParam.h"

namespace o2::gtrk
{
// Each detector's instance, touched here so that it is in the parameter
// database and --configKeyValues finds its keys.
static auto& sRecoParamITS = RecoParam<o2::detectors::DetID::ITS>::Instance();
static auto& sRecoParamMFT = RecoParam<o2::detectors::DetID::MFT>::Instance();
#ifdef ENABLE_UPGRADES
static auto& sRecoParamTRK = RecoParam<o2::detectors::DetID::TRK>::Instance();
#endif
} // namespace o2::gtrk
