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

#ifndef GTRK_WORKFLOW_ITSMFTOPTIONS_H_
#define GTRK_WORKFLOW_ITSMFTOPTIONS_H_

#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "DataFormatsITSMFT/DPLAlpideParamInitializer.h"
#include "Framework/ConfigContext.h"

#include <concepts>

namespace o2::gtrk::workflow
{

template <DetectorLayout Det>
[[nodiscard]] bool staggeredInput(const o2::framework::ConfigContext& context)
{
  if constexpr (std::same_as<Det, detectors::ITS>) {
    return o2::itsmft::DPLAlpideParamInitializer::isITSStaggeringEnabled(context);
  } else if constexpr (std::same_as<Det, detectors::MFT>) {
    return o2::itsmft::DPLAlpideParamInitializer::isMFTStaggeringEnabled(context);
  } else {
    return false;
  }
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_ITSMFTOPTIONS_H_
