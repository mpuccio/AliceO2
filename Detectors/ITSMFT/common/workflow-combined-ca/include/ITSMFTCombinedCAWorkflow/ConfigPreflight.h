// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file ConfigPreflight.h
/// \brief Driver-level configuration preflight for the combined
///        ITS+MFT common-CA tracker workflow.
///
/// This implementation remains independent of the single-detector workflows;
/// its fatal conditions match the shared Sync/useDiamond policy where relevant.
///
/// Called, in order, from defineDataProcessing() -- before any
/// DataProcessorSpec/device is constructed -- so a rejected configuration
/// never reaches a running device:
///
///   1. requireCombinedTrackingEnabledOrFatal()
///   2. requireSyncTrackingModeOrFatal(trackingMode)
///   3. requireDiamondVertexConstraintOrFatal()
///   4. requireNoMFTIRFrameConfigOrFatal()
///
/// The MFT continuous-readout GRP fact is CCDB-sourced and checked at runtime.

#ifndef ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_CONFIGPREFLIGHT_H_
#define ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_CONFIGPREFLIGHT_H_

#include "ITSMFTTracking/Configuration.h"

namespace o2::itsmft::combined
{

/// Fatals unless ITSMFTCombinedCATrackerParam::enabled is set.
void requireCombinedTrackingEnabledOrFatal();

/// Fatals unless mode == Sync, naming the rejected mode explicitly.
void requireSyncTrackingModeOrFatal(o2::itsmft::TrackingMode::Type mode);

/// Fatals unless ITSCommonCATrackerParam::useDiamond is set.
void requireDiamondVertexConstraintOrFatal();

/// Fatals if MFTTrackingParam::irFramesOnly is set; combined tracking has no
/// IR-frame input or mask.
void requireNoMFTIRFrameConfigOrFatal();

/// Fatals unless the fetched MFT GRP fact reports continuous readout. This is
/// called before the first TimeFrame is processed.
void requireContinuousMFTReadoutOrFatal(bool continuousReadout);

} // namespace o2::itsmft::combined

#endif // ALICEO2_ITSMFT_COMBINEDCAWORKFLOW_CONFIGPREFLIGHT_H_
