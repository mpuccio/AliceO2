// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTCombinedCAWorkflow/ConfigPreflight.h"

#include "Framework/Logger.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "MFTTracking/MFTTrackingParam.h"

namespace o2::itsmft::combined
{

void requireCombinedTrackingEnabledOrFatal()
{
  if (!o2::itsmft::ITSMFTCombinedCATrackerParam::Instance().enabled) {
    LOGP(fatal,
         "o2-itsmft-combined-ca-tracker-workflow requires ITSMFTCombinedCATrackerParam.enabled=true; "
         "refusing to run a combined ITS+MFT tracking pass without an explicit opt-in");
  }
}

void requireSyncTrackingModeOrFatal(o2::itsmft::TrackingMode::Type mode)
{
  if (mode != o2::itsmft::TrackingMode::Sync) {
    LOGP(fatal,
         "Combined ITS+MFT common-CA tracker workflow supports only tracking-mode 'sync'; '{}' is not supported",
         o2::itsmft::TrackingMode::toString(mode));
  }
}

void requireDiamondVertexConstraintOrFatal()
{
  const auto& tc = o2::itsmft::ITSCommonCATrackerParam::Instance();
  if (!tc.useDiamond) {
    LOGP(fatal,
         "Combined ITS+MFT common-CA tracker workflow requires ITSCommonCATrackerParam.useDiamond=true: the "
         "combined coordinator's ITS leg has no real per-event vertexing capability yet, so any other "
         "vertex/beam constraint mode would silently run tracklet/cell finding against an always-empty "
         "per-ROF primary vertex table instead of failing loudly. Set "
         "ITSCommonCATrackerParam.useDiamond=true (and optionally .diamondPos[0..2]/.pvRes)");
  }
}

void requireNoMFTIRFrameConfigOrFatal()
{
  if (o2::mft::MFTTrackingParam::Instance().irFramesOnly) {
    LOGP(fatal,
         "Combined ITS+MFT common-CA tracker workflow does not support MFTTrackingParam.irFramesOnly=true: "
         "the combined ITS+MFT common-CA tracking flow has no IR-frame parameter and always enables every MFT "
         "ROF unmasked, so an IR-frame-filtered configuration would silently run unmasked triggered MFT "
         "data instead of failing loudly. Run the single-detector o2-mft-ca-tracker-workflow for "
         "IR-frame/triggered MFT data instead");
  }
}

void requireContinuousMFTReadoutOrFatal(bool continuousReadout)
{
  if (!continuousReadout) {
    LOGP(fatal,
         "Combined ITS+MFT common-CA tracker workflow does not support triggered MFT readout: GRPECS "
         "reports MFT is not in continuous readout for this run, but the combined ITS+MFT common-CA tracking flow "
         "always enables every MFT ROF unmasked (no IR-frame/trigger-window filtering). Run the "
         "single-detector o2-mft-ca-tracker-workflow for triggered MFT data instead");
  }
}

} // namespace o2::itsmft::combined
