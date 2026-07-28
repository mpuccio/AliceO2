// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "ITSCAWorkflow/ConfigPreflight.h"

#include "CommonUtils/ConfigurableParam.h"
#include "Framework/Logger.h"
#include "ITSMFTTracking/ConfigKeyValuesPreflight.h"
#include "ITSMFTTracking/TrackingConfigParam.h"

namespace o2::its::ca
{

void applyConfigKeyValuesOrFatal(const std::string& configKeyValues)
{
  const auto preflight = o2::itsmft::tracking::checkITSCommonCAConfigKeyValues(configKeyValues);
  if (preflight.outcome == o2::itsmft::tracking::ConfigKeyValuesPreflightOutcome::RejectedLegacyITSNamespace) {
    LOGP(fatal,
        "ITS common-CA tracker workflow rejects legacy '{}' --configKeyValues override ('{}'); "
        "use the dedicated 'ITSCommonCATrackerParam' namespace instead",
        preflight.offendingNamespace, preflight.offendingToken);
  }
  o2::conf::ConfigurableParam::updateFromString(configKeyValues);
}

void requireSyncTrackingModeOrFatal(o2::itsmft::TrackingMode::Type mode)
{
  if (mode != o2::itsmft::TrackingMode::Sync) {
    LOGP(fatal,
        "ITS common-CA tracker workflow supports only tracking-mode 'sync'; '{}' is not supported yet",
        o2::itsmft::TrackingMode::toString(mode));
  }
}

void requireDiamondVertexConstraintOrFatal()
{
  const auto& tc = o2::itsmft::ITSCommonCATrackerParam::Instance();
  if (!tc.useDiamond) {
    LOGP(fatal,
        "ITS common-CA tracker workflow requires ITSCommonCATrackerParam.useDiamond=true: the common "
        "tracker has no real per-event vertexing capability for ITS yet, so any other vertex/beam "
        "constraint mode (the legacy real-vertexer default, or a bare CCDB MeanVertex beam-position "
        "override) would silently run tracklet/cell finding against an always-empty per-ROF primary "
        "vertex table instead of failing loudly. Set ITSCommonCATrackerParam.useDiamond=true (and "
        "optionally .diamondPos[0..2]/.pvRes) to select the one Sync vertex/beam-constraint mode this "
        "workflow can reproduce faithfully");
  }
}

} // namespace o2::its::ca
