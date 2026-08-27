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

#include "ITSCAWorkflow/ConfigPreflight.h"

#include <string_view>

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/StringUtils.h"
#include "Framework/Logger.h"
#include "ITSMFTTracking/TrackingConfigParam.h"

namespace o2::its::ca
{

namespace
{
// This is an ITS common-CA workflow policy, not a generic tracking-parameter
// validation. Keep the legacy namespace spelling local to the workflow that
// rejects it, so the shared parameter library does not expose a workflow API.
constexpr std::string_view kLegacyITSNamespace = "ITSCATrackerParam";
} // namespace

void applyConfigKeyValuesOrFatal(const std::string& configKeyValues)
{
  // Mirror ConfigurableParam::updateFromString()'s tokenization: split on
  // ';', trim each token, skip empty tokens, and split at the first '='.
  // Malformed tokens remain the configurator's responsibility.
  const auto tokens = o2::utils::Str::tokenize(configKeyValues, ';', true);
  for (const auto& token : tokens) {
    const auto eq = token.find('=');
    if (eq == std::string::npos || eq == 0 || eq == token.size() - 1) {
      continue;
    }
    const auto key = token.substr(0, eq);
    const auto dot = key.find('.');
    const auto ns = dot == std::string::npos ? key : key.substr(0, dot);
    if (ns == kLegacyITSNamespace) {
      LOGP(fatal,
           "ITS common-CA tracker workflow rejects legacy '{}' --configKeyValues override ('{}'); "
           "use the dedicated 'ITSCommonCATrackerParam' namespace instead",
           ns, token);
    }
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
