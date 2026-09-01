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

void requireSupportedTrackingModeOrFatal(o2::itsmft::TrackingMode::Type mode)
{
  if (mode != o2::itsmft::TrackingMode::Sync && mode != o2::itsmft::TrackingMode::Async) {
    LOGP(fatal,
         "ITS common-CA tracker workflow supports tracking-mode 'sync' and 'async'; '{}' is not supported",
         o2::itsmft::TrackingMode::toString(mode));
  }
}

void requireVertexConstraintOrFatal()
{
  const auto& tc = o2::itsmft::ITSCommonCATrackerParam::Instance();
  const bool useTruth = o2::its::VertexerParamConfig::Instance().useTruthSeeding;
  if (tc.useDiamond == useTruth) {
    LOGP(fatal,
         "ITS common-CA tracker workflow requires exactly one supported vertex source: set either "
         "ITSCommonCATrackerParam.useDiamond=true for a static constraint or "
         "ITSVertexerParam.useTruthSeeding=true for per-collision MC vertices");
  }
}

} // namespace o2::its::ca
