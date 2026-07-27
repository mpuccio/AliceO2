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

#include "ITSMFTTracking/ConfigKeyValuesPreflight.h"

#include "CommonUtils/StringUtils.h"

namespace o2::itsmft::tracking
{

namespace
{
// Matches o2::its::TrackerParamConfig's O2ParamDef name
// (ITStracking/TrackingConfigParam.h) exactly -- the frozen legacy ITS
// tracking parameter namespace this preflight rejects for the opt-in
// common-CA workflow.
constexpr std::string_view kLegacyITSNamespace = "ITSCATrackerParam";
} // namespace

ConfigKeyValuesPreflightResult checkITSCommonCAConfigKeyValues(std::string_view configKeyValues)
{
  ConfigKeyValuesPreflightResult result;

  // Same grammar as o2::conf::ConfigurableParam::updateFromString(): split on
  // ';', trim each token, drop empty tokens (Str::tokenize's default
  // skipEmpty=true, matching updateFromString's own call).
  const auto tokens = o2::utils::Str::tokenize(std::string(configKeyValues), ';', true);

  for (auto const& token : tokens) {
    const auto eq = token.find('=');
    if (eq == std::string::npos || eq == 0 || eq == token.size() - 1) {
      // Malformed key=value token: not this helper's concern.
      continue;
    }
    const auto key = token.substr(0, eq);
    const auto dot = key.find('.');
    const auto ns = dot == std::string::npos ? key : key.substr(0, dot);
    if (ns == kLegacyITSNamespace) {
      result.outcome = ConfigKeyValuesPreflightOutcome::RejectedLegacyITSNamespace;
      result.offendingToken = token;
      result.offendingNamespace = ns;
      return result;
    }
  }

  return result;
}

} // namespace o2::itsmft::tracking
