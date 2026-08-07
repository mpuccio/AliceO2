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
///
/// \file ConfigKeyValuesPreflight.h
/// \brief Pure --configKeyValues namespace preflight for the opt-in ITS
///        common-CA workflow.

#ifndef ALICEO2_ITSMFT_TRACKING_CONFIG_KEY_VALUES_PREFLIGHT_H_
#define ALICEO2_ITSMFT_TRACKING_CONFIG_KEY_VALUES_PREFLIGHT_H_

#include <string>
#include <string_view>

namespace o2::itsmft::tracking
{

enum class ConfigKeyValuesPreflightOutcome {
  Accepted,                  // no legacy-ITS-namespace override present
  RejectedLegacyITSNamespace // at least one "ITSCATrackerParam.*" override present
};

struct ConfigKeyValuesPreflightResult {
  ConfigKeyValuesPreflightOutcome outcome = ConfigKeyValuesPreflightOutcome::Accepted;
  std::string offendingToken;     // full "key=value" token that triggered rejection (empty if Accepted)
  std::string offendingNamespace; // namespace portion of the offending key (empty if Accepted)
};

/// Pure, DPL-free preflight. Tokens are semicolon-separated, trimmed, and
/// split at the first '='; malformed tokens are left to the configurator.
/// The first key in namespace `ITSCATrackerParam` is reported. The helper
/// neither logs nor aborts; the caller decides how to report rejection.
ConfigKeyValuesPreflightResult checkITSCommonCAConfigKeyValues(std::string_view configKeyValues);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_CONFIG_KEY_VALUES_PREFLIGHT_H_ */
