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
///        common-CA workflow (workflow-onboarding Slice 1).

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

/// Pure, DPL-free namespace preflight for a raw --configKeyValues string.
///
/// Mirrors the exact tokenization grammar
/// o2::conf::ConfigurableParam::updateFromString() uses (semicolon-separated
/// via o2::utils::Str::tokenize(..., ';', true) -- each token trimmed, empty
/// tokens skipped -- then each token's key/value split on the first '='), so
/// a string this helper accepts is guaranteed to be tokenized identically
/// once actually applied.
///
/// Flags the first token whose key's namespace (the substring before the
/// first '.', or the whole key if there is no '.') is exactly the legacy ITS
/// namespace "ITSCATrackerParam" -- the name registered by the frozen
/// legacy o2::its::TrackerParamConfig in O2::ITStracking. The dedicated
/// ITSCommonCATrackerParam namespace and every other namespace (unrelated
/// globally registered parameters such as MatLUT or HBFUtils) pass through
/// unchanged. A malformed token (no '=', or '=' at the start/end) is not
/// this helper's concern -- it is left for
/// o2::conf::ConfigurableParam::updateFromString() to reject when the
/// string is actually applied -- so it is skipped here without effect.
///
/// Never itself logs or aborts: the caller (a future workflow spec) is
/// responsible for turning a rejection into its own fatal diagnostic.
ConfigKeyValuesPreflightResult checkITSCommonCAConfigKeyValues(std::string_view configKeyValues);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_CONFIG_KEY_VALUES_PREFLIGHT_H_ */
