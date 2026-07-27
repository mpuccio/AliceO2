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

// Workflow-onboarding Slice 1: focused tests for the pure
// checkITSCommonCAConfigKeyValues() namespace preflight
// (ConfigKeyValuesPreflight.h). This helper is deliberately additive and
// unwired: no workflow spec exists yet, and it never itself logs or aborts
// -- these tests only exercise the pure string-parsing contract. Linked only
// against the lightweight O2::ITSMFTTrackingParams library (not the full
// O2::ITSMFTTracking core), demonstrating the helper does not need it.

#define BOOST_TEST_MODULE ITSMFT ITSCommonCAConfigKeyValuesPreflight
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/ConfigKeyValuesPreflight.h"

using namespace o2::itsmft::tracking;

BOOST_AUTO_TEST_CASE(EmptyStringIsAccepted)
{
  const auto result = checkITSCommonCAConfigKeyValues("");
  BOOST_CHECK(result.outcome == ConfigKeyValuesPreflightOutcome::Accepted);
  BOOST_CHECK(result.offendingToken.empty());
  BOOST_CHECK(result.offendingNamespace.empty());
}

BOOST_AUTO_TEST_CASE(DedicatedNamespaceIsAccepted)
{
  const auto result = checkITSCommonCAConfigKeyValues("ITSCommonCATrackerParam.dropTFUponFailure=true");
  BOOST_CHECK(result.outcome == ConfigKeyValuesPreflightOutcome::Accepted);
}

BOOST_AUTO_TEST_CASE(UnrelatedNamespacesAreAccepted)
{
  for (auto const* s : {"MatLUT.x=1", "HBFUtils.nHBFPerTF=32", "ITSCommonCATrackerParam.useDiamond=1;MatLUT.x=1",
                        "SomeFlagWithNoDot=1"}) {
    const auto result = checkITSCommonCAConfigKeyValues(s);
    BOOST_CHECK_MESSAGE(result.outcome == ConfigKeyValuesPreflightOutcome::Accepted, "unexpected rejection for: " << s);
  }
}

BOOST_AUTO_TEST_CASE(LegacyNamespaceIsRejected)
{
  const auto result = checkITSCommonCAConfigKeyValues("ITSCATrackerParam.trackFollowerTop=1");
  BOOST_CHECK(result.outcome == ConfigKeyValuesPreflightOutcome::RejectedLegacyITSNamespace);
  BOOST_CHECK_EQUAL(result.offendingNamespace, "ITSCATrackerParam");
  BOOST_CHECK_EQUAL(result.offendingToken, "ITSCATrackerParam.trackFollowerTop=1");
}

BOOST_AUTO_TEST_CASE(LegacyNamespaceWithNoFieldStillHasThatNamespace)
{
  // "ITSCATrackerParam" alone (no '.') is its own namespace too -- still
  // must be rejected, since it is exactly the legacy name.
  const auto result = checkITSCommonCAConfigKeyValues("ITSCATrackerParam=1");
  BOOST_CHECK(result.outcome == ConfigKeyValuesPreflightOutcome::RejectedLegacyITSNamespace);
  BOOST_CHECK_EQUAL(result.offendingNamespace, "ITSCATrackerParam");
}

BOOST_AUTO_TEST_CASE(MixedInputRejectsWhenLegacyNamespacePresentAnywhere)
{
  // Order matters only for which offending token is reported, not whether
  // rejection happens: a dedicated-namespace token earlier in the string
  // must not mask a legacy-namespace token later in the same string.
  const auto result = checkITSCommonCAConfigKeyValues(
    "ITSCommonCATrackerParam.dropTFUponFailure=true;ITSCATrackerParam.trackFollowerTop=1;MatLUT.x=1");
  BOOST_CHECK(result.outcome == ConfigKeyValuesPreflightOutcome::RejectedLegacyITSNamespace);
  BOOST_CHECK_EQUAL(result.offendingNamespace, "ITSCATrackerParam");
  BOOST_CHECK_EQUAL(result.offendingToken, "ITSCATrackerParam.trackFollowerTop=1");
}

BOOST_AUTO_TEST_CASE(WhitespaceBearingInputIsHandledLikeConfigurableParam)
{
  // Str::tokenize(..., ';', true) trims each token around ';' boundaries,
  // matching o2::conf::ConfigurableParam::updateFromString()'s own call.
  const auto result = checkITSCommonCAConfigKeyValues("  ITSCommonCATrackerParam.useDiamond=1 ; ITSCATrackerParam.trackFollowerTop=1  ");
  BOOST_CHECK(result.outcome == ConfigKeyValuesPreflightOutcome::RejectedLegacyITSNamespace);
  BOOST_CHECK_EQUAL(result.offendingToken, "ITSCATrackerParam.trackFollowerTop=1");
}

BOOST_AUTO_TEST_CASE(InternalWhitespaceAroundEqualsStillIdentifiesNamespace)
{
  // Internal whitespace between the key and '=' is not stripped by
  // tokenize's outer trim; the namespace (before the first '.') is
  // unaffected by it either way.
  const auto result = checkITSCommonCAConfigKeyValues("ITSCATrackerParam.trackFollowerTop = 1");
  BOOST_CHECK(result.outcome == ConfigKeyValuesPreflightOutcome::RejectedLegacyITSNamespace);
  BOOST_CHECK_EQUAL(result.offendingNamespace, "ITSCATrackerParam");
}

BOOST_AUTO_TEST_CASE(EmptyEntriesBetweenSemicolonsAreSkipped)
{
  const auto result = checkITSCommonCAConfigKeyValues(";;ITSCommonCATrackerParam.useDiamond=1;;;MatLUT.x=1;;");
  BOOST_CHECK(result.outcome == ConfigKeyValuesPreflightOutcome::Accepted);
}

BOOST_AUTO_TEST_CASE(MultipleEntriesAllDedicatedOrUnrelatedAreAccepted)
{
  const auto result = checkITSCommonCAConfigKeyValues(
    "ITSCommonCATrackerParam.dropTFUponFailure=true;ITSCommonCATrackerParam.useDiamond=1;MatLUT.x=1;HBFUtils.y=2");
  BOOST_CHECK(result.outcome == ConfigKeyValuesPreflightOutcome::Accepted);
}

BOOST_AUTO_TEST_CASE(MalformedTokensAreNotThisHelpersConcern)
{
  // No '=', or '=' at the very start/end: skipped, not flagged as a legacy
  // rejection -- o2::conf::ConfigurableParam::updateFromString() rejects
  // these itself once the string is actually applied.
  for (auto const* s : {"ITSCATrackerParamNoEquals", "=ITSCATrackerParam.x", "ITSCATrackerParam.x="}) {
    const auto result = checkITSCommonCAConfigKeyValues(s);
    BOOST_CHECK_MESSAGE(result.outcome == ConfigKeyValuesPreflightOutcome::Accepted, "unexpected rejection for: " << s);
  }
}

BOOST_AUTO_TEST_CASE(RepeatedCallsAreDeterministic)
{
  for (int i = 0; i < 5; ++i) {
    BOOST_CHECK(checkITSCommonCAConfigKeyValues("ITSCATrackerParam.x=1").outcome ==
                ConfigKeyValuesPreflightOutcome::RejectedLegacyITSNamespace);
    BOOST_CHECK(checkITSCommonCAConfigKeyValues("ITSCommonCATrackerParam.x=1").outcome ==
                ConfigKeyValuesPreflightOutcome::Accepted);
  }
}
