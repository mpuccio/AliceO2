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

// Gate 3 workflow-onboarding Slice 2: focused tests for the driver-level
// preflight in ITSCAWorkflow/ConfigPreflight.h -- the three checks
// its-ca-tracker-workflow.cxx runs, in order, before constructing any
// DataProcessorSpec/device.

#define BOOST_TEST_MODULE ITSMFT ITSCAConfigPreflight
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <stdexcept>

#include <boost/property_tree/ptree.hpp>
#include <fairlogger/Logger.h>

#include "CommonUtils/ConfigurableParam.h"
#include "ITSCAWorkflow/ConfigPreflight.h"
#include "ITSMFTTracking/TrackingConfigParam.h"

using namespace o2::its::ca;

namespace
{
struct FatalToExceptionFixture {
  FatalToExceptionFixture()
  {
    fair::Logger::OnFatal([]() { throw std::runtime_error("fatal"); });
  }
};
} // namespace

// --- applyConfigKeyValuesOrFatal(): preflight runs before the update -------

BOOST_FIXTURE_TEST_CASE(LegacyNamespaceIsRejectedBeforeAnyUpdate, FatalToExceptionFixture)
{
  // Sentinel: if the rejection did not actually run before
  // ConfigurableParam::updateFromString(), this legacy-namespace string
  // would still throw from updateFromString() itself (unknown param), so
  // this alone would not distinguish "preflight fired first" from "update
  // itself fatal'd" -- the meaningful assertion is in the next test, which
  // confirms the dedicated param was NOT mutated by the rejected string.
  BOOST_CHECK_THROW(applyConfigKeyValuesOrFatal("ITSCATrackerParam.trackFollowerTop=1"), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(RejectedStringNeverReachesConfigurableParamUpdate, FatalToExceptionFixture)
{
  // A malicious/confused string mixing a real dedicated-namespace override
  // with an offending legacy one must not have its dedicated part applied
  // either -- the whole string is rejected pre-update, atomically.
  o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "useDiamond", false);
  BOOST_CHECK_THROW(
    applyConfigKeyValuesOrFatal("ITSCommonCATrackerParam.useDiamond=true;ITSCATrackerParam.trackFollowerTop=1"),
    std::runtime_error);
  BOOST_CHECK_EQUAL(o2::itsmft::ITSCommonCATrackerParam::Instance().useDiamond, false);
}

BOOST_FIXTURE_TEST_CASE(DedicatedNamespaceIsAcceptedAndApplied, FatalToExceptionFixture)
{
  o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "useDiamond", false);
  BOOST_CHECK_NO_THROW(applyConfigKeyValuesOrFatal("ITSCommonCATrackerParam.useDiamond=true"));
  BOOST_CHECK_EQUAL(o2::itsmft::ITSCommonCATrackerParam::Instance().useDiamond, true);
  o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "useDiamond", false);
}

BOOST_FIXTURE_TEST_CASE(EmptyStringIsAcceptedAndApplied, FatalToExceptionFixture)
{
  BOOST_CHECK_NO_THROW(applyConfigKeyValuesOrFatal(""));
}

BOOST_FIXTURE_TEST_CASE(LegacyNamespaceWithoutFieldIsRejected, FatalToExceptionFixture)
{
  BOOST_CHECK_THROW(applyConfigKeyValuesOrFatal("ITSCATrackerParam=1"), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(MixedInputRejectsLegacyNamespaceInEitherPosition, FatalToExceptionFixture)
{
  for (const auto* config : {"ITSCATrackerParam.trackFollowerTop=1;ITSCommonCATrackerParam.useDiamond=true",
                             "ITSCommonCATrackerParam.useDiamond=true;ITSCATrackerParam.trackFollowerTop=1"}) {
    o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "useDiamond", false);
    BOOST_CHECK_THROW(applyConfigKeyValuesOrFatal(config), std::runtime_error);
    BOOST_CHECK_EQUAL(o2::itsmft::ITSCommonCATrackerParam::Instance().useDiamond, false);
  }
}

BOOST_FIXTURE_TEST_CASE(OuterWhitespaceAndInternalKeyWhitespaceKeepNamespace, FatalToExceptionFixture)
{
  o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "useDiamond", false);
  BOOST_CHECK_THROW(
    applyConfigKeyValuesOrFatal("  ITSCommonCATrackerParam.useDiamond=true ; ITSCATrackerParam.trackFollowerTop=1  "),
    std::runtime_error);
  BOOST_CHECK_EQUAL(o2::itsmft::ITSCommonCATrackerParam::Instance().useDiamond, false);

  BOOST_CHECK_THROW(applyConfigKeyValuesOrFatal("ITSCATrackerParam.trackFollowerTop = 1"), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(EmptyEntriesAndUnrelatedNamespacesAreAccepted, FatalToExceptionFixture)
{
  o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "dropTFUponFailure", false);
  o2::conf::ConfigurableParam::setValue<int>("ITSVertexerParam", "nIterations", 1);
  BOOST_CHECK_NO_THROW(applyConfigKeyValuesOrFatal(
    ";;ITSCommonCATrackerParam.dropTFUponFailure=true;;;ITSVertexerParam.nIterations=2;;"));
  BOOST_CHECK_EQUAL(o2::itsmft::ITSCommonCATrackerParam::Instance().dropTFUponFailure, true);
  BOOST_CHECK_EQUAL(o2::its::VertexerParamConfig::Instance().nIterations, 2);
  o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "dropTFUponFailure", false);
  o2::conf::ConfigurableParam::setValue<int>("ITSVertexerParam", "nIterations", 1);
}

BOOST_FIXTURE_TEST_CASE(MalformedTokensRemainConfiguratorErrors, FatalToExceptionFixture)
{
  for (const auto* config : {"ITSCATrackerParamNoEquals", "=ITSCATrackerParam.x", "ITSCATrackerParam.x="}) {
    BOOST_CHECK_THROW(applyConfigKeyValuesOrFatal(config), std::runtime_error);
  }
}

BOOST_FIXTURE_TEST_CASE(RepeatedAcceptedAndRejectedCallsRemainDeterministic, FatalToExceptionFixture)
{
  for (int i = 0; i < 5; ++i) {
    o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "useDiamond", false);
    BOOST_CHECK_THROW(applyConfigKeyValuesOrFatal("ITSCATrackerParam.trackFollowerTop=1"), std::runtime_error);
    BOOST_CHECK_NO_THROW(applyConfigKeyValuesOrFatal("ITSCommonCATrackerParam.useDiamond=true"));
    BOOST_CHECK_EQUAL(o2::itsmft::ITSCommonCATrackerParam::Instance().useDiamond, true);
  }
  o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "useDiamond", false);
}

// --- requireSyncTrackingModeOrFatal(): every non-Sync mode fails closed ----

BOOST_FIXTURE_TEST_CASE(SyncModeIsAccepted, FatalToExceptionFixture)
{
  BOOST_CHECK_NO_THROW(requireSyncTrackingModeOrFatal(o2::itsmft::TrackingMode::Sync));
}

BOOST_FIXTURE_TEST_CASE(EveryNonSyncModeFailsClosed, FatalToExceptionFixture)
{
  const std::array<o2::itsmft::TrackingMode::Type, 4> rejected{
    o2::itsmft::TrackingMode::Off, o2::itsmft::TrackingMode::Unset,
    o2::itsmft::TrackingMode::Async, o2::itsmft::TrackingMode::Cosmics};
  for (const auto mode : rejected) {
    BOOST_CHECK_THROW(requireSyncTrackingModeOrFatal(mode), std::runtime_error);
  }
}

// --- requireDiamondVertexConstraintOrFatal(): useDiamond is mandatory -----

BOOST_FIXTURE_TEST_CASE(MissingDiamondConstraintFailsClosed, FatalToExceptionFixture)
{
  o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "useDiamond", false);
  BOOST_CHECK_THROW(requireDiamondVertexConstraintOrFatal(), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(ConfiguredDiamondConstraintIsAccepted, FatalToExceptionFixture)
{
  o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "useDiamond", true);
  BOOST_CHECK_NO_THROW(requireDiamondVertexConstraintOrFatal());
  o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "useDiamond", false);
}
