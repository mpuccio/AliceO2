// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 4 C4: focused tests for the driver-level preflight in
// ITSMFTCombinedCAWorkflow/ConfigPreflight.h -- the four checks
// itsmft-combined-ca-tracker-workflow.cxx runs, in order, before
// constructing any DataProcessorSpec/device, plus the fifth (runtime,
// CCDB-sourced) check called from CombinedCATrackerSpec.cxx.

#define BOOST_TEST_MODULE ITSMFT ITSMFTCombinedCAConfigPreflight
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <stdexcept>

#include <boost/property_tree/ptree.hpp>
#include <fairlogger/Logger.h>

#include "CommonUtils/ConfigurableParam.h"
#include "ITSMFTCombinedCAWorkflow/ConfigPreflight.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "MFTTracking/MFTTrackingParam.h"

using namespace o2::itsmft::combined;

namespace
{
struct FatalToExceptionFixture {
  FatalToExceptionFixture()
  {
    fair::Logger::OnFatal([]() { throw std::runtime_error("fatal"); });
  }
};
} // namespace

// --- requireCombinedTrackingEnabledOrFatal(): default-false opt-in --------

BOOST_FIXTURE_TEST_CASE(DisabledByDefaultFailsClosed, FatalToExceptionFixture)
{
  o2::conf::ConfigurableParam::setValue<bool>("ITSMFTCombinedCATrackerParam", "enabled", false);
  BOOST_CHECK_THROW(requireCombinedTrackingEnabledOrFatal(), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(ExplicitlyEnabledIsAccepted, FatalToExceptionFixture)
{
  o2::conf::ConfigurableParam::setValue<bool>("ITSMFTCombinedCATrackerParam", "enabled", true);
  BOOST_CHECK_NO_THROW(requireCombinedTrackingEnabledOrFatal());
  o2::conf::ConfigurableParam::setValue<bool>("ITSMFTCombinedCATrackerParam", "enabled", false);
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

// --- requireNoMFTIRFrameConfigOrFatal(): irFramesOnly is rejected ---------

BOOST_FIXTURE_TEST_CASE(IRFramesOnlyConfigFailsClosed, FatalToExceptionFixture)
{
  o2::conf::ConfigurableParam::setValue<bool>("MFTTracking", "irFramesOnly", true);
  BOOST_CHECK_THROW(requireNoMFTIRFrameConfigOrFatal(), std::runtime_error);
  o2::conf::ConfigurableParam::setValue<bool>("MFTTracking", "irFramesOnly", false);
}

BOOST_FIXTURE_TEST_CASE(DefaultIRFramesConfigIsAccepted, FatalToExceptionFixture)
{
  o2::conf::ConfigurableParam::setValue<bool>("MFTTracking", "irFramesOnly", false);
  BOOST_CHECK_NO_THROW(requireNoMFTIRFrameConfigOrFatal());
}

// --- requireContinuousMFTReadoutOrFatal(): triggered MFT readout is rejected

BOOST_FIXTURE_TEST_CASE(ContinuousReadoutIsAccepted, FatalToExceptionFixture)
{
  BOOST_CHECK_NO_THROW(requireContinuousMFTReadoutOrFatal(true));
}

BOOST_FIXTURE_TEST_CASE(TriggeredReadoutFailsClosed, FatalToExceptionFixture)
{
  BOOST_CHECK_THROW(requireContinuousMFTReadoutOrFatal(false), std::runtime_error);
}
