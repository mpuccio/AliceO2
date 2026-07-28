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

// ITSCommonCATrackerParam.nThreads: focused tests for the common ITS
// interface's own dedicated thread-count knob (TrackingInterface.cxx's
// initialiseTracker()). No catalog/decoder/geometry fixture is needed here
// -- ITSMFTTrackingInterface<7>::initialise() (resolveTrackingParameters() +
// initialiseMemoryPool() + initialiseTracker()) touches neither, and the ITS
// Sync branch of TrackingMode::getTrackingParameters() (unlike MFT's) never
// dereferences o2::base::Propagator::Instance(), so the bare 3-arg
// constructor is enough to reach a real, fully-constructed TrackerTraits<7>
// and its tbb::task_arena.

#define BOOST_TEST_MODULE ITSMFT ITSCommonCANThreads
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <fairlogger/Logger.h>
#include <stdexcept>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITSMFTTracking/TrackingInterface.h"
#include "ITStracking/TrackingConfigParam.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{
// ConfigurableParamHelper<P>::Instance() deliberately returns const P& to
// discourage casual mutation from production code; tests that need to stage
// a specific singleton value use the same const_cast pattern as e.g.
// testTrackingInterfaceLoadFailureContract.cxx's mutableMFTAlpideParam().
ITSCommonCATrackerParam& mutableITSCommonCATrackerParam()
{
  return const_cast<ITSCommonCATrackerParam&>(ITSCommonCATrackerParam::Instance());
}

o2::its::TrackerParamConfig& mutableLegacyITSTrackerParamConfig()
{
  return const_cast<o2::its::TrackerParamConfig&>(o2::its::TrackerParamConfig::Instance());
}

// Restores every singleton field this test file stages, on scope exit
// (including a thrown/failed check), so no test case here can leak state
// into another test case in this binary.
struct ScopedParams {
  ScopedParams()
    : originalUseDiamond(mutableITSCommonCATrackerParam().useDiamond),
      originalNThreads(mutableITSCommonCATrackerParam().nThreads),
      originalLegacyNThreads(mutableLegacyITSTrackerParamConfig().nThreads)
  {
    // ITS Sync requires useDiamond=true or TrackingMode::getTrackingParameters()
    // fatals (workflow-onboarding Slice 1); every test case in this file
    // needs a successful resolveTrackingParameters() to reach initialiseTracker().
    mutableITSCommonCATrackerParam().useDiamond = true;
  }
  ~ScopedParams()
  {
    mutableITSCommonCATrackerParam().useDiamond = originalUseDiamond;
    mutableITSCommonCATrackerParam().nThreads = originalNThreads;
    mutableLegacyITSTrackerParamConfig().nThreads = originalLegacyNThreads;
  }
  bool originalUseDiamond;
  int originalNThreads;
  int originalLegacyNThreads;
};

// LOGP(fatal, ...) normally terminates the process (FairLogger default). A
// process-local OnFatal handler converts it into a catchable exception, the
// same pattern testITSCommonCATrackingModeConfiguration.cxx already uses.
// Each ITSMFT test source file builds its own executable (o2_add_test ==
// one binary per SOURCES file), so this handler cannot leak into unrelated
// test binaries.
struct FatalToExceptionFixture {
  FatalToExceptionFixture() { fair::Logger::OnFatal([]() { throw std::runtime_error("fatal"); }); }
};
} // namespace

BOOST_FIXTURE_TEST_CASE(DefaultNThreadsIsConsumedByTheITSInterface, ScopedParams)
{
  BOOST_REQUIRE_EQUAL(mutableITSCommonCATrackerParam().nThreads, 1);

  ITSMFTTrackingInterface<ITSNLayers> interface{false, o2::itsmft::TrackingMode::Sync, false};
  interface.initialise();
  BOOST_REQUIRE(interface.isActive());
  BOOST_CHECK_EQUAL(interface.getTrackerNThreads(), 1);
}

BOOST_FIXTURE_TEST_CASE(OverriddenNThreadsIsConsumedNotTheLegacyNamespace, ScopedParams)
{
  // Deliberately different values: if the interface mistakenly consumed
  // TrackerParamRef<ITS>::get() (the frozen legacy o2::its::TrackerParamConfig,
  // "ITSCATrackerParam") instead of its own ITSCommonCATrackerParam, the
  // observed thread count below would be 7, not 4.
  mutableITSCommonCATrackerParam().nThreads = 4;
  mutableLegacyITSTrackerParamConfig().nThreads = 7;

  ITSMFTTrackingInterface<ITSNLayers> interface{false, o2::itsmft::TrackingMode::Sync, false};
  interface.initialise();
  BOOST_REQUIRE(interface.isActive());
  BOOST_CHECK_EQUAL(interface.getTrackerNThreads(), 4);
}

BOOST_FIXTURE_TEST_CASE(NonPositiveNThreadsFailsSafelyBeforeTracking, FatalToExceptionFixture)
{
  ScopedParams guard;
  for (const int invalid : {0, -1, -8}) {
    mutableITSCommonCATrackerParam().nThreads = invalid;
    ITSMFTTrackingInterface<ITSNLayers> interface{false, o2::itsmft::TrackingMode::Sync, false};
    BOOST_CHECK_THROW(interface.initialise(), std::runtime_error);
  }
}
