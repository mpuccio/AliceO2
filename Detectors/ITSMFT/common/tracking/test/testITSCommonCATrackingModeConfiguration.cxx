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

// Workflow-onboarding Slice 1: focused tests for the dedicated
// ITSCommonCATrackerParam configuration type (TrackingConfigParam.h) and the
// real ITS Sync branch of TrackingMode::getTrackingParameters()
// (Configuration.cxx). No workflow spec exists yet -- these tests call the
// common-tracking library directly, the same way
// testTrackingInterfaceLoadFailureContract.cxx already documents that the
// ITS branch of getTrackingParameters() used to unconditionally
// LOGP(fatal, ...) regardless of mode; that fatal is now real per-mode
// behaviour instead, exercised here.

#define BOOST_TEST_MODULE ITSMFT ITSCommonCATrackingModeConfiguration
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/property_tree/ptree.hpp>
#include <fairlogger/Logger.h>
#include <stdexcept>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/TrackingConfigParam.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

// --- Dedicated name is distinct from every other registered CA param name --

BOOST_AUTO_TEST_CASE(DedicatedNameIsDistinctFromLegacyAndMFTNames)
{
  const auto& itsCommonCA = ITSCommonCATrackerParam::Instance();
  const auto& itsLegacy = o2::its::TrackerParamConfig::Instance();
  const auto& mftCommonCA = TrackerParamConfig<o2::detectors::DetID::MFT>::Instance();

  BOOST_CHECK_EQUAL(itsCommonCA.getName(), "ITSCommonCATrackerParam");
  BOOST_CHECK_EQUAL(itsLegacy.getName(), "ITSCATrackerParam");
  BOOST_CHECK_EQUAL(mftCommonCA.getName(), "MFTCATrackerParam");

  BOOST_CHECK(itsCommonCA.getName() != itsLegacy.getName());
  BOOST_CHECK(itsCommonCA.getName() != mftCommonCA.getName());
  BOOST_CHECK(itsLegacy.getName() != mftCommonCA.getName());
}

// --- ITSCommonCATrackerParam defaults match the documented Sync baseline ---

BOOST_AUTO_TEST_CASE(DedicatedDefaultsMatchDocumentedSyncBaseline)
{
  const auto& tc = ITSCommonCATrackerParam::Instance();
  BOOST_CHECK_EQUAL(tc.dropTFUponFailure, false);
  BOOST_CHECK_EQUAL(tc.printMemory, false);
  BOOST_CHECK_EQUAL(tc.maxMemory, std::numeric_limits<size_t>::max());
  BOOST_CHECK_EQUAL(tc.saveTimeBenchmarks, false);
  BOOST_CHECK_EQUAL(tc.useDiamond, false);
  BOOST_CHECK_EQUAL(tc.diamondPos[0], 0.f);
  BOOST_CHECK_EQUAL(tc.diamondPos[1], 0.f);
  BOOST_CHECK_EQUAL(tc.diamondPos[2], 0.f);
  BOOST_CHECK_EQUAL(tc.pvRes, -1.f);
  BOOST_CHECK_EQUAL(tc.nThreads, 1);
}

// --- ITS Sync construction is valid, one-iteration, with expected values ---

BOOST_AUTO_TEST_CASE(ITSSyncTrackingParametersAreValidOneIteration)
{
  const auto trackParams = TrackingMode::getTrackingParameters(o2::detectors::DetID::ITS, TrackingMode::Sync);

  BOOST_REQUIRE_EQUAL(trackParams.size(), 1u);
  const auto& p = trackParams[0];

  BOOST_CHECK_EQUAL(p.NLayers, tracking::ITSNLayers);
  BOOST_CHECK_EQUAL(p.MinTrackLength, TrackerParamConfig<o2::detectors::DetID::ITS>::MinTrackLength);
  BOOST_CHECK_EQUAL(p.MinPt.size(), static_cast<size_t>(TrackerParamConfig<o2::detectors::DetID::ITS>::MaxTrackLength -
                                                        TrackerParamConfig<o2::detectors::DetID::ITS>::MinTrackLength + 1));
  BOOST_CHECK_EQUAL(p.StartLayerMask.count(), tracking::ITSNLayers); // default mask: all 7 barrel layers active

  // Administrative fields wired straight from the dedicated config's defaults.
  BOOST_CHECK_EQUAL(p.DropTFUponFailure, false);
  BOOST_CHECK_EQUAL(p.PrintMemory, false);
  BOOST_CHECK_EQUAL(p.MaxMemory, std::numeric_limits<size_t>::max());
  BOOST_CHECK_EQUAL(p.SaveTimeBenchmarks, false);
  BOOST_CHECK_EQUAL(p.UseDiamond, false);

  // resetDetectorDefaults(..., DetID::ITS) supplies real barrel geometry
  // defaults (TrackingParameters' own struct defaults); confirm they were
  // not clobbered.
  BOOST_CHECK_EQUAL(p.LayerRadii.size(), static_cast<size_t>(tracking::ITSNLayers));
  BOOST_CHECK_EQUAL(p.LayerZ.size(), static_cast<size_t>(tracking::ITSNLayers));
}

BOOST_AUTO_TEST_CASE(ITSSyncTrackingParametersAreDeterministic)
{
  const auto a = TrackingMode::getTrackingParameters(o2::detectors::DetID::ITS, TrackingMode::Sync);
  const auto b = TrackingMode::getTrackingParameters(o2::detectors::DetID::ITS, TrackingMode::Sync);
  BOOST_REQUIRE_EQUAL(a.size(), b.size());
  BOOST_CHECK_EQUAL(a[0].MinTrackLength, b[0].MinTrackLength);
  BOOST_CHECK_EQUAL(a[0].NLayers, b[0].NLayers);
  BOOST_CHECK(a[0].LayerRadii == b[0].LayerRadii);
}

// --- Every unsupported TrackingMode fails closed, none silently mapped -----
//
// LOGP(fatal, ...) normally terminates the process (FairLogger default). A
// process-local OnFatal handler converts it into a catchable exception so
// this remains a normal, non-crashing ctest case. Each ITSMFT test source
// file builds its own executable (o2_add_test == one binary per SOURCES
// file), so this handler cannot leak into unrelated test binaries.

namespace
{
struct FatalToExceptionFixture {
  FatalToExceptionFixture()
  {
    fair::Logger::OnFatal([]() { throw std::runtime_error("fatal"); });
  }
};
} // namespace

BOOST_FIXTURE_TEST_CASE(EveryUnsupportedITSModeFailsClosed, FatalToExceptionFixture)
{
  const std::array<TrackingMode::Type, 4> unsupported{
    TrackingMode::Off, TrackingMode::Unset, TrackingMode::Async, TrackingMode::Cosmics};

  for (const auto mode : unsupported) {
    BOOST_CHECK_THROW(TrackingMode::getTrackingParameters(o2::detectors::DetID::ITS, mode), std::runtime_error);
  }
}

BOOST_FIXTURE_TEST_CASE(SyncStillSucceedsAfterFatalHandlerInstalled, FatalToExceptionFixture)
{
  // The OnFatal fixture above must not turn the *supported* path into a
  // false failure: Sync should still construct normally.
  BOOST_CHECK_NO_THROW(TrackingMode::getTrackingParameters(o2::detectors::DetID::ITS, TrackingMode::Sync));
}

// --- MFT configuration behaviour: no regression -----------------------------
//
// This slice only adds an early-return ITS branch ahead of the pre-existing
// MFT body in getTrackingParameters(); every line reachable for
// detId==MFT (the "if (detId != detectors::DetID::MFT)" check onward) is
// byte-for-byte unchanged by this slice -- confirmed by inspection, not just
// by this test.
//
// Sync/Async/Cosmics are deliberately NOT exercised here for MFT: their
// shared tail (after the per-mode branches) calls
// o2::base::Propagator::Instance(), which -- with no GRP/magnetic-field
// singleton set up, as in a bare focused unit test -- segfaults inside
// PropagatorImpl::updateField() (Detectors/Base/src/Propagator.cxx). This is
// a pre-existing property of the untouched MFT tail, reproducible on the
// pre-Slice-1 code exactly the same way, unrelated to and unaffected by this
// slice; giving it real GRP/field setup is out of scope here. TrackingMode
// ::Off is the one mode whose early return (before the Propagator-dependent
// tail) is safe to call bare, and is what the MFT no-regression check below
// uses instead.

BOOST_FIXTURE_TEST_CASE(MFTOffStillReturnsEmptyNotFatal, FatalToExceptionFixture)
{
  // MFT's own Off semantics (legitimate "tracking disabled", not a failure)
  // are unrelated to and unaffected by ITS's stricter Slice 1 fail-closed
  // Off handling.
  BOOST_CHECK_NO_THROW({
    const auto trackParams = TrackingMode::getTrackingParameters(o2::detectors::DetID::MFT, TrackingMode::Off);
    BOOST_CHECK(trackParams.empty());
  });
}

// --- workflow-onboarding Slice 2: diamondPos/pvRes are wired through -------
//
// Mutates the global ITSCommonCATrackerParam singleton via
// ConfigurableParam::setValue -- deliberately placed last in this
// translation unit so no other test observes the mutated state.

BOOST_AUTO_TEST_CASE(DiamondPosAndPVresAreWiredIntoITSSyncTrackingParameters)
{
  o2::conf::ConfigurableParam::setValue<float>("ITSCommonCATrackerParam", "diamondPos[0]", 1.5f);
  o2::conf::ConfigurableParam::setValue<float>("ITSCommonCATrackerParam", "diamondPos[1]", -2.5f);
  o2::conf::ConfigurableParam::setValue<float>("ITSCommonCATrackerParam", "diamondPos[2]", 3.5f);
  o2::conf::ConfigurableParam::setValue<float>("ITSCommonCATrackerParam", "pvRes", 0.25f);
  o2::conf::ConfigurableParam::setValue<bool>("ITSCommonCATrackerParam", "useDiamond", true);

  const auto trackParams = TrackingMode::getTrackingParameters(o2::detectors::DetID::ITS, TrackingMode::Sync);
  BOOST_REQUIRE_EQUAL(trackParams.size(), 1u);
  const auto& p = trackParams[0];

  BOOST_CHECK_EQUAL(p.UseDiamond, true);
  BOOST_CHECK_EQUAL(p.Diamond[0], 1.5f);
  BOOST_CHECK_EQUAL(p.Diamond[1], -2.5f);
  BOOST_CHECK_EQUAL(p.Diamond[2], 3.5f);
  BOOST_CHECK_EQUAL(p.PVres, 0.25f);
}

// The shared ITS/MFT common-CA output route is explicitly opt-in.  Keep this
// mutation last and restore the default so a future in-process workflow test
// cannot inherit a route selection from this configuration test.
BOOST_AUTO_TEST_CASE(CommonTrackOutputSelectionIsDefaultFalseAndRuntimeVisible)
{
  auto& output = ITSMFTCommonCAOutputParam::Instance();
  BOOST_CHECK_EQUAL(output.getName(), "ITSMFTCommonCAOutputParam");
  BOOST_CHECK_EQUAL(output.useCommonTrackOutput, false);

  o2::conf::ConfigurableParam::setValue<bool>("ITSMFTCommonCAOutputParam", "useCommonTrackOutput", true);
  BOOST_CHECK_EQUAL(output.useCommonTrackOutput, true);
  o2::conf::ConfigurableParam::setValue<bool>("ITSMFTCommonCAOutputParam", "useCommonTrackOutput", false);
  BOOST_CHECK_EQUAL(output.useCommonTrackOutput, false);
}
