// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 3 ITS common-CA config-ownership correction: TrackingInterface.cxx's
// resolveTrackingParameters() must not consult TrackerParamRef<ITS>::get()
// (o2::its::TrackerParamConfig, the frozen legacy "ITSCATrackerParam"
// namespace) for its tracking mode -- the common workflow
// constructor/setTrackingMode() is the sole owner of mTrackingMode for ITS.
// MFT is unchanged: TrackerParamRef<MFT>::get() (TrackerParamConfig<MFT>,
// "MFTCATrackerParam") is MFT's own live configuration, not a legacy
// namespace, and stays authoritative for MFT's tracking mode exactly as
// before this correction.

#define BOOST_TEST_MODULE ITSMFT TrackingModeOwnership
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <TGeoGlobalMagField.h>
#include <fairlogger/Logger.h>
#include <stdexcept>

#include "DataFormatsParameters/GRPECSObject.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "Field/MagneticField.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InputSpec.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITSMFTTracking/TrackingInterface.h"
#include "ITStracking/TrackingConfigParam.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

// --- Global fixtures needed only for the MFT regression case: MFT's shared
// tail in TrackingMode::getTrackingParameters() dereferences
// o2::base::Propagator::Instance() and configureROFLookupTables() needs a
// GRPECS singleton (same requirement/pattern as
// testTrackingInterfaceLoadFailureContract.cxx). Not needed for ITS: the ITS
// Sync branch touches neither. ---

struct FieldFixture {
  FieldFixture()
  {
    if (!TGeoGlobalMagField::Instance()->GetField()) {
      TGeoGlobalMagField::Instance()->SetField(o2::field::MagneticField::createNominalField(5, true));
      TGeoGlobalMagField::Instance()->Lock();
    }
  }
};

struct GRPECSFixture {
  static o2::parameters::GRPECSObject& grpEcs()
  {
    static o2::parameters::GRPECSObject obj;
    return obj;
  }

  GRPECSFixture()
  {
    auto& obj = grpEcs();
    obj.addDetContinuousReadOut(o2::detectors::DetID::MFT);
    obj.setNHBFPerTF(128);

    std::vector<o2::framework::InputSpec> inputs;
    auto req = std::make_shared<o2::base::GRPGeomRequest>(
      false, true, false, false, false, o2::base::GRPGeomRequest::None, inputs);
    o2::base::GRPGeomHelper::instance().setRequest(req);
    o2::framework::ConcreteDataMatcher matcher{"GLO", "GRPECS", 0};
    o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, &obj);
  }
};

BOOST_GLOBAL_FIXTURE(FieldFixture);
BOOST_GLOBAL_FIXTURE(GRPECSFixture);

// ConfigurableParamHelper<P>::Instance() deliberately returns const P& to
// discourage casual mutation from production code; tests that need to stage
// a specific singleton value use the same const_cast pattern as e.g.
// testITSCommonCANThreads.cxx's mutable*() helpers.
ITSCommonCATrackerParam& mutableITSCommonCATrackerParam()
{
  return const_cast<ITSCommonCATrackerParam&>(ITSCommonCATrackerParam::Instance());
}

o2::its::TrackerParamConfig& mutableLegacyITSTrackerParamConfig()
{
  return const_cast<o2::its::TrackerParamConfig&>(o2::its::TrackerParamConfig::Instance());
}

TrackerParamConfig<o2::detectors::DetID::MFT>& mutableMFTTrackerParamConfig()
{
  return const_cast<TrackerParamConfig<o2::detectors::DetID::MFT>&>(
    TrackerParamConfig<o2::detectors::DetID::MFT>::Instance());
}

// Restores every singleton field this test file stages, on scope exit
// (including a thrown/failed check), so no test case here can leak state
// into another test case in this binary. RAII, matching
// testITSCommonCANThreads.cxx's ScopedParams.
struct ScopedParams {
  ScopedParams()
    : originalITSUseDiamond(mutableITSCommonCATrackerParam().useDiamond),
      originalLegacyITSTrackingMode(mutableLegacyITSTrackerParamConfig().trackingMode),
      originalMFTTrackingMode(mutableMFTTrackerParamConfig().trackingMode)
  {
    // ITS Sync's getTrackingParameters() itself does not require
    // useDiamond=true, but matches the established convention (see
    // testITSCommonCANThreads.cxx) so this file stays consistent with it.
    mutableITSCommonCATrackerParam().useDiamond = true;
  }
  ~ScopedParams()
  {
    mutableITSCommonCATrackerParam().useDiamond = originalITSUseDiamond;
    mutableLegacyITSTrackerParamConfig().trackingMode = originalLegacyITSTrackingMode;
    mutableMFTTrackerParamConfig().trackingMode = originalMFTTrackingMode;
  }
  bool originalITSUseDiamond;
  int originalLegacyITSTrackingMode;
  int originalMFTTrackingMode;
};

// LOGP(fatal, ...) normally terminates the process (FairLogger default). A
// process-local OnFatal handler converts it into a catchable exception, the
// same pattern testITSCommonCATrackingModeConfiguration.cxx already uses.
// Each ITSMFT test source file builds its own executable (o2_add_test == one
// binary per SOURCES file), so this handler cannot leak into unrelated test
// binaries.
struct FatalToExceptionFixture {
  FatalToExceptionFixture() { fair::Logger::OnFatal([]() { throw std::runtime_error("fatal"); }); }
};

} // namespace

// --- ITS: a legacy trackingMode override must not reach the common ITS
// interface. TrackingMode::getTrackingParameters(ITS, mode) LOGP(fatal, ...)s
// for every mode except Sync (workflow-onboarding Slice 1) -- if
// resolveTrackingParameters() leaked the legacy Async/Cosmics value through,
// initialise() would throw instead of succeeding. ---

BOOST_FIXTURE_TEST_CASE(ITSCommonSyncInitialisationIgnoresLegacyAsyncMode, ScopedParams)
{
  mutableLegacyITSTrackerParamConfig().trackingMode = static_cast<int>(o2::itsmft::TrackingMode::Async);

  ITSMFTTrackingInterface<ITSNLayers> interface{false, o2::itsmft::TrackingMode::Sync, false};
  BOOST_CHECK_NO_THROW(interface.initialise());
  BOOST_CHECK(interface.isActive());
  BOOST_REQUIRE_EQUAL(interface.getTrackingParameters().size(), 1u);
}

BOOST_FIXTURE_TEST_CASE(ITSCommonSyncInitialisationIgnoresLegacyCosmicsMode, ScopedParams)
{
  mutableLegacyITSTrackerParamConfig().trackingMode = static_cast<int>(o2::itsmft::TrackingMode::Cosmics);

  ITSMFTTrackingInterface<ITSNLayers> interface{false, o2::itsmft::TrackingMode::Sync, false};
  BOOST_CHECK_NO_THROW(interface.initialise());
  BOOST_CHECK(interface.isActive());
  BOOST_REQUIRE_EQUAL(interface.getTrackingParameters().size(), 1u);
}

BOOST_FIXTURE_TEST_CASE(ITSCommonSyncInitialisationIgnoresLegacyAsyncModeEvenWhenModeIsUnset, ScopedParams)
{
  // The Unset-default constructor path is the one that used to defer to the
  // legacy configurable param (see resolveTrackingParameters()'s pre-fix
  // "mode == Unset" branch); confirm it now resolves to Sync unconditionally
  // for ITS instead of adopting the legacy override.
  mutableLegacyITSTrackerParamConfig().trackingMode = static_cast<int>(o2::itsmft::TrackingMode::Async);

  ITSMFTTrackingInterface<ITSNLayers> interface{false, o2::itsmft::TrackingMode::Unset, false};
  BOOST_CHECK_NO_THROW(interface.initialise());
  BOOST_CHECK(interface.isActive());
  BOOST_REQUIRE_EQUAL(interface.getTrackingParameters().size(), 1u);
}

BOOST_FIXTURE_TEST_CASE(EveryUnsupportedITSModeStillFailsClosedRegardlessOfLegacyParam, FatalToExceptionFixture)
{
  ScopedParams guard;
  // Confirms the isolation fix did not accidentally make ITS silently
  // tolerant of unsupported modes: a directly requested unsupported mode
  // must still fail, matching testITSCommonCATrackingModeConfiguration.cxx's
  // EveryUnsupportedITSModeFailsClosed at the getTrackingParameters() level.
  mutableLegacyITSTrackerParamConfig().trackingMode = static_cast<int>(o2::itsmft::TrackingMode::Sync);
  ITSMFTTrackingInterface<ITSNLayers> interface{false, o2::itsmft::TrackingMode::Cosmics, false};
  BOOST_CHECK_THROW(interface.initialise(), std::runtime_error);
}

// --- MFT: no regression. TrackerParamRef<MFT>::get().trackingMode
// ("MFTCATrackerParam", MFT's own live configuration, not a legacy
// namespace) must remain authoritative for MFT exactly as before this
// correction. Cosmics vs. the Sync/default baseline is distinguished via
// PVres (Cosmics sets 1.e5f; TrackingParameters' own struct default is
// 1.e-2f) -- confirmed empirically to survive
// getTrackingParameters()'s shared per-mode tail (unlike e.g. ColBins/
// RowBins, which that tail unconditionally overwrites from
// TrackerParamConfig<MFT>'s own LUTbinsU/LUTbinsV defaults regardless of
// mode). Not distinguished via iteration count either:
// TrackerParamConfig<MFT>::nIterations defaults to 1, so Cosmics' single
// iteration is not distinguishable from Sync's by size() alone. ---

BOOST_FIXTURE_TEST_CASE(MFTStillConsultsItsOwnLiveTrackingModeConfig, ScopedParams)
{
  mutableMFTTrackerParamConfig().trackingMode = static_cast<int>(o2::itsmft::TrackingMode::Cosmics);

  ITSMFTTrackingInterface<o2::mft::constants::mft::LayersNumber> interface{
    false, o2::itsmft::TrackingMode::Unset, false};
  BOOST_CHECK_NO_THROW(interface.initialise());
  BOOST_REQUIRE(interface.isActive());
  BOOST_REQUIRE_EQUAL(interface.getTrackingParameters().size(), 1u);
  BOOST_CHECK_CLOSE(interface.getTrackingParameters()[0].PVres, 1.e5f, 1e-4); // Cosmics-specific, not the Sync/default 1.e-2f
}

BOOST_FIXTURE_TEST_CASE(MFTExplicitModeStillOverridesLegacyDefaultUnset, ScopedParams)
{
  // mTrackingMode == Sync (explicit, not Unset): the legacy param is left at
  // its Unset(-1) default, so this exercises the same "mode already decided,
  // legacy has nothing to add" path the MFT branch has always had.
  BOOST_REQUIRE_EQUAL(mutableMFTTrackerParamConfig().trackingMode, -1);

  ITSMFTTrackingInterface<o2::mft::constants::mft::LayersNumber> interface{
    false, o2::itsmft::TrackingMode::Sync, false};
  BOOST_CHECK_NO_THROW(interface.initialise());
  BOOST_REQUIRE(interface.isActive());
  BOOST_REQUIRE_EQUAL(interface.getTrackingParameters().size(), 1u);
  BOOST_CHECK_CLOSE(interface.getTrackingParameters()[0].PVres, 1.e-2f, 1e-4); // Sync/default, not Cosmics' 1.e5f
}
