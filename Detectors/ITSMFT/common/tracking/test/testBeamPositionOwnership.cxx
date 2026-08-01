// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 3 ITS common-CA config-ownership correction: DetectorTraits.cxx's
// TrackingLoadPolicy<ITS, NLayers>::configureBeamPosition() must not consult
// TrackerParamRef<ITS>::get() (o2::its::TrackerParamConfig, the frozen legacy
// "ITSCATrackerParam" namespace) for overrideBeamEstimation -- the common
// workflow constructor/configuration (the overrideBeamEstimation argument
// itself, plus TrackingParameters::UseDiamond, both ultimately sourced from
// ITSCommonCATrackerParam via TrackingMode::getTrackingParameters()) is the
// sole owner of the selected static-diamond constraint. MFT is unchanged:
// its branch never reads TrackerParamRef<MFT>::get() at all (it always uses
// the diamond), so it needs no isolation and this file only re-confirms that
// with a no-regression check.
//
// TrackingLoadPolicy<DetId, NLayers>::configureBeamPosition() is already a
// public static member (DetectorTraits.h); called directly here with a bare
// TimeFrame (non-templated, Gate 4 B3.1), no geometry/dictionary/GRP fixture
// required.

#define BOOST_TEST_MODULE ITSMFT BeamPositionOwnership
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/DetectorTraits.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITStracking/TrackingConfigParam.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

o2::its::TrackerParamConfig& mutableLegacyITSTrackerParamConfig()
{
  return const_cast<o2::its::TrackerParamConfig&>(o2::its::TrackerParamConfig::Instance());
}

// Restores the legacy singleton field this test file stages, on scope exit
// (including a thrown/failed check).
struct ScopedLegacyOverrideBeamEstimation {
  ScopedLegacyOverrideBeamEstimation()
    : original(mutableLegacyITSTrackerParamConfig().overrideBeamEstimation)
  {
  }
  ~ScopedLegacyOverrideBeamEstimation()
  {
    mutableLegacyITSTrackerParamConfig().overrideBeamEstimation = original;
  }
  bool original;
};

TrackingParameters diamondParameters(float dx, float dy, float dz)
{
  TrackingParameters p;
  p.UseDiamond = true;
  p.Diamond[0] = dx;
  p.Diamond[1] = dy;
  p.Diamond[2] = dz;
  return p;
}

// Distinct from the diamond above so a test can tell which source actually
// won.
o2::dataformats::MeanVertexObject meanVertexAt(float x, float y)
{
  return o2::dataformats::MeanVertexObject{x, y, 0.f, 0.02f, 0.02f, 6.f, 0.f, 0.f};
}

} // namespace

// --- ITS: a legacy overrideBeamEstimation=true must not, by itself, steer
// the common ITS beam position onto the MeanVertex source. Only the
// constructor/configuration-supplied overrideBeamEstimation argument may. ---

BOOST_FIXTURE_TEST_CASE(LegacyOverrideBeamEstimationAloneDoesNotSelectMeanVertex, ScopedLegacyOverrideBeamEstimation)
{
  mutableLegacyITSTrackerParamConfig().overrideBeamEstimation = true;

  TimeFrame tf;
  const auto p = diamondParameters(1.f, 2.f, 3.f);
  const auto meanVertex = meanVertexAt(-9.f, -8.f);

  // overrideBeamEstimation argument itself is false: if the legacy field
  // above still leaked in (the pre-fix "overrideBeamEstimation || tc.overrideBeamEstimation"
  // condition), the beam position would come from meanVertex (-9, -8)
  // instead of the diamond (1, 2).
  TrackingLoadPolicyN<ITSNLayers>::configureBeamPosition(tf, p, &meanVertex, /*overrideBeamEstimation=*/false);

  BOOST_CHECK_CLOSE(tf.getBeamX(), 1.f, 1e-4);
  BOOST_CHECK_CLOSE(tf.getBeamY(), 2.f, 1e-4);
}

BOOST_FIXTURE_TEST_CASE(ConstructorOverrideBeamEstimationStillSelectsMeanVertex, ScopedLegacyOverrideBeamEstimation)
{
  // The legacy field is left at its default (false) here: the common
  // workflow's own overrideBeamEstimation argument is the sole owner and
  // must still work on its own, unassisted by the legacy field.
  BOOST_REQUIRE_EQUAL(mutableLegacyITSTrackerParamConfig().overrideBeamEstimation, false);

  TimeFrame tf;
  const auto p = diamondParameters(1.f, 2.f, 3.f);
  const auto meanVertex = meanVertexAt(-9.f, -8.f);

  TrackingLoadPolicyN<ITSNLayers>::configureBeamPosition(tf, p, &meanVertex, /*overrideBeamEstimation=*/true);

  BOOST_CHECK_CLOSE(tf.getBeamX(), -9.f, 1e-4);
  BOOST_CHECK_CLOSE(tf.getBeamY(), -8.f, 1e-4);
}

BOOST_FIXTURE_TEST_CASE(ITSFallsBackToDiamondWhenNeitherOverrideIsSet, ScopedLegacyOverrideBeamEstimation)
{
  mutableLegacyITSTrackerParamConfig().overrideBeamEstimation = false;

  TimeFrame tf;
  const auto p = diamondParameters(1.f, 2.f, 3.f);
  const auto meanVertex = meanVertexAt(-9.f, -8.f);

  TrackingLoadPolicyN<ITSNLayers>::configureBeamPosition(tf, p, &meanVertex, /*overrideBeamEstimation=*/false);

  BOOST_CHECK_CLOSE(tf.getBeamX(), 1.f, 1e-4);
  BOOST_CHECK_CLOSE(tf.getBeamY(), 2.f, 1e-4);
}

// --- MFT: no regression. The MFT branch always uses the diamond and never
// reads overrideBeamEstimation (constructor argument or legacy field). ---

BOOST_FIXTURE_TEST_CASE(MFTAlwaysUsesDiamondRegardlessOfOverrideArguments, ScopedLegacyOverrideBeamEstimation)
{
  mutableLegacyITSTrackerParamConfig().overrideBeamEstimation = true; // must have no effect on MFT either

  TimeFrame tf;
  TrackingParameters p;
  p.Diamond[0] = 4.f;
  p.Diamond[1] = 5.f;
  p.Diamond[2] = 6.f;
  const auto meanVertex = meanVertexAt(-9.f, -8.f);

  TrackingLoadPolicyN<o2::mft::constants::mft::LayersNumber>::configureBeamPosition(tf, p, &meanVertex, /*overrideBeamEstimation=*/true);

  BOOST_CHECK_CLOSE(tf.getBeamX(), 4.f, 1e-4);
  BOOST_CHECK_CLOSE(tf.getBeamY(), 5.f, 1e-4);
}
