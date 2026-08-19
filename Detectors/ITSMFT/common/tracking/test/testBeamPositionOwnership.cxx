// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#define BOOST_TEST_MODULE ITSMFT BeamPositionOwnership
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <memory>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITStracking/TrackingConfigParam.h"

#include "TraversalTestSupport.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

o2::its::TrackerParamConfig& mutableLegacyITSTrackerParamConfig()
{
  return const_cast<o2::its::TrackerParamConfig&>(o2::its::TrackerParamConfig::Instance());
}

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

struct TrackerRig {
  explicit TrackerRig(bool useDiamond)
  {
    resetDetectorDefaults(parameters, o2::detectors::DetID::ITS);
    parameters.UseDiamond = useDiamond;
    parameters.Diamond[0] = 1.f;
    parameters.Diamond[1] = 2.f;
    parameters.DiamondCov[3] = 0.25f;

    TrackerInitialization configuration;
    configuration.catalog = {kITSStaticSurfaceCatalog.data(), static_cast<uint32_t>(kITSStaticSurfaceCatalog.size())};
    configuration.layout = makeDetectorLayout();
    configuration.parameters.push_back(parameters);
    configuration.memoryPool = std::make_shared<BoundedMemoryResource>();
    BOOST_REQUIRE(tracker.initialize(frame, configuration).ok());
  }

  TrackingParameters parameters;
  TimeFrame frame;
  Tracker tracker;
};

} // namespace

BOOST_FIXTURE_TEST_CASE(TrackerAppliesItsOwnedDiamondConfiguration, ScopedLegacyOverrideBeamEstimation)
{
  mutableLegacyITSTrackerParamConfig().overrideBeamEstimation = true;
  TrackerRig rig{/*useDiamond=*/true};
  rig.frame.setBeamPosition(7.f, 8.f, 0.5f);

  TrackerTestAccess::configureBeamPosition(rig.tracker, rig.frame);

  BOOST_CHECK_CLOSE(rig.frame.getBeamX(), 1.f, 1e-4);
  BOOST_CHECK_CLOSE(rig.frame.getBeamY(), 2.f, 1e-4);
  BOOST_CHECK_EQUAL(rig.frame.getBeamPositionVariance(), 0.25f);
}

BOOST_AUTO_TEST_CASE(TrackerReappliesDiamondForEveryTimeFrame)
{
  TrackerRig rig{/*useDiamond=*/true};
  TrackerTestAccess::configureBeamPosition(rig.tracker, rig.frame);
  rig.frame.setBeamPosition(-9.f, -8.f, 0.5f);

  TrackerTestAccess::configureBeamPosition(rig.tracker, rig.frame);

  BOOST_CHECK_CLOSE(rig.frame.getBeamX(), 1.f, 1e-4);
  BOOST_CHECK_CLOSE(rig.frame.getBeamY(), 2.f, 1e-4);
  BOOST_CHECK_EQUAL(rig.frame.getBeamPositionVariance(), 0.25f);
}

BOOST_AUTO_TEST_CASE(TrackerPreservesEstimatedBeamWhenDiamondIsDisabled)
{
  TrackerRig rig{/*useDiamond=*/false};
  rig.frame.setBeamPosition(7.f, 8.f, 0.5f);

  TrackerTestAccess::configureBeamPosition(rig.tracker, rig.frame);

  BOOST_CHECK_CLOSE(rig.frame.getBeamX(), 7.f, 1e-4);
  BOOST_CHECK_CLOSE(rig.frame.getBeamY(), 8.f, 1e-4);
  BOOST_CHECK_EQUAL(rig.frame.getBeamPositionVariance(), 0.5f);
}
