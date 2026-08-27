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

// Gate 3 common-CA failure contract: unit test for the pure publication
// decision `decideCATrackerPublicationAction()` factored out of
// CATrackerDPL::run() in CATrackerSpec.cxx. This is the narrowest feasible
// test seam around the DPL publication contract: CATrackerDPL::run() itself
// takes a framework::ProcessingContext, and this repository has no harness
// that constructs one outside a running DPL driver (no workflow/test
// directory anywhere under Detectors/ITSMFT exercises a DataProcessorSpec's
// run() directly; the two precedents in the whole Detectors tree --
// Detectors/TPC/workflow/test/test_TPCWorkflow.cxx and
// Detectors/TRD/workflow -- either drive a full external process or a shell
// script, not a Boost.Test call into run()). What the actual publication
// contract for CATrackerDPL::run() depends on -- once inputs are pulled and
// tracking has returned -- reduces to exactly two pieces of information,
// tracker-active and the tracking result, so factoring that decision out
// into a pure, directly testable function is the narrowest seam available
// without adding a ProcessingContext test harness to the framework itself.
// The end-to-end publication behavior (no output on a dropped TF, a live
// non-fatal device, no accepted tracks) is additionally covered by the
// workflow-level negative replay documented in
// Detectors/ITSMFT/common/tracking/test/gate0-baseline/README.md.

#define BOOST_TEST_MODULE MFT CA tracker publication decision
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/Tracker.h"
#include "MFTWorkflow/CATrackerSpec.h"

using namespace o2::mft;

BOOST_AUTO_TEST_CASE(InactiveTrackerAlwaysPublishesEmptyRegardlessOfResultValue)
{
  BOOST_CHECK(decideCATrackerPublicationAction(false, o2::itsmft::tracking::TrackingOutcome::Success) == CATrackerPublicationAction::PublishInactiveEmpty);
  BOOST_CHECK(decideCATrackerPublicationAction(false, o2::itsmft::tracking::TrackingOutcome::RecoverableDropped) == CATrackerPublicationAction::PublishInactiveEmpty);
  BOOST_CHECK(decideCATrackerPublicationAction(false, o2::itsmft::tracking::TrackingOutcome::Structural) == CATrackerPublicationAction::PublishInactiveEmpty);
}

BOOST_AUTO_TEST_CASE(ActiveTrackerWithRecoverableDropSkipsPublication)
{
  BOOST_CHECK(decideCATrackerPublicationAction(true, o2::itsmft::tracking::TrackingOutcome::RecoverableDropped) == CATrackerPublicationAction::SkipDroppedTimeFrame);
}

BOOST_AUTO_TEST_CASE(ActiveTrackerWithNonDroppedResultPublishes)
{
  BOOST_CHECK(decideCATrackerPublicationAction(true, o2::itsmft::tracking::TrackingOutcome::Success) == CATrackerPublicationAction::PublishActiveResult);
  BOOST_CHECK(decideCATrackerPublicationAction(true, o2::itsmft::tracking::TrackingOutcome::Structural) == CATrackerPublicationAction::PublishActiveResult);
}
