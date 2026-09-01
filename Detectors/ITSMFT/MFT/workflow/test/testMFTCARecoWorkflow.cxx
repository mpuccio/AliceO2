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

#define BOOST_TEST_MODULE MFTCARecoWorkflow
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <string_view>

#include "MFTWorkflow/CARecoWorkflow.h"

namespace
{
bool hasDevice(const o2::framework::WorkflowSpec& workflow, std::string_view name)
{
  return std::any_of(workflow.begin(), workflow.end(), [name](const auto& spec) { return spec.name == name; });
}
} // namespace

BOOST_AUTO_TEST_CASE(DefaultWorkflowIsMonolithic)
{
  const auto workflow = o2::mft::ca_reco_workflow::getWorkflow(
    false, false, false, false, false, false, false, false, false, true, true,
    o2::itsmft::TrackingMode::Sync, false);

  BOOST_CHECK(hasDevice(workflow, "mft-digit-reader"));
  BOOST_CHECK(hasDevice(workflow, "mft-clusterer"));
  BOOST_CHECK(hasDevice(workflow, "mft-cluster-writer"));
  BOOST_CHECK(hasDevice(workflow, "mft-ca-tracker"));
  BOOST_CHECK(hasDevice(workflow, "mft-track-writer"));
  BOOST_CHECK(!hasDevice(workflow, "mft-tracker"));
}

BOOST_AUTO_TEST_CASE(UpstreamClustersCanRunTrackerOnly)
{
  const auto workflow = o2::mft::ca_reco_workflow::getWorkflow(
    false, false, false, false, false, true, false, true, false, true, true,
    o2::itsmft::TrackingMode::Sync, false);

  BOOST_REQUIRE_EQUAL(workflow.size(), 1);
  BOOST_CHECK_EQUAL(workflow.front().name, "mft-ca-tracker");
}
