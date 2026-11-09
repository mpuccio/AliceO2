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

/// \file testClusterShapeBand.cxx
/// \brief A cluster-shape band's slope interval per cluster shape, and a
///        line's local slope against one.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK cluster shape band
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/ShapeCut.h"
#include "GenericTracking/Workflow/ClusterShapeBand.h"

#include <cstddef>
#include <vector>

using o2::gtrk::workflow::ClusterShapeBand;

namespace
{

/// A band of three layers, layer 1 calibrated: at most two rows; one column
/// at slopes [0, 1), two [0, 2) (the hull of two bins), five [1, 2), three
/// none.
ClusterShapeBand band()
{
  ClusterShapeBand result{
    .maxRows = {0, 2, 0},
    .slopeLow = std::vector<float>(3 * ClusterShapeBand::kSpans, ClusterShapeBand::kNone.first),
    .slopeHigh = std::vector<float>(3 * ClusterShapeBand::kSpans, ClusterShapeBand::kNone.second),
  };
  const auto set = [&](std::size_t span, float low, float high) {
    result.slopeLow[ClusterShapeBand::kSpans + span] = low;
    result.slopeHigh[ClusterShapeBand::kSpans + span] = high;
  };
  set(1, 0.F, 1.F);
  set(2, 0.F, 2.F);
  set(5, 1.F, 2.F);
  return result;
}

} // namespace

// A cluster's shape gives its interval on a calibrated layer; more rows than
// the band's, a span no slope gives, a cluster without a pattern and an
// uncalibrated layer, or one beyond the band's, give none or every one.
BOOST_AUTO_TEST_CASE(AClusterShapeGivesItsSlopeInterval)
{
  const auto shapes = band();
  BOOST_TEST(shapes.fits(3U));
  BOOST_TEST(!shapes.fits(2U));
  BOOST_TEST((shapes.interval(5U, 1U, 1U) == ClusterShapeBand::kEvery));
  BOOST_TEST((shapes.interval(1U, 1U, 1U) == ClusterShapeBand::Interval(0.F, 1.F)));
  BOOST_TEST((shapes.interval(1U, 2U, 2U) == ClusterShapeBand::Interval(0.F, 2.F)));
  BOOST_TEST((shapes.interval(1U, 1U, 5U) == ClusterShapeBand::Interval(1.F, 2.F)));
  BOOST_TEST((shapes.interval(1U, 3U, 1U) == ClusterShapeBand::kNone));
  BOOST_TEST((shapes.interval(1U, 1U, 3U) == ClusterShapeBand::kNone));
  BOOST_TEST((shapes.interval(1U, 1U, 0U) == ClusterShapeBand::kEvery));
  BOOST_TEST((shapes.interval(0U, 1U, 2U) == ClusterShapeBand::kEvery));
}
