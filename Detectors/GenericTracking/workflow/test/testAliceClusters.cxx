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

/// \file testAliceClusters.cxx
/// \brief The dictionary-free half of the O2 cluster conversion: readout
///        frames, times and the sensor frame (GenericTracking/Alice/Clusters.h).
///
/// The dictionary half is held to O2 on real clusters by
/// `o2-gtrk-cluster-closure`.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK ALICE clusters
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Base/TimeEstPS.h"

#include "CommonConstants/LHCConstants.h"
#include "DataFormatsITSMFT/ROFRecord.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace o2::gtrk
{

namespace
{

using o2::itsmft::ROFRecord;

/// \a count frames of \a lengthBC from \a firstBC of orbit \a orbit, with
/// \a perFrame clusters each.
std::vector<ROFRecord> frames(int count, int lengthBC, int firstBC, std::uint32_t orbit, int perFrame)
{
  std::vector<ROFRecord> result;
  for (int i = 0; i < count; ++i) {
    const int bc = firstBC + (i * lengthBC);
    const o2::InteractionRecord ir{static_cast<std::uint16_t>(bc % o2::constants::lhc::LHCMaxBunches),
                                   orbit + static_cast<std::uint32_t>(bc / o2::constants::lhc::LHCMaxBunches)};
    result.emplace_back(ir, i, i * perFrame, perFrame);
  }
  return result;
}

} // namespace

BOOST_AUTO_TEST_CASE(FrameLengthIsTheRecordsSpacing)
{
  BOOST_TEST(alice::readoutFrames(frames(40, 198, 64, 0U, 0)).lengthBC == 198);
  BOOST_TEST(alice::readoutFrames(frames(3, 891, 0, 7U, 0)).lengthBC == 891);
  BOOST_CHECK_THROW(static_cast<void>(alice::readoutFrames(frames(1, 198, 0, 0U, 0))), std::runtime_error);
  auto uneven = frames(4, 198, 0, 0U, 0);
  uneven.back().getBCData() += 1;
  BOOST_CHECK_THROW(static_cast<void>(alice::readoutFrames(uneven)), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ClustersTakeTheirFrameAsTimeBoxAndOurSensorFrame)
{
  // Two frames of 198 BC in a timeframe starting at orbit 32, the second
  // across the orbit boundary.
  const auto rofs = frames(2, 198, 3500, 32U, 2);
  const std::vector<alice::LocalCluster> local{
    {.chip = 3, .x = 0.1F, .z = -0.5F, .sigma2X = 4e-8F, .sigma2Z = 9e-8F, .nPixels = 2U},
    {.chip = 4, .x = -0.3F, .z = 1.2F, .sigma2X = 1e-8F, .sigma2Z = 1e-8F, .nPixels = 1U},
    {.chip = 5, .x = 0.F, .z = 0.F, .sigma2X = 1e-8F, .sigma2Z = 1e-8F, .nPixels = 7U},
    {.chip = 6, .x = 0.6F, .z = -1.4F, .sigma2X = 1e-8F, .sigma2Z = 1e-8F, .nPixels = 3U},
  };
  const auto clusters = alice::toClusters(local, rofs, {.lengthBC = 198}, 32U);
  BOOST_REQUIRE(clusters.size() == local.size());

  BOOST_TEST(clusters[0].chipID == 3U);
  BOOST_TEST(std::abs(static_cast<double>(clusters[0].xLocal) - (0.1 - alice::kActiveCentreX)) < 1e-7);
  BOOST_TEST(clusters[0].yLocal == -0.5F);
  BOOST_TEST(std::abs(clusters[0].sigmaXLocal - 2e-4F) < 1e-9F);
  BOOST_TEST(std::abs(clusters[0].sigmaYLocal - 3e-4F) < 1e-9F);

  // Each cluster's readout frame as its time box, its start and its length
  // rounded down; the second frame across the orbit boundary.
  BOOST_TEST(clusters[0].timePS == static_cast<std::uint32_t>(psOfBC(3500)));
  BOOST_TEST(clusters[0].widthPS == static_cast<std::uint32_t>(psOfBC(198)));
  BOOST_TEST(clusters[1].timePS == clusters[0].timePS);
  BOOST_TEST(clusters[2].timePS == static_cast<std::uint32_t>(psOfBC(3500 + 198)));
  BOOST_TEST(clusters[2].widthPS == clusters[0].widthPS);
  BOOST_TEST(clusters[3].timePS == clusters[2].timePS);
}

BOOST_AUTO_TEST_CASE(FramesMustCoverTheClustersInOrder)
{
  const std::vector<alice::LocalCluster> local(4);
  BOOST_CHECK_THROW(static_cast<void>(alice::toClusters(local, frames(1, 198, 0, 0U, 2), {.lengthBC = 198}, 0U)), std::runtime_error);
  BOOST_CHECK_THROW(static_cast<void>(alice::toClusters(local, frames(3, 198, 0, 0U, 2), {.lengthBC = 198}, 0U)), std::runtime_error);
  // A frame before the timeframe's first orbit.
  BOOST_CHECK_THROW(static_cast<void>(alice::toClusters(local, frames(2, 198, 0, 0U, 2), {.lengthBC = 198}, 1U)), std::runtime_error);
  BOOST_CHECK_NO_THROW(static_cast<void>(alice::toClusters(local, frames(2, 198, 0, 1U, 2), {.lengthBC = 198}, 1U)));
}

} // namespace o2::gtrk
