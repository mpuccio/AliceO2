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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE IOTOFClusterTiming
#include <boost/test/unit_test.hpp>

#include "CommonConstants/LHCConstants.h"
#include "IOTOFReconstruction/Clusterer.h"
#include <TMemFile.h>

#include <cmath>
#include <vector>

using namespace o2::iotof;

BOOST_AUTO_TEST_CASE(FineTimeSurvivesSecondsAndROOTIO)
{
  // Well beyond a second; 1 ps still separates two clusters on the same BC.
  constexpr uint64_t bc = 8000000000ULL;
  Cluster first(1, 2, 1, 1, 1, 0, 3, 24.123F, bc);
  auto second = first;
  second.setTime(first.getTime() + 0.001F);
  BOOST_TEST(first.getBc() == bc);
  BOOST_TEST(std::abs((second.getTime() - first.getTime()) - 0.001) < 0.000002);
  second.setBc(bc + 1);
  BOOST_CHECK(!(first == second));
  TMemFile file("timing.root", "RECREATE");
  file.WriteObject(&first, "cluster");
  auto* read = file.Get<Cluster>("cluster");
  BOOST_REQUIRE(read != nullptr);
  BOOST_CHECK(*read == first);
}

BOOST_AUTO_TEST_CASE(ClustererUsesFineTDCAtLargeBC)
{
  Clusterer parent;
  Clusterer::ClustererThread thread(&parent);
  constexpr uint64_t bc = 8000000000ULL;
  const double bin = DPLDigitizerParam::Instance().tdcBin;
  // Deliberately unusable absolute digit times: BC/TDC is authoritative.
  std::vector<Digit> digits{Digit(0, 3, 4, 100, 0., bc, 123)};
  thread.findClustersSingleHit(digits, 0, nullptr, nullptr);
  BOOST_REQUIRE_EQUAL(thread.mClusters.size(), 1);
  BOOST_TEST(thread.mClusters[0].getBc() == bc);
  BOOST_TEST(std::abs(thread.mClusters[0].getTime() - 123 * bin) < 0.000001);
}

BOOST_AUTO_TEST_CASE(ClusterMeanAndNeighbourGateCrossBCBoundary)
{
  Clusterer parent;
  Clusterer::ClustererThread thread(&parent);
  constexpr uint64_t bc = 8000000000ULL;
  constexpr double spacing = o2::constants::lhc::LHCBunchSpacingNS;
  const double bin = DPLDigitizerParam::Instance().tdcBin;
  const auto last = static_cast<int>(std::floor(spacing / bin));
  std::vector<Digit> digits{Digit(0, 3, 4, 100, 0., bc, last),
                            Digit(0, 3, 5, 100, 0., bc + 1, 1)};
  std::vector<uint32_t> indices{0, 1};
  const double expected = (last * bin + spacing + bin) / 2.;
  thread.findClustersMultipleHits(digits, indices, nullptr, nullptr);
  BOOST_REQUIRE_EQUAL(thread.mClusters.size(), 1);
  const auto& cluster = thread.mClusters.front();
  const double actual = static_cast<double>(cluster.getBc() - bc) * spacing + cluster.getTime();
  BOOST_TEST(std::abs(actual - expected) < 0.000002);
  BOOST_TEST(cluster.getTime() >= 0.F);
  BOOST_TEST(cluster.getTime() < spacing);
  // Same pixels in the next BC, far from the boundary, are separate hits.
  digits = {Digit(0, 3, 4, 100, 0., bc, 10), Digit(0, 3, 5, 100, 0., bc + 1, 10)};
  auto groups = thread.buildPreclusters(digits, indices, 3, 0.02F);
  BOOST_TEST(groups.size() == 2);
}
