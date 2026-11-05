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

/// \file testVertexQA.cxx
/// \brief Vertex truth assignment, coverage, splits, fakes and merges.
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK VertexQA
#include <boost/test/unit_test.hpp>
#include "GenericTracking/Study/VertexQA.h"
#include "ToyBarrel8.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <vector>

namespace o2::gtrk::study
{
BOOST_AUTO_TEST_CASE(ClaimsDetermineCoverageSplitsFakesAndMerges)
{
  std::vector<TruthRow> rows;
  ParticleTable particles;
  // Two findable primaries in separate collisions, each measured on four
  // layers. The second collision contributes to the first's vertex.
  for (std::uint64_t collision = 0; collision < 2; ++collision) {
    const auto particle = (collision << o2::MCCompLabel::nbitsTrackID) | 1U;
    particles[particle] = {.pt = 1.F, .primary = true, .vz = 0.F};
    for (LayerID layer = 0; layer < 4; ++layer) {
      TruthRow row;
      row.particles[0] = particle;
      row.count = 1U;
      row.layer = layer;
      rows.push_back(row);
    }
  }
  rows.emplace_back(); // noise
  const std::vector<TripletClaim> claims{
    {.rows = {0, 1, 2}, .vertex = 0}, {.rows = {1, 2, 3}, .vertex = 0}, {.rows = {0, 1, 3}, .vertex = 1}, {.rows = {4, 5, 6}, .vertex = 0}, {.rows = {8, 8, 8}, .vertex = 2}};
  VertexerResult found;
  found.vertices = {{.z = 0.01F, .varZ = 0.0001F}, {.z = 0.02F, .varZ = 0.0001F}, {.z = 3.F}};
  found.contributors = {3U, 1U, 1U};
  const auto score = scoreVertices<detectors::ToyBarrel8>(rows, particles, claims, found);
  BOOST_TEST(score.vertices == 3U);
  BOOST_TEST(score.collisions == 2U);
  BOOST_TEST(score.found == 1U);
  BOOST_TEST(score.findable == 2U);
  BOOST_TEST(score.covered == 1U);
  BOOST_TEST(score.fakes == 1U);
  BOOST_TEST(score.splits == 1U);
  BOOST_TEST(score.merged == 1U);
  BOOST_REQUIRE(score.residualsUM.size() == 1U);
  BOOST_TEST(score.residualsUM[0] == 100., boost::test_tools::tolerance(1e-6));
  BOOST_REQUIRE(score.pulls.size() == 1U);
  BOOST_TEST(score.pulls[0] == 1., boost::test_tools::tolerance(1e-6));
  // Of two vertices of a collision the contributor maximum supplies its
  // residual; this split has more contributors than the original vertex.
  found.contributors[1] = 4U;
  const auto splitWins = scoreVertices<detectors::ToyBarrel8>(rows, particles, claims, found);
  BOOST_TEST(splitWins.residualsUM[0] == 200., boost::test_tools::tolerance(1e-6));
}

BOOST_AUTO_TEST_CASE(TiesKeepTheFirstCollisionAndFirstVertex)
{
  const std::vector<std::map<std::int64_t, std::uint32_t>> labels{{{0, 2U}, {1, 2U}}, {{0, 1U}}, {{-1, 1U}, {1, 1U}}};
  const std::array<std::uint32_t, 3> contributors{4U, 4U, 2U};
  const std::vector<std::map<std::int64_t, std::uint32_t>> wentTo{{{0, 2U}, {1, 1U}}, {{0, 2U}, {2, 1U}}};
  const auto score = scoreVertexClaims(labels, contributors, wentTo);
  BOOST_TEST(score.collision == (std::vector<std::int64_t>{0, 0, -1}), boost::test_tools::per_element());
  BOOST_TEST(score.best == (std::vector<std::int64_t>{0, -1}), boost::test_tools::per_element());
  BOOST_TEST(score.fakes == 1U);
  BOOST_TEST(score.splits == 1U);
  BOOST_TEST(score.merged == 1U);
}
} // namespace o2::gtrk::study
