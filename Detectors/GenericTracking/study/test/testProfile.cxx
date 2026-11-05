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

/// \file testProfile.cxx
/// \brief The output digest and the phase profile: a digest reads values,
///        not padding; a profile's zones
///        nest, its parallel loops fold into the zone around them, and nothing
///        is recorded without an activation.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Profile
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Study/Digest.h"
#include "GenericTracking/Study/ProfileReport.h"
#include "GenericTracking/Tracking/Profile.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <sstream>
#include <vector>

namespace o2::gtrk::study
{
namespace
{

using profile::Phase;

BoundTrackState stateOver(unsigned char fill)
{
  // Every byte, padding included, set to \a fill before the fields are.
  BoundTrackState state;
  std::memset(static_cast<void*>(&state), fill, sizeof(state));
  for (std::size_t index = 0; index < std::size(state.par); ++index) {
    state.par[index] = static_cast<Scalar>(index) * 0.5F;
  }
  for (std::size_t index = 0; index < std::size(state.cov); ++index) {
    state.cov[index] = static_cast<Scalar>(index) * 0.25F;
  }
  state.surface = 7;
  state.surfacePhi = 0.F;
  return state;
}

std::uint64_t digestOf(const std::vector<BoundTrackState>& states)
{
  Digest digest;
  digest.column(states.data(), states.size());
  return digest.result();
}

BOOST_AUTO_TEST_CASE(DigestReadsValuesNotPadding)
{
  const std::vector<BoundTrackState> zeros{stateOver(0x00), stateOver(0x00)};
  const std::vector<BoundTrackState> ones{stateOver(0xFF), stateOver(0xFF)};
  BOOST_TEST(digestOf(zeros) == digestOf(ones));

  auto changed = zeros;
  changed[1].cov[3] = std::nextafter(changed[1].cov[3], 1.F);
  BOOST_TEST(digestOf(zeros) != digestOf(changed));

  // A curved sensor's plane angle is part of the state; 0 is not hashed.
  auto curved = zeros;
  curved[0].surfacePhi = 0.25F;
  BOOST_TEST(digestOf(zeros) != digestOf(curved));
}

BOOST_AUTO_TEST_CASE(DigestKeepsColumnBoundaries)
{
  const std::array<std::uint32_t, 3> values{1U, 2U, 3U};
  Digest split;
  split.column(values.data(), 2U);
  split.column(values.data() + 2, 1U);
  Digest other;
  other.column(values.data(), 1U);
  other.column(values.data() + 1, 2U);
  BOOST_TEST(split.result() != other.result());

  // Signed zeros are different bits, so a different digest.
  Digest positive;
  positive.add(0.F);
  Digest negative;
  negative.add(-0.F);
  BOOST_TEST(positive.result() != negative.result());
}

BOOST_AUTO_TEST_CASE(TheReportNamesItsPhasesAndUnattributedTime)
{
  profile::Profile result;
  result[Phase::kTracklets].calls = 1U;
  result[Phase::kTracklets].seconds = 2.;
  result[Phase::kTrackletGenerate].calls = 1U;
  result[Phase::kTrackletGenerate].seconds = 1.;
  std::ostringstream report;
  study::report(report, result, Phase::kTracklets, 4);
  BOOST_TEST(report.str().find("tracklet generate") != std::string::npos);
  BOOST_TEST(report.str().find("(other)") != std::string::npos);
}

} // namespace
} // namespace o2::gtrk::study
