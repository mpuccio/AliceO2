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

/// \file testDeviceStorage.cxx
/// \brief The device algorithms and storage every pass shares, against the
///        obvious host loops.
///
/// Built once per backend. Needs the device.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK device storage
#include <boost/test/unit_test.hpp>

#include "DeviceTestFixture.h"
#include "DeviceTestKernels.h"

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/Tracklet.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace o2::gtrk
{
namespace
{

using namespace device_test;

bool sameTracklet(const Tracklet& lhs, const Tracklet& rhs)
{
  return lhs.source == rhs.source && lhs.target == rhs.target && lhs.link == rhs.link &&
         std::bit_cast<std::uint32_t>(lhs.lambda) == std::bit_cast<std::uint32_t>(rhs.lambda) &&
         std::bit_cast<std::uint32_t>(lhs.phi) == std::bit_cast<std::uint32_t>(rhs.phi);
}

} // namespace

BOOST_AUTO_TEST_CASE(SelectIndicesKeepsTheOrder)
{
  Device device;
  std::mt19937 random{5U};
  for (const std::size_t count : {std::size_t{0}, std::size_t{1}, std::size_t{1000}, std::size_t{300000}}) {
    std::vector<std::uint8_t> mask(count);
    std::ranges::generate(mask, [&] { return static_cast<std::uint8_t>(random() % 3U == 0U ? 1U : 0U); });
    std::vector<std::uint32_t> expected;
    for (std::uint32_t index = 0; index < count; ++index) {
      if (mask[index] != 0U) {
        expected.push_back(index);
      }
    }
    std::vector<std::uint32_t> selected;
    BOOST_REQUIRE((kernels().selectIndices(*device.runtime, device.stream, *device.arena, mask, selected) == StageStatus::kComplete));
    BOOST_TEST(selected == expected, boost::test_tools::per_element());
  }
  BOOST_TEST(device.arena->used(ArenaEnd::kScratch) == 0U);
}

BOOST_AUTO_TEST_CASE(ExclusiveScanIn32And64Bits)
{
  Device device;
  std::mt19937 random{6U};
  std::vector<std::uint32_t> narrow(100001U);
  std::ranges::generate(narrow, [&] { return static_cast<std::uint32_t>(random() % 50U); });
  std::vector<std::uint32_t> narrowExpected(narrow.size());
  std::exclusive_scan(narrow.begin(), narrow.end(), narrowExpected.begin(), std::uint32_t{0U});
  std::vector<std::uint32_t> narrowResult;
  BOOST_REQUIRE((kernels().scan32(*device.runtime, device.stream, *device.arena, narrow, narrowResult) == StageStatus::kComplete));
  BOOST_TEST(narrowResult == narrowExpected, boost::test_tools::per_element());

  // Counts whose total passes 32 bits: the wide scan keeps it.
  std::vector<std::uint64_t> wide(4097U, std::uint64_t{std::numeric_limits<std::uint32_t>::max()} / 1000U);
  std::vector<std::uint64_t> wideExpected(wide.size());
  std::exclusive_scan(wide.begin(), wide.end(), wideExpected.begin(), std::uint64_t{0U});
  BOOST_REQUIRE(wideExpected.back() > std::numeric_limits<std::uint32_t>::max());
  std::vector<std::uint64_t> wideResult;
  BOOST_REQUIRE((kernels().scan64(*device.runtime, device.stream, *device.arena, wide, wideResult) == StageStatus::kComplete));
  BOOST_TEST(wideResult == wideExpected, boost::test_tools::per_element());
  BOOST_TEST(device.arena->used(ArenaEnd::kScratch) == 0U);
}

/// A chain reaches its fixed point in at most its length's passes, whatever
/// the schedule; the second flag is sticky; a bound too tight throws.
BOOST_AUTO_TEST_CASE(SweepToFixedPointStopsAndChecksItsBound)
{
  Device device;
  constexpr std::uint32_t kLength = 17U;
  std::vector<std::uint32_t> values;
  std::array<std::uint32_t, 2> flags{};
  std::uint32_t changing = 0U;
  BOOST_REQUIRE((kernels().sweepChain(*device.runtime, device.stream, *device.arena, kLength, kLength + 5U, kLength - 1U, values, flags, changing) ==
                 StageStatus::kComplete));
  for (std::uint32_t index = 0; index < kLength; ++index) {
    BOOST_TEST(values[index] == index + 1U);
  }
  BOOST_TEST(flags[0] == 0U);
  BOOST_TEST(flags[1] == 0U);
  BOOST_TEST(changing >= 1U);
  BOOST_TEST(changing <= kLength - 1U);

  // A value above the threshold, early in the chain: the second flag stays up.
  BOOST_REQUIRE(
    (kernels().sweepChain(*device.runtime, device.stream, *device.arena, kLength, 3U, kLength - 1U, values, flags, changing) == StageStatus::kComplete));
  BOOST_TEST(flags[0] == 0U);
  BOOST_TEST(flags[1] == 1U);

  // No changing pass allowed: a chain that needs one throws.
  BOOST_CHECK_THROW(static_cast<void>(kernels().sweepChain(*device.runtime, device.stream, *device.arena, kLength, kLength + 5U, 0U, values, flags, changing)),
                    std::logic_error);
  // A chain of one value never changes: no changing pass needed.
  BOOST_REQUIRE((kernels().sweepChain(*device.runtime, device.stream, *device.arena, 1U, 5U, 0U, values, flags, changing) == StageStatus::kComplete));
  BOOST_TEST(changing == 0U);
}

BOOST_AUTO_TEST_CASE(CopyColumnsAtOffsets)
{
  Device device;
  std::vector<Tracklet> records(1000U);
  for (std::uint32_t i = 0; i < records.size(); ++i) {
    records[i] = Tracklet{.source = i, .target = 2U * i, .link = static_cast<std::uint16_t>(i % 7U), .lambda = 0.001F * static_cast<float>(i), .phi = -0.5F};
  }
  const auto source = recordsOf<TrackletT>(records);
  std::vector<Tracklet> destination;
  BOOST_REQUIRE((kernels().copyColumns(*device.runtime, device.stream, *device.arena, *source, 250U, 600U, 100U, 400U, destination) == StageStatus::kComplete));
  BOOST_REQUIRE(destination.size() == 600U);
  const Tracklet untouched{.source = ~0U, .target = ~0U, .link = 0xFFFFU, .lambda = std::bit_cast<float>(~0U), .phi = std::bit_cast<float>(~0U)};
  for (std::size_t row = 0; row < destination.size(); ++row) {
    const bool copied = row >= 100U && row < 500U;
    BOOST_TEST(sameTracklet(destination[row], copied ? records[row - 100U + 250U] : untouched), "row " << row);
  }
  // Past either end is refused.
  BOOST_CHECK_THROW(static_cast<void>(kernels().copyColumns(*device.runtime, device.stream, *device.arena, *source, 900U, 600U, 0U, 200U, destination)),
                    std::out_of_range);
  BOOST_CHECK_THROW(static_cast<void>(kernels().copyColumns(*device.runtime, device.stream, *device.arena, *source, 0U, 600U, 500U, 200U, destination)),
                    std::out_of_range);
}

/// Chunks of several sizes, one empty, gathered back in order; an arena too
/// small refuses a chunk and is left as it was.
BOOST_AUTO_TEST_CASE(ColumnChunksGatherInOrder)
{
  std::vector<Tracklet> records(1000U);
  for (std::uint32_t i = 0; i < records.size(); ++i) {
    records[i] = Tracklet{.source = i, .target = 3U * i, .link = static_cast<std::uint16_t>(i % 11U), .lambda = 0.01F * static_cast<float>(i), .phi = 1.F};
  }
  const auto source = recordsOf<TrackletT>(records);
  const std::vector<std::size_t> chunkRows{100U, 0U, 350U, 1U, 549U};
  {
    Device device;
    std::vector<Tracklet> gathered;
    BOOST_REQUIRE(kernels().gatherChunks(*device.runtime, device.stream, *device.arena, *source, chunkRows, gathered));
    BOOST_REQUIRE(gathered.size() == records.size());
    for (std::size_t row = 0; row < records.size(); ++row) {
      BOOST_TEST(sameTracklet(gathered[row], records[row]), "row " << row);
    }
    BOOST_TEST(device.arena->used() == 0U);
  }
  {
    Device small{std::size_t{40} * 1024U};
    std::vector<Tracklet> gathered;
    BOOST_TEST(!kernels().gatherChunks(*small.runtime, small.stream, *small.arena, *source, chunkRows, gathered));
    BOOST_TEST(small.arena->used() == 0U);
  }
}

/// A flat object lands usable: the device reads its values through the
/// pointer relocated for the device buffer, also where a second upload of
/// the kept image lands elsewhere. An arena too small refuses it.
BOOST_AUTO_TEST_CASE(AFlatObjectIsUsableOnTheDevice)
{
  std::mt19937 random{9U};
  std::vector<std::int32_t> values(10000U);
  std::ranges::generate(values, [&] { return static_cast<std::int32_t>(random()); });
  ToyFlatObject object;
  object.construct(values);
  {
    Device device;
    std::vector<std::int32_t> readBack;
    std::vector<std::int32_t> moved;
    BOOST_REQUIRE(kernels().flatObject(*device.runtime, device.stream, *device.arena, object, readBack, moved));
    BOOST_TEST(readBack == values, boost::test_tools::per_element());
    BOOST_TEST(moved == values, boost::test_tools::per_element());
    BOOST_TEST(device.arena->used() == 0U);
  }
  {
    Device small{std::size_t{4} * 1024U};
    std::vector<std::int32_t> readBack;
    std::vector<std::int32_t> moved;
    BOOST_TEST(!kernels().flatObject(*small.runtime, small.stream, *small.arena, object, readBack, moved));
    BOOST_TEST(small.arena->used() == 0U);
  }
}

} // namespace o2::gtrk
