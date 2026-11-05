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

/// \file testHostAlgorithms.cxx
/// \brief The host's ordered compaction and exclusive scan against the serial
///        loops they replace, at one and eight threads and several chunk sizes.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK HostAlgorithms
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Base/HostAlgorithms.h"
#include "GenericTracking/Tracking/Parallel.h"

#include "SerialReferences.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <limits>
#include <memory_resource>
#include <vector>

namespace o2::gtrk
{

namespace
{

using host_algorithms::Index;

/// body(i) for i in [0, count) on \a arena's threads, one index per task.
struct ArenaForEach {
  detail::TaskArena* arena = nullptr;
  template <typename Body>
  void operator()(Index count, const Body& body) const
  {
    detail::parallelFor(*arena, 0, count, 1, [&](std::size_t first, std::size_t last) {
      for (std::size_t i = first; i < last; ++i) {
        body(static_cast<Index>(i));
      }
    });
  }
};

constexpr std::array<Index, 6> kSizes{0U, 1U, 5U, 256U, 1000U, 100003U};
constexpr std::array<Index, 4> kChunkSizes{1U, 3U, 256U, 100000U};

/// A pseudo-random bit per index, about one in three set.
bool scattered(Index i) { return ((i * 2654435761U) >> 7U) % 3U == 0U; }

template <typename ForEach>
void checkCompaction(const ForEach& forEach)
{
  const std::array<bool (*)(Index), 3> predicates{[](Index) { return true; }, [](Index) { return false; }, scattered};
  for (const Index n : kSizes) {
    for (const Index chunkSize : kChunkSizes) {
      for (const auto keep : predicates) {
        std::vector<Index> expected;
        for (Index i = 0; i < n; ++i) {
          if (keep(i)) {
            expected.push_back(i);
          }
        }
        BoundedVector<Index> out{std::pmr::new_delete_resource()};
        out.assign(7U, 99U); // stale contents and size are replaced
        host_algorithms::compactIndices(forEach, n, chunkSize, keep, out);
        BOOST_TEST_CONTEXT("n " << n << ", chunks of " << chunkSize)
        {
          BOOST_TEST(std::vector<Index>(out.begin(), out.end()) == expected, boost::test_tools::per_element());
        }
      }
    }
  }
}

template <typename ForEach>
void checkScan(const ForEach& forEach)
{
  const auto value = [](Index i) { return (i * 7U) % 11U; };
  for (const Index n : kSizes) {
    for (const Index chunkSize : kChunkSizes) {
      std::vector<Index> expected(std::size_t{n} + 1U);
      Index running = 0U;
      for (Index i = 0; i < n; ++i) {
        expected[i] = running;
        running += value(i);
      }
      expected[n] = running;
      std::vector<Index> offsets(std::size_t{n} + 1U, 12345U);
      BOOST_TEST_CONTEXT("n " << n << ", chunks of " << chunkSize)
      {
        BOOST_TEST(host_algorithms::exclusiveScan(forEach, n, chunkSize, value, offsets, std::pmr::new_delete_resource()));
        BOOST_TEST(offsets == expected, boost::test_tools::per_element());
      }
    }
  }
  // A total past the index's range fails, as the serial loop's running check.
  std::vector<Index> offsets(4U);
  const auto large = [](Index) { return std::numeric_limits<Index>::max() / 2U; };
  BOOST_TEST(host_algorithms::exclusiveScan(forEach, 2U, 1U, large, offsets, std::pmr::new_delete_resource()));
  BOOST_TEST(!host_algorithms::exclusiveScan(forEach, 3U, 1U, large, offsets, std::pmr::new_delete_resource()));
}

} // namespace

/// Each index of [begin, end) runs exactly once, in grains of at most the
/// grain size starting at begin plus a multiple of it, with the calling
/// thread one of the workers.
BOOST_AUTO_TEST_CASE(ParallelForRunsEachIndexOnce)
{
  for (const int threads : {1, 3, 8}) {
    detail::TaskArena arena{threads};
    for (const std::size_t begin : {std::size_t{0}, std::size_t{7}}) {
      for (const std::size_t count : {std::size_t{0}, std::size_t{1}, std::size_t{5}, std::size_t{1000}, std::size_t{100003}}) {
        for (const std::size_t grain : {std::size_t{1}, std::size_t{3}, std::size_t{64}, std::size_t{200000}}) {
          std::vector<std::atomic<int>> runs(count);
          std::atomic<int> badGrains{0};
          detail::parallelFor(arena, begin, begin + count, grain, [&](std::size_t first, std::size_t last) {
            if (first >= last || last - first > grain || (first - begin) % grain != 0U) {
              ++badGrains;
            }
            for (std::size_t i = first; i < last; ++i) {
              ++runs[i - begin];
            }
          });
          BOOST_TEST(badGrains.load() == 0);
          BOOST_TEST(std::ranges::all_of(runs, [](const std::atomic<int>& run) { return run.load() == 1; }));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(CompactionIsTheSerialLoops)
{
  checkCompaction(SerialForEach{});
  detail::TaskArena one{1};
  checkCompaction(ArenaForEach{&one});
  detail::TaskArena eight{8};
  checkCompaction(ArenaForEach{&eight});
}

BOOST_AUTO_TEST_CASE(ScanIsTheSerialLoops)
{
  checkScan(SerialForEach{});
  detail::TaskArena one{1};
  checkScan(ArenaForEach{&one});
  detail::TaskArena eight{8};
  checkScan(ArenaForEach{&eight});
}

} // namespace o2::gtrk
