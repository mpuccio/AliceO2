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

/// \file testTripletEndIndex.cxx
/// \brief The index from a tracklet to the triplets ending with it, against a
///        brute-force scan of the triplet arena.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TripletEndIndex
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/TripletEndIndex.h"

#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/Neighbour.h"
#include "SerialReferences.h"

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <memory>
#include <random>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace o2::gtrk
{

namespace
{

using Triplets = StageOutput<TripletT>;
using Index = TripletEndIndex::Index;
constexpr Index kAnyChunks = std::numeric_limits<Index>::max();

/// The triplets whose second tracklet is \a tracklet, in arena order.
std::span<const Index> ending(const TripletEndIndex& index, Index tracklet) { return {index.triplets() + index.offsets()[tracklet], index.offsets()[tracklet + 1U] - index.offsets()[tracklet]}; }

/// A triplet arena whose element `e` emits the triplets `(first, second)` of
/// \a records[e].
void fill(Triplets& triplets, const std::vector<std::vector<std::pair<Index, Index>>>& records)
{
  triplets.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  triplets.prepare(static_cast<Index>(records.size()));
  const auto status = runSerial(triplets, [&](Index element, auto& emitter) {
    for (const auto& [first, second] : records[element]) {
      emitter.emit(Triplet{first, second, 0.F, 0.F, 1U, 0U});
    }
  });
  BOOST_REQUIRE(status == StageStatus::kComplete);
}

/// Every tracklet's range holds exactly the triplets ending with it, in
/// arena order.
void checkAgainstBruteForce(const TripletEndIndex& index, Index trackletCount, const Triplets& triplets)
{
  const auto view = triplets.view();
  BOOST_REQUIRE(index.offsets()[trackletCount] == triplets.size());
  for (Index tracklet = 0; tracklet < trackletCount; ++tracklet) {
    std::vector<Index> expected;
    for (Index triplet = 0; triplet < triplets.size(); ++triplet) {
      if (view.second[triplet] == tracklet) {
        expected.push_back(triplet);
      }
    }
    const auto range = ending(index, tracklet);
    BOOST_TEST_CONTEXT("tracklet " << tracklet)
    {
      BOOST_TEST(std::vector<Index>(range.begin(), range.end()) == expected, boost::test_tools::per_element());
      BOOST_TEST(index.offsets()[tracklet + 1U] - index.offsets()[tracklet] == range.size());
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(RangesHoldTheTripletsEndingWithEachTracklet)
{
  // Tracklet 4 ends three triplets, from two elements; 0, 2, 5 and 6 end none.
  Triplets triplets;
  fill(triplets, {{{0, 4}, {0, 3}}, {}, {{1, 4}, {2, 1}, {3, 4}}, {{5, 3}}});
  TripletEndIndex index;
  index.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  BOOST_REQUIRE(index.build(7U, triplets, SerialForEach{}, kAnyChunks) == StageStatus::kComplete);
  checkAgainstBruteForce(index, 7U, triplets);
  const auto four = ending(index, 4U);
  BOOST_TEST(std::vector<Index>(four.begin(), four.end()) == (std::vector<Index>{0U, 2U, 4U}), boost::test_tools::per_element());
  BOOST_TEST(ending(index, 6U).empty());
}

BOOST_AUTO_TEST_CASE(ALargerArenaAndARebuild)
{
  constexpr Index kTracklets = 97U;
  std::vector<std::vector<std::pair<Index, Index>>> records(40U);
  for (Index element = 0; element < records.size(); ++element) {
    for (Index k = 0; k < element % 7U; ++k) {
      records[element].emplace_back(element, ((element * 31U) + (k * 17U)) % kTracklets);
    }
  }
  Triplets triplets;
  fill(triplets, records);
  TripletEndIndex index;
  index.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  BOOST_REQUIRE(index.build(kTracklets, triplets, SerialForEach{}, kAnyChunks) == StageStatus::kComplete);
  checkAgainstBruteForce(index, kTracklets, triplets);

  // A second build over a smaller slice replaces the first.
  Triplets smaller;
  fill(smaller, {{{0, 2}}, {{1, 0}, {1, 2}}});
  BOOST_REQUIRE(index.build(3U, smaller, SerialForEach{}, kAnyChunks) == StageStatus::kComplete);
  checkAgainstBruteForce(index, 3U, smaller);
}

BOOST_AUTO_TEST_CASE(NoTriplets)
{
  Triplets triplets;
  fill(triplets, {{}, {}});
  TripletEndIndex index;
  index.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  BOOST_REQUIRE(index.build(5U, triplets, SerialForEach{}, kAnyChunks) == StageStatus::kComplete);
  BOOST_TEST(index.offsets()[5U] == 0U);
  checkAgainstBruteForce(index, 5U, triplets);
}

/// The two-level sort on the host threads gives the one-thread build's index
/// byte for byte, at 1, 3 and 8 threads, over arenas from empty to several
/// chunks and buckets, each checked against the brute force where that is
/// quick.
BOOST_AUTO_TEST_CASE(ParallelBuildEqualsSerialBuild)
{
  std::mt19937 engine{12345U};
  for (const auto& [elements, perElement, tracklets] : std::initializer_list<std::tuple<Index, Index, Index>>{{0U, 0U, 1U}, {1U, 1U, 1U}, {300U, 3U, 50U}, {700U, 4U, TripletEndIndex::kBucketTracklets}, {5000U, 7U, 9000U}, {40000U, 5U, 12345U}}) {
    std::vector<std::vector<std::pair<Index, Index>>> records(elements);
    std::uniform_int_distribution<Index> second{0U, tracklets - 1U};
    for (Index element = 0; element < elements; ++element) {
      const auto emitted = static_cast<Index>(engine() % (perElement + 1U));
      for (Index k = 0; k < emitted; ++k) {
        records[element].emplace_back(element, second(engine));
      }
    }
    Triplets triplets;
    fill(triplets, records);
    TripletEndIndex serial;
    serial.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    BOOST_REQUIRE(serial.build(tracklets, triplets, SerialForEach{}, kAnyChunks) == StageStatus::kComplete);
    if (static_cast<std::uint64_t>(tracklets) * triplets.size() < 20000000U) {
      checkAgainstBruteForce(serial, tracklets, triplets);
    }
    // The chunks bounded as the tracker bounds them, to one, or not at all.
    for (const auto& [threads, maxChunks] : std::initializer_list<std::pair<int, Index>>{{1, 0U}, {3, 0U}, {8, 0U}, {3, 1U}, {8, 64U}}) {
      BOOST_TEST_CONTEXT("elements " << elements << ", tracklets " << tracklets << ", threads " << threads << ", chunks at most " << maxChunks)
      {
        detail::TaskArena cpuArena{threads};
        CpuStage<NeighbourT> cpu{cpuArena};
        TripletEndIndex parallel;
        parallel.setMemoryPool(std::make_shared<BoundedMemoryResource>());
        const auto forEach = [&](Index count, const auto& body) { cpu.forEachElement(count, 1U, body); };
        BOOST_REQUIRE(parallel.build(tracklets, triplets, forEach, maxChunks == 0U ? kAnyChunks : maxChunks) == StageStatus::kComplete);
        BOOST_TEST(parallel.offsets()[tracklets] == serial.offsets()[tracklets]);
        BOOST_TEST(std::vector<Index>(parallel.offsets(), parallel.offsets() + tracklets + 1U) == std::vector<Index>(serial.offsets(), serial.offsets() + tracklets + 1U), boost::test_tools::per_element());
        BOOST_TEST(std::vector<Index>(parallel.triplets(), parallel.triplets() + parallel.offsets()[tracklets]) == std::vector<Index>(serial.triplets(), serial.triplets() + serial.offsets()[tracklets]), boost::test_tools::per_element());
      }
    }
  }
}

} // namespace o2::gtrk
