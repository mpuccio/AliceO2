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

/// \file testSelection.cxx
/// \brief The global selection's rounds of bidding against hand-made cases
///        and the sequential greedy walk, with and without shared clusters.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Selection
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/Selection.h"
#include "SerialReferences.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <span>
#include <vector>

namespace o2::gtrk
{

namespace
{

constexpr std::size_t kCapacity = 8U;
using Selection = GlobalSelection<kCapacity>;
using Records = Selection::Records;
using Record = BoundedColumns<SelectionRecordOf<kCapacity>::template RecordT>::Record;

struct Road {
  std::vector<std::uint32_t> rows; ///< innermost first
  Scalar chi2 = 0.F;
  std::uint8_t kinks = 0U; ///< the outward pass's
};

Record record(const Road& road)
{
  Record result;
  for (std::size_t i = 0; i < kCapacity; ++i) {
    result.clusters.index[i] = i < road.rows.size() ? road.rows[i] : ClusterList<kCapacity>::kNoCluster;
  }
  result.count = static_cast<std::uint8_t>(road.rows.size());
  result.chi2 = road.chi2;
  result.kinks = road.kinks;
  return result;
}

/// One slice's records.
void append(Selection& selection, const std::vector<Road>& roads)
{
  Records records;
  records.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  records.prepare(static_cast<Records::Index>(roads.size()));
  BOOST_REQUIRE(runSerial(records, [&](Records::Index element, auto& emitter) { emitter.emit(record(roads[element])); }) == StageStatus::kComplete);
  BOOST_REQUIRE(selection.append(records, SerialForEach{}) == StageStatus::kComplete);
}

const auto kSerial = [](std::uint32_t elements, const auto& body) {
  for (std::uint32_t element = 0; element < elements; ++element) {
    body(element);
  }
};

const auto kSerialSort = [](std::span<std::uint32_t> indices, const auto& less) { std::ranges::sort(indices, less); };

/// The selected roads, each as its rows.
std::vector<std::vector<std::uint32_t>> selectedRows(const Selection& selection)
{
  std::vector<std::vector<std::uint32_t>> result;
  const auto records = selection.records();
  for (const auto index : selection.selected()) {
    result.emplace_back(records.clusters[index].index, records.clusters[index].index + records.count[index]);
  }
  return result;
}

struct Outcome {
  std::vector<std::vector<std::uint32_t>> selected;
  std::vector<std::uint8_t> mask;
  std::uint32_t rounds = 0U;
};

Outcome run(const std::vector<std::vector<Road>>& slices, std::uint32_t sharedInnermost, std::vector<std::uint8_t> mask)
{
  Selection selection;
  selection.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  for (const auto& slice : slices) {
    append(selection, slice);
  }
  BOOST_REQUIRE(selection.select(kSerial, kSerialSort, sharedInnermost, mask) == StageStatus::kComplete);
  BOOST_TEST(selection.mismatches() == 0U);
  return {.selected = selectedRows(selection), .mask = mask, .rounds = selection.rounds()};
}

} // namespace

BOOST_AUTO_TEST_CASE(TheBestRecordWinsAndItsLoserFreesTheNext)
{
  // A beats B on row 3; B would beat C on row 6, but B is out, so C is in:
  // the greedy result, in two rounds.
  const std::vector<Road> roads{{.rows = {6, 7, 8, 9}, .chi2 = 3.F}, {.rows = {3, 4, 5, 6}, .chi2 = 2.F}, {.rows = {0, 1, 2, 3}, .chi2 = 1.F}};
  const auto outcome = run({roads}, 0U, std::vector<std::uint8_t>(12U, kFree));
  BOOST_REQUIRE(outcome.selected.size() == 2U);
  BOOST_TEST(outcome.selected[0] == roads[2].rows, boost::test_tools::per_element());
  BOOST_TEST(outcome.selected[1] == roads[0].rows, boost::test_tools::per_element());
  BOOST_TEST(outcome.rounds == 2U);
  const std::vector<std::uint8_t> mask{2, 2, 2, 2, 0, 0, 2, 2, 2, 2, 0, 0};
  BOOST_TEST(outcome.mask == mask, boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(MoreClustersBeatALowerChi2)
{
  const std::vector<Road> roads{{.rows = {0, 1, 2, 3}, .chi2 = 1.F}, {.rows = {3, 4, 5, 6, 7}, .chi2 = 10.F}};
  const auto outcome = run({roads}, 0U, std::vector<std::uint8_t>(8U, kFree));
  BOOST_REQUIRE(outcome.selected.size() == 1U);
  BOOST_TEST(outcome.selected[0] == roads[1].rows, boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(ARoadFoundByTwoSlicesIsOneRecord)
{
  const Road road{.rows = {0, 1, 2, 3}, .chi2 = 1.F};
  Selection selection;
  selection.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  append(selection, {road});
  append(selection, {road});
  std::vector<std::uint8_t> mask(4U, kFree);
  BOOST_REQUIRE(selection.select(kSerial, kSerialSort, 0U, mask) == StageStatus::kComplete);
  BOOST_TEST(selection.duplicates() == 1U);
  BOOST_TEST(selection.mismatches() == 0U);
  BOOST_TEST(selection.selected().size() == 1U);

  // The same rows with another chi2: a road that depends on its slice.
  append(selection, {{.rows = {0, 1, 2, 3}, .chi2 = 2.F}});
  std::ranges::fill(mask, kFree);
  BOOST_REQUIRE(selection.select(kSerial, kSerialSort, 0U, mask) == StageStatus::kComplete);
  BOOST_TEST(selection.mismatches() == 1U);
}

BOOST_AUTO_TEST_CASE(TheInnermostClustersAreSharedWhenConfigured)
{
  // A and B share their first cluster; C has it third.
  const std::vector<Road> roads{{.rows = {0, 1, 2, 3}, .chi2 = 1.F}, {.rows = {0, 4, 5, 6}, .chi2 = 2.F}, {.rows = {7, 8, 0, 9}, .chi2 = 3.F}};
  const auto strict = run({roads}, 0U, std::vector<std::uint8_t>(10U, kFree));
  BOOST_TEST(strict.selected.size() == 1U);
  const auto shared = run({roads}, 1U, std::vector<std::uint8_t>(10U, kFree));
  BOOST_REQUIRE(shared.selected.size() == 2U);
  BOOST_TEST(shared.selected[1] == roads[1].rows, boost::test_tools::per_element());
  BOOST_TEST(shared.mask[0] == kShared);
  BOOST_TEST(shared.mask[1] == kExclusive);
  // With three shareable, C's third cluster may be shared too.
  BOOST_TEST(run({roads}, 3U, std::vector<std::uint8_t>(10U, kFree)).selected.size() == 3U);
}

BOOST_AUTO_TEST_CASE(EarlierIterationsTracksAreRespected)
{
  const std::vector<Road> roads{{.rows = {0, 1, 2, 3}, .chi2 = 1.F}, {.rows = {4, 5, 6, 7}, .chi2 = 1.F}};
  std::vector<std::uint8_t> mask(8U, kFree);
  mask[2] = kExclusive; // the first road's third cluster is taken
  mask[4] = kShared;    // the second road's first is shared
  BOOST_TEST(run({roads}, 0U, mask).selected.empty());
  const auto shared = run({roads}, 1U, mask);
  BOOST_REQUIRE(shared.selected.size() == 1U);
  BOOST_TEST(shared.selected[0] == roads[1].rows, boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(TheRoundsAreTheSequentialGreedyWalk)
{
  // Dense conflicts: 3000 roads of 4 to 8 clusters on 1500 rows, duplicates
  // among them, some rows already taken; for k = 0, 1 and 2, serially and on
  // 1 to 4 threads, in the given order and shuffled over two slices.
  // A fixed seed: the outcome is compared, not a rate.
  std::mt19937_64 generator(20260926U); // NOLINT(cert-msc32-c,cert-msc51-cpp)
  constexpr std::uint32_t kRows = 1500U;
  std::vector<Road> roads;
  std::uniform_int_distribution<std::uint32_t> row(0U, kRows - 1U);
  std::uniform_int_distribution<std::uint32_t> length(4U, 8U);
  std::uniform_real_distribution<float> chi2(0.F, 20.F);
  for (int i = 0; i < 3000; ++i) {
    Road road{.rows = {}, .chi2 = chi2(generator)};
    const auto count = length(generator);
    while (road.rows.size() < count) {
      const auto candidate = row(generator);
      if (std::ranges::find(road.rows, candidate) == road.rows.end()) {
        road.rows.push_back(candidate);
      }
    }
    roads.push_back(road);
    if (i % 17 == 0) {
      roads.push_back(road);
    }
  }
  std::vector<std::uint8_t> initial(kRows, kFree);
  for (std::uint32_t r = 0; r < kRows; r += 23U) {
    initial[r] = r % 2U == 0U ? kShared : kExclusive;
  }

  for (const std::uint32_t sharedInnermost : {0U, 1U, 2U}) {
    BOOST_TEST_CONTEXT("k " << sharedInnermost)
    {
      Selection selection;
      selection.setMemoryPool(std::make_shared<BoundedMemoryResource>());
      append(selection, roads);
      auto mask = initial;
      BOOST_REQUIRE(selection.select(kSerial, kSerialSort, sharedInnermost, mask) == StageStatus::kComplete);
      BOOST_TEST(selection.duplicates() == 177U);
      BOOST_TEST(selection.mismatches() == 0U);
      auto oracleMask = initial;
      const auto oracle = selectSequentially<kCapacity>(selection.records(), selection.ranked(), sharedInnermost, oracleMask);
      const std::vector<std::uint32_t> selected(selection.selected().begin(), selection.selected().end());
      BOOST_TEST(selected == oracle, boost::test_tools::per_element());
      BOOST_TEST(mask == oracleMask, boost::test_tools::per_element());
      BOOST_TEST_MESSAGE("k " << sharedInnermost << ": " << selected.size() << " of " << selection.ranked().size() << " in " << selection.rounds() << " rounds");
      BOOST_TEST(selection.rounds() > 1U);

      for (const int threads : {1, 2, 4}) {
        for (const std::uint32_t chunk : {1U, 7U, 64U}) {
          detail::TaskArena cpuArena{threads};
          CpuStage<SelectionRecordOf<kCapacity>::template RecordT> cpu{cpuArena};
          cpu.setMemoryPool(std::make_shared<BoundedMemoryResource>());
          auto parallelMask = initial;
          BOOST_REQUIRE(selection.select([&](std::uint32_t elements, const auto& body) { cpu.forEachElement(elements, chunk, body); }, [&](std::span<std::uint32_t> indices, const auto& less) { cpu.sortIndices(indices, less); }, sharedInnermost, parallelMask) == StageStatus::kComplete);
          const std::vector<std::uint32_t> parallel(selection.selected().begin(), selection.selected().end());
          BOOST_TEST_CONTEXT("threads " << threads << ", chunk " << chunk)
          {
            BOOST_TEST(parallel == selected, boost::test_tools::per_element());
            BOOST_TEST(parallelMask == mask, boost::test_tools::per_element());
          }
        }
      }

      auto shuffled = roads;
      std::ranges::shuffle(shuffled, generator);
      const std::vector<Road> first(shuffled.begin(), shuffled.begin() + 1000);
      const std::vector<Road> second(shuffled.begin() + 1000, shuffled.end());
      const auto reordered = run({second, first}, sharedInnermost, initial);
      auto expected = selectedRows(selection);
      BOOST_TEST((reordered.selected == expected));
      BOOST_TEST(reordered.mask == mask, boost::test_tools::per_element());
    }
  }
}

} // namespace o2::gtrk
