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

/// \file testSelectionGlue.cxx
/// \brief An iteration's global selection on the device against
///        GlobalSelection::select on the host: the ranked and selected
///        records, the duplicates and
///        mismatches, the rounds and the mask, all equal.
///
/// Generated records with ties in every part of the quality order,
/// duplicates the same and differing, sharing on and off, and a mask that
/// arrives with an earlier iteration's tracks. Built once per backend. Needs
/// the device.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK selection glue
#include <boost/test/unit_test.hpp>

#include "DeviceTestFixture.h"
#include "DeviceTestKernels.h"

#include "GenericTracking/Tracking/Definitions.h"
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

using namespace device_test;
constexpr std::size_t kK = kToyBarrel8Clusters;
using Record = ToyBarrel8Records::Record;
using Selection = GlobalSelection<kK>;

const auto kSerial = [](std::uint32_t elements, const auto& body) {
  for (std::uint32_t element = 0; element < elements; ++element) {
    body(element);
  }
};

/// \a count random records over \a rows rows: 3 to K clusters, the chi2 from
/// four values and the kinks from two, so that every tie of the quality
/// order occurs; every seventh an exact duplicate of an earlier record, and
/// with \a differing every thirteenth one with the same clusters in another
/// order, which the host counts as a mismatch.
std::vector<Record> generate(std::uint32_t count, std::uint32_t rows, bool differing, std::uint32_t seed)
{
  std::mt19937 random{seed};
  std::uniform_int_distribution<std::uint32_t> row{0U, rows - 1U};
  std::uniform_int_distribution<std::uint32_t> length{3U, static_cast<std::uint32_t>(kK)};
  std::uniform_int_distribution<std::uint32_t> pick{0U, 3U};
  constexpr std::array<Scalar, 4> kChi2{0.F, 1.5F, 1.5F, 7.25F};
  std::vector<Record> records;
  records.reserve(count);
  for (std::uint32_t i = 0; i < count; ++i) {
    if (i > 0U && i % 7U == 0U) {
      std::uniform_int_distribution<std::uint32_t> earlier{0U, i - 1U};
      records.push_back(records[earlier(random)]);
      continue;
    }
    if (differing && i > 0U && i % 13U == 0U) {
      auto copy = records[i - 1U];
      std::reverse(copy.clusters.index, copy.clusters.index + copy.count);
      records.push_back(copy);
      continue;
    }
    Record record;
    const auto clusters = length(random);
    std::vector<std::uint32_t> chosen;
    while (chosen.size() < clusters) {
      const auto candidate = row(random);
      if (std::ranges::find(chosen, candidate) == chosen.end()) {
        chosen.push_back(candidate);
      }
    }
    for (std::size_t k = 0; k < kK; ++k) {
      record.clusters.index[k] = k < chosen.size() ? chosen[k] : ClusterList<kK>::kNoCluster;
    }
    record.count = static_cast<std::uint8_t>(clusters);
    record.chi2 = kChi2[pick(random)];
    record.kinks = static_cast<std::uint8_t>(pick(random) % 2U);
    records.push_back(record);
  }
  return records;
}

std::unique_ptr<ToyBarrel8Records> columnsOf(const std::vector<Record>& records, const std::shared_ptr<BoundedMemoryResource>& pool)
{
  auto columns = std::make_unique<ToyBarrel8Records>();
  columns->setMemoryResource(pool.get());
  columns->resize(records.size());
  const auto view = columns->view();
  for (std::size_t i = 0; i < records.size(); ++i) {
    view.clusters[i] = records[i].clusters;
    view.count[i] = records[i].count;
    view.chi2[i] = records[i].chi2;
    view.kinks[i] = records[i].kinks;
  }
  return columns;
}

/// The host's selection of the same records, into its own copy of the mask.
struct HostOutcome {
  std::vector<std::uint32_t> ranked;
  std::vector<std::uint32_t> selected;
  std::vector<std::uint8_t> mask;
  std::uint32_t duplicates = 0U;
  std::uint32_t mismatches = 0U;
  std::uint32_t rounds = 0U;
};

HostOutcome host(const std::vector<Record>& records, std::vector<std::uint8_t> mask, std::uint32_t sharedInnermost)
{
  Selection selection;
  selection.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  selection.clear();
  Selection::Records output;
  output.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  output.prepare(static_cast<std::uint32_t>(records.size()));
  BOOST_REQUIRE(runSerial(output, [&](std::uint32_t element, auto& emitter) { emitter.emit(records[element]); }) == StageStatus::kComplete);
  BOOST_REQUIRE(selection.append(output, SerialForEach{}) == StageStatus::kComplete);
  BOOST_REQUIRE(selection.select(kSerial, [](std::span<std::uint32_t> indices, const auto& less) { std::ranges::sort(indices, less); }, sharedInnermost, mask) == StageStatus::kComplete);
  return {
    .ranked = {selection.ranked().begin(), selection.ranked().end()},
    .selected = {selection.selected().begin(), selection.selected().end()},
    .mask = mask,
    .duplicates = selection.duplicates(),
    .mismatches = selection.mismatches(),
    .rounds = selection.rounds(),
  };
}

/// \a rows rows, every eleventh shared and every seventeenth exclusive: an
/// earlier iteration's tracks.
std::vector<std::uint8_t> earlierMask(std::uint32_t rows)
{
  std::vector<std::uint8_t> mask(rows, kFree);
  for (std::uint32_t row = 0; row < rows; ++row) {
    if (row % 17U == 0U) {
      mask[row] = kExclusive;
    } else if (row % 11U == 0U) {
      mask[row] = kShared;
    }
  }
  return mask;
}

void compare(const std::vector<Record>& records, const std::vector<std::uint8_t>& mask, std::uint32_t sharedInnermost)
{
  const auto expected = host(records, mask, sharedInnermost);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto columns = columnsOf(records, pool);
  Device device;
  SelectionGlueResult result;
  BOOST_REQUIRE((kernels().selectionGlue(*device.runtime, device.stream, *device.arena, *columns, mask, sharedInnermost, result) == StageStatus::kComplete));
  BOOST_TEST(result.ranked == expected.ranked, boost::test_tools::per_element());
  BOOST_TEST(result.selected == expected.selected, boost::test_tools::per_element());
  BOOST_TEST(result.mask == expected.mask, boost::test_tools::per_element());
  BOOST_TEST(result.duplicates == expected.duplicates);
  BOOST_TEST(result.mismatches == expected.mismatches);
  BOOST_TEST(result.rounds == expected.rounds);
  BOOST_TEST(result.scratchAfter == 0U);
  // The rounds end in the sequential walk's result.
  std::vector<std::uint8_t> walked = mask;
  const auto sequential = selectSequentially<kK>(std::as_const(*columns).view(), expected.ranked, sharedInnermost, walked);
  BOOST_TEST(result.selected == sequential, boost::test_tools::per_element());
  BOOST_TEST_MESSAGE(runtime()->backend() << ": " << records.size() << " records, " << result.ranked.size() << " unique, " << result.duplicates << " duplicates, " << result.mismatches << " mismatches, " << result.rounds << " rounds, " << result.selected.size() << " selected, " << result.readbacks << " readbacks");
}

} // namespace

BOOST_AUTO_TEST_CASE(TheDeviceSelectsAsTheHost)
{
  constexpr std::uint32_t kRows = 3000U;
  for (const std::uint32_t shared : {0U, 2U}) {
    for (const bool earlier : {false, true}) {
      BOOST_TEST_CONTEXT("sharedInnermost " << shared << (earlier ? ", with an earlier iteration's mask" : ", free mask"))
      {
        const auto mask = earlier ? earlierMask(kRows) : std::vector<std::uint8_t>(kRows, kFree);
        compare(generate(5000U, kRows, false, 7U + shared), mask, shared);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(DuplicatesThatDifferAreCountedAsTheHostCountsThem)
{
  constexpr std::uint32_t kRows = 500U;
  const auto records = generate(2000U, kRows, true, 23U);
  const auto expected = host(records, std::vector<std::uint8_t>(kRows, kFree), 0U);
  BOOST_TEST(expected.mismatches > 0U);
  compare(records, std::vector<std::uint8_t>(kRows, kFree), 0U);
}

BOOST_AUTO_TEST_CASE(NoRecordsSelectNothing) { compare({}, std::vector<std::uint8_t>(10U, kFree), 0U); }

/// A refusal of the arena returns kOutOfMemory and leaves its scratch end as
/// it was.
BOOST_AUTO_TEST_CASE(ARefusalLeavesTheArena)
{
  constexpr std::uint32_t kRows = 3000U;
  const auto records = generate(5000U, kRows, false, 3U);
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto columns = columnsOf(records, pool);
  // Room for the records and the mask, not for the keys: K words a record.
  const std::size_t recordBytes = records.size() * (sizeof(ClusterList<kK>) + sizeof(std::uint8_t) + sizeof(Scalar) + sizeof(std::uint8_t));
  Device device{recordBytes + kRows + (std::size_t{64} * 1024U)};
  SelectionGlueResult result;
  const std::vector<std::uint8_t> mask(kRows, kFree);
  BOOST_TEST((kernels().selectionGlue(*device.runtime, device.stream, *device.arena, *columns, mask, 0U, result) == StageStatus::kOutOfMemory));
  BOOST_TEST(result.scratchAfter == 0U);
  BOOST_TEST(device.arena->used(ArenaEnd::kBase) == 0U);
}

} // namespace o2::gtrk
