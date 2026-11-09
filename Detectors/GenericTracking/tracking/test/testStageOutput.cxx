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

/// \file testStageOutput.cxx
/// \brief Generate, scan and place: the serial reference and the CPU backend
///        produce the same canonical output.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK StageOutput
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "SyntheticStage.h"
#include "GenericTracking/Tracking/Tracklet.h"

#include "StageOutputComparison.h"
#include "SerialReferences.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <numeric>
#include <random>
#include <thread>
#include <utility>
#include <vector>

namespace o2::gtrk
{

namespace
{

using Output = StageOutput<TrackletT>;
using Cpu = CpuStage<TrackletT>;
using Index = Output::Index;

template <typename Counts>
auto generator(Counts counts)
{
  return synthetic::CountGenerator{counts};
}

/// Dense then empty, as barrel sources next to disc ones: the case a staging
/// split by element count gets wrong.
const auto kUneven = [](Index element) { return element < 40U ? 6U : 0U; };
/// Irregular, including zero.
const auto kIrregular = [](Index element) { return (element * 7U) % 11U; };

std::shared_ptr<BoundedMemoryResource> pool() { return std::make_shared<BoundedMemoryResource>(); }

struct LayeredDependencyGenerator {
  std::array<Index, 5> layerBoundaries{0U, 37U, 1000U, 1000U, 1127U};
  /// The layers by depth, as `layerDepthGroups` lists a layout's: layer 3
  /// (read by layer 1) and the empty layer 2 at depth 0, layer 1 (read by
  /// layer 0) at depth 1, layer 0 at depth 2.
  struct Groups {
    std::array<Index, 4> order{3U, 2U, 1U, 0U};
    std::array<std::size_t, 5> groups{0U, 2U, 3U, 4U, 4U};
    std::size_t count = 3U;
  } layerGroups{};
  std::span<Index> values;
  template <typename Emitter>
  void operator()(Index element, Emitter& emitter) const
  {
    const auto value = [&] {
      if (element >= 1000U) {
        return 1U;
      }
      if (element >= 37U) {
        return values[1000U + element % 127U] + 1U;
      }
      return values[37U + element % 963U] + 1U;
    }();
    values[element] = value;
    if (element % 5U != 0U) {
      for (Index ordinal = 0U; ordinal < 4U; ++ordinal) {
        emitter.emit(Tracklet{.source = element, .target = value, .link = 0U, .lambda = 0.F, .phi = static_cast<Scalar>(ordinal)});
      }
    }
  }
};

} // namespace

BOOST_AUTO_TEST_CASE(EveryColumnReadsBackThroughItsOwnField)
{
  // One sentinel per column, each unlike the others in every byte, so a view
  // whose pointers were bound in an order other than the declared one reads a
  // wrong value rather than a plausible one.
  const Tracklet sentinel{
    .source = 0xA1A1A1A1U,
    .target = 0xB2B2B2B2U,
    .link = 0xC3C3U,
    .lambda = std::bit_cast<Scalar>(0x3ED4D4D4U),
    .phi = std::bit_cast<Scalar>(0xBEE5E5E5U),
  };
  const auto memory = pool();
  BoundedColumns<TrackletT> columns;
  columns.setMemoryResource(memory.get());
  columns.resize(3U);
  for (std::size_t i = 0; i < columns.size(); ++i) {
    Tracklet record = sentinel;
    record.source += static_cast<std::uint32_t>(i);
    record.link = static_cast<std::uint16_t>(record.link + i);
    columns.write(i, record);
  }

  const TrackletView view = columns.view();
  const TrackletConstView constView = std::as_const(columns).view();
  for (std::size_t i = 0; i < columns.size(); ++i) {
    BOOST_TEST(view.source[i] == sentinel.source + i);
    BOOST_TEST(view.target[i] == sentinel.target);
    BOOST_TEST(view.link[i] == sentinel.link + i);
    BOOST_TEST(std::bit_cast<std::uint32_t>(view.lambda[i]) == std::bit_cast<std::uint32_t>(sentinel.lambda));
    BOOST_TEST(std::bit_cast<std::uint32_t>(view.phi[i]) == std::bit_cast<std::uint32_t>(sentinel.phi));

    BOOST_TEST(constView.source[i] == view.source[i]);
    BOOST_TEST(constView.target[i] == view.target[i]);
    BOOST_TEST(constView.link[i] == view.link[i]);
    BOOST_TEST(std::bit_cast<std::uint32_t>(constView.lambda[i]) == std::bit_cast<std::uint32_t>(view.lambda[i]));
    BOOST_TEST(std::bit_cast<std::uint32_t>(constView.phi[i]) == std::bit_cast<std::uint32_t>(view.phi[i]));

    const Tracklet back = columns.read(i);
    BOOST_TEST(back.source == view.source[i]);
    BOOST_TEST(back.link == view.link[i]);
    BOOST_TEST(std::bit_cast<std::uint32_t>(back.phi) == std::bit_cast<std::uint32_t>(sentinel.phi));
  }
}

BOOST_AUTO_TEST_CASE(TheSerialReferenceDefinesTheCanonicalOutput)
{
  Output output;
  output.setMemoryPool(pool());
  output.prepare(5U);
  const auto counts = [](Index element) { return std::array<Index, 5>{1, 2, 0, 2, 0}[element]; };
  BOOST_REQUIRE((runSerial(output, generator(counts)) == StageStatus::kComplete));

  const std::array<Index, 6> offsets{0, 1, 3, 3, 5, 5};
  BOOST_TEST(std::ranges::equal(output.offsets(), offsets));
  BOOST_TEST(output.size() == 5U);
  const std::array<std::uint32_t, 5> targets{0, 1000, 1001, 3000, 3001};
  const auto view = std::as_const(output).view();
  for (std::size_t i = 0; i < targets.size(); ++i) {
    BOOST_TEST(view[i].target == targets[i]);
  }
  BOOST_TEST(std::bit_cast<std::uint32_t>(static_cast<Scalar>(view[2].lambda)) == std::bit_cast<std::uint32_t>(1.125F));
}

BOOST_AUTO_TEST_CASE(TheCpuBackendMatchesTheSerialReference)
{
  constexpr Index kElements = 97;
  for (const bool uneven : {true, false}) {
    Output serial;
    serial.setMemoryPool(pool());
    serial.prepare(kElements);
    const StageStatus serialStatus = uneven ? runSerial(serial, generator(kUneven)) : runSerial(serial, generator(kIrregular));
    BOOST_REQUIRE((serialStatus == StageStatus::kComplete));

    for (const int threads : {1, 2, 4}) {
      for (const Index chunk : {1U, 3U, 8U, 64U, 1000U}) {
        Output parallel;
        const auto parallelPool = pool();
        parallel.setMemoryPool(parallelPool);
        parallel.prepare(kElements);
        detail::TaskArena cpuArena{threads};
        Cpu cpu{cpuArena};
        BOOST_TEST(cpu.maxConcurrency() == threads);
        cpu.setMemoryPool(parallelPool);
        const StageStatus status = uneven ? cpu.run(parallel, chunk, generator(kUneven)) : cpu.run(parallel, chunk, generator(kIrregular));
        BOOST_REQUIRE_MESSAGE(status == StageStatus::kComplete, "threads " << threads << ", chunk " << chunk << ", uneven " << uneven);
        BOOST_TEST(sameStageOutput(parallel, serial), "threads " << threads << ", chunk " << chunk << ", uneven " << uneven);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(AStagingOverflowGivesTheExactRetrySize)
{
  constexpr Index kElements = 60;
  const auto generate = generator(kUneven);

  Output generous;
  generous.setMemoryPool(pool());
  generous.prepare(kElements);
  BOOST_REQUIRE((runSerial(generous, generate) == StageStatus::kComplete));

  Output tight;
  tight.setMemoryPool(pool());
  tight.prepare(kElements);
  SerialStage<TrackletT> one{1U};
  BOOST_TEST((one.run(tight, generate) == StageStatus::kStagingOverflow));
  BOOST_TEST(one.stagedCount() == generous.size());

  const std::size_t exactCapacity = one.stagedCount();
  BOOST_REQUIRE((SerialStage<TrackletT>{exactCapacity}.run(tight, generate) == StageStatus::kComplete));
  BOOST_TEST(sameStageOutput(tight, generous));
}

BOOST_AUTO_TEST_CASE(TheCpuStagingPredictionAvoidsGrowthAndAMissStaysExact)
{
  constexpr Index kElements = 60U;
  constexpr std::size_t kRecords = 240U;

  Output predicted;
  const auto predictedPool = pool();
  predicted.setMemoryPool(predictedPool);
  predicted.prepare(kElements);
  detail::TaskArena preparedArena{1};
  Cpu prepared{preparedArena};
  prepared.setMemoryPool(predictedPool);
  prepared.prepareStagingPerThread(kRecords);
  BOOST_REQUIRE((prepared.run(predicted, 8U, generator(kUneven)) == StageStatus::kComplete));
  BOOST_TEST(!prepared.grewSinceReset());
  BOOST_TEST(prepared.stagingCapacity() == kRecords);

  Output missed;
  const auto missedPool = pool();
  missed.setMemoryPool(missedPool);
  missed.prepare(kElements);
  detail::TaskArena smallArena{1};
  Cpu small{smallArena};
  small.setMemoryPool(missedPool);
  small.prepareStagingPerThread(1U);
  BOOST_REQUIRE((small.run(missed, 8U, generator(kUneven)) == StageStatus::kComplete));
  BOOST_TEST(small.grewSinceReset());
  BOOST_TEST(sameStageOutput(missed, predicted));
}

BOOST_AUTO_TEST_CASE(ResetKeepsTheStorageAndReleaseReturnsIt)
{
  const auto memory = pool();
  Output output;
  output.setMemoryPool(memory);
  output.prepare(50U);
  BOOST_REQUIRE((runSerial(output, generator(kIrregular)) == StageStatus::kComplete));
  const auto used = memory->getUsedMemory();
  BOOST_TEST(used > 0U);

  resetOutput(output);
  BOOST_TEST(memory->getUsedMemory() == used);
  BOOST_TEST(output.size() == 0U);
  BOOST_TEST(std::ranges::all_of(output.counts(), [](Index count) { return count == 0U; }));

  output.clearStorage();
  BOOST_TEST(memory->getUsedMemory() == 0U);
}

BOOST_AUTO_TEST_CASE(TooManyRecordsForOneElementFailAlikeOnBothPaths)
{
  const auto count = [](Index element) { return element == 1U ? static_cast<Index>(Output::kMaxPerElement) + 1U : 1U; };
  const auto generate = generator(count);

  Output serial;
  serial.setMemoryPool(pool());
  serial.prepare(3U);
  BOOST_TEST((runSerial(serial, generate) == StageStatus::kOrdinalOverflow));

  Output parallel;
  const auto parallelPool = pool();
  parallel.setMemoryPool(parallelPool);
  parallel.prepare(3U);
  detail::TaskArena cpuArena;
  Cpu cpu{cpuArena};
  cpu.setMemoryPool(parallelPool);
  BOOST_TEST((cpu.run(parallel, 1U, generate) == StageStatus::kOrdinalOverflow));
}

BOOST_AUTO_TEST_CASE(AnExhaustedPoolIsReportedNotThrown)
{
  const auto small = std::make_shared<BoundedMemoryResource>(std::size_t{256} * 1024U);
  Output output;
  output.setMemoryPool(small);
  output.prepare(1000U);
  detail::TaskArena cpuArena;
  Cpu cpu{cpuArena};
  cpu.setMemoryPool(small);
  const auto many = [](Index /*element*/) { return 50U; };
  BOOST_TEST((cpu.run(output, 16U, generator(many)) == StageStatus::kOutOfMemory));
}

// A failed rebuild leaves the vector empty and valid: its owner reuses it, or
// destroys it once, and the pool's accounting stays exact.
// A column set grows all or nothing: when a later column's allocation is
// refused, the earlier ones go back to the old size, so `size()`, which reads
// the first, never promises rows another column lacks.
BOOST_AUTO_TEST_CASE(AColumnSetGrowsAllOrNothing)
{
  // 4 + 4 + 2 + 4 + 4 bytes a row: 100 rows of the first two columns fit,
  // the third does not.
  BoundedMemoryResource bounded{900U};
  {
    BoundedColumns<TrackletT> columns;
    columns.setMemoryResource(&bounded);
    BOOST_CHECK_THROW(columns.resize(100U), std::bad_alloc);
    BOOST_TEST(columns.size() == 0U);
    columns.resize(10U);
    BOOST_TEST(columns.size() == 10U);
    columns.write(9U, Tracklet{.source = 1U, .target = 2U, .link = 3U, .lambda = 0.5F, .phi = 0.25F});
    BOOST_TEST(columns.read(9U).link == 3U);
  }
  BOOST_TEST(bounded.getUsedMemory() == 0U);
}

// A stage's sort runs on its task arena's threads alone: one thread, the
// caller's, not the default arena's every core, whatever the tracker's
// thread count.
BOOST_AUTO_TEST_CASE(AStagesSortRunsOnItsArenaAlone)
{
  detail::TaskArena cpuArena{1};
  Cpu cpu{cpuArena};
  std::vector<std::uint32_t> indices(std::size_t{1} << 20U);
  std::iota(indices.begin(), indices.end(), 0U);
  std::ranges::shuffle(indices, std::mt19937{44U});
  const auto caller = std::this_thread::get_id();
  std::atomic<bool> elsewhere{false};
  cpu.sortIndices(indices, [&](std::uint32_t lhs, std::uint32_t rhs) {
    if (std::this_thread::get_id() != caller) {
      elsewhere.store(true, std::memory_order_relaxed);
    }
    return lhs < rhs;
  });
  BOOST_TEST(!elsewhere.load());
  BOOST_TEST(std::ranges::is_sorted(indices));

  // More threads give the same order.
  detail::TaskArena wideArena{4};
  Cpu wide{wideArena};
  std::ranges::shuffle(indices, std::mt19937{45U});
  wide.sortIndices(indices, [](std::uint32_t lhs, std::uint32_t rhs) { return lhs < rhs; });
  BOOST_TEST(std::ranges::is_sorted(indices));
}

BOOST_AUTO_TEST_CASE(AnEmptyStageIsAValidEmptyOutput)
{
  Output output;
  const auto outputPool = pool();
  output.setMemoryPool(outputPool);
  output.prepare(0U);
  bool visited = false;
  const auto visit = [&visited](Index, auto&) { visited = true; };
  BOOST_TEST((runSerial(output, visit) == StageStatus::kComplete));

  detail::TaskArena cpuArena;
  Cpu cpu{cpuArena};
  cpu.setMemoryPool(outputPool);
  BOOST_TEST((cpu.run(output, 8U, visit) == StageStatus::kComplete));
  BOOST_TEST(!visited);
  BOOST_TEST(output.size() == 0U);
  BOOST_TEST(output.offsets().size() == 1U);
  BOOST_TEST(output.offsets().front() == 0U);
}

} // namespace o2::gtrk

namespace o2::gtrk
{
BOOST_AUTO_TEST_CASE(LayeredGenerationHonoursDependenciesAndCanonicalPlacement)
{
  constexpr Index kElements = 1127U;
  for (const auto threads : {1, 8}) {
    for (const auto chunk : {1U, 7U, 256U, 1000U}) {
      detail::TaskArena cpuArena{threads};
      Cpu cpu{cpuArena};
      cpu.setMemoryPool(pool());
      cpu.prepareStagingPerThread(1U);
      Output output;
      output.setMemoryPool(pool());
      output.prepare(kElements);
      std::array<Index, kElements> values{};
      const LayeredDependencyGenerator generate{.values = values};
      BOOST_REQUIRE(cpu.run(output, chunk, generate) == StageStatus::kComplete);
      const auto view = output.view();
      Index placed = 0U;
      for (Index element = 0U; element < kElements; ++element) {
        const auto expected = [&] {
          if (element >= 1000U) {
            return 1U;
          }
          if (element >= 37U) {
            return 2U;
          }
          return 3U;
        }();
        BOOST_TEST(values[element] == expected);
        BOOST_TEST(output.offsets()[element] == placed);
        if (element % 5U != 0U) {
          for (Index ordinal = 0U; ordinal < 4U; ++ordinal) {
            BOOST_TEST(view.source[placed] == element);
            BOOST_TEST(view.target[placed] == expected);
            BOOST_TEST(view.phi[placed] == static_cast<Scalar>(ordinal));
            ++placed;
          }
        }
      }
      BOOST_TEST(output.size() == placed);
      BOOST_TEST(output.offsets().back() == placed);
    }
  }
}
} // namespace o2::gtrk
