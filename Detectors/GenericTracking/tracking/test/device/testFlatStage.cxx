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

/// \file testFlatStage.cxx
/// \brief The flat candidate list against the serial reference:
///        `flat::runFlat` over the
///        synthetic candidates, its output downloaded and compared byte for
///        byte with `StageOutput::runSerial` of the same walk and test --
///        counts, offsets and every column -- over whole runs, listed runs,
///        forced chunks, a source past the ordinal's reach and empty runs.
///        Built once per backend; needs the device.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK FlatStage
#include <boost/test/unit_test.hpp>

#include "DeviceTestFixture.h"
#include "DeviceTestKernels.h"
#include "StageOutputComparison.h"
#include "SyntheticStage.h"

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "SerialReferences.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

namespace o2::gtrk
{
namespace
{

using namespace device_test;
using Output = StageOutput<TrackletT>;
using synthetic::Kind;

/// The host's serial reference: the same candidates, the same test.
std::unique_ptr<Output> serial(Kind kind, std::uint32_t divisor, std::uint32_t elements, StageStatus expected = StageStatus::kComplete)
{
  auto output = std::make_unique<Output>();
  output->setMemoryPool(std::make_shared<BoundedMemoryResource>());
  output->prepare(elements);
  BOOST_REQUIRE((runSerial(*output, synthetic::FlatGenerator{.count = {.kind = kind}, .accept = {.divisor = divisor}}) == expected));
  return output;
}

struct DeviceRun {
  StageStatus status = StageStatus::kComplete;
  FlatRunInfo info;
  Output output;
};

std::unique_ptr<DeviceRun> onDevice(DeviceRuntime& device, DeviceArena& arena, Kind kind, std::uint32_t divisor, std::uint32_t elements, bool listed = false,
                                    std::size_t maxChunkPairs = 0U)
{
  auto stream = device.createStream();
  auto run = std::make_unique<DeviceRun>();
  run->output.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  run->status = kernels().runFlatSynthetic(device, stream, arena, kind, divisor, elements, listed, maxChunkPairs, run->output, run->info);
  stream.synchronize();
  return run;
}

} // namespace

/// With every candidate accepted the flat run is the staged run: the
/// canonical sample output, record for record.
BOOST_AUTO_TEST_CASE(TheSampleGivesTheCanonicalOutput)
{
  auto device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  const auto run = onDevice(*device, arena, Kind::kSample, 0U, 5U);
  BOOST_REQUIRE((run->status == StageStatus::kComplete));
  const std::array<std::uint32_t, 6> offsets{0, 1, 3, 3, 5, 5};
  BOOST_TEST(std::ranges::equal(run->output.offsets(), offsets));
  BOOST_TEST(run->output.size() == 5U);
  BOOST_TEST(run->info.report.candidates == 5U);
  BOOST_TEST(run->info.report.chunks == 1U);
  BOOST_TEST(run->info.report.largest == 2U);
  BOOST_TEST(sameStageOutput(run->output, *serial(Kind::kSample, 0U, 5U)));
}

/// The test refuses a third of the candidates: the records, their order and
/// the counts are the serial loop's, for each kind of element count.
BOOST_AUTO_TEST_CASE(TheDeviceMatchesTheSerialReference)
{
  auto device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  for (const Kind kind : {Kind::kSample, Kind::kUneven, Kind::kIrregular}) {
    for (const std::uint32_t divisor : {0U, 2U, 3U}) {
      const auto reference = serial(kind, divisor, 97U);
      const auto run = onDevice(*device, arena, kind, divisor, 97U);
      BOOST_REQUIRE((run->status == StageStatus::kComplete));
      BOOST_TEST(sameStageOutput(run->output, *reference), "kind " << static_cast<int>(kind) << " divisor " << divisor);
    }
  }
  constexpr std::uint32_t kLarge = 200000U;
  const auto reference = serial(Kind::kIrregular, 3U, kLarge);
  const auto run = onDevice(*device, arena, Kind::kIrregular, 3U, kLarge);
  BOOST_REQUIRE((run->status == StageStatus::kComplete));
  BOOST_TEST(run->output.size() == reference->size());
  BOOST_TEST(sameStageOutput(run->output, *reference));
  BOOST_TEST_MESSAGE(device->backend() << ": " << kLarge << " elements, " << run->info.report.candidates << " candidates, " << run->output.size()
                                       << " records, " << run->info.report.chunks << " chunk(s), " << run->info.report.readbacks << " readbacks");
}

/// Chunks of one, seven and a thousand pairs cut the list inside elements;
/// the chunks' outputs concatenate into the one-chunk run's.
BOOST_AUTO_TEST_CASE(ChunksConcatenate)
{
  auto device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  constexpr std::uint32_t kElements = 2000U;
  const auto reference = serial(Kind::kIrregular, 3U, kElements);
  const auto whole = onDevice(*device, arena, Kind::kIrregular, 3U, kElements);
  BOOST_REQUIRE((whole->status == StageStatus::kComplete));
  BOOST_TEST(whole->info.report.chunks == 1U);
  for (const std::size_t chunk : {std::size_t{1U}, std::size_t{7U}, std::size_t{1000U}}) {
    const auto run = onDevice(*device, arena, Kind::kIrregular, 3U, kElements, false, chunk);
    BOOST_REQUIRE((run->status == StageStatus::kComplete));
    BOOST_TEST(run->info.report.chunks > 1U);
    BOOST_TEST(run->info.report.candidates == whole->info.report.candidates);
    BOOST_TEST(sameStageOutput(run->output, *reference), "chunks of " << chunk);
    BOOST_TEST(sameStageOutput(run->output, whole->output), "chunks of " << chunk);
  }
}

/// Over a list -- every element that emits and every third of the rest --
/// the output is the full run's: the rest's counts are 0.
BOOST_AUTO_TEST_CASE(AListedRunGivesTheFullOutput)
{
  auto device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  constexpr std::uint32_t kElements = 2000U;
  for (const Kind kind : {Kind::kUneven, Kind::kIrregular}) {
    const auto reference = serial(kind, 3U, kElements);
    const auto run = onDevice(*device, arena, kind, 3U, kElements, true);
    BOOST_REQUIRE((run->status == StageStatus::kComplete));
    BOOST_TEST(run->info.listed < kElements);
    BOOST_TEST(sameStageOutput(run->output, *reference), "kind " << static_cast<int>(kind));
    const auto chunked = onDevice(*device, arena, kind, 3U, kElements, true, 5U);
    BOOST_REQUIRE((chunked->status == StageStatus::kComplete));
    BOOST_TEST(sameStageOutput(chunked->output, *reference), "kind " << static_cast<int>(kind) << " in chunks");
  }
}

/// An element with more candidates than an ordinal numbers is no concern of
/// the list; its accepted records are, as the staging's limit: past it the
/// run reports the overflow the serial reference reports.
BOOST_AUTO_TEST_CASE(AnOversizedElementIsCountedNotNumbered)
{
  auto device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  constexpr std::uint32_t kElements = 4U;
  const auto reference = serial(Kind::kOverflow, 3U, kElements);
  const auto run = onDevice(*device, arena, Kind::kOverflow, 3U, kElements);
  BOOST_REQUIRE((run->status == StageStatus::kComplete));
  BOOST_TEST(run->info.report.largest == Output::kMaxPerElement + 1U);
  BOOST_TEST(sameStageOutput(run->output, *reference));
  const auto chunked = onDevice(*device, arena, Kind::kOverflow, 3U, kElements, false, 4096U);
  BOOST_REQUIRE((chunked->status == StageStatus::kComplete));
  BOOST_TEST(chunked->info.report.chunks > 1U);
  BOOST_TEST(sameStageOutput(chunked->output, *reference));

  serial(Kind::kOverflow, 0U, kElements, StageStatus::kOrdinalOverflow);
  const auto overflow = onDevice(*device, arena, Kind::kOverflow, 0U, kElements);
  BOOST_TEST((overflow->status == StageStatus::kOrdinalOverflow));
}

/// No elements, and elements with no candidates: a complete, empty output.
BOOST_AUTO_TEST_CASE(EmptyRunsComplete)
{
  auto device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  const auto none = onDevice(*device, arena, Kind::kIrregular, 3U, 0U);
  BOOST_REQUIRE((none->status == StageStatus::kComplete));
  BOOST_TEST(none->output.size() == 0U);
  BOOST_TEST(none->info.report.candidates == 0U);
  // kUneven's elements past the 40th have no candidates: from the 50th on, none at all.
  const auto reference = serial(Kind::kIrregular, 3U, 11U);
  const auto run = onDevice(*device, arena, Kind::kIrregular, 3U, 11U);
  BOOST_REQUIRE((run->status == StageStatus::kComplete));
  BOOST_TEST(sameStageOutput(run->output, *reference));
  BOOST_TEST(run->output.counts()[0] == 0U); // element 0 has (0 * 7) mod 11 = 0 candidates
}

} // namespace o2::gtrk
