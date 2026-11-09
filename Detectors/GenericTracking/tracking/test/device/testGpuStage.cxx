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

/// \file testGpuStage.cxx
/// \brief GpuStage against the serial reference: the device half of the
///        stages' output equality test.
///
/// Built once per backend. Every device result is downloaded into a host
/// StageOutput and compared with runSerial byte for byte: counts, offsets and
/// every column. Needs the device.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK GpuStage
#include <boost/test/unit_test.hpp>

#include "DeviceTestFixture.h"
#include "DeviceTestKernels.h"
#include "StageOutputComparison.h"
#include "SyntheticStage.h"

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "SerialReferences.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk
{
namespace
{

using namespace device_test;

using Output = StageOutput<TrackletT>;
using synthetic::Kind;

/// The host's serial reference for \a kind over \a elements elements.
std::unique_ptr<Output> serial(Kind kind, std::uint32_t elements, StageStatus expected = StageStatus::kComplete)
{
  auto output = std::make_unique<Output>();
  output->setMemoryPool(std::make_shared<BoundedMemoryResource>());
  output->prepare(elements);
  BOOST_REQUIRE((runSerial(*output, synthetic::Generator{.count = {.kind = kind}}) == expected));
  return output;
}

/// A device run of \a kind, downloaded into a host StageOutput.
struct DeviceRun {
  StageStatus status = StageStatus::kComplete;
  SyntheticRunInfo info;
  Output output;
};

std::unique_ptr<DeviceRun> onDevice(DeviceRuntime& device, DeviceArena& arena, Kind kind, std::uint32_t elements, std::size_t capacity)
{
  auto stream = device.createStream();
  auto run = std::make_unique<DeviceRun>();
  run->output.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  run->status = kernels().runSynthetic(device, stream, arena, kind, elements, capacity, run->output, run->info);
  stream.synchronize();
  return run;
}

} // namespace

BOOST_AUTO_TEST_CASE(TheSampleGivesTheCanonicalOutput)
{
  auto* device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  const auto run = onDevice(*device, arena, Kind::kSample, 5U, 16U);
  BOOST_REQUIRE((run->status == StageStatus::kComplete));
  const std::array<std::uint32_t, 6> offsets{0, 1, 3, 3, 5, 5};
  BOOST_TEST(std::ranges::equal(run->output.offsets(), offsets));
  BOOST_TEST(run->output.size() == 5U);
  const std::array<std::uint32_t, 5> targets{0, 1000, 1001, 3000, 3001};
  const auto view = std::as_const(run->output).view();
  for (std::size_t i = 0; i < targets.size(); ++i) {
    BOOST_TEST(view[i].target == targets[i]);
  }
  BOOST_TEST(std::bit_cast<std::uint32_t>(static_cast<Scalar>(view[2].lambda)) == std::bit_cast<std::uint32_t>(1.125F));
  BOOST_TEST(sameStageOutput(run->output, *serial(Kind::kSample, 5U)));
}

/// Thousands of threads race for staging slots; the placement is still the
/// serial loop's, record for record.
BOOST_AUTO_TEST_CASE(TheDeviceMatchesTheSerialReference)
{
  auto* device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  for (const Kind kind : {Kind::kUneven, Kind::kIrregular}) {
    const auto reference = serial(kind, 97U);
    const auto run = onDevice(*device, arena, kind, 97U, 4096U);
    BOOST_REQUIRE((run->status == StageStatus::kComplete));
    BOOST_TEST(!run->info.predictionMissed);
    BOOST_TEST(sameStageOutput(run->output, *reference), "kind " << static_cast<int>(kind));
  }

  constexpr std::uint32_t kLarge = 200000U;
  const auto reference = serial(Kind::kIrregular, kLarge);
  const auto run = onDevice(*device, arena, Kind::kIrregular, kLarge, std::size_t{kLarge} * 6U);
  BOOST_REQUIRE((run->status == StageStatus::kComplete));
  BOOST_TEST(run->output.size() == reference->size());
  BOOST_TEST(sameStageOutput(run->output, *reference));
}

/// Over a list -- every element that emits, and every third of the rest --
/// the stage gives the full run's output: the rest's counts are 0.
BOOST_AUTO_TEST_CASE(AListedRunGivesTheFullOutput)
{
  auto* device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  for (const Kind kind : {Kind::kUneven, Kind::kIrregular}) {
    constexpr std::uint32_t kElements = 2000U;
    const auto reference = serial(kind, kElements);
    auto stream = device->createStream();
    DeviceRun run;
    run.output.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    run.status = kernels().runSyntheticListed(*device, stream, arena, kind, kElements, std::size_t{kElements} * 12U, run.output, run.info);
    stream.synchronize();
    BOOST_REQUIRE((run.status == StageStatus::kComplete));
    BOOST_TEST(run.info.listed < kElements);
    BOOST_TEST(sameStageOutput(run.output, *reference), "kind " << static_cast<int>(kind));
  }
}

/// A staging capacity of one overflows, the retry runs at the exact
/// size and keeps it, and the output does not depend on the mispredict.
BOOST_AUTO_TEST_CASE(AMispredictRetriesOnceAtTheExactSize)
{
  auto* device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  const auto reference = serial(Kind::kUneven, 60U);
  const auto run = onDevice(*device, arena, Kind::kUneven, 60U, 1U);
  BOOST_REQUIRE((run->status == StageStatus::kComplete));
  BOOST_TEST(run->info.predictionMissed);
  BOOST_TEST(run->info.capacityAfter == reference->size());
  BOOST_TEST(sameStageOutput(run->output, *reference));
}

BOOST_AUTO_TEST_CASE(TooManyRecordsForOneElementFailAsOnTheHost)
{
  auto* device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  static_cast<void>(serial(Kind::kOverflow, 3U, StageStatus::kOrdinalOverflow));
  const auto run = onDevice(*device, arena, Kind::kOverflow, 3U, Output::kMaxPerElement + 8U);
  BOOST_TEST((run->status == StageStatus::kOrdinalOverflow));
}

/// An arena too small is a refusal, not a throw, and leaves the arena as the
/// run found it: counts and offsets on the base end, nothing on the scratch
/// end. Once the staging is refused, once the output after a staging that fit.
BOOST_AUTO_TEST_CASE(AnExhaustedArenaIsOutOfMemory)
{
  auto* device = runtime();
  // 1001 counts and 1001 offsets, 4004 bytes each, the second 256-aligned.
  constexpr std::size_t kCountsAndOffsets = 4096U + 4004U;
  for (const std::size_t bytes : {std::size_t{16} * 1024U, std::size_t{150} * 1024U}) {
    const SliceScope slices;
    DeviceArena arena{*device, block(bytes)};
    const auto run = onDevice(*device, arena, Kind::kIrregular, 1000U, 5000U);
    BOOST_TEST((run->status == StageStatus::kOutOfMemory), "arena " << bytes);
    BOOST_TEST(arena.used(ArenaEnd::kScratch) == 0U);
    BOOST_TEST(arena.used(ArenaEnd::kBase) == kCountsAndOffsets, "arena " << bytes);
    if (bytes > std::size_t{100} * 1024U) {
      BOOST_TEST(arena.peak() > std::size_t{100} * 1024U); // the staging did fit: the output was refused
    }
  }
}

BOOST_AUTO_TEST_CASE(AnEmptyStageIsAValidEmptyOutput)
{
  auto* device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  const auto run = onDevice(*device, arena, Kind::kIrregular, 0U, 16U);
  BOOST_REQUIRE((run->status == StageStatus::kComplete));
  BOOST_TEST(run->output.size() == 0U);
  BOOST_TEST(run->output.offsets().size() == 1U);
  BOOST_TEST(run->output.offsets().front() == 0U);
}

/// Staging and scan storage are gone when run returns; counts, offsets and
/// the placed output stay on the base end, in the caller's scope.
BOOST_AUTO_TEST_CASE(TheScratchIsReleasedAndTheOutputStays)
{
  auto* device = runtime();
  DeviceArena arena{*device, block(kArenaBytes)};
  const auto run = onDevice(*device, arena, Kind::kIrregular, 5000U, 40000U);
  BOOST_REQUIRE((run->status == StageStatus::kComplete));
  BOOST_TEST(run->info.scratchAfter == 0U);
  const std::size_t kept = arena.used(ArenaEnd::kBase);
  constexpr std::size_t kTrackletBytes = 4U + 4U + 2U + 4U + 4U;
  const std::size_t minimum = (std::size_t{2} * (std::size_t{5000} + 1U) * sizeof(std::uint32_t)) + (std::size_t{run->output.size()} * kTrackletBytes);
  BOOST_TEST(kept >= minimum);
  BOOST_TEST(arena.peak() > kept); // the staging was there, on the other end
}

} // namespace o2::gtrk
