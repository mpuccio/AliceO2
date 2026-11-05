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
/// \brief Profile activation, nested zones, pool peaks and stage loops.
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK tracking profile
#include <boost/test/unit_test.hpp>
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/Profile.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "SyntheticStage.h"
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
namespace o2::gtrk
{
using profile::Phase;
BOOST_AUTO_TEST_CASE(NothingWithoutActivation)
{
  BOOST_TEST(profile::current() == nullptr);
  const profile::Zone zone{Phase::kRun};
  profile::count(Phase::kRun, 1U, 1U);
  BOOST_TEST(profile::current() == nullptr);
}

BOOST_AUTO_TEST_CASE(ZonesNest)
{
  profile::Profile result;
  {
    const profile::Activation activation{result, nullptr};
    const profile::Zone run{Phase::kRun};
    for (int repeat = 0; repeat < 3; ++repeat) {
      const profile::Zone begin{Phase::kBegin};
      BOOST_TEST((profile::current()->phase == Phase::kBegin));
    }
    BOOST_TEST((profile::current()->phase == Phase::kRun));
    {
      // An inner activation profiles into its own profile, and restores.
      profile::Profile inner;
      const profile::Activation nested{inner, nullptr};
      const profile::Zone iteration{Phase::kIteration};
      BOOST_TEST(profile::current()->profile == &inner);
    }
    BOOST_TEST(profile::current()->profile == &result);
  }
  BOOST_TEST(profile::current() == nullptr);
  BOOST_TEST(result[Phase::kRun].calls == 1U);
  BOOST_TEST(result[Phase::kBegin].calls == 3U);
  BOOST_TEST(result[Phase::kIteration].calls == 0U);
  BOOST_TEST(result[Phase::kBegin].seconds <= result[Phase::kRun].seconds);
  BOOST_TEST((profile::info(Phase::kBegin).parent == Phase::kRun));
}

BOOST_AUTO_TEST_CASE(ZonesSeeThePoolsPeakInTheirWindow)
{
  BoundedMemoryResource pool;
  void* held = pool.allocate(100U);
  profile::Profile result;
  {
    const profile::Activation activation{result, &pool};
    const profile::Zone run{Phase::kRun};
    {
      // A phase's peak counts what the pool held when it began.
      const profile::Zone begin{Phase::kBegin};
      void* bytes = pool.allocate(1000U);
      pool.deallocate(bytes, 1000U);
    }
    BOOST_TEST(pool.getPeakMemory() == 1100U);
    {
      // A later, smaller phase sees its own peak, not the earlier one's.
      const profile::Zone iteration{Phase::kIteration};
      void* bytes = pool.allocate(10U);
      pool.deallocate(bytes, 10U);
    }
    // Closing a window keeps the larger of the two.
    BOOST_TEST(pool.getPeakMemory() == 1100U);
  }
  pool.deallocate(held, 100U);
  BOOST_TEST(result[Phase::kBegin].peakBytes == 1100U);
  BOOST_TEST(result[Phase::kIteration].peakBytes == 110U);
  BOOST_TEST(result[Phase::kRun].peakBytes == 1100U);
  BOOST_TEST(result[Phase::kBegin].peakGrowthBytes == 1000U);
  BOOST_TEST(result[Phase::kIteration].peakGrowthBytes == 10U);
  BOOST_TEST(result[Phase::kRun].peakGrowthBytes == 1000U);
  BOOST_TEST(pool.getPeakMemory() == 1100U);

  profile::Profile other;
  other[Phase::kRun].peakBytes = 2000U;
  result += other;
  BOOST_TEST(result[Phase::kRun].peakBytes == 2000U);
}

BOOST_AUTO_TEST_CASE(StageLoopsFoldIntoTheirZones)
{
  using Output = StageOutput<TrackletT>;
  constexpr std::uint32_t kElements = 1000U;
  constexpr std::uint32_t kChunk = 64U;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  Output output;
  output.setMemoryPool(pool);
  output.prepare(kElements);
  detail::TaskArena cpuArena{4};
  CpuStage<TrackletT> cpu{cpuArena};
  cpu.setMemoryPool(pool);
  cpu.setPhases({.reset = Phase::kTrackletReset, .generate = Phase::kTrackletGenerate, .scan = Phase::kTrackletScan, .place = Phase::kTrackletPlace});

  profile::Profile result;
  profile::Trace trace{4U};
  {
    const profile::Activation activation{result, pool.get(), &trace};
    const profile::Zone zone{Phase::kTracklets};
    BOOST_REQUIRE(cpu.run(output, kChunk, synthetic::Generator{}) == StageStatus::kComplete);
    const profile::Zone sort{Phase::kTrackletSort};
    cpu.forEachElement(kElements, kChunk, [](std::uint32_t) {});
  }

  const auto chunks = (kElements + kChunk - 1U) / kChunk;
  const auto& generate = result[Phase::kTrackletGenerate];
  BOOST_TEST(generate.calls == 1U);
  BOOST_TEST(generate.regions == 1U);
  BOOST_TEST(generate.elements == kElements);
  BOOST_TEST(generate.outputs == output.size());
  BOOST_TEST(generate.chunks == chunks);
  BOOST_TEST(generate.busySeconds > 0.);
  BOOST_TEST(generate.busySeconds <= generate.capacitySeconds);
  BOOST_TEST(generate.regionSeconds <= generate.seconds);
  const auto& place = result[Phase::kTrackletPlace];
  BOOST_TEST(place.elements == output.size());
  // The scan's loop over the chunks is a region of its phase too.
  const auto& scan = result[Phase::kTrackletScan];
  BOOST_TEST(scan.regions == 1U);
  BOOST_TEST(scan.elements == kElements);
  BOOST_TEST(scan.chunks == chunks);
  BOOST_TEST(result[Phase::kTrackletReset].calls == 1U);
  const auto& sort = result[Phase::kTrackletSort];
  BOOST_TEST(sort.elements == kElements);
  BOOST_TEST(sort.chunks == chunks);
  // The loops were in their phases, not in the zone around them.
  BOOST_TEST(result[Phase::kTracklets].regions == 0U);
  const auto& tracklets = result[Phase::kTracklets];
  BOOST_TEST(tracklets.seconds >= generate.seconds + sort.seconds);

  std::ostringstream timeline;
  trace.write(timeline);
  BOOST_TEST(timeline.str().find(R"("name":"tracklet generate","cat":"chunk")") != std::string::npos);
  BOOST_TEST(timeline.str().find(R"("name":"tracklets","cat":"zone")") != std::string::npos);
}
} // namespace o2::gtrk
