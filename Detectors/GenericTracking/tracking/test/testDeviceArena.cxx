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

/// \file testDeviceArena.cxx
/// \brief The device arena and the device storage's allocation and liveness:
///        host code, run here on a host buffer as the block.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Device Arena
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/Tracklet.h"

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>

namespace o2::gtrk
{
namespace
{

/// The arena only keeps its runtime: this one throws if it is ever called.
class NoRuntime final : public DeviceRuntime
{
 public:
  [[nodiscard]] const char* backend() const noexcept final { return "none"; }
  [[nodiscard]] const DeviceProperties& properties() const noexcept final { return mProperties; }
  void copyToDevice(void* /*to*/, const void* /*from*/, std::size_t /*bytes*/, DeviceStream& /*stream*/) final { fail(); }
  void copyToHost(void* /*to*/, const void* /*from*/, std::size_t /*bytes*/, DeviceStream& /*stream*/) final { fail(); }
  void copyOnDevice(void* /*to*/, const void* /*from*/, std::size_t /*bytes*/, DeviceStream& /*stream*/) final { fail(); }
  void fill(void* /*to*/, std::uint8_t /*value*/, std::size_t /*bytes*/, DeviceStream& /*stream*/) final { fail(); }
  void registerHost(void* /*host*/, std::size_t /*bytes*/) final { fail(); }
  void unregisterHost(void* /*host*/) noexcept final {}

 protected:
  [[nodiscard]] void* createStreamHandle() final { fail(); }

 private:
  [[noreturn]] static void fail() { throw std::logic_error("the arena called its runtime"); }
  void destroyStreamHandle(void* /*handle*/) noexcept final {}
  void synchronizeStreamHandle(void* /*handle*/) final { fail(); }
  [[nodiscard]] std::unique_ptr<DeviceTracker<detectors::ITS, o2::base::MatLayerCylSet>> makeITSTracker(DeviceBlock /*block*/) final { fail(); }
  [[nodiscard]] std::unique_ptr<DeviceTracker<detectors::MFT, o2::base::MatLayerCylSet>> makeMFTTracker(DeviceBlock /*block*/) final { fail(); }
  [[nodiscard]] std::unique_ptr<DeviceTracker<detectors::TRKFT3, o2::base::MatLayerCylSet>> makeTRKFT3Tracker(DeviceBlock /*block*/) final { fail(); }

  DeviceProperties mProperties;
};

/// The runtime and host memory as the arena's block, aligned as the device's
/// blocks are.
struct Host {
  NoRuntime runtime;
  alignas(4096) std::array<std::byte, std::size_t{16} * 1024U> memory{};

  DeviceBlock block(std::size_t bytes)
  {
    BOOST_REQUIRE(bytes <= memory.size());
    return {.data = memory.data(), .bytes = bytes};
  }
};

} // namespace

BOOST_AUTO_TEST_CASE(ArenaAlignmentMetricsResetAndRefusal)
{
  Host host;
  DeviceArena arena{host.runtime, host.block(2048U)};
  const auto first = arena.allocate(7U, 16U);
  const auto second = arena.allocate(19U, 256U);
  BOOST_REQUIRE(first);
  BOOST_REQUIRE(second);
  BOOST_TEST(static_cast<std::byte*>(second.data) - static_cast<std::byte*>(first.data) == 256);
  BOOST_TEST(arena.used() == 275U);
  BOOST_TEST(arena.peak() == 275U);

  const auto used = arena.used();
  const auto peak = arena.peak();
  const auto refused = arena.allocate(arena.capacity());
  BOOST_TEST(!static_cast<bool>(refused));
  BOOST_TEST(arena.used() == used);
  BOOST_TEST(arena.peak() == peak);
  arena.reset();
  const auto again = arena.allocate(7U, 16U);
  BOOST_TEST(again.data == first.data);
  arena.resetPeak();
  BOOST_TEST(arena.peak() == arena.used());
}

BOOST_AUTO_TEST_CASE(ScopesInvalidateOnlyWhatTheyRewind)
{
  Host host;
  DeviceArena arena{host.runtime, host.block(4096U)};
  DeviceSpan<std::uint32_t> timeframe;
  BOOST_REQUIRE(timeframe.allocate(arena, 8U));
  DeviceSpan<std::uint32_t> slice;
  {
    auto scope = arena.scope();
    BOOST_REQUIRE(slice.allocate(arena, 8U));
    BOOST_CHECK_NO_THROW(static_cast<void>(slice.data()));
  }
  BOOST_CHECK_THROW(static_cast<void>(slice.data()), std::logic_error);
  BOOST_CHECK_NO_THROW(static_cast<void>(timeframe.data()));

  arena.reset();
  BOOST_CHECK_THROW(static_cast<void>(timeframe.data()), std::logic_error);
}

BOOST_AUTO_TEST_CASE(ColumnAllocationIsTransactional)
{
  Host host;
  DeviceArena arena{host.runtime, host.block(800U)};
  const auto used = arena.used();
  DeviceColumns<TrackletT> columns;
  BOOST_TEST(!columns.allocate(arena, 1U));
  BOOST_TEST(arena.used() == used);
}

/// A transactional allocation opens its own scope and, on success, hands its
/// columns to the scope around it: they outlive an inner scope opened later and
/// die with the one they were allocated in.
BOOST_AUTO_TEST_CASE(CommittedColumnsLiveAndDieWithTheirScope)
{
  Host host;
  DeviceArena arena{host.runtime, host.block(std::size_t{16} * 1024U)};
  DeviceColumns<TrackletT> columns;
  {
    auto slice = arena.scope();
    const auto before = arena.used();
    BOOST_REQUIRE(columns.allocate(arena, 4U));
    BOOST_TEST(arena.used() > before);
    BOOST_CHECK_NO_THROW(static_cast<void>(columns.view()));
    {
      auto inner = arena.scope();
      DeviceSpan<std::uint32_t> scratch;
      BOOST_REQUIRE(scratch.allocate(arena, 16U));
    }
    BOOST_CHECK_NO_THROW(static_cast<void>(columns.view()));
    BOOST_TEST(columns.size() == 4U);
  }
  BOOST_CHECK_THROW(static_cast<void>(columns.view()), std::logic_error);
}

/// The scratch end grows down from the top: both ends allocate without
/// touching each other, meet but never cross, and each end's scopes rewind
/// only that end.
BOOST_AUTO_TEST_CASE(TheScratchEndIsIndependentOfTheBase)
{
  Host host;
  constexpr std::size_t kCapacity = 4096U;
  DeviceArena arena{host.runtime, host.block(kCapacity)};
  DeviceSpan<std::uint8_t> kept;
  BOOST_REQUIRE(kept.allocate(arena, 1000U));
  // The first allocation sits at the arena's base (the block is 256-aligned).
  const auto* origin = reinterpret_cast<const std::byte*>(kept.data());
  BOOST_TEST(arena.used(ArenaEnd::kBase) == 1000U);
  BOOST_TEST(arena.used(ArenaEnd::kScratch) == 0U);

  DeviceSpan<std::byte> staging;
  {
    auto scratch = arena.scope(ArenaEnd::kScratch);
    BOOST_REQUIRE(staging.allocate(arena, 1000U, ArenaEnd::kScratch));
    const std::byte* at = staging.data();
    BOOST_TEST(at >= origin + arena.used(ArenaEnd::kBase));
    BOOST_TEST(at + 1000 <= origin + kCapacity);
    BOOST_TEST(std::bit_cast<std::uintptr_t>(at) % 256U == 0U);
    BOOST_TEST(arena.used() == arena.used(ArenaEnd::kBase) + arena.used(ArenaEnd::kScratch));

    // The ends meet: what does not fit between them is refused and changes nothing.
    const auto base = arena.used(ArenaEnd::kBase);
    const auto top = arena.used(ArenaEnd::kScratch);
    BOOST_TEST(!static_cast<bool>(arena.allocate(kCapacity - base - top + 1U, 1U, ArenaEnd::kScratch)));
    BOOST_TEST(!static_cast<bool>(arena.allocate(kCapacity - base - top + 1U, 1U, ArenaEnd::kBase)));
    BOOST_TEST(arena.used(ArenaEnd::kBase) == base);
    BOOST_TEST(arena.used(ArenaEnd::kScratch) == top);

    // A base-end scope opened and closed meanwhile leaves the scratch alone.
    {
      auto slice = arena.scope();
      BOOST_REQUIRE(arena.allocate(16U));
    }
    BOOST_TEST(arena.used(ArenaEnd::kScratch) == top);
    BOOST_CHECK_NO_THROW(static_cast<void>(staging.data()));
  }
  BOOST_TEST(arena.used(ArenaEnd::kScratch) == 0U);
  BOOST_CHECK_THROW(static_cast<void>(staging.data()), std::logic_error);
  BOOST_CHECK_NO_THROW(static_cast<void>(kept.data()));
  BOOST_TEST(arena.peak() >= 2000U);

  // Transactional columns on the scratch end, and a reset that clears both ends.
  DeviceColumns<TrackletT> columns;
  BOOST_REQUIRE(columns.allocate(arena, 8U, ArenaEnd::kScratch));
  BOOST_TEST(arena.used(ArenaEnd::kScratch) > 0U);
  BOOST_CHECK_NO_THROW(static_cast<void>(columns.view()));
  arena.reset();
  BOOST_TEST(arena.used() == 0U);
  BOOST_CHECK_THROW(static_cast<void>(columns.view()), std::logic_error);
}

} // namespace o2::gtrk
