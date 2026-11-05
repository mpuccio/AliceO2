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

/// \file DeviceTestFixture.h
/// \brief What every device test needs: the backend it is built for, its
///        test kernels, the device on O2's GPU route with arenas on it and a
///        stream, and host stage outputs made from record lists.
///
/// A device test is built once per backend with GTRK_TEST_HIP or without it,
/// for CUDA. Include it after <boost/test/unit_test.hpp>, with the test's
/// BOOST_TEST_MODULE defined, in the test's one source: it holds the test
/// process's device, a global fixture.

#ifndef GTRK_TEST_DEVICE_DEVICETESTFIXTURE_H_
#define GTRK_TEST_DEVICE_DEVICETESTFIXTURE_H_

#include "DeviceTestKernels.h"

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceContext.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStream.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "SerialReferences.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/tree/observer.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace o2::gtrk::device_test
{

/// That backend's test kernels.
inline const DeviceTestEntry& kernels()
{
#ifdef GTRK_TEST_HIP
  return *testKernelsHIP();
#else
  return *testKernelsCUDA();
#endif
}

inline constexpr std::size_t kArenaBytes = std::size_t{256} * 1024U * 1024U;
/// The test process's one block, of which every arena of a test case takes
/// a slice: what a case holds at once.
inline constexpr std::size_t kBlockBytes = std::size_t{1} << 30U;

/// The test process's device: the backend's first, on O2's GPU route.
struct Process final : public boost::unit_test::test_observer {
  Process()
  {
#ifdef GTRK_TEST_HIP
    constexpr const char* kBackend = "HIP";
#else
    constexpr const char* kBackend = "CUDA";
#endif
    context = std::make_unique<DeviceContext>(kBackend, 0, std::array{kBlockBytes});
    instance() = this;
    boost::unit_test::framework::register_observer(*this);
  }
  ~Process() override
  {
    boost::unit_test::framework::deregister_observer(*this);
    instance() = nullptr;
  }
  Process(const Process&) = delete;
  Process& operator=(const Process&) = delete;
  Process(Process&&) = delete;
  Process& operator=(Process&&) = delete;

  /// Every test case starts with the whole block.
  void test_unit_start(const boost::unit_test::test_unit& /*unit*/) override { used = 0U; }

  static Process*& instance()
  {
    static Process* process = nullptr;
    return process;
  }

  std::unique_ptr<DeviceContext> context;
  std::size_t used = 0U;
};
BOOST_TEST_GLOBAL_FIXTURE(Process);

/// The device's runtime.
inline DeviceRuntime* runtime() { return &Process::instance()->context->runtime(); }

/// The next \a bytes of the test case's block, 4 kB aligned.
inline DeviceBlock block(std::size_t bytes)
{
  constexpr std::size_t kAlignment = 4096U;
  auto& process = *Process::instance();
  const auto whole = process.context->block(0U);
  const std::size_t at = (process.used + kAlignment - 1U) / kAlignment * kAlignment;
  BOOST_REQUIRE_MESSAGE(at + bytes <= whole.bytes, "the test case's arenas exceed the test block");
  process.used = at + bytes;
  return {.data = static_cast<std::byte*>(whole.data) + at, .bytes = bytes};
}

/// Gives the slices taken inside it back when it closes.
class SliceScope
{
 public:
  SliceScope() : mUsed(Process::instance()->used) {}
  ~SliceScope() { Process::instance()->used = mUsed; }
  SliceScope(const SliceScope&) = delete;
  SliceScope& operator=(const SliceScope&) = delete;
  SliceScope(SliceScope&&) = delete;
  SliceScope& operator=(SliceScope&&) = delete;

 private:
  std::size_t mUsed;
};

/// An arena on the device and a stream.
struct Device {
  explicit Device(std::size_t arenaBytes = kArenaBytes)
    : runtime(device_test::runtime()), arena(std::make_unique<DeviceArena>(*runtime, block(arenaBytes))), stream(runtime->createStream())
  {
  }

  SliceScope slices;
  DeviceRuntime* runtime;
  std::unique_ptr<DeviceArena> arena;
  DeviceStream stream;
};

/// A host stage output of one record per element, made by the serial
/// reference.
template <template <template <class> class> class S, typename Record>
std::unique_ptr<StageOutput<S>> recordsOf(const std::vector<Record>& records)
{
  auto output = std::make_unique<StageOutput<S>>();
  output->setMemoryPool(std::make_shared<BoundedMemoryResource>());
  output->prepare(static_cast<std::uint32_t>(records.size()));
  BOOST_REQUIRE(runSerial(*output, [&](std::uint32_t element, auto& emitter) { emitter.emit(records[element]); }) == StageStatus::kComplete);
  return output;
}

} // namespace o2::gtrk::device_test

#endif // GTRK_TEST_DEVICE_DEVICETESTFIXTURE_H_
