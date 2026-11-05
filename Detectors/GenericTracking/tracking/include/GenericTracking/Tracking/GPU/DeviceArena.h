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

/// \file DeviceArena.h
/// \brief A bounded, scoped bump allocator over one block of device memory.
///
/// Two ends: what lives as long as a scope of
/// the run, the timeframe, the iteration or the slice - parameters, the
/// timeframe's columns, a stage's output - grows up from the base; what must
/// outlive a base-end scope closed before it - the timeframe's track chunks,
/// the iteration's record chunks - or is needed only while a stage runs -
/// staging, scan storage - grows down from the top. Each end has its own
/// nested scopes; the two ends meet and never cross.

#ifndef GTRK_TRACKING_GPU_DEVICEARENA_H_
#define GTRK_TRACKING_GPU_DEVICEARENA_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace o2::gtrk
{

class DeviceRuntime;
class DeviceArena;
template <template <template <class> class> class S>
class DeviceColumns;

/// Which end of the arena an allocation or a scope belongs to.
enum class ArenaEnd : std::uint8_t {
  kBase = 0,    ///< grows up: storage of the scopes the run nests
  kScratch = 1, ///< grows down: chunks that outlive a base scope, and a stage's staging
};

/// Device memory an arena works in, owned by whoever hands it out (O2's
/// GPU reconstruction's pool, through GPUChainGenericTracking).
struct DeviceBlock {
  void* data = nullptr;
  std::size_t bytes = 0;
};

struct ArenaTag {
  const DeviceArena* arena = nullptr;
  ArenaEnd end = ArenaEnd::kBase;
  std::size_t depth = 0;
  std::uint64_t serial = 0;
};

struct DeviceAllocation {
  void* data = nullptr;
  std::size_t bytes = 0;
  ArenaTag tag;

  [[nodiscard]] explicit operator bool() const noexcept { return data != nullptr || bytes == 0U; }
};

class ArenaScope
{
 public:
  ArenaScope() = default;
  ~ArenaScope();
  ArenaScope(const ArenaScope&) = delete;
  ArenaScope& operator=(const ArenaScope&) = delete;
  ArenaScope(ArenaScope&& other) noexcept;
  ArenaScope& operator=(ArenaScope&& other) noexcept;

 private:
  friend class DeviceArena;
  template <template <template <class> class> class S>
  friend class DeviceColumns;
  ArenaScope(DeviceArena& arena, ArenaEnd end, std::size_t depth, std::uint64_t serial) noexcept : mArena(&arena), mEnd(end), mDepth(depth), mSerial(serial) {}
  /// Keep a transactional allocation and let its descriptor take the parent
  /// scope's tag. Ordinary arena scopes always rewind on destruction.
  void commit();
  void close() noexcept;

  DeviceArena* mArena = nullptr;
  ArenaEnd mEnd = ArenaEnd::kBase;
  std::size_t mDepth = 0;
  std::uint64_t mSerial = 0;
};

class DeviceArena
{
 public:
  DeviceArena(DeviceRuntime& runtime, DeviceBlock block);
  ~DeviceArena() = default;
  DeviceArena(const DeviceArena&) = delete;
  DeviceArena& operator=(const DeviceArena&) = delete;
  DeviceArena(DeviceArena&&) = delete;
  DeviceArena& operator=(DeviceArena&&) = delete;

  [[nodiscard]] DeviceAllocation allocate(std::size_t bytes, std::size_t alignment = 256U, ArenaEnd end = ArenaEnd::kBase);
  [[nodiscard]] ArenaScope scope(ArenaEnd end = ArenaEnd::kBase);
  /// Rewind both ends and close every scope; the block stays.
  void reset() noexcept;
  void resetPeak() noexcept { mPeak = used(); }
  /// The peak becomes at least \a peak: closes a window `resetPeak` opened.
  void raisePeak(std::size_t peak) noexcept { mPeak = std::max(mPeak, peak); }

  /// Both ends together, and each end on its own.
  [[nodiscard]] std::size_t used() const noexcept { return mUsed[0] + mUsed[1]; }
  [[nodiscard]] std::size_t used(ArenaEnd end) const noexcept { return mUsed[static_cast<std::size_t>(end)]; }
  [[nodiscard]] std::size_t peak() const noexcept { return mPeak; }
  [[nodiscard]] std::size_t capacity() const noexcept { return mCapacity; }
  [[nodiscard]] DeviceRuntime& runtime() const noexcept { return *mRuntime; }
  [[nodiscard]] ArenaTag currentTag(ArenaEnd end = ArenaEnd::kBase) const noexcept;
  void requireLive(const ArenaTag& tag) const;

 private:
  friend class ArenaScope;
  [[nodiscard]] bool live(const ArenaTag& tag) const noexcept;
  struct ScopeState {
    std::size_t offset = 0;
    std::uint64_t serial = 0;
  };

  [[nodiscard]] DeviceAllocation refuse(std::size_t bytes, ArenaEnd end);
  void closeScope(ArenaEnd end, std::size_t depth, std::uint64_t serial, bool commit) noexcept;

  DeviceRuntime* mRuntime = nullptr;
  void* mBase = nullptr;
  std::size_t mCapacity = 0;
  /// Bytes in use at each end: from the base up, and from the top down.
  std::array<std::size_t, 2> mUsed{};
  std::size_t mPeak = 0;
  std::uint64_t mNextSerial = 1;
  std::array<std::vector<ScopeState>, 2> mScopes;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_DEVICEARENA_H_
