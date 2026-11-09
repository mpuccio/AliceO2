// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "GenericTracking/Base/BoundedAllocator.h"

#include <format>

#include "GPUCommonLogger.h"

namespace o2::gtrk
{
namespace
{
constexpr double kGigabyte = 1024. * 1024. * 1024.;
}

BoundedMemoryResource::MemoryLimitExceeded::MemoryLimitExceeded(size_t attempted, size_t used, size_t max)
{
  if (attempted != 0) {
    mMsg = std::format("Reached set memory limit (attempted: {}, used: {}, max: {})", attempted, used, max);
  } else {
    mMsg = std::format("New set maximum below current used (newMax: {}, used: {})", max, used);
  }
}

const char* BoundedMemoryResource::MemoryLimitExceeded::what() const noexcept { return mMsg.c_str(); }

std::pmr::memory_resource* BoundedMemoryResource::cachingUpstream()
{
  static std::pmr::synchronized_pool_resource pool{std::pmr::get_default_resource()};
  return &pool;
}

BoundedMemoryResource::BoundedMemoryResource(size_t maxBytes, std::pmr::memory_resource* upstream)
  : mMaxMemory(maxBytes),
    mUpstream(upstream != nullptr ? upstream : cachingUpstream())
{
}

BoundedMemoryResource::BoundedMemoryResource(std::unique_ptr<std::pmr::memory_resource> upstream, size_t maxBytes)
  : mMaxMemory(maxBytes),
    mOwnedUpstream(std::move(upstream)),
    mUpstream(mOwnedUpstream.get())
{
}

void* BoundedMemoryResource::do_allocate(size_t bytes, size_t alignment)
{
  size_t newUsed{0};
  size_t currentUsed{mUsedMemory.load(std::memory_order_relaxed)};
  for (;;) {
    const size_t limit = mMaxMemory.load(std::memory_order_relaxed);
    // Test before adding, so even SIZE_MAX requests cannot wrap under budget.
    if (currentUsed > limit || bytes > limit - currentUsed) {
      const size_t attempted = bytes > std::numeric_limits<size_t>::max() - currentUsed ? std::numeric_limits<size_t>::max() : currentUsed + bytes;
      mCountThrow.fetch_add(1, std::memory_order_relaxed);
      throw MemoryLimitExceeded(attempted, currentUsed, limit);
    }
    newUsed = currentUsed + bytes;
    if (mUsedMemory.compare_exchange_weak(currentUsed, newUsed, std::memory_order_acq_rel, std::memory_order_relaxed)) {
      break;
    }
  }

  void* p{nullptr};
  try {
    p = mUpstream->allocate(bytes, alignment);
  } catch (...) {
    mUsedMemory.fetch_sub(bytes, std::memory_order_relaxed);
    throw;
  }

  mAllocations.fetch_add(1, std::memory_order_relaxed);
  mAllocatedBytes.fetch_add(bytes, std::memory_order_relaxed);

  size_t peak = mPeakUsedMemory.load(std::memory_order_relaxed);
  while (newUsed > peak && !mPeakUsedMemory.compare_exchange_weak(peak, newUsed, std::memory_order_relaxed)) {
  }

  return p;
}

void BoundedMemoryResource::do_deallocate(void* p, size_t bytes, size_t alignment)
{
  mUpstream->deallocate(p, bytes, alignment);
  mUsedMemory.fetch_sub(bytes, std::memory_order_relaxed);
}

bool BoundedMemoryResource::do_is_equal(const std::pmr::memory_resource& other) const noexcept { return this == &other; }

size_t BoundedMemoryResource::getUsedMemory() const noexcept { return mUsedMemory.load(std::memory_order_relaxed); }

size_t BoundedMemoryResource::getMaxMemory() const noexcept { return mMaxMemory.load(std::memory_order_relaxed); }

size_t BoundedMemoryResource::getThrowCount() const noexcept { return mCountThrow.load(std::memory_order_relaxed); }

size_t BoundedMemoryResource::getPeakMemory() const noexcept { return mPeakUsedMemory.load(std::memory_order_relaxed); }

size_t BoundedMemoryResource::getPeakMemoryDelta() const noexcept
{
  const size_t peak = mPeakUsedMemory.load(std::memory_order_relaxed);
  const size_t baseline = mPeakBaselineMemory.load(std::memory_order_relaxed);
  return peak > baseline ? peak - baseline : 0;
}

size_t BoundedMemoryResource::resetPeakMemory() noexcept
{
  const size_t used = mUsedMemory.load(std::memory_order_acquire);
  mPeakBaselineMemory.store(used, std::memory_order_release);
  return mPeakUsedMemory.exchange(used, std::memory_order_acq_rel);
}

void BoundedMemoryResource::raisePeakMemory(size_t peak) noexcept
{
  size_t current = mPeakUsedMemory.load(std::memory_order_relaxed);
  while (peak > current && !mPeakUsedMemory.compare_exchange_weak(current, peak, std::memory_order_relaxed)) {
  }
}

void BoundedMemoryResource::setMaxMemory(size_t max)
{
  size_t current = mMaxMemory.load(std::memory_order_relaxed);
  if (max == current) {
    return;
  }
  for (;;) {
    const size_t used = mUsedMemory.load(std::memory_order_acquire);
    if (used > max) {
      mCountThrow.fetch_add(1, std::memory_order_relaxed);
      throw MemoryLimitExceeded(0, used, max);
    }
    if (mMaxMemory.compare_exchange_weak(current, max, std::memory_order_release, std::memory_order_relaxed)) {
      return;
    }
    if (current == max) {
      return;
    }
  }
}

std::string BoundedMemoryResource::asString() const
{
  const auto throwCount = mCountThrow.load(std::memory_order_relaxed);
  const auto used = static_cast<double>(mUsedMemory.load(std::memory_order_relaxed));
  const auto peak = static_cast<double>(mPeakUsedMemory.load(std::memory_order_relaxed));
  const auto peakDelta = static_cast<double>(getPeakMemoryDelta());
  const auto maxMemory = mMaxMemory.load(std::memory_order_relaxed);
  std::string result;
  if (maxMemory == std::numeric_limits<size_t>::max()) {
    result += std::format("maxthrow={} maxmem=unbounded used={:.2f} GB stagepeak={:.2f} GB stagealloc={:.2f} GB", throwCount, used / kGigabyte, peak / kGigabyte, peakDelta / kGigabyte);
  } else {
    result += std::format("maxthrow={} maxmem={:.2f} GB used={:.2f} GB ({:.2f}%) stagepeak={:.2f} GB stagealloc={:.2f} GB", throwCount, static_cast<double>(maxMemory) / kGigabyte, used / kGigabyte, 100.0 * used / static_cast<double>(maxMemory), peak / kGigabyte, peakDelta / kGigabyte);
  }
  return result;
}

void BoundedMemoryResource::print() const { LOGP(info, "{}", asString()); }

} // namespace o2::gtrk
