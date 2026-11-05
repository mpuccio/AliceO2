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

#include "GenericTracking/Tracking/GPU/DeviceArena.h"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <utility>

namespace o2::gtrk
{

DeviceArena::DeviceArena(DeviceRuntime& runtime, DeviceBlock block) : mRuntime(&runtime), mBase(block.data), mCapacity(block.bytes)
{
  for (auto& scopes : mScopes) {
    scopes.push_back({.offset = 0U, .serial = mNextSerial++});
  }
}

DeviceAllocation DeviceArena::refuse(std::size_t bytes, ArenaEnd end)
{
  return {.data = nullptr, .bytes = bytes, .tag = currentTag(end)};
}

DeviceAllocation DeviceArena::allocate(std::size_t bytes, std::size_t alignment, ArenaEnd end)
{
  if (!std::has_single_bit(alignment)) {
    throw std::invalid_argument("device arena alignment must be a non-zero power of two");
  }
  auto* const origin = static_cast<std::byte*>(mBase);
  const std::size_t free = mCapacity - mUsed[0] - mUsed[1];
  if (bytes == 0U) {
    std::byte* at = nullptr;
    if (origin != nullptr) {
      at = end == ArenaEnd::kBase ? origin + mUsed[0] : origin + (mCapacity - mUsed[1]);
    }
    return {.data = at, .bytes = 0U, .tag = currentTag(end)};
  }
  if (origin == nullptr || bytes > free) {
    return refuse(bytes, end);
  }

  // Addresses, not offsets, are aligned.
  const std::size_t mask = alignment - 1U;
  if (end == ArenaEnd::kBase) {
    std::byte* const low = origin + mUsed[0];
    const std::size_t padding = (alignment - (std::bit_cast<std::uintptr_t>(low) & mask)) & mask;
    if (padding > free - bytes) {
      return refuse(bytes, end);
    }
    mUsed[0] += padding + bytes;
    mPeak = std::max(mPeak, used());
    return {.data = low + padding, .bytes = bytes, .tag = currentTag(end)};
  }
  std::byte* const unaligned = origin + (mCapacity - mUsed[1] - bytes);
  const std::size_t padding = std::bit_cast<std::uintptr_t>(unaligned) & mask;
  if (padding > free - bytes) {
    return refuse(bytes, end);
  }
  mUsed[1] += bytes + padding;
  mPeak = std::max(mPeak, used());
  return {.data = unaligned - padding, .bytes = bytes, .tag = currentTag(end)};
}

ArenaScope DeviceArena::scope(ArenaEnd end)
{
  auto& scopes = mScopes[static_cast<std::size_t>(end)];
  scopes.push_back({.offset = mUsed[static_cast<std::size_t>(end)], .serial = mNextSerial++});
  return ArenaScope{*this, end, scopes.size() - 1U, scopes.back().serial};
}

void DeviceArena::reset() noexcept
{
  mUsed = {};
  for (auto& scopes : mScopes) {
    scopes.clear();
    scopes.push_back({.offset = 0U, .serial = mNextSerial++});
  }
}

ArenaTag DeviceArena::currentTag(ArenaEnd end) const noexcept
{
  const auto& scopes = mScopes[static_cast<std::size_t>(end)];
  return {.arena = this, .end = end, .depth = scopes.size() - 1U, .serial = scopes.back().serial};
}

bool DeviceArena::live(const ArenaTag& tag) const noexcept
{
  if (tag.arena != this) {
    return false;
  }
  const auto& scopes = mScopes[static_cast<std::size_t>(tag.end)];
  return tag.depth < scopes.size() && scopes[tag.depth].serial == tag.serial;
}

void DeviceArena::requireLive(const ArenaTag& tag) const
{
  if (!live(tag)) {
    throw std::logic_error("device storage belongs to a rewound arena scope");
  }
}

void DeviceArena::closeScope(ArenaEnd end, std::size_t depth, std::uint64_t serial, bool commit) noexcept
{
  auto& scopes = mScopes[static_cast<std::size_t>(end)];
  // reset() deliberately makes outstanding ArenaScope objects inert.
  if (depth >= scopes.size() || scopes[depth].serial != serial) {
    return;
  }
  if (depth + 1U != scopes.size() || depth == 0U) {
    std::terminate();
  }
  if (!commit) {
    mUsed[static_cast<std::size_t>(end)] = scopes.back().offset;
  }
  scopes.pop_back();
}

ArenaScope::~ArenaScope() { close(); }

ArenaScope::ArenaScope(ArenaScope&& other) noexcept
  : mArena(std::exchange(other.mArena, nullptr)), mEnd(other.mEnd), mDepth(other.mDepth), mSerial(other.mSerial)
{
}

ArenaScope& ArenaScope::operator=(ArenaScope&& other) noexcept
{
  if (this != &other) {
    close();
    mArena = std::exchange(other.mArena, nullptr);
    mEnd = other.mEnd;
    mDepth = other.mDepth;
    mSerial = other.mSerial;
  }
  return *this;
}

void ArenaScope::commit()
{
  if (mArena == nullptr) {
    throw std::logic_error("cannot commit an empty device arena scope");
  }
  mArena->closeScope(mEnd, mDepth, mSerial, true);
  mArena = nullptr;
}

void ArenaScope::close() noexcept
{
  if (mArena != nullptr) {
    mArena->closeScope(mEnd, mDepth, mSerial, false);
    mArena = nullptr;
  }
}

} // namespace o2::gtrk
