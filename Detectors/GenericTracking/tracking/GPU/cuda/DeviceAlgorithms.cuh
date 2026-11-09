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

/// \file DeviceAlgorithms.cuh
/// \brief Device algorithms the tracking's passes share.
///
///   forEach            a device functor for every element, one thread each
///   forEachGroup       a device functor for every item, a group of lanes each
///   Gather             one array gathered through an index array
///   keyBits            the bits a radix sort needs for keys up to a largest
///   selectIndices      the indices a predicate holds for, in order
///   exclusiveScan      an exclusive sum with its storage on the scratch end
///   inclusiveScan      an inclusive scan by any operator, the same
///   sortPairs          a stable radix sort of key-value pairs, the same
///   sortKeys           a stable radix sort of keys alone, the same
///   segmentedStableSortPairs  key-value pairs sorted within each segment
///   sweepToFixedPoint  a pass repeated until it changes nothing, its flags
///                      read back once per pass
///
/// Each takes the runtime, arena and stream; temporaries go to the scratch
/// end and are released on return. Released scratch is only reused by work
/// queued after it on the same stream, as in GpuStage. Included by the .cu
/// sources of a backend library or its test library, in the backend's
/// namespace.

#ifndef GTRK_TRACKING_GPU_DEVICEALGORITHMS_CUH_
#define GTRK_TRACKING_GPU_DEVICEALGORITHMS_CUH_

#include "LaneGroup.cuh"
#include "Platform.cuh"

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/StageOutput.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::algorithms
{

namespace detail
{

/// \a body(element) for every element in [0, \a count): the launch covers
/// them all, one thread each (launchFor).
template <typename Body>
GPUg() void GTRK_KERNEL_BOUNDS forEachKernel(Body body, std::uint32_t count)
{
  const std::uint32_t element = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (element < count) {
    body(element);
  }
}

/// \a body.together(item, group) for every item in [0, \a count), each by a
/// group of `L` consecutive threads (`WarpGroup`).
/// No grid-stride loop: the launch covers every item, and a group past the
/// last does nothing as a whole.
template <std::uint32_t L, typename Body>
GPUg() void GTRK_KERNEL_BOUNDS forEachGroupKernel(Body body, std::uint32_t count)
{
  const std::uint64_t item = ((std::uint64_t{blockIdx.x} * blockDim.x) + threadIdx.x) / L;
  if (item < count) {
    body.together(static_cast<std::uint32_t>(item), WarpGroup<L>{});
  }
}

template <typename Predicate>
GPUg() void GTRK_KERNEL_BOUNDS flagIndicesKernel(Predicate predicate, std::uint32_t count, std::uint32_t* indices, std::uint8_t* flags)
{
  const std::uint32_t index = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (index < count) {
    indices[index] = index;
    flags[index] = predicate(index) ? 1U : 0U;
  }
}

/// CUB's and hipCUB's item counts are ints.
inline void requireItems(std::size_t items)
{
  if (items > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("a device algorithm's item count exceeds INT_MAX");
  }
}

} // namespace detail

/// \a body(element) -- a device functor -- for every element in [0, \a count),
/// one thread each; queued on \a stream.
template <typename Body>
void forEach(DeviceRuntime& runtime, DeviceStream& stream, std::uint32_t count, const Body& body)
{
  if (count == 0U) {
    return;
  }
  const auto shape = launchFor(count);
  detail::forEachKernel<<<shape.blocks, shape.threads, 0U, nativeStream(runtime, stream)>>>(body, count);
  checkLaunch();
}

/// \a body.together(item, group) -- a device functor -- for every item in
/// [0, \a count), `L` lanes each walking it in step (`WarpGroup`); queued on
/// \a stream.
template <std::uint32_t L, typename Body>
void forEachGroup(DeviceRuntime& runtime, DeviceStream& stream, std::uint32_t count, const Body& body)
{
  if (count == 0U) {
    return;
  }
  if (std::uint64_t{count} * L > std::numeric_limits<std::uint32_t>::max()) {
    throw std::length_error("forEachGroup: more threads than a launch counts");
  }
  const auto shape = launchFor(count * L);
  detail::forEachGroupKernel<L><<<shape.blocks, shape.threads, 0U, nativeStream(runtime, stream)>>>(body, count);
  checkLaunch();
}

/// to[i] = from[indices[i]]: one array gathered through an index array.
template <typename T>
struct Gather {
  const T* from = nullptr;
  const std::uint32_t* indices = nullptr;
  T* to = nullptr;
  GPUdi() void operator()(std::uint32_t i) const { to[i] = from[indices[i]]; }
};

/// The bits a radix sort needs for keys up to and including \a largest.
[[nodiscard]] inline int keyBits(std::uint32_t largest) { return std::max(1, static_cast<int>(std::bit_width(largest))); }

namespace detail
{
/// A library call with its temporary storage: \a call(temporary, bytes) once
/// with no storage, which sets the bytes it needs, then with that much of
/// the arena's scratch. kOutOfMemory if the scratch has not got it.
template <typename Call>
[[nodiscard]] StageStatus withTemporary(DeviceArena& arena, const Call& call)
{
  std::size_t temporaryBytes = 0U;
  GPUChkErrS(call(nullptr, temporaryBytes));
  DeviceSpan<std::byte> temporary;
  if (!temporary.allocate(arena, temporaryBytes, ArenaEnd::kScratch)) {
    return StageStatus::kOutOfMemory;
  }
  GPUChkErrS(call(temporary.data(), temporaryBytes));
  return StageStatus::kComplete;
}
} // namespace detail

/// The indices in [0, \a count) for which \a predicate -- a device functor
/// `bool(std::uint32_t)` -- holds, in increasing order, into \a output
/// (\a count entries on \a outputEnd); \a selected is how many, read back.
template <typename Predicate>
[[nodiscard]] StageStatus selectIndices(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, std::uint32_t count, const Predicate& predicate, DeviceSpan<std::uint32_t>& output, std::uint32_t& selected, ArenaEnd outputEnd = ArenaEnd::kBase)
{
  selected = 0U;
  detail::requireItems(count);
  if (!output.allocate(arena, count, outputEnd)) {
    return StageStatus::kOutOfMemory;
  }
  if (count == 0U) {
    return StageStatus::kComplete;
  }
  auto scratch = arena.scope(ArenaEnd::kScratch);
  DeviceSpan<std::uint32_t> indices;
  DeviceSpan<std::uint8_t> flags;
  DeviceSpan<std::uint32_t> selectedCount;
  if (!indices.allocate(arena, count, ArenaEnd::kScratch) || !flags.allocate(arena, count, ArenaEnd::kScratch) || !selectedCount.allocate(arena, 1U, ArenaEnd::kScratch)) {
    return StageStatus::kOutOfMemory;
  }
  const auto native = nativeStream(runtime, stream);
  const auto shape = launchFor(count);
  detail::flagIndicesKernel<<<shape.blocks, shape.threads, 0U, native>>>(predicate, count, indices.data(), flags.data());
  checkLaunch();
  const auto flagged = detail::withTemporary(arena, [&](void* temporary, std::size_t& bytes) { return cubapi::DeviceSelect::Flagged(temporary, bytes, indices.data(), flags.data(), output.data(), selectedCount.data(), static_cast<int>(count), native); });
  if (flagged != StageStatus::kComplete) {
    return flagged;
  }
  download(selectedCount, &selected, 1U, stream);
  stream.synchronize();
  return StageStatus::kComplete;
}

/// `output[i]` = the sum of `input[0..i)`, over \a items entries of any
/// integer type; queued on \a stream, nothing read back.
template <typename T>
[[nodiscard]] StageStatus exclusiveScan(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, const T* input, T* output, std::size_t items)
{
  detail::requireItems(items);
  if (items == 0U) {
    return StageStatus::kComplete;
  }
  auto scratch = arena.scope(ArenaEnd::kScratch);
  const auto native = nativeStream(runtime, stream);
  return detail::withTemporary(arena, [&](void* temporary, std::size_t& bytes) { return cubapi::DeviceScan::ExclusiveSum(temporary, bytes, input, output, static_cast<int>(items), native); });
}

/// `output[i]` = \a op over `input[0..i]`, over \a items entries; queued on
/// \a stream, nothing read back.
template <typename T, typename Op>
[[nodiscard]] StageStatus inclusiveScan(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, const T* input, T* output, std::size_t items, const Op& op)
{
  detail::requireItems(items);
  if (items == 0U) {
    return StageStatus::kComplete;
  }
  auto scratch = arena.scope(ArenaEnd::kScratch);
  const auto native = nativeStream(runtime, stream);
  return detail::withTemporary(arena, [&](void* temporary, std::size_t& bytes) { return cubapi::DeviceScan::InclusiveScan(temporary, bytes, input, output, op, static_cast<int>(items), native); });
}

/// Sort \a count pairs by \a keys' low \a endBit bits, stably: \a keys into
/// \a sortedKeys, \a values into \a sortedValues. Queued on \a stream.
template <typename Key, typename Value>
[[nodiscard]] StageStatus sortPairs(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, const Key* keys, Key* sortedKeys, const Value* values, Value* sortedValues, std::uint32_t count, int endBit)
{
  detail::requireItems(count);
  auto scratch = arena.scope(ArenaEnd::kScratch);
  const auto native = nativeStream(runtime, stream);
  return detail::withTemporary(arena, [&](void* temporary, std::size_t& bytes) { return cubapi::DeviceRadixSort::SortPairs(temporary, bytes, keys, sortedKeys, values, sortedValues, static_cast<int>(count), 0, endBit, native); });
}

/// Sort \a count keys by their low \a endBit bits, stably, \a keys into
/// \a sortedKeys. Queued on \a stream.
template <typename Key>
[[nodiscard]] StageStatus sortKeys(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, const Key* keys, Key* sortedKeys, std::uint32_t count, int endBit)
{
  detail::requireItems(count);
  auto scratch = arena.scope(ArenaEnd::kScratch);
  const auto native = nativeStream(runtime, stream);
  return detail::withTemporary(arena, [&](void* temporary, std::size_t& bytes) { return cubapi::DeviceRadixSort::SortKeys(temporary, bytes, keys, sortedKeys, static_cast<int>(count), 0, endBit, native); });
}

/// Sort the items of every segment [\a begin[s], \a end[s]) by key, stably,
/// keys and values together, through the double buffers: afterwards each
/// buffer's `Current()` holds the sorted items. The library chooses how by
/// the segments' lengths. Queued on \a stream.
template <typename Key, typename Value, typename Offsets>
[[nodiscard]] StageStatus segmentedStableSortPairs(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, cubapi::DoubleBuffer<Key>& keys, cubapi::DoubleBuffer<Value>& values, std::uint32_t items, std::uint32_t segments, Offsets begin, Offsets end)
{
  detail::requireItems(items);
  detail::requireItems(segments);
  if (items == 0U) {
    return StageStatus::kComplete;
  }
  auto scratch = arena.scope(ArenaEnd::kScratch);
  const auto native = nativeStream(runtime, stream);
  return detail::withTemporary(arena, [&](void* temporary, std::size_t& bytes) { return cubapi::DeviceSegmentedSort::StableSortPairs(temporary, bytes, keys, values, static_cast<int>(items), static_cast<int>(segments), begin, end, native); });
}

/// What a sweep to a fixed point ended with.
template <std::size_t Flags>
struct SweepResult {
  StageStatus status = StageStatus::kComplete; ///< kOutOfMemory if the flags could not be had
  std::array<std::uint32_t, Flags> flags{};    ///< [0] is the change flag, 0 at the end; the others are sticky
  std::uint32_t changingSweeps = 0U;           ///< passes that raised the change flag; one readback per pass, the last unchanged
};

/// Repeat \a sweep until a pass leaves the change flag down. \a sweep is a
/// host callable `void(std::uint32_t* flags)` that queues one pass on
/// \a stream, the pass raising `flags[0]` when it changed anything and the
/// others as it needs. `flags[0]` is cleared before every pass, the others
/// only before the first, so they tell whether any pass raised them. All of
/// them come back in one copy per pass. More than \a maxChangingSweeps
/// changing passes is a broken invariant of the caller's -- its fixed point
/// is known to be that close -- and throws std::logic_error.
template <std::size_t Flags, typename Sweep>
[[nodiscard]] SweepResult<Flags> sweepToFixedPoint(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, std::uint32_t maxChangingSweeps, const Sweep& sweep)
{
  static_assert(Flags >= 1U, "a sweep needs its change flag");
  SweepResult<Flags> result;
  auto scratch = arena.scope(ArenaEnd::kScratch);
  DeviceSpan<std::uint32_t> flags;
  if (!flags.allocate(arena, Flags, ArenaEnd::kScratch)) {
    result.status = StageStatus::kOutOfMemory;
    return result;
  }
  runtime.fill(flags.data(), 0U, Flags * sizeof(std::uint32_t), stream);
  while (true) {
    runtime.fill(flags.data(), 0U, sizeof(std::uint32_t), stream);
    sweep(flags.data());
    download(flags, result.flags.data(), Flags, stream);
    stream.synchronize();
    if (result.flags[0] == 0U) {
      return result;
    }
    if (++result.changingSweeps > maxChangingSweeps) {
      throw std::logic_error("a sweep did not reach its fixed point within its bound: a broken invariant");
    }
  }
}

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::algorithms

#endif // GTRK_TRACKING_GPU_DEVICEALGORITHMS_CUH_
