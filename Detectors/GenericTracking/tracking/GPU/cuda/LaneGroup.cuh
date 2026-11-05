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

/// \file LaneGroup.cuh
/// \brief Lanes of one warp that walk in step: the
///        devices' group of Lanes.h's contract.

#ifndef GTRK_TRACKING_GPU_LANEGROUP_CUH_
#define GTRK_TRACKING_GPU_LANEGROUP_CUH_

#include "Platform.cuh"

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
{

/// `L` consecutive threads of a block, aligned to `L` inside one warp: lane
/// `i` is the thread `i` places after the group's first. Its collectives
/// (`SoloGroup`'s contract, Lanes.h) are the warp's cross-lane operations on
/// the group's lanes alone, so the other groups of the warp may be anywhere
/// else meanwhile: CUDA's `_sync` forms with the group's mask and width `L`,
/// HIP's with width `L` and its ballot masked to the group. Every lane of the
/// group must reach each call.
template <std::uint32_t L>
struct WarpGroup {
  static_assert(L >= 2U && (L & (L - 1U)) == 0U && static_cast<std::uint32_t>(kWarpSize) % L == 0U, "a group is a power-of-two share of one warp");
  static constexpr std::uint32_t kLanes = L;

  [[nodiscard]] GPUdii() std::uint32_t lane() const { return threadIdx.x % L; }

  /// Bit `i` for lane `i` of the group where \a set holds.
  [[nodiscard]] GPUdii() std::uint64_t ballot(bool set) const
  {
#ifdef __HIPCC__
    const std::uint64_t all = __ballot(set);
    if constexpr (L == 64U) {
      return all;
    } else {
      return (all >> first()) & ((std::uint64_t{1} << L) - 1U);
    }
#else
    return (__ballot_sync(mask(), set) & mask()) >> first();
#endif
  }

  /// Lane \a from's \a value on every lane of the group, word by word.
  template <typename T>
  [[nodiscard]] GPUdii() T broadcast(const T& value, std::uint32_t from) const
  {
    static_assert(std::is_trivially_copyable_v<T> && sizeof(T) % sizeof(std::uint32_t) == 0U, "a broadcast value is whole 32-bit words");
    constexpr std::size_t kWords = sizeof(T) / sizeof(std::uint32_t);
    std::uint32_t words[kWords];
    __builtin_memcpy(words, &value, sizeof(T));
    for (std::size_t i = 0; i < kWords; ++i) {
#ifdef __HIPCC__
      words[i] = __shfl(words[i], static_cast<int>(from), static_cast<int>(L));
#else
      words[i] = __shfl_sync(mask(), words[i], static_cast<int>(from), static_cast<int>(L));
#endif
    }
    T out;
    __builtin_memcpy(&out, words, sizeof(T));
    return out;
  }

 private:
  /// The group's first lane in its warp.
  [[nodiscard]] GPUdii() std::uint32_t first() const { return (threadIdx.x % static_cast<std::uint32_t>(kWarpSize)) - lane(); }

#ifndef __HIPCC__
  /// The group's lanes in the warp's mask.
  [[nodiscard]] GPUdii() std::uint32_t mask() const
  {
    if constexpr (L == 32U) {
      return 0xFFFFFFFFU;
    } else {
      return ((1U << L) - 1U) << first();
    }
  }
#endif
};

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE

#endif // GTRK_TRACKING_GPU_LANEGROUP_CUH_
