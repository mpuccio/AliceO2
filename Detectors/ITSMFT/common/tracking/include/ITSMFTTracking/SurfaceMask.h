// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEMASK_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEMASK_H_

#include <cstdint>
#include <type_traits>

#ifndef GPUCA_GPUCODE
#include <fmt/format.h>
#include <string>

#include <gsl/span>
#endif

#include "GPUCommonDef.h"
#include "GPUCommonMath.h"
#include "ITStracking/Constants.h"
#include "ITSMFTTracking/SurfaceId.h"

namespace o2::itsmft::tracking
{

struct LayerMask {
  GPUhdDefault() constexpr LayerMask() noexcept = default;
  GPUhdDefault() constexpr LayerMask(uint16_t mask) noexcept : mBits{mask} {}
  GPUhdDefault() constexpr LayerMask(int layer0, int layer1, int layer2) noexcept
    : mBits{static_cast<uint16_t>((uint16_t(1) << layer0) | (uint16_t(1) << layer1) | (uint16_t(1) << layer2))}
  {
  }
  GPUhdi() constexpr operator uint16_t() const noexcept { return mBits; }
  GPUhdi() constexpr uint16_t value() const noexcept { return mBits; }
  GPUhdi() constexpr void set(int layer) noexcept { mBits |= (uint16_t(1) << layer); }

  GPUhdi() LayerMask operator~() const noexcept { return LayerMask{static_cast<uint16_t>(~mBits)}; }
  GPUhdi() LayerMask operator&(LayerMask other) const noexcept { return LayerMask{static_cast<uint16_t>(mBits & other.mBits)}; }
  GPUhdi() LayerMask operator|(LayerMask other) const noexcept { return LayerMask{static_cast<uint16_t>(mBits | other.mBits)}; }
  GPUhdi() LayerMask& operator&=(LayerMask other) noexcept
  {
    mBits &= other.mBits;
    return *this;
  }
  GPUhdi() LayerMask& operator|=(LayerMask other) noexcept
  {
    mBits |= other.mBits;
    return *this;
  }

  GPUhdi() bool empty() const noexcept { return mBits == 0; }
  GPUhdi() bool has(int layer) const noexcept { return mBits & (uint16_t(1) << layer); }
  GPUhdi() bool isSubsetOf(LayerMask allowed) const noexcept { return (*this & ~allowed).empty(); }
  GPUhdi() bool isAllowedHoleMask(int maxHoles, LayerMask allowedHoleMask) const noexcept
  {
    const int allowedHoles = maxHoles > 0 ? maxHoles : 0;
    return count() <= allowedHoles && isSubsetOf(allowedHoleMask);
  }
  GPUhdi() bool isAllowed(int maxHoles, LayerMask allowedHoleMask) const noexcept
  {
    return holeMask().isAllowedHoleMask(maxHoles, allowedHoleMask);
  }
  GPUhdi() int length() const noexcept { return empty() ? 0 : last() - first() + 1; }
  GPUhdi() int count() const noexcept { return static_cast<int>(o2::gpu::GPUCommonMath::Popcount(mBits)); }
  GPUhdi() int first() const noexcept { return mBits ? static_cast<int>(o2::gpu::GPUCommonMath::Ctz(mBits)) : o2::its::constants::UnusedIndex; }
  GPUhdi() int last() const noexcept { return mBits ? 31 - static_cast<int>(o2::gpu::GPUCommonMath::Clz(mBits)) : o2::its::constants::UnusedIndex; }
  GPUhdi() LayerMask holeMask() const noexcept
  {
    return empty() ? LayerMask{0} : (span(first(), last()) & ~(*this));
  }

  GPUhdi() int slot(int layer) const noexcept
  {
    if (!has(layer)) {
      return o2::its::constants::UnusedIndex;
    }
    const uint32_t lowerLayers = (uint32_t(1) << layer) - 1;
    return static_cast<int>(o2::gpu::GPUCommonMath::Popcount(static_cast<uint32_t>(mBits) & lowerLayers));
  }

  static GPUhdi() LayerMask span(int fromLayer, int toLayer) noexcept
  {
    if (fromLayer > toLayer) {
      return 0;
    }
    const uint32_t upper = (uint32_t(1) << (toLayer + 1)) - 1;
    const uint32_t lower = (uint32_t(1) << fromLayer) - 1;
    return static_cast<uint16_t>(upper & ~lower);
  }

  static GPUhdi() LayerMask skipped(int fromLayer, int toLayer) noexcept
  {
    return (toLayer - fromLayer <= 1) ? LayerMask{0} : span(fromLayer + 1, toLayer - 1);
  }

#ifndef GPUCA_GPUCODE
  std::string asString() const { return fmt::format("{:016b}", mBits); }
#endif

 private:
  uint16_t mBits{0};
};

static_assert(std::is_standard_layout_v<LayerMask>);
static_assert(std::is_trivially_copyable_v<LayerMask>);
static_assert(sizeof(LayerMask) == sizeof(uint16_t));
static_assert(alignof(LayerMask) == alignof(uint16_t));

struct SurfaceMask {
  using ValueType = uint32_t;

  GPUhdDefault() constexpr SurfaceMask() noexcept = default;
  GPUhdDefault() explicit constexpr SurfaceMask(ValueType bits) noexcept : mBits{bits} {}

  GPUhdi() constexpr ValueType value() const noexcept { return mBits; }
  GPUhdi() constexpr bool empty() const noexcept { return mBits == 0; }

  GPUhdi() constexpr bool has(SurfaceId surface) const noexcept
  {
    return surface.isValid() && surface.value() < MaxLayoutSurfaces && (mBits & (ValueType{1} << surface.value()));
  }
  GPUhdi() constexpr void set(SurfaceId surface) noexcept
  {
    if (surface.isValid() && surface.value() < MaxLayoutSurfaces) {
      mBits |= ValueType{1} << surface.value();
    }
  }
  GPUhdi() constexpr void reset(SurfaceId surface) noexcept
  {
    if (surface.isValid() && surface.value() < MaxLayoutSurfaces) {
      mBits &= ~(ValueType{1} << surface.value());
    }
  }

  GPUhdi() SurfaceMask operator~() const noexcept { return SurfaceMask{~mBits}; }
  GPUhdi() SurfaceMask operator&(SurfaceMask other) const noexcept { return SurfaceMask{mBits & other.mBits}; }
  GPUhdi() SurfaceMask operator|(SurfaceMask other) const noexcept { return SurfaceMask{mBits | other.mBits}; }
  GPUhdi() friend constexpr bool operator==(SurfaceMask lhs, SurfaceMask rhs) noexcept { return lhs.mBits == rhs.mBits; }
  GPUhdi() friend constexpr bool operator!=(SurfaceMask lhs, SurfaceMask rhs) noexcept { return !(lhs == rhs); }
  GPUhdi() SurfaceMask& operator&=(SurfaceMask other) noexcept
  {
    mBits &= other.mBits;
    return *this;
  }
  GPUhdi() SurfaceMask& operator|=(SurfaceMask other) noexcept
  {
    mBits |= other.mBits;
    return *this;
  }

  GPUhdi() bool isSubsetOf(SurfaceMask allowed) const noexcept { return (*this & ~allowed).empty(); }
  GPUhdi() int count() const noexcept { return static_cast<int>(o2::gpu::GPUCommonMath::Popcount(mBits)); }
  GPUhdi() int first() const noexcept { return mBits ? static_cast<int>(o2::gpu::GPUCommonMath::Ctz(mBits)) : -1; }
  GPUhdi() int last() const noexcept { return mBits ? 31 - static_cast<int>(o2::gpu::GPUCommonMath::Clz(mBits)) : -1; }

#ifndef GPUCA_GPUCODE
  std::string asString() const { return fmt::format("{:032b}", mBits); }
#endif

 private:
  ValueType mBits{0};
};

static_assert(std::is_standard_layout_v<SurfaceMask>);
static_assert(std::is_trivially_copyable_v<SurfaceMask>);
static_assert(sizeof(SurfaceMask) == sizeof(uint32_t));
static_assert(alignof(SurfaceMask) == alignof(uint32_t));

#ifndef GPUCA_GPUCODE
// Converts a positional LayerMask -- each set bit is a *position* in
// `orderedSurfaces`, never a numeric comparison against the SurfaceId values
// it contains -- into the corresponding global SurfaceMask. Only the first
// `activeCount` positions are considered; `activeCount` must not exceed
// `orderedSurfaces.size()` (unchecked, as with the two call sites this
// consolidates).
inline SurfaceMask positionalSurfaceMask(LayerMask layerMask, gsl::span<const SurfaceId> orderedSurfaces, uint32_t activeCount) noexcept
{
  SurfaceMask result;
  for (uint32_t position = 0; position < activeCount; ++position) {
    if (layerMask.has(position)) {
      result.set(orderedSurfaces[position]);
    }
  }
  return result;
}
#endif

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACEMASK_H_ */
