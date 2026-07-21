// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENTINDEX_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENTINDEX_H_

#include <cstdint>
#include <limits>
#include <type_traits>

#include "GPUCommonDef.h"

namespace o2::itsmft::tracking
{

// Position in one surface's complete TimeFrame measurement span.
class SurfaceMeasurementIndex
{
 public:
  using ValueType = uint32_t;
  static constexpr ValueType InvalidValue = std::numeric_limits<ValueType>::max();

  GPUhdDefault() constexpr SurfaceMeasurementIndex() noexcept = default;
  GPUhdDefault() explicit constexpr SurfaceMeasurementIndex(ValueType value) noexcept : mValue{value} {}

  GPUhdi() constexpr ValueType value() const noexcept { return mValue; }
  GPUhdi() constexpr bool isValid() const noexcept { return mValue != InvalidValue; }
  GPUhdi() static constexpr SurfaceMeasurementIndex invalid() noexcept { return SurfaceMeasurementIndex{InvalidValue}; }

  GPUhdi() friend constexpr bool operator==(SurfaceMeasurementIndex lhs, SurfaceMeasurementIndex rhs) noexcept
  {
    return lhs.mValue == rhs.mValue;
  }
  GPUhdi() friend constexpr bool operator!=(SurfaceMeasurementIndex lhs, SurfaceMeasurementIndex rhs) noexcept { return !(lhs == rhs); }

 private:
  ValueType mValue{InvalidValue};
};

static_assert(std::is_standard_layout_v<SurfaceMeasurementIndex>);
static_assert(std::is_trivially_copyable_v<SurfaceMeasurementIndex>);
static_assert(sizeof(SurfaceMeasurementIndex) == sizeof(uint32_t));
static_assert(alignof(SurfaceMeasurementIndex) == alignof(uint32_t));

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENTINDEX_H_ */
