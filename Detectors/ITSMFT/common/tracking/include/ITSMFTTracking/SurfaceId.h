// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEID_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEID_H_

#include <cstdint>
#include <limits>
#include <type_traits>

#include "GPUCommonDef.h"

namespace o2::itsmft::tracking
{

namespace detail
{
template <typename Tag>
class Identifier16
{
 public:
  using ValueType = uint16_t;
  static constexpr ValueType InvalidValue = std::numeric_limits<ValueType>::max();

  GPUhdDefault() constexpr Identifier16() noexcept = default;
  GPUhdDefault() explicit constexpr Identifier16(ValueType value) noexcept : mValue{value} {}

  GPUhdi() constexpr ValueType value() const noexcept { return mValue; }
  GPUhdi() constexpr bool isValid() const noexcept { return mValue != InvalidValue; }
  GPUhdi() static constexpr Identifier16 invalid() noexcept { return Identifier16{InvalidValue}; }

  GPUhdi() friend constexpr bool operator==(Identifier16 lhs, Identifier16 rhs) noexcept { return lhs.mValue == rhs.mValue; }
  GPUhdi() friend constexpr bool operator!=(Identifier16 lhs, Identifier16 rhs) noexcept { return !(lhs == rhs); }
  GPUhdi() friend constexpr bool operator<(Identifier16 lhs, Identifier16 rhs) noexcept { return lhs.mValue < rhs.mValue; }

 private:
  ValueType mValue{InvalidValue};
};
} // namespace detail

struct SurfaceIdTag;
struct TransitionIdTag;
struct CellTopologyIdTag;
struct ClusterSourceIdTag;

using SurfaceId = detail::Identifier16<SurfaceIdTag>;
using TransitionId = detail::Identifier16<TransitionIdTag>;
using CellTopologyId = detail::Identifier16<CellTopologyIdTag>;
// Dense, TimeFrame-local input-stream ID; all-ones is invalid.
using ClusterSourceId = detail::Identifier16<ClusterSourceIdTag>;

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

inline constexpr uint32_t MaxLayoutSurfaces = 32;
inline constexpr uint32_t MaxLayoutTransitions = MaxLayoutSurfaces * (MaxLayoutSurfaces - 1);
inline constexpr uint32_t MaxLayoutCellTopologies = MaxLayoutSurfaces * (MaxLayoutSurfaces - 1) * (MaxLayoutSurfaces - 1);

static_assert(MaxLayoutTransitions < TransitionId::InvalidValue);
static_assert(MaxLayoutCellTopologies < CellTopologyId::InvalidValue);
static_assert(std::is_standard_layout_v<SurfaceId> && std::is_trivially_copyable_v<SurfaceId>);
static_assert(std::is_standard_layout_v<TransitionId> && std::is_trivially_copyable_v<TransitionId>);
static_assert(std::is_standard_layout_v<CellTopologyId> && std::is_trivially_copyable_v<CellTopologyId>);
static_assert(std::is_standard_layout_v<ClusterSourceId> && std::is_trivially_copyable_v<ClusterSourceId>);
static_assert(std::is_standard_layout_v<SurfaceMeasurementIndex>);
static_assert(std::is_trivially_copyable_v<SurfaceMeasurementIndex>);
static_assert(sizeof(SurfaceId) == sizeof(uint16_t));
static_assert(sizeof(TransitionId) == sizeof(uint16_t));
static_assert(sizeof(CellTopologyId) == sizeof(uint16_t));
static_assert(sizeof(ClusterSourceId) == sizeof(uint16_t));
static_assert(sizeof(SurfaceMeasurementIndex) == sizeof(uint32_t));
static_assert(alignof(SurfaceMeasurementIndex) == alignof(uint32_t));

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACEID_H_ */
