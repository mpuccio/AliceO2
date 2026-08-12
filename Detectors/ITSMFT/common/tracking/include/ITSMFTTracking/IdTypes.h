// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_IDTYPES_H_
#define ALICEO2_ITSMFT_TRACKING_IDTYPES_H_

#include <cstdint>
#include <limits>
#include <type_traits>

#include "GPUCommonDef.h"

namespace o2::itsmft::tracking
{

// The coordinate convention of a surface and every state defined on it.
enum class SurfaceKind : uint8_t {
  Undefined,
  Cylinder,
  Disk
};

static_assert(std::is_same_v<std::underlying_type_t<SurfaceKind>, uint8_t>);
static_assert(sizeof(SurfaceKind) == sizeof(uint8_t));

namespace detail
{
template <typename Tag, typename ValueType>
class Identifier
{
 public:
  static constexpr ValueType InvalidValue = std::numeric_limits<ValueType>::max();

  GPUhdDefault() constexpr Identifier() noexcept = default;
  GPUhdDefault() explicit constexpr Identifier(ValueType value) noexcept : mValue{value} {}

  GPUhdi() constexpr ValueType value() const noexcept { return mValue; }
  GPUhdi() constexpr bool isValid() const noexcept { return mValue != InvalidValue; }
  GPUhdi() static constexpr Identifier invalid() noexcept { return Identifier{InvalidValue}; }

  GPUhdi() friend constexpr bool operator==(Identifier lhs, Identifier rhs) noexcept { return lhs.mValue == rhs.mValue; }
  GPUhdi() friend constexpr bool operator!=(Identifier lhs, Identifier rhs) noexcept { return !(lhs == rhs); }
  GPUhdi() friend constexpr bool operator<(Identifier lhs, Identifier rhs) noexcept { return lhs.mValue < rhs.mValue; }

 private:
  ValueType mValue{InvalidValue};
};
} // namespace detail

struct SurfaceIdTag;
struct TransitionIdTag;
struct CellTopologyIdTag;
struct ClusterSourceIdTag;
struct ClusterIndexTag;

using SurfaceId = detail::Identifier<SurfaceIdTag, uint16_t>;
using TransitionId = detail::Identifier<TransitionIdTag, uint16_t>;
using CellTopologyId = detail::Identifier<CellTopologyIdTag, uint16_t>;
using ClusterSourceId = detail::Identifier<ClusterSourceIdTag, uint16_t>;
using SurfaceMeasurementIndex = detail::Identifier<ClusterIndexTag, uint32_t>;

GPUhdi() constexpr bool isRecognizedSurfaceKind(SurfaceKind kind) noexcept
{
  return kind == SurfaceKind::Cylinder || kind == SurfaceKind::Disk;
}

inline constexpr uint32_t MaxLayoutSurfaces = 32;
inline constexpr uint32_t MaxLayoutTransitions = MaxLayoutSurfaces * (MaxLayoutSurfaces - 1);
inline constexpr uint32_t MaxLayoutCellTopologies = MaxLayoutSurfaces * (MaxLayoutSurfaces - 1) * (MaxLayoutSurfaces - 1);

static_assert(MaxLayoutTransitions < TransitionId::InvalidValue);
static_assert(MaxLayoutCellTopologies < CellTopologyId::InvalidValue);

} // namespace o2::itsmft::tracking

#endif
