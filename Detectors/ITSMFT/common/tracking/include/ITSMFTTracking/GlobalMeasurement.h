// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_GLOBALMEASUREMENT_H_
#define ALICEO2_ITSMFT_TRACKING_GLOBALMEASUREMENT_H_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/SurfaceId.h"

namespace o2::itsmft::tracking
{

struct ClusterRef {
  GPUhdDefault() constexpr ClusterRef() noexcept = default;
  GPUhdDefault() constexpr ClusterRef(ClusterSourceId sourceValue, uint32_t indexValue, uint16_t flagsValue = 0) noexcept
    : source{sourceValue}, flags{flagsValue}, index{indexValue}
  {
  }

  ClusterSourceId source{};
  uint16_t flags{0};
  uint32_t index{std::numeric_limits<uint32_t>::max()};

  GPUhdi() constexpr bool isValid() const noexcept { return source.isValid() && index != std::numeric_limits<uint32_t>::max(); }
  GPUhdi() friend constexpr bool operator==(ClusterRef lhs, ClusterRef rhs) noexcept { return lhs.source == rhs.source && lhs.index == rhs.index; }
  GPUhdi() friend constexpr bool operator!=(ClusterRef lhs, ClusterRef rhs) noexcept { return !(lhs == rhs); }
};

struct DetectorSensorId {
  uint32_t detector{std::numeric_limits<uint32_t>::max()};
  uint32_t sensor{std::numeric_limits<uint32_t>::max()};

  GPUhdi() constexpr bool isValid() const noexcept
  {
    return detector != std::numeric_limits<uint32_t>::max() && sensor != std::numeric_limits<uint32_t>::max();
  }
  GPUhdi() friend constexpr bool operator==(DetectorSensorId, DetectorSensorId) noexcept = default;
  GPUhdi() friend constexpr bool operator!=(DetectorSensorId, DetectorSensorId) noexcept = default;
};

struct GlobalPoint3F {
  float x{0.f};
  float y{0.f};
  float z{0.f};
};

struct GlobalCovariance3F {
  float xx{0.f};
  float xy{0.f};
  float xz{0.f};
  float yy{0.f};
  float yz{0.f};
  float zz{0.f};
};

struct ClusterShape {
  uint32_t nPixels{0};
  uint16_t rowSpan{0};
  uint16_t columnSpan{0};
};

struct GlobalMeasurement {
  GlobalPoint3F position{};
  float radius{0.f};
  GlobalCovariance3F covariance{};
  DetectorSensorId sensor{};
  ClusterRef cluster{};
  ClusterShape shape{};
  uint32_t sourceROF{std::numeric_limits<uint32_t>::max()};
  SurfaceId surface{};
  uint16_t flags{0};
};

#define O2_ITSMFT_ASSERT_GLOBAL_TYPE(Type, Size)     \
  static_assert(std::is_standard_layout_v<Type>);    \
  static_assert(std::is_trivially_copyable_v<Type>); \
  static_assert(sizeof(Type) == Size)

O2_ITSMFT_ASSERT_GLOBAL_TYPE(ClusterRef, 8);
O2_ITSMFT_ASSERT_GLOBAL_TYPE(DetectorSensorId, 8);
O2_ITSMFT_ASSERT_GLOBAL_TYPE(GlobalPoint3F, 12);
O2_ITSMFT_ASSERT_GLOBAL_TYPE(GlobalCovariance3F, 24);
O2_ITSMFT_ASSERT_GLOBAL_TYPE(ClusterShape, 8);
O2_ITSMFT_ASSERT_GLOBAL_TYPE(GlobalMeasurement, 72);

#undef O2_ITSMFT_ASSERT_GLOBAL_TYPE

static_assert(alignof(GlobalMeasurement) == 4);
static_assert(offsetof(GlobalMeasurement, position) == 0);
static_assert(offsetof(GlobalMeasurement, radius) == 12);
static_assert(offsetof(GlobalMeasurement, covariance) == 16);
static_assert(offsetof(GlobalMeasurement, sensor) == 40);
static_assert(offsetof(GlobalMeasurement, cluster) == 48);
static_assert(offsetof(GlobalMeasurement, shape) == 56);
static_assert(offsetof(GlobalMeasurement, sourceROF) == 64);
static_assert(offsetof(GlobalMeasurement, surface) == 68);
static_assert(offsetof(GlobalMeasurement, flags) == 70);

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_GLOBALMEASUREMENT_H_
