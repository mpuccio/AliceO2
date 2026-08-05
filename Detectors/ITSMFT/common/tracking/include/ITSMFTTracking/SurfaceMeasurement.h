// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENT_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENT_H_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include <gsl/span>

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
  GPUhdi() friend constexpr bool operator==(DetectorSensorId lhs, DetectorSensorId rhs) noexcept
  {
    return lhs.detector == rhs.detector && lhs.sensor == rhs.sensor;
  }
  GPUhdi() friend constexpr bool operator!=(DetectorSensorId lhs, DetectorSensorId rhs) noexcept { return !(lhs == rhs); }
};

struct GlobalPoint3F {
  float x{0.f};
  float y{0.f};
  float z{0.f};
};

// q is normal to the surface. The measured coordinates are always (u, v).
struct SurfaceFramePoint {
  float q{0.f};
  float u{0.f};
  float v{0.f};
  float frameAngle{0.f};
};

// Packed symmetric covariance of the measured (u, v) coordinates.
struct SurfaceCovariance2F {
  float uu{0.f};
  float uv{0.f};
  float vv{0.f};
};

struct ClusterShape {
  uint32_t nPixels{0};
  uint16_t rowSpan{0};
  uint16_t columnSpan{0};
};

struct SurfaceMeasurement {
  GlobalPoint3F global{};
  SurfaceFramePoint frame{};
  SurfaceCovariance2F covariance{};
  DetectorSensorId sensor{};
  ClusterRef cluster{};
  ClusterShape shape{};
  uint32_t sourceROF{std::numeric_limits<uint32_t>::max()};
  SurfaceId surface{};
  uint16_t flags{0};
};

#define O2_ITSMFT_ASSERT_DEVICE_TYPE(Type, Size)     \
  static_assert(std::is_standard_layout_v<Type>);    \
  static_assert(std::is_trivially_copyable_v<Type>); \
  static_assert(sizeof(Type) == Size)

O2_ITSMFT_ASSERT_DEVICE_TYPE(ClusterRef, 8);
O2_ITSMFT_ASSERT_DEVICE_TYPE(DetectorSensorId, 8);
O2_ITSMFT_ASSERT_DEVICE_TYPE(GlobalPoint3F, 12);
O2_ITSMFT_ASSERT_DEVICE_TYPE(SurfaceFramePoint, 16);
O2_ITSMFT_ASSERT_DEVICE_TYPE(SurfaceCovariance2F, 12);
O2_ITSMFT_ASSERT_DEVICE_TYPE(ClusterShape, 8);
O2_ITSMFT_ASSERT_DEVICE_TYPE(SurfaceMeasurement, 72);
static_assert(offsetof(ClusterRef, source) == 0);
static_assert(offsetof(ClusterRef, flags) == 2);
static_assert(offsetof(ClusterRef, index) == 4);
static_assert(alignof(SurfaceMeasurement) == 4);
static_assert(offsetof(SurfaceMeasurement, global) == 0);
static_assert(offsetof(SurfaceMeasurement, frame) == 12);
static_assert(offsetof(SurfaceMeasurement, covariance) == 28);
static_assert(offsetof(SurfaceMeasurement, sensor) == 40);
static_assert(offsetof(SurfaceMeasurement, cluster) == 48);
static_assert(offsetof(SurfaceMeasurement, shape) == 56);
static_assert(offsetof(SurfaceMeasurement, sourceROF) == 64);
static_assert(offsetof(SurfaceMeasurement, surface) == 68);
static_assert(offsetof(SurfaceMeasurement, flags) == 70);

#undef O2_ITSMFT_ASSERT_DEVICE_TYPE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACEMEASUREMENT_H_ */
