// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_MEASUREMENTVIEW_H_
#define ALICEO2_ITSMFT_TRACKING_MEASUREMENTVIEW_H_

#include <cstdint>
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/GlobalMeasurement.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/SurfaceTiming.h"

namespace o2::itsmft::tracking
{

struct SurfaceMeasurementSpan {
  const GlobalMeasurement* global{nullptr};
  const SurfaceMeasurement* data{nullptr};
  uint32_t count{0};
};

static_assert(std::is_standard_layout_v<SurfaceMeasurementSpan>);
static_assert(std::is_trivially_copyable_v<SurfaceMeasurementSpan>);

// GPU-safe read-only view over TimeFrame's event measurement storage.
struct MeasurementView {
  const SurfaceMeasurementSpan* surfaces{nullptr};
  uint32_t nSurfaces{0};
  const ROFIntervalBC* rofIntervals{nullptr};
  const uint32_t* sourceROFOffsets{nullptr};
  uint32_t nSources{0};

  GPUhdi() const SurfaceMeasurement* getSurfaceMeasurements(SurfaceId surface) const noexcept
  {
    const auto& span = surfaces[surface.value()];
    return span.count == 0 ? nullptr : span.data;
  }
  GPUhdi() const GlobalMeasurement* getGlobalMeasurements(SurfaceId surface) const noexcept
  {
    const auto& span = surfaces[surface.value()];
    return span.count == 0 ? nullptr : span.global;
  }
  GPUhdi() uint32_t getSurfaceMeasurementCount(SurfaceId surface) const noexcept { return surfaces[surface.value()].count; }
  GPUhdi() const GlobalMeasurement* getGlobalMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept
  {
    if (!surface.isValid() || surface.value() >= nSurfaces || !index.isValid()) {
      return nullptr;
    }
    const auto& span = surfaces[surface.value()];
    return index.value() < span.count ? span.global + index.value() : nullptr;
  }
  GPUhdi() const SurfaceMeasurement* getSurfaceMeasurement(SurfaceId surface, SurfaceMeasurementIndex index) const noexcept
  {
    if (!surface.isValid() || surface.value() >= nSurfaces || !index.isValid()) {
      return nullptr;
    }
    const auto& span = surfaces[surface.value()];
    return index.value() < span.count ? span.data + index.value() : nullptr;
  }
  GPUhdi() const ROFIntervalBC* getSourceIntervals(ClusterSourceId source) const noexcept
  {
    const auto first = sourceROFOffsets[source.value()];
    const auto last = sourceROFOffsets[source.value() + 1];
    return first == last ? nullptr : rofIntervals + first;
  }
  GPUhdi() uint32_t getSourceIntervalCount(ClusterSourceId source) const noexcept
  {
    return sourceROFOffsets[source.value() + 1] - sourceROFOffsets[source.value()];
  }
};

static_assert(std::is_standard_layout_v<MeasurementView>);
static_assert(std::is_trivially_copyable_v<MeasurementView>);

} // namespace o2::itsmft::tracking

#endif
