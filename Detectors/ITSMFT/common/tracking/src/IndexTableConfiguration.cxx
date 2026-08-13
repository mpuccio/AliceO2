// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/IndexTableConfiguration.h"

#include <cmath>
#include <cstdint>
#include <limits>

#include "CommonConstants/MathConstants.h"

namespace o2::itsmft::tracking
{

using o2::itsmft::IndexTableCoordType;

IndexTableConfigError bindIndexTableConfiguration(o2::itsmft::IndexTableUtilsCore& staged,
                                                  const TrackingParameters& params,
                                                  int activeSurfaceCount,
                                                  SurfaceKind kind) noexcept
{
  if (kind != SurfaceKind::Cylinder && kind != SurfaceKind::Disk) {
    return IndexTableConfigError::InvalidSurfaceKind;
  }
  if (!(activeSurfaceCount > 0 && activeSurfaceCount <= o2::itsmft::IndexTableUtilsCore::MaxLayers) ||
      params.NLayers != activeSurfaceCount) {
    return IndexTableConfigError::InvalidActiveLayerCount;
  }
  if (params.RowBins <= 0) {
    return IndexTableConfigError::NonPositiveRowBins;
  }
  if (params.ColBins <= 0) {
    return IndexTableConfigError::NonPositiveColBins;
  }

  const std::uint64_t binCount = static_cast<std::uint64_t>(params.RowBins) * static_cast<std::uint64_t>(params.ColBins);
  if (binCount > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    return IndexTableConfigError::RowColBinCountExceedsIndexRange;
  }

  float rowMin = 0.f;
  float rowMax = o2::constants::math::TwoPI;
  if (kind == SurfaceKind::Disk) {
    constexpr float defaultRowMin{-20.f};
    constexpr float defaultRowMax{20.f};
    const bool hasRowRange = params.IndexRowMax != 0.f;
    rowMin = hasRowRange ? params.IndexRowMin : defaultRowMin;
    rowMax = hasRowRange ? params.IndexRowMax : defaultRowMax;
    if (!std::isfinite(rowMin) || !std::isfinite(rowMax)) {
      return IndexTableConfigError::NonFiniteRowRange;
    }
    if (!(rowMax > rowMin)) {
      return IndexTableConfigError::DegenerateRowRange;
    }
  }

  const bool useColHalfExtent = !params.LayerColHalfExtent.empty();
  const auto& source = useColHalfExtent ? params.LayerColHalfExtent : params.LayerZ;
  if (source.size() < static_cast<std::size_t>(activeSurfaceCount)) {
    return IndexTableConfigError::InsufficientLayerColHalfExtent;
  }

  for (int iLayer = 0; iLayer < activeSurfaceCount; ++iLayer) {
    const float extent = source[iLayer];
    if (!std::isfinite(extent)) {
      return IndexTableConfigError::NonFiniteColHalfExtent;
    }
    if (!(extent > 0.f)) {
      return IndexTableConfigError::NonPositiveColHalfExtent;
    }
  }

  // Mutate `staged` only after validation. Use the low-level setter because
  // the convenience wrappers silently zero-fill short extent sources instead
  // of rejecting them; `source` is already known to cover the active layers.
  staged.setIndexTableParams(kind == SurfaceKind::Disk ? IndexTableCoordType::XY : IndexTableCoordType::PhiZ,
                             params.RowBins, params.ColBins, rowMin, rowMax,
                             gsl::span<const float>{source.data(), static_cast<std::size_t>(activeSurfaceCount)});
  return IndexTableConfigError::None;
}

} // namespace o2::itsmft::tracking
