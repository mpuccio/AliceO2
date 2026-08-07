// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_INDEXTABLECONFIGURATION_H_
#define ALICEO2_ITSMFT_TRACKING_INDEXTABLECONFIGURATION_H_

#include <cstdint>

// Host-only: TrackingParameters (Configuration.h) owns std::vector members
// and is not itself device-compatible -- same reasoning as
// the host material/configuration binding. Kept as its own boundary, separate
// from the host binding, so that header's existing consumers
// (bindAttachHitConfig, bindLayerGeometryConfig,
// host kernel-parameter conversion) do not transitively pick up
// IndexTableUtils.h's own ITStracking/Cluster.h and MFTTracking/Constants.h
// dependencies merely because this binder also lives under
// host binding's include.
#ifndef GPUCA_GPUCODE

#include <cmath>
#include <cstddef>
#include <limits>

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"

namespace o2::itsmft::tracking
{

enum class IndexTableConfigError : uint8_t {
  None,
  NonPositiveRowBins,
  NonPositiveColBins,
  RowColBinCountExceedsIndexRange, // RowBins*ColBins would exceed int, the type used for bin indices in the unchanged hot loop
  InvalidActiveLayerCount,         // !(0 < activeSurfaceCount <= MaxLayoutSurfaces) or params.NLayers != activeSurfaceCount
  InsufficientLayerColHalfExtent,  // selected extent source shorter than activeSurfaceCount
  NonFiniteColHalfExtent,          // NaN or +/-Inf
  NonPositiveColHalfExtent,        // finite but <= 0
  NonFiniteRowRange,               // XY only: rowMin or rowMax non-finite
  DegenerateRowRange,              // XY only: finite but rowMax <= rowMin
  InvalidSurfaceKind,              // kind is neither Cylinder nor Disk
};

/// Binds and validates one iteration's TrackingParameters into `staged`,
/// keyed by the active endpoint SurfaceKind. The caller must resolve `kind`
/// from the validated SurfaceGraph/SurfacePlanBinding for this iteration --
/// never from NLayers or DetId. Every field of `params` is validated before
/// `staged` is mutated; on any error `staged` is left completely untouched.
/// Must be called once per iteration, outside any candidate/neighbour/road loop.
IndexTableConfigError bindIndexTableConfiguration(o2::itsmft::IndexTableUtilsCore& staged,
                                                  const TrackingParameters& params,
                                                  int activeSurfaceCount,
                                                  SurfaceKind kind) noexcept;

/// True iff every field setIndexTableParams stores agrees exactly between
/// `a` and `b`. Used to enforce that a non-FirstPass iteration's freshly
/// bound configuration matches the TimeFrame-owned configuration and LUT
/// content it intends to reuse or resort, before that TimeFrame is touched.
inline bool indexTableConfigurationsMatch(const o2::itsmft::IndexTableUtilsCore& a,
                                          const o2::itsmft::IndexTableUtilsCore& b,
                                          int activeSurfaceCount) noexcept
{
  if (a.getCoordType() != b.getCoordType() ||
      a.getNrowBins() != b.getNrowBins() ||
      a.getNcolBins() != b.getNcolBins() ||
      a.getRowOrigin() != b.getRowOrigin() ||
      a.getRowCoordinateSpan() != b.getRowCoordinateSpan()) {
    return false;
  }
  if (activeSurfaceCount <= 0 || activeSurfaceCount > o2::itsmft::IndexTableUtilsCore::MaxLayers) {
    return false;
  }
  for (int iLayer = 0; iLayer < activeSurfaceCount; ++iLayer) {
    if (a.getLayerColHalfExtent(iLayer) != b.getLayerColHalfExtent(iLayer)) {
      return false;
    }
  }
  return true;
}

/// Checked size_t multiplication for index-table allocation sizes: returns
/// false (leaving `result` unset) iff `a * b` would overflow size_t. A
/// division-guard, valid for any operand values -- never relies on operands
/// being "realistically" bounded.
inline bool checkedIndexTableSizeProduct(std::size_t a, std::size_t b, std::size_t& result) noexcept
{
  if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) {
    return false;
  }
  result = a * b;
  return true;
}

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_INDEXTABLECONFIGURATION_H_ */
