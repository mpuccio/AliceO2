// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_COMBINEDSURFACECATALOGBUILDER_H_
#define ALICEO2_ITSMFT_TRACKING_COMBINEDSURFACECATALOGBUILDER_H_

#ifndef GPUCA_GPUCODE

#include <cstddef>
#include <limits>
#include <vector>

#include <gsl/span>

#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"
#include "ITSMFTTracking/SurfaceId.h"

namespace o2::itsmft::tracking
{

// One detector's contribution to a combined, multi-detector surface catalog:
// its own existing DetectorSurfaceCatalogProvider (used exactly as it is
// used standalone today -- always queried at firstSurface=0, its own
// detector and count) plus that detector's declared surface count. Neither
// this type nor buildCombinedSurfaceCatalog() below is "owned" by any one
// detector: a combined catalog of N slices treats every slice identically,
// regardless of its position in `slices`.
struct CombinedCatalogSlice {
  const DetectorSurfaceCatalogProvider* provider{nullptr};
  o2::detectors::DetID::ID detector{o2::detectors::DetID::ITS};
  uint32_t count{0};
};

struct CombinedSurfaceCatalogResult {
  std::vector<SurfaceDescriptor> catalog{};
  DetectorSurfaceCatalogError error{DetectorSurfaceCatalogError::None};
  // Index into the `slices` span passed to buildCombinedSurfaceCatalog() that
  // failed. Meaningful only when !ok(); std::numeric_limits<size_t>::max()
  // (the default) otherwise.
  size_t failedSlice{std::numeric_limits<size_t>::max()};

  bool ok() const noexcept { return error == DetectorSurfaceCatalogError::None; }
};

// Builds one global, densely-numbered SurfaceDescriptor catalog spanning
// multiple detectors' own catalogs, by calling each slice's provider
// independently -- at firstSurface=0, exactly as that provider is invoked by
// any single-detector caller today -- and concatenating the results with a
// per-slice global SurfaceId offset equal to the running sum of the counts
// of every preceding slice. This does not modify, and is not a variant of,
// any existing DetectorSurfaceCatalogProvider::buildCatalog() call: every
// slice's provider sees the identical request shape (firstSurface=0, its own
// detector and count) it would see if queried alone; only this function's
// own output is subsequently re-numbered.
//
// Transactional: on any slice's failure (null provider, provider error, or a
// slice whose returned catalog size disagrees with its declared count), the
// whole result is empty and `error`/`failedSlice` identify the first
// offending slice; no partial catalog is ever returned. Slices whose
// cumulative count would exceed MaxLayoutSurfaces are rejected the same way,
// before any provider beyond the offending one is even queried.
//
// Only SurfaceDescriptor::id (the global identity) is shifted by the running
// offset. SurfaceDescriptor::detectorSurfaceIndex (the detector-local index)
// is left exactly as each provider returned it -- it is not, and must never
// become, a function of slice position.
CombinedSurfaceCatalogResult buildCombinedSurfaceCatalog(gsl::span<const CombinedCatalogSlice> slices);

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_COMBINEDSURFACECATALOGBUILDER_H_ */
