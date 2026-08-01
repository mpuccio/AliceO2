// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B1 Slice 1: the only file in this validator target that touches
// GeometryTGeo. Uses the same GeometryTGeo API calls
// (getNumberOfChips/getLayer/getMatrixL2G) and the same
// o2::itsmft::SegmentationAlpide active-area constants the now-deleted
// ITSSurfaceCatalogProvider.cxx/MFTSurfaceCatalogProvider.cxx/
// GeometrySurfaceCatalogProvider.cxx (Gate 4 B2 Slice 3) used for both
// detectors -- just outside O2::ITSMFTTracking's link graph and without
// depending on any of that library's catalog/provider types. No CCDB, no
// alignment: callers are
// required to have already loaded geometry via
// GeometryManager::loadGeometry(prefix, /*applyMisalignment=*/false,
// /*preferAlignedFile=*/false) and filled the selected detector's L2G matrix
// cache before calling resolveGeometrySurfaceSource().

#ifndef ALICEO2_ITSMFT_VALIDATION_GEOMETRYTGEOSURFACESOURCE_H_
#define ALICEO2_ITSMFT_VALIDATION_GEOMETRYTGEOSURFACESOURCE_H_

#include <cstddef>
#include <cstdint>

#include "DetectorSurfaceCatalogAggregation.h"
#include "ValidatorArgs.h"

namespace o2::itsmft::validation
{

struct GeometrySurfaceSource {
  size_t chipCount{0};
  o2::itsmft::tracking::detail::SurfaceForChip surfaceForChip{};
  o2::itsmft::tracking::detail::LocalToGlobal localToGlobal{};
  o2::itsmft::tracking::detail::SurfaceReferenceCoordinate referenceCoordinate{
    o2::itsmft::tracking::detail::SurfaceReferenceCoordinate::MeanRadius};
};

enum class GeometrySourceError : uint8_t {
  None,
  GeometryNotInitialized,
  L2GCacheNotFilled
};

struct GeometrySurfaceSourceResult {
  GeometrySurfaceSource source{};
  GeometrySourceError error{GeometrySourceError::None};

  bool ok() const noexcept { return error == GeometrySourceError::None; }
};

// Reads the already-populated GeometryTGeo singleton for `detector` (ITS or
// MFT); does not itself load geometry or fill the matrix cache.
GeometrySurfaceSourceResult resolveGeometrySurfaceSource(DetectorSelection detector);

// The standard ALPIDE active footprint (o2::itsmft::SegmentationAlpide
// constants), the same one the now-deleted GeometrySurfaceCatalogProvider.cxx
// used for both ITS and MFT -- not a new or invented value.
o2::itsmft::tracking::detail::LocalActiveArea alpideActiveArea();

} // namespace o2::itsmft::validation

#endif /* ALICEO2_ITSMFT_VALIDATION_GEOMETRYTGEOSURFACESOURCE_H_ */
