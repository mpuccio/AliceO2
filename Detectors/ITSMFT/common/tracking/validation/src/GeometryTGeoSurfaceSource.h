// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Accesses GeometryTGeo directly for ITS and MFT validation. Callers must
// load geometry without alignment and fill the selected detector's L2G cache
// before calling resolveGeometrySurfaceSource().

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

// Reads the populated GeometryTGeo singleton; does not load geometry or fill
// the matrix cache.
GeometrySurfaceSourceResult resolveGeometrySurfaceSource(DetectorSelection detector);

// Returns the standard ALPIDE active footprint from SegmentationAlpide.
o2::itsmft::tracking::detail::LocalActiveArea alpideActiveArea();

} // namespace o2::itsmft::validation

#endif /* ALICEO2_ITSMFT_VALIDATION_GEOMETRYTGEOSURFACESOURCE_H_ */
