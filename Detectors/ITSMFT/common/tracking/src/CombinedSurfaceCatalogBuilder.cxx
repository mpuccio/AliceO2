// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/CombinedSurfaceCatalogBuilder.h"

namespace o2::itsmft::tracking
{

CombinedSurfaceCatalogResult buildCombinedSurfaceCatalog(gsl::span<const CombinedCatalogSlice> slices)
{
  CombinedSurfaceCatalogResult result;
  uint64_t offset = 0;
  for (size_t i = 0; i < slices.size(); ++i) {
    const auto& slice = slices[i];
    if (slice.provider == nullptr) {
      return {{}, DetectorSurfaceCatalogError::InvalidRequest, i};
    }
    if (offset + slice.count > MaxLayoutSurfaces) {
      return {{}, DetectorSurfaceCatalogError::InvalidRequest, i};
    }

    const DetectorSurfaceCatalogRequest request{slice.detector, SurfaceId{0}, slice.count};
    auto sliceResult = slice.provider->buildCatalog(request);
    if (!sliceResult.ok()) {
      return {{}, sliceResult.error, i};
    }
    if (sliceResult.catalog.size() != slice.count) {
      return {{}, DetectorSurfaceCatalogError::InvalidRequest, i};
    }

    for (auto descriptor : sliceResult.catalog) {
      descriptor.id = SurfaceId{static_cast<uint16_t>(descriptor.id.value() + offset)};
      result.catalog.push_back(descriptor);
    }
    offset += slice.count;
  }
  return result;
}

} // namespace o2::itsmft::tracking
