// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_GEOMETRYSURFACECATALOGPROVIDER_H_
#define ALICEO2_ITSMFT_TRACKING_GEOMETRYSURFACECATALOGPROVIDER_H_

#include <cstddef>
#include <functional>
#include <optional>

#include "DetectorSurfaceCatalogAggregation.h"
#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"

namespace o2::itsmft::tracking::detail
{

struct DetectorGeometryCatalogSpec {
  o2::detectors::DetID::ID detector;
  uint32_t surfaceCount;
  SurfaceKind kind;
  SurfaceReferenceCoordinate referenceCoordinate;
};

struct DetectorGeometryView {
  size_t chipCount{0};
  bool l2gCacheFilled{false};
  size_t l2gCacheSize{0};
  SurfaceForChip surfaceForChip{};
  LocalToGlobal localToGlobal{};
};

using DetectorGeometryLookup = std::function<std::optional<DetectorGeometryView>()>;

DetectorSurfaceCatalogResult buildGeometrySurfaceCatalog(const DetectorSurfaceCatalogRequest& request,
                                                         const DetectorGeometryCatalogSpec& spec,
                                                         const DetectorGeometryLookup& geometryLookup);

} // namespace o2::itsmft::tracking::detail

#endif /* ALICEO2_ITSMFT_TRACKING_GEOMETRYSURFACECATALOGPROVIDER_H_ */
