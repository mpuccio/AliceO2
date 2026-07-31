// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "GeometrySurfaceCatalogProvider.h"

#include "ITSMFTBase/SegmentationAlpide.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking::detail
{
namespace
{
DetectorSurfaceCatalogError mapAggregationError(SurfaceGeometryAggregationError error)
{
  switch (error) {
    case SurfaceGeometryAggregationError::SurfaceOutOfRange:
    case SurfaceGeometryAggregationError::EmptySurface:
      return DetectorSurfaceCatalogError::SurfaceLookupFailure;
    case SurfaceGeometryAggregationError::NonFiniteCoordinate:
    case SurfaceGeometryAggregationError::InvalidDerivedGeometry:
      return DetectorSurfaceCatalogError::InvalidSurfaceGeometry;
    case SurfaceGeometryAggregationError::None:
      return DetectorSurfaceCatalogError::None;
  }
  return DetectorSurfaceCatalogError::InvalidSurfaceGeometry;
}
} // namespace

DetectorSurfaceCatalogResult buildGeometrySurfaceCatalog(const DetectorSurfaceCatalogRequest& request,
                                                         const DetectorGeometryCatalogSpec& spec,
                                                         const DetectorGeometryLookup& geometryLookup)
{
  if (request.detector != spec.detector) {
    return {{}, DetectorSurfaceCatalogError::UnsupportedDetector};
  }
  if (request.firstSurface != SurfaceId{0} || request.detectorSurfaceCount != spec.surfaceCount) {
    return {{}, DetectorSurfaceCatalogError::InvalidRequest};
  }

  const auto geometry = geometryLookup();
  if (!geometry) {
    return {{}, DetectorSurfaceCatalogError::GeometryNotInitialized};
  }
  if (!geometry->l2gCacheFilled || geometry->l2gCacheSize < geometry->chipCount) {
    return {{}, DetectorSurfaceCatalogError::GeometryUnavailable};
  }

  using Segmentation = o2::itsmft::SegmentationAlpide;
  const LocalActiveArea activeArea{-0.5 * Segmentation::SensorSizeRows + Segmentation::PassiveEdgeReadOut,
                                   0.5 * Segmentation::SensorSizeRows - Segmentation::PassiveEdgeTop,
                                   -0.5 * Segmentation::ActiveMatrixSizeCols,
                                   0.5 * Segmentation::ActiveMatrixSizeCols};
  auto aggregation = aggregateSurfaceGeometry(geometry->chipCount, spec.surfaceCount, activeArea,
                                              spec.referenceCoordinate, geometry->surfaceForChip,
                                              geometry->localToGlobal);
  if (!aggregation.ok()) {
    return {{}, mapAggregationError(aggregation.error)};
  }

  // Before indexing: the nominal-material span supplied by the concrete
  // provider must cover every catalog surface.
  if (spec.nominalLayerX0.size() != spec.surfaceCount) {
    return {{}, DetectorSurfaceCatalogError::InvalidRequest};
  }

  DetectorSurfaceCatalogResult result;
  result.catalog.reserve(spec.surfaceCount);
  for (uint16_t surface = 0; surface < spec.surfaceCount; ++surface) {
    const auto& geometryValues = aggregation.surfaces[surface];
    SurfaceDescriptor descriptor{SurfaceId{surface}, surface,
                                 static_cast<uint8_t>(spec.detector), spec.kind, 0,
                                 geometryValues.referenceCoordinate};
    descriptor.material.xOverX0 = spec.nominalLayerX0[surface];
    descriptor.material.arealDensityGPerCm2 = spec.nominalLayerX0[surface] *
                                              o2::its::constants::Radl *
                                              o2::its::constants::Rho;
    result.catalog.push_back(descriptor);
  }
  return result;
}

} // namespace o2::itsmft::tracking::detail
