// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "GeometryTGeoSurfaceSource.h"

#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "MFTBase/GeometryTGeo.h"
#include "MathUtils/Cartesian.h"

namespace o2::itsmft::validation
{

using o2::itsmft::tracking::detail::GeometryPoint;
using o2::itsmft::tracking::detail::LocalActiveArea;
using o2::itsmft::tracking::detail::SurfaceReferenceCoordinate;

LocalActiveArea alpideActiveArea()
{
  using Segmentation = o2::itsmft::SegmentationAlpide;
  return LocalActiveArea{-0.5 * Segmentation::SensorSizeRows + Segmentation::PassiveEdgeReadOut,
                         0.5 * Segmentation::SensorSizeRows - Segmentation::PassiveEdgeTop,
                         -0.5 * Segmentation::ActiveMatrixSizeCols,
                         0.5 * Segmentation::ActiveMatrixSizeCols};
}

namespace
{
template <typename GeometryTGeoT>
GeometrySurfaceSourceResult resolveFrom(SurfaceReferenceCoordinate referenceCoordinate)
{
  if (!GeometryTGeoT::instanceExist()) {
    return {{}, GeometrySourceError::GeometryNotInitialized};
  }
  const auto* geometry = GeometryTGeoT::Instance();
  const auto& cache = geometry->getCacheL2G();
  const auto chipCount = static_cast<size_t>(geometry->getNumberOfChips());
  if (!cache.isFilled() || static_cast<size_t>(cache.getSize()) < chipCount) {
    return {{}, GeometrySourceError::L2GCacheNotFilled};
  }

  GeometrySurfaceSource source;
  source.chipCount = chipCount;
  source.referenceCoordinate = referenceCoordinate;
  source.surfaceForChip = [geometry](size_t chip) { return geometry->getLayer(static_cast<int>(chip)); };
  source.localToGlobal = [geometry](size_t chip, const GeometryPoint& local) {
    const auto global = geometry->getMatrixL2G(static_cast<int>(chip))(o2::math_utils::Point3D<double>{local.x, local.y, local.z});
    return GeometryPoint{global.X(), global.Y(), global.Z()};
  };
  return {source, GeometrySourceError::None};
}
} // namespace

GeometrySurfaceSourceResult resolveGeometrySurfaceSource(DetectorSelection detector)
{
  if (detector == DetectorSelection::ITS) {
    return resolveFrom<o2::its::GeometryTGeo>(SurfaceReferenceCoordinate::MeanRadius);
  }
  return resolveFrom<o2::mft::GeometryTGeo>(SurfaceReferenceCoordinate::MeanZ);
}

} // namespace o2::itsmft::validation
