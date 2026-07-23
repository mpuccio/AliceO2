// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/MFTSurfaceCatalogProvider.h"

#include "GeometrySurfaceCatalogProvider.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "MathUtils/Cartesian.h"
#include "MFTBase/GeometryTGeo.h"

namespace o2::itsmft::tracking
{

DetectorSurfaceCatalogResult MFTSurfaceCatalogProvider::buildCatalog(const DetectorSurfaceCatalogRequest& request) const
{
  const detail::DetectorGeometryCatalogSpec spec{o2::detectors::DetID::MFT, MFTNLayers, SurfaceKind::Disk,
                                                 detail::SurfaceReferenceCoordinate::MeanZ,
                                                 gsl::span<const float>{kNominalMFTLayerX0.data(), kNominalMFTLayerX0.size()}};
  return detail::buildGeometrySurfaceCatalog(request, spec, []() -> std::optional<detail::DetectorGeometryView> {
    if (!o2::mft::GeometryTGeo::instanceExist()) {
      return std::nullopt;
    }
    const auto* geometry = o2::mft::GeometryTGeo::Instance();
    const auto& cache = geometry->getCacheL2G();
    return detail::DetectorGeometryView{
      static_cast<size_t>(geometry->getNumberOfChips()), cache.isFilled(), static_cast<size_t>(cache.getSize()),
      [geometry](size_t chip) { return geometry->getLayer(static_cast<int>(chip)); },
      [geometry](size_t chip, const detail::GeometryPoint& local) {
        const auto global = geometry->getMatrixL2G(static_cast<int>(chip))(o2::math_utils::Point3D<double>{local.x, local.y, local.z});
        return detail::GeometryPoint{global.X(), global.Y(), global.Z()};
      }};
  });
}

} // namespace o2::itsmft::tracking
