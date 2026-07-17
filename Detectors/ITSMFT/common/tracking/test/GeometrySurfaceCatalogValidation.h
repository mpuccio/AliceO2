// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_GEOMETRYSURFACECATALOGVALIDATION_H_
#define ALICEO2_ITSMFT_TRACKING_GEOMETRYSURFACECATALOGVALIDATION_H_

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "DetectorsBase/GeometryManager.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTTracking/ITSSurfaceCatalogProvider.h"
#include "ITSMFTTracking/MFTSurfaceCatalogProvider.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"
#include "MFTBase/GeometryTGeo.h"

namespace o2::itsmft::tracking::test
{
namespace detail
{
inline void require(bool condition, const std::string& message)
{
  if (!condition) {
    throw std::runtime_error(message);
  }
}

inline const char* catalogErrorName(DetectorSurfaceCatalogError error)
{
  switch (error) {
    case DetectorSurfaceCatalogError::None:
      return "None";
    case DetectorSurfaceCatalogError::InvalidRequest:
      return "InvalidRequest";
    case DetectorSurfaceCatalogError::UnsupportedDetector:
      return "UnsupportedDetector";
    case DetectorSurfaceCatalogError::GeometryUnavailable:
      return "GeometryUnavailable";
    case DetectorSurfaceCatalogError::GeometryNotInitialized:
      return "GeometryNotInitialized";
    case DetectorSurfaceCatalogError::SurfaceLookupFailure:
      return "SurfaceLookupFailure";
    case DetectorSurfaceCatalogError::InvalidSurfaceGeometry:
      return "InvalidSurfaceGeometry";
  }
  return "Unknown";
}

inline void printAndValidate(const char* name,
                             const DetectorSurfaceCatalogResult& result,
                             o2::detectors::DetID::ID detector,
                             SurfaceKind kind,
                             size_t expectedSize)
{
  require(result.ok(), std::string{name} + " provider error " + catalogErrorName(result.error));
  require(result.catalog.size() == expectedSize, std::string{name} + " catalog size mismatch");
  for (size_t i = 0; i < result.catalog.size(); ++i) {
    const auto& surface = result.catalog[i];
    require(surface.id == SurfaceId{static_cast<uint16_t>(i)}, std::string{name} + " non-dense surface id");
    require(surface.detectorSurfaceIndex == i, std::string{name} + " detector surface index mismatch");
    require(surface.detectorId == static_cast<uint8_t>(detector), std::string{name} + " detector mismatch");
    require(surface.kind == kind, std::string{name} + " surface-kind mismatch");
    require(surface.flags == 0, std::string{name} + " nonzero flags");
    require(std::isfinite(surface.referenceCoordinate) && std::isfinite(surface.radialMin) &&
              std::isfinite(surface.radialMax),
            std::string{name} + " non-finite surface geometry");
    require(surface.radialMin >= 0.f && surface.radialMin <= surface.radialMax,
            std::string{name} + " unordered radial bounds");
    std::cout << "SURFACE " << name << ' ' << i << " reference=" << surface.referenceCoordinate
              << " radialMin=" << surface.radialMin << " radialMax=" << surface.radialMax << '\n';
  }
}
} // namespace detail

inline void validateGeometrySurfaceCatalogs(const std::string& geometryPrefix)
{
  detail::require(!geometryPrefix.empty(), "geometry prefix is empty");
  o2::base::GeometryManager::loadGeometry(geometryPrefix);

  constexpr int l2gMask = o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G);
  auto* itsGeometry = o2::its::GeometryTGeo::Instance();
  auto* mftGeometry = o2::mft::GeometryTGeo::Instance();
  detail::require(itsGeometry != nullptr, "ITS geometry singleton is unavailable");
  detail::require(mftGeometry != nullptr, "MFT geometry singleton is unavailable");
  itsGeometry->fillMatrixCache(l2gMask);
  mftGeometry->fillMatrixCache(l2gMask);

  detail::require(itsGeometry->getNumberOfChips() > 0, "ITS geometry contains no chips");
  detail::require(mftGeometry->getNumberOfChips() > 0, "MFT geometry contains no chips");
  detail::require(itsGeometry->getCacheL2G().getSize() >= itsGeometry->getNumberOfChips(), "ITS L2G cache undersized");
  detail::require(mftGeometry->getCacheL2G().getSize() >= mftGeometry->getNumberOfChips(), "MFT L2G cache undersized");

  ITSSurfaceCatalogProvider itsProvider;
  MFTSurfaceCatalogProvider mftProvider;
  const auto its = itsProvider.buildCatalog(
    {o2::detectors::DetID::ITS, SurfaceId{0}, static_cast<uint32_t>(ITSNLayers)});
  const auto mft = mftProvider.buildCatalog(
    {o2::detectors::DetID::MFT, SurfaceId{0}, static_cast<uint32_t>(MFTNLayers)});
  detail::printAndValidate("ITS", its, o2::detectors::DetID::ITS, SurfaceKind::Cylinder, ITSNLayers);
  detail::printAndValidate("MFT", mft, o2::detectors::DetID::MFT, SurfaceKind::Disk, MFTNLayers);

  const std::array<float, ITSNLayers> itsRadiusMin{1.5f, 2.5f, 3.2f, 18.f, 23.f, 33.f, 38.f};
  const std::array<float, ITSNLayers> itsRadiusMax{3.f, 4.f, 5.f, 21.f, 27.f, 36.f, 42.f};
  for (size_t i = 0; i < its.catalog.size(); ++i) {
    detail::require(its.catalog[i].referenceCoordinate > itsRadiusMin[i] &&
                      its.catalog[i].referenceCoordinate < itsRadiusMax[i],
                    "ITS layer reference radius outside documented plausibility window");
  }
  for (size_t i = 1; i < mft.catalog.size(); ++i) {
    detail::require(mft.catalog[i].referenceCoordinate < mft.catalog[i - 1].referenceCoordinate,
                    "MFT disk-plane z coordinates are not ordered away from the interaction point");
  }
  for (const auto& surface : mft.catalog) {
    detail::require(surface.referenceCoordinate < -40.f && surface.referenceCoordinate > -85.f,
                    "MFT disk-plane z outside documented plausibility window");
    detail::require(surface.radialMin > 1.f && surface.radialMax < 25.f,
                    "MFT radial envelope outside documented plausibility window");
  }

  // MFT's geometry construction and L2G extraction both use the common
  // SegmentationAlpide sensor definition, so the standard asymmetric ALPIDE
  // active footprint used by the provider is the detector's actual footprint.
  std::cout << "ALPIDE_FOOTPRINT MFT=standard-common-SegmentationAlpide\n";
  std::cout << "VALIDATION_OK ITS=" << ITSNLayers << " MFT=" << MFTNLayers << '\n';
}

} // namespace o2::itsmft::tracking::test

#endif /* ALICEO2_ITSMFT_TRACKING_GEOMETRYSURFACECATALOGVALIDATION_H_ */
