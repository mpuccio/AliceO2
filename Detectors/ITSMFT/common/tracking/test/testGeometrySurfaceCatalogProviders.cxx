// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT geometry surface catalog providers
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <array>
#include <limits>

#include "../src/GeometrySurfaceCatalogProvider.h"
#include "ITSMFTTracking/ITSSurfaceCatalogProvider.h"
#include "ITSMFTTracking/MFTSurfaceCatalogProvider.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/Constants.h"

using namespace o2::itsmft::tracking;
using namespace o2::itsmft::tracking::detail;

namespace
{
const gsl::span<const float> itsNominalLayerX0{kNominalITSLayerX0.data(), kNominalITSLayerX0.size()};
const gsl::span<const float> mftNominalLayerX0{kNominalMFTLayerX0.data(), kNominalMFTLayerX0.size()};
const DetectorGeometryCatalogSpec itsSpec{o2::detectors::DetID::ITS, ITSNLayers, SurfaceKind::Cylinder,
                                          SurfaceReferenceCoordinate::MeanRadius, itsNominalLayerX0};
const DetectorGeometryCatalogSpec mftSpec{o2::detectors::DetID::MFT, MFTNLayers, SurfaceKind::Disk,
                                          SurfaceReferenceCoordinate::MeanZ, mftNominalLayerX0};

DetectorSurfaceCatalogRequest request(o2::detectors::DetID::ID detector, uint16_t first, uint32_t count)
{
  return {detector, SurfaceId{first}, count};
}

DetectorGeometryView geometry(size_t chips, size_t surfaces)
{
  return DetectorGeometryView{
    chips, true, chips,
    [surfaces](size_t chip) { return static_cast<int>(chip % surfaces); },
    [](size_t chip, const GeometryPoint& local) {
      return GeometryPoint{10. + static_cast<double>(chip) + local.x, local.z,
                           -50. - static_cast<double>(chip) + 0.1 * local.x};
    }};
}

void checkEmptyFailure(const DetectorSurfaceCatalogResult& result, DetectorSurfaceCatalogError error)
{
  BOOST_CHECK(result.error == error);
  BOOST_CHECK(result.catalog.empty());
}
} // namespace

BOOST_AUTO_TEST_CASE(fixed_request_contracts_precede_geometry_lookup)
{
  bool lookedUp = false;
  const DetectorGeometryLookup lookup = [&]() -> std::optional<DetectorGeometryView> {
    lookedUp = true;
    return std::nullopt;
  };
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::MFT, 0, ITSNLayers), itsSpec, lookup),
                    DetectorSurfaceCatalogError::UnsupportedDetector);
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 1, ITSNLayers), itsSpec, lookup),
                    DetectorSurfaceCatalogError::InvalidRequest);
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, ITSNLayers - 1), itsSpec, lookup),
                    DetectorSurfaceCatalogError::InvalidRequest);
  BOOST_CHECK(!lookedUp);

  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, MFTNLayers), mftSpec, lookup),
                    DetectorSurfaceCatalogError::UnsupportedDetector);
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::MFT, 1, MFTNLayers), mftSpec, lookup),
                    DetectorSurfaceCatalogError::InvalidRequest);
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::MFT, 0, MFTNLayers - 1), mftSpec, lookup),
                    DetectorSurfaceCatalogError::InvalidRequest);
  BOOST_CHECK(!lookedUp);
}

BOOST_AUTO_TEST_CASE(geometry_initialization_and_cache_failures)
{
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, ITSNLayers), itsSpec,
                                                []() -> std::optional<DetectorGeometryView> { return std::nullopt; }),
                    DetectorSurfaceCatalogError::GeometryNotInitialized);

  auto absent = geometry(ITSNLayers, ITSNLayers);
  absent.l2gCacheFilled = false;
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, ITSNLayers), itsSpec,
                                                [absent] { return absent; }),
                    DetectorSurfaceCatalogError::GeometryUnavailable);

  auto undersized = geometry(ITSNLayers, ITSNLayers);
  undersized.l2gCacheSize = ITSNLayers - 1;
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, ITSNLayers), itsSpec,
                                                [undersized] { return undersized; }),
                    DetectorSurfaceCatalogError::GeometryUnavailable);
}

BOOST_AUTO_TEST_CASE(surface_lookup_failures)
{
  auto outOfRange = geometry(ITSNLayers, ITSNLayers);
  outOfRange.surfaceForChip = [](size_t) { return ITSNLayers; };
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, ITSNLayers), itsSpec,
                                                [outOfRange] { return outOfRange; }),
                    DetectorSurfaceCatalogError::SurfaceLookupFailure);

  const auto emptyBucket = geometry(ITSNLayers - 1, ITSNLayers);
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, ITSNLayers), itsSpec,
                                                [emptyBucket] { return emptyBucket; }),
                    DetectorSurfaceCatalogError::SurfaceLookupFailure);
}

BOOST_AUTO_TEST_CASE(invalid_surface_geometry_is_mapped_and_transactional)
{
  auto invalid = geometry(ITSNLayers, ITSNLayers);
  invalid.localToGlobal = [](size_t, const GeometryPoint&) {
    return GeometryPoint{std::numeric_limits<double>::quiet_NaN(), 0., 0.};
  };
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, ITSNLayers), itsSpec,
                                                [invalid] { return invalid; }),
                    DetectorSurfaceCatalogError::InvalidSurfaceGeometry);
}

BOOST_AUTO_TEST_CASE(successful_its_and_mft_catalogs)
{
  const auto itsGeometry = geometry(ITSNLayers, ITSNLayers);
  const auto its = buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, ITSNLayers), itsSpec,
                                               [itsGeometry] { return itsGeometry; });
  BOOST_REQUIRE(its.ok());
  BOOST_REQUIRE_EQUAL(its.catalog.size(), ITSNLayers);
  for (uint16_t i = 0; i < its.catalog.size(); ++i) {
    BOOST_CHECK(its.catalog[i].id == SurfaceId{i});
    BOOST_CHECK_EQUAL(its.catalog[i].detectorSurfaceIndex, i);
    BOOST_CHECK_EQUAL(its.catalog[i].detectorId, static_cast<uint8_t>(o2::detectors::DetID::ITS));
    BOOST_CHECK(its.catalog[i].kind == SurfaceKind::Cylinder);
    BOOST_CHECK_EQUAL(its.catalog[i].flags, 0);
    BOOST_CHECK_LE(its.catalog[i].radialMin, its.catalog[i].radialMax);
    BOOST_CHECK_CLOSE(its.catalog[i].material.xOverX0, kNominalITSLayerX0[i], 1.e-6f);
    BOOST_CHECK_CLOSE(its.catalog[i].material.arealDensityGPerCm2,
                      kNominalITSLayerX0[i] * o2::its::constants::Radl * o2::its::constants::Rho, 1.e-6f);
  }

  const auto mftGeometry = geometry(MFTNLayers, MFTNLayers);
  const auto mft = buildGeometrySurfaceCatalog(request(o2::detectors::DetID::MFT, 0, MFTNLayers), mftSpec,
                                               [mftGeometry] { return mftGeometry; });
  BOOST_REQUIRE(mft.ok());
  BOOST_REQUIRE_EQUAL(mft.catalog.size(), MFTNLayers);
  for (uint16_t i = 0; i < mft.catalog.size(); ++i) {
    BOOST_CHECK(mft.catalog[i].id == SurfaceId{i});
    BOOST_CHECK_EQUAL(mft.catalog[i].detectorSurfaceIndex, i);
    BOOST_CHECK_EQUAL(mft.catalog[i].detectorId, static_cast<uint8_t>(o2::detectors::DetID::MFT));
    BOOST_CHECK(mft.catalog[i].kind == SurfaceKind::Disk);
    BOOST_CHECK_EQUAL(mft.catalog[i].flags, 0);
    BOOST_CHECK_LT(mft.catalog[i].referenceCoordinate, 0.f);
    BOOST_CHECK_CLOSE(mft.catalog[i].material.xOverX0, kNominalMFTLayerX0[i], 1.e-6f);
    BOOST_CHECK_CLOSE(mft.catalog[i].material.arealDensityGPerCm2,
                      kNominalMFTLayerX0[i] * o2::its::constants::Radl * o2::its::constants::Rho, 1.e-6f);
  }
}

BOOST_AUTO_TEST_CASE(its_and_mft_use_the_same_xoverx0_radl_rho_conversion)
{
  const auto itsGeometry = geometry(ITSNLayers, ITSNLayers);
  const auto its = buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, ITSNLayers), itsSpec,
                                               [itsGeometry] { return itsGeometry; });
  BOOST_REQUIRE(its.ok());
  const auto mftGeometry = geometry(MFTNLayers, MFTNLayers);
  const auto mft = buildGeometrySurfaceCatalog(request(o2::detectors::DetID::MFT, 0, MFTNLayers), mftSpec,
                                               [mftGeometry] { return mftGeometry; });
  BOOST_REQUIRE(mft.ok());

  for (const auto& surface : its.catalog) {
    BOOST_CHECK_CLOSE(surface.material.arealDensityGPerCm2,
                      surface.material.xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho, 1.e-6f);
  }
  for (const auto& surface : mft.catalog) {
    BOOST_CHECK_CLOSE(surface.material.arealDensityGPerCm2,
                      surface.material.xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho, 1.e-6f);
  }
}

BOOST_AUTO_TEST_CASE(provider_rejects_wrong_sized_nominal_material_span_before_indexing)
{
  const std::array<float, ITSNLayers - 1> undersized{};
  const DetectorGeometryCatalogSpec badSpec{o2::detectors::DetID::ITS, ITSNLayers, SurfaceKind::Cylinder,
                                            SurfaceReferenceCoordinate::MeanRadius,
                                            gsl::span<const float>{undersized.data(), undersized.size()}};
  bool lookedUp = false;
  const auto itsGeometry = geometry(ITSNLayers, ITSNLayers);
  const auto result = buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, ITSNLayers), badSpec,
                                                  [&]() -> std::optional<DetectorGeometryView> {
                                                    lookedUp = true;
                                                    return itsGeometry;
                                                  });
  BOOST_CHECK(lookedUp);
  checkEmptyFailure(result, DetectorSurfaceCatalogError::InvalidRequest);
}

BOOST_AUTO_TEST_CASE(public_providers_do_not_create_missing_singletons_for_bad_requests)
{
  ITSSurfaceCatalogProvider its;
  checkEmptyFailure(its.buildCatalog(request(o2::detectors::DetID::MFT, 0, ITSNLayers)),
                    DetectorSurfaceCatalogError::UnsupportedDetector);
  checkEmptyFailure(its.buildCatalog(request(o2::detectors::DetID::ITS, 1, ITSNLayers)),
                    DetectorSurfaceCatalogError::InvalidRequest);
  checkEmptyFailure(its.buildCatalog(request(o2::detectors::DetID::ITS, 0, ITSNLayers)),
                    DetectorSurfaceCatalogError::GeometryNotInitialized);

  MFTSurfaceCatalogProvider mft;
  checkEmptyFailure(mft.buildCatalog(request(o2::detectors::DetID::ITS, 0, MFTNLayers)),
                    DetectorSurfaceCatalogError::UnsupportedDetector);
  checkEmptyFailure(mft.buildCatalog(request(o2::detectors::DetID::MFT, 0, MFTNLayers - 1)),
                    DetectorSurfaceCatalogError::InvalidRequest);
  checkEmptyFailure(mft.buildCatalog(request(o2::detectors::DetID::MFT, 0, MFTNLayers)),
                    DetectorSurfaceCatalogError::GeometryNotInitialized);
}
