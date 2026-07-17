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

#include <limits>

#include "../src/GeometrySurfaceCatalogProvider.h"
#include "ITSMFTTracking/ITSSurfaceCatalogProvider.h"
#include "ITSMFTTracking/MFTSurfaceCatalogProvider.h"

using namespace o2::itsmft::tracking;
using namespace o2::itsmft::tracking::detail;

namespace
{
const DetectorGeometryCatalogSpec itsSpec{o2::detectors::DetID::ITS, 7, SurfaceKind::Cylinder,
                                          SurfaceReferenceCoordinate::MeanRadius};
const DetectorGeometryCatalogSpec mftSpec{o2::detectors::DetID::MFT, 10, SurfaceKind::Disk,
                                          SurfaceReferenceCoordinate::MeanZ};

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
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::MFT, 0, 7), itsSpec, lookup),
                    DetectorSurfaceCatalogError::UnsupportedDetector);
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 1, 7), itsSpec, lookup),
                    DetectorSurfaceCatalogError::InvalidRequest);
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, 6), itsSpec, lookup),
                    DetectorSurfaceCatalogError::InvalidRequest);
  BOOST_CHECK(!lookedUp);

  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, 10), mftSpec, lookup),
                    DetectorSurfaceCatalogError::UnsupportedDetector);
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::MFT, 1, 10), mftSpec, lookup),
                    DetectorSurfaceCatalogError::InvalidRequest);
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::MFT, 0, 9), mftSpec, lookup),
                    DetectorSurfaceCatalogError::InvalidRequest);
  BOOST_CHECK(!lookedUp);
}

BOOST_AUTO_TEST_CASE(geometry_initialization_and_cache_failures)
{
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, 7), itsSpec,
                                                []() -> std::optional<DetectorGeometryView> { return std::nullopt; }),
                    DetectorSurfaceCatalogError::GeometryNotInitialized);

  auto absent = geometry(7, 7);
  absent.l2gCacheFilled = false;
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, 7), itsSpec,
                                                [absent] { return absent; }),
                    DetectorSurfaceCatalogError::GeometryUnavailable);

  auto undersized = geometry(7, 7);
  undersized.l2gCacheSize = 6;
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, 7), itsSpec,
                                                [undersized] { return undersized; }),
                    DetectorSurfaceCatalogError::GeometryUnavailable);
}

BOOST_AUTO_TEST_CASE(surface_lookup_failures)
{
  auto outOfRange = geometry(7, 7);
  outOfRange.surfaceForChip = [](size_t) { return 7; };
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, 7), itsSpec,
                                                [outOfRange] { return outOfRange; }),
                    DetectorSurfaceCatalogError::SurfaceLookupFailure);

  const auto emptyBucket = geometry(6, 7);
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, 7), itsSpec,
                                                [emptyBucket] { return emptyBucket; }),
                    DetectorSurfaceCatalogError::SurfaceLookupFailure);
}

BOOST_AUTO_TEST_CASE(invalid_surface_geometry_is_mapped_and_transactional)
{
  auto invalid = geometry(7, 7);
  invalid.localToGlobal = [](size_t, const GeometryPoint&) {
    return GeometryPoint{std::numeric_limits<double>::quiet_NaN(), 0., 0.};
  };
  checkEmptyFailure(buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, 7), itsSpec,
                                                [invalid] { return invalid; }),
                    DetectorSurfaceCatalogError::InvalidSurfaceGeometry);
}

BOOST_AUTO_TEST_CASE(successful_its_and_mft_catalogs)
{
  const auto itsGeometry = geometry(7, 7);
  const auto its = buildGeometrySurfaceCatalog(request(o2::detectors::DetID::ITS, 0, 7), itsSpec,
                                               [itsGeometry] { return itsGeometry; });
  BOOST_REQUIRE(its.ok());
  BOOST_REQUIRE_EQUAL(its.catalog.size(), 7);
  for (uint16_t i = 0; i < its.catalog.size(); ++i) {
    BOOST_CHECK(its.catalog[i].id == SurfaceId{i});
    BOOST_CHECK_EQUAL(its.catalog[i].detectorSurfaceIndex, i);
    BOOST_CHECK_EQUAL(its.catalog[i].detectorId, static_cast<uint8_t>(o2::detectors::DetID::ITS));
    BOOST_CHECK(its.catalog[i].kind == SurfaceKind::Cylinder);
    BOOST_CHECK_EQUAL(its.catalog[i].flags, 0);
    BOOST_CHECK_LE(its.catalog[i].radialMin, its.catalog[i].radialMax);
  }

  const auto mftGeometry = geometry(10, 10);
  const auto mft = buildGeometrySurfaceCatalog(request(o2::detectors::DetID::MFT, 0, 10), mftSpec,
                                               [mftGeometry] { return mftGeometry; });
  BOOST_REQUIRE(mft.ok());
  BOOST_REQUIRE_EQUAL(mft.catalog.size(), 10);
  for (uint16_t i = 0; i < mft.catalog.size(); ++i) {
    BOOST_CHECK(mft.catalog[i].id == SurfaceId{i});
    BOOST_CHECK_EQUAL(mft.catalog[i].detectorSurfaceIndex, i);
    BOOST_CHECK_EQUAL(mft.catalog[i].detectorId, static_cast<uint8_t>(o2::detectors::DetID::MFT));
    BOOST_CHECK(mft.catalog[i].kind == SurfaceKind::Disk);
    BOOST_CHECK_EQUAL(mft.catalog[i].flags, 0);
    BOOST_CHECK_LT(mft.catalog[i].referenceCoordinate, 0.f);
  }
}

BOOST_AUTO_TEST_CASE(public_providers_do_not_create_missing_singletons_for_bad_requests)
{
  ITSSurfaceCatalogProvider its;
  checkEmptyFailure(its.buildCatalog(request(o2::detectors::DetID::MFT, 0, 7)),
                    DetectorSurfaceCatalogError::UnsupportedDetector);
  checkEmptyFailure(its.buildCatalog(request(o2::detectors::DetID::ITS, 1, 7)),
                    DetectorSurfaceCatalogError::InvalidRequest);
  checkEmptyFailure(its.buildCatalog(request(o2::detectors::DetID::ITS, 0, 7)),
                    DetectorSurfaceCatalogError::GeometryNotInitialized);

  MFTSurfaceCatalogProvider mft;
  checkEmptyFailure(mft.buildCatalog(request(o2::detectors::DetID::ITS, 0, 10)),
                    DetectorSurfaceCatalogError::UnsupportedDetector);
  checkEmptyFailure(mft.buildCatalog(request(o2::detectors::DetID::MFT, 0, 9)),
                    DetectorSurfaceCatalogError::InvalidRequest);
  checkEmptyFailure(mft.buildCatalog(request(o2::detectors::DetID::MFT, 0, 10)),
                    DetectorSurfaceCatalogError::GeometryNotInitialized);
}
