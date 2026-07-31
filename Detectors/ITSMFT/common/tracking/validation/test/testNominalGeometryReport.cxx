// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B1 Slice 1: report-building/formatting-determinism tests for the
// detached nominal-geometry validator, using synthetic surfaceForChip/
// localToGlobal callbacks -- never GeometryTGeo, never a geometry file. The
// production entry point (nominal-geometry-validator.cxx) supplies
// GeometryTGeo-backed callbacks instead of these synthetic ones;
// buildValidationReport() itself cannot tell the difference.

#define BOOST_TEST_MODULE Test ITSMFTTrackingValidation NominalGeometryReport
#include <boost/test/unit_test.hpp>

#include <limits>
#include <string>

#include "NominalGeometryReport.h"

using namespace o2::itsmft::validation;
using o2::itsmft::tracking::detail::GeometryPoint;
using o2::itsmft::tracking::detail::LocalActiveArea;
using o2::itsmft::tracking::detail::SurfaceGeometryAggregationError;
using o2::itsmft::tracking::detail::SurfaceReferenceCoordinate;

namespace
{
// Two chips on one surface, unit-square active area, chip N translated by
// 10*N along x with no rotation: chip centers land at global x=0 and x=10,
// so the surface's mean radius (MeanRadius reference, hypot(x,y) of each
// chip center) is exactly (0+10)/2 = 5, hand-checkable without relying on
// GeometryTGeo or any detector constant.
GeometryProvenance syntheticProvenance()
{
  GeometryProvenance provenance;
  provenance.geometryFilePath = "synthetic://two-chip-one-surface";
  provenance.detectorLabel = "SYN";
  return provenance;
}

LocalActiveArea unitSquare()
{
  return LocalActiveArea{-0.5, 0.5, -0.5, 0.5};
}

ValidationReport buildTwoChipReport()
{
  return buildValidationReport(
    syntheticProvenance(), /*chipCount=*/2, /*surfaceCount=*/1, unitSquare(), SurfaceReferenceCoordinate::MeanRadius,
    [](size_t /*chip*/) { return 0; },
    [](size_t chip, const GeometryPoint& local) {
      return GeometryPoint{local.x + 10.0 * static_cast<double>(chip), local.y, local.z};
    });
}
} // namespace

BOOST_AUTO_TEST_CASE(SyntheticTwoChipReportIsOkAndHandCheckable)
{
  const auto report = buildTwoChipReport();
  BOOST_REQUIRE(report.ok());
  BOOST_REQUIRE_EQUAL(report.surfaces.size(), 1u);
  BOOST_CHECK_CLOSE(report.surfaces[0].aggregated.referenceCoordinate, 5.f, 1e-4f);
  BOOST_CHECK_EQUAL(report.surfaces[0].surfaceIndex, 0u);
  BOOST_CHECK(report.status == ValidatorStatus::Ok);
}

BOOST_AUTO_TEST_CASE(ReportBuildingIsDeterministic)
{
  const auto first = buildTwoChipReport();
  const auto second = buildTwoChipReport();
  BOOST_REQUIRE_EQUAL(first.surfaces.size(), second.surfaces.size());
  BOOST_CHECK_EQUAL(first.surfaces[0].aggregated.referenceCoordinate, second.surfaces[0].aggregated.referenceCoordinate);
  BOOST_CHECK_EQUAL(first.surfaces[0].aggregated.radialMin, second.surfaces[0].aggregated.radialMin);
  BOOST_CHECK_EQUAL(first.surfaces[0].aggregated.radialMax, second.surfaces[0].aggregated.radialMax);
}

BOOST_AUTO_TEST_CASE(HumanReadableFormattingIsDeterministic)
{
  const auto first = formatHumanReadable(buildTwoChipReport());
  const auto second = formatHumanReadable(buildTwoChipReport());
  BOOST_CHECK_EQUAL(first, second);
  BOOST_CHECK(first.find("status=Ok") != std::string::npos);
  BOOST_CHECK(first.find("surface_count=1") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(MachineReadableFormattingIsDeterministic)
{
  const auto first = formatMachineReadable(buildTwoChipReport());
  const auto second = formatMachineReadable(buildTwoChipReport());
  BOOST_CHECK_EQUAL(first, second);
  BOOST_CHECK(first.find("\"status\":\"Ok\"") != std::string::npos);
  BOOST_CHECK(first.find("\"surfaces\":[") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(SurfaceOutOfRangeIsReportedAsAggregationFailed)
{
  const auto report = buildValidationReport(
    syntheticProvenance(), /*chipCount=*/1, /*surfaceCount=*/1, unitSquare(), SurfaceReferenceCoordinate::MeanRadius,
    [](size_t /*chip*/) { return 5; }, // out of range: only surface 0 exists
    [](size_t /*chip*/, const GeometryPoint& local) { return local; });

  BOOST_REQUIRE(!report.ok());
  BOOST_CHECK(report.status == ValidatorStatus::AggregationFailed);
  BOOST_CHECK(report.aggregationError == SurfaceGeometryAggregationError::SurfaceOutOfRange);
  BOOST_CHECK(report.surfaces.empty());
}

BOOST_AUTO_TEST_CASE(NonFiniteCoordinateIsReportedAsAggregationFailed)
{
  const auto report = buildValidationReport(
    syntheticProvenance(), /*chipCount=*/1, /*surfaceCount=*/1, unitSquare(), SurfaceReferenceCoordinate::MeanRadius,
    [](size_t /*chip*/) { return 0; },
    [](size_t /*chip*/, const GeometryPoint& local) {
      return GeometryPoint{std::numeric_limits<double>::quiet_NaN(), local.y, local.z};
    });

  BOOST_REQUIRE(!report.ok());
  BOOST_CHECK(report.status == ValidatorStatus::AggregationFailed);
  BOOST_CHECK(report.aggregationError == SurfaceGeometryAggregationError::NonFiniteCoordinate);
}

BOOST_AUTO_TEST_CASE(FailedReportsBeforeAggregationCarryProvenanceAndEmptySurfaces)
{
  for (auto status : {ValidatorStatus::InvalidArguments, ValidatorStatus::GeometryLoadFailed,
                      ValidatorStatus::GeometryNotInitialized}) {
    const auto report = buildFailedReport(syntheticProvenance(), status);
    BOOST_CHECK(!report.ok());
    BOOST_CHECK(report.status == status);
    BOOST_CHECK(report.surfaces.empty());
    BOOST_CHECK(report.aggregationError == SurfaceGeometryAggregationError::None);
    BOOST_CHECK_EQUAL(report.provenance.geometryFilePath, "synthetic://two-chip-one-surface");
  }
}

BOOST_AUTO_TEST_CASE(EveryStatusAndErrorHasANonUnknownName)
{
  for (auto status : {ValidatorStatus::Ok, ValidatorStatus::InvalidArguments, ValidatorStatus::GeometryLoadFailed,
                      ValidatorStatus::GeometryNotInitialized, ValidatorStatus::AggregationFailed}) {
    BOOST_CHECK(std::string{toString(status)} != "Unknown");
  }
  for (auto error : {SurfaceGeometryAggregationError::None, SurfaceGeometryAggregationError::SurfaceOutOfRange,
                     SurfaceGeometryAggregationError::EmptySurface, SurfaceGeometryAggregationError::NonFiniteCoordinate,
                     SurfaceGeometryAggregationError::InvalidDerivedGeometry}) {
    BOOST_CHECK(std::string{toString(error)} != "Unknown");
  }
}
