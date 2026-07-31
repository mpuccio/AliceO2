// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B1 Slice 1: detached nominal-geometry validator infrastructure.
//
// This header is deliberately independent of GeometryTGeo, CCDB, and
// O2::ITSMFTTracking: buildValidationReport() below is a pure function of
// its explicit chip/surface/callback arguments, exactly like
// aggregateSurfaceGeometry() itself (DetectorSurfaceCatalogAggregation.h,
// compiled directly into this validator target -- see this directory's
// CMakeLists.txt; the real GeometryTGeo-backed callbacks are supplied only
// by GeometryTGeoSurfaceSource.{h,cxx} and nominal-geometry-validator.cxx,
// never by this file). This split is what lets testNominalGeometryReport.cxx
// exercise report building and formatting deterministically, with synthetic
// chip layouts, and no geometry file at all.
//
// Slice 1 produces no ITSSurfaceSpec/MFTSurfaceSpec comparison: there is
// nothing authored yet to compare against (Gate 4 B1 report). This tool
// reports the geometry-derived per-surface values and their provenance; a
// later slice adds comparison against the authored static tables once they
// exist and are signed off.

#ifndef ALICEO2_ITSMFT_VALIDATION_NOMINALGEOMETRYREPORT_H_
#define ALICEO2_ITSMFT_VALIDATION_NOMINALGEOMETRYREPORT_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "DetectorSurfaceCatalogAggregation.h"

namespace o2::itsmft::validation
{

// Provenance is recorded, never inferred: the caller (nominal-geometry-
// validator.cxx's main()) states exactly which file it asked
// GeometryManager::loadGeometry() to load and with what arguments, so the
// report is self-describing regardless of how it is later archived or
// forwarded for detector-owner sign-off.
struct GeometryProvenance {
  std::string geometryFilePath{};
  std::string detectorLabel{};
  bool preferAlignedFile{false};
  bool misalignmentApplied{false};
};

enum class ValidatorStatus : uint8_t {
  Ok,
  InvalidArguments,
  GeometryLoadFailed,
  GeometryNotInitialized,
  AggregationFailed
};

const char* toString(ValidatorStatus status) noexcept;
const char* toString(o2::itsmft::tracking::detail::SurfaceGeometryAggregationError error) noexcept;

struct SurfaceDiagnostic {
  uint16_t surfaceIndex{0};
  o2::itsmft::tracking::detail::AggregatedSurfaceGeometry aggregated{};
};

struct ValidationReport {
  GeometryProvenance provenance{};
  ValidatorStatus status{ValidatorStatus::Ok};
  o2::itsmft::tracking::detail::SurfaceGeometryAggregationError aggregationError{
    o2::itsmft::tracking::detail::SurfaceGeometryAggregationError::None};
  std::vector<SurfaceDiagnostic> surfaces{};

  bool ok() const noexcept { return status == ValidatorStatus::Ok; }
};

// Pure aggregation-and-report step. Never touches GeometryTGeo, CCDB, or any
// alignment object itself -- surfaceForChip/localToGlobal are supplied by
// the caller, synthetic in tests, GeometryTGeo-backed in production. Reports
// AggregationFailed (with the specific aggregationError populated) on any
// aggregateSurfaceGeometry() failure; never throws.
ValidationReport buildValidationReport(
  GeometryProvenance provenance,
  size_t chipCount,
  size_t surfaceCount,
  const o2::itsmft::tracking::detail::LocalActiveArea& activeArea,
  o2::itsmft::tracking::detail::SurfaceReferenceCoordinate referenceCoordinate,
  const o2::itsmft::tracking::detail::SurfaceForChip& surfaceForChip,
  const o2::itsmft::tracking::detail::LocalToGlobal& localToGlobal);

// A report built from an InvalidArguments/GeometryLoadFailed/
// GeometryNotInitialized status before any aggregation attempt (surfaces
// stays empty; aggregationError stays None).
ValidationReport buildFailedReport(GeometryProvenance provenance, ValidatorStatus status);

// The exact JSON-number token formatMachineReadable() embeds for each
// geometry float field (Gate 4 acceptance-cleanup C1): std::to_chars'
// shortest round-tripping representation, always valid JSON number syntax.
// Exposed on its own -- not just exercised indirectly through a full
// report -- so its round-trip guarantee can be proven directly against
// std::from_chars (the exact inverse it is built on) without going through
// a third-party JSON parser's own numeric-type classification, which is a
// separate concern this function makes no claim about.
std::string formatLosslessFloat(float value);

// Both deterministic given identical input: no timestamps, no
// locale-dependent formatting, no pointer/address-derived content.
//
// formatHumanReadable(): float values use a fixed, explicit six-decimal
// precision -- readable, not intended as a numeric-provenance source.
//
// formatMachineReadable(): geometry float fields (referenceCoordinate,
// radialMin, radialMax) use std::to_chars' shortest round-tripping
// representation instead (Gate 4 acceptance-cleanup C1) -- the printed
// token, re-parsed with std::from_chars, reproduces the exact source
// float32 bit pattern. Still deterministic: to_chars' shortest-round-trip
// output is itself a pure function of the input value.
std::string formatHumanReadable(const ValidationReport& report);
std::string formatMachineReadable(const ValidationReport& report);

} // namespace o2::itsmft::validation

#endif /* ALICEO2_ITSMFT_VALIDATION_NOMINALGEOMETRYREPORT_H_ */
