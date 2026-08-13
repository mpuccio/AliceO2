// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Detached nominal-geometry validator infrastructure. It does not depend on
// GeometryTGeo or CCDB: callers provide the geometry callbacks explicitly,
// allowing deterministic tests with synthetic chip layouts. The report
// contains geometry-derived surface values and their provenance; comparison
// with authored static tables is outside this slice.

#ifndef ALICEO2_ITSMFT_VALIDATION_NOMINALGEOMETRYREPORT_H_
#define ALICEO2_ITSMFT_VALIDATION_NOMINALGEOMETRYREPORT_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "DetectorSurfaceCatalogAggregation.h"

namespace o2::itsmft::validation
{

// Records the geometry file and loading options supplied by the caller.
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

// Aggregates caller-supplied geometry callbacks into a report. Reports
// AggregationFailed with the underlying error when aggregation fails.
ValidationReport buildValidationReport(
  GeometryProvenance provenance,
  size_t chipCount,
  size_t surfaceCount,
  const o2::itsmft::tracking::detail::LocalActiveArea& activeArea,
  o2::itsmft::tracking::detail::SurfaceReferenceCoordinate referenceCoordinate,
  const o2::itsmft::tracking::detail::SurfaceForChip& surfaceForChip,
  const o2::itsmft::tracking::detail::LocalToGlobal& localToGlobal);

// Builds a report for a failure that occurred before aggregation.
ValidationReport buildFailedReport(GeometryProvenance provenance, ValidatorStatus status);

// Returns the shortest std::to_chars representation that round-trips through
// std::from_chars and is valid as a JSON number.
std::string formatLosslessFloat(float value);

// Both formatters are deterministic for identical input. Human-readable
// output uses six decimal places; machine-readable geometry floats use the
// shortest round-tripping representation.
std::string formatHumanReadable(const ValidationReport& report);
std::string formatMachineReadable(const ValidationReport& report);

} // namespace o2::itsmft::validation

#endif /* ALICEO2_ITSMFT_VALIDATION_NOMINALGEOMETRYREPORT_H_ */
