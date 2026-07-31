// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "NominalGeometryReport.h"

#include <iomanip>
#include <sstream>

namespace o2::itsmft::validation
{

using o2::itsmft::tracking::detail::aggregateSurfaceGeometry;
using o2::itsmft::tracking::detail::SurfaceGeometryAggregationError;

const char* toString(ValidatorStatus status) noexcept
{
  switch (status) {
    case ValidatorStatus::Ok:
      return "Ok";
    case ValidatorStatus::InvalidArguments:
      return "InvalidArguments";
    case ValidatorStatus::GeometryLoadFailed:
      return "GeometryLoadFailed";
    case ValidatorStatus::GeometryNotInitialized:
      return "GeometryNotInitialized";
    case ValidatorStatus::AggregationFailed:
      return "AggregationFailed";
  }
  return "Unknown";
}

const char* toString(SurfaceGeometryAggregationError error) noexcept
{
  switch (error) {
    case SurfaceGeometryAggregationError::None:
      return "None";
    case SurfaceGeometryAggregationError::SurfaceOutOfRange:
      return "SurfaceOutOfRange";
    case SurfaceGeometryAggregationError::EmptySurface:
      return "EmptySurface";
    case SurfaceGeometryAggregationError::NonFiniteCoordinate:
      return "NonFiniteCoordinate";
    case SurfaceGeometryAggregationError::InvalidDerivedGeometry:
      return "InvalidDerivedGeometry";
  }
  return "Unknown";
}

ValidationReport buildFailedReport(GeometryProvenance provenance, ValidatorStatus status)
{
  ValidationReport report;
  report.provenance = std::move(provenance);
  report.status = status;
  return report;
}

ValidationReport buildValidationReport(
  GeometryProvenance provenance,
  size_t chipCount,
  size_t surfaceCount,
  const o2::itsmft::tracking::detail::LocalActiveArea& activeArea,
  o2::itsmft::tracking::detail::SurfaceReferenceCoordinate referenceCoordinate,
  const o2::itsmft::tracking::detail::SurfaceForChip& surfaceForChip,
  const o2::itsmft::tracking::detail::LocalToGlobal& localToGlobal)
{
  ValidationReport report;
  report.provenance = std::move(provenance);

  const auto aggregation = aggregateSurfaceGeometry(chipCount, surfaceCount, activeArea, referenceCoordinate,
                                                     surfaceForChip, localToGlobal);
  if (!aggregation.ok()) {
    report.status = ValidatorStatus::AggregationFailed;
    report.aggregationError = aggregation.error;
    return report;
  }

  report.status = ValidatorStatus::Ok;
  report.surfaces.reserve(aggregation.surfaces.size());
  for (size_t i = 0; i < aggregation.surfaces.size(); ++i) {
    report.surfaces.push_back(SurfaceDiagnostic{static_cast<uint16_t>(i), aggregation.surfaces[i]});
  }
  return report;
}

namespace
{
// Fixed precision, classic ("C") locale, no grouping: the same
// ValidationReport always formats to the same string, independent of the
// environment's locale/thousands-separator configuration.
std::ostringstream freshStream()
{
  std::ostringstream stream;
  stream.imbue(std::locale::classic());
  stream << std::fixed << std::setprecision(6);
  return stream;
}

std::string escapeJson(const std::string& value)
{
  std::string escaped;
  escaped.reserve(value.size());
  for (char c : value) {
    if (c == '"' || c == '\\') {
      escaped.push_back('\\');
    }
    escaped.push_back(c);
  }
  return escaped;
}
} // namespace

std::string formatHumanReadable(const ValidationReport& report)
{
  auto stream = freshStream();
  stream << "NOMINAL_GEOMETRY_VALIDATION_REPORT\n"
         << "geometry_file=" << report.provenance.geometryFilePath << '\n'
         << "detector=" << report.provenance.detectorLabel << '\n'
         << "prefer_aligned_file=" << (report.provenance.preferAlignedFile ? "true" : "false") << '\n'
         << "misalignment_applied=" << (report.provenance.misalignmentApplied ? "true" : "false") << '\n'
         << "status=" << toString(report.status) << '\n';
  if (report.status == ValidatorStatus::AggregationFailed) {
    stream << "aggregation_error=" << toString(report.aggregationError) << '\n';
  }
  stream << "surface_count=" << report.surfaces.size() << '\n';
  for (const auto& surface : report.surfaces) {
    stream << "surface[" << surface.surfaceIndex << "] reference=" << surface.aggregated.referenceCoordinate
           << " radialMin=" << surface.aggregated.radialMin << " radialMax=" << surface.aggregated.radialMax
           << '\n';
  }
  return stream.str();
}

std::string formatMachineReadable(const ValidationReport& report)
{
  auto stream = freshStream();
  stream << "{"
         << "\"geometryFile\":\"" << escapeJson(report.provenance.geometryFilePath) << "\","
         << "\"detector\":\"" << escapeJson(report.provenance.detectorLabel) << "\","
         << "\"preferAlignedFile\":" << (report.provenance.preferAlignedFile ? "true" : "false") << ","
         << "\"misalignmentApplied\":" << (report.provenance.misalignmentApplied ? "true" : "false") << ","
         << "\"status\":\"" << toString(report.status) << "\",";
  if (report.status == ValidatorStatus::AggregationFailed) {
    stream << "\"aggregationError\":\"" << toString(report.aggregationError) << "\",";
  }
  stream << "\"surfaces\":[";
  for (size_t i = 0; i < report.surfaces.size(); ++i) {
    if (i != 0) {
      stream << ",";
    }
    const auto& surface = report.surfaces[i];
    stream << "{\"index\":" << surface.surfaceIndex << ",\"referenceCoordinate\":" << surface.aggregated.referenceCoordinate
           << ",\"radialMin\":" << surface.aggregated.radialMin << ",\"radialMax\":" << surface.aggregated.radialMax << "}";
  }
  stream << "]}";
  return stream.str();
}

} // namespace o2::itsmft::validation
