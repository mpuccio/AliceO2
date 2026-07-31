// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B1 Slice 1: standalone offline nominal-geometry validator.
//
// Outside O2::ITSMFTTracking's and every workflow's link graph: this
// executable links only O2::DetectorsBase (GeometryManager), O2::ITSBase /
// O2::MFTBase (GeometryTGeo), O2::ITSMFTBase (SegmentationAlpide), and
// O2::CommonUtils (NameConf); it never links O2::ITSMFTTracking and never
// runs from any DPL workflow.
//
// Loads geometry with GeometryManager::loadGeometry(prefix,
// /*applyMisalignment=*/false, /*preferAlignedFile=*/false) -- the
// unaligned-file, no-misalignment-applied path (GeometryManager.cxx),
// which never touches CCDB or o2::base::Aligner. Slice 1 reports the
// geometry-derived per-surface values only; it does not yet compare against
// any authored ITSSurfaceSpec/MFTSurfaceSpec table (none exists yet -- Gate
// 4 B1 report).
//
// Usage: o2-itsmft-nominal-geometry-validator --geometry <prefix-or-path>
//        --detector <ITS|MFT> --surfaces <N> [--format text|json]

#include <cstdlib>
#include <filesystem>
#include <iostream>

#include "CommonUtils/NameConf.h"
#include "DetectorsBase/GeometryManager.h"
#include "GeometryTGeoSurfaceSource.h"
#include "ITSBase/GeometryTGeo.h"
#include "MFTBase/GeometryTGeo.h"
#include "MathUtils/Cartesian.h"
#include "NominalGeometryReport.h"
#include "ValidatorArgs.h"

namespace
{
using namespace o2::itsmft::validation;

int exitCodeFor(ValidatorStatus status)
{
  switch (status) {
    case ValidatorStatus::Ok:
      return 0;
    case ValidatorStatus::InvalidArguments:
      return 2;
    case ValidatorStatus::GeometryLoadFailed:
      return 3;
    case ValidatorStatus::GeometryNotInitialized:
      return 4;
    case ValidatorStatus::AggregationFailed:
      return 5;
  }
  return 1;
}

void print(const ValidationReport& report, OutputFormat format)
{
  std::cout << (format == OutputFormat::Json ? formatMachineReadable(report) : formatHumanReadable(report)) << '\n';
}
} // namespace

int main(int argc, char** argv)
{
  const auto parsed = parseValidatorArgs(argc, argv);
  if (!parsed.ok()) {
    std::cerr << "argument error: " << parsed.diagnostic << " (" << toString(parsed.status) << ")\n"
              << "usage: " << (argc > 0 ? argv[0] : "o2-itsmft-nominal-geometry-validator")
              << " --geometry <prefix-or-path> --detector <ITS|MFT> --surfaces <N> [--format text|json]\n";
    GeometryProvenance provenance;
    print(buildFailedReport(provenance, ValidatorStatus::InvalidArguments), OutputFormat::Text);
    return exitCodeFor(ValidatorStatus::InvalidArguments);
  }

  const auto& args = parsed.args;
  GeometryProvenance provenance;
  provenance.detectorLabel = args.detectorLabel;
  provenance.preferAlignedFile = false;
  provenance.misalignmentApplied = false;

  // Preflight the file's existence ourselves: GeometryManager::loadGeometry()
  // fatals (via FairLogger) on a missing/unreadable/malformed geometry
  // object rather than returning an error, which this tool cannot catch as a
  // typed ValidatorStatus. Checking existence first turns the common "wrong
  // path" mistake into a clean, typed GeometryLoadFailed report instead of a
  // framework-level abort; a file that exists but is not a valid geometry
  // object still goes through that uncatchable fatal path -- documented, not
  // solved, here.
  provenance.geometryFilePath = o2::base::NameConf::getGeomFileName(args.geometryPrefixOrPath);
  if (!std::filesystem::exists(provenance.geometryFilePath)) {
    std::cerr << "geometry file not found: " << provenance.geometryFilePath << '\n';
    print(buildFailedReport(provenance, ValidatorStatus::GeometryLoadFailed), args.format);
    return exitCodeFor(ValidatorStatus::GeometryLoadFailed);
  }

  o2::base::GeometryManager::loadGeometry(args.geometryPrefixOrPath, /*applyMisalignment=*/false,
                                          /*preferAlignedFile=*/false);

  const auto l2gMask = o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G);
  if (args.detector == DetectorSelection::ITS) {
    o2::its::GeometryTGeo::Instance()->fillMatrixCache(l2gMask);
  } else {
    o2::mft::GeometryTGeo::Instance()->fillMatrixCache(l2gMask);
  }

  const auto sourceResult = resolveGeometrySurfaceSource(args.detector);
  if (!sourceResult.ok()) {
    std::cerr << "geometry source unavailable after load\n";
    print(buildFailedReport(provenance, ValidatorStatus::GeometryNotInitialized), args.format);
    return exitCodeFor(ValidatorStatus::GeometryNotInitialized);
  }

  const auto report = buildValidationReport(provenance, sourceResult.source.chipCount, args.surfaceCount,
                                            alpideActiveArea(), sourceResult.source.referenceCoordinate,
                                            sourceResult.source.surfaceForChip, sourceResult.source.localToGlobal);
  print(report, args.format);
  return exitCodeFor(report.status);
}
