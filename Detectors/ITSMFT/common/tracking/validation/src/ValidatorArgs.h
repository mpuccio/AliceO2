// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B1 Slice 1: CLI argument parsing/validation for the detached
// nominal-geometry validator. Pure argc/argv logic -- never touches
// GeometryTGeo, GeometryManager, or any file I/O -- so
// testValidatorArgs.cxx exercises it without any geometry file or
// synthetic-seam plumbing at all.
//
// Surface count is a required argument, not a compiled-in ITS=7/MFT=10
// assumption: Gate 4 B1's accepted design forbids inventing detector
// constants in this slice, and topology surface counts are exactly the kind
// of fact the future ITSSurfaceSpec/MFTSurfaceSpec tables (not this tool)
// are meant to become the authority on.

#ifndef ALICEO2_ITSMFT_VALIDATION_VALIDATORARGS_H_
#define ALICEO2_ITSMFT_VALIDATION_VALIDATORARGS_H_

#include <cstddef>
#include <cstdint>
#include <string>

namespace o2::itsmft::validation
{

enum class DetectorSelection : uint8_t {
  ITS,
  MFT
};

enum class OutputFormat : uint8_t {
  Text,
  Json
};

struct ValidatorArgs {
  std::string geometryPrefixOrPath{};
  DetectorSelection detector{DetectorSelection::ITS};
  std::string detectorLabel{}; // as given on the command line, preserved verbatim for report provenance
  size_t surfaceCount{0};
  OutputFormat format{OutputFormat::Text};
};

enum class ArgParseStatus : uint8_t {
  Ok,
  MissingRequiredArgument,
  InvalidValue,
  UnknownDetector,
  UnknownOption,
  UnknownFormat
};

struct ArgParseResult {
  ValidatorArgs args{};
  ArgParseStatus status{ArgParseStatus::Ok};
  std::string diagnostic{};

  bool ok() const noexcept { return status == ArgParseStatus::Ok; }
};

// argv[0] (the program name) is expected at index 0, exactly as passed to
// main(); parsing starts at argv[1].
ArgParseResult parseValidatorArgs(int argc, const char* const* argv);

const char* toString(ArgParseStatus status) noexcept;

} // namespace o2::itsmft::validation

#endif /* ALICEO2_ITSMFT_VALIDATION_VALIDATORARGS_H_ */
