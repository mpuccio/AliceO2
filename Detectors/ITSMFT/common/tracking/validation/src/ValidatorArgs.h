// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// CLI argument parsing for the detached nominal-geometry validator.
// This layer only processes argc/argv and does not access geometry or files.

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
  std::string detectorLabel{}; // original command-line value for report provenance
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

// argv[0] is the program name; parsing starts at argv[1].
ArgParseResult parseValidatorArgs(int argc, const char* const* argv);

const char* toString(ArgParseStatus status) noexcept;

} // namespace o2::itsmft::validation

#endif /* ALICEO2_ITSMFT_VALIDATION_VALIDATORARGS_H_ */
