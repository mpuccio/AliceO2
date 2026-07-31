// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ValidatorArgs.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstring>

namespace o2::itsmft::validation
{

const char* toString(ArgParseStatus status) noexcept
{
  switch (status) {
    case ArgParseStatus::Ok:
      return "Ok";
    case ArgParseStatus::MissingRequiredArgument:
      return "MissingRequiredArgument";
    case ArgParseStatus::InvalidValue:
      return "InvalidValue";
    case ArgParseStatus::UnknownDetector:
      return "UnknownDetector";
    case ArgParseStatus::UnknownOption:
      return "UnknownOption";
    case ArgParseStatus::UnknownFormat:
      return "UnknownFormat";
  }
  return "Unknown";
}

namespace
{
// Carries forward whatever was already successfully parsed into `args`
// (notably --format, which may have appeared before the option that
// ultimately fails) rather than discarding it: main()'s --format json
// stdout contract must hold on a validation failure too, and it can only
// do that if a --format seen before the failing option survives into the
// returned ArgParseResult.
ArgParseResult fail(ValidatorArgs args, ArgParseStatus status, std::string diagnostic)
{
  ArgParseResult result;
  result.args = std::move(args);
  result.status = status;
  result.diagnostic = std::move(diagnostic);
  return result;
}

std::string toUpper(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return value;
}
} // namespace

ArgParseResult parseValidatorArgs(int argc, const char* const* argv)
{
  if (argv == nullptr) {
    return fail(ValidatorArgs{}, ArgParseStatus::MissingRequiredArgument, "argv is null");
  }

  ValidatorArgs args;
  bool haveGeometry = false;
  bool haveDetector = false;
  bool haveSurfaces = false;

  for (int i = 1; i < argc; ++i) {
    if (argv[i] == nullptr) {
      return fail(args, ArgParseStatus::InvalidValue, "argv contains a null entry");
    }
    const std::string option{argv[i]};
    const auto needValue = [&](const char* name) -> const char* {
      if (i + 1 >= argc || argv[i + 1] == nullptr) {
        return nullptr;
      }
      ++i;
      return argv[i];
    };

    if (option == "--geometry") {
      const char* value = needValue("--geometry");
      if (value == nullptr) {
        return fail(args, ArgParseStatus::MissingRequiredArgument, "--geometry requires a value");
      }
      args.geometryPrefixOrPath = value;
      haveGeometry = true;
    } else if (option == "--detector") {
      const char* value = needValue("--detector");
      if (value == nullptr) {
        return fail(args, ArgParseStatus::MissingRequiredArgument, "--detector requires a value");
      }
      args.detectorLabel = value;
      const auto upper = toUpper(args.detectorLabel);
      if (upper == "ITS") {
        args.detector = DetectorSelection::ITS;
      } else if (upper == "MFT") {
        args.detector = DetectorSelection::MFT;
      } else {
        return fail(args, ArgParseStatus::UnknownDetector, "unknown --detector value: " + args.detectorLabel);
      }
      haveDetector = true;
    } else if (option == "--surfaces") {
      const char* value = needValue("--surfaces");
      if (value == nullptr) {
        return fail(args, ArgParseStatus::MissingRequiredArgument, "--surfaces requires a value");
      }
      const std::string_view text{value};
      size_t parsed = 0;
      const auto convResult = std::from_chars(text.data(), text.data() + text.size(), parsed);
      if (convResult.ec != std::errc{} || convResult.ptr != text.data() + text.size() || parsed == 0) {
        return fail(args, ArgParseStatus::InvalidValue, "--surfaces must be a positive integer, got: " + std::string{text});
      }
      args.surfaceCount = parsed;
      haveSurfaces = true;
    } else if (option == "--format") {
      const char* value = needValue("--format");
      if (value == nullptr) {
        return fail(args, ArgParseStatus::MissingRequiredArgument, "--format requires a value");
      }
      const auto lower = std::string{value};
      if (lower == "text") {
        args.format = OutputFormat::Text;
      } else if (lower == "json") {
        args.format = OutputFormat::Json;
      } else {
        return fail(args, ArgParseStatus::UnknownFormat, "unknown --format value: " + lower);
      }
    } else {
      return fail(args, ArgParseStatus::UnknownOption, "unknown option: " + option);
    }
  }

  if (!haveGeometry) {
    return fail(args, ArgParseStatus::MissingRequiredArgument, "--geometry is required");
  }
  if (!haveDetector) {
    return fail(args, ArgParseStatus::MissingRequiredArgument, "--detector is required");
  }
  if (!haveSurfaces) {
    return fail(args, ArgParseStatus::MissingRequiredArgument, "--surfaces is required");
  }

  ArgParseResult result;
  result.args = std::move(args);
  result.status = ArgParseStatus::Ok;
  return result;
}

} // namespace o2::itsmft::validation
