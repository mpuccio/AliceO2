// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B1 Slice 1: argument-validation tests for the detached
// nominal-geometry validator's CLI. Pure argc/argv logic -- no geometry
// file, no GeometryTGeo, no synthetic seam even needed here.

#define BOOST_TEST_MODULE Test ITSMFTTrackingValidation ValidatorArgs
#include <boost/test/unit_test.hpp>

#include <vector>

#include "ValidatorArgs.h"

using namespace o2::itsmft::validation;

namespace
{
ArgParseResult parse(std::vector<const char*> args)
{
  args.insert(args.begin(), "o2-itsmft-nominal-geometry-validator");
  return parseValidatorArgs(static_cast<int>(args.size()), args.data());
}
} // namespace

BOOST_AUTO_TEST_CASE(ValidFullParse)
{
  const auto result = parse({"--geometry", "o2sim", "--detector", "ITS", "--surfaces", "7"});
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result.args.geometryPrefixOrPath, "o2sim");
  BOOST_CHECK(result.args.detector == DetectorSelection::ITS);
  BOOST_CHECK_EQUAL(result.args.detectorLabel, "ITS");
  BOOST_CHECK_EQUAL(result.args.surfaceCount, 7u);
  BOOST_CHECK(result.args.format == OutputFormat::Text);
}

BOOST_AUTO_TEST_CASE(DetectorLabelIsCaseInsensitiveButPreservedVerbatim)
{
  const auto result = parse({"--geometry", "o2sim", "--detector", "mft", "--surfaces", "10"});
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(result.args.detector == DetectorSelection::MFT);
  BOOST_CHECK_EQUAL(result.args.detectorLabel, "mft");
}

BOOST_AUTO_TEST_CASE(JsonFormatParses)
{
  const auto result = parse({"--geometry", "o2sim", "--detector", "ITS", "--surfaces", "7", "--format", "json"});
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(result.args.format == OutputFormat::Json);
}

BOOST_AUTO_TEST_CASE(MissingGeometryIsRejected)
{
  const auto result = parse({"--detector", "ITS", "--surfaces", "7"});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.status == ArgParseStatus::MissingRequiredArgument);
}

BOOST_AUTO_TEST_CASE(MissingDetectorIsRejected)
{
  const auto result = parse({"--geometry", "o2sim", "--surfaces", "7"});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.status == ArgParseStatus::MissingRequiredArgument);
}

BOOST_AUTO_TEST_CASE(MissingSurfacesIsRejected)
{
  const auto result = parse({"--geometry", "o2sim", "--detector", "ITS"});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.status == ArgParseStatus::MissingRequiredArgument);
}

BOOST_AUTO_TEST_CASE(UnknownDetectorIsRejected)
{
  const auto result = parse({"--geometry", "o2sim", "--detector", "TPC", "--surfaces", "7"});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.status == ArgParseStatus::UnknownDetector);
}

BOOST_AUTO_TEST_CASE(ZeroSurfacesIsRejected)
{
  const auto result = parse({"--geometry", "o2sim", "--detector", "ITS", "--surfaces", "0"});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.status == ArgParseStatus::InvalidValue);
}

BOOST_AUTO_TEST_CASE(NonNumericSurfacesIsRejected)
{
  const auto result = parse({"--geometry", "o2sim", "--detector", "ITS", "--surfaces", "seven"});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.status == ArgParseStatus::InvalidValue);
}

BOOST_AUTO_TEST_CASE(NegativeSurfacesIsRejected)
{
  const auto result = parse({"--geometry", "o2sim", "--detector", "ITS", "--surfaces", "-1"});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.status == ArgParseStatus::InvalidValue);
}

BOOST_AUTO_TEST_CASE(UnknownFormatIsRejected)
{
  const auto result = parse({"--geometry", "o2sim", "--detector", "ITS", "--surfaces", "7", "--format", "xml"});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.status == ArgParseStatus::UnknownFormat);
}

BOOST_AUTO_TEST_CASE(UnknownOptionIsRejected)
{
  const auto result = parse({"--geometry", "o2sim", "--detector", "ITS", "--surfaces", "7", "--bogus", "value"});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.status == ArgParseStatus::UnknownOption);
}

BOOST_AUTO_TEST_CASE(DanglingValuelessOptionIsRejected)
{
  const auto result = parse({"--geometry", "o2sim", "--detector", "ITS", "--surfaces"});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.status == ArgParseStatus::MissingRequiredArgument);
}

BOOST_AUTO_TEST_CASE(NoArgumentsAtAllIsRejected)
{
  const auto result = parse({});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.status == ArgParseStatus::MissingRequiredArgument);
}
