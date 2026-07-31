// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 acceptance-cleanup C1: proves formatLosslessFloat() -- the
// encoder formatMachineReadable() uses for referenceCoordinate, radialMin,
// and radialMax -- round-trips exactly: original float32 bit pattern ->
// printed token -> re-parsed float32 bit pattern, for known nontrivial
// positive and negative bit patterns, including cases where equality alone
// would be misleading (+0.0f vs -0.0f compare equal as floats but have
// different bit patterns).
//
// Decoding uses std::strtof, not std::from_chars, for two independent
// reasons found while building this test, both worth recording:
//
// 1. std::from_chars(float) is unavailable at this toolchain's deployment
//    target: Apple's libc++ marks the floating-point overload
//    availability(macos,introduced=26.0), and this build compiles with
//    -mmacosx-version-min=14.0 (confirmed by compiling a minimal probe
//    with the production compile line's exact flags -- see the commit
//    message). std::to_chars(float) carries no such restriction and is
//    what NominalGeometryReport.cxx's encoder actually uses; only the
//    decode side needed an alternative here.
// 2. RapidJSON's own GetFloat() is not that alternative either: parsing a
//    full report payload and reading surfaces[i].referenceCoordinate back
//    through RapidJSON was tried first and found to silently drop -0's
//    sign bit (a token with no '.' or 'e' is integer-classified internally
//    and widened through a signless path) and to classify some
//    legitimately-shortest float tokens (e.g. "2160324864" -- the shortest
//    round-tripping form for a float value whose magnitude exceeds
//    float32's integer-exactness range) as an integer type rather than a
//    float type. Both are real, reproducible RapidJSON behaviors, not
//    defects in the encoder under test; they are simply not part of what
//    this encoder promises, so they should not gate this proof.
//
// std::strtof is C, locale-sensitive only in its decimal-point character,
// and carries no deployment-target restriction; LC_NUMERIC is pinned to
// "C" for this file's process for that reason, matching to_chars' own
// locale-independent '.' output.

#define BOOST_TEST_MODULE Test ITSMFTTrackingValidation LosslessFloatRoundTrip
#include <boost/test/unit_test.hpp>

#include <clocale>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ios>
#include <string>

#include <rapidjson/document.h>

#include "NominalGeometryReport.h"

using namespace o2::itsmft::validation;

namespace
{
struct PinCLocale {
  PinCLocale() { std::setlocale(LC_NUMERIC, "C"); }
};
BOOST_GLOBAL_FIXTURE(PinCLocale);

uint32_t bitsOf(float value)
{
  uint32_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

float floatFromBits(uint32_t bits)
{
  float value{};
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

void checkRoundTrips(float value)
{
  const std::string token = formatLosslessFloat(value);
  BOOST_REQUIRE_MESSAGE(!token.empty(), "formatLosslessFloat produced an empty token");

  char* endptr = nullptr;
  const float parsed = std::strtof(token.c_str(), &endptr);
  BOOST_REQUIRE_MESSAGE(endptr == token.c_str() + token.size(),
                        "strtof left unparsed trailing characters in \"" << token << "\"");

  BOOST_CHECK_MESSAGE(bitsOf(parsed) == bitsOf(value), "bit mismatch: original=0x" << std::hex << bitsOf(value)
                                                                                    << " parsed=0x" << bitsOf(parsed)
                                                                                    << std::dec << " token=\"" << token << "\"");
}
} // namespace

BOOST_AUTO_TEST_CASE(RealGeometryValuePositive)
{
  // ITS surface[6]'s actual reference coordinate from the Slice 1
  // real-fixture run -- a genuine, previously-observed geometry value, not
  // an arbitrary literal.
  checkRoundTrips(39.310642f);
}

BOOST_AUTO_TEST_CASE(RealGeometryValueNegative)
{
  // MFT surface[9]'s actual reference coordinate from the same run.
  checkRoundTrips(-77.511101f);
}

BOOST_AUTO_TEST_CASE(PositiveAndNegativeZeroAreDistinguishedByBits)
{
  // +0.0f and -0.0f compare equal as floats but differ in their sign bit;
  // a plain BOOST_CHECK_EQUAL(parsed, value) would not catch a formatter
  // that silently normalized -0 to "0". checkRoundTrips()'s bitwise
  // comparison does.
  checkRoundTrips(0.0f);
  checkRoundTrips(-0.0f);
}

BOOST_AUTO_TEST_CASE(ValuesNeedingMaximumSignificantDigits)
{
  // Adjacent-representable-value bit patterns: 1.0f's immediate successor
  // and predecessor need every one of float's max_digits10 (9) significant
  // digits to disambiguate from 1.0f itself -- the kind of value most
  // likely to expose a formatter that silently truncates precision.
  checkRoundTrips(floatFromBits(0x3F800001u)); // just above 1.0f
  checkRoundTrips(floatFromBits(0x3F7FFFFFu)); // just below 1.0f
  checkRoundTrips(floatFromBits(0xBF800001u)); // just below -1.0f (more negative)
}

BOOST_AUTO_TEST_CASE(ExtremeMagnitudesAndSubnormals)
{
  checkRoundTrips(1.0e30f);
  checkRoundTrips(-1.0e30f);
  checkRoundTrips(1.0e-30f);
  checkRoundTrips(-1.234567e-5f);
  checkRoundTrips(floatFromBits(0x00000001u)); // smallest positive subnormal
  checkRoundTrips(floatFromBits(0x80000001u)); // smallest-magnitude negative subnormal
  checkRoundTrips(floatFromBits(0x7F7FFFFFu)); // largest finite positive
  checkRoundTrips(floatFromBits(0xFF7FFFFFu)); // largest finite negative (most negative)
}

BOOST_AUTO_TEST_CASE(SweepOfArbitraryBitPatternsRoundTrips)
{
  // A deterministic, non-cherry-picked sweep across the finite bit-pattern
  // space (skipping NaN/Inf exponent-all-ones patterns, which
  // aggregateSurfaceGeometry() never returns on a successful aggregation
  // and which formatLosslessFloat() is not specified to handle). Includes
  // magnitudes where float32's ULP exceeds 1 -- the shortest round-trip
  // token there is a bare integer with no '.' or 'e', exactly the case
  // that exposed RapidJSON's integer/float type-classification split
  // documented at the top of this file; strtof has no such split.
  for (uint32_t bits = 0x00000001u; bits < 0x7F800000u; bits += 104729u /* prime stride, ~20500 points */) {
    const float value = floatFromBits(bits);
    checkRoundTrips(value);
    checkRoundTrips(-value);
  }
}

BOOST_AUTO_TEST_CASE(FullReportPayloadIsStillSyntacticallyValidJson)
{
  // Narrow, honest complement to the round-trip proof above: confirms the
  // full report (with formatLosslessFloat()'s tokens embedded) still
  // parses as one JSON document -- already the subject of
  // testJsonStdoutContract.cxx's own coverage, kept here too since it is
  // the one property this file can still correctly attribute to RapidJSON
  // after the findings in this file's header comment.
  ValidationReport report;
  report.provenance.geometryFilePath = "synthetic://lossless-float-round-trip";
  report.provenance.detectorLabel = "SYN";
  report.status = ValidatorStatus::Ok;
  report.surfaces.push_back(SurfaceDiagnostic{
    0, o2::itsmft::tracking::detail::AggregatedSurfaceGeometry{-0.0f, floatFromBits(0x00000001u), 2160324864.0f}});

  const auto payload = formatMachineReadable(report);
  rapidjson::Document document;
  document.Parse(payload.c_str());
  BOOST_CHECK_MESSAGE(!document.HasParseError(), "payload did not parse as JSON: " << payload);
}
