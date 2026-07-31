// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B1 Slice 1 correction: proves --format json's stdout contract --
// exactly one JSON document, on both the success and a typed-failure path
// -- by actually capturing the process's real stdout file descriptor
// (the same mechanism nominal-geometry-validator.cxx's main() uses via
// JsonStdoutGuard) around a synthetic workload that emits incidental
// [INFO]-style noise the way GeometryManager::loadGeometry()/GeometryTGeo
// do in production, and parsing the captured bytes with RapidJSON. No
// GeometryTGeo, geometry file, or O2::ITSMFTTracking dependency.

#define BOOST_TEST_MODULE Test ITSMFTTrackingValidation JsonStdoutContract
#include <boost/test/unit_test.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

#include <rapidjson/document.h>

#include "JsonStdoutGuard.h"
#include "NominalGeometryReport.h"

using namespace o2::itsmft::validation;

namespace
{
// Redirects the *test process's own* real stdout (fd 1) to a fresh temp
// file for the guarded scope, so this test can inspect exactly what
// JsonStdoutGuard::emit() writes to "the real stdout" it captured at its
// own construction. Restoration in the destructor is unconditional (runs
// even if a BOOST_REQUIRE inside the guarded scope aborts the test case),
// which is required here: without it, a failed assertion would leave every
// later test case's own console output redirected into this temp file.
class CapturedStdout
{
 public:
  // Deliberately uses plain exceptions here, never a BOOST_* check macro:
  // Boost.Test's own "check has passed" logging (at --log_level=all and
  // above) writes through std::cout/fd 1 at the point each check
  // evaluates -- so a BOOST_REQUIRE* wrapping the dup2() call below would
  // have its own pass/fail message land *in the just-redirected capture
  // file* instead of the real console, corrupting the very capture this
  // class exists to keep clean. A thrown exception here produces no
  // console output at all until it is caught well outside any redirected
  // scope (Boost's own top-level handler, after this object's destructor
  // has already restored fd 1).
  CapturedStdout()
  {
    std::string pathTemplate = (std::filesystem::temp_directory_path() /
                                "o2-itsmft-validator-stdout-capture-XXXXXX")
                                 .string();
    std::vector<char> buffer(pathTemplate.begin(), pathTemplate.end());
    buffer.push_back('\0');
    mCapturedFd = ::mkstemp(buffer.data());
    if (mCapturedFd < 0) {
      throw std::runtime_error("mkstemp() failed for stdout-capture temp file");
    }
    mPath = buffer.data();

    mSavedStdoutFd = ::dup(STDOUT_FILENO);
    if (mSavedStdoutFd < 0) {
      throw std::runtime_error("dup(STDOUT_FILENO) failed");
    }
    if (::dup2(mCapturedFd, STDOUT_FILENO) < 0) {
      throw std::runtime_error("dup2() onto stdout failed");
    }
  }

  ~CapturedStdout()
  {
    if (mSavedStdoutFd >= 0) {
      ::dup2(mSavedStdoutFd, STDOUT_FILENO);
      ::close(mSavedStdoutFd);
    }
    if (mCapturedFd >= 0) {
      ::close(mCapturedFd);
    }
    if (!mPath.empty()) {
      std::remove(mPath.c_str());
    }
  }

  CapturedStdout(const CapturedStdout&) = delete;
  CapturedStdout& operator=(const CapturedStdout&) = delete;

  std::string readAll() const
  {
    std::ifstream stream(mPath, std::ios::binary);
    std::ostringstream buffer;
    buffer << stream.rdbuf();
    return buffer.str();
  }

 private:
  std::string mPath;
  int mCapturedFd{-1};
  int mSavedStdoutFd{-1};
};

// Exactly one line ending in '\n', with no embedded newline before it --
// i.e. one, and only one, complete text record was written.
std::string requireSingleLine(const std::string& captured)
{
  BOOST_REQUIRE(!captured.empty());
  BOOST_REQUIRE_EQUAL(captured.back(), '\n');
  const auto payload = captured.substr(0, captured.size() - 1);
  BOOST_REQUIRE(payload.find('\n') == std::string::npos);
  return payload;
}

rapidjson::Document parseOrFail(const std::string& payload)
{
  rapidjson::Document document;
  document.Parse(payload.c_str());
  BOOST_REQUIRE_MESSAGE(!document.HasParseError(),
                        "captured stdout payload did not parse as JSON: " << payload);
  return document;
}
} // namespace

BOOST_AUTO_TEST_CASE(SuccessfulSyntheticReportProducesExactlyOneJsonDocumentDespiteNoise)
{
  std::string captured;
  {
    CapturedStdout capture;
    JsonStdoutGuard guard;
    BOOST_REQUIRE(guard.ok());

    // Simulate the incidental [INFO]-style stdout writes
    // GeometryManager::loadGeometry()/GeometryTGeo::fillMatrixCache() make
    // in production, via both C++ iostreams and C stdio -- dup2 must catch
    // either, and does, since it operates below both abstractions.
    std::cout << "[INFO] pretend geometry loading noise via std::cout\n";
    std::fprintf(stdout, "[INFO] pretend geometry loading noise via stdio\n");
    std::fflush(stdout);

    GeometryProvenance provenance;
    provenance.geometryFilePath = "synthetic://json-stdout-contract";
    provenance.detectorLabel = "SYN";
    const auto report = buildValidationReport(
      provenance, /*chipCount=*/1, /*surfaceCount=*/1,
      o2::itsmft::tracking::detail::LocalActiveArea{-0.5, 0.5, -0.5, 0.5},
      o2::itsmft::tracking::detail::SurfaceReferenceCoordinate::MeanRadius,
      [](size_t /*chip*/) { return 0; },
      [](size_t /*chip*/, const o2::itsmft::tracking::detail::GeometryPoint& local) { return local; });
    BOOST_REQUIRE(report.ok());

    guard.emit(formatMachineReadable(report));
    captured = capture.readAll();
  }

  BOOST_CHECK(captured.find("noise") == std::string::npos);
  const auto payload = requireSingleLine(captured);
  const auto document = parseOrFail(payload);
  BOOST_REQUIRE(document.HasMember("status"));
  BOOST_CHECK_EQUAL(std::string{document["status"].GetString()}, "Ok");
  BOOST_REQUIRE(document.HasMember("surfaces"));
  BOOST_REQUIRE(document["surfaces"].IsArray());
  BOOST_CHECK_EQUAL(document["surfaces"].Size(), 1u);
}

BOOST_AUTO_TEST_CASE(InvalidArgumentTypedFailureProducesExactlyOneJsonDocument)
{
  std::string captured;
  {
    CapturedStdout capture;
    JsonStdoutGuard guard;
    BOOST_REQUIRE(guard.ok());

    // Mirrors main()'s pre-geometry argument-validation failure path: no
    // geometry loading is ever reached, so this also proves the guard
    // itself (not just what it redirects around) produces clean output.
    GeometryProvenance provenance;
    const auto report = buildFailedReport(provenance, ValidatorStatus::InvalidArguments);
    guard.emit(formatMachineReadable(report));
    captured = capture.readAll();
  }

  const auto payload = requireSingleLine(captured);
  const auto document = parseOrFail(payload);
  BOOST_REQUIRE(document.HasMember("status"));
  BOOST_CHECK_EQUAL(std::string{document["status"].GetString()}, "InvalidArguments");
  BOOST_REQUIRE(document.HasMember("surfaces"));
  BOOST_REQUIRE(document["surfaces"].IsArray());
  BOOST_CHECK_EQUAL(document["surfaces"].Size(), 0u);
}

BOOST_AUTO_TEST_CASE(MissingGeometryTypedFailureProducesExactlyOneJsonDocument)
{
  std::string captured;
  {
    CapturedStdout capture;
    JsonStdoutGuard guard;
    BOOST_REQUIRE(guard.ok());

    GeometryProvenance provenance;
    provenance.geometryFilePath = "/nonexistent/o2sim_geometry.root";
    provenance.detectorLabel = "ITS";
    const auto report = buildFailedReport(provenance, ValidatorStatus::GeometryLoadFailed);
    guard.emit(formatMachineReadable(report));
    captured = capture.readAll();
  }

  const auto payload = requireSingleLine(captured);
  const auto document = parseOrFail(payload);
  BOOST_REQUIRE(document.HasMember("status"));
  BOOST_CHECK_EQUAL(std::string{document["status"].GetString()}, "GeometryLoadFailed");
  BOOST_REQUIRE(document.HasMember("geometryFile"));
  BOOST_CHECK_EQUAL(std::string{document["geometryFile"].GetString()}, "/nonexistent/o2sim_geometry.root");
}
