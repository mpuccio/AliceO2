// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// Gate 4 M6f source/dependency guard. The common-CA production include/src
// tree has one workspace, one plan binding, and one whole-track seed. This
// guard scans that tree directly so a compatibility template, bridge type,
// or stale production reference cannot be reintroduced unnoticed.
//
// The scan is intentionally scoped to common/tracking/include and
// common/tracking/src. Frozen legacy ITS implementation outside this common
// tree is not part of M6f, and ITSMFTLegacyParticipantSet remains the one
// explicitly permitted M6g-only application-adapter class whose name still
// contains "Legacy".

#define BOOST_TEST_MODULE ITSMFT M6f legacy bridge guard
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <regex>
#include <string>

namespace fs = std::filesystem;

namespace
{
fs::path trackingRoot()
{
  const fs::path testFile{__FILE__};
  return testFile.parent_path().parent_path();
}

bool isProductionSource(const fs::path& path)
{
  static constexpr std::array<const char*, 5> extensions{ ".h", ".hpp", ".c", ".cc", ".cxx" };
  const auto extension = path.extension().string();
  for (const auto* candidate : extensions) {
    if (extension == candidate) {
      return true;
    }
  }
  return false;
}
} // namespace

BOOST_AUTO_TEST_CASE(CommonCAProductionHasNoRetiredWorkspaceBindingOrSeedSpellings)
{
  const auto root = trackingRoot();
  const auto includeRoot = root / "include";
  const auto sourceRoot = root / "src";
  BOOST_REQUIRE_MESSAGE(fs::is_directory(includeRoot), "cannot inspect " << includeRoot.string());
  BOOST_REQUIRE_MESSAGE(fs::is_directory(sourceRoot), "cannot inspect " << sourceRoot.string());

  const std::array<std::regex, 6> forbidden{
    std::regex{"\\bLegacyTrackerScratch\\b"},
    std::regex{"\\bDetectorTraversalBinding\\b"},
    std::regex{"\\bLegacyCATrackingParticipant\\b"},
    std::regex{"\\bTrackSeedTpl\\b"},
    std::regex{"\\bScratchT\\b"},
    std::regex{"\\bBindingT\\b"}};

  for (const auto& directory : {includeRoot, sourceRoot}) {
    for (const auto& entry : fs::recursive_directory_iterator(directory)) {
      if (!entry.is_regular_file() || !isProductionSource(entry.path())) {
        continue;
      }
      std::ifstream input{entry.path()};
      BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << entry.path().string());
      const std::string text{std::istreambuf_iterator<char>{input}, {}};
      for (const auto& pattern : forbidden) {
        BOOST_CHECK_MESSAGE(!std::regex_search(text, pattern),
                            "retired M6f bridge spelling remains in " << entry.path().string());
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(RetiredBridgeFilesAreGoneAndPlanParticipantIsNamedNarrowly)
{
  const auto root = trackingRoot();
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/LegacyTrackerScratch.h"));
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/detail/DetectorTraversalBinding.h"));
  BOOST_CHECK(!fs::exists(root / "src/LegacyTrackerScratch.cxx"));
  BOOST_CHECK(fs::is_regular_file(root / "include/ITSMFTTracking/SurfacePlanTrackingParticipant.h"));
  BOOST_CHECK(fs::is_regular_file(root / "src/SurfacePlanTrackingParticipant.cxx"));
}
