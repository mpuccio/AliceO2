// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// M6b (Detectors/ITSMFT/common/tracking/doc/design/0002-m6-generic-workspace-migration.md
// Sec 3.2, 7, 9) source-level/include guard: SurfacePlanBinding.h itself must
// never reference ITS/MFT detector identity, the fixed 7/10 layer counts, or
// the fixed ClusterSourceId{0}/{1} combined-load positions -- all of that
// stays inside DetectorTraversalBinding (unmodified) and, further up the
// stack, inside whatever adapter eventually calls SurfacePlanBinding::build()
// (M6d/M6e; today, nothing). Scans only detail/SurfacePlanBinding.h itself
// (the "new binding implementation" the task names), not the whole
// production tree -- DetectorTraversalBinding.h and TrackerTraits.cxx
// legitimately still use every one of these tokens and are unaffected by
// this guard. Mirrors the scan/allowlist idiom of
// testNoLegacyFittingDependency.cxx (ADR 0008) and
// testTransitionPolicyTagContainment.cxx (M4), scoped to a single file since
// this guard has no exceptions to allowlist.

#define BOOST_TEST_MODULE ITSMFT SurfacePlanBindingNoDetectorDependency
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <regex>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace
{

struct TokenCheck {
  std::regex pattern;
  std::string description;
};

const std::vector<TokenCheck>& forbiddenTokens()
{
  static const std::vector<TokenCheck> checks = {
    {std::regex{"DetID\\s*::\\s*ITS"}, "hardcoded ITS detector identity"},
    {std::regex{"DetID\\s*::\\s*MFT"}, "hardcoded MFT detector identity"},
    {std::regex{"\\bDetID\\s*::"}, "any DetID enumerator (this type accepts no detector identity at all)"},
    {std::regex{"\\bITSNLayers\\b"}, "fixed ITS layer count (7)"},
    {std::regex{"\\bMFTNLayers\\b"}, "fixed MFT layer count (10)"},
    {std::regex{"ClusterSourceId\\s*\\{\\s*0\\s*\\}"}, "fixed source-0 fact (the combined-load ITS position)"},
    {std::regex{"ClusterSourceId\\s*\\{\\s*1\\s*\\}"}, "fixed source-1 fact (the combined-load MFT position)"},
  };
  return checks;
}

fs::path surfacePlanBindingHeader()
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  return fs::path(testDirectory) / ".." / "include" / "ITSMFTTracking" / "detail" / "SurfacePlanBinding.h";
}

} // namespace

BOOST_AUTO_TEST_CASE(SurfacePlanBindingHeaderMentionsNoDetectorIdentitySourceOrFixedLayerCountFacts)
{
  const auto path = surfacePlanBindingHeader();
  BOOST_REQUIRE_MESSAGE(fs::is_regular_file(path), "cannot find " << path.string());

  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path.string());

  std::vector<std::string> lines;
  std::string line;
  while (std::getline(input, line)) {
    lines.push_back(line);
  }
  BOOST_CHECK_GT(lines.size(), 0u);

  for (const auto& check : forbiddenTokens()) {
    for (size_t i = 0; i < lines.size(); ++i) {
      BOOST_CHECK_MESSAGE(!std::regex_search(lines[i], check.pattern),
                          "SurfacePlanBinding.h:" << (i + 1) << " mentions a forbidden token ("
                                                  << check.description << "): " << lines[i]);
    }
  }
}

BOOST_AUTO_TEST_CASE(SurfacePlanBindingHeaderDoesNotIncludeTheDetIdHeader)
{
  const auto path = surfacePlanBindingHeader();
  BOOST_REQUIRE_MESSAGE(fs::is_regular_file(path), "cannot find " << path.string());

  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path.string());

  const std::regex includeRegex{"#include\\s*\"DetectorsCommonDataFormats/DetID\\.h\""};
  std::string line;
  size_t lineNumber = 0;
  bool found = false;
  while (std::getline(input, line)) {
    ++lineNumber;
    if (std::regex_search(line, includeRegex)) {
      found = true;
      BOOST_CHECK_MESSAGE(false, "SurfacePlanBinding.h:" << lineNumber << " includes the DetID header, which this type must not depend on: " << line);
    }
  }
  BOOST_CHECK(!found);
}
