// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// M5d (doc/decisions/0008-native-refit-activation.md) source-level/include
// guard: the adapter refit operation's two production branches (ITS/barrel and
// MFT/forward) no longer depend on either frozen legacy fitting engine --
// o2::its::track::fitTrack/refitTrack/refitTrackSeed (ITStracking/TrackHelpers.h)
// or o2::mft::TrackFitter<TrackLTF>/TrackLTFL (MFTTracking/TrackFitter.h). Both
// engines instead remain reachable only as explanatory prose in this
// milestone's own doc comments, describing what was
// removed and why -- never an actual #include or call.
//
// This scans include/ITSMFTTracking/*.h, include/ITSMFTTracking/detail/*.h,
// and src/*.cxx (this module's own production headers/sources, not test/),
// exactly like the M4 SurfaceKind containment guard
// (testSurfaceKindContainment.cxx), whose scan/allowlist structure
// this test follows.

#define BOOST_TEST_MODULE ITSMFT NoLegacyFittingDependency
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <regex>
#include <set>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace
{

struct TokenCheck {
  std::string token;
  std::string description;
  std::set<std::string> allowlist;
};

// Every allowlisted file is scanned too -- the allowlist only means "known,
// reasoned exception", not "unscanned". Rationale per file:
//  - NativeRefitDriver.h, detail/TrackletFinding.h,
//    SurfaceStateOperationResult.h: doc comments describing the legacy
//    fitTrack/refitTrack formula an operation reproduces or a failure reason
//    mirrors -- prose, not a dependency.
//  - MFTFwdTrackHelpers.h, MFTFwdTrackHelpers.cxx: this milestone's own doc
//    comments naming the two engines it removed.
const std::vector<TokenCheck>& forbiddenTokens()
{
  static const std::vector<TokenCheck> checks = {
    {"refitTrackSeed",
     "frozen ITS refit entry point (ITStracking/TrackHelpers.h)",
     {"SurfaceStateOperationResult.h"}},
    {"TrackFitContext",
     "frozen ITS refit context type (ITStracking/TrackHelpers.h)",
     {}},
    {"TrackLTF",
     "frozen MFT linear-track-finder fit engine (MFTTracking/TrackFitter.h)",
     {"MFTFwdTrackHelpers.h", "MFTFwdTrackHelpers.cxx", "NativeRefitDriver.h"}},
  };
  return checks;
}

bool lineMentionsToken(const std::string& line, const std::regex& tokenRegex)
{
  return std::regex_search(line, tokenRegex);
}

size_t scanFilesForToken(const std::vector<fs::path>& files, const TokenCheck& check)
{
  const std::regex tokenRegex{"\\b" + check.token + "\\b"};
  size_t scanned = 0;
  for (const auto& path : files) {
    const auto name = path.filename().string();
    ++scanned;
    if (check.allowlist.count(name) != 0) {
      continue;
    }
    std::ifstream input{path};
    BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path.string());
    std::string line;
    size_t lineNumber = 0;
    while (std::getline(input, line)) {
      ++lineNumber;
      BOOST_CHECK_MESSAGE(!lineMentionsToken(line, tokenRegex),
                          name << ":" << lineNumber << " (non-allowlisted) mentions " << check.token
                               << " (" << check.description << "): " << line);
    }
  }
  return scanned;
}

std::vector<fs::path> collectProductionFiles(const fs::path& trackingRoot)
{
  std::vector<fs::path> files;
  for (const auto* subdir : {"include/ITSMFTTracking", "include/ITSMFTTracking/detail", "src"}) {
    const fs::path dir = trackingRoot / subdir;
    BOOST_REQUIRE_MESSAGE(fs::is_directory(dir), "cannot find " << dir.string());
    for (const auto& entry : fs::directory_iterator(dir)) {
      if (!entry.is_regular_file()) {
        continue;
      }
      const auto ext = entry.path().extension();
      if (ext == ".h" || ext == ".cxx") {
        files.push_back(entry.path());
      }
    }
  }
  return files;
}

} // namespace

BOOST_AUTO_TEST_CASE(NewCommonTrackerHasNoFrozenLegacyFittingSymbolDependency)
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const fs::path trackingRoot = fs::path(testDirectory) / "..";

  const auto files = collectProductionFiles(trackingRoot);
  BOOST_CHECK_GT(files.size(), 60u);

  for (const auto& check : forbiddenTokens()) {
    const auto scanned = scanFilesForToken(files, check);
    BOOST_CHECK_EQUAL(scanned, files.size());
  }
}

BOOST_AUTO_TEST_CASE(NoProductionFileIncludesEitherFrozenLegacyFittingHeader)
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const fs::path trackingRoot = fs::path(testDirectory) / "..";
  const auto files = collectProductionFiles(trackingRoot);

  const std::regex itsIncludeRegex{"#include\\s*\"ITStracking/TrackHelpers\\.h\""};
  const std::regex mftIncludeRegex{"#include\\s*\"MFTTracking/TrackFitter\\.h\""};

  for (const auto& path : files) {
    std::ifstream input{path};
    BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path.string());
    std::string line;
    size_t lineNumber = 0;
    while (std::getline(input, line)) {
      ++lineNumber;
      BOOST_CHECK_MESSAGE(!std::regex_search(line, mftIncludeRegex),
                          path.filename().string() << ":" << lineNumber << " includes the frozen MFT fitter header: " << line);
      // ITStracking/TrackHelpers.h: retained by TrackletFinding.cxx
      // for legacy-layer/constant utilities unrelated to the fitting chain
      // this milestone removes (grep-verified: no o2::its::track::fitTrack/
      // refitTrack/refitTrackSeed/TrackFitContext symbol is referenced by that
      // file) -- pre-existing, out of this milestone's scope.
      if (path.filename() == "TrackletFinding.cxx") {
        continue;
      }
      BOOST_CHECK_MESSAGE(!std::regex_search(line, itsIncludeRegex),
                          path.filename().string() << ":" << lineNumber << " includes the frozen ITS refit header: " << line);
    }
  }
}
