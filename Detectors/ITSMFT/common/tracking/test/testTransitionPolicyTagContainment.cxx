// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// M4 (GenericTrackingEngineMigration.md; ADR 0007 decisions 7-8) compile/grep
// boundary guard: TransitionPolicyTag is a temporary legacy hot-loop-dispatch
// implementation detail, confined to ITSMFTTracking/detail/ (headers whose
// consumers are, by construction, only other ITSMFTTracking translation
// units -- never a DPL workflow spec or another detector). No public,
// non-detail header under ITSMFTTracking/ may name TransitionPolicyTag.
//
// This scans every *.h file directly under include/ITSMFTTracking/ (not
// recursing into detail/, which is expected to mention the token) so the
// guard also catches a *future* header that starts referencing the tag --
// not just the files known at the time this test was written.

#define BOOST_TEST_MODULE ITSMFT TransitionPolicyTagContainment
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <set>
#include <string>

namespace fs = std::filesystem;

namespace
{

// Explicit, reasoned allowlist: every one of these is confirmed (repo-wide
// grep, M4 audit) to have zero consumers outside this same library -- never
// included by a DPL workflow spec, another detector, or any other component
// -- so the tag/family-templated internals they still declare (all
// unaffected leaf hot-loop machinery this milestone leaves behind for M5,
// per GenericTrackingEngineMigration.md's M4 entry) are not part of any
// adapter-facing contract, even though the file itself is not physically
// under detail/.
const std::set<std::string>& allowlist()
{
  static const std::set<std::string> files = {
    // TrackerTraits<NLayers> is the legacy CA participant implementation
    // (M2/M3): every tag-templated member left here is private
    // (processNeighbours<Tag> included, moved private by this milestone);
    // the public surface uses only StateFamily (getPolicyBindingCount).
    "TrackerTraits.h",
    // Tag-templated legacy index-table binder (bindIndexTableConfiguration<Tag>);
    // declaration-only, host-only, called exclusively from TrackerTraits.cxx.
    "IndexTableConfiguration.h",
    // Unwired native-refit driver; its one internal `constexpr auto Tag =
    // TransitionPolicyTag::CylinderCylinder` selects a detail/
    // TransitionPolicyOperations.h specialization, never exposed in this
    // file's own public signatures.
    "NativeCylinderCylinderRefitDriver.h",
  };
  return files;
}

} // namespace

BOOST_AUTO_TEST_CASE(NoPublicNonDetailHeaderExposesTheLegacyTransitionPolicyTag)
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const fs::path includeDir = fs::path(testDirectory) / ".." / "include" / "ITSMFTTracking";
  BOOST_REQUIRE_MESSAGE(fs::is_directory(includeDir), "cannot find " << includeDir.string());

  // detail/ is expected to still declare/define TransitionPolicyTag and its
  // dispatch/grouping/tag-keyed-traits machinery -- that is exactly what
  // "confined to detail/" means -- so it is deliberately not scanned here.
  const fs::path detailDir = includeDir / "detail";
  BOOST_REQUIRE_MESSAGE(fs::is_directory(detailDir), "detail/ boundary directory is missing: " << detailDir.string());

  size_t scanned = 0;
  for (const auto& entry : fs::directory_iterator(includeDir)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".h") {
      continue;
    }
    const auto name = entry.path().filename().string();
    ++scanned;
    if (allowlist().count(name) != 0) {
      continue;
    }
    std::ifstream input{entry.path()};
    BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << entry.path().string());
    std::string line;
    size_t lineNumber = 0;
    while (std::getline(input, line)) {
      ++lineNumber;
      BOOST_CHECK_MESSAGE(line.find("TransitionPolicyTag") == std::string::npos,
                          name << ":" << lineNumber << " (non-detail, non-allowlisted) mentions TransitionPolicyTag: " << line);
    }
  }
  // Sanity: this must actually have scanned a realistic number of headers --
  // an empty/near-empty loop (e.g. from a resolved-but-wrong includeDir)
  // would make every BOOST_CHECK_MESSAGE above vacuously pass.
  BOOST_CHECK_GT(scanned, 40u);
}
