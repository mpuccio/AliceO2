// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// Gate 4 M6g source/dependency guard. The common tracking library contains
// only detector-neutral engine/participant primitives and concrete plan-driven
// participant legs. The combined DPL task is the sole owner of the ITS/MFT
// application schedule and event publication context.
//
// The scan is intentionally scoped to common/tracking/include+src and the
// combined workflow include+src. Frozen legacy ITS code elsewhere in the
// repository is explicitly outside this guard. No M6g coordinator spelling is
// exempted here; ITSMFTLegacyParticipantSet was the final temporary exception
// and is retired by this gate.

#define BOOST_TEST_MODULE ITSMFT M6g participant ownership guard
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <regex>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace
{
fs::path trackingRoot()
{
  const fs::path testFile{__FILE__};
  return testFile.parent_path().parent_path();
}

fs::path workflowRoot()
{
  return trackingRoot().parent_path() / "workflow-combined-ca";
}

bool isProductionSource(const fs::path& path)
{
  static constexpr std::array<const char*, 5> extensions{".h", ".hpp", ".c", ".cc", ".cxx"};
  const auto extension = path.extension().string();
  for (const auto* candidate : extensions) {
    if (extension == candidate) {
      return true;
    }
  }
  return false;
}

std::string readFile(const fs::path& path)
{
  std::ifstream input{path};
  if (!input.good()) {
    return {};
  }
  return {std::istreambuf_iterator<char>{input}, {}};
}

void checkTreeHasNoForbiddenSpellings(const fs::path& root, const std::vector<std::regex>& forbidden)
{
  BOOST_REQUIRE_MESSAGE(fs::is_directory(root), "cannot inspect " << root.string());
  for (const auto& entry : fs::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file() || !isProductionSource(entry.path())) {
      continue;
    }
    const auto text = readFile(entry.path());
    BOOST_REQUIRE_MESSAGE(!text.empty(), "cannot inspect " << entry.path().string());
    for (const auto& pattern : forbidden) {
      BOOST_CHECK_MESSAGE(!std::regex_search(text, pattern),
                          "retired M6g spelling remains in " << entry.path().string());
    }
  }
}
} // namespace

BOOST_AUTO_TEST_CASE(CommonCAAndCombinedWorkflowHaveNoRetiredCoordinatorOrBridgeSpellings)
{
  const std::vector<std::regex> forbidden{
    std::regex{"\\bLegacyTrackerScratch\\b"},
    std::regex{"\\bDetectorTraversalBinding\\b"},
    std::regex{"\\bLegacyCATrackingParticipant\\b"},
    std::regex{"\\bTrackSeedTpl\\b"},
    std::regex{"\\bITSMFTLegacyParticipantSet\\b"},
    std::regex{"\\bScratchT\\b"},
    std::regex{"\\bBindingT\\b"}};
  const auto root = trackingRoot();
  checkTreeHasNoForbiddenSpellings(root / "include", forbidden);
  checkTreeHasNoForbiddenSpellings(root / "src", forbidden);
  checkTreeHasNoForbiddenSpellings(workflowRoot() / "include", forbidden);
  checkTreeHasNoForbiddenSpellings(workflowRoot() / "src", forbidden);
}

BOOST_AUTO_TEST_CASE(EventPublicationStateIsOwnedByTheCombinedDPLTask)
{
  const auto tracking = trackingRoot();
  const auto workflow = workflowRoot();
  const auto trackingText = readFile(tracking / "include/ITSMFTTracking/TrackingEngine.h") +
                            readFile(tracking / "include/ITSMFTTracking/TrackingParticipant.h") +
                            readFile(tracking / "include/ITSMFTTracking/SurfacePlanTrackingParticipant.h") +
                            readFile(tracking / "src/TrackingEngine.cxx") +
                            readFile(tracking / "src/SurfacePlanTrackingParticipant.cxx");
  BOOST_CHECK(trackingText.find("mITSClock") == std::string::npos);
  BOOST_CHECK(trackingText.find("mMFTClock") == std::string::npos);
  BOOST_CHECK(trackingText.find("mPublicationValid") == std::string::npos);

  const auto workflowText = readFile(workflow / "include/ITSMFTCombinedCAWorkflow/CombinedCATrackerSpec.h") +
                            readFile(workflow / "src/CombinedCATrackerSpec.cxx");
  for (const auto* spelling : {"mITSClock", "mMFTClock", "mPublicationValid", "clearPublicationSidecars",
                               "invalidatePublication", "markPublicationValid"}) {
    BOOST_CHECK_MESSAGE(workflowText.find(spelling) != std::string::npos,
                        "combined DPL task no longer visibly owns " << spelling);
  }
}

BOOST_AUTO_TEST_CASE(RetiredParticipantSetFilesAreGoneAndPlanDrivenParticipantsRemain)
{
  const auto root = trackingRoot();
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/ITSMFTLegacyParticipantSet.h"));
  BOOST_CHECK(!fs::exists(root / "src/ITSMFTLegacyParticipantSet.cxx"));
  BOOST_CHECK(!fs::exists(root / "test/testITSMFTLegacyParticipantSet.cxx"));
  BOOST_CHECK(fs::is_regular_file(root / "include/ITSMFTTracking/SurfacePlanTrackingParticipant.h"));
  BOOST_CHECK(fs::is_regular_file(root / "src/SurfacePlanTrackingParticipant.cxx"));
}
