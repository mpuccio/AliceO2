// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#define BOOST_TEST_MODULE ITSMFT L8 Tracker orchestration guard
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

namespace
{
namespace fs = std::filesystem;

std::string read(const fs::path& path)
{
  std::ifstream input{path};
  return {std::istreambuf_iterator<char>{input}, {}};
}

bool isSource(const fs::path& path)
{
  static constexpr std::array<const char*, 5> extensions{".h", ".hpp", ".c", ".cc", ".cxx"};
  const auto suffix = path.extension().string();
  for (const auto* extension : extensions) {
    if (suffix == extension) {
      return true;
    }
  }
  return false;
}

void assertNoSpellings(const fs::path& root, const std::vector<std::regex>& forbidden)
{
  BOOST_REQUIRE(fs::is_directory(root));
  for (const auto& entry : fs::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file() || !isSource(entry.path())) {
      continue;
    }
    const auto text = read(entry.path());
    for (const auto& pattern : forbidden) {
      BOOST_CHECK_MESSAGE(!std::regex_search(text, pattern), "retired orchestration spelling in " << entry.path());
    }
  }
}
} // namespace

BOOST_AUTO_TEST_CASE(RetiredEngineAndParticipantCompositionIsAbsent)
{
  const fs::path trackingRoot{fs::path{__FILE__}.parent_path().parent_path()};
  const fs::path workflowRoot = trackingRoot.parent_path() / "workflow-combined-ca";
  const std::vector<std::regex> forbidden{
    std::regex{"\\bTrackingEngine\\b"},
    std::regex{"\\bTrackingParticipant\\b"},
    std::regex{"\\bSurfacePlanTrackingParticipant\\b"},
    std::regex{"\\bParticipantId\\b"},
    std::regex{"\\bmEngine\\b"},
    std::regex{"\\bmSchedule\\b"},
    std::regex{"\\bexecuteEvent\\b"}};
  assertNoSpellings(trackingRoot / "include", forbidden);
  assertNoSpellings(trackingRoot / "src", forbidden);
  assertNoSpellings(workflowRoot / "include", forbidden);
  assertNoSpellings(workflowRoot / "src", forbidden);
}

BOOST_AUTO_TEST_CASE(TrackerIsStatelessAndOwnsOnlyTheOperationEdge)
{
  const fs::path trackingRoot{fs::path{__FILE__}.parent_path().parent_path()};
  const auto header = read(trackingRoot / "include/ITSMFTTracking/Tracker.h");
  BOOST_CHECK(header.find("TrackingResult run(TimeFrame& frame, TrackerTraits& traits)") != std::string::npos);
  BOOST_CHECK(header.find("TimeFrame*") == std::string::npos);
  BOOST_CHECK(header.find("TrackerTraits*") == std::string::npos);
  BOOST_CHECK(header.find("TimeFrameScratch m") == std::string::npos);
  BOOST_CHECK(header.find("std::vector<" + std::string{"Surface"} + "Graph>") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(TraversalTestAccessIsConfinedToTestSupport)
{
  const fs::path trackingRoot{fs::path{__FILE__}.parent_path().parent_path()};
  const auto tracker = read(trackingRoot / "include/ITSMFTTracking/Tracker.h");
  const auto traits = read(trackingRoot / "include/ITSMFTTracking/TrackerTraits.h");
  const auto trackerFriend = tracker.find("friend struct TrackerTestAccess");
  const auto traitsFriend = traits.find("friend struct TrackerTestAccess");
  BOOST_REQUIRE_NE(trackerFriend, std::string::npos);
  BOOST_REQUIRE_NE(traitsFriend, std::string::npos);
  BOOST_CHECK_EQUAL(tracker.find("friend struct TrackerTestAccess", trackerFriend + 1), std::string::npos);
  BOOST_CHECK_EQUAL(traits.find("friend struct TrackerTestAccess", traitsFriend + 1), std::string::npos);

  const std::regex access{"\\bTrackerTestAccess\\b"};
  for (const auto& root : {trackingRoot / "include", trackingRoot / "src"}) {
    BOOST_REQUIRE(fs::is_directory(root));
    for (const auto& entry : fs::recursive_directory_iterator(root)) {
      if (!entry.is_regular_file() || !isSource(entry.path())) {
        continue;
      }
      if (entry.path() == trackingRoot / "include/ITSMFTTracking/Tracker.h" ||
          entry.path() == trackingRoot / "include/ITSMFTTracking/TrackerTraits.h") {
        continue;
      }
      const auto text = read(entry.path());
      BOOST_CHECK_MESSAGE(!std::regex_search(text, access), "test-only traversal access leaked into " << entry.path());
    }
  }
}

BOOST_AUTO_TEST_CASE(TraversalTraitsKeepsNoAdoptedTraversalState)
{
  const fs::path trackingRoot{fs::path{__FILE__}.parent_path().parent_path()};
  const auto traits = read(trackingRoot / "include/ITSMFTTracking/TrackerTraits.h");
  const std::vector<std::regex> retired{
    std::regex{"\\badopt(Frame|Scratch)\\b"},
    std::regex{"\\b(updateTrackingParameters|setBz|setMemoryPool|resetTraversalCache|hasTraversalCache)\\b"},
    std::regex{"\\b(acceptedTracksForSharedStatus|clearAcceptedTracksForSharedStatus|initialiseTimeFrame)\\b"},
    std::regex{"\\bm(Frame|Scratch|Binding|TrkParams|Bz|MemoryPool|TraversalGraph|TraversalCacheValid|KernelParameters|AttachHitConfig|LayerMaterial|LayerMeasurements|LayerGlobalMeasurements|DiskLayerReferenceZ|AcceptedTracksForSharedStatus)\\b"}};
  for (const auto& pattern : retired) {
    BOOST_CHECK_MESSAGE(!std::regex_search(traits, pattern), "retired TrackerTraits state/API remains in its header");
  }
  BOOST_CHECK(traits.find("TraversalWorkspaceView") != std::string::npos);
  BOOST_CHECK(traits.find("mTaskArena") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(CombinedWorkflowComposesDirectly)
{
  const fs::path trackingRoot{fs::path{__FILE__}.parent_path().parent_path()};
  const fs::path workflowRoot = trackingRoot.parent_path() / "workflow-combined-ca";
  const auto header = read(workflowRoot / "include/ITSMFTCombinedCAWorkflow/CombinedCATrackerSpec.h");
  const auto source = read(workflowRoot / "src/CombinedCATrackerSpec.cxx");
  BOOST_CHECK(header.find("std::unique_ptr<o2::itsmft::tracking::Tracker>") != std::string::npos);
  BOOST_CHECK(header.find("std::unique_ptr<o2::itsmft::tracking::TrackerTraits>") != std::string::npos);
  BOOST_CHECK(source.find("loadTimeFrameSources") != std::string::npos);
  BOOST_CHECK(source.find("->run(mFrame") != std::string::npos);
  BOOST_CHECK(source.find("mEngine") == std::string::npos);
}
