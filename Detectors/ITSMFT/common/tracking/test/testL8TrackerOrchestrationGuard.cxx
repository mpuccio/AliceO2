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
  const auto source = read(trackingRoot / "src/Tracker.cxx");
  BOOST_CHECK(header.find("TrackingResult run(TimeFrame& frame, TrackerTraits& traits)") != std::string::npos);
  BOOST_CHECK(header.find("TimeFrame*") == std::string::npos);
  BOOST_CHECK(header.find("TrackerTraits*") == std::string::npos);
  BOOST_CHECK(header.find("SurfaceTrackingScratch m") == std::string::npos);
  BOOST_CHECK(header.find("std::vector<SurfaceGraph>") == std::string::npos);
  BOOST_CHECK(source.find("mFrame") == std::string::npos);
  BOOST_CHECK(source.find("mTraits") == std::string::npos);
  BOOST_CHECK(source.find("mGraphs") == std::string::npos);
  BOOST_CHECK(source.find("mMemoryPool") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(CombinedWorkflowComposesDirectly)
{
  const fs::path trackingRoot{fs::path{__FILE__}.parent_path().parent_path()};
  const fs::path workflowRoot = trackingRoot.parent_path() / "workflow-combined-ca";
  const auto header = read(workflowRoot / "include/ITSMFTCombinedCAWorkflow/CombinedCATrackerSpec.h");
  const auto source = read(workflowRoot / "src/CombinedCATrackerSpec.cxx");
  BOOST_CHECK(header.find("std::unique_ptr<o2::itsmft::tracking::Tracker>") != std::string::npos);
  BOOST_CHECK(header.find("std::unique_ptr<o2::itsmft::tracking::TrackerTraits>") != std::string::npos);
  BOOST_CHECK(source.find("MultiSourceTimeFrameLoader::load") != std::string::npos);
  BOOST_CHECK(source.find("->run(mFrame") != std::string::npos);
  BOOST_CHECK(source.find("mEngine") == std::string::npos);
}
