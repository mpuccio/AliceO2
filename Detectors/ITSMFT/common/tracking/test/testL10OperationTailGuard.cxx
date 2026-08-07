// Gate 4 L10: the operation seam is refit-only and retired compatibility
// paths cannot re-enter the common production/test boundary.

#define BOOST_TEST_MODULE ITSMFT L10 operation tail guard
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>

#include <boost/test/unit_test.hpp>

namespace fs = std::filesystem;

namespace
{
fs::path trackingRoot()
{
  return fs::path{__FILE__}.parent_path().parent_path();
}

std::string readFile(const fs::path& path)
{
  std::ifstream input{path};
  return {std::istreambuf_iterator<char>{input}, {}};
}

bool isSource(const fs::path& path)
{
  const auto extension = path.extension().string();
  return extension == ".h" || extension == ".hpp" || extension == ".c" || extension == ".cc" || extension == ".cxx";
}

void checkNoTokenInTree(const fs::path& root, std::string_view token)
{
  for (const auto& entry : fs::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file() || !isSource(entry.path())) {
      continue;
    }
    if (entry.path().filename() == "testL10OperationTailGuard.cxx") {
      continue;
    }
    const auto source = readFile(entry.path());
    BOOST_CHECK_MESSAGE(source.find(token) == std::string::npos,
                        entry.path().lexically_relative(root) << " retains retired L10 token " << token);
  }
}

BOOST_AUTO_TEST_CASE(operation_seam_contains_only_refit)
{
  const auto root = trackingRoot();
  const auto adapter = readFile(root / "include/ITSMFTTracking/TrackingOperationAdapter.h");
  BOOST_REQUIRE(!adapter.empty());
  BOOST_CHECK(adapter.find("refitSeed") != std::string::npos);
  BOOST_CHECK(adapter.find("completeAccepted") == std::string::npos);
  BOOST_CHECK(adapter.find("resetAdapterState") == std::string::npos);

  for (const auto& relative : {"src/Tracker.cxx", "src/TrackerTraits.cxx", "include/ITSMFTTracking/Tracker.h",
                               "include/ITSMFTTracking/TrackerTraits.h"}) {
    const auto source = readFile(root / relative);
    BOOST_REQUIRE(!source.empty());
    BOOST_CHECK(source.find("completeAccepted") == std::string::npos);
    BOOST_CHECK(source.find("resetAdapterState") == std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(retired_public_compatibility_paths_are_absent)
{
  const auto root = trackingRoot();
  checkNoTokenInTree(root / "include", "TrackingFrameInfoAdapters");
  checkNoTokenInTree(root / "src", "TrackingFrameInfoAdapters");
  checkNoTokenInTree(root / "test", "TrackingFrameInfoAdapters");
  checkNoTokenInTree(root / "include", "loadClusterTrackingFrameInfo");
  checkNoTokenInTree(root / "src", "loadClusterTrackingFrameInfo");
  checkNoTokenInTree(root / "test", "loadClusterTrackingFrameInfo");
  checkNoTokenInTree(root / "include", "NativeCylinderRefitDriver");
  checkNoTokenInTree(root / "src", "NativeCylinderRefitDriver");
  checkNoTokenInTree(root / "test", "NativeCylinderRefitDriver");
  checkNoTokenInTree(root / "include", "getDevicePropagator");
  checkNoTokenInTree(root / "src", "getDevicePropagator");
  checkNoTokenInTree(root / "test", "getDevicePropagator");
  checkNoTokenInTree(root / "include", "setDevicePropagator");
  checkNoTokenInTree(root / "src", "setDevicePropagator");
  checkNoTokenInTree(root / "test", "setDevicePropagator");
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/TrackingFrameInfoAdapters.h"));
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/NativeCylinderRefitDriver.h"));
}

BOOST_AUTO_TEST_CASE(publication_lifecycle_is_at_the_application_edge)
{
  const auto root = trackingRoot();
  const auto its = readFile(root.parent_path().parent_path() / "ITS/workflow-ca/src/CATrackerSpec.cxx");
  const auto mft = readFile(root.parent_path().parent_path() / "MFT/workflow/src/CATrackerSpec.cxx");
  const auto combined = readFile(root.parent_path().parent_path() / "common/workflow-combined-ca/src/CombinedCATrackerSpec.cxx");
  for (const auto source : {its, mft, combined}) {
    BOOST_CHECK(source.find("completeAccepted") != std::string::npos);
    BOOST_CHECK(source.find("PublicationAdapter.reset") != std::string::npos);
  }
}
} // namespace
