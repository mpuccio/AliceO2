// Gate 4 L3: one public owning graph and one device-facing graph view.

#define BOOST_TEST_MODULE ITSMFT L3 SurfaceGraph guard
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
  return path.extension() == ".h" || path.extension() == ".hpp" || path.extension() == ".c" ||
         path.extension() == ".cc" || path.extension() == ".cxx";
}
} // namespace

BOOST_AUTO_TEST_CASE(OnlySurfaceGraphPublicGraphModelRemains)
{
  const auto root = trackingRoot();
  BOOST_REQUIRE(fs::exists(root / "include/ITSMFTTracking/SurfaceGraph.h"));
  BOOST_REQUIRE(fs::exists(root / "include/ITSMFTTracking/SurfaceGraphBuilder.h"));
  for (const auto oldName : {"DetectorLayout.h", "DetectorLayoutBuilder.h", "DetectorLayoutSet.h", "SparseTrackingTopology.h"}) {
    BOOST_CHECK_MESSAGE(!fs::exists(root / ("include/ITSMFTTracking/" + std::string{oldName})),
                        "deleted graph representation was recreated: " << oldName);
  }
}

BOOST_AUTO_TEST_CASE(NoCompetingGraphVocabularyInCommonProduction)
{
  const auto root = trackingRoot();
  static constexpr std::array<std::string_view, 7> forbidden{
    "DetectorLayout", "DetectorLayoutView", "DetectorLayoutConfigurationKey", "TrackingTopology<",
    "SparseTrackingTopology", "using CATracker", "SurfaceGraphSet"};
  for (const auto directory : {root / "include", root / "src"}) {
    for (const auto& entry : fs::recursive_directory_iterator(directory)) {
      if (!entry.is_regular_file() || !isSource(entry.path())) {
        continue;
      }
      const auto source = readFile(entry.path());
      BOOST_REQUIRE_MESSAGE(!source.empty(), "cannot inspect " << entry.path().string());
      for (const auto token : forbidden) {
        BOOST_CHECK_MESSAGE(source.find(token) == std::string::npos,
                            "obsolete graph representation in " << entry.path().string() << ": " << token);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(SurfaceGraphViewIsTheOnlyDeviceView)
{
  const auto source = readFile(trackingRoot() / "include/ITSMFTTracking/SurfaceGraph.h");
  BOOST_CHECK(source.find("struct SurfaceGraphView") != std::string::npos);
  BOOST_CHECK(source.find("class SurfaceGraph") != std::string::npos);
  BOOST_CHECK(source.find("struct DetectorLayoutView") == std::string::npos);
  BOOST_CHECK(source.find("struct SparseTrackingTopologyView") == std::string::npos);
}
