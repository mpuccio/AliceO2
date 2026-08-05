// Gate 4 M7e: typed refit/output compatibility is adapter-owned, while the
// Tracker/TrackerTraits execution path remains detector-neutral.

#define BOOST_TEST_MODULE ITSMFT M7e adapter boundary guard
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
  if (!input.good()) {
    return {};
  }
  return {std::istreambuf_iterator<char>{input}, {}};
}

BOOST_AUTO_TEST_CASE(core_headers_and_sources_have_no_typed_output_dependency)
{
  const auto root = trackingRoot();
  static constexpr std::array<std::string_view, 6> coreFiles{
    "include/ITSMFTTracking/Tracker.h",
    "include/ITSMFTTracking/TrackerTraits.h",
    "include/ITSMFTTracking/CATracker.h",
    "include/ITSMFTTracking/TrackingOperationAdapter.h",
    "src/CATracker.cxx",
    "src/TrackerTraits.cxx",
  };
  static constexpr std::array<std::string_view, 13> forbidden{
    "DetectorTraits",
    "CATrackType",
    "LayerMeasurementSpans",
    "TrackITSExt",
    "TrackMFT",
    "MFTCATrack",
    "DetID",
    "o2::detectors",
    "DPL",
    "Workflow",
    "TrackWriter",
    "sharedClusters",
    "hasSharedClusters",
  };

  for (const auto relative : coreFiles) {
    const auto path = root / relative;
    const auto source = readFile(path);
    BOOST_REQUIRE_MESSAGE(!source.empty(), "cannot inspect " << path.string());
    for (const auto token : forbidden) {
      BOOST_CHECK_MESSAGE(source.find(token) == std::string::npos,
                          relative << " names adapter/output dependency " << token);
    }
  }
}

BOOST_AUTO_TEST_CASE(deleted_typed_bridges_have_no_production_definition)
{
  const auto root = trackingRoot();
  static constexpr std::array<std::string_view, 5> forbidden{
    "DetectorTraits<",
    "CATrackType<",
    "LayerMeasurementSpans<",
    "Tracker<7",
    "Tracker<10",
  };
  for (const auto& directory : {root / "include", root / "src"}) {
    for (const auto& entry : fs::recursive_directory_iterator(directory)) {
      if (!entry.is_regular_file()) {
        continue;
      }
      const auto extension = entry.path().extension().string();
      if (extension != ".h" && extension != ".hpp" && extension != ".c" && extension != ".cc" && extension != ".cxx") {
        continue;
      }
      const auto source = readFile(entry.path());
      BOOST_REQUIRE_MESSAGE(!source.empty(), "cannot inspect " << entry.path().string());
      for (const auto token : forbidden) {
        BOOST_CHECK_MESSAGE(source.find(token) == std::string::npos,
                            entry.path().lexically_relative(root).generic_string() << " retains " << token);
      }
    }
  }
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/DetectorTraits.h"));
}

BOOST_AUTO_TEST_CASE(adapter_files_are_the_only_typed_compatibility_seam)
{
  const auto root = trackingRoot();
  const auto refit = readFile(root / "include/ITSMFTTracking/MFTFwdTrackHelpers.h");
  const auto adapterRefit = readFile(root / "include/ITSMFTTracking/MFTAdapterRefit.h");
  const auto support = readFile(root / "include/ITSMFTTracking/DetectorTrackingOperationAdapterSupport.h");
  const auto publication = readFile(root / "include/ITSMFTTracking/DetectorPublicationAdapter.h");
  BOOST_CHECK(refit.find("MFTCATrack") == std::string::npos);
  BOOST_CHECK(adapterRefit.find("MFTCATrack") != std::string::npos);
  BOOST_CHECK(support.find("DetectorTraits") == std::string::npos);
  BOOST_CHECK(publication.find("TrackITSExt") == std::string::npos);
  BOOST_CHECK(publication.find("TrackMFT") == std::string::npos);
}

} // namespace
