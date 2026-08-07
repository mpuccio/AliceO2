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

bool containsIdentifier(std::string_view source, std::string_view identifier)
{
  const auto isIdentifierCharacter = [](char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_'; };
  for (std::size_t offset = source.find(identifier); offset != std::string_view::npos; offset = source.find(identifier, offset + 1)) {
    const bool startsIdentifier = offset == 0 || !isIdentifierCharacter(source[offset - 1]);
    const auto end = offset + identifier.size();
    const bool endsIdentifier = end == source.size() || !isIdentifierCharacter(source[end]);
    if (startsIdentifier && endsIdentifier) {
      return true;
    }
  }
  return false;
}

BOOST_AUTO_TEST_CASE(core_headers_and_sources_have_no_typed_output_dependency)
{
  const auto root = trackingRoot();
  static constexpr std::array<std::string_view, 4> coreFiles{
    "include/ITSMFTTracking/Tracker.h",
    "include/ITSMFTTracking/TrackerTraits.h",
    "src/Tracker.cxx",
    "src/TrackerTraits.cxx",
  };
  static constexpr std::array<std::string_view, 13> forbidden{
    "DetectorTraits",
    "CATrackType",
    "LayerMeasurementSpans",
    "TrackITSExt",
    "TrackMFT",
    "MFT"
    "CATrack",
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
  static constexpr std::array<std::string_view, 7> forbidden{
    "DetectorTraits<",
    "CATrackType<",
    "LayerMeasurementSpans<",
    "Tracker<7",
    "Tracker<10",
    "MFT"
    "AdapterRefit",
    "MFT"
    "CATrack",
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
        const auto found = token == ("MFT"
                                     "CATrack")
                             ? containsIdentifier(source, token)
                             : source.find(token) != std::string::npos;
        BOOST_CHECK_MESSAGE(!found,
                            entry.path().lexically_relative(root).generic_string() << " retains " << token);
      }
    }
  }
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/DetectorTraits.h"));
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/MFT"
                                 "AdapterRefit.h"));
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/MFT"
                                 "CATrack.h"));
  const auto trackingCMake = readFile(root / "CMakeLists.txt");
  BOOST_CHECK(trackingCMake.find("MFT"
                                 "AdapterRefit") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(adapter_files_are_the_only_typed_compatibility_seam)
{
  const auto root = trackingRoot();
  const auto refit = readFile(root / "include/ITSMFTTracking/detail/MFTFwdTrackHelpers.h");
  const auto support = readFile(root / "include/ITSMFTTracking/detail/DetectorTrackingOperationAdapterSupport.h");
  const auto publication = readFile(root / "include/ITSMFTTracking/detail/DetectorPublicationAdapter.h");
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/TrackingOperationAdapter.h"));
  BOOST_CHECK(!containsIdentifier(refit,
                                  "MFT"
                                  "CATrack"));
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/MFT"
                                 "AdapterRefit.h"));
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/MFT"
                                 "CATrack.h"));
  BOOST_CHECK(support.find("DetectorTraits") == std::string::npos);
  BOOST_CHECK(publication.find("TrackITSExt") == std::string::npos);
  BOOST_CHECK(publication.find("TrackMFT") == std::string::npos);
}

} // namespace
