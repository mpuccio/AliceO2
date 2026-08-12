// Gate 4 M7f: one runtime-plan core, with only named adapter/device
// compatibility boundaries remaining in common tracking production sources.

#define BOOST_TEST_MODULE ITSMFT M7f runtime core guard
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>

#include <boost/test/unit_test.hpp>

namespace fs = std::filesystem;

namespace
{

enum class ResidualKind { FixedDeviceCapacity,
                          AdapterCompatibility };

struct ResidualClassification {
  ResidualKind kind;
  std::string_view role;
};

const char* kindName(ResidualKind kind)
{
  switch (kind) {
    case ResidualKind::FixedDeviceCapacity:
      return "fixed device ABI/capacity";
    case ResidualKind::AdapterCompatibility:
      return "narrow adapter-owned compatibility";
  }
  return "unclassified";
}

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
  const auto isIdentifierCharacter = [](char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_';
  };
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

// Migration comments deliberately document the deleted bridges. They are not
// production dependencies, so remove them before classifying the remaining
// code-level NLayers uses.
std::string withoutComments(std::string_view source)
{
  std::string result;
  result.reserve(source.size());
  bool blockComment = false;
  for (std::size_t i = 0; i < source.size();) {
    if (blockComment) {
      if (i + 1 < source.size() && source[i] == '*' && source[i + 1] == '/') {
        blockComment = false;
        result += "  ";
        i += 2;
      } else {
        result += source[i++] == '\n' ? '\n' : ' ';
      }
      continue;
    }
    if (i + 1 < source.size() && source[i] == '/' && source[i + 1] == '*') {
      blockComment = true;
      result += "  ";
      i += 2;
    } else if (i + 1 < source.size() && source[i] == '/' && source[i + 1] == '/') {
      while (i < source.size() && source[i] != '\n') {
        result += ' ';
        ++i;
      }
    } else {
      result += source[i++];
    }
  }
  return result;
}

std::optional<ResidualClassification> classifyNLayers(std::string_view relative)
{
  // These are the only common-tracking production files in which a live
  // compile-time layer parameter remains. Every one is an application-edge
  // bridge: frozen ROF/configuration construction or participant/publication
  // compatibility. The generic Tracker/TrackerTraits path is not listed.
  if (relative == "include/ITSMFTTracking/Configuration.h") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "configuration compatibility accessors and DetID-to-layer adapter mapping"};
  }
  if (relative == "include/ITSMFTTracking/TrackingConfigParam.h" ||
      relative == "include/ITSMFTTracking/ITSMFTDetectorDefinitions.h" ||
      relative == "include/ITSMFTTracking/ITSMFTDetectorDefinitions.h") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "detector application constants and static surface descriptors"};
  }
  if (relative == "src/Configuration.cxx") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "configuration serialization/default construction"};
  }
  if (relative == "src/SurfaceGraphBuilder.cxx") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "adapter plan validation against TrackingParameters"};
  }
  if (relative == "src/IndexTableConfiguration.cxx") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "adapter-edge active-count validation"};
  }
  if (relative == "src/TrackerTraits.cxx") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "one adapter-edge NLayers/active-surface validation"};
  }
  if (relative == "src/Tracker.cxx") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "adapter-edge NLayers/active-surface validation"};
  }
  if (relative == "include/ITSMFTTracking/detail/DetectorPublicationAdapter.h") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "typed publication sidecar adapter"};
  }
  return std::nullopt;
}

void scanResidualNLayers(const fs::path& root, const fs::path& directory)
{
  BOOST_REQUIRE_MESSAGE(fs::is_directory(directory), "cannot inspect " << directory.string());
  for (const auto& entry : fs::recursive_directory_iterator(directory)) {
    if (!entry.is_regular_file()) {
      continue;
    }
    const auto extension = entry.path().extension().string();
    if (extension != ".h" && extension != ".hpp" && extension != ".c" && extension != ".cc" && extension != ".cxx") {
      continue;
    }
    const auto raw = readFile(entry.path());
    BOOST_REQUIRE_MESSAGE(!raw.empty(), "cannot inspect " << entry.path().string());
    const auto code = withoutComments(raw);
    if (code.find("NLayers") == std::string::npos) {
      continue;
    }
    const auto relative = entry.path().lexically_relative(root).generic_string();
    const auto classification = classifyNLayers(relative);
    BOOST_REQUIRE_MESSAGE(classification.has_value(),
                          "unclassified common-production NLayers use in " << relative);
    BOOST_TEST_MESSAGE(relative << ": " << kindName(classification->kind) << " (" << classification->role << ")");
  }
}

void scanForbiddenProductionVocabulary(const fs::path& root)
{
  static constexpr std::array<std::string_view, 16> forbidden{
    "TrackSeedTpl",
    "DetectorTraversalBinding",
    "LegacyTrackerScratch",
    "DetectorTraits<",
    "CATrackType<",
    "LayerMeasurementSpans<",
    "Tracker<",
    "TrackerTraits<",
    "Tracker<7",
    "Tracker<10",
    "TrackerTraits<7",
    "TrackerTraits<10",
    "CATrackerITS",
    "CATrackerMFT",
    "TrackerITS",
    "TrackerMFT",
  };
  static constexpr std::array<std::string_view, 7> deletedBridges{
    "IndexTableUtilsN",
    "mSurfaceToSlot",
    "stagedSurfaceToSlot",
    "kInvalidSurfaceSlot",
    "surfaceToLegacyLayer",
    "legacyLayerToSurface",
    "exportNativeRefitToTrackITSExt",
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
      const auto raw = readFile(entry.path());
      BOOST_REQUIRE_MESSAGE(!raw.empty(), "cannot inspect " << entry.path().string());
      const auto code = withoutComments(raw);
      const auto relative = entry.path().lexically_relative(root).generic_string();
      for (const auto token : forbidden) {
        BOOST_CHECK_MESSAGE(code.find(token) == std::string::npos,
                            relative << " retains deleted common-core vocabulary " << token);
      }
      for (const auto token : deletedBridges) {
        BOOST_CHECK_MESSAGE(code.find(token) == std::string::npos,
                            relative << " retains deleted bridge " << token);
      }
    }
  }
}

void scanCoreBoundary(const fs::path& root)
{
  static constexpr std::array<std::string_view, 6> coreFiles{
    "include/ITSMFTTracking/Tracker.h",
    "include/ITSMFTTracking/TrackerTraits.h",
    "include/ITSMFTTracking/detail/SurfacePlanBinding.h",
    "include/ITSMFTTracking/SurfaceGraph.h",
    "include/ITSMFTTracking/IndexTableConfiguration.h",
    "src/TrackerTraits.cxx",
  };
  static constexpr std::array<std::string_view, 13> forbidden{
    "o2::detectors::DetID",
    "ITSNLayers",
    "MFTNLayers",
    "ROFOverlapTable",
    "ROFVertexLookupTable",
    "ROFMaskTable",
    "TrackITSExt",
    "CATrack",
    "DetectorTraits",
    "CATrackType",
    "LayerMeasurementSpans",
    "DPL",
    "Writer",
  };

  for (const auto relative : coreFiles) {
    const auto path = root / relative;
    const auto raw = readFile(path);
    BOOST_REQUIRE_MESSAGE(!raw.empty(), "cannot inspect " << path.string());
    const auto code = withoutComments(raw);
    for (const auto token : forbidden) {
      BOOST_CHECK_MESSAGE(code.find(token) == std::string::npos,
                          relative << " selects detector/ROF/output behavior through " << token);
    }
  }
}

void scanCanonicalTrackerFiles(const fs::path& root)
{
  const auto trackerHeader = root / "include/ITSMFTTracking/Tracker.h";
  const auto trackerSource = root / "src/Tracker.cxx";
  BOOST_REQUIRE(fs::exists(trackerHeader));
  BOOST_REQUIRE(fs::exists(trackerSource));
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/CATracker.h"));
  BOOST_CHECK(!fs::exists(root / "src/CATracker.cxx"));

  const auto trackingCMake = readFile(root / "CMakeLists.txt");
  BOOST_CHECK(trackingCMake.find("src/Tracker.cxx") != std::string::npos);
  BOOST_CHECK(trackingCMake.find("CATracker.h") == std::string::npos);
  BOOST_CHECK(trackingCMake.find("CATracker.cxx") == std::string::npos);

  for (const auto& directory : {root / "include", root / "src"}) {
    for (const auto& entry : fs::recursive_directory_iterator(directory)) {
      if (!entry.is_regular_file()) {
        continue;
      }
      BOOST_CHECK_MESSAGE(entry.path().filename().string().find("CATracker") == std::string::npos,
                          "obsolete CATracker filename remains: " << entry.path().string());
      const auto extension = entry.path().extension().string();
      if (extension != ".h" && extension != ".hpp" && extension != ".c" && extension != ".cc" && extension != ".cxx") {
        continue;
      }
      const auto code = withoutComments(readFile(entry.path()));
      BOOST_CHECK_MESSAGE(code.find("CATracker.h") == std::string::npos,
                          entry.path().string() << " retains the deleted CATracker header include/name");
      BOOST_CHECK_MESSAGE(code.find("CATracker.cxx") == std::string::npos,
                          entry.path().string() << " retains the deleted CATracker source name");
      BOOST_CHECK_MESSAGE(!containsIdentifier(code, "CATracker"),
                          entry.path().string() << " retains the obsolete CATracker identifier");
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(AllResidualNLayersUsesHaveAnExactNamedException)
{
  const auto root = trackingRoot();
  scanResidualNLayers(root, root / "include");
  scanResidualNLayers(root, root / "src");
}

BOOST_AUTO_TEST_CASE(DeletedBridgesAndOldCoreVocabularyAreAbsent)
{
  scanForbiddenProductionVocabulary(trackingRoot());
}

BOOST_AUTO_TEST_CASE(GenericCoreHasNoDetectorLayerOrROFCompatibilityAuthority)
{
  scanCoreBoundary(trackingRoot());
}

BOOST_AUTO_TEST_CASE(TrackerFilesAreCanonicalAndCATrackerIsDeleted)
{
  scanCanonicalTrackerFiles(trackingRoot());
}

BOOST_AUTO_TEST_CASE(RuntimePlanAndFixedCapacityAuthoritiesRemainExplicit)
{
  const auto root = trackingRoot();
  const auto traits = readFile(root / "src/TrackerTraits.cxx");
  const auto scratch = readFile(root / "include/ITSMFTTracking/detail/SurfaceTrackingScratch.h");
  const auto binding = readFile(root / "include/ITSMFTTracking/detail/SurfacePlanBinding.h");
  const auto seed = readFile(root / "include/ITSMFTTracking/Cell.h");
  BOOST_CHECK(traits.find("getNOwnedSurfaces()") != std::string::npos);
  BOOST_CHECK(traits.find("getOrderedSurfaces()") != std::string::npos);
  BOOST_CHECK(scratch.find("getNOwnedSurfaces()") != std::string::npos);
  BOOST_CHECK(binding.find("getOwnedSurfaceIndex") != std::string::npos);
  BOOST_CHECK(seed.find("MaxLayoutSurfaces") != std::string::npos);
  BOOST_CHECK(seed.find("class CellSeed") != std::string::npos);
}
