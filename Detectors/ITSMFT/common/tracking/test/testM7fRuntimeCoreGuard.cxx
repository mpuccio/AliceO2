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
      relative == "include/ITSMFTTracking/ITSSurfaceSpec.h" ||
      relative == "include/ITSMFTTracking/MFTSurfaceSpec.h" ||
      relative == "include/ITSMFTTracking/StaticDetectorCatalogs.h" ||
      relative == "include/ITSMFTTracking/NominalSurfaceMaterialDefaults.h") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "detector application constants and static surface descriptors"};
  }
  if (relative == "src/Configuration.cxx") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "configuration serialization/default construction"};
  }
  if (relative == "src/DetectorLayoutSet.cxx") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "adapter plan validation against TrackingParameters"};
  }
  if (relative == "src/IndexTableConfiguration.cxx") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "adapter-edge active-count validation"};
  }
  if (relative == "src/TrackerTraits.cxx") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "one adapter-edge NLayers/active-surface validation"};
  }
  if (relative == "include/ITSMFTTracking/DetectorPublicationAdapter.h") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "typed publication sidecar adapter"};
  }
  if (relative == "include/ITSMFTTracking/SurfacePlanTrackingParticipant.h" ||
      relative == "src/SurfacePlanTrackingParticipant.cxx") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "ITS/MFT participant, ROF, refit, and publication edge"};
  }
  if (relative == "include/ITSMFTTracking/TrackingInterface.h" ||
      relative == "src/TrackingInterface.cxx") {
    return ResidualClassification{ResidualKind::AdapterCompatibility, "ITS/MFT workflow-facing interface and frozen ROF-table builder"};
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
  static constexpr std::array<std::string_view, 10> coreFiles{
    "include/ITSMFTTracking/Tracker.h",
    "include/ITSMFTTracking/TrackerTraits.h",
    "include/ITSMFTTracking/CATracker.h",
    "include/ITSMFTTracking/TrackingEngine.h",
    "include/ITSMFTTracking/TrackingParticipant.h",
    "include/ITSMFTTracking/TrackingOperationAdapter.h",
    "include/ITSMFTTracking/detail/SurfacePlanBinding.h",
    "include/ITSMFTTracking/SparseTrackingTopology.h",
    "include/ITSMFTTracking/IndexTableUtils.h",
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
    "MFTCATrack",
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

BOOST_AUTO_TEST_CASE(RuntimePlanAndFixedCapacityAuthoritiesRemainExplicit)
{
  const auto root = trackingRoot();
  const auto traits = readFile(root / "src/TrackerTraits.cxx");
  const auto scratch = readFile(root / "include/ITSMFTTracking/SurfaceTrackingScratch.h");
  const auto binding = readFile(root / "include/ITSMFTTracking/detail/SurfacePlanBinding.h");
  const auto seed = readFile(root / "include/ITSMFTTracking/Cell.h");
  BOOST_CHECK(traits.find("getNOwnedSurfaces()") != std::string::npos);
  BOOST_CHECK(traits.find("getOrderedSurfaces()") != std::string::npos);
  BOOST_CHECK(scratch.find("getNOwnedSurfaces()") != std::string::npos);
  BOOST_CHECK(binding.find("getOwnedSurfaceIndex") != std::string::npos);
  BOOST_CHECK(seed.find("MaxLayoutSurfaces") != std::string::npos);
  BOOST_CHECK(seed.find("class CellSeed") != std::string::npos);
}
