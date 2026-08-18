// Gate 4 M7b source/dependency guard.
//
// This guard makes the remaining adapter/device boundaries reviewable and
// rejects a new semantically runtime-sized generic-core use.
// The categories are deliberately narrow:
//   fixed-device-abi      fixed value capacity/representation;
//   private-operation     a temporary typed operation boundary;
//   adapter-edge          detector/output/configuration compatibility;
// It scans common production include/src only.  Tests and documentation are
// not authority sources and are intentionally outside this dependency guard.

#define BOOST_TEST_MODULE ITSMFT M7b runtime count authority guard
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace fs = std::filesystem;

namespace
{
enum class CountBoundary { FixedDeviceABI,
                           PrivateOperation,
                           AdapterEdge };

const char* boundaryName(CountBoundary boundary)
{
  switch (boundary) {
    case CountBoundary::FixedDeviceABI:
      return "fixed device ABI/capacity";
    case CountBoundary::PrivateOperation:
      return "temporary private operation implementation";
    case CountBoundary::AdapterEdge:
      return "adapter edge";
  }
  return "unclassified";
}

// Remove comments before looking for code-level NLayers uses.  The source
// files contain extensive migration notes which are useful to humans but are
// not dependencies or count authorities.
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

std::optional<CountBoundary> classify(const fs::path& path, std::string_view codeLine)
{
  const auto name = path.filename().string();
  if (name == "Cell.h") {
    return CountBoundary::FixedDeviceABI;
  }
  if ((name == "TrackerTraits.cxx" || name == "Tracker.cxx" || name == "TraversalTopology.cxx") && codeLine.find(".NLayers") != std::string_view::npos) {
    return CountBoundary::AdapterEdge;
  }
  if (name == "IndexTableConfiguration.h" || name == "IndexTableConfiguration.cxx" ||
      name == "RefitDriver.h" ||
      name == "CandidateFinding.h" ||
      name == "CandidateFinding.cxx") {
    // The old native-cylinder output loop is an adapter edge; its remaining
    // NLayers template uses are otherwise private operation plumbing.
    if ((name == "IndexTableConfiguration.h" || name == "IndexTableConfiguration.cxx") &&
        codeLine.find("params.NLayers") != std::string_view::npos) {
      return CountBoundary::AdapterEdge;
    }
    return CountBoundary::PrivateOperation;
  }
  if (name == "Configuration.h" || name == "MFTFwdTrackHelpers.h" ||
      name == "MFTFwdTrackHelpers.cxx" || name == "SurfaceMeasurement.h" ||
      name == "Tracker.h" || name == "Tracker.cxx" ||
      name == "DetectorPublicationAdapter.h" ||
      name == "DetectorRefitSupport.h" ||
      name == "Configuration.cxx" ||
      name == "ITSMFTDetectorDefinitions.h" ||
      name == "TrackingConfigParam.h" || name == "ITSMFTDetectorDefinitions.h") {
    return CountBoundary::AdapterEdge;
  }
  return std::nullopt;
}

void checkNoRuntimeSizedTemplateUse(const fs::path& path, std::string_view code)
{
  const auto text = std::string{code};
  if (path.filename() == "TrackerTraits.cxx") {
    BOOST_CHECK_MESSAGE(text.find("mTrkParams[iteration].NLayers != activeSurfaceCount") != std::string::npos,
                        "TrackerTraits must retain only the explicit adapter-edge NLayers/count validation");
    BOOST_CHECK_MESSAGE(text.find("std::array<NominalSurfaceMaterial, NLayers>") == std::string::npos,
                        "TrackerTraits retains a plan-sized material array");
    BOOST_CHECK_MESSAGE(text.find("std::array<float, NLayers>") == std::string::npos,
                        "TrackerTraits retains a plan-sized operation array");
    BOOST_CHECK_MESSAGE(text.find("< NLayers") == std::string::npos,
                        "TrackerTraits contains an NLayers-bounded hot loop");
  }
  if (path.filename() == "TimeFrameScratch.cxx") {
    BOOST_CHECK_MESSAGE(text.find(".NLayers") == std::string::npos,
                        "TimeFrameScratch uses TrackingParameters::NLayers as a count authority");
  }
  BOOST_CHECK_MESSAGE(text.find("std::array<SurfaceMeasurement, NLayers>") == std::string::npos,
                      "common refit production code retains a plan-sized host measurement buffer");
}

void scanProductionTree(const fs::path& root)
{
  BOOST_REQUIRE_MESSAGE(fs::is_directory(root), "cannot inspect " << root.string());
  for (const auto& entry : fs::recursive_directory_iterator(root)) {
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
    checkNoRuntimeSizedTemplateUse(entry.path(), code);
    std::size_t lineStart = 0;
    while (lineStart < code.size()) {
      const auto lineEnd = code.find('\n', lineStart);
      const auto line = code.substr(lineStart, lineEnd == std::string::npos ? code.size() - lineStart : lineEnd - lineStart);
      if (line.find("NLayers") != std::string::npos) {
        const auto category = classify(entry.path(), line);
        BOOST_REQUIRE_MESSAGE(category.has_value(),
                              "unclassified common production NLayers use in " << entry.path().string() << ": " << line);
        BOOST_TEST_MESSAGE(entry.path().filename().string() << ": " << boundaryName(*category));
      }
      if (lineEnd == std::string::npos) {
        break;
      }
      lineStart = lineEnd + 1;
    }
  }
}
} // namespace

BOOST_AUTO_TEST_CASE(AllRemainingNLayersUsesAreClassifiedAndRuntimeHotLoopsAreClean)
{
  const auto root = trackingRoot();
  scanProductionTree(root / "include");
  scanProductionTree(root / "src");
}

BOOST_AUTO_TEST_CASE(RuntimePlanAndFixedCapacityContractsAreVisibleInProductionSources)
{
  const auto root = trackingRoot();
  const auto traits = readFile(root / "src/TrackerTraits.cxx");
  const auto scratch = readFile(root / "src/TimeFrameScratch.cxx");
  const auto refit = readFile(root / "include/ITSMFTTracking/RefitDriver.h");
  BOOST_CHECK(traits.find("context.configuration.topology.orderedSurfaces") != std::string::npos);
  BOOST_CHECK(traits.find("context.configuration.edges") != std::string::npos);
  BOOST_CHECK(scratch.find("mNOwnedSurfaces") == std::string::npos);
  BOOST_CHECK(refit.find("std::vector<detail::RefitMeasurementSlot>") != std::string::npos);
  BOOST_CHECK(refit.find("std::array<SurfaceMeasurement, NLayers>") == std::string::npos);
}
