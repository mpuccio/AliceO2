// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Source guards for shared tracklet/cell orchestration and the remaining
// neighbour/road operation boundary.

#define BOOST_TEST_MODULE ITSMFT TraversalOperationBindingContainment
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace
{
namespace fs = std::filesystem;

fs::path commonTrackingRoot()
{
  return fs::path{__FILE__}.parent_path().parent_path();
}

std::string readFile(const std::string& path)
{
  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path);
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

std::string readTrackerTraitsSource()
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  return readFile(testDirectory + "/../src/TrackerTraits.cxx");
}

std::string readTrackerTraitsHeader()
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  return readFile(testDirectory + "/../include/ITSMFTTracking/TrackerTraits.h");
}

std::string readTrackletFindingSource()
{
  return readFile((commonTrackingRoot() / "src/TrackletFinding.cxx").string());
}

/// Extracts the body (inclusive of the opening/closing braces) of the first
/// `void TrackerTraits::<methodName>(` definition found in `source`,
/// by counting braces from that signature's opening `{` to its matching
/// closing `}`. `methodName` must be immediately followed by `(` in the
/// match (a plain string search on the qualified name, not a regex), so
/// `computeLayerTracklets` never matches `computeLayerTrackletsForKind`.
std::string extractMethodBody(const std::string& source, const std::string& methodName)
{
  const std::string signature = "TrackerTraits::" + methodName + "(";
  const auto signaturePos = source.find(signature);
  BOOST_REQUIRE_MESSAGE(signaturePos != std::string::npos, "cannot find " << signature << " in TrackerTraits.cxx");
  const auto openBrace = source.find('{', signaturePos);
  BOOST_REQUIRE_MESSAGE(openBrace != std::string::npos, "cannot find opening brace for " << methodName);

  int depth = 0;
  size_t pos = openBrace;
  for (; pos < source.size(); ++pos) {
    if (source[pos] == '{') {
      ++depth;
    } else if (source[pos] == '}') {
      --depth;
      if (depth == 0) {
        break;
      }
    }
  }
  BOOST_REQUIRE_MESSAGE(depth == 0, "unbalanced braces while extracting " << methodName);
  return source.substr(openBrace, pos - openBrace + 1);
}

/// Strips `//`-style line comments (this file's own convention throughout
/// TrackerTraits.cxx -- no `/* */` block comments are used there), so a
/// token search below only ever matches actual code, never explanatory prose
/// about what a method used to do.
std::string stripLineComments(const std::string& body)
{
  std::istringstream lines{body};
  std::ostringstream stripped;
  std::string line;
  while (std::getline(lines, line)) {
    const auto commentPos = line.find("//");
    stripped << (commentPos == std::string::npos ? line : line.substr(0, commentPos)) << '\n';
  }
  return stripped.str();
}

bool mentionsToken(const std::string& body, const std::string& token)
{
  const std::regex tokenRegex{"\\b" + token + "\\b"};
  return std::regex_search(body, tokenRegex);
}

size_t countOccurrences(const std::string& text, const std::string& needle)
{
  size_t count = 0;
  for (size_t pos = 0; (pos = text.find(needle, pos)) != std::string::npos; pos += needle.size()) {
    ++count;
  }
  return count;
}

} // namespace

BOOST_AUTO_TEST_CASE(TrackletAndCellEntryPointsUseOneGlobalTraversal)
{
  const auto source = readTrackerTraitsSource();
  const std::array<std::tuple<std::string, std::string, std::string>, 2> methods{{
    {"computeLayerTracklets", "computeLayerTrackletsImpl", "getGlobalTransitions"},
    {"computeLayerCells", "computeLayerCellsImpl", "getGlobalCells"},
  }};
  for (const auto& [method, implementation, globalIds] : methods) {
    const auto code = stripLineComments(extractMethodBody(source, method));
    BOOST_CHECK_EQUAL(countOccurrences(code, implementation + "("), 1u);
    BOOST_CHECK_EQUAL(countOccurrences(code, globalIds + "()"), 1u);
    BOOST_CHECK(!mentionsToken(code, "SurfaceKind"));
  }
}

BOOST_AUTO_TEST_CASE(NeighbourAndRoadEntryPointsKeepTheirExistingLeafBoundary)
{
  const auto source = readTrackerTraitsSource();
  const std::array<std::tuple<std::string, std::string, std::string>, 2> methods{{
    {"findCellsNeighbours", "findCellsNeighboursCylinder", "findCellsNeighboursDisk"},
    {"findRoads", "findRoadsCylinder", "findRoadsDisk"},
  }};
  const std::array<std::string, 3> forbiddenTokens{"SurfaceKind", "StateFamily", "constexpr"};

  for (const auto& [method, cylinderLeaf, diskLeaf] : methods) {
    const auto body = extractMethodBody(source, method);
    BOOST_REQUIRE_GT(body.size(), 0u);
    const auto code = stripLineComments(body);
    for (const auto& token : forbiddenTokens) {
      BOOST_CHECK_MESSAGE(!mentionsToken(code, token),
                          "TrackerTraits::" << method << "() still mentions " << token
                                            << " in code -- a SurfaceKind/StateFamily dispatch branch leaked back"
                                            << " into the shared entry point");
    }
    BOOST_CHECK(code.find(cylinderLeaf) != std::string::npos);
    BOOST_CHECK(code.find(diskLeaf) != std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(RetiredTraversalOperationAdapterDoesNotReturn)
{
  const auto source = readTrackerTraitsSource();
  const auto header = readTrackerTraitsHeader();
  for (const auto token : {"TraversalOperationBinding", "mTraversalOperation", "bindTraversalOperation"}) {
    BOOST_CHECK(header.find(token) == std::string::npos);
    BOOST_CHECK(source.find(token) == std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(OperationPartitionsComeFromSurfaceDescriptorsOnly)
{
  const auto source = readTrackerTraitsSource();
  const auto code = stripLineComments(extractMethodBody(source, "initialiseTimeFrame"));
  BOOST_CHECK(code.find("layout.getSurface(layout.getTransition(transitionId).from).kind") != std::string::npos);
  BOOST_CHECK(code.find("mTransitionsByKind[kindIndex(kind)].push_back(transitionId)") != std::string::npos);
  BOOST_CHECK(code.find("mCellsByKind") == std::string::npos);
  BOOST_CHECK(code.find("mScheduledCellsByKind[kindIndex(kind)].push_back(cellId)") != std::string::npos);
  BOOST_CHECK(code.find("mRoadStartCellsByKind[kindIndex(kind)].push_back(cellId)") != std::string::npos);
  for (const auto token : {"ClusterSourceId", "DetID", "parametersByKind", "parametersForKind"}) {
    BOOST_CHECK_MESSAGE(!mentionsToken(code, token),
                        "operation partitioning depends on " << token << " instead of SurfaceDescriptor::kind");
  }
}

BOOST_AUTO_TEST_CASE(RemainingNonTemplateWrapperTargetsForwardToKindSpecificLeafImplementations)
{
  const auto source = readTrackerTraitsSource();
  const std::array<std::pair<std::string, std::string>, 4> wrapperToLeaf{{
    {"findCellsNeighboursCylinder", "findCellsNeighboursForKind<SurfaceKind::Cylinder>"},
    {"findCellsNeighboursDisk", "findCellsNeighboursForKind<SurfaceKind::Disk>"},
    {"findRoadsCylinder", "findRoadsForKind<SurfaceKind::Cylinder>"},
    {"findRoadsDisk", "findRoadsForKind<SurfaceKind::Disk>"},
  }};
  for (const auto& [wrapper, leafCall] : wrapperToLeaf) {
    const auto body = extractMethodBody(source, wrapper);
    BOOST_REQUIRE_GT(body.size(), 0u);
    BOOST_CHECK_MESSAGE(body.find(leafCall) != std::string::npos,
                        "TrackerTraits::" << wrapper << "() does not forward to " << leafCall);
  }
  for (const auto retired : {"computeLayerTrackletsCylinder", "computeLayerTrackletsDisk",
                             "computeLayerCellsCylinder", "computeLayerCellsDisk",
                             "computeLayerTrackletsForKind", "computeLayerCellsForKind"}) {
    BOOST_CHECK(source.find(retired) == std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(ProcessNeighboursKeepsSurfaceKindSelectionCompileTime)
{
  const auto source = readTrackerTraitsSource();
  BOOST_CHECK(source.find("template <SurfaceKind Kind, typename InputSeed>") != std::string::npos);
  BOOST_CHECK(source.find("params.kind") == std::string::npos);
  BOOST_CHECK(source.find("processNeighbours<Kind>(") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(RetiredCoordinateCutsAreAbsentFromCommonProductionSources)
{
  const std::array<fs::path, 2> productionRoots{
    commonTrackingRoot() / "include",
    commonTrackingRoot() / "src"};
  const std::array<std::string, 6> retired{
    "TrackletMinAbsX", "trackletMinAbsX", "CellRoadRCut", "cellRoadRCut",
    "passesCylinderCellRoadPrecut", "passesDiskCellRoadPrecut"};

  for (const auto& root : productionRoots) {
    for (const auto& entry : fs::recursive_directory_iterator{root}) {
      if (!entry.is_regular_file() || (entry.path().extension() != ".h" && entry.path().extension() != ".cxx")) {
        continue;
      }
      const auto source = readFile(entry.path().string());
      for (const auto& token : retired) {
        BOOST_CHECK_MESSAGE(source.find(token) == std::string::npos,
                            entry.path() << " still contains retired coordinate cut " << token);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(CellCandidateLoopHasOneDescriptorSelectedLeafBoundary)
{
  const auto source = readTrackerTraitsSource();
  const auto code = stripLineComments(extractMethodBody(source, "computeLayerCellsImpl"));

  BOOST_CHECK_EQUAL(countOccurrences(code, "buildCellSeed("), 1u);
  BOOST_CHECK(code.find("buildCylinderCellSeed(") == std::string::npos);
  BOOST_CHECK(code.find("buildDiskCellSeed(") == std::string::npos);
  BOOST_CHECK(code.find("SurfaceKind::Cylinder") == std::string::npos);
  BOOST_CHECK(code.find("SurfaceKind::Disk") == std::string::npos);
  for (const auto token : {"DetID", "ClusterSourceId", "ROAD", "Precut"}) {
    BOOST_CHECK_MESSAGE(!mentionsToken(code, token),
                        "cell orchestration contains detector/source/cut dispatch token " << token);
  }
}

BOOST_AUTO_TEST_CASE(TrackletLoopUsesOnlyCoordinateNeutralLeafFacades)
{
  const auto source = readTrackerTraitsSource();
  const auto code = stripLineComments(extractMethodBody(source, "computeLayerTrackletsImpl"));
  for (const auto required : {"bindTrackletProjectionState(", "projectTrackletSearchWindow(",
                              "trackletSearchRowBin(", "acceptTrackletCandidate("}) {
    BOOST_CHECK_EQUAL(countOccurrences(code, required), 1u);
  }
  for (const auto forbidden : {"projectCylinderSearchWindow(", "projectDiskSearchWindow(",
                               "CylinderTrackletSearchWindow", "DiskTrackletSearchWindow",
                               "SurfaceKind::Cylinder", "SurfaceKind::Disk",
                               "DetID", "ClusterSourceId"}) {
    BOOST_CHECK(code.find(forbidden) == std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(KernelPolicyIsOneSurfaceNeutralRecord)
{
  const auto header = readTrackerTraitsHeader();
  const auto kernelHeader = readFile((commonTrackingRoot() / "include/ITSMFTTracking/detail/TrackingKernelParameters.h").string());
  BOOST_CHECK(header.find("TrackingKernelParameters mKernelParameters") != std::string::npos);
  BOOST_CHECK(header.find("array<TrackingKernelParameters") == std::string::npos);
  BOOST_CHECK(kernelHeader.find("SurfaceKind") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(SurfaceSelectionLivesInTrackletAndCellLeafSource)
{
  const auto source = readTrackletFindingSource();
  BOOST_CHECK(source.find("case SurfaceKind::Cylinder:") != std::string::npos);
  BOOST_CHECK(source.find("case SurfaceKind::Disk:") != std::string::npos);
  BOOST_CHECK(source.find("projectCylinderSearchWindow(") != std::string::npos);
  BOOST_CHECK(source.find("projectDiskSearchWindow(") != std::string::npos);
  BOOST_CHECK(source.find("buildCylinderCellSeed(") != std::string::npos);
  BOOST_CHECK(source.find("buildDiskCellSeed(") != std::string::npos);
}
