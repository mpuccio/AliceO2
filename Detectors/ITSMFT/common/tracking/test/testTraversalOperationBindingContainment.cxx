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

std::string readTrackerSource()
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  return readFile(testDirectory + "/../src/Tracker.cxx");
}

std::string readTrackerTraitsHeader()
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  return readFile(testDirectory + "/../include/ITSMFTTracking/TrackerTraits.h");
}

std::string readCandidateFindingSource()
{
  return readFile((commonTrackingRoot() / "src/CandidateFinding.cxx").string());
}

/// Extracts the body (inclusive of the opening/closing braces) of the first
/// `void TrackerTraits::<methodName>(` definition found in `source`,
/// by counting braces from that signature's opening `{` to its matching
/// closing `}`. `methodName` must be immediately followed by `(` in the
/// match (a plain string search on the qualified name, not a regex), so
/// `computeLayerTracklets` never matches `computeLayerTrackletsForKind`.
std::string extractMethodBody(const std::string& source, const std::string& methodName,
                              const std::string& owner = "TrackerTraits")
{
  const std::string signature = owner + "::" + methodName + "(";
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
    {"computeLayerTracklets", "computeLayerTrackletsImpl", "context.configuration.edges"},
    {"computeLayerCells", "computeLayerCellsImpl", "context.configuration.cells"},
  }};
  for (const auto& [method, implementation, globalIds] : methods) {
    const auto code = stripLineComments(extractMethodBody(source, method));
    BOOST_CHECK_EQUAL(countOccurrences(code, implementation + "("), 1u);
    BOOST_CHECK_EQUAL(countOccurrences(code, globalIds), 1u);
    BOOST_CHECK(!mentionsToken(code, "SurfaceKind"));
  }
}

BOOST_AUTO_TEST_CASE(NeighbourEntryPointUsesOneCoordinateNeutralCompatibilityPath)
{
  const auto source = readTrackerTraitsSource();
  const std::array<std::string, 3> forbiddenTokens{"SurfaceKind", "SurfaceKind", "constexpr"};
  const auto entryPoint = stripLineComments(extractMethodBody(source, "findCellsNeighbours"));
  BOOST_CHECK_EQUAL(countOccurrences(entryPoint, "findCellsNeighboursForSchedule("), 1u);
  for (const auto& token : forbiddenTokens) {
    BOOST_CHECK_MESSAGE(!mentionsToken(entryPoint, token),
                        "TrackerTraits::findCellsNeighbours() still mentions " << token);
  }

  const auto compatibility = stripLineComments(extractMethodBody(source, "findCellsNeighboursForSchedule"));
  BOOST_REQUIRE_GT(compatibility.size(), 0u);
  BOOST_CHECK_EQUAL(countOccurrences(compatibility, "fitAdjacentTripletFactors("), 1u);
  for (const auto& token : {"SurfaceKind", "SurfaceKind", "ClusterSourceId", "DetID", "constexpr"}) {
    BOOST_CHECK_MESSAGE(!mentionsToken(compatibility, token),
                        "cell-neighbour compatibility depends on " << token);
  }
}

BOOST_AUTO_TEST_CASE(RoadEntryPointUsesOneGraphSchedule)
{
  const auto source = readTrackerTraitsSource();
  const auto code = stripLineComments(extractMethodBody(source, "findRoads"));
  BOOST_REQUIRE_GT(code.size(), 0u);
  for (const auto token : {"SurfaceKind", "SurfaceKind", "constexpr"}) {
    BOOST_CHECK(!mentionsToken(code, token));
  }
  BOOST_CHECK(code.find("findRoadsImpl") != std::string::npos);
  BOOST_CHECK(code.find("findRoadsCylinder") == std::string::npos);
  BOOST_CHECK(code.find("findRoadsDisk") == std::string::npos);
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

BOOST_AUTO_TEST_CASE(EdgePreparationUsesOneGraphSchedule)
{
  const auto source = readTrackerSource();
  const auto code = stripLineComments(extractMethodBody(source, "initializeIteration", "Tracker"));
  BOOST_CHECK(code.find("prepareTraversalEdgeTolerances") != std::string::npos);
  BOOST_CHECK(code.find("mLinksByKind") == std::string::npos);
  BOOST_CHECK(code.find("mRoadStartCellsByKind") == std::string::npos);
  for (const auto token : {"ClusterSourceId", "DetID", "parametersByKind", "parametersForKind"}) {
    BOOST_CHECK_MESSAGE(!mentionsToken(code, token),
                        "operation partitioning depends on " << token << " instead of SurfaceDescriptor::kind");
  }
}

BOOST_AUTO_TEST_CASE(RetiredRoadWrappersDoNotReturn)
{
  const auto source = readTrackerTraitsSource();
  for (const auto retired : {"findRoadsCylinder", "findRoadsDisk", "findRoadsForKind",
                             "computeLayerTrackletsCylinder", "computeLayerTrackletsDisk",
                             "computeLayerCellsCylinder", "computeLayerCellsDisk"}) {
    BOOST_CHECK(source.find(retired) == std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(ProcessNeighboursUsesTheStateDrivenPropagator)
{
  const auto source = readTrackerTraitsSource();
  BOOST_CHECK(source.find("template <SurfaceKind Kind, typename InputSeed>") == std::string::npos);
  BOOST_CHECK(source.find("params.kind") == std::string::npos);
  const auto code = stripLineComments(extractMethodBody(source, "processNeighbours"));
  BOOST_CHECK(code.find("Propagator::attachMeasurement") != std::string::npos);
  BOOST_CHECK(code.find("SurfaceKind::Cylinder") == std::string::npos);
  BOOST_CHECK(code.find("SurfaceKind::Disk") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(RefitWorkersUseTheDescriptorDrivenBoundary)
{
  const auto source = readTrackerTraitsSource();
  const auto code = stripLineComments(extractMethodBody(source, "findRoadsImpl"));
  BOOST_REQUIRE_GT(code.size(), 0u);
  BOOST_CHECK(code.find("refitSources") == std::string::npos);
  BOOST_CHECK(code.find("TraversalFailureReason::SurfaceKindMismatch") == std::string::npos);
  BOOST_CHECK(code.find("context.configuration.topology.orderedSurfaces") != std::string::npos);
  BOOST_CHECK(code.find("mLayerGlobalMeasurements[position].front()") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(RetiredTraversalForwardersDoNotReturn)
{
  const auto source = readTrackerTraitsSource();
  const auto header = readTrackerTraitsHeader();
  for (const auto token : {"findRoadsForGraph", "findRoadsForSchedule", "LayerGeometryConfigView", "bindLayerGeometryConfig"}) {
    BOOST_CHECK(source.find(token) == std::string::npos);
    BOOST_CHECK(header.find(token) == std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(RetiredCoordinateCutsAreAbsentFromCommonProductionSources)
{
  const std::array<fs::path, 2> productionRoots{
    commonTrackingRoot() / "include",
    commonTrackingRoot() / "src"};
  const std::array<std::string, 10> retired{
    "TrackletMinAbsX", "trackletMinAbsX", "CellRoadRCut", "cellRoadRCut",
    "passesCylinderCellRoadPrecut", "passesDiskCellRoadPrecut",
    "CellDeltaTanLambdaSigma", "cellDeltaTanLambdaSigma",
    "computeDiskCellDirectionCompatibilityChi2",
    "cellDirectionsAreCompatible(SurfaceKind"};

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

  BOOST_CHECK_EQUAL(countOccurrences(code, "cellDirectionsAreCompatible("), 1u);
  BOOST_CHECK_EQUAL(countOccurrences(code, "makeDirectionObservation("), 3u);
  BOOST_CHECK_EQUAL(countOccurrences(code, "trackletDirectionsAreTransverselyCompatible("), 1u);
  BOOST_CHECK_EQUAL(countOccurrences(code, "makeTransverseDirectionObservation("), 3u);
  BOOST_CHECK_EQUAL(countOccurrences(code, "getEdgeMSAngle(secondEdgeId)"), 1u);
  BOOST_CHECK_EQUAL(countOccurrences(code, "DirectionProcessNoise"), 1u);
  BOOST_CHECK_EQUAL(countOccurrences(code, "buildCellSeed("), 1u);
  BOOST_CHECK(code.find("CellDeltaTanLambdaSigma") == std::string::npos);
  BOOST_CHECK(code.find("deltaTanLambdaSigma") == std::string::npos);
  BOOST_CHECK(code.find("buildCylinderCellSeed(") == std::string::npos);
  BOOST_CHECK(code.find("buildDiskCellSeed(") == std::string::npos);
  BOOST_CHECK(code.find("SurfaceKind::Cylinder") == std::string::npos);
  BOOST_CHECK(code.find("SurfaceKind::Disk") == std::string::npos);
  BOOST_CHECK(code.find(".tanLambda") == std::string::npos);
  BOOST_CHECK(code.find("switch (kind)") == std::string::npos);
  for (const auto token : {"DetID", "ClusterSourceId", "ROAD", "Precut"}) {
    BOOST_CHECK_MESSAGE(!mentionsToken(code, token),
                        "cell orchestration contains detector/source/cut dispatch token " << token);
  }
}

BOOST_AUTO_TEST_CASE(TransverseTrackletDirectionUsesOneFamilyNeutralAlgorithm)
{
  const auto source = readCandidateFindingSource();
  const auto begin = source.find("bool trackletDirectionsAreTransverselyCompatible(");
  const auto end = source.find("bool makeDirectionObservation(", begin);
  BOOST_REQUIRE(begin != std::string::npos);
  BOOST_REQUIRE(end != std::string::npos);
  const auto compatibility = source.substr(begin, end - begin);
  for (const auto forbidden : {"SurfaceKind", "switch", "if constexpr",
                               "Cylinder", "Disk", "DetID", "ClusterSourceId"}) {
    BOOST_CHECK_MESSAGE(compatibility.find(forbidden) == std::string::npos,
                        "shared transverse-direction algorithm contains " << forbidden);
  }
  for (const auto required : {"trackletMinPt", "maximumCurvature", "maximumBending",
                              "deltaPhi", "covarianceXY", "angularVariance", "chi2"}) {
    BOOST_CHECK_MESSAGE(compatibility.find(required) != std::string::npos,
                        "shared transverse-direction algorithm lacks " << required);
  }
}

BOOST_AUTO_TEST_CASE(CellDirectionUsesOneFamilyNeutralAlgorithmAndConsolidatedBoundary)
{
  const auto root = commonTrackingRoot();
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/detail/CellFinding.h"));
  BOOST_CHECK(!fs::exists(root / "include/ITSMFTTracking/detail/TrackletFinding.h"));
  BOOST_CHECK(fs::exists(root / "include/ITSMFTTracking/detail/CandidateFinding.h"));

  const auto source = readCandidateFindingSource();
  const auto begin = source.find("bool cellDirectionsAreCompatible(");
  const auto end = source.find("bool buildCylinderCellSeed(", begin);
  BOOST_REQUIRE(begin != std::string::npos);
  BOOST_REQUIRE(end != std::string::npos);
  const auto compatibility = source.substr(begin, end - begin);
  for (const auto forbidden : {"SurfaceKind", "switch", "if constexpr",
                               "Cylinder", "Disk", "DetID", "ClusterSourceId"}) {
    BOOST_CHECK_MESSAGE(compatibility.find(forbidden) == std::string::npos,
                        "shared cell-direction algorithm contains " << forbidden);
  }
  for (const auto required : {"covarianceRZ", "deltaZ01", "deltaZ12",
                              "deltaR01", "deltaR12", "variance", "chi2"}) {
    BOOST_CHECK_MESSAGE(compatibility.find(required) != std::string::npos,
                        "shared cell-direction algorithm lacks " << required);
  }
}

BOOST_AUTO_TEST_CASE(TrackletLoopUsesOnlyCoordinateNeutralLeafFacades)
{
  const auto source = readTrackerTraitsSource();
  const auto code = stripLineComments(extractMethodBody(source, "computeLayerTrackletsImpl"));
  BOOST_CHECK_EQUAL(countOccurrences(code, "bindTrackletProjectionCache("), 2u);
  for (const auto required : {"projectTrackletSearchWindow(", "searchWindow.bins", "acceptTrackletCandidate("}) {
    BOOST_CHECK_EQUAL(countOccurrences(code, required), 1u);
  }
  for (const auto forbidden : {"projectCylinderSearchWindow(", "projectDiskSearchWindow(",
                               "CylinderTrackletSearchWindow", "DiskTrackletSearchWindow",
                               "trackletSearchBins(", "trackletSearchRowCount(", "trackletSearchRowBin(",
                               "SurfaceKind::Cylinder", "SurfaceKind::Disk",
                               "DetID", "ClusterSourceId", "makeLinkState"}) {
    BOOST_CHECK(code.find(forbidden) == std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(TrackletSearchWindowIsOneDataOnlyCovarianceContract)
{
  const auto header = readFile((commonTrackingRoot() / "include/ITSMFTTracking/detail/CandidateFinding.h").string());
  const auto begin = header.find("struct TrackletSearchWindow");
  const auto end = header.find("bool bindTrackletProjectionCache", begin);
  BOOST_REQUIRE(begin != std::string::npos);
  BOOST_REQUIRE(end != std::string::npos);
  const auto window = header.substr(begin, end - begin);
  for (const auto required : {"int4 bins", "float prediction[2]", "float variance[3]"}) {
    BOOST_CHECK(window.find(required) != std::string::npos);
  }
  for (const auto forbidden : {"std::variant", "CylinderTrackletSearchWindow", "DiskTrackletSearchWindow", "nSigmaCut",
                               "projectCylinderSearchWindow", "projectDiskSearchWindow",
                               "buildCylinderCellSeed", "buildDiskCellSeed"}) {
    BOOST_CHECK(window.find(forbidden) == std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(KernelPolicyIsOneSurfaceNeutralRecord)
{
  const auto header = readTrackerTraitsHeader();
  const auto kernelHeader = readFile((commonTrackingRoot() / "include/ITSMFTTracking/detail/TrackingKernelParameters.h").string());
  BOOST_CHECK(header.find("TrackingKernelParameters mKernelParameters") == std::string::npos);
  BOOST_CHECK(header.find("array<TrackingKernelParameters") == std::string::npos);
  BOOST_CHECK(kernelHeader.find("SurfaceKind") == std::string::npos);
}

BOOST_AUTO_TEST_CASE(SurfaceSelectionLivesInCandidateFindingLeaves)
{
  const auto source = readCandidateFindingSource();
  BOOST_CHECK(source.find("case SurfaceKind::Cylinder:") != std::string::npos);
  BOOST_CHECK(source.find("case SurfaceKind::Disk:") != std::string::npos);
  BOOST_CHECK(source.find("projectCylinderSearchWindow(") != std::string::npos);
  BOOST_CHECK(source.find("projectDiskSearchWindow(") != std::string::npos);
  BOOST_CHECK(source.find("buildCylinderCellSeed(") != std::string::npos);
  BOOST_CHECK(source.find("buildDiskCellSeed(") != std::string::npos);
  BOOST_CHECK(source.find("makeDirectionObservation(") != std::string::npos);
  BOOST_CHECK(source.find("makeTransverseDirectionObservation(") != std::string::npos);
  BOOST_CHECK(source.find("makeCylinderDirectionObservation(") == std::string::npos);
  BOOST_CHECK(source.find("makeDiskDirectionObservation(") == std::string::npos);
  BOOST_CHECK(source.find("makeCylinderTransverseDirectionObservation(") == std::string::npos);
  BOOST_CHECK(source.find("makeDiskTransverseDirectionObservation(") == std::string::npos);
}
