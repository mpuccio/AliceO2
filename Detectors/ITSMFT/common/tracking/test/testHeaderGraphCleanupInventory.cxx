// Source guard for the final common-tracking header boundary.

#define BOOST_TEST_MODULE ITSMFT HeaderGraphCleanupInventory
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>

namespace
{
namespace fs = std::filesystem;

fs::path trackingRoot()
{
  return fs::path{__FILE__}.parent_path().parent_path();
}

fs::path itsmftRoot()
{
  return trackingRoot().parent_path().parent_path();
}

std::string readFile(const fs::path& path)
{
  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path.string());
  return {std::istreambuf_iterator<char>{input}, {}};
}

bool isSourceFile(const fs::path& path)
{
  const auto extension = path.extension().string();
  return extension == ".h" || extension == ".hpp" || extension == ".c" || extension == ".cc" ||
         extension == ".cxx" || extension == ".cpp" || extension == ".cu";
}

std::vector<fs::path> productionFiles()
{
  const auto tracking = trackingRoot();
  const auto itsmft = itsmftRoot();
  const std::vector<fs::path> roots{
    tracking / "include",
    tracking / "src",
    itsmft / "ITS/workflow-ca/include",
    itsmft / "ITS/workflow-ca/src",
    itsmft / "MFT/workflow/include",
    itsmft / "MFT/workflow/src",
    itsmft / "common/workflow-combined-ca/include",
    itsmft / "common/workflow-combined-ca/src",
  };

  std::vector<fs::path> files;
  for (const auto& root : roots) {
    BOOST_REQUIRE_MESSAGE(fs::is_directory(root), "cannot inspect " << root.string());
    for (const auto& entry : fs::recursive_directory_iterator(root)) {
      if (entry.is_regular_file() && isSourceFile(entry.path())) {
        files.push_back(entry.path());
      }
    }
  }
  return files;
}

std::string withoutComments(std::string_view source)
{
  std::string result;
  result.reserve(source.size());
  bool blockComment = false;
  for (std::size_t index = 0; index < source.size();) {
    if (blockComment) {
      if (index + 1 < source.size() && source[index] == '*' && source[index + 1] == '/') {
        blockComment = false;
        result += "  ";
        index += 2;
      } else {
        result += source[index++] == '\n' ? '\n' : ' ';
      }
      continue;
    }
    if (index + 1 < source.size() && source[index] == '/' && source[index + 1] == '*') {
      blockComment = true;
      result += "  ";
      index += 2;
    } else if (index + 1 < source.size() && source[index] == '/' && source[index + 1] == '/') {
      while (index < source.size() && source[index] != '\n') {
        result += ' ';
        ++index;
      }
    } else {
      result += source[index++];
    }
  }
  return result;
}

std::size_t directHeaderCount(const fs::path& directory)
{
  BOOST_REQUIRE_MESSAGE(fs::is_directory(directory), "cannot inspect " << directory.string());
  std::size_t count = 0;
  for (const auto& entry : fs::directory_iterator(directory)) {
    count += entry.is_regular_file() && entry.path().extension() == ".h";
  }
  return count;
}

} // namespace

BOOST_AUTO_TEST_CASE(FinalHeaderInventoryIsExact)
{
  const auto include = trackingRoot() / "include/ITSMFTTracking";
  BOOST_CHECK_EQUAL(directHeaderCount(include), 31U);
  BOOST_CHECK_EQUAL(directHeaderCount(include / "detail"), 13U);
  BOOST_CHECK(fs::is_regular_file(include / "TripletFitting.h"));

  const auto surfaceLayout = include / "SurfaceLayout.h";
  BOOST_REQUIRE_MESSAGE(fs::is_regular_file(surfaceLayout), "SurfaceLayout must remain a direct public configuration header");
  const auto layoutText = withoutComments(readFile(surfaceLayout));
  BOOST_CHECK(layoutText.find("class SurfaceLayout") != std::string::npos);
  BOOST_CHECK(layoutText.find("SurfaceLayoutDefinition") != std::string::npos);
  const auto retiredGraph = std::string{"Surface"} + "Graph";
  BOOST_CHECK(layoutText.find(retiredGraph) == std::string::npos);

  const auto traversalTopology = include / "TraversalTopology.h";
  BOOST_REQUIRE_MESSAGE(fs::is_regular_file(traversalTopology), "TraversalTopology must remain a direct topology header");
  const auto topologyText = withoutComments(readFile(traversalTopology));
  BOOST_CHECK(topologyText.find("struct CellPath") != std::string::npos);
  BOOST_CHECK(topologyText.find("struct TraversalTopologyView") != std::string::npos);
  const auto retiredCell = std::string{"Surface"} + "CellTopology";
  BOOST_CHECK(topologyText.find(retiredCell) == std::string::npos);
}

BOOST_AUTO_TEST_CASE(ObsoleteTopologyModelIsAbsentFromLiveSource)
{
  const auto root = trackingRoot();
  const std::vector<std::string> retiredFiles{
    std::string{"Surface"} + "Graph.h",
    std::string{"Surface"} + "GraphBuilder.h",
    std::string{"Surface"} + "GraphBuilder.cxx"};
  for (const auto& relative : retiredFiles) {
    BOOST_CHECK_MESSAGE(!fs::exists(root / "include/ITSMFTTracking" / relative) &&
                          !fs::exists(root / "src" / relative),
                        "obsolete topology file remains: " << relative);
  }

  const std::vector<std::string> retiredSpellings{
    std::string{"Surface"} + "Graph",
    std::string{"Surface"} + "GraphBuilder",
    std::string{"Surface"} + "GraphView",
    std::string{"Surface"} + "GraphDefinition",
    std::string{"Surface"} + "GraphBuildResult",
    std::string{"Surface"} + "GraphBatchResult",
    std::string{"Surface"} + "CellTopology",
    std::string{"Cell"} + "TopologyId"};
  for (const auto& directory : {root / "include", root / "src", root / "test"}) {
    for (const auto& entry : fs::recursive_directory_iterator(directory)) {
      if (!entry.is_regular_file() || !isSourceFile(entry.path())) {
        continue;
      }
      const auto source = readFile(entry.path());
      for (const auto& spelling : retiredSpellings) {
        BOOST_CHECK_MESSAGE(source.find(spelling) == std::string::npos,
                            entry.path().string() << " retains obsolete topology spelling " << spelling);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(RetiredAndRelocatedPublicPathsAreAbsent)
{
  const auto include = trackingRoot() / "include/ITSMFTTracking";
  const std::vector<std::string_view> retiredPublic{
    "SurfaceMeasurementIndex.h",
    "SurfaceId.h",
    "LayerMask.h",
    "SurfaceLinearizationReference.h",
    "RefitLegAssembly.h",
    "BarrelSurfaceStateOperations.h",
    "ForwardSurfaceStateOperations.h",
    "SeedAnchor.h",
    "SurfaceCatalogView.h",
    "StaticSurfaceDescriptor.h",
    "ITSSurfaceSpec.h",
    "MFTSurfaceSpec.h",
    "ClockTimingPublicationView.h",
    "ROFTimingUniformity.h",
    "ClusterDecoder.h",
    "DecodedCluster.h",
    "SurfaceMeasurementAdapters.h",
    "ClusterSource.h",
    "MultiSourceLoading.h",
    "TimeFrameLoadFailure.h",
    "CommonTrackShadow.h",
    "GenericTrackShadow.h",
    "SurfaceTrackingScratch.h",
    "DetectorPublicationAdapter.h",
    "DetectorTrackingOperationAdapterSupport.h",
    "ITSSharedClusterCompatibility.h",
    "MFTPublicationCompatibility.h",
    "MFTFwdTrackHelpers.h",
    "SurfaceKinematicStateLegacyAdapters.h",
    "ConfigKeyValuesPreflight.h",
    "SurfaceKind.h",
    "TrackingOperationAdapter.h",
    "CommonTrack.h",
    "CommonTrackOutputAdapter.h",
    "NativeRefitDriver.h",
  };
  for (const auto header : retiredPublic) {
    BOOST_CHECK_MESSAGE(!fs::exists(include / header), "retired public path remains: " << header);
  }

  const std::string policy = std::string{"Edge"} + "Policy";
  for (const auto suffix : {".h", "Binding.h", "Dispatch.h", "Operations.h", "State.h"}) {
    BOOST_CHECK_MESSAGE(!fs::exists(include / "detail" / (policy + suffix)),
                        "retired policy header remains: " << policy + suffix);
  }

  const std::vector<std::string> forbiddenSpellings{
    policy,
    std::string{"edge"} + "Policy",
    std::string{"For"} + "Policy",
    std::string{"bindLink"} + "PolicyParams",
    "PropagationModel",
  };
  for (const auto& path : productionFiles()) {
    const auto text = readFile(path);
    for (const auto& spelling : forbiddenSpellings) {
      BOOST_CHECK_MESSAGE(text.find(spelling) == std::string::npos,
                          path.string() << " retains retired policy spelling " << spelling);
    }
    for (const auto header : retiredPublic) {
      const auto oldInclude = std::string{"ITSMFTTracking/"} + std::string{header};
      BOOST_CHECK_MESSAGE(text.find(oldInclude) == std::string::npos,
                          path.string() << " retains retired include path " << oldInclude);
    }
  }
}

BOOST_AUTO_TEST_CASE(NoRetiredShadowPublicationPathOrPublisherExists)
{
  for (const auto& path : productionFiles()) {
    const auto text = readFile(path);
    BOOST_CHECK_MESSAGE(text.find("CommonTrackShadow") == std::string::npos,
                        path.string() << " retains the retired CommonTrack shadow identifier");
    BOOST_CHECK_MESSAGE(text.find("GenericTrackShadow") == std::string::npos,
                        path.string() << " recreates a GenericTrack shadow identifier");
    BOOST_CHECK_MESSAGE(text.find("ShadowPublisher") == std::string::npos,
                        path.string() << " recreates a shadow publication layer");
  }
}

BOOST_AUTO_TEST_CASE(GenericTrackAndRefitDriverUseOnlyCurrentProductNames)
{
  const auto include = trackingRoot() / "include/ITSMFTTracking";
  BOOST_CHECK(fs::is_regular_file(include / "GenericTrack.h"));
  BOOST_CHECK(fs::is_regular_file(include / "GenericTrackOutputAdapter.h"));
  BOOST_CHECK(fs::is_regular_file(include / "RefitDriver.h"));

  for (const auto& path : productionFiles()) {
    const auto text = readFile(path);
    BOOST_CHECK_MESSAGE(text.find("ITSMFTTracking/NativeRefitDriver.h") == std::string::npos,
                        path.string() << " retains the superseded RefitDriver include");
    BOOST_CHECK_MESSAGE(text.find("ITSMFTTracking/CommonTrack.h") == std::string::npos,
                        path.string() << " retains the superseded GenericTrack model include");
    BOOST_CHECK_MESSAGE(text.find("ITSMFTTracking/CommonTrackOutputAdapter.h") == std::string::npos,
                        path.string() << " retains the superseded GenericTrack adapter include");
  }
}

BOOST_AUTO_TEST_CASE(WorkflowsUseOnlyTheDetailAdapterFacades)
{
  const auto itsmft = itsmftRoot();
  const std::vector<fs::path> workflowRoots{
    itsmft / "ITS/workflow-ca",
    itsmft / "MFT/workflow",
    itsmft / "common/workflow-combined-ca",
  };
  const std::vector<std::string_view> forbiddenDetailIncludes{
    "ITSMFTTracking/detail/CommonTrackShadow.h",
    "ITSMFTTracking/detail/GenericTrackShadow.h",
    "ITSMFTTracking/detail/SurfaceTrackingScratch.h",
    "ITSMFTTracking/detail/ITSSharedClusterCompatibility.h",
    "ITSMFTTracking/detail/MFTPublicationCompatibility.h",
    "ITSMFTTracking/detail/MFTFwdTrackHelpers.h",
    "ITSMFTTracking/detail/SurfaceKinematicStateLegacyAdapters.h",
    "ITSMFTTracking/detail/TrackingKernelParameters.h",
    "ITSMFTTracking/detail/CandidateFinding.h",
  };

  for (const auto& root : workflowRoots) {
    for (const auto& entry : fs::recursive_directory_iterator(root)) {
      if (!entry.is_regular_file() || !isSourceFile(entry.path())) {
        continue;
      }
      const auto text = readFile(entry.path());
      for (const auto include : forbiddenDetailIncludes) {
        BOOST_CHECK_MESSAGE(text.find(include) == std::string::npos,
                            entry.path().string() << " bypasses the workflow detail facade with " << include);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(GenericPlanTypesHaveNoDetectorOrLayerCountAuthority)
{
  const auto root = trackingRoot();
  const std::vector<std::string_view> genericPlanFiles{
    "include/ITSMFTTracking/SurfaceDescriptor.h",
    "include/ITSMFTTracking/IdTypes.h",
    "include/ITSMFTTracking/SurfaceKinematicState.h",
    "include/ITSMFTTracking/SurfaceMeasurement.h",
    "include/ITSMFTTracking/SurfaceSpec.h",
    "include/ITSMFTTracking/detail/TrackingKernelParameters.h",
  };
  const std::vector<std::string_view> forbiddenAuthorities{
    "o2::detectors::DetID",
    "DetID::ITS",
    "DetID::MFT",
    "ITSNLayers",
    "MFTNLayers",
    "LayersNumber",
    "LayerZCoordinate",
    "o2::its::",
    "o2::mft::",
    "#include \"ITStracking/",
    "#include \"MFTTracking/",
    "NLayers",
  };

  for (const auto relative : genericPlanFiles) {
    const auto path = root / relative;
    const auto code = withoutComments(readFile(path));
    for (const auto authority : forbiddenAuthorities) {
      BOOST_CHECK_MESSAGE(code.find(authority) == std::string::npos,
                          relative << " acquires detector/layer-count authority through " << authority);
    }
  }
}

BOOST_AUTO_TEST_CASE(TraversalWorkspaceOwnsPassTopologyAndTrackerBuildsIt)
{
  const auto root = trackingRoot();
  const auto timeFrame = withoutComments(readFile(root / "include/ITSMFTTracking/TimeFrame.h"));
  const auto workspace = withoutComments(readFile(root / "include/ITSMFTTracking/detail/SurfaceTrackingScratch.h"));
  const auto tracker = withoutComments(readFile(root / "include/ITSMFTTracking/Tracker.h"));
  const auto trackerSource = withoutComments(readFile(root / "src/Tracker.cxx"));

  BOOST_CHECK(timeFrame.find("SurfacePlanBinding") == std::string::npos);
  BOOST_CHECK(timeFrame.find("mBindings") == std::string::npos);
  BOOST_CHECK(workspace.find("std::vector<EdgeId> edges") != std::string::npos);
  BOOST_CHECK(workspace.find("std::vector<CellPathId> cells") != std::string::npos);
  BOOST_CHECK(workspace.find("roadStartCells") != std::string::npos);
  BOOST_CHECK(tracker.find("buildTraversalPlan") != std::string::npos);
  BOOST_CHECK(trackerSource.find("buildTraversalPlan(context.workspace") != std::string::npos);
}
