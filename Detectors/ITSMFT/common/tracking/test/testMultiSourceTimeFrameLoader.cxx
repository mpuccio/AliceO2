// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.

#define BOOST_TEST_MODULE ITSMFT MultiSourceTimeFrameLoader
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{
class FakeClusterDecoder final : public ClusterDecoder
{
 public:
  explicit FakeClusterDecoder(o2::detectors::DetID::ID detector) : mDetector(detector) {}

  o2::itsmft::tracking::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster, BoundedPatternCursor& patterns, const TopologyDictionary* dictionary,
    gsl::span<const SurfaceId> layerToSurface, ClusterSourceId source, uint32_t externalIndex,
    uint32_t sourceROF, bool) const override
  {
    const auto decodedPattern = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dictionary);
    o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
    if (!decodedPattern.ok()) {
      result.error = decodedPattern.error;
      return result;
    }
    const auto layer = static_cast<int>(cluster.getSensorID());
    if (layer < 0 || static_cast<size_t>(layer) >= layerToSurface.size()) {
      return result;
    }
    result.layer = layer;
    result.layerMapped = true;
    result.kind = SurfaceKind::Cylinder;
    DecodedCluster decoded{};
    decoded.global = {static_cast<float>(layer), static_cast<float>(cluster.getRow()), static_cast<float>(cluster.getCol())};
    decoded.cylinderFrame = {10.f + layer, 1.f, 2.f, 0.1f};
    decoded.rowColumnCovariance = {decodedPattern.sig2Row, 0.f, decodedPattern.sig2Col};
    decoded.shape = decodedPattern.shape;
    decoded.sensor = static_cast<uint32_t>(layer);
    decoded.layer = layer;
    result = makeCylinderMeasurementDecodeResult(
      decoded, DetectorSensorId{static_cast<uint32_t>(mDetector), decoded.sensor}, layerToSurface[layer],
      ClusterRef{source, externalIndex}, sourceROF);
    return result;
  }

 private:
  o2::detectors::DetID::ID mDetector;
};

constexpr std::array<unsigned char, 3> onePixelPattern{1, 1, 0x80};

const TopologyDictionary& dictionary()
{
  static const TopologyDictionary value;
  return value;
}

struct OneClusterSource {
  std::vector<CompClusterExt> clusters{{0, 1, CompCluster::InvalidPatternID, 0}};
  std::vector<unsigned char> patterns{onePixelPattern.begin(), onePixelPattern.end()};
  std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, 1}};
  std::array<SurfaceId, 2> layerToSurface;
  FakeClusterDecoder decoder;

  explicit OneClusterSource(SurfaceId surface) : layerToSurface{surface, SurfaceId{static_cast<uint16_t>(surface.value() + 1)}}, decoder{o2::detectors::DetID::ITS} {}

  ClusterSourceInput input(ClusterSourceId source) const
  {
    ClusterSourceInput value{};
    value.id = source;
    value.detector = o2::detectors::DetID::ITS;
    value.clusters = clusters;
    value.patterns = patterns;
    value.rofs = rofs;
    value.dictionary = &dictionary();
    value.layerToSurface = layerToSurface;
    value.timing = ROFTimingConfig{40, 0, 0, 0};
    value.decoder = &decoder;
    return value;
  }
};

struct ThreeSourceConfiguration {
  SurfaceGraph graph;
  std::vector<SurfaceDescriptor> catalog;
};

ThreeSourceConfiguration makeConfiguration()
{
  std::vector<SurfaceDescriptor> catalog;
  for (uint16_t id = 0; id < 6; ++id) {
    catalog.push_back(SurfaceDescriptor{SurfaceId{id}, static_cast<uint8_t>(id % 2), static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  }
  SurfaceGraph topology{6};
  topology.finalize();
  return {SurfaceGraph{catalog, std::move(topology)}, std::move(catalog)};
}

void configureFrame(TimeFrame& frame, const SurfaceGraph& graph)
{
  const auto view = graph.getView();
  std::vector<SurfaceId> ordered;
  SurfaceMask owned;
  for (uint16_t id = 0; id < 3; ++id) {
    for (const auto surface : {SurfaceId{static_cast<uint16_t>(id * 2)}, SurfaceId{static_cast<uint16_t>(id * 2 + 1)}}) {
      owned.set(surface);
      ordered.push_back(surface);
    }
  }
  auto built = SurfacePlanBinding::build(view, owned, ordered);
  BOOST_REQUIRE_MESSAGE(built.ok(), "binding error=" << static_cast<int>(built.error));
  std::vector<SurfaceGraph> graphs;
  graphs.push_back(graph);
  std::vector<TrackingParameters> parameters(1);
  std::vector<std::unique_ptr<SurfacePlanBinding>> bindings;
  bindings.push_back(std::move(built.binding));
  std::vector<TrackingWorkspaceCapacity> capacities{{ordered.size(), 0, 0}};
  BOOST_REQUIRE(frame.commitConfiguration(std::move(graphs), std::move(parameters), std::move(bindings),
                                          std::move(capacities), std::make_shared<BoundedMemoryResource>()));
}

std::vector<ClusterSourceInput> makeSources(const std::array<OneClusterSource, 3>& sources)
{
  return {sources[0].input(ClusterSourceId{0}), sources[1].input(ClusterSourceId{1}), sources[2].input(ClusterSourceId{2})};
}
} // namespace

BOOST_AUTO_TEST_CASE(DirectThreeSourceTransactionInstallsAllSources)
{
  auto configuration = makeConfiguration();
  TimeFrame frame;
  configureFrame(frame, configuration.graph);
  const std::array<OneClusterSource, 3> inputs{OneClusterSource{SurfaceId{0}}, OneClusterSource{SurfaceId{2}}, OneClusterSource{SurfaceId{4}}};
  const auto sources = makeSources(inputs);

  const auto result = MultiSourceTimeFrameLoader::load(frame, sources, configuration.graph.getView().getSurfaceCatalogView(), {50, 5});
  BOOST_REQUIRE_MESSAGE(result.ok(), "load error=" << static_cast<int>(result.error) << " source=" << result.source.value());
  BOOST_CHECK_EQUAL(frame.getEventResetCount(), 1u);
  for (uint16_t id = 0; id < 3; ++id) {
    BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{static_cast<uint16_t>(id * 2)}).size(), 1u);
  }
  BOOST_CHECK_EQUAL(frame.getWorkspace().getTotalClusters(), 3);
}

BOOST_AUTO_TEST_CASE(FailedSourcePartitionLeavesPriorEventAndRetrySucceeds)
{
  auto configuration = makeConfiguration();
  TimeFrame frame;
  configureFrame(frame, configuration.graph);
  const std::array<OneClusterSource, 3> inputs{OneClusterSource{SurfaceId{0}}, OneClusterSource{SurfaceId{2}}, OneClusterSource{SurfaceId{4}}};
  auto sources = makeSources(inputs);
  const auto baseline = MultiSourceTimeFrameLoader::load(frame, sources, configuration.graph.getView().getSurfaceCatalogView(), {50, 5});
  BOOST_REQUIRE_MESSAGE(baseline.ok(), "baseline load error=" << static_cast<int>(baseline.error) << " source=" << baseline.source.value());
  const auto resetCount = frame.getEventResetCount();

  auto malformedInputs = inputs;
  malformedInputs[2].layerToSurface[0] = SurfaceId{0};
  const auto malformedSources = makeSources(malformedInputs);
  const auto failed = MultiSourceTimeFrameLoader::load(frame, malformedSources, configuration.graph.getView().getSurfaceCatalogView(), {50, 5});
  BOOST_CHECK(failed.error == MultiSourceLoadError::InvalidLayerMapping);
  BOOST_CHECK_EQUAL(frame.getEventResetCount(), resetCount);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{4}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getWorkspace().getTotalClusters(), 3);

  BOOST_REQUIRE(MultiSourceTimeFrameLoader::load(frame, sources, configuration.graph.getView().getSurfaceCatalogView(), {50, 5}).ok());
  BOOST_CHECK_EQUAL(frame.getEventResetCount(), resetCount + 1);
}

BOOST_AUTO_TEST_CASE(UnconfiguredFrameAndSourceQualificationFailBeforeCommit)
{
  auto configuration = makeConfiguration();
  TimeFrame frame;
  const std::array<OneClusterSource, 3> inputs{OneClusterSource{SurfaceId{0}}, OneClusterSource{SurfaceId{2}}, OneClusterSource{SurfaceId{4}}};
  const auto sources = makeSources(inputs);
  BOOST_CHECK(MultiSourceTimeFrameLoader::load(frame, sources, configuration.graph.getView().getSurfaceCatalogView(), {50, 5}).error ==
              MultiSourceLoadError::FrameNotConfigured);

  configureFrame(frame, configuration.graph);
  auto wrong = sources;
  wrong[1].id = ClusterSourceId{5};
  const auto result = MultiSourceTimeFrameLoader::load(frame, wrong, configuration.graph.getView().getSurfaceCatalogView(), {50, 5});
  BOOST_CHECK(result.error == MultiSourceLoadError::NonDenseSourceIds);
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
  BOOST_CHECK_EQUAL(frame.getEventResetCount(), 0u);
}

namespace
{
std::vector<std::string> readLines(const std::string& path)
{
  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path);
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(input, line)) {
    lines.push_back(line);
  }
  return lines;
}
} // namespace

BOOST_AUTO_TEST_CASE(DirectLoaderHierarchyGuard)
{
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const std::string root = testDirectory + "/..";
  const std::vector<std::string> forbidden = {
    "Load"
    "Target",
    "load"
    "Target(",
    "Atomic"
    "LoadBinding",
    "m"
    "Load"
    "Target",
    "friend class Multi"
    "SourceTimeFrameLoader"};
  for (const auto& directory : {root + "/include", root + "/src", root + "/test"}) {
    for (const auto& entry : std::filesystem::recursive_directory_iterator{directory}) {
      if (!entry.is_regular_file()) {
        continue;
      }
      const auto extension = entry.path().extension();
      if (extension != ".h" && extension != ".cxx" && extension != ".cmake") {
        continue;
      }
      const auto file = entry.path().string();
      for (const auto& line : readLines(file)) {
        for (const auto& token : forbidden) {
          BOOST_CHECK_MESSAGE(line.find(token) == std::string::npos,
                              file << " contains deleted L6 token '" << token << "': " << line);
        }
      }
    }
  }
}
