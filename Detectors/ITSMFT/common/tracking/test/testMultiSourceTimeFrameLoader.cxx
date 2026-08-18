// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.

#define BOOST_TEST_MODULE ITSMFT TimeFrame source loading
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
#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/SurfaceLayout.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/TimeFrame.h"
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

  o2::itsmft::tracking::ClusterDecodeResult decode(
    const CompClusterExt& cluster, BoundedPatternCursor& patterns, const TopologyDictionary* dictionary,
    uint32_t, bool) const override
  {
    const auto decodedPattern = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dictionary);
    o2::itsmft::tracking::ClusterDecodeResult result;
    if (!decodedPattern.ok()) {
      result.error = decodedPattern.error;
      return result;
    }
    const auto layer = static_cast<int>(cluster.getSensorID());
    auto& decoded = result.decoded;
    decoded.global = {static_cast<float>(layer), static_cast<float>(cluster.getRow()), static_cast<float>(cluster.getCol())};
    decoded.cylinderFrame = {10.f + layer, 1.f, 2.f, 0.1f};
    decoded.rowColumnCovariance = {decodedPattern.sig2Row, 0.f, decodedPattern.sig2Col};
    decoded.shape = decodedPattern.shape;
    decoded.layer = layer;
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
  std::array<LayerId, 2> layerToSurface;
  FakeClusterDecoder decoder;

  explicit OneClusterSource(LayerId surface) : layerToSurface{surface, LayerId{static_cast<uint16_t>(surface.value() + 1)}}, decoder{o2::detectors::DetID::ITS} {}

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
  SurfaceLayout layout;
  std::vector<SurfaceDescriptor> catalog;
};

ThreeSourceConfiguration makeConfiguration()
{
  std::vector<SurfaceDescriptor> catalog;
  for (uint16_t id = 0; id < 6; ++id) {
    catalog.push_back(SurfaceDescriptor{LayerId{id}, static_cast<uint8_t>(id % 2), static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  }
  std::vector<LayerId> ordered;
  for (uint16_t id = 0; id < 6; ++id) {
    ordered.push_back(LayerId{id});
  }
  return {SurfaceLayout{gsl::span<const SurfaceDescriptor>{catalog.data(), catalog.size()}, makeSurfaceLayoutChain(ordered)}, std::move(catalog)};
}

void configureFrame(TimeFrame& frame, const SurfaceLayout& layout)
{
  std::vector<LayerId> ordered;
  for (uint16_t id = 0; id < 3; ++id) {
    for (const auto surface : {LayerId{static_cast<uint16_t>(id * 2)}, LayerId{static_cast<uint16_t>(id * 2 + 1)}}) {
      ordered.push_back(surface);
    }
  }
  auto configuredLayout = layout;
  BOOST_REQUIRE(frame.configure(std::move(configuredLayout), 0, 0,
                                std::make_shared<BoundedMemoryResource>()));
}

std::vector<ClusterSourceInput> makeSources(const std::array<OneClusterSource, 3>& sources)
{
  return {sources[0].input(ClusterSourceId{0}), sources[1].input(ClusterSourceId{1}), sources[2].input(ClusterSourceId{2})};
}
} // namespace

BOOST_AUTO_TEST_CASE(DirectThreeSourceLoadInstallsAllSources)
{
  auto configuration = makeConfiguration();
  TimeFrame frame;
  configureFrame(frame, configuration.layout);
  const std::array<OneClusterSource, 3> inputs{OneClusterSource{LayerId{0}}, OneClusterSource{LayerId{2}}, OneClusterSource{LayerId{4}}};
  const auto sources = makeSources(inputs);

  const auto result = loadTimeFrameSources(frame, sources, configuration.layout.getSurfaceCatalog(), {50, 5});
  BOOST_REQUIRE_MESSAGE(result.ok(), "load error=" << static_cast<int>(result.error) << " source=" << result.source.value());
  for (uint16_t id = 0; id < 3; ++id) {
    BOOST_CHECK_EQUAL(frame.getGlobalMeasurements(LayerId{static_cast<uint16_t>(id * 2)}).size(), 1u);
  }
  BOOST_CHECK_EQUAL(frame.getTotalClusters(), 3);
}

BOOST_AUTO_TEST_CASE(FailedSourcePartitionClearsTimeFrameAndRetrySucceeds)
{
  auto configuration = makeConfiguration();
  TimeFrame frame;
  configureFrame(frame, configuration.layout);
  const std::array<OneClusterSource, 3> inputs{OneClusterSource{LayerId{0}}, OneClusterSource{LayerId{2}}, OneClusterSource{LayerId{4}}};
  auto sources = makeSources(inputs);
  const auto baseline = loadTimeFrameSources(frame, sources, configuration.layout.getSurfaceCatalog(), {50, 5});
  BOOST_REQUIRE_MESSAGE(baseline.ok(), "baseline load error=" << static_cast<int>(baseline.error) << " source=" << baseline.source.value());
  auto malformedInputs = inputs;
  malformedInputs[2].layerToSurface[0] = LayerId{0};
  const auto malformedSources = makeSources(malformedInputs);
  const auto failed = loadTimeFrameSources(frame, malformedSources, configuration.layout.getSurfaceCatalog(), {50, 5});
  BOOST_CHECK(failed.error == MultiSourceLoadError::InvalidLayerMapping);
  BOOST_CHECK_EQUAL(frame.getTotalClusters(), 0);

  BOOST_REQUIRE(loadTimeFrameSources(frame, sources, configuration.layout.getSurfaceCatalog(), {50, 5}).ok());
  BOOST_CHECK_EQUAL(frame.getTotalClusters(), 3);
}

BOOST_AUTO_TEST_CASE(UnconfiguredFrameAndSourceQualificationFailBeforeCommit)
{
  auto configuration = makeConfiguration();
  TimeFrame frame;
  const std::array<OneClusterSource, 3> inputs{OneClusterSource{LayerId{0}}, OneClusterSource{LayerId{2}}, OneClusterSource{LayerId{4}}};
  const auto sources = makeSources(inputs);
  BOOST_CHECK(loadTimeFrameSources(frame, sources, configuration.layout.getSurfaceCatalog(), {50, 5}).error ==
              MultiSourceLoadError::FrameNotConfigured);

  configureFrame(frame, configuration.layout);
  auto wrong = sources;
  wrong[1].id = ClusterSourceId{5};
  const auto result = loadTimeFrameSources(frame, wrong, configuration.layout.getSurfaceCatalog(), {50, 5});
  BOOST_CHECK(result.error == MultiSourceLoadError::NonDenseSourceIds);
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
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
