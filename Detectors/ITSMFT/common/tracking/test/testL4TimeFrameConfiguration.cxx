// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details.

#define BOOST_TEST_MODULE ITSMFT L4 TimeFrame configuration
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/Tracker.h"

using namespace o2::itsmft::tracking;
using o2::itsmft::TrackingParameters;

namespace
{
std::vector<SurfaceDescriptor> makeCatalog()
{
  std::vector<SurfaceDescriptor> catalog;
  for (uint16_t layer = 0; layer < 3; ++layer) {
    SurfaceDescriptor descriptor{layer, 0, SurfaceKind::Cylinder};
    descriptor.referenceCoordinate = static_cast<float>(layer + 1);
    descriptor.chartRange = {-10.f, 10.f};
    catalog.push_back(descriptor);
  }
  return catalog;
}

TrackerInitialization makeConfiguration(const std::vector<SurfaceDescriptor>& catalog,
                                        std::shared_ptr<BoundedMemoryResource> pool)
{
  TrackerInitialization configuration;
  configuration.catalog = SurfaceCatalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  configuration.memoryPool = std::move(pool);
  configuration.layout = makeDetectorLayout(LayerMask{1u << 1});
  for (const auto seeds : {uint32_t{1u << 2}, uint32_t{1u << 0}}) {
    TrackingParameters parameters;
    parameters.NLayers = 3;
    parameters.LayerxX0.assign(3, 0.f);
    parameters.LayerRadii = {1.f, 2.f, 3.f};
    parameters.LayerResolution.assign(3, 0.001f);
    parameters.SystError2Row.assign(3, 0.f);
    parameters.SystError2Col.assign(3, 0.f);
    parameters.MaxHoles = 1;
    parameters.StartLayerMask = seeds;
    configuration.parameters.push_back(parameters);
  }
  return configuration;
}

bool contains(const std::filesystem::path& path, const std::string& text)
{
  std::ifstream input(path);
  const std::string content{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
  return content.find(text) != std::string::npos;
}
} // namespace

BOOST_AUTO_TEST_CASE(ConfigurationIsInstalledOnce)
{
  const auto catalog = makeCatalog();
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  Tracker tracker;

  auto configuration = makeConfiguration(catalog, pool);
  const auto initialResult = tracker.initialize(frame, configuration);
  BOOST_REQUIRE_MESSAGE(initialResult.ok(), "initial configuration error=" << static_cast<int>(initialResult.error)
                                                                           << " layout=" << static_cast<int>(initialResult.layoutError));
  BOOST_REQUIRE(frame.isConfigured());
  BOOST_REQUIRE(tracker.isConfiguredFor(frame));
  BOOST_REQUIRE_EQUAL(tracker.getIterationConfigurations().size(), 2u);
  BOOST_CHECK_EQUAL(frame.getLayout()[LayerId{0}].detectorSurfaceIndex, 0u);
  BOOST_CHECK(tracker.getIterationConfigurations()[0].parameters.StartLayerMask.has(2));
  BOOST_CHECK(tracker.getIterationConfigurations()[1].parameters.StartLayerMask.has(0));
  BOOST_CHECK(frame.getGenericTracks().empty());
  const auto* oldLayout = &frame.getLayout();
  const auto* scratch = &frame.getScratch();
  const auto* oldPool = frame.getMemoryPool().get();

  const auto replacementPool = std::make_shared<BoundedMemoryResource>();
  auto replacement = makeConfiguration(catalog, replacementPool);
  const auto replacementResult = tracker.initialize(frame, replacement);
  BOOST_CHECK(!replacementResult.ok());
  BOOST_CHECK(replacementResult.error == TrackerInitializationError::FrameAlreadyConfigured);
  BOOST_CHECK(&frame.getLayout() == oldLayout);
  BOOST_CHECK(&frame.getScratch() == scratch);
  BOOST_CHECK(frame.getMemoryPool().get() == oldPool);
  BOOST_CHECK(tracker.isConfiguredFor(frame));
}

BOOST_AUTO_TEST_CASE(ResetPreservesStaticConfigurationAndCapacity)
{
  const auto catalog = makeCatalog();
  TimeFrame frame;
  Tracker tracker;
  auto configuration = makeConfiguration(catalog, std::make_shared<BoundedMemoryResource>());
  const auto initialResult = tracker.initialize(frame, configuration);
  BOOST_REQUIRE_MESSAGE(initialResult.ok(), "initial configuration error=" << static_cast<int>(initialResult.error)
                                                                           << " layout=" << static_cast<int>(initialResult.layoutError));
  const auto* layout = &frame.getLayout();
  const auto cellCapacity = frame.getScratch().getNCells();
  frame.getGenericTracks().push_back(GenericTrack{});
  frame.resetTimeFrame();
  BOOST_CHECK(frame.isConfigured());
  BOOST_CHECK(&frame.getLayout() == layout);
  BOOST_CHECK_EQUAL(frame.getScratch().getNCells(), cellCapacity);
  BOOST_CHECK(frame.getGenericTracks().empty());
}

BOOST_AUTO_TEST_CASE(ConfigurationOwnershipGuard)
{
  const auto trackingRoot = std::filesystem::path{__FILE__}.parent_path().parent_path();
  const std::vector<std::filesystem::path> sources{
    trackingRoot / "include/ITSMFTTracking/TimeFrame.h",
    trackingRoot / "include/ITSMFTTracking/Tracker.h",
    trackingRoot / "include/ITSMFTTracking/TrackerTraits.h",
    trackingRoot / "src/Tracker.cxx",
    trackingRoot / "src/TrackerTraits.cxx"};
  for (const auto& source : sources) {
    BOOST_REQUIRE_MESSAGE(std::filesystem::exists(source), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "mGraphs;"), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "mTrackParams;"), source.string());
    if (source.filename() != "TrackerTraits.h" && source.filename() != "TimeFrame.h") {
      BOOST_CHECK_MESSAGE(!contains(source, "mMemoryPool;"), source.string());
    }
    if (source.filename() == "TimeFrame.h") {
      BOOST_CHECK_MESSAGE(!contains(source, "unique_ptr<TimeFrameScratch"), source.string());
      BOOST_CHECK_MESSAGE(!contains(source, "publishConfiguration"), source.string());
    }
  }
}
