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
  return {SurfaceDescriptor{SurfaceId{0}, 0, 0, SurfaceKind::Cylinder},
          SurfaceDescriptor{SurfaceId{1}, 1, 0, SurfaceKind::Cylinder},
          SurfaceDescriptor{SurfaceId{2}, 2, 0, SurfaceKind::Cylinder}};
}

TrackerInitialization makeConfiguration(const std::vector<SurfaceDescriptor>& catalog,
                                        std::shared_ptr<BoundedMemoryResource> pool)
{
  const std::vector<SurfaceId> ordered{SurfaceId{2}, SurfaceId{0}, SurfaceId{1}};
  TrackerInitialization configuration;
  configuration.catalog = SurfaceCatalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  configuration.memoryPool = std::move(pool);
  for (const auto [holes, seeds] : {std::pair{uint32_t{1u << 1}, uint32_t{1u << 2}},
                                    std::pair{uint32_t{1u << 2}, uint32_t{1u << 0}}}) {
    TrackerIterationConfiguration iteration;
    TrackingParameters parameters;
    parameters.NLayers = 0;
    iteration.layout = makeSurfaceLayoutChain(ordered, 1, SurfaceMask{holes}, SurfaceMask{seeds});
    iteration.parameters = parameters;
    configuration.iterations.push_back(std::move(iteration));
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

BOOST_AUTO_TEST_CASE(ConfigurationCommitIsAtomic)
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
  BOOST_REQUIRE_EQUAL(frame.getNIterations(), 2u);
  BOOST_CHECK(frame.getLayout(0).getOrderedSurfaces()[0] == SurfaceId{2});
  BOOST_CHECK(frame.getLayout(1).getOrderedSurfaces()[0] == SurfaceId{2});
  BOOST_CHECK(frame.getLayout(0).getSeedingSurfaces().has(SurfaceId{2}));
  BOOST_CHECK(frame.getLayout(1).getSeedingSurfaces().has(SurfaceId{0}));
  BOOST_CHECK_EQUAL(frame.getWorkspace().getNTraversalWorkspaces(), 2u);
  BOOST_CHECK(!frame.getWorkspace().getTraversalWorkspace(0).valid);
  BOOST_CHECK(!frame.getWorkspace().getTraversalWorkspace(1).valid);
  const auto* oldLayout = &frame.getLayout(0);
  const auto* oldPool = frame.getMemoryPool().get();

  TrackerInitialization invalid;
  invalid.catalog = configuration.catalog;
  invalid.memoryPool = std::make_shared<BoundedMemoryResource>();
  invalid.iterations.push_back(TrackerIterationConfiguration{});
  BOOST_CHECK(!tracker.initialize(frame, invalid).ok());
  BOOST_CHECK(frame.isConfigured());
  BOOST_CHECK(&frame.getLayout(0) == oldLayout);
  BOOST_CHECK(frame.getMemoryPool().get() == oldPool);

  auto replacement = makeConfiguration(catalog, std::make_shared<BoundedMemoryResource>());
  BOOST_REQUIRE(tracker.initialize(frame, replacement).ok());
  BOOST_CHECK(&frame.getLayout(0) != oldLayout);
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
  const auto* layout = &frame.getLayout(0);
  const auto capacity = *frame.getWorkspaceCapacity(0);
  BOOST_CHECK_EQUAL(frame.getWorkspace().getNTraversalWorkspaces(), 2u);
  frame.getGenericTracks().push_back(GenericTrack{});
  frame.resetEvent();
  BOOST_CHECK(frame.isConfigured());
  BOOST_CHECK(&frame.getLayout(0) == layout);
  BOOST_CHECK_EQUAL(frame.getWorkspaceCapacity(0)->cells, capacity.cells);
  BOOST_CHECK(frame.getGenericTracks().empty());
}

BOOST_AUTO_TEST_CASE(ConfigurationOwnershipGuard)
{
  const auto trackingRoot = std::filesystem::path{__FILE__}.parent_path().parent_path();
  const std::vector<std::filesystem::path> sources{
    trackingRoot / "include/ITSMFTTracking/Tracker.h",
    trackingRoot / "include/ITSMFTTracking/TrackerTraits.h",
    trackingRoot / "src/Tracker.cxx",
    trackingRoot / "src/TrackerTraits.cxx"};
  for (const auto& source : sources) {
    BOOST_REQUIRE_MESSAGE(std::filesystem::exists(source), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "mGraphs;"), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "mTrackParams;"), source.string());
    if (source.filename() != "TrackerTraits.h") {
      BOOST_CHECK_MESSAGE(!contains(source, "mMemoryPool;"), source.string());
    }
  }
}
