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
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
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
                                        std::shared_ptr<BoundedMemoryResource> pool,
                                        ClusterSourceId source = ClusterSourceId{3})
{
  const std::vector<SurfaceId> ordered{SurfaceId{2}, SurfaceId{0}, SurfaceId{1}};
  SurfaceMask owned;
  for (const auto surface : ordered) {
    owned.set(surface);
  }

  TrackerInitialization configuration;
  configuration.catalog = SurfaceCatalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  configuration.memoryPool = std::move(pool);
  for (const auto [holes, seeds] : {std::pair{uint32_t{1u << 1}, uint32_t{1u << 2}},
                                    std::pair{uint32_t{1u << 2}, uint32_t{1u << 0}}}) {
    TrackerIterationConfiguration iteration;
    iteration.graphSubgraphs.push_back(SurfaceGraphSubgraph{ordered, 1, SurfaceMask{holes}, SurfaceMask{seeds}});
    TrackingParameters parameters;
    parameters.NLayers = 0;
    iteration.parameters.push_back(parameters);
    iteration.bindings.push_back(SurfacePlanBinding::Declaration{source, owned, ordered,
                                                                 SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder});
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

BOOST_AUTO_TEST_CASE(ConfigurationCommitIsAtomicAndSourceQualified)
{
  const auto catalog = makeCatalog();
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  TrackerTraits traits;
  Tracker tracker{&traits};
  tracker.setSource(ClusterSourceId{3});

  auto configuration = makeConfiguration(catalog, pool);
  const auto initialResult = tracker.initialize(frame, configuration);
  BOOST_REQUIRE_MESSAGE(initialResult.ok(), "initial configuration error=" << static_cast<int>(initialResult.error)
                                                                           << " graph=" << static_cast<int>(initialResult.graphError)
                                                                           << " binding=" << static_cast<int>(initialResult.bindingError));
  BOOST_REQUIRE(frame.isConfigured());
  BOOST_REQUIRE_EQUAL(frame.getNIterations(), 2u);
  BOOST_REQUIRE(frame.getBinding(0, ClusterSourceId{3}) != nullptr);
  BOOST_REQUIRE(frame.getBinding(1, ClusterSourceId{3}) != nullptr);
  BOOST_CHECK(frame.getBinding(0, ClusterSourceId{3})->getOrderedSurfaces()[0] == SurfaceId{2});
  BOOST_CHECK(frame.getBinding(1, ClusterSourceId{3})->getOrderedSurfaces()[0] == SurfaceId{2});
  BOOST_CHECK(frame.getGraph(0).getView().seedingSurfaces.has(SurfaceId{2}));
  BOOST_CHECK(frame.getGraph(1).getView().seedingSurfaces.has(SurfaceId{0}));
  const auto* oldBinding = frame.getBinding(0, ClusterSourceId{3});
  const auto* oldPool = frame.getMemoryPool().get();

  TrackerInitialization invalid;
  invalid.catalog = configuration.catalog;
  invalid.memoryPool = std::make_shared<BoundedMemoryResource>();
  invalid.iterations.push_back(TrackerIterationConfiguration{});
  BOOST_CHECK(!tracker.initialize(frame, invalid).ok());
  BOOST_CHECK(frame.isConfigured());
  BOOST_CHECK(frame.getBinding(0, ClusterSourceId{3}) == oldBinding);
  BOOST_CHECK(frame.getMemoryPool().get() == oldPool);

  auto replacement = makeConfiguration(catalog, std::make_shared<BoundedMemoryResource>(), ClusterSourceId{9});
  tracker.setSource(ClusterSourceId{9});
  BOOST_REQUIRE(tracker.initialize(frame, replacement).ok());
  BOOST_CHECK(frame.getBinding(0, ClusterSourceId{9}) != nullptr);
  BOOST_CHECK(frame.getBinding(0, ClusterSourceId{3}) == nullptr);
}

BOOST_AUTO_TEST_CASE(ResetPreservesStaticConfigurationAndCapacity)
{
  const auto catalog = makeCatalog();
  TimeFrame frame;
  TrackerTraits traits;
  Tracker tracker{&traits};
  tracker.setSource(ClusterSourceId{3});
  auto configuration = makeConfiguration(catalog, std::make_shared<BoundedMemoryResource>());
  const auto initialResult = tracker.initialize(frame, configuration);
  BOOST_REQUIRE_MESSAGE(initialResult.ok(), "initial configuration error=" << static_cast<int>(initialResult.error)
                                                                           << " graph=" << static_cast<int>(initialResult.graphError)
                                                                           << " binding=" << static_cast<int>(initialResult.bindingError));
  const auto* binding = frame.getBinding(0, ClusterSourceId{3});
  const auto capacity = *frame.getWorkspaceCapacity(0, ClusterSourceId{3});
  frame.getCommonTracks().push_back(CommonTrack{});
  frame.wipe();
  BOOST_CHECK(frame.isConfigured());
  BOOST_CHECK(frame.getBinding(0, ClusterSourceId{3}) == binding);
  BOOST_CHECK_EQUAL(frame.getWorkspaceCapacity(0, ClusterSourceId{3})->cells, capacity.cells);
  BOOST_CHECK(frame.getCommonTracks().empty());
}

BOOST_AUTO_TEST_CASE(ConfigurationOwnershipGuard)
{
  const auto trackingRoot = std::filesystem::path{__FILE__}.parent_path().parent_path();
  const std::vector<std::filesystem::path> sources{
    trackingRoot / "include/ITSMFTTracking/Tracker.h",
    trackingRoot / "include/ITSMFTTracking/TrackingInterface.h",
    trackingRoot / "include/ITSMFTTracking/SurfacePlanTrackingParticipant.h",
    trackingRoot / "src/Tracker.cxx",
    trackingRoot / "src/TrackingInterface.cxx",
    trackingRoot / "src/SurfacePlanTrackingParticipant.cxx"};
  for (const auto& source : sources) {
    BOOST_REQUIRE_MESSAGE(std::filesystem::exists(source), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "mGraphs;"), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "mTrackParams;"), source.string());
    if (source.filename() != "TrackerTraits.h") {
      BOOST_CHECK_MESSAGE(!contains(source, "mMemoryPool;"), source.string());
    }
  }
}
