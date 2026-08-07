// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details.

#define BOOST_TEST_MODULE ITSMFT L5 TimeFrame workspace ownership
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

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
                                        const std::shared_ptr<BoundedMemoryResource>& pool,
                                        ClusterSourceId source)
{
  const std::vector<SurfaceId> ordered{SurfaceId{2}, SurfaceId{0}, SurfaceId{1}};
  SurfaceMask owned;
  for (const auto surface : ordered) {
    owned.set(surface);
  }

  TrackerInitialization configuration;
  configuration.catalog = SurfaceCatalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  configuration.memoryPool = pool;
  TrackerIterationConfiguration iteration;
  iteration.graphSubgraphs.push_back(SurfaceGraphSubgraph{ordered, 1, SurfaceMask{1u << 1}, SurfaceMask{1u << 2}});
  TrackingParameters parameters;
  parameters.NLayers = 0;
  iteration.parameters.push_back(parameters);
  iteration.bindings.push_back(SurfacePlanBinding::Declaration{source, owned, ordered,
                                                               SurfaceKind::Cylinder});
  configuration.iterations.push_back(std::move(iteration));
  return configuration;
}

bool contains(const std::filesystem::path& path, const std::string& text)
{
  std::ifstream input{path};
  const std::string content{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
  return content.find(text) != std::string::npos;
}
} // namespace

BOOST_AUTO_TEST_CASE(TimeFrameOwnsWorkspaceAndResetsItOnce)
{
  const auto catalog = makeCatalog();
  const auto pool = std::make_shared<BoundedMemoryResource>();
  constexpr ClusterSourceId source{3};
  TimeFrame frame;
  Tracker tracker;
  tracker.setSource(source);
  BOOST_REQUIRE(tracker.initialize(frame, makeConfiguration(catalog, pool, source)).ok());

  auto& workspace = frame.getWorkspace(source);
  auto* workspaceAddress = &workspace;
  const auto* allocator = frame.getMemoryPool().get();
  const auto capacity = workspace.getNOwnedSurfaces();
  BOOST_REQUIRE(!workspace.getTracklets().empty());
  workspace.getTracklets().front().push_back(o2::its::Tracklet{});
  frame.getCommonTracks().push_back(CommonTrack{});

  frame.resetEvent();
  BOOST_CHECK_EQUAL(frame.getEventResetCount(), 1u);
  BOOST_CHECK_EQUAL(&frame.getWorkspace(source), workspaceAddress);
  BOOST_CHECK_EQUAL(frame.getMemoryPool().get(), allocator);
  BOOST_CHECK_EQUAL(frame.getWorkspace(source).getNOwnedSurfaces(), capacity);
  BOOST_CHECK(frame.getWorkspace(source).getTracklets().front().empty());
  BOOST_CHECK(frame.getCommonTracks().empty());
  BOOST_CHECK(frame.isConfigured());

  frame.resetEvent();
  BOOST_CHECK_EQUAL(frame.getEventResetCount(), 2u);
  BOOST_CHECK(frame.isConfigured());
}

BOOST_AUTO_TEST_CASE(FailedConfigurationDoesNotReplaceWorkspace)
{
  const auto catalog = makeCatalog();
  TimeFrame frame;
  Tracker tracker;
  constexpr ClusterSourceId source{3};
  tracker.setSource(source);
  BOOST_REQUIRE(tracker.initialize(frame, makeConfiguration(catalog, std::make_shared<BoundedMemoryResource>(), source)).ok());
  auto* oldWorkspace = &frame.getWorkspace(source);
  const auto* oldAllocator = frame.getMemoryPool().get();
  const auto oldCapacity = oldWorkspace->getNOwnedSurfaces();

  TrackerInitialization invalid;
  invalid.catalog = SurfaceCatalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  invalid.memoryPool = std::make_shared<BoundedMemoryResource>();
  invalid.iterations.push_back(TrackerIterationConfiguration{});
  BOOST_CHECK(!tracker.initialize(frame, invalid).ok());
  BOOST_CHECK(frame.isConfigured());
  BOOST_CHECK_EQUAL(&frame.getWorkspace(source), oldWorkspace);
  BOOST_CHECK_EQUAL(frame.getMemoryPool().get(), oldAllocator);
  BOOST_CHECK_EQUAL(frame.getWorkspace(source).getNOwnedSurfaces(), oldCapacity);
}

BOOST_AUTO_TEST_CASE(ProductionHasOneFrameResetAndNoIndependentLiveScratchOwner)
{
  const auto trackingRoot = std::filesystem::path{__FILE__}.parent_path().parent_path();
  const std::vector<std::filesystem::path> sources{
    trackingRoot / "include/ITSMFTTracking/Tracker.h",
    trackingRoot / "src/Tracker.cxx"};
  for (const auto& source : sources) {
    BOOST_REQUIRE_MESSAGE(std::filesystem::exists(source), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "SurfaceTrackingScratch mScratch"), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "SurfaceTrackingScratch* mScratch"), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "resetTimeFrameEvent"), source.string());
  }
  BOOST_CHECK(contains(trackingRoot / "src/Tracker.cxx", "frame.resetEvent();"));
  BOOST_CHECK(contains(trackingRoot / "include/ITSMFTTracking/TimeFrame.h", "mWorkspaces"));
  BOOST_CHECK(!contains(trackingRoot / "include/ITSMFTTracking/TimeFrame.h", "virtual void wipe"));
}
