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
  std::vector<SurfaceDescriptor> catalog;
  for (uint16_t layer = 0; layer < 3; ++layer) {
    SurfaceDescriptor descriptor{LayerId{layer}, layer, 0, SurfaceKind::Cylinder};
    descriptor.referenceCoordinate = static_cast<float>(layer + 1);
    descriptor.chartRange = {-10.f, 10.f};
    catalog.push_back(descriptor);
  }
  return catalog;
}

TrackerInitialization makeConfiguration(const std::vector<SurfaceDescriptor>& catalog,
                                        const std::shared_ptr<BoundedMemoryResource>& pool)
{
  const std::vector<LayerId> ordered{LayerId{2}, LayerId{0}, LayerId{1}};
  TrackerInitialization configuration;
  configuration.catalog = SurfaceCatalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  configuration.memoryPool = pool;
  TrackingParameters parameters;
  parameters.NLayers = 3;
  parameters.LayerxX0.assign(3, 0.f);
  parameters.LayerRadii = {1.f, 2.f, 3.f};
  parameters.LayerResolution.assign(3, 0.001f);
  parameters.SystError2Row.assign(3, 0.f);
  parameters.SystError2Col.assign(3, 0.f);
  parameters.MaxHoles = 1;
  parameters.StartLayerMask = 1u << 2;
  configuration.layout = makeSurfaceLayoutChain(ordered, LayerMask{1u << 1});
  configuration.parameters.push_back(parameters);
  return configuration;
}

bool contains(const std::filesystem::path& path, const std::string& text)
{
  std::ifstream input{path};
  const std::string content{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
  return content.find(text) != std::string::npos;
}
} // namespace

BOOST_AUTO_TEST_CASE(TimeFrameOwnsAndResetsWorkspace)
{
  const auto catalog = makeCatalog();
  const auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  Tracker tracker;
  BOOST_REQUIRE(tracker.initialize(frame, makeConfiguration(catalog, pool)).ok());

  auto& workspace = frame.getScratch();
  auto* workspaceAddress = &workspace;
  const auto* allocator = frame.getMemoryPool().get();
  const auto edgeCapacity = workspace.getNEdges();
  BOOST_REQUIRE(!workspace.getTracklets().empty());
  workspace.getTracklets().front().push_back(Tracklet{});
  frame.getGenericTracks().push_back(GenericTrack{});

  frame.resetTimeFrame();
  BOOST_CHECK_EQUAL(&frame.getScratch(), workspaceAddress);
  BOOST_CHECK_EQUAL(frame.getMemoryPool().get(), allocator);
  BOOST_CHECK_EQUAL(frame.getScratch().getNEdges(), edgeCapacity);
  BOOST_CHECK(frame.getScratch().getTracklets().front().empty());
  BOOST_CHECK(frame.getGenericTracks().empty());
  BOOST_CHECK(frame.isConfigured());

  frame.resetTimeFrame();
  BOOST_CHECK(frame.isConfigured());
}

BOOST_AUTO_TEST_CASE(FailedConfigurationDoesNotReplaceWorkspace)
{
  const auto catalog = makeCatalog();
  TimeFrame frame;
  Tracker tracker;
  BOOST_REQUIRE(tracker.initialize(frame, makeConfiguration(catalog, std::make_shared<BoundedMemoryResource>())).ok());
  auto* oldWorkspace = &frame.getScratch();
  const auto* oldAllocator = frame.getMemoryPool().get();
  const auto oldEdgeCapacity = oldWorkspace->getNEdges();

  TrackerInitialization invalid;
  invalid.catalog = SurfaceCatalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  invalid.memoryPool = std::make_shared<BoundedMemoryResource>();
  invalid.parameters.push_back(TrackingParameters{});
  BOOST_CHECK(!tracker.initialize(frame, invalid).ok());
  BOOST_CHECK(frame.isConfigured());
  BOOST_CHECK_EQUAL(&frame.getScratch(), oldWorkspace);
  BOOST_CHECK_EQUAL(frame.getMemoryPool().get(), oldAllocator);
  BOOST_CHECK_EQUAL(frame.getScratch().getNEdges(), oldEdgeCapacity);
}

BOOST_AUTO_TEST_CASE(ProductionHasOneFrameResetAndNoIndependentLiveScratchOwner)
{
  const auto trackingRoot = std::filesystem::path{__FILE__}.parent_path().parent_path();
  const std::vector<std::filesystem::path> sources{
    trackingRoot / "include/ITSMFTTracking/Tracker.h",
    trackingRoot / "src/Tracker.cxx"};
  for (const auto& source : sources) {
    BOOST_REQUIRE_MESSAGE(std::filesystem::exists(source), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "TimeFrameScratch mScratch"), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "TimeFrameScratch* mScratch"), source.string());
    BOOST_CHECK_MESSAGE(!contains(source, "resetTimeFrameEvent"), source.string());
  }
  BOOST_CHECK(contains(trackingRoot / "src/Tracker.cxx", "frame.resetTimeFrame();"));
  BOOST_CHECK(contains(trackingRoot / "include/ITSMFTTracking/TimeFrame.h", "mScratch"));
  BOOST_CHECK(!contains(trackingRoot / "include/ITSMFTTracking/TimeFrame.h", "virtual void wipe"));
}
