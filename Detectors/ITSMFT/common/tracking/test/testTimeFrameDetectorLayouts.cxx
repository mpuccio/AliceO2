// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#define BOOST_TEST_MODULE ITSMFT TimeFrame detector layouts
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <memory>
#include <vector>

#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITStracking/Constants.h"

#include "TraversalTestSupport.h"

using namespace o2::itsmft::tracking;
using o2::itsmft::TrackingParameters;

namespace
{
bool noopSeedRefit(const TrackSeed&, const TrackingParameters&, float, SurfaceTrackingScratch&,
                   gsl::span<const gsl::span<const GlobalMeasurement>>,
                   gsl::span<const gsl::span<const SurfaceMeasurement>>, SurfaceCatalogView,
                   gsl::span<const SurfaceId>, TrackingCandidate&)
{
  return false;
}

std::vector<SurfaceId> order(size_t count)
{
  std::vector<SurfaceId> result;
  for (uint16_t id = 0; id < count; ++id) {
    result.emplace_back(id);
  }
  return result;
}

float nominalXOverX0(o2::detectors::DetID::ID detector, uint16_t surface)
{
  return detector == o2::detectors::DetID::MFT ? kNominalMFTLayerX0[surface % MFTNLayers] : kNominalITSLayerX0[surface % ITSNLayers];
}

std::vector<SurfaceDescriptor> catalog(size_t count, SurfaceKind kind, o2::detectors::DetID::ID detector)
{
  std::vector<SurfaceDescriptor> result;
  result.reserve(count);
  for (uint16_t id = 0; id < count; ++id) {
    const float xOverX0 = nominalXOverX0(detector, id);
    SurfaceDescriptor descriptor{SurfaceId{id}, id, static_cast<uint8_t>(detector), kind, 0, static_cast<float>(id + 1), 0.f, 100.f};
    descriptor.material.xOverX0 = xOverX0;
    descriptor.material.arealDensityGPerCm2 = xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
    result.push_back(descriptor);
  }
  return result;
}

TrackingParameters parameters(o2::detectors::DetID::ID detector)
{
  TrackingParameters result;
  resetDetectorDefaults(result, detector);
  return result;
}

void configure(TimeFrame& frame, Tracker& tracker, const std::shared_ptr<BoundedMemoryResource>& pool,
               gsl::span<const SurfaceDescriptor> surfaces, gsl::span<const SurfaceId> ordered,
               const std::vector<TrackingParameters>& parameters)
{
  TrackerInitialization configuration;
  configuration.catalog = {surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  configuration.memoryPool = pool;
  for (const auto& params : parameters) {
    TrackerIterationConfiguration iteration;
    iteration.graph = makeSurfaceChain(ordered, params.MaxHoles,
                                       positionalSurfaceMask(params.HoleLayerMask, ordered, params.NLayers),
                                       positionalSurfaceMask(params.StartLayerMask, ordered, params.NLayers));
    iteration.parameters = params;
    configuration.iterations.push_back(std::move(iteration));
  }
  BOOST_REQUIRE(tracker.initialize(frame, configuration).ok());
}
} // namespace

BOOST_AUTO_TEST_CASE(buildSurfaceGraphsRejectsInvalidDefinitions)
{
  const auto surfaces = catalog(7, SurfaceKind::Cylinder, o2::detectors::DetID::ITS);
  const SurfaceCatalogView view{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  const auto ordered = order(7);
  auto invalidCount = parameters(o2::detectors::DetID::ITS);
  invalidCount.NLayers = 8;
  const std::vector<TrackingParameters> countParameters{invalidCount};
  const auto countResult = buildSurfaceGraphs(view, ordered, countParameters);
  BOOST_CHECK(!countResult.ok());
  BOOST_CHECK_EQUAL(countResult.failedIteration, 0u);

  auto invalidHoles = parameters(o2::detectors::DetID::ITS);
  invalidHoles.NLayers = 7;
  invalidHoles.MaxHoles = -1;
  const std::vector<TrackingParameters> holeParameters{invalidHoles};
  const auto holeResult = buildSurfaceGraphs(view, ordered, holeParameters);
  BOOST_CHECK(!holeResult.ok());
  BOOST_CHECK_EQUAL(holeResult.failedIteration, 0u);
}

BOOST_AUTO_TEST_CASE(nonidentity_surface_order_builds_the_expected_binding)
{
  const std::vector<SurfaceId> ordered{SurfaceId{3}, SurfaceId{0}, SurfaceId{6}, SurfaceId{2}, SurfaceId{5}, SurfaceId{1}, SurfaceId{4}};
  auto params = parameters(o2::detectors::DetID::ITS);
  params.NLayers = 7;
  params.MaxHoles = 1;
  params.HoleLayerMask = uint16_t{1} << 1;
  params.StartLayerMask = (uint16_t{1} << 0) | (uint16_t{1} << 4);
  const auto surfaces = catalog(7, SurfaceKind::Cylinder, o2::detectors::DetID::ITS);
  const SurfaceCatalogView view{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  const std::vector<TrackingParameters> graphParameters{params};
  const auto result = buildSurfaceGraphs(view, ordered, graphParameters);
  BOOST_REQUIRE(result.ok());
  const auto& graph = result.graphs.front();
  SurfaceMask owned;
  for (const auto surface : ordered) {
    owned.set(surface);
  }
  const auto binding = SurfacePlanBinding::build(graph.getView(), owned, ordered);
  BOOST_REQUIRE(binding.ok());
  BOOST_CHECK_EQUAL(binding.binding->getOwnedSurfaces().count(), 7);
  BOOST_CHECK(std::is_sorted(binding.binding->getGlobalEdges().begin(), binding.binding->getGlobalEdges().end()));
  BOOST_CHECK(!binding.binding->getGlobalRoadStartCells().empty());
}

BOOST_AUTO_TEST_CASE(traversal_configuration_allocates_one_workspace_per_iteration)
{
  const auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  Tracker tracker{noopSeedRefit};
  TrackerTraits traits;
  auto first = parameters(o2::detectors::DetID::MFT);
  auto second = first;
  second.MaxHoles = 1;
  second.HoleLayerMask = uint16_t{1} << 3;
  const auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  const auto ordered = order(10);
  configure(frame, tracker, pool, surfaces, ordered, {first, second});

  const auto& workspace = frame.getWorkspace();
  BOOST_CHECK_EQUAL(frame.getNIterations(), 2u);
  BOOST_CHECK(!workspace.getTraversalWorkspace(0).valid);
  BOOST_CHECK(!workspace.getTraversalWorkspace(1).valid);
  BOOST_CHECK_NE(&workspace.getTraversalWorkspace(0), &workspace.getTraversalWorkspace(1));
  BOOST_CHECK_NE(frame.getGraph(0).getView().nEdges, frame.getGraph(1).getView().nEdges);
}

BOOST_AUTO_TEST_CASE(configuration_retains_the_selected_binding)
{
  const auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  Tracker tracker{noopSeedRefit};
  auto params = parameters(o2::detectors::DetID::MFT);
  const auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  const auto ordered = order(10);
  configure(frame, tracker, pool, surfaces, ordered, {params});
  BOOST_REQUIRE(frame.getBinding(0) != nullptr);
  BOOST_CHECK_EQUAL(frame.getBinding(0)->getOwnedSurfaces().count(), 10);
  BOOST_CHECK(frame.getGenericTracks().empty());
}

BOOST_AUTO_TEST_CASE(empty_road_start_is_represented_by_the_binding)
{
  const auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  Tracker tracker{noopSeedRefit};
  auto params = parameters(o2::detectors::DetID::MFT);
  params.StartLayerMask = 0;
  const auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  const auto ordered = order(10);
  configure(frame, tracker, pool, surfaces, ordered, {params});
  BOOST_REQUIRE(frame.getBinding(0) != nullptr);
  BOOST_CHECK(frame.getBinding(0)->getGlobalRoadStartCells().empty());
  BOOST_CHECK(frame.getGenericTracks().empty());
}

BOOST_AUTO_TEST_CASE(surface_plan_binding_rejects_cycles_and_accepts_disconnected_kinds)
{
  SurfaceGraph cycle{3};
  BOOST_REQUIRE(cycle.addEdge({SurfaceId{0}, SurfaceId{1}, {}, 0}).isValid());
  BOOST_REQUIRE(cycle.addEdge({SurfaceId{1}, SurfaceId{2}, {}, 0}).isValid());
  BOOST_REQUIRE(cycle.addEdge({SurfaceId{2}, SurfaceId{0}, {}, 0}).isValid());
  BOOST_REQUIRE(cycle.finalize());
  const auto disk = catalog(3, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  SurfaceGraph cyclicGraph{disk, std::move(cycle)};
  const auto cyclicBinding = SurfacePlanBinding::build(cyclicGraph.getView(), SurfaceMask{0x7}, cyclicGraph.getOrderedSurfaces());
  BOOST_CHECK(!cyclicBinding.ok());
}
