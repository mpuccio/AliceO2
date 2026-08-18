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
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TraversalTopology.h"
#include "ITStracking/Constants.h"

#include "TraversalTestSupport.h"

using namespace o2::itsmft::tracking;
using o2::itsmft::TrackingParameters;

namespace
{
bool noopSeedRefit(const TrackSeed&, const TimeFrame&, const TrackingParameters&, float,
                   gsl::span<const gsl::span<const GlobalMeasurement>>,
                   SurfaceCatalogView, gsl::span<const LayerId>, TrackingCandidate&)
{
  return false;
}

std::vector<LayerId> order(size_t count)
{
  std::vector<LayerId> result;
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
    SurfaceDescriptor descriptor{LayerId{id}, id, static_cast<uint8_t>(detector), kind, 0, static_cast<float>(id + 1), 0.f, 100.f};
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
               gsl::span<const SurfaceDescriptor> surfaces, gsl::span<const LayerId> ordered,
               const std::vector<TrackingParameters>& parameters)
{
  TrackerInitialization configuration;
  configuration.catalog = {surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  configuration.memoryPool = pool;
  for (const auto& params : parameters) {
    TrackerIterationConfiguration iteration;
    iteration.layout = makeSurfaceLayoutChain(ordered, params.MaxHoles,
                                              positionalSurfaceMask(params.HoleLayerMask, ordered, params.NLayers),
                                              positionalSurfaceMask(params.StartLayerMask, ordered, params.NLayers));
    iteration.parameters = params;
    configuration.iterations.push_back(std::move(iteration));
  }
  BOOST_REQUIRE(tracker.initialize(frame, configuration).ok());
}
} // namespace

BOOST_AUTO_TEST_CASE(surfaceLayoutsRejectInvalidDefinitions)
{
  const auto surfaces = catalog(7, SurfaceKind::Cylinder, o2::detectors::DetID::ITS);
  const SurfaceCatalogView view{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  const gsl::span<const SurfaceDescriptor> catalogSpan{surfaces.data(), surfaces.size()};
  const auto ordered = order(7);
  auto invalidOrder = ordered;
  invalidOrder.push_back(LayerId{7});
  const auto invalidSurface = SurfaceLayout{catalogSpan, makeSurfaceLayoutChain(invalidOrder)};
  BOOST_CHECK(!invalidSurface.valid());
  BOOST_CHECK_EQUAL(static_cast<int>(invalidSurface.getError()), static_cast<int>(SurfaceLayoutError::InvalidSurface));

  const auto invalidHoleLayout = SurfaceLayout{catalogSpan, makeSurfaceLayoutChain(ordered, -1)};
  BOOST_CHECK(!invalidHoleLayout.valid());
  BOOST_CHECK_EQUAL(static_cast<int>(invalidHoleLayout.getError()), static_cast<int>(SurfaceLayoutError::NegativeMaxHoles));
}

BOOST_AUTO_TEST_CASE(nonidentity_surface_order_builds_the_expected_topology)
{
  const std::vector<LayerId> ordered{LayerId{3}, LayerId{0}, LayerId{6}, LayerId{2}, LayerId{5}, LayerId{1}, LayerId{4}};
  auto params = parameters(o2::detectors::DetID::ITS);
  params.NLayers = 7;
  params.MaxHoles = 1;
  params.HoleLayerMask = uint16_t{1} << 1;
  params.StartLayerMask = (uint16_t{1} << 0) | (uint16_t{1} << 4);
  const auto surfaces = catalog(7, SurfaceKind::Cylinder, o2::detectors::DetID::ITS);
  const SurfaceCatalogView view{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  const auto layout = SurfaceLayout{gsl::span<const SurfaceDescriptor>{surfaces.data(), surfaces.size()},
                                    makeSurfaceLayoutChain(ordered, params.MaxHoles,
                                                           positionalSurfaceMask(params.HoleLayerMask, ordered, params.NLayers),
                                                           positionalSurfaceMask(params.StartLayerMask, ordered, params.NLayers))};
  BOOST_REQUIRE(layout.valid());
  const auto result = deriveTraversalTopology(layout);
  BOOST_REQUIRE(result.ok());
  const auto& topology = *result.topology;
  BOOST_CHECK_EQUAL(topology.orderedSurfaces.size(), 7u);
  for (size_t position = 0; position < ordered.size(); ++position) {
    BOOST_CHECK_EQUAL(topology.orderedSurfaces[position].value(), ordered[position].value());
  }
  BOOST_CHECK(!topology.edges.empty());
  BOOST_CHECK(!topology.paths.empty());
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

  auto& workspace = frame.getWorkspace();
  BOOST_CHECK_EQUAL(frame.getNIterations(), 2u);
  BOOST_CHECK(!workspace.getTraversalWorkspace(0).valid);
  BOOST_CHECK(!workspace.getTraversalWorkspace(1).valid);
  BOOST_CHECK_NE(&workspace.getTraversalWorkspace(0), &workspace.getTraversalWorkspace(1));
  TrackerTestAccess::preparePlan(tracker, workspace.getTraversalWorkspace(0), frame.getLayout(0));
  TrackerTestAccess::preparePlan(tracker, workspace.getTraversalWorkspace(1), frame.getLayout(1));
  BOOST_CHECK_NE(workspace.getTraversalWorkspace(0).getTopologyView().nEdges,
                 workspace.getTraversalWorkspace(1).getTopologyView().nEdges);
}

BOOST_AUTO_TEST_CASE(configuration_retains_the_selected_workspace_plan)
{
  const auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  Tracker tracker{noopSeedRefit};
  auto params = parameters(o2::detectors::DetID::MFT);
  const auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  const auto ordered = order(10);
  configure(frame, tracker, pool, surfaces, ordered, {params});
  auto& workspace = frame.getWorkspace().getTraversalWorkspace(0);
  TrackerTestAccess::preparePlan(tracker, workspace, frame.getLayout(0));
  BOOST_CHECK_EQUAL(workspace.orderedSurfaces.size(), 10u);
  for (size_t position = 0; position < ordered.size(); ++position) {
    BOOST_CHECK_EQUAL(workspace.orderedSurfaces[position].value(), ordered[position].value());
  }
  BOOST_CHECK(frame.getGenericTracks().empty());
}

BOOST_AUTO_TEST_CASE(empty_road_start_is_represented_by_the_workspace_plan)
{
  const auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  Tracker tracker{noopSeedRefit};
  auto params = parameters(o2::detectors::DetID::MFT);
  params.StartLayerMask = 0;
  const auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  const auto ordered = order(10);
  configure(frame, tracker, pool, surfaces, ordered, {params});
  auto& workspace = frame.getWorkspace().getTraversalWorkspace(0);
  TrackerTestAccess::preparePlan(tracker, workspace, frame.getLayout(0));
  BOOST_CHECK(workspace.roadStartCells.empty());
  BOOST_CHECK(frame.getGenericTracks().empty());
}
