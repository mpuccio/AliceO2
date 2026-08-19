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
using o2::itsmft::IterationParameters;
using o2::itsmft::TrackingParameters;

namespace
{
bool noopSeedRefit(const TrackSeed&, const TimeFrame&, const IterationParameters&, float,
                   gsl::span<const gsl::span<const GlobalMeasurement>>,
                   SurfaceCatalogView, TrackingCandidate&)
{
  return false;
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
    SurfaceDescriptor descriptor{id, static_cast<uint8_t>(detector), kind};
    descriptor.referenceCoordinate = static_cast<float>(id + 1);
    descriptor.chartRange = {0.f, 100.f};
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
               gsl::span<const SurfaceDescriptor> surfaces,
               const std::vector<TrackingParameters>& parameters, LayerMask holeLayers = {})
{
  TrackerInitialization configuration;
  configuration.catalog = {surfaces.data(), static_cast<uint32_t>(surfaces.size())};
  configuration.memoryPool = pool;
  configuration.layout = makeDetectorLayout(holeLayers);
  configuration.parameters = parameters;
  BOOST_REQUIRE(tracker.initialize(frame, configuration).ok());
}
} // namespace

BOOST_AUTO_TEST_CASE(detectorLayoutsRejectInvalidDefinitions)
{
  const auto surfaces = catalog(7, SurfaceKind::Cylinder, o2::detectors::DetID::ITS);
  const gsl::span<const SurfaceDescriptor> catalogSpan{surfaces.data(), surfaces.size()};
  const auto emptyLayout = DetectorLayout{};
  BOOST_CHECK(!emptyLayout.valid());
  BOOST_CHECK_EQUAL(static_cast<int>(emptyLayout.getError()), static_cast<int>(DetectorLayoutError::EmptyCatalog));

  const auto invalidHoleLayout = DetectorLayout{catalogSpan, makeDetectorLayout(LayerMask{1u << 7})};
  BOOST_CHECK(!invalidHoleLayout.valid());
  BOOST_CHECK_EQUAL(static_cast<int>(invalidHoleLayout.getError()), static_cast<int>(DetectorLayoutError::HoleLayersOutsideLayout));
}

BOOST_AUTO_TEST_CASE(dense_layout_positions_build_the_expected_topology)
{
  auto params = parameters(o2::detectors::DetID::ITS);
  params.NLayers = 7;
  params.MaxHoles = 1;
  params.StartLayerMask = (uint16_t{1} << 0) | (uint16_t{1} << 4);
  params.InactiveLayerMask = uint16_t{1} << 1;
  const auto surfaces = catalog(7, SurfaceKind::Cylinder, o2::detectors::DetID::ITS);
  const auto layout = DetectorLayout{gsl::span<const SurfaceDescriptor>{surfaces.data(), surfaces.size()},
                                     makeDetectorLayout(LayerMask{1u << 1})};
  BOOST_REQUIRE(layout.valid());
  const auto result = deriveTraversalTopology(layout, params);
  BOOST_REQUIRE(result.ok());
  const auto& topology = *result.topology;
  BOOST_CHECK_EQUAL(topology.nLayers, 7u);
  BOOST_REQUIRE_EQUAL(topology.activeSurfaceList.size(), 6u);
  BOOST_CHECK_EQUAL(topology.activeSurfaceList[0].value(), 0u);
  BOOST_CHECK_EQUAL(topology.activeSurfaceList[1].value(), 2u);
  BOOST_CHECK_EQUAL(topology.activeSurfaceList.back().value(), 6u);
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
  const auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  configure(frame, tracker, pool, surfaces, {first, second}, LayerMask{1u << 3});

  auto& scratch = frame.getScratch();
  BOOST_CHECK_EQUAL(tracker.getIterationConfigurations().size(), 2u);
  BOOST_CHECK(frame.getGenericTracks().empty());
  BOOST_CHECK_NE(tracker.getIterationConfigurations()[0].topology.edges.size(),
                 tracker.getIterationConfigurations()[1].topology.edges.size());
}

BOOST_AUTO_TEST_CASE(configuration_retains_the_selected_workspace_plan)
{
  const auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  Tracker tracker{noopSeedRefit};
  auto params = parameters(o2::detectors::DetID::MFT);
  const auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  configure(frame, tracker, pool, surfaces, {params});
  const auto& configuration = tracker.getIterationConfigurations()[0];
  BOOST_CHECK_EQUAL(configuration.topology.nLayers, 10u);
  BOOST_REQUIRE_EQUAL(configuration.topology.activeSurfaceList.size(), 10u);
  for (uint16_t position = 0; position < 10; ++position) {
    BOOST_CHECK_EQUAL(configuration.topology.activeSurfaceList[position].value(), position);
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
  configure(frame, tracker, pool, surfaces, {params});
  const auto& configuration = tracker.getIterationConfigurations()[0];
  BOOST_CHECK(configuration.topology.roadStartPaths.empty());
  BOOST_CHECK(frame.getGenericTracks().empty());
}
