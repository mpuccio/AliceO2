// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details.

#define BOOST_TEST_MODULE ITSMFT TraversalWorkspacePlan
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "ITSMFTTracking/IterationConfiguration.h"

using namespace o2::itsmft::tracking;
using o2::itsmft::TrackingParameters;

namespace
{
DetectorLayout makeLayout()
{
  std::array<SurfaceDescriptor, 4> catalog{};
  for (uint16_t id = 0; id < catalog.size(); ++id) {
    catalog[id] = SurfaceDescriptor{id, 0, SurfaceKind::Cylinder};
  }
  DetectorLayoutDefinition definition;
  definition.holeLayers = LayerMask{uint32_t{1} << 1};
  return DetectorLayout{gsl::span<const SurfaceDescriptor>{catalog.data(), catalog.size()}, std::move(definition)};
}

IterationConfiguration makeConfiguration(const DetectorLayout& layout, const TrackingParameters& parameters)
{
  const auto result = deriveTraversalTopology(layout, parameters);
  if (!result.ok()) {
    throw std::runtime_error{"cannot derive test iteration configuration"};
  }
  IterationConfiguration configuration;
  configuration.parameters = parameters;
  configuration.topology = *result.topology;
  for (uint16_t edge = 0; edge < configuration.topology.edges.size(); ++edge) {
    configuration.edges.emplace_back(edge);
  }
  for (uint16_t cell = 0; cell < configuration.topology.paths.size(); ++cell) {
    configuration.cells.emplace_back(cell);
  }
  return configuration;
}
} // namespace

BOOST_AUTO_TEST_CASE(IterationConfigurationDerivesSelectedTopologyFromTheStaticLayout)
{
  const auto layout = makeLayout();
  o2::itsmft::TrackingParameters firstParameters;
  firstParameters.NLayers = 4;
  firstParameters.MaxHoles = 1;
  firstParameters.StartLayerMask = uint32_t{1} << 3;
  auto secondParameters = firstParameters;
  secondParameters.InactiveLayerMask = uint32_t{1} << 1;

  const auto first = makeConfiguration(layout, firstParameters);
  const auto second = makeConfiguration(layout, secondParameters);

  BOOST_CHECK_EQUAL(first.topology.nLayers, 4u);
  BOOST_CHECK_EQUAL(first.edges.size(), 4u);
  BOOST_CHECK_EQUAL(first.cells.size(), 3u);
  BOOST_CHECK_EQUAL(second.topology.nLayers, 4u);
  BOOST_CHECK_EQUAL(second.topology.activeLayers.count(), 3);
  BOOST_CHECK_EQUAL(second.edges.size(), 2u);
  BOOST_CHECK_EQUAL(second.cells.size(), 1u);
  BOOST_CHECK(second.hasLayer(LayerId{1}));
  BOOST_CHECK(second.getEdgeSlot(EdgeId{1}).has_value()); // 0 -> 2, skipping disabled surface 1
  BOOST_CHECK_EQUAL(second.topology.scheduledPaths.size(), second.cells.size());
}

BOOST_AUTO_TEST_CASE(InvalidLayoutCannotProduceAnIterationConfiguration)
{
  const DetectorLayout invalid{};
  TrackingParameters parameters;
  const auto result = deriveTraversalTopology(invalid, parameters);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == TraversalTopologyError::InvalidLayout);
}

BOOST_AUTO_TEST_CASE(KernelViewBorrowsIterationConfigurationTopology)
{
  const auto layout = makeLayout();
  TrackingParameters parameters;
  parameters.NLayers = 4;
  const auto configuration = makeConfiguration(layout, parameters);
  const auto view = configuration.getTopologyView(layout.getSurfaceCatalog());
  BOOST_CHECK_EQUAL(view.nLayers, configuration.topology.nLayers);
  BOOST_CHECK(view.edges == configuration.topology.edges.data());
  BOOST_CHECK(view.paths == configuration.topology.paths.data());
  BOOST_CHECK(view.pathsByFirstEdgeOffsets == configuration.topology.pathsByFirstEdgeOffsets.data());
  BOOST_CHECK(view.pathsByFirstEdge == configuration.topology.pathsByFirstEdge.data());
  BOOST_CHECK(view.scheduledPaths == configuration.topology.scheduledPaths.data());
  BOOST_CHECK(view.roadStartPaths == configuration.topology.roadStartPaths.data());
  BOOST_CHECK_EQUAL(view.nEdges, 3u);
  BOOST_CHECK_EQUAL(view.nPaths, 2u);
}

BOOST_AUTO_TEST_CASE(ExecutionSourcesUseCurrentTopology)
{
  namespace fs = std::filesystem;
  const auto trackingRoot = fs::path{__FILE__}.parent_path().parent_path();
  const std::array<fs::path, 7> sources{
    trackingRoot / "include/ITSMFTTracking/TrackerTraits.h",
    trackingRoot / "src/TrackerTraits.cxx",
    trackingRoot / "include/ITSMFTTracking/detail/TimeFrameScratch.h",
    trackingRoot / "src/TimeFrameScratch.cxx",
    trackingRoot / "include/ITSMFTTracking/TraversalTopology.h",
    trackingRoot / "src/TraversalTopology.cxx",
    trackingRoot / "src/Tracker.cxx"};
  for (const auto& path : sources) {
    std::ifstream input{path};
    const std::string source{std::istreambuf_iterator<char>{input}, {}};
    BOOST_REQUIRE_MESSAGE(!source.empty(), "cannot read execution source " << path.string());
    BOOST_CHECK_MESSAGE(source.find("template <typename LayoutView>") == std::string::npos,
                        "workspace view retains ignored LayoutView constructor: " << path.string());
    BOOST_CHECK_MESSAGE(source.find("LayerMask skippedSurfaces") == std::string::npos,
                        "lean Edge contract retains skipped-surface storage: " << path.string());
    BOOST_CHECK_MESSAGE(source.find("uint16_t flags") == std::string::npos,
                        "lean Edge contract retains unused flags storage: " << path.string());
  }
}
