// Copyright 2019-2026 CERN and copyright holders of ALICE O2.

#define BOOST_TEST_MODULE ITSMFT CombinedStaticSurfaceCatalogTopology
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/TraversalTopology.h"

namespace
{
using namespace o2::itsmft::tracking;
using o2::itsmft::TrackingParameters;

} // namespace

BOOST_AUTO_TEST_CASE(CombinedStaticCatalogHasDenseGlobalIdsAndQualifiedIdentity)
{
  BOOST_REQUIRE_EQUAL(kITSMFTCombinedStaticSurfaceCatalog.size(), static_cast<std::size_t>(ITSNLayers + MFTNLayers));
  for (uint16_t i = 0; i < kITSMFTCombinedStaticSurfaceCatalog.size(); ++i) {
    const auto& surface = kITSMFTCombinedStaticSurfaceCatalog[i];
    BOOST_CHECK(&surface == &kITSMFTCombinedStaticSurfaceCatalog[i]);
    if (i < ITSNLayers) {
      BOOST_CHECK_EQUAL(surface.detectorId, static_cast<uint8_t>(o2::detectors::DetID::ITS));
      BOOST_CHECK(surface.kind == SurfaceKind::Cylinder);
      BOOST_CHECK_EQUAL(surface.detectorSurfaceIndex, i);
    } else {
      BOOST_CHECK_EQUAL(surface.detectorId, static_cast<uint8_t>(o2::detectors::DetID::MFT));
      BOOST_CHECK(surface.kind == SurfaceKind::Disk);
      BOOST_CHECK_EQUAL(surface.detectorSurfaceIndex, i - ITSNLayers);
    }
  }
}
BOOST_AUTO_TEST_CASE(CombinedStaticCatalogDerivesDisconnectedHoleTopology)
{
  DetectorLayoutDefinition definition;
  definition.componentOffsets = {0, ITSNLayers};
  definition.holeLayers.set(3);
  definition.holeLayers.set(static_cast<uint16_t>(ITSNLayers + 5));

  const auto layout = DetectorLayout{kITSMFTCombinedStaticSurfaceCatalog, std::move(definition)};
  TrackingParameters parameters;
  parameters.NLayers = static_cast<int>(layout.size());
  parameters.MaxHoles = 1;
  const auto result = deriveTraversalTopology(layout, parameters);
  BOOST_REQUIRE(result.ok());
  const auto& topology = *result.topology;
  BOOST_CHECK_EQUAL(topology.edges.size(), 17u);
  BOOST_CHECK_EQUAL(topology.paths.size(), 17u);

  for (const auto& edge : topology.edges) {
    const bool its = edge.from.value() < ITSNLayers && edge.to.value() < ITSNLayers;
    const bool mft = edge.from.value() >= ITSNLayers && edge.to.value() >= ITSNLayers;
    BOOST_CHECK(its || mft);
  }
  const auto view = topology.getView(layout.getSurfaceCatalog());
  uint32_t csrEntries = 0;
  for (uint32_t edge = 0; edge < view.nEdges; ++edge) {
    csrEntries += view.getPathsStartingWithEdge(EdgeId{static_cast<uint16_t>(edge)}).getEntries();
  }
  BOOST_CHECK_EQUAL(csrEntries, topology.paths.size());
}
