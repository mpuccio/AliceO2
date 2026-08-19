// Copyright 2019-2026 CERN and copyright holders of ALICE O2.

#define BOOST_TEST_MODULE ITSMFT DetectorLayout contract
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <type_traits>
#include <filesystem>
#include <fstream>
#include <iterator>

#include "ITSMFTTracking/DetectorLayout.h"

namespace
{
using namespace o2::itsmft::tracking;

std::vector<SurfaceDescriptor> catalog(uint16_t count)
{
  std::vector<SurfaceDescriptor> result;
  for (uint16_t id = 0; id < count; ++id) {
    result.push_back({id < 7 ? id : static_cast<uint16_t>(id - 7), static_cast<uint8_t>(id < 7 ? 0 : 8),
                      id < 7 ? SurfaceKind::Cylinder : SurfaceKind::Disk});
  }
  return result;
}

} // namespace

BOOST_AUTO_TEST_CASE(ValidLayoutProvidesExactlyDensePositionIds)
{
  const auto surfaces = catalog(4);
  DetectorLayoutDefinition definition;
  definition.holeLayers.set(1);
  const DetectorLayout layout{surfaces, definition};
  BOOST_REQUIRE(layout.valid());
  BOOST_CHECK_EQUAL(layout.size(), 4u);
  for (uint16_t id = 0; id < layout.size(); ++id) {
    BOOST_CHECK(&layout[LayerId{id}] == &layout.getLayers()[id]);
  }
  BOOST_CHECK(!layout.getSurfaceCatalog().hasSurface(LayerId{4}));
  BOOST_CHECK(layout.sameComponent(0, 3));
}

BOOST_AUTO_TEST_CASE(ComponentOffsetsExpressDisconnectedLayouts)
{
  const auto surfaces = catalog(17);
  DetectorLayoutDefinition definition;
  definition.componentOffsets = {0, 7};
  const DetectorLayout layout{surfaces, definition};
  BOOST_REQUIRE(layout.valid());
  BOOST_CHECK(layout.sameComponent(0, 6));
  BOOST_CHECK(layout.sameComponent(7, 16));
  BOOST_CHECK(!layout.sameComponent(6, 7));
}

BOOST_AUTO_TEST_CASE(RejectsInvalidBoundaries)
{
  const auto surfaces = catalog(4);
  DetectorLayoutDefinition missingZero;
  missingZero.componentOffsets = {1};
  BOOST_CHECK((DetectorLayout{surfaces, missingZero}.getError() == DetectorLayoutError::InvalidComponentBoundary));
  DetectorLayoutDefinition repeated;
  repeated.componentOffsets = {0, 2, 2};
  BOOST_CHECK((DetectorLayout{surfaces, repeated}.getError() == DetectorLayoutError::InvalidComponentBoundary));
}

BOOST_AUTO_TEST_CASE(HolesMustBeASubsetOfTheLayout)
{
  const auto surfaces = catalog(3);
  DetectorLayoutDefinition holes;
  holes.holeLayers.set(7);
  BOOST_CHECK((DetectorLayout{surfaces, holes}.getError() == DetectorLayoutError::HoleLayersOutsideLayout));
}

BOOST_AUTO_TEST_CASE(LayoutHasNoStaticExpandedTopology)
{
  BOOST_CHECK(std::is_standard_layout_v<DetectorLayoutDefinition>);
  namespace fs = std::filesystem;
  const auto root = fs::path{__FILE__}.parent_path().parent_path();
  std::ifstream input{root / "include/ITSMFTTracking/DetectorLayout.h"};
  const std::string source{std::istreambuf_iterator<char>{input}, {}};
  BOOST_REQUIRE(!source.empty());
  BOOST_CHECK(source.find("std::vector<Edge>") == std::string::npos);
  const auto retiredCell = std::string{"Surface"} + "CellTopology";
  BOOST_CHECK(source.find(retiredCell) == std::string::npos);
  BOOST_CHECK(source.find("pathsByFirstEdge") == std::string::npos);
  BOOST_CHECK(source.find("roadStart") == std::string::npos);
  BOOST_CHECK(source.find("orderedSurfaces") == std::string::npos);
  BOOST_CHECK(source.find("mSurfaceIndexById") == std::string::npos);
  BOOST_CHECK(source.find("surfaceIndicesById") == std::string::npos);
}
