// Copyright 2019-2026 CERN and copyright holders of ALICE O2.

#define BOOST_TEST_MODULE ITSMFT SurfaceLayout contract
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <type_traits>
#include <filesystem>
#include <fstream>
#include <iterator>

#include "ITSMFTTracking/SurfaceLayout.h"

namespace
{
using namespace o2::itsmft::tracking;

std::vector<SurfaceDescriptor> catalog(uint16_t count)
{
  std::vector<SurfaceDescriptor> result;
  for (uint16_t id = 0; id < count; ++id) {
    result.push_back({LayerId{id}, id, 0, id < 7 ? SurfaceKind::Cylinder : SurfaceKind::Disk});
  }
  return result;
}

SurfaceLayoutDefinition ordered(uint16_t count)
{
  SurfaceLayoutDefinition result;
  for (uint16_t id = 0; id < count; ++id) {
    result.orderedSurfaces.push_back(LayerId{id});
  }
  return result;
}
} // namespace

BOOST_AUTO_TEST_CASE(ValidOrderedLayoutProvidesValidatedCatalog)
{
  const auto surfaces = catalog(4);
  auto definition = ordered(4);
  definition.holeLayers.set(1);
  const SurfaceLayout layout{surfaces, definition};
  BOOST_REQUIRE(layout.valid());
  BOOST_CHECK_EQUAL(layout.getOrderedSurfaces().size(), 4u);
  BOOST_CHECK_EQUAL(layout.getSurfaceCatalog().getSurface(LayerId{3}).id.value(), 3u);
  BOOST_CHECK(layout.sameComponent(0, 3));
}

BOOST_AUTO_TEST_CASE(ComponentOffsetsExpressDisconnectedLayouts)
{
  const auto surfaces = catalog(17);
  auto definition = ordered(17);
  definition.componentOffsets = {0, 7};
  const SurfaceLayout layout{surfaces, definition};
  BOOST_REQUIRE(layout.valid());
  BOOST_CHECK(layout.sameComponent(0, 6));
  BOOST_CHECK(layout.sameComponent(7, 16));
  BOOST_CHECK(!layout.sameComponent(6, 7));
}

BOOST_AUTO_TEST_CASE(RejectsInvalidSurfacesAndBoundaries)
{
  const auto surfaces = catalog(4);
  auto duplicate = ordered(4);
  duplicate.orderedSurfaces[3] = LayerId{2};
  BOOST_CHECK((SurfaceLayout{surfaces, duplicate}.getError() == SurfaceLayoutError::DuplicateSurface));
  auto outside = ordered(4);
  outside.orderedSurfaces[3] = LayerId{7};
  BOOST_CHECK((SurfaceLayout{surfaces, outside}.getError() == SurfaceLayoutError::InvalidSurface));
  auto missingZero = ordered(4);
  missingZero.componentOffsets = {1};
  BOOST_CHECK((SurfaceLayout{surfaces, missingZero}.getError() == SurfaceLayoutError::InvalidComponentBoundary));
  auto repeated = ordered(4);
  repeated.componentOffsets = {0, 2, 2};
  BOOST_CHECK((SurfaceLayout{surfaces, repeated}.getError() == SurfaceLayoutError::InvalidComponentBoundary));
}

BOOST_AUTO_TEST_CASE(HolesMustBeASubsetOfTheLayout)
{
  const auto surfaces = catalog(3);
  auto holes = ordered(3);
  holes.holeLayers.set(7);
  BOOST_CHECK((SurfaceLayout{surfaces, holes}.getError() == SurfaceLayoutError::HoleLayersOutsideLayout));
}

BOOST_AUTO_TEST_CASE(LayoutHasNoStaticExpandedTopology)
{
  BOOST_CHECK(std::is_standard_layout_v<SurfaceLayoutDefinition>);
  namespace fs = std::filesystem;
  const auto root = fs::path{__FILE__}.parent_path().parent_path();
  std::ifstream input{root / "include/ITSMFTTracking/SurfaceLayout.h"};
  const std::string source{std::istreambuf_iterator<char>{input}, {}};
  BOOST_REQUIRE(!source.empty());
  BOOST_CHECK(source.find("std::vector<Edge>") == std::string::npos);
  const auto retiredCell = std::string{"Surface"} + "CellTopology";
  BOOST_CHECK(source.find(retiredCell) == std::string::npos);
  BOOST_CHECK(source.find("pathsByFirstEdge") == std::string::npos);
  BOOST_CHECK(source.find("roadStart") == std::string::npos);
}
