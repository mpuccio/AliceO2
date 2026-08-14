#define BOOST_TEST_MODULE ITSMFT LookupCoordinates
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/detail/LookupCoordinates.h"
#include "MFTTracking/Constants.h"

using namespace o2::itsmft::tracking;

namespace
{
GlobalMeasurement measurement(float x, float y, float z, GlobalCovariance3F covariance)
{
  GlobalMeasurement value;
  value.position = {x, y, z};
  value.covariance = covariance;
  return value;
}
} // namespace

BOOST_AUTO_TEST_CASE(DescriptorChartRangesAreFiniteForBothSurfaceKinds)
{
  for (const auto& surface : kITSStaticSurfaceCatalog) {
    BOOST_CHECK(surface.kind == SurfaceKind::Cylinder);
    BOOST_CHECK(surface.chartRange.isValid());
    BOOST_CHECK_EQUAL(surface.chartRange.min, -kITSLookupZHalfExtent[surface.detectorSurfaceIndex]);
    BOOST_CHECK_EQUAL(surface.chartRange.max, kITSLookupZHalfExtent[surface.detectorSurfaceIndex]);
  }
  for (const auto& surface : kMFTStaticSurfaceCatalog) {
    BOOST_CHECK(surface.kind == SurfaceKind::Disk);
    BOOST_CHECK(surface.chartRange.isValid());
    BOOST_CHECK_GT(surface.chartRange.min, 0.f);
    BOOST_CHECK_EQUAL(surface.chartRange.min, o2::mft::constants::index_table::RMin[surface.detectorSurfaceIndex]);
    BOOST_CHECK_EQUAL(surface.chartRange.max, o2::mft::constants::index_table::RMax[surface.detectorSurfaceIndex]);
  }
}

BOOST_AUTO_TEST_CASE(DiskTransformKeepsPhiRadialCorrelation)
{
  const auto input = measurement(3.f, 4.f, -50.f, {4.f, 1.f, 0.f, 9.f, 0.f, 16.f});
  LookupCoordinates coordinates;
  BOOST_REQUIRE(makeLookupCoordinates(kMFTStaticSurfaceCatalog[0], input, coordinates));
  BOOST_CHECK_CLOSE_FRACTION(coordinates.transverse, 5.f, 1.e-6f);
  BOOST_CHECK_CLOSE_FRACTION(coordinates.phi, std::atan2(4.f, 3.f), 1.e-6f);
  BOOST_CHECK_NE(coordinates.covariance[1], 0.f);
  BOOST_CHECK_GE(coordinates.covariance[0] * coordinates.covariance[2] -
                   coordinates.covariance[1] * coordinates.covariance[1],
                 0.f);
}

BOOST_AUTO_TEST_CASE(CylinderTransformUsesPhiAndZWithCorrelation)
{
  const auto input = measurement(3.f, 4.f, 7.f, {4.f, 1.f, 2.f, 9.f, 3.f, 16.f});
  LookupCoordinates coordinates;
  BOOST_REQUIRE(makeLookupCoordinates(kITSStaticSurfaceCatalog[0], input, coordinates));
  BOOST_CHECK_CLOSE_FRACTION(coordinates.transverse, 7.f, 1.e-6f);
  BOOST_CHECK_NE(coordinates.covariance[1], 0.f);
}

BOOST_AUTO_TEST_CASE(ZeroRadiusAndUndefinedSurfaceFailClosed)
{
  LookupCoordinates coordinates;
  BOOST_CHECK(!makeLookupCoordinates(kMFTStaticSurfaceCatalog[0], measurement(0.f, 0.f, 0.f, {}), coordinates));
  auto surface = kMFTStaticSurfaceCatalog[0];
  surface.kind = SurfaceKind::Undefined;
  BOOST_CHECK(!makeLookupCoordinates(surface, measurement(1.f, 0.f, 0.f, {}), coordinates));
}

BOOST_AUTO_TEST_CASE(WindowWrapsPhiAndClampsToDescriptorRange)
{
  LookupCoordinates coordinates{0.01f, 2.15f, {0.01f, 0.f, 0.04f}};
  LookupWindow window;
  BOOST_REQUIRE(makeLookupWindow(coordinates, kMFTStaticSurfaceCatalog[0].chartRange, 1.f, window));
  BOOST_CHECK(window.wrapsPhi);
  BOOST_CHECK_EQUAL(window.transverseMin, kMFTStaticSurfaceCatalog[0].chartRange.min);
  BOOST_CHECK_LT(window.transverseMax, kMFTStaticSurfaceCatalog[0].chartRange.max);
}

BOOST_AUTO_TEST_CASE(LookupTransformDoesNotDispatchOnDetectorIdentity)
{
  const auto header = std::filesystem::path{__FILE__}.parent_path().parent_path() /
                      "include" / "ITSMFTTracking" / "detail" / "LookupCoordinates.h";
  std::ifstream stream{header};
  BOOST_REQUIRE(stream);
  const std::string source{std::istreambuf_iterator<char>{stream}, std::istreambuf_iterator<char>{}};
  BOOST_CHECK(source.find("detectorId") == std::string::npos);
  BOOST_CHECK(source.find("ClusterSourceId") == std::string::npos);
}
