// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file testPlacedLayouts.cxx
/// \brief Placed-layout extents, empty tracking and curved-sensor geometry.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK placed layouts
#include <boost/test/unit_test.hpp>
#include "GenericTracking/Base/LayerEnvelopes.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/Tracker.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GridTestSupport.h"
#include "UniformMaterial.h"
#include "ToyBarrel3.h"
#include "ToyBarrel4Disc8.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>
namespace o2::gtrk
{
using detectors::ITS;
using detectors::MFT;
using detectors::ToyBarrel3;
using detectors::ToyBarrel4Disc8;
constexpr double kPositionTolerance = 1e-9;
static_assert(DetectorLayout<ToyBarrel3>);
static_assert(DetectorLayout<ToyBarrel4Disc8>);

// The topology each toy exists to represent.
static_assert(Detector<ToyBarrel3>::hasBarrel && !Detector<ToyBarrel3>::hasDiscs);
static_assert(Detector<ToyBarrel4Disc8>::hasBarrel && Detector<ToyBarrel4Disc8>::hasDiscs);
static_assert(Detector<ToyBarrel3>::nLayers == 3);
static_assert(Detector<ToyBarrel4Disc8>::nBarrelLayers == 4);
static_assert(Detector<ToyBarrel4Disc8>::nDiscLayers == 8);
static_assert(Detector<ToyBarrel4Disc8>::isBarrelLayer(0));
static_assert(Detector<ToyBarrel4Disc8>::isBarrelLayer(3));
static_assert(!Detector<ToyBarrel4Disc8>::isBarrelLayer(4));
static_assert(Detector<ToyBarrel4Disc8>::discIndex(4) == 0U);

// Sensor counts are derived, never typed in twice, so this cannot drift.
static_assert(ToyBarrel3::nSensors == 12);

namespace
{
static_assert(Detector<ITS>::nBarrelLayers == 7 && Detector<ITS>::nDiscLayers == 0);
static_assert(Detector<MFT>::nBarrelLayers == 0 && Detector<MFT>::nDiscLayers == 10);
static_assert(ITS::nSensors == 24120 && MFT::nSensors == 936);

/// Every link joins two layers of its detector, and none joins a layer to itself.
template <DetectorLayout Det>
void checkLinks()
{
  for (const auto& link : Det::layerLinks) {
    BOOST_TEST(link.from < Detector<Det>::nLayers);
    BOOST_TEST(link.to < Detector<Det>::nLayers);
    BOOST_TEST(link.from != link.to);
    BOOST_TEST(link.sourceMin < link.sourceMax);
  }
}

/// Where a stand-in layer sits: the ITS's radii and the MFT's disc z, roughly.
template <DetectorLayout Det>
double standInPosition(LayerID layer)
{
  return Detector<Det>::isBarrelLayer(layer) ? 2.3 + (6.2 * static_cast<double>(layer)) : -45.3 - (3.6 * static_cast<double>(layer));
}

/// A stand-in placement: every chip on its layer at the same place, facing
/// along x on a barrel layer and along z on a disc, with the disc chips 8 cm
/// from the beam. Enough to build the tracking geometry; not a detector.
template <DetectorLayout Det>
std::vector<Sensor> standInPlacements(double shiftZ = 0.)
{
  std::vector<Sensor> sensors;
  for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    for (int chip = 0; chip < Detector<Det>::sensorCount(id); ++chip) {
      Sensor sensor{};
      if (Detector<Det>::isBarrelLayer(id)) {
        sensor.rot[3] = sensor.rot[7] = sensor.rot[2] = 1.; // local x, y, z along global y, z, x
        sensor.trans[0] = standInPosition<Det>(id);
      } else {
        sensor.rot[0] = sensor.rot[4] = sensor.rot[8] = 1.; // local axes along the global ones
        sensor.trans[0] = 8.;
        sensor.trans[2] = standInPosition<Det>(id);
      }
      sensor.trans[2] += shiftZ;
      sensor.layer = id;
      sensors.push_back(sensor);
    }
  }
  return sensors;
}

/// The whole tracker, run on an empty timeframe: that it compiles and runs
/// for a barrel-only and a disc-only placed layout is the point.
template <DetectorLayout Det>
void runEmpty()
{
  const auto sensors = standInPlacements<Det>();
  const auto geometry = std::make_unique<const Geometry<Det>>(sensors);
  const auto tracking = std::make_unique<const TrackingGeometry<Det>>(*geometry);
  BOOST_TEST(tracking->surface(0).bounds.halfX == static_cast<Scalar>(detectors::alpide::kSensor.halfX));

  const SearchWindowParams windows{.minPtGeV = 0.3, .bzKiloGauss = 5.};
  const VertexConstraint beamSpot{.z = 0., .sigmaZ = 6., .sigmaTransverse = 0.002};
  const GridGeometry<Det> grid{20'000'000, test::withToyFrames(GridParams{}), test::oneIterationWindows<Det>(windows, beamSpot)};
  std::vector<std::uint32_t> offsets(static_cast<std::size_t>(grid.bucketCount()) + 1U, 0U);
  const MeasurementStore columns;
  const TimeFrame<Det> frame{measurementView(columns), grid, std::span<const std::uint32_t>{offsets}};

  const IterationSettings settings{.windows = windows, .vertexSource = VertexSource::kBeamLine};
  Tracker<Det, UniformMaterial> tracker{settings, beamSpot, scatteringTable<Det>(uniformLayerMaterial<Det>(0.01)),
                                        std::make_shared<BoundedMemoryResource>(), 1U};
  tracker.setTrackModel(*tracking, UniformMaterial{}, Propagator{.bzKiloGauss = 5.F});
  const auto result = tracker.run(frame, tracker.planTimeframe(frame, {}, 1, 0.));
  BOOST_TEST((result.status == StageStatus::kComplete));
  BOOST_REQUIRE(result.tracks != nullptr);
  BOOST_TEST(result.tracks->size() == 0U);
}

} // namespace

BOOST_AUTO_TEST_CASE(ExtentsNeedAGeometry)
{
  BOOST_CHECK_THROW(static_cast<void>(layerDistance<ITS>(0)), std::logic_error);
  BOOST_CHECK_THROW(static_cast<void>(discZ<MFT>(0)), std::logic_error);
}

BOOST_AUTO_TEST_CASE(LinkTablesJoinLayersOfTheirDetector)
{
  checkLinks<ITS>();
  checkLinks<MFT>();
}

BOOST_AUTO_TEST_CASE(ExtentsAreMeasuredFromThePlacements)
{
  const auto& half = detectors::alpide::kSensor;
  {
    const Geometry<ITS> geometry{standInPlacements<ITS>()};
    for (LayerID layer = 0; layer < Detector<ITS>::nLayers; ++layer) {
      const auto& e = PlacedEnvelopes<ITS>::all()[layer];
      const double radius = standInPosition<ITS>(layer);
      BOOST_TEST(layerDistance<ITS>(layer) == e.position);
      BOOST_TEST(e.position == radius, boost::test_tools::tolerance(1e-9));
      // The chip faces along x: its thickness spans r, its columns z.
      BOOST_TEST(radiusMin<ITS>(layer) == e.rMin);
      BOOST_TEST(e.rMin <= radius - half.halfThickness);
      BOOST_TEST(e.rMin > radius - half.halfThickness - 2. * kEnvelopeResolution);
      BOOST_TEST(radiusMax<ITS>(layer) == e.rMax);
      BOOST_TEST(e.rMax >= std::hypot(radius + half.halfThickness, half.halfX));
      BOOST_TEST(secondMin<ITS>(layer) == e.zMin);
      BOOST_TEST(e.zMin <= -half.halfY);
      BOOST_TEST(secondMax<ITS>(layer) == e.zMax);
      BOOST_TEST(e.zMax >= half.halfY);
      BOOST_TEST(discZ<ITS>(layer) == 0.);
    }
  }
  {
    const Geometry<MFT> geometry{standInPlacements<MFT>()};
    for (LayerID layer = 0; layer < Detector<MFT>::nLayers; ++layer) {
      const auto& e = PlacedEnvelopes<MFT>::all()[layer];
      const double z = standInPosition<MFT>(layer);
      BOOST_TEST(discZ<MFT>(layer) == e.position);
      BOOST_TEST(e.position == z, boost::test_tools::tolerance(1e-9));
      BOOST_TEST(layerDistance<MFT>(layer) == -e.position); // the MFT is on the C side
      BOOST_TEST(secondMin<MFT>(layer) == e.rMin);
      BOOST_TEST(e.rMin <= 8. - half.halfX);
      BOOST_TEST(secondMax<MFT>(layer) == e.rMax);
      BOOST_TEST(e.rMax >= std::hypot(8. + half.halfX, half.halfY));
      BOOST_TEST(zMin<MFT>(layer) == e.zMin);
      BOOST_TEST(e.zMin <= z - half.halfThickness);
      BOOST_TEST(zMax<MFT>(layer) == e.zMax);
      BOOST_TEST(e.zMax >= z + half.halfThickness);
    }
  }
  // A later geometry of the same layout replaces the table: shifted by 1 cm
  // along z, every z extent follows to the micron.
  {
    const Geometry<MFT> shifted{standInPlacements<MFT>(1.)};
    for (LayerID layer = 0; layer < Detector<MFT>::nLayers; ++layer) {
      BOOST_TEST(discZ<MFT>(layer) == standInPosition<MFT>(layer) + 1., boost::test_tools::tolerance(1e-9));
      BOOST_TEST(zMin<MFT>(layer) == PlacedEnvelopes<MFT>::all()[layer].zMin);
    }
  }
}

BOOST_AUTO_TEST_CASE(GridPointsAreAboutTheBeamLine)
{
  constexpr double kBeamX = -0.057;
  constexpr double kBeamY = -0.0565;
  auto moved = standInPlacements<ITS>();
  for (auto& sensor : moved) {
    sensor.trans[0] += kBeamX;
    sensor.trans[1] += kBeamY;
  }
  const Geometry<ITS> reference{standInPlacements<ITS>()};
  const Geometry<ITS> shifted{moved};
  const auto tolerance = boost::test_tools::tolerance(1e-5F);
  for (const SensorID chip : {SensorID{0}, SensorID{107}, SensorID{24119}}) {
    for (const auto& [loc0, loc1] : {std::pair{0.3, -1.1}, std::pair{-0.6, 0.9}}) {
      const auto a = gridPoint<ITS>(reference, chip, loc0, loc1);
      const auto b = gridPoint<ITS>(shifted, chip, loc0, loc1, kBeamX, kBeamY);
      BOOST_TEST(b.layer == a.layer);
      BOOST_TEST(b.x == a.x, tolerance);
      BOOST_TEST(b.y == a.y, tolerance);
      BOOST_TEST(b.z == a.z, tolerance);
      BOOST_TEST(b.r == a.r, tolerance);
      BOOST_TEST(b.phi == a.phi, tolerance);
      BOOST_TEST(b.second == a.second, tolerance);
      const auto va = pointVariances<ITS>(reference, chip, loc0, a, 1e-6F, 0.F, 4e-6F);
      const auto vb = pointVariances<ITS>(shifted, chip, loc0, b, 1e-6F, 0.F, 4e-6F);
      BOOST_TEST(vb.z == va.z, tolerance);
      BOOST_TEST(vb.r == va.r, tolerance);
      BOOST_TEST(vb.rPhi == va.rPhi, tolerance);
      // And without the subtraction the moved geometry is off by the beam.
      const auto c = gridPoint<ITS>(shifted, chip, loc0, loc1);
      BOOST_TEST(static_cast<double>(c.x - a.x) == kBeamX, boost::test_tools::tolerance(1e-3));
      BOOST_TEST(static_cast<double>(c.y - a.y) == kBeamY, boost::test_tools::tolerance(1e-3));
    }
  }
}

BOOST_AUTO_TEST_CASE(PlacedGeometryRejectsAWrongCount)
{
  auto sensors = standInPlacements<MFT>();
  sensors.pop_back();
  BOOST_CHECK_THROW(Geometry<MFT>{sensors}, std::invalid_argument);
  sensors = standInPlacements<MFT>();
  sensors.back().layer = static_cast<LayerID>(Detector<MFT>::nLayers);
  BOOST_CHECK_THROW(Geometry<MFT>{sensors}, std::invalid_argument);
  // Every chip moved onto layer 0: the count agrees, but a layer is empty.
  sensors = standInPlacements<MFT>();
  for (auto& sensor : sensors) {
    sensor.layer = 0;
  }
  BOOST_CHECK_THROW(Geometry<MFT>{sensors}, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TrackerRunsOnTheIts) { runEmpty<ITS>(); }

BOOST_AUTO_TEST_CASE(TrackerRunsOnTheMft) { runEmpty<MFT>(); }
namespace
{

/// A curved sensor like the ALICE 3 VD's: a 71 degree arc of a 0.5 cm cylinder,
/// 20 um thick and 50 cm long.
constexpr SensorSpec
  kCurved{.shape = ShapeType::kCylinder, .halfX = 0.3112, .halfY = 25., .halfThickness = 0.001, .pitchX = 0.001, .pitchY = 0.001, .radius = 0.5};

/// A sensor rotated by \a angle about the beam axis and moved by \a t.
Sensor rotatedAboutBeam(double angle, std::array<double, 3> t)
{
  Sensor sensor{};
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  const double rot[9] = {c, -s, 0., s, c, 0., 0., 0., 1.};
  std::ranges::copy(rot, std::begin(sensor.rot));
  std::ranges::copy(t, std::begin(sensor.trans));
  return sensor;
}

} // namespace

BOOST_AUTO_TEST_CASE(CurvedSensorUnrollsOntoItsCylinder)
{
  // The middle of the arc is on the local +x axis, the unrolled z outwards.
  double frame[3];
  kCurved.unrolledToFrame(0., 3., 0.001, frame);
  BOOST_TEST(std::abs(frame[0] - 0.501) < kPositionTolerance);
  BOOST_TEST(std::abs(frame[1]) < kPositionTolerance);
  BOOST_TEST(std::abs(frame[2] - 3.) < kPositionTolerance);
  // The unrolled x is the arc length on the inner surface, counter-clockwise.
  kCurved.unrolledToFrame(0.5 * std::numbers::pi / 4., 0., 0., frame);
  BOOST_TEST(std::abs(frame[0] - frame[1]) < kPositionTolerance);
  BOOST_TEST(frame[1] > 0.);
  static_assert(exactlyEqual(kCurved.surfaceMin(), 0.) && exactlyEqual(kCurved.surfaceMax(), 2. * kCurved.halfThickness));
}

BOOST_AUTO_TEST_CASE(CurvedEnvelopeIsTheHullOfItsSilicon)
{
  // On the beam axis, and moved off it so far that the nearest and farthest
  // points lie inside the arc rather than at its ends.
  for (const auto& t : {std::array<double, 3>{0., 0., 0.}, std::array<double, 3>{0.03, -0.02, 1.}, std::array<double, 3>{0.2, 0.25, -2.}}) {
    const auto sensor = rotatedAboutBeam(std::numbers::pi / 4., t);
    const std::array<SensorSpec, 1> specs{kCurved};
    const auto envelope = measureEnvelopes(std::span(&sensor, 1), specs, 1).front();
    // The silicon sampled densely lies within the hull and reaches it.
    double rMin = 1e9;
    double rMax = 0.;
    for (int i = 0; i <= 2000; ++i) {
      const double x = -kCurved.halfX + (2. * kCurved.halfX * i / 2000.);
      for (const double z : {kCurved.surfaceMin(), kCurved.surfaceMax()}) {
        double frame[3];
        kCurved.unrolledToFrame(x, 0., z, frame);
        double global[3];
        sensor.toGlobal(frame, global);
        rMin = std::min(rMin, std::hypot(global[0], global[1]));
        rMax = std::max(rMax, std::hypot(global[0], global[1]));
      }
    }
    BOOST_TEST(envelope.rMin <= rMin + kPositionTolerance);
    BOOST_TEST(envelope.rMax >= rMax - kPositionTolerance);
    // 2000 steps over the arc miss an extremum by at most a micron.
    BOOST_TEST(rMin - envelope.rMin < 1e-4);
    BOOST_TEST(envelope.rMax - rMax < 1e-4);
    BOOST_TEST(std::abs(envelope.zMin - (t[2] - 25.)) < kPositionTolerance);
    BOOST_TEST(std::abs(envelope.zMax - (t[2] + 25.)) < kPositionTolerance);
  }
  // On the axis the position is the middle of the silicon.
  const auto onAxis = rotatedAboutBeam(1., {0., 0., 0.});
  const std::array<SensorSpec, 1> specs{kCurved};
  BOOST_TEST(std::abs(measureEnvelopes(std::span(&onAxis, 1), specs, 1).front().position - 0.501) < kPositionTolerance);
}

BOOST_AUTO_TEST_CASE(CurvedSensorMustBeAboutTheBeam)
{
  auto tilted = rotatedAboutBeam(0., {0., 0., 0.});
  const double c = std::cos(0.01);
  const double s = std::sin(0.01);
  const double rot[9] = {1., 0., 0., 0., c, -s, 0., s, c};
  std::ranges::copy(rot, std::begin(tilted.rot));
  const std::array<SensorSpec, 1> specs{kCurved};
  BOOST_CHECK_THROW(static_cast<void>(measureEnvelopes(std::span(&tilted, 1), specs, 1)), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TangentPlaneTouchesItsCylinder)
{
  TrackingSurface cylinder{};
  // A petal rotated by 135 degrees about a beam moved off the origin.
  const Scalar c = std::cos(0.75F * std::numbers::pi_v<Scalar>);
  const Scalar s = std::sin(0.75F * std::numbers::pi_v<Scalar>);
  cylinder.axis0[0] = c;
  cylinder.axis0[1] = s;
  cylinder.axis1[0] = -s;
  cylinder.axis1[1] = c;
  cylinder.normal[2] = 1.F;
  cylinder.origin[0] = 0.01F;
  cylinder.origin[1] = -0.02F;
  cylinder.origin[2] = 0.5F;
  cylinder.radius = 1.2F;
  for (const Scalar angle : {-0.7F, -0.1F, 0.F, 0.3F, 0.73F}) {
    const auto plane = tangentPlane(cylinder, angle);
    BOOST_TEST(plane.surfacePhi == angle);
    // (axis0, axis1, normal) is right-handed and orthonormal.
    const Scalar cross[3]{(plane.axis0[1] * plane.axis1[2]) - (plane.axis0[2] * plane.axis1[1]),
                          (plane.axis0[2] * plane.axis1[0]) - (plane.axis0[0] * plane.axis1[2]),
                          (plane.axis0[0] * plane.axis1[1]) - (plane.axis0[1] * plane.axis1[0])};
    for (int i = 0; i < 3; ++i) {
      BOOST_TEST(std::abs(cross[i] - plane.normal[i]) < 1e-6F);
    }
    // The touching point is on the cylinder at the angle, and its plane
    // coordinates are its unrolled ones.
    Scalar touch[3];
    surfaceToGlobal(plane, cylinder.radius * angle, 0.2F, touch);
    const Scalar local[3]{touch[0] - cylinder.origin[0], touch[1] - cylinder.origin[1], touch[2] - cylinder.origin[2]};
    const Scalar x = (local[0] * c) + (local[1] * s);
    const Scalar y = -(local[0] * s) + (local[1] * c);
    BOOST_TEST(std::abs(std::hypot(x, y) - cylinder.radius) < 1e-6F);
    BOOST_TEST(std::abs(std::atan2(y, x) - angle) < 1e-6F);
    BOOST_TEST(std::abs(local[2] - 0.2F) < 1e-6F);
  }
}

} // namespace o2::gtrk
