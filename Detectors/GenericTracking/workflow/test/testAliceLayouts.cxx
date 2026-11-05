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

/// \file testAliceLayouts.cxx
/// \brief O2 sensor placements against L2G matrices and ideal envelopes.
/// Set GTRK_ALICE_GEOMETRY and GTRK_ALICE_GEOMETRY_IDEAL for both checks.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK ALICE layouts
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Base/LayerEnvelopes.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"

#include "DetectorsBase/GeometryManager.h"
#include "ITSBase/GeometryTGeo.h"
#include "MFTBase/GeometryTGeo.h"
#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk
{

namespace
{

using detectors::ITS;
using detectors::MFT;

/// Our sensor frame against O2's L2G: ALICE's local (x, 0, z) is our
/// (x - kActiveCentreX, z, 0) on the same chip.
template <typename GeometryTGeo>
void checkAgainstL2G(GeometryTGeo& o2Geometry, std::span<const Sensor> sensors)
{
  o2Geometry.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  double worst = 0.;
  for (std::size_t chip = 0; chip < sensors.size(); chip += 7U) {
    for (const double xa : {-0.6, 0.0581, 0.7}) {
      for (const double za : {-1.4, 0., 1.4}) {
        const auto expected = o2Geometry.getMatrixL2G(static_cast<int>(chip))(o2::math_utils::Point3D<double>{xa, 0., za});
        const double local[3] = {xa - alice::kActiveCentreX, za, 0.};
        double global[3];
        sensors[chip].toGlobal(local, global);
        worst = std::max({worst, std::abs(global[0] - expected.X()), std::abs(global[1] - expected.Y()), std::abs(global[2] - expected.Z())});
      }
    }
    BOOST_TEST(static_cast<int>(sensors[chip].layer) == o2Geometry.getLayer(static_cast<int>(chip)));
  }
  BOOST_TEST(worst < 1e-9);
}

/// The reference envelopes of the ITS and the MFT, measured by
/// o2-gtrk-import-geometry from an ideal o2sim_geometry-aligned.root. The
/// measurement at load must give them back to the bit on that geometry.
// Measured values that happen to lie near mathematical constants are not them.
// NOLINTBEGIN(modernize-use-std-numbers)
constexpr std::array<LayerEnvelope, 7> kReferenceIts{{
  {.position = 2.3260, .rMin = 2.2331, .rMax = 2.6036, .zMin = -13.5571, .zMax = 13.5571},
  {.position = 3.1354, .rMin = 3.0066, .rMax = 3.4012, .zMin = -13.5571, .zMax = 13.5571},
  {.position = 3.9162, .rMin = 3.7820, .rMax = 4.1614, .zMin = -13.5571, .zMax = 13.5571},
  {.position = 19.5882, .rMin = 19.3518, .rMax = 19.8957, .zMin = -42.2071, .zMax = 42.2071},
  {.position = 24.5272, .rMin = 24.3018, .rMax = 24.8096, .zMin = -42.2071, .zMax = 42.2071},
  {.position = 34.3546, .rMin = 34.1418, .rMax = 34.6084, .zMin = -73.8721, .zMax = 73.8721},
  {.position = 39.3106, .rMin = 39.1018, .rMax = 39.5554, .zMin = -73.8721, .zMax = 73.8721},
}};
constexpr std::array<LayerEnvelope, 10> kReferenceMft{{
  {.position = -45.2889, .rMin = 2.3029, .rMax = 12.2259, .zMin = -45.2903, .zMax = -45.2874},
  {.position = -46.7111, .rMin = 2.3029, .rMax = 12.2259, .zMin = -46.7126, .zMax = -46.7097},
  {.position = -48.5889, .rMin = 2.3029, .rMax = 12.2259, .zMin = -48.5903, .zMax = -48.5874},
  {.position = -50.0111, .rMin = 2.3029, .rMax = 12.2259, .zMin = -50.0126, .zMax = -50.0097},
  {.position = -52.3889, .rMin = 2.3029, .rMax = 13.6984, .zMin = -52.3903, .zMax = -52.3874},
  {.position = -53.8111, .rMin = 2.3029, .rMax = 13.6984, .zMin = -53.8126, .zMax = -53.8097},
  {.position = -67.6889, .rMin = 3.4029, .rMax = 16.7089, .zMin = -67.6903, .zMax = -67.6874},
  {.position = -69.1111, .rMin = 3.4029, .rMax = 16.7089, .zMin = -69.1126, .zMax = -69.1097},
  {.position = -76.0889, .rMin = 3.8029, .rMax = 17.3105, .zMin = -76.0903, .zMax = -76.0874},
  {.position = -77.5111, .rMin = 3.8029, .rMax = 17.3105, .zMin = -77.5126, .zMax = -77.5097},
}};
// NOLINTEND(modernize-use-std-numbers)

template <DetectorLayout Det, std::size_t N>
void checkReference(const std::array<LayerEnvelope, N>& reference)
{
  static_assert(N == Detector<Det>::nLayers);
  for (std::size_t layer = 0; layer < N; ++layer) {
    const auto& e = PlacedEnvelopes<Det>::all()[layer];
    BOOST_TEST(e.position == reference[layer].position);
    BOOST_TEST(e.rMin == reference[layer].rMin);
    BOOST_TEST(e.rMax == reference[layer].rMax);
    BOOST_TEST(e.zMin == reference[layer].zMin);
    BOOST_TEST(e.zMax == reference[layer].zMax);
  }
}

const char* geometryFile() { return std::getenv("GTRK_ALICE_GEOMETRY"); }

} // namespace

BOOST_AUTO_TEST_CASE(PlacementsFollowO2)
{
  if (geometryFile() == nullptr) {
    BOOST_TEST_MESSAGE("GTRK_ALICE_GEOMETRY not set: O2's placements are not checked");
    return;
  }
  o2::base::GeometryManager::loadGeometry(geometryFile());
  const auto its = alice::placeChecked<ITS>();
  const auto mft = alice::placeChecked<MFT>();
  checkAgainstL2G(*o2::its::GeometryTGeo::Instance(), its);
  checkAgainstL2G(*o2::mft::GeometryTGeo::Instance(), mft);

  // The normal convention the frame mapping gives: away from the beam on the
  // ITS, alternating along z on the MFT.
  for (const auto& sensor : its) {
    const double outward = (sensor.rot[2] * sensor.trans[0]) + (sensor.rot[5] * sensor.trans[1]);
    BOOST_TEST_REQUIRE(outward > 0.);
  }
  for (const auto& sensor : mft) {
    BOOST_TEST_REQUIRE(std::abs(sensor.rot[8]) > 0.999999);
  }

  // The envelopes measured at load, against the reference tables. Only the
  // ideal geometry reproduces them; an aligned one does not, and is not meant
  // to.
  if (std::getenv("GTRK_ALICE_GEOMETRY_IDEAL") != nullptr) {
    [[maybe_unused]] const Geometry<ITS> itsGeometry{its};
    [[maybe_unused]] const Geometry<MFT> mftGeometry{mft};
    checkReference<ITS>(kReferenceIts);
    checkReference<MFT>(kReferenceMft);
  } else {
    BOOST_TEST_MESSAGE("GTRK_ALICE_GEOMETRY_IDEAL not set: the measured envelopes are not held to the reference tables");
  }
}

} // namespace o2::gtrk
