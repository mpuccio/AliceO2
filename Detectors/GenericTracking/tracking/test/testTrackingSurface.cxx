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

/// \file testTrackingSurface.cxx
/// \brief The host-to-tracking boundary: packed covariance indexing, and the
///        surface records the propagator will be written against.
///
/// The load-bearing test here is the one against Geometry<Det>::toGlobal.  A
/// surface built from the *rows* of Sensor::rot instead of its columns is a
/// transpose that no round trip through TrackingSurface itself can detect,
/// because both directions would use the same wrong matrix.  Only an
/// independent implementation of the same transform sees it.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TrackingSurface
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Base/Geometry.h"
#include "ToyBarrel4Disc8.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/ShapeCut.h"

#include <cmath>
#include <limits>
#include <set>

namespace o2::gtrk
{

using Mixed = detectors::ToyBarrel4Disc8;

/// The position budget: 1 um, against a 14 um single-pixel cluster
/// RMS.  It has to hold at |z| = 152 cm, where a float ulp is already
/// 1.2e-5 cm, so this is a real constraint on the narrowing and not a
/// formality.
constexpr double kMaxPosError = 1e-4; // cm

BOOST_AUTO_TEST_CASE(PackedCovarianceIndexingIsUnambiguous)
{
  static_assert(kNBoundParams == 5);
  static_assert(kBoundCovSize == 15);

  std::set<int> seen;
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j <= i; ++j) {
      const int idx = covIndex(i, j);
      BOOST_TEST(idx >= 0);
      BOOST_TEST(idx < kBoundCovSize);
      // Symmetric in its arguments: this is the property that keeps a caller
      // from having to order the indices, which is the usual way a packed
      // covariance goes quietly wrong.
      BOOST_TEST(covIndex(j, i) == idx);
      seen.insert(idx);
    }
  }
  // Every slot used exactly once -- no collisions, no gaps.
  BOOST_TEST(seen.size() == static_cast<std::size_t>(kBoundCovSize));

  BOOST_TEST(varIndex(kLoc0) == 0);
  BOOST_TEST(varIndex(kQOverP) == kBoundCovSize - 1);
}

BOOST_AUTO_TEST_CASE(SurfaceAgreesWithTheGeometryItWasBuiltFrom)
{
  const auto geo = toyGeometry<Mixed>();
  const TrackingGeometry<Mixed> tgeo(geo);

  // One barrel sensor and one disc sensor: the two orientations the single
  // bound model has to serve.
  const SensorID barrelId = 0;
  const SensorID discId = firstSensorOfLayer<Mixed>(4); // first disc layer

  for (const SensorID id : {barrelId, discId}) {
    const TrackingSurface& s = tgeo.surface(id);

    for (const double loc0 : {-0.3, 0.0, 0.5}) {
      for (const double loc1 : {-1.2, 0.0, 2.0}) {
        const double local[3] = {loc0, loc1, 0.};
        double expected[3] = {};
        geo.toGlobal(id, local, expected);

        Scalar got[3] = {};
        surfaceToGlobal(s, static_cast<Scalar>(loc0), static_cast<Scalar>(loc1), got);

        for (int k = 0; k < 3; ++k) {
          BOOST_TEST(std::abs(static_cast<double>(got[k]) - expected[k]) < kMaxPosError);
        }

        // ...and back, including the out-of-plane component, which must be
        // zero for a point that was constructed on the plane.
        Scalar back[3] = {};
        surfaceToLocal(s, got, back);
        BOOST_TEST(std::abs(static_cast<double>(back[0]) - loc0) < kMaxPosError);
        BOOST_TEST(std::abs(static_cast<double>(back[1]) - loc1) < kMaxPosError);
        BOOST_TEST(std::abs(static_cast<double>(back[2])) < kMaxPosError);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(AxesAreOrthonormalAndRightHanded)
{
  const auto geo = toyGeometry<Mixed>();
  const TrackingGeometry<Mixed> tgeo(geo);

  for (SensorID id = 0; id < static_cast<SensorID>(Mixed::nSensors); id += 97) {
    const TrackingSurface& s = tgeo.surface(id);
    const Scalar* ax[3] = {s.axis0, s.axis1, s.normal};

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const Scalar dot = ax[i][0] * ax[j][0] + ax[i][1] * ax[j][1] + ax[i][2] * ax[j][2];
        BOOST_TEST(std::abs(dot - (i == j ? 1.F : 0.F)) < 1e-5F);
      }
    }

    // normal == axis0 x axis1, so the frame is right-handed and the sign of
    // the normal is not a free choice.
    const Scalar cx = s.axis0[1] * s.axis1[2] - s.axis0[2] * s.axis1[1];
    const Scalar cy = s.axis0[2] * s.axis1[0] - s.axis0[0] * s.axis1[2];
    const Scalar cz = s.axis0[0] * s.axis1[1] - s.axis0[1] * s.axis1[0];
    BOOST_TEST(std::abs(cx - s.normal[0]) < 1e-5F);
    BOOST_TEST(std::abs(cy - s.normal[1]) < 1e-5F);
    BOOST_TEST(std::abs(cz - s.normal[2]) < 1e-5F);
  }
}

BOOST_AUTO_TEST_CASE(NormalsFallIntoTheTwoFamiliesTheClosedFormIntersectionNeeds)
{
  // The propagator's two closed forms are only legitimate because every
  // sensor normal is either purely transverse (barrel: the plane is a line in
  // the transverse projection, so the intersection is a quadratic) or exactly
  // along z (disc: z is linear in path length). This pins that precondition.
  const auto geo = toyGeometry<Mixed>();
  const TrackingGeometry<Mixed> tgeo(geo);

  int barrel = 0;
  int disc = 0;
  for (SensorID id = 0; id < static_cast<SensorID>(Mixed::nSensors); ++id) {
    const TrackingSurface& s = tgeo.surface(id);
    if (s.layer < Detector<Mixed>::nBarrelLayers) {
      BOOST_TEST(std::abs(s.normal[2]) < 1e-6F);
      ++barrel;
    } else {
      BOOST_TEST(std::abs(s.normal[0]) < 1e-6F);
      BOOST_TEST(std::abs(s.normal[1]) < 1e-6F);
      BOOST_TEST(std::abs(std::abs(s.normal[2]) - 1.F) < 1e-6F);
      ++disc;
    }
  }
  BOOST_TEST(barrel > 0);
  BOOST_TEST(disc > 0);
  BOOST_TEST(barrel + disc == Mixed::nSensors);
}

BOOST_AUTO_TEST_CASE(BoundsAreARectangleThatRejectsNaN)
{
  const SurfaceBounds b{.halfX = 0.85F, .halfY = 3.4F};
  BOOST_TEST(b.isInside(0.F, 0.F));
  BOOST_TEST(b.isInside(0.84F, 3.3F));
  BOOST_TEST(!b.isInside(0.9F, 0.F));
  BOOST_TEST(!b.isInside(0.F, 3.5F));
  BOOST_TEST(!b.isInside(std::numeric_limits<Scalar>::quiet_NaN(), 0.F));
  BOOST_TEST(!b.isInside(0.F, std::numeric_limits<Scalar>::quiet_NaN()));
}

BOOST_AUTO_TEST_CASE(ALineIsHeldToTheSlopeAlongTheNormal)
{
  BOOST_TEST(slopeAdmits(1.F, 3.F, 2.F, 0.F, 2.F, 1.5F, 0.7F, -3.F));
  BOOST_TEST(slopeAdmits(2.F, 3.F, 2.F, 0.F, 2.F, 1.5F, 0.7F, 3.F));
  BOOST_TEST(!slopeAdmits(1.F, 2.F, 2.F, 0.F, 2.F, 1.5F, 0.7F, 3.F));
  BOOST_TEST(!slopeAdmits(std::numeric_limits<Scalar>::infinity(), 0.F, 2.F, 0.F, 2.F, 1.5F, 0.F, 3.F));
  BOOST_TEST(slopeAdmits(0.F, kSlopeUnbounded, 2.F, 0.F, 2.F, 1.5F, 0.F, 300.F));
}

} // namespace o2::gtrk
