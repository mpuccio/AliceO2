// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT DetectorTraitsTrackRepresentation
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/DetectorTraits.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TimeFrame.h"

using namespace o2::itsmft::tracking;

namespace
{

o2::track::TrackParCovF makeBarrelState(float q2pt)
{
  const o2::track::TrackParCovF::params_t params{1.f, 2.f, 0.1f, 0.2f, q2pt};
  o2::track::TrackParCovF::covMat_t covariance{};
  for (size_t index = 0; index < covariance.size(); ++index) {
    covariance[index] = static_cast<float>(index + 1);
  }
  return o2::track::TrackParCovF{3.f, 0.4f, params, covariance};
}

o2::track::TrackParCovFwd makeForwardState(double invQPt)
{
  const o2::track::SMatrix5 params{1., 2., 0.3, 0.4, invQPt};
  o2::track::SMatrix55Sym covariance{};
  for (size_t row = 0; row < 5; ++row) {
    for (size_t column = 0; column < 5; ++column) {
      covariance(row, column) = 10. * row + column;
    }
  }
  return o2::track::TrackParCovFwd{5., params, covariance, 6.};
}

void checkBarrelStateEqual(const o2::track::TrackParCovF& lhs, const o2::track::TrackParCovF& rhs)
{
  BOOST_CHECK_EQUAL(lhs.getX(), rhs.getX());
  BOOST_CHECK_EQUAL(lhs.getAlpha(), rhs.getAlpha());
  for (int parameter = 0; parameter < 5; ++parameter) {
    BOOST_CHECK_EQUAL(lhs.getParam(parameter), rhs.getParam(parameter));
  }
  const auto& lhsCovariance = lhs.getCov();
  const auto& rhsCovariance = rhs.getCov();
  for (size_t element = 0; element < lhsCovariance.size(); ++element) {
    BOOST_CHECK_EQUAL(lhsCovariance[element], rhsCovariance[element]);
  }
  BOOST_CHECK_EQUAL(lhs.getUserField(), rhs.getUserField());
}

void checkForwardStateEqual(const o2::track::TrackParCovFwd& lhs, const o2::track::TrackParCovFwd& rhs)
{
  BOOST_CHECK_EQUAL(lhs.getZ(), rhs.getZ());
  for (size_t parameter = 0; parameter < 5; ++parameter) {
    BOOST_CHECK_EQUAL(lhs.getParameters()(parameter), rhs.getParameters()(parameter));
  }
  BOOST_CHECK_EQUAL(lhs.getTrackChi2(), rhs.getTrackChi2());
  const auto& lhsCovariance = lhs.getCovariances();
  const auto& rhsCovariance = rhs.getCovariances();
  for (size_t row = 0; row < 5; ++row) {
    for (size_t column = 0; column < 5; ++column) {
      BOOST_CHECK_EQUAL(lhsCovariance(row, column), rhsCovariance(row, column));
    }
  }
}

void checkITSTrackEqual(const o2::its::TrackITSExt& lhs, const o2::its::TrackITSExt& rhs)
{
  checkBarrelStateEqual(lhs, rhs);
  checkBarrelStateEqual(lhs.getParamOut(), rhs.getParamOut());
  BOOST_CHECK_EQUAL(lhs.getChi2(), rhs.getChi2());
  BOOST_CHECK_EQUAL(lhs.getPattern(), rhs.getPattern());
  BOOST_CHECK_EQUAL(lhs.getClusterSizes(), rhs.getClusterSizes());
  BOOST_CHECK_EQUAL(lhs.getTimeStamp().getTimeStamp(), rhs.getTimeStamp().getTimeStamp());
  BOOST_CHECK_EQUAL(lhs.getTimeStamp().getTimeStampError(), rhs.getTimeStamp().getTimeStampError());
  for (int layer = 0; layer < o2::its::TrackITSExt::MaxClusters; ++layer) {
    BOOST_CHECK_EQUAL(lhs.getClusterIndex(layer), rhs.getClusterIndex(layer));
  }
}

void checkMFTTrackEqual(const MFTCATrack& lhs, const MFTCATrack& rhs)
{
  checkForwardStateEqual(lhs.getTrack(), rhs.getTrack());
  BOOST_CHECK_EQUAL(lhs.getPattern(), rhs.getPattern());
  BOOST_CHECK_EQUAL(lhs.getSeedPattern(), rhs.getSeedPattern());
  BOOST_CHECK_EQUAL(lhs.hasSharedClusters(), rhs.hasSharedClusters());
  BOOST_CHECK_EQUAL(lhs.getTimeStamp().getTimeStamp(), rhs.getTimeStamp().getTimeStamp());
  BOOST_CHECK_EQUAL(lhs.getTimeStamp().getTimeStampError(), rhs.getTimeStamp().getTimeStampError());
  for (int layer = 0; layer < MFTCATrack::MaxClusters; ++layer) {
    BOOST_CHECK_EQUAL(lhs.getClusterIndex(layer), rhs.getClusterIndex(layer));
    BOOST_CHECK_EQUAL(lhs.getClusterSize(layer), rhs.getClusterSize(layer));
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(copy_seed_pattern_to_track)
{
  TrackSeed mftSeed;
  mftSeed.setSurfaceMask(SurfaceMask{0x02a5});
  DetectorTraits<10>::TrackType mftTrack;
  DetectorTraits<10>::copySeedPatternToTrack(mftTrack, mftSeed);
  BOOST_CHECK_EQUAL(mftTrack.getSeedPattern(), 0x02a5);

  TrackSeed itsSeed;
  itsSeed.setSurfaceMask(SurfaceMask{0x007f});
  DetectorTraits<7>::TrackType itsTrack{makeBarrelState(0.5f)};
  itsTrack.getParamOut() = makeBarrelState(-0.25f);
  itsTrack.setChi2(7.f);
  itsTrack.setPattern(0x1234567u);
  itsTrack.setClusterSize(0, 3);
  itsTrack.setExternalClusterIndex(0, 42, true);
  const auto before = itsTrack;
  DetectorTraits<7>::copySeedPatternToTrack(itsTrack, itsSeed);
  checkITSTrackEqual(itsTrack, before);
}

BOOST_AUTO_TEST_CASE(clear_transient_layer_pattern)
{
  DetectorTraits<7>::TrackType itsTrack{makeBarrelState(0.5f)};
  itsTrack.getParamOut() = makeBarrelState(-0.25f);
  itsTrack.setExtendedLayerPattern<7>(0x55);
  itsTrack.getParamOut().setUserField(0xaa);
  const auto equivalentTrack = itsTrack;
  auto expectedTrack = equivalentTrack;
  expectedTrack.clearExtendedLayerPattern();
  DetectorTraits<7>::clearTransientLayerPattern(itsTrack);
  checkITSTrackEqual(itsTrack, expectedTrack);
  BOOST_CHECK_EQUAL(itsTrack.getUserField(), 0);
  BOOST_CHECK_EQUAL(itsTrack.getParamOut().getUserField(), 0);

  DetectorTraits<10>::TrackType mftTrack;
  mftTrack.getTrack() = o2::mft::TrackMFT{};
  static_cast<o2::track::TrackParCovFwd&>(mftTrack.getTrack()) = makeForwardState(-0.5);
  mftTrack.setSeedPattern(0x01f3);
  mftTrack.setClusterIndex(3, 17);
  mftTrack.setClusterSize(3, 12);
  mftTrack.setSharedClusters();
  const auto before = mftTrack;
  DetectorTraits<10>::clearTransientLayerPattern(mftTrack);
  checkMFTTrackEqual(mftTrack, before);
}

BOOST_AUTO_TEST_CASE(have_same_polarity)
{
  DetectorTraits<7>::TrackType itsPositive{makeBarrelState(0.5f)};
  DetectorTraits<7>::TrackType itsSame{makeBarrelState(0.25f)};
  DetectorTraits<7>::TrackType itsOpposite{makeBarrelState(-0.25f)};
  BOOST_CHECK(DetectorTraits<7>::haveSamePolarity(itsPositive, itsSame));
  BOOST_CHECK(!DetectorTraits<7>::haveSamePolarity(itsPositive, itsOpposite));

  DetectorTraits<10>::TrackType mftPositive;
  DetectorTraits<10>::TrackType mftSame;
  DetectorTraits<10>::TrackType mftOpposite;
  static_cast<o2::track::TrackParCovFwd&>(mftPositive.getTrack()) = makeForwardState(0.5);
  static_cast<o2::track::TrackParCovFwd&>(mftSame.getTrack()) = makeForwardState(0.25);
  static_cast<o2::track::TrackParCovFwd&>(mftOpposite.getTrack()) = makeForwardState(-0.25);
  BOOST_CHECK(DetectorTraits<10>::haveSamePolarity(mftPositive, mftSame));
  BOOST_CHECK(!DetectorTraits<10>::haveSamePolarity(mftPositive, mftOpposite));
}

// --- M5d: single-shot barrel refitter boundary ------------------------------
//
// DetectorTraits<7>::refitSeed's barrel/ITS branch (refitSeedITS,
// DetectorTraits.cxx) now goes through the shared native driver
// (fitTrackSeedLegs, NativeRefitDriver.h) instead of the frozen legacy
// o2::its::track::refitTrackSeed chain -- see
// doc/decisions/0008-native-refit-activation.md. This test proves that an
// invalid (StateFamily::Invalid) seed state -- a genuinely corrupted/
// caller-bug state, distinct from a merely cross-family-convertible Forward
// state, which Propagator::convertFamily now legitimately supports -- makes
// the whole refit fail cleanly (return false, track untouched) via
// makeLinearizationReference's own hasRecognizedFamily() check, before any
// per-layer measurement is ever indexed (so an intentionally empty
// layerMeasurements/surfaceCatalog is safe to pass here).
BOOST_AUTO_TEST_CASE(RefitSeedFailsCleanlyForInvalidFamilySeedState)
{
  TrackSeed seed;
  seed.state().family = StateFamily::Invalid; // corrupted/unrecognized family
  seed.state().parameters[4] = 0.2f;
  seed.setSurfaceMask(SurfaceMask{0x007f});

  DetectorTraits<7>::TrackType track;
  o2::itsmft::TrackingParameters params;
  SurfaceTrackingScratch tf;
  const std::vector<gsl::span<const SurfaceMeasurement>> layerMeasurements(7);
  const SurfaceCatalogView surfaceCatalog{};

  const bool refitSuccess = DetectorTraits<7>::refitSeed(seed, track, params, 0.5f, tf, layerMeasurements, surfaceCatalog, ClusterSourceId{0});
  BOOST_CHECK(!refitSuccess);
}
