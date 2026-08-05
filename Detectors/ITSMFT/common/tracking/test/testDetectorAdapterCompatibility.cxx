// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#define BOOST_TEST_MODULE ITSMFT DetectorAdapterCompatibility
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <vector>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/DetectorTrackingOperationAdapterSupport.h"
#include "ITSMFTTracking/SurfaceKinematicStateLegacyAdapters.h"

using namespace o2::itsmft::tracking;

namespace
{

o2::track::TrackParCovF makeBarrelState(float q2pt)
{
  const o2::track::TrackParCovF::params_t params{1.f, 2.f, 0.1f, 0.2f, q2pt};
  o2::track::TrackParCovF::covMat_t covariance{};
  return o2::track::TrackParCovF{3.f, 0.4f, params, covariance};
}

o2::track::TrackParCovFwd makeForwardState(double invQPt)
{
  const o2::track::SMatrix5 params{1., 2., 0.3, 0.4, invQPt};
  o2::track::SMatrix55Sym covariance{};
  return o2::track::TrackParCovFwd{5., params, covariance, 6.};
}

} // namespace

BOOST_AUTO_TEST_CASE(generic_candidate_kinematics_are_adapter_edge_only)
{
  TrackingCandidate itsCandidate;
  BOOST_REQUIRE(legacy::importBarrelTrackParCov(makeBarrelState(0.5f), itsCandidate.track.innerState));
  BOOST_REQUIRE(detail::fillCandidateKinematics(itsCandidate));
  BOOST_CHECK_EQUAL(itsCandidate.charge, makeBarrelState(0.5f).getSign());
  BOOST_CHECK_EQUAL(itsCandidate.phi, makeBarrelState(0.5f).getPhi());
  BOOST_CHECK_EQUAL(itsCandidate.eta, makeBarrelState(0.5f).getEta());

  TrackingCandidate mftCandidate;
  BOOST_REQUIRE(legacy::importLegacyForwardTrackParCov(makeForwardState(-0.25), mftCandidate.track.innerState));
  BOOST_REQUIRE(detail::fillCandidateKinematics(mftCandidate));
  BOOST_CHECK_EQUAL(mftCandidate.charge, makeForwardState(-0.25).getCharge());
  BOOST_CHECK_CLOSE(mftCandidate.phi, makeForwardState(-0.25).getPhi(), 1e-4);
  BOOST_CHECK_CLOSE(mftCandidate.eta, makeForwardState(-0.25).getEta(), 1e-4);
}

BOOST_AUTO_TEST_CASE(its_compatibility_consumes_generic_results)
{
  TrackingCandidate first;
  first.commonTrackIndex = 4;
  TrackingCandidate second;
  second.commonTrackIndex = 9;

  std::vector<TrackingCandidate> results{first, second};
  ITSSharedClusterCompatibility sidecar;
  DetectorPublicationAdapter<ITSNLayers> adapter;
  adapter.adoptITSSharedClusterCompatibility(&sidecar);
  o2::itsmft::TrackingParameters params;
  SurfaceTrackingScratch scratch;
  BOOST_REQUIRE(adapter.completeAccepted(results, params, scratch, true));
  BOOST_REQUIRE_EQUAL(sidecar.entries().size(), 2);
  BOOST_CHECK_EQUAL(sidecar.entries()[0].commonTrackIndex, 4);
  BOOST_CHECK(!sidecar.entries()[0].hasSharedClusters);
  BOOST_CHECK(!sidecar.entries()[1].hasSharedClusters);

  std::vector<uint8_t> sharedFlags(10, 0);
  sharedFlags[4] = 1;
  BOOST_REQUIRE(sidecar.replaceFromAcceptedResults(results, sharedFlags));
  BOOST_CHECK(sidecar.entries()[0].hasSharedClusters);
  BOOST_CHECK(!sidecar.entries()[1].hasSharedClusters);
}

BOOST_AUTO_TEST_CASE(mft_compatibility_consumes_seed_and_common_result)
{
  TrackingCandidate result;
  result.commonTrackIndex = 7;
  result.seed.state().parameters[4] = -0.25f;
  result.seed.setSurfaceMask(SurfaceMask{0x02a5});

  std::vector<TrackingCandidate> results{result};
  MFTPublicationCompatibility sidecar;
  DetectorPublicationAdapter<MFTNLayers> adapter;
  adapter.adoptMFTPublicationCompatibility(&sidecar);
  o2::itsmft::TrackingParameters params;
  SurfaceTrackingScratch scratch;
  BOOST_REQUIRE(adapter.completeAccepted(results, params, scratch, true));
  BOOST_REQUIRE_EQUAL(sidecar.entries().size(), 1);
  BOOST_CHECK_EQUAL(sidecar.entries()[0].commonTrackIndex, 7);
  BOOST_CHECK_EQUAL(sidecar.entries()[0].invQPtSeed, -0.25);
  BOOST_CHECK_EQUAL(sidecar.entries()[0].chi2QPtSeed, 0.);
  BOOST_CHECK_EQUAL(sidecar.entries()[0].seedPattern, 0x02a5);
}

BOOST_AUTO_TEST_CASE(native_refit_rejects_invalid_generic_state)
{
  TrackSeed seed;
  seed.state().family = StateFamily::Invalid;
  seed.state().parameters[4] = 0.2f;
  seed.setSurfaceMask(SurfaceMask{0x007f});

  TrackingCandidate candidate;
  o2::itsmft::TrackingParameters params;
  const std::vector<gsl::span<const SurfaceMeasurement>> layerMeasurements(7);
  const SurfaceCatalogView surfaceCatalog{};
  BOOST_CHECK(!detail::refitITSSeed<ITSNLayers>(seed, params, 0.5f, layerMeasurements, surfaceCatalog, candidate));
}
