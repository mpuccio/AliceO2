// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#define BOOST_TEST_MODULE ITSMFT DetectorAdapterCompatibility
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/detail/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/detail/DetectorRefitSupport.h"
#include "ITSMFTTracking/detail/SurfaceKinematicStateLegacyAdapters.h"

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

struct MixedRefitFixture {
  SurfaceTrackingScratch scratch;
  std::array<std::vector<SurfaceMeasurement>, 3> measurementsStorage;
  std::array<std::vector<GlobalMeasurement>, 3> globalsStorage;
  std::vector<gsl::span<const SurfaceMeasurement>> measurements{3};
  std::vector<gsl::span<const GlobalMeasurement>> globals{3};
  std::array<SurfaceDescriptor, 3> surfaces{};
  std::array<SurfaceId, 3> ordered{SurfaceId{0}, SurfaceId{1}, SurfaceId{2}};
  TrackSeed seed{};
  o2::itsmft::TrackingParameters parameters{};

  MixedRefitFixture()
  {
    scratch.setMemoryPool(std::make_shared<o2::its::BoundedMemoryResource>());
    scratch.adoptPlan(3, 0, 0);
    const std::array<ClusterSourceId, 3> sources{ClusterSourceId{0}, ClusterSourceId{1}, ClusterSourceId{0}};
    BOOST_REQUIRE(scratch.setSurfaceSources(sources));
    for (int position = 0; position < 3; ++position) {
      surfaces[position].id = ordered[position];
      surfaces[position].kind = position == 1 ? SurfaceKind::Disk : SurfaceKind::Cylinder;
      surfaces[position].referenceCoordinate = position == 1 ? -10.f : 2.5f;
      surfaces[position].material = NominalSurfaceMaterial{0.f, 0.f};
      SurfaceMeasurement measurement{};
      measurement.frame.q = surfaces[position].referenceCoordinate;
      measurement.frame.frameAngle = 0.3f;
      measurement.frame.u = position == 1 ? 5.f : 0.8f;
      measurement.frame.v = position == 1 ? -5.f : -0.45f;
      measurement.covariance = {10.f, 0.f, 10.f};
      GlobalMeasurement global{};
      global.position = {measurement.frame.u, measurement.frame.v, measurement.frame.q};
      global.radius = std::hypot(global.position.x, global.position.y);
      global.covariance = {10.f, 0.f, 0.f, 10.f, 0.f, 10.f};
      global.cluster = ClusterRef{sources[position], static_cast<uint32_t>(100 + position)};
      global.surface = ordered[position];
      measurementsStorage[position].push_back(measurement);
      globalsStorage[position].push_back(global);
      measurements[position] = measurementsStorage[position];
      globals[position] = globalsStorage[position];
      scratch.addClusterExternalIndexToLayer(position, 100 + position);
      seed.getClusters()[position] = 0;
    }
    seed.setSurfaceMask(SurfaceMask{0x7});
    auto& state = seed.state();
    state.parameters[0] = 1.25f;
    state.parameters[1] = -0.75f;
    state.parameters[2] = 0.2f;
    state.parameters[3] = -0.35f;
    state.parameters[4] = 0.1f;
    state.referenceCoordinate = 4.f;
    state.alpha = 0.3f;
    state.kind = SurfaceKind::Cylinder;
    state.absCharge = 1;
    for (uint8_t i = 0; i < 5; ++i) {
      state.covariance[packedCovarianceIndex(i, i)] = 0.1f;
    }
    parameters.MinTrackLength = 3;
    parameters.MinPt.assign(4, 0.f);
    parameters.MaxChi2ClusterAttachment = 1.e8f;
    parameters.MaxChi2NDF = 1.e8f;
  }

  SurfaceCatalogView catalog() { return {surfaces.data(), static_cast<uint32_t>(surfaces.size())}; }
};

} // namespace

BOOST_AUTO_TEST_CASE(generic_candidate_kinematics_are_adapter_edge_only)
{
  TrackingCandidate itsCandidate;
  BOOST_REQUIRE(legacy::importBarrelTrackParCov(makeBarrelState(0.5f), itsCandidate.track.innerState));
  BOOST_REQUIRE(detail::fillCandidateKinematics(itsCandidate));
  BOOST_CHECK_EQUAL(itsCandidate.charge, makeBarrelState(0.5f).getSign());
  BOOST_CHECK_EQUAL(itsCandidate.phi, makeBarrelState(0.5f).getPhi());
  BOOST_CHECK_CLOSE(itsCandidate.eta, makeBarrelState(0.5f).getEta(), 1e-4);

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
  first.genericTrackIndex = 4;
  TrackingCandidate second;
  second.genericTrackIndex = 9;

  std::vector<TrackingCandidate> results{first, second};
  ITSSharedClusterCompatibility sidecar;
  DetectorPublicationAdapter<ITSNLayers> adapter;
  adapter.adoptITSSharedClusterCompatibility(&sidecar);
  o2::itsmft::TrackingParameters params;
  SurfaceTrackingScratch scratch;
  BOOST_REQUIRE(adapter.completeAccepted(results, params, scratch, true));
  BOOST_REQUIRE_EQUAL(sidecar.entries().size(), 2);
  BOOST_CHECK_EQUAL(sidecar.entries()[0].genericTrackIndex, 4);
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
  result.genericTrackIndex = 7;
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
  BOOST_CHECK_EQUAL(sidecar.entries()[0].genericTrackIndex, 7);
  BOOST_CHECK_EQUAL(sidecar.entries()[0].invQPtSeed, -0.25);
  BOOST_CHECK_EQUAL(sidecar.entries()[0].chi2QPtSeed, 0.);
  BOOST_CHECK_EQUAL(sidecar.entries()[0].seedPattern, 0x02a5);
}

BOOST_AUTO_TEST_CASE(native_refit_rejects_invalid_generic_state)
{
  TrackSeed seed;
  seed.state().kind = SurfaceKind::Undefined;
  seed.state().parameters[4] = 0.2f;
  seed.setSurfaceMask(SurfaceMask{0x007f});

  TrackingCandidate candidate;
  o2::itsmft::TrackingParameters params;
  const std::vector<gsl::span<const GlobalMeasurement>> layerGlobals(7);
  const std::vector<gsl::span<const SurfaceMeasurement>> layerMeasurements(7);
  const SurfaceCatalogView surfaceCatalog{};
  SurfaceTrackingScratch scratch;
  const std::vector<SurfaceId> orderedSurfaces(7);
  BOOST_CHECK(!detail::refitSurfaceSeed(seed, params, 0.5f, scratch, layerGlobals, layerMeasurements, surfaceCatalog, orderedSurfaces, candidate));
}

BOOST_AUTO_TEST_CASE(generic_refit_accepts_mixed_family_and_source_seed)
{
  MixedRefitFixture fixture;
  TrackingCandidate candidate;
  BOOST_REQUIRE(detail::refitSurfaceSeed(fixture.seed, fixture.parameters, 0.5f, fixture.scratch,
                                         fixture.globals, fixture.measurements, fixture.catalog(), fixture.ordered, candidate));
  BOOST_CHECK(candidate.track.innerState.hasRecognizedKind());
  BOOST_CHECK(candidate.track.outerState.hasRecognizedKind());
  BOOST_CHECK_EQUAL(candidate.getNumberOfClusters(), 3);
}

BOOST_AUTO_TEST_CASE(generic_refit_rejects_foreign_per_surface_source_before_commit)
{
  MixedRefitFixture fixture;
  fixture.globalsStorage[1][0].cluster.source = ClusterSourceId{0};
  fixture.globals[1] = fixture.globalsStorage[1];
  TrackingCandidate before{};
  TrackingCandidate candidate = before;
  BOOST_CHECK(!detail::refitSurfaceSeed(fixture.seed, fixture.parameters, 0.5f, fixture.scratch,
                                        fixture.globals, fixture.measurements, fixture.catalog(), fixture.ordered, candidate));
  BOOST_CHECK_EQUAL(candidate.track.chi2, before.track.chi2);
  BOOST_CHECK_EQUAL(candidate.seed.getSurfaceMask().value(), before.seed.getSurfaceMask().value());
}

BOOST_AUTO_TEST_CASE(generic_refit_validates_each_measurement_before_commit)
{
  const auto rejects = [](MixedRefitFixture& fixture) {
    TrackingCandidate candidate;
    const TrackingCandidate before = candidate;
    BOOST_CHECK(!detail::refitSurfaceSeed(fixture.seed, fixture.parameters, 0.5f, fixture.scratch,
                                          fixture.globals, fixture.measurements, fixture.catalog(), fixture.ordered, candidate));
    BOOST_CHECK_EQUAL(candidate.track.chi2, before.track.chi2);
    BOOST_CHECK_EQUAL(candidate.seed.getSurfaceMask().value(), before.seed.getSurfaceMask().value());
  };

  {
    MixedRefitFixture fixture;
    fixture.globalsStorage[1][0].surface = SurfaceId{0};
    fixture.globals[1] = fixture.globalsStorage[1];
    rejects(fixture);
  }
  {
    MixedRefitFixture fixture;
    fixture.measurementsStorage[1][0].covariance.uu = std::numeric_limits<float>::quiet_NaN();
    fixture.measurements[1] = fixture.measurementsStorage[1];
    rejects(fixture);
  }
  {
    MixedRefitFixture fixture;
    fixture.seed.getClusters()[1] = 1;
    rejects(fixture);
  }
  {
    MixedRefitFixture fixture;
    SurfaceTrackingScratch missingMapping;
    missingMapping.setMemoryPool(std::make_shared<o2::its::BoundedMemoryResource>());
    missingMapping.adoptPlan(3, 0, 0);
    for (int position = 0; position < 3; ++position) {
      missingMapping.addClusterExternalIndexToLayer(position, 100 + position);
    }
    TrackingCandidate candidate;
    BOOST_CHECK(!detail::refitSurfaceSeed(fixture.seed, fixture.parameters, 0.5f, missingMapping,
                                          fixture.globals, fixture.measurements, fixture.catalog(), fixture.ordered, candidate));
  }
}

BOOST_AUTO_TEST_CASE(common_refit_boundary_has_no_first_hit_detector_dispatch)
{
  namespace fs = std::filesystem;
  const auto testPath = fs::path{__FILE__};
  const auto header = testPath.parent_path().parent_path() / "include/ITSMFTTracking/detail/DetectorRefitSupport.h";
  std::ifstream input{header};
  BOOST_REQUIRE(input);
  std::ostringstream contents;
  contents << input.rdbuf();
  const auto source = contents.str();
  const auto begin = source.find("inline bool refitSurfaceSeed");
  BOOST_REQUIRE(begin != std::string::npos);
  const auto body = source.substr(begin);
  BOOST_CHECK(body.find("refitITSSeed") == std::string::npos);
  BOOST_CHECK(body.find("refitMFTSeed") == std::string::npos);
  BOOST_CHECK(body.find("getSurface(surface).kind") == std::string::npos);
}
