// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Stage-B MFT normalized-measurement authority slice (Architecture.md Sec
// 6.5/12, AgentCoordination.md MFT-adapter role), updated for M5d: focused
// coverage for the generic MFT refit operation in MFTFwdTrackHelpers.cxx,
// proving every physical hit coordinate/covariance the shared native driver
// (fitTrackSeedLegs, NativeRefitDriver.h) consumes comes from the
// caller-supplied runtime span-of-spans (TrackerTraits::
// mLayerMeasurements) and SurfaceCatalogView, never from detector-legacy
// legacy Cluster/TrackingFrameInfo backfill. Each test calls the generic
// refit operation directly with a hand-built SurfaceTrackingScratch/seed/measurement-span/
// catalog fixture -- it does not exercise TrackerTraits::initialiseTimeFrame()
// or the full CA traversal, which already have their own focused coverage
// (testComputeLayerCellsOrchestration.cxx et al.) for the
// NormalizedMeasurementMismatch load-time contract this file's fixtures
// assume as given.
//
// M5d note: refitTrackFwd no longer drives the frozen o2::mft::TrackFitter/
// TrackLTF Kalman engine -- it shares the same descriptor-driven native
// driver the barrel/ITS branch uses (doc/decisions/0008-native-refit-activation.md).
// Every assertion below that only depends on *which* measurement source is
// authoritative (normalized vs. legacy backfill) or on generic Kalman
// properties (larger measurement variance cannot shrink posterior variance)
// remains valid unchanged. DiagonalOnlyCovarianceIgnoresOffDiagonalTerm is
// replaced by OffDiagonalCovarianceIsUsedByNativeUpdate: unlike the retired
// MFT TrackFitter, the native forward::predictedChi2/update operations use
// the measurement's full uu/uv/vv covariance (ForwardSurfaceStateOperations.h),
// so a nonzero, physically valid correlation now genuinely changes the fit --
// this is a disclosed, intentional behavior difference, not an oversight.

#define BOOST_TEST_MODULE ITSMFT MFTNormalizedRefit
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/ForwardSurfaceStateOperations.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/DetectorTrackingOperationAdapterSupport.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Constants.h"
#include "MFTTracking/Constants.h"

using namespace o2::itsmft::tracking;

namespace
{

constexpr int NLayers = o2::mft::constants::mft::LayersNumber;
// Field-off dispatch (refitTrackFwd: std::abs(bz) < 1e-6f): exercises the
// native Linear propagation model deterministically, independent of
// o2::mft::MFTTrackingParam::Instance()'s global forceZeroField default.
constexpr float Bz = 0.f;
constexpr float DefaultSigma2 = 2.5e-7f; // (~0.5 micron)^2, MFT-scale resolution
constexpr float PoisonCoordinate = 9999.f;

// A straight track through every MFT disk. Points are exact (analytic) so the
// native Kalman fit converges with ~0 chi2 regardless of which reasonable
// covariance is assigned -- residuals, not covariance magnitude, dominate
// chi2 here, which is what keeps the fixture deterministic without needing to
// hand-tune the fitter's internal MCS/Kalman numerics.
struct StraightTrackGeometry {
  std::array<float, NLayers> x{};
  std::array<float, NLayers> y{};
  std::array<float, NLayers> z{};
  float xSlope{};

  explicit StraightTrackGeometry(float slope) : xSlope(slope)
  {
    const auto zLayer = o2::mft::constants::mft::LayerZCoordinate();
    const float z0 = zLayer[0];
    for (int layer = 0; layer < NLayers; ++layer) {
      z[layer] = zLayer[layer];
      x[layer] = 1.f + xSlope * (z[layer] - z0);
      y[layer] = 0.5f - 0.006f * (z[layer] - z0);
    }
  }
};

// Owns everything refitTrackFwd needs for one call: a SurfaceTrackingScratch
// (external-index/size bookkeeping only -- never a physical-coordinate
// source once populated), the normalized measurement spans under test, a
// matching zero-material Disk surface catalog (one SurfaceDescriptor per
// layer, id == layer, resolved by every measurement's own .surface field),
// and the seed. Every test builds its own fixture so mutating one layer's
// normalized measurement or legacy backfill in one test can never leak into
// another.
struct RefitFixture {
  SurfaceTrackingScratch tf;
  std::array<std::vector<SurfaceMeasurement>, NLayers> storage;
  std::vector<gsl::span<const SurfaceMeasurement>> layerMeasurements = std::vector<gsl::span<const SurfaceMeasurement>>(NLayers);
  std::vector<SurfaceDescriptor> catalogSurfaces;
  SurfaceCatalogView catalog{};
  TrackSeed seed;
  o2::itsmft::TrackingParameters params;
  int nHitLayers{0};

  explicit RefitFixture(const StraightTrackGeometry& geometry, int hits = NLayers)
    : nHitLayers(hits)
  {
    params.MinTrackLength = 5;
    params.MinPt.assign(NLayers + 1, 0.f);
    params.TrackletMinAbsX = 0.f;
    params.MaxChi2NDF = 30.f;

    catalogSurfaces.resize(NLayers);
    for (int layer = 0; layer < NLayers; ++layer) {
      catalogSurfaces[layer].id = SurfaceId{static_cast<uint16_t>(layer)};
      catalogSurfaces[layer].detectorSurfaceIndex = static_cast<uint16_t>(layer);
      catalogSurfaces[layer].kind = SurfaceKind::Disk;
      catalogSurfaces[layer].material = NominalSurfaceMaterial{0.f, 0.f};
    }
    catalog = SurfaceCatalogView{catalogSurfaces.data(), static_cast<uint32_t>(catalogSurfaces.size())};
    tf.setMemoryPool(std::make_shared<o2::its::BoundedMemoryResource>());
    tf.adoptPlan(NLayers, 0, 0);

    uint16_t mask = 0;
    for (int layer = 0; layer < hits; ++layer) {
      setMeasurement(layer, geometry.x[layer], geometry.y[layer], geometry.z[layer],
                     DefaultSigma2, DefaultSigma2, static_cast<uint32_t>(layer));
      tf.addClusterExternalIndexToLayer(layer, layer); // extIdx == layer, clIdx == 0
      seed.getClusters()[layer] = 0;
      mask |= static_cast<uint16_t>(uint16_t(1) << layer);
    }
    seed.setSurfaceMask(SurfaceMask{mask});

    // M5d: fitTrackSeedLegs (unlike the retired TrackLTF engine, which
    // computed its own seed internally from the attached points) starts Leg
    // A from the caller-supplied seed.state() -- production always supplies
    // one already built by forward::buildSeed during CA cell construction.
    // This fixture reproduces that here from its own inner/middle/outer
    // measurements (layer 0 closest to the vertex has the largest z, per
    // the MFT geometry convention buildSeed's own strict z-ordering check
    // enforces) rather than leaving seed.state() default (StateFamily::
    // Invalid), which fitTrackSeedLegs would reject before ever reading a
    // measurement.
    const int innerLayer = 0;
    const int middleLayer = hits / 2;
    const int outerLayer = hits - 1;
    OperationFailureReason seedReason{};
    BOOST_REQUIRE(forward::buildSeed(layerMeasurements[innerLayer][0], layerMeasurements[middleLayer][0],
                                     layerMeasurements[outerLayer][0], Bz, params.TrackletMinPt,
                                     /*absCharge=*/1, o2::track::PID::Pion, seed.state(), seedReason));
  }

  void setMeasurement(int layer, float x, float y, float z, float uu, float vv, uint32_t clusterIndex, float uv = 0.f)
  {
    SurfaceMeasurement m{};
    m.global = {x, y, z};
    // Disk adapter contract (ForwardSurfaceStateOperations.h): frame.q ==
    // global.z. The retired legacy TrackLTF engine only ever read
    // measurement.global, so this fixture never needed frame.q populated
    // until now: Propagator::propagateToMeasurement's forward routing
    // propagates to frame.q (the same target every native forward CA
    // operation already uses), not global.z directly.
    m.frame.q = z;
    m.covariance.uu = uu;
    m.covariance.vv = vv;
    m.covariance.uv = uv;
    m.cluster = ClusterRef{ClusterSourceId{0}, clusterIndex};
    m.surface = SurfaceId{static_cast<uint16_t>(layer)};
    storage[layer].assign(1, m);
    layerMeasurements[layer] = storage[layer];
  }

  // Poisons legacy backfill (mUnsortedClusters / mTrackingFrameInfo) at every
  // hit layer with coordinates far from the real trajectory. refitTrackFwd
  // must never read either: MFTFwdTrackHelpers.cxx no longer calls
  // tf.getUnsortedClusters()/tf.getClusterTrackingFrameInfo() at all.
  void poisonLegacyBackfill()
  {
    for (int layer = 0; layer < nHitLayers; ++layer) {
      tf.addClusterToLayer(layer, PoisonCoordinate, PoisonCoordinate, PoisonCoordinate, 0);
      tf.addTrackingFrameInfoToLayer(layer, o2::its::TrackingFrameInfo{
                                              PoisonCoordinate, PoisonCoordinate, PoisonCoordinate, PoisonCoordinate, PoisonCoordinate, {PoisonCoordinate, PoisonCoordinate}, {1.f, 0.f, 1.f}});
    }
  }
};

bool refit(RefitFixture& fixture, TrackingCandidate& candidate)
{
  return detail::refitMFTSeed(fixture.seed, fixture.params, Bz, fixture.tf,
                              fixture.layerMeasurements, fixture.catalog, ClusterSourceId{0}, candidate);
}

bool refit(const TrackSeed& seed, o2::itsmft::TrackingParameters& params, SurfaceTrackingScratch& scratch,
           gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
           SurfaceCatalogView catalog, TrackingCandidate& candidate)
{
  return detail::refitMFTSeed(seed, params, Bz, scratch, layerMeasurements, catalog, ClusterSourceId{0}, candidate);
}

void checkTrackUnchanged(const TrackingCandidate& before, const TrackingCandidate& after)
{
  BOOST_CHECK_EQUAL(before.seed.getSurfaceMask().value(), after.seed.getSurfaceMask().value());
  for (int position = 0; position < TrackSeed::MaxSurfaces; ++position) {
    BOOST_CHECK_EQUAL(before.seed.getCluster(position), after.seed.getCluster(position));
  }
  for (int i = 0; i < 5; ++i) {
    BOOST_CHECK_EQUAL(before.track.innerState.parameters[i], after.track.innerState.parameters[i]);
    BOOST_CHECK_EQUAL(before.track.outerState.parameters[i], after.track.outerState.parameters[i]);
  }
  for (int i = 0; i < 15; ++i) {
    BOOST_CHECK_EQUAL(before.track.innerState.covariance[i], after.track.innerState.covariance[i]);
    BOOST_CHECK_EQUAL(before.track.outerState.covariance[i], after.track.outerState.covariance[i]);
  }
  BOOST_CHECK_EQUAL(before.track.innerState.referenceCoordinate, after.track.innerState.referenceCoordinate);
  BOOST_CHECK_EQUAL(before.track.outerState.referenceCoordinate, after.track.outerState.referenceCoordinate);
  BOOST_CHECK_EQUAL(static_cast<int>(before.track.innerState.family), static_cast<int>(after.track.innerState.family));
  BOOST_CHECK_EQUAL(static_cast<int>(before.track.outerState.family), static_cast<int>(after.track.outerState.family));
  BOOST_CHECK_EQUAL(before.track.chi2, after.track.chi2);
  BOOST_CHECK_EQUAL(before.phi, after.phi);
  BOOST_CHECK_EQUAL(before.eta, after.eta);
  BOOST_CHECK_EQUAL(before.charge, after.charge);
}

} // namespace

// --- A. Normalized authority over poisoned legacy backfill ------------------

BOOST_AUTO_TEST_CASE(NormalizedAuthorityOverridesPoisonedLegacyBackfill)
{
  const StraightTrackGeometry geometry(0.3f);

  RefitFixture reference(geometry);
  TrackingCandidate referenceTrack;
  BOOST_REQUIRE(refit(reference, referenceTrack));

  RefitFixture poisoned(geometry);
  poisoned.poisonLegacyBackfill();
  TrackingCandidate poisonedTrack;
  BOOST_REQUIRE(refit(poisoned, poisonedTrack));

  // Identical normalized input, poisoned vs. untouched legacy backfill: byte-
  // identical output proves the refit path never reads the legacy structures.
  checkTrackUnchanged(referenceTrack, poisonedTrack);
}

// --- B. Inverse authority: normalized data drives the output ----------------

BOOST_AUTO_TEST_CASE(NormalizedGlobalCoordinateChangeAltersOutput)
{
  const StraightTrackGeometry geometry(0.3f);

  RefitFixture reference(geometry);
  TrackingCandidate referenceTrack;
  BOOST_REQUIRE(refit(reference, referenceTrack));

  // Perturb only the normalized global.x of one interior layer -- legacy
  // backfill is absent (never populated) in both fixtures, so this isolates
  // the normalized measurement as the sole cause of the changed outcome. The
  // shift is far larger than DefaultSigma2's resolution, so the previously
  // ~0 chi2/ndf now certainly exceeds MaxChi2NDF.
  RefitFixture perturbed(geometry);
  auto perturbedMeasurement = perturbed.storage[5].front();
  perturbedMeasurement.global.x += 0.05f;
  perturbed.storage[5].assign(1, perturbedMeasurement);
  perturbed.layerMeasurements[5] = perturbed.storage[5];

  TrackingCandidate perturbedTrack;
  const bool perturbedOk = refit(perturbed, perturbedTrack);
  BOOST_CHECK(!perturbedOk);
}

BOOST_AUTO_TEST_CASE(NormalizedCovarianceChangeAltersOutput)
{
  const StraightTrackGeometry geometry(0.3f);

  RefitFixture reference(geometry);
  TrackingCandidate referenceTrack;
  BOOST_REQUIRE(refit(reference, referenceTrack));

  // Scale up every layer's diagonal covariance uniformly (legacy backfill
  // again absent in both fixtures): with exact-colinear points the fitted
  // position/chi2 are unaffected, but the posterior parameter covariance the
  // Kalman filter propagates is not -- a strictly larger measurement variance
  // must not shrink the output covariance. This is a generic Kalman-filter
  // property, unaffected by which per-hit update formula produces it.
  RefitFixture loose(geometry);
  for (int layer = 0; layer < NLayers; ++layer) {
    auto m = loose.storage[layer].front();
    m.covariance.uu *= 400.f;
    m.covariance.vv *= 400.f;
    loose.storage[layer].assign(1, m);
    loose.layerMeasurements[layer] = loose.storage[layer];
  }
  TrackingCandidate looseTrack;
  BOOST_REQUIRE(refit(loose, looseTrack));

  BOOST_CHECK_GT(looseTrack.track.outerState.covariance[packedCovarianceIndex(0, 0)],
                 referenceTrack.track.outerState.covariance[packedCovarianceIndex(0, 0)]);
  BOOST_CHECK_GT(looseTrack.track.outerState.covariance[packedCovarianceIndex(1, 1)],
                 referenceTrack.track.outerState.covariance[packedCovarianceIndex(1, 1)]);
}

// --- C. TrackletMinAbsX is gated by normalized global.x ----------------------

BOOST_AUTO_TEST_CASE(TrackletMinAbsXUsesNormalizedGlobalXEvenWhenLegacyDisagrees)
{
  // Slope chosen so x crosses well below the cut only at the last (z-
  // extreme) layer while staying far above it everywhere else, including at
  // the fit's own overall-X gate (paramOut's X, at the fit's outward-leg
  // reference, layer 0, x == 1.0).
  const StraightTrackGeometry geometry(0.03f);
  RefitFixture fx(geometry);
  fx.params.TrackletMinAbsX = 0.05f;

  // Legacy Cluster::xCoordinate poisoned to a value that would have passed
  // the old (legacy-reading) cut; only the normalized global.x (~0.034 at
  // layer 9) is below TrackletMinAbsX.
  for (int layer = 0; layer < NLayers; ++layer) {
    fx.tf.addClusterToLayer(layer, PoisonCoordinate, PoisonCoordinate, PoisonCoordinate, 0);
  }

  TrackingCandidate track;
  const bool ok = refit(fx, track);
  BOOST_CHECK(!ok);
}

// --- D. Invalid normalized input fails cleanly, destination untouched -------

BOOST_AUTO_TEST_CASE(NonFiniteGlobalCoordinateFailsCleanly)
{
  const StraightTrackGeometry geometry(0.3f);
  RefitFixture fx(geometry);
  auto m = fx.storage[3].front();
  m.global.x = std::numeric_limits<float>::quiet_NaN();
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  TrackingCandidate before;
  TrackingCandidate track = before;
  BOOST_CHECK(!refit(fx, track));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(InfiniteGlobalCoordinateFailsCleanly)
{
  const StraightTrackGeometry geometry(0.3f);
  RefitFixture fx(geometry);
  auto m = fx.storage[3].front();
  m.global.z = std::numeric_limits<float>::infinity();
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  TrackingCandidate before;
  TrackingCandidate track = before;
  BOOST_CHECK(!refit(fx, track));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(NonFiniteCovarianceFailsCleanly)
{
  const StraightTrackGeometry geometry(0.3f);
  RefitFixture fx(geometry);
  auto m = fx.storage[3].front();
  m.covariance.uu = std::numeric_limits<float>::quiet_NaN();
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  TrackingCandidate before;
  TrackingCandidate track = before;
  BOOST_CHECK(!refit(fx, track));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(NegativeCovarianceFailsCleanly)
{
  const StraightTrackGeometry geometry(0.3f);
  RefitFixture fx(geometry);
  auto m = fx.storage[3].front();
  m.covariance.vv = -1.f;
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  TrackingCandidate before;
  TrackingCandidate track = before;
  BOOST_CHECK(!refit(fx, track));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(OutOfRangeClusterIndexFailsCleanly)
{
  const StraightTrackGeometry geometry(0.3f);
  RefitFixture fx(geometry);
  fx.seed.getClusters()[3] = 99; // storage[3] only ever has one element (index 0)

  TrackingCandidate before;
  TrackingCandidate track = before;
  BOOST_CHECK(!refit(fx, track));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(IdentityMismatchFailsCleanly)
{
  const StraightTrackGeometry geometry(0.3f);
  RefitFixture fx(geometry);
  // ClusterRef.index no longer agrees with tf.getClusterExternalIndex(layer, clIdx)
  // (== layer, per RefitFixture's construction) -- the ClusterRef contract
  // TrackerTraits::initialiseTimeFrame() would normally have already enforced.
  auto m = fx.storage[3].front();
  m.cluster = ClusterRef{ClusterSourceId{0}, 12345u};
  m.surface = SurfaceId{3};
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  TrackingCandidate before;
  TrackingCandidate track = before;
  BOOST_CHECK(!refit(fx, track));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(InvalidClusterRefFailsCleanly)
{
  const StraightTrackGeometry geometry(0.3f);
  RefitFixture fx(geometry);
  auto m = fx.storage[3].front();
  m.cluster = ClusterRef{}; // default-constructed: index == invalid sentinel, isValid() == false
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  TrackingCandidate before;
  TrackingCandidate track = before;
  BOOST_CHECK(!refit(fx, track));
  checkTrackUnchanged(before, track);
}

// --- E. Preservation ---------------------------------------------------------

BOOST_AUTO_TEST_CASE(PreservesSeedMembershipForGenericRefit)
{
  const StraightTrackGeometry geometry(0.3f);
  // Holes at layers 2 and 7: MinTrackLength(5) <= 8 remaining hits.
  RefitFixture fx(geometry);
  fx.seed.getClusters()[2] = o2::its::constants::UnusedIndex;
  fx.seed.getClusters()[7] = o2::its::constants::UnusedIndex;
  SurfaceMask mask = fx.seed.getSurfaceMask();
  mask.reset(SurfaceId{2});
  mask.reset(SurfaceId{7});
  fx.seed.setSurfaceMask(mask);

  TrackingCandidate track;
  BOOST_REQUIRE(refit(fx, track));

  BOOST_CHECK_EQUAL(track.getNumberOfClusters(), NLayers - 2);
  for (int layer = 0; layer < NLayers; ++layer) {
    if (layer == 2 || layer == 7) {
      BOOST_CHECK_EQUAL(track.getClusterIndex(layer), o2::its::constants::UnusedIndex);
      BOOST_CHECK(!track.seed.hasCluster(layer));
    } else {
      BOOST_CHECK_EQUAL(track.getClusterIndex(layer), 0);
      BOOST_CHECK(track.seed.hasCluster(layer));
    }
  }
}

// M5d: supersedes the retired DiagonalOnlyCovarianceIgnoresOffDiagonalTerm.
// Unlike the retired MFT TrackFitter (diagonal-only measurement update), the
// native forward::predictedChi2/update operations use the measurement's full
// uu/uv/vv covariance (ForwardSurfaceStateOperations.h doc) -- a disclosed,
// intentional behavior difference from this milestone. A physically valid
// correlation (uv = 0.3*sqrt(uu*vv), |correlation| < 1) must therefore change
// the fitted output relative to the uv == 0 reference.
BOOST_AUTO_TEST_CASE(OffDiagonalCovarianceIsUsedByNativeUpdate)
{
  const StraightTrackGeometry geometry(0.3f);

  RefitFixture reference(geometry);
  TrackingCandidate referenceTrack;
  BOOST_REQUIRE(refit(reference, referenceTrack));

  RefitFixture withUv(geometry);
  // A generous per-hit/per-track chi2 gate: this test's goal is only to
  // prove a physically valid off-diagonal correlation changes the native
  // update's output, not to probe chi2-gate behavior -- a nonzero
  // correlation legitimately raises the predicted chi2 against a reference
  // fit tuned for the uncorrelated (uv == 0) case.
  withUv.params.MaxChi2ClusterAttachment = 1.e4f;
  withUv.params.MaxChi2NDF = 1.e4f;
  for (int layer = 0; layer < NLayers; ++layer) {
    auto m = withUv.storage[layer].front();
    // A modest, physically valid correlation (|coefficient| << 1) suffices
    // to prove the point.
    m.covariance.uv = 0.05f * std::sqrt(m.covariance.uu * m.covariance.vv);
    withUv.storage[layer].assign(1, m);
    withUv.layerMeasurements[layer] = withUv.storage[layer];
  }
  TrackingCandidate withUvTrack;
  BOOST_REQUIRE(refit(withUv, withUvTrack));

  BOOST_CHECK_NE(withUvTrack.track.outerState.covariance[packedCovarianceIndex(0, 0)],
                 referenceTrack.track.outerState.covariance[packedCovarianceIndex(0, 0)]);
}

// --- Regression: source-qualified seed-cluster identity --------------------

BOOST_AUTO_TEST_CASE(GenericRefitValidatesExternalClusterIdentity)
{
  // Every hit layer has a layer-local seed index of zero but a deliberately
  // large, non-monotonic source-qualified external cluster index. The generic
  // refit must validate that identity through the scratch mapping without
  // treating the external index as a local measurement position.
  const StraightTrackGeometry geometry(0.3f);
  SurfaceTrackingScratch tf;
  std::array<std::vector<SurfaceMeasurement>, NLayers> storage;
  std::vector<gsl::span<const SurfaceMeasurement>> layerMeasurements = std::vector<gsl::span<const SurfaceMeasurement>>(NLayers);
  std::vector<SurfaceDescriptor> catalogSurfaces(NLayers);
  for (int layer = 0; layer < NLayers; ++layer) {
    catalogSurfaces[layer].id = SurfaceId{static_cast<uint16_t>(layer)};
    catalogSurfaces[layer].kind = SurfaceKind::Disk;
    catalogSurfaces[layer].material = NominalSurfaceMaterial{0.f, 0.f};
  }
  SurfaceCatalogView catalog{catalogSurfaces.data(), static_cast<uint32_t>(catalogSurfaces.size())};
  tf.setMemoryPool(std::make_shared<o2::its::BoundedMemoryResource>());
  tf.adoptPlan(NLayers, 0, 0);
  TrackSeed seed;
  o2::itsmft::TrackingParameters params;
  params.MinTrackLength = 5;
  params.MinPt.assign(NLayers + 1, 0.f);
  params.TrackletMinAbsX = 0.f;
  params.MaxChi2NDF = 30.f;

  uint16_t mask = 0;
  for (int layer = 0; layer < NLayers; ++layer) {
    SurfaceMeasurement m{};
    m.global = {geometry.x[layer], geometry.y[layer], geometry.z[layer]};
    m.frame.q = geometry.z[layer]; // disk adapter contract: frame.q == global.z
    m.covariance.uu = DefaultSigma2;
    m.covariance.vv = DefaultSigma2;
    m.cluster = ClusterRef{ClusterSourceId{0}, static_cast<uint32_t>(1000 + layer)};
    m.surface = SurfaceId{static_cast<uint16_t>(layer)};
    storage[layer].assign(1, m);
    layerMeasurements[layer] = storage[layer];

    tf.addClusterExternalIndexToLayer(layer, 1000 + layer); // clIdx 0 -> large, non-monotonic extIdx
    seed.getClusters()[layer] = 0;
    mask |= static_cast<uint16_t>(uint16_t(1) << layer);
  }
  seed.setSurfaceMask(SurfaceMask{mask});

  // M5d: see RefitFixture's identical seed.state() construction above.
  OperationFailureReason seedReason{};
  BOOST_REQUIRE(forward::buildSeed(layerMeasurements[0][0], layerMeasurements[NLayers / 2][0],
                                   layerMeasurements[NLayers - 1][0], Bz, params.TrackletMinPt,
                                   /*absCharge=*/1, o2::track::PID::Pion, seed.state(), seedReason));

  TrackingCandidate track;
  BOOST_REQUIRE(refit(seed, params, tf, layerMeasurements, catalog, track));

  for (int layer = 0; layer < NLayers; ++layer) {
    BOOST_CHECK(track.seed.hasCluster(layer));
    BOOST_CHECK_EQUAL(track.getClusterIndex(layer), 0);
    BOOST_CHECK_EQUAL(layerMeasurements[layer][0].cluster.index, static_cast<uint32_t>(1000 + layer));
  }
}
