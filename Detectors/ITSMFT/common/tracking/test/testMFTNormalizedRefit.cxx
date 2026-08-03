// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Stage-B MFT normalized-measurement authority slice (Architecture.md Sec
// 6.5/12, AgentCoordination.md MFT-adapter role), updated for M5d: focused
// coverage for refitTrackFwd/refitTrackFwdImpl (MFTFwdTrackHelpers.cxx),
// proving every physical hit coordinate/covariance the shared native driver
// (fitTrackSeedLegs, NativeRefitDriver.h) consumes comes from the
// caller-supplied LayerMeasurementSpans<NLayers> (TrackerTraits::
// mLayerMeasurements) and SurfaceCatalogView, never from LegacyTrackerScratch's
// legacy Cluster/TrackingFrameInfo backfill. Each test calls refitTrackFwd
// directly with a hand-built LegacyTrackerScratch/seed/measurement-span/
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
#include "ITSMFTTracking/LegacyTrackerScratch.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/SurfaceCatalogView.h"
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

// Owns everything refitTrackFwd needs for one call: a LegacyTrackerScratch
// (external-index/size bookkeeping only -- never a physical-coordinate
// source once populated), the normalized measurement spans under test, a
// matching zero-material Disk surface catalog (one SurfaceDescriptor per
// layer, id == layer, resolved by every measurement's own .surface field),
// and the seed. Every test builds its own fixture so mutating one layer's
// normalized measurement or legacy backfill in one test can never leak into
// another.
struct RefitFixture {
  LegacyTrackerScratch<NLayers> tf;
  std::array<std::vector<SurfaceMeasurement>, NLayers> storage;
  LayerMeasurementSpans<NLayers> layerMeasurements{};
  std::vector<SurfaceDescriptor> catalogSurfaces;
  SurfaceCatalogView catalog{};
  TrackSeedN<NLayers> seed;
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

    uint16_t mask = 0;
    for (int layer = 0; layer < hits; ++layer) {
      setMeasurement(layer, geometry.x[layer], geometry.y[layer], geometry.z[layer],
                     DefaultSigma2, DefaultSigma2, static_cast<uint32_t>(layer));
      tf.addClusterExternalIndexToLayer(layer, layer); // extIdx == layer, clIdx == 0
      // Per-layer mClusterSize, indexed by the layer-local clIdx (0 here) --
      // mirrors loadNormalizedSource()'s real per-layer layout. A single
      // flat mClusterSize[0] sized NLayers and indexed by extIdx (this
      // fixture's old setup) only matched the pre-fix bug's own
      // (incorrect) convention, not TimeFrame::getClusterSize()'s actual
      // per-layer contract.
      bounded_vector<uint8_t> layerSizes;
      layerSizes.resize(1, uint8_t{1});
      tf.setClusterSize(layer, layerSizes);
      seed.getClusters()[layer] = 0;
      mask |= static_cast<uint16_t>(uint16_t(1) << layer);
    }
    seed.setHitLayerMask(LayerMask{mask});

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

void checkTrackUnchanged(const MFTCATrack& before, const MFTCATrack& after)
{
  BOOST_CHECK_EQUAL(before.getPattern(), after.getPattern());
  BOOST_CHECK_EQUAL(before.getSeedPattern(), after.getSeedPattern());
  BOOST_CHECK_EQUAL(before.hasSharedClusters(), after.hasSharedClusters());
  BOOST_CHECK_EQUAL(before.getTrack().getX(), after.getTrack().getX());
  BOOST_CHECK_EQUAL(before.getTrack().getY(), after.getTrack().getY());
  BOOST_CHECK_EQUAL(before.getTrack().getZ(), after.getTrack().getZ());
  BOOST_CHECK_EQUAL(before.getTrack().getTanl(), after.getTrack().getTanl());
  BOOST_CHECK_EQUAL(before.getTrack().getInvQPt(), after.getTrack().getInvQPt());
  BOOST_CHECK_EQUAL(before.getTrack().getTrackChi2(), after.getTrack().getTrackChi2());
  BOOST_CHECK_EQUAL(before.getTrack().getNumberOfPoints(), after.getTrack().getNumberOfPoints());
  for (int layer = 0; layer < MFTCATrack::MaxClusters; ++layer) {
    BOOST_CHECK_EQUAL(before.getClusterIndex(layer), after.getClusterIndex(layer));
    BOOST_CHECK_EQUAL(before.getClusterSize(layer), after.getClusterSize(layer));
  }
}

} // namespace

// --- A. Normalized authority over poisoned legacy backfill ------------------

BOOST_AUTO_TEST_CASE(NormalizedAuthorityOverridesPoisonedLegacyBackfill)
{
  const StraightTrackGeometry geometry(0.3f);

  RefitFixture reference(geometry);
  MFTCATrack referenceTrack;
  BOOST_REQUIRE(refitTrackFwd(reference.seed, referenceTrack, reference.tf, reference.params, Bz, reference.layerMeasurements, reference.catalog, ClusterSourceId{0}));

  RefitFixture poisoned(geometry);
  poisoned.poisonLegacyBackfill();
  MFTCATrack poisonedTrack;
  BOOST_REQUIRE(refitTrackFwd(poisoned.seed, poisonedTrack, poisoned.tf, poisoned.params, Bz, poisoned.layerMeasurements, poisoned.catalog, ClusterSourceId{0}));

  // Identical normalized input, poisoned vs. untouched legacy backfill: byte-
  // identical output proves the refit path never reads the legacy structures.
  checkTrackUnchanged(referenceTrack, poisonedTrack);
}

// --- B. Inverse authority: normalized data drives the output ----------------

BOOST_AUTO_TEST_CASE(NormalizedGlobalCoordinateChangeAltersOutput)
{
  const StraightTrackGeometry geometry(0.3f);

  RefitFixture reference(geometry);
  MFTCATrack referenceTrack;
  BOOST_REQUIRE(refitTrackFwd(reference.seed, referenceTrack, reference.tf, reference.params, Bz, reference.layerMeasurements, reference.catalog, ClusterSourceId{0}));

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

  MFTCATrack perturbedTrack;
  const bool perturbedOk = refitTrackFwd(perturbed.seed, perturbedTrack, perturbed.tf, perturbed.params, Bz, perturbed.layerMeasurements, perturbed.catalog, ClusterSourceId{0});
  BOOST_CHECK(!perturbedOk);
}

BOOST_AUTO_TEST_CASE(NormalizedCovarianceChangeAltersOutput)
{
  const StraightTrackGeometry geometry(0.3f);

  RefitFixture reference(geometry);
  MFTCATrack referenceTrack;
  BOOST_REQUIRE(refitTrackFwd(reference.seed, referenceTrack, reference.tf, reference.params, Bz, reference.layerMeasurements, reference.catalog, ClusterSourceId{0}));

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
  MFTCATrack looseTrack;
  BOOST_REQUIRE(refitTrackFwd(loose.seed, looseTrack, loose.tf, loose.params, Bz, loose.layerMeasurements, loose.catalog, ClusterSourceId{0}));

  BOOST_CHECK_GT(looseTrack.getTrack().getSigma2X(), referenceTrack.getTrack().getSigma2X());
  BOOST_CHECK_GT(looseTrack.getTrack().getSigma2Y(), referenceTrack.getTrack().getSigma2Y());
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

  MFTCATrack track;
  const bool ok = refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements, fx.catalog, ClusterSourceId{0});
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

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements, fx.catalog, ClusterSourceId{0}));
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

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements, fx.catalog, ClusterSourceId{0}));
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

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements, fx.catalog, ClusterSourceId{0}));
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

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements, fx.catalog, ClusterSourceId{0}));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(OutOfRangeClusterIndexFailsCleanly)
{
  const StraightTrackGeometry geometry(0.3f);
  RefitFixture fx(geometry);
  fx.seed.getClusters()[3] = 99; // storage[3] only ever has one element (index 0)

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements, fx.catalog, ClusterSourceId{0}));
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

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements, fx.catalog, ClusterSourceId{0}));
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

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements, fx.catalog, ClusterSourceId{0}));
  checkTrackUnchanged(before, track);
}

// --- E. Preservation ---------------------------------------------------------

BOOST_AUTO_TEST_CASE(PreservesHitPatternClusterIndicesAndSizes)
{
  const StraightTrackGeometry geometry(0.3f);
  // Holes at layers 2 and 7: MinTrackLength(5) <= 8 remaining hits.
  RefitFixture fx(geometry);
  fx.seed.getClusters()[2] = o2::its::constants::UnusedIndex;
  fx.seed.getClusters()[7] = o2::its::constants::UnusedIndex;
  LayerMask mask = fx.seed.getHitLayerMask();
  mask &= ~LayerMask{static_cast<uint16_t>(uint16_t(1) << 2)};
  mask &= ~LayerMask{static_cast<uint16_t>(uint16_t(1) << 7)};
  fx.seed.setHitLayerMask(mask);

  MFTCATrack track;
  BOOST_REQUIRE(refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements, fx.catalog, ClusterSourceId{0}));

  BOOST_CHECK_EQUAL(track.getTrack().getNumberOfPoints(), NLayers - 2);
  for (int layer = 0; layer < NLayers; ++layer) {
    if (layer == 2 || layer == 7) {
      BOOST_CHECK_EQUAL(track.getClusterIndex(layer), o2::its::constants::UnusedIndex);
      BOOST_CHECK(!track.hasHitOnLayer(layer));
    } else {
      BOOST_CHECK_EQUAL(track.getClusterIndex(layer), 0);
      BOOST_CHECK_EQUAL(track.getClusterSize(layer), 1);
      BOOST_CHECK(track.hasHitOnLayer(layer));
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
  MFTCATrack referenceTrack;
  BOOST_REQUIRE(refitTrackFwd(reference.seed, referenceTrack, reference.tf, reference.params, Bz, reference.layerMeasurements, reference.catalog, ClusterSourceId{0}));

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
  MFTCATrack withUvTrack;
  BOOST_REQUIRE(refitTrackFwd(withUv.seed, withUvTrack, withUv.tf, withUv.params, Bz, withUv.layerMeasurements, withUv.catalog, ClusterSourceId{0}));

  BOOST_CHECK_NE(withUvTrack.getTrack().getSigma2X(), referenceTrack.getTrack().getSigma2X());
}

// --- Regression: local-vs-external cluster-size index domain ---------------

BOOST_AUTO_TEST_CASE(ClusterSizeIsReadFromItsOwnLayerNotFromLayerZeroByExternalIndex)
{
  // Regression test for the bug fixed in MFTFwdTrackHelpers.cxx:
  // refitTrackFwdImpl must read tf.getClusterSize(layer, clIdx) -- clIdx
  // being THIS layer's own layer-local cluster identity -- and must never
  // read tf.getClusterSize(0, extIdx) using the external/global identity:
  // mClusterSize is a per-layer vector (see
  // LegacyTrackerScratch::loadNormalizedSource() and
  // LegacyTrackerScratch::getClusterSize()), not a flat array addressable by
  // an external id.
  //
  // Every hit layer gets a real, distinctive size (10 + layer) at its own
  // layer-local clIdx == 0, and a deliberately large, non-monotonic external
  // cluster index (1000 + layer) -- proving the fix does not depend on
  // external indices being small or increasing with layer. Layer 0's own
  // cluster-size vector is then independently given poisoned entries at
  // exactly the slots (1000 + layer) that the old buggy code would have
  // read via tf.getClusterSize(0, 1000 + layer): under the fix that value
  // must never leak into any other layer's published size.
  const StraightTrackGeometry geometry(0.3f);
  LegacyTrackerScratch<NLayers> tf;
  std::array<std::vector<SurfaceMeasurement>, NLayers> storage;
  LayerMeasurementSpans<NLayers> layerMeasurements{};
  std::vector<SurfaceDescriptor> catalogSurfaces(NLayers);
  for (int layer = 0; layer < NLayers; ++layer) {
    catalogSurfaces[layer].id = SurfaceId{static_cast<uint16_t>(layer)};
    catalogSurfaces[layer].kind = SurfaceKind::Disk;
    catalogSurfaces[layer].material = NominalSurfaceMaterial{0.f, 0.f};
  }
  SurfaceCatalogView catalog{catalogSurfaces.data(), static_cast<uint32_t>(catalogSurfaces.size())};
  TrackSeedN<NLayers> seed;
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
    bounded_vector<uint8_t> sizes;
    sizes.resize(1, static_cast<uint8_t>(10 + layer));
    tf.setClusterSize(layer, sizes);

    seed.getClusters()[layer] = 0;
    mask |= static_cast<uint16_t>(uint16_t(1) << layer);
  }
  seed.setHitLayerMask(LayerMask{mask});

  // M5d: see RefitFixture's identical seed.state() construction above.
  OperationFailureReason seedReason{};
  BOOST_REQUIRE(forward::buildSeed(layerMeasurements[0][0], layerMeasurements[NLayers / 2][0],
                                   layerMeasurements[NLayers - 1][0], Bz, params.TrackletMinPt,
                                   /*absCharge=*/1, o2::track::PID::Pion, seed.state(), seedReason));

  // Poison layer 0's own vector at every slot the old (buggy) code would
  // have queried for another layer.
  bounded_vector<uint8_t> poisonedLayer0;
  poisonedLayer0.resize(1000 + NLayers, uint8_t{0});
  poisonedLayer0[0] = uint8_t{10}; // layer 0's real value, at its own clIdx == 0
  for (int layer = 1; layer < NLayers; ++layer) {
    poisonedLayer0[1000 + layer] = uint8_t{250}; // distinct poison value, never a real size here
  }
  tf.setClusterSize(0, poisonedLayer0);

  MFTCATrack track;
  BOOST_REQUIRE(refitTrackFwd(seed, track, tf, params, Bz, layerMeasurements, catalog, ClusterSourceId{0}));

  for (int layer = 0; layer < NLayers; ++layer) {
    BOOST_CHECK(track.hasHitOnLayer(layer));
    BOOST_CHECK_EQUAL(track.getClusterSize(layer), 10 + layer);
    BOOST_CHECK_NE(track.getClusterSize(layer), 250);
  }
}
