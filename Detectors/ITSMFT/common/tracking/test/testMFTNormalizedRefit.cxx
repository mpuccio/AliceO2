// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Stage-B MFT retained-refitter normalized-measurement authority slice
// (Architecture.md Sec 6.5/12, AgentCoordination.md MFT-adapter role): focused
// coverage for refitTrackFwd/refitTrackFwdImpl (MFTFwdTrackHelpers.cxx),
// proving every physical hit coordinate/covariance the frozen
// o2::mft::TrackFitter/TrackLTF consumes comes from the caller-supplied
// LayerMeasurementSpans<NLayers> (TrackerTraits::mLayerMeasurements), never
// from TimeFrame's legacy Cluster/TrackingFrameInfo backfill. Each test calls
// refitTrackFwd directly with a hand-built TimeFrame/seed/measurement-span
// fixture -- it does not exercise TrackerTraits::initialiseTimeFrame() or the
// full CA traversal, which already have their own focused coverage
// (testComputeLayerCellsOrchestration.cxx et al.) for the
// NormalizedMeasurementMismatch load-time contract this file's fixtures
// assume as given.

#define BOOST_TEST_MODULE ITSMFT MFTNormalizedRefit
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Constants.h"
#include "MFTTracking/Constants.h"

using namespace o2::itsmft::tracking;

namespace
{

constexpr int NLayers = o2::mft::constants::mft::LayersNumber;
// Field-off dispatch (refitTrackFwd: std::abs(bz) < 1e-6f): exercises the
// TrackLTFL linear-model fit deterministically, independent of
// o2::mft::MFTTrackingParam::Instance()'s global forceZeroField default.
constexpr float Bz = 0.f;
constexpr float DefaultSigma2 = 2.5e-7f; // (~0.5 micron)^2, MFT-scale resolution
constexpr float PoisonCoordinate = 9999.f;

// A straight track through every MFT disk. Points are exact (analytic) so the
// retained Kalman fit converges with ~0 chi2 regardless of which reasonable
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

// Owns everything refitTrackFwd needs for one call: a TimeFrame (external-
// index/size bookkeeping only -- never a physical-coordinate source once
// populated), the normalized measurement spans under test, and the seed.
// Every test builds its own fixture so mutating one layer's normalized
// measurement or legacy backfill in one test can never leak into another.
struct RefitFixture {
  TimeFrame<NLayers> tf;
  std::array<std::vector<SurfaceMeasurement>, NLayers> storage;
  LayerMeasurementSpans<NLayers> layerMeasurements{};
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
  }

  void setMeasurement(int layer, float x, float y, float z, float uu, float vv, uint32_t clusterIndex, float uv = 0.f)
  {
    SurfaceMeasurement m{};
    m.global = {x, y, z};
    m.covariance.uu = uu;
    m.covariance.vv = vv;
    m.covariance.uv = uv;
    m.cluster = ClusterRef{ClusterSourceId{0}, clusterIndex};
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
                                              PoisonCoordinate, PoisonCoordinate, PoisonCoordinate, PoisonCoordinate, PoisonCoordinate,
                                              {PoisonCoordinate, PoisonCoordinate}, {1.f, 0.f, 1.f}});
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
  const StraightTrackGeometry geometry(0.01f);

  RefitFixture reference(geometry);
  MFTCATrack referenceTrack;
  BOOST_REQUIRE(refitTrackFwd(reference.seed, referenceTrack, reference.tf, reference.params, Bz, reference.layerMeasurements));

  RefitFixture poisoned(geometry);
  poisoned.poisonLegacyBackfill();
  MFTCATrack poisonedTrack;
  BOOST_REQUIRE(refitTrackFwd(poisoned.seed, poisonedTrack, poisoned.tf, poisoned.params, Bz, poisoned.layerMeasurements));

  // Identical normalized input, poisoned vs. untouched legacy backfill: byte-
  // identical output proves the refit path never reads the legacy structures.
  checkTrackUnchanged(referenceTrack, poisonedTrack);
}

// --- B. Inverse authority: normalized data drives the output ----------------

BOOST_AUTO_TEST_CASE(NormalizedGlobalCoordinateChangeAltersOutput)
{
  const StraightTrackGeometry geometry(0.01f);

  RefitFixture reference(geometry);
  MFTCATrack referenceTrack;
  BOOST_REQUIRE(refitTrackFwd(reference.seed, referenceTrack, reference.tf, reference.params, Bz, reference.layerMeasurements));

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
  const bool perturbedOk = refitTrackFwd(perturbed.seed, perturbedTrack, perturbed.tf, perturbed.params, Bz, perturbed.layerMeasurements);
  BOOST_CHECK(!perturbedOk);
}

BOOST_AUTO_TEST_CASE(NormalizedCovarianceChangeAltersOutput)
{
  const StraightTrackGeometry geometry(0.01f);

  RefitFixture reference(geometry);
  MFTCATrack referenceTrack;
  BOOST_REQUIRE(refitTrackFwd(reference.seed, referenceTrack, reference.tf, reference.params, Bz, reference.layerMeasurements));

  // Scale up every layer's diagonal covariance uniformly (legacy backfill
  // again absent in both fixtures): with exact-colinear points the fitted
  // position/chi2 are unaffected, but the posterior parameter covariance the
  // Kalman filter propagates is not -- a strictly larger measurement variance
  // must not shrink the output covariance.
  RefitFixture loose(geometry);
  for (int layer = 0; layer < NLayers; ++layer) {
    auto m = loose.storage[layer].front();
    m.covariance.uu *= 400.f;
    m.covariance.vv *= 400.f;
    loose.storage[layer].assign(1, m);
    loose.layerMeasurements[layer] = loose.storage[layer];
  }
  MFTCATrack looseTrack;
  BOOST_REQUIRE(refitTrackFwd(loose.seed, looseTrack, loose.tf, loose.params, Bz, loose.layerMeasurements));

  BOOST_CHECK_GT(looseTrack.getTrack().getSigma2X(), referenceTrack.getTrack().getSigma2X());
  BOOST_CHECK_GT(looseTrack.getTrack().getSigma2Y(), referenceTrack.getTrack().getSigma2Y());
}

// --- C. TrackletMinAbsX is gated by normalized global.x ----------------------

BOOST_AUTO_TEST_CASE(TrackletMinAbsXUsesNormalizedGlobalXEvenWhenLegacyDisagrees)
{
  // Slope chosen so x crosses well below the cut only at the last (z-
  // extreme) layer while staying far above it everywhere else, including at
  // the fit's own overall-X gate (ltf.getX(), evaluated at the innermost
  // point, layer 0, x == 1.0).
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
  const bool ok = refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements);
  BOOST_CHECK(!ok);
}

// --- D. Invalid normalized input fails cleanly, destination untouched -------

BOOST_AUTO_TEST_CASE(NonFiniteGlobalCoordinateFailsCleanly)
{
  const StraightTrackGeometry geometry(0.01f);
  RefitFixture fx(geometry);
  auto m = fx.storage[3].front();
  m.global.x = std::numeric_limits<float>::quiet_NaN();
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(InfiniteGlobalCoordinateFailsCleanly)
{
  const StraightTrackGeometry geometry(0.01f);
  RefitFixture fx(geometry);
  auto m = fx.storage[3].front();
  m.global.z = std::numeric_limits<float>::infinity();
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(NonFiniteCovarianceFailsCleanly)
{
  const StraightTrackGeometry geometry(0.01f);
  RefitFixture fx(geometry);
  auto m = fx.storage[3].front();
  m.covariance.uu = std::numeric_limits<float>::quiet_NaN();
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(NegativeCovarianceFailsCleanly)
{
  const StraightTrackGeometry geometry(0.01f);
  RefitFixture fx(geometry);
  auto m = fx.storage[3].front();
  m.covariance.vv = -1.f;
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(OutOfRangeClusterIndexFailsCleanly)
{
  const StraightTrackGeometry geometry(0.01f);
  RefitFixture fx(geometry);
  fx.seed.getClusters()[3] = 99; // storage[3] only ever has one element (index 0)

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(IdentityMismatchFailsCleanly)
{
  const StraightTrackGeometry geometry(0.01f);
  RefitFixture fx(geometry);
  // ClusterRef.index no longer agrees with tf.getClusterExternalIndex(layer, clIdx)
  // (== layer, per RefitFixture's construction) -- the ClusterRef contract
  // TrackerTraits::initialiseTimeFrame() would normally have already enforced.
  auto m = fx.storage[3].front();
  m.cluster = ClusterRef{ClusterSourceId{0}, 12345u};
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements));
  checkTrackUnchanged(before, track);
}

BOOST_AUTO_TEST_CASE(InvalidClusterRefFailsCleanly)
{
  const StraightTrackGeometry geometry(0.01f);
  RefitFixture fx(geometry);
  auto m = fx.storage[3].front();
  m.cluster = ClusterRef{}; // default-constructed: index == invalid sentinel, isValid() == false
  fx.storage[3].assign(1, m);
  fx.layerMeasurements[3] = fx.storage[3];

  MFTCATrack before;
  MFTCATrack track = before;
  BOOST_CHECK(!refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements));
  checkTrackUnchanged(before, track);
}

// --- E. Preservation ---------------------------------------------------------

BOOST_AUTO_TEST_CASE(PreservesHitPatternClusterIndicesAndSizes)
{
  const StraightTrackGeometry geometry(0.01f);
  // Holes at layers 2 and 7: MinTrackLength(5) <= 8 remaining hits.
  RefitFixture fx(geometry);
  fx.seed.getClusters()[2] = o2::its::constants::UnusedIndex;
  fx.seed.getClusters()[7] = o2::its::constants::UnusedIndex;
  LayerMask mask = fx.seed.getHitLayerMask();
  mask &= ~LayerMask{static_cast<uint16_t>(uint16_t(1) << 2)};
  mask &= ~LayerMask{static_cast<uint16_t>(uint16_t(1) << 7)};
  fx.seed.setHitLayerMask(mask);

  MFTCATrack track;
  BOOST_REQUIRE(refitTrackFwd(fx.seed, track, fx.tf, fx.params, Bz, fx.layerMeasurements));

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

BOOST_AUTO_TEST_CASE(DiagonalOnlyCovarianceIgnoresOffDiagonalTerm)
{
  const StraightTrackGeometry geometry(0.01f);

  RefitFixture reference(geometry);
  MFTCATrack referenceTrack;
  BOOST_REQUIRE(refitTrackFwd(reference.seed, referenceTrack, reference.tf, reference.params, Bz, reference.layerMeasurements));

  // The retained fitter is diagonal-only: measurement.covariance.uv must be
  // read by nothing on this path. A wild, physically nonsensical uv (this
  // large a correlation is invalid for any real 2x2 covariance) must not
  // change the output at all versus uv == 0.
  RefitFixture withUv(geometry);
  for (int layer = 0; layer < NLayers; ++layer) {
    auto m = withUv.storage[layer].front();
    m.covariance.uv = 1.e6f;
    withUv.storage[layer].assign(1, m);
    withUv.layerMeasurements[layer] = withUv.storage[layer];
  }
  MFTCATrack withUvTrack;
  BOOST_REQUIRE(refitTrackFwd(withUv.seed, withUvTrack, withUv.tf, withUv.params, Bz, withUv.layerMeasurements));

  checkTrackUnchanged(referenceTrack, withUvTrack);
}

// --- Regression: local-vs-external cluster-size index domain ---------------

BOOST_AUTO_TEST_CASE(ClusterSizeIsReadFromItsOwnLayerNotFromLayerZeroByExternalIndex)
{
  // Regression test for the bug fixed in MFTFwdTrackHelpers.cxx:
  // refitTrackFwdImpl must read tf.getClusterSize(layer, clIdx) -- clIdx
  // being THIS layer's own layer-local cluster identity -- and must never
  // read tf.getClusterSize(0, extIdx) using the external/global identity:
  // mClusterSize is a per-layer vector (see TimeFrame::loadNormalizedSource()
  // and TimeFrame::getClusterSize()), not a flat array addressable by an
  // external id.
  //
  // Every hit layer gets a real, distinctive size (10 + layer) at its own
  // layer-local clIdx == 0, and a deliberately large, non-monotonic external
  // cluster index (1000 + layer) -- proving the fix does not depend on
  // external indices being small or increasing with layer. Layer 0's own
  // cluster-size vector is then independently given poisoned entries at
  // exactly the slots (1000 + layer) that the old buggy code would have
  // read via tf.getClusterSize(0, 1000 + layer): under the fix that value
  // must never leak into any other layer's published size.
  const StraightTrackGeometry geometry(0.01f);
  TimeFrame<NLayers> tf;
  std::array<std::vector<SurfaceMeasurement>, NLayers> storage;
  LayerMeasurementSpans<NLayers> layerMeasurements{};
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
    m.covariance.uu = DefaultSigma2;
    m.covariance.vv = DefaultSigma2;
    m.cluster = ClusterRef{ClusterSourceId{0}, static_cast<uint32_t>(1000 + layer)};
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
  BOOST_REQUIRE(refitTrackFwd(seed, track, tf, params, Bz, layerMeasurements));

  for (int layer = 0; layer < NLayers; ++layer) {
    BOOST_CHECK(track.hasHitOnLayer(layer));
    BOOST_CHECK_EQUAL(track.getClusterSize(layer), 10 + layer);
    BOOST_CHECK_NE(track.getClusterSize(layer), 250);
  }
}
