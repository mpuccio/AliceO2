// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 2 correction, updated for Gate 4 B3.1's TimeFrame /
// LegacyTrackerScratch<NLayers> split: TimeFrame lifecycle and
// transactional legacy backfill.
//
// A. Wipe lifecycle: TimeFrame::wipe() must unconditionally clear the
//    normalized owner (mNormalizedFrame) associated by
//    LegacyTrackerScratch<NLayers>::loadNormalizedSource() -- every
//    normalized accessor obtained *after* wipe() must report empty/zero
//    content. This test never dereferences a view obtained before wipe():
//    every post-wipe check re-obtains its accessor. (Gate 4 B3.1: neither
//    TimeFrame nor LegacyTrackerScratch<NLayers> stores mDetId at all any
//    more -- callers pass the detector explicitly to every call that needs
//    it -- so there is nothing detector-identity-shaped left for wipe() to
//    preserve or clear.)
//
// B. Strong exception transactionality: the owner-level load operation
//    LegacyTrackerScratch<NLayers>::loadNormalizedSource(TimeFrame&, ...)
//    stages both the shared TimeFrame's normalized update and its own
//    legacy backfill (unsorted clusters, TrackingFrameInfo, external
//    indices, cluster sizes, ROF boundaries, label pointers) before
//    touching either live owner. If a BoundedMemoryResource allocation
//    fails while building the staged legacy backfill -- after loadSources()
//    has already decoded a valid normalized frame into local, not-yet-
//    committed storage -- the thrown MemoryLimitExceeded must leave both
//    the TimeFrame's normalized owner and every legacy layer's
//    compatibility structures on the scratch exactly at their pre-call
//    baseline: a failed stage commits to neither owner.

#define BOOST_TEST_MODULE ITSMFT TimeFrame lifecycle
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <memory>
#include <optional>
#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/LegacyTrackerScratch.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

// Deterministic, geometry-free stand-in for GeometryClusterDecoder<DetId>
// (same construction as testTimeFrameNormalizedSource.cxx / testMultiSourceLoading.cxx):
// sensorID is used directly as the detector-local layer, global/frame
// coordinates are pure functions of (sensorID, row, col), and pattern
// consumption goes through the real production helper so cursor bookkeeping
// is exercised identically to GeometryClusterDecoder.
class LegacyLikeDecoder final : public ClusterDecoder
{
 public:
  explicit LegacyLikeDecoder(o2::detectors::DetID::ID detector) : mDetector(detector) {}

  o2::itsmft::ioutils::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const TopologyDictionary* dict,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool applySysErrors) const override
  {
    const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dict);
    if (!clusterData.ok()) {
      o2::itsmft::ioutils::SurfaceMeasurementDecodeResult result;
      result.error = clusterData.error;
      return result;
    }

    o2::itsmft::ioutils::SurfaceMeasurementDecodeResult result;
    const int sensorID = cluster.getSensorID();
    const int layer = sensorID;
    result.layer = layer;
    if (layer < 0 || static_cast<size_t>(layer) >= layerToSurface.size()) {
      return result;
    }
    result.layerMapped = true;
    result.kind = SurfaceKind::Cylinder;

    DecodedCluster decoded{};
    decoded.global = {static_cast<float>(sensorID) * 10.f, static_cast<float>(cluster.getRow()), static_cast<float>(cluster.getCol())};
    decoded.cylinderFrame = {static_cast<float>(sensorID) + 100.f, static_cast<float>(cluster.getRow()) + 1.f, static_cast<float>(cluster.getCol()) + 2.f, 0.01f * sensorID};
    decoded.rowColumnCovariance = {clusterData.sig2Row, 0.f, clusterData.sig2Col};
    decoded.shape = clusterData.shape;
    decoded.sensor = static_cast<uint32_t>(sensorID);
    decoded.layer = layer;

    const auto surface = layerToSurface[layer];
    const DetectorSensorId sensor{static_cast<uint32_t>(mDetector), decoded.sensor};
    const ClusterRef clusterRef{source, externalIndex};
    result.measurement = makeCylinderSurfaceMeasurement(decoded, sensor, surface, clusterRef, sourceROF);
    // Counts only clusters this decoder actually turned into a measurement
    // (the early-return failure paths above never reach here), so a test can
    // prove every cluster of a given input was successfully decoded by
    // checking how much this counter advanced across that call.
    ++decodeCount;
    return result;
  }

  mutable int decodeCount{0};

 private:
  o2::detectors::DetID::ID mDetector;
};

const TopologyDictionary& dict()
{
  static const TopologyDictionary d;
  return d;
}

constexpr std::array<unsigned char, 3> onePixelPattern{1, 1, 0x80};   // 1x1, 1 pixel
constexpr std::array<unsigned char, 3> threePixelPattern{1, 3, 0xE0}; // 1x3, 3 pixels

std::vector<unsigned char> concatPatterns(std::initializer_list<gsl::span<const unsigned char>> parts)
{
  std::vector<unsigned char> bytes;
  for (const auto& p : parts) {
    bytes.insert(bytes.end(), p.begin(), p.end());
  }
  return bytes;
}

std::vector<SurfaceDescriptor> makeITSTestCatalog()
{
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.reserve(ITSNLayers);
  for (uint16_t i = 0; i < ITSNLayers; ++i) {
    surfaces.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  }
  return surfaces;
}

std::vector<SurfaceId> identitySurfaces(uint16_t nLayers)
{
  std::vector<SurfaceId> mapping;
  mapping.reserve(nLayers);
  for (uint16_t i = 0; i < nLayers; ++i) {
    mapping.push_back(SurfaceId{i});
  }
  return mapping;
}

GlobalPoint3F expectedGlobal(int sensorID, int row, int col)
{
  return {static_cast<float>(sensorID) * 10.f, static_cast<float>(row), static_cast<float>(col)};
}

struct Fixture {
  std::vector<CompClusterExt> clusters;
  std::vector<unsigned char> patterns;
  std::vector<ROFRecord> rofs;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
};

// 4 clusters on layers {0,1,0,2}, partitioned into 3 ROFs: ROF0={c0,c1},
// ROF1={c2}, ROF2={c3}. Identical shape to testTimeFrameNormalizedSource.cxx's
// fixture, so parity with that accepted test coverage is preserved.
Fixture makeFixture()
{
  Fixture f;
  f.clusters = {
    CompClusterExt{10, 20, CompCluster::InvalidPatternID, 0}, // sensor 0 -> layer 0
    CompClusterExt{11, 21, CompCluster::InvalidPatternID, 1}, // sensor 1 -> layer 1
    CompClusterExt{12, 22, CompCluster::InvalidPatternID, 0}, // sensor 0 -> layer 0
    CompClusterExt{13, 23, CompCluster::InvalidPatternID, 2}, // sensor 2 -> layer 2
  };
  f.patterns = concatPatterns({onePixelPattern, threePixelPattern, onePixelPattern, threePixelPattern});
  f.rofs = {
    ROFRecord{{100, 5}, 0, 0, 2},
    ROFRecord{{140, 5}, 1, 2, 1},
    ROFRecord{{1000, 6}, 2, 3, 1}};
  for (uint32_t i = 0; i < f.clusters.size(); ++i) {
    f.labels.addElement(i, o2::MCCompLabel{static_cast<int>(i) + 1, 0, 0});
  }
  return f;
}

// A second, distinct, independently valid fixture: different sensors/layers
// (3,4,3,5,3 instead of 0,1,0,2), different rows/columns, a different
// pattern arrangement, a different ROF partition (3 ROFs over 5 clusters
// instead of 4), and its own separate MCTruthContainer with different label
// values. Used as the *replacement* load in the strong-exception-safety
// test, so that if any partial commit ever leaked through, it would be
// observable as foreign data (wrong layer, wrong coordinates, wrong label)
// rather than being masked by coincidentally reloading the same values.
Fixture makeReplacementFixture()
{
  Fixture f;
  f.clusters = {
    CompClusterExt{50, 60, CompCluster::InvalidPatternID, 3}, // sensor 3 -> layer 3
    CompClusterExt{51, 61, CompCluster::InvalidPatternID, 4}, // sensor 4 -> layer 4
    CompClusterExt{52, 62, CompCluster::InvalidPatternID, 3}, // sensor 3 -> layer 3
    CompClusterExt{53, 63, CompCluster::InvalidPatternID, 5}, // sensor 5 -> layer 5
    CompClusterExt{54, 64, CompCluster::InvalidPatternID, 3}, // sensor 3 -> layer 3
  };
  f.patterns = concatPatterns({threePixelPattern, threePixelPattern, onePixelPattern, threePixelPattern, onePixelPattern});
  f.rofs = {
    ROFRecord{{500, 1}, 0, 0, 3},
    ROFRecord{{540, 1}, 1, 3, 1},
    ROFRecord{{2000, 2}, 2, 4, 1}};
  for (uint32_t i = 0; i < f.clusters.size(); ++i) {
    f.labels.addElement(i, o2::MCCompLabel{static_cast<int>(i) + 101, 1, 1});
  }
  return f;
}

struct Expected {
  uint32_t externalIndex;
  int layer;
  int sensorID;
  int row, col;
  uint32_t sourceROF;
  uint32_t nPixels;
};

const std::vector<Expected> expectedClusters{
  {0, 0, 0, 10, 20, 0, 1},
  {1, 1, 1, 11, 21, 0, 3},
  {2, 0, 0, 12, 22, 1, 1},
  {3, 2, 2, 13, 23, 2, 3},
};

// Full parity check between the Fixture and every observable legacy/normalized
// accessor -- reused both right after a successful load and after a later
// call that must have left both owners untouched. `frame` owns the
// normalized measurements; `tf` owns every legacy per-layer compatibility
// structure (Gate 4 B3.1 split).
void verifyFixtureLoaded(const TimeFrame& frame, LegacyTrackerScratch<ITSNLayers>& tf, const Fixture& f, const o2::InteractionRecord& origin, const ROFTimingConfig& timing)
{
  constexpr ClusterSourceId kSourceId{0};

  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{0}).size(), 2u);
  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{1}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{2}).size(), 1u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 0).size() + tf.getUnsortedClustersOnLayer(1, 0).size(), 2u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 1).size(), 1u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(2, 2).size(), 1u);
  for (int l = 3; l < ITSNLayers; ++l) {
    BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{static_cast<uint16_t>(l)}).size(), 0u);
    BOOST_CHECK_EQUAL(tf.getNrof(l), static_cast<int>(f.rofs.size()));
  }

  BOOST_CHECK_EQUAL(tf.getNrof(0), static_cast<int>(f.rofs.size()));
  BOOST_REQUIRE_EQUAL(frame.getNormalizedFrame().getSources().size(), 1u);
  BOOST_CHECK(frame.getNormalizedFrame().getSources()[0].detector == o2::detectors::DetID::ITS);
  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getSources()[0].nROFs, f.rofs.size());
  const auto intervals = frame.getNormalizedFrame().getSourceIntervals(kSourceId);
  BOOST_REQUIRE_EQUAL(intervals.size(), f.rofs.size());
  for (uint32_t r = 0; r < f.rofs.size(); ++r) {
    BOOST_CHECK_EQUAL(intervals[r].sourceROF, r);
    BOOST_CHECK_EQUAL(intervals[r].begin, f.rofs[r].getBCData().differenceInBC(origin));
    BOOST_CHECK_EQUAL(intervals[r].length(), static_cast<uint64_t>(timing.rofLength));
  }

  for (const auto& e : expectedClusters) {
    const tracking::Cluster* legacyCluster = nullptr;
    int clId = -1;
    for (int rof = 0; rof < tf.getNrof(e.layer); ++rof) {
      for (const auto& c : tf.getUnsortedClustersOnLayer(rof, e.layer)) {
        if (tf.getClusterExternalIndex(e.layer, c.clusterId) == static_cast<int>(e.externalIndex)) {
          legacyCluster = &c;
          clId = c.clusterId;
          break;
        }
      }
      if (legacyCluster != nullptr) {
        break;
      }
    }
    BOOST_REQUIRE(legacyCluster != nullptr);

    const SurfaceMeasurement* measurement = nullptr;
    for (const auto& m : frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{static_cast<uint16_t>(e.layer)})) {
      if (m.cluster.index == e.externalIndex) {
        measurement = &m;
        break;
      }
    }
    BOOST_REQUIRE(measurement != nullptr);

    const auto g = expectedGlobal(e.sensorID, e.row, e.col);
    BOOST_CHECK_EQUAL(legacyCluster->xCoordinate, g.x);
    BOOST_CHECK_EQUAL(legacyCluster->yCoordinate, g.y);
    BOOST_CHECK_EQUAL(legacyCluster->zCoordinate, g.z);
    BOOST_CHECK_EQUAL(measurement->global.x, g.x);
    BOOST_CHECK_EQUAL(measurement->global.y, g.y);
    BOOST_CHECK_EQUAL(measurement->global.z, g.z);

    const auto& tfInfo = tf.getClusterTrackingFrameInfo(e.layer, *legacyCluster);
    BOOST_CHECK_EQUAL(tfInfo.xCoordinate, g.x);
    BOOST_CHECK_EQUAL(tfInfo.yCoordinate, g.y);
    BOOST_CHECK_EQUAL(tfInfo.zCoordinate, g.z);
    BOOST_CHECK_EQUAL(tfInfo.xTrackingFrame, static_cast<float>(e.sensorID) + 100.f);
    BOOST_CHECK_EQUAL(tfInfo.alphaTrackingFrame, 0.01f * e.sensorID);
    BOOST_CHECK_EQUAL(tfInfo.positionTrackingFrame[0], static_cast<float>(e.row) + 1.f);
    BOOST_CHECK_EQUAL(tfInfo.positionTrackingFrame[1], static_cast<float>(e.col) + 2.f);
    BOOST_CHECK_EQUAL(measurement->frame.q, tfInfo.xTrackingFrame);
    BOOST_CHECK_EQUAL(measurement->frame.u, tfInfo.positionTrackingFrame[0]);
    BOOST_CHECK_EQUAL(measurement->frame.v, tfInfo.positionTrackingFrame[1]);
    BOOST_CHECK_EQUAL(measurement->frame.frameAngle, tfInfo.alphaTrackingFrame);

    BOOST_CHECK_EQUAL(tfInfo.covarianceTrackingFrame[0], o2::itsmft::ioutils::DefClusError2Row);
    BOOST_CHECK_EQUAL(tfInfo.covarianceTrackingFrame[1], 0.f);
    BOOST_CHECK_EQUAL(tfInfo.covarianceTrackingFrame[2], o2::itsmft::ioutils::DefClusError2Col);
    BOOST_CHECK_EQUAL(measurement->covariance.uu, o2::itsmft::ioutils::DefClusError2Row);
    BOOST_CHECK_EQUAL(measurement->covariance.uv, 0.f);
    BOOST_CHECK_EQUAL(measurement->covariance.vv, o2::itsmft::ioutils::DefClusError2Col);

    BOOST_CHECK_EQUAL(tf.getClusterExternalIndex(e.layer, clId), static_cast<int>(e.externalIndex));
    BOOST_CHECK_EQUAL(measurement->cluster.index, e.externalIndex);
    BOOST_CHECK(measurement->cluster.source == kSourceId);
    BOOST_CHECK(measurement->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::ITS));
    BOOST_CHECK_EQUAL(measurement->sensor.sensor, static_cast<uint32_t>(e.sensorID));
    BOOST_CHECK(measurement->surface == SurfaceId{static_cast<uint16_t>(e.layer)});
    BOOST_CHECK_EQUAL(measurement->sourceROF, e.sourceROF);

    BOOST_CHECK_EQUAL(static_cast<uint32_t>(tf.getClusterSize(e.layer, clId)), e.nPixels);
    BOOST_CHECK_EQUAL(measurement->shape.nPixels, e.nPixels);

    const auto legacyLabels = tf.getClusterLabels(e.layer, clId);
    const auto normalizedLabels = frame.getNormalizedFrame().getLabels(ClusterRef{kSourceId, e.externalIndex});
    BOOST_REQUIRE_EQUAL(legacyLabels.size(), 1u);
    BOOST_REQUIRE_EQUAL(normalizedLabels.size(), 1u);
    BOOST_CHECK(legacyLabels[0] == normalizedLabels[0]);
    BOOST_CHECK(legacyLabels[0] == o2::MCCompLabel(static_cast<int>(e.externalIndex) + 1, 0, 0));
  }
}

} // namespace

// --- A. Wipe lifecycle -------------------------------------------------

BOOST_AUTO_TEST_CASE(WipeClearsNormalizedFrameButPreservesDetId)
{
  const auto catalog = makeITSTestCatalog();
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0};

  TimeFrame frame;
  LegacyTrackerScratch<ITSNLayers> tf;
  std::vector<TrackingParameters> noIterations;
  auto planResult = buildDetectorLayoutSet(catalogView, orderedSurfaces, noIterations);
  BOOST_REQUIRE(planResult.ok());
  const auto plan = std::move(*planResult.layout);

  const auto f = makeFixture();
  const auto result = tf.loadNormalizedSource(frame, decoder, origin, timing, f.clusters, f.patterns, f.rofs, &dict(), &f.labels, o2::detectors::DetID::ITS,
                                              gsl::span<const SurfaceId>{plan.getConfigurationKey().orderedSurfaces}, plan.getSurfaceCatalog());
  BOOST_REQUIRE(result.ok());
  // Sanity: the successful load itself has the expected content, matching
  // the accepted parity coverage in testTimeFrameNormalizedSource.cxx.
  verifyFixtureLoaded(frame, tf, f, origin, timing);

  frame.wipe();

  // --- inspect only freshly obtained normalized accessors/views ---
  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getTotalMeasurements(), 0u);
  BOOST_CHECK(frame.getNormalizedFrame().getSources().empty());
  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getNSurfaces(), 0u);
  for (uint16_t s = 0; s < ITSNLayers; ++s) {
    BOOST_CHECK(frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{s}).empty());
  }
  BOOST_CHECK(frame.getNormalizedFrame().getSourceIntervals(ClusterSourceId{0}).empty());
  BOOST_CHECK(frame.getNormalizedFrame().getLabels(ClusterRef{ClusterSourceId{0}, 0}).empty());
  const auto freshView = frame.getNormalizedFrameView();
  BOOST_CHECK_EQUAL(freshView.nSurfaces, 0u);
  BOOST_CHECK_EQUAL(freshView.nSources, 0u);

  // Gate 4 B3.1: neither owner stores mDetId any more -- the plan lives on
  // `plan` above, entirely outside both TimeFrame and LegacyTrackerScratch,
  // so wipe() has no detector-identity state to preserve or clear.
}

// --- B. Strong exception transactionality ------------------------------

BOOST_AUTO_TEST_CASE(BackfillAllocationFailureLeavesNormalizedAndLegacyStateAtBaseline)
{
  const auto catalog = makeITSTestCatalog();
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0};

  // Declared before `frame`/`tf`, so it is destroyed after both at scope
  // exit: local destruction order is the reverse of declaration order, and
  // the pool-backed bounded_vector members on both owners must be destroyed
  // (returning their memory to the pool) before the pool itself can safely
  // go away.
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  LegacyTrackerScratch<ITSNLayers> tf;
  frame.setMemoryPool(pool);
  tf.setMemoryPool(pool);

  std::vector<TrackingParameters> noIterations;
  auto planResult = buildDetectorLayoutSet(catalogView, orderedSurfaces, noIterations);
  BOOST_REQUIRE(planResult.ok());
  const auto plan = std::move(*planResult.layout);
  const gsl::span<const SurfaceId> planOrderedSurfaces{plan.getConfigurationKey().orderedSurfaces};

  const auto f = makeFixture();
  const auto replacement = makeReplacementFixture();

  // Baseline: a successful normalized load, with real content in both the
  // normalized owner and every legacy compatibility structure.
  const auto baseline = tf.loadNormalizedSource(frame, decoder, origin, timing, f.clusters, f.patterns, f.rofs, &dict(), &f.labels, o2::detectors::DetID::ITS,
                                                planOrderedSurfaces, plan.getSurfaceCatalog());
  BOOST_REQUIRE(baseline.ok());
  BOOST_CHECK_EQUAL(decoder.decodeCount, static_cast<int>(f.clusters.size()));
  verifyFixtureLoaded(frame, tf, f, origin, timing);

  // Remove all BoundedMemoryResource headroom: any further allocation from
  // this pool must now throw.
  pool->setMaxMemory(pool->getUsedMemory());

  // A second, independently valid normalized load of a *distinct* fixture
  // (different sensors/layers, coordinates, pattern shapes, ROF partition
  // and MC labels -- see makeReplacementFixture()): loadSources() itself
  // runs on plain heap storage and succeeds, decoding every one of this
  // fixture's clusters, before the load deterministically throws while
  // building the staged legacy backfill, which allocates from the
  // now-exhausted bounded pool.
  bool threw = false;
  try {
    tf.loadNormalizedSource(frame, decoder, origin, timing, replacement.clusters, replacement.patterns, replacement.rofs, &dict(), &replacement.labels, o2::detectors::DetID::ITS,
                            planOrderedSurfaces, plan.getSurfaceCatalog());
  } catch (const BoundedMemoryResource::MemoryLimitExceeded&) {
    threw = true;
  }
  BOOST_REQUIRE(threw);
  // All of the replacement fixture's clusters were successfully decoded
  // (loadSources() ran to completion) before the staged legacy backfill
  // allocation failed: the decoder was called exactly once per baseline
  // cluster plus once per replacement cluster, never partially.
  BOOST_CHECK_EQUAL(decoder.decodeCount, static_cast<int>(f.clusters.size() + replacement.clusters.size()));

  // Every normalized measurement/source/timing/label on `frame`, and every
  // legacy layer's clusters/tracking information/external indices/sizes/ROF
  // boundaries/label pointers on `tf`, remain exactly at their baseline
  // state -- not even partially replaced by the distinct replacement
  // fixture's data. This is the owner-level load's all-or-nothing contract:
  // a failed stage leaves both live owners unchanged.
  verifyFixtureLoaded(frame, tf, f, origin, timing);
}

// --- C. LegacyTrackerScratch as sole owner of its BoundedMemoryResource -
//
// Exercises the member-destruction-order contract directly (see the
// mExtMemoryPool/mMemoryPool declaration-order comment in
// LegacyTrackerScratch.h): the pool owner is declared before every
// pmr/bounded_vector member, so LegacyTrackerScratch destroys every
// pool-backed vector -- returning its memory to the pool -- before
// releasing its own shared_ptr to that pool. Here the caller's shared_ptr
// is released first, making the scratch the *sole* remaining owner of the
// BoundedMemoryResource while its pool-backed vectors (mUnsortedClusters,
// mTrackingFrameInfo, mClusterExternalIndices, mClusterSize,
// mROFramesClusters) are still populated; a regression in that member order
// would have the scratch free the pool while those vectors still reference
// it, then crash or corrupt memory when they are destroyed. A sanitizer
// build (e.g. ASan) would catch such a regression far more reliably than
// this plain host run can -- a use-after-free here is not guaranteed to
// crash every time -- but this test still documents and exercises the
// ordering contract this correction depends on. (Gate 4 B3.1: these
// pool-backed legacy vectors moved from TimeFrame<NLayers> to
// LegacyTrackerScratch<NLayers>, so this test now targets the scratch's own
// memory-pool ownership; TimeFrame itself no longer holds any of them.)
BOOST_AUTO_TEST_CASE(ScratchOutlivesSoleOwnershipOfItsMemoryPool)
{
  const auto catalog = makeITSTestCatalog();
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0};
  const auto f = makeFixture();

  std::vector<TrackingParameters> noIterations;
  auto planResult = buildDetectorLayoutSet(catalogView, orderedSurfaces, noIterations);
  BOOST_REQUIRE(planResult.ok());
  const auto plan = std::move(*planResult.layout);

  {
    auto pool = std::make_shared<BoundedMemoryResource>();
    TimeFrame frame;
    LegacyTrackerScratch<ITSNLayers> tf;
    tf.setMemoryPool(pool);

    // Populate every pool-backed vector (mUnsortedClusters,
    // mTrackingFrameInfo, mClusterExternalIndices, mClusterSize,
    // mROFramesClusters) with real content.
    const auto result = tf.loadNormalizedSource(frame, decoder, origin, timing, f.clusters, f.patterns, f.rofs, &dict(), &f.labels, o2::detectors::DetID::ITS,
                                                gsl::span<const SurfaceId>{plan.getConfigurationKey().orderedSurfaces}, plan.getSurfaceCatalog());
    BOOST_REQUIRE(result.ok());
    BOOST_REQUIRE(tf.getUnsortedClustersOnLayer(0, 0).size() + tf.getUnsortedClustersOnLayer(1, 0).size() > 0);

    // tf.mMemoryPool holds a second reference at this point (frame never
    // received a copy of `pool` in this test, since the vectors under test
    // all live on the scratch).
    BOOST_REQUIRE_EQUAL(pool.use_count(), 2);
    // Release the caller's reference: tf is now the sole owner of this
    // BoundedMemoryResource, while its pool-backed vectors above are still
    // fully populated and alive.
    pool.reset();
    // Falling off this scope destroys `tf` -- and, per the member-order
    // contract, every pool-backed vector before the pool's own last
    // shared_ptr -- without any of them touching already-freed memory.
    // `frame` is declared before `tf`, so it is destroyed after `tf`, per
    // the Gate 4 B3.1 lifetime contract (TimeFrame outlives every scratch
    // bound to it).
  }
}
