// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 2 correction: TimeFrame lifecycle and transactional legacy backfill.
//
// A. Wipe lifecycle: TimeFrame<NLayers>::wipe() must unconditionally clear
//    the normalized owner (mNormalizedFrame) associated by
//    loadNormalizedSource() -- every normalized accessor obtained *after*
//    wipe() must report empty/zero content -- while leaving catalog/layout
//    ownership, the required layout configuration, the required geometry
//    epoch and mDetId completely untouched. This test never dereferences a
//    view obtained before wipe(): every post-wipe check re-obtains its
//    accessor.
//
// B. Strong exception transactionality: loadNormalizedSource() stages its
//    entire legacy backfill (unsorted clusters, TrackingFrameInfo, external
//    indices, cluster sizes, ROF boundaries, label pointers) on local
//    owners before touching any live TimeFrame member. If a
//    BoundedMemoryResource allocation fails while building that staged
//    backfill -- after loadSources() has already committed a valid
//    normalized decode into its own *local* scratch owner -- the thrown
//    MemoryLimitExceeded must leave mDetId, the normalized owner and every
//    legacy layer's compatibility structures exactly at their pre-call
//    baseline.

#define BOOST_TEST_MODULE ITSMFT TimeFrame lifecycle
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <memory>
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
#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"
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
    return result;
  }

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

class FakeCatalogProvider final : public DetectorSurfaceCatalogProvider
{
 public:
  explicit FakeCatalogProvider(std::vector<SurfaceDescriptor> catalog) : mCatalog{std::move(catalog)} {}

  DetectorSurfaceCatalogResult buildCatalog(const DetectorSurfaceCatalogRequest&) const final
  {
    return {mCatalog, DetectorSurfaceCatalogError::None};
  }

  std::vector<SurfaceDescriptor> mCatalog;
};

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
// call that must have left the TimeFrame untouched.
void verifyFixtureLoaded(TimeFrame<ITSNLayers>& tf, const Fixture& f, const o2::InteractionRecord& origin, const ROFTimingConfig& timing)
{
  constexpr ClusterSourceId kSourceId{0};
  BOOST_CHECK(tf.getDetId() == o2::detectors::DetID::ITS);

  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{0}).size(), 2u);
  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{1}).size(), 1u);
  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{2}).size(), 1u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 0).size() + tf.getUnsortedClustersOnLayer(1, 0).size(), 2u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 1).size(), 1u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(2, 2).size(), 1u);
  for (int l = 3; l < ITSNLayers; ++l) {
    BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{static_cast<uint16_t>(l)}).size(), 0u);
    BOOST_CHECK_EQUAL(tf.getNrof(l), static_cast<int>(f.rofs.size()));
  }

  BOOST_CHECK_EQUAL(tf.getNrof(0), static_cast<int>(f.rofs.size()));
  BOOST_REQUIRE_EQUAL(tf.getNormalizedFrame().getSources().size(), 1u);
  BOOST_CHECK(tf.getNormalizedFrame().getSources()[0].detector == o2::detectors::DetID::ITS);
  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSources()[0].nROFs, f.rofs.size());
  const auto intervals = tf.getNormalizedFrame().getSourceIntervals(kSourceId);
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
    for (const auto& m : tf.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{static_cast<uint16_t>(e.layer)})) {
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
    const auto normalizedLabels = tf.getNormalizedFrame().getLabels(ClusterRef{kSourceId, e.externalIndex});
    BOOST_REQUIRE_EQUAL(legacyLabels.size(), 1u);
    BOOST_REQUIRE_EQUAL(normalizedLabels.size(), 1u);
    BOOST_CHECK(legacyLabels[0] == normalizedLabels[0]);
    BOOST_CHECK(legacyLabels[0] == o2::MCCompLabel(static_cast<int>(e.externalIndex) + 1, 0, 0));
  }
}

} // namespace

// --- A. Wipe lifecycle -------------------------------------------------

BOOST_AUTO_TEST_CASE(WipeClearsNormalizedFrameButPreservesLayoutConfigurationEpochAndDetId)
{
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  FakeCatalogProvider provider{makeITSTestCatalog()};
  const DetectorSurfaceCatalogRequest catalogRequest{o2::detectors::DetID::ITS, SurfaceId{0}, ITSNLayers};
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0};

  TimeFrame<ITSNLayers> tf;
  std::vector<TrackingParameters> noIterations;
  BOOST_REQUIRE(tf.ensureDetectorLayouts(&provider, catalogRequest, orderedSurfaces, TransitionPolicyTag::CylinderCylinder, noIterations).ok());

  const auto f = makeFixture();
  const auto result = tf.loadNormalizedSource(decoder, origin, timing, f.clusters, f.patterns, f.rofs, &dict(), &f.labels, o2::detectors::DetID::ITS);
  BOOST_REQUIRE(result.ok());
  // Sanity: the successful load itself has the expected content, matching
  // the accepted parity coverage in testTimeFrameNormalizedSource.cxx.
  verifyFixtureLoaded(tf, f, origin, timing);

  // Record catalog/layout/configuration/epoch/mDetId *before* wipe().
  BOOST_REQUIRE(tf.getSurfaceCatalog() != nullptr);
  const auto catalogSizeBefore = tf.getSurfaceCatalog()->size();
  BOOST_REQUIRE(tf.getDetectorLayouts() != nullptr);
  const auto configKeyBefore = tf.getDetectorLayouts()->getConfigurationKey();
  const auto epochBefore = tf.getRequiredDetectorGeometryEpoch();
  const auto detIdBefore = tf.getDetId();
  BOOST_REQUIRE(tf.detectorLayoutsCurrent());

  tf.wipe();

  // --- inspect only freshly obtained normalized accessors/views ---
  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getTotalMeasurements(), 0u);
  BOOST_CHECK(tf.getNormalizedFrame().getSources().empty());
  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getNSurfaces(), 0u);
  for (uint16_t s = 0; s < ITSNLayers; ++s) {
    BOOST_CHECK(tf.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{s}).empty());
  }
  BOOST_CHECK(tf.getNormalizedFrame().getSourceIntervals(ClusterSourceId{0}).empty());
  BOOST_CHECK(tf.getNormalizedFrame().getLabels(ClusterRef{ClusterSourceId{0}, 0}).empty());
  const auto freshView = tf.getNormalizedFrameView();
  BOOST_CHECK_EQUAL(freshView.nMeasurements, 0u);
  BOOST_CHECK_EQUAL(freshView.nSurfaces, 0u);
  BOOST_CHECK_EQUAL(freshView.nSources, 0u);

  // --- catalog/layout ownership remains stored, current and unchanged ---
  BOOST_CHECK(tf.hasStoredDetectorLayouts());
  BOOST_CHECK(tf.detectorLayoutsCurrent());
  BOOST_REQUIRE(tf.getSurfaceCatalog() != nullptr);
  BOOST_CHECK_EQUAL(tf.getSurfaceCatalog()->size(), catalogSizeBefore);
  BOOST_REQUIRE(tf.getDetectorLayouts() != nullptr);
  BOOST_CHECK(tf.getDetectorLayouts()->getConfigurationKey() == configKeyBefore);

  // --- configuration, epoch and mDetId remain unchanged ---
  BOOST_CHECK_EQUAL(tf.getRequiredDetectorGeometryEpoch(), epochBefore);
  BOOST_CHECK(tf.getDetId() == detIdBefore);
}

// --- B. Strong exception transactionality ------------------------------

BOOST_AUTO_TEST_CASE(BackfillAllocationFailureLeavesNormalizedAndLegacyStateAtBaseline)
{
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  FakeCatalogProvider provider{makeITSTestCatalog()};
  const DetectorSurfaceCatalogRequest catalogRequest{o2::detectors::DetID::ITS, SurfaceId{0}, ITSNLayers};
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0};

  // Declared before `tf`, so it is destroyed after `tf` at scope exit: local
  // destruction order is the reverse of declaration order, and TimeFrame's
  // pool-backed bounded_vector members must be destroyed (returning their
  // memory to the pool) before the pool itself can safely go away.
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame<ITSNLayers> tf;
  tf.setMemoryPool(pool);

  std::vector<TrackingParameters> noIterations;
  BOOST_REQUIRE(tf.ensureDetectorLayouts(&provider, catalogRequest, orderedSurfaces, TransitionPolicyTag::CylinderCylinder, noIterations).ok());

  const auto f = makeFixture();

  // Baseline: a successful normalized load, with real content in both the
  // normalized owner and every legacy compatibility structure.
  const auto baseline = tf.loadNormalizedSource(decoder, origin, timing, f.clusters, f.patterns, f.rofs, &dict(), &f.labels, o2::detectors::DetID::ITS);
  BOOST_REQUIRE(baseline.ok());
  verifyFixtureLoaded(tf, f, origin, timing);

  // Remove all BoundedMemoryResource headroom: any further allocation from
  // this pool must now throw.
  pool->setMaxMemory(pool->getUsedMemory());

  // A second, otherwise-valid normalized load (loadSources() itself runs on
  // plain heap storage and succeeds) that deterministically throws while
  // building the staged legacy backfill, which allocates from the
  // now-exhausted bounded pool.
  bool threw = false;
  try {
    tf.loadNormalizedSource(decoder, origin, timing, f.clusters, f.patterns, f.rofs, &dict(), &f.labels, o2::detectors::DetID::ITS);
  } catch (const BoundedMemoryResource::MemoryLimitExceeded&) {
    threw = true;
  }
  BOOST_REQUIRE(threw);

  // mDetId, every normalized measurement/source/timing/label, and every
  // legacy layer's clusters/tracking information/external indices/sizes/ROF
  // boundaries/label pointers remain exactly at their baseline state.
  verifyFixtureLoaded(tf, f, origin, timing);
}
