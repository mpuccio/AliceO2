// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 1 compatibility-boundary parity tests (AgentCoordination.md / Wave 1
// exit criteria: "Current single-detector TimeFrames still load and expose
// equivalent clusters"). These instantiate the actual common
// o2::itsmft::tracking::TimeFrame<ITSNLayers>/<MFTNLayers> (the TimeFrame
// already shared with MFT production via CATrackerSpec.cxx, and compiled --
// though not yet workflow-wired -- for ITS) and exercise
// TimeFrame<NLayers>::loadNormalizedSource(), which loads one single-detector
// cluster stream through the normalized MultiSourceFrame owner
// (ITSMFTTracking/MultiSourceLoading.h) and then backfills this same
// TimeFrame's existing legacy compatibility structures (unsorted clusters,
// TrackingFrameInfo, external indices, cluster sizes, ROF boundaries, label
// lookup) from the committed normalized measurements.
//
// Deterministic adapter parity (this test) vs. real decode parity: like
// every other Gate 1 test in this directory (see testMultiSourceLoading.cxx),
// this uses a host-only decoder instead of the production
// GeometryClusterDecoder<DetId>, because no unit test in this component
// depends on a loaded ITS/MFT TGeo geometry singleton. It still reuses the
// real o2::itsmft::ioutils::extractClusterData pattern/dictionary
// consumption path, so cursor bookkeeping is exercised exactly as production
// would; only the geometry transform itself is a deterministic stand-in.
//
// Real decode parity is not re-validated here: it was already confirmed by
// a separate, accepted real-geometry validation that ran
// ITSGeometryClusterDecoder, MFTGeometryClusterDecoder and loadSources() over
// all 7,057 clusters of the retained Gate 0 fixture with zero discrepancies,
// against real ITS/MFT geometry and the recorded CCDB dictionaries. That
// validation covered common/group pattern consumption, ITS legacy-frame
// parity, the explicit separation of MFT's legacy synthetic coordinates from
// normalized disk coordinates, and ITS-only, MFT-only and combined
// normalized loading. This test suite instead covers what that validation
// does not: the common TimeFrame<NLayers> compatibility boundary itself
// (loadNormalizedSource() and its backfill of legacy compatibility
// structures from the normalized owner), which is unaffected by which
// decoder produced the measurements.

#define BOOST_TEST_MODULE ITSMFT TimeFrameNormalizedSource
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
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
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/SparseTrackingTopology.h"
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

// Deterministic, geometry-free stand-in for GeometryClusterDecoder<DetId>:
// sensorID is used directly as the detector-local layer (an identity mapping
// chosen only to keep this fixture simple), and global/frame coordinates are
// pure, independently-recomputable functions of (sensorID, row, col). Pattern
// consumption goes through the real production helper
// (o2::itsmft::ioutils::extractClusterData), so cursor bookkeeping is
// exercised identically to GeometryClusterDecoder. Also records the
// `applySysErrors` value it was called with, so tests can verify
// loadNormalizedSource() propagates it correctly.
class LegacyLikeDecoder final : public ClusterDecoder
{
 public:
  LegacyLikeDecoder(o2::detectors::DetID::ID detector, bool disk) : mDetector(detector), mDisk(disk) {}

  o2::itsmft::ioutils::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    gsl::span<const unsigned char>::iterator& pattIt,
    const TopologyDictionary* dict,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool applySysErrors) const override
  {
    lastApplySysErrors = applySysErrors;
    float sigma2Row{0.f};
    float sigma2Col{0.f};
    ClusterShape shape{};
    o2::itsmft::ioutils::extractClusterData(cluster, pattIt, dict, sigma2Row, sigma2Col, nullptr, &shape);

    o2::itsmft::ioutils::SurfaceMeasurementDecodeResult result;
    const int sensorID = cluster.getSensorID();
    const int layer = sensorID;
    result.layer = layer;
    if (layer < 0 || static_cast<size_t>(layer) >= layerToSurface.size()) {
      return result;
    }
    result.layerMapped = true;
    result.kind = mDisk ? SurfaceKind::Disk : SurfaceKind::Cylinder;

    DecodedCluster decoded{};
    decoded.global = {static_cast<float>(sensorID) * 10.f, static_cast<float>(cluster.getRow()), static_cast<float>(cluster.getCol())};
    decoded.cylinderFrame = {static_cast<float>(sensorID) + 100.f, static_cast<float>(cluster.getRow()) + 1.f, static_cast<float>(cluster.getCol()) + 2.f, 0.01f * sensorID};
    decoded.rowColumnCovariance = {sigma2Row, 0.f, sigma2Col};
    decoded.shape = shape;
    decoded.sensor = static_cast<uint32_t>(sensorID);
    decoded.layer = layer;

    const auto surface = layerToSurface[layer];
    const DetectorSensorId sensor{static_cast<uint32_t>(mDetector), decoded.sensor};
    const ClusterRef clusterRef{source, externalIndex};
    result.measurement = mDisk ? makeDiskSurfaceMeasurement(decoded, sensor, surface, clusterRef, sourceROF)
                               : makeCylinderSurfaceMeasurement(decoded, sensor, surface, clusterRef, sourceROF);
    return result;
  }

  mutable bool lastApplySysErrors{false};

 private:
  o2::detectors::DetID::ID mDetector;
  bool mDisk;
};

const TopologyDictionary& dict()
{
  static const TopologyDictionary d;
  return d;
}

// Explicit (non-grouped, InvalidPatternID) patterns: header {rowSpan,columnSpan}
// plus one bitmap byte. Distinct shapes let a parity check catch a
// misaligned pattern cursor (a shift would corrupt every following cluster's
// shape).
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

// Explicit, hard-coded 7-surface ITS (all cylinder) layout: no
// detector-to-kind inference, just the literal SurfaceKind this test wants.
DetectorLayout makeITSTestLayout()
{
  SparseTrackingTopology topology{static_cast<uint32_t>(ITSNLayers)};
  topology.finalize();
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.reserve(ITSNLayers);
  for (uint16_t i = 0; i < ITSNLayers; ++i) {
    surfaces.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  }
  return DetectorLayout{std::move(surfaces), std::move(topology)};
}

// Explicit, hard-coded 10-surface MFT (all disk) layout.
DetectorLayout makeMFTTestLayout()
{
  SparseTrackingTopology topology{static_cast<uint32_t>(MFTNLayers)};
  topology.finalize();
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.reserve(MFTNLayers);
  for (uint16_t i = 0; i < MFTNLayers; ++i) {
    surfaces.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
  }
  return DetectorLayout{std::move(surfaces), std::move(topology)};
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
  o2::detectors::DetID::ID detector;
  bool disk;
  std::vector<CompClusterExt> clusters;
  std::vector<unsigned char> patterns;
  std::vector<ROFRecord> rofs;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
};

// 4 clusters on layers {0,1,0,2}, partitioned into 3 ROFs: ROF0={c0,c1},
// ROF1={c2}, ROF2={c3}. Row/col distinguish every cluster; patterns
// alternate 1x1/1x3 shapes; one label per cluster.
Fixture makeFixture(o2::detectors::DetID::ID detector, bool disk)
{
  Fixture f;
  f.detector = detector;
  f.disk = disk;
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
  uint32_t nPixels, rowSpan, columnSpan;
};

const std::vector<Expected> expectedClusters{
  {0, 0, 0, 10, 20, 0, 1, 1, 1},
  {1, 1, 1, 11, 21, 0, 3, 1, 3},
  {2, 0, 0, 12, 22, 1, 1, 1, 1},
  {3, 2, 2, 13, 23, 2, 3, 1, 3},
};

template <int NLayers>
void checkParity(const DetectorLayout& layout, const Fixture& f)
{
  BOOST_REQUIRE(layout.valid());
  const auto layerToSurface = identitySurfaces(static_cast<uint16_t>(NLayers));

  LegacyLikeDecoder decoder{f.detector, f.disk};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0}; // rofLength > 0: required, no unusable zero-length default.
  // loadNormalizedSource() always submits exactly one source and fixes its
  // ID to ClusterSourceId{0} internally (loadSources() requires dense,
  // zero-based IDs, so no other value could ever succeed for a single
  // source); this local constant mirrors that fixed value for building
  // expected ClusterRefs below, it is not passed to loadNormalizedSource().
  constexpr ClusterSourceId kSourceId{0};

  TimeFrame<NLayers> tf;
  const auto result = tf.loadNormalizedSource(layout.getView(), layerToSurface, decoder, origin, timing,
                                              f.clusters, f.patterns, f.rofs, &dict(), &f.labels, f.detector);
  BOOST_REQUIRE(result.ok());

  // --- cluster counts per legacy layer == normalized surface (identity layout) ---
  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{0}).size(), 2u); // clusters 0,2
  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{1}).size(), 1u); // cluster 1
  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{2}).size(), 1u); // cluster 3
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 0).size() + tf.getUnsortedClustersOnLayer(1, 0).size(), 2u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 1).size(), 1u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(2, 2).size(), 1u);
  for (int l = 3; l < NLayers; ++l) {
    BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{static_cast<uint16_t>(l)}).size(), 0u);
    BOOST_CHECK_EQUAL(tf.getNrof(l), static_cast<int>(f.rofs.size()));
  }

  // --- per-ROF counts: legacy cumulative table vs. per-source ROF count ---
  BOOST_CHECK_EQUAL(tf.getNrof(0), static_cast<int>(f.rofs.size()));
  BOOST_REQUIRE_EQUAL(tf.getNormalizedFrame().getSources().size(), 1u);
  BOOST_CHECK(tf.getNormalizedFrame().getSources()[0].detector == f.detector);
  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSources()[0].nROFs, f.rofs.size());
  const auto intervals = tf.getNormalizedFrame().getSourceIntervals(kSourceId);
  BOOST_REQUIRE_EQUAL(intervals.size(), f.rofs.size());
  for (uint32_t r = 0; r < f.rofs.size(); ++r) {
    BOOST_CHECK_EQUAL(intervals[r].sourceROF, r);
    BOOST_CHECK_EQUAL(intervals[r].begin, f.rofs[r].getBCData().differenceInBC(origin));
    BOOST_CHECK_EQUAL(intervals[r].length(), static_cast<uint64_t>(timing.rofLength));
  }

  for (const auto& e : expectedClusters) {
    // Find this cluster on its expected layer via the legacy unsorted-cluster
    // accessor (using cumulative ROF boundaries to scan every ROF on that layer).
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

    // Find the matching normalized measurement independently.
    const SurfaceMeasurement* measurement = nullptr;
    for (const auto& m : tf.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{static_cast<uint16_t>(e.layer)})) {
      if (m.cluster.index == e.externalIndex) {
        measurement = &m;
        break;
      }
    }
    BOOST_REQUIRE(measurement != nullptr);

    // --- global position: legacy accessor vs. normalized owner vs. independently expected ---
    const auto g = expectedGlobal(e.sensorID, e.row, e.col);
    BOOST_CHECK_EQUAL(legacyCluster->xCoordinate, g.x);
    BOOST_CHECK_EQUAL(legacyCluster->yCoordinate, g.y);
    BOOST_CHECK_EQUAL(legacyCluster->zCoordinate, g.z);
    BOOST_CHECK_EQUAL(measurement->global.x, g.x);
    BOOST_CHECK_EQUAL(measurement->global.y, g.y);
    BOOST_CHECK_EQUAL(measurement->global.z, g.z);

    // --- tracking-frame coordinates and covariance: legacy TrackingFrameInfo vs. normalized owner ---
    const auto& tfInfo = tf.getClusterTrackingFrameInfo(e.layer, *legacyCluster);
    BOOST_CHECK_EQUAL(tfInfo.xCoordinate, g.x);
    BOOST_CHECK_EQUAL(tfInfo.yCoordinate, g.y);
    BOOST_CHECK_EQUAL(tfInfo.zCoordinate, g.z);
    if (f.disk) {
      // Established synthetic legacy MFT representation: global position +
      // row/column covariance, never disk-frame (z,x,y) coordinates.
      BOOST_CHECK_EQUAL(tfInfo.xTrackingFrame, g.x);
      BOOST_CHECK_EQUAL(tfInfo.alphaTrackingFrame, 0.f);
      BOOST_CHECK_EQUAL(tfInfo.positionTrackingFrame[0], g.y);
      BOOST_CHECK_EQUAL(tfInfo.positionTrackingFrame[1], g.z);
    } else {
      BOOST_CHECK_EQUAL(tfInfo.xTrackingFrame, static_cast<float>(e.sensorID) + 100.f);
      BOOST_CHECK_EQUAL(tfInfo.alphaTrackingFrame, 0.01f * e.sensorID);
      BOOST_CHECK_EQUAL(tfInfo.positionTrackingFrame[0], static_cast<float>(e.row) + 1.f);
      BOOST_CHECK_EQUAL(tfInfo.positionTrackingFrame[1], static_cast<float>(e.col) + 2.f);
      // The normalized owner's cylinder frame carries the identical q/u/v/frameAngle.
      BOOST_CHECK_EQUAL(measurement->frame.q, tfInfo.xTrackingFrame);
      BOOST_CHECK_EQUAL(measurement->frame.u, tfInfo.positionTrackingFrame[0]);
      BOOST_CHECK_EQUAL(measurement->frame.v, tfInfo.positionTrackingFrame[1]);
      BOOST_CHECK_EQUAL(measurement->frame.frameAngle, tfInfo.alphaTrackingFrame);
    }
    // With CompCluster::InvalidPatternID, extractClusterData always returns
    // the fixed default half-pixel covariance regardless of pattern shape or
    // detector -- the same constants the legacy per-detector loaders fall
    // back to.
    BOOST_CHECK_EQUAL(tfInfo.covarianceTrackingFrame[0], o2::itsmft::ioutils::DefClusError2Row);
    BOOST_CHECK_EQUAL(tfInfo.covarianceTrackingFrame[1], 0.f);
    BOOST_CHECK_EQUAL(tfInfo.covarianceTrackingFrame[2], o2::itsmft::ioutils::DefClusError2Col);
    BOOST_CHECK_EQUAL(measurement->covariance.uu, o2::itsmft::ioutils::DefClusError2Row);
    BOOST_CHECK_EQUAL(measurement->covariance.uv, 0.f);
    BOOST_CHECK_EQUAL(measurement->covariance.vv, o2::itsmft::ioutils::DefClusError2Col);

    // --- external indices and source-qualified references ---
    BOOST_CHECK_EQUAL(tf.getClusterExternalIndex(e.layer, clId), static_cast<int>(e.externalIndex));
    BOOST_CHECK_EQUAL(measurement->cluster.index, e.externalIndex);
    BOOST_CHECK(measurement->cluster.source == kSourceId);
    BOOST_CHECK(measurement->sensor.detector == static_cast<uint32_t>(f.detector));
    BOOST_CHECK_EQUAL(measurement->sensor.sensor, static_cast<uint32_t>(e.sensorID));
    BOOST_CHECK(measurement->surface == SurfaceId{static_cast<uint16_t>(e.layer)});
    BOOST_CHECK_EQUAL(measurement->sourceROF, e.sourceROF);

    // --- cluster shape / sizes: legacy vs. normalized, explicit pattern consumed exactly once ---
    BOOST_CHECK_EQUAL(static_cast<uint32_t>(tf.getClusterSize(e.layer, clId)), e.nPixels);
    BOOST_CHECK_EQUAL(measurement->shape.nPixels, e.nPixels);
    BOOST_CHECK_EQUAL(measurement->shape.rowSpan, e.rowSpan);
    BOOST_CHECK_EQUAL(measurement->shape.columnSpan, e.columnSpan);

    // --- labels: legacy lookup vs. normalized ClusterRef lookup ---
    const auto legacyLabels = tf.getClusterLabels(e.layer, clId);
    const auto normalizedLabels = tf.getNormalizedFrame().getLabels(ClusterRef{kSourceId, e.externalIndex});
    BOOST_REQUIRE_EQUAL(legacyLabels.size(), 1u);
    BOOST_REQUIRE_EQUAL(normalizedLabels.size(), 1u);
    BOOST_CHECK(legacyLabels[0] == normalizedLabels[0]);
    BOOST_CHECK(legacyLabels[0] == o2::MCCompLabel(static_cast<int>(e.externalIndex) + 1, 0, 0));
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(ITSTimeFrameNormalizedSourceParity)
{
  checkParity<ITSNLayers>(makeITSTestLayout(), makeFixture(o2::detectors::DetID::ITS, false));
}

BOOST_AUTO_TEST_CASE(MFTTimeFrameNormalizedSourceParity)
{
  checkParity<MFTNLayers>(makeMFTTestLayout(), makeFixture(o2::detectors::DetID::MFT, true));
}

BOOST_AUTO_TEST_CASE(EmptyInputsAreLegalForBothDetectors)
{
  {
    const auto layout = makeITSTestLayout();
    BOOST_REQUIRE(layout.valid());
    const auto layerToSurface = identitySurfaces(ITSNLayers);
    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
    TimeFrame<ITSNLayers> tf;
    const auto result = tf.loadNormalizedSource(layout.getView(), layerToSurface, decoder, {0, 0},
                                                ROFTimingConfig{40, 0, 0, 0}, {}, {}, {}, &dict(), nullptr, o2::detectors::DetID::ITS);
    BOOST_CHECK(result.ok());
    BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getTotalMeasurements(), 0u);
    BOOST_REQUIRE_EQUAL(tf.getNormalizedFrame().getSources().size(), 1u);
    BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSources()[0].nROFs, 0u);
    BOOST_CHECK_EQUAL(tf.getNrof(0), 0);
  }
  {
    const auto layout = makeMFTTestLayout();
    BOOST_REQUIRE(layout.valid());
    const auto layerToSurface = identitySurfaces(MFTNLayers);
    LegacyLikeDecoder decoder{o2::detectors::DetID::MFT, true};
    TimeFrame<MFTNLayers> tf;
    const auto result = tf.loadNormalizedSource(layout.getView(), layerToSurface, decoder, {0, 0},
                                                ROFTimingConfig{40, 0, 0, 0}, {}, {}, {}, &dict(), nullptr, o2::detectors::DetID::MFT);
    BOOST_CHECK(result.ok());
    BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getTotalMeasurements(), 0u);
    BOOST_CHECK_EQUAL(tf.getNrof(0), 0);
  }
}

BOOST_AUTO_TEST_CASE(FailedNormalizedLoadLeavesBothRepresentationsUnchanged)
{
  const auto layout = makeITSTestLayout();
  BOOST_REQUIRE(layout.valid());
  const auto layerToSurface = identitySurfaces(ITSNLayers);
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const ROFTimingConfig timing{40, 0, 0, 0};

  // First, a valid baseline load to give the TimeFrame real content in both
  // the normalized owner and the legacy compatibility structures.
  const std::vector<CompClusterExt> goodClusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto goodPatterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> goodRofs{ROFRecord{{0, 0}, 0, 0, 1}};

  TimeFrame<ITSNLayers> tf;
  const auto baseline = tf.loadNormalizedSource(layout.getView(), layerToSurface, decoder, {0, 0},
                                                timing, goodClusters, goodPatterns, goodRofs, &dict(), nullptr, o2::detectors::DetID::ITS);
  BOOST_REQUIRE(baseline.ok());
  BOOST_REQUIRE_EQUAL(tf.getNormalizedFrame().getTotalMeasurements(), 1u);
  BOOST_REQUIRE_EQUAL(tf.getUnsortedClustersOnLayer(0, 0).size(), 1u);
  BOOST_REQUIRE_EQUAL(tf.getNrof(0), 1);

  // Then a malformed ROF partition (leading gap): firstEntry != 0.
  const std::vector<CompClusterExt> badClusters{{2, 2, CompCluster::InvalidPatternID, 0}};
  const auto badPatterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> badRofs{ROFRecord{{0, 0}, 0, 1, 1}}; // gap: cluster 0 unreferenced

  const auto failed = tf.loadNormalizedSource(layout.getView(), layerToSurface, decoder, {0, 0},
                                              timing, badClusters, badPatterns, badRofs, &dict(), nullptr, o2::detectors::DetID::ITS);
  BOOST_CHECK(!failed.ok());
  BOOST_CHECK(failed.error == MultiSourceLoadError::InvalidROFRange);

  // Both the normalized owner and the legacy compatibility structures retain
  // exactly their pre-failure (baseline) content.
  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getTotalMeasurements(), 1u);
  BOOST_REQUIRE_EQUAL(tf.getNormalizedFrame().getSources().size(), 1u);
  BOOST_CHECK_EQUAL(tf.getNormalizedFrame().getSources()[0].nROFs, 1u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 0).size(), 1u);
  BOOST_CHECK_EQUAL(tf.getNrof(0), 1);
  BOOST_CHECK_EQUAL(tf.getClusterExternalIndex(0, 0), 0);
}

// TimeFrame<NLayers>::loadROFrameData() calls loadClusterTrackingFrameInfo()
// with its own default applySysErrors=true; loadNormalizedSource() must match
// that default (covariance compatibility) and must also honor and propagate
// an explicit override, since callers may deliberately want the
// GeometryClusterDecoder sys-error convention turned off.
BOOST_AUTO_TEST_CASE(ApplySysErrorsDefaultsTrueAndPropagatesToTheDecoder)
{
  const auto layout = makeITSTestLayout();
  BOOST_REQUIRE(layout.valid());
  const auto layerToSurface = identitySurfaces(ITSNLayers);
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
  const ROFTimingConfig timing{40, 0, 0, 0};

  {
    // No explicit argument: must match loadROFrameData()'s own default.
    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
    TimeFrame<ITSNLayers> tf;
    const auto result = tf.loadNormalizedSource(layout.getView(), layerToSurface, decoder, {0, 0}, timing,
                                                clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS);
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK(decoder.lastApplySysErrors);
  }
  {
    // Explicit false is honored and reaches the decoder unchanged.
    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
    TimeFrame<ITSNLayers> tf;
    const auto result = tf.loadNormalizedSource(layout.getView(), layerToSurface, decoder, {0, 0}, timing,
                                                clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS, false);
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK(!decoder.lastApplySysErrors);
  }
}
