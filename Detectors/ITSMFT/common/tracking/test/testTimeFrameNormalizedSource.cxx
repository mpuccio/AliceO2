// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 1/2 compatibility-boundary parity tests (AgentCoordination.md / Wave 1
// exit criteria: "Current single-detector TimeFrames still load and expose
// equivalent clusters"), updated for the Gate 4 B3.1 TimeFrame /
// TimeFrameScratch split. These instantiate the actual common
// o2::itsmft::tracking::TimeFrame (the detector-neutral normalized owner)
// together with o2::itsmft::tracking::TimeFrameScratch/
// <MFTNLayers> (the legacy per-detector compatibility scratch, already
// shared with MFT production via CATrackerSpec.cxx, and compiled -- though
// not yet workflow-wired -- for ITS) and exercise
// TimeFrameScratch::loadNormalizedSource(TimeFrame&, ...), the
// owner-level load operation which loads one single-detector cluster stream
// through TimeFrame event measurement storage (ITSMFTTracking/
// IOUtils.h) and then backfills the scratch's existing legacy
// compatibility structures (unsorted clusters, TrackingFrameInfo, external
// indices, cluster sizes, ROF boundaries, label lookup) from the committed
// normalized measurements.
//
// Gate 2 (Slice B3): loadNormalizedSource() no longer accepts an externally
// supplied topology view or layer-to-surface mapping. TimeFrame owns no
// catalog/plan of its own; every test passes the immutable catalog/order
// explicitly to loadNormalizedSource().
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
// does not: the common TimeFrame/TimeFrameScratch compatibility
// boundary itself (loadNormalizedSource() and its backfill of legacy
// compatibility structures from the normalized owner), which is unaffected
// by which decoder produced the measurements.

#define BOOST_TEST_MODULE ITSMFT TimeFrameNormalizedSource
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <memory>
#include <optional>
#include <tuple>
#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/SurfaceLayout.h"
#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{
SurfaceLayout catalogGraph(SurfaceCatalogView catalog, gsl::span<const SurfaceId> ordered)
{
  return SurfaceLayout{gsl::span<const SurfaceDescriptor>{catalog.surfaces, catalog.nSurfaces},
                       makeSurfaceLayoutChain(ordered)};
}

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

  o2::itsmft::tracking::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const TopologyDictionary* dict,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool applySysErrors) const override
  {
    lastApplySysErrors = applySysErrors;
    const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dict);
    if (!clusterData.ok()) {
      o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
      result.error = clusterData.error;
      return result;
    }

    o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
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
    decoded.rowColumnCovariance = {clusterData.sig2Row, 0.f, clusterData.sig2Col};
    decoded.shape = clusterData.shape;
    decoded.sensor = static_cast<uint32_t>(sensorID);
    decoded.layer = layer;

    const auto surface = layerToSurface[layer];
    const DetectorSensorId sensor{static_cast<uint32_t>(mDetector), decoded.sensor};
    const ClusterRef clusterRef{source, externalIndex};
    result = mDisk ? makeDiskMeasurementDecodeResult(decoded, sensor, surface, clusterRef, sourceROF)
                   : makeCylinderMeasurementDecodeResult(decoded, sensor, surface, clusterRef, sourceROF);
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

void configureScratchFromPlan(TimeFrameScratch& scratch, std::size_t nOwnedSurfaces)
{
  scratch.setMemoryPool(std::make_shared<o2::its::BoundedMemoryResource>());
  scratch.adoptPlan(nOwnedSurfaces, 0, 0);
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

// Explicit, hard-coded 7-surface ITS (all cylinder) canonical catalog: no
// detector-to-kind inference, just the literal SurfaceKind this test wants.
std::vector<SurfaceDescriptor> makeITSTestCatalog()
{
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.reserve(ITSNLayers);
  for (uint16_t i = 0; i < ITSNLayers; ++i) {
    surfaces.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  }
  return surfaces;
}

// Explicit, hard-coded 10-surface MFT (all disk) canonical catalog.
std::vector<SurfaceDescriptor> makeMFTTestCatalog()
{
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.reserve(MFTNLayers);
  for (uint16_t i = 0; i < MFTNLayers; ++i) {
    surfaces.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
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

BOOST_AUTO_TEST_CASE(combined_owner_load_keeps_detector_backfills_separate)
{
  auto its = makeFixture(o2::detectors::DetID::ITS, false);
  auto mft = makeFixture(o2::detectors::DetID::MFT, true);
  std::vector<SurfaceDescriptor> catalog;
  for (uint16_t layer = 0; layer < ITSNLayers; ++layer) {
    catalog.push_back({SurfaceId{layer}, layer, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  }
  for (uint16_t layer = 0; layer < MFTNLayers; ++layer) {
    catalog.push_back({SurfaceId{static_cast<uint16_t>(ITSNLayers + layer)}, layer,
                       static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
  }
  const std::array<SurfaceId, ITSNLayers> itsSurfaces{SurfaceId{0}, SurfaceId{1}, SurfaceId{2}, SurfaceId{3}, SurfaceId{4}, SurfaceId{5}, SurfaceId{6}};
  std::array<SurfaceId, MFTNLayers> mftSurfaces{};
  for (uint16_t layer = 0; layer < MFTNLayers; ++layer) {
    mftSurfaces[layer] = SurfaceId{static_cast<uint16_t>(ITSNLayers + layer)};
  }
  LegacyLikeDecoder itsDecoder{o2::detectors::DetID::ITS, false};
  LegacyLikeDecoder mftDecoder{o2::detectors::DetID::MFT, true};
  const ROFTimingConfig timing{40, 0, 0, 0};
  ClusterSourceInput itsSource{ClusterSourceId{0}, o2::detectors::DetID::ITS, its.clusters, its.patterns, its.rofs, &dict(), &its.labels,
                               itsSurfaces, timing, &itsDecoder};
  ClusterSourceInput mftSource{ClusterSourceId{1}, o2::detectors::DetID::MFT, mft.clusters, mft.patterns, mft.rofs, &dict(), &mft.labels,
                               mftSurfaces, timing, &mftDecoder};
  TimeFrame frame;
  const SurfaceCatalogView view{catalog.data(), static_cast<uint32_t>(catalog.size())};
  std::vector<SurfaceId> combinedSurfaces{itsSurfaces.begin(), itsSurfaces.end()};
  combinedSurfaces.insert(combinedSurfaces.end(), mftSurfaces.begin(), mftSurfaces.end());
  auto layout = catalogGraph(view, combinedSurfaces);
  std::vector<SurfaceLayout> layouts;
  layouts.push_back(std::move(layout));
  std::vector<TrackingParameters> parameters(1);
  std::vector<TrackingWorkspaceCapacity> capacities{{combinedSurfaces.size(), 0, 0}};
  BOOST_REQUIRE(frame.commitConfiguration(std::move(layouts), std::move(parameters),
                                          std::move(capacities), std::make_shared<BoundedMemoryResource>()));
  const std::array<ClusterSourceInput, 2> sources{itsSource, mftSource};
  BOOST_REQUIRE(loadTimeFrameSources(frame, gsl::span<const ClusterSourceInput>{sources}, view, {50, 5}).ok());
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0}).size(), 2u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{static_cast<uint16_t>(ITSNLayers)}).size(), 2u);
  BOOST_CHECK_EQUAL(frame.getWorkspace().getTotalClusters(), static_cast<int>(its.clusters.size() + mft.clusters.size()));
  BOOST_CHECK(frame.getGlobalMeasurements(SurfaceId{0})[0].cluster.source == ClusterSourceId{0});
  BOOST_CHECK(frame.getGlobalMeasurements(SurfaceId{static_cast<uint16_t>(ITSNLayers)})[0].cluster.source == ClusterSourceId{1});
  BOOST_CHECK(frame.getWorkspace().getSurfaceSource(0) == ClusterSourceId{0});
  BOOST_CHECK(frame.getWorkspace().getSurfaceSource(ITSNLayers) == ClusterSourceId{1});

  // Frame reset clears the workspace and normalized ownership together.
  frame.resetTimeFrame();
  BOOST_CHECK(frame.getWorkspace().empty());
  BOOST_CHECK(!frame.getWorkspace().getSurfaceSource(0));

  // A malformed replacement is rejected before the no-throw three-owner
  // commit; the still-live MFT scratch and shared normalized owner survive.
  auto malformedMFT = mftSource;
  malformedMFT.patterns = {};
  const std::array<ClusterSourceInput, 2> retrySources{itsSource, malformedMFT};
  BOOST_CHECK(!loadTimeFrameSources(frame, gsl::span<const ClusterSourceInput>{retrySources}, view, {50, 5}).ok());
  BOOST_CHECK(frame.getWorkspace().empty());
  frame.resetTimeFrame();
  BOOST_CHECK(frame.getWorkspace().empty());
}

template <int NLayers>
void checkParity(std::vector<SurfaceDescriptor> catalog, const Fixture& f)
{
  const auto orderedSurfaces = identitySurfaces(static_cast<uint16_t>(NLayers));
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};

  LegacyLikeDecoder decoder{f.detector, f.disk};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0}; // rofLength > 0: required, no unusable zero-length default.
  // loadNormalizedSource() always submits exactly one source and fixes its
  // ID to ClusterSourceId{0} internally (loadSources() requires dense,
  // zero-based IDs, so no other value could ever succeed for a single
  // source); this local constant mirrors that fixed value for building
  // expected ClusterRefs below, it is not passed to loadNormalizedSource().
  constexpr ClusterSourceId kSourceId{0};

  TimeFrame frame;

  TimeFrameScratch tf;
  std::vector<TrackingParameters> noIterations;
  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());

  const auto result = tf.loadNormalizedSource(frame, decoder, origin, timing,
                                              f.clusters, f.patterns, f.rofs, &dict(), &f.labels, f.detector,
                                              gsl::span<const SurfaceId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog());
  BOOST_REQUIRE(result.ok());

  // --- cluster counts per legacy layer == normalized surface (identity layout) ---
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0}).size(), 2u); // clusters 0,2
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{1}).size(), 1u); // cluster 1
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{2}).size(), 1u); // cluster 3
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 0).size() + tf.getUnsortedClustersOnLayer(1, 0).size(), 2u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 1).size(), 1u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(2, 2).size(), 1u);

  // loadROFrameData() has always sized mNTrackletsPerCluster/
  // mNTrackletsPerClusterSum (both 2-element, not per-layer, consumed by
  // computeTrackletsPerROFScans() on the CA hot path) from layer 1's total
  // cluster count; loadNormalizedSource() must do the same, or a TF loaded
  // through this path would silently carry stale/empty tracklet-per-cluster
  // scratch into tracking.
  const auto nClustersLayer1 = tf.mUnsortedClusters[1].size();
  for (int i = 0; i < 2; ++i) {
    BOOST_CHECK_EQUAL(tf.mNTrackletsPerCluster[i].size(), nClustersLayer1);
    BOOST_CHECK_EQUAL(tf.mNTrackletsPerClusterSum[i].size(), nClustersLayer1 + 1);
  }
  for (int l = 3; l < NLayers; ++l) {
    BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{static_cast<uint16_t>(l)}).size(), 0u);
    BOOST_CHECK_EQUAL(tf.getNrof(l), static_cast<int>(f.rofs.size()));
  }

  // --- per-ROF counts: legacy cumulative table vs. per-source ROF count ---
  BOOST_CHECK_EQUAL(tf.getNrof(0), static_cast<int>(f.rofs.size()));
  BOOST_CHECK_EQUAL(frame.getSourceIntervals(kSourceId).size(), f.rofs.size());
  const auto intervals = frame.getSourceIntervals(kSourceId);
  BOOST_REQUIRE_EQUAL(intervals.size(), f.rofs.size());
  for (uint32_t r = 0; r < f.rofs.size(); ++r) {
    BOOST_CHECK_EQUAL(intervals[r].sourceROF, r);
    BOOST_CHECK_EQUAL(intervals[r].begin, f.rofs[r].getBCData().differenceInBC(origin));
    BOOST_CHECK_EQUAL(intervals[r].length(), static_cast<uint64_t>(timing.rofLength));
  }

  for (const auto& e : expectedClusters) {
    // Find this cluster on its expected layer via the legacy unsorted-cluster
    // accessor (using cumulative ROF boundaries to scan every ROF on that layer).
    const o2::its::Cluster* legacyCluster = nullptr;
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
    const GlobalMeasurement* globalMeasurement = nullptr;
    const SurfaceMeasurement* measurement = nullptr;
    const auto surface = SurfaceId{static_cast<uint16_t>(e.layer)};
    const auto globals = frame.getGlobalMeasurements(surface);
    const auto locals = frame.getSurfaceMeasurements(surface);
    for (size_t index = 0; index < globals.size(); ++index) {
      if (globals[index].cluster.index == e.externalIndex) {
        globalMeasurement = &globals[index];
        measurement = &locals[index];
        break;
      }
    }
    BOOST_REQUIRE(globalMeasurement != nullptr);
    BOOST_REQUIRE(measurement != nullptr);

    // --- global position: legacy accessor vs. normalized owner vs. independently expected ---
    const auto g = expectedGlobal(e.sensorID, e.row, e.col);
    BOOST_CHECK_EQUAL(legacyCluster->xCoordinate, g.x);
    BOOST_CHECK_EQUAL(legacyCluster->yCoordinate, g.y);
    BOOST_CHECK_EQUAL(legacyCluster->zCoordinate, g.z);
    BOOST_CHECK_EQUAL(globalMeasurement->position.x, g.x);
    BOOST_CHECK_EQUAL(globalMeasurement->position.y, g.y);
    BOOST_CHECK_EQUAL(globalMeasurement->position.z, g.z);

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
    BOOST_CHECK_EQUAL(globalMeasurement->cluster.index, e.externalIndex);
    BOOST_CHECK(globalMeasurement->cluster.source == kSourceId);
    BOOST_CHECK(globalMeasurement->sensor.detector == static_cast<uint32_t>(f.detector));
    BOOST_CHECK_EQUAL(globalMeasurement->sensor.sensor, static_cast<uint32_t>(e.sensorID));
    BOOST_CHECK(globalMeasurement->surface == SurfaceId{static_cast<uint16_t>(e.layer)});
    BOOST_CHECK_EQUAL(globalMeasurement->sourceROF, e.sourceROF);

    // --- cluster shape / sizes: legacy vs. normalized, explicit pattern consumed exactly once ---
    BOOST_CHECK_EQUAL(static_cast<uint32_t>(tf.getClusterSize(e.layer, clId)), e.nPixels);
    BOOST_CHECK_EQUAL(globalMeasurement->shape.nPixels, e.nPixels);
    BOOST_CHECK_EQUAL(globalMeasurement->shape.rowSpan, e.rowSpan);
    BOOST_CHECK_EQUAL(globalMeasurement->shape.columnSpan, e.columnSpan);

    // --- labels: legacy lookup vs. normalized ClusterRef lookup ---
    const auto legacyLabels = tf.getClusterLabels(e.layer, clId);
    const auto normalizedLabels = frame.getLabels(ClusterRef{kSourceId, e.externalIndex});
    BOOST_REQUIRE_EQUAL(legacyLabels.size(), 1u);
    BOOST_REQUIRE_EQUAL(normalizedLabels.size(), 1u);
    BOOST_CHECK(legacyLabels[0] == normalizedLabels[0]);
    BOOST_CHECK(legacyLabels[0] == o2::MCCompLabel(static_cast<int>(e.externalIndex) + 1, 0, 0));
  }
}

// A caller can no longer supply an unrelated topology view or
// layer-to-surface mapping by API construction: taking the current member
// function's address is unambiguous (there is exactly one overload), so
// checking its invocability against the removed Gate 1 argument list proves
// the parameters are gone from the signature itself, not merely unused.
static_assert(!std::is_invocable_v<decltype(&TimeFrameScratch::loadNormalizedSource),
                                   TimeFrameScratch&,
                                   TimeFrame&,
                                   const TraversalTopologyView&,
                                   gsl::span<const SurfaceId>,
                                   const ClusterDecoder&,
                                   const o2::InteractionRecord&,
                                   const ROFTimingConfig&,
                                   gsl::span<const o2::itsmft::CompClusterExt>,
                                   gsl::span<const unsigned char>,
                                   gsl::span<const o2::itsmft::ROFRecord>,
                                   const o2::itsmft::TopologyDictionary*,
                                   const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*,
                                   o2::detectors::DetID::ID>,
              "loadNormalizedSource must no longer accept an externally supplied layout or layer-to-surface mapping");

} // namespace

BOOST_AUTO_TEST_CASE(ITSTimeFrameNormalizedSourceParity)
{
  checkParity<ITSNLayers>(makeITSTestCatalog(), makeFixture(o2::detectors::DetID::ITS, false));
}

BOOST_AUTO_TEST_CASE(MFTTimeFrameNormalizedSourceParity)
{
  checkParity<MFTNLayers>(makeMFTTestCatalog(), makeFixture(o2::detectors::DetID::MFT, true));
}

BOOST_AUTO_TEST_CASE(EmptyInputsAreLegalForBothDetectors)
{
  {
    const auto orderedSurfaces = identitySurfaces(ITSNLayers);
    const auto catalog = makeITSTestCatalog();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
    TimeFrame frame;
    TimeFrameScratch tf;
    std::vector<TrackingParameters> noIterations;
    const auto plan = catalogGraph(catalogView, identitySurfaces(ITSNLayers));
    configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());
    const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0},
                                                ROFTimingConfig{40, 0, 0, 0}, {}, {}, {}, &dict(), nullptr, o2::detectors::DetID::ITS,
                                                gsl::span<const SurfaceId>{orderedSurfaces}, plan.getSurfaceCatalog());
    BOOST_CHECK(result.ok());
    BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
    BOOST_CHECK_EQUAL(frame.getSourceIntervals(ClusterSourceId{0}).size(), 0u);
    BOOST_CHECK_EQUAL(tf.getNrof(0), 0);
  }
  {
    const auto orderedSurfaces = identitySurfaces(MFTNLayers);
    const auto catalog = makeMFTTestCatalog();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    LegacyLikeDecoder decoder{o2::detectors::DetID::MFT, true};
    TimeFrame frame;
    TimeFrameScratch tf;
    std::vector<TrackingParameters> noIterations;
    const auto plan = catalogGraph(catalogView, identitySurfaces(MFTNLayers));
    configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());
    const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0},
                                                ROFTimingConfig{40, 0, 0, 0}, {}, {}, {}, &dict(), nullptr, o2::detectors::DetID::MFT,
                                                gsl::span<const SurfaceId>{orderedSurfaces}, plan.getSurfaceCatalog());
    BOOST_CHECK(result.ok());
    BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
    BOOST_CHECK_EQUAL(tf.getNrof(0), 0);
  }
}

// Focused test: a canonical catalog configured with zero tracking iterations
// (no tracking layout ever committed) must still support normalized loading --
// loadNormalizedSource() never selects or requires any tracking-iteration
// layout.
BOOST_AUTO_TEST_CASE(ZeroIterationCatalogOnlyLoadingSucceeds)
{
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const auto catalog = makeITSTestCatalog();
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};

  TimeFrame frame;

  TimeFrameScratch tf;
  std::vector<TrackingParameters> noIterations;
  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());
  BOOST_CHECK_EQUAL(frame.getNIterations(), 0u);

  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                              clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                              gsl::span<const SurfaceId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog());
  BOOST_CHECK(result.ok());
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 1u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 0).size(), 1u);
}

BOOST_AUTO_TEST_CASE(NeverConfiguredCatalogIsRejected)
{
  // Gate 4 B2 Slice 2: TimeFrame owns no catalog of its own -- "never
  // configured" is now expressed by the caller passing an empty/default
  // SurfaceCatalogView explicitly, not by TimeFrame's own internal state.
  TimeFrame frame;
  TimeFrameScratch tf;
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                              clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                              gsl::span<const SurfaceId>{}, SurfaceCatalogView{});
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::SurfaceCatalogNotConfigured);
  BOOST_CHECK(result.source == ClusterSourceId{0});
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
}

// Gate 4 B2 Slice 2 removed the StaleCatalogAfterInvalidationIsRejected test
// that used to live here (TimeFrame graph invalidation +
// MultiSourceLoadError::SurfaceCatalogStale): loadNormalizedSource() now
// receives its SurfaceCatalogView explicitly from the caller on every call,
// with no TimeFrame-owned currency concept to invalidate -- SurfaceCatalogStale
// is retained in the enum only for wire/enum-value compatibility (see
// IOUtils.h) and is no longer producible through this API.

// TPC (or any non-ITS/MFT detector) is rejected as UnsupportedDetector before
// the catalog is ever consulted -- see UnsupportedDetectorWinsOverSharedNLayers
// below. To exercise the catalog/detId mismatch specifically, both values
// must be a supported detector (ITS or MFT); the catalog here is a 7-entry
// catalog labeled MFT (an artificial, but validation-passing, construction)
// while the load is probed with detId=ITS, so detId is supported and shares
// NLayers with the TimeFrame, yet every mapped SurfaceDescriptor.detectorId
// still disagrees with it (Gate 4 B2 Slice 2: this per-surface check, against
// the caller-supplied SurfaceCatalogView, is what DetectorSurfaceMismatch
// means now -- there is no separate catalogRequest.detector to disagree with
// any more).
BOOST_AUTO_TEST_CASE(CatalogRequestDetectorMismatchIsRejected)
{
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  auto catalog = makeITSTestCatalog();
  for (auto& surface : catalog) {
    surface.detectorId = static_cast<uint8_t>(o2::detectors::DetID::MFT);
  }
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};

  TimeFrame frame;

  TimeFrameScratch tf;
  std::vector<TrackingParameters> noIterations;
  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());

  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                              clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                              gsl::span<const SurfaceId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog());
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::DetectorSurfaceMismatch);
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
}

// A TimeFrameScratch probed with detId=TPC shares NLayers==7 with the
// scratch. The detector-identity preflight must still reject it as
// UnsupportedDetector, and it must do so before catalog ownership is ever
// inspected: no catalog is configured here at all, so a result of
// SurfaceCatalogNotConfigured (or anything catalog-related) would prove the
// unsupported-detector check ran too late or not at all.
BOOST_AUTO_TEST_CASE(UnsupportedDetectorWinsOverSharedNLayers)
{
  TimeFrame frame;
  TimeFrameScratch tf;
  LegacyLikeDecoder decoder{o2::detectors::DetID::TPC, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                              clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::TPC,
                                              gsl::span<const SurfaceId>{}, SurfaceCatalogView{});
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::UnsupportedDetector);
  BOOST_CHECK(result.source == ClusterSourceId{0});
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
}

// A supported detector still requires a caller-adopted surface plan. The
// runtime-sized scratch has no detector-specific layer-count identity of its
// own, so an unconfigured MFT load is a catalog configuration failure.
BOOST_AUTO_TEST_CASE(SupportedDetectorRequiresConfiguredSurfacePlan)
{
  TimeFrame frame;
  TimeFrameScratch tf;
  LegacyLikeDecoder decoder{o2::detectors::DetID::MFT, true};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                              clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::MFT,
                                              gsl::span<const SurfaceId>{}, SurfaceCatalogView{});
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::SurfaceCatalogNotConfigured);
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
}

BOOST_AUTO_TEST_CASE(WrongMappingCardinalityIsRejected)
{
  // Ordered surfaces shorter than NLayers: the catalog itself is still a
  // valid, self-consistent 7-surface ITS catalog, but only the first 6 are
  // designated as the per-layer mapping.
  const std::vector<SurfaceId> shortOrderedSurfaces{SurfaceId{0}, SurfaceId{1}, SurfaceId{2}, SurfaceId{3}, SurfaceId{4}, SurfaceId{5}};
  const auto catalog = makeITSTestCatalog();
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};

  TimeFrame frame;

  TimeFrameScratch tf;
  std::vector<TrackingParameters> noIterations;
  const auto plan = catalogGraph(catalogView, shortOrderedSurfaces);
  // The scratch is intentionally configured for the canonical seven-surface
  // ITS plan; the six-surface mapping below must therefore be rejected.
  configureScratchFromPlan(tf, catalog.size());

  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                              clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                              gsl::span<const SurfaceId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog());
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);
}

BOOST_AUTO_TEST_CASE(InvalidOrOutOfRangeMappedSurfaceIsRejected)
{
  {
    // Explicitly invalid SurfaceId in the mapping.
    auto orderedSurfaces = identitySurfaces(ITSNLayers);
    orderedSurfaces[3] = SurfaceId::invalid();
    const auto catalog = makeITSTestCatalog();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    TimeFrame frame;
    TimeFrameScratch tf;
    std::vector<TrackingParameters> noIterations;
    const auto plan = catalogGraph(catalogView, identitySurfaces(ITSNLayers));
    configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());

    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
    const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
    const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
    const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
    const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                                clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                                gsl::span<const SurfaceId>{orderedSurfaces}, plan.getSurfaceCatalog());
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);
  }
  {
    // Out-of-range SurfaceId: beyond the 7-entry catalog.
    auto orderedSurfaces = identitySurfaces(ITSNLayers);
    orderedSurfaces[3] = SurfaceId{100};
    const auto catalog = makeITSTestCatalog();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    TimeFrame frame;
    TimeFrameScratch tf;
    std::vector<TrackingParameters> noIterations;
    const auto plan = catalogGraph(catalogView, identitySurfaces(ITSNLayers));
    configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());

    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
    const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
    const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
    const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
    const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                                clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                                gsl::span<const SurfaceId>{orderedSurfaces}, plan.getSurfaceCatalog());
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);
  }
}

BOOST_AUTO_TEST_CASE(DuplicateMappedSurfaceIsRejected)
{
  auto orderedSurfaces = identitySurfaces(ITSNLayers);
  orderedSurfaces[3] = orderedSurfaces[2]; // duplicate: surface 2 mapped twice, surface 3 unmapped
  const auto catalog = makeITSTestCatalog();
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};

  TimeFrame frame;

  TimeFrameScratch tf;
  std::vector<TrackingParameters> noIterations;
  const auto plan = catalogGraph(catalogView, identitySurfaces(ITSNLayers));
  configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());

  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                              clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                              gsl::span<const SurfaceId>{orderedSurfaces}, plan.getSurfaceCatalog());
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);
}

// A mapped SurfaceDescriptor.detectorId that disagrees with detId for only
// ONE mapped surface out of NLayers, proving the DetectorSurfaceMismatch
// check in loadNormalizedSource() is genuinely per-surface (against the
// caller-supplied SurfaceCatalogView), not a blanket catalog-level check.
// Gate 4 B2 Slice 2 removed the runtime catalog-request validation window
// this test used to route around (layout derivation performs no
// runtime catalog validation at all -- see its own doc comment), so the
// combined MFT+ITS catalog below needs no firstSurface/count bookkeeping any
// more; it is built and used directly.
BOOST_AUTO_TEST_CASE(MappedDescriptorDetectorMismatchIsRejected)
{
  auto mftHead = makeMFTTestCatalog(); // occupies global SurfaceIds [0, MFTNLayers)
  auto itsTail = makeITSTestCatalog();
  for (uint16_t i = 0; i < itsTail.size(); ++i) {
    itsTail[i].id = SurfaceId{static_cast<uint16_t>(MFTNLayers + i)};
  }
  std::vector<SurfaceDescriptor> combinedCatalog{mftHead};
  combinedCatalog.insert(combinedCatalog.end(), itsTail.begin(), itsTail.end());
  const SurfaceCatalogView catalogView{combinedCatalog.data(), static_cast<uint32_t>(combinedCatalog.size())};

  auto orderedSurfaces = identitySurfaces(ITSNLayers);
  for (auto& surfaceId : orderedSurfaces) {
    surfaceId = SurfaceId{static_cast<uint16_t>(MFTNLayers + surfaceId.value())};
  }
  orderedSurfaces[3] = SurfaceId{2}; // an MFT surface from the unrelated prefix, wrong detector for an ITS mapping

  TimeFrame frame;

  TimeFrameScratch tf;
  std::vector<TrackingParameters> noIterations;
  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());

  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                              clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                              gsl::span<const SurfaceId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog());
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::DetectorSurfaceMismatch);
}

BOOST_AUTO_TEST_CASE(FailedNormalizedLoadLeavesBothRepresentationsUnchanged)
{
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const auto catalog = makeITSTestCatalog();
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const ROFTimingConfig timing{40, 0, 0, 0};

  TimeFrame frame;

  TimeFrameScratch tf;
  std::vector<TrackingParameters> noIterations;
  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());
  const gsl::span<const SurfaceId> planOrderedSurfaces{plan.getOrderedSurfaces()};

  // First, a valid baseline load to give the TimeFrame real content in both
  // the normalized owner and the legacy compatibility structures.
  const std::vector<CompClusterExt> goodClusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto goodPatterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> goodRofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto baseline = tf.loadNormalizedSource(frame, decoder, {0, 0},
                                                timing, goodClusters, goodPatterns, goodRofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                                planOrderedSurfaces, plan.getSurfaceCatalog());
  BOOST_REQUIRE(baseline.ok());
  BOOST_REQUIRE_EQUAL(frame.getTotalMeasurements(), 1u);
  BOOST_REQUIRE_EQUAL(tf.getUnsortedClustersOnLayer(0, 0).size(), 1u);
  BOOST_REQUIRE_EQUAL(tf.getNrof(0), 1);
  BOOST_REQUIRE(tf.getSurfaceSource(0) == ClusterSourceId{0});

  // Then a malformed ROF partition (leading gap): firstEntry != 0. This
  // fails inside loadSources() itself, after the preflight above has passed.
  const std::vector<CompClusterExt> badClusters{{2, 2, CompCluster::InvalidPatternID, 0}};
  const auto badPatterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> badRofs{ROFRecord{{0, 0}, 0, 1, 1}}; // gap: cluster 0 unreferenced

  const auto failed = tf.loadNormalizedSource(frame, decoder, {0, 0},
                                              timing, badClusters, badPatterns, badRofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                              planOrderedSurfaces, plan.getSurfaceCatalog());
  BOOST_CHECK(!failed.ok());
  BOOST_CHECK(failed.error == MultiSourceLoadError::InvalidROFRange);

  // Both the normalized owner and the legacy compatibility structures retain
  // exactly their pre-failure (baseline) content.
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 1u);
  BOOST_CHECK_EQUAL(frame.getSourceIntervals(ClusterSourceId{0}).size(), 1u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 0).size(), 1u);
  BOOST_CHECK_EQUAL(tf.getNrof(0), 1);
  BOOST_CHECK_EQUAL(tf.getClusterExternalIndex(0, 0), 0);
  BOOST_CHECK(tf.getSurfaceSource(0) == ClusterSourceId{0});
}

// Every preflight failure -- not only a loadSources()-internal one -- must
// also preserve a prior successful normalized and legacy load. Gate 4 B2
// Slice 2 removed this test's original mechanism (invalidate the
// TimeFrame-owned catalog between the baseline load and the failing one, to
// produce SurfaceCatalogStale): loadNormalizedSource() now receives its
// SurfaceCatalogView explicitly on every call, with no TimeFrame-owned
// currency to invalidate. The second call below instead passes an
// empty/unconfigured view directly, producing SurfaceCatalogNotConfigured --
// still a preflight-only failure, caught entirely in loadNormalizedSource()'s
// own preflight, before loadSources() is ever called, exactly like the
// removed mechanism.
BOOST_AUTO_TEST_CASE(PreflightFailureAfterBaselineLoadPreservesState)
{
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const auto catalog = makeITSTestCatalog();
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const ROFTimingConfig timing{40, 0, 0, 0};

  TimeFrame frame;

  TimeFrameScratch tf;
  std::vector<TrackingParameters> noIterations;
  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());

  const std::vector<CompClusterExt> goodClusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto goodPatterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> goodRofs{ROFRecord{{0, 0}, 0, 0, 1}};
  const auto baseline = tf.loadNormalizedSource(frame, decoder, {0, 0},
                                                timing, goodClusters, goodPatterns, goodRofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                                gsl::span<const SurfaceId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog());
  BOOST_REQUIRE(baseline.ok());
  BOOST_REQUIRE_EQUAL(frame.getTotalMeasurements(), 1u);
  BOOST_REQUIRE_EQUAL(tf.getUnsortedClustersOnLayer(0, 0).size(), 1u);

  const auto failed = tf.loadNormalizedSource(frame, decoder, {0, 0},
                                              timing, goodClusters, goodPatterns, goodRofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                              gsl::span<const SurfaceId>{}, SurfaceCatalogView{});
  BOOST_CHECK(!failed.ok());
  BOOST_CHECK(failed.error == MultiSourceLoadError::SurfaceCatalogNotConfigured);

  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 1u);
  BOOST_CHECK_EQUAL(tf.getUnsortedClustersOnLayer(0, 0).size(), 1u);
  BOOST_CHECK_EQUAL(tf.getNrof(0), 1);
  BOOST_CHECK_EQUAL(tf.getClusterExternalIndex(0, 0), 0);
}

// Normalized loading must match the configured covariance default and honor
// an explicit override when callers deliberately disable systematic errors.
BOOST_AUTO_TEST_CASE(ApplySysErrorsDefaultsTrueAndPropagatesToTheDecoder)
{
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
  const ROFTimingConfig timing{40, 0, 0, 0};

  {
    // No explicit argument: must match loadROFrameData()'s own default.
    const auto catalog = makeITSTestCatalog();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
    TimeFrame frame;
    TimeFrameScratch tf;
    std::vector<TrackingParameters> noIterations;
    const auto plan = catalogGraph(catalogView, orderedSurfaces);
    configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());
    const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, timing,
                                                clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                                gsl::span<const SurfaceId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog());
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK(decoder.lastApplySysErrors);
  }
  {
    // Explicit false is honored and reaches the decoder unchanged.
    const auto catalog = makeITSTestCatalog();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
    TimeFrame frame;
    TimeFrameScratch tf;
    std::vector<TrackingParameters> noIterations;
    const auto plan = catalogGraph(catalogView, orderedSurfaces);
    configureScratchFromPlan(tf, plan.getOrderedSurfaces().size());
    const auto result = tf.loadNormalizedSource(frame, decoder, {0, 0}, timing,
                                                clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                                gsl::span<const SurfaceId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog(), false);
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK(!decoder.lastApplySysErrors);
  }
}
