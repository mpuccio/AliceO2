// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 1 compatibility-bridge parity tests (AgentCoordination.md / Wave 1 exit
// criteria: "Current single-detector TimeFrames still load and expose
// equivalent clusters"). These exercise o2::itsmft::tracking::bridge::
// loadLegacySource() -- the facade that wraps one legacy single-detector
// cluster stream (the same raw ROFRecord/CompClusterExt/pattern/dictionary/
// label inputs that o2::its::TimeFrame<NLayers>::loadROFrameData and
// o2::mft::ioutils::loadROFrameData already consume) into loadSources().
//
// Like every other Gate 1 test in this directory (see
// testMultiSourceLoading.cxx), these tests use a host-only decoder instead of
// the production GeometryClusterDecoder<DetId>, because no unit test in this
// component depends on a loaded ITS/MFT TGeo geometry singleton (that is only
// exercised by the heavier Gate 0 fixture/replay protocol under
// test/gate0-baseline). The decoder still reuses the real
// o2::itsmft::ioutils::extractClusterData pattern/dictionary consumption
// path, so pattern-cursor bookkeeping is exercised exactly as production
// would; only the geometry transform itself is a deterministic stand-in.
// Real-geometry decode parity for both ITS and MFT (T2L/L2G matrices,
// covariance ordering) was already established by the accepted D005
// SurfaceMeasurementAdapters work and is not re-derived here.

#define BOOST_TEST_MODULE ITSMFT LegacyTimeFrameBridge
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
#include "ITSMFTTracking/LegacyTimeFrameBridge.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
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
// exercised identically to GeometryClusterDecoder.
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
    bool /*applySysErrors*/) const override
  {
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

 private:
  o2::detectors::DetID::ID mDetector;
  bool mDisk;
};

const TopologyDictionary& dict()
{
  static const TopologyDictionary d;
  return d;
}

// Explicit (non-grouped, InvalidPatternID) patterns: header {rowSpan,colSpan}
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

// Expected values, computed independently of LegacyLikeDecoder's formulas by
// re-deriving them from the same (sensorID,row,col) inputs, so the test does
// not just echo the decoder's own arithmetic back at itself for the
// higher-level bridge-wiring checks (identity/timing/labels/counts) below.
GlobalPoint3F expectedGlobal(int sensorID, int row, int col)
{
  return {static_cast<float>(sensorID) * 10.f, static_cast<float>(row), static_cast<float>(col)};
}

struct Fixture {
  o2::detectors::DetID::ID detector;
  uint16_t nLayers;
  bool disk;
  std::vector<CompClusterExt> clusters;
  std::vector<unsigned char> patterns;
  std::vector<ROFRecord> rofs;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
};

// 4 clusters on layers {0,1,0,2}, partitioned into 3 ROFs: ROF0={c0,c1},
// ROF1={c2}, ROF2={c3}. Row/col distinguish every cluster; patterns
// alternate 1x1/1x3 shapes; one label per cluster.
Fixture makeFixture(o2::detectors::DetID::ID detector, uint16_t nLayers, bool disk)
{
  Fixture f;
  f.detector = detector;
  f.nLayers = nLayers;
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

void checkParity(const Fixture& f)
{
  const auto layout = bridge::makeSingleDetectorLayout(f.detector, f.nLayers);
  BOOST_REQUIRE(layout.valid());
  const auto layerToSurface = bridge::identityLayerToSurface(f.nLayers);

  LegacyLikeDecoder decoder{f.detector, f.disk};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0};

  MultiSourceFrame frame;
  const auto result = bridge::loadLegacySource(frame, layout.getView(), f.detector, layerToSurface, ClusterSourceId{0},
                                               f.rofs, f.clusters, f.patterns, &dict(), &f.labels, decoder, origin, timing);
  BOOST_REQUIRE(result.ok());

  // --- cluster counts per legacy layer == normalized surface (identity layout) ---
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0}).size(), 2u); // clusters 0,2
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{1}).size(), 1u); // cluster 1
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{2}).size(), 1u); // cluster 3
  for (uint16_t l = 3; l < f.nLayers; ++l) {
    BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{l}).size(), 0u);
  }
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), f.clusters.size());

  const auto findOn = [&](SurfaceId s, uint32_t externalIndex) -> const SurfaceMeasurement& {
    for (const auto& m : frame.getSurfaceMeasurements(s)) {
      if (m.cluster.index == externalIndex) {
        return m;
      }
    }
    BOOST_FAIL("measurement not found");
    static SurfaceMeasurement dummy{};
    return dummy;
  };

  struct Expected {
    uint32_t externalIndex;
    SurfaceId surface;
    int sensorID;
    int row, col;
    uint32_t sourceROF;
    uint32_t nPixels, rowSpan, columnSpan;
  };
  const std::vector<Expected> expected{
    {0, SurfaceId{0}, 0, 10, 20, 0, 1, 1, 1},
    {1, SurfaceId{1}, 1, 11, 21, 0, 3, 1, 3},
    {2, SurfaceId{0}, 0, 12, 22, 1, 1, 1, 1},
    {3, SurfaceId{2}, 2, 13, 23, 2, 3, 1, 3},
  };

  for (const auto& e : expected) {
    const auto& m = findOn(e.surface, e.externalIndex);

    // --- global position ---
    const auto g = expectedGlobal(e.sensorID, e.row, e.col);
    BOOST_CHECK_EQUAL(m.global.x, g.x);
    BOOST_CHECK_EQUAL(m.global.y, g.y);
    BOOST_CHECK_EQUAL(m.global.z, g.z);

    // --- tracking/surface-frame coordinates and covariance semantics ---
    if (f.disk) {
      // makeDiskSurfaceMeasurement's accepted convention (D005): disk frame is
      // {global.z, global.x, global.y, 0}, never the legacy synthetic
      // TrackingFrameInfo projection.
      BOOST_CHECK_EQUAL(m.frame.q, g.z);
      BOOST_CHECK_EQUAL(m.frame.u, g.x);
      BOOST_CHECK_EQUAL(m.frame.v, g.y);
      BOOST_CHECK_EQUAL(m.frame.frameAngle, 0.f);
    } else {
      BOOST_CHECK_EQUAL(m.frame.q, static_cast<float>(e.sensorID) + 100.f);
      BOOST_CHECK_EQUAL(m.frame.u, static_cast<float>(e.row) + 1.f);
      BOOST_CHECK_EQUAL(m.frame.v, static_cast<float>(e.col) + 2.f);
      BOOST_CHECK_EQUAL(m.frame.frameAngle, 0.01f * e.sensorID);
    }
    // With CompCluster::InvalidPatternID, extractClusterData always returns
    // the fixed default half-pixel covariance regardless of pattern shape or
    // detector -- the same constants o2::its::TimeFrame::loadROFrameData and
    // o2::mft::ioutils::loadROFrameData fall back to.
    BOOST_CHECK_EQUAL(m.covariance.uu, o2::itsmft::ioutils::DefClusError2Row);
    BOOST_CHECK_EQUAL(m.covariance.uv, 0.f);
    BOOST_CHECK_EQUAL(m.covariance.vv, o2::itsmft::ioutils::DefClusError2Col);

    // --- sensor, source, external-index, surface and source-ROF metadata ---
    BOOST_CHECK_EQUAL(m.sensor.detector, static_cast<uint32_t>(f.detector));
    BOOST_CHECK_EQUAL(m.sensor.sensor, static_cast<uint32_t>(e.sensorID));
    BOOST_CHECK(m.cluster.source == ClusterSourceId{0});
    BOOST_CHECK_EQUAL(m.cluster.index, e.externalIndex);
    BOOST_CHECK(m.surface == e.surface);
    BOOST_CHECK_EQUAL(m.sourceROF, e.sourceROF);

    // --- cluster shape: explicit pattern consumed exactly once, in order ---
    BOOST_CHECK_EQUAL(m.shape.nPixels, e.nPixels);
    BOOST_CHECK_EQUAL(m.shape.rowSpan, e.rowSpan);
    BOOST_CHECK_EQUAL(m.shape.columnSpan, e.columnSpan);

    // --- labels resolved through ClusterRef ---
    const auto labelSpan = frame.getLabels(ClusterRef{ClusterSourceId{0}, e.externalIndex});
    BOOST_REQUIRE_EQUAL(labelSpan.size(), 1u);
    BOOST_CHECK(labelSpan[0] == o2::MCCompLabel(static_cast<int>(e.externalIndex) + 1, 0, 0));
  }

  // --- ROF / timing association ---
  BOOST_REQUIRE_EQUAL(frame.getSources().size(), 1u);
  BOOST_CHECK(frame.getSources()[0].detector == f.detector);
  BOOST_CHECK_EQUAL(frame.getSources()[0].nROFs, f.rofs.size());
  const auto intervals = frame.getSourceIntervals(ClusterSourceId{0});
  BOOST_REQUIRE_EQUAL(intervals.size(), f.rofs.size());
  for (uint32_t r = 0; r < f.rofs.size(); ++r) {
    BOOST_CHECK_EQUAL(intervals[r].sourceROF, r);
    BOOST_CHECK_EQUAL(intervals[r].begin, f.rofs[r].getBCData().differenceInBC(origin));
    BOOST_CHECK_EQUAL(intervals[r].length(), static_cast<uint64_t>(timing.rofLength));
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(ITSOnlyLegacyLoadingBridgesToEquivalentSurfaceMeasurements)
{
  checkParity(makeFixture(o2::detectors::DetID::ITS, ITSNLayers, false));
}

BOOST_AUTO_TEST_CASE(MFTOnlyLegacyLoadingBridgesToEquivalentSurfaceMeasurements)
{
  checkParity(makeFixture(o2::detectors::DetID::MFT, MFTNLayers, true));
}

// The current MFT production loader (o2::mft::ioutils::loadROFrameData)
// keeps only the first MC label per cluster
// (`*(mcLabels->getLabels(...).begin())`); the normalized owner intentionally
// exposes the full label span instead. This is a deliberate, documented
// asymmetry (Wave 1 goal: richer normalized exposure), not a regression,
// since the legacy MFT loader itself is untouched.
BOOST_AUTO_TEST_CASE(MFTBridgeExposesFullLabelSpanUnlikeLegacyFirstLabelOnly)
{
  const auto layout = bridge::makeSingleDetectorLayout(o2::detectors::DetID::MFT, MFTNLayers);
  BOOST_REQUIRE(layout.valid());
  const auto layerToSurface = bridge::identityLayerToSurface(MFTNLayers);

  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
  labels.addElement(0, o2::MCCompLabel{1, 0, 0});
  labels.addElement(0, o2::MCCompLabel{2, 0, 0});

  LegacyLikeDecoder decoder{o2::detectors::DetID::MFT, true};
  MultiSourceFrame frame;
  const auto result = bridge::loadLegacySource(frame, layout.getView(), o2::detectors::DetID::MFT, layerToSurface, ClusterSourceId{0},
                                               rofs, clusters, patterns, &dict(), &labels, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0});
  BOOST_REQUIRE(result.ok());

  const auto labelSpan = frame.getLabels(ClusterRef{ClusterSourceId{0}, 0});
  BOOST_REQUIRE_EQUAL(labelSpan.size(), 2u);
  BOOST_CHECK(labelSpan[0] == o2::MCCompLabel(1, 0, 0));
  BOOST_CHECK(labelSpan[1] == o2::MCCompLabel(2, 0, 0));
}

BOOST_AUTO_TEST_CASE(EmptyInputsAreLegalForBothDetectors)
{
  for (const auto detector : {o2::detectors::DetID::ITS, o2::detectors::DetID::MFT}) {
    const uint16_t nLayers = (detector == o2::detectors::DetID::MFT) ? MFTNLayers : ITSNLayers;
    const auto layout = bridge::makeSingleDetectorLayout(detector, nLayers);
    BOOST_REQUIRE(layout.valid());
    const auto layerToSurface = bridge::identityLayerToSurface(nLayers);

    const std::vector<CompClusterExt> clusters{};
    const std::vector<unsigned char> patterns{};
    const std::vector<ROFRecord> rofs{};
    LegacyLikeDecoder decoder{detector, detector == o2::detectors::DetID::MFT};

    MultiSourceFrame frame;
    const auto result = bridge::loadLegacySource(frame, layout.getView(), detector, layerToSurface, ClusterSourceId{0},
                                                 rofs, clusters, patterns, &dict(), nullptr, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0});
    BOOST_CHECK(result.ok());
    BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
    BOOST_REQUIRE_EQUAL(frame.getSources().size(), 1u);
    BOOST_CHECK_EQUAL(frame.getSources()[0].nROFs, 0u);
  }
}

BOOST_AUTO_TEST_CASE(FailedLoadLeavesFrameUnmodifiedTransactionally)
{
  const auto layout = bridge::makeSingleDetectorLayout(o2::detectors::DetID::ITS, ITSNLayers);
  BOOST_REQUIRE(layout.valid());
  const auto layerToSurface = bridge::identityLayerToSurface(ITSNLayers);
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};

  // First, a valid baseline load to give the frame real content.
  const std::vector<CompClusterExt> goodClusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto goodPatterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> goodRofs{ROFRecord{{0, 0}, 0, 0, 1}};

  MultiSourceFrame frame;
  const auto baseline = bridge::loadLegacySource(frame, layout.getView(), o2::detectors::DetID::ITS, layerToSurface, ClusterSourceId{0},
                                                 goodRofs, goodClusters, goodPatterns, &dict(), nullptr, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0});
  BOOST_REQUIRE(baseline.ok());
  BOOST_REQUIRE_EQUAL(frame.getTotalMeasurements(), 1u);

  // Then a malformed ROF partition (leading gap): firstEntry != 0.
  const std::vector<CompClusterExt> badClusters{{2, 2, CompCluster::InvalidPatternID, 0}};
  const auto badPatterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> badRofs{ROFRecord{{0, 0}, 0, 1, 1}}; // gap: cluster 0 unreferenced

  const auto failed = bridge::loadLegacySource(frame, layout.getView(), o2::detectors::DetID::ITS, layerToSurface, ClusterSourceId{0},
                                               badRofs, badClusters, badPatterns, &dict(), nullptr, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0});
  BOOST_CHECK(!failed.ok());
  BOOST_CHECK(failed.error == MultiSourceLoadError::InvalidROFRange);

  // Frame retains exactly its pre-failure (baseline) content.
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 1u);
  BOOST_REQUIRE_EQUAL(frame.getSources().size(), 1u);
  BOOST_CHECK_EQUAL(frame.getSources()[0].nROFs, 1u);
}
