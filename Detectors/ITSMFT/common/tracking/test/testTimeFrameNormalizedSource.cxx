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
// through TimeFrame measurement storage (ITSMFTTracking/IOUtils.h) and then
// backfills the detector adapter's compatibility views from the loaded
// measurements.
//
// Gate 2 (Slice B3): loadNormalizedSource() no longer accepts an externally
// supplied topology view or layer-to-surface mapping. The immutable
// catalog/order is passed explicitly at the adapter boundary.
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
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{
DetectorLayout catalogGraph(SurfaceCatalogView catalog, gsl::span<const LayerId> ordered)
{
  return DetectorLayout{gsl::span<const SurfaceDescriptor>{catalog.surfaces, catalog.nSurfaces},
                        makeDetectorLayout()};
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

  o2::itsmft::tracking::ClusterDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const TopologyDictionary* dict,
    uint32_t,
    bool applySysErrors) const override
  {
    lastApplySysErrors = applySysErrors;
    const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dict);
    if (!clusterData.ok()) {
      o2::itsmft::tracking::ClusterDecodeResult result;
      result.error = clusterData.error;
      return result;
    }

    o2::itsmft::tracking::ClusterDecodeResult result;
    const int sensorID = cluster.getSensorID();
    auto& decoded = result.decoded;
    decoded.global = {static_cast<float>(sensorID) * 10.f, static_cast<float>(cluster.getRow()), static_cast<float>(cluster.getCol())};
    decoded.cylinderFrame = {static_cast<float>(sensorID) + 100.f, static_cast<float>(cluster.getRow()) + 1.f, static_cast<float>(cluster.getCol()) + 2.f, 0.01f * sensorID};
    decoded.rowColumnCovariance = {clusterData.sig2Row, 0.f, clusterData.sig2Col};
    decoded.shape = clusterData.shape;
    decoded.layer = sensorID;
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

void configureFrame(TimeFrame& frame, SurfaceCatalogView catalog, gsl::span<const LayerId> orderedSurfaces)
{
  auto layout = catalogGraph(catalog, orderedSurfaces);
  BOOST_REQUIRE(frame.configure(std::move(layout), 0, 0,
                                std::make_shared<o2::its::BoundedMemoryResource>()));
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
    surfaces.push_back(SurfaceDescriptor{i, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  }
  return surfaces;
}

// Explicit, hard-coded 10-surface MFT (all disk) canonical catalog.
std::vector<SurfaceDescriptor> makeMFTTestCatalog()
{
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.reserve(MFTNLayers);
  for (uint16_t i = 0; i < MFTNLayers; ++i) {
    surfaces.push_back(SurfaceDescriptor{i, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
  }
  return surfaces;
}

std::vector<LayerId> identitySurfaces(uint16_t nLayers)
{
  std::vector<LayerId> mapping;
  mapping.reserve(nLayers);
  for (uint16_t i = 0; i < nLayers; ++i) {
    mapping.push_back(LayerId{i});
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
    catalog.push_back({layer, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  }
  for (uint16_t layer = 0; layer < MFTNLayers; ++layer) {
    catalog.push_back({layer, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
  }
  const std::array<LayerId, ITSNLayers> itsSurfaces{LayerId{0}, LayerId{1}, LayerId{2}, LayerId{3}, LayerId{4}, LayerId{5}, LayerId{6}};
  std::array<LayerId, MFTNLayers> mftSurfaces{};
  for (uint16_t layer = 0; layer < MFTNLayers; ++layer) {
    mftSurfaces[layer] = LayerId{static_cast<uint16_t>(ITSNLayers + layer)};
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
  std::vector<LayerId> combinedSurfaces{itsSurfaces.begin(), itsSurfaces.end()};
  combinedSurfaces.insert(combinedSurfaces.end(), mftSurfaces.begin(), mftSurfaces.end());
  auto layout = catalogGraph(view, combinedSurfaces);
  BOOST_REQUIRE(frame.configure(std::move(layout), 0, 0,
                                std::make_shared<BoundedMemoryResource>()));
  const std::array<ClusterSourceInput, 2> sources{itsSource, mftSource};
  std::vector<std::vector<uint32_t>> externalIndicesBySurface;
  std::vector<std::vector<uint32_t>> clusterSizesBySurface;
  BOOST_REQUIRE(loadTimeFrameSources(frame, gsl::span<const ClusterSourceInput>{sources}, view, {50, 5},
                                     &externalIndicesBySurface, &clusterSizesBySurface)
                  .ok());
  BOOST_CHECK_EQUAL(frame.getGlobalMeasurements(LayerId{0}).size(), 2u);
  BOOST_CHECK_EQUAL(frame.getGlobalMeasurements(LayerId{static_cast<uint16_t>(ITSNLayers)}).size(), 2u);
  BOOST_CHECK_EQUAL(frame.getTotalClusters(), static_cast<int>(its.clusters.size() + mft.clusters.size()));
  BOOST_CHECK(frame.getGlobalMeasurements(LayerId{0})[0].hasValidClusterId());
  BOOST_CHECK(frame.getGlobalMeasurements(LayerId{static_cast<uint16_t>(ITSNLayers)})[0].hasValidClusterId());
  BOOST_CHECK_EQUAL(externalIndicesBySurface[0][0], 0u);
  BOOST_CHECK_EQUAL(externalIndicesBySurface[ITSNLayers][0], 0u);
  BOOST_CHECK_EQUAL(clusterSizesBySurface[0].size(), 2u);

  // Frame reset clears only normalized event data; adapter-owned index
  // translation has an independent lifecycle.
  frame.resetTimeFrame();
  BOOST_CHECK(frame.empty());
  BOOST_CHECK(!externalIndicesBySurface.empty());
  BOOST_CHECK(!clusterSizesBySurface.empty());
  externalIndicesBySurface.clear();
  clusterSizesBySurface.clear();

  // A malformed replacement is rejected before the no-throw three-owner
  // commit; the still-live MFT scratch and shared normalized owner survive.
  auto malformedMFT = mftSource;
  malformedMFT.patterns = {};
  const std::array<ClusterSourceInput, 2> retrySources{itsSource, malformedMFT};
  BOOST_CHECK(!loadTimeFrameSources(frame, gsl::span<const ClusterSourceInput>{retrySources}, view, {50, 5}).ok());
  BOOST_CHECK(frame.empty());
  frame.resetTimeFrame();
  BOOST_CHECK(frame.empty());
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
  TimeFrame frame;
  std::vector<std::vector<uint32_t>> externalIndicesBySurface;
  std::vector<std::vector<uint32_t>> clusterSizesBySurface;

  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  configureFrame(frame, catalogView, orderedSurfaces);

  const auto result = loadTimeFrameSource(frame, decoder, origin, timing,
                                          f.clusters, f.patterns, f.rofs, &dict(), &f.labels, f.detector,
                                          gsl::span<const LayerId>{identitySurfaces(static_cast<uint16_t>(plan.size()))}, plan.getSurfaceCatalog(), true,
                                          &externalIndicesBySurface, &clusterSizesBySurface);
  BOOST_REQUIRE(result.ok());

  // --- cluster counts per normalized surface (identity layout) ---
  BOOST_CHECK_EQUAL(frame.getGlobalMeasurements(LayerId{0}).size(), 2u); // clusters 0,2
  BOOST_CHECK_EQUAL(frame.getGlobalMeasurements(LayerId{1}).size(), 1u); // cluster 1
  BOOST_CHECK_EQUAL(frame.getGlobalMeasurements(LayerId{2}).size(), 1u); // cluster 3
  for (int l = 3; l < NLayers; ++l) {
    BOOST_CHECK_EQUAL(frame.getGlobalMeasurements(LayerId{static_cast<uint16_t>(l)}).size(), 0u);
    BOOST_CHECK_EQUAL(frame.getNrof(l), static_cast<int>(f.rofs.size()));
  }

  BOOST_CHECK_EQUAL(frame.getNrof(0), static_cast<int>(f.rofs.size()));

  for (const auto& e : expectedClusters) {
    // Find the matching normalized measurement.
    const GlobalMeasurement* globalMeasurement = nullptr;
    const SurfaceMeasurement* measurement = nullptr;
    uint32_t localClusterId = std::numeric_limits<uint32_t>::max();
    const auto surface = LayerId{static_cast<uint16_t>(e.layer)};
    const auto globals = frame.getGlobalMeasurements(surface);
    for (size_t index = 0; index < globals.size(); ++index) {
      if (externalIndicesBySurface[surface.value()][globals[index].clusterId] == e.externalIndex) {
        globalMeasurement = &globals[index];
        localClusterId = globals[index].clusterId;
        measurement = frame.getSurfaceMeasurement(surface, globals[index].clusterId);
        break;
      }
    }
    BOOST_REQUIRE(globalMeasurement != nullptr);
    BOOST_REQUIRE(measurement != nullptr);

    // --- global position ---
    const auto g = expectedGlobal(e.sensorID, e.row, e.col);
    BOOST_CHECK_EQUAL(globalMeasurement->position.x, g.x);
    BOOST_CHECK_EQUAL(globalMeasurement->position.y, g.y);
    BOOST_CHECK_EQUAL(globalMeasurement->position.z, g.z);

    // --- tracking-frame coordinates and covariance ---
    if (f.disk) {
      BOOST_CHECK_EQUAL(measurement->frame.q, g.z);
      BOOST_CHECK_EQUAL(measurement->frame.u, g.x);
      BOOST_CHECK_EQUAL(measurement->frame.v, g.y);
    } else {
      BOOST_CHECK_EQUAL(measurement->frame.q, static_cast<float>(e.sensorID) + 100.f);
      BOOST_CHECK_EQUAL(measurement->frame.u, static_cast<float>(e.row) + 1.f);
      BOOST_CHECK_EQUAL(measurement->frame.v, static_cast<float>(e.col) + 2.f);
      BOOST_CHECK_EQUAL(measurement->frame.frameAngle, 0.01f * e.sensorID);
    }
    // With CompCluster::InvalidPatternID, extractClusterData always returns
    // the fixed default half-pixel covariance regardless of pattern shape or
    // detector -- the same constants the legacy per-detector loaders fall
    // back to.
    BOOST_CHECK_EQUAL(measurement->covariance.uu, o2::itsmft::ioutils::DefClusError2Row);
    BOOST_CHECK_EQUAL(measurement->covariance.uv, 0.f);
    BOOST_CHECK_EQUAL(measurement->covariance.vv, o2::itsmft::ioutils::DefClusError2Col);

    // --- stable TimeFrame-local identity and adapter-owned external identity ---
    BOOST_CHECK_EQUAL(externalIndicesBySurface[surface.value()][localClusterId], e.externalIndex);

    // --- cluster size: explicit pattern consumed exactly once ---
    BOOST_CHECK_EQUAL(clusterSizesBySurface[surface.value()][localClusterId], e.nPixels);

    // --- labels ---
    const auto normalizedLabels = frame.getLabels(surface, localClusterId);
    BOOST_REQUIRE_EQUAL(normalizedLabels.size(), 1u);
    BOOST_CHECK(normalizedLabels[0] == o2::MCCompLabel(static_cast<int>(e.externalIndex) + 1, 0, 0));
  }
}

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
    const auto plan = catalogGraph(catalogView, identitySurfaces(ITSNLayers));
    configureFrame(frame, catalogView, identitySurfaces(static_cast<uint16_t>(plan.size())));
    const auto result = loadTimeFrameSource(frame, decoder, {0, 0},
                                            ROFTimingConfig{40, 0, 0, 0}, {}, {}, {}, &dict(), nullptr, o2::detectors::DetID::ITS,
                                            gsl::span<const LayerId>{orderedSurfaces}, plan.getSurfaceCatalog());
    BOOST_CHECK(result.ok());
    BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
    BOOST_CHECK_EQUAL(frame.getNrof(0), 0);
  }
  {
    const auto orderedSurfaces = identitySurfaces(MFTNLayers);
    const auto catalog = makeMFTTestCatalog();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    LegacyLikeDecoder decoder{o2::detectors::DetID::MFT, true};
    TimeFrame frame;
    const auto plan = catalogGraph(catalogView, identitySurfaces(MFTNLayers));
    configureFrame(frame, catalogView, identitySurfaces(static_cast<uint16_t>(plan.size())));
    const auto result = loadTimeFrameSource(frame, decoder, {0, 0},
                                            ROFTimingConfig{40, 0, 0, 0}, {}, {}, {}, &dict(), nullptr, o2::detectors::DetID::MFT,
                                            gsl::span<const LayerId>{orderedSurfaces}, plan.getSurfaceCatalog());
    BOOST_CHECK(result.ok());
    BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
    BOOST_CHECK_EQUAL(frame.getNrof(0), 0);
  }
}

// A configured frame owns both its catalog/layout and normalized event data.
BOOST_AUTO_TEST_CASE(ConfiguredCatalogLoadingSucceeds)
{
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const auto catalog = makeITSTestCatalog();
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};

  TimeFrame frame;

  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  configureFrame(frame, catalogView, identitySurfaces(static_cast<uint16_t>(plan.size())));

  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                          clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                          gsl::span<const LayerId>{identitySurfaces(static_cast<uint16_t>(plan.size()))}, plan.getSurfaceCatalog());
  BOOST_CHECK(result.ok());
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 1u);
  BOOST_CHECK_EQUAL(frame.getGlobalMeasurements(LayerId{0}).size(), 1u);
}

BOOST_AUTO_TEST_CASE(NeverConfiguredCatalogIsRejected)
{
  // Gate 4 B2 Slice 2: TimeFrame owns no catalog of its own -- "never
  // configured" is now expressed by the caller passing an empty/default
  // SurfaceCatalogView explicitly, not by TimeFrame's own internal state.
  TimeFrame frame;
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                          clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                          gsl::span<const LayerId>{}, SurfaceCatalogView{});
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

  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  configureFrame(frame, catalogView, identitySurfaces(static_cast<uint16_t>(plan.size())));

  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                          clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                          gsl::span<const LayerId>{identitySurfaces(static_cast<uint16_t>(plan.size()))}, plan.getSurfaceCatalog());
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
  LegacyLikeDecoder decoder{o2::detectors::DetID::TPC, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                          clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::TPC,
                                          gsl::span<const LayerId>{}, SurfaceCatalogView{});
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
  LegacyLikeDecoder decoder{o2::detectors::DetID::MFT, true};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                          clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::MFT,
                                          gsl::span<const LayerId>{}, SurfaceCatalogView{});
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::SurfaceCatalogNotConfigured);
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
}

BOOST_AUTO_TEST_CASE(WrongMappingCardinalityIsRejected)
{
  // Ordered surfaces shorter than NLayers: the catalog itself is still a
  // valid, self-consistent 7-surface ITS catalog, but only the first 6 are
  // designated as the per-layer mapping.
  const std::vector<LayerId> shortLayerMapping{LayerId{0}, LayerId{1}, LayerId{2}, LayerId{3}, LayerId{4}, LayerId{5}};
  const auto catalog = makeITSTestCatalog();
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};

  TimeFrame frame;

  const auto plan = catalogGraph(catalogView, shortLayerMapping);
  // The scratch is intentionally configured for the canonical seven-surface
  // ITS plan; the six-surface mapping below must therefore be rejected.
  configureFrame(frame, catalogView, identitySurfaces(ITSNLayers));

  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                          clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                          gsl::span<const LayerId>{shortLayerMapping}, plan.getSurfaceCatalog());
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);
}

BOOST_AUTO_TEST_CASE(InvalidOrOutOfRangeMappedSurfaceIsRejected)
{
  {
    // Explicitly invalid LayerId in the mapping.
    auto orderedSurfaces = identitySurfaces(ITSNLayers);
    orderedSurfaces[3] = LayerId::invalid();
    const auto catalog = makeITSTestCatalog();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    TimeFrame frame;
    const auto plan = catalogGraph(catalogView, identitySurfaces(ITSNLayers));
    configureFrame(frame, catalogView, identitySurfaces(static_cast<uint16_t>(plan.size())));

    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
    const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
    const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
    const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
    const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                            clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                            gsl::span<const LayerId>{orderedSurfaces}, plan.getSurfaceCatalog());
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);
  }
  {
    // Out-of-range LayerId: beyond the 7-entry catalog.
    auto orderedSurfaces = identitySurfaces(ITSNLayers);
    orderedSurfaces[3] = LayerId{100};
    const auto catalog = makeITSTestCatalog();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    TimeFrame frame;
    const auto plan = catalogGraph(catalogView, identitySurfaces(ITSNLayers));
    configureFrame(frame, catalogView, identitySurfaces(static_cast<uint16_t>(plan.size())));

    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
    const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
    const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
    const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
    const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                            clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                            gsl::span<const LayerId>{orderedSurfaces}, plan.getSurfaceCatalog());
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

  const auto plan = catalogGraph(catalogView, identitySurfaces(ITSNLayers));
  configureFrame(frame, catalogView, identitySurfaces(static_cast<uint16_t>(plan.size())));

  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                          clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                          gsl::span<const LayerId>{orderedSurfaces}, plan.getSurfaceCatalog());
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
  auto mftHead = makeMFTTestCatalog(); // occupies global LayerIds [0, MFTNLayers)
  auto itsTail = makeITSTestCatalog();
  for (uint16_t i = 0; i < itsTail.size(); ++i) {
  }
  std::vector<SurfaceDescriptor> combinedCatalog{mftHead};
  combinedCatalog.insert(combinedCatalog.end(), itsTail.begin(), itsTail.end());
  const SurfaceCatalogView catalogView{combinedCatalog.data(), static_cast<uint32_t>(combinedCatalog.size())};

  auto orderedSurfaces = identitySurfaces(ITSNLayers);
  for (auto& layerId : orderedSurfaces) {
    layerId = LayerId{static_cast<uint16_t>(MFTNLayers + layerId.value())};
  }
  orderedSurfaces[3] = LayerId{2}; // an MFT surface from the unrelated prefix, wrong detector for an ITS mapping

  TimeFrame frame;

  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  configureFrame(frame, catalogView, identitySurfaces(static_cast<uint16_t>(plan.size())));

  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, ROFTimingConfig{40, 0, 0, 0},
                                          clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                          gsl::span<const LayerId>{identitySurfaces(static_cast<uint16_t>(plan.size()))}, plan.getSurfaceCatalog());
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::DetectorSurfaceMismatch);
}

BOOST_AUTO_TEST_CASE(FailedNormalizedLoadClearsTheTimeFrame)
{
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const auto catalog = makeITSTestCatalog();
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const ROFTimingConfig timing{40, 0, 0, 0};

  TimeFrame frame;

  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  const auto planLayerMapping = identitySurfaces(static_cast<uint16_t>(plan.size()));
  configureFrame(frame, catalogView, planLayerMapping);
  const gsl::span<const LayerId> planMapping{planLayerMapping};

  // First, a valid baseline load to give the TimeFrame real content in both
  // the normalized owner and the legacy compatibility structures.
  const std::vector<CompClusterExt> goodClusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto goodPatterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> goodRofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const auto baseline = loadTimeFrameSource(frame, decoder, {0, 0},
                                            timing, goodClusters, goodPatterns, goodRofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                            planMapping, plan.getSurfaceCatalog());
  BOOST_REQUIRE(baseline.ok());
  BOOST_REQUIRE_EQUAL(frame.getTotalMeasurements(), 1u);
  BOOST_REQUIRE_EQUAL(frame.getGlobalMeasurements(LayerId{0}).size(), 1u);
  BOOST_REQUIRE_EQUAL(frame.getNrof(0), 1);

  // Then a malformed ROF partition (leading gap): firstEntry != 0. This
  // fails inside loadSources() itself, after the preflight above has passed.
  const std::vector<CompClusterExt> badClusters{{2, 2, CompCluster::InvalidPatternID, 0}};
  const auto badPatterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> badRofs{ROFRecord{{0, 0}, 0, 1, 1}}; // gap: cluster 0 unreferenced

  const auto failed = loadTimeFrameSource(frame, decoder, {0, 0},
                                          timing, badClusters, badPatterns, badRofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                          planMapping, plan.getSurfaceCatalog());
  BOOST_CHECK(!failed.ok());
  BOOST_CHECK(failed.error == MultiSourceLoadError::InvalidROFRange);

  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
  BOOST_CHECK(frame.getGlobalMeasurements(LayerId{0}).empty());
}

// Every preflight failure -- not only a loadSources()-internal one -- must
// A preflight failure must clear a prior successful load as well as failures
// encountered after decoding begins.
BOOST_AUTO_TEST_CASE(PreflightFailureAfterBaselineLoadClearsState)
{
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const auto catalog = makeITSTestCatalog();
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
  const ROFTimingConfig timing{40, 0, 0, 0};

  TimeFrame frame;

  const auto plan = catalogGraph(catalogView, orderedSurfaces);
  configureFrame(frame, catalogView, identitySurfaces(static_cast<uint16_t>(plan.size())));

  const std::vector<CompClusterExt> goodClusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto goodPatterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> goodRofs{ROFRecord{{0, 0}, 0, 0, 1}};
  const auto baseline = loadTimeFrameSource(frame, decoder, {0, 0},
                                            timing, goodClusters, goodPatterns, goodRofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                            gsl::span<const LayerId>{identitySurfaces(static_cast<uint16_t>(plan.size()))}, plan.getSurfaceCatalog());
  BOOST_REQUIRE(baseline.ok());
  BOOST_REQUIRE_EQUAL(frame.getTotalMeasurements(), 1u);
  BOOST_REQUIRE_EQUAL(frame.getGlobalMeasurements(LayerId{0}).size(), 1u);

  const auto failed = loadTimeFrameSource(frame, decoder, {0, 0},
                                          timing, goodClusters, goodPatterns, goodRofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                          gsl::span<const LayerId>{}, SurfaceCatalogView{});
  BOOST_CHECK(!failed.ok());
  BOOST_CHECK(failed.error == MultiSourceLoadError::SurfaceCatalogNotConfigured);

  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
  BOOST_CHECK(frame.getGlobalMeasurements(LayerId{0}).empty());
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
    const auto plan = catalogGraph(catalogView, orderedSurfaces);
    configureFrame(frame, catalogView, identitySurfaces(static_cast<uint16_t>(plan.size())));
    const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, timing,
                                            clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                            gsl::span<const LayerId>{identitySurfaces(static_cast<uint16_t>(plan.size()))}, plan.getSurfaceCatalog());
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK(decoder.lastApplySysErrors);
  }
  {
    // Explicit false is honored and reaches the decoder unchanged.
    const auto catalog = makeITSTestCatalog();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS, false};
    TimeFrame frame;
    const auto plan = catalogGraph(catalogView, orderedSurfaces);
    configureFrame(frame, catalogView, identitySurfaces(static_cast<uint16_t>(plan.size())));
    const auto result = loadTimeFrameSource(frame, decoder, {0, 0}, timing,
                                            clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                            gsl::span<const LayerId>{identitySurfaces(static_cast<uint16_t>(plan.size()))}, plan.getSurfaceCatalog(), false);
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK(!decoder.lastApplySysErrors);
  }
}
