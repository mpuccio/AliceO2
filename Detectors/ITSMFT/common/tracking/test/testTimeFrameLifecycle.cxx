// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// TimeFrame lifecycle and transactional configuration/load publication.
//
// A. Reset lifecycle: TimeFrame::resetTimeFrame() must unconditionally clear the
//    normalized owner (mNormalizedFrame) associated by
//    TimeFrameScratch::loadNormalizedSource() -- every
//    normalized accessor obtained *after* resetTimeFrame() must report empty/zero
//    content. This test never dereferences a view obtained before resetTimeFrame():
//    every post-wipe check re-obtains its accessor. (Gate 4 B3.1: neither
//    TimeFrame nor TimeFrameScratch stores mDetId at all any
//    more -- callers pass the detector explicitly to every call that needs
//    it -- so there is nothing detector-identity-shaped left for resetTimeFrame() to
//    preserve or clear.)
//
// B. Strong configuration transactionality: a BoundedMemoryResource failure
//    while staging a valid replacement must preserve the live configuration,
//    workspace, allocator and capacities, as well as an already loaded event,
//    its allocator-backed storage, navigation, and results.
//
// C. Event loading decodes into heap-backed staging and publishes through the
//    no-throw TimeFrame::commitLoadedEvent(). Malformed replacement input is
//    therefore validation coverage: rejection before publication must leave
//    the loaded event unchanged.

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

  o2::itsmft::tracking::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const TopologyDictionary* dict,
    gsl::span<const LayerId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool applySysErrors) const override
  {
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
    result = makeCylinderMeasurementDecodeResult(decoded, sensor, surface, clusterRef, sourceROF);
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
    surfaces.push_back(SurfaceDescriptor{LayerId{i}, i, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
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

SurfaceLayout catalogLayout(SurfaceCatalogView catalog, gsl::span<const LayerId> ordered)
{
  return SurfaceLayout{gsl::span<const SurfaceDescriptor>{catalog.surfaces, catalog.nSurfaces},
                       makeSurfaceLayoutChain(ordered)};
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

void verifyFixtureLoaded(const TimeFrame& frame, const Fixture& f)
{
  constexpr ClusterSourceId kSourceId{0};
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(LayerId{0}).size(), 2u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(LayerId{1}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(LayerId{2}).size(), 1u);
  for (int l = 3; l < ITSNLayers; ++l) {
    BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(LayerId{static_cast<uint16_t>(l)}).size(), 0u);
    BOOST_CHECK_EQUAL(frame.getNrof(l), static_cast<int>(f.rofs.size()));
  }

  BOOST_CHECK_EQUAL(frame.getNrof(0), static_cast<int>(f.rofs.size()));

  for (const auto& e : expectedClusters) {
    const GlobalMeasurement* globalMeasurement = nullptr;
    const SurfaceMeasurement* measurement = nullptr;
    const auto surface = LayerId{static_cast<uint16_t>(e.layer)};
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

    const auto g = expectedGlobal(e.sensorID, e.row, e.col);
    BOOST_CHECK_EQUAL(globalMeasurement->position.x, g.x);
    BOOST_CHECK_EQUAL(globalMeasurement->position.y, g.y);
    BOOST_CHECK_EQUAL(globalMeasurement->position.z, g.z);

    BOOST_CHECK_EQUAL(measurement->frame.q, static_cast<float>(e.sensorID) + 100.f);
    BOOST_CHECK_EQUAL(measurement->frame.u, static_cast<float>(e.row) + 1.f);
    BOOST_CHECK_EQUAL(measurement->frame.v, static_cast<float>(e.col) + 2.f);
    BOOST_CHECK_EQUAL(measurement->frame.frameAngle, 0.01f * e.sensorID);

    BOOST_CHECK_EQUAL(measurement->covariance.uu, o2::itsmft::ioutils::DefClusError2Row);
    BOOST_CHECK_EQUAL(measurement->covariance.uv, 0.f);
    BOOST_CHECK_EQUAL(measurement->covariance.vv, o2::itsmft::ioutils::DefClusError2Col);

    BOOST_CHECK_EQUAL(globalMeasurement->cluster.index, e.externalIndex);
    BOOST_CHECK(globalMeasurement->cluster.source == kSourceId);
    BOOST_CHECK(globalMeasurement->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::ITS));
    BOOST_CHECK_EQUAL(globalMeasurement->sensor.sensor, static_cast<uint32_t>(e.sensorID));
    BOOST_CHECK(globalMeasurement->surface == LayerId{static_cast<uint16_t>(e.layer)});
    BOOST_CHECK_EQUAL(globalMeasurement->sourceROF, e.sourceROF);

    BOOST_CHECK_EQUAL(globalMeasurement->shape.nPixels, e.nPixels);

    const auto normalizedLabels = frame.getLabels(ClusterRef{kSourceId, e.externalIndex});
    BOOST_REQUIRE_EQUAL(normalizedLabels.size(), 1u);
    BOOST_CHECK(normalizedLabels[0] == o2::MCCompLabel(static_cast<int>(e.externalIndex) + 1, 0, 0));
  }
}

void configureFrame(TimeFrame& frame, SurfaceCatalogView catalog, gsl::span<const LayerId> orderedSurfaces,
                    std::shared_ptr<BoundedMemoryResource> pool = std::make_shared<BoundedMemoryResource>())
{
  std::vector<SurfaceLayout> layouts;
  layouts.push_back(catalogLayout(catalog, orderedSurfaces));
  std::vector<TrackingParameters> parameters(1);
  std::vector<TrackingWorkspaceCapacity> capacities{{orderedSurfaces.size(), 0, 0}};
  BOOST_REQUIRE(frame.commitConfiguration(std::move(layouts), std::move(parameters), std::move(capacities), std::move(pool)));
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
  const auto plan = catalogLayout(catalogView, orderedSurfaces);
  configureFrame(frame, catalogView, orderedSurfaces);

  const auto f = makeFixture();
  const auto result = loadTimeFrameSource(frame, decoder, origin, timing, f.clusters, f.patterns, f.rofs, &dict(), &f.labels, o2::detectors::DetID::ITS,
                                          gsl::span<const LayerId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog());
  BOOST_REQUIRE(result.ok());
  // Sanity: the successful load itself has the expected content, matching
  // the accepted parity coverage in testTimeFrameNormalizedSource.cxx.
  verifyFixtureLoaded(frame, f);

  frame.resetTimeFrame();

  // --- inspect only freshly obtained normalized accessors/views ---
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
  BOOST_CHECK_EQUAL(frame.getNMeasurementSurfaces(), 0u);
  for (uint16_t s = 0; s < ITSNLayers; ++s) {
    BOOST_CHECK(frame.getSurfaceMeasurements(LayerId{s}).empty());
  }
  BOOST_CHECK(frame.getLabels(ClusterRef{ClusterSourceId{0}, 0}).empty());

  // Gate 4 B3.1: neither owner stores mDetId any more -- the plan lives on
  // `plan` above, entirely outside both TimeFrame and LegacyTrackerScratch,
  // so resetTimeFrame() has no detector-identity state to preserve or clear.
}

BOOST_AUTO_TEST_CASE(FailedConfigurationAllocationPreservesLoadedFrame)
{
  const auto catalog = makeITSTestCatalog();
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0};
  const auto fixture = makeFixture();
  const auto plan = catalogLayout(catalogView, orderedSurfaces);
  auto livePool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  configureFrame(frame, catalogView, orderedSurfaces, livePool);

  const auto loaded = loadTimeFrameSource(frame, decoder, origin, timing, fixture.clusters,
                                          fixture.patterns, fixture.rofs, &dict(),
                                          &fixture.labels, o2::detectors::DetID::ITS,
                                          gsl::span<const LayerId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog());
  BOOST_REQUIRE(loaded.ok());

  frame.getWorkspace().getPositionResolutions().resize(ITSNLayers);
  frame.getWorkspace().getPositionResolutions()[0] = 17.f;
  BOOST_REQUIRE_EQUAL(frame.getClusters().size(), ITSNLayers);
  for (uint16_t layer = 0; layer < ITSNLayers; ++layer) {
    const auto globals = frame.getGlobalMeasurements(LayerId{layer});
    for (std::size_t index = 0; index < globals.size(); ++index) {
      const auto& position = globals[index].position;
      frame.getClusters()[layer].emplace_back(position.x, position.y, position.z, static_cast<int>(index));
    }
  }
  BOOST_REQUIRE(!frame.getClusters()[0].empty());

  GenericTrack result{};
  result.chi2 = 42.f;
  result.firstClusterRef = 0;
  result.clusterRefEnd = 1;
  frame.getGenericTracks().push_back(result);
  frame.getTrackClusterIndices().push_back(TrackClusterReference{LayerId{0}, MeasurementIndex{0}});
  frame.addPrimaryVertex(Vertex{});

  verifyFixtureLoaded(frame, fixture);
  const auto* const liveLayout = &frame.getLayout(0);
  const auto* const liveParameters = frame.getTrackingParameters(0);
  const auto liveCapacity = *frame.getWorkspaceCapacity(0);
  const auto* const liveWorkspace = &frame.getWorkspace();
  const auto* const liveWorkspacePool = frame.getWorkspace().getMemoryPool().get();
  const auto liveWorkspaceSurfaces = frame.getWorkspace().getNOwnedSurfaces();
  const auto liveWorkspaceEdges = frame.getWorkspace().getNEdges();
  const auto liveWorkspaceCells = frame.getWorkspace().getNCells();
  const auto liveTraversalWorkspaces = frame.getWorkspace().getNTraversalWorkspaces();
  const auto livePositionResolution = frame.getWorkspace().getPositionResolution(0);
  const auto livePoolUsed = livePool->getUsedMemory();
  const auto livePoolMax = livePool->getMaxMemory();
  const auto liveResetCount = frame.getEventResetCount();
  const auto* const liveGlobalMeasurements = frame.getGlobalMeasurements(LayerId{0}).data();
  const auto* const liveSurfaceMeasurements = frame.getSurfaceMeasurements(LayerId{0}).data();
  const auto* const liveClusters = frame.getClusters()[0].data();
  const auto liveClusterCount = frame.getClusters()[0].size();
  const auto liveClusterId = frame.getClusters()[0][0].clusterId;
  const auto liveUsedClusters = frame.getUsedClusters(0);
  const auto* const liveUsedClusterData = liveUsedClusters.data();
  const auto liveROFBoundaries = frame.getROFrameClusters(0);
  const std::vector<int> expectedROFBoundaries{liveROFBoundaries.begin(), liveROFBoundaries.end()};
  const auto* const liveROFBoundaryData = liveROFBoundaries.data();
  const auto* const liveGenericTracks = frame.getGenericTracks().data();
  const auto* const liveTrackClusterIndices = frame.getTrackClusterIndices().data();
  const auto* const livePrimaryVertices = frame.getPrimaryVertices().data();

  std::vector<SurfaceLayout> replacementLayouts;
  replacementLayouts.push_back(catalogLayout(catalogView, orderedSurfaces));
  std::vector<TrackingParameters> replacementParameters(1);
  std::vector<TrackingWorkspaceCapacity> replacementCapacities{{ITSNLayers, 1, 1}};
  auto failingPool = std::make_shared<BoundedMemoryResource>(0);

  BOOST_CHECK(!frame.commitConfiguration(std::move(replacementLayouts), std::move(replacementParameters),
                                         std::move(replacementCapacities), failingPool));
  BOOST_CHECK_EQUAL(failingPool->getThrowCount(), 1u);
  BOOST_CHECK_EQUAL(failingPool->getUsedMemory(), 0u);

  BOOST_CHECK(frame.isConfigured());
  BOOST_CHECK_EQUAL(frame.getNIterations(), 1u);
  BOOST_CHECK(&frame.getLayout(0) == liveLayout);
  BOOST_CHECK(frame.getTrackingParameters(0) == liveParameters);
  BOOST_REQUIRE(frame.getWorkspaceCapacity(0) != nullptr);
  BOOST_CHECK_EQUAL(frame.getWorkspaceCapacity(0)->ownedSurfaces, liveCapacity.ownedSurfaces);
  BOOST_CHECK_EQUAL(frame.getWorkspaceCapacity(0)->edges, liveCapacity.edges);
  BOOST_CHECK_EQUAL(frame.getWorkspaceCapacity(0)->cells, liveCapacity.cells);
  BOOST_CHECK(&frame.getWorkspace() == liveWorkspace);
  BOOST_CHECK(frame.getMemoryPool().get() == livePool.get());
  BOOST_CHECK(frame.getWorkspace().getMemoryPool().get() == liveWorkspacePool);
  BOOST_CHECK_EQUAL(frame.getWorkspace().getNOwnedSurfaces(), liveWorkspaceSurfaces);
  BOOST_CHECK_EQUAL(frame.getWorkspace().getNEdges(), liveWorkspaceEdges);
  BOOST_CHECK_EQUAL(frame.getWorkspace().getNCells(), liveWorkspaceCells);
  BOOST_CHECK_EQUAL(frame.getWorkspace().getNTraversalWorkspaces(), liveTraversalWorkspaces);
  BOOST_CHECK_EQUAL(frame.getWorkspace().getPositionResolution(0), livePositionResolution);
  BOOST_CHECK_EQUAL(livePool->getUsedMemory(), livePoolUsed);
  BOOST_CHECK_EQUAL(livePool->getMaxMemory(), livePoolMax);
  BOOST_CHECK_EQUAL(frame.getEventResetCount(), liveResetCount);

  verifyFixtureLoaded(frame, fixture);
  BOOST_CHECK(frame.getGlobalMeasurements(LayerId{0}).data() == liveGlobalMeasurements);
  BOOST_CHECK(frame.getSurfaceMeasurements(LayerId{0}).data() == liveSurfaceMeasurements);
  BOOST_CHECK(frame.getClusters()[0].data() == liveClusters);
  BOOST_CHECK_EQUAL(frame.getClusters()[0].size(), liveClusterCount);
  BOOST_CHECK_EQUAL(frame.getClusters()[0][0].clusterId, liveClusterId);
  BOOST_CHECK(frame.getUsedClusters(0).data() == liveUsedClusterData);
  BOOST_CHECK(frame.getUsedClusters(0).empty());
  const auto preservedROFBoundaries = frame.getROFrameClusters(0);
  BOOST_CHECK(preservedROFBoundaries.data() == liveROFBoundaryData);
  BOOST_CHECK_EQUAL_COLLECTIONS(preservedROFBoundaries.begin(), preservedROFBoundaries.end(),
                                expectedROFBoundaries.begin(), expectedROFBoundaries.end());
  BOOST_REQUIRE_EQUAL(frame.getGenericTracks().size(), 1u);
  BOOST_CHECK(frame.getGenericTracks().data() == liveGenericTracks);
  BOOST_CHECK_EQUAL(frame.getGenericTracks()[0].chi2, 42.f);
  BOOST_REQUIRE_EQUAL(frame.getTrackClusterIndices().size(), 1u);
  BOOST_CHECK(frame.getTrackClusterIndices().data() == liveTrackClusterIndices);
  BOOST_CHECK(frame.getTrackClusterIndices()[0].surface == LayerId{0});
  BOOST_CHECK(frame.getTrackClusterIndices()[0].index == MeasurementIndex{0});
  BOOST_REQUIRE_EQUAL(frame.getPrimaryVerticesNum(), 1u);
  BOOST_CHECK(frame.getPrimaryVertices().data() == livePrimaryVertices);
}

BOOST_AUTO_TEST_CASE(MalformedReplacementValidationLeavesLoadedEventUnchanged)
{
  const auto catalog = makeITSTestCatalog();
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0};
  const auto baselineFixture = makeFixture();
  auto malformedReplacement = makeReplacementFixture();
  const auto plan = catalogLayout(catalogView, orderedSurfaces);
  TimeFrame frame;
  configureFrame(frame, catalogView, orderedSurfaces);
  const auto baseline = loadTimeFrameSource(frame, decoder, origin, timing, baselineFixture.clusters,
                                            baselineFixture.patterns, baselineFixture.rofs, &dict(),
                                            &baselineFixture.labels, o2::detectors::DetID::ITS,
                                            gsl::span<const LayerId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog());
  BOOST_REQUIRE(baseline.ok());
  verifyFixtureLoaded(frame, baselineFixture);

  malformedReplacement.rofs.front().setFirstEntry(1);
  const auto failed = loadTimeFrameSource(frame, decoder, origin, timing, malformedReplacement.clusters,
                                          malformedReplacement.patterns, malformedReplacement.rofs, &dict(),
                                          &malformedReplacement.labels, o2::detectors::DetID::ITS,
                                          gsl::span<const LayerId>{plan.getOrderedSurfaces()}, plan.getSurfaceCatalog());
  BOOST_CHECK(!failed.ok());
  BOOST_CHECK(failed.error == MultiSourceLoadError::InvalidROFRange);
  verifyFixtureLoaded(frame, baselineFixture);
}
