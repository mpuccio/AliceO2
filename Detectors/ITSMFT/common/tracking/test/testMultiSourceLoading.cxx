// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT MultiSourceLoading
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

// Host-only test decoder (no geometry singletons): maps a chip ID to a
// detector-local layer via an explicit table, and reuses the same pattern
// consumption path (extractClusterData) that the production decoder uses,
// so pattern-cursor bookkeeping is exercised identically.
class FakeClusterDecoder final : public ClusterDecoder
{
 public:
  FakeClusterDecoder(o2::detectors::DetID::ID detector, std::vector<int> sensorToLayer, bool disk)
    : mDetector(detector), mSensorToLayer(std::move(sensorToLayer)), mDisk(disk)
  {
  }

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
    float sigma2Row{0.f};
    float sigma2Col{0.f};
    ClusterShape shape{};
    o2::itsmft::ioutils::extractClusterData(cluster, pattIt, dict, sigma2Row, sigma2Col, nullptr, &shape);

    o2::itsmft::ioutils::SurfaceMeasurementDecodeResult result;
    const auto sensorID = cluster.getSensorID();
    const int layer = (sensorID >= 0 && static_cast<size_t>(sensorID) < mSensorToLayer.size()) ? mSensorToLayer[sensorID] : -1;
    result.layer = layer;
    if (layer < 0 || static_cast<size_t>(layer) >= layerToSurface.size()) {
      return result;
    }
    result.layerMapped = true;

    DecodedCluster decoded{};
    decoded.global = {static_cast<float>(sensorID), static_cast<float>(cluster.getRow()), static_cast<float>(cluster.getCol())};
    decoded.cylinderFrame = {10.f + sensorID, 1.f, 2.f, 0.1f};
    decoded.rowColumnCovariance = {sigma2Row, 0.f, sigma2Col};
    decoded.shape = shape;
    decoded.sensor = static_cast<uint32_t>(sensorID);
    decoded.layer = layer;

    const auto surface = layerToSurface[layer];
    const DetectorSensorId sensor{static_cast<uint32_t>(mDetector), decoded.sensor};
    const ClusterRef clusterRef{source, externalIndex};
    if (mDisk) {
      result.measurement = makeDiskSurfaceMeasurement(decoded, sensor, surface, clusterRef, sourceROF);
    } else {
      result.measurement = makeCylinderSurfaceMeasurement(decoded, sensor, surface, clusterRef, sourceROF);
    }
    return result;
  }

 private:
  o2::detectors::DetID::ID mDetector;
  std::vector<int> mSensorToLayer;
  bool mDisk;
};

// 4-surface disconnected ITS(cylinder)+MFT(disk) layout: surfaces {0,1} are
// ITS layers 0/1, surfaces {2,3} are MFT layers 0/1. No transitions are
// needed to exercise loading.
DetectorLayout makeCombinedLayout()
{
  SparseTrackingTopology topology{4};
  topology.finalize();
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.push_back(SurfaceDescriptor{SurfaceId{0}, 0, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{1}, 1, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{2}, 0, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{3}, 1, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
  return DetectorLayout{std::move(surfaces), std::move(topology)};
}

// One explicit (non-grouped) 1-pixel pattern: rowSpan=1, colSpan=1, one
// bitmap byte. Three bytes are consumed per cluster.
constexpr std::array<unsigned char, 3> onePixelPattern{1, 1, 0x80};

std::vector<unsigned char> makePatternBytes(size_t nClusters)
{
  std::vector<unsigned char> bytes;
  bytes.reserve(nClusters * onePixelPattern.size());
  for (size_t i = 0; i < nClusters; ++i) {
    bytes.insert(bytes.end(), onePixelPattern.begin(), onePixelPattern.end());
  }
  return bytes;
}

const TopologyDictionary& dict()
{
  static const TopologyDictionary d;
  return d;
}

const std::array<SurfaceId, 2> itsLayerToSurface{SurfaceId{0}, SurfaceId{1}};
const std::array<SurfaceId, 2> mftLayerToSurface{SurfaceId{2}, SurfaceId{3}};

} // namespace

BOOST_AUTO_TEST_CASE(SingleITSSourceLoadsIntoExpectedSurfaces)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> clusters{
    {10, 20, CompCluster::InvalidPatternID, 0}, // sensor 0 -> layer 0
    {11, 21, CompCluster::InvalidPatternID, 1}, // sensor 1 -> layer 1
  };
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 2}};

  FakeClusterDecoder decoder{o2::detectors::DetID::ITS, {0, 1}, false};
  ClusterSourceInput src;
  src.id = ClusterSourceId{0};
  src.detector = o2::detectors::DetID::ITS;
  src.clusters = clusters;
  src.patterns = patterns;
  src.rofs = rofs;
  src.dictionary = &dict();
  src.layerToSurface = itsLayerToSurface;
  src.timing = ROFTimingConfig{40, 0, 0, 0};
  src.decoder = &decoder;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_REQUIRE(result.ok());

  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{1}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{2}).size(), 0u);
  BOOST_REQUIRE_EQUAL(frame.getSources().size(), 1u);
  BOOST_CHECK(frame.getSources()[0].detector == o2::detectors::DetID::ITS);
  BOOST_CHECK_EQUAL(frame.getSources()[0].nROFs, 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0})[0].sensor.detector, static_cast<uint32_t>(o2::detectors::DetID::ITS));
}

BOOST_AUTO_TEST_CASE(SingleMFTSourceLoadsIntoExpectedSurfaces)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> clusters{
    {5, 6, CompCluster::InvalidPatternID, 0},
    {7, 8, CompCluster::InvalidPatternID, 1},
  };
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 2}};

  FakeClusterDecoder decoder{o2::detectors::DetID::MFT, {0, 1}, true};
  ClusterSourceInput src;
  src.id = ClusterSourceId{0};
  src.detector = o2::detectors::DetID::MFT;
  src.clusters = clusters;
  src.patterns = patterns;
  src.rofs = rofs;
  src.dictionary = &dict();
  src.layerToSurface = mftLayerToSurface;
  src.timing = ROFTimingConfig{40, 0, 0, 0};
  src.decoder = &decoder;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_REQUIRE(result.ok());

  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{2}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{3}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{2})[0].sensor.detector, static_cast<uint32_t>(o2::detectors::DetID::MFT));
}

BOOST_AUTO_TEST_CASE(CombinedITSAndMFTSourcesLoadTogether)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> itsClusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto itsPatterns = makePatternBytes(itsClusters.size());
  const std::vector<ROFRecord> itsRofs{ROFRecord{{0, 0}, 0, 0, 1}};
  FakeClusterDecoder itsDecoder{o2::detectors::DetID::ITS, {0}, false};

  const std::vector<CompClusterExt> mftClusters{{2, 2, CompCluster::InvalidPatternID, 0}};
  const auto mftPatterns = makePatternBytes(mftClusters.size());
  const std::vector<ROFRecord> mftRofs{ROFRecord{{0, 0}, 0, 0, 1}};
  FakeClusterDecoder mftDecoder{o2::detectors::DetID::MFT, {1}, true}; // sensor 0 -> layer 1 -> surface 3

  std::array<ClusterSourceInput, 2> sources{};
  sources[0].id = ClusterSourceId{0};
  sources[0].detector = o2::detectors::DetID::ITS;
  sources[0].clusters = itsClusters;
  sources[0].patterns = itsPatterns;
  sources[0].rofs = itsRofs;
  sources[0].dictionary = &dict();
  sources[0].layerToSurface = itsLayerToSurface;
  sources[0].timing = ROFTimingConfig{40, 0, 0, 0};
  sources[0].decoder = &itsDecoder;

  sources[1].id = ClusterSourceId{1};
  sources[1].detector = o2::detectors::DetID::MFT;
  sources[1].clusters = mftClusters;
  sources[1].patterns = mftPatterns;
  sources[1].rofs = mftRofs;
  sources[1].dictionary = &dict();
  sources[1].layerToSurface = mftLayerToSurface;
  sources[1].timing = ROFTimingConfig{50, 0, 0, 0};
  sources[1].decoder = &mftDecoder;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
  BOOST_REQUIRE(result.ok());

  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{3}).size(), 1u);
  BOOST_REQUIRE_EQUAL(frame.getSources().size(), 2u);
}

BOOST_AUTO_TEST_CASE(TwoSourcesOfSameDetectorBothAppendToOneSurface)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> clustersA{{1, 1, CompCluster::InvalidPatternID, 0}};
  const std::vector<CompClusterExt> clustersB{{2, 2, CompCluster::InvalidPatternID, 0}};
  const auto patternsA = makePatternBytes(clustersA.size());
  const auto patternsB = makePatternBytes(clustersB.size());
  const std::vector<ROFRecord> rofsA{ROFRecord{{0, 0}, 0, 0, 1}};
  const std::vector<ROFRecord> rofsB{ROFRecord{{0, 0}, 0, 0, 1}};

  FakeClusterDecoder decoderA{o2::detectors::DetID::ITS, {0}, false};
  FakeClusterDecoder decoderB{o2::detectors::DetID::ITS, {0}, false};

  std::array<ClusterSourceInput, 2> sources{};
  sources[0].id = ClusterSourceId{0};
  sources[0].detector = o2::detectors::DetID::ITS;
  sources[0].clusters = clustersA;
  sources[0].patterns = patternsA;
  sources[0].rofs = rofsA;
  sources[0].dictionary = &dict();
  sources[0].layerToSurface = itsLayerToSurface;
  sources[0].timing = ROFTimingConfig{40, 0, 0, 0};
  sources[0].decoder = &decoderA;

  sources[1].id = ClusterSourceId{1};
  sources[1].detector = o2::detectors::DetID::ITS;
  sources[1].clusters = clustersB;
  sources[1].patterns = patternsB;
  sources[1].rofs = rofsB;
  sources[1].dictionary = &dict();
  sources[1].layerToSurface = itsLayerToSurface;
  sources[1].timing = ROFTimingConfig{40, 0, 0, 0};
  sources[1].decoder = &decoderB;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
  BOOST_REQUIRE(result.ok());

  const auto onSurfaceZero = frame.getSurfaceMeasurements(SurfaceId{0});
  BOOST_REQUIRE_EQUAL(onSurfaceZero.size(), 2u);
  BOOST_CHECK(onSurfaceZero[0].cluster.source != onSurfaceZero[1].cluster.source);
}

BOOST_AUTO_TEST_CASE(IdenticalExternalIndicesInDifferentSourcesDoNotCollide)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> clustersA{{1, 1, CompCluster::InvalidPatternID, 0}}; // external index 0
  const std::vector<CompClusterExt> clustersB{{2, 2, CompCluster::InvalidPatternID, 0}}; // external index 0 too
  const auto patternsA = makePatternBytes(clustersA.size());
  const auto patternsB = makePatternBytes(clustersB.size());
  const std::vector<ROFRecord> rofsA{ROFRecord{{0, 0}, 0, 0, 1}};
  const std::vector<ROFRecord> rofsB{ROFRecord{{0, 0}, 0, 0, 1}};

  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labelsA;
  labelsA.addElement(0, o2::MCCompLabel{1, 0, 0});
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labelsB;
  labelsB.addElement(0, o2::MCCompLabel{2, 0, 0});

  FakeClusterDecoder decoderA{o2::detectors::DetID::ITS, {0}, false};
  FakeClusterDecoder decoderB{o2::detectors::DetID::ITS, {0}, false};

  std::array<ClusterSourceInput, 2> sources{};
  sources[0].id = ClusterSourceId{0};
  sources[0].detector = o2::detectors::DetID::ITS;
  sources[0].clusters = clustersA;
  sources[0].patterns = patternsA;
  sources[0].rofs = rofsA;
  sources[0].dictionary = &dict();
  sources[0].labels = &labelsA;
  sources[0].layerToSurface = itsLayerToSurface;
  sources[0].timing = ROFTimingConfig{40, 0, 0, 0};
  sources[0].decoder = &decoderA;

  sources[1].id = ClusterSourceId{1};
  sources[1].detector = o2::detectors::DetID::ITS;
  sources[1].clusters = clustersB;
  sources[1].patterns = patternsB;
  sources[1].rofs = rofsB;
  sources[1].dictionary = &dict();
  sources[1].labels = &labelsB;
  sources[1].layerToSurface = itsLayerToSurface;
  sources[1].timing = ROFTimingConfig{40, 0, 0, 0};
  sources[1].decoder = &decoderB;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
  BOOST_REQUIRE(result.ok());

  const auto onSurfaceZero = frame.getSurfaceMeasurements(SurfaceId{0});
  BOOST_REQUIRE_EQUAL(onSurfaceZero.size(), 2u);
  for (const auto& m : onSurfaceZero) {
    BOOST_CHECK_EQUAL(m.cluster.index, 0u);
  }

  const ClusterRef refA{ClusterSourceId{0}, 0};
  const ClusterRef refB{ClusterSourceId{1}, 0};
  const auto labelSpanA = frame.getLabels(refA);
  const auto labelSpanB = frame.getLabels(refB);
  BOOST_REQUIRE_EQUAL(labelSpanA.size(), 1u);
  BOOST_REQUIRE_EQUAL(labelSpanB.size(), 1u);
  BOOST_CHECK(labelSpanA[0] != labelSpanB[0]);
}

BOOST_AUTO_TEST_CASE(IndependentROFCountsAcrossSourcesAreAllowed)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  // Source A: 3 ROFs of 1 cluster each. Source B: 1 ROF of 1 cluster.
  const std::vector<CompClusterExt> clustersA{
    {1, 1, CompCluster::InvalidPatternID, 0},
    {2, 2, CompCluster::InvalidPatternID, 0},
    {3, 3, CompCluster::InvalidPatternID, 0}};
  const auto patternsA = makePatternBytes(clustersA.size());
  const std::vector<ROFRecord> rofsA{
    ROFRecord{{0, 0}, 0, 0, 1},
    ROFRecord{{40, 0}, 1, 1, 1},
    ROFRecord{{80, 0}, 2, 2, 1}};

  const std::vector<CompClusterExt> clustersB{{4, 4, CompCluster::InvalidPatternID, 1}};
  const auto patternsB = makePatternBytes(clustersB.size());
  const std::vector<ROFRecord> rofsB{ROFRecord{{0, 0}, 0, 0, 1}};

  FakeClusterDecoder decoderA{o2::detectors::DetID::ITS, {0}, false};
  FakeClusterDecoder decoderB{o2::detectors::DetID::ITS, {-1, 1}, false}; // sensor 1 -> layer 1

  std::array<ClusterSourceInput, 2> sources{};
  sources[0].id = ClusterSourceId{0};
  sources[0].detector = o2::detectors::DetID::ITS;
  sources[0].clusters = clustersA;
  sources[0].patterns = patternsA;
  sources[0].rofs = rofsA;
  sources[0].dictionary = &dict();
  sources[0].layerToSurface = itsLayerToSurface;
  sources[0].timing = ROFTimingConfig{40, 0, 0, 0};
  sources[0].decoder = &decoderA;

  sources[1].id = ClusterSourceId{1};
  sources[1].detector = o2::detectors::DetID::ITS;
  sources[1].clusters = clustersB;
  sources[1].patterns = patternsB;
  sources[1].rofs = rofsB;
  sources[1].dictionary = &dict();
  sources[1].layerToSurface = itsLayerToSurface;
  sources[1].timing = ROFTimingConfig{100, 0, 0, 0};
  sources[1].decoder = &decoderB;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
  BOOST_REQUIRE(result.ok());

  BOOST_CHECK_EQUAL(frame.getSources()[0].nROFs, 3u);
  BOOST_CHECK_EQUAL(frame.getSources()[1].nROFs, 1u);
  BOOST_CHECK_EQUAL(frame.getSourceIntervals(ClusterSourceId{0}).size(), 3u);
  BOOST_CHECK_EQUAL(frame.getSourceIntervals(ClusterSourceId{1}).size(), 1u);
}

BOOST_AUTO_TEST_CASE(OverlappingAndNonOverlappingSourceTimingIntervals)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> clustersA{{1, 1, CompCluster::InvalidPatternID, 0}};
  const std::vector<CompClusterExt> clustersB{{2, 2, CompCluster::InvalidPatternID, 0}};
  const auto patternsA = makePatternBytes(clustersA.size());
  const auto patternsB = makePatternBytes(clustersB.size());
  // Source A ROF at BC 0..40 (TF-relative); source B ROF at real BC 30 -> its
  // own interval overlaps A's despite a different, unrelated ROF ordinal.
  const std::vector<ROFRecord> rofsA{ROFRecord{{0, 0}, 0, 0, 1}};
  const std::vector<ROFRecord> rofsB{ROFRecord{{30, 0}, 0, 0, 1}};
  // Source C ROF at real BC 1000: far away, must not overlap A.
  const std::vector<CompClusterExt> clustersC{{3, 3, CompCluster::InvalidPatternID, 0}};
  const auto patternsC = makePatternBytes(clustersC.size());
  const std::vector<ROFRecord> rofsC{ROFRecord{{1000, 0}, 0, 0, 1}};

  FakeClusterDecoder decoderA{o2::detectors::DetID::ITS, {0}, false};
  FakeClusterDecoder decoderB{o2::detectors::DetID::ITS, {0}, false};
  FakeClusterDecoder decoderC{o2::detectors::DetID::ITS, {0}, false};

  std::array<ClusterSourceInput, 3> sources{};
  sources[0].id = ClusterSourceId{0};
  sources[0].detector = o2::detectors::DetID::ITS;
  sources[0].clusters = clustersA;
  sources[0].patterns = patternsA;
  sources[0].rofs = rofsA;
  sources[0].dictionary = &dict();
  sources[0].layerToSurface = itsLayerToSurface;
  sources[0].timing = ROFTimingConfig{40, 0, 0, 0};
  sources[0].decoder = &decoderA;

  sources[1].id = ClusterSourceId{1};
  sources[1].detector = o2::detectors::DetID::ITS;
  sources[1].clusters = clustersB;
  sources[1].patterns = patternsB;
  sources[1].rofs = rofsB;
  sources[1].dictionary = &dict();
  sources[1].layerToSurface = itsLayerToSurface;
  sources[1].timing = ROFTimingConfig{40, 0, 0, 0};
  sources[1].decoder = &decoderB;

  sources[2].id = ClusterSourceId{2};
  sources[2].detector = o2::detectors::DetID::ITS;
  sources[2].clusters = clustersC;
  sources[2].patterns = patternsC;
  sources[2].rofs = rofsC;
  sources[2].dictionary = &dict();
  sources[2].layerToSurface = itsLayerToSurface;
  sources[2].timing = ROFTimingConfig{40, 0, 0, 0};
  sources[2].decoder = &decoderC;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
  BOOST_REQUIRE(result.ok());

  const auto a = frame.getSourceIntervals(ClusterSourceId{0})[0];
  const auto b = frame.getSourceIntervals(ClusterSourceId{1})[0];
  const auto c = frame.getSourceIntervals(ClusterSourceId{2})[0];
  BOOST_CHECK(intersects(a, b));
  BOOST_CHECK(!intersects(a, c));
}

BOOST_AUTO_TEST_CASE(SourceSpecificPatternCursorsAreIndependent)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> clustersA{
    {1, 1, CompCluster::InvalidPatternID, 0},
    {2, 2, CompCluster::InvalidPatternID, 0}};
  const std::vector<CompClusterExt> clustersB{
    {3, 3, CompCluster::InvalidPatternID, 0},
    {4, 4, CompCluster::InvalidPatternID, 0}};
  const auto patternsA = makePatternBytes(clustersA.size());
  const auto patternsB = makePatternBytes(clustersB.size());
  const std::vector<ROFRecord> rofsA{ROFRecord{{0, 0}, 0, 0, 2}};
  const std::vector<ROFRecord> rofsB{ROFRecord{{0, 0}, 0, 0, 2}};

  FakeClusterDecoder decoderA{o2::detectors::DetID::ITS, {0}, false};
  FakeClusterDecoder decoderB{o2::detectors::DetID::ITS, {0}, false};

  std::array<ClusterSourceInput, 2> sources{};
  sources[0].id = ClusterSourceId{0};
  sources[0].detector = o2::detectors::DetID::ITS;
  sources[0].clusters = clustersA;
  sources[0].patterns = patternsA;
  sources[0].rofs = rofsA;
  sources[0].dictionary = &dict();
  sources[0].layerToSurface = itsLayerToSurface;
  sources[0].timing = ROFTimingConfig{40, 0, 0, 0};
  sources[0].decoder = &decoderA;

  sources[1].id = ClusterSourceId{1};
  sources[1].detector = o2::detectors::DetID::ITS;
  sources[1].clusters = clustersB;
  sources[1].patterns = patternsB;
  sources[1].rofs = rofsB;
  sources[1].dictionary = &dict();
  sources[1].layerToSurface = itsLayerToSurface;
  sources[1].timing = ROFTimingConfig{40, 0, 0, 0};
  sources[1].decoder = &decoderB;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
  BOOST_REQUIRE(result.ok());

  // Every cluster consumed exactly one 1-pixel pattern regardless of source.
  for (const auto& m : frame.getSurfaceMeasurements(SurfaceId{0})) {
    BOOST_CHECK_EQUAL(m.shape.nPixels, 1u);
  }
}

BOOST_AUTO_TEST_CASE(AbsentLabelsAreLegal)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  FakeClusterDecoder decoder{o2::detectors::DetID::ITS, {0}, false};
  ClusterSourceInput src;
  src.id = ClusterSourceId{0};
  src.detector = o2::detectors::DetID::ITS;
  src.clusters = clusters;
  src.patterns = patterns;
  src.rofs = rofs;
  src.dictionary = &dict();
  src.labels = nullptr; // no MC labels for this source
  src.layerToSurface = itsLayerToSurface;
  src.timing = ROFTimingConfig{40, 0, 0, 0};
  src.decoder = &decoder;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_REQUIRE(result.ok());

  BOOST_CHECK(frame.getLabels(ClusterRef{ClusterSourceId{0}, 0}).empty());
  BOOST_CHECK(frame.getLabels(ClusterRef{}).empty()); // invalid ref is also legal
}

BOOST_AUTO_TEST_CASE(NonDenseAndDuplicateAndInvalidSourceIdsAreRejected)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
  FakeClusterDecoder decoderA{o2::detectors::DetID::ITS, {0}, false};
  FakeClusterDecoder decoderB{o2::detectors::DetID::ITS, {0}, false};

  auto makeSource = [&](ClusterSourceId id, FakeClusterDecoder& decoder) {
    ClusterSourceInput src;
    src.id = id;
    src.detector = o2::detectors::DetID::ITS;
    src.clusters = clusters;
    src.patterns = patterns;
    src.rofs = rofs;
    src.dictionary = &dict();
    src.layerToSurface = itsLayerToSurface;
    src.timing = ROFTimingConfig{40, 0, 0, 0};
    src.decoder = &decoder;
    return src;
  };

  {
    // Non-dense: ids {0, 2} for two sources.
    std::array<ClusterSourceInput, 2> sources{makeSource(ClusterSourceId{0}, decoderA), makeSource(ClusterSourceId{2}, decoderB)};
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::NonDenseSourceIds);
  }
  {
    // Duplicate ids {0, 0}.
    std::array<ClusterSourceInput, 2> sources{makeSource(ClusterSourceId{0}, decoderA), makeSource(ClusterSourceId{0}, decoderB)};
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::DuplicateSourceId);
  }
  {
    // Explicitly invalid id.
    std::array<ClusterSourceInput, 1> sources{makeSource(ClusterSourceId::invalid(), decoderA)};
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::NonDenseSourceIds);
  }
}

BOOST_AUTO_TEST_CASE(InvalidROFClusterRangesAreRejected)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());
  const std::vector<CompClusterExt> clusters{
    {1, 1, CompCluster::InvalidPatternID, 0},
    {2, 2, CompCluster::InvalidPatternID, 0}};
  const auto patterns = makePatternBytes(clusters.size());
  FakeClusterDecoder decoder{o2::detectors::DetID::ITS, {0}, false};

  auto makeSrc = [&](const std::vector<ROFRecord>& rofs) {
    ClusterSourceInput src;
    src.id = ClusterSourceId{0};
    src.detector = o2::detectors::DetID::ITS;
    src.clusters = clusters;
    src.patterns = patterns;
    src.rofs = rofs;
    src.dictionary = &dict();
    src.layerToSurface = itsLayerToSurface;
    src.timing = ROFTimingConfig{40, 0, 0, 0};
    src.decoder = &decoder;
    return src;
  };

  {
    // Out of bounds: firstEntry+nEntries exceeds the cluster span.
    const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 5}};
    auto src = makeSrc(rofs);
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidROFRange);
  }
  {
    // Overlapping ranges.
    const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 2}, ROFRecord{{40, 0}, 1, 1, 1}};
    auto src = makeSrc(rofs);
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidROFRange);
  }
}

BOOST_AUTO_TEST_CASE(InvalidLayerToSurfaceMappingIsRejected)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 1}}; // sensor 1 -> layer 1
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
  FakeClusterDecoder decoder{o2::detectors::DetID::ITS, {-1, 1}, false};

  ClusterSourceInput src;
  src.id = ClusterSourceId{0};
  src.detector = o2::detectors::DetID::ITS;
  src.clusters = clusters;
  src.patterns = patterns;
  src.rofs = rofs;
  src.dictionary = &dict();
  src.layerToSurface = gsl::span<const SurfaceId>(itsLayerToSurface.data(), 1); // too short: only covers layer 0
  src.timing = ROFTimingConfig{40, 0, 0, 0};
  src.decoder = &decoder;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);
}

BOOST_AUTO_TEST_CASE(DetectorSurfaceMismatchIsRejected)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
  FakeClusterDecoder decoder{o2::detectors::DetID::ITS, {0}, false};

  ClusterSourceInput src;
  src.id = ClusterSourceId{0};
  src.detector = o2::detectors::DetID::ITS;
  src.clusters = clusters;
  src.patterns = patterns;
  src.rofs = rofs;
  src.dictionary = &dict();
  // Deliberately mapped to an MFT surface: ITS source, MFT surface.
  const std::array<SurfaceId, 1> wrongMapping{SurfaceId{2}};
  src.layerToSurface = wrongMapping;
  src.timing = ROFTimingConfig{40, 0, 0, 0};
  src.decoder = &decoder;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::DetectorSurfaceMismatch);
}

BOOST_AUTO_TEST_CASE(FailedLoadLeavesNoPartialState)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
  FakeClusterDecoder decoder{o2::detectors::DetID::ITS, {0}, false};

  ClusterSourceInput goodSrc;
  goodSrc.id = ClusterSourceId{0};
  goodSrc.detector = o2::detectors::DetID::ITS;
  goodSrc.clusters = clusters;
  goodSrc.patterns = patterns;
  goodSrc.rofs = rofs;
  goodSrc.dictionary = &dict();
  goodSrc.layerToSurface = itsLayerToSurface;
  goodSrc.timing = ROFTimingConfig{40, 0, 0, 0};
  goodSrc.decoder = &decoder;

  MultiSourceFrame frame;
  BOOST_REQUIRE(loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(&goodSrc, 1), {0, 0}).ok());
  const auto totalBefore = frame.getTotalMeasurements();
  BOOST_REQUIRE_EQUAL(totalBefore, 1u);

  // Now attempt an invalid load (duplicate ids) on the SAME frame.
  std::array<ClusterSourceInput, 2> badSources{goodSrc, goodSrc}; // both id==0
  const auto result = loadSources(frame, layout.getView(), gsl::span<const ClusterSourceInput>(badSources), {0, 0});
  BOOST_REQUIRE(!result.ok());

  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), totalBefore);
  BOOST_REQUIRE_EQUAL(frame.getSources().size(), 1u);
  BOOST_CHECK(frame.getSources()[0].detector == o2::detectors::DetID::ITS);
}

BOOST_AUTO_TEST_CASE(ViewsAreStandardLayoutAndTriviallyCopyable)
{
  static_assert(std::is_standard_layout_v<MultiSourceFrameView>);
  static_assert(std::is_trivially_copyable_v<MultiSourceFrameView>);
  static_assert(std::is_standard_layout_v<SurfaceMeasurementRange>);
  static_assert(std::is_trivially_copyable_v<SurfaceMeasurementRange>);
  BOOST_CHECK(std::is_standard_layout_v<MultiSourceFrameView>);
  BOOST_CHECK(std::is_trivially_copyable_v<MultiSourceFrameView>);
}
