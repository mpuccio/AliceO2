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
#include "ITSMFTTracking/SurfaceGraph.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

// Host-only test decoder (no geometry singletons): maps a chip ID to a
// detector-local layer via an explicit table, and reuses the same pattern
// consumption path (extractClusterData) that the production decoder uses,
// so pattern-cursor bookkeeping is exercised identically. An optional
// `Corruption` deliberately breaks one authoritative-identity field, to
// prove that loadSources() itself catches a buggy host adapter rather than
// trusting the decoder's output blindly.
enum class Corruption {
  None,
  WrongSurface,
  WrongSource,
  WrongIndex,
  WrongSourceROF,
  WrongSensorDetector,
  WrongKind,
  // These two report layerMapped=true while `layer` itself is not a usable
  // index into layerToSurface: loadSources() must not trust layerMapped
  // alone before indexing with `layer`.
  LayerMappedTrueWithNegativeLayer,
  LayerMappedTrueWithLayerOutOfRange
};

class FakeClusterDecoder final : public ClusterDecoder
{
 public:
  FakeClusterDecoder(o2::detectors::DetID::ID detector, std::vector<int> sensorToLayer, bool disk, Corruption corruption = Corruption::None)
    : mDetector(detector), mSensorToLayer(std::move(sensorToLayer)), mDisk(disk), mCorruption(corruption)
  {
  }

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
    // Simulate a decoder that lies about layerMapped without ever touching
    // the pattern cursor or geometry decode: loadSources() must catch this
    // from layerMapped/layer alone, before indexing layerToSurface.
    if (mCorruption == Corruption::LayerMappedTrueWithNegativeLayer) {
      o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
      result.layerMapped = true;
      result.layer = -1;
      return result;
    }
    if (mCorruption == Corruption::LayerMappedTrueWithLayerOutOfRange) {
      o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
      result.layerMapped = true;
      result.layer = static_cast<int>(layerToSurface.size());
      return result;
    }

    const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dict);
    if (!clusterData.ok()) {
      o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
      result.error = clusterData.error;
      return result;
    }

    o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
    const auto sensorID = cluster.getSensorID();
    const int layer = (sensorID >= 0 && static_cast<size_t>(sensorID) < mSensorToLayer.size()) ? mSensorToLayer[sensorID] : -1;
    result.layer = layer;
    if (layer < 0 || static_cast<size_t>(layer) >= layerToSurface.size()) {
      return result;
    }
    result.layerMapped = true;
    result.kind = mDisk ? SurfaceKind::Disk : SurfaceKind::Cylinder;

    DecodedCluster decoded{};
    decoded.global = {static_cast<float>(sensorID), static_cast<float>(cluster.getRow()), static_cast<float>(cluster.getCol())};
    decoded.cylinderFrame = {10.f + sensorID, 1.f, 2.f, 0.1f};
    decoded.rowColumnCovariance = {clusterData.sig2Row, 0.f, clusterData.sig2Col};
    decoded.shape = clusterData.shape;
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

    switch (mCorruption) {
      case Corruption::WrongSurface:
        result.measurement.surface = SurfaceId{static_cast<uint16_t>(surface.value() == 0 ? 1 : 0)};
        break;
      case Corruption::WrongSource:
        result.measurement.cluster.source = ClusterSourceId{static_cast<uint16_t>(source.value() + 7)};
        break;
      case Corruption::WrongIndex:
        result.measurement.cluster.index = externalIndex + 1;
        break;
      case Corruption::WrongSourceROF:
        result.measurement.sourceROF = sourceROF + 1;
        break;
      case Corruption::WrongSensorDetector:
        result.measurement.sensor.detector = static_cast<uint32_t>(
          mDetector == o2::detectors::DetID::ITS ? o2::detectors::DetID::MFT : o2::detectors::DetID::ITS);
        break;
      case Corruption::WrongKind:
        result.kind = (result.kind == SurfaceKind::Cylinder) ? SurfaceKind::Disk : SurfaceKind::Cylinder;
        break;
      case Corruption::LayerMappedTrueWithNegativeLayer:
      case Corruption::LayerMappedTrueWithLayerOutOfRange:
        // Handled by the early returns above, before normal decoding.
        break;
      case Corruption::None:
        break;
    }
    return result;
  }

 private:
  o2::detectors::DetID::ID mDetector;
  std::vector<int> mSensorToLayer;
  bool mDisk;
  Corruption mCorruption;
};

// Geometry-free decoder used only to exercise the normalized loader's
// dictionary/common/group/explicit pattern contract. Pattern ID 0 represents
// a common dictionary entry (no explicit bytes), pattern ID 1 represents a
// grouped dictionary entry (explicit bytes required), and InvalidPatternID
// represents an ordinary explicit pattern.
class PatternContractDecoder final : public ClusterDecoder
{
 public:
  o2::itsmft::tracking::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const TopologyDictionary* dictionary,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool) const override
  {
    o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
    if (dictionary == nullptr) {
      result.error = ClusterDecodeError::MissingDictionary;
      return result;
    }
    if (layerToSurface.empty()) {
      result.error = ClusterDecodeError::InvalidLayerMapping;
      return result;
    }

    ClusterShape shape{1, 1, 1};
    if (cluster.getPatternID() != 0) {
      ClusterPattern pattern;
      result.error = patterns.acquirePattern(pattern);
      if (!result.ok()) {
        return result;
      }
      shape = {static_cast<uint32_t>(pattern.getNPixels()),
               static_cast<uint16_t>(pattern.getRowSpan()),
               static_cast<uint16_t>(pattern.getColumnSpan())};
    }

    DecodedCluster decoded{};
    decoded.global = {1.f, 2.f, 3.f};
    decoded.cylinderFrame = {4.f, 5.f, 6.f, 0.f};
    decoded.rowColumnCovariance = {0.1f, 0.f, 0.2f};
    decoded.shape = shape;
    decoded.sensor = static_cast<uint32_t>(cluster.getSensorID());
    decoded.layer = 0;
    result.layer = 0;
    result.layerMapped = true;
    result.kind = SurfaceKind::Cylinder;
    result.measurement = makeCylinderSurfaceMeasurement(
      decoded, {o2::detectors::DetID::ITS, decoded.sensor}, layerToSurface[0],
      {source, externalIndex}, sourceROF);
    return result;
  }
};

/// SurfaceGraph no longer owns a surface copy (Slice 3, shared ownership):
/// it borrows a caller-supplied catalog span only for construction/view
/// assembly. This fixture keeps its own catalog alongside it so `.getView()`
/// keeps working as a zero-argument call at every existing call site below.
struct BuiltLayout {
  SurfaceGraph layout;
  std::vector<SurfaceDescriptor> surfaces;

  bool valid() const noexcept { return layout.valid(); }
  SurfaceGraphView getView() const noexcept
  {
    const auto masks = computeSurfaceKindMasks(surfaces);
    return layout.getView();
  }
};

// 4-surface disconnected ITS(cylinder)+MFT(disk) layout: surfaces {0,1} are
// ITS layers 0/1, surfaces {2,3} are MFT layers 0/1. No transitions are
// needed to exercise loading.
BuiltLayout makeCombinedLayout()
{
  SurfaceGraph topology{4};
  topology.finalize();
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.push_back(SurfaceDescriptor{SurfaceId{0}, 0, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{1}, 1, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{2}, 0, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{3}, 1, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
  return BuiltLayout{SurfaceGraph{surfaces, std::move(topology)}, std::move(surfaces)};
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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_REQUIRE(result.ok());
  // A success result must retain the timingDetail default: it is only ever
  // meaningful when error == TimingError.
  BOOST_CHECK(result.timingDetail == TimingBuildError::None);

  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{1}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{2}).size(), 0u);
  BOOST_REQUIRE_EQUAL(frame.getSources().size(), 1u);
  BOOST_CHECK(frame.getSources()[0].detector == o2::detectors::DetID::ITS);
  BOOST_CHECK_EQUAL(frame.getSources()[0].nROFs, 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0})[0].sensor.detector, static_cast<uint32_t>(o2::detectors::DetID::ITS));
}

BOOST_AUTO_TEST_CASE(InvalidTimingConfigurationIsReportedWithBuildErrorDetail)
{
  // computeROFIntervalBC()'s own exhaustive TimingBuildError coverage lives
  // in testSurfaceTiming.cxx (InvalidROFLengthIsRejected, OverflowIsDetected
  // AndChecked, InvalidSourceROFIsRejected); this test only proves that
  // loadSources() actually plumbs that detail into LoadSourcesResult rather
  // than discarding it. InvalidROFLength (rofLength <= 0) is the only one of
  // the three practically reachable through loadSources() itself:
  // InvalidSourceROF would require a source ROF count exceeding UINT32_MAX,
  // and Overflow requires contrived BC values already covered directly at
  // the computeROFIntervalBC() level.
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
  src.layerToSurface = itsLayerToSurface;
  src.timing = ROFTimingConfig{0, 0, 0, 0}; // rofLength <= 0
  src.decoder = &decoder;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_CHECK(result.error == MultiSourceLoadError::TimingError);
  BOOST_CHECK(result.timingDetail == TimingBuildError::InvalidROFLength);
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
  BOOST_REQUIRE(result.ok());

  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{3}).size(), 1u);
  BOOST_REQUIRE_EQUAL(frame.getSources().size(), 2u);
}

// Hardening regression test (MultiSourceFrame cached-span safety): proves
// MultiSourceFrameView's cached per-surface spans (mSurfaceSpans) resolve
// correctly, on more than one populated surface, immediately after a real,
// successful load -- not just via the owner-side accessors every other
// test above already exercises. Also cross-checks the view's resolution
// against the owner's own per-surface accessor, proving the cache and the
// underlying storage never diverge after a real commit.
BOOST_AUTO_TEST_CASE(ViewResolvesMeasurementsOnMultipleSurfacesAfterSuccessfulLoad)
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
  BOOST_REQUIRE(loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0}).ok());

  const auto view = frame.getView();
  BOOST_REQUIRE_EQUAL(view.nSurfaces, 4u);
  BOOST_CHECK_EQUAL(view.getSurfaceMeasurementCount(SurfaceId{0}), 1u);
  BOOST_CHECK_EQUAL(view.getSurfaceMeasurementCount(SurfaceId{3}), 1u);

  const auto* onITS = view.getMeasurement(SurfaceId{0}, SurfaceMeasurementIndex{0});
  const auto* onMFT = view.getMeasurement(SurfaceId{3}, SurfaceMeasurementIndex{0});
  BOOST_REQUIRE(onITS != nullptr);
  BOOST_REQUIRE(onMFT != nullptr);
  BOOST_CHECK(onITS->surface == SurfaceId{0});
  BOOST_CHECK(onMFT->surface == SurfaceId{3});
  BOOST_CHECK(onITS->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::ITS));
  BOOST_CHECK(onMFT->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::MFT));

  BOOST_CHECK_EQUAL(onITS->global.x, frame.getSurfaceMeasurements(SurfaceId{0})[0].global.x);
  BOOST_CHECK_EQUAL(onMFT->global.x, frame.getSurfaceMeasurements(SurfaceId{3})[0].global.x);
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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
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

BOOST_AUTO_TEST_CASE(ClusterRefFlagsDoNotAffectIdentityOrLabelLookup)
{
  // {ClusterSourceId, external index} is the whole identity; flags are a
  // side channel that operator== and label lookup must ignore entirely.
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
  labels.addElement(0, o2::MCCompLabel{1, 0, 0});

  FakeClusterDecoder decoder{o2::detectors::DetID::ITS, {0}, false};
  ClusterSourceInput src;
  src.id = ClusterSourceId{0};
  src.detector = o2::detectors::DetID::ITS;
  src.clusters = clusters;
  src.patterns = patterns;
  src.rofs = rofs;
  src.dictionary = &dict();
  src.labels = &labels;
  src.layerToSurface = itsLayerToSurface;
  src.timing = ROFTimingConfig{40, 0, 0, 0};
  src.decoder = &decoder;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_REQUIRE(result.ok());

  const ClusterRef plain{ClusterSourceId{0}, 0};
  const ClusterRef flagged{ClusterSourceId{0}, 0, 0xdead};
  BOOST_CHECK(plain == flagged);
  BOOST_CHECK(!(plain != flagged));

  const auto labelPlain = frame.getLabels(plain);
  const auto labelFlagged = frame.getLabels(flagged);
  BOOST_REQUIRE_EQUAL(labelPlain.size(), 1u);
  BOOST_REQUIRE_EQUAL(labelFlagged.size(), 1u);
  BOOST_CHECK(labelPlain[0] == labelFlagged[0]);

  // Two refs that differ only in flags but disagree on source or index are
  // still, correctly, not the same identity.
  BOOST_CHECK(plain != ClusterRef(ClusterSourceId{0}, 1, 0xdead));
  BOOST_CHECK(plain != ClusterRef(ClusterSourceId{1}, 0, 0xdead));

  // The measurement's own stored cluster ref -- produced by the production
  // decode path -- must also compare equal regardless of any flags a caller
  // later probes it with.
  const auto measurement = frame.getSurfaceMeasurements(SurfaceId{0})[0];
  BOOST_CHECK(measurement.cluster == flagged);
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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
  BOOST_REQUIRE(result.ok());

  const auto a = frame.getSourceIntervals(ClusterSourceId{0})[0];
  const auto b = frame.getSourceIntervals(ClusterSourceId{1})[0];
  const auto c = frame.getSourceIntervals(ClusterSourceId{2})[0];
  BOOST_CHECK(intersects(a, b));
  BOOST_CHECK(!intersects(a, c));
}

BOOST_AUTO_TEST_CASE(TriggeredAndContinuousReadoutAreBothSupportedTogether)
{
  // Continuous source: ROFs sit at a fixed cadence equal to the readout
  // length, so consecutive interval begins are regularly spaced by
  // rofLength (mirrors a periodic strobe). Triggered source: ROFs sit at
  // sparse, irregular real interaction records (individual triggers) with a
  // short trigger-specific window, so consecutive interval begins follow the
  // trigger BCs exactly rather than any ordinal*rofLength formula. Both must
  // load into the same frame and their intervals must remain independently
  // and correctly comparable via intersection.
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> continuousClusters{
    {1, 1, CompCluster::InvalidPatternID, 0},
    {2, 2, CompCluster::InvalidPatternID, 0},
    {3, 3, CompCluster::InvalidPatternID, 0}};
  const auto continuousPatterns = makePatternBytes(continuousClusters.size());
  const std::vector<ROFRecord> continuousRofs{
    ROFRecord{{0, 0}, 0, 0, 1},
    ROFRecord{{40, 0}, 1, 1, 1},
    ROFRecord{{80, 0}, 2, 2, 1}};
  constexpr TFBC continuousRofLength = 40;

  const std::vector<CompClusterExt> triggeredClusters{
    {4, 4, CompCluster::InvalidPatternID, 0},
    {5, 5, CompCluster::InvalidPatternID, 0},
    {6, 6, CompCluster::InvalidPatternID, 0}};
  const auto triggeredPatterns = makePatternBytes(triggeredClusters.size());
  // Sparse, irregular trigger BCs; a short single-BC-scale trigger window.
  const std::vector<ROFRecord> triggeredRofs{
    ROFRecord{{5, 0}, 0, 0, 1},
    ROFRecord{{137, 0}, 1, 1, 1},
    ROFRecord{{812, 0}, 2, 2, 1}};
  constexpr TFBC triggeredRofLength = 4;

  FakeClusterDecoder continuousDecoder{o2::detectors::DetID::ITS, {0}, false};
  FakeClusterDecoder triggeredDecoder{o2::detectors::DetID::ITS, {0}, false};

  std::array<ClusterSourceInput, 2> sources{};
  sources[0].id = ClusterSourceId{0};
  sources[0].detector = o2::detectors::DetID::ITS;
  sources[0].clusters = continuousClusters;
  sources[0].patterns = continuousPatterns;
  sources[0].rofs = continuousRofs;
  sources[0].dictionary = &dict();
  sources[0].layerToSurface = itsLayerToSurface;
  sources[0].timing = ROFTimingConfig{continuousRofLength, 0, 0, 0};
  sources[0].decoder = &continuousDecoder;

  sources[1].id = ClusterSourceId{1};
  sources[1].detector = o2::detectors::DetID::ITS;
  sources[1].clusters = triggeredClusters;
  sources[1].patterns = triggeredPatterns;
  sources[1].rofs = triggeredRofs;
  sources[1].dictionary = &dict();
  sources[1].layerToSurface = itsLayerToSurface;
  sources[1].timing = ROFTimingConfig{triggeredRofLength, 0, 0, 0};
  sources[1].decoder = &triggeredDecoder;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
  BOOST_REQUIRE(result.ok());

  const auto continuousIntervals = frame.getSourceIntervals(ClusterSourceId{0});
  BOOST_REQUIRE_EQUAL(continuousIntervals.size(), 3u);
  for (size_t i = 1; i < continuousIntervals.size(); ++i) {
    BOOST_CHECK_EQUAL(continuousIntervals[i].begin - continuousIntervals[i - 1].begin, continuousRofLength);
  }

  const auto triggeredIntervals = frame.getSourceIntervals(ClusterSourceId{1});
  BOOST_REQUIRE_EQUAL(triggeredIntervals.size(), 3u);
  // The triggered ROFs' own irregular real BCs are reproduced exactly (5,
  // 137, 812), not forced onto a uniform cadence.
  BOOST_CHECK_EQUAL(triggeredIntervals[0].begin, 5);
  BOOST_CHECK_EQUAL(triggeredIntervals[1].begin, 137);
  BOOST_CHECK_EQUAL(triggeredIntervals[2].begin, 812);
  for (const auto& interval : triggeredIntervals) {
    BOOST_CHECK_EQUAL(interval.length(), triggeredRofLength);
  }

  // Cross-checked overlap: the continuous source's second ROF spans [40,80);
  // it must not spuriously overlap the nearby-but-disjoint triggered ROF at
  // [137,141), and must correctly not overlap the far-away one at [812,816).
  BOOST_CHECK(!intersects(continuousIntervals[1], triggeredIntervals[1]));
  BOOST_CHECK(!intersects(continuousIntervals[0], triggeredIntervals[2]));
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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
  BOOST_REQUIRE(result.ok());

  // Every cluster consumed exactly one 1-pixel pattern regardless of source.
  for (const auto& m : frame.getSurfaceMeasurements(SurfaceId{0})) {
    BOOST_CHECK_EQUAL(m.shape.nPixels, 1u);
  }
}

BOOST_AUTO_TEST_CASE(CommonDictionaryPatternDoesNotConsumeExplicitBytes)
{
  const auto layout = makeCombinedLayout();
  const std::vector<CompClusterExt> clusters{
    {1, 1, 0, 0},                              // common dictionary pattern
    {2, 2, CompCluster::InvalidPatternID, 0}}; // explicit pattern
  const std::vector<unsigned char> patterns{onePixelPattern.begin(), onePixelPattern.end()};
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 2}};
  PatternContractDecoder decoder;

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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0}).size(), 2u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0})[0].shape.nPixels, 1u);
  BOOST_CHECK_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0})[1].shape.nPixels, 1u);
}

BOOST_AUTO_TEST_CASE(ExplicitAndGroupedPatternTruncationIsTypedAndContextual)
{
  const auto layout = makeCombinedLayout();
  PatternContractDecoder decoder;
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
  constexpr std::array<unsigned char, 4> encoded{3, 3, 0x80, 0x80};

  for (const auto patternID : {CompCluster::InvalidPatternID, static_cast<unsigned short>(1)}) {
    const std::vector<CompClusterExt> clusters{{1, 1, patternID, 0}};
    for (size_t available = 0; available < encoded.size(); ++available) {
      ClusterSourceInput src;
      src.id = ClusterSourceId{0};
      src.detector = o2::detectors::DetID::ITS;
      src.clusters = clusters;
      src.patterns = gsl::span<const unsigned char>{encoded.data(), available};
      src.rofs = rofs;
      src.dictionary = &dict();
      src.layerToSurface = itsLayerToSurface;
      src.timing = ROFTimingConfig{40, 0, 0, 0};
      src.decoder = &decoder;

      MultiSourceFrame frame;
      const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
      BOOST_CHECK(result.error == MultiSourceLoadError::TruncatedExplicitPattern);
      BOOST_CHECK(result.source == ClusterSourceId{0});
      BOOST_CHECK_EQUAL(result.rof, 0u);
      BOOST_CHECK_EQUAL(result.clusterIndex, 0u);
      BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
    }
  }

  const std::vector<CompClusterExt> malformedClusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const std::array<unsigned char, 2> malformedPattern{0, 1};
  ClusterSourceInput malformedSource;
  malformedSource.id = ClusterSourceId{0};
  malformedSource.detector = o2::detectors::DetID::ITS;
  malformedSource.clusters = malformedClusters;
  malformedSource.patterns = malformedPattern;
  malformedSource.rofs = rofs;
  malformedSource.dictionary = &dict();
  malformedSource.layerToSurface = itsLayerToSurface;
  malformedSource.timing = ROFTimingConfig{40, 0, 0, 0};
  malformedSource.decoder = &decoder;
  MultiSourceFrame frame;
  const auto malformed = loadSources(
    frame, layout.getView().getSurfaceCatalogView(),
    gsl::span<const ClusterSourceInput>(&malformedSource, 1), {0, 0});
  BOOST_CHECK(malformed.error == MultiSourceLoadError::MalformedExplicitPattern);
  BOOST_CHECK_EQUAL(malformed.rof, 0u);
  BOOST_CHECK_EQUAL(malformed.clusterIndex, 0u);
}

BOOST_AUTO_TEST_CASE(ExactPatternConsumptionSucceedsAndTrailingBytesAreRejected)
{
  const auto layout = makeCombinedLayout();
  PatternContractDecoder decoder;
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  auto makeSource = [&](gsl::span<const unsigned char> patterns) {
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

  const std::vector<unsigned char> exact{onePixelPattern.begin(), onePixelPattern.end()};
  auto exactSource = makeSource(exact);
  MultiSourceFrame frame;
  BOOST_REQUIRE(loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&exactSource, 1), {0, 0}).ok());

  const std::vector<unsigned char> trailing{1, 1, 0x80, 0xff};
  auto trailingSource = makeSource(trailing);
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&trailingSource, 1), {0, 0});
  BOOST_CHECK(result.error == MultiSourceLoadError::TrailingPatternData);
  BOOST_CHECK_EQUAL(result.rof, 1u);
  BOOST_CHECK_EQUAL(result.clusterIndex, 1u);
  // The failed replacement load is transactional.
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 1u);

  auto missingDictionarySource = makeSource(exact);
  missingDictionarySource.dictionary = nullptr;
  const auto missingDictionary = loadSources(
    frame, layout.getView().getSurfaceCatalogView(),
    gsl::span<const ClusterSourceInput>(&missingDictionarySource, 1), {0, 0});
  BOOST_CHECK(missingDictionary.error == MultiSourceLoadError::MissingDictionary);
  BOOST_CHECK(missingDictionary.source == ClusterSourceId{0});
  BOOST_CHECK_EQUAL(missingDictionary.rof, 0u);
  BOOST_CHECK_EQUAL(missingDictionary.clusterIndex, 0u);
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 1u);
}

BOOST_AUTO_TEST_CASE(MissingDictionaryIsTypedBeforeProductionGeometryDecode)
{
  ITSGeometryClusterDecoder decoder;
  const CompClusterExt cluster{1, 1, CompCluster::InvalidPatternID, 0};
  BoundedPatternCursor patterns{onePixelPattern};
  const std::array<SurfaceId, 1> mapping{SurfaceId{0}};
  const auto decoded = decoder.decode(cluster, patterns, nullptr, mapping, ClusterSourceId{0}, 0, 0, false);
  BOOST_CHECK(decoded.error == ClusterDecodeError::MissingDictionary);
  BOOST_CHECK_EQUAL(patterns.consumed(), 0u);
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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
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
    const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::NonDenseSourceIds);
  }
  {
    // Duplicate ids {0, 0}.
    std::array<ClusterSourceInput, 2> sources{makeSource(ClusterSourceId{0}, decoderA), makeSource(ClusterSourceId{0}, decoderB)};
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::DuplicateSourceId);
  }
  {
    // Explicitly invalid id.
    std::array<ClusterSourceInput, 1> sources{makeSource(ClusterSourceId::invalid(), decoderA)};
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
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
    const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidROFRange);
  }
  {
    // Overlapping ranges.
    const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 2}, ROFRecord{{40, 0}, 1, 1, 1}};
    auto src = makeSrc(rofs);
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidROFRange);
  }
  {
    // Leading gap: first ROF does not begin at cluster index 0.
    const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 1, 1}};
    auto src = makeSrc(rofs);
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidROFRange);
  }
  {
    // Internal gap: rof0 covers [0,1), rof1 covers [2,2) i.e. starts at 2
    // while only cluster index 1 is unreferenced in between (2 clusters
    // total, so this leaves cluster 1 outside any ROF).
    const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}, ROFRecord{{40, 0}, 1, 2, 0}};
    auto src = makeSrc(rofs);
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidROFRange);
  }
  {
    // Trailing cluster: the ROFs cover only the first cluster, leaving the
    // second cluster unreferenced by any ROF.
    const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
    auto src = makeSrc(rofs);
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidROFRange);
  }
  {
    // Clusters without ROFs: zero ROFs is only valid when clusters is also
    // empty, but this source has two clusters.
    const std::vector<ROFRecord> rofs{};
    auto src = makeSrc(rofs);
    MultiSourceFrame frame;
    const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidROFRange);
  }
}

BOOST_AUTO_TEST_CASE(ZeroROFsIsValidWithZeroClusters)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());
  const std::vector<CompClusterExt> clusters{};
  const std::vector<unsigned char> patterns{};
  const std::vector<ROFRecord> rofs{};
  FakeClusterDecoder decoder{o2::detectors::DetID::ITS, {0}, false};

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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_CHECK(result.ok());
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
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
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::DetectorSurfaceMismatch);
}

BOOST_AUTO_TEST_CASE(InconsistentDecoderMetadataIsRejected)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const std::array<Corruption, 5> corruptions{
    Corruption::WrongSurface, Corruption::WrongSource, Corruption::WrongIndex,
    Corruption::WrongSourceROF, Corruption::WrongSensorDetector};
  for (const auto corruption : corruptions) {
    FakeClusterDecoder decoder{o2::detectors::DetID::ITS, {0}, false, corruption};
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
    const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InconsistentDecoderMetadata);
  }
}

BOOST_AUTO_TEST_CASE(LayerMappedTrueWithUnsafeLayerIsRejected)
{
  // A decoder reporting layerMapped=true must never be trusted to also have
  // produced a `layer` that is safe to index into layerToSurface with:
  // negative, or equal to/beyond its size, must both be rejected the same
  // way as layerMapped=false, without ever indexing layerToSurface[layer].
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());
  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};

  const std::array<Corruption, 2> corruptions{
    Corruption::LayerMappedTrueWithNegativeLayer, Corruption::LayerMappedTrueWithLayerOutOfRange};
  for (const auto corruption : corruptions) {
    FakeClusterDecoder decoder{o2::detectors::DetID::ITS, {0}, false, corruption};
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
    const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.error == MultiSourceLoadError::InvalidLayerMapping);
    BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
  }
}

BOOST_AUTO_TEST_CASE(SurfaceKindMismatchIsRejected)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  // An ITS decoder that reports it decoded a disk measurement while mapping
  // onto a cylinder surface must be rejected without inferring geometry
  // from surface count.
  const std::vector<CompClusterExt> itsClusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto itsPatterns = makePatternBytes(itsClusters.size());
  const std::vector<ROFRecord> itsRofs{ROFRecord{{0, 0}, 0, 0, 1}};
  FakeClusterDecoder wrongKindDecoder{o2::detectors::DetID::ITS, {0}, false, Corruption::WrongKind};

  ClusterSourceInput src;
  src.id = ClusterSourceId{0};
  src.detector = o2::detectors::DetID::ITS;
  src.clusters = itsClusters;
  src.patterns = itsPatterns;
  src.rofs = itsRofs;
  src.dictionary = &dict();
  src.layerToSurface = itsLayerToSurface;
  src.timing = ROFTimingConfig{40, 0, 0, 0};
  src.decoder = &wrongKindDecoder;

  MultiSourceFrame frame;
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0});
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::SurfaceKindMismatch);
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
  BOOST_REQUIRE(loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&goodSrc, 1), {0, 0}).ok());
  const auto totalBefore = frame.getTotalMeasurements();
  BOOST_REQUIRE_EQUAL(totalBefore, 1u);

  // Now attempt an invalid load (duplicate ids) on the SAME frame.
  std::array<ClusterSourceInput, 2> badSources{goodSrc, goodSrc}; // both id==0
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(badSources), {0, 0});
  BOOST_REQUIRE(!result.ok());

  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), totalBefore);
  BOOST_REQUIRE_EQUAL(frame.getSources().size(), 1u);
  BOOST_CHECK(frame.getSources()[0].detector == o2::detectors::DetID::ITS);
}

BOOST_AUTO_TEST_CASE(FailedLoadAfterFirstSourceStagedLeavesNoPartialState)
{
  // Unlike FailedLoadLeavesNoPartialState (which fails during up-front
  // source-id validation, before any source is even staged), this exercises
  // failure during decode/validation of the SECOND source, after the first
  // source has already been fully decoded into local staging storage.
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());

  const std::vector<CompClusterExt> clusters{{1, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 1}};
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
  labels.addElement(0, o2::MCCompLabel{1, 0, 0});
  FakeClusterDecoder decoder{o2::detectors::DetID::ITS, {0}, false};

  ClusterSourceInput goodSrc;
  goodSrc.id = ClusterSourceId{0};
  goodSrc.detector = o2::detectors::DetID::ITS;
  goodSrc.clusters = clusters;
  goodSrc.patterns = patterns;
  goodSrc.rofs = rofs;
  goodSrc.dictionary = &dict();
  goodSrc.labels = &labels;
  goodSrc.layerToSurface = itsLayerToSurface;
  goodSrc.timing = ROFTimingConfig{40, 0, 0, 0};
  goodSrc.decoder = &decoder;

  MultiSourceFrame frame;
  BOOST_REQUIRE(loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&goodSrc, 1), {0, 0}).ok());

  // Baseline: content and view pointer identity before the failing call.
  const std::vector<SurfaceMeasurement> baselineMeasurements(
    frame.getSurfaceMeasurements(SurfaceId{0}).begin(), frame.getSurfaceMeasurements(SurfaceId{0}).end());
  const std::vector<ROFIntervalBC> baselineIntervals(
    frame.getSourceIntervals(ClusterSourceId{0}).begin(), frame.getSourceIntervals(ClusterSourceId{0}).end());
  const std::vector<SourceMetadata> baselineSources = frame.getSources();
  const auto baselineLabel = frame.getLabels(ClusterRef{ClusterSourceId{0}, 0});
  BOOST_REQUIRE_EQUAL(baselineLabel.size(), 1u);
  const auto baselineView = frame.getView();

  // Second source: dense/unique id (so id-level validation passes and the
  // decoder actually runs for source 0), but fails once ITS is asked to map
  // onto an MFT surface -- i.e. only after source 0 has already been staged.
  FakeClusterDecoder decoderA{o2::detectors::DetID::ITS, {0}, false};
  FakeClusterDecoder decoderB{o2::detectors::DetID::ITS, {0}, false};

  ClusterSourceInput srcA = goodSrc;
  srcA.decoder = &decoderA;

  ClusterSourceInput srcB;
  srcB.id = ClusterSourceId{1};
  srcB.detector = o2::detectors::DetID::ITS;
  srcB.clusters = clusters;
  srcB.patterns = patterns;
  srcB.rofs = rofs;
  srcB.dictionary = &dict();
  const std::array<SurfaceId, 1> wrongMapping{SurfaceId{2}}; // MFT surface for an ITS source
  srcB.layerToSurface = wrongMapping;
  srcB.timing = ROFTimingConfig{40, 0, 0, 0};
  srcB.decoder = &decoderB;

  std::array<ClusterSourceInput, 2> sources{srcA, srcB};
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0});
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error == MultiSourceLoadError::DetectorSurfaceMismatch);
  BOOST_CHECK(result.source == ClusterSourceId{1});

  // The frame must be exactly as it was: same view pointers (proving the
  // owning vectors were never touched, not just left with equal content).
  // `surfaces` is the per-surface pointer/count span cache; its own address
  // only changes when assignLoadedData()/clear() actually reassigns the
  // underlying per-surface storage, so pointer identity here is exactly as
  // strong a proof as the old flattened `measurements` pointer used to be.
  const auto afterView = frame.getView();
  BOOST_CHECK(afterView.surfaces == baselineView.surfaces);
  BOOST_CHECK(afterView.rofIntervals == baselineView.rofIntervals);
  BOOST_CHECK(afterView.sourceROFOffsets == baselineView.sourceROFOffsets);
  BOOST_CHECK_EQUAL(afterView.nSurfaces, baselineView.nSurfaces);
  BOOST_CHECK_EQUAL(afterView.nSources, baselineView.nSources);

  // ...and identical measurements, identities, timing intervals, source
  // metadata, and label lookup.
  const auto afterMeasurements = frame.getSurfaceMeasurements(SurfaceId{0});
  BOOST_REQUIRE_EQUAL(afterMeasurements.size(), baselineMeasurements.size());
  for (size_t i = 0; i < afterMeasurements.size(); ++i) {
    BOOST_CHECK(afterMeasurements[i].cluster == baselineMeasurements[i].cluster);
    BOOST_CHECK(afterMeasurements[i].surface == baselineMeasurements[i].surface);
    BOOST_CHECK_EQUAL(afterMeasurements[i].sourceROF, baselineMeasurements[i].sourceROF);
  }
  const auto afterIntervals = frame.getSourceIntervals(ClusterSourceId{0});
  BOOST_REQUIRE_EQUAL(afterIntervals.size(), baselineIntervals.size());
  for (size_t i = 0; i < afterIntervals.size(); ++i) {
    BOOST_CHECK_EQUAL(afterIntervals[i].begin, baselineIntervals[i].begin);
    BOOST_CHECK_EQUAL(afterIntervals[i].end, baselineIntervals[i].end);
  }
  BOOST_REQUIRE_EQUAL(frame.getSources().size(), baselineSources.size());
  BOOST_CHECK(frame.getSources()[0].detector == baselineSources[0].detector);
  BOOST_CHECK_EQUAL(frame.getSources()[0].nROFs, baselineSources[0].nROFs);
  const auto afterLabel = frame.getLabels(ClusterRef{ClusterSourceId{0}, 0});
  BOOST_REQUIRE_EQUAL(afterLabel.size(), 1u);
  BOOST_CHECK(afterLabel[0] == baselineLabel[0]);
  // The second source's label container must not have been staged either.
  BOOST_CHECK(frame.getLabels(ClusterRef{ClusterSourceId{1}, 0}).empty());
}

BOOST_AUTO_TEST_CASE(EmptyFrameAccessorsAvoidNullPointerArithmetic)
{
  MultiSourceFrame frame;

  BOOST_CHECK(frame.getSurfaceMeasurements(SurfaceId{0}).empty());
  BOOST_CHECK(frame.getSourceIntervals(ClusterSourceId{0}).empty());
  BOOST_CHECK(frame.getLabels(ClusterRef{ClusterSourceId{0}, 0}).empty());

  // Loading zero sources into a layout with surfaces is legal and must
  // leave every per-surface bucket and per-source interval list empty.
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.valid());
  const auto result = loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>{}, {0, 0});
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);

  const auto view = frame.getView();
  BOOST_CHECK_EQUAL(view.nSources, 0u);
  BOOST_CHECK_EQUAL(view.getSurfaceMeasurementCount(SurfaceId{0}), 0u);
  BOOST_CHECK(view.getSurfaceMeasurements(SurfaceId{0}) == nullptr);
  BOOST_CHECK(view.getMeasurement(SurfaceId{0}, SurfaceMeasurementIndex{0}) == nullptr);
  BOOST_CHECK(frame.getSurfaceMeasurements(SurfaceId{0}).empty());
  BOOST_CHECK(frame.getMeasurement(SurfaceId{0}, SurfaceMeasurementIndex{0}) == nullptr);
}

BOOST_AUTO_TEST_CASE(EmptyLayoutWithZeroSourcesLoadsSuccessfully)
{
  // A layout with no surfaces at all, combined with zero sources, is the
  // most degenerate legal input: nothing to validate, nothing to decode,
  // nothing to commit.
  SurfaceGraph emptyTopology{0};
  emptyTopology.finalize();
  SurfaceGraph emptyLayout{{}, std::move(emptyTopology)};
  BOOST_REQUIRE(emptyLayout.valid());
  BOOST_CHECK_EQUAL(emptyLayout.getView().nSurfaces, 0u);

  MultiSourceFrame frame;
  const auto result = loadSources(frame, emptyLayout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>{}, {0, 0});
  BOOST_CHECK(result.ok());
  BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
  BOOST_CHECK_EQUAL(frame.getNSurfaces(), 0u);
  BOOST_CHECK_EQUAL(frame.getSources().size(), 0u);
}

BOOST_AUTO_TEST_CASE(ViewsAreStandardLayoutAndTriviallyCopyable)
{
  static_assert(std::is_standard_layout_v<MultiSourceFrameView>);
  static_assert(std::is_trivially_copyable_v<MultiSourceFrameView>);
  static_assert(std::is_standard_layout_v<SurfaceMeasurementSpan>);
  static_assert(std::is_trivially_copyable_v<SurfaceMeasurementSpan>);
  BOOST_CHECK(std::is_standard_layout_v<MultiSourceFrameView>);
  BOOST_CHECK(std::is_trivially_copyable_v<MultiSourceFrameView>);
}
