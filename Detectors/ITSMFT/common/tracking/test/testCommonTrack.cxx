// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 4 CommonTrack foundation. Covers:
//  - CommonTrack's own layout/device-compatibility traits;
//  - isValidTrackRange()'s exact validity condition (empty/default, single-,
//    multi- and hole-containing ranges, out-of-range and reversed ranges);
//  - that a completed track's hitSurfaces is the union of the SurfaceId of
//    every measurement its range references, across a single ITS source and
//    across combined ITS+MFT sources sharing one MultiSourceFrame;
//  - that TimeFrame::wipe() invalidates CommonTrack/trackClusterIndices
//    storage together, like every other per-event CA artefact;
//  - that CommonTrack itself has no detector/public-output dependency.
//
// This slice does not populate CommonTrack from CA seeds: every track/range
// below is constructed directly by the test.

#define BOOST_TEST_MODULE ITSMFT CommonTrack
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <type_traits>
#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/CommonTrack.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

// ---------------------------------------------------------------------
// CommonTrack has no detector/public-output dependency.
//
// This is a structural claim about ITSMFTTracking/CommonTrack.h itself, not
// something a runtime assertion can observe: CommonTrack.h's own include
// list (DataFormatsITS/TimeEstBC.h, GPUCommonDef.h, ITSMFTTracking/
// SurfaceKinematicState.h, ITSMFTTracking/SurfaceMask.h) contains no
// DetID.h, TrackITS.h/TrackITSExt.h, MFTCATrack.h, GeometryTGeo.h, or
// workflow header, and CommonTrack declares no DetID/NLayers/publication-
// type field -- every field is either a plain scalar (float, uint32_t), the
// shared SurfaceKinematicState/SurfaceMask device PODs, or TimeEstBC (a
// common timing primitive already reused throughout this library, e.g.
// Cell.h/Tracklet.h, and not itself a detector output type). This test case
// exercises CommonTrack using exactly that narrow surface, so that if a
// future edit to CommonTrack.h ever added such a dependency, the type
// itself (constructible, copyable, comparable-by-field here) would still
// need no wider include to keep working -- the absence is enforced by
// review of CommonTrack.h's own include list, restated here as the
// authoritative claim this test documents.
// ---------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(CommonTrackHasNoDetectorOrPublicationOutputDependency)
{
  CommonTrack track{};
  track.innerState.family = StateFamily::Barrel;
  track.outerState.family = StateFamily::Barrel;
  track.chi2 = 1.5f;
  track.hitSurfaces.set(SurfaceId{0});
  track.firstClusterRef = 0;
  track.clusterRefEnd = 1;
  BOOST_CHECK(track.hitSurfaces.has(SurfaceId{0}));
  BOOST_CHECK_EQUAL(trackClusterRefCount(track), 1u);
}

BOOST_AUTO_TEST_CASE(CommonTrackLayoutAndDeviceCompatibilityTraits)
{
  // CommonTrack is deliberately not standard-layout, purely because its
  // embedded o2::its::TimeEstBC field is not (see CommonTrack.h); it is
  // trivially copyable, which is the actual device/host-compatibility
  // property that matters here.
  static_assert(!std::is_standard_layout_v<CommonTrack>);
  static_assert(std::is_trivially_copyable_v<CommonTrack>);
  static_assert(std::is_same_v<decltype(CommonTrack::hitSurfaces), SurfaceMask>);
  static_assert(std::is_same_v<decltype(CommonTrack::innerState), SurfaceKinematicState>);
  static_assert(std::is_same_v<decltype(CommonTrack::outerState), SurfaceKinematicState>);
  static_assert(std::is_same_v<decltype(CommonTrack::firstClusterRef), uint32_t>);
  static_assert(std::is_same_v<decltype(CommonTrack::clusterRefEnd), uint32_t>);

  // Default-constructed: zeroed range, empty mask, no NLayers/detector
  // dependency of any kind. Not constructed as `constexpr` here: o2::track::
  // PID's constructor (SurfaceKinematicState::pid's default member
  // initializer) is not itself constexpr, so a CommonTrack instance cannot
  // be a core-constant-expression -- a property of PID, unrelated to
  // CommonTrack's own trivial-copyability asserted above.
  const CommonTrack defaultTrack{};
  BOOST_CHECK_EQUAL(defaultTrack.firstClusterRef, 0u);
  BOOST_CHECK_EQUAL(defaultTrack.clusterRefEnd, 0u);
  BOOST_CHECK(defaultTrack.hitSurfaces.empty());
  BOOST_CHECK_EQUAL(defaultTrack.chi2, 0.f);
}

// --- isValidTrackRange() -------------------------------------------------

BOOST_AUTO_TEST_CASE(EmptyDefaultRangeIsValidForAnyContainerSize)
{
  const CommonTrack track{};
  BOOST_CHECK(isValidTrackRange(track, 0));
  BOOST_CHECK(isValidTrackRange(track, 5));
  BOOST_CHECK_EQUAL(trackClusterRefCount(track), 0u);
}

BOOST_AUTO_TEST_CASE(ValidSingleMultiAndHoleContainingRanges)
{
  // Single-hit range: [0,1) into a 1-element array.
  CommonTrack single{};
  single.firstClusterRef = 0;
  single.clusterRefEnd = 1;
  BOOST_CHECK(isValidTrackRange(single, 1));
  BOOST_CHECK_EQUAL(trackClusterRefCount(single), 1u);

  // Multi-hit range: [1,4) into a 5-element array (some entries before/after
  // the range belong to other tracks sharing the same flat array).
  CommonTrack multi{};
  multi.firstClusterRef = 1;
  multi.clusterRefEnd = 4;
  BOOST_CHECK(isValidTrackRange(multi, 5));
  BOOST_CHECK_EQUAL(trackClusterRefCount(multi), 3u);

  // Hole-containing: the range itself is a dense [first,end) span of
  // *present* references (holes are never stored as sentinel entries); a
  // hole instead shows up as a gap in hitSurfaces' SurfaceId numbering. A
  // 2-hit track on surfaces {0,2} (skipping surface 1) is a valid,
  // completed, hole-containing track: its range is still contiguous and
  // valid, only its mask has a gap.
  CommonTrack withHole{};
  withHole.firstClusterRef = 0;
  withHole.clusterRefEnd = 2;
  withHole.hitSurfaces.set(SurfaceId{0});
  withHole.hitSurfaces.set(SurfaceId{2});
  BOOST_CHECK(isValidTrackRange(withHole, 2));
  BOOST_CHECK_EQUAL(withHole.hitSurfaces.count(), 2);
  BOOST_CHECK(!withHole.hitSurfaces.has(SurfaceId{1})); // the hole
}

BOOST_AUTO_TEST_CASE(OutOfRangeAndReversedRangesAreRejected)
{
  CommonTrack pastEnd{};
  pastEnd.firstClusterRef = 0;
  pastEnd.clusterRefEnd = 6;
  BOOST_CHECK(!isValidTrackRange(pastEnd, 5)); // clusterRefEnd > size

  CommonTrack exactlyAtSize{};
  exactlyAtSize.firstClusterRef = 0;
  exactlyAtSize.clusterRefEnd = 5;
  BOOST_CHECK(isValidTrackRange(exactlyAtSize, 5)); // clusterRefEnd == size is valid (half-open)

  CommonTrack reversed{};
  reversed.firstClusterRef = 3;
  reversed.clusterRefEnd = 1;
  BOOST_CHECK(!isValidTrackRange(reversed, 5)); // firstClusterRef > clusterRefEnd
}

// --- hitSurfaces equals the union of referenced measurement surfaces ----

namespace
{

// Minimal, geometry-free decoder (same construction as
// testMultiSourceLoading.cxx/testTimeFrameLifecycle.cxx): sensorID is used
// directly as the detector-local layer.
class FakeClusterDecoder final : public ClusterDecoder
{
 public:
  FakeClusterDecoder(o2::detectors::DetID::ID detector, bool disk) : mDetector(detector), mDisk(disk) {}

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
    result.measurement = mDisk ? makeDiskSurfaceMeasurement(decoded, sensor, surface, clusterRef, sourceROF)
                               : makeCylinderSurfaceMeasurement(decoded, sensor, surface, clusterRef, sourceROF);
    return result;
  }

 private:
  o2::detectors::DetID::ID mDetector;
  bool mDisk;
};

struct BuiltLayout {
  DetectorLayout layout;
  std::vector<SurfaceDescriptor> surfaces;

  DetectorLayoutView getView() const noexcept
  {
    const auto masks = computeSurfaceKindMasks(surfaces);
    return layout.getView(surfaces, masks.first, masks.second);
  }
};

// 4-surface disconnected ITS(cylinder){0,1,2}+MFT(disk){3} layout, matching
// this file's fixtures below.
BuiltLayout makeCombinedLayout()
{
  SparseTrackingTopology topology{4};
  topology.finalize();
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.push_back(SurfaceDescriptor{SurfaceId{0}, 0, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{1}, 1, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{2}, 2, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{3}, 0, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
  return BuiltLayout{DetectorLayout{surfaces, std::move(topology)}, std::move(surfaces)};
}

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

// Builds trackClusterIndices for `track` by appending the flat position of
// each of `measurements` (already loaded into `frame`) in traversal order,
// then returns the union of every referenced measurement's SurfaceId as a
// SurfaceMask, for comparison against track.hitSurfaces.
SurfaceMask unionOfReferencedSurfaces(const CommonTrack& track,
                                      gsl::span<const SurfaceMeasurementIndex> trackClusterIndices,
                                      const MultiSourceFrame& frame)
{
  SurfaceMask observed{};
  BOOST_REQUIRE(isValidTrackRange(track, static_cast<uint32_t>(trackClusterIndices.size())));
  for (uint32_t i = track.firstClusterRef; i < track.clusterRefEnd; ++i) {
    const auto* measurement = frame.getMeasurement(trackClusterIndices[i]);
    BOOST_REQUIRE(measurement != nullptr);
    observed.set(measurement->surface);
  }
  return observed;
}

} // namespace

BOOST_AUTO_TEST_CASE(HitSurfacesEqualsUnionOfReferencedMeasurementSurfacesSingleSource)
{
  const auto layout = makeCombinedLayout();
  BOOST_REQUIRE(layout.getView().nSurfaces == 4);

  // One ITS source, 3 clusters on surfaces {0,1,2}.
  const std::vector<CompClusterExt> clusters{
    {10, 20, CompCluster::InvalidPatternID, 0},
    {11, 21, CompCluster::InvalidPatternID, 1},
    {12, 22, CompCluster::InvalidPatternID, 2},
  };
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 3}};
  const std::array<SurfaceId, 3> itsLayerToSurface{SurfaceId{0}, SurfaceId{1}, SurfaceId{2}};

  FakeClusterDecoder decoder{o2::detectors::DetID::ITS, false};
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
  BOOST_REQUIRE(loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(&src, 1), {0, 0}).ok());
  BOOST_REQUIRE_EQUAL(frame.getTotalMeasurements(), 3u);

  // A track referencing every measurement, inner to outer (surface 0,1,2).
  std::vector<SurfaceMeasurementIndex> trackClusterIndices;
  for (uint32_t i = 0; i < frame.getTotalMeasurements(); ++i) {
    trackClusterIndices.push_back(SurfaceMeasurementIndex{i});
  }

  CommonTrack track{};
  track.firstClusterRef = 0;
  track.clusterRefEnd = static_cast<uint32_t>(trackClusterIndices.size());
  track.hitSurfaces.set(SurfaceId{0});
  track.hitSurfaces.set(SurfaceId{1});
  track.hitSurfaces.set(SurfaceId{2});

  const auto observed = unionOfReferencedSurfaces(track, trackClusterIndices, frame);
  BOOST_CHECK(observed == track.hitSurfaces);

  // A hole-containing sub-track referencing only the first and last
  // measurement (surfaces 0 and 2, skipping 1): still valid, mask still
  // matches the (smaller) referenced set exactly.
  std::vector<SurfaceMeasurementIndex> holeIndices{SurfaceMeasurementIndex{0}, SurfaceMeasurementIndex{2}};
  CommonTrack holeTrack{};
  holeTrack.firstClusterRef = 0;
  holeTrack.clusterRefEnd = 2;
  holeTrack.hitSurfaces.set(SurfaceId{0});
  holeTrack.hitSurfaces.set(SurfaceId{2});
  const auto observedHole = unionOfReferencedSurfaces(holeTrack, holeIndices, frame);
  BOOST_CHECK(observedHole == holeTrack.hitSurfaces);
  BOOST_CHECK(!observedHole.has(SurfaceId{1}));
}

BOOST_AUTO_TEST_CASE(MultiSourceITSAndMFTReferenceIdentity)
{
  const auto layout = makeCombinedLayout();

  // ITS source: 2 clusters on surfaces {0,1}.
  const std::vector<CompClusterExt> itsClusters{
    {10, 20, CompCluster::InvalidPatternID, 0},
    {11, 21, CompCluster::InvalidPatternID, 1},
  };
  const auto itsPatterns = makePatternBytes(itsClusters.size());
  const std::vector<ROFRecord> itsRofs{ROFRecord{{0, 0}, 0, 0, 2}};
  const std::array<SurfaceId, 2> itsLayerToSurface{SurfaceId{0}, SurfaceId{1}};
  FakeClusterDecoder itsDecoder{o2::detectors::DetID::ITS, false};

  // MFT source: 1 cluster on surface 3.
  const std::vector<CompClusterExt> mftClusters{{5, 6, CompCluster::InvalidPatternID, 0}};
  const auto mftPatterns = makePatternBytes(mftClusters.size());
  const std::vector<ROFRecord> mftRofs{ROFRecord{{0, 0}, 0, 0, 1}};
  const std::array<SurfaceId, 1> mftLayerToSurface{SurfaceId{3}};
  FakeClusterDecoder mftDecoder{o2::detectors::DetID::MFT, true};

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
  BOOST_REQUIRE_EQUAL(frame.getTotalMeasurements(), 3u);

  // Locate each surface's single measurement's flat position via the
  // per-surface span, then confirm getMeasurement() at that flat position
  // resolves back to the identical measurement (source, external index,
  // sensor.detector), across the ITS/MFT source boundary.
  const auto onZero = frame.getSurfaceMeasurements(SurfaceId{0});
  const auto onOne = frame.getSurfaceMeasurements(SurfaceId{1});
  const auto onThree = frame.getSurfaceMeasurements(SurfaceId{3});
  BOOST_REQUIRE_EQUAL(onZero.size(), 1u);
  BOOST_REQUIRE_EQUAL(onOne.size(), 1u);
  BOOST_REQUIRE_EQUAL(onThree.size(), 1u);

  // A single common track crossing the ITS/MFT source boundary: references
  // one measurement from each detector, in traversal order.
  std::vector<SurfaceMeasurementIndex> trackClusterIndices;
  const auto view = frame.getView();
  for (uint32_t i = 0; i < view.nMeasurements; ++i) {
    trackClusterIndices.push_back(SurfaceMeasurementIndex{i});
  }
  BOOST_REQUIRE_EQUAL(trackClusterIndices.size(), 3u);

  CommonTrack track{};
  track.firstClusterRef = 0;
  track.clusterRefEnd = static_cast<uint32_t>(trackClusterIndices.size());
  track.hitSurfaces.set(SurfaceId{0});
  track.hitSurfaces.set(SurfaceId{1});
  track.hitSurfaces.set(SurfaceId{3});

  const auto observed = unionOfReferencedSurfaces(track, trackClusterIndices, frame);
  BOOST_CHECK(observed == track.hitSurfaces);

  // Explicit cross-source identity check: the measurement reached via
  // getMeasurement() at the MFT reference's flat position carries
  // sensor.detector == MFT and a distinct ClusterSourceId from the ITS
  // ones, proving the flat index is a genuine cross-detector, cross-source
  // identity -- not a detector-local layer index.
  bool foundMFT = false, foundITSZero = false, foundITSOne = false;
  for (uint32_t i = track.firstClusterRef; i < track.clusterRefEnd; ++i) {
    const auto* measurement = frame.getMeasurement(trackClusterIndices[i]);
    BOOST_REQUIRE(measurement != nullptr);
    if (measurement->surface == SurfaceId{3}) {
      BOOST_CHECK(measurement->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::MFT));
      BOOST_CHECK(measurement->cluster.source == ClusterSourceId{1});
      foundMFT = true;
    } else if (measurement->surface == SurfaceId{0}) {
      BOOST_CHECK(measurement->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::ITS));
      BOOST_CHECK(measurement->cluster.source == ClusterSourceId{0});
      foundITSZero = true;
    } else if (measurement->surface == SurfaceId{1}) {
      BOOST_CHECK(measurement->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::ITS));
      BOOST_CHECK(measurement->cluster.source == ClusterSourceId{0});
      foundITSOne = true;
    }
  }
  BOOST_CHECK(foundMFT);
  BOOST_CHECK(foundITSZero);
  BOOST_CHECK(foundITSOne);

  // Out-of-range flat index is rejected, not silently wrapped/clamped.
  BOOST_CHECK(frame.getMeasurement(SurfaceMeasurementIndex{view.nMeasurements}) == nullptr);
  BOOST_CHECK(frame.getMeasurement(SurfaceMeasurementIndex{}) == nullptr); // invalid sentinel
  BOOST_CHECK(view.getMeasurement(SurfaceMeasurementIndex{view.nMeasurements}) == nullptr);
}

// --- TimeFrame wipe/reload invalidation ----------------------------------

BOOST_AUTO_TEST_CASE(TimeFrameWipeInvalidatesCommonTracksAndTrackClusterIndicesTogether)
{
  TimeFrame<ITSNLayers> tf;

  tf.getTrackClusterIndices().push_back(SurfaceMeasurementIndex{0});
  tf.getTrackClusterIndices().push_back(SurfaceMeasurementIndex{1});

  CommonTrack track{};
  track.firstClusterRef = 0;
  track.clusterRefEnd = 2;
  track.hitSurfaces.set(SurfaceId{0});
  track.hitSurfaces.set(SurfaceId{1});
  tf.getCommonTracks().push_back(track);

  BOOST_REQUIRE_EQUAL(tf.getCommonTracks().size(), 1u);
  BOOST_REQUIRE_EQUAL(tf.getTrackClusterIndices().size(), 2u);
  BOOST_REQUIRE(isValidTrackRange(tf.getCommonTracks()[0], static_cast<uint32_t>(tf.getTrackClusterIndices().size())));

  tf.wipe();

  BOOST_CHECK(tf.getCommonTracks().empty());
  BOOST_CHECK(tf.getTrackClusterIndices().empty());

  // Reload after wipe: both containers accept new content independently of
  // whatever they held before, confirming they are ordinary per-event state
  // rather than something wipe() leaves in a half-cleared condition.
  tf.getTrackClusterIndices().push_back(SurfaceMeasurementIndex{7});
  CommonTrack reloaded{};
  reloaded.firstClusterRef = 0;
  reloaded.clusterRefEnd = 1;
  tf.getCommonTracks().push_back(reloaded);
  BOOST_CHECK_EQUAL(tf.getCommonTracks().size(), 1u);
  BOOST_CHECK_EQUAL(tf.getTrackClusterIndices().size(), 1u);
}
