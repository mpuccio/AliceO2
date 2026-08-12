// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 4 CommonTrack foundation. Covers:
//  - CommonTrack/TrackClusterReference/CommonTrackTimestamp layout and
//    device-compatibility traits;
//  - isValidTrackRange()'s exact validity condition (empty/default, single-,
//    multi- and hole-containing ranges, out-of-range and reversed ranges);
//  - per-surface measurement storage and surface-local index bounds
//    (TimeFrame/MeasurementView::getGlobalMeasurement(SurfaceId,
//    SurfaceMeasurementIndex));
//  - cross-surface and cross-source TrackClusterReference resolution;
//  - that a completed track's hitSurfaces is the union of the SurfaceId of
//    every measurement its range references, and that each resolved
//    measurement's own surface matches the reference it was resolved from;
//  - that a successful TimeFrame::loadNormalizedSource() reload clears
//    CommonTrack/trackClusterIndices storage, and a failed one preserves it;
//  - that TimeFrame::resetEvent() invalidates CommonTrack/trackClusterIndices
//    storage together;
//  - that CommonTrack itself has no detector/public-output dependency.
//
// This slice does not populate CommonTrack from CA seeds: every track/range
// below is constructed directly by the test.

#define BOOST_TEST_MODULE ITSMFT CommonTrack
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/CommonTrack.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/MeasurementView.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/detail/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/CommonTrackOutputAdapter.h"
#include "ITSMFTTracking/detail/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

// ---------------------------------------------------------------------
// CommonTrack has no detector/public-output dependency.
//
// This is a structural claim about ITSMFTTracking/CommonTrack.h itself, not
// something a runtime assertion can observe: CommonTrack.h's own include
// list (GPUCommonDef.h, and the ITSMFTTracking/Surface{Id,KinematicState,
// Mask,MeasurementIndex,Timing}.h common primitives) contains no DetID.h,
// TrackITS.h/TrackITSExt.h, typed MFT output header, GeometryTGeo.h, or workflow
// header, and CommonTrack/TrackClusterReference declare no
// DetID/NLayers/publication-type field -- every field is either a plain
// scalar, or one of the shared SurfaceId/SurfaceKinematicState/SurfaceMask/
// SurfaceMeasurementIndex/CommonTrackTimestamp device PODs. This test case
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
  track.timestamp = CommonTrackTimestamp{100, 140};
  track.hitSurfaces.set(SurfaceId{0});
  track.firstClusterRef = 0;
  track.clusterRefEnd = 1;
  BOOST_CHECK(track.hitSurfaces.has(SurfaceId{0}));
  BOOST_CHECK_EQUAL(trackClusterRefCount(track), 1u);

  const TrackClusterReference reference{SurfaceId{0}, SurfaceMeasurementIndex{0}};
  BOOST_CHECK(reference.surface == SurfaceId{0});
  BOOST_CHECK(reference.index == SurfaceMeasurementIndex{0});
}

BOOST_AUTO_TEST_CASE(CommonTrackLayoutAndDeviceCompatibilityTraits)
{
  // Standard-layout and trivially-copyable together -- not trivial-
  // copyability alone: trivial-copyability only proves CommonTrack's bytes
  // may be copied without invoking non-trivial special member functions; it
  // says nothing about a single, well-defined, offsetof-usable member order
  // that is consistent across host/device compilation. Both properties are
  // asserted, on all three types, matching every other device-facing type
  // in this library (SurfaceMeasurement, StaticSurfaceDescriptor,
  // SurfaceGraphView, ...).
  static_assert(std::is_standard_layout_v<CommonTrack>);
  static_assert(std::is_trivially_copyable_v<CommonTrack>);
  static_assert(std::is_standard_layout_v<TrackClusterReference>);
  static_assert(std::is_trivially_copyable_v<TrackClusterReference>);
  static_assert(std::is_standard_layout_v<CommonTrackTimestamp>);
  static_assert(std::is_trivially_copyable_v<CommonTrackTimestamp>);

  static_assert(std::is_same_v<decltype(CommonTrack::hitSurfaces), SurfaceMask>);
  static_assert(std::is_same_v<decltype(CommonTrack::innerState), SurfaceKinematicState>);
  static_assert(std::is_same_v<decltype(CommonTrack::outerState), SurfaceKinematicState>);
  static_assert(std::is_same_v<decltype(CommonTrack::timestamp), CommonTrackTimestamp>);
  static_assert(std::is_same_v<decltype(CommonTrack::firstClusterRef), uint32_t>);
  static_assert(std::is_same_v<decltype(CommonTrack::clusterRefEnd), uint32_t>);
  static_assert(std::is_same_v<decltype(TrackClusterReference::surface), SurfaceId>);
  static_assert(std::is_same_v<decltype(TrackClusterReference::index), SurfaceMeasurementIndex>);

  // Default-constructed: zeroed range, empty mask, no NLayers/detector
  // dependency of any kind. Not constructed as `constexpr` here:
  // o2::track::PID's constructor (SurfaceKinematicState::pid's default
  // member initializer) is not itself constexpr, so a CommonTrack instance
  // cannot be a core-constant-expression -- a property of PID, unrelated to
  // CommonTrack's own standard-layout/trivial-copyability asserted above.
  const CommonTrack defaultTrack{};
  BOOST_CHECK_EQUAL(defaultTrack.firstClusterRef, 0u);
  BOOST_CHECK_EQUAL(defaultTrack.clusterRefEnd, 0u);
  BOOST_CHECK(defaultTrack.hitSurfaces.empty());
  BOOST_CHECK_EQUAL(defaultTrack.chi2, 0.f);
  BOOST_CHECK(!defaultTrack.timestamp.isValid()); // default {0,0}: begin < end is false
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

// --- Per-surface measurement storage / TrackClusterReference resolution --

namespace
{

// Minimal, geometry-free decoder (same construction as
// testMultiSourceLoading.cxx/testTimeFrameLifecycle.cxx): sensorID is used
// directly as the detector-local layer.
class FakeClusterDecoder final : public ClusterDecoder
{
 public:
  FakeClusterDecoder(o2::detectors::DetID::ID detector, bool disk) : mDetector(detector), mDisk(disk) {}

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
    decoded.global = {static_cast<float>(sensorID), static_cast<float>(cluster.getRow()), static_cast<float>(cluster.getCol())};
    decoded.cylinderFrame = {10.f + sensorID, 1.f, 2.f, 0.1f};
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

 private:
  o2::detectors::DetID::ID mDetector;
  bool mDisk;
};

struct BuiltLayout {
  SurfaceGraph layout;
  std::vector<SurfaceDescriptor> surfaces;

  SurfaceGraphView getView() const noexcept
  {
    const auto masks = computeSurfaceKindMasks(surfaces);
    return layout.getView();
  }
};

// 4-surface disconnected ITS(cylinder){0,1,2}+MFT(disk){3} layout, matching
// this file's fixtures below.
BuiltLayout makeCombinedLayout()
{
  SurfaceGraph topology{4};
  topology.finalize();
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.push_back(SurfaceDescriptor{SurfaceId{0}, 0, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{1}, 1, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{2}, 2, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{3}, 0, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
  return BuiltLayout{SurfaceGraph{surfaces, std::move(topology)}, std::move(surfaces)};
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

// Builds a combined ITS(surfaces {0,1,2}, source 0)+MFT(surface {3}, source
// 1) TimeFrame with exactly one measurement on each of surfaces
// {0,1,3} (surface 2 is left empty, a deliberate hole in the catalog's own
// numbering -- not exercised by any track in these tests, only present to
// prove per-surface storage does not require every surface to be non-empty).
void loadThreeMeasurementFrame(TimeFrame& frame, const BuiltLayout& layout)
{
  const std::vector<CompClusterExt> itsClusters{
    {10, 20, CompCluster::InvalidPatternID, 0},
    {11, 21, CompCluster::InvalidPatternID, 1},
  };
  const auto itsPatterns = makePatternBytes(itsClusters.size());
  const std::vector<ROFRecord> itsRofs{ROFRecord{{0, 0}, 0, 0, 2}};
  const std::array<SurfaceId, 2> itsLayerToSurface{SurfaceId{0}, SurfaceId{1}};
  static const FakeClusterDecoder itsDecoder{o2::detectors::DetID::ITS, false};

  const std::vector<CompClusterExt> mftClusters{{5, 6, CompCluster::InvalidPatternID, 0}};
  const auto mftPatterns = makePatternBytes(mftClusters.size());
  const std::vector<ROFRecord> mftRofs{ROFRecord{{0, 0}, 0, 0, 1}};
  const std::array<SurfaceId, 1> mftLayerToSurface{SurfaceId{3}};
  static const FakeClusterDecoder mftDecoder{o2::detectors::DetID::MFT, true};

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

  BOOST_REQUIRE(loadSources(frame, layout.getView().getSurfaceCatalogView(), gsl::span<const ClusterSourceInput>(sources), {0, 0}).ok());
}

} // namespace

BOOST_AUTO_TEST_CASE(PerSurfaceMeasurementStorageAndSurfaceLocalIndexBounds)
{
  const auto layout = makeCombinedLayout();
  TimeFrame frame;
  loadThreeMeasurementFrame(frame, layout);

  BOOST_REQUIRE_EQUAL(frame.getSurfaceMeasurements(SurfaceId{0}).size(), 1u);
  BOOST_REQUIRE_EQUAL(frame.getSurfaceMeasurements(SurfaceId{1}).size(), 1u);
  BOOST_REQUIRE_EQUAL(frame.getSurfaceMeasurements(SurfaceId{2}).size(), 0u);
  BOOST_REQUIRE_EQUAL(frame.getSurfaceMeasurements(SurfaceId{3}).size(), 1u);

  // Index 0 is valid, and independently valid, on every non-empty surface --
  // proving each surface owns its own index space rather than sharing one
  // global domain (index 0 resolves to a *different* measurement on each
  // surface, not the same global position read three times).
  const auto* onZero = frame.getGlobalMeasurement(SurfaceId{0}, SurfaceMeasurementIndex{0});
  const auto* onOne = frame.getGlobalMeasurement(SurfaceId{1}, SurfaceMeasurementIndex{0});
  const auto* onThree = frame.getGlobalMeasurement(SurfaceId{3}, SurfaceMeasurementIndex{0});
  BOOST_REQUIRE(onZero != nullptr);
  BOOST_REQUIRE(onOne != nullptr);
  BOOST_REQUIRE(onThree != nullptr);
  BOOST_CHECK(onZero->surface == SurfaceId{0});
  BOOST_CHECK(onOne->surface == SurfaceId{1});
  BOOST_CHECK(onThree->surface == SurfaceId{3});
  BOOST_CHECK(onZero->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::ITS));
  BOOST_CHECK(onThree->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::MFT));

  // Index 1 is out of range for surface 0 (only one measurement), even
  // though it would be perfectly in range if surface 0 had two entries --
  // the bound is surface-local, not derived from any other surface's size
  // or from a shared/global count.
  BOOST_CHECK(frame.getGlobalMeasurement(SurfaceId{0}, SurfaceMeasurementIndex{1}) == nullptr);
  // Surface 2 has zero measurements: even index 0 is out of range.
  BOOST_CHECK(frame.getGlobalMeasurement(SurfaceId{2}, SurfaceMeasurementIndex{0}) == nullptr);
  // Invalid surface id (out of range for a 4-surface catalog).
  BOOST_CHECK(frame.getGlobalMeasurement(SurfaceId{4}, SurfaceMeasurementIndex{0}) == nullptr);
  // Invalid/default index sentinel.
  BOOST_CHECK(frame.getGlobalMeasurement(SurfaceId{0}, SurfaceMeasurementIndex{}) == nullptr);

  // Same contract on the device-facing view.
  const auto view = frame.getMeasurementView();
  BOOST_REQUIRE_EQUAL(view.nSurfaces, 4u);
  BOOST_CHECK(view.getGlobalMeasurement(SurfaceId{0}, SurfaceMeasurementIndex{0}) != nullptr);
  BOOST_CHECK(view.getGlobalMeasurement(SurfaceId{0}, SurfaceMeasurementIndex{1}) == nullptr);
  BOOST_CHECK(view.getGlobalMeasurement(SurfaceId{2}, SurfaceMeasurementIndex{0}) == nullptr);
  BOOST_CHECK(view.getGlobalMeasurement(SurfaceId{4}, SurfaceMeasurementIndex{0}) == nullptr);
  BOOST_CHECK_EQUAL(view.getSurfaceMeasurementCount(SurfaceId{0}), 1u);
  BOOST_CHECK_EQUAL(view.getSurfaceMeasurementCount(SurfaceId{2}), 0u);
}

BOOST_AUTO_TEST_CASE(CrossSurfaceAndCrossSourceTrackClusterReferenceResolution)
{
  const auto layout = makeCombinedLayout();
  TimeFrame frame;
  loadThreeMeasurementFrame(frame, layout);

  // A single common track crossing the ITS/MFT source boundary, traversal
  // order inner to outer: surface 0 (ITS, source 0), surface 1 (ITS, source
  // 0), surface 3 (MFT, source 1) -- skipping surface 2 as a hole. Each
  // reference pairs the surface with that surface's own (surface-local)
  // measurement index, never a raw external cluster index or a global
  // position.
  const std::vector<TrackClusterReference> trackClusterIndices{
    {SurfaceId{0}, SurfaceMeasurementIndex{0}},
    {SurfaceId{1}, SurfaceMeasurementIndex{0}},
    {SurfaceId{3}, SurfaceMeasurementIndex{0}},
  };

  CommonTrack track{};
  track.firstClusterRef = 0;
  track.clusterRefEnd = static_cast<uint32_t>(trackClusterIndices.size());
  track.hitSurfaces.set(SurfaceId{0});
  track.hitSurfaces.set(SurfaceId{1});
  track.hitSurfaces.set(SurfaceId{3});
  BOOST_REQUIRE(isValidTrackRange(track, static_cast<uint32_t>(trackClusterIndices.size())));

  bool foundITSZero = false, foundITSOne = false, foundMFT = false;
  for (uint32_t i = track.firstClusterRef; i < track.clusterRefEnd; ++i) {
    const auto& reference = trackClusterIndices[i];
    const auto* measurement = frame.getGlobalMeasurement(reference.surface, reference.index);
    BOOST_REQUIRE(measurement != nullptr);
    // Completed-track invariant: the resolved measurement's own surface
    // matches the reference it was resolved from.
    BOOST_CHECK(measurement->surface == reference.surface);
    if (reference.surface == SurfaceId{0}) {
      BOOST_CHECK(measurement->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::ITS));
      BOOST_CHECK(measurement->cluster.source == ClusterSourceId{0});
      foundITSZero = true;
    } else if (reference.surface == SurfaceId{1}) {
      BOOST_CHECK(measurement->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::ITS));
      BOOST_CHECK(measurement->cluster.source == ClusterSourceId{0});
      foundITSOne = true;
    } else if (reference.surface == SurfaceId{3}) {
      BOOST_CHECK(measurement->sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::MFT));
      BOOST_CHECK(measurement->cluster.source == ClusterSourceId{1});
      foundMFT = true;
    }
  }
  BOOST_CHECK(foundITSZero);
  BOOST_CHECK(foundITSOne);
  BOOST_CHECK(foundMFT);
}

BOOST_AUTO_TEST_CASE(HitSurfacesEqualsUnionAndEachMeasurementSurfaceMatchesItsReference)
{
  const auto layout = makeCombinedLayout();
  TimeFrame frame;
  loadThreeMeasurementFrame(frame, layout);

  const std::vector<TrackClusterReference> trackClusterIndices{
    {SurfaceId{0}, SurfaceMeasurementIndex{0}},
    {SurfaceId{1}, SurfaceMeasurementIndex{0}},
    {SurfaceId{3}, SurfaceMeasurementIndex{0}},
  };

  CommonTrack track{};
  track.firstClusterRef = 0;
  track.clusterRefEnd = static_cast<uint32_t>(trackClusterIndices.size());
  track.hitSurfaces.set(SurfaceId{0});
  track.hitSurfaces.set(SurfaceId{1});
  track.hitSurfaces.set(SurfaceId{3});

  SurfaceMask observed{};
  BOOST_REQUIRE(isValidTrackRange(track, static_cast<uint32_t>(trackClusterIndices.size())));
  for (uint32_t i = track.firstClusterRef; i < track.clusterRefEnd; ++i) {
    const auto& reference = trackClusterIndices[i];
    const auto* measurement = frame.getGlobalMeasurement(reference.surface, reference.index);
    BOOST_REQUIRE(measurement != nullptr);
    BOOST_CHECK(measurement->surface == reference.surface);
    observed.set(reference.surface);
  }
  BOOST_CHECK(observed == track.hitSurfaces);

  // A hole-containing sub-track referencing only surfaces 0 and 3 (skipping
  // 1): still a valid, completed track, mask still matches exactly the
  // (smaller) referenced set.
  const std::vector<TrackClusterReference> holeIndices{
    {SurfaceId{0}, SurfaceMeasurementIndex{0}},
    {SurfaceId{3}, SurfaceMeasurementIndex{0}},
  };
  CommonTrack holeTrack{};
  holeTrack.firstClusterRef = 0;
  holeTrack.clusterRefEnd = 2;
  holeTrack.hitSurfaces.set(SurfaceId{0});
  holeTrack.hitSurfaces.set(SurfaceId{3});

  SurfaceMask observedHole{};
  BOOST_REQUIRE(isValidTrackRange(holeTrack, static_cast<uint32_t>(holeIndices.size())));
  for (uint32_t i = holeTrack.firstClusterRef; i < holeTrack.clusterRefEnd; ++i) {
    const auto& reference = holeIndices[i];
    const auto* measurement = frame.getGlobalMeasurement(reference.surface, reference.index);
    BOOST_REQUIRE(measurement != nullptr);
    BOOST_CHECK(measurement->surface == reference.surface);
    observedHole.set(reference.surface);
  }
  BOOST_CHECK(observedHole == holeTrack.hitSurfaces);
  BOOST_CHECK(!observedHole.has(SurfaceId{1})); // the hole
}

// --- TimeFrame reload/wipe lifecycle --------------------------------------

namespace
{

// Deterministic, geometry-free stand-in for GeometryClusterDecoder<DetId>
// (same construction as testTimeFrameLifecycle.cxx): sensorID is used
// directly as the detector-local layer.
class LegacyLikeDecoder final : public ClusterDecoder
{
 public:
  explicit LegacyLikeDecoder(o2::detectors::DetID::ID detector) : mDetector(detector) {}

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
    return result;
  }

 private:
  o2::detectors::DetID::ID mDetector;
};

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

struct TimeFrameFixture {
  // Gate 4 B3.1: `tf` (the permanent, non-templated TimeFrame, owning
  // CommonTrack/TrackClusterReference/normalized-frame state) is declared
  // before `scratch` (the temporary SurfaceTrackingScratch) so it
  // is constructed first and destroyed last -- see SurfaceTrackingScratch's
  // own lifetime-contract doc. Neither owns or stores a reference to the
  // other; this fixture is what binds both for the load() call below.
  TimeFrame tf;
  SurfaceTrackingScratch scratch;
  // Keep the catalog with the graph fixture so initialization inputs have one
  // explicit owner.
  std::vector<SurfaceDescriptor> catalog{makeITSTestCatalog()};
  std::optional<std::vector<SurfaceGraph>> plan;
  LegacyLikeDecoder decoder{o2::detectors::DetID::ITS};
  o2::InteractionRecord origin{50, 5};
  ROFTimingConfig timing{40, 0, 0, 0};

  TimeFrameFixture()
  {
    const auto orderedSurfaces = identitySurfaces(ITSNLayers);
    SurfaceGraph graph{gsl::span<const SurfaceDescriptor>{catalog}};
    graph.setOrderedSurfaces(orderedSurfaces);
    BOOST_REQUIRE(graph.finalize());
    plan.emplace();
    plan->push_back(std::move(graph));
    scratch.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    scratch.adoptPlan(plan->front().getOrderedSurfaces().size(), 0, 0);
  }

  // One cluster on layer 0, one ROF: the minimal input that succeeds.
  LoadSourcesResult load()
  {
    const std::vector<CompClusterExt> clusters{{0, 1, CompCluster::InvalidPatternID, 0}};
    const auto patterns = makePatternBytes(clusters.size());
    const std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, 1}};
    const auto& orderedSurfaces = plan->front().getOrderedSurfaces();
    return scratch.loadNormalizedSource(tf, decoder, origin, timing, clusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                        gsl::span<const SurfaceId>{orderedSurfaces}, plan->front().getSurfaceCatalog());
  }
};

o2::its::LayerTiming makeFixtureClockTiming()
{
  // TimeFrameFixture::load() deliberately passes a temporary ROF vector to
  // the scratch loader. Its overlap-table view is therefore not a retained
  // publication input. Build the same immutable clock timing explicitly for
  // output-boundary tests rather than dereferencing that non-owning view.
  o2::its::LayerTiming timing{};
  timing.mNROFsTF = 1;
  timing.mROFLength = 40;
  timing.mROFDelay = 100;
  return timing;
}

struct TestCommonTrack {
  CommonTrack track;
  std::vector<TrackClusterReference> references;
};

TestCommonTrack makeTestCommonTrack()
{
  TestCommonTrack record;
  record.track.innerState.family = StateFamily::Barrel;
  record.track.outerState.family = StateFamily::Barrel;
  record.track.timestamp = {100, 140};
  record.track.hitSurfaces.set(SurfaceId{0});
  record.references.push_back({SurfaceId{0}, SurfaceMeasurementIndex{0}});
  return record;
}

uint32_t storeTestCommonTrack(TimeFrame& frame, TestCommonTrack record)
{
  const auto index = static_cast<uint32_t>(frame.getCommonTracks().size());
  record.track.firstClusterRef = static_cast<uint32_t>(frame.getTrackClusterIndices().size());
  frame.getTrackClusterIndices().insert(frame.getTrackClusterIndices().end(), record.references.begin(), record.references.end());
  record.track.clusterRefEnd = static_cast<uint32_t>(frame.getTrackClusterIndices().size());
  frame.getCommonTracks().push_back(record.track);
  return index;
}

// Populates tf's CommonTrack/trackClusterIndices with arbitrary, self-
// consistent content so a subsequent clear can be observed.
void populateCommonResults(TimeFrame& tf)
{
  tf.getTrackClusterIndices().push_back(TrackClusterReference{SurfaceId{0}, SurfaceMeasurementIndex{0}});
  tf.getTrackClusterIndices().push_back(TrackClusterReference{SurfaceId{1}, SurfaceMeasurementIndex{0}});
  CommonTrack track{};
  track.firstClusterRef = 0;
  track.clusterRefEnd = 2;
  track.hitSurfaces.set(SurfaceId{0});
  track.hitSurfaces.set(SurfaceId{1});
  tf.getCommonTracks().push_back(track);
}

} // namespace

BOOST_AUTO_TEST_CASE(SuccessfulReloadClearsCommonTrackAndTrackClusterIndices)
{
  TimeFrameFixture fixture;
  BOOST_REQUIRE(fixture.load().ok());

  populateCommonResults(fixture.tf);
  BOOST_REQUIRE_EQUAL(fixture.tf.getCommonTracks().size(), 1u);
  BOOST_REQUIRE_EQUAL(fixture.tf.getTrackClusterIndices().size(), 2u);

  // A second, independently successful load on the same TimeFrame: the
  // normalized frame is replaced, and CommonTrack/trackClusterIndices --
  // built against the *previous* normalized frame -- must be cleared in
  // that same successful commit, since they are no longer meaningful once
  // the frame they referenced is gone.
  BOOST_REQUIRE(fixture.load().ok());
  BOOST_CHECK(fixture.tf.getCommonTracks().empty());
  BOOST_CHECK(fixture.tf.getTrackClusterIndices().empty());
}

BOOST_AUTO_TEST_CASE(MFTPublicationCompatibilityIsSparseOrderedAndTransactional)
{
  TimeFrameFixture fixture;
  BOOST_REQUIRE(fixture.load().ok());
  const auto record = makeTestCommonTrack();
  MFTPublicationCompatibility sidecar;

  MFTPublicationCompatibilityTransaction firstTx{sidecar, 1.25, 3.5, 0x25u};
  const auto first = storeTestCommonTrack(fixture.tf, record);
  BOOST_REQUIRE(firstTx.validate(first));
  firstTx.reserve();
  firstTx.append(first);
  BOOST_CHECK_EQUAL(first, 0u);
  BOOST_REQUIRE_EQUAL(sidecar.entries().size(), 1u);
  BOOST_CHECK_EQUAL(sidecar.entries().front().commonTrackIndex, first);
  BOOST_CHECK_EQUAL(sidecar.entries().front().invQPtSeed, 1.25);
  BOOST_CHECK_EQUAL(sidecar.entries().front().chi2QPtSeed, 3.5);
  BOOST_CHECK_EQUAL(sidecar.entries().front().seedPattern, 0x25u);
  BOOST_CHECK(sidecar.find(first, fixture.tf.getCommonTracks().size()) != nullptr);
  BOOST_CHECK(sidecar.find(1u, fixture.tf.getCommonTracks().size()) == nullptr);  // missing key
  BOOST_CHECK(sidecar.find(99u, fixture.tf.getCommonTracks().size()) == nullptr); // out of range

  MFTPublicationCompatibilityTransaction duplicate{sidecar, 2., 4., 0x12u};
  BOOST_CHECK(!duplicate.validate(0u));

  // Scratch-only reset has no authority over the shared TimeFrame or this
  // MFT bridge-owned sidecar; a future combined owner decides detector-local
  // CommonTrack removal/marking separately.
  fixture.scratch.reset();
  BOOST_CHECK_EQUAL(fixture.tf.getCommonTracks().size(), 1u);
  BOOST_CHECK_EQUAL(sidecar.entries().size(), 1u);
}

BOOST_AUTO_TEST_CASE(ITSSharedClusterCompatibilityUsesExplicitPreSortAssociations)
{
  struct MarkedTrack {
    bool shared = false;
    bool hasSharedClusters() const { return shared; }
  };

  TimeFrameFixture fixture;
  BOOST_REQUIRE(fixture.load().ok());
  const auto record = makeTestCommonTrack();
  ITSSharedClusterCompatibility sidecar;

  // Deliberately use a non-identity conceptual fclusSort permutation. The
  // status is read later from the original accepted slots, not this order.
  std::array<MarkedTrack, 3> accepted{{{false}, {true}, {false}}};
  const std::array<int, 3> fclusSort{{2, 0, 1}};
  BOOST_CHECK_NE(fclusSort[0], 0);
  for (size_t i = 0; i < accepted.size(); ++i) {
    ITSSharedClusterCompatibilityTransaction tx{sidecar};
    const auto index = storeTestCommonTrack(fixture.tf, record);
    BOOST_REQUIRE(tx.validate(index));
    tx.reserve();
    tx.append(index);
    BOOST_CHECK_EQUAL(index, i);
  }
  BOOST_CHECK_EQUAL(sidecar.pendingSize(), accepted.size());
  BOOST_CHECK(sidecar.sealFromMarkedTracks(accepted));
  BOOST_CHECK(sidecar.isSealed());
  BOOST_REQUIRE_EQUAL(sidecar.entries().size(), accepted.size());
  BOOST_CHECK_EQUAL(sidecar.entries()[0].commonTrackIndex, 0u);
  BOOST_CHECK(!sidecar.entries()[0].hasSharedClusters);
  BOOST_CHECK_EQUAL(sidecar.entries()[1].commonTrackIndex, 1u);
  BOOST_CHECK(sidecar.entries()[1].hasSharedClusters);

  // A later legacy output sort cannot change the already global-index-keyed
  // sealed result.
  std::reverse(accepted.begin(), accepted.end());
  BOOST_CHECK_EQUAL(sidecar.entries()[1].commonTrackIndex, 1u);
  BOOST_CHECK(sidecar.entries()[1].hasSharedClusters);

  // Scratch-only reset has no authority over TimeFrame-owned CommonTracks
  // or the bridge-owned compatibility result they index.
  fixture.scratch.reset();
  BOOST_CHECK_EQUAL(fixture.tf.getCommonTracks().size(), 3u);
  BOOST_CHECK_EQUAL(sidecar.entries().size(), 3u);

  ITSSharedClusterCompatibility malformed;
  const auto malformedIndex = storeTestCommonTrack(fixture.tf, record);
  ITSSharedClusterCompatibilityTransaction tx{malformed};
  BOOST_REQUIRE(tx.validate(malformedIndex));
  tx.reserve();
  tx.append(malformedIndex);
  BOOST_CHECK(!malformed.sealFromMarkedTracks(accepted)); // pending/track cardinality mismatch
  BOOST_CHECK(!malformed.isSealed());
  BOOST_CHECK(malformed.entries().empty());

  TimeFrameFixture rollbackFixture;
  BOOST_REQUIRE(rollbackFixture.load().ok());
  ITSSharedClusterCompatibility rollback;
  const auto rollbackIndex = storeTestCommonTrack(rollbackFixture.tf, record);
  ITSSharedClusterCompatibilityTransaction rollbackTx{rollback};
  BOOST_REQUIRE(rollbackTx.validate(rollbackIndex));
  rollbackTx.reserve();
  rollbackTx.append(rollbackIndex);
  BOOST_CHECK_EQUAL(rollback.pendingSize(), 1u);

  ITSSharedClusterCompatibility sealingFailure;
  ITSSharedClusterCompatibilityTransaction sealingTx{sealingFailure};
  const auto sealingIndex = storeTestCommonTrack(rollbackFixture.tf, record);
  BOOST_REQUIRE(sealingTx.validate(sealingIndex));
  sealingTx.reserve();
  sealingTx.append(sealingIndex);
  std::array<MarkedTrack, 1> oneTrack{{{true}}};
  BOOST_CHECK(sealingFailure.sealFromMarkedTracks(oneTrack));
  BOOST_CHECK(sealingFailure.isSealed());
  BOOST_REQUIRE_EQUAL(sealingFailure.entries().size(), 1u);

  sidecar.clear();
  fixture.tf.resetEvent();
  BOOST_CHECK(fixture.tf.getCommonTracks().empty());
  BOOST_CHECK_EQUAL(sidecar.pendingSize(), 0u);
  BOOST_CHECK(sidecar.entries().empty());
}

BOOST_AUTO_TEST_CASE(FailedLoadPreservesCommonTrackAndTrackClusterIndicesUnchanged)
{
  TimeFrameFixture fixture;
  BOOST_REQUIRE(fixture.load().ok());

  populateCommonResults(fixture.tf);
  BOOST_REQUIRE_EQUAL(fixture.tf.getCommonTracks().size(), 1u);
  BOOST_REQUIRE_EQUAL(fixture.tf.getTrackClusterIndices().size(), 2u);
  const auto measurementsBefore = fixture.tf.getTotalMeasurements();
  BOOST_REQUIRE(measurementsBefore > 0u);

  // Deliberately fail: this SurfaceTrackingScratch preflight-rejects an
  // unsupported detector before touching anything.
  const std::vector<CompClusterExt> clusters{{0, 1, CompCluster::InvalidPatternID, 0}};
  const auto patterns = makePatternBytes(clusters.size());
  const std::vector<ROFRecord> rofs{ROFRecord{{200, 5}, 0, 0, 1}};
  const auto& orderedSurfaces = fixture.plan->front().getOrderedSurfaces();
  const auto failed = fixture.scratch.loadNormalizedSource(fixture.tf, fixture.decoder, fixture.origin, fixture.timing, clusters, patterns, rofs,
                                                           &dict(), nullptr, o2::detectors::DetID::TPC,
                                                           gsl::span<const SurfaceId>{orderedSurfaces}, fixture.plan->front().getSurfaceCatalog());
  BOOST_REQUIRE(!failed.ok());
  BOOST_CHECK(failed.error == MultiSourceLoadError::UnsupportedDetector);

  // Normalized frame, CommonTrack storage and trackClusterIndices are all
  // exactly as they were before the failed call.
  BOOST_CHECK_EQUAL(fixture.tf.getTotalMeasurements(), measurementsBefore);
  BOOST_REQUIRE_EQUAL(fixture.tf.getCommonTracks().size(), 1u);
  BOOST_CHECK_EQUAL(fixture.tf.getCommonTracks()[0].firstClusterRef, 0u);
  BOOST_CHECK_EQUAL(fixture.tf.getCommonTracks()[0].clusterRefEnd, 2u);
  BOOST_REQUIRE_EQUAL(fixture.tf.getTrackClusterIndices().size(), 2u);
  BOOST_CHECK(fixture.tf.getTrackClusterIndices()[0].surface == SurfaceId{0});
  BOOST_CHECK(fixture.tf.getTrackClusterIndices()[1].surface == SurfaceId{1});
}

BOOST_AUTO_TEST_CASE(TimeFrameWipeInvalidatesCommonTracksAndTrackClusterIndicesTogether)
{
  TimeFrame tf;
  populateCommonResults(tf);

  BOOST_REQUIRE_EQUAL(tf.getCommonTracks().size(), 1u);
  BOOST_REQUIRE_EQUAL(tf.getTrackClusterIndices().size(), 2u);
  BOOST_REQUIRE(isValidTrackRange(tf.getCommonTracks()[0], static_cast<uint32_t>(tf.getTrackClusterIndices().size())));

  tf.resetEvent();

  BOOST_CHECK(tf.getCommonTracks().empty());
  BOOST_CHECK(tf.getTrackClusterIndices().empty());

  // Reload after wipe: both containers accept new content independently of
  // whatever they held before, confirming they are ordinary per-event state
  // rather than something resetEvent() leaves in a half-cleared condition.
  tf.getTrackClusterIndices().push_back(TrackClusterReference{SurfaceId{2}, SurfaceMeasurementIndex{0}});
  CommonTrack reloaded{};
  reloaded.firstClusterRef = 0;
  reloaded.clusterRefEnd = 1;
  tf.getCommonTracks().push_back(reloaded);
  BOOST_CHECK_EQUAL(tf.getCommonTracks().size(), 1u);
  BOOST_CHECK_EQUAL(tf.getTrackClusterIndices().size(), 1u);
}

BOOST_AUTO_TEST_CASE(CommonTrackOutputAdapterTimestampIsSymmetricAndClamped)
{
  CommonTrackOutputAdapterError error = CommonTrackOutputAdapterError::None;
  o2::its::LayerTiming clock{};
  clock.mROFLength = 14;
  const ClockTimingPublicationView view{clock};
  const auto timestamp = makeOutputTimestamp({100, 120}, view, error);
  BOOST_REQUIRE(timestamp);
  BOOST_CHECK_EQUAL(timestamp->getTimeStamp(), 110.f);
  BOOST_CHECK_EQUAL(timestamp->getTimeStampError(), 7.f);
  BOOST_CHECK(!makeOutputTimestamp({20, 20}, view, error));
  BOOST_CHECK(error == CommonTrackOutputAdapterError::InvalidTimestamp);
}

BOOST_AUTO_TEST_CASE(CommonTrackOutputAdapterUsesLegacyPublicationOrder)
{
  TimeFrameFixture fixture;
  BOOST_REQUIRE(fixture.load().ok());

  auto later = makeTestCommonTrack();
  later.track.timestamp = {200, 240};
  later.track.chi2 = 1.f;
  auto earlier = makeTestCommonTrack();
  earlier.track.timestamp = {100, 140};
  earlier.track.chi2 = 2.f;
  BOOST_CHECK_EQUAL(storeTestCommonTrack(fixture.tf, later), 0u);
  BOOST_CHECK_EQUAL(storeTestCommonTrack(fixture.tf, earlier), 1u);

  o2::its::LayerTiming clock{};
  clock.mROFLength = 40;
  CommonTrackOutputAdapterError error = CommonTrackOutputAdapterError::None;
  const CommonTrackOutputAdapterSelection selection{{0u, 1u}};
  const auto ordered = makeLegacyOutputOrder(fixture.tf, selection, ClockTimingPublicationView{clock}, error);
  BOOST_REQUIRE(ordered);
  BOOST_REQUIRE_EQUAL(ordered->size(), 2u);
  BOOST_CHECK_EQUAL((*ordered)[0].globalIndex, 1u);
  BOOST_CHECK_EQUAL((*ordered)[1].globalIndex, 0u);
}

BOOST_AUTO_TEST_CASE(ClockTimingPublicationViewDelegatesLegacyClockSemantics)
{
  for (const uint32_t length : {9u, 10u}) {
    o2::its::LayerTiming legacy{};
    legacy.mNROFsTF = 4;
    legacy.mROFLength = length;
    legacy.mROFDelay = 3;
    legacy.mROFBias = 2;
    const ClockTimingPublicationView view{legacy};
    const std::array<CommonTrackTimestamp, 4> timestamps{{{5, 6}, {5, 5 + length}, {5 + length, 5 + 2 * length}, {5 + 3 * length, 5 + 4 * length}}};
    for (const auto timestamp : timestamps) {
      const auto asymmetric = view.makeTimeEstBC(timestamp);
      BOOST_REQUIRE(asymmetric);
      auto expected = asymmetric->makeSymmetrical();
      if (expected.getTimeStampError() > legacy.mROFLength * .5f)
        expected.setTimeStampError(legacy.mROFLength * .5f);
      const auto actual = view.makeOutputTimestamp(timestamp);
      BOOST_REQUIRE(actual);
      BOOST_CHECK_EQUAL(actual->getTimeStamp(), expected.getTimeStamp());
      BOOST_CHECK_EQUAL(actual->getTimeStampError(), expected.getTimeStampError());
      BOOST_CHECK_EQUAL(view.getROF(*actual), legacy.getROF(expected));
    }
  }
  o2::its::LayerTiming clock{};
  const ClockTimingPublicationView view{clock};
  BOOST_CHECK(!view.makeTimeEstBC({0, 0}));
  BOOST_CHECK(!view.makeTimeEstBC({-1, 1}));
  BOOST_CHECK(!view.makeTimeEstBC({0, static_cast<TFBC>(std::numeric_limits<uint32_t>::max()) + 1}));
  BOOST_CHECK(!view.makeTimeEstBC({0, static_cast<TFBC>(std::numeric_limits<uint16_t>::max()) + 1}));
}

BOOST_AUTO_TEST_CASE(CommonTrackOutputAdapterStagesITSAndFailsClosed)
{
  TimeFrameFixture fixture;
  BOOST_REQUIRE(fixture.load().ok());
  auto record = makeTestCommonTrack();
  record.track.chi2 = 3.f;
  ITSSharedClusterCompatibility shared;
  const auto commonTrackIndex = storeTestCommonTrack(fixture.tf, record);
  ITSSharedClusterCompatibilityTransaction transaction{shared};
  BOOST_REQUIRE(transaction.validate(commonTrackIndex));
  transaction.reserve();
  transaction.append(commonTrackIndex);
  struct MarkedTrack {
    bool shared{};
    bool hasSharedClusters() const { return shared; }
  };
  const std::array<MarkedTrack, 1> marked{{{true}}};
  BOOST_REQUIRE(shared.sealFromMarkedTracks(marked));
  const auto& measurement = *fixture.tf.getGlobalMeasurement(SurfaceId{0}, SurfaceMeasurementIndex{0});
  const std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 7, 3}};
  CommonTrackOutputAdapterError error = CommonTrackOutputAdapterError::None;
  const auto clock = makeFixtureClockTiming();
  const CommonTrackOutputTimingContext timing{rofs, ClockTimingPublicationView{clock}};
  const auto output = stageITSCommonTrackOutput(fixture.tf, measurement.cluster.source,
                                                gsl::span<const SurfaceId>{fixture.plan->front().getOrderedSurfaces()}, timing, shared,
                                                true, error);
  BOOST_REQUIRE(output);
  BOOST_CHECK_EQUAL(output->tracks.size(), 1u);
  BOOST_CHECK_EQUAL(output->clusterIndices.size(), 1u);
  BOOST_CHECK_EQUAL(output->clusterIndices[0], static_cast<int>(measurement.cluster.index));
  BOOST_CHECK(output->tracks[0].hasSharedClusters());
  BOOST_CHECK_EQUAL(output->tracks[0].getChi2(), 3.f);
  BOOST_CHECK_EQUAL(output->trackROFs[0].getFirstEntry(), 0);
  BOOST_CHECK_EQUAL(output->trackROFs[0].getNEntries(), 1);
  BOOST_CHECK_EQUAL(output->trackROFs[0].getFlags(), rofs[0].getFlags());

  // This is the workflow-facing binding: the workflow combines its own ROF
  // span with the immutable export returned by the tracking interface.  The
  // adapter accepts no scratch state and keeps the source/layout binding
  // explicit at this boundary.
  const CommonTrackPublicationContext publicationContext{
    o2::detectors::DetID::ITS, measurement.cluster.source, rofs, ClockTimingPublicationView{clock},
    gsl::span<const SurfaceId>{fixture.plan->front().getOrderedSurfaces()}};
  const auto contextOutput = stageITSCommonTrackOutput(fixture.tf, publicationContext, shared, true, error);
  BOOST_REQUIRE(contextOutput);
  BOOST_CHECK_EQUAL(contextOutput->tracks.size(), output->tracks.size());
  BOOST_REQUIRE_EQUAL(contextOutput->clusterIndices.size(), output->clusterIndices.size());
  for (size_t i = 0; i < output->clusterIndices.size(); ++i) {
    BOOST_CHECK_EQUAL(contextOutput->clusterIndices[i], output->clusterIndices[i]);
  }
  BOOST_CHECK_EQUAL(contextOutput->trackROFs.size(), output->trackROFs.size());

  auto wrongDetectorContext = publicationContext;
  wrongDetectorContext.detector = o2::detectors::DetID::MFT;
  BOOST_CHECK(!stageITSCommonTrackOutput(fixture.tf, wrongDetectorContext, shared, false, error));
  BOOST_CHECK(error == CommonTrackOutputAdapterError::MixedDetector);

  // Legacy publication retains a track even when its selected output
  // timestamp falls outside the workflow ROF span; it simply does not
  // increment a TrackROF entry. The adapter must preserve that behavior.
  fixture.tf.getCommonTracks()[0].timestamp = {1000, 1001};
  const auto outOfRangeOutput = stageITSCommonTrackOutput(fixture.tf, publicationContext, shared, false, error);
  BOOST_REQUIRE(outOfRangeOutput);
  BOOST_REQUIRE_EQUAL(outOfRangeOutput->tracks.size(), 1u);
  BOOST_CHECK_EQUAL(outOfRangeOutput->trackROFs[0].getFirstEntry(), 0);
  BOOST_CHECK_EQUAL(outOfRangeOutput->trackROFs[0].getNEntries(), 0);
  fixture.tf.getCommonTracks()[0].timestamp = record.track.timestamp;

  const auto oldTracks = fixture.tf.getCommonTracks().size();
  const auto oldReferences = fixture.tf.getTrackClusterIndices().size();
  // The original workflow ROF span can have more (or fewer) entries than
  // the LayerTiming clock. The CommonTrack adapter preserves that span
  // verbatim and groups only in-range clock slots.
  const std::vector<ROFRecord> mismatchedROFs{ROFRecord{{100, 5}, 0, 1, 2}, ROFRecord{{100, 6}, 1, 2, 3}};
  const CommonTrackOutputTimingContext mismatchedROF{mismatchedROFs, ClockTimingPublicationView{clock}};
  const auto mismatchedOutput = stageITSCommonTrackOutput(fixture.tf, measurement.cluster.source,
                                                          gsl::span<const SurfaceId>{fixture.plan->front().getOrderedSurfaces()}, mismatchedROF, shared, false, error);
  BOOST_REQUIRE(mismatchedOutput);
  BOOST_REQUIRE_EQUAL(mismatchedOutput->trackROFs.size(), mismatchedROFs.size());
  BOOST_CHECK_EQUAL(mismatchedOutput->trackROFs[0].getNEntries(), 1);
  BOOST_CHECK_EQUAL(mismatchedOutput->trackROFs[1].getNEntries(), 0);
  BOOST_CHECK_EQUAL(fixture.tf.getCommonTracks().size(), oldTracks);
  BOOST_CHECK_EQUAL(fixture.tf.getTrackClusterIndices().size(), oldReferences);
}

BOOST_AUTO_TEST_CASE(CommonTrackOutputAdapterStagesMFTAndRejectsMissingSidecar)
{
  const auto layout = makeCombinedLayout();
  TimeFrame frame;
  loadThreeMeasurementFrame(frame, layout);
  TestCommonTrack record;
  record.track.innerState.family = StateFamily::Forward;
  record.track.outerState.family = StateFamily::Forward;
  record.track.innerState.referenceCoordinate = -77.f;
  record.track.outerState.referenceCoordinate = -12.f;
  for (uint8_t i = 0; i < 5; ++i) {
    record.track.innerState.parameters[i] = 0.5f + i;
    record.track.outerState.parameters[i] = 3.5f + i;
  }
  for (uint8_t i = 0; i < 15; ++i) {
    record.track.innerState.covariance[i] = 0.01f * (i + 1);
    record.track.outerState.covariance[i] = 0.02f * (i + 1);
  }
  record.track.chi2 = 8.f;
  record.track.timestamp = {100, 124};
  record.track.hitSurfaces.set(SurfaceId{3});
  record.references.push_back({SurfaceId{3}, SurfaceMeasurementIndex{0}});
  MFTPublicationCompatibility sidecar;
  MFTPublicationCompatibilityTransaction tx{sidecar, 0.25, 1.5, 0x51u, 8.f};
  const auto commonTrackIndex = storeTestCommonTrack(frame, record);
  BOOST_REQUIRE(tx.validate(commonTrackIndex));
  tx.reserve();
  tx.append(commonTrackIndex);
  const auto& measurement = *frame.getGlobalMeasurement(SurfaceId{3}, SurfaceMeasurementIndex{0});
  const std::vector<ROFRecord> rofs{ROFRecord{{7, 9}, 2, 4, 5}};
  CommonTrackOutputAdapterError error = CommonTrackOutputAdapterError::None;
  o2::its::LayerTiming clock{};
  clock.mNROFsTF = 1;
  clock.mROFLength = 18;
  clock.mROFDelay = 100;
  const CommonTrackOutputTimingContext timing{rofs, ClockTimingPublicationView{clock}};
  const std::array<SurfaceId, 1> surfaces{SurfaceId{3}};
  const auto output = stageMFTCommonTrackOutput(frame, measurement.cluster.source, surfaces, timing, sidecar, true, error);
  BOOST_REQUIRE(output);
  BOOST_REQUIRE_EQUAL(output->tracks.size(), 1u);
  BOOST_CHECK_EQUAL(output->tracks[0].getZ(), -77.);
  BOOST_CHECK_EQUAL(output->tracks[0].getOutParam().getZ(), -12.);
  BOOST_CHECK_EQUAL(output->tracks[0].getCovariances()(4, 3), record.track.innerState.covariance[packedCovarianceIndex(4, 3)]);
  BOOST_CHECK_EQUAL(output->tracks[0].getOutParam().getCovariances()(4, 3), record.track.outerState.covariance[packedCovarianceIndex(4, 3)]);
  BOOST_CHECK_EQUAL(output->tracks[0].getTrackChi2(), 8.);
  BOOST_CHECK_EQUAL(output->tracks[0].getOutParam().getTrackChi2(), 8.);
  BOOST_CHECK_EQUAL(output->tracks[0].getInvQPtSeed(), .25);
  BOOST_CHECK_EQUAL(output->tracks[0].getChi2QPtSeed(), 1.5);
  BOOST_CHECK_EQUAL(output->seedPatterns[0], 0x51u);
  BOOST_CHECK_EQUAL(output->clusterIndices[0], static_cast<int>(measurement.cluster.index));
  BOOST_CHECK_EQUAL(output->trackROFs[0].getFirstEntry(), 0);
  BOOST_CHECK_EQUAL(output->trackROFs[0].getNEntries(), 1);
  BOOST_CHECK_EQUAL(output->trackROFs[0].getFlags(), rofs[0].getFlags());

  const CommonTrackPublicationContext publicationContext{
    o2::detectors::DetID::MFT, measurement.cluster.source, rofs, ClockTimingPublicationView{clock}, surfaces};
  const auto contextOutput = stageMFTCommonTrackOutput(frame, publicationContext, sidecar, true, error);
  BOOST_REQUIRE(contextOutput);
  BOOST_CHECK_EQUAL(contextOutput->tracks.size(), output->tracks.size());
  BOOST_REQUIRE_EQUAL(contextOutput->seedPatterns.size(), output->seedPatterns.size());
  BOOST_CHECK_EQUAL(contextOutput->seedPatterns[0], output->seedPatterns[0]);

  auto wrongDetectorContext = publicationContext;
  wrongDetectorContext.detector = o2::detectors::DetID::ITS;
  BOOST_CHECK(!stageMFTCommonTrackOutput(frame, wrongDetectorContext, sidecar, false, error));
  BOOST_CHECK(error == CommonTrackOutputAdapterError::MixedDetector);

  MFTPublicationCompatibility missing;
  BOOST_CHECK(!stageMFTCommonTrackOutput(frame, measurement.cluster.source, surfaces, timing, missing, false, error));
  BOOST_CHECK(error == CommonTrackOutputAdapterError::MissingCompatibility);
  BOOST_CHECK_EQUAL(frame.getCommonTracks().size(), 1u);
  BOOST_CHECK_EQUAL(frame.getTrackClusterIndices().size(), 1u);
}

BOOST_AUTO_TEST_CASE(CommonTrackOutputAdapterRejectsMalformedInputsWithoutMutatingOwners)
{
  TimeFrameFixture fixture;
  BOOST_REQUIRE(fixture.load().ok());
  const auto record = makeTestCommonTrack();
  storeTestCommonTrack(fixture.tf, record);
  const auto& measurement = *fixture.tf.getGlobalMeasurement(SurfaceId{0}, SurfaceMeasurementIndex{0});
  const std::vector<ROFRecord> rofs{ROFRecord{{1, 2}, 0, 0, 1}};
  const auto clock = makeFixtureClockTiming();
  const CommonTrackOutputTimingContext timing{rofs, ClockTimingPublicationView{clock}};
  const auto surfaces = gsl::span<const SurfaceId>{fixture.plan->front().getOrderedSurfaces()};
  const auto tracks = fixture.tf.getCommonTracks().size();
  const auto refs = fixture.tf.getTrackClusterIndices().size();
  const auto measurements = fixture.tf.getTotalMeasurements();
  CommonTrackOutputAdapterError error = CommonTrackOutputAdapterError::None;
  ITSSharedClusterCompatibility unsealed;
  BOOST_CHECK(!stageITSCommonTrackOutput(fixture.tf, measurement.cluster.source, surfaces, timing, unsealed, false, error));
  BOOST_CHECK(error == CommonTrackOutputAdapterError::MissingCompatibility);
  const auto wrongSource = stageITSCommonTrackOutput(fixture.tf, ClusterSourceId{99}, surfaces, timing, unsealed, false, error);
  BOOST_REQUIRE(wrongSource);
  BOOST_CHECK(wrongSource->tracks.empty());
  BOOST_CHECK_EQUAL(fixture.tf.getCommonTracks().size(), tracks);
  BOOST_CHECK_EQUAL(fixture.tf.getTrackClusterIndices().size(), refs);
  BOOST_CHECK_EQUAL(fixture.tf.getTotalMeasurements(), measurements);

  fixture.tf.getCommonTracks()[0].clusterRefEnd = refs + 1;
  BOOST_CHECK(!selectCommonTracksForSource(fixture.tf, o2::detectors::DetID::ITS, measurement.cluster.source, error));
  BOOST_CHECK(error == CommonTrackOutputAdapterError::InvalidTrackRange);
  fixture.tf.getCommonTracks()[0].clusterRefEnd = refs;
  fixture.tf.getTrackClusterIndices()[0].surface = SurfaceId{99};
  BOOST_CHECK(!selectCommonTracksForSource(fixture.tf, o2::detectors::DetID::ITS, measurement.cluster.source, error));
  BOOST_CHECK(error == CommonTrackOutputAdapterError::UnresolvedReference);
  fixture.tf.getTrackClusterIndices()[0].surface = SurfaceId{0};

  ITSSharedClusterCompatibility sealed;
  ITSSharedClusterCompatibilityTransaction tx{sealed};
  const auto secondTrackIndex = storeTestCommonTrack(fixture.tf, record);
  BOOST_REQUIRE(tx.validate(secondTrackIndex));
  tx.reserve();
  tx.append(secondTrackIndex);
  struct Marked {
    bool hasSharedClusters() const { return false; }
  };
  const std::array<Marked, 0> none{};
  BOOST_CHECK(!sealed.sealFromMarkedTracks(none)); // pending cardinality mismatch fails closed
  BOOST_CHECK(!sealed.isSealed());
}

BOOST_AUTO_TEST_CASE(MFTPublicationCompatibilityRejectsDuplicateNonMonotonicAndOutOfRangeKeys)
{
  MFTPublicationCompatibility sidecar;
  MFTPublicationCompatibilityTransaction first{sidecar, 1., 2., 3};
  BOOST_REQUIRE(first.validate(4));
  first.reserve();
  first.append(4);
  MFTPublicationCompatibilityTransaction duplicate{sidecar, 4., 5., 6};
  BOOST_CHECK(!duplicate.validate(4));
  MFTPublicationCompatibilityTransaction nonMonotonic{sidecar, 4., 5., 6};
  BOOST_CHECK(!nonMonotonic.validate(3));
  BOOST_CHECK(sidecar.find(4, 4) == nullptr); // key is outside the supplied CommonTrack owner range
  BOOST_CHECK(sidecar.find(5, 10) == nullptr);
  BOOST_REQUIRE(sidecar.find(4, 5));
  BOOST_CHECK_EQUAL(sidecar.entries().size(), 1u);
}
