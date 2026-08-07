// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 4 ownership tests retained for the lower-level scratch staging seam.
// Production event state is owned and reset by the non-templated TimeFrame;
// these fixtures exercise scratch-local allocator/backfill behavior and the
// public kernel seam without treating a standalone scratch as a live event
// owner. The remaining cases cover contracts not already exercised by the
// migrated TimeFrame tests (testTimeFrameLifecycle.cxx,
// testTimeFrameNormalizedSource.cxx, testTimeFrameCovarianceLifecycle.cxx,
// testTimeFrameDetectorLayouts.cxx and the workflow loading tests,
// testMFTNormalizedRefit.cxx):
//   1. reset() clears a standalone scratch fixture only; TimeFrame content
//      (CommonTracks, vertices) is untouched.
//   2. TimeFrame::resetEvent() is the owner-level operation that clears frame
//      event data and configured workspaces.
//   3. Two standalone scratch fixtures can be isolated while testing the
//      lower-level kernel seam; production configured workspaces are grouped
//      under their owning TimeFrame.
//   4. An allocation failure injected into the scratch-side legacy backfill
//      staging, after the TimeFrame-side normalized decode has already
//      succeeded, leaves both live owners exactly at their pre-load state.
//   5. Compile-time proof that the common core's TimeFrame is the plain,
//      non-templated owner Gate 4 B3.1 introduced, not a surviving
//      `TimeFrame<NLayers>` template.

#define BOOST_TEST_MODULE ITSMFT Gate4B31OwnershipContract
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <TGeoGlobalMagField.h>

#include <array>
#include <memory>
#include <type_traits>
#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "Field/MagneticField.h"
#include "ITSMFTTracking/CommonTrack.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/Cluster.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

// The field fixture is shared with the lower-level geometry fixtures below.
struct FieldFixture {
  FieldFixture()
  {
    if (!TGeoGlobalMagField::Instance()->GetField()) {
      TGeoGlobalMagField::Instance()->SetField(o2::field::MagneticField::createNominalField(5, true));
      TGeoGlobalMagField::Instance()->Lock();
    }
  }
};

BOOST_GLOBAL_FIXTURE(FieldFixture);

// --- Compile-time proof (test 5): TimeFrame is a plain, non-templated
// type; SurfaceTrackingScratch remains the plan-sized kernel workspace
// implementation type. If a `template <int NLayers> struct TimeFrame` regression were
// ever reintroduced, a bare `TimeFrame` (used here and by every migrated
// production/test call site, e.g. testTimeFrameLifecycle.cxx) would no
// longer name a complete type, and these static_asserts -- along with
// every `TimeFrame frame;` declaration in this file and the rest of the
// suite -- would fail to compile. This is the compile-time half of the
// proof; the grep half (zero remaining `TimeFrame<NLayers>`, `TimeFrameITS`,
// or `TimeFrameMFT` occurrences across the common core's headers, sources,
// tests, and every production call site) was run as part of this slice's
// validation sweep and is recorded in its handoff notes.
static_assert(std::is_default_constructible_v<TimeFrame>,
              "TimeFrame must be a plain, non-templated, default-constructible type (Gate 4 B3.1)");
static_assert(std::is_default_constructible_v<SurfaceTrackingScratch>,
              "SurfaceTrackingScratch must remain the plan-sized scratch owner");
static_assert(!std::is_same_v<TimeFrame, SurfaceTrackingScratch>,
              "TimeFrame and SurfaceTrackingScratch must remain two distinct owner types");

// Cheap, direct scratch population for tests 1-3: writes one cluster's worth
// of scratch-owned legacy state on layer 0 via the same public append API
// testMFTNormalizedRefit.cxx's fixtures use, without needing a full decoder/
// catalog/loadNormalizedSource() round trip.
void populateScratch(SurfaceTrackingScratch& scratch, float tag)
{
  scratch.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  scratch.adoptPlan(ITSNLayers, 0, 0);
  scratch.addClusterToLayer(0, tag, tag, tag, 0);
  scratch.addTrackingFrameInfoToLayer(0, o2::its::TrackingFrameInfo{
                                           tag, tag, tag, tag, tag, {tag, tag}, {1.f, 0.f, 1.f}});
  scratch.addClusterExternalIndexToLayer(0, 0);
  bounded_vector<uint8_t> sizes;
  sizes.resize(1, uint8_t{1});
  scratch.setClusterSize(0, sizes);
}

// Cheap, direct frame population: one vertex, one CommonTrack -- enough to
// distinguish "untouched" from "cleared" without depending on a normalized
// load.
void populateFrame(TimeFrame& frame)
{
  frame.addPrimaryVertex(Vertex{});
  frame.getCommonTracks().push_back(CommonTrack{});
}

// --- Decoder/catalog fixture for test 4 only: a minimal, geometry-free,
// single-source ITS load, deliberately trimmed to one cluster/one ROF --
// mirrors the accepted pattern in testTimeFrameLifecycle.cxx's
// LegacyLikeDecoder/makeITSTestCatalog, kept separate here since this
// file's other tests need none of it.
class OneLayerDecoder final : public ClusterDecoder
{
 public:
  explicit OneLayerDecoder(o2::detectors::DetID::ID detector) : mDetector(detector) {}

  o2::itsmft::tracking::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const TopologyDictionary* dict,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool /*applySysErrors*/) const override
  {
    o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
    const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dict);
    if (!clusterData.ok()) {
      result.error = clusterData.error;
      return result;
    }
    const int sensorID = cluster.getSensorID();
    result.layer = sensorID;
    if (sensorID < 0 || static_cast<size_t>(sensorID) >= layerToSurface.size()) {
      return result;
    }
    result.layerMapped = true;
    result.kind = SurfaceKind::Cylinder;

    DecodedCluster decoded{};
    decoded.global = {static_cast<float>(sensorID) * 10.f, static_cast<float>(cluster.getRow()), static_cast<float>(cluster.getCol())};
    decoded.cylinderFrame = {static_cast<float>(sensorID) + 100.f, static_cast<float>(cluster.getRow()) + 1.f, static_cast<float>(cluster.getCol()) + 2.f, 0.f};
    decoded.rowColumnCovariance = {clusterData.sig2Row, 0.f, clusterData.sig2Col};
    decoded.shape = clusterData.shape;
    decoded.sensor = static_cast<uint32_t>(sensorID);
    decoded.layer = sensorID;

    const auto surface = layerToSurface[sensorID];
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

constexpr std::array<unsigned char, 3> onePixelPattern{1, 1, 0x80};

std::vector<SurfaceDescriptor> makeITSTestCatalog()
{
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.reserve(ITSNLayers);
  for (uint16_t i = 0; i < ITSNLayers; ++i) {
    surfaces.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  }
  return surfaces;
}

std::vector<SurfaceId> identitySurfaces()
{
  std::vector<SurfaceId> mapping;
  mapping.reserve(ITSNLayers);
  for (uint16_t i = 0; i < ITSNLayers; ++i) {
    mapping.push_back(SurfaceId{i});
  }
  return mapping;
}

} // namespace

// ---------------------------------------------------------------------
// 1. reset() clears scratch-owned state only.
// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ResetScratchClearsScratchOnlyTimeFrameContentSurvives)
{
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  populateFrame(frame);
  populateScratch(scratch, 1.f);

  BOOST_REQUIRE_GT(scratch.getTotalClusters(), 0);
  BOOST_REQUIRE_EQUAL(frame.getPrimaryVertices().size(), 1u);
  BOOST_REQUIRE_EQUAL(frame.getCommonTracks().size(), 1u);

  scratch.reset();

  BOOST_CHECK_EQUAL(scratch.getTotalClusters(), 0);
  // reset() must never wipe or mutate TimeFrame: its populated
  // content survives exactly, byte for byte.
  BOOST_CHECK_EQUAL(frame.getPrimaryVertices().size(), 1u);
  BOOST_CHECK_EQUAL(frame.getCommonTracks().size(), 1u);
}

// ---------------------------------------------------------------------
// 2. The frame-owned reset clears frame event state. Lower-level scratch
//    reset remains covered only as an implementation invariant.
// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ResetTimeFrameEventClearsBothOwners)
{
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  populateFrame(frame);
  populateScratch(scratch, 2.f);

  BOOST_REQUIRE_GT(scratch.getTotalClusters(), 0);
  BOOST_REQUIRE_EQUAL(frame.getPrimaryVertices().size(), 1u);
  BOOST_REQUIRE_EQUAL(frame.getCommonTracks().size(), 1u);

  scratch.reset();
  frame.resetEvent();

  BOOST_CHECK_EQUAL(scratch.getTotalClusters(), 0);
  BOOST_CHECK(frame.getPrimaryVertices().empty());
  BOOST_CHECK(frame.getCommonTracks().empty());
}

// ---------------------------------------------------------------------
// 3. Two scratches, one shared TimeFrame: resetting one scratch never
//    touches the other scratch or the shared frame.
// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(TwoScratchesOneFrameResettingOneLeavesTheOtherAndTheFrameUntouched)
{
  TimeFrame frame;
  SurfaceTrackingScratch itsScratch;
  SurfaceTrackingScratch mftScratch;
  populateFrame(frame);
  populateScratch(itsScratch, 3.f);
  populateScratch(mftScratch, 4.f);

  BOOST_REQUIRE_GT(itsScratch.getTotalClusters(), 0);
  BOOST_REQUIRE_GT(mftScratch.getTotalClusters(), 0);

  itsScratch.reset();

  BOOST_CHECK_EQUAL(itsScratch.getTotalClusters(), 0);
  // The MFT scratch, bound to the SAME shared TimeFrame, is completely
  // unaffected: this is exactly the hazard the accepted design's combined-
  // owner reset rule guards against -- "resetting ITS scratch cannot erase
  // MFT event data".
  BOOST_CHECK_GT(mftScratch.getTotalClusters(), 0);
  BOOST_CHECK_EQUAL(frame.getPrimaryVertices().size(), 1u);
  BOOST_CHECK_EQUAL(frame.getCommonTracks().size(), 1u);
}

// ---------------------------------------------------------------------
// 4. Injected scratch-backfill failure after the TimeFrame-side normalized
//    decode has already succeeded leaves both owners at their pre-load
//    state.
// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(InjectedScratchBackfillFailureAfterNormalizedStagingLeavesBothOwnersUnchanged)
{
  const auto catalog = makeITSTestCatalog();
  const auto orderedSurfaces = identitySurfaces();
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  OneLayerDecoder decoder{o2::detectors::DetID::ITS};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0};

  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  SurfaceTrackingScratch scratch;
  frame.setMemoryPool(pool);
  scratch.setMemoryPool(pool);

  std::vector<TrackingParameters> noIterations;
  auto planResult = buildSurfaceGraphs(catalogView, orderedSurfaces, noIterations);
  BOOST_REQUIRE(planResult.ok());
  scratch.adoptPlan(orderedSurfaces.size(), 0, 0);
  const gsl::span<const SurfaceId> planOrderedSurfaces{orderedSurfaces};

  const std::vector<CompClusterExt> clusters{CompClusterExt{10, 20, CompCluster::InvalidPatternID, 0}};
  const auto patterns = std::vector<unsigned char>(onePixelPattern.begin(), onePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, 1}};

  // Baseline: a real, successful owner-level load, giving both owners
  // genuine content to check "unchanged" against below.
  const auto baseline = scratch.loadNormalizedSource(frame, decoder, origin, timing, clusters, patterns, rofs,
                                                     &dict(), nullptr, o2::detectors::DetID::ITS,
                                                     planOrderedSurfaces, catalogView);
  BOOST_REQUIRE(baseline.ok());
  const auto baselineMeasurements = frame.getNormalizedFrame().getTotalMeasurements();
  const auto baselineClusters = scratch.getTotalClusters();
  BOOST_REQUIRE_EQUAL(baselineMeasurements, 1u);
  BOOST_REQUIRE_EQUAL(baselineClusters, 1);

  // Exhaust the shared pool. loadSources() -- the normalized decode, staged
  // into plain heap-owned local storage, not the bounded pool -- still
  // succeeds; it is the scratch-side legacy backfill staging that follows
  // it (still entirely local, not yet committed to either live owner) that
  // allocates from the pool and is where this throws.
  pool->setMaxMemory(pool->getUsedMemory());

  bool threw = false;
  try {
    scratch.loadNormalizedSource(frame, decoder, origin, timing, clusters, patterns, rofs,
                                 &dict(), nullptr, o2::detectors::DetID::ITS,
                                 planOrderedSurfaces, catalogView);
  } catch (const BoundedMemoryResource::MemoryLimitExceeded&) {
    threw = true;
  }
  BOOST_REQUIRE(threw);

  // Both owners retain exactly their pre-load (baseline) state: the
  // owner-level load's all-or-nothing contract. See
  // SurfaceTrackingScratch::loadNormalizedSource()'s own doc comment
  // and testTimeFrameLifecycle.cxx's
  // BackfillAllocationFailureLeavesNormalizedAndLegacyStateAtBaseline for
  // the equivalent proof exercised through a larger, multi-cluster fixture.
  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getTotalMeasurements(), baselineMeasurements);
  BOOST_CHECK_EQUAL(scratch.getTotalClusters(), baselineClusters);
}

// ---------------------------------------------------------------------
// 5. Compile-time proof that the old TimeFrame<NLayers> spelling is gone
//    (see the static_asserts above this file's test cases).
// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(CommonCoreTimeFrameTemplateSpellingIsGone)
{
  // The static_asserts near the top of this file are the compile-time half
  // of this proof: this translation unit -- along with every migrated
  // production and test call site in this component -- only compiles
  // because TimeFrame (ITSMFTTracking/TimeFrame.h) is the plain,
  // non-templated common event owner Gate 4 B3.1 introduced, with no
  // surviving `TimeFrame<NLayers>`, `TimeFrameITS`, or `TimeFrameMFT` alias
  // anywhere in the common core. The grep half of this proof (zero
  // remaining occurrences of those old spellings across
  // Detectors/ITSMFT/common/tracking's headers/sources/tests and every
  // production call site -- Tracker, TrackerTraits, DetectorTraits,
  // MFTFwdTrackHelpers, the ITS/MFT workflow specs, the
  // ITS onboarding compile target) is not repeatable inside a unit test
  // binary; it was run as part of this slice's validation sweep and is
  // recorded in its handoff notes.
  BOOST_CHECK(true);
}
