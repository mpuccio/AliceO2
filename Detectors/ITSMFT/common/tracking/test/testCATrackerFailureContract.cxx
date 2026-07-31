// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 3 common-CA failure contract: Tracker<NLayers>::clustersToTracks()
// exception classification, wipe-on-every-failure, and the exact drop
// sentinel.
//
// Contract under test (see CATracker.h/CATracker.cxx):
//  - TraversalException (structural/configuration failure): TimeFrame is
//    wiped, then the exception always rethrows, regardless of
//    DropTFUponFailure.
//  - BoundedMemoryResource::MemoryLimitExceeded and std::bad_alloc
//    (recoverable, per-TF resource failures): TimeFrame is wiped;
//    DropTFUponFailure=true returns the exact kDroppedTimeFrameResult
//    sentinel, DropTFUponFailure=false rethrows.
//  - Any other std::exception (e.g. std::runtime_error): treated as
//    unclassified/structural, wiped, always rethrows regardless of the
//    flag -- it must never be silently converted into a dropped-TF result.
//  - Valid empty input (a real layout/topology with zero loaded clusters)
//    completes without throwing and returns a non-negative, non-sentinel
//    result.
//  - A tracker instance that dropped one TimeFrame can immediately process a
//    following one successfully.
//
// The std::bad_alloc and unclassified-std::exception cases are exercised
// through InjectingTrackerTraits, a test-only TrackerTraits<ITSNLayers>
// subclass that overrides the virtual computeLayerTracklets() traversal-stage
// boundary to throw on demand, deterministically, without provoking real
// host OOM or needing to reach genuine tracklet/cell/road computation.
//
// Every fixture below establishes a real layout/plan (buildDetectorLayoutSet())
// and then loads a normalized source -- even the structural-failure cases,
// and even when that source carries zero clusters/ROFs -- before running
// tracking. This is load-bearing, not incidental: TimeFrame::initialise()
// unconditionally calls getNrof(layer) = mROFramesClusters[layer].size()-1
// on every layer, and a never-loaded (default-constructed, size-0)
// mROFramesClusters underflows that subtraction, corrupting memory deep
// inside prepareClusters() rather than throwing a clean exception.
// loadNormalizedSource() sizes mROFramesClusters[layer] to rofs.size()+1 for
// every layer regardless of whether clusters/rofs are empty, which is what
// makes that call, and every "iterate 0..getNrof()" loop reached afterward,
// safe. The structural-failure cases below produce their TraversalException
// through an invalid TrackingParameters/index-table configuration, not
// through a missing/stale plan: Gate 4 B2 Slice 2 removed the plan-currency
// concept entirely (initialiseTimeFrame() now takes the plan as an explicit
// `const DetectorLayoutSet&` parameter, so "no plan" is no longer a state a
// caller can even construct) -- see the removed
// StructuralFailureViaStaleLayoutAlwaysRethrowsAndWipes test's replacement
// note below for what covers the "always rethrows and wipes" contract now.
//
// The recoverable-failure fixtures trigger BoundedMemoryResource::
// MemoryLimitExceeded through TrackingParameters::MaxMemory (via
// Rig::forceMemoryLimitBelowCurrentUsage()), not by calling
// BoundedMemoryResource::setMaxMemory() directly on the pool: Tracker<N>::
// clustersToTracks() itself calls mMemoryPool->setMaxMemory(mTrkParams[
// iteration].MaxMemory) as the first statement in its try block on every
// call, which would silently undo a limit set directly on the pool object
// beforehand.

#define BOOST_TEST_MODULE ITSMFT CATracker failure contract
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

#include <gsl/gsl>
#include <oneapi/tbb/task_arena.h>

#include <TGeoGlobalMagField.h>
#include "Field/MagneticField.h"

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/CATracker.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/Constants.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

// Deterministic, geometry-free stand-in for GeometryClusterDecoder<DetId>,
// identical construction to testTimeFrameLifecycle.cxx /
// testTimeFrameNormalizedSource.cxx / testMultiSourceLoading.cxx.
class LegacyLikeDecoder final : public ClusterDecoder
{
 public:
  explicit LegacyLikeDecoder(o2::detectors::DetID::ID detector) : mDetector(detector) {}

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

// TrackerTraits::findRoads() unconditionally touches the global
// o2::base::Propagator singleton on first use, which in turn requires
// TGeoGlobalMagField to already hold a real o2::field::MagneticField
// object -- with none set (the state of every other test in this suite,
// none of which calls clustersToTracks() end to end), Propagator falls
// back to a legacy FairRunAna singleton that also does not exist in this
// process and segfaults dereferencing it. Only the tests that expect a
// genuinely successful clustersToTracks() run (valid empty input,
// continued processing after a drop) reach findRoads(); the
// structural/recoverable-failure tests throw/return before ever getting
// there and do not need this. A trivial default-constructed
// MagneticField (no field map file, zero solenoid current) is sufficient
// -- these tests never fit or propagate an actual trajectory since there
// are no clusters. TGeoGlobalMagField::Instance()->Lock() only allows one
// SetField() call per process, so this must run at most once.
void ensureTrivialMagneticFieldIsSet()
{
  static const bool done = [] {
    TGeoGlobalMagField::Instance()->SetField(new o2::field::MagneticField());
    TGeoGlobalMagField::Instance()->Lock();
    return true;
  }();
  (void)done;
}

constexpr std::array<unsigned char, 3> onePixelPattern{1, 1, 0x80};
constexpr std::array<unsigned char, 3> threePixelPattern{1, 3, 0xE0};

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
    surfaces.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
    // Matches o2::itsmft::resetDetectorDefaults(..., DetID::ITS)'s LayerxX0
    // default, so TrackerTraits::initialiseTimeFrame()'s LegacyMaterialMismatch
    // compatibility check passes for these unperturbed fixtures.
    const float xOverX0 = kNominalITSLayerX0[i];
    surfaces.back().material.xOverX0 = xOverX0;
    surfaces.back().material.arealDensityGPerCm2 = xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
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

struct Fixture {
  std::vector<CompClusterExt> clusters;
  std::vector<unsigned char> patterns;
  std::vector<ROFRecord> rofs;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
};

// 4 clusters on layers {0,1,0,2}, partitioned into 3 ROFs. Same shape as
// testTimeFrameLifecycle.cxx's fixture -- only needed to give the
// recoverable-failure fixture genuine per-event content to wipe.
Fixture makeFixture()
{
  Fixture f;
  f.clusters = {
    CompClusterExt{10, 20, CompCluster::InvalidPatternID, 0},
    CompClusterExt{11, 21, CompCluster::InvalidPatternID, 1},
    CompClusterExt{12, 22, CompCluster::InvalidPatternID, 0},
    CompClusterExt{13, 23, CompCluster::InvalidPatternID, 2},
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

std::vector<TrackingParameters> makeOneIterationITSParams(bool dropTFUponFailure, size_t maxMemory = std::numeric_limits<size_t>::max())
{
  std::vector<TrackingParameters> params(1);
  resetDetectorDefaults(params[0], o2::detectors::DetID::ITS);
  params[0].DropTFUponFailure = dropTFUponFailure;
  params[0].MaxMemory = maxMemory;
  return params;
}

// A valid FirstPass iteration 0 followed by a non-FirstPass (RebuildClusterLUT
// only, matching the legacy ITS async-iteration-3 shape) iteration 1, both ITS
// defaults -- callers mutate params[1]'s index-table fields to construct a
// deliberate mismatch against the configuration iteration 0 will commit.
std::vector<TrackingParameters> makeTwoIterationITSParams(bool dropTFUponFailure)
{
  std::vector<TrackingParameters> params(2);
  resetDetectorDefaults(params[0], o2::detectors::DetID::ITS);
  resetDetectorDefaults(params[1], o2::detectors::DetID::ITS);
  params[1].PassFlags = IterationSteps{IterationStep::RebuildClusterLUT};
  for (auto& p : params) {
    p.DropTFUponFailure = dropTFUponFailure;
  }
  return params;
}

// Deterministic injection seam for the std::bad_alloc and
// unclassified-std::exception cases: TrackerTraits::computeLayerTracklets()
// is virtual and is the first traversal stage Tracker<N>::clustersToTracks()
// calls after initialiseTimeFrame() succeeds (see CATracker.cxx's do/while
// loop), so overriding it to throw immediately exercises CATracker.cxx's
// catch chain deterministically -- without provoking real host OOM, and
// without ever reaching genuine tracklet/cell/road computation (so, unlike
// ValidEmptyInputCompletesWithoutErrorAndProducesNoTracks /
// TrackerRemainsUsableAfterADroppedTimeFrame, these cases do not need
// ensureTrivialMagneticFieldIsSet(): findRoads() is never reached).
enum class InjectedFailure { None,
                             BadAlloc,
                             UnclassifiedRuntimeError };

class InjectingTrackerTraits final : public TrackerTraits<ITSNLayers>
{
 public:
  InjectedFailure failure = InjectedFailure::None;

  void computeLayerTracklets(const int iteration, int iVertex) override
  {
    switch (failure) {
      case InjectedFailure::BadAlloc:
        throw std::bad_alloc{};
      case InjectedFailure::UnclassifiedRuntimeError:
        throw std::runtime_error{"injected unclassified failure"};
      case InjectedFailure::None:
        break;
    }
    TrackerTraits<ITSNLayers>::computeLayerTracklets(iteration, iVertex);
  }
};

// Bundles a TimeFrame<ITSNLayers>, a TraitsT/Tracker<ITSNLayers> pair, and a
// bounded memory pool -- the minimal wiring Tracker<N>::clustersToTracks()
// needs to run at all (task arena included: TrackerTraits::
// computeLayerTracklets() dereferences it unconditionally, even though the
// structural-failure cases never reach that far). TraitsT defaults to the
// real TrackerTraits<ITSNLayers>; RigT<InjectingTrackerTraits> (aliased
// ThrowingRig below) is used by the injected-failure tests.
template <class TraitsT = TrackerTraits<ITSNLayers>>
struct RigT {
  explicit RigT(bool dropTFUponFailure, size_t maxMemory = std::numeric_limits<size_t>::max())
    : pool(std::make_shared<BoundedMemoryResource>()),
      params(makeOneIterationITSParams(dropTFUponFailure, maxMemory)),
      tracker(&traits)
  {
    tf.setMemoryPool(pool);
    traits.setMemoryPool(pool);
    traits.setNThreads(1, arena);
    tracker.adoptTimeFrame(tf);
    tracker.setParameters(params);
    tracker.setMemoryPool(pool);
    tracker.setBz(0.5f);
  }

  std::shared_ptr<BoundedMemoryResource> pool;
  std::vector<TrackingParameters> params;
  TimeFrame<ITSNLayers> tf;
  TraitsT traits;
  Tracker<ITSNLayers> tracker;
  std::shared_ptr<tbb::task_arena> arena;
  // Must outlive `plan` (DetectorLayoutSet borrows a SurfaceCatalogView into
  // it, Gate 4 B2 Slice 2) -- declared before `plan` so it is constructed
  // first and destroyed last. Owned by this Rig, not by TimeFrame: wipe()
  // cannot touch either, by construction.
  std::vector<SurfaceDescriptor> catalog;
  std::optional<DetectorLayoutSet> plan;

  // Establishes a real, valid detector layout/plan + topology, without
  // loading any clusters, and binds it into the tracker (mirroring
  // ITSMFTTrackingInterface::runTracking()'s adoptDetectorLayoutSet() call).
  void establishValidLayout()
  {
    catalog = makeITSTestCatalog();
    const auto orderedSurfaces = identitySurfaces(ITSNLayers);
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    auto result = buildDetectorLayoutSet(catalogView, orderedSurfaces, TransitionPolicyTag::CylinderCylinder, params);
    BOOST_REQUIRE(result.ok());
    plan.emplace(std::move(*result.layout));
    tracker.adoptDetectorLayoutSet(*plan);
    tf.initTrackerTopologies(params);
  }

  // Loads clusters (or, with an empty Fixture, zero clusters -- still a
  // valid load that sizes every per-layer ROF boundary table to a real,
  // if trivial, state) through the same normalized-loading path production
  // code uses. This sizing is load-bearing: TimeFrame::initialise() calls
  // getNrof(layer) = mROFramesClusters[layer].size() - 1 unconditionally,
  // and a never-loaded (default-constructed, size-0) mROFramesClusters
  // underflows that subtraction, crashing deep inside prepareClusters()
  // before any failure-contract check ever runs. loadNormalizedSource()
  // sizes mROFramesClusters[layer] to rofs.size()+1 for every layer even
  // when rofs/clusters are empty, so calling it with an empty Fixture is
  // the only proven-safe way to reach a genuinely valid, still-empty
  // TimeFrame state.
  void loadSource(const Fixture& f)
  {
    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS};
    const o2::InteractionRecord origin{50, 5};
    const ROFTimingConfig timing{40, 0, 0, 0};
    const auto& orderedSurfaces = plan->getConfigurationKey().orderedSurfaces;
    const auto result = tf.loadNormalizedSource(decoder, origin, timing, f.clusters, f.patterns, f.rofs, &dict(),
                                                f.labels.getIndexedSize() > 0 ? &f.labels : nullptr, o2::detectors::DetID::ITS,
                                                gsl::span<const SurfaceId>{orderedSurfaces}, plan->getSurfaceCatalog());
    BOOST_REQUIRE(result.ok());

    // TrackerTraits::computeLayerTracklets() reads per-layer ROF counts
    // from mROFOverlapTableView (o2::its::LayerTiming), a separate table
    // from mROFramesClusters/getNrof() -- it is never populated by
    // loadNormalizedSource() and defaults to an unconfigured/garbage view.
    // A traversal that reaches computeLayerTracklets() without this being
    // set derives its ROF loop bound from that garbage view and walks out
    // of bounds. Mirrors TrackingInterface::configureROFLookupTables()'s
    // shape, but with every layer given the same trivial timing matching
    // this fixture's single combined ROF stream (real production input has
    // per-detector-param ROF length/delay/bias; none of that is exercised
    // by the failure-contract cases here, only the ROF *count* is load
    // -bearing).
    o2::its::LayerTiming timing2{};
    timing2.mNROFsTF = static_cast<unsigned int>(f.rofs.size());
    timing2.mROFLength = 40;
    typename TimeFrame<ITSNLayers>::ROFOverlapTableN rofTable;
    for (int iLayer = 0; iLayer < ITSNLayers; ++iLayer) {
      rofTable.defineLayer(iLayer, timing2);
    }
    rofTable.init();
    tf.setROFOverlapTable(rofTable);

    typename TimeFrame<ITSNLayers>::ROFMaskTableN mask{rofTable};
    mask.resetMask();
    for (int iLayer = 0; iLayer < ITSNLayers; ++iLayer) {
      mask.setROFsEnabled(iLayer, 0, timing2.mNROFsTF, 1);
    }
    tf.setMultiplicityCutMask(std::move(mask));
  }

  // Tracker<N>::clustersToTracks() calls mMemoryPool->setMaxMemory(params[0].
  // MaxMemory) as the very first statement in its try block, on every call
  // -- so a memory pool limit tightened directly on the pool object (rather
  // than through TrackingParameters::MaxMemory) is silently undone the
  // instant clustersToTracks() runs. Setting MaxMemory here below the
  // pool's current usage makes that same setMaxMemory() call throw
  // BoundedMemoryResource::MemoryLimitExceeded immediately, before
  // initialiseTimeFrame()/prepareClusters() or any traversal code runs.
  // Tracker<N> keeps its own copy of the parameters vector (setParameters()
  // copies by value), so it must be re-applied after this mutation.
  void forceMemoryLimitBelowCurrentUsage()
  {
    const auto used = pool->getUsedMemory();
    params[0].MaxMemory = used > 0 ? used - 1 : 0;
    tracker.setParameters(params);
  }

  void restoreUnboundedMemory()
  {
    params[0].MaxMemory = std::numeric_limits<size_t>::max();
    tracker.setParameters(params);
  }
};

using Rig = RigT<>;
using ThrowingRig = RigT<InjectingTrackerTraits>;

Fixture emptyFixture()
{
  return Fixture{};
}

} // namespace

// --- Sentinel exactness --------------------------------------------------

BOOST_AUTO_TEST_CASE(SentinelIsExactMatchNotSignCheck)
{
  BOOST_CHECK(isDroppedTimeFrame(kDroppedTimeFrameResult));
  BOOST_CHECK(!isDroppedTimeFrame(0.f));
  BOOST_CHECK(!isDroppedTimeFrame(3.5f));
  BOOST_CHECK(!isDroppedTimeFrame(-2.f));
  BOOST_CHECK(!isDroppedTimeFrame(std::numeric_limits<float>::infinity()));
  BOOST_CHECK(!isDroppedTimeFrame(std::numeric_limits<float>::quiet_NaN()));
}

// --- Structural failure: always rethrows, always wipes -------------------
//
// Gate 4 B2 Slice 2 removed this section's original mechanism
// (StructuralFailureViaStaleLayoutAlwaysRethrowsAndWipes: establish a valid
// layout, then TimeFrame::invalidateDetectorLayouts() right before running
// tracking to deterministically produce TraversalException{StaleLayout}).
// Neither invalidateDetectorLayouts() nor TraversalFailureReason::StaleLayout
// is reachable any more: initialiseTimeFrame() now takes the plan as an
// explicit `const DetectorLayoutSet&` parameter with no TimeFrame-owned
// currency concept to invalidate. The "TraversalException (structural/
// configuration failure): TimeFrame is wiped, then the exception always
// rethrows, regardless of DropTFUponFailure" contract this test protected is
// still covered below, through a different structural-failure reason
// (InvalidIndexTableConfigurationAlwaysRethrowsAndWipesRegardlessOfFlag /
// IndexTableConfigurationMismatchAlwaysRethrowsAndWipesRegardlessOfFlag): the
// contract under test is about TraversalException as a *category*, not about
// any one specific TraversalFailureReason value.

// --- Recoverable failure: DropTFUponFailure decides, always wipes --------

BOOST_AUTO_TEST_CASE(RecoverableFailureDroppedReturnsExactSentinelAndWipes)
{
  Rig rig{/*dropTFUponFailure=*/true};
  rig.establishValidLayout();
  rig.loadSource(makeFixture());
  BOOST_REQUIRE(rig.tf.getNormalizedFrame().getTotalMeasurements() > 0u);

  rig.forceMemoryLimitBelowCurrentUsage();

  const float result = rig.tracker.clustersToTracks();

  BOOST_CHECK(isDroppedTimeFrame(result));
  BOOST_CHECK_EQUAL(rig.tf.getNormalizedFrame().getTotalMeasurements(), 0u);
  BOOST_CHECK(rig.tf.getTracks().empty());
  // The plan lives on `rig`, not on TimeFrame (Gate 4 B2 Slice 2): wipe()
  // cannot touch it, by construction -- nothing left to assert here.
}

BOOST_AUTO_TEST_CASE(RecoverableFailureNotDroppedRethrowsButStillWipesFirst)
{
  Rig rig{/*dropTFUponFailure=*/false};
  rig.establishValidLayout();
  rig.loadSource(makeFixture());
  BOOST_REQUIRE(rig.tf.getNormalizedFrame().getTotalMeasurements() > 0u);

  rig.forceMemoryLimitBelowCurrentUsage();

  BOOST_CHECK_THROW(rig.tracker.clustersToTracks(), BoundedMemoryResource::MemoryLimitExceeded);

  // Wipe must have already happened before the exception propagated -- not
  // "the process is going down anyway".
  BOOST_CHECK_EQUAL(rig.tf.getNormalizedFrame().getTotalMeasurements(), 0u);
  BOOST_CHECK(rig.tf.getTracks().empty());
}

// --- std::bad_alloc: recoverable, same drop-or-rethrow policy ------------
//
// Injected via InjectingTrackerTraits rather than genuine host OOM (see the
// class comment above): computeLayerTracklets() throws std::bad_alloc
// immediately, after initialiseTimeFrame() has already established a real,
// valid state from establishValidLayout()+loadSource(makeFixture()).

BOOST_AUTO_TEST_CASE(BadAllocDroppedReturnsExactSentinelAndWipes)
{
  ThrowingRig rig{/*dropTFUponFailure=*/true};
  rig.establishValidLayout();
  rig.loadSource(makeFixture());
  BOOST_REQUIRE(rig.tf.getNormalizedFrame().getTotalMeasurements() > 0u);

  rig.traits.failure = InjectedFailure::BadAlloc;
  const float result = rig.tracker.clustersToTracks();

  BOOST_CHECK(isDroppedTimeFrame(result));
  BOOST_CHECK_EQUAL(rig.tf.getNormalizedFrame().getTotalMeasurements(), 0u);
  BOOST_CHECK(rig.tf.getTracks().empty());
}

BOOST_AUTO_TEST_CASE(BadAllocNotDroppedRethrowsButStillWipesFirst)
{
  ThrowingRig rig{/*dropTFUponFailure=*/false};
  rig.establishValidLayout();
  rig.loadSource(makeFixture());
  BOOST_REQUIRE(rig.tf.getNormalizedFrame().getTotalMeasurements() > 0u);

  rig.traits.failure = InjectedFailure::BadAlloc;
  BOOST_CHECK_THROW(rig.tracker.clustersToTracks(), std::bad_alloc);

  BOOST_CHECK_EQUAL(rig.tf.getNormalizedFrame().getTotalMeasurements(), 0u);
  BOOST_CHECK(rig.tf.getTracks().empty());
}

// --- Unclassified std::exception: always structural, never a sentinel ----
//
// A plain std::runtime_error (or any std::exception that is neither
// TraversalException, BoundedMemoryResource::MemoryLimitExceeded, nor
// std::bad_alloc) must always rethrow and never be silently converted into
// a dropped-TF result, regardless of DropTFUponFailure.

BOOST_AUTO_TEST_CASE(UnclassifiedExceptionAlwaysRethrowsAndWipesRegardlessOfFlag)
{
  for (const bool dropFlag : {false, true}) {
    ThrowingRig rig{dropFlag};
    rig.establishValidLayout();
    rig.loadSource(makeFixture());
    BOOST_REQUIRE(rig.tf.getNormalizedFrame().getTotalMeasurements() > 0u);

    rig.traits.failure = InjectedFailure::UnclassifiedRuntimeError;
    BOOST_CHECK_THROW(rig.tracker.clustersToTracks(), std::runtime_error);

    BOOST_CHECK_EQUAL(rig.tf.getNormalizedFrame().getTotalMeasurements(), 0u);
    BOOST_CHECK(rig.tf.getTracks().empty());
  }
}

// --- Index-table configuration failures: structural, always rethrow -------
//
// Both new TraversalFailureReason values (InvalidIndexTableConfiguration,
// IndexTableConfigurationMismatch; TrackerTraits.cxx::initialiseTimeFrame())
// are TraversalException, the same structural-failure category the removed
// StaleLayout test above used to cover -- so they must follow the identical
// always-rethrow-and-wipe contract, regardless of DropTFUponFailure.

BOOST_AUTO_TEST_CASE(InvalidIndexTableConfigurationAlwaysRethrowsAndWipesRegardlessOfFlag)
{
  for (const bool dropFlag : {false, true}) {
    Rig rig{dropFlag};
    rig.params[0].RowBins = 0; // structurally invalid
    rig.tracker.setParameters(rig.params);
    rig.establishValidLayout();
    rig.loadSource(emptyFixture());

    rig.tf.getTracks().push_back(CATrackType<ITSNLayers>{});
    BOOST_REQUIRE(!rig.tf.getTracks().empty());

    bool threw = false;
    try {
      rig.tracker.clustersToTracks();
    } catch (const TraversalException& e) {
      threw = true;
      BOOST_CHECK(e.getReason() == TraversalFailureReason::InvalidIndexTableConfiguration);
    }
    BOOST_CHECK(threw);

    BOOST_CHECK(rig.tf.getTracks().empty());
  }
}

BOOST_AUTO_TEST_CASE(IndexTableConfigurationMismatchAlwaysRethrowsAndWipesRegardlessOfFlag)
{
  ensureTrivialMagneticFieldIsSet(); // iteration 0 runs to completion (findRoads()) before iteration 1 fails
  for (const bool dropFlag : {false, true}) {
    Rig rig{dropFlag};
    rig.params = makeTwoIterationITSParams(dropFlag);
    rig.params[1].RowBins = rig.params[0].RowBins + 1; // mismatched against what iteration 0 will commit
    rig.tracker.setParameters(rig.params);
    rig.establishValidLayout();
    rig.loadSource(emptyFixture());

    bool threw = false;
    try {
      rig.tracker.clustersToTracks();
    } catch (const TraversalException& e) {
      threw = true;
      BOOST_CHECK(e.getReason() == TraversalFailureReason::IndexTableConfigurationMismatch);
    }
    BOOST_CHECK(threw);

    BOOST_CHECK(rig.tf.getTracks().empty());
  }
}

// --- Valid empty input -----------------------------------------------------

BOOST_AUTO_TEST_CASE(ValidEmptyInputCompletesWithoutErrorAndProducesNoTracks)
{
  ensureTrivialMagneticFieldIsSet();
  Rig rig{/*dropTFUponFailure=*/false};
  rig.establishValidLayout();
  rig.loadSource(emptyFixture());
  BOOST_REQUIRE_EQUAL(rig.tf.getNormalizedFrame().getTotalMeasurements(), 0u);

  float result = std::numeric_limits<float>::quiet_NaN();
  BOOST_CHECK_NO_THROW(result = rig.tracker.clustersToTracks());

  BOOST_CHECK(!isDroppedTimeFrame(result));
  BOOST_CHECK(result >= 0.f);
  BOOST_CHECK_EQUAL(rig.tf.getNumberOfTracks(), 0u);
}

// --- Continued processing after a drop ------------------------------------

BOOST_AUTO_TEST_CASE(TrackerRemainsUsableAfterADroppedTimeFrame)
{
  ensureTrivialMagneticFieldIsSet();
  Rig rig{/*dropTFUponFailure=*/true};
  rig.establishValidLayout();
  rig.loadSource(makeFixture());

  rig.forceMemoryLimitBelowCurrentUsage();
  const float dropped = rig.tracker.clustersToTracks();
  BOOST_REQUIRE(isDroppedTimeFrame(dropped));

  // Restore headroom and process a fresh (here, empty) TimeFrame on the
  // SAME Tracker/TrackerTraits instance -- proving the tracker/device stays
  // usable after a drop, matching the DPL device staying alive.
  rig.restoreUnboundedMemory();
  rig.loadSource(emptyFixture());

  float result = std::numeric_limits<float>::quiet_NaN();
  BOOST_CHECK_NO_THROW(result = rig.tracker.clustersToTracks());
  BOOST_CHECK(!isDroppedTimeFrame(result));
  BOOST_CHECK(result >= 0.f);
}
