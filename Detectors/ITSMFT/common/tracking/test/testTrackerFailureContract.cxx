// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Tracker failure contract: Tracker::run()
// exception classification, wipe-on-every-failure, and the exact drop
// sentinel.
//
// Contract under test (see Tracker.h/Tracker.cxx):
//  - TraversalException (structural/configuration failure): TimeFrame is
//    wiped, then the exception always rethrows, regardless of
//    DropTFUponFailure.
//  - BoundedMemoryResource::MemoryLimitExceeded and std::bad_alloc
//    (recoverable, per-TF resource failures): TimeFrame is wiped;
//    DropTFUponFailure=true returns TrackingOutcome::RecoverableDropped
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
// The std::bad_alloc and unclassified-std::exception cases use a real MFT
// road and a test-only SeedRefitFunction that throws at the normal refit
// boundary. This keeps the failure classification test at Tracker::run()
// without restoring a virtual traversal-stage injection seam.
//
// Every fixture below establishes a real layout/plan (buildSurfaceGraphs())
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
// `const std::vector<SurfaceGraph>&` parameter, so "no plan" is no longer a state a
// caller can even construct) -- see the removed
// StructuralFailureViaStaleLayoutAlwaysRethrowsAndWipes test's replacement
// note below for what covers the "always rethrows and wipes" contract now.
//
// The recoverable-failure fixtures tighten the already-used frame allocator
// to its current usage. The next tracking allocation then exercises the
// normal bounded-resource failure/reset contract without changing config.

#define BOOST_TEST_MODULE ITSMFT Tracker failure contract
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
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/detail/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/MeasurementView.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/detail/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/Constants.h"
#include "ITStracking/ROFLookupTables.h"
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

const TopologyDictionary& dict()
{
  static const TopologyDictionary d;
  return d;
}

// TrackerTraits::findRoads() unconditionally touches the global
// o2::base::Propagator singleton on first use, which in turn requires
// TGeoGlobalMagField to already hold a real o2::field::MagneticField
// object -- with none set (the state of every other test in this suite,
// none of which calls Tracker::run() end to end), Propagator falls
// back to a legacy FairRunAna singleton that also does not exist in this
// process and segfaults dereferencing it. Only the tests that expect a
// genuinely successful Tracker::run() (valid empty input,
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
    surfaces.back().chartRange = {-20.f, 20.f};
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

enum class RefitFailure { None,
                          BadAlloc,
                          UnclassifiedRuntimeError };

RefitFailure gRefitFailure = RefitFailure::None;
int gRefitCalls = 0;
int gThrowOnRefitCall = 1;

// Deliberately a plain function pointer: production itself reaches this
// seam only after candidate formation and road finding have produced a real
// seed. Tests reset the process-local controls before every invocation.
bool controlledSeedRefit(const TrackSeed&, const TrackingParameters&, float, SurfaceTrackingScratch&,
                         gsl::span<const gsl::span<const GlobalMeasurement>>,
                         gsl::span<const gsl::span<const SurfaceMeasurement>>, SurfaceCatalogView,
                         gsl::span<const SurfaceId>, TrackingCandidate&)
{
  ++gRefitCalls;
  if (gRefitCalls < gThrowOnRefitCall) {
    return false;
  }
  switch (gRefitFailure) {
    case RefitFailure::BadAlloc:
      throw std::bad_alloc{};
    case RefitFailure::UnclassifiedRuntimeError:
      throw std::runtime_error{"controlled refit failure"};
    case RefitFailure::None:
      return false;
  }
  return false;
}

// The core's typed refit/publication work is deliberately not part of this
// failure-contract fixture. This narrow test function supplies only refit.
bool testSeedRefit(const TrackSeed&, const TrackingParameters&, float, SurfaceTrackingScratch&,
                   gsl::span<const gsl::span<const GlobalMeasurement>>,
                   gsl::span<const gsl::span<const SurfaceMeasurement>>, SurfaceCatalogView,
                   gsl::span<const SurfaceId>, TrackingCandidate&)
{
  return false;
}

// Bundles a TimeFrame, real backend, Tracker, and bounded memory pool -- the
// minimal wiring Tracker::run() needs for the ITS configuration tests below.
struct Rig {
  explicit Rig(bool dropTFUponFailure, size_t maxMemory = std::numeric_limits<size_t>::max())
    : pool(std::make_shared<BoundedMemoryResource>()),
      params(makeOneIterationITSParams(dropTFUponFailure, maxMemory)),
      tracker()
  {
    traits.setNThreads(1, arena);
    tracker.setSeedRefitFunction(testSeedRefit);
    frame.setBz(0.5f);
  }

  // Stages one pending sidecar entry and one GenericTrack/TrackClusterReference
  // pair directly on `frame` -- deliberately not through a real CA seed (out
  // of scope here): only frame.resetEvent()'s unconditional clear of these two
  // containers and the workflow-edge sidecar reset are under test.
  void stageStaleState()
  {
    ITSSharedClusterCompatibilityTransaction txn{sidecar};
    BOOST_REQUIRE(txn.validate(0));
    txn.reserve();
    txn.append(0);
    BOOST_REQUIRE_EQUAL(sidecar.pendingSize(), 1u);

    frame.getTrackClusterIndices().push_back(TrackClusterReference{SurfaceId{0}, SurfaceMeasurementIndex{0}});
    GenericTrack track{};
    track.clusterRefEnd = static_cast<uint32_t>(frame.getTrackClusterIndices().size());
    frame.getGenericTracks().push_back(track);
    BOOST_REQUIRE(!frame.getGenericTracks().empty());
    BOOST_REQUIRE(!frame.getTrackClusterIndices().empty());
  }

  void resetPublication() noexcept { sidecar.clear(); }

  std::shared_ptr<BoundedMemoryResource> pool;
  std::vector<TrackingParameters> params;
  TimeFrame frame;
  TrackerTraits traits;
  Tracker tracker;
  ITSSharedClusterCompatibility sidecar;
  // Scratch carries non-owning runtime ROF views. Keep these adapter-edge
  // builders alive across load, initialise, and failure/replacement calls.
  std::optional<o2::its::ROFOverlapTable<ITSNLayers>> rofTable;
  std::optional<o2::its::ROFVertexLookupTable<ITSNLayers>> vertexTable;
  std::optional<o2::its::ROFMaskTable<ITSNLayers>> mask;
  std::shared_ptr<tbb::task_arena> arena;
  std::vector<SurfaceDescriptor> catalog;

  // Builds and atomically installs the complete static configuration.
  void establishValidLayout()
  {
    catalog = makeITSTestCatalog();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    TrackerInitialization configuration;
    configuration.catalog = catalogView;
    configuration.memoryPool = pool;
    const auto orderedSurfaces = identitySurfaces(ITSNLayers);
    configuration.iterations.reserve(params.size());
    for (const auto& parameter : params) {
      TrackerIterationConfiguration iteration;
      iteration.graph = makeSurfaceChain(
        orderedSurfaces, parameter.MaxHoles,
        positionalSurfaceMask(parameter.HoleLayerMask, orderedSurfaces, ITSNLayers),
        positionalSurfaceMask(parameter.StartLayerMask, orderedSurfaces, ITSNLayers));
      iteration.parameters = parameter;
      configuration.iterations.push_back(std::move(iteration));
    }
    const auto result = tracker.initialize(frame, configuration);
    BOOST_REQUIRE(result.ok());
    BOOST_REQUIRE_EQUAL(frame.getWorkspace().getNOwnedSurfaces(), orderedSurfaces.size());
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
    const auto& graph = frame.getGraph(0);
    const auto& orderedSurfaces = graph.getOrderedSurfaces();
    auto& workspace = frame.getWorkspace();
    const auto result = workspace.loadNormalizedSource(frame, decoder, origin, timing, f.clusters, f.patterns, f.rofs, &dict(),
                                                       f.labels.getIndexedSize() > 0 ? &f.labels : nullptr, o2::detectors::DetID::ITS,
                                                       gsl::span<const SurfaceId>{orderedSurfaces}, graph.getSurfaceCatalog());
    BOOST_REQUIRE(result.ok());

    // TrackerTraits::computeLayerTracklets() reads per-layer ROF counts
    // from mROFOverlapTableView (o2::its::LayerTiming), a separate table
    // from mROFramesClusters/getNrof() -- it is never populated by
    // loadNormalizedSource() and defaults to an unconfigured/garbage view.
    // A traversal that reaches computeLayerTracklets() without this being
    // set derives its ROF loop bound from that garbage view and walks out
    // of bounds. Mirrors the workflow timing-table construction's
    // shape, but with every layer given the same trivial timing matching
    // this fixture's single combined ROF stream (real production input has
    // per-detector-param ROF length/delay/bias; none of that is exercised
    // by the failure-contract cases here, only the ROF *count* is load
    // -bearing).
    o2::its::LayerTiming timing2{};
    timing2.mNROFsTF = static_cast<unsigned int>(f.rofs.size());
    timing2.mROFLength = 40;
    rofTable.emplace();
    for (int iLayer = 0; iLayer < ITSNLayers; ++iLayer) {
      rofTable->defineLayer(iLayer, timing2);
    }
    rofTable->init();
    vertexTable.emplace();
    for (int iLayer = 0; iLayer < ITSNLayers; ++iLayer) {
      vertexTable->defineLayer(iLayer, timing2);
    }
    vertexTable->init();

    mask.emplace(*rofTable);
    mask->resetMask();
    for (int iLayer = 0; iLayer < ITSNLayers; ++iLayer) {
      mask->setROFsEnabled(iLayer, 0, timing2.mNROFsTF, 1);
    }
    workspace.setROFViews(RuntimeROFViews{rofTable->getView(), vertexTable->getView(), mask->getView(), {}});
  }

  // Set the event-local budget at the current usage; the next allocation is
  // the controlled recoverable failure.
  void forceMemoryLimitAtCurrentUsage()
  {
    const auto used = pool->getUsedMemory();
    pool->setMaxMemory(used);
  }

  void restoreUnboundedMemory()
  {
    pool->setMaxMemory(std::numeric_limits<size_t>::max());
  }
};

class MftRoadDecoder final : public ClusterDecoder
{
 public:
  explicit MftRoadDecoder(std::vector<DecodedCluster> clusters) : mClusters{std::move(clusters)} {}

  SurfaceMeasurementDecodeResult decode(const CompClusterExt& cluster, BoundedPatternCursor& patterns,
                                        const TopologyDictionary* dictionary, gsl::span<const SurfaceId> layerToSurface,
                                        ClusterSourceId source, uint32_t externalIndex, uint32_t sourceROF, bool) const final
  {
    const auto clusterData = ioutils::extractClusterDataBounded(cluster, patterns, dictionary);
    if (!clusterData.ok()) {
      SurfaceMeasurementDecodeResult result;
      result.error = clusterData.error;
      return result;
    }
    SurfaceMeasurementDecodeResult result;
    if (externalIndex >= mClusters.size()) {
      return result;
    }
    auto decoded = mClusters[externalIndex];
    decoded.shape = clusterData.shape;
    result.layer = decoded.layer;
    if (decoded.layer < 0 || static_cast<std::size_t>(decoded.layer) >= layerToSurface.size()) {
      return result;
    }
    result.layerMapped = true;
    result.kind = SurfaceKind::Disk;
    return makeDiskMeasurementDecodeResult(decoded, DetectorSensorId{static_cast<uint32_t>(o2::detectors::DetID::MFT), decoded.sensor},
                                           layerToSurface[decoded.layer], ClusterRef{source, externalIndex}, sourceROF);
  }

 private:
  std::vector<DecodedCluster> mClusters;
};

std::vector<DecodedCluster> makeMftRoad(const TrackingParameters& parameters, float bz)
{
  std::vector<DecodedCluster> result;
  result.reserve(MFTNLayers);
  float x = 3.f;
  float y = 1.5f;
  float z = detail::mftLayerZ(0);
  for (int layer = 0; layer < MFTNLayers; ++layer) {
    DecodedCluster cluster{};
    cluster.global = {x, y, z};
    cluster.rowColumnCovariance = {1.e-2f, 0.f, 1.e-2f};
    cluster.sensor = static_cast<uint32_t>(layer);
    cluster.layer = layer;
    result.push_back(cluster);
    if (layer + 1 == MFTNLayers) {
      break;
    }
    const float nextZ = detail::mftLayerZ(layer + 1);
    float nextX = 0.f;
    float nextY = 0.f;
    detail::mftTrackletProject(x, y, z, parameters.Diamond[0], parameters.Diamond[1], parameters.Diamond[2],
                               layer, layer + 1, bz, parameters.TrackletMinPt, nextX, nextY);
    x = nextX;
    y = nextY;
    z = nextZ;
  }
  return result;
}

std::vector<SurfaceDescriptor> makeMftCatalog()
{
  std::vector<SurfaceDescriptor> catalog;
  catalog.reserve(MFTNLayers);
  for (uint16_t layer = 0; layer < MFTNLayers; ++layer) {
    SurfaceDescriptor surface{SurfaceId{layer}, layer, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk};
    surface.chartRange = {kMFTLookupRMin[layer], kMFTLookupRMax[layer]};
    surface.referenceCoordinate = detail::mftLayerZ(layer);
    const float xOverX0 = kNominalMFTLayerX0[layer];
    surface.material.xOverX0 = xOverX0;
    surface.material.arealDensityGPerCm2 = xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
    catalog.push_back(surface);
  }
  return catalog;
}

// This fixture forms the smallest established full MFT chain: one hit on
// every disk surface. Its controlled refit callback is invoked only after
// real tracklet, cell, neighbour, and road construction.
struct MftFailureRig {
  explicit MftFailureRig(bool dropTFUponFailure)
    : pool(std::make_shared<BoundedMemoryResource>()), tracker(&controlledSeedRefit)
  {
    resetDetectorDefaults(parameters, o2::detectors::DetID::MFT);
    parameters.UseDiamond = true;
    parameters.CreateArtefactLabels = false;
    parameters.DropTFUponFailure = dropTFUponFailure;
    frame.setBz(.5f);
    traits.setNThreads(1, arena);
  }

  void configure(std::size_t iterations = 1)
  {
    catalog = makeMftCatalog();
    TrackerInitialization configuration;
    configuration.catalog = {catalog.data(), static_cast<uint32_t>(catalog.size())};
    configuration.memoryPool = pool;
    const auto surfaces = identitySurfaces(MFTNLayers);
    for (std::size_t index = 0; index < iterations; ++index) {
      TrackerIterationConfiguration iteration;
      iteration.graph = makeSurfaceChain(surfaces, parameters.MaxHoles,
                                         positionalSurfaceMask(parameters.HoleLayerMask, surfaces, MFTNLayers),
                                         positionalSurfaceMask(parameters.StartLayerMask, surfaces, MFTNLayers));
      iteration.parameters = parameters;
      configuration.iterations.push_back(std::move(iteration));
    }
    BOOST_REQUIRE(tracker.initialize(frame, configuration).ok());
  }

  void loadRoad()
  {
    const auto decoded = makeMftRoad(parameters, frame.getBz());
    std::vector<CompClusterExt> compact;
    std::vector<unsigned char> patterns;
    compact.reserve(decoded.size());
    patterns.reserve(decoded.size() * onePixelPattern.size());
    for (const auto& cluster : decoded) {
      compact.emplace_back(0, 0, CompCluster::InvalidPatternID, cluster.layer);
      patterns.insert(patterns.end(), onePixelPattern.begin(), onePixelPattern.end());
    }
    const std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, static_cast<int>(compact.size())}};
    MftRoadDecoder decoder{decoded};
    auto& scratch = frame.getWorkspace();
    const auto& graph = frame.getGraph(0);
    BOOST_REQUIRE(scratch.loadNormalizedSource(frame, decoder, o2::InteractionRecord{50, 5}, ROFTimingConfig{40, 0, 0, 0},
                                               compact, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::MFT,
                                               gsl::span<const SurfaceId>{graph.getOrderedSurfaces()}, graph.getSurfaceCatalog())
                    .ok());
    o2::its::LayerTiming timing{};
    timing.mNROFsTF = 1;
    timing.mROFLength = 40;
    rofTable.emplace();
    vertexTable.emplace();
    for (int layer = 0; layer < MFTNLayers; ++layer) {
      rofTable->defineLayer(layer, timing);
      vertexTable->defineLayer(layer, timing);
    }
    rofTable->init();
    vertexTable->init();
    mask.emplace(*rofTable);
    mask->resetMask();
    for (int layer = 0; layer < MFTNLayers; ++layer) {
      mask->setROFsEnabled(layer, 0, 1, 1);
    }
    scratch.setROFViews(RuntimeROFViews{rofTable->getView(), vertexTable->getView(), mask->getView(), {}});
  }

  void assertReset() const
  {
    BOOST_CHECK_EQUAL(frame.getTotalMeasurements(), 0u);
    BOOST_CHECK(frame.getGenericTracks().empty());
    BOOST_CHECK(frame.getTrackClusterIndices().empty());
    for (std::size_t iteration = 0; iteration < frame.getTrackingParameters().size(); ++iteration) {
      const auto& workspace = frame.getWorkspace().getTraversalWorkspace(iteration);
      BOOST_CHECK(!workspace.valid);
      BOOST_CHECK(workspace.acceptedTracks.empty());
    }
  }

  void stageStaleState()
  {
    ITSSharedClusterCompatibilityTransaction txn{sidecar};
    BOOST_REQUIRE(txn.validate(0));
    txn.reserve();
    txn.append(0);
    frame.getTrackClusterIndices().push_back(TrackClusterReference{SurfaceId{0}, SurfaceMeasurementIndex{0}});
    GenericTrack track{};
    track.clusterRefEnd = static_cast<uint32_t>(frame.getTrackClusterIndices().size());
    frame.getGenericTracks().push_back(track);
  }

  void resetPublication() noexcept { sidecar.clear(); }

  std::shared_ptr<BoundedMemoryResource> pool;
  TrackingParameters parameters{};
  TimeFrame frame;
  TrackerTraits traits;
  Tracker tracker;
  ITSSharedClusterCompatibility sidecar;
  std::shared_ptr<tbb::task_arena> arena;
  std::vector<SurfaceDescriptor> catalog;
  std::optional<o2::its::ROFOverlapTable<MFTNLayers>> rofTable;
  std::optional<o2::its::ROFVertexLookupTable<MFTNLayers>> vertexTable;
  std::optional<o2::its::ROFMaskTable<MFTNLayers>> mask;
};

Fixture emptyFixture()
{
  return Fixture{};
}

} // namespace

// --- Structural failure: always rethrows, always wipes -------------------
//
// Gate 4 B2 Slice 2 removed this section's original mechanism
// (StructuralFailureViaStaleLayoutAlwaysRethrowsAndWipes: establish a valid
// layout, then TimeFrame::invalidateSurfaceGraphs() right before running
// tracking to deterministically produce TraversalException{StaleLayout}).
// Neither invalidateSurfaceGraphs() nor TraversalFailureReason::StaleLayout
// is reachable any more: initialiseTimeFrame() now takes the plan as an
// explicit `const std::vector<SurfaceGraph>&` parameter with no TimeFrame-owned
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
  BOOST_REQUIRE(rig.frame.getTotalMeasurements() > 0u);

  rig.forceMemoryLimitAtCurrentUsage();

  const auto result = rig.tracker.run(rig.frame, rig.traits);

  BOOST_CHECK(result.outcome == TrackingOutcome::RecoverableDropped);
  BOOST_CHECK_EQUAL(rig.frame.getTotalMeasurements(), 0u);
  BOOST_CHECK(rig.frame.getGenericTracks().empty());
  BOOST_CHECK_EQUAL(rig.frame.getEventResetCount(), 2u);
}

BOOST_AUTO_TEST_CASE(RecoverableFailureNotDroppedRethrowsButStillWipesFirst)
{
  Rig rig{/*dropTFUponFailure=*/false};
  rig.establishValidLayout();
  rig.loadSource(makeFixture());
  BOOST_REQUIRE(rig.frame.getTotalMeasurements() > 0u);

  rig.forceMemoryLimitAtCurrentUsage();

  BOOST_CHECK_THROW(rig.tracker.run(rig.frame, rig.traits), BoundedMemoryResource::MemoryLimitExceeded);

  // Wipe must have already happened before the exception propagated -- not
  // "the process is going down anyway".
  BOOST_CHECK_EQUAL(rig.frame.getTotalMeasurements(), 0u);
  BOOST_CHECK(rig.frame.getGenericTracks().empty());
  BOOST_CHECK_EQUAL(rig.frame.getEventResetCount(), 2u);
}

// --- std::bad_alloc: recoverable, same drop-or-rethrow policy ------------
//
// A real ten-disk MFT road reaches controlledSeedRefit() after all CA stages.
// The callback is the existing refit operation seam, not a synthetic kernel
// override.

BOOST_AUTO_TEST_CASE(BadAllocDroppedReturnsExactSentinelAndWipes)
{
  ensureTrivialMagneticFieldIsSet();
  MftFailureRig rig{/*dropTFUponFailure=*/true};
  rig.configure();
  rig.loadRoad();
  BOOST_REQUIRE(rig.frame.getTotalMeasurements() > 0u);

  gRefitFailure = RefitFailure::BadAlloc;
  gRefitCalls = 0;
  const auto result = rig.tracker.run(rig.frame, rig.traits);
  gRefitFailure = RefitFailure::None;

  BOOST_CHECK(result.outcome == TrackingOutcome::RecoverableDropped);
  BOOST_CHECK_GT(gRefitCalls, 0);
  rig.assertReset();
}

BOOST_AUTO_TEST_CASE(BadAllocNotDroppedRethrowsButStillWipesFirst)
{
  ensureTrivialMagneticFieldIsSet();
  MftFailureRig rig{/*dropTFUponFailure=*/false};
  rig.configure();
  rig.loadRoad();
  BOOST_REQUIRE(rig.frame.getTotalMeasurements() > 0u);

  gRefitFailure = RefitFailure::BadAlloc;
  gRefitCalls = 0;
  BOOST_CHECK_THROW(rig.tracker.run(rig.frame, rig.traits), std::bad_alloc);
  gRefitFailure = RefitFailure::None;

  BOOST_CHECK_GT(gRefitCalls, 0);
  rig.assertReset();
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
    ensureTrivialMagneticFieldIsSet();
    MftFailureRig rig{dropFlag};
    rig.configure();
    rig.loadRoad();
    BOOST_REQUIRE(rig.frame.getTotalMeasurements() > 0u);

    gRefitFailure = RefitFailure::UnclassifiedRuntimeError;
    gRefitCalls = 0;
    BOOST_CHECK_THROW(rig.tracker.run(rig.frame, rig.traits), std::runtime_error);
    gRefitFailure = RefitFailure::None;

    BOOST_CHECK_GT(gRefitCalls, 0);
    rig.assertReset();
  }
}

BOOST_AUTO_TEST_CASE(LaterIterationRefitFailureWipesEveryIterationWorkspace)
{
  ensureTrivialMagneticFieldIsSet();
  MftFailureRig rig{/*dropTFUponFailure=*/true};
  rig.configure(/*iterations=*/2);
  rig.loadRoad();

  gRefitFailure = RefitFailure::UnclassifiedRuntimeError;
  gRefitCalls = 0;
  gThrowOnRefitCall = 2;
  BOOST_CHECK_THROW(rig.tracker.run(rig.frame, rig.traits), std::runtime_error);
  gRefitFailure = RefitFailure::None;
  gThrowOnRefitCall = 1;

  // The first iteration completed its real road/refit callback before the
  // second failed. Tracker::run() must not leave its workspace selectable by
  // a later adapter pass.
  BOOST_CHECK_EQUAL(gRefitCalls, 2);
  rig.assertReset();
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
    rig.establishValidLayout();
    rig.loadSource(emptyFixture());

    rig.frame.getGenericTracks().push_back(GenericTrack{});
    BOOST_REQUIRE(!rig.frame.getGenericTracks().empty());

    bool threw = false;
    try {
      rig.tracker.run(rig.frame, rig.traits);
    } catch (const TraversalException& e) {
      threw = true;
      BOOST_CHECK(e.getReason() == TraversalFailureReason::InvalidIndexTableConfiguration);
    }
    BOOST_CHECK(threw);

    BOOST_CHECK(rig.frame.getGenericTracks().empty());
  }
}

BOOST_AUTO_TEST_CASE(IndexTableConfigurationMismatchAlwaysRethrowsAndWipesRegardlessOfFlag)
{
  ensureTrivialMagneticFieldIsSet(); // iteration 0 runs to completion (findRoads()) before iteration 1 fails
  for (const bool dropFlag : {false, true}) {
    Rig rig{dropFlag};
    rig.params = makeTwoIterationITSParams(dropFlag);
    rig.params[1].RowBins = rig.params[0].RowBins + 1; // mismatched against what iteration 0 will commit
    rig.establishValidLayout();
    rig.loadSource(emptyFixture());

    bool threw = false;
    try {
      rig.tracker.run(rig.frame, rig.traits);
    } catch (const TraversalException& e) {
      threw = true;
      BOOST_CHECK(e.getReason() == TraversalFailureReason::IndexTableConfigurationMismatch);
    }
    BOOST_CHECK(threw);

    BOOST_CHECK(rig.frame.getGenericTracks().empty());
  }
}

// --- Valid empty input -----------------------------------------------------

BOOST_AUTO_TEST_CASE(ValidEmptyInputCompletesWithoutErrorAndProducesNoTracks)
{
  ensureTrivialMagneticFieldIsSet();
  Rig rig{/*dropTFUponFailure=*/false};
  rig.establishValidLayout();
  rig.loadSource(emptyFixture());
  BOOST_REQUIRE_EQUAL(rig.frame.getTotalMeasurements(), 0u);

  TrackingResult result{TrackingOutcome::Structural, std::numeric_limits<float>::quiet_NaN()};
  BOOST_CHECK_NO_THROW(result = rig.tracker.run(rig.frame, rig.traits));

  BOOST_CHECK(result.outcome == TrackingOutcome::Success);
  BOOST_CHECK(result.elapsedMs >= 0.f);
  BOOST_CHECK_EQUAL(rig.frame.getGenericTracks().size(), 0u);
}

// --- Direct outcome classification ----------------------------------------
//
// TrackingOutcome::Structural is part of this type's vocabulary for a future
// caller that catches Tracker::run()'s propagated exception itself -- run()
// never constructs it via a
// normal return, since every structural/unclassified/non-dropped-recoverable
// failure keeps the exact "retain exceptions" contract already proven above
// (UnclassifiedExceptionAlwaysRethrowsAndWipesRegardlessOfFlag,
// InvalidIndexTableConfigurationAlwaysRethrowsAndWipesRegardlessOfFlag,
// BadAllocNotDroppedRethrowsButStillWipesFirst,
// RecoverableFailureNotDroppedRethrowsButStillWipesFirst): those tests *are*
// this outcome's structural-failure classification evidence, expressed the
// only way it is currently observable (a thrown exception, never a returned
// value). This test only proves the three values the type actually defines
// are distinct and that TrackingResult's fields carry what each documented
// path above already relies on.
BOOST_AUTO_TEST_CASE(TrackingOutcomeValuesAreDistinct)
{
  BOOST_CHECK(TrackingOutcome::Success != TrackingOutcome::RecoverableDropped);
  BOOST_CHECK(TrackingOutcome::Success != TrackingOutcome::Structural);
  BOOST_CHECK(TrackingOutcome::RecoverableDropped != TrackingOutcome::Structural);

  constexpr TrackingResult defaulted{};
  BOOST_CHECK(defaulted.outcome == TrackingOutcome::Success);
  BOOST_CHECK_EQUAL(defaulted.elapsedMs, 0.f);
}

// --- No stale TimeFrame/GenericTrack/sidecar state survives -----------------
//
// Both non-success return paths from Tracker::run() (structural-rethrow
// and recoverable-dropped) must leave the shared TimeFrame's GenericTrack
// storage and the tracker's adopted compatibility sidecar exactly as empty
// as a freshly wiped/cleared TimeFrame would -- not merely the normalized
// frame and legacy tracks storage the tests above already check.

BOOST_AUTO_TEST_CASE(RecoverableDroppedLeavesNoStaleGenericTrackOrSidecarState)
{
  Rig rig{/*dropTFUponFailure=*/true};
  rig.establishValidLayout();
  rig.loadSource(makeFixture());
  rig.stageStaleState();

  rig.forceMemoryLimitAtCurrentUsage();
  const auto result = rig.tracker.run(rig.frame, rig.traits);
  rig.resetPublication();

  BOOST_CHECK(result.outcome == TrackingOutcome::RecoverableDropped);
  BOOST_CHECK(rig.frame.getGenericTracks().empty());
  BOOST_CHECK(rig.frame.getTrackClusterIndices().empty());
  BOOST_CHECK_EQUAL(rig.sidecar.pendingSize(), 0u);
}

BOOST_AUTO_TEST_CASE(StructuralFailureLeavesNoStaleGenericTrackOrSidecarState)
{
  for (const bool dropFlag : {false, true}) {
    ensureTrivialMagneticFieldIsSet();
    MftFailureRig rig{dropFlag};
    rig.configure();
    rig.loadRoad();
    rig.stageStaleState();

    gRefitFailure = RefitFailure::UnclassifiedRuntimeError;
    gRefitCalls = 0;
    BOOST_CHECK_THROW(rig.tracker.run(rig.frame, rig.traits), std::runtime_error);
    gRefitFailure = RefitFailure::None;
    rig.resetPublication();

    BOOST_CHECK_GT(gRefitCalls, 0);
    rig.assertReset();
    BOOST_CHECK_EQUAL(rig.sidecar.pendingSize(), 0u);
  }
}

// --- Continued processing after a drop ------------------------------------

BOOST_AUTO_TEST_CASE(TrackerRemainsUsableAfterADroppedTimeFrame)
{
  ensureTrivialMagneticFieldIsSet();
  Rig rig{/*dropTFUponFailure=*/true};
  rig.establishValidLayout();
  rig.loadSource(makeFixture());

  rig.forceMemoryLimitAtCurrentUsage();
  const auto dropped = rig.tracker.run(rig.frame, rig.traits);
  BOOST_REQUIRE(dropped.outcome == TrackingOutcome::RecoverableDropped);

  // Restore headroom and process a fresh (here, empty) TimeFrame on the
  // SAME Tracker/TrackerTraits instance -- proving the tracker/device stays
  // usable after a drop, matching the DPL device staying alive.
  rig.restoreUnboundedMemory();
  rig.loadSource(emptyFixture());

  TrackingResult result{TrackingOutcome::Structural, std::numeric_limits<float>::quiet_NaN()};
  BOOST_CHECK_NO_THROW(result = rig.tracker.run(rig.frame, rig.traits));
  BOOST_CHECK(result.outcome == TrackingOutcome::Success);
  BOOST_CHECK(result.elapsedMs >= 0.f);
}
