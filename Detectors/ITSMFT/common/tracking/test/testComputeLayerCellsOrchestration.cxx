// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Orchestration coverage for TrackerTraits<NLayers>::computeLayerCells after
// its detector-family branch was replaced by a one-shot outer dispatch to
// cell-seed leaves (Architecture.md Sec 10/10.1). cell-seed leaves's
// numerical parity with the legacy inline formulas is already proven by
// testTrackletFinding.cxx; this file does not re-derive or
// duplicate that formula. It proves instead that the real public
// computeLayerCells() entry point:
//  - resolves the three clusters/hits/material values for a candidate in
//    strict {inner, middle, outer} order and hands them to cell-seed leaves
//    unchanged (checked against an oracle call to cell-seed leaves itself,
//    built from the same values refetched through the TimeFrame, not
//    re-typed literals);
//  - exercises cylinder and disk cells through the same public orchestration
//    entry point, with coordinate differences confined to cell-seed leaves;
//  - leaves cellTopologyId indexing, the LUT, MC-label construction, and
//    one-pass/two-pass ordering untouched;
//  - fails closed (TraversalException::InvalidTraversalSchedule) through the
//    existing public API alone, with no test-only seam into private
//    traversal-cache state.

#define BOOST_TEST_MODULE ITSMFT ComputeLayerCells orchestration
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <oneapi/tbb/task_arena.h>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITSMFTTracking/detail/CandidateFinding.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Constants.h"
#include "ITStracking/Tracklet.h"
#include "MFTTracking/Constants.h"

#include "TraversalTestSupport.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{
constexpr uint8_t kCompatibilityAbsCharge = 1;
const o2::track::PID kCompatibilityPID = o2::track::PID::Pion;
} // namespace

namespace
{

constexpr float Bz = 0.5f;

// Preflight-only fixtures (Rig::establishLayout()) load zero real clusters --
// this decoder's decode() is never actually invoked there. It exists only to
// satisfy loadNormalizedSource()'s interface, mirroring
// testTrackerFailureContract.cxx's LegacyLikeDecoder.
class NeverDecodedDecoder final : public ClusterDecoder
{
 public:
  explicit NeverDecodedDecoder(o2::detectors::DetID::ID detector) : mDetector(detector) {}

  o2::itsmft::tracking::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt&, BoundedPatternCursor&, const TopologyDictionary*,
    gsl::span<const SurfaceId>, ClusterSourceId, uint32_t, uint32_t, bool) const override
  {
    return {};
  }

 private:
  o2::detectors::DetID::ID mDetector;
};

// Stage-B normalized-CA-measurements slice: computeLayerCells() now reads
// TrackerTraits::mLayerMeasurements, resolved once per initialiseTimeFrame()
// from the TimeFrame's already-loaded normalized frame. Candidate fixtures
// therefore load their three clusters through the real loadNormalizedSource()
// path -- backfilling both the normalized frame and every legacy
// compatibility structure (unsorted clusters, TrackingFrameInfo, external
// indices, ROF boundaries) together, in lockstep -- rather than poking legacy
// structures directly. This decoder returns exactly the caller-supplied
// SurfaceMeasurement for a given detector-local layer (encoded as the
// synthetic CompClusterExt's chipID/sensorID), verbatim except for the
// surface/sensor/cluster/sourceROF identity fields loadSources() itself
// validates against its own call parameters.
class FixedMeasurementDecoder final : public ClusterDecoder
{
 public:
  struct MeasurementPair {
    GlobalMeasurement global{};
    SurfaceMeasurement local{};
  };

  FixedMeasurementDecoder(o2::detectors::DetID::ID detector, SurfaceKind kind) : mDetector(detector), mKind(kind) {}

  void setMeasurement(int layer, const MeasurementPair& measurement) { mByLayer[layer] = measurement; }

  o2::itsmft::tracking::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor&,
    const TopologyDictionary*,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool) const override
  {
    o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
    const int layer = cluster.getSensorID();
    result.layer = layer;
    if (layer < 0 || static_cast<size_t>(layer) >= layerToSurface.size()) {
      return result;
    }
    result.layerMapped = true;
    result.kind = mKind;
    const auto it = mByLayer.find(layer);
    BOOST_REQUIRE(it != mByLayer.end());
    auto measurement = it->second;
    measurement.global.surface = layerToSurface[layer];
    measurement.global.sensor = DetectorSensorId{static_cast<uint32_t>(mDetector), static_cast<uint32_t>(layer)};
    measurement.global.cluster = ClusterRef{source, externalIndex};
    measurement.global.sourceROF = sourceROF;
    result.global = measurement.global;
    result.measurement = measurement.local;
    return result;
  }

 private:
  o2::detectors::DetID::ID mDetector;
  SurfaceKind mKind;
  std::map<int, MeasurementPair> mByLayer;
};

const TopologyDictionary& dict()
{
  static const TopologyDictionary d;
  return d;
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

// `layerxX0` supplies each surface's nominal material to embed in the
// catalog, in legacy layer order -- callers pass the exact TrackingParameters::
// LayerxX0 values they intend to run with, so TrackerTraits::
// initialiseTimeFrame()'s LegacyMaterialMismatch compatibility check (which
// compares the two) always passes here; this file exercises computeLayerCells()
// orchestration, not that compatibility check.
std::vector<SurfaceDescriptor> makeCatalog(uint16_t nLayers, o2::detectors::DetID::ID det, SurfaceKind kind, gsl::span<const float> layerxX0)
{
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.reserve(nLayers);
  for (uint16_t i = 0; i < nLayers; ++i) {
    surfaces.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(det), kind});
    surfaces.back().chartRange = kind == SurfaceKind::Disk ? SurfaceChartRange{0.1f, 20.f} : SurfaceChartRange{-20.f, 20.f};
    surfaces.back().referenceCoordinate = kind == SurfaceKind::Cylinder
                                            ? 3.f + static_cast<float>(i)
                                            : -0.4f - 0.2f * static_cast<float>(i);
    const float xOverX0 = layerxX0[i];
    surfaces.back().material.xOverX0 = xOverX0;
    surfaces.back().material.arealDensityGPerCm2 = xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
  }
  return surfaces;
}

NominalSurfaceMaterial toMaterial(float xOverX0)
{
  return NominalSurfaceMaterial{xOverX0, xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho};
}

std::array<NominalSurfaceMaterial, 3> toMaterial(const std::array<float, 3>& xOverX0)
{
  return {toMaterial(xOverX0[0]), toMaterial(xOverX0[1]), toMaterial(xOverX0[2])};
}

// Same construction as testTrackletFinding.cxx's helpers -- plain
// input-struct builders, not a reimplementation of any fit formula.
o2::its::Cluster makeGlobalCluster(float x, float y, float z, int id = 0)
{
  return o2::its::Cluster{x, y, z, id};
}

o2::its::TrackingFrameInfo makeBarrelHit(float xTF, float alpha, float y, float z, float sigma2Y = 1.e-4f, float sigma2Z = 1.e-4f)
{
  return o2::its::TrackingFrameInfo{xTF, y, z, xTF, alpha, {y, z}, {sigma2Y, 0.f, sigma2Z}};
}

o2::its::TrackingFrameInfo makeDiskHit(float z, float x, float y, float sigma2X = 1.e-2f, float sigma2Y = 1.e-2f)
{
  return o2::its::TrackingFrameInfo{x, y, z, 0.f, 0.f, {x, y}, {sigma2X, 0.f, sigma2Y}};
}

// Test-local field-mapping helpers (not a production API), matching the same
// Cylinder/Disk field mapping used by the production migration and by
// testCellFinding.cxx: the single SurfaceMeasurement now
// standing in for the retired {Cluster, TrackingFrameInfo} pair at each
// candidate position.
FixedMeasurementDecoder::MeasurementPair barrelMeasurementFor(const o2::its::Cluster& cluster, const o2::its::TrackingFrameInfo& hit)
{
  FixedMeasurementDecoder::MeasurementPair measurement{};
  measurement.global.position = {cluster.xCoordinate, cluster.yCoordinate, cluster.zCoordinate};
  measurement.global.radius = std::hypot(cluster.xCoordinate, cluster.yCoordinate);
  const float sine = std::sin(hit.alphaTrackingFrame);
  const float cosine = std::cos(hit.alphaTrackingFrame);
  const float varianceU = hit.covarianceTrackingFrame[0];
  const float covarianceUV = hit.covarianceTrackingFrame[1];
  measurement.global.covariance = {sine * sine * varianceU,
                                   -sine * cosine * varianceU,
                                   -sine * covarianceUV,
                                   cosine * cosine * varianceU,
                                   cosine * covarianceUV,
                                   hit.covarianceTrackingFrame[2]};
  measurement.local.frame.q = hit.xTrackingFrame;
  measurement.local.frame.frameAngle = hit.alphaTrackingFrame;
  measurement.local.frame.u = hit.positionTrackingFrame[0];
  measurement.local.frame.v = hit.positionTrackingFrame[1];
  measurement.local.covariance.uu = hit.covarianceTrackingFrame[0];
  measurement.local.covariance.uv = hit.covarianceTrackingFrame[1];
  measurement.local.covariance.vv = hit.covarianceTrackingFrame[2];
  return measurement;
}

FixedMeasurementDecoder::MeasurementPair diskMeasurementFor(const o2::its::Cluster& cluster, const o2::its::TrackingFrameInfo& hit)
{
  FixedMeasurementDecoder::MeasurementPair measurement{};
  measurement.global.position = {cluster.xCoordinate, cluster.yCoordinate, cluster.zCoordinate};
  measurement.global.radius = std::hypot(cluster.xCoordinate, cluster.yCoordinate);
  measurement.global.covariance = {hit.covarianceTrackingFrame[0], 0.f, 0.f,
                                   hit.covarianceTrackingFrame[2], 0.f, 0.f};
  measurement.local.frame = {cluster.zCoordinate, cluster.xCoordinate, cluster.yCoordinate, 0.f};
  measurement.local.covariance.uu = hit.covarianceTrackingFrame[0];
  measurement.local.covariance.uv = 0.f;
  measurement.local.covariance.vv = hit.covarianceTrackingFrame[2];
  return measurement;
}

// Stage-B activation: the produced Cell no longer inherits a track
// parametrization, so the oracle comparison is done directly on
// SurfaceKinematicState (both the produced cell's own .state() and the
// oracle's native cell-seed leaves output are this type).
void checkSurfaceKinematicStateEqual(const SurfaceKinematicState& lhs, const SurfaceKinematicState& rhs)
{
  for (int i = 0; i < 5; ++i) {
    BOOST_CHECK_EQUAL(lhs.parameters[i], rhs.parameters[i]);
  }
  for (int i = 0; i < 15; ++i) {
    BOOST_CHECK_EQUAL(lhs.covariance[i], rhs.covariance[i]);
  }
  BOOST_CHECK_EQUAL(lhs.referenceCoordinate, rhs.referenceCoordinate);
  BOOST_CHECK_EQUAL(lhs.alpha, rhs.alpha);
}

// Minimal wiring TrackerTraits<NLayers>::computeLayerCells() needs: a real
// layout/topology (so initialiseTimeFrame() genuinely binds
// the link/cell schedule and tracking parameters
// -- computeLayerCells()'s own private caches, never poked directly), and a
// validly-sized-but-empty normalized load (proven pattern from
// testTrackerFailureContract.cxx: TimeFrame::initialise() unconditionally
// reads mROFramesClusters sizes, which only loadNormalizedSource() sets up
// safely, even for zero clusters).
struct RigFrameStorage {
  RigFrameStorage() : pool(std::make_shared<BoundedMemoryResource>()) { frame.setMemoryPool(pool); }

  std::shared_ptr<BoundedMemoryResource> pool;
  TimeFrame frame;
};

template <int NLayers>
struct Rig : RigFrameStorage {

  Rig(o2::detectors::DetID::ID det, SurfaceKind kind, int nThreads = 1)
    : params(1),
      mDet(det),
      mKind(kind)
  {
    resetDetectorDefaults(params[0], det);
    // This file bypasses computeLayerTracklets()'s phi/z/index-table cuts
    // entirely (candidates are injected directly, see
    // injectCandidateTracklets() below): clearing RebuildClusterLUT keeps
    // TimeFrame::initialise() from also exercising prepareClusters()'s
    // index-table row/col binning and ROF-mask lookup on the synthetic
    // candidate positions/ROF this file uses -- that out-of-scope subsystem
    // is not configured for this file's candidates (in particular, no
    // multiplicity/UPC ROF mask is ever loaded, so its default view is
    // never a valid one to index once real clusters are present, unlike
    // when this file loaded zero real clusters).
    params[0].PassFlags.reset(IterationStep::RebuildClusterLUT);
    traits.setNThreads(nThreads, arena);
    frame.setBz(Bz);
  }

  // Establishes the catalog/layout and loads a (zero-cluster) normalized
  // source. Deliberately not run by the constructor: it builds the catalog's
  // nominal material from the *current* params[0].LayerxX0, so callers must
  // finish any TrackingParameters::LayerxX0 override first -- matching
  // TrackerTraits::initialiseTimeFrame()'s LegacyMaterialMismatch
  // compatibility check, which this file does not itself exercise.
  void establishLayout()
  {
    catalog = makeCatalog(static_cast<uint16_t>(NLayers), mDet, mKind, gsl::span<const float>(params[0].LayerxX0));
    const auto orderedSurfaces = identitySurfaces(static_cast<uint16_t>(NLayers));
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    TrackerInitialization configuration;
    configuration.catalog = catalogView;
    configuration.memoryPool = pool;
    TrackerIterationConfiguration iteration;
    iteration.graph = makeSurfaceChain(
      orderedSurfaces, params[0].MaxHoles,
      positionalSurfaceMask(params[0].HoleLayerMask, orderedSurfaces, NLayers),
      positionalSurfaceMask(params[0].StartLayerMask, orderedSurfaces, NLayers));
    iteration.parameters = params[0];
    configuration.iterations.push_back(std::move(iteration));
    BOOST_REQUIRE(tracker.initialize(frame, configuration).ok());
    tf = &frame.getWorkspace();
    const auto& graph = frame.getGraph(0);

    NeverDecodedDecoder decoder{mDet};
    const o2::InteractionRecord origin{50, 5};
    const ROFTimingConfig timing{40, 0, 0, 0};
    const std::vector<CompClusterExt> noClusters;
    const std::vector<unsigned char> noPatterns;
    const std::vector<ROFRecord> noRofs;
    const auto& loadOrderedSurfaces = graph.getOrderedSurfaces();
    const auto loadResult = tf->loadNormalizedSource(frame, decoder, origin, timing, noClusters, noPatterns, noRofs, &dict(), nullptr, mDet,
                                                     gsl::span<const SurfaceId>{loadOrderedSurfaces}, graph.getSurfaceCatalog());
    BOOST_REQUIRE(loadResult.ok());
  }

  o2::detectors::DetID::ID detector() const noexcept { return mDet; }
  SurfaceKind kind() const noexcept { return mKind; }

  std::vector<TrackingParameters> params;
  // Gate 4 B3.1: `frame` declared before `tf` so it is constructed first and
  // destroyed last (see SurfaceTrackingScratch's own lifetime-contract doc).
  SurfaceTrackingScratch* tf{nullptr};
  Tracker tracker;
  TrackerTraits traits;
  std::shared_ptr<tbb::task_arena> arena;
  // Must outlive `plan` (std::vector<SurfaceGraph> borrows a SurfaceCatalogView into
  // it, Gate 4 B2 Slice 2) -- declared before `plan` so it is constructed
  // first and destroyed last.
  std::vector<SurfaceDescriptor> catalog;

 private:
  o2::detectors::DetID::ID mDet;
  SurfaceKind mKind;
};

template <int NLayers>
TraversalWorkspaceView prepare(Rig<NLayers>& rig)
{
  return TrackerTestAccess::prepare(rig.tracker, rig.frame, 0);
}

// Loads exactly the three supplied {cluster, hit} candidates at legacy
// layers {0, 1, 2} (every test in this file locates its candidate cell via
// findCellTopologyId(topology, 0, 1, 2), so the layer mapping is always this
// identity triple) through the real loadNormalizedSource() path, via
// FixedMeasurementDecoder -- so the normalized frame and every legacy
// compatibility structure are populated together, in lockstep, exactly as
// TrackerTraits::initialiseTimeFrame()'s one-time normalized-measurement
// binding requires. Must be called after Rig::establishLayout() (which needs
// the catalog/topology first) and before TrackerTraits::initialiseTimeFrame()
// (which validates the normalized frame against the legacy structures this
// call also populates).
template <int NLayers>
void loadCandidateClusters(Rig<NLayers>& rig,
                           const std::array<o2::its::Cluster, 3>& clusters,
                           const std::array<o2::its::TrackingFrameInfo, 3>& hits)
{
  const bool isDisk = rig.kind() == SurfaceKind::Disk;
  FixedMeasurementDecoder decoder{rig.detector(), rig.kind()};
  std::vector<CompClusterExt> compClusters;
  compClusters.reserve(3);
  for (int layer = 0; layer < 3; ++layer) {
    compClusters.emplace_back(0, 0, CompCluster::InvalidPatternID, static_cast<uint16_t>(layer));
    const auto measurement = isDisk ? diskMeasurementFor(clusters[layer], hits[layer]) : barrelMeasurementFor(clusters[layer], hits[layer]);
    decoder.setMeasurement(layer, measurement);
  }
  const std::vector<unsigned char> noPatterns;
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, 3}};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0};
  const auto& orderedSurfaces = rig.frame.getGraph(0).getOrderedSurfaces();
  const auto result = rig.tf->loadNormalizedSource(rig.frame, decoder, origin, timing, compClusters, noPatterns, rofs, &dict(), nullptr, rig.detector(),
                                                   gsl::span<const SurfaceId>{orderedSurfaces}, rig.frame.getGraph(0).getSurfaceCatalog());
  BOOST_REQUIRE(result.ok());
}

// Finds the cellTopologyId whose two links span exactly
// inner->middle->outer, without assuming any particular enumeration order
// out of the sparse topology's builder enumeration.
template <typename TopologyView>
int findCellTopologyId(const TopologyView& topology, int inner, int middle, int outer)
{
  for (int i = 0; i < topology.nCells; ++i) {
    const auto& cell = topology.getCell(CellTopologyId{static_cast<uint16_t>(i)});
    const auto& first = topology.getLink(cell.firstLink);
    const auto& second = topology.getLink(cell.secondLink);
    if (first.from.value() == inner && first.to.value() == middle && second.from.value() == middle && second.to.value() == outer) {
      return i;
    }
  }
  return -1;
}

// Bypasses the real (untouched, out-of-scope-for-this-change)
// computeLayerTracklets() phi/z/index-table cuts entirely: RebuildClusterLUT
// is cleared (Rig's constructor), so TimeFrame::initialise() resizes
// mClusters[layer] to match the single already-loaded unsorted cluster but
// leaves it default-constructed (clusterId == UnusedIndex) rather than
// sorting real content into it. This helper assigns the real sorted-cluster
// entry at index 0 directly (matching loadCandidateClusters()'s own
// unsorted-index-0 cluster on each layer) and injects exactly one synthetic
// tracklet per link of cellTopologyId, wired so
// computeLayerCellsForKind<Tag>'s tracklet-pairing loop finds exactly one
// candidate pair.
template <int NLayers>
void injectCandidateTracklets(Rig<NLayers>& rig, int cellTopologyId, const std::array<o2::its::Cluster, 3>& clusters)
{
  const auto topology = rig.frame.getGraph(0).getView();
  const auto& cell = topology.getCell(CellTopologyId{static_cast<uint16_t>(cellTopologyId)});
  const auto& first = topology.getLink(cell.firstLink);
  const auto& second = topology.getLink(cell.secondLink);
  const int layers[3] = {first.from.value(), first.to.value(), second.to.value()};

  for (int i = 0; i < 3; ++i) {
    BOOST_REQUIRE_EQUAL(rig.tf->getClusters()[layers[i]].size(), 1u);
    rig.tf->getClusters()[layers[i]][0] = clusters[i];
  }

  const o2::its::TimeEstBC ts{static_cast<uint32_t>(0), static_cast<uint16_t>(1)};
  const float firstPhi = std::atan2(clusters[0].yCoordinate - clusters[1].yCoordinate,
                                    clusters[0].xCoordinate - clusters[1].xCoordinate);
  const float secondPhi = std::atan2(clusters[1].yCoordinate - clusters[2].yCoordinate,
                                     clusters[1].xCoordinate - clusters[2].xCoordinate);
  rig.tf->getTracklets()[cell.firstLink.value()].push_back(o2::its::Tracklet{0, 0, 0.f, firstPhi, ts});
  rig.tf->getTracklets()[cell.secondLink.value()].push_back(o2::its::Tracklet{0, 0, 0.f, secondPhi, ts});

  auto& secondLUT = rig.tf->getTrackletsLookupTable()[cell.secondLink.value()];
  secondLUT.resize(2);
  secondLUT[0] = 0;
  secondLUT[1] = 1;
}

// Gate 4 Slice 0b additions below: multi-cell parity coverage for the
// migrated computeLayerCells()/computeLayerCellsForKind(), extending this
// file's existing single-cell (always layers {0,1,2}) machinery to an
// arbitrary ordered set of N>=3 layers so several simultaneously-populated
// cells (sharing links between adjacent triples) can be checked in one
// run.

// Same technique as loadCandidateClusters() (real loadNormalizedSource()
// path via FixedMeasurementDecoder), generalized to N candidate layers
// instead of the fixed {0,1,2} triple.
template <int NLayers, size_t N>
void loadCandidateClustersAtLayers(Rig<NLayers>& rig,
                                   const std::array<int, N>& layers,
                                   const std::array<o2::its::Cluster, N>& clusters,
                                   const std::array<o2::its::TrackingFrameInfo, N>& hits)
{
  const bool isDisk = rig.kind() == SurfaceKind::Disk;
  FixedMeasurementDecoder decoder{rig.detector(), rig.kind()};
  std::vector<CompClusterExt> compClusters;
  compClusters.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    compClusters.emplace_back(0, 0, CompCluster::InvalidPatternID, static_cast<uint16_t>(layers[i]));
    const auto measurement = isDisk ? diskMeasurementFor(clusters[i], hits[i]) : barrelMeasurementFor(clusters[i], hits[i]);
    decoder.setMeasurement(layers[i], measurement);
  }
  const std::vector<unsigned char> noPatterns;
  const std::vector<ROFRecord> rofs{ROFRecord{{0, 0}, 0, 0, static_cast<int>(N)}};
  const o2::InteractionRecord origin{50, 5};
  const ROFTimingConfig timing{40, 0, 0, 0};
  const auto& orderedSurfaces = rig.frame.getGraph(0).getOrderedSurfaces();
  const auto result = rig.tf->loadNormalizedSource(rig.frame, decoder, origin, timing, compClusters, noPatterns, rofs, &dict(), nullptr, rig.detector(),
                                                   gsl::span<const SurfaceId>{orderedSurfaces}, rig.frame.getGraph(0).getSurfaceCatalog());
  BOOST_REQUIRE(result.ok());
}

// Finds the linkId spanning exactly from->to, mirroring
// findCellTopologyId()'s linear-search style over the legacy view.
template <typename TopologyView>
int findLinkId(const TopologyView& topology, int from, int to)
{
  for (int i = 0; i < topology.nLinks; ++i) {
    const auto& t = topology.getLink(LinkId{static_cast<uint16_t>(i)});
    if (t.from.value() == from && t.to.value() == to) {
      return i;
    }
  }
  return -1;
}

// Generalizes injectCandidateTracklets() to an ordered chain of N>=3 layers:
// writes each physical layer's single cluster exactly once, then touches
// each of the N-1 adjacent-pair links exactly once (one synthetic
// tracklet + one LUT {0,1}), regardless of how many downstream cells in the
// chain share that link. Naively calling the single-cell
// injectCandidateTracklets() once per overlapping cell would instead
// double-write any shared link (extra duplicate tracklet, and a LUT
// left however the last call set it) and silently clobber a shared physical
// layer's cluster across calls -- this helper touches every physical layer
// and every link exactly once, by construction.
template <int NLayers, size_t N>
void injectChainCandidateTracklets(Rig<NLayers>& rig, const std::array<int, N>& layers, const std::array<o2::its::Cluster, N>& clusters)
{
  static_assert(N >= 3, "a chain needs at least 3 layers to form one cell");
  const auto topology = rig.frame.getGraph(0).getView();
  for (size_t i = 0; i < N; ++i) {
    BOOST_REQUIRE_EQUAL(rig.tf->getClusters()[layers[i]].size(), 1u);
    rig.tf->getClusters()[layers[i]][0] = clusters[i];
  }

  const o2::its::TimeEstBC ts{static_cast<uint32_t>(0), static_cast<uint16_t>(1)};
  for (size_t i = 0; i + 1 < N; ++i) {
    const int linkId = findLinkId(topology, layers[i], layers[i + 1]);
    BOOST_REQUIRE_GE(linkId, 0);
    rig.tf->getTracklets()[linkId].push_back(o2::its::Tracklet{0, 0, 0.f, 0.f, ts});
    auto& lut = rig.tf->getTrackletsLookupTable()[linkId];
    lut.resize(2);
    lut[0] = 0;
    lut[1] = 1;
  }
}

} // namespace

// --- Barrel: real orchestration matches the cell-seed leaves oracle -----

BOOST_AUTO_TEST_CASE(CylinderComputeLayerCellsMatchesBuildCellSeedOracle)
{
  Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder};
  rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
  rig.params[0].LayerxX0[0] = 0.005f; // inner
  rig.params[0].LayerxX0[1] = 0.005f; // middle
  rig.params[0].LayerxX0[2] = 0.f;    // outer: contractually unused by Cylinder
  rig.establishLayout();

  const std::array<o2::its::Cluster, 3> clusters{makeGlobalCluster(3.0f, 0.100f, 0.9f, 0),
                                                 makeGlobalCluster(4.0f, 0.150f, 1.05f, 0),
                                                 makeGlobalCluster(5.0f, 0.201f, 1.25f, 0)};
  loadCandidateClusters(rig, clusters,
                        {makeBarrelHit(3.f, 0.f, 0.100f, 0.9f),
                         makeBarrelHit(4.f, 0.f, 0.150f, 1.05f),
                         makeBarrelHit(5.f, 0.f, 0.201f, 1.25f)});

  auto view = TrackerTestAccess::prepare(rig.tracker, rig.frame, 0);

  const auto topology = rig.frame.getGraph(0).getView();
  const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
  BOOST_REQUIRE_GE(cellTopologyId, 0);

  injectCandidateTracklets(rig, cellTopologyId, clusters);

  // Any other cellTopologyId keeps its empty-link early-continue
  // path: cleared once up front, never touched again.
  int otherCellTopologyId = -1;
  for (int i = 0; i < topology.nCells; ++i) {
    if (i != cellTopologyId) {
      otherCellTopologyId = i;
      break;
    }
  }
  BOOST_REQUIRE_GE(otherCellTopologyId, 0);

  TrackerTestAccess::computeCells(rig.traits, view);

  BOOST_CHECK(rig.tf->getCells()[otherCellTopologyId].empty());
  BOOST_CHECK(rig.tf->getCellsLookupTable()[otherCellTopologyId].empty());
  BOOST_CHECK(rig.tf->getCellsLabel(otherCellTopologyId).empty());

  BOOST_REQUIRE_EQUAL(rig.tf->getCells()[cellTopologyId].size(), 1u);
  const auto& producedCell = rig.tf->getCells()[cellTopologyId][0];
  BOOST_CHECK(producedCell.tripletFactor().isValid());
  for (int slot = 0; slot < 3; ++slot) {
    const auto reference = producedCell.getClusterReference(slot);
    BOOST_CHECK_EQUAL(reference.surfacePosition, slot);
    BOOST_CHECK_EQUAL(reference.clusterIndex, producedCell.getClusters()[slot]);
  }

  BOOST_REQUIRE_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId].size(), 2u);
  BOOST_CHECK_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId][0], 0);
  BOOST_CHECK_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId][1], 1);

  // hasMCinformation() is false (no labels were loaded), so label
  // construction is skipped, exactly as before this change.
  BOOST_CHECK(rig.tf->getCellsLabel(cellTopologyId).empty());

  // Oracle: independently refetch the same measurements through
  // TrackerTraits::getLayerMeasurements() (never re-typed literals) and use
  // the same generic seed construction contract as production.
  const auto layerMeasurements = gsl::span<const gsl::span<const SurfaceMeasurement>>{view.workspace.layerMeasurements};
  const auto layerGlobalMeasurements = gsl::span<const gsl::span<const GlobalMeasurement>>{view.workspace.layerGlobalMeasurements};
  const auto& oracleGlobalInner = layerGlobalMeasurements[0][producedCell.getFirstClusterIndex()];
  const auto& oracleGlobalMiddle = layerGlobalMeasurements[1][producedCell.getSecondClusterIndex()];
  const auto& oracleMeasurementInner = layerMeasurements[0][producedCell.getFirstClusterIndex()];
  const auto& oracleMeasurementMiddle = layerMeasurements[1][producedCell.getSecondClusterIndex()];
  const auto& oracleMeasurementOuter = layerMeasurements[2][producedCell.getThirdClusterIndex()];
  const std::array<float, 3> xOverX0{rig.params[0].LayerxX0[0], rig.params[0].LayerxX0[1], rig.params[0].LayerxX0[2]};
  const auto material = toMaterial(xOverX0);

  TrackingKernelParameters trackingParams;
  trackingParams.maxChi2ClusterAttachment = rig.params[0].MaxChi2ClusterAttachment;

  SurfaceKinematicState oracleState{};
  float oracleChi2 = 0.f;
  OperationFailureReason oracleReason{};
  BOOST_REQUIRE(buildCellSeed(
    SurfaceKind::Cylinder, oracleGlobalInner, oracleGlobalMiddle, {}, oracleMeasurementInner, oracleMeasurementMiddle, oracleMeasurementOuter,
    material, Bz, kCompatibilityAbsCharge, kCompatibilityPID, oracleState, oracleChi2, trackingParams, oracleReason));

  checkSurfaceKinematicStateEqual(producedCell.state(), oracleState);
  BOOST_CHECK_EQUAL(producedCell.getChi2(), oracleChi2);
}

BOOST_AUTO_TEST_CASE(CylinderCellCombinationUsesTrackletMinPtScattering)
{
  auto acceptedCells = [](float trackletMinPt) {
    Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder};
    rig.params[0].TrackletMinPt = trackletMinPt;
    rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
    rig.params[0].LayerxX0[0] = 0.005f;
    rig.params[0].LayerxX0[1] = 0.01f;
    rig.params[0].LayerxX0[2] = 0.f;
    rig.establishLayout();

    const std::array<o2::its::Cluster, 3> clusters{
      makeGlobalCluster(3.f, 0.100f, 0.9f),
      makeGlobalCluster(4.f, 0.150f, 1.05f),
      makeGlobalCluster(5.f, 0.201f, 1.22f)};
    loadCandidateClusters(rig, clusters,
                          {makeBarrelHit(3.f, 0.f, 0.100f, 0.9f, 1.e-6f, 1.e-6f),
                           makeBarrelHit(4.f, 0.f, 0.150f, 1.05f, 1.e-6f, 1.e-6f),
                           makeBarrelHit(5.f, 0.f, 0.201f, 1.22f, 1.e-6f, 1.e-6f)});

    auto view = prepare(rig);
    const auto topology = rig.frame.getGraph(0).getView();
    const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
    BOOST_REQUIRE_GE(cellTopologyId, 0);
    injectCandidateTracklets(rig, cellTopologyId, clusters);
    TrackerTestAccess::computeCells(rig.traits, view);
    return rig.tf->getCells()[cellTopologyId].size();
  };

  BOOST_CHECK_EQUAL(acceptedCells(0.3f), 1u);
  BOOST_CHECK_EQUAL(acceptedCells(1.f), 0u);
}

// --- Disk: real orchestration matches the generic cell-seed oracle -------

BOOST_AUTO_TEST_CASE(DiskComputeLayerCellsMatchesBuildCellSeedOracle)
{
  Rig<MFTNLayers> rig{o2::detectors::DetID::MFT, SurfaceKind::Disk};
  rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
  rig.params[0].TrackletMinPt = 0.3f;
  rig.params[0].LayerxX0[0] = 0.015f; // inner
  rig.params[0].LayerxX0[1] = 0.017f; // middle
  rig.params[0].LayerxX0[2] = 0.02f;  // outer
  rig.establishLayout();

  const std::array<o2::its::Cluster, 3> clusters{makeGlobalCluster(1.0f, 0.5f, -0.4f, 0),
                                                 makeGlobalCluster(1.3f, 0.62f, -0.6f, 0),
                                                 makeGlobalCluster(1.7f, 0.78f, -0.9f, 0)};
  loadCandidateClusters(rig, clusters,
                        {makeDiskHit(-0.4f, 1.0f, 0.5f),
                         makeDiskHit(-0.6f, 1.3f, 0.62f),
                         makeDiskHit(-0.9f, 1.7f, 0.78f)});

  auto view = prepare(rig);

  const auto topology = rig.frame.getGraph(0).getView();
  const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
  BOOST_REQUIRE_GE(cellTopologyId, 0);

  injectCandidateTracklets(rig, cellTopologyId, clusters);

  TrackerTestAccess::computeCells(rig.traits, view);

  BOOST_REQUIRE_EQUAL(rig.tf->getCells()[cellTopologyId].size(), 1u);
  const auto& producedCell = rig.tf->getCells()[cellTopologyId][0];
  BOOST_CHECK(producedCell.tripletFactor().isValid());
  for (int slot = 0; slot < 3; ++slot) {
    const auto reference = producedCell.getClusterReference(slot);
    BOOST_CHECK_EQUAL(reference.surfacePosition, slot);
    BOOST_CHECK_EQUAL(reference.clusterIndex, producedCell.getClusters()[slot]);
  }

  const auto layerMeasurements = gsl::span<const gsl::span<const SurfaceMeasurement>>{view.workspace.layerMeasurements};
  const auto& oracleMeasurementInner = layerMeasurements[0][producedCell.getFirstClusterIndex()];
  const auto& oracleMeasurementMiddle = layerMeasurements[1][producedCell.getSecondClusterIndex()];
  const auto& oracleMeasurementOuter = layerMeasurements[2][producedCell.getThirdClusterIndex()];
  const std::array<float, 3> xOverX0{rig.params[0].LayerxX0[0], rig.params[0].LayerxX0[1], rig.params[0].LayerxX0[2]};
  const auto material = toMaterial(xOverX0);

  TrackingKernelParameters trackingParams;
  trackingParams.maxChi2ClusterAttachment = rig.params[0].MaxChi2ClusterAttachment;
  trackingParams.trackletMinPt = rig.params[0].TrackletMinPt;

  SurfaceKinematicState oracleState{};
  float oracleChi2 = 0.f;
  OperationFailureReason oracleReason{};
  BOOST_REQUIRE(buildCellSeed(
    SurfaceKind::Disk, {}, {}, {}, oracleMeasurementInner, oracleMeasurementMiddle, oracleMeasurementOuter,
    material, Bz, kCompatibilityAbsCharge, kCompatibilityPID, oracleState, oracleChi2, trackingParams, oracleReason));

  checkSurfaceKinematicStateEqual(producedCell.state(), oracleState);
  BOOST_CHECK_EQUAL(producedCell.getChi2(), oracleChi2);
}

// --- One-pass vs two-pass: identical result regardless of thread count ----

BOOST_AUTO_TEST_CASE(CylinderComputeLayerCellsOnePassAndTwoPassAgree)
{
  struct Result {
    int cellTopologyId{-1};
    std::vector<int> lut;
    float chi2{0.f};
    int cl0{-1}, cl1{-1}, cl2{-1};
  };

  auto run = [](int nThreads) {
    Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, nThreads};
    rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
    rig.params[0].LayerxX0[0] = 0.005f;
    rig.params[0].LayerxX0[1] = 0.005f;
    rig.establishLayout();

    const std::array<o2::its::Cluster, 3> clusters{makeGlobalCluster(3.0f, 0.100f, 0.9f, 0),
                                                   makeGlobalCluster(4.0f, 0.150f, 1.05f, 0),
                                                   makeGlobalCluster(5.0f, 0.201f, 1.25f, 0)};
    loadCandidateClusters(rig, clusters,
                          {makeBarrelHit(3.f, 0.f, 0.100f, 0.9f),
                           makeBarrelHit(4.f, 0.f, 0.150f, 1.05f),
                           makeBarrelHit(5.f, 0.f, 0.201f, 1.25f)});

    auto view = prepare(rig);

    const auto topology = rig.frame.getGraph(0).getView();
    const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
    BOOST_REQUIRE_GE(cellTopologyId, 0);

    injectCandidateTracklets(rig, cellTopologyId, clusters);

    TrackerTestAccess::computeCells(rig.traits, view);

    Result r;
    r.cellTopologyId = cellTopologyId;
    const auto& lut = rig.tf->getCellsLookupTable()[cellTopologyId];
    r.lut.assign(lut.begin(), lut.end());
    BOOST_REQUIRE_EQUAL(rig.tf->getCells()[cellTopologyId].size(), 1u);
    const auto& cell = rig.tf->getCells()[cellTopologyId][0];
    r.chi2 = cell.getChi2();
    r.cl0 = cell.getFirstClusterIndex();
    r.cl1 = cell.getSecondClusterIndex();
    r.cl2 = cell.getThirdClusterIndex();
    return r;
  };

  const auto onePass = run(1);
  const auto twoPass = run(4);

  BOOST_CHECK_EQUAL(onePass.cellTopologyId, twoPass.cellTopologyId);
  BOOST_CHECK_EQUAL_COLLECTIONS(onePass.lut.begin(), onePass.lut.end(), twoPass.lut.begin(), twoPass.lut.end());
  BOOST_CHECK_EQUAL(onePass.chi2, twoPass.chi2);
  BOOST_CHECK_EQUAL(onePass.cl0, twoPass.cl0);
  BOOST_CHECK_EQUAL(onePass.cl1, twoPass.cl1);
  BOOST_CHECK_EQUAL(onePass.cl2, twoPass.cl2);
}

BOOST_AUTO_TEST_CASE(DiskComputeLayerCellsOnePassAndTwoPassAgree)
{
  struct Result {
    int cellTopologyId{-1};
    std::vector<int> lut;
    float chi2{0.f};
    int cl0{-1}, cl1{-1}, cl2{-1};
  };

  auto run = [](int nThreads) {
    Rig<MFTNLayers> rig{o2::detectors::DetID::MFT, SurfaceKind::Disk, nThreads};
    rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
    rig.params[0].TrackletMinPt = 0.3f;
    rig.establishLayout();

    const std::array<o2::its::Cluster, 3> clusters{makeGlobalCluster(1.0f, 0.5f, -0.4f, 0),
                                                   makeGlobalCluster(1.3f, 0.62f, -0.6f, 0),
                                                   makeGlobalCluster(1.7f, 0.78f, -0.9f, 0)};
    loadCandidateClusters(rig, clusters,
                          {makeDiskHit(-0.4f, 1.0f, 0.5f),
                           makeDiskHit(-0.6f, 1.3f, 0.62f),
                           makeDiskHit(-0.9f, 1.7f, 0.78f)});

    auto view = prepare(rig);

    const auto topology = rig.frame.getGraph(0).getView();
    const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
    BOOST_REQUIRE_GE(cellTopologyId, 0);

    injectCandidateTracklets(rig, cellTopologyId, clusters);

    TrackerTestAccess::computeCells(rig.traits, view);

    Result r;
    r.cellTopologyId = cellTopologyId;
    const auto& lut = rig.tf->getCellsLookupTable()[cellTopologyId];
    r.lut.assign(lut.begin(), lut.end());
    BOOST_REQUIRE_EQUAL(rig.tf->getCells()[cellTopologyId].size(), 1u);
    const auto& cell = rig.tf->getCells()[cellTopologyId][0];
    r.chi2 = cell.getChi2();
    r.cl0 = cell.getFirstClusterIndex();
    r.cl1 = cell.getSecondClusterIndex();
    r.cl2 = cell.getThirdClusterIndex();
    return r;
  };

  const auto onePass = run(1);
  const auto twoPass = run(4);

  BOOST_CHECK_EQUAL(onePass.cellTopologyId, twoPass.cellTopologyId);
  BOOST_CHECK_EQUAL_COLLECTIONS(onePass.lut.begin(), onePass.lut.end(), twoPass.lut.begin(), twoPass.lut.end());
  BOOST_CHECK_EQUAL(onePass.chi2, twoPass.chi2);
  BOOST_CHECK_EQUAL(onePass.cl0, twoPass.cl0);
  BOOST_CHECK_EQUAL(onePass.cl1, twoPass.cl1);
  BOOST_CHECK_EQUAL(onePass.cl2, twoPass.cl2);
}

BOOST_AUTO_TEST_CASE(RepeatedComputeLayerCellsCallsDoNotRebindOrIncreaseCounts)
{
  Rig<MFTNLayers> rig{o2::detectors::DetID::MFT, SurfaceKind::Disk};
  rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
  rig.params[0].TrackletMinPt = 0.3f;
  rig.establishLayout();

  const std::array<o2::its::Cluster, 3> clusters{makeGlobalCluster(1.0f, 0.5f, -0.4f, 0),
                                                 makeGlobalCluster(1.3f, 0.62f, -0.6f, 0),
                                                 makeGlobalCluster(1.7f, 0.78f, -0.9f, 0)};
  loadCandidateClusters(rig, clusters,
                        {makeDiskHit(-0.4f, 1.0f, 0.5f),
                         makeDiskHit(-0.6f, 1.3f, 0.62f),
                         makeDiskHit(-0.9f, 1.7f, 0.78f)});

  auto view = prepare(rig);

  const auto topology = rig.frame.getGraph(0).getView();
  const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
  BOOST_REQUIRE_GE(cellTopologyId, 0);

  injectCandidateTracklets(rig, cellTopologyId, clusters);

  TrackerTestAccess::computeCells(rig.traits, view);
  BOOST_REQUIRE_EQUAL(rig.tf->getCells()[cellTopologyId].size(), 1u);
  const auto firstChi2 = rig.tf->getCells()[cellTopologyId][0].getChi2();

  TrackerTestAccess::computeCells(rig.traits, view);
  TrackerTestAccess::computeCells(rig.traits, view);

  // Re-inject tracklets and recompute (the underlying candidate clusters/
  // measurements loaded above are untouched -- reloading them here would
  // invalidate TrackerTraits::mLayerMeasurements without a fresh
  // initialiseTimeFrame() call to re-resolve it, which is not what this test
  // checks): a fresh call after the tracklets were consumed must still
  // reproduce the identical chi2 through the same cache.
  injectCandidateTracklets(rig, cellTopologyId, clusters);
  TrackerTestAccess::computeCells(rig.traits, view);

  BOOST_REQUIRE_EQUAL(rig.tf->getCells()[cellTopologyId].size(), 1u);
  BOOST_CHECK_EQUAL(rig.tf->getCells()[cellTopologyId][0].getChi2(), firstChi2);
}

// Material-correction preflight has its own focused test target; this file
// covers only direct cell-stage orchestration.

BOOST_AUTO_TEST_CASE(CylinderComputeLayerCellsMultiCellChainProducesCorrectCellsAndOrder)
{
  // 5-layer chain at global X = 3..7 (small Y, alpha=0.f, matching the
  // single-cell oracle tests' convention above): proves link-level/
  // cell-level parity across three simultaneously-populated cells (0,1,2),
  // (1,2,3), (2,3,4) -- each resolved through the migrated
  // computeLayerCellsForKind() via a fresh mSurfaceToLegacyLayer lookup
  // per CellTopologyId -- not just the single cell the tests above check,
  // while every non-participating cellTopologyId stays empty.
  Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder};
  rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
  for (int layer = 0; layer < ITSNLayers; ++layer) {
    rig.params[0].LayerxX0[layer] = 0.005f;
  }
  rig.establishLayout();

  // Y values lie exactly on one real circle (center (0,50), radius 50,
  // through the origin) rather than an ad hoc linear Y(X): a single
  // physically consistent curvature across all 5 points avoids the
  // rotation-boundary edge cases (BarrelSurfaceStateOperations.cxx's
  // csp*ca+snp*sa<0 checks) an inconsistent, near-degenerate linear Y(X)
  // can trip for some sub-triples but not others.
  constexpr std::array<int, 5> layers{0, 1, 2, 3, 4};
  constexpr std::array<float, 5> xs{3.f, 4.f, 5.f, 6.f, 7.f};
  constexpr std::array<float, 5> ys{0.0902f, 0.1603f, 0.2506f, 0.3614f, 0.4924f};
  constexpr std::array<float, 5> zs{0.90f, 1.05f, 1.20f, 1.35f, 1.50f};
  std::array<o2::its::Cluster, 5> clusters;
  std::array<o2::its::TrackingFrameInfo, 5> hits;
  for (size_t i = 0; i < 5; ++i) {
    clusters[i] = makeGlobalCluster(xs[i], ys[i], zs[i], 0);
    hits[i] = makeBarrelHit(xs[i], 0.f, ys[i], zs[i]);
  }
  loadCandidateClustersAtLayers(rig, layers, clusters, hits);

  auto view = prepare(rig);

  const auto topology = rig.frame.getGraph(0).getView();
  injectChainCandidateTracklets(rig, layers, clusters);

  TrackerTestAccess::computeCells(rig.traits, view);

  const std::array<std::array<int, 3>, 3> triples{{{0, 1, 2}, {1, 2, 3}, {2, 3, 4}}};
  std::array<int, 3> topologyIds{};
  std::vector<bool> participating(topology.nCells, false);
  for (size_t i = 0; i < triples.size(); ++i) {
    const auto& triple = triples[i];
    const int cellTopologyId = findCellTopologyId(topology, triple[0], triple[1], triple[2]);
    BOOST_REQUIRE_GE(cellTopologyId, 0);
    topologyIds[i] = cellTopologyId;
    participating[cellTopologyId] = true;

    BOOST_REQUIRE_EQUAL(rig.tf->getCells()[cellTopologyId].size(), 1u);
    const auto& producedCell = rig.tf->getCells()[cellTopologyId][0];
    BOOST_CHECK_EQUAL(producedCell.getFirstClusterIndex(), 0);
    BOOST_CHECK_EQUAL(producedCell.getSecondClusterIndex(), 0);
    BOOST_CHECK_EQUAL(producedCell.getThirdClusterIndex(), 0);
    BOOST_CHECK_EQUAL(producedCell.getHitLayerMask().value(), LayerMask(triple[0], triple[1], triple[2]).value());

    BOOST_REQUIRE_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId].size(), 2u);
    BOOST_CHECK_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId][0], 0);
    BOOST_CHECK_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId][1], 1);
  }

  for (int i = 0; i < topology.nCells; ++i) {
    if (!participating[i]) {
      BOOST_CHECK(rig.tf->getCells()[i].empty());
    }
  }

  TrackerTestAccess::findNeighbours(rig.traits, view);
  for (size_t i = 0; i < topologyIds.size(); ++i) {
    BOOST_CHECK_EQUAL(rig.tf->getCells()[topologyIds[i]][0].getLevel(), static_cast<int>(i + 1));
    if (i == 0) {
      BOOST_CHECK(rig.tf->getCellsNeighbours()[topologyIds[i]].empty());
      continue;
    }
    BOOST_REQUIRE_EQUAL(rig.tf->getCellsNeighbours()[topologyIds[i]].size(), 1u);
    BOOST_CHECK_EQUAL(rig.tf->getCellsNeighbours()[topologyIds[i]][0], 0);
    BOOST_CHECK_EQUAL(rig.tf->getCellsNeighboursTopology()[topologyIds[i]][0], topologyIds[i - 1]);
  }
}

BOOST_AUTO_TEST_CASE(DiskComputeLayerCellsMultiCellChainProducesCorrectCellsAndOrder)
{
  // Same multi-cell parity property for the Disk/forward family:
  // cell-seed leaves genuinely branches per family (Cylinder
  // reads [1] then [0]; Disk reads [2],[1],[0] -- see the comment on
  // that call in computeLayerCellsForKind()), so multi-link
  // cell-chaining for Disk is real, otherwise-unproven coverage.
  Rig<MFTNLayers> rig{o2::detectors::DetID::MFT, SurfaceKind::Disk};
  rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
  rig.params[0].TrackletMinPt = 0.3f;
  rig.establishLayout();

  constexpr std::array<int, 5> layers{0, 1, 2, 3, 4};
  constexpr std::array<float, 5> xs{1.0f, 1.3f, 1.6f, 1.9f, 2.2f};
  constexpr std::array<float, 5> ys{0.50f, 0.62f, 0.74f, 0.86f, 0.98f};
  constexpr std::array<float, 5> zs{-0.40f, -0.60f, -0.80f, -1.00f, -1.20f};
  std::array<o2::its::Cluster, 5> clusters;
  std::array<o2::its::TrackingFrameInfo, 5> hits;
  for (size_t i = 0; i < 5; ++i) {
    clusters[i] = makeGlobalCluster(xs[i], ys[i], zs[i], 0);
    hits[i] = makeDiskHit(zs[i], xs[i], ys[i]);
  }
  loadCandidateClustersAtLayers(rig, layers, clusters, hits);

  auto view = prepare(rig);

  const auto topology = rig.frame.getGraph(0).getView();
  injectChainCandidateTracklets(rig, layers, clusters);

  TrackerTestAccess::computeCells(rig.traits, view);

  const std::array<std::array<int, 3>, 3> triples{{{0, 1, 2}, {1, 2, 3}, {2, 3, 4}}};
  std::array<int, 3> topologyIds{};
  std::vector<bool> participating(topology.nCells, false);
  for (size_t i = 0; i < triples.size(); ++i) {
    const auto& triple = triples[i];
    const int cellTopologyId = findCellTopologyId(topology, triple[0], triple[1], triple[2]);
    BOOST_REQUIRE_GE(cellTopologyId, 0);
    topologyIds[i] = cellTopologyId;
    participating[cellTopologyId] = true;

    BOOST_REQUIRE_EQUAL(rig.tf->getCells()[cellTopologyId].size(), 1u);
    const auto& producedCell = rig.tf->getCells()[cellTopologyId][0];
    BOOST_CHECK_EQUAL(producedCell.getHitLayerMask().value(), LayerMask(triple[0], triple[1], triple[2]).value());

    BOOST_REQUIRE_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId].size(), 2u);
    BOOST_CHECK_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId][0], 0);
    BOOST_CHECK_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId][1], 1);
  }

  for (int i = 0; i < topology.nCells; ++i) {
    if (!participating[i]) {
      BOOST_CHECK(rig.tf->getCells()[i].empty());
    }
  }

  TrackerTestAccess::findNeighbours(rig.traits, view);
  for (size_t i = 0; i < topologyIds.size(); ++i) {
    BOOST_CHECK_EQUAL(rig.tf->getCells()[topologyIds[i]][0].getLevel(), static_cast<int>(i + 1));
    if (i == 0) {
      BOOST_CHECK(rig.tf->getCellsNeighbours()[topologyIds[i]].empty());
      continue;
    }
    BOOST_REQUIRE_EQUAL(rig.tf->getCellsNeighbours()[topologyIds[i]].size(), 1u);
    BOOST_CHECK_EQUAL(rig.tf->getCellsNeighbours()[topologyIds[i]][0], 0);
    BOOST_CHECK_EQUAL(rig.tf->getCellsNeighboursTopology()[topologyIds[i]][0], topologyIds[i - 1]);
  }
}

BOOST_AUTO_TEST_CASE(CylinderComputeLayerCellsHoleCellReconstructsCorrectLayerMask)
{
  // MaxHoles=1 with layer 1 an allowed hole introduces a (0,2)-skip-1
  // link; combined with the adjacent (2,3) link this forms cell
  // (0,2,3) -- a direct, non-adjacent exercise of resolveCellHitLayers()
  // (mSurfaceToLegacyLayer) resolving a cell's endpoints correctly, and of
  // hole/skipped-surface behaviour staying identical to the pre-migration
  // code (which read the same fromLayer/toLayer straight off the legacy
  // view). No cluster is placed on layer 1 at all.
  Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder};
  rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
  rig.params[0].MaxHoles = 1;
  rig.params[0].HoleLayerMask = LayerMask{static_cast<uint16_t>(1u << 1)};
  rig.establishLayout();

  constexpr std::array<int, 3> layers{0, 2, 3};
  const std::array<o2::its::Cluster, 3> clusters{
    makeGlobalCluster(3.f, 0.10f, 0.90f, 0),
    makeGlobalCluster(5.f, 0.20f, 1.20f, 0),
    makeGlobalCluster(6.f, 0.25f, 1.35f, 0)};
  const std::array<o2::its::TrackingFrameInfo, 3> hits{
    makeBarrelHit(3.f, 0.f, 0.10f, 0.90f),
    makeBarrelHit(5.f, 0.f, 0.20f, 1.20f),
    makeBarrelHit(6.f, 0.f, 0.25f, 1.35f)};
  loadCandidateClustersAtLayers(rig, layers, clusters, hits);

  auto view = prepare(rig);

  const auto topology = rig.frame.getGraph(0).getView();
  injectChainCandidateTracklets(rig, layers, clusters);

  TrackerTestAccess::computeCells(rig.traits, view);

  const int cellTopologyId = findCellTopologyId(topology, 0, 2, 3);
  BOOST_REQUIRE_GE(cellTopologyId, 0);
  BOOST_REQUIRE_EQUAL(rig.tf->getCells()[cellTopologyId].size(), 1u);
  const auto& producedCell = rig.tf->getCells()[cellTopologyId][0];
  BOOST_CHECK_EQUAL(producedCell.getHitLayerMask().value(), LayerMask(0, 2, 3).value());

  BOOST_REQUIRE_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId].size(), 2u);
  BOOST_CHECK_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId][0], 0);
  BOOST_CHECK_EQUAL(rig.tf->getCellsLookupTable()[cellTopologyId][1], 1);

  for (int i = 0; i < topology.nCells; ++i) {
    if (i != cellTopologyId) {
      BOOST_CHECK(rig.tf->getCells()[i].empty());
    }
  }
}

BOOST_AUTO_TEST_CASE(ComputeLayerCellsFailsClosedOnCellStorageSizeMismatch)
{
  // Gate 4 Slice 0b defensive check: computeLayerCells() now compares
  // mTimeFrame->getCells().size()/getCellsLookupTable().size() against the
  // sparse cell count before its own clear loop touches either container
  // (that loop is no longer definitionally safe once its bound is
  // sparse-derived rather than mTimeFrame's own legacy-shaped storage size).
  // Desyncing a public container after a genuinely successful
  // initialiseTimeFrame() -- mirroring testTimeFrameSurfaceGraphs.cxx's
  // established legacy-cell-container-size-mismatch seam -- exercises this
  // exactly, through computeLayerCells()'s own public contract, with no
  // private-state bypass.
  Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder};
  rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
  rig.establishLayout();

  const std::array<o2::its::Cluster, 3> clusters{makeGlobalCluster(3.0f, 0.100f, 0.9f, 0),
                                                 makeGlobalCluster(4.0f, 0.150f, 1.05f, 0),
                                                 makeGlobalCluster(5.0f, 0.201f, 1.25f, 0)};
  loadCandidateClusters(rig, clusters,
                        {makeBarrelHit(3.f, 0.f, 0.100f, 0.9f),
                         makeBarrelHit(4.f, 0.f, 0.150f, 1.05f),
                         makeBarrelHit(5.f, 0.f, 0.201f, 1.25f)});

  auto view = prepare(rig);

  rig.tf->getCells().pop_back();

  BOOST_CHECK_EXCEPTION(TrackerTestAccess::computeCells(rig.traits, view), TraversalException, [](const TraversalException& error) {
    return error.getReason() == TraversalFailureReason::SparseTopologyMismatch;
  });
}
