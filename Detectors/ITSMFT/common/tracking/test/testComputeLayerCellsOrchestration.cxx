// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Orchestration coverage for TrackerTraits<NLayers>::computeLayerCells after
// its detector-family branch was replaced by a one-shot outer dispatch to
// buildCellSeed<Tag> (Architecture.md Sec 10/10.1). buildCellSeed<Tag>'s
// numerical parity with the legacy inline formulas is already proven by
// testTransitionPolicyOperations.cxx; this file does not re-derive or
// duplicate that formula. It proves instead that the real public
// computeLayerCells() entry point:
//  - resolves the three clusters/hits/material values for a candidate in
//    strict {inner, middle, outer} order and hands them to buildCellSeed<Tag>
//    unchanged (checked against an oracle call to buildCellSeed<Tag> itself,
//    built from the same values refetched through the TimeFrame, not
//    re-typed literals);
//  - keeps the MFT road pre-cut (passesCellRoadPrecut<DiskDisk>, formerly
//    detail::validateMFTCellClusters) outside buildCellSeed, called
//    unconditionally for both families with no detector-ID branch in the
//    candidate loop, and driven by the bound DiskDiskPolicyParams::cellRoadRCut
//    and the once-per-iteration cached legacy reference-z span
//    (TrackerTraits::mDiskLayerReferenceZ);
//  - leaves cellTopologyId indexing, the LUT, MC-label construction, and
//    one-pass/two-pass ordering untouched;
//  - fails closed (TraversalException::InvalidTraversalSchedule) through the
//    existing public API alone, with no test-only seam into private
//    traversal-cache state.

#define BOOST_TEST_MODULE ITSMFT ComputeLayerCells orchestration
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
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
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITSMFTTracking/TransitionPolicyOperations.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Constants.h"
#include "ITStracking/Tracklet.h"
#include "MFTTracking/Constants.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

constexpr float Bz = 0.5f;

// Preflight-only fixtures (Rig::establishLayout()) load zero real clusters --
// this decoder's decode() is never actually invoked there. It exists only to
// satisfy loadNormalizedSource()'s interface, mirroring
// testCATrackerFailureContract.cxx's LegacyLikeDecoder.
class NeverDecodedDecoder final : public ClusterDecoder
{
 public:
  explicit NeverDecodedDecoder(o2::detectors::DetID::ID detector) : mDetector(detector) {}

  o2::itsmft::ioutils::SurfaceMeasurementDecodeResult decode(
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
  FixedMeasurementDecoder(o2::detectors::DetID::ID detector, SurfaceKind kind) : mDetector(detector), mKind(kind) {}

  void setMeasurement(int layer, const SurfaceMeasurement& measurement) { mByLayer[layer] = measurement; }

  o2::itsmft::ioutils::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor&,
    const TopologyDictionary*,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool) const override
  {
    o2::itsmft::ioutils::SurfaceMeasurementDecodeResult result;
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
    measurement.surface = layerToSurface[layer];
    measurement.sensor = DetectorSensorId{static_cast<uint32_t>(mDetector), static_cast<uint32_t>(layer)};
    measurement.cluster = ClusterRef{source, externalIndex};
    measurement.sourceROF = sourceROF;
    result.measurement = measurement;
    return result;
  }

 private:
  o2::detectors::DetID::ID mDetector;
  SurfaceKind mKind;
  std::map<int, SurfaceMeasurement> mByLayer;
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

class FakeCatalogProvider final : public DetectorSurfaceCatalogProvider
{
 public:
  explicit FakeCatalogProvider(std::vector<SurfaceDescriptor> catalog) : mCatalog{std::move(catalog)} {}

  DetectorSurfaceCatalogResult buildCatalog(const DetectorSurfaceCatalogRequest&) const final
  {
    return {mCatalog, DetectorSurfaceCatalogError::None};
  }

  std::vector<SurfaceDescriptor> mCatalog;
};

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

// Same construction as testTransitionPolicyOperations.cxx's helpers -- plain
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
// testTransitionPolicyOperationsNative.cxx: the single SurfaceMeasurement now
// standing in for the retired {Cluster, TrackingFrameInfo} pair at each
// candidate position.
SurfaceMeasurement barrelMeasurementFor(const o2::its::Cluster& cluster, const o2::its::TrackingFrameInfo& hit)
{
  SurfaceMeasurement measurement{};
  measurement.global = {cluster.xCoordinate, cluster.yCoordinate, cluster.zCoordinate};
  measurement.frame.q = hit.xTrackingFrame;
  measurement.frame.frameAngle = hit.alphaTrackingFrame;
  measurement.frame.u = hit.positionTrackingFrame[0];
  measurement.frame.v = hit.positionTrackingFrame[1];
  measurement.covariance.uu = hit.covarianceTrackingFrame[0];
  measurement.covariance.uv = hit.covarianceTrackingFrame[1];
  measurement.covariance.vv = hit.covarianceTrackingFrame[2];
  return measurement;
}

SurfaceMeasurement diskMeasurementFor(const o2::its::Cluster& cluster, const o2::its::TrackingFrameInfo& hit)
{
  SurfaceMeasurement measurement{};
  measurement.global = {cluster.xCoordinate, cluster.yCoordinate, cluster.zCoordinate};
  measurement.frame.q = cluster.zCoordinate; // disk adapter contract: frame.q == global.z
  measurement.covariance.uu = hit.covarianceTrackingFrame[0];
  measurement.covariance.uv = 0.f;
  measurement.covariance.vv = hit.covarianceTrackingFrame[2];
  return measurement;
}

// Stage-B activation: the produced Cell no longer inherits a track
// parametrization, so the oracle comparison is done directly on
// SurfaceKinematicState (both the produced cell's own .state() and the
// oracle's native buildCellSeed<Tag> output are this type).
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
// mTraversalGrouping/mCylinderPolicyParams|mDiskPolicyParams/mAttachHitConfig
// -- computeLayerCells()'s own private caches, never poked directly), and a
// validly-sized-but-empty normalized load (proven pattern from
// testCATrackerFailureContract.cxx: TimeFrame::initialise() unconditionally
// reads mROFramesClusters sizes, which only loadNormalizedSource() sets up
// safely, even for zero clusters).
template <int NLayers>
struct Rig {
  Rig(o2::detectors::DetID::ID det, SurfaceKind kind, TransitionPolicyTag tag, int nThreads = 1)
    : pool(std::make_shared<BoundedMemoryResource>()),
      params(1),
      mDet(det),
      mKind(kind),
      mTag(tag)
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
    tf.setMemoryPool(pool);
    traits.setMemoryPool(pool);
    traits.setNThreads(nThreads, arena);
    traits.adoptTimeFrame(&tf);
    traits.updateTrackingParameters(params);
    traits.setBz(Bz);
  }

  // Establishes the catalog/layout and loads a (zero-cluster) normalized
  // source. Deliberately not run by the constructor: it builds the catalog's
  // nominal material from the *current* params[0].LayerxX0, so callers must
  // finish any TrackingParameters::LayerxX0 override first -- matching
  // TrackerTraits::initialiseTimeFrame()'s LegacyMaterialMismatch
  // compatibility check, which this file does not itself exercise.
  void establishLayout()
  {
    const auto orderedSurfaces = identitySurfaces(static_cast<uint16_t>(NLayers));
    FakeCatalogProvider provider{makeCatalog(static_cast<uint16_t>(NLayers), mDet, mKind,
                                             gsl::span<const float>(params[0].LayerxX0))};
    const DetectorSurfaceCatalogRequest catalogRequest{mDet, SurfaceId{0}, static_cast<uint32_t>(NLayers)};
    BOOST_REQUIRE(tf.ensureDetectorLayouts(&provider, catalogRequest, orderedSurfaces, mTag, params).ok());
    tf.initTrackerTopologies(params);

    NeverDecodedDecoder decoder{mDet};
    const o2::InteractionRecord origin{50, 5};
    const ROFTimingConfig timing{40, 0, 0, 0};
    const std::vector<CompClusterExt> noClusters;
    const std::vector<unsigned char> noPatterns;
    const std::vector<ROFRecord> noRofs;
    const auto result = tf.loadNormalizedSource(decoder, origin, timing, noClusters, noPatterns, noRofs, &dict(), nullptr, mDet);
    BOOST_REQUIRE(result.ok());
  }

  o2::detectors::DetID::ID detector() const noexcept { return mDet; }
  SurfaceKind kind() const noexcept { return mKind; }

  std::shared_ptr<BoundedMemoryResource> pool;
  std::vector<TrackingParameters> params;
  TimeFrame<NLayers> tf;
  TrackerTraits<NLayers> traits;
  std::shared_ptr<tbb::task_arena> arena;

 private:
  o2::detectors::DetID::ID mDet;
  SurfaceKind mKind;
  TransitionPolicyTag mTag;
};

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
  const auto result = rig.tf.loadNormalizedSource(decoder, origin, timing, compClusters, noPatterns, rofs, &dict(), nullptr, rig.detector());
  BOOST_REQUIRE(result.ok());
}

// Finds the cellTopologyId whose two transitions span exactly
// inner->middle->outer, without assuming any particular enumeration order
// out of TrackingTopology<NLayers>::init().
template <typename TopologyView>
int findCellTopologyId(const TopologyView& topology, int inner, int middle, int outer)
{
  for (int i = 0; i < topology.nCells; ++i) {
    const auto& cell = topology.getCell(i);
    const auto& first = topology.getTransition(cell.firstTransition);
    const auto& second = topology.getTransition(cell.secondTransition);
    if (first.fromLayer == inner && first.toLayer == middle && second.fromLayer == middle && second.toLayer == outer) {
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
// tracklet per transition of cellTopologyId, wired so
// computeLayerCellsForPolicy<Tag>'s tracklet-pairing loop finds exactly one
// candidate pair.
template <int NLayers>
void injectCandidateTracklets(Rig<NLayers>& rig, int cellTopologyId, const std::array<o2::its::Cluster, 3>& clusters)
{
  const auto topology = rig.tf.getTrackingTopologyView();
  const auto& cell = topology.getCell(cellTopologyId);
  const auto& first = topology.getTransition(cell.firstTransition);
  const auto& second = topology.getTransition(cell.secondTransition);
  const int layers[3] = {first.fromLayer, first.toLayer, second.toLayer};

  for (int i = 0; i < 3; ++i) {
    BOOST_REQUIRE_EQUAL(rig.tf.getClusters()[layers[i]].size(), 1u);
    rig.tf.getClusters()[layers[i]][0] = clusters[i];
  }

  const o2::its::TimeEstBC ts{static_cast<uint32_t>(0), static_cast<uint16_t>(1)};
  rig.tf.getTracklets()[cell.firstTransition].push_back(o2::its::Tracklet{0, 0, 0.f, 0.f, ts});
  rig.tf.getTracklets()[cell.secondTransition].push_back(o2::its::Tracklet{0, 0, 0.f, 0.f, ts});

  auto& secondLUT = rig.tf.getTrackletsLookupTable()[cell.secondTransition];
  secondLUT.resize(2);
  secondLUT[0] = 0;
  secondLUT[1] = 1;
}

} // namespace

// --- Barrel: real orchestration matches the buildCellSeed<Tag> oracle -----

BOOST_AUTO_TEST_CASE(CylinderComputeLayerCellsMatchesBuildCellSeedOracle)
{
  Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder};
  rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
  rig.params[0].LayerxX0[0] = 0.005f; // inner
  rig.params[0].LayerxX0[1] = 0.005f; // middle
  rig.params[0].LayerxX0[2] = 0.f;    // outer: contractually unused by CylinderCylinder
  rig.establishLayout();

  const std::array<o2::its::Cluster, 3> clusters{makeGlobalCluster(3.0f, 0.100f, 0.9f, 0),
                                                 makeGlobalCluster(4.0f, 0.150f, 1.05f, 0),
                                                 makeGlobalCluster(5.0f, 0.201f, 1.25f, 0)};
  loadCandidateClusters(rig, clusters,
                        {makeBarrelHit(3.f, 0.f, 0.100f, 0.9f),
                         makeBarrelHit(4.f, 0.f, 0.150f, 1.05f),
                         makeBarrelHit(5.f, 0.f, 0.201f, 1.25f)});

  rig.traits.updateTrackingParameters(rig.params);
  rig.traits.initialiseTimeFrame(0);
  BOOST_REQUIRE(rig.traits.hasTraversalCache());

  const auto topology = rig.tf.getTrackingTopologyView();
  const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
  BOOST_REQUIRE_GE(cellTopologyId, 0);

  injectCandidateTracklets(rig, cellTopologyId, clusters);

  // Any other cellTopologyId keeps its empty-transition early-continue
  // path: cleared once up front, never touched again.
  int otherCellTopologyId = -1;
  for (int i = 0; i < topology.nCells; ++i) {
    if (i != cellTopologyId) {
      otherCellTopologyId = i;
      break;
    }
  }
  BOOST_REQUIRE_GE(otherCellTopologyId, 0);

  rig.traits.computeLayerCells(0);

  BOOST_CHECK(rig.tf.getCells()[otherCellTopologyId].empty());
  BOOST_CHECK(rig.tf.getCellsLookupTable()[otherCellTopologyId].empty());
  BOOST_CHECK(rig.tf.getCellsLabel(otherCellTopologyId).empty());

  BOOST_REQUIRE_EQUAL(rig.tf.getCells()[cellTopologyId].size(), 1u);
  const auto& producedCell = rig.tf.getCells()[cellTopologyId][0];

  BOOST_REQUIRE_EQUAL(rig.tf.getCellsLookupTable()[cellTopologyId].size(), 2u);
  BOOST_CHECK_EQUAL(rig.tf.getCellsLookupTable()[cellTopologyId][0], 0);
  BOOST_CHECK_EQUAL(rig.tf.getCellsLookupTable()[cellTopologyId][1], 1);

  // hasMCinformation() is false (no labels were loaded), so label
  // construction is skipped, exactly as before this change.
  BOOST_CHECK(rig.tf.getCellsLabel(cellTopologyId).empty());

  // Oracle: independently refetch the same measurements through
  // TrackerTraits::getLayerMeasurements() (never re-typed literals) and call
  // buildCellSeed<CylinderCylinder> directly.
  const auto layerMeasurements = rig.traits.getLayerMeasurements();
  const auto& oracleMeasurementInner = layerMeasurements[0][producedCell.getFirstClusterIndex()];
  const auto& oracleMeasurementMiddle = layerMeasurements[1][producedCell.getSecondClusterIndex()];
  const auto& oracleMeasurementOuter = layerMeasurements[2][producedCell.getThirdClusterIndex()];
  const std::array<float, 3> xOverX0{rig.params[0].LayerxX0[0], rig.params[0].LayerxX0[1], rig.params[0].LayerxX0[2]};
  const auto material = toMaterial(xOverX0);

  CylinderCylinderPolicyParams policyParams;
  policyParams.maxChi2ClusterAttachment = rig.params[0].MaxChi2ClusterAttachment;

  SurfaceKinematicState oracleState{};
  float oracleChi2 = 0.f;
  OperationFailureReason oracleReason{};
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::CylinderCylinder>(
    oracleMeasurementInner, oracleMeasurementMiddle, oracleMeasurementOuter,
    material, Bz, kCompatibilityAbsCharge, kCompatibilityPID, oracleState, oracleChi2, policyParams, oracleReason));

  checkSurfaceKinematicStateEqual(producedCell.state(), oracleState);
  BOOST_CHECK_EQUAL(producedCell.getChi2(), oracleChi2);
}

// --- Disk: real orchestration matches the buildCellSeed<Tag> oracle -------

BOOST_AUTO_TEST_CASE(DiskComputeLayerCellsMatchesBuildCellSeedOracle)
{
  Rig<MFTNLayers> rig{o2::detectors::DetID::MFT, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk};
  rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
  rig.params[0].TrackletMinPt = 0.3f;
  rig.params[0].CellRoadRCut = 1000.f; // generous: road pre-cut must pass here
  rig.params[0].LayerxX0[0] = 0.015f;  // inner
  rig.params[0].LayerxX0[1] = 0.017f;  // middle
  rig.params[0].LayerxX0[2] = 0.02f;   // outer
  rig.establishLayout();

  const std::array<o2::its::Cluster, 3> clusters{makeGlobalCluster(1.0f, 0.5f, -0.4f, 0),
                                                 makeGlobalCluster(1.3f, 0.62f, -0.6f, 0),
                                                 makeGlobalCluster(1.7f, 0.78f, -0.9f, 0)};
  loadCandidateClusters(rig, clusters,
                        {makeDiskHit(-0.4f, 1.0f, 0.5f),
                         makeDiskHit(-0.6f, 1.3f, 0.62f),
                         makeDiskHit(-0.9f, 1.7f, 0.78f)});

  rig.traits.updateTrackingParameters(rig.params);
  rig.traits.initialiseTimeFrame(0);
  BOOST_REQUIRE(rig.traits.hasTraversalCache());

  const auto topology = rig.tf.getTrackingTopologyView();
  const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
  BOOST_REQUIRE_GE(cellTopologyId, 0);

  injectCandidateTracklets(rig, cellTopologyId, clusters);

  rig.traits.computeLayerCells(0);

  BOOST_REQUIRE_EQUAL(rig.tf.getCells()[cellTopologyId].size(), 1u);
  const auto& producedCell = rig.tf.getCells()[cellTopologyId][0];

  const auto layerMeasurements = rig.traits.getLayerMeasurements();
  const auto& oracleMeasurementInner = layerMeasurements[0][producedCell.getFirstClusterIndex()];
  const auto& oracleMeasurementMiddle = layerMeasurements[1][producedCell.getSecondClusterIndex()];
  const auto& oracleMeasurementOuter = layerMeasurements[2][producedCell.getThirdClusterIndex()];
  const std::array<float, 3> xOverX0{rig.params[0].LayerxX0[0], rig.params[0].LayerxX0[1], rig.params[0].LayerxX0[2]};
  const auto material = toMaterial(xOverX0);

  DiskDiskPolicyParams policyParams;
  policyParams.maxChi2ClusterAttachment = rig.params[0].MaxChi2ClusterAttachment;
  policyParams.trackletMinPt = rig.params[0].TrackletMinPt;

  SurfaceKinematicState oracleState{};
  float oracleChi2 = 0.f;
  OperationFailureReason oracleReason{};
  BOOST_REQUIRE(buildCellSeed<TransitionPolicyTag::DiskDisk>(
    oracleMeasurementInner, oracleMeasurementMiddle, oracleMeasurementOuter,
    material, Bz, kCompatibilityAbsCharge, kCompatibilityPID, oracleState, oracleChi2, policyParams, oracleReason));

  checkSurfaceKinematicStateEqual(producedCell.state(), oracleState);
  BOOST_CHECK_EQUAL(producedCell.getChi2(), oracleChi2);
}

BOOST_AUTO_TEST_CASE(DiskComputeLayerCellsRoadPreCutRejectsBeforeBuildCellSeed)
{
  Rig<MFTNLayers> rig{o2::detectors::DetID::MFT, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk};
  rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
  rig.params[0].TrackletMinPt = 0.3f;
  // Same geometry as DiskComputeLayerCellsMatchesBuildCellSeedOracle, but a
  // vanishingly small road cut: proves passesCellRoadPrecut<DiskDisk> still
  // runs before buildCellSeed<DiskDisk> and is driven by the bound
  // DiskDiskPolicyParams::cellRoadRCut, not a stale/uninitialized value.
  rig.params[0].CellRoadRCut = 1.e-6f;
  rig.establishLayout();

  const std::array<o2::its::Cluster, 3> clusters{makeGlobalCluster(1.0f, 0.5f, -0.4f, 0),
                                                 makeGlobalCluster(1.3f, 0.62f, -0.6f, 0),
                                                 makeGlobalCluster(1.7f, 0.78f, -0.9f, 0)};
  loadCandidateClusters(rig, clusters,
                        {makeDiskHit(-0.4f, 1.0f, 0.5f),
                         makeDiskHit(-0.6f, 1.3f, 0.62f),
                         makeDiskHit(-0.9f, 1.7f, 0.78f)});

  rig.traits.updateTrackingParameters(rig.params);
  rig.traits.initialiseTimeFrame(0);

  const auto topology = rig.tf.getTrackingTopologyView();
  const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
  BOOST_REQUIRE_GE(cellTopologyId, 0);

  injectCandidateTracklets(rig, cellTopologyId, clusters);

  rig.traits.computeLayerCells(0);

  BOOST_CHECK(rig.tf.getCells()[cellTopologyId].empty());
  BOOST_REQUIRE_EQUAL(rig.tf.getCellsLookupTable()[cellTopologyId].size(), 2u);
  BOOST_CHECK_EQUAL(rig.tf.getCellsLookupTable()[cellTopologyId][0], 0);
  BOOST_CHECK_EQUAL(rig.tf.getCellsLookupTable()[cellTopologyId][1], 0);
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
    Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder, nThreads};
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

    rig.traits.updateTrackingParameters(rig.params);
    rig.traits.initialiseTimeFrame(0);

    const auto topology = rig.tf.getTrackingTopologyView();
    const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
    BOOST_REQUIRE_GE(cellTopologyId, 0);

    injectCandidateTracklets(rig, cellTopologyId, clusters);

    rig.traits.computeLayerCells(0);

    Result r;
    r.cellTopologyId = cellTopologyId;
    const auto& lut = rig.tf.getCellsLookupTable()[cellTopologyId];
    r.lut.assign(lut.begin(), lut.end());
    BOOST_REQUIRE_EQUAL(rig.tf.getCells()[cellTopologyId].size(), 1u);
    const auto& cell = rig.tf.getCells()[cellTopologyId][0];
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
    Rig<MFTNLayers> rig{o2::detectors::DetID::MFT, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk, nThreads};
    rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
    rig.params[0].TrackletMinPt = 0.3f;
    rig.params[0].CellRoadRCut = 1000.f; // generous: road pre-cut must pass here
    rig.establishLayout();

    const std::array<o2::its::Cluster, 3> clusters{makeGlobalCluster(1.0f, 0.5f, -0.4f, 0),
                                                   makeGlobalCluster(1.3f, 0.62f, -0.6f, 0),
                                                   makeGlobalCluster(1.7f, 0.78f, -0.9f, 0)};
    loadCandidateClusters(rig, clusters,
                          {makeDiskHit(-0.4f, 1.0f, 0.5f),
                           makeDiskHit(-0.6f, 1.3f, 0.62f),
                           makeDiskHit(-0.9f, 1.7f, 0.78f)});

    rig.traits.updateTrackingParameters(rig.params);
    rig.traits.initialiseTimeFrame(0);

    const auto topology = rig.tf.getTrackingTopologyView();
    const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
    BOOST_REQUIRE_GE(cellTopologyId, 0);

    injectCandidateTracklets(rig, cellTopologyId, clusters);

    rig.traits.computeLayerCells(0);

    Result r;
    r.cellTopologyId = cellTopologyId;
    const auto& lut = rig.tf.getCellsLookupTable()[cellTopologyId];
    r.lut.assign(lut.begin(), lut.end());
    BOOST_REQUIRE_EQUAL(rig.tf.getCells()[cellTopologyId].size(), 1u);
    const auto& cell = rig.tf.getCells()[cellTopologyId][0];
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

// --- Gate 3 cell-road pre-cut cache: empty span safety and one-time binding

BOOST_AUTO_TEST_CASE(CylinderComputeLayerCellsSafeWithEmptyDiskReferenceSpan)
{
  // A CylinderCylinder iteration never binds a legacy MFT reference-z span:
  // TrackerTraits::mDiskLayerReferenceZ stays default-constructed (empty) for
  // the whole traversal, exactly as it does today. passesCellRoadPrecut<
  // CylinderCylinder> is still called unconditionally by the shared candidate
  // loop and must never read it. This is already incidentally exercised by
  // every Cylinder test above; this case documents and asserts it explicitly.
  Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder};
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

  rig.traits.updateTrackingParameters(rig.params);
  rig.traits.initialiseTimeFrame(0);
  BOOST_REQUIRE(rig.traits.hasTraversalCache());

  const auto topology = rig.tf.getTrackingTopologyView();
  const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
  BOOST_REQUIRE_GE(cellTopologyId, 0);

  injectCandidateTracklets(rig, cellTopologyId, clusters);

  rig.traits.computeLayerCells(0);

  BOOST_REQUIRE_EQUAL(rig.tf.getCells()[cellTopologyId].size(), 1u);
}

BOOST_AUTO_TEST_CASE(RepeatedComputeLayerCellsCallsDoNotRebindOrIncreaseCounts)
{
  // getTraversalGroupingCount()/getPolicyBindingCount() only increment inside
  // initialiseTimeFrame() (Gate 2 counters, unchanged by this slice).
  // computeLayerCells() has no code path that rebinds mDiskLayerReferenceZ or
  // any other cached policy/geometry state -- calling it repeatedly for the
  // same iteration must leave both counters exactly as they were after the
  // single initialiseTimeFrame() call, the closest observable proxy, through
  // the public API alone, for "the legacy reference-z span is bound once per
  // iteration, not once per candidate". (computeLayerCells() itself clears
  // and consumes tracklets as an existing, unrelated post-step, so a second
  // call with no freshly injected tracklets legitimately produces zero cells
  // -- that is not what this test checks.)
  Rig<MFTNLayers> rig{o2::detectors::DetID::MFT, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk};
  rig.params[0].MaxChi2ClusterAttachment = 1.e6f;
  rig.params[0].TrackletMinPt = 0.3f;
  rig.params[0].CellRoadRCut = 1000.f;
  rig.establishLayout();

  const std::array<o2::its::Cluster, 3> clusters{makeGlobalCluster(1.0f, 0.5f, -0.4f, 0),
                                                 makeGlobalCluster(1.3f, 0.62f, -0.6f, 0),
                                                 makeGlobalCluster(1.7f, 0.78f, -0.9f, 0)};
  loadCandidateClusters(rig, clusters,
                        {makeDiskHit(-0.4f, 1.0f, 0.5f),
                         makeDiskHit(-0.6f, 1.3f, 0.62f),
                         makeDiskHit(-0.9f, 1.7f, 0.78f)});

  rig.traits.updateTrackingParameters(rig.params);
  rig.traits.initialiseTimeFrame(0);
  BOOST_REQUIRE(rig.traits.hasTraversalCache());

  const int groupingCountAfterInit = rig.traits.getTraversalGroupingCount();
  const int diskBindingCountAfterInit = rig.traits.getPolicyBindingCount(TransitionPolicyTag::DiskDisk);
  const int cylinderBindingCountAfterInit = rig.traits.getPolicyBindingCount(TransitionPolicyTag::CylinderCylinder);

  const auto topology = rig.tf.getTrackingTopologyView();
  const int cellTopologyId = findCellTopologyId(topology, 0, 1, 2);
  BOOST_REQUIRE_GE(cellTopologyId, 0);

  injectCandidateTracklets(rig, cellTopologyId, clusters);

  rig.traits.computeLayerCells(0);
  BOOST_REQUIRE_EQUAL(rig.tf.getCells()[cellTopologyId].size(), 1u);
  const auto firstChi2 = rig.tf.getCells()[cellTopologyId][0].getChi2();

  rig.traits.computeLayerCells(0);
  rig.traits.computeLayerCells(0);

  BOOST_CHECK_EQUAL(rig.traits.getTraversalGroupingCount(), groupingCountAfterInit);
  BOOST_CHECK_EQUAL(rig.traits.getPolicyBindingCount(TransitionPolicyTag::DiskDisk), diskBindingCountAfterInit);
  BOOST_CHECK_EQUAL(rig.traits.getPolicyBindingCount(TransitionPolicyTag::CylinderCylinder), cylinderBindingCountAfterInit);

  // Re-inject tracklets and recompute (the underlying candidate clusters/
  // measurements loaded above are untouched -- reloading them here would
  // invalidate TrackerTraits::mLayerMeasurements without a fresh
  // initialiseTimeFrame() call to re-resolve it, which is not what this test
  // checks): a fresh call after the tracklets were consumed must still
  // reproduce the identical chi2 through the same, never-rebound
  // mDiskLayerReferenceZ/mDiskPolicyParams cache.
  injectCandidateTracklets(rig, cellTopologyId, clusters);
  rig.traits.computeLayerCells(0);

  BOOST_CHECK_EQUAL(rig.traits.getTraversalGroupingCount(), groupingCountAfterInit);
  BOOST_CHECK_EQUAL(rig.traits.getPolicyBindingCount(TransitionPolicyTag::DiskDisk), diskBindingCountAfterInit);
  BOOST_REQUIRE_EQUAL(rig.tf.getCells()[cellTopologyId].size(), 1u);
  BOOST_CHECK_EQUAL(rig.tf.getCells()[cellTopologyId][0].getChi2(), firstChi2);
}

// --- Fail-closed: no test-only seam, only the existing public entry point -

BOOST_AUTO_TEST_CASE(ComputeLayerCellsFailsClosedWithoutInitialiseTimeFrame)
{
  TimeFrame<ITSNLayers> tf;
  TrackerTraits<ITSNLayers> traits;
  traits.adoptTimeFrame(&tf);

  BOOST_CHECK_EXCEPTION(traits.computeLayerCells(0), TraversalException, [](const TraversalException& e) {
    return e.getReason() == TraversalFailureReason::InvalidTraversalSchedule;
  });
}

// --- Stage-B activation: material-correction-mode preflight wiring --------
//
// checkMaterialCorrectionModeSupport<Tag>() itself is unit-tested in
// isolation (testMaterialCorrectionModePreflight.cxx); these tests instead
// prove it is genuinely wired into TrackerTraits<NLayers>::initialiseTimeFrame()
// through the real public API, using the same Rig harness as the rest of
// this file.

BOOST_AUTO_TEST_CASE(CylinderNoneCorrectionInitialisesSuccessfully)
{
  Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder};
  rig.params[0].CorrType = o2::base::PropagatorF::MatCorrType::USEMatCorrNONE;
  rig.establishLayout();
  rig.traits.updateTrackingParameters(rig.params);
  BOOST_CHECK_NO_THROW(rig.traits.initialiseTimeFrame(0));
  BOOST_CHECK(rig.traits.hasTraversalCache());
}

BOOST_AUTO_TEST_CASE(CylinderUnsupportedMaterialCorrectionModeThrowsBeforeTimeFrameMutation)
{
  for (const auto corrType : {o2::base::PropagatorF::MatCorrType::USEMatCorrLUT, o2::base::PropagatorF::MatCorrType::USEMatCorrTGeo}) {
    Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder};
    rig.params[0].CorrType = corrType;
    rig.establishLayout();
    rig.traits.updateTrackingParameters(rig.params);

    BOOST_CHECK_EXCEPTION(rig.traits.initialiseTimeFrame(0), TraversalException, [](const TraversalException& e) {
      return e.getReason() == TraversalFailureReason::UnsupportedMaterialCorrectionMode;
    });
    // Structural failure: the traversal cache stays exactly as
    // resetTraversalCache() left it -- never partially populated -- and a
    // later call with a supported mode still succeeds, proving no lingering
    // corrupted state from the rejected attempt.
    BOOST_CHECK(!rig.traits.hasTraversalCache());
    BOOST_CHECK_EQUAL(rig.traits.getPolicyBindingCount(TransitionPolicyTag::CylinderCylinder), 0);

    rig.params[0].CorrType = o2::base::PropagatorF::MatCorrType::USEMatCorrNONE;
    rig.traits.updateTrackingParameters(rig.params);
    BOOST_CHECK_NO_THROW(rig.traits.initialiseTimeFrame(0));
    BOOST_CHECK(rig.traits.hasTraversalCache());
  }
}

BOOST_AUTO_TEST_CASE(InvalidCorrTypeRetainsExistingInvalidPolicyParametersReason)
{
  // A CorrType value isRecognizedMatCorrType() does not recognize must
  // surface as the pre-existing InvalidPolicyParameters failure (from
  // AttachHitPolicyConfigView::isValid()), never as
  // UnsupportedMaterialCorrectionMode -- the preflight explicitly defers to
  // that established classification for this case.
  Rig<ITSNLayers> rig{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder};
  rig.params[0].CorrType = static_cast<o2::base::PropagatorF::MatCorrType>(99);
  rig.establishLayout();
  rig.traits.updateTrackingParameters(rig.params);

  BOOST_CHECK_EXCEPTION(rig.traits.initialiseTimeFrame(0), TraversalException, [](const TraversalException& e) {
    return e.getReason() == TraversalFailureReason::InvalidPolicyParameters;
  });
  BOOST_CHECK(!rig.traits.hasTraversalCache());
}

BOOST_AUTO_TEST_CASE(DiskAcceptsAllRecognizedMaterialCorrectionModes)
{
  // DiskDisk's native path always uses descriptor-based nominal material
  // regardless of CorrType, so it must never be rejected by this preflight
  // for any recognized mode.
  for (const auto corrType : {o2::base::PropagatorF::MatCorrType::USEMatCorrNONE,
                              o2::base::PropagatorF::MatCorrType::USEMatCorrLUT,
                              o2::base::PropagatorF::MatCorrType::USEMatCorrTGeo}) {
    Rig<MFTNLayers> rig{o2::detectors::DetID::MFT, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk};
    rig.params[0].CorrType = corrType;
    rig.establishLayout();
    rig.traits.updateTrackingParameters(rig.params);
    BOOST_CHECK_NO_THROW(rig.traits.initialiseTimeFrame(0));
    BOOST_CHECK(rig.traits.hasTraversalCache());
  }
}
