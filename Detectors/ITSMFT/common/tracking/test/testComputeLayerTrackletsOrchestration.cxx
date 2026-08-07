// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT ComputeLayerTracklets orchestration
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <oneapi/tbb/task_arena.h>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/detail/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/detail/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/Constants.h"
#include "ITStracking/MathUtils.h"
#include "ITStracking/ROFLookupTables.h"
#include "ITStracking/Tracklet.h"
#include "MFTTracking/Constants.h"
#include "CommonConstants/MathConstants.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

constexpr float Bz = 0.5f;
constexpr std::array<unsigned char, 3> OnePixelPattern{1, 1, 0x80};

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

std::vector<SurfaceDescriptor> makeCatalog(uint16_t nLayers, o2::detectors::DetID::ID detector, SurfaceKind kind)
{
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.reserve(nLayers);
  for (uint16_t i = 0; i < nLayers; ++i) {
    surfaces.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(detector), kind});
    // Matches o2::itsmft::resetDetectorDefaults()'s per-detector LayerxX0
    // default, so TrackerTraits::initialiseTimeFrame()'s LegacyMaterialMismatch
    // compatibility check passes for these unperturbed fixtures.
    const float xOverX0 = detector == o2::detectors::DetID::MFT ? kNominalMFTLayerX0[i % MFTNLayers] : kNominalITSLayerX0[i % ITSNLayers];
    surfaces.back().material.xOverX0 = xOverX0;
    surfaces.back().material.arealDensityGPerCm2 = xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
  }
  return surfaces;
}

class PrescribedDecoder final : public ClusterDecoder
{
 public:
  PrescribedDecoder(o2::detectors::DetID::ID detector, SurfaceKind kind, std::vector<DecodedCluster> clusters)
    : mDetector{detector}, mKind{kind}, mClusters{std::move(clusters)}
  {
  }

  o2::itsmft::tracking::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const TopologyDictionary* dictionary,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool) const final
  {
    const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dictionary);
    if (!clusterData.ok()) {
      o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
      result.error = clusterData.error;
      return result;
    }

    o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
    if (externalIndex >= mClusters.size()) {
      return result;
    }
    auto decoded = mClusters[externalIndex];
    decoded.shape = clusterData.shape;
    const int layer = decoded.layer;
    result.layer = layer;
    if (layer < 0 || static_cast<size_t>(layer) >= layerToSurface.size()) {
      return result;
    }
    result.layerMapped = true;
    result.kind = mKind;
    const DetectorSensorId sensor{static_cast<uint32_t>(mDetector), decoded.sensor};
    const ClusterRef clusterRef{source, externalIndex};
    result.measurement = mKind == SurfaceKind::Disk
                           ? makeDiskSurfaceMeasurement(decoded, sensor, layerToSurface[layer], clusterRef, sourceROF)
                           : makeCylinderSurfaceMeasurement(decoded, sensor, layerToSurface[layer], clusterRef, sourceROF);
    return result;
  }

 private:
  o2::detectors::DetID::ID mDetector;
  SurfaceKind mKind;
  std::vector<DecodedCluster> mClusters;
};

struct TrackletSnapshot {
  int transitionId{-1};
  std::vector<o2::its::Tracklet> tracklets;
  std::vector<int> lookup;
  o2::its::TimeEstBC expectedTimestamp;
  bool nonparticipatingTransitionsEmpty{false};
  // Gate 4 Slice 0a additions: full per-(legacy-transitionId) tracklet/LUT
  // content and (fromLayer,toLayer) identity, for multi-transition
  // candidate-set/order/LUT parity checks that go beyond the single
  // `transitionId` above. Indices across these three vectors correspond
  // 1:1, in ascending legacy transitionId order.
  std::vector<int> allTransitionFromLayer;
  std::vector<int> allTransitionToLayer;
  std::vector<std::vector<o2::its::Tracklet>> allTracklets;
  std::vector<std::vector<int>> allLookups;
};

/// Independent acceptance oracle for the Gate 3 transition-preparation slice
/// (layerMultipleScatteringAngle<Tag>, clampTransitionCurvature<Tag>,
/// prepareTransitionScatteringAndBending, relocated into
/// TrackerTraits::initialiseTimeFrame()). Re-derives the frozen legacy
/// per-layer/per-transition formula directly -- from math_utils::MSangle
/// (barrel) or detail::mftLayerMSAngle (disk, which itself still calls the
/// legacy mftLayerZ()/LayerZCoordinate() constants internally, exactly as
/// production did before this migration) and the exact former
/// TimeFrame::initialise() transition loop -- and deliberately never calls
/// layerMultipleScatteringAngle<Tag>, clampTransitionCurvature<Tag>, or
/// prepareTransitionScatteringAndBending, so this is a genuine external
/// oracle for those operations rather than a caller of them. Preserves the
/// half-open [fromLayer, toLayer) MS accumulation range, threads oneOverR in
/// increasing legacy transitionId order exactly as the production loop
/// does, and uses the literal matching each family (`isDisk`selects `0.5f`
/// float for Disk vs `0.5` double-promoted for Cylinder, per the
/// integration review finding preserved -- not canonicalized -- in part 1/4
/// of this slice).
template <int NLayers>
void computeLegacyTransitionMSAndPhiCut(const TrackingParameters& trkParam, float bz, bool isDisk,
                                        const SurfaceGraphView& topology,
                                        gsl::span<const float> positionResolution,
                                        std::vector<float>& msAnglesOut, std::vector<float>& phiCutsOut)
{
  std::array<float, NLayers> msAngles{};
  for (unsigned int iLayer{0}; iLayer < NLayers; ++iLayer) {
    msAngles[iLayer] = isDisk ? detail::mftLayerMSAngle(iLayer, trkParam)
                              : o2::its::math_utils::MSangle(0.14f, trkParam.TrackletMinPt, trkParam.LayerxX0[iLayer]);
  }

  msAnglesOut.assign(topology.nTransitions, 0.f);
  phiCutsOut.assign(topology.nTransitions, 0.f);
  float oneOverR{0.001f * 0.3f * std::abs(bz) / trkParam.TrackletMinPt};
  for (int transitionId{0}; transitionId < static_cast<int>(topology.nTransitions); ++transitionId) {
    const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(transitionId)});
    const int from = transition.from.value();
    const int to = transition.to.value();
    float ms2 = 0.f;
    for (int layer = from; layer < to; ++layer) {
      ms2 += o2::its::math_utils::Sq(msAngles[layer]);
    }
    const float msAngle = o2::gpu::CAMath::Sqrt(ms2);
    const float r1 = trkParam.LayerRadii[from];
    const float r2 = trkParam.LayerRadii[to];
    if (isDisk) {
      oneOverR = (0.5f * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
    } else {
      oneOverR = (0.5 * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
    }
    const float res1 = o2::gpu::CAMath::Hypot(trkParam.PVres, positionResolution[from]);
    const float res2 = o2::gpu::CAMath::Hypot(trkParam.PVres, positionResolution[to]);
    const float cosTheta1half = o2::gpu::CAMath::Sqrt(1.f - o2::its::math_utils::Sq(0.5f * r1 * oneOverR));
    const float cosTheta2half = o2::gpu::CAMath::Sqrt(1.f - o2::its::math_utils::Sq(0.5f * r2 * oneOverR));
    const float x = (r2 * cosTheta1half) - (r1 * cosTheta2half);
    const float delta = o2::gpu::CAMath::Sqrt(1.f / (1.f - 0.25f * o2::its::math_utils::Sq(x * oneOverR)) *
                                              (o2::its::math_utils::Sq((0.25f * r1 * r2 * o2::its::math_utils::Sq(oneOverR) / cosTheta2half) + cosTheta1half) * o2::its::math_utils::Sq(res1) +
                                               o2::its::math_utils::Sq((0.25f * r1 * r2 * o2::its::math_utils::Sq(oneOverR) / cosTheta1half) + cosTheta2half) * o2::its::math_utils::Sq(res2)));
    msAnglesOut[transitionId] = msAngle;
    phiCutsOut[transitionId] = o2::gpu::CAMath::Min(o2::gpu::CAMath::ASin(0.5f * x * oneOverR) + 2.f * msAngle + delta, o2::constants::math::PI * 0.5f);
  }
}

template <int NLayers>
TrackletSnapshot runFixture(o2::detectors::DetID::ID detector,
                            SurfaceKind kind,
                            SurfaceKind tag,
                            std::vector<DecodedCluster> decoded,
                            int nThreads,
                            std::function<void(TrackingParameters&)> customizeParams = {})
{
  auto pool = std::make_shared<BoundedMemoryResource>();
  // Gate 4 B3.1: `frame` declared before `tf` so it is constructed first and
  // destroyed last (see SurfaceTrackingScratch's own lifetime-contract doc).
  TimeFrame frame;
  SurfaceTrackingScratch tf;
  TrackerTraits traits;
  std::shared_ptr<tbb::task_arena> arena;
  std::vector<TrackingParameters> params(1);
  resetDetectorDefaults(params[0], detector);
  params[0].UseDiamond = true;
  params[0].CreateArtefactLabels = false;
  params[0].PassFlags.reset();
  params[0].PassFlags.set(IterationStep::FirstPass, IterationStep::RebuildClusterLUT);
  if (customizeParams) {
    customizeParams(params[0]);
  }

  frame.setMemoryPool(pool);
  tf.setMemoryPool(pool);
  traits.setMemoryPool(pool);
  traits.setNThreads(nThreads, arena);
  traits.adoptScratch(&tf);
  traits.adoptFrame(&frame);
  traits.updateTrackingParameters(params);
  traits.setBz(Bz);

  const auto orderedSurfaces = identitySurfaces(static_cast<uint16_t>(NLayers));
  const auto catalog = makeCatalog(static_cast<uint16_t>(NLayers), detector, kind);
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  auto planResult = buildSurfaceGraphs(catalogView, orderedSurfaces, params);
  BOOST_REQUIRE(planResult.ok());
  const auto graph = std::move(planResult.graphs.front());
  const std::vector<SurfaceGraph> plan{graph};
  const auto layoutView = graph.getView();
  tf.adoptPlan(orderedSurfaces.size(), layoutView.nTransitions, layoutView.nCells);

  std::vector<CompClusterExt> compactClusters;
  std::vector<unsigned char> patterns;
  compactClusters.reserve(decoded.size());
  patterns.reserve(decoded.size() * OnePixelPattern.size());
  for (const auto& cluster : decoded) {
    compactClusters.emplace_back(0, 0, CompCluster::InvalidPatternID, cluster.layer);
    patterns.insert(patterns.end(), OnePixelPattern.begin(), OnePixelPattern.end());
  }
  const std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, static_cast<int>(compactClusters.size())}};
  PrescribedDecoder decoder{detector, kind, std::move(decoded)};
  const auto load = tf.loadNormalizedSource(frame, decoder, o2::InteractionRecord{50, 5}, ROFTimingConfig{40, 0, 0, 0},
                                            compactClusters, patterns, rofs, &dict(), nullptr, detector,
                                            gsl::span<const SurfaceId>{graph.getOrderedSurfaces()}, graph.getSurfaceCatalog());
  BOOST_REQUIRE(load.ok());

  o2::its::LayerTiming layerTiming{};
  layerTiming.mNROFsTF = 1;
  layerTiming.mROFLength = 40;
  o2::its::ROFOverlapTable<NLayers> rofTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    rofTable.defineLayer(layer, layerTiming);
  }
  rofTable.init();
  // Real production workflow timing construction
  // always builds and sets this alongside the ROFOverlapTable above, from
  // the same per-layer LayerTiming, regardless of UseDiamond -- the diamond
  // vertex derived per-ROF for tracklet finding (TrackerTraits.cxx) is
  // checked through the genuine isVertexCompatible() on this table, not a
  // useDiamond-skipped shortcut, so this fixture needs it populated too.
  o2::its::ROFVertexLookupTable<NLayers> vtxTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    vtxTable.defineLayer(layer, layerTiming);
  }
  vtxTable.init();
  o2::its::ROFMaskTable<NLayers> mask{rofTable};
  mask.resetMask();
  for (int layer = 0; layer < NLayers; ++layer) {
    mask.setROFsEnabled(layer, 0, 1, 1);
  }
  tf.setROFViews(RuntimeROFViews{rofTable.getView(), vtxTable.getView(), mask.getView(), {}});

  SurfaceMask owned;
  for (const auto surface : graph.getOrderedSurfaces()) {
    owned.set(surface);
  }
  auto bindingResult = SurfacePlanBinding::build(graph.getView(), ClusterSourceId{0}, owned,
                                                 graph.getOrderedSurfaces(), kind);
  BOOST_REQUIRE(bindingResult.ok());
  traits.adoptSurfacePlanBinding(bindingResult.binding.get());
  traits.initialiseTimeFrame(0, plan);
  BOOST_CHECK_EQUAL(traits.getTraversalGroupingCount(), 1);

  // Gate 3 transition-preparation slice: successful initialisation must fill
  // every transition entry (relocated from TimeFrame::initialise() into
  // TrackerTraits::initialiseTimeFrame(), see TrackletFinding.h).
  // Exercised here for both Cylinder and Disk through the
  // existing fixture rather than a separate harness. Beyond finiteness, each
  // entry is checked bit-for-bit against computeLegacyTransitionMSAndPhiCut's
  // independent oracle -- the only replay-grade acceptance evidence for the
  // common Cylinder path, since no real-geometry common-CA ITS
  // replay exists yet.
  {
    const auto preparedTopology = layoutView;
    const auto& msAngles = tf.getTransitionMSAngles();
    const auto& phiCuts = tf.getTransitionPhiCuts();
    BOOST_REQUIRE_EQUAL(msAngles.size(), static_cast<size_t>(preparedTopology.nTransitions));
    BOOST_REQUIRE_EQUAL(phiCuts.size(), static_cast<size_t>(preparedTopology.nTransitions));
    for (int id = 0; id < preparedTopology.nTransitions; ++id) {
      BOOST_CHECK(std::isfinite(msAngles[id]));
      BOOST_CHECK(std::isfinite(phiCuts[id]));
    }

    std::vector<float> positionResolution(NLayers);
    for (int layer = 0; layer < NLayers; ++layer) {
      positionResolution[layer] = tf.getPositionResolution(layer);
    }
    std::vector<float> expectedMSAngles;
    std::vector<float> expectedPhiCuts;
    computeLegacyTransitionMSAndPhiCut<NLayers>(params[0], Bz, kind == SurfaceKind::Disk, preparedTopology,
                                                gsl::span<const float>(positionResolution.data(), positionResolution.size()),
                                                expectedMSAngles, expectedPhiCuts);
    BOOST_REQUIRE_EQUAL(expectedMSAngles.size(), msAngles.size());
    BOOST_REQUIRE_EQUAL(expectedPhiCuts.size(), phiCuts.size());
    for (int id = 0; id < preparedTopology.nTransitions; ++id) {
      BOOST_CHECK_EQUAL(msAngles[id], expectedMSAngles[id]);
      BOOST_CHECK_EQUAL(phiCuts[id], expectedPhiCuts[id]);
    }
  }

  const auto topology = layoutView;
  int transitionId = -1;
  for (int id = 0; id < topology.nTransitions; ++id) {
    const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(id)});
    if (transition.from.value() == 0 && transition.to.value() == 1) {
      transitionId = id;
      break;
    }
  }
  BOOST_REQUIRE_GE(transitionId, 0);

  // Repeated per-vertex processing must reuse the initialization-time
  // traversal grouping and typed parameter binding.
  traits.computeLayerTracklets(0, 0);
  BOOST_CHECK_EQUAL(traits.getTraversalGroupingCount(), 1);
  traits.computeLayerTracklets(0, 0);
  BOOST_CHECK_EQUAL(traits.getTraversalGroupingCount(), 1);

  TrackletSnapshot result;
  result.transitionId = transitionId;
  result.expectedTimestamp = tf.getROFOverlapView().getTimeStamp(0, 0, 1, 0);
  const auto& tracklets = tf.getTracklets()[transitionId];
  result.tracklets.assign(tracklets.begin(), tracklets.end());
  const auto& lookup = tf.getTrackletsLookupTable()[transitionId];
  result.lookup.assign(lookup.begin(), lookup.end());
  result.nonparticipatingTransitionsEmpty = true;
  for (int id = 0; id < topology.nTransitions; ++id) {
    if (id != transitionId && !tf.getTracklets()[id].empty()) {
      result.nonparticipatingTransitionsEmpty = false;
      break;
    }
  }

  // Gate 4 Slice 0a: full per-transition snapshot, ascending legacy
  // transitionId order, for multi-transition candidate-set/order/LUT parity
  // checks (see e.g. ItsIdentityLayoutTrackletsSpanMultipleAdjacentTransitionsInOrder).
  for (int id = 0; id < topology.nTransitions; ++id) {
    const auto& transition = topology.getTransition(TransitionId{static_cast<uint16_t>(id)});
    result.allTransitionFromLayer.push_back(transition.from.value());
    result.allTransitionToLayer.push_back(transition.to.value());
    const auto& idTracklets = tf.getTracklets()[id];
    result.allTracklets.emplace_back(idTracklets.begin(), idTracklets.end());
    const auto& idLookup = tf.getTrackletsLookupTable()[id];
    result.allLookups.emplace_back(idLookup.begin(), idLookup.end());
  }
  return result;
}

void checkSame(const TrackletSnapshot& serial, const TrackletSnapshot& parallel)
{
  BOOST_CHECK_EQUAL(serial.transitionId, parallel.transitionId);
  BOOST_REQUIRE_EQUAL(serial.tracklets.size(), parallel.tracklets.size());
  BOOST_CHECK_EQUAL_COLLECTIONS(serial.lookup.begin(), serial.lookup.end(), parallel.lookup.begin(), parallel.lookup.end());
  for (size_t i = 0; i < serial.tracklets.size(); ++i) {
    BOOST_CHECK(serial.tracklets[i] == parallel.tracklets[i]);
    BOOST_CHECK_EQUAL(serial.tracklets[i].tanLambda, parallel.tracklets[i].tanLambda);
    BOOST_CHECK_EQUAL(serial.tracklets[i].phi, parallel.tracklets[i].phi);
    BOOST_CHECK_EQUAL(serial.tracklets[i].getTimeStamp().getTimeStamp(), parallel.tracklets[i].getTimeStamp().getTimeStamp());
    BOOST_CHECK_EQUAL(serial.tracklets[i].getTimeStamp().getTimeStampError(), parallel.tracklets[i].getTimeStamp().getTimeStampError());
  }
}

void checkExactTracklet(const TrackletSnapshot& snapshot, float expectedTanLambda, float expectedPhi)
{
  BOOST_REQUIRE_EQUAL(snapshot.tracklets.size(), 1u);
  const auto& tracklet = snapshot.tracklets.front();
  BOOST_CHECK_EQUAL(tracklet.firstClusterIndex, 0);
  BOOST_CHECK_EQUAL(tracklet.secondClusterIndex, 0);
  BOOST_CHECK_EQUAL(tracklet.tanLambda, expectedTanLambda);
  BOOST_CHECK_EQUAL(tracklet.phi, expectedPhi);
  BOOST_CHECK_EQUAL(tracklet.getTimeStamp().getTimeStamp(), snapshot.expectedTimestamp.getTimeStamp());
  BOOST_CHECK_EQUAL(tracklet.getTimeStamp().getTimeStampError(), snapshot.expectedTimestamp.getTimeStampError());
  const std::vector<int> expectedLookup{0, 1};
  BOOST_CHECK_EQUAL_COLLECTIONS(snapshot.lookup.begin(), snapshot.lookup.end(), expectedLookup.begin(), expectedLookup.end());
  BOOST_CHECK(snapshot.nonparticipatingTransitionsEmpty);
}

DecodedCluster cylinderCluster(float radius, float z, int layer)
{
  DecodedCluster cluster{};
  cluster.global = {radius, 0.f, z};
  cluster.cylinderFrame = {radius, 0.f, z, 0.f};
  cluster.rowColumnCovariance = {1.e-4f, 0.f, 1.e-4f};
  cluster.sensor = static_cast<uint32_t>(layer);
  cluster.layer = layer;
  return cluster;
}

DecodedCluster diskCluster(float x, float y, float z, int layer)
{
  DecodedCluster cluster{};
  cluster.global = {x, y, z};
  cluster.rowColumnCovariance = {1.e-2f, 0.f, 1.e-2f};
  cluster.sensor = static_cast<uint32_t>(layer);
  cluster.layer = layer;
  return cluster;
}

/// A chain of `nHops + 1` disk clusters (layers 0..nHops) consistent with a
/// single forward trajectory: each hop's target position is computed by
/// projecting from the previous hop's own cluster position via
/// detail::mftTrackletProject -- the same primitive
/// projectDiskSearchWindow itself uses internally -- so every adjacent
/// pair in the chain is a genuine geometric match, not just the first one.
std::vector<DecodedCluster> buildMftChainClusters(const TrackingParameters& params, float bz, int nHops)
{
  std::vector<DecodedCluster> clusters;
  float x = 1.f, y = 0.5f;
  float z = detail::mftLayerZ(0);
  clusters.push_back(diskCluster(x, y, z, 0));
  for (int hop = 0; hop < nHops; ++hop) {
    const float nextZ = detail::mftLayerZ(hop + 1);
    float targetX = 0.f, targetY = 0.f;
    detail::mftTrackletProject(x, y, z, params.Diamond[0], params.Diamond[1], params.Diamond[2],
                               hop, hop + 1, bz, params.TrackletMinPt, targetX, targetY);
    clusters.push_back(diskCluster(targetX, targetY, nextZ, hop + 1));
    x = targetX;
    y = targetY;
    z = nextZ;
  }
  return clusters;
}

/// Gate 4 Slice 0a fail-closed coverage, revised for Gate 4 B2 Slice 2: a
/// production plan always builds a single-subgraph, single-kind layout
/// from a duplicate-free orderedSurfaces span (SurfaceGraphBuilder already
/// rejects a duplicate SurfaceId within one subgraph, and
/// buildSurfaceGraphs() has no way to author a combined/mixed-kind
/// layout at all). To directly exercise TrackerTraits::initialiseTimeFrame()'s
/// own fail-closed checks (SurfaceLayerMappingMismatch, MixedSurfaceKindLayout)
/// against inputs that cannot arise through that production path, the two
/// tests below construct a deliberately-corrupted std::vector<SurfaceGraph> directly
/// and pass it to initialiseTimeFrame() as its explicit plan parameter -- no
/// TimeFrame-subclass injection hack is needed any more, since the plan is no
/// longer TimeFrame-owned state to smuggle in.

/// A valid, identity-ordered chain SurfaceGraph for `detector`/`kind`
/// (same shape buildSurfaceGraphs() would build for that detector), built
/// directly via SurfaceGraphBuilder so it can be installed with a
/// deliberately-corrupted SurfaceGraphConfigurationKey::orderedSurfaces
/// afterwards.
std::pair<SurfaceGraph, std::vector<SurfaceDescriptor>> buildIdentityChainLayout(uint16_t nSurfaces, o2::detectors::DetID::ID detector, SurfaceKind kind)
{
  auto surfaces = makeCatalog(nSurfaces, detector, kind);
  SurfaceGraphBuilder builder{SurfaceCatalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())}};
  SurfaceGraphSubgraph subgraph;
  subgraph.orderedSurfaces = identitySurfaces(nSurfaces);
  subgraph.maxHoles = 0;
  auto result = builder.addSubgraph(std::move(subgraph)).build();
  BOOST_REQUIRE(result.ok());
  return {std::move(*result.graph), std::move(surfaces)};
}

std::vector<SurfaceId> rangeSurfaces(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> ids;
  ids.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    ids.push_back(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return ids;
}

/// Disconnected catalog spanning [0, nCylinders) as Cylinder/ITS surfaces and
/// [nCylinders, nCylinders + nDisks) as Disk/MFT surfaces, in one shared
/// global SurfaceId space -- material is left default-initialized since
/// CombinedCylinderAndDiskLayoutIsRejectedBeforeTrackletProcessing never
/// reaches the material-validation step.
std::vector<SurfaceDescriptor> combinedCatalog(uint16_t nCylinders, uint16_t nDisks)
{
  std::vector<SurfaceDescriptor> surfaces;
  for (uint16_t s = 0; s < nCylinders; ++s) {
    surfaces.push_back(SurfaceDescriptor{SurfaceId{s}, s, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  }
  for (uint16_t s = 0; s < nDisks; ++s) {
    const auto id = SurfaceId{static_cast<uint16_t>(nCylinders + s)};
    surfaces.push_back(SurfaceDescriptor{id, s, static_cast<uint8_t>(o2::detectors::DetID::MFT), SurfaceKind::Disk});
  }
  return surfaces;
}

} // namespace

BOOST_AUTO_TEST_CASE(CylinderOnePassAndTwoPassProduceIdenticalTracklets)
{
  const std::vector<DecodedCluster> clusters{
    cylinderCluster(3.f, 0.3f, 0),
    cylinderCluster(4.f, 0.4f, 1)};
  const auto serial = runFixture<ITSNLayers>(o2::detectors::DetID::ITS, SurfaceKind::Cylinder,
                                             SurfaceKind::Cylinder, clusters, 1);
  const auto parallel = runFixture<ITSNLayers>(o2::detectors::DetID::ITS, SurfaceKind::Cylinder,
                                               SurfaceKind::Cylinder, clusters, 4);
  checkExactTracklet(serial, (0.3f - 0.4f) / (3.f - 4.f), o2::gpu::CAMath::ATan2(0.f, -1.f));
  checkExactTracklet(parallel, (0.3f - 0.4f) / (3.f - 4.f), o2::gpu::CAMath::ATan2(0.f, -1.f));
  checkSame(serial, parallel);
}

BOOST_AUTO_TEST_CASE(DiskOnePassAndTwoPassProduceIdenticalTracklets)
{
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::MFT);
  const float fromZ = detail::mftLayerZ(0);
  const float toZ = detail::mftLayerZ(1);
  float targetX = 0.f;
  float targetY = 0.f;
  detail::mftTrackletProject(1.f, 0.5f, fromZ,
                             params.Diamond[0], params.Diamond[1], params.Diamond[2],
                             0, 1, Bz, params.TrackletMinPt, targetX, targetY);
  const std::vector<DecodedCluster> clusters{
    diskCluster(1.f, 0.5f, fromZ, 0),
    diskCluster(targetX, targetY, toZ, 1)};
  const auto serial = runFixture<MFTNLayers>(o2::detectors::DetID::MFT, SurfaceKind::Disk,
                                             SurfaceKind::Disk, clusters, 1);
  const auto parallel = runFixture<MFTNLayers>(o2::detectors::DetID::MFT, SurfaceKind::Disk,
                                               SurfaceKind::Disk, clusters, 4);
  const float expectedTanLambda = (fromZ - toZ) / (toZ - fromZ);
  const float expectedPhi = o2::gpu::CAMath::ATan2(0.5f - targetY, 1.f - targetX);
  checkExactTracklet(serial, expectedTanLambda, expectedPhi);
  checkExactTracklet(parallel, expectedTanLambda, expectedPhi);
  checkSame(serial, parallel);
}

BOOST_AUTO_TEST_CASE(ComputeLayerTrackletsFailsClosedWithoutInitialiseTimeFrame)
{
  SurfaceTrackingScratch tf;
  TrackerTraits traits;
  traits.adoptScratch(&tf);

  BOOST_CHECK_EXCEPTION(traits.computeLayerTracklets(0, -1), TraversalException, [](const TraversalException& error) {
    return error.getReason() == TraversalFailureReason::InvalidTraversalSchedule;
  });
}

BOOST_AUTO_TEST_CASE(InitialiseTimeFrameFailureLeavesTransitionArraysZeroFilledNotPartial)
{
  // Gate 3 transition-preparation slice failure contract: TimeFrame::initialise()
  // already clears/resizes mTransitionMSAngles/mTransitionPhiCuts to
  // nTransitions before any surface-kind/geometry validation runs (unchanged by
  // this slice). A later fallible check (here: an invalid CorrType, so
  // AttachHitConfigView::isValid() fails) must leave those arrays
  // exactly zero-filled at the correct size -- never a mixture of computed
  // and zero entries -- because the (non-throwing) value-computation loop
  // this slice added never starts until every fallible check has succeeded.
  // Deliberately not a corrupted LayerxX0 value: since the compatibility-only
  // LegacyMaterialMismatch check (TrackerTraits::initialiseTimeFrame()) now
  // runs before TimeFrame::initialise(), a corrupted LayerxX0 would fail
  // there instead, before the arrays are ever resized -- this test wants a
  // fallible check that fires strictly after TimeFrame::initialise().
  auto pool = std::make_shared<BoundedMemoryResource>();
  // Gate 4 B3.1: `frame` declared before `tf` so it is constructed first and
  // destroyed last (see SurfaceTrackingScratch's own lifetime-contract doc).
  TimeFrame frame;
  SurfaceTrackingScratch tf;
  TrackerTraits traits;
  std::shared_ptr<tbb::task_arena> arena;
  std::vector<TrackingParameters> params(1);
  resetDetectorDefaults(params[0], o2::detectors::DetID::ITS);
  params[0].PassFlags.reset();
  params[0].PassFlags.set(IterationStep::FirstPass, IterationStep::RebuildClusterLUT);
  params[0].CorrType = static_cast<o2::base::PropagatorImpl<float>::MatCorrType>(99); // invalid: AttachHitConfigView rejects it

  frame.setMemoryPool(pool);
  tf.setMemoryPool(pool);
  traits.setMemoryPool(pool);
  traits.setNThreads(1, arena);
  traits.adoptScratch(&tf);
  traits.adoptFrame(&frame);
  traits.updateTrackingParameters(params);
  traits.setBz(Bz);

  const auto orderedSurfaces = identitySurfaces(static_cast<uint16_t>(ITSNLayers));
  const auto catalog = makeCatalog(static_cast<uint16_t>(ITSNLayers), o2::detectors::DetID::ITS, SurfaceKind::Cylinder);
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  auto planResult = buildSurfaceGraphs(catalogView, orderedSurfaces, params);
  BOOST_REQUIRE(planResult.ok());
  const auto graph = std::move(planResult.graphs.front());
  const std::vector<SurfaceGraph> plan{graph};
  const auto layoutView = graph.getView();
  tf.adoptPlan(orderedSurfaces.size(), layoutView.nTransitions, layoutView.nCells);

  // Same minimal cluster/ROF/mask setup as runFixture(): TimeFrame::initialise()
  // (called unconditionally, before any of this test's induced failure) needs
  // it to size mIndexTables/mClusters correctly, regardless of what this test
  // is actually probing.
  const std::vector<DecodedCluster> decoded{cylinderCluster(3.f, 0.3f, 0), cylinderCluster(4.f, 0.4f, 1)};
  std::vector<CompClusterExt> compactClusters;
  std::vector<unsigned char> patterns;
  compactClusters.reserve(decoded.size());
  patterns.reserve(decoded.size() * OnePixelPattern.size());
  for (const auto& cluster : decoded) {
    compactClusters.emplace_back(0, 0, CompCluster::InvalidPatternID, cluster.layer);
    patterns.insert(patterns.end(), OnePixelPattern.begin(), OnePixelPattern.end());
  }
  const std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, static_cast<int>(compactClusters.size())}};
  PrescribedDecoder decoder{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, decoded};
  const auto load = tf.loadNormalizedSource(frame, decoder, o2::InteractionRecord{50, 5}, ROFTimingConfig{40, 0, 0, 0},
                                            compactClusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS,
                                            gsl::span<const SurfaceId>{graph.getOrderedSurfaces()}, graph.getSurfaceCatalog());
  BOOST_REQUIRE(load.ok());

  o2::its::LayerTiming layerTiming{};
  layerTiming.mNROFsTF = 1;
  layerTiming.mROFLength = 40;
  o2::its::ROFOverlapTable<ITSNLayers> rofTable;
  for (int layer = 0; layer < ITSNLayers; ++layer) {
    rofTable.defineLayer(layer, layerTiming);
  }
  rofTable.init();
  o2::its::ROFVertexLookupTable<ITSNLayers> vtxTable;
  for (int layer = 0; layer < ITSNLayers; ++layer) {
    vtxTable.defineLayer(layer, layerTiming);
  }
  vtxTable.init();
  o2::its::ROFMaskTable<ITSNLayers> mask{rofTable};
  mask.resetMask();
  for (int layer = 0; layer < ITSNLayers; ++layer) {
    mask.setROFsEnabled(layer, 0, 1, 1);
  }
  tf.setROFViews(RuntimeROFViews{rofTable.getView(), vtxTable.getView(), mask.getView(), {}});

  SurfaceMask owned;
  for (const auto surface : graph.getOrderedSurfaces()) {
    owned.set(surface);
  }
  auto bindingResult = SurfacePlanBinding::build(graph.getView(), ClusterSourceId{0}, owned,
                                                 graph.getOrderedSurfaces(), SurfaceKind::Cylinder);
  BOOST_REQUIRE(bindingResult.ok());
  traits.adoptSurfacePlanBinding(bindingResult.binding.get());

  BOOST_CHECK_EXCEPTION(traits.initialiseTimeFrame(0, plan), TraversalException, [](const TraversalException& error) {
    return error.getReason() == TraversalFailureReason::InvalidSurfaceParameters;
  });

  const auto topology = layoutView;
  const auto& msAngles = tf.getTransitionMSAngles();
  const auto& phiCuts = tf.getTransitionPhiCuts();
  BOOST_REQUIRE_EQUAL(msAngles.size(), static_cast<size_t>(topology.nTransitions));
  BOOST_REQUIRE_EQUAL(phiCuts.size(), static_cast<size_t>(topology.nTransitions));
  for (int id = 0; id < topology.nTransitions; ++id) {
    BOOST_CHECK_EQUAL(msAngles[id], 0.f);
    BOOST_CHECK_EQUAL(phiCuts[id], 0.f);
  }
}

// ---------------------------------------------------------------------------
// Gate 4 Slice 0a (sparse-topology tracklet migration) additions below.
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ItsIdentityLayoutTrackletsSpanMultipleAdjacentTransitionsInOrder)
{
  // Collinear track across 4 barrel layers (z = 0.1 * r for every cluster).
  // Under ITS's default MaxHoles=0 only strictly-adjacent transitions exist
  // at all, so this directly proves transition-level tracklet/LUT/order
  // parity across three distinct transitions simultaneously -- each
  // resolved through the migrated computeLayerTrackletsForKind() via a
  // fresh mSurfaceToLegacyLayer lookup -- not just the single transition the
  // tests above check, while every non-participating transition (touching
  // layers 4/5/6) stays empty.
  const std::vector<DecodedCluster> clusters{
    cylinderCluster(3.f, 0.3f, 0),
    cylinderCluster(4.f, 0.4f, 1),
    cylinderCluster(5.f, 0.5f, 2),
    cylinderCluster(6.f, 0.6f, 3)};
  const auto snapshot = runFixture<ITSNLayers>(o2::detectors::DetID::ITS, SurfaceKind::Cylinder,
                                               SurfaceKind::Cylinder, clusters, 1);
  // Each transition's expected tanLambda is computed from its own specific
  // (radius, z) pair rather than one shared constant: although every pair
  // shares the same nominal slope (z = 0.1 * r), float subtraction/division
  // of different operand pairs does not generally round to the identical
  // bit pattern even when the mathematical result is the same value.
  constexpr std::array<float, 4> radii{3.f, 4.f, 5.f, 6.f};
  constexpr std::array<float, 4> zs{0.3f, 0.4f, 0.5f, 0.6f};
  const float expectedPhi = o2::gpu::CAMath::ATan2(0.f, -1.f);
  const std::vector<int> expectedLookup{0, 1};

  BOOST_REQUIRE_EQUAL(snapshot.allTransitionFromLayer.size(), snapshot.allTracklets.size());
  BOOST_REQUIRE_EQUAL(snapshot.allTransitionFromLayer.size(), snapshot.allLookups.size());
  bool sawTransition01 = false, sawTransition12 = false, sawTransition23 = false;
  for (size_t id = 0; id < snapshot.allTransitionFromLayer.size(); ++id) {
    const int from = snapshot.allTransitionFromLayer[id];
    const int to = snapshot.allTransitionToLayer[id];
    const bool participates = (from == 0 && to == 1) || (from == 1 && to == 2) || (from == 2 && to == 3);
    if (participates) {
      BOOST_REQUIRE_EQUAL(snapshot.allTracklets[id].size(), 1u);
      const auto& tracklet = snapshot.allTracklets[id].front();
      BOOST_CHECK_EQUAL(tracklet.firstClusterIndex, 0);
      BOOST_CHECK_EQUAL(tracklet.secondClusterIndex, 0);
      const float expectedTanLambda = (zs[from] - zs[to]) / (radii[from] - radii[to]);
      BOOST_CHECK_EQUAL(tracklet.tanLambda, expectedTanLambda);
      BOOST_CHECK_EQUAL(tracklet.phi, expectedPhi);
      BOOST_CHECK_EQUAL_COLLECTIONS(snapshot.allLookups[id].begin(), snapshot.allLookups[id].end(), expectedLookup.begin(), expectedLookup.end());
      sawTransition01 |= (from == 0 && to == 1);
      sawTransition12 |= (from == 1 && to == 2);
      sawTransition23 |= (from == 2 && to == 3);
    } else {
      BOOST_CHECK(snapshot.allTracklets[id].empty());
    }
  }
  BOOST_CHECK(sawTransition01);
  BOOST_CHECK(sawTransition12);
  BOOST_CHECK(sawTransition23);
}

BOOST_AUTO_TEST_CASE(MftIdentityLayoutTrackletsSpanMultipleAdjacentTransitionsInOrder)
{
  // Same multi-transition parity property for the Disk/forward family:
  // a 4-disk chain built hop-by-hop with detail::mftTrackletProject (the
  // same primitive projectDiskSearchWindow uses internally), proving
  // transitions (0,1),(1,2),(2,3) each get exactly one correctly-ordered
  // tracklet and every other transition stays empty.
  TrackingParameters params;
  resetDetectorDefaults(params, o2::detectors::DetID::MFT);
  const auto clusters = buildMftChainClusters(params, Bz, 3);
  BOOST_REQUIRE_EQUAL(clusters.size(), 4u);
  const auto snapshot = runFixture<MFTNLayers>(o2::detectors::DetID::MFT, SurfaceKind::Disk,
                                               SurfaceKind::Disk, clusters, 1);
  const std::vector<int> expectedLookup{0, 1};

  BOOST_REQUIRE_EQUAL(snapshot.allTransitionFromLayer.size(), snapshot.allTracklets.size());
  BOOST_REQUIRE_EQUAL(snapshot.allTransitionFromLayer.size(), snapshot.allLookups.size());
  bool sawTransition01 = false, sawTransition12 = false, sawTransition23 = false;
  for (size_t id = 0; id < snapshot.allTransitionFromLayer.size(); ++id) {
    const int from = snapshot.allTransitionFromLayer[id];
    const int to = snapshot.allTransitionToLayer[id];
    const bool participates = (from == 0 && to == 1) || (from == 1 && to == 2) || (from == 2 && to == 3);
    if (participates) {
      BOOST_REQUIRE_EQUAL(snapshot.allTracklets[id].size(), 1u);
      const auto& tracklet = snapshot.allTracklets[id].front();
      BOOST_CHECK_EQUAL(tracklet.firstClusterIndex, 0);
      BOOST_CHECK_EQUAL(tracklet.secondClusterIndex, 0);
      // (fromZ - toZ) / (toZ - fromZ) is exactly -1 for any two distinct
      // z values -- the same identity DiskOnePassAndTwoPassProduceIdenticalTracklets
      // above already relies on for its single-hop expectedTanLambda.
      BOOST_CHECK_EQUAL(tracklet.tanLambda, -1.f);
      BOOST_CHECK_EQUAL_COLLECTIONS(snapshot.allLookups[id].begin(), snapshot.allLookups[id].end(), expectedLookup.begin(), expectedLookup.end());
      sawTransition01 |= (from == 0 && to == 1);
      sawTransition12 |= (from == 1 && to == 2);
      sawTransition23 |= (from == 2 && to == 3);
    } else {
      BOOST_CHECK(snapshot.allTracklets[id].empty());
    }
  }
  BOOST_CHECK(sawTransition01);
  BOOST_CHECK(sawTransition12);
  BOOST_CHECK(sawTransition23);
}

BOOST_AUTO_TEST_CASE(ItsHoleTransitionTrackletResolvesCorrectLegacyLayerEndpoints)
{
  // MaxHoles=1 with layer 1 an allowed hole introduces a (0,2)-skip-1
  // transition whose sparse SurfaceTransition endpoints are SurfaceId{0}/
  // SurfaceId{2} -- a direct, non-adjacent exercise of mSurfaceToLegacyLayer
  // resolving a transition's endpoints correctly, and of hole/skipped-surface
  // behaviour staying identical to the pre-migration code (which read the
  // same fromLayer/toLayer straight off the legacy view). No cluster is
  // placed on layer 1 at all, so only the hole transition can produce a
  // tracklet.
  const std::vector<DecodedCluster> clusters{
    cylinderCluster(3.f, 0.3f, 0),
    cylinderCluster(5.f, 0.5f, 2)};
  const auto snapshot = runFixture<ITSNLayers>(
    o2::detectors::DetID::ITS, SurfaceKind::Cylinder, SurfaceKind::Cylinder, clusters, 1,
    [](TrackingParameters& p) {
      p.MaxHoles = 1;
      p.HoleLayerMask = LayerMask{static_cast<uint16_t>(1u << 1)};
    });

  const float expectedTanLambda = (0.3f - 0.5f) / (3.f - 5.f);
  const float expectedPhi = o2::gpu::CAMath::ATan2(0.f, -1.f);
  bool sawHoleTransition = false;
  BOOST_REQUIRE_EQUAL(snapshot.allTransitionFromLayer.size(), snapshot.allTracklets.size());
  for (size_t id = 0; id < snapshot.allTransitionFromLayer.size(); ++id) {
    const int from = snapshot.allTransitionFromLayer[id];
    const int to = snapshot.allTransitionToLayer[id];
    if (from == 0 && to == 2) {
      sawHoleTransition = true;
      BOOST_REQUIRE_EQUAL(snapshot.allTracklets[id].size(), 1u);
      const auto& tracklet = snapshot.allTracklets[id].front();
      BOOST_CHECK_EQUAL(tracklet.tanLambda, expectedTanLambda);
      BOOST_CHECK_EQUAL(tracklet.phi, expectedPhi);
      const std::vector<int> expectedLookup{0, 1};
      BOOST_CHECK_EQUAL_COLLECTIONS(snapshot.allLookups[id].begin(), snapshot.allLookups[id].end(), expectedLookup.begin(), expectedLookup.end());
    } else {
      BOOST_CHECK(snapshot.allTracklets[id].empty());
    }
  }
  BOOST_CHECK(sawHoleTransition);
}

BOOST_AUTO_TEST_CASE(DuplicateSurfaceIdMappingFailsClosedBeforeTrackletProcessing)
{
  // buildSurfaceGraphs() -- the production path -- always builds its
  // subgraph from the caller-supplied orderedSurfaces, and SurfaceGraphBuilder
  // already rejects a duplicate SurfaceId within that subgraph
  // (DuplicateSurfaceInSubgraph), so a duplicate mapping can never reach
  // TrackerTraits::initialiseTimeFrame() through that path. This test
  // constructs an otherwise-valid, identity-topology std::vector<SurfaceGraph>
  // directly, with a SurfaceGraphConfigurationKey::orderedSurfaces that
  // duplicates SurfaceId{0} at legacy layers 0 and 1, then passes it to
  // initialiseTimeFrame() as its explicit plan argument -- exercising
  // exactly (and only) the mSurfaceToLegacyLayer bijectivity preflight,
  // through that method's own public contract.
  auto built = buildIdentityChainLayout(static_cast<uint16_t>(ITSNLayers), o2::detectors::DetID::ITS, SurfaceKind::Cylinder);
  const std::vector<SurfaceDescriptor> surfaces = std::move(built.second);

  auto planGraph = std::move(built.first);
  planGraph.setOrderedSurfaces({SurfaceId{0}, SurfaceId{0}, SurfaceId{2}, SurfaceId{3}, SurfaceId{4}, SurfaceId{5}, SurfaceId{6}});
  std::vector<SurfaceGraph> plan;
  plan.push_back(std::move(planGraph));

  // Gate 4 B3.1: `frame` declared before `tf` so it is constructed first and
  // destroyed last (see SurfaceTrackingScratch's own lifetime-contract doc).
  TimeFrame frame;
  SurfaceTrackingScratch tf;
  auto pool = std::make_shared<BoundedMemoryResource>();
  std::shared_ptr<tbb::task_arena> arena;
  TrackerTraits traits;
  std::vector<TrackingParameters> params(1);
  resetDetectorDefaults(params[0], o2::detectors::DetID::ITS);
  frame.setMemoryPool(pool);
  tf.setMemoryPool(pool);
  traits.setMemoryPool(pool);
  traits.setNThreads(1, arena);
  traits.adoptScratch(&tf);
  traits.adoptFrame(&frame);
  traits.updateTrackingParameters(params);
  traits.setBz(Bz);
  const auto layoutView = plan.front().getView();
  tf.adoptPlan(plan.front().getOrderedSurfaces().size(), layoutView.nTransitions, layoutView.nCells);

  SurfaceMask owned;
  for (const auto surface : plan.front().getOrderedSurfaces()) {
    owned.set(surface);
  }
  const auto bindingResult = SurfacePlanBinding::build(plan.front().getView(), ClusterSourceId{0}, owned,
                                                       plan.front().getOrderedSurfaces(), SurfaceKind::Cylinder);
  BOOST_CHECK(!bindingResult.ok());
}

BOOST_AUTO_TEST_CASE(CombinedCylinderAndDiskLayoutIsRejectedBeforeTrackletProcessing)
{
  // Gate 4 Slice 0a scoping note (accepted plan revision): a layout
  // combining both surface kinds cannot be authored through the production
  // buildSurfaceGraphs() path at all -- it always builds a single
  // subgraph over one detector's own homogeneous-kind orderedSurfaces,
  // matching the architecture note that no cross-detector edges exist
  // because none are authored. This test
  // constructs one directly to prove TrackerTraits::initialiseTimeFrame()
  // still fails closed (MixedSurfaceKindLayout) before the plan binding is
  // committed and before computeLayerTracklets() could process anything --
  // i.e. that no duplicate/cross-tag candidate processing is possible even
  // if such a layout existed. Per-kind span exactness/disjointness itself is
  // already proven by the SurfacePlanBinding schedule tests
  // level by the surface-plan schedule tests; this test is the
  // TrackerTraits-level complement covering the "rejected" case.
  const auto nCylinders = static_cast<uint16_t>(ITSNLayers);
  const auto nDisks = static_cast<uint16_t>(MFTNLayers);
  const std::vector<SurfaceDescriptor> surfaces = combinedCatalog(nCylinders, nDisks);
  const SurfaceCatalogView catalogView{surfaces.data(), static_cast<uint32_t>(surfaces.size())};

  SurfaceGraphBuilder builder{catalogView};
  SurfaceGraphSubgraph cylinderSubgraph;
  cylinderSubgraph.orderedSurfaces = rangeSurfaces(0, nCylinders);
  SurfaceGraphSubgraph diskSubgraph;
  diskSubgraph.orderedSurfaces = rangeSurfaces(nCylinders, nDisks);
  auto result = builder.addSubgraph(std::move(cylinderSubgraph)).addSubgraph(std::move(diskSubgraph)).build();
  BOOST_REQUIRE(result.ok());

  std::vector<SurfaceGraph> plan;
  plan.push_back(std::move(*result.graph));

  // Gate 4 B3.1: `frame` declared before `tf` so it is constructed first and
  // destroyed last (see SurfaceTrackingScratch's own lifetime-contract doc).
  TimeFrame frame;
  SurfaceTrackingScratch tf;
  auto pool = std::make_shared<BoundedMemoryResource>();
  std::shared_ptr<tbb::task_arena> arena;
  TrackerTraits traits;
  std::vector<TrackingParameters> params(1);
  resetDetectorDefaults(params[0], o2::detectors::DetID::ITS);
  frame.setMemoryPool(pool);
  tf.setMemoryPool(pool);
  traits.setMemoryPool(pool);
  traits.setNThreads(1, arena);
  traits.adoptScratch(&tf);
  traits.adoptFrame(&frame);
  traits.updateTrackingParameters(params);
  traits.setBz(Bz);

  SurfaceMask owned;
  for (const auto& surface : surfaces) {
    owned.set(surface.id);
  }
  const auto bindingResult = SurfacePlanBinding::build(plan.front().getView(), ClusterSourceId{0}, owned,
                                                       plan.front().getOrderedSurfaces(), SurfaceKind::Cylinder);
  BOOST_CHECK(!bindingResult.ok());
}
