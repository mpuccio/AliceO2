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
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/Constants.h"
#include "ITStracking/MathUtils.h"
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

class FakeCatalogProvider final : public DetectorSurfaceCatalogProvider
{
 public:
  explicit FakeCatalogProvider(std::vector<SurfaceDescriptor> catalog) : mCatalog{std::move(catalog)} {}

  DetectorSurfaceCatalogResult buildCatalog(const DetectorSurfaceCatalogRequest&) const final
  {
    return {mCatalog, DetectorSurfaceCatalogError::None};
  }

 private:
  std::vector<SurfaceDescriptor> mCatalog;
};

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

  o2::itsmft::ioutils::SurfaceMeasurementDecodeResult decode(
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
      o2::itsmft::ioutils::SurfaceMeasurementDecodeResult result;
      result.error = clusterData.error;
      return result;
    }

    o2::itsmft::ioutils::SurfaceMeasurementDecodeResult result;
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
/// float for DiskDisk vs `0.5` double-promoted for CylinderCylinder, per the
/// integration review finding preserved -- not canonicalized -- in part 1/4
/// of this slice).
template <int NLayers>
void computeLegacyTransitionMSAndPhiCut(const TrackingParameters& trkParam, float bz, bool isDisk,
                                        const typename TimeFrame<NLayers>::TrackingTopologyN::View& topology,
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
    const auto& transition = topology.getTransition(transitionId);
    float ms2 = 0.f;
    for (int layer = transition.fromLayer; layer < transition.toLayer; ++layer) {
      ms2 += o2::its::math_utils::Sq(msAngles[layer]);
    }
    const float msAngle = o2::gpu::CAMath::Sqrt(ms2);
    const float r1 = trkParam.LayerRadii[transition.fromLayer];
    const float r2 = trkParam.LayerRadii[transition.toLayer];
    if (isDisk) {
      oneOverR = (0.5f * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
    } else {
      oneOverR = (0.5 * oneOverR >= 1.f / r2) ? (2.f / r2) - o2::constants::math::Almost0 : oneOverR;
    }
    const float res1 = o2::gpu::CAMath::Hypot(trkParam.PVres, positionResolution[transition.fromLayer]);
    const float res2 = o2::gpu::CAMath::Hypot(trkParam.PVres, positionResolution[transition.toLayer]);
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
                            TransitionPolicyTag tag,
                            std::vector<DecodedCluster> decoded,
                            int nThreads)
{
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame<NLayers> tf;
  TrackerTraits<NLayers> traits;
  std::shared_ptr<tbb::task_arena> arena;
  std::vector<TrackingParameters> params(1);
  resetDetectorDefaults(params[0], detector);
  params[0].UseDiamond = true;
  params[0].CreateArtefactLabels = false;
  params[0].PassFlags.reset();
  params[0].PassFlags.set(IterationStep::FirstPass, IterationStep::RebuildClusterLUT);

  tf.setMemoryPool(pool);
  traits.setMemoryPool(pool);
  traits.setNThreads(nThreads, arena);
  traits.adoptTimeFrame(&tf);
  traits.updateTrackingParameters(params);
  traits.setBz(Bz);

  const auto orderedSurfaces = identitySurfaces(static_cast<uint16_t>(NLayers));
  FakeCatalogProvider provider{makeCatalog(static_cast<uint16_t>(NLayers), detector, kind)};
  const DetectorSurfaceCatalogRequest catalogRequest{detector, SurfaceId{0}, static_cast<uint32_t>(NLayers)};
  BOOST_REQUIRE(tf.ensureDetectorLayouts(&provider, catalogRequest, orderedSurfaces, tag, params).ok());
  tf.initTrackerTopologies(params);

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
  const auto load = tf.loadNormalizedSource(decoder, o2::InteractionRecord{50, 5}, ROFTimingConfig{40, 0, 0, 0},
                                            compactClusters, patterns, rofs, &dict(), nullptr, detector);
  BOOST_REQUIRE(load.ok());

  o2::its::LayerTiming layerTiming{};
  layerTiming.mNROFsTF = 1;
  layerTiming.mROFLength = 40;
  typename TimeFrame<NLayers>::ROFOverlapTableN rofTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    rofTable.defineLayer(layer, layerTiming);
  }
  rofTable.init();
  tf.setROFOverlapTable(rofTable);
  typename TimeFrame<NLayers>::ROFMaskTableN mask{rofTable};
  mask.resetMask();
  for (int layer = 0; layer < NLayers; ++layer) {
    mask.setROFsEnabled(layer, 0, 1, 1);
  }
  tf.setMultiplicityCutMask(std::move(mask));

  traits.initialiseTimeFrame(0);
  BOOST_CHECK_EQUAL(traits.getTraversalGroupingCount(), 1);
  BOOST_CHECK_EQUAL(traits.getPolicyBindingCount(tag), 1);
  const auto otherTag = tag == TransitionPolicyTag::CylinderCylinder ? TransitionPolicyTag::DiskDisk : TransitionPolicyTag::CylinderCylinder;
  BOOST_CHECK_EQUAL(traits.getPolicyBindingCount(otherTag), 0);

  // Gate 3 transition-preparation slice: successful initialisation must fill
  // every transition entry (relocated from TimeFrame::initialise() into
  // TrackerTraits::initialiseTimeFrame(), see TransitionPolicyOperations.h).
  // Exercised here for both CylinderCylinder and DiskDisk through the
  // existing fixture rather than a separate harness. Beyond finiteness, each
  // entry is checked bit-for-bit against computeLegacyTransitionMSAndPhiCut's
  // independent oracle -- the only replay-grade acceptance evidence for the
  // common CylinderCylinder path, since no real-geometry common-CA ITS
  // replay exists yet.
  {
    const auto preparedTopology = tf.getTrackingTopologyView();
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

  const auto topology = tf.getTrackingTopologyView();
  int transitionId = -1;
  for (int id = 0; id < topology.nTransitions; ++id) {
    const auto& transition = topology.getTransition(id);
    if (transition.fromLayer == 0 && transition.toLayer == 1) {
      transitionId = id;
      break;
    }
  }
  BOOST_REQUIRE_GE(transitionId, 0);

  // Repeated per-vertex processing must reuse the initialization-time
  // traversal grouping and typed parameter binding.
  traits.computeLayerTracklets(0, 0);
  BOOST_CHECK_EQUAL(traits.getTraversalGroupingCount(), 1);
  BOOST_CHECK_EQUAL(traits.getPolicyBindingCount(tag), 1);
  traits.computeLayerTracklets(0, 0);
  BOOST_CHECK_EQUAL(traits.getTraversalGroupingCount(), 1);
  BOOST_CHECK_EQUAL(traits.getPolicyBindingCount(tag), 1);

  TrackletSnapshot result;
  result.transitionId = transitionId;
  result.expectedTimestamp = tf.getROFOverlapTableView().getTimeStamp(0, 0, 1, 0);
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
  return result;
}

void checkSame(const TrackletSnapshot& serial, const TrackletSnapshot& parallel)
{
  BOOST_CHECK_EQUAL(serial.transitionId, parallel.transitionId);
  BOOST_REQUIRE_EQUAL(serial.tracklets.size(), parallel.tracklets.size());
  BOOST_CHECK_EQUAL_COLLECTIONS(serial.lookup.begin(), serial.lookup.end(), parallel.lookup.begin(), parallel.lookup.end());
  for (size_t i = 0; i < serial.tracklets.size(); ++i) {
    BOOST_CHECK(serial.tracklets[i] == parallel.tracklets[i]);
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

} // namespace

BOOST_AUTO_TEST_CASE(CylinderOnePassAndTwoPassProduceIdenticalTracklets)
{
  const std::vector<DecodedCluster> clusters{
    cylinderCluster(3.f, 0.3f, 0),
    cylinderCluster(4.f, 0.4f, 1)};
  const auto serial = runFixture<ITSNLayers>(o2::detectors::DetID::ITS, SurfaceKind::Cylinder,
                                             TransitionPolicyTag::CylinderCylinder, clusters, 1);
  const auto parallel = runFixture<ITSNLayers>(o2::detectors::DetID::ITS, SurfaceKind::Cylinder,
                                               TransitionPolicyTag::CylinderCylinder, clusters, 4);
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
                                             TransitionPolicyTag::DiskDisk, clusters, 1);
  const auto parallel = runFixture<MFTNLayers>(o2::detectors::DetID::MFT, SurfaceKind::Disk,
                                               TransitionPolicyTag::DiskDisk, clusters, 4);
  const float expectedTanLambda = (fromZ - toZ) / (toZ - fromZ);
  const float expectedPhi = o2::gpu::CAMath::ATan2(0.5f - targetY, 1.f - targetX);
  checkExactTracklet(serial, expectedTanLambda, expectedPhi);
  checkExactTracklet(parallel, expectedTanLambda, expectedPhi);
  checkSame(serial, parallel);
}

BOOST_AUTO_TEST_CASE(ComputeLayerTrackletsFailsClosedWithoutInitialiseTimeFrame)
{
  TimeFrame<ITSNLayers> tf;
  TrackerTraits<ITSNLayers> traits;
  traits.adoptTimeFrame(&tf);

  BOOST_CHECK_EXCEPTION(traits.computeLayerTracklets(0, -1), TraversalException, [](const TraversalException& error) {
    return error.getReason() == TraversalFailureReason::InvalidTraversalSchedule;
  });
}

BOOST_AUTO_TEST_CASE(InitialiseTimeFrameFailureLeavesTransitionArraysZeroFilledNotPartial)
{
  // Gate 3 transition-preparation slice failure contract: TimeFrame::initialise()
  // already clears/resizes mTransitionMSAngles/mTransitionPhiCuts to
  // nTransitions before any policy/geometry validation runs (unchanged by
  // this slice). A later fallible check (here: an invalid CorrType, so
  // AttachHitPolicyConfigView::isValid() fails) must leave those arrays
  // exactly zero-filled at the correct size -- never a mixture of computed
  // and zero entries -- because the (non-throwing) value-computation loop
  // this slice added never starts until every fallible check has succeeded.
  // Deliberately not a corrupted LayerxX0 value: since the compatibility-only
  // LegacyMaterialMismatch check (TrackerTraits::initialiseTimeFrame()) now
  // runs before TimeFrame::initialise(), a corrupted LayerxX0 would fail
  // there instead, before the arrays are ever resized -- this test wants a
  // fallible check that fires strictly after TimeFrame::initialise().
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame<ITSNLayers> tf;
  TrackerTraits<ITSNLayers> traits;
  std::shared_ptr<tbb::task_arena> arena;
  std::vector<TrackingParameters> params(1);
  resetDetectorDefaults(params[0], o2::detectors::DetID::ITS);
  params[0].PassFlags.reset();
  params[0].PassFlags.set(IterationStep::FirstPass, IterationStep::RebuildClusterLUT);
  params[0].CorrType = static_cast<o2::base::PropagatorImpl<float>::MatCorrType>(99); // invalid: AttachHitPolicyConfigView rejects it

  tf.setMemoryPool(pool);
  traits.setMemoryPool(pool);
  traits.setNThreads(1, arena);
  traits.adoptTimeFrame(&tf);
  traits.updateTrackingParameters(params);
  traits.setBz(Bz);

  const auto orderedSurfaces = identitySurfaces(static_cast<uint16_t>(ITSNLayers));
  FakeCatalogProvider provider{makeCatalog(static_cast<uint16_t>(ITSNLayers), o2::detectors::DetID::ITS, SurfaceKind::Cylinder)};
  const DetectorSurfaceCatalogRequest catalogRequest{o2::detectors::DetID::ITS, SurfaceId{0}, static_cast<uint32_t>(ITSNLayers)};
  BOOST_REQUIRE(tf.ensureDetectorLayouts(&provider, catalogRequest, orderedSurfaces, TransitionPolicyTag::CylinderCylinder, params).ok());
  tf.initTrackerTopologies(params);

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
  const auto load = tf.loadNormalizedSource(decoder, o2::InteractionRecord{50, 5}, ROFTimingConfig{40, 0, 0, 0},
                                            compactClusters, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::ITS);
  BOOST_REQUIRE(load.ok());

  o2::its::LayerTiming layerTiming{};
  layerTiming.mNROFsTF = 1;
  layerTiming.mROFLength = 40;
  TimeFrame<ITSNLayers>::ROFOverlapTableN rofTable;
  for (int layer = 0; layer < ITSNLayers; ++layer) {
    rofTable.defineLayer(layer, layerTiming);
  }
  rofTable.init();
  tf.setROFOverlapTable(rofTable);
  TimeFrame<ITSNLayers>::ROFMaskTableN mask{rofTable};
  mask.resetMask();
  for (int layer = 0; layer < ITSNLayers; ++layer) {
    mask.setROFsEnabled(layer, 0, 1, 1);
  }
  tf.setMultiplicityCutMask(std::move(mask));

  BOOST_CHECK_EXCEPTION(traits.initialiseTimeFrame(0), TraversalException, [](const TraversalException& error) {
    return error.getReason() == TraversalFailureReason::InvalidPolicyParameters;
  });

  const auto topology = tf.getTrackingTopologyView();
  const auto& msAngles = tf.getTransitionMSAngles();
  const auto& phiCuts = tf.getTransitionPhiCuts();
  BOOST_REQUIRE_EQUAL(msAngles.size(), static_cast<size_t>(topology.nTransitions));
  BOOST_REQUIRE_EQUAL(phiCuts.size(), static_cast<size_t>(topology.nTransitions));
  for (int id = 0; id < topology.nTransitions; ++id) {
    BOOST_CHECK_EQUAL(msAngles[id], 0.f);
    BOOST_CHECK_EQUAL(phiCuts[id], 0.f);
  }
}
