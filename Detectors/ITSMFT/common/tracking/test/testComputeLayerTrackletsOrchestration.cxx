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
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/Tracklet.h"
#include "MFTTracking/Constants.h"

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
  // existing fixture rather than a separate harness.
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
  // this slice). A later fallible check (here: LayerxX0 corrupted so
  // AttachHitPolicyConfigView::isValid() fails) must leave those arrays
  // exactly zero-filled at the correct size -- never a mixture of computed
  // and zero entries -- because the (non-throwing) value-computation loop
  // this slice added never starts until every fallible check has succeeded.
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame<ITSNLayers> tf;
  TrackerTraits<ITSNLayers> traits;
  std::shared_ptr<tbb::task_arena> arena;
  std::vector<TrackingParameters> params(1);
  resetDetectorDefaults(params[0], o2::detectors::DetID::ITS);
  params[0].PassFlags.reset();
  params[0].PassFlags.set(IterationStep::FirstPass, IterationStep::RebuildClusterLUT);
  params[0].LayerxX0[2] = -1.f; // invalid: AttachHitPolicyConfigView rejects negative xX0

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
