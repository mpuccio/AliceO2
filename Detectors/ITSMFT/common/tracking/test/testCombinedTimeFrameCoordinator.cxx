// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file testCombinedTimeFrameCoordinator.cxx
/// \brief Gate 4 C3: focused tests for the host-only combined ITS+MFT
/// disconnected-tracking coordinator.

#define BOOST_TEST_MODULE ITSMFT CombinedTimeFrameCoordinator
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <cmath>
#include <memory>
#include <optional>
#include <vector>

#include <boost/test/unit_test.hpp>

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
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/CombinedTimeFrameCoordinator.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/LegacyTrackerScratch.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITSMFTTracking/TransitionPolicyOperations.h"
#include "ITStracking/Constants.h"
#include "ReconstructionDataFormats/Track.h"

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

void ensureTrivialMagneticFieldIsSet()
{
  static const bool done = [] {
    TGeoGlobalMagField::Instance()->SetField(new o2::field::MagneticField());
    TGeoGlobalMagField::Instance()->Lock();
    return true;
  }();
  (void)done;
}

std::vector<SurfaceId> ordered(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
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

DecodedCluster diskCluster(float x, float y, float z, int layer)
{
  DecodedCluster cluster{};
  cluster.global = {x, y, z};
  cluster.rowColumnCovariance = {1.e-2f, 0.f, 1.e-2f};
  cluster.sensor = static_cast<uint32_t>(layer);
  cluster.layer = layer;
  return cluster;
}

DecodedCluster cylinderCluster(float radius, float phi, float tanLambda, int layer)
{
  DecodedCluster cluster{};
  cluster.global = {radius * std::cos(phi), radius * std::sin(phi), radius * tanLambda};
  cluster.rowColumnCovariance = {1.e-2f, 0.f, 1.e-2f};
  cluster.sensor = static_cast<uint32_t>(layer);
  cluster.layer = layer;
  return cluster;
}

/// Same chained-projection construction as
/// testComputeLayerTrackletsOrchestration.cxx's buildMftChainClusters() /
/// testDetectorTraversalBindingOrchestration.cxx's identically-named helper:
/// each hop's target is a genuine geometric match via
/// detail::mftTrackletProject, so every adjacent pair in the chain produces
/// a real tracklet, and a full-length chain reaches acceptance.
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

/// A genuine, low-but-nonzero-curvature helical ITS barrel trajectory,
/// sampled at each nominal layer radius via the same standard O2 barrel-
/// propagation utility production ITS/TPC-matching code already uses
/// (o2::track::TrackPar::getXatLabR() to find the local x where the helix
/// crosses a given lab radius, then getXYZGloAt() to read the global point
/// there -- both const, no incremental state mutation between layers).
///
/// A perfectly collinear ("infinite pT" / zero-curvature) triple is a
/// genuine, deliberate rejection of TransitionPolicyOperations.cxx's
/// barrel::buildSeed circle fit -- see testBarrelSurfaceStateOperations.cxx's
/// own BuildSeedDegenerateZeroFieldGeometryRejectsViaNonFiniteOutput and its
/// "three well-separated, non-collinear points" fixture comment. A first
/// attempt at this fixture used exactly such a collinear radial line and
/// reproduced that same rejection (RotationFailure) on one specific
/// three-layer cell, confirming it is not usable as an ITS road fixture;
/// this helix construction is deliberately non-degenerate instead.
std::vector<DecodedCluster> buildItsHelixChainClusters(const std::vector<float>& radii, float bz, float pt, float phi0, float tanl)
{
  const float px = pt * std::cos(phi0);
  const float py = pt * std::sin(phi0);
  const float pz = pt * tanl;
  o2::track::TrackPar seed(std::array<float, 3>{0.f, 0.f, 0.f}, std::array<float, 3>{px, py, pz}, 1, true);

  std::vector<DecodedCluster> clusters;
  clusters.reserve(radii.size());
  for (size_t layer = 0; layer < radii.size(); ++layer) {
    float xAtR = 0.f;
    if (!seed.getXatLabR(radii[layer], xAtR, bz, o2::track::DirType::DirOutward)) {
      return {};
    }
    bool ok = false;
    const auto point = seed.getXYZGloAt(xAtR, bz, ok);
    if (!ok) {
      return {};
    }
    DecodedCluster cluster{};
    cluster.global = {static_cast<float>(point.X()), static_cast<float>(point.Y()), static_cast<float>(point.Z())};
    cluster.rowColumnCovariance = {1.e-2f, 0.f, 1.e-2f};
    cluster.sensor = static_cast<uint32_t>(layer);
    cluster.layer = static_cast<int>(layer);
    clusters.push_back(cluster);
  }
  return clusters;
}

TrackingParameters makeItsParams()
{
  TrackingParameters p;
  resetDetectorDefaults(p, o2::detectors::DetID::ITS);
  // Tracklet formation needs a primary vertex to seed the search window
  // (TrackerTraits.cxx's forTracklets()): with UseDiamond=false (ITS's own
  // default) that must come from TimeFrame::getPrimaryVertices(), which
  // these focused fixtures never populate. UseDiamond=true instead uses the
  // fixed Diamond{0,0,0} vertex every synthetic radial chain below is built
  // through, with no TimeFrame vertex needed -- the same knob
  // buildMftChainClusters()'s MFT fixtures already rely on.
  p.UseDiamond = true;
  return p;
}

TrackingParameters makeMftParams()
{
  TrackingParameters p;
  resetDetectorDefaults(p, o2::detectors::DetID::MFT);
  p.UseDiamond = true;
  p.CreateArtefactLabels = false;
  return p;
}

/// Encodes `decoded` as compact/pattern input and returns a
/// ClusterSourceInput referencing `decoder`/`compactOut`/`patternsOut`/
/// `rofsOut` (kept alive by the caller for the lifetime of every process()
/// call that uses it).
ClusterSourceInput makeSource(ClusterSourceId id, o2::detectors::DetID::ID det, const std::vector<SurfaceId>& surfaces,
                              const PrescribedDecoder& decoder, std::vector<CompClusterExt>& compactOut,
                              std::vector<unsigned char>& patternsOut, std::vector<ROFRecord>& rofsOut,
                              const std::vector<DecodedCluster>& decoded)
{
  compactOut.reserve(decoded.size());
  patternsOut.reserve(decoded.size() * OnePixelPattern.size());
  for (const auto& cluster : decoded) {
    compactOut.emplace_back(0, 0, CompCluster::InvalidPatternID, cluster.layer);
    patternsOut.insert(patternsOut.end(), OnePixelPattern.begin(), OnePixelPattern.end());
  }
  rofsOut = {ROFRecord{{100, 5}, 0, 0, static_cast<int>(compactOut.size())}};

  ClusterSourceInput source{};
  source.id = id;
  source.detector = det;
  source.clusters = compactOut;
  source.patterns = patternsOut;
  source.rofs = rofsOut;
  source.dictionary = &dict();
  source.layerToSurface = surfaces;
  source.timing = ROFTimingConfig{40, 0, 0, 0};
  source.decoder = &decoder;
  return source;
}

/// A source that is valid (dense-empty ROF, zero clusters) but describes no
/// hits at all -- the coordinator's own required "the other detector may be
/// empty" shape, matching ITSMFTTrackingInterface's own zero-cluster path.
ClusterSourceInput makeEmptySource(ClusterSourceId id, o2::detectors::DetID::ID det, const std::vector<SurfaceId>& surfaces,
                                   const PrescribedDecoder& decoder)
{
  ClusterSourceInput source{};
  source.id = id;
  source.detector = det;
  source.dictionary = &dict();
  source.layerToSurface = surfaces;
  source.timing = ROFTimingConfig{40, 0, 0, 0};
  source.decoder = &decoder;
  return source;
}

/// Independent, non-combined, single-detector reference run: the same shape
/// ITSMFTTrackingInterface<NLayers>/testDetectorTraversalBindingOrchestration
/// .cxx's StandaloneMftRun already use -- global SurfaceIds equal compact
/// scratch slots, no DetectorTraversalBinding adopted. Used as the "reproduce
/// the standalone oracle count" reference for the combined coordinator.
template <int NLayers>
struct StandaloneRun {
  TimeFrame frame;
  LegacyTrackerScratch<NLayers> scratch;
  TrackerTraits<NLayers> traits;
  Tracker<NLayers> tracker{&traits};
  std::shared_ptr<tbb::task_arena> arena;
  std::shared_ptr<BoundedMemoryResource> pool = std::make_shared<BoundedMemoryResource>();
  std::vector<TrackingParameters> params;
  std::vector<SurfaceDescriptor> catalog;
  std::optional<DetectorLayoutSet> plan;
  ITSSharedClusterCompatibility itsSidecar;
  MFTPublicationCompatibility mftSidecar;
  TrackingResult result;

  StandaloneRun(o2::detectors::DetID::ID det, SurfaceKind kind, TransitionPolicyTag policy,
                const TrackingParameters& singleParams, const std::vector<DecodedCluster>& decoded)
    : params{singleParams}
  {
    frame.setMemoryPool(pool);
    scratch.setMemoryPool(pool);
    traits.setMemoryPool(pool);
    traits.setNThreads(1, arena);
    traits.adoptScratch(&scratch);
    traits.adoptFrame(&frame);
    traits.updateTrackingParameters(params);
    traits.setBz(Bz);

    const auto orderedSurfaces = ordered(0, NLayers);
    catalog.reserve(NLayers);
    for (uint16_t i = 0; i < NLayers; ++i) {
      SurfaceDescriptor surface{SurfaceId{i}, i, static_cast<uint8_t>(det), kind};
      const float xOverX0 = det == o2::detectors::DetID::MFT ? kNominalMFTLayerX0[i] : kNominalITSLayerX0[i];
      surface.material.xOverX0 = xOverX0;
      surface.material.arealDensityGPerCm2 = xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
      catalog.push_back(surface);
    }
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    auto planResult = buildDetectorLayoutSet(catalogView, orderedSurfaces, policy, params);
    BOOST_REQUIRE(planResult.ok());
    plan.emplace(std::move(*planResult.layout));

    std::vector<CompClusterExt> compact;
    std::vector<unsigned char> patterns;
    for (const auto& cluster : decoded) {
      compact.emplace_back(0, 0, CompCluster::InvalidPatternID, cluster.layer);
      patterns.insert(patterns.end(), OnePixelPattern.begin(), OnePixelPattern.end());
    }
    const std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, static_cast<int>(compact.size())}};
    PrescribedDecoder decoder{det, kind, decoded};
    const auto load = scratch.loadNormalizedSource(frame, decoder, o2::InteractionRecord{50, 5}, ROFTimingConfig{40, 0, 0, 0},
                                                   compact, patterns, rofs, &dict(), nullptr, det,
                                                   gsl::span<const SurfaceId>{plan->getConfigurationKey().orderedSurfaces},
                                                   plan->getSurfaceCatalog());
    BOOST_REQUIRE(load.ok());

    o2::its::LayerTiming layerTiming{};
    layerTiming.mNROFsTF = 1;
    layerTiming.mROFLength = 40;
    typename LegacyTrackerScratch<NLayers>::ROFOverlapTableN rofTable;
    for (int layer = 0; layer < NLayers; ++layer) {
      rofTable.defineLayer(layer, layerTiming);
    }
    rofTable.init();
    scratch.setROFOverlapTable(rofTable);
    typename LegacyTrackerScratch<NLayers>::ROFVertexLookupTableN vtxTable;
    for (int layer = 0; layer < NLayers; ++layer) {
      vtxTable.defineLayer(layer, layerTiming);
    }
    vtxTable.init();
    scratch.setROFVertexLookupTable(vtxTable);
    typename LegacyTrackerScratch<NLayers>::ROFMaskTableN mask{rofTable};
    mask.resetMask();
    for (int layer = 0; layer < NLayers; ++layer) {
      mask.setROFsEnabled(layer, 0, 1, 1);
    }
    scratch.setMultiplicityCutMask(std::move(mask));
    scratch.initTrackerTopologies(params);

    if (det == o2::detectors::DetID::ITS) {
      tracker.adoptITSSharedClusterCompatibility(itsSidecar);
    } else {
      tracker.adoptMFTPublicationCompatibility(mftSidecar);
    }
    tracker.adoptScratch(scratch);
    tracker.adoptFrame(frame);
    tracker.adoptDetectorLayoutSet(*plan);
    tracker.setParameters(params);
    tracker.setMemoryPool(pool);
    tracker.setBz(Bz);
    result = tracker.clustersToTracks();
  }
};

CombinedTimeFrameCoordinator makeCoordinator(const TrackingParameters& itsParams, const TrackingParameters& mftParams)
{
  return CombinedTimeFrameCoordinator{std::vector<TrackingParameters>{itsParams}, std::vector<TrackingParameters>{mftParams}};
}

} // namespace

BOOST_AUTO_TEST_CASE(CombinedLoadingBackfillsIndependentCompactScratches)
{
  // TrackerTraits::findRoads() unconditionally touches the global
  // o2::base::Propagator singleton on first use, regardless of whether any
  // road is actually found -- required before any clustersToTracks() call.
  ensureTrivialMagneticFieldIsSet();
  const auto itsSurfaces = ordered(0, ITSNLayers);
  const auto mftSurfaces = ordered(ITSNLayers, MFTNLayers);
  const auto itsClusters = std::vector<DecodedCluster>{cylinderCluster(3.f, 0.2f, 0.1f, 0), cylinderCluster(4.f, 0.2f, 0.1f, 1)};
  const auto mftClusters = std::vector<DecodedCluster>{diskCluster(1.f, 0.5f, detail::mftLayerZ(0), 0), diskCluster(1.f, 0.5f, detail::mftLayerZ(1), 1)};

  PrescribedDecoder itsDecoder{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, itsClusters};
  PrescribedDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, mftClusters};
  std::vector<CompClusterExt> itsCompact, mftCompact;
  std::vector<unsigned char> itsPatterns, mftPatterns;
  std::vector<ROFRecord> itsRofs, mftRofs;
  const auto itsSource = makeSource(ClusterSourceId{0}, o2::detectors::DetID::ITS, itsSurfaces, itsDecoder, itsCompact, itsPatterns, itsRofs, itsClusters);
  const auto mftSource = makeSource(ClusterSourceId{1}, o2::detectors::DetID::MFT, mftSurfaces, mftDecoder, mftCompact, mftPatterns, mftRofs, mftClusters);

  auto coordinator = makeCoordinator(makeItsParams(), makeMftParams());
  TimeFrame frame;
  coordinator.adoptFrame(frame);
  coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  coordinator.setBz(Bz);
  coordinator.setNThreads(1);

  const auto result = coordinator.process(itsSource, mftSource, o2::InteractionRecord{50, 5});
  BOOST_REQUIRE(result.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::Success);

  BOOST_CHECK_EQUAL(coordinator.getITSScratch().getTotalClusters(), static_cast<int>(itsClusters.size()));
  BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getTotalClusters(), static_cast<int>(mftClusters.size()));
  // Independent compact backfills: each scratch only ever sees its own
  // detector's layer count worth of per-layer arrays, and only its own
  // clusters landed there.
  BOOST_CHECK_EQUAL(coordinator.getITSScratch().getNrof(0), 1);
  BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getNrof(0), 1);
}

BOOST_AUTO_TEST_CASE(MftGlobalIdsWorkEndToEndThroughRefitAndReproducesStandaloneCount)
{
  ensureTrivialMagneticFieldIsSet();
  const auto itsSurfaces = ordered(0, ITSNLayers);
  const auto mftSurfaces = ordered(ITSNLayers, MFTNLayers);

  const auto mftParams = makeMftParams();
  const auto mftClusters = buildMftChainClusters(mftParams, Bz, MFTNLayers - 1);
  BOOST_REQUIRE_EQUAL(mftClusters.size(), static_cast<size_t>(MFTNLayers));

  StandaloneRun<MFTNLayers> standalone{o2::detectors::DetID::MFT, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk, mftParams, mftClusters};
  BOOST_REQUIRE(standalone.result.outcome == TrackingOutcome::Success);
  BOOST_REQUIRE_GT(standalone.scratch.getNumberOfTracks(), 0u);

  PrescribedDecoder itsDecoder{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, {}};
  PrescribedDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, mftClusters};
  std::vector<CompClusterExt> mftCompact;
  std::vector<unsigned char> mftPatterns;
  std::vector<ROFRecord> mftRofs;
  const auto itsSource = makeEmptySource(ClusterSourceId{0}, o2::detectors::DetID::ITS, itsSurfaces, itsDecoder);
  const auto mftSource = makeSource(ClusterSourceId{1}, o2::detectors::DetID::MFT, mftSurfaces, mftDecoder, mftCompact, mftPatterns, mftRofs, mftClusters);

  auto coordinator = makeCoordinator(makeItsParams(), mftParams);
  TimeFrame frame;
  coordinator.adoptFrame(frame);
  coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  coordinator.setBz(Bz);
  coordinator.setNThreads(1);

  const auto result = coordinator.process(itsSource, mftSource, o2::InteractionRecord{50, 5});
  BOOST_REQUIRE(result.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::Success);
  BOOST_CHECK_EQUAL(result.nITSTracks, 0u);
  // Global MFT SurfaceIds 7..16 plus source 1 worked end to end through
  // refit: the combined pass through the coordinator reproduces the
  // standalone (global==compact, unbound) oracle count exactly.
  BOOST_CHECK_EQUAL(result.nMFTTracks, standalone.scratch.getNumberOfTracks());
  BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getNumberOfTracks(), standalone.scratch.getNumberOfTracks());
}

BOOST_AUTO_TEST_CASE(ITSAndMFTAcceptedResultsReproduceStandaloneCountsInOneCombinedPass)
{
  ensureTrivialMagneticFieldIsSet();
  const auto itsSurfaces = ordered(0, ITSNLayers);
  const auto mftSurfaces = ordered(ITSNLayers, MFTNLayers);

  const auto itsParams = makeItsParams();
  const auto mftParams = makeMftParams();
  const auto itsClusters = buildItsHelixChainClusters(itsParams.LayerRadii, Bz, 1.f, 0.4f, 0.3f);
  BOOST_REQUIRE_EQUAL(itsClusters.size(), static_cast<size_t>(ITSNLayers));
  const auto mftClusters = buildMftChainClusters(mftParams, Bz, MFTNLayers - 1);
  BOOST_REQUIRE_EQUAL(mftClusters.size(), static_cast<size_t>(MFTNLayers));

  StandaloneRun<ITSNLayers> standaloneIts{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder, itsParams, itsClusters};
  BOOST_REQUIRE(standaloneIts.result.outcome == TrackingOutcome::Success);
  // A genuine full 7-layer road (MinTrackLength=7, MaxHoles=0): the helix
  // fixture above is a real, non-degenerate curved trajectory, so this is a
  // nonzero accepted-track oracle, not a 0==0 parity check.
  BOOST_REQUIRE_GT(standaloneIts.scratch.getNumberOfTracks(), 0u);
  StandaloneRun<MFTNLayers> standaloneMft{o2::detectors::DetID::MFT, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk, mftParams, mftClusters};
  BOOST_REQUIRE(standaloneMft.result.outcome == TrackingOutcome::Success);
  BOOST_REQUIRE_GT(standaloneMft.scratch.getNumberOfTracks(), 0u);

  PrescribedDecoder itsDecoder{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, itsClusters};
  PrescribedDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, mftClusters};
  std::vector<CompClusterExt> itsCompact, mftCompact;
  std::vector<unsigned char> itsPatterns, mftPatterns;
  std::vector<ROFRecord> itsRofs, mftRofs;
  const auto itsSource = makeSource(ClusterSourceId{0}, o2::detectors::DetID::ITS, itsSurfaces, itsDecoder, itsCompact, itsPatterns, itsRofs, itsClusters);
  const auto mftSource = makeSource(ClusterSourceId{1}, o2::detectors::DetID::MFT, mftSurfaces, mftDecoder, mftCompact, mftPatterns, mftRofs, mftClusters);

  auto coordinator = makeCoordinator(itsParams, mftParams);
  TimeFrame frame;
  coordinator.adoptFrame(frame);
  coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  coordinator.setBz(Bz);
  coordinator.setNThreads(1);

  const auto result = coordinator.process(itsSource, mftSource, o2::InteractionRecord{50, 5});
  BOOST_REQUIRE(result.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::Success);

  // ITS and MFT accepted results reproduce their standalone oracle counts in
  // one combined pass -- nonzero on both sides, not a 0==0 check.
  BOOST_CHECK_GT(result.nITSTracks, 0u);
  BOOST_CHECK_GT(result.nMFTTracks, 0u);
  BOOST_CHECK_EQUAL(result.nITSTracks, standaloneIts.scratch.getNumberOfTracks());
  BOOST_CHECK_EQUAL(result.nMFTTracks, standaloneMft.scratch.getNumberOfTracks());

  // No cross-detector topology element reached either tracker: had the
  // adopted bindings leaked a foreign transition/cell, the combined cell
  // counts below (the real, non-degenerate evidence -- getNumberOfTracklets()
  // is always 0 by the time clustersToTracks() returns: computeLayerCells()
  // clears the per-transition tracklet arrays once it has consumed them,
  // TrackerTraits.cxx) would diverge from the independently-built,
  // single-detector-catalog standalone references.
  BOOST_CHECK_EQUAL(coordinator.getITSScratch().getNumberOfTracklets(), standaloneIts.scratch.getNumberOfTracklets());
  BOOST_CHECK_EQUAL(coordinator.getITSScratch().getNumberOfCells(), standaloneIts.scratch.getNumberOfCells());
  BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getNumberOfTracklets(), standaloneMft.scratch.getNumberOfTracklets());
  BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getNumberOfCells(), standaloneMft.scratch.getNumberOfCells());
  BOOST_CHECK_GT(coordinator.getITSScratch().getNumberOfCells(), 0u);
  BOOST_CHECK_GT(coordinator.getMFTScratch().getNumberOfCells(), 0u);

  // CommonTrack global references resolve correctly and ordering is ITS
  // then MFT: every accepted track's hitSurfaces mask stays within exactly
  // one detector's own global range, and every ITS-range entry precedes
  // every MFT-range entry (shared TimeFrame, append-only, ITS run first).
  const auto itsMask = SurfaceMask{uint32_t{(1u << ITSNLayers) - 1u}};
  const auto mftMask = SurfaceMask{static_cast<uint32_t>(((1u << MFTNLayers) - 1u) << ITSNLayers)};
  const auto& commonTracks = frame.getCommonTracks();
  BOOST_REQUIRE_EQUAL(commonTracks.size(), result.nITSTracks + result.nMFTTracks);
  bool seenMft = false;
  for (size_t i = 0; i < commonTracks.size(); ++i) {
    const auto& track = commonTracks[i];
    BOOST_REQUIRE(track.hitSurfaces.isSubsetOf(itsMask) || track.hitSurfaces.isSubsetOf(mftMask));
    const bool isMft = track.hitSurfaces.isSubsetOf(mftMask) && !track.hitSurfaces.empty();
    if (isMft) {
      seenMft = true;
    } else {
      BOOST_CHECK_MESSAGE(!seenMft, "ITS CommonTrack at index " << i << " appeared after an MFT one");
    }
    for (uint32_t ref = track.firstClusterRef; ref < track.clusterRefEnd; ++ref) {
      const auto& reference = frame.getTrackClusterIndices()[ref];
      const auto* measurement = frame.getNormalizedFrame().getMeasurement(reference.surface, reference.index);
      BOOST_REQUIRE(measurement != nullptr);
      BOOST_CHECK(measurement->surface == reference.surface);
      BOOST_CHECK(isMft ? mftMask.has(reference.surface) : itsMask.has(reference.surface));
    }
  }
  BOOST_CHECK_EQUAL(seenMft, result.nMFTTracks > 0);

  // Publication exports are valid after success, source-qualified, and
  // carry each detector's own ordered-surface span.
  const auto itsExport = coordinator.getITSPublicationExport();
  const auto mftExport = coordinator.getMFTPublicationExport();
  BOOST_REQUIRE(itsExport.has_value());
  BOOST_REQUIRE(mftExport.has_value());
  BOOST_CHECK(itsExport->detector == o2::detectors::DetID::ITS);
  BOOST_CHECK(itsExport->source == ClusterSourceId{0});
  BOOST_CHECK_EQUAL(itsExport->orderedSurfaces.size(), static_cast<size_t>(ITSNLayers));
  BOOST_CHECK(itsExport->orderedSurfaces[0] == SurfaceId{0});
  BOOST_CHECK(mftExport->detector == o2::detectors::DetID::MFT);
  BOOST_CHECK(mftExport->source == ClusterSourceId{1});
  BOOST_CHECK_EQUAL(mftExport->orderedSurfaces.size(), static_cast<size_t>(MFTNLayers));
  BOOST_CHECK(mftExport->orderedSurfaces[0] == SurfaceId{ITSNLayers});
}

BOOST_AUTO_TEST_CASE(LoadFailureResetsWholeCombinedTFExactlyOnceAndInvalidatesPublication)
{
  ensureTrivialMagneticFieldIsSet();
  const auto itsSurfaces = ordered(0, ITSNLayers);
  const auto mftSurfaces = ordered(ITSNLayers, MFTNLayers);
  const auto itsClusters = std::vector<DecodedCluster>{cylinderCluster(3.f, 0.2f, 0.1f, 0), cylinderCluster(4.f, 0.2f, 0.1f, 1)};
  const auto mftClusters = std::vector<DecodedCluster>{diskCluster(1.f, 0.5f, detail::mftLayerZ(0), 0), diskCluster(1.f, 0.5f, detail::mftLayerZ(1), 1)};

  PrescribedDecoder itsDecoder{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, itsClusters};
  PrescribedDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, mftClusters};
  std::vector<CompClusterExt> itsCompact, mftCompact;
  std::vector<unsigned char> itsPatterns, mftPatterns;
  std::vector<ROFRecord> itsRofs, mftRofs;
  const auto itsSource = makeSource(ClusterSourceId{0}, o2::detectors::DetID::ITS, itsSurfaces, itsDecoder, itsCompact, itsPatterns, itsRofs, itsClusters);
  auto mftSource = makeSource(ClusterSourceId{1}, o2::detectors::DetID::MFT, mftSurfaces, mftDecoder, mftCompact, mftPatterns, mftRofs, mftClusters);

  auto coordinator = makeCoordinator(makeItsParams(), makeMftParams());
  TimeFrame frame;
  coordinator.adoptFrame(frame);
  coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  coordinator.setBz(Bz);
  coordinator.setNThreads(1);

  // First pass genuinely succeeds, so there is real state (scratches,
  // CommonTracks, publication exports) for the second, failing pass to
  // actually have to clear.
  const auto first = coordinator.process(itsSource, mftSource, o2::InteractionRecord{50, 5});
  BOOST_REQUIRE(first.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::Success);
  BOOST_REQUIRE(coordinator.getITSPublicationExport().has_value());
  BOOST_REQUIRE(coordinator.getMFTPublicationExport().has_value());

  // Malformed MFT ROF partition (a gap before the second cluster): a
  // structural load failure MultiSourceTimeFrameLoader::loadITSAndMFT()
  // must reject before touching either scratch or the shared TimeFrame.
  std::vector<ROFRecord> malformedMftRofs{ROFRecord{{100, 5}, 0, 0, 1}, ROFRecord{{140, 5}, 0, 2, 1}};
  mftSource.rofs = malformedMftRofs;

  const auto second = coordinator.process(itsSource, mftSource, o2::InteractionRecord{50, 5});
  // MFT's own DropTFUponFailure defaults false (makeMftParams() never sets
  // it), so this recoverable InvalidROFRange load error is still classified
  // Structural -- see CombinedOutcome's own doc.
  BOOST_CHECK(second.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::Structural);
  BOOST_CHECK_EQUAL(second.nITSTracks, 0u);
  BOOST_CHECK_EQUAL(second.nMFTTracks, 0u);

  BOOST_CHECK_EQUAL(coordinator.getITSScratch().getTotalClusters(), 0);
  BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getTotalClusters(), 0);
  BOOST_CHECK_EQUAL(coordinator.getITSScratch().getNumberOfTracks(), 0u);
  BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getNumberOfTracks(), 0u);
  BOOST_CHECK(frame.getCommonTracks().empty());
  BOOST_CHECK(frame.getTrackClusterIndices().empty());
  BOOST_CHECK(!coordinator.getITSPublicationExport().has_value());
  BOOST_CHECK(!coordinator.getMFTPublicationExport().has_value());
}

BOOST_AUTO_TEST_CASE(MFTTrackingFailureAfterITSSuccessStillResetsBothScratches)
{
  ensureTrivialMagneticFieldIsSet();
  const auto itsSurfaces = ordered(0, ITSNLayers);
  const auto mftSurfaces = ordered(ITSNLayers, MFTNLayers);
  const auto itsClusters = std::vector<DecodedCluster>{cylinderCluster(3.f, 0.2f, 0.1f, 0), cylinderCluster(4.f, 0.2f, 0.1f, 1)};
  const auto mftClusters = std::vector<DecodedCluster>{diskCluster(1.f, 0.5f, detail::mftLayerZ(0), 0), diskCluster(1.f, 0.5f, detail::mftLayerZ(1), 1)};

  PrescribedDecoder itsDecoder{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, itsClusters};
  PrescribedDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, mftClusters};
  std::vector<CompClusterExt> itsCompact, mftCompact;
  std::vector<unsigned char> itsPatterns, mftPatterns;
  std::vector<ROFRecord> itsRofs, mftRofs;
  const auto itsSource = makeSource(ClusterSourceId{0}, o2::detectors::DetID::ITS, itsSurfaces, itsDecoder, itsCompact, itsPatterns, itsRofs, itsClusters);
  const auto mftSource = makeSource(ClusterSourceId{1}, o2::detectors::DetID::MFT, mftSurfaces, mftDecoder, mftCompact, mftPatterns, mftRofs, mftClusters);

  // ITS runs (and would succeed) first; MFT's own MaxMemory is exhausted
  // immediately, so its Tracker<MFTNLayers>::clustersToTracks() returns
  // RecoverableDropped. The coordinator must still reset ITS's own scratch
  // -- Tracker<MFTNLayers>'s own internal recovery only resets its own
  // scratch plus the shared TimeFrame, never a sibling detector's scratch.
  auto mftParams = makeMftParams();
  mftParams.MaxMemory = 1;
  mftParams.DropTFUponFailure = true;

  auto coordinator = makeCoordinator(makeItsParams(), mftParams);
  TimeFrame frame;
  coordinator.adoptFrame(frame);
  coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  coordinator.setBz(Bz);
  coordinator.setNThreads(1);

  const auto result = coordinator.process(itsSource, mftSource, o2::InteractionRecord{50, 5});
  // Tracker<NLayers>::clustersToTracks() already gated this resource
  // exhaustion on MFT's own DropTFUponFailure=true internally, so it
  // returned TrackingOutcome::RecoverableDropped rather than throwing; the
  // coordinator's non-Success branch for a tracking-phase outcome is always
  // RecoverableDropped (see CombinedOutcome's own doc).
  BOOST_CHECK(result.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::RecoverableDropped);
  BOOST_CHECK_EQUAL(coordinator.getITSScratch().getTotalClusters(), 0);
  BOOST_CHECK_EQUAL(coordinator.getITSScratch().getNumberOfTracks(), 0u);
  BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getTotalClusters(), 0);
  BOOST_CHECK(frame.getCommonTracks().empty());
  BOOST_CHECK(!coordinator.getITSPublicationExport().has_value());
  BOOST_CHECK(!coordinator.getMFTPublicationExport().has_value());
}

BOOST_AUTO_TEST_CASE(ConstructorRejectsMultiIterationParameters)
{
  std::vector<TrackingParameters> twoIterations(2);
  BOOST_CHECK_THROW((CombinedTimeFrameCoordinator{twoIterations, {makeMftParams()}}), std::invalid_argument);
  BOOST_CHECK_THROW((CombinedTimeFrameCoordinator{{makeItsParams()}, twoIterations}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(BothBindingsDeriveFromOneIdenticalGlobalTopology)
{
  // Single authoritative combined topology: the coordinator builds its one
  // shared ITS+MFT DetectorLayout exactly once (CombinedTimeFrameCoordinator
  // .cxx's buildCombinedLayout(), called once from the constructor) and both
  // DetectorLayoutSets/DetectorTraversalBindings only ever receive a passive
  // copy of that one built object (ownDetectorPlan()). This test proves the
  // two copies are byte-identical in content, including global
  // TransitionId/CellTopologyId identity, so neither can have diverged into
  // an independent authority.
  auto coordinator = makeCoordinator(makeItsParams(), makeMftParams());
  const auto itsView = coordinator.getITSLayoutView();
  const auto mftView = coordinator.getMFTLayoutView();

  // Both bindings were built against the same 17-surface catalog.
  BOOST_REQUIRE_EQUAL(itsView.nSurfaces, mftView.nSurfaces);
  BOOST_REQUIRE_EQUAL(itsView.nSurfaces, static_cast<uint32_t>(ITSNLayers + MFTNLayers));
  BOOST_CHECK(itsView.surfaces == mftView.surfaces); // pointer identity: same static catalog storage
  BOOST_CHECK(itsView.cylinderSurfaces == mftView.cylinderSurfaces);
  BOOST_CHECK(itsView.diskSurfaces == mftView.diskSurfaces);

  // Both bindings were built against content-identical topologies: same
  // transition/cell counts, and every global TransitionId/CellTopologyId
  // resolves to the exact same (from, to, policyTag, skippedSurfaces) /
  // (firstTransition, secondTransition, hitSurfaces) content on both sides.
  const auto& itsTopo = itsView.topology;
  const auto& mftTopo = mftView.topology;
  BOOST_REQUIRE_EQUAL(itsTopo.nTransitions, mftTopo.nTransitions);
  BOOST_REQUIRE_EQUAL(itsTopo.nCells, mftTopo.nCells);
  for (uint32_t t = 0; t < itsTopo.nTransitions; ++t) {
    const auto& a = itsTopo.getTransition(TransitionId{static_cast<uint16_t>(t)});
    const auto& b = mftTopo.getTransition(TransitionId{static_cast<uint16_t>(t)});
    BOOST_CHECK(a.from == b.from);
    BOOST_CHECK(a.to == b.to);
    BOOST_CHECK(a.policyTag == b.policyTag);
    BOOST_CHECK(a.skippedSurfaces == b.skippedSurfaces);
  }
  for (uint32_t c = 0; c < itsTopo.nCells; ++c) {
    const auto& a = itsTopo.getCell(CellTopologyId{static_cast<uint16_t>(c)});
    const auto& b = mftTopo.getCell(CellTopologyId{static_cast<uint16_t>(c)});
    BOOST_CHECK(a.firstTransition == b.firstTransition);
    BOOST_CHECK(a.secondTransition == b.secondTransition);
    BOOST_CHECK(a.hitSurfaces == b.hitSurfaces);
  }

  // No detector-local topology can diverge from the shared one: every
  // surface/transition/cell either binding owns is a genuine subset of this
  // one identical 17-surface topology, and the two owned subsets are
  // disjoint (proving both scope into the same global space rather than
  // each defining its own local one).
  const auto itsMask = SurfaceMask{uint32_t{(1u << ITSNLayers) - 1u}};
  const auto mftMask = SurfaceMask{static_cast<uint32_t>(((1u << MFTNLayers) - 1u) << ITSNLayers)};
  BOOST_CHECK((itsMask & mftMask).empty());
  for (uint32_t t = 0; t < itsTopo.nTransitions; ++t) {
    const auto& transition = itsTopo.getTransition(TransitionId{static_cast<uint16_t>(t)});
    BOOST_CHECK(itsMask.has(transition.from) == itsMask.has(transition.to));
    BOOST_CHECK(mftMask.has(transition.from) == mftMask.has(transition.to));
  }
}

// Gate 4 C4: focused coverage for CombinedOutcome's three-way classification
// (RecoverableDropped / Structural), the always-valid ordered-surface
// getters, the compatibility sidecar getters, and sequential-TF state
// replacement -- the coordinator-side contract the opt-in combined DPL
// workflow (Detectors/ITSMFT/common/workflow-combined-ca/) relies on.
namespace
{

/// A minimal, always-valid ITS+MFT source pair sharing the two-cluster
/// fixture already used by CombinedLoadingBackfillsIndependentCompactScratches.
struct MinimalFixture {
  std::vector<SurfaceId> itsSurfaces = ordered(0, ITSNLayers);
  std::vector<SurfaceId> mftSurfaces = ordered(ITSNLayers, MFTNLayers);
  std::vector<DecodedCluster> itsClusters{cylinderCluster(3.f, 0.2f, 0.1f, 0), cylinderCluster(4.f, 0.2f, 0.1f, 1)};
  std::vector<DecodedCluster> mftClusters{diskCluster(1.f, 0.5f, detail::mftLayerZ(0), 0), diskCluster(1.f, 0.5f, detail::mftLayerZ(1), 1)};
  PrescribedDecoder itsDecoder{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, itsClusters};
  PrescribedDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, mftClusters};
  std::vector<CompClusterExt> itsCompact, mftCompact;
  std::vector<unsigned char> itsPatterns, mftPatterns;
  std::vector<ROFRecord> itsRofs, mftRofs;
  ClusterSourceInput itsSource;
  ClusterSourceInput mftSource;

  MinimalFixture()
  {
    itsSource = makeSource(ClusterSourceId{0}, o2::detectors::DetID::ITS, itsSurfaces, itsDecoder, itsCompact, itsPatterns, itsRofs, itsClusters);
    mftSource = makeSource(ClusterSourceId{1}, o2::detectors::DetID::MFT, mftSurfaces, mftDecoder, mftCompact, mftPatterns, mftRofs, mftClusters);
  }
};

/// A malformed (gap-before-second-cluster) ROF partition for one detector's
/// source, reproducing MultiSourceLoadError::InvalidROFRange -- a
/// *recoverable* per-TF data error under isRecoverableLoadError()
/// (TimeFrameLoadFailure.cxx) -- without touching the other detector's
/// (still valid) source.
void makeRofGap(std::vector<ROFRecord>& rofs)
{
  rofs = {ROFRecord{{100, 5}, 0, 0, 1}, ROFRecord{{140, 5}, 0, 2, 1}};
}

} // namespace

BOOST_AUTO_TEST_CASE(RecoverableITSLoadFailureIsDroppedOnlyWhenITSDropTFAllows)
{
  ensureTrivialMagneticFieldIsSet();

  for (const bool itsDropTF : {true, false}) {
    MinimalFixture fixture;
    makeRofGap(fixture.itsRofs);
    fixture.itsSource.rofs = fixture.itsRofs;

    auto itsParams = makeItsParams();
    itsParams.DropTFUponFailure = itsDropTF;
    auto coordinator = makeCoordinator(itsParams, makeMftParams());
    TimeFrame frame;
    coordinator.adoptFrame(frame);
    coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    coordinator.setBz(Bz);
    coordinator.setNThreads(1);

    const auto result = coordinator.process(fixture.itsSource, fixture.mftSource, o2::InteractionRecord{50, 5});
    const auto expected = itsDropTF ? CombinedTimeFrameCoordinator::CombinedOutcome::RecoverableDropped
                                    : CombinedTimeFrameCoordinator::CombinedOutcome::Structural;
    BOOST_CHECK_MESSAGE(result.outcome == expected, "ITS DropTFUponFailure=" << itsDropTF);
    // Every non-success path still performs exactly one whole reset:
    // both scratches, the shared TimeFrame's CommonTracks, and both
    // publication exports are empty/invalid regardless of classification.
    BOOST_CHECK_EQUAL(coordinator.getITSScratch().getTotalClusters(), 0);
    BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getTotalClusters(), 0);
    BOOST_CHECK(frame.getCommonTracks().empty());
    BOOST_CHECK(!coordinator.getITSPublicationExport().has_value());
    BOOST_CHECK(!coordinator.getMFTPublicationExport().has_value());
  }
}

BOOST_AUTO_TEST_CASE(RecoverableMFTLoadFailureIsDroppedOnlyWhenMFTDropTFAllows)
{
  ensureTrivialMagneticFieldIsSet();

  for (const bool mftDropTF : {true, false}) {
    MinimalFixture fixture;
    makeRofGap(fixture.mftRofs);
    fixture.mftSource.rofs = fixture.mftRofs;

    auto mftParams = makeMftParams();
    mftParams.DropTFUponFailure = mftDropTF;
    auto coordinator = makeCoordinator(makeItsParams(), mftParams);
    TimeFrame frame;
    coordinator.adoptFrame(frame);
    coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    coordinator.setBz(Bz);
    coordinator.setNThreads(1);

    const auto result = coordinator.process(fixture.itsSource, fixture.mftSource, o2::InteractionRecord{50, 5});
    const auto expected = mftDropTF ? CombinedTimeFrameCoordinator::CombinedOutcome::RecoverableDropped
                                    : CombinedTimeFrameCoordinator::CombinedOutcome::Structural;
    BOOST_CHECK_MESSAGE(result.outcome == expected, "MFT DropTFUponFailure=" << mftDropTF);
    BOOST_CHECK_EQUAL(coordinator.getITSScratch().getTotalClusters(), 0);
    BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getTotalClusters(), 0);
    BOOST_CHECK(frame.getCommonTracks().empty());
    BOOST_CHECK(!coordinator.getITSPublicationExport().has_value());
    BOOST_CHECK(!coordinator.getMFTPublicationExport().has_value());
  }
}

BOOST_AUTO_TEST_CASE(StructuralLoadErrorIsAlwaysStructuralRegardlessOfDropTF)
{
  ensureTrivialMagneticFieldIsSet();

  // A missing dictionary is MultiSourceLoadError::MissingDictionary, never
  // recoverable under isRecoverableLoadError() -- DropTFUponFailure=true
  // must not turn this into a dropped TF.
  MinimalFixture fixture;
  fixture.itsSource.dictionary = nullptr;

  auto itsParams = makeItsParams();
  itsParams.DropTFUponFailure = true;
  auto coordinator = makeCoordinator(itsParams, makeMftParams());
  TimeFrame frame;
  coordinator.adoptFrame(frame);
  coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  coordinator.setBz(Bz);
  coordinator.setNThreads(1);

  const auto result = coordinator.process(fixture.itsSource, fixture.mftSource, o2::InteractionRecord{50, 5});
  BOOST_CHECK(result.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::Structural);
  BOOST_CHECK(frame.getCommonTracks().empty());
  BOOST_CHECK(!coordinator.getITSPublicationExport().has_value());
}

BOOST_AUTO_TEST_CASE(UnrecognizedLoadSourceIsAlwaysStructural)
{
  ensureTrivialMagneticFieldIsSet();

  // loadITSAndMFT() rejects any id other than its own fixed ITS=0/MFT=1
  // contract as MultiSourceLoadError::UnsupportedDetector before ever
  // calling loadSources() -- LoadSourcesResult::source then carries the
  // caller's own (unrecognized) id verbatim. Even if a future loader
  // variant ever reported a recoverable error against such an id, this
  // boundary must still classify Structural: an unidentifiable source is
  // never eligible for recoverable/DropTFUponFailure treatment.
  MinimalFixture fixture;
  fixture.itsSource.id = ClusterSourceId{5};

  auto coordinator = makeCoordinator(makeItsParams(), makeMftParams());
  TimeFrame frame;
  coordinator.adoptFrame(frame);
  coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  coordinator.setBz(Bz);
  coordinator.setNThreads(1);

  const auto result = coordinator.process(fixture.itsSource, fixture.mftSource, o2::InteractionRecord{50, 5});
  BOOST_CHECK(result.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::Structural);
  BOOST_CHECK(frame.getCommonTracks().empty());
}

BOOST_AUTO_TEST_CASE(StructuralTrackingExceptionIsClassifiedStructuralAfterWholeReset)
{
  ensureTrivialMagneticFieldIsSet();

  // MaxMemory=1 with DropTFUponFailure left at its false default: MFT's
  // Tracker<MFTNLayers>::clustersToTracks() throws
  // BoundedMemoryResource::MemoryLimitExceeded (a std::bad_alloc, hence a
  // std::exception) rather than returning RecoverableDropped, since that
  // internal gate is only "drop" when DropTFUponFailure is explicitly true
  // (see MFTTrackingFailureAfterITSSuccessStillResetsBothScratches for the
  // DropTFUponFailure=true counterpart).
  MinimalFixture fixture;
  auto mftParams = makeMftParams();
  mftParams.MaxMemory = 1;

  auto coordinator = makeCoordinator(makeItsParams(), mftParams);
  TimeFrame frame;
  coordinator.adoptFrame(frame);
  coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  coordinator.setBz(Bz);
  coordinator.setNThreads(1);

  const auto result = coordinator.process(fixture.itsSource, fixture.mftSource, o2::InteractionRecord{50, 5});
  BOOST_CHECK(result.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::Structural);
  BOOST_CHECK_EQUAL(coordinator.getITSScratch().getTotalClusters(), 0);
  BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getTotalClusters(), 0);
  BOOST_CHECK(frame.getCommonTracks().empty());
  BOOST_CHECK(!coordinator.getITSPublicationExport().has_value());
  BOOST_CHECK(!coordinator.getMFTPublicationExport().has_value());
}

BOOST_AUTO_TEST_CASE(SequentialSuccessfulTFsReplaceStateWithoutStaleAccumulation)
{
  ensureTrivialMagneticFieldIsSet();
  const auto itsSurfaces = ordered(0, ITSNLayers);
  const auto mftSurfaces = ordered(ITSNLayers, MFTNLayers);
  const auto itsParams = makeItsParams();
  const auto mftParams = makeMftParams();
  // A genuine nonzero-track fixture (same construction as
  // ITSAndMFTAcceptedResultsReproduceStandaloneCountsInOneCombinedPass): if
  // CommonTrack/TrackClusterIndices storage ever accumulated across TFs
  // instead of being replaced, the second TF's count below would silently
  // double rather than reproduce the same per-TF value.
  const auto itsClusters = buildItsHelixChainClusters(itsParams.LayerRadii, Bz, 1.f, 0.4f, 0.3f);
  BOOST_REQUIRE_EQUAL(itsClusters.size(), static_cast<size_t>(ITSNLayers));
  const auto mftClusters = buildMftChainClusters(mftParams, Bz, MFTNLayers - 1);
  BOOST_REQUIRE_EQUAL(mftClusters.size(), static_cast<size_t>(MFTNLayers));

  PrescribedDecoder itsDecoder{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, itsClusters};
  PrescribedDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, mftClusters};
  std::vector<CompClusterExt> itsCompact, mftCompact;
  std::vector<unsigned char> itsPatterns, mftPatterns;
  std::vector<ROFRecord> itsRofs, mftRofs;
  const auto itsSource = makeSource(ClusterSourceId{0}, o2::detectors::DetID::ITS, itsSurfaces, itsDecoder, itsCompact, itsPatterns, itsRofs, itsClusters);
  const auto mftSource = makeSource(ClusterSourceId{1}, o2::detectors::DetID::MFT, mftSurfaces, mftDecoder, mftCompact, mftPatterns, mftRofs, mftClusters);

  auto coordinator = makeCoordinator(itsParams, mftParams);
  TimeFrame frame;
  coordinator.adoptFrame(frame);
  coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  coordinator.setBz(Bz);
  coordinator.setNThreads(1);

  const auto firstResult = coordinator.process(itsSource, mftSource, o2::InteractionRecord{50, 5});
  BOOST_REQUIRE(firstResult.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::Success);
  BOOST_REQUIRE_GT(firstResult.nITSTracks + firstResult.nMFTTracks, 0u);
  const auto firstCommonTrackCount = frame.getCommonTracks().size();
  BOOST_REQUIRE_EQUAL(firstCommonTrackCount, firstResult.nITSTracks + firstResult.nMFTTracks);

  // No explicit reset between successful TFs: MultiSourceTimeFrameLoader::
  // loadITSAndMFT()'s commitNormalizedFrame() atomically replaces the
  // normalized frame and clears mCommonTracks/mTrackClusterIndices in the
  // same commit (TimeFrame.h), so the second process() call alone -- on the
  // identical fixture again -- must reproduce the same per-TF count, not
  // the first TF's count plus the second's.
  const auto secondResult = coordinator.process(itsSource, mftSource, o2::InteractionRecord{60, 6});
  BOOST_REQUIRE(secondResult.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::Success);

  BOOST_CHECK_EQUAL(secondResult.nITSTracks, firstResult.nITSTracks);
  BOOST_CHECK_EQUAL(secondResult.nMFTTracks, firstResult.nMFTTracks);
  BOOST_CHECK_EQUAL(frame.getCommonTracks().size(), firstCommonTrackCount);
  BOOST_CHECK_EQUAL(coordinator.getITSScratch().getTotalClusters(), static_cast<int>(itsClusters.size()));
  BOOST_CHECK_EQUAL(coordinator.getMFTScratch().getTotalClusters(), static_cast<int>(mftClusters.size()));
}

BOOST_AUTO_TEST_CASE(OrderedSurfaceGettersAreAlwaysValidUnlikePublicationExports)
{
  auto coordinator = makeCoordinator(makeItsParams(), makeMftParams());

  // Valid immediately after construction -- no adoptFrame()/process() call
  // needed, unlike getITSPublicationExport()/getMFTPublicationExport().
  const auto itsSurfacesBefore = coordinator.getITSOrderedSurfaces();
  const auto mftSurfacesBefore = coordinator.getMFTOrderedSurfaces();
  BOOST_REQUIRE_EQUAL(itsSurfacesBefore.size(), static_cast<size_t>(ITSNLayers));
  BOOST_REQUIRE_EQUAL(mftSurfacesBefore.size(), static_cast<size_t>(MFTNLayers));
  BOOST_CHECK(itsSurfacesBefore[0] == SurfaceId{0});
  BOOST_CHECK(mftSurfacesBefore[0] == SurfaceId{ITSNLayers});
  BOOST_CHECK(!coordinator.getITSPublicationExport().has_value());
  BOOST_CHECK(!coordinator.getMFTPublicationExport().has_value());

  // Still identical after a failure (which invalidates the publication
  // exports but must never move the fixed catalog-offset spans).
  ensureTrivialMagneticFieldIsSet();
  MinimalFixture fixture;
  makeRofGap(fixture.mftRofs);
  fixture.mftSource.rofs = fixture.mftRofs;
  TimeFrame frame;
  coordinator.adoptFrame(frame);
  coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  coordinator.setBz(Bz);
  coordinator.setNThreads(1);
  const auto failed = coordinator.process(fixture.itsSource, fixture.mftSource, o2::InteractionRecord{50, 5});
  BOOST_REQUIRE(failed.outcome != CombinedTimeFrameCoordinator::CombinedOutcome::Success);
  BOOST_CHECK(coordinator.getITSOrderedSurfaces().data() == itsSurfacesBefore.data());
  BOOST_CHECK(coordinator.getMFTOrderedSurfaces().data() == mftSurfacesBefore.data());
}

BOOST_AUTO_TEST_CASE(CompatibilitySidecarGettersReflectSealAndReset)
{
  ensureTrivialMagneticFieldIsSet();
  MinimalFixture fixture;
  auto coordinator = makeCoordinator(makeItsParams(), makeMftParams());
  TimeFrame frame;
  coordinator.adoptFrame(frame);
  coordinator.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  coordinator.setBz(Bz);
  coordinator.setNThreads(1);

  // Not yet sealed before any process() call.
  BOOST_CHECK(!coordinator.getITSSharedClusterCompatibility().isSealed());

  const auto result = coordinator.process(fixture.itsSource, fixture.mftSource, o2::InteractionRecord{50, 5});
  BOOST_REQUIRE(result.outcome == CombinedTimeFrameCoordinator::CombinedOutcome::Success);
  // A successful run always seals the ITS sidecar (Tracker<ITSNLayers>::
  // clustersToTracks() -> markTracks() -> sealFromMarkedTracks()), which is
  // exactly what stageITSCommonTrackOutput() requires
  // (CommonTrackOutputAdapter.h).
  BOOST_CHECK(coordinator.getITSSharedClusterCompatibility().isSealed());

  // A whole reset clears both sidecars back to their pre-process() state.
  makeRofGap(fixture.mftRofs);
  fixture.mftSource.rofs = fixture.mftRofs;
  const auto failed = coordinator.process(fixture.itsSource, fixture.mftSource, o2::InteractionRecord{60, 6});
  BOOST_REQUIRE(failed.outcome != CombinedTimeFrameCoordinator::CombinedOutcome::Success);
  BOOST_CHECK(!coordinator.getITSSharedClusterCompatibility().isSealed());
  BOOST_CHECK(coordinator.getITSSharedClusterCompatibility().entries().empty());
  BOOST_CHECK(coordinator.getMFTPublicationCompatibility().find(0, 0) == nullptr);
}
