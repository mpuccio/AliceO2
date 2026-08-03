// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file testDetectorTraversalBindingOrchestration.cxx
/// \brief Gate 4 C2 Slice 1: TrackerTraits<MFTNLayers> hot-loop orchestration
/// bound to a DetectorTraversalBinding built against the combined 17-surface
/// ITS+MFT catalog (global MFT SurfaceIds 7..16, ClusterSourceId{1}).
///
/// Proves: (1) tracklet/cell content is byte-identical between an MFT-only
/// production-shaped run (global == compact, today's Gate 3 behavior) and
/// the same clusters run through a binding-scoped combined layout (global
/// offset by ITSNLayers, translated to the same compact scratch slots) --
/// this is the slot-parity and no-boundary-traversal evidence; (2) a
/// mismatched binding/instance pairing fails closed before any scratch
/// access, never silently misreads or corrupts state.

#define BOOST_TEST_MODULE ITSMFT DetectorTraversalBinding orchestration
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <array>
#include <memory>
#include <utility>
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
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/detail/DetectorTraversalBinding.h"
#include "ITSMFTTracking/MFTFwdTrackHelpers.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "ITSMFTTracking/LegacyTrackerScratch.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/Constants.h"

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

// Every layer eligible as a road start, matching resetDetectorDefaults()'s
// own StartLayerMask default ((1u << nLayers) - 1u, i.e. "any layer may
// seed a road") -- StandaloneMftRun gets this for free through
// buildDetectorLayoutSet()'s positionalSurfaceMask(parameters.StartLayerMask,
// ...); a manually-built combined DetectorLayoutBuilder subgraph has no such
// implicit default and must pass it explicitly, or findRoads()'s
// roadStartCells span is permanently empty (topology.seedingSurfaces.has(...)
// can never be true) regardless of chain length or cuts.
SurfaceMask itsSeedMask() { return SurfaceMask{uint32_t{0x7f}}; }
SurfaceMask mftSeedMask() { return SurfaceMask{uint32_t{0x1ff80}}; }

// TrackerTraits::findRoads() unconditionally touches the global
// o2::base::Propagator singleton on first use, which requires
// TGeoGlobalMagField to already hold a real o2::field::MagneticField object
// -- see testCATrackerFailureContract.cxx's identical helper for the full
// rationale. Only the findRoads()-reaching tests below need this.
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

// Same per-surface xOverX0 convention as testComputeLayerTrackletsOrchestration.cxx's
// makeCatalog(), generalized to an arbitrary global SurfaceId offset so ITS
// and MFT can share one combined catalog while each individually still
// matches o2::itsmft::resetDetectorDefaults()'s per-detector LayerxX0.
void appendMaterialCatalog(std::vector<SurfaceDescriptor>& out, uint16_t firstId, uint16_t nLayers,
                           o2::detectors::DetID::ID detector, SurfaceKind kind)
{
  for (uint16_t i = 0; i < nLayers; ++i) {
    SurfaceDescriptor surface{SurfaceId{static_cast<uint16_t>(firstId + i)}, i, static_cast<uint8_t>(detector), kind};
    const float xOverX0 = detector == o2::detectors::DetID::MFT ? kNominalMFTLayerX0[i % MFTNLayers] : kNominalITSLayerX0[i % ITSNLayers];
    surface.material.xOverX0 = xOverX0;
    surface.material.arealDensityGPerCm2 = xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
    out.push_back(surface);
  }
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

/// Same chained-projection construction as testComputeLayerTrackletsOrchestration.cxx's
/// buildMftChainClusters(): each hop's target is a genuine geometric match
/// via detail::mftTrackletProject, so every adjacent pair in the chain
/// produces a real tracklet.
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

ClusterSourceInput emptyITSSource(const std::vector<SurfaceId>& itsSurfaces, const PrescribedDecoder& decoder)
{
  ClusterSourceInput source{};
  source.id = ClusterSourceId{0};
  source.detector = o2::detectors::DetID::ITS;
  source.dictionary = &dict();
  source.layerToSurface = itsSurfaces;
  source.timing = ROFTimingConfig{40, 0, 0, 0};
  source.decoder = &decoder;
  return source;
}

/// Encodes `decoded` as compact/pattern input and returns a ClusterSourceInput
/// referencing `decoder`/`compactOut`/`patternsOut` (kept alive by the
/// caller for the lifetime of the load call).
ClusterSourceInput mftSource(const std::vector<SurfaceId>& mftSurfaces, const PrescribedDecoder& decoder,
                             std::vector<CompClusterExt>& compactOut, std::vector<unsigned char>& patternsOut,
                             std::vector<ROFRecord>& rofsOut, const std::vector<DecodedCluster>& decoded)
{
  compactOut.reserve(decoded.size());
  patternsOut.reserve(decoded.size() * OnePixelPattern.size());
  for (const auto& cluster : decoded) {
    compactOut.emplace_back(0, 0, CompCluster::InvalidPatternID, cluster.layer);
    patternsOut.insert(patternsOut.end(), OnePixelPattern.begin(), OnePixelPattern.end());
  }
  rofsOut = {ROFRecord{{100, 5}, 0, 0, static_cast<int>(compactOut.size())}};

  ClusterSourceInput source{};
  source.id = ClusterSourceId{1};
  source.detector = o2::detectors::DetID::MFT;
  source.clusters = compactOut;
  source.patterns = patternsOut;
  source.rofs = rofsOut;
  source.dictionary = &dict();
  source.layerToSurface = mftSurfaces;
  source.timing = ROFTimingConfig{40, 0, 0, 0};
  source.decoder = &decoder;
  return source;
}

void setUpMftScratchTables(LegacyTrackerScratchMFT& scratch, const std::vector<TrackingParameters>& params)
{
  o2::its::LayerTiming layerTiming{};
  layerTiming.mNROFsTF = 1;
  layerTiming.mROFLength = 40;
  LegacyTrackerScratchMFT::ROFOverlapTableN rofTable;
  for (int layer = 0; layer < MFTNLayers; ++layer) {
    rofTable.defineLayer(layer, layerTiming);
  }
  rofTable.init();
  scratch.setROFOverlapTable(rofTable);
  LegacyTrackerScratchMFT::ROFVertexLookupTableN vtxTable;
  for (int layer = 0; layer < MFTNLayers; ++layer) {
    vtxTable.defineLayer(layer, layerTiming);
  }
  vtxTable.init();
  scratch.setROFVertexLookupTable(vtxTable);
  LegacyTrackerScratchMFT::ROFMaskTableN mask{rofTable};
  mask.resetMask();
  for (int layer = 0; layer < MFTNLayers; ++layer) {
    mask.setROFsEnabled(layer, 0, 1, 1);
  }
  scratch.setMultiplicityCutMask(std::move(mask));
  scratch.initTrackerTopologies(params);
}

std::vector<TrackingParameters> mftParams()
{
  std::vector<TrackingParameters> params(1);
  resetDetectorDefaults(params[0], o2::detectors::DetID::MFT);
  params[0].UseDiamond = true;
  params[0].CreateArtefactLabels = false;
  params[0].PassFlags.reset();
  params[0].PassFlags.set(IterationStep::FirstPass, IterationStep::RebuildClusterLUT);
  return params;
}

/// Standalone (no binding), production-shaped MFT-only baseline: global ==
/// compact, today's Gate 3 behavior unchanged by this slice.
struct StandaloneMftRun {
  TimeFrame frame;
  LegacyTrackerScratchMFT scratch;
  TrackerTraits<MFTNLayers> traits;
  std::shared_ptr<tbb::task_arena> arena;
  std::shared_ptr<BoundedMemoryResource> pool = std::make_shared<BoundedMemoryResource>();
  std::vector<TrackingParameters> params = mftParams();
  // Declared before `plan`: DetectorLayoutSet borrows a non-owning
  // SurfaceCatalogView into this vector (never a copy), so it must outlive
  // `plan` -- C++ destroys members in reverse declaration order, so this
  // member is destroyed last. A constructor-local vector here would leave
  // `plan`'s view dangling the instant the constructor returns (only latent
  // as long as nothing re-reads the catalog after construction, e.g. a
  // second initialiseTimeFrame() call from Tracker::clustersToTracks()).
  std::vector<SurfaceDescriptor> catalog;

  StandaloneMftRun(const std::vector<DecodedCluster>& decoded)
  {
    frame.setMemoryPool(pool);
    scratch.setMemoryPool(pool);
    traits.setMemoryPool(pool);
    traits.setNThreads(1, arena);
    traits.adoptScratch(&scratch);
    traits.adoptFrame(&frame);
    traits.updateTrackingParameters(params);
    traits.setBz(Bz);

    const auto orderedSurfaces = ordered(0, MFTNLayers);
    appendMaterialCatalog(catalog, 0, MFTNLayers, o2::detectors::DetID::MFT, SurfaceKind::Disk);
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    auto planResult = buildDetectorLayoutSet(catalogView, orderedSurfaces, params);
    BOOST_REQUIRE(planResult.ok());
    plan.emplace(std::move(*planResult.layout));

    std::vector<CompClusterExt> compact;
    std::vector<unsigned char> patterns;
    compact.reserve(decoded.size());
    patterns.reserve(decoded.size() * OnePixelPattern.size());
    for (const auto& cluster : decoded) {
      compact.emplace_back(0, 0, CompCluster::InvalidPatternID, cluster.layer);
      patterns.insert(patterns.end(), OnePixelPattern.begin(), OnePixelPattern.end());
    }
    const std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, static_cast<int>(compact.size())}};
    PrescribedDecoder decoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, decoded};
    const auto load = scratch.loadNormalizedSource(frame, decoder, o2::InteractionRecord{50, 5}, ROFTimingConfig{40, 0, 0, 0},
                                                    compact, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::MFT,
                                                    gsl::span<const SurfaceId>{plan->getConfigurationKey().orderedSurfaces}, plan->getSurfaceCatalog());
    BOOST_REQUIRE(load.ok());

    setUpMftScratchTables(scratch, params);
    traits.initialiseTimeFrame(0, *plan);
  }

  std::optional<DetectorLayoutSet> plan;
};

/// Combined-catalog run: ITS occupies global SurfaceIds [0, ITSNLayers), MFT
/// occupies [ITSNLayers, ITSNLayers+MFTNLayers) in one shared 17-surface
/// DetectorLayout/topology, exactly as a future combined coordinator would
/// build it. The MFT TrackerTraits adopts a DetectorTraversalBinding scoping
/// it to its own owned global surfaces/transitions/cells; ITS contributes an
/// empty source (zero clusters) so this test exercises the MFT/binding path
/// in isolation without needing a second TrackerTraits instance.
struct CombinedMftRun {
  TimeFrame frame;
  LegacyTrackerScratchITS itsScratch;
  LegacyTrackerScratchMFT mftScratch;
  TrackerTraits<MFTNLayers> traits;
  std::shared_ptr<tbb::task_arena> arena;
  std::shared_ptr<BoundedMemoryResource> pool = std::make_shared<BoundedMemoryResource>();
  std::vector<TrackingParameters> params = mftParams();
  std::vector<SurfaceDescriptor> catalog;
  std::unique_ptr<DetectorTraversalBinding> binding;
  std::optional<DetectorLayoutSet> plan;

  CombinedMftRun(const std::vector<DecodedCluster>& decoded)
  {
    frame.setMemoryPool(pool);
    itsScratch.setMemoryPool(pool);
    mftScratch.setMemoryPool(pool);
    traits.setMemoryPool(pool);
    traits.setNThreads(1, arena);
    traits.adoptScratch(&mftScratch);
    traits.adoptFrame(&frame);
    traits.updateTrackingParameters(params);
    traits.setBz(Bz);

    appendMaterialCatalog(catalog, 0, ITSNLayers, o2::detectors::DetID::ITS, SurfaceKind::Cylinder);
    appendMaterialCatalog(catalog, ITSNLayers, MFTNLayers, o2::detectors::DetID::MFT, SurfaceKind::Disk);
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};

    const auto itsSurfaces = ordered(0, ITSNLayers);
    const auto mftSurfaces = ordered(ITSNLayers, MFTNLayers);
    DetectorLayoutBuilder builder{catalogView};
    builder.addSubgraph({itsSurfaces, 0, SurfaceMask{}, itsSeedMask()});
    builder.addSubgraph({mftSurfaces, 0, SurfaceMask{}, mftSeedMask()});
    auto built = builder.build();
    BOOST_REQUIRE(built.ok());
    DetectorLayout combinedLayout = std::move(*built.layout);
    const gsl::span<const SurfaceDescriptor> catalogSpan{catalogView.surfaces, catalogView.nSurfaces};
    const auto masks = computeSurfaceKindMasks(catalogSpan);
    const auto combinedView = combinedLayout.getView(catalogSpan, masks.first, masks.second);

    auto bindingResult = DetectorTraversalBinding::build(combinedView, o2::detectors::DetID::MFT, ClusterSourceId{1},
                                                          SurfaceMask{uint32_t{0x1ff80}}, mftSurfaces);
    BOOST_REQUIRE(bindingResult.ok());
    binding = std::move(bindingResult.binding);
    traits.adoptDetectorTraversalBinding(binding.get());

    DetectorLayoutConfigurationKey key;
    key.orderedSurfaces = mftSurfaces;
    key.iterations.push_back(DetectorLayoutIterationConfiguration{static_cast<uint32_t>(MFTNLayers), 0, LayerMask{0}, LayerMask{0}});
    std::vector<DetectorLayout> layouts;
    layouts.push_back(std::move(combinedLayout));
    plan.emplace(std::move(key), catalogView, std::move(layouts));

    PrescribedDecoder itsDecoder{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, {}};
    PrescribedDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, decoded};
    std::vector<CompClusterExt> compact;
    std::vector<unsigned char> patterns;
    std::vector<ROFRecord> rofs;
    const auto itsSource = emptyITSSource(itsSurfaces, itsDecoder);
    const auto mftSrc = mftSource(mftSurfaces, mftDecoder, compact, patterns, rofs, decoded);
    const auto load = MultiSourceTimeFrameLoader::loadITSAndMFT(frame, itsScratch, mftScratch, itsSource, mftSrc,
                                                                 catalogView, o2::InteractionRecord{50, 5});
    BOOST_REQUIRE(load.ok());

    setUpMftScratchTables(mftScratch, params);
    traits.initialiseTimeFrame(0, *plan);
  }
};

} // namespace

BOOST_AUTO_TEST_CASE(MftBoundBindingReproducesMftOnlyCatalogSingleTracklet)
{
  TrackingParameters probeParams;
  resetDetectorDefaults(probeParams, o2::detectors::DetID::MFT);
  const float fromZ = detail::mftLayerZ(0);
  const float toZ = detail::mftLayerZ(1);
  float targetX = 0.f, targetY = 0.f;
  detail::mftTrackletProject(1.f, 0.5f, fromZ, probeParams.Diamond[0], probeParams.Diamond[1], probeParams.Diamond[2],
                             0, 1, Bz, probeParams.TrackletMinPt, targetX, targetY);
  const std::vector<DecodedCluster> clusters{diskCluster(1.f, 0.5f, fromZ, 0), diskCluster(targetX, targetY, toZ, 1)};

  StandaloneMftRun standalone{clusters};
  CombinedMftRun combined{clusters};

  standalone.traits.computeLayerTracklets(0, 0);
  combined.traits.computeLayerTracklets(0, 0);

  // Global SurfaceIds 7,8 (MFT layers 0,1 inside the combined 17-surface
  // catalog) must resolve, through the adopted binding, to exactly the same
  // compact scratch slot 0 that the MFT-only catalog (global == compact)
  // already uses -- this is the slot-parity/no-boundary-traversal proof:
  // if the binding leaked a foreign id or misassigned a slot, these two
  // independently-built LegacyTrackerScratch<MFTNLayers> instances would
  // disagree.
  const auto& standaloneTracklets = standalone.scratch.getTracklets();
  const auto& combinedTracklets = combined.mftScratch.getTracklets();
  BOOST_REQUIRE_EQUAL(standaloneTracklets.size(), combinedTracklets.size());
  BOOST_REQUIRE_GE(standaloneTracklets.size(), 1u);
  BOOST_REQUIRE_EQUAL(standaloneTracklets[0].size(), 1u);
  BOOST_REQUIRE_EQUAL(combinedTracklets[0].size(), 1u);
  BOOST_CHECK(standaloneTracklets[0][0] == combinedTracklets[0][0]);
  BOOST_CHECK_EQUAL(standaloneTracklets[0][0].tanLambda, combinedTracklets[0][0].tanLambda);
  BOOST_CHECK_EQUAL(standaloneTracklets[0][0].phi, combinedTracklets[0][0].phi);
  for (size_t id = 1; id < combinedTracklets.size(); ++id) {
    BOOST_CHECK(combinedTracklets[id].empty());
  }
}

BOOST_AUTO_TEST_CASE(MftBoundBindingMultiHopChainReproducesCellsAndNeighbours)
{
  auto params = mftParams();
  const auto clusters = buildMftChainClusters(params[0], Bz, 3);
  BOOST_REQUIRE_EQUAL(clusters.size(), 4u);

  StandaloneMftRun standalone{clusters};
  CombinedMftRun combined{clusters};

  standalone.traits.computeLayerTracklets(0, 0);
  standalone.traits.computeLayerCells(0);
  standalone.traits.findCellsNeighbours(0);
  combined.traits.computeLayerTracklets(0, 0);
  combined.traits.computeLayerCells(0);
  combined.traits.findCellsNeighbours(0);

  BOOST_REQUIRE_EQUAL(standalone.scratch.getCells().size(), combined.mftScratch.getCells().size());
  size_t standaloneNonEmptyCells = 0, combinedNonEmptyCells = 0;
  size_t standaloneCellCount = 0, combinedCellCount = 0;
  for (size_t id = 0; id < standalone.scratch.getCells().size(); ++id) {
    standaloneNonEmptyCells += !standalone.scratch.getCells()[id].empty();
    combinedNonEmptyCells += !combined.mftScratch.getCells()[id].empty();
    standaloneCellCount += standalone.scratch.getCells()[id].size();
    combinedCellCount += combined.mftScratch.getCells()[id].size();
  }
  // Non-degenerate: the 3-hop chain must actually produce cells (two
  // transitions per cell, so at least 2 of this fixture's 3 transitions
  // participate) for this to be a meaningful cross-check at all.
  BOOST_CHECK_GT(standaloneCellCount, 0u);
  BOOST_CHECK_EQUAL(standaloneNonEmptyCells, combinedNonEmptyCells);
  BOOST_CHECK_EQUAL(standaloneCellCount, combinedCellCount);

  // findCellsNeighbours' dynamically-discovered-neighbour translation site
  // (TrackerTraits.cxx, the nextTopologyId -> requireScratchCellSlot() site):
  // a 3-hop chain produces at least one cell-to-cell neighbour relation
  // (cell on layers 0-1-2 feeding cell on layers 1-2-3), reflected in a
  // non-degenerate CellsNeighboursLUT for whichever compact slot the source
  // cell occupies. Cross-checked in aggregate between the two runs, exactly
  // as the cell counts above.
  size_t standaloneNeighbourEntries = 0, combinedNeighbourEntries = 0;
  for (size_t id = 0; id < standalone.scratch.getCellsNeighbours().size(); ++id) {
    standaloneNeighbourEntries += standalone.scratch.getCellsNeighbours()[id].size();
    combinedNeighbourEntries += combined.mftScratch.getCellsNeighbours()[id].size();
  }
  BOOST_CHECK_GT(standaloneNeighbourEntries, 0u);
  BOOST_CHECK_EQUAL(standaloneNeighbourEntries, combinedNeighbourEntries);
}

BOOST_AUTO_TEST_CASE(MismatchedBindingFailsClosedBeforeTrackletProcessing)
{
  // A binding built for ITS (CylinderCylinder, global surfaces 0..6) adopted
  // into a TrackerTraits<MFTNLayers> instance (DiskDisk-shaped) is a
  // construction-time misuse -- the same class of bug
  // TraversalFailureReason::TraversalBindingMismatch exists to catch -- and
  // must fail closed before any scratch container is touched, never
  // silently misread a foreign transition/cell.
  std::vector<SurfaceDescriptor> catalog;
  appendMaterialCatalog(catalog, 0, ITSNLayers, o2::detectors::DetID::ITS, SurfaceKind::Cylinder);
  appendMaterialCatalog(catalog, ITSNLayers, MFTNLayers, o2::detectors::DetID::MFT, SurfaceKind::Disk);
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  const auto itsSurfaces = ordered(0, ITSNLayers);
  const auto mftSurfaces = ordered(ITSNLayers, MFTNLayers);

  DetectorLayoutBuilder builder{catalogView};
  builder.addSubgraph({itsSurfaces, 0, SurfaceMask{}, itsSeedMask()});
  builder.addSubgraph({mftSurfaces, 0, SurfaceMask{}, mftSeedMask()});
  auto built = builder.build();
  BOOST_REQUIRE(built.ok());
  DetectorLayout combinedLayout = std::move(*built.layout);
  const gsl::span<const SurfaceDescriptor> catalogSpan{catalogView.surfaces, catalogView.nSurfaces};
  const auto masks = computeSurfaceKindMasks(catalogSpan);
  const auto combinedView = combinedLayout.getView(catalogSpan, masks.first, masks.second);

  auto itsBindingResult = DetectorTraversalBinding::build(combinedView, o2::detectors::DetID::ITS, ClusterSourceId{0},
                                                           SurfaceMask{uint32_t{0x7f}}, itsSurfaces);
  BOOST_REQUIRE(itsBindingResult.ok());

  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  LegacyTrackerScratchMFT scratch;
  TrackerTraits<MFTNLayers> traits;
  std::shared_ptr<tbb::task_arena> arena;
  auto params = mftParams();
  frame.setMemoryPool(pool);
  scratch.setMemoryPool(pool);
  traits.setMemoryPool(pool);
  traits.setNThreads(1, arena);
  traits.adoptScratch(&scratch);
  traits.adoptFrame(&frame);
  traits.updateTrackingParameters(params);
  traits.setBz(Bz);
  // Deliberate misuse: an ITS-shaped binding adopted onto an MFT-shaped
  // TrackerTraits.
  traits.adoptDetectorTraversalBinding(itsBindingResult.binding.get());

  DetectorLayoutConfigurationKey key;
  key.orderedSurfaces = mftSurfaces;
  key.iterations.push_back(DetectorLayoutIterationConfiguration{static_cast<uint32_t>(MFTNLayers), 0, LayerMask{0}, LayerMask{0}});
  std::vector<DetectorLayout> layouts;
  layouts.push_back(std::move(combinedLayout));
  const DetectorLayoutSet plan{std::move(key), catalogView, std::move(layouts)};

  PrescribedDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, {diskCluster(1.f, 0.5f, detail::mftLayerZ(0), 0)}};
  std::vector<CompClusterExt> compact{CompClusterExt(0, 0, CompCluster::InvalidPatternID, 0)};
  std::vector<unsigned char> patterns(OnePixelPattern.begin(), OnePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, 1}};
  const auto load = scratch.loadNormalizedSource(frame, mftDecoder, o2::InteractionRecord{50, 5}, ROFTimingConfig{40, 0, 0, 0},
                                                 compact, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::MFT,
                                                 gsl::span<const SurfaceId>{plan.getConfigurationKey().orderedSurfaces}, plan.getSurfaceCatalog());
  BOOST_REQUIRE(load.ok());
  setUpMftScratchTables(scratch, params);

  BOOST_CHECK_EXCEPTION(traits.initialiseTimeFrame(0, plan), TraversalException, [](const TraversalException&) { return true; });
  BOOST_CHECK(!traits.hasTraversalCache());
}

// Gate 4 C2 Slice 2: the same construction-time-misuse scenario as
// MismatchedBindingFailsClosedBeforeTrackletProcessing above, but driven
// through Tracker<MFTNLayers>::clustersToTracks() with
// DropTFUponFailure=true -- proving the binding-mismatch failure (a
// TraversalException, structural by TrackerTraits.cxx's own classification)
// is never reclassified as TrackingOutcome::RecoverableDropped merely
// because the policy flag happens to be set. clustersToTracks()'s
// TraversalException catch block never inspects DropTFUponFailure at all
// (CATracker.cxx), so this is a direct exercise of that existing,
// unconditional behavior against the one new failure reason Slice 1 added.
BOOST_AUTO_TEST_CASE(MismatchedBindingNeverReclassifiesStructuralAsRecoverableDropped)
{
  std::vector<SurfaceDescriptor> catalog;
  appendMaterialCatalog(catalog, 0, ITSNLayers, o2::detectors::DetID::ITS, SurfaceKind::Cylinder);
  appendMaterialCatalog(catalog, ITSNLayers, MFTNLayers, o2::detectors::DetID::MFT, SurfaceKind::Disk);
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  const auto itsSurfaces = ordered(0, ITSNLayers);
  const auto mftSurfaces = ordered(ITSNLayers, MFTNLayers);

  DetectorLayoutBuilder builder{catalogView};
  builder.addSubgraph({itsSurfaces, 0, SurfaceMask{}, itsSeedMask()});
  builder.addSubgraph({mftSurfaces, 0, SurfaceMask{}, mftSeedMask()});
  auto built = builder.build();
  BOOST_REQUIRE(built.ok());
  DetectorLayout combinedLayout = std::move(*built.layout);
  const gsl::span<const SurfaceDescriptor> catalogSpan{catalogView.surfaces, catalogView.nSurfaces};
  const auto masks = computeSurfaceKindMasks(catalogSpan);
  const auto combinedView = combinedLayout.getView(catalogSpan, masks.first, masks.second);

  auto itsBindingResult = DetectorTraversalBinding::build(combinedView, o2::detectors::DetID::ITS, ClusterSourceId{0},
                                                          SurfaceMask{uint32_t{0x7f}}, itsSurfaces);
  BOOST_REQUIRE(itsBindingResult.ok());

  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame frame;
  LegacyTrackerScratchMFT scratch;
  TrackerTraits<MFTNLayers> traits;
  Tracker<MFTNLayers> tracker(&traits);
  std::shared_ptr<tbb::task_arena> arena;
  auto params = mftParams();
  params[0].DropTFUponFailure = true;
  frame.setMemoryPool(pool);
  scratch.setMemoryPool(pool);
  traits.setMemoryPool(pool);
  traits.setNThreads(1, arena);
  traits.adoptScratch(&scratch);
  traits.adoptFrame(&frame);
  traits.updateTrackingParameters(params);
  traits.setBz(Bz);
  // Deliberate misuse: an ITS-shaped binding adopted onto an MFT-shaped
  // TrackerTraits, exactly as the initialiseTimeFrame()-level negative test
  // above.
  traits.adoptDetectorTraversalBinding(itsBindingResult.binding.get());

  tracker.adoptScratch(scratch);
  tracker.adoptFrame(frame);
  tracker.setParameters(params);
  tracker.setMemoryPool(pool);
  tracker.setBz(Bz);

  DetectorLayoutConfigurationKey key;
  key.orderedSurfaces = mftSurfaces;
  key.iterations.push_back(DetectorLayoutIterationConfiguration{static_cast<uint32_t>(MFTNLayers), 0, LayerMask{0}, LayerMask{0}});
  std::vector<DetectorLayout> layouts;
  layouts.push_back(std::move(combinedLayout));
  const DetectorLayoutSet plan{std::move(key), catalogView, std::move(layouts)};
  tracker.adoptDetectorLayoutSet(plan);

  PrescribedDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, {diskCluster(1.f, 0.5f, detail::mftLayerZ(0), 0)}};
  std::vector<CompClusterExt> compact{CompClusterExt(0, 0, CompCluster::InvalidPatternID, 0)};
  std::vector<unsigned char> patterns(OnePixelPattern.begin(), OnePixelPattern.end());
  const std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, 1}};
  const auto load = scratch.loadNormalizedSource(frame, mftDecoder, o2::InteractionRecord{50, 5}, ROFTimingConfig{40, 0, 0, 0},
                                                 compact, patterns, rofs, &dict(), nullptr, o2::detectors::DetID::MFT,
                                                 gsl::span<const SurfaceId>{plan.getConfigurationKey().orderedSurfaces}, plan.getSurfaceCatalog());
  BOOST_REQUIRE(load.ok());
  setUpMftScratchTables(scratch, params);

  BOOST_CHECK_EXCEPTION(tracker.clustersToTracks(), TraversalException, [](const TraversalException&) { return true; });
}

// Gate 4 C2 Slice 2: closes the Gate 4 C2 Slice 1 handoff's stated gap --
// findRoads()'s road-start translation site
// (mBinding->getGlobalRoadStartCells(), TrackerTraits.cxx), and the complete
// Tracker<MFTNLayers>::clustersToTracks() traversal through the combined-
// catalog binding, were previously only exercised structurally (through
// DetectorTraversalBinding's own getGlobalRoadStartCells() test in
// testDetectorTraversalBinding.cxx), never end to end. This drives the same
// 3-hop chain fixture MftBoundBindingMultiHopChainReproducesCellsAndNeighbours
// uses all the way through findRoads()'s refit/road-building stage via
// Tracker<MFTNLayers>::clustersToTracks(), and cross-checks the resulting
// track count against the identical standalone MFT-only-catalog run.
BOOST_AUTO_TEST_CASE(MftBoundBindingFullClustersToTracksReproducesStandaloneTrackCount)
{
  ensureTrivialMagneticFieldIsSet();
  auto chainParams = mftParams();
  // A full-length, all-10-layer chain -- unlike
  // MftBoundBindingMultiHopChainReproducesCellsAndNeighbours's 3-hop probe
  // (sufficient for computeLayerCells()/findCellsNeighbours(), which do not
  // filter by level/length), findRoads() only accepts a road whose level
  // reaches TrackingParameters::CellMinimumLevel() (5 for MFT's own
  // unmodified default MinTrackLength=7, Configuration.h) -- reachable only
  // with enough chained cells, never by lowering the (untouched, default)
  // cuts themselves.
  const auto clusters = buildMftChainClusters(chainParams[0], Bz, MFTNLayers - 1);
  BOOST_REQUIRE_EQUAL(clusters.size(), static_cast<size_t>(MFTNLayers));

  StandaloneMftRun standalone{clusters};
  CombinedMftRun combined{clusters};

  // Production always adopts this before an accepted track can be shadow-
  // published (ITSMFTTrackingInterface<MFTNLayers>'s always-owned
  // MFTPublicationCompatibilityOwner sidecar) -- acceptTracks() treats a
  // missing sidecar as a fatal invariant violation the instant it actually
  // has a track to publish (AcceptedTrackShadowPublisher<MFTNLayers>::
  // publish() returns nullopt with mSidecar == nullptr), which every earlier
  // test in this file never reached (no test before this one produced a
  // real accepted track).
  MFTPublicationCompatibility standaloneSidecar;
  MFTPublicationCompatibility combinedSidecar;

  Tracker<MFTNLayers> standaloneTracker(&standalone.traits);
  standaloneTracker.adoptScratch(standalone.scratch);
  standaloneTracker.adoptFrame(standalone.frame);
  standaloneTracker.adoptDetectorLayoutSet(*standalone.plan);
  standaloneTracker.adoptMFTPublicationCompatibility(standaloneSidecar);
  standaloneTracker.setParameters(standalone.params);
  standaloneTracker.setMemoryPool(standalone.pool);
  standaloneTracker.setBz(Bz);
  const auto standaloneResult = standaloneTracker.clustersToTracks();
  BOOST_REQUIRE(standaloneResult.outcome == TrackingOutcome::Success);

  Tracker<MFTNLayers> combinedTracker(&combined.traits);
  combinedTracker.adoptScratch(combined.mftScratch);
  combinedTracker.adoptFrame(combined.frame);
  combinedTracker.adoptDetectorLayoutSet(*combined.plan);
  combinedTracker.adoptMFTPublicationCompatibility(combinedSidecar);
  combinedTracker.setParameters(combined.params);
  combinedTracker.setMemoryPool(combined.pool);
  combinedTracker.setBz(Bz);
  const auto combinedResult = combinedTracker.clustersToTracks();
  BOOST_REQUIRE(combinedResult.outcome == TrackingOutcome::Success);

  BOOST_CHECK_GT(standalone.scratch.getNumberOfTracks(), 0u);
  // Gate 4 C2 source-identity correction: findRoadsForPolicy() now resolves
  // expectedSource once from mBinding->getSource() (ClusterSourceId{1} for
  // this combined-catalog binding) and threads it through
  // DetectorTraits::refitSeed()/refitTrackFwd() into
  // MFTFwdTrackHelpers.cxx::refitTrackFwdImpl()'s normalized-measurement
  // identity re-check, replacing the previous hardcoded ClusterSourceId{0}
  // comparison that rejected every genuinely-source-1 candidate seed. With
  // that fixed, the binding-translation path (already proven correct through
  // road-building by every earlier test in this file) now also reproduces
  // the standalone track count end to end.
  BOOST_CHECK_EQUAL(standalone.scratch.getNumberOfTracks(), combined.mftScratch.getNumberOfTracks());
}

// Gate 4 C2 source-identity correction: a binding whose expected source does
// not match the genuinely-loaded measurement source (ClusterSourceId{1} for
// combined MFT, per MultiSourceTimeFrameLoader::loadITSAndMFT) must fail
// closed before any track is published. In practice this is caught even
// earlier than the refit boundary this correction touches:
// initialiseTimeFrame()'s own normalized-measurement validation (Gate 3
// Stage-B / Gate 4 C2 Slice 1, TrackerTraits.cxx) already compares every
// mLayerMeasurements entry's source against this same
// mBinding->getSource() -- since refitTrackFwdImpl()'s expectedSource
// (this correction) is resolved from that identical binding, a mismatch
// can never reach the refit stage undetected; it surfaces as a
// NormalizedMeasurementMismatch TraversalException out of
// clustersToTracks() first. That TraversalException path is structural
// (CATracker.cxx never reclassifies it as TrackingOutcome::
// RecoverableDropped, regardless of DropTFUponFailure), so this still
// proves the required end-to-end property: a foreign-source binding never
// reaches AcceptedTrackShadowPublisher, whichever validation layer catches
// it first. Built for ClusterSourceId{2} (a source neither ITS nor MFT
// ever legitimately uses in this fixture) instead of the correct
// ClusterSourceId{1}.
BOOST_AUTO_TEST_CASE(ForeignExpectedSourceBindingFailsClosedBeforePublication)
{
  ensureTrivialMagneticFieldIsSet();
  auto chainParams = mftParams();
  const auto clusters = buildMftChainClusters(chainParams[0], Bz, MFTNLayers - 1);
  BOOST_REQUIRE_EQUAL(clusters.size(), static_cast<size_t>(MFTNLayers));

  std::vector<SurfaceDescriptor> catalog;
  appendMaterialCatalog(catalog, 0, ITSNLayers, o2::detectors::DetID::ITS, SurfaceKind::Cylinder);
  appendMaterialCatalog(catalog, ITSNLayers, MFTNLayers, o2::detectors::DetID::MFT, SurfaceKind::Disk);
  const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
  const auto itsSurfaces = ordered(0, ITSNLayers);
  const auto mftSurfaces = ordered(ITSNLayers, MFTNLayers);

  DetectorLayoutBuilder builder{catalogView};
  builder.addSubgraph({itsSurfaces, 0, SurfaceMask{}, itsSeedMask()});
  builder.addSubgraph({mftSurfaces, 0, SurfaceMask{}, mftSeedMask()});
  auto built = builder.build();
  BOOST_REQUIRE(built.ok());
  DetectorLayout combinedLayout = std::move(*built.layout);
  const gsl::span<const SurfaceDescriptor> catalogSpan{catalogView.surfaces, catalogView.nSurfaces};
  const auto masks = computeSurfaceKindMasks(catalogSpan);
  const auto combinedView = combinedLayout.getView(catalogSpan, masks.first, masks.second);

  // Deliberately foreign expected source: ClusterSourceId{2}, not the
  // genuinely-loaded MFT source (ClusterSourceId{1}) this fixture uses below.
  auto bindingResult = DetectorTraversalBinding::build(combinedView, o2::detectors::DetID::MFT, ClusterSourceId{2},
                                                       SurfaceMask{uint32_t{0x1ff80}}, mftSurfaces);
  BOOST_REQUIRE(bindingResult.ok());
  auto binding = std::move(bindingResult.binding);

  TimeFrame frame;
  LegacyTrackerScratchITS itsScratch;
  LegacyTrackerScratchMFT mftScratch;
  TrackerTraits<MFTNLayers> traits;
  std::shared_ptr<tbb::task_arena> arena;
  auto pool = std::make_shared<BoundedMemoryResource>();
  auto params = mftParams();
  frame.setMemoryPool(pool);
  itsScratch.setMemoryPool(pool);
  mftScratch.setMemoryPool(pool);
  traits.setMemoryPool(pool);
  traits.setNThreads(1, arena);
  traits.adoptScratch(&mftScratch);
  traits.adoptFrame(&frame);
  traits.updateTrackingParameters(params);
  traits.setBz(Bz);
  traits.adoptDetectorTraversalBinding(binding.get());

  DetectorLayoutConfigurationKey key;
  key.orderedSurfaces = mftSurfaces;
  key.iterations.push_back(DetectorLayoutIterationConfiguration{static_cast<uint32_t>(MFTNLayers), 0, LayerMask{0}, LayerMask{0}});
  std::vector<DetectorLayout> layouts;
  layouts.push_back(std::move(combinedLayout));
  DetectorLayoutSet plan{std::move(key), catalogView, std::move(layouts)};

  PrescribedDecoder itsDecoder{o2::detectors::DetID::ITS, SurfaceKind::Cylinder, {}};
  PrescribedDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk, clusters};
  std::vector<CompClusterExt> compact;
  std::vector<unsigned char> patterns;
  std::vector<ROFRecord> rofs;
  const auto itsSource = emptyITSSource(itsSurfaces, itsDecoder);
  const auto mftSrc = mftSource(mftSurfaces, mftDecoder, compact, patterns, rofs, clusters);
  const auto load = MultiSourceTimeFrameLoader::loadITSAndMFT(frame, itsScratch, mftScratch, itsSource, mftSrc,
                                                              catalogView, o2::InteractionRecord{50, 5});
  BOOST_REQUIRE(load.ok());

  setUpMftScratchTables(mftScratch, params);

  MFTPublicationCompatibility sidecar;
  Tracker<MFTNLayers> tracker(&traits);
  tracker.adoptScratch(mftScratch);
  tracker.adoptFrame(frame);
  tracker.adoptDetectorLayoutSet(plan);
  tracker.adoptMFTPublicationCompatibility(sidecar);
  tracker.setParameters(params);
  tracker.setMemoryPool(pool);
  tracker.setBz(Bz);

  // The genuine source (1) never matches the binding's foreign expected
  // source (2): fails closed with a structural TraversalException before
  // any track reaches publication, never a crash or a silently accepted
  // foreign-source track.
  BOOST_CHECK_EXCEPTION(tracker.clustersToTracks(), TraversalException, [](const TraversalException&) { return true; });
  BOOST_CHECK_EQUAL(mftScratch.getNumberOfTracks(), 0u);
}
