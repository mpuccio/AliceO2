// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 3 index-table activation: TrackerTraits<NLayers>::initialiseTimeFrame()
// binding/commit/reuse lifecycle, exercised directly (no Tracker<N>, no task
// arena, no magnetic field -- initialiseTimeFrame() never reaches
// computeLayerTracklets()/findRoads()). The DropTFUponFailure/wipe gating
// contract for the new TraversalFailureReason values is covered separately
// in testTrackerFailureContract.cxx, which needs the full Tracker<N>::
// clustersToTracks() path.
//
// Contract under test (IndexTableConfiguration.h, TrackerTraits.cxx,
// TimeFrame.cxx):
//  - A FirstPass iteration binds+validates this iteration's TrackingParameters
//    into a local scratch IndexTableUtils and, only if every check passes,
//    commits it into the TimeFrame by value and (re)allocates the LUT from
//    it.
//  - An invalid binding throws TraversalException{InvalidIndexTableConfiguration}
//    before TimeFrame::initialise() is ever called, so the TimeFrame-owned
//    configuration (and any previously allocated LUT) is left byte-identical.
//  - A non-FirstPass iteration does not recommit or reallocate; its freshly
//    bound configuration must match the TimeFrame-owned one exactly, or the
//    call throws TraversalException{IndexTableConfigurationMismatch} before
//    TimeFrame is touched.
//  - A FirstPass iteration may legitimately change RowBins/ColBins/extents
//    relative to a previous iteration; that recommits and reallocates.
//  - Layout/binding are validated before TimeFrame state is committed.

#define BOOST_TEST_MODULE ITSMFT IndexTableActivation
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/IndexTableConfiguration.h"
#include "ITSMFTTracking/IndexTableUtils.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/detail/TimeFrameScratch.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "TraversalTestSupport.h"
#include "ITStracking/Constants.h"
#include "ITStracking/ROFLookupTables.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

// Deterministic, geometry-free stand-in for GeometryClusterDecoder<DetId>;
// identical construction to testTimeFrameLifecycle.cxx / testTrackerFailureContract.cxx.
class LegacyLikeDecoder final : public ClusterDecoder
{
 public:
  explicit LegacyLikeDecoder(o2::detectors::DetID::ID detector) : mDetector(detector) {}

  o2::itsmft::tracking::ClusterDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const TopologyDictionary* dict,
    uint32_t,
    bool applySysErrors) const override
  {
    const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dict);
    if (!clusterData.ok()) {
      o2::itsmft::tracking::ClusterDecodeResult result;
      result.error = clusterData.error;
      return result;
    }

    o2::itsmft::tracking::ClusterDecodeResult result;
    const int sensorID = cluster.getSensorID();
    auto& decoded = result.decoded;
    decoded.global = {static_cast<float>(sensorID) * 10.f, static_cast<float>(cluster.getRow()), static_cast<float>(cluster.getCol())};
    decoded.cylinderFrame = {static_cast<float>(sensorID) + 100.f, static_cast<float>(cluster.getRow()) + 1.f, static_cast<float>(cluster.getCol()) + 2.f, 0.01f * sensorID};
    decoded.rowColumnCovariance = {clusterData.sig2Row, 0.f, clusterData.sig2Col};
    decoded.shape = clusterData.shape;
    decoded.layer = sensorID;
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
    surfaces.push_back(SurfaceDescriptor{i, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
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

std::vector<LayerId> identitySurfaces(uint16_t nLayers)
{
  std::vector<LayerId> mapping;
  mapping.reserve(nLayers);
  for (uint16_t i = 0; i < nLayers; ++i) {
    mapping.push_back(LayerId{i});
  }
  return mapping;
}

struct Fixture {
  std::vector<CompClusterExt> clusters;
  std::vector<unsigned char> patterns;
  std::vector<ROFRecord> rofs;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
};

Fixture emptyFixture()
{
  return Fixture{};
}

// 4 clusters on layers {0,1,0,2}, partitioned into 3 ROFs -- same shape as
// testTrackerFailureContract.cxx's fixture.
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

std::vector<TrackingParameters> makeOneIterationITSParams()
{
  std::vector<TrackingParameters> params(1);
  resetDetectorDefaults(params[0], o2::detectors::DetID::ITS);
  return params;
}

// Two identical FirstPass ITS iterations, except PassFlags: iteration 1 drops
// FirstPass but keeps RebuildClusterLUT -- the shape TimeFrame::initialise()
// documents for legacy ITS's async iteration 3 (RebuildClusterLUT without
// FirstPass, reusing/resorting into whatever configuration the TimeFrame
// already owns). Callers mutate params[1] further to construct a mismatch.
std::vector<TrackingParameters> makeTwoIterationITSParams()
{
  std::vector<TrackingParameters> params(2);
  resetDetectorDefaults(params[0], o2::detectors::DetID::ITS);
  resetDetectorDefaults(params[1], o2::detectors::DetID::ITS);
  params[1].PassFlags = IterationSteps{IterationStep::RebuildClusterLUT};
  return params;
}

struct Rig {
  Rig() : pool(std::make_shared<BoundedMemoryResource>())
  {
  }

  std::shared_ptr<BoundedMemoryResource> pool;
  TimeFrame frame;
  TimeFrameScratch* tf{nullptr};
  Tracker tracker;
  TrackerTraits traits;
  // Scratch carries non-owning runtime ROF views. Keep these adapter-edge
  // builders alive across every load/initialise call on this fixture.
  std::optional<o2::its::ROFOverlapTable<ITSNLayers>> rofTable;
  std::optional<o2::its::ROFVertexLookupTable<ITSNLayers>> vertexTable;
  std::optional<o2::its::ROFMaskTable<ITSNLayers>> mask;
  std::vector<SurfaceDescriptor> catalog;

  void establishValidLayout(gsl::span<const TrackingParameters> params)
  {
    catalog = makeITSTestCatalog();
    const auto orderedSurfaces = identitySurfaces(ITSNLayers);
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    TrackerInitialization configuration;
    configuration.catalog = catalogView;
    configuration.memoryPool = pool;
    configuration.layout = makeDetectorLayout();
    configuration.parameters.assign(params.begin(), params.end());
    BOOST_REQUIRE(tracker.initialize(frame, configuration).ok());
    tf = &frame.getScratch();
  }

  // See testTrackerFailureContract.cxx's identical helper for why loading a
  // (possibly empty) normalized source is load-bearing before initialise()
  // is allowed to run at all: it sizes every per-layer ROF boundary table
  // that getNrof() unconditionally reads. The ROF overlap table and
  // multiplicity mask are equally load-bearing here: prepareClusters()
  // (reached from initialise() whenever RebuildClusterLUT is set, which is
  // every iteration exercised by this file) reads mROFMaskView::isROFEnabled()
  // unconditionally for every (layer, rof); an unset mask view is not merely
  // "everything enabled", it is a garbage/unconfigured view that crashes.
  void loadSource(const Fixture& f)
  {
    LegacyLikeDecoder decoder{o2::detectors::DetID::ITS};
    const o2::InteractionRecord origin{50, 5};
    const ROFTimingConfig timing{40, 0, 0, 0};
    const auto& layout = frame.getLayout();
    const auto layerMapping = identitySurfaces(ITSNLayers);
    const auto result = loadTimeFrameSource(frame, decoder, origin, timing, f.clusters, f.patterns, f.rofs, &dict(),
                                            f.labels.getIndexedSize() > 0 ? &f.labels : nullptr, o2::detectors::DetID::ITS,
                                            gsl::span<const LayerId>{layerMapping}, layout.getSurfaceCatalog());
    BOOST_REQUIRE(result.ok());

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
    frame.setROFViews(RuntimeROFViews{rofTable->getView(), vertexTable->getView(), mask->getView(), {}});
  }
};

std::vector<int> snapshotIndexTable(TimeFrame& tf, int layer)
{
  std::vector<int> flat;
  for (int rof = 0; rof < tf.getNrof(layer); ++rof) {
    const auto span = tf.getIndexTable(rof, layer);
    flat.insert(flat.end(), span.begin(), span.end());
  }
  return flat;
}

} // namespace

// --- FirstPass commit ------------------------------------------------------

BOOST_AUTO_TEST_CASE(FirstPassCommitsValidatedConfigurationIntoTimeFrame)
{
  Rig rig;
  auto params = makeOneIterationITSParams();
  rig.establishValidLayout(params);
  rig.loadSource(makeFixture());
  BOOST_TEST(rig.frame.getLayout().getSurfaceCatalog().getSurface(LayerId{0}).chartRange.min == -20.f);
  BOOST_TEST(rig.frame.getLayout().getSurfaceCatalog().getSurface(LayerId{0}).chartRange.max == 20.f);
  BOOST_CHECK_NO_THROW(TrackerTestAccess::prepare(rig.tracker, rig.frame, 0));

  IndexTableUtilsCore expected;
  std::array<SurfaceChartRange, ITSNLayers> chartRanges{};
  for (auto& range : chartRanges) {
    range = {-20.f, 20.f};
  }
  BOOST_REQUIRE(bindIndexTableConfiguration(expected, params[0], ITSNLayers, SurfaceKind::Cylinder, chartRanges) == IndexTableConfigError::None);
  BOOST_CHECK(indexTableConfigurationsMatch(rig.frame.getIndexTableUtils(), expected, ITSNLayers));
  BOOST_CHECK_GT(rig.frame.getIndexTableUtils().getNrowBins(), 0);
}

BOOST_AUTO_TEST_CASE(LoadingNormalizesGlobalsToBeamCoordinatesExactlyOnce)
{
  const auto params = makeOneIterationITSParams();

  Rig originReference;
  originReference.establishValidLayout(params);
  originReference.loadSource(makeFixture());

  Rig rig;
  rig.establishValidLayout(params);
  rig.frame.setBeamPosition(1.f, 2.f, 0.04f);
  rig.loadSource(makeFixture());

  const auto& referenceClusters = originReference.frame.getClusters();
  const auto& normalizedClusters = rig.frame.getClusters();
  BOOST_REQUIRE_EQUAL(normalizedClusters.size(), referenceClusters.size());
  for (size_t layer = 0; layer < normalizedClusters.size(); ++layer) {
    BOOST_REQUIRE_EQUAL(normalizedClusters[layer].size(), referenceClusters[layer].size());
    for (size_t cluster = 0; cluster < normalizedClusters[layer].size(); ++cluster) {
      const auto& reference = referenceClusters[layer][cluster];
      const auto& normalized = normalizedClusters[layer][cluster];
      BOOST_CHECK_EQUAL(normalized.clusterId, reference.clusterId);
      BOOST_CHECK_EQUAL(normalized.x, reference.x - 1.f);
      BOOST_CHECK_EQUAL(normalized.y, reference.y - 2.f);
      BOOST_CHECK_EQUAL(normalized.z, reference.z);
      BOOST_CHECK_CLOSE_FRACTION(normalized.radius, std::hypot(normalized.x, normalized.y), 1.e-6f);
      BOOST_CHECK_CLOSE_FRACTION(normalized.phi, o2::its::math_utils::computePhi(normalized.x, normalized.y), 1.e-6f);
    }
  }

  const auto beforePrepare = normalizedClusters;
  BOOST_REQUIRE_NO_THROW(TrackerTestAccess::prepare(rig.tracker, rig.frame, 0));

  const auto checkUnchanged = [&] {
    const auto& current = rig.frame.getClusters();
    BOOST_REQUIRE_EQUAL(current.size(), beforePrepare.size());
    for (size_t layer = 0; layer < current.size(); ++layer) {
      BOOST_REQUIRE_EQUAL(current[layer].size(), beforePrepare[layer].size());
      for (size_t cluster = 0; cluster < current[layer].size(); ++cluster) {
        BOOST_CHECK_EQUAL(current[layer][cluster].x, beforePrepare[layer][cluster].x);
        BOOST_CHECK_EQUAL(current[layer][cluster].y, beforePrepare[layer][cluster].y);
        BOOST_CHECK_EQUAL(current[layer][cluster].z, beforePrepare[layer][cluster].z);
      }
    }
  };
  checkUnchanged();
  BOOST_CHECK_EQUAL(rig.frame.getBeamPositionVariance(), 0.04f);

  BOOST_REQUIRE_NO_THROW(TrackerTestAccess::prepare(rig.tracker, rig.frame, 0));
  checkUnchanged();
}

// --- Invalid binding leaves TimeFrame untouched -----------------------------

BOOST_AUTO_TEST_CASE(InvalidBindingLeavesTimeFrameOwnedConfigurationUnchanged)
{
  Rig rig;
  auto params = makeOneIterationITSParams();
  params[0].RowBins = 0; // structurally invalid
  rig.catalog = makeITSTestCatalog();
  const auto orderedSurfaces = identitySurfaces(ITSNLayers);
  TrackerInitialization configuration;
  configuration.catalog = {rig.catalog.data(), static_cast<uint32_t>(rig.catalog.size())};
  configuration.memoryPool = rig.pool;
  configuration.layout = makeDetectorLayout();
  configuration.parameters = params;
  BOOST_CHECK(!rig.tracker.initialize(rig.frame, configuration).ok());
  BOOST_CHECK(!rig.frame.isConfigured());
}

// --- Non-FirstPass reuse -----------------------------------------------------

BOOST_AUTO_TEST_CASE(NonFirstPassMatchingReuseSucceedsWithoutRecommit)
{
  Rig rig;
  auto params = makeTwoIterationITSParams(); // params[1] is a byte-identical, non-FirstPass reuse
  rig.establishValidLayout(params);
  rig.loadSource(makeFixture());
  BOOST_REQUIRE_NO_THROW(TrackerTestAccess::prepare(rig.tracker, rig.frame, 0));
  const IndexTableUtilsCore committed = rig.frame.getIndexTableUtils();

  BOOST_CHECK_NO_THROW(TrackerTestAccess::prepare(rig.tracker, rig.frame, 1));
  BOOST_CHECK(indexTableConfigurationsMatch(rig.frame.getIndexTableUtils(), committed, ITSNLayers));
}

BOOST_AUTO_TEST_CASE(IterationSpecificRowColBinsAreRejectedAtInitialization)
{
  Rig rig;
  auto params = makeTwoIterationITSParams();
  params[1].RowBins = params[0].RowBins + 1;
  rig.catalog = makeITSTestCatalog();
  TrackerInitialization configuration;
  configuration.catalog = {rig.catalog.data(), static_cast<uint32_t>(rig.catalog.size())};
  configuration.memoryPool = rig.pool;
  configuration.layout = makeDetectorLayout();
  configuration.parameters = params;
  const auto result = rig.tracker.initialize(rig.frame, configuration);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK_EQUAL(result.failedIteration, 1u);
  BOOST_CHECK(!rig.frame.isConfigured());
}

BOOST_AUTO_TEST_CASE(ParameterSideExtentDoesNotOverrideDescriptorChartRange)
{
  Rig rig;
  auto params = makeTwoIterationITSParams();
  // The descriptor is now the authoritative chart-range source: changing a
  // retired parameter-side extent cannot alter a configured LUT.
  BOOST_REQUIRE(!params[1].LayerZ.empty());
  params[1].LayerZ[0] += 1.f;
  rig.establishValidLayout(params);
  rig.loadSource(makeFixture());
  BOOST_REQUIRE_NO_THROW(TrackerTestAccess::prepare(rig.tracker, rig.frame, 0));
  const IndexTableUtilsCore committedBefore = rig.frame.getIndexTableUtils();

  bool threw = false;
  try {
    TrackerTestAccess::prepare(rig.tracker, rig.frame, 1);
  } catch (const TraversalException& e) {
    threw = true;
    BOOST_CHECK(e.getReason() == TraversalFailureReason::IndexTableConfigurationMismatch);
  }
  BOOST_CHECK(!threw);
  BOOST_CHECK(indexTableConfigurationsMatch(rig.frame.getIndexTableUtils(), committedBefore, ITSNLayers));
}

BOOST_AUTO_TEST_CASE(FirstPassCannotChangeInvariantIndexTableConfiguration)
{
  Rig rig;
  auto params = makeOneIterationITSParams();
  params.push_back(params[0]);
  params[1].RowBins = params[0].RowBins + 8;
  rig.catalog = makeITSTestCatalog();
  TrackerInitialization configuration;
  configuration.catalog = {rig.catalog.data(), static_cast<uint32_t>(rig.catalog.size())};
  configuration.memoryPool = rig.pool;
  configuration.layout = makeDetectorLayout();
  configuration.parameters = params;
  const auto result = rig.tracker.initialize(rig.frame, configuration);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK_EQUAL(result.failedIteration, 1u);
  BOOST_CHECK(!rig.frame.isConfigured());
}

// Gate 4 B2 Slice 2 removed the MissingLayout reordering-regression test that
// used to live here: initialiseTimeFrame() now takes the plan as an explicit
// layout/topology view, so "no layout established"
// is no longer a state a caller can even construct -- there is no longer a
// TimeFrame-owned "missing layout" condition for TraversalFailureReason::
// MissingLayout to classify. The reordering property this test protected
// (TimeFrame::initialise() runs only after every structural check, never as
// initialiseTimeFrame()'s first statement) still holds structurally, since
// initialiseTimeFrame()'s only remaining precondition check on `layouts` is
// the IterationOutOfRange bounds check, still ahead of every other check.
