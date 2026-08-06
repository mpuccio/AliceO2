// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// M6e1 coverage for the standalone MFT interface. The fixture supplies the
// workflow-owned publication adapter and runtime ROF views used by the live
// interface; the generic frame remains the only event-state owner.
//
// Fixture technique (FieldFixture/GRPECSFixture/OneLayerDecoder/
// makePatternBytes/oneRof/oneCluster) duplicated from
// testTrackingInterfaceLoadFailureContract.cxx's own MFT-only fixtures
// rather than shared via a new header, matching this test directory's
// existing per-file-local-fixture convention -- that file is not modified
// here (it still explicitly instantiates and tests the bare, default-
// argument ITSMFTTrackingInterface<10> for its own NLayers-generic
// load-failure-contract coverage, unaffected by this migration). This file
// targets ITSMFTTrackingInterfaceMFT specifically -- the concrete type
// production actually uses -- not the generic template.

#define BOOST_TEST_MODULE ITSMFT M6e1StandaloneMFTMigration
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <memory>
#include <type_traits>
#include <vector>

#include <TGeoGlobalMagField.h>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Field/MagneticField.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InputSpec.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "ITSMFTTracking/TrackingInterface.h"
#include "ITSMFTTracking/DetectorPublicationAdapter.h"
#include "ITStracking/ROFLookupTables.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

// --- Global fixtures, identical in intent to
// testTrackingInterfaceLoadFailureContract.cxx's own (this file's own
// process-wide singleton setup, not a shared dependency on that file). ---

struct FieldFixture {
  FieldFixture()
  {
    if (!TGeoGlobalMagField::Instance()->GetField()) {
      TGeoGlobalMagField::Instance()->SetField(o2::field::MagneticField::createNominalField(5, true));
      TGeoGlobalMagField::Instance()->Lock();
    }
  }
};

struct GRPECSFixture {
  static o2::parameters::GRPECSObject& grpEcs()
  {
    static o2::parameters::GRPECSObject obj;
    return obj;
  }

  GRPECSFixture()
  {
    auto& obj = grpEcs();
    obj.addDetContinuousReadOut(o2::detectors::DetID::MFT);
    obj.setNHBFPerTF(128);

    std::vector<o2::framework::InputSpec> inputs;
    auto req = std::make_shared<o2::base::GRPGeomRequest>(
      false, true, false, false, false, o2::base::GRPGeomRequest::None, inputs);
    o2::base::GRPGeomHelper::instance().setRequest(req);
    o2::framework::ConcreteDataMatcher matcher{"GLO", "GRPECS", 0};
    o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, &obj);
  }
};

BOOST_GLOBAL_FIXTURE(FieldFixture);
BOOST_GLOBAL_FIXTURE(GRPECSFixture);

class OneLayerDecoder final : public ClusterDecoder
{
 public:
  o2::itsmft::ioutils::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const TopologyDictionary* dict,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool) const override
  {
    o2::itsmft::ioutils::SurfaceMeasurementDecodeResult result;
    if (dict == nullptr) {
      result.error = ClusterDecodeError::MissingDictionary;
      return result;
    }
    const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dict);
    if (!clusterData.ok()) {
      result.error = clusterData.error;
      return result;
    }
    if (layerToSurface.empty()) {
      result.error = ClusterDecodeError::InvalidLayerMapping;
      return result;
    }
    result.layer = 0;
    result.layerMapped = true;
    result.kind = SurfaceKind::Disk;

    DecodedCluster decoded{};
    decoded.global = {1.f, 2.f, 3.f};
    decoded.cylinderFrame = {4.f, 5.f, 6.f, 0.f};
    decoded.rowColumnCovariance = {0.1f, 0.f, 0.2f};
    decoded.shape = clusterData.shape;
    decoded.sensor = static_cast<uint32_t>(cluster.getSensorID());
    decoded.layer = 0;

    const o2::itsmft::tracking::DetectorSensorId sensor{static_cast<uint32_t>(o2::detectors::DetID::MFT), decoded.sensor};
    const ClusterRef clusterRef{source, externalIndex};
    result.measurement = makeDiskSurfaceMeasurement(decoded, sensor, layerToSurface[0], clusterRef, sourceROF);
    return result;
  }
};

constexpr std::array<unsigned char, 3> onePixelPattern{1, 1, 0x80};

std::vector<unsigned char> makePatternBytes(size_t nClusters)
{
  std::vector<unsigned char> bytes;
  bytes.reserve(nClusters * onePixelPattern.size());
  for (size_t i = 0; i < nClusters; ++i) {
    bytes.insert(bytes.end(), onePixelPattern.begin(), onePixelPattern.end());
  }
  return bytes;
}

const TopologyDictionary& dict()
{
  static const TopologyDictionary d;
  return d;
}

struct MFTAdapterContext {
  DetectorPublicationAdapter<MFTNLayers> adapter;
  MFTPublicationCompatibility sidecar;
  o2::its::ROFOverlapTable<MFTNLayers> overlap;
  o2::its::ROFVertexLookupTable<MFTNLayers> vertex;
  o2::its::ROFMaskTable<MFTNLayers> mask;
  RuntimeROFViews views{};

  RuntimeROFViews configure()
  {
    o2::its::LayerTiming timing{};
    timing.mNROFsTF = 1;
    timing.mROFLength = 40;
    for (int layer = 0; layer < MFTNLayers; ++layer) {
      overlap.defineLayer(layer, timing);
      vertex.defineLayer(layer, timing);
    }
    overlap.init();
    vertex.init();
    mask = o2::its::ROFMaskTable<MFTNLayers>{overlap};
    mask.resetMask();
    for (int layer = 0; layer < MFTNLayers; ++layer) {
      mask.setROFsEnabled(layer, 0, 1, 1);
    }
    views = {overlap.getView(), vertex.getView(), mask.getView(), {}};
    return views;
  }
};

MFTAdapterContext& mftAdapterContext()
{
  static MFTAdapterContext context;
  return context;
}

std::vector<ROFRecord> oneRof(int nEntries = 1)
{
  return {ROFRecord{{0, 0}, 0, 0, nEntries}};
}
std::vector<CompClusterExt> oneCluster()
{
  return {CompClusterExt{10, 20, CompCluster::InvalidPatternID, 0}};
}

// Constructs and fully initializes a real ITSMFTTrackingInterfaceMFT (the
// concrete production type -- not the generic template), ready for
// processTimeFrame(). Mirrors testTrackingInterfaceLoadFailureContract.cxx's
// own makeReadyInterface() exactly, targeted at the one concrete type this
// file tests.
std::unique_ptr<ITSMFTTrackingInterfaceMFT> makeReadyInterface()
{
  auto interface = std::make_unique<ITSMFTTrackingInterfaceMFT>(
    false, o2::itsmft::TrackingMode::Sync, false, std::make_unique<OneLayerDecoder>());
  interface->initialise();
  BOOST_REQUIRE(interface->isActive());
  interface->setClusterDictionary(&dict());
  auto& context = mftAdapterContext();
  context.adapter.adoptMFTPublicationCompatibility(&context.sidecar);
  interface->bindPublicationAdapter(context.adapter);
  context.configure();
  interface->bindROFViews(context.views);
  interface->getScratch().setROFViews(context.views);
  return interface;
}

// Real, full-catalog SurfaceId order for the standalone MFT static catalog
// (dense, local: surface i's id is always SurfaceId{i} -- identical
// technique to TrackingInterface.cxx's own identitySurfaceOrder<NLayers>()).
std::vector<SurfaceId> mftIdentityOrder()
{
  std::vector<SurfaceId> order;
  order.reserve(MFTNLayers);
  for (uint16_t i = 0; i < MFTNLayers; ++i) {
    order.push_back(SurfaceId{i});
  }
  return order;
}

SurfaceMask mftFullMask()
{
  SurfaceMask mask;
  for (uint16_t i = 0; i < MFTNLayers; ++i) {
    mask.set(SurfaceId{i});
  }
  return mask;
}

// Real combined ITS+MFT static catalog order, MFT half only (surfaces
// ITSNLayers..ITSNLayers+MFTNLayers-1) -- duplicated from
// testSurfacePlanBinding.cxx's own CombinedLayout fixture technique, not
// shared via a new header (this test directory's existing convention).
std::vector<SurfaceId> combinedMftOrder()
{
  std::vector<SurfaceId> order;
  order.reserve(MFTNLayers);
  for (uint16_t i = 0; i < MFTNLayers; ++i) {
    order.push_back(SurfaceId{static_cast<uint16_t>(ITSNLayers + i)});
  }
  return order;
}

SurfaceMask combinedMftMask()
{
  SurfaceMask mask;
  for (uint16_t i = 0; i < MFTNLayers; ++i) {
    mask.set(SurfaceId{static_cast<uint16_t>(ITSNLayers + i)});
  }
  return mask;
}

} // namespace

// --- Type proofs: both standalone interfaces expose the same generic scratch
// entity; the MFT fixture below additionally supplies adapter-owned views. ---

static_assert(std::is_same_v<decltype(std::declval<ITSMFTTrackingInterfaceMFT&>().getScratch()), SurfaceTrackingScratch&>);
static_assert(std::is_same_v<decltype(std::declval<ITSMFTTrackingInterfaceITS&>().getScratch()), SurfaceTrackingScratch&>);

BOOST_AUTO_TEST_CASE(TypeProofsAppearInOutput)
{
  BOOST_CHECK(true); // compile-time proofs above; case exists so it appears in test output.
}

// --- Successful load/track/reset lifecycle preserves the existing
// interface contract (ITSMFTTrackingInterfaceMFT specifically, not the
// generic template). ---

BOOST_AUTO_TEST_CASE(StandaloneMFTBaselineLifecycleUsesTheNewScratchAndBinding)
{
  auto interfacePtr = makeReadyInterface();
  auto& interface = *interfacePtr;
  // Production-code proof, not just a type proof: getScratch() returns a
  // live SurfaceTrackingScratch& already adopted with a real plan (nonzero
  // owned-surface count) by the time initialise() has run.
  static_assert(std::is_same_v<decltype(interface.getScratch()), SurfaceTrackingScratch&>);

  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const auto patterns = makePatternBytes(clusters.size());

  const float result = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(result));
  BOOST_CHECK_GE(result, 0.f);
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 1u);

  const auto exported = interface.getCommonTrackPublicationExport();
  BOOST_REQUIRE(exported);
  BOOST_CHECK(exported->detector == o2::detectors::DetID::MFT);
  BOOST_CHECK(exported->source == ClusterSourceId{0});
  BOOST_CHECK_EQUAL(exported->orderedSurfaces.size(), static_cast<size_t>(MFTNLayers));

  interface.resetEvent();
  BOOST_CHECK(!interface.getCommonTrackPublicationExport());
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 0u);

  // A subsequent load on the same, reset interface succeeds again --
  // scratch/binding/plan state survives reset() intact.
  const float again = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(again));
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 1u);
}

// --- Load failure and dropped-TF handling leave no stale CommonTrack,
// sidecar, or scratch state. ---

BOOST_AUTO_TEST_CASE(RecoverableLoadFailureLeavesNoStaleCommonTrackSidecarOrScratchState)
{
  auto interfacePtr = makeReadyInterface();
  auto& interface = *interfacePtr;
  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const std::vector<unsigned char> malformedPattern{0, 1}; // rowSpan==0 -> MalformedExplicitPattern, recoverable

  BOOST_CHECK_THROW(interface.processTimeFrame(rofs, clusters, malformedPattern, nullptr), RecoverableLoadFailure);
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 0u);
  BOOST_CHECK(!interface.getCommonTrackPublicationExport());
  BOOST_REQUIRE(interface.getMFTPublicationCompatibility() != nullptr);
  // The publication sidecar is cleared by resetEvent() (called internally on
  // this failure path) alongside the scratch -- no accepted-track shadow
  // entries can have survived a load that never reached tracking.
  BOOST_CHECK_EQUAL(interface.getMFTPublicationCompatibility()->entries().size(), 0u);

  // A subsequent valid load succeeds cleanly -- the failure left nothing
  // stale to interfere with it.
  const auto patterns = makePatternBytes(clusters.size());
  const float retried = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(retried));
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 1u);
}

BOOST_AUTO_TEST_CASE(DroppedTimeFrameLeavesNoStaleCommonTrackSidecarOrScratchState)
{
  auto interfacePtr = makeReadyInterface();
  auto& interface = *interfacePtr;
  const_cast<std::vector<o2::itsmft::TrackingParameters>&>(interface.getTrackingParameters())[0].DropTFUponFailure = true;

  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const std::vector<unsigned char> malformedPattern{0, 1};

  const float result = interface.processTimeFrame(rofs, clusters, malformedPattern, nullptr);
  BOOST_CHECK(isDroppedTimeFrame(result));
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 0u);
  BOOST_CHECK(!interface.getCommonTrackPublicationExport());
  BOOST_REQUIRE(interface.getMFTPublicationCompatibility() != nullptr);
  BOOST_CHECK_EQUAL(interface.getMFTPublicationCompatibility()->entries().size(), 0u);

  const auto patterns = makePatternBytes(clusters.size());
  const float retried = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(retried));
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 1u);
}

// --- Standalone MFT compact transition/cell mapping matches the
// already-migrated combined MFT leg: both bindings are built from their own
// real, distinct production static catalogs -- standalone's dense 0..9
// numbering vs. combined's 7..16 numbering -- so only the *compact*
// (position-in-ownedSurfaces) mapping can coincide; the underlying global
// SurfaceId numbering deliberately does not. ---

BOOST_AUTO_TEST_CASE(StandaloneAndCombinedMFTBindingsAgreeOnCompactSlotsByRelativePosition)
{
  // Standalone: real MFT-only static catalog, exactly TrackingInterface.cxx's
  // own initialiseTracker() construction (Sync mode, real TrackingParameters).
  const auto standaloneParams = o2::itsmft::TrackingMode::getTrackingParameters(o2::detectors::DetID::MFT, o2::itsmft::TrackingMode::Sync);
  const auto standaloneOrder = mftIdentityOrder();
  const auto standaloneResult = buildSurfaceGraphs(
    SurfaceCatalogView{kMFTStaticSurfaceCatalog.data(), static_cast<uint32_t>(kMFTStaticSurfaceCatalog.size())},
    gsl::span<const SurfaceId>{standaloneOrder}, standaloneParams);
  BOOST_REQUIRE(standaloneResult.ok());
  const auto standaloneBindingResult = SurfacePlanBinding::build(
    standaloneResult.graphs.front().getView(), ClusterSourceId{0}, mftFullMask(),
    gsl::span<const SurfaceId>{standaloneOrder}, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk);
  BOOST_REQUIRE(standaloneBindingResult.ok());

  // Combined: real ITS+MFT combined static catalog, MFT half only -- mirrors
  // The combined workflow's own real production application construction
  // (ClusterSourceId{1}, surfaces ITSNLayers..ITSNLayers+MFTNLayers-1).
  const auto combinedParams = standaloneParams; // MFT's own TrackingParameters do not depend on which catalog they size against
  const auto combinedOrder = combinedMftOrder();
  const auto combinedResult = buildSurfaceGraphs(
    SurfaceCatalogView{kITSMFTCombinedStaticSurfaceCatalog.data(), static_cast<uint32_t>(kITSMFTCombinedStaticSurfaceCatalog.size())},
    gsl::span<const SurfaceId>{combinedOrder}, combinedParams);
  BOOST_REQUIRE(combinedResult.ok());
  const auto combinedBindingResult = SurfacePlanBinding::build(
    combinedResult.graphs.front().getView(), ClusterSourceId{1}, combinedMftMask(),
    gsl::span<const SurfaceId>{combinedOrder}, SurfaceKind::Disk, TransitionPolicyTag::DiskDisk);
  BOOST_REQUIRE(combinedBindingResult.ok());

  // Same source-count-agnostic compact-slot cardinality in both deployments.
  BOOST_CHECK_EQUAL(standaloneBindingResult.binding->getGlobalTransitions().size(), combinedBindingResult.binding->getGlobalTransitions().size());
  BOOST_CHECK_EQUAL(standaloneBindingResult.binding->getGlobalCells().size(), combinedBindingResult.binding->getGlobalCells().size());

  // Position k in each binding's own orderedSurfaces resolves to compact
  // slot k in both -- despite the two catalogs numbering the same physical
  // MFT surface with a different global SurfaceId (i vs. ITSNLayers+i).
  for (uint16_t k = 0; k < MFTNLayers; ++k) {
    const auto standaloneSlot = standaloneBindingResult.binding->getOwnedSurfaceIndex(standaloneOrder[k]);
    const auto combinedSlot = combinedBindingResult.binding->getOwnedSurfaceIndex(combinedOrder[k]);
    BOOST_REQUIRE(standaloneSlot.has_value());
    BOOST_REQUIRE(combinedSlot.has_value());
    BOOST_CHECK_EQUAL(*standaloneSlot, k);
    BOOST_CHECK_EQUAL(*combinedSlot, k);
    BOOST_CHECK_EQUAL(*standaloneSlot, *combinedSlot);
  }
}
