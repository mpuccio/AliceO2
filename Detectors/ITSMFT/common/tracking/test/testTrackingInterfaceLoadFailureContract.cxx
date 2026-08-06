// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Interface-level loading failure-contract tests. No real ITS/MFT geometry
// singleton is required: a ClusterDecoder and workflow-owned runtime ROF
// context are injected, while GRPECS is supplied through the same CCDB hook
// used by a DPL device.
//
// The fixture is exercised with the live MFT interface instantiation. Generic
// loader and frame failure contracts are covered independently by the
// NLayers-neutral TimeFrame tests.

#define BOOST_TEST_MODULE ITSMFT TrackingInterfaceLoadFailureContract
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <memory>
#include <stdexcept>
#include <vector>

#include <TGeoGlobalMagField.h>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InputSpec.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "ITSMFTTracking/TimeFrameLoadFailure.h"
#include "ITSMFTTracking/TrackingInterface.h"
#include "ITStracking/Constants.h"
#include "ITStracking/ROFLookupTables.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

// --- Global fixtures: field + a minimal GRPECSObject, both process-wide
// singletons that must be configured exactly once for this test binary. ---

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
  // Must outlive every test case: GRPGeomHelper only stores the raw pointer
  // finaliseCCDB() was given, exactly like a real CCDB-fetched object would.
  static o2::parameters::GRPECSObject& grpEcs()
  {
    static o2::parameters::GRPECSObject obj;
    return obj;
  }

  GRPECSFixture()
  {
    auto& obj = grpEcs();
    obj.addDetContinuousReadOut(o2::detectors::DetID::MFT); // MFT runtime timing context uses this condition.
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

// --- Fake decoder: geometry-free, deterministic single-layer mapping
// (sensorID -> layer 0), reusing the real bounded pattern-consumption path
// (extractClusterDataBounded) so malformed/truncated pattern bytes exercise
// exactly the same typed ClusterDecodeError classification production does.
// A separate ThrowingDecoder below covers the "decoder throws" case. ---

class OneLayerDecoder final : public ClusterDecoder
{
 public:
  OneLayerDecoder(o2::detectors::DetID::ID detector, SurfaceKind kind) : mDetector(detector), mKind(kind) {}

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
    result.kind = mKind;

    DecodedCluster decoded{};
    decoded.global = {1.f, 2.f, 3.f};
    decoded.cylinderFrame = {4.f, 5.f, 6.f, 0.f};
    decoded.rowColumnCovariance = {0.1f, 0.f, 0.2f};
    decoded.shape = clusterData.shape;
    decoded.sensor = static_cast<uint32_t>(cluster.getSensorID());
    decoded.layer = 0;

    const o2::itsmft::tracking::DetectorSensorId sensor{static_cast<uint32_t>(mDetector), decoded.sensor};
    const ClusterRef clusterRef{source, externalIndex};
    result.measurement = (mKind == SurfaceKind::Disk)
                           ? makeDiskSurfaceMeasurement(decoded, sensor, layerToSurface[0], clusterRef, sourceROF)
                           : makeCylinderSurfaceMeasurement(decoded, sensor, layerToSurface[0], clusterRef, sourceROF);
    return result;
  }

 private:
  o2::detectors::DetID::ID mDetector;
  SurfaceKind mKind;
};

class ThrowingDecoder final : public ClusterDecoder
{
 public:
  o2::itsmft::ioutils::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt&, BoundedPatternCursor&, const TopologyDictionary*,
    gsl::span<const SurfaceId>, ClusterSourceId, uint32_t, uint32_t, bool) const override
  {
    throw std::runtime_error{"ThrowingDecoder: simulated decoder bug"};
  }
};

// One explicit (non-grouped) 1-pixel pattern: rowSpan=1, colSpan=1, one
// bitmap byte -- three bytes consumed per cluster (see testMultiSourceLoading.cxx).
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

template <int NLayers>
struct TestTraits;
template <>
struct TestTraits<7> {
  static constexpr auto detId = o2::detectors::DetID::ITS;
  static constexpr auto kind = SurfaceKind::Cylinder;
  static constexpr auto policyTag = TransitionPolicyTag::CylinderCylinder;
};
template <>
struct TestTraits<10> {
  static constexpr auto detId = o2::detectors::DetID::MFT;
  static constexpr auto kind = SurfaceKind::Disk;
  static constexpr auto policyTag = TransitionPolicyTag::DiskDisk;
};

// Constructs and fully initializes an interface (real Sync-mode tracking
// parameters, real static-catalog plan construction -- Gate 4 B2 Slice 2:
// ITSMFTTrackingInterface::initialise() builds its one immutable plan
// unconditionally from the compile-time-selected static per-detector
// catalog, so no separate layout-configuration call is needed or possible
// any more) ready for processTimeFrame(). decoderOut is a non-owning
// observer into the interface's injected decoder, passed via the 4-arg
// constructor overload (Gate 4 B2 Slice 3 removed the catalogProvider
// parameter the 5-arg overload used to take).
// Returned via unique_ptr, not by value: ITSMFTTrackingInterface has
// unique_ptr members (mTrackerTraits, mTracker, ...), so its copy
// constructor is deleted and its implicit move constructor is not
// guaranteed to be usable either (SurfaceTrackingScratch is not
// trivially movable); returning it by value would require relying on unguaranteed
// NRVO through several mutating calls in between.
template <int NLayers>
struct TestAdapterContext {
  DetectorPublicationAdapter<NLayers> adapter;
  ITSSharedClusterCompatibility itsSidecar;
  MFTPublicationCompatibility mftSidecar;
  o2::its::ROFOverlapTable<NLayers> overlap;
  o2::its::ROFVertexLookupTable<NLayers> vertex;
  o2::its::ROFMaskTable<NLayers> mask;
  RuntimeROFViews views{};

  RuntimeROFViews configure()
  {
    o2::its::LayerTiming timing{};
    timing.mNROFsTF = 1;
    timing.mROFLength = 40;
    for (int layer = 0; layer < NLayers; ++layer) {
      overlap.defineLayer(layer, timing);
      vertex.defineLayer(layer, timing);
    }
    overlap.init();
    vertex.init();
    mask = o2::its::ROFMaskTable<NLayers>{overlap};
    mask.resetMask();
    for (int layer = 0; layer < NLayers; ++layer) {
      mask.setROFsEnabled(layer, 0, 1, 1);
    }
    views = {overlap.getView(), vertex.getView(), mask.getView(), {}};
    return views;
  }
};

template <int NLayers>
TestAdapterContext<NLayers>& testAdapterContext()
{
  static TestAdapterContext<NLayers> context;
  return context;
}

template <int NLayers, typename DecoderT = OneLayerDecoder>
std::unique_ptr<ITSMFTTrackingInterface<NLayers>> makeReadyInterface(DecoderT*& decoderOut)
{
  using Traits = TestTraits<NLayers>;
  std::unique_ptr<DecoderT> decoder;
  if constexpr (std::is_same_v<DecoderT, OneLayerDecoder>) {
    decoder = std::make_unique<OneLayerDecoder>(Traits::detId, Traits::kind);
  } else {
    decoder = std::make_unique<DecoderT>();
  }
  decoderOut = decoder.get();

  auto interface = std::make_unique<ITSMFTTrackingInterface<NLayers>>(
    false, o2::itsmft::TrackingMode::Sync, false, std::move(decoder));
  interface->initialise();
  BOOST_REQUIRE(interface->isActive());
  interface->setClusterDictionary(&dict());
  auto& context = testAdapterContext<NLayers>();
  if constexpr (NLayers == ITSNLayers) {
    context.adapter.adoptITSSharedClusterCompatibility(&context.itsSidecar);
  } else {
    context.adapter.adoptMFTPublicationCompatibility(&context.mftSidecar);
  }
  interface->bindPublicationAdapter(context.adapter);
  context.configure();
  interface->bindROFViews(context.views);
  interface->getScratch().setROFViews(context.views);
  return interface;
}

// One ROF, one cluster on sensor 0, everything the OneLayerDecoder needs.
std::vector<ROFRecord> oneRof(int nEntries = 1)
{
  return {ROFRecord{{0, 0}, 0, 0, nEntries}};
}
std::vector<CompClusterExt> oneCluster()
{
  return {CompClusterExt{10, 20, CompCluster::InvalidPatternID, 0}};
}

// RAII guard for the process-wide DPLAlpideParam<MFT> singleton: stashes one
// layer's roFrameLayerLengthInBC on construction and restores it on
// destruction (including on a thrown/failed check), so a test that stages a
// non-uniform per-layer configuration can never leak that state into a
// later test case in this same test binary.
o2::itsmft::DPLAlpideParam<o2::detectors::DetID::MFT>& mutableMFTAlpideParam()
{
  // ConfigurableParamHelper<P>::Instance() deliberately returns const P& to
  // discourage casual mutation from production code; tests that need to
  // stage a specific singleton value use the same const_cast pattern as
  // e.g. testCommonTrackingParameters.cxx (o2::its::TrackerParamConfig) and
  // several ZDC RecoParamZDC call sites.
  return const_cast<o2::itsmft::DPLAlpideParam<o2::detectors::DetID::MFT>&>(
    o2::itsmft::DPLAlpideParam<o2::detectors::DetID::MFT>::Instance());
}

class ScopedMFTLayerROFLengthOverride
{
 public:
  ScopedMFTLayerROFLengthOverride(int layer, int overrideValue)
    : mLayer(layer), mOriginal(mutableMFTAlpideParam().roFrameLayerLengthInBC[layer])
  {
    mutableMFTAlpideParam().roFrameLayerLengthInBC[mLayer] = overrideValue;
  }
  ~ScopedMFTLayerROFLengthOverride()
  {
    mutableMFTAlpideParam().roFrameLayerLengthInBC[mLayer] = mOriginal;
  }
  ScopedMFTLayerROFLengthOverride(const ScopedMFTLayerROFLengthOverride&) = delete;
  ScopedMFTLayerROFLengthOverride& operator=(const ScopedMFTLayerROFLengthOverride&) = delete;

 private:
  int mLayer;
  int mOriginal;
};

} // namespace

// ---------------------------------------------------------------------
// Baseline success
// ---------------------------------------------------------------------

template <int NLayers>
void checkBaselineSuccess()
{
  OneLayerDecoder* decoder = nullptr;
  auto interfacePtr = makeReadyInterface<NLayers>(decoder);
  auto& interface = *interfacePtr;
  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const auto patterns = makePatternBytes(clusters.size());

  const float result = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(result));
  BOOST_CHECK_GE(result, 0.f);
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 1u);
}

BOOST_AUTO_TEST_CASE(MFT_BaselineSuccess) { checkBaselineSuccess<10>(); }

template <int NLayers>
void checkCommonTrackPublicationExportLifecycle()
{
  OneLayerDecoder* decoder = nullptr;
  auto interface = makeReadyInterface<NLayers>(decoder);
  BOOST_CHECK(!interface->getCommonTrackPublicationExport());
  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const auto patterns = makePatternBytes(clusters.size());
  BOOST_REQUIRE_GE(interface->processTimeFrame(rofs, clusters, patterns, nullptr), 0.f);
  const auto exported = interface->getCommonTrackPublicationExport();
  BOOST_REQUIRE(exported);
  BOOST_CHECK(exported->detector == TestTraits<NLayers>::detId);
  BOOST_CHECK(exported->source == ClusterSourceId{0});
  BOOST_CHECK_EQUAL(exported->orderedSurfaces.size(), static_cast<size_t>(NLayers));
  const auto direct = exported->clock.getLegacyClockLayer();
  const CommonTrackTimestamp timestamp{0, 1};
  BOOST_REQUIRE(exported->clock.makeOutputTimestamp(timestamp));
  BOOST_CHECK_EQUAL(exported->clock.getROF(*exported->clock.makeOutputTimestamp(timestamp)), direct.getROF(*exported->clock.makeOutputTimestamp(timestamp)));
  interface->resetEvent();
  BOOST_CHECK(!interface->getCommonTrackPublicationExport());
}

BOOST_AUTO_TEST_CASE(MFT_CommonTrackPublicationExportLifecycle) { checkCommonTrackPublicationExportLifecycle<10>(); }

// ---------------------------------------------------------------------
// Recoverable: malformed pattern, gated by DropTFUponFailure
// ---------------------------------------------------------------------

template <int NLayers>
void checkRecoverableMalformedPatternIsGated()
{
  OneLayerDecoder* decoder = nullptr;
  auto interfacePtr = makeReadyInterface<NLayers>(decoder);
  auto& interface = *interfacePtr;
  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const std::vector<unsigned char> malformedPattern{0, 1}; // rowSpan==0 -> MalformedExplicitPattern

  BOOST_CHECK_THROW(interface.processTimeFrame(rofs, clusters, malformedPattern, nullptr), RecoverableLoadFailure);
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 0u);

  // A subsequent valid load on the same interface succeeds (state was
  // preserved/restored, catalog/layout/detId untouched by the failure).
  const auto patterns = makePatternBytes(clusters.size());
  const float retried = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(retried));
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 1u);
}

BOOST_AUTO_TEST_CASE(MFT_RecoverableMalformedPatternThrowsWhenNotDropped) { checkRecoverableMalformedPatternIsGated<10>(); }

template <int NLayers>
void checkRecoverableFailureIsDroppedWhenConfigured()
{
  OneLayerDecoder* decoder = nullptr;
  auto interfacePtr = makeReadyInterface<NLayers>(decoder);
  auto& interface = *interfacePtr;
  const_cast<std::vector<o2::itsmft::TrackingParameters>&>(interface.getTrackingParameters())[0].DropTFUponFailure = true;

  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const std::vector<unsigned char> malformedPattern{0, 1};

  const float result = interface.processTimeFrame(rofs, clusters, malformedPattern, nullptr);
  BOOST_CHECK(isDroppedTimeFrame(result));
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 0u);

  const auto patterns = makePatternBytes(clusters.size());
  const float retried = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(retried));
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 1u);
}

BOOST_AUTO_TEST_CASE(MFT_RecoverableFailureDroppedWhenConfigured) { checkRecoverableFailureIsDroppedWhenConfigured<10>(); }

// ---------------------------------------------------------------------
// Recoverable: invalid ROF range, never confused with a structural failure
// ---------------------------------------------------------------------

template <int NLayers>
void checkInvalidROFRangeIsRecoverable()
{
  OneLayerDecoder* decoder = nullptr;
  auto interfacePtr = makeReadyInterface<NLayers>(decoder);
  auto& interface = *interfacePtr;
  // Two clusters declared, but the single ROF only claims one: an internal
  // gap the transactional loader must reject before any decode.
  const std::vector<CompClusterExt> clusters{
    CompClusterExt{10, 20, CompCluster::InvalidPatternID, 0},
    CompClusterExt{11, 21, CompCluster::InvalidPatternID, 0}};
  const auto rofs = oneRof(1); // claims only 1 of the 2 clusters
  const auto patterns = makePatternBytes(clusters.size());

  BOOST_CHECK_THROW(interface.processTimeFrame(rofs, clusters, patterns, nullptr), RecoverableLoadFailure);
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 0u);
}

BOOST_AUTO_TEST_CASE(MFT_InvalidROFRangeIsRecoverable) { checkInvalidROFRangeIsRecoverable<10>(); }

// ---------------------------------------------------------------------
// Structural: dictionary never configured -- distinct from a malformed
// pattern/dictionary-content failure on an otherwise-configured source.
// ---------------------------------------------------------------------

template <int NLayers>
void checkDictionaryNotConfiguredIsStructural()
{
  OneLayerDecoder* decoder = nullptr;
  using Traits = TestTraits<NLayers>;
  auto decoderOwner = std::make_unique<OneLayerDecoder>(Traits::detId, Traits::kind);
  decoder = decoderOwner.get();
  ITSMFTTrackingInterface<NLayers> interface{false, o2::itsmft::TrackingMode::Sync, false, std::move(decoderOwner)};
  interface.initialise();
  auto& context = testAdapterContext<NLayers>();
  if constexpr (NLayers == ITSNLayers) {
    context.adapter.adoptITSSharedClusterCompatibility(&context.itsSidecar);
  } else {
    context.adapter.adoptMFTPublicationCompatibility(&context.mftSidecar);
  }
  interface.bindPublicationAdapter(context.adapter);
  context.configure();
  interface.bindROFViews(context.views);
  interface.getScratch().setROFViews(context.views);
  // setClusterDictionary() deliberately never called.

  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const auto patterns = makePatternBytes(clusters.size());

  // Structural: never gated by DropTFUponFailure, even when set.
  const_cast<std::vector<o2::itsmft::TrackingParameters>&>(interface.getTrackingParameters())[0].DropTFUponFailure = true;
  BOOST_CHECK_EXCEPTION(interface.processTimeFrame(rofs, clusters, patterns, nullptr), TimeFrameLoadException,
                        [](const TimeFrameLoadException& e) { return e.reason() == TimeFrameLoadFailureReason::DictionaryNotConfigured; });

  interface.setClusterDictionary(&dict());
  const float retried = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(retried));
}

BOOST_AUTO_TEST_CASE(MFT_DictionaryNotConfiguredIsStructural) { checkDictionaryNotConfiguredIsStructural<10>(); }

// ---------------------------------------------------------------------
// Distinguish DictionaryNotConfigured from a malformed-content failure on
// an otherwise-configured (non-null) dictionary.
// ---------------------------------------------------------------------

template <int NLayers>
void checkMalformedPatternWithConfiguredDictionaryIsRecoverableNotMissingDictionary()
{
  OneLayerDecoder* decoder = nullptr;
  auto interfacePtr = makeReadyInterface<NLayers>(decoder); // setClusterDictionary() already called with a non-null dictionary
  auto& interface = *interfacePtr;
  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const std::vector<unsigned char> malformedPattern{0, 1};

  bool sawRecoverable = false;
  try {
    interface.processTimeFrame(rofs, clusters, malformedPattern, nullptr);
    BOOST_FAIL("expected RecoverableLoadFailure");
  } catch (const RecoverableLoadFailure& err) {
    sawRecoverable = true;
    BOOST_CHECK(err.error() != MultiSourceLoadError::MissingDictionary);
    BOOST_CHECK(err.error() == MultiSourceLoadError::MalformedExplicitPattern);
  }
  BOOST_CHECK(sawRecoverable);
}

BOOST_AUTO_TEST_CASE(MFT_MalformedPatternWithConfiguredDictionaryIsRecoverable) { checkMalformedPatternWithConfiguredDictionaryIsRecoverableNotMissingDictionary<10>(); }

// Gate 4 B2 Slice 2 removed the checkUnconfiguredCatalogIsStructural test
// that used to live here (skip configureSurfaceGraphs() to produce
// MultiSourceLoadError::SurfaceCatalogNotConfigured): configureSurfaceGraphs()
// no longer exists, and ITSMFTTrackingInterface::initialise() now
// unconditionally builds its one immutable plan from the static per-detector
// catalog before returning -- there is no longer a way to construct a ready,
// active interface with no plan at all. This scenario is not merely
// untested, it is unconstructible.

// ---------------------------------------------------------------------
// Recoverable: configured memory exhaustion, restores accounting, never
// leaves the pool's used-memory elevated after wipe.
// ---------------------------------------------------------------------

template <int NLayers>
void checkMemoryLimitIsRecoverableAndRestoresAccounting()
{
  OneLayerDecoder* decoder = nullptr;
  auto interfacePtr = makeReadyInterface<NLayers>(decoder);
  auto& interface = *interfacePtr;
  auto& pool = *interface.getTimeFrame().getMemoryPool();
  // Release plan-adoption allocations so the test starts at the same empty
  // per-event baseline as the load attempt being bounded. The adopted plan's
  // runtime count remains attached to the scratch for the retry.
  interface.resetEvent();
  const size_t usedBeforeFailure = pool.getUsedMemory();
  // Plan adoption reserves the fixed scratch topology before event loading.
  // Cap at the already-accounted baseline so the first event allocation is
  // still rejected without attempting to set a maximum below current use.
  pool.setMaxMemory(usedBeforeFailure);

  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const auto patterns = makePatternBytes(clusters.size());

  bool threw = false;
  try {
    interface.processTimeFrame(rofs, clusters, patterns, nullptr);
    BOOST_FAIL("expected a memory-limit exception");
  } catch (const BoundedMemoryResource::MemoryLimitExceeded&) {
    threw = true;
  }
  BOOST_CHECK(threw);
  BOOST_CHECK_EQUAL(pool.getUsedMemory(), usedBeforeFailure);
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 0u);

  // A subsequent valid load succeeds once the bound is restored.
  pool.setMaxMemory(std::numeric_limits<size_t>::max());
  const float retried = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(retried));
}

BOOST_AUTO_TEST_CASE(MFT_MemoryLimitIsRecoverableAndRestoresAccounting) { checkMemoryLimitIsRecoverableAndRestoresAccounting<10>(); }

// ---------------------------------------------------------------------
// Unclassified: a decoder bug propagates by its original type, never
// dropped regardless of DropTFUponFailure.
// ---------------------------------------------------------------------

template <int NLayers>
void checkDecoderExceptionPropagatesByOriginalType()
{
  ThrowingDecoder* decoder = nullptr;
  auto interfacePtr = makeReadyInterface<NLayers, ThrowingDecoder>(decoder);
  auto& interface = *interfacePtr;
  const_cast<std::vector<o2::itsmft::TrackingParameters>&>(interface.getTrackingParameters())[0].DropTFUponFailure = true;

  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const auto patterns = makePatternBytes(clusters.size());

  BOOST_CHECK_EXCEPTION(interface.processTimeFrame(rofs, clusters, patterns, nullptr), std::runtime_error,
                        [](const std::runtime_error& e) { return std::string(e.what()).find("ThrowingDecoder") != std::string::npos; });
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 0u);
}

BOOST_AUTO_TEST_CASE(MFT_DecoderExceptionPropagatesByOriginalType) { checkDecoderExceptionPropagatesByOriginalType<10>(); }

// ---------------------------------------------------------------------
// No tracking callback runs past a load failure.
// ---------------------------------------------------------------------

template <int NLayers>
class CountingInterface : public ITSMFTTrackingInterface<NLayers>
{
 public:
  using ITSMFTTrackingInterface<NLayers>::ITSMFTTrackingInterface;
  int loadedCount = 0;
  int finishedCount = 0;

 protected:
  void onTimeFrameLoaded() override { ++loadedCount; }
  void onTrackingFinished(float) override { ++finishedCount; }
};

template <int NLayers>
void checkNoCallbacksPastLoadFailure()
{
  using Traits = TestTraits<NLayers>;
  auto decoderOwner = std::make_unique<OneLayerDecoder>(Traits::detId, Traits::kind);
  CountingInterface<NLayers> interface{false, o2::itsmft::TrackingMode::Sync, false, std::move(decoderOwner)};
  interface.initialise();
  interface.setClusterDictionary(&dict());
  auto& context = testAdapterContext<NLayers>();
  if constexpr (NLayers == ITSNLayers) {
    context.adapter.adoptITSSharedClusterCompatibility(&context.itsSidecar);
  } else {
    context.adapter.adoptMFTPublicationCompatibility(&context.mftSidecar);
  }
  interface.bindPublicationAdapter(context.adapter);
  context.configure();
  interface.bindROFViews(context.views);
  interface.getScratch().setROFViews(context.views);

  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const std::vector<unsigned char> malformedPattern{0, 1};

  BOOST_CHECK_THROW(interface.processTimeFrame(rofs, clusters, malformedPattern, nullptr), RecoverableLoadFailure);
  BOOST_CHECK_EQUAL(interface.loadedCount, 0);
  BOOST_CHECK_EQUAL(interface.finishedCount, 0);

  const auto patterns = makePatternBytes(clusters.size());
  const float result = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(result));
  BOOST_CHECK_EQUAL(interface.loadedCount, 1);
  BOOST_CHECK_EQUAL(interface.finishedCount, 1);
}

BOOST_AUTO_TEST_CASE(MFT_NoCallbacksPastLoadFailure) { checkNoCallbacksPastLoadFailure<10>(); }

// ---------------------------------------------------------------------
// MFT thread configuration is unaffected by ITS's dedicated nThreads.
//
// ITSMFTTrackingInterface<NLayers>::initialiseTracker() now branches on
// DetId to source nThreads (TrackingInterface.cxx): ITS reads its own
// ITSCommonCATrackerParam.nThreads (testITSCommonCANThreads.cxx), MFT is
// untouched and still reads TrackerParamConfig<MFT>.nThreads via
// TrackerParamRef<MFT>::get() -- confirmed here through the real interface,
// not just by inspection.
// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(MFT_ThreadConfigurationUnaffectedByITSDedicatedNThreads)
{
  auto& mftParam = const_cast<TrackerParamConfig<o2::detectors::DetID::MFT>&>(
    TrackerParamConfig<o2::detectors::DetID::MFT>::Instance());
  const int original = mftParam.nThreads;
  struct RestoreNThreads {
    int& field;
    int original;
    ~RestoreNThreads() { field = original; }
  } restoreGuard{mftParam.nThreads, original};
  mftParam.nThreads = 3;

  OneLayerDecoder* decoder = nullptr;
  auto interfacePtr = makeReadyInterface<10>(decoder);
  BOOST_CHECK_EQUAL(interfacePtr->getTrackerNThreads(), 3);
}
