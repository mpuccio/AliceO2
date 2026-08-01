// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Interface-level (ITSMFTTrackingInterface::processTimeFrame()) loading
// failure-contract tests. No real ITS/MFT geometry singleton is required: a
// ClusterDecoder is injected via the Slice 5 constructor overload instead of
// the production GeometryClusterDecoder<DetId>. A minimal, self-contained
// GRPECSObject is injected through GRPGeomHelper::finaliseCCDB() (the same
// entry point a real DPL device uses from its own finaliseCCDB() method) so
// that MFT's configureROFLookupTables() -- the only per-detector code path
// that dereferences GRPGeomHelper::instance().getGRPECS() -- does not need a
// running DPL topology or CCDB access.
//
// MFT only, not ITS: every check function below is templated over NLayers
// and works unchanged for NLayers=7, but ITSMFTTrackingInterface<7>::
// initialise() cannot succeed on this branch (or on main) for anyone, test
// or production. o2::itsmft::TrackingMode::getTrackingParameters() --
// Configuration.cxx, pre-existing, unrelated to this loading-boundary
// correction -- unconditionally LOGP(fatal, ...)s for
// detId==o2::detectors::DetID::ITS ("ITS CA tracking via O2::ITSMFTTracking
// is not enabled yet; use O2::ITStracking"), regardless of tracking mode.
// This is confirmed empirically: instantiating every check<7>() here (kept
// in an earlier revision of this file) reliably hit that fatal through
// resolveTrackingParameters(). ITS opt-in onboarding is a separate,
// not-yet-started piece of work (AgentCoordination.md Gate 3 status: "ITS
// parameter/workflow onboarding ... follow[s]"), out of scope for this
// bounded correction per binding requirement #10 ("do not start any other
// feature"). The loading-boundary machinery this file exercises is already
// covered for LegacyTrackerScratch<7> independently of
// ITSMFTTrackingInterface: see testTimeFrameNormalizedSource.cxx's
// ITSTimeFrameNormalizedSourceParity
// (LegacyTrackerScratch<7>::loadNormalizedSource() directly) and
// testTimeFrameLoadFailure.cxx's exhaustiveness/exception tests (NLayers-
// independent). Only MFT has a real opt-in ITSMFTTrackingInterface consumer
// today (CATrackerSpec.cxx / o2-mft-ca-tracker-workflow), so MFT-only
// coverage here matches this design's actual current deployment scope.

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
    obj.addDetContinuousReadOut(o2::detectors::DetID::MFT); // configureROFLookupTables()'s MFT-only continuous-readout branch
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
// guaranteed to be usable either (LegacyTrackerScratch<NLayers> is not
// trivially movable); returning it by value would require relying on unguaranteed
// NRVO through several mutating calls in between.
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
// that used to live here (skip configureDetectorLayouts() to produce
// MultiSourceLoadError::SurfaceCatalogNotConfigured): configureDetectorLayouts()
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
  const size_t usedBeforeFailure = pool.getUsedMemory();
  pool.setMaxMemory(1); // far below what even one cluster's backfill needs

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
// Non-uniform per-layer ROF timing (configureROFLookupTables()) is
// structural: it must throw TimeFrameLoadException{NonUniformROFTiming},
// never a dropped-TF sentinel, even with DropTFUponFailure=true; it must
// wipe existing TimeFrame event state; and a subsequent valid load must
// succeed once the configuration is restored. ITS is not covered here (see
// this file's header comment): only MFT has a real opt-in
// ITSMFTTrackingInterface consumer today, and DPLAlpideParam<MFT> is the
// singleton this test perturbs.
// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(MFT_NonUniformROFTimingIsStructuralNeverDroppedAndWipes)
{
  OneLayerDecoder* decoder = nullptr;
  auto interfacePtr = makeReadyInterface<10>(decoder);
  auto& interface = *interfacePtr;
  // Structural failures are never gated by DropTFUponFailure: set it so a
  // dropped-sentinel return (rather than a thrown exception) would be an
  // observable, wrong outcome if this were misclassified as recoverable.
  const_cast<std::vector<o2::itsmft::TrackingParameters>&>(interface.getTrackingParameters())[0].DropTFUponFailure = true;

  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const auto patterns = makePatternBytes(clusters.size());

  // Baseline valid load first, so the later wipe check proves state was
  // actually cleared, not merely "still empty".
  const float baseline = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_REQUIRE(!isDroppedTimeFrame(baseline));
  BOOST_REQUIRE_EQUAL(interface.getScratch().getTotalClusters(), 1u);

  {
    // MFT default roFrameLengthInBC is LHCMaxBunches/18 = 198; overriding
    // one layer to 202 changes that layer's own nROFsPerOrbit (3564/202=17)
    // away from every other layer's (3564/198=18), so mNROFsTF -- not only
    // mROFLength -- differs across layers, exactly the case the removed
    // per-layer mNROFsTF fatal check used to catch, now reachable only
    // through deriveUniformROFTimingConfig()'s mROFLength comparison.
    constexpr int overriddenLayer = 3;
    constexpr int overrideROFLengthInBC = 202;
    ScopedMFTLayerROFLengthOverride guard{overriddenLayer, overrideROFLengthInBC};

    bool threw = false;
    try {
      interface.processTimeFrame(rofs, clusters, patterns, nullptr);
      BOOST_FAIL("expected TimeFrameLoadException{NonUniformROFTiming}");
    } catch (const TimeFrameLoadException& err) {
      threw = true;
      BOOST_CHECK(err.reason() == TimeFrameLoadFailureReason::NonUniformROFTiming);
    }
    BOOST_CHECK(threw);
    // Wiped: the baseline TF's clusters are gone, not merely never-replaced.
    BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 0u);
  } // guard restores DPLAlpideParam<MFT>::Instance().roFrameLayerLengthInBC[3] here, even though the checks above never threw past it

  // Catalog/layout/detId are untouched by the failure (Slice 6/7 contract);
  // a subsequent valid load succeeds on the same interface once the
  // configuration is uniform again.
  const float retried = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(retried));
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 1u);
}

BOOST_AUTO_TEST_CASE(MFT_NonPositiveROFLengthIsStructuralBeforeDivision)
{
  OneLayerDecoder* decoder = nullptr;
  auto interfacePtr = makeReadyInterface<10>(decoder);
  auto& interface = *interfacePtr;

  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const auto patterns = makePatternBytes(clusters.size());

  {
    // roFrameLayerLengthInBC[layer] == 0 makes getROFLengthInBC(layer) fall
    // back to the shared default (it is a "use the global value" sentinel,
    // not literally zero), so a non-positive *effective* ROF length can
    // only be reached by overriding the shared default itself. This proves
    // configureROFLookupTables() checks positivity before ever computing
    // LHCMaxBunches/rofLengthInBC (a division that would otherwise be
    // undefined behavior for a zero divisor) rather than crashing.
    auto& par = mutableMFTAlpideParam();
    const int originalShared = par.roFrameLengthInBC;
    struct RestoreShared {
      int& field;
      int original;
      ~RestoreShared() { field = original; }
    } restoreGuard{par.roFrameLengthInBC, originalShared};
    par.roFrameLengthInBC = 0;

    bool threw = false;
    try {
      interface.processTimeFrame(rofs, clusters, patterns, nullptr);
      BOOST_FAIL("expected TimeFrameLoadException{NonUniformROFTiming}");
    } catch (const TimeFrameLoadException& err) {
      threw = true;
      BOOST_CHECK(err.reason() == TimeFrameLoadFailureReason::NonUniformROFTiming);
    }
    BOOST_CHECK(threw);
  }

  const float retried = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(retried));
}

// ---------------------------------------------------------------------
// A per-layer ROF length exceeding LHCMaxBunches (3564) makes
// nROFsPerOrbit (== LHCMaxBunches / rofLengthInBC) integer-divide to 0,
// so mNROFsTF == 0: a real, malformed per-TF timing configuration with no
// ROF at all to anchor a diamond-vertex TF interval envelope on (see
// TrackerTraits::computeLayerTrackletsForPolicy). This must fail
// structurally (TimeFrameLoadException{ZeroROFCount}) before any
// tracklet/cell code ever indexes ROF 0 or ROF mNROFsTF-1, exactly like
// every other malformed per-layer timing case above -- never silently
// treated as "0 ROFs to track", and never a crash.
// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(MFT_ZeroROFCountFromOversizedROFLengthIsStructural)
{
  OneLayerDecoder* decoder = nullptr;
  auto interfacePtr = makeReadyInterface<10>(decoder);
  auto& interface = *interfacePtr;

  const auto rofs = oneRof();
  const auto clusters = oneCluster();
  const auto patterns = makePatternBytes(clusters.size());

  {
    constexpr int overriddenLayer = 0;
    constexpr int oversizedROFLengthInBC = 4000; // > LHCMaxBunches (3564)
    ScopedMFTLayerROFLengthOverride guard{overriddenLayer, oversizedROFLengthInBC};

    bool threw = false;
    try {
      interface.processTimeFrame(rofs, clusters, patterns, nullptr);
      BOOST_FAIL("expected TimeFrameLoadException{ZeroROFCount}");
    } catch (const TimeFrameLoadException& err) {
      threw = true;
      BOOST_CHECK(err.reason() == TimeFrameLoadFailureReason::ZeroROFCount);
    }
    BOOST_CHECK(threw);
    BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 0u);
  }

  const float retried = interface.processTimeFrame(rofs, clusters, patterns, nullptr);
  BOOST_CHECK(!isDroppedTimeFrame(retried));
  BOOST_CHECK_EQUAL(interface.getScratch().getTotalClusters(), 1u);
}

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
