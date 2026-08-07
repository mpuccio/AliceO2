// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Stage-B Slice A: exactly-once systematic covariance in the common
// normalized-loading path. SurfaceMeasurement covariance is authoritative;
// compatibility TrackingFrameInfo copies it, and TimeFrame initialization
// never mutates either representation.

#define BOOST_TEST_MODULE ITSMFT TimeFrameCovarianceLifecycle
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <bit>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "GPUCommonMath.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITStracking/Constants.h"
#include "ITStracking/ROFLookupTables.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

constexpr std::array<unsigned char, 3> OnePixelPattern{1, 1, 0x80};
constexpr float ConfiguredRowIncrement{0.125f};
constexpr float ConfiguredColumnIncrement{0.25f};

const TopologyDictionary& dictionary()
{
  static const TopologyDictionary value;
  return value;
}

uint32_t floatBits(float value)
{
  return std::bit_cast<uint32_t>(value);
}

void checkBitIdentical(float lhs, float rhs)
{
  BOOST_CHECK_EQUAL(floatBits(lhs), floatBits(rhs));
}

// Geometry-free decoder seam with the production covariance contract:
// extract the real base covariance, then apply the configured detector
// increments once if and only if loadNormalizedSource() requests it.
class SystematicContractDecoder final : public ClusterDecoder
{
 public:
  SystematicContractDecoder(o2::detectors::DetID::ID detector, SurfaceKind kind,
                            float rowIncrement, float columnIncrement)
    : mDetector(detector), mKind(kind), mRowIncrement(rowIncrement), mColumnIncrement(columnIncrement)
  {
  }

  o2::itsmft::tracking::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const TopologyDictionary* dict,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool applySysErrors) const override
  {
    lastApplySysErrors = applySysErrors;
    const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dict);
    if (!clusterData.ok()) {
      o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
      result.error = clusterData.error;
      return result;
    }

    o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
    const int layer = cluster.getSensorID();
    result.layer = layer;
    if (layer < 0 || static_cast<size_t>(layer) >= layerToSurface.size()) {
      result.error = ClusterDecodeError::InvalidLayerMapping;
      return result;
    }
    result.layerMapped = true;
    result.kind = mKind;

    const float rowCovariance = clusterData.sig2Row + (applySysErrors ? mRowIncrement : 0.f);
    const float columnCovariance = clusterData.sig2Col + (applySysErrors ? mColumnIncrement : 0.f);
    DecodedCluster decoded{
      {3.f, 4.f, 5.f},
      {6.f, 7.f, 8.f, 0.125f},
      {rowCovariance, 0.f, columnCovariance},
      clusterData.shape,
      static_cast<uint32_t>(cluster.getSensorID()),
      layer};

    const DetectorSensorId sensor{static_cast<uint32_t>(mDetector), decoded.sensor};
    const ClusterRef clusterRef{source, externalIndex};
    result.measurement = mKind == SurfaceKind::Disk
                           ? makeDiskSurfaceMeasurement(decoded, sensor, layerToSurface[layer], clusterRef, sourceROF)
                           : makeCylinderSurfaceMeasurement(decoded, sensor, layerToSurface[layer], clusterRef, sourceROF);
    return result;
  }

  mutable bool lastApplySysErrors{false};

 private:
  o2::detectors::DetID::ID mDetector;
  SurfaceKind mKind;
  float mRowIncrement;
  float mColumnIncrement;
};

template <int NLayers>
std::vector<SurfaceId> identitySurfaces()
{
  std::vector<SurfaceId> result;
  result.reserve(NLayers);
  for (uint16_t layer = 0; layer < NLayers; ++layer) {
    result.emplace_back(layer);
  }
  return result;
}

template <int NLayers>
std::vector<SurfaceDescriptor> makeCatalog(o2::detectors::DetID::ID detector, SurfaceKind kind)
{
  std::vector<SurfaceDescriptor> result;
  result.reserve(NLayers);
  for (uint16_t layer = 0; layer < NLayers; ++layer) {
    result.push_back(SurfaceDescriptor{
      SurfaceId{layer}, layer, static_cast<uint8_t>(detector), kind, 0,
      static_cast<float>(layer + 1), 0.f, 100.f});
    const float xOverX0 = detector == o2::detectors::DetID::ITS
                            ? kNominalITSLayerX0[layer]
                            : kNominalMFTLayerX0[layer];
    result.back().material.xOverX0 = xOverX0;
    result.back().material.arealDensityGPerCm2 =
      xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
  }
  return result;
}

template <int NLayers>
std::vector<TrackingParameters> makeLifecycleParameters(o2::detectors::DetID::ID detector,
                                                        float rowIncrement,
                                                        float columnIncrement)
{
  std::vector<TrackingParameters> result(3);
  for (auto& parameters : result) {
    resetDetectorDefaults(parameters, detector);
    parameters.SystError2Row.assign(NLayers, rowIncrement);
    parameters.SystError2Col.assign(NLayers, columnIncrement);
  }
  result[1].PassFlags.reset(); // later iteration, LUT reuse
  result[2].PassFlags = IterationSteps{IterationStep::RebuildClusterLUT};
  return result;
}

struct CovarianceSnapshot {
  std::array<uint32_t, 3> normalized{};
  std::array<uint32_t, 3> compatibility{};

  friend bool operator==(const CovarianceSnapshot&, const CovarianceSnapshot&) = default;
};

CovarianceSnapshot snapshotCovariance(const TimeFrame& frame, const SurfaceTrackingScratch& tf)
{
  const auto measurements = frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{0});
  BOOST_REQUIRE_EQUAL(measurements.size(), 1u);
  const auto& compatibility = tf.getTrackingFrameInfoOnLayer(0);
  BOOST_REQUIRE_EQUAL(compatibility.size(), 1u);
  return {
    {floatBits(measurements[0].covariance.uu),
     floatBits(measurements[0].covariance.uv),
     floatBits(measurements[0].covariance.vv)},
    {floatBits(compatibility[0].covarianceTrackingFrame[0]),
     floatBits(compatibility[0].covarianceTrackingFrame[1]),
     floatBits(compatibility[0].covarianceTrackingFrame[2])}};
}

void checkRepresentationsAligned(const TimeFrame& frame, const SurfaceTrackingScratch& tf)
{
  const auto snapshot = snapshotCovariance(frame, tf);
  BOOST_CHECK(snapshot.normalized == snapshot.compatibility);
}

template <int NLayers>
struct Rig {
  Rig(o2::detectors::DetID::ID detectorValue, SurfaceKind kindValue,
      float rowIncrement, float columnIncrement)
    : detector(detectorValue), kind(kindValue), decoder(detectorValue, kindValue, rowIncrement, columnIncrement), pool(std::make_shared<BoundedMemoryResource>())
  {
    frame.setMemoryPool(pool);
    tf.setMemoryPool(pool);
    traits.setMemoryPool(pool);
    traits.adoptScratch(&tf);
    traits.adoptFrame(&frame);
    traits.setBz(5.f);
  }

  void configure(gsl::span<const TrackingParameters> parameters)
  {
    catalog = makeCatalog<NLayers>(detector, kind);
    const auto orderedSurfaces = identitySurfaces<NLayers>();
    const SurfaceCatalogView catalogView{catalog.data(), static_cast<uint32_t>(catalog.size())};
    auto result = buildSurfaceGraphs(catalogView, orderedSurfaces, parameters);
    BOOST_REQUIRE(result.ok());
    plan.emplace(std::move(result.graphs));
    const auto layoutView = plan->front().getView();
    tf.adoptPlan(plan->front().getOrderedSurfaces().size(), layoutView.nTransitions, layoutView.nCells);
    SurfaceMask owned;
    for (const auto surface : plan->front().getOrderedSurfaces()) {
      owned.set(surface);
    }
    auto bindingResult = SurfacePlanBinding::build(layoutView, ClusterSourceId{0}, owned,
                                                   plan->front().getOrderedSurfaces(), kind);
    BOOST_REQUIRE(bindingResult.ok());
    traits.adoptSurfacePlanBinding(bindingResult.binding.get());
    binding = std::move(bindingResult.binding);
  }

  void load(bool applySysErrors)
  {
    const std::vector<CompClusterExt> clusters{
      CompClusterExt{10, 20, CompCluster::InvalidPatternID, 0}};
    const std::vector<unsigned char> patterns{OnePixelPattern.begin(), OnePixelPattern.end()};
    const std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, 1}};
    const auto& orderedSurfaces = plan->front().getOrderedSurfaces();
    const auto result = tf.loadNormalizedSource(
      frame, decoder, o2::InteractionRecord{50, 5}, ROFTimingConfig{40, 0, 0, 0},
      clusters, patterns, rofs, &dictionary(), nullptr, detector,
      gsl::span<const SurfaceId>{orderedSurfaces}, plan->front().getSurfaceCatalog(), applySysErrors);
    BOOST_REQUIRE(result.ok());

    o2::its::LayerTiming timing{};
    timing.mNROFsTF = 1;
    timing.mROFLength = 40;
    rofTable.emplace();
    for (int layer = 0; layer < NLayers; ++layer) {
      rofTable->defineLayer(layer, timing);
    }
    rofTable->init();
    vertexTable.emplace();
    for (int layer = 0; layer < NLayers; ++layer) {
      vertexTable->defineLayer(layer, timing);
    }
    vertexTable->init();

    mask.emplace(*rofTable);
    mask->resetMask();
    for (int layer = 0; layer < NLayers; ++layer) {
      mask->setROFsEnabled(layer, 0, 1, 1);
    }
    tf.setROFViews(RuntimeROFViews{rofTable->getView(), vertexTable->getView(), mask->getView(), {}});
  }

  o2::detectors::DetID::ID detector;
  SurfaceKind kind;
  SystematicContractDecoder decoder;
  std::shared_ptr<BoundedMemoryResource> pool;
  // Gate 4 B3.1: `frame` declared before `tf` so it is constructed first and
  // destroyed last (see SurfaceTrackingScratch's own lifetime-contract doc).
  TimeFrame frame;
  SurfaceTrackingScratch tf;
  TrackerTraits traits;
  // Scratch carries non-owning runtime ROF views. Keep these adapter-edge
  // builders alive for every subsequent tracking call.
  std::optional<o2::its::ROFOverlapTable<NLayers>> rofTable;
  std::optional<o2::its::ROFVertexLookupTable<NLayers>> vertexTable;
  std::optional<o2::its::ROFMaskTable<NLayers>> mask;
  // Must outlive `plan` (std::vector<SurfaceGraph> borrows a SurfaceCatalogView into
  // it, Gate 4 B2 Slice 2) -- declared before `plan` so it is constructed
  // first and destroyed last.
  std::vector<SurfaceDescriptor> catalog;
  std::optional<std::vector<SurfaceGraph>> plan;
  std::unique_ptr<SurfacePlanBinding> binding;
};

template <int NLayers>
void checkPositionResolution(const Rig<NLayers>& rig,
                             const TrackingParameters& parameters,
                             float rowIncrement,
                             float columnIncrement)
{
  const float expected = o2::gpu::CAMath::Sqrt(
    0.5f * (columnIncrement + rowIncrement) +
    parameters.LayerResolution[0] * parameters.LayerResolution[0]);
  checkBitIdentical(rig.tf.getPositionResolution(0), expected);
}

template <int NLayers>
void exerciseDisabled(o2::detectors::DetID::ID detector, SurfaceKind kind)
{
  auto parameters = makeLifecycleParameters<NLayers>(detector, 0.f, 0.f);
  Rig<NLayers> rig{detector, kind, 0.f, 0.f};
  rig.configure(parameters);
  rig.load(true);
  rig.traits.updateTrackingParameters(parameters);

  const auto before = snapshotCovariance(rig.frame, rig.tf);
  checkBitIdentical(std::bit_cast<float>(before.normalized[0]), o2::itsmft::ioutils::DefClusError2Row);
  checkBitIdentical(std::bit_cast<float>(before.normalized[2]), o2::itsmft::ioutils::DefClusError2Col);
  checkRepresentationsAligned(rig.frame, rig.tf);

  BOOST_REQUIRE_NO_THROW(rig.traits.initialiseTimeFrame(0, *rig.plan));
  BOOST_CHECK(snapshotCovariance(rig.frame, rig.tf) == before);
  checkRepresentationsAligned(rig.frame, rig.tf);
}

template <int NLayers>
void exerciseEnabledLifecycle(o2::detectors::DetID::ID detector, SurfaceKind kind)
{
  auto parameters = makeLifecycleParameters<NLayers>(
    detector, ConfiguredRowIncrement, ConfiguredColumnIncrement);
  Rig<NLayers> rig{
    detector, kind, ConfiguredRowIncrement, ConfiguredColumnIncrement};
  rig.configure(parameters);
  rig.load(true);
  BOOST_REQUIRE(rig.decoder.lastApplySysErrors);
  rig.traits.updateTrackingParameters(parameters);

  const float expectedRow = o2::itsmft::ioutils::DefClusError2Row + ConfiguredRowIncrement;
  const float expectedColumn = o2::itsmft::ioutils::DefClusError2Col + ConfiguredColumnIncrement;
  const auto loaded = snapshotCovariance(rig.frame, rig.tf);
  checkBitIdentical(std::bit_cast<float>(loaded.normalized[0]), expectedRow);
  checkBitIdentical(std::bit_cast<float>(loaded.normalized[2]), expectedColumn);
  checkRepresentationsAligned(rig.frame, rig.tf);

  // FirstPass + LUT rebuild.
  BOOST_REQUIRE_NO_THROW(rig.traits.initialiseTimeFrame(0, *rig.plan));
  BOOST_CHECK(snapshotCovariance(rig.frame, rig.tf) == loaded);
  checkPositionResolution(rig, parameters[0], ConfiguredRowIncrement, ConfiguredColumnIncrement);

  // Later iteration + LUT reuse.
  BOOST_REQUIRE_NO_THROW(rig.traits.initialiseTimeFrame(1, *rig.plan));
  BOOST_CHECK(snapshotCovariance(rig.frame, rig.tf) == loaded);

  // Later iteration + LUT rebuild.
  BOOST_REQUIRE_NO_THROW(rig.traits.initialiseTimeFrame(2, *rig.plan));
  BOOST_CHECK(snapshotCovariance(rig.frame, rig.tf) == loaded);

  // Repeated FirstPass initialization is bit-identical.
  BOOST_REQUIRE_NO_THROW(rig.traits.initialiseTimeFrame(0, *rig.plan));
  BOOST_CHECK(snapshotCovariance(rig.frame, rig.tf) == loaded);
  BOOST_REQUIRE_NO_THROW(rig.traits.initialiseTimeFrame(0, *rig.plan));
  BOOST_CHECK(snapshotCovariance(rig.frame, rig.tf) == loaded);
  checkRepresentationsAligned(rig.frame, rig.tf);

  // A rejected later iteration leaves normalized and compatibility
  // covariance untouched.
  auto invalidParameters = parameters;
  ++invalidParameters[1].RowBins;
  rig.traits.updateTrackingParameters(invalidParameters);
  bool rejected = false;
  try {
    rig.traits.initialiseTimeFrame(1, *rig.plan);
  } catch (const TraversalException& error) {
    rejected = true;
    BOOST_CHECK(error.getReason() == TraversalFailureReason::IndexTableConfigurationMismatch);
  }
  BOOST_REQUIRE(rejected);
  BOOST_CHECK(snapshotCovariance(rig.frame, rig.tf) == loaded);
  checkRepresentationsAligned(rig.frame, rig.tf);

  // The frame-owned reset clears event state but preserves the configured
  // layout. Reloading and initializing the same event starts again from one
  // decoded increment, never from the previous compatibility copy.
  rig.tf.reset();
  rig.frame.resetEvent();
  rig.load(true);
  rig.traits.updateTrackingParameters(parameters);
  const auto reloaded = snapshotCovariance(rig.frame, rig.tf);
  BOOST_CHECK(reloaded == loaded);
  BOOST_REQUIRE_NO_THROW(rig.traits.initialiseTimeFrame(0, *rig.plan));
  BOOST_CHECK(snapshotCovariance(rig.frame, rig.tf) == loaded);
  checkRepresentationsAligned(rig.frame, rig.tf);
}

template <int NLayers>
void exerciseExplicitFalse(o2::detectors::DetID::ID detector, SurfaceKind kind)
{
  auto parameters = makeLifecycleParameters<NLayers>(
    detector, ConfiguredRowIncrement, ConfiguredColumnIncrement);
  Rig<NLayers> rig{
    detector, kind, ConfiguredRowIncrement, ConfiguredColumnIncrement};
  rig.configure(parameters);
  rig.load(false);
  BOOST_REQUIRE(!rig.decoder.lastApplySysErrors);
  rig.traits.updateTrackingParameters(parameters);

  const auto loaded = snapshotCovariance(rig.frame, rig.tf);
  checkBitIdentical(std::bit_cast<float>(loaded.normalized[0]), o2::itsmft::ioutils::DefClusError2Row);
  checkBitIdentical(std::bit_cast<float>(loaded.normalized[2]), o2::itsmft::ioutils::DefClusError2Col);
  checkRepresentationsAligned(rig.frame, rig.tf);

  BOOST_REQUIRE_NO_THROW(rig.traits.initialiseTimeFrame(0, *rig.plan));
  BOOST_CHECK(snapshotCovariance(rig.frame, rig.tf) == loaded);
  BOOST_REQUIRE_NO_THROW(rig.traits.initialiseTimeFrame(1, *rig.plan));
  BOOST_CHECK(snapshotCovariance(rig.frame, rig.tf) == loaded);
  BOOST_REQUIRE_NO_THROW(rig.traits.initialiseTimeFrame(2, *rig.plan));
  BOOST_CHECK(snapshotCovariance(rig.frame, rig.tf) == loaded);

  // Position-resolution preparation intentionally still consumes configured
  // systematic contributions even when this explicit loading override keeps
  // cluster covariance at its base value.
  checkPositionResolution(rig, parameters[2], ConfiguredRowIncrement, ConfiguredColumnIncrement);
  checkRepresentationsAligned(rig.frame, rig.tf);
}

} // namespace

BOOST_AUTO_TEST_CASE(SystematicsDisabledPreservesBaseCovarianceForITSAndMFT)
{
  exerciseDisabled<ITSNLayers>(o2::detectors::DetID::ITS, SurfaceKind::Cylinder);
  exerciseDisabled<MFTNLayers>(o2::detectors::DetID::MFT, SurfaceKind::Disk);
}

BOOST_AUTO_TEST_CASE(SystematicsEnabledRemainExactlyOnceAcrossTheFullLifecycle)
{
  exerciseEnabledLifecycle<ITSNLayers>(o2::detectors::DetID::ITS, SurfaceKind::Cylinder);
  exerciseEnabledLifecycle<MFTNLayers>(o2::detectors::DetID::MFT, SurfaceKind::Disk);
}

BOOST_AUTO_TEST_CASE(ExplicitApplySysErrorsFalseSurvivesInitialization)
{
  exerciseExplicitFalse<ITSNLayers>(o2::detectors::DetID::ITS, SurfaceKind::Cylinder);
  exerciseExplicitFalse<MFTNLayers>(o2::detectors::DetID::MFT, SurfaceKind::Disk);
}
