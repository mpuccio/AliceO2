// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT TimeFrame detector layouts
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <TGeoGlobalMagField.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"
#include "Field/MagneticField.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingInterface.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITStracking/Constants.h"

using namespace o2::itsmft::tracking;
using o2::itsmft::TrackingParameters;

struct TraversalPropagatorFieldFixture {
  TraversalPropagatorFieldFixture()
  {
    if (!TGeoGlobalMagField::Instance()->GetField()) {
      TGeoGlobalMagField::Instance()->SetField(o2::field::MagneticField::createNominalField(5, true));
      TGeoGlobalMagField::Instance()->Lock();
    }
  }
};

BOOST_GLOBAL_FIXTURE(TraversalPropagatorFieldFixture);

namespace
{
class FakeCatalogProvider final : public DetectorSurfaceCatalogProvider
{
 public:
  explicit FakeCatalogProvider(std::vector<SurfaceDescriptor> catalog) : mCatalog{std::move(catalog)} {}

  DetectorSurfaceCatalogResult buildCatalog(const DetectorSurfaceCatalogRequest& request) const final
  {
    ++calls;
    requests.push_back(request);
    return fail ? DetectorSurfaceCatalogResult{{}, DetectorSurfaceCatalogError::GeometryUnavailable}
                : DetectorSurfaceCatalogResult{mCatalog, DetectorSurfaceCatalogError::None};
  }

  mutable int calls{0};
  mutable std::vector<DetectorSurfaceCatalogRequest> requests;
  bool fail{false};
  std::vector<SurfaceDescriptor> mCatalog;
};

struct EpochTestTimeFrame : TimeFrame<7> {
  void setRequiredEpoch(DetectorGeometryEpoch epoch) { mRequiredDetectorGeometryEpoch = epoch; }
  const DetectorLayoutSet* getStoredDetectorLayouts() const { return mDetectorLayouts ? &*mDetectorLayouts : nullptr; }
};

/// DetectorLayout no longer owns a surface copy (Slice 3, shared ownership):
/// test fixtures that build one in isolation keep the surfaces alongside it.
struct BuiltLayout {
  DetectorLayout layout;
  std::vector<SurfaceDescriptor> surfaces;
};

struct TraversalTestTimeFrame : TimeFrame<10> {
  void installLayout(BuiltLayout built)
  {
    DetectorLayoutConfigurationKey key;
    key.geometryEpoch = mRequiredDetectorGeometryEpoch;
    std::vector<DetectorLayout> layouts;
    layouts.push_back(std::move(built.layout));
    mRequiredDetectorLayoutConfiguration = key;
    mDetectorLayouts.emplace(std::move(key), std::move(built.surfaces), std::move(layouts));
  }
};

static_assert(std::is_nothrow_move_constructible_v<DetectorLayoutSet>);

DetectorSurfaceCatalogRequest request(uint32_t count,
                                      o2::detectors::DetID::ID detector = o2::detectors::DetID::ITS,
                                      uint16_t firstSurface = 0)
{
  return DetectorSurfaceCatalogRequest{detector, SurfaceId{firstSurface}, count};
}

// Nominal material matching this detector's default TrackingParameters::LayerxX0
// (o2::itsmft::resetDetectorDefaults()), so fixtures that reach
// TrackerTraits::initialiseTimeFrame() with unperturbed parameters satisfy the
// LegacyMaterialMismatch compatibility check by construction. Indexed by
// global surface id, which equals the legacy layer index for every identity-
// ordered fixture in this file.
float nominalXOverX0(o2::detectors::DetID::ID detector, uint16_t surfaceIndex)
{
  if (detector == o2::detectors::DetID::MFT) {
    return kNominalMFTLayerX0[surfaceIndex % MFTNLayers];
  }
  return kNominalITSLayerX0[surfaceIndex % ITSNLayers];
}

std::vector<SurfaceDescriptor> catalog(size_t count, SurfaceKind kind = SurfaceKind::Cylinder,
                                       o2::detectors::DetID::ID detector = o2::detectors::DetID::ITS)
{
  std::vector<SurfaceDescriptor> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(detector), kind, 0, static_cast<float>(i + 1), 0.f, 100.f});
    const float xOverX0 = nominalXOverX0(detector, i);
    result.back().material.xOverX0 = xOverX0;
    result.back().material.arealDensityGPerCm2 = xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
  }
  return result;
}

std::vector<SurfaceDescriptor> catalog(const DetectorSurfaceCatalogRequest& catalogRequest,
                                       SurfaceKind kind = SurfaceKind::Cylinder)
{
  const auto size = catalogRequest.firstSurface.value() + catalogRequest.detectorSurfaceCount;
  auto result = catalog(size, kind, catalogRequest.detector);
  for (uint32_t localIndex = 0; localIndex < catalogRequest.detectorSurfaceCount; ++localIndex) {
    auto& surface = result[catalogRequest.firstSurface.value() + localIndex];
    surface.detectorId = static_cast<uint8_t>(catalogRequest.detector);
    surface.detectorSurfaceIndex = localIndex;
  }
  return result;
}

std::vector<SurfaceId> order(size_t count)
{
  std::vector<SurfaceId> result;
  for (uint16_t i = 0; i < count; ++i) {
    result.emplace_back(i);
  }
  return result;
}

BuiltLayout cyclicDiskLayout()
{
  SparseTrackingTopology topology{10};
  BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{0}, SurfaceId{1}, {}, TransitionPolicyTag::DiskDisk, 0}).isValid());
  BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{1}, SurfaceId{2}, {}, TransitionPolicyTag::DiskDisk, 0}).isValid());
  BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{2}, SurfaceId{0}, {}, TransitionPolicyTag::DiskDisk, 0}).isValid());
  BOOST_REQUIRE(topology.finalize());
  auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  return BuiltLayout{DetectorLayout{surfaces, std::move(topology)}, std::move(surfaces)};
}

BuiltLayout mixedDisconnectedLayout()
{
  SparseTrackingTopology topology{10};
  for (uint16_t id = 0; id < 4; ++id) {
    BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{id}, SurfaceId{static_cast<uint16_t>(id + 1)}, {}, TransitionPolicyTag::CylinderCylinder, 0}).isValid());
  }
  for (uint16_t id = 5; id < 9; ++id) {
    BOOST_REQUIRE(topology.addTransition(SurfaceTransition{SurfaceId{id}, SurfaceId{static_cast<uint16_t>(id + 1)}, {}, TransitionPolicyTag::DiskDisk, 0}).isValid());
  }
  BOOST_REQUIRE(topology.finalize());
  auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  for (uint16_t id = 0; id < 5; ++id) {
    surfaces[id].kind = SurfaceKind::Cylinder;
  }
  return BuiltLayout{DetectorLayout{surfaces, std::move(topology)}, std::move(surfaces)};
}

TrackingParameters parameters(int activeCount, int maxHoles = 0, uint16_t holes = 0, uint16_t starts = 0xffff)
{
  TrackingParameters result;
  result.NLayers = activeCount;
  result.MaxHoles = maxHoles;
  result.HoleLayerMask = holes;
  result.StartLayerMask = starts;
  return result;
}

template <int NLayers>
void prepareTraversalFrame(TimeFrame<NLayers>& frame,
                           TrackerTraits<NLayers>& traits,
                           const std::shared_ptr<BoundedMemoryResource>& pool,
                           const std::vector<TrackingParameters>& params)
{
  frame.setMemoryPool(pool);
  for (auto& rofOffsets : frame.mROFramesClusters) {
    rofOffsets.resize(1, 0);
  }
  frame.initTrackerTopologies(params);
  traits.setMemoryPool(pool);
  traits.adoptTimeFrame(&frame);
  traits.updateTrackingParameters(params);
}

std::vector<TrackingParameters> mftTraversalParameters()
{
  std::vector<TrackingParameters> params(1);
  o2::itsmft::resetDetectorDefaults(params.front(), o2::detectors::DetID::MFT);
  return params;
}

template <int NLayers>
void checkLegacyParity(SurfaceKind kind, TransitionPolicyTag policyTag, uint16_t startMask)
{
  TimeFrame<NLayers> frame;
  const auto catalogRequest = request(NLayers, kind == SurfaceKind::Disk ? o2::detectors::DetID::MFT : o2::detectors::DetID::ITS);
  FakeCatalogProvider provider{catalog(catalogRequest, kind)};
  auto ordered = order(NLayers);
  std::vector<TrackingParameters> params{parameters(NLayers, 1, uint16_t{1} << (NLayers / 2), startMask)};
  const auto result = frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, policyTag, params);
  BOOST_REQUIRE(result.ok());

  TrackingTopology<NLayers> legacy;
  legacy.init(NLayers, params[0].MaxHoles, params[0].HoleLayerMask);
  const auto legacyView = legacy.getView();
  const auto layoutView = frame.getDetectorLayoutView(0);
  const auto sparse = layoutView.topology;
  BOOST_REQUIRE_EQUAL(sparse.nTransitions, legacyView.nTransitions);
  BOOST_REQUIRE_EQUAL(sparse.nCells, legacyView.nCells);
  for (uint32_t i = 0; i < sparse.nTransitions; ++i) {
    BOOST_CHECK_EQUAL(sparse.transitions[i].from.value(), legacyView.transitions[i].fromLayer);
    BOOST_CHECK_EQUAL(sparse.transitions[i].to.value(), legacyView.transitions[i].toLayer);
  }
  // Independent legacy oracle for which CellTopologyIds the *former*
  // findRoadsForPolicy predicate (hitLayerMask.last() + StartLayerMask.has())
  // would have selected as road starts, recomputed here from the frozen
  // legacy TrackingTopology<NLayers> view -- not by calling into any
  // TransitionPolicyGrouping/TrackerTraits production code.
  std::vector<CellTopologyId> legacyOracleStarts;
  for (uint32_t i = 0; i < sparse.nCells; ++i) {
    const auto& sparseCell = sparse.cells[i];
    const auto& legacyCell = legacyView.cells[i];
    BOOST_CHECK_EQUAL(sparseCell.firstTransition.value(), legacyCell.firstTransition);
    BOOST_CHECK_EQUAL(sparseCell.secondTransition.value(), legacyCell.secondTransition);
    BOOST_CHECK_EQUAL(sparseCell.hitSurfaces.value(), static_cast<uint16_t>(legacyCell.hitLayerMask));
    const auto legacyStartLayer = legacyCell.hitLayerMask.last();
    const auto sparseStartSurface = SurfaceId{static_cast<uint16_t>(sparseCell.hitSurfaces.last())};
    BOOST_CHECK_EQUAL(sparse.seedingSurfaces.has(sparseStartSurface), params[0].StartLayerMask.has(legacyStartLayer));
    if (params[0].StartLayerMask.has(legacyStartLayer)) {
      legacyOracleStarts.push_back(CellTopologyId{static_cast<uint16_t>(i)});
    }
  }

  // Item 7: exact identity-layout parity between roadStartCellsForTag() and
  // the former StartLayerMask/hitLayerMask.last() selection, for both
  // TrackerTraits<7> (ITS-like) and TrackerTraits<10> (MFT-like) call sites
  // of this helper.
  TransitionPolicyGrouping grouping{layoutView};
  BOOST_REQUIRE(grouping.valid());
  const auto starts = grouping.roadStartCellsForTag(policyTag);
  BOOST_CHECK(std::is_sorted(starts.begin(), starts.end()));
  BOOST_REQUIRE_EQUAL(starts.size(), legacyOracleStarts.size());
  for (size_t i = 0; i < starts.size(); ++i) {
    BOOST_CHECK(starts[i] == legacyOracleStarts[i]);
  }
}
} // namespace

BOOST_AUTO_TEST_CASE(initial_dirty_missing_provider_and_bounds)
{
  TimeFrame<7> frame;
  const auto catalogRequest = request(7);
  auto ordered = order(7);
  std::vector<TrackingParameters> params{parameters(7)};
  BOOST_CHECK(!frame.detectorLayoutsCurrent());
  BOOST_CHECK(!frame.hasStoredDetectorLayouts());
  BOOST_CHECK(frame.getDetectorLayouts() == nullptr);
  BOOST_CHECK(frame.getSurfaceCatalog() == nullptr);
  BOOST_CHECK(frame.getSurfaceCatalogView().empty());
  BOOST_CHECK(frame.getDetectorLayout(0) == nullptr);
  BOOST_CHECK(frame.getDetectorLayoutView(0).surfaces == nullptr);
  const auto missing = frame.ensureDetectorLayouts(nullptr, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params);
  BOOST_CHECK(missing.error == DetectorLayoutSetBuildError::MissingProvider);
  BOOST_CHECK(!frame.detectorLayoutsCurrent());
}

BOOST_AUTO_TEST_CASE(initial_build_current_stale_bounds_and_wipe_preservation)
{
  TimeFrame<7> frame;
  const auto catalogRequest = request(7);
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{parameters(7), parameters(5)};
  const auto first = frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params);
  BOOST_REQUIRE(first.ok());
  BOOST_CHECK(first.rebuilt);
  BOOST_REQUIRE_EQUAL(provider.requests.size(), 1);
  BOOST_CHECK(provider.requests.front() == catalogRequest);
  BOOST_CHECK(frame.detectorLayoutsCurrent());
  BOOST_REQUIRE(frame.getSurfaceCatalog() != nullptr);
  BOOST_CHECK_EQUAL(frame.getSurfaceCatalog()->size(), 7);
  BOOST_CHECK_EQUAL(frame.getSurfaceCatalogView().size(), 7);
  BOOST_REQUIRE(frame.getDetectorLayouts() != nullptr);
  BOOST_CHECK_EQUAL(frame.getDetectorLayouts()->size(), 2);
  BOOST_CHECK(frame.getDetectorLayout(1) != nullptr);
  BOOST_CHECK(frame.getDetectorLayout(2) == nullptr);
  BOOST_CHECK(frame.getDetectorLayoutView(2).surfaces == nullptr);
  const auto current = frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params);
  BOOST_CHECK(current.ok());
  BOOST_CHECK(!current.rebuilt);
  BOOST_CHECK_EQUAL(provider.calls, 1);
  frame.wipe();
  BOOST_CHECK(frame.detectorLayoutsCurrent());
  BOOST_REQUIRE(frame.getSurfaceCatalog() != nullptr);
  BOOST_CHECK_EQUAL(frame.getSurfaceCatalogView().size(), 7);
  BOOST_CHECK(frame.getDetectorLayout(0) != nullptr);
  frame.invalidateDetectorLayouts();
  BOOST_CHECK(!frame.detectorLayoutsCurrent());
  BOOST_CHECK(frame.hasStoredDetectorLayouts());
  BOOST_CHECK(frame.getSurfaceCatalog() == nullptr);
  BOOST_CHECK(frame.getSurfaceCatalogView().empty());
  BOOST_CHECK(frame.getDetectorLayouts() == nullptr);
  BOOST_CHECK(frame.getDetectorLayoutView(0).surfaces == nullptr);
}

BOOST_AUTO_TEST_CASE(epoch_wrap_contract)
{
  BOOST_CHECK_EQUAL(nextDetectorGeometryEpoch(InitialDetectorGeometryEpoch), InitialDetectorGeometryEpoch + 1);
  BOOST_CHECK_EQUAL(nextDetectorGeometryEpoch(std::numeric_limits<DetectorGeometryEpoch>::max()), InitialDetectorGeometryEpoch);

  EpochTestTimeFrame frame;
  const auto catalogRequest = request(7);
  frame.setRequiredEpoch(std::numeric_limits<DetectorGeometryEpoch>::max());
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{parameters(7)};
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params).ok());
  BOOST_REQUIRE(frame.hasStoredDetectorLayouts());
  frame.invalidateDetectorLayouts();
  BOOST_CHECK_EQUAL(frame.getRequiredDetectorGeometryEpoch(), InitialDetectorGeometryEpoch);
  BOOST_CHECK(!frame.hasStoredDetectorLayouts());
  BOOST_CHECK(!frame.detectorLayoutsCurrent());
}

BOOST_AUTO_TEST_CASE(invalidation_replacement_and_provider_failure_are_transactional)
{
  EpochTestTimeFrame frame;
  const auto catalogRequest = request(7);
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{parameters(7)};
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params).ok());
  const auto* storedOwner = frame.getStoredDetectorLayouts();
  BOOST_REQUIRE(storedOwner != nullptr);
  const auto storedEpoch = storedOwner->getConfigurationKey().geometryEpoch;
  const auto firstEpoch = frame.getRequiredDetectorGeometryEpoch();
  BOOST_CHECK_EQUAL((*frame.getSurfaceCatalog())[0].referenceCoordinate, 1.f);
  frame.invalidateDetectorLayouts();
  BOOST_CHECK_EQUAL(frame.getRequiredDetectorGeometryEpoch(), firstEpoch + 1);
  provider.fail = true;
  const auto failed = frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params);
  BOOST_CHECK(failed.error == DetectorLayoutSetBuildError::CatalogProviderFailure);
  BOOST_CHECK(failed.catalogError == DetectorSurfaceCatalogError::GeometryUnavailable);
  BOOST_CHECK(frame.hasStoredDetectorLayouts());
  BOOST_CHECK(!frame.detectorLayoutsCurrent());
  BOOST_CHECK(frame.getStoredDetectorLayouts() == storedOwner);
  BOOST_CHECK_EQUAL(frame.getStoredDetectorLayouts()->getConfigurationKey().geometryEpoch, storedEpoch);
  provider.fail = false;
  provider.mCatalog[0].referenceCoordinate = 42.f;
  const auto replaced = frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params);
  BOOST_REQUIRE(replaced.ok());
  BOOST_CHECK(replaced.rebuilt);
  BOOST_CHECK(frame.detectorLayoutsCurrent());
  BOOST_CHECK_EQUAL((*frame.getSurfaceCatalog())[0].referenceCoordinate, 42.f);
}

BOOST_AUTO_TEST_CASE(request_change_forces_rebuild)
{
  TimeFrame<7> frame;
  auto catalogRequest = request(7);
  FakeCatalogProvider provider{catalog(catalogRequest)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{parameters(7)};
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params).ok());
  BOOST_CHECK_EQUAL(provider.calls, 1);

  catalogRequest = request(7, o2::detectors::DetID::MFT);
  provider.mCatalog = catalog(catalogRequest);
  const auto rebuilt = frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params);
  BOOST_REQUIRE(rebuilt.ok());
  BOOST_CHECK(rebuilt.rebuilt);
  BOOST_CHECK_EQUAL(provider.calls, 2);
  BOOST_CHECK(frame.getDetectorLayouts()->getConfigurationKey().catalogRequest == catalogRequest);
  BOOST_CHECK_EQUAL(frame.getSurfaceCatalogView().size(), 7);
}

BOOST_AUTO_TEST_CASE(catalog_request_validation)
{
  auto checkFailure = [](const DetectorSurfaceCatalogRequest& catalogRequest,
                         std::vector<SurfaceDescriptor> providedCatalog,
                         DetectorSurfaceCatalogValidationError expected) {
    TimeFrame<7> frame;
    FakeCatalogProvider provider{std::move(providedCatalog)};
    auto ordered = order(7);
    std::vector<TrackingParameters> noIterations;
    const auto result = frame.ensureDetectorLayouts(&provider, catalogRequest, ordered,
                                                    TransitionPolicyTag::CylinderCylinder, noIterations);
    BOOST_CHECK(result.error == DetectorLayoutSetBuildError::InvalidCatalog);
    BOOST_CHECK(result.catalogValidationError == expected);
    BOOST_CHECK(!frame.hasStoredDetectorLayouts());
  };

  auto wrongDetector = catalog(7, SurfaceKind::Cylinder, o2::detectors::DetID::MFT);
  checkFailure(request(7), std::move(wrongDetector), DetectorSurfaceCatalogValidationError::DetectorMismatch);
  checkFailure(request(6), catalog(7), DetectorSurfaceCatalogValidationError::SizeMismatch);
  checkFailure(request(6, o2::detectors::DetID::ITS, 1), catalog(7),
               DetectorSurfaceCatalogValidationError::DetectorSurfaceIndexOutOfRange);
  checkFailure(request(1, o2::detectors::DetID::ITS, MaxLayoutSurfaces), {},
               DetectorSurfaceCatalogValidationError::TooManySurfaces);

  auto duplicateIndex = catalog(7);
  duplicateIndex[6].detectorSurfaceIndex = 5;
  checkFailure(request(7), std::move(duplicateIndex), DetectorSurfaceCatalogValidationError::DuplicateDetectorSurfaceIndex);

  auto missingIndex = catalog(7);
  missingIndex[6].detectorSurfaceIndex = std::numeric_limits<uint16_t>::max();
  checkFailure(request(7), std::move(missingIndex), DetectorSurfaceCatalogValidationError::MissingDetectorSurfaceIndex);

  auto negativeXOverX0 = catalog(7);
  negativeXOverX0[2].material.xOverX0 = -1.f;
  checkFailure(request(7), std::move(negativeXOverX0), DetectorSurfaceCatalogValidationError::InvalidMaterial);

  auto nanXOverX0 = catalog(7);
  nanXOverX0[2].material.xOverX0 = std::numeric_limits<float>::quiet_NaN();
  checkFailure(request(7), std::move(nanXOverX0), DetectorSurfaceCatalogValidationError::InvalidMaterial);

  auto infAreal = catalog(7);
  infAreal[5].material.arealDensityGPerCm2 = std::numeric_limits<float>::infinity();
  checkFailure(request(7), std::move(infAreal), DetectorSurfaceCatalogValidationError::InvalidMaterial);

  auto negativeAreal = catalog(7);
  negativeAreal[5].material.arealDensityGPerCm2 = -0.1f;
  checkFailure(request(7), std::move(negativeAreal), DetectorSurfaceCatalogValidationError::InvalidMaterial);
}

BOOST_AUTO_TEST_CASE(malformed_catalog_rejected_with_zero_iterations)
{
  TimeFrame<7> frame;
  const auto catalogRequest = request(7);
  auto malformed = catalog(7);
  malformed[3].id = SurfaceId{4};
  FakeCatalogProvider provider{std::move(malformed)};
  auto ordered = order(7);
  std::vector<TrackingParameters> noIterations;
  const auto result = frame.ensureDetectorLayouts(&provider, catalogRequest, ordered,
                                                  TransitionPolicyTag::CylinderCylinder, noIterations);
  BOOST_CHECK(result.error == DetectorLayoutSetBuildError::InvalidCatalog);
  BOOST_CHECK(result.catalogValidationError == DetectorSurfaceCatalogValidationError::NonDenseGlobalSurfaceIds);
  BOOST_CHECK(!frame.hasStoredDetectorLayouts());
}

BOOST_AUTO_TEST_CASE(interface_owns_provider_from_construction)
{
  const auto catalogRequest = request(7);
  auto provider = std::make_unique<FakeCatalogProvider>(catalog(catalogRequest));
  auto* providerObserver = provider.get();
  ITSMFTTrackingInterface<7> interface{false, o2::itsmft::TrackingMode::Unset, false, std::move(provider)};
  auto ordered = order(7);
  const auto result = interface.configureDetectorLayouts(catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(result.rebuilt);
  BOOST_CHECK_EQUAL(providerObserver->calls, 1);
  BOOST_REQUIRE(interface.getTimeFrame().getSurfaceCatalog() != nullptr);
  BOOST_CHECK_EQUAL(interface.getTimeFrame().getSurfaceCatalogView().size(), 7);
  BOOST_CHECK(interface.getTimeFrame().getDetectorLayout(0) == nullptr);
}

BOOST_AUTO_TEST_CASE(final_iteration_builder_failure_commits_nothing)
{
  TimeFrame<7> frame;
  const auto catalogRequest = request(7);
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> valid{parameters(7)};
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, valid).ok());
  frame.invalidateDetectorLayouts();
  std::vector<TrackingParameters> failing{parameters(7), parameters(6), parameters(5, -1)};
  const auto result = frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, failing);
  BOOST_CHECK(result.error == DetectorLayoutSetBuildError::LayoutBuilderFailure);
  BOOST_CHECK_EQUAL(result.failedIteration, 2);
  BOOST_CHECK(result.layoutBuildError == DetectorLayoutBuildError::NegativeMaxHoles);
  BOOST_CHECK(frame.hasStoredDetectorLayouts());
  BOOST_CHECK(!frame.detectorLayoutsCurrent());
}

BOOST_AUTO_TEST_CASE(catalog_identity_active_count_and_mask_mapping)
{
  TimeFrame<7> frame;
  const auto catalogRequest = request(7);
  FakeCatalogProvider provider{catalog(7)};
  const std::vector<SurfaceId> ordered{SurfaceId{3}, SurfaceId{0}, SurfaceId{6}, SurfaceId{2}, SurfaceId{5}, SurfaceId{1}, SurfaceId{4}};
  std::vector<TrackingParameters> params{
    parameters(7, 1, uint16_t{1} << 1, (uint16_t{1} << 0) | (uint16_t{1} << 4)),
    parameters(4, 1, uint16_t{1} << 1, (uint16_t{1} << 0) | (uint16_t{1} << 4) | (uint16_t{1} << 6))};
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params).ok());
  const auto* full = frame.getDetectorLayout(0);
  const auto* reduced = frame.getDetectorLayout(1);
  BOOST_REQUIRE(full != nullptr && reduced != nullptr);
  // Slice 3: the surface catalog is now owned exactly once by
  // DetectorLayoutSet, shared by every iteration -- there is no longer a
  // separate per-DetectorLayout copy to compare element-by-element. Both
  // iterations trivially observe the same single catalog by construction.
  const auto* sharedCatalog = frame.getSurfaceCatalog();
  BOOST_REQUIRE(sharedCatalog != nullptr);
  BOOST_REQUIRE_EQUAL(sharedCatalog->size(), 7);
  BOOST_CHECK_LT(reduced->getTopology().getTransitions().size(), full->getTopology().getTransitions().size());
  BOOST_CHECK(full->getTopology().getView().seedingSurfaces.has(SurfaceId{3}));
  BOOST_CHECK(full->getTopology().getView().seedingSurfaces.has(SurfaceId{5}));
  BOOST_CHECK(!full->getTopology().getView().seedingSurfaces.has(SurfaceId{0}));
  BOOST_CHECK(reduced->getTopology().getView().seedingSurfaces.has(SurfaceId{3}));
  BOOST_CHECK_EQUAL(reduced->getTopology().getView().seedingSurfaces.count(), 1);
  const auto& reducedTransitions = reduced->getTopology().getTransitions();
  BOOST_REQUIRE(!reducedTransitions.empty());
  BOOST_CHECK(reducedTransitions.front().from == SurfaceId{3});
  BOOST_CHECK(reducedTransitions.front().to == SurfaceId{0});
  const auto skipped = std::find_if(reducedTransitions.begin(), reducedTransitions.end(), [](const auto& transition) {
    return transition.from == SurfaceId{3} && transition.to == SurfaceId{6};
  });
  BOOST_REQUIRE(skipped != reducedTransitions.end());
  BOOST_CHECK(skipped->skippedSurfaces.has(SurfaceId{0}));

  // Selected-iteration isolation (item 6): iterations 0 and 1 use different
  // StartLayerMask-derived seeding masks over the same catalog/ordering
  // (`ordered`: position 0 = SurfaceId 3, position 4 = SurfaceId 5).
  // Position 0 is a seeding surface in *both* iterations' masks but can
  // never be a cell's transition endpoint (transitions only run from an
  // earlier to a later position) -- "a seeded surface with no terminating
  // cell" in both layouts simultaneously. Iteration 0's other seeding
  // surface (5, position 4) is reachable within the full 7-active-surface
  // layout; iteration 1's activeCount=4 mask only ever tests position 0
  // (positions 4 and 6 of its starts bitset are out of range and ignored by
  // positionalSurfaceMask), so iteration 1's own layout must select no road
  // starts at all -- proving each iteration's roadStartCellsForTag reflects
  // only its own layout's seedingSurfaces.
  const auto fullView = frame.getDetectorLayoutView(0);
  const auto reducedView = frame.getDetectorLayoutView(1);
  TransitionPolicyGrouping fullGrouping{fullView};
  TransitionPolicyGrouping reducedGrouping{reducedView};
  BOOST_REQUIRE(fullGrouping.valid());
  BOOST_REQUIRE(reducedGrouping.valid());
  const auto fullStarts = fullGrouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  const auto reducedStarts = reducedGrouping.roadStartCellsForTag(TransitionPolicyTag::CylinderCylinder);
  BOOST_CHECK(std::is_sorted(fullStarts.begin(), fullStarts.end()));
  BOOST_CHECK(reducedStarts.empty());
  BOOST_REQUIRE_GT(fullStarts.size(), 0u);
  for (const auto id : fullStarts) {
    const auto endpoint = fullView.topology.getTransition(fullView.topology.getCell(id).secondTransition).to;
    BOOST_CHECK(endpoint == SurfaceId{5});
  }
}

BOOST_AUTO_TEST_CASE(traversal_initialisation_classifies_missing_and_stale_layouts)
{
  auto params = mftTraversalParameters();
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame<10> frame;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, traits, pool, params);

  try {
    traits.initialiseTimeFrame(0);
    BOOST_FAIL("missing layout must throw");
  } catch (const TraversalException& error) {
    BOOST_CHECK_EQUAL(error.getIteration(), 0);
    BOOST_CHECK(error.getReason() == TraversalFailureReason::MissingLayout);
  }
  BOOST_CHECK(!traits.hasTraversalCache());

  const auto catalogRequest = request(10, o2::detectors::DetID::MFT);
  FakeCatalogProvider provider{catalog(catalogRequest, SurfaceKind::Disk)};
  const auto ordered = order(10);
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered,
                                            TransitionPolicyTag::DiskDisk, params)
                  .ok());
  frame.invalidateDetectorLayouts();
  try {
    traits.initialiseTimeFrame(0);
    BOOST_FAIL("stale layout must throw");
  } catch (const TraversalException& error) {
    BOOST_CHECK(error.getReason() == TraversalFailureReason::StaleLayout);
  }
  BOOST_CHECK(!traits.hasTraversalCache());

  auto twoIterations = mftTraversalParameters();
  twoIterations.push_back(twoIterations.front());
  TimeFrame<10> shortLayoutFrame;
  TrackerTraits<10> shortLayoutTraits;
  prepareTraversalFrame(shortLayoutFrame, shortLayoutTraits, pool, twoIterations);
  std::vector<TrackingParameters> oneLayout{twoIterations.front()};
  BOOST_REQUIRE(shortLayoutFrame.ensureDetectorLayouts(&provider, catalogRequest, ordered,
                                                       TransitionPolicyTag::DiskDisk, oneLayout)
                  .ok());
  try {
    shortLayoutTraits.initialiseTimeFrame(1);
    BOOST_FAIL("iteration beyond the configured layout set must throw");
  } catch (const TraversalException& error) {
    BOOST_CHECK_EQUAL(error.getIteration(), 1);
    BOOST_CHECK(error.getReason() == TraversalFailureReason::IterationOutOfRange);
  }
  BOOST_CHECK(!shortLayoutTraits.hasTraversalCache());
}

BOOST_AUTO_TEST_CASE(traversal_cache_groups_and_binds_once_across_repeated_neighbour_and_road_calls)
{
  auto params = mftTraversalParameters();
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame<10> frame;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, traits, pool, params);
  const auto catalogRequest = request(10, o2::detectors::DetID::MFT);
  FakeCatalogProvider provider{catalog(catalogRequest, SurfaceKind::Disk)};
  const auto ordered = order(10);
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered,
                                            TransitionPolicyTag::DiskDisk, params)
                  .ok());

  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);
  traits.initialiseTimeFrame(0);
  BOOST_REQUIRE(traits.hasTraversalCache());
  BOOST_CHECK_EQUAL(traits.getTraversalGroupingCount(), 1);
  BOOST_CHECK_EQUAL(traits.getPolicyBindingCount(TransitionPolicyTag::CylinderCylinder), 0);
  BOOST_CHECK_EQUAL(traits.getPolicyBindingCount(TransitionPolicyTag::DiskDisk), 1);

  traits.findCellsNeighbours(0);
  traits.findCellsNeighbours(0);
  traits.findRoads(0);
  traits.findRoads(0);
  BOOST_CHECK_EQUAL(traits.getTraversalGroupingCount(), 1);
  BOOST_CHECK_EQUAL(traits.getPolicyBindingCount(TransitionPolicyTag::DiskDisk), 1);

  std::vector<TrackingParameters> itsParams{parameters(7, 0, 0, 0x7f)};
  TimeFrame<7> itsFrame;
  TrackerTraits<7> itsTraits;
  prepareTraversalFrame(itsFrame, itsTraits, pool, itsParams);
  const auto itsCatalogRequest = request(7, o2::detectors::DetID::ITS);
  FakeCatalogProvider itsProvider{catalog(itsCatalogRequest, SurfaceKind::Cylinder)};
  BOOST_REQUIRE(itsFrame.ensureDetectorLayouts(&itsProvider, itsCatalogRequest, order(7),
                                               TransitionPolicyTag::CylinderCylinder, itsParams)
                  .ok());
  itsTraits.setNThreads(1, arena);
  itsTraits.initialiseTimeFrame(0);
  itsTraits.findRoads(0);
  itsTraits.findRoads(0);
  BOOST_CHECK_EQUAL(itsTraits.getTraversalGroupingCount(), 1);
  BOOST_CHECK_EQUAL(itsTraits.getPolicyBindingCount(TransitionPolicyTag::CylinderCylinder), 1);
  BOOST_CHECK_EQUAL(itsTraits.getPolicyBindingCount(TransitionPolicyTag::DiskDisk), 0);
}

BOOST_AUTO_TEST_CASE(traversal_empty_road_start_span_is_valid_and_produces_no_tracks)
{
  // StartLayerMask=0 -> an empty seeding mask -> an empty roadStartCellsForTag
  // span (Architecture.md Sec 10, item 7: "empty road-start span is valid").
  // Unlike testCATrackerFailureContract.cxx's
  // ValidEmptyInputCompletesWithoutErrorAndProducesNoTracks (full StartLayerMask,
  // no cluster data), this exercises the *topologically empty* road-start case
  // through the same initialiseTimeFrame()/findRoads() pair.
  auto params = mftTraversalParameters();
  params[0].StartLayerMask = 0;
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame<10> frame;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, traits, pool, params);
  const auto catalogRequest = request(10, o2::detectors::DetID::MFT);
  FakeCatalogProvider provider{catalog(catalogRequest, SurfaceKind::Disk)};
  const auto ordered = order(10);
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered,
                                            TransitionPolicyTag::DiskDisk, params)
                  .ok());

  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);
  BOOST_CHECK_NO_THROW(traits.initialiseTimeFrame(0));
  BOOST_REQUIRE(traits.hasTraversalCache());
  BOOST_CHECK_NO_THROW(traits.findRoads(0));
  BOOST_CHECK_EQUAL(frame.getNumberOfTracks(), 0u);
}

BOOST_AUTO_TEST_CASE(traversal_legacy_cell_container_size_mismatch_fails_before_indexing)
{
  // Item 4/7: findRoads() indexes mTimeFrame->getCells() with sparse
  // CellTopologyId values; a desync between that legacy container and the
  // cached sparse layout must fail with LegacyIndexMismatch rather than
  // index out of bounds. Reached here through TimeFrame::getCells(), the
  // existing public, non-const production accessor (TimeFrame.h) -- no new
  // mutation API is added for this test.
  auto params = mftTraversalParameters();
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame<10> frame;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, traits, pool, params);
  const auto catalogRequest = request(10, o2::detectors::DetID::MFT);
  FakeCatalogProvider provider{catalog(catalogRequest, SurfaceKind::Disk)};
  const auto ordered = order(10);
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered,
                                            TransitionPolicyTag::DiskDisk, params)
                  .ok());

  std::shared_ptr<tbb::task_arena> arena;
  traits.setNThreads(1, arena);
  traits.initialiseTimeFrame(0);
  BOOST_REQUIRE(traits.hasTraversalCache());
  BOOST_REQUIRE(!frame.getCells().empty());
  frame.getCells().pop_back();

  try {
    traits.findRoads(0);
    BOOST_FAIL("legacy cell-container size mismatch must throw before indexing");
  } catch (const TraversalException& error) {
    BOOST_CHECK(error.getReason() == TraversalFailureReason::LegacyIndexMismatch);
  }
}

BOOST_AUTO_TEST_CASE(traversal_preflight_rejects_legacy_mismatch_state_mismatch_and_bad_parameters)
{
  auto checkFailure = [](std::vector<TrackingParameters> params,
                         std::vector<SurfaceDescriptor> surfaces,
                         std::vector<SurfaceId> ordered,
                         TransitionPolicyTag tag,
                         TraversalFailureReason expected) {
    auto pool = std::make_shared<BoundedMemoryResource>();
    TimeFrame<10> frame;
    TrackerTraits<10> traits;
    prepareTraversalFrame(frame, traits, pool, params);
    const auto catalogRequest = request(10, o2::detectors::DetID::MFT);
    FakeCatalogProvider provider{std::move(surfaces)};
    BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, tag, params).ok());
    try {
      traits.initialiseTimeFrame(0);
      BOOST_FAIL("invalid traversal preflight must throw");
    } catch (const TraversalException& error) {
      BOOST_CHECK(error.getReason() == expected);
    }
    BOOST_CHECK(!traits.hasTraversalCache());
  };

  auto params = mftTraversalParameters();
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT),
               {SurfaceId{3}, SurfaceId{0}, SurfaceId{6}, SurfaceId{2}, SurfaceId{9}, SurfaceId{5}, SurfaceId{1}, SurfaceId{8}, SurfaceId{4}, SurfaceId{7}},
               TransitionPolicyTag::DiskDisk, TraversalFailureReason::LegacyIndexMismatch);
  checkFailure(params, catalog(10, SurfaceKind::Cylinder, o2::detectors::DetID::MFT), order(10),
               TransitionPolicyTag::CylinderCylinder, TraversalFailureReason::StateFamilyMismatch);
  params[0].MaxChi2ClusterAttachment = -1.f;
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk,
               TraversalFailureReason::InvalidPolicyParameters);

  // Perturbing the temporary legacy LayerxX0 away from the catalog's
  // authoritative material now fails the material-compatibility check
  // (LegacyMaterialMismatch) before ever reaching attachHitConfig's own
  // finite/non-negative validation -- the catalog's material is unperturbed
  // and finite, so the mismatch is purely numeric.
  params = mftTraversalParameters();
  params[0].LayerxX0[3] = std::numeric_limits<float>::quiet_NaN();
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk,
               TraversalFailureReason::LegacyMaterialMismatch);

  params = mftTraversalParameters();
  params[0].LayerxX0[7] = std::numeric_limits<float>::infinity();
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk,
               TraversalFailureReason::LegacyMaterialMismatch);

  params = mftTraversalParameters();
  params[0].CorrType = static_cast<o2::base::PropagatorF::MatCorrType>(99);
  checkFailure(params, catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT), order(10), TransitionPolicyTag::DiskDisk,
               TraversalFailureReason::InvalidPolicyParameters);
}

BOOST_AUTO_TEST_CASE(every_iteration_resolves_identical_authoritative_material)
{
  auto params = mftTraversalParameters();
  params.push_back(params.front());
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame<10> frame;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, traits, pool, params);
  const auto catalogRequest = request(10, o2::detectors::DetID::MFT);
  FakeCatalogProvider provider{catalog(catalogRequest, SurfaceKind::Disk)};
  const auto ordered = order(10);
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered,
                                            TransitionPolicyTag::DiskDisk, params)
                  .ok());

  traits.initialiseTimeFrame(0);
  const auto firstIterationMaterial = traits.getLayerMaterial();
  std::vector<NominalSurfaceMaterial> firstIteration(firstIterationMaterial.begin(), firstIterationMaterial.end());

  traits.initialiseTimeFrame(1);
  const auto secondIteration = traits.getLayerMaterial();
  BOOST_REQUIRE_EQUAL(secondIteration.size(), firstIteration.size());
  for (size_t layer = 0; layer < firstIteration.size(); ++layer) {
    BOOST_CHECK_EQUAL(secondIteration[layer].xOverX0, firstIteration[layer].xOverX0);
    BOOST_CHECK_EQUAL(secondIteration[layer].arealDensityGPerCm2, firstIteration[layer].arealDensityGPerCm2);
  }
}

BOOST_AUTO_TEST_CASE(rejected_initialisation_does_not_mutate_surface_descriptor_material)
{
  auto params = mftTraversalParameters();
  params[0].LayerxX0[4] = std::numeric_limits<float>::quiet_NaN();
  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame<10> frame;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, traits, pool, params);
  const auto catalogRequest = request(10, o2::detectors::DetID::MFT);
  FakeCatalogProvider provider{catalog(catalogRequest, SurfaceKind::Disk)};
  const auto ordered = order(10);
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered,
                                            TransitionPolicyTag::DiskDisk, params)
                  .ok());
  const NominalSurfaceMaterial materialBefore = frame.getSurfaceCatalog()->at(4).material;

  try {
    traits.initialiseTimeFrame(0);
    BOOST_FAIL("perturbed LayerxX0 must throw");
  } catch (const TraversalException& error) {
    BOOST_CHECK(error.getReason() == TraversalFailureReason::LegacyMaterialMismatch);
  }

  const auto& materialAfter = frame.getSurfaceCatalog()->at(4).material;
  BOOST_CHECK_EQUAL(materialAfter.xOverX0, materialBefore.xOverX0);
  BOOST_CHECK_EQUAL(materialAfter.arealDensityGPerCm2, materialBefore.arealDensityGPerCm2);
}

BOOST_AUTO_TEST_CASE(non_monotonic_ordered_surfaces_maps_legacy_layers_to_correct_surface_ids)
{
  // Distinct-per-surface material (not the uniform MFT nominal default) so an
  // identity-assuming mapping bug (mLayerMaterial[legacyLayer] read from
  // catalog[legacyLayer] instead of catalog[orderedSurfaces[legacyLayer]])
  // would be observable as a numeric difference, not masked by every layer
  // happening to carry the same value.
  auto surfaces = catalog(10, SurfaceKind::Disk, o2::detectors::DetID::MFT);
  for (uint16_t id = 0; id < surfaces.size(); ++id) {
    const float xOverX0 = 0.001f * static_cast<float>(id + 1);
    surfaces[id].material.xOverX0 = xOverX0;
    surfaces[id].material.arealDensityGPerCm2 = xOverX0 * o2::its::constants::Radl * o2::its::constants::Rho;
  }
  const std::vector<SurfaceId> nonMonotonicOrder{
    SurfaceId{3}, SurfaceId{0}, SurfaceId{6}, SurfaceId{2}, SurfaceId{9},
    SurfaceId{5}, SurfaceId{1}, SurfaceId{8}, SurfaceId{4}, SurfaceId{7}};

  auto params = mftTraversalParameters();
  for (size_t legacyLayer = 0; legacyLayer < nonMonotonicOrder.size(); ++legacyLayer) {
    params[0].LayerxX0[legacyLayer] = surfaces[nonMonotonicOrder[legacyLayer].value()].material.xOverX0;
  }

  auto pool = std::make_shared<BoundedMemoryResource>();
  TimeFrame<10> frame;
  TrackerTraits<10> traits;
  prepareTraversalFrame(frame, traits, pool, params);
  const auto catalogRequest = request(10, o2::detectors::DetID::MFT);
  FakeCatalogProvider provider{surfaces};
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, nonMonotonicOrder,
                                            TransitionPolicyTag::DiskDisk, params)
                  .ok());

  // The non-identity mapping is structurally incompatible with the separate
  // legacy-topology-parity check (validateLegacyParity), so the overall call
  // still fails -- with LegacyIndexMismatch, not a material reason -- but
  // only after mLayerMaterial has already been resolved and committed from
  // the correct (non-identity) mapping, which is what this test verifies.
  try {
    traits.initialiseTimeFrame(0);
    BOOST_FAIL("non-monotonic ordering must fail legacy topology parity");
  } catch (const TraversalException& error) {
    BOOST_CHECK(error.getReason() == TraversalFailureReason::LegacyIndexMismatch);
  }

  const auto resolvedMaterial = traits.getLayerMaterial();
  BOOST_REQUIRE_EQUAL(resolvedMaterial.size(), nonMonotonicOrder.size());
  for (size_t legacyLayer = 0; legacyLayer < nonMonotonicOrder.size(); ++legacyLayer) {
    const auto& expected = surfaces[nonMonotonicOrder[legacyLayer].value()].material;
    BOOST_CHECK_EQUAL(resolvedMaterial[legacyLayer].xOverX0, expected.xOverX0);
    BOOST_CHECK_EQUAL(resolvedMaterial[legacyLayer].arealDensityGPerCm2, expected.arealDensityGPerCm2);
    // Confirms the mapping is genuinely position-driven, not an identity
    // shortcut, at every position where this permutation actually moves the
    // surface (nonMonotonicOrder[legacyLayer] != legacyLayer): the
    // identity-mapped value would differ from `expected` there.
    if (nonMonotonicOrder[legacyLayer].value() != legacyLayer) {
      BOOST_CHECK_NE(resolvedMaterial[legacyLayer].xOverX0, surfaces[legacyLayer].material.xOverX0);
    }
  }
}

BOOST_AUTO_TEST_CASE(traversal_preflight_reports_invalid_schedule_and_mixed_policy_layout)
{
  auto checkInstalledLayout = [](BuiltLayout layout, TraversalFailureReason expected) {
    auto params = mftTraversalParameters();
    auto pool = std::make_shared<BoundedMemoryResource>();
    TraversalTestTimeFrame frame;
    TrackerTraits<10> traits;
    prepareTraversalFrame(frame, traits, pool, params);
    frame.installLayout(std::move(layout));
    try {
      traits.initialiseTimeFrame(0);
      BOOST_FAIL("invalid installed layout must throw");
    } catch (const TraversalException& error) {
      BOOST_CHECK(error.getReason() == expected);
    }
    BOOST_CHECK(!traits.hasTraversalCache());
  };

  checkInstalledLayout(cyclicDiskLayout(), TraversalFailureReason::InvalidTraversalSchedule);
  checkInstalledLayout(mixedDisconnectedLayout(), TraversalFailureReason::MixedPolicyLayout);
}

BOOST_AUTO_TEST_CASE(malformed_catalog_failure_preserves_stale_owner)
{
  TimeFrame<7> frame;
  const auto catalogRequest = request(7);
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{parameters(7)};
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params).ok());
  frame.invalidateDetectorLayouts();
  provider.mCatalog[3].id = SurfaceId{6};
  const auto failed = frame.ensureDetectorLayouts(&provider, catalogRequest, ordered, TransitionPolicyTag::CylinderCylinder, params);
  BOOST_CHECK(failed.error == DetectorLayoutSetBuildError::InvalidCatalog);
  BOOST_CHECK(failed.catalogValidationError == DetectorSurfaceCatalogValidationError::NonDenseGlobalSurfaceIds);
  BOOST_CHECK(frame.hasStoredDetectorLayouts());
  BOOST_CHECK(!frame.detectorLayoutsCurrent());
}

BOOST_AUTO_TEST_CASE(its_legacy_topology_and_road_start_parity)
{
  checkLegacyParity<7>(SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder,
                       (uint16_t{1} << 6) | (uint16_t{1} << 3));
}

BOOST_AUTO_TEST_CASE(mft_legacy_topology_and_road_start_parity)
{
  checkLegacyParity<10>(SurfaceKind::Disk, TransitionPolicyTag::DiskDisk,
                        (uint16_t{1} << 9) | (uint16_t{1} << 5));
}
