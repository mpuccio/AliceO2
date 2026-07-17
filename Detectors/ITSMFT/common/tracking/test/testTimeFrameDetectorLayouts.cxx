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

#include <algorithm>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingInterface.h"

using namespace o2::itsmft::tracking;
using o2::itsmft::TrackingParameters;

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

static_assert(std::is_nothrow_move_constructible_v<DetectorLayoutSet>);

DetectorSurfaceCatalogRequest request(uint32_t count,
                                      o2::detectors::DetID::ID detector = o2::detectors::DetID::ITS,
                                      uint16_t firstSurface = 0)
{
  return DetectorSurfaceCatalogRequest{detector, SurfaceId{firstSurface}, count};
}

std::vector<SurfaceDescriptor> catalog(size_t count, SurfaceKind kind = SurfaceKind::Cylinder,
                                       o2::detectors::DetID::ID detector = o2::detectors::DetID::ITS)
{
  std::vector<SurfaceDescriptor> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(detector), kind, 0, static_cast<float>(i + 1), 0.f, 100.f});
  }
  return result;
}

std::vector<SurfaceDescriptor> catalog(const DetectorSurfaceCatalogRequest& catalogRequest,
                                       SurfaceKind kind = SurfaceKind::Cylinder)
{
  const auto size = catalogRequest.firstSurface.value() + catalogRequest.detectorSurfaceCount;
  auto result = catalog(size, kind, o2::detectors::DetID::ITS);
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
  const auto sparse = frame.getDetectorLayoutView(0).topology;
  BOOST_REQUIRE_EQUAL(sparse.nTransitions, legacyView.nTransitions);
  BOOST_REQUIRE_EQUAL(sparse.nCells, legacyView.nCells);
  for (uint32_t i = 0; i < sparse.nTransitions; ++i) {
    BOOST_CHECK_EQUAL(sparse.transitions[i].from.value(), legacyView.transitions[i].fromLayer);
    BOOST_CHECK_EQUAL(sparse.transitions[i].to.value(), legacyView.transitions[i].toLayer);
  }
  for (uint32_t i = 0; i < sparse.nCells; ++i) {
    const auto& sparseCell = sparse.cells[i];
    const auto& legacyCell = legacyView.cells[i];
    BOOST_CHECK_EQUAL(sparseCell.firstTransition.value(), legacyCell.firstTransition);
    BOOST_CHECK_EQUAL(sparseCell.secondTransition.value(), legacyCell.secondTransition);
    BOOST_CHECK_EQUAL(sparseCell.hitSurfaces.value(), static_cast<uint16_t>(legacyCell.hitLayerMask));
    const auto legacyStartLayer = legacyCell.hitLayerMask.last();
    const auto sparseStartSurface = SurfaceId{static_cast<uint16_t>(sparseCell.hitSurfaces.last())};
    BOOST_CHECK_EQUAL(sparse.seedingSurfaces.has(sparseStartSurface), params[0].StartLayerMask.has(legacyStartLayer));
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
  BOOST_CHECK_EQUAL(frame.getDetectorLayout(0)->getSurfaces()[0].referenceCoordinate, 1.f);
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
  BOOST_CHECK_EQUAL(frame.getDetectorLayout(0)->getSurfaces()[0].referenceCoordinate, 42.f);
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
  BOOST_REQUIRE_EQUAL(full->getSurfaces().size(), 7);
  BOOST_REQUIRE_EQUAL(reduced->getSurfaces().size(), 7);
  for (size_t i = 0; i < 7; ++i) {
    const auto& lhs = full->getSurfaces()[i];
    const auto& rhs = reduced->getSurfaces()[i];
    BOOST_CHECK_EQUAL(lhs.id.value(), rhs.id.value());
    BOOST_CHECK_EQUAL(lhs.detectorSurfaceIndex, rhs.detectorSurfaceIndex);
    BOOST_CHECK_EQUAL(lhs.detectorId, rhs.detectorId);
    BOOST_CHECK(lhs.kind == rhs.kind);
    BOOST_CHECK_EQUAL(lhs.flags, rhs.flags);
    BOOST_CHECK_EQUAL(lhs.referenceCoordinate, rhs.referenceCoordinate);
    BOOST_CHECK_EQUAL(lhs.radialMin, rhs.radialMin);
    BOOST_CHECK_EQUAL(lhs.radialMax, rhs.radialMax);
  }
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
