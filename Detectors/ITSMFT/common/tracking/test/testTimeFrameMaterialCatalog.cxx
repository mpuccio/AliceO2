// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT TimeFrame material catalog
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <limits>
#include <vector>

#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"
#include "ITSMFTTracking/TimeFrame.h"

using namespace o2::itsmft::tracking;
using o2::itsmft::TrackingParameters;

namespace
{
class FakeCatalogProvider final : public DetectorSurfaceCatalogProvider
{
 public:
  explicit FakeCatalogProvider(std::vector<SurfaceDescriptor> catalog) : mCatalog{std::move(catalog)} {}

  DetectorSurfaceCatalogResult buildCatalog(const DetectorSurfaceCatalogRequest&) const final
  {
    return DetectorSurfaceCatalogResult{mCatalog, DetectorSurfaceCatalogError::None};
  }

  std::vector<SurfaceDescriptor> mCatalog;
};

DetectorSurfaceCatalogRequest request(uint32_t count)
{
  return DetectorSurfaceCatalogRequest{o2::detectors::DetID::ITS, SurfaceId{0}, count};
}

std::vector<SurfaceDescriptor> catalog(size_t count)
{
  std::vector<SurfaceDescriptor> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder, 0, static_cast<float>(i + 1), 0.f, 100.f});
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

std::vector<NominalSurfaceMaterialEntry> matchingEntries(const std::vector<SurfaceDescriptor>& c, float xOverX0 = 0.01f, float arealDensity = 0.02f)
{
  std::vector<NominalSurfaceMaterialEntry> entries;
  entries.reserve(c.size());
  for (const auto& s : c) {
    entries.push_back(NominalSurfaceMaterialEntry{s.id, NominalSurfaceMaterialBudget{xOverX0, arealDensity}});
  }
  return entries;
}

struct MaterialTestTimeFrame : TimeFrame<7> {
  void setRequiredGeometryEpoch(DetectorGeometryEpoch epoch) { mRequiredDetectorGeometryEpoch = epoch; }
  void setRequiredMaterialEpoch(MaterialCatalogEpoch epoch) { mRequiredMaterialCatalogEpoch = epoch; }
  const DetectorLayoutSet* getStoredDetectorLayouts() const { return mDetectorLayouts ? &*mDetectorLayouts : nullptr; }
};
} // namespace

// -------------------------------------------------------------------------
// Compatibility overload (no materialEntries parameter at all): auto-filled,
// zero, still validated. This is the temporary compatibility boundary --
// see ensureDetectorLayouts()'s declaration comment in TimeFrame.h.
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(CompatibilityOverloadAutoFillsZeroMaterialAndSucceeds)
{
  TimeFrame<7> frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  const auto result = frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params);
  BOOST_REQUIRE(result.ok());
  const auto* layouts = frame.getDetectorLayouts();
  BOOST_REQUIRE(layouts != nullptr);
  BOOST_REQUIRE_EQUAL(layouts->getNominalMaterial().size(), 7u);
  for (const auto& budget : layouts->getNominalMaterial()) {
    BOOST_CHECK_EQUAL(budget.normalXOverX0, 0.f);
    BOOST_CHECK_EQUAL(budget.normalArealDensityGPerCm2, 0.f);
  }
}

// -------------------------------------------------------------------------
// Explicit overload: validates exactly what is supplied, no auto-fill.
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ExplicitOverloadWithEmptySpanAndNonEmptyCatalogueFailsSizeMismatch)
{
  TimeFrame<7> frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  // Explicitly empty, not omitted: this must NOT auto-fill (that is the
  // compatibility overload's job), so it fails SizeMismatch against the
  // resolved 7-surface catalog.
  const gsl::span<const NominalSurfaceMaterialEntry> emptyEntries{};
  const auto result = frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, emptyEntries);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutSetBuildError::InvalidMaterial);
  BOOST_CHECK(result.materialValidationError == DetectorLayoutMaterialValidationError::SizeMismatch);
  BOOST_CHECK(!frame.hasStoredDetectorLayouts());
}

// A zero-count DetectorSurfaceCatalogRequest is already rejected by the
// pre-existing, unrelated geometry validation (EmptyDetector, checked
// before material validation ever runs) -- so "empty surface catalogue"
// is not a state ensureDetectorLayouts() can reach via its public,
// request-based provider path. The material-validation contract for a
// genuinely empty geometry/material pair is exercised directly at the
// DetectorLayout::getView() layer instead:
// testDetectorLayoutView.cxx::GetViewValidForLegitimateEmptyGeometryAndEmptyMaterial.

BOOST_AUTO_TEST_CASE(ExplicitOverloadWithCorrectlySizedEntriesSucceeds)
{
  TimeFrame<7> frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  const auto entries = matchingEntries(catalog(7), 0.03f, 0.04f);
  const auto result = frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries);
  BOOST_REQUIRE(result.ok());
  const auto& stored = frame.getDetectorLayouts()->getNominalMaterial();
  BOOST_REQUIRE_EQUAL(stored.size(), 7u);
  for (const auto& budget : stored) {
    BOOST_CHECK_EQUAL(budget.normalXOverX0, 0.03f);
    BOOST_CHECK_EQUAL(budget.normalArealDensityGPerCm2, 0.04f);
  }
}

BOOST_AUTO_TEST_CASE(ExplicitOverloadEmptySpanFailureIsTransactional)
{
  // A prior successful build (via the compatibility overload) must survive,
  // stale-but-stored, after a failed explicit-overload call.
  TimeFrame<7> frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params).ok());
  BOOST_REQUIRE(frame.hasStoredDetectorLayouts());

  frame.invalidateNominalMaterial();
  const gsl::span<const NominalSurfaceMaterialEntry> emptyEntries{};
  const auto failed = frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, emptyEntries);
  BOOST_CHECK(!failed.ok());
  BOOST_CHECK(frame.hasStoredDetectorLayouts()); // old owner retained
  BOOST_CHECK(!frame.detectorLayoutsCurrent());  // but exposed as stale
}

// -------------------------------------------------------------------------
// Transactional validation precedence
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(MaterialSizeMismatchRejected)
{
  TimeFrame<7> frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  auto entries = matchingEntries(catalog(7));
  entries.pop_back();
  const auto result = frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutSetBuildError::InvalidMaterial);
  BOOST_CHECK(result.materialValidationError == DetectorLayoutMaterialValidationError::SizeMismatch);
}

BOOST_AUTO_TEST_CASE(ReorderedSurfaceIdRejected)
{
  TimeFrame<7> frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  auto entries = matchingEntries(catalog(7));
  std::swap(entries[0], entries[1]);
  const auto result = frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error == DetectorLayoutSetBuildError::InvalidMaterial);
  BOOST_CHECK(result.materialValidationError == DetectorLayoutMaterialValidationError::SurfaceIdMismatch);
}

BOOST_AUTO_TEST_CASE(DuplicateSurfaceIdRejected)
{
  TimeFrame<7> frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  auto entries = matchingEntries(catalog(7));
  entries[1].surface = entries[0].surface;
  const auto result = frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.materialValidationError == DetectorLayoutMaterialValidationError::SurfaceIdMismatch);
}

BOOST_AUTO_TEST_CASE(NegativeNaNInfiniteXOverX0RejectedIndependently)
{
  TimeFrame<7> frame;
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  for (const float bad : {-0.1f, std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity()}) {
    FakeCatalogProvider provider{catalog(7)};
    auto entries = matchingEntries(catalog(7));
    entries[3].budget.normalXOverX0 = bad;
    const auto result = frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries);
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.materialValidationError == DetectorLayoutMaterialValidationError::InvalidBudget);
    frame.invalidateDetectorLayouts();
  }
}

BOOST_AUTO_TEST_CASE(NegativeNaNInfiniteArealDensityRejectedIndependently)
{
  TimeFrame<7> frame;
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  for (const float bad : {-0.1f, std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity()}) {
    FakeCatalogProvider provider{catalog(7)};
    auto entries = matchingEntries(catalog(7));
    entries[3].budget.normalArealDensityGPerCm2 = bad;
    const auto result = frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries);
    BOOST_CHECK(!result.ok());
    BOOST_CHECK(result.materialValidationError == DetectorLayoutMaterialValidationError::InvalidBudget);
    frame.invalidateDetectorLayouts();
  }
}

BOOST_AUTO_TEST_CASE(IndependentlyZeroXOverX0AndArealDensityAccepted)
{
  TimeFrame<7> frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  auto entries = matchingEntries(catalog(7), 0.f, 0.03f); // zero xOverX0, nonzero areal density
  entries[2].budget = NominalSurfaceMaterialBudget{0.05f, 0.f}; // nonzero xOverX0, zero areal density
  const auto result = frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries);
  BOOST_REQUIRE(result.ok());
  const auto& stored = frame.getDetectorLayouts()->getNominalMaterial();
  BOOST_CHECK_EQUAL(stored[0].normalXOverX0, 0.f);
  BOOST_CHECK_EQUAL(stored[0].normalArealDensityGPerCm2, 0.03f);
  BOOST_CHECK_EQUAL(stored[2].normalXOverX0, 0.05f);
  BOOST_CHECK_EQUAL(stored[2].normalArealDensityGPerCm2, 0.f);
}

// -------------------------------------------------------------------------
// Epoch semantics
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(GeometryOnlyEpochChangeInvalidatesIndependentlyOfMaterial)
{
  MaterialTestTimeFrame frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  auto entries = matchingEntries(catalog(7));
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries).ok());
  const auto materialEpochBefore = frame.getRequiredMaterialCatalogEpoch();
  frame.invalidateDetectorLayouts(); // geometry only
  BOOST_CHECK(!frame.detectorLayoutsCurrent());
  BOOST_CHECK_EQUAL(frame.getRequiredMaterialCatalogEpoch(), materialEpochBefore); // unchanged
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries).ok());
  BOOST_CHECK(frame.detectorLayoutsCurrent());
}

BOOST_AUTO_TEST_CASE(MaterialOnlyEpochChangeInvalidatesIndependentlyOfGeometry)
{
  MaterialTestTimeFrame frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  auto entries = matchingEntries(catalog(7));
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries).ok());
  const auto geometryEpochBefore = frame.getRequiredDetectorGeometryEpoch();
  frame.invalidateNominalMaterial(); // material only
  BOOST_CHECK(!frame.detectorLayoutsCurrent());
  BOOST_CHECK_EQUAL(frame.getRequiredDetectorGeometryEpoch(), geometryEpochBefore); // unchanged
  entries[0].budget.normalXOverX0 = 0.5f; // documented contract: content changed alongside the epoch bump
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries).ok());
  BOOST_CHECK(frame.detectorLayoutsCurrent());
  BOOST_CHECK_EQUAL(frame.getDetectorLayouts()->getNominalMaterial()[0].normalXOverX0, 0.5f);
}

BOOST_AUTO_TEST_CASE(MaterialEpochWrapDropsStoredOwner)
{
  MaterialTestTimeFrame frame;
  frame.setRequiredMaterialEpoch(std::numeric_limits<MaterialCatalogEpoch>::max());
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  auto entries = matchingEntries(catalog(7));
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries).ok());
  BOOST_REQUIRE(frame.hasStoredDetectorLayouts());
  frame.invalidateNominalMaterial();
  BOOST_CHECK_EQUAL(frame.getRequiredMaterialCatalogEpoch(), InitialMaterialCatalogEpoch);
  BOOST_CHECK(!frame.hasStoredDetectorLayouts());
  BOOST_CHECK(!frame.detectorLayoutsCurrent());
}

// -------------------------------------------------------------------------
// Transactional failed rebuild
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(FailedMaterialRebuildLeavesOldOwnerStoredButStale)
{
  MaterialTestTimeFrame frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  auto goodEntries = matchingEntries(catalog(7));
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, goodEntries).ok());
  const auto* storedOwner = frame.getStoredDetectorLayouts();
  BOOST_REQUIRE(storedOwner != nullptr);
  const auto storedMaterialEpoch = storedOwner->getConfigurationKey().materialEpoch;

  frame.invalidateNominalMaterial();
  auto badEntries = matchingEntries(catalog(7));
  badEntries[0].budget.normalXOverX0 = -1.f;
  const auto failed = frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, badEntries);
  BOOST_CHECK(!failed.ok());
  BOOST_CHECK(failed.error == DetectorLayoutSetBuildError::InvalidMaterial);
  // Old owner retained (not reset), but now exposed as stale.
  BOOST_CHECK(frame.getStoredDetectorLayouts() == storedOwner);
  BOOST_CHECK_EQUAL(frame.getStoredDetectorLayouts()->getConfigurationKey().materialEpoch, storedMaterialEpoch);
  BOOST_CHECK(frame.hasStoredDetectorLayouts());
  BOOST_CHECK(!frame.detectorLayoutsCurrent());
  BOOST_CHECK(frame.getDetectorLayouts() == nullptr); // current-only accessor hides the stale owner
}

// -------------------------------------------------------------------------
// wipe() preservation
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(WipePreservesLayoutMaterialOwnershipAndBothEpochs)
{
  TimeFrame<7> frame;
  FakeCatalogProvider provider{catalog(7)};
  auto ordered = order(7);
  std::vector<TrackingParameters> params{};
  auto entries = matchingEntries(catalog(7));
  BOOST_REQUIRE(frame.ensureDetectorLayouts(&provider, request(7), ordered, TransitionPolicyTag::CylinderCylinder, params, entries).ok());
  const auto geometryEpochBefore = frame.getRequiredDetectorGeometryEpoch();
  const auto materialEpochBefore = frame.getRequiredMaterialCatalogEpoch();
  BOOST_REQUIRE(frame.detectorLayoutsCurrent());

  frame.wipe();

  BOOST_CHECK(frame.detectorLayoutsCurrent());
  BOOST_CHECK_EQUAL(frame.getRequiredDetectorGeometryEpoch(), geometryEpochBefore);
  BOOST_CHECK_EQUAL(frame.getRequiredMaterialCatalogEpoch(), materialEpochBefore);
  BOOST_REQUIRE(frame.getDetectorLayouts() != nullptr);
  BOOST_CHECK_EQUAL(frame.getDetectorLayouts()->getNominalMaterial().size(), 7u);
}
