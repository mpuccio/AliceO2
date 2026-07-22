// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE ITSMFT DetectorLayoutView
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <type_traits>
#include <vector>

#include "ITSMFTTracking/DetectorLayoutSet.h"

using namespace o2::itsmft::tracking;

namespace
{
SurfaceDescriptor surface(uint16_t id, SurfaceKind kind = SurfaceKind::Cylinder)
{
  return SurfaceDescriptor{SurfaceId{id}, id, 0, kind};
}

DetectorLayoutSet buildSet(std::vector<SurfaceDescriptor> catalog, std::vector<NominalSurfaceMaterialBudget> material, int nIterations)
{
  DetectorLayoutConfigurationKey key;
  key.geometryEpoch = 1;
  std::vector<DetectorLayout> layouts;
  layouts.reserve(nIterations);
  for (int i = 0; i < nIterations; ++i) {
    SparseTrackingTopology topology{static_cast<uint32_t>(catalog.size())};
    topology.finalize();
    layouts.emplace_back(catalog, std::move(topology));
  }
  return DetectorLayoutSet{std::move(key), std::move(catalog), std::move(material), std::move(layouts)};
}

// Structural proof that no per-iteration surface ownership remains:
// DetectorLayout no longer exposes getSurfaces()/getCylinderSurfaces()/
// getDiskSurfaces() at all (compile-time checked, not just "unused").
template <typename T>
concept HasGetSurfaces = requires(const T& t) { t.getSurfaces(); };
template <typename T>
concept HasGetCylinderSurfaces = requires(const T& t) { t.getCylinderSurfaces(); };
template <typename T>
concept HasGetDiskSurfaces = requires(const T& t) { t.getDiskSurfaces(); };
static_assert(!HasGetSurfaces<DetectorLayout>);
static_assert(!HasGetCylinderSurfaces<DetectorLayout>);
static_assert(!HasGetDiskSurfaces<DetectorLayout>);
} // namespace

// -------------------------------------------------------------------------
// POD / ABI lock (compiler-verified, see DetectorLayout.h)
// -------------------------------------------------------------------------

static_assert(std::is_standard_layout_v<DetectorLayoutView>);
static_assert(std::is_trivially_copyable_v<DetectorLayoutView>);
static_assert(sizeof(DetectorLayoutView) == 88);
static_assert(alignof(DetectorLayoutView) == 8);
static_assert(offsetof(DetectorLayoutView, surfaces) == 0);
static_assert(offsetof(DetectorLayoutView, nSurfaces) == 8);
static_assert(offsetof(DetectorLayoutView, cylinderSurfaces) == 12);
static_assert(offsetof(DetectorLayoutView, diskSurfaces) == 16);
static_assert(offsetof(DetectorLayoutView, topology) == 24);
static_assert(offsetof(DetectorLayoutView, nominalMaterial) == 72);
static_assert(offsetof(DetectorLayoutView, status) == 80);
static_assert(static_cast<uint8_t>(DetectorLayoutViewStatus::Invalid) == 0);
static_assert(static_cast<uint8_t>(DetectorLayoutViewStatus::Valid) == 1);

// -------------------------------------------------------------------------
// Valid / invalid status contract
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(DefaultConstructedViewIsInvalid)
{
  DetectorLayoutView view{};
  BOOST_CHECK(!view.isValid());
  BOOST_CHECK(view.status == DetectorLayoutViewStatus::Invalid);
  BOOST_CHECK(view.surfaces == nullptr);
  BOOST_CHECK(view.nominalMaterial == nullptr);
}

BOOST_AUTO_TEST_CASE(OutOfRangeIterationReturnsInvalid)
{
  auto set = buildSet({surface(0), surface(1)}, {{}, {}}, 1);
  BOOST_CHECK(!set.getLayoutView(5).isValid());
  BOOST_CHECK(set.getLayoutView(5).status == DetectorLayoutViewStatus::Invalid);
}

BOOST_AUTO_TEST_CASE(ZeroIterationsHasNoValidView)
{
  auto set = buildSet({surface(0)}, {{}}, 0);
  BOOST_CHECK_EQUAL(set.size(), 0u);
  BOOST_CHECK(!set.getLayoutView(0).isValid());
}

BOOST_AUTO_TEST_CASE(ValidIterationIsValidIndependentOfZeroCounts)
{
  // Zero surfaces, one (trivially valid) iteration: a legitimate
  // empty-but-valid layout, distinguished from the invalid default purely
  // via status, not via nSurfaces/pointer nullness.
  auto set = buildSet({}, {}, 1);
  const auto view = set.getLayoutView(0);
  BOOST_CHECK(view.isValid());
  BOOST_CHECK(view.status == DetectorLayoutViewStatus::Valid);
  BOOST_CHECK_EQUAL(view.nSurfaces, 0u);
}

// -------------------------------------------------------------------------
// Shared storage: identical pointers across multiple iteration views
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(MultipleIterationViewsShareTheSameSurfaceAndMaterialPointers)
{
  std::vector<SurfaceDescriptor> catalog{surface(0, SurfaceKind::Cylinder), surface(1, SurfaceKind::Cylinder), surface(2, SurfaceKind::Disk)};
  std::vector<NominalSurfaceMaterialBudget> material{{0.01f, 0.f}, {0.02f, 0.f}, {0.03f, 0.5f}};
  auto set = buildSet(catalog, material, 3);

  const auto v0 = set.getLayoutView(0);
  const auto v1 = set.getLayoutView(1);
  const auto v2 = set.getLayoutView(2);
  BOOST_REQUIRE(v0.isValid() && v1.isValid() && v2.isValid());

  // Same shared arrays, not per-iteration copies: identical pointers.
  BOOST_CHECK(v0.surfaces == v1.surfaces);
  BOOST_CHECK(v1.surfaces == v2.surfaces);
  BOOST_CHECK(v0.nominalMaterial == v1.nominalMaterial);
  BOOST_CHECK(v1.nominalMaterial == v2.nominalMaterial);
  BOOST_CHECK_EQUAL(v0.nSurfaces, 3u);
  BOOST_CHECK_EQUAL(v0.cylinderSurfaces.value(), v1.cylinderSurfaces.value());
  BOOST_CHECK_EQUAL(v0.diskSurfaces.value(), v1.diskSurfaces.value());
  BOOST_CHECK_EQUAL(v0.cylinderSurfaces.value(), 0x3u);
  BOOST_CHECK_EQUAL(v0.diskSurfaces.value(), 0x4u);

  BOOST_CHECK_EQUAL(v0.getNominalMaterial(SurfaceId{2}).normalArealDensityGPerCm2, 0.5f);
}

BOOST_AUTO_TEST_CASE(SharedCatalogAndMasksAreComputedOnceFromDetectorLayoutSet)
{
  std::vector<SurfaceDescriptor> catalog{surface(0, SurfaceKind::Disk), surface(1, SurfaceKind::Cylinder)};
  auto set = buildSet(catalog, {{}, {}}, 2);
  BOOST_CHECK_EQUAL(set.getCylinderSurfaces().value(), 0x2u);
  BOOST_CHECK_EQUAL(set.getDiskSurfaces().value(), 0x1u);
  BOOST_CHECK_EQUAL(&set.getSurfaceCatalog(), &set.getSurfaceCatalog()); // same object, not recomputed
}

// -------------------------------------------------------------------------
// DetectorLayout::getView() alignment requirement (nominalMaterial.size()
// must equal surfaces.size(); no default argument -- see DetectorLayout.h)
// -------------------------------------------------------------------------

namespace
{
DetectorLayout buildLayout(const std::vector<SurfaceDescriptor>& surfaces)
{
  SparseTrackingTopology topology{static_cast<uint32_t>(surfaces.size())};
  BOOST_REQUIRE(topology.finalize());
  return DetectorLayout{surfaces, std::move(topology)};
}
} // namespace

BOOST_AUTO_TEST_CASE(GetViewValidForMatchingNonEmptyGeometryAndMaterial)
{
  const std::vector<SurfaceDescriptor> surfaces{surface(0), surface(1)};
  auto layout = buildLayout(surfaces);
  BOOST_REQUIRE(layout.valid());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const std::vector<NominalSurfaceMaterialBudget> material{{0.01f, 0.f}, {0.02f, 0.5f}};
  const auto view = layout.getView(surfaces, masks.first, masks.second, material);
  BOOST_CHECK(view.isValid());
  BOOST_CHECK_EQUAL(view.nSurfaces, 2u);
  BOOST_CHECK_EQUAL(view.getNominalMaterial(SurfaceId{1}).normalArealDensityGPerCm2, 0.5f);
}

BOOST_AUTO_TEST_CASE(GetViewInvalidForNonEmptyGeometryWithEmptyMaterial)
{
  const std::vector<SurfaceDescriptor> surfaces{surface(0), surface(1)};
  auto layout = buildLayout(surfaces);
  BOOST_REQUIRE(layout.valid());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = layout.getView(surfaces, masks.first, masks.second, {});
  BOOST_CHECK(!view.isValid());
  BOOST_CHECK(view.status == DetectorLayoutViewStatus::Invalid);
}

BOOST_AUTO_TEST_CASE(GetViewInvalidForMismatchedNonEmptySizes)
{
  const std::vector<SurfaceDescriptor> surfaces{surface(0), surface(1)};
  auto layout = buildLayout(surfaces);
  BOOST_REQUIRE(layout.valid());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const std::vector<NominalSurfaceMaterialBudget> tooFewMaterial{{0.01f, 0.f}}; // 1 entry for 2 surfaces
  const auto view = layout.getView(surfaces, masks.first, masks.second, tooFewMaterial);
  BOOST_CHECK(!view.isValid());
}

BOOST_AUTO_TEST_CASE(GetViewValidForLegitimateEmptyGeometryAndEmptyMaterial)
{
  auto layout = buildLayout({});
  BOOST_REQUIRE(layout.valid());
  const auto view = layout.getView({}, {}, {}, {});
  BOOST_CHECK(view.isValid());
  BOOST_CHECK_EQUAL(view.nSurfaces, 0u);
}

BOOST_AUTO_TEST_CASE(DetectorLayoutSetProducedViewsRemainValid)
{
  std::vector<SurfaceDescriptor> catalog{surface(0), surface(1), surface(2)};
  std::vector<NominalSurfaceMaterialBudget> material{{0.01f, 0.f}, {0.02f, 0.f}, {0.03f, 0.f}};
  auto set = buildSet(catalog, material, 2);
  BOOST_CHECK(set.getLayoutView(0).isValid());
  BOOST_CHECK(set.getLayoutView(1).isValid());
}
