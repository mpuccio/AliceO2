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
SurfaceDescriptor surface(uint16_t id, SurfaceKind kind = SurfaceKind::Cylinder, NominalSurfaceMaterial material = {})
{
  SurfaceDescriptor result{SurfaceId{id}, id, 0, kind};
  result.material = material;
  return result;
}

// DetectorLayoutSet borrows (SurfaceCatalogView): `catalog` is the caller's
// own named local, which must outlive the returned DetectorLayoutSet -- this
// helper never takes ownership of it.
DetectorLayoutSet buildSet(const std::vector<SurfaceDescriptor>& catalog, int nIterations)
{
  DetectorLayoutConfigurationKey key;
  const SurfaceCatalogView view{catalog.data(), static_cast<uint32_t>(catalog.size())};
  std::vector<DetectorLayout> layouts;
  layouts.reserve(nIterations);
  for (int i = 0; i < nIterations; ++i) {
    SparseTrackingTopology topology{static_cast<uint32_t>(catalog.size())};
    topology.finalize();
    layouts.emplace_back(catalog, std::move(topology));
  }
  return DetectorLayoutSet{std::move(key), view, std::move(layouts)};
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
static_assert(sizeof(DetectorLayoutView) == 72);
static_assert(alignof(DetectorLayoutView) == 8);
static_assert(offsetof(DetectorLayoutView, surfaces) == 0);
static_assert(offsetof(DetectorLayoutView, nSurfaces) == 8);
static_assert(offsetof(DetectorLayoutView, cylinderSurfaces) == 12);
static_assert(offsetof(DetectorLayoutView, diskSurfaces) == 16);
static_assert(offsetof(DetectorLayoutView, topology) == 24);

// -------------------------------------------------------------------------
// Default / out-of-range sentinel
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(DefaultConstructedViewIsTheSentinel)
{
  DetectorLayoutView view{};
  BOOST_CHECK(view.surfaces == nullptr);
  BOOST_CHECK_EQUAL(view.nSurfaces, 0u);
}

BOOST_AUTO_TEST_CASE(OutOfRangeIterationReturnsSentinel)
{
  const std::vector<SurfaceDescriptor> catalog{surface(0), surface(1)};
  auto set = buildSet(catalog, 1);
  BOOST_CHECK(set.getLayoutView(5).surfaces == nullptr);
}

// -------------------------------------------------------------------------
// Shared storage: identical pointers across multiple iteration views, and
// material values visible identically from every iteration view.
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(MultipleIterationViewsShareTheSameSurfaceCatalog)
{
  std::vector<SurfaceDescriptor> catalog{surface(0, SurfaceKind::Cylinder, {0.01f, 0.f}),
                                         surface(1, SurfaceKind::Cylinder, {0.02f, 0.f}),
                                         surface(2, SurfaceKind::Disk, {0.03f, 0.5f})};
  auto set = buildSet(catalog, 3);

  const auto v0 = set.getLayoutView(0);
  const auto v1 = set.getLayoutView(1);
  const auto v2 = set.getLayoutView(2);
  BOOST_REQUIRE(v0.surfaces != nullptr && v1.surfaces != nullptr && v2.surfaces != nullptr);

  // Same shared catalog, not per-iteration copies: identical pointers.
  BOOST_CHECK(v0.surfaces == v1.surfaces);
  BOOST_CHECK(v1.surfaces == v2.surfaces);
  BOOST_CHECK_EQUAL(v0.nSurfaces, 3u);
  BOOST_CHECK_EQUAL(v0.cylinderSurfaces.value(), v1.cylinderSurfaces.value());
  BOOST_CHECK_EQUAL(v0.diskSurfaces.value(), v1.diskSurfaces.value());
  BOOST_CHECK_EQUAL(v0.cylinderSurfaces.value(), 0x3u);
  BOOST_CHECK_EQUAL(v0.diskSurfaces.value(), 0x4u);

  // Material is reached through the shared SurfaceDescriptor, not a
  // parallel array -- and is identical from every iteration's view.
  BOOST_CHECK_EQUAL(v0.getSurface(SurfaceId{2}).material.arealDensityGPerCm2, 0.5f);
  BOOST_CHECK_EQUAL(v0.getSurface(SurfaceId{2}).material.arealDensityGPerCm2,
                    v1.getSurface(SurfaceId{2}).material.arealDensityGPerCm2);
  BOOST_CHECK_EQUAL(v1.getSurface(SurfaceId{2}).material.arealDensityGPerCm2,
                    v2.getSurface(SurfaceId{2}).material.arealDensityGPerCm2);
  BOOST_CHECK(&v0.getSurface(SurfaceId{0}) == &v1.getSurface(SurfaceId{0}));
}

BOOST_AUTO_TEST_CASE(SharedCatalogAndMasksAreComputedOnceFromDetectorLayoutSet)
{
  std::vector<SurfaceDescriptor> catalog{surface(0, SurfaceKind::Disk), surface(1, SurfaceKind::Cylinder)};
  auto set = buildSet(catalog, 2);
  BOOST_CHECK_EQUAL(set.getCylinderSurfaces().value(), 0x2u);
  BOOST_CHECK_EQUAL(set.getDiskSurfaces().value(), 0x1u);
  // getSurfaceCatalog() returns a SurfaceCatalogView by value now (Gate 4
  // B2 Slice 2: a borrowed view, not an owned container), so "not
  // recomputed" is checked via the borrowed pointer's identity instead of
  // the accessor's own (nonexistent) object identity.
  BOOST_CHECK(set.getSurfaceCatalog().surfaces == set.getSurfaceCatalog().surfaces);
}

// -------------------------------------------------------------------------
// DetectorLayout::getView()
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

BOOST_AUTO_TEST_CASE(GetViewValidForNonEmptyGeometry)
{
  const std::vector<SurfaceDescriptor> surfaces{surface(0, SurfaceKind::Cylinder, {0.01f, 0.f}),
                                                surface(1, SurfaceKind::Cylinder, {0.02f, 0.5f})};
  auto layout = buildLayout(surfaces);
  BOOST_REQUIRE(layout.valid());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = layout.getView(surfaces, masks.first, masks.second);
  BOOST_CHECK(view.surfaces != nullptr);
  BOOST_CHECK_EQUAL(view.nSurfaces, 2u);
  BOOST_CHECK_EQUAL(view.getSurface(SurfaceId{1}).material.arealDensityGPerCm2, 0.5f);
}

BOOST_AUTO_TEST_CASE(GetViewIsSentinelForInvalidLayout)
{
  // Mismatched topology surface count makes the layout invalid.
  SparseTrackingTopology topology{3};
  BOOST_REQUIRE(topology.finalize());
  const std::vector<SurfaceDescriptor> surfaces{surface(0), surface(1)};
  DetectorLayout layout{surfaces, std::move(topology)};
  BOOST_CHECK(!layout.valid());
  const auto masks = computeSurfaceKindMasks(surfaces);
  const auto view = layout.getView(surfaces, masks.first, masks.second);
  BOOST_CHECK(view.surfaces == nullptr);
  BOOST_CHECK_EQUAL(view.nSurfaces, 0u);
}

BOOST_AUTO_TEST_CASE(GetViewValidForLegitimateEmptyGeometry)
{
  auto layout = buildLayout({});
  BOOST_REQUIRE(layout.valid());
  const auto view = layout.getView({}, {}, {});
  BOOST_CHECK_EQUAL(view.nSurfaces, 0u);
}

BOOST_AUTO_TEST_CASE(DetectorLayoutSetProducedViewsRemainValid)
{
  std::vector<SurfaceDescriptor> catalog{surface(0), surface(1), surface(2)};
  auto set = buildSet(catalog, 2);
  BOOST_CHECK(set.getLayoutView(0).surfaces != nullptr);
  BOOST_CHECK(set.getLayoutView(1).surfaces != nullptr);
}
