// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#define BOOST_TEST_MODULE Test ITSMFTTracking SurfaceSpec
#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "ITSMFTTracking/IdTypes.h"
#include "ITSMFTTracking/SurfaceSpec.h"

using namespace o2::itsmft::tracking;

namespace
{
constexpr StaticSurfaceDescriptor cylinder(uint16_t id, uint8_t detector, uint16_t local, float radius = 2.f,
                                           float material = 0.01f, float arealDensity = 0.f)
{
  return {LayerId{id},
          {detector, local},
          SurfaceKind::Cylinder,
          radius,
          {material, arealDensity}};
}

constexpr StaticSurfaceDescriptor disk(uint16_t id, uint8_t detector, uint16_t local, float z = -40.f,
                                       float material = 0.02f, float arealDensity = 0.f)
{
  return {LayerId{id},
          {detector, local},
          SurfaceKind::Disk,
          z,
          {material, arealDensity}};
}

struct Cylinders {
  inline static constexpr std::array surfaces{cylinder(0, 37, 0, 2.f), cylinder(1, 37, 1, 4.f)};
};

struct Disks {
  inline static constexpr std::array surfaces{disk(0, 201, 0, -40.f), disk(1, 201, 1, -50.f)};
};

using Combined = ConcatenatedSurfaceSpec<Cylinders, Disks>;

template <typename T>
concept HasEdges = requires { T::edges; };

template <typename T>
concept HasTopology = requires { T::topology; };

template <std::size_t N>
consteval auto cylindersFor(uint8_t detector)
{
  std::array<StaticSurfaceDescriptor, N> result{};
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = cylinder(static_cast<uint16_t>(i), detector, static_cast<uint16_t>(i));
  }
  return result;
}

struct ThirtyTwoSurfaces {
  inline static constexpr auto surfaces = cylindersFor<32>(93);
};

struct ThirtyOneSurfaces {
  inline static constexpr auto surfaces = cylindersFor<31>(93);
};

struct ThirtyThreeSurfaces {
  inline static constexpr auto surfaces = cylindersFor<33>(93);
};

struct OneSurface {
  inline static constexpr std::array surfaces{cylinder(0, 94, 0)};
};

template <auto Mutator>
struct MutatedCylinders {
  inline static constexpr auto surfaces = [] {
    auto result = Cylinders::surfaces;
    Mutator(result);
    return result;
  }();
};

constexpr auto duplicateId = [](auto& surfaces) { surfaces[1].id = LayerId{0}; };
constexpr auto sparseId = [](auto& surfaces) { surfaces[1].id = LayerId{2}; };
constexpr auto duplicateIdentity = [](auto& surfaces) { surfaces[1].identity = surfaces[0].identity; };
constexpr auto sparseLocalIdentity = [](auto& surfaces) { surfaces[1].identity.detectorSurfaceIndex = 2; };
constexpr auto invalidKind = [](auto& surfaces) { surfaces[0].kind = static_cast<SurfaceKind>(0xff); };
constexpr auto nanCoordinate = [](auto& surfaces) {
  surfaces[0].nominalReferenceCoordinate = std::numeric_limits<float>::quiet_NaN();
};
constexpr auto infiniteCoordinate = [](auto& surfaces) {
  surfaces[0].nominalReferenceCoordinate = std::numeric_limits<float>::infinity();
};
constexpr auto zeroCylinderRadius = [](auto& surfaces) { surfaces[0].nominalReferenceCoordinate = 0.f; };
constexpr auto negativeCylinderRadius = [](auto& surfaces) { surfaces[0].nominalReferenceCoordinate = -1.f; };
constexpr auto zeroMaterial = [](auto& surfaces) { surfaces[0].material.xOverX0 = 0.f; };
constexpr auto negativeMaterial = [](auto& surfaces) { surfaces[0].material.xOverX0 = -0.01f; };
constexpr auto nanMaterial = [](auto& surfaces) {
  surfaces[0].material.xOverX0 = std::numeric_limits<float>::quiet_NaN();
};
constexpr auto infiniteMaterial = [](auto& surfaces) {
  surfaces[0].material.xOverX0 = std::numeric_limits<float>::infinity();
};
constexpr auto zeroArealDensity = [](auto& surfaces) { surfaces[0].material.arealDensityGPerCm2 = 0.f; };
constexpr auto negativeArealDensity = [](auto& surfaces) { surfaces[0].material.arealDensityGPerCm2 = -0.02f; };
constexpr auto nanArealDensity = [](auto& surfaces) {
  surfaces[0].material.arealDensityGPerCm2 = std::numeric_limits<float>::quiet_NaN();
};
constexpr auto infiniteArealDensity = [](auto& surfaces) {
  surfaces[0].material.arealDensityGPerCm2 = std::numeric_limits<float>::infinity();
};
constexpr auto zeroBothMaterialFields = [](auto& surfaces) {
  surfaces[0].material.xOverX0 = 0.f;
  surfaces[0].material.arealDensityGPerCm2 = 0.f;
};

struct CrossBoundaryIdentityCollision {
  inline static constexpr std::array surfaces{cylinder(0, 37, 0)};
};

using ThirtyTwoCombined = ConcatenatedSurfaceSpec<ThirtyOneSurfaces, OneSurface>;

static_assert(SurfaceSpec<Cylinders> && SurfaceSpec<Disks> && SurfaceSpec<Combined>);
static_assert(SurfaceSpecDefinition<MutatedCylinders<sparseId>>);
static_assert(!SurfaceSpec<MutatedCylinders<sparseId>>);
static_assert(SurfaceSpecDefinition<ThirtyThreeSurfaces>);
static_assert(!SurfaceSpec<ThirtyThreeSurfaces>);
static_assert(validateSurfaceSpec<Cylinders>());
static_assert(validateSurfaceSpec<Disks>());
static_assert(validateSurfaceSpec<Combined>());
static_assert(SurfaceCount<Cylinders> == 2);
static_assert(SurfaceCount<Disks> == 2);
static_assert(SurfaceCount<Combined> == 4);
static_assert(Cylinders::surfaces.data() == &Cylinders::surfaces[0]);
static_assert(!HasEdges<Cylinders> && !HasTopology<Cylinders>);
static_assert(!HasEdges<Combined> && !HasTopology<Combined>);

static_assert(!validateSurfaceSpec<MutatedCylinders<duplicateId>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<sparseId>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<duplicateIdentity>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<sparseLocalIdentity>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<invalidKind>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<nanCoordinate>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<infiniteCoordinate>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<zeroCylinderRadius>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<negativeCylinderRadius>>());
// Zero is legal, independently, for either material field: a surface with no
// material configured yet is a valid SurfaceSpec, not a rejected one.
static_assert(validateSurfaceSpec<MutatedCylinders<zeroMaterial>>());
static_assert(validateSurfaceSpec<MutatedCylinders<zeroArealDensity>>());
static_assert(validateSurfaceSpec<MutatedCylinders<zeroBothMaterialFields>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<negativeMaterial>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<nanMaterial>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<infiniteMaterial>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<negativeArealDensity>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<nanArealDensity>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<infiniteArealDensity>>());
static_assert(validateSurfaceSpec<ThirtyTwoSurfaces>());
static_assert(!validateSurfaceSpec<ThirtyThreeSurfaces>());
static_assert(SurfaceSpecsCanBeConcatenated<Cylinders, Disks>);
static_assert(!SurfaceSpecsCanBeConcatenated<MutatedCylinders<sparseId>, Disks>);
static_assert(!SurfaceSpecsCanBeConcatenated<Disks, MutatedCylinders<sparseId>>);
static_assert(!SurfaceSpecsCanBeConcatenated<Cylinders, CrossBoundaryIdentityCollision>);
static_assert(SurfaceSpecsCanBeConcatenated<ThirtyOneSurfaces, OneSurface>);
static_assert(!SurfaceSpecsCanBeConcatenated<ThirtyTwoSurfaces, OneSurface>);
static_assert(SurfaceSpec<ThirtyTwoCombined>);
static_assert(SurfaceCount<ThirtyTwoCombined> == 32);

static_assert(std::is_standard_layout_v<DetectorLayerIdentity>);
static_assert(std::is_trivially_copyable_v<DetectorLayerIdentity>);
static_assert(std::is_standard_layout_v<StaticSurfaceDescriptor>);
static_assert(std::is_trivially_copyable_v<StaticSurfaceDescriptor>);

} // namespace

// -------------------------------------------------------------------------
// Runtime ABI/layout lock, alongside the compile-time static_asserts in
// SurfaceDescriptor.h/SurfaceSpec.h: nominal material is a field
// directly on both descriptor types, not a parallel array or pointer.
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(SurfaceDescriptorRuntimeAbiLock)
{
  BOOST_CHECK_EQUAL(sizeof(NominalSurfaceMaterial), 8u);
  BOOST_CHECK_EQUAL(alignof(NominalSurfaceMaterial), 4u);
  BOOST_CHECK_EQUAL(offsetof(NominalSurfaceMaterial, xOverX0), 0u);
  BOOST_CHECK_EQUAL(offsetof(NominalSurfaceMaterial, arealDensityGPerCm2), 4u);

  BOOST_CHECK_EQUAL(sizeof(SurfaceDescriptor), 28u);
  BOOST_CHECK_EQUAL(alignof(SurfaceDescriptor), 4u);
  BOOST_CHECK_EQUAL(offsetof(SurfaceDescriptor, material), 12u);

  BOOST_CHECK_EQUAL(sizeof(StaticSurfaceDescriptor), 28u);
  BOOST_CHECK_EQUAL(alignof(StaticSurfaceDescriptor), 4u);
  BOOST_CHECK_EQUAL(offsetof(StaticSurfaceDescriptor, material), 12u);

  // Default-constructed material on both descriptor types is zero on both
  // fields independently -- not just "some default", but exactly zero.
  constexpr SurfaceDescriptor defaultRuntime{};
  BOOST_CHECK_EQUAL(defaultRuntime.material.xOverX0, 0.f);
  BOOST_CHECK_EQUAL(defaultRuntime.material.arealDensityGPerCm2, 0.f);
  constexpr StaticSurfaceDescriptor defaultStatic{};
  BOOST_CHECK_EQUAL(defaultStatic.material.xOverX0, 0.f);
  BOOST_CHECK_EQUAL(defaultStatic.material.arealDensityGPerCm2, 0.f);
}

BOOST_AUTO_TEST_CASE(ConcatenationRebasesAndPreservesFields)
{
  BOOST_CHECK_EQUAL(Combined::surfaces[0].id.value(), 0);
  BOOST_CHECK_EQUAL(Combined::surfaces[1].id.value(), 1);
  BOOST_CHECK_EQUAL(Combined::surfaces[2].id.value(), 2);
  BOOST_CHECK_EQUAL(Combined::surfaces[3].id.value(), 3);
  BOOST_CHECK((Combined::surfaces[0].identity == DetectorLayerIdentity{37, 0}));
  BOOST_CHECK((Combined::surfaces[1].identity == DetectorLayerIdentity{37, 1}));
  BOOST_CHECK((Combined::surfaces[2].identity == DetectorLayerIdentity{201, 0}));
  BOOST_CHECK((Combined::surfaces[3].identity == DetectorLayerIdentity{201, 1}));
  BOOST_CHECK(Combined::surfaces[2].kind == Disks::surfaces[0].kind);
  BOOST_CHECK_EQUAL(Combined::surfaces[2].nominalReferenceCoordinate, Disks::surfaces[0].nominalReferenceCoordinate);
  BOOST_CHECK_EQUAL(Combined::surfaces[2].material.xOverX0, Disks::surfaces[0].material.xOverX0);
}

BOOST_AUTO_TEST_CASE(RuntimeProjectionHasIdealStaticSemantics)
{
  const auto projectedCylinder = toRuntimeSurfaceDescriptor(Cylinders::surfaces[1]);
  BOOST_CHECK(projectedCylinder.id == LayerId{1});
  BOOST_CHECK_EQUAL(projectedCylinder.detectorId, 37);
  BOOST_CHECK_EQUAL(projectedCylinder.detectorSurfaceIndex, 1);
  BOOST_CHECK(projectedCylinder.kind == SurfaceKind::Cylinder);
  BOOST_CHECK_EQUAL(projectedCylinder.flags, 0);
  BOOST_CHECK_EQUAL(projectedCylinder.referenceCoordinate, 4.f);
  BOOST_CHECK_EQUAL(projectedCylinder.material.xOverX0, Cylinders::surfaces[1].material.xOverX0);
  BOOST_CHECK_EQUAL(projectedCylinder.material.arealDensityGPerCm2, Cylinders::surfaces[1].material.arealDensityGPerCm2);

  const auto projectedDisk = toRuntimeSurfaceDescriptor(Disks::surfaces[0]);
  BOOST_CHECK(projectedDisk.id == LayerId{0});
  BOOST_CHECK_EQUAL(projectedDisk.detectorId, 201);
  BOOST_CHECK_EQUAL(projectedDisk.detectorSurfaceIndex, 0);
  BOOST_CHECK(projectedDisk.kind == SurfaceKind::Disk);
  BOOST_CHECK_EQUAL(projectedDisk.flags, 0);
  BOOST_CHECK_EQUAL(projectedDisk.referenceCoordinate, -40.f);
  BOOST_CHECK_EQUAL(projectedDisk.material.xOverX0, Disks::surfaces[0].material.xOverX0);
  BOOST_CHECK_EQUAL(projectedDisk.material.arealDensityGPerCm2, Disks::surfaces[0].material.arealDensityGPerCm2);
}

BOOST_AUTO_TEST_CASE(StaticToRuntimeMaterialProjectionCopiesBothFieldsIndependently)
{
  constexpr auto surface = cylinder(0, 37, 0, 2.f, 0.03f, 0.04f);
  const auto projected = toRuntimeSurfaceDescriptor(surface);
  BOOST_CHECK_EQUAL(projected.material.xOverX0, 0.03f);
  BOOST_CHECK_EQUAL(projected.material.arealDensityGPerCm2, 0.04f);

  constexpr auto zeroXOverX0 = cylinder(0, 37, 0, 2.f, 0.f, 0.05f);
  const auto projectedZeroXOverX0 = toRuntimeSurfaceDescriptor(zeroXOverX0);
  BOOST_CHECK_EQUAL(projectedZeroXOverX0.material.xOverX0, 0.f);
  BOOST_CHECK_EQUAL(projectedZeroXOverX0.material.arealDensityGPerCm2, 0.05f);

  constexpr auto zeroArealDensityValue = cylinder(0, 37, 0, 2.f, 0.06f, 0.f);
  const auto projectedZeroArealDensity = toRuntimeSurfaceDescriptor(zeroArealDensityValue);
  BOOST_CHECK_EQUAL(projectedZeroArealDensity.material.xOverX0, 0.06f);
  BOOST_CHECK_EQUAL(projectedZeroArealDensity.material.arealDensityGPerCm2, 0.f);
}
