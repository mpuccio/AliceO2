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

#include "ITSMFTTracking/SurfaceMeasurementIndex.h"
#include "ITSMFTTracking/SurfaceSpec.h"

using namespace o2::itsmft::tracking;

namespace
{
constexpr StaticSurfaceDescriptor cylinder(uint16_t id, uint8_t detector, uint16_t local, float radius = 2.f,
                                           float zMin = -10.f, float zMax = 10.f, float material = 0.01f)
{
  return {SurfaceId{id},
          {detector, local},
          SurfaceKind::Cylinder,
          radius,
          SurfaceAcceptance::fromCylinder({zMin, zMax}),
          {material},
          SurfaceIndexingFamily::CylindricalPhiZ};
}

constexpr StaticSurfaceDescriptor disk(uint16_t id, uint8_t detector, uint16_t local, float z = -40.f,
                                       float radiusMin = 1.f, float radiusMax = 20.f, float material = 0.02f)
{
  return {SurfaceId{id},
          {detector, local},
          SurfaceKind::Disk,
          z,
          SurfaceAcceptance::fromDisk({radiusMin, radiusMax}),
          {material},
          SurfaceIndexingFamily::CartesianXY};
}

struct Cylinders {
  inline static constexpr std::array surfaces{cylinder(0, 37, 0, 2.f), cylinder(1, 37, 1, 4.f)};
};

struct Disks {
  inline static constexpr std::array surfaces{disk(0, 201, 0, -40.f), disk(1, 201, 1, -50.f)};
};

using Combined = ConcatenatedSurfaceSpec<Cylinders, Disks>;

template <typename T>
concept HasTransitions = requires { T::transitions; };

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

constexpr auto duplicateId = [](auto& surfaces) { surfaces[1].id = SurfaceId{0}; };
constexpr auto sparseId = [](auto& surfaces) { surfaces[1].id = SurfaceId{2}; };
constexpr auto duplicateIdentity = [](auto& surfaces) { surfaces[1].identity = surfaces[0].identity; };
constexpr auto sparseLocalIdentity = [](auto& surfaces) { surfaces[1].identity.detectorSurfaceIndex = 2; };
constexpr auto invalidKind = [](auto& surfaces) { surfaces[0].kind = static_cast<SurfaceKind>(0xff); };
constexpr auto invalidIndexing = [](auto& surfaces) { surfaces[0].indexingFamily = SurfaceIndexingFamily::Invalid; };
constexpr auto mismatchedAcceptance = [](auto& surfaces) {
  surfaces[0].nominalTrackingAcceptance = SurfaceAcceptance::fromDisk({1.f, 2.f});
};
constexpr auto invalidAcceptance = [](auto& surfaces) { surfaces[0].nominalTrackingAcceptance = {}; };
constexpr auto nanCoordinate = [](auto& surfaces) {
  surfaces[0].nominalReferenceCoordinate = std::numeric_limits<float>::quiet_NaN();
};
constexpr auto infiniteCoordinate = [](auto& surfaces) {
  surfaces[0].nominalReferenceCoordinate = std::numeric_limits<float>::infinity();
};
constexpr auto zeroCylinderRadius = [](auto& surfaces) { surfaces[0].nominalReferenceCoordinate = 0.f; };
constexpr auto negativeCylinderRadius = [](auto& surfaces) { surfaces[0].nominalReferenceCoordinate = -1.f; };
constexpr auto nanAcceptance = [](auto& surfaces) {
  surfaces[0].nominalTrackingAcceptance = SurfaceAcceptance::fromCylinder(
    {std::numeric_limits<float>::quiet_NaN(), 1.f});
};
constexpr auto infiniteAcceptance = [](auto& surfaces) {
  surfaces[0].nominalTrackingAcceptance = SurfaceAcceptance::fromCylinder(
    {-1.f, std::numeric_limits<float>::infinity()});
};
constexpr auto reversedAcceptance = [](auto& surfaces) {
  surfaces[0].nominalTrackingAcceptance = SurfaceAcceptance::fromCylinder({2.f, 1.f});
};
constexpr auto zeroMaterial = [](auto& surfaces) { surfaces[0].material.xOverX0 = 0.f; };
constexpr auto negativeMaterial = [](auto& surfaces) { surfaces[0].material.xOverX0 = -0.01f; };
constexpr auto nanMaterial = [](auto& surfaces) {
  surfaces[0].material.xOverX0 = std::numeric_limits<float>::quiet_NaN();
};
constexpr auto infiniteMaterial = [](auto& surfaces) {
  surfaces[0].material.xOverX0 = std::numeric_limits<float>::infinity();
};

struct NegativeDiskRadius {
  inline static constexpr std::array surfaces{disk(0, 201, 0, -40.f, -1.f, 20.f)};
};

struct ReversedDiskBounds {
  inline static constexpr std::array surfaces{disk(0, 201, 0, -40.f, 20.f, 1.f)};
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
static_assert(!HasTransitions<Cylinders> && !HasTopology<Cylinders>);
static_assert(!HasTransitions<Combined> && !HasTopology<Combined>);

static_assert(!validateSurfaceSpec<MutatedCylinders<duplicateId>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<sparseId>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<duplicateIdentity>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<sparseLocalIdentity>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<invalidKind>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<invalidIndexing>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<mismatchedAcceptance>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<invalidAcceptance>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<nanCoordinate>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<infiniteCoordinate>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<zeroCylinderRadius>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<negativeCylinderRadius>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<nanAcceptance>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<infiniteAcceptance>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<reversedAcceptance>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<zeroMaterial>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<negativeMaterial>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<nanMaterial>>());
static_assert(!validateSurfaceSpec<MutatedCylinders<infiniteMaterial>>());
static_assert(!validateSurfaceSpec<NegativeDiskRadius>());
static_assert(!validateSurfaceSpec<ReversedDiskBounds>());
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

static_assert(std::is_standard_layout_v<DetectorSurfaceIdentity>);
static_assert(std::is_trivially_copyable_v<DetectorSurfaceIdentity>);
static_assert(std::is_standard_layout_v<SurfaceAcceptance>);
static_assert(std::is_trivially_copyable_v<SurfaceAcceptance>);
static_assert(std::is_standard_layout_v<StaticSurfaceDescriptor>);
static_assert(std::is_trivially_copyable_v<StaticSurfaceDescriptor>);
static_assert(std::is_standard_layout_v<SurfaceMeasurementIndex>);
static_assert(std::is_trivially_copyable_v<SurfaceMeasurementIndex>);

static_assert(!std::is_convertible_v<uint32_t, SurfaceMeasurementIndex>);
static_assert(!SurfaceMeasurementIndex{}.isValid());
static_assert(SurfaceMeasurementIndex{0}.isValid());
static_assert(SurfaceMeasurementIndex{0xfffffffeu}.isValid());
static_assert(!SurfaceMeasurementIndex{0xffffffffu}.isValid());
static_assert(SurfaceMeasurementIndex{7} == SurfaceMeasurementIndex{7});
static_assert(SurfaceMeasurementIndex{7} != SurfaceMeasurementIndex{8});

} // namespace

BOOST_AUTO_TEST_CASE(ConcatenationRebasesAndPreservesFields)
{
  BOOST_CHECK_EQUAL(Combined::surfaces[0].id.value(), 0);
  BOOST_CHECK_EQUAL(Combined::surfaces[1].id.value(), 1);
  BOOST_CHECK_EQUAL(Combined::surfaces[2].id.value(), 2);
  BOOST_CHECK_EQUAL(Combined::surfaces[3].id.value(), 3);
  BOOST_CHECK((Combined::surfaces[0].identity == DetectorSurfaceIdentity{37, 0}));
  BOOST_CHECK((Combined::surfaces[1].identity == DetectorSurfaceIdentity{37, 1}));
  BOOST_CHECK((Combined::surfaces[2].identity == DetectorSurfaceIdentity{201, 0}));
  BOOST_CHECK((Combined::surfaces[3].identity == DetectorSurfaceIdentity{201, 1}));
  BOOST_CHECK(Combined::surfaces[2].kind == Disks::surfaces[0].kind);
  BOOST_CHECK_EQUAL(Combined::surfaces[2].nominalReferenceCoordinate, Disks::surfaces[0].nominalReferenceCoordinate);
  BOOST_CHECK_EQUAL(Combined::surfaces[2].nominalTrackingAcceptance.min, Disks::surfaces[0].nominalTrackingAcceptance.min);
  BOOST_CHECK_EQUAL(Combined::surfaces[2].material.xOverX0, Disks::surfaces[0].material.xOverX0);
  BOOST_CHECK(Combined::surfaces[2].indexingFamily == Disks::surfaces[0].indexingFamily);
}

BOOST_AUTO_TEST_CASE(RuntimeProjectionHasIdealStaticSemantics)
{
  const auto projectedCylinder = toRuntimeSurfaceDescriptor(Cylinders::surfaces[1]);
  BOOST_CHECK(projectedCylinder.id == SurfaceId{1});
  BOOST_CHECK_EQUAL(projectedCylinder.detectorId, 37);
  BOOST_CHECK_EQUAL(projectedCylinder.detectorSurfaceIndex, 1);
  BOOST_CHECK(projectedCylinder.kind == SurfaceKind::Cylinder);
  BOOST_CHECK_EQUAL(projectedCylinder.flags, 0);
  BOOST_CHECK_EQUAL(projectedCylinder.referenceCoordinate, 4.f);
  BOOST_CHECK_EQUAL(projectedCylinder.radialMin, 4.f);
  BOOST_CHECK_EQUAL(projectedCylinder.radialMax, 4.f);

  const auto projectedDisk = toRuntimeSurfaceDescriptor(Disks::surfaces[0]);
  BOOST_CHECK(projectedDisk.id == SurfaceId{0});
  BOOST_CHECK_EQUAL(projectedDisk.detectorId, 201);
  BOOST_CHECK_EQUAL(projectedDisk.detectorSurfaceIndex, 0);
  BOOST_CHECK(projectedDisk.kind == SurfaceKind::Disk);
  BOOST_CHECK_EQUAL(projectedDisk.flags, 0);
  BOOST_CHECK_EQUAL(projectedDisk.referenceCoordinate, -40.f);
  BOOST_CHECK_EQUAL(projectedDisk.radialMin, 1.f);
  BOOST_CHECK_EQUAL(projectedDisk.radialMax, 20.f);
}

BOOST_AUTO_TEST_CASE(SurfaceMeasurementIndexBoundaries)
{
  BOOST_CHECK(!SurfaceMeasurementIndex{}.isValid());
  BOOST_CHECK_EQUAL(SurfaceMeasurementIndex{}.value(), SurfaceMeasurementIndex::InvalidValue);
  BOOST_CHECK(SurfaceMeasurementIndex{0}.isValid());
  BOOST_CHECK(SurfaceMeasurementIndex{0xfffffffeu}.isValid());
  BOOST_CHECK(!SurfaceMeasurementIndex{0xffffffffu}.isValid());
}
