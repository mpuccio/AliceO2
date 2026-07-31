// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_STATICSURFACEDESCRIPTOR_H_
#define ALICEO2_ITSMFT_TRACKING_STATICSURFACEDESCRIPTOR_H_

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"

namespace o2::itsmft::tracking
{

// Detector-qualified identity of a logical tracking surface. The detector ID
// remains open to every detector representable by the existing uint8_t field.
struct DetectorSurfaceIdentity {
  uint8_t detectorId{0};
  uint16_t detectorSurfaceIndex{0};

  GPUhdi() friend constexpr bool operator==(DetectorSurfaceIdentity lhs, DetectorSurfaceIdentity rhs) noexcept
  {
    return lhs.detectorId == rhs.detectorId && lhs.detectorSurfaceIndex == rhs.detectorSurfaceIndex;
  }
  GPUhdi() friend constexpr bool operator!=(DetectorSurfaceIdentity lhs, DetectorSurfaceIdentity rhs) noexcept { return !(lhs == rhs); }
};

// NominalSurfaceMaterial itself is defined once in SurfaceDescriptor.h and
// shared by both descriptor types; no duplicate definition here.

enum class SurfaceIndexingFamily : uint8_t {
  Invalid,
  CylindricalPhiZ,
  CartesianXY
};

struct StaticSurfaceDescriptor {
  SurfaceId id{};
  DetectorSurfaceIdentity identity{};
  SurfaceKind kind{SurfaceKind::Cylinder};
  float nominalReferenceCoordinate{0.f};
  NominalSurfaceMaterial material{};
  SurfaceIndexingFamily indexingFamily{SurfaceIndexingFamily::Invalid};
};

// Precondition: source belongs to a validated SurfaceSpec. This is an
// ideal/static layout projection only; it neither validates nor repairs an
// arbitrary descriptor. Runtime geometry observations are validation data and
// do not define a second Stage-B surface catalogue.
GPUhdi() constexpr SurfaceDescriptor toRuntimeSurfaceDescriptor(const StaticSurfaceDescriptor& source) noexcept
{
  return SurfaceDescriptor{source.id,
                           source.identity.detectorSurfaceIndex,
                           source.identity.detectorId,
                           source.kind,
                           0,
                           source.nominalReferenceCoordinate,
                           source.material};
}

#define O2_ITSMFT_ASSERT_STATIC_SURFACE_TYPE(Type, Size, Alignment) \
  static_assert(std::is_standard_layout_v<Type>);                   \
  static_assert(std::is_trivially_copyable_v<Type>);                \
  static_assert(sizeof(Type) == Size);                              \
  static_assert(alignof(Type) == Alignment)

O2_ITSMFT_ASSERT_STATIC_SURFACE_TYPE(DetectorSurfaceIdentity, 4, 2);
// NominalSurfaceMaterial's own layout is asserted once, in
// SurfaceDescriptor.h, where the shared type is defined.
O2_ITSMFT_ASSERT_STATIC_SURFACE_TYPE(StaticSurfaceDescriptor, 24, 4);

static_assert(offsetof(DetectorSurfaceIdentity, detectorId) == 0);
static_assert(offsetof(DetectorSurfaceIdentity, detectorSurfaceIndex) == 2);
static_assert(offsetof(StaticSurfaceDescriptor, id) == 0);
static_assert(offsetof(StaticSurfaceDescriptor, identity) == 2);
static_assert(offsetof(StaticSurfaceDescriptor, kind) == 6);
static_assert(offsetof(StaticSurfaceDescriptor, nominalReferenceCoordinate) == 8);
static_assert(offsetof(StaticSurfaceDescriptor, material) == 12);
static_assert(offsetof(StaticSurfaceDescriptor, indexingFamily) == 20);

#undef O2_ITSMFT_ASSERT_STATIC_SURFACE_TYPE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_STATICSURFACEDESCRIPTOR_H_ */
