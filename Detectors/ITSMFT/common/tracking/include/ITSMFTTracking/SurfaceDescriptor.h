// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEDESCRIPTOR_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEDESCRIPTOR_H_

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/IdTypes.h"

namespace o2::itsmft::tracking
{

// Nominal normal-incidence material; zero denotes material not configured.
struct NominalSurfaceMaterial {
  float xOverX0{0.f};
  float arealDensityGPerCm2{0.f};
};

struct SurfaceChartRange {
  float min{0.f};
  float max{0.f};

  GPUhdi() constexpr bool isValid() const noexcept { return min < max; }
};

static_assert(std::is_standard_layout_v<NominalSurfaceMaterial>);
static_assert(std::is_trivially_copyable_v<NominalSurfaceMaterial>);
static_assert(sizeof(NominalSurfaceMaterial) == 8);
static_assert(alignof(NominalSurfaceMaterial) == 4);
static_assert(offsetof(NominalSurfaceMaterial, xOverX0) == 0);
static_assert(offsetof(NominalSurfaceMaterial, arealDensityGPerCm2) == 4);
static_assert(std::is_standard_layout_v<SurfaceChartRange>);
static_assert(std::is_trivially_copyable_v<SurfaceChartRange>);
static_assert(sizeof(SurfaceChartRange) == 8);

// Immutable surface identity, geometry and nominal material.
struct SurfaceDescriptor {
  SurfaceId id{};
  uint16_t detectorSurfaceIndex{0};
  uint8_t detectorId{0};
  SurfaceKind kind{SurfaceKind::Undefined};
  uint16_t flags{0};
  float referenceCoordinate{0.f}; // nominal radius for cylinders, z for disks
  NominalSurfaceMaterial material{};
  SurfaceChartRange chartRange{};
};

static_assert(std::is_standard_layout_v<SurfaceDescriptor>);
static_assert(std::is_trivially_copyable_v<SurfaceDescriptor>);
static_assert(sizeof(SurfaceDescriptor) == 28);
static_assert(alignof(SurfaceDescriptor) == 4);
static_assert(offsetof(SurfaceDescriptor, id) == 0);
static_assert(offsetof(SurfaceDescriptor, detectorSurfaceIndex) == 2);
static_assert(offsetof(SurfaceDescriptor, detectorId) == 4);
static_assert(offsetof(SurfaceDescriptor, kind) == 5);
static_assert(offsetof(SurfaceDescriptor, flags) == 6);
static_assert(offsetof(SurfaceDescriptor, referenceCoordinate) == 8);
static_assert(offsetof(SurfaceDescriptor, material) == 12);
static_assert(offsetof(SurfaceDescriptor, chartRange) == 20);

// Non-owning surface-catalog view. Topology, timing and measurements stay
// outside this POD so loading and propagation do not depend on them.
struct SurfaceCatalogView {
  const SurfaceDescriptor* surfaces{nullptr};
  uint32_t nSurfaces{0};
  // Optional lookup for non-dense IDs; null keeps dense catalog indexing.
  const uint8_t* surfaceIndicesById{nullptr};

  GPUhdi() uint32_t getSurfaceIndex(SurfaceId id) const
  {
    if (!id.isValid() || id.value() >= MaxLayoutSurfaces) {
      return nSurfaces;
    }
    if (surfaceIndicesById != nullptr) {
      return surfaceIndicesById[id.value()];
    }
    if (id.value() < nSurfaces && surfaces[id.value()].id == id) {
      return id.value();
    }
    for (uint32_t i = 0; i < nSurfaces; ++i) {
      if (surfaces[i].id == id) {
        return i;
      }
    }
    return nSurfaces;
  }

  GPUhdi() bool hasSurface(SurfaceId id) const { return getSurfaceIndex(id) < nSurfaces; }
  GPUhdi() const SurfaceDescriptor& getSurface(SurfaceId id) const { return surfaces[getSurfaceIndex(id)]; }
};

static_assert(std::is_standard_layout_v<SurfaceCatalogView>);
static_assert(std::is_trivially_copyable_v<SurfaceCatalogView>);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACEDESCRIPTOR_H_ */
