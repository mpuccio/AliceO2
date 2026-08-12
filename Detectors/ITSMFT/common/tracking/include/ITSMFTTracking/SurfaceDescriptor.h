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

// Nominal (normal-incidence) per-surface material. An immutable property of
// the surface itself, alongside its geometry -- not a separate catalogue
// keyed by identity. Both fields are independently legal at zero (a surface
// with no material yet configured); shared by SurfaceDescriptor and
// StaticSurfaceDescriptor so there is exactly one definition of this type.
struct NominalSurfaceMaterial {
  float xOverX0{0.f};
  float arealDensityGPerCm2{0.f};
};

static_assert(std::is_standard_layout_v<NominalSurfaceMaterial>);
static_assert(std::is_trivially_copyable_v<NominalSurfaceMaterial>);
static_assert(sizeof(NominalSurfaceMaterial) == 8);
static_assert(alignof(NominalSurfaceMaterial) == 4);
static_assert(offsetof(NominalSurfaceMaterial, xOverX0) == 0);
static_assert(offsetof(NominalSurfaceMaterial, arealDensityGPerCm2) == 4);

// Geometry identity shared by host and device code, plus its nominal
// material. Measurement, timing and indexing descriptors are deliberately
// composed separately.
struct SurfaceDescriptor {
  SurfaceId id{};
  uint16_t detectorSurfaceIndex{0};
  uint8_t detectorId{0};
  SurfaceKind kind{SurfaceKind::Undefined};
  uint16_t flags{0};
  float referenceCoordinate{0.f}; // nominal radius for cylinders, z for disks
  NominalSurfaceMaterial material{};
};

static_assert(std::is_standard_layout_v<SurfaceDescriptor>);
static_assert(std::is_trivially_copyable_v<SurfaceDescriptor>);
static_assert(sizeof(SurfaceDescriptor) == 20);
static_assert(alignof(SurfaceDescriptor) == 4);
static_assert(offsetof(SurfaceDescriptor, id) == 0);
static_assert(offsetof(SurfaceDescriptor, detectorSurfaceIndex) == 2);
static_assert(offsetof(SurfaceDescriptor, detectorId) == 4);
static_assert(offsetof(SurfaceDescriptor, kind) == 5);
static_assert(offsetof(SurfaceDescriptor, flags) == 6);
static_assert(offsetof(SurfaceDescriptor, referenceCoordinate) == 8);
static_assert(offsetof(SurfaceDescriptor, material) == 12);

// Minimal, topology-free and mask-free view over a canonical
// surface catalog: a pointer to const SurfaceDescriptor plus a count,
// nothing else. Deliberately carries no topology, masks, transition dispatch,
// STL ownership or detector dependency, so consumers that only need the
// surface descriptions (e.g. loading) do not have to depend on graph
// adjacency or transition dispatch.
//
// A SurfaceDescriptor includes immutable identity, geometry and nominal
// material, so this view is not geometry-only: consumers that only care about
// identity/geometry may simply ignore the `material` field on each returned
// SurfaceDescriptor, but the view does not hide it. No separate material
// pointer or accessor is added here for that -- it is naturally reachable
// through the returned SurfaceDescriptor itself. What keeps this view narrow
// is what it deliberately excludes: topology, masks, transition policies,
// timing, measurements, and any runtime (endpoint-dependent) material-query
// result.
struct SurfaceCatalogView {
  const SurfaceDescriptor* surfaces{nullptr};
  uint32_t nSurfaces{0};
  // Optional compact lookup for catalogs whose global ids are not dense.
  // A null pointer preserves direct indexing for the canonical dense catalogs.
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
