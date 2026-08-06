// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACECATALOGVIEW_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACECATALOGVIEW_H_

#include <cstdint>
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceId.h"

namespace o2::itsmft::tracking
{

// Minimal, topology-free, mask-free and policy-free view over a canonical
// surface catalog: a pointer to const SurfaceDescriptor plus a count,
// nothing else. Deliberately carries no topology, masks, transition-policy,
// STL ownership or detector dependency, so consumers that only need the
// surface descriptions (e.g. loading) do not have to depend on graph
// adjacency or transition policy.
//
// A SurfaceDescriptor now includes immutable identity, geometry and nominal
// material (see SurfaceDescriptor.h), so this view is not geometry-only:
// consumers that only care about identity/geometry may simply ignore the
// `material` field on each returned SurfaceDescriptor, but the view does
// not hide it. No separate material pointer or accessor is added here for
// that -- it is naturally reachable through the returned SurfaceDescriptor
// itself. What keeps this view narrow is what it deliberately excludes:
// topology, masks, transition policies, timing, measurements, and any
// runtime (endpoint-dependent) material-query result.
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

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACECATALOGVIEW_H_ */
