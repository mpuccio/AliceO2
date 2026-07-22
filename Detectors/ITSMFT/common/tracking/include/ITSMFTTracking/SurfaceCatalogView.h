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
// surface descriptions (e.g. loading) do not have to depend on
// DetectorLayout/SparseTrackingTopology/TransitionPolicy.
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

  // Bounds-unchecked, matching DetectorLayoutView::getSurface().
  GPUhdi() const SurfaceDescriptor& getSurface(SurfaceId id) const { return surfaces[id.value()]; }
};

static_assert(std::is_standard_layout_v<SurfaceCatalogView>);
static_assert(std::is_trivially_copyable_v<SurfaceCatalogView>);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACECATALOGVIEW_H_ */
