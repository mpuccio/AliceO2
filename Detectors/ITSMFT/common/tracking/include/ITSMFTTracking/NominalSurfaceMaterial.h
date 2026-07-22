// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_NOMINALSURFACEMATERIAL_H_
#define ALICEO2_ITSMFT_TRACKING_NOMINALSURFACEMATERIAL_H_

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "ITSMFTTracking/SurfaceId.h"

// This header is deliberately narrow: no STL container, gsl::span, owner,
// provider, material-query API, or TGeo/geometry dependency. It defines only
// the per-surface material payload and the identity-bearing construction
// input used to validate it positionally against a dense surface catalog
// (see DetectorLayoutSet).

namespace o2::itsmft::tracking
{

// Nominal (perpendicular-incidence) per-surface material budget. Both
// fields are non-negative by construction; converting this into a
// direction/incidence-corrected budget is out of scope for this slice.
struct NominalSurfaceMaterialBudget {
  float normalXOverX0{0.f};
  float normalArealDensityGPerCm2{0.f};
};

// One positional construction input: the SurfaceId an entry claims to
// describe, paired with its nominal budget. The owner validates this claim
// positionally against a dense surface catalog; the dense catalog order,
// not this struct alone, is the source of truth for identity.
struct NominalSurfaceMaterialEntry {
  SurfaceId surface{};
  NominalSurfaceMaterialBudget budget{};
};

#define O2_ITSMFT_ASSERT_NOMINAL_MATERIAL_TYPE(Type, Size, Alignment) \
  static_assert(std::is_standard_layout_v<Type>);                     \
  static_assert(std::is_trivially_copyable_v<Type>);                  \
  static_assert(sizeof(Type) == Size);                                \
  static_assert(alignof(Type) == Alignment)

O2_ITSMFT_ASSERT_NOMINAL_MATERIAL_TYPE(NominalSurfaceMaterialBudget, 8, 4);
O2_ITSMFT_ASSERT_NOMINAL_MATERIAL_TYPE(NominalSurfaceMaterialEntry, 12, 4);

static_assert(offsetof(NominalSurfaceMaterialBudget, normalXOverX0) == 0);
static_assert(offsetof(NominalSurfaceMaterialBudget, normalArealDensityGPerCm2) == 4);
static_assert(offsetof(NominalSurfaceMaterialEntry, surface) == 0);
static_assert(offsetof(NominalSurfaceMaterialEntry, budget) == 4);

#undef O2_ITSMFT_ASSERT_NOMINAL_MATERIAL_TYPE

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_NOMINALSURFACEMATERIAL_H_ */
