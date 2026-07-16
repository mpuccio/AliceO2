// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEDESCRIPTOR_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEDESCRIPTOR_H_

#include <cstdint>
#include <type_traits>

#include "ITSMFTTracking/SurfaceId.h"

namespace o2::itsmft::tracking
{

enum class SurfaceKind : uint8_t {
  Cylinder,
  Disk
};

// Geometry identity shared by host and device code. Measurement, timing,
// material and indexing descriptors are deliberately composed separately.
struct SurfaceDescriptor {
  SurfaceId id{};
  uint16_t detectorSurfaceIndex{0};
  uint8_t detectorId{0};
  SurfaceKind kind{SurfaceKind::Cylinder};
  uint16_t flags{0};
  float referenceCoordinate{0.f}; // nominal radius for cylinders, z for disks
  float radialMin{0.f};
  float radialMax{0.f};
};

static_assert(std::is_standard_layout_v<SurfaceDescriptor>);
static_assert(std::is_trivially_copyable_v<SurfaceDescriptor>);
static_assert(sizeof(SurfaceDescriptor) == 20);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACEDESCRIPTOR_H_ */
