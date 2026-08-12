// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Static detector definitions with process-lifetime storage.

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORDEFINITIONS_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORDEFINITIONS_H_

#include <array>
#include <cstddef>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/SurfaceSpec.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking
{

static_assert(MFTNLayers % 2 == 0);
inline constexpr int MFTDisks = MFTNLayers / 2;
inline constexpr std::array<float, ITSNLayers> kNominalITSLayerX0{
  5.e-3f, 5.e-3f, 5.e-3f, 1.e-2f, 1.e-2f, 1.e-2f, 1.e-2f};
inline constexpr float kMFTNominalRadLength = 0.042f;

constexpr std::array<float, MFTNLayers> makeNominalMFTLayerX0()
{
  std::array<float, MFTNLayers> values{};
  for (auto& value : values) {
    value = kMFTNominalRadLength / static_cast<float>(MFTDisks);
  }
  return values;
}

inline constexpr std::array<float, MFTNLayers> kNominalMFTLayerX0 = makeNominalMFTLayerX0();

constexpr NominalSurfaceMaterial itsLayerMaterial(std::size_t layer) noexcept
{
  const float x0 = kNominalITSLayerX0[layer];
  return {x0, x0 * o2::its::constants::Radl * o2::its::constants::Rho};
}

struct ITSSurfaceSpec {
  inline static constexpr std::array<StaticSurfaceDescriptor, ITSNLayers> surfaces{
    StaticSurfaceDescriptor{SurfaceId{0}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 0}, SurfaceKind::Cylinder, 2.3259652f, itsLayerMaterial(0)},
    StaticSurfaceDescriptor{SurfaceId{1}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 1}, SurfaceKind::Cylinder, 3.1353536f, itsLayerMaterial(1)},
    StaticSurfaceDescriptor{SurfaceId{2}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 2}, SurfaceKind::Cylinder, 3.9162421f, itsLayerMaterial(2)},
    StaticSurfaceDescriptor{SurfaceId{3}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 3}, SurfaceKind::Cylinder, 19.58824f, itsLayerMaterial(3)},
    StaticSurfaceDescriptor{SurfaceId{4}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 4}, SurfaceKind::Cylinder, 24.527159f, itsLayerMaterial(4)},
    StaticSurfaceDescriptor{SurfaceId{5}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 5}, SurfaceKind::Cylinder, 34.354595f, itsLayerMaterial(5)},
    StaticSurfaceDescriptor{SurfaceId{6}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 6}, SurfaceKind::Cylinder, 39.310642f, itsLayerMaterial(6)},
  };
};

static_assert(SurfaceSpec<ITSSurfaceSpec>);
static_assert(SurfaceCount<ITSSurfaceSpec> == ITSNLayers);

constexpr NominalSurfaceMaterial mftLayerMaterial(std::size_t layer) noexcept
{
  const float x0 = kNominalMFTLayerX0[layer];
  return {x0, x0 * o2::its::constants::Radl * o2::its::constants::Rho};
}

struct MFTSurfaceSpec {
  inline static constexpr std::array<StaticSurfaceDescriptor, MFTNLayers> surfaces{
    StaticSurfaceDescriptor{SurfaceId{0}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 0}, SurfaceKind::Disk, -45.2889f, mftLayerMaterial(0)},
    StaticSurfaceDescriptor{SurfaceId{1}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 1}, SurfaceKind::Disk, -46.7111f, mftLayerMaterial(1)},
    StaticSurfaceDescriptor{SurfaceId{2}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 2}, SurfaceKind::Disk, -48.5889f, mftLayerMaterial(2)},
    StaticSurfaceDescriptor{SurfaceId{3}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 3}, SurfaceKind::Disk, -50.0111f, mftLayerMaterial(3)},
    StaticSurfaceDescriptor{SurfaceId{4}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 4}, SurfaceKind::Disk, -52.3889f, mftLayerMaterial(4)},
    StaticSurfaceDescriptor{SurfaceId{5}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 5}, SurfaceKind::Disk, -53.8111f, mftLayerMaterial(5)},
    StaticSurfaceDescriptor{SurfaceId{6}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 6}, SurfaceKind::Disk, -67.6889f, mftLayerMaterial(6)},
    StaticSurfaceDescriptor{SurfaceId{7}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 7}, SurfaceKind::Disk, -69.1111f, mftLayerMaterial(7)},
    StaticSurfaceDescriptor{SurfaceId{8}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 8}, SurfaceKind::Disk, -76.0889f, mftLayerMaterial(8)},
    StaticSurfaceDescriptor{SurfaceId{9}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 9}, SurfaceKind::Disk, -77.5111f, mftLayerMaterial(9)},
  };
};

static_assert(SurfaceSpec<MFTSurfaceSpec>);
static_assert(SurfaceCount<MFTSurfaceSpec> == MFTNLayers);

template <SurfaceSpec Spec>
consteval std::array<SurfaceDescriptor, SurfaceCount<Spec>> projectStaticSurfaceCatalog() noexcept
{
  std::array<SurfaceDescriptor, SurfaceCount<Spec>> result{};
  for (std::size_t i = 0; i < SurfaceCount<Spec>; ++i) {
    result[i] = toRuntimeSurfaceDescriptor(Spec::surfaces[i]);
  }
  return result;
}

inline constexpr auto kITSStaticSurfaceCatalog = projectStaticSurfaceCatalog<ITSSurfaceSpec>();
inline constexpr auto kMFTStaticSurfaceCatalog = projectStaticSurfaceCatalog<MFTSurfaceSpec>();

static_assert(kITSStaticSurfaceCatalog.size() == ITSNLayers);
static_assert(kMFTStaticSurfaceCatalog.size() == MFTNLayers);

// One combined ITS+MFT global id space, compile-time concatenated via the
// existing ConcatenatedSurfaceSpec mechanism -- not a new abstraction.
// concatenateAndRebase() rebases `id` densely across the two specs while
// leaving each StaticSurfaceDescriptor::identity (detector-qualified: own
// detectorId, own local detectorSurfaceIndex) untouched, so every surface's
// detector-local identity survives global rebasing unchanged.
using CombinedITSMFTSurfaceSpec = ConcatenatedSurfaceSpec<ITSSurfaceSpec, MFTSurfaceSpec>;
static_assert(SurfaceSpecsCanBeConcatenated<ITSSurfaceSpec, MFTSurfaceSpec>);
static_assert(SurfaceSpec<CombinedITSMFTSurfaceSpec>);

inline constexpr auto kITSMFTCombinedStaticSurfaceCatalog = projectStaticSurfaceCatalog<CombinedITSMFTSurfaceSpec>();

static_assert(kITSMFTCombinedStaticSurfaceCatalog.size() == ITSNLayers + MFTNLayers);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTORDEFINITIONS_H_ */
