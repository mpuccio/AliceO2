// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B1 Slice C3a: candidate static ITS DetectorSpec table.
//
// Existing mechanism only -- StaticSurfaceDescriptor, the SurfaceSpec
// concept, and toRuntimeSurfaceDescriptor() -- no new wrapper type, no
// runtime catalog, no acceptance/radial-envelope field (removed in Gate 4
// acceptance-cleanup C2; ADR 0003).
//
// nominalReferenceCoordinate literals below are std::to_chars' shortest
// round-tripping tokens, transcribed verbatim (not re-rounded) from the
// Gate 4 B1 Slice 1 detached validator's --format json output
// (NominalGeometryReport.cxx's formatLosslessFloat(), Gate 4
// acceptance-cleanup C1) against the checksum-pinned unaligned
// o2sim_geometry.root. Exact source SHA-256, package, and artifact paths
// are recorded in doc/decisions/0004-its-mft-static-surface-spec-tables.md;
// testITSMFTSurfaceSpecProjection.cxx proves each literal here reproduces
// the exact source JSON token's float32 bit pattern.
//
// material reuses kNominalITSLayerX0 by reference, not by re-literaling
// it, so the reused behavioral baseline can never silently drift from that
// array (Gate 4 B1 Slice 2 report).
//
// SurfaceIds are dense and local to this spec (0..ITSNLayers-1); global
// rebasing across detectors, if ever needed, is ConcatenatedSurfaceSpec's
// job, not this file's.

#ifndef ALICEO2_ITSMFT_TRACKING_ITSSURFACESPEC_H_
#define ALICEO2_ITSMFT_TRACKING_ITSSURFACESPEC_H_

#include <array>
#include <cstddef>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/StaticSurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceSpec.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking
{

constexpr NominalSurfaceMaterial itsLayerMaterial(std::size_t layer) noexcept
{
  const float x0 = kNominalITSLayerX0[layer];
  return {x0, x0 * o2::its::constants::Radl * o2::its::constants::Rho};
}

struct ITSSurfaceSpec {
  inline static constexpr std::array<StaticSurfaceDescriptor, ITSNLayers> surfaces{
    StaticSurfaceDescriptor{SurfaceId{0}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 0}, SurfaceKind::Cylinder,
                            2.3259652f, itsLayerMaterial(0), SurfaceIndexingFamily::CylindricalPhiZ},
    StaticSurfaceDescriptor{SurfaceId{1}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 1}, SurfaceKind::Cylinder,
                            3.1353536f, itsLayerMaterial(1), SurfaceIndexingFamily::CylindricalPhiZ},
    StaticSurfaceDescriptor{SurfaceId{2}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 2}, SurfaceKind::Cylinder,
                            3.9162421f, itsLayerMaterial(2), SurfaceIndexingFamily::CylindricalPhiZ},
    StaticSurfaceDescriptor{SurfaceId{3}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 3}, SurfaceKind::Cylinder,
                            19.58824f, itsLayerMaterial(3), SurfaceIndexingFamily::CylindricalPhiZ},
    StaticSurfaceDescriptor{SurfaceId{4}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 4}, SurfaceKind::Cylinder,
                            24.527159f, itsLayerMaterial(4), SurfaceIndexingFamily::CylindricalPhiZ},
    StaticSurfaceDescriptor{SurfaceId{5}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 5}, SurfaceKind::Cylinder,
                            34.354595f, itsLayerMaterial(5), SurfaceIndexingFamily::CylindricalPhiZ},
    StaticSurfaceDescriptor{SurfaceId{6}, {static_cast<uint8_t>(o2::detectors::DetID::ITS), 6}, SurfaceKind::Cylinder,
                            39.310642f, itsLayerMaterial(6), SurfaceIndexingFamily::CylindricalPhiZ},
  };
};

static_assert(SurfaceSpec<ITSSurfaceSpec>);
static_assert(SurfaceCount<ITSSurfaceSpec> == ITSNLayers);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_ITSSURFACESPEC_H_ */
