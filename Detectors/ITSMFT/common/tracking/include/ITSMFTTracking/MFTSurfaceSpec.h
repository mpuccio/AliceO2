// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B1 Slice C3a: candidate static MFT DetectorSpec table.
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
// material reuses kNominalMFTLayerX0 by reference, not by re-literaling
// it, so the reused behavioral baseline can never silently drift from that
// array (Gate 4 B1 Slice 2 report).
//
// SurfaceIds are dense and local to this spec (0..MFTNLayers-1); global
// rebasing across detectors, if ever needed, is ConcatenatedSurfaceSpec's
// job, not this file's.

#ifndef ALICEO2_ITSMFT_TRACKING_MFTSURFACESPEC_H_
#define ALICEO2_ITSMFT_TRACKING_MFTSURFACESPEC_H_

#include <array>
#include <cstddef>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/NominalSurfaceMaterialDefaults.h"
#include "ITSMFTTracking/StaticSurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceSpec.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking
{

constexpr NominalSurfaceMaterial mftLayerMaterial(std::size_t layer) noexcept
{
  const float x0 = kNominalMFTLayerX0[layer];
  return {x0, x0 * o2::its::constants::Radl * o2::its::constants::Rho};
}

struct MFTSurfaceSpec {
  inline static constexpr std::array<StaticSurfaceDescriptor, MFTNLayers> surfaces{
    StaticSurfaceDescriptor{SurfaceId{0}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 0}, SurfaceKind::Disk,
                            -45.2889f, mftLayerMaterial(0), SurfaceIndexingFamily::CartesianXY},
    StaticSurfaceDescriptor{SurfaceId{1}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 1}, SurfaceKind::Disk,
                            -46.7111f, mftLayerMaterial(1), SurfaceIndexingFamily::CartesianXY},
    StaticSurfaceDescriptor{SurfaceId{2}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 2}, SurfaceKind::Disk,
                            -48.5889f, mftLayerMaterial(2), SurfaceIndexingFamily::CartesianXY},
    StaticSurfaceDescriptor{SurfaceId{3}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 3}, SurfaceKind::Disk,
                            -50.0111f, mftLayerMaterial(3), SurfaceIndexingFamily::CartesianXY},
    StaticSurfaceDescriptor{SurfaceId{4}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 4}, SurfaceKind::Disk,
                            -52.3889f, mftLayerMaterial(4), SurfaceIndexingFamily::CartesianXY},
    StaticSurfaceDescriptor{SurfaceId{5}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 5}, SurfaceKind::Disk,
                            -53.8111f, mftLayerMaterial(5), SurfaceIndexingFamily::CartesianXY},
    StaticSurfaceDescriptor{SurfaceId{6}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 6}, SurfaceKind::Disk,
                            -67.6889f, mftLayerMaterial(6), SurfaceIndexingFamily::CartesianXY},
    StaticSurfaceDescriptor{SurfaceId{7}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 7}, SurfaceKind::Disk,
                            -69.1111f, mftLayerMaterial(7), SurfaceIndexingFamily::CartesianXY},
    StaticSurfaceDescriptor{SurfaceId{8}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 8}, SurfaceKind::Disk,
                            -76.0889f, mftLayerMaterial(8), SurfaceIndexingFamily::CartesianXY},
    StaticSurfaceDescriptor{SurfaceId{9}, {static_cast<uint8_t>(o2::detectors::DetID::MFT), 9}, SurfaceKind::Disk,
                            -77.5111f, mftLayerMaterial(9), SurfaceIndexingFamily::CartesianXY},
  };
};

static_assert(SurfaceSpec<MFTSurfaceSpec>);
static_assert(SurfaceCount<MFTSurfaceSpec> == MFTNLayers);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_MFTSURFACESPEC_H_ */
