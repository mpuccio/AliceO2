// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B2 Slice 1: constexpr, static-storage-duration projections of the
// candidate ITSSurfaceSpec/MFTSurfaceSpec tables (Gate 4 B1 Slice C3a) into
// runtime SurfaceDescriptor catalogs, plus their compile-time concatenation
// into one combined ITS+MFT global id space.
//
// Existing mechanism only: toRuntimeSurfaceDescriptor() (already used by
// testITSMFTSurfaceSpecProjection.cxx per-surface) and ConcatenatedSurfaceSpec
// (SurfaceSpec.h, already proven concatenable for ITSSurfaceSpec/MFTSurfaceSpec
// in that same test). No new catalog abstraction is introduced here.
//
// `inline constexpr` gives each array external linkage and static storage
// duration, so a SurfaceCatalogView (SurfaceCatalogView.h) borrowing a
// pointer into one of these arrays is valid for the process's entire
// lifetime -- never dangling, never requiring a copy.
//
// kITSStaticSurfaceCatalog/kMFTStaticSurfaceCatalog remain the sole plan
// source for single-detector production tracking as of Gate 4 B2 Slice 2:
// ITSMFTTrackingInterface::initialiseTracker() builds its one immutable
// DetectorLayoutSet from whichever of the two matches DetId, once, via
// buildDetectorLayoutSet() (DetectorLayoutSet.h). They are adapter/
// application data for that single-detector path, not a combined-tracking
// input.
//
// kITSMFTCombinedStaticSurfaceCatalog is, as of Gate 4 C2/C3, the
// authoritative combined-catalog source for combined disconnected tracking:
// ITSMFTLegacyParticipantSet (ITSMFTLegacyParticipantSet.h) builds its one
// 17-surface DetectorLayout directly from this catalog and derives both
// detectors' the retired traversal binding from that single build, so ITS and
// MFT traversal always share one identical global id space. Mixed-detector
// track publication stays rejected/diagnostic outside that set
// (TrackerTraits::initialiseTimeFrame()'s existing MixedPolicyLayout gate
// already covers that, unchanged).

#ifndef ALICEO2_ITSMFT_TRACKING_STATICDETECTORCATALOGS_H_
#define ALICEO2_ITSMFT_TRACKING_STATICDETECTORCATALOGS_H_

#include <array>
#include <cstddef>

#include "ITSMFTTracking/ITSSurfaceSpec.h"
#include "ITSMFTTracking/MFTSurfaceSpec.h"
#include "ITSMFTTracking/StaticSurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceSpec.h"
#include "ITSMFTTracking/TrackingConfigParam.h"

namespace o2::itsmft::tracking
{

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

#endif /* ALICEO2_ITSMFT_TRACKING_STATICDETECTORCATALOGS_H_ */
