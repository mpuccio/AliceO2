// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_GENERICTRACK_H_
#define ALICEO2_ITSMFT_TRACKING_GENERICTRACK_H_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "GPUCommonDef.h"
#ifndef GPUCA_GPUCODE
#include "ITSMFTTracking/Cell.h"
#endif
#include "ITSMFTTracking/IdTypes.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMask.h"
#include "ITSMFTTracking/SurfaceTiming.h"

namespace o2::itsmft::tracking
{

// A measurement index is local to its surface. Resolve a reference through
// TimeFrame::getGlobalMeasurement(surface, index); raw external cluster IDs
// remain source-qualified ClusterRef data.
struct TrackClusterReference {
  SurfaceId surface{};
  SurfaceMeasurementIndex index{};
};

static_assert(std::is_standard_layout_v<TrackClusterReference>);
static_assert(std::is_trivially_copyable_v<TrackClusterReference>);
static_assert(sizeof(TrackClusterReference) == 8);
static_assert(alignof(TrackClusterReference) == 4);
static_assert(offsetof(TrackClusterReference, surface) == 0);
static_assert(offsetof(TrackClusterReference, index) == 4);

// Detector-neutral result owned by TimeFrame. Its reference range is
// inner-to-outer in the frame-owned trackClusterIndices array; hitSurfaces is
// their union. The result and its references are valid only with the same
// normalized frame and are cleared atomically on event reset or replacement.
struct GenericTrack {
  SurfaceKinematicState innerState{};
  SurfaceKinematicState outerState{};
  float chi2{0.f};
  GenericTrackTimestamp timestamp{};
  SurfaceMask hitSurfaces{};
  uint32_t firstClusterRef{0};
  uint32_t clusterRefEnd{0};
};

#ifndef GPUCA_GPUCODE

// Detector-neutral result of one successful seed refit. Typed output tracks
// and publication sidecars remain outside the common tracker.
struct TrackingCandidate {
  TrackSeed seed;
  GenericTrack track{};
  float phi{0.f};
  float eta{0.f};
  double charge{0.};
  uint32_t genericTrackIndex{std::numeric_limits<uint32_t>::max()};

  int getNumberOfClusters() const noexcept { return seed.getActiveSurfaceCount(); }
  int getClusterIndex(int position) const noexcept { return seed.getCluster(position); }
  int getFirstClusterLayer() const noexcept { return seed.getSurfaceMask().first(); }
};

#endif

// Both properties are required for a portable device-facing value: bytes are
// copyable and member order is stable for host/device compilation.
static_assert(std::is_standard_layout_v<GenericTrack>);
static_assert(std::is_trivially_copyable_v<GenericTrack>);

// A track range is valid iff 0 <= firstClusterRef <= clusterRefEnd <=
// trackClusterIndicesSize (firstClusterRef's lower bound is automatic: both
// fields are unsigned). `trackClusterIndicesSize` is the caller's current
// TimeFrame::getTrackClusterIndices().size(), never inferred from the track
// itself.
GPUhdi() constexpr bool isValidTrackRange(const GenericTrack& track, uint32_t trackClusterIndicesSize) noexcept
{
  return track.firstClusterRef <= track.clusterRefEnd && track.clusterRefEnd <= trackClusterIndicesSize;
}

GPUhdi() constexpr uint32_t trackClusterRefCount(const GenericTrack& track) noexcept
{
  return track.clusterRefEnd - track.firstClusterRef;
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_GENERICTRACK_H_ */
