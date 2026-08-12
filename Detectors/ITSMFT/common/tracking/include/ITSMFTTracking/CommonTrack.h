// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_COMMONTRACK_H_
#define ALICEO2_ITSMFT_TRACKING_COMMONTRACK_H_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "GPUCommonDef.h"
#ifndef GPUCA_GPUCODE
#include "ITSMFTTracking/Cell.h"
#endif
#include "ITSMFTTracking/SurfaceId.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMask.h"
#include "ITSMFTTracking/SurfaceTiming.h"

namespace o2::itsmft::tracking
{

// Detector-neutral element of the TimeFrame-owned flat track-membership
// array (TimeFrame::getTrackClusterIndices()). `surface` identifies which
// surface's own measurement array `index` is a position into -- measurements
// are owned per surface by TimeFrame, never as one
// flattened global array, so a bare SurfaceMeasurementIndex alone is never a
// complete identity: it must always be paired with the SurfaceId it is local
// to. A reference resolves through exactly one call:
//   normalizedFrame.getGlobalMeasurement(reference.surface, reference.index)
// Never a detector-local raw cluster index or external cluster ID: those
// remain ClusterRef's job (source-qualified, ITSMFTTracking/
// SurfaceMeasurement.h), a different identity axis entirely.
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

// Detector-neutral common-CA result owned by the TimeFrame. Publication
// adapters remain outside the tracking core.
//
// Cluster membership is a half-open range of positions in the TimeFrame-owned
// trackClusterIndices array. References resolve through the normalized frame
// and are stored in inner-to-outer traversal order.
//
// innerState/outerState use the family-selected SurfaceKinematicState;
// hitSurfaces is the global 32-bit SurfaceMask. A completed track resolves
// every reference in the range and has their union in hitSurfaces.
//
// timestamp is CommonTrackTimestamp (ITSMFTTracking/SurfaceTiming.h), the
// common TFBC-based half-open BC interval -- not o2::its::TimeEstBC, which
// is not standard-layout (its TimeStampWithError/TimeStamp base hierarchy
// declares non-static data members at more than one level) and is an
// ITS-namespaced type despite being reused elsewhere in this library.
//
// Lifetime: track-cluster-reference storage and CommonTrack storage are
// TimeFrame event data. Both are cleared together with every other
// per-event CA artefact by TimeFrame::resetEvent(), and are also cleared together,
// in the same successful commit, whenever TimeFrame::loadNormalizedSource()
// replaces the normalized frame they were built against (a stale
// CommonTrack/trackClusterIndices pair referencing measurements from a
// now-replaced normalized frame is not meaningful and must not survive a
// reload). A *failed* loadNormalizedSource() call leaves the normalized
// frame, CommonTrack storage and trackClusterIndices all completely
// unchanged, matching that call's existing transactional contract for every
// other TimeFrame member. A CommonTrack's range, and every
// TrackClusterReference it reaches, are only meaningful together with the
// TimeFrame's normalized frame that was current when the track was built;
// none of the three is individually meaningful once any of the other two
// has been wiped or reloaded.
struct CommonTrack {
  SurfaceKinematicState innerState{};
  SurfaceKinematicState outerState{};
  float chi2{0.f};
  CommonTrackTimestamp timestamp{};
  SurfaceMask hitSurfaces{};
  uint32_t firstClusterRef{0};
  uint32_t clusterRefEnd{0};
};

#ifndef GPUCA_GPUCODE

// Detector-neutral result of one successful seed refit. Typed output tracks
// and publication sidecars remain outside the common tracker.
struct TrackingCandidate {
  TrackSeed seed;
  CommonTrack track{};
  float phi{0.f};
  float eta{0.f};
  double charge{0.};
  uint32_t commonTrackIndex{std::numeric_limits<uint32_t>::max()};

  int getNumberOfClusters() const noexcept { return seed.getActiveSurfaceCount(); }
  int getClusterIndex(int position) const noexcept { return seed.getCluster(position); }
  int getFirstClusterLayer() const noexcept { return seed.getSurfaceMask().first(); }
};

#endif

// Both properties are required for a portable device-facing value: bytes are
// copyable and member order is stable for host/device compilation.
static_assert(std::is_standard_layout_v<CommonTrack>);
static_assert(std::is_trivially_copyable_v<CommonTrack>);

// A track range is valid iff 0 <= firstClusterRef <= clusterRefEnd <=
// trackClusterIndicesSize (firstClusterRef's lower bound is automatic: both
// fields are unsigned). `trackClusterIndicesSize` is the caller's current
// TimeFrame::getTrackClusterIndices().size(), never inferred from the track
// itself.
GPUhdi() constexpr bool isValidTrackRange(const CommonTrack& track, uint32_t trackClusterIndicesSize) noexcept
{
  return track.firstClusterRef <= track.clusterRefEnd && track.clusterRefEnd <= trackClusterIndicesSize;
}

GPUhdi() constexpr uint32_t trackClusterRefCount(const CommonTrack& track) noexcept
{
  return track.clusterRefEnd - track.firstClusterRef;
}

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_COMMONTRACK_H_ */
