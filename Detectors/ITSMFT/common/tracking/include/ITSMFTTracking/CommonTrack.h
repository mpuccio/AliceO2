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
#include <type_traits>

#include "GPUCommonDef.h"
#include "ITSMFTTracking/SurfaceId.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMask.h"
#include "ITSMFTTracking/SurfaceMeasurementIndex.h"
#include "ITSMFTTracking/SurfaceTiming.h"

namespace o2::itsmft::tracking
{

// Detector-neutral element of the TimeFrame-owned flat track-membership
// array (TimeFrame::getTrackClusterIndices()). `surface` identifies which
// surface's own measurement array `index` is a position into -- measurements
// are owned per surface (MultiSourceFrame/MultiSourceFrameView), never as one
// flattened global array, so a bare SurfaceMeasurementIndex alone is never a
// complete identity: it must always be paired with the SurfaceId it is local
// to. A reference resolves through exactly one call:
//   normalizedFrame.getMeasurement(reference.surface, reference.index)
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

// Detector-neutral common-CA tracking result (Architecture.md Sec 12: "The
// TimeFrame stores a generic internal result"). Owned by the common
// TimeFrame (see TimeFrame::getCommonTracks()/getTrackClusterIndices()).
// Deliberately carries no ITS/MFT/ALICE-3 publication type, detector ID,
// compile-time layer parameter, or workflow dependency: detector output
// adapters (TrackITSExt, TrackMFT, a future ALICE-3 type, ...) are built
// from this type by code outside the common core, never the reverse.
//
// Cluster membership is not stored inline. firstClusterRef/clusterRefEnd is
// a half-open [firstClusterRef, clusterRefEnd) range of *positions* into the
// TimeFrame-owned flat trackClusterIndices array (see
// TimeFrame::getTrackClusterIndices(), an array of TrackClusterReference);
// see isValidTrackRange() below for the exact validity condition. Each
// referenced TrackClusterReference resolves through
// normalizedFrame.getMeasurement(reference.surface, reference.index) -- see
// that struct's own doc above. References are stored in traversal order,
// inner to outer.
//
// innerState/outerState reuse SurfaceKinematicState directly (Barrel or
// Forward family, selected by SurfaceKinematicState::family): this is
// deliberately not a second five-parameter/covariance representation.
// hitSurfaces uses the global 32-bit SurfaceMask (ITSMFTTracking/
// SurfaceMask.h), never a legacy fixed-layer LayerMask. For a valid,
// completed track: every TrackClusterReference in
// [firstClusterRef, clusterRefEnd) resolves to an existing measurement whose
// own SurfaceMeasurement::surface equals that reference's surface field, and
// hitSurfaces is the union of those surfaces. This is a consumer-side
// invariant (see testCommonTrack.cxx), not something this struct enforces by
// construction, since populating a CommonTrack from CA seeds is out of scope
// for this slice.
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

// Standard-layout and trivially-copyable together -- not trivial-copyability
// alone. Trivial-copyability only guarantees CommonTrack's bytes may be
// copied (e.g. memcpy'd to a device buffer) without invoking non-trivial
// special member functions; it says nothing about a consistent, portable
// field layout. Standard-layout is the property that additionally
// guarantees a single, well-defined member order usable with `offsetof` and
// consistent across host and device compilation -- both properties are
// required together for the same reason every other device-facing type in
// this library asserts both (SurfaceMeasurement, StaticSurfaceDescriptor,
// SurfaceGraphView, ...), and neither one substitutes for the other.
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
