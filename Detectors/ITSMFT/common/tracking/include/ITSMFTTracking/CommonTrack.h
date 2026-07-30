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

#include "DataFormatsITS/TimeEstBC.h"
#include "GPUCommonDef.h"
#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ITSMFTTracking/SurfaceMask.h"

namespace o2::itsmft::tracking
{

// Detector-neutral common-CA tracking result (Architecture.md Sec 12: "The
// TimeFrame stores a generic internal result"). Owned by the common
// TimeFrame (see TimeFrame::getCommonTracks()/getTrackClusterIndices()).
// Deliberately carries no ITS/MFT/ALICE-3 publication type, detector ID,
// NLayers template parameter, or workflow dependency: detector output
// adapters (TrackITSExt, TrackMFT, a future ALICE-3 type, ...) are built
// from this type by code outside the common core, never the reverse.
//
// Cluster membership is not stored inline. firstClusterRef/clusterRefEnd is
// a half-open [firstClusterRef, clusterRefEnd) range of *positions* into the
// TimeFrame-owned flat trackClusterIndices array (see
// TimeFrame::getTrackClusterIndices()); see isValidTrackRange() below for the
// exact validity condition. Each element of that array is, in turn, a
// SurfaceMeasurementIndex -- a canonical position into the flattened,
// normalized SurfaceMeasurement array owned by MultiSourceFrame
// (MultiSourceFrame::getMeasurement()/MultiSourceFrameView::getMeasurement())
// -- never a detector-local layer index or a raw external cluster index.
// References are stored in traversal order, inner to outer.
//
// innerState/outerState reuse SurfaceKinematicState directly (Barrel or
// Forward family, selected by SurfaceKinematicState::family): this is
// deliberately not a second five-parameter/covariance representation.
// hitSurfaces uses the global 32-bit SurfaceMask (ITSMFTTracking/
// SurfaceMask.h), never the legacy 16-bit per-NLayers LayerMask. For a
// valid, completed track, hitSurfaces is the union of the SurfaceId of every
// measurement referenced by [firstClusterRef, clusterRefEnd) -- this is a
// consumer-side invariant (see testCommonTrack.cxx), not something this
// struct enforces by construction, since populating a CommonTrack from CA
// seeds is out of scope for this slice.
//
// Lifetime: track-cluster-index storage and CommonTrack storage are
// TimeFrame event data, invalidated together by TimeFrame::wipe() (and by
// any subsequent reload), exactly like every other per-event CA artefact
// (mCells, mTracks, ...). A CommonTrack's firstClusterRef/clusterRefEnd
// range, and every SurfaceMeasurementIndex it reaches through
// trackClusterIndices, are only meaningful together with the TimeFrame's
// normalized frame (MultiSourceFrame) that was current when the track was
// built; none of the three are individually meaningful once any of the
// other two has been wiped or reloaded.
struct CommonTrack {
  SurfaceKinematicState innerState{};
  SurfaceKinematicState outerState{};
  float chi2{0.f};
  o2::its::TimeEstBC timestamp{};
  SurfaceMask hitSurfaces{};
  uint32_t firstClusterRef{0};
  uint32_t clusterRefEnd{0};
};

// Not standard-layout: o2::its::TimeEstBC's own base hierarchy
// (TimeStampWithError<uint32_t,uint16_t> -> TimeStamp<uint32_t>) declares
// non-static data members at more than one level, which by itself makes any
// aggregate containing it non-standard-layout -- regardless of CommonTrack's
// own (single-level, non-virtual, all-public) structure. This does not
// affect device/host compatibility in practice: CommonTrack has no virtual
// functions or bases of its own and every member (including TimeEstBC) is
// trivially copyable, which is what device views and bounded_vector storage
// actually require.
static_assert(!std::is_standard_layout_v<o2::its::TimeEstBC>);
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
