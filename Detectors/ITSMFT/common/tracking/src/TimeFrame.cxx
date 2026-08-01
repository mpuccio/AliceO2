// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file TimeFrame.cxx
/// \brief
///

#include "ITSMFTTracking/TimeFrame.h"

namespace o2::itsmft::tracking
{

using o2::its::deepVectorClear;

void TimeFrame::addPrimaryVertex(const Vertex& vert)
{
  mPrimaryVertices.emplace_back(vert);
  if (!isBeamPositionOverridden) {
    const float w = vert.getNContributors();
    mBeamPos[0] = (mBeamPos[0] * mBeamPosWeight + vert.getX() * w) / (mBeamPosWeight + w);
    mBeamPos[1] = (mBeamPos[1] * mBeamPosWeight + vert.getY() * w) / (mBeamPosWeight + w);
    mBeamPosWeight += w;
  }
}

void TimeFrame::resetBeamXY(const float x, const float y, const float w)
{
  mBeamPos[0] = x;
  mBeamPos[1] = y;
  mBeamPosWeight = w;
}

void TimeFrame::commitNormalizedFrame(MultiSourceFrame&& staged) noexcept
{
  mNormalizedFrame = std::move(staged);
  // Gate 4 CommonTrack foundation: mNormalizedFrame above has just been
  // replaced, so any CommonTrack/TrackClusterReference built against the
  // previous normalized frame no longer refers to anything meaningful and
  // must not survive this commit -- cleared here, in the same successful
  // commit, exactly like every other member the owner-level load operation
  // replaces.
  deepVectorClear(mCommonTracks);
  deepVectorClear(mTrackClusterIndices);
}

void TimeFrame::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
{
  mMemoryPool = pool;

  auto initVector = [&]<typename T>(bounded_vector<T>& vec) {
    deepVectorClear(vec, mMemoryPool.get());
  };

  initVector(mPrimaryVertices);
  initVector(mPrimaryVerticesLabels);
  initVector(mCommonTracks);
  initVector(mTrackClusterIndices);
}

void TimeFrame::wipe()
{
  deepVectorClear(mPrimaryVertices);
  deepVectorClear(mPrimaryVerticesLabels);
  // Gate 4 CommonTrack foundation: common-CA result storage is per-event
  // data, invalidated together -- a CommonTrack's cluster-reference range
  // and every SurfaceMeasurementIndex it reaches are only meaningful
  // alongside the normalized frame current when it was built (CommonTrack.h),
  // so both must be cleared here, unconditionally.
  deepVectorClear(mCommonTracks);
  deepVectorClear(mTrackClusterIndices);
  // Event-owned normalized data, cleared unconditionally (MultiSourceFrame is
  // host-only and never framework/GPU-managed). This clears mNormalizedFrame
  // in place: a getNormalizedFrame() reference obtained before this call
  // still refers to that same live member and remains safe to use, now
  // observing its cleared (empty) state. Any MultiSourceFrameView or
  // gsl::span obtained before this call (from getNormalizedFrameView(),
  // getSurfaceMeasurements(), getSourceIntervals(), getLabels()) is
  // invalidated by it, since clear() may reallocate/free the buffers those
  // point into.
  mNormalizedFrame.clear();
}

} // namespace o2::itsmft::tracking
