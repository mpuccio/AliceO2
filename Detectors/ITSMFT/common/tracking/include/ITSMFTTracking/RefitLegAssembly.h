// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_REFITLEGASSEMBLY_H_
#define ALICEO2_ITSMFT_TRACKING_REFITLEGASSEMBLY_H_

#include "GPUCommonDef.h"

#ifndef GPUCA_GPUCODE

#include <gsl/span>

#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking
{

/// Gate 3 Slice A (native ITS refit driver prerequisite): builds one
/// traversal-ordered leg of `SurfaceMeasurement` slots for `driveRefitLeg<Tag>`
/// from a `TrackSeed`'s already-attached, layer-indexed
/// cluster bookkeeping -- exactly the mapping `computeLayerCellsForPolicy`/
/// `processNeighbours` already use per candidate (TrackerTraits.cxx), walked
/// once here across a whole seed's leg instead of once per candidate.
///
/// `[start, end)` stepping by `step` is the same legacy layer-index range the
/// frozen `o2::its::track::fitTrack` loop iterates (ITStracking/TrackHelpers.h):
/// `step == +1` for an inward-index (increasing-layer) leg, `step == -1` for an
/// outward-index (decreasing-layer) leg. This function does not choose or
/// validate that range; the caller (the future native refit driver) supplies
/// it exactly as the frozen `refitTrack` driver does for each of its legs.
///
/// Per visited plan position, `seed.getCluster(position)` (Cell.h) is either
/// `o2::its::constants::UnusedIndex` -- a hole, written as a default-
/// constructed `SurfaceMeasurement` (`cluster.isValid() == false`), matching
/// `driveRefitLeg`'s own hole contract exactly -- or a valid index into
/// `layerMeasurements[position]`, the same authoritative, per-iteration
/// validated span `TrackerTraits::mLayerMeasurements` already is (its
/// correspondence to every attached cluster index is established once per
/// iteration by `initialiseTimeFrame()`'s `NormalizedMeasurementMismatch`
/// validation block, TrackerTraits.cxx, before any candidate is built). This
/// function trusts that established mapping and does not re-validate cluster
/// index bounds itself, exactly as `attachHit`/`buildCellSeed`'s own
/// production callers already do for the identical lookup.
///
/// Slots are written at `out[position++]` while iterating in the caller's
/// traversal order (`position` running `start, start+step, ...`), never at
/// `out[position]` from the source coordinate -- the array index and traversal
/// order coincide only for an increasing (`step == +1`) leg. Writing by source
/// `position` would silently reorder a decreasing leg's slots back into
/// ascending-layer order, breaking `driveRefitLeg`'s hole/gate/direction
/// semantics for that leg without any crash or bounds violation to reveal it.
/// The returned span is the populated contiguous prefix `[0, position)`; for
/// every call this function currently expects (the frozen `refitTrack`
/// driver's three legs, each spanning the complete runtime surface range),
/// the returned span is the populated prefix and the caller must use it.
///
/// Unwired: no production call site uses this in this slice.
inline gsl::span<const SurfaceMeasurement> assembleRefitLegSlots(
  const TrackSeed& seed,
  gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
  int start, int end, int step,
  gsl::span<SurfaceMeasurement> out) noexcept
{
  int position = 0;
  for (int surfacePosition = start; surfacePosition != end && position < static_cast<int>(out.size()); surfacePosition += step) {
    const int clsIdx = seed.getCluster(surfacePosition);
    out[position++] = (clsIdx == o2::its::constants::UnusedIndex)
                        ? SurfaceMeasurement{}
                        : layerMeasurements[surfacePosition][clsIdx];
  }
  return gsl::span<const SurfaceMeasurement>(out.data(), position);
}

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_REFITLEGASSEMBLY_H_ */
