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

#include <array>

#include <gsl/span>

#include "ITSMFTTracking/Cell.h"
#include "ITSMFTTracking/SurfaceMeasurement.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITStracking/Constants.h"

namespace o2::itsmft::tracking
{

/// Gate 3 Slice A (native ITS refit driver prerequisite): builds one
/// traversal-ordered leg of `SurfaceMeasurement` slots for `driveRefitLeg<Tag>`
/// from a `TrackSeedN<NLayers>`'s already-attached, legacy-layer-indexed
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
/// Per visited legacy layer, `seed.getCluster(legacyLayer)` (Cell.h) is either
/// `o2::its::constants::UnusedIndex` -- a hole, written as a default-
/// constructed `SurfaceMeasurement` (`cluster.isValid() == false`), matching
/// `driveRefitLeg`'s own hole contract exactly -- or a valid index into
/// `layerMeasurements[legacyLayer]`, the same authoritative, per-iteration
/// validated span `TrackerTraits::mLayerMeasurements` already is (its
/// correspondence to every attached cluster index is established once per
/// iteration by `initialiseTimeFrame()`'s `NormalizedMeasurementMismatch`
/// validation block, TrackerTraits.cxx, before any candidate is built). This
/// function trusts that established mapping and does not re-validate cluster
/// index bounds itself, exactly as `attachHit`/`buildCellSeed`'s own
/// production callers already do for the identical lookup.
///
/// Slots are written at `out[position++]` while iterating in the caller's
/// traversal order (`legacyLayer` running `start, start+step, ...`), never at
/// `out[legacyLayer]` -- the array index and the traversal order coincide only
/// for an increasing (`step == +1`) leg. Writing by `legacyLayer` instead of by
/// `position` would silently reorder a decreasing leg's slots back into
/// ascending-layer order, breaking `driveRefitLeg`'s hole/gate/direction
/// semantics for that leg without any crash or bounds violation to reveal it.
/// The returned span is the populated contiguous prefix `[0, position)`; for
/// every call this function currently expects (the frozen `refitTrack`
/// driver's three legs, each spanning the complete `[0, NLayers)` or
/// `[NLayers - 1, -1)` range), `position` ends at `NLayers`, but the caller
/// must use the returned span, not assume `NLayers` directly.
///
/// Unwired: no production call site uses this in this slice.
template <int NLayers>
gsl::span<const SurfaceMeasurement> assembleRefitLegSlots(
  const TrackSeedN<NLayers>& seed,
  const LayerMeasurementSpans<NLayers>& layerMeasurements,
  int start, int end, int step,
  std::array<SurfaceMeasurement, NLayers>& out) noexcept
{
  int position = 0;
  for (int legacyLayer = start; legacyLayer != end; legacyLayer += step) {
    const int clsIdx = seed.getCluster(legacyLayer);
    out[position++] = (clsIdx == o2::its::constants::UnusedIndex)
                        ? SurfaceMeasurement{}
                        : layerMeasurements[legacyLayer][clsIdx];
  }
  return gsl::span<const SurfaceMeasurement>(out.data(), position);
}

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_REFITLEGASSEMBLY_H_ */
