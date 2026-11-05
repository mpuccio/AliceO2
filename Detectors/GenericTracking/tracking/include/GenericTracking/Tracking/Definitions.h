// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Definitions.h
/// \brief Scalars, parameter indices and packed-covariance indexing for the
///        track model, the tracking's statuses and kinds, and the dispatch
///        and sort helpers the stages share.
///
/// Nothing here may depend on ROOT or on anything a device compiler cannot
/// see.
///
/// The enums the tracking's algorithms share. A status is compared where it
/// is read, never wrapped in a predicate. Enums of one piece of machinery stay
/// with it.

#ifndef GTRK_TRACKING_DEFINITIONS_H_
#define GTRK_TRACKING_DEFINITIONS_H_

#include "GenericTracking/Base/Definitions.h"

#include <cfloat>
#include <cstdint>
#include <type_traits>

namespace o2::gtrk
{

/// The scalar used for the track parameters, covariance and their arithmetic.
/// Time is not a track parameter: a measurement's and a track's time is a
/// box in integer picoseconds (`TimeEstPS`).
using Scalar = float;

/// The scalar used where the geometry's dynamic range must survive: surface
/// intersection and the path-length solve, which run at |z| up to 300 cm and
/// have to land inside a 1 um budget.  Results narrow to Scalar on return, so
/// this type never appears in a stored record.
using ExactScalar = double;

/// Indices into the bound parameter vector.
///
///   loc0, loc1  sensor-local coordinates on the reference plane   [cm]
///   phi, theta  global momentum direction                         [rad]
///   qOverP      signed charge over momentum                       [(GeV/c)^-1]
///
/// The order is part of the record layout and of every Jacobian written
/// against it; changing it is not a refactor.
enum BoundIndex : std::uint8_t { kLoc0 = 0,
                                 kLoc1 = 1,
                                 kPhi = 2,
                                 kTheta = 3,
                                 kQOverP = 4,
                                 kNBoundParams = 5 };

/// A symmetric 5x5 needs 15 numbers, not 25.
inline constexpr int kBoundCovSize = kNBoundParams * (kNBoundParams + 1) / 2;

/// The strict upper triangle of a 5x5, the U of a U D U^T factorisation.
inline constexpr int kBoundUpperSize = kBoundCovSize - kNBoundParams;

/// Index of element (row, col) in the packed lower triangle, row-major:
///
///     0
///     1  2
///     3  4  5
///     6  7  8  9
///    10 11 12 13 14
///
/// This is the layout O2 uses for its 5x5: o2::track::CovarMap is
/// {{0,1,3,6,10},...} with DiagMap {0,2,5,9,14}, which is what this formula
/// produces.  O2 tabulates it and we compute it.
///
/// Symmetric in its arguments, so a caller never has to know which of the two
/// indices is the larger one -- which is exactly the mistake that makes a
/// packed covariance silently wrong rather than loudly broken.
GPUhdi() GPUAlwaysInline() constexpr int covIndex(int row, int col)
{
  const int hi = row > col ? row : col;
  const int lo = row > col ? col : row;
  return (hi * (hi + 1)) / 2 + lo;
}

/// Index of the variance of parameter \a i.
GPUhdi() GPUAlwaysInline() constexpr int varIndex(int i) { return covIndex(i, i); }

// Propagation and fits.

/// What propagation did, in place of an exception, which device code cannot
/// throw.  A caller has to look at this: several of these are ordinary
/// physics, not faults.
enum class PropagationStatus : std::uint8_t {
  kOk = 0,
  /// The trajectory never reaches the target plane.  A track whose transverse
  /// circle closes inside the target radius is the common case, and it is a
  /// fact about the track rather than an error.
  kNoIntersection,
  /// The plane is only reachable against the requested direction of travel.
  kWrongDirection,
  /// The intersection lies outside the sensor's active area.  Reported, not
  /// enforced -- a fit may want the point anyway, a search will not.
  kOutOfBounds,
  /// The track grazes the plane: |normal . direction| is below threshold, so
  /// the bound representation is singular and the covariance would be
  /// meaningless.  Real at large |eta| on the outer barrel.
  kGrazing,
  /// A non-finite parameter arrived or was produced.
  kNotFinite,
  /// The crossing lies more than half a turn of the helix along the track: a
  /// looper coming round again, not the next layer of the same arc.
  kTooManyTurns,
  /// The mass hypothesis's energy loss stopped the particle in the step's
  /// material, or left it below 10 MeV (`energyLoss`).
  kStopped
};

/// Which way along the trajectory the target is expected to lie.
enum class PropagationDirection : std::int8_t {
  kForward = 1,  ///< along the momentum
  kBackward = -1 ///< against it, as the inward refit needs
};

/// Which of a barrel plane's two crossings a propagation takes
/// (`propagate`). A compile-time choice
/// (`withCrossingRule`); `kFirstAhead` is the default.
enum class CrossingRule : std::uint8_t {
  /// The crossing of smaller offset along the plane's line whenever it lies
  /// ahead; the other only where it does not, and only nearer the centre.
  kFirstAhead = 0,
  /// Of the crossings ahead, the one nearer the plane's centre; a tie to the
  /// first. With the first behind, as the default.
  kNearerCentre,
  /// Of the crossings ahead, one inside the sensor's bounds, the first if
  /// both are; where neither is, the one nearer the plane's centre, a tie to
  /// the first. With the first behind, as the default.
  kInsideFirst
};

/// Where a failed fit stopped: the outward filter or the inward one.
enum class FitPass : std::uint8_t { kNone = 0,
                                    kOutward,
                                    kInward };

/// The road fit's outcome (TrackFitter.h, `fitOutward`): accepted, refused by
/// its cuts, or a propagation failed on the way.
enum class RoadFitStatus : std::uint8_t { kAccepted = 0,
                                          kRejected,
                                          kPropagationFailed };

// The stages.

/// Free; shared, taken by tracks only at positions among their `k` innermost;
/// exclusive, taken at any other. Written by the selection alone, between
/// iterations; read by the next iteration's tracklet stage, which skips the
/// exclusive rows.
enum UsedState : std::uint8_t { kFree = 0U,
                                kShared = 1U,
                                kExclusive = 2U };

/// Where an iteration's vertex seeds come from.
enum class VertexSource : std::uint8_t {
  kBeamLine,   ///< one seed, the beam spot, compatible with every time
  kVertexList, ///< the timeframe's vertex seeds, sorted by time
};

/// How a pair's second-coordinate window follows from a seed.
enum class PairKind : std::uint8_t {
  kBarrel,       ///< barrel -> barrel: linear in the seed's z
  kDisc,         ///< disc -> disc: linear in the seed's r
  kBarrelToDisc, ///< evaluated per seed, see `crossingPairWindow`
  kDiscToBarrel  ///< evaluated per seed, see `crossingPairWindow`
};

/// What a step did with a row it visited, in the order its tests run: held
/// by a track, outside the step's window, crossing another layer
/// (LinkHoles.h), no usable update, above the chi2 gate, or passing.
enum class StepOutcome : std::uint8_t { kUsed,
                                        kOutsideWindow,
                                        kHole,
                                        kNoUpdate,
                                        kChi2,
                                        kPasses };

/// Why a direction of a follow ended: no layer to search from where it
/// stands; a sensor crossed without a row past the iteration's holes; every
/// layer searched without a row within the holes (crossed in gaps, or not
/// reached); or the track as long as a chain can be.
enum class FollowEnd : std::uint8_t { kNoTarget,
                                      kHoles,
                                      kNoRow,
                                      kCapacity };

/// What the map makes of one triplet, and why it is not a candidate.
enum class VertexCandidateStatus : std::uint8_t { kAccepted,
                                                  kChi2RZ,
                                                  kDca,
                                                  kDegenerate };

/// Whether a value is finite, in every build: O2's `GPUCommonMath::Finite`
/// and `IsNaN` are constants outside the deterministic mode, on the host as
/// on the device. NaN fails both comparisons, an
/// infinity one of them. The bounds are the C macros, which a device compiler
/// accepts where a `std::numeric_limits` call may not.
GPUhdi() GPUAlwaysInline() constexpr bool isFiniteScalar(float x) { return x >= -FLT_MAX && x <= FLT_MAX; }

/// `body(std::bool_constant<b>{}...)` for the run-time flags: an
/// iteration's choices, each a compile-time variant of a stage, dispatched
/// once per combination; the host's and a device's one switch.
template <typename Body>
decltype(auto) withBools(const Body& body)
{
  return body();
}
template <typename Body, typename... Rest>
decltype(auto) withBools(const Body& body, bool first, Rest... rest)
{
  if (first) {
    return withBools([&](auto... fixed) { return body(std::true_type{}, fixed...); }, rest...);
  }
  return withBools([&](auto... fixed) { return body(std::false_type{}, fixed...); }, rest...);
}

/// `body(std::integral_constant<CrossingRule, rule>{})` for the run-time \a rule.
template <typename Body>
decltype(auto) withCrossingRule(CrossingRule rule, const Body& body)
{
  switch (rule) {
    case CrossingRule::kNearerCentre:
      return body(std::integral_constant<CrossingRule, CrossingRule::kNearerCentre>{});
    case CrossingRule::kInsideFirst:
      return body(std::integral_constant<CrossingRule, CrossingRule::kInsideFirst>{});
    case CrossingRule::kFirstAhead:
      break;
  }
  return body(std::integral_constant<CrossingRule, CrossingRule::kFirstAhead>{});
}

/// \a count values sorted in place by \a less, stable: an insertion sort,
/// for the few entries one thread holds. Stable, and the same order on host
/// and device; O2's `sort` is neither.
template <typename T, typename Less>
GPUhdi() GPUAlwaysInline() void insertionSort(T* values, std::uint32_t count, const Less& less)
{
  for (std::uint32_t i = 1; i < count; ++i) {
    const T value = values[i];
    std::uint32_t slot = i;
    while (slot > 0U && less(value, values[slot - 1U])) {
      values[slot] = values[slot - 1U];
      --slot;
    }
    values[slot] = value;
  }
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_DEFINITIONS_H_
