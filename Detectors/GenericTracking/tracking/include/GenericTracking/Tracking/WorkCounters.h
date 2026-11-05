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

/// \file WorkCounters.h
/// \brief Counts of the work inside the stages' per-element loops.
///
/// The candidates a stage looks at, what each cut
/// rejects, the fits it reaches and the time they take.
/// A count touches a hot loop, so only a build with `GTRK_WORK_COUNTERS`
/// counts (O2's CMake option of that name, which defines it for every unit
/// that links the tracking), and only on the host: elsewhere `countWork` and
/// `WorkTimer` are empty. Each thread counts into its own block; `workCounts` sums them, and
/// is read between parallel loops, whose join orders the counts before it.
/// Nothing here reads or writes what the tracking computes.

#ifndef GTRK_TRACKING_WORKCOUNTERS_H_
#define GTRK_TRACKING_WORKCOUNTERS_H_

#include "GenericTracking/Base/Definitions.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#if defined(GTRK_WORK_COUNTERS) && !defined(__CUDACC__) && !defined(__HIPCC__)
#include <chrono>
#endif

namespace o2::gtrk
{

/// Every count, in report order.
enum class WorkCounter : std::uint8_t {
  kTripletFirsts,       ///< first tracklets visited
  kTripletNoMiddle,     ///< of them, whose middle is not a source of the slice
  kTripletWindow,       ///< second tracklets inside the dip-angle window
  kTripletTime,         ///< of them, rejected by the time cut
  kTripletCurvature,    ///< rejected by the curvature cut
  kTripletFits,         ///< reaching the fit
  kTripletFitFailed,    ///< the fit failed
  kTripletChi2,         ///< rejected by the fit's chi2
  kTripletOwnMomentum,  ///< rejected by the fit's chi2 with the scattering at the triplet's own momentum
  kTripletKinked,       ///< kept only with the kink allowance (counted in the rejections above as well, once each)
  kTripletFitNS,        ///< time inside the fit
  kNeighbourTriplets,   ///< triplets visited
  kNeighbourCandidates, ///< triplets ending on their first tracklet
  kNeighbourBothKinked, ///< of them, both triplets kinked: never neighbours
  kNeighbourCurvature,  ///< rejected by the curvature agreement
  kNeighbourTime,       ///< rejected by the time cut
  kNeighbourFits,       ///< reaching the fit
  kNeighbourFitFailed,  ///< the fit failed
  kNeighbourChi2,       ///< rejected by the fit's chi2
  kNeighbourFitNS,      ///< time inside the fit
  kRoadTime,            ///< rejected before the fit: their boxes share no time
  kRoadFits,            ///< roads fitted
  kRoadStepsPossible,   ///< the filter steps their complete outward passes would take
  kRoadSteps,           ///< the steps taken, the chi2 cut stopping a pass early
  kRoadChi2,            ///< rejected by the chi2 cut, early or at the end
  kRoadPropagation,     ///< the propagation failed
  kRoadFloor,           ///< below the pT floor for their cluster count
  kRoadBeam,            ///< refused by the beam line
  kRoadAccepted,        ///< selection records emitted
  kRoadTrimFits,        ///< roads fitted again without their innermost clusters (innerTrim)
  kRoadTrimAccepted,    ///< of the roads, the ones a fit without their innermost clusters took
  kTrackletSources,     ///< sources visited that no earlier track holds exclusively
  kTrackletQueries,     ///< their out-links whose window they reach: one grid query each
  kTrackletUsed,        ///< candidates an earlier track holds exclusively
  kTrackletAccepted,    ///< tracklets emitted
  kPairTests,           ///< candidates reaching the pair cut's transverse test
  kPairTransverse,      ///< of them, passing it
  kGridQueries,         ///< the tracklet stage's grid queries (the follow's are kFollowGridQueries)
  kGridBlocks,          ///< time blocks they span
  kGridRanges,          ///< row ranges they read, each a run of consecutive cells
  kGridCells,           ///< cells those runs cover
  kGridSearchedCells,   ///< of them, in runs dense enough for the per-cell binary searches
  kGridRows,            ///< rows in the runs
  kGridInBound,         ///< of them, inside the query's time bound
  kGridVisits,          ///< of those, passing the exact time cut: handed to the caller
  kPropagations,        ///< plane crossings solved (the fits' and the follow's)
  kPropagationsGrazing, ///< of them, at |n . d| <= 0.1: the ones a float solve would not serve
  kSortSegments,        ///< tracklet sort: sources with two tracklets or more
  kSortSegmentsOver8,   ///< of them, with more than 8
  kSortSegmentsOver32,  ///< with more than 32
  kSortSegmentsOver128, ///< with more than 128
  kSortSegmentsOver512, ///< with more than 512
  kSortTracklets,       ///< the tracklets of the sources counted in kSortSegments
  kSortTrackletsOver32, ///< of them, in sources with more than 32
  kPairLongitudinal,    ///< pair tests passing the longitudinal (second-coordinate) test, after the transverse one
  // The grid's counts of the follow's queries, in kGridQueries..kGridVisits's
  // order: the grid counts there while a FollowGridScope is open.
  kFollowGridQueries,
  kFollowGridBlocks,
  kFollowGridRanges,
  kFollowGridCells,
  kFollowGridSearchedCells,
  kFollowGridRows,
  kFollowGridInBound,
  kFollowGridVisits,
  kTrackletHoles,       ///< candidates whose segment crosses another layer's envelope, in an iteration taking no skips
  kTrackletHoleQueries, ///< queries not made: every pair of the source on the link would cross one
  kPairBetweenBound,    ///< candidates the transverse test refuses at the link's bound on the kink between the layers, before their own
  // The follow: its searches, its passes, and the
  // passes of the earlier iterations' follow apart.
  kFollowUpdates,              ///< its candidates reaching the update's chi2: free and in the step's window
  kFollowProposed,             ///< rows the follow proposes, before the bidding
  kFollowOutbid,               ///< of them, cut by the bidding: on or after a row another track won
  kFollowRefitRejected,        ///< of the rest, dropped with their track's failed refit
  kFollowKinkSteps,            ///< steps no row passed that tried the kink allowance
  kFollowKinksTaken,           ///< of them, taking a confirmed row
  kFollowLookaheads,           ///< steps more than one row passed that ranked them by their confirmation
  kFollowKinkUpdates,          ///< the kink passes' candidates reaching the update's chi2
  kFollowWitnessUpdates,       ///< the kink rows' confirmations' candidates reaching the update's chi2
  kFollowWitnessTrimmed,       ///< rows the bidding left that are dropped with the witness of a kink step
  kFollowHoleTests,            ///< layers without a passing row tested for the track crossing a sensor
  kFollowHoles,                ///< of them, holes: a sensor crossed within its bounds
  kFollowUnconfirmed,          ///< last rows dropped, no row confirming them where a sensor was crossed
  kFollowDensityKept,          ///< of the unconfirmed last rows, kept on their density (followKeepDensity)
  kFollowEndRefused,           ///< inward last rows refused at the layout's end on sqrt|V| (followEndMaxSqrtDet)
  kFollowTrimFits,             ///< tracks refitted without their last outward row (followTrimMargin)
  kFollowTrimmed,              ///< of them, the row trimmed: two or more passed its gate, the track kept the cuts
  kFollowTrimUpdates,          ///< the updates counting the rows that pass the trimmed row's gate
  kFollowEndFits,              ///< roads fitted again without their outermost row (followEnd)
  kFollowEndReplaced,          ///< of them, the follow's first row another than that row
  kFollowEarlierProposed,      ///< rows the earlier iterations' follow proposes (`followEarlier`)
  kFollowEarlierOutbid,        ///< of them, cut by the bidding
  kFollowEarlierRefitRejected, ///< of the rest, dropped with their track's failed refit
  kTrackletInWindow,           ///< candidates inside their query's own phi and second-coordinate window
  kTrackletAcceptedOutside,    ///< tracklets emitted whose target lies outside that window
  kPairBetweenBoundIn,         ///< of kPairBetweenBound, candidates inside the window
  kPairTestsIn,                ///< of kPairTests, inside the window
  kPairLongitudinalTestsIn,    ///< candidates reaching the longitudinal test (kPairTransverse) that are inside the window
  kFollowBeamEndRefused,       ///< follow: inward last rows refused against the beam line (followBeamEnd)
  kTrackletShape,              ///< candidates whose line crosses a measurement outside its slope interval (shapeTest)
  kTripletShape,               ///< thirds whose triplet's line crosses a measurement outside its slope interval
  kCount
};

inline constexpr std::size_t kWorkCounterCount = static_cast<std::size_t>(WorkCounter::kCount);
using WorkCounts = std::array<std::uint64_t, kWorkCounterCount>;

/// The tracklet stage's counts per link, indexed by the layout's link.
enum class LinkCounter : std::uint8_t {
  kQueries,    ///< grid queries, as kTrackletQueries
  kCandidates, ///< candidates reaching the pair cut: not held by an earlier track, not crossing another layer
  kTracklets,  ///< tracklets emitted
  kCount
};
inline constexpr std::size_t kLinkCounterCount = static_cast<std::size_t>(LinkCounter::kCount);
using LinkWorkCounts = std::vector<std::array<std::uint64_t, kLinkCounterCount>>;

#if defined(GTRK_WORK_COUNTERS) && !defined(__CUDACC__) && !defined(__HIPCC__)

/// Whether this build counts.
inline constexpr bool kWorkCounters = true;

namespace detail
{
/// The calling thread's block, registered on first use.
[[nodiscard]] WorkCounts& threadWorkCounts();
/// The calling thread's per-link counts, in the same block.
[[nodiscard]] LinkWorkCounts& threadLinkWorkCounts();
} // namespace detail

/// Every thread's counts since the process began, summed.
[[nodiscard]] WorkCounts workCounts();
/// The same per link.
[[nodiscard]] LinkWorkCounts linkWorkCounts();

inline void countWork(WorkCounter counter, std::uint64_t amount = 1U) { detail::threadWorkCounts()[static_cast<std::size_t>(counter)] += amount; }

inline void countLinkWork(std::size_t link, LinkCounter counter, std::uint64_t amount = 1U)
{
  auto& counts = detail::threadLinkWorkCounts();
  if (link >= counts.size()) {
    counts.resize(link + 1U);
  }
  counts[link][static_cast<std::size_t>(counter)] += amount;
}

namespace detail
{
/// Whether the calling thread's grid queries are the follow's.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
inline thread_local bool tFollowGrid = false;
} // namespace detail

/// A grid counter, kGridQueries to kGridVisits: the follow's own while a
/// FollowGridScope is open on the calling thread.
inline void countGridWork(WorkCounter counter, std::uint64_t amount = 1U)
{
  if (detail::tFollowGrid) {
    counter = static_cast<WorkCounter>(static_cast<std::size_t>(WorkCounter::kFollowGridQueries) + static_cast<std::size_t>(counter) -
                                       static_cast<std::size_t>(WorkCounter::kGridQueries));
  }
  countWork(counter, amount);
}

/// The grid's counts on the calling thread go to the follow's counters
/// until destroyed.
class FollowGridScope
{
 public:
  FollowGridScope() : mOuter{detail::tFollowGrid} { detail::tFollowGrid = true; }
  FollowGridScope(const FollowGridScope&) = delete;
  FollowGridScope& operator=(const FollowGridScope&) = delete;
  FollowGridScope(FollowGridScope&&) = delete;
  FollowGridScope& operator=(FollowGridScope&&) = delete;
  ~FollowGridScope() { detail::tFollowGrid = mOuter; }

 private:
  bool mOuter;
};

/// Adds the nanoseconds of its lifetime to a counter.
class WorkTimer
{
 public:
  explicit WorkTimer(WorkCounter counter) : mCounter{counter}, mStart{std::chrono::steady_clock::now()} {}
  WorkTimer(const WorkTimer&) = delete;
  WorkTimer& operator=(const WorkTimer&) = delete;
  WorkTimer(WorkTimer&&) = delete;
  WorkTimer& operator=(WorkTimer&&) = delete;
  ~WorkTimer()
  {
    const auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - mStart).count();
    countWork(mCounter, static_cast<std::uint64_t>(nanoseconds));
  }

 private:
  WorkCounter mCounter;
  std::chrono::steady_clock::time_point mStart;
};

#else

inline constexpr bool kWorkCounters = false;

[[nodiscard]] inline WorkCounts workCounts() { return {}; }
[[nodiscard]] inline LinkWorkCounts linkWorkCounts() { return {}; }

GPUhdi() GPUAlwaysInline() void countWork(WorkCounter /*counter*/, std::uint64_t /*amount*/ = 1U) {}
GPUhdi() GPUAlwaysInline() void countLinkWork(std::size_t /*link*/, LinkCounter /*counter*/, std::uint64_t /*amount*/ = 1U) {}
GPUhdi() GPUAlwaysInline() void countGridWork(WorkCounter /*counter*/, std::uint64_t /*amount*/ = 1U) {}

class FollowGridScope
{
 public:
  // User-provided, so that a const scope -- the form its users declare -- is
  // well-formed without members to initialise.
  GPUhdi() GPUAlwaysInline() FollowGridScope() {} // NOLINT(modernize-use-equals-default)
};

class WorkTimer
{
 public:
  GPUhdi() GPUAlwaysInline() explicit WorkTimer(WorkCounter /*counter*/) {}
};

#endif

} // namespace o2::gtrk

#endif // GTRK_TRACKING_WORKCOUNTERS_H_
