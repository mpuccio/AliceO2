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

/// \file VertexGlue.cuh
/// \brief The vertexer on the device: the
///        vertexing pass's triplets to the timeframe's vertex list, as
///        `Vertexer::find` makes it on the host.
///
///   Candidates: one candidate per triplet: `vertexCandidate`, its r-z chi2,
///       its box; the accepted ones made cluster-disjoint in rounds of bids -- each
///       live one bids (chi2, place) for its rows, the lowest everywhere is
///       kept, those on its rows dropped -- the host's greedy pass in chi2
///       order
///   Order: the kept ones in triplet order, sorted stably by (frame, z)
///   Rounds: under `sweepToFixedPoint`, every frame at once: scans of the
///       free candidates' counts and z; each candidate's window and
///       sidebands by binary search in its frame; the frame's peak in two
///       stages -- the largest excess, then the lowest index with it -- over
///       the whole 64-bit excess; the centre from the scans; the claims, their
///       z's sums and the box by integer atomics, one set per frame and round;
///       then the round's vertices own the measurements of their frames'
///       triplets within the claim -- the cut ones too -- and the free
///       candidates with an owned one are dropped, each
///       measurement a slot per frame, (frame << 32) | row, sorted and unique
///   Seeds: the sums read back; the host makes the seeds (`vertexSeedOf`) and
///       orders each frame's vertices by z, ties by round, as the host does
///
/// Everything after the map is integer, or float on integers, so the list
/// does not depend on the schedule; the map itself is float, and a device's
/// candidates can differ from the host's where its contractions round
/// otherwise. Included by DeviceTracker.cuh, in the backend's namespace.

#ifndef GTRK_TRACKING_GPU_VERTEXGLUE_CUH_
#define GTRK_TRACKING_GPU_VERTEXGLUE_CUH_

#include "DeviceAlgorithms.cuh"
#include "Platform.cuh"

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/Triplet.h"
#include "GenericTracking/Tracking/VertexStage.h"

#include "GPUCommonAlgorithm.h"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::vertex_glue
{

using Index = std::uint32_t;

namespace detail
{
namespace
{

constexpr Index kNoIndex = std::numeric_limits<Index>::max();
constexpr std::uint64_t kNoBid = std::numeric_limits<std::uint64_t>::max();

/// Candidates: each triplet's candidate key (0 for none), its key as an ownership
/// source (0 for a degenerate one), its r-z chi2, its weight
/// from its z0 error -- the pass's material and the beam's transverse spread,
/// `ScatteringErrors` -- and its box.
template <DetectorLayout Det>
struct CandidateMap {
  MeasurementView measurements;
  TrackletConstView tracklets;
  TripletConstView triplets;
  VertexerSettings settings;
  const IterationParameters<Det>* pass = nullptr;
  float varTransverse = 0.F;
  std::uint64_t* keys = nullptr;
  std::uint64_t* sourceKeys = nullptr;
  float* chi2 = nullptr;
  std::uint64_t* weights = nullptr;
  Index* lower = nullptr;
  Index* upper = nullptr;
  GPUdi() void operator()(Index triplet) const
  {
    const ScatteringErrors<Det> errors{.table = &pass->scattering, .cut = pass->pairCut, .varTransverse = varTransverse};
    VertexCandidateFit fit;
    TimeEstPS box;
    const auto status = measureCandidate(measurements, tracklets, triplets, triplet, settings, errors, fit, weights[triplet], box, nullptr);
    keys[triplet] = status == VertexCandidateStatus::kAccepted ? fit.key : 0U;
    sourceKeys[triplet] = status != VertexCandidateStatus::kDegenerate ? fit.key : 0U;
    chi2[triplet] = fit.chi2RZ;
    lower[triplet] = box.lower();
    upper[triplet] = box.upper();
  }
};

struct Accepted {
  const std::uint64_t* keys = nullptr;
  GPUdi() bool operator()(Index triplet) const { return keys[triplet] != 0U; }
};

/// The deduplication's state of an accepted candidate.
enum : std::uint8_t { kLive = 0U,
                      kKept = 1U,
                      kDropped = 2U };

/// What the deduplication's rounds read and write: the accepted candidates
/// in triplet order, their rows and bids, and per row its lowest bid and
/// whether a kept candidate holds it.
struct Deduplication {
  const Index* accepted = nullptr;
  TrackletConstView tracklets;
  TripletConstView triplets;
  const float* chi2 = nullptr;
  std::uint8_t* state = nullptr;
  std::uint64_t* slots = nullptr;
  std::uint8_t* taken = nullptr;
  GPUdi() std::array<Index, 3> rows(Index candidate) const { return tripletRows(tracklets, triplets, accepted[candidate]); }
  /// The bid: the chi2's bits -- non-negative, so they order as the values
  /// -- then the candidate's place, the lower triplet first on ties.
  GPUdi() std::uint64_t bid(Index candidate) const { return (static_cast<std::uint64_t>(std::bit_cast<std::uint32_t>(chi2[accepted[candidate]])) << 32U) | candidate; }
};

/// Every live candidate's rows open for this round's bids.
struct Reopen {
  Deduplication d;
  GPUdi() void operator()(Index candidate) const
  {
    if (d.state[candidate] != kLive) {
      return;
    }
    const auto rows = d.rows(candidate);
    for (const Index row : rows) {
      d.slots[row] = kNoBid;
    }
  }
};

/// Every live candidate bids for its rows.
struct Bid {
  Deduplication d;
  GPUdi() void operator()(Index candidate) const
  {
    if (d.state[candidate] != kLive) {
      return;
    }
    const auto rows = d.rows(candidate);
    const auto bid = static_cast<unsigned long long>(d.bid(candidate));
    for (const Index row : rows) {
      atomicMin(reinterpret_cast<unsigned long long*>(d.slots) + row, bid);
    }
  }
};

/// A live candidate holding the lowest bid on all its rows is kept: nothing
/// before it in the order is left on any of them. The round's change flag
/// goes up.
struct Decide {
  Deduplication d;
  Index* changed = nullptr;
  GPUdi() void operator()(Index candidate) const
  {
    if (d.state[candidate] != kLive) {
      return;
    }
    const auto rows = d.rows(candidate);
    const std::uint64_t bid = d.bid(candidate);
    for (const Index row : rows) {
      if (d.slots[row] != bid) {
        return;
      }
    }
    d.state[candidate] = kKept;
    for (const Index row : rows) {
      d.taken[row] = 1U;
    }
    *changed = 1U;
  }
};

/// A live candidate on a row a kept one holds is dropped.
struct Drop {
  Deduplication d;
  GPUdi() void operator()(Index candidate) const
  {
    if (d.state[candidate] != kLive) {
      return;
    }
    const auto rows = d.rows(candidate);
    for (const Index row : rows) {
      if (d.taken[row] != 0U) {
        d.state[candidate] = kDropped;
        return;
      }
    }
  }
};

/// The kept candidates of the owned frames (`OwnedFrames`).
struct KeptOnes {
  const std::uint8_t* state = nullptr;
  const Index* accepted = nullptr;
  const std::uint64_t* keys = nullptr;
  OwnedFrames owned;
  GPUdi() bool operator()(Index candidate) const { return state[candidate] == kKept && owned.owns(vertexKeyFrame(keys[accepted[candidate]])); }
};

/// The kept candidates' keys and triplets, in triplet order, for the sort.
struct GatherKept {
  const Index* kept = nullptr;
  const Index* accepted = nullptr;
  const std::uint64_t* keys = nullptr;
  std::uint64_t* candidateKeys = nullptr;
  Index* candidateTriplets = nullptr;
  GPUdi() void operator()(Index candidate) const
  {
    const Index triplet = accepted[kept[candidate]];
    candidateKeys[candidate] = keys[triplet];
    candidateTriplets[candidate] = triplet;
  }
};

/// The sorted candidates' z [um] and box, and whether each begins a frame.
struct SortedColumns {
  const std::uint64_t* sortedKeys = nullptr;
  const Index* sortedTriplets = nullptr;
  const Index* lower = nullptr;
  const Index* upper = nullptr;
  const std::uint64_t* weights = nullptr;
  long long* z = nullptr;
  std::uint64_t* weight = nullptr;
  Index* boxLower = nullptr;
  Index* boxUpper = nullptr;
  Index* frameStart = nullptr;
  GPUdi() void operator()(Index candidate) const
  {
    const std::uint64_t key = sortedKeys[candidate];
    z[candidate] = vertexKeyZUM(key);
    const Index triplet = sortedTriplets[candidate];
    weight[candidate] = weights[triplet];
    boxLower[candidate] = lower[triplet];
    boxUpper[candidate] = upper[triplet];
    frameStart[candidate] = candidate == 0U || vertexKeyFrame(key) != vertexKeyFrame(sortedKeys[candidate - 1U]) ? 1U : 0U;
  }
};

struct FrameStarts {
  const Index* frameStart = nullptr;
  GPUdi() bool operator()(Index candidate) const { return frameStart[candidate] != 0U; }
};

/// Each candidate's frame, from the starts' exclusive scan.
struct FrameOf {
  const Index* frameStart = nullptr;
  const Index* startsBefore = nullptr;
  Index* frame = nullptr;
  GPUdi() void operator()(Index candidate) const { frame[candidate] = startsBefore[candidate] + frameStart[candidate] - 1U; }
};

/// The ownership's sources: every triplet of the owned
/// frames but the degenerate.
struct Sources {
  const std::uint64_t* sourceKeys = nullptr;
  OwnedFrames owned;
  GPUdi() bool operator()(Index triplet) const { return sourceKeys[triplet] != 0U && owned.owns(vertexKeyFrame(sourceKeys[triplet])); }
};

/// Each source's three slot keys, (frame << 32) | row (`ownedSlotKey`).
struct SlotKeys {
  const Index* sources = nullptr;
  const std::uint64_t* sourceKeys = nullptr;
  TrackletConstView tracklets;
  TripletConstView triplets;
  std::uint64_t* slotKeys = nullptr;
  GPUdi() void operator()(Index source) const
  {
    const Index triplet = sources[source];
    const auto rows = tripletRows(tracklets, triplets, triplet);
    for (Index k = 0; k < 3U; ++k) {
      slotKeys[(3U * source) + k] = ownedSlotKey(sourceKeys[triplet], rows[k]);
    }
  }
};

/// The first of each run of equal sorted slot keys.
struct SlotStarts {
  const std::uint64_t* sorted = nullptr;
  GPUdi() bool operator()(Index i) const { return i == 0U || sorted[i] != sorted[i - 1U]; }
};

struct GatherSlots {
  const Index* starts = nullptr;
  const std::uint64_t* sorted = nullptr;
  std::uint64_t* unique = nullptr;
  GPUdi() void operator()(Index slot) const { unique[slot] = sorted[starts[slot]]; }
};

/// The slots of a triplet's measurements: their keys' places among the
/// unique ones.
struct SlotsOf {
  const std::uint64_t* unique = nullptr;
  Index slots = 0U;
  TrackletConstView tracklets;
  TripletConstView triplets;
  GPUdi() void operator()(std::uint64_t key, Index triplet, Index* out) const
  {
    const auto rows = tripletRows(tracklets, triplets, triplet);
    for (Index k = 0; k < 3U; ++k) {
      const std::uint64_t slotKey = ownedSlotKey(key, rows[k]);
      out[k] = o2::gpu::CAAlgo::partitionPoint(Index{0}, slots, [&](Index i) { return unique[i] < slotKey; });
    }
  }
};

/// Each source's slots, z [um] and frame among the candidates' frames
/// (kNoIndex if its frame has no candidate).
struct SourceColumns {
  const Index* sources = nullptr;
  const std::uint64_t* sourceKeys = nullptr;
  SlotsOf slotsOf;
  const std::uint64_t* sortedKeys = nullptr;
  const Index* frameBegin = nullptr;
  Index frames = 0U;
  Index* sourceSlots = nullptr;
  long long* sourceZ = nullptr;
  Index* sourceFrame = nullptr;
  GPUdi() void operator()(Index source) const
  {
    const Index triplet = sources[source];
    const std::uint64_t key = sourceKeys[triplet];
    slotsOf(key, triplet, sourceSlots + (3U * source));
    sourceZ[source] = vertexKeyZUM(key);
    const std::int64_t frame = vertexKeyFrame(key);
    const Index f = o2::gpu::CAAlgo::partitionPoint(Index{0}, frames, [&](Index g) { return vertexKeyFrame(sortedKeys[frameBegin[g]]) < frame; });
    sourceFrame[source] = f < frames && vertexKeyFrame(sortedKeys[frameBegin[f]]) == frame ? f : kNoIndex;
  }
};

/// Each sorted candidate's slots.
struct CandidateSlots {
  const std::uint64_t* sortedKeys = nullptr;
  const Index* sortedTriplets = nullptr;
  SlotsOf slotsOf;
  Index* candidateSlots = nullptr;
  GPUdi() void operator()(Index candidate) const { slotsOf(sortedKeys[candidate], sortedTriplets[candidate], candidateSlots + (3U * candidate)); }
};

/// What the rounds read of the candidates and their frames.
struct Candidates {
  const long long* z = nullptr;
  const std::uint64_t* weight = nullptr;
  const Index* boxLower = nullptr;
  const Index* boxUpper = nullptr;
  const Index* frame = nullptr;
  const Index* frameBegin = nullptr;
  Index frames = 0U;
  Index count = 0U;
  GPUdi() Index frameEnd(Index f) const { return f + 1U < frames ? frameBegin[f + 1U] : count; }
  /// The first candidate of [\a begin, \a end) whose z is at least \a value:
  /// the host's `std::ranges::lower_bound` in the frame.
  GPUdi() Index at(Index begin, Index end, long long value) const
  {
    return o2::gpu::CAAlgo::partitionPoint(begin, end - begin, [&](Index i) { return z[i] < value; });
  }
};

/// One round's per-frame state: the peak and the vertex it makes.
struct Round {
  Counter* bestExcess = nullptr;
  Index* bestIndex = nullptr;
  float* centre = nullptr;
  Index* contributors = nullptr;
  Index* earliest = nullptr;
  Index* latest = nullptr;
};

/// The free candidates' count and z, for the round's scans; the entry past
/// the last is zero, so the scans' last entry is the total.
struct FreeColumns {
  Candidates candidates;
  const std::int32_t* claims = nullptr;
  Index* freeCount = nullptr;
  long long* freeZ = nullptr;
  GPUdi() void operator()(Index candidate) const
  {
    const bool free = candidate < candidates.count && claims[candidate] == kClaimFree;
    freeCount[candidate] = free ? 1U : 0U;
    freeZ[candidate] = free ? candidates.z[candidate] : 0;
  }
};

/// Each free candidate's window, its sidebands and its excess 2W - S if it
/// may be a peak (0 otherwise); the frame's largest excess by atomic max.
struct Score {
  Candidates candidates;
  const std::int32_t* claims = nullptr;
  const Index* countBefore = nullptr;
  long long window = 0;
  long long claim = 0;
  Index minContributors = 0U;
  float k2 = 0.F;
  long long* excess = nullptr;
  Round round;
  GPUdi() void operator()(Index i) const
  {
    excess[i] = 0;
    if (claims[i] != kClaimFree) {
      return;
    }
    const Index f = candidates.frame[i];
    const Index begin = candidates.frameBegin[f];
    const Index end = candidates.frameEnd(f);
    const long long z = candidates.z[i];
    const auto countIn = [&](std::int64_t low, std::int64_t high) {
      return countBefore[candidates.at(begin, end, high)] - countBefore[candidates.at(begin, end, low)];
    };
    const std::int64_t score = peakExcess(z, window, claim, minContributors, k2, countIn);
    if (score <= 0) {
      return;
    }
    excess[i] = score;
    atomicMax(round.bestExcess + f, static_cast<Counter>(score));
  }
};

/// The lowest index among the frame's candidates with its largest excess.
struct Pick {
  Candidates candidates;
  const long long* excess = nullptr;
  Round round;
  GPUdi() void operator()(Index i) const
  {
    const Index f = candidates.frame[i];
    if (excess[i] > 0 && static_cast<Counter>(excess[i]) == round.bestExcess[f]) {
      atomicMin(round.bestIndex + f, i);
    }
  }
};

/// Each frame with a peak: the centre, the mean of the window's free
/// candidates' z, and the round's change flag raised.
struct Centre {
  Candidates candidates;
  const Index* countBefore = nullptr;
  const long long* zBefore = nullptr;
  long long window = 0;
  Round round;
  Index* changed = nullptr;
  GPUdi() void operator()(Index f) const
  {
    const Index best = round.bestIndex[f];
    if (best == kNoIndex) {
      return;
    }
    const Index windowEnd = candidates.at(candidates.frameBegin[f], candidates.frameEnd(f), candidates.z[best] + window);
    const long long sum = zBefore[windowEnd] - zBefore[best];
    const Index members = countBefore[windowEnd] - countBefore[best];
    round.centre[f] = static_cast<float>(sum) / static_cast<float>(members);
    *changed = 1U;
  }
};

/// Each free candidate of a frame with a peak, within the claim of its
/// centre: claimed for the round's vertex, counted, its box taken.
struct Claim {
  Candidates candidates;
  std::int32_t* claims = nullptr;
  float claim = 0.F;
  std::int32_t vertexBase = 0; ///< the round's first vertex number: round * frames
  Round round;
  GPUdi() void operator()(Index i) const
  {
    if (claims[i] != kClaimFree) {
      return;
    }
    const Index f = candidates.frame[i];
    if (round.bestIndex[f] == kNoIndex) {
      return;
    }
    if (!withinZ(candidates.z[i], round.centre[f], claim)) {
      return;
    }
    claims[i] = vertexBase + static_cast<std::int32_t>(f);
    atomicAdd(round.contributors + f, 1U);
    atomicMin(round.earliest + f, candidates.boxLower[i]);
    atomicMax(round.latest + f, candidates.boxUpper[i]);
  }
};

/// Each source of a frame with a vertex this round, within its claim: its
/// measurements owned.
struct Own {
  const Index* sourceSlots = nullptr;
  const long long* sourceZ = nullptr;
  const Index* sourceFrame = nullptr;
  std::uint8_t* ownedSlots = nullptr;
  float claim = 0.F;
  Round round;
  GPUdi() void operator()(Index source) const
  {
    const Index f = sourceFrame[source];
    if (f == kNoIndex || round.bestIndex[f] == kNoIndex || !withinZ(sourceZ[source], round.centre[f], claim)) {
      return;
    }
    for (Index k = 0; k < 3U; ++k) {
      ownedSlots[sourceSlots[(3U * source) + k]] = 1U;
    }
  }
};

/// Each free candidate with an owned measurement is dropped.
struct DropOwned {
  const Index* candidateSlots = nullptr;
  const std::uint8_t* ownedSlots = nullptr;
  std::int32_t* claims = nullptr;
  GPUdi() void operator()(Index i) const
  {
    if (claims[i] != kClaimFree) {
      return;
    }
    for (Index k = 0; k < 3U; ++k) {
      if (ownedSlots[candidateSlots[(3U * i) + k]] != 0U) {
        claims[i] = kClaimDropped;
        return;
      }
    }
  }
};

/// What a vertex's z and error come from (`VertexClaim`): its selected
/// candidates, grown from the core, one slot
/// per vertex number -- round * frames + frame, as the claims hold them.
struct Selection {
  Candidates candidates;
  const std::int32_t* claims = nullptr;
  std::uint8_t* selected = nullptr;
  Counter* weightSum = nullptr;
  Counter* weightedZSum = nullptr;
  Counter* scatterSum = nullptr;
  Counter* count = nullptr;
  float* mean = nullptr;
  float* varMean = nullptr;
  GPUdi() void select(Index i, Index v) const
  {
    selected[i] = 1U;
    const std::uint64_t w = candidates.weight[i];
    atomicAdd(weightSum + v, static_cast<Counter>(w));
    atomicAdd(weightedZSum + v, static_cast<Counter>(static_cast<long long>(w) * candidates.z[i]));
    atomicAdd(count + v, Counter{1U});
  }
};

/// The core: each claimed candidate within the window of its vertex's
/// centre starts the selection.
struct CoreStart {
  Selection selection;
  const float* centre = nullptr;
  float window = 0.F;
  GPUdi() void operator()(Index i) const
  {
    const std::int32_t v = selection.claims[i];
    if (v >= 0 && withinZ(selection.candidates.z[i], centre[v], window)) {
      selection.select(i, static_cast<Index>(v));
    }
  }
};

/// Each vertex's weighted mean and its variance from its selection.
struct Moments {
  Selection selection;
  GPUdi() void operator()(Index v) const
  {
    if (selection.count[v] == 0U) {
      return;
    }
    const auto moments = zMoments(static_cast<std::uint64_t>(selection.weightSum[v]), static_cast<std::int64_t>(selection.weightedZSum[v]));
    selection.mean[v] = moments.mean;
    selection.varMean[v] = moments.varMean;
  }
};

/// Each claimed candidate not yet selected joins if it lies within nSigma of
/// its own and the mean's errors of its vertex's mean; the change flag goes up.
struct Grow {
  Selection selection;
  float nSigma = 0.F;
  Index* changed = nullptr;
  GPUdi() void operator()(Index i) const
  {
    const std::int32_t v = selection.claims[i];
    if (v < 0 || selection.selected[i] != 0U) {
      return;
    }
    const auto at = static_cast<Index>(v);
    if (withinError(selection.candidates.z[i], selection.candidates.weight[i], ZMoments{.mean = selection.mean[at], .varMean = selection.varMean[at]},
                    nSigma)) {
      selection.select(i, at);
      *changed = 1U;
    }
  }
};

/// Each selected candidate's share of its vertex's scatter about the settled
/// weighted mean (`candidateScatter`), for the sandwich error.
struct Scatter {
  Selection selection;
  GPUdi() void operator()(Index i) const
  {
    const std::int32_t v = selection.claims[i];
    if (v < 0 || selection.selected[i] == 0U) {
      return;
    }
    atomicAdd(selection.scatterSum + v, static_cast<Counter>(candidateScatter(selection.candidates.z[i], selection.candidates.weight[i], selection.mean[v])));
  }
};

/// A round's per-frame arrays on the base end.
struct RoundStorage {
  DeviceSpan<Counter> bestExcess;
  DeviceSpan<Index> bestIndex;
  DeviceSpan<float> centre;
  DeviceSpan<Index> contributors;
  DeviceSpan<Index> earliest;
  DeviceSpan<Index> latest;

  [[nodiscard]] bool allocate(DeviceArena& arena, Index frames)
  {
    return bestExcess.allocate(arena, frames) && bestIndex.allocate(arena, frames) && centre.allocate(arena, frames) && contributors.allocate(arena, frames) &&
           earliest.allocate(arena, frames) && latest.allocate(arena, frames);
  }
  void reset(DeviceRuntime& runtime, DeviceStream& stream, Index frames)
  {
    runtime.fill(bestExcess.data(), 0U, frames * sizeof(Counter), stream);
    runtime.fill(bestIndex.data(), 0xFFU, frames * sizeof(Index), stream);
    runtime.fill(contributors.data(), 0U, frames * sizeof(Index), stream);
    runtime.fill(earliest.data(), 0xFFU, frames * sizeof(Index), stream);
    runtime.fill(latest.data(), 0U, frames * sizeof(Index), stream);
  }
  [[nodiscard]] Round view()
  {
    return {.bestExcess = bestExcess.data(),
            .bestIndex = bestIndex.data(),
            .centre = centre.data(),
            .contributors = contributors.data(),
            .earliest = earliest.data(),
            .latest = latest.data()};
  }
};

} // namespace
} // namespace detail

/// The vertices of the pass's \a tripletCount triplets on \a measurements,
/// \a rows of them, into \a result (`DeviceTracker::findVertices`); \a pass
/// the pass's parameters on the device, for the candidates' errors; the
/// vertices of \a owned's frames, a slice's, as `Vertexer::find` finds them.
/// Its arrays on the base end, in the caller's scope; kOutOfMemory if the
/// arena refuses.
template <DetectorLayout Det>
[[nodiscard]] StageStatus findVertices(DeviceRuntime& runtime, DeviceArena& arena, DeviceStream& stream, MeasurementView measurements, Index rows,
                                       TrackletConstView tracklets, TripletConstView triplets, Index tripletCount, const IterationParameters<Det>* pass,
                                       const VertexerSettings& settings, const VertexConstraint& beamSpot, const OwnedFrames& owned, VertexerResult& result)
{
  using namespace detail;
  result = {};
  result.triplets = tripletCount;
  if (tripletCount == 0U) {
    return StageStatus::kComplete;
  }
  // Candidates.
  DeviceSpan<std::uint64_t> keys;
  DeviceSpan<std::uint64_t> sourceKeys;
  DeviceSpan<float> chi2;
  DeviceSpan<std::uint64_t> weights;
  DeviceSpan<Index> lower;
  DeviceSpan<Index> upper;
  if (!keys.allocate(arena, tripletCount) || !sourceKeys.allocate(arena, tripletCount) || !chi2.allocate(arena, tripletCount) ||
      !weights.allocate(arena, tripletCount) || !lower.allocate(arena, tripletCount) || !upper.allocate(arena, tripletCount)) {
    return StageStatus::kOutOfMemory;
  }
  algorithms::forEach(runtime, stream, tripletCount,
                      CandidateMap<Det>{.measurements = measurements,
                                        .tracklets = tracklets,
                                        .triplets = triplets,
                                        .settings = settings,
                                        .pass = pass,
                                        .varTransverse = static_cast<float>(beamSpot.sigmaTransverse * beamSpot.sigmaTransverse),
                                        .keys = keys.data(),
                                        .sourceKeys = sourceKeys.data(),
                                        .chi2 = chi2.data(),
                                        .weights = weights.data(),
                                        .lower = lower.data(),
                                        .upper = upper.data()});
  DeviceSpan<Index> accepted;
  Index acceptedCount = 0U;
  StageStatus status = algorithms::selectIndices(runtime, arena, stream, tripletCount, Accepted{keys.data()}, accepted, acceptedCount);
  if (status != StageStatus::kComplete) {
    return status;
  }
  result.accepted = acceptedCount;
  if (acceptedCount == 0U) {
    return StageStatus::kComplete;
  }

  // The deduplication: rounds of bids until
  // every accepted candidate is kept or dropped -- the host's greedy pass in
  // chi2 order, since a candidate is kept only when nothing before it is left
  // on any of its rows.
  DeviceSpan<std::uint8_t> state;
  DeviceSpan<std::uint64_t> slots;
  DeviceSpan<std::uint8_t> taken;
  if (!state.allocate(arena, acceptedCount) || !slots.allocate(arena, rows) || !taken.allocate(arena, rows)) {
    return StageStatus::kOutOfMemory;
  }
  runtime.fill(state.data(), kLive, acceptedCount, stream);
  runtime.fill(taken.data(), 0U, rows, stream);
  const Deduplication deduplication{.accepted = accepted.data(),
                                    .tracklets = tracklets,
                                    .triplets = triplets,
                                    .chi2 = chi2.data(),
                                    .state = state.data(),
                                    .slots = slots.data(),
                                    .taken = taken.data()};
  const auto deduplicated = algorithms::sweepToFixedPoint<1U>(runtime, arena, stream, acceptedCount + 1U, [&](Index* changed) {
    algorithms::forEach(runtime, stream, acceptedCount, Reopen{deduplication});
    algorithms::forEach(runtime, stream, acceptedCount, Bid{deduplication});
    algorithms::forEach(runtime, stream, acceptedCount, Decide{deduplication, changed});
    algorithms::forEach(runtime, stream, acceptedCount, Drop{deduplication});
  });
  if (deduplicated.status != StageStatus::kComplete) {
    return deduplicated.status;
  }
  DeviceSpan<Index> kept;
  Index count = 0U;
  status = algorithms::selectIndices(runtime, arena, stream, acceptedCount, KeptOnes{state.data(), accepted.data(), keys.data(), owned}, kept, count);
  if (status != StageStatus::kComplete) {
    return status;
  }
  result.candidates = count;
  if (count == 0U) {
    return StageStatus::kComplete;
  }

  // Order: the kept ones in triplet order, then stably by key.
  DeviceSpan<std::uint64_t> candidateKeys;
  DeviceSpan<Index> candidateTriplets;
  DeviceSpan<std::uint64_t> sortedKeys;
  DeviceSpan<Index> sortedTriplets;
  if (!candidateKeys.allocate(arena, count) || !candidateTriplets.allocate(arena, count) || !sortedKeys.allocate(arena, count) ||
      !sortedTriplets.allocate(arena, count)) {
    return StageStatus::kOutOfMemory;
  }
  algorithms::forEach(runtime, stream, count, GatherKept{kept.data(), accepted.data(), keys.data(), candidateKeys.data(), candidateTriplets.data()});
  status = algorithms::sortPairs(runtime, arena, stream, std::as_const(candidateKeys).data(), sortedKeys.data(), std::as_const(candidateTriplets).data(),
                                 sortedTriplets.data(), count, 64);
  if (status != StageStatus::kComplete) {
    return status;
  }
  DeviceSpan<long long> z;
  DeviceSpan<Index> boxLower;
  DeviceSpan<Index> boxUpper;
  DeviceSpan<Index> frameStart;
  DeviceSpan<Index> startsBefore;
  DeviceSpan<Index> frameOf;
  DeviceSpan<std::uint64_t> weight;
  if (!weight.allocate(arena, count) || !z.allocate(arena, count) || !boxLower.allocate(arena, count) || !boxUpper.allocate(arena, count) ||
      !frameStart.allocate(arena, count) || !startsBefore.allocate(arena, count) || !frameOf.allocate(arena, count)) {
    return StageStatus::kOutOfMemory;
  }
  algorithms::forEach(runtime, stream, count,
                      SortedColumns{.sortedKeys = sortedKeys.data(),
                                    .sortedTriplets = sortedTriplets.data(),
                                    .lower = lower.data(),
                                    .upper = upper.data(),
                                    .weights = weights.data(),
                                    .z = z.data(),
                                    .weight = weight.data(),
                                    .boxLower = boxLower.data(),
                                    .boxUpper = boxUpper.data(),
                                    .frameStart = frameStart.data()});
  DeviceSpan<Index> frameBegin;
  Index frames = 0U;
  status = algorithms::selectIndices(runtime, arena, stream, count, FrameStarts{frameStart.data()}, frameBegin, frames);
  if (status != StageStatus::kComplete) {
    return status;
  }
  status = algorithms::exclusiveScan(runtime, arena, stream, std::as_const(frameStart).data(), startsBefore.data(), count);
  if (status != StageStatus::kComplete) {
    return status;
  }
  algorithms::forEach(runtime, stream, count, FrameOf{frameStart.data(), startsBefore.data(), frameOf.data()});
  const Candidates candidates{.z = z.data(),
                              .weight = weight.data(),
                              .boxLower = boxLower.data(),
                              .boxUpper = boxUpper.data(),
                              .frame = frameOf.data(),
                              .frameBegin = frameBegin.data(),
                              .frames = frames,
                              .count = count};

  // The ownership's sources and slots.
  DeviceSpan<Index> sources;
  Index sourceCount = 0U;
  status = algorithms::selectIndices(runtime, arena, stream, tripletCount, Sources{std::as_const(sourceKeys).data(), owned}, sources, sourceCount);
  if (status != StageStatus::kComplete) {
    return status;
  }
  const Index slotEntries = 3U * sourceCount;
  DeviceSpan<std::uint64_t> slotKeys;
  DeviceSpan<std::uint64_t> sortedSlotKeys;
  if (!slotKeys.allocate(arena, slotEntries) || !sortedSlotKeys.allocate(arena, slotEntries)) {
    return StageStatus::kOutOfMemory;
  }
  algorithms::forEach(runtime, stream, sourceCount,
                      SlotKeys{.sources = std::as_const(sources).data(),
                               .sourceKeys = std::as_const(sourceKeys).data(),
                               .tracklets = tracklets,
                               .triplets = triplets,
                               .slotKeys = slotKeys.data()});
  status = algorithms::sortKeys(runtime, arena, stream, std::as_const(slotKeys).data(), sortedSlotKeys.data(), slotEntries, 64);
  if (status != StageStatus::kComplete) {
    return status;
  }
  DeviceSpan<Index> slotStarts;
  Index slotCount = 0U;
  status = algorithms::selectIndices(runtime, arena, stream, slotEntries, SlotStarts{std::as_const(sortedSlotKeys).data()}, slotStarts, slotCount);
  if (status != StageStatus::kComplete) {
    return status;
  }
  DeviceSpan<std::uint64_t> uniqueSlots;
  DeviceSpan<std::uint8_t> ownedSlots;
  DeviceSpan<Index> sourceSlots;
  DeviceSpan<long long> sourceZ;
  DeviceSpan<Index> sourceFrame;
  DeviceSpan<Index> candidateSlots;
  if (!uniqueSlots.allocate(arena, slotCount) || !ownedSlots.allocate(arena, slotCount) || !sourceSlots.allocate(arena, slotEntries) ||
      !sourceZ.allocate(arena, sourceCount) || !sourceFrame.allocate(arena, sourceCount) || !candidateSlots.allocate(arena, 3U * count)) {
    return StageStatus::kOutOfMemory;
  }
  algorithms::forEach(runtime, stream, slotCount, GatherSlots{std::as_const(slotStarts).data(), std::as_const(sortedSlotKeys).data(), uniqueSlots.data()});
  runtime.fill(ownedSlots.data(), 0U, slotCount, stream);
  const SlotsOf slotsOf{.unique = std::as_const(uniqueSlots).data(), .slots = slotCount, .tracklets = tracklets, .triplets = triplets};
  algorithms::forEach(runtime, stream, sourceCount,
                      SourceColumns{.sources = std::as_const(sources).data(),
                                    .sourceKeys = std::as_const(sourceKeys).data(),
                                    .slotsOf = slotsOf,
                                    .sortedKeys = std::as_const(sortedKeys).data(),
                                    .frameBegin = std::as_const(frameBegin).data(),
                                    .frames = frames,
                                    .sourceSlots = sourceSlots.data(),
                                    .sourceZ = sourceZ.data(),
                                    .sourceFrame = sourceFrame.data()});
  algorithms::forEach(runtime, stream, count,
                      CandidateSlots{.sortedKeys = std::as_const(sortedKeys).data(),
                                     .sortedTriplets = std::as_const(sortedTriplets).data(),
                                     .slotsOf = slotsOf,
                                     .candidateSlots = candidateSlots.data()});

  // Rounds: every frame at once.
  DeviceSpan<std::int32_t> claims;
  DeviceSpan<Index> freeCount;
  DeviceSpan<long long> freeZ;
  DeviceSpan<Index> countBefore;
  DeviceSpan<long long> zBefore;
  DeviceSpan<long long> excess;
  const Index scanned = count + 1U;
  if (!claims.allocate(arena, count) || !freeCount.allocate(arena, scanned) || !freeZ.allocate(arena, scanned) || !countBefore.allocate(arena, scanned) ||
      !zBefore.allocate(arena, scanned) || !excess.allocate(arena, count)) {
    return StageStatus::kOutOfMemory;
  }
  runtime.fill(claims.data(), 0xFFU, std::size_t{count} * sizeof(std::int32_t), stream);
  const auto window = static_cast<long long>(settings.windowUM);
  const auto claim = static_cast<long long>(settings.claimUM);
  const Scalar k2 = 2.F * settings.significance * settings.significance;
  std::vector<RoundStorage> rounds;
  bool outOfMemory = false;
  const auto swept = algorithms::sweepToFixedPoint<1U>(runtime, arena, stream, count + 1U, [&](Index* changed) {
    if (outOfMemory) {
      return;
    }
    auto& round = rounds.emplace_back();
    if (!round.allocate(arena, frames)) {
      outOfMemory = true;
      return;
    }
    round.reset(runtime, stream, frames);
    const auto roundView = round.view();
    algorithms::forEach(runtime, stream, scanned, FreeColumns{candidates, std::as_const(claims).data(), freeCount.data(), freeZ.data()});
    if (algorithms::exclusiveScan(runtime, arena, stream, std::as_const(freeCount).data(), countBefore.data(), scanned) != StageStatus::kComplete ||
        algorithms::exclusiveScan(runtime, arena, stream, std::as_const(freeZ).data(), zBefore.data(), scanned) != StageStatus::kComplete) {
      outOfMemory = true;
      return;
    }
    algorithms::forEach(runtime, stream, count,
                        Score{.candidates = candidates,
                              .claims = std::as_const(claims).data(),
                              .countBefore = std::as_const(countBefore).data(),
                              .window = window,
                              .claim = claim,
                              .minContributors = settings.minContributors,
                              .k2 = k2,
                              .excess = excess.data(),
                              .round = roundView});
    algorithms::forEach(runtime, stream, count, Pick{candidates, std::as_const(excess).data(), roundView});
    algorithms::forEach(runtime, stream, frames,
                        Centre{candidates, std::as_const(countBefore).data(), std::as_const(zBefore).data(), window, roundView, changed});
    algorithms::forEach(runtime, stream, count,
                        Claim{.candidates = candidates,
                              .claims = claims.data(),
                              .claim = static_cast<float>(claim),
                              .vertexBase = static_cast<std::int32_t>((rounds.size() - 1U) * frames),
                              .round = roundView});
    algorithms::forEach(runtime, stream, sourceCount,
                        Own{.sourceSlots = std::as_const(sourceSlots).data(),
                            .sourceZ = std::as_const(sourceZ).data(),
                            .sourceFrame = std::as_const(sourceFrame).data(),
                            .ownedSlots = ownedSlots.data(),
                            .claim = static_cast<float>(claim),
                            .round = roundView});
    algorithms::forEach(runtime, stream, count, DropOwned{std::as_const(candidateSlots).data(), std::as_const(ownedSlots).data(), claims.data()});
  });
  if (swept.status != StageStatus::kComplete) {
    return swept.status;
  }
  if (outOfMemory) {
    return StageStatus::kOutOfMemory;
  }

  // Each vertex's selection, grown from its core: one slot per vertex number,
  // round * frames + frame.
  const std::size_t found = rounds.size() - 1U; // the last round found nothing
  const auto vertexSlots = static_cast<Index>(found * frames);
  DeviceSpan<float> centres;
  DeviceSpan<std::uint8_t> selected;
  DeviceSpan<Counter> selectedWeightSum;
  DeviceSpan<Counter> selectedWeightedZSum;
  DeviceSpan<Counter> selectedScatterSum;
  DeviceSpan<Counter> selectedCount;
  DeviceSpan<float> mean;
  DeviceSpan<float> varMean;
  if (vertexSlots != 0U) {
    if (!centres.allocate(arena, vertexSlots) || !selected.allocate(arena, count) || !selectedWeightSum.allocate(arena, vertexSlots) ||
        !selectedWeightedZSum.allocate(arena, vertexSlots) || !selectedScatterSum.allocate(arena, vertexSlots) || !selectedCount.allocate(arena, vertexSlots) ||
        !mean.allocate(arena, vertexSlots) || !varMean.allocate(arena, vertexSlots)) {
      return StageStatus::kOutOfMemory;
    }
    for (std::size_t r = 0; r < found; ++r) {
      runtime.copyOnDevice(centres.data() + (r * frames), std::as_const(rounds[r].centre).data(), frames * sizeof(float), stream);
    }
    runtime.fill(selected.data(), 0U, count, stream);
    runtime.fill(selectedWeightSum.data(), 0U, vertexSlots * sizeof(Counter), stream);
    runtime.fill(selectedWeightedZSum.data(), 0U, vertexSlots * sizeof(Counter), stream);
    runtime.fill(selectedScatterSum.data(), 0U, vertexSlots * sizeof(Counter), stream);
    runtime.fill(selectedCount.data(), 0U, vertexSlots * sizeof(Counter), stream);
    const Selection selection{.candidates = candidates,
                              .claims = std::as_const(claims).data(),
                              .selected = selected.data(),
                              .weightSum = selectedWeightSum.data(),
                              .weightedZSum = selectedWeightedZSum.data(),
                              .scatterSum = selectedScatterSum.data(),
                              .count = selectedCount.data(),
                              .mean = mean.data(),
                              .varMean = varMean.data()};
    algorithms::forEach(runtime, stream, count, CoreStart{selection, std::as_const(centres).data(), static_cast<float>(window)});
    algorithms::forEach(runtime, stream, vertexSlots, Moments{selection});
    const auto grown = algorithms::sweepToFixedPoint<1U>(runtime, arena, stream, count + 1U, [&](Index* changed) {
      algorithms::forEach(runtime, stream, count, Grow{selection, settings.clipNSigma, changed});
      algorithms::forEach(runtime, stream, vertexSlots, Moments{selection});
    });
    if (grown.status != StageStatus::kComplete) {
      return grown.status;
    }
    // The sweep ends on Moments, so the means are the settled selection's.
    algorithms::forEach(runtime, stream, count, Scatter{selection});
  }

  // Seeds: the sums back, the seeds and each frame's order on the host.
  std::vector<Index> bestIndex(found * frames);
  std::vector<Index> contributors(found * frames);
  std::vector<Counter> weightSum(found * frames);
  std::vector<Counter> weightedZSum(found * frames);
  std::vector<Counter> scatterSum(found * frames);
  std::vector<Counter> members(found * frames);
  std::vector<Index> earliest(found * frames);
  std::vector<Index> latest(found * frames);
  for (std::size_t r = 0; r < found; ++r) {
    const std::size_t at = r * frames;
    download(rounds[r].bestIndex, bestIndex.data() + at, frames, stream);
    download(rounds[r].contributors, contributors.data() + at, frames, stream);
    download(rounds[r].earliest, earliest.data() + at, frames, stream);
    download(rounds[r].latest, latest.data() + at, frames, stream);
  }
  if (vertexSlots != 0U) {
    download(selectedWeightSum, weightSum.data(), vertexSlots, stream);
    download(selectedWeightedZSum, weightedZSum.data(), vertexSlots, stream);
    download(selectedScatterSum, scatterSum.data(), vertexSlots, stream);
    download(selectedCount, members.data(), vertexSlots, stream);
  }
  stream.synchronize();
  std::uint32_t mostInFrame = 0U;
  for (Index f = 0; f < frames; ++f) {
    struct Found {
      VertexSeed seed;
      Index contributors = 0U;
    };
    std::vector<Found> inFrame;
    for (std::size_t r = 0; r < found; ++r) {
      const std::size_t at = (r * frames) + f;
      if (bestIndex[at] == kNoIndex) {
        continue;
      }
      const VertexClaim sums{.weightSum = static_cast<std::uint64_t>(weightSum[at]),
                             .weightedZSum = static_cast<std::int64_t>(weightedZSum[at]),
                             .scatterSum = static_cast<std::uint64_t>(scatterSum[at]),
                             .count = static_cast<std::int64_t>(members[at]),
                             .earliestPS = earliest[at],
                             .latestPS = latest[at]};
      inFrame.push_back({.seed = vertexSeedOf(sums, beamSpot), .contributors = contributors[at]});
    }
    mostInFrame = std::max(mostInFrame, static_cast<std::uint32_t>(inFrame.size()));
    std::ranges::stable_sort(inFrame, {}, [](const Found& vertex) { return vertex.seed.z; });
    for (const auto& vertex : inFrame) {
      result.vertices.push_back(vertex.seed);
      result.contributors.push_back(vertex.contributors);
    }
  }
  result.rounds = mostInFrame + 1U;
  return StageStatus::kComplete;
}

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE::vertex_glue

#endif // GTRK_TRACKING_GPU_VERTEXGLUE_CUH_
