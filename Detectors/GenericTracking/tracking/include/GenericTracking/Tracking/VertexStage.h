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

/// \file VertexStage.h
/// \brief The vertexer: the vertexing pass's triplets to the timeframe's
///        vertex list -- the candidate map and the peak, mean and error
///        arithmetic shared with the device (VertexGlue.cuh), and the host's
///        rounds (`Vertexer`).
///
/// The vertexing pass is an iteration at the beam line whose links stay
/// within the innermost layers (`IterationSettings::layerMask`); it stops at
/// its triplets. Each triplet becomes a candidate: the z at the beam
/// line of the straight line through its points in r-z, if that line's chi2
/// and the transverse distance of the circle through its points to the beam
/// axis pass their cuts. The candidates are made cluster-disjoint:
/// greedily by their r-z chi2, the lower
/// triplet first on ties, a candidate is kept only if none of its clusters
/// is already kept's -- one measurement counts once, so a particle's second
/// triplet where sensors overlap, and the many triplets one fake tracklet
/// makes, do not pile into a peak or shrink a vertex's error. The kept
/// ones are sorted by (frame, z), a frame being the box start of the
/// middle measurement. Each frame then finds its vertices in rounds: the
/// free window [z_i, z_i + w) of the most candidates over its sidebands, beyond
/// the claim's reach on both sides, wins if it holds `minContributors` (a
/// frame where none does: `aloneMinContributors`, its rounds rerun) and
/// stands out of the sidebands by `significance` of their Poisson spread;
/// it claims every free candidate within the claim of its centre, and its z
/// is their mean. A found vertex owns the measurements of every triplet of
/// its frame within the claim of its centre -- those the cuts rejected too --
/// and the frame's later rounds drop the free candidates with an owned
/// measurement before they count anything: a dense
/// collision's measurements, recombined into triplets of several of its
/// particles that point elsewhere, make no second peak. Rounds repeat until
/// no window qualifies.
///
/// Everything after the candidate map is integer arithmetic -- z in
/// micrometres, counts -- or one IEEE operation on integers, so
/// the list does not depend on the order in which the work is done: not on
/// the threads here, nor on a device's schedule (VertexGlue.cuh). The map
/// itself is GPUhdi and float: the r-z line
/// fitted about its weighted mean radius, the circle's distance to the axis
/// from the axis's power with respect to it, so that neither subtracts two
/// numbers of the circle's size.

#ifndef GTRK_TRACKING_VERTEXSTAGE_H_
#define GTRK_TRACKING_VERTEXSTAGE_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/TimeEstPS.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/PairCut.h"

#ifndef GPUCA_GPUCODE_DEVICE
#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Tracking/Parallel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory_resource>
#endif
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace o2::gtrk
{

/// The vertexer's knobs.
struct VertexerSettings {
  /// The r-z straight line's chi2 on the points' z errors (1 dof) above
  /// which a triplet is no candidate; 0 for no cut.
  Scalar maxChi2RZ = 5.F;
  /// The transverse distance of the circle through a triplet's points to
  /// the beam axis above which it is no candidate [cm]; 0 for no cut.
  Scalar maxDcaCM = 0.02F;
  std::int32_t windowUM = 150;       ///< w: a peak's window [um]
  std::int32_t claimUM = 2000;       ///< the half-width around a peak's centre its vertex claims [um]
  std::uint32_t minContributors = 4; ///< the least candidates a peak's window must hold
  /// The least candidates a peak's window must hold in a frame whose rounds
  /// at `minContributors` found no vertex: such a frame runs its rounds again
  /// at this minimum, from its untouched candidates. A frame with a vertex is
  /// never rerun, so a collision's debris makes no vertex at the lower
  /// minimum. 0, or `minContributors` and above: no rerun.
  std::uint32_t aloneMinContributors = 0;
  Scalar significance = 5.F; ///< k: a peak's excess over the sidebands' mean b must reach k sqrt(b)
  /// A vertex's z and its error from its claimed candidates within this many
  /// of their own and the mean's errors of their weighted mean, grown from the
  /// core within the window: its own tails
  /// in, a neighbouring collision's candidates in the claim out. The run's
  /// `windowNSigma`.
  Scalar clipNSigma = 3.F;
};

/// z in micrometres is stored offset by this, so that the key sorts as an
/// unsigned integer.
inline constexpr std::int64_t kVertexZOffsetUM = std::int64_t{1} << 31U;

/// What the map measures of one triplet: its sort key (accepted ones), its
/// r-z line's chi2 -- the deduplication's order -- and z at the beam line
/// [cm], that z's variance from the points' z errors alone, the weighted
/// least squares' intercept [cm^2], and its circle's radius [cm], infinite
/// (1e30) for a straight line: the error model's momentum.
struct VertexCandidateFit {
  std::uint64_t key = 0U;
  Scalar chi2RZ = 0.F;
  Scalar z0 = 0.F;
  Scalar varZ0Line = 0.F;
  Scalar radius = 1.e30F;
};

/// One triplet's constraint on the beam point: the circle through its points comes closest to the beam axis at a
/// signed distance \a d along the unit direction (\a nx, \a ny), so a beam
/// point p moves that distance to d - n . p, and the point where every
/// triplet's distance vanishes obeys n . p = d. The linear term is exact for
/// a straight line and good to (|p| / R)^2 for a circle of radius R.
struct BeamTerm {
  Scalar nx = 0.F;
  Scalar ny = 0.F;
  Scalar d = 0.F;
  bool valid = false; ///< the geometry was computed: the triplet passed the r-z cut and is not degenerate
};

/// The candidate of the triplet \a a, \a b, \a c (the middle \a b's box
/// starting at \a middlePS): its sort key, (middlePS << 32) | (z in um +
/// offset), into \a key unless degenerate -- the candidates of one readout frame
/// share the frame's start -- with what else it
/// measures (`VertexCandidateFit`) into \a fit.
/// \a term, if given, receives the triplet's BeamTerm whenever its circle
/// is computed, before the transverse cut decides: a beam line off by more
/// than the cut still needs every triplet's say.
[[nodiscard]] GPUhdi() GPUAlwaysInline() VertexCandidateStatus vertexCandidate(const PairPoint& a, const PairPoint& b, const PairPoint& c, std::uint32_t middlePS, const VertexerSettings& settings, VertexCandidateFit& fit, BeamTerm* term = nullptr)
{
  using o2::gpu::GPUCommonMath;
  const PairPoint points[3] = {a, b, c};
  // The straight line z = z0 + slope r, each point weighted by 1 / its z
  // variance, fitted about the weighted means.
  Scalar weights[3] = {};
  Scalar sw = 0.F;
  Scalar swr = 0.F;
  Scalar swz = 0.F;
  for (int i = 0; i < 3; ++i) {
    weights[i] = 1.F / (points[i].varZ > 1.e-12F ? points[i].varZ : 1.e-12F);
    sw += weights[i];
    swr += weights[i] * points[i].r;
    swz += weights[i] * points[i].z;
  }
  const Scalar rMean = swr / sw;
  const Scalar zMean = swz / sw;
  Scalar srr = 0.F;
  Scalar srz = 0.F;
  for (int i = 0; i < 3; ++i) {
    const Scalar dr = points[i].r - rMean;
    srr += weights[i] * dr * dr;
    srz += weights[i] * dr * (points[i].z - zMean);
  }
  if (!(srr > 0.F)) {
    return VertexCandidateStatus::kDegenerate;
  }
  fit.varZ0Line = (1.F / sw) + (rMean * rMean / srr);
  const Scalar slope = srz / srr;
  const Scalar z0 = zMean - (slope * rMean);
  // The key, also for a candidate the cuts reject: its z owns measurements.
  // Rounded by adding one half away from zero and
  // truncating, the same on host and device (not `llround`: a value just
  // under a half can round up).
  const Scalar zUM = z0 * 1.e4F;
  const auto zIndex = static_cast<std::int64_t>(zUM < 0.F ? zUM - 0.5F : zUM + 0.5F);
  const std::int64_t zKey = zIndex + kVertexZOffsetUM;
  if (zKey < 0 || zKey >= (std::int64_t{1} << 32U)) {
    return VertexCandidateStatus::kDegenerate;
  }
  fit.key = (static_cast<std::uint64_t>(middlePS) << 32U) | static_cast<std::uint64_t>(zKey);
  Scalar chi2 = 0.F;
  for (int i = 0; i < 3; ++i) {
    const Scalar residual = points[i].z - zMean - (slope * (points[i].r - rMean));
    chi2 += weights[i] * residual * residual;
  }
  fit.chi2RZ = chi2;
  fit.z0 = z0;
  if (settings.maxChi2RZ > 0.F && chi2 > settings.maxChi2RZ) {
    return VertexCandidateStatus::kChi2RZ;
  }
  {
    // The circle through the three points in x-y, about the middle one: A and
    // C the outer points from it, U the centre from it; a straight line when
    // they are collinear.
    const Scalar ax = a.x - b.x;
    const Scalar ay = a.y - b.y;
    const Scalar cx = c.x - b.x;
    const Scalar cy = c.y - b.y;
    const Scalar twiceArea = 2.F * ((ax * cy) - (ay * cx));
    Scalar dca = 0.F;
    Scalar nx = 0.F;
    Scalar ny = 0.F;
    if (twiceArea > 1.e-12F || twiceArea < -1.e-12F) {
      const Scalar a2 = (ax * ax) + (ay * ay);
      const Scalar c2 = (cx * cx) + (cy * cy);
      const Scalar ux = ((cy * a2) - (ay * c2)) / twiceArea;
      const Scalar uy = ((ax * c2) - (cx * a2)) / twiceArea;
      const Scalar radius = GPUCommonMath::Sqrt((ux * ux) + (uy * uy));
      fit.radius = radius;
      // The centre from the axis, b + U; the axis's distance to the circle
      // |b + U| - R as its power |b + U|^2 - R^2 = |b|^2 + 2 b . U over
      // |b + U| + R.
      const Scalar centreX = b.x + ux;
      const Scalar centreY = b.y + uy;
      const Scalar centre = GPUCommonMath::Sqrt((centreX * centreX) + (centreY * centreY));
      const Scalar power = (b.x * b.x) + (b.y * b.y) + (2.F * ((b.x * ux) + (b.y * uy)));
      dca = power / (centre + radius);
      // The circle comes closest to the axis along the centre's direction.
      nx = centre > 0.F ? centreX / centre : 0.F;
      ny = centre > 0.F ? centreY / centre : 0.F;
    } else {
      const Scalar dx = c.x - a.x;
      const Scalar dy = c.y - a.y;
      const Scalar length = GPUCommonMath::Sqrt((dx * dx) + (dy * dy));
      dca = ((a.x * c.y) - (a.y * c.x)) / length;
      // A line's signed distance from a point p is d - (p x t): along (t_y, -t_x).
      nx = length > 0.F ? dy / length : 0.F;
      ny = length > 0.F ? -dx / length : 0.F;
    }
    if (term != nullptr) {
      term->nx = nx;
      term->ny = ny;
      term->d = dca;
      term->valid = ((nx * nx) + (ny * ny)) > 0.5F; // a unit direction, or none
    }
    if (settings.maxDcaCM > 0.F && (dca > settings.maxDcaCM || dca < -settings.maxDcaCM)) {
      return VertexCandidateStatus::kDca;
    }
  }
  return VertexCandidateStatus::kAccepted;
}

/// The excess of a peak starting at the free candidate at \a z [um] -- its
/// window [z, z + \a window) and the two sidebands \a claim beyond it, their
/// free candidates counted by \a countIn(low, high) -- or 0 where it is no
/// peak: fewer than \a minContributors in the window, or the excess 2W - S
/// over the sidebands' mean not positive or short of \a k2 (= 2 k^2) in
/// (2W - S)^2 >= k2 S. The host's and the device's search
/// (`Vertexer::find`, `vertex_glue::Score`).
template <typename CountIn>
[[nodiscard]] GPUhdi() GPUAlwaysInline() std::int64_t peakExcess(std::int64_t z, std::int64_t window, std::int64_t claim, std::uint32_t minContributors, Scalar k2, const CountIn& countIn)
{
  const auto inWindow = countIn(z, z + window);
  if (inWindow < minContributors) {
    return 0;
  }
  const auto sidebands = countIn(z - claim - window, z - claim) + countIn(z + window + claim, z + (2 * window) + claim);
  const auto excess = (2 * static_cast<std::int64_t>(inWindow)) - static_cast<std::int64_t>(sidebands);
  if (excess <= 0) {
    return 0;
  }
  const auto excessF = static_cast<Scalar>(excess);
  return excessF * excessF < k2 * static_cast<Scalar>(sidebands) ? 0 : excess;
}

/// The vertexer's candidate errors with the scattering of the pass's
/// material and the beam's transverse spread.
template <DetectorLayout Det>
struct ScatteringErrors {
  const ScatteringTable<Det>* table = nullptr;
  PairCut cut;
  Scalar varTransverse = 0.F;
  /// A candidate's z0 variance [cm^2]: its r-z line's (\a fit's
  /// `varZ0Line`), the scattering before and at its first point \a a at its
  /// own momentum -- the pair cut's `sourceScattering` from the fit's z0 at
  /// the pass's pT floor, scaled by the curvature: the transverse momentum
  /// goes as the circle's radius, so (R_min / R)^2 of the floor's, none for a
  /// straight line; without a field the momentum is unknown and the floor's
  /// is taken -- and the beam's transverse spread: the collision lies off the
  /// axis the line is extrapolated to, which moves z0 by the slope times
  /// that. The scattering's displacement across the track becomes one along
  /// z at the beam line by g = sqrt(1 + slope^2).
  ///
  /// \a a is taken on layer 0: the vertexing pass's links join the three
  /// innermost layers, so every triplet starts there. A pass on other layers
  /// has to take its first layer here, or its own layer's scattering is left
  /// out.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar operator()(const PairPoint& a, const VertexCandidateFit& fit) const
  {
    const VertexSeed from{.z = fit.z0};
    const Scalar atFloor = sourceScattering<Det>(*table, cut, from, LayerID{0}, a);
    const Scalar halfRMin = 0.5F * cut.twoRMin;
    const Scalar scale = halfRMin > 0.F && fit.radius > 0.F ? (halfRMin / fit.radius) * (halfRMin / fit.radius) : 1.F;
    const Scalar slope = a.r > 0.F ? (a.z - fit.z0) / a.r : 0.F;
    return fit.varZ0Line + (atFloor * scale * (1.F + (slope * slope))) + (slope * slope * varTransverse);
  }
};

/// A candidate's weight in its vertex's mean: 2^32 / its z0 variance in um^2,
/// truncated to an integer, so that the mean's sums are exact in any order.
inline constexpr Scalar kVertexWeightScale = 4294967296.F;

[[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint64_t candidateWeight(Scalar varZ0CM2)
{
  const Scalar varUM2 = varZ0CM2 * 1.e8F;
  const Scalar weight = kVertexWeightScale / (varUM2 > 1.F ? varUM2 : 1.F);
  return static_cast<std::uint64_t>(weight);
}

/// A candidate's frame and z from its key.
[[nodiscard]] GPUhdi() GPUAlwaysInline() std::int64_t vertexKeyFrame(std::uint64_t key) { return static_cast<std::int64_t>(key >> 32U); }
[[nodiscard]] GPUhdi() GPUAlwaysInline() std::int64_t vertexKeyZUM(std::uint64_t key) { return static_cast<std::int64_t>(key & 0xFFFFFFFFULL) - kVertexZOffsetUM; }

/// The ownership slot's key of measurement \a row of a triplet with
/// candidate key \a key: its frame and the row, so that a measurement in
/// triplets of two frames is owned in each apart.
[[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint64_t ownedSlotKey(std::uint64_t key, std::uint32_t row) { return (key & ~0xFFFFFFFFULL) | row; }

/// Whether a candidate at \a z [um] lies within \a halfWidth [um] of
/// \a centre: a vertex's claim, or its core.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool withinZ(std::int64_t z, Scalar centre, Scalar halfWidth)
{
  const auto zf = static_cast<Scalar>(z);
  const Scalar distance = zf > centre ? zf - centre : centre - zf;
  return distance <= halfWidth;
}

/// The weighted mean [um] of a vertex's selected candidates and its variance
/// [um^2], 1 / the sum of their weights, from their integer sums: one float
/// form for the host's rounds and the device's.
struct ZMoments {
  Scalar mean = 0.F;
  Scalar varMean = 0.F;
};

[[nodiscard]] GPUhdi() GPUAlwaysInline() ZMoments zMoments(std::uint64_t weightSum, std::int64_t weightedZSum)
{
  ZMoments moments;
  if (weightSum != 0U) {
    moments.mean = static_cast<Scalar>(weightedZSum) / static_cast<Scalar>(weightSum);
    moments.varMean = kVertexWeightScale / static_cast<Scalar>(weightSum);
  }
  return moments;
}

/// Whether a candidate at \a z [um] with \a weight lies within \a nSigma of
/// its own and the mean's errors of the mean.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool withinError(std::int64_t z, std::uint64_t weight, const ZMoments& moments, Scalar nSigma)
{
  const auto zf = static_cast<Scalar>(z);
  const Scalar distance = zf > moments.mean ? zf - moments.mean : moments.mean - zf;
  const Scalar var = (kVertexWeightScale / static_cast<Scalar>(weight > 0U ? weight : 1U)) + moments.varMean;
  return distance <= nSigma * o2::gpu::GPUCommonMath::Sqrt(var);
}

/// A selected candidate's share of its vertex's scatter about the weighted
/// \a mean [um]: w^2 (z - mean)^2 / 2^32, its weight times its chi2, as an
/// integer, so that the sum is exact in any order.
[[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint64_t candidateScatter(std::int64_t z, std::uint64_t weight, Scalar mean)
{
  const Scalar residual = static_cast<Scalar>(z) - mean;
  const auto w = static_cast<Scalar>(weight);
  return static_cast<std::uint64_t>(w * (w / kVertexWeightScale) * residual * residual);
}

/// What a candidate's claim holds besides its vertex: free, or dropped -- a
/// measurement of it owned by a vertex of its frame.
inline constexpr std::int32_t kClaimFree = -1;
inline constexpr std::int32_t kClaimDropped = -2;

/// What a vertex's z and error come from: its selected candidates' weights
/// and weighted z [um] as integer sums, their scatter about the weighted mean
/// (`candidateScatter`), their count; and its claimed candidates' boxes'
/// extent.
struct VertexClaim {
  std::uint64_t weightSum = 0U;
  std::int64_t weightedZSum = 0;
  std::uint64_t scatterSum = 0U;
  std::int64_t count = 0;
  std::uint32_t earliestPS = 0U;
  std::uint32_t latestPS = 0U;
};

/// The seed of a vertex from its \a claim: z the selected candidates'
/// weighted mean, each weight the inverse of the candidate's own z0 variance,
/// and its variance the sandwich one, sum w^2 (z - mean)^2 / (sum w)^2 times
/// n / (n - 1): the error from the candidates' scatter, honest whatever the
/// model's tails -- its box the claimed
/// ones' extent, its transverse spread \a beamSpot's. The host's rounds and
/// the device's read-back sums both end here.
[[nodiscard]] inline VertexSeed vertexSeedOf(const VertexClaim& claim, const VertexConstraint& beamSpot)
{
  const auto moments = zMoments(claim.weightSum, claim.weightedZSum);
  Scalar varMeanUM2 = moments.varMean;
  if (claim.count > 1 && claim.weightSum != 0U) {
    const auto weightSum = static_cast<Scalar>(claim.weightSum);
    const auto n = static_cast<Scalar>(claim.count);
    varMeanUM2 = (kVertexWeightScale / weightSum) * (static_cast<Scalar>(claim.scatterSum) / weightSum) * (n / (n - 1.F));
  }
  return {
    .timePS = claim.earliestPS,
    .widthPS = claim.latestPS - claim.earliestPS,
    .z = moments.mean * 1.e-4F,
    .varZ = varMeanUM2 * 1.e-8F,
    .varTransverse = static_cast<Scalar>(beamSpot.sigmaTransverse * beamSpot.sigmaTransverse),
  };
}

/// Triplet \a triplet's three rows, innermost first; the column sets are
/// anything with the tracklets' `source`/`target` and the triplets'
/// `first`/`second`, as `Vertexer::find` takes them.
template <typename Tracklets, typename Triplets>
[[nodiscard]] GPUhdi() GPUAlwaysInline() std::array<std::uint32_t, 3> tripletRows(const Tracklets& tracklets, const Triplets& triplets, std::uint32_t triplet)
{
  const auto first = triplets.first[triplet];
  return {tracklets.source[first], tracklets.target[first], tracklets.target[triplets.second[triplet]]};
}

/// One triplet as a vertex candidate: its fit (`vertexCandidate`) into
/// \a fit, its weight from \a errors (`candidateWeight` of `errors(a, fit)`,
/// \a a its first point) into \a weight and
/// its rows' time box into \a box; \a term as `vertexCandidate` takes it.
/// The host's map and a device's are this.
template <typename Measurements, typename Tracklets, typename Triplets, typename Errors>
GPUhdi() GPUAlwaysInline() VertexCandidateStatus measureCandidate(const Measurements& measurements, const Tracklets& tracklets, const Triplets& triplets, std::uint32_t triplet, const VertexerSettings& settings, const Errors& errors, VertexCandidateFit& fit, std::uint64_t& weight, TimeEstPS& box, BeamTerm* term)
{
  const auto rows = tripletRows(tracklets, triplets, triplet);
  const auto pa = pairPointAt(measurements, rows[0]);
  const auto pb = pairPointAt(measurements, rows[1]);
  const auto pc = pairPointAt(measurements, rows[2]);
  const auto status = vertexCandidate(pa, pb, pc, measurements.timePS[rows[1]], settings, fit, term);
  weight = candidateWeight(errors(pa, fit));
  box = timeBoxAt(measurements, rows[0]);
  box += timeBoxAt(measurements, rows[1]);
  box += timeBoxAt(measurements, rows[2]);
  return status;
}

/// The frames a slice's vertexer finds vertices in: those whose start -- a
/// candidate's middle measurement's box start, its key's frame -- lies in
/// `[beginPS, endPS)`, the frames whose block is in the slice's core
/// (`ownedFrames`). Every frame starts in one core, so the slices' lists
/// add up to the timeframe's, each vertex once; the frame's triplets are all
/// in the slice that owns it, its gate reaching past every link.
struct OwnedFrames {
  std::int64_t beginPS = 0;
  std::int64_t endPS = std::int64_t{1} << 32U; ///< past every box start

  [[nodiscard]] GPUhdi() GPUAlwaysInline() constexpr bool owns(std::int64_t framePS) const { return framePS >= beginPS && framePS < endPS; }
};

/// The frames whose start's block (`blockOfPS`) lies in \a core of
/// \a geometry: from the first time in its first block to the first in the
/// block after it, to the end with the timeframe's last block.
template <typename Geometry, typename Blocks>
[[nodiscard]] constexpr OwnedFrames ownedFrames(const Geometry& geometry, const Blocks& core)
{
  const auto firstPS = [&](std::int32_t block) -> std::int64_t {
    if (block >= geometry.blockCount()) {
      return OwnedFrames{}.endPS;
    }
    const std::int64_t start = geometry.blockStart(block);
    return block <= 0 || start <= 0 ? 0 : (start + geometry.partsPerPS() - 1) / geometry.partsPerPS();
  };
  return {.beginPS = firstPS(core.begin), .endPS = firstPS(core.end)};
}

/// What the vertexer found in one timeframe, or in one slice's frames.
struct VertexerResult {
  std::vector<VertexSeed> vertices;        ///< sorted by time, then z
  std::vector<std::uint32_t> contributors; ///< per vertex: the candidates it claimed
  /// `tripletVertex`'s codes other than a vertex's index: a free candidate;
  /// none (cut, degenerate or deduplicated); another slice's (its frame not
  /// owned).
  static constexpr std::int32_t kFree = -1;
  static constexpr std::int32_t kNone = -2;
  static constexpr std::int32_t kNotOwned = -3;
  /// Per triplet: the vertex that claimed it, else kFree, kNone or
  /// kNotOwned. Filled when asked for.
  std::vector<std::int32_t> tripletVertex;
  /// The triplets given and the accepted ones count every triplet of a
  /// slice, its halo's too: over slices they add up the slices' work.
  std::size_t triplets = 0U;
  std::size_t accepted = 0U;   ///< the candidates past the cuts
  std::size_t candidates = 0U; ///< those kept, cluster-disjoint, of the frames owned
  std::uint32_t rounds = 0U;   ///< the most rounds any frame took, the last finding nothing included
};

/// A triplet's rows -- its first, middle and last measurement -- and the
/// vertex that claimed it, else `VertexerResult::kFree` or `kNone`: the
/// vertexing pass's claims over all its slices.
struct TripletClaim {
  std::array<std::uint32_t, 3> rows{};
  std::int32_t vertex = VertexerResult::kNone;
};

/// A slice's vertices and counts after \a into's, as the pass gathers the
/// timeframe's slice by slice (the slices in time order, so the list stays in
/// it); the per-triplet claims are each slice's own and stay out.
inline void appendVertices(VertexerResult& into, const VertexerResult& slice)
{
  into.vertices.insert(into.vertices.end(), slice.vertices.begin(), slice.vertices.end());
  into.contributors.insert(into.contributors.end(), slice.contributors.begin(), slice.contributors.end());
  into.triplets += slice.triplets;
  into.accepted += slice.accepted;
  into.candidates += slice.candidates;
  into.rounds = into.rounds > slice.rounds ? into.rounds : slice.rounds;
}

#ifndef GPUCA_GPUCODE_DEVICE

/// The vertexer on the host. Its whole work is the candidate map, one sort
/// and the rounds of each frame, the frames in parallel in \a arena if given.
class Vertexer
{
 public:
  explicit Vertexer(const VertexerSettings& settings, const VertexConstraint& beamSpot) : mSettings(settings), mBeamSpot(beamSpot) {}

  [[nodiscard]] const VertexerSettings& settings() const { return mSettings; }

  /// The vertices of one timeframe from the vertexing pass's triplets:
  /// \a measurements the frame's columns, \a tracklets and \a triplets the
  /// pass's stage outputs (a slice's, or one over the timeframe), each with its
  /// element count. A frame is the middle measurements' box start; a vertex's
  /// box spans its claimed triplets'.
  /// \a errors gives a candidate's z0 variance from its first point and its
  /// fit (`ScatteringErrors`), its weight in its vertex's mean. The work arrays,
  /// as long as the triplets, live on \a pool (the pass's; the default
  /// resource without one), so that its limit and its peak see them;
  /// the result is the
  /// caller's. The vertices are \a owned's frames' only: a slice's, whose
  /// deduplication still sees every candidate it has, so that a frame's
  /// comes out as in one slice over the timeframe.
  template <typename Measurements, typename Tracklets, typename Triplets, typename Errors>
  [[nodiscard]] VertexerResult find(const Measurements& measurements, const Tracklets& tracklets, const Triplets& triplets, std::size_t tripletCount, detail::TaskArena* arena, bool keepTripletVertex, const Errors& errors, std::pmr::memory_resource* pool = nullptr, const OwnedFrames& owned = {}) const
  {
    std::pmr::memory_resource* const resource = pool != nullptr ? pool : std::pmr::get_default_resource();
    VertexerResult result;
    result.triplets = tripletCount;
    // Candidates: in triplet order; the key 0 marks a degenerate
    // triplet (no frame 0 candidate has z at -2^31 um), a cut one has its
    // key too. A triplet's box is its three boxes' intersection, not empty by
    // the triplet cut.
    BoundedVector<std::uint64_t> keys(tripletCount, 0U, resource);
    BoundedVector<Scalar> chi2(tripletCount, 0.F, resource);
    BoundedVector<std::uint64_t> weights(tripletCount, 0U, resource);
    BoundedVector<TimeEstPS> boxes(tripletCount, resource);
    BoundedVector<VertexCandidateStatus> status(tripletCount, VertexCandidateStatus::kDegenerate, resource);
    // body(begin, end) over [0, n) in \a arena, in pieces of \a grain, or at once.
    const auto parallel = [&](std::size_t n, std::size_t grain, const auto& body) {
      if (arena != nullptr) {
        detail::parallelFor(*arena, 0U, n, grain, body);
      } else {
        body(0U, n);
      }
    };
    const auto map = [&](std::size_t begin, std::size_t end) {
      for (auto t = begin; t < end; ++t) {
        VertexCandidateFit fit;
        status[t] = measureCandidate(measurements, tracklets, triplets, static_cast<std::uint32_t>(t), mSettings, errors, fit, weights[t], boxes[t], nullptr);
        keys[t] = fit.key;
        chi2[t] = fit.chi2RZ;
      }
    };
    parallel(tripletCount, 4096U, map);
    BoundedVector<std::uint32_t> accepted(resource);
    accepted.reserve(tripletCount);
    for (std::size_t t = 0; t < tripletCount; ++t) {
      if (status[t] == VertexCandidateStatus::kAccepted) {
        accepted.push_back(static_cast<std::uint32_t>(t));
      }
    }
    result.accepted = accepted.size();
    // The deduplication: in increasing chi2, the lower triplet first on
    // ties, each accepted candidate kept if none of its clusters is a kept
    // one's; kept ones in triplet order.
    BoundedVector<std::uint32_t> byChi2(accepted, resource);
    std::ranges::stable_sort(byChi2, {}, [&](std::uint32_t t) { return chi2[t]; });
    std::uint32_t rowEnd = 0U;
    for (const auto t : accepted) {
      rowEnd = std::max(rowEnd, std::ranges::max(tripletRows(tracklets, triplets, t)) + 1U);
    }
    BoundedVector<std::uint8_t> kept(tripletCount, 0U, resource);
    BoundedVector<std::uint8_t> taken(rowEnd, 0U, resource);
    for (const auto t : byChi2) {
      const auto rows = tripletRows(tracklets, triplets, t);
      if (std::ranges::none_of(rows, [&](std::uint32_t row) { return taken[row] != 0U; })) {
        kept[t] = 1U;
        for (const auto row : rows) {
          taken[row] = 1U;
        }
      }
    }
    // The ownership's sources: every triplet of the owned
    // frames but the degenerate, by key, and each measurement of a frame's
    // triplets one slot, (frame << 32) | row, so that frames own apart.
    BoundedVector<Source> sources(resource);
    for (std::size_t t = 0; t < tripletCount; ++t) {
      if (status[t] != VertexCandidateStatus::kDegenerate && owned.owns(vertexKeyFrame(keys[t]))) {
        sources.push_back({.key = keys[t], .triplet = static_cast<std::uint32_t>(t)});
      }
    }
    std::ranges::stable_sort(sources, {}, &Source::key);
    BoundedVector<std::uint64_t> slotKeys(resource);
    slotKeys.reserve(3U * sources.size());
    for (const auto& source : sources) {
      for (const auto row : tripletRows(tracklets, triplets, source.triplet)) {
        slotKeys.push_back(ownedSlotKey(source.key, row));
      }
    }
    std::ranges::sort(slotKeys);
    const auto repeated = std::ranges::unique(slotKeys);
    slotKeys.erase(repeated.begin(), repeated.end());
    const auto slotsOf = [&](std::uint64_t key, std::uint32_t triplet) {
      std::array<std::uint32_t, 3> slots{};
      const auto rows = tripletRows(tracklets, triplets, triplet);
      for (std::size_t k = 0; k < 3U; ++k) {
        slots[k] = static_cast<std::uint32_t>(std::ranges::lower_bound(slotKeys, ownedSlotKey(key, rows[k])) - slotKeys.begin());
      }
      return slots;
    };
    parallel(sources.size(), 4096U, [&](std::size_t begin, std::size_t end) {
      for (auto i = begin; i < end; ++i) {
        sources[i].slots = slotsOf(sources[i].key, sources[i].triplet);
      }
    });
    BoundedVector<std::uint8_t> ownedSlots(slotKeys.size(), 0U, resource);
    BoundedVector<Candidate> candidates(resource);
    candidates.reserve(accepted.size());
    for (const auto t : accepted) {
      if (kept[t] != 0U && owned.owns(vertexKeyFrame(keys[t]))) {
        candidates.push_back({.key = keys[t], .triplet = t, .box = boxes[t], .weight = weights[t], .slots = slotsOf(keys[t], t)});
      }
    }
    result.candidates = candidates.size();
    // Order: one stable sort, so that equal keys keep the triplets' order.
    std::ranges::stable_sort(candidates, {}, &Candidate::key);
    // Rounds: the frames' segments, and each frame's rounds.
    BoundedVector<std::size_t> frameBegin(resource);
    for (std::size_t i = 0; i < candidates.size(); ++i) {
      if (i == 0U || vertexKeyFrame(candidates[i].key) != vertexKeyFrame(candidates[i - 1U].key)) {
        frameBegin.push_back(i);
      }
    }
    frameBegin.push_back(candidates.size());
    const std::size_t frames = frameBegin.size() - 1U;
    BoundedVector<BoundedVector<Found>> found(frames, BoundedVector<Found>(resource), resource);
    BoundedVector<std::int32_t> claimedBy(candidates.size(), kClaimFree, resource);
    BoundedVector<std::uint32_t> rounds(frames, 0U, resource);
    const auto rounder = [&](std::size_t begin, std::size_t end) {
      for (auto f = begin; f < end; ++f) {
        const std::span<const Candidate> segment{candidates.data() + frameBegin[f], frameBegin[f + 1U] - frameBegin[f]};
        const std::span<std::int32_t> claims{claimedBy.data() + frameBegin[f], segment.size()};
        const auto frame = vertexKeyFrame(segment.front().key);
        const auto first = std::ranges::lower_bound(sources, frame, {}, [](const Source& source) { return vertexKeyFrame(source.key); });
        const auto last = std::ranges::upper_bound(sources, frame, {}, [](const Source& source) { return vertexKeyFrame(source.key); });
        rounds[f] = findInFrame(segment, claims, found[f], {first, last}, ownedSlots);
      }
    };
    parallel(frames, 1U, rounder);
    // Seeds: the list, in frame order and each frame's in z; the claims renumbered.
    BoundedVector<std::int32_t> firstOfFrame(frames, 0, resource);
    for (std::size_t f = 0; f < frames; ++f) {
      firstOfFrame[f] = static_cast<std::int32_t>(result.vertices.size());
      result.rounds = std::max(result.rounds, rounds[f]);
      for (const auto& vertex : found[f]) {
        result.vertices.push_back(vertex.seed);
        result.contributors.push_back(vertex.contributors);
      }
    }
    // Each candidate's vertex in the list: free, a dropped one too.
    BoundedVector<std::int32_t> vertexOf(candidates.size(), VertexerResult::kFree, resource);
    for (std::size_t f = 0; f < frames; ++f) {
      for (auto i = frameBegin[f]; i < frameBegin[f + 1U]; ++i) {
        if (claimedBy[i] >= 0) {
          vertexOf[i] = firstOfFrame[f] + claimedBy[i];
        }
      }
    }
    if (keepTripletVertex) {
      result.tripletVertex.assign(tripletCount, VertexerResult::kNone);
      for (std::size_t t = 0; t < tripletCount; ++t) {
        if (!owned.owns(measurements.timePS[tracklets.target[triplets.first[t]]])) {
          result.tripletVertex[t] = VertexerResult::kNotOwned;
        }
      }
      for (std::size_t i = 0; i < candidates.size(); ++i) {
        result.tripletVertex[candidates[i].triplet] = vertexOf[i];
      }
    }
    return result;
  }

 private:
  struct Candidate {
    std::uint64_t key = 0U;
    std::uint32_t triplet = 0U;
    TimeEstPS box;                        ///< the triplet's, its boxes' intersection
    std::uint64_t weight = 0U;            ///< in its vertex's mean (`candidateWeight`)
    std::array<std::uint32_t, 3> slots{}; ///< its measurements' ownership slots
  };

  /// A triplet whose measurements a vertex can own: every one of the owned
  /// frames but the degenerate, the cut ones too.
  struct Source {
    std::uint64_t key = 0U;
    std::uint32_t triplet = 0U;
    std::array<std::uint32_t, 3> slots{};
  };

  struct Found {
    VertexSeed seed;
    std::uint32_t contributors = 0U;
  };

  /// The rounds of one frame's \a segment, sorted by z, claiming into
  /// \a claims (kClaimFree, kClaimDropped, else the vertex's index in \a found);
  /// \a found sorted by z at the end. After each vertex, the measurements of
  /// the frame's \a sources (sorted by z) within its claim become owned in
  /// \a ownedSlots, and the free candidates with an owned one are dropped.
  /// A frame whose rounds find no vertex runs them again at
  /// `aloneMinContributors` when that is the lower minimum; nothing was
  /// claimed, owned or dropped, so it starts from the same state. The rounds
  /// taken, the last included -- of the rerun alone where it ran.
  std::uint32_t findInFrame(std::span<const Candidate> segment, std::span<std::int32_t> claims, BoundedVector<Found>& found, std::span<const Source> sources, std::span<std::uint8_t> ownedSlots) const
  {
    const std::size_t n = segment.size();
    std::pmr::memory_resource* const resource = found.get_allocator().resource();
    BoundedVector<std::int64_t> z(n, resource);
    for (std::size_t i = 0; i < n; ++i) {
      z[i] = vertexKeyZUM(segment[i].key);
    }
    const auto at = [&](std::int64_t value) { return static_cast<std::size_t>(std::ranges::lower_bound(z, value) - z.begin()); };
    const std::int64_t window = mSettings.windowUM;
    const std::int64_t claim = mSettings.claimUM;
    // The test (2W - S)^2 >= 2 k^2 S, with W the window's count and S the
    // two sidebands': the excess over the sidebands' mean against k times its
    // root, squared.
    const Scalar k2 = 2.F * mSettings.significance * mSettings.significance;
    BoundedVector<std::uint32_t> count(n + 1U, 0U, resource);
    std::uint32_t minContributors = mSettings.minContributors;
    std::uint32_t rounds = 0U;
    while (true) {
      ++rounds;
      for (std::size_t i = 0; i < n; ++i) {
        count[i + 1U] = count[i] + (claims[i] == kClaimFree ? 1U : 0U);
      }
      const auto countIn = [&](std::int64_t low, std::int64_t high) { return count[at(high)] - count[at(low)]; };
      std::int64_t bestExcess = 0;
      std::size_t best = n;
      for (std::size_t i = 0; i < n; ++i) {
        if (claims[i] != kClaimFree) {
          continue;
        }
        const auto excess = peakExcess(z[i], window, claim, minContributors, k2, countIn);
        if (excess <= bestExcess) {
          continue;
        }
        bestExcess = excess;
        best = i;
      }
      if (best == n) {
        if (found.empty() && mSettings.aloneMinContributors != 0U && mSettings.aloneMinContributors < minContributors) {
          minContributors = mSettings.aloneMinContributors;
          rounds = 0U;
          continue;
        }
        break;
      }
      // The centre: the mean of the window's free candidates.
      const auto windowEnd = at(z[best] + window);
      std::int64_t sum = 0;
      for (auto i = best; i < windowEnd; ++i) {
        sum += claims[i] == kClaimFree ? z[i] : 0;
      }
      const auto members = static_cast<std::int64_t>(count[windowEnd] - count[best]);
      const Scalar centre = static_cast<Scalar>(sum) / static_cast<Scalar>(members);
      // Claim within the claim's reach.
      const auto id = static_cast<std::int32_t>(found.size());
      Found vertex;
      BoundedVector<std::int64_t> claimed(resource);
      BoundedVector<std::uint64_t> claimedWeight(resource);
      std::uint32_t earliest = std::numeric_limits<std::uint32_t>::max();
      std::uint32_t latest = 0U;
      for (auto i = at(static_cast<std::int64_t>(centre) - claim - 1); i < n && static_cast<Scalar>(z[i]) <= centre + static_cast<Scalar>(claim); ++i) {
        if (claims[i] != kClaimFree || !withinZ(z[i], centre, static_cast<Scalar>(claim))) {
          continue;
        }
        claims[i] = id;
        ++vertex.contributors;
        earliest = std::min(earliest, segment[i].box.lower());
        latest = std::max(latest, segment[i].box.upper());
        claimed.push_back(z[i]);
        claimedWeight.push_back(segment[i].weight);
      }
      // The vertex's z and error: the weighted mean of a selection grown from
      // the core within the window, by the claimed candidates within
      // clipNSigma of their own and the mean's errors of it, until none joins.
      BoundedVector<std::uint8_t> selected(claimed.size(), 0U, resource);
      std::uint64_t weightSum = 0U;
      std::int64_t weightedZSum = 0;
      std::int64_t selectedCount = 0;
      const auto select = [&](std::size_t k) {
        selected[k] = 1U;
        weightSum += claimedWeight[k];
        weightedZSum += static_cast<std::int64_t>(claimedWeight[k]) * claimed[k];
        ++selectedCount;
      };
      for (std::size_t k = 0; k < claimed.size(); ++k) {
        if (withinZ(claimed[k], centre, static_cast<Scalar>(window))) {
          select(k);
        }
      }
      for (bool grew = true; grew;) {
        grew = false;
        const auto moments = zMoments(weightSum, weightedZSum);
        for (std::size_t k = 0; k < claimed.size(); ++k) {
          if (selected[k] == 0U && withinError(claimed[k], claimedWeight[k], moments, mSettings.clipNSigma)) {
            select(k);
            grew = true;
          }
        }
      }
      const auto settled = zMoments(weightSum, weightedZSum);
      std::uint64_t scatterSum = 0U;
      for (std::size_t k = 0; k < claimed.size(); ++k) {
        if (selected[k] != 0U) {
          scatterSum += candidateScatter(claimed[k], claimedWeight[k], settled.mean);
        }
      }
      vertex.seed = vertexSeedOf({.weightSum = weightSum, .weightedZSum = weightedZSum, .scatterSum = scatterSum, .count = selectedCount, .earliestPS = earliest, .latestPS = latest}, mBeamSpot);
      found.push_back(vertex);
      // The vertex owns the measurements of the frame's triplets within
      // its claim; the free candidates with an owned one are dropped before
      // the next round counts.
      const auto sourceZ = [](const Source& source) { return vertexKeyZUM(source.key); };
      const auto reach = static_cast<std::int64_t>(centre);
      for (auto it = std::ranges::lower_bound(sources, reach - claim - 1, {}, sourceZ); it != sources.end() && sourceZ(*it) <= reach + claim + 1; ++it) {
        if (withinZ(sourceZ(*it), centre, static_cast<Scalar>(claim))) {
          for (const auto slot : it->slots) {
            ownedSlots[slot] = 1U;
          }
        }
      }
      for (std::size_t i = 0; i < n; ++i) {
        if (claims[i] == kClaimFree && std::ranges::any_of(segment[i].slots, [&](std::uint32_t slot) { return ownedSlots[slot] != 0U; })) {
          claims[i] = kClaimDropped;
        }
      }
    }
    // The frame's vertices by z, ties by the order found, and the claims
    // renumbered to match.
    BoundedVector<std::uint32_t> order(found.size(), resource);
    for (std::size_t v = 0; v < order.size(); ++v) {
      order[v] = static_cast<std::uint32_t>(v);
    }
    std::ranges::stable_sort(order, {}, [&](std::uint32_t v) { return found[v].seed.z; });
    BoundedVector<std::int32_t> placeOf(found.size(), resource);
    BoundedVector<Found> sorted(resource);
    sorted.reserve(found.size());
    for (std::size_t place = 0; place < order.size(); ++place) {
      placeOf[order[place]] = static_cast<std::int32_t>(place);
      sorted.push_back(found[order[place]]);
    }
    for (auto& owner : claims) {
      owner = owner < 0 ? owner : placeOf[static_cast<std::size_t>(owner)];
    }
    found.swap(sorted);
    return rounds;
  }

  VertexerSettings mSettings;
  VertexConstraint mBeamSpot;
};

#endif

} // namespace o2::gtrk

#endif // GTRK_TRACKING_VERTEXSTAGE_H_
