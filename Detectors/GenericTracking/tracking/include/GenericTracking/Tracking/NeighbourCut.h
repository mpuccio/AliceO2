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

/// \file NeighbourCut.h
/// \brief The cuts two neighbouring triplets pass.
///
/// Triplets `a b c` and `b c d`, on links `a -> b`, `b -> c` and `c -> d`,
/// are neighbours when they pass, cheapest first:
///
/// 1. the agreement of the two stored curvatures;
/// 2. the boxes of `a` and `d` overlap;
/// 3. the adjacent-triplet fit's chi2, three degrees of freedom.
///
/// The curvature agreement is the caller's (`neighbourCurvatureCompatible`),
/// the others are `neighbourCompatible`'s.
///
/// Every tolerance follows from the scattering, the measurements' errors and
/// pT_min, with the iteration's `nSigma`.

#ifndef GTRK_TRACKING_NEIGHBOURCUT_H_
#define GTRK_TRACKING_NEIGHBOURCUT_H_

#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/TripletCut.h"
#include "GenericTracking/Tracking/TripletFit.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonMath.h"

namespace o2::gtrk
{

/// Both triplets fit the 3D curvature of one helix, so the two agree within
/// `n` sigma of their difference. The shared `b` and `c` correlate them
/// negatively -- `c` ends the first triplet and bends the second, `b` the
/// reverse -- so the variances' sum is not a bound; `(sigma_1 + sigma_2)^2`,
/// correlation -1, is. The adjacent fit then takes the correlation exactly.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool neighbourCurvatureCompatible(const IterationParameters<Det>& parameters, Scalar firstCurvature, Scalar firstVariance, Scalar secondCurvature, Scalar secondVariance)
{
  const Scalar difference = firstCurvature - secondCurvature;
  const Scalar nSigma = parameters.pairCut.nSigma;
  const Scalar bound = firstVariance + secondVariance + (2.F * o2::gpu::CAMath::Sqrt(firstVariance * secondVariance));
  return difference * difference <= nSigma * nSigma * bound;
}

/// The first and the fourth measurement's boxes overlap: with the triplets'
/// pairs, the four boxes share a time.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool neighbourTimeCompatible(const TimedPoint& first, const TimedPoint& fourth) { return first.time.isCompatible(fourth.time); }

namespace detail
{

/// The adjacent fit's outcome against its chi2 cut; \a fit receives it when it
/// passes.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool neighbourChi2Compatible(const IterationParameters<Det>& parameters, bool fitted, const TripletFitResult& result, TripletFitResult& fit)
{
  if (!fitted) {
    countWork(WorkCounter::kNeighbourFitFailed);
    return false;
  }
  if (!(result.chi2 <= parameters.neighbourChi2Max)) {
    countWork(WorkCounter::kNeighbourChi2);
    return false;
  }
  fit = result;
  return true;
}

} // namespace detail

/// Every cut after the curvature agreement; \a fit receives the adjacent fit
/// when all pass. `b` and `c` lie on \a secondLayer and \a thirdLayer. \a
/// secondScale and \a thirdScale multiply the scattering at `b` and `c`: the
/// iteration's kink allowance where the triplet with that middle point is
/// kinked, 1 otherwise. Both triplets' kink terms are made from the points,
/// past the time cut: `a b c`'s for this
/// candidate, `b c d`'s into \a bcd at the first candidate that needs them
/// (\a made, and \a valid with `makeKinkTerms`' result), which the caller
/// keeps for the others of the same triplet. The terms are the triplet fit's
/// own: the same function of the same points.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool neighbourCompatible(const IterationParameters<Det>& parameters, LayerID secondLayer, LayerID thirdLayer, const TimedPoint& a, const TimedPoint& b, const TimedPoint& c, const TimedPoint& d, KinkTerms& bcd, bool& made, bool& valid, TripletFitResult& fit, Scalar secondScale, Scalar thirdScale)
{
  if (!neighbourTimeCompatible(a, d)) {
    countWork(WorkCounter::kNeighbourTime);
    return false;
  }
  countWork(WorkCounter::kNeighbourFits);
  TripletFitResult result;
  bool fitted = false;
  {
    const WorkTimer timer{WorkCounter::kNeighbourFitNS};
    if (!made) {
      valid = makeKinkTerms(b.point, c.point, d.point, bcd);
      made = true;
    }
    KinkTerms abc;
    fitted = valid && makeKinkTerms(a.point, b.point, c.point, abc) && fitAdjacentTriplets(abc, bcd, a.point, b.point, c.point, d.point, secondLayer, thirdLayer, parameters.scattering, parameters.pairCut, result, secondScale, thirdScale);
  }
  return detail::neighbourChi2Compatible(parameters, fitted, result, fit);
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_NEIGHBOURCUT_H_
