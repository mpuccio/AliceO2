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

/// \file LinkHoles.h
/// \brief Whether a pair of measurements leaves a hole: the segment between
///        them crosses another layer's envelope.
///
/// A link joins two layers some trajectory crosses in a row, and one link
/// carries both kinds of pairs: ITS 0 -> 3 holds the steep tracks past the
/// ends of layers 1 and 2, and the central ones that lost their hits there.
/// Only the pair's own measurements tell which. The segment between them in
/// the (r, z) plane stands for the helix: the transverse arc outgrows r by
/// about r^3 / (24 R^2), below a millimetre over ITS 0 -> 5 at 0.1 GeV.
///
/// An iteration that takes no skips makes no tracklet that crosses a layer,
/// and the findable counts a step of a particle's chain only when it crosses
/// none: one test for both.

#ifndef GTRK_TRACKING_LINKHOLES_H_
#define GTRK_TRACKING_LINKHOLES_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/SearchWindow.h"

#include "GPUCommonMath.h"

#include <cstdint>
#include <limits>

namespace o2::gtrk
{

/// A layer's envelope as the hole test reads it: a barrel layer's radius and
/// z range, or a disc's signed z and r range [cm]; and the range of the
/// position its measurements take, radii on a barrel layer, z on a disc.
struct LayerSpan {
  Scalar position = 0.F;
  Scalar low = 0.F;
  Scalar high = 0.F;
  Scalar positionMin = 0.F;
  Scalar positionMax = 0.F;
};

/// Layer \a layer's span, from the extents the windows and the scattering
/// table read (SearchWindow.h).
template <DetectorLayout Det>
[[nodiscard]] inline LayerSpan layerSpan(LayerID layer)
{
  const bool barrel = Detector<Det>::isBarrelLayer(layer);
  return {.position = static_cast<Scalar>(barrel ? layerDistance<Det>(layer) : discZ<Det>(layer)),
          .low = static_cast<Scalar>(secondMin<Det>(layer)),
          .high = static_cast<Scalar>(secondMax<Det>(layer)),
          .positionMin = static_cast<Scalar>(barrel ? radiusMin<Det>(layer) : zMin<Det>(layer)),
          .positionMax = static_cast<Scalar>(barrel ? radiusMax<Det>(layer) : zMax<Det>(layer))};
}

/// Whether the segment from (\a ra, \a za) to (\a rb, \a zb) crosses
/// \a span, a barrel layer's when \a barrel, a disc's otherwise: it passes
/// the layer's position strictly between its ends, inside its range there.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool segmentCrosses(bool barrel, const LayerSpan& span, Scalar ra, Scalar za, Scalar rb, Scalar zb)
{
  const Scalar a = barrel ? ra : za;
  const Scalar b = barrel ? rb : zb;
  const bool between = (a < span.position && span.position < b) || (b < span.position && span.position < a);
  if (!between) {
    return false;
  }
  const Scalar t = (span.position - a) / (b - a);
  const Scalar along = barrel ? za + (t * (zb - za)) : ra + (t * (rb - ra));
  return along >= span.low && along <= span.high;
}

/// Whether the segment crosses any layer of \a holes, a bit per layer, with
/// \a barrels the barrel layers' bits and \a spans each layer's span.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool segmentLeavesHole(std::uint64_t holes, std::uint64_t barrels, const LayerSpan* spans, Scalar ra, Scalar za,
                                                                Scalar rb, Scalar zb)
{
  for (; holes != 0U; holes &= holes - 1U) {
    const std::uint32_t layer = o2::gpu::CAMath::Ctz64(holes);
    if (segmentCrosses(((barrels >> layer) & 1U) != 0U, spans[layer], ra, za, rb, zb)) {
      return true;
    }
  }
  return false;
}

namespace detail
{

/// A range of target second coordinates; empty if low > high.
struct Shadow {
  Scalar low;
  Scalar high;
};

/// The target second coordinates a segment from \a a, in the target's frame
/// (P its position axis, Q its second), to a target point at position \a w
/// crosses \a span with: \a span's shadow seen from \a a, an interval,
/// empty (low > high) when there is none. \a sameAxis: the span lies across
/// P, as a barrel layer before a barrel target; otherwise across Q.
[[nodiscard]] GPUhdi() GPUAlwaysInline() Shadow shadow(bool sameAxis, const LayerSpan& span, Scalar aP, Scalar aQ, Scalar w)
{
  constexpr Shadow kNone{.low = 1.F, .high = 0.F};
  constexpr Scalar kInfinity = std::numeric_limits<Scalar>::infinity();
  const Scalar toTarget = w - aP;
  if (sameAxis) {
    // The span at P = position over Q in [low, high], crossed when it lies
    // strictly between the source and the target: projected from the source.
    const Scalar toSpan = span.position - aP;
    const bool between = (toSpan > 0.F && toSpan < toTarget) || (toSpan < 0.F && toSpan > toTarget);
    if (!between) {
      return kNone;
    }
    const Scalar scale = toTarget / toSpan;
    return {.low = aQ + (scale * (span.low - aQ)), .high = aQ + (scale * (span.high - aQ))};
  }
  // The span at Q = position over P in [low, high]: a segment reaching the
  // target at Q = u crosses it at fraction s = (position - aQ) / (u - aQ) of
  // the way, at P = aP + s (w - aP), which has to lie in the range.
  const Scalar toSpan = span.position - aQ;
  if (exactlyEqual(toTarget, 0.F) || exactlyEqual(toSpan, 0.F)) {
    return kNone;
  }
  const Scalar s1 = (span.low - aP) / toTarget;
  const Scalar s2 = (span.high - aP) / toTarget;
  const Scalar sMin = s1 < s2 ? s1 : s2;
  const Scalar sWide = s1 < s2 ? s2 : s1;
  const Scalar sMax = sWide < 1.F ? sWide : 1.F;
  if (!(sMax > 0.F) || sMin > sMax) {
    return kNone;
  }
  const Scalar near = aQ + (toSpan / sMax);
  const Scalar unbounded = toSpan > 0.F ? kInfinity : -kInfinity;
  const Scalar far = sMin > 0.F ? aQ + (toSpan / sMin) : unbounded;
  return near < far ? Shadow{.low = near, .high = far} : Shadow{.low = far, .high = near};
}

} // namespace detail

/// The range of the second coordinate on layer \a target, within [\a low,
/// \a high] -- the target's extent --, where a pair from source (\a ra,
/// \a za) can leave no hole: the range with every shadow of a layer of
/// \a holes that covers one of its ends, at every position the target's
/// measurements take, cut off; empty (low > high) when nothing is left. A
/// hull: a shadow inside it stays, and `segmentLeavesHole` decides per pair.
/// It spares the query of a source no direct pair can leave, and the cells
/// of the others beyond it, and changes no tracklet.
[[nodiscard]] GPUhdi() GPUAlwaysInline() LayerSpan directPairRange(std::uint64_t holes, std::uint64_t barrels, const LayerSpan* spans, LayerID target,
                                                                   Scalar ra, Scalar za, Scalar low, Scalar high)
{
  const bool targetBarrel = ((barrels >> target) & 1U) != 0U;
  const Scalar aP = targetBarrel ? ra : za;
  const Scalar aQ = targetBarrel ? za : ra;
  const Scalar wLow = spans[target].positionMin;
  const Scalar wHigh = spans[target].positionMax;
  // Each layer's shadow over every target position: the shadows at the two
  // ends of the position range intersected, the shadow's ends being
  // monotonic in the position. The range's ends are cut back while a shadow
  // covers one.
  for (bool changed = true; changed && low <= high;) {
    changed = false;
    for (std::uint64_t left = holes; left != 0U; left &= left - 1U) {
      const std::uint32_t layer = o2::gpu::CAMath::Ctz64(left);
      const bool sameAxis = (((barrels >> layer) & 1U) != 0U) == targetBarrel;
      const auto atLow = detail::shadow(sameAxis, spans[layer], aP, aQ, wLow);
      const auto atHigh = detail::shadow(sameAxis, spans[layer], aP, aQ, wHigh);
      const Scalar shadowLow = atLow.low > atHigh.low ? atLow.low : atHigh.low;
      const Scalar shadowHigh = atLow.high < atHigh.high ? atLow.high : atHigh.high;
      if (atLow.low > atLow.high || atHigh.low > atHigh.high || shadowLow > shadowHigh) {
        continue;
      }
      if (shadowLow <= low && shadowHigh > low) {
        low = shadowHigh;
        changed = true;
      }
      if (shadowHigh >= high && shadowLow < high) {
        high = shadowLow;
        changed = true;
      }
    }
  }
  return {.position = spans[target].position, .low = low, .high = high, .positionMin = wLow, .positionMax = wHigh};
}

/// The layers a pair on link \a from -> \a to can cross, a bit per layer: a
/// barrel layer whose radius lies inside the radii the two layers' measurements
/// span, a disc whose z lies inside their z. A superset: `segmentCrosses`
/// decides per pair. 0 for a link with none, as l -> l + 1 on a barrel.
template <DetectorLayout Det>
[[nodiscard]] inline std::uint64_t linkHoleLayers(LayerID from, LayerID to)
{
  static_assert(Detector<Det>::nLayers <= std::numeric_limits<std::uint64_t>::digits);
  std::uint64_t holes = 0U;
  for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    if (id == from || id == to) {
      continue;
    }
    const bool barrel = Detector<Det>::isBarrelLayer(id);
    const double position = barrel ? layerDistance<Det>(id) : discZ<Det>(id);
    holes |= positionBetweenLayers<Det>(barrel, position, from, to) ? std::uint64_t{1} << layer : 0U;
  }
  return holes;
}

/// The barrel layers of \a Det, a bit per layer.
template <DetectorLayout Det>
[[nodiscard]] constexpr std::uint64_t barrelLayerBits()
{
  std::uint64_t bits = 0U;
  for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
    bits |= Detector<Det>::isBarrelLayer(static_cast<LayerID>(layer)) ? std::uint64_t{1} << layer : 0U;
  }
  return bits;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_LINKHOLES_H_
