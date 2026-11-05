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

/// \file MeasurementGrid.h
/// \brief The adjacent-layer candidate index.
///
/// Each layer's SoA rows are ordered by `(time block, phi cell, second cell,
/// box start)` and a prefix-sum offset table over `(block, cell)` travels with
/// them.  The second coordinate is z on a barrel layer and r on a disc, by
/// `if constexpr`.  Three keys, each doing the job it is good at:
///
///  - **the block key keeps slicing contiguous.**  A slice is a union of whole
///    blocks, so its core is one row range per layer and its halo is more
///    blocks on each side -- still one range.  That is what makes a slice's
///    bounded working set affordable, and it is why the block is the *outer* key: one
///    linear order cannot be contiguous in both time and space, and the slice
///    is the thing that needs it.
///  - **the cell key makes the spatial cut a lookup.**  A cell's rows are a
///    contiguous range found by offset, so a phi or second-coordinate window
///    selects cells instead of scanning rows.
///  - **the box start orders rows within a cell,** so the exact pairwise time
///    cut -- the boxes overlap -- is a binary search over the rows' own boxes.
///    The index is never keyed on a readout frame:
///    a frame is a property of a detector's readout, and an index that
///    assumed one would be wrong the moment two layers are staggered.
///
/// Cells are a conservative pre-filter and nothing more.  The query applies
/// the exact time cut; the exact phi and second-coordinate cuts are applied to
/// what it returns, from the azimuth, radius and z the producer ships with
/// every row (`spatiallyCompatible`).  That separation is what makes the cell
/// size a tuning parameter with no effect on the result, and it is why a
/// densely populated cell means the cells are too coarse for speed rather
/// than that the window is wrong.
///
/// The layout is built once per timeframe, by the producer, and reused across
/// slices and iterations.  An iteration with a looser minimum pT does not
/// rebuild it: a wider window simply touches more cells.

#ifndef GTRK_TRACKING_MEASUREMENTGRID_H_
#define GTRK_TRACKING_MEASUREMENTGRID_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Base/Lanes.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/MeasurementRange.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Tracking/Trig.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonAlgorithm.h"
#include "GPUCommonMath.h"
#include "MathUtils/detail/trigonometric.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>

namespace o2::gtrk
{

/// How many cells one query window spans, per axis: the derivation's
/// constants.
///
/// A window W lands anywhere on cells W / c wide: it touches about c + 1 cells
/// and reads rows over (1 + 1/c) W. The best c balances the cost of a cell
/// against the cost of a row read, which depends on the rows per window -- the
/// occupancy the grid is never sized from -- and on the backend; so c is a
/// constant, not derived: ACTS's convention, 2 in phi (one cell per one-sided
/// deflection), 1 along the second coordinate, which the 6 cm beam spot makes
/// coarse.
inline constexpr double kPhiCellsPerWindow = 2.;
inline constexpr double kSecondCellsPerWindow = 1.;

/// What a layer's cells are sized from: the widest and the narrowest window,
/// per axis, of the queries the configured iterations make into it
/// (GridSizing.h); a layer no query enters has one cell per axis. A narrowest
/// width of 0 is the widest.
struct GridLayerWindow {
  double phiHalfWidth = 0.;    ///< [rad]
  double secondHalfWidth = 0.; ///< [cm]
  bool queried = false;
  double narrowPhiHalfWidth = 0.;    ///< [rad]
  double narrowSecondHalfWidth = 0.; ///< [cm]
};

/// The grid's parameters beside the layer windows.
struct GridParams {
  /// The block edges' clock: `partsPerPS` parts make a picosecond, and the
  /// two quantities below are in parts. A bunch crossing is not a whole
  /// number of picoseconds; with 400789 parts per ps it is 10^10 parts, so
  /// edges on readout frames are exact in integers.
  std::int64_t partsPerPS = 1;

  /// Width of one time block [parts], and so the quantum of a slice; 0 for
  /// one block over the timeframe.
  ///
  /// The one parameter here that is *budgeted* rather than derived: it follows
  /// from the scratch a worker may hold, not from the detector. It is bounded
  /// above by that scratch and below by the offset table, which scales
  /// linearly in the block count. Bisection cannot go below one block, so
  /// this also sets how often the overflow ladder reaches the spill rung, and a halo rounded up to whole blocks
  /// reads more rows than the time cut needs.
  std::int64_t blockWidth = 10'000'000;

  /// Where block 0 starts [parts from the timeframe's start]; it may lie
  /// before it. A detector with readout frames puts it a picosecond before
  /// its first frame -- a box's start is its frame's rounded down -- and
  /// takes whole frames per block, so that a frame's rows are binned in one
  /// block and a zero-margin query from a block's second frame reads that
  /// block alone.
  std::int64_t blockOrigin = 0;

  /// Each layer's own cell counts, 0 for the derived count: a grid of a given
  /// shape, for the tests; no configuration
  /// key sets them.
  static constexpr std::size_t kMaxLayers = 64;
  std::int32_t phiBins[kMaxLayers]{};
  std::int32_t secondBins[kMaxLayers]{};

  /// Each layer's readout frames, of a positive length: what the slices'
  /// reach reads to tell which links cannot drift in time (SlicePlan.h's
  /// `sliceReach`).
  FrameTiming frames[kMaxLayers]{};
};

/// A run of whole time blocks, half-open.
struct BlockRange {
  std::int32_t begin = 0;
  std::int32_t end = 0;

  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::int32_t size() const { return end - begin; }
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool empty() const { return end <= begin; }
};

/// A time slice as the grid expresses it: whole blocks, core and halo.
///
/// This is the tracker's slice with its width quantised to the block grid.
/// The *result* must be invariant under slice width, not the width
/// continuous, and the halo covers cross-boundary candidates either way -- but
/// it does mean the block width is chosen when the payload is written rather
/// than when the tracker runs.
struct GridSlice {
  BlockRange core;
  BlockRange halo;
};

/// Azimuth of a global position, in [0, 2pi).
///
/// This and the transverse radius are computed once per measurement, by the
/// producer from the float position it ships, and shipped beside it: the
/// producer bins by them and the tracker cuts on them, so a row's bucket and
/// the exact cut applied to it read the same numbers, and no row lands on the
/// far side of a cell edge from where its query looks. Nothing downstream
/// evaluates either again.
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar azimuth(Scalar x, Scalar y)
{
  const Scalar phi = atan2Of(y, x);
  return o2::math_utils::detail::to02Pi<Scalar>(phi);
}

/// The second binning coordinate, picked from the cached ones: z on a barrel
/// layer, r on a disc.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar secondCoordinate(LayerID layer, Scalar r, Scalar z)
{
  return Detector<Det>::isBarrelLayer(layer) ? z : r;
}

/// The bin counts one layer's axes carry, and the affine maps onto them.
///
/// Public aggregate and trivially copyable on purpose: this is uploaded to a
/// device as-is.
struct LayerAxes {
  std::int32_t nPhi = 1;
  std::int32_t nSecond = 1;
  Scalar secondMin = 0.F;      ///< lower edge of the second-coordinate axis [cm]
  Scalar secondInvWidth = 0.F; ///< cells per cm along it
  Scalar phiInvWidth = 0.F;    ///< cells per radian

  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::int32_t cellCount() const { return nPhi * nSecond; }
};

/// Bin count for an axis: as many cells as it takes to cover \a extent with
/// cells \a windowWidth / \a cellsPerWindow wide, or 0 when that is more
/// than an int32 holds -- checked before the conversion, so it never wraps;
/// the geometry then refuses the layout.
///
/// A window wider than the axis, or an axis that a detector does not resolve
/// at all, collapses to one cell rather than to none -- the axis always
/// exists, and a weak constraint yields few cells rather than no binning.
[[nodiscard]] inline std::int32_t derivedBinCount(double extent, double windowWidth, double cellsPerWindow)
{
  if (!(extent > 0.) || !(windowWidth > 0.) || !(cellsPerWindow > 0.)) {
    return 1;
  }
  const double count = std::ceil(extent * cellsPerWindow / windowWidth);
  if (!(count > 1.)) {
    return 1;
  }
  return count <= static_cast<double>(INT32_MAX) ? static_cast<std::int32_t>(count) : 0;
}

/// Bin count for an axis searched by windows from \a narrowWidth to
/// \a wideWidth: one cell width for both
/// ends, each held to the `cellsPerWindow` convention. The narrowest window
/// reads rows over its width plus a cell, (1 + w / a) against the
/// convention's (1 + 1 / c); the widest spans b / w + 1 cells against c + 1.
/// The width at which both exceed the convention by the same factor solves
/// (c / a) w^2 + (c - 1) w - b = 0; with one width, a = b, it is a / c, the
/// convention itself. Every window between exceeds it by less.
[[nodiscard]] inline std::int32_t balancedBinCount(double extent, double narrowWidth, double wideWidth, double cellsPerWindow)
{
  if (!(narrowWidth > 0.) || !(narrowWidth < wideWidth) || !(cellsPerWindow > 0.)) {
    return derivedBinCount(extent, wideWidth, cellsPerWindow);
  }
  const double c = cellsPerWindow;
  const double a = narrowWidth;
  const double b = wideWidth;
  const double width = a * (std::sqrt(((c - 1.) * (c - 1.)) + (4. * c * b / a)) - (c - 1.)) / (2. * c);
  return derivedBinCount(extent, width * c, c);
}

/// The bin of \a scaled -- a coordinate already divided by the bin width -- on
/// an axis of \a count bins: the first bin for NaN or anything not above zero,
/// the last for anything at or past the end. The range is checked before the
/// conversion because converting a float outside int32's range is undefined:
/// arm64 and the GPUs saturate, while x86 returns INT_MIN, which passes a
/// `bin < count` clamp and indexes far outside the grid.
template <typename T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() constexpr std::int32_t clampedBin(T scaled, std::int32_t count)
{
  if (!(scaled > T{0})) {
    return 0;
  }
  if (!(scaled < static_cast<T>(count))) {
    return count - 1;
  }
  return static_cast<std::int32_t>(scaled);
}

/// Where every cell of every layer sits, and how a `(time, phi, second)`
/// triple becomes a bucket number.
///
/// Built on the host from `Det`, the timeframe duration and the parameters
/// above; a plain array of PODs afterwards, so it copies to a device with the
/// rest of the tracker's configuration.
///
/// Bucket numbering is layer-major, then block, then phi, then the second
/// coordinate. That order is the whole point:
///
///  - layer-major means the grid's prefix offsets also provide each whole
///    layer's extent;
///  - block before cell means a slice's rows are one contiguous range;
///  - the second coordinate innermost means a query's `(block, phi)` pair
///    contributes one contiguous run of cells rather than a stride.
template <DetectorLayout Det>
class GridGeometry
{
 public:
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);

  GridGeometry() = default;

  /// \param durationPS  the timeframe's span; the block grid covers exactly it
  /// \param params      the block clock, and explicit counts for tests
  /// \param windows     per layer, the widest and the narrowest window of the
  ///                    queries into it (GridSizing.h's `gridLayerWindows`)
  GridGeometry(std::int64_t durationPS, const GridParams& params, std::span<const GridLayerWindow> windows)
    : mPartsPerPS(params.partsPerPS),
      mBlockOrigin(params.blockOrigin),
      mBlockWidth(chosenBlockWidth((durationPS * params.partsPerPS) - params.blockOrigin, params)),
      mBlockCount(blocksCovering((durationPS * params.partsPerPS) - params.blockOrigin, mBlockWidth))
  {
    static_assert(kLayerCount <= GridParams::kMaxLayers);
    assert(windows.size() == kLayerCount);
    // The boxes are 32-bit picoseconds from the timeframe's start (TimeEstPS):
    // a longer timeframe would wrap them silently.
    if (std::cmp_greater(durationPS, std::numeric_limits<std::uint32_t>::max())) {
      throw std::invalid_argument("the timeframe is longer than the 32-bit picosecond time box holds (4.29 ms)");
    }
    std::uint64_t base = 0;
    for (std::size_t layer = 0; layer < kLayerCount; ++layer) {
      if (params.frames[layer].lengthBC <= 0) {
        throw std::invalid_argument("every layer of the grid has readout frames of a positive length");
      }
      mFrames[layer] = params.frames[layer];
      const auto id = static_cast<LayerID>(layer);
      const auto& window = windows[layer];
      const double secondMinValue = secondMin<Det>(id);
      const double extent = secondMax<Det>(id) - secondMinValue;
      // A layer no query enters needs no cells: one per axis.
      const auto derived = [&](std::int32_t set, double axisExtent, double halfWidth, double cellsPerWindow) {
        if (set > 0) {
          return set;
        }
        return window.queried ? derivedBinCount(axisExtent, 2. * halfWidth, cellsPerWindow) : 1;
      };

      // Phi: balanced between the narrowest and the widest query, whose phi
      // cells are a range each. The second coordinate: by the narrowest, as
      // its cells are the innermost key -- a wider query spans more of them in
      // the same range per phi cell.
      const double narrowPhi = window.narrowPhiHalfWidth > 0. ? window.narrowPhiHalfWidth : window.phiHalfWidth;
      const double narrowSecond = window.narrowSecondHalfWidth > 0. ? window.narrowSecondHalfWidth : window.secondHalfWidth;
      auto& axes = mAxes[layer];
      axes.nPhi = params.phiBins[layer];
      if (axes.nPhi <= 0) {
        axes.nPhi = window.queried ? balancedBinCount(2 * std::numbers::pi_v<double>, 2. * narrowPhi, 2. * window.phiHalfWidth, kPhiCellsPerWindow) : 1;
      }
      axes.nSecond = derived(params.secondBins[layer], extent, narrowSecond, kSecondCellsPerWindow);
      // A count past int32 comes back as 0: the layout is refused, and the
      // axis kept at one cell so nothing below divides by or wraps it.
      if (axes.nPhi <= 0 || axes.nSecond <= 0) {
        mOverflow = true;
        axes.nPhi = std::max(axes.nPhi, 1);
        axes.nSecond = std::max(axes.nSecond, 1);
      }
      axes.phiInvWidth = static_cast<Scalar>(static_cast<double>(axes.nPhi) / (2 * std::numbers::pi_v<double>));
      axes.secondMin = static_cast<Scalar>(secondMinValue);
      axes.secondInvWidth = extent > 0. ? static_cast<Scalar>(static_cast<double>(axes.nSecond) / extent) : 0.F;

      mLayerBase[layer] = static_cast<std::uint32_t>(base);
      // A layer's buckets, checked before they are formed: each count is at
      // most 2^31, so the cells of a block fit 64 bits, and if they fit a
      // signed 32-bit offset so do the blocks' times them; bucket() computes
      // block * cells + cell in 32-bit signed arithmetic.
      const std::uint64_t cells = static_cast<std::uint64_t>(axes.nPhi) * static_cast<std::uint64_t>(axes.nSecond);
      const std::uint64_t layerBuckets = cells <= std::uint64_t{INT32_MAX} ? static_cast<std::uint64_t>(mBlockCount) * cells : std::uint64_t{INT32_MAX} + 1U;
      mOverflow = mOverflow || layerBuckets > std::uint64_t{INT32_MAX};
      base = std::min(base + layerBuckets, std::uint64_t{UINT32_MAX} + 1U);
    }
    mLayerBase[kLayerCount] = static_cast<std::uint32_t>(std::min(base, std::uint64_t{UINT32_MAX}));
    mOverflow = mOverflow || base > std::uint64_t{UINT32_MAX};
  }

  /// Whether a layer's buckets outnumber a 32-bit signed bucket offset, or all
  /// of them a 32-bit bucket number: the layout is then unusable, and
  /// makeGridGeometry refuses it.
  [[nodiscard]] bool overflows() const { return mOverflow; }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::int32_t blockCount() const { return mBlockCount; }
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::int64_t partsPerPS() const { return mPartsPerPS; }
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::int64_t blockWidth() const { return mBlockWidth; }
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::int64_t blockOrigin() const { return mBlockOrigin; }
  /// Layer \a layer's readout frames; none without.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() const FrameTiming& frames(LayerID layer) const { return mFrames[layer]; }
  /// Where block \a block starts [parts from the timeframe's start].
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::int64_t blockStart(std::int32_t block) const
  {
    return mBlockOrigin + (static_cast<std::int64_t>(block) * mBlockWidth);
  }
  [[nodiscard]] GPUhdi() GPUAlwaysInline() const LayerAxes& axes(LayerID layer) const { return mAxes[layer]; }
  [[nodiscard]] GPUhdi() GPUAlwaysInline() BlockRange allBlocks() const { return {.begin = 0, .end = mBlockCount}; }

  /// Number of `(block, cell)` buckets over the whole detector. The offset
  /// table has one more entry than this.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t bucketCount() const { return mLayerBase[kLayerCount]; }

  /// First bucket of a layer. `offsets[layerBucketBase(l)]` is that layer's
  /// first row, which is how a whole-layer range is recovered.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t layerBucketBase(LayerID layer) const { return mLayerBase[layer]; }

  /// The block of a time box's start [ps from the timeframe's start], exact;
  /// a query's bound before the first block lands in it, one past the last
  /// in the last.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::int32_t blockOfPS(std::int64_t timePS) const
  {
    const std::int64_t parts = (timePS * mPartsPerPS) - mBlockOrigin;
    if (parts < 0) {
      return 0;
    }
    const std::int64_t block = parts / mBlockWidth;
    return block < mBlockCount ? static_cast<std::int32_t>(block) : mBlockCount - 1;
  }

  /// \a phi must already be in [0, 2pi).
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::int32_t phiBinOf(LayerID layer, Scalar phi) const
  {
    const auto& axes = mAxes[layer]; // NOLINT(clang-analyzer-security.ArrayBound): layer < kLayerCount, the geometry's
    return clampedBin(phi * axes.phiInvWidth, axes.nPhi);
  }

  /// Clamped, so a coordinate off the end of a layer lands in its edge cell.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::int32_t secondBinOf(LayerID layer, Scalar second) const
  {
    const auto& axes = mAxes[layer];
    return clampedBin((second - axes.secondMin) * axes.secondInvWidth, axes.nSecond);
  }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t bucket(LayerID layer, std::int32_t block, std::int32_t phiBin, std::int32_t secondBin) const
  {
    const auto& axes = mAxes[layer];
    const auto cell = phiBin * axes.nSecond + secondBin;
    return mLayerBase[layer] + static_cast<std::uint32_t>(block * axes.cellCount() + cell);
  }

  /// The bucket a measurement belongs to, from its own coordinates: the
  /// block of its time box's start.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t bucketAt(LayerID layer, std::uint32_t timePS, Scalar phi, Scalar second) const
  {
    return bucket(layer, blockOfPS(timePS), phiBinOf(layer, phi), secondBinOf(layer, second));
  }

  /// Blocks a slice of \a widthPS covers, 0 for one slice over the
  /// timeframe; at least one block.
  [[nodiscard]] std::int32_t blocksPerSlice(std::int64_t widthPS) const
  {
    if (widthPS <= 0) {
      return mBlockCount;
    }
    return static_cast<std::int32_t>(std::clamp<std::int64_t>(ceilDiv(widthPS * mPartsPerPS, mBlockWidth), 1, mBlockCount));
  }

  /// Widen a core block range by a halo of \a haloPS, in whole blocks,
  /// clipped at the ends of the timeframe -- the hard input boundary of the
  /// slicing.
  [[nodiscard]] GridSlice slice(BlockRange core, std::int64_t haloPS) const
  {
    const auto pad = static_cast<std::int32_t>(std::clamp<std::int64_t>(ceilDiv(std::max<std::int64_t>(haloPS, 0) * mPartsPerPS, mBlockWidth), 0, mBlockCount));
    const auto clamp = [this](std::int32_t block) {
      if (block < 0) {
        return 0;
      }
      return block > mBlockCount ? mBlockCount : block;
    };
    return {.core = {clamp(core.begin), clamp(core.end)}, .halo = {clamp(core.begin - pad), clamp(core.end + pad)}};
  }

 private:
  /// \a a / \a b rounded up, \a a >= 0, \a b > 0.
  [[nodiscard]] static constexpr std::int64_t ceilDiv(std::int64_t a, std::int64_t b) { return (a + b - 1) / b; }

  /// A block width of zero means one block over the whole timeframe, which is
  /// what a tracker that does not want to slice at all asks for.
  [[nodiscard]] static std::int64_t chosenBlockWidth(std::int64_t span, const GridParams& params)
  {
    const std::int64_t width = params.blockWidth > 0 ? params.blockWidth : span;
    return width > 0 ? width : 1;
  }

  [[nodiscard]] static std::int32_t blocksCovering(std::int64_t span, std::int64_t width)
  {
    return span > width ? static_cast<std::int32_t>(ceilDiv(span, width)) : 1;
  }

  // Declaration order is initialisation order, and the block count is derived
  // from the width -- so the width comes first here, not for tidiness.
  std::int64_t mPartsPerPS = 1;
  std::int64_t mBlockOrigin = 0;
  std::int64_t mBlockWidth = 1;
  std::int32_t mBlockCount = 1;
  LayerAxes mAxes[kLayerCount]{};
  std::uint32_t mLayerBase[kLayerCount + 1]{};
  FrameTiming mFrames[kLayerCount]{};
  bool mOverflow = false;
};

/// One seed's question, in the coordinates the grid bins.
///
/// Everything here is arithmetic on a `LayerPairWindow` the host computed once
/// for this layer pair and this iteration: the seed contributes its time box,
/// an azimuth and its r and z, and nothing transcendental happens per seed.
///
/// The time is a closed range of the candidates' box starts, `timePS`, and a
/// time the candidates' boxes must end after: the rows a box overlaps.
/// The range is signed and may begin before
/// the timeframe.
struct CandidateWindow {
  std::int64_t earliestPS = 0;
  std::int64_t latestPS = 0;
  std::uint32_t afterPS = 0U; ///< a candidate's box must end after this
  Scalar phiCentre = 0.F;     ///< in [0, 2pi)
  Scalar phiHalfWidth = 0.F;
  Scalar secondCentre = 0.F;
  Scalar secondHalfWidth = 0.F;
};

/// A range of the target's second coordinate, as the grid queries it.
struct SecondRange {
  Scalar centre = 0.F;
  Scalar halfWidth = 0.F;
};

struct SecondRun {
  Scalar low = 0.F;
  Scalar high = 0.F;
};

/// The target's second coordinate a seed can reach across a barrel <-> disc
/// pair, from the fields `crossingPairWindow`
/// stored.
///
/// The landing point is monotonic in the vertex z, in the target's position
/// and in the seed's own coordinate separately, so its extremes lie at the
/// corners of those three ranges. Barrel -> disc divides by `z_i - z_v`: a
/// vertex level with the seed, or beyond it as seen from the disc, sends no
/// track through the seed to the disc, and the landing radius grows without
/// bound as the vertex approaches the seed. Such a corner opens the window to
/// the disc's outer edge. With no corner valid, the whole extent is returned,
/// which the seed's source interval normally prevents from being asked.
[[nodiscard]] GPUhdi() GPUAlwaysInline() SecondRange crossingRange(const LayerPairWindow& pair, Scalar seedR, Scalar seedZ, Scalar vertexLow, Scalar vertexHigh,
                                                                   Scalar padding)
{
  const Scalar vertices[2] = {vertexLow, vertexHigh};
  const Scalar targets[2] = {pair.targetLow, pair.targetHigh};
  const Scalar spreads[2] = {-pair.seedSpread, pair.seedSpread};
  const bool toDisc = pair.kind == PairKind::kBarrelToDisc;
  Scalar low = pair.targetSecondMax;
  Scalar high = pair.targetSecondMin;
  bool open = false;
  for (const Scalar spread : spreads) {
    for (const Scalar vertex : vertices) {
      for (const Scalar target : targets) {
        Scalar landing = 0.F;
        if (toDisc) {
          const Scalar z = seedZ + spread;
          const Scalar lever = z - vertex;
          if (!(lever * (target - z) > 0.F)) {
            open = true;
            continue;
          }
          landing = seedR * (target - vertex) / lever;
        } else {
          landing = vertex + ((seedZ - vertex) * target / (seedR + spread));
        }
        low = landing < low ? landing : low;
        high = landing > high ? landing : high;
      }
    }
  }
  if (open) {
    high = pair.targetSecondMax;
  }
  if (high < low) {
    low = pair.targetSecondMin;
    high = pair.targetSecondMax;
  }
  return {.centre = 0.5F * (low + high), .halfWidth = (0.5F * (high - low)) + padding};
}

[[nodiscard]] GPUhdi() GPUAlwaysInline() SecondRange crossingRange(const LayerPairWindow& pair, Scalar seedR, Scalar seedZ)
{
  return crossingRange(pair, seedR, seedZ, pair.vertexLow, pair.vertexHigh, pair.secondHalfWidth);
}

/// One source/link/vertex interval. Its host-derived half-width is shared by
/// the timeframe's vertices; only the centre and the vertex's own z reach are
/// evaluated here.
[[nodiscard]] GPUhdi() GPUAlwaysInline() SecondRange vertexSecondRange(const LayerPairWindow& pair, const VertexPairWindow& bounds, Scalar seedR, Scalar seedZ,
                                                                       Scalar vertexZ, Scalar vertexVarZ, Scalar nSigma)
{
  if (pair.kind != PairKind::kBarrel && pair.kind != PairKind::kDisc) {
    const Scalar reach = nSigma * o2::gpu::CAMath::Sqrt(vertexVarZ);
    return crossingRange(pair, seedR, seedZ, vertexZ - reach, vertexZ + reach, bounds.secondPadding);
  }
  const Scalar vertexTerm = bounds.vertexLever * bounds.vertexLever * vertexVarZ;
  const Scalar halfWidth = (nSigma * o2::gpu::CAMath::Sqrt(vertexTerm + bounds.secondVariance)) + bounds.secondPadding;
  if (pair.kind == PairKind::kBarrel) {
    return {.centre = vertexZ + ((seedZ - vertexZ) * bounds.secondScale), .halfWidth = halfWidth};
  }
  const Scalar centre = seedR * (bounds.targetPosition - vertexZ) / (bounds.sourcePosition - vertexZ);
  return {.centre = centre, .halfWidth = halfWidth};
}

[[nodiscard]] GPUhdi() GPUAlwaysInline() bool secondRangeContains(const SecondRange& range, Scalar second)
{
  return o2::gpu::CAMath::Abs(second - range.centre) <= range.halfWidth;
}

/// Assemble one from a precomputed pair window and a seed.
///
/// \a seedR and \a seedZ are the seed's shipped radius and z: a same-kind pair
/// reads the seed's second coordinate from them, a crossing pair both.
/// \a maxCandidateWidthPS must bound the boxes' widths over the searched rows:
/// a box overlapping the seed's starts less than its own width before it.
[[nodiscard]] GPUhdi() GPUAlwaysInline() CandidateWindow candidateWindow(const LayerPairWindow& pair, const TimeEstPS& seed, Scalar seedPhi, Scalar seedR,
                                                                         Scalar seedZ, std::uint32_t maxCandidateWidthPS)
{
  SecondRange second;
  if (pair.kind == PairKind::kBarrel || pair.kind == PairKind::kDisc) {
    const Scalar seedSecond = pair.kind == PairKind::kBarrel ? seedZ : seedR;
    second = {.centre = (pair.secondScale * seedSecond) + pair.secondOffset, .halfWidth = pair.secondHalfWidth};
  } else {
    second = crossingRange(pair, seedR, seedZ);
  }
  return {.earliestPS = static_cast<std::int64_t>(seed.lower()) - static_cast<std::int64_t>(maxCandidateWidthPS) + 1,
          .latestPS = static_cast<std::int64_t>(seed.upper()) - 1,
          .afterPS = seed.lower(),
          .phiCentre = seedPhi,
          .phiHalfWidth = pair.phiHalfWidth,
          .secondCentre = second.centre,
          .secondHalfWidth = second.halfWidth};
}

/// Apply the exact spatial cuts after the grid cells have supplied a
/// conservative candidate set. Both azimuths must be in [0, 2pi).
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool spatiallyCompatible(const CandidateWindow& window, Scalar phi, Scalar second)
{
  const Scalar deltaPhi = o2::gpu::GPUCommonMath::Abs(o2::math_utils::detail::toPMPi<Scalar>(phi - window.phiCentre));
  const Scalar deltaSecond = o2::gpu::GPUCommonMath::Abs(second - window.secondCentre);
  return deltaPhi <= window.phiHalfWidth && deltaSecond <= window.secondHalfWidth;
}

/// The cells one query visits: a run of whole blocks, at most two runs of phi
/// cells split where the window crosses zero, and one run of second-coordinate
/// cells, plus the box-start bracket searched inside each.
struct CellSpan {
  std::int64_t earliestPS = 0;
  std::int64_t latestPS = 0;
  std::uint32_t afterPS = 0U;
  BlockRange blocks;
  std::int32_t firstPhi = 0; ///< first run: [firstPhi, phiEnd)
  std::int32_t phiEnd = 0;
  std::int32_t wrapEnd = 0; ///< second run, after the wrap: [0, wrapEnd)
  std::int32_t firstSecond = 0;
  std::int32_t lastSecond = -1; ///< inclusive

  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::int32_t cellCount() const
  {
    return blocks.empty() ? 0 : blocks.size() * (phiEnd - firstPhi + wrapEnd) * (lastSecond - firstSecond + 1);
  }
};

/// The index, bound to the columns whose order it describes.
///
/// It borrows both: the measurement columns and the offset table belong to the
/// incoming message, exactly as `TimeFrame<Det>` borrows its view.
template <DetectorLayout Det>
class MeasurementGrid
{
 public:
  using Index = std::uint32_t;

  MeasurementGrid(MeasurementView measurements, const GridGeometry<Det>& grid, std::span<const Index> offsets)
    : mMeasurements(measurements), mGrid(grid), mOffsets(offsets.data())
  {
    assert(offsets.size() == static_cast<std::size_t>(grid.bucketCount()) + 1U);
  }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() const GridGeometry<Det>& geometry() const { return mGrid; }
  /// The prefix offsets the grid borrows, `bucketCount() + 1` of them.
  [[nodiscard]] std::span<const Index> offsets() const { return {mOffsets, static_cast<std::size_t>(mGrid.bucketCount()) + 1U}; }

  /// The rows of a run of whole blocks: one range, which is the property that
  /// makes the block the outer key.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() MeasurementRange rows(LayerID layer, BlockRange blocks) const
  {
    if (blocks.empty()) {
      const auto at = mGrid.bucket(layer, blocks.begin, 0, 0);
      return {mOffsets[at], mOffsets[at]};
    }
    const auto first = mGrid.bucket(layer, blocks.begin, 0, 0);
    const auto lastBucket = mGrid.layerBucketBase(layer) + static_cast<std::uint32_t>(blocks.end) * static_cast<std::uint32_t>(mGrid.axes(layer).cellCount());
    return {mOffsets[first], mOffsets[lastBucket]};
  }

  /// The whole layer, recovered from the grid prefix offsets.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() MeasurementRange layerRows(LayerID layer) const { return rows(layer, mGrid.allBlocks()); }
  /// The measurement columns the grid orders.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() MeasurementView measurements() const { return mMeasurements; }

  /// Visit every row of \a target that could continue the seed \a window
  /// describes, restricted to \a blocks -- a slice's halo range.
  ///
  /// What the visitor gets is a global row index, which is the tracking's
  /// internal measurement ID; `Measurement::clusterIndex` remains the
  /// back-reference to the input cluster.
  ///
  /// The exact time cut -- the boxes overlap -- is applied here. The exact phi and
  /// second-coordinate cuts are the caller's, on the row's shipped azimuth and
  /// second coordinate: what the cells did is narrow the rows it has to
  /// evaluate.
  ///
  /// Rows are stored in the order this visits them -- block, phi cell,
  /// second-coordinate cell, box start -- so a run of consecutive cells is one
  /// range of rows. Where the window spans the whole second axis, a block's
  /// phi run is one range; otherwise each phi cell's run of second cells is.
  /// Either way the rows and their order are those of a cell-by-cell walk.
  /// With \a lanes (`Lanes`), this lane's share of them, in the same order,
  /// for a visitor without collectives.
  ///
  /// With a group (`WarpGroup`) the lanes walk in
  /// step: lane `i` of `L` gets the rows `Lanes{i, L}` would, in the same
  /// order, and every lane calls \a visit(row, has) at every step of the
  /// group -- \a has whether the lane has a row there, \a row meaningful only
  /// then -- the same number of times in every range and cell, so that
  /// \a visit may call the group's collectives.
  template <typename Visitor, typename Share = Lanes>
  GPUhdi() void visitCandidates(LayerID target, const CandidateWindow& window, BlockRange blocks, Visitor visit, Share lanes = {}) const
  {
    const auto span = querySpan(target, window, blocks);
    const auto lastSecondBin = mGrid.axes(target).nSecond - 1;
    const bool wholeSecond = span.firstSecond == 0 && span.lastSecond == lastSecondBin;
    for (auto block = span.blocks.begin; block < span.blocks.end; ++block) {
      if (wholeSecond) {
        if (span.firstPhi < span.phiEnd) {
          visitBuckets(mGrid.bucket(target, block, span.firstPhi, 0), mGrid.bucket(target, block, span.phiEnd - 1, lastSecondBin) + 1U, span, visit, lanes);
        }
        if (span.wrapEnd > 0) {
          visitBuckets(mGrid.bucket(target, block, 0, 0), mGrid.bucket(target, block, span.wrapEnd - 1, lastSecondBin) + 1U, span, visit, lanes);
        }
        continue;
      }
      forEachPhiBin(span, [&](std::int32_t phiBin) { visitSecondRun(target, block, phiBin, span.firstSecond, span.lastSecond, span, visit, lanes); });
    }
  }

  /// Visit ordered, disjoint second-coordinate runs without visiting a row
  /// twice when different runs clamp to or overlap the same grid cell.
  template <typename Visitor>
  GPUhdi() void visitCandidates(LayerID target, const CandidateWindow& window, const SecondRun* secondRuns, std::uint32_t secondRunCount, BlockRange blocks,
                                Visitor visit, Lanes lanes = {}) const
  {
    const auto span = querySpan(target, window, blocks);
    for (auto block = span.blocks.begin; block < span.blocks.end; ++block) {
      forEachPhiBin(span, [&](std::int32_t phiBin) { visitRuns(target, block, phiBin, span, secondRuns, secondRunCount, visit, lanes); });
    }
  }

  /// The cells `visitCandidates` walks for this window.
  ///
  /// Public so that the study can count a query's cost (`CompleteAnatomy.h`);
  /// the result never depends on the cells.
  [[nodiscard]] GPUhdi() CellSpan cellSpan(LayerID target, const CandidateWindow& window, BlockRange blocks) const
  {
    const auto& axes = mGrid.axes(target);

    CellSpan span;
    span.earliestPS = window.earliestPS;
    span.latestPS = window.latestPS;
    span.afterPS = window.afterPS;

    const auto firstBlock = mGrid.blockOfPS(span.earliestPS);
    const auto lastBlock = mGrid.blockOfPS(span.latestPS);
    span.blocks.begin = firstBlock > blocks.begin ? firstBlock : blocks.begin;
    span.blocks.end = lastBlock < blocks.end - 1 ? lastBlock + 1 : blocks.end;
    if (span.blocks.empty()) {
      return span;
    }

    // No early return for a window past either end of the axis: the edge
    // cells are clamped, not bounded. A disc's outer cell holds clusters out
    // to its sensors' corners, and a window lying wholly beyond the nominal
    // edge can still contain them. The exact cut rejects what does not fit.
    span.firstSecond = mGrid.secondBinOf(target, window.secondCentre - window.secondHalfWidth);
    span.lastSecond = mGrid.secondBinOf(target, window.secondCentre + window.secondHalfWidth);

    // The phi cells as at most two runs, split where the window crosses zero:
    // [firstPhi, phiEnd) and then [0, wrapEnd). Where the wrap falls is known
    // before the loop, so the walk needs neither a modulo nor a per-cell
    // wrap test.
    span.phiEnd = axes.nPhi;
    if (2.F * window.phiHalfWidth < (2 * std::numbers::pi_v<Scalar>)) {
      Scalar low = window.phiCentre - window.phiHalfWidth;
      Scalar high = window.phiCentre + window.phiHalfWidth;
      if (low < 0.F) {
        low += 2 * std::numbers::pi_v<Scalar>;
      }
      if (high >= 2 * std::numbers::pi_v<Scalar>) {
        high -= 2 * std::numbers::pi_v<Scalar>;
      }
      const auto firstPhi = mGrid.phiBinOf(target, low);
      const auto lastPhi = mGrid.phiBinOf(target, high);
      if (!(high < low)) {
        span.firstPhi = firstPhi;
        span.phiEnd = lastPhi + 1;
      } else if (lastPhi < firstPhi) {
        span.firstPhi = firstPhi;
        span.wrapEnd = lastPhi + 1;
      }
      // Otherwise the window crosses zero with both ends in one cell: the gap
      // it leaves lies inside that cell, so every cell is visited.
    }
    return span;
  }

 private:
  /// A range scanned row by row costs what it holds, a binary search per cell
  /// what the cell's time bound selects. Up to this many rows per cell the scan
  /// is the cheaper; denser cells -- many collisions in a block and a narrow
  /// time bound -- keep the searches. The rows and their order are the same
  /// either way, so this changes the cost and never the result.
  static constexpr Index kScanRowsPerCell = 8U;

  /// A query's cells, counted.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() CellSpan querySpan(LayerID target, const CandidateWindow& window, BlockRange blocks) const
  {
    const auto span = cellSpan(target, window, blocks);
    countGridWork(WorkCounter::kGridQueries);
    countGridWork(WorkCounter::kGridBlocks, static_cast<std::uint64_t>(span.blocks.empty() ? 0 : span.blocks.size()));
    return span;
  }

  /// \a body(phiBin) over the span's phi cells: the run up to the axis's end,
  /// then the one after the wrap.
  template <typename Body>
  GPUhdi() GPUAlwaysInline() static void forEachPhiBin(const CellSpan& span, const Body& body)
  {
    for (auto phiBin = span.firstPhi; phiBin < span.phiEnd; ++phiBin) {
      body(phiBin);
    }
    for (std::int32_t phiBin = 0; phiBin < span.wrapEnd; ++phiBin) {
      body(phiBin);
    }
  }

  template <typename Visitor>
  GPUhdi() GPUAlwaysInline() void visitRuns(LayerID target, std::int32_t block, std::int32_t phiBin, const CellSpan& span, const SecondRun* secondRuns,
                                            std::uint32_t secondRunCount, Visitor& visit, Lanes lanes) const
  {
    std::int32_t previousLast = -1;
    for (std::uint32_t i = 0; i < secondRunCount; ++i) {
      const auto& run = secondRuns[i];
      auto first = mGrid.secondBinOf(target, run.low);
      const auto last = mGrid.secondBinOf(target, run.high);
      first = first > previousLast ? first : previousLast + 1;
      if (first <= last) {
        visitSecondRun(target, block, phiBin, first, last, span, visit, lanes);
      }
      previousLast = last > previousLast ? last : previousLast;
    }
  }

  /// One `(block, phi)` pair's second cells `[firstSecond, lastSecond]`:
  /// consecutive buckets, because the second coordinate is the innermost key.
  template <typename Visitor, typename Share>
  GPUhdi() GPUAlwaysInline() void visitSecondRun(LayerID target, std::int32_t block, std::int32_t phiBin, std::int32_t firstSecond, std::int32_t lastSecond,
                                                 const CellSpan& span, Visitor& visit, Share lanes) const
  {
    visitBuckets(mGrid.bucket(target, block, phiBin, firstSecond), mGrid.bucket(target, block, phiBin, lastSecond) + 1U, span, visit, lanes);
  }

  /// The rows of the consecutive buckets `[firstBucket, endBucket)` whose box
  /// starts lie in the span's time bound and whose boxes end after its
  /// `afterPS`, in storage order.
  ///
  /// Only a cell's rows are ordered by box start -- a run of cells is
  /// contiguous in memory but not in time -- so the bound is either tested row
  /// by row over the whole range, or cell by cell: a binary search for the
  /// first row at or after `earliestPS`, and a walk from there to the first
  /// row after `latestPS`. Either way a row is kept exactly when it lies at or
  /// after the one and not after the other. The walk rather than a second
  /// search, because a dense cell's time bound holds a fraction of a row on
  /// average.
  ///
  /// With several \a lanes each takes every `count`-th row of the range, or
  /// of each cell's run after its search, from its own `lane` on; a lane
  /// past the time bound stops, as every row after it in the cell is too.
  template <typename Visitor>
  GPUhdi() void visitBuckets(std::uint32_t firstBucket, std::uint32_t endBucket, const CellSpan& span, Visitor& visit, Lanes lanes) const
  {
    const Index begin = mOffsets[firstBucket];
    const Index end = mOffsets[endBucket];
    countGridWork(WorkCounter::kGridRanges);
    countGridWork(WorkCounter::kGridCells, endBucket - firstBucket);
    countGridWork(WorkCounter::kGridRows, end - begin);
    if (end - begin <= kScanRowsPerCell * (endBucket - firstBucket)) {
      for (auto row = begin + lanes.lane; row < end; row += lanes.count) {
        const std::int64_t timePS = mMeasurements.timePS[row];
        if (timePS < span.earliestPS || span.latestPS < timePS) {
          continue;
        }
        countGridWork(WorkCounter::kGridInBound);
        visitIfCompatible(row, span, visit);
      }
      return;
    }
    countGridWork(WorkCounter::kGridSearchedCells, endBucket - firstBucket);
    for (auto bucket = firstBucket; bucket < endBucket; ++bucket) {
      const MeasurementRange range{.begin = mOffsets[bucket], .end = mOffsets[bucket + 1]};
      for (auto row = lowerBound(range, span.earliestPS) + lanes.lane; row < range.end; row += lanes.count) {
        if (span.latestPS < static_cast<std::int64_t>(mMeasurements.timePS[row])) {
          break;
        }
        countGridWork(WorkCounter::kGridInBound);
        visitIfCompatible(row, span, visit);
      }
    }
  }

  /// `visitBuckets` for a \a group in step: the group's lanes take the rows
  /// `Lanes{lane, kLanes}` would, `kLanes` at a time from the range's first
  /// row (scanned) or from each cell's first row in the bound (searched), and
  /// every lane calls \a visit(row, has) at every step. A searched cell ends
  /// at the first step whose first row -- the same row for every lane -- is
  /// past the time bound: the cell's rows are ordered by box start, so the
  /// lanes past the bound within a step are its last ones. A scanned range
  /// spans cells, whose times start again at each, so it goes to its end.
  template <typename Visitor, typename Group>
  GPUhdi() void visitBuckets(std::uint32_t firstBucket, std::uint32_t endBucket, const CellSpan& span, Visitor& visit, const Group& group) const
  {
    constexpr Index kLanes = Group::kLanes;
    const Index lane = group.lane();
    const Index begin = mOffsets[firstBucket];
    const Index end = mOffsets[endBucket];
    countGridWork(WorkCounter::kGridRanges);
    countGridWork(WorkCounter::kGridCells, endBucket - firstBucket);
    countGridWork(WorkCounter::kGridRows, end - begin);
    if (end - begin <= kScanRowsPerCell * (endBucket - firstBucket)) {
      for (Index first = begin; first < end; first += kLanes) {
        const Index row = first + lane;
        bool has = false;
        if (row < end) {
          const std::int64_t timePS = mMeasurements.timePS[row];
          has = span.earliestPS <= timePS && timePS <= span.latestPS && compatibleInBound(row, span);
        }
        visit(row, has);
      }
      return;
    }
    countGridWork(WorkCounter::kGridSearchedCells, endBucket - firstBucket);
    for (auto bucket = firstBucket; bucket < endBucket; ++bucket) {
      const MeasurementRange range{.begin = mOffsets[bucket], .end = mOffsets[bucket + 1]};
      for (Index first = lowerBound(range, span.earliestPS); first < range.end && !(span.latestPS < std::int64_t{mMeasurements.timePS[first]});
           first += kLanes) {
        const Index row = first + lane;
        const bool has = row < range.end && !(span.latestPS < std::int64_t{mMeasurements.timePS[row]}) && compatibleInBound(row, span);
        visit(row, has);
      }
    }
  }

  /// A row in the time bound, counted, and whether its box ends after the
  /// span's `afterPS` (`visitIfCompatible`'s test), counted where it does.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool compatibleInBound(Index row, const CellSpan& span) const
  {
    countGridWork(WorkCounter::kGridInBound);
    if (timeBoxAt(mMeasurements, row).upper() > span.afterPS) {
      countGridWork(WorkCounter::kGridVisits);
      return true;
    }
    return false;
  }

  /// A row in the bound whose box ends after the span's `afterPS`: with the
  /// bound's end, the boxes overlap.
  template <typename Visitor>
  GPUhdi() void visitIfCompatible(Index row, const CellSpan& span, Visitor& visit) const
  {
    if (timeBoxAt(mMeasurements, row).upper() > span.afterPS) {
      countGridWork(WorkCounter::kGridVisits);
      visit(row);
    }
  }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() Index lowerBound(MeasurementRange range, std::int64_t timePS) const
  {
    return o2::gpu::CAAlgo::partitionPoint(range.begin, range.size(), [&](Index row) {
      const std::int64_t start = mMeasurements.timePS[row];
      return start < timePS;
    });
  }

  MeasurementView mMeasurements;
  GridGeometry<Det> mGrid;
  const Index* mOffsets = nullptr;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_MEASUREMENTGRID_H_
