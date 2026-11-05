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

/// \file SlicePlan.h
/// \brief A timeframe's slices, planned once on the host from the frame and
///        the iterations, and run by the host and the device trackers alike.
///
/// A slice's gate window and the halo width depend on the measurements'
/// time variances, not only on the grid offsets; the tracklet stage keeps
/// them as values. The plan computes them once, with the source maps, the
/// vertex windows and the extension's per-layer bound, and records what it
/// was planned from, so that a tracker can refuse a plan of another frame or
/// configuration.

#ifndef GTRK_TRACKING_SLICEPLAN_H_
#define GTRK_TRACKING_SLICEPLAN_H_

#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LayerRowMap.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TimeFrameScratch.h"
#include "GenericTracking/Tracking/TrackletStage.h"
#include "GenericTracking/Tracking/TripletCut.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <stdexcept>
#include <vector>

namespace o2::gtrk
{

/// Every layer's largest values over a range of rows: the time box's width
/// the tracklet stage's and the extension's windows and the slices' reach are
/// bounded with, and the z, r and r-phi variances the triplet stage's
/// dip-angle window is; and whether every row's box is exactly one of the
/// layer's readout frames (`FrameTiming::isFrame`), which the slices' reach
/// reads (`sliceReach`).
template <DetectorLayout Det>
struct HaloBounds {
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);
  std::array<std::uint32_t, kLayerCount> widthPS{};
  std::array<TripletTargetBounds, kLayerCount> target{};
  std::array<bool, kLayerCount> rowsAreFrames{};
};

template <DetectorLayout Det>
[[nodiscard]] HaloBounds<Det> boundsOver(const TimeFrame<Det>& frame, std::span<const MeasurementRange, HaloBounds<Det>::kLayerCount> rows)
{
  HaloBounds<Det> result;
  const auto columns = frame.measurements();
  const auto& geometry = frame.grid().geometry();
  for (std::size_t layer = 0; layer < HaloBounds<Det>::kLayerCount; ++layer) {
    auto& target = result.target[layer];
    const auto& frames = geometry.frames(static_cast<LayerID>(layer));
    result.rowsAreFrames[layer] = true;
    for (auto row = rows[layer].begin; row < rows[layer].end; ++row) {
      result.widthPS[layer] = std::max(result.widthPS[layer], columns.widthPS[row]);
      result.rowsAreFrames[layer] = result.rowsAreFrames[layer] && frames.isFrame({.start = columns.timePS[row], .width = columns.widthPS[row]});
      target.varZ = std::max(target.varZ, columns.varZ[row]);
      target.varR = std::max(target.varR, columns.varR[row]);
      target.varRPhi = std::max(target.varRPhi, columns.varRPhi[row]);
    }
  }
  return result;
}

/// The bounds over the whole frame: what every slice's stages take, so that
/// no window depends on the slicing. The dip-angle window is a cut at its
/// margin, not only a pre-filter, so a bound taken over a slice's halo would
/// make the candidates depend on which rows the halo held.
template <DetectorLayout Det>
[[nodiscard]] HaloBounds<Det> frameBounds(const TimeFrame<Det>& frame)
{
  std::array<MeasurementRange, HaloBounds<Det>::kLayerCount> layers{};
  for (std::size_t layer = 0; layer < layers.size(); ++layer) {
    layers[layer] = frame.grid().layerRows(static_cast<LayerID>(layer));
  }
  return boundsOver(frame, std::span<const MeasurementRange, HaloBounds<Det>::kLayerCount>{layers});
}

/// How far in time a slice's stages reach past its core, in ps, from each
/// link's pair reaches: how much later and
/// how much earlier the box of a pair's target can start than its source's.
///
/// Two overlapping boxes `a` (the source) and `b` (the target) have
/// `b.start - a.start` in `(-width_b, width_a)`: a link reaches its source
/// layer's widest box forward and its target layer's backward. When the two
/// layers have readout frames, the same ones, and the rows of both are exactly
/// those frames, two boxes are equal or disjoint, so a pair lies in one frame
/// and the link reaches 0 both ways. Pairwise overlap is not transitive once
/// boxes differ, so a chain drifts by its links' reaches: `R+` and `R-` are
/// the largest sums of the forward and of the backward reaches along a path
/// through the links.
struct SliceReach {
  std::int64_t pairForwardPS = 0;  ///< the largest link's forward reach
  std::int64_t pairBackwardPS = 0; ///< the largest link's backward reach
  std::int64_t forwardPS = 0;      ///< `R+`
  std::int64_t backwardPS = 0;     ///< `R-`

  /// How far from the core every stage works: `R+ + R-` on both sides. A road
  /// the core owns starts in it and its measurements start within `R-` before
  /// and `R+` after its first; the levels on it depend on the chains into its
  /// measurements, which start within `R+` before and `R-` after them. Within
  /// the gate everything is built, so the output does not depend on the
  /// slicing. With the same widths on every link, `2R`.
  [[nodiscard]] std::int64_t gatePS() const { return forwardPS + backwardPS; }

  /// The time halo a slice needs: the gate, plus one link's reach for the
  /// targets of the tracklets starting at its edges (the backward one before
  /// it, the forward one after; the halo is whole blocks, so the larger).
  [[nodiscard]] std::int64_t haloPS() const { return gatePS() + std::max(pairForwardPS, pairBackwardPS); }
};

/// The reach of a frame whose rows have \a bounds, on \a geometry's frames.
template <DetectorLayout Det>
[[nodiscard]] SliceReach sliceReach(const HaloBounds<Det>& bounds, const GridGeometry<Det>& geometry)
{
  const auto sameFrames = [&](const LayerLink& link) {
    const auto& frames = geometry.frames(link.from);
    return frames == geometry.frames(link.to) && bounds.rowsAreFrames[link.from] && bounds.rowsAreFrames[link.to];
  };
  SliceReach result;
  std::array<std::int64_t, HaloBounds<Det>::kLayerCount> forwardTo{};
  std::array<std::int64_t, HaloBounds<Det>::kLayerCount> backwardTo{};
  for (const auto layer : layerOrder<Det>) {
    result.forwardPS = std::max(result.forwardPS, forwardTo[layer]);
    result.backwardPS = std::max(result.backwardPS, backwardTo[layer]);
    for (auto i = outLinks<Det>.begin(layer); i < outLinks<Det>.end(layer); ++i) {
      const auto& link = Det::layerLinks[outLinks<Det>.links[i]];
      const bool zero = sameFrames(link);
      const std::int64_t forward = zero ? 0 : bounds.widthPS[link.from];
      const std::int64_t backward = zero ? 0 : bounds.widthPS[link.to];
      result.pairForwardPS = std::max(result.pairForwardPS, forward);
      result.pairBackwardPS = std::max(result.pairBackwardPS, backward);
      forwardTo[link.to] = std::max(forwardTo[link.to], forwardTo[layer] + forward);
      backwardTo[link.to] = std::max(backwardTo[link.to], backwardTo[layer] + backward);
    }
  }
  return result;
}

/// One slice of a timeframe: its blocks and rows, its tracklet sources, and
/// the gate's window of box starts around its core.
template <DetectorLayout Det>
struct SlicePlan {
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);
  using Slice = typename TimeFrame<Det>::TimeSlice;

  Slice slice;
  LayerRowMap<Det> sources;
  std::int64_t gateEarliestPS = 0;
  std::int64_t gateLatestPS = 0;

  /// The tracklet stage's inputs, with the frame's widest candidate boxes per
  /// layer \a maxCandidateWidthPS (`TimeframePlan::layerMaxWidthPS`) and the
  /// gate's core ranges pointing into this plan, which must outlive them.
  [[nodiscard]] TrackletInputs<Det> inputs(const IterationParameters<Det>& parameters, const std::array<std::uint32_t, kLayerCount>& maxCandidateWidthPS) const
  {
    return {.parameters = &parameters,
            .sources = &sources,
            .targetBlocks = slice.haloBlocks,
            .maxCandidateWidthPS = maxCandidateWidthPS,
            .gate = {.core = slice.core.data(), .earliestPS = gateEarliestPS, .latestPS = gateLatestPS},
            .halo = slice.halo};
  }
};

/// The plan of \a slice, its gate reaching \a gatePS past the core
/// (`SliceReach::gatePS` of the frame's bounds).
template <DetectorLayout Det>
[[nodiscard]] SlicePlan<Det> planSlice(const TimeFrame<Det>& frame, const typename TimeFrame<Det>::TimeSlice& slice, std::int64_t gatePS)
{
  SlicePlan<Det> plan{.slice = slice, .sources = TimeFrameScratch<Det>::makeSources(slice.halo)};
  // The core's box starts, widened by the gate: a halo source beyond it lies
  // on no chain the core's output depends on. The reach is the frame's rows',
  // as the halo's; the core's edges are rounded outward to picoseconds.
  const auto& geometry = frame.grid().geometry();
  const std::int64_t parts = geometry.partsPerPS();
  const auto floorDiv = [](std::int64_t a, std::int64_t b) { return (a / b) - ((a % b) < 0 ? 1 : 0); };
  plan.gateEarliestPS = floorDiv(geometry.blockStart(slice.coreBlocks.begin), parts) - gatePS;
  plan.gateLatestPS = -floorDiv(-geometry.blockStart(slice.coreBlocks.end), parts) + gatePS;
  return plan;
}

/// A timeframe's plan for every iteration of a run -- its slices, the same
/// for every iteration, and each iteration's vertex windows -- with what it
/// was planned from.
template <DetectorLayout Det>
struct TimeframePlan {
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);

  // What it was planned from.
  std::uint32_t rows = 0U;
  std::int32_t blockCount = 0;
  std::int64_t blockWidth = 0; ///< [parts of a ps]
  std::vector<IterationParameters<Det>> parameters;

  VertexSeed beamLine;
  std::vector<VertexSeed> vertices; ///< the timeframe's vertex list, sorted by time
  /// The widest candidate box per layer, over the whole frame (`frameBounds`).
  std::array<std::uint32_t, kLayerCount> layerMaxWidthPS{};
  /// The triplet stage's target bounds, over the whole frame (`frameBounds`).
  std::array<TripletTargetBounds, kLayerCount> layerTargetBounds{};
  std::int64_t haloPS = 0;
  std::vector<SlicePlan<Det>> slices;                     ///< in order
  std::vector<VertexWindowParameters<Det>> vertexWindows; ///< per iteration

  /// Throw std::invalid_argument unless this plan was planned from \a frame
  /// -- its rows and grid -- and from \a configured. The parameters compare
  /// as bytes: every copy of them comes from the same objects by trivial
  /// copies, padding included.
  void check(const TimeFrame<Det>& frame, std::span<const IterationParameters<Det>> configured) const
  {
    const auto& geometry = frame.grid().geometry();
    if (rows != frame.size() || blockCount != geometry.blockCount() || blockWidth != geometry.blockWidth()) {
      throw std::invalid_argument("a timeframe plan of another frame");
    }
    if (configured.size() != parameters.size() || vertexWindows.size() != parameters.size()) {
      throw std::invalid_argument("a timeframe plan of another number of iterations");
    }
    for (std::size_t i = 0; i < parameters.size(); ++i) {
      // Padding included, deliberately: both sides are copies of the same
      // objects, made by std::copy or assign, which move their bytes as they are.
      // NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
      if (std::memcmp(&configured[i], &parameters[i], sizeof(IterationParameters<Det>)) != 0) {
        throw std::invalid_argument("a timeframe plan of other iteration parameters");
      }
    }
  }
};

/// Plan every iteration of \a iterations over \a frame in slices of
/// \a sliceWidthBlocks grid blocks (`GridGeometry::blocksPerSlice` of the
/// configured width) with at least \a configuredHaloPS of halo, more when
/// the frame's reach needs it -- the vertexing pass too, each
/// slice finding its core's frames' vertices (`OwnedFrames`). \a
/// windowSettings are the settings each iteration's parameters were derived
/// from; \a beamLine the seed a beam-line iteration reads; \a vertices the
/// timeframe's vertex list, sorted by time.
template <DetectorLayout Det>
[[nodiscard]] TimeframePlan<Det> planTimeframe(const TimeFrame<Det>& frame, std::span<const IterationParameters<Det>> iterations,
                                               std::span<const SearchWindowParams> windowSettings, const VertexSeed& beamLine,
                                               std::span<const VertexSeed> vertices, std::int32_t sliceWidthBlocks, std::int64_t configuredHaloPS)
{
  if (windowSettings.size() != iterations.size()) {
    throw std::invalid_argument("planTimeframe needs one window setting per iteration");
  }
  TimeframePlan<Det> plan;
  plan.rows = frame.size();
  plan.blockCount = frame.grid().geometry().blockCount();
  plan.blockWidth = frame.grid().geometry().blockWidth();
  plan.parameters.assign(iterations.begin(), iterations.end());
  plan.beamLine = beamLine;
  plan.vertices.assign(vertices.begin(), vertices.end());
  const auto bounds = frameBounds(frame);
  plan.layerMaxWidthPS = bounds.widthPS;
  plan.layerTargetBounds = bounds.target;
  const auto reach = sliceReach(bounds, frame.grid().geometry());
  plan.haloPS = std::max(configuredHaloPS, reach.haloPS());
  const auto width = std::max<std::int32_t>(1, sliceWidthBlocks);
  for (std::int32_t begin = 0; begin < plan.blockCount; begin += width) {
    const typename TimeFrame<Det>::TimeSlice slice{frame, {.begin = begin, .end = std::min(plan.blockCount, begin + width)}, plan.haloPS};
    plan.slices.push_back(planSlice(frame, slice, reach.gatePS()));
  }
  plan.vertexWindows.reserve(iterations.size());
  for (std::size_t iteration = 0; iteration < iterations.size(); ++iteration) {
    plan.vertexWindows.push_back(vertexWindowParameters<Det>(iterations[iteration], windowSettings[iteration], vertices));
  }
  return plan;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SLICEPLAN_H_
