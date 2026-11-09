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

/// \file GridTestSupport.h
/// \brief Test grids of a known shape, and the toys' time boxes and frames.
///
/// The grid has no caps; tests that need a small grid set its counts
/// explicitly. `cappedGrid` gives
/// the layout one iteration with the given windows derives, each count held
/// to a ceiling.

#ifndef GTRK_TEST_GRIDTESTSUPPORT_H_
#define GTRK_TEST_GRIDTESTSUPPORT_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/TimeEstPS.h"
#include "GenericTracking/Tracking/GridSizing.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Tracking/SlicePlan.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <memory>
#include <vector>

namespace o2::gtrk::test
{

/// The toys' clock: 25 ns ticks, a tick's box widened by 8 ns at both ends,
/// so that two hits one tick apart are compatible
/// and two ticks apart are not.
inline constexpr double kToyTickNS = 25.;
inline constexpr double kToyMarginNS = 8.;
/// The width of every toy box not clamped at the timeframe's start [ps].
inline constexpr std::uint32_t kToyBoxWidthPS = 41'000U;
/// The time of a toy track whose hits share one tick [ns].
inline constexpr double kToyTrackTimeNS = 10.;

/// The box of the toy tick holding \a timeNS, clamped at the timeframe's start.
[[nodiscard]] inline TimeEstPS toyTimeBox(double timeNS)
{
  constexpr double kPSPerNS = 1000.;
  const double tick = std::floor(timeNS / kToyTickNS);
  const double lower = std::max((tick * kToyTickNS) - kToyMarginNS, 0.);
  const double upper = ((tick + 1.) * kToyTickNS) + kToyMarginNS;
  const auto start = std::llround(lower * kPSPerNS);
  return {.start = static_cast<std::uint32_t>(start), .width = static_cast<std::uint32_t>(std::llround(upper * kPSPerNS) - start)};
}

/// The toys' readout frames: one BC on every layer. Their boxes, widened
/// ticks, are not frames.
inline constexpr FrameTiming kToyFrames{.lengthBC = 1, .firstStartBC = 0};

/// \a params with the toys' frames on every layer.
[[nodiscard]] inline GridParams withToyFrames(GridParams params)
{
  std::ranges::fill(params.frames, kToyFrames);
  return params;
}

/// The windows of one iteration with \a windows, taking every link (chains
/// of two clusters, `takenLinks`) and extending: every layer a link touches
/// is sized.
template <DetectorLayout Det>
[[nodiscard]] auto oneIterationWindows(const SearchWindowParams& windows, const VertexConstraint& beamSpot)
{
  const IterationSettings settings{.windows = windows, .minClusters = 2};
  return gridLayerWindows<Det>(std::span<const IterationSettings>{&settings, 1U}, beamSpot);
}

/// The derived layout of one iteration with \a windows, each layer's counts
/// at most \a maxPhi and \a maxSecond.
template <DetectorLayout Det>
[[nodiscard]] GridGeometry<Det> cappedGrid(double durationNS, GridParams params, const SearchWindowParams& windows, const VertexConstraint& beamSpot, std::int32_t maxPhi, std::int32_t maxSecond)
{
  const auto layers = oneIterationWindows<Det>(windows, beamSpot);
  // The toys' times are in ns; the grid's in ps.
  const auto durationPS = static_cast<std::int64_t>(std::llround(durationNS * 1.e3));
  const GridGeometry<Det> derived{durationPS, params, layers};
  for (std::size_t layer = 0; layer < GridGeometry<Det>::kLayerCount; ++layer) {
    const auto& axes = derived.axes(static_cast<LayerID>(layer));
    params.phiBins[layer] = std::min(axes.nPhi, maxPhi);
    params.secondBins[layer] = std::min(axes.nSecond, maxSecond);
  }
  return GridGeometry<Det>{durationPS, params, layers};
}

/// Stable grid placement shared by frame fixtures and grid-only queries.
/// Callers supply the bucket of each input and fill its placed row.
template <typename BucketOf, typename Fill>
void buildGridRows(std::size_t count, std::span<std::uint32_t> offsets, const BucketOf& bucketOf, const Fill& fill)
{
  MeasurementGridBuilder builder{offsets};
  std::vector<std::uint32_t> buckets(count);
  for (std::size_t i = 0; i < count; ++i) {
    buckets[i] = bucketOf(i);
    builder.count(buckets[i]);
  }
  builder.prefixSum();
  for (std::size_t i = 0; i < count; ++i) {
    fill(i, builder.place(buckets[i]));
  }
  builder.finish();
}

struct FramePoint {
  LayerID layer = kInvalidLayerID;
  Scalar x = 0.F;
  Scalar y = 0.F;
  Scalar z = 0.F;
  Scalar r = 0.F;
  Scalar phi = 0.F;
};

struct PointVariances {
  Scalar z = 0.F;
  Scalar r = 0.F;
  Scalar rPhi = 0.F;
};

struct NoLocalMeasurement {
  void operator()(std::size_t /*input*/, std::uint32_t /*row*/) const {}
};

/// Build a frame in caller-owned columns and offsets. Logical input indices
/// become cluster indices; optional local measurements follow the same placement.
template <DetectorLayout Det, typename PointOf, typename BoxOf, typename VarianceOf, typename FillLocal = NoLocalMeasurement>
std::unique_ptr<TimeFrame<Det>> frameOf(const GridGeometry<Det>& geometry, std::size_t size, const PointOf& pointOf, const BoxOf& boxOf, const VarianceOf& varianceOf, MeasurementStore& columns, std::vector<std::uint32_t>& offsets, const FillLocal& fillLocal = {})
{
  resizeMeasurementColumns(columns, size);
  offsets.resize(static_cast<std::size_t>(geometry.bucketCount()) + 1U);
  buildGridRows(
    size, offsets,
    [&](std::size_t i) {
      const auto p = pointOf(i);
      return geometry.bucketAt(p.layer, boxOf(i).start, p.phi, secondCoordinate<Det>(p.layer, p.r, p.z));
    },
    [&](std::size_t i, std::uint32_t row) {
      const auto p = pointOf(i);
      const auto box = boxOf(i);
      const auto variance = varianceOf(i);
      columns.timePS[row] = box.start;
      columns.widthPS[row] = box.width;
      columns.clusterIndex[row] = static_cast<std::uint32_t>(i);
      columns.x[row] = p.x;
      columns.y[row] = p.y;
      columns.z[row] = p.z;
      columns.r[row] = p.r;
      columns.phi[row] = p.phi;
      columns.varZ[row] = variance.z;
      columns.varR[row] = variance.r;
      columns.varRPhi[row] = variance.rPhi;
      fillLocal(i, row);
    });
  return std::make_unique<TimeFrame<Det>>(measurementView(columns), geometry, std::span<const std::uint32_t>{offsets});
}

/// Iteration \a iteration of \a tracker over \a slice alone
/// (`Tracker::runSlice`), with the frame's bounds and the windows of
/// \a vertices as `planTimeframe` plans them, the slice planned as it plans
/// its slices.
template <typename Tracker, DetectorLayout Det>
[[nodiscard]] auto runSlice(Tracker& tracker, std::size_t iteration, const TimeFrame<Det>& frame, const typename TimeFrame<Det>::TimeSlice& slice, std::span<const VertexSeed> vertices)
{
  const auto plan = tracker.planTimeframe(frame, vertices, frame.grid().geometry().blockCount(), 0);
  const auto slicePlan = planSlice(frame, slice, sliceReach(frameBounds(frame), frame.grid().geometry()).gatePS());
  return tracker.runSlice(iteration, frame, plan, slicePlan, vertices, plan.vertexWindows[iteration]);
}

} // namespace o2::gtrk::test

#endif // GTRK_TEST_GRIDTESTSUPPORT_H_
