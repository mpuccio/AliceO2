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

/// \file TrackletStage.h
/// \brief The tracklet stage's per-source function.

#ifndef GTRK_TRACKING_TRACKLETSTAGE_H_
#define GTRK_TRACKING_TRACKLETSTAGE_H_

#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/DisplacedPair.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LayerRowMap.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/ShapeCut.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/Trig.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonAlgorithm.h"
#include "GPUCommonMath.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace o2::gtrk
{

/// Which of a slice's source rows the tracklet stage works on: every core
/// row, and a halo row only if its box starts within the gate
/// (`SliceReach::gatePS`) of the core's: everything a chain the output
/// depends on can reach. The halo itself is whole time blocks,
/// far wider. With no core ranges every source is worked.
struct SourceGate {
  const MeasurementRange* core = nullptr; ///< per layer
  std::int64_t earliestPS = 0;
  std::int64_t latestPS = 0;
};

/// What one slice's tracklet stage is built from besides the timeframe, the
/// vertices and the used mask. The tracker derives it once per slice, so that
/// a device run of a slice uses exactly what the host run did.
/// Host memory throughout; the gate's core
/// ranges point into the slice, which must outlive it.
template <DetectorLayout Det>
struct TrackletInputs {
  const IterationParameters<Det>* parameters = nullptr;
  const LayerRowMap<Det>* sources = nullptr;
  BlockRange targetBlocks{};
  std::array<std::uint32_t, LayerRowMap<Det>::kLayerCount> maxCandidateWidthPS{};
  SourceGate gate{};
  /// The slice's halo rows on every layer, source or not: what the later
  /// stages' variance bounds cover.
  std::array<MeasurementRange, LayerRowMap<Det>::kLayerCount> halo{};
};

/// Where a tracklet stage's inputs are for a run somewhere else: the same
/// objects copied to a device, each address in that device's memory.
/// `TrackletStage::bind` swaps them in.
template <DetectorLayout Det>
struct TrackletBindings {
  const IterationParameters<Det>* parameters = nullptr;
  const LayerRowMap<Det>* sources = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const std::uint8_t* used = nullptr;
  const VertexWindowParameters<Det>* vertexWindows = nullptr;
  const VertexSeed* vertices = nullptr;       ///< a copy of the vertices the stage reads: the beam line or the list
  const MeasurementRange* gateCore = nullptr; ///< per layer, if the stage is gated
};

/// A displaced iteration's pair test:
/// `displacedPairWithin` with the two points' errors, \a scattering the
/// source's (`sourceScattering`), \a between the passives' between the two
/// (`betweenScattering`), the beam line's width (\a vertex, the beam line),
/// the pair cut's curvature bound and sigma, and the iteration's impact bound
/// and z bound about the beam spot's centre. The tracklet stage's, the stage
/// funnel's and the sink's brute force's, one test.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() displaced::DisplacedPairResult displacedPairOf(const IterationParameters<Det>& parameters, const VertexSeed& vertex,
                                                                                        const PairPoint& source, Scalar scattering, const PairPoint& target,
                                                                                        Scalar between)
{
  const displaced::PairErrors<Scalar> errors{.sourceVarR = source.varR,
                                             .sourceVarRPhi = source.varRPhi,
                                             .sourceVarZ = source.varZ,
                                             .targetVarR = target.varR,
                                             .targetVarRPhi = target.varRPhi,
                                             .targetVarZ = target.varZ,
                                             .sourceScatter = scattering,
                                             .betweenScatter = between,
                                             .beamVarTransverse = vertex.varTransverse};
  const Scalar centre = parameters.vertexZCentre;
  return displaced::displacedPairWithin<Scalar>(source.x, source.y, source.z, target.x, target.y, target.z, errors, parameters.pairCut.nSigma,
                                                parameters.pairCut.maxCurvature, parameters.impactMax, centre - parameters.impactZMax,
                                                centre + parameters.impactZMax);
}

/// Generate every accepted outgoing tracklet of one source measurement.
///
/// The call operator has the common StageOutput/CpuStage generator shape.
/// Sources, outgoing links and grid candidates are visited in their canonical
/// order; the first compatible vertex whose pair cut passes accepts a target.
///
/// With \a BeamLine the stage serves only a beam-line iteration: one vertex,
/// no per-vertex intervals and no fallback past the cache, so neither is
/// compiled. A beam-line iteration makes the same tracklets with either; the
/// trackers take this one for it.
/// With \a Displaced, on the beam line only, the stage serves a displaced
/// iteration (`IterationParameters::impactMax`): its pair test is `displacedPairWithin`, a circle above pT_min
/// passing within the impact bound of the beam line at a z within the z
/// bound of the beam spot's centre, in place of the beam point's.
///
/// With \a Shape, the iteration's `shapeTest`: a candidate's line must cross
/// both its measurements at a local slope in their slope intervals
/// (ShapeCut.h), tested first.
template <DetectorLayout Det, std::uint32_t VertexCapacity = 8U, bool BeamLine = false, bool Displaced = false, bool Shape = false>
class TrackletStage
{
 public:
  using Index = typename LayerRowMap<Det>::Index;
  static constexpr std::size_t kLayerCount = LayerRowMap<Det>::kLayerCount;
  static_assert(!BeamLine || VertexCapacity == 1U, "the beam line is one vertex");
  static_assert(!Displaced || BeamLine, "a displaced iteration pairs through the beam line");
  // How many compatible vertices get their own interval. A source with more
  // takes the beam-spot window with every one of them tested; both windows
  // hold every vertex `vertexUsable` admits, so the tracklets do not depend on
  // it.
  static constexpr std::uint32_t kVertexCapacity = VertexCapacity;
  static_assert(kVertexCapacity > 0U);
  static_assert(nLinks<Det> <= std::size_t{std::numeric_limits<std::uint16_t>::max()} + 1U);

  TrackletStage(const IterationParameters<Det>& parameters, const LayerRowMap<Det>& sources, const TimeFrame<Det>& frame, BlockRange targetBlocks,
                const VertexSeed& beamLine, std::span<const VertexSeed> vertices, const std::array<std::uint32_t, kLayerCount>& maxCandidateWidthPS,
                const VertexWindowParameters<Det>& vertexWindows, const std::uint8_t* used, const SourceGate& gate = {})
    : mParameters(&parameters),
      mSources(&sources),
      mFrame(&frame),
      mTargetBlocks(targetBlocks),
      mGate(gate),
      mUsed(used),
      mVertexWindows(&vertexWindows),
      mMaxCandidateWidthPS(maxCandidateWidthPS)
  {
    if (BeamLine && parameters.vertexSource != VertexSource::kBeamLine) {
      throw std::invalid_argument("the beam-line tracklet stage serves only a beam-line iteration");
    }
    if (Displaced != (parameters.impactMax > 0.F)) {
      throw std::invalid_argument("the displaced tracklet stage serves exactly the displaced iterations");
    }
    if (Shape != parameters.shapeTest) {
      throw std::invalid_argument("the shape-testing tracklet stage serves exactly the iterations with shapeTest");
    }
    if (parameters.vertexSource == VertexSource::kBeamLine) {
      mVertices = &beamLine;
      mVertexCount = 1U;
    } else {
      mVertices = vertices.data();
      mVertexCount = static_cast<std::uint32_t>(vertices.size());
    }
    for (std::uint32_t i = 0; i < mVertexCount; ++i) {
      mMaxVertexWidthPS = std::max(mMaxVertexWidthPS, mVertices[i].widthPS);
    }
  }

  /// This stage reading its inputs at \a bindings instead: every pointer
  /// swapped, every value the constructor derived on the host -- the vertex
  /// choice and count, the variance bounds, the gate's window -- kept, so
  /// that the bound stage runs on exactly the numbers this one does. Host
  /// code; the result is what a kernel receives.
  [[nodiscard]] TrackletStage bind(const TrackletBindings<Det>& bindings) const
  {
    const bool complete = bindings.parameters != nullptr && bindings.sources != nullptr && bindings.frame != nullptr && bindings.used != nullptr &&
                          bindings.vertexWindows != nullptr && (bindings.vertices != nullptr || mVertexCount == 0U) &&
                          (bindings.gateCore != nullptr || mGate.core == nullptr);
    if (!complete) {
      throw std::invalid_argument("a tracklet stage binding is missing an input the stage reads");
    }
    TrackletStage bound = *this;
    bound.mParameters = bindings.parameters;
    bound.mSources = bindings.sources;
    bound.mFrame = bindings.frame;
    bound.mUsed = bindings.used;
    bound.mVertexWindows = bindings.vertexWindows;
    bound.mVertices = bindings.vertices;
    bound.mGate.core = mGate.core == nullptr ? nullptr : bindings.gateCore;
    return bound;
  }

  using Between = BetweenSource<Det>;
  /// The measurement columns, which a device's candidate tests read
  /// (`testCandidate`, `record`).
  [[nodiscard]] GPUhdi() GPUAlwaysInline() MeasurementView columns() const { return mFrame->measurements(); }

  /// The vertices a source can pair through, cached with their scattering at
  /// the source.
  struct CompatibleVertex {
    std::uint32_t vertexIndex;
    Scalar scattering;
  };

  struct VertexRange {
    std::uint32_t begin = 0U;
    std::uint32_t end = 0U;
  };

  /// What one source brings to each of its queries (`prepare`), besides its
  /// cached vertices: an array apart, so that a compiler can keep the
  /// context's scalars in registers while the array is indexed at run time.
  struct SourceContext {
    typename LayerRowMap<Det>::Row at{};
    TimeEstPS box{};
    PairPoint point{}; ///< its position as the tests read it
    Scalar phi = 0.F;
    Scalar second = 0.F;
    VertexRange vertexRange{};
    VertexRange uncachedRange{};
    std::uint32_t cachedCount = 0U;
  };

  /// One query of a source: a link and the window it searches (`query`).
  struct QueryContext {
    std::size_t linkIndex = 0U;
    LayerID targetLayer = kInvalidLayerID;
    std::uint64_t holes = 0U;
    CandidateWindow window{};
    Scalar directLow = -std::numeric_limits<Scalar>::infinity();
    Scalar directHigh = std::numeric_limits<Scalar>::infinity();
  };

  /// With \a lanes, this lane's share of every query's candidates (`Lanes`):
  /// lanes walking one source together emit its tracklets between them, in
  /// another order than one thread would. The order does not matter: each
  /// source's tracklets are sorted by (link, lambda, target) afterwards, a
  /// key no two of them share.
  ///
  /// The source's work: `prepare` once, a query per out-link, and per
  /// candidate `cheap`, then `expensive`; a device runs the candidates apart
  /// (`forEachCandidate`, `testCandidate`, `record`).
  template <typename Emitter>
  GPUhdi() void operator()(Index sourceElement, Emitter& emitter, Lanes lanes = {}) const
  {
    const auto test = [&](const SourceContext& source, const CompatibleVertex* cachedVertices, const BetweenSource<Det>& between, const QueryContext& query,
                          const SecondRange* intervals, const MeasurementView& columns,
                          Index targetRow) { expensive(source, cachedVertices, between, query, intervals, columns, targetRow, emitter); };
    forEachCandidate(sourceElement, test, lanes);
  }

  /// The candidates of \a sourceElement that pass `cheap`, counted without
  /// the pair tests: what a flat candidate list holds for the source.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t candidateCount(Index sourceElement) const
  {
    std::uint32_t count = 0U;
    const auto one = [&](const SourceContext&, const CompatibleVertex*, const BetweenSource<Det>&, const QueryContext&, const SecondRange*,
                         const MeasurementView&, Index) { ++count; };
    forEachCandidate(sourceElement, one);
    return count;
  }

  /// \a body(source, cachedVertices, between, query, intervals, columns,
  /// targetRow) for every candidate of \a sourceElement that passes `cheap`,
  /// in the serial order -- the out-links in order, each query's rows in the
  /// grid walk's order: the one traversal the stage's forms share.
  /// With \a lanes, this lane's share of it.
  template <typename Body>
  GPUhdi() void forEachCandidate(Index sourceElement, Body& body, Lanes lanes = {}) const
  {
    SourceContext source;
    CompatibleVertex cachedVertices[kVertexCapacity];
    if (!prepare(sourceElement, source, cachedVertices)) {
      return;
    }
    // Apart from the context, as the vertex cache (its passives are indexed at run time).
    BetweenSource<Det> between;
    if constexpr (BeamLine) {
      between = betweenSource<Det>(mParameters->scattering, mVertices[cachedVertices[0].vertexIndex], source.point);
    }
    const MeasurementView columns = mFrame->measurements();
    forEachQuery(source, [&](const QueryContext& query) {
      // One candidate of the query, against the vertices' second-coordinate
      // \a intervals, or nullptr when the query took the whole window.
      const auto visitTarget = [&](const SecondRange* intervals) {
        return [&, intervals](Index targetRow) {
          if (cheap(targetRow)) {
            body(source, cachedVertices, between, query, intervals, columns, targetRow);
          }
        };
      };
      if constexpr (BeamLine) {
        mFrame->grid().visitCandidates(query.targetLayer, query.window, mTargetBlocks, visitTarget(nullptr), lanes);
      } else {
        const bool fallback = source.uncachedRange.begin != source.vertexRange.end;
        if (mParameters->vertexSource == VertexSource::kBeamLine || fallback) {
          mFrame->grid().visitCandidates(query.targetLayer, query.window, mTargetBlocks, visitTarget(nullptr), lanes);
          return false;
        }

        SecondRange intervals[kVertexCapacity];
        SecondRun runs[kVertexCapacity];
        const std::uint32_t runCount = vertexRuns(source, cachedVertices, query, intervals, runs);
        mFrame->grid().visitCandidates(query.targetLayer, query.window, runs, runCount, mTargetBlocks, visitTarget(intervals), lanes);
      }
      return false;
    });
  }

  /// Whether \a sourceElement has a query: `prepare` and some out-link's
  /// `prepareQuery`, as `operator()` decides them, without reading a target
  /// row. A source without one emits nothing, so a device launches the
  /// stage over the sources with one only.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool hasQuery(Index sourceElement) const
  {
    SourceContext source;
    CompatibleVertex cachedVertices[kVertexCapacity];
    return prepare(sourceElement, source, cachedVertices) && forEachQuery(source, [](const QueryContext& /*query*/) { return true; });
  }

 private:
  /// \a body(query) for each out-link of \a source that makes a query
  /// (`prepareQuery`), in the links' order, until it returns true; whether
  /// it did.
  template <typename Body>
  GPUhdi() GPUAlwaysInline() bool forEachQuery(const SourceContext& source, const Body& body) const
  {
    for (auto at = mParameters->outLinkOffsets[source.at.layer]; at < mParameters->outLinkOffsets[static_cast<std::size_t>(source.at.layer) + 1U]; ++at) {
      QueryContext query;
      if (prepareQuery(source, mParameters->outLinkIndices[at], query) && body(query)) {
        return true;
      }
    }
    return false;
  }

  /// The source of \a sourceElement in \a source and its vertices in
  /// \a cachedVertices (`kVertexCapacity` of them): false when it has no
  /// query -- held exclusively by an earlier track, outside the slice's gate,
  /// or with no vertex to pair through.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool prepare(Index sourceElement, SourceContext& source, CompatibleVertex* cachedVertices) const
  {
    source.at = mSources->row(sourceElement);
    // A row an earlier iteration's track holds exclusively is on no road of
    // this one; a shared one may still be.
    if (mUsed[source.at.row] == kExclusive) {
      return false;
    }
    countWork(WorkCounter::kTrackletSources);
    const auto measurement = mFrame->measurement(source.at.row);
    source.box = timeBox(measurement);
    if (mGate.core != nullptr) {
      const auto& core = mGate.core[source.at.layer];
      const bool inCore = source.at.row >= core.begin && source.at.row < core.end;
      const std::int64_t start = source.box.lower();
      if (!inCore && (start < mGate.earliestPS || start > mGate.latestPS)) {
        return false;
      }
    }
    source.point = pairPoint(measurement);
    source.phi = measurement.phi;
    source.vertexRange = compatibleVertexRange(source.box);
    source.cachedCount = 0U;
    source.uncachedRange = {.begin = source.vertexRange.end, .end = source.vertexRange.end};
    for (std::uint32_t i = source.vertexRange.begin; i < source.vertexRange.end; ++i) {
      const auto& vertex = mVertices[i];
      if (!vertexUsable(*mParameters, source.box, vertex)) {
        continue;
      }
      if (source.cachedCount == kVertexCapacity) {
        source.uncachedRange.begin = i;
        source.uncachedRange.end = source.vertexRange.end;
        break;
      }
      cachedVertices[source.cachedCount] = {
        .vertexIndex = i, .scattering = sourceScattering<Det>(mParameters->scattering, mParameters->pairCut, vertex, source.at.layer, source.point)};
      ++source.cachedCount;
    }
    if (source.cachedCount == 0U) {
      return false;
    }
    source.second = secondCoordinate<Det>(source.at.layer, source.point.r, source.point.z);
    return true;
  }

  /// The query of \a source on link \a linkIndex in \a query: false when it
  /// makes none -- the source outside the link's reach, or no direct pair.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool prepareQuery(const SourceContext& source, std::size_t linkIndex, QueryContext& query) const
  {
    query.linkIndex = linkIndex;
    query.targetLayer = mParameters->linkTargets[linkIndex];
    if (source.second < mParameters->sourceMin[linkIndex] || source.second > mParameters->sourceMax[linkIndex]) {
      return false;
    }
    query.holes = mParameters->holeLayers[linkIndex];
    const auto& at = source.point;
    query.window = candidateWindow(mParameters->windows[linkIndex], source.box, source.phi, at.r, at.z, mMaxCandidateWidthPS[query.targetLayer]);
    // In an iteration without skips, a source queries only the part of its
    // window a pair of it can reach without crossing another layer, and
    // nothing when there is none (LinkHoles.h); every pair is still tested.
    query.directLow = -std::numeric_limits<Scalar>::infinity();
    query.directHigh = std::numeric_limits<Scalar>::infinity();
    if (query.holes != 0U) {
      const auto& target = mParameters->layerSpans[query.targetLayer];
      const auto direct =
        directPairRange(query.holes, mParameters->barrelLayers, mParameters->layerSpans, query.targetLayer, at.r, at.z, target.low, target.high);
      const Scalar low = o2::gpu::CAMath::Max(direct.low, query.window.secondCentre - query.window.secondHalfWidth);
      const Scalar high = o2::gpu::CAMath::Min(direct.high, query.window.secondCentre + query.window.secondHalfWidth);
      if (!(low <= high)) {
        countWork(WorkCounter::kTrackletHoleQueries);
        return false;
      }
      query.directLow = direct.low;
      query.directHigh = direct.high;
      query.window.secondCentre = 0.5F * (low + high);
      query.window.secondHalfWidth = 0.5F * (high - low);
    }
    countWork(WorkCounter::kTrackletQueries);
    countLinkWork(linkIndex, LinkCounter::kQueries);
    return true;
  }

  /// A candidate's first test: false when an earlier track holds the row
  /// exclusively.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool cheap(Index targetRow) const
  {
    if (mUsed[targetRow] == kExclusive) {
      countWork(WorkCounter::kTrackletUsed);
      return false;
    }
    return true;
  }

  /// The rest of a candidate's tests -- the pair cut through every vertex the
  /// source can pair with, then the hole -- and its tracklet if they pass.
  /// \a between: the source's part of `betweenScattering` on the beam line
  /// (`betweenSource`), unread on a vertex list. \a intervals: the vertices'
  /// second-coordinate intervals of the query, or nullptr when it took the
  /// whole window.
  ///
  /// The hole test last: both tests have to pass and neither has a side
  /// effect, and the pair test refuses most rows.
  template <typename Columns, typename Emitter>
  GPUhdi() void expensive(const SourceContext& source, const CompatibleVertex* cachedVertices, const BetweenSource<Det>& between, const QueryContext& query,
                          const SecondRange* intervals, const Columns& columns, Index targetRow, Emitter& emitter) const
  {
    const auto target = pairPointAt(columns, targetRow);
    if (!testCandidate(source, cachedVertices, between, query, intervals, columns, targetRow, target)) {
      return;
    }
    emit(source.at.row, targetRow, query.linkIndex, source.point, target, emitter);
  }

 public:
  /// `expensive`'s tests on \a target, the candidate's point: true where the
  /// candidate is a tracklet (the test apart from the record).
  template <typename Columns>
  [[nodiscard]] GPUhdi() bool testCandidate(const SourceContext& source, const CompatibleVertex* cachedVertices, const BetweenSource<Det>& between,
                                            const QueryContext& query, const SecondRange* intervals, const Columns& columns, Index targetRow,
                                            const PairPoint& target) const
  {
    countLinkWork(query.linkIndex, LinkCounter::kCandidates);
    if constexpr (Shape) {
      if (!shapeAdmits(columns, source.at.row, source.point.x, source.point.y, source.point.z, target.x, target.y, target.z) ||
          !shapeAdmits(columns, targetRow, source.point.x, source.point.y, source.point.z, target.x, target.y, target.z)) {
        countWork(WorkCounter::kTrackletShape);
        return false;
      }
    }
    // Counting only: the candidate against the query's own window.
    [[maybe_unused]] bool inWindow = true; // NOLINT(misc-const-correctness): written with the work counters
    if constexpr (kWorkCounters) {
      inWindow = spatiallyCompatible(query.window, columns.phi[targetRow], secondCoordinate<Det>(query.targetLayer, target.r, target.z));
      if (inWindow) {
        countWork(WorkCounter::kTrackletInWindow);
      }
    }
    if (!accepted(source.at.layer, source.box, source.point, target, cachedVertices, source.cachedCount, source.uncachedRange, intervals,
                  secondCoordinate<Det>(query.targetLayer, target.r, target.z), mParameters->betweenBound[query.linkIndex],
                  mParameters->betweenPassives[query.linkIndex], between, inWindow)) {
      return false;
    }
    if (query.holes != 0U &&
        segmentLeavesHole(query.holes, mParameters->barrelLayers, mParameters->layerSpans, source.point.r, source.point.z, target.r, target.z)) {
      countWork(WorkCounter::kTrackletHoles);
      return false;
    }
    if (!inWindow) {
      countWork(WorkCounter::kTrackletAcceptedOutside);
    }
    return true;
  }

  /// The tracklet of an accepted candidate (`emit`'s record).
  template <typename Source, typename Target>
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Tracklet record(Index sourceRow, Index targetRow, std::size_t linkIndex, const Source& source,
                                                           const Target& target) const
  {
    const Scalar dx = target.x - source.x;
    const Scalar dy = target.y - source.y;
    const Scalar dz = target.z - source.z;
    const Scalar transverse = o2::gpu::CAMath::Hypot(dx, dy);
    return Tracklet{
      .source = sourceRow, .target = targetRow, .link = static_cast<std::uint16_t>(linkIndex), .lambda = atan2Of(dz, transverse), .phi = atan2Of(dy, dx)};
  }

 private:
  /// A vertex-list query's second-coordinate runs: the cached vertices'
  /// intervals (\a intervals), sorted, merged and clipped to the direct part
  /// of the target; their count.
  GPUhdi() GPUAlwaysInline() std::uint32_t vertexRuns(const SourceContext& source, const CompatibleVertex* cachedVertices, const QueryContext& query,
                                                      SecondRange* intervals, SecondRun* runs) const
  {
    for (std::uint32_t i = 0; i < source.cachedCount; ++i) {
      const auto& vertex = mVertices[cachedVertices[i].vertexIndex];
      intervals[i] = vertexSecondRange(mParameters->windows[query.linkIndex], mVertexWindows->windows[query.linkIndex], source.point.r, source.point.z,
                                       vertex.z, vertex.varZ, mParameters->pairCut.nSigma);
      runs[i] = {.low = intervals[i].centre - intervals[i].halfWidth, .high = intervals[i].centre + intervals[i].halfWidth};
    }
    // With one slot there is nothing to sort; compiling the loop out also
    // spares GCC's -Warray-bounds a subscript it cannot prove unreachable.
    if constexpr (kVertexCapacity > 1U) {
      insertionSort(runs, source.cachedCount, [](const SecondRun& a, const SecondRun& b) { return a.low < b.low; });
    }
    std::uint32_t runCount = 0U;
    for (std::uint32_t i = 0; i < source.cachedCount; ++i) {
      if (runCount == 0U || runs[i].low > runs[runCount - 1U].high) {
        runs[runCount++] = runs[i];
      } else if (runs[i].high > runs[runCount - 1U].high) {
        runs[runCount - 1U].high = runs[i].high;
      }
    }
    // The direct part of the target, as for the window (`prepareQuery`).
    if (query.holes != 0U) {
      std::uint32_t kept = 0U;
      for (std::uint32_t i = 0; i < runCount; ++i) {
        const SecondRun clipped{.low = o2::gpu::CAMath::Max(runs[i].low, query.directLow), .high = o2::gpu::CAMath::Min(runs[i].high, query.directHigh)};
        if (clipped.low <= clipped.high) {
          runs[kept++] = clipped;
        }
      }
      runCount = kept;
    }
    return runCount;
  }

  /// \a betweenBound: the link's bound on the kink of a passive scatterer
  /// between its layers (`IterationParameters::betweenBound`), 0 if none can
  /// lie there. With one, the pair's transverse test runs at the bound first,
  /// and the pair's own kink (`betweenScattering`, a loop over the passives)
  /// is computed only if it passes: the test only widens with the kink, so
  /// the same pairs pass.
  /// \a betweenPassives: the passives a pair on the link can cross.
  [[nodiscard]] GPUhdi() bool accepted(LayerID sourceLayer, const TimeEstPS& sourceBox, const PairPoint& source, const PairPoint& target,
                                       const CompatibleVertex* cachedVertices, std::uint32_t cachedCount, VertexRange uncachedRange,
                                       const SecondRange* intervals, Scalar targetSecond, Scalar betweenBound, std::uint64_t betweenPassives,
                                       const BetweenSource<Det>& between, [[maybe_unused]] bool inWindow) const
  {
    // The pair's own kink: on the beam line from the source's part, computed
    // once per source; else from the pair, per vertex.
    const auto betweenOf = [&](const VertexSeed& vertex) {
      if constexpr (BeamLine) {
        return betweenScattering<Det>(between, mParameters->pairCut, vertex, target, betweenPassives);
      } else {
        return betweenScattering<Det>(mParameters->scattering, mParameters->pairCut, vertex, source, target, betweenPassives);
      }
    };
    // What the pair tests compute of the two positions alone, once.
    const auto chord = pairChord(source, target);
    const auto compatible = [&](const VertexSeed& vertex, Scalar scattering) {
      const auto& cut = mParameters->pairCut;
      if constexpr (Displaced) {
        return displacedCompatible(vertex, source, scattering, target, betweenBound > 0.F ? betweenOf(vertex) : 0.F);
      }
      // Counting only: the cascade's stages split by the window.
      const auto countStages = [&](Scalar betweenScatter) {
        if constexpr (kWorkCounters) {
          if (inWindow) {
            countWork(WorkCounter::kPairTestsIn);
            if (transverseCompatible(cut, transversePairTerms(cut, vertex, source, scattering, target, chord, betweenScatter))) {
              countWork(WorkCounter::kPairLongitudinalTestsIn);
            }
          }
        }
      };
      if (!(betweenBound > 0.F)) {
        countStages(0.F);
        return pairCompatible(cut, vertex, source, scattering, target, chord);
      }
      if (!transverseCompatible(cut, transversePairTerms(cut, vertex, source, scattering, target, chord, betweenBound))) {
        countWork(WorkCounter::kPairBetweenBound);
        if constexpr (kWorkCounters) {
          if (inWindow) {
            countWork(WorkCounter::kPairBetweenBoundIn);
          }
        }
        return false;
      }
      const Scalar own = betweenOf(vertex);
      countStages(own);
      return pairCompatible(cut, vertex, source, scattering, target, chord, own);
    };
    for (std::uint32_t i = 0; i < cachedCount; ++i) {
      if constexpr (!BeamLine) {
        if (intervals != nullptr && !secondRangeContains(intervals[i], targetSecond)) {
          continue;
        }
      }
      const auto& cached = cachedVertices[i];
      const auto& vertex = mVertices[cached.vertexIndex];
      if (compatible(vertex, cached.scattering)) {
        return true;
      }
    }
    if constexpr (!BeamLine) {
      for (std::uint32_t i = uncachedRange.begin; i < uncachedRange.end; ++i) {
        const auto& vertex = mVertices[i];
        if (!vertexUsable(*mParameters, sourceBox, vertex)) {
          continue;
        }
        const Scalar scattering = sourceScattering<Det>(mParameters->scattering, mParameters->pairCut, vertex, sourceLayer, source);
        if (compatible(vertex, scattering)) {
          return true;
        }
      }
    }
    return false;
  }

  /// A displaced iteration's pair test (`displacedPairOf`), counted.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool displacedCompatible(const VertexSeed& vertex, const PairPoint& source, Scalar scattering,
                                                                    const PairPoint& target, Scalar between) const
  {
    countWork(WorkCounter::kPairTests);
    const auto result = displacedPairOf(*mParameters, vertex, source, scattering, target, between);
    if (!result.transverse) {
      return false;
    }
    countWork(WorkCounter::kPairTransverse);
    if (!result.longitudinal) {
      return false;
    }
    countWork(WorkCounter::kPairLongitudinal);
    return true;
  }

  template <typename Source, typename Target, typename Emitter>
  GPUhdi() GPUAlwaysInline() void emit(Index sourceRow, Index targetRow, std::size_t linkIndex, const Source& source, const Target& target,
                                       Emitter& emitter) const
  {
    countWork(WorkCounter::kTrackletAccepted);
    countLinkWork(linkIndex, LinkCounter::kTracklets);
    emitter.emit(record(sourceRow, targetRow, linkIndex, source, target));
  }

  /// The vertices, sorted by their boxes' starts, that can overlap \a box:
  /// starting before its end, and no earlier than the widest vertex box
  /// before its start.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() VertexRange compatibleVertexRange(const TimeEstPS& box) const
  {
    if (mVertexCount == 0U) {
      return {};
    }
    const std::int64_t low = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(mMaxVertexWidthPS);
    return {.begin = startingBefore(low), .end = startingBefore(box.upper())};
  }

  /// The first vertex whose box starts at or after \a timePS.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t startingBefore(std::int64_t timePS) const
  {
    return o2::gpu::CAAlgo::partitionPoint(0U, mVertexCount, [&](std::uint32_t vertex) {
      const std::int64_t start = mVertices[vertex].timePS;
      return start < timePS;
    });
  }

  const IterationParameters<Det>* mParameters = nullptr;
  const LayerRowMap<Det>* mSources = nullptr;
  const TimeFrame<Det>* mFrame = nullptr;
  BlockRange mTargetBlocks{};
  SourceGate mGate{};
  const std::uint8_t* mUsed = nullptr; ///< the used mask, per row
  const VertexWindowParameters<Det>* mVertexWindows = nullptr;
  const VertexSeed* mVertices = nullptr;
  std::uint32_t mVertexCount = 0U;
  std::uint32_t mMaxVertexWidthPS = 0U;
  std::array<std::uint32_t, kLayerCount> mMaxCandidateWidthPS{};
};

/// The tracklet stage a tracker runs for an iteration: the beam-line one
/// for a beam-line iteration, its displaced one for a displaced iteration,
/// the general one otherwise; each testing the shapes with \a Shape.
template <DetectorLayout Det, bool BeamLine, bool Displaced = false, bool Shape = false>
using TrackletStageFor = std::conditional_t<BeamLine, TrackletStage<Det, 1U, true, Displaced, Shape>, TrackletStage<Det, 8U, false, false, Shape>>;

/// \a body called with two `std::bool_constant`s: whether the iteration is a
/// beam-line one (\a source), and whether it is displaced (\a displaced, on
/// the beam line only): the choice of `TrackletStageFor`, host code.
template <typename Body>
decltype(auto) withTrackletVertices(VertexSource source, bool displaced, const Body& body)
{
  if (source == VertexSource::kBeamLine) {
    return displaced ? body(std::true_type{}, std::true_type{}) : body(std::true_type{}, std::false_type{});
  }
  return body(std::false_type{}, std::false_type{});
}

/// \a body called with a `std::bool_constant`: whether the iteration tests
/// the shapes (\a shapeTest), the choice of the stages' Shape, host code.
template <typename Body>
decltype(auto) withShapeTest(bool shapeTest, const Body& body)
{
  return shapeTest ? body(std::true_type{}) : body(std::false_type{});
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKLETSTAGE_H_
