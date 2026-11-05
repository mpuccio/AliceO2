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

/// \file TripletStage.h
/// \brief The triplet stage's per-first-tracklet function.

#ifndef GTRK_TRACKING_TRIPLETSTAGE_H_
#define GTRK_TRACKING_TRIPLETSTAGE_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LayerRowMap.h"
#include "GenericTracking/Tracking/ShapeCut.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/TrackletPaths.h"
#include "GenericTracking/Tracking/Triplet.h"
#include "GenericTracking/Tracking/TripletCut.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonAlgorithm.h"

#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace o2::gtrk
{

/// Generate every accepted continuation of one first tracklet.
///
/// The device path reads the iteration's plain link-topology inputs, not host
/// constexpr storage. Within each continuing link, the middle measurement's
/// tracklets are searched by their sorted (link, lambda, target) key. Elements
/// index a permutation of the first-tracklet arena in arena order
/// (TripletPermutation.h).
///
/// The work in the pieces a flat run keeps apart:
/// `prepareFirst` once per first tracklet, `forEachCandidate` its dip-angle
/// windows' candidates in the serial order, `testCandidate` one candidate's
/// cuts and fit, `record` its triplet. `operator()` runs them in sequence.
///
/// With \a Shape, the iteration's `shapeTest`: the line from the first
/// measurement to the third must cross all three at a local slope in their
/// slope intervals (ShapeCut.h), before the fit; the tracklet stage held each pair to its own line.
template <DetectorLayout Det, bool Shape = false>
class TripletStage
{
 public:
  using Index = typename StageOutput<TripletT>::Index;
  using Sources = LayerRowMap<Det>;
  using Tracklets = StageOutput<TrackletT>;

  TripletStage(const IterationParameters<Det>& parameters, const Sources& sources, const TimeFrame<Det>& frame, const Tracklets& tracklets,
               const TripletTargetBounds* targetBounds, const Index* permutation, TrackletPathBounds paths)
    : TripletStage(&parameters, &sources, &frame, tracklets.view(), tracklets.offsets().data(), targetBounds, permutation, paths)
  {
    if (Shape != parameters.shapeTest) {
      throw std::invalid_argument("the shape-testing triplet stage serves exactly the iterations with shapeTest");
    }
  }

  /// The same stage over inputs wherever they are -- a device's copies:
  /// every argument is only stored.
  TripletStage(const IterationParameters<Det>* parameters, const Sources* sources, const TimeFrame<Det>* frame, TrackletConstView tracklets,
               const Index* trackletOffsets, const TripletTargetBounds* targetBounds, const Index* permutation, TrackletPathBounds paths)
    : mParameters(parameters),
      mSources(sources),
      mFrame(frame),
      mTracklets(tracklets),
      mOffsets(trackletOffsets),
      mTargetBounds(targetBounds),
      mPermutation(permutation),
      mPaths(paths)
  {
  }

  /// What a first tracklet's tests read: the first and middle points, the
  /// middle's layer and the path prefix; and what its windows need: the
  /// middle's tracklet range and the first's lambda.
  struct FirstContext {
    TimedPoint first{};
    TimedPoint middle{};
    Index firstIndex = 0U;
    Index begin = 0U; ///< the middle's tracklets [begin, end)
    Index end = 0U;
    Scalar lambda = 0.F;
    std::uint32_t prefix = 0U;
    LayerID middleLayer = kInvalidLayerID;
  };

  template <typename Emitter>
  GPUhdi() void operator()(Index element, Emitter& emitter) const
  {
    FirstContext context;
    if (!prepareFirst(element, context)) {
      return;
    }
    const auto one = [&](Index candidateIndex) {
      TripletFitResult fit;
      bool kinked = false;
      if (testCandidate(context, candidateIndex, fit, kinked)) {
        emitter.emit(record(context, candidateIndex, fit, kinked));
      }
    };
    forEachWindowCandidate(context, one);
  }

  /// \a element's context: false when it has no continuation -- its middle
  /// outside the slice, or without a tracklet.
  [[nodiscard]] GPUhdi() bool prepareFirst(Index element, FirstContext& context) const
  {
    context.firstIndex = mPermutation[element];
    const std::size_t firstLink = mTracklets.link[context.firstIndex];
    context.middleLayer = mParameters->linkTargets[firstLink];
    const Index middleRow = mTracklets.target[context.firstIndex];
    countWork(WorkCounter::kTripletFirsts);
    if (!mSources->contains(context.middleLayer, middleRow)) {
      countWork(WorkCounter::kTripletNoMiddle);
      return false;
    }
    const Index middleElement = mSources->element(context.middleLayer, middleRow);
    context.begin = mOffsets[middleElement];
    context.end = mOffsets[middleElement + 1U];
    if (context.begin == context.end) {
      return false;
    }
    // The first tracklet is on a path of minClusters (TripletPermutation);
    // each third must be too.
    const Index sourceRow = mTracklets.source[context.firstIndex];
    context.prefix = mPaths.prefix[mSources->element(mParameters->linkSources[firstLink], sourceRow)];
    const auto columns = mFrame->measurements();
    context.first = timedPointAt(columns, sourceRow);
    context.middle = timedPointAt(columns, middleRow);
    context.lambda = mTracklets.lambda[context.firstIndex];
    return true;
  }

  /// \a body(candidateIndex) for every tracklet of the middle in the first
  /// tracklet's dip-angle window on each of the middle's out-links, in the
  /// serial order: the links in `outLinkIndices` order, each window's
  /// indices ascending from the lower bound of `lambda - halfWidth` while
  /// `lambda <= high` (a scan: cheaper than a search for the end on the host,
  /// the same on a device).
  template <typename Body>
  GPUhdi() void forEachWindowCandidate(const FirstContext& context, Body& body) const
  {
    // The dip-angle window's first-tracklet part, once for every out-link.
    const auto firstTerms = tripletFirstTerms(*mParameters, context.middleLayer, context.first.point, context.middle.point);
    for (auto out = mParameters->outLinkOffsets[context.middleLayer]; out < mParameters->outLinkOffsets[static_cast<std::size_t>(context.middleLayer) + 1U];
         ++out) {
      const std::size_t secondLink = mParameters->outLinkIndices[out];
      const Index linkBegin = lowerBoundLink(context.begin, context.end, secondLink);
      const Index linkEnd = upperBoundLink(linkBegin, context.end, secondLink);
      if (linkBegin == linkEnd) {
        continue;
      }
      const LayerID targetLayer = mParameters->linkTargets[secondLink];
      const Scalar halfWidth = tripletLambdaHalfWidth(*mParameters, secondLink, firstTerms, mTargetBounds[targetLayer]);
      const Scalar high = context.lambda + halfWidth;
      for (Index candidateIndex = lowerBoundLambda(linkBegin, linkEnd, context.lambda - halfWidth);
           candidateIndex < linkEnd && mTracklets.lambda[candidateIndex] <= high; ++candidateIndex) {
        body(candidateIndex);
      }
    }
  }

  /// The candidates of \a element's windows, counted without the tests: what
  /// a flat candidate list holds for the first tracklet.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t candidateCount(Index element) const
  {
    FirstContext context;
    if (!prepareFirst(element, context)) {
      return 0U;
    }
    std::uint32_t count = 0U;
    const auto one = [&](Index) { ++count; };
    forEachWindowCandidate(context, one);
    return count;
  }

  /// One candidate's cuts: the path, the time, the curvature and the fit
  /// (`tripletCompatible`, with the kinked retry); \a fit and \a kinked where
  /// it passes.
  [[nodiscard]] GPUhdi() bool testCandidate(const FirstContext& context, Index candidateIndex, TripletFitResult& fit, bool& kinked) const
  {
    countWork(WorkCounter::kTripletWindow);
    const LayerID targetLayer = mParameters->linkTargets[mTracklets.link[candidateIndex]];
    const auto suffix = measurementSuffix(mPaths, *mSources, mParameters->roadEndLayers, targetLayer, mTracklets.target[candidateIndex]);
    if (suffix == 0U || context.prefix + 1U + suffix < mPaths.minClusters) {
      return false;
    }
    const auto third = timedPointAt(mFrame->measurements(), mTracklets.target[candidateIndex]);
    if constexpr (Shape) {
      const auto columns = mFrame->measurements();
      const auto& from = context.first.point;
      const auto& to = third.point;
      if (!shapeAdmits(columns, mTracklets.source[context.firstIndex], from.x, from.y, from.z, to.x, to.y, to.z) ||
          !shapeAdmits(columns, mTracklets.target[context.firstIndex], from.x, from.y, from.z, to.x, to.y, to.z) ||
          !shapeAdmits(columns, mTracklets.target[candidateIndex], from.x, from.y, from.z, to.x, to.y, to.z)) {
        countWork(WorkCounter::kTripletShape);
        return false;
      }
    }
    KinkTerms terms;
    kinked = false;
    return tripletCompatible(*mParameters, context.middleLayer, context.first, context.middle, third, fit, terms, kinked);
  }

  /// The triplet of an accepted candidate.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Triplet record(const FirstContext& context, Index candidateIndex, const TripletFitResult& fit, bool kinked) const
  {
    return Triplet{.first = context.firstIndex,
                   .second = candidateIndex,
                   .curvature = fit.curvature,
                   .varCurvature = fit.varCurvature,
                   .level = 1U,
                   .kinked = kinked ? std::uint8_t{1} : std::uint8_t{0}};
  }

 private:
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Index lowerBoundLink(Index searchBegin, Index searchEnd, std::size_t value) const
  {
    return o2::gpu::CAAlgo::partitionPoint(searchBegin, searchEnd - searchBegin, [&](Index tracklet) { return mTracklets.link[tracklet] < value; });
  }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() Index upperBoundLink(Index searchBegin, Index searchEnd, std::size_t value) const
  {
    return o2::gpu::CAAlgo::partitionPoint(searchBegin, searchEnd - searchBegin, [&](Index tracklet) { return mTracklets.link[tracklet] <= value; });
  }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() Index lowerBoundLambda(Index first, Index last, Scalar lambda) const
  {
    return o2::gpu::CAAlgo::partitionPoint(first, last - first, [&](Index tracklet) { return mTracklets.lambda[tracklet] < lambda; });
  }

  const IterationParameters<Det>* mParameters = nullptr;
  const Sources* mSources = nullptr;
  const TimeFrame<Det>* mFrame = nullptr;
  TrackletConstView mTracklets;
  const Index* mOffsets = nullptr;
  const TripletTargetBounds* mTargetBounds = nullptr;
  const Index* mPermutation = nullptr;
  TrackletPathBounds mPaths;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRIPLETSTAGE_H_
